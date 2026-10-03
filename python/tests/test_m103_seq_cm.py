"""M103: SEQ_CM (codec id 19), read bases without a reference.

The ctypes wrapper round-trips every byte value and zero-length reads,
codes the same reads to the same bytes, and rejects lengths that do not
match. Through the writers it is the blocks_v1 sequences default when a
run has no reference, and an accepted sequences override on the
whole-channel layout; both read back exactly.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

pytest.importorskip("h5py")
import h5py

from ttio import SpectralDataset
from ttio.codecs import seq_cm
from ttio.enums import Compression
from ttio.written_genomic_run import WrittenGenomicRun

needs_native = pytest.mark.skipif(not seq_cm.HAVE_NATIVE_LIB,
                                  reason="libttio_rans not loaded (TTIO_RANS_LIB_PATH)")


def _reads(seed: int, n: int, *, extra: bytes = b""):
    """Overlapping reads from a random genome, with N and lower case."""
    rng = np.random.default_rng(seed)
    genome = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), 20_000)
    lengths = rng.integers(40, 160, n).astype(np.uint64)
    lengths[3] = 0
    parts = []
    for i, L in enumerate(lengths):
        s = int(rng.integers(0, genome.size - 200))
        r = bytearray(genome[s:s + int(L)].tobytes())
        if L and i % 7 == 0:
            r[int(L) // 2] = ord("N")
        if L and i % 11 == 0:
            r[0] = ord("a")
        parts.append(bytes(r))
    seq = b"".join(parts) + extra
    if extra:
        lengths[-1] += len(extra)
    return seq, lengths


@needs_native
def test_round_trip_every_byte_and_empty_reads():
    seq, lengths = _reads(3, 400, extra=bytes(range(256)))
    blob = seq_cm.encode(seq, lengths)
    assert blob[:4] == b"SQC1"
    assert seq_cm.decode(blob, lengths) == seq
    assert seq_cm.encode(seq, lengths) == blob  # deterministic


@needs_native
def test_empty_block():
    blob = seq_cm.encode(b"", np.zeros(0, dtype=np.uint64))
    assert seq_cm.decode(blob, np.zeros(0, dtype=np.uint64)) == b""


@needs_native
def test_lengths_must_match():
    seq, lengths = _reads(5, 50)
    with pytest.raises(ValueError):
        seq_cm.encode(seq, lengths[:-1])
    blob = seq_cm.encode(seq, lengths)
    with pytest.raises(ValueError):
        seq_cm.decode(blob, lengths[:-1])


def _run(seed: int, n: int, **kw) -> WrittenGenomicRun:
    seq, lengths = _reads(seed, n)
    rng = np.random.default_rng(seed + 1)
    lengths = lengths.astype(np.uint32)
    offsets = np.zeros(n, dtype=np.uint64)
    offsets[1:] = np.cumsum(lengths[:-1])
    return WrittenGenomicRun(
        acquisition_mode=7, reference_uri="", platform="ILLUMINA", sample_name="M103",
        positions=np.zeros(n, dtype=np.int64),
        mapping_qualities=np.zeros(n, dtype=np.uint8),
        flags=np.full(n, 0x4, dtype=np.uint32),
        sequences=np.frombuffer(seq, dtype=np.uint8).copy(),
        qualities=rng.integers(33, 73, len(seq)).astype(np.uint8),
        offsets=offsets, lengths=lengths,
        cigars=["*"] * n,
        read_names=[f"r{i:05d}" for i in range(n)],
        mate_chromosomes=["*"] * n,
        mate_positions=np.full(n, -1, dtype=np.int64),
        template_lengths=np.zeros(n, dtype=np.int32),
        chromosomes=["*"] * n,
        **kw,
    )


def _expected(run: WrittenGenomicRun) -> list[str]:
    seq = bytes(run.sequences.tobytes())
    out, o = [], 0
    for L in run.lengths:
        out.append(seq[o:o + int(L)].decode("ascii"))
        o += int(L)
    return out


@needs_native
def test_blocks_v1_default_without_reference(tmp_path):
    from ttio.genomic._blocks import slice_run
    from ttio.genomic.stream_writer import GenomicStreamWriter

    run = _run(7, 600)
    path = tmp_path / "blocks.tio"
    SpectralDataset.write_minimal(str(path), title="m103", isa_investigation_id="i", runs={})
    ds = SpectralDataset.open(str(path), writable=True)
    w = GenomicStreamWriter(ds.study_group, "g", acquisition_mode=run.acquisition_mode,
                            reference_uri=run.reference_uri, platform=run.platform,
                            sample_name=run.sample_name, block_reads=200)
    with ds, w:
        for s in range(0, 600, 100):
            w.append_batch(slice_run(run, s, s + 100))
    with h5py.File(path, "r") as f:
        idx = f["study/genomic_runs/g/blocks/index"][()]
        assert len(idx) == 3
        assert set(idx["sequences_codec"].tolist()) == {int(Compression.SEQ_CM)}
    with SpectralDataset.open(str(path)) as ds:
        got = [r.sequence for r in ds.genomic_runs["g"].iter_reads()]
    assert got == _expected(run)


@needs_native
def test_whole_channel_override(tmp_path):
    run = dataclasses.replace(_run(11, 300), opt_legacy_whole_channel=True,
                              signal_codec_overrides={"sequences": Compression.SEQ_CM})
    path = tmp_path / "whole.tio"
    SpectralDataset.write_minimal(str(path), title="m103", isa_investigation_id="i",
                                  runs={}, genomic_runs={"g": run})
    with h5py.File(path, "r") as f:
        ds = f["study/genomic_runs/g/signal_channels/sequences"]
        assert int(ds.attrs["compression"]) == int(Compression.SEQ_CM)
        assert bytes(ds[:4].tobytes()) == b"SQC1"
    with SpectralDataset.open(str(path)) as ds:
        got = [r.sequence for r in ds.genomic_runs["g"].iter_reads()]
    assert got == _expected(run)
