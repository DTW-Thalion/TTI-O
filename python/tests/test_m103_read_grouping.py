"""M103: reads grouped by sequence before blocking (layout blocks_v1_grouped).

The writer reorders an unaligned run's reads with the native kernel, cuts
blocks over the reordered run and stores genomic_index/input_index; every
read accessor presents input order, block iteration stays in stored order,
and the grouping survives per-AU encryption, the encrypted transport
(BlockSidecar input_index, transport-spec v0.13) and signatures. Plaintext
transport delivers input order and the receiver writes it ungrouped.
"""
from __future__ import annotations

import dataclasses
import shutil

import numpy as np
import pytest

pytest.importorskip("h5py")
import h5py

from ttio import SpectralDataset
from ttio.codecs import seq_group
from ttio.genomic._blocks import slice_run
from ttio.genomic.stream_writer import GenomicStreamWriter
from ttio.written_genomic_run import WrittenGenomicRun

needs_native = pytest.mark.skipif(not seq_group.HAVE_NATIVE_LIB,
                                  reason="libttio_rans not loaded (TTIO_RANS_LIB_PATH)")
KEY = bytes(range(32))


def _comp(seq: bytes) -> bytes:
    t = bytes.maketrans(b"ACGT", b"TGCA")
    return seq.translate(t)[::-1]


def _run(seed: int = 7, n: int = 3000, L: int = 120, G: int = 30_000) -> WrittenGenomicRun:
    """Paired reads from both strands of a random genome, shuffled, with a
    few N bases and a zero-length read: about 12x coverage."""
    rng = np.random.default_rng(seed)
    genome = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), G).tobytes()
    seqs, names = [], []
    for i in range(n // 2):
        p = int(rng.integers(0, G - 500))
        r1 = genome[p:p + L]
        r2 = _comp(genome[p + 300:p + 300 + L])
        if i % 37 == 0:
            r1 = r1[:50] + b"N" + r1[51:]
        seqs += [r1, r2]
        names += [f"frag{i}/1", f"frag{i}/2"]
    seqs[5] = b""
    perm = rng.permutation(len(seqs))
    seqs = [seqs[i] for i in perm]
    names = [names[i] for i in perm]
    lengths = np.array([len(s) for s in seqs], dtype=np.uint32)
    offsets = np.zeros(len(seqs), dtype=np.uint64)
    offsets[1:] = np.cumsum(lengths[:-1])
    seq = b"".join(seqs)
    m = len(seqs)
    return WrittenGenomicRun(
        acquisition_mode=7, reference_uri="", platform="ILLUMINA", sample_name="M103G",
        positions=np.arange(m, dtype=np.int64),            # distinct per read, to check order
        mapping_qualities=(np.arange(m) % 61).astype(np.uint8),
        flags=np.full(m, 0x4, dtype=np.uint32),
        sequences=np.frombuffer(seq, dtype=np.uint8).copy(),
        qualities=rng.integers(33, 73, len(seq)).astype(np.uint8),
        offsets=offsets, lengths=lengths,
        cigars=["*"] * m, read_names=names,
        mate_chromosomes=["*"] * m, mate_positions=np.full(m, -1, dtype=np.int64),
        template_lengths=np.zeros(m, dtype=np.int32), chromosomes=["*"] * m,
    )


def _write(path, run, *, group: bool, block_reads: int = 700, **kw):
    SpectralDataset.write_minimal(str(path), title="m103", isa_investigation_id="i", runs={})
    ds = SpectralDataset.open(str(path), writable=True)
    w = GenomicStreamWriter(ds.study_group, "g", acquisition_mode=run.acquisition_mode,
                            reference_uri=run.reference_uri, platform=run.platform,
                            sample_name=run.sample_name, block_reads=block_reads,
                            group_reads=group, **kw)
    n = len(run.lengths)
    with ds, w:
        for s in range(0, n, 500):
            w.append_batch(slice_run(run, s, min(s + 500, n)))
    return path


def _expected(run):
    seq = bytes(run.sequences.tobytes())
    out, o = [], 0
    for L in run.lengths:
        out.append(seq[o:o + int(L)].decode("ascii"))
        o += int(L)
    return out


def _reads(path):
    with SpectralDataset.open(str(path)) as ds:
        g = ds.genomic_runs["g"]
        return [(r.read_name, r.sequence, int(r.position)) for r in g.iter_reads()]


@needs_native
def test_grouped_layout_and_input_order(tmp_path):
    run = _run()
    path = _write(tmp_path / "grp.tio", run, group=True)
    n = len(run.lengths)
    with h5py.File(path, "r") as f:
        rg = f["study/genomic_runs/g"]
        layout = rg.attrs["layout"]
        assert (layout.decode() if isinstance(layout, bytes) else layout) == "blocks_v1_grouped"
        ii = rg["genomic_index/input_index"][()]
        assert ii.dtype == np.uint32 and sorted(ii.tolist()) == list(range(n))
        assert ii.tolist() != list(range(n)), "the reads were reordered"
        stored_pos = rg["genomic_index/positions"][()]
        assert stored_pos.tolist() == ii.tolist(), "every index column moved with its read"
        assert len(rg["blocks/index"]) > 1
    want = list(zip(run.read_names, _expected(run), range(n)))
    assert _reads(path) == want
    with SpectralDataset.open(str(path)) as ds:
        g = ds.genomic_runs["g"]
        assert g.layout == "blocks_v1_grouped"
        assert len(g) == n
        for i in (0, 1, 5, n // 2, n - 1, -1):
            r = g[i]
            j = i % n
            assert (r.read_name, r.sequence) == (run.read_names[j], want[j][1])
        assert g.index.positions.tolist() == list(range(n)), "index in input order"
        assert g.index.lengths.tolist() == run.lengths.tolist()
        part = [r.read_name for r in g.iter_reads(1000, 1100)]
        assert part == run.read_names[1000:1100]
        # Block iteration is in stored order; input_index maps its rows.
        seen = {}

        def visit(view, view_start, first_read, n_reads):
            for k in range(n_reads):
                seen[int(g.input_index[first_read + k])] = view[view_start + k].read_name

        g.for_each_block(visit, threads=1)
        assert [seen[i] for i in range(n)] == run.read_names


@needs_native
def test_grouping_shrinks_sequences(tmp_path):
    run = _run(n=6000)
    a = _write(tmp_path / "plain.tio", run, group=False, block_reads=1000)
    b = _write(tmp_path / "grp.tio", run, group=True, block_reads=1000)
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        sa = fa["study/genomic_runs/g/signal_channels/sequences/data"].shape[0]
        sb = fb["study/genomic_runs/g/signal_channels/sequences/data"].shape[0]
    assert sb < sa * 0.8, (sa, sb)


@needs_native
def test_write_minimal_opt_group_reads(tmp_path):
    run = dataclasses.replace(_run(n=800), opt_group_reads=True)
    path = tmp_path / "wm.tio"
    SpectralDataset.write_minimal(str(path), title="m103", isa_investigation_id="i",
                                  runs={}, genomic_runs={"g": run})
    with SpectralDataset.open(str(path)) as ds:
        g = ds.genomic_runs["g"]
        assert g.layout == "blocks_v1_grouped"
        assert [r.read_name for r in g.iter_reads()] == run.read_names


@needs_native
def test_grouping_refuses_aligned_and_legacy(tmp_path):
    SpectralDataset.write_minimal(str(tmp_path / "x.tio"), title="m103",
                                  isa_investigation_id="i", runs={})
    with SpectralDataset.open(str(tmp_path / "x.tio"), writable=True) as ds:
        common = dict(acquisition_mode=7, reference_uri="r", platform="ILLUMINA",
                      sample_name="s", group_reads=True)
        with pytest.raises(ValueError):
            GenomicStreamWriter(ds.study_group, "a", reference_chrom_seqs={"chr1": b"ACGT"},
                                **common)
        with pytest.raises(ValueError):
            GenomicStreamWriter(ds.study_group, "b", opt_legacy_whole_channel=True, **common)


@needs_native
def test_per_au_restore_keeps_grouping(tmp_path):
    from ttio.encryption_per_au import decrypt_per_au_in_place, encrypt_per_au
    run = _run(n=1500)
    path = _write(tmp_path / "grp.tio", run, group=True, block_reads=400)
    pristine = tmp_path / "pristine.tio"
    shutil.copyfile(path, pristine)
    encrypt_per_au(str(path), KEY)
    decrypt_per_au_in_place(str(path), KEY)
    with h5py.File(path, "r") as f, h5py.File(pristine, "r") as p:
        for name in ("blocks/index", "genomic_index/input_index",
                     "signal_channels/sequences/data", "signal_channels/qualities"):
            assert np.array_equal(f[f"study/genomic_runs/g/{name}"][()],
                                  p[f"study/genomic_runs/g/{name}"][()]), name
    assert [r[0] for r in _reads(path)] == run.read_names


@needs_native
def test_encrypted_transport_carries_input_index(tmp_path):
    from ttio.encryption_per_au import decrypt_per_au_in_place, encrypt_per_au
    from ttio.transport.codec import TransportWriter
    from ttio.transport.encrypted import read_encrypted_to_file, write_encrypted_dataset
    run = _run(n=1500)
    path = _write(tmp_path / "grp.tio", run, group=True, block_reads=400)
    with h5py.File(path, "r") as f:
        want_ii = f["study/genomic_runs/g/genomic_index/input_index"][()]
    encrypt_per_au(str(path), KEY)
    stream = tmp_path / "s.tis"
    with TransportWriter(str(stream)) as writer:
        write_encrypted_dataset(writer, str(path))
    got = tmp_path / "got.tio"
    read_encrypted_to_file(str(stream), str(got))
    with h5py.File(got, "r") as f:
        rg = f["study/genomic_runs/g"]
        layout = rg.attrs["layout"]
        assert (layout.decode() if isinstance(layout, bytes) else layout) == "blocks_v1_grouped"
        assert np.array_equal(rg["genomic_index/input_index"][()], want_ii)
        feats = f.attrs.get("ttio_features")
        assert b"transport_blocks_v1_grouped" not in (feats if isinstance(feats, bytes)
                                                      else str(feats).encode())
    decrypt_per_au_in_place(str(got), KEY)
    assert [r[0] for r in _reads(got)] == run.read_names


@needs_native
def test_plaintext_transport_delivers_input_order(tmp_path):
    from ttio.transport.codec import file_to_transport, transport_to_file
    run = _run(n=1200)
    path = _write(tmp_path / "grp.tio", run, group=True, block_reads=300)
    stream = tmp_path / "p.tis"
    out = tmp_path / "out.tio"
    file_to_transport(str(path), str(stream))
    transport_to_file(str(stream), str(out))
    with SpectralDataset.open(str(out)) as ds:
        g = next(iter(ds.genomic_runs.values()))
        assert g.layout != "blocks_v1_grouped"
        assert [r.read_name for r in g.iter_reads()] == run.read_names


@needs_native
def test_signatures_cover_input_index(tmp_path):
    from ttio.signatures import sign_genomic_run, verify_genomic_run
    path = _write(tmp_path / "grp.tio", _run(n=600), group=True, block_reads=200)
    key = b"k" * 32
    with h5py.File(path, "r+") as f:
        rg = f["study/genomic_runs/g"]
        assert "genomic_index/input_index" in sign_genomic_run(rg, key)
        assert verify_genomic_run(rg, key)
        ii = rg["genomic_index/input_index"][...]
        ii[[0, 1]] = ii[[1, 0]]
        rg["genomic_index/input_index"][...] = ii
        assert not verify_genomic_run(rg, key)


@needs_native
def test_corrupt_input_index_is_refused(tmp_path):
    path = _write(tmp_path / "grp.tio", _run(n=600), group=True, block_reads=200)
    with h5py.File(path, "r+") as f:
        ds = f["study/genomic_runs/g/genomic_index/input_index"]
        ii = ds[...]
        ii[0] = ii[1]
        ds[...] = ii
    with pytest.raises(ValueError, match="permutation"):
        with SpectralDataset.open(str(path)) as ds:
            ds.genomic_runs["g"]
