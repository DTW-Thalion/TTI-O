"""Block encoder: a block's blobs equal the whole-run writer's bytes."""
from __future__ import annotations

import numpy as np

from _genomic_fixture import make_written_genomic_run
from ttio.genomic import _blocks
from ttio.providers.memory import MemoryProvider
from ttio._dataset_write_genomic import _write_genomic_run


def test_slice_run_rebases_offsets():
    run = make_written_genomic_run(n_reads=10, read_len=7)
    s = _blocks.slice_run(run, 3, 6)
    assert len(s.lengths) == 3 and int(s.offsets[0]) == 0
    assert bytes(s.sequences) == bytes(run.sequences[21:42])
    assert s.read_names == run.read_names[3:6]
    assert s.chromosomes == run.chromosomes[3:6]


def _open_channel(sc, ch, seq_layout):
    if ch == "sequences":
        try:
            g = sc.open_group("sequences")
        except Exception:
            return sc.open_dataset("sequences")
        return g.open_dataset(seq_layout)
    if ch == "mate_info":
        return sc.open_group("mate_info").open_dataset("inline_v2")
    return sc.open_dataset(ch)


def test_encode_block_equals_whole_run_writer_bytes():
    run = make_written_genomic_run(n_reads=40, read_len=50, with_reference=True, paired=True)
    blobs = _blocks.encode_block(run)
    # blocks_v1 defaults cigars to RANS_ORDER0 (the M82 compound VL-string
    # default has no blob form); the whole-run reference write gets the
    # same override so the byte comparison is like for like.
    from ttio.enums import Compression
    import dataclasses
    ref = dataclasses.replace(run, signal_codec_overrides={"cigars": Compression.RANS_ORDER0,
                                                           "qualities": Compression.FQZCOMP_NX16_Z})
    prov = MemoryProvider.open("memory://blocks-enc-test", mode="w")
    root = prov.root_group()
    _write_genomic_run(root, "r", ref)
    sc = root.open_group("r").open_group("signal_channels")
    for ch in _blocks.REQUIRED_BLOCK_CHANNELS:      # the fixture carries no tags
        ds = _open_channel(sc, ch, blobs.seq_layout)
        assert bytes(np.asarray(ds.read(), dtype=np.uint8).tobytes()) == blobs.blobs[ch], ch
        assert int(ds.get_attribute("compression")) == blobs.compression[ch], ch
    assert blobs.seq_layout == "refdiff_v2"
    assert blobs.n_reads == 40 and blobs.n_bases == 2000


def test_encode_block_uses_preassigned_chrom_ids():
    run = make_written_genomic_run(n_reads=8, read_len=10,
                                   chromosomes=["chrB"] * 4 + ["chrA"] * 4)
    run.chrom_name_to_id = {"chrA": 0}
    _blocks.encode_block(run)
    assert run.chrom_name_to_id == {"chrA": 0, "chrB": 1}


def test_encode_block_forwards_qual_strategy_hint():
    from ttio.codecs import fqzcomp_nx16_z as fq
    run = make_written_genomic_run(n_reads=40, read_len=50)
    default = _blocks.encode_block(run)
    # Small block: the size floor keeps auto on V4.
    assert fq.stream_strategy(default.blobs["qualities"]) == 4
    forced = _blocks.encode_block(run, qual_strategy_hint=5)
    assert fq.stream_strategy(forced.blobs["qualities"]) == 5
    # An explicit -1 is the default path, byte for byte.
    explicit = _blocks.encode_block(run, qual_strategy_hint=-1)
    assert explicit.blobs["qualities"] == default.blobs["qualities"]
