"""Read-side support for the ``blocks_v1`` genomic layout.

A block is materialised as a v1.8-shaped in-memory run group (the
block's blobs with their attributes, the matching slice of the index
arrays, and the run-level name tables), and the existing
:class:`~ttio.genomic_run.GenomicRun` decode path runs over it. The
codecs never learn about blocks.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .. import _hdf5_io as io
from ..enums import Precision
from . import _blocks


@dataclass
class BlockTable:
    """The ``blocks/index`` rows as column arrays."""
    read_start: np.ndarray      # uint64
    n_reads: np.ndarray         # uint32
    base_start: np.ndarray      # uint64
    n_bases: np.ndarray         # uint64
    ranges: dict[str, tuple[np.ndarray, np.ndarray]]   # channel -> (off, len)
    codecs: dict[str, np.ndarray] | None = None        # channel -> codec id per block
    # The channels whose columns the index carries: the required five,
    # plus "tags" when the run has the tags triple (M101).
    channels: tuple = _blocks.REQUIRED_BLOCK_CHANNELS

    @property
    def count(self) -> int:
        return int(len(self.read_start))

    @property
    def read_count(self) -> int:
        if self.count == 0:
            return 0
        return int(self.read_start[-1]) + int(self.n_reads[-1])

    def block_for(self, i: int) -> int:
        b = int(np.searchsorted(self.read_start, i, side="right")) - 1
        if b < 0 or i >= int(self.read_start[b]) + int(self.n_reads[b]):
            raise IndexError(f"read index {i} out of range [0, {self.read_count})")
        return b

    @classmethod
    def read(cls, run_group) -> "BlockTable":
        rows = run_group.open_group("blocks").open_dataset("index").read_rows()
        def col(name, dt):
            return np.asarray([r[name] for r in rows], dtype=dt)

        # An optional channel's columns are absent from runs that do
        # not carry it: every block then has an empty range.
        def present(ch):
            return ch in _blocks.REQUIRED_BLOCK_CHANNELS or not rows or f"{ch}_off" in rows[0]

        def zeros(dt):
            return np.zeros(len(rows), dtype=dt)
        ranges = {ch: ((col(f"{ch}_off", np.uint64), col(f"{ch}_len", np.uint64))
                       if present(ch) else (zeros(np.uint64), zeros(np.uint64)))
                  for ch in _blocks.BLOCK_CHANNELS}
        codecs = None
        if rows and f"{_blocks.BLOCK_CHANNELS[0]}_codec" in rows[0]:
            codecs = {ch: col(f"{ch}_codec", np.uint32) if present(ch) else zeros(np.uint32)
                      for ch in _blocks.BLOCK_CHANNELS}
        channels = tuple(ch for ch in _blocks.BLOCK_CHANNELS
                         if ch in _blocks.REQUIRED_BLOCK_CHANNELS
                         or (rows and f"{ch}_off" in rows[0]))
        return cls(read_start=col("read_start", np.uint64), n_reads=col("n_reads", np.uint32),
                   base_start=col("base_start", np.uint64), n_bases=col("n_bases", np.uint64),
                   ranges=ranges, codecs=codecs, channels=channels)


_INDEX_ARRAYS = (("lengths", Precision.UINT32), ("positions", Precision.INT64),
                 ("mapping_qualities", Precision.UINT8), ("flags", Precision.UINT32),
                 ("chromosome_ids", Precision.UINT16))


def _copy_attrs(src, dst) -> None:
    for k in src.attribute_names():
        dst.set_attribute(k, src.get_attribute(k))


def _rows_of(group, name) -> list[dict]:
    return [{"name": (r["name"].decode() if isinstance(r["name"], bytes) else r["name"])}
            for r in io.read_compound_dataset(group, name)]


def materialise_block(run_group, table: BlockTable, b: int, *,
                      chrom_name_rows: list[dict] | None = None,
                      mate_chrom_rows: list[dict] | None = None,
                      skip_channels: tuple = ()):
    """Return an in-memory v1.8-shaped run group holding block ``b``.

    ``skip_channels`` names blob channels to leave out of the view —
    the per-AU decrypt walker skips the encrypted (deleted) channels
    and injects their decrypted raw bytes itself.
    """
    from ..providers.memory import MemoryProvider

    root = MemoryProvider.open(f"memory://ttio-block-view-{id(run_group)}-{b}", mode="w").root_group()
    view = root.create_group("run")
    for k in run_group.attribute_names():
        if k in ("layout", "block_policy", "base_count"):
            continue          # the view is a v1.8-shaped whole-channel run
        view.set_attribute(k, run_group.get_attribute(k))
    r0, n = int(table.read_start[b]), int(table.n_reads[b])
    io.write_int_attr(view, "read_count", n)

    # genomic_index slice
    src_idx = run_group.open_group("genomic_index")
    dst_idx = view.create_group("genomic_index")
    for name, prec in _INDEX_ARRAYS:
        src = src_idx.open_dataset(name)
        arr = src.read(r0, n)
        ds = dst_idx.create_dataset(name, prec, n)
        ds.write(arr)
        _copy_attrs(src, ds)
    if chrom_name_rows is None:
        chrom_name_rows = _rows_of(src_idx, "chromosome_names")
    io.write_compound_dataset(dst_idx, "chromosome_names", chrom_name_rows,
                              [("name", io.vl_str())])

    # signal channels
    src_sc = run_group.open_group("signal_channels")
    dst_sc = view.create_group("signal_channels")
    from ..enums import Compression
    for ch in _blocks.BLOCK_CHANNELS:
        if ch in skip_channels:
            continue
        off, ln = int(table.ranges[ch][0][b]), int(table.ranges[ch][1][b])
        if ln == 0:
            continue
        codec = int(table.codecs[ch][b]) if table.codecs is not None else None
        if ch == "sequences":
            g = src_sc.open_group("sequences")
            src = g.open_dataset("data")
            if codec is None:
                codec = int(src.get_attribute("compression"))
            if codec == int(Compression.REF_DIFF_V2):
                dst_parent, dst_name = dst_sc.create_group("sequences"), "refdiff_v2"
            else:
                dst_parent, dst_name = dst_sc, "sequences"
        elif ch == "mate_info":
            src = src_sc.open_group("mate_info").open_dataset("inline_v2")
            dst_parent, dst_name = dst_sc.create_group("mate_info"), "inline_v2"
        else:
            src = src_sc.open_dataset(ch)
            dst_parent, dst_name = dst_sc, ch
        ds = dst_parent.create_dataset(dst_name, Precision.UINT8, ln)
        ds.write(np.asarray(src.read(off, ln), dtype=np.uint8))
        _copy_attrs(src, ds)
        if codec is not None:
            io.write_int_attr(ds, "compression", codec, dtype="<u1")
    if src_sc.has_child("mate_info"):
        if mate_chrom_rows is None:
            mate_chrom_rows = _rows_of(src_sc.open_group("mate_info"), "chrom_names")
        if not dst_sc.has_child("mate_info"):
            dst_sc.create_group("mate_info")
        io.write_compound_dataset(dst_sc.open_group("mate_info"), "chrom_names",
                                  mate_chrom_rows, [("name", io.vl_str())])
    return view


class LazyGenomicIndex:
    """:class:`~ttio.genomic_index.GenomicIndex` look-alike for
    ``blocks_v1`` runs: ``count`` comes from the block table, the
    per-read arrays load from disk on first use."""

    def __init__(self, idx_group, table: BlockTable):
        self._idx_group = idx_group
        self._table = table
        self._loaded = None

    @property
    def count(self) -> int:
        return self._table.read_count

    def _load(self):
        if self._loaded is None:
            from ..genomic_index import GenomicIndex
            self._loaded = GenomicIndex.read(self._idx_group)
        return self._loaded

    def __getattr__(self, name):
        # offsets, lengths, chromosomes, positions, mapping_qualities,
        # flags, chromosome_ids, chromosome_names, indices_for_*
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._load(), name)
