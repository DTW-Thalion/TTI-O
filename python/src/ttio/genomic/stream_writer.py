"""``GenomicStreamWriter``: write a genomic run as ``blocks_v1`` with
bounded memory.

Reads are buffered until a block is full (``block_reads`` reads or
``block_bytes`` sequence bytes, whichever first), encoded through
:func:`ttio.genomic._blocks.encode_block`, and appended to extendable
per-channel datasets; ``blocks/index`` records where each block's
blob lives. Format: ``docs/format-spec.md`` section 10.12.

Cross-language equivalents
--------------------------
Java: ``global.thalion.ttio.genomics.GenomicStreamWriter`` ·
Objective-C: ``TTIOGenomicStreamWriter`` (both follow this module).

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import concurrent.futures as _cf
import dataclasses
import json
import os
from typing import Any

import numpy as np

from .. import _hdf5_io as io
from ..codecs import quality
from ..enums import Compression, Precision
from ..providers.base import CompoundField, CompoundFieldKind
from ..written_genomic_run import WrittenGenomicRun
from .._threads import (resolve_threads, pool_context,
                        resolve_memory_budget, apply_v6_segment_threads)
from . import _blocks

LAYOUT = "blocks_v1"
#: blocks_v1 with the reads reordered by sequence and
#: ``genomic_index/input_index`` restoring input order (M103). A distinct
#: value so readers that predate it refuse the run rather than read it
#: in stored order (format-spec 10.12.7).
LAYOUT_GROUPED = "blocks_v1_grouped"
DEFAULT_BLOCK_READS = 1_000_000
# A block is the unit of encode concurrency and of residency, so a
# big one costs both: the pool cannot start a block that is still
# filling, and every block in flight is held whole. At 256 MiB a
# 2x250 import of 10.6 M reads was 11 blocks, which left 7 workers
# at 3.05x of one thread and 10.4 GiB resident. At 64 MiB the same
# import is 43.0 s against 70.3 s and 4.50 GiB against 10.38 GiB.
# The cost is context the codecs no longer share: 0.16% on 2x250,
# 0.05% on 150 bp FASTQ, and 0.46% on HiFi, where reads are long
# enough that a block holds few of them.
DEFAULT_BLOCK_BYTES = 64 << 20
# Unfiltered chunked datasets store whole chunks, so the chunk is kept
# moderate: a 10-read run costs 5 x 256 KiB, a 100 GB channel ~400 k chunks.
CHANNEL_CHUNK = 256 << 10

#: Block index schema, in the fixed column order of format-spec 10.12.
INDEX_FIELDS: list[CompoundField] = (
    [CompoundField("read_start", CompoundFieldKind.UINT64),
     CompoundField("n_reads", CompoundFieldKind.UINT32),
     CompoundField("base_start", CompoundFieldKind.UINT64),
     CompoundField("n_bases", CompoundFieldKind.UINT64)]
    + [CompoundField(f"{ch}_{k}", CompoundFieldKind.UINT64)
       for ch in _blocks.REQUIRED_BLOCK_CHANNELS for k in ("off", "len")]
    + [CompoundField(f"{ch}_codec", CompoundFieldKind.UINT32)
       for ch in _blocks.REQUIRED_BLOCK_CHANNELS]
)
#: The schema of a run that carries the tags channel (M101): the
#: tags (off, len, codec) triple follows the required columns.
INDEX_FIELDS_WITH_TAGS: list[CompoundField] = INDEX_FIELDS + [
    CompoundField("tags_off", CompoundFieldKind.UINT64),
    CompoundField("tags_len", CompoundFieldKind.UINT64),
    CompoundField("tags_codec", CompoundFieldKind.UINT32),
]


def index_fields_for(rows_or_channels) -> list[CompoundField]:
    """The index schema for a run: with the tags triple when the given
    index rows (dicts) carry ``tags_off``, or the channel names include
    ``tags``."""
    items = list(rows_or_channels)
    if items and isinstance(items[0], dict):
        has = "tags_off" in items[0]
    else:
        has = "tags" in items
    return INDEX_FIELDS_WITH_TAGS if has else INDEX_FIELDS

_INDEX_ARRAYS = (
    ("lengths", Precision.UINT32, np.uint32),
    ("positions", Precision.INT64, np.int64),
    ("mapping_qualities", Precision.UINT8, np.uint8),
    ("flags", Precision.UINT32, np.uint32),
    ("chromosome_ids", Precision.UINT16, np.uint16),
)
#: The grouped layout's extra column: row j holds the input index of the
#: read stored at row j (M103).
_INPUT_INDEX = (("input_index", Precision.UINT32, np.uint32),)


def register_block_chromosomes(block: WrittenGenomicRun, chrom_map: dict) -> None:
    """Assign ids for every chromosome name the block introduces, in the
    order the block encoder assigns them (own names in read order, '*'
    included, then mate names that are not '', '*' or '='), so the
    encoder only reads the map and blocks can encode concurrently."""
    for name in block.chromosomes:
        if name not in chrom_map:
            chrom_map[name] = len(chrom_map)
    for name in block.mate_chromosomes:
        if name and name not in ("*", "=") and name not in chrom_map:
            chrom_map[name] = len(chrom_map)


class GenomicStreamWriter:
    """Append reads to one genomic run of an open-for-write dataset.

    With ``threads`` > 1 (the ``threads`` argument, else ``TTIO_THREADS``,
    else cores minus 2) completed blocks are encoded on a pool and written
    in order by the caller's thread; at most ``threads + 1`` blocks are
    pending or encoding, so memory is about that many block working sets
    (about 1.8 GB for a 1 M-read block of 100 bp reads, 10 GB for
    10.6 M reads of 250 bp). The file is byte for byte the one thread's.
    ``read_count`` and ``block_count`` count written blocks; blocks still
    in flight are counted after ``close``."""

    def __init__(self, study_group, run_name: str, *,
                 acquisition_mode: int, reference_uri: str, platform: str,
                 sample_name: str, reference_chrom_seqs=None,
                 embed_reference: bool = False,
                 block_reads: int = DEFAULT_BLOCK_READS,
                 block_bytes: int = DEFAULT_BLOCK_BYTES,
                 opt_disable_qualities_v5: bool = False,
                 signal_codec_overrides: dict | None = None,
                 signal_compression: str = "gzip",
                 opt_legacy_whole_channel: bool = False,
                 provenance_records=None,
                 threads: int | None = None,
                 memory_budget_bytes: int | None = None,
                 read_role: str | None = None,
                 ref_diff_slice_bytes: int = 0,
                 group_reads: bool = False):
        self._study = study_group
        self._provenance = list(provenance_records or [])
        self._name = run_name
        # M97: reject the mistuned lossy table before anything is
        # written; write_minimal re-checks per run.
        overrides = dict(signal_codec_overrides or {})
        if (
            overrides.get("qualities") == Compression.QUALITY_BINNED
            and not quality.binned_allowed_for_platform(platform)
        ):
            raise ValueError(
                "signal_codec_overrides['qualities']: QUALITY_BINNED "
                "applies the fixed Illumina-8 bin table, which does "
                "not fit the quality distribution of platform "
                f"{platform!r} (M97)."
            )
        self._meta: dict[str, Any] = dict(
            acquisition_mode=int(acquisition_mode), reference_uri=reference_uri,
            platform=platform, sample_name=sample_name,
            reference_chrom_seqs=reference_chrom_seqs,
            embed_reference=bool(embed_reference),
            opt_disable_qualities_v5=bool(opt_disable_qualities_v5),
            signal_codec_overrides=overrides,
            signal_compression=signal_compression,
            read_role=read_role,
            ref_diff_slice_bytes=int(ref_diff_slice_bytes or 0),
        )
        if block_reads < 1 or block_bytes < 1:
            raise ValueError("block_reads and block_bytes must be >= 1")
        self._block_reads = int(block_reads)
        self._block_bytes = int(block_bytes)
        self._legacy = bool(opt_legacy_whole_channel)
        # M103: buffer the whole run, reorder it by sequence at close,
        # then block it as usual (format-spec 10.12.7).
        self._group = bool(group_reads)
        if self._group and self._legacy:
            raise ValueError("group_reads needs the blocks_v1 layout, not the legacy whole-channel one")
        if self._group and reference_chrom_seqs is not None:
            raise ValueError("group_reads applies to runs without a reference; aligned "
                             "runs keep coordinate order for REF_DIFF_V2")
        if self._group:
            from ..codecs import seq_group as _sg
            if not _sg.HAVE_NATIVE_LIB:
                raise RuntimeError("group_reads requires libttio_rans (set TTIO_RANS_LIB_PATH)")
        self._group_parts: list[WrittenGenomicRun] = []
        self._input_index: np.ndarray | None = None
        self._pending: list[WrittenGenomicRun] = []
        self._pending_reads = 0
        self._pending_bytes = 0
        self._chrom_map: dict[str, int] = {}
        self._pending_chrom: str | None = None
        self._reference_md5: bytes | None = None
        self._read_count = 0
        self._base_count = 0
        self._block_count = 0
        self._with_tags = False   # the index carries the tags columns (M101)
        # Per-run sticky qualities strategy: block 0 auto-tunes, the
        # winner is pinned for the rest of the run.
        self._qual_hint = -1
        self._qual_exhaustive = os.environ.get("TTIO_M94Z_EXHAUSTIVE") == "1"
        # A positive TTIO_M94Z_HINT pins the strategy before block 0
        # instead, skipping the tune. V6 is reachable no other way, so
        # measuring it through the writer needs this.
        if not self._qual_exhaustive:
            raw = os.environ.get("TTIO_M94Z_HINT")
            if raw:
                try:
                    v = int(raw.strip())
                except ValueError:
                    v = 0
                if v > 0:
                    self._qual_hint = v
        self._rg = None
        self._ds: dict[str, Any] = {}
        self._idx: dict[str, Any] = {}
        self._index = None
        self._seq_layout: str | None = None
        self._embedded = False
        self._closed = False
        self._legacy_parts: list[WrittenGenomicRun] = []
        self._threads = resolve_threads(threads)
        self._window = self._threads + 1
        self._budget = resolve_memory_budget(
            memory_budget_bytes, self._threads, self._block_bytes)
        self._inflight_bytes = 0
        self._max_inflight_bytes = 0
        self._pool = None
        self._pool_ctx = None
        self._inflight: list = []
        if self._threads > 1 and not self._legacy:
            self._pool_ctx = pool_context(self._threads)
            self._pool_ctx.__enter__()
            self._pool = _cf.ThreadPoolExecutor(max_workers=self._threads,
                                                thread_name_prefix="ttio-block-encode")

    # ------------------------------------------------------------------
    @property
    def threads(self) -> int:
        return self._threads

    @property
    def max_inflight_bytes_observed(self) -> int:
        """High-water mark of the estimated bytes of in-flight
        (submitted, not yet written) blocks. The writer stalls so
        this stays at or under half the resolved byte budget
        (``memory_budget_bytes`` / ``TTIO_MEMORY_BUDGET``)."""
        return self._max_inflight_bytes

    @property
    def read_count(self) -> int:
        return self._read_count

    @property
    def block_count(self) -> int:
        return self._block_count

    def append(self, read) -> None:
        """Append one :class:`~ttio.aligned_read.AlignedRead`."""
        self.append_batch(_single_read_run(read, self._meta))

    def append_batch(self, batch: WrittenGenomicRun) -> None:
        """Append the reads of ``batch`` (its run-level metadata is
        ignored; the writer's own applies)."""
        if self._closed:
            raise RuntimeError("writer is closed")
        self._append(batch)

    def _append(self, batch: WrittenGenomicRun) -> None:
        n = int(len(batch.lengths))
        if n == 0:
            return
        if self._legacy:
            self._legacy_parts.append(batch)
            return
        if self._group and self._input_index is None:
            self._group_parts.append(batch)
            return
        # A block never spans two chromosomes: REF_DIFF_V2 encodes one
        # chromosome per blob, and a coordinate-sorted BAM streams
        # through it as long same-chromosome stretches. Unmapped reads
        # ('*') form their own blocks (BASE_PACK fallback per block).
        chroms = batch.chromosomes
        start = 0
        while start < n:
            chrom = chroms[start]
            seg_end = start + 1
            while seg_end < n and chroms[seg_end] == chrom:
                seg_end += 1
            if self._pending and self._pending_chrom != chrom:
                self._cut_block()
            self._pending_chrom = chrom
            while start < seg_end:
                room_reads = self._block_reads - self._pending_reads
                room_bytes = self._block_bytes - self._pending_bytes
                stop = min(seg_end, start + max(room_reads, 1))
                cum = np.cumsum(np.asarray(batch.lengths[start:stop], dtype=np.int64))
                fit = int(np.searchsorted(cum, room_bytes, side="right"))
                if fit < stop - start:
                    stop = start + max(fit, 1)
                part = batch if (start, stop) == (0, n) else _blocks.slice_run(batch, start, stop)
                self._pending.append(part)
                self._pending_reads += stop - start
                self._pending_bytes += int(np.asarray(part.lengths, dtype=np.int64).sum())
                if self._pending_reads >= self._block_reads or self._pending_bytes >= self._block_bytes:
                    self._cut_block()
                start = stop

    def flush(self) -> None:
        """Cut the pending reads into a block and write every block:
        after ``flush`` returns nothing is in flight. The writer cuts
        blocks on its own at block boundaries without waiting."""
        self._cut_block()
        self._drain(block_until=0)

    def _cut_block(self) -> None:
        """Encode the pending reads as one block: inline with one thread,
        else on the pool (written, in order, by a later drain)."""
        if self._legacy or not self._pending:
            return
        block = _blocks.concat_runs(self._pending)
        self._pending = []
        self._pending_reads = 0
        self._pending_bytes = 0
        if self._reference_md5 is None and self._meta["reference_chrom_seqs"] is not None:
            from .._dataset_write_genomic import _reference_md5_for_run
            probe = _apply_meta(block, self._meta, None)
            self._reference_md5 = _reference_md5_for_run(probe)
        register_block_chromosomes(block, self._chrom_map)
        block = _apply_meta(block, self._meta, self._chrom_map, self._reference_md5)
        if not self._embedded and self._meta["embed_reference"]:
            from .._dataset_write_genomic import _embed_references_for_runs
            _embed_references_for_runs(self._study, {self._name: block})
            self._embedded = True
        if self._pool is None:
            self._write_encoded(block, _blocks.encode_block(
                block, qual_strategy_hint=self._qual_hint))
            return
        if (not self._qual_exhaustive and self._qual_hint == -1
                and (self._inflight or self._block_count)):
            # Block-0 gate: the pin comes from the earliest written
            # M94Z block; hold later submissions until it is known so
            # the choice is timing-independent.
            self._drain(block_until=0)
        self._drain(block_until=self._window - 1)
        est = int(block.sequences.nbytes + block.qualities.nbytes
                  + block.offsets.nbytes * 3) * 2
        half = self._budget // 2
        while self._inflight and self._inflight_bytes + est > half:
            self._drain(block_until=len(self._inflight) - 1)
        self._inflight_bytes += est
        self._max_inflight_bytes = max(self._max_inflight_bytes, self._inflight_bytes)
        # Segment threads follow the blocks actually in flight, not
        # the pool's size; see resolve_v6_segment_threads.
        apply_v6_segment_threads(len(self._inflight) + 1)
        self._inflight.append((block, est,
                               self._pool.submit(_blocks.encode_block, block,
                                                 qual_strategy_hint=self._qual_hint)))

    def _drain(self, block_until: int) -> None:
        """Write completed blocks in sequence order; wait on the oldest
        until at most ``block_until`` remain in flight."""
        while self._inflight and (len(self._inflight) > block_until
                                  or self._inflight[0][2].done()):
            block, est, fut = self._inflight.pop(0)
            self._inflight_bytes -= est
            self._write_encoded(block, fut.result())

    def _shutdown_pool(self) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=True)
            self._pool = None
        if self._pool_ctx is not None:
            self._pool_ctx.__exit__(None, None, None)
            self._pool_ctx = None

    def _write_encoded(self, block: WrittenGenomicRun, blobs: _blocks.BlockBlobs) -> None:
        self._ensure_layout(blobs)
        if blobs.blobs.get("tags") and not self._with_tags:
            self._add_tags_columns()
        row = {"read_start": self._read_count, "n_reads": blobs.n_reads,
               "base_start": self._base_count, "n_bases": blobs.n_bases}
        for ch in _blocks.BLOCK_CHANNELS:
            if ch == "tags" and not self._with_tags:
                continue
            data = blobs.blobs[ch]
            ds = self._ds.get(ch)
            row[f"{ch}_codec"] = int(blobs.compression[ch])
            if ds is None:
                if data:
                    ds = self._create_channel(ch, blobs)
                else:
                    row[f"{ch}_off"] = 0
                    row[f"{ch}_len"] = 0
                    continue
            row[f"{ch}_off"] = int(ds.length)
            row[f"{ch}_len"] = len(data)
            if data:
                ds.append(np.frombuffer(data, dtype=np.uint8))
        self._index.append([row])
        self._append_index_arrays(block)
        self._read_count += blobs.n_reads
        self._base_count += blobs.n_bases
        self._block_count += 1
        if (self._qual_hint == -1 and not self._qual_exhaustive
                and int(blobs.compression.get("qualities", 0))
                    == int(Compression.FQZCOMP_NX16_Z)):
            from ..codecs import fqzcomp_nx16_z as _fq
            s_ = _fq.stream_strategy(blobs.blobs["qualities"])
            self._qual_hint = _fq.HINT_V4_AUTO if s_ == 4 else s_
        io.write_int_attr(self._rg, "read_count", self._read_count)
        io.write_int_attr(self._rg, "base_count", self._base_count)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self._close_inner()
        finally:
            self._shutdown_pool()

    def _close_inner(self) -> None:
        if self._legacy:
            if self._legacy_parts:
                from .._dataset_write_genomic import (
                    _embed_references_for_runs, _write_genomic_run)
                whole = _apply_meta(_blocks.concat_runs(self._legacy_parts), self._meta, None)
                if self._meta["embed_reference"]:
                    _embed_references_for_runs(self._study, {self._name: whole})
                whole = dataclasses.replace(whole, provenance_records=list(self._provenance))
                _write_genomic_run(self._runs_group(), self._name, whole)
                self._read_count = int(len(whole.lengths))
            self._legacy_parts = []
            return
        if self._group:
            parts, self._group_parts = self._group_parts, []
            whole = _blocks.concat_runs(parts) if parts else None
            order = _blocks.group_order(whole) if whole is not None else np.zeros(0, np.uint32)
            self._input_index = order
            if whole is not None:
                self._append(_blocks.take_run(whole, order))
        self.flush()
        if self._rg is None:
            self._ensure_layout(None)
        self._write_close_tables()
        if self._provenance:
            from .._dataset_write_metadata import _write_provenance
            prov = self._rg.create_group("provenance")
            _write_provenance(prov, self._provenance, dataset_name="steps")

    def __enter__(self) -> "GenomicStreamWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ------------------------------------------------------------------
    def _runs_group(self):
        if self._study.has_child("genomic_runs"):
            g = self._study.open_group("genomic_runs")
        else:
            g = self._study.create_group("genomic_runs")
            io.write_fixed_string_attr(g, "_run_names", "")
        names = [n for n in (io.read_string_attr(g, "_run_names", default="") or "").split(",") if n]
        if self._name not in names:
            names.append(self._name)
            io.write_fixed_string_attr(g, "_run_names", ",".join(names))
        return g

    def _ensure_layout(self, blobs) -> None:
        if self._rg is not None:
            return
        g = self._runs_group()
        if g.has_child(self._name):
            raise ValueError(f"genomic run {self._name!r} already exists")
        rg = g.create_group(self._name)
        m = self._meta
        io.write_int_attr(rg, "acquisition_mode", m["acquisition_mode"])
        io.write_fixed_string_attr(rg, "modality", "genomic_sequencing")
        io.write_int_attr(rg, "spectrum_class", 5)
        io.write_fixed_string_attr(rg, "reference_uri", m["reference_uri"])
        io.write_fixed_string_attr(rg, "platform", m["platform"])
        if m.get("read_role"):
            io.write_fixed_string_attr(rg, "read_role", m["read_role"])
        io.write_fixed_string_attr(rg, "sample_name", m["sample_name"])
        io.write_int_attr(rg, "read_count", 0)
        io.write_int_attr(rg, "base_count", 0)
        io.write_fixed_string_attr(rg, "layout", LAYOUT_GROUPED if self._group else LAYOUT)
        io.write_fixed_string_attr(
            rg, "block_policy", f"reads={self._block_reads},bytes={self._block_bytes}")
        # Persist non-default writer policy that shapes the coded
        # blobs, so a per-AU restore re-encodes with the same policy
        # (format-spec 9.1.1).
        if m["ref_diff_slice_bytes"]:
            io.write_int_attr(rg, "ref_diff_slice_bytes",
                              m["ref_diff_slice_bytes"], dtype="<u8")
        if m["opt_disable_qualities_v5"]:
            io.write_int_attr(rg, "opt_disable_qualities_v5", 1)
        if m["reference_chrom_seqs"] is not None:
            # The chromosome set and the reference-set digest (the
            # md5 the REF_DIFF blob headers carry and the resolver
            # verifies), so a restore can rebuild the reference from
            # /study/references or REF_PATH (format-spec 9.1.1).
            seqs = m["reference_chrom_seqs"]
            md5 = self._reference_md5
            if md5 is None:
                import hashlib
                set_md5 = getattr(seqs, "set_md5", None)
                if callable(set_md5):
                    md5 = set_md5()
                else:
                    h = hashlib.md5()
                    for c in sorted(seqs):
                        h.update(seqs[c])
                    md5 = h.digest()
                self._reference_md5 = md5
            io.write_fixed_string_attr(
                rg, "reference_md5s",
                json.dumps({c: md5.hex() for c in sorted(seqs)},
                           separators=(",", ":")))
        blocks = rg.create_group("blocks")
        self._with_tags = bool(blobs is not None and blobs.blobs.get("tags"))
        self._index = blocks.create_compound_dataset(
            "index", INDEX_FIELDS_WITH_TAGS if self._with_tags else INDEX_FIELDS,
            0, extendable=True, chunk_rows=1024)
        idx_group = rg.create_group("genomic_index")
        for name, prec, _ in _INDEX_ARRAYS + (_INPUT_INDEX if self._group else ()):
            self._idx[name] = idx_group.create_dataset(
                name, prec, 0, chunk_size=io.DEFAULT_SIGNAL_CHUNK,
                compression=Compression.ZLIB, compression_level=6, extendable=True)
        rg.create_group("signal_channels")
        self._rg = rg

    def _add_tags_columns(self) -> None:
        """The first block with tags arrives after blocks without: give
        the index the tags columns, the earlier blocks an empty range."""
        blocks = self._rg.open_group("blocks")
        rows = self._index.read_rows() if self._block_count else []
        for row in rows:
            row.update(tags_off=0, tags_len=0, tags_codec=0)
        blocks.delete_child("index")
        self._index = blocks.create_compound_dataset(
            "index", INDEX_FIELDS_WITH_TAGS, 0, extendable=True, chunk_rows=1024)
        if rows:
            self._index.append(rows)
        self._with_tags = True

    def _create_channel(self, ch: str, blobs: _blocks.BlockBlobs):
        sc = self._rg.open_group("signal_channels")
        if ch == "sequences":
            # blocks_v1 keeps every block's sequences blob, whatever its
            # codec (REF_DIFF_V2, or the BASE_PACK/rANS fallback of an
            # unmapped or reference-less block), in sequences/data; the
            # per-block codec is the sequences_codec index column.
            parent, name = sc.create_group("sequences"), "data"
        elif ch == "mate_info":
            parent, name = sc.create_group("mate_info"), "inline_v2"
        else:
            parent, name = sc, ch
        # Codec output is high-entropy: no filter. A raw channel (codec 0)
        # keeps the v1.8 zlib filter so it is not stored uncompressed.
        raw = int(blobs.compression[ch]) == 0
        ds = parent.create_dataset(name, Precision.UINT8, 0, chunk_size=CHANNEL_CHUNK,
                                   compression=(Compression.ZLIB if raw else Compression.NONE),
                                   compression_level=6, extendable=True)
        io.write_int_attr(ds, "compression", int(blobs.compression[ch]), dtype="<u1")
        for k, v in blobs.extra_attrs.get(ch, {}).items():
            ds.set_attribute(k, v)
        self._ds[ch] = ds
        return ds

    def _append_index_arrays(self, block: WrittenGenomicRun) -> None:
        # The block encoder already registered this block's chromosome
        # names in the shared map (genomic_index assigns slots to every
        # name, '*' included), so every name resolves here.
        m = self._chrom_map
        ids = np.fromiter((m[c] for c in block.chromosomes), dtype=np.uint16,
                          count=len(block.chromosomes))
        arrays = {
            "lengths": np.asarray(block.lengths, dtype=np.uint32),
            "positions": np.asarray(block.positions, dtype=np.int64),
            "mapping_qualities": np.asarray(block.mapping_qualities, dtype=np.uint8),
            "flags": np.asarray(block.flags, dtype=np.uint32),
            "chromosome_ids": ids,
        }
        for name, _, dt in _INDEX_ARRAYS:
            self._idx[name].append(arrays[name].astype(dt, copy=False))
        if self._group:
            # Blocks are written in order, so this block's rows start at
            # the reads written so far.
            r0 = self._read_count
            self._idx["input_index"].append(
                self._input_index[r0:r0 + len(block.lengths)].astype(np.uint32, copy=False))

    def _write_close_tables(self) -> None:
        names = sorted(self._chrom_map, key=self._chrom_map.get)
        rows = [{"name": n} for n in names]
        idx_group = self._rg.open_group("genomic_index")
        io.write_compound_dataset(idx_group, "chromosome_names", rows,
                                  [("name", io.vl_str())])
        sc = self._rg.open_group("signal_channels")
        if sc.has_child("mate_info"):
            mate = sc.open_group("mate_info")
        else:
            mate = sc.create_group("mate_info")
        if not mate.has_child("chrom_names"):
            io.write_compound_dataset(mate, "chrom_names", rows,
                                      [("name", io.vl_str())])


# ----------------------------------------------------------------------
def _apply_meta(run: WrittenGenomicRun, meta: dict, chrom_map,
                reference_md5: bytes | None = None) -> WrittenGenomicRun:
    return dataclasses.replace(
        run,
        reference_md5=reference_md5,
        acquisition_mode=meta["acquisition_mode"], reference_uri=meta["reference_uri"],
        platform=meta["platform"], sample_name=meta["sample_name"],
        reference_chrom_seqs=meta["reference_chrom_seqs"],
        embed_reference=meta["embed_reference"],
        opt_disable_qualities_v5=meta["opt_disable_qualities_v5"],
        signal_codec_overrides=dict(meta["signal_codec_overrides"]),
        signal_compression=meta["signal_compression"],
        chrom_name_to_id=chrom_map,
        read_role=meta.get("read_role"),
        ref_diff_slice_bytes=int(meta.get("ref_diff_slice_bytes") or 0),
    )


def _single_read_run(read, meta: dict) -> WrittenGenomicRun:
    seq = read.sequence
    seq_b = seq.encode("ascii") if isinstance(seq, str) else bytes(seq)
    qual = np.asarray(read.qualities, dtype=np.uint8)
    n = len(seq_b)
    return WrittenGenomicRun(
        acquisition_mode=meta["acquisition_mode"], reference_uri=meta["reference_uri"],
        platform=meta["platform"], sample_name=meta["sample_name"],
        positions=np.array([int(read.position)], dtype=np.int64),
        mapping_qualities=np.array([int(read.mapping_quality)], dtype=np.uint8),
        flags=np.array([int(read.flags)], dtype=np.uint32),
        sequences=np.frombuffer(seq_b, dtype=np.uint8),
        qualities=qual,
        offsets=np.zeros(1, dtype=np.uint64),
        lengths=np.array([n], dtype=np.uint32),
        cigars=[read.cigar],
        read_names=[read.read_name],
        mate_chromosomes=[read.mate_chromosome if read.mate_chromosome else "*"],
        mate_positions=np.array([int(read.mate_position)], dtype=np.int64),
        template_lengths=np.array([int(read.template_length)], dtype=np.int32),
        chromosomes=[read.chromosome],
        tags=[getattr(read, "tags", "") or ""],
    )
