"""GenomicRun — lazy view over /study/genomic_runs/<name>/.

Materialises :class:`ttio.aligned_read.AlignedRead` instances on demand
from the signal channels stored under ``signal_channels/``.  The
:class:`ttio.genomic_index.GenomicIndex` is loaded eagerly at open time
for cheap filtering and offset lookups.

Genomic analogue of :class:`ttio.acquisition_run.AcquisitionRun`.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Iterator, TYPE_CHECKING

import numpy as np

from typing import Any

from .aligned_read import AlignedRead
from .enums import AcquisitionMode, Compression, Precision
from .genomic_index import GenomicIndex
from . import _hdf5_io as io

# Hoist codec imports out of per-read accessor hot paths. Without
# these, every per-read lazy decode call invokes
# ``importlib._handle_fromlist`` — measured at ~6% of decode wall
# (~9s on chr22) when these were inside the per-read methods.
#
# v1.0 reset (Phase 2c): the v1 ``name_tokenizer`` codec was removed.
# The ULEB128 varint primitives (used by the cigars rANS schema-lift
# reader to walk the length-prefix-concat byte stream) live in the
# shared :mod:`ttio.codecs._varint` helper.
from .codecs.rans import decode as _rans_decode
from .codecs._varint import varint_decode as _varint_decode


# v1.6 (L4): _INTEGER_CHANNEL_DTYPES removed. The dict only ever
# contained positions/flags/mapping_qualities, all dropped from
# signal_channels/ in v1.6 — see docs/format-spec.md §10.7.
# The companion _int_channel_array helper has been removed below.

if TYPE_CHECKING:
    from .providers.base import StorageGroup
    from .codecs._context import CodecContext


def _wrap_hdf5_group(obj: object) -> "StorageGroup":
    """Adapt an h5py.Group to a StorageGroup; pass-through for StorageGroup."""
    from .providers.base import StorageGroup as _SG
    if isinstance(obj, _SG):
        return obj
    from .providers.hdf5 import _Group as _Hdf5Group
    return _Hdf5Group(obj)  # type: ignore[arg-type]


#: Decoded blocks a threaded sequential walk holds at most (the
#: decode-ahead window; memory is about this many block working sets).
_READ_AHEAD_BLOCKS = 4


def _read_ahead_blocks() -> int:
    """The decode-ahead window, with ``TTIO_READ_AHEAD_BLOCKS`` allowed
    to override it.

    Each block in flight stays resident until it is consumed, so this is
    a memory setting as much as a latency one, and the two do not pull
    the same way. The override exists so the trade can be measured
    rather than argued; see ``ttio.tools.genomic_read_bench``. Java:
    ``GenomicRun.readAheadBlocks``; Objective-C:
    ``TTIOReadAheadBlocks``.
    """
    raw = os.environ.get("TTIO_READ_AHEAD_BLOCKS")
    if raw:
        try:
            v = int(raw.strip())
        except ValueError:
            return _READ_AHEAD_BLOCKS
        if v > 0:
            return v
    return _READ_AHEAD_BLOCKS


@dataclass(slots=True)
class GenomicRun:
    """Lazy view over one /study/genomic_runs/<name>/ group.

    The :class:`GenomicIndex` and run-level attributes are loaded eagerly
    at :meth:`open` time so that ``len(gr)`` and region filtering are
    cheap.  Signal channel datasets are opened lazily and cached on first
    access.  Compound datasets (``cigars``, ``read_names``, ``mate_info``)
    are read whole-dataset once and cached.
    """

    name: str
    acquisition_mode: AcquisitionMode
    modality: str
    reference_uri: str
    platform: str
    sample_name: str
    index: GenomicIndex
    group: "StorageGroup"
    channel_names: list[str]  # populated for introspection / future tooling; not read by __getitem__

    # M97 — the ``@read_role`` run attribute (``hifi``, ``ont_ul``,
    # ``hic_r1``, ``hic_r2``, ``parental_maternal``,
    # ``parental_paternal``, ``illumina_wgs``, or any other stored
    # string); None when the attribute is absent. ObjC: readRole;
    # Java: getReadRole().
    read_role: str | None = None

    _signal_cache: dict = field(default_factory=dict, repr=False, compare=False)
    # Cached handle to the ``signal_channels`` child group. Opened once on
    # first per-record signal access and reused thereafter (PT2, P1.4):
    # the group is otherwise re-opened on every uncached channel access.
    # ``None`` until first use. Lifetime is the GenomicRun instance, so the
    # handle GCs with the run — same lifecycle as the eager open in open().
    _signal_group: object | None = field(default=None, repr=False, compare=False)
    _compound_cache: dict[str, list[dict]] = field(default_factory=dict, repr=False, compare=False)
    # lazy whole-channel decode cache for byte channels whose
    # @compression attribute names a TTIO codec (rANS / BASE_PACK).
    # Codec output is byte-stream non-sliceable, so the whole channel
    # is decoded once on first access and the decoded buffer is
    # sliced from memory thereafter (). Cache
    # lifetime is the GenomicRun instance — re-opening the file
    # incurs the decode cost again (Gotcha §101).
    _decoded_byte_channels: dict[str, bytes] = field(
        default_factory=dict, repr=False, compare=False,
    )
    # lazy decode cache for the read_names channel when
    # it carries a NAME_TOKENIZED codec override. Held as a
    # ``list[str]`` because the codec returns the decoded names
    # already split into per-read entries (different value type from
    # ``_decoded_byte_channels``, which holds raw concatenated
    # bytes). Per the two caches are kept
    # separate. ``None`` until first access; when populated, the
    # whole list is materialised in memory (a few hundred MB for
    # 10M reads — acceptable for typical genomic workloads, see
    # Gotcha §125).
    _decoded_read_names: list[str] | None = field(
        default=None, repr=False, compare=False,
    )
    # lazy decode cache for the cigars channel when it
    # carries a TTIO codec override (RANS_ORDER0 / RANS_ORDER1 /
    # NAME_TOKENIZED). Held as a ``list[str]`` because all three
    # codec paths produce per-read string entries — the rANS path
    # walks varint-length-prefix entries inside the decoded byte
    # buffer, and NAME_TOKENIZED returns ``list[str]`` directly.
    # Per (Option A from §2.3) this cache is
    # **separate** from ``_decoded_read_names`` — the lower-risk
    # choice that does not touch shipped Phase E code. A future
    # generalisation (Option B) could fold both into a
    # ``dict[str, list[str]]`` if a third list-of-strings channel
    # appears. Cache lifetime is the GenomicRun instance (Gotcha
    # §138).
    _decoded_cigars: list[str] | None = field(
        default=None, repr=False, compare=False,
    )
    # M101: decoded SAM_TAGS channel, one string per read; [] when the
    # run has no tags channel. None until first access.
    _decoded_tags: list[str] | None = field(
        default=None, repr=False, compare=False,
    )
    # combined per-field cache for the mate_info subgroup
    # layout (Gotcha §144). Held as a single
    # ``dict[str, Any]`` keyed by field name (``"chrom"`` →
    # list[str]; ``"pos"`` → np.ndarray int64; ``"tlen"`` →
    # np.ndarray int32) because the three fields have three different
    # value types — a typed-per-field cache would force a union or
    # three separate fields. Used only for the Phase F subgroup
    # layout; the M82 compound layout still uses the existing
    # ``_compound_cache`` via :meth:`_compound`. ``_mate_<field>_at``
    # populates the corresponding key on first access.
    _decoded_mate_info: dict[str, Any] = field(
        default_factory=dict, repr=False, compare=False,
    )

    # Cached result of `_mate_info_is_subgroup()`. Without this, the
    # method does an HDF5 link-type probe on every call, and the
    # per-read decode path calls it 3x per read — at chr22 scale
    # (1.77M reads) that's 5.3M redundant probes resolving to the
    # same answer, dominating decode wall-time. None = not yet
    # computed; True/False = cached result.
    _mate_info_subgroup_cached: "bool | None" = field(
        default=None, repr=False, compare=False,
    )

    # cached whole-sequence decode from the refdiff_v2 blob.
    # None = not yet probed; b"" = probed and found to be v1/BASE_PACK;
    # non-empty bytes = decoded concatenated sequence bytes (total_bases long).
    # Populated on first access by _sequences_is_ref_diff_v2() +
    # _decode_ref_diff_v2_sequences(). Cache lifetime = GenomicRun instance.
    _decoded_ref_diff_v2: "bytes | None" = field(
        default=None, repr=False, compare=False,
    )
    # None = not yet probed; True/False = cached probe result.
    _sequences_is_v2_cached: "bool | None" = field(
        default=None, repr=False, compare=False,
    )

    # Lazily-built run-derived CodecContext shared across all registry
    # decode calls (codec-registry refactor). None until first
    # _codec_context() call; immutable thereafter for the run lifetime.
    _codec_ctx_cache: "CodecContext | None" = field(
        default=None, repr=False, compare=False,
    )

    # The /study/references group (StorageGroup) threaded in by
    # SpectralDataset._from_provider so REF_DIFF decode can resolve
    # embedded reference chromosomes through the storage protocol rather
    # than a raw h5py handle (P3.9). None when the file carries no
    # embedded references; resolve() then falls back to REF_PATH or
    # raises RefMissingError.
    _references_group: "StorageGroup | None" = field(
        default=None, repr=False, compare=False,
    )

    # When True (local files), :meth:`_byte_channel_slice` bulk-reads an
    # uncompressed byte channel once and slices in memory. When False
    # (remote fsspec-backed runs), it keeps the per-record hyperslab read
    # so a single random access does not pull the whole column over the
    # network. Codec-compressed channels are always whole-channel decoded
    # (the codec output is not sliceable) regardless of this flag. Set by
    # SpectralDataset._from_provider. The returned bytes are identical
    # either way.
    _bulk_read: bool = field(default=True, repr=False, compare=False)
    # blocks_v1 (format-spec 10.12): the block table, the last
    # materialised block view (block number, GenomicRun over the view)
    # and a counter for tests. ``_layout`` is "whole" for v1.8 files.
    _layout: str = field(default="whole", repr=False, compare=False)
    _block_table: "Any" = field(default=None, repr=False, compare=False)
    _block_cache: "tuple[int, GenomicRun] | None" = field(default=None, repr=False, compare=False)
    _blocks_materialised: int = field(default=0, repr=False, compare=False)
    _name_tables: "tuple[list, list | None] | None" = field(default=None, repr=False, compare=False)

    # ------------------------------------------------------------------
    # Sequence protocol
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self.index.count

    def __iter__(self) -> Iterator[AlignedRead]:
        return self.iter_reads()

    @property
    def layout(self) -> str:
        """``"blocks_v1"`` or ``"whole"`` (the v1.8 whole-channel layout)."""
        return self._layout

    @property
    def block_count(self) -> int:
        """Number of blocks (1 for a whole-channel run)."""
        return self._block_table.count if self._block_table is not None else 1

    def iter_reads(self, start: int = 0, stop: int | None = None, *,
                   threads: int | None = None) -> Iterator[AlignedRead]:
        """Yield reads ``[start, stop)`` in order. Under ``blocks_v1``
        the next few blocks decode ahead on a pool, at most
        ``_READ_AHEAD_BLOCKS`` of them; one thread keeps the one-block
        path, holding one decoded block at a time. Each block in
        flight stays resident, so the window costs memory as much as
        it saves latency."""
        n = len(self)
        if stop is None or stop > n:
            stop = n
        if start < 0:
            start += n
        if self._layout != "blocks_v1":
            for i in range(max(start, 0), stop):
                yield self[i]
            return
        from ._threads import resolve_threads, pool_context
        t = self._block_table
        nthreads = resolve_threads(threads)
        i = max(start, 0)
        if nthreads <= 1 or i >= stop:
            while i < stop:
                b = t.block_for(i)
                r0 = int(t.read_start[b])
                b_end = min(r0 + int(t.n_reads[b]), stop)
                view = self._block_view(b)
                for j in range(i, b_end):
                    yield view[j - r0]
                i = b_end
            return
        import concurrent.futures as _cf
        b_first, b_last = t.block_for(i), t.block_for(stop - 1)
        # The serial consumer only needs enough decode-ahead to never
        # stall; more in flight is pure memory, since the window is how
        # many decoded blocks stay resident.
        window = min(nthreads, _read_ahead_blocks())
        # pool_context sizes V6 segment threads from the number of
        # workers, so it takes the window and not nthreads: only this
        # many blocks decode at once, and the rest of the machine is
        # what the segments of each one are free to use.
        with pool_context(window), _cf.ThreadPoolExecutor(
                max_workers=window, thread_name_prefix="ttio-block-decode") as pool:
            pending: dict[int, "_cf.Future"] = {}

            def submit(b: int) -> None:
                if b <= b_last and b not in pending:
                    pending[b] = pool.submit(self._warm, self._prefetch_view(b))

            for b in range(b_first, min(b_last, b_first + window - 1) + 1):
                submit(b)
            b = b_first
            while i < stop:
                view = pending.pop(b).result()
                submit(b + window)
                r0 = int(t.read_start[b])
                b_end = min(r0 + int(t.n_reads[b]), stop)
                for j in range(i, b_end):
                    yield view[j - r0]
                i = b_end
                b += 1

    def _block_window(self, nthreads: int, b_first: int, b_last: int) -> int:
        """How many blocks may be in flight at once, given how many
        threads want one each.

        A block in flight is resident for as long as the caller is
        inside it, so this is a memory setting before it is a
        concurrency one, and the thread count is only an upper bound. A
        block is charged 8 times its sequence bytes, the figure
        ``resolve_memory_budget`` uses on the writing side, so
        ``TTIO_MEMORY_BUDGET`` and the half-of-physical cap mean the
        same thing to a reader as to a writer. The widest block in the
        range sets the size, not the mean.

        Java: ``GenomicRun.blockWindow``; Objective-C:
        ``-_blockWindowForThreads:first:last:``.
        """
        from ._threads import resolve_memory_budget
        t = self._block_table
        widest = max(int(t.n_bases[b]) for b in range(b_first, b_last + 1))
        if widest == 0:
            return 1
        budget = resolve_memory_budget(None, nthreads, widest)
        admits = max(1, budget // (widest * 8))
        blocks = b_last - b_first + 1
        return max(1, min(nthreads, admits, blocks))

    def for_each_block(self, fn, start: int = 0, stop: int | None = None, *,
                       threads: int | None = None) -> None:
        """Call ``fn(view, view_start, first_read, n_reads)`` once per
        decoded block, on the pool.

        Read record ``k`` of the delivered range as
        ``view[view_start + k]``; its index in the whole run is
        ``first_read + k``. The two differ whenever a range starts
        part-way into a block, so a caller that indexes the view by
        ``first_read`` reads the wrong records.

        ``fn`` is called from several threads at once and in no
        particular order, and must be safe to call that way; ``view`` is
        a run over one block's reads, indexed from 0, valid only for the
        duration of the call. That relaxed ordering is the whole
        difference from :meth:`iter_reads`, whose in-order delivery on
        the calling thread is what bounds it.

        ⚠️ **In CPython this buys nothing, and past a few workers it
        costs.** The callback holds the GIL, so the consumers serialise
        against each other while still paying for the contention;
        measured over 1.5M reads: serial ``iter_reads`` 6.2 MB/s, this
        method 6.1 at one worker, 6.8 at four, and **4.5 at eight and
        above** — slower than not threading at all. Only the decode
        overlaps, because it releases the GIL inside the native call.
        The compiled SDKs behave the other way round, which is why the
        default is capped here and not there.

        So an unset ``threads`` resolves to at most
        ``_read_ahead_blocks()`` workers rather than to ``TTIO_THREADS``:
        the default has to be the best measured setting, not the worst.
        An explicit ``threads`` is honoured as given. The same
        divergence already applies to the decode-ahead window, which is
        a per-SDK setting for the same reason.

        Java: ``GenomicRun.iterBlocks``; Objective-C:
        ``-iterBlocksFrom:to:threads:error:usingBlock:``.
        """
        n = len(self)
        if stop is None or stop > n:
            stop = n
        lo = max(start + n if start < 0 else start, 0)
        if lo >= stop:
            return
        if self._layout != "blocks_v1":
            fn(self, lo, lo, stop - lo)
            return
        import concurrent.futures as _cf
        from ._threads import resolve_threads, pool_context
        t = self._block_table
        nthreads = (resolve_threads(threads) if threads is not None
                    else min(resolve_threads(None), _read_ahead_blocks()))
        b_first, b_last = t.block_for(lo), t.block_for(stop - 1)
        window = self._block_window(nthreads, b_first, b_last)
        blocks = b_last - b_first + 1
        ahead_n = min(blocks, max(2, _read_ahead_blocks()))

        # Below the decode-ahead depth, one consumer fed by that many
        # decoders beats the same number of threads each decoding and
        # then consuming its own block: decoding is the slower half, so
        # until there are threads enough to cover it, pipelining wins
        # over dividing.
        if window < ahead_n:
            self._for_each_block_pipelined(
                fn, lo, stop, b_first, b_last, min(ahead_n, max(window, 1) * 4))
            return

        def one(b: int, view: "GenomicRun") -> None:
            view = self._warm(view)
            r0 = int(t.read_start[b])
            frm, to = max(r0, lo), min(r0 + int(t.n_reads[b]), stop)
            if to > frm:
                fn(view, frm - r0, frm, to - frm)

        with pool_context(window), _cf.ThreadPoolExecutor(
                max_workers=window, thread_name_prefix="ttio-block") as pool:
            pending: dict[int, "_cf.Future"] = {}

            def submit(b: int) -> None:
                # _prefetch_view runs here, on the calling thread: those
                # are storage reads, and only the decode is safe to hand
                # to the pool.
                if b <= b_last and b not in pending:
                    pending[b] = pool.submit(one, b, self._prefetch_view(b))

            for b in range(b_first, min(b_last, b_first + window - 1) + 1):
                submit(b)
            for b in range(b_first, b_last + 1):
                f = pending.pop(b, None)
                if f is None:
                    continue
                f.result()
                submit(b + window)

    def _for_each_block_pipelined(self, fn, lo: int, stop: int, b_first: int,
                                  b_last: int, ahead_n: int) -> None:
        """One consumer on this thread, ``ahead_n`` decoders behind it.
        The view is built lazily, so :meth:`_warm` has to run on the
        pool: without it the channel decode lands on the consumer's
        thread at its first read and the overlap is lost entirely."""
        import concurrent.futures as _cf
        from ._threads import pool_context
        t = self._block_table
        with pool_context(ahead_n), _cf.ThreadPoolExecutor(
                max_workers=max(2, ahead_n),
                thread_name_prefix="ttio-block-decode") as pool:
            pending: dict[int, "_cf.Future"] = {}

            def submit(b: int) -> None:
                if b <= b_last and b not in pending:
                    pending[b] = pool.submit(self._warm, self._prefetch_view(b))

            for b in range(b_first, min(b_last, b_first + ahead_n - 1) + 1):
                submit(b)
            for b in range(b_first, b_last + 1):
                view = pending.pop(b).result()
                submit(b + ahead_n)
                r0 = int(t.read_start[b])
                frm, to = max(r0, lo), min(r0 + int(t.n_reads[b]), stop)
                if to > frm:
                    fn(view, frm - r0, frm, to - frm)

    def _prefetch_view(self, b: int) -> "GenomicRun":
        """The block view for ``b``, built on the caller's thread (storage
        reads) but not yet decoded; :meth:`_warm` decodes it and is safe
        to run on a pool thread."""
        from .genomic._block_view import materialise_block, _rows_of
        if self._name_tables is None:
            idx = self.group.open_group("genomic_index")
            sc = self.group.open_group("signal_channels")
            mate = _rows_of(sc.open_group("mate_info"), "chrom_names") if sc.has_child("mate_info") else None
            self._name_tables = (_rows_of(idx, "chromosome_names"), mate)
        grp = materialise_block(self.group, self._block_table, b,
                                chrom_name_rows=self._name_tables[0],
                                mate_chrom_rows=self._name_tables[1])
        self._blocks_materialised += 1
        return GenomicRun.open(grp, self.name, references_group=self._references_group,
                               bulk_read=self._bulk_read)

    @staticmethod
    def _warm(view: "GenomicRun") -> "GenomicRun":
        if len(view):
            view[0]          # decodes every channel of the block into its caches
        return view

    def _block_view(self, b: int) -> "GenomicRun":
        if self._block_cache is not None and self._block_cache[0] == b:
            return self._block_cache[1]
        sub = self._prefetch_view(b)
        self._block_cache = (b, sub)
        return sub

    def provenance_chain(self) -> "list":
        """Return per-run provenance records in insertion order.

        Closes the M91 read-side gap: until Phase 1 the lazy
        ``GenomicRun`` view didn't expose provenance, so cross-
        modality queries had to fall back to the ``sample_name``
        attribute. Now both run types share the same accessor.

        Reads from the ``<run>/provenance/steps`` compound dataset
        written by :func:`spectral_dataset._write_genomic_run`.
        Returns ``[]`` for runs that carry no provenance.
        """
        # Navigate the compound-provenance subgroup through the
        # StorageGroup protocol. The provenance compound layout is
        # identical to the MS path (see
        # acquisition_run.AcquisitionRun.provenance) and is decoded by
        # the same helper.
        from .acquisition_run import _decode_provenance_compound
        if self.group.has_child("provenance"):
            prov = self.group.open_group("provenance")
            if prov.has_child("steps"):
                return _decode_provenance_compound(prov, "steps")
        return []

    def __getitem__(self, i: int) -> AlignedRead:
        if i < 0:
            i += len(self)
        if not 0 <= i < len(self):
            raise IndexError(
                f"read index {i} out of range [0, {len(self)})"
            )

        if self._layout == "blocks_v1":
            b = self._block_table.block_for(i)
            return self._block_view(b)[i - int(self._block_table.read_start[b])]
        offset = int(self.index.offsets[i])
        length = int(self.index.lengths[i])

        # Per-read scalar fields come straight from the index.
        position = int(self.index.positions[i])
        mapq = int(self.index.mapping_qualities[i])
        flag = int(self.index.flags[i])
        chrom = self.index.chromosomes[i]

        # Sequence and qualities — read a slice of the per-base channels.
        # routed through _byte_channel_slice so that channels
        # written with a TTIO codec override (@compression > 0) are
        # decoded transparently before slicing.
        seq_bytes = self._byte_channel_slice("sequences", offset, length)
        sequence = seq_bytes.decode("ascii")
        qualities = self._byte_channel_slice("qualities", offset, length)

        # Compound / codec-lifted channels — dispatch on dataset
        # shape (Phases C and E). The compound path delegates
        # to ``_compound`` (whole-dataset cached), the codec path
        # decodes once and caches ``list[str]``.
        cigar = self._cigar_at(i)

        read_name = self._read_name_at(i)

        # dispatch on HDF5 link type (compound dataset =
        # M82 path; subgroup = Phase F per-field path). The three
        # helpers each open the bare ``mate_info`` link, detect the
        # layout, and route to either the existing ``_compound``
        # cache (M82) or the per-field codec/natural-dtype dispatch
        # (Phase F).
        mate_chromosome = self._mate_chrom_at(i)
        mate_position = self._mate_pos_at(i)
        template_length = self._mate_tlen_at(i)
        tags = self._all_tags()

        return AlignedRead(
            read_name=read_name,
            chromosome=chrom,
            position=position,
            mapping_quality=mapq,
            cigar=cigar,
            sequence=sequence,
            qualities=qualities,
            flags=flag,
            mate_chromosome=mate_chromosome,
            mate_position=mate_position,
            template_length=template_length,
            tags=tags[i] if tags else "",
        )

    # ------------------------------------------------------------------
    # Region query
    # ------------------------------------------------------------------

    def reads_in_region(
        self, chromosome: str, start: int, end: int
    ) -> list[AlignedRead]:
        """Return reads on ``chromosome`` whose mapping position is in ``[start, end)``.

        Note: filters by mapping position only, not by read end coordinate.
        A read whose start lies outside the window but extends into it
        will NOT be returned. Use SAM-style overlap semantics in a future
        enhancement if needed.
        """
        return [
            self[i]
            for i in self.index.indices_for_region(chromosome, start, end)
        ]

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def open(cls, group, name: str, *, references_group=None,
             bulk_read: bool = True) -> "GenomicRun":
        """Open an existing genomic_runs/<name>/ group.

        Mirrors :meth:`ttio.acquisition_run.AcquisitionRun.open`: the
        caller resolves the child group before calling this classmethod.
        The genomic index and run-level attributes are loaded eagerly;
        signal channel datasets remain closed until first access.

        Args:
            references_group: the ``/study/references`` group as a
                :class:`~ttio.providers.base.StorageGroup`, threaded
                through to the REF_DIFF :class:`ReferenceResolver`. May
                be ``None`` when the file carries no embedded references.
        """

        sgroup = _wrap_hdf5_group(group)

        layout = io.read_string_attr(sgroup, "layout") or "whole"
        block_table = None
        idx_group = sgroup.open_group("genomic_index")
        if layout == "blocks_v1":
            # blocks_v1: the block table is small; the per-read index
            # arrays load lazily on first use (format-spec 10.12).
            from .genomic._block_view import BlockTable, LazyGenomicIndex
            block_table = BlockTable.read(sgroup)
            index = LazyGenomicIndex(idx_group, block_table)
        elif layout != "whole":
            raise ValueError(
                f"genomic run {name!r}: unsupported layout {layout!r} "
                "(this reader knows the whole-channel layout and blocks_v1)")
        else:
            # Eager: load the genomic index.
            index = GenomicIndex.read(idx_group)

        # Eager: list signal channel names.
        sig = sgroup.open_group("signal_channels")
        channel_names = list(sig.child_names())

        # Eager: run-level attributes written by _write_genomic_run.
        acq_mode_raw = io.read_int_attr(sgroup, "acquisition_mode")
        modality = io.read_string_attr(sgroup, "modality") or ""
        reference_uri = io.read_string_attr(sgroup, "reference_uri") or ""
        platform = io.read_string_attr(sgroup, "platform") or ""
        sample_name = io.read_string_attr(sgroup, "sample_name") or ""
        # @read_role is absent on pre-M97 files; None, not "".
        read_role = io.read_string_attr(sgroup, "read_role") or None

        return cls(
            name=name,
            acquisition_mode=AcquisitionMode(int(acq_mode_raw)),
            modality=modality,
            reference_uri=reference_uri,
            platform=platform,
            sample_name=sample_name,
            read_role=read_role,
            index=index,
            group=sgroup,
            channel_names=channel_names,
            _references_group=references_group,
            _bulk_read=bulk_read,
            _layout=layout,
            _block_table=block_table,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _signal_channels_group(self):
        """Open the signal_channels group once and cache the handle."""
        if self._signal_group is None:
            self._signal_group = self.group.open_group("signal_channels")
        return self._signal_group

    def _signal_dataset(self, name: str):
        """Open a primitive signal-channel dataset and cache the handle."""
        if name not in self._signal_cache:
            sig = self._signal_channels_group()
            self._signal_cache[name] = sig.open_dataset(name)
        return self._signal_cache[name]

    def _codec_context(self) -> "CodecContext":
        """Build (and cache) the run-derived CodecContext for registry decode.

        Every registry decode call routed from this run shares one
        immutable context. ``own_chrom_ids`` is sourced with the same
        encounter-order uint16 assignment that ``_decode_mate_inline_v2``
        uses, so the MATE_INLINE_V2 decode adapter sees identical ids.
        """
        from .codecs._context import CodecContext
        if self._codec_ctx_cache is not None:
            return self._codec_ctx_cache
        import numpy as np
        idx = self.index
        flags = np.asarray(idx.flags, dtype=np.uint32)
        revcomp = ((flags & 16) != 0).astype(np.uint8)  # vectorized

        # own_chrom_ids (uint16): encounter-order id assignment over the
        # index chromosomes — verbatim mirror of _decode_mate_inline_v2's
        # writer-matching derivation so MATE_INLINE_V2 decode is identical.
        n = int(idx.count)
        # The writer's id assignment is the mate_info/chrom_names table
        # (row index = id): encounter order over own chromosomes plus
        # mate-only names for a whole-channel run, and the run-wide map
        # for a blocks_v1 run (whose per-block encounter order differs).
        # Seed from that table when present; otherwise rebuild encounter
        # order (files without mate_info never need own ids).
        name_to_id: dict[str, int] = {}
        try:
            sig = self._signal_channels_group()
            if sig.has_child("mate_info"):
                mate_group = sig.open_group("mate_info")
                if mate_group.has_child("chrom_names"):
                    for row_i, row in enumerate(io.read_compound_dataset(mate_group, "chrom_names")):
                        v = row["name"]
                        name_to_id[v.decode("utf-8") if isinstance(v, bytes) else v] = row_i
        except KeyError:
            name_to_id = {}
        own_chrom_ids = np.empty(n, dtype=np.uint16)
        for i, cname in enumerate(idx.chromosomes):
            if cname == "*" or not cname:
                own_chrom_ids[i] = 0xFFFF
            else:
                if cname not in name_to_id:
                    name_to_id[cname] = len(name_to_id)
                own_chrom_ids[i] = name_to_id[cname]

        from .genomic.reference_resolver import ReferenceResolver
        # ReferenceResolver navigates the /study/references group through
        # the StorageGroup protocol (P3.9). references_group is None when
        # the file has no embedded refs; resolve() then uses REF_PATH or
        # raises RefMissingError — the same clear failure the old
        # non-HDF5 path produced.
        resolver = ReferenceResolver(self._references_group)
        ctx = CodecContext(
            read_lengths=np.asarray(idx.lengths, dtype=np.uint32),
            revcomp_flags=revcomp,
            read_count=int(idx.count),
            positions=np.asarray(idx.positions, dtype=np.int64),
            cigars_provider=self._all_cigars,
            total_bases=int(sum(idx.lengths)),
            chromosomes=list(idx.chromosomes),
            own_positions=np.asarray(idx.positions, dtype=np.int64),
            n_records=int(idx.count),
            own_chrom_ids=own_chrom_ids,
            reference_resolver=resolver,
            # fqzcomp V5: the qualities decoder needs the decoded
            # sequences channel; route through the decode-once byte
            # cache so ref_diff_v2 and plain layouts both work and the
            # provider fires only for version-5 streams.
            sequences_provider=lambda: self._byte_channel_slice(
                "sequences", 0, int(sum(idx.lengths))),
            tag_references_provider=lambda: self._tag_references(name_to_id, resolver),
        )
        self._codec_ctx_cache = ctx
        return ctx

    def _tag_references(self, name_to_id: dict[str, int], resolver) -> "list[bytes | None]":
        """References for SAM_TAGS MD/NM derivation, indexed like
        ``own_chrom_ids``. The writer derives only on REF_DIFF_V2 runs
        (binding decision §96), whose one chromosome resolves the same
        way the sequences decode does; any other run gets none."""
        if not self._sequences_is_ref_diff_v2():
            return []
        from .codecs import ref_diff_v2
        ds = self._signal_channels_group().open_group("sequences").open_dataset("refdiff_v2")
        # The outer header is 38 bytes plus a uint16-length URI.
        header = ref_diff_v2.parse_blob_header(
            bytes(ds.read(offset=0, count=min(int(ds.length), 38 + 0xFFFF))))
        names = sorted(set(self.index.chromosomes) - {"*", ""})
        refs: list[bytes | None] = [None] * len(name_to_id)
        if len(names) == 1 and names[0] in name_to_id:
            refs[name_to_id[names[0]]] = resolver.resolve(
                uri=header.reference_uri, expected_md5=header.reference_md5, chromosome=names[0])
        return refs

    def _all_tags(self) -> list[str]:
        """Every read's SAM tag text (M101), decoded once; [] when the
        run has no ``tags`` channel."""
        if self._decoded_tags is not None:
            return self._decoded_tags
        sig = self._signal_channels_group()
        if not sig.has_child("tags"):
            self._decoded_tags = []
            return self._decoded_tags
        ds = sig.open_dataset("tags")
        codec_id = io.read_int_attr(ds, "compression", default=0) or 0
        if codec_id != int(Compression.SAM_TAGS):
            raise ValueError(
                f"signal_channel 'tags': @compression={codec_id} is not "
                "SAM_TAGS (18), the only codec of the tags channel")
        from .codecs._registry import CODEC_REGISTRY
        from .codecs._context import ChannelPayload
        blob = bytes(ds.read(offset=0, count=int(ds.length)))
        self._decoded_tags = CODEC_REGISTRY[Compression.SAM_TAGS].decode(
            ChannelPayload.of_bytes(blob), self._codec_context()).as_str_list()
        return self._decoded_tags

    def _byte_channel_slice(self, name: str, offset: int, count: int) -> bytes:
        """Return bytes ``[offset, offset+count)`` for a uint8 byte channel.

        M86 dispatch: for codec-compressed channels (``@compression > 0``)
        the whole channel is decoded once on first access, the decoded
        buffer is cached on this :class:`GenomicRun` instance, and
        subsequent slices are taken from the cached bytes. For
        uncompressed channels (no attribute or value 0) the existing
        per-slice :meth:`StorageDataset.read` path is used unchanged
        — no whole-channel materialisation, no behaviour change vs M82.

        v1.8 extension: when ``name == "sequences"`` and the link at
        signal_channels/sequences is a GROUP (refdiff_v2 layout), the
        v2 decode path is used and the result is stored in
        ``_decoded_ref_diff_v2`` / ``_decoded_byte_channels["sequences"]``.
        """
        cached = self._decoded_byte_channels.get(name)
        if cached is not None:
            return cached[offset:offset + count]

        # v1.8 probe: for sequences, check for the group layout first.
        if name == "sequences" and self._sequences_is_ref_diff_v2():
            from .codecs._registry import CODEC_REGISTRY
            from .codecs._context import ChannelPayload
            sig = self._signal_channels_group()
            decoded = CODEC_REGISTRY[Compression.REF_DIFF_V2].decode(
                ChannelPayload.of_group(sig.open_group("sequences")),
                self._codec_context(),
            ).as_bytes()
            self._decoded_byte_channels[name] = decoded
            self._decoded_ref_diff_v2 = decoded
            return decoded[offset:offset + count]

        ds = self._signal_dataset(name)
        # The @compression attr is probed only on a cache MISS — once the
        # channel is cached below, the top-of-function ``cached`` check
        # returns the slice without re-reading the attr, so this is a
        # one-time-per-channel determination (was previously re-read on
        # every slice for the uncompressed path).
        codec_id = io.read_int_attr(ds, "compression", default=0) or 0
        if codec_id == 0:
            if not self._bulk_read:
                # Remote (fsspec) path: keep the per-record hyperslab read
                # so a single random access only pulls its own slice over
                # the network, not the whole column. Byte-identical to the
                # bulk-then-slice path.
                return bytes(ds.read(offset=offset, count=count))
            # Uncompressed: bulk-read the WHOLE byte channel once, cache,
            # and slice in memory thereafter — mirroring the compressed
            # branch below and the ObjC/Java fast path. The cached slice
            # ``[offset:offset+count]`` is byte-identical to the old
            # per-call ``ds.read(offset=offset, count=count)``.
            all_bytes = bytes(ds.read(offset=0, count=int(ds.length)))
            self._decoded_byte_channels[name] = all_bytes
            return all_bytes[offset:offset + count]

        # Compressed: read all bytes, decode, cache for subsequent slices.
        all_bytes = bytes(ds.read(offset=0, count=int(ds.length)))
        from .codecs._registry import CODEC_REGISTRY
        from .codecs._context import ChannelPayload
        try:
            codec = CODEC_REGISTRY[Compression(codec_id)]
        except (KeyError, ValueError):
            raise ValueError(
                f"signal_channel '{name}': @compression={codec_id} "
                "is not a supported TTIO codec id")
        decoded = codec.decode(ChannelPayload.of_bytes(all_bytes), self._codec_context()).as_bytes()
        self._decoded_byte_channels[name] = decoded
        return decoded[offset:offset + count]

    def _sequences_is_ref_diff_v2(self) -> bool:
        """True iff signal_channels/sequences is a GROUP containing refdiff_v2.

        v1.8+ layout: sequences is a group with a refdiff_v2 child dataset
        (@compression = 14). Result is cached via ``_sequences_is_v2_cached``.
        """
        if self._sequences_is_v2_cached is not None:
            return self._sequences_is_v2_cached
        sig = self._signal_channels_group()
        try:
            seq_grp = sig.open_group("sequences")
            # It's a group — check for the refdiff_v2 child dataset.
            try:
                seq_grp.open_dataset("refdiff_v2")
                result = True
            except KeyError:
                result = False
        except KeyError:
            # sequences is a dataset (v1) — not v2.
            result = False
        self._sequences_is_v2_cached = result
        return result

    def _all_cigars(self) -> list[str]:
        """Return the full list of CIGAR strings for this run.

        Honours the M86 Phase C codec dispatch on the ``cigars``
        channel (RANS / NAME_TOKENIZED override → uint8 dataset; no
        override → M82 compound dataset). Caches the result on
        ``self._decoded_cigars`` so subsequent ``_cigar_at`` calls
        hit the cache.
        """
        if self._decoded_cigars is not None:
            return self._decoded_cigars
        # Trigger the existing per-read decode path once; the helper
        # populates ``self._decoded_cigars`` for the rANS / tokenised
        # paths. The compound (M82) path doesn't populate the field,
        # so fall back to a manual walk in that case.
        _ = self._cigar_at(0) if len(self) > 0 else None
        if self._decoded_cigars is not None:
            return self._decoded_cigars
        # M82 compound layout — gather from the compound cache.
        cigars = self._compound("cigars")
        out = [c["value"] for c in cigars]
        self._decoded_cigars = out
        return out

    # v1.6 (L4): _int_channel_array removed. The helper supported
    # reading positions/flags/mapping_qualities from signal_channels/
    # via codec dispatch — but the per-record reader path (__getitem__)
    # always used genomic_index/ via self.index.{positions,...}. With
    # v1.6 dropping the signal_channels/ duplicates, no caller has a
    # reason to use this helper. See docs/format-spec.md §10.7.

    def _compound(self, name: str) -> list[dict]:
        """Read a compound dataset whole and cache it.

        ``read_compound_dataset`` already decodes VL bytes to ``str``, so
        callers never need to check ``isinstance(v, bytes)``.
        """
        if name not in self._compound_cache:
            sig = self._signal_channels_group()
            self._compound_cache[name] = io.read_compound_dataset(sig, name)
        return self._compound_cache[name]

    def _read_name_at(self, i: int) -> str:
        """Return the read name at index ``i``.

        v1.0 reset (Phase 2c): read_names is always a flat 1-D uint8
        dataset encoded with NAME_TOKENIZED_V2 (codec id 15). The v1
        NAME_TOKENIZED (codec id 8) decoder and the M82 compound
        fallback were both removed; readers reject the v1 layout
        with a clear migration error. Decoded names are cached as a
        ``list[str]`` on this :class:`GenomicRun` instance.
        """
        cached = self._decoded_read_names
        if cached is not None:
            return cached[i]

        sig = self._signal_channels_group()
        ds = sig.open_dataset("read_names")

        if ds.precision != Precision.UINT8:
            raise ValueError(
                "signal_channel 'read_names': dataset is not a flat "
                "uint8 codec stream — the v1.0 reader requires the "
                "NAME_TOKENIZED_V2 layout. The legacy M82 compound "
                "VL-string layout is no longer supported; this file "
                "was written with an older TTI-O version. Re-encode "
                "with v1.0+."
            )

        codec_id = io.read_int_attr(ds, "compression", default=0) or 0
        if codec_id == int(Compression.NAME_TOKENIZED_V2):
            from .codecs._registry import CODEC_REGISTRY
            from .codecs._context import ChannelPayload
            all_bytes = bytes(ds.read(offset=0, count=int(ds.length)))
            self._decoded_read_names = CODEC_REGISTRY[
                Compression.NAME_TOKENIZED_V2
            ].decode(
                ChannelPayload.of_bytes(all_bytes), self._codec_context()
            ).as_str_list()
            return self._decoded_read_names[i]
        raise ValueError(
            f"signal_channel 'read_names': @compression={codec_id} "
            "is not a supported TTIO codec id for the read_names "
            "channel (only NAME_TOKENIZED_V2 = 15 is recognised)"
        )

    def _cigar_at(self, i: int) -> str:
        """Return the cigar string at index ``i``, dispatching on shape.

        M86 Phase C: cigars has two on-disk layouts (Binding
        Decisions §120-§123):

        - **M82 compound** (no override): VL_STRING-in-compound
          dataset, read whole-and-cache via :meth:`_compound`.
        - **rANS codec** (override active): flat 1-D uint8 dataset
          carrying a length-prefix-concat byte stream
          (``varint(len) + bytes`` per CIGAR — §2.5 of the Phase C
          plan / format-spec §10.6 extended). Decoded once on first
          access and cached as a ``list[str]`` per Binding Decision
          §123. Two codec ids are recognised: ``RANS_ORDER0`` (4)
          and ``RANS_ORDER1`` (5).

        Dispatch is on dataset shape — a 1-D uint8 dataset routes
        through the codec path; anything else (compound) falls
        through to the M82 path. The :attr:`_decoded_cigars` cache
        holds the entire decoded list across calls.
        """
        cached = self._decoded_cigars
        if cached is not None:
            return cached[i]

        sig = self._signal_channels_group()
        ds = sig.open_dataset("cigars")

        if ds.precision == Precision.UINT8:
            codec_id = io.read_int_attr(ds, "compression", default=0) or 0
            all_bytes = bytes(ds.read(offset=0, count=int(ds.length)))
            if codec_id in (
                int(Compression.RANS_ORDER0),
                int(Compression.RANS_ORDER1),
            ):
                decoded = _rans_decode(all_bytes)
                _vd = _varint_decode
                # Walk the length-prefix-concat byte stream — the
                # mirror of the writer's serialisation contract
                # (§2.5). Each entry is varint(len) + len bytes
                # of ASCII payload.
                out: list[str] = []
                offset = 0
                n = len(decoded)
                while offset < n:
                    length, offset = _vd(decoded, offset)
                    if offset + length > n:
                        raise ValueError(
                            "cigars rANS stream: length-prefix-concat "
                            f"entry runs off end of decoded buffer "
                            f"(offset={offset}, length={length}, "
                            f"buffer_size={n})"
                        )
                    payload = decoded[offset:offset + length]
                    offset += length
                    try:
                        out.append(payload.decode("ascii"))
                    except UnicodeDecodeError as exc:
                        raise ValueError(
                            "cigars rANS stream: entry contains "
                            "non-ASCII bytes"
                        ) from exc
                self._decoded_cigars = out
                return out[i]
            raise ValueError(
                f"signal_channel 'cigars': @compression={codec_id} "
                "is not a supported TTIO codec id for the cigars "
                "channel (only RANS_ORDER0 = 4 and RANS_ORDER1 = 5 "
                "are recognised)"
            )

        # Compound path (M82, no override). Materialise the whole
        # list on first call and cache in self._decoded_cigars —
        # without this, per-record _cigar_at(i) goes back through
        # the structured-record decode each call, dominating the
        # per-record time on the genomic transport encode path
        # (mirrors Java fix / ObjC parity).
        cigars = self._compound("cigars")
        out: list[str] = []
        for row in cigars:
            v = row["value"]
            if isinstance(v, bytes):
                out.append(v.decode("utf-8"))
            else:
                out.append(v)
        self._decoded_cigars = out
        return out[i]

    # ------------------------------------------------------------------
    # M86 Phase F — mate_info per-field dispatch
    # ------------------------------------------------------------------

    def _mate_info_is_subgroup(self) -> bool:
        """True iff ``signal_channels/mate_info`` is a group (Phase F).

        Per / Gotcha §141, dispatch is on HDF5
        link type, NOT on ``@compression`` attribute presence on the
        bare link. The StorageGroup protocol's ``open_group`` raises
        ``KeyError`` when the named child is a dataset (verified in
        :class:`ttio.providers.hdf5._Group.open_group`); we use that
        as the link-type query.

        Result is cached on the instance — the file structure is
        immutable for the lifetime of an open run, so the link-type
        probe only runs once.
        """
        if self._mate_info_subgroup_cached is not None:
            return self._mate_info_subgroup_cached
        sig = self._signal_channels_group()
        try:
            sig.open_group("mate_info")
            self._mate_info_subgroup_cached = True
        except KeyError:
            self._mate_info_subgroup_cached = False
        return self._mate_info_subgroup_cached

    def _mate_info_is_inline_v2(self) -> bool:
        """True iff signal_channels/mate_info/inline_v2 exists.

        v1.7+ inline-v2 layout. Implies _mate_info_is_subgroup() == True.
        Result is cached via _decoded_mate_info sentinel key.
        """
        if not self._mate_info_is_subgroup():
            return False
        cached = self._decoded_mate_info.get("_is_inline_v2")
        if cached is not None:
            return cached
        sig = self._signal_channels_group()
        mate_group = sig.open_group("mate_info")
        try:
            mate_group.open_dataset("inline_v2")
            result = True
        except KeyError:
            result = False
        self._decoded_mate_info["_is_inline_v2"] = result
        return result

    def _decode_mate_inline_v2(self) -> "dict[str, Any]":
        """Decode the inline_v2 blob; cache all three fields together.

        Returns the _decoded_mate_info dict populated with keys:
        'chrom' (list[str]), 'pos' (np.int64 array), 'tlen' (np.int32
        array). Reader-side dependency: own_positions and own_chrom_ids
        reconstructed from genomic_index; chrom name resolution uses
        the mate_info/chrom_names table written alongside the blob.
        """
        if "inline_v2" in self._decoded_mate_info:
            return self._decoded_mate_info

        from .codecs._registry import CODEC_REGISTRY
        from .codecs._context import ChannelPayload

        sig = self._signal_channels_group()
        mate_group = sig.open_group("mate_info")
        ds = mate_group.open_dataset("inline_v2")
        blob = bytes(np.asarray(ds.read(offset=0, count=int(ds.length))).tobytes())

        decoded = CODEC_REGISTRY[Compression.MATE_INLINE_V2].decode(
            ChannelPayload.of_bytes(blob), self._codec_context()
        ).as_mate_info()
        mc = decoded["mate_chrom_ids"]
        mp = decoded["mate_positions"]
        ts = decoded["template_lengths"]

        # Read the full chrom_id → name table written by the writer.
        # This covers mate-only chromosomes absent from genomic_index.
        chrom_name_rows = io.read_compound_dataset(mate_group, "chrom_names")
        chrom_names_by_id: list[str] = []
        for row in chrom_name_rows:
            v = row["name"]
            chrom_names_by_id.append(
                v.decode("utf-8") if isinstance(v, bytes) else v
            )

        # Convert mc (int32 ids, -1 = unmapped) back to mate_chromosomes list[str].
        mate_chromosomes: list[str] = []
        for v in mc:
            iv = int(v)
            if iv == -1:
                mate_chromosomes.append("*")
            elif 0 <= iv < len(chrom_names_by_id):
                mate_chromosomes.append(chrom_names_by_id[iv])
            else:
                # Should not happen in well-formed files.
                mate_chromosomes.append(f"chr_id_{iv}")

        self._decoded_mate_info["chrom"] = mate_chromosomes
        self._decoded_mate_info["pos"] = mp
        self._decoded_mate_info["tlen"] = ts
        self._decoded_mate_info["inline_v2"] = True  # marker for round-trip cache
        return self._decoded_mate_info

    def _raise_unsupported_mate_layout(self) -> None:
        """Raise a clear migration error for any non-inline_v2 mate_info layout.

        v1.0 reset (Phase 2c): the M86 Phase F per-field subgroup
        layout (with chrom / pos / tlen child datasets) and the
        legacy M82 compound dataset layout were both removed. Only
        the v1.7+ inline_v2 BLOB path under
        ``signal_channels/mate_info/inline_v2`` (codec id 13) is
        decoded; everything else surfaces a clear error so callers
        learn they need to re-encode with v1.0+.
        """
        raise ValueError(
            "signal_channels/mate_info: legacy layout detected — "
            "the v1.0 reader requires the inline_v2 blob "
            "(signal_channels/mate_info/inline_v2 with @compression=13). "
            "The M86 Phase F per-field subgroup (chrom/pos/tlen) and "
            "the M82 compound dataset layouts were removed in v1.0. "
            "This file was written with an older TTI-O version; "
            "re-encode with v1.0+ to use MATE_INLINE_V2."
        )

    def _mate_chrom_at(self, i: int) -> str:
        """Return the mate chromosome at index ``i``.

        v1.0 reset (Phase 2c): only the v1.7+ inline_v2 layout is
        supported. Any other layout raises ValueError.
        """
        if self._mate_info_is_inline_v2():
            return self._decode_mate_inline_v2()["chrom"][i]
        self._raise_unsupported_mate_layout()
        raise AssertionError("unreachable")  # pragma: no cover

    def _mate_pos_at(self, i: int) -> int:
        """Return the mate position at index ``i``.

        v1.0 reset (Phase 2c): inline_v2 is the only supported layout.
        """
        if self._mate_info_is_inline_v2():
            return int(self._decode_mate_inline_v2()["pos"][i])
        self._raise_unsupported_mate_layout()
        raise AssertionError("unreachable")  # pragma: no cover

    def _mate_tlen_at(self, i: int) -> int:
        """Return the template length at index ``i``.

        v1.0 reset (Phase 2c): inline_v2 is the only supported layout.
        """
        if self._mate_info_is_inline_v2():
            return int(self._decode_mate_inline_v2()["tlen"][i])
        self._raise_unsupported_mate_layout()
        raise AssertionError("unreachable")  # pragma: no cover
