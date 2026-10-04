"""Block encoder for the ``blocks_v1`` genomic layout.

A block is a run consisting of a contiguous range of reads. Its
channel blobs are produced by running the whole-run writer
(:func:`ttio._dataset_write_genomic._write_genomic_run`) against an
in-memory provider group and harvesting each channel's dataset bytes
and ``@compression``; the codecs and their wire formats are untouched
and a block's blob is byte-identical to what a v1.8 whole-run write of
those reads would produce (format-spec 10.12).

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np

from ..written_genomic_run import WrittenGenomicRun

#: Blob channels every blocks_v1 run indexes, in block-index column order.
REQUIRED_BLOCK_CHANNELS = ("sequences", "qualities", "read_names", "cigars", "mate_info")
#: Channels whose index columns exist only in runs that carry them (M101:
#: the SAM tags). Their (off, len, codec) columns follow the required
#: channels' codec columns; a reader treats missing columns as an empty
#: channel (format-spec 10.12.2).
OPTIONAL_BLOCK_CHANNELS = ("tags",)
#: Every blob channel of a block.
BLOCK_CHANNELS = REQUIRED_BLOCK_CHANNELS + OPTIONAL_BLOCK_CHANNELS


@dataclass
class BlockBlobs:
    blobs: dict[str, bytes]
    compression: dict[str, int]
    seq_layout: str                       # "refdiff_v2" | "raw"
    extra_attrs: dict[str, dict] = field(default_factory=dict)
    n_reads: int = 0
    n_bases: int = 0


def slice_run(run: WrittenGenomicRun, start: int, stop: int) -> WrittenGenomicRun:
    """Reads ``[start, stop)`` of ``run`` as a run of their own, offsets
    rebased to 0. Run-level metadata is shared by reference."""
    if stop <= start:
        b0 = b1 = 0
    else:
        b0 = int(run.offsets[start])
        b1 = int(run.offsets[stop - 1]) + int(run.lengths[stop - 1])
    return dataclasses.replace(
        run,
        positions=run.positions[start:stop],
        mapping_qualities=run.mapping_qualities[start:stop],
        flags=run.flags[start:stop],
        sequences=run.sequences[b0:b1],
        qualities=run.qualities[b0:b1],
        offsets=(np.asarray(run.offsets[start:stop], dtype=np.uint64) - np.uint64(b0)),
        lengths=run.lengths[start:stop],
        cigars=list(run.cigars[start:stop]),
        read_names=list(run.read_names[start:stop]),
        mate_chromosomes=list(run.mate_chromosomes[start:stop]),
        mate_positions=run.mate_positions[start:stop],
        template_lengths=run.template_lengths[start:stop],
        chromosomes=list(run.chromosomes[start:stop]),
        tags=None if run.tags is None else list(run.tags[start:stop]),
        provenance_records=[],
    )


def concat_runs(parts: list[WrittenGenomicRun]) -> WrittenGenomicRun:
    """The inverse of :func:`slice_run` for a list of consecutive parts."""
    if len(parts) == 1:
        return parts[0]
    first = parts[0]
    seqs = np.concatenate([p.sequences for p in parts])
    quals = np.concatenate([p.qualities for p in parts])
    lengths = np.concatenate([p.lengths for p in parts]).astype(np.uint32)
    offsets = np.zeros(len(lengths), dtype=np.uint64)
    if len(lengths) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.uint64)
    return dataclasses.replace(
        first,
        positions=np.concatenate([p.positions for p in parts]),
        mapping_qualities=np.concatenate([p.mapping_qualities for p in parts]),
        flags=np.concatenate([p.flags for p in parts]),
        sequences=seqs, qualities=quals, offsets=offsets, lengths=lengths,
        cigars=[c for p in parts for c in p.cigars],
        read_names=[n for p in parts for n in p.read_names],
        mate_chromosomes=[m for p in parts for m in p.mate_chromosomes],
        mate_positions=np.concatenate([p.mate_positions for p in parts]),
        template_lengths=np.concatenate([p.template_lengths for p in parts]),
        chromosomes=[c for p in parts for c in p.chromosomes],
        tags=(None if all(p.tags is None for p in parts)
              else [t for p in parts for t in (p.tags or [""] * len(p.cigars))]),
        provenance_records=[],
    )


def take_run(run: WrittenGenomicRun, order: np.ndarray) -> WrittenGenomicRun:
    """The reads of ``run`` in the order ``order`` (read ``order[j]`` at row
    ``j``), offsets recomputed. Run-level metadata is shared by reference."""
    order = np.asarray(order, dtype=np.int64)
    lengths = np.asarray(run.lengths, dtype=np.uint32)[order]
    starts = np.asarray(run.offsets, dtype=np.int64)[order]
    offsets = np.zeros(len(order), dtype=np.uint64)
    if len(order) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.uint64)
    total = int(lengths.sum(dtype=np.uint64))
    if total:
        # Byte gather: for each output byte, its source byte.
        rep = np.repeat(starts - offsets.astype(np.int64), lengths.astype(np.int64))
        src = np.arange(total, dtype=np.int64) + rep
        seqs = np.asarray(run.sequences)[src]
        quals = np.asarray(run.qualities)[src] if len(run.qualities) else run.qualities[:0]
    else:
        seqs = np.asarray(run.sequences)[:0]
        quals = np.asarray(run.qualities)[:0]
    idx = order.tolist()
    return dataclasses.replace(
        run,
        positions=np.asarray(run.positions)[order],
        mapping_qualities=np.asarray(run.mapping_qualities)[order],
        flags=np.asarray(run.flags)[order],
        sequences=seqs, qualities=quals, offsets=offsets, lengths=lengths,
        cigars=[run.cigars[i] for i in idx],
        read_names=[run.read_names[i] for i in idx],
        mate_chromosomes=[run.mate_chromosomes[i] for i in idx],
        mate_positions=np.asarray(run.mate_positions)[order],
        template_lengths=np.asarray(run.template_lengths)[order],
        chromosomes=[run.chromosomes[i] for i in idx],
        tags=None if run.tags is None else [run.tags[i] for i in idx],
        provenance_records=[],
    )


def group_order(run: WrittenGenomicRun) -> np.ndarray:
    """The read-grouping permutation of ``run`` (M103): reads are grouped
    by sequence within each chromosome label, labels in first-seen order,
    so a block still never spans two labels. Entry ``j`` is the input
    index of the read stored at row ``j`` (uint32)."""
    from ..codecs import seq_group
    n = len(run.lengths)
    if n == 0:
        return np.zeros(0, dtype=np.uint32)
    labels: dict[str, list[int]] = {}
    for i, c in enumerate(run.chromosomes):
        labels.setdefault(c, []).append(i)
    parts = []
    for members in labels.values():
        sub = run if len(labels) == 1 else take_run(run, np.asarray(members))
        local = seq_group.group(bytes(np.asarray(sub.sequences, dtype=np.uint8).tobytes()),
                                np.asarray(sub.lengths, dtype=np.uint64), list(sub.read_names))
        parts.append(np.asarray(members, dtype=np.uint32)[local])
    return np.concatenate(parts).astype(np.uint32)


def _try_group(parent, name: str):
    try:
        return parent.open_group(name)
    except Exception:
        return None


def _harvest(ds) -> tuple[bytes, int, dict]:
    attrs = {k: ds.get_attribute(k) for k in ds.attribute_names() if k != "compression"}
    comp = int(ds.get_attribute("compression")) if ds.has_attribute("compression") else 0
    return bytes(np.asarray(ds.read(), dtype=np.uint8).tobytes()), comp, attrs


def encode_block(block: WrittenGenomicRun,
                 qual_strategy_hint: int = -1) -> BlockBlobs:
    """Encode one block's channels through the whole-run writer.

    ``qual_strategy_hint``: -1 auto (default), 5/6 forced V5, 7 V4 with
    internal preset selection; the stream writer passes its per-run pin.
    """
    from .._dataset_write_genomic import _write_genomic_run
    from ..providers.memory import MemoryProvider

    from ..enums import Compression

    # Under blocks_v1 every blob channel is a flat byte stream. The
    # v1.8 default for cigars is the M82 compound VL-string dataset,
    # which has no blob form, so blocks_v1 defaults cigars to the
    # existing RANS_ORDER0 wiring (format-spec 10.8) unless the caller
    # chose a cigars codec.
    overrides = dict(block.signal_codec_overrides)
    if "cigars" not in overrides:
        overrides["cigars"] = Compression.RANS_ORDER0
    # blocks_v1 always codes qualities with FQZCOMP_NX16_Z (the v1.8
    # writer only does so when another v1.5 codec is active on the run,
    # else it leaves zlib-filtered raw bytes, which have no blob form).
    if "qualities" not in overrides:
        # FQZCOMP_NX16_Z: the v1.8 kernel fails to decode a stream in
        # which a zero-length read (SEQ '*', e.g. a secondary
        # alignment) precedes another read (rc=-6); such blocks use
        # RANS_ORDER0 for qualities. Recorded per block in the
        # qualities_codec index column.
        lens = np.asarray(block.lengths)
        if len(lens) and (lens == 0).any():
            overrides["qualities"] = Compression.RANS_ORDER0
        else:
            overrides["qualities"] = Compression.FQZCOMP_NX16_Z
    # Sequences without a reference: SEQ_CM (M103), or RANS_ORDER1 when
    # libttio_rans is not loaded. With a reference the REF_DIFF_V2
    # default applies; unmapped reads inside a mapped block are carried
    # by the codec (UL substream).
    if "sequences" not in overrides and block.reference_chrom_seqs is None:
        from ..codecs import seq_cm
        overrides["sequences"] = (Compression.SEQ_CM if seq_cm.HAVE_NATIVE_LIB
                                  else Compression.RANS_ORDER1)
    if overrides != block.signal_codec_overrides:
        block = dataclasses.replace(block, signal_codec_overrides=overrides)
    root = MemoryProvider.open("memory://ttio-block-encode", mode="w").root_group()
    _write_genomic_run(root, "b", block, qual_strategy_hint)
    sc = root.open_group("b").open_group("signal_channels")
    out = BlockBlobs(blobs={}, compression={}, seq_layout="raw",
                     n_reads=int(len(block.lengths)),
                     n_bases=int(np.asarray(block.lengths, dtype=np.uint64).sum()))
    for ch in BLOCK_CHANNELS:
        ds = None
        if ch == "sequences":
            g = _try_group(sc, "sequences")
            if g is not None and g.has_child("refdiff_v2"):
                ds = g.open_dataset("refdiff_v2")
                out.seq_layout = "refdiff_v2"
            elif g is None and sc.has_child("sequences"):
                ds = sc.open_dataset("sequences")
        elif ch == "mate_info":
            g = _try_group(sc, "mate_info")
            if g is not None and g.has_child("inline_v2"):
                ds = g.open_dataset("inline_v2")
        elif sc.has_child(ch) and _try_group(sc, ch) is None:
            ds = sc.open_dataset(ch)
        if ds is None:
            out.blobs[ch] = b""
            out.compression[ch] = 0
            out.extra_attrs[ch] = {}
            continue
        out.blobs[ch], out.compression[ch], out.extra_attrs[ch] = _harvest(ds)
    return out
