"""Shared transport-codec helpers used by both writer and reader.

Pure code-movement split of ``codec.py`` (OO-assessment P3.10).
Holds the wire-mapping constants and the genomic/codec helpers that
:mod:`ttio.transport._writer` and :mod:`ttio.transport._reader` both
depend on. Importing this module must not pull in either of those
(one-way dependency: writer/reader -> _common).
"""
from __future__ import annotations

import logging
import struct
from typing import Iterator

from ..enums import Compression, Polarity, Precision
from .._hdf5_io import read_int_attr as io_attr_int  # M90.10 wire codec probe
from .packets import (
    AccessUnit,
    ChannelData,
)

_CHECKSUM_STRUCT = struct.Struct("<I")

_LOG = logging.getLogger("ttio.transport.codec")

# ---------------------------------------------------------- wire mappings

# Wire polarity (0=positive, 1=negative, 2=unknown) vs Polarity enum
# (UNKNOWN=0, POSITIVE=1, NEGATIVE=-1). The transport spec uses a
# nonneg-only layout for portability across languages that can't
# round-trip negative uint8 values.
_POLARITY_TO_WIRE = {
    Polarity.POSITIVE: 0,
    Polarity.NEGATIVE: 1,
    Polarity.UNKNOWN: 2,
}
_WIRE_TO_POLARITY = {v: k for k, v in _POLARITY_TO_WIRE.items()}

_SPECTRUM_CLASS_TO_WIRE = {
    "TTIOMassSpectrum": 0,
    "TTIONMRSpectrum": 1,
    "TTIONMR2DSpectrum": 2,
    "TTIOFreeInductionDecay": 3,
    "TTIOMSImagePixel": 4,
    "TTIOGenomicRead": 5,  # M89.2
}
_WIRE_TO_SPECTRUM_CLASS = {v: k for k, v in _SPECTRUM_CLASS_TO_WIRE.items()}


# M86 codec dispatch for genomic UINT8 channels on the wire.

_RANS_ORDER0_WIRE = int(Compression.RANS_ORDER0)
_RANS_ORDER1_WIRE = int(Compression.RANS_ORDER1)
_BASE_PACK_WIRE = int(Compression.BASE_PACK)


def _read_mate_chrom_names_table(mate_grp) -> list[str]:
    """Read the ``mate_info/chrom_names`` compound dataset into a list
    of names ordered by row index (matches the chrom_id encoding used
    by MATE_INLINE_V2). Returns ``[]`` when the table is missing
    (mate_info/inline_v2 is allowed without a name table on empty
    runs)."""
    from .. import _hdf5_io as _io
    try:
        records = _io.read_compound_dataset(mate_grp, "chrom_names")
    except (KeyError, ValueError):
        return []
    return [str(r.get("name", "")) for r in records]


_WIRE_CODEC_IDS = frozenset({0, 4, 5, 6})   # NONE, RANS_ORDER0, RANS_ORDER1, BASE_PACK


def _apply_wire_codec(plaintext: bytes, codec: int) -> bytes:
    """Encode ``plaintext`` with the wire codec id (NONE → identity)."""
    if codec == 0:  # NONE
        return plaintext
    if codec == _RANS_ORDER0_WIRE:
        from ..codecs import rans
        return rans.encode(plaintext, order=0)
    if codec == _RANS_ORDER1_WIRE:
        from ..codecs import rans
        return rans.encode(plaintext, order=1)
    if codec == _BASE_PACK_WIRE:
        from ..codecs import base_pack
        return base_pack.encode(plaintext)
    # Other compression ids (zlib for MS, etc.) take the existing
    # paths in this module; this helper is genomic-channel-only.
    raise NotImplementedError(
        f"_apply_wire_codec: codec id {codec} not supported for genomic UINT8"
    )


def _decode_wire_codec(payload: bytes, codec: int) -> bytes:
    """Decode a payload encoded by :func:`_apply_wire_codec`."""
    if codec == 0:
        return payload
    if codec == _RANS_ORDER0_WIRE or codec == _RANS_ORDER1_WIRE:
        from ..codecs import rans
        return rans.decode(payload)
    if codec == _BASE_PACK_WIRE:
        from ..codecs import base_pack
        return base_pack.decode(payload)
    raise NotImplementedError(
        f"_decode_wire_codec: codec id {codec} not supported for genomic UINT8"
    )


def _run_has_tags(run) -> bool:
    """True when a stored genomic run carries the SAM tags channel (M101)."""
    try:
        return run.group.open_group("signal_channels").has_child("tags")
    except Exception:
        return False


def genomic_channel_names(run) -> list[str]:
    """The DatasetHeader channel list of a genomic run: the five M89
    channels, plus ``tags`` (transport-spec 4.3.1) when the run carries
    SAM tags, so a tag-less run's stream is unchanged."""
    names = ["sequences", "qualities", "cigar", "read_name", "mate_chromosome"]
    if _run_has_tags(run):
        names.append("tags")
    return names


def _iter_genomic_run_access_units(run) -> Iterator[tuple[int, "AccessUnit"]]:
    """Yield ``(au_sequence, AccessUnit)`` tuples for every AlignedRead
    in ``run``.

    Body extracted from :meth:`TransportWriter._emit_genomic_run_access_units`
    in #141 so :func:`ttio.transport.walker.walk_dataset` can emit
    genomic AUs without duplicating the per-read construction. Both
    callers see byte-identical AUs.
    """
    n_reads = len(run)
    precision_uint8 = int(Precision.UINT8) & 0xFF
    compression_none = int(Compression.NONE) & 0xFF
    acq_mode = int(run.acquisition_mode) & 0xFF
    # M90.10 wire codec: a whole-channel run whose sequences/qualities
    # dataset names an M86 codec is re-encoded per AU with that codec.
    # A blocks_v1 run (format-spec 10.12) keeps its blobs per block and
    # is sent with the plain per-AU channel bytes.
    seq_codec = qual_codec = compression_none
    blocks = getattr(run, "layout", "whole") in ("blocks_v1", "blocks_v1_grouped")
    try:
        sig_group = run.group.open_group("signal_channels")
        if sig_group.has_child("sequences"):
            if blocks:
                seq_ds = sig_group.open_group("sequences").open_dataset("data")
            else:
                seq_ds = sig_group.open_dataset("sequences")
            seq_codec = (io_attr_int(seq_ds, "compression",
                                        default=0) or 0) & 0xFF
        if sig_group.has_child("qualities"):
            qual_ds = sig_group.open_dataset("qualities")
            qual_codec = (io_attr_int(qual_ds, "compression",
                                         default=0) or 0) & 0xFF
    except Exception:
        seq_codec = qual_codec = compression_none
    # Only the per-AU wire codecs propagate (RANS_ORDER0/1, BASE_PACK);
    # a whole-run codec (REF_DIFF_V2, FQZCOMP_NX16_Z, ...) has no per-AU
    # form and those channels are sent as plain bytes.
    if seq_codec not in _WIRE_CODEC_IDS:
        seq_codec = compression_none
    if qual_codec not in _WIRE_CODEC_IDS:
        qual_codec = compression_none
    has_tags = _run_has_tags(run)

    # iter_reads holds one decoded block at a time for blocks_v1 and
    # walks the whole-channel caches for the v1.8 layout.
    for i, r in enumerate(run.iter_reads()):
        seq_bytes = (r.sequence or "").encode("ascii")
        qual_bytes = bytes(r.qualities)
        length = len(seq_bytes)
        seq_payload = _apply_wire_codec(seq_bytes, seq_codec)
        qual_payload = _apply_wire_codec(qual_bytes, qual_codec)
        cigar_bytes = (r.cigar or "").encode("utf-8")
        name_bytes = (r.read_name or "").encode("utf-8")
        mate_chr_bytes = (r.mate_chromosome or "").encode("utf-8")
        channels = [
            ChannelData("sequences", precision_uint8,
                        seq_codec, length, seq_payload),
            ChannelData("qualities", precision_uint8,
                        qual_codec, length, qual_payload),
            ChannelData("cigar", precision_uint8,
                        compression_none, len(cigar_bytes), cigar_bytes),
            ChannelData("read_name", precision_uint8,
                        compression_none, len(name_bytes), name_bytes),
            ChannelData("mate_chromosome", precision_uint8,
                        compression_none, len(mate_chr_bytes),
                        mate_chr_bytes),
        ]
        if has_tags:
            tag_bytes = (r.tags or "").encode("utf-8")
            channels.append(ChannelData("tags", precision_uint8, compression_none,
                                        len(tag_bytes), tag_bytes))
        au = AccessUnit(
            spectrum_class=5,
            acquisition_mode=acq_mode,
            ms_level=0,
            polarity=2,
            retention_time=0.0,
            precursor_mz=0.0,
            precursor_charge=0,
            ion_mobility=0.0,
            base_peak_intensity=0.0,
            channels=channels,
            chromosome=r.chromosome,
            position=int(r.position),
            mapping_quality=int(r.mapping_quality),
            flags=int(r.flags) & 0xFFFF,
            mate_position=int(r.mate_position),
            template_length=int(r.template_length),
        )
        yield i, au
