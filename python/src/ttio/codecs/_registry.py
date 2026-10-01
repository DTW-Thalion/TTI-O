"""Codec registry: maps Compression ids to Codec adapters (codec-registry refactor).

Context-aware codecs are added in Task 3. No wire change — adapters wrap the
existing codec functions verbatim.
"""
from __future__ import annotations

from typing import Protocol

from ..enums import Compression
from . import base_pack, delta_rans, quality, rans
from . import fqzcomp_nx16_z, mate_info_v2, name_tokenizer_v2, ref_diff_v2, sam_tags
from ._context import ChannelPayload, CodecContext, DecodedChannel, EncodedChannel


class Codec(Protocol):
    id: Compression
    is_context_aware: bool
    needs_embedded_reference: bool

    def decode(self, payload: ChannelPayload, ctx: CodecContext) -> DecodedChannel: ...
    def encode(self, value: DecodedChannel, ctx: CodecContext) -> EncodedChannel: ...


class _RansCodec:
    """rANS O0/O1. decode() is order-agnostic (order is in the stream)."""
    def __init__(self, cid: Compression, order: int) -> None:
        self.id = cid
        self._order = order
        self.is_context_aware = False
        self.needs_embedded_reference = False

    def decode(self, payload, ctx):
        return DecodedChannel.of_bytes(rans.decode(payload.as_bytes()))

    def encode(self, value, ctx):
        return EncodedChannel.of_dataset(rans.encode(value.as_bytes(), order=self._order))


class _BasePackCodec:
    id = Compression.BASE_PACK
    is_context_aware = False
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        return DecodedChannel.of_bytes(base_pack.decode(payload.as_bytes()))

    def encode(self, value, ctx):
        return EncodedChannel.of_dataset(base_pack.encode(value.as_bytes()))


class _QualityBinnedCodec:
    id = Compression.QUALITY_BINNED
    is_context_aware = False
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        return DecodedChannel.of_bytes(quality.decode(payload.as_bytes()))

    def encode(self, value, ctx):
        return EncodedChannel.of_dataset(quality.encode(value.as_bytes()))


class _DeltaRansCodec:
    id = Compression.DELTA_RANS_ORDER0
    is_context_aware = False
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        return DecodedChannel.of_bytes(delta_rans.decode(payload.as_bytes()))

    def encode(self, value, ctx):
        if ctx.element_size is None:
            raise ValueError("DELTA_RANS encode requires CodecContext.element_size")
        return EncodedChannel.of_dataset(
            delta_rans.encode(value.as_bytes(), ctx.element_size)
        )


class _NameTokenizedV2Codec:
    id = Compression.NAME_TOKENIZED_V2
    is_context_aware = False  # str-list domain, but no run context needed
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        return DecodedChannel.of_str_list(name_tokenizer_v2.decode(payload.as_bytes()))

    def encode(self, value, ctx):
        return EncodedChannel.of_dataset(name_tokenizer_v2.encode(value.as_str_list()))


class _FqzcompNx16ZCodec:
    id = Compression.FQZCOMP_NX16_Z
    is_context_aware = True
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        flags = None
        if ctx.revcomp_flags is not None:
            flags = [int(x) for x in ctx.revcomp_flags]
        qualities, _read_lengths, _rc = fqzcomp_nx16_z.decode_with_metadata(
            payload.as_bytes(), flags,
            sequences_provider=ctx.sequences_provider,
        )
        return DecodedChannel.of_bytes(qualities)

    def encode(self, value, ctx):
        if ctx.read_lengths is None or ctx.revcomp_flags is None:
            raise ValueError(
                "FQZCOMP_NX16_Z encode requires CodecContext.read_lengths + revcomp_flags"
            )
        blob = fqzcomp_nx16_z.encode(
            value.as_bytes(),
            [int(x) for x in ctx.read_lengths],
            [int(x) for x in ctx.revcomp_flags],
            sequences=ctx.sequences,
            v4_strategy_hint=ctx.qual_strategy_hint,
        )
        return EncodedChannel.of_dataset(blob)


class _MateInlineV2Codec:
    id = Compression.MATE_INLINE_V2
    is_context_aware = True
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        if ctx.own_chrom_ids is None or ctx.own_positions is None or ctx.n_records is None:
            raise ValueError(
                "MATE_INLINE_V2 decode requires CodecContext.own_chrom_ids/own_positions/n_records"
            )
        mc, mp, tl = mate_info_v2.decode(
            payload.as_bytes(), ctx.own_chrom_ids, ctx.own_positions, ctx.n_records
        )
        return DecodedChannel.of_mate_info(
            {"mate_chrom_ids": mc, "mate_positions": mp, "template_lengths": tl}
        )

    def encode(self, value, ctx):
        if ctx.own_chrom_ids is None or ctx.own_positions is None:
            raise ValueError(
                "MATE_INLINE_V2 encode requires CodecContext.own_chrom_ids/own_positions"
            )
        d = value.as_mate_info()
        blob = mate_info_v2.encode(
            d["mate_chrom_ids"], d["mate_positions"], d["template_lengths"],
            ctx.own_chrom_ids, ctx.own_positions,
        )
        return EncodedChannel.of_dataset(blob)


class _RefDiffV2Codec:
    id = Compression.REF_DIFF_V2
    is_context_aware = True
    needs_embedded_reference = True

    def decode(self, payload, ctx):
        # Relocated from GenomicRun._decode_ref_diff_v2_sequences: parse the blob
        # header, resolve the reference via ctx, then decode.
        ds = payload.group().open_dataset("refdiff_v2")
        blob = bytes(ds.read(offset=0, count=int(ds.length)))
        header = ref_diff_v2.parse_blob_header(blob)
        if ctx.reference_resolver is None or ctx.chromosomes is None:
            raise ValueError("REF_DIFF_V2 decode requires CodecContext.reference_resolver + chromosomes")
        unique = set(ctx.chromosomes)
        if len(unique) == 0:
            chrom = ""
        elif len(unique) > 1:
            raise RuntimeError(
                "REF_DIFF_V2 supports single-chromosome runs only; "
                f"this run carries {sorted(unique)}.")
        else:
            chrom = next(iter(unique))
        chrom_seq = ctx.reference_resolver.resolve(
            uri=header.reference_uri, expected_md5=header.reference_md5, chromosome=chrom)
        cigars = ctx.cigars_provider() if ctx.cigars_provider else []
        out_seq, _off = ref_diff_v2.decode(
            blob, ctx.positions, cigars, chrom_seq, ctx.read_count, ctx.total_bases)
        return DecodedChannel.of_bytes(bytes(out_seq))

    def encode(self, value, ctx):
        # Task 5c: encode a single-chromosome refdiff_v2 blob. The
        # reference inputs (offsets/reference/md5/uri/cigars/positions)
        # live on the encode-only CodecContext fields because they are
        # written *into* the blob header. Returns a GROUP layout: the
        # writer materialises a ``sequences`` group with one
        # ``refdiff_v2`` child carrying these bytes + @compression. The
        # bytes are byte-identical to the prior direct ref_diff_v2.encode
        # call (same args, same order).
        if (
            ctx.offsets is None
            or ctx.positions is None
            or ctx.cigar_strings is None
            or ctx.reference is None
            or ctx.reference_md5 is None
            or ctx.reference_uri is None
        ):
            raise ValueError(
                "REF_DIFF_V2 encode requires CodecContext.offsets/positions/"
                "cigar_strings/reference/reference_md5/reference_uri")
        import numpy as _np
        blob = ref_diff_v2.encode(
            _np.frombuffer(value.as_bytes(), dtype=_np.uint8),
            _np.asarray(ctx.offsets, dtype=_np.uint64),
            _np.asarray(ctx.positions, dtype=_np.int64),
            list(ctx.cigar_strings),
            ctx.reference,
            ctx.reference_md5,
            ctx.reference_uri,
            reads_per_slice=(
                ctx.reads_per_slice if ctx.reads_per_slice is not None else 10_000
            ),
            slice_bytes=(ctx.slice_bytes if ctx.slice_bytes is not None else 0),
        )
        return EncodedChannel.of_group({"refdiff_v2": blob}, {})


class _SamTagsCodec:
    """SAM_TAGS (M101): per-read tag text. MD/NM derivation needs the
    reads' sequences, CIGARs, positions and chromosome ids plus the
    reference bases; without ``tag_references`` none of that is read."""
    id = Compression.SAM_TAGS
    is_context_aware = True
    # Derivation reuses the reference REF_DIFF_V2 already requires; the
    # codec never makes a writer embed one (binding decision §98).
    needs_embedded_reference = False

    def decode(self, payload, ctx):
        import numpy as _np
        n = int(ctx.read_count or 0)
        refs = ctx.tag_references_provider() if ctx.tag_references_provider else None
        kw = {}
        if refs and any(r is not None for r in refs):
            lens = _np.asarray(ctx.read_lengths, dtype=_np.uint64)
            off = _np.zeros(n + 1, dtype=_np.uint64)
            _np.cumsum(lens, out=off[1:])
            kw = dict(sequences=ctx.sequences_provider(), seq_offsets=off,
                      cigars=ctx.cigars_provider(), positions=ctx.positions,
                      chrom_ids=ctx.own_chrom_ids, references=refs)
        return DecodedChannel.of_str_list(sam_tags.decode(payload.as_bytes(), n, **kw))

    def encode(self, value, ctx):
        kw = {}
        if ctx.tag_references and any(r is not None for r in ctx.tag_references):
            kw = dict(sequences=ctx.sequences, seq_offsets=ctx.offsets,
                      cigars=ctx.cigar_strings, positions=ctx.positions,
                      chrom_ids=ctx.own_chrom_ids, references=ctx.tag_references)
        return EncodedChannel.of_dataset(sam_tags.encode(value.as_str_list(), **kw))


CODEC_REGISTRY: "dict[Compression, Codec]" = {
    Compression.RANS_ORDER0: _RansCodec(Compression.RANS_ORDER0, 0),
    Compression.RANS_ORDER1: _RansCodec(Compression.RANS_ORDER1, 1),
    Compression.BASE_PACK: _BasePackCodec(),
    Compression.QUALITY_BINNED: _QualityBinnedCodec(),
    Compression.DELTA_RANS_ORDER0: _DeltaRansCodec(),
    Compression.NAME_TOKENIZED_V2: _NameTokenizedV2Codec(),
    Compression.FQZCOMP_NX16_Z: _FqzcompNx16ZCodec(),
    Compression.MATE_INLINE_V2: _MateInlineV2Codec(),
    Compression.REF_DIFF_V2: _RefDiffV2Codec(),
    Compression.SAM_TAGS: _SamTagsCodec(),
}
