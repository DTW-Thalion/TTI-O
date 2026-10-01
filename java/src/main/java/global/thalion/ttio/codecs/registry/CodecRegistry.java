package global.thalion.ttio.codecs.registry;

import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.codecs.BasePack;
import global.thalion.ttio.codecs.DeltaRans;
import global.thalion.ttio.codecs.FqzcompNx16Z;
import global.thalion.ttio.codecs.MateInfoV2;
import global.thalion.ttio.codecs.NameTokenizerV2;
import global.thalion.ttio.codecs.Quality;
import global.thalion.ttio.codecs.Rans;
import global.thalion.ttio.codecs.RefDiffV2;
import global.thalion.ttio.codecs.SamTags;
import java.util.EnumMap;
import java.util.Map;

/** Maps Compression ids to Codec adapters. Adapters wrap the existing static
 *  codec classes verbatim — no wire change. */
public final class CodecRegistry {
    private CodecRegistry() {}

    public static final Map<Compression, Codec> CODEC_REGISTRY = build();

    private static byte[] bytes(DecodedChannel v) {
        return ((DecodedChannel.Bytes) v).data();
    }
    private static byte[] payloadBytes(ChannelPayload p) {
        return ((ChannelPayload.BytesPayload) p).bytes();
    }

    static final class RansCodec implements Codec {
        private final Compression id;
        private final int order;
        RansCodec(Compression id, int order) { this.id = id; this.order = order; }
        public Compression id() { return id; }
        public boolean isContextAware() { return false; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            return new DecodedChannel.Bytes(Rans.decode(payloadBytes(p)));
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            return new EncodedChannel.DatasetBytes(Rans.encode(bytes(v), order));
        }
    }

    static final class BasePackCodec implements Codec {
        public Compression id() { return Compression.BASE_PACK; }
        public boolean isContextAware() { return false; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            return new DecodedChannel.Bytes(BasePack.decode(payloadBytes(p)));
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            return new EncodedChannel.DatasetBytes(BasePack.encode(bytes(v)));
        }
    }

    static final class QualityCodec implements Codec {
        public Compression id() { return Compression.QUALITY_BINNED; }
        public boolean isContextAware() { return false; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            return new DecodedChannel.Bytes(Quality.decode(payloadBytes(p)));
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            return new EncodedChannel.DatasetBytes(Quality.encode(bytes(v)));
        }
    }

    static final class DeltaRansCodec implements Codec {
        public Compression id() { return Compression.DELTA_RANS_ORDER0; }
        public boolean isContextAware() { return false; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            return new DecodedChannel.Bytes(DeltaRans.decode(payloadBytes(p)));
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            if (ctx.elementSize() == null) {
                throw new IllegalArgumentException(
                    "DELTA_RANS encode requires CodecContext.elementSize");
            }
            return new EncodedChannel.DatasetBytes(DeltaRans.encode(bytes(v), ctx.elementSize()));
        }
    }

    static final class NameTokenizedCodec implements Codec {
        public Compression id() { return Compression.NAME_TOKENIZED_V2; }
        public boolean isContextAware() { return false; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            return new DecodedChannel.StrList(NameTokenizerV2.decode(payloadBytes(p)));
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            return new EncodedChannel.DatasetBytes(
                NameTokenizerV2.encode(((DecodedChannel.StrList) v).names()));
        }
    }

    static final class FqzcompCodec implements Codec {
        public Compression id() { return Compression.FQZCOMP_NX16_Z; }
        public boolean isContextAware() { return true; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            FqzcompNx16Z.DecodeResult dr =
                FqzcompNx16Z.decode(payloadBytes(p), ctx.revcompFlags(),
                                    ctx.sequencesProvider());
            return new DecodedChannel.Bytes(dr.qualities());
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            if (ctx.readLengths() == null || ctx.revcompFlags() == null) {
                throw new IllegalArgumentException(
                    "FQZCOMP_NX16_Z encode requires CodecContext.readLengths + revcompFlags");
            }
            int hint = ctx.qualStrategyHint() != null
                ? ctx.qualStrategyHint() : -1;
            var opts = hint != -1
                ? new FqzcompNx16Z.EncodeOptions().v4StrategyHint(hint)
                : null;
            return new EncodedChannel.DatasetBytes(
                FqzcompNx16Z.encode(bytes(v), ctx.readLengths(),
                                    ctx.revcompFlags(), ctx.sequences(),
                                    opts));
        }
    }

    static final class MateInfoCodec implements Codec {
        public Compression id() { return Compression.MATE_INLINE_V2; }
        public boolean isContextAware() { return true; }
        public boolean needsEmbeddedReference() { return false; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            if (ctx.ownChromIds() == null || ctx.ownPositions() == null || ctx.nRecords() == null) {
                throw new IllegalArgumentException(
                    "MATE_INLINE_V2 decode requires ownChromIds/ownPositions/nRecords");
            }
            MateInfoV2.Triple t = MateInfoV2.decode(
                payloadBytes(p), ctx.ownChromIds(), ctx.ownPositions(), ctx.nRecords());
            return new DecodedChannel.MateInfo(
                t.mateChromIds, t.matePositions, t.templateLengths);
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            if (ctx.ownChromIds() == null || ctx.ownPositions() == null) {
                throw new IllegalArgumentException(
                    "MATE_INLINE_V2 encode requires ownChromIds/ownPositions");
            }
            DecodedChannel.MateInfo mi = (DecodedChannel.MateInfo) v;
            return new EncodedChannel.DatasetBytes(MateInfoV2.encode(
                mi.mateChromIds(), mi.matePositions(), mi.templateLengths(),
                ctx.ownChromIds(), ctx.ownPositions()));
        }
    }

    static final class RefDiffCodec implements Codec {
        public Compression id() { return Compression.REF_DIFF_V2; }
        public boolean isContextAware() { return true; }
        public boolean needsEmbeddedReference() { return true; }
        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            var group = ((ChannelPayload.GroupPayload) p).group();
            var ds = group.openDataset("refdiff_v2");
            byte[] blob = (byte[]) ds.readSlice(0L, ds.shape()[0]);
            RefDiffV2.BlobHeader header = RefDiffV2.parseBlobHeader(blob);
            if (ctx.referenceResolver() == null || ctx.chromosomes() == null) {
                throw new IllegalArgumentException(
                    "REF_DIFF_V2 decode requires referenceResolver + chromosomes");
            }
            java.util.LinkedHashSet<String> uniq =
                new java.util.LinkedHashSet<>(java.util.Arrays.asList(ctx.chromosomes()));
            String chrom;
            if (uniq.isEmpty()) chrom = "";
            else if (uniq.size() > 1) throw new IllegalStateException(
                "REF_DIFF_V2 supports single-chromosome runs only; this run carries " + uniq);
            else chrom = uniq.iterator().next();
            byte[] reference = ctx.referenceResolver().resolve(
                header.referenceUri(), header.referenceMd5(), chrom);
            String[] cigars = ctx.cigarsProvider() != null
                ? ctx.cigarsProvider().get() : new String[0];
            RefDiffV2.Pair out = RefDiffV2.decode(
                blob, ctx.positions(), cigars, reference, ctx.readCount(), ctx.totalBases());
            return new DecodedChannel.Bytes(out.sequences);
        }
        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            if (ctx.offsets() == null || ctx.positions() == null
                    || ctx.reference() == null || ctx.referenceMd5() == null
                    || ctx.cigarsProvider() == null) {
                throw new IllegalArgumentException(
                    "REF_DIFF_V2 encode requires offsets/positions/reference/"
                    + "referenceMd5/cigarsProvider");
            }
            int readsPerSlice = ctx.readsPerSlice() != null
                ? ctx.readsPerSlice() : 10_000;
            long sliceBytes = ctx.sliceBytes() != null ? ctx.sliceBytes() : 0L;
            byte[] blob = RefDiffV2.encode(
                bytes(v), ctx.offsets(), ctx.positions(),
                ctx.cigarsProvider().get(), ctx.reference(),
                ctx.referenceMd5(), ctx.referenceUri(), readsPerSlice,
                sliceBytes);
            return new EncodedChannel.GroupLayout(
                new java.util.LinkedHashMap<>(
                    java.util.Map.of("refdiff_v2", blob)),
                new java.util.LinkedHashMap<>());
        }
    }

    /** SAM_TAGS (M101): per-read tag text. MD/NM derivation needs the
     *  reads' sequences, CIGARs, positions and chromosome ids plus the
     *  reference bases; without tag references none of that is read. */
    static final class SamTagsCodec implements Codec {
        public Compression id() { return Compression.SAM_TAGS; }
        public boolean isContextAware() { return true; }
        public boolean needsEmbeddedReference() { return false; }

        private static boolean anyRef(java.util.List<byte[]> refs) {
            if (refs == null) return false;
            for (byte[] r : refs) if (r != null) return true;
            return false;
        }

        public DecodedChannel decode(ChannelPayload p, CodecContext ctx) {
            int n = ctx.readCount() != null ? ctx.readCount() : 0;
            java.util.List<byte[]> refs = ctx.tagReferencesProvider() != null
                ? ctx.tagReferencesProvider().get() : null;
            SamTags.Context tc = SamTags.Context.none();
            if (anyRef(refs)) {
                long[] off = new long[n + 1];
                for (int i = 0; i < n; i++) off[i + 1] = off[i] + ctx.readLengths()[i];
                tc = new SamTags.Context(ctx.sequencesProvider().get(), off,
                    java.util.Arrays.asList(ctx.cigarsProvider().get()),
                    ctx.positions(), ctx.ownChromIds(), refs);
            }
            return new DecodedChannel.StrList(SamTags.decode(payloadBytes(p), n, tc));
        }

        public EncodedChannel encode(DecodedChannel v, CodecContext ctx) {
            SamTags.Context tc = SamTags.Context.none();
            if (anyRef(ctx.tagReferences())) {
                tc = new SamTags.Context(ctx.sequences(), ctx.offsets(),
                    java.util.Arrays.asList(ctx.cigarsProvider().get()),
                    ctx.positions(), ctx.ownChromIds(), ctx.tagReferences());
            }
            return new EncodedChannel.DatasetBytes(
                SamTags.encode(((DecodedChannel.StrList) v).names(), tc));
        }
    }

    private static Map<Compression, Codec> build() {
        EnumMap<Compression, Codec> m = new EnumMap<>(Compression.class);
        m.put(Compression.RANS_ORDER0, new RansCodec(Compression.RANS_ORDER0, 0));
        m.put(Compression.RANS_ORDER1, new RansCodec(Compression.RANS_ORDER1, 1));
        m.put(Compression.BASE_PACK, new BasePackCodec());
        m.put(Compression.QUALITY_BINNED, new QualityCodec());
        m.put(Compression.DELTA_RANS_ORDER0, new DeltaRansCodec());
        m.put(Compression.NAME_TOKENIZED_V2, new NameTokenizedCodec());
        m.put(Compression.FQZCOMP_NX16_Z, new FqzcompCodec());
        m.put(Compression.MATE_INLINE_V2, new MateInfoCodec());
        m.put(Compression.REF_DIFF_V2, new RefDiffCodec());
        m.put(Compression.SAM_TAGS, new SamTagsCodec());
        return m;
    }
}
