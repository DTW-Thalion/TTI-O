/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.codecs;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * SAM_TAGS — SAM optional fields as a column codec (codec id 18, M101).
 *
 * <p>Per-read tag text (SAM columns 12+, tab-joined as samtools prints
 * them) coded as a tag-line dictionary plus one column per tag key,
 * with MD:Z and NM:i recomputed from the reference when they match.
 * The blob is written as {@code signal_channels/tags} with
 * {@code @compression = 18}. Spec: {@code docs/codecs/sam_tags.md}.</p>
 *
 * <p>Encode and decode take the same {@link Context}: the reads'
 * sequences, CIGARs, 1-based positions and chromosome ids, and one
 * reference sequence per chromosome id ({@code null} where there is
 * none). Without references MD and NM are stored like any other tag.
 * Delegates to {@code ttio_sam_tags_encode/_decode} through
 * {@link TtioRansNative}.</p>
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.codecs.sam_tags}, Objective-C {@code TTIOSamTags}.</p>
 */
public final class SamTags {

    private SamTags() {}

    /** @return {@code true} iff the native JNI library loaded. */
    public static boolean isAvailable() { return TtioRansNative.isAvailable(); }

    /** MD/NM derivation context.
     *
     *  @param sequences  concatenated read bases
     *  @param seqOffsets {@code n + 1} byte offsets into {@code sequences}
     *  @param cigars     one CIGAR per read
     *  @param positions  1-based POS per read (0 = unmapped)
     *  @param chromIds   chromosome id per read; {@code (short) 0xFFFF} = none
     *  @param references bases per chromosome id ({@code null} entries
     *                    allowed); {@code null} or all-null disables
     *                    derivation and the other fields are ignored */
    public record Context(byte[] sequences, long[] seqOffsets, List<String> cigars,
                          long[] positions, short[] chromIds, List<byte[]> references) {
        /** No derivation: MD and NM stored verbatim. */
        public static Context none() { return new Context(null, null, null, null, null, null); }

        boolean derives() {
            if (references == null) return false;
            for (byte[] r : references) if (r != null) return true;
            return false;
        }
    }

    /** Encode one blob's per-read tag text. */
    public static byte[] encode(List<String> tags, Context ctx) {
        Joined t = join(tags);
        Native c = Native.of(tags.size(), ctx);
        return TtioRansNative.encodeSamTags(t.bytes, t.offsets, c.seq, c.seqOffsets,
            c.cigars, c.cigarOffsets, c.positions, c.chromIds, c.refs);
    }

    /** Decode a blob to {@code nReads} per-read tag strings (same
     *  context as {@link #encode}). */
    public static List<String> decode(byte[] blob, int nReads, Context ctx) {
        Native c = Native.of(nReads, ctx);
        Object[] out = TtioRansNative.decodeSamTags(blob, nReads, c.seq, c.seqOffsets,
            c.cigars, c.cigarOffsets, c.positions, c.chromIds, c.refs);
        return split((byte[]) out[0], (long[]) out[1], nReads);
    }

    /** {@code flat[offsets[i] .. offsets[i+1])} as UTF-8 strings. */
    static List<String> split(byte[] flat, long[] offsets, int n) {
        List<String> res = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            int a = (int) offsets[i], b = (int) offsets[i + 1];
            res.add(new String(flat, a, b - a, StandardCharsets.UTF_8));
        }
        return res;
    }

    private record Joined(byte[] bytes, long[] offsets) {}

    private static Joined join(List<String> values) {
        int n = values.size();
        byte[][] enc = new byte[n][];
        long[] off = new long[n + 1];
        long total = 0;
        for (int i = 0; i < n; i++) {
            String v = values.get(i);
            enc[i] = v == null ? new byte[0] : v.getBytes(StandardCharsets.UTF_8);
            total += enc[i].length;
            off[i + 1] = total;
        }
        byte[] flat = new byte[(int) total];
        int p = 0;
        for (byte[] e : enc) {
            System.arraycopy(e, 0, flat, p, e.length);
            p += e.length;
        }
        return new Joined(flat, off);
    }

    /** The context marshalled for the JNI entry points. */
    private record Native(byte[] seq, long[] seqOffsets, byte[] cigars, long[] cigarOffsets,
                          long[] positions, short[] chromIds, byte[][] refs) {
        static Native of(int n, Context ctx) {
            if (ctx == null || !ctx.derives() || n == 0) {
                return new Native(null, null, null, null, null, null, null);
            }
            if ((ctx.sequences() == null || ctx.seqOffsets() == null
                    || ctx.cigars() == null || ctx.positions() == null
                    || ctx.chromIds() == null)) {
                throw new IllegalArgumentException(
                    "MD/NM derivation needs sequences, cigars, positions and chrom ids");
            }
            for (Object[] check : new Object[][]{
                    {"seqOffsets", (long) ctx.seqOffsets().length, n + 1L},
                    {"cigars", (long) ctx.cigars().size(), (long) n},
                    {"positions", (long) ctx.positions().length, (long) n},
                    {"chromIds", (long) ctx.chromIds().length, (long) n}}) {
                if (!check[1].equals(check[2])) {
                    throw new IllegalArgumentException(
                        check[0] + " must have " + check[2] + " entries, has " + check[1]);
                }
            }
            Joined cig = join(ctx.cigars());
            return new Native(ctx.sequences(), ctx.seqOffsets(), cig.bytes, cig.offsets,
                ctx.positions(), ctx.chromIds(), ctx.references().toArray(new byte[0][]));
        }
    }
}
