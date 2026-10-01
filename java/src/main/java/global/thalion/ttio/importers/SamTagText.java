/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.SAMRecord;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;

/**
 * SAM optional fields as {@code samtools view} prints them (M101).
 *
 * <p>Python and ObjC read a BAM/SAM through {@code samtools view} and
 * keep columns 12+ verbatim, so the Java importer must produce the same
 * text from an htsjdk record: fields in record order, tab-joined, each
 * {@code KEY:TYPE:VALUE} with the rules of htslib's
 * {@code sam_format_aux1}:</p>
 * <ul>
 *   <li>every integer type ({@code c C s S i I}) prints as {@code i}
 *       with its decimal value (unsigned types unsigned);</li>
 *   <li>{@code A}, {@code Z} and {@code H} print verbatim;</li>
 *   <li>{@code f} (and {@code d}) print as C {@code %g}
 *       ({@link #formatDouble});</li>
 *   <li>{@code B} arrays print {@code B:<subtype>} followed by
 *       {@code ,value} per element, integers in decimal and floats
 *       as C {@code %g}.</li>
 * </ul>
 *
 * <p>htsjdk keeps a record's attributes sorted by tag, so record order
 * is not recoverable from {@link SAMRecord#getAttributes()}: a BAM
 * record is formatted from its raw aux bytes
 * ({@link BAMRecord#getVariableBinaryRepresentation()}), a SAM text
 * record from the line's own columns ({@link #normalizeSamText}). A
 * CRAM record has neither, so it falls back to the sorted attributes
 * ({@link #fromAttributes}), where an {@code H} value reads as
 * {@code Z}.</p>
 */
public final class SamTagText {

    private SamTagText() {}

    /** The tag text of {@code rec}: from its raw aux bytes when it is an
     *  unmodified {@link BAMRecord}, else from its attributes. */
    public static String fromRecord(SAMRecord rec) {
        if (rec instanceof BAMRecord b) {
            String s = fromBamRecord(b);
            if (s != null) return s;
        }
        return fromAttributes(rec);
    }

    /** Format a BAM record's aux block in record order; {@code null}
     *  when the raw bytes are no longer available. */
    static String fromBamRecord(BAMRecord rec) {
        byte[] raw = rec.getVariableBinaryRepresentation();
        int auxLen = rec.getAttributesBinarySize();
        if (raw == null || auxLen < 0) return null;
        int start = raw.length - auxLen;
        // htslib moves a long CIGAR out of the CG tag into the record
        // (bam_tag2cigar) when the stored CIGAR is the <l_seq>S<n>N
        // placeholder, so samtools does not print that CG tag.
        boolean dropCg = false;
        if (rec.getCigarLength() == 2) {
            int nameLen = rec.getReadNameLength() + 1;
            ByteBuffer cb = ByteBuffer.wrap(raw, nameLen, 8).order(ByteOrder.LITTLE_ENDIAN);
            int op0 = cb.getInt(), op1 = cb.getInt();
            dropCg = (op0 & 0xF) == 4 && (op0 >>> 4) == rec.getReadLength() && (op1 & 0xF) == 3;
        }
        return formatAux(raw, start, raw.length, dropCg);
    }

    /** Format a BAM aux block {@code buf[from, to)} as samtools prints it. */
    public static String formatAux(byte[] buf, int from, int to, boolean dropCg) {
        ByteBuffer bb = ByteBuffer.wrap(buf, from, to - from).order(ByteOrder.LITTLE_ENDIAN);
        StringBuilder sb = new StringBuilder();
        while (bb.remaining() >= 3) {
            char k0 = (char) (bb.get() & 0xFF), k1 = (char) (bb.get() & 0xFF);
            char type = (char) (bb.get() & 0xFF);
            int mark = sb.length();
            if (sb.length() > 0) sb.append('\t');
            sb.append(k0).append(k1).append(':');
            switch (type) {
                case 'A' -> sb.append("A:").append((char) (bb.get() & 0xFF));
                case 'c' -> sb.append("i:").append(bb.get());
                case 'C' -> sb.append("i:").append(bb.get() & 0xFF);
                case 's' -> sb.append("i:").append(bb.getShort());
                case 'S' -> sb.append("i:").append(bb.getShort() & 0xFFFF);
                case 'i' -> sb.append("i:").append(bb.getInt());
                case 'I' -> sb.append("i:").append(bb.getInt() & 0xFFFFFFFFL);
                case 'f' -> sb.append("f:").append(formatDouble(bb.getFloat()));
                case 'd' -> sb.append("d:").append(formatDouble(bb.getDouble()));
                case 'Z', 'H' -> {
                    int s = bb.position();
                    int e = s;
                    while (e < from + (to - from) && buf[e] != 0) e++;
                    sb.append(type).append(':')
                      .append(new String(buf, s, e - s, StandardCharsets.UTF_8));
                    bb.position(Math.min(e + 1, to));
                }
                case 'B' -> {
                    char sub = (char) (bb.get() & 0xFF);
                    int n = bb.getInt();
                    sb.append("B:").append(sub);
                    for (int i = 0; i < n; i++) {
                        sb.append(',');
                        switch (sub) {
                            case 'c' -> sb.append(bb.get());
                            case 'C' -> sb.append(bb.get() & 0xFF);
                            case 's' -> sb.append(bb.getShort());
                            case 'S' -> sb.append(bb.getShort() & 0xFFFF);
                            case 'i' -> sb.append(bb.getInt());
                            case 'I' -> sb.append(bb.getInt() & 0xFFFFFFFFL);
                            case 'f' -> sb.append(formatDouble(bb.getFloat()));
                            default -> throw new IllegalArgumentException(
                                "BAM aux tag " + k0 + k1 + ": unknown B subtype '" + sub + "'");
                        }
                    }
                }
                default -> throw new IllegalArgumentException(
                    "BAM aux tag " + k0 + k1 + ": unknown type '" + type + "'");
            }
            if (dropCg && k0 == 'C' && k1 == 'G' && type == 'B') sb.setLength(mark);
        }
        return sb.toString();
    }

    /** The tag text of {@code rec} from its htsjdk attributes, in the
     *  (tag-sorted) order htsjdk keeps them. */
    public static String fromAttributes(SAMRecord rec) {
        StringBuilder sb = new StringBuilder();
        for (SAMRecord.SAMTagAndValue tv : rec.getAttributes()) {
            if (sb.length() > 0) sb.append('\t');
            sb.append(tv.tag).append(':');
            Object v = tv.value;
            boolean unsigned = rec.isUnsignedArrayAttribute(tv.tag);
            if (v instanceof String s) {
                sb.append("Z:").append(s);
            } else if (v instanceof Character c) {
                sb.append("A:").append(c.charValue());
            } else if (v instanceof Float f) {
                sb.append("f:").append(formatDouble(f));
            } else if (v instanceof Double d) {
                sb.append("d:").append(formatDouble(d));
            } else if (v instanceof Number num) {
                sb.append("i:").append(num.longValue());
            } else if (v instanceof byte[] a) {
                sb.append("B:").append(unsigned ? 'C' : 'c');
                for (byte x : a) sb.append(',').append(unsigned ? x & 0xFF : x);
            } else if (v instanceof short[] a) {
                sb.append("B:").append(unsigned ? 'S' : 's');
                for (short x : a) sb.append(',').append(unsigned ? x & 0xFFFF : x);
            } else if (v instanceof int[] a) {
                sb.append("B:").append(unsigned ? 'I' : 'i');
                for (int x : a) sb.append(',').append(unsigned ? x & 0xFFFFFFFFL : x);
            } else if (v instanceof float[] a) {
                sb.append("B:f");
                for (float x : a) sb.append(',').append(formatDouble(x));
            } else {
                throw new IllegalArgumentException(
                    "SAM tag " + tv.tag + ": unsupported value type " + v.getClass().getName());
            }
        }
        return sb.toString();
    }

    /** Normalise SAM text optional fields (columns 12+, tab-joined) the
     *  way {@code samtools view} re-prints them after parsing: integers
     *  in canonical decimal, floats as C {@code %g}, {@code A /
     *  Z / H} verbatim, {@code B} arrays element by element. */
    public static String normalizeSamText(String fields) {
        if (fields == null || fields.isEmpty()) return "";
        StringBuilder sb = new StringBuilder(fields.length());
        for (String f : fields.split("\t", -1)) {
            if (sb.length() > 0) sb.append('\t');
            if (f.length() < 5 || f.charAt(2) != ':' || f.charAt(4) != ':') {
                throw new IllegalArgumentException("malformed SAM optional field '" + f + "'");
            }
            String key = f.substring(0, 2);
            char type = f.charAt(3);
            String val = f.substring(5);
            sb.append(key).append(':');
            switch (type) {
                case 'A', 'a', 'c', 'C' -> sb.append("A:").append(val);
                case 'i', 'I' -> sb.append("i:").append(Long.parseLong(stripPlus(val.trim())));
                case 'f' -> sb.append("f:").append(formatDouble((float) Double.parseDouble(val.trim())));
                case 'd' -> sb.append("d:").append(formatDouble(Double.parseDouble(val.trim())));
                case 'Z', 'H' -> sb.append(type).append(':').append(val);
                case 'B' -> {
                    String[] parts = val.split(",", -1);
                    char sub = parts[0].charAt(0);
                    sb.append("B:").append(sub);
                    for (int i = 1; i < parts.length; i++) {
                        String p = parts[i].trim();
                        sb.append(',');
                        if (sub == 'f') {
                            sb.append(formatDouble((float) Double.parseDouble(p)));
                        } else {
                            sb.append(Long.parseLong(stripPlus(p)));
                        }
                    }
                }
                default -> throw new IllegalArgumentException(
                    "SAM optional field '" + f + "': unknown type '" + type + "'");
            }
        }
        return sb.toString();
    }

    private static String stripPlus(String s) {
        return s.startsWith("+") ? s.substring(1) : s;
    }

    /** A float tag value ({@code f}, {@code d}, {@code B:f} element) as
     *  samtools 1.22 prints it: C {@code %g} (six significant digits,
     *  glibc rounding of the exact binary value, so ties go to even:
     *  755088.5 prints {@code 755088}); zero prints {@code 0}, or
     *  {@code -0} with the sign bit. Verified against
     *  {@code samtools view} over a seeded sweep of values. */
    public static String formatDouble(double d) {
        if (d == 0) return (Double.doubleToRawLongBits(d) < 0) ? "-0" : "0";
        if (Double.isNaN(d)) return cG(d);
        return d < 0 ? "-" + cG(-d) : cG(d);
    }

    /** C {@code printf("%g", d)} for finite positive {@code d} (glibc
     *  rounds the exact binary value), {@code inf} / {@code nan}
     *  otherwise. */
    static String cG(double d) {
        if (Double.isNaN(d)) return Double.doubleToRawLongBits(d) < 0 ? "-nan" : "nan";
        if (Double.isInfinite(d)) return d < 0 ? "-inf" : "inf";
        BigDecimal r = new BigDecimal(d).round(new MathContext(6, RoundingMode.HALF_EVEN));
        String digits = r.unscaledValue().abs().toString();
        int x = digits.length() - 1 - r.scale();          // decimal exponent
        // Pad the significand to 6 digits.
        StringBuilder sig = new StringBuilder(digits);
        while (sig.length() < 6) sig.append('0');
        if (x < -4 || x >= 6) {
            String mant = sig.charAt(0) + "." + sig.substring(1);
            mant = stripZeros(mant);
            String exp = Integer.toString(Math.abs(x));
            if (exp.length() < 2) exp = "0" + exp;
            return mant + "e" + (x < 0 ? "-" : "+") + exp;
        }
        String fixed = r.setScale(Math.max(5 - x, 0), RoundingMode.HALF_EVEN).toPlainString();
        return stripZeros(fixed);
    }

    private static String stripZeros(String s) {
        if (s.indexOf('.') < 0) return s;
        int e = s.length();
        while (s.charAt(e - 1) == '0') e--;
        if (s.charAt(e - 1) == '.') e--;
        return s.substring(0, e);
    }
}
