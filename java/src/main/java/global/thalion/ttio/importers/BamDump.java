/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import global.thalion.ttio.genomics.WrittenGenomicRun;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.List;
import java.util.TreeMap;

/**
 * BamDump — canonical-JSON dump of a SAM/BAM file for the M87
 * cross-language conformance harness.
 *
 * <p>Usage:
 * <pre>
 *   java -cp ... global.thalion.ttio.importers.BamDump &lt;bam_or_sam_path&gt; [--name NAME]
 * </pre>
 *
 * <p>Reads the file via {@link BamReader} and emits a canonical JSON
 * document on stdout matching the schema documented in
 * {@code HANDOFF.md} §7. The output is byte-identical to the Python
 * {@code python -m ttio.importers.bam_dump} reference.</p>
 *
 * <p>Format: keys sorted alphabetically; 2-space indent; each list
 * element on its own indented line (matches Python
 * {@code json.dumps(sort_keys=True, indent=2)}); trailing newline.</p>
 *
 * (M87)
 */
public final class BamDump {

    private BamDump() {}

    public static int run(String[] args, Writer out) throws IOException {
        String pathStr = null;
        String name = "genomic_0001";
        String reference = null;
        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if ("--name".equals(a) && i + 1 < args.length) {
                name = args[++i];
            } else if (a.startsWith("--name=")) {
                name = a.substring("--name=".length());
            } else if ("--reference".equals(a) && i + 1 < args.length) {
                reference = args[++i];
            } else if (a.startsWith("--reference=")) {
                reference = a.substring("--reference=".length());
            } else if ("-h".equals(a) || "--help".equals(a)) {
                out.write("Usage: BamDump <path> [--reference FASTA] [--name NAME]\n");
                return 0;
            } else if (pathStr == null) {
                pathStr = a;
            } else {
                throw new IllegalArgumentException(
                    "Unexpected positional argument: " + a);
            }
        }
        if (pathStr == null) {
            throw new IllegalArgumentException(
                "Usage: BamDump <path> [--reference FASTA] [--name NAME]");
        }

        // dispatch on .cram extension to CramReader. For
        // BAM/SAM paths --reference is accepted but unused.
        BamReader reader;
        if (pathStr.toLowerCase().endsWith(".cram")) {
            if (reference == null) {
                System.err.println(
                    "error: --reference <fasta> is required for .cram input");
                return 2;
            }
            reader = new CramReader(Paths.get(pathStr), Paths.get(reference));
        } else {
            reader = new BamReader(Paths.get(pathStr));
        }
        WrittenGenomicRun run = reader.toGenomicRun(name);

        TreeMap<String, Object> payload = new TreeMap<>();
        payload.put("name", name);
        payload.put("read_count", run.readNames().size());
        payload.put("sample_name", run.sampleName());
        payload.put("platform", run.platform());
        payload.put("reference_uri", run.referenceUri());
        payload.put("read_names", run.readNames());
        payload.put("positions", longArray(run.positions()));
        payload.put("chromosomes", run.chromosomes());
        payload.put("flags", intArray(run.flags()));
        payload.put("mapping_qualities", byteArrayUnsigned(run.mappingQualities()));
        payload.put("cigars", run.cigars());
        payload.put("mate_chromosomes", run.mateChromosomes());
        payload.put("mate_positions", longArray(run.matePositions()));
        payload.put("template_lengths", intArray(run.templateLengths()));
        payload.put("sequences_md5", md5Hex(run.sequences()));
        payload.put("qualities_md5", md5Hex(run.qualities()));
        payload.put("provenance_count", reader.lastProvenance().size());
        // M101: the SAM optional fields of each read, as samtools
        // prints them ("" when a read has none).
        payload.put("tags", run.tags() != null ? run.tags()
            : java.util.Collections.nCopies(run.readNames().size(), ""));

        StringBuilder sb = new StringBuilder(4096);
        writeJson(sb, payload, 0);
        sb.append('\n');
        out.write(sb.toString());
        out.flush();
        return 0;
    }

    public static void main(String[] args) throws IOException {
        Writer w = new OutputStreamWriter(System.out, StandardCharsets.UTF_8);
        int code = run(args, w);
        w.flush();
        if (code != 0) System.exit(code);
    }

    // ── Helpers ─────────────────────────────────────────────────────

    private static List<Long> longArray(long[] a) {
        List<Long> out = new java.util.ArrayList<>(a.length);
        for (long v : a) out.add(v);
        return out;
    }

    private static List<Integer> intArray(int[] a) {
        List<Integer> out = new java.util.ArrayList<>(a.length);
        for (int v : a) out.add(v);
        return out;
    }

    private static List<Integer> byteArrayUnsigned(byte[] a) {
        List<Integer> out = new java.util.ArrayList<>(a.length);
        for (byte v : a) out.add(v & 0xFF);
        return out;
    }

    static String md5Hex(byte[] data) {
        try {
            MessageDigest md = MessageDigest.getInstance("MD5");
            byte[] digest = md.digest(data);
            StringBuilder sb = new StringBuilder(digest.length * 2);
            for (byte b : digest) {
                sb.append(String.format("%02x", b & 0xFF));
            }
            return sb.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("MD5 not available", e);
        }
    }

    /**
     * Hand-rolled JSON emitter matching Python's
     * {@code json.dumps(sort_keys=True, indent=2)} byte-for-byte.
     *
     * <p>Object: each key/value on its own line, trailing key has no
     * comma, closing brace dedented; empty object renders as
     * {@code "{}"}. Array: each element on its own line, trailing
     * element has no comma, closing bracket dedented; empty array
     * renders as {@code "[]"}. Both objects and arrays open with the
     * brace immediately, then a newline, then indented children.</p>
     */
    @SuppressWarnings("unchecked")
    static void writeJson(StringBuilder sb, Object value, int indent) {
        if (value == null) {
            sb.append("null");
        } else if (value instanceof Boolean b) {
            sb.append(b.booleanValue() ? "true" : "false");
        } else if (value instanceof Number) {
            sb.append(value.toString());
        } else if (value instanceof CharSequence) {
            writeJsonString(sb, value.toString());
        } else if (value instanceof java.util.Map<?, ?> m) {
            // TreeMap iteration is sorted; LinkedHashMap is insertion
            // order. For the canonical dump we always pass a TreeMap.
            if (m.isEmpty()) {
                sb.append("{}");
                return;
            }
            sb.append('{');
            boolean first = true;
            for (var entry : m.entrySet()) {
                if (!first) sb.append(',');
                first = false;
                sb.append('\n');
                appendIndent(sb, indent + 1);
                writeJsonString(sb, entry.getKey().toString());
                sb.append(": ");
                writeJson(sb, entry.getValue(), indent + 1);
            }
            sb.append('\n');
            appendIndent(sb, indent);
            sb.append('}');
        } else if (value instanceof java.util.List<?> l) {
            if (l.isEmpty()) {
                sb.append("[]");
                return;
            }
            sb.append('[');
            boolean first = true;
            for (Object item : l) {
                if (!first) sb.append(',');
                first = false;
                sb.append('\n');
                appendIndent(sb, indent + 1);
                writeJson(sb, item, indent + 1);
            }
            sb.append('\n');
            appendIndent(sb, indent);
            sb.append(']');
        } else {
            // Fallback: stringify.
            writeJsonString(sb, value.toString());
        }
    }

    private static void appendIndent(StringBuilder sb, int level) {
        for (int i = 0; i < level; i++) sb.append("  ");
    }

    /**
     * Emit a JSON string literal matching Python's default
     * {@code json.dumps} encoding: ASCII output (non-ASCII codepoints
     * escaped as {@code \\uXXXX}), control characters escaped, plus
     * {@code "}, {@code \\}, and the C0 escapes for newline / tab /
     * etc.
     */
    private static void writeJsonString(StringBuilder sb, String s) {
        sb.append('"');
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '\\': sb.append("\\\\"); break;
                case '"':  sb.append("\\\""); break;
                case '\b': sb.append("\\b"); break;
                case '\f': sb.append("\\f"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < 0x20 || c > 0x7F) {
                        sb.append(String.format("\\u%04x", (int) c));
                    } else {
                        sb.append(c);
                    }
            }
        }
        sb.append('"');
    }
}
