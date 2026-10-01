/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.codecs.SamTags;
import global.thalion.ttio.exporters.BamWriter;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.importers.BamReader;
import global.thalion.ttio.importers.ImporterRegistry;
import global.thalion.ttio.importers.SamTagText;
import global.thalion.ttio.protection.EncryptedTransport;
import global.thalion.ttio.protection.PerAUFile;
import global.thalion.ttio.protection.SignatureManager;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import global.thalion.ttio.transport.AccessUnit;
import global.thalion.ttio.transport.ChannelData;
import global.thalion.ttio.transport.PacketType;
import global.thalion.ttio.transport.TransportReader;
import global.thalion.ttio.transport.TransportWriter;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.junit.jupiter.api.Assertions.*;

/** M101: SAM optional tags through import, storage, read and export.
 *  Mirrors python/tests/test_m101_sam_tags.py. The SAM fixture is built
 *  in the test (random reference; reads with mismatches, deletions,
 *  insertions, soft clips and an unmapped read) with MD/NM written by
 *  {@code samtools calmd} as the independent ground truth. */
class M101SamTagsTest {

    @TempDir Path tmp;

    private static final int REF_LEN = 20_000;
    private static final String RUN = "genomic_0001";

    // ── samtools helpers ─────────────────────────────────────────────

    private static boolean samtools() {
        try {
            Process p = new ProcessBuilder("samtools", "--version")
                .redirectErrorStream(true).start();
            p.getInputStream().readAllBytes();
            return p.waitFor() == 0;
        } catch (IOException | InterruptedException e) {
            return false;
        }
    }

    private static void needs() {
        Assumptions.assumeTrue(SamTags.isAvailable(), "libttio_rans_jni");
        Assumptions.assumeTrue(samtools(), "samtools on PATH");
    }

    private static byte[] run(Path stdout, String... cmd) throws Exception {
        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectError(ProcessBuilder.Redirect.DISCARD);
        if (stdout != null) pb.redirectOutput(stdout.toFile());
        Process p = pb.start();
        byte[] out = stdout == null ? p.getInputStream().readAllBytes() : new byte[0];
        assertEquals(0, p.waitFor(), String.join(" ", cmd));
        return out;
    }

    /** {@code samtools view} alignment lines (no header). */
    private static List<String> view(Path p) throws Exception {
        String s = new String(run(null, "samtools", "view", p.toString()), StandardCharsets.UTF_8);
        List<String> out = new ArrayList<>();
        for (String line : s.split("\n")) if (!line.isEmpty()) out.add(line);
        return out;
    }

    /** Columns 12+ of each alignment line, as samtools prints them. */
    private static List<String> viewTags(Path p) throws Exception {
        List<String> out = new ArrayList<>();
        for (String line : view(p)) {
            String[] c = line.split("\t", 12);
            out.add(c.length > 11 ? c[11] : "");
        }
        return out;
    }

    /** Write ref.fa and in.sam (coordinate-sorted, one chromosome);
     *  the first {@code taglessPrefix} reads carry no optional fields. */
    private Path[] makeSam(Path dir, int n, boolean tags, int taglessPrefix, long seed)
            throws Exception {
        Files.createDirectories(dir);
        Random rng = new Random(seed);
        StringBuilder rs = new StringBuilder(REF_LEN);
        for (int i = 0; i < REF_LEN; i++) rs.append("ACGT".charAt(rng.nextInt(4)));
        String ref = rs.toString();
        Path fa = dir.resolve("ref.fa");
        StringBuilder fs = new StringBuilder(">chr1\n");
        for (int i = 0; i < REF_LEN; i += 60) fs.append(ref, i, Math.min(i + 60, REF_LEN)).append('\n');
        Files.writeString(fa, fs.toString());
        run(null, "samtools", "faidx", fa.toString());

        List<String> lines = new ArrayList<>(List.of("@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:" + REF_LEN, "@RG\tID:g1\tSM:S1\tPL:ILLUMINA"));
        String comp = "ACGT", sub = "CGTA";
        for (int k = 0; k <= n; k++) {
            List<String> r = new ArrayList<>();
            if (k == n) {
                r.addAll(List.of("u0", "4", "*", "0", "0", "*", "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII"));
            } else {
                int pos = 1 + (k * (REF_LEN - 400)) / n;
                StringBuilder q = new StringBuilder();
                for (int j = 0; j < 100; j++) q.append("#+5:?FI".charAt(rng.nextInt(7)));
                String seq = ref.substring(pos - 1, pos - 1 + 100);
                String cigar = "100M";
                switch (k % 6) {
                    case 1 -> {
                        char[] s = seq.toCharArray();
                        for (int m = 0; m < 3; m++) {
                            int j = rng.nextInt(100);
                            s[j] = sub.charAt(comp.indexOf(s[j]));
                        }
                        seq = new String(s);
                    }
                    case 2 -> {
                        seq = ref.substring(pos - 1, pos + 39) + ref.substring(pos + 44, pos + 104);
                        cigar = "40M5D60M";
                    }
                    case 3 -> {
                        seq = ref.substring(pos - 1, pos + 49) + "TTT" + ref.substring(pos + 49, pos + 96);
                        cigar = "50M3I47M";
                    }
                    case 4 -> {
                        seq = "GGGGG" + ref.substring(pos - 1, pos + 94);
                        cigar = "5S95M";
                    }
                    default -> { }
                }
                int flag = k % 2 == 0 ? 99 : 147;
                r.addAll(List.of(String.format("r%05d", k), Integer.toString(flag), "chr1",
                    Integer.toString(pos), "60", cigar, "=", Integer.toString(pos + 150),
                    flag == 99 ? "250" : "-250", seq, q.toString()));
            }
            if (tags && k >= taglessPrefix) {
                int a = rng.nextInt(300);
                r.addAll(List.of("AS:i:" + a, "UQ:i:" + a, "RG:Z:g1", "XS:i:" + (rng.nextInt(10) - 5)));
                if (k % 7 == 0) {
                    r.add("XA:Z:chr1,+" + (1 + rng.nextInt(9998)) + ",100M," + rng.nextInt(4) + ";");
                }
                if (k % 11 == 0) r.add("BC:B:c,-3,0,7");
                if (k % 13 == 0) r.add("XF:f:0.25");
            }
            lines.add(String.join("\t", r));
        }
        Path raw = dir.resolve("raw.sam");
        Files.writeString(raw, String.join("\n", lines) + "\n");
        Path sam = dir.resolve("in.sam");
        if (!tags) {
            Files.copy(raw, sam);
            return new Path[]{ sam, fa };
        }
        // calmd appends NM/MD where they are missing (all tagged reads).
        run(sam, "samtools", "calmd", raw.toString(), fa.toString());
        if (taglessPrefix > 0) {
            List<String> out = new ArrayList<>();
            for (String line : Files.readAllLines(sam)) {
                String[] c = line.split("\t");
                if (!line.startsWith("@") && c[0].startsWith("r")
                        && Integer.parseInt(c[0].substring(1)) < taglessPrefix) {
                    line = String.join("\t", java.util.Arrays.copyOf(c, 11));
                }
                out.add(line);
            }
            Files.writeString(sam, String.join("\n", out) + "\n");
        }
        return new Path[]{ sam, fa };
    }

    private Path importSam(Path sam, Path ref, Integer blockReads, String name) throws Exception {
        Path out = tmp.resolve(name);
        Map<String, Object> opts = new LinkedHashMap<>();
        if (ref != null) {
            opts.put("reference", ref.toString());
            opts.put("embed_reference", true);
        }
        if (blockReads != null) opts.put("block_reads", blockReads);
        ImporterRegistry.encode("sam", List.of(sam.toString()), out, opts);
        return out;
    }

    private static List<String> storedTags(Path tio) {
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            GenomicRun g = ds.genomicRuns().values().iterator().next();
            List<String> out = new ArrayList<>();
            var it = g.iterReads();
            while (it.hasNext()) out.add(it.next().tags());
            return out;
        }
    }

    /** {@code samtools view} lines of the run exported to BAM. */
    private List<String> exported(Path tio, String name) throws Exception {
        Path bam = tmp.resolve(name);
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            GenomicRun g = ds.genomicRuns().values().iterator().next();
            new BamWriter(bam).write(g, List.of(), false);
        }
        return view(bam);
    }

    private static boolean hasFlag(Path tio, String flag) {
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            return FeatureFlags.readFrom(sp.rootGroup()).features().contains(flag);
        }
    }

    private static List<Map<String, Object>> indexRows(Path tio) {
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            StorageGroup rg = sp.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup(RUN);
            try (StorageDataset ds = rg.openGroup("blocks").openDataset("index")) {
                return ds.readRows();
            }
        }
    }

    /** (index rows, tags blob, sequences blob, qualities blob). */
    private static List<Object> channelState(Path tio) {
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            StorageGroup rg = sp.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup(RUN);
            StorageGroup sc = rg.openGroup("signal_channels");
            List<Object> out = new ArrayList<>();
            try (StorageDataset ds = rg.openGroup("blocks").openDataset("index")) {
                out.add(ds.readRows());
            }
            try (StorageDataset ds = sc.openDataset("tags")) { out.add(List.of((byte[]) ds.readAll())); }
            try (StorageDataset ds = sc.openGroup("sequences").openDataset("data")) {
                out.add(List.of((byte[]) ds.readAll()));
            }
            try (StorageDataset ds = sc.openDataset("qualities")) { out.add(List.of((byte[]) ds.readAll())); }
            return out;
        }
    }

    private static void assertStateEquals(List<Object> a, List<Object> b, String what) {
        assertEquals(a.get(0), b.get(0), what + ": blocks/index");
        for (int i = 1; i < a.size(); i++) {
            assertArrayEquals((byte[]) ((List<?>) a.get(i)).get(0),
                (byte[]) ((List<?>) b.get(i)).get(0), what + ": channel " + i);
        }
    }

    // ── codec ────────────────────────────────────────────────────────

    @Test
    void codecRoundTripWithDerivation() {
        Assumptions.assumeTrue(SamTags.isAvailable(), "libttio_rans_jni");
        Random rng = new Random(5);
        byte[] ref = new byte[4000];
        for (int i = 0; i < ref.length; i++) ref[i] = (byte) "ACGT".charAt(rng.nextInt(4));
        int n = 50;
        byte[] seq = new byte[n * 50];
        long[] off = new long[n + 1];
        long[] pos = new long[n];
        short[] chr = new short[n];
        List<String> cig = new ArrayList<>();
        List<String> tags = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            pos[i] = 1 + i * 70L;
            System.arraycopy(ref, (int) pos[i] - 1, seq, i * 50, 50);
            int mm = 10 + i % 30;
            byte orig = seq[i * 50 + mm];
            seq[i * 50 + mm] = (byte) (orig == 'A' ? 'C' : 'A');
            off[i + 1] = off[i] + 50;
            cig.add("50M");
            String md = mm + String.valueOf((char) orig) + (49 - mm);
            tags.add("NM:i:1\tMD:Z:" + md + "\tAS:i:" + (40 + i % 5) + "\tXS:i:" + (40 + i % 5));
        }
        SamTags.Context ctx = new SamTags.Context(seq, off, cig, pos, chr, List.of(ref));
        byte[] derived = SamTags.encode(tags, ctx);
        byte[] verbatim = SamTags.encode(tags, SamTags.Context.none());
        assertEquals(tags, SamTags.decode(derived, n, ctx));
        assertEquals(tags, SamTags.decode(verbatim, n, SamTags.Context.none()));
        assertTrue(derived.length < verbatim.length,
            "MD/NM recomputed: " + derived.length + " < " + verbatim.length);
        // Empty strings and non-canonical text round-trip verbatim.
        List<String> odd = List.of("", "zz:Z:not a tag", "AS:i:007", "é:Z:x");
        assertEquals(odd, SamTags.decode(SamTags.encode(odd, SamTags.Context.none()), 4,
            SamTags.Context.none()));
    }

    // ── samtools-text formatter (item 3) ─────────────────────────────

    @Test
    void formatDoubleMatchesKnownHtslibOutput() {
        assertEquals("0.25", SamTagText.formatDouble(0.25f));
        assertEquals("0.1", SamTagText.formatDouble(0.1f));
        assertEquals("3.14159", SamTagText.formatDouble(3.14159274f));
        assertEquals("-2.5", SamTagText.formatDouble(-2.5f));
        assertEquals("123457", SamTagText.formatDouble(123456.7f));
        assertEquals("1e-05", SamTagText.formatDouble(1e-05f));
        assertEquals("1e+06", SamTagText.formatDouble(1e6f));
        assertEquals("1.23457e+07", SamTagText.formatDouble(12345678f));
        assertEquals("0", SamTagText.formatDouble(0f));
        assertEquals("-0", SamTagText.formatDouble(-0f));
        assertEquals("755088", SamTagText.formatDouble(755088.5f));
        assertEquals("XF:f:0.25\tXI:i:4294967295\tXB:B:f,0.1,-2.5",
            SamTagText.normalizeSamText("XF:f:0.250\tXI:i:4294967295\tXB:B:f,0.1,-2.50"));
    }

    /** Every tag type, BAM and SAM input, against {@code samtools view}:
     *  record order kept, integers as {@code i}, {@code A/Z/H} verbatim,
     *  floats and {@code B} arrays formatted as htslib does (plus a
     *  seeded sweep of float values across the whole range). */
    @Test
    void importerTagTextEqualsSamtoolsView() throws Exception {
        needs();
        List<String> lines = new ArrayList<>(List.of("@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:1000"));
        String base = "\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII";
        lines.add("t0" + base + "\tZZ:A:x\tXC:i:-128\tXD:i:255\tXE:i:-32768\tXG:i:65535"
            + "\tXH:i:-2147483648\tXI:i:4294967295\tXJ:i:2147483648\tXK:i:0"
            + "\tXF:f:0.25\tXZ:Z:hello world\tXX:H:1AE301\tAB:B:c,-1,2,-3\tXU:B:C,0,255"
            + "\tXS:B:s,-300,300\tXT:B:S,65535\tXV:B:i,-70000,70000\tXW:B:I,4294967295"
            + "\tXY:B:f,0.1,-2.5,1e-05,123456.7\tXQ:B:c");
        lines.add("t1" + base);
        lines.add("t2" + base + "\tNM:i:+3\tXF:f:1.50\tAA:i:007");
        // A rounding tie: %g rounds the exact binary value half to even.
        lines.add("t3" + base + "\tXP:f:755088.5\tXR:B:f,755088.5,-0,2.5e-05");
        Random rng = new Random(101);
        for (int r = 0; r < 3000; r++) {
            StringBuilder sb = new StringBuilder("f" + r + base);
            float v;
            if (r % 3 == 0) {
                do { v = Float.intBitsToFloat(rng.nextInt()); } while (!Float.isFinite(v));
            } else {
                v = (float) (Math.pow(10, rng.nextDouble() * 14 - 6) * (rng.nextBoolean() ? 1 : -1));
            }
            sb.append("\tFA:f:").append(v).append("\tFB:B:f");
            for (int k = 0; k < 6; k++) {
                float w = (float) (Math.pow(10, rng.nextDouble() * 16 - 7) * (rng.nextBoolean() ? 1 : -1));
                sb.append(',').append(w);
            }
            lines.add(sb.toString());
        }
        Path sam = tmp.resolve("types.sam");
        Files.writeString(sam, String.join("\n", lines) + "\n");
        Path bam = tmp.resolve("types.bam");
        run(null, "samtools", "view", "-b", "-o", bam.toString(), sam.toString());

        List<String> want = viewTags(sam);
        assertEquals(want, viewTags(bam), "samtools prints SAM and BAM alike");
        assertSameTags(want, new BamReader(sam).toGenomicRun("s").tags(), "SAM input");
        assertSameTags(want, new BamReader(bam).toGenomicRun("b").tags(), "BAM input");
        assertEquals("XP:f:755088\tXR:B:f,755088,-0,2.5e-05", want.get(3));
        // H is kept as H (htsjdk would read it as Z).
        assertTrue(want.get(0).contains("\tXX:H:1AE301\t"));

        // Streaming batches carry the same text.
        List<String> streamed = new ArrayList<>();
        var it = new BamReader(bam).iterBatches("b", null, null, 7);
        while (it.hasNext()) streamed.addAll(it.next().tags());
        assertEquals(want, streamed);
    }

    /** Equal lists, reporting the first differing field. */
    private static void assertSameTags(List<String> want, List<String> got, String what) {
        assertEquals(want.size(), got.size(), what + ": read count");
        for (int i = 0; i < want.size(); i++) {
            String[] a = want.get(i).split("\t", -1), b = got.get(i).split("\t", -1);
            for (int k = 0; k < Math.min(a.length, b.length); k++) {
                assertEquals(a[k], b[k], what + ": read " + i + " field " + k);
            }
            assertEquals(want.get(i), got.get(i), what + ": read " + i);
        }
    }

    @Test
    void bamDumpCarriesTags() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("dump"), 30, true, 0, 3);
        java.io.StringWriter w = new java.io.StringWriter();
        assertEquals(0, global.thalion.ttio.importers.BamDump.run(new String[]{ f[0].toString() }, w));
        String json = w.toString();
        assertTrue(json.contains("\"tags\": ["), json);
        assertTrue(json.contains("\"AS:i:"), json);
    }

    // ── import / store / read / export ───────────────────────────────

    @Test
    void importWithReferenceDerivesMdNm() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("a"), 600, true, 0, 11);
        Path tio = importSam(f[0], f[1], 200, "a.tio");
        assertEquals(viewTags(f[0]), storedTags(tio));
        assertTrue(hasFlag(tio, FeatureFlags.OPT_SAM_TAGS));
        List<Map<String, Object>> rows = indexRows(tio);
        assertTrue(rows.get(0).containsKey("tags_off"));
        for (Map<String, Object> r : rows) {
            if (((Number) r.get("tags_len")).longValue() > 0) {
                assertEquals(Compression.SAM_TAGS.ordinal(), ((Number) r.get("tags_codec")).intValue());
            }
        }
        long stored = 0;
        for (Map<String, Object> r : rows) stored += ((Number) r.get("tags_len")).longValue();
        long text = 0;
        for (String t : viewTags(f[0])) text += t.length();
        assertTrue(stored * 8 < text, "MD/NM recomputed: " + stored + " vs " + text);
        // Full records survive import -> write -> read -> export.
        assertEquals(view(f[0]), exported(tio, "a.bam"));
    }

    @Test
    void importWithoutReference() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("b"), 600, true, 0, 11);
        Path tio = importSam(f[0], null, 250, "b.tio");
        assertEquals(viewTags(f[0]), storedTags(tio));
        assertEquals(view(f[0]), exported(tio, "b.bam"));
    }

    @Test
    void taglessInputWritesNoTagsChannel() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("c"), 300, false, 0, 11);
        Path tio = importSam(f[0], f[1], 200, "c.tio");
        assertFalse(hasFlag(tio, FeatureFlags.OPT_SAM_TAGS));
        assertFalse(indexRows(tio).get(0).containsKey("tags_off"));
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            assertFalse(sp.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup(RUN)
                .openGroup("signal_channels").hasChild("tags"));
        }
        for (String t : storedTags(tio)) assertEquals("", t);
        assertEquals(view(f[0]), exported(tio, "c.bam"));
    }

    @Test
    void tagsFirstAppearingInALaterBlock() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("d"), 600, true, 250, 11);
        Path tio = importSam(f[0], f[1], 100, "d.tio");
        List<Map<String, Object>> rows = indexRows(tio);
        assertEquals(0L, ((Number) rows.get(0).get("tags_len")).longValue());
        assertEquals(0L, ((Number) rows.get(1).get("tags_len")).longValue());
        assertTrue(((Number) rows.get(rows.size() - 1).get("tags_len")).longValue() > 0);
        assertEquals(viewTags(f[0]), storedTags(tio));
        assertTrue(hasFlag(tio, FeatureFlags.OPT_SAM_TAGS));
    }

    @Test
    void wholeChannelLayout() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("e"), 120, true, 0, 11);
        WrittenGenomicRun run = new BamReader(f[0]).toGenomicRun("g").withOptLegacyWholeChannel(true);
        Path tio = tmp.resolve("e.tio");
        SpectralDataset.create(tio.toString(), "t", "i", List.of(), List.of(run),
            List.of(), List.of(), List.of(), FeatureFlags.defaultCurrent()).close();
        assertTrue(hasFlag(tio, FeatureFlags.OPT_SAM_TAGS));
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            StorageGroup sc = sp.rootGroup().openGroup("study").openGroup("genomic_runs")
                .openGroup(RUN).openGroup("signal_channels");
            try (StorageDataset ds = sc.openDataset("tags")) {
                assertEquals(Compression.SAM_TAGS.ordinal(),
                    ((Number) ds.getAttribute("compression")).intValue());
            }
        }
        assertEquals(viewTags(f[0]), storedTags(tio));
        assertEquals(view(f[0]), exported(tio, "e.bam"));
        // The text test seam appends the tags after column 11.
        String text = buildSamText(run);
        assertTrue(text.contains("\tAS:i:"), text);
    }

    private static String buildSamText(WrittenGenomicRun run) throws Exception {
        var m = BamWriter.class.getDeclaredMethod("buildSamText", WrittenGenomicRun.class,
            List.class, boolean.class);
        m.setAccessible(true);
        return (String) m.invoke(new BamWriter(Path.of("x.bam")), run, List.of(), false);
    }

    // ── protection ───────────────────────────────────────────────────

    private static byte[] key() {
        byte[] k = new byte[32];
        java.util.Arrays.fill(k, (byte) 0x65);
        return k;
    }

    @Test
    void perAuEncryptsTagsAndRestoresByteIdentical() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("p"), 600, true, 0, 11);
        Path tio = importSam(f[0], f[1], 200, "p.tio");
        List<Object> before = channelState(tio);
        List<String> want = storedTags(tio);

        PerAUFile.encryptFile(tio.toString(), key(), false, "hdf5");
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            StorageGroup sc = sp.rootGroup().openGroup("study").openGroup("genomic_runs")
                .openGroup(RUN).openGroup("signal_channels");
            assertFalse(sc.hasChild("tags"));
            assertTrue(sc.hasChild("tags_segments"));
        }
        assertEquals(want, PerAUFile.decryptFile(tio.toString(), key(), "hdf5").get(RUN).tags());

        PerAUFile.decryptFileInPlace(tio.toString(), key(), "hdf5");
        assertStateEquals(before, channelState(tio), "restore");
        assertEquals(view(f[0]), exported(tio, "p.bam"));
    }

    /** HDF5 filter count of a dataset (0 = unfiltered). */
    private static int filterCount(Path tio, String dataset) throws Exception {
        long f = hdf.hdf5lib.H5.H5Fopen(tio.toString(), hdf.hdf5lib.HDF5Constants.H5F_ACC_RDONLY,
            hdf.hdf5lib.HDF5Constants.H5P_DEFAULT);
        try {
            long d = hdf.hdf5lib.H5.H5Dopen(f, dataset, hdf.hdf5lib.HDF5Constants.H5P_DEFAULT);
            try {
                long p = hdf.hdf5lib.H5.H5Dget_create_plist(d);
                try {
                    return hdf.hdf5lib.H5.H5Pget_nfilters(p);
                } finally {
                    hdf.hdf5lib.H5.H5Pclose(p);
                }
            } finally {
                hdf.hdf5lib.H5.H5Dclose(d);
            }
        } finally {
            hdf.hdf5lib.H5.H5Fclose(f);
        }
    }

    /** Tags first appear after untagged blocks: the restore creates the
     *  tags dataset at the first non-empty blob (codec 18, unfiltered),
     *  as the stream writer did, so the file is byte-identical. */
    @Test
    void perAuRestoreOfLateTagsIsByteIdentical() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("pl"), 600, true, 250, 11);
        Path tio = importSam(f[0], f[1], 100, "pl.tio");
        String tagsPath = "/study/genomic_runs/" + RUN + "/signal_channels/tags";
        assertEquals(0, filterCount(tio, tagsPath), "stream writer: tags unfiltered");
        List<Object> before = channelState(tio);

        PerAUFile.encryptFile(tio.toString(), key(), false, "hdf5");
        PerAUFile.decryptFileInPlace(tio.toString(), key(), "hdf5");

        assertStateEquals(before, channelState(tio), "late-tags restore");
        assertEquals(0, filterCount(tio, tagsPath), "restored tags unfiltered");
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ, "hdf5")) {
            try (StorageDataset ds = sp.rootGroup().openGroup("study").openGroup("genomic_runs")
                    .openGroup(RUN).openGroup("signal_channels").openDataset("tags")) {
                assertEquals(Compression.SAM_TAGS.ordinal(),
                    ((Number) ds.getAttribute("compression")).intValue());
            }
        }
        assertEquals(viewTags(f[0]), storedTags(tio));
    }

    @Test
    void perAuRefusesWholeChannelAndRegionEncryptionOfTags() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("w"), 60, true, 0, 11);
        WrittenGenomicRun run = new BamReader(f[0]).toGenomicRun("g").withOptLegacyWholeChannel(true);
        Path tio = tmp.resolve("w.tio");
        SpectralDataset.create(tio.toString(), "t", "i", List.of(), List.of(run),
            List.of(), List.of(), List.of(), FeatureFlags.defaultCurrent()).close();
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
            () -> PerAUFile.encryptFile(tio.toString(), key(), false, "hdf5"));
        assertTrue(e.getMessage().contains("SAM tags"), e.getMessage());

        Path blocks = importSam(f[0], f[1], 20, "w2.tio");
        IllegalArgumentException r = assertThrows(IllegalArgumentException.class,
            () -> PerAUFile.encryptByRegion(blocks.toString(), Map.of("chr1", key()), "hdf5"));
        assertTrue(r.getMessage().contains("SAM tags"), r.getMessage());
    }

    @Test
    void signaturesCoverTheTags() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("s"), 120, true, 0, 11);
        Path tio = importSam(f[0], f[1], null, "s.tio");
        byte[] k = "kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk".getBytes(StandardCharsets.US_ASCII);
        try (StorageProvider sp = ProviderRegistry.open(tio.toString(), StorageProvider.Mode.READ_WRITE, "hdf5")) {
            StorageGroup rg = sp.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup(RUN);
            assertTrue(SignatureManager.signGenomicRun(rg, k).containsKey("signal_channels/tags"));
            assertTrue(SignatureManager.verifyGenomicRun(rg, k));
            try (StorageDataset ds = rg.openGroup("signal_channels").openDataset("tags")) {
                byte[] data = (byte[]) ds.readAll();
                data[data.length - 1] ^= 1;
                ds.writeAll(data);
            }
            assertFalse(SignatureManager.verifyGenomicRun(rg, k));
        }
    }

    // ── transport ────────────────────────────────────────────────────

    private static byte[] plaintextStream(Path tio) throws Exception {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (SpectralDataset ds = SpectralDataset.open(tio.toString());
             TransportWriter tw = new TransportWriter(bos)) {
            tw.writeDataset(ds);
        }
        return bos.toByteArray();
    }

    @Test
    void plaintextTransportCarriesTags() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("t"), 600, true, 0, 11);
        Path src = importSam(f[0], f[1], 200, "t.tio");
        byte[] stream = plaintextStream(src);
        Path out = tmp.resolve("t-out.tio");
        try (TransportReader tr = new TransportReader(stream)) {
            tr.materializeTo(out.toString()).close();
        }
        assertEquals(viewTags(f[0]), storedTags(out));
        assertTrue(hasFlag(out, FeatureFlags.OPT_SAM_TAGS));
    }

    @Test
    void taglessPlaintextStreamIsUnchanged() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("u"), 60, false, 0, 11);
        Path src = importSam(f[0], null, null, "u.tio");
        byte[] stream = plaintextStream(src);
        int aus = 0;
        try (TransportReader tr = new TransportReader(stream)) {
            for (TransportReader.PacketRecord p : tr.readAllPackets()) {
                if (p.header.packetType != PacketType.ACCESS_UNIT) continue;
                aus++;
                for (ChannelData c : AccessUnit.decode(p.payload).channels) {
                    assertNotEquals("tags", c.name);
                }
            }
        }
        assertEquals(61, aus);
    }

    @Test
    void encryptedTransportRestoresTags() throws Exception {
        needs();
        Path[] f = makeSam(tmp.resolve("x"), 600, true, 0, 11);
        Path tio = importSam(f[0], f[1], 200, "x.tio");
        List<Object> before = channelState(tio);
        Map<String, byte[]> refSeqs = readFasta(f[1]);
        String refUri;
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            refUri = ds.genomicRuns().get(RUN).referenceUri();
        }
        PerAUFile.encryptFile(tio.toString(), key(), false, "hdf5");
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (TransportWriter w = new TransportWriter(bos)) {
            EncryptedTransport.writeEncryptedDataset(tio.toString(), w, "hdf5");
        }
        Path rx = tmp.resolve("x-rx.tio");
        EncryptedTransport.readEncryptedToPath(rx.toString(), bos.toByteArray(), "hdf5");
        assertTrue(indexRows(rx).get(0).containsKey("tags_off"), "receiver index carries the triple");
        // The stream carries no reference bytes (the M99.1 limit): embed
        // the reference in the received container instead of REF_PATH.
        try (StorageProvider sp = ProviderRegistry.open(rx.toString(), StorageProvider.Mode.READ_WRITE, "hdf5")) {
            StorageGroup study = sp.rootGroup().openGroup("study");
            WrittenGenomicRun carrier = new WrittenGenomicRun(
                Enums.AcquisitionMode.GENOMIC_WGS, refUri, "", "",
                new long[0], new byte[0], new int[0], new byte[0], new byte[0], new long[0], new int[0],
                List.of(), List.of(), List.of(), new long[0], new int[0], List.of(),
                Compression.ZLIB, Map.of(), List.of(), true, refSeqs, null);
            SpectralDatasetGenomicWriter.embedReferencesForRuns(study, List.of(carrier));
        }
        PerAUFile.decryptFileInPlace(rx.toString(), key(), "hdf5");
        List<Object> after = channelState(rx);
        assertEquals(before.get(0), after.get(0), "blocks/index");
        assertArrayEquals((byte[]) ((List<?>) before.get(1)).get(0),
            (byte[]) ((List<?>) after.get(1)).get(0), "tags blob");
        assertEquals(viewTags(f[0]), storedTags(rx));
    }

    private static Map<String, byte[]> readFasta(Path fa) throws IOException {
        Map<String, byte[]> out = new LinkedHashMap<>();
        String name = null;
        StringBuilder sb = new StringBuilder();
        for (String line : Files.readAllLines(fa)) {
            if (line.startsWith(">")) {
                if (name != null) out.put(name, sb.toString().getBytes(StandardCharsets.US_ASCII));
                name = line.substring(1).trim();
                sb.setLength(0);
            } else {
                sb.append(line.trim());
            }
        }
        if (name != null) out.put(name, sb.toString().getBytes(StandardCharsets.US_ASCII));
        return out;
    }

    @Test
    void alignedReadKeepsTheOldConstructor() {
        AlignedRead r = new AlignedRead("q", "chr1", 1, 60, "4M", "ACGT", new byte[4], 0, "", -1, 0);
        assertEquals("", r.tags());
    }
}
