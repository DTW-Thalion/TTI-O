/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.codecs.SeqGroup;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicBlocks;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.GenomicStreamWriter;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.protection.EncryptedTransport;
import global.thalion.ttio.protection.PerAUFile;
import global.thalion.ttio.protection.SignatureManager;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import global.thalion.ttio.transport.PacketType;
import global.thalion.ttio.transport.TransportReader;
import global.thalion.ttio.transport.TransportWriter;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Consumer;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.junit.jupiter.api.Assertions.*;

/**
 * M103: reads grouped by sequence before blocking (layout
 * {@code blocks_v1_grouped}).
 *
 * <p>The writer reorders an unaligned run's reads with the native kernel,
 * cuts blocks over the reordered run and stores
 * {@code genomic_index/input_index}; every read accessor presents input
 * order, block iteration stays in stored order, and the grouping survives
 * per-AU encryption, the encrypted transport (BlockSidecar input_index,
 * transport-spec v0.13) and signatures. Plaintext transport delivers input
 * order and the receiver writes it ungrouped.</p>
 *
 * <p>Mirrors {@code python/tests/test_m103_read_grouping.py}.</p>
 */
final class M103ReadGroupingTest {

    @TempDir Path tmp;

    private static byte[] key() {
        byte[] k = new byte[32];
        for (int i = 0; i < 32; i++) k[i] = (byte) i;
        return k;
    }

    @BeforeEach
    void needsNative() {
        Assumptions.assumeTrue(SeqGroup.isAvailable(), "libttio_rans_jni not loaded");
    }

    private static byte[] revComp(byte[] s) {
        byte[] out = new byte[s.length];
        for (int i = 0; i < s.length; i++) {
            byte b = s[s.length - 1 - i];
            out[i] = switch (b) {
                case 'A' -> (byte) 'T';
                case 'C' -> (byte) 'G';
                case 'G' -> (byte) 'C';
                case 'T' -> (byte) 'A';
                default -> b;
            };
        }
        return out;
    }

    /** Paired reads from both strands of a random genome, shuffled, with
     *  a few N bases and a zero-length read: about 12x coverage. */
    private static WrittenGenomicRun run(int seed, int n) {
        int len = 120, g = 30_000;
        Random rng = new Random(seed);
        byte[] acgt = {'A', 'C', 'G', 'T'};
        byte[] genome = new byte[g];
        for (int i = 0; i < g; i++) genome[i] = acgt[rng.nextInt(4)];
        List<byte[]> seqs = new ArrayList<>();
        List<String> names = new ArrayList<>();
        for (int i = 0; i < n / 2; i++) {
            int p = rng.nextInt(g - 500);
            byte[] r1 = Arrays.copyOfRange(genome, p, p + len);
            byte[] r2 = revComp(Arrays.copyOfRange(genome, p + 300, p + 300 + len));
            if (i % 37 == 0) r1[50] = 'N';
            seqs.add(r1);
            seqs.add(r2);
            names.add("frag" + i + "/1");
            names.add("frag" + i + "/2");
        }
        seqs.set(5, new byte[0]);
        int m = seqs.size();
        Integer[] perm = new Integer[m];
        for (int i = 0; i < m; i++) perm[i] = i;
        java.util.Collections.shuffle(Arrays.asList(perm), rng);
        List<byte[]> sh = new ArrayList<>(m);
        List<String> shNames = new ArrayList<>(m);
        for (int i : perm) { sh.add(seqs.get(i)); shNames.add(names.get(i)); }
        int[] lengths = new int[m];
        int total = 0;
        for (int i = 0; i < m; i++) { lengths[i] = sh.get(i).length; total += lengths[i]; }
        byte[] seq = new byte[total];
        byte[] qual = new byte[total];
        long[] offsets = new long[m];
        int o = 0;
        for (int i = 0; i < m; i++) {
            offsets[i] = o;
            System.arraycopy(sh.get(i), 0, seq, o, lengths[i]);
            o += lengths[i];
        }
        for (int i = 0; i < total; i++) qual[i] = (byte) (33 + rng.nextInt(40));
        long[] positions = new long[m];
        byte[] mapqs = new byte[m];
        int[] flags = new int[m];
        long[] matePos = new long[m];
        for (int i = 0; i < m; i++) {
            positions[i] = i;                     // distinct per read, to check order
            mapqs[i] = (byte) (i % 61);
            flags[i] = 0x4;
            matePos[i] = -1;
        }
        List<String> stars = java.util.Collections.nCopies(m, "*");
        return new WrittenGenomicRun(
            AcquisitionMode.GENOMIC_WGS, "", "ILLUMINA", "M103G",
            positions, mapqs, flags, seq, qual, offsets, lengths,
            stars, shNames, stars, matePos, new int[m], stars,
            Compression.ZLIB, Map.of(), List.of(), false, null, null,
            null, false, false, null, 0L, null);
    }

    private Path write(String name, WrittenGenomicRun run, boolean group, int blockReads) {
        Path path = tmp.resolve(name);
        SpectralDataset.create(path.toString(), "m103", "i",
            List.of(), List.of(), List.of(), List.of(), List.of(),
            FeatureFlags.defaultCurrent()).close();
        try (StorageProvider sp = ProviderRegistry.open(path.toString(),
                StorageProvider.Mode.READ_WRITE, "hdf5");
             StorageGroup study = sp.rootGroup().openGroup("study");
             GenomicStreamWriter w = new GenomicStreamWriter(study, "g",
                 GenomicStreamWriter.Options.fromRun(run)
                     .withBlockPolicy(blockReads, Long.MAX_VALUE).withGroupReads(group))) {
            int n = run.readCount();
            for (int s = 0; s < n; s += 500) {
                w.appendBatch(GenomicBlocks.sliceRun(run, s, Math.min(s + 500, n)));
            }
        }
        return path;
    }

    private static List<String> expectedSeqs(WrittenGenomicRun run) {
        List<String> out = new ArrayList<>();
        for (int i = 0; i < run.readCount(); i++) {
            out.add(new String(run.sequences(), (int) run.offsets()[i], run.lengths()[i],
                StandardCharsets.US_ASCII));
        }
        return out;
    }

    private static void withRun(Path path, StorageProvider.Mode mode, Consumer<StorageGroup> fn) {
        try (StorageProvider sp = ProviderRegistry.open(path.toString(), mode, "hdf5");
             StorageGroup rg = sp.rootGroup().openGroup("study").openGroup("genomic_runs")
                 .openGroup("g")) {
            fn.accept(rg);
        }
    }

    private static int[] ints(Object raw) {
        if (raw instanceof int[] a) return a;
        long[] l = (long[]) raw;
        int[] out = new int[l.length];
        for (int i = 0; i < l.length; i++) out[i] = (int) l[i];
        return out;
    }

    private static int[] readIndexColumn(StorageGroup rg, String name) {
        try (StorageGroup idx = rg.openGroup("genomic_index");
             StorageDataset ds = idx.openDataset(name)) {
            return ints(ds.readAll());
        }
    }

    private static long[] readPositions(StorageGroup rg) {
        try (StorageGroup idx = rg.openGroup("genomic_index");
             StorageDataset ds = idx.openDataset("positions")) {
            return (long[]) ds.readAll();
        }
    }

    private static String layoutOf(StorageGroup rg) {
        Object v = rg.getAttribute("layout");
        return v instanceof byte[] b ? new String(b, StandardCharsets.UTF_8) : v.toString();
    }

    private static List<String> namesInOrder(Path path) {
        try (SpectralDataset ds = SpectralDataset.open(path.toString())) {
            GenomicRun g = ds.genomicRuns().values().iterator().next();
            List<String> out = new ArrayList<>();
            Iterator<AlignedRead> it = g.iterReads();
            while (it.hasNext()) out.add(it.next().readName());
            return out;
        }
    }

    @Test
    void groupedLayoutAndInputOrder() {
        WrittenGenomicRun run = run(7, 3000);
        int n = run.readCount();
        Path path = write("grp.tio", run, true, 700);
        withRun(path, StorageProvider.Mode.READ, rg -> {
            assertEquals("blocks_v1_grouped", layoutOf(rg));
            int[] ii = readIndexColumn(rg, "input_index");
            int[] sorted = ii.clone();
            Arrays.sort(sorted);
            for (int i = 0; i < n; i++) assertEquals(i, sorted[i]);
            boolean identity = true;
            for (int i = 0; i < n; i++) identity &= ii[i] == i;
            assertFalse(identity, "the reads were reordered");
            long[] storedPos = readPositions(rg);
            for (int j = 0; j < n; j++) {
                assertEquals(ii[j], storedPos[j], "every index column moved with its read");
            }
            try (StorageGroup blocks = rg.openGroup("blocks");
                 StorageDataset ds = blocks.openDataset("index")) {
                assertTrue(ds.readRows().size() > 1);
            }
        });
        List<String> seqs = expectedSeqs(run);
        try (SpectralDataset ds = SpectralDataset.open(path.toString())) {
            GenomicRun g = ds.genomicRuns().get("g");
            assertEquals("blocks_v1_grouped", g.layout());
            assertEquals(n, g.readCount());
            Iterator<AlignedRead> it = g.iterReads();
            for (int i = 0; i < n; i++) {
                AlignedRead r = it.next();
                assertEquals(run.readNames().get(i), r.readName());
                assertEquals(seqs.get(i), r.sequence());
                assertEquals(i, r.position());
            }
            assertFalse(it.hasNext());
            for (int i : new int[]{0, 1, 5, n / 2, n - 1}) {
                AlignedRead r = g.readAt(i);
                assertEquals(run.readNames().get(i), r.readName());
                assertEquals(seqs.get(i), r.sequence());
                assertEquals(run.readNames().get(i), g.readNameAt(i));
            }
            for (int i = 0; i < n; i++) {
                assertEquals(i, g.index().positionAt(i), "index in input order");
                assertEquals(run.lengths()[i], g.index().lengthAt(i));
            }
            assertEquals(run.readNames(), g.readNamesAll());
            assertArrayEquals(run.sequences(), g.sequencesFull());
            assertArrayEquals(run.qualities(), g.qualitiesFull());
            List<String> part = new ArrayList<>();
            Iterator<AlignedRead> pit = g.iterReads(1000, 1100);
            while (pit.hasNext()) part.add(pit.next().readName());
            assertEquals(run.readNames().subList(1000, 1100), part);
            // The threaded iterator also presents input order.
            Iterator<AlignedRead> tit = g.iterReads(0, n, 4);
            for (int i = 0; i < n; i++) assertEquals(run.readNames().get(i), tit.next().readName());
            // Block iteration is in stored order; inputIndex maps its rows.
            int[] ii = g.inputIndex();
            Map<Integer, String> seen = new HashMap<>();
            g.iterBlocks(0, n, 1, (view, viewStart, firstRead, nReads) -> {
                for (int k = 0; k < nReads; k++) {
                    seen.put(ii[firstRead + k], view.readAt(viewStart + k).readName());
                }
            });
            for (int i = 0; i < n; i++) assertEquals(run.readNames().get(i), seen.get(i));
        }
    }

    @Test
    void groupingShrinksSequences() {
        WrittenGenomicRun run = run(7, 6000);
        Path a = write("plain.tio", run, false, 1000);
        Path b = write("grp.tio", run, true, 1000);
        long[] sizes = new long[2];
        int k = 0;
        for (Path p : new Path[]{a, b}) {
            final int slot = k++;
            withRun(p, StorageProvider.Mode.READ, rg -> {
                try (StorageGroup sc = rg.openGroup("signal_channels");
                     StorageDataset ds = sc.openGroup("sequences").openDataset("data")) {
                    sizes[slot] = ds.shape()[0];
                }
            });
        }
        assertTrue(sizes[1] < sizes[0] * 0.8, sizes[0] + " -> " + sizes[1]);
    }

    @Test
    void defaultWriteWithOptGroupReads() {
        WrittenGenomicRun run = run(7, 800).withOptGroupReads(true);
        assertTrue(run.optGroupReads());
        Path path = tmp.resolve("wm.tio");
        SpectralDataset.create(path.toString(), "m103", "i", List.of(), List.of(run),
            List.of(), List.of(), List.of(), FeatureFlags.defaultCurrent()).close();
        try (SpectralDataset ds = SpectralDataset.open(path.toString())) {
            GenomicRun g = ds.genomicRuns().values().iterator().next();
            assertEquals("blocks_v1_grouped", g.layout());
        }
        assertEquals(run.readNames(), namesInOrder(path));
    }

    @Test
    void groupingRefusesAlignedAndLegacy() {
        Path path = tmp.resolve("x.tio");
        SpectralDataset.create(path.toString(), "m103", "i",
            List.of(), List.of(), List.of(), List.of(), List.of(),
            FeatureFlags.defaultCurrent()).close();
        WrittenGenomicRun run = run(7, 20);
        try (StorageProvider sp = ProviderRegistry.open(path.toString(),
                StorageProvider.Mode.READ_WRITE, "hdf5");
             StorageGroup study = sp.rootGroup().openGroup("study")) {
            var o = GenomicStreamWriter.Options.fromRun(run).withGroupReads(true);
            assertThrows(IllegalArgumentException.class, () -> new GenomicStreamWriter(study, "a",
                o.withReference(Map.of("chr1", "ACGT".getBytes(StandardCharsets.US_ASCII)), false)));
            assertThrows(IllegalArgumentException.class,
                () -> new GenomicStreamWriter(study, "b", o.withLegacy(true)));
        }
    }

    @Test
    void perAuRestoreKeepsGrouping() throws Exception {
        WrittenGenomicRun run = run(7, 1500);
        Path path = write("grp.tio", run, true, 400);
        Path pristine = tmp.resolve("pristine.tio");
        Files.copy(path, pristine);
        PerAUFile.encryptFile(path.toString(), key(), false, "hdf5");
        PerAUFile.decryptFileInPlace(path.toString(), key(), "hdf5");
        Object[][] got = snapshot(path), want = snapshot(pristine);
        assertEquals(want[0][0], got[0][0], "blocks/index");
        for (int k = 1; k < 4; k++) {
            assertTrue(java.util.Objects.deepEquals(want[0][k], got[0][k]), "column " + k);
        }
        assertEquals(run.readNames(), namesInOrder(path));
    }

    /** blocks/index rows, input_index, sequences/data, qualities. */
    private static Object[][] snapshot(Path p) {
        Object[][] out = new Object[1][4];
        withRun(p, StorageProvider.Mode.READ, rg -> {
            try (StorageGroup blocks = rg.openGroup("blocks");
                 StorageDataset ds = blocks.openDataset("index")) {
                out[0][0] = ds.readRows();
            }
            out[0][1] = readIndexColumn(rg, "input_index");
            try (StorageGroup sc = rg.openGroup("signal_channels")) {
                try (StorageDataset ds = sc.openGroup("sequences").openDataset("data")) {
                    out[0][2] = ds.readAll();
                }
                try (StorageDataset ds = sc.openDataset("qualities")) {
                    out[0][3] = ds.readAll();
                }
            }
        });
        return out;
    }

    @Test
    void encryptedTransportCarriesInputIndex() throws Exception {
        WrittenGenomicRun run = run(7, 1500);
        Path path = write("grp.tio", run, true, 400);
        int[][] wantIi = new int[1][];
        withRun(path, StorageProvider.Mode.READ, rg -> wantIi[0] = readIndexColumn(rg, "input_index"));
        PerAUFile.encryptFile(path.toString(), key(), false, "hdf5");
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (TransportWriter writer = new TransportWriter(bos)) {
            EncryptedTransport.writeEncryptedDataset(path.toString(), writer, "hdf5");
        }
        byte[] stream = bos.toByteArray();
        boolean announced = false;
        try (TransportReader reader = new TransportReader(stream)) {
            for (TransportReader.PacketRecord r : reader.readAllPackets()) {
                if (r.header.packetType == PacketType.STREAM_HEADER) {
                    announced = new String(r.payload, StandardCharsets.UTF_8)
                        .contains(PacketType.TRANSPORT_BLOCKS_V1_GROUPED_FEATURE);
                }
            }
        }
        assertTrue(announced, "StreamHeader carries transport_blocks_v1_grouped");
        Path got = tmp.resolve("got.tio");
        EncryptedTransport.readEncryptedToPath(got.toString(), stream, "hdf5");
        withRun(got, StorageProvider.Mode.READ, rg -> {
            assertEquals("blocks_v1_grouped", layoutOf(rg));
            assertArrayEquals(wantIi[0], readIndexColumn(rg, "input_index"));
        });
        try (StorageProvider sp = ProviderRegistry.open(got.toString(),
                StorageProvider.Mode.READ, "hdf5")) {
            java.util.Set<String> feats = FeatureFlags.readFrom(sp.rootGroup()).features();
            assertFalse(feats.contains(PacketType.TRANSPORT_BLOCKS_V1_GROUPED_FEATURE));
            assertFalse(feats.contains(PacketType.TRANSPORT_BLOCKS_V1_FEATURE));
        }
        PerAUFile.decryptFileInPlace(got.toString(), key(), "hdf5");
        assertEquals(run.readNames(), namesInOrder(got));
    }

    @Test
    void plaintextTransportDeliversInputOrder() throws Exception {
        WrittenGenomicRun run = run(7, 1200);
        Path path = write("grp.tio", run, true, 300);
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (SpectralDataset ds = SpectralDataset.open(path.toString());
             TransportWriter tw = new TransportWriter(bos)) {
            tw.writeDataset(ds);
        }
        Path out = tmp.resolve("out.tio");
        try (TransportReader tr = new TransportReader(bos.toByteArray())) {
            tr.materializeTo(out.toString()).close();
        }
        try (SpectralDataset ds = SpectralDataset.open(out.toString())) {
            GenomicRun g = ds.genomicRuns().values().iterator().next();
            assertNotEquals("blocks_v1_grouped", g.layout());
            List<String> seqs = expectedSeqs(run);
            Iterator<AlignedRead> it = g.iterReads();
            for (int i = 0; i < run.readCount(); i++) {
                AlignedRead r = it.next();
                assertEquals(run.readNames().get(i), r.readName());
                assertEquals(seqs.get(i), r.sequence());
            }
        }
    }

    @Test
    void signaturesCoverInputIndex() {
        Path path = write("grp.tio", run(7, 600), true, 200);
        byte[] k = new byte[32];
        Arrays.fill(k, (byte) 'k');
        withRun(path, StorageProvider.Mode.READ_WRITE, rg -> {
            assertTrue(SignatureManager.signGenomicRun(rg, k).containsKey("genomic_index/input_index"));
            assertTrue(SignatureManager.verifyGenomicRun(rg, k));
            swapFirstTwo(rg);
            assertFalse(SignatureManager.verifyGenomicRun(rg, k));
        });
    }

    private static void swapFirstTwo(StorageGroup rg) {
        int[] ii = readIndexColumn(rg, "input_index");
        try (StorageGroup idx = rg.openGroup("genomic_index");
             StorageDataset ds = idx.openDataset("input_index")) {
            ds.writeSlice(0, new int[]{ii[1], ii[0]});
        }
    }

    @Test
    void corruptInputIndexIsRefused() {
        Path path = write("grp.tio", run(7, 600), true, 200);
        withRun(path, StorageProvider.Mode.READ_WRITE, rg -> {
            int[] ii = readIndexColumn(rg, "input_index");
            try (StorageGroup idx = rg.openGroup("genomic_index");
                 StorageDataset ds = idx.openDataset("input_index")) {
                ds.writeSlice(0, new int[]{ii[1]});
            }
        });
        IllegalStateException e = assertThrows(IllegalStateException.class, () -> {
            try (SpectralDataset ds = SpectralDataset.open(path.toString())) {
                ds.genomicRuns().get("g");
            }
        });
        assertTrue(String.valueOf(e.getMessage()).contains("permutation"), e.getMessage());
    }

    @Test
    void groupOrderKeepsLabelsApart() {
        // Two labels interleaved: the order groups within each label,
        // labels in first-seen order (a block never spans two).
        WrittenGenomicRun base = run(11, 400);
        int n = base.readCount();
        List<String> chroms = new ArrayList<>(n);
        for (int i = 0; i < n; i++) chroms.add(i % 3 == 0 ? "b" : "a");
        WrittenGenomicRun r = new WrittenGenomicRun(
            base.acquisitionMode(), base.referenceUri(), base.platform(), base.sampleName(),
            base.positions(), base.mappingQualities(), base.flags(), base.sequences(),
            base.qualities(), base.offsets(), base.lengths(), base.cigars(), base.readNames(),
            base.mateChromosomes(), base.matePositions(), base.templateLengths(), chroms,
            Compression.ZLIB, Map.of(), List.of(), false, null, null,
            null, false, false, null, 0L, null);
        int[] order = GenomicBlocks.groupOrder(r);
        assertEquals(n, order.length);
        // "b" is seen first (read 0), so its reads come first.
        int nb = 0;
        for (String c : chroms) if (c.equals("b")) nb++;
        for (int j = 0; j < n; j++) {
            assertEquals(j < nb ? "b" : "a", chroms.get(order[j]), "row " + j);
        }
    }
}
