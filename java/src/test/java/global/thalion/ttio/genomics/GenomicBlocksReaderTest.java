/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.importers.BamReader;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/** Reading {@code blocks_v1} runs through {@link GenomicRun}. */
class GenomicBlocksReaderTest {

    static WrittenGenomicRun m87() throws IOException {
        return new BamReader(Paths.get("src", "test", "resources", "ttio",
                "fixtures", "genomic", "m87_test.bam")).toGenomicRun("genomic_0001");
    }

    static Path writeBlocks(Path tmp, int blockReads) throws Exception {
        WrittenGenomicRun run = m87();
        Path out = tmp.resolve("b" + blockReads + ".tio");
        try (StorageProvider p = ProviderRegistry.open(out.toString(), StorageProvider.Mode.CREATE, "hdf5")) {
            StorageGroup study = p.rootGroup().createGroup("study");
            try (GenomicStreamWriter w = new GenomicStreamWriter(study, "genomic_0001",
                    GenomicStreamWriter.Options.fromRun(run).withBlockPolicy(blockReads, Long.MAX_VALUE))) {
                w.appendBatch(run);
            }
        }
        return out;
    }

    @ParameterizedTest
    @ValueSource(ints = {1, 3, 1_000_000})
    void readsAgreeWithWholeRunDecode(int blockReads, @TempDir Path tmp) throws Exception {
        WrittenGenomicRun run = m87();
        try (SpectralDataset ds = SpectralDataset.open(writeBlocks(tmp, blockReads).toString())) {
            GenomicRun g = ds.genomicRuns().get("genomic_0001");
            assertEquals("blocks_v1", g.layout());
            assertEquals(run.readCount(), g.readCount());
            if (blockReads >= run.readCount()) {
                // one block per chromosome run in a coordinate-sorted BAM
                int runs = 1;
                for (int i = 1; i < run.readCount(); i++) {
                    if (!run.chromosomes().get(i).equals(run.chromosomes().get(i - 1))) runs++;
                }
                assertEquals(runs, g.blockCount());
            } else {
                assertTrue(g.blockCount() >= (run.readCount() + blockReads - 1) / blockReads);
            }
            for (int i = 0; i < run.readCount(); i++) {
                AlignedRead r = g.readAt(i);
                assertEquals(run.readNames().get(i), r.readName(), "read " + i);
                assertEquals(run.cigars().get(i), r.cigar(), "read " + i);
                assertEquals(run.chromosomes().get(i), r.chromosome(), "read " + i);
                assertEquals(run.positions()[i], r.position());
                assertEquals(run.mappingQualities()[i] & 0xFF, r.mappingQuality());
                assertEquals(run.flags()[i], r.flags());
                String mate = run.mateChromosomes().get(i);
                assertEquals(mate.equals("=") ? run.chromosomes().get(i) : mate, r.mateChromosome(), "read " + i);
                assertEquals(run.matePositions()[i], r.matePosition());
                assertEquals(run.templateLengths()[i], r.templateLength());
                int o = (int) run.offsets()[i], l = run.lengths()[i];
                assertEquals(new String(run.sequences(), o, l, StandardCharsets.US_ASCII), r.sequence());
                assertArrayEquals(Arrays.copyOfRange(run.qualities(), o, o + l), r.qualities());
            }
            int n = 0;
            Iterator<AlignedRead> it = g.iterReads(0, g.readCount());
            while (it.hasNext()) { it.next(); n++; }
            assertEquals(run.readCount(), n);
            List<AlignedRead> region = g.readsInRegion(run.chromosomes().get(0), 0, Long.MAX_VALUE);
            assertFalse(region.isEmpty());
            assertArrayEquals(run.sequences(), g.sequencesFull());
            assertArrayEquals(run.qualities(), g.qualitiesFull());
            assertEquals(run.readNames(), g.readNamesAll());
            assertEquals(new java.util.LinkedHashSet<>(run.chromosomes()).size(), g.chromosomeNames().size());
        }
    }

    @Test
    void unknownLayoutIsRejected(@TempDir Path tmp) throws Exception {
        Path f = writeBlocks(tmp, 3);
        try (StorageProvider p = ProviderRegistry.open(f.toString(), StorageProvider.Mode.READ_WRITE, "hdf5")) {
            p.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup("genomic_0001")
                .setAttribute("layout", "blocks_v9");
        }
        assertThrows(IllegalStateException.class, () -> SpectralDataset.open(f.toString()));
    }

    @Test
    void signaturesCoverTheBlockLayout(@TempDir Path tmp) throws Exception {
        Path f = writeBlocks(tmp, 3);
        byte[] key = "0123456789abcdef0123456789abcdef".getBytes(StandardCharsets.US_ASCII);
        try (StorageProvider p = ProviderRegistry.open(f.toString(), StorageProvider.Mode.READ_WRITE, "hdf5")) {
            StorageGroup rg = p.rootGroup().openGroup("study").openGroup("genomic_runs").openGroup("genomic_0001");
            var out = global.thalion.ttio.protection.SignatureManager.signGenomicRun(rg, key);
            assertTrue(out.containsKey("signal_channels/sequences/data"), out.keySet().toString());
            assertTrue(out.containsKey("signal_channels/qualities"));
            assertTrue(out.containsKey("blocks/index"));
            assertTrue(out.containsKey("genomic_index/lengths"));
            assertTrue(global.thalion.ttio.protection.SignatureManager.verifyGenomicRun(rg, key));
            try (var ds = rg.openGroup("signal_channels").openGroup("sequences").openDataset("data")) {
                byte[] first = (byte[]) ds.readSlice(0, 1);
                ds.writeSlice(0, new byte[]{ (byte) (first[0] ^ 0x55) });
            }
            assertFalse(global.thalion.ttio.protection.SignatureManager.verifyGenomicRun(rg, key));
        }
    }

    @Test
    void partialFileReadsUpToLastBlock(@TempDir Path tmp) throws Exception {
        WrittenGenomicRun run = m87();
        Path out = tmp.resolve("partial.tio");
        try (StorageProvider p = ProviderRegistry.open(out.toString(), StorageProvider.Mode.CREATE, "hdf5")) {
            StorageGroup study = p.rootGroup().createGroup("study");
            GenomicStreamWriter w = new GenomicStreamWriter(study, "genomic_0001",
                    GenomicStreamWriter.Options.fromRun(run).withBlockPolicy(2, Long.MAX_VALUE), 1);
            w.appendBatch(GenomicBlocks.sliceRun(run, 0, 4));   // two blocks flushed
            w.appendBatch(GenomicBlocks.sliceRun(run, 4, 5));   // pending; no close()
        }
        try (StorageProvider p = ProviderRegistry.open(out.toString(), StorageProvider.Mode.READ, "hdf5")) {
            BlockTable t = BlockTable.read(p.rootGroup().openGroup("study")
                    .openGroup("genomic_runs").openGroup("genomic_0001"));
            assertEquals(2, t.count());
            assertEquals(4, t.readCount());
        }
    }

    @org.junit.jupiter.api.Test
    void iterReadsThreadedMatchesSerial() throws Exception {
        // reference-less (SEQ_CM sequences, M103): a memory-provider readFrom has no resolver
        WrittenGenomicRun run = GenomicStreamWriterTest.bigSyntheticRun(30_000, 11).withReference(false, null, null);
        StorageGroup study = GenomicStreamWriterTest.writeWithThreads("memory://gbrt-th", run, 1, 5_000);
        try (GenomicRun g = GenomicRun.readFrom(study.openGroup("genomic_runs").openGroup("g"), "g")) {
            java.util.List<String> serial = new java.util.ArrayList<>();
            for (var it = g.iterReads(0, g.readCount(), 1); it.hasNext(); ) {
                var r = it.next();
                serial.add(r.readName() + "|" + r.sequence() + "|" + r.cigar() + "|" + r.mateChromosome());
            }
            java.util.List<String> threaded = new java.util.ArrayList<>();
            for (var it = g.iterReads(0, g.readCount(), 4); it.hasNext(); ) {
                var r = it.next();
                threaded.add(r.readName() + "|" + r.sequence() + "|" + r.cigar() + "|" + r.mateChromosome());
            }
            org.junit.jupiter.api.Assertions.assertEquals(serial, threaded);
            java.util.List<String> part = new java.util.ArrayList<>();
            for (var it = g.iterReads(12_345, 17_890, 3); it.hasNext(); ) part.add(it.next().readName());
            org.junit.jupiter.api.Assertions.assertEquals(serial.subList(12_345, 17_890).stream()
                .map(x -> x.split("\\|")[0]).toList(), part);
        }
    }
}
