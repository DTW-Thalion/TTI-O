/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.codecs.FqzcompNx16Z;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicBlocks;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.GenomicWriteContext;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.importers.BamReader;
import global.thalion.ttio.providers.MemoryProvider;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/** The whole-channel genomic writer over a memory-provider root, and the
 *  block encoder built on it. */
class GenomicBlocksTest {

    static WrittenGenomicRun m87() throws IOException {
        return new BamReader(Paths.get("src", "test", "resources", "ttio",
                "fixtures", "genomic", "m87_test.bam")).toGenomicRun("genomic_0001");
    }

    static StorageGroup memRoot(String url) {
        return new MemoryProvider().open(url, StorageProvider.Mode.CREATE).rootGroup();
    }

    @Test
    void v18WriterRoundTripsThroughMemoryProvider() throws Exception {
        WrittenGenomicRun run = m87();
        StorageGroup root = memRoot("memory://gb-roundtrip");
        SpectralDatasetGenomicWriter.writeGenomicRunSubtree(root, "r", run,
                GenomicWriteContext.none());
        try (GenomicRun g = GenomicRun.readFrom(root.openGroup("r"), "r")) {
            assertEquals(run.readCount(), g.readCount());
            for (int i = 0; i < run.readCount(); i++) {
                AlignedRead r = g.readAt(i);
                assertEquals(run.readNames().get(i), r.readName());
                assertEquals(run.cigars().get(i), r.cigar());
                assertEquals(run.chromosomes().get(i), r.chromosome());
                assertEquals(run.positions()[i], r.position());
            }
        }
        MemoryProvider.discardStore("memory://gb-roundtrip");
    }

    @Test
    void sharedChromMapGivesStableIdsAcrossWrites() throws Exception {
        WrittenGenomicRun run = m87();
        Map<String, Integer> shared = new LinkedHashMap<>();
        shared.put("chrZ", 0);
        StorageGroup root = memRoot("memory://gb-shared");
        SpectralDatasetGenomicWriter.writeGenomicRunSubtree(root, "a", run,
                new GenomicWriteContext(shared, null));
        short[] ids = (short[]) root.openGroup("a").openGroup("genomic_index")
                .openDataset("chromosome_ids").readAll();
        String first = run.chromosomes().get(0);
        assertTrue(shared.get(first) >= 1, "pre-seeded id 0 must survive");
        assertEquals(shared.get(first).intValue(), ids[0]);
        // the run-level name tables list names in id order, chrZ first
        var names = root.openGroup("a").openGroup("genomic_index")
                .openDataset("chromosome_names").readRows();
        assertEquals("chrZ", names.get(0).get("name").toString());
        MemoryProvider.discardStore("memory://gb-shared");
    }

    @Test
    void sliceAndConcatAreInverse() throws Exception {
        WrittenGenomicRun run = m87();
        WrittenGenomicRun a = GenomicBlocks.sliceRun(run, 0, 4);
        WrittenGenomicRun b = GenomicBlocks.sliceRun(run, 4, run.readCount());
        assertEquals(4, a.readCount());
        assertEquals(0L, a.offsets()[0]);
        assertEquals(0L, b.offsets()[0]);
        assertEquals(run.lengths()[4], b.lengths()[0]);
        WrittenGenomicRun back = GenomicBlocks.concatRuns(java.util.List.of(a, b));
        assertArrayEquals(run.sequences(), back.sequences());
        assertArrayEquals(run.qualities(), back.qualities());
        assertArrayEquals(run.offsets(), back.offsets());
        assertArrayEquals(run.lengths(), back.lengths());
        assertEquals(run.readNames(), back.readNames());
        assertEquals(run.mateChromosomes(), back.mateChromosomes());
        assertEquals(run.chromosomes(), back.chromosomes());
    }

    @Test
    void encodeBlockMatchesTheV18WriterBytes() throws Exception {
        WrittenGenomicRun run = m87();
        String chr = run.chromosomes().get(0);
        int stop = 0;
        while (stop < run.readCount() && run.chromosomes().get(stop).equals(chr)) stop++;
        WrittenGenomicRun block = GenomicBlocks.sliceRun(run, 0, stop);
        GenomicBlocks.BlockBlobs blobs = GenomicBlocks.encodeBlock(block,
                new GenomicWriteContext(new LinkedHashMap<>(), null));
        assertEquals(stop, blobs.nReads());
        assertEquals(Enums.Compression.RANS_ORDER0.ordinal(), blobs.codecs().get("cigars"));
        assertEquals(Enums.Compression.FQZCOMP_NX16_Z.ordinal(), blobs.codecs().get("qualities"));
        // No reference: SEQ_CM when the native library loaded (M103), else RANS_ORDER1.
        Enums.Compression expSeq = global.thalion.ttio.codecs.SeqCm.isAvailable()
            ? Enums.Compression.SEQ_CM : Enums.Compression.RANS_ORDER1;
        assertEquals(expSeq.ordinal(), blobs.codecs().get("sequences"));
        assertEquals(Enums.Compression.NAME_TOKENIZED_V2.ordinal(), blobs.codecs().get("read_names"));
        assertEquals(Enums.Compression.MATE_INLINE_V2.ordinal(), blobs.codecs().get("mate_info"));

        Map<String, Enums.Compression> ov = new LinkedHashMap<>();
        ov.put("cigars", Enums.Compression.RANS_ORDER0);
        ov.put("qualities", Enums.Compression.FQZCOMP_NX16_Z);
        ov.put("sequences", expSeq);
        WrittenGenomicRun same = block.withSignalCodecOverrides(ov);
        StorageGroup root = memRoot("memory://gb-cmp");
        SpectralDatasetGenomicWriter.writeGenomicRunSubtree(root, "r", same, GenomicWriteContext.none());
        StorageGroup sc = root.openGroup("r").openGroup("signal_channels");
        assertArrayEquals((byte[]) sc.openDataset("qualities").readAll(), blobs.blobs().get("qualities"));
        assertArrayEquals((byte[]) sc.openDataset("sequences").readAll(), blobs.blobs().get("sequences"));
        assertArrayEquals((byte[]) sc.openDataset("cigars").readAll(), blobs.blobs().get("cigars"));
        assertArrayEquals((byte[]) sc.openDataset("read_names").readAll(), blobs.blobs().get("read_names"));
        assertArrayEquals((byte[]) sc.openGroup("mate_info").openDataset("inline_v2").readAll(),
                blobs.blobs().get("mate_info"));
        MemoryProvider.discardStore("memory://gb-cmp");
    }

    @Test
    void encodeBlockForwardsQualStrategyHint() throws Exception {
        WrittenGenomicRun run = m87();
        String chr = run.chromosomes().get(0);
        int stop = 0;
        while (stop < run.readCount() && run.chromosomes().get(stop).equals(chr)) stop++;
        WrittenGenomicRun block = GenomicBlocks.sliceRun(run, 0, stop);
        GenomicBlocks.BlockBlobs def = GenomicBlocks.encodeBlock(block,
                new GenomicWriteContext(new LinkedHashMap<>(), null));
        // Small block: the size floor keeps auto on V4.
        assertEquals(4, FqzcompNx16Z.streamStrategy(def.blobs().get("qualities")));
        GenomicBlocks.BlockBlobs forced = GenomicBlocks.encodeBlock(block,
                new GenomicWriteContext(new LinkedHashMap<>(), null, 5));
        assertEquals(5, FqzcompNx16Z.streamStrategy(forced.blobs().get("qualities")));
        // An explicit -1 is the default path, byte for byte.
        GenomicBlocks.BlockBlobs explicit = GenomicBlocks.encodeBlock(block,
                new GenomicWriteContext(new LinkedHashMap<>(), null, -1));
        assertArrayEquals(def.blobs().get("qualities"),
                explicit.blobs().get("qualities"));
    }

    @Test
    void zeroLengthReadForcesRansQualities() throws Exception {
        WrittenGenomicRun run = m87();
        int z = -1;
        for (int i = 0; i < run.readCount(); i++) if (run.lengths()[i] == 0) { z = i; break; }
        org.junit.jupiter.api.Assumptions.assumeTrue(z >= 0, "m87 has a SEQ '*' read");
        WrittenGenomicRun block = GenomicBlocks.sliceRun(run, z, z + 1);
        GenomicBlocks.BlockBlobs blobs = GenomicBlocks.encodeBlock(block,
                new GenomicWriteContext(new LinkedHashMap<>(), null));
        assertEquals(Enums.Compression.RANS_ORDER0.ordinal(), blobs.codecs().get("qualities"));
    }
}
