/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import global.thalion.ttio.genomics.GenomicStreamWriter;
import global.thalion.ttio.genomics.LazyReference;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.io.ProgressSink;
import global.thalion.ttio.providers.StorageGroup;

import java.nio.file.Path;
import java.util.Iterator;
import java.util.function.Supplier;

/**
 * A genomic run delivered as consecutive {@link WrittenGenomicRun}
 * batches (see {@link BamReader#iterBatches} and
 * {@link FastqReader#iterBatches}), written by
 * {@link ImportedDataset#write} through {@link GenomicStreamWriter} with
 * bounded memory. {@code referenceFasta} enables REF_DIFF_V2 through a
 * {@link LazyReference}. Python: {@code ttio.importers.GenomicStreamSource}.
 *
 * @param batches      supplies a fresh iterator over the batches
 * @param blockReads   writer block policy; {@code null} = writer default
 * @param blockBytes   writer block policy; {@code null} = writer default
 * @param groupReads   M103: group the reads by sequence before blocking
 *                     (layout {@code blocks_v1_grouped}; no reference)
 */
public record GenomicStreamSource(String name, Supplier<Iterator<WrittenGenomicRun>> batches,
                                  Path referenceFasta, boolean embedReference,
                                  Integer blockReads, Long blockBytes,
                                  boolean optLegacyWholeChannel, boolean groupReads) {

    /** Pre-M103 signature: reads are not grouped. */
    public GenomicStreamSource(String name, Supplier<Iterator<WrittenGenomicRun>> batches,
                               Path referenceFasta, boolean embedReference,
                               Integer blockReads, Long blockBytes,
                               boolean optLegacyWholeChannel) {
        this(name, batches, referenceFasta, embedReference, blockReads, blockBytes,
            optLegacyWholeChannel, false);
    }

    /** Same source with a writer block policy and layout choice
     *  ({@code null} keeps the writer defaults). */
    public GenomicStreamSource withPolicy(Integer blockReads, Long blockBytes, boolean legacy) {
        return new GenomicStreamSource(name, batches, referenceFasta, embedReference,
            blockReads, blockBytes, legacy, groupReads);
    }

    /** As {@link #withPolicy(Integer, Long, boolean)}, also choosing read
     *  grouping (M103). */
    public GenomicStreamSource withPolicy(Integer blockReads, Long blockBytes, boolean legacy,
                                          boolean group) {
        return new GenomicStreamSource(name, batches, referenceFasta, embedReference,
            blockReads, blockBytes, legacy, group);
    }

    /** Write the run into {@code /study} {@code study}; returns the reads
     *  written. {@code progress} sees {@code (readsSoFar, -1)} per batch
     *  and {@code (n, n)} at the end. */
    public long writeInto(StorageGroup study, ProgressSink progress) {
        ProgressSink sink = progress == null ? ProgressSink.discard() : progress;
        LazyReference ref = referenceFasta != null ? new LazyReference(referenceFasta) : null;
        GenomicStreamWriter writer = null;
        long n = 0;
        Iterator<WrittenGenomicRun> it = batches.get();
        try {
            while (it.hasNext()) {
                WrittenGenomicRun batch = it.next();
                if (writer == null) {
                    GenomicStreamWriter.Options o = GenomicStreamWriter.Options.fromRun(batch);
                    if (ref != null) o = o.withReference(ref, embedReference || batch.embedReference());
                    if (blockReads != null || blockBytes != null) {
                        o = o.withBlockPolicy(blockReads != null ? blockReads : o.blockReads(),
                                              blockBytes != null ? blockBytes : o.blockBytes());
                    }
                    if (optLegacyWholeChannel) o = o.withLegacy(true);
                    if (groupReads) o = o.withGroupReads(true);
                    // Same rule as the producer side, so the two halves
                    // of the pipeline budget are sized off one count.
                    writer = new GenomicStreamWriter(
                        study, name, o, global.thalion.ttio.Threads.resolveImportThreads());
                }
                writer.appendBatch(batch);
                n += batch.readCount();
                sink.onProgress(n, -1L);
            }
        } finally {
            if (it instanceof AutoCloseable c) {
                try { c.close(); } catch (Exception ignored) { }
            }
        }
        if (writer != null) writer.close();
        sink.onProgress(n, n);
        return n;
    }
}
