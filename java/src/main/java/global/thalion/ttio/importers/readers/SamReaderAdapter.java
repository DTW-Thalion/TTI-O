/*
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */
package global.thalion.ttio.importers.readers;

import global.thalion.ttio.importers.ImportedDataset;
import global.thalion.ttio.importers.Reader;
import global.thalion.ttio.importers.SamReader;
import global.thalion.ttio.io.ProgressSink;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/** {@link Reader} adapter for SAM. Mirrors the GUI
 *  {@code ImportTask.importBamLike("SAM", ...)}: builds a
 *  {@code WrittenGenomicRun} via {@link SamReader#toGenomicRun}
 *  ({@link SamReader} extends {@code BamReader}). Same opts as
 *  {@link BamReaderAdapter}: {@code name}, {@code region}, {@code sample}. */
public final class SamReaderAdapter implements Reader {
    @Override
    public ImportedDataset read(List<String> inputs, Map<String, Object> opts,
                                ProgressSink progress) throws IOException {
        String name = BamReaderAdapter.optString(opts, "name", "genomic_0001");
        String region = BamReaderAdapter.optString(opts, "region", null);
        String sample = BamReaderAdapter.optString(opts, "sample", null);

        ImportedDataset d = new ImportedDataset();
        SamReader r = new SamReader(Path.of(inputs.get(0)));
        d.genomicStreams.put(name, r.stream(name, region, sample, StreamOpts.referencePath(opts),
            StreamOpts.flag(opts, "embed_reference"), StreamOpts.batchReads(opts))
            .withPolicy(StreamOpts.blockReads(opts), StreamOpts.blockBytes(opts),
                        StreamOpts.flag(opts, "legacy_whole_channel"),
                        StreamOpts.groupReads(opts)));
        return d;
    }
}
