/*
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */
package global.thalion.ttio.importers.readers;

import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.importers.BamReader;
import global.thalion.ttio.importers.ImportedDataset;
import global.thalion.ttio.importers.Reader;
import global.thalion.ttio.io.ProgressSink;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/** {@link Reader} adapter for BAM. Mirrors the GUI
 *  {@code ImportTask.importBamLike("BAM", ...)}: builds a
 *  {@link WrittenGenomicRun} via {@link BamReader#toGenomicRun}.
 *
 *  <p>Opts (matching the Python genomic readers): {@code name}
 *  (default {@code "genomic_0001"}), {@code region}, {@code sample}.</p> */
public final class BamReaderAdapter implements Reader {
    @Override
    public ImportedDataset read(List<String> inputs, Map<String, Object> opts,
                                ProgressSink progress) throws IOException {
        String name = optString(opts, "name", "genomic_0001");
        String region = optString(opts, "region", null);
        String sample = optString(opts, "sample", null);

        ImportedDataset d = new ImportedDataset();
        BamReader r = new BamReader(Path.of(inputs.get(0)));
        d.genomicStreams.put(name, r.stream(name, region, sample, StreamOpts.referencePath(opts),
            StreamOpts.flag(opts, "embed_reference"), StreamOpts.batchReads(opts))
            .withPolicy(StreamOpts.blockReads(opts), StreamOpts.blockBytes(opts),
                        StreamOpts.flag(opts, "legacy_whole_channel"),
                        StreamOpts.groupReads(opts)));
        return d;
    }

    static String optString(Map<String, Object> opts, String key, String dflt) {
        Object v = opts.get(key);
        return (v instanceof String s) ? s : dflt;
    }
}
