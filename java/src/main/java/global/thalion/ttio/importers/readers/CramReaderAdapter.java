/*
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */
package global.thalion.ttio.importers.readers;

import global.thalion.ttio.importers.CramReader;
import global.thalion.ttio.importers.ImportedDataset;
import global.thalion.ttio.importers.Reader;
import global.thalion.ttio.io.ProgressSink;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/** {@link Reader} adapter for CRAM. Mirrors the GUI
 *  {@code ImportTask.importBamLike("CRAM", ...)}: a CRAM import requires
 *  a reference FASTA (taken from {@code opts.get("reference")}), then
 *  builds a {@code WrittenGenomicRun} via {@link CramReader#toGenomicRun}.
 *  Same {@code name} / {@code region} / {@code sample} opts as
 *  {@link BamReaderAdapter}. */
public final class CramReaderAdapter implements Reader {
    @Override
    public ImportedDataset read(List<String> inputs, Map<String, Object> opts,
                                ProgressSink progress) throws IOException {
        Path reference = referencePath(opts.get("reference"));
        if (reference == null) {
            throw new IllegalArgumentException(
                "CRAM import requires a reference FASTA (opts \"reference\")");
        }

        String name = BamReaderAdapter.optString(opts, "name", "genomic_0001");
        String region = BamReaderAdapter.optString(opts, "region", null);
        String sample = BamReaderAdapter.optString(opts, "sample", null);

        ImportedDataset d = new ImportedDataset();
        CramReader r = new CramReader(Path.of(inputs.get(0)), reference);
        d.genomicStreams.put(name, r.stream(name, region, sample, reference,
            StreamOpts.flag(opts, "embed_reference"), StreamOpts.batchReads(opts))
            .withPolicy(StreamOpts.blockReads(opts), StreamOpts.blockBytes(opts),
                        StreamOpts.flag(opts, "legacy_whole_channel"),
                        StreamOpts.groupReads(opts)));
        return d;
    }

    private static Path referencePath(Object ref) {
        if (ref instanceof Path p) return p;
        if (ref instanceof String s && !s.isEmpty()) return Path.of(s);
        return null;
    }
}
