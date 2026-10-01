package global.thalion.ttio.importers;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Identification;
import global.thalion.ttio.IRImage;
import global.thalion.ttio.MSImage;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.Quantification;
import global.thalion.ttio.RamanImage;
import global.thalion.ttio.Sample;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.Subject;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.io.ProgressSink;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Normalized in-memory draft every importer produces; {@link #write} is
 *  the single dataset-write call site. Cross-language equivalent: Python
 *  {@code ttio.importers.imported_dataset.ImportedDataset}. */
public final class ImportedDataset {
    public String title = "";
    public String isaInvestigationId = "";
    public final List<AcquisitionRun> runs = new ArrayList<>();
    public final List<WrittenGenomicRun> genomicRuns = new ArrayList<>();
    public final List<Identification> identifications = new ArrayList<>();
    public final List<Quantification> quantifications = new ArrayList<>();
    public final List<ProvenanceRecord> provenance = new ArrayList<>();
    public final List<Subject> subjects = new ArrayList<>();
    public final List<Sample> samples = new ArrayList<>();
    public MSImage image;
    public RamanImage ramanImage;
    public IRImage irImage;
    /** Genomic runs delivered as batch streams; written after the
     *  in-memory content through {@link GenomicStreamSource#writeInto}. */
    public final Map<String, GenomicStreamSource> genomicStreams = new LinkedHashMap<>();
    /** Spectral runs delivered as batch streams; written after the
     *  in-memory content through {@link SpectralStreamSource#writeInto}. */
    public final Map<String, SpectralStreamSource> spectralStreams = new LinkedHashMap<>();

    /** Optional write-through delegate. When non-null, {@link #write}
     *  routes the write through it instead of {@link SpectralDataset#create}.
     *  Used by CLI-delegated (write-through) importers — e.g.
     *  {@link BrukerTDFReader} — whose {@code .tio} is produced by an
     *  external process at write time, not assembled in memory. */
    @FunctionalInterface
    public interface WriteDelegate {
        java.nio.file.Path write(java.nio.file.Path output,
                global.thalion.ttio.io.ProgressSink progress)
                throws java.io.IOException;
    }

    /** When non-null, {@link #write} routes through this delegate. */
    public WriteDelegate writeDelegate;

    /** Create a delegate-backed draft for a write-through importer. */
    public static ImportedDataset delegated(WriteDelegate delegate) {
        ImportedDataset d = new ImportedDataset();
        d.writeDelegate = delegate;
        return d;
    }

    /** Write this draft to {@code output} with no progress reporting. */
    public Path write(Path output) throws java.io.IOException {
        return write(output, null);
    }

    /** Write this draft to {@code output} via the single image-aware
     *  {@link SpectralDataset#create} call site, optionally reporting
     *  per-section progress. Returns the written {@code .tio} path.
     *  When a {@link #writeDelegate} is set, the write is routed
     *  through it (write-through importers). */
    public Path write(Path output, ProgressSink progress) throws java.io.IOException {
        if (writeDelegate != null) {
            return writeDelegate.write(output, progress);
        }
        Path written = SpectralDataset.create(
            output.toString(),
            title.isEmpty() ? "imported" : title,
            isaInvestigationId,
            runs, genomicRuns,
            identifications, quantifications, provenance,
            subjects, samples,
            image, ramanImage, irImage,
            progress != null ? progress : ProgressSink.discard());
        if (!genomicStreams.isEmpty() || !spectralStreams.isEmpty()) {
            try (global.thalion.ttio.providers.StorageProvider p =
                     global.thalion.ttio.providers.ProviderRegistry.open(
                         written.toString(), global.thalion.ttio.providers.StorageProvider.Mode.READ_WRITE, "hdf5")) {
                global.thalion.ttio.providers.StorageGroup root = p.rootGroup();
                global.thalion.ttio.providers.StorageGroup study = root.openGroup("study");
                if (!genomicStreams.isEmpty()) {
                    global.thalion.ttio.FeatureFlags flags = global.thalion.ttio.FeatureFlags.readFrom(root);
                    if (!flags.has(global.thalion.ttio.FeatureFlags.OPT_GENOMIC)) {
                        flags.with(global.thalion.ttio.FeatureFlags.OPT_GENOMIC).writeTo(root);
                    }
                }
                for (GenomicStreamSource src : genomicStreams.values()) src.writeInto(study, progress);
                for (SpectralStreamSource src : spectralStreams.values()) src.writeInto(study, progress);
                if (!genomicStreams.isEmpty()) flagSamTags(root, study);
            }
        }
        return written;
    }

    /** Add {@code opt_sam_tags} once a streamed genomic run turned out
     *  to carry SAM tags (M101); the static write could not know. */
    private void flagSamTags(global.thalion.ttio.providers.StorageGroup root,
                             global.thalion.ttio.providers.StorageGroup study) {
        if (!study.hasChild("genomic_runs")) return;
        boolean any = false;
        try (var runs = study.openGroup("genomic_runs")) {
            for (GenomicStreamSource src : genomicStreams.values()) {
                if (!runs.hasChild(src.name())) continue;
                try (var rg = runs.openGroup(src.name())) {
                    if (rg.hasChild("signal_channels")) {
                        try (var sc = rg.openGroup("signal_channels")) {
                            any |= sc.hasChild("tags");
                        }
                    }
                }
            }
        }
        if (!any) return;
        global.thalion.ttio.FeatureFlags flags = global.thalion.ttio.FeatureFlags.readFrom(root);
        if (!flags.has(global.thalion.ttio.FeatureFlags.OPT_SAM_TAGS)) {
            flags.with(global.thalion.ttio.FeatureFlags.OPT_SAM_TAGS).writeTo(root);
        }
    }
}
