/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.assembly.AssemblyGraph;
import global.thalion.ttio.assembly.WrittenAssemblyGraph;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.hdf5.Hdf5File;
import global.thalion.ttio.hdf5.Hdf5Group;
import global.thalion.ttio.io.ProgressSink;
import global.thalion.ttio.providers.Hdf5Provider;
import global.thalion.ttio.providers.StorageProvider;

import java.util.*;

/**
 * Root reader/writer for TTI-O {@code .tio} files. Implements
 * {@link AutoCloseable} for try-with-resources.
 *
 * <p>HDF5 layout: root group has {@code @ttio_format_version} and
 * {@code @ttio_features} attributes. The {@code /study/} group contains
 * runs, identifications, quantifications, and provenance.</p>
 *
 * <p><b>Compound metadata (§6 of format-spec).</b> Writing emits the
 * native HDF5 compound dataset <em>and</em> a JSON attribute mirror
 * (transition-window behavior, §6.4). Reading prefers the JSON attribute
 * because JHI5 1.10 cannot marshal VL-string fields out of a compound;
 * when only the native compound is present, primitive fields are
 * recovered via type projection and VL-string fields decode as empty
 * strings. The mirror is emitted to keep Java-written files fully
 * round-trippable by every implementation.</p>
 *
 * <p><b>API status:</b> Stable. {@code Encryptable} conformance is
 * delivered in slice 41.5 when the encryption manager lands in Java.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOSpectralDataset}, Python
 * {@code ttio.spectral_dataset.SpectralDataset}.</p>
 *
 *
 */
public class SpectralDataset implements
        global.thalion.ttio.protocols.Encryptable,
        AutoCloseable {

    private final StorageProvider provider;  // owning provider
    private final Hdf5File file;             // native handle (kept for
                                              // signature/encryption paths)
    private final FeatureFlags featureFlags;
    private final String title;
    private final String isaInvestigationId;
    private final Map<String, AcquisitionRun> msRuns;
    private final Map<String, GenomicRun> genomicRuns;
    // M98: assembly graphs at /study/assembly_graphs/<name>/,
    // populated post-construction by the open + create paths (the
    // ctor chain predates the field; the encryptedAlgorithm precedent).
    private Map<String, AssemblyGraph> assemblyGraphs = Map.of();
    private final Map<String, global.thalion.ttio.genomics.ReferenceImport> references;
    private final MSImage image;  // null when /study/image_cube absent
    private final RamanImage ramanImage;  // null when /study/raman_image_cube absent
    private final IRImage irImage;  // null when /study/ir_image_cube absent
    private final List<Identification> identifications;
    private final List<Quantification> quantifications;
    private final List<ProvenanceRecord> provenanceRecords;
    // Stage 6 (transport-spec v0.11, Deferral 2): first-class
    // Subject + Sample lists. Eagerly populated from
    // /study/subjects/ + /study/samples/ on open(); empty by default
    // for back-compat with pre-Stage-6 files.
    private final List<Subject> subjects;
    private final List<Sample> samples;
    // Encryptable conformance.
    private global.thalion.ttio.protection.AccessPolicy accessPolicy;
    // root-level encryption state that survives close/reopen.
    // Empty string when the dataset carries no @encrypted root attribute;
    // "aes-256-gcm" when it does. Updated by encryptWithKey and by both
    // readers.
    private String encryptedAlgorithm;

    private SpectralDataset(StorageProvider provider, Hdf5File file,
                            FeatureFlags featureFlags,
                            String title, String isaInvestigationId,
                            Map<String, AcquisitionRun> msRuns,
                            Map<String, GenomicRun> genomicRuns,
                            Map<String, global.thalion.ttio.genomics.ReferenceImport> references,
                            MSImage image,
                            RamanImage ramanImage,
                            IRImage irImage,
                            List<Identification> identifications,
                            List<Quantification> quantifications,
                            List<ProvenanceRecord> provenanceRecords,
                            List<Subject> subjects,
                            List<Sample> samples,
                            String encryptedAlgorithm) {
        this.provider = provider;
        this.file = file;
        this.featureFlags = featureFlags;
        this.title = title;
        this.isaInvestigationId = isaInvestigationId;
        this.msRuns = msRuns;
        this.genomicRuns = genomicRuns != null ? genomicRuns : Map.of();
        this.references = references != null ? references : Map.of();
        this.image = image;
        this.ramanImage = ramanImage;
        this.irImage = irImage;
        this.identifications = identifications;
        this.quantifications = quantifications;
        this.provenanceRecords = provenanceRecords;
        this.subjects = subjects != null ? subjects : List.of();
        this.samples = samples != null ? samples : List.of();
        this.encryptedAlgorithm = encryptedAlgorithm != null ? encryptedAlgorithm : "";
    }

    // Pre-Stage-6 long-form constructor (no subjects/samples). Kept for
    // call-sites that have not yet been widened; forwards with empty
    // Subject + Sample lists.
    private SpectralDataset(StorageProvider provider, Hdf5File file,
                            FeatureFlags featureFlags,
                            String title, String isaInvestigationId,
                            Map<String, AcquisitionRun> msRuns,
                            Map<String, GenomicRun> genomicRuns,
                            Map<String, global.thalion.ttio.genomics.ReferenceImport> references,
                            MSImage image,
                            RamanImage ramanImage,
                            IRImage irImage,
                            List<Identification> identifications,
                            List<Quantification> quantifications,
                            List<ProvenanceRecord> provenanceRecords,
                            String encryptedAlgorithm) {
        this(provider, file, featureFlags, title, isaInvestigationId, msRuns,
                genomicRuns, references, image, ramanImage, irImage,
                identifications, quantifications, provenanceRecords,
                List.of(), List.of(), encryptedAlgorithm);
    }

    // Forwarding constructor for callers that don't pass references (write
    // paths that haven't yet read /study/references/ from disk). The
    // Phase 0 tio-browser wiring populates references in the open paths;
    // create paths leave it empty until the file is reopened.
    private SpectralDataset(StorageProvider provider, Hdf5File file,
                            FeatureFlags featureFlags,
                            String title, String isaInvestigationId,
                            Map<String, AcquisitionRun> msRuns,
                            Map<String, GenomicRun> genomicRuns,
                            List<Identification> identifications,
                            List<Quantification> quantifications,
                            List<ProvenanceRecord> provenanceRecords,
                            String encryptedAlgorithm) {
        this(provider, file, featureFlags, title, isaInvestigationId, msRuns,
                genomicRuns, Map.of(), null, null, null, identifications, quantifications,
                provenanceRecords, encryptedAlgorithm);
    }

    // Pre-M82.3 constructors (kept for callers that don't use genomic_runs).
    private SpectralDataset(StorageProvider provider, Hdf5File file,
                            FeatureFlags featureFlags,
                            String title, String isaInvestigationId,
                            Map<String, AcquisitionRun> msRuns,
                            List<Identification> identifications,
                            List<Quantification> quantifications,
                            List<ProvenanceRecord> provenanceRecords,
                            String encryptedAlgorithm) {
        this(provider, file, featureFlags, title, isaInvestigationId, msRuns,
                Map.of(), Map.of(), null, null, null, identifications, quantifications, provenanceRecords,
                encryptedAlgorithm);
    }

    private SpectralDataset(StorageProvider provider, Hdf5File file,
                            FeatureFlags featureFlags,
                            String title, String isaInvestigationId,
                            Map<String, AcquisitionRun> msRuns,
                            List<Identification> identifications,
                            List<Quantification> quantifications,
                            List<ProvenanceRecord> provenanceRecords) {
        this(provider, file, featureFlags, title, isaInvestigationId, msRuns,
                Map.of(), Map.of(), null, null, null, identifications, quantifications, provenanceRecords, "");
    }

    /** The absolute path of the underlying .tio file (null for in-memory datasets). */
    public String filePath() {
        return file != null ? file.getPath() : null;
    }

    /** the owning storage provider. New call sites should reach
     *  for this instead of the native {@link Hdf5File}. */
    public StorageProvider provider() { return provider; }

    // ── Accessors ───────────────────────────────────────────────────

    /** @return The format feature flags advertised by the file (version + opt-in features). */
    public FeatureFlags featureFlags() { return featureFlags; }

    /** @return The dataset's free-form title attribute; empty when unset. */
    public String title() { return title; }

    /** @return The ISA investigation identifier linking this dataset to an ISA-Tab/JSON bundle; empty when unset. */
    public String isaInvestigationId() { return isaInvestigationId; }

    /** @return Unmodifiable map of mass-spectrometry / vibrational / NMR acquisition runs keyed by run name. */
    public Map<String, AcquisitionRun> msRuns() { return msRuns; }
    /** zero or more named genomic runs. Empty for pre-M82
     *  files; populated when {@code /study/genomic_runs/} is present. */
    public Map<String, GenomicRun> genomicRuns() { return genomicRuns; }

    /** Assembly graphs stored under {@code /study/assembly_graphs/},
     *  keyed by name; empty for files without the
     *  {@code opt_assembly_graph} feature (M98). */
    public Map<String, AssemblyGraph> assemblyGraphs() {
        return assemblyGraphs;
    }

    /**
     * Returns embedded references discovered under
     * {@code /study/references/} on this dataset.
     *
     * <p>Keys are reference URIs (the same string returned by
     * {@link global.thalion.ttio.genomics.GenomicRun#referenceUri()}).
     * Values are fully-materialized {@link
     * global.thalion.ttio.genomics.ReferenceImport} instances ready
     * for diff-based codecs and for inspection in user-facing
     * tooling.
     *
     * <p>Datasets written without embedded references (writer flag
     * {@code embedReference=false}) return an empty map even if
     * {@link global.thalion.ttio.genomics.GenomicRun#referenceUri()}
     * is non-null on individual runs.
     *
     * @return unmodifiable map; never null
     * @since 1.1.0
     */
    public Map<String, global.thalion.ttio.genomics.ReferenceImport> references() {
        return Collections.unmodifiableMap(references);
    }

    /**
     * The embedded {@link Image} for the requested {@code kind}, or
     * {@code null} when that modality's cube is absent.
     *
     * <p>JIT2: replaces the typed {@code image()} / {@code ramanImage()}
     * / {@code irImage()} accessors with a single uniform lookup keyed
     * by {@link Enums.ImageKind}. Callers needing a typed handle cast
     * the result (e.g. {@code (MSImage) ds.imageForKind(ImageKind.MS)}).
     * Returns the same eagerly-materialised, cached value object the
     * removed accessors returned — no re-read occurs per call.</p>
     *
     * @param kind which modality (MS / Raman / IR) to fetch
     * @return the materialised image, or {@code null} when absent
     * @since 1.2.0
     */
    public Image imageForKind(Enums.ImageKind kind) {
        switch (kind) {
            case MS:    return image;
            case RAMAN: return ramanImage;
            case IR:    return irImage;
            default:
                throw new IllegalArgumentException("unknown ImageKind: " + kind);
        }
    }

    /**
     * The embedded images on this dataset keyed by modality, containing
     * only the kinds actually present (non-null).
     *
     * @return an {@link java.util.EnumMap} over the present {@link
     *         Enums.ImageKind}s; empty when the dataset carries no images
     * @since 1.2.0
     */
    public java.util.Map<Enums.ImageKind, Image> images() {
        java.util.Map<Enums.ImageKind, Image> out =
            new java.util.EnumMap<>(Enums.ImageKind.class);
        for (Enums.ImageKind k : Enums.ImageKind.values()) {
            Image img = imageForKind(k);
            if (img != null) out.put(k, img);
        }
        return out;
    }

    // ── Phase 2 (post-M91) — canonical unified runs accessor ────────

    /** Phase 2: canonical mapping over every run in the file (MS +
     *  genomic), keyed by run name. Values conform to the
     *  {@link global.thalion.ttio.protocols.Run} interface so callers
     *  can iterate uniformly without knowing the underlying modality:
     *
     *  <pre>{@code
     *  for (var entry : ds.runs().entrySet()) {
     *      Run run = entry.getValue();
     *      System.out.println(run.name() + ": " + run.count() + " measurements");
     *  }
     *  }</pre>
     *
     *  <p>Use {@link #runsOfModality(Class)} to narrow by class, or
     *  {@link #runsForSample(String)} to filter by provenance sample
     *  URI. Phase 2 promotes this to the canonical access pattern;
     *  the legacy {@link #msRuns()} / {@link #genomicRuns()} maps
     *  continue to work, but new code should prefer {@code runs()}.</p>
     *
     *  <p>NMR runs are reported alongside MS runs because the Java
     *  implementation does not split them on disk —
     *  {@link AcquisitionRun} carries both modalities, and
     *  {@code msRuns} already covers both.</p> */
    public Map<String, global.thalion.ttio.protocols.Run> runs() {
        Map<String, global.thalion.ttio.protocols.Run> merged =
            new LinkedHashMap<>();
        for (var entry : msRuns.entrySet()) {
            merged.put(entry.getKey(), entry.getValue());
        }
        for (var entry : genomicRuns.entrySet()) {
            // First-write-wins, matching Python's ``setdefault`` semantics.
            merged.putIfAbsent(entry.getKey(), entry.getValue());
        }
        return merged;
    }

    /** Phase 1 (post-M91): every run associated with {@code sampleUri}.
     *  A run is considered associated when its
     *  {@link global.thalion.ttio.protocols.Run#provenanceChain
     *  provenanceChain} carries {@code sampleUri} in any record's
     *  {@link ProvenanceRecord#inputRefs}. Walks all modalities (MS,
     *  NMR, genomic) uniformly via the Run interface — closes the M91
     *  cross-modality query gap that previously had to fork on
     *  access pattern.
     *
     *  <p>Returns a map keyed by run name; empty when no run matches.
     *  Iteration order is the unified order of {@link #runs()}.</p> */
    public Map<String, global.thalion.ttio.protocols.Run> runsForSample(
            String sampleUri) {
        Map<String, global.thalion.ttio.protocols.Run> out =
            new LinkedHashMap<>();
        for (var entry : runs().entrySet()) {
            global.thalion.ttio.protocols.Run run = entry.getValue();
            List<ProvenanceRecord> chain;
            try {
                chain = run.provenanceChain();
            } catch (Exception e) {
                continue;
            }
            if (chain == null) continue;
            for (ProvenanceRecord r : chain) {
                if (r.inputRefs().contains(sampleUri)) {
                    out.put(entry.getKey(), run);
                    break;
                }
            }
        }
        return out;
    }

    /** Phase 1 (post-M91): every run whose value is an instance of
     *  {@code runType}. Pass {@link AcquisitionRun}{@code .class} to
     *  get the union of MS + NMR runs (any spectrum-class subtype);
     *  pass {@link GenomicRun}{@code .class} to get genomic only. The
     *  return is a thin filter over {@link #runs()}. */
    public Map<String, global.thalion.ttio.protocols.Run> runsOfModality(
            Class<? extends global.thalion.ttio.protocols.Run> runType) {
        Map<String, global.thalion.ttio.protocols.Run> out =
            new LinkedHashMap<>();
        for (var entry : runs().entrySet()) {
            if (runType.isInstance(entry.getValue())) {
                out.put(entry.getKey(), entry.getValue());
            }
        }
        return out;
    }

    /** @return Unmodifiable list of dataset-level identification records. */
    public List<Identification> identifications() { return identifications; }

    /** @return Unmodifiable list of dataset-level quantification records. */
    public List<Quantification> quantifications() { return quantifications; }

    /** @return Unmodifiable list of dataset-level provenance records (per-run provenance lives on each {@link AcquisitionRun}). */
    public List<ProvenanceRecord> provenanceRecords() { return provenanceRecords; }

    /** Stage 6 (transport-spec v0.11, Deferral 2): every {@link Subject}
     *  persisted under {@code /study/subjects/} on this dataset, in
     *  on-disk iteration order. Empty list when no Subjects were
     *  written (which is most pre-Stage-6 files). Backed by an
     *  eagerly-read snapshot taken when {@link #open} runs.
     *  @since 1.4.0 */
    public List<Subject> subjects() { return Collections.unmodifiableList(subjects); }

    /** Stage 6 (transport-spec v0.11, Deferral 2): every {@link Sample}
     *  persisted under {@code /study/samples/} on this dataset, in
     *  on-disk iteration order. Empty list when no Samples were
     *  written. Backed by an eagerly-read snapshot taken when
     *  {@link #open} runs.
     *  @since 1.4.0 */
    public List<Sample> samples() { return Collections.unmodifiableList(samples); }

    /** {@code true} iff this dataset carries an {@code @encrypted} root
     *  attribute. Survives close/reopen because the value is read back
     *  from disk by {@link #open}. Mirrors Python
     *  {@code SpectralDataset.is_encrypted} and Objective-C
     *  {@code -[TTIOSpectralDataset isEncrypted]}. */
    public boolean isEncrypted() { return !encryptedAlgorithm.isEmpty(); }

    /** Algorithm string stored in the root {@code @encrypted} attribute,
     *  or the empty string when the dataset is not encrypted. Typical
     *  value is {@code "aes-256-gcm"}. Mirrors Python
     *  {@code SpectralDataset.encrypted_algorithm} and Objective-C
     *  {@code -[TTIOSpectralDataset encryptedAlgorithm]}. */
    public String encryptedAlgorithm() { return encryptedAlgorithm; }

    // ── Open (read) ─────────────────────────────────────────────────

    /** Open an existing .tio file for reading. (Java):
     *  a URL scheme ({@code memory://}, {@code sqlite://},
     *  {@code zarr://}) dispatches to the matching StorageProvider
     *  and reads the whole dataset through the protocol
     *  (StorageGroup-based); bare paths (and {@code file://} URLs)
     *  stay on the HDF5 fast path for byte parity with pre-M64.5
     *  files and the cross-language smoke suite. */
    public static SpectralDataset open(String pathOrUrl) {
        return open(pathOrUrl, false);
    }

    /** Open an existing .tio file, read-write when {@code writable}
     *  is {@code true} so in-place writers such as
     *  {@code ReferenceImport.writeToDataset} can reach a writable
     *  provider through the open dataset. The file must exist in
     *  either mode. Mirrors Python
     *  {@code SpectralDataset.open(path, writable=True)} and
     *  Objective-C {@code +readFromFilePath:writable:error:}.
     *  @since 1.9.0 */
    public static SpectralDataset open(String pathOrUrl, boolean writable) {
        StorageProvider.Mode mode = writable
                ? StorageProvider.Mode.READ_WRITE
                : StorageProvider.Mode.READ;
        if (pathOrUrl != null && isNonHdf5Url(pathOrUrl)) {
            return openViaProvider(pathOrUrl, mode);
        }
        Hdf5Provider provider = (Hdf5Provider) new Hdf5Provider()
                .open(pathOrUrl, mode);
        Hdf5File file = (Hdf5File) provider.nativeHandle();
        try (Hdf5Group root = file.rootGroup()) {
            FeatureFlags flags = FeatureFlags.readFrom(root);
            String encryptedAlg = root.hasAttribute("encrypted")
                    ? root.readStringAttribute("encrypted") : "";

            String title = null;
            String isaId = null;
            Map<String, AcquisitionRun> runs = new LinkedHashMap<>();
            Map<String, GenomicRun> genomicRuns = new LinkedHashMap<>();
            Map<String, AssemblyGraph> assemblyGraphsMap = new LinkedHashMap<>();  // M98
            Map<String, global.thalion.ttio.genomics.ReferenceImport> references =
                    new LinkedHashMap<>();
            List<Identification> idents = List.of();
            List<Quantification> quants = List.of();
            List<ProvenanceRecord> prov = List.of();
            List<Subject> subjects = List.of();
            List<Sample> samples = List.of();
            MSImage image = null;
            RamanImage ramanImage = null;
            IRImage irImage = null;

            if (root.hasChild("study")) {
                try (Hdf5Group study = root.openGroup("study")) {
                    if (study.hasAttribute("title"))
                        title = study.readStringAttribute("title");
                    if (study.hasAttribute("isa_investigation_id"))
                        isaId = study.readStringAttribute("isa_investigation_id");

                    // Read MS runs
                    if (study.hasChild("ms_runs")) {
                        try (Hdf5Group msRunsGroup = study.openGroup("ms_runs")) {
                            if (msRunsGroup.hasAttribute("_run_names")) {
                                String names = msRunsGroup.readStringAttribute("_run_names");
                                for (String rn : names.split(",")) {
                                    String name = rn.strip();
                                    if (!name.isEmpty() && msRunsGroup.hasChild(name)) {
                                        // AcquisitionRun.readFrom takes
                                        // StorageGroup; wrap the raw Hdf5Group.
                                        AcquisitionRun run = AcquisitionRun.readFrom(
                                                Hdf5Provider.adapterForGroup(msRunsGroup), name);
                                        run.setPersistenceContext(pathOrUrl, name);
                                        runs.put(name, run);
                                    }
                                }
                            }
                        }
                    }

                    // read genomic_runs/ when present.
                    if (study.hasChild("genomic_runs")) {
                        try (Hdf5Group gG = study.openGroup("genomic_runs")) {
                            if (gG.hasAttribute("_run_names")) {
                                String names = gG.readStringAttribute("_run_names");
                                var gAdapter = Hdf5Provider.adapterForGroup(gG);
                                for (String rn : names.split(",")) {
                                    String name = rn.strip();
                                    if (!name.isEmpty() && gG.hasChild(name)) {
                                        var rgGroup = gAdapter.openGroup(name);
                                        genomicRuns.put(name,
                                            GenomicRun.readFrom(rgGroup, name));
                                    }
                                }
                            }
                        }
                    }

                    // /study/references/<uri>/ each holds an embedded reference
                    // (1.1.0 — Phase 0 tio-browser read-back).
                    if (study.hasChild("references")) {
                        try (Hdf5Group refsGroup = study.openGroup("references")) {
                            var refsAdapter = Hdf5Provider.adapterForGroup(refsGroup);
                            for (String uri : refsGroup.childNames()) {
                                try (var oneRef = refsAdapter.openGroup(uri)) {
                                    references.put(uri,
                                        global.thalion.ttio.genomics.ReferenceImport
                                            .readFromGroup(oneRef));
                                }
                            }
                        }
                    }

                    // /study/image_cube — eagerly materialise into a value object (1.2.0).
                    if (study.hasChild("image_cube")) {
                        image = MSImage.readFrom(
                            global.thalion.ttio.providers.Hdf5Provider
                                .adapterForGroup(study));
                    }
                    // /study/raman_image_cube — eagerly materialise (1.2.0).
                    if (study.hasChild("raman_image_cube")) {
                        ramanImage = RamanImage.readFrom(
                            global.thalion.ttio.providers.Hdf5Provider
                                .adapterForGroup(study));
                    }
                    // /study/ir_image_cube — eagerly materialise (1.2.0).
                    if (study.hasChild("ir_image_cube")) {
                        irImage = IRImage.readFrom(
                            global.thalion.ttio.providers.Hdf5Provider
                                .adapterForGroup(study));
                    }
                    idents = SpectralDatasetMetadataIO.readIdentifications(study);
                    quants = SpectralDatasetMetadataIO.readQuantifications(study);
                    prov = SpectralDatasetMetadataIO.readProvenance(study);
                    // Stage 6: per-row subject + sample groups
                    // ({@code /study/subjects/<external_id>/} +
                    // {@code /study/samples/<sample_id>/}). Mirrors how
                    // {@code /study/references/<uri>/} is read above.
                    subjects = SpectralDatasetMetadataIO.readSubjects(study);
                    samples = SpectralDatasetMetadataIO.readSamples(study);
                    // M98: provider-agnostic assembly_graphs read.
                    assemblyGraphsMap.putAll(AssemblyGraph.readAll(
                        Hdf5Provider.adapterForGroup(study)));
                }
            }

            SpectralDataset ds = new SpectralDataset(provider, file, flags, title, isaId, runs,
                    genomicRuns, references, image, ramanImage, irImage,
                    idents, quants, prov, subjects, samples, encryptedAlg);
            ds.assemblyGraphs = assemblyGraphsMap;
            return ds;
        }
    }

    // ── URL-scheme detection () ───────────────────────────

    private static final java.util.regex.Pattern NON_HDF5_URL =
            java.util.regex.Pattern.compile("^(memory|sqlite|zarr)://.*");

    private static boolean isNonHdf5Url(String pathOrUrl) {
        return NON_HDF5_URL.matcher(pathOrUrl).matches();
    }

    // ── Provider-aware read path () ───────────────────────

    private static SpectralDataset openViaProvider(String url,
            StorageProvider.Mode mode) {
        StorageProvider provider = global.thalion.ttio.providers
                .ProviderRegistry.open(url, mode);
        try (global.thalion.ttio.providers.StorageGroup root =
                provider.rootGroup()) {
            FeatureFlags flags = FeatureFlags.readFrom(root);
            String encryptedAlg = "";
            if (root.hasAttribute("encrypted")) {
                Object v = root.getAttribute("encrypted");
                if (v != null) encryptedAlg = v.toString();
            }
            String title = null, isaId = null;
            Map<String, AcquisitionRun> runs = new LinkedHashMap<>();
            Map<String, GenomicRun> genomicRuns = new LinkedHashMap<>();
            Map<String, AssemblyGraph> assemblyGraphsMap = new LinkedHashMap<>();  // M98
            Map<String, global.thalion.ttio.genomics.ReferenceImport> references =
                    new LinkedHashMap<>();
            List<Identification> idents = List.of();
            List<Quantification> quants = List.of();
            List<ProvenanceRecord> prov = List.of();
            List<Subject> subjects = List.of();
            List<Sample> samples = List.of();

            if (root.hasChild("study")) {
                try (global.thalion.ttio.providers.StorageGroup study =
                        root.openGroup("study")) {
                    if (study.hasAttribute("title")) {
                        Object v = study.getAttribute("title");
                        title = v != null ? v.toString() : null;
                    }
                    if (study.hasAttribute("isa_investigation_id")) {
                        Object v = study.getAttribute("isa_investigation_id");
                        isaId = v != null ? v.toString() : null;
                    }
                    if (study.hasChild("ms_runs")) {
                        try (global.thalion.ttio.providers.StorageGroup ms =
                                study.openGroup("ms_runs")) {
                            if (ms.hasAttribute("_run_names")) {
                                Object names = ms.getAttribute("_run_names");
                                String csv = names != null ? names.toString() : "";
                                for (String rn : csv.split(",")) {
                                    String name = rn.strip();
                                    if (!name.isEmpty() && ms.hasChild(name)) {
                                        AcquisitionRun run =
                                                AcquisitionRun.readFrom(ms, name);
                                        run.setPersistenceContext(url, name);
                                        runs.put(name, run);
                                    }
                                }
                            }
                        }
                    }
                    // read genomic_runs/ from any provider.
                    if (study.hasChild("genomic_runs")) {
                        try (var gG = study.openGroup("genomic_runs")) {
                            if (gG.hasAttribute("_run_names")) {
                                Object n = gG.getAttribute("_run_names");
                                String csv = n != null ? n.toString() : "";
                                for (String rn : csv.split(",")) {
                                    String name = rn.strip();
                                    if (!name.isEmpty() && gG.hasChild(name)) {
                                        var rgGroup = gG.openGroup(name);
                                        genomicRuns.put(name,
                                            GenomicRun.readFrom(rgGroup, name));
                                    }
                                }
                            }
                        }
                    }
                    // /study/references/<uri>/ each holds an embedded reference
                    // (1.1.0 — Phase 0 tio-browser read-back).
                    if (study.hasChild("references")) {
                        try (var refsGroup = study.openGroup("references")) {
                            for (String uri : refsGroup.childNames()) {
                                try (var oneRef = refsGroup.openGroup(uri)) {
                                    references.put(uri,
                                        global.thalion.ttio.genomics.ReferenceImport
                                            .readFromGroup(oneRef));
                                }
                            }
                        }
                    }
                    idents = SpectralDatasetMetadataIO.readIdentificationsFromJson(study);
                    quants = SpectralDatasetMetadataIO.readQuantificationsFromJson(study);
                    prov = SpectralDatasetMetadataIO.readProvenanceFromJson(study);
                    // Stage 6: per-row subject + sample groups.
                    subjects = SpectralDatasetMetadataIO.readSubjectsFromProvider(study);
                    samples = SpectralDatasetMetadataIO.readSamplesFromProvider(study);
                    // M98: provider-agnostic assembly_graphs read.
                    assemblyGraphsMap.putAll(AssemblyGraph.readAll(study));
                }
            }
            SpectralDataset ds = new SpectralDataset(provider, null, flags, title, isaId, runs,
                    genomicRuns, references, null, null, null,
                    idents, quants, prov, subjects, samples, encryptedAlg);
            ds.assemblyGraphs = assemblyGraphsMap;
            return ds;
        }
    }

    private static SpectralDataset createViaProviderMixed(
            String url, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            List<String> genomicRunNames,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            List<Subject> subjects,
            List<Sample> samples,
            FeatureFlags featureFlags) {
        return createViaProviderMixed(url, title, isaInvestigationId,
                runs, genomicRuns, genomicRunNames, Map.of(),
                identifications, quantifications, provenanceRecords,
                subjects, samples, featureFlags, () -> {});
    }

    /** Stage D ProgressSink-aware variant of {@link #createViaProviderMixed}.
     *  {@code bumpSection} runs once after each §5.4-ordered section
     *  finishes writing (caller has already emitted the (0, total)
     *  baseline). */
    private static SpectralDataset createViaProviderMixed(
            String url, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            List<String> genomicRunNames,
            Map<String, WrittenAssemblyGraph> assemblyGraphs,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            List<Subject> subjects,
            List<Sample> samples,
            FeatureFlags featureFlags,
            Runnable bumpSection) {
        // Issue #251: genomic writes request the large amortising meta
        // block; pure-spectral writes take HDF5 defaults. Non-HDF5
        // providers (memory/sqlite/zarr) ignore the hint.
        boolean hasGenomic = genomicRuns != null && !genomicRuns.isEmpty();
        StorageProvider provider = global.thalion.ttio.providers
                .ProviderRegistry.open(url, StorageProvider.Mode.CREATE, hasGenomic);
        // Batch all create-time writes into a single provider transaction so
        // SQLite doesn't fsync per group/dataset/attribute. No-op for
        // providers without explicit transactions (default StorageProvider
        // impl).
        provider.beginTransaction();
        try (global.thalion.ttio.providers.StorageGroup root =
                provider.rootGroup()) {
            featureFlags.writeTo(root);
            try (global.thalion.ttio.providers.StorageGroup study =
                    root.createGroup("study")) {
                if (title != null) study.setAttribute("title", title);
                if (isaInvestigationId != null)
                    study.setAttribute("isa_investigation_id", isaInvestigationId);

                Map<String, AcquisitionRun> runMap = new LinkedHashMap<>();
                if (runs != null && !runs.isEmpty()) {
                    try (global.thalion.ttio.providers.StorageGroup ms =
                            study.createGroup("ms_runs")) {
                        StringBuilder names = new StringBuilder();
                        for (int i = 0; i < runs.size(); i++) {
                            AcquisitionRun run = runs.get(i);
                            if (i > 0) names.append(",");
                            names.append(run.name());
                            run.writeTo(ms);
                            runMap.put(run.name(), run);
                        }
                        ms.setAttribute("_run_names", names.toString());
                    }
                    bumpSection.run();  // ms_runs done
                }
                if (identifications != null && !identifications.isEmpty()) {
                    study.setAttribute("identifications_json",
                            SpectralDatasetMetadataIO.buildIdentificationsJson(identifications));
                    bumpSection.run();  // identifications done
                }
                if (quantifications != null && !quantifications.isEmpty()) {
                    study.setAttribute("quantifications_json",
                            buildQuantificationsJson(quantifications));
                    bumpSection.run();  // quantifications done
                }
                if (provenanceRecords != null && !provenanceRecords.isEmpty()) {
                    study.setAttribute("provenance_json",
                            buildProvenanceJson(provenanceRecords));
                    bumpSection.run();  // provenance done
                }

                // genomic_runs subtree (provider-agnostic).
                Map<String, GenomicRun> genomicMap = new LinkedHashMap<>();
                if (genomicRuns != null && !genomicRuns.isEmpty()) {
                    // M93 v1.2: embed references at /study/references/
                    // before writing genomic_runs (provider-agnostic
                    // mirror of the HDF5 fast path).
                    SpectralDatasetGenomicWriter.embedReferencesForRuns(study, genomicRuns);
                    bumpSection.run();  // references done
                    try (var gG = study.createGroup("genomic_runs")) {
                        StringBuilder names = new StringBuilder();
                        for (int i = 0; i < genomicRuns.size(); i++) {
                            WrittenGenomicRun gr = genomicRuns.get(i);
                            String gname = genomicRunNames.get(i);
                            if (i > 0) names.append(",");
                            names.append(gname);
                            writeGenomicRunDefault(study, gG, gname, gr);
                            try (var rgGroup = gG.openGroup(gname)) {
                                genomicMap.put(gname, GenomicRun.readFrom(rgGroup, gname));
                            }
                        }
                        gG.setAttribute("_run_names", names.toString());
                    }
                    bumpSection.run();  // genomic_runs done
                }

                // M98: assembly_graphs subtree (study is already a
                // StorageGroup here). Names are written in sorted
                // order, matching the HDF5 fast path.
                Map<String, AssemblyGraph> assemblyMap = new LinkedHashMap<>();
                if (assemblyGraphs != null && !assemblyGraphs.isEmpty()) {
                    for (String agName : new TreeMap<>(assemblyGraphs).keySet()) {
                        AssemblyGraph.write(assemblyGraphs.get(agName),
                            agName, study);
                    }
                    assemblyMap = AssemblyGraph.readAll(study);
                }

                // Stage 6: per-row subject + sample groups (provider-
                // agnostic). Mirrors the HDF5 fast path; validation
                // already ran upstream in createMixed.
                SpectralDatasetMetadataIO.writeSubjectsViaProvider(study, subjects);
                if (subjects != null && !subjects.isEmpty()) bumpSection.run();
                SpectralDatasetMetadataIO.writeSamplesViaProvider(study, samples);
                if (samples != null && !samples.isEmpty()) bumpSection.run();

                SpectralDataset out = new SpectralDataset(provider, null,
                        featureFlags, title, isaInvestigationId, runMap,
                        genomicMap, Map.of(), null, null, null,
                        identifications != null ? identifications : List.of(),
                        quantifications != null ? quantifications : List.of(),
                        provenanceRecords != null ? provenanceRecords : List.of(),
                        subjects, samples,
                        "");
                out.assemblyGraphs = assemblyMap;
                provider.commitTransaction();
                return out;
            }
        }
    }

    // ── Create (write) ──────────────────────────────────────────────

    /** Sections in §5.4 on-disk write order. Used by the Stage D
     *  per-section ProgressSink reporting. The on-disk write order
     *  (which the writer follows below) puts {@code ms_runs} before
     *  {@code references} so spectral-only files keep their existing
     *  byte layout; the per-section sink reports section indices in
     *  that same emit order so the UI label tracks the actual write
     *  step. */
    private static final List<String> CREATE_SECTION_LABELS = List.of(
        "ms_runs", "references", "genomic_runs",
        "identifications", "quantifications", "provenance",
        "subjects", "samples");

    /** Create a new .tio file with the given content. */
    public static SpectralDataset create(String path, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords) {
        return create(path, title, isaInvestigationId, runs,
                identifications, quantifications, provenanceRecords,
                autoFeatureFlags(runs));
    }

    /**
     * Stage D overload of
     * {@link #create(String, String, String, List, List, List, List)}
     * that fires {@code progress.onProgress(sectionIdx, sectionCount)}
     * once per §5.4-ordered section as it is materialised. UI can drive
     * a determinate progress bar showing "writing identifications..."
     * with section-count progress.
     *
     * <p>The total {@code sectionCount} reflects only sections that
     * will actually be written (empty collections are skipped, so a
     * spectral-only .tio with no identifications / no provenance fires
     * fewer reports than one with the full set).</p>
     *
     * <p>The outer sink does not see byte counts (HDF5 does not expose
     * those per-section cleanly); for per-record progress within a
     * section, consumers should wire a sink directly into the relevant
     * writer (e.g. {@link MzMLWriter#write(AcquisitionRun, String,
     * boolean, ProgressSink)}) — Stage D's writer side has those.</p>
     *
     * @since 1.5.0
     */
    public static SpectralDataset create(String path, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          ProgressSink progress) {
        return createMixed(path, title, isaInvestigationId,
                runs != null ? runs : List.of(),
                List.of(),
                new java.util.ArrayList<String>(),
                identifications, quantifications, provenanceRecords,
                List.of(), List.of(),
                autoFeatureFlags(runs),
                progress != null ? progress : ProgressSink.discard(),
                null, null, null);
    }

    /** Stage 6 (transport-spec v0.11, Deferral 2): create a .tio file
     *  with first-class {@link Subject} + {@link Sample} lists.
     *  Subjects are persisted as {@code /study/subjects/<external_id>/}
     *  per-row groups; Samples as {@code /study/samples/<sample_id>/}.
     *
     *  <p>Validation per spec §4.4:</p>
     *  <ul>
     *    <li>Duplicate {@code Subject.externalId} or
     *        {@code Sample.sampleId} raises
     *        {@link IllegalArgumentException}.</li>
     *    <li>{@code Sample.subjectExternalId} that does not match any
     *        Subject in the same dataset logs a WARNING but does not
     *        fail (anonymous / cross-dataset samples are valid).</li>
     *  </ul>
     *
     *  <p>{@code AcquisitionRun.sampleName} remains the canonical
     *  run → sample link string. Adding Sample rows does not change
     *  that contract.</p>
     *
     *  @since 1.4.0
     */
    public static SpectralDataset create(String path, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          List<Subject> subjects,
                                          List<Sample> samples) {
        return createMixed(path, title, isaInvestigationId,
                runs != null ? runs : List.of(),
                List.of(),
                new java.util.ArrayList<String>(),
                identifications, quantifications, provenanceRecords,
                subjects != null ? subjects : List.of(),
                samples != null ? samples : List.of(),
                autoFeatureFlags(runs));
    }

    /** M74 Slice E: default feature flags upgraded with
     *  {@code opt_ms2_activation_detail} + format version bumped to
     *  {@code "1.3"} whenever any run's {@link SpectrumIndex} carries
     *  the four optional activation/isolation columns. Files without
     *  M74 content keep the legacy 1.1 layout so existing byte-parity
     *  tests continue to pass. */
    private static FeatureFlags autoFeatureFlags(List<AcquisitionRun> runs) {
        FeatureFlags base = FeatureFlags.defaultCurrent();
        if (runs == null) return base;
        boolean anyM74 = false;
        for (AcquisitionRun r : runs) {
            if (r.spectrumIndex() != null
                    && r.spectrumIndex().activationMethods() != null) {
                anyM74 = true;
                break;
            }
        }
        if (!anyM74) return base;
        java.util.Set<String> withFlag = new java.util.LinkedHashSet<>(base.features());
        withFlag.add(FeatureFlags.OPT_MS2_ACTIVATION_DETAIL);
        return new FeatureFlags("1.3", withFlag);
    }

    /** M102: {@code true} when any run's {@link SpectrumIndex} carries
     *  the pixel coordinate columns. */
    private static boolean anyPixelCoordinates(List<AcquisitionRun> runs) {
        if (runs == null) return false;
        for (AcquisitionRun r : runs) {
            if (r != null && r.spectrumIndex() != null
                    && r.spectrumIndex().hasPixelCoordinates()) {
                return true;
            }
        }
        return false;
    }

    /** Convenience overload that delegates to the
     *  {@code genomicRuns}-aware variant with an empty genomic list. */
    public static SpectralDataset create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          FeatureFlags featureFlags) {
        return create(pathOrUrl, title, isaInvestigationId, runs, List.of(),
                identifications, quantifications, provenanceRecords, featureFlags);
    }

    /** Phase 2 (post-M91): mixed-Map create. The {@code runs} map may
     *  carry both {@link AcquisitionRun} (MS / NMR) and
     *  {@link WrittenGenomicRun} (genomic) values; this overload
     *  dispatches by {@code instanceof} on each value and forwards to
     *  the typed-list create API. Mirrors the Python
     *  {@code SpectralDataset.write_minimal} mixed-dict path.
     *
     *  <p>Run-name collision between an MS entry and a genomic entry
     *  raises {@link IllegalArgumentException}. Names are preserved on
     *  disk as-is — the genomic entries no longer get an automatic
     *  {@code genomic_NNNN} prefix when supplied via this overload, so
     *  callers control the storage name.</p>
     *
     *  <p>{@code values} may be empty. Acquired ordering is preserved
     *  (use {@link java.util.LinkedHashMap}). Other parameters mirror
     *  the typed-list overload. */
    public static SpectralDataset create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          Map<String, Object> runs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          FeatureFlags featureFlags) {
        if (runs == null) runs = Map.of();
        List<AcquisitionRun> msList = new ArrayList<>();
        List<WrittenGenomicRun> gList = new ArrayList<>();
        java.util.Set<String> msNames = new java.util.LinkedHashSet<>();
        java.util.Set<String> gNames = new java.util.LinkedHashSet<>();
        for (var entry : runs.entrySet()) {
            String name = entry.getKey();
            Object value = entry.getValue();
            if (value instanceof AcquisitionRun ar) {
                if (gNames.contains(name)) {
                    throw new IllegalArgumentException(
                        "Phase 2 mixed runs map: name '" + name
                        + "' appears as both AcquisitionRun and "
                        + "WrittenGenomicRun");
                }
                msNames.add(name);
                // The on-disk name comes from the AcquisitionRun's own
                // ``name()`` field; reject mismatches early so the
                // caller doesn't silently get a different on-disk name.
                if (!name.equals(ar.name())) {
                    throw new IllegalArgumentException(
                        "Phase 2 mixed runs map: key '" + name
                        + "' does not match AcquisitionRun.name() = '"
                        + ar.name() + "'");
                }
                msList.add(ar);
            } else if (value instanceof WrittenGenomicRun gr) {
                if (msNames.contains(name)) {
                    throw new IllegalArgumentException(
                        "Phase 2 mixed runs map: name '" + name
                        + "' appears as both AcquisitionRun and "
                        + "WrittenGenomicRun");
                }
                gNames.add(name);
                gList.add(gr);
            } else if (value == null) {
                throw new IllegalArgumentException(
                    "Phase 2 mixed runs map: value for '" + name
                    + "' is null");
            } else {
                throw new IllegalArgumentException(
                    "Phase 2 mixed runs map: value for '" + name
                    + "' has unsupported type "
                    + value.getClass().getName()
                    + " (expected AcquisitionRun or WrittenGenomicRun)");
            }
        }
        // Phase 2: the mixed-Map path uses the caller-supplied genomic
        // run names verbatim, bypassing the ``genomic_NNNN`` auto-
        // naming used by the typed-list create. Forward through a
        // private helper so we keep the existing list-based factory
        // intact for back-compat.
        return createMixed(pathOrUrl, title, isaInvestigationId,
                           msList, gList, gNames,
                           identifications, quantifications,
                           provenanceRecords, featureFlags);
    }

    /** full create signature accepting genomic runs
     *  alongside MS runs. When {@code genomicRuns} is non-empty,
     *  {@link FeatureFlags#OPT_GENOMIC} is added (idempotent if the
     *  caller-supplied {@code featureFlags} already includes it) and
     *  the format version is bumped to {@code "1.4"}. */
    public static SpectralDataset create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<WrittenGenomicRun> genomicRuns,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          FeatureFlags featureFlags) {
        return create(pathOrUrl, title, isaInvestigationId,
                runs, genomicRuns, identifications, quantifications,
                provenanceRecords, featureFlags, ProgressSink.discard());
    }

    /**
     * Stage D overload of the genomic-aware
     * {@link #create(String, String, String, List, List, List, List, List, FeatureFlags)}
     * that fires {@code progress.onProgress(sectionIdx, sectionCount)}
     * once per §5.4-ordered section as it is materialised.
     *
     * @since 1.5.0
     */
    public static SpectralDataset create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<WrittenGenomicRun> genomicRuns,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          FeatureFlags featureFlags,
                                          ProgressSink progress) {
        // Phase 2: forward through the names-aware helper with the
        // legacy auto-naming scheme (genomic_NNNN). The mixed-Map create
        // overload calls createMixed directly with caller-supplied names.
        java.util.List<String> autoNames = new java.util.ArrayList<>();
        if (genomicRuns != null) {
            for (int i = 0; i < genomicRuns.size(); i++) {
                autoNames.add("genomic_" + String.format("%04d", i + 1));
            }
        }
        return createMixed(pathOrUrl, title, isaInvestigationId,
                runs != null ? runs : List.of(),
                genomicRuns != null ? genomicRuns : List.of(),
                autoNames,
                identifications, quantifications, provenanceRecords,
                List.of(), List.of(),
                featureFlags,
                progress != null ? progress : ProgressSink.discard(),
                null, null, null);
    }

    /** M98 assembly-graphs-aware create. When {@code assemblyGraphs}
     *  is non-empty, {@link FeatureFlags#OPT_ASSEMBLY_GRAPH} is added
     *  (idempotent) and each graph is written under
     *  {@code /study/assembly_graphs/<name>/} in sorted-name order.
     *  Genomic runs use the legacy {@code genomic_NNNN} auto-naming
     *  like the other typed-list overloads. */
    public static SpectralDataset create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<WrittenGenomicRun> genomicRuns,
                                          Map<String, WrittenAssemblyGraph> assemblyGraphs,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          FeatureFlags featureFlags) {
        java.util.List<String> autoNames = new java.util.ArrayList<>();
        if (genomicRuns != null) {
            for (int i = 0; i < genomicRuns.size(); i++) {
                autoNames.add("genomic_" + String.format("%04d", i + 1));
            }
        }
        return createMixed(pathOrUrl, title, isaInvestigationId,
                runs != null ? runs : List.of(),
                genomicRuns != null ? genomicRuns : List.of(),
                autoNames,
                assemblyGraphs != null ? assemblyGraphs : Map.of(),
                identifications, quantifications, provenanceRecords,
                List.of(), List.of(),
                featureFlags,
                ProgressSink.discard(),
                null, null, null);
    }

    /** JT2: one-shot, image-aware create used by the importer/exporter
     *  registry. This is the single overload that carries the full
     *  normalized draft an importer produces: MS {@code runs} +
     *  {@code genomicRuns} + identifications / quantifications /
     *  provenance + subjects / samples + the three optional embedded
     *  images (MS / Raman / IR) + a {@link ProgressSink}. It forwards
     *  straight to {@link #createMixed} (the widest backend, with image
     *  embedding wired by JT1) using the legacy {@code genomic_NNNN}
     *  auto-naming scheme — matching the other typed-list {@code create}
     *  overloads. The created dataset is closed (flushed to disk) before
     *  returning so the caller may immediately {@link #open(String)} it.
     *
     *  <p>{@link ImportedDataset#write} is the sole call site.
     *  Cross-language equivalent: the Python importer draft's
     *  {@code write(...)}.</p>
     *
     *  @return the {@link java.nio.file.Path} of the written {@code .tio}.
     *  @since 1.7.0 */
    public static java.nio.file.Path create(String pathOrUrl, String title,
                                          String isaInvestigationId,
                                          List<AcquisitionRun> runs,
                                          List<WrittenGenomicRun> genomicRuns,
                                          List<Identification> identifications,
                                          List<Quantification> quantifications,
                                          List<ProvenanceRecord> provenanceRecords,
                                          List<Subject> subjects,
                                          List<Sample> samples,
                                          MSImage image, RamanImage ramanImage,
                                          IRImage irImage,
                                          ProgressSink progress) {
        java.util.List<String> autoNames = new java.util.ArrayList<>();
        if (genomicRuns != null) {
            for (int i = 0; i < genomicRuns.size(); i++) {
                autoNames.add("genomic_" + String.format("%04d", i + 1));
            }
        }
        try (SpectralDataset ds = createMixed(pathOrUrl, title, isaInvestigationId,
                runs != null ? runs : List.of(),
                genomicRuns != null ? genomicRuns : List.of(),
                autoNames,
                identifications, quantifications, provenanceRecords,
                subjects != null ? subjects : List.of(),
                samples != null ? samples : List.of(),
                autoFeatureFlags(runs),
                progress != null ? progress : ProgressSink.discard(),
                image, ramanImage, irImage)) {
            // try-with-resources closes (flushes) the dataset.
        }
        return java.nio.file.Path.of(pathOrUrl);
    }

    /** Back-compat overload for callers that don't pass Subject /
     *  Sample lists. Forwards with empty Stage-6 collections. */
    private static SpectralDataset createMixed(
            String pathOrUrl, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            java.util.Collection<String> genomicRunNames,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            FeatureFlags featureFlags) {
        return createMixed(pathOrUrl, title, isaInvestigationId,
                runs, genomicRuns, genomicRunNames,
                identifications, quantifications, provenanceRecords,
                List.of(), List.of(),
                featureFlags,
                ProgressSink.discard(),
                null, null, null);
    }

    /** Phase 2 (post-M91): names-aware backend used by both the
     *  typed-list {@link #create(String, String, String, List, List,
     *  List, List, List, FeatureFlags) create} (auto-named genomic
     *  runs) and the mixed-Map {@link #create(String, String, String,
     *  Map, List, List, List, FeatureFlags) create} (caller-supplied
     *  genomic names). Stage 6 (transport-spec v0.11) widens with
     *  {@code subjects} + {@code samples}. Kept private — callers go
     *  through one of the public {@code create} overloads. */
    private static SpectralDataset createMixed(
            String pathOrUrl, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            java.util.Collection<String> genomicRunNames,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            List<Subject> subjects,
            List<Sample> samples,
            FeatureFlags featureFlags) {
        return createMixed(pathOrUrl, title, isaInvestigationId,
                runs, genomicRuns, genomicRunNames,
                identifications, quantifications, provenanceRecords,
                subjects, samples, featureFlags,
                ProgressSink.discard(),
                null, null, null);
    }

    /** Stage D ProgressSink-aware overload of {@link #createMixed} that
     *  fires {@code progress.onProgress(sectionIdx, sectionCount)} once
     *  per §5.4-ordered section as it is materialised. Kept with the
     *  pre-M98 signature; forwards with no assembly graphs. */
    private static SpectralDataset createMixed(
            String pathOrUrl, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            java.util.Collection<String> genomicRunNames,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            List<Subject> subjects,
            List<Sample> samples,
            FeatureFlags featureFlags,
            ProgressSink progress,
            MSImage image, RamanImage ramanImage, IRImage irImage) {
        return createMixed(pathOrUrl, title, isaInvestigationId,
                runs, genomicRuns, genomicRunNames, Map.of(),
                identifications, quantifications, provenanceRecords,
                subjects, samples, featureFlags, progress,
                image, ramanImage, irImage);
    }

    /** M98 assembly-graphs-aware backend all createMixed overloads
     *  funnel into. When {@code assemblyGraphs} is non-empty,
     *  {@link FeatureFlags#OPT_ASSEMBLY_GRAPH} is added and the
     *  graphs are written under {@code /study/assembly_graphs/} in
     *  sorted-name order after the genomic subtree. */
    private static SpectralDataset createMixed(
            String pathOrUrl, String title, String isaInvestigationId,
            List<AcquisitionRun> runs,
            List<WrittenGenomicRun> genomicRuns,
            java.util.Collection<String> genomicRunNames,
            Map<String, WrittenAssemblyGraph> assemblyGraphs,
            List<Identification> identifications,
            List<Quantification> quantifications,
            List<ProvenanceRecord> provenanceRecords,
            List<Subject> subjects,
            List<Sample> samples,
            FeatureFlags featureFlags,
            ProgressSink progress,
            MSImage image, RamanImage ramanImage, IRImage irImage) {
        if (progress == null) progress = ProgressSink.discard();
        // v1.0 single format-version stamp. Readers gate optional
        // features by the feature-flag list (opt_*), not by version
        // equality.
        boolean hasGenomic = genomicRuns != null && !genomicRuns.isEmpty();
        boolean hasGraphs = assemblyGraphs != null && !assemblyGraphs.isEmpty();
        String targetVersion = "1.0";
        if (hasGenomic) {
            java.util.Set<String> withFlags =
                new java.util.LinkedHashSet<>(featureFlags.features());
            if (!withFlags.contains(FeatureFlags.OPT_GENOMIC)) {
                withFlags.add(FeatureFlags.OPT_GENOMIC);
            }
            featureFlags = new FeatureFlags(targetVersion, withFlags);
        } else if (!targetVersion.equals(featureFlags.formatVersion())) {
            featureFlags = new FeatureFlags(targetVersion, featureFlags.features());
        }
        // M98: opt_assembly_graph advertises assembly-graph content
        // the same way, added only when a graph is present so ms-only
        // files stay byte-identical.
        if (hasGraphs
                && !featureFlags.features().contains(FeatureFlags.OPT_ASSEMBLY_GRAPH)) {
            java.util.Set<String> withFlags =
                new java.util.LinkedHashSet<>(featureFlags.features());
            withFlags.add(FeatureFlags.OPT_ASSEMBLY_GRAPH);
            featureFlags = new FeatureFlags(targetVersion, withFlags);
        }
        // M101: opt_sam_tags when any genomic run carries SAM optional
        // fields (a tags channel); tag-less files are unchanged.
        boolean hasSamTags = false;
        if (hasGenomic) {
            for (WrittenGenomicRun g : genomicRuns) hasSamTags |= g.hasTags();
        }
        if (hasSamTags
                && !featureFlags.features().contains(FeatureFlags.OPT_SAM_TAGS)) {
            java.util.Set<String> withFlags =
                new java.util.LinkedHashSet<>(featureFlags.features());
            withFlags.add(FeatureFlags.OPT_SAM_TAGS);
            featureFlags = new FeatureFlags(targetVersion, withFlags);
        }
        // M102: opt_pixel_coordinates when any run's spectrum_index
        // carries the pixel_x/y/z columns (format-spec §4b). Stamped
        // here, the backend every create path funnels into, so runs
        // without the columns leave the flag list unchanged.
        if (anyPixelCoordinates(runs)
                && !featureFlags.features().contains(FeatureFlags.OPT_PIXEL_COORDINATES)) {
            java.util.Set<String> withFlags =
                new java.util.LinkedHashSet<>(featureFlags.features());
            withFlags.add(FeatureFlags.OPT_PIXEL_COORDINATES);
            featureFlags = new FeatureFlags(targetVersion, withFlags);
        }

        java.util.List<String> gNamesList = genomicRunNames != null
            ? new java.util.ArrayList<>(genomicRunNames) : new java.util.ArrayList<>();
        if (hasGenomic && gNamesList.size() != genomicRuns.size()) {
            throw new IllegalStateException(
                "createMixed: genomicRunNames (" + gNamesList.size()
                + ") does not match genomicRuns (" + genomicRuns.size() + ")");
        }

        // Stage 6: normalise + validate Subject + Sample lists early so
        // that the writer fails fast (before any HDF5 mutation) on
        // duplicate IDs. Soft-FK warnings are emitted after both lists
        // are known so we can compare them.
        List<Subject> subjectsList = subjects != null ? subjects : List.of();
        List<Sample> samplesList = samples != null ? samples : List.of();
        SpectralDatasetMetadataIO.validateSubjectsAndSamples(subjectsList, samplesList);

        // Stage D: pre-compute which §5.4 sections will be written so
        // the per-section ProgressSink reports a determinate total.
        // Emit order matches the writer's actual on-disk emission:
        //   ms_runs → references → genomic_runs → identifications →
        //   quantifications → provenance → subjects → samples.
        // (references is bundled with genomic_runs because the writer
        // calls embedReferencesForRuns immediately before genomic_runs.)
        long sectionTotal = 0L;
        if (runs != null && !runs.isEmpty()) sectionTotal++;
        if (hasGenomic) sectionTotal += 2L;  // references + genomic_runs
        if (identifications != null && !identifications.isEmpty()) sectionTotal++;
        if (quantifications != null && !quantifications.isEmpty()) sectionTotal++;
        if (provenanceRecords != null && !provenanceRecords.isEmpty()) sectionTotal++;
        if (!subjectsList.isEmpty()) sectionTotal++;
        if (!samplesList.isEmpty()) sectionTotal++;
        final long sectionTotalFinal = sectionTotal;
        final ProgressSink sink = progress;
        // Initial fire at (0, total) so listeners can establish a
        // determinate baseline before the first section completes.
        sink.onProgress(0L, sectionTotalFinal);
        final long[] sectionDone = { 0L };
        Runnable bumpSection = () -> {
            sectionDone[0]++;
            sink.onProgress(sectionDone[0], sectionTotalFinal);
        };

        if (pathOrUrl != null && isNonHdf5Url(pathOrUrl)) {
            // JT1: image embedding is currently only wired on the HDF5
            // fast path. Fail fast rather than silently drop images on a
            // memory:// / sqlite:// / zarr:// target.
            if (image != null || ramanImage != null || irImage != null) {
                throw new UnsupportedOperationException(
                    "createMixed: image embedding is only supported for "
                    + "local .tio (HDF5) targets, not provider URL "
                    + pathOrUrl);
            }
            return createViaProviderMixed(pathOrUrl, title, isaInvestigationId,
                    runs, genomicRuns, gNamesList,
                    assemblyGraphs != null ? assemblyGraphs : Map.of(),
                    identifications, quantifications,
                    provenanceRecords, subjectsList, samplesList, featureFlags,
                    bumpSection);
        }
        // Issue #251: only genomic writes need the 8 MB amortising meta
        // block; pure-spectral files take the HDF5 default (small) blocks.
        Hdf5Provider provider = (Hdf5Provider) new Hdf5Provider()
                .open(pathOrUrl, StorageProvider.Mode.CREATE, hasGenomic);
        Hdf5File file = (Hdf5File) provider.nativeHandle();
        try (Hdf5Group root = file.rootGroup()) {
            featureFlags.writeTo(root);

            try (Hdf5Group study = root.createGroup("study")) {
                if (title != null) study.setStringAttribute("title", title);
                if (isaInvestigationId != null)
                    study.setStringAttribute("isa_investigation_id", isaInvestigationId);

                // JT1: write any embedded images via the same writeTo(...)
                // path the tio-browser GUI uses, immediately after the
                // study-level title/isa attributes and BEFORE the §5.4
                // run sections, so image-free datasets are byte-identical
                // (these branches are skipped when all images are null).
                if (image != null)
                    image.writeTo(Hdf5Provider.adapterForGroup(study));
                if (ramanImage != null)
                    ramanImage.writeTo(Hdf5Provider.adapterForGroup(study));
                if (irImage != null)
                    irImage.writeTo(Hdf5Provider.adapterForGroup(study));

                Map<String, AcquisitionRun> runMap = new LinkedHashMap<>();
                if (runs != null && !runs.isEmpty()) {
                    try (Hdf5Group msRunsGroup = study.createGroup("ms_runs")) {
                        StringBuilder names = new StringBuilder();
                        for (int i = 0; i < runs.size(); i++) {
                            AcquisitionRun run = runs.get(i);
                            if (i > 0) names.append(",");
                            names.append(run.name());
                            // writeTo takes StorageGroup; wrap the
                            // raw Hdf5Group via the provider adapter.
                            run.writeTo(Hdf5Provider.adapterForGroup(msRunsGroup));
                            runMap.put(run.name(), run);
                        }
                        msRunsGroup.setStringAttribute("_run_names", names.toString());
                    }
                    bumpSection.run();  // ms_runs done
                }

                // genomic_runs subtree (only when non-empty).
                Map<String, GenomicRun> genomicMap = new LinkedHashMap<>();
                if (hasGenomic) {
                    // M93 v1.2: embed referenced chromosome sequences at
                    // /study/references/<uri>/ before writing genomic
                    // runs so the writer's REF_DIFF dispatch can resolve
                    // the md5 attribute back from disk if needed.
                    SpectralDatasetGenomicWriter.embedReferencesForRuns(
                        Hdf5Provider.adapterForGroup(study), genomicRuns);
                    bumpSection.run();  // references done
                    try (Hdf5Group gRunsGroup = study.createGroup("genomic_runs")) {
                        StringBuilder names = new StringBuilder();
                        for (int i = 0; i < genomicRuns.size(); i++) {
                            WrittenGenomicRun gr = genomicRuns.get(i);
                            String gname = gNamesList.get(i);
                            if (i > 0) names.append(",");
                            names.append(gname);
                            writeGenomicRunDefault(Hdf5Provider.adapterForGroup(study),
                                Hdf5Provider.adapterForGroup(gRunsGroup), gname, gr);
                            // Open a read-side handle to populate genomicMap.
                            try (var gAdapter = Hdf5Provider.adapterForGroup(gRunsGroup);
                                 var rgGroup = gAdapter.openGroup(gname)) {
                                genomicMap.put(gname, GenomicRun.readFrom(rgGroup, gname));
                            }
                        }
                        gRunsGroup.setStringAttribute("_run_names", names.toString());
                    }
                    bumpSection.run();  // genomic_runs done
                }

                // M98: assembly_graphs subtree (only when non-empty),
                // written through the StorageGroup adapter so the
                // layout is identical to the provider path. Names go
                // in sorted order.
                Map<String, AssemblyGraph> assemblyMap = new LinkedHashMap<>();
                if (hasGraphs) {
                    var studyAdapter = Hdf5Provider.adapterForGroup(study);
                    for (String agName : new TreeMap<>(assemblyGraphs).keySet()) {
                        AssemblyGraph.write(assemblyGraphs.get(agName),
                            agName, studyAdapter);
                    }
                    assemblyMap = AssemblyGraph.readAll(studyAdapter);
                }

                if (identifications != null && !identifications.isEmpty()) {
                    SpectralDatasetMetadataIO.writeIdentifications(study, identifications);
                    bumpSection.run();  // identifications done
                }
                if (quantifications != null && !quantifications.isEmpty()) {
                    SpectralDatasetMetadataIO.writeQuantifications(study, quantifications);
                    bumpSection.run();  // quantifications done
                }
                if (provenanceRecords != null && !provenanceRecords.isEmpty()) {
                    SpectralDatasetMetadataIO.writeProvenance(study, provenanceRecords);
                    bumpSection.run();  // provenance done
                }

                // Stage 6 (transport-spec v0.11, Deferral 2): per-row
                // subject + sample groups under /study/subjects/ +
                // /study/samples/. Validation already ran upstream
                // (duplicate-ID raise + soft-FK warning), so the writer
                // here just emits the typed attributes.
                SpectralDatasetMetadataIO.writeSubjects(study, subjectsList);
                if (!subjectsList.isEmpty()) bumpSection.run();  // subjects done
                SpectralDatasetMetadataIO.writeSamples(study, samplesList);
                if (!samplesList.isEmpty()) bumpSection.run();  // samples done

                SpectralDataset ds = new SpectralDataset(provider, file, featureFlags, title, isaInvestigationId,
                        runMap, genomicMap, Map.of(), null, null, null,
                        identifications != null ? identifications : List.of(),
                        quantifications != null ? quantifications : List.of(),
                        provenanceRecords != null ? provenanceRecords : List.of(),
                        subjectsList, samplesList,
                        "");
                ds.assemblyGraphs = assemblyMap;
                return ds;
            }
        }
    }

    /** JT1: one-shot writer for an image-bearing dataset. Writes the
     *  given (nullable) {@link MSImage} / {@link RamanImage} /
     *  {@link IRImage} into {@code /study/} via the same
     *  {@code writeTo(StorageGroup)} path the tio-browser GUI uses, with
     *  no run / metadata sections. The created dataset is closed (and the
     *  HDF5 file flushed to disk) before returning so the caller may
     *  immediately {@link #open(String)} it.
     *
     *  @return the {@link Path} of the written {@code .tio} file.
     *  @since 1.7.0 */
    public static java.nio.file.Path createWithImages(String path, String title,
                                        String isaInvestigationId,
                                        MSImage image, RamanImage ramanImage,
                                        IRImage irImage) {
        FeatureFlags flags = FeatureFlags.defaultCurrent();
        try (SpectralDataset ds = createMixed(path, title, isaInvestigationId,
                List.of(), List.of(), new java.util.ArrayList<String>(),
                List.of(), List.of(), List.of(),
                List.of(), List.of(),
                flags,
                ProgressSink.discard(),
                image, ramanImage, irImage)) {
            // try-with-resources closes (flushes) the dataset.
        }
        return java.nio.file.Path.of(path);
    }


    /** Phase 1: build the JSON array attribute carrying per-run
     *  provenance for a genomic run. Same shape as
     *  {@link global.thalion.ttio.AcquisitionRun#writeProvenance}. */
    /** Write one genomic run: {@code blocks_v1} through the stream writer
     *  unless the run asks for the v1.8 whole-channel layout. */
    static void writeGenomicRunDefault(
            global.thalion.ttio.providers.StorageGroup study,
            global.thalion.ttio.providers.StorageGroup runsGroup,
            String name, WrittenGenomicRun run) {
        if (run.optLegacyWholeChannel()) {
            SpectralDatasetGenomicWriter.writeGenomicRunSubtree(runsGroup, name, run);
            return;
        }
        try (global.thalion.ttio.genomics.GenomicStreamWriter w =
                 new global.thalion.ttio.genomics.GenomicStreamWriter(study, name,
                     global.thalion.ttio.genomics.GenomicStreamWriter.Options.fromRun(run))) {
            w.appendBatch(run);
        }
    }

    static String buildProvenanceJsonArray(
            List<ProvenanceRecord> records) {
        StringBuilder json = new StringBuilder("[");
        for (int i = 0; i < records.size(); i++) {
            if (i > 0) json.append(",");
            ProvenanceRecord r = records.get(i);
            json.append("{\"timestamp_unix\":").append(r.timestampUnix())
                .append(",\"software\":\"").append(
                    r.software().replace("\"", "\\\""))
                .append("\"")
                .append(",\"parameters\":").append(r.parametersJson())
                .append(",\"input_refs\":").append(r.inputRefsJson())
                .append(",\"output_refs\":").append(r.outputRefsJson())
                .append("}");
        }
        return json.append("]").toString();
    }

    /** Thin delegator preserving the {@code SpectralDataset.buildQuantificationsJson}
     *  call surface after the metadata-IO statics moved to
     *  {@link SpectralDatasetMetadataIO} (P3.10). */
    static String buildQuantificationsJson(List<Quantification> quants) {
        return SpectralDatasetMetadataIO.buildQuantificationsJson(quants);
    }

    /** Thin delegator preserving the {@code SpectralDataset.buildProvenanceJson}
     *  call surface after the metadata-IO statics moved to
     *  {@link SpectralDatasetMetadataIO} (P3.10). */
    static String buildProvenanceJson(List<ProvenanceRecord> records) {
        return SpectralDatasetMetadataIO.buildProvenanceJson(records);
    }


    // writeMateInfoSubgroup + writeMateIntField
    // were removed — the mate_info v1 per-field subgroup writer no
    // longer has a code path. mate_info is v2-only (inline_v2 blob)
    // in v1.0+; non-empty runs without the native lib raise
    // IllegalStateException at the call site (see writeGenomicRunSubtree).


















    // ── Compound metadata: identifications ──────────────────────────



    // ── Compound metadata: quantifications ──────────────────────────



    // ── Compound metadata: provenance ───────────────────────────────



    // ── StorageGroup-based JSON metadata () ───────────────







    // ── JSON parsing (attribute fallback path) ──────────────────────





    // ── Stage 6 (transport-spec v0.11, Deferral 2): Subjects + Samples ──













    // ---- Encryptable conformance ----

    /**
     * Encrypt every MS run's intensity channel in place on disk under
     * AES-256-GCM, then mark the root {@code @encrypted} attribute so
     * the encrypted state survives close/reopen.
     *
     * @param key   32-byte AES-256 key material
     * @param level Encryption granularity (per-run, per-dataset, per-AU)
     * @throws Exception on I/O or cipher failure
     */
    @Override
    public void encryptWithKey(byte[] key, global.thalion.ttio.Enums.EncryptionLevel level)
            throws Exception {
        for (var run : msRuns.values()) run.encryptWithKey(key, level);
        markRootEncrypted();
    }

    /** write the root {@code @encrypted} attribute so
     *  {@link #isEncrypted} / {@link #encryptedAlgorithm} survive
     *  close/reopen. For HDF5-backed datasets (the only backend with
     *  per-run on-disk encryption today) this opens the file R/W after
     *  each run has finished its own encrypt pass. Callers must treat
     *  the dataset as logically closed after this — matching the
     *  Objective-C contract where {@code closeFile} precedes encrypt. */
    private void markRootEncrypted() {
        String path = file != null ? file.getPath() : null;
        if (path != null) {
            try (global.thalion.ttio.hdf5.Hdf5File f =
                         global.thalion.ttio.hdf5.Hdf5File.open(path);
                 Hdf5Group root = f.rootGroup()) {
                root.setStringAttribute("encrypted", "aes-256-gcm");
            }
        }
        this.encryptedAlgorithm = "aes-256-gcm";
    }

    /**
     * Decrypt every MS run's intensity channel into an in-memory
     * overlay. <b>Read-only</b>: the on-disk file is NOT modified, the
     * root {@code @encrypted} attribute is left in place, and
     * {@link #isEncrypted()} continues to return {@code true} on this
     * instance and on any reopen.
     *
     * <p>This is the <b>asymmetric</b> counterpart to
     * {@link #encryptWithKey(byte[], global.thalion.ttio.Enums.EncryptionLevel)}
     * (which IS persistent + flag-flipping) by design: in-memory rehydration
     * lets a process read encrypted data without rewriting the file. To
     * fully reverse encryption on disk and clear the {@code @encrypted}
     * attribute, use the static {@link #decryptInPlace(String, byte[])}
     * — close any open instance first.
     *
     * <p><b>Cross-language equivalents:</b> Python
     * {@code SpectralDataset.decrypt_with_key} (same in-memory-only
     * semantics; returns {@code dict[str, bytes]}), Objective-C
     * {@code -[TTIOSpectralDataset decryptWithKey:error:]} (same).
     */
    @Override
    public void decryptWithKey(byte[] key) throws Exception {
        for (var run : msRuns.values()) run.decryptWithKey(key);
    }

    /**
     * v1.1.1: persist-to-disk decrypt. Strips AES-256-GCM encryption
     * from the {@code .tio} file at {@code path}: for every MS run
     * with an encrypted intensity channel, writes plaintext back as
     * {@code intensity_values} and removes the encrypted siblings.
     * Finally clears the root {@code @encrypted} attribute so
     * {@link #isEncrypted} returns {@code false} when the file is
     * reopened.
     *
     * <p>Symmetric with {@link #encryptWithKey(byte[],
     * global.thalion.ttio.Enums.EncryptionLevel)} (which leaves the root
     * attribute set). After this call the file is byte-compatible with
     * the pre-encryption layout.</p>
     *
     * <p>The file must not be held open by another writer.</p>
     *
     * <p><b>Cross-language equivalents:</b> Python
     * {@code SpectralDataset.decrypt_in_place}, Objective-C
     * {@code +[TTIOSpectralDataset decryptInPlaceAtPath:withKey:error:]}.</p>
     */
    public static void decryptInPlace(String path, byte[] key) {
        if (key == null || key.length != 32) {
            throw new IllegalArgumentException("Key must be exactly 32 bytes");
        }
        java.util.List<String> runNames = new java.util.ArrayList<>();
        try (global.thalion.ttio.hdf5.Hdf5File f =
                     global.thalion.ttio.hdf5.Hdf5File.openReadOnly(path);
             Hdf5Group root = f.rootGroup()) {
            if (root.hasChild("study")) {
                try (Hdf5Group study = root.openGroup("study")) {
                    if (study.hasChild("ms_runs")) {
                        try (Hdf5Group msRunsG = study.openGroup("ms_runs")) {
                            runNames.addAll(msRunsG.childNames());
                        }
                    }
                }
            }
        }

        for (String runName : runNames) {
            global.thalion.ttio.protection.EncryptionManager
                    .decryptIntensityChannelInRunInPlace(path, runName, key);
        }

        try (global.thalion.ttio.hdf5.Hdf5File f =
                     global.thalion.ttio.hdf5.Hdf5File.open(path);
             Hdf5Group root = f.rootGroup()) {
            if (root.hasAttribute("encrypted")) {
                root.deleteAttribute("encrypted");
            }
        }
    }

    /**
     * @return The dataset's attached access policy, or {@code null}
     *         when none is set. Typed as {@code Object} to keep the
     *         {@code Encryptable} protocol free of the protection
     *         package dependency.
     */
    @Override
    public Object accessPolicy() { return accessPolicy; }

    /**
     * Attach an access policy to the dataset.
     *
     * @param policy a {@code global.thalion.ttio.protection.AccessPolicy}
     *               instance, or {@code null} to clear
     * @throws ClassCastException when {@code policy} is non-null and
     *         not an {@code AccessPolicy}
     */
    @Override
    public void setAccessPolicy(Object policy) {
        this.accessPolicy = (global.thalion.ttio.protection.AccessPolicy) policy;
    }

    /**
     * Release the underlying storage handles ({@link StorageProvider}
     * preferred, falling back to the {@link Hdf5File} for legacy
     * callers). Idempotent. Wired through {@link AutoCloseable} so
     * try-with-resources works.
     */
    @Override
    public void close() {
        // Runs hold open channel datasets and block views; release them
        // before the file so the native handle count drops to zero.
        for (AcquisitionRun r : msRuns.values()) {
            try { r.close(); } catch (RuntimeException ignored) { }
        }
        for (GenomicRun g : genomicRuns.values()) {
            try { g.close(); } catch (RuntimeException ignored) { }
        }
        // Prefer closing via the provider (owns the native handle); fall
        // back to direct file close for legacy callers that didn't go
        // through Hdf5Provider.
        if (provider != null) {
            provider.close();
        } else if (file != null) {
            file.close();
        }
    }
}

