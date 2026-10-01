/* TTI-O Java Implementation / Copyright (c) 2026 The Thalion Initiative / SPDX-License-Identifier: LGPL-3.0-or-later */
package global.thalion.ttio.protection;

import global.thalion.ttio.Enums;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.Enums.Precision;
import global.thalion.ttio.FeatureFlags;
import global.thalion.ttio.protection.PerAUEncryption.AUHeaderPlaintext;
import global.thalion.ttio.protection.PerAUEncryption.ChannelSegment;
import global.thalion.ttio.protection.PerAUEncryption.GcmResult;
import global.thalion.ttio.protection.PerAUEncryption.HeaderSegment;
import global.thalion.ttio.providers.CompoundField;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import global.thalion.ttio.SpectralDatasetGenomicWriter;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.BlockTable;
import global.thalion.ttio.genomics.BlockView;
import global.thalion.ttio.genomics.GenomicBlocks;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.GenomicStreamWriter;
import global.thalion.ttio.genomics.GenomicWriteContext;
import global.thalion.ttio.genomics.PackedReference;
import global.thalion.ttio.genomics.WrittenGenomicRun;

import java.io.ByteArrayOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * v1.0 file-level per-Access-Unit encryption orchestrator.
 *
 * <p>Reads plaintext {@code <channel>_values} datasets from an
 * TTI-O file, encrypts each spectrum independently with
 * {@link PerAUEncryption}, and rewrites the file's
 * {@code signal_channels} groups with the
 * {@code <channel>_segments} compound layout from
 * {@code docs/format-spec.md} §9.1. Routes through the
 * {@link StorageProvider} abstraction so any backend supporting
 * VL_BYTES compound fields (HDF5 + Memory today) works.
 *
 * <p>When {@code encryptHeaders} is true, also encrypts the six
 * semantic index arrays into {@code spectrum_index/au_header_segments}
 * and deletes the plaintext children. Offsets + lengths stay
 * plaintext (structural framing, not semantic PHI).
 *
 * <p>Sets {@code opt_per_au_encryption} (and
 * {@code opt_encrypted_au_headers} when applicable) on the root
 * group.
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.encryption_per_au.encrypt_per_au / decrypt_per_au},
 * Objective-C {@code TTIOPerAUFile}.
 *
 *
 */
public final class PerAUFile {

    private PerAUFile() {}

    /** Result of {@link #decryptFile}: per-run, per-channel plaintext
     *  bytes; optional {@code auHeaders} list when
     *  {@code opt_encrypted_au_headers} is set. M90.11: optional
     *  {@code indexPlain} carries the four genomic_index columns
     *  recovered when the file was encrypted with the reserved
     *  {@code "_headers"} key. M90.12: {@code isGenomic} tags the
     *  run as genomic-uint8 (vs MS-float64) so the CLI can emit the
     *  right MPAD v1 dtype code without having to re-open the file. */
    public record DecryptedRun(Map<String, byte[]> channels,
                                 List<AUHeaderPlaintext> auHeaders,
                                 GenomicIndexPlain indexPlain,
                                 boolean isGenomic,
                                 /** M101: a genomic run's SAM tag text, one
                                  *  string per read; null when the run
                                  *  carries no tags. */
                                 List<String> tags) {
        public DecryptedRun(Map<String, byte[]> channels,
                              List<AUHeaderPlaintext> auHeaders,
                              GenomicIndexPlain indexPlain,
                              boolean isGenomic) {
            this(channels, auHeaders, indexPlain, isGenomic, null);
        }

        public DecryptedRun(Map<String, byte[]> channels,
                              List<AUHeaderPlaintext> auHeaders) {
            this(channels, auHeaders, null, false);
        }

        public DecryptedRun(Map<String, byte[]> channels,
                              List<AUHeaderPlaintext> auHeaders,
                              GenomicIndexPlain indexPlain) {
            this(channels, auHeaders, indexPlain, false);
        }
    }

    /** plaintext genomic_index columns recovered from a file
     *  encrypted with the reserved {@code "_headers"} key. */
    public record GenomicIndexPlain(List<String> chromosomes,
                                      long[] positions,
                                      byte[] mappingQualities,
                                      int[] flags) {}

    /** reserved key name in the {@code keyMap} that signals
     *  the caller wants the genomic_index columns encrypted. */
    public static final String HEADERS_KEY_NAME = "_headers";

    /** Encrypt {@code path} in place. */
    public static void encryptFile(String path, byte[] key,
                                     boolean encryptHeaders,
                                     String providerName) {
        if (key.length != 32) {
            throw new IllegalArgumentException(
                "AES-256-GCM key must be 32 bytes, got " + key.length);
        }
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ_WRITE, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);

            int datasetId = 1;
            try (StorageGroup study = root.openGroup("study")) {
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        for (String runName : runNames(msRuns)) {
                            encryptOneRun(msRuns, runName, datasetId, key,
                                            encryptHeaders);
                            datasetId++;
                        }
                    }
                }
                // continue dataset_id_counter into genomic_runs.
                if (study.hasChild("genomic_runs")) {
                    try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                        for (String runName : runNames(gRuns)) {
                            if (isBlocksV1(gRuns, runName)) {
                                encryptBlocksV1Run(study, gRuns, runName,
                                                   datasetId, key);
                            } else {
                                encryptOneGenomicRun(gRuns, runName,
                                                     datasetId, key);
                            }
                            datasetId++;
                        }
                    }
                }
                // M98: assembly graphs; datasetId continues after
                // the genomic runs.
                if (study.hasChild("assembly_graphs")) {
                    try (StorageGroup agRoot =
                            study.openGroup("assembly_graphs")) {
                        for (String agName : runNames(agRoot)) {
                            encryptOneAssemblyGraph(agRoot, agName,
                                                    datasetId, key);
                            datasetId++;
                        }
                    }
                }
            }

            List<String> updatedFeatures = new ArrayList<>(flags.features());
            if (!updatedFeatures.contains(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                updatedFeatures.add(FeatureFlags.OPT_PER_AU_ENCRYPTION);
            }
            if (encryptHeaders
                    && !updatedFeatures.contains(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS)) {
                updatedFeatures.add(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);
            }
            java.util.Collections.sort(updatedFeatures);
            new FeatureFlags(flags.formatVersion(), updatedFeatures).writeTo(root);
        } finally {
            sp.close();
        }
    }

    /** Read-only decrypt. Returns a map keyed by run name. */
    public static Map<String, DecryptedRun> decryptFile(String path, byte[] key,
                                                          String providerName) {
        if (key.length != 32) {
            throw new IllegalArgumentException(
                "AES-256-GCM key must be 32 bytes, got " + key.length);
        }
        Map<String, DecryptedRun> out = new LinkedHashMap<>();
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);
            if (!flags.has(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                throw new IllegalStateException(
                    "file at " + path + " does not carry opt_per_au_encryption");
            }
            boolean headersEncrypted = flags.has(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);

            int datasetId = 1;
            try (StorageGroup study = root.openGroup("study")) {
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        for (String runName : runNames(msRuns)) {
                            out.put(runName, decryptOneRun(msRuns, runName,
                                                             datasetId, key,
                                                             headersEncrypted));
                            datasetId++;
                        }
                    }
                }
                // dataset_id_counter continues into genomic_runs so
                // AAD reconstruction matches the encrypt path exactly.
                if (study.hasChild("genomic_runs")) {
                    try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                        for (String runName : runNames(gRuns)) {
                            out.put(runName, decryptOneGenomicRun(
                                gRuns, runName, datasetId, key));
                            datasetId++;
                        }
                    }
                }
            }
            return out;
        } finally {
            sp.close();
        }
    }

    /**
     * Persist-to-disk per-AU decrypt counterpart to
     * {@link #encryptFile(String, byte[], boolean, String)}.
     *
     * <p>The legacy
     * {@link global.thalion.ttio.SpectralDataset#decryptInPlace(String, byte[])}
     * (mirror of ObjC {@code +[TTIOSpectralDataset decryptInPlaceAtPath:withKey:]})
     * only handles the {@code intensity_values_encrypted} single-dataset layout
     * and is a silent idempotent no-op on per-AU files. This method is the
     * per-AU equivalent.</p>
     *
     * <p>For each MS run with {@code <channel>_segments} under
     * {@code signal_channels}, decrypts each spectrum row with the per-AU GCM
     * scheme (AAD = {@code dataset_id || au_sequence || channel_name}), writes
     * the concatenated plaintext back as {@code <channel>_values} (float64),
     * removes {@code <channel>_segments} and the {@code <channel>_algorithm}
     * attribute. When the file carries {@code opt_encrypted_au_headers} the
     * six plaintext index datasets are restored from
     * {@code au_header_segments} and the encrypted compound is removed.</p>
     *
     * <p>Genomic runs are handled the same way:
     * {@code study/genomic_runs/<name>/signal_channels/<sequences|qualities>_segments}
     * is decrypted (uint8) and written back as the bare
     * {@code <sequences|qualities>} dataset. {@code dataset_id} continues
     * from where the MS loop left off, mirroring {@link #encryptFile}'s
     * AAD numbering exactly.</p>
     *
     * <p>Strips the {@code opt_per_au_encryption} /
     * {@code opt_encrypted_au_headers} feature flags and the root
     * {@code @encrypted} attribute on completion, so subsequent readers
     * see a fully unprotected file. Idempotent on already-plaintext files
     * (returns without throwing).</p>
     *
     * <p>Cross-language equivalents:
     * Python {@code ttio.encryption_per_au.decrypt_per_au_in_place};
     * ObjC {@code +[TTIOPerAUFile decryptFilePathInPlace:withKey:providerName:error:]}.</p>
     */
    public static void decryptFileInPlace(String path, byte[] key,
                                            String providerName) {
        if (key.length != 32) {
            throw new IllegalArgumentException(
                "AES-256-GCM key must be 32 bytes, got " + key.length);
        }
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ_WRITE, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);
            if (!flags.has(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                return; // idempotent: already plaintext
            }
            boolean headersEncrypted = flags.has(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);

            try (StorageGroup study = root.openGroup("study")) {
                int datasetId = 1;
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        for (String runName : runNames(msRuns)) {
                            decryptOneRunInPlace(msRuns, runName, datasetId,
                                                  key, headersEncrypted);
                            datasetId++;
                        }
                    }
                }
                // datasetId continues into genomic_runs (IDs N+1..N+M),
                // matching encryptFile's AAD numbering.
                if (study.hasChild("genomic_runs")) {
                    try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                        for (String runName : runNames(gRuns)) {
                            if (isBlocksV1(gRuns, runName)) {
                                decryptBlocksV1RunInPlace(study, gRuns,
                                                          runName,
                                                          datasetId, key);
                            } else {
                                decryptOneGenomicRunInPlace(gRuns, runName,
                                                              datasetId, key);
                            }
                            datasetId++;
                        }
                    }
                }
                // M98: assembly graphs; datasetId numbering mirrors
                // encryptFile exactly.
                if (study.hasChild("assembly_graphs")) {
                    try (StorageGroup agRoot =
                            study.openGroup("assembly_graphs")) {
                        for (String agName : runNames(agRoot)) {
                            decryptOneAssemblyGraphInPlace(agRoot, agName,
                                                           datasetId, key);
                            datasetId++;
                        }
                    }
                }
            }

            // All MS + genomic segments are now plaintext: strip the per-AU
            // feature flags + the root @encrypted attribute.
            List<String> updated = new ArrayList<>(flags.features());
            updated.remove(FeatureFlags.OPT_PER_AU_ENCRYPTION);
            updated.remove(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);
            java.util.Collections.sort(updated);
            new FeatureFlags(flags.formatVersion(), updated).writeTo(root);
            if (root.hasAttribute("encrypted")) {
                root.deleteAttribute("encrypted");
            }
        } finally {
            sp.close();
        }
    }

    /** Per-AU decrypt one genomic run in place: replace each
     *  {@code <sequences|qualities>_segments} with a bare uint8 dataset
     *  under the channel name (genomic layout has no _values suffix). */
    private static void decryptOneGenomicRunInPlace(StorageGroup gRuns,
                                                      String runName,
                                                      int datasetId,
                                                      byte[] key) {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels")) {
            for (String cname : new String[]{"sequences", "qualities"}) {
                String segName = cname + "_segments";
                if (!sig.hasChild(segName)) continue;
                List<ChannelSegment> segs = readChannelSegments(sig, segName);
                byte[] plain = PerAUEncryption.decryptChannelFromSegments(
                    segs, datasetId, cname, key, 1);
                if (sig.hasChild(cname)) sig.deleteChild(cname);
                try (StorageDataset ds = sig.createDataset(
                        cname, Precision.UINT8, plain.length, 0,
                        Compression.NONE, 0)) {
                    ds.writeAll(plain);
                }
                sig.deleteChild(segName);
                String algAttr = cname + "_algorithm";
                if (sig.hasAttribute(algAttr)) sig.deleteAttribute(algAttr);
            }
        }
    }

    /** Per-AU decrypt one MS run in place: replace each
     *  {@code <channel>_segments} with a plaintext {@code <channel>_values}
     *  and restore the 6 header columns when applicable. */
    private static void decryptOneRunInPlace(StorageGroup msRuns,
                                              String runName,
                                              int datasetId, byte[] key,
                                              boolean headersEncrypted) {
        try (StorageGroup run = msRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("spectrum_index")) {

            String rawNames = (String) sig.getAttribute("channel_names");
            for (String cname : splitNames(rawNames)) {
                String segName = cname + "_segments";
                if (!sig.hasChild(segName)) continue;
                List<ChannelSegment> segs = readChannelSegments(sig, segName);
                byte[] plain = PerAUEncryption.decryptChannelFromSegments(
                    segs, datasetId, cname, key);
                double[] values = leBytesToDoubles(plain);

                String valuesName = cname + "_values";
                if (sig.hasChild(valuesName)) sig.deleteChild(valuesName);
                try (StorageDataset ds = sig.createDataset(
                        valuesName, Precision.FLOAT64, values.length, 0,
                        Compression.NONE, 0)) {
                    ds.writeAll(values);
                }
                sig.deleteChild(segName);
                String algAttr = cname + "_algorithm";
                if (sig.hasAttribute(algAttr)) sig.deleteAttribute(algAttr);
            }

            if (headersEncrypted && idx.hasChild("au_header_segments")) {
                List<HeaderSegment> hdrSegs =
                    readHeaderSegments(idx, "au_header_segments");
                List<AUHeaderPlaintext> rows =
                    PerAUEncryption.decryptHeaderSegments(hdrSegs, datasetId, key);
                int n = rows.size();
                int[] msLevels = new int[n];
                int[] polarities = new int[n];
                int[] pcs = new int[n];
                double[] rts = new double[n];
                double[] pmzs = new double[n];
                double[] bpis = new double[n];
                for (int i = 0; i < n; i++) {
                    AUHeaderPlaintext h = rows.get(i);
                    msLevels[i] = h.msLevel();
                    polarities[i] = h.polarity();
                    pcs[i] = h.precursorCharge();
                    rts[i] = h.retentionTime();
                    pmzs[i] = h.precursorMz();
                    bpis[i] = h.basePeakIntensity();
                }
                writeIndexColumn(idx, "ms_levels", Precision.INT32, msLevels);
                writeIndexColumn(idx, "polarities", Precision.INT32, polarities);
                writeIndexColumn(idx, "precursor_charges", Precision.INT32, pcs);
                writeIndexColumn(idx, "retention_times", Precision.FLOAT64, rts);
                writeIndexColumn(idx, "precursor_mzs", Precision.FLOAT64, pmzs);
                writeIndexColumn(idx, "base_peak_intensities", Precision.FLOAT64, bpis);
                idx.deleteChild("au_header_segments");
            }
        }
    }

    private static void writeIndexColumn(StorageGroup idx, String name,
                                          Precision precision, Object data) {
        if (idx.hasChild(name)) idx.deleteChild(name);
        int length = (data instanceof int[]) ? ((int[]) data).length
                                              : ((double[]) data).length;
        try (StorageDataset ds = idx.createDataset(name, precision, length, 0,
                                                     Compression.NONE, 0)) {
            ds.writeAll(data);
        }
    }

    private static double[] leBytesToDoubles(byte[] b) {
        double[] out = new double[b.length / 8];
        java.nio.ByteBuffer.wrap(b).order(java.nio.ByteOrder.LITTLE_ENDIAN)
            .asDoubleBuffer().get(out);
        return out;
    }

    // ─────────────────────────────────────────── M90.4 region encryption

    /** encrypt genomic signal channels with a per-chromosome
     *  key map. Reads whose chromosome appears in {@code keyMap} get
     *  AES-256-GCM encrypted with that key; reads on chromosomes not
     *  in the map are stored as clear segments (empty IV/tag,
     *  plaintext rides in the ciphertext slot — see
     *  {@link PerAUEncryption#encryptChannelByRegion}).
     *
     *  <p>MS runs are NOT touched — chromosome is a genomic concept.
     *  Use {@link #encryptFile} for MS encryption.
     *
     *  <p>Sets both {@link FeatureFlags#OPT_PER_AU_ENCRYPTION} and
     *  {@link FeatureFlags#OPT_REGION_KEYED_ENCRYPTION} on the root. */
    public static void encryptByRegion(String path,
                                         Map<String, byte[]> keyMap,
                                         String providerName) {
        for (Map.Entry<String, byte[]> e : keyMap.entrySet()) {
            if (e.getValue().length != 32) {
                throw new IllegalArgumentException(
                    "AES-256-GCM key for chromosome '" + e.getKey()
                    + "' must be 32 bytes, got " + e.getValue().length);
            }
        }
        // split off the reserved "_headers" entry. The
        // remaining map drives per-AU signal-channel dispatch.
        Map<String, byte[]> chromosomeKeys = new LinkedHashMap<>();
        byte[] headersKey = null;
        for (Map.Entry<String, byte[]> e : keyMap.entrySet()) {
            if (HEADERS_KEY_NAME.equals(e.getKey())) {
                headersKey = e.getValue();
            } else {
                chromosomeKeys.put(e.getKey(), e.getValue());
            }
        }
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ_WRITE, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);

            try (StorageGroup study = root.openGroup("study")) {
                if (!study.hasChild("genomic_runs")) {
                    return;  // no genomic data — nothing to encrypt
                }
                // Match the dataset_id_counter convention from the MS
                // path: MS runs occupy 1..N, genomic N+1..N+M. Region
                // encryption only touches genomic, but we still walk MS
                // first to advance the counter.
                int nMs = 0;
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        nMs = runNames(msRuns).size();
                    }
                }
                int datasetId = nMs + 1;
                // Signal-channel encryption runs in two cases (M90.11
                // semantics, mirroring Python's run_signal_encrypt):
                //   (a) caller supplied chromosome keys (path)
                //   (b) caller supplied an empty key_map (no-op)
                // The only path that SKIPS signal-channel encryption
                // is the headers-only case (key_map == {"_headers": K}).
                boolean runSignalEncrypt =
                    !chromosomeKeys.isEmpty() || headersKey == null;
                try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                    for (String runName : runNames(gRuns)) {
                        encryptOneGenomicRunByRegion(gRuns, runName,
                                                       datasetId,
                                                       chromosomeKeys,
                                                       headersKey,
                                                       runSignalEncrypt);
                        datasetId++;
                    }
                }
            }

            List<String> updatedFeatures = new ArrayList<>(flags.features());
            // Feature-flag set rules (mirror Python):
            //  * OPT_PER_AU_ENCRYPTION — set whenever signal-channel
            //    encryption ran (chromosome keys present OR empty
            //    key_map no-op path) OR when headers_key is provided.
            //  * OPT_REGION_KEYED_ENCRYPTION — only when at least one
            //    chromosome key was provided.
            //  * OPT_ENCRYPTED_AU_HEADERS — set when "_headers" key
            //    was used (M90.11).
            // The two Python predicates collapse: "chromosomeKeys
            // present OR headersKey is null" covers the M90.4 path
            // (incl. empty key_map no-op) and "headersKey != null"
            // covers the M90.11 headers-only path. Their union is
            // always true once we reach this point, so we
            // unconditionally add the flag.
            if (!updatedFeatures.contains(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                updatedFeatures.add(FeatureFlags.OPT_PER_AU_ENCRYPTION);
            }
            if (!chromosomeKeys.isEmpty()
                    && !updatedFeatures.contains(FeatureFlags.OPT_REGION_KEYED_ENCRYPTION)) {
                updatedFeatures.add(FeatureFlags.OPT_REGION_KEYED_ENCRYPTION);
            }
            if (headersKey != null
                    && !updatedFeatures.contains(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS)) {
                updatedFeatures.add(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);
            }
            java.util.Collections.sort(updatedFeatures);
            new FeatureFlags(flags.formatVersion(), updatedFeatures).writeTo(root);
        } finally {
            sp.close();
        }
    }

    /** decrypt a region-encrypted file using a per-chromosome
     *  key map. Caller may supply a subset of the keys used at
     *  encryption time — clear segments decode without any key, but
     *  encrypted segments whose chromosome key isn't in {@code keyMap}
     *  raise {@link IllegalStateException}.
     *
     *  <p>Returns {@code {runName -> DecryptedRun}}. The MS runs are
     *  decrypted via the standard single-key path inside this function
     *  iff the file also carries MS encryption under the supplied
     *  key — the M90.4 convention is that MS encryption (if any)
     *  uses the standard {@link #encryptFile} entry point first, and
     *  region encryption layers on top for genomic only. */
    public static Map<String, DecryptedRun> decryptByRegion(String path,
            Map<String, byte[]> keyMap, String providerName) {
        Map<String, DecryptedRun> out = new LinkedHashMap<>();
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);
            if (!flags.has(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                throw new IllegalStateException(
                    "file at " + path + " does not carry opt_per_au_encryption");
            }
            // when the file carries opt_encrypted_au_headers,
            // decrypt requires the reserved "_headers" key. Without
            // it we can't even reconstruct the chromosomes column
            // needed for per-AU dispatch on signal channels.
            boolean headersEncrypted =
                flags.has(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);
            byte[] headersKey = keyMap.get(HEADERS_KEY_NAME);
            if (headersEncrypted && headersKey == null) {
                throw new IllegalStateException(
                    "file at " + path + " carries opt_encrypted_au_headers; "
                    + "caller must provide a '_headers' entry in keyMap "
                    + "to decrypt the genomic_index columns");
            }
            Map<String, byte[]> chromosomeKeys = new LinkedHashMap<>();
            for (Map.Entry<String, byte[]> e : keyMap.entrySet()) {
                if (!HEADERS_KEY_NAME.equals(e.getKey())) {
                    chromosomeKeys.put(e.getKey(), e.getValue());
                }
            }

            try (StorageGroup study = root.openGroup("study")) {
                int nMs = 0;
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        nMs = runNames(msRuns).size();
                    }
                }
                if (!study.hasChild("genomic_runs")) {
                    return out;
                }
                int datasetId = nMs + 1;
                try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                    for (String runName : runNames(gRuns)) {
                        out.put(runName, decryptOneGenomicRunByRegion(
                            gRuns, runName, datasetId, chromosomeKeys,
                            headersEncrypted, headersKey));
                        datasetId++;
                    }
                }
            }
            return out;
        } finally {
            sp.close();
        }
    }

    // ────────────────────────────────────────────────── encrypt helpers

    private static void encryptOneRun(StorageGroup msRuns, String runName,
                                        int datasetId, byte[] key,
                                        boolean encryptHeaders) {
        try (StorageGroup run = msRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("spectrum_index")) {

            int[] lengths = readInts(idx, "lengths");
            // offsets is no longer stored on disk by default;
            // synthesize from cumsum(lengths). Pre-v1.10 files have it.
            long[] offsets = idx.hasChild("offsets")
                ? readLongs(idx, "offsets")
                : global.thalion.ttio.genomics.GenomicIndex.offsetsFromLengths(lengths);

            String rawNames = (String) sig.getAttribute("channel_names");
            List<String> channelNames = splitNames(rawNames);

            for (String cname : channelNames) {
                String valuesName = cname + "_values";
                if (!sig.hasChild(valuesName)) continue;
                // FLOAT_DELTA_ZSTD (codec id 17, the MS default since
                // Phase 2): decode to float64 before slicing — the
                // per-AU segment contract is per-spectrum float64,
                // and decrypt writes plain float64 back.
                double[] values;
                try (StorageDataset ds = sig.openDataset(valuesName)) {
                    long codecId = 0L;
                    if (ds.hasAttribute("compression")) {
                        Object v = ds.getAttribute("compression");
                        if (v instanceof Number num) codecId = num.longValue();
                    }
                    if (codecId == global.thalion.ttio.Enums.Compression
                            .FLOAT_DELTA_ZSTD.ordinal()) {
                        values = global.thalion.ttio.codecs.FloatDeltaZstd
                                .decode((byte[]) ds.readAll());
                    } else {
                        values = (double[]) ds.readAll();
                    }
                }
                byte[] bytes = doublesToLeBytes(values);
                List<ChannelSegment> segs = PerAUEncryption.encryptChannelToSegments(
                    bytes, offsets, lengths, datasetId, cname, key);
                writeChannelSegments(sig, cname + "_segments", segs);
                sig.deleteChild(valuesName);
                sig.setAttribute(cname + "_algorithm", "aes-256-gcm");
            }

            if (encryptHeaders) {
                int acqMode = ((Number) getAttrOr(run, "acquisition_mode", 0L)).intValue();
                double[] rts = readDoubles(idx, "retention_times");
                int[] msLevels = readInts(idx, "ms_levels");
                int[] pols = readInts(idx, "polarities");
                double[] pmzs = readDoubles(idx, "precursor_mzs");
                int[] pcs = readInts(idx, "precursor_charges");
                double[] bpis = readDoubles(idx, "base_peak_intensities");

                List<AUHeaderPlaintext> rows = new ArrayList<>(rts.length);
                for (int i = 0; i < rts.length; i++) {
                    rows.add(new AUHeaderPlaintext(acqMode, msLevels[i], pols[i],
                                                     rts[i], pmzs[i], pcs[i], 0.0,
                                                     bpis[i]));
                }
                List<HeaderSegment> segs =
                    PerAUEncryption.encryptHeaderSegments(rows, datasetId, key);
                writeHeaderSegments(idx, "au_header_segments", segs);

                for (String name : new String[]{"retention_times", "ms_levels",
                                                  "polarities", "precursor_mzs",
                                                  "precursor_charges",
                                                  "base_peak_intensities"}) {
                    if (idx.hasChild(name)) idx.deleteChild(name);
                }
            }
        }
    }

    // ────────────────────────────────────────────────── decrypt helpers

    private static DecryptedRun decryptOneRun(StorageGroup msRuns, String runName,
                                                int datasetId, byte[] key,
                                                boolean headersEncrypted) {
        try (StorageGroup run = msRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("spectrum_index")) {
            Map<String, byte[]> channels = new LinkedHashMap<>();
            String rawNames = (String) sig.getAttribute("channel_names");
            for (String cname : splitNames(rawNames)) {
                String segName = cname + "_segments";
                if (!sig.hasChild(segName)) continue;
                List<ChannelSegment> segs = readChannelSegments(sig, segName);
                channels.put(cname,
                    PerAUEncryption.decryptChannelFromSegments(segs, datasetId,
                                                                  cname, key));
            }

            List<AUHeaderPlaintext> auHeaders = null;
            if (headersEncrypted && idx.hasChild("au_header_segments")) {
                List<HeaderSegment> segs = readHeaderSegments(idx,
                                                                "au_header_segments");
                auHeaders = PerAUEncryption.decryptHeaderSegments(segs, datasetId,
                                                                     key);
            }
            return new DecryptedRun(channels, auHeaders);
        }
    }

    // ─────────────────────────────────────── genomic encrypt / decrypt

    /** encrypt one {@code /study/genomic_runs/<name>/} subtree.
     *  Sequences and qualities are uint8 (one byte per logical
     *  element), AAD reuses the standard
     *  {@code dataset_id || au_sequence || channel_name} layout. */
    private static void encryptOneGenomicRun(StorageGroup gRuns, String runName,
                                               int datasetId, byte[] key) {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("genomic_index")) {
            if (sig.hasChild("tags")) {
                // M101: tags are encrypted with the bases; only the
                // blocks_v1 walker carries them, and leaving them in
                // plaintext would expose MD strings.
                throw new IllegalArgumentException(
                    "per-AU encryption: genomic run '" + runName + "' uses the "
                    + "whole-channel layout and carries SAM tags; per-AU "
                    + "protection of tags needs the blocks_v1 layout. Rewrite "
                    + "the run without optLegacyWholeChannel.");
            }

            int[] lengths = readInts(idx, "lengths");
            // offsets is no longer stored on disk by default;
            // synthesize from cumsum(lengths). Pre-v1.10 files have it.
            long[] offsets = idx.hasChild("offsets")
                ? readLongs(idx, "offsets")
                : global.thalion.ttio.genomics.GenomicIndex.offsetsFromLengths(lengths);

            for (String cname : new String[]{"sequences", "qualities"}) {
                if (!sig.hasChild(cname)) continue;
                byte[] plaintext;
                try (StorageDataset ds = sig.openDataset(cname)) {
                    plaintext = (byte[]) ds.readAll();
                }
                List<ChannelSegment> segs =
                    PerAUEncryption.encryptChannelToSegments(
                        plaintext, offsets, lengths, datasetId, cname, key, 1);
                writeChannelSegments(sig, cname + "_segments", segs);
                sig.deleteChild(cname);
                sig.setAttribute(cname + "_algorithm", "aes-256-gcm");
            }
        }
    }

    /** M98: per-AU encrypt one assembly graph's sequences channel
     *  in place. One AU per segment record; offsets / lengths come
     *  from {@code segments/records} ({@code seq_missing} rows have
     *  length 0 and encrypt to empty ciphertext). The stored channel
     *  is codec-encoded ({@code @compression}), so it is decoded to
     *  raw bytes before slicing; {@code decryptOneAssemblyGraphInPlace}
     *  writes the raw channel back. */
    private static void encryptOneAssemblyGraph(StorageGroup agRoot,
                                                 String graphName,
                                                 int datasetId,
                                                 byte[] key) {
        try (StorageGroup g = agRoot.openGroup(graphName)) {
            if (!g.hasChild("segments")) return;
            try (StorageGroup seg = g.openGroup("segments")) {
                if (!seg.hasChild("sequences")
                        || !seg.hasChild("records")) {
                    return;
                }
                byte[] raw;
                try (StorageDataset ds = seg.openDataset("sequences")) {
                    raw = global.thalion.ttio.assembly.AssemblyGraph
                        .decodeBytesChannel(ds);
                }
                List<Map<String, Object>> rows;
                try (StorageDataset recs = seg.openDataset("records")) {
                    rows = recs.readRows();
                }
                long[] offsets = new long[rows.size()];
                int[] lengths = new int[rows.size()];
                for (int i = 0; i < rows.size(); i++) {
                    offsets[i] = ((Number) rows.get(i)
                        .get("seq_offset")).longValue();
                    lengths[i] = (int) ((Number) rows.get(i)
                        .get("length")).longValue();
                }
                List<ChannelSegment> segs =
                    PerAUEncryption.encryptChannelToSegments(
                        raw, offsets, lengths, datasetId, "sequences",
                        key, 1);
                writeChannelSegments(seg, "sequences_segments", segs);
                seg.deleteChild("sequences");
                seg.setAttribute("sequences_algorithm", "aes-256-gcm");
            }
        }
    }

    /** M98: per-AU decrypt one assembly graph in place. The raw
     *  sequences bytes come back as a plain uint8 dataset with no
     *  {@code @compression}; the graph re-emits byte-exactly from
     *  raw bytes. */
    private static void decryptOneAssemblyGraphInPlace(StorageGroup agRoot,
                                                        String graphName,
                                                        int datasetId,
                                                        byte[] key) {
        try (StorageGroup g = agRoot.openGroup(graphName)) {
            if (!g.hasChild("segments")) return;
            try (StorageGroup seg = g.openGroup("segments")) {
                if (!seg.hasChild("sequences_segments")) return;
                List<ChannelSegment> segs =
                    readChannelSegments(seg, "sequences_segments");
                byte[] plain = PerAUEncryption.decryptChannelFromSegments(
                    segs, datasetId, "sequences", key, 1);
                if (seg.hasChild("sequences")) seg.deleteChild("sequences");
                try (StorageDataset ds = seg.createDataset(
                        "sequences", Precision.UINT8, plain.length, 0,
                        Compression.NONE, 0)) {
                    ds.writeAll(plain);
                }
                seg.deleteChild("sequences_segments");
                if (seg.hasAttribute("sequences_algorithm")) {
                    seg.deleteAttribute("sequences_algorithm");
                }
            }
        }
    }

    /** decrypt one genomic run subtree. Returns a
     *  {@link DecryptedRun} whose {@code channels} map carries
     *  {@code "sequences"} and {@code "qualities"} as flat uint8
     *  byte arrays (no element-width unpacking). */
    private static DecryptedRun decryptOneGenomicRun(StorageGroup gRuns,
            String runName, int datasetId, byte[] key) {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels")) {
            Map<String, byte[]> channels = new LinkedHashMap<>();
            for (String cname : new String[]{"sequences", "qualities"}) {
                String segName = cname + "_segments";
                if (!sig.hasChild(segName)) continue;
                List<ChannelSegment> segs = readChannelSegments(sig, segName);
                channels.put(cname,
                    PerAUEncryption.decryptChannelFromSegments(
                        segs, datasetId, cname, key, 1));
            }
            // M101: tags come back as one string per read.
            List<String> tags = null;
            if (sig.hasChild("tags_segments")) {
                List<ChannelSegment> segs = readChannelSegments(sig, "tags_segments");
                tags = splitTagTexts(PerAUEncryption.decryptChannelFromSegments(
                    segs, datasetId, "tags", key, 1), segs);
            }
            return new DecryptedRun(channels, null, null, /* isGenomic */ true, tags);
        }
    }

    /** M90.4 + M90.11: encrypt one genomic run with per-chromosome
     *  dispatch on signal channels and optional reserved-{@code "_headers"}
     *  encryption of the genomic_index columns.
     *
     *  <p>{@code runSignalEncrypt} is {@code false} only in the
     *  M90.11 headers-only case (key_map == {"_headers": K}) where
     *  the caller wants the index columns encrypted but the signal
     *  channels left untouched. */
    private static void encryptOneGenomicRunByRegion(StorageGroup gRuns,
            String runName, int datasetId,
            Map<String, byte[]> chromosomeKeys,
            byte[] headersKey,
            boolean runSignalEncrypt) {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("genomic_index")) {
            if (sig.hasChild("tags")) {
                // M101: region keys do not cover the SAM tags; refuse
                // rather than leave MD strings in plaintext.
                throw new IllegalArgumentException(
                    "per-AU region encryption: genomic run '" + runName + "' carries "
                    + "SAM tags, which region keys do not protect; use encryptFile");
            }

            int[] lengths = readInts(idx, "lengths");
            // offsets is no longer stored on disk by default;
            // synthesize from cumsum(lengths). Pre-v1.10 files have it.
            long[] offsets = idx.hasChild("offsets")
                ? readLongs(idx, "offsets")
                : global.thalion.ttio.genomics.GenomicIndex.offsetsFromLengths(lengths);
            List<String> chromosomes = readChromosomes(idx);

            if (runSignalEncrypt) {
                for (String cname : new String[]{"sequences", "qualities"}) {
                    if (!sig.hasChild(cname)) continue;
                    byte[] plaintext;
                    try (StorageDataset ds = sig.openDataset(cname)) {
                        plaintext = (byte[]) ds.readAll();
                    }
                    List<ChannelSegment> segs =
                        PerAUEncryption.encryptChannelByRegion(
                            plaintext, offsets, lengths, chromosomes,
                            datasetId, cname, chromosomeKeys);
                    writeChannelSegments(sig, cname + "_segments", segs);
                    sig.deleteChild(cname);
                    sig.setAttribute(cname + "_algorithm",
                                      "aes-256-gcm-by-region");
                }
            }

            // encrypt genomic_index columns under the
            // reserved _headers key.
            if (headersKey != null) {
                encryptGenomicIndex(idx, datasetId, headersKey, chromosomes);
            }
        }
    }

    /** M90.4 + M90.11: decrypt one region-encrypted genomic run. */
    private static DecryptedRun decryptOneGenomicRunByRegion(
            StorageGroup gRuns, String runName, int datasetId,
            Map<String, byte[]> chromosomeKeys,
            boolean headersEncrypted, byte[] headersKey) {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("genomic_index")) {
            // decrypt the genomic_index columns first so the
            // per-AU signal-channel dispatch (which needs chromosomes)
            // can proceed even when the source columns were encrypted.
            List<String> chromosomes;
            GenomicIndexPlain indexPlain = null;
            if (headersEncrypted) {
                indexPlain = decryptGenomicIndex(idx, datasetId, headersKey);
                chromosomes = indexPlain.chromosomes();
            } else {
                chromosomes = readChromosomes(idx);
            }
            Map<String, byte[]> channels = new LinkedHashMap<>();
            for (String cname : new String[]{"sequences", "qualities"}) {
                String segName = cname + "_segments";
                if (!sig.hasChild(segName)) continue;
                List<ChannelSegment> segs = readChannelSegments(sig, segName);
                channels.put(cname,
                    PerAUEncryption.decryptChannelByRegion(
                        segs, chromosomes, datasetId, cname, chromosomeKeys));
            }
            return new DecryptedRun(channels, null, indexPlain,
                                      /* isGenomic */ true);
        }
    }

    /** encrypt the four genomic_index columns
     *  (chromosomes, positions, mapping_qualities, flags) and replace
     *  the plaintext datasets with {@code <column>_encrypted} blobs
     *  containing {@code iv || tag || ciphertext}. {@code offsets} /
     *  {@code lengths} stay plaintext (structural framing).
     *
     *  <p>Per-column AES-GCM with AAD =
     *  {@code "genomic_headers:" + datasetId + ":" + column_name}.
     *  Chromosomes (a VL compound) is JSON-serialised before
     *  encryption to match Python's
     *  {@code json.dumps(chromosomes).encode("utf-8")} byte form. */
    private static void encryptGenomicIndex(StorageGroup idx, int datasetId,
                                              byte[] key,
                                              List<String> chromosomes) {
        long[] positions = readLongs(idx, "positions");
        byte[] mapqs;
        try (StorageDataset ds = idx.openDataset("mapping_qualities")) {
            mapqs = (byte[]) ds.readAll();
        }
        int[] flags = readInts(idx, "flags");

        Map<String, byte[]> columns = new LinkedHashMap<>();
        columns.put("chromosomes",
            chromosomesJson(chromosomes).getBytes(StandardCharsets.UTF_8));
        columns.put("positions", longsToLeBytes(positions));
        columns.put("mapping_qualities", mapqs.clone());
        columns.put("flags", intsToLeBytes(flags));

        for (Map.Entry<String, byte[]> e : columns.entrySet()) {
            String colName = e.getKey();
            byte[] plaintext = e.getValue();
            byte[] aad = ("genomic_headers:" + datasetId + ":" + colName)
                .getBytes(StandardCharsets.US_ASCII);
            GcmResult r = PerAUEncryption.encryptWithAad(plaintext, key, aad, null);
            byte[] blob = new byte[r.iv().length + r.tag().length
                                     + r.ciphertext().length];
            System.arraycopy(r.iv(), 0, blob, 0, r.iv().length);
            System.arraycopy(r.tag(), 0, blob, r.iv().length, r.tag().length);
            System.arraycopy(r.ciphertext(), 0, blob,
                              r.iv().length + r.tag().length,
                              r.ciphertext().length);
            if (idx.hasChild(colName)) {
                idx.deleteChild(colName);
            }
            // L1 (Task #82 Phase B.1): the on-disk chromosomes column
            // is now decomposed into chromosome_ids + chromosome_names
            // — also delete those when encrypting the logical
            // "chromosomes" column.
            if ("chromosomes".equals(colName)) {
                if (idx.hasChild("chromosome_ids")) {
                    idx.deleteChild("chromosome_ids");
                }
                if (idx.hasChild("chromosome_names")) {
                    idx.deleteChild("chromosome_names");
                }
            }
            String encName = colName + "_encrypted";
            if (idx.hasChild(encName)) {
                idx.deleteChild(encName);
            }
            try (StorageDataset out = idx.createDataset(encName,
                    Enums.Precision.UINT8, blob.length, 0,
                    Enums.Compression.NONE, 0)) {
                out.writeAll(blob);
            }
        }
    }

    /** inverse of {@link #encryptGenomicIndex}. Returns the
     *  four plaintext columns. */
    private static GenomicIndexPlain decryptGenomicIndex(
            StorageGroup idx, int datasetId, byte[] key) {
        List<String> chromosomes = null;
        long[] positions = null;
        byte[] mapqs = null;
        int[] flags = null;
        for (String colName : new String[]{
                "chromosomes", "positions", "mapping_qualities", "flags"}) {
            String encName = colName + "_encrypted";
            if (!idx.hasChild(encName)) {
                throw new IllegalStateException(
                    "genomic_index/" + encName + " missing — file does not "
                    + "appear to carry M90.11 encrypted headers");
            }
            byte[] blob;
            try (StorageDataset ds = idx.openDataset(encName)) {
                blob = (byte[]) ds.readAll();
            }
            if (blob.length < 12 + 16) {
                throw new IllegalStateException(
                    "genomic_index/" + encName + " too short for IV+TAG");
            }
            byte[] iv = Arrays.copyOfRange(blob, 0, 12);
            byte[] tag = Arrays.copyOfRange(blob, 12, 28);
            byte[] ciphertext = Arrays.copyOfRange(blob, 28, blob.length);
            byte[] aad = ("genomic_headers:" + datasetId + ":" + colName)
                .getBytes(StandardCharsets.US_ASCII);
            byte[] plain = PerAUEncryption.decryptWithAad(iv, tag, ciphertext,
                                                            key, aad);
            switch (colName) {
                case "chromosomes":
                    chromosomes = chromosomesFromJson(
                        new String(plain, StandardCharsets.UTF_8));
                    break;
                case "positions":
                    positions = leBytesToLongs(plain);
                    break;
                case "mapping_qualities":
                    mapqs = plain;
                    break;
                case "flags":
                    flags = leBytesToInts(plain);
                    break;
            }
        }
        return new GenomicIndexPlain(chromosomes, positions, mapqs, flags);
    }

    /** Match Python {@code json.dumps(list_of_str)} with default
     *  separators ({@code ", "} between items). */
    static String chromosomesJson(List<String> chromosomes) {
        StringBuilder sb = new StringBuilder();
        sb.append('[');
        for (int i = 0; i < chromosomes.size(); i++) {
            if (i > 0) sb.append(", ");
            sb.append('"');
            // Chromosome names are simple ASCII identifiers (chr1,
            // chr6, chrX, ...) — no escaping needed in practice. We
            // still escape backslash + double-quote for safety.
            String s = chromosomes.get(i);
            for (int j = 0; j < s.length(); j++) {
                char c = s.charAt(j);
                if (c == '"' || c == '\\') sb.append('\\');
                sb.append(c);
            }
            sb.append('"');
        }
        sb.append(']');
        return sb.toString();
    }

    /** Tiny JSON parser for the chromosomes column — accepts the
     *  flat string-array shape produced by {@link #chromosomesJson}
     *  or by Python's {@code json.dumps}. */
    static List<String> chromosomesFromJson(String json) {
        List<String> out = new ArrayList<>();
        int n = json.length();
        int i = 0;
        // Skip leading whitespace and the opening '['.
        while (i < n && Character.isWhitespace(json.charAt(i))) i++;
        if (i >= n || json.charAt(i) != '[') {
            throw new IllegalStateException(
                "chromosomes JSON must start with '[': " + json);
        }
        i++;
        while (i < n) {
            while (i < n && (Character.isWhitespace(json.charAt(i))
                              || json.charAt(i) == ',')) i++;
            if (i >= n) break;
            char c = json.charAt(i);
            if (c == ']') break;
            if (c != '"') {
                throw new IllegalStateException(
                    "chromosomes JSON expected '\"' at " + i + ": " + json);
            }
            i++;  // past opening quote
            StringBuilder sb = new StringBuilder();
            while (i < n) {
                char ch = json.charAt(i);
                if (ch == '\\' && i + 1 < n) {
                    sb.append(json.charAt(i + 1));
                    i += 2;
                    continue;
                }
                if (ch == '"') {
                    i++;
                    break;
                }
                sb.append(ch);
                i++;
            }
            out.add(sb.toString());
        }
        return out;
    }

    /** Read the genomic_index chromosome columns into a
     *  {@code List<String>}. L1 (Task #82 Phase B.1, 2026-05-01):
     *  chromosomes are stored as {@code chromosome_ids} (uint16) +
     *  {@code chromosome_names} (compound) instead of a single
     *  VL-string compound. */
    @SuppressWarnings("unchecked")
    private static List<String> readChromosomes(StorageGroup idx) {
        short[] ids;
        try (StorageDataset ds = idx.openDataset("chromosome_ids")) {
            ids = (short[]) ds.readAll();
        }
        List<Object[]> nameRows;
        try (StorageDataset ds = idx.openDataset("chromosome_names")) {
            nameRows = (List<Object[]>) ds.readAll();
        }
        List<String> nameTable = new ArrayList<>(nameRows.size());
        for (Object[] r : nameRows) {
            Object v = r[0];
            if (v == null) {
                nameTable.add("");
            } else if (v instanceof byte[] b) {
                nameTable.add(new String(b, java.nio.charset.StandardCharsets.UTF_8));
            } else {
                nameTable.add(v.toString());
            }
        }
        List<String> out = new ArrayList<>(ids.length);
        for (short id : ids) {
            int idx2 = Short.toUnsignedInt(id);
            out.add(idx2 < nameTable.size() ? nameTable.get(idx2) : "");
        }
        return out;
    }

    // ────────────────────────────────────── M99: blocks_v1 walkers
    //
    // The default genomic layout stores codec-coded per-block blobs,
    // so the walkers stream block by block: decode one block, slice
    // its reads, encrypt one AU per read with GLOBAL AU numbering,
    // append to extendable segments tables. Restore re-encodes each
    // block with the stream writer's machinery, honouring the writer
    // policy the run group persists (opt_disable_qualities_v5,
    // ref_diff_slice_bytes); when a re-encoded blob still lands on
    // different ranges (older files without the attrs), the block
    // index is rewritten to the ranges actually written instead of
    // refusing, keeping the file consistent and readable.

    // M101: the SAM tags are encrypted with the bases they describe (an
    // MD string lists the read's mismatches). One AU per read, the
    // plaintext being the read's UTF-8 tag text; "tags" is handled as
    // text, not bytes, and is never injected into a block view.
    private static final List<String> BLOCKS_V1_CHANNELS =
        List.of("sequences", "qualities", "tags");

    /** Split decrypted tag text by the segments' plaintext lengths. */
    private static List<String> splitTagTexts(byte[] flat, List<ChannelSegment> segs) {
        List<String> out = new ArrayList<>(segs.size());
        int p = 0;
        for (ChannelSegment s : segs) {
            out.add(new String(flat, p, s.length(), StandardCharsets.UTF_8));
            p += s.length();
        }
        return out;
    }

    private static boolean isBlocksV1(StorageGroup gRuns, String runName) {
        try (StorageGroup run = gRuns.openGroup(runName)) {
            if (!run.hasAttribute("layout")) return false;
            Object layout = run.getAttribute("layout");
            return layout != null && "blocks_v1".equals(layout.toString());
        }
    }

    /** {@code {chromosome: bytes}} of an embedded reference; null when
     *  absent. */
    private static Map<String, byte[]> embeddedReferenceSeqs(
            StorageGroup study, String uri) {
        if (uri == null || uri.isEmpty()
                || !study.hasChild("references")) return null;
        try (StorageGroup refs = study.openGroup("references")) {
            if (!refs.hasChild(uri)) return null;
            try (StorageGroup ref = refs.openGroup(uri)) {
                if (!ref.hasChild("chromosomes")) return null;
                try (StorageGroup chroms = ref.openGroup("chromosomes")) {
                    Map<String, byte[]> out = new LinkedHashMap<>();
                    for (String cname : chroms.childNames()) {
                        try (StorageGroup cg = chroms.openGroup(cname)) {
                            out.put(cname,
                                    PackedReference.readChromosomeBytes(cg));
                        }
                    }
                    return out;
                }
            }
        }
    }

    /** The run-wide chromosome-id map the writer accumulated: the
     *  {@code mate_info/chrom_names} table in row order. */
    private static Map<String, Integer> blocksV1ChromMap(StorageGroup sig) {
        Map<String, Integer> map = new LinkedHashMap<>();
        if (sig.hasChild("mate_info")) {
            try (StorageGroup mate = sig.openGroup("mate_info")) {
                List<String> names = BlockView.readNames(mate, "chrom_names");
                for (int i = 0; i < names.size(); i++) map.put(names.get(i), i);
            }
        }
        return map;
    }

    /** {@code {chromosome: bytes}} for re-encoding a run's REF_DIFF
     *  blocks. With the {@code @reference_md5s} run attribute the
     *  reference is rebuilt through the resolver (embedded or
     *  {@code REF_PATH}); older files without the attribute fall back
     *  to the embedded-directory walk. Null when neither works. */
    private static Map<String, byte[]> blocksV1ReferenceSeqs(
            StorageGroup runGroup, StorageGroup study, String uri) {
        String raw = null;
        if (runGroup.hasAttribute("reference_md5s")) {
            Object v = runGroup.getAttribute("reference_md5s");
            raw = v == null ? null
                : v instanceof byte[] b
                    ? new String(b, java.nio.charset.StandardCharsets.UTF_8)
                    : v.toString();
        }
        if (raw != null && !raw.isEmpty()) {
            Map<String, String> md5s;
            try {
                md5s = new com.fasterxml.jackson.databind.ObjectMapper()
                    .readValue(raw,
                        new com.fasterxml.jackson.core.type
                            .TypeReference<Map<String, String>>() {});
            } catch (com.fasterxml.jackson.core.JacksonException e) {
                throw new IllegalStateException(
                    "per-AU blocks_v1: run attribute reference_md5s "
                    + "is not a JSON object of chromosome to md5 hex: "
                    + raw, e);
            }
            global.thalion.ttio.hdf5.Hdf5Group h5g =
                global.thalion.ttio.providers.Hdf5Provider
                    .tryUnwrapHdf5Group(runGroup);
            if (h5g == null) return null;
            var resolver = new global.thalion.ttio.codecs
                .ReferenceResolver(h5g.owningFile());
            Map<String, byte[]> out = new LinkedHashMap<>();
            for (Map.Entry<String, String> e : md5s.entrySet()) {
                byte[] md5 = new byte[16];
                String hex = e.getValue();
                for (int i = 0; i < 16; i++) {
                    md5[i] = (byte) Integer.parseInt(
                        hex.substring(i * 2, i * 2 + 2), 16);
                }
                out.put(e.getKey(),
                        resolver.resolve(uri, md5, e.getKey()));
            }
            return out;
        }
        return embeddedReferenceSeqs(study, uri);
    }

    private record BlockRun(WrittenGenomicRun run, byte[] seq, byte[] qual) {}

    private static BlockRun blocksV1BlockRun(GenomicRun rd,
            StorageGroup runGroup, StorageGroup study, BlockTable t,
            int b, long indexBase) {
        return blocksV1BlockRun(rd, runGroup, study, t, b, indexBase, null);
    }

    /** Collect reads {@code [indexBase, indexBase+nn)} from an open
     *  reader into a per-block {@link WrittenGenomicRun}. The encrypt
     *  walker reads block {@code b} through the run's own reader
     *  ({@code indexBase = readStartAt(b)}); the decrypt walker reads
     *  a materialised one-block view ({@code indexBase = 0}). */
    /** {@code tags}, when given, are the block's decrypted tag texts
     *  (M101); otherwise each read's tags come from the reader. */
    private static BlockRun blocksV1BlockRun(GenomicRun rd,
            StorageGroup runGroup, StorageGroup study, BlockTable t,
            int b, long indexBase, List<String> tags) {
        int nn = t.nReadsAt(b);
        List<String> readTags = new ArrayList<>(nn);
        long[] positions = new long[nn];
        byte[] mapqs = new byte[nn];
        int[] flags = new int[nn];
        long[] offsets = new long[nn];
        int[] lengths = new int[nn];
        long[] matePos = new long[nn];
        int[] tlens = new int[nn];
        List<String> cigars = new ArrayList<>(nn);
        List<String> names = new ArrayList<>(nn);
        List<String> mateChroms = new ArrayList<>(nn);
        List<String> chroms = new ArrayList<>(nn);
        ByteArrayOutputStream seq = new ByteArrayOutputStream();
        ByteArrayOutputStream qual = new ByteArrayOutputStream();
        for (int i = 0; i < nn; i++) {
            AlignedRead r = rd.readAt((int) (indexBase + i));
            byte[] sb = r.sequence().getBytes(StandardCharsets.US_ASCII);
            offsets[i] = seq.size();
            seq.writeBytes(sb);
            qual.writeBytes(r.qualities());
            lengths[i] = sb.length;
            positions[i] = r.position();
            mapqs[i] = (byte) r.mappingQuality();
            flags[i] = r.flags();
            matePos[i] = r.matePosition();
            tlens[i] = r.templateLength();
            cigars.add(r.cigar());
            names.add(r.readName());
            mateChroms.add(r.mateChromosome());
            chroms.add(r.chromosome());
            readTags.add(r.tags());
        }

        String refUri = rd.referenceUri();
        int seqCodec = t.hasCodecs() ? t.codecOf("sequences", b) : 0;
        int qualCodec = t.hasCodecs() ? t.codecOf("qualities", b) : 0;
        Map<String, byte[]> refSeqs = null;
        Map<String, Compression> overrides = new LinkedHashMap<>();
        if (seqCodec == Compression.REF_DIFF_V2.ordinal()) {
            refSeqs = blocksV1ReferenceSeqs(runGroup, study, refUri);
            if (refSeqs == null) {
                throw new IllegalStateException(
                    "per-AU blocks_v1: block " + b + " codes sequences "
                    + "with REF_DIFF_V2 but reference '" + refUri
                    + "' is not embedded in /study/references and not "
                    + "resolvable via REF_PATH; restoring the blob "
                    + "needs the reference bytes");
            }
        } else if (seqCodec != 0) {
            overrides.put("sequences", Compression.values()[seqCodec]);
        }
        if (qualCodec != 0) {
            overrides.put("qualities", Compression.values()[qualCodec]);
        }

        byte[] seqBytes = seq.toByteArray();
        byte[] qualBytes = qual.toByteArray();
        boolean disableV5 =
            readLongRunAttr(runGroup, "opt_disable_qualities_v5") != 0;
        long sliceBytes = readLongRunAttr(runGroup, "ref_diff_slice_bytes");
        WrittenGenomicRun block = new WrittenGenomicRun(
            rd.acquisitionMode(), refUri, rd.platform(), rd.sampleName(),
            positions, mapqs, flags, seqBytes, qualBytes,
            offsets, lengths, cigars, names, mateChroms, matePos, tlens,
            chroms, Compression.ZLIB, overrides,
            List.of(), false, refSeqs, null, null,
            disableV5, false, rd.getReadRole(), sliceBytes,
            tags != null ? tags : readTags);
        return new BlockRun(block, seqBytes, qualBytes);
    }

    private static long readLongRunAttr(StorageGroup g, String name) {
        if (!g.hasAttribute(name)) return 0L;
        Object v = g.getAttribute(name);
        return v instanceof Number n ? n.longValue() : 0L;
    }

    /** The stream writer's sticky qualities discipline: after the
     *  first FQZCOMP_NX16_Z block, read the winning strategy back from
     *  the encoded stream and pin it for the rest of the run. */
    private static int blocksV1DeriveQualHint(
            GenomicBlocks.BlockBlobs blobs, int current) {
        if (current != -1) return current;
        Integer qc = blobs.codecs().get("qualities");
        if (qc == null || qc != Compression.FQZCOMP_NX16_Z.ordinal()) {
            return current;
        }
        int strat = global.thalion.ttio.codecs.FqzcompNx16Z
            .streamStrategy(blobs.blobs().get("qualities"));
        if (strat <= 0) return current;
        return strat == 4
            ? global.thalion.ttio.codecs.FqzcompNx16Z.HINT_V4_AUTO : strat;
    }

    private static void encryptBlocksV1Run(StorageGroup study,
            StorageGroup gRuns, String runName, int datasetId,
            byte[] key) {
        try (StorageGroup runGroup = gRuns.openGroup(runName);
             StorageGroup sig = runGroup.openGroup("signal_channels")) {
            List<String> channels = new ArrayList<>();
            for (String ch : BLOCKS_V1_CHANNELS) {
                if (sig.hasChild(ch)) channels.add(ch);
            }
            if (channels.isEmpty()) return;
            BlockTable t = BlockTable.read(runGroup);
            GenomicRun rd = GenomicRun.readFrom(runGroup, runName);

            Map<String, StorageDataset> segDs = new LinkedHashMap<>();
            long tagOffset = 0;   // global plaintext offset of the tags channel
            try {
                for (String ch : channels) {
                    String segName = ch + "_segments";
                    if (sig.hasChild(segName)) sig.deleteChild(segName);
                    segDs.put(ch, sig.createCompoundDataset(
                        segName, CHANNEL_SEG_FIELDS, 0, true, 1024));
                }
                for (int b = 0; b < t.count(); b++) {
                    BlockRun br = blocksV1BlockRun(rd, runGroup, study,
                                                   t, b, t.readStartAt(b));
                    long[] local = new long[t.nReadsAt(b)];
                    int[] blkLens = br.run().lengths();
                    long cum = 0;
                    for (int i = 0; i < blkLens.length; i++) {
                        local[i] = cum;
                        cum += blkLens[i];
                    }
                    for (String ch : channels) {
                        byte[] plain;
                        long[] chLocal = local;
                        int[] chLens = blkLens;
                        long base = t.baseStartAt(b);
                        if (ch.equals("tags")) {
                            List<String> texts = br.run().tags();
                            int nn = t.nReadsAt(b);
                            chLocal = new long[nn];
                            chLens = new int[nn];
                            java.io.ByteArrayOutputStream tb =
                                new java.io.ByteArrayOutputStream();
                            for (int i = 0; i < nn; i++) {
                                byte[] tx = (texts == null ? "" : texts.get(i))
                                    .getBytes(StandardCharsets.UTF_8);
                                chLocal[i] = tb.size();
                                chLens[i] = tx.length;
                                tb.writeBytes(tx);
                            }
                            plain = tb.toByteArray();
                            base = tagOffset;
                            tagOffset += plain.length;
                        } else {
                            plain = ch.equals("sequences") ? br.seq() : br.qual();
                        }
                        List<ChannelSegment> segs =
                            PerAUEncryption.encryptChannelToSegments(
                                plain, chLocal, chLens, datasetId, ch,
                                key, 1, (int) t.readStartAt(b), base);
                        List<Object[]> rows =
                            new ArrayList<>(segs.size());
                        for (ChannelSegment s : segs) {
                            rows.add(new Object[]{ s.offset(), s.length(),
                                s.iv(), s.tag(), s.ciphertext() });
                        }
                        segDs.get(ch).append(rows);
                    }
                }
            } catch (RuntimeException e) {
                for (String ch : channels) {
                    String segName = ch + "_segments";
                    if (sig.hasChild(segName)) sig.deleteChild(segName);
                }
                throw e;
            } finally {
                for (StorageDataset ds : segDs.values()) ds.close();
            }
            for (String ch : channels) {
                sig.deleteChild(ch);
                sig.setAttribute(ch + "_algorithm", "aes-256-gcm");
            }
        }
    }

    private static void decryptBlocksV1RunInPlace(StorageGroup study,
            StorageGroup gRuns, String runName, int datasetId,
            byte[] key) {
        try (StorageGroup runGroup = gRuns.openGroup(runName);
             StorageGroup sig = runGroup.openGroup("signal_channels")) {
            List<String> channels = new ArrayList<>();
            for (String ch : BLOCKS_V1_CHANNELS) {
                if (sig.hasChild(ch + "_segments")) channels.add(ch);
            }
            if (channels.isEmpty()) return;
            BlockTable t = BlockTable.read(runGroup);
            Map<String, Integer> chromMap = blocksV1ChromMap(sig);
            List<String> chromNames;
            try (StorageGroup idx = runGroup.openGroup("genomic_index")) {
                chromNames = BlockView.readNames(idx, "chromosome_names");
            }
            List<String> mateChromNames = List.of();
            if (sig.hasChild("mate_info")) {
                try (StorageGroup mate = sig.openGroup("mate_info")) {
                    mateChromNames = BlockView.readNames(mate, "chrom_names");
                }
            }
            Set<String> skip = Set.copyOf(channels);

            Map<String, StorageDataset> newDs = new LinkedHashMap<>();
            Map<String, Long> written = new LinkedHashMap<>();
            Map<String, List<long[]>> newRanges = new LinkedHashMap<>();
            for (String ch : channels) {
                written.put(ch, 0L);
                newRanges.put(ch, new ArrayList<>());
            }
            boolean mismatch = false;
            int qualHint = -1;
            byte[] refMd5 = null;
            try {
                for (int b = 0; b < t.count(); b++) {
                    long r0 = t.readStartAt(b);
                    int nn = t.nReadsAt(b);
                    Map<String, byte[]> decrypted = new LinkedHashMap<>();
                    List<String> tagTexts = null;
                    for (String ch : channels) {
                        List<ChannelSegment> segs =
                            readChannelSegmentsSlice(sig, ch + "_segments",
                                                     r0, nn);
                        byte[] flat = PerAUEncryption.decryptChannelFromSegments(
                            segs, datasetId, ch, key, 1, (int) r0);
                        if (ch.equals("tags")) {
                            // M101: one string per read, straight into
                            // the block run (never a codec-0 dataset).
                            tagTexts = splitTagTexts(flat, segs);
                        } else {
                            decrypted.put(ch, flat);
                        }
                    }

                    BlockView.Handle view = BlockView.materialise(
                        runGroup, t, b, chromNames, mateChromNames, skip);
                    BlockRun br;
                    try {
                        StorageGroup viewSig =
                            view.group().openGroup("signal_channels");
                        for (String ch : decrypted.keySet()) {
                            byte[] raw = decrypted.get(ch);
                            try (StorageDataset ds = viewSig.createDataset(
                                    ch, Precision.UINT8, raw.length, 0,
                                    Compression.NONE, 0)) {
                                ds.writeAll(raw);
                                ds.setAttribute("compression", 0);
                            }
                        }
                        GenomicRun rd =
                            GenomicRun.readFrom(view.group(), "block");
                        br = blocksV1BlockRun(rd, runGroup, study, t, b, 0,
                                              tagTexts);
                    } finally {
                        view.discard();
                    }
                    if (refMd5 == null
                            && br.run().referenceChromSeqs() != null) {
                        refMd5 = SpectralDatasetGenomicWriter
                            .referenceMd5ForRun(br.run());
                    }
                    GenomicWriteContext ctx =
                        new GenomicWriteContext(chromMap, refMd5, qualHint);
                    GenomicBlocks.BlockBlobs blobs =
                        GenomicBlocks.encodeBlock(br.run(), ctx);
                    qualHint = blocksV1DeriveQualHint(blobs, qualHint);

                    for (String ch : channels) {
                        byte[] got = blobs.blobs()
                            .getOrDefault(ch, new byte[0]);
                        long off = t.offsetOf(ch, b);
                        long ln = t.lengthOf(ch, b);
                        long pos = written.get(ch);
                        if (got.length != ln || pos != off) {
                            mismatch = true;
                        }
                        newRanges.get(ch).add(new long[]{
                            pos, got.length,
                            blobs.codecs().getOrDefault(ch, 0) });
                        StorageDataset ds = newDs.get(ch);
                        // The channel dataset is created at its first
                        // non-empty blob, with that blob's codec, as the
                        // stream writer does (M101: tags first appearing
                        // after untagged blocks).
                        if (ds == null && got.length == 0) continue;
                        if (ds == null) {
                            StorageGroup parent;
                            String dsName;
                            if (ch.equals("sequences")) {
                                parent = sig.createGroup("sequences");
                                dsName = "data";
                            } else {
                                parent = sig;
                                dsName = ch;
                            }
                            int codec = blobs.codecs().getOrDefault(ch, 0);
                            ds = parent.createDataset(dsName,
                                Precision.UINT8, 0, 256 << 10,
                                codec == 0 ? Compression.ZLIB
                                           : Compression.NONE,
                                6, true);
                            ds.setAttribute("compression", codec);
                            Map<String, Object> extra =
                                blobs.extraAttrs().getOrDefault(ch, Map.of());
                            for (Map.Entry<String, Object> e
                                    : extra.entrySet()) {
                                ds.setAttribute(e.getKey(), e.getValue());
                            }
                            newDs.put(ch, ds);
                        }
                        if (got.length > 0) {
                            ds.append(got);
                        }
                        written.put(ch, pos + got.length);
                    }
                }
            } finally {
                for (StorageDataset ds : newDs.values()) ds.close();
            }
            if (mismatch) {
                rewriteBlocksIndex(runGroup, channels, newRanges);
            }
            for (String ch : channels) {
                sig.deleteChild(ch + "_segments");
                String algAttr = ch + "_algorithm";
                if (sig.hasAttribute(algAttr)) sig.deleteAttribute(algAttr);
            }
        }
    }

    /** Fallback when a re-encoded blob does not land on the ranges the
     *  block index records: patch off/len/codec for the re-encoded
     *  channels and recreate {@code blocks/index}; the other columns
     *  and channels are carried over unchanged. */
    @SuppressWarnings("unchecked")
    private static void rewriteBlocksIndex(StorageGroup runGroup,
            List<String> channels, Map<String, List<long[]>> newRanges) {
        try (StorageGroup blocks = runGroup.openGroup("blocks")) {
            List<Object[]> rows;
            try (StorageDataset ds = blocks.openDataset("index")) {
                rows = (List<Object[]>) ds.readAll();
            }
            // M101: an index with the tags triple keeps it.
            boolean withTags = !rows.isEmpty()
                && rows.get(0).length == GenomicStreamWriter.INDEX_FIELDS_WITH_TAGS.size();
            for (int b = 0; b < rows.size(); b++) {
                Object[] row = rows.get(b);
                for (String ch : channels) {
                    long[] r = newRanges.get(ch).get(b);
                    int oc = GenomicStreamWriter.offsetColumn(ch);
                    row[oc] = r[0];
                    row[oc + 1] = r[1];
                    row[GenomicStreamWriter.codecColumn(ch)] = (int) r[2];
                }
            }
            blocks.deleteChild("index");
            try (StorageDataset out = blocks.createCompoundDataset(
                    "index", GenomicStreamWriter.indexFields(withTags), 0,
                    true, 1024)) {
                if (!rows.isEmpty()) out.append(rows);
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static List<ChannelSegment> readChannelSegmentsSlice(
            StorageGroup parent, String name, long offset, int count) {
        try (StorageDataset ds = parent.openDataset(name)) {
            List<Object[]> rows = (List<Object[]>) ds.readSlice(offset, count);
            List<ChannelSegment> out = new ArrayList<>(rows.size());
            for (Object[] r : rows) {
                out.add(new ChannelSegment(
                    ((Number) r[0]).longValue(),
                    ((Number) r[1]).intValue(),
                    (byte[]) r[2], (byte[]) r[3], (byte[]) r[4]));
            }
            return out;
        }
    }

    // ────────────────────────────────────────────── compound I/O helpers

    private static final List<CompoundField> CHANNEL_SEG_FIELDS = List.of(
        new CompoundField("offset", CompoundField.Kind.INT64),
        new CompoundField("length", CompoundField.Kind.UINT32),
        new CompoundField("iv", CompoundField.Kind.VL_BYTES),
        new CompoundField("tag", CompoundField.Kind.VL_BYTES),
        new CompoundField("ciphertext", CompoundField.Kind.VL_BYTES));

    private static final List<CompoundField> HEADER_SEG_FIELDS = List.of(
        new CompoundField("iv", CompoundField.Kind.VL_BYTES),
        new CompoundField("tag", CompoundField.Kind.VL_BYTES),
        new CompoundField("ciphertext", CompoundField.Kind.VL_BYTES));

    static void writeChannelSegments(StorageGroup parent, String name,
                                      List<ChannelSegment> segments) {
        if (parent.hasChild(name)) parent.deleteChild(name);
        List<Object[]> rows = new ArrayList<>(segments.size());
        for (ChannelSegment s : segments) {
            rows.add(new Object[]{ s.offset(), s.length(),
                                    s.iv(), s.tag(), s.ciphertext() });
        }
        try (StorageDataset ds = parent.createCompoundDataset(name,
                CHANNEL_SEG_FIELDS, rows.size())) {
            ds.writeAll(rows);
        }
    }

    @SuppressWarnings("unchecked")
    static List<ChannelSegment> readChannelSegments(StorageGroup parent,
                                                      String name) {
        try (StorageDataset ds = parent.openDataset(name)) {
            List<Object[]> rows = (List<Object[]>) ds.readAll();
            List<ChannelSegment> out = new ArrayList<>(rows.size());
            for (Object[] r : rows) {
                out.add(new ChannelSegment(
                    ((Number) r[0]).longValue(),
                    ((Number) r[1]).intValue(),
                    (byte[]) r[2], (byte[]) r[3], (byte[]) r[4]));
            }
            return out;
        }
    }

    static void writeHeaderSegments(StorageGroup parent, String name,
                                     List<HeaderSegment> segments) {
        if (parent.hasChild(name)) parent.deleteChild(name);
        List<Object[]> rows = new ArrayList<>(segments.size());
        for (HeaderSegment s : segments) {
            rows.add(new Object[]{ s.iv(), s.tag(), s.ciphertext() });
        }
        try (StorageDataset ds = parent.createCompoundDataset(name,
                HEADER_SEG_FIELDS, rows.size())) {
            ds.writeAll(rows);
        }
    }

    @SuppressWarnings("unchecked")
    static List<HeaderSegment> readHeaderSegments(StorageGroup parent,
                                                    String name) {
        try (StorageDataset ds = parent.openDataset(name)) {
            List<Object[]> rows = (List<Object[]>) ds.readAll();
            List<HeaderSegment> out = new ArrayList<>(rows.size());
            for (Object[] r : rows) {
                out.add(new HeaderSegment((byte[]) r[0], (byte[]) r[1],
                                            (byte[]) r[2]));
            }
            return out;
        }
    }

    // ─────────────────────────────────────────────────────── misc helpers

    private static List<String> runNames(StorageGroup msRuns) {
        List<String> out = new ArrayList<>();
        for (String n : msRuns.childNames()) {
            if (!n.startsWith("_") && msRuns.hasChild(n)) out.add(n);
        }
        return out;
    }

    private static List<String> splitNames(String raw) {
        if (raw == null || raw.isEmpty()) return List.of();
        String[] parts = raw.split(",");
        List<String> out = new ArrayList<>(parts.length);
        for (String p : parts) {
            String t = p.strip();
            if (!t.isEmpty()) out.add(t);
        }
        return out;
    }

    private static Object getAttrOr(StorageGroup g, String name, Object fallback) {
        if (!g.hasAttribute(name)) return fallback;
        Object v = g.getAttribute(name);
        return v == null ? fallback : v;
    }

    private static long[] readLongs(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (long[]) ds.readAll();
        }
    }

    private static int[] readInts(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (int[]) ds.readAll();
        }
    }

    private static double[] readDoubles(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (double[]) ds.readAll();
        }
    }

    private static byte[] doublesToLeBytes(double[] v) {
        ByteBuffer bb = ByteBuffer.allocate(v.length * 8).order(ByteOrder.LITTLE_ENDIAN);
        for (double d : v) bb.putDouble(d);
        return bb.array();
    }

    private static byte[] longsToLeBytes(long[] v) {
        ByteBuffer bb = ByteBuffer.allocate(v.length * 8).order(ByteOrder.LITTLE_ENDIAN);
        for (long l : v) bb.putLong(l);
        return bb.array();
    }

    private static byte[] intsToLeBytes(int[] v) {
        ByteBuffer bb = ByteBuffer.allocate(v.length * 4).order(ByteOrder.LITTLE_ENDIAN);
        for (int i : v) bb.putInt(i);
        return bb.array();
    }

    private static long[] leBytesToLongs(byte[] b) {
        if ((b.length & 7) != 0) {
            throw new IllegalStateException(
                "leBytesToLongs: length " + b.length + " not multiple of 8");
        }
        ByteBuffer bb = ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN);
        long[] out = new long[b.length / 8];
        for (int i = 0; i < out.length; i++) out[i] = bb.getLong();
        return out;
    }

    private static int[] leBytesToInts(byte[] b) {
        if ((b.length & 3) != 0) {
            throw new IllegalStateException(
                "leBytesToInts: length " + b.length + " not multiple of 4");
        }
        ByteBuffer bb = ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN);
        int[] out = new int[b.length / 4];
        for (int i = 0; i < out.length; i++) out[i] = bb.getInt();
        return out;
    }
}
