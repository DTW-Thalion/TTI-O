/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.hdf5.Hdf5Group;
import java.util.*;

/**
 * Reader/writer for TTI-O feature flags stored as root group attributes.
 * {@code @ttio_format_version = "1.1"}
 * {@code @ttio_features = JSON array of feature strings}
 *
 *
 */
public final class FeatureFlags {
    // Required features (must refuse if unrecognized)
    public static final String BASE_V1 = "base_v1";
    public static final String COMPOUND_IDENTIFICATIONS = "compound_identifications";
    public static final String COMPOUND_QUANTIFICATIONS = "compound_quantifications";
    public static final String COMPOUND_PROVENANCE = "compound_provenance";
    public static final String COMPOUND_PER_RUN_PROVENANCE = "compound_per_run_provenance";

    // Optional features (ignorable if unrecognized)
    public static final String OPT_COMPOUND_HEADERS = "opt_compound_headers";
    public static final String OPT_NATIVE_2D_NMR = "opt_native_2d_nmr";
    public static final String OPT_NATIVE_MSIMAGE_CUBE = "opt_native_msimage_cube";
    public static final String OPT_DATASET_ENCRYPTION = "opt_dataset_encryption";
    public static final String OPT_DIGITAL_SIGNATURES = "opt_digital_signatures";
    public static final String OPT_CANONICAL_SIGNATURES = "opt_canonical_signatures";
    public static final String OPT_KEY_ROTATION = "opt_key_rotation";
    public static final String OPT_ANONYMIZED = "opt_anonymized";
    /** file uses post-quantum crypto (ML-KEM-1024 and/or
     *  ML-DSA-87). Opt-flag — a reader without PQC can still open the
     *  file and read unencrypted datasets. */
    public static final String OPT_PQC_PREVIEW = "opt_pqc_preview";
    /** channels encrypted per Access Unit. */
    public static final String OPT_PER_AU_ENCRYPTION = "opt_per_au_encryption";
    /** AU semantic header encrypted. */
    public static final String OPT_ENCRYPTED_AU_HEADERS = "opt_encrypted_au_headers";
    /** per-AU encryption keyed by chromosome (region-based).
     *  Set in addition to {@link #OPT_PER_AU_ENCRYPTION} when a file
     *  was encrypted via the per-chromosome key-map dispatch path.
     *  Clear AUs encode as a {@code <channel>_segments} row with
     *  empty IV/tag and plaintext bytes in the ciphertext slot;
     *  encrypted AUs use the standard 12/16-byte AES-GCM layout.
     * M90.4 */
    public static final String OPT_REGION_KEYED_ENCRYPTION = "opt_region_keyed_encryption";
    /** {@code spectrum_index} carries four optional
     *  parallel columns recording MS/MS activation method and precursor
     *  isolation window. Files that set this flag are written with
     *  format version {@code "1.3"}; pre-M74 readers that ignore the
     *  flag still parse the base columns. */
    public static final String OPT_MS2_ACTIVATION_DETAIL = "opt_ms2_activation_detail";
    /** file contains one or more genomic runs under
     *  {@code /study/genomic_runs/}. Files that set this flag are
     *  written with format version {@code "1.4"}; pre-M82 readers
     *  that ignore the flag still parse the MS pipeline normally
     *  (genomic runs are a separate group hierarchy). M82 */
    public static final String OPT_GENOMIC = "opt_genomic";
    /** v1.2 L4 : genomic runs in this file do NOT carry the
     *  {@code signal_channels/{positions,flags,mapping_qualities}}
     *  duplicates that v1.5 and earlier wrote alongside the canonical
     *  {@code genomic_index/} copies. Tooling that wants to detect
     *  v1.6 layout (without enumerating the {@code signal_channels/}
     *  group) can check for this flag. Strictly additive: pre-v1.6
     *  readers ignore unknown opt_* flags and continue to read from
     *  {@code genomic_index/} correctly. */
    public static final String OPT_NO_SIGNAL_INT_DUPS = "opt_no_signal_int_dups";
    /** file contains one or more assembly graphs under
     *  {@code /study/assembly_graphs/} (M98, format-spec §11a).
     *  Strictly additive: readers that ignore the flag parse the rest
     *  of the file normally (the graphs are a separate group
     *  hierarchy). */
    public static final String OPT_ASSEMBLY_GRAPH = "opt_assembly_graph";
    /** one or more genomic runs carry SAM optional fields as a
     *  {@code signal_channels/tags} channel (SAM_TAGS, codec id 18;
     *  M101, format-spec §10.13). Tag-less files never carry it. */
    public static final String OPT_SAM_TAGS = "opt_sam_tags";
    /** at least one run's {@code spectrum_index/} carries the
     *  {@code pixel_x}, {@code pixel_y}, {@code pixel_z} int32 columns:
     *  each spectrum's position on an imaging grid (M102,
     *  format-spec §4b). Written only when a run carries the columns,
     *  so other files stay byte-identical. */
    public static final String OPT_PIXEL_COORDINATES = "opt_pixel_coordinates";

    private static final Set<String> REQUIRED = Set.of(
        BASE_V1, COMPOUND_IDENTIFICATIONS, COMPOUND_QUANTIFICATIONS,
        COMPOUND_PROVENANCE, COMPOUND_PER_RUN_PROVENANCE
    );

    private static final Set<String> KNOWN_OPTIONAL = Set.of(
        OPT_COMPOUND_HEADERS, OPT_NATIVE_2D_NMR, OPT_NATIVE_MSIMAGE_CUBE,
        OPT_DATASET_ENCRYPTION, OPT_DIGITAL_SIGNATURES, OPT_CANONICAL_SIGNATURES,
        OPT_KEY_ROTATION, OPT_ANONYMIZED, OPT_PQC_PREVIEW,
        OPT_PER_AU_ENCRYPTION, OPT_ENCRYPTED_AU_HEADERS,
        OPT_REGION_KEYED_ENCRYPTION,
        OPT_MS2_ACTIVATION_DETAIL,
        OPT_PIXEL_COORDINATES
    );

    private final String formatVersion;
    private final Set<String> features;

    /**
     * Build a {@code FeatureFlags} from a format version label and a
     * collection of feature strings. The collection is copied into a
     * {@link LinkedHashSet} preserving iteration order.
     *
     * @param formatVersion file format version label (e.g. {@code "1.5"})
     * @param features      collection of feature flag names
     */
    public FeatureFlags(String formatVersion, Collection<String> features) {
        this.formatVersion = formatVersion;
        this.features = new LinkedHashSet<>(features);
    }

    /** @return The format version label (e.g. {@code "1.5"}). */
    public String formatVersion() { return formatVersion; }

    /** @return Unmodifiable view of the feature flag set, iteration order preserved. */
    public Set<String> features() { return Collections.unmodifiableSet(features); }

    /**
     * @param flag feature flag name to test
     * @return     {@code true} when the flag is present in this set
     */
    public boolean has(String flag) { return features.contains(flag); }

    /** Check if this is a v0.1 file (no features attribute). */
    public boolean isV1Legacy() {
        return features.isEmpty() || "1.0.0".equals(formatVersion);
    }

    /** Default features for a new v0.4+ file. */
    public static FeatureFlags defaultCurrent() {
        List<String> flags = new ArrayList<>(REQUIRED);
        return new FeatureFlags("1.1", flags);
    }

    /** Add an optional feature flag. */
    public FeatureFlags with(String flag) {
        Set<String> updated = new LinkedHashSet<>(features);
        updated.add(flag);
        return new FeatureFlags(formatVersion, updated);
    }

    /** Read feature flags from an HDF5 root group. */
    public static FeatureFlags readFrom(Hdf5Group root) {
        String version = "1.0.0";
        Set<String> flags = new LinkedHashSet<>();

        if (root.hasAttribute("ttio_format_version")) {
            version = root.readStringAttribute("ttio_format_version");
        }
        if (root.hasAttribute("ttio_features")) {
            String json = root.readStringAttribute("ttio_features");
            // Parse simple JSON array: ["flag1","flag2",...]
            flags = parseJsonArray(json);
        }
        return new FeatureFlags(version, flags);
    }

    /** Read feature flags from any provider's root group . */
    public static FeatureFlags readFrom(global.thalion.ttio.providers.StorageGroup root) {
        String version = "1.0.0";
        Set<String> flags = new LinkedHashSet<>();
        if (root.hasAttribute("ttio_format_version")) {
            Object v = root.getAttribute("ttio_format_version");
            if (v != null) version = v.toString();
        }
        if (root.hasAttribute("ttio_features")) {
            Object v = root.getAttribute("ttio_features");
            if (v != null) flags = parseJsonArray(v.toString());
        }
        return new FeatureFlags(version, flags);
    }

    /** Write feature flags to an HDF5 root group. */
    public void writeTo(Hdf5Group root) {
        root.setStringAttribute("ttio_format_version", formatVersion);
        root.setStringAttribute("ttio_features", toJsonArray());
    }

    /** Write feature flags to any provider's root group . */
    public void writeTo(global.thalion.ttio.providers.StorageGroup root) {
        root.setAttribute("ttio_format_version", formatVersion);
        root.setAttribute("ttio_features", toJsonArray());
    }

    private String toJsonArray() {
        StringBuilder sb = new StringBuilder("[");
        boolean first = true;
        for (String f : features) {
            if (!first) sb.append(",");
            sb.append("\"").append(f).append("\"");
            first = false;
        }
        sb.append("]");
        return sb.toString();
    }

    static Set<String> parseJsonArray(String json) {
        Set<String> result = new LinkedHashSet<>();
        if (json == null || json.isBlank()) return result;
        // Strip brackets and split by comma
        String inner = json.strip();
        if (inner.startsWith("[")) inner = inner.substring(1);
        if (inner.endsWith("]")) inner = inner.substring(0, inner.length() - 1);
        for (String part : inner.split(",")) {
            String trimmed = part.strip();
            if (trimmed.startsWith("\"") && trimmed.endsWith("\"")) {
                trimmed = trimmed.substring(1, trimmed.length() - 1);
            }
            if (!trimmed.isEmpty()) result.add(trimmed);
        }
        return result;
    }
}
