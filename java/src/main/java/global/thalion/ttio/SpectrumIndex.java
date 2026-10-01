/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.Enums.ActivationMethod;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.Enums.Polarity;
import global.thalion.ttio.Enums.Precision;
import global.thalion.ttio.genomics.GenomicIndex;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

/**
 * Compressed-domain query index for spectra in an acquisition run.
 *
 * <p>Parallel arrays: offsets, lengths, retention times, MS levels,
 * polarities, precursor m/z, precursor charges, base peak intensities.
 * Kept entirely in memory; signal channels remain lazy on disk.</p>
 *
 * <p>HDF5 layout: {@code <run>/spectrum_index/} group with named datasets.</p>
 *
 * <p><b>API status:</b> Stable.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOSpectrumIndex}, Python
 * {@code ttio.acquisition_run.SpectrumIndex}.</p>
 *
 *
 */
public class SpectrumIndex {

    private final int count;
    private final long[] offsets;
    private final int[] lengths;
    private final double[] retentionTimes;
    private final int[] msLevels;
    private final int[] polarities;
    private final double[] precursorMzs;
    private final int[] precursorCharges;
    private final double[] basePeakIntensities;
    // M74 : optional parallel columns. Non-null iff the file was
    // written with the opt_ms2_activation_detail feature flag set.
    // All four are null-or-all-populated (enforced by constructors).
    private final int[] activationMethods;
    private final double[] isolationTargetMzs;
    private final double[] isolationLowerOffsets;
    private final double[] isolationUpperOffsets;
    // Optional per-spectrum centroided flag. {@code null} for files
    // written before the centroided column was added; otherwise
    // {@code 0 = profile, 1 = centroided} per spectrum. Mirrors mzML
    // CV terms MS:1000127 (centroid) and MS:1000128 (profile).
    private final int[] centroideds;
    // M102: optional imaging pixel coordinates (format-spec §4b).
    // All three are null or all three are populated with length
    // {@code count}; set only through {@link #withPixelCoordinates}.
    private int[] pixelX;
    private int[] pixelY;
    private int[] pixelZ;

    /** Pre-M74 legacy constructor; defaults the four M74 columns to null. */
    public SpectrumIndex(int count, long[] offsets, int[] lengths,
                         double[] retentionTimes, int[] msLevels, int[] polarities,
                         double[] precursorMzs, int[] precursorCharges,
                         double[] basePeakIntensities) {
        this(count, offsets, lengths, retentionTimes, msLevels, polarities,
             precursorMzs, precursorCharges, basePeakIntensities,
             null, null, null, null);
    }

    /** M74 constructor (pre-centroided overload retained for source-compat). */
    public SpectrumIndex(int count, long[] offsets, int[] lengths,
                         double[] retentionTimes, int[] msLevels, int[] polarities,
                         double[] precursorMzs, int[] precursorCharges,
                         double[] basePeakIntensities,
                         int[] activationMethods,
                         double[] isolationTargetMzs,
                         double[] isolationLowerOffsets,
                         double[] isolationUpperOffsets) {
        this(count, offsets, lengths, retentionTimes, msLevels, polarities,
             precursorMzs, precursorCharges, basePeakIntensities,
             activationMethods, isolationTargetMzs,
             isolationLowerOffsets, isolationUpperOffsets, null);
    }

    /** Full constructor including the optional {@code centroideds} column.
     *
     *  <p>The four M74 activation / isolation arrays must be either all-null
     *  or all non-null with length equal to {@code count}. The
     *  {@code centroideds} array is independently optional: {@code null}
     *  for legacy files; otherwise length must equal {@code count}.</p>
     */
    public SpectrumIndex(int count, long[] offsets, int[] lengths,
                         double[] retentionTimes, int[] msLevels, int[] polarities,
                         double[] precursorMzs, int[] precursorCharges,
                         double[] basePeakIntensities,
                         int[] activationMethods,
                         double[] isolationTargetMzs,
                         double[] isolationLowerOffsets,
                         double[] isolationUpperOffsets,
                         int[] centroideds) {
        boolean anyNull = activationMethods == null
                || isolationTargetMzs == null
                || isolationLowerOffsets == null
                || isolationUpperOffsets == null;
        boolean allNull = activationMethods == null
                && isolationTargetMzs == null
                && isolationLowerOffsets == null
                && isolationUpperOffsets == null;
        if (anyNull && !allNull) {
            throw new IllegalArgumentException(
                "M74 columns must be all-null or all-populated");
        }
        if (centroideds != null && centroideds.length != count) {
            throw new IllegalArgumentException(
                "centroideds length " + centroideds.length
                + " does not match spectrum count " + count);
        }
        this.count = count;
        this.offsets = offsets;
        this.lengths = lengths;
        this.retentionTimes = retentionTimes;
        this.msLevels = msLevels;
        this.polarities = polarities;
        this.precursorMzs = precursorMzs;
        this.precursorCharges = precursorCharges;
        this.basePeakIntensities = basePeakIntensities;
        this.activationMethods = activationMethods;
        this.isolationTargetMzs = isolationTargetMzs;
        this.isolationLowerOffsets = isolationLowerOffsets;
        this.isolationUpperOffsets = isolationUpperOffsets;
        this.centroideds = centroideds;
    }

    /** @return Number of spectra indexed; equals the length of every parallel column. */
    public int count() { return count; }

    /** @return Per-spectrum element offsets into the run's concatenated channel buffers. */
    public long[] offsets() { return offsets; }

    /** @return Per-spectrum element counts (peaks) parallel to {@link #offsets()}. */
    public int[] lengths() { return lengths; }

    /** @return Per-spectrum retention times in seconds. */
    public double[] retentionTimes() { return retentionTimes; }

    /** @return Per-spectrum MS levels (1 for MS1, 2 for MS2, ...). */
    public int[] msLevels() { return msLevels; }

    /** @return Per-spectrum polarity codes; decode via {@link Polarity#fromInt(int)}. */
    public int[] polarities() { return polarities; }

    /** @return Per-spectrum precursor m/z; 0.0 for MS1 entries. */
    public double[] precursorMzs() { return precursorMzs; }

    /** @return Per-spectrum precursor charge states; 0 for MS1 entries. */
    public int[] precursorCharges() { return precursorCharges; }

    /** @return Per-spectrum base-peak intensities. */
    public double[] basePeakIntensities() { return basePeakIntensities; }

    /** @param i spectrum index. @return element offset of spectrum {@code i} into the channel buffers. */
    public long offsetAt(int i) { return offsets[i]; }

    /** @param i spectrum index. @return element count (peaks) of spectrum {@code i}. */
    public int lengthAt(int i) { return lengths[i]; }

    /** @param i spectrum index. @return retention time of spectrum {@code i} in seconds. */
    public double retentionTimeAt(int i) { return retentionTimes[i]; }

    /** @param i spectrum index. @return MS level of spectrum {@code i}. */
    public int msLevelAt(int i) { return msLevels[i]; }

    /** @param i spectrum index. @return polarity enum value of spectrum {@code i}. */
    public Polarity polarityAt(int i) { return Polarity.fromInt(polarities[i]); }

    /** @param i spectrum index. @return precursor m/z of spectrum {@code i} (0.0 for MS1). */
    public double precursorMzAt(int i) { return precursorMzs[i]; }

    /** @param i spectrum index. @return precursor charge state of spectrum {@code i} (0 for MS1). */
    public int precursorChargeAt(int i) { return precursorCharges[i]; }

    /** @param i spectrum index. @return base-peak intensity of spectrum {@code i}. */
    public double basePeakIntensityAt(int i) { return basePeakIntensities[i]; }

    /** Optional parallel columns; {@code null} when the file was
     *  written without {@code opt_ms2_activation_detail}.
     *  @return per-spectrum activation method codes, or {@code null} */
    public int[] activationMethods() { return activationMethods; }

    /** @return per-spectrum precursor isolation target m/z, or {@code null} when absent. */
    public double[] isolationTargetMzs() { return isolationTargetMzs; }

    /** @return per-spectrum lower-side isolation offset in m/z, or {@code null} when absent. */
    public double[] isolationLowerOffsets() { return isolationLowerOffsets; }

    /** @return per-spectrum upper-side isolation offset in m/z, or {@code null} when absent. */
    public double[] isolationUpperOffsets() { return isolationUpperOffsets; }

    /** Optional per-spectrum centroided column; {@code null} when the
     *  file was written before centroided tracking. */
    public int[] centroideds() { return centroideds; }

    /** Returns whether spectrum {@code i} is centroided (mzML
     *  MS:1000127). Returns {@code false} when the column is absent —
     *  callers that need to distinguish "unknown" from "profile"
     *  should check {@link #centroideds()} for {@code null}.
     *
     *  <p><b>Cross-language equivalents:</b> Python
     *  {@code SpectrumIndex.centroided_at(i)}, Objective-C
     *  {@code -[TTIOSpectrumIndex centroidedAt:]}.</p>
     */
    public boolean centroidedAt(int i) {
        return centroideds != null && centroideds[i] != 0;
    }

    /** (M74) Returns the activation method at spectrum {@code i};
     *  {@link ActivationMethod#NONE} when the M74 column is absent or
     *  the stored value is 0. */
    public ActivationMethod activationMethodAt(int i) {
        if (activationMethods == null) return ActivationMethod.NONE;
        return ActivationMethod.fromInt(activationMethods[i]);
    }

    /** (M74) Returns the isolation window at spectrum {@code i}, or
     *  {@code null} when the M74 columns are absent or the stored
     *  target+offsets are all zero (MS1 sentinel). */
    public IsolationWindow isolationWindowAt(int i) {
        if (isolationTargetMzs == null || isolationLowerOffsets == null
                || isolationUpperOffsets == null) return null;
        double t = isolationTargetMzs[i];
        double l = isolationLowerOffsets[i];
        double u = isolationUpperOffsets[i];
        if (t == 0.0 && l == 0.0 && u == 0.0) return null;
        return new IsolationWindow(t, l, u);
    }

    /** (M102) @return per-spectrum pixel x positions, or {@code null}
     *  when the run carries no pixel coordinates. */
    public int[] pixelX() { return pixelX; }

    /** (M102) @return per-spectrum pixel y positions, or {@code null}. */
    public int[] pixelY() { return pixelY; }

    /** (M102) @return per-spectrum pixel z positions ({@code 1} for a
     *  2-D acquisition), or {@code null}. */
    public int[] pixelZ() { return pixelZ; }

    /** (M102) @return {@code true} when the three
     *  {@code pixel_x/pixel_y/pixel_z} columns are present. */
    public boolean hasPixelCoordinates() { return pixelX != null; }

    /** (M102) @param i spectrum index.
     *  @return {@code {x, y, z}} of spectrum {@code i}, or {@code null}
     *  when the run carries no pixel coordinates. */
    public int[] pixelCoordinatesAt(int i) {
        if (pixelX == null) return null;
        return new int[]{ pixelX[i], pixelY[i], pixelZ[i] };
    }

    /** (M102) Returns a copy of this index carrying the given pixel
     *  coordinate columns (format-spec §4b). The three arrays must be
     *  all {@code null} (the copy carries no coordinates) or all
     *  non-null with length {@link #count()}; the other columns are
     *  shared with this index.
     *
     *  <p><b>Cross-language equivalents:</b> Python
     *  {@code SpectrumIndex.pixel_x/pixel_y/pixel_z}, Objective-C
     *  {@code TTIOSpectrumIndex pixelX/pixelY/pixelZ}.</p>
     *
     *  @throws IllegalArgumentException on partial or mis-sized columns */
    public SpectrumIndex withPixelCoordinates(int[] x, int[] y, int[] z) {
        boolean anyNull = x == null || y == null || z == null;
        boolean allNull = x == null && y == null && z == null;
        if (anyNull && !allNull) {
            throw new IllegalArgumentException(
                "pixel coordinate columns must be all-null or all-populated");
        }
        if (!allNull && (x.length != count || y.length != count || z.length != count)) {
            throw new IllegalArgumentException(
                "pixel coordinate column lengths (" + x.length + ", " + y.length
                + ", " + z.length + ") do not match spectrum count " + count);
        }
        SpectrumIndex copy = new SpectrumIndex(count, offsets, lengths,
                retentionTimes, msLevels, polarities, precursorMzs,
                precursorCharges, basePeakIntensities,
                activationMethods, isolationTargetMzs,
                isolationLowerOffsets, isolationUpperOffsets, centroideds);
        copy.pixelX = x;
        copy.pixelY = y;
        copy.pixelZ = z;
        return copy;
    }

    /**
     * @return indices whose retention time lies within
     *         {@code [range.minimum(), range.maximum()]}.
     */
    public java.util.List<Integer> indicesInRetentionTimeRange(ValueRange range) {
        java.util.List<Integer> out = new java.util.ArrayList<>();
        for (int i = 0; i < count; i++) {
            double t = retentionTimes[i];
            if (t >= range.minimum() && t <= range.maximum()) out.add(i);
        }
        return out;
    }

    /** @return indices whose {@code msLevel} equals {@code msLevel}. */
    public java.util.List<Integer> indicesForMsLevel(int msLevel) {
        java.util.List<Integer> out = new java.util.ArrayList<>();
        for (int i = 0; i < count; i++) if (msLevels[i] == msLevel) out.add(i);
        return out;
    }

    /** Write this index to a storage group (creates spectrum_index/ subgroup).
     *
     *  <p>parameter type relaxed to {@link StorageGroup} so the
     *  index can be written through any provider (HDF5, SQLite, Memory).</p>
     *
     *  <p>The mathematically redundant {@code offsets} column is omitted;
     *  readers synthesize it from {@code cumsum(lengths)}.</p> */
    public void writeTo(StorageGroup runGroup) {
        try (StorageGroup idx = runGroup.createGroup("spectrum_index")) {
            idx.setAttribute("count", (long) count);

            writeDataset(idx, "lengths", Precision.UINT32, lengths);
            writeDataset(idx, "retention_times", Precision.FLOAT64, retentionTimes);
            writeDataset(idx, "ms_levels", Precision.INT32, msLevels);
            writeDataset(idx, "polarities", Precision.INT32, polarities);
            writeDataset(idx, "precursor_mzs", Precision.FLOAT64, precursorMzs);
            writeDataset(idx, "precursor_charges", Precision.INT32, precursorCharges);
            writeDataset(idx, "base_peak_intensities", Precision.FLOAT64, basePeakIntensities);
            // M74 schema-gating: emit the four optional columns only
            // when they were supplied. Constructor already enforces
            // all-or-nothing, so checking one covers all four.
            if (activationMethods != null) {
                writeDataset(idx, "activation_methods", Precision.INT32, activationMethods);
                writeDataset(idx, "isolation_target_mzs", Precision.FLOAT64, isolationTargetMzs);
                writeDataset(idx, "isolation_lower_offsets", Precision.FLOAT64, isolationLowerOffsets);
                writeDataset(idx, "isolation_upper_offsets", Precision.FLOAT64, isolationUpperOffsets);
            }
            // Optional centroided column (independent of M74 gating).
            if (centroideds != null) {
                writeDataset(idx, "centroideds", Precision.INT32, centroideds);
            }
            // M102: pixel coordinates, all three or none (format-spec §4b).
            if (pixelX != null) {
                writeDataset(idx, "pixel_x", Precision.INT32, pixelX);
                writeDataset(idx, "pixel_y", Precision.INT32, pixelY);
                writeDataset(idx, "pixel_z", Precision.INT32, pixelZ);
            }
        }
    }

    /** Read spectrum index from an existing run group.
     *
     *  <p>parameter type relaxed to {@link StorageGroup}.</p> */
    public static SpectrumIndex readFrom(StorageGroup runGroup) {
        try (StorageGroup idx = runGroup.openGroup("spectrum_index")) {
            int count = ((Number) idx.getAttribute("count")).intValue();

            // offsets is omitted from disk by default and
            // computed from cumsum(lengths) at read time. Pre-v1.10
            // files have it on disk (read directly).
            int[] lengths = readInts(idx, "lengths");
            long[] offsets = idx.hasChild("offsets")
                ? readLongs(idx, "offsets")
                : GenomicIndex.offsetsFromLengths(lengths);
            double[] retentionTimes = readDoubles(idx, "retention_times");
            int[] msLevels = readInts(idx, "ms_levels");
            int[] polarities = readInts(idx, "polarities");
            double[] precursorMzs = readDoubles(idx, "precursor_mzs");
            int[] precursorCharges = readInts(idx, "precursor_charges");
            double[] basePeakIntensities = readDoubles(idx, "base_peak_intensities");

            // M74 schema-gating: probe for the four optional columns.
            // Present-all or absent-all is the contract; partial
            // presence indicates a malformed file and is flagged.
            boolean hasAct = idx.hasChild("activation_methods");
            boolean hasTgt = idx.hasChild("isolation_target_mzs");
            boolean hasLo = idx.hasChild("isolation_lower_offsets");
            boolean hasHi = idx.hasChild("isolation_upper_offsets");
            if (hasAct != hasTgt || hasAct != hasLo || hasAct != hasHi) {
                throw new IllegalStateException(
                    "spectrum_index is malformed: partial M74 columns present");
            }
            int[] activationMethods = hasAct ? readInts(idx, "activation_methods") : null;
            double[] isolationTargetMzs = hasAct ? readDoubles(idx, "isolation_target_mzs") : null;
            double[] isolationLowerOffsets = hasAct ? readDoubles(idx, "isolation_lower_offsets") : null;
            double[] isolationUpperOffsets = hasAct ? readDoubles(idx, "isolation_upper_offsets") : null;

            // Optional centroided column — independent of M74 gating.
            int[] centroideds = idx.hasChild("centroideds")
                    ? readInts(idx, "centroideds") : null;

            SpectrumIndex out = new SpectrumIndex(count, offsets, lengths, retentionTimes,
                    msLevels, polarities, precursorMzs, precursorCharges,
                    basePeakIntensities,
                    activationMethods, isolationTargetMzs,
                    isolationLowerOffsets, isolationUpperOffsets,
                    centroideds);

            // M102: optional pixel coordinates; partial presence is a
            // malformed file, like the M74 check above.
            boolean hasPx = idx.hasChild("pixel_x");
            boolean hasPy = idx.hasChild("pixel_y");
            boolean hasPz = idx.hasChild("pixel_z");
            if (hasPx != hasPy || hasPx != hasPz) {
                throw new IllegalStateException(
                    "spectrum_index is malformed: partial pixel coordinate columns present");
            }
            if (hasPx) {
                int[] px = readInts(idx, "pixel_x");
                int[] py = readInts(idx, "pixel_y");
                int[] pz = readInts(idx, "pixel_z");
                if (px.length != count || py.length != count || pz.length != count) {
                    throw new IllegalStateException(
                        "spectrum_index is malformed: pixel coordinate column length "
                        + "does not match count " + count);
                }
                out = out.withPixelCoordinates(px, py, pz);
            }
            return out;
        }
    }

    // Format parity: Python and ObjC chunk + zlib-compress the 8
    // parallel index datasets (offsets, lengths, retention_times,
    // ms_levels, polarities, precursor_mzs, precursor_charges,
    // base_peak_intensities). Previously Java wrote them contiguous +
    // uncompressed (chunkSize=0, Compression.NONE) which saved a
    // zlib pass but inflated files by ~4.8 MB at 100 K spectra and
    // broke bit-level parity with the other two writers.
    static final int INDEX_CHUNK_SIZE = 4096;

    private static void writeDataset(StorageGroup group, String name,
                                     Precision precision, Object data) {
        int len = java.lang.reflect.Array.getLength(data);
        // ZarrProvider (Java v0.8) throws UnsupportedOperationException
        // on compressed datasets — probe and fall back to contiguous
        // uncompressed for providers that don't implement zlib. HDF5
        // and Memory providers take the compressed path.
        StorageDataset ds;
        try {
            ds = group.createDataset(name, precision, len,
                    INDEX_CHUNK_SIZE, Compression.ZLIB, 6);
        } catch (UnsupportedOperationException e) {
            ds = group.createDataset(name, precision, len,
                    0, Compression.NONE, 0);
        }
        try (StorageDataset closeMe = ds) {
            closeMe.writeAll(data);
        }
    }

    private static double[] readDoubles(StorageGroup group, String name) {
        try (StorageDataset ds = group.openDataset(name)) {
            return (double[]) ds.readAll();
        }
    }

    private static int[] readInts(StorageGroup group, String name) {
        try (StorageDataset ds = group.openDataset(name)) {
            return (int[]) ds.readAll();
        }
    }

    private static long[] readLongs(StorageGroup group, String name) {
        try (StorageDataset ds = group.openDataset(name)) {
            return (long[]) ds.readAll();
        }
    }
}
