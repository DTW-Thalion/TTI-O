/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.FeatureFlags;
import global.thalion.ttio.InstrumentConfig;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.SpectrumIndex;
import global.thalion.ttio.exporters.ImzMLWriter;
import global.thalion.ttio.exporters.writers.ImzMLWriterAdapter;
import global.thalion.ttio.importers.ImzMLReader.ImzMLImport;
import global.thalion.ttio.importers.ImzMLReader.PixelSpectrum;
import global.thalion.ttio.providers.MemoryProvider;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * M102 — imzML pixel coordinates stored as {@code spectrum_index/pixel_x/y/z}
 * (format-spec §4b, feature flag {@code opt_pixel_coordinates}).
 */
class ImzMLPixelCoordinatesTest {

    private static final String UUID_HEX = "0123456789abcdef0123456789abcdef";

    // ── helpers ─────────────────────────────────────────────────────

    /** Processed-mode pixels on a {@code w x h} grid (1-based), 3-5
     *  peaks per pixel with pixel-dependent m/z. */
    private static List<PixelSpectrum> processedPixels(int w, int h) {
        List<PixelSpectrum> out = new ArrayList<>(w * h);
        for (int y = 1; y <= h; y++) {
            for (int x = 1; x <= w; x++) {
                int k = 3 + ((x + y) % 3);
                double[] mz = new double[k];
                double[] in = new double[k];
                for (int j = 0; j < k; j++) {
                    mz[j] = 100.0 + 50.0 * j + 0.001 * x + 0.0001 * y;
                    in[j] = 10.0 * (j + 1) + x + 0.5 * y;
                }
                out.add(new PixelSpectrum(x, y, 1, mz, in));
            }
        }
        return out;
    }

    private static SpectrumIndex plainIndex(int n) {
        long[] off = new long[n];
        int[] len = new int[n];
        for (int i = 0; i < n; i++) { off[i] = 2L * i; len[i] = 2; }
        int[] ms = new int[n];
        java.util.Arrays.fill(ms, 1);
        return new SpectrumIndex(n, off, len, new double[n], ms, new int[n],
            new double[n], new int[n], new double[n]);
    }

    // ── SpectrumIndex unit coverage ─────────────────────────────────

    @Test
    void withPixelCoordinatesValidatesAllOrNoneAndLength() {
        SpectrumIndex idx = plainIndex(3);
        assertFalse(idx.hasPixelCoordinates());
        assertNull(idx.pixelX());
        assertNull(idx.pixelCoordinatesAt(0));

        int[] a = {1, 2, 3};
        assertThrows(IllegalArgumentException.class,
            () -> idx.withPixelCoordinates(a, a, null));
        assertThrows(IllegalArgumentException.class,
            () -> idx.withPixelCoordinates(a, new int[]{1, 2}, a));

        SpectrumIndex with = idx.withPixelCoordinates(
            new int[]{1, 2, 3}, new int[]{4, 5, 6}, new int[]{1, 1, 1});
        assertTrue(with.hasPixelCoordinates());
        assertArrayEquals(new int[]{2, 5, 1}, with.pixelCoordinatesAt(1));
        assertFalse(idx.hasPixelCoordinates(), "original index is unchanged");
        assertFalse(with.withPixelCoordinates(null, null, null).hasPixelCoordinates());
    }

    @Test
    void columnsRoundTripAndPartialColumnsAreRejected() {
        String url = "memory://m102-spectrum-index";
        try {
            StorageGroup root = new MemoryProvider()
                .open(url, StorageProvider.Mode.CREATE).rootGroup();
            SpectrumIndex idx = plainIndex(2).withPixelCoordinates(
                new int[]{7, 8}, new int[]{9, 10}, new int[]{1, 2});
            try (StorageGroup run = root.createGroup("run")) {
                idx.writeTo(run);
                SpectrumIndex back = SpectrumIndex.readFrom(run);
                assertArrayEquals(new int[]{7, 8}, back.pixelX());
                assertArrayEquals(new int[]{9, 10}, back.pixelY());
                assertArrayEquals(new int[]{1, 2}, back.pixelZ());

                try (StorageGroup si = run.openGroup("spectrum_index")) {
                    si.deleteChild("pixel_z");
                }
                IllegalStateException e = assertThrows(IllegalStateException.class,
                    () -> SpectrumIndex.readFrom(run));
                assertTrue(e.getMessage().contains("pixel"), e.getMessage());
            }
            try (StorageGroup run2 = root.createGroup("run2")) {
                plainIndex(2).writeTo(run2);
                try (StorageGroup si = run2.openGroup("spectrum_index")) {
                    assertFalse(si.hasChild("pixel_x"));
                }
                assertNull(SpectrumIndex.readFrom(run2).pixelX());
            }
        } finally {
            MemoryProvider.discardStore(url);
        }
    }

    @Test
    void partialColumnsOnDiskFailTheOpen(@TempDir Path tmp) throws Exception {
        Path imzml = tmp.resolve("small.imzML");
        ImzMLWriter.write(processedPixels(3, 2), imzml, null, "processed",
            0, 0, 0, 0, 0, "flyback", UUID_HEX);
        Path tio = tmp.resolve("small.tio");
        ImzMLReader.read(imzml).toImportedDataset(null).write(tio);
        try (StorageProvider p = ProviderRegistry.open(tio.toString(),
                StorageProvider.Mode.READ_WRITE, "hdf5")) {
            try (StorageGroup study = p.rootGroup().openGroup("study");
                 StorageGroup runs = study.openGroup("ms_runs");
                 StorageGroup run = runs.openGroup(ImzMLReader.PIXEL_RUN_NAME);
                 StorageGroup si = run.openGroup("spectrum_index")) {
                si.deleteChild("pixel_y");
            }
        }
        Exception e = assertThrows(Exception.class, () -> {
            try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
                ds.msRuns().get(ImzMLReader.PIXEL_RUN_NAME).spectrumIndex();
            }
        });
        Throwable t = e;
        boolean found = false;
        while (t != null) {
            if (t instanceof IllegalStateException
                    && String.valueOf(t.getMessage()).contains("pixel")) found = true;
            t = t.getCause();
        }
        assertTrue(found, "expected a partial-pixel-column error, got " + e);
    }

    @Test
    void runsWithoutPixelColumnsCarryNoFlag(@TempDir Path tmp) {
        AcquisitionRun run = new AcquisitionRun("run_0001", AcquisitionMode.MS1_DDA,
            plainIndex(2), new InstrumentConfig("", "", "", "", "", ""),
            Map.of("mz", new double[]{1, 2, 3, 4}, "intensity", new double[]{5, 6, 7, 8}),
            List.of(), List.of(), null, 0.0);
        Path tio = tmp.resolve("plain.tio");
        SpectralDataset.create(tio.toString(), "plain", "", List.of(run),
            List.of(), List.of(), List.of()).close();
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            assertFalse(ds.featureFlags().has(FeatureFlags.OPT_PIXEL_COORDINATES));
            AcquisitionRun back = ds.msRuns().get("run_0001");
            assertFalse(back.spectrumIndex().hasPixelCoordinates());
            assertNull(ImzMLReader.pixelCoordinatesOf(back, ds.provenanceRecords()));
        }
    }

    // ── Large processed-mode import: the pre-M102 64 KB failure ─────

    @Test
    void largeProcessedImportStoresCoordinatesAsColumns(@TempDir Path tmp) throws Exception {
        // 15,000 pixels: the pre-M102 CSV (~123 KB) overflowed the 64 KB
        // @provenance_json attribute from ~9,600 pixels on a 100-wide grid.
        int w = 150, h = 100;
        List<PixelSpectrum> src = processedPixels(w, h);
        Path imzml = tmp.resolve("big.imzML");
        Path ibd = tmp.resolve("big.ibd");
        ImzMLWriter.write(src, imzml, ibd, "processed",
            w, h, 1, 25.0, 30.0, "meander", UUID_HEX);

        ImportedDataset imp = ImporterRegistry.specFor("imzml").reader()
            .read(List.of(imzml.toString(), ibd.toString()), Map.of(), null);
        assertNull(imp.image, "processed mode imports as a pixel run, not a cube");
        assertEquals(1, imp.runs.size());
        assertEquals("imzML import: big.imzML", imp.title);
        assertEquals(1, imp.provenance.size());

        Path tio = tmp.resolve("big.tio");
        imp.write(tio);

        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            assertTrue(ds.featureFlags().has(FeatureFlags.OPT_PIXEL_COORDINATES));
            AcquisitionRun run = ds.msRuns().get(ImzMLReader.PIXEL_RUN_NAME);
            assertNotNull(run);
            SpectrumIndex idx = run.spectrumIndex();
            assertEquals(w * h, idx.count());
            assertTrue(idx.hasPixelCoordinates());
            for (int i = 0; i < src.size(); i++) {
                PixelSpectrum p = src.get(i);
                assertEquals(p.x(), idx.pixelX()[i]);
                assertEquals(p.y(), idx.pixelY()[i]);
                assertEquals(p.z(), idx.pixelZ()[i]);
                assertEquals(1, idx.msLevelAt(i));
                assertEquals(0.0, idx.retentionTimeAt(i));
            }
            assertEquals(src.get(5).intensity()[src.get(5).intensity().length - 1],
                idx.basePeakIntensityAt(5));

            List<ProvenanceRecord> all = new ArrayList<>(ds.provenanceRecords());
            all.addAll(run.provenanceRecords());
            assertFalse(all.isEmpty());
            for (ProvenanceRecord r : all) {
                assertFalse(r.parameters().containsKey(
                    ImzMLReader.PARAM_LEGACY_COORDINATES_CSV),
                    "coordinates must not be written to provenance");
            }
            ProvenanceRecord runProv = run.provenanceRecords().get(0);
            assertEquals(ImzMLReader.IMPORTER_SOFTWARE, runProv.software());
            assertEquals("processed", runProv.parameters().get(ImzMLReader.PARAM_MODE));
            assertEquals(UUID_HEX, runProv.parameters().get(ImzMLReader.PARAM_UUID_HEX));
            assertEquals("150", runProv.parameters().get(ImzMLReader.PARAM_GRID_MAX_X));
            assertEquals(List.of(imzml.toString(), ibd.toString()), runProv.inputRefs());

            int[][] coords = ImzMLReader.pixelCoordinatesOf(run, ds.provenanceRecords());
            assertNotNull(coords);
            assertArrayEquals(new int[]{src.get(12345).x(), src.get(12345).y(), 1}, coords[12345]);
        }

        // Export through the registry writer adapter and read it back.
        Path out = tmp.resolve("exported.imzML");
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            new ImzMLWriterAdapter().write(ds, null, out, Map.of());
        }
        ImzMLImport back = ImzMLReader.read(out, tmp.resolve("exported.ibd"));
        assertEquals("processed", back.mode());
        assertEquals(UUID_HEX, back.uuidHex());
        assertEquals(w, back.gridMaxX());
        assertEquals(h, back.gridMaxY());
        assertEquals(25.0, back.pixelSizeX());
        assertEquals(30.0, back.pixelSizeY());
        assertEquals(src.size(), back.spectra().size());
        for (int i = 0; i < src.size(); i++) {
            PixelSpectrum a = src.get(i), b = back.spectra().get(i);
            assertEquals(a.x(), b.x());
            assertEquals(a.y(), b.y());
            assertEquals(a.z(), b.z());
            assertArrayEquals(a.mz(), b.mz());
            assertArrayEquals(a.intensity(), b.intensity());
        }
    }

    // ── Legacy: coordinates in the pre-M102 provenance parameter ────

    @Test
    void legacyCsvParameterIsReadAndExported(@TempDir Path tmp) throws Exception {
        int n = 4;
        int[][] coords = {{1, 1, 1}, {2, 1, 1}, {1, 2, 1}, {2, 2, 1}};
        double[] mzAxis = {150.0, 250.0, 350.0};
        double[] mz = new double[n * 3];
        double[] in = new double[n * 3];
        long[] off = new long[n];
        int[] len = new int[n];
        int[] ms = new int[n];
        double[] bp = new double[n];
        StringBuilder csv = new StringBuilder();
        for (int i = 0; i < n; i++) {
            off[i] = 3L * i;
            len[i] = 3;
            ms[i] = 1;
            for (int j = 0; j < 3; j++) {
                mz[3 * i + j] = mzAxis[j];
                in[3 * i + j] = 100.0 * i + j + 1;
            }
            bp[i] = 100.0 * i + 3;
            if (i > 0) csv.append(';');
            csv.append(coords[i][0]).append(',').append(coords[i][1])
               .append(',').append(coords[i][2]);
        }
        Map<String, String> params = new LinkedHashMap<>();
        params.put(ImzMLReader.PARAM_MODE, "continuous");
        params.put(ImzMLReader.PARAM_UUID_HEX, UUID_HEX);
        params.put(ImzMLReader.PARAM_GRID_MAX_X, "2");
        params.put(ImzMLReader.PARAM_GRID_MAX_Y, "2");
        params.put(ImzMLReader.PARAM_GRID_MAX_Z, "1");
        params.put(ImzMLReader.PARAM_PIXEL_SIZE_X, "12.5");
        params.put(ImzMLReader.PARAM_PIXEL_SIZE_Y, "12.5");
        params.put(ImzMLReader.PARAM_SCAN_PATTERN, "flyback");
        params.put(ImzMLReader.PARAM_LEGACY_COORDINATES_CSV, csv.toString());
        ProvenanceRecord prov = new ProvenanceRecord(1700000000L,
            ImzMLReader.IMPORTER_SOFTWARE, params,
            List.of("legacy.imzML", "legacy.ibd"), List.of());

        Map<String, double[]> ch = new LinkedHashMap<>();
        ch.put("mz", mz);
        ch.put("intensity", in);
        AcquisitionRun run = new AcquisitionRun(ImzMLReader.PIXEL_RUN_NAME,
            AcquisitionMode.MS1_DDA,
            new SpectrumIndex(n, off, len, new double[n], ms, new int[n],
                new double[n], new int[n], bp),
            new InstrumentConfig("", "", "", "", "", ""),
            ch, List.of(), List.of(prov), null, 0.0);

        Path tio = tmp.resolve("legacy.tio");
        SpectralDataset.create(tio.toString(), "imzML import: legacy.imzML", "",
            List.of(run), List.of(), List.of(), List.of(prov)).close();

        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            assertFalse(ds.featureFlags().has(FeatureFlags.OPT_PIXEL_COORDINATES),
                "no columns, no flag");
            AcquisitionRun back = ds.msRuns().get(ImzMLReader.PIXEL_RUN_NAME);
            assertFalse(back.spectrumIndex().hasPixelCoordinates());
            int[][] got = ImzMLReader.pixelCoordinatesOf(back, ds.provenanceRecords());
            assertNotNull(got);
            for (int i = 0; i < n; i++) assertArrayEquals(coords[i], got[i]);
            // Dataset-level fallback when the run record lacks the CSV.
            AcquisitionRun noRunProv = new AcquisitionRun("r", AcquisitionMode.MS1_DDA,
                back.spectrumIndex(), null, Map.of(), List.of(), List.of(), null, 0.0);
            int[][] fromDs = ImzMLReader.pixelCoordinatesOf(noRunProv, ds.provenanceRecords());
            assertNotNull(fromDs);
            assertArrayEquals(coords[3], fromDs[3]);
        }

        Path out = tmp.resolve("legacy_out.imzML");
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            new ImzMLWriterAdapter().write(ds, null, out, Map.of());
        }
        ImzMLImport exported = ImzMLReader.read(out, tmp.resolve("legacy_out.ibd"));
        assertEquals("continuous", exported.mode());
        assertEquals(UUID_HEX, exported.uuidHex());
        assertEquals(12.5, exported.pixelSizeX());
        assertEquals(n, exported.spectra().size());
        for (int i = 0; i < n; i++) {
            PixelSpectrum p = exported.spectra().get(i);
            assertEquals(coords[i][0], p.x());
            assertEquals(coords[i][1], p.y());
            assertEquals(coords[i][2], p.z());
            assertArrayEquals(mzAxis, p.mz());
            assertArrayEquals(java.util.Arrays.copyOfRange(in, 3 * i, 3 * i + 3),
                p.intensity());
        }
    }

    @Test
    void legacyCsvWithWrongTripleCountIsIgnored() {
        int[][] ok = ImzMLReader.parseLegacyCsv("1,1,1;2,1,1", 2);
        assertNotNull(ok);
        assertNull(ImzMLReader.parseLegacyCsv("1,1,1;2,1,1", 3));
        assertNull(ImzMLReader.parseLegacyCsv("1,1;2,1,1", 2));
        assertNull(ImzMLReader.parseLegacyCsv("1,a,1;2,1,1", 2));
    }

    @Test
    void datasetWithNeitherImageNorPixelRunStillFailsExport(@TempDir Path tmp) {
        AcquisitionRun run = new AcquisitionRun("run_0001", AcquisitionMode.MS1_DDA,
            plainIndex(2), new InstrumentConfig("", "", "", "", "", ""),
            Map.of("mz", new double[]{1, 2, 3, 4}, "intensity", new double[]{5, 6, 7, 8}),
            List.of(), List.of(), null, 0.0);
        Path tio = tmp.resolve("noimg.tio");
        SpectralDataset.create(tio.toString(), "plain", "", List.of(run),
            List.of(), List.of(), List.of()).close();
        try (SpectralDataset ds = SpectralDataset.open(tio.toString())) {
            assertThrows(IllegalArgumentException.class, () ->
                new ImzMLWriterAdapter().write(ds, null, tmp.resolve("x.imzML"), Map.of()));
        }
    }
}
