/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.exporters;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.SpectrumIndex;
import global.thalion.ttio.importers.ImzMLReader;
import global.thalion.ttio.importers.ImzMLReader.ImzMLImport;
import global.thalion.ttio.importers.ImzMLReader.PixelSpectrum;
import global.thalion.ttio.io.ProgressSink;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.UUID;

/**
 * imzML + .ibd exporter — v0.9+.
 *
 * <p>Reverses {@link ImzMLReader}: takes a list of pixel spectra plus
 * grid metadata and emits the canonical paired {@code .imzML} +
 * {@code .ibd} files. Supports both continuous mode (shared m/z axis
 * once at the head of the .ibd + per-pixel intensity arrays) and
 * processed mode (per-pixel m/z + intensity arrays).</p>
 *
 * <p>The emitted XML uses the canonical IMS accessions that real-world
 * imzML files (pyimzML test corpus, imzML 1.1 spec) use:
 * {@code IMS:1000080} for UUID, {@code IMS:1000042 / 1000043} for max
 * pixel counts. Output round-trips through our {@link ImzMLReader}
 * bit-identically and is readable by pyimzml.</p>
 *
 * <p><b>Cross-language equivalents:</b><br>
 * Python: {@code ttio.exporters.imzml} &middot;
 * Objective-C: {@code TTIOImzMLWriter}</p>
 *
 *
 */
public final class ImzMLWriter {

    /** Return value of {@link #write}: resolved paths + UUID. */
    public record WriteResult(
        Path imzmlPath,
        Path ibdPath,
        String uuidHex,
        String mode,
        int nPixels
    ) {}

    /** Stage D: emit-every-N cadence for {@link ProgressSink} callbacks
     *  during imzML pixel serialisation. Mirrors
     *  {@code ImzMLReader.PROGRESS_INTERVAL_SPECTRA = 100}. */
    public static final int PROGRESS_INTERVAL_SPECTRA = 100;

    private ImzMLWriter() {}

    /**
     * Write {@code pixels} as an imzML + .ibd pair.
     *
     * @param pixels       List of {@link PixelSpectrum}.
     * @param imzmlPath    Destination .imzML path.
     * @param ibdPath      Optional .ibd path; null => derive by swapping extension.
     * @param mode         {@code "continuous"} or {@code "processed"}.
     * @param gridMaxX/Y/Z Pixel grid extents (0 => derive from pixel coordinates).
     * @param pixelSizeX/Y Pixel size in micrometres; 0 => omit cvParam.
     * @param scanPattern  Free-text scan pattern ("flyback", etc.).
     * @param uuidHex      Optional UUID string (dashes/braces tolerated); null => random UUID4.
     */
    public static WriteResult write(
            List<PixelSpectrum> pixels,
            Path imzmlPath,
            Path ibdPath,
            String mode,
            int gridMaxX, int gridMaxY, int gridMaxZ,
            double pixelSizeX, double pixelSizeY,
            String scanPattern,
            String uuidHex
    ) {
        return write(pixels, imzmlPath, ibdPath, mode,
            gridMaxX, gridMaxY, gridMaxZ,
            pixelSizeX, pixelSizeY,
            scanPattern, uuidHex, ProgressSink.discard());
    }

    /**
     * Stage D overload of the main {@link #write} entry that fires
     * {@code progress.onProgress(pixelsDone, totalPixels)} every
     * {@link #PROGRESS_INTERVAL_SPECTRA} pixels during the .ibd
     * encode loop and a final fire once both .imzML + .ibd are
     * flushed to disk.
     *
     * @since 1.5.0
     */
    public static WriteResult write(
            List<PixelSpectrum> pixels,
            Path imzmlPath,
            Path ibdPath,
            String mode,
            int gridMaxX, int gridMaxY, int gridMaxZ,
            double pixelSizeX, double pixelSizeY,
            String scanPattern,
            String uuidHex,
            ProgressSink progress
    ) {
        if (progress == null) progress = ProgressSink.discard();
        if (!"continuous".equals(mode) && !"processed".equals(mode)) {
            throw new IllegalArgumentException(
                    "mode must be 'continuous' or 'processed', got: " + mode);
        }
        if (pixels == null || pixels.isEmpty()) {
            throw new IllegalArgumentException(
                    "at least one pixel spectrum is required");
        }

        Path ibd = ibdPath != null ? ibdPath : deriveIbdPath(imzmlPath);

        String uuid = uuidHex == null ? UUID.randomUUID().toString()
                                       : uuidHex;
        uuid = normaliseUuid(uuid);
        if (uuid.length() != 32) {
            throw new IllegalArgumentException(
                    "uuidHex must be 32 hex chars after normalisation, got: " + uuid.length());
        }

        int gx = gridMaxX, gy = gridMaxY, gz = gridMaxZ;
        if (gx == 0) gx = pixels.stream().mapToInt(p -> p.x()).max().orElse(0);
        if (gy == 0) gy = pixels.stream().mapToInt(p -> p.y()).max().orElse(0);
        if (gz == 0) gz = 1;

        // ── .ibd assembly ──────────────────────────────────────────
        List<byte[]> ibdChunks = new ArrayList<>();
        ibdChunks.add(hexToBytes(uuid));   // 16-byte header
        long cursor = 16;

        int[][] offsets = new int[pixels.size()][4];  // {mzOff, mzLen, inOff, inLen}
        long totalPixels = pixels.size();

        if ("continuous".equals(mode)) {
            double[] sharedMz = pixels.get(0).mz();
            byte[] sharedMzBytes = doubleArrayToLE(sharedMz);
            long mzOffset = cursor;
            int mzLen = sharedMz.length;
            ibdChunks.add(sharedMzBytes);
            cursor += sharedMzBytes.length;
            for (int i = 0; i < pixels.size(); i++) {
                PixelSpectrum p = pixels.get(i);
                if (!Arrays.equals(p.mz(), sharedMz)) {
                    throw new IllegalArgumentException(
                            "continuous-mode imzML requires all pixels to share the same m/z axis; " +
                            "pixel " + i + " (x=" + p.x() + ", y=" + p.y() + ") differs");
                }
                byte[] intenBytes = doubleArrayToLE(p.intensity());
                long intOffset = cursor;
                ibdChunks.add(intenBytes);
                cursor += intenBytes.length;
                offsets[i] = new int[] {
                    (int) mzOffset, mzLen,
                    (int) intOffset, p.intensity().length
                };
                long done = (long) (i + 1);
                if (done % PROGRESS_INTERVAL_SPECTRA == 0 && done < totalPixels) {
                    progress.onProgress(done, totalPixels);
                }
            }
        } else {
            for (int i = 0; i < pixels.size(); i++) {
                PixelSpectrum p = pixels.get(i);
                if (p.mz().length != p.intensity().length) {
                    throw new IllegalArgumentException(
                            "processed-mode pixel (x=" + p.x() + ", y=" + p.y() +
                            "): mz and intensity arrays must be the same length");
                }
                byte[] mzBytes = doubleArrayToLE(p.mz());
                long mzOffset = cursor;
                ibdChunks.add(mzBytes);
                cursor += mzBytes.length;
                byte[] intenBytes = doubleArrayToLE(p.intensity());
                long intOffset = cursor;
                ibdChunks.add(intenBytes);
                cursor += intenBytes.length;
                offsets[i] = new int[] {
                    (int) mzOffset, p.mz().length,
                    (int) intOffset, p.intensity().length
                };
                long done = (long) (i + 1);
                if (done % PROGRESS_INTERVAL_SPECTRA == 0 && done < totalPixels) {
                    progress.onProgress(done, totalPixels);
                }
            }
        }

        byte[] ibdBytes = concatBytes(ibdChunks);
        try {
            Files.write(ibd, ibdBytes, StandardOpenOption.CREATE,
                    StandardOpenOption.TRUNCATE_EXISTING,
                    StandardOpenOption.WRITE);
        } catch (IOException e) {
            throw new UncheckedIOException("Failed to write .ibd: " + ibd, e);
        }
        String ibdSha1 = sha1Hex(ibdBytes);

        // ── .imzML XML ─────────────────────────────────────────────
        String xml = buildXml(uuid, ibdSha1, mode,
                gx, gy, gz, pixelSizeX, pixelSizeY,
                scanPattern, pixels, offsets);
        try {
            Files.writeString(imzmlPath, xml);
        } catch (IOException e) {
            throw new UncheckedIOException("Failed to write .imzML: " + imzmlPath, e);
        }

        progress.onProgress(totalPixels, totalPixels);
        return new WriteResult(imzmlPath, ibd, uuid, mode, pixels.size());
    }

    /** Round-trip helper: re-emit an {@link ImzMLImport}. */
    public static WriteResult writeFromImport(ImzMLImport imp, Path imzmlPath, Path ibdPath) {
        return write(
                imp.spectra(), imzmlPath, ibdPath, imp.mode(),
                imp.gridMaxX(), imp.gridMaxY(), imp.gridMaxZ(),
                imp.pixelSizeX(), imp.pixelSizeY(),
                imp.scanPattern(), imp.uuidHex());
    }

    /**
     * (M102) First mass-spectrum run of {@code ds} that carries pixel
     * coordinates ({@link ImzMLReader#pixelCoordinatesOf}: the
     * {@code spectrum_index/pixel_*} columns, else the legacy
     * {@code imzml_pixel_coordinates_csv} provenance parameter), in run
     * order; {@code null} when there is none.
     */
    public static AcquisitionRun findPixelRun(SpectralDataset ds) {
        for (AcquisitionRun run : ds.msRuns().values()) {
            if (!run.hasChannel("mz") || !run.hasChannel("intensity")) continue;
            if (ImzMLReader.pixelCoordinatesOf(run, ds.provenanceRecords()) != null) {
                return run;
            }
        }
        return null;
    }

    /**
     * (M102) Write a pixel run — one spectrum per pixel, as the imzML
     * importer stores a processed-mode file — as an imzML + .ibd pair.
     *
     * <p>Pixels are written in spectrum order with the coordinates from
     * {@link ImzMLReader#pixelCoordinatesOf} and each spectrum's
     * {@code mz} / {@code intensity} arrays. The mode is {@code mode}
     * when non-null, else the {@code imzml_mode} provenance parameter,
     * else {@code "continuous"} when every spectrum's m/z array is
     * identical and {@code "processed"} otherwise. UUID, grid extents,
     * pixel size and scan pattern come from the {@code imzml_*}
     * provenance parameters (run-level first, then
     * {@code datasetProvenance}) when present, else the defaults of
     * {@link #write} (random UUID, extents derived from the
     * coordinates, no pixel size, {@code "flyback"}).</p>
     *
     * @param run               the pixel run
     * @param datasetProvenance dataset-level provenance records (may be {@code null})
     * @param imzmlPath         destination {@code .imzML}
     * @param ibdPath           destination {@code .ibd}; {@code null} derives it
     * @param mode              {@code "continuous"}, {@code "processed"} or {@code null}
     * @param progress          optional progress sink
     * @throws IllegalArgumentException when the run has no pixel coordinates
     */
    public static WriteResult writeFromRun(AcquisitionRun run,
                                           List<ProvenanceRecord> datasetProvenance,
                                           Path imzmlPath, Path ibdPath,
                                           String mode, ProgressSink progress) {
        int[][] coords = ImzMLReader.pixelCoordinatesOf(run, datasetProvenance);
        if (coords == null) {
            throw new IllegalArgumentException(
                "run '" + run.name() + "' carries no pixel coordinates");
        }
        SpectrumIndex idx = run.spectrumIndex();
        Map<String, double[]> channels = run.channels();
        double[] mzAll = channels.get("mz");
        double[] inAll = channels.get("intensity");
        if (mzAll == null || inAll == null) {
            throw new IllegalArgumentException(
                "run '" + run.name() + "' has no mz/intensity channels");
        }
        int n = idx.count();
        List<PixelSpectrum> pixels = new ArrayList<>(n);
        boolean sameMz = true;
        double[] firstMz = null;
        for (int i = 0; i < n; i++) {
            int off = (int) idx.offsetAt(i);
            int len = idx.lengthAt(i);
            double[] mz = Arrays.copyOfRange(mzAll, off, off + len);
            double[] in = Arrays.copyOfRange(inAll, off, off + len);
            if (firstMz == null) firstMz = mz;
            else if (sameMz && !Arrays.equals(firstMz, mz)) sameMz = false;
            pixels.add(new PixelSpectrum(coords[i][0], coords[i][1], coords[i][2], mz, in));
        }

        String effMode = mode;
        if (effMode == null) {
            String p = ImzMLReader.imzmlParameter(run, datasetProvenance, ImzMLReader.PARAM_MODE);
            if ("processed".equals(p)) effMode = p;
            // A recorded "continuous" holds only while the arrays still
            // share one m/z axis; otherwise processed is the only
            // faithful encoding.
            else if ("continuous".equals(p)) effMode = sameMz ? p : "processed";
        }
        if (effMode == null) effMode = sameMz ? "continuous" : "processed";

        String uuid = ImzMLReader.imzmlParameter(run, datasetProvenance, ImzMLReader.PARAM_UUID_HEX);
        if (uuid != null && normaliseUuid(uuid).length() != 32) uuid = null;
        int gx = intParam(run, datasetProvenance, ImzMLReader.PARAM_GRID_MAX_X);
        int gy = intParam(run, datasetProvenance, ImzMLReader.PARAM_GRID_MAX_Y);
        int gz = intParam(run, datasetProvenance, ImzMLReader.PARAM_GRID_MAX_Z);
        double psx = doubleParam(run, datasetProvenance, ImzMLReader.PARAM_PIXEL_SIZE_X);
        double psy = doubleParam(run, datasetProvenance, ImzMLReader.PARAM_PIXEL_SIZE_Y);
        String scan = ImzMLReader.imzmlParameter(run, datasetProvenance, ImzMLReader.PARAM_SCAN_PATTERN);
        if (scan == null || scan.isEmpty()) scan = "flyback";

        return write(pixels, imzmlPath, ibdPath, effMode,
            gx, gy, gz, psx, psy, scan, uuid, progress);
    }

    private static int intParam(AcquisitionRun run, List<ProvenanceRecord> ds, String key) {
        String v = ImzMLReader.imzmlParameter(run, ds, key);
        if (v == null || v.isBlank()) return 0;
        try {
            return Integer.parseInt(v.trim());
        } catch (NumberFormatException e) {
            try {
                return (int) Double.parseDouble(v.trim());
            } catch (NumberFormatException e2) {
                return 0;
            }
        }
    }

    private static double doubleParam(AcquisitionRun run, List<ProvenanceRecord> ds, String key) {
        String v = ImzMLReader.imzmlParameter(run, ds, key);
        if (v == null || v.isBlank()) return 0.0;
        try {
            return Double.parseDouble(v.trim());
        } catch (NumberFormatException e) {
            return 0.0;
        }
    }

    // ── helpers ────────────────────────────────────────────────────

    private static Path deriveIbdPath(Path imzmlPath) {
        String name = imzmlPath.getFileName().toString();
        int dot = name.lastIndexOf('.');
        String stem = dot > 0 ? name.substring(0, dot) : name;
        return imzmlPath.resolveSibling(stem + ".ibd");
    }

    private static String normaliseUuid(String raw) {
        StringBuilder sb = new StringBuilder(raw.length());
        for (int i = 0; i < raw.length(); i++) {
            char c = raw.charAt(i);
            if (c == '{' || c == '}' || c == '-' || Character.isWhitespace(c)) continue;
            sb.append(Character.toLowerCase(c));
        }
        return sb.toString();
    }

    private static byte[] hexToBytes(String hex) {
        byte[] out = new byte[hex.length() / 2];
        for (int i = 0; i < out.length; i++) {
            out[i] = (byte) Integer.parseInt(hex.substring(i * 2, i * 2 + 2), 16);
        }
        return out;
    }

    private static byte[] doubleArrayToLE(double[] values) {
        ByteBuffer buf = ByteBuffer.allocate(values.length * 8)
                                    .order(ByteOrder.LITTLE_ENDIAN);
        for (double v : values) buf.putDouble(v);
        return buf.array();
    }

    private static byte[] concatBytes(List<byte[]> chunks) {
        int total = 0;
        for (byte[] c : chunks) total += c.length;
        byte[] out = new byte[total];
        int pos = 0;
        for (byte[] c : chunks) {
            System.arraycopy(c, 0, out, pos, c.length);
            pos += c.length;
        }
        return out;
    }

    private static String sha1Hex(byte[] data) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-1");
            byte[] digest = md.digest(data);
            StringBuilder sb = new StringBuilder(digest.length * 2);
            for (byte b : digest) sb.append(String.format("%02x", b));
            return sb.toString();
        } catch (NoSuchAlgorithmException e) {
            throw new AssertionError("SHA-1 must be available in every JDK", e);
        }
    }

    private static String xmlEscape(String s) {
        StringBuilder sb = new StringBuilder(s.length());
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '&': sb.append("&amp;"); break;
                case '<': sb.append("&lt;"); break;
                case '>': sb.append("&gt;"); break;
                case '"': sb.append("&quot;"); break;
                case '\'': sb.append("&apos;"); break;
                default: sb.append(c);
            }
        }
        return sb.toString();
    }

    private static String buildXml(
            String uuid, String ibdSha1, String mode,
            int gx, int gy, int gz,
            double pixelSizeX, double pixelSizeY,
            String scanPattern,
            List<PixelSpectrum> pixels,
            int[][] offsets
    ) {
        String modeAcc = "continuous".equals(mode) ? "IMS:1000030" : "IMS:1000031";
        String modeName = mode;

        StringBuilder sb = new StringBuilder(8192);
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"");
        sb.append(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"");
        sb.append(" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml");
        sb.append(" http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" version=\"1.1\">\n");

        sb.append("  <cvList count=\"3\">\n");
        sb.append("    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\"");
        sb.append(" version=\"4.1.0\"");
        sb.append(" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n");
        sb.append("    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"2020-03-10\"");
        sb.append(" URI=\"http://ontologies.berkeleybop.org/uo.obo\"/>\n");
        sb.append("    <cv id=\"IMS\" fullName=\"Mass Spectrometry Imaging Ontology\" version=\"1.1.0\"");
        sb.append(" URI=\"https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo\"/>\n");
        sb.append("  </cvList>\n");
        sb.append("  <fileDescription>\n");
        sb.append("    <fileContent>\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000080\"");
        sb.append(" name=\"universally unique identifier\" value=\"").append(uuid).append("\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000091\" name=\"ibd SHA-1\"");
        sb.append(" value=\"").append(ibdSha1).append("\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"").append(modeAcc);
        sb.append("\" name=\"").append(modeName).append("\" value=\"\"/>\n");
        sb.append("    </fileContent>\n");
        sb.append("  </fileDescription>\n");

        sb.append("  <referenceableParamGroupList count=\"2\">\n");
        sb.append("    <referenceableParamGroup id=\"mzArray\">\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"");
        sb.append(" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\"");
        sb.append(" name=\"external data\" value=\"true\"/>\n");
        sb.append("    </referenceableParamGroup>\n");
        sb.append("    <referenceableParamGroup id=\"intensityArray\">\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"");
        sb.append(" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\"");
        sb.append(" name=\"external data\" value=\"true\"/>\n");
        sb.append("    </referenceableParamGroup>\n");
        sb.append("  </referenceableParamGroupList>\n");

        sb.append("  <softwareList count=\"1\">\n");
        sb.append("    <software id=\"ttio\" version=\"0.9.0\">\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000799\"");
        sb.append(" name=\"custom unreleased software tool\" value=\"ttio\"/>\n");
        sb.append("    </software>\n");
        sb.append("  </softwareList>\n");

        sb.append("  <scanSettingsList count=\"1\">\n");
        sb.append("    <scanSettings id=\"scansettings1\">\n");
        sb.append("      <userParam name=\"scan pattern\" value=\"").append(xmlEscape(scanPattern));
        sb.append("\" type=\"xsd:string\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000040\" name=\"linescan sequence\"");
        sb.append(" value=\"").append(xmlEscape(scanPattern)).append("\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000042\"");
        sb.append(" name=\"max count of pixels x\" value=\"").append(gx).append("\"/>\n");
        sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000043\"");
        sb.append(" name=\"max count of pixels y\" value=\"").append(gy).append("\"/>\n");
        if (pixelSizeX > 0.0) {
            sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000046\"");
            sb.append(" name=\"pixel size (x)\" value=\"").append(pixelSizeX).append("\"");
            sb.append(" unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n");
        }
        if (pixelSizeY > 0.0) {
            sb.append("      <cvParam cvRef=\"IMS\" accession=\"IMS:1000047\"");
            sb.append(" name=\"pixel size y\" value=\"").append(pixelSizeY).append("\"");
            sb.append(" unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n");
        }
        sb.append("    </scanSettings>\n");
        sb.append("  </scanSettingsList>\n");

        sb.append("  <instrumentConfigurationList count=\"1\">\n");
        sb.append("    <instrumentConfiguration id=\"IC1\">\n");
        sb.append("      <cvParam cvRef=\"MS\" accession=\"MS:1000031\"");
        sb.append(" name=\"instrument model\" value=\"\"/>\n");
        sb.append("    </instrumentConfiguration>\n");
        sb.append("  </instrumentConfigurationList>\n");

        sb.append("  <dataProcessingList count=\"1\">\n");
        sb.append("    <dataProcessing id=\"dp_export\">\n");
        sb.append("      <processingMethod order=\"0\" softwareRef=\"ttio\">\n");
        sb.append("        <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\n");
        sb.append("      </processingMethod>\n");
        sb.append("    </dataProcessing>\n");
        sb.append("  </dataProcessingList>\n");

        sb.append("  <run id=\"ttio_imzml_export\" defaultInstrumentConfigurationRef=\"IC1\">\n");
        sb.append("    <spectrumList count=\"").append(pixels.size());
        sb.append("\" defaultDataProcessingRef=\"dp_export\">\n");

        for (int i = 0; i < pixels.size(); i++) {
            PixelSpectrum p = pixels.get(i);
            int mzOff = offsets[i][0], mzLen = offsets[i][1];
            int inOff = offsets[i][2], inLen = offsets[i][3];
            int mzEnc = mzLen * 8, inEnc = inLen * 8;

            sb.append("      <spectrum id=\"Scan=").append(i + 1);
            sb.append("\" index=\"").append(i).append("\" defaultArrayLength=\"0\">\n");
            sb.append("        <cvParam cvRef=\"MS\" accession=\"MS:1000579\"");
            sb.append(" name=\"MS1 spectrum\" value=\"\"/>\n");
            sb.append("        <cvParam cvRef=\"MS\" accession=\"MS:1000511\"");
            sb.append(" name=\"ms level\" value=\"1\"/>\n");
            sb.append("        <scanList count=\"1\">\n");
            sb.append("          <cvParam cvRef=\"MS\" accession=\"MS:1000795\"");
            sb.append(" name=\"no combination\" value=\"\"/>\n");
            sb.append("          <scan instrumentConfigurationRef=\"IC1\">\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000050\"");
            sb.append(" name=\"position x\" value=\"").append(p.x()).append("\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000051\"");
            sb.append(" name=\"position y\" value=\"").append(p.y()).append("\"/>\n");
            if (p.z() != 1) {
                sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000052\"");
                sb.append(" name=\"position z\" value=\"").append(p.z()).append("\"/>\n");
            }
            sb.append("          </scan>\n");
            sb.append("        </scanList>\n");
            sb.append("        <binaryDataArrayList count=\"2\">\n");

            // m/z array
            sb.append("          <binaryDataArray encodedLength=\"0\">\n");
            sb.append("            <referenceableParamGroupRef ref=\"mzArray\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000523\"");
            sb.append(" name=\"64-bit float\" value=\"\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000576\"");
            sb.append(" name=\"no compression\" value=\"\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000514\"");
            sb.append(" name=\"m/z array\" value=\"\"");
            sb.append(" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\"");
            sb.append(" name=\"external data\" value=\"true\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000102\"");
            sb.append(" name=\"external offset\" value=\"").append(mzOff).append("\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000103\"");
            sb.append(" name=\"external array length\" value=\"").append(mzLen).append("\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000104\"");
            sb.append(" name=\"external encoded length\" value=\"").append(mzEnc).append("\"/>\n");
            sb.append("            <binary/>\n");
            sb.append("          </binaryDataArray>\n");

            // intensity array
            sb.append("          <binaryDataArray encodedLength=\"0\">\n");
            sb.append("            <referenceableParamGroupRef ref=\"intensityArray\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000523\"");
            sb.append(" name=\"64-bit float\" value=\"\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000576\"");
            sb.append(" name=\"no compression\" value=\"\"/>\n");
            sb.append("            <cvParam cvRef=\"MS\" accession=\"MS:1000515\"");
            sb.append(" name=\"intensity array\" value=\"\"");
            sb.append(" unitCvRef=\"MS\" unitAccession=\"MS:1000131\"");
            sb.append(" unitName=\"number of detector counts\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\"");
            sb.append(" name=\"external data\" value=\"true\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000102\"");
            sb.append(" name=\"external offset\" value=\"").append(inOff).append("\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000103\"");
            sb.append(" name=\"external array length\" value=\"").append(inLen).append("\"/>\n");
            sb.append("            <cvParam cvRef=\"IMS\" accession=\"IMS:1000104\"");
            sb.append(" name=\"external encoded length\" value=\"").append(inEnc).append("\"/>\n");
            sb.append("            <binary/>\n");
            sb.append("          </binaryDataArray>\n");
            sb.append("        </binaryDataArrayList>\n");
            sb.append("      </spectrum>\n");
        }

        sb.append("    </spectrumList>\n");
        sb.append("  </run>\n");
        sb.append("</mzML>\n");
        return sb.toString();
    }
}
