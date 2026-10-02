/*
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */
package global.thalion.ttio.importers;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums;
import global.thalion.ttio.InstrumentConfig;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.SpectrumIndex;
import global.thalion.ttio.io.ProgressSink;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.time.Instant;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * imzML + .ibd importer
 *
 * <p>imzML is the dominant interchange format for mass-spectrometry
 * imaging. The format is a pair of files:
 * <ul>
 *   <li>{@code <stem>.imzML} — XML metadata mirroring mzML, with
 *       each {@code <spectrum>} carrying an external offset /
 *       external array length cvParam triple
 *       ({@code IMS:1000102} / {@code IMS:1000103} /
 *       {@code IMS:1000104}).</li>
 *   <li>{@code <stem>.ibd} — concatenated binary mass / intensity
 *       arrays prefixed by a 16-byte UUID that must match the
 *       {@code IMS:1000042 universally unique identifier} cvParam
 *       in the metadata </li>
 * </ul>
 *
 * <p>Two storage modes:
 * <ul>
 *   <li><b>continuous</b> ({@code IMS:1000030}) — single shared m/z
 *       array stored once; per-pixel intensity arrays follow.</li>
 *   <li><b>processed</b> ({@code IMS:1000031}) — per-pixel m/z +
 *       intensity arrays.</li>
 * </ul>
 *
 * <p>API status: Provisional ().</p>
 *
 * <p>Cross-language equivalents: Python
 * {@code ttio.importers.imzml}, Objective-C
 * {@code TTIOImzMLReader}.</p>
 *
 *
 */
public final class ImzMLReader {

    /** Result of parsing an imzML + .ibd pair. */
    public record ImzMLImport(
        String mode,
        String uuidHex,
        int gridMaxX,
        int gridMaxY,
        int gridMaxZ,
        double pixelSizeX,
        double pixelSizeY,
        String scanPattern,
        List<PixelSpectrum> spectra,
        String sourceImzML,
        String sourceIbd
    ) {
        /** The imzML provenance record shared by the pixel run and the
         *  dataset: the scalar {@code imzml_*} parameters (M102 — the
         *  coordinates themselves live in {@code spectrum_index/pixel_*}),
         *  software {@code "ttio imzml importer v0.9"}, input refs
         *  {@code [imzML, ibd]}.
         *
         *  @param outputRef recorded as the only output ref when
         *                   non-null; {@code null} records none. */
        public ProvenanceRecord provenanceRecord(String outputRef) {
            Map<String, String> params = new LinkedHashMap<>();
            params.put(PARAM_MODE, mode);
            params.put(PARAM_UUID_HEX, uuidHex);
            params.put(PARAM_GRID_MAX_X, Integer.toString(gridMaxX));
            params.put(PARAM_GRID_MAX_Y, Integer.toString(gridMaxY));
            params.put(PARAM_GRID_MAX_Z, Integer.toString(gridMaxZ));
            params.put(PARAM_PIXEL_SIZE_X, Double.toString(pixelSizeX));
            params.put(PARAM_PIXEL_SIZE_Y, Double.toString(pixelSizeY));
            params.put(PARAM_SCAN_PATTERN, scanPattern);
            return new ProvenanceRecord(Instant.now().getEpochSecond(),
                IMPORTER_SOFTWARE, params,
                List.of(sourceImzML, sourceIbd),
                outputRef != null ? List.of(outputRef) : List.of());
        }

        /** Build the {@code imzml_pixels} run: one mass spectrum per
         *  pixel in file order, {@code mz} + {@code intensity} channels
         *  concatenated, retention time 0, MS level 1, polarity unknown,
         *  no precursor, base peak = the pixel's maximum intensity, and
         *  the pixel coordinates as {@code spectrum_index/pixel_x/y/z}
         *  (M102). Mirrors the run built by Python
         *  {@code ImzMLImport.to_imported_dataset}.
         *
         *  @param provenance the run-level provenance record (normally
         *                    {@link #provenanceRecord(String)}). */
        public AcquisitionRun toPixelRun(ProvenanceRecord provenance) {
            int n = spectra.size();
            if (n == 0) {
                throw new IllegalStateException(sourceImzML + ": no spectra parsed");
            }
            long[] offsets = new long[n];
            int[] lengths = new int[n];
            long total = 0;
            for (int i = 0; i < n; i++) {
                offsets[i] = total;
                lengths[i] = spectra.get(i).mz().length;
                total += lengths[i];
            }
            if (total > Integer.MAX_VALUE) {
                throw new IllegalStateException(
                    sourceImzML + ": too many peaks for one in-memory run (" + total + ")");
            }
            double[] mzBuf = new double[(int) total];
            double[] itBuf = new double[(int) total];
            double[] basePeaks = new double[n];
            int[] px = new int[n], py = new int[n], pz = new int[n];
            int[] msLevels = new int[n];
            for (int i = 0; i < n; i++) {
                PixelSpectrum s = spectra.get(i);
                System.arraycopy(s.mz(), 0, mzBuf, (int) offsets[i], lengths[i]);
                System.arraycopy(s.intensity(), 0, itBuf, (int) offsets[i], lengths[i]);
                double bp = 0.0;
                if (s.intensity().length > 0) {
                    bp = Double.NEGATIVE_INFINITY;
                    for (double v : s.intensity()) if (v > bp) bp = v;
                }
                basePeaks[i] = bp;
                px[i] = s.x();
                py[i] = s.y();
                pz[i] = s.z();
                msLevels[i] = 1;
            }
            SpectrumIndex index = new SpectrumIndex(n, offsets, lengths,
                new double[n], msLevels, new int[n],
                new double[n], new int[n], basePeaks)
                .withPixelCoordinates(px, py, pz);
            Map<String, double[]> channels = new LinkedHashMap<>();
            channels.put("mz", mzBuf);
            channels.put("intensity", itBuf);
            return new AcquisitionRun(PIXEL_RUN_NAME,
                Enums.AcquisitionMode.MS1_DDA, index,
                new InstrumentConfig("", "", "", "", "", ""),
                channels, List.of(),
                provenance != null ? List.of(provenance) : List.of(),
                null, 0.0);
        }

        /** Build the normalized import draft: the {@code imzml_pixels}
         *  run, title {@code "imzML import: <file name>"}, and the same
         *  provenance record at run and dataset level. Mirrors Python
         *  {@code ImzMLImport.to_imported_dataset}.
         *
         *  @param outputRef optional output ref for the provenance record. */
        public ImportedDataset toImportedDataset(String outputRef) {
            ProvenanceRecord prov = provenanceRecord(outputRef);
            ImportedDataset d = new ImportedDataset();
            d.title = "imzML import: " + Path.of(sourceImzML).getFileName();
            d.runs.add(toPixelRun(prov));
            d.provenance.add(prov);
            return d;
        }
    }

    /** Name of the run the importer builds from an imzML file's pixels. */
    public static final String PIXEL_RUN_NAME = "imzml_pixels";
    /** Provenance {@code software} string of the imzML importer. */
    public static final String IMPORTER_SOFTWARE = "ttio imzml importer v0.9";
    /** Provenance parameter keys (format-spec §4b). */
    public static final String PARAM_MODE = "imzml_mode";
    public static final String PARAM_UUID_HEX = "imzml_uuid_hex";
    public static final String PARAM_GRID_MAX_X = "imzml_grid_max_x";
    public static final String PARAM_GRID_MAX_Y = "imzml_grid_max_y";
    public static final String PARAM_GRID_MAX_Z = "imzml_grid_max_z";
    public static final String PARAM_PIXEL_SIZE_X = "imzml_pixel_size_x";
    public static final String PARAM_PIXEL_SIZE_Y = "imzml_pixel_size_y";
    public static final String PARAM_SCAN_PATTERN = "imzml_scan_pattern";
    /** Pre-M102 coordinate parameter ({@code "x,y,z;x,y,z;..."}); read
     *  for older files only, never written. */
    public static final String PARAM_LEGACY_COORDINATES_CSV = "imzml_pixel_coordinates_csv";

    /**
     * (M102) Pixel coordinates of {@code run}'s spectra, or {@code null}
     * when it has none.
     *
     * <p>The run's {@code spectrum_index/pixel_x/y/z} columns win. For a
     * run without them (files written before M102) the legacy
     * {@code imzml_pixel_coordinates_csv} provenance parameter is
     * looked up in the run's own provenance records, then in
     * {@code datasetProvenance}, and used only when it parses to exactly
     * one triple per spectrum (format-spec §4b).</p>
     *
     * @param run               the run to inspect
     * @param datasetProvenance the dataset-level provenance records; may be {@code null}
     * @return an {@code N x 3} array, row {@code i} = {@code {x, y, z}}
     *         of spectrum {@code i}; or {@code null}
     */
    public static int[][] pixelCoordinatesOf(AcquisitionRun run,
                                             List<ProvenanceRecord> datasetProvenance) {
        if (run == null || run.spectrumIndex() == null) return null;
        SpectrumIndex idx = run.spectrumIndex();
        int n = idx.count();
        if (idx.hasPixelCoordinates()) {
            int[][] out = new int[n][];
            for (int i = 0; i < n; i++) out[i] = idx.pixelCoordinatesAt(i);
            return out;
        }
        int[][] fromRun = legacyCoordinates(run.provenanceRecords(), n);
        if (fromRun != null) return fromRun;
        return legacyCoordinates(datasetProvenance, n);
    }

    /** (M102) Value of the imzML provenance parameter {@code key} for
     *  {@code run}: the run's own records first, then
     *  {@code datasetProvenance}; {@code null} when absent. */
    public static String imzmlParameter(AcquisitionRun run,
                                        List<ProvenanceRecord> datasetProvenance,
                                        String key) {
        if (run != null) {
            for (ProvenanceRecord r : run.provenanceRecords()) {
                String v = r.parameters().get(key);
                if (v != null) return v;
            }
        }
        if (datasetProvenance != null) {
            for (ProvenanceRecord r : datasetProvenance) {
                String v = r.parameters().get(key);
                if (v != null) return v;
            }
        }
        return null;
    }

    private static int[][] legacyCoordinates(List<ProvenanceRecord> records, int n) {
        if (records == null) return null;
        for (ProvenanceRecord r : records) {
            String csv = r.parameters().get(PARAM_LEGACY_COORDINATES_CSV);
            if (csv == null) continue;
            int[][] parsed = parseLegacyCsv(csv, n);
            if (parsed != null) return parsed;
        }
        return null;
    }

    /** Parse {@code "x,y,z;x,y,z;..."}; {@code null} unless it holds
     *  exactly {@code n} well-formed triples. */
    static int[][] parseLegacyCsv(String csv, int n) {
        String trimmed = csv.trim();
        if (trimmed.isEmpty()) return n == 0 ? new int[0][] : null;
        String[] triples = trimmed.split(";", -1);
        if (triples.length != n) return null;
        int[][] out = new int[n][];
        for (int i = 0; i < n; i++) {
            String[] parts = triples[i].split(",", -1);
            if (parts.length != 3) return null;
            try {
                out[i] = new int[]{
                    Integer.parseInt(parts[0].trim()),
                    Integer.parseInt(parts[1].trim()),
                    Integer.parseInt(parts[2].trim())
                };
            } catch (NumberFormatException e) {
                return null;
            }
        }
        return out;
    }

    /** One pixel's spatial coordinates plus its m/z + intensity arrays. */
    public record PixelSpectrum(
        int x, int y, int z,
        double[] mz,
        double[] intensity
    ) {}

    /** Raised when the imzML XML is structurally invalid for our needs. */
    public static final class ImzMLParseException extends IOException {
        private static final long serialVersionUID = 1L;
        public ImzMLParseException(String msg) { super(msg); }
        public ImzMLParseException(String msg, Throwable cause) { super(msg, cause); }
    }

    /** Raised when the .ibd binary disagrees with the .imzML metadata
     *  (UUID mismatch or offset overflow). */
    public static final class ImzMLBinaryException extends IOException {
        private static final long serialVersionUID = 1L;
        public ImzMLBinaryException(String msg) { super(msg); }
    }

    /** Default emit-every-N cadence for {@link ProgressSink} callbacks
     *  during the .ibd materialisation loop. Imaging-MS files often
     *  have tens of thousands of pixels, so reporting per-pixel would
     *  drown the UI thread; 100 pixels keeps overhead negligible while
     *  still giving sub-second updates on typical files. */
    public static final int PROGRESS_INTERVAL_SPECTRA = 100;

    private ImzMLReader() {}

    /** Parse an imzML metadata file, locating the sibling .ibd via
     *  filename rewriting. */
    public static ImzMLImport read(Path imzmlPath) throws IOException {
        return read(imzmlPath, null, ProgressSink.discard());
    }

    /** Parse an imzML metadata file with an explicit .ibd location. */
    public static ImzMLImport read(Path imzmlPath, Path ibdPath) throws IOException {
        return read(imzmlPath, ibdPath, ProgressSink.discard());
    }

    /** Overload of {@link #read(Path)} that accepts a
     *  {@link ProgressSink}. */
    public static ImzMLImport read(Path imzmlPath, ProgressSink progress)
            throws IOException {
        return read(imzmlPath, null, progress);
    }

    /**
     * Parse an imzML metadata file with an explicit .ibd location,
     * firing {@code progress.onProgress(pixelsDone, totalPixels)} every
     * {@link #PROGRESS_INTERVAL_SPECTRA} pixels during the .ibd
     * materialisation phase. The total is known up front (one stub per
     * {@code <spectrum>}), so callbacks carry a real {@code total} from
     * the first fire instead of {@code -1L}.
     *
     * @since 1.5.0
     */
    public static ImzMLImport read(Path imzmlPath, Path ibdPath,
                                    ProgressSink progress) throws IOException {
        if (progress == null) progress = ProgressSink.discard();
        if (!Files.isRegularFile(imzmlPath)) {
            throw new ImzMLParseException("imzML metadata not found: " + imzmlPath);
        }
        Path resolvedIbd = (ibdPath != null) ? ibdPath : siblingIbd(imzmlPath);
        if (!Files.isRegularFile(resolvedIbd)) {
            throw new ImzMLParseException("imzML binary not found: " + resolvedIbd);
        }

        Handler handler = new Handler();
        try {
            SAXParserFactory factory = SAXParserFactory.newInstance();
            factory.setNamespaceAware(true);
            SAXParser parser = factory.newSAXParser();
            parser.parse(imzmlPath.toFile(), handler);
        } catch (SAXException | javax.xml.parsers.ParserConfigurationException e) {
            throw new ImzMLParseException(
                "failed to parse imzML '" + imzmlPath + "': " + e.getMessage(), e);
        }

        if (handler.mode.isEmpty()) {
            throw new ImzMLParseException(
                imzmlPath + ": no continuous/processed mode CV term found");
        }
        if (handler.uuidHex.isEmpty()) {
            throw new ImzMLParseException(
                imzmlPath + ": missing IMS:1000042 universally unique identifier");
        }
        if (handler.stubs.isEmpty()) {
            throw new ImzMLParseException(imzmlPath + ": no <spectrum> elements parsed");
        }

        long ibdSize = Files.size(resolvedIbd);
        if (ibdSize < 16) {
            throw new ImzMLBinaryException(
                resolvedIbd + ": shorter than the 16-byte UUID header");
        }

        try (RandomAccessFile raf = new RandomAccessFile(resolvedIbd.toFile(), "r")) {
            byte[] uuidBytes = new byte[16];
            raf.readFully(uuidBytes);
            String ibdUuidHex = HexFormat.of().formatHex(uuidBytes);
            if (!ibdUuidHex.equals(handler.uuidHex)) {
                throw new ImzMLBinaryException(
                    "UUID mismatch: imzML declares " + handler.uuidHex
                  + " but .ibd header is " + ibdUuidHex);
            }

            List<PixelSpectrum> pixels = new ArrayList<>(handler.stubs.size());
            long totalPixels = handler.stubs.size();
            double[] sharedMz = null;
            boolean continuous = "continuous".equals(handler.mode);
            long pixelsDone = 0L;
            for (Stub stub : handler.stubs) {
                double[] mz = readArray(raf, stub.mzOffset, stub.mzLength,
                                         stub.mzPrecision, ibdSize, resolvedIbd, "m/z");
                double[] intensity = readArray(raf, stub.intOffset, stub.intLength,
                                                stub.intPrecision, ibdSize, resolvedIbd, "intensity");
                if (mz.length != intensity.length) {
                    throw new ImzMLBinaryException(
                        resolvedIbd + ": pixel (" + stub.x + "," + stub.y
                      + ") mz/intensity size mismatch (" + mz.length
                      + " vs " + intensity.length + ")");
                }
                double[] effectiveMz;
                if (continuous) {
                    if (sharedMz == null) sharedMz = mz;
                    effectiveMz = sharedMz;
                } else {
                    effectiveMz = mz;
                }
                pixels.add(new PixelSpectrum(stub.x, stub.y, stub.z,
                                              effectiveMz, intensity));
                pixelsDone++;
                if (pixelsDone % PROGRESS_INTERVAL_SPECTRA == 0) {
                    progress.onProgress(pixelsDone, totalPixels);
                }
            }
            // Final fire — total is known up front for imzML.
            progress.onProgress(pixelsDone, totalPixels);

            return new ImzMLImport(
                handler.mode, handler.uuidHex,
                handler.gridMaxX, handler.gridMaxY, handler.gridMaxZ,
                handler.pixelSizeX, handler.pixelSizeY,
                handler.scanPattern,
                List.copyOf(pixels),
                imzmlPath.toString(), resolvedIbd.toString()
            );
        }
    }

    // ────────────────────────────────────────────────────────────────────
    // Internal: SAX handler + binary materialisation.
    // ────────────────────────────────────────────────────────────────────

    private static Path siblingIbd(Path imzmlPath) {
        String name = imzmlPath.getFileName().toString();
        int dot = name.lastIndexOf('.');
        String stem = (dot > 0) ? name.substring(0, dot) : name;
        Path parent = imzmlPath.getParent();
        return (parent != null) ? parent.resolve(stem + ".ibd")
                                : Paths.get(stem + ".ibd");
    }

    private static String normaliseUuid(String value) {
        return value.replace("{", "").replace("}", "")
                    .replace("-", "").trim().toLowerCase();
    }

    private static double[] readArray(RandomAccessFile raf,
                                       long offset, long length,
                                       String precision,
                                       long ibdSize, Path ibdPath, String label)
            throws IOException
    {
        if (offset < 0 || length < 0) {
            throw new ImzMLBinaryException(
                ibdPath + ": negative offset/length for " + label + " array");
        }
        if (length == 0) return new double[0];
        int bytesPer = "64".equals(precision) ? 8 : 4;
        long nbytes = length * bytesPer;
        if (offset + nbytes > ibdSize) {
            throw new ImzMLBinaryException(
                ibdPath + ": " + label + " array reads past end of file "
              + "(offset=" + offset + ", bytes=" + nbytes + ", size=" + ibdSize + ")");
        }
        byte[] raw = new byte[(int) nbytes];
        raf.seek(offset);
        raf.readFully(raw);
        ByteBuffer buf = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        double[] out = new double[(int) length];
        if (bytesPer == 8) {
            for (int i = 0; i < length; i++) out[i] = buf.getDouble();
        } else {
            for (int i = 0; i < length; i++) out[i] = buf.getFloat();
        }
        return out;
    }

    private static final class Stub {
        int x = 0, y = 0, z = 1;
        long mzOffset = -1, mzLength = 0;
        long intOffset = -1, intLength = 0;
        String mzPrecision = "64";
        String intPrecision = "64";
    }

    private static final class Handler extends DefaultHandler {
        String mode = "";
        String uuidHex = "";
        int gridMaxX = 0, gridMaxY = 0, gridMaxZ = 1;
        double pixelSizeX = 0.0, pixelSizeY = 0.0;
        String scanPattern = "";
        final List<Stub> stubs = new ArrayList<>();

        private Stub current;
        private boolean inBinaryArray;
        private boolean inScan;
        private String arrayKind = "";

        @Override
        public void startElement(String uri, String localName,
                                  String qName, Attributes attrs)
        {
            String local = localName(localName, qName);
            switch (local) {
                case "spectrum":
                    current = new Stub();
                    break;
                case "binaryDataArray":
                    inBinaryArray = true;
                    arrayKind = "";
                    break;
                case "scan":
                    inScan = true;
                    break;
                case "cvParam":
                    handleCvParam(attrs);
                    break;
                default:
                    break;
            }
        }

        @Override
        public void endElement(String uri, String localName, String qName) {
            String local = localName(localName, qName);
            switch (local) {
                case "spectrum":
                    if (current != null) stubs.add(current);
                    current = null;
                    break;
                case "binaryDataArray":
                    inBinaryArray = false;
                    arrayKind = "";
                    break;
                case "scan":
                    inScan = false;
                    break;
                default:
                    break;
            }
        }

        private static String localName(String localName, String qName) {
            if (localName != null && !localName.isEmpty()) return localName;
            int colon = qName.indexOf(':');
            return (colon >= 0) ? qName.substring(colon + 1) : qName;
        }

        private void handleCvParam(Attributes attrs) {
            String acc = attrs.getValue("accession");
            String value = attrs.getValue("value");
            if (acc == null) return;
            if (value == null) value = "";

            switch (acc) {
                // imzML storage mode: only the IMS-namespaced forms are
                // real. MS:1000030 = "vendor processing software",
                // MS:1000031 = "instrument model" — unrelated terms.
                case "IMS:1000030": mode = "continuous"; return;
                case "IMS:1000031": mode = "processed"; return;
                // Canonical IMS accessions (real-world imzML 1.1).
                case "IMS:1000080":
                    if (!value.isEmpty()) uuidHex = normaliseUuid(value);
                    return;
                case "IMS:1000042":
                    if (!value.isEmpty()) {
                        // Canonical max count of pixels x. Legacy TTIO
                        // pre-v0.9 used IMS:1000042 for UUID with a
                        // hex-string value; if the value isn't a plain
                        // integer and we haven't seen a UUID yet, treat
                        // it as the legacy UUID accession.
                        try {
                            gridMaxX = Integer.parseInt(value);
                        } catch (NumberFormatException ignored) {
                            if (uuidHex.isEmpty()) {
                                String cand = normaliseUuid(value);
                                if (cand.length() == 32) uuidHex = cand;
                            }
                        }
                    }
                    return;
                case "IMS:1000043":
                    if (!value.isEmpty()) gridMaxY = Integer.parseInt(value); return;
                case "IMS:1000003":
                    if (!value.isEmpty()) gridMaxX = Integer.parseInt(value); return;
                case "IMS:1000004":
                    if (!value.isEmpty()) gridMaxY = Integer.parseInt(value); return;
                case "IMS:1000005":
                    if (!value.isEmpty()) gridMaxZ = Integer.parseInt(value); return;
                case "IMS:1000046":
                    if (!value.isEmpty()) pixelSizeX = Double.parseDouble(value); return;
                case "IMS:1000047":
                    if (!value.isEmpty()) pixelSizeY = Double.parseDouble(value); return;
                case "IMS:1000040": case "IMS:1000048":
                    if (!value.isEmpty() && scanPattern.isEmpty()) scanPattern = value;
                    return;
                default:
                    break;
            }

            if (inScan && current != null) {
                switch (acc) {
                    case "IMS:1000050":
                        if (!value.isEmpty()) current.x = Integer.parseInt(value); return;
                    case "IMS:1000051":
                        if (!value.isEmpty()) current.y = Integer.parseInt(value); return;
                    case "IMS:1000052":
                        if (!value.isEmpty()) current.z = Integer.parseInt(value); return;
                    default: break;
                }
            }
            if (inBinaryArray && current != null) {
                switch (acc) {
                    case "MS:1000514":
                        arrayKind = "mz"; return;
                    case "MS:1000515":
                        arrayKind = "intensity"; return;
                    case "MS:1000523":
                        if ("mz".equals(arrayKind)) current.mzPrecision = "64";
                        else if ("intensity".equals(arrayKind)) current.intPrecision = "64";
                        return;
                    case "MS:1000521":
                        if ("mz".equals(arrayKind)) current.mzPrecision = "32";
                        else if ("intensity".equals(arrayKind)) current.intPrecision = "32";
                        return;
                    case "IMS:1000102":
                        if (value.isEmpty()) return;
                        if ("mz".equals(arrayKind)) current.mzOffset = Long.parseLong(value);
                        else if ("intensity".equals(arrayKind)) current.intOffset = Long.parseLong(value);
                        return;
                    case "IMS:1000103":
                        if (value.isEmpty()) return;
                        if ("mz".equals(arrayKind)) current.mzLength = Long.parseLong(value);
                        else if ("intensity".equals(arrayKind)) current.intLength = Long.parseLong(value);
                        return;
                    default: break;
                }
            }
        }
    }
}
