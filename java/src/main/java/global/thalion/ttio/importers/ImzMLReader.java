/*
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */
package global.thalion.ttio.importers;

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
import java.util.HashMap;
import java.util.HexFormat;
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
    ) {}

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
        // <referenceableParamGroup id=...> cvParams (accession, value),
        // applied wherever a <referenceableParamGroupRef ref=...>
        // appears. Writers commonly declare the array kind / precision /
        // compression of every binaryDataArray this way (e.g. the
        // HR2MSI PXD001283 files).
        private final Map<String, List<String[]>> paramGroups = new HashMap<>();
        private List<String[]> currentGroup;

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
                case "referenceableParamGroup":
                    currentGroup = new ArrayList<>();
                    paramGroups.put(nonNull(attrs.getValue("id")), currentGroup);
                    break;
                case "referenceableParamGroupRef":
                    for (String[] cv : paramGroups.getOrDefault(
                            nonNull(attrs.getValue("ref")), List.of())) {
                        handleCv(cv[0], cv[1]);
                    }
                    break;
                case "cvParam": {
                    String acc = attrs.getValue("accession");
                    if (acc == null) break;
                    String value = nonNull(attrs.getValue("value"));
                    if (currentGroup != null) currentGroup.add(new String[] {acc, value});
                    handleCv(acc, value);
                    break;
                }
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
                case "referenceableParamGroup":
                    currentGroup = null;
                    break;
                default:
                    break;
            }
        }

        private static String nonNull(String s) {
            return (s == null) ? "" : s;
        }

        private static String localName(String localName, String qName) {
            if (localName != null && !localName.isEmpty()) return localName;
            int colon = qName.indexOf(':');
            return (colon >= 0) ? qName.substring(colon + 1) : qName;
        }

        private void handleCv(String acc, String value) {
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
