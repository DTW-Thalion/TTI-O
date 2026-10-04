/*
 * TTI-O Java Implementation
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.transport;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums;
import global.thalion.ttio.Enums.IRMode;
import global.thalion.ttio.IRImage;
import global.thalion.ttio.Identification;
import global.thalion.ttio.InstrumentConfig;
import global.thalion.ttio.MSImage;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.Quantification;
import global.thalion.ttio.RamanImage;
import global.thalion.ttio.Sample;
import global.thalion.ttio.Subject;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.Spectrum;
import global.thalion.ttio.MassSpectrum;
import global.thalion.ttio.SpectrumIndex;
import global.thalion.ttio.codecs.BasePack;
import global.thalion.ttio.codecs.Rans;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.ReferenceImport;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.io.IOException;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Serialises an {@link SpectralDataset} as a transport byte stream.
 *
 * <p>Cross-language equivalents: Python
 * {@code ttio.transport.codec.TransportWriter}, Objective-C
 * {@code TTIOTransportWriter}.</p>
 */
public final class TransportWriter implements AutoCloseable {

    private final OutputStream out;
    private final boolean ownsStream;
    private boolean useChecksum = false;
    private boolean useCompression = false;
    /** Default AU channel codec when compression is on: codec 17,
     *  one FDZ1 stream per channel (wire id 17). */
    public static final String DEFAULT_COMPRESSION_CODEC = "float_delta_zstd";
    /** AU channel codec when {@link #useCompression} is on:
     *  "float_delta_zstd" (wire id 17, default), "zstd" (wire id 16,
     *  level 3) or "zlib" (wire id 1). */
    private String compressionCodec = DEFAULT_COMPRESSION_CODEC;
    /** Phase 2c-T: when true, probe each genomic run for v2 codec
     *  blobs and emit BlobV2* packets carrying them verbatim. */
    private boolean useBulkMode = false;

    public TransportWriter(OutputStream out) {
        this.out = out;
        this.ownsStream = false;
    }

    public TransportWriter(Path path) throws IOException {
        this.out = new FileOutputStream(path.toFile());
        this.ownsStream = true;
    }

    public void setUseChecksum(boolean v) { this.useChecksum = v; }
    public void setUseCompression(boolean v) { this.useCompression = v; }
    public void setUseBulkMode(boolean v) { this.useBulkMode = v; }
    public boolean useCompression() { return useCompression; }
    public boolean useBulkMode() { return useBulkMode; }

    /** Select the AU channel codec used when compression is on:
     *  "float_delta_zstd" (default), "zstd" or "zlib". Readers older
     *  than a codec's addition reject its id (id 16: pre-1.8.0; id 17:
     *  1.8.0 and older), so name "zlib" explicitly until a
     *  deployment's readers are current. */
    public void setCompressionCodec(String codec) {
        if (!"float_delta_zstd".equals(codec) && !"zstd".equals(codec)
                && !"zlib".equals(codec)) {
            throw new IllegalArgumentException(
                "unsupported compressionCodec '" + codec
                + "'; expected \"float_delta_zstd\", \"zstd\" or \"zlib\"");
        }
        this.compressionCodec = codec;
    }
    public String compressionCodec() { return compressionCodec; }

    @Override
    public void close() throws IOException {
        if (ownsStream) out.close();
    }

    // ---------------------------------------------------------- primitives

    private static long nowNs() { return System.currentTimeMillis() * 1_000_000L; }

    private void emit(PacketType type, byte[] payload, int datasetId, long auSequence)
            throws IOException {
        int flags = useChecksum ? PacketHeader.FLAG_HAS_CHECKSUM : 0;
        emitRawPacket(type, flags, datasetId, auSequence, payload);
    }

    /** Emit a packet with arbitrary flag bits. Used by
     *  {@link global.thalion.ttio.protection.EncryptedTransport} so it
     *  can set FLAG_ENCRYPTED / FLAG_ENCRYPTED_HEADER on AUs.
     *
     *
     */
    public void emitRawPacket(PacketType type, int flags, int datasetId,
                                long auSequence, byte[] payload) throws IOException {
        int finalFlags = useChecksum
            ? (flags | PacketHeader.FLAG_HAS_CHECKSUM)
            : flags;
        PacketHeader h = new PacketHeader(type, finalFlags, datasetId, auSequence,
                payload.length, nowNs());
        out.write(h.encode());
        out.write(payload);
        if (useChecksum) {
            int crc = Crc32c.compute(payload);
            ByteBuffer crcBuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
            crcBuf.putInt(crc);
            out.write(crcBuf.array());
        }
    }

    private static void appendLEString(ByteBuffer buf, String s, int width) {
        byte[] b = s == null
                ? new byte[0]
                : s.getBytes(StandardCharsets.UTF_8);
        if (width == 2) buf.putShort((short) (b.length & 0xFFFF));
        else            buf.putInt(b.length);
        buf.put(b);
    }

    private static int sizeLEString(String s, int width) {
        int b = s == null ? 0 : s.getBytes(StandardCharsets.UTF_8).length;
        return width + b;
    }

    // ---------------------------------------------------------- packets

    public void writeStreamHeader(String formatVersion, String title,
                                    String isaInvestigation, List<String> features,
                                    int nDatasets) throws IOException {
        int size = sizeLEString(formatVersion, 2)
                 + sizeLEString(title, 2)
                 + sizeLEString(isaInvestigation, 2)
                 + 2;
        for (String f : features) size += sizeLEString(f, 2);
        size += 2;
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        appendLEString(buf, formatVersion, 2);
        appendLEString(buf, title, 2);
        appendLEString(buf, isaInvestigation, 2);
        buf.putShort((short) (features.size() & 0xFFFF));
        for (String f : features) appendLEString(buf, f, 2);
        buf.putShort((short) (nDatasets & 0xFFFF));
        emit(PacketType.STREAM_HEADER, buf.array(), 0, 0);
    }

    public void writeDatasetHeader(int datasetId, String name, int acquisitionMode,
                                     String spectrumClass, List<String> channelNames,
                                     String instrumentJSON, long expectedAUCount)
            throws IOException {
        int size = 2
                 + sizeLEString(name, 2)
                 + 1
                 + sizeLEString(spectrumClass, 2)
                 + 1;
        for (String c : channelNames) size += sizeLEString(c, 2);
        size += sizeLEString(instrumentJSON, 4) + 4;
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (datasetId & 0xFFFF));
        appendLEString(buf, name, 2);
        buf.put((byte) (acquisitionMode & 0xFF));
        appendLEString(buf, spectrumClass, 2);
        buf.put((byte) (channelNames.size() & 0xFF));
        for (String c : channelNames) appendLEString(buf, c, 2);
        appendLEString(buf, instrumentJSON, 4);
        buf.putInt((int) (expectedAUCount & 0xFFFFFFFFL));
        emit(PacketType.DATASET_HEADER, buf.array(), datasetId, 0);
    }

    public void writeAccessUnit(int datasetId, long auSequence, AccessUnit au)
            throws IOException {
        emit(PacketType.ACCESS_UNIT, au.encode(), datasetId, auSequence);
    }

    /** Phase 2c-T (transport-spec §4.10). */
    public void writeBlobV2MateInfo(int datasetId, List<String> chromNames,
                                      byte[] blob) throws IOException {
        int size = 5;  // dataset_id(2) + codec_id(1) + n_chrom_names(2)
        for (String n : chromNames) size += sizeLEString(n, 2);
        size += 4 + blob.length;
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (datasetId & 0xFFFF));
        buf.put((byte) PacketType.CODEC_ID_MATE_INLINE_V2);
        buf.putShort((short) (chromNames.size() & 0xFFFF));
        for (String n : chromNames) appendLEString(buf, n, 2);
        buf.putInt(blob.length);
        buf.put(blob);
        emit(PacketType.BLOB_V2_MATE_INFO, buf.array(), datasetId, 0);
    }

    /** Phase 2c-T (transport-spec §4.11). */
    public void writeBlobV2RefDiff(int datasetId, String referenceUri,
                                     byte[] blob) throws IOException {
        int size = 3 + sizeLEString(referenceUri, 2) + 4 + blob.length;
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (datasetId & 0xFFFF));
        buf.put((byte) PacketType.CODEC_ID_REF_DIFF_V2);
        appendLEString(buf, referenceUri, 2);
        buf.putInt(blob.length);
        buf.put(blob);
        emit(PacketType.BLOB_V2_REF_DIFF, buf.array(), datasetId, 0);
    }

    /** Phase 2c-T (transport-spec §4.12). */
    public void writeBlobV2NameTok(int datasetId, byte[] blob) throws IOException {
        int size = 3 + 4 + blob.length;
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (datasetId & 0xFFFF));
        buf.put((byte) PacketType.CODEC_ID_NAME_TOKENIZED_V2);
        buf.putInt(blob.length);
        buf.put(blob);
        emit(PacketType.BLOB_V2_NAME_TOK, buf.array(), datasetId, 0);
    }

    // ----------------------------------------------- v0.11 §4.23

    /**
     * v0.11 Task 1.5: emit an {@code ENCRYPTION_ALGORITHM} (0x1B)
     * packet carrying the dataset-level {@code @encrypted} algorithm
     * name (e.g. {@code "aes-256-gcm"}). Wire layout matches
     * transport-spec §4.23:
     *
     * <pre>
     * algorithm_length:  uint16
     * algorithm_utf8:    bytes[algorithm_length]
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7.</p>
     *
     * <p>This packet only conveys the algorithm-name string; per-AU
     * key material continues to ride on {@code ProtectionMetadata}
     * (0x04).</p>
     */
    public void writeEncryptionAlgorithm(String algorithm) throws IOException {
        if (algorithm == null) {
            throw new IllegalArgumentException(
                "writeEncryptionAlgorithm: algorithm must not be null");
        }
        byte[] algoBytes = algorithm.getBytes(StandardCharsets.UTF_8);
        if (algoBytes.length > 0xFFFF) {
            throw new IOException(
                "ENCRYPTION_ALGORITHM: algorithm name " + algoBytes.length
                + " bytes exceeds uint16 max");
        }
        ByteBuffer buf = ByteBuffer.allocate(2 + algoBytes.length)
            .order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (algoBytes.length & 0xFFFF));
        buf.put(algoBytes);
        emit(PacketType.ENCRYPTION_ALGORITHM, buf.array(), 0, 0);
    }

    // ----------------------------------------------- v0.11 §4.21

    /**
     * v0.11 Task 1.6: emit a {@code DATASET_PROVENANCE} (0x18) packet
     * carrying the dataset-level provenance chain (format-spec §6.3).
     * One packet carries all records — see transport-spec §4.21 for
     * the wire layout:
     *
     * <pre>
     * record_count:        uint32
     * # repeated record_count times:
     * timestamp_unix:      int64
     * software_length:     uint16, software bytes[..]   (UTF-8)
     * parameters_length:   uint16, parameters_json[..]  (UTF-8 JSON)
     * input_refs_length:   uint16, input_refs_csv[..]   (UTF-8 CSV)
     * output_refs_length:  uint16, output_refs_csv[..]  (UTF-8 CSV)
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7. The
     * input_refs / output_refs lists ride as comma-joined UTF-8 — a
     * single empty string for an empty list (no separators).</p>
     *
     * <p>Distinct from the per-run {@code Provenance} (0x06) packet,
     * which carries one JSON record per packet.</p>
     */
    public void writeDatasetProvenance(List<ProvenanceRecord> records)
            throws IOException {
        if (records == null) {
            throw new IllegalArgumentException(
                "writeDatasetProvenance: records must not be null");
        }
        // Pre-compute UTF-8 byte arrays so we can size the buffer
        // exactly. Mirrors the StreamHeader/DatasetHeader emit pattern.
        byte[][] softwareBytes = new byte[records.size()][];
        byte[][] paramsBytes   = new byte[records.size()][];
        byte[][] inputsBytes   = new byte[records.size()][];
        byte[][] outputsBytes  = new byte[records.size()][];
        int size = 4;  // record_count
        for (int i = 0; i < records.size(); i++) {
            ProvenanceRecord r = records.get(i);
            softwareBytes[i] = nz(r.software())
                .getBytes(StandardCharsets.UTF_8);
            paramsBytes[i]   = r.parametersJson()
                .getBytes(StandardCharsets.UTF_8);
            inputsBytes[i]   = csvJoin(r.inputRefs())
                .getBytes(StandardCharsets.UTF_8);
            outputsBytes[i]  = csvJoin(r.outputRefs())
                .getBytes(StandardCharsets.UTF_8);
            for (byte[] b : new byte[][]{softwareBytes[i], paramsBytes[i],
                                          inputsBytes[i], outputsBytes[i]}) {
                if (b.length > 0xFFFF) {
                    throw new IOException(
                        "DATASET_PROVENANCE: per-field length " + b.length
                        + " exceeds uint16 max");
                }
            }
            size += 8                          // timestamp_unix
                  + 2 + softwareBytes[i].length
                  + 2 + paramsBytes[i].length
                  + 2 + inputsBytes[i].length
                  + 2 + outputsBytes[i].length;
        }
        ByteBuffer buf = ByteBuffer.allocate(size).order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(records.size());
        for (int i = 0; i < records.size(); i++) {
            ProvenanceRecord r = records.get(i);
            buf.putLong(r.timestampUnix());
            buf.putShort((short) (softwareBytes[i].length & 0xFFFF));
            buf.put(softwareBytes[i]);
            buf.putShort((short) (paramsBytes[i].length & 0xFFFF));
            buf.put(paramsBytes[i]);
            buf.putShort((short) (inputsBytes[i].length & 0xFFFF));
            buf.put(inputsBytes[i]);
            buf.putShort((short) (outputsBytes[i].length & 0xFFFF));
            buf.put(outputsBytes[i]);
        }
        emit(PacketType.DATASET_PROVENANCE, buf.array(), 0, 0);
    }

    private static String csvJoin(List<String> refs) {
        if (refs == null || refs.isEmpty()) return "";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < refs.size(); i++) {
            if (i > 0) sb.append(',');
            sb.append(refs.get(i));
        }
        return sb.toString();
    }

    // ----------------------------------------------- v0.11 §4.13–§4.15

    /**
     * Threshold below which a chromosome rides as raw UINT8
     * (encoding=0). Mirrors transport-spec §4.14: ZLIB framing costs
     * dominate short sequences, so the writer skips compression below
     * 4 KiB and lets the reader handle both encodings.
     */
    static final int REFERENCE_CHROMOSOME_ZLIB_THRESHOLD = 4096;

    /**
     * v0.11 Stage 1: emit a {@link ReferenceImport} as the packet
     * sequence
     * {@code REFERENCE_GROUP_HEADER (0x10) → N × REFERENCE_CHROMOSOME (0x11)
     *  → END_OF_REFERENCE_GROUP (0x12)}.
     *
     * <p>Wire layout matches transport-spec §4.13–§4.15. All
     * multi-byte integers are LITTLE-ENDIAN (spec §1.7). The
     * chromosome index rides in the packet header's
     * {@code auSequence} field (0-based). The MD5 hex string from
     * {@link ReferenceImport#md5Hex()} is emitted verbatim as 32
     * ASCII bytes.</p>
     *
     * <p>The encoding byte on each chromosome record is 0
     * (uncompressed UINT8) when the raw sequence is shorter than
     * {@link #REFERENCE_CHROMOSOME_ZLIB_THRESHOLD}, otherwise 1
     * (zlib via {@link java.util.zip.Deflater} with default
     * settings).</p>
     *
     * <p>Reader-side materialisation is added by Task 1.3; this
     * method only emits the wire bytes.</p>
     */
    public void writeReferenceGroup(ReferenceImport ref) throws IOException {
        List<String> chromNames = ref.chromosomes();
        List<byte[]> seqs = ref.sequences();
        int chromCount = chromNames.size();
        long totalBases = ref.totalBases();
        String md5Hex = ref.md5Hex();
        if (md5Hex == null || md5Hex.length() != 32) {
            throw new IOException(
                "ReferenceImport.md5Hex() must be 32 hex chars, got "
                + (md5Hex == null ? "null" : md5Hex.length()));
        }

        // -- REFERENCE_GROUP_HEADER (0x10) -------------------------------
        byte[] uriBytes = ref.uri().getBytes(StandardCharsets.UTF_8);
        byte[] md5HexBytes = md5Hex.getBytes(StandardCharsets.US_ASCII);
        int hdrSize = 2 + uriBytes.length + 4 + 8 + 32;
        ByteBuffer hbuf = ByteBuffer.allocate(hdrSize).order(ByteOrder.LITTLE_ENDIAN);
        hbuf.putShort((short) (uriBytes.length & 0xFFFF));
        hbuf.put(uriBytes);
        hbuf.putInt(chromCount);
        hbuf.putLong(totalBases);
        hbuf.put(md5HexBytes);
        emit(PacketType.REFERENCE_GROUP_HEADER, hbuf.array(), 0, 0);

        // -- REFERENCE_CHROMOSOME (0x11) — one per contig ----------------
        for (int i = 0; i < chromCount; i++) {
            String name = chromNames.get(i);
            byte[] seq = seqs.get(i);
            byte[] nameBytes = name.getBytes(StandardCharsets.UTF_8);
            byte encoding;
            byte[] payload;
            if (seq.length < REFERENCE_CHROMOSOME_ZLIB_THRESHOLD) {
                encoding = 0;
                payload = seq;
            } else {
                encoding = 1;
                payload = zlibDeflate(seq);
            }
            int recSize = 2 + nameBytes.length + 8 + 1 + 4 + payload.length;
            ByteBuffer cbuf = ByteBuffer.allocate(recSize).order(ByteOrder.LITTLE_ENDIAN);
            cbuf.putShort((short) (nameBytes.length & 0xFFFF));
            cbuf.put(nameBytes);
            cbuf.putLong((long) seq.length);
            cbuf.put(encoding);
            cbuf.putInt(payload.length);
            cbuf.put(payload);
            emit(PacketType.REFERENCE_CHROMOSOME, cbuf.array(), 0, i);
        }

        // -- END_OF_REFERENCE_GROUP (0x12) -------------------------------
        ByteBuffer ebuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ebuf.putInt(chromCount);
        emit(PacketType.END_OF_REFERENCE_GROUP, ebuf.array(), 0, 0);
    }

    // ----------------------------------------------- v0.11 §4.16-§4.18

    /**
     * v0.11 Task 1.7: emit an {@link MSImage} as the packet sequence
     * {@code IMAGE_HEADER (0x13) → N × IMAGE_PIXEL (0x14)
     *  → END_OF_IMAGE (0x15)}, where {@code N = width * height}.
     *
     * <p>Wire layout matches transport-spec §4.16-§4.18. All
     * multi-byte integers are LITTLE-ENDIAN (spec §1.7). Each pixel
     * rides as a continuous-mode IMAGE_PIXEL — the shared m/z axis
     * lives on the IMAGE_HEADER, and every pixel carries only its
     * intensities (FLOAT64, uncompressed). The pixel index rides in
     * the packet header's {@code auSequence} field
     * ({@code y * width + x}; 0-based).</p>
     *
     * <p>Processed-mode (per-pixel axis, signalled by
     * {@code is_continuous == 0}) is not yet emitted by this writer;
     * the corresponding decoder path in {@link TransportReader} is
     * also continuous-only at Task 1.7. Adding processed-mode is
     * deferred to a follow-up task once a Java MSImage value type
     * exposes per-pixel axes (the current {@link MSImage} model is
     * continuous by construction: one {@code mzAxis} shared by all
     * pixels).</p>
     *
     * <p>Reader-side materialisation is added by Task 1.7 in
     * {@link TransportReader}; this method only emits the wire
     * bytes.</p>
     */
    public void writeImage(MSImage img) throws IOException {
        if (img == null) {
            throw new IllegalArgumentException(
                "writeImage: image must not be null");
        }
        int width  = img.width();
        int height = img.height();
        int bins   = img.spectralPoints();
        double[] axis = img.mzAxis();  // length 0 or == bins
        byte axisKind = 0;  // mz
        byte scanPatternByte = scanPatternToByte(img.scanPattern());
        // Modality=0 (MS) has no modality-specific extras (spec §4.16).
        byte[] extras = new byte[0];

        // -- IMAGE_HEADER (0x13) -----------------------------------------
        emitImageHeader((byte) 0, width, height, bins,
            img.pixelSizeX(), img.pixelSizeY(),
            scanPatternByte, axisKind, axis != null ? axis : new double[0],
            (byte) 1,  // is_continuous = 1 (shared axis)
            img.title(), img.isaInvestigationId(),
            extras);

        // -- IMAGE_PIXEL (0x14) — one per pixel --------------------------
        // Continuous-mode payload: x(u32) + y(u32) + precision(u8) +
        // compression(u8) + payload_length(u32) + intensities[..].
        // We always emit FLOAT64 (precision=1) uncompressed (compression=0)
        // intensities mirroring the on-disk MSImage cube precision.
        byte precision = 1;       // FLOAT64
        byte compression = 0;     // NONE
        long pixelIndex = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int payloadLen = 8 * bins;
                int recSize = 4 + 4 + 1 + 1 + 4 + payloadLen;
                ByteBuffer pbuf = ByteBuffer.allocate(recSize)
                    .order(ByteOrder.LITTLE_ENDIAN);
                pbuf.putInt(x);
                pbuf.putInt(y);
                pbuf.put(precision);
                pbuf.put(compression);
                pbuf.putInt(payloadLen);
                // MSImage.spectrumAt uses (row, col) ordering — row=y, col=x.
                double[] spec = img.spectrumAt(y, x);
                for (int k = 0; k < bins; k++) pbuf.putDouble(spec[k]);
                emit(PacketType.IMAGE_PIXEL, pbuf.array(), 0, pixelIndex);
                pixelIndex++;
            }
        }

        // -- END_OF_IMAGE (0x15) -----------------------------------------
        // pixel_count_seen: uint32 per spec §4.18.
        ByteBuffer ebuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ebuf.putInt((int) (pixelIndex & 0xFFFFFFFFL));
        emit(PacketType.END_OF_IMAGE, ebuf.array(), 0, 0);
    }

    /**
     * v0.11 Task 5.3 (Deferral 1): emit a {@link RamanImage} as the
     * packet sequence {@code IMAGE_HEADER (0x13) → N × IMAGE_PIXEL
     * (0x14) → END_OF_IMAGE (0x15)} with {@code modality=1}.
     *
     * <p>Wire layout per transport-spec §4.16. The shared axis on
     * the IMAGE_HEADER carries the Raman wavenumbers vector
     * ({@code axis_kind = 1 = wavenumber}). The
     * {@code modality_extras} slot at the tail of the IMAGE_HEADER
     * carries the Raman-specific fields:</p>
     * <pre>
     *   excitation_wavelength_nm:  float64
     *   laser_power_mw:            float64
     * </pre>
     * <p>(16 bytes total.)</p>
     *
     * <p>Each pixel rides as a continuous-mode IMAGE_PIXEL whose
     * {@code payload_bytes} is a dense vector of
     * {@code spectrum_bins} FLOAT64 intensities at the shared
     * wavenumber axis.</p>
     */
    public void writeRamanImage(RamanImage img) throws IOException {
        if (img == null) {
            throw new IllegalArgumentException(
                "writeRamanImage: image must not be null");
        }
        int width  = img.width();
        int height = img.height();
        int bins   = img.spectralPoints();
        double[] axis = img.wavenumbers();
        byte axisKind = 1;  // wavenumber
        byte scanPatternByte = scanPatternToByte(img.scanPattern());
        // Raman modality_extras: 8B excitation_wavelength_nm + 8B laser_power_mw.
        ByteBuffer extrasBuf = ByteBuffer.allocate(16)
            .order(ByteOrder.LITTLE_ENDIAN);
        extrasBuf.putDouble(img.excitationWavelengthNm());
        extrasBuf.putDouble(img.laserPowerMw());

        emitImageHeader((byte) 1, width, height, bins,
            img.pixelSizeX(), img.pixelSizeY(),
            scanPatternByte, axisKind, axis != null ? axis : new double[0],
            (byte) 1,
            img.title(), img.isaInvestigationId(),
            extrasBuf.array());

        byte precision = 1;
        byte compression = 0;
        long pixelIndex = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int payloadLen = 8 * bins;
                int recSize = 4 + 4 + 1 + 1 + 4 + payloadLen;
                ByteBuffer pbuf = ByteBuffer.allocate(recSize)
                    .order(ByteOrder.LITTLE_ENDIAN);
                pbuf.putInt(x);
                pbuf.putInt(y);
                pbuf.put(precision);
                pbuf.put(compression);
                pbuf.putInt(payloadLen);
                double[] spec = img.spectrumAt(y, x);
                for (int k = 0; k < bins; k++) pbuf.putDouble(spec[k]);
                emit(PacketType.IMAGE_PIXEL, pbuf.array(), 0, pixelIndex);
                pixelIndex++;
            }
        }

        ByteBuffer ebuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ebuf.putInt((int) (pixelIndex & 0xFFFFFFFFL));
        emit(PacketType.END_OF_IMAGE, ebuf.array(), 0, 0);
    }

    /**
     * v0.11 Task 5.3 (Deferral 1): emit an {@link IRImage} as the
     * packet sequence {@code IMAGE_HEADER (0x13) → N × IMAGE_PIXEL
     * (0x14) → END_OF_IMAGE (0x15)} with {@code modality=2}.
     *
     * <p>Wire layout per transport-spec §4.16. The shared axis on
     * the IMAGE_HEADER carries the IR wavenumbers vector
     * ({@code axis_kind = 1 = wavenumber}). The
     * {@code modality_extras} slot at the tail of the IMAGE_HEADER
     * carries the IR-specific fields:</p>
     * <pre>
     *   ir_mode:            uint8   # 0=transmittance, 1=absorbance
     *   resolution_cm_inv:  float64
     * </pre>
     * <p>(9 bytes total.)</p>
     */
    public void writeIRImage(IRImage img) throws IOException {
        if (img == null) {
            throw new IllegalArgumentException(
                "writeIRImage: image must not be null");
        }
        int width  = img.width();
        int height = img.height();
        int bins   = img.spectralPoints();
        double[] axis = img.wavenumbers();
        byte axisKind = 1;  // wavenumber
        byte scanPatternByte = scanPatternToByte(img.scanPattern());
        // IR modality_extras: u8 ir_mode + f64 resolution_cm_inv = 9B.
        ByteBuffer extrasBuf = ByteBuffer.allocate(9)
            .order(ByteOrder.LITTLE_ENDIAN);
        extrasBuf.put((byte) (img.mode() == IRMode.ABSORBANCE ? 1 : 0));
        extrasBuf.putDouble(img.resolutionCmInv());

        emitImageHeader((byte) 2, width, height, bins,
            img.pixelSizeX(), img.pixelSizeY(),
            scanPatternByte, axisKind, axis != null ? axis : new double[0],
            (byte) 1,
            img.title(), img.isaInvestigationId(),
            extrasBuf.array());

        byte precision = 1;
        byte compression = 0;
        long pixelIndex = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int payloadLen = 8 * bins;
                int recSize = 4 + 4 + 1 + 1 + 4 + payloadLen;
                ByteBuffer pbuf = ByteBuffer.allocate(recSize)
                    .order(ByteOrder.LITTLE_ENDIAN);
                pbuf.putInt(x);
                pbuf.putInt(y);
                pbuf.put(precision);
                pbuf.put(compression);
                pbuf.putInt(payloadLen);
                double[] spec = img.spectrumAt(y, x);
                for (int k = 0; k < bins; k++) pbuf.putDouble(spec[k]);
                emit(PacketType.IMAGE_PIXEL, pbuf.array(), 0, pixelIndex);
                pixelIndex++;
            }
        }

        ByteBuffer ebuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ebuf.putInt((int) (pixelIndex & 0xFFFFFFFFL));
        emit(PacketType.END_OF_IMAGE, ebuf.array(), 0, 0);
    }

    /**
     * v0.11 Task 5.3: shared IMAGE_HEADER (0x13) packing routine.
     * Used by {@link #writeImage}, {@link #writeRamanImage},
     * {@link #writeIRImage} (and the processed-mode variant) so the
     * common shape stays byte-stable across modalities and the
     * {@code modality_extras} slot is appended once per call.
     */
    private void emitImageHeader(byte modality,
                                 int width, int height, int bins,
                                 double pixelSizeX, double pixelSizeY,
                                 byte scanPatternByte, byte axisKind,
                                 double[] axis,
                                 byte isContinuous,
                                 String title, String isaId,
                                 byte[] extras) throws IOException {
        if (axis == null) axis = new double[0];
        if (extras == null) extras = new byte[0];
        byte[] titleBytes = (title == null ? "" : title)
            .getBytes(StandardCharsets.UTF_8);
        byte[] isaBytes = (isaId == null ? "" : isaId)
            .getBytes(StandardCharsets.UTF_8);
        if (titleBytes.length > 0xFFFF) {
            throw new IOException("IMAGE_HEADER: title " + titleBytes.length
                + " bytes exceeds uint16 max");
        }
        if (isaBytes.length > 0xFFFF) {
            throw new IOException("IMAGE_HEADER: isa_id " + isaBytes.length
                + " bytes exceeds uint16 max");
        }
        if (extras.length > 0xFFFF) {
            throw new IOException("IMAGE_HEADER: modality_extras "
                + extras.length + " bytes exceeds uint16 max");
        }
        int axisLength = axis.length;
        int hdrSize = 1                  // modality
                    + 4                  // width
                    + 4                  // height
                    + 4                  // spectrum_bins
                    + 8                  // pixel_size_x
                    + 8                  // pixel_size_y
                    + 1                  // scan_pattern
                    + 1                  // axis_kind
                    + 4                  // axis_length
                    + 8 * axisLength     // axis values
                    + 1                  // is_continuous
                    + 2 + titleBytes.length
                    + 2 + isaBytes.length
                    + 2 + extras.length;  // modality_extras
        ByteBuffer hbuf = ByteBuffer.allocate(hdrSize).order(ByteOrder.LITTLE_ENDIAN);
        hbuf.put(modality);
        hbuf.putInt(width);
        hbuf.putInt(height);
        hbuf.putInt(bins);
        hbuf.putDouble(pixelSizeX);
        hbuf.putDouble(pixelSizeY);
        hbuf.put(scanPatternByte);
        hbuf.put(axisKind);
        hbuf.putInt(axisLength);
        for (int i = 0; i < axisLength; i++) hbuf.putDouble(axis[i]);
        hbuf.put(isContinuous);
        hbuf.putShort((short) (titleBytes.length & 0xFFFF));
        hbuf.put(titleBytes);
        hbuf.putShort((short) (isaBytes.length & 0xFFFF));
        hbuf.put(isaBytes);
        hbuf.putShort((short) (extras.length & 0xFFFF));
        hbuf.put(extras);
        emit(PacketType.IMAGE_HEADER, hbuf.array(), 0, 0);
    }

    /**
     * v0.11 Task 5.1 (Deferral 1): emit an {@link MSImage} as the
     * packet sequence {@code IMAGE_HEADER (0x13) → N × IMAGE_PIXEL
     * (0x14) → END_OF_IMAGE (0x15)} in <b>processed mode</b>
     * (sparse), where each pixel carries only its nonzero
     * {@code (channel_index, intensity)} pairs indexed into the
     * shared {@code mzAxis}. The dense cube is reconstructed by the
     * reader.
     *
     * <p>Wire layout per transport-spec §4.17 (LITTLE-ENDIAN). The
     * IMAGE_HEADER is identical to {@link #writeImage} except for
     * {@code is_continuous = 0}; each IMAGE_PIXEL payload is</p>
     * <pre>
     *   x(u32) + y(u32) + precision(u8) + compression(u8)
     *     + payload_length(u32)
     *     + payload_bytes = nonzero_count(u32)
     *         + nonzero_count × { channel_index(u32) + intensity(f64) }
     * </pre>
     *
     * <p>Nonzero is defined strictly as {@code v != 0.0}; NaN is
     * preserved (NaN compares unequal to 0.0, so it is emitted
     * verbatim as a nonzero entry). The MSImage data model stays
     * dense; processed mode is purely a wire optimisation for sparse
     * cubes.</p>
     *
     * <p>This is an opt-in sibling of {@link #writeImage}. Callers
     * pick continuous vs processed mode explicitly today; the
     * automatic heuristic (emit whichever is smaller) lands in a
     * follow-up task.</p>
     */
    public void writeImageProcessed(MSImage img) throws IOException {
        if (img == null) {
            throw new IllegalArgumentException(
                "writeImageProcessed: image must not be null");
        }
        int width  = img.width();
        int height = img.height();
        int bins   = img.spectralPoints();
        double[] axis = img.mzAxis();
        byte axisKind = 0;  // mz
        byte scanPatternByte = scanPatternToByte(img.scanPattern());
        // Modality=0 (MS) has no modality-specific extras (spec §4.16).
        byte[] extras = new byte[0];

        // -- IMAGE_HEADER (0x13) — shared layout, is_continuous=0 -------
        emitImageHeader((byte) 0, width, height, bins,
            img.pixelSizeX(), img.pixelSizeY(),
            scanPatternByte, axisKind, axis != null ? axis : new double[0],
            (byte) 0,  // is_continuous = 0 (processed)
            img.title(), img.isaInvestigationId(),
            extras);

        // -- IMAGE_PIXEL (0x14) — sparse per spec §4.17 -----------------
        // We always emit FLOAT64 (precision=1) uncompressed
        // (compression=0) intensities — mirrors writeImage above and
        // keeps the wire round-trip byte-exact with the cube.
        byte precision = 1;       // FLOAT64
        byte compression = 0;     // NONE
        long pixelIndex = 0;
        double[] cube = img.intensityCube();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int base = (y * width + x) * bins;
                // First pass: count nonzeros (v != 0.0; NaN counts as
                // nonzero because NaN != 0.0, preserving NaN on the wire).
                int nonzero = 0;
                for (int k = 0; k < bins; k++) {
                    if (cube[base + k] != 0.0) nonzero++;
                }
                int payloadLen = 4 + nonzero * (4 + 8);
                int recSize = 4 + 4 + 1 + 1 + 4 + payloadLen;
                ByteBuffer pbuf = ByteBuffer.allocate(recSize)
                    .order(ByteOrder.LITTLE_ENDIAN);
                pbuf.putInt(x);
                pbuf.putInt(y);
                pbuf.put(precision);
                pbuf.put(compression);
                pbuf.putInt(payloadLen);
                pbuf.putInt(nonzero);
                for (int k = 0; k < bins; k++) {
                    double v = cube[base + k];
                    if (v != 0.0) {
                        pbuf.putInt(k);
                        pbuf.putDouble(v);
                    }
                }
                emit(PacketType.IMAGE_PIXEL, pbuf.array(), 0, pixelIndex);
                pixelIndex++;
            }
        }

        // -- END_OF_IMAGE (0x15) -----------------------------------------
        ByteBuffer ebuf = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ebuf.putInt((int) (pixelIndex & 0xFFFFFFFFL));
        emit(PacketType.END_OF_IMAGE, ebuf.array(), 0, 0);
    }

    /** Map the MSImage scan-pattern string to the wire byte per spec
     *  §4.16 ({@code 0=flyback, 1=meander, 2=random}). The on-disk
     *  format uses "raster" as the default name for the flyback
     *  pattern. Unknown values map to 0 (flyback) defensively. */
    static byte scanPatternToByte(String scanPattern) {
        if (scanPattern == null) return 0;
        return switch (scanPattern) {
            case "raster", "flyback" -> 0;
            case "meander"           -> 1;
            case "random"            -> 2;
            default                  -> 0;
        };
    }

    // ----------------------------------------------- v0.11 §4.19 / §4.20

    /**
     * v0.11 Task 1.8: emit an {@code IDENTIFICATIONS_TABLE} (0x16)
     * packet carrying the full identifications table as a single
     * length-prefixed Apache Arrow IPC stream. Wire layout per
     * transport-spec §4.19:
     *
     * <pre>
     * arrow_ipc_length:    uint32
     * arrow_ipc:           bytes[arrow_ipc_length]   # self-describing IPC stream
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7. The Arrow
     * IPC stream carries its own schema, row count, and null bitmaps,
     * so no per-row TLV envelope is needed. Empty lists MUST NOT
     * invoke this method — callers (see {@link #writeDataset}) gate on
     * {@code !list.isEmpty()} per §5.4 step 6 ("zero or more").</p>
     */
    public void writeIdentifications(List<Identification> identifications)
            throws IOException {
        if (identifications == null) {
            throw new IllegalArgumentException(
                "writeIdentifications: identifications must not be null");
        }
        byte[] ipc = ArrowIpcCodec.encodeIdentifications(identifications);
        ByteBuffer buf = ByteBuffer.allocate(4 + ipc.length)
            .order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(ipc.length);
        buf.put(ipc);
        emit(PacketType.IDENTIFICATIONS_TABLE, buf.array(), 0, 0);
    }

    /**
     * v0.11 Task 1.8: emit a {@code QUANTIFICATIONS_TABLE} (0x17)
     * packet carrying the full quantifications table as a single
     * length-prefixed Apache Arrow IPC stream. Wire layout per
     * transport-spec §4.20 — identical shape to §4.19 but with a
     * distinct packet type so receivers can dispatch without parsing
     * the IPC payload first.
     *
     * <pre>
     * arrow_ipc_length:    uint32
     * arrow_ipc:           bytes[arrow_ipc_length]
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7. Empty
     * lists MUST NOT invoke this method (§5.4 step 6).</p>
     */
    public void writeQuantifications(List<Quantification> quantifications)
            throws IOException {
        if (quantifications == null) {
            throw new IllegalArgumentException(
                "writeQuantifications: quantifications must not be null");
        }
        byte[] ipc = ArrowIpcCodec.encodeQuantifications(quantifications);
        ByteBuffer buf = ByteBuffer.allocate(4 + ipc.length)
            .order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(ipc.length);
        buf.put(ipc);
        emit(PacketType.QUANTIFICATIONS_TABLE, buf.array(), 0, 0);
    }

    // ----------------------------------------------- v0.11 §4.22 (Stage 6)

    /**
     * Stage 6 / Task 6.2: emit a {@code SUBJECT_METADATA} (0x19) packet
     * carrying the full {@link Subject} table as a single length-
     * prefixed Apache Arrow IPC stream. Wire layout per transport-spec
     * §4.22:
     *
     * <pre>
     * arrow_ipc_length:    uint32
     * arrow_ipc:           bytes[arrow_ipc_length]   # self-describing IPC stream
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7. The
     * payload schema and null-handling conventions live in
     * {@link ArrowIpcCodec}. Empty lists MUST NOT invoke this method —
     * callers (see {@link #writeDataset}) gate on
     * {@code !list.isEmpty()} per §5.4 step 5 ("zero or more").</p>
     */
    public void writeSubjectMetadata(List<Subject> subjects)
            throws IOException {
        if (subjects == null) {
            throw new IllegalArgumentException(
                "writeSubjectMetadata: subjects must not be null");
        }
        byte[] ipc = ArrowIpcCodec.encodeSubjects(subjects);
        ByteBuffer buf = ByteBuffer.allocate(4 + ipc.length)
            .order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(ipc.length);
        buf.put(ipc);
        emit(PacketType.SUBJECT_METADATA, buf.array(), 0, 0);
    }

    /**
     * Stage 6 / Task 6.2: emit a {@code SAMPLE_METADATA} (0x1A) packet
     * carrying the full {@link Sample} table as a single length-prefixed
     * Apache Arrow IPC stream. Wire layout per transport-spec §4.22 —
     * identical shape to the SUBJECT_METADATA framing with a distinct
     * packet type so receivers can dispatch without parsing the IPC
     * payload first.
     *
     * <pre>
     * arrow_ipc_length:    uint32
     * arrow_ipc:           bytes[arrow_ipc_length]
     * </pre>
     *
     * <p>All multi-byte integers LITTLE-ENDIAN per spec §1.7. Empty
     * lists MUST NOT invoke this method (§5.4 step 5).</p>
     */
    public void writeSampleMetadata(List<Sample> samples)
            throws IOException {
        if (samples == null) {
            throw new IllegalArgumentException(
                "writeSampleMetadata: samples must not be null");
        }
        byte[] ipc = ArrowIpcCodec.encodeSamples(samples);
        ByteBuffer buf = ByteBuffer.allocate(4 + ipc.length)
            .order(ByteOrder.LITTLE_ENDIAN);
        buf.putInt(ipc.length);
        buf.put(ipc);
        emit(PacketType.SAMPLE_METADATA, buf.array(), 0, 0);
    }

    public void writeEndOfDataset(int datasetId, long finalAUSequence) throws IOException {
        ByteBuffer buf = ByteBuffer.allocate(6).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) (datasetId & 0xFFFF));
        buf.putInt((int) (finalAUSequence & 0xFFFFFFFFL));
        emit(PacketType.END_OF_DATASET, buf.array(), datasetId, 0);
    }

    public void writeEndOfStream() throws IOException {
        emit(PacketType.END_OF_STREAM, new byte[0], 0, 0);
    }

    // ---------------------------------------------------------- high-level

    public void writeDataset(SpectralDataset dataset) throws IOException {
        Map<String, AcquisitionRun> runs = dataset.msRuns();
        Map<String, GenomicRun> genomicRuns = dataset.genomicRuns();
        List<String> features = new ArrayList<>();
        for (String f : dataset.featureFlags().features()) features.add(f);

        // Phase 2c-T: declare bulk-mode in StreamHeader features when
        // the writer is bulk-enabled AND there is at least one
        // genomic run with v2 blobs to carry.
        boolean bulkActive = false;
        if (useBulkMode) {
            for (GenomicRun g : genomicRuns.values()) {
                if (bulkCarriable(g)) { bulkActive = true; break; }
            }
        }
        if (bulkActive && !features.contains(PacketType.BULK_MODE_V2_BLOBS_FEATURE)) {
            features.add(PacketType.BULK_MODE_V2_BLOBS_FEATURE);
        }

        // Task 1.4/1.5/1.6/1.7/1.8/6.2: detect v0.11 content
        // (references + encryption algorithm + dataset provenance +
        // image cube + identifications/quantifications tables +
        // subjects + samples — all the v0.11 prelude content types
        // ship under the same TRANSPORT_V0_11_FEATURE flag).
        // JIT2: capture the (cached) typed image handles once and reuse
        // them for both the v0.11 gating check and emission below — the
        // accessor reads a cached field, so this preserves prior behavior.
        MSImage msImage = (MSImage) dataset.imageForKind(Enums.ImageKind.MS);
        RamanImage ramanImage = (RamanImage) dataset.imageForKind(Enums.ImageKind.RAMAN);
        IRImage irImage = (IRImage) dataset.imageForKind(Enums.ImageKind.IR);
        boolean v011 = !dataset.references().isEmpty()
                    || dataset.isEncrypted()
                    || !dataset.provenanceRecords().isEmpty()
                    || msImage != null
                    || ramanImage != null
                    || irImage != null
                    || !dataset.identifications().isEmpty()
                    || !dataset.quantifications().isEmpty()
                    || !dataset.subjects().isEmpty()
                    || !dataset.samples().isEmpty();
        if (v011 && !features.contains(PacketType.TRANSPORT_V0_11_FEATURE)) {
            features.add(PacketType.TRANSPORT_V0_11_FEATURE);
        }

        writeStreamHeader("1.2", dataset.title(), dataset.isaInvestigationId(),
                features, runs.size() + genomicRuns.size());

        // Task 1.4/1.5/1.6/1.7/1.8/6.2: v0.11 prelude -- per §5.4
        // ordering, v0.11 sections come BEFORE the v0.10 dataset/run
        // sections, and the sub-sections appear in this order:
        //   §5.4.1 ENCRYPTION_ALGORITHM
        //   §5.4.2 DATASET_PROVENANCE
        //   §5.4.3 SUBJECT_METADATA / SAMPLE_METADATA  (subjects first)
        //   §5.4.4 reference groups
        //   §5.4.5 image cubes
        //   §5.4.6 IDENTIFICATIONS_TABLE / QUANTIFICATIONS_TABLE
        if (v011) {
            if (dataset.isEncrypted()) {
                writeEncryptionAlgorithm(dataset.encryptedAlgorithm());
            }
            if (!dataset.provenanceRecords().isEmpty()) {
                writeDatasetProvenance(dataset.provenanceRecords());
            }
            // §5.4.3 Stage 6 / Task 6.2: SUBJECT_METADATA (0x19) emits
            // before SAMPLE_METADATA (0x1A) so a downstream reader sees
            // subjects ahead of any samples that soft-FK into them.
            // Empty lists emit NO packet (spec §5.4 step 5 says "zero
            // or more").
            if (!dataset.subjects().isEmpty()) {
                writeSubjectMetadata(dataset.subjects());
            }
            if (!dataset.samples().isEmpty()) {
                writeSampleMetadata(dataset.samples());
            }
            for (ReferenceImport ref : dataset.references().values()) {
                writeReferenceGroup(ref);
            }
            // Task 5.3 (Deferral 1): image cubes emit in MS → Raman → IR
            // order so a downstream reader sees a deterministic sequence
            // when more than one modality is present on the same dataset.
            if (msImage != null) {
                writeImage(msImage);
            }
            if (ramanImage != null) {
                writeRamanImage(ramanImage);
            }
            if (irImage != null) {
                writeIRImage(irImage);
            }
            // §5.4 step 6: identifications first, then quantifications.
            // Empty lists emit NO packet (spec says "zero or more").
            if (!dataset.identifications().isEmpty()) {
                writeIdentifications(dataset.identifications());
            }
            if (!dataset.quantifications().isEmpty()) {
                writeQuantifications(dataset.quantifications());
            }
        }

        // Spectral dataset headers: ids 1..N.
        int id = 1;
        for (Map.Entry<String, AcquisitionRun> e : runs.entrySet()) {
            AcquisitionRun run = e.getValue();
            List<String> channelNames = new ArrayList<>(run.channels().keySet());
            String instrumentJSON = instrumentConfigJson(run.instrumentConfig());
            writeDatasetHeader(id, e.getKey(),
                    run.acquisitionMode().ordinal(),
                    run.spectrumClassName(),
                    channelNames,
                    instrumentJSON,
                    run.spectrumCount());
            id++;
        }

        // M89.2/M89.4: Genomic dataset headers: ids N+1..N+M (contiguous).
        // now lists 5 channels (sequences, qualities + the 3
        // per-AU compound strings cigar/read_name/mate_chromosome).
        for (Map.Entry<String, GenomicRun> e : genomicRuns.entrySet()) {
            GenomicRun grun = e.getValue();
            writeDatasetHeader(id, e.getKey(),
                    grun.acquisitionMode().ordinal(),
                    "TTIOGenomicRead",
                    genomicChannelNames(grun),
                    genomicRunMetadataJson(grun),
                    grun.readCount());
            id++;
        }

        // Spectral AUs first.
        id = 1;
        for (Map.Entry<String, AcquisitionRun> e : runs.entrySet()) {
            AcquisitionRun run = e.getValue();
            int n = run.spectrumCount();
            List<String> channelNames = new ArrayList<>(run.channels().keySet());
            for (int i = 0; i < n; i++) {
                AccessUnit au = spectrumToAccessUnit(run, i, channelNames,
                        useCompression, compressionCodec);
                writeAccessUnit(id, i, au);
            }
            writeEndOfDataset(id, n);
            id++;
        }

        // Then genomic AUs. Phase 2c-T: in bulk mode, emit
        // verbatim v2 codec blobs first, then per-AU AccessUnits.
        for (Map.Entry<String, GenomicRun> e : genomicRuns.entrySet()) {
            if (useBulkMode && bulkCarriable(e.getValue())) {
                emitGenomicRunV2Blobs(id, e.getValue());
            }
            emitGenomicRunAccessUnits(id, e.getValue());
            writeEndOfDataset(id, e.getValue().readCount());
            id++;
        }
        writeEndOfStream();
    }

    /** Whole-channel runs and one-block {@code blocks_v1} runs carry
     *  their v2 blobs verbatim; a multi-block run goes per AU, and so
     *  does a grouped run (M103): its blobs are in stored order while its
     *  access units are in input order. */
    static boolean bulkCarriable(GenomicRun run) {
        if ("blocks_v1_grouped".equals(run.layout())) return false;
        return !"blocks_v1".equals(run.layout()) || run.blockCount() == 1;
    }

    /** Phase 2c-T: probe a {@link GenomicRun} for verbatim v2
     *  codec blobs and emit the matching BlobV2* packets. Each
     *  blob is independently optional. */
    private void emitGenomicRunV2Blobs(int datasetId, GenomicRun run)
            throws IOException {
        byte[] mateBlob = run.readMateInfoInlineV2BlobBytes();
        if (mateBlob != null) {
            List<String> names = run.readMateInfoChromNamesTable();
            writeBlobV2MateInfo(datasetId, names, mateBlob);
        }
        byte[] nameBlob = run.readNameTokV2BlobBytes();
        if (nameBlob != null) {
            writeBlobV2NameTok(datasetId, nameBlob);
        }
        byte[] refDiffBlob = run.readRefDiffV2BlobBytes();
        if (refDiffBlob != null) {
            String refUri = run.referenceUri() == null ? "" : run.referenceUri();
            writeBlobV2RefDiff(datasetId, refUri, refDiffBlob);
        }
    }

    /** Write a single {@link GenomicRun} as a stream segment.
     *
     *  <p>Used by callers that drive emission manually (multiplexed
     *  streams, M89.4). The dataset header + AUs + end-of-dataset are
     *  emitted; the caller is responsible for stream framing
     *  (writeStreamHeader / writeEndOfStream).</p>
     *
     * (M89.2)
     */
    public void writeGenomicRun(int datasetId, String name, GenomicRun run)
            throws IOException {
        writeDatasetHeader(datasetId, name,
                run.acquisitionMode().ordinal(),
                "TTIOGenomicRead",
                genomicChannelNames(run),
                genomicRunMetadataJson(run),
                run.readCount());
        emitGenomicRunAccessUnits(datasetId, run);
        writeEndOfDataset(datasetId, run.readCount());
    }

    /** M89.2/M90.9: emit one ACCESS_UNIT packet per AlignedRead in
     *  {@code run}.
     *
     *  <p>M89.2: per-read fixed fields go into the AU's genomic suffix
     *  (chromosome / position / mapping_quality / flags). The
     *  variable-length sequences and qualities ride as two UINT8
     *  channels with the per-read slice as data.</p>
     *
     *  <p>M90.9: compound fields now also round-trip on the wire.
     *  cigar, read_name, mate_chromosome ride as additional UINT8
     *  string channels (one per AU). mate_position + template_length
     *  live in the M90.9 mate extension at the end of the AU genomic
     *  suffix.</p>
     *
     *  <p>M90.10: when the source channel carries an {@code @compression}
     *  attribute naming an M86 codec (RANS_ORDER0/1, BASE_PACK), the
     *  writer re-encodes each per-AU slice with the same codec on
     *  the wire. The wire ChannelData.compression byte tells the
     *  reader which decoder to dispatch. The 3 string channels
     *  (cigar / read_name / mate_chromosome) ALWAYS ride uncompressed —
     *  per-AU codec framing dominates short strings.</p>
     */
    private void emitGenomicRunAccessUnits(int datasetId, GenomicRun run)
            throws IOException {
        // Delegate to the package-private helper so DatasetWalker (#141)
        // can yield the same AUs without re-implementing per-read
        // construction. Both callers see byte-identical AUs.
        List<AccessUnit> aus = genomicRunAccessUnits(run);
        for (int i = 0; i < aus.size(); i++) {
            writeAccessUnit(datasetId, i, aus.get(i));
        }
    }

    /** The DatasetHeader channel list of a genomic run: the five M89
     *  channels, plus {@code tags} (transport-spec 4.3.1) when the run
     *  carries SAM tags (M101), so a tag-less run's stream is unchanged. */
    static List<String> genomicChannelNames(GenomicRun run) {
        List<String> names = new ArrayList<>(List.of("sequences", "qualities",
            "cigar", "read_name", "mate_chromosome"));
        if (run.hasTagsChannel()) names.add("tags");
        return names;
    }

    /** Build the per-read {@link AccessUnit} list for {@code run},
     *  matching the layout {@link #emitGenomicRunAccessUnits} writes
     *  to the wire. Extracted in #141 so {@link DatasetWalker} can
     *  emit identical AUs without duplicating the construction. */
    static List<AccessUnit> genomicRunAccessUnits(GenomicRun run) {
        int n = run.readCount();
        int precisionUint8 = Enums.Precision.UINT8.ordinal();
        int compressionNone = Enums.Compression.NONE.ordinal();
        int acqMode = run.acquisitionMode().ordinal() & 0xFF;
        // probe source @compression on sequences + qualities
        // so the wire codec mirrors the file's codec choice. The
        // string channels (cigar/read_name/mate_chromosome) always
        // ride uncompressed.
        int seqCodec = run.signalChannelCompressionCode("sequences");
        int qualCodec = run.signalChannelCompressionCode("qualities");
        // Bulk-fetch the byte channels + read-names list once. Mirrors
        // the Python encoder's pattern in transport/codec.py and the
        // FastqWriter bulk path. Without this pre-fetch, per-record
        // time is dominated by the objectAtIndex String roundtrip +
        // AlignedRead allocation; pre-fetching skips both for the seq
        // and name paths.
        // A grouped run (M103) is walked through its input-order read
        // iterator, which decodes each block once per pass; its per-read
        // accessors would hop between blocks.
        java.util.Iterator<global.thalion.ttio.genomics.AlignedRead> grouped =
            "blocks_v1_grouped".equals(run.layout()) ? run.iterReads() : null;
        byte[] seqAll = n > 0 && grouped == null ? run.sequencesFull() : new byte[0];
        byte[] qualAll = n > 0 && grouped == null ? run.qualitiesFull() : new byte[0];
        java.util.List<String> namesAll = grouped == null ? run.readNamesAll() : List.of();
        global.thalion.ttio.genomics.GenomicIndex idx = grouped == null ? run.index() : null;
        boolean hasTags = run.hasTagsChannel();
        List<AccessUnit> out = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            if (grouped != null) {
                out.add(genomicReadAccessUnit(grouped.next(), acqMode, seqCodec,
                                              qualCodec, hasTags));
                continue;
            }
            long offset = idx.offsetAt(i);
            int length = idx.lengthAt(i);
            byte[] seqBytes = new byte[length];
            if (length > 0) System.arraycopy(seqAll, (int) offset, seqBytes, 0, length);
            byte[] qualBytes;
            if (qualAll.length >= offset + length) {
                qualBytes = new byte[length];
                if (length > 0) System.arraycopy(qualAll, (int) offset, qualBytes, 0, length);
            } else {
                qualBytes = new byte[0];
            }
            // Bypass run.readAt(i) entirely: pull cigar / mateChrom /
            // matePos / tlen via the per-field cached accessors and
            // chrom/pos/mapq/flags from the index. All cached after
            // first call, so per-record cost is dominated by the
            // .getBytes() calls on the three string channels.
            String cigar = run.cigarAt(i);
            String name = i < namesAll.size() ? namesAll.get(i) : "";
            String mateChrom = run.mateChromAt(i);
            long matePos = run.matePosAt(i);
            int tlen = run.mateTlenAt(i);
            String chromosome = idx.chromosomeAt(i);
            long position = idx.positionAt(i);
            int mappingQuality = idx.mappingQualityAt(i);
            int flagsValue = idx.flagsAt(i);
            // re-encode per-AU slice with the M86 codec when
            // the source channel had an @compression attribute set.
            byte[] seqPayload = applyWireCodec(seqBytes, seqCodec);
            byte[] qualPayload = applyWireCodec(qualBytes, qualCodec);
            byte[] cigarBytes = (cigar == null ? "" : cigar)
                .getBytes(StandardCharsets.UTF_8);
            byte[] nameBytes = name.getBytes(StandardCharsets.UTF_8);
            byte[] mateChrBytes = (mateChrom == null ? "" : mateChrom)
                .getBytes(StandardCharsets.UTF_8);
            List<ChannelData> channels = new ArrayList<>(5);
            channels.add(new ChannelData("sequences", precisionUint8,
                    seqCodec, length, seqPayload));
            channels.add(new ChannelData("qualities", precisionUint8,
                    qualCodec, qualBytes.length, qualPayload));
            channels.add(new ChannelData("cigar", precisionUint8,
                    compressionNone, cigarBytes.length, cigarBytes));
            channels.add(new ChannelData("read_name", precisionUint8,
                    compressionNone, nameBytes.length, nameBytes));
            channels.add(new ChannelData("mate_chromosome", precisionUint8,
                    compressionNone, mateChrBytes.length, mateChrBytes));
            if (hasTags) {
                // M101: the read's SAM tag text, UTF-8, uncompressed.
                byte[] tagBytes = run.tagsAt(i).getBytes(StandardCharsets.UTF_8);
                channels.add(new ChannelData("tags", precisionUint8,
                        compressionNone, tagBytes.length, tagBytes));
            }
            AccessUnit au = new AccessUnit(
                    5,                  // spectrum_class GenomicRead
                    acqMode,
                    0,                  // ms_level
                    2,                  // polarity = unknown (wire)
                    0.0, 0.0, 0,        // rt, precursor_mz, precursor_charge
                    0.0, 0.0,           // ion_mobility, base_peak_intensity
                    channels,
                    0L, 0L, 0L,         // pixel_x/y/z (unused for class==5)
                    chromosome,
                    position,
                    mappingQuality,
                    flagsValue & 0xFFFF,
                    matePos,
                    tlen);
            out.add(au);
        }
        return out;
    }

    /** One read's AU, from the read itself: the same bytes
     *  {@link #genomicRunAccessUnits} builds from the bulk channels. */
    private static AccessUnit genomicReadAccessUnit(
            global.thalion.ttio.genomics.AlignedRead r, int acqMode,
            int seqCodec, int qualCodec, boolean hasTags) {
        int precisionUint8 = Enums.Precision.UINT8.ordinal();
        int compressionNone = Enums.Compression.NONE.ordinal();
        byte[] seqBytes = (r.sequence() == null ? "" : r.sequence())
            .getBytes(StandardCharsets.US_ASCII);
        byte[] qualBytes = r.qualities() == null ? new byte[0] : r.qualities();
        byte[] cigarBytes = (r.cigar() == null ? "" : r.cigar())
            .getBytes(StandardCharsets.UTF_8);
        byte[] nameBytes = (r.readName() == null ? "" : r.readName())
            .getBytes(StandardCharsets.UTF_8);
        byte[] mateChrBytes = (r.mateChromosome() == null ? "" : r.mateChromosome())
            .getBytes(StandardCharsets.UTF_8);
        List<ChannelData> channels = new ArrayList<>(5);
        channels.add(new ChannelData("sequences", precisionUint8,
                seqCodec, seqBytes.length, applyWireCodec(seqBytes, seqCodec)));
        channels.add(new ChannelData("qualities", precisionUint8,
                qualCodec, qualBytes.length, applyWireCodec(qualBytes, qualCodec)));
        channels.add(new ChannelData("cigar", precisionUint8,
                compressionNone, cigarBytes.length, cigarBytes));
        channels.add(new ChannelData("read_name", precisionUint8,
                compressionNone, nameBytes.length, nameBytes));
        channels.add(new ChannelData("mate_chromosome", precisionUint8,
                compressionNone, mateChrBytes.length, mateChrBytes));
        if (hasTags) {
            byte[] tagBytes = (r.tags() == null ? "" : r.tags()).getBytes(StandardCharsets.UTF_8);
            channels.add(new ChannelData("tags", precisionUint8,
                    compressionNone, tagBytes.length, tagBytes));
        }
        return new AccessUnit(
                5, acqMode, 0, 2,
                0.0, 0.0, 0,
                0.0, 0.0,
                channels,
                0L, 0L, 0L,
                r.chromosome(),
                r.position(),
                r.mappingQuality(),
                r.flags() & 0xFFFF,
                r.matePosition(),
                r.templateLength());
    }

    /** encode {@code plaintext} with the given wire codec id.
     *  NONE → identity. Other ids dispatch to the matching M86 codec.
     *  Mirrors Python {@code _apply_wire_codec}. */
    private static byte[] applyWireCodec(byte[] plaintext, int codecId) {
        if (codecId == 0) return plaintext;  // NONE
        if (codecId == Enums.Compression.RANS_ORDER0.ordinal()) {
            return Rans.encode(plaintext, 0);
        }
        if (codecId == Enums.Compression.RANS_ORDER1.ordinal()) {
            return Rans.encode(plaintext, 1);
        }
        if (codecId == Enums.Compression.BASE_PACK.ordinal()) {
            return BasePack.encode(plaintext);
        }
        throw new UnsupportedOperationException(
            "applyWireCodec: codec id " + codecId
            + " not supported for genomic UINT8");
    }

    /** Per-genomic-run metadata serialised into the
     *  {@code instrument_json} slot of the dataset header. Mirrors
     *  Python {@code _genomic_run_metadata_json}: JSON object with
     *  reference_uri, platform, sample_name, modality, sort_keys=true.
     */
    static String genomicRunMetadataJson(GenomicRun run) {
        StringBuilder sb = new StringBuilder(96);
        sb.append('{');
        appendJsonField(sb, "modality",      nz(run.modality()),     false);
        appendJsonField(sb, "platform",      nz(run.platform()),     true);
        // M97: always present, "" when the run has no read role, in
        // sorted-key position so the JSON stays byte-identical with
        // the Python and ObjC writers.
        appendJsonField(sb, "read_role",     nz(run.getReadRole()),  true);
        appendJsonField(sb, "reference_uri", nz(run.referenceUri()), true);
        appendJsonField(sb, "sample_name",   nz(run.sampleName()),   true);
        sb.append('}');
        return sb.toString();
    }

    private static String nz(String s) { return s == null ? "" : s; }

    private static void appendJsonField(StringBuilder sb, String key,
                                          String value, boolean needsComma) {
        if (needsComma) sb.append(", ");
        sb.append('"').append(key).append("\": \"");
        for (int i = 0; i < value.length(); i++) {
            char c = value.charAt(i);
            if (c == '"' || c == '\\') sb.append('\\').append(c);
            else sb.append(c);
        }
        sb.append('"');
    }

    static String instrumentConfigJson(InstrumentConfig cfg) {
        if (cfg == null) return "{}";
        // Minimal sorted-key JSON emitter to match Python / ObjC output
        // without pulling a full JSON library in.
        StringBuilder sb = new StringBuilder(128);
        sb.append('{');
        appendField(sb, "analyzer_type", cfg.analyzerType(), false);
        appendField(sb, "detector_type", cfg.detectorType(), true);
        appendField(sb, "manufacturer",  cfg.manufacturer(),  true);
        appendField(sb, "model",         cfg.model(),         true);
        appendField(sb, "serial_number", cfg.serialNumber(),  true);
        appendField(sb, "source_type",   cfg.sourceType(),    true);
        sb.append('}');
        return sb.toString();
    }

    private static void appendField(StringBuilder sb, String key, String value,
                                      boolean needsComma) {
        if (needsComma) sb.append(", ");
        sb.append('"').append(key).append("\": \"");
        String v = value == null ? "" : value;
        for (int i = 0; i < v.length(); i++) {
            char c = v.charAt(i);
            if (c == '"' || c == '\\') sb.append('\\').append(c);
            else sb.append(c);
        }
        sb.append('"');
    }

    static AccessUnit spectrumToAccessUnit(AcquisitionRun run, int i,
                                              List<String> channelNames) {
        return spectrumToAccessUnit(run, i, channelNames, false,
                DEFAULT_COMPRESSION_CODEC);
    }

    static AccessUnit spectrumToAccessUnit(AcquisitionRun run, int i,
                                              List<String> channelNames,
                                              boolean useCompression) {
        return spectrumToAccessUnit(run, i, channelNames, useCompression,
                DEFAULT_COMPRESSION_CODEC);
    }

    static AccessUnit spectrumToAccessUnit(AcquisitionRun run, int i,
                                              List<String> channelNames,
                                              boolean useCompression,
                                              String compressionCodec) {
        Spectrum sp = run.objectAtIndex(i);
        int wireClass = wireFromSpectrumClassName(run.spectrumClassName());
        int msLevel = 0;
        int polarityWire = 2;
        if (sp instanceof MassSpectrum ms) {
            msLevel = ms.msLevel();
            polarityWire = wireFromPolarity(ms.polarity());
        }

        double bpi = 0.0;
        SpectrumIndex idx = run.spectrumIndex();
        if (idx != null && i < idx.count()) {
            bpi = idx.basePeakIntensities()[i];
        }

        List<ChannelData> channels = new ArrayList<>();
        for (String cname : channelNames) {
            int off = (int) idx.offsetAt(i);
            int len = idx.lengthAt(i);
            double[] vals = run.channelRange(cname, off, len);
            if (vals == null) continue;
            byte[] payload;
            int compressionCode;
            if (useCompression && "float_delta_zstd".equals(compressionCodec)) {
                payload = global.thalion.ttio.codecs.FloatDeltaZstd.encode(vals);
                compressionCode = Enums.Compression.FLOAT_DELTA_ZSTD.ordinal();
            } else {
                byte[] raw = new byte[len * 8];
                ByteBuffer buf = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
                for (int k = 0; k < len; k++) buf.putDouble(vals[k]);
                payload = raw;
                compressionCode = Enums.Compression.NONE.ordinal();
                if (useCompression && "zstd".equals(compressionCodec)) {
                    payload = zstdCompress(raw);
                    compressionCode = Enums.Compression.ZSTD.ordinal();
                } else if (useCompression) {
                    payload = zlibDeflate(raw);
                    compressionCode = Enums.Compression.ZLIB.ordinal();
                }
            }
            channels.add(new ChannelData(cname,
                    Enums.Precision.FLOAT64.ordinal(),
                    compressionCode,
                    len, payload));
        }

        return new AccessUnit(
                wireClass,
                run.acquisitionMode().ordinal(),
                msLevel,
                polarityWire,
                sp.scanTimeSeconds(),
                sp.precursorMz(),
                sp.precursorCharge(),
                0.0,
                bpi,
                channels,
                0, 0, 0);
    }

    private static byte[] zstdCompress(byte[] input) {
        io.airlift.compress.zstd.ZstdCompressor c =
                new io.airlift.compress.zstd.ZstdCompressor();
        byte[] buf = new byte[c.maxCompressedLength(input.length)];
        int n = c.compress(input, 0, input.length, buf, 0, buf.length);
        return java.util.Arrays.copyOf(buf, n);
    }

    private static byte[] zlibDeflate(byte[] input) {
        java.util.zip.Deflater def = new java.util.zip.Deflater();
        def.setInput(input);
        def.finish();
        byte[] buf = new byte[Math.max(64, input.length)];
        java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
        while (!def.finished()) {
            int n = def.deflate(buf);
            out.write(buf, 0, n);
        }
        def.end();
        return out.toByteArray();
    }

    private static int wireFromSpectrumClassName(String name) {
        if (name == null) return 0;
        return switch (name) {
            case "TTIOMassSpectrum"       -> 0;
            case "TTIONMRSpectrum"        -> 1;
            case "TTIONMR2DSpectrum"      -> 2;
            case "TTIOFreeInductionDecay" -> 3;
            case "TTIOMSImagePixel"       -> 4;
            default -> 0;
        };
    }

    private static int wireFromPolarity(Enums.Polarity p) {
        return switch (p) {
            case POSITIVE -> 0;
            case NEGATIVE -> 1;
            case UNKNOWN  -> 2;
        };
    }
}
