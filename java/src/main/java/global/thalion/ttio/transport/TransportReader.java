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
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.Subject;
import global.thalion.ttio.SpectrumIndex;
import global.thalion.ttio.MiniJson;
import global.thalion.ttio.codecs.BasePack;
import global.thalion.ttio.codecs.Rans;
import global.thalion.ttio.genomics.BulkV2Blobs;
import global.thalion.ttio.genomics.ReferenceImport;
import global.thalion.ttio.genomics.WrittenGenomicRun;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Deserialises a transport byte stream into {@link PacketRecord} values
 * or materializes the stream into a {@link SpectralDataset}.
 *
 * <p>Cross-language equivalents: Python
 * {@code ttio.transport.codec.TransportReader}, Objective-C
 * {@code TTIOTransportReader}.</p>
 */
public final class TransportReader implements AutoCloseable {

    private static final Logger LOG =
        Logger.getLogger(TransportReader.class.getName());

    private final InputStream in;
    private final boolean ownsStream;

    // v0.11 Stage 1 / Task 1.3: per-stream accumulator state for the
    // REFERENCE_GROUP_HEADER → N × REFERENCE_CHROMOSOME → END_OF_REFERENCE_GROUP
    // packet sequence. Reset at the start of every materializeTo() call.
    // Reads of these fields outside materializeTo() are undefined; the
    // reader is single-use (AutoCloseable) by contract.
    private String currentRefUri;
    private final List<String> currentChromNames = new ArrayList<>();
    private final List<byte[]> currentChromSeqs  = new ArrayList<>();
    private final List<ReferenceImport> collectedRefs = new ArrayList<>();
    // v0.11 Task 1.5: dataset-level @encrypted algorithm string carried
    // by ENCRYPTION_ALGORITHM (0x1B) packets. null when no such packet
    // appeared in the stream. Reset at the start of every
    // materializeTo() call.
    private String collectedEncryptionAlgorithm;
    // v0.11 Task 1.6: dataset-level provenance chain decoded from
    // DATASET_PROVENANCE (0x18) packets. Reset at the start of every
    // materializeTo() call. Passed into SpectralDataset.create as the
    // provenanceRecords arg so the on-disk /study/provenance_json
    // attribute round-trips.
    private final List<ProvenanceRecord> collectedProvenance = new ArrayList<>();
    // v0.11 Task 1.7: image-cube accumulator state for the
    // IMAGE_HEADER (0x13) → N × IMAGE_PIXEL (0x14) → END_OF_IMAGE
    // (0x15) packet sequence. Reset at the start of every
    // materializeTo() call. `currentImageBuilder` is non-null between
    // IMAGE_HEADER and END_OF_IMAGE.
    private ImageBuilder currentImageBuilder;
    private MSImage collectedImage;
    // v0.11 Task 5.3 (Deferral 1): per-modality image accumulators.
    // The IMAGE_HEADER's modality byte selects which builder type is
    // primed by startImage(); appendPixel + finishImage forward to the
    // same builder. Unknown modalities prime `currentImageSkipping`
    // (forward-compat skip-unknown path) and no image is materialised.
    private RamanImageBuilder currentRamanBuilder;
    private RamanImage collectedRamanImage;
    private IRImageBuilder currentIrBuilder;
    private IRImage collectedIrImage;
    // True between IMAGE_HEADER (unknown modality) and END_OF_IMAGE.
    // appendPixel + finishImage skip work; finishImage clears the flag.
    private boolean currentImageSkipping;
    // v0.11 Task 1.8: identification/quantification rows decoded from
    // IDENTIFICATIONS_TABLE (0x16) / QUANTIFICATIONS_TABLE (0x17)
    // packets. Reset at the start of every materializeTo() call.
    // Passed into SpectralDataset.create as the identifications /
    // quantifications args so the on-disk
    // /study/identifications_json + /study/quantifications_json
    // attributes round-trip.
    private final List<Identification> collectedIdentifications = new ArrayList<>();
    private final List<Quantification> collectedQuantifications = new ArrayList<>();
    // Stage 6 / Task 6.2: subject + sample rows decoded from
    // SUBJECT_METADATA (0x19) / SAMPLE_METADATA (0x1A) packets. Reset
    // at the start of every materializeTo() call. Persisted to the
    // resulting HDF5 file's /study/subjects/ + /study/samples/ groups
    // before the close+reopen dance so SpectralDataset.open() surfaces
    // them on the returned handle.
    private final List<Subject> collectedSubjects = new ArrayList<>();
    private final List<Sample> collectedSamples = new ArrayList<>();

    public TransportReader(InputStream in) {
        this.in = in;
        this.ownsStream = false;
    }

    public TransportReader(Path path) throws IOException {
        this.in = new FileInputStream(path.toFile());
        this.ownsStream = true;
    }

    public TransportReader(byte[] data) {
        this.in = new ByteArrayInputStream(data);
        this.ownsStream = true;
    }

    @Override
    public void close() throws IOException {
        if (ownsStream) in.close();
    }

    // ---------------------------------------------------------- record

    /** One parsed packet as header + payload. */
    public static final class PacketRecord {
        public final PacketHeader header;
        public final byte[] payload;
        public PacketRecord(PacketHeader header, byte[] payload) {
            this.header = header;
            this.payload = payload;
        }
    }

    // ---------------------------------------------------------- iteration

    public List<PacketRecord> readAllPackets() throws IOException {
        List<PacketRecord> out = new ArrayList<>();
        while (true) {
            byte[] headerBytes = in.readNBytes(PacketHeader.HEADER_SIZE);
            if (headerBytes.length == 0) return out;
            if (headerBytes.length < PacketHeader.HEADER_SIZE) {
                throw new IOException("truncated header: " + headerBytes.length);
            }
            PacketHeader header = PacketHeader.decode(headerBytes);
            byte[] payload = in.readNBytes((int) header.payloadLength);
            if (payload.length != header.payloadLength) {
                throw new IOException("truncated payload: "
                        + payload.length + "/" + header.payloadLength);
            }
            if ((header.flags & PacketHeader.FLAG_HAS_CHECKSUM) != 0) {
                byte[] crcBytes = in.readNBytes(4);
                if (crcBytes.length != 4) throw new IOException("truncated CRC-32C");
                int expected = ByteBuffer.wrap(crcBytes).order(ByteOrder.LITTLE_ENDIAN).getInt();
                int actual = Crc32c.compute(payload);
                if (expected != actual) {
                    throw new IOException("CRC-32C mismatch on packet type 0x"
                            + Integer.toHexString(header.packetTypeByte())
                            + ": expected=" + Integer.toHexString(expected)
                            + ", got=" + Integer.toHexString(actual));
                }
            }
            // Forward-compat (transport-spec §6): tolerate unknown
            // packet types so v0.10 readers can ingest v0.11+ streams.
            // The header was length-prefixed so we already consumed
            // the payload (and CRC if present) correctly — just log
            // and continue iterating past it.
            if (header.packetType == null) {
                LOG.fine("skipping unknown packet type 0x"
                    + Integer.toHexString(header.packetTypeByte()));
            }
            out.add(new PacketRecord(header, payload));
            if (header.packetType == PacketType.END_OF_STREAM) return out;
        }
    }

    /** Test-only accessor mirroring {@link #readAllPackets} for the
     *  skip-unknown forward-compat tests. Package-private. */
    List<PacketRecord> recordsForTest() throws IOException {
        return readAllPackets();
    }

    // ---------------------------------------------------------- materialize

    public SpectralDataset materializeTo(String outputPath) throws IOException {
        // Reset Task 1.3 reference-group accumulator state so a single
        // reader instance can drive multiple materializeTo() calls
        // safely (defensive — the AutoCloseable contract limits us to
        // one in practice, but the reset keeps the invariants explicit).
        currentRefUri = null;
        currentChromNames.clear();
        currentChromSeqs.clear();
        collectedRefs.clear();
        collectedEncryptionAlgorithm = null;
        collectedProvenance.clear();
        currentImageBuilder = null;
        collectedImage = null;
        currentRamanBuilder = null;
        collectedRamanImage = null;
        currentIrBuilder = null;
        collectedIrImage = null;
        currentImageSkipping = false;
        collectedIdentifications.clear();
        collectedQuantifications.clear();
        collectedSubjects.clear();
        collectedSamples.clear();

        List<PacketRecord> packets = readAllPackets();
        String title = "";
        String isa = "";
        Map<Integer, DatasetMeta> datasetMetas = new LinkedHashMap<>();
        Map<Integer, RunAccumulator> runAccs = new LinkedHashMap<>();
        // parallel accumulator for genomic datasets.
        Map<Integer, GenomicAccumulator> genomicAccs = new LinkedHashMap<>();
        // Phase 2c-T: per-genomic-dataset_id buffers for verbatim
        // v2 blobs. Populated when BlobV2* packets arrive.
        Map<Integer, BulkV2BlobsBuilder> bulkBlobs = new LinkedHashMap<>();
        Map<Integer, Long> lastSeq = new LinkedHashMap<>();
        boolean sawStreamHeader = false;
        boolean bulkModeRequired = false;

        for (PacketRecord rec : packets) {
            PacketHeader h = rec.header;
            // Forward-compat: skip packets whose type byte wasn't a
            // known PacketType. readAllPackets already consumed the
            // bytes; we just ignore the record for materialization.
            if (h.packetType == null) continue;
            ByteBuffer buf = ByteBuffer.wrap(rec.payload).order(ByteOrder.LITTLE_ENDIAN);

            if (h.packetType == PacketType.STREAM_HEADER) {
                if (sawStreamHeader) continue;
                sawStreamHeader = true;
                readLEString(buf, 2); // format_version
                title = readLEString(buf, 2);
                isa = readLEString(buf, 2);
                int nFeatures = buf.getShort() & 0xFFFF;
                for (int i = 0; i < nFeatures; i++) {
                    String f = readLEString(buf, 2);
                    if (PacketType.BULK_MODE_V2_BLOBS_FEATURE.equals(f)) {
                        bulkModeRequired = true;
                    }
                }
                // n_datasets - we don't need it; the headers carry their own ids.
                continue;
            }
            if (!sawStreamHeader) {
                throw new IOException("first packet must be StreamHeader; got " + h.packetType);
            }

            if (h.packetType == PacketType.DATASET_HEADER) {
                int datasetId = buf.getShort() & 0xFFFF;
                String name = readLEString(buf, 2);
                int acqMode = buf.get() & 0xFF;
                String spectrumClass = readLEString(buf, 2);
                int nch = buf.get() & 0xFF;
                List<String> channelNames = new ArrayList<>(nch);
                for (int i = 0; i < nch; i++) channelNames.add(readLEString(buf, 2));
                String instrumentJson = readLEString(buf, 4);
                long expected = buf.getInt() & 0xFFFFFFFFL;
                datasetMetas.put(datasetId, new DatasetMeta(
                        datasetId, name, acqMode, spectrumClass, channelNames,
                        instrumentJson, expected));
                // genomic datasets get a parallel accumulator.
                if ("TTIOGenomicRead".equals(spectrumClass)) {
                    genomicAccs.put(datasetId, new GenomicAccumulator());
                } else {
                    runAccs.put(datasetId, new RunAccumulator(channelNames));
                }
                continue;
            }
            if (h.packetType == PacketType.ACCESS_UNIT) {
                DatasetMeta meta = datasetMetas.get(h.datasetId);
                if (meta == null) {
                    throw new IOException("AccessUnit before DatasetHeader for id " + h.datasetId);
                }
                Long prev = lastSeq.get(h.datasetId);
                if (prev != null && h.auSequence <= prev) {
                    throw new IOException("non-monotonic au_sequence in dataset "
                            + h.datasetId + ": prev=" + prev + ", got=" + h.auSequence);
                }
                lastSeq.put(h.datasetId, h.auSequence);
                AccessUnit au = AccessUnit.decode(rec.payload);
                if (genomicAccs.containsKey(h.datasetId)) {
                    genomicAccs.get(h.datasetId).ingest(au);
                } else {
                    runAccs.get(h.datasetId).ingest(au);
                }
                continue;
            }
            if (h.packetType == PacketType.BLOB_V2_MATE_INFO) {
                int dsId = buf.getShort() & 0xFFFF;
                int codecId = buf.get() & 0xFF;
                if (codecId != PacketType.CODEC_ID_MATE_INLINE_V2) {
                    throw new IOException("BlobV2MateInfo bad codec_id " + codecId);
                }
                if (dsId != h.datasetId) {
                    throw new IOException("BlobV2MateInfo dataset_id " + dsId
                        + " != header.datasetId " + h.datasetId);
                }
                int nNames = buf.getShort() & 0xFFFF;
                List<String> names = new ArrayList<>(nNames);
                for (int i = 0; i < nNames; i++) names.add(readLEString(buf, 2));
                int blobLength = buf.getInt();
                byte[] blob = new byte[blobLength];
                buf.get(blob);
                BulkV2BlobsBuilder slot = bulkBlobs.computeIfAbsent(
                    dsId, k -> new BulkV2BlobsBuilder());
                if (slot.mateBlob != null) {
                    throw new IOException("duplicate BlobV2MateInfo for dataset_id " + dsId);
                }
                slot.mateBlob = blob;
                slot.mateChromNames = names;
                continue;
            }
            if (h.packetType == PacketType.BLOB_V2_REF_DIFF) {
                int dsId = buf.getShort() & 0xFFFF;
                int codecId = buf.get() & 0xFF;
                if (codecId != PacketType.CODEC_ID_REF_DIFF_V2) {
                    throw new IOException("BlobV2RefDiff bad codec_id " + codecId);
                }
                if (dsId != h.datasetId) {
                    throw new IOException("BlobV2RefDiff dataset_id mismatch");
                }
                String refUri = readLEString(buf, 2);
                int blobLength = buf.getInt();
                byte[] blob = new byte[blobLength];
                buf.get(blob);
                BulkV2BlobsBuilder slot = bulkBlobs.computeIfAbsent(
                    dsId, k -> new BulkV2BlobsBuilder());
                if (slot.refDiffBlob != null) {
                    throw new IOException("duplicate BlobV2RefDiff for dataset_id " + dsId);
                }
                slot.refDiffBlob = blob;
                slot.refDiffReferenceUri = refUri;
                continue;
            }
            if (h.packetType == PacketType.BLOB_V2_NAME_TOK) {
                int dsId = buf.getShort() & 0xFFFF;
                int codecId = buf.get() & 0xFF;
                if (codecId != PacketType.CODEC_ID_NAME_TOKENIZED_V2) {
                    throw new IOException("BlobV2NameTok bad codec_id " + codecId);
                }
                if (dsId != h.datasetId) {
                    throw new IOException("BlobV2NameTok dataset_id mismatch");
                }
                int blobLength = buf.getInt();
                byte[] blob = new byte[blobLength];
                buf.get(blob);
                BulkV2BlobsBuilder slot = bulkBlobs.computeIfAbsent(
                    dsId, k -> new BulkV2BlobsBuilder());
                if (slot.nameTokBlob != null) {
                    throw new IOException("duplicate BlobV2NameTok for dataset_id " + dsId);
                }
                slot.nameTokBlob = blob;
                continue;
            }
            if (h.packetType == PacketType.ENCRYPTION_ALGORITHM) {
                decodeEncryptionAlgorithm(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.DATASET_PROVENANCE) {
                decodeDatasetProvenance(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.REFERENCE_GROUP_HEADER) {
                startReferenceGroup(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.REFERENCE_CHROMOSOME) {
                appendChromosome(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.END_OF_REFERENCE_GROUP) {
                finishReferenceGroup();
                continue;
            }
            if (h.packetType == PacketType.IMAGE_HEADER) {
                startImage(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.IMAGE_PIXEL) {
                appendPixel(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.END_OF_IMAGE) {
                finishImage(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.IDENTIFICATIONS_TABLE) {
                decodeIdentificationsTable(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.QUANTIFICATIONS_TABLE) {
                decodeQuantificationsTable(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.SUBJECT_METADATA) {
                decodeSubjectMetadata(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.SAMPLE_METADATA) {
                decodeSampleMetadata(rec.payload);
                continue;
            }
            if (h.packetType == PacketType.END_OF_DATASET) continue;
            if (h.packetType == PacketType.END_OF_STREAM) break;
            // Annotation / Provenance / Chromatogram / Protection: skipped in M67.
        }

        // Phase 2c-T: a stream that declared bulk_mode_v2_blobs but
        // shipped zero blob packets is malformed.
        if (bulkModeRequired && bulkBlobs.isEmpty()) {
            throw new IOException(
                "StreamHeader declared " + PacketType.BULK_MODE_V2_BLOBS_FEATURE
                + " but no BlobV2* packets were received");
        }

        List<AcquisitionRun> runs = new ArrayList<>();
        for (Map.Entry<Integer, DatasetMeta> e : datasetMetas.entrySet()) {
            DatasetMeta meta = e.getValue();
            if (genomicAccs.containsKey(e.getKey())) continue;
            RunAccumulator acc = runAccs.get(e.getKey());
            SpectrumIndex idx = acc.toSpectrumIndex();
            Map<String, double[]> channelMap = acc.toChannelMap();
            InstrumentConfig cfg = new InstrumentConfig("", "", "", "", "", "");
            Enums.AcquisitionMode acqMode =
                    Enums.AcquisitionMode.values()[
                            Math.min(meta.acquisitionMode,
                                    Enums.AcquisitionMode.values().length - 1)];
            runs.add(new AcquisitionRun(meta.name, acqMode, idx, cfg,
                    channelMap, List.of(), List.of(), "", 0.0));
        }

        // build WrittenGenomicRun for each genomic dataset.
        // Phase 2c-T: attach any verbatim v2 blobs collected for this
        // dataset_id so write_minimal skips the v2 codec encode for
        // those channels.
        List<WrittenGenomicRun> genomicRuns = new ArrayList<>();
        for (Map.Entry<Integer, GenomicAccumulator> e : genomicAccs.entrySet()) {
            DatasetMeta meta = datasetMetas.get(e.getKey());
            WrittenGenomicRun gr = e.getValue().toWrittenGenomicRun(meta);
            BulkV2BlobsBuilder slot = bulkBlobs.get(e.getKey());
            if (slot != null && slot.hasAny()) {
                gr = gr.withBulkV2Blobs(new BulkV2Blobs(
                    slot.mateBlob, slot.mateChromNames,
                    slot.nameTokBlob,
                    slot.refDiffBlob, slot.refDiffReferenceUri));
            }
            genomicRuns.add(gr);
        }

        // Create the file then re-open so the returned dataset's
        // genomic StorageGroup handles are live (the create() call
        // closes its read-side handles after writing - GenomicRun
        // would then fail to open signal_channels lazily).
        SpectralDataset created = SpectralDataset.create(
            outputPath, title, isa, runs, genomicRuns,
            // v0.11 Task 1.8: pass collected IDENTIFICATIONS_TABLE
            // and QUANTIFICATIONS_TABLE rows into create so the
            // resulting on-disk file carries the round-tripped
            // identifications + quantifications tables.
            new ArrayList<>(collectedIdentifications),
            new ArrayList<>(collectedQuantifications),
            // v0.11 Task 1.6: pass collected DATASET_PROVENANCE records
            // into create so the resulting on-disk file carries the
            // round-tripped provenance chain.
            new ArrayList<>(collectedProvenance),
            global.thalion.ttio.FeatureFlags.defaultCurrent());
        // v0.11 Stage 1 / Task 1.3: embed any reference groups decoded
        // from the stream's REFERENCE_* packets. ReferenceImport.write-
        // ToDataset(ds) requires an open writable provider, so this
        // must run before the genomic close+reopen dance below.
        if (!collectedRefs.isEmpty()) {
            for (ReferenceImport ref : collectedRefs) {
                ref.writeToDataset(created);
            }
        }
        // v0.11 Task 1.7: embed any image cube decoded from the
        // stream's IMAGE_* packets. MSImage.writeTo(StorageGroup)
        // requires the /study/ group of an open writable provider,
        // so this runs before the close+reopen dance below.
        // Task 5.3 (Deferral 1): each modality lives in its own
        // /study/{image_cube,raman_image_cube,ir_image_cube} sub-group,
        // so all three can coexist on the same dataset.
        if (collectedImage != null
                || collectedRamanImage != null
                || collectedIrImage != null) {
            try (var studyGrp =
                    created.provider().rootGroup().openGroup("study")) {
                if (collectedImage != null) {
                    collectedImage.writeTo(studyGrp);
                }
                if (collectedRamanImage != null) {
                    collectedRamanImage.writeTo(studyGrp);
                }
                if (collectedIrImage != null) {
                    collectedIrImage.writeTo(studyGrp);
                }
            }
        }
        // Stage 6 / Task 6.2: persist Subject / Sample rows decoded
        // from SUBJECT_METADATA (0x19) / SAMPLE_METADATA (0x1A) packets
        // as per-row groups under /study/subjects/ + /study/samples/.
        // SpectralDataset.create()'s public-overloads do not yet
        // expose a "runs + genomicRuns + subjects + samples" form, so
        // we layer them onto the open file via StorageGroup directly —
        // mirrors the encryption-algorithm + image layering pattern
        // immediately above and below. The close+reopen at the bottom
        // of this method then surfaces them via SpectralDataset.open's
        // readSubjects + readSamples on the returned handle. The HDF5
        // layout matches SpectralDataset.writeSubjectsViaProvider /
        // writeSamplesViaProvider exactly (per-row groups with typed
        // attributes; absent group when the list is empty).
        if (!collectedSubjects.isEmpty() || !collectedSamples.isEmpty()) {
            try (var studyGrp =
                    created.provider().rootGroup().openGroup("study")) {
                if (!collectedSubjects.isEmpty()) {
                    try (var subjectsGroup = studyGrp.createGroup("subjects")) {
                        for (Subject s : collectedSubjects) {
                            try (var row = subjectsGroup.createGroup(s.externalId())) {
                                row.setAttribute("external_id", s.externalId());
                                if (!s.project().isEmpty()) {
                                    row.setAttribute("project", s.project());
                                }
                                if (!s.sex().isEmpty()) {
                                    row.setAttribute("sex", s.sex());
                                }
                                row.setAttribute("birth_year",
                                    Long.valueOf(s.birthYear()));
                                row.setAttribute("attributes_json",
                                    s.attributesJson());
                            }
                        }
                    }
                }
                if (!collectedSamples.isEmpty()) {
                    try (var samplesGroup = studyGrp.createGroup("samples")) {
                        for (Sample s : collectedSamples) {
                            try (var row = samplesGroup.createGroup(s.sampleId())) {
                                row.setAttribute("sample_id", s.sampleId());
                                if (!s.subjectExternalId().isEmpty()) {
                                    row.setAttribute("subject_external_id",
                                        s.subjectExternalId());
                                }
                                if (!s.sampleKind().isEmpty()) {
                                    row.setAttribute("sample_kind",
                                        s.sampleKind());
                                }
                                row.setAttribute("collected_at",
                                    Long.valueOf(s.collectedAt()));
                                row.setAttribute("attributes_json",
                                    s.attributesJson());
                            }
                        }
                    }
                }
            }
        }
        // v0.11 Task 1.5: persist the dataset-level @encrypted root
        // attribute so the materialised file reports isEncrypted() ==
        // true on reopen. SpectralDataset caches encryptedAlgorithm in
        // an instance field at construction time, so we must close +
        // reopen to surface the value on the returned dataset (even
        // for non-genomic streams).
        if (collectedEncryptionAlgorithm != null) {
            created.provider().rootGroup()
                .setAttribute("encrypted", collectedEncryptionAlgorithm);
        }
        // v0.11 Task 1.7 + Stage 6 / Task 6.2: when an image, subject,
        // or sample was embedded after create(), the returned `created`
        // dataset still has image()/subjects()/samples() == null /
        // empty (SpectralDataset caches them at construction). Force a
        // close+reopen so callers see the round-tripped content on the
        // returned handle.
        if (genomicRuns.isEmpty()
                && collectedEncryptionAlgorithm == null
                && collectedImage == null
                && collectedRamanImage == null
                && collectedIrImage == null
                && collectedSubjects.isEmpty()
                && collectedSamples.isEmpty()) {
            return created;
        }
        created.close();
        return SpectralDataset.open(outputPath);
    }

    // ---------------------------------------------------------- v0.11 §4.21 / §4.23

    /** v0.11 Task 1.5: decode an ENCRYPTION_ALGORITHM (0x1B) payload
     *  and stash the algorithm string for application at
     *  materialize-time. Multiple 0x1B packets are tolerated — the
     *  last-write-wins (spec §5.4 says "zero or more"; in practice
     *  the writer emits exactly one). */
    private void decodeEncryptionAlgorithm(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int len = bb.getShort() & 0xFFFF;
        byte[] algoBytes = new byte[len];
        bb.get(algoBytes);
        collectedEncryptionAlgorithm =
            new String(algoBytes, StandardCharsets.UTF_8);
    }

    /** v0.11 Task 1.6: decode a DATASET_PROVENANCE (0x18) payload per
     *  transport-spec §4.21 and append each record to
     *  {@link #collectedProvenance}. Multiple 0x18 packets MAY appear
     *  in a stream (spec §5.4 says "zero or more"); each carries its
     *  own record_count + records and they accumulate in emission
     *  order. */
    private void decodeDatasetProvenance(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int recordCount = bb.getInt();
        for (int i = 0; i < recordCount; i++) {
            long timestamp = bb.getLong();
            String software   = readLEString(bb, 2);
            String paramsJson = readLEString(bb, 2);
            String inputsCsv  = readLEString(bb, 2);
            String outputsCsv = readLEString(bb, 2);
            Map<String, String> params = parseParametersJson(paramsJson);
            List<String> inputs  = parseCsv(inputsCsv);
            List<String> outputs = parseCsv(outputsCsv);
            collectedProvenance.add(new ProvenanceRecord(
                timestamp, software, params, inputs, outputs));
        }
    }

    /** Parse a {@code {"k":"v", ...}} string back into a Map<String,
     *  String>, matching the format emitted by
     *  {@link ProvenanceRecord#parametersJson()}. Empty / "{}" → empty
     *  map. */
    private static Map<String, String> parseParametersJson(String json) {
        if (json == null || json.isEmpty() || "{}".equals(json)) {
            return java.util.Map.of();
        }
        // Defer to the shared MiniJson parser so we tolerate any
        // future field-set extensions byte-for-byte the same way the
        // dataset open path does (see SpectralDataset.parseProvenance-
        // Json). The format is "{...}" with string values only —
        // anything else gets coerced via toString(). LinkedHashMap to
        // preserve insertion order across the round-trip.
        Map<String, String> out = new java.util.LinkedHashMap<>();
        try {
            Object parsed = global.thalion.ttio.MiniJson.parse(json);
            if (parsed instanceof Map<?, ?> mraw) {
                for (var e : mraw.entrySet()) {
                    Object k = e.getKey();
                    Object v = e.getValue();
                    if (k != null && v != null) {
                        out.put(k.toString(), v.toString());
                    }
                }
            }
        } catch (Exception ignore) {
            // best-effort: leave map empty on parse failure
        }
        return out;
    }

    /** Parse a comma-joined UTF-8 ref list. Empty string → empty list. */
    private static List<String> parseCsv(String csv) {
        if (csv == null || csv.isEmpty()) return List.of();
        // No quoting / escaping in v0.11 — URIs are URL-encoded so
        // they cannot themselves contain commas. Plain split.
        String[] parts = csv.split(",", -1);
        return java.util.Arrays.asList(parts);
    }

    // ---------------------------------------------------------- v0.11 §4.19 / §4.20

    /** v0.11 Task 1.8: decode an IDENTIFICATIONS_TABLE (0x16) payload
     *  per transport-spec §4.19 — uint32 IPC length, then that many
     *  bytes of an Apache Arrow IPC stream. The decoded rows append to
     *  {@link #collectedIdentifications} so multiple 0x16 packets
     *  accumulate in emission order (spec §5.4 step 6 says "zero or
     *  more"). */
    private void decodeIdentificationsTable(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        long ipcLenLong = bb.getInt() & 0xFFFFFFFFL;
        int ipcLen = (int) ipcLenLong;
        byte[] ipc = new byte[ipcLen];
        bb.get(ipc);
        collectedIdentifications.addAll(ArrowIpcCodec.decodeIdentifications(ipc));
    }

    /** v0.11 Task 1.8: decode a QUANTIFICATIONS_TABLE (0x17) payload
     *  per transport-spec §4.20 — identical wire shape to §4.19 with
     *  a distinct dispatch. Rows append to
     *  {@link #collectedQuantifications}. */
    private void decodeQuantificationsTable(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        long ipcLenLong = bb.getInt() & 0xFFFFFFFFL;
        int ipcLen = (int) ipcLenLong;
        byte[] ipc = new byte[ipcLen];
        bb.get(ipc);
        collectedQuantifications.addAll(ArrowIpcCodec.decodeQuantifications(ipc));
    }

    // ---------------------------------------------------------- v0.11 §4.22 (Stage 6 / Task 6.2)

    /** Stage 6 / Task 6.2: decode a SUBJECT_METADATA (0x19) payload per
     *  transport-spec §4.22 — uint32 IPC length, then that many bytes
     *  of an Apache Arrow IPC stream. Rows append to
     *  {@link #collectedSubjects} so multiple 0x19 packets accumulate
     *  in emission order (spec §5.4 step 5 says "zero or more"). */
    private void decodeSubjectMetadata(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        long ipcLenLong = bb.getInt() & 0xFFFFFFFFL;
        int ipcLen = (int) ipcLenLong;
        byte[] ipc = new byte[ipcLen];
        bb.get(ipc);
        collectedSubjects.addAll(ArrowIpcCodec.decodeSubjects(ipc));
    }

    /** Stage 6 / Task 6.2: decode a SAMPLE_METADATA (0x1A) payload per
     *  transport-spec §4.22 — identical wire shape to the
     *  SUBJECT_METADATA framing with a distinct dispatch. Rows append
     *  to {@link #collectedSamples}. */
    private void decodeSampleMetadata(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        long ipcLenLong = bb.getInt() & 0xFFFFFFFFL;
        int ipcLen = (int) ipcLenLong;
        byte[] ipc = new byte[ipcLen];
        bb.get(ipc);
        collectedSamples.addAll(ArrowIpcCodec.decodeSamples(ipc));
    }


    // ---------------------------------------------------------- reference-group helpers

    /** v0.11 Stage 1 / Task 1.3: decode a REFERENCE_GROUP_HEADER
     *  (0x10) payload and prime the per-group accumulator. The
     *  {@code chromosome_count} field is parsed for cross-check
     *  purposes; the actual count is enforced by
     *  {@link #finishReferenceGroup()} when the EOR packet arrives. */
    private void startReferenceGroup(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int uriLen = bb.getShort() & 0xFFFF;
        byte[] uri = new byte[uriLen];
        bb.get(uri);
        currentRefUri = new String(uri, StandardCharsets.UTF_8);
        // chromosome_count (uint32) — used as a hint; the actual count
        // is the size of currentChromNames after all REFERENCE_CHROMOSOME
        // packets land. Parsed here to advance the buffer position.
        bb.getInt();
        // total_bases (uint64) — informational; the actual total comes
        // from summing the per-chromosome seq lengths. Parsed to advance.
        bb.getLong();
        // md5_hex[32] — informational; ReferenceImport recomputes from
        // the chromosome bytes (case-preserving, sort-by-name). Parsed
        // to consume the trailing 32 bytes of the payload.
        byte[] md5Hex = new byte[32];
        bb.get(md5Hex);
        currentChromNames.clear();
        currentChromSeqs.clear();
    }

    /** v0.11 Stage 1 / Task 1.3: decode a REFERENCE_CHROMOSOME (0x11)
     *  payload and append it to the current group. */
    private void appendChromosome(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int nameLen = bb.getShort() & 0xFFFF;
        byte[] name = new byte[nameLen];
        bb.get(name);
        long seqLength = bb.getLong();
        int encoding = bb.get() & 0xFF;
        long dataLenLong = bb.getInt() & 0xFFFFFFFFL;
        int dataLen = (int) dataLenLong;
        byte[] raw = new byte[dataLen];
        bb.get(raw);
        byte[] seq;
        if (encoding == 0) {
            seq = raw;
        } else if (encoding == 1) {
            seq = zlibInflate(raw);
            if (seq.length != seqLength) {
                throw new IllegalStateException(
                    "REFERENCE_CHROMOSOME zlib payload inflated to "
                    + seq.length + " bytes; expected " + seqLength);
            }
        } else {
            throw new IllegalStateException(
                "unknown REFERENCE_CHROMOSOME encoding: " + encoding);
        }
        currentChromNames.add(new String(name, StandardCharsets.UTF_8));
        currentChromSeqs.add(seq);
    }

    /** v0.11 Stage 1 / Task 1.3: close out the current reference
     *  group on END_OF_REFERENCE_GROUP (0x12). Builds a {@link
     *  ReferenceImport} from the accumulated chromosomes and stages it
     *  for {@code writeToDataset} after the {@code SpectralDataset} is
     *  created. */
    private void finishReferenceGroup() {
        if (currentRefUri == null) {
            throw new IllegalStateException(
                "END_OF_REFERENCE_GROUP without prior REFERENCE_GROUP_HEADER");
        }
        collectedRefs.add(new ReferenceImport(
            currentRefUri,
            new ArrayList<>(currentChromNames),
            new ArrayList<>(currentChromSeqs)));
        currentRefUri = null;
        currentChromNames.clear();
        currentChromSeqs.clear();
    }

    // ---------------------------------------------------------- v0.11 §4.16-§4.18

    /** v0.11 Task 1.7 / 5.1 / 5.3 (Deferral 1): decode an
     *  IMAGE_HEADER (0x13) payload and prime the per-image
     *  accumulator. Wire layout matches transport-spec §4.16. Both
     *  continuous-mode ({@code is_continuous == 1}) and
     *  processed-mode ({@code is_continuous == 0}, sparse
     *  {channel,intensity} pairs indexed into the shared axis) are
     *  supported; the mode flag is cached on the builder and read by
     *  {@link #appendPixel}.
     *
     *  <p>Modality dispatch (Task 5.3): the IMAGE_HEADER ends with a
     *  {@code u16 modality_extras_length + extras[length]} slot
     *  carrying modality-specific fields. Layout per modality:</p>
     *  <ul>
     *    <li>{@code modality == 0} (MS): extras is empty.</li>
     *    <li>{@code modality == 1} (Raman):
     *        {@code f64 excitation_wavelength_nm + f64 laser_power_mw}
     *        (16 bytes).</li>
     *    <li>{@code modality == 2} (IR):
     *        {@code u8 ir_mode + f64 resolution_cm_inv}
     *        (9 bytes).</li>
     *  </ul>
     *  Unknown modalities are logged + skipped — the wire-format
     *  {@code modality_extras_length} field is self-describing so the
     *  reader still advances past the header, and a flag is set so
     *  the following IMAGE_PIXEL + END_OF_IMAGE packets are dropped
     *  without aborting the stream (forward compat per §4.16). */
    private void startImage(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int modality = bb.get() & 0xFF;
        int width    = bb.getInt();
        int height   = bb.getInt();
        int bins     = bb.getInt();
        double pxX   = bb.getDouble();
        double pxY   = bb.getDouble();
        int scanPat  = bb.get() & 0xFF;
        int axisKind = bb.get() & 0xFF;
        int axisLen  = bb.getInt();
        double[] axis = new double[axisLen];
        for (int i = 0; i < axisLen; i++) axis[i] = bb.getDouble();
        int isContinuous = bb.get() & 0xFF;
        int titleLen = bb.getShort() & 0xFFFF;
        byte[] titleBytes = new byte[titleLen];
        bb.get(titleBytes);
        int isaLen = bb.getShort() & 0xFFFF;
        byte[] isaBytes = new byte[isaLen];
        bb.get(isaBytes);
        int extrasLen = bb.getShort() & 0xFFFF;
        byte[] extras = new byte[extrasLen];
        bb.get(extras);

        if (isContinuous != 0 && isContinuous != 1) {
            throw new IllegalStateException(
                "IMAGE_HEADER: is_continuous must be 0 or 1; got "
                + isContinuous);
        }
        String title = new String(titleBytes, StandardCharsets.UTF_8);
        String isa   = new String(isaBytes,   StandardCharsets.UTF_8);
        String scanPattern = scanPatternFromByte(scanPat);

        if (modality == 0) {
            currentImageBuilder = new ImageBuilder(
                width, height, bins,
                pxX, pxY, scanPattern, axis,
                title, isa,
                isContinuous == 1);
        } else if (modality == 1) {
            // Raman extras: 16 bytes (f64 excitation + f64 laser power).
            if (extrasLen != 16) {
                throw new IllegalStateException(
                    "IMAGE_HEADER (modality=1, Raman) expects 16-byte "
                    + "modality_extras (excitation + laser_power); got "
                    + extrasLen);
            }
            ByteBuffer eb = ByteBuffer.wrap(extras).order(ByteOrder.LITTLE_ENDIAN);
            double exc = eb.getDouble();
            double laser = eb.getDouble();
            currentRamanBuilder = new RamanImageBuilder(
                width, height, bins,
                pxX, pxY, scanPattern, axis,
                title, isa,
                exc, laser,
                isContinuous == 1);
        } else if (modality == 2) {
            // IR extras: 9 bytes (u8 ir_mode + f64 resolution).
            if (extrasLen != 9) {
                throw new IllegalStateException(
                    "IMAGE_HEADER (modality=2, IR) expects 9-byte "
                    + "modality_extras (ir_mode + resolution); got "
                    + extrasLen);
            }
            ByteBuffer eb = ByteBuffer.wrap(extras).order(ByteOrder.LITTLE_ENDIAN);
            int irModeByte = eb.get() & 0xFF;
            double resolution = eb.getDouble();
            IRMode mode = (irModeByte == 1) ? IRMode.ABSORBANCE
                                            : IRMode.TRANSMITTANCE;
            currentIrBuilder = new IRImageBuilder(
                width, height, bins,
                pxX, pxY, scanPattern, axis,
                title, isa,
                mode, resolution,
                isContinuous == 1);
        } else {
            // Unknown modality — skip the entire image block.
            LOG.warning("IMAGE_HEADER: unknown modality=" + modality
                + "; skipping image block (extrasLen=" + extrasLen
                + ", width=" + width + ", height=" + height + ")");
            currentImageSkipping = true;
        }
    }

    /** v0.11 Task 1.7 / 5.1: decode an IMAGE_PIXEL (0x14) payload
     *  per transport-spec §4.17 and stash the intensities at the
     *  pixel's row/col slot. The wire shape inside
     *  {@code payload_bytes} branches on the cached
     *  {@code isContinuous} from the IMAGE_HEADER:
     *  <ul>
     *    <li>continuous: dense {@code spectrum_bins} intensities</li>
     *    <li>processed: {@code u32 nonzero_count} + that many
     *        {@code u32 channel_index + fXX intensity} pairs;
     *        unmentioned channels stay at 0.0</li>
     *  </ul>
     */
    private void appendPixel(byte[] payload) {
        // Task 5.3: unknown-modality stream — silently drop the pixel.
        if (currentImageSkipping) return;
        if (currentImageBuilder == null
                && currentRamanBuilder == null
                && currentIrBuilder == null) {
            throw new IllegalStateException(
                "IMAGE_PIXEL received before IMAGE_HEADER");
        }
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        int x = bb.getInt();
        int y = bb.getInt();
        int precision   = bb.get() & 0xFF;
        int compression = bb.get() & 0xFF;
        long payloadLenLong = bb.getInt() & 0xFFFFFFFFL;
        int payloadLen = (int) payloadLenLong;
        byte[] raw = new byte[payloadLen];
        bb.get(raw);
        // §4.17 pixel enum: 0=none, 1=zstd, 2=zlib. Writers emit 0
        // today; the inflate paths make the enum real.
        if (compression == 1) {
            io.airlift.compress.zstd.ZstdDecompressor zd =
                    new io.airlift.compress.zstd.ZstdDecompressor();
            long plain = io.airlift.compress.zstd.ZstdDecompressor
                    .getDecompressedSize(raw, 0, raw.length);
            if (plain < 0 || plain > (1 << 27)) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL zstd frame without a sane content size: "
                    + plain);
            }
            byte[] out = new byte[(int) plain];
            zd.decompress(raw, 0, raw.length, out, 0, out.length);
            raw = out;
        } else if (compression == 2) {
            raw = zlibInflate(raw);
        } else if (compression != 0) {
            throw new IllegalStateException(
                "IMAGE_PIXEL compression=" + compression
                + " not supported (0=none, 1=zstd, 2=zlib)");
        }
        if (precision != 0 && precision != 1) {
            throw new IllegalStateException(
                "IMAGE_PIXEL precision=" + precision
                + " not supported (expected 0=float32 or 1=float64)");
        }
        // Resolve the active builder (one of MS/Raman/IR) and its
        // wire-mode flag. Each builder exposes the same setPixel
        // signature + spectralPoints field so the parse logic below
        // is identical across modalities.
        boolean isContinuous;
        int specBins;
        if (currentImageBuilder != null) {
            isContinuous = currentImageBuilder.isContinuous;
            specBins = currentImageBuilder.spectralPoints;
        } else if (currentRamanBuilder != null) {
            isContinuous = currentRamanBuilder.isContinuous;
            specBins = currentRamanBuilder.spectralPoints;
        } else {
            isContinuous = currentIrBuilder.isContinuous;
            specBins = currentIrBuilder.spectralPoints;
        }
        ByteBuffer ibuf = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        double[] intensities;
        if (isContinuous) {
            int n;
            if (precision == 1) {
                n = payloadLen / 8;
                intensities = new double[n];
                for (int k = 0; k < n; k++) intensities[k] = ibuf.getDouble();
            } else {
                n = payloadLen / 4;
                intensities = new double[n];
                for (int k = 0; k < n; k++) intensities[k] = ibuf.getFloat();
            }
        } else {
            // Processed-mode (sparse): u32 nonzero_count + entries.
            int nonzero = ibuf.getInt();
            intensities = new double[specBins];
            for (int k = 0; k < nonzero; k++) {
                int ch = ibuf.getInt();
                double v = (precision == 1) ? ibuf.getDouble()
                                            : (double) ibuf.getFloat();
                if (ch < 0 || ch >= specBins) {
                    throw new IllegalStateException(
                        "IMAGE_PIXEL (processed) channel_index " + ch
                        + " out of range [0, " + specBins + ") at pixel ("
                        + x + ", " + y + ")");
                }
                intensities[ch] = v;
            }
        }
        if (currentImageBuilder != null) {
            currentImageBuilder.setPixel(x, y, intensities);
        } else if (currentRamanBuilder != null) {
            currentRamanBuilder.setPixel(x, y, intensities);
        } else {
            currentIrBuilder.setPixel(x, y, intensities);
        }
    }

    /** v0.11 Task 1.7: close out the current image cube on
     *  END_OF_IMAGE (0x15). Verifies the {@code pixel_count_seen}
     *  field matches the per-pixel ingest count and stages the
     *  built {@link MSImage} for write-out after
     *  {@link SpectralDataset#create} returns. */
    private void finishImage(byte[] payload) {
        // Task 5.3: drain a skipped-modality block. The
        // pixel_count_seen field is still consumed for stream
        // hygiene but not validated against any per-pixel count
        // (we never accumulated one).
        if (currentImageSkipping) {
            currentImageSkipping = false;
            return;
        }
        if (currentImageBuilder == null
                && currentRamanBuilder == null
                && currentIrBuilder == null) {
            throw new IllegalStateException(
                "END_OF_IMAGE without prior IMAGE_HEADER");
        }
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        long declared = bb.getInt() & 0xFFFFFFFFL;
        if (currentImageBuilder != null) {
            long actual = currentImageBuilder.pixelsSeen();
            long expected = (long) currentImageBuilder.width
                          * currentImageBuilder.height;
            if (declared != actual) {
                throw new IllegalStateException(
                    "END_OF_IMAGE pixel_count_seen mismatch: declared="
                    + declared + ", actual=" + actual + " (width*height="
                    + expected + ")");
            }
            if (actual != expected) {
                throw new IllegalStateException(
                    "END_OF_IMAGE pixel count " + actual
                    + " does not equal width*height=" + expected);
            }
            collectedImage = currentImageBuilder.build();
            currentImageBuilder = null;
        } else if (currentRamanBuilder != null) {
            long actual = currentRamanBuilder.pixelsSeen();
            long expected = (long) currentRamanBuilder.width
                          * currentRamanBuilder.height;
            if (declared != actual) {
                throw new IllegalStateException(
                    "END_OF_IMAGE (Raman) pixel_count_seen mismatch: "
                    + "declared=" + declared + ", actual=" + actual
                    + " (width*height=" + expected + ")");
            }
            if (actual != expected) {
                throw new IllegalStateException(
                    "END_OF_IMAGE (Raman) pixel count " + actual
                    + " does not equal width*height=" + expected);
            }
            collectedRamanImage = currentRamanBuilder.build();
            currentRamanBuilder = null;
        } else {
            long actual = currentIrBuilder.pixelsSeen();
            long expected = (long) currentIrBuilder.width
                          * currentIrBuilder.height;
            if (declared != actual) {
                throw new IllegalStateException(
                    "END_OF_IMAGE (IR) pixel_count_seen mismatch: "
                    + "declared=" + declared + ", actual=" + actual
                    + " (width*height=" + expected + ")");
            }
            if (actual != expected) {
                throw new IllegalStateException(
                    "END_OF_IMAGE (IR) pixel count " + actual
                    + " does not equal width*height=" + expected);
            }
            collectedIrImage = currentIrBuilder.build();
            currentIrBuilder = null;
        }
    }

    /** Inverse of {@link TransportWriter#scanPatternToByte}. */
    private static String scanPatternFromByte(int b) {
        return switch (b) {
            case 0 -> "raster";
            case 1 -> "meander";
            case 2 -> "random";
            default -> "raster";
        };
    }

    /** Per-image accumulator. Fills the intensity cube as pixels
     *  arrive, validated against width/height/spectrum_bins from the
     *  IMAGE_HEADER. */
    private static final class ImageBuilder {
        final int width;
        final int height;
        final int spectralPoints;
        final double pixelSizeX;
        final double pixelSizeY;
        final String scanPattern;
        final double[] axis;
        final String title;
        final String isaInvestigationId;
        final double[] cube;
        // Bitset of seen pixels to guarantee width*height unique pixels
        // (catches duplicate (x,y) writes that would silently overwrite).
        final boolean[] seen;
        int seenCount;
        // Wire-mode marker forwarded from the IMAGE_HEADER so
        // appendPixel knows whether to parse a dense intensity vector
        // (continuous) or sparse channel/intensity pairs (processed).
        final boolean isContinuous;

        ImageBuilder(int width, int height, int spectralPoints,
                     double pixelSizeX, double pixelSizeY,
                     String scanPattern, double[] axis,
                     String title, String isaInvestigationId,
                     boolean isContinuous) {
            this.width = width;
            this.height = height;
            this.spectralPoints = spectralPoints;
            this.pixelSizeX = pixelSizeX;
            this.pixelSizeY = pixelSizeY;
            this.scanPattern = scanPattern;
            this.axis = axis != null ? axis : new double[0];
            this.title = title;
            this.isaInvestigationId = isaInvestigationId;
            this.cube = new double[(long) width * height
                                   * spectralPoints > Integer.MAX_VALUE
                ? Integer.MAX_VALUE
                : width * height * spectralPoints];
            this.seen = new boolean[width * height];
            this.isContinuous = isContinuous;
        }

        void setPixel(int x, int y, double[] intensities) {
            if (x < 0 || x >= width || y < 0 || y >= height) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL coordinates out of bounds: x=" + x
                    + ", y=" + y + " (width=" + width + ", height=" + height + ")");
            }
            if (intensities.length != spectralPoints) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL intensity count " + intensities.length
                    + " does not match IMAGE_HEADER.spectrum_bins="
                    + spectralPoints);
            }
            int pixelIdx = y * width + x;
            if (seen[pixelIdx]) {
                throw new IllegalStateException(
                    "duplicate IMAGE_PIXEL at (x=" + x + ", y=" + y + ")");
            }
            seen[pixelIdx] = true;
            seenCount++;
            int base = (y * width + x) * spectralPoints;
            System.arraycopy(intensities, 0, cube, base, spectralPoints);
        }

        long pixelsSeen() { return seenCount; }

        MSImage build() {
            return new MSImage(width, height, spectralPoints, 0,
                pixelSizeX, pixelSizeY, scanPattern,
                cube, axis,
                title, isaInvestigationId,
                List.of(), List.of(), List.of());
        }
    }

    /** v0.11 Task 5.3 (Deferral 1): Raman counterpart to
     *  {@link ImageBuilder}. The pixel-ingest logic is identical
     *  (continuous or processed mode, same axis dispatch); the
     *  only differences are the modality-specific extras
     *  ({@code excitationWavelengthNm}, {@code laserPowerMw}) and
     *  the {@link RamanImage} construction at build-time. */
    private static final class RamanImageBuilder {
        final int width;
        final int height;
        final int spectralPoints;
        final double pixelSizeX;
        final double pixelSizeY;
        final String scanPattern;
        final double[] axis;
        final String title;
        final String isaInvestigationId;
        final double excitationWavelengthNm;
        final double laserPowerMw;
        final double[] cube;
        final boolean[] seen;
        int seenCount;
        final boolean isContinuous;

        RamanImageBuilder(int width, int height, int spectralPoints,
                          double pixelSizeX, double pixelSizeY,
                          String scanPattern, double[] axis,
                          String title, String isaInvestigationId,
                          double excitationWavelengthNm, double laserPowerMw,
                          boolean isContinuous) {
            this.width = width;
            this.height = height;
            this.spectralPoints = spectralPoints;
            this.pixelSizeX = pixelSizeX;
            this.pixelSizeY = pixelSizeY;
            this.scanPattern = scanPattern;
            this.axis = axis != null ? axis : new double[0];
            this.title = title;
            this.isaInvestigationId = isaInvestigationId;
            this.excitationWavelengthNm = excitationWavelengthNm;
            this.laserPowerMw = laserPowerMw;
            this.cube = new double[width * height * spectralPoints];
            this.seen = new boolean[width * height];
            this.isContinuous = isContinuous;
        }

        void setPixel(int x, int y, double[] intensities) {
            if (x < 0 || x >= width || y < 0 || y >= height) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL (Raman) coordinates out of bounds: x="
                    + x + ", y=" + y);
            }
            if (intensities.length != spectralPoints) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL (Raman) intensity count "
                    + intensities.length + " does not match spectrum_bins="
                    + spectralPoints);
            }
            int pixelIdx = y * width + x;
            if (seen[pixelIdx]) {
                throw new IllegalStateException(
                    "duplicate IMAGE_PIXEL (Raman) at (x=" + x + ", y="
                    + y + ")");
            }
            seen[pixelIdx] = true;
            seenCount++;
            int base = (y * width + x) * spectralPoints;
            System.arraycopy(intensities, 0, cube, base, spectralPoints);
        }

        long pixelsSeen() { return seenCount; }

        RamanImage build() {
            return new RamanImage(width, height, spectralPoints, 0,
                pixelSizeX, pixelSizeY, scanPattern,
                excitationWavelengthNm, laserPowerMw,
                cube, axis,
                title, isaInvestigationId,
                List.of(), List.of(), List.of());
        }
    }

    /** v0.11 Task 5.3 (Deferral 1): IR counterpart to
     *  {@link ImageBuilder}. Modality-specific extras carried at
     *  build time are {@code mode} (TRANSMITTANCE / ABSORBANCE)
     *  and {@code resolutionCmInv}. */
    private static final class IRImageBuilder {
        final int width;
        final int height;
        final int spectralPoints;
        final double pixelSizeX;
        final double pixelSizeY;
        final String scanPattern;
        final double[] axis;
        final String title;
        final String isaInvestigationId;
        final IRMode mode;
        final double resolutionCmInv;
        final double[] cube;
        final boolean[] seen;
        int seenCount;
        final boolean isContinuous;

        IRImageBuilder(int width, int height, int spectralPoints,
                       double pixelSizeX, double pixelSizeY,
                       String scanPattern, double[] axis,
                       String title, String isaInvestigationId,
                       IRMode mode, double resolutionCmInv,
                       boolean isContinuous) {
            this.width = width;
            this.height = height;
            this.spectralPoints = spectralPoints;
            this.pixelSizeX = pixelSizeX;
            this.pixelSizeY = pixelSizeY;
            this.scanPattern = scanPattern;
            this.axis = axis != null ? axis : new double[0];
            this.title = title;
            this.isaInvestigationId = isaInvestigationId;
            this.mode = mode;
            this.resolutionCmInv = resolutionCmInv;
            this.cube = new double[width * height * spectralPoints];
            this.seen = new boolean[width * height];
            this.isContinuous = isContinuous;
        }

        void setPixel(int x, int y, double[] intensities) {
            if (x < 0 || x >= width || y < 0 || y >= height) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL (IR) coordinates out of bounds: x="
                    + x + ", y=" + y);
            }
            if (intensities.length != spectralPoints) {
                throw new IllegalStateException(
                    "IMAGE_PIXEL (IR) intensity count "
                    + intensities.length + " does not match spectrum_bins="
                    + spectralPoints);
            }
            int pixelIdx = y * width + x;
            if (seen[pixelIdx]) {
                throw new IllegalStateException(
                    "duplicate IMAGE_PIXEL (IR) at (x=" + x + ", y="
                    + y + ")");
            }
            seen[pixelIdx] = true;
            seenCount++;
            int base = (y * width + x) * spectralPoints;
            System.arraycopy(intensities, 0, cube, base, spectralPoints);
        }

        long pixelsSeen() { return seenCount; }

        IRImage build() {
            return new IRImage(width, height, spectralPoints, 0,
                pixelSizeX, pixelSizeY, scanPattern,
                mode, resolutionCmInv,
                cube, axis,
                title, isaInvestigationId,
                List.of(), List.of(), List.of());
        }
    }

    // ---------------------------------------------------------- helpers

    private static String readLEString(ByteBuffer buf, int widthBytes) {
        int len;
        if (widthBytes == 2) len = buf.getShort() & 0xFFFF;
        else                  len = buf.getInt();
        byte[] b = new byte[len];
        buf.get(b);
        return new String(b, StandardCharsets.UTF_8);
    }

    // ---------------------------------------------------------- accumulators

    /** Phase 2c-T: per-dataset verbatim v2 blob accumulator. */
    private static final class BulkV2BlobsBuilder {
        byte[] mateBlob;
        List<String> mateChromNames;
        byte[] nameTokBlob;
        byte[] refDiffBlob;
        String refDiffReferenceUri;

        boolean hasAny() {
            return mateBlob != null || nameTokBlob != null || refDiffBlob != null;
        }
    }

    private static final class DatasetMeta {
        final int datasetId;
        final String name;
        final int acquisitionMode;
        final String spectrumClass;
        final List<String> channelNames;
        final String instrumentJson;
        final long expectedAUCount;
        DatasetMeta(int id, String n, int mode, String cls, List<String> ch,
                    String instrumentJson, long exp) {
            datasetId = id; name = n; acquisitionMode = mode; spectrumClass = cls;
            channelNames = ch; this.instrumentJson = instrumentJson;
            expectedAUCount = exp;
        }
    }

    /** M89.2/M90.9: per-dataset accumulator for genomic AUs. Mirrors
     *  the Python {@code _new_genomic_accumulator} dict.
     *
     *  <p>M89.2: sequences and qualities ride as UINT8 channels; the
     *  suffix carries chromosome / position / mapq / flags.</p>
     *
     *  <p>M90.9: cigar / read_name / mate_chromosome ride as 3
     *  additional UINT8 string channels (per-AU). mate_position +
     *  template_length ride on the M90.9 mate extension at the end of
     *  the AU genomic suffix; the {@link AccessUnit#decode} path
     *  defaults them to -1 / 0 when absent (fixtures).</p>
     *
     *  <p>M90.10: dispatches on the wire {@code compression} byte to
     *  pick the M86 codec decoder (rANS / BASE_PACK) for the
     *  sequences + qualities channels. The 3 string channels are
     *  always uncompressed.</p> */
    private static final class GenomicAccumulator {
        final List<String> chromosomes = new ArrayList<>();
        final List<Long> positions = new ArrayList<>();
        final List<Integer> mappingQualities = new ArrayList<>();
        final List<Integer> flags = new ArrayList<>();
        final java.io.ByteArrayOutputStream sequences = new java.io.ByteArrayOutputStream();
        final java.io.ByteArrayOutputStream qualities = new java.io.ByteArrayOutputStream();
        final List<Long> offsets = new ArrayList<>();
        final List<Integer> lengths = new ArrayList<>();
        // M90.9 compound-field accumulators.
        final List<String> cigars = new ArrayList<>();
        final List<String> readNames = new ArrayList<>();
        final List<String> mateChroms = new ArrayList<>();
        final List<Long> matePositions = new ArrayList<>();
        final List<Integer> templateLengths = new ArrayList<>();
        // M101: per-read SAM tag text; anyTags marks a stream that
        // carried the tags channel at all.
        final List<String> tags = new ArrayList<>();
        boolean anyTags = false;
        long runningOffset = 0L;
        int acquisitionMode = 0;

        void ingest(AccessUnit au) {
            if (au.spectrumClass != 5) {
                throw new IllegalStateException(
                    "genomic accumulator received spectrum_class " + au.spectrumClass);
            }
            chromosomes.add(au.chromosome);
            positions.add(au.position);
            mappingQualities.add(au.mappingQuality);
            flags.add(au.flags);
            // mate extension fields ride on the AU genomic suffix.
            matePositions.add(au.matePosition);
            templateLengths.add(au.templateLength);
            int length = 0;
            // compound-string channels default to "" if absent
            // (an M89.2-era AU). Channel-name dispatch covers both
            // layouts.
            String cigarStr = "";
            String nameStr = "";
            String mateChrStr = "";
            String tagsStr = "";
            for (ChannelData ch : au.channels) {
                if (ch.precision != Enums.Precision.UINT8.ordinal()) {
                    throw new IllegalStateException(
                        "genomic channel precision " + ch.precision
                        + " not yet supported (UINT8 only)");
                }
                // dispatch on wire compression byte (NONE /
                // RANS_* / BASE_PACK). See decodeWireCodec.
                byte[] decoded = decodeWireCodec(ch.data, ch.compression);
                if ("sequences".equals(ch.name)) {
                    try { sequences.write(decoded); }
                    catch (java.io.IOException e) { throw new IllegalStateException(e); }
                    length = decoded.length;
                } else if ("qualities".equals(ch.name)) {
                    try { qualities.write(decoded); }
                    catch (java.io.IOException e) { throw new IllegalStateException(e); }
                    if (length == 0) length = decoded.length;
                } else if ("cigar".equals(ch.name)) {
                    cigarStr = new String(decoded, StandardCharsets.UTF_8);
                } else if ("read_name".equals(ch.name)) {
                    nameStr = new String(decoded, StandardCharsets.UTF_8);
                } else if ("mate_chromosome".equals(ch.name)) {
                    mateChrStr = new String(decoded, StandardCharsets.UTF_8);
                } else if ("tags".equals(ch.name)) {
                    tagsStr = new String(decoded, StandardCharsets.UTF_8);
                    anyTags = true;
                }
            }
            tags.add(tagsStr);
            cigars.add(cigarStr);
            readNames.add(nameStr);
            mateChroms.add(mateChrStr);
            offsets.add(runningOffset);
            lengths.add(length);
            runningOffset += length;
        }

        /** decode a wire payload encoded by
         *  {@code TransportWriter.applyWireCodec}. NONE → identity. */
        private static byte[] decodeWireCodec(byte[] payload, int codecId) {
            if (codecId == 0) return payload;  // NONE
            if (codecId == Enums.Compression.RANS_ORDER0.ordinal()
                    || codecId == Enums.Compression.RANS_ORDER1.ordinal()) {
                return Rans.decode(payload);
            }
            if (codecId == Enums.Compression.BASE_PACK.ordinal()) {
                return BasePack.decode(payload);
            }
            throw new UnsupportedOperationException(
                "decodeWireCodec: codec id " + codecId
                + " not supported for genomic UINT8");
        }

        WrittenGenomicRun toWrittenGenomicRun(DatasetMeta meta) {
            int n = chromosomes.size();
            long[] offsetsArr = new long[n];
            int[] lengthsArr = new int[n];
            long[] positionsArr = new long[n];
            byte[] mqArr = new byte[n];
            int[] flagsArr = new int[n];
            for (int i = 0; i < n; i++) {
                offsetsArr[i] = offsets.get(i);
                lengthsArr[i] = lengths.get(i);
                positionsArr[i] = positions.get(i);
                mqArr[i] = (byte) (mappingQualities.get(i) & 0xFF);
                flagsArr[i] = flags.get(i);
            }
            // compound fields now round-trip on the wire. When
            // the source is an M89.2-era stream the per-AU decoders
            // default the missing strings to "" and the mate scalars
            // to -1 / 0 (preserved by AccessUnit.decode + the
            // accumulator defaults).
            long[] mateP = new long[n];
            int[] tlens = new int[n];
            for (int i = 0; i < n; i++) {
                mateP[i] = matePositions.get(i);
                tlens[i] = templateLengths.get(i);
            }
            List<String> cigarsOut = new ArrayList<>(cigars);
            List<String> readNamesOut = new ArrayList<>(readNames);
            List<String> mateChromsOut = new ArrayList<>(mateChroms);

            // Decode instrument_json metadata.
            String referenceUri = "", platform = "", sampleName = "", modality = "";
            String readRole = "";
            try {
                Object parsed = MiniJson.parse(meta.instrumentJson);
                if (parsed instanceof Map<?, ?> mraw) {
                    @SuppressWarnings("unchecked")
                    Map<String, Object> m = (Map<String, Object>) mraw;
                    Object ru = m.get("reference_uri"); if (ru != null) referenceUri = ru.toString();
                    Object pl = m.get("platform"); if (pl != null) platform = pl.toString();
                    Object sn = m.get("sample_name"); if (sn != null) sampleName = sn.toString();
                    Object md = m.get("modality"); if (md != null) modality = md.toString();
                    // M97: "" on the wire means no role.
                    Object rr = m.get("read_role"); if (rr != null) readRole = rr.toString();
                }
            } catch (Exception ignore) {
                // Unparseable instrument_json - leave defaults.
            }
            // modality not currently surfaced on WrittenGenomicRun's
            // constructor; the GenomicRun reader pulls it from the
            // file-level modality attribute (defaulted at write time).
            Enums.AcquisitionMode acqMode;
            try {
                acqMode = Enums.AcquisitionMode.values()[
                    Math.min(meta.acquisitionMode,
                             Enums.AcquisitionMode.values().length - 1)];
            } catch (Exception e) {
                acqMode = Enums.AcquisitionMode.GENOMIC_WGS;
            }
            WrittenGenomicRun decoded = new WrittenGenomicRun(
                acqMode, referenceUri, platform, sampleName,
                positionsArr, mqArr, flagsArr,
                sequences.toByteArray(), qualities.toByteArray(),
                offsetsArr, lengthsArr,
                cigarsOut, readNamesOut, mateChromsOut, mateP, tlens,
                new ArrayList<>(chromosomes),
                Enums.Compression.ZLIB);
            if (anyTags) decoded = decoded.withTags(new ArrayList<>(tags));
            return readRole.isEmpty() ? decoded
                                      : decoded.withReadRole(readRole);
        }
    }

    private static final class RunAccumulator {
        final List<String> channelNames;
        final Map<String, List<double[]>> perSpectrumChannels = new LinkedHashMap<>();
        long runningOffset = 0;
        final List<Long> offsets = new ArrayList<>();
        final List<Integer> lengths = new ArrayList<>();
        final List<Double> retentionTimes = new ArrayList<>();
        final List<Integer> msLevels = new ArrayList<>();
        final List<Integer> polarities = new ArrayList<>();
        final List<Double> precursorMzs = new ArrayList<>();
        final List<Integer> precursorCharges = new ArrayList<>();
        final List<Double> basePeakIntensities = new ArrayList<>();

        RunAccumulator(List<String> channelNames) {
            this.channelNames = channelNames;
            for (String c : channelNames) perSpectrumChannels.put(c, new ArrayList<>());
        }

        void ingest(AccessUnit au) {
            int length = 0;
            Map<String, double[]> perAu = new LinkedHashMap<>();
            for (ChannelData ch : au.channels) {
                if (ch.precision != Enums.Precision.FLOAT64.ordinal()) {
                    throw new IllegalStateException(
                            "reader supports FLOAT64 precision only");
                }
                double[] arr;
                if (ch.compression == Enums.Compression.FLOAT_DELTA_ZSTD.ordinal()) {
                    arr = global.thalion.ttio.codecs.FloatDeltaZstd.decode(ch.data);
                    if (arr.length != ch.nElements) {
                        throw new IllegalStateException(
                                "FLOAT_DELTA_ZSTD channel decoded " + arr.length
                                + " values but the channel header declares n_elements="
                                + ch.nElements);
                    }
                } else {
                    byte[] raw;
                    if (ch.compression == Enums.Compression.NONE.ordinal()) {
                        raw = ch.data;
                    } else if (ch.compression == Enums.Compression.ZLIB.ordinal()) {
                        raw = zlibInflate(ch.data);
                    } else if (ch.compression == Enums.Compression.ZSTD.ordinal()) {
                        // The channel header carries the exact plaintext
                        // element count, so decompress into an exact buffer.
                        raw = zstdDecompress(ch.data, ch.nElements * 8);
                    } else {
                        throw new IllegalStateException(
                                "reader supports NONE/ZLIB/ZSTD/FLOAT_DELTA_ZSTD "
                                + "compression only, got " + ch.compression);
                    }
                    int n = raw.length / 8;
                    arr = new double[n];
                    ByteBuffer buf = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
                    for (int k = 0; k < n; k++) arr[k] = buf.getDouble();
                }
                int n = arr.length;
                perAu.put(ch.name, arr);
                if (length != 0 && length != n) {
                    throw new IllegalStateException(
                            "channels in one AU have mismatched lengths");
                }
                length = n;
            }
            offsets.add(runningOffset);
            lengths.add(length);
            runningOffset += length;
            for (String c : channelNames) {
                double[] arr = perAu.get(c);
                if (arr == null) arr = new double[length];
                perSpectrumChannels.get(c).add(arr);
            }
            retentionTimes.add(au.retentionTime);
            msLevels.add(au.msLevel);
            polarities.add(wireToPolarityInt(au.polarity));
            precursorMzs.add(au.precursorMz);
            precursorCharges.add(au.precursorCharge);
            basePeakIntensities.add(au.basePeakIntensity);
        }

        SpectrumIndex toSpectrumIndex() {
            int n = offsets.size();
            long[] off = new long[n];
            int[] len = new int[n];
            double[] rt = new double[n];
            int[] ms = new int[n];
            int[] pol = new int[n];
            double[] pmz = new double[n];
            int[] pc = new int[n];
            double[] bpi = new double[n];
            for (int i = 0; i < n; i++) {
                off[i] = offsets.get(i);
                len[i] = lengths.get(i);
                rt[i] = retentionTimes.get(i);
                ms[i] = msLevels.get(i);
                pol[i] = polarities.get(i);
                pmz[i] = precursorMzs.get(i);
                pc[i] = precursorCharges.get(i);
                bpi[i] = basePeakIntensities.get(i);
            }
            return new SpectrumIndex(n, off, len, rt, ms, pol, pmz, pc, bpi);
        }

        Map<String, double[]> toChannelMap() {
            Map<String, double[]> out = new LinkedHashMap<>();
            for (String c : channelNames) {
                List<double[]> chunks = perSpectrumChannels.get(c);
                int total = 0;
                for (double[] d : chunks) total += d.length;
                double[] flat = new double[total];
                int off = 0;
                for (double[] d : chunks) {
                    System.arraycopy(d, 0, flat, off, d.length);
                    off += d.length;
                }
                out.put(c, flat);
            }
            return out;
        }
    }

    private static byte[] zstdDecompress(byte[] input, int plainLength) {
        io.airlift.compress.zstd.ZstdDecompressor d =
                new io.airlift.compress.zstd.ZstdDecompressor();
        byte[] out = new byte[plainLength];
        int n = d.decompress(input, 0, input.length, out, 0, out.length);
        if (n != plainLength) {
            throw new IllegalStateException(
                    "zstd payload inflated to " + n + " bytes; expected "
                    + plainLength);
        }
        return out;
    }

    private static byte[] zlibInflate(byte[] input) {
        java.util.zip.Inflater inf = new java.util.zip.Inflater();
        inf.setInput(input);
        byte[] buf = new byte[Math.max(64, input.length * 4)];
        java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
        try {
            while (!inf.finished()) {
                int n = inf.inflate(buf);
                if (n == 0) {
                    if (inf.needsInput() || inf.needsDictionary()) {
                        throw new IllegalStateException("zlib underflow");
                    }
                }
                out.write(buf, 0, n);
            }
        } catch (java.util.zip.DataFormatException e) {
            throw new IllegalStateException("zlib inflate failed", e);
        } finally {
            inf.end();
        }
        return out.toByteArray();
    }

    private static int wireToPolarityInt(int wire) {
        return switch (wire) {
            case 0 -> Enums.Polarity.POSITIVE.intValue();
            case 1 -> Enums.Polarity.NEGATIVE.intValue();
            default -> Enums.Polarity.UNKNOWN.intValue();
        };
    }
}
