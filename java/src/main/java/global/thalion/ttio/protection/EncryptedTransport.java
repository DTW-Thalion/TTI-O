/* TTI-O Java Implementation / Copyright (c) 2026 The Thalion Initiative / SPDX-License-Identifier: LGPL-3.0-or-later */
package global.thalion.ttio.protection;

import global.thalion.ttio.Enums;
import global.thalion.ttio.FeatureFlags;
import global.thalion.ttio.protection.PerAUEncryption.ChannelSegment;
import global.thalion.ttio.protection.PerAUEncryption.HeaderSegment;
import global.thalion.ttio.providers.CompoundField;
import global.thalion.ttio.providers.ProviderRegistry;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;
import global.thalion.ttio.transport.PacketHeader;
import global.thalion.ttio.transport.PacketType;
import global.thalion.ttio.transport.TransportReader;
import global.thalion.ttio.transport.TransportWriter;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * v1.0 transport layer for per-AU encrypted files. Writer pushes
 * ciphertext bytes from the {@code <channel>_segments} compound onto
 * the wire unmodified (the server never decrypts in transit, per
 * {@code docs/transport-spec.md} §6.2). Reader materialises an
 * encrypted {@code .tio} on the receiver side preserving ciphertext
 * bytes verbatim.
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.transport.encrypted}, Objective-C
 * {@code TTIOEncryptedTransport}.
 *
 *
 */
public final class EncryptedTransport {

    private EncryptedTransport() {}

    // ────────────────────────────────────────────────────────── probe

    /** {@code true} iff the file at {@code path} carries
     *  {@code opt_per_au_encryption}. */
    public static boolean isPerAUEncrypted(String path, String providerName) {
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);
            return flags.has(FeatureFlags.OPT_PER_AU_ENCRYPTION);
        } finally {
            sp.close();
        }
    }

    // ────────────────────────────────────────────────────────── writer

    /** Emit the full transport stream from a per-AU-encrypted file.
     *
     *  <p>M90.8: also walks {@code /study/genomic_runs/} after MS runs.
     *  Genomic dataset_id continues from MS (1..N MS, N+1..N+M genomic)
     *  so AAD reconstruction matches the per-AU encrypt path
     *  (M90.1).</p> */
    public static void writeEncryptedDataset(String ttioPath,
                                               TransportWriter writer,
                                               String providerName)
            throws IOException {
        StorageProvider sp = ProviderRegistry.open(ttioPath,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            FeatureFlags flags = FeatureFlags.readFrom(root);
            if (!flags.has(FeatureFlags.OPT_PER_AU_ENCRYPTION)) {
                throw new IllegalStateException(
                    ttioPath + " does not carry opt_per_au_encryption");
            }
            boolean headersEncrypted = flags.has(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);

            try (StorageGroup study = root.openGroup("study")) {
                String title = attrStr(study, "title", "");
                String isa = attrStr(study, "isa_investigation_id", "");

                List<String> msRunNames = new ArrayList<>();
                if (study.hasChild("ms_runs")) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        for (String n : msRuns.childNames()) {
                            if (!n.startsWith("_") && msRuns.hasChild(n)) {
                                msRunNames.add(n);
                            }
                        }
                    }
                }
                List<String> genomicRunNames = new ArrayList<>();
                Map<String, String> genomicLayouts = new LinkedHashMap<>();
                if (study.hasChild("genomic_runs")) {
                    try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                        for (String n : gRuns.childNames()) {
                            if (!n.startsWith("_") && gRuns.hasChild(n)) {
                                genomicRunNames.add(n);
                            }
                        }
                        for (String n : genomicRunNames) {
                            try (StorageGroup g = gRuns.openGroup(n)) {
                                genomicLayouts.put(n,
                                    g.hasAttribute("layout")
                                        ? g.getAttribute("layout").toString()
                                        : "");
                            }
                        }
                    }
                }

                // blocks_v1 runs ride the same per-read AU stream plus
                // the M99.1 sidecar packets (GenomicRunSidecar + one
                // BlockSidecar per block); the StreamHeader announces
                // them with a wire feature token so receivers know the
                // stream needs them.
                List<String> features = new ArrayList<>(flags.features());
                if (genomicLayouts.containsValue("blocks_v1")
                        || genomicLayouts.containsValue("blocks_v1_grouped")) {
                    features.add(PacketType.TRANSPORT_BLOCKS_V1_FEATURE);
                }
                // v0.13 (M103): a grouped run's BlockSidecars end with the
                // block's slice of genomic_index/input_index.
                if (genomicLayouts.containsValue("blocks_v1_grouped")) {
                    features.add(PacketType.TRANSPORT_BLOCKS_V1_GROUPED_FEATURE);
                }
                writer.writeStreamHeader("1.2", title, isa, features,
                                          msRunNames.size() + genomicRunNames.size());
                int did = 1;
                if (!msRunNames.isEmpty()) {
                    try (StorageGroup msRuns = study.openGroup("ms_runs")) {
                        for (String runName : msRunNames) {
                            try (StorageGroup run = msRuns.openGroup(runName);
                                 StorageGroup sig = run.openGroup("signal_channels")) {
                                List<String> channelNames = splitNames(
                                    attrStr(sig, "channel_names", ""));
                                String firstCh = channelNames.isEmpty()
                                    ? "intensity" : channelNames.get(0);
                                String cipherSuite = attrStr(sig,
                                    firstCh + "_algorithm", "aes-256-gcm");
                                String kek = attrStr(sig,
                                    firstCh + "_kek_algorithm", "");
                                byte[] wrapped = readWrappedDekAttr(sig,
                                    firstCh + "_wrapped_dek");
                                List<Recipient> extra = decodeRecipientBlockBytes(
                                    readWrappedDekAttr(sig,
                                        firstCh + "_wrapped_dek_recipients"));

                                writer.emitRawPacket(PacketType.PROTECTION_METADATA, 0,
                                    did, 0, encodeProtection(cipherSuite, kek, wrapped, extra,
                                        attrStr(sig, firstCh + "_server_kek_id", "")));

                                long expectedAUs = firstChannelSegmentCount(sig, firstCh);
                                int acqMode = intAttr(run, "acquisition_mode", 0);
                                String spectrumClass = attrStr(run, "spectrum_class",
                                                                 "TTIOMassSpectrum");
                                writer.writeDatasetHeader(did, runName, acqMode,
                                    spectrumClass, channelNames, "{}", expectedAUs);
                            }
                            did++;
                        }

                        did = 1;
                        for (String runName : msRunNames) {
                            long n = emitRunAUs(writer, msRuns, runName, did,
                                                  headersEncrypted);
                            writer.writeEndOfDataset(did, n);
                            did++;
                        }
                    }
                }

                if (!genomicRunNames.isEmpty()) {
                    try (StorageGroup gRuns = study.openGroup("genomic_runs")) {
                        int firstGenomicDid = did;
                        for (String runName : genomicRunNames) {
                            emitGenomicDatasetHeader(writer, gRuns, runName, did);
                            if (global.thalion.ttio.genomics.GenomicStreamWriter
                                    .isBlocksLayout(genomicLayouts.get(runName))) {
                                emitBlocksV1Sidecars(writer, gRuns,
                                                     runName, did);
                            }
                            did++;
                        }
                        did = firstGenomicDid;
                        for (String runName : genomicRunNames) {
                            long n = emitGenomicRunAUs(writer, gRuns, runName, did);
                            writer.writeEndOfDataset(did, n);
                            did++;
                        }
                    }
                }
            }
            writer.writeEndOfStream();
        } finally {
            sp.close();
        }
    }

    /** emit ProtectionMetadata + DatasetHeader for one genomic
     *  run. Genomic only encrypts {@code sequences} + {@code qualities}
     *  (per M90.1); other channels (cigars, read_names, mate_info)
     *  stay plaintext on the source file and are not part of the
     *  encrypted-transport contract. */
    private static void emitGenomicDatasetHeader(TransportWriter writer,
                                                   StorageGroup gRuns,
                                                   String runName,
                                                   int datasetId)
            throws IOException {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels")) {
            // Genomic encrypts sequences + qualities (M90.1) and, when
            // the run carries them, the SAM tags (M101).
            List<String> channelNames = new ArrayList<>(3);
            for (String c : new String[]{"sequences", "qualities", "tags"}) {
                if (sig.hasChild(c + "_segments")) channelNames.add(c);
            }
            String firstCh = channelNames.isEmpty()
                ? "sequences" : channelNames.get(0);
            String cipherSuite = attrStr(sig,
                firstCh + "_algorithm", "aes-256-gcm");
            String kek = attrStr(sig, firstCh + "_kek_algorithm", "");
            byte[] wrapped = readWrappedDekAttr(sig, firstCh + "_wrapped_dek");
            List<Recipient> extra = decodeRecipientBlockBytes(
                readWrappedDekAttr(sig, firstCh + "_wrapped_dek_recipients"));
            writer.emitRawPacket(PacketType.PROTECTION_METADATA, 0,
                datasetId, 0, encodeProtection(cipherSuite, kek, wrapped, extra,
                    attrStr(sig, firstCh + "_server_kek_id", "")));

            int acqMode = intAttr(run, "acquisition_mode", 0);
            String metadataJson = genomicRunMetadataJson(
                attrStr(run, "modality", ""),
                attrStr(run, "platform", ""),
                attrStr(run, "read_role", ""),
                attrStr(run, "reference_uri", ""),
                attrStr(run, "sample_name", ""));
            long expectedAUs = channelNames.isEmpty()
                ? 0L
                : PerAUFile.readChannelSegments(sig,
                      channelNames.get(0) + "_segments").size();
            writer.writeDatasetHeader(datasetId, runName, acqMode,
                "TTIOGenomicRead", channelNames, metadataJson, expectedAUs);
        }
    }

    /** blocks_v1 channels whose coded blobs stay plaintext under
     *  per-AU protection and ride in BlockSidecar packets. */
    private static final List<String> SIDECAR_BLOB_CHANNELS =
        List.of("read_names", "cigars", "mate_info");

    private static void putLEString(java.io.ByteArrayOutputStream out,
                                    String s) {
        byte[] d = (s == null ? "" : s).getBytes(StandardCharsets.UTF_8);
        out.write(d.length & 0xFF);
        out.write((d.length >> 8) & 0xFF);
        out.writeBytes(d);
    }

    private static void putLE(java.io.ByteArrayOutputStream out,
                              long v, int width) {
        for (int i = 0; i < width; i++) {
            out.write((int) ((v >>> (8 * i)) & 0xFF));
        }
    }

    /** One GenomicRunSidecar, then one BlockSidecar per block, for a
     *  blocks_v1 per-AU-encrypted genomic run (transport-spec
     *  §4.24). */
    private static void emitBlocksV1Sidecars(TransportWriter writer,
                                             StorageGroup gRuns,
                                             String runName,
                                             int datasetId)
            throws IOException {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels")) {
            global.thalion.ttio.genomics.BlockTable t =
                global.thalion.ttio.genomics.BlockTable.read(run);

            Map<String, Object> attrs = new TreeMap<>();
            long slice = longAttr(run, "ref_diff_slice_bytes");
            if (slice != 0) attrs.put("ref_diff_slice_bytes", slice);
            if (longAttr(run, "opt_disable_qualities_v5") != 0) {
                attrs.put("opt_disable_qualities_v5", 1);
            }
            String md5s = attrStr(run, "reference_md5s", "");
            if (!md5s.isEmpty()) attrs.put("reference_md5s", md5s);

            List<Map<String, Object>> channels = new ArrayList<>();
            Map<String, StorageDataset> blobDs = new LinkedHashMap<>();
            try {
                for (String ch : SIDECAR_BLOB_CHANNELS) {
                    StorageDataset ds = null;
                    if (ch.equals("mate_info")) {
                        if (sig.hasChild("mate_info")) {
                            try (StorageGroup mg =
                                    sig.openGroup("mate_info")) {
                                if (mg.hasChild("inline_v2")) {
                                    ds = mg.openDataset("inline_v2");
                                }
                            }
                        }
                    } else if (sig.hasChild(ch)) {
                        ds = sig.openDataset(ch);
                    }
                    if (ds == null) continue;
                    blobDs.put(ch, ds);
                    Map<String, Object> centry = new TreeMap<>();
                    centry.put("name", ch);
                    centry.put("compression", ds.hasAttribute("compression")
                        ? ((Number) ds.getAttribute("compression"))
                              .longValue()
                        : 0L);
                    Map<String, Object> extra = new TreeMap<>();
                    for (String k : ds.attributeNames()) {
                        if (k.equals("compression")) continue;
                        Object v = ds.getAttribute(k);
                        if (v instanceof byte[] bb) {
                            v = new String(bb, StandardCharsets.UTF_8);
                        }
                        if (v != null) extra.put(k, v);
                    }
                    if (!extra.isEmpty()) centry.put("extra_attrs", extra);
                    channels.add(centry);
                }

                List<String> chromNames;
                try (StorageGroup idx = run.openGroup("genomic_index")) {
                    chromNames = global.thalion.ttio.genomics.BlockView
                        .readNames(idx, "chromosome_names");
                }
                List<String> mateNames = List.of();
                if (sig.hasChild("mate_info")) {
                    try (StorageGroup mg = sig.openGroup("mate_info")) {
                        if (mg.hasChild("chrom_names")) {
                            mateNames = global.thalion.ttio.genomics
                                .BlockView.readNames(mg, "chrom_names");
                        }
                    }
                }
                // M103: a grouped run carries each block's input_index
                // rows at the end of its BlockSidecar (v0.13).
                long[] inputIndex = null;
                if ("blocks_v1_grouped".equals(attrStr(run, "layout", ""))) {
                    try (StorageGroup idx = run.openGroup("genomic_index");
                         StorageDataset ds = idx.openDataset("input_index")) {
                        inputIndex = toUnsignedLongs(ds.readAll());
                    }
                }

                com.fasterxml.jackson.databind.ObjectMapper om =
                    new com.fasterxml.jackson.databind.ObjectMapper();
                byte[] attrsJson = om.writeValueAsBytes(attrs);
                byte[] channelsJson = om.writeValueAsBytes(channels);

                java.io.ByteArrayOutputStream p =
                    new java.io.ByteArrayOutputStream();
                putLEString(p, attrStr(run, "layout", ""));
                putLEString(p, attrStr(run, "block_policy", ""));
                putLE(p, longAttr(run, "read_count"), 8);
                putLE(p, longAttr(run, "base_count"), 8);
                putLE(p, attrsJson.length, 4);
                p.writeBytes(attrsJson);
                putLE(p, channelsJson.length, 4);
                p.writeBytes(channelsJson);
                putLE(p, chromNames.size(), 4);
                for (String n : chromNames) putLEString(p, n);
                putLE(p, mateNames.size(), 4);
                for (String n : mateNames) putLEString(p, n);
                writer.emitRawPacket(PacketType.GENOMIC_RUN_SIDECAR, 0,
                                     datasetId, 0, p.toByteArray());

                // The index's own channel set: a run without tags sends
                // the five required entries, exactly as before M101.
                List<String> blockChannels = t.channels();
                for (int b = 0; b < t.count(); b++) {
                    java.io.ByteArrayOutputStream bp =
                        new java.io.ByteArrayOutputStream();
                    putLE(bp, b, 4);
                    putLE(bp, t.readStartAt(b), 8);
                    putLE(bp, t.nReadsAt(b), 4);
                    putLE(bp, t.baseStartAt(b), 8);
                    putLE(bp, t.nBasesAt(b), 8);
                    bp.write(blockChannels.size() & 0xFF);
                    for (String ch : blockChannels) {
                        putLEString(bp, ch);
                        putLE(bp, t.offsetOf(ch, b), 8);
                        putLE(bp, t.lengthOf(ch, b), 8);
                        putLE(bp, t.hasCodecs() ? t.codecOf(ch, b) : 0, 4);
                    }
                    bp.write(blobDs.size() & 0xFF);
                    for (Map.Entry<String, StorageDataset> e
                            : blobDs.entrySet()) {
                        String ch = e.getKey();
                        long off = t.offsetOf(ch, b);
                        long ln = t.lengthOf(ch, b);
                        byte[] data = ln > 0
                            ? (byte[]) e.getValue().readSlice(off, ln)
                            : new byte[0];
                        putLEString(bp, ch);
                        putLE(bp, data.length, 4);
                        bp.writeBytes(data);
                    }
                    if (inputIndex != null) {
                        // v0.13: the block's rows of genomic_index/input_index.
                        int r0 = (int) t.readStartAt(b);
                        for (int k = 0; k < t.nReadsAt(b); k++) {
                            putLE(bp, inputIndex[r0 + k], 4);
                        }
                    }
                    writer.emitRawPacket(PacketType.BLOCK_SIDECAR, 0,
                                         datasetId, b, bp.toByteArray());
                }
            } finally {
                for (StorageDataset ds : blobDs.values()) ds.close();
            }
        }
    }

    /** A uint32 dataset's values, whatever array type the provider
     *  hands back. */
    private static long[] toUnsignedLongs(Object raw) {
        if (raw instanceof int[] a) {
            long[] out = new long[a.length];
            for (int i = 0; i < a.length; i++) out[i] = Integer.toUnsignedLong(a[i]);
            return out;
        }
        if (raw instanceof long[] l) return l;
        throw new IllegalStateException("input_index has an unexpected type");
    }

    private static long longAttr(StorageGroup g, String name) {
        if (!g.hasAttribute(name)) return 0L;
        Object v = g.getAttribute(name);
        return v instanceof Number n ? n.longValue() : 0L;
    }

    /** Decoded GenomicRunSidecar packet (transport-spec §4.24). */
    static final class RunSidecar {
        String layout = "";
        String blockPolicy = "";
        long readCount;
        long baseCount;
        Map<String, Object> attrs = Map.of();
        List<Map<String, Object>> channels = List.of();
        List<String> chromNames = List.of();
        List<String> mateChromNames = List.of();
    }

    /** Decoded BlockSidecar packet (transport-spec §4.24). */
    static final class BlockSidecar {
        int blockIndex;
        long readStart;
        int nReads;
        long baseStart;
        long nBases;
        Map<String, long[]> channels = new LinkedHashMap<>();
        Map<String, byte[]> blobs = new LinkedHashMap<>();
        /** v0.13 (M103): the block's input_index rows of a
         *  blocks_v1_grouped run; null when the sidecar carries none. */
        int[] inputIndex;
    }

    private static String getLEString(ByteBuffer bb) {
        int n = bb.getShort() & 0xFFFF;
        byte[] d = new byte[n];
        bb.get(d);
        return new String(d, StandardCharsets.UTF_8);
    }

    @SuppressWarnings("unchecked")
    private static RunSidecar parseRunSidecar(byte[] payload)
            throws IOException {
        ByteBuffer bb = ByteBuffer.wrap(payload)
            .order(ByteOrder.LITTLE_ENDIAN);
        RunSidecar sc = new RunSidecar();
        sc.layout = getLEString(bb);
        sc.blockPolicy = getLEString(bb);
        sc.readCount = bb.getLong();
        sc.baseCount = bb.getLong();
        com.fasterxml.jackson.databind.ObjectMapper om =
            new com.fasterxml.jackson.databind.ObjectMapper();
        int jl = bb.getInt();
        byte[] jd = new byte[jl];
        bb.get(jd);
        sc.attrs = om.readValue(jd, Map.class);
        jl = bb.getInt();
        jd = new byte[jl];
        bb.get(jd);
        sc.channels = om.readValue(jd, List.class);
        int n = bb.getInt();
        List<String> chroms = new ArrayList<>(n);
        for (int i = 0; i < n; i++) chroms.add(getLEString(bb));
        sc.chromNames = chroms;
        n = bb.getInt();
        List<String> mates = new ArrayList<>(n);
        for (int i = 0; i < n; i++) mates.add(getLEString(bb));
        sc.mateChromNames = mates;
        return sc;
    }

    private static BlockSidecar parseBlockSidecar(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload)
            .order(ByteOrder.LITTLE_ENDIAN);
        BlockSidecar bs = new BlockSidecar();
        bs.blockIndex = bb.getInt();
        bs.readStart = bb.getLong();
        bs.nReads = bb.getInt();
        bs.baseStart = bb.getLong();
        bs.nBases = bb.getLong();
        int nch = bb.get() & 0xFF;
        for (int i = 0; i < nch; i++) {
            String name = getLEString(bb);
            long off = bb.getLong();
            long len = bb.getLong();
            int codec = bb.getInt();
            bs.channels.put(name, new long[]{ off, len, codec });
        }
        int nb = bb.get() & 0xFF;
        for (int i = 0; i < nb; i++) {
            String name = getLEString(bb);
            int bl = bb.getInt();
            byte[] data = new byte[bl];
            bb.get(data);
            bs.blobs.put(name, data);
        }
        // v0.13 (M103): a blocks_v1_grouped run's sidecar ends with the
        // block's input_index rows, n_reads x uint32.
        if (bb.hasRemaining()) {
            long want = 4L * Integer.toUnsignedLong(bs.nReads);
            if (bb.remaining() != want) {
                throw new IllegalStateException("BlockSidecar " + bs.blockIndex + ": "
                    + bb.remaining() + " trailing bytes, expected " + want
                    + " of input_index");
            }
            bs.inputIndex = new int[bs.nReads];
            for (int k = 0; k < bs.nReads; k++) bs.inputIndex[k] = bb.getInt();
        }
        return bs;
    }

    /** emit one ENCRYPTED ACCESS_UNIT packet per read.
     *  Each AU carries {@code spectrum_class=5}, the M89.1 chromosome
     *  + position + mapq + flags suffix (sourced from the plaintext
     *  {@code genomic_index/}), and UINT8 ChannelData with
     *  {@code IV || TAG || ciphertext} for each encrypted channel. */
    private static long emitGenomicRunAUs(TransportWriter writer,
                                            StorageGroup gRuns,
                                            String runName,
                                            int datasetId)
            throws IOException {
        try (StorageGroup run = gRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("genomic_index")) {
            // Genomic encrypts sequences + qualities (M90.1) and, when
            // the run carries them, the SAM tags (M101).
            List<String> channelNames = new ArrayList<>(3);
            for (String c : new String[]{"sequences", "qualities", "tags"}) {
                if (sig.hasChild(c + "_segments")) channelNames.add(c);
            }
            int acqMode = intAttr(run, "acquisition_mode", 0);
            Map<String, List<ChannelSegment>> segsByCh = new LinkedHashMap<>();
            for (String c : channelNames) {
                segsByCh.put(c, PerAUFile.readChannelSegments(sig, c + "_segments"));
            }
            long[] positions;
            byte[] mapqs;
            int[] flagsArr;
            try (StorageDataset d = idx.openDataset("positions")) {
                positions = (long[]) d.readAll();
            }
            try (StorageDataset d = idx.openDataset("mapping_qualities")) {
                mapqs = (byte[]) d.readAll();
            }
            try (StorageDataset d = idx.openDataset("flags")) {
                flagsArr = (int[]) d.readAll();
            }
            List<String> chromosomes = readGenomicChromosomes(idx);
            int n = channelNames.isEmpty() ? 0 : segsByCh.get(channelNames.get(0)).size();
            int uint8Precision = global.thalion.ttio.Enums.Precision.UINT8.ordinal();

            for (int i = 0; i < n; i++) {
                List<global.thalion.ttio.transport.ChannelData> channels =
                    new ArrayList<>(channelNames.size());
                for (String c : channelNames) {
                    ChannelSegment seg = segsByCh.get(c).get(i);
                    byte[] data = concat(seg.iv(), seg.tag(), seg.ciphertext());
                    channels.add(new global.thalion.ttio.transport.ChannelData(
                        c, uint8Precision, 0, seg.length(), data));
                }
                global.thalion.ttio.transport.AccessUnit au =
                    new global.thalion.ttio.transport.AccessUnit(
                        5, acqMode, 0, 2,
                        0.0, 0.0, 0,
                        0.0, 0.0,
                        channels,
                        0L, 0L, 0L,
                        chromosomes.get(i),
                        positions[i],
                        mapqs[i] & 0xFF,
                        flagsArr[i] & 0xFFFF);
                writer.emitRawPacket(PacketType.ACCESS_UNIT,
                    PacketHeader.FLAG_ENCRYPTED, datasetId, i, au.encode());
            }
            return n;
        }
    }

    @SuppressWarnings("unchecked")
    private static List<String> readGenomicChromosomes(StorageGroup idx) {
        // L1 (Task #82 Phase B.1): chromosomes are stored as
        // chromosome_ids (uint16) + chromosome_names (compound).
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
            if (v == null) nameTable.add("");
            else if (v instanceof byte[] b) nameTable.add(new String(b, StandardCharsets.UTF_8));
            else nameTable.add(v.toString());
        }
        {
            List<String> out = new ArrayList<>(ids.length);
            for (short id : ids) {
                int idx2 = Short.toUnsignedInt(id);
                out.add(idx2 < nameTable.size() ? nameTable.get(idx2) : "");
            }
            return out;
        }
    }

    private static String genomicRunMetadataJson(String modality, String platform,
                                                   String readRole,
                                                   String referenceUri, String sampleName) {
        char QT = (char) 34;
        char OB = (char) 123;
        char CB = (char) 125;
        StringBuilder sb = new StringBuilder(96);
        sb.append(OB);
        sb.append(QT).append("modality").append(QT).append(": ").append(QT).append(esc(modality)).append(QT);
        sb.append(", ");
        sb.append(QT).append("platform").append(QT).append(": ").append(QT).append(esc(platform)).append(QT);
        sb.append(", ");
        // M97: always present, "" when the run has no @read_role.
        sb.append(QT).append("read_role").append(QT).append(": ").append(QT).append(esc(readRole)).append(QT);
        sb.append(", ");
        sb.append(QT).append("reference_uri").append(QT).append(": ").append(QT).append(esc(referenceUri)).append(QT);
        sb.append(", ");
        sb.append(QT).append("sample_name").append(QT).append(": ").append(QT).append(esc(sampleName)).append(QT);
        sb.append(CB);
        return sb.toString();
    }

    private static String esc(String value) {
        if (value == null) return "";
        char QT = (char) 34;
        char BS = (char) 92;
        StringBuilder out = new StringBuilder(value.length());
        for (int i = 0; i < value.length(); i++) {
            char c = value.charAt(i);
            if (c == QT || c == BS) out.append(BS);
            out.append(c);
        }
        return out.toString();
    }

    /** @return number of AU packets emitted for this run. */
    private static long emitRunAUs(TransportWriter writer, StorageGroup msRuns,
                                     String runName, int datasetId,
                                     boolean headersEncrypted) throws IOException {
        try (StorageGroup run = msRuns.openGroup(runName);
             StorageGroup sig = run.openGroup("signal_channels");
             StorageGroup idx = run.openGroup("spectrum_index")) {
        List<String> channelNames = splitNames(
            attrStr(sig, "channel_names", ""));
        int acqMode = intAttr(run, "acquisition_mode", 0);
        int wireClass = wireFromSpectrumClass(
            attrStr(run, "spectrum_class", "TTIOMassSpectrum"));

        Map<String, List<ChannelSegment>> segsByCh = new LinkedHashMap<>();
        for (String c : channelNames) {
            segsByCh.put(c, PerAUFile.readChannelSegments(sig, c + "_segments"));
        }
        List<HeaderSegment> hdrSegs = null;
        if (headersEncrypted) {
            hdrSegs = PerAUFile.readHeaderSegments(idx, "au_header_segments");
        }

        List<ChannelSegment> firstSegs = segsByCh.get(channelNames.get(0));
        int n = firstSegs.size();

        // Plaintext filter arrays (only read when headers are NOT encrypted).
        double[] rts = null, pmzs = null, bpis = null;
        int[] msLevels = null, pols = null, pcs = null;
        if (!headersEncrypted) {
            try (StorageDataset d = idx.openDataset("retention_times")) {
                rts = (double[]) d.readAll(); }
            try (StorageDataset d = idx.openDataset("ms_levels")) {
                msLevels = (int[]) d.readAll(); }
            try (StorageDataset d = idx.openDataset("polarities")) {
                pols = (int[]) d.readAll(); }
            try (StorageDataset d = idx.openDataset("precursor_mzs")) {
                pmzs = (double[]) d.readAll(); }
            try (StorageDataset d = idx.openDataset("precursor_charges")) {
                pcs = (int[]) d.readAll(); }
            try (StorageDataset d = idx.openDataset("base_peak_intensities")) {
                bpis = (double[]) d.readAll(); }
        }

        for (int i = 0; i < n; i++) {
            List<byte[]> channelPayloads = new ArrayList<>(channelNames.size());
            for (String c : channelNames) {
                ChannelSegment seg = segsByCh.get(c).get(i);
                byte[] data = concat(seg.iv(), seg.tag(), seg.ciphertext());
                channelPayloads.add(encodeChannelData(c, seg.length(), data));
            }
            int flags = PacketHeader.FLAG_ENCRYPTED;
            byte[] payload;
            if (headersEncrypted) {
                flags |= PacketHeader.FLAG_ENCRYPTED_HEADER;
                payload = encodeEncryptedHeaderAU(wireClass, channelNames.size(),
                    hdrSegs.get(i), channelPayloads);
            } else {
                int polWire = polToWire(pols[i]);
                payload = encodePlaintextHeaderAU(wireClass, acqMode, msLevels[i],
                    polWire, rts[i], pmzs[i], pcs[i], 0.0, bpis[i],
                    channelPayloads);
            }
            writer.emitRawPacket(PacketType.ACCESS_UNIT, flags, datasetId, i,
                                    payload);
        }
        return n;
        }
    }

    // ────────────────────────────────────────────────────────── reader

    /** Materialise an encrypted transport stream into {@code outputPath}. */
    public static void readEncryptedToPath(String outputPath, byte[] streamData,
                                             String providerName)
            throws IOException {
        TransportReader reader = new TransportReader(streamData);
        List<TransportReader.PacketRecord> packets = reader.readAllPackets();
        reader.close();

        String title = "";
        String isa = "";
        List<String> features = new ArrayList<>();
        Map<Integer, DatasetAccumulator> datasets = new TreeMap<>();
        Map<Integer, ProtectionMeta> protection = new TreeMap<>();

        for (TransportReader.PacketRecord rec : packets) {
            PacketType t = rec.header.packetType;
            // Forward-compat: skip unknown packet types (their bytes
            // were already consumed by the reader's outer loop).
            if (t == null) continue;
            switch (t) {
                case STREAM_HEADER -> {
                    ParsedStreamHeader h = parseStreamHeader(rec.payload);
                    title = h.title;
                    isa = h.isa;
                    features.addAll(h.features);
                }
                case PROTECTION_METADATA -> {
                    protection.put(rec.header.datasetId,
                                    parseProtection(rec.payload));
                }
                case DATASET_HEADER -> {
                    ParsedDatasetHeader h = parseDatasetHeader(rec.payload);
                    boolean isGenomic = "TTIOGenomicRead".equals(h.spectrumClass);
                    DatasetAccumulator acc = new DatasetAccumulator(
                        h.name, h.acqMode, h.spectrumClass, h.channelNames,
                        h.instrumentJson, isGenomic);
                    datasets.put(h.datasetId, acc);
                }
                case GENOMIC_RUN_SIDECAR -> {
                    DatasetAccumulator acc = datasets.get(rec.header.datasetId);
                    if (acc == null) {
                        throw new IllegalStateException(
                            "GenomicRunSidecar for unknown dataset_id "
                            + rec.header.datasetId);
                    }
                    acc.runSidecar = parseRunSidecar(rec.payload);
                }
                case BLOCK_SIDECAR -> {
                    DatasetAccumulator acc = datasets.get(rec.header.datasetId);
                    if (acc == null) {
                        throw new IllegalStateException(
                            "BlockSidecar for unknown dataset_id "
                            + rec.header.datasetId);
                    }
                    acc.blockSidecars.add(parseBlockSidecar(rec.payload));
                }
                case ACCESS_UNIT -> {
                    DatasetAccumulator acc = datasets.get(rec.header.datasetId);
                    if (acc == null) {
                        throw new IllegalStateException(
                            "AU for unknown dataset_id " + rec.header.datasetId);
                    }
                    ingestAU(rec.header.flags, rec.payload, acc);
                }
                default -> { /* END_OF_DATASET / END_OF_STREAM: ignore */ }
            }
        }

        TreeSet<String> featureSet = new TreeSet<>(features);
        // The blocks_v1 wire token is transport-scoped and never a
        // container feature flag.
        featureSet.remove(PacketType.TRANSPORT_BLOCKS_V1_FEATURE);
        featureSet.remove(PacketType.TRANSPORT_BLOCKS_V1_GROUPED_FEATURE);
        featureSet.add(FeatureFlags.OPT_PER_AU_ENCRYPTION);
        boolean anyHeaderEnc = datasets.values().stream()
            .anyMatch(d -> d.usedEncryptedHeaders);
        if (anyHeaderEnc) featureSet.add(FeatureFlags.OPT_ENCRYPTED_AU_HEADERS);

        writeEncryptedFile(outputPath, providerName, title, isa,
                            new ArrayList<>(featureSet), protection, datasets);
    }

    // ──────────────────────────────────── wrapped-DEK stamp / read (PQC)

    /** Stamp a wrapped DEK + KEK algorithm onto every run's
     *  {@code signal_channels} so {@link #writeEncryptedDataset} emits it
     *  in the ProtectionMetadata packet.
     *
     *  <p>Used for envelope / PQC per-AU upload, where the per-run DEK is
     *  wrapped under the recipient's KEK ({@code aes-256-gcm}) or KEM
     *  public key ({@code ml-kem-1024}). All channels in a run share the
     *  DEK (v1.0 design), so the wrapper is stamped on each run's first
     *  channel -- the same attribute the writer reads.</p>
     *
     *  <p>Cross-language equivalent: Python
     *  {@code transport.encrypted.stamp_transport_wrapped_dek}.</p> */
    public static void stampTransportWrappedDek(String ttioPath,
                                                  byte[] wrappedDek,
                                                  String kekAlgorithm,
                                                  String providerName)
            throws IOException {
        stampTransportWrappedDek(ttioPath, wrappedDek, kekAlgorithm,
            java.util.List.of(), providerName);
    }

    /** Multi-recipient overload (FD-1 Phase A): the same DEK wrapped under
     *  ``additionalRecipients``' keys is stamped as the
     *  {@code <channel>_wrapped_dek_recipients} attribute and emitted as
     *  the packet's trailing block. Empty list keeps single-recipient runs
     *  byte-identical to pre-Phase-A. */
    public static void stampTransportWrappedDek(String ttioPath,
                                                  byte[] wrappedDek,
                                                  String kekAlgorithm,
                                                  List<Recipient> additionalRecipients,
                                                  String providerName)
            throws IOException {
        stampTransportWrappedDek(ttioPath, wrappedDek, kekAlgorithm,
            additionalRecipients, null, providerName);
    }

    /** FD-1 C-2a overload: also stamps {@code <channel>_server_kek_id} (the
     *  opaque label naming the KEK under which the primary wrapped_dek is
     *  wrapped), so {@code writeEncryptedDataset} emits it in the packet and
     *  the daemon can decide server-processability. {@code serverKekId} null
     *  marks the container BYOK / not server-processable. */
    public static void stampTransportWrappedDek(String ttioPath,
                                                  byte[] wrappedDek,
                                                  String kekAlgorithm,
                                                  List<Recipient> additionalRecipients,
                                                  String serverKekId,
                                                  String providerName)
            throws IOException {
        byte[] recipientBlock = encodeRecipientBlock(additionalRecipients);
        StorageProvider sp = ProviderRegistry.open(ttioPath,
            StorageProvider.Mode.READ_WRITE, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            if (!root.hasChild("study")) return;
            try (StorageGroup study = root.openGroup("study")) {
                for (String parent : new String[]{"ms_runs", "genomic_runs"}) {
                    if (!study.hasChild(parent)) continue;
                    try (StorageGroup g = study.openGroup(parent)) {
                        for (String n : g.childNames()) {
                            if (n.startsWith("_") || !g.hasChild(n)) continue;
                            try (StorageGroup run = g.openGroup(n)) {
                                if (!run.hasChild("signal_channels")) continue;
                                try (StorageGroup sig =
                                         run.openGroup("signal_channels")) {
                                    List<String> names = splitNames(
                                        attrStr(sig, "channel_names", ""));
                                    if (names.isEmpty()) continue;
                                    String fc = names.get(0);
                                    setWrappedDekAttr(sig, fc + "_wrapped_dek",
                                                       wrappedDek);
                                    sig.setAttribute(fc + "_kek_algorithm",
                                                      kekAlgorithm);
                                    setWrappedDekAttr(sig,
                                        fc + "_wrapped_dek_recipients",
                                        recipientBlock);
                                    if (serverKekId != null && !serverKekId.isEmpty()) {
                                        sig.setAttribute(fc + "_server_kek_id",
                                                          serverKekId);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } finally {
            sp.close();
        }
    }

    /** FD-1 C-2a: read the {@code server_kek_id} stamped on the first run's
     *  signal_channels (the KEK the daemon resolves to process the container
     *  server-side), or null if absent (BYOK / not server-processable). */
    public static String readTransportServerKekId(String ttioPath,
                                                  String providerName)
            throws IOException {
        StorageProvider sp = ProviderRegistry.open(ttioPath,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            if (!root.hasChild("study")) return null;
            try (StorageGroup study = root.openGroup("study")) {
                for (String parent : new String[]{"ms_runs", "genomic_runs"}) {
                    if (!study.hasChild(parent)) continue;
                    try (StorageGroup g = study.openGroup(parent)) {
                        for (String n : g.childNames()) {
                            if (n.startsWith("_") || !g.hasChild(n)) continue;
                            try (StorageGroup run = g.openGroup(n)) {
                                if (!run.hasChild("signal_channels")) continue;
                                try (StorageGroup sig =
                                         run.openGroup("signal_channels")) {
                                    List<String> names = splitNames(
                                        attrStr(sig, "channel_names", ""));
                                    if (names.isEmpty()) continue;
                                    String v = attrStr(sig,
                                        names.get(0) + "_server_kek_id", "");
                                    return v.isEmpty() ? null : v;
                                }
                            }
                        }
                    }
                }
            }
            return null;
        } finally {
            sp.close();
        }
    }

    /** A wrapped DEK + the KEK algorithm that unwraps it. */
    public record WrappedDek(byte[] wrappedDek, String kekAlgorithm) {}

    /** Return the wrapped DEK + KEK algorithm stamped on the first run's
     *  {@code signal_channels}, or an empty {@code WrappedDek} if none.
     *  The receiver unwraps the DEK with its KEK / KEM private key, then
     *  calls {@link PerAUFile#decryptFile}.
     *
     *  <p>Cross-language equivalent: Python
     *  {@code transport.encrypted.read_transport_wrapped_dek}.</p> */
    public static WrappedDek readTransportWrappedDek(String ttioPath,
                                                       String providerName)
            throws IOException {
        StorageProvider sp = ProviderRegistry.open(ttioPath,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            if (!root.hasChild("study")) return new WrappedDek(new byte[0], "");
            try (StorageGroup study = root.openGroup("study")) {
                for (String parent : new String[]{"ms_runs", "genomic_runs"}) {
                    if (!study.hasChild(parent)) continue;
                    try (StorageGroup g = study.openGroup(parent)) {
                        for (String n : g.childNames()) {
                            if (n.startsWith("_") || !g.hasChild(n)) continue;
                            try (StorageGroup run = g.openGroup(n)) {
                                if (!run.hasChild("signal_channels")) continue;
                                try (StorageGroup sig =
                                         run.openGroup("signal_channels")) {
                                    List<String> names = splitNames(
                                        attrStr(sig, "channel_names", ""));
                                    if (names.isEmpty()) continue;
                                    String fc = names.get(0);
                                    if (!sig.hasAttribute(fc + "_wrapped_dek")) {
                                        continue;
                                    }
                                    return new WrappedDek(
                                        readWrappedDekAttr(sig, fc + "_wrapped_dek"),
                                        attrStr(sig, fc + "_kek_algorithm", ""));
                                }
                            }
                        }
                    }
                }
            }
            return new WrappedDek(new byte[0], "");
        } finally {
            sp.close();
        }
    }

    /** Return the full DEK recipient list stamped on the first run's
     *  {@code signal_channels} (FD-1 Phase A): index 0 is the primary
     *  (recipient id {@code ""}), followed by any additional recipients.
     *  Empty list when no wrapped DEK is present. {@link
     *  #readTransportWrappedDek} stays the single-recipient accessor. */
    public static List<Recipient> readTransportRecipients(String ttioPath,
                                                            String providerName)
            throws IOException {
        StorageProvider sp = ProviderRegistry.open(ttioPath,
            StorageProvider.Mode.READ, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            if (!root.hasChild("study")) return new ArrayList<>();
            try (StorageGroup study = root.openGroup("study")) {
                for (String parent : new String[]{"ms_runs", "genomic_runs"}) {
                    if (!study.hasChild(parent)) continue;
                    try (StorageGroup g = study.openGroup(parent)) {
                        for (String n : g.childNames()) {
                            if (n.startsWith("_") || !g.hasChild(n)) continue;
                            try (StorageGroup run = g.openGroup(n)) {
                                if (!run.hasChild("signal_channels")) continue;
                                try (StorageGroup sig =
                                         run.openGroup("signal_channels")) {
                                    List<String> names = splitNames(
                                        attrStr(sig, "channel_names", ""));
                                    if (names.isEmpty()
                                        || !sig.hasAttribute(
                                            names.get(0) + "_wrapped_dek")) {
                                        continue;
                                    }
                                    String fc = names.get(0);
                                    List<Recipient> out = new ArrayList<>();
                                    out.add(new Recipient("",
                                        attrStr(sig, fc + "_kek_algorithm", ""),
                                        readWrappedDekAttr(sig, fc + "_wrapped_dek")));
                                    out.addAll(decodeRecipientBlockBytes(
                                        readWrappedDekAttr(sig,
                                            fc + "_wrapped_dek_recipients")));
                                    return out;
                                }
                            }
                        }
                    }
                }
            }
            return new ArrayList<>();
        } finally {
            sp.close();
        }
    }

    // ────────────────────────────────────────────── payload encoders

    // Package-private (not private) so the FD-1 Phase A-4 conformance test
    // can pin its output against the shared golden vectors.
    static byte[] encodeProtection(String cipherSuite, String kek,
                                             byte[] wrapped,
                                             List<Recipient> additional) {
        return encodeProtection(cipherSuite, kek, wrapped, additional, null);
    }

    static byte[] encodeProtection(String cipherSuite, String kek,
                                             byte[] wrapped,
                                             List<Recipient> additional,
                                             String serverKekId) {
        byte[] cs = cipherSuite.getBytes(StandardCharsets.UTF_8);
        byte[] kk = kek.getBytes(StandardCharsets.UTF_8);
        byte[] trailing = encodeProtectionTrailing(additional, serverKekId);
        int len = 2 + cs.length + 2 + kk.length + 4 + wrapped.length
                + 2 + 4 + trailing.length;
        ByteBuffer buf = ByteBuffer.allocate(len).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) cs.length); buf.put(cs);
        buf.putShort((short) kk.length); buf.put(kk);
        buf.putInt(wrapped.length); buf.put(wrapped);
        buf.putShort((short) 0); // signature_algorithm
        buf.putInt(0);            // public_key length
        buf.put(trailing);
        return buf.array();
    }

    /** FD-1 C-2a trailing section after the five §4.4 fields. Emitted ONLY
     *  when there are additional recipients OR a {@code serverKekId}, so
     *  pure BYOK / single-recipient packets stay byte-identical to §4.4. A
     *  single-recipient server-processable container emits
     *  {@code additional_recipient_count = 0} followed by server_kek_id. */
    private static byte[] encodeProtectionTrailing(List<Recipient> additional,
                                                   String serverKekId) {
        boolean haveAdd = additional != null && !additional.isEmpty();
        boolean haveKid = serverKekId != null && !serverKekId.isEmpty();
        if (!haveAdd && !haveKid) return new byte[0];
        byte[] block = haveAdd ? encodeRecipientBlock(additional)
                               : new byte[]{0, 0};   // additional_recipient_count = 0 (LE u16)
        byte[] kid = haveKid ? serverKekId.getBytes(StandardCharsets.UTF_8)
                             : new byte[0];
        int len = block.length + (haveKid ? 2 + kid.length : 0);
        ByteBuffer buf = ByteBuffer.allocate(len).order(ByteOrder.LITTLE_ENDIAN);
        buf.put(block);
        if (haveKid) { buf.putShort((short) kid.length); buf.put(kid); }
        return buf.array();
    }

    /** FD-1 Phase A append-only trailing block:
     *  {@code additional_recipient_count}(u16) + per recipient
     *  {@code {recipient_id, kek_algorithm, wrapped_dek}}. Empty bytes for
     *  no additional recipients, so single-recipient packets stay
     *  byte-identical to transport-spec §4.4. Byte-for-byte identical to
     *  the Python `_encode_recipient_block` (the cross-language contract).
     *  Package-private for the FD-1 Phase A-4 conformance test. */
    static byte[] encodeRecipientBlock(List<Recipient> additional) {
        if (additional == null || additional.isEmpty()) return new byte[0];
        int len = 2;
        for (Recipient r : additional) {
            len += 2 + r.recipientId().getBytes(StandardCharsets.UTF_8).length
                 + 2 + r.kekAlgorithm().getBytes(StandardCharsets.UTF_8).length
                 + 4 + r.wrappedDek().length;
        }
        ByteBuffer buf = ByteBuffer.allocate(len).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) additional.size());
        for (Recipient r : additional) {
            byte[] rid = r.recipientId().getBytes(StandardCharsets.UTF_8);
            byte[] ka = r.kekAlgorithm().getBytes(StandardCharsets.UTF_8);
            buf.putShort((short) rid.length); buf.put(rid);
            buf.putShort((short) ka.length); buf.put(ka);
            buf.putInt(r.wrappedDek().length); buf.put(r.wrappedDek());
        }
        return buf.array();
    }

    /** Decode the trailing recipient block from a buffer positioned at its
     *  start; returns empty when fewer than 2 bytes remain. */
    private static List<Recipient> decodeRecipientBlock(ByteBuffer bb) {
        List<Recipient> out = new ArrayList<>();
        if (bb.remaining() < 2) return out;
        int n = bb.getShort() & 0xFFFF;
        for (int i = 0; i < n; i++) {
            String rid = readLEString(bb, 2);
            String ka = readLEString(bb, 2);
            int wl = bb.getInt();
            byte[] wd = new byte[wl]; bb.get(wd);
            out.add(new Recipient(rid, ka, wd));
        }
        return out;
    }

    static List<Recipient> decodeRecipientBlockBytes(byte[] block) {
        if (block == null || block.length == 0) return new ArrayList<>();
        return decodeRecipientBlock(
            ByteBuffer.wrap(block).order(ByteOrder.LITTLE_ENDIAN));
    }

    private static byte[] encodeChannelData(String cname, int nElements,
                                              byte[] data) {
        byte[] nameUtf = cname.getBytes(StandardCharsets.UTF_8);
        int len = 2 + nameUtf.length + 1 + 1 + 4 + 4 + data.length;
        ByteBuffer buf = ByteBuffer.allocate(len).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort((short) nameUtf.length); buf.put(nameUtf);
        buf.put((byte) Enums.Precision.FLOAT64.ordinal());
        buf.put((byte) Enums.Compression.NONE.ordinal());
        buf.putInt(nElements);
        buf.putInt(data.length);
        buf.put(data);
        return buf.array();
    }

    private static byte[] encodeEncryptedHeaderAU(int wireClass, int nChannels,
                                                    HeaderSegment hdrSeg,
                                                    List<byte[]> channelPayloads) {
        int total = 1 + 1 + 12 + 16 + 36;
        for (byte[] cp : channelPayloads) total += cp.length;
        ByteBuffer buf = ByteBuffer.allocate(total).order(ByteOrder.LITTLE_ENDIAN);
        buf.put((byte) wireClass);
        buf.put((byte) nChannels);
        buf.put(hdrSeg.iv());
        buf.put(hdrSeg.tag());
        buf.put(hdrSeg.ciphertext());
        for (byte[] cp : channelPayloads) buf.put(cp);
        return buf.array();
    }

    private static byte[] encodePlaintextHeaderAU(int wireClass, int acqMode,
                                                    int msLevel, int polWire,
                                                    double rt, double pmz,
                                                    int pc, double ionMob,
                                                    double bpi,
                                                    List<byte[]> channelPayloads) {
        int total = 1 + 1 + 1 + 1 + 8 + 8 + 1 + 8 + 8 + 1;
        for (byte[] cp : channelPayloads) total += cp.length;
        ByteBuffer buf = ByteBuffer.allocate(total).order(ByteOrder.LITTLE_ENDIAN);
        buf.put((byte) wireClass);
        buf.put((byte) acqMode);
        buf.put((byte) (msLevel & 0xFF));
        buf.put((byte) polWire);
        buf.putDouble(rt);
        buf.putDouble(pmz);
        buf.put((byte) (pc & 0xFF));
        buf.putDouble(ionMob);
        buf.putDouble(bpi);
        buf.put((byte) channelPayloads.size());
        for (byte[] cp : channelPayloads) buf.put(cp);
        return buf.array();
    }

    // ────────────────────────────────────────────── payload decoders

    private static final class ParsedStreamHeader {
        String title = "", isa = "";
        List<String> features = new ArrayList<>();
    }

    private static ParsedStreamHeader parseStreamHeader(byte[] payload) {
        ParsedStreamHeader out = new ParsedStreamHeader();
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        readLEString(bb, 2);   // format_version — ignored here
        out.title = readLEString(bb, 2);
        out.isa = readLEString(bb, 2);
        int nFeat = bb.getShort() & 0xFFFF;
        for (int i = 0; i < nFeat; i++) out.features.add(readLEString(bb, 2));
        // n_datasets follows; discovered via DATASET_HEADER packets.
        return out;
    }

    private static final class ParsedDatasetHeader {
        int datasetId;
        String name;
        int acqMode;
        String spectrumClass;
        List<String> channelNames;
        String instrumentJson = "";
    }

    private static ParsedDatasetHeader parseDatasetHeader(byte[] payload) {
        ParsedDatasetHeader out = new ParsedDatasetHeader();
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        out.datasetId = bb.getShort() & 0xFFFF;
        out.name = readLEString(bb, 2);
        out.acqMode = bb.get() & 0xFF;
        out.spectrumClass = readLEString(bb, 2);
        int nch = bb.get() & 0xFF;
        out.channelNames = new ArrayList<>(nch);
        for (int i = 0; i < nch; i++) out.channelNames.add(readLEString(bb, 2));
        // instrument_json carries the genomic-run metadata for the
        // reader to rebuild modality/platform/reference_uri/sample_name.
        if (bb.remaining() >= 4) {
            out.instrumentJson = readLEString(bb, 4);
        } else {
            out.instrumentJson = "";
        }
        return out;
    }

    static ProtectionMeta parseProtection(byte[] payload) {
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);
        String cs = readLEString(bb, 2);
        String kek = readLEString(bb, 2);
        int wLen = bb.getInt();
        byte[] wrapped = new byte[wLen];
        bb.get(wrapped);
        // Consume signature_algorithm + public_key to reach the FD-1
        // trailing recipient block (both empty on pre-Phase-A packets, in
        // which case no bytes remain -> no additional recipients).
        List<Recipient> additional = new ArrayList<>();
        String serverKekId = null;
        if (bb.remaining() > 0) {
            readLEString(bb, 2);                     // signature_algorithm
            int pkLen = bb.getInt();
            bb.position(bb.position() + pkLen);      // public_key
            additional = decodeRecipientBlock(bb);
            // FD-1 C-2a: anything after the recipient block is server_kek_id.
            if (bb.remaining() > 0) {
                serverKekId = readLEString(bb, 2);
            }
        }
        return new ProtectionMeta(cs, kek, wrapped, additional, serverKekId);
    }

    private static void ingestAU(int flags, byte[] payload,
                                   DatasetAccumulator acc) {
        boolean encHeader = (flags & PacketHeader.FLAG_ENCRYPTED_HEADER) != 0;
        boolean encChannel = (flags & PacketHeader.FLAG_ENCRYPTED) != 0;
        if (!encChannel) {
            throw new IllegalStateException("ingestAU on plaintext AU");
        }
        acc.usedEncryptedHeaders = encHeader;
        ByteBuffer bb = ByteBuffer.wrap(payload).order(ByteOrder.LITTLE_ENDIAN);

        if (encHeader) {
            bb.get();
            int nChannels = bb.get() & 0xFF;
            byte[] iv = new byte[12]; bb.get(iv);
            byte[] tag = new byte[16]; bb.get(tag);
            byte[] ct = new byte[36]; bb.get(ct);
            acc.headerSegments.add(new HeaderSegment(iv, tag, ct));
            readEncryptedChannels(bb, nChannels, acc);
            return;
        }
        bb.get();
        int acq = bb.get() & 0xFF;
        int msLevel = bb.get() & 0xFF;
        int polWire = bb.get() & 0xFF;
        double rt = bb.getDouble();
        double pmz = bb.getDouble();
        int pc = bb.get() & 0xFF;
        double ionMob = bb.getDouble();
        double bpi = bb.getDouble();
        int nChannels = bb.get() & 0xFF;
        if (acc.acquisitionMode == 0) acc.acquisitionMode = acq;
        if (acc.isGenomic) {
            readEncryptedChannels(bb, nChannels, acc);
            global.thalion.ttio.transport.AccessUnit au =
                global.thalion.ttio.transport.AccessUnit.decode(payload);
            acc.genomicChromosomes.add(au.chromosome);
            acc.genomicPositions.add(au.position);
            acc.genomicMapqs.add(au.mappingQuality);
            acc.genomicFlags.add(au.flags & 0xFFFFFFFFL);
        } else {
            acc.plaintextRts.add(rt);
            acc.plaintextMsLevels.add(msLevel);
            int polInt = 0;
            if (polWire == 0) polInt = 1;
            else if (polWire == 1) polInt = -1;
            acc.plaintextPolarities.add(polInt);
            acc.plaintextPmzs.add(pmz);
            acc.plaintextPcs.add(pc);
            acc.plaintextBpis.add(bpi);
            readEncryptedChannels(bb, nChannels, acc);
        }
        @SuppressWarnings("unused") double _im = ionMob;
    }

    private static void readEncryptedChannels(ByteBuffer bb, int nChannels,
                                                 DatasetAccumulator acc) {
        for (int c = 0; c < nChannels; c++) {
            int nameLen = bb.getShort() & 0xFFFF;
            byte[] nameBytes = new byte[nameLen];
            bb.get(nameBytes);
            String cname = new String(nameBytes, StandardCharsets.UTF_8);
            bb.get();                           // precision
            bb.get();                           // compression
            int nElements = bb.getInt();
            int dataLen = bb.getInt();
            if (dataLen < 28) {
                throw new IllegalStateException(
                    "encrypted channel data shorter than IV+TAG");
            }
            byte[] iv = new byte[12]; bb.get(iv);
            byte[] tag = new byte[16]; bb.get(tag);
            byte[] ct = new byte[dataLen - 28]; bb.get(ct);

            List<ChannelSegment> segs = acc.channelSegments
                .computeIfAbsent(cname, k -> new ArrayList<>());
            if (!acc.channelNames.contains(cname)) acc.channelNames.add(cname);
            long prior = 0;
            for (ChannelSegment s : segs) prior += s.length();
            segs.add(new ChannelSegment(prior, nElements, iv, tag, ct));
        }
    }

    // ──────────────────────────────────────────── file materialisation

    private static void writeEncryptedFile(String path, String providerName,
                                             String title, String isa,
                                             List<String> features,
                                             Map<Integer, ProtectionMeta> protection,
                                             Map<Integer, DatasetAccumulator> datasets) {
        StorageProvider sp = ProviderRegistry.open(path,
            StorageProvider.Mode.CREATE, providerName);
        try {
            StorageGroup root = sp.rootGroup();
            // bump format_version to 1.4 + add opt_genomic when any
            // genomic dataset came through the stream (matches the
            // SpectralDataset.create heuristic).
            boolean anyGenomic = false;
            for (DatasetAccumulator d : datasets.values()) {
                if (d.isGenomic) { anyGenomic = true; break; }
            }
            String formatVersion = anyGenomic ? "1.4" : "1.1";
            List<String> writtenFeatures = new ArrayList<>(features);
            if (anyGenomic && !writtenFeatures.contains(FeatureFlags.OPT_GENOMIC)) {
                writtenFeatures.add(FeatureFlags.OPT_GENOMIC);
                java.util.Collections.sort(writtenFeatures);
            }
            new FeatureFlags(formatVersion, writtenFeatures).writeTo(root);

            try (StorageGroup study = root.createGroup("study")) {
                study.setAttribute("title", title == null ? "" : title);
                study.setAttribute("isa_investigation_id",
                                     isa == null ? "" : isa);

                // split datasets into MS vs genomic.
                Map<Integer, DatasetAccumulator> msDs = new TreeMap<>();
                Map<Integer, DatasetAccumulator> gDs = new TreeMap<>();
                for (Map.Entry<Integer, DatasetAccumulator> e : datasets.entrySet()) {
                    if (e.getValue().isGenomic) gDs.put(e.getKey(), e.getValue());
                    else msDs.put(e.getKey(), e.getValue());
                }

                // Always present, even with zero MS datasets: the
                // Python and ObjC receivers write it, and per-AU
                // decrypt walkers expect it.
                try (StorageGroup msRuns = study.createGroup("ms_runs")) {
                    StringBuilder names = new StringBuilder();
                    boolean first = true;
                    for (DatasetAccumulator acc : msDs.values()) {
                        if (!first) names.append(',');
                        names.append(acc.name);
                        first = false;
                    }
                    msRuns.setAttribute("_run_names", names.toString());
                    for (Map.Entry<Integer, DatasetAccumulator> e : msDs.entrySet()) {
                        materialiseRun(msRuns, e.getValue(), protection.get(e.getKey()));
                    }
                }

                if (!gDs.isEmpty()) {
                    try (StorageGroup gRuns = study.createGroup("genomic_runs")) {
                        StringBuilder gnames = new StringBuilder();
                        boolean first = true;
                        for (DatasetAccumulator acc : gDs.values()) {
                            if (!first) gnames.append(',');
                            gnames.append(acc.name);
                            first = false;
                        }
                        gRuns.setAttribute("_run_names", gnames.toString());
                        for (Map.Entry<Integer, DatasetAccumulator> e : gDs.entrySet()) {
                            if (e.getValue().runSidecar != null) {
                                materialiseBlocksV1GenomicRun(gRuns,
                                    e.getValue(),
                                    protection.get(e.getKey()));
                            } else {
                                materialiseGenomicRun(gRuns, e.getValue(),
                                    protection.get(e.getKey()));
                            }
                        }
                    }
                }
            }
        } finally {
            sp.close();
        }
    }

    private static void materialiseRun(StorageGroup msRuns,
                                         DatasetAccumulator acc,
                                         ProtectionMeta pm) {
        try (StorageGroup run = msRuns.createGroup(acc.name)) {
        run.setAttribute("acquisition_mode", (long) acc.acquisitionMode);
        run.setAttribute("spectrum_class",
                          acc.spectrumClass == null
                            ? "TTIOMassSpectrum" : acc.spectrumClass);
        int spectrumCount = acc.headerSegments.isEmpty()
            ? acc.channelSegments.get(acc.channelNames.get(0)).size()
            : acc.headerSegments.size();
        run.setAttribute("spectrum_count", (long) spectrumCount);

        try (StorageGroup cfg = run.createGroup("instrument_config")) {
            for (String f : new String[]{"manufacturer", "model", "serial_number",
                                           "source_type", "analyzer_type",
                                           "detector_type"}) {
                cfg.setAttribute(f, "");
            }
        }

        try (StorageGroup sig = run.createGroup("signal_channels")) {
        sig.setAttribute("channel_names", String.join(",", acc.channelNames));
        for (String cname : acc.channelNames) {
            List<ChannelSegment> segs = acc.channelSegments.get(cname);
            PerAUFile.writeChannelSegments(sig, cname + "_segments", segs);
            sig.setAttribute(cname + "_algorithm",
                pm != null && pm.cipherSuite != null
                    ? pm.cipherSuite : "aes-256-gcm");
            if (pm != null && pm.wrappedDek != null && pm.wrappedDek.length > 0) {
                setWrappedDekAttr(sig, cname + "_wrapped_dek", pm.wrappedDek);
                sig.setAttribute(cname + "_kek_algorithm",
                                  pm.kekAlgorithm == null ? "" : pm.kekAlgorithm);
                // FD-1 Phase A: persist additional recipients (no-op when single).
                setWrappedDekAttr(sig, cname + "_wrapped_dek_recipients",
                                  encodeRecipientBlock(pm.additionalRecipients));
                // FD-1 C-2a: persist server_kek_id (no-op for BYOK).
                if (pm.serverKekId != null && !pm.serverKekId.isEmpty()) {
                    sig.setAttribute(cname + "_server_kek_id", pm.serverKekId);
                }
            }
        }

        try (StorageGroup idx = run.createGroup("spectrum_index")) {
            List<ChannelSegment> firstSegs = acc.channelSegments.get(acc.channelNames.get(0));
            idx.setAttribute("count", (long) firstSegs.size());
            writePrimitiveArray(idx, "offsets", Enums.Precision.INT64,
                firstSegs.stream().mapToLong(ChannelSegment::offset).toArray());
            writePrimitiveArray(idx, "lengths", Enums.Precision.UINT32,
                firstSegs.stream().mapToInt(ChannelSegment::length).toArray());

            if (acc.usedEncryptedHeaders) {
                PerAUFile.writeHeaderSegments(idx, "au_header_segments",
                                                 acc.headerSegments);
            } else {
                writePrimitiveArray(idx, "retention_times",
                    Enums.Precision.FLOAT64, toDoubleArr(acc.plaintextRts));
                writePrimitiveArray(idx, "ms_levels",
                    Enums.Precision.INT32, toIntArr(acc.plaintextMsLevels));
                writePrimitiveArray(idx, "polarities",
                    Enums.Precision.INT32, toIntArr(acc.plaintextPolarities));
                writePrimitiveArray(idx, "precursor_mzs",
                    Enums.Precision.FLOAT64, toDoubleArr(acc.plaintextPmzs));
                writePrimitiveArray(idx, "precursor_charges",
                    Enums.Precision.INT32, toIntArr(acc.plaintextPcs));
                writePrimitiveArray(idx, "base_peak_intensities",
                    Enums.Precision.FLOAT64, toDoubleArr(acc.plaintextBpis));
            }
        }
        }  // run try-with-resources
        }  // sig try-with-resources
    }

    /** materialise one genomic dataset into the destination .tio.
     *  Mirrors materialiseRun for genomic_runs/ subtree: writes the run
     *  group with modality/platform/reference_uri/sample_name attrs
     *  (parsed from instrument_json), the encrypted signal_channels
     *  segments, and the genomic_index columns + chromosomes compound. */
    private static void materialiseGenomicRun(StorageGroup gRuns,
                                                DatasetAccumulator acc,
                                                ProtectionMeta pm) {
        try (StorageGroup run = gRuns.createGroup(acc.name)) {
            run.setAttribute("acquisition_mode", (long) acc.acquisitionMode);
            run.setAttribute("spectrum_class", acc.spectrumClass);
            run.setAttribute("modality", extractJsonField(acc.instrumentJson, "modality"));
            run.setAttribute("platform", extractJsonField(acc.instrumentJson, "platform"));
            run.setAttribute("reference_uri", extractJsonField(acc.instrumentJson, "reference_uri"));
            run.setAttribute("sample_name", extractJsonField(acc.instrumentJson, "sample_name"));
            // M97: absent on the container when "" on the wire.
            String rxRole = extractJsonField(acc.instrumentJson, "read_role");
            if (!rxRole.isEmpty()) run.setAttribute("read_role", rxRole);
            run.setAttribute("read_count", (long) acc.genomicChromosomes.size());

            try (StorageGroup sig = run.createGroup("signal_channels")) {
                StringBuilder cn = new StringBuilder();
                boolean first = true;
                for (String c : acc.channelNames) {
                    if (!first) cn.append(',');
                    cn.append(c);
                    first = false;
                }
                sig.setAttribute("channel_names", cn.toString());
                for (String cname : acc.channelNames) {
                    List<ChannelSegment> segs = acc.channelSegments.get(cname);
                    if (segs == null) continue;
                    PerAUFile.writeChannelSegments(sig, cname + "_segments", segs);
                    sig.setAttribute(cname + "_algorithm",
                        pm != null && pm.cipherSuite != null
                            ? pm.cipherSuite : "aes-256-gcm");
                    if (pm != null && pm.wrappedDek != null && pm.wrappedDek.length > 0) {
                        setWrappedDekAttr(sig, cname + "_wrapped_dek", pm.wrappedDek);
                        sig.setAttribute(cname + "_kek_algorithm",
                                          pm.kekAlgorithm == null ? "" : pm.kekAlgorithm);
                        setWrappedDekAttr(sig, cname + "_wrapped_dek_recipients",
                                          encodeRecipientBlock(pm.additionalRecipients));
                        if (pm.serverKekId != null && !pm.serverKekId.isEmpty()) {
                            sig.setAttribute(cname + "_server_kek_id", pm.serverKekId);
                        }
                    }
                }
            }

            try (StorageGroup idx = run.createGroup("genomic_index")) {
                int n = acc.genomicChromosomes.size();
                idx.setAttribute("count", (long) n);
                long[] offsetsArr = new long[n];
                int[] lengthsArr = new int[n];
                long[] positionsArr = new long[n];
                byte[] mqArr = new byte[n];
                int[] flagsArr = new int[n];
                List<ChannelSegment> firstSegs = acc.channelNames.isEmpty()
                    ? java.util.Collections.<ChannelSegment>emptyList()
                    : acc.channelSegments.get(acc.channelNames.get(0));
                for (int i = 0; i < n; i++) {
                    offsetsArr[i] = firstSegs != null && i < firstSegs.size() ? firstSegs.get(i).offset() : 0L;
                    lengthsArr[i] = firstSegs != null && i < firstSegs.size() ? firstSegs.get(i).length() : 0;
                    positionsArr[i] = acc.genomicPositions.get(i);
                    mqArr[i] = (byte) (acc.genomicMapqs.get(i) & 0xFF);
                    flagsArr[i] = (int) (acc.genomicFlags.get(i) & 0xFFFFFFFFL);
                }
                writePrimitiveArray(idx, "offsets", Enums.Precision.UINT64, offsetsArr);
                writePrimitiveArray(idx, "lengths", Enums.Precision.UINT32, lengthsArr);
                writePrimitiveArray(idx, "positions", Enums.Precision.INT64, positionsArr);
                writePrimitiveArray(idx, "mapping_qualities", Enums.Precision.UINT8, mqArr);
                writePrimitiveArray(idx, "flags", Enums.Precision.UINT32, flagsArr);
                writeChromosomesCompound(idx, acc.genomicChromosomes);
            }
        }
    }

    /** Materialise one blocks_v1 per-AU-encrypted genomic run from
     *  its sidecar packets and AU stream, in the shape the stream
     *  writer creates it, so decrypt-in-place restores it
     *  (transport-spec §4.24). */
    private static void materialiseBlocksV1GenomicRun(StorageGroup gRuns,
                                                       DatasetAccumulator acc,
                                                       ProtectionMeta pm) {
        RunSidecar sc = acc.runSidecar;
        try (StorageGroup run = gRuns.createGroup(acc.name)) {
            run.setAttribute("acquisition_mode", (long) acc.acquisitionMode);
            run.setAttribute("modality",
                extractJsonField(acc.instrumentJson, "modality"));
            run.setAttribute("spectrum_class", 5L);
            run.setAttribute("reference_uri",
                extractJsonField(acc.instrumentJson, "reference_uri"));
            run.setAttribute("platform",
                extractJsonField(acc.instrumentJson, "platform"));
            String rxRole = extractJsonField(acc.instrumentJson, "read_role");
            if (!rxRole.isEmpty()) run.setAttribute("read_role", rxRole);
            run.setAttribute("sample_name",
                extractJsonField(acc.instrumentJson, "sample_name"));
            run.setAttribute("read_count", sc.readCount);
            run.setAttribute("base_count", sc.baseCount);
            run.setAttribute("layout", sc.layout);
            run.setAttribute("block_policy", sc.blockPolicy);
            Object slice = sc.attrs.get("ref_diff_slice_bytes");
            if (slice instanceof Number sn && sn.longValue() != 0) {
                run.setAttribute("ref_diff_slice_bytes", sn.longValue());
            }
            Object noV5 = sc.attrs.get("opt_disable_qualities_v5");
            if (noV5 instanceof Number vn && vn.longValue() != 0) {
                run.setAttribute("opt_disable_qualities_v5", 1L);
            }
            Object md5s = sc.attrs.get("reference_md5s");
            if (md5s instanceof String ms && !ms.isEmpty()) {
                run.setAttribute("reference_md5s", ms);
            }

            List<BlockSidecar> sidecars = new ArrayList<>(acc.blockSidecars);
            sidecars.sort(java.util.Comparator.comparingInt(
                b -> b.blockIndex));
            // The index carries the tags triple when the sender's did
            // (M101): its block sidecars then list a tags entry.
            boolean withTags = false;
            for (BlockSidecar bs : sidecars) withTags |= bs.channels.containsKey("tags");
            List<String> order = new ArrayList<>(
                global.thalion.ttio.genomics.GenomicBlocks.REQUIRED_BLOCK_CHANNELS);
            if (withTags) order.add("tags");
            var fields = global.thalion.ttio.genomics.GenomicStreamWriter
                .indexFields(withTags);
            try (StorageGroup blocks = run.createGroup("blocks");
                 StorageDataset idxDs = blocks.createCompoundDataset(
                     "index", fields, 0, true, 1024)) {
                List<Object[]> rows = new ArrayList<>(sidecars.size());
                for (BlockSidecar bs : sidecars) {
                    Object[] row = new Object[fields.size()];
                    row[0] = bs.readStart;
                    row[1] = bs.nReads;
                    row[2] = bs.baseStart;
                    row[3] = bs.nBases;
                    for (String ch : order) {
                        long[] triple = bs.channels.get(ch);
                        int oc = global.thalion.ttio.genomics
                            .GenomicStreamWriter.offsetColumn(ch);
                        row[oc] = triple == null ? 0L : triple[0];
                        row[oc + 1] = triple == null ? 0L : triple[1];
                        row[global.thalion.ttio.genomics.GenomicStreamWriter
                            .codecColumn(ch)] = triple == null
                            ? 0 : (int) triple[2];
                    }
                    rows.add(row);
                }
                if (!rows.isEmpty()) idxDs.append(rows);
            }

            List<ChannelSegment> firstSegs = acc.channelNames.isEmpty()
                ? List.of()
                : acc.channelSegments.get(acc.channelNames.get(0));
            int n = firstSegs.size();
            Map<String, Integer> nameToId = new LinkedHashMap<>();
            for (int i = 0; i < sc.chromNames.size(); i++) {
                nameToId.put(sc.chromNames.get(i), i);
            }
            int[] lengthsArr = new int[n];
            long[] positionsArr = new long[n];
            byte[] mqArr = new byte[n];
            int[] flagsArr = new int[n];
            short[] cidArr = new short[n];
            for (int i = 0; i < n; i++) {
                lengthsArr[i] = firstSegs.get(i).length();
                positionsArr[i] = acc.genomicPositions.get(i);
                mqArr[i] = (byte) (acc.genomicMapqs.get(i) & 0xFF);
                flagsArr[i] = (int) (acc.genomicFlags.get(i) & 0xFFFFFFFFL);
                Integer slot = nameToId.get(acc.genomicChromosomes.get(i));
                if (slot == null) {
                    throw new IllegalStateException(
                        "AU chromosome '" + acc.genomicChromosomes.get(i)
                        + "' is not in the sidecar chromosome_names "
                        + "table");
                }
                cidArr[i] = slot.shortValue();
            }
            try (StorageGroup idx = run.createGroup("genomic_index")) {
                // Extendable, chunked and zlib-6 like the stream
                // writer's genomic_index arrays (GenomicIndex chunk).
                Object[][] arrays = {
                    {"lengths", Enums.Precision.UINT32, lengthsArr},
                    {"positions", Enums.Precision.INT64, positionsArr},
                    {"mapping_qualities", Enums.Precision.UINT8, mqArr},
                    {"flags", Enums.Precision.UINT32, flagsArr},
                    {"chromosome_ids", Enums.Precision.UINT16, cidArr},
                };
                for (Object[] a : arrays) {
                    try (StorageDataset ds = idx.createDataset(
                            (String) a[0], (Enums.Precision) a[1], 0,
                            65536, Enums.Compression.ZLIB, 6, true)) {
                        ds.append(a[2]);
                    }
                }
                if ("blocks_v1_grouped".equals(sc.layout)) {
                    // M103: the permutation back to input order, from the
                    // sidecars' slices in block order.
                    int total = 0;
                    for (BlockSidecar bs : sidecars) {
                        if (bs.inputIndex == null) {
                            throw new IllegalStateException("genomic run '" + acc.name
                                + "': blocks_v1_grouped stream without input_index "
                                + "in its BlockSidecars");
                        }
                        total += bs.inputIndex.length;
                    }
                    int[] ii = new int[total];
                    int w = 0;
                    for (BlockSidecar bs : sidecars) {
                        System.arraycopy(bs.inputIndex, 0, ii, w, bs.inputIndex.length);
                        w += bs.inputIndex.length;
                    }
                    try (StorageDataset ds = idx.createDataset(
                            "input_index", Enums.Precision.UINT32, 0,
                            65536, Enums.Compression.ZLIB, 6, true)) {
                        ds.append(ii);
                    }
                }
                writeNameTable(idx, "chromosome_names", sc.chromNames);
            }

            try (StorageGroup sig = run.createGroup("signal_channels")) {
                for (String cname : acc.channelNames) {
                    List<ChannelSegment> segs =
                        acc.channelSegments.get(cname);
                    if (segs == null) continue;
                    PerAUFile.writeChannelSegments(sig,
                        cname + "_segments", segs);
                    sig.setAttribute(cname + "_algorithm",
                        pm != null && pm.cipherSuite != null
                            ? pm.cipherSuite : "aes-256-gcm");
                    if (pm != null && pm.wrappedDek != null
                            && pm.wrappedDek.length > 0) {
                        setWrappedDekAttr(sig, cname + "_wrapped_dek",
                                          pm.wrappedDek);
                        sig.setAttribute(cname + "_kek_algorithm",
                            pm.kekAlgorithm == null ? "" : pm.kekAlgorithm);
                        setWrappedDekAttr(sig,
                            cname + "_wrapped_dek_recipients",
                            encodeRecipientBlock(pm.additionalRecipients));
                        if (pm.serverKekId != null
                                && !pm.serverKekId.isEmpty()) {
                            sig.setAttribute(cname + "_server_kek_id",
                                             pm.serverKekId);
                        }
                    }
                }

                StorageGroup mateGroup = null;
                try {
                    for (Map<String, Object> centry : sc.channels) {
                        String ch = (String) centry.get("name");
                        StorageGroup parent;
                        String dsName;
                        if ("mate_info".equals(ch)) {
                            mateGroup = sig.createGroup("mate_info");
                            parent = mateGroup;
                            dsName = "inline_v2";
                        } else {
                            parent = sig;
                            dsName = ch;
                        }
                        long codec = centry.get("compression")
                                instanceof Number cn ? cn.longValue() : 0L;
                        try (StorageDataset ds = parent.createDataset(
                                dsName, Enums.Precision.UINT8, 0,
                                global.thalion.ttio.genomics
                                    .GenomicStreamWriter.CHANNEL_CHUNK,
                                codec == 0 ? Enums.Compression.ZLIB
                                           : Enums.Compression.NONE,
                                6, true)) {
                            ds.setAttribute("compression", (int) codec);
                            Object extra = centry.get("extra_attrs");
                            if (extra instanceof Map<?, ?> em) {
                                for (Map.Entry<?, ?> e : em.entrySet()) {
                                    ds.setAttribute((String) e.getKey(),
                                                    e.getValue());
                                }
                            }
                            for (BlockSidecar bs : sidecars) {
                                byte[] data = bs.blobs.get(ch);
                                if (data != null && data.length > 0) {
                                    ds.append(data);
                                }
                            }
                        }
                    }
                    // The stream writer creates mate_info and its
                    // chrom_names table at close even when no mate
                    // blob was written; mirror that.
                    if (!sc.mateChromNames.isEmpty()) {
                        if (mateGroup == null) {
                            mateGroup = sig.createGroup("mate_info");
                        }
                        writeNameTable(mateGroup, "chrom_names",
                                       sc.mateChromNames);
                    }
                } finally {
                    if (mateGroup != null) mateGroup.close();
                }
            }
        }
    }

    private static void writeNameTable(StorageGroup parent, String name,
                                       List<String> names) {
        List<CompoundField> fields = List.of(new CompoundField(
            "name", CompoundField.Kind.VL_STRING));
        List<Object[]> rows = new ArrayList<>(names.size());
        for (String s : names) rows.add(new Object[]{ s });
        try (StorageDataset ds = parent.createCompoundDataset(
                name, fields, rows.size())) {
            ds.writeAll(rows);
        }
    }

    /** minimal JSON value extractor for the metadata JSON we
     *  emit on the wire. Looks for the literal pattern "key": "value"
     *  and returns the captured value; returns the empty string when
     *  the key is not present or the JSON cannot be parsed. */
    private static String extractJsonField(String json, String key) {
        if (json == null || json.isEmpty()) return "";
        char QT = (char) 34;
        char BS = (char) 92;
        String needle = QT + key + QT + ": " + QT;
        int start = json.indexOf(needle);
        if (start < 0) return "";
        start += needle.length();
        // Walk forward looking for an unescaped closing quote.
        StringBuilder out = new StringBuilder();
        int i = start;
        while (i < json.length()) {
            char c = json.charAt(i);
            if (c == BS && i + 1 < json.length()) {
                out.append(json.charAt(i + 1));
                i += 2;
                continue;
            }
            if (c == QT) break;
            out.append(c);
            i++;
        }
        return out.toString();
    }

    /** write a chromosomes compound dataset (single VL_STRING
     *  field {@code value}) mirroring the GenomicIndex layout.
     *  L1 (Task #82 Phase B.1, 2026-05-01): write chromosomes as
     *  {@code chromosome_ids} (uint16) + {@code chromosome_names}
     *  (compound) instead of a single VL-string compound. */
    private static void writeChromosomesCompound(StorageGroup idx,
                                                    List<String> chromosomes) {
        java.util.LinkedHashMap<String, Integer> nameToId = new java.util.LinkedHashMap<>();
        short[] ids = new short[chromosomes.size()];
        for (int i = 0; i < chromosomes.size(); i++) {
            String name = chromosomes.get(i);
            if (name == null) name = "";
            Integer slot = nameToId.get(name);
            if (slot == null) {
                if (nameToId.size() > 65535) {
                    throw new IllegalStateException(
                        "encrypted-transport: > 65,535 unique chromosome names");
                }
                slot = nameToId.size();
                nameToId.put(name, slot);
            }
            ids[i] = slot.shortValue();
        }
        try (StorageDataset cids = idx.createDataset(
                "chromosome_ids", Enums.Precision.UINT16, ids.length,
                0, Enums.Compression.NONE, 0)) {
            cids.writeAll(ids);
        }
        java.util.List<global.thalion.ttio.providers.CompoundField> fields =
            java.util.List.of(new global.thalion.ttio.providers.CompoundField(
                "name",
                global.thalion.ttio.providers.CompoundField.Kind.VL_STRING));
        java.util.List<Object[]> rows = new java.util.ArrayList<>(nameToId.size());
        for (String n : nameToId.keySet()) rows.add(new Object[]{ n });
        try (StorageDataset ds = idx.createCompoundDataset(
                "chromosome_names", fields, rows.size())) {
            ds.writeAll(rows);
        }
    }

    private static void writePrimitiveArray(StorageGroup parent, String name,
                                              Enums.Precision p, Object data) {
        int len = java.lang.reflect.Array.getLength(data);
        try (StorageDataset ds = parent.createDataset(name, p, len,
                0, Enums.Compression.NONE, 0)) {
            ds.writeAll(data);
        }
    }

    // ──────────────────────────────────────────────────────── helpers

    private static byte[] concat(byte[]... parts) {
        int n = 0;
        for (byte[] p : parts) n += p.length;
        byte[] out = new byte[n];
        int off = 0;
        for (byte[] p : parts) {
            System.arraycopy(p, 0, out, off, p.length);
            off += p.length;
        }
        return out;
    }

    private static String readLEString(ByteBuffer bb, int widthBytes) {
        int len = widthBytes == 2 ? (bb.getShort() & 0xFFFF) : bb.getInt();
        byte[] data = new byte[len];
        bb.get(data);
        return new String(data, StandardCharsets.UTF_8);
    }

    private static String attrStr(StorageGroup g, String name, String fallback) {
        if (!g.hasAttribute(name)) return fallback;
        Object v = g.getAttribute(name);
        return v == null ? fallback : v.toString();
    }

    private static byte[] attrBytes(StorageGroup g, String name) {
        if (!g.hasAttribute(name)) return new byte[0];
        Object v = g.getAttribute(name);
        if (v instanceof byte[] b) return b;
        if (v instanceof String s) return s.getBytes(StandardCharsets.UTF_8);
        return new byte[0];
    }

    // The wrapped DEK is binary (v1.2 / ML-KEM blob with embedded NULs).
    // Storage providers only support String / Number attributes (no
    // byte[]), so it is base64-encoded at the attribute boundary. The
    // byte[] flows through the rest of the code unchanged; this is the
    // Java analogue of the Python uint8-array attribute (each language
    // round-trips its own .tio, so the on-disk encodings need not match
    // — the cross-language interchange is the transport packet bytes).
    private static void setWrappedDekAttr(StorageGroup g, String name,
                                            byte[] wrapped) {
        if (wrapped == null || wrapped.length == 0) return;
        g.setAttribute(name, Base64.getEncoder().encodeToString(wrapped));
    }

    private static byte[] readWrappedDekAttr(StorageGroup g, String name) {
        if (!g.hasAttribute(name)) return new byte[0];
        Object v = g.getAttribute(name);
        if (v instanceof String s && !s.isEmpty()) {
            return Base64.getDecoder().decode(s);
        }
        if (v instanceof byte[] b) return b;  // legacy / direct byte[]
        return new byte[0];
    }

    private static int intAttr(StorageGroup g, String name, int fallback) {
        if (!g.hasAttribute(name)) return fallback;
        Object v = g.getAttribute(name);
        if (v instanceof Number n) return n.intValue();
        return fallback;
    }

    private static long firstChannelSegmentCount(StorageGroup sig, String firstCh) {
        List<ChannelSegment> segs = PerAUFile.readChannelSegments(
            sig, firstCh + "_segments");
        return segs.size();
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

    private static int wireFromSpectrumClass(String s) {
        return switch (s) {
            case "TTIONMRSpectrum" -> 1;
            case "TTIONMR2DSpectrum" -> 2;
            case "TTIOFreeInductionDecay" -> 3;
            case "TTIOMSImagePixel" -> 4;
            default -> 0;
        };
    }

    private static int polToWire(int polInt) {
        if (polInt == 1) return 0;
        if (polInt == -1) return 1;
        return 2;
    }

    private static double[] toDoubleArr(List<Double> xs) {
        double[] out = new double[xs.size()];
        for (int i = 0; i < xs.size(); i++) out[i] = xs.get(i);
        return out;
    }

    private static int[] toIntArr(List<Integer> xs) {
        int[] out = new int[xs.size()];
        for (int i = 0; i < xs.size(); i++) out[i] = xs.get(i);
        return out;
    }

    // ─────────────────────────────────────────── accumulator types

    static final class DatasetAccumulator {
        String name;
        int acquisitionMode;
        String spectrumClass;
        List<String> channelNames;
        Map<String, List<ChannelSegment>> channelSegments = new LinkedHashMap<>();
        List<HeaderSegment> headerSegments = new ArrayList<>();
        boolean usedEncryptedHeaders;
        List<Double> plaintextRts = new ArrayList<>();
        List<Integer> plaintextMsLevels = new ArrayList<>();
        List<Integer> plaintextPolarities = new ArrayList<>();
        List<Double> plaintextPmzs = new ArrayList<>();
        List<Integer> plaintextPcs = new ArrayList<>();
        List<Double> plaintextBpis = new ArrayList<>();
        // genomic accumulator fields
        boolean isGenomic;
        String instrumentJson = "";
        List<String> genomicChromosomes = new ArrayList<>();
        List<Long> genomicPositions = new ArrayList<>();
        List<Integer> genomicMapqs = new ArrayList<>();
        List<Long> genomicFlags = new ArrayList<>();
        // M99.1 blocks_v1 sidecars; runSidecar stays null for
        // legacy-shaped genomic runs.
        RunSidecar runSidecar;
        List<BlockSidecar> blockSidecars = new ArrayList<>();

        DatasetAccumulator(String name, int acqMode, String spectrumClass,
                             List<String> channelNames) {
            this(name, acqMode, spectrumClass, channelNames, "", false);
        }

        DatasetAccumulator(String name, int acqMode, String spectrumClass,
                             List<String> channelNames, String instrumentJson,
                             boolean isGenomic) {
            this.name = name;
            this.acquisitionMode = acqMode;
            this.spectrumClass = spectrumClass;
            this.channelNames = new ArrayList<>(channelNames);
            this.instrumentJson = instrumentJson == null ? "" : instrumentJson;
            this.isGenomic = isGenomic;
            for (String c : channelNames) {
                channelSegments.put(c, new ArrayList<>());
            }
        }
    }

    /** One DEK recipient: the same per-run DEK wrapped under this
     *  recipient's KEK. FD-1 Phase A multi-recipient support. */
    public record Recipient(String recipientId, String kekAlgorithm,
                            byte[] wrappedDek) {}

    record ProtectionMeta(String cipherSuite, String kekAlgorithm,
                            byte[] wrappedDek,
                            List<Recipient> additionalRecipients,
                            String serverKekId) {}
}
