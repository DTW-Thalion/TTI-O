/*
 * TTI-O Java Implementation
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.transport;

/**
 * Transport packet types. See {@code docs/transport-spec.md} §3.2.
 *
 * <p>Cross-language equivalents: Python
 * {@code ttio.transport.packets.PacketType}, Objective-C
 * {@code TTIOTransportPacketType}.</p>
 */
public enum PacketType {
    STREAM_HEADER       (0x01),
    DATASET_HEADER      (0x02),
    ACCESS_UNIT         (0x03),
    PROTECTION_METADATA (0x04),
    ANNOTATION          (0x05),
    PROVENANCE          (0x06),
    CHROMATOGRAM        (0x07),
    END_OF_DATASET      (0x08),
    /**
     * Phase 2c-T bulk mode (transport-spec §4.10): verbatim
     * {@code mate_info/inline_v2} blob carriage. Emitted only when
     * the StreamHeader features list contains
     * {@code "bulk_mode_v2_blobs"}.
     */
    BLOB_V2_MATE_INFO   (0x09),
    /** Phase 2c-T bulk mode (transport-spec §4.11): verbatim
     * {@code sequences/refdiff_v2} blob carriage. */
    BLOB_V2_REF_DIFF    (0x0A),
    /** Phase 2c-T bulk mode (transport-spec §4.12): verbatim
     * {@code read_names/name_tok_v2} blob carriage. */
    BLOB_V2_NAME_TOK    (0x0B),
    // -------- v0.11 (transport-spec-complete-coverage 2026-05-25) --------
    /** Per-reference group header. See transport-spec §4.13. */
    REFERENCE_GROUP_HEADER (0x10),
    /** Per-chromosome record inside a reference group. §4.14. */
    REFERENCE_CHROMOSOME   (0x11),
    /** Terminator of a reference group. §4.15. */
    END_OF_REFERENCE_GROUP (0x12),
    /** Image-cube header (MSImage / vibrational / UV-Vis imaging). §4.16. */
    IMAGE_HEADER           (0x13),
    /** Single pixel of an image cube. §4.17. */
    IMAGE_PIXEL            (0x14),
    /** Terminator of an image cube. §4.18. */
    END_OF_IMAGE           (0x15),
    /** Identifications table — Arrow IPC stream. §4.19. */
    IDENTIFICATIONS_TABLE  (0x16),
    /** Quantifications table — Arrow IPC stream. §4.20. */
    QUANTIFICATIONS_TABLE  (0x17),
    /** Dataset-level provenance chain. §4.21. */
    DATASET_PROVENANCE     (0x18),
    /** Subjects metadata — Arrow IPC stream. §4.22. */
    SUBJECT_METADATA       (0x19),
    /** Samples metadata — Arrow IPC stream. §4.22. */
    SAMPLE_METADATA        (0x1A),
    /** Dataset-level @encrypted algorithm name. §4.23. */
    ENCRYPTION_ALGORITHM   (0x1B),
    /** M99.1 blocks_v1 run-level sidecar. §4.24. */
    GENOMIC_RUN_SIDECAR    (0x1C),
    /** M99.1 blocks_v1 per-block sidecar. §4.24. */
    BLOCK_SIDECAR          (0x1D),
    END_OF_STREAM       (0xFF);

    /** Phase 2c-T feature flag in StreamHeader features list. */
    public static final String BULK_MODE_V2_BLOBS_FEATURE =
        "bulk_mode_v2_blobs";

    /** v0.11 feature flag in StreamHeader features list. */
    public static final String TRANSPORT_V0_11_FEATURE =
        "transport_v0_11";

    /** M99.1 feature token in StreamHeader features list when
     *  blocks_v1 sidecar packets ride on the wire. Wire-scoped:
     *  receivers strip it before writing container feature flags. */
    public static final String TRANSPORT_BLOCKS_V1_FEATURE =
        "transport_blocks_v1";

    /** Required StreamHeader token (transport-spec v0.13) when a run in
     *  the stream has the {@code blocks_v1_grouped} layout: its
     *  BlockSidecars end with the block's slice of
     *  {@code genomic_index/input_index} (M103). Wire-scoped, like
     *  {@link #TRANSPORT_BLOCKS_V1_FEATURE}. */
    public static final String TRANSPORT_BLOCKS_V1_GROUPED_FEATURE =
        "transport_blocks_v1_grouped";

    /** Phase 2c-T codec id constants (mirror enums.Compression). */
    public static final int CODEC_ID_MATE_INLINE_V2    = 13;
    public static final int CODEC_ID_REF_DIFF_V2       = 14;
    public static final int CODEC_ID_NAME_TOKENIZED_V2 = 15;

    private final int wire;
    PacketType(int wire) { this.wire = wire; }

    /** Wire byte value for this packet type. */
    public int wire() { return wire; }

    public static PacketType fromWire(int v) {
        PacketType t = fromWireOrNull(v);
        if (t == null) {
            throw new IllegalArgumentException("unknown packet type: 0x"
                    + Integer.toHexString(v));
        }
        return t;
    }

    /** Like {@link #fromWire} but returns {@code null} for unknown
     *  packet type bytes. Used by the reader's skip-unknown path so
     *  v0.10 readers can tolerate v0.11+ packet types they don't
     *  recognise (see transport-spec §6 / spec coverage v0.11). */
    public static PacketType fromWireOrNull(int v) {
        for (PacketType t : values()) if (t.wire == v) return t;
        return null;
    }
}
