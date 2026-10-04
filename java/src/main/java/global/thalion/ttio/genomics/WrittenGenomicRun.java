/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.ProvenanceRecord;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Write-side container for a single genomic run, passed to
 * {@link global.thalion.ttio.SpectralDataset#writeMinimal SpectralDataset.writeMinimal}.
 *
 * <p>Genomic analogue of {@link global.thalion.ttio.WrittenRun}. Pure
 * data — no methods beyond accessors and the canonical record
 * components.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOWrittenGenomicRun}, Python
 * {@code ttio.written_genomic_run.WrittenGenomicRun}.</p>
 *
 * @param acquisitionMode    {@link AcquisitionMode#GENOMIC_WGS} or
 *                           {@link AcquisitionMode#GENOMIC_WES}.
 * @param referenceUri       e.g. {@code "GRCh38.p14"}.
 * @param platform           e.g. {@code "ILLUMINA"}.
 * @param sampleName         e.g. {@code "NA12878"}.
 * @param positions          int64 per-read 0-based mapping positions.
 * @param mappingQualities   uint8 per-read mapping qualities.
 * @param flags              uint32 per-read SAM flags.
 * @param sequences          uint8 concatenated bases (one ASCII byte
 *                           per base in M82; base-packing deferred to
 *                           the codec milestone).
 * @param qualities          uint8 concatenated Phred quality scores.
 * @param offsets            uint64 per-read byte offset into
 *                           sequences/qualities.
 * @param lengths            uint32 per-read length in bases.
 * @param cigars             one CIGAR string per read.
 * @param readNames          one read name per read.
 * @param mateChromosomes    one mate chromosome per read; {@code ""}
 *                           if unpaired.
 * @param matePositions      int64 per-read mate position; {@code -1}
 *                           if unpaired.
 * @param templateLengths    int32 per-read template length; {@code 0}
 *                           if unpaired.
 * @param chromosomes        one chromosome per read (for the index).
 * @param signalCompression  codec applied to typed signal channels.
 *                           Defaults to {@link Compression#ZLIB}; pass
 *                           {@link Compression#NONE} to skip.
 * @param signalCodecOverrides M86: per-channel TTI-O codec opt-in.
 *                           Maps channel name (only {@code "sequences"}
 *                           and {@code "qualities"} accepted) to a
 *                           codec id (only {@link Compression#RANS_ORDER0},
 *                           {@link Compression#RANS_ORDER1}, or
 *                           {@link Compression#BASE_PACK} accepted).
 *                           Channels not in this map use the
 *                           {@link #signalCompression} HDF5-filter
 *                           dispatch path. Defaults to
 *                           {@link Map#of() empty}; never {@code null}.
 * @param provenanceRecords  Phase 1 (post-M91): per-run provenance
 *                           chain in insertion order. Defaults to
 *                           {@link List#of() empty}; never {@code null}.
 *                           Round-trips through the
 *                           {@code <run>/provenance_json} attribute on
 *                           the genomic run group, mirroring
 *                           {@link global.thalion.ttio.AcquisitionRun}'s
 *                           layout.
 * @param embedReference     M93 : when true (default), the writer
 *                           embeds {@link #referenceChromSeqs} at
 *                           {@code /study/references/<reference_uri>/}.
 *                           Set to false when the reference is supplied
 *                           externally and the file should not duplicate
 *                           the bytes. Has no effect when
 *                           {@link #referenceChromSeqs} is {@code null}.
 * @param referenceChromSeqs M93 : per-chromosome reference
 *                           sequence bytes (uppercase ACGTN), keyed by
 *                           chromosome name. Required input for the
 *                           REF_DIFF context-aware codec on the
 *                           {@code sequences} channel; absent values
 *                           cause REF_DIFF to fall back to BASE_PACK
 *                           silently (Q5b=C). May be {@code null}.
 * @param externalReferencePath M93 : optional explicit external
 *                           FASTA path that REF_DIFF readers can fall
 *                           back to when the reference is not embedded
 *                           and {@code REF_PATH} is unset. May be
 *                           {@code null}; not yet honoured by Java
 *                           writers (kept for symmetry with Python /
 *                           ObjC).
 */
public record WrittenGenomicRun(
    AcquisitionMode acquisitionMode,
    String referenceUri,
    String platform,
    String sampleName,
    long[] positions,
    byte[] mappingQualities,
    int[]  flags,
    byte[] sequences,
    byte[] qualities,
    long[] offsets,
    int[]  lengths,
    List<String> cigars,
    List<String> readNames,
    List<String> mateChromosomes,
    long[] matePositions,
    int[]  templateLengths,
    List<String> chromosomes,
    Compression signalCompression,
    Map<String, Compression> signalCodecOverrides,
    List<ProvenanceRecord> provenanceRecords,
    boolean embedReference,
    Map<String, byte[]> referenceChromSeqs,
    Path externalReferencePath,
    /** Phase 2c-T verbatim v2 codec blobs; null disables. */
    BulkV2Blobs bulkV2Blobs,
    /** Removes the V5 sequence-context strategies from the qualities
     *  auto-tune set (spec 2.4). Python:
     *  {@code opt_disable_qualities_v5}; ObjC:
     *  {@code optDisableQualitiesV5}. */
    boolean optDisableQualitiesV5,
    /** Write this run in the v1.8 whole-channel layout instead of
     *  {@code blocks_v1} (format-spec 10.12). Python:
     *  {@code opt_legacy_whole_channel}. */
    boolean optLegacyWholeChannel,
    /** M97 — assembly read role, persisted as the {@code @read_role}
     *  UTF-8 run attribute when non-null. Recognised values:
     *  {@code hifi}, {@code ont_ul}, {@code hic_r1}, {@code hic_r2},
     *  {@code parental_maternal}, {@code parental_paternal},
     *  {@code illumina_wgs}; other strings are stored unchecked.
     *  Python: {@code read_role}; ObjC: {@code readRole}. */
    String readRole,
    /** M97 — REF_DIFF_V2 slice byte budget: a slice closes before the
     *  read that would push it past this many bases (the 10,000-read
     *  cap still applies). 0 keeps the fixed-count rule. Writer policy
     *  only — the wire format and decoder are unchanged. Python:
     *  {@code ref_diff_slice_bytes}; ObjC: {@code refDiffSliceBytes}. */
    long refDiffSliceBytes,
    /** M101 — SAM optional fields per read (columns 12+, tab-joined as
     *  samtools prints them; {@code ""} for a read without tags).
     *  {@code null}, or a list whose entries are all {@code ""}, writes
     *  no {@code tags} channel. Stored with SAM_TAGS (codec id 18);
     *  MD:Z / NM:i are recomputed from the reference when the
     *  sequences channel is REF_DIFF_V2. Python: {@code tags};
     *  ObjC: {@code tags}. */
    List<String> tags,
    /** M103 — reorder the reads by sequence before blocking and store
     *  the permutation (layout {@code blocks_v1_grouped}, format-spec
     *  10.12.7). Only for runs without a reference; readers present input
     *  order. Memory is then unbounded in the run size. Python:
     *  {@code opt_group_reads}; ObjC: {@code optGroupReads}. */
    boolean optGroupReads
) {
    /** Pre-M103 canonical signature (29 components); reads are not
     *  grouped. */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath,
        BulkV2Blobs bulkV2Blobs,
        boolean optDisableQualitiesV5,
        boolean optLegacyWholeChannel,
        String readRole,
        long refDiffSliceBytes,
        List<String> tags
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
             readRole, refDiffSliceBytes, tags, false);
    }

    /** Pre-M101 canonical signature (28 components); no SAM tags. */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath,
        BulkV2Blobs bulkV2Blobs,
        boolean optDisableQualitiesV5,
        boolean optLegacyWholeChannel,
        String readRole,
        long refDiffSliceBytes
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
             readRole, refDiffSliceBytes, null);
    }

    /** Previous canonical signature (26 components); no read role, the
     *  fixed REF_DIFF_V2 slice rule. */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath,
        BulkV2Blobs bulkV2Blobs,
        boolean optDisableQualitiesV5,
        boolean optLegacyWholeChannel
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
             null, 0L);
    }

    /** Pre-M97 signature (25 components before the layout flag); the
     *  run is written as {@code blocks_v1}. */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath,
        BulkV2Blobs bulkV2Blobs,
        boolean optDisableQualitiesV5
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             bulkV2Blobs, optDisableQualitiesV5, false);
    }

    /** Signature before qualities V5 (25 components); qualities V5
     *  stays enabled. */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath,
        BulkV2Blobs bulkV2Blobs
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             bulkV2Blobs, false);
    }

    public WrittenGenomicRun {
        Objects.requireNonNull(acquisitionMode);
        Objects.requireNonNull(referenceUri);
        Objects.requireNonNull(platform);
        Objects.requireNonNull(sampleName);
        Objects.requireNonNull(signalCompression);
        Objects.requireNonNull(signalCodecOverrides,
            "signalCodecOverrides must not be null; pass Map.of() for none");
        Objects.requireNonNull(provenanceRecords,
            "provenanceRecords must not be null; pass List.of() for none");
        cigars                = List.copyOf(cigars);
        readNames             = List.copyOf(readNames);
        mateChromosomes       = List.copyOf(mateChromosomes);
        chromosomes           = List.copyOf(chromosomes);
        signalCodecOverrides  = Map.copyOf(signalCodecOverrides);
        provenanceRecords     = List.copyOf(provenanceRecords);
        if (tags != null) {
            List<String> t = new java.util.ArrayList<>(tags.size());
            for (String s : tags) t.add(s == null ? "" : s);
            tags = java.util.Collections.unmodifiableList(t);
        }
        // referenceChromSeqs / externalReferencePath are intentionally
        // not deep-copied — byte[] values are large and the writer
        // does not mutate them. The map reference is kept verbatim
        // (may be null).
    }

    /**
     * Backwards-compatible constructor (pre-M86) that defaults
     * {@link #signalCodecOverrides} to {@link Map#of() empty}. Existing
     * callers that build a run without per-channel codec overrides
     * continue to work unchanged.
     */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, Map.of(), List.of(),
             false, null, null, null, false);
    }

    /**
     * Backwards-compatible constructor (era, 19 components) that
     * defaults {@link #provenanceRecords} to {@link List#of() empty}.
     * Existing callers built before Phase 1 (post-M91) continue to
     * work unchanged.
     */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, List.of(),
             false, null, null, null, false);
    }

    /**
     * Backwards-compatible constructor (Phase 1 / post-M91, 20 components)
     * that defaults the M93 reference-related fields to their natural
     * "no reference embedded" values. Existing callers that don't yet
     * carry a reference continue to work unchanged.
     */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             false, null, null, null, false);
    }

    /**
     * Backwards-compatible constructor (M93 v1.5 era, 23 components)
     * that defaults the Phase 2c-T {@link #bulkV2Blobs} to {@code null}.
     * Callers that pre-date the 2026-05-05 verbatim-v2-blob carriage
     * additive continue to work unchanged.
     */
    public WrittenGenomicRun(
        AcquisitionMode acquisitionMode,
        String referenceUri,
        String platform,
        String sampleName,
        long[] positions,
        byte[] mappingQualities,
        int[]  flags,
        byte[] sequences,
        byte[] qualities,
        long[] offsets,
        int[]  lengths,
        List<String> cigars,
        List<String> readNames,
        List<String> mateChromosomes,
        long[] matePositions,
        int[]  templateLengths,
        List<String> chromosomes,
        Compression signalCompression,
        Map<String, Compression> signalCodecOverrides,
        List<ProvenanceRecord> provenanceRecords,
        boolean embedReference,
        Map<String, byte[]> referenceChromSeqs,
        Path externalReferencePath
    ) {
        this(acquisitionMode, referenceUri, platform, sampleName,
             positions, mappingQualities, flags, sequences, qualities,
             offsets, lengths, cigars, readNames, mateChromosomes,
             matePositions, templateLengths, chromosomes,
             signalCompression, signalCodecOverrides, provenanceRecords,
             embedReference, referenceChromSeqs, externalReferencePath,
             null);
    }

    /**
     * M93 full-fat builder. Returns a new instance with the same
     * payload but with the M93 reference fields replaced. Builder-style
     * convenience for callers that already have a base run and want to
     * attach a reference for REF_DIFF dispatch.
     */
    public WrittenGenomicRun withReference(
        boolean embed, Map<String, byte[]> chromSeqs, Path externalPath) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embed, chromSeqs, externalPath, bulkV2Blobs,
            optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Phase 2c-T builder: returns a new instance with the given
     *  verbatim v2 blobs attached. Used by the transport bulk-mode
     *  receiver to bypass the v2 codec encode in the writer. */
    public WrittenGenomicRun withBulkV2Blobs(BulkV2Blobs blobs) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Same run written in the v1.8 whole-channel layout when
     *  {@code legacy} is true, as {@code blocks_v1} otherwise. */
    public WrittenGenomicRun withOptLegacyWholeChannel(boolean legacy) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, legacy,
            readRole, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Same run with the given per-channel codec overrides. */
    public WrittenGenomicRun withSignalCodecOverrides(Map<String, Compression> overrides) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, overrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Same run with the given provenance records. */
    public WrittenGenomicRun withProvenance(List<ProvenanceRecord> records) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, records,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Same run with the given {@code @read_role} value (M97). */
    public WrittenGenomicRun withReadRole(String role) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            role, refDiffSliceBytes, tags, optGroupReads);
    }

    /** Same run with the given REF_DIFF_V2 slice byte budget (M97). */
    public WrittenGenomicRun withRefDiffSliceBytes(long sliceBytes) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, sliceBytes, tags, optGroupReads);
    }

    /** Same run with the given per-read SAM tag text (M101); {@code null}
     *  for none. */
    public WrittenGenomicRun withTags(List<String> newTags) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, newTags, optGroupReads);
    }

    /** Same run with read grouping (M103, layout
     *  {@code blocks_v1_grouped}) switched on or off. */
    public WrittenGenomicRun withOptGroupReads(boolean group) {
        return new WrittenGenomicRun(
            acquisitionMode, referenceUri, platform, sampleName,
            positions, mappingQualities, flags, sequences, qualities,
            offsets, lengths, cigars, readNames, mateChromosomes,
            matePositions, templateLengths, chromosomes,
            signalCompression, signalCodecOverrides, provenanceRecords,
            embedReference, referenceChromSeqs, externalReferencePath,
            bulkV2Blobs, optDisableQualitiesV5, optLegacyWholeChannel,
            readRole, refDiffSliceBytes, tags, group);
    }

    /** True when some read carries non-empty SAM tag text (M101): the
     *  run then writes a {@code tags} channel. */
    public boolean hasTags() {
        if (tags == null) return false;
        for (String t : tags) if (!t.isEmpty()) return true;
        return false;
    }

    /** Number of reads (derived from {@link #offsets} length). */
    public int readCount() { return offsets.length; }
}
