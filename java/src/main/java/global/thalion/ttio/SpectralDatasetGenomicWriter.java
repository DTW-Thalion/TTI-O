/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.genomics.BulkV2Blobs;          // Phase 2c-T
import global.thalion.ttio.genomics.GenomicIndex;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.hdf5.Hdf5CompoundIO;
import global.thalion.ttio.hdf5.Hdf5Group;

import java.util.*;

/** Package-private helper extracted from {@link SpectralDataset} (P3.10).
 *  Pure code movement; no behavior change. */
public final class SpectralDatasetGenomicWriter {

    private SpectralDatasetGenomicWriter() { }

    /** SAM REVERSE flag bit (0x10). */

    static final int SAM_REVERSE_FLAG = 16;

    /** write one {@code /study/genomic_runs/<name>/}
     *  subtree via the StorageGroup protocol. Provider-agnostic — used
     *  by both the HDF5 fast path and the {@code memory://} /
     *  {@code sqlite://} / {@code zarr://} paths.
     *
     *  <p>M86: validates {@link WrittenGenomicRun#signalCodecOverrides}
     *  before any file mutation (: only sequences
     *  / qualities accept overrides). Phase A allowed RANS_ORDER0 /
     *  RANS_ORDER1 / BASE_PACK on either byte channel. Phase D
     *  () adds QUALITY_BINNED but restricts it to
     *  the {@code qualities} channel only — applying it to ACGT bytes
     *  would silently destroy the sequence via Phred-bin quantisation.
     *  Validation throws {@link IllegalArgumentException} so the caller
     *  sees the failure immediately and the file is left untouched.</p> */
    static void writeGenomicRunSubtree(
            global.thalion.ttio.providers.StorageGroup parent,
            String name,
            WrittenGenomicRun run) {
        writeGenomicRunSubtree(parent, name, run,
            global.thalion.ttio.genomics.GenomicWriteContext.none());
    }

    /** Whole-channel writer entry point; {@code ctx} carries the state a
     *  {@code blocks_v1} writer shares across blocks. */
    public static void writeGenomicRunSubtree(
            global.thalion.ttio.providers.StorageGroup parent,
            String name,
            WrittenGenomicRun run,
            global.thalion.ttio.genomics.GenomicWriteContext ctx) {
        // M86 Phase D/E/B: per-channel allowed-codec map (Gotcha §119).
        // Sequences accepts RANS+BASE_PACK; qualities additionally
        // accepts QUALITY_BINNED. Phase B adds
        // positions/flags/mapping_qualities which accept ONLY the rANS
        // codecs (— BASE_PACK / QUALITY_BINNED
        // would silently corrupt the integer values).
        // Runs BEFORE any group/dataset is created so the file is
        // untouched on a bad override (Gotcha §96 — no half-written run).
        //
        // read_names + mate_info_* + cigars
        // NAME_TOKENIZED entries removed. read_names is now v2-only
        // (no override surface — v2 is the auto-default and only path).
        // mate_info is v2-only (inline_v2 blob; the v1 per-field
        // subgroup writer is gone). cigars retains rANS only; the v1
        // NAME_TOKENIZED codec was deleted. REF_DIFF (id 9) override
        // for sequences was also deleted — REF_DIFF_V2 is the
        // auto-default when a reference is available.
        java.util.Map<String, java.util.Set<Enums.Compression>>
            allowedCodecsByChannel = java.util.Map.of(
                "sequences", java.util.Set.of(
                    Enums.Compression.RANS_ORDER0,
                    Enums.Compression.RANS_ORDER1,
                    Enums.Compression.BASE_PACK),
                "qualities", java.util.Set.of(
                    Enums.Compression.RANS_ORDER0,
                    Enums.Compression.RANS_ORDER1,
                    Enums.Compression.BASE_PACK,
                    Enums.Compression.QUALITY_BINNED,
                    // M94.Z v1.2: CRAM-mimic rANS-Nx16 quality codec.
                    Enums.Compression.FQZCOMP_NX16_Z),
                "cigars", java.util.Set.of(
                    Enums.Compression.RANS_ORDER0,
                    Enums.Compression.RANS_ORDER1));
                // positions / flags / mapping_qualities REMOVED
                // from the override surface. These per-record integer
                // fields live only in genomic_index/ now (mirroring
                // MS's spectrum_index/ pattern). The droppedIntChannels
                // set below catches the keys with a clear v1.6 error
                // pointing at genomic_index/.
                // read_names + mate_info_* keys
                // removed entirely from the override surface (rejected
                // below by the unconditional reject branches and the
                // generic "channel not supported" branch).
        // per-record integer metadata channels removed from the
        // signal_channels/ override surface. They live exclusively
        // under genomic_index/ now. Hard-error so callers with stale
        // code learn immediately.
        java.util.Set<String> droppedIntChannels = java.util.Set.of(
            "positions", "flags", "mapping_qualities");
        for (var entry : run.signalCodecOverrides().entrySet()) {
            String chName = entry.getKey();
            Enums.Compression codec = entry.getValue();
            if (droppedIntChannels.contains(chName)) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides[\"" + chName + "\"]: removed "
                    + "in v1.6 — per-record integer metadata fields "
                    + "(positions, flags, mapping_qualities) are stored "
                    + "only under genomic_index/, not signal_channels/. "
                    + "The override no longer applies. See "
                    + "docs/format-spec.md §4 and §10.7.");
            }
            // mate_info v2 (inline_v2) is the only
            // path; the v1 per-field override surface is gone. Reject
            // mate_info_* keys unconditionally — there is no longer a
            // writer code path that honours them.
            if (chName.equals("mate_info_chrom")
                    || chName.equals("mate_info_pos")
                    || chName.equals("mate_info_tlen")) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides['" + chName + "']: "
                    + "mate_info v2 (inline_v2 blob) is the only "
                    + "supported path in v1.0+; the v1 per-field "
                    + "override keys (mate_info_chrom / mate_info_pos "
                    + "/ mate_info_tlen) were removed in Phase 2c.");
            }
            // read_names is v2-only. Reject any
            // override on read_names (the v1 NAME_TOKENIZED writer
            // dispatch was removed; v2 is the auto-default and only
            // supported codec).
            if ("read_names".equals(chName)) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides['read_names']: read_names is "
                    + "v2-only in v1.0+ (NAME_TOKENIZED_V2 = 15, the "
                    + "auto-default). The v1 NAME_TOKENIZED override "
                    + "(codec id 8) was removed in Phase 2c. Build "
                    + "with the native libttio_rans library to use "
                    + "the v2 codec.");
            }
            // the bare "mate_info" key is
            // rejected; mate_info is v2-only (inline_v2 blob), with
            // no per-field override surface.
            if ("mate_info".equals(chName)) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides['mate_info']: the 'mate_info' "
                    + "key is rejected — mate_info is v2-only "
                    + "(inline_v2 blob, codec id 13) in v1.0+. The v1 "
                    + "per-field override surface (mate_info_chrom, "
                    + "mate_info_pos, mate_info_tlen) was removed in "
                    + "Phase 2c.");
            }
            if (!allowedCodecsByChannel.containsKey(chName)) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides: channel '" + chName
                    + "' not supported (in v1.0+ only sequences, "
                    + "qualities, and cigars accept overrides; "
                    + "read_names + mate_info are v2-only with no "
                    + "override surface)");
            }
            java.util.Set<Enums.Compression> allowed =
                allowedCodecsByChannel.get(chName);
            if (codec == null || !allowed.contains(codec)) {
                // Phase D : explicit message for
                // the (sequences, QUALITY_BINNED) category error —
                // names the codec, the channel, and the
                // lossy-quantisation rationale.
                if (codec == Enums.Compression.QUALITY_BINNED
                        && "sequences".equals(chName)) {
                    throw new IllegalArgumentException(
                        "signalCodecOverrides['" + chName + "']: codec "
                        + "QUALITY_BINNED is not valid on the '"
                        + chName + "' channel — quality binning is "
                        + "lossy and only applies to Phred quality "
                        + "scores. Applying it to ACGT sequence bytes "
                        + "would silently destroy the sequence via "
                        + "Phred-bin quantisation. Use the 'qualities' "
                        + "channel for QUALITY_BINNED, or RANS_ORDER0/"
                        + "RANS_ORDER1/BASE_PACK on sequences.");
                }
                // Phase C : explicit messages for
                // wrong-content codecs on the cigars channel. CIGAR
                // strings are 7-bit ASCII per the SAM spec; BASE_PACK
                // assumes ACGT bytes and QUALITY_BINNED assumes Phred
                // values, so either would silently corrupt the CIGARs.
                if ("cigars".equals(chName)) {
                    if (codec == Enums.Compression.BASE_PACK) {
                        throw new IllegalArgumentException(
                            "signalCodecOverrides['" + chName + "']: codec "
                            + "BASE_PACK is not valid on the '"
                            + chName + "' channel — BASE_PACK 2-bit-packs "
                            + "ACGT sequence bytes and would silently "
                            + "corrupt the structured ASCII strings stored "
                            + "on this channel. Use RANS_ORDER0 or "
                            + "RANS_ORDER1 on '" + chName + "'.");
                    }
                    if (codec == Enums.Compression.QUALITY_BINNED) {
                        throw new IllegalArgumentException(
                            "signalCodecOverrides['" + chName + "']: codec "
                            + "QUALITY_BINNED is not valid on the '"
                            + chName + "' channel — QUALITY_BINNED "
                            + "quantises Phred quality scores onto an "
                            + "8-bin centre table and would silently "
                            + "destroy the structured ASCII strings stored "
                            + "on this channel. Use RANS_ORDER0 or "
                            + "RANS_ORDER1 on '" + chName + "'.");
                    }
                }
                throw new IllegalArgumentException(
                    "signalCodecOverrides['" + chName + "']: codec "
                    + codec + " not supported on the '" + chName
                    + "' channel (allowed: " + allowed + ")");
            }
            // M97: QUALITY_BINNED is a fixed Illumina-8 bin table;
            // reject it on long-read platforms rather than write a
            // mistuned lossy file.
            if ("qualities".equals(chName)
                    && codec == Enums.Compression.QUALITY_BINNED
                    && !global.thalion.ttio.codecs.Quality
                        .binnedAllowedForPlatform(run.platform())) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides['qualities']: QUALITY_BINNED "
                    + "applies the fixed Illumina-8 bin table, which "
                    + "does not fit the quality distribution of "
                    + "platform '" + run.platform() + "' (M97).");
            }
        }
        try (var rg = parent.createGroup(name)) {
            // Run-level attributes.
            rg.setAttribute("acquisition_mode",
                (long) run.acquisitionMode().ordinal());
            rg.setAttribute("modality", "genomic_sequencing");
            rg.setAttribute("spectrum_class", 5L);
            rg.setAttribute("reference_uri", run.referenceUri());
            rg.setAttribute("platform", run.platform());
            if (run.readRole() != null && !run.readRole().isEmpty())
                rg.setAttribute("read_role", run.readRole());
            rg.setAttribute("sample_name", run.sampleName());
            rg.setAttribute("read_count", (long) run.readCount());

            // genomic_index (delegates to GenomicIndex.writeTo).
            GenomicIndex idx = new GenomicIndex(
                run.offsets(), run.lengths(), run.chromosomes(),
                run.positions(), run.mappingQualities(), run.flags());
            try (var ig = rg.createGroup("genomic_index")) {
                idx.writeTo(ig, ctx.chromNameToId());
            }

            // signal_channels: 5 typed channels + 3 compound datasets.
            try (var sc = rg.createGroup("signal_channels")) {
                // positions / flags / mapping_qualities are NOT
                // written under signal_channels/. They live exclusively
                // in genomic_index/, mirroring MS's spectrum_index/
                // pattern. See docs/format-spec.md §4 and §10.7.
                // Override-validation rejects these channel names.
                // sequences/qualities go through the codec
                // dispatch helper; absent from the override map →
                // existing HDF5-filter path with @compression unset.
                // ref-diff path: writeSequencesRefDiff handles both
                // the v2 fast path (when the native lib is available
                // and the run is eligible) and the BASE_PACK fallback
                // (no reference or native lib absent).
                // The path is selected when the caller has not provided
                // an explicit sequences codec, signalCompression is the
                // default ZLIB, and referenceChromSeqs is supplied.
                Enums.Compression seqCodec =
                    run.signalCodecOverrides().get("sequences");
                BulkV2Blobs bulkBlobs = run.bulkV2Blobs();
                // usesRefDiffDefaultPath captures (ZLIB default + no
                // "sequences" override); seqCodec == null is equivalent to
                // !containsKey("sequences") because signalCodecOverrides is
                // a Map.copyOf (no null values). The referenceChromSeqs gate
                // is extra to this site only.
                boolean useRefDiffPath =
                    usesRefDiffDefaultPath(run)
                    && run.referenceChromSeqs() != null;
                if (bulkBlobs != null && bulkBlobs.refDiffBlob() != null) {
                    // Phase 2c-T: skip codec encode and write the
                    // verbatim ref-diff blob.
                    if (!java.util.Objects.equals(
                            bulkBlobs.refDiffReferenceUri(),
                            run.referenceUri())) {
                        throw new IllegalArgumentException(
                            "BulkV2Blobs.refDiffReferenceUri "
                            + bulkBlobs.refDiffReferenceUri()
                            + " != run.referenceUri " + run.referenceUri());
                    }
                    writeBulkSequencesRefDiff(sc, bulkBlobs.refDiffBlob());
                } else if (useRefDiffPath) {
                    writeSequencesRefDiff(sc, run, ctx.referenceMd5());
                } else {
                    writeByteChannelWithCodec(sc, "sequences",
                        run.sequences(), run.signalCompression(),
                        seqCodec);
                }
                // M94.Z v1.2: FQZCOMP_NX16_Z is the auto-default quality
                // codec. Apply ONLY when the run is already on a v1.5
                // path (ref-diff selected for sequences OR an explicit
                // v1.5 codec override is active). This gate preserves
                // M82 byte-parity for pure-baseline writes (no
                // reference, no v1.5 overrides).
                Enums.Compression qualCodec =
                    run.signalCodecOverrides().get("qualities");
                if (qualCodec == null
                    && run.signalCompression() == Enums.Compression.ZLIB) {
                    boolean isV1_5Candidate = useRefDiffPath;
                    if (!isV1_5Candidate) {
                        for (Enums.Compression ovr
                                : run.signalCodecOverrides().values()) {
                            if (ovr == Enums.Compression.FQZCOMP_NX16_Z
                                || ovr == Enums.Compression.DELTA_RANS_ORDER0) {
                                isV1_5Candidate = true;
                                break;
                            }
                        }
                    }
                    if (isV1_5Candidate) {
                        qualCodec = Enums.Compression.FQZCOMP_NX16_Z;
                    }
                }
                if (qualCodec == Enums.Compression.FQZCOMP_NX16_Z) {
                    writeQualitiesFqzcompNx16Z(sc, run, ctx.qualStrategyHint());
                } else {
                    writeByteChannelWithCodec(sc, "qualities",
                        run.qualities(), run.signalCompression(),
                        qualCodec);
                }
                // (positions / flags / mapping_qualities removed in
                // v1.6 — see comment above and genomic_index/ writer.)

                // Compound datasets: cigars + read_names (single
                // VL_STRING). M82.4: Java now reads VL_STRING in
                // compounds correctly via Unsafe-based char* deref;
                // wire format matches Python and ObjC.
                List<global.thalion.ttio.providers.CompoundField> vlField = List.of(
                    new global.thalion.ttio.providers.CompoundField("value",
                        global.thalion.ttio.providers.CompoundField.Kind.VL_STRING));
                // schema lift for cigars. When an override
                // is present, replace the M82 compound dataset with a
                // flat 1-D uint8 dataset of the same name carrying the
                // codec output, plus an @compression attribute (Binding
                // Decisions §120-§122). only rANS
                // codecs accepted now (the v1 NAME_TOKENIZED branch was
                // removed). Two codec ids supported:
                //   * RANS_ORDER0 / RANS_ORDER1: serialise the CIGAR
                //     list[String] via length-prefix-concat
                //     (varint(asciiLen) + asciiBytes per CIGAR — §2.5
                //     of the Phase C plan), then encode through M83
                //     rANS.
                // The two layouts (override vs M82 compound) are
                // mutually exclusive within a single run; readers
                // dispatch on dataset shape and the @compression
                // attribute. No HDF5 filter is applied — codec output
                // is high-entropy ().
                if (run.signalCodecOverrides().containsKey("cigars")) {
                    Enums.Compression cigarsCodec =
                        run.signalCodecOverrides().get("cigars");
                    byte[] encoded = encodeCigars(run.cigars(), cigarsCodec);
                    global.thalion.ttio.providers.StorageDataset cgDs;
                    try {
                        cgDs = sc.createDataset("cigars",
                            Enums.Precision.UINT8, encoded.length,
                            65536, Enums.Compression.NONE, 0);
                    } catch (UnsupportedOperationException e) {
                        cgDs = sc.createDataset("cigars",
                            Enums.Precision.UINT8, encoded.length,
                            0, Enums.Compression.NONE, 0);
                    }
                    try (var closeMe = cgDs) {
                        closeMe.writeAll(encoded);
                        closeMe.setAttribute("compression",
                            codecIdFor(cigarsCodec));
                    }
                } else {
                    writeCompoundOneCol(sc, "cigars", vlField, run.cigars());
                }
                // read_names is v2-only.
                // - readCount == 0: short-circuit, write a zero-length
                //   uint8 dataset with @compression=15 so readers see a
                //   present-but-empty channel without needing the native
                //   library. Cross-language convention shared with Python
                //   and ObjC writers (uniform with the regular layout —
                //   just length 0).
                // - readCount > 0 && native lib available: encode via
                //   NameTokenizerV2 → uint8 dataset with @compression == 15.
                // - readCount > 0 && native lib unavailable: throw
                //   IllegalStateException (no fallback in v1.0+; the v1
                //   M82-compound and v1 NAME_TOKENIZED paths were removed
                //   in Phase 2c).
                if (bulkBlobs != null && bulkBlobs.nameTokBlob() != null) {
                    // Phase 2c-T: skip codec encode.
                    writeBulkReadNames(sc, bulkBlobs.nameTokBlob());
                } else if (run.readCount() == 0) {
                    global.thalion.ttio.providers.StorageDataset emptyRn;
                    try {
                        emptyRn = sc.createDataset("read_names",
                            Enums.Precision.UINT8, 0,
                            0, Enums.Compression.NONE, 0);
                    } catch (UnsupportedOperationException e) {
                        emptyRn = sc.createDataset("read_names",
                            Enums.Precision.UINT8, 0,
                            0, Enums.Compression.NONE, 0);
                    }
                    try (var closeMe = emptyRn) {
                        closeMe.setAttribute("compression",
                            codecIdFor(Enums.Compression.NAME_TOKENIZED_V2));
                    }
                } else if (!global.thalion.ttio.codecs.NameTokenizerV2.isAvailable()) {
                    throw new IllegalStateException(
                        "NAME_TOKENIZED_V2 codec requires the native "
                        + "libttio_rans library to be linked. Build with "
                        + "-Dttio.native=true or install the native "
                        + "package. (The v1 NAME_TOKENIZED codec and the "
                        + "M82 compound fallback were both removed in "
                        + "Phase 2c — there is no non-native code path "
                        + "for read_names with readCount > 0.)");
                } else {
                    byte[] encoded =
                        ((global.thalion.ttio.codecs.registry.EncodedChannel.DatasetBytes)
                            global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                                .get(Enums.Compression.NAME_TOKENIZED_V2)
                                .encode(new global.thalion.ttio.codecs.registry.DecodedChannel.StrList(
                                    run.readNames()),
                                    global.thalion.ttio.codecs.registry.CodecContext.empty())).bytes();
                    global.thalion.ttio.providers.StorageDataset rnDs;
                    try {
                        rnDs = sc.createDataset("read_names",
                            Enums.Precision.UINT8, encoded.length,
                            65536, Enums.Compression.NONE, 0);
                    } catch (UnsupportedOperationException e) {
                        rnDs = sc.createDataset("read_names",
                            Enums.Precision.UINT8, encoded.length,
                            0, Enums.Compression.NONE, 0);
                    }
                    try (var closeMe = rnDs) {
                        closeMe.writeAll(encoded);
                        closeMe.setAttribute("compression",
                            codecIdFor(Enums.Compression.NAME_TOKENIZED_V2));
                    }
                }

                // mate_info is v2-only (inline_v2
                // blob). The v1 per-field subgroup writer and the M82
                // compound fallback were removed. Empty runs OMIT the
                // mate_info group entirely (cross-language convention
                // shared with Python; ObjC was reconciled to the same).
                // Readers treat an absent group as "no mate info".
                if (bulkBlobs != null && bulkBlobs.mateInfoBlob() != null) {
                    // Phase 2c-T: skip codec encode.
                    writeBulkMateInfo(sc, bulkBlobs.mateInfoBlob(),
                        bulkBlobs.mateInfoChromNames());
                } else if (run.readCount() == 0) {
                    // Omit the group — readers handle absence as no mates.
                } else if (!global.thalion.ttio.codecs.MateInfoV2.isAvailable()) {
                    throw new IllegalStateException(
                        "MATE_INLINE_V2 codec requires the native "
                        + "libttio_rans library to be linked. Build with "
                        + "-Dttio.native=true or install the native "
                        + "package. (The v1 mate_info per-field subgroup "
                        + "and M82 compound paths were both removed in "
                        + "Phase 2c — there is no non-native code path "
                        + "for mate_info with readCount > 0.)");
                } else {
                    writeMateInfoV2(sc, run, ctx.chromNameToId());
                }

                // M101 — SAM optional fields. A run whose reads carry
                // no tags writes no tags channel, so tag-less files are
                // unchanged.
                if (run.hasTags()) {
                    writeTags(sc, run);
                }
            }

            // Phase 2 (post-M91): per-run provenance, mirroring
            // AcquisitionRun.writeProvenance. On the HDF5 fast path
            // we write the canonical compound ``provenance/steps``
            // (the same layout Python writes). The JSON attribute is
            // also emitted so non-HDF5 providers
            // (memory/sqlite/zarr) and legacy Java readers still see
            // the chain.
            if (!run.provenanceRecords().isEmpty()) {
                writeRunProvenance(rg, run.provenanceRecords());
            }
        }
    }

    /** Write {@code signal_channels/tags} through SAM_TAGS (codec id 18).
     *
     *  <p>MD:Z / NM:i are recomputed from the reference only when the
     *  sequences channel went through REF_DIFF_V2, so a reader that can
     *  decode the sequences can always decode the tags (binding decision
     *  §96). REF_DIFF_V2 runs hold one chromosome: reads on it get the
     *  reference, every other read (unmapped {@code *}) none. Mirrors
     *  Python {@code _write_tags}.</p> */
    static void writeTags(global.thalion.ttio.providers.StorageGroup sc,
                          WrittenGenomicRun run) {
        int n = run.cigars().size();
        if (run.tags().size() != n) {
            throw new IllegalArgumentException("WrittenGenomicRun.tags has "
                + run.tags().size() + " entries for " + n + " reads");
        }
        boolean derive = false;
        var seqGroup = global.thalion.ttio.genomics.GenomicBlocks.tryGroup(sc, "sequences");
        if (seqGroup != null) derive = seqGroup.hasChild("refdiff_v2");
        List<byte[]> refs = new ArrayList<>();
        short[] chromIds = new short[n];
        Arrays.fill(chromIds, (short) 0xFFFF);
        if (derive && run.referenceChromSeqs() != null) {
            TreeSet<String> names = new TreeSet<>(run.chromosomes());
            names.remove("*");
            names.remove("");
            if (names.size() == 1 && run.referenceChromSeqs().get(names.first()) != null) {
                String chrom = names.first();
                refs.add(run.referenceChromSeqs().get(chrom));
                for (int i = 0; i < n; i++) {
                    if (chrom.equals(run.chromosomes().get(i))) chromIds[i] = 0;
                }
            }
        }
        long[] offsets = new long[n + 1];
        for (int i = 0; i < n; i++) offsets[i + 1] = offsets[i] + run.lengths()[i];
        String[] cigarArr = run.cigars().toArray(new String[0]);
        var tagCtx = global.thalion.ttio.codecs.registry.CodecContext.builder()
            .sequences(refs.isEmpty() ? null : run.sequences())
            .offsets(offsets)
            .cigarsProvider(() -> cigarArr)
            .positions(run.positions())
            .ownChromIds(chromIds)
            .tagReferences(refs)
            .build();
        byte[] blob = ((global.thalion.ttio.codecs.registry.EncodedChannel.DatasetBytes)
            global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                .get(Enums.Compression.SAM_TAGS)
                .encode(new global.thalion.ttio.codecs.registry.DecodedChannel.StrList(
                    run.tags()), tagCtx)).bytes();
        global.thalion.ttio.providers.StorageDataset ds;
        try {
            ds = sc.createDataset("tags", Enums.Precision.UINT8, blob.length,
                65536, Enums.Compression.NONE, 0);
        } catch (UnsupportedOperationException e) {
            ds = sc.createDataset("tags", Enums.Precision.UINT8, blob.length,
                0, Enums.Compression.NONE, 0);
        }
        try (var closeMe = ds) {
            closeMe.writeAll(blob);
            closeMe.setAttribute("compression",
                codecIdFor(Enums.Compression.SAM_TAGS));
        }
    }

    /** Per-run provenance: the compound {@code provenance/steps} on HDF5
     *  and the {@code provenance_json} attribute on every provider. */
    public static void writeRunProvenance(
            global.thalion.ttio.providers.StorageGroup rg,
            List<ProvenanceRecord> recs) {
        try (var prov = rg.createGroup("provenance")) {
            Hdf5Group h5 = global.thalion.ttio.providers.Hdf5Provider
                .tryUnwrapHdf5Group(prov);
            if (h5 != null) {
                Hdf5CompoundIO.writeCompoundDataset(h5, "steps",
                    Hdf5CompoundIO.provenanceSchema(),
                    recs.size(),
                    (row, pool) -> {
                        ProvenanceRecord r = recs.get(row);
                        return new Object[]{
                            r.timestampUnix(),
                            pool.addString(r.software()),
                            pool.addString(r.parametersJson()),
                            pool.addString(r.inputRefsJson()),
                            pool.addString(r.outputRefsJson())
                        };
                    });
            }
        }
        rg.setAttribute("provenance_json",
            SpectralDataset.buildProvenanceJsonArray(recs));
    }

    /** Task 13 (mate_info v2): write the CRAM-style inline_v2 blob.
     *
     *  <p>Creates a subgroup {@code signal_channels/mate_info/} containing:
     *  <ul>
     *    <li>{@code inline_v2} — uint8 1-D dataset (the encoded blob),
     *        {@code @compression = 13} (MATE_INLINE_V2).</li>
     *    <li>{@code chrom_names} — compound[(name, VL_STRING)] sidecar
     *        mapping chrom_id integers (used inside the blob) back to
     *        chromosome name strings. One row per chrom in encounter
     *        order (own chroms first, then mate-only chroms). Row index
     *        == chrom_id used in the blob. {@code mate_chrom_id == -1}
     *        means unmapped (RNEXT='*'); no sidecar row for that sentinel.</li>
     *  </ul>
     *
     *  <p>Chrom-id encoding: own chromosomes are indexed by encounter
     *  order from {@code run.chromosomes()} (same as the genomic_index
     *  chromosome_ids). Mate chromosomes that reference a chrom not in
     *  the own set extend the table; {@code "*"} is mapped to -1 and
     *  never gets a sidecar row; {@code ""} (unpaired) is also -1.
     *
     *  <p>Own chrom ids come from the GenomicIndex encounter-order map
     *  and are passed as {@code short[]} (Java's closest to uint16;
     *  (short)0xFFFF for the unmapped sentinel). */
    static void writeMateInfoV2(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run) {
        writeMateInfoV2(sc, run, null);
    }

    /** {@code shared}, when given, is the run-wide chromosome id map of a
     *  {@code blocks_v1} run; it is extended in place with mate-only
     *  names so ids stay stable across blocks. */
    static void writeMateInfoV2(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run,
            java.util.Map<String, Integer> shared) {
        int n = run.readCount();

        // Build encounter-order chrom table, starting from own chroms.
        // The GenomicIndex writer uses the same encounter order; we
        // replicate it here so chrom_ids are consistent on the read side.
        java.util.Map<String, Integer> chromToId =
            shared != null ? shared : new java.util.LinkedHashMap<>();
        for (String chr : run.chromosomes()) {
            if (!chromToId.containsKey(chr)) {
                chromToId.put(chr, chromToId.size());
            }
        }
        // Extend with mate-only chroms (non-'*', non-'', non-'=').
        for (String mc : run.mateChromosomes()) {
            if (mc == null || mc.isEmpty() || "*".equals(mc) || "=".equals(mc)) continue;
            if (!chromToId.containsKey(mc)) {
                chromToId.put(mc, chromToId.size());
            }
        }

        // Build typed arrays for encode.
        short[] ownChromIds   = new short[n];
        long[]  ownPositions  = new long[n];
        int[]   mateChromIds  = new int[n];
        long[]  matePositions = new long[n];
        int[]   templateLens  = new int[n];

        for (int i = 0; i < n; i++) {
            // Own chrom id from table (unmapped own → 0xFFFF).
            String ownChr = run.chromosomes().get(i);
            Integer ownId = chromToId.get(ownChr);
            ownChromIds[i] = (ownId == null) ? (short) 0xFFFF
                           : ownId.shortValue();
            ownPositions[i] = run.positions()[i];

            // Mate chrom id: '*'/'' → -1; '=' → own chrom id.
            String mc = run.mateChromosomes().get(i);
            if (mc == null || mc.isEmpty() || "*".equals(mc)) {
                mateChromIds[i] = -1;
            } else if ("=".equals(mc)) {
                mateChromIds[i] = (ownId == null) ? -1 : ownId;
            } else {
                Integer mcId = chromToId.get(mc);
                mateChromIds[i] = (mcId == null) ? -1 : mcId;
            }
            matePositions[i] = run.matePositions()[i];
            templateLens[i]  = run.templateLengths()[i];
        }

        // Encode to the inline_v2 blob via the native JNI library.
        var mateCtx = global.thalion.ttio.codecs.registry.CodecContext.builder()
            .ownChromIds(ownChromIds)
            .ownPositions(ownPositions)
            .build();
        byte[] blob = ((global.thalion.ttio.codecs.registry.EncodedChannel.DatasetBytes)
            global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                .get(Enums.Compression.MATE_INLINE_V2)
                .encode(new global.thalion.ttio.codecs.registry.DecodedChannel.MateInfo(
                    mateChromIds, matePositions, templateLens), mateCtx)).bytes();

        // Write the mate_info group with inline_v2 + chrom_names.
        try (var mateGroup = sc.createGroup("mate_info")) {
            // Write inline_v2 blob dataset.
            global.thalion.ttio.providers.StorageDataset blobDs;
            try {
                blobDs = mateGroup.createDataset("inline_v2",
                    Enums.Precision.UINT8, blob.length,
                    65536, Enums.Compression.NONE, 0);
            } catch (UnsupportedOperationException e) {
                blobDs = mateGroup.createDataset("inline_v2",
                    Enums.Precision.UINT8, blob.length,
                    0, Enums.Compression.NONE, 0);
            }
            try (var closeMe = blobDs) {
                closeMe.writeAll(blob);
                closeMe.setAttribute("compression",
                    codecIdFor(Enums.Compression.MATE_INLINE_V2));
            }

            // Write chrom_names sidecar compound[(name, VL_STRING)].
            // One row per chrom in encounter order (row index == chrom_id).
            List<global.thalion.ttio.providers.CompoundField> nameFields = List.of(
                new global.thalion.ttio.providers.CompoundField("name",
                    global.thalion.ttio.providers.CompoundField.Kind.VL_STRING));
            List<Object[]> nameRows = new ArrayList<>(chromToId.size());
            for (String chromName : GenomicIndex.namesInIdOrder(chromToId)) {
                nameRows.add(new Object[]{ chromName });
            }
            try (var ds = mateGroup.createCompoundDataset(
                    "chrom_names", nameFields, nameRows.size())) {
                ds.writeAll(nameRows);
            }
        }
    }

    /** Phase 2c-T : write a verbatim {@code mate_info/inline_v2}
     *  blob plus the {@code chrom_names} sidecar table, BYPASSING the
     *  v2 codec encode. Used by the transport bulk-mode receiver. */
    static void writeBulkMateInfo(
            global.thalion.ttio.providers.StorageGroup sc,
            byte[] blob, List<String> chromNames) {
        try (var mateGroup = sc.createGroup("mate_info")) {
            global.thalion.ttio.providers.StorageDataset blobDs;
            try {
                blobDs = mateGroup.createDataset("inline_v2",
                    Enums.Precision.UINT8, blob.length,
                    65536, Enums.Compression.NONE, 0);
            } catch (UnsupportedOperationException e) {
                blobDs = mateGroup.createDataset("inline_v2",
                    Enums.Precision.UINT8, blob.length,
                    0, Enums.Compression.NONE, 0);
            }
            try (var closeMe = blobDs) {
                if (blob.length > 0) closeMe.writeAll(blob);
                closeMe.setAttribute("compression",
                    codecIdFor(Enums.Compression.MATE_INLINE_V2));
            }
            List<global.thalion.ttio.providers.CompoundField> nameFields = List.of(
                new global.thalion.ttio.providers.CompoundField("name",
                    global.thalion.ttio.providers.CompoundField.Kind.VL_STRING));
            List<Object[]> nameRows = new ArrayList<>(chromNames.size());
            for (String chromName : chromNames) {
                nameRows.add(new Object[]{ chromName });
            }
            try (var ds = mateGroup.createCompoundDataset(
                    "chrom_names", nameFields, nameRows.size())) {
                ds.writeAll(nameRows);
            }
        }
    }

    /** Phase 2c-T: write a verbatim {@code read_names} blob, bypassing
     *  the NameTokenizerV2 codec encode. */
    static void writeBulkReadNames(
            global.thalion.ttio.providers.StorageGroup sc, byte[] blob) {
        global.thalion.ttio.providers.StorageDataset rnDs;
        try {
            rnDs = sc.createDataset("read_names",
                Enums.Precision.UINT8, blob.length,
                65536, Enums.Compression.NONE, 0);
        } catch (UnsupportedOperationException e) {
            rnDs = sc.createDataset("read_names",
                Enums.Precision.UINT8, blob.length,
                0, Enums.Compression.NONE, 0);
        }
        try (var closeMe = rnDs) {
            if (blob.length > 0) closeMe.writeAll(blob);
            closeMe.setAttribute("compression",
                codecIdFor(Enums.Compression.NAME_TOKENIZED_V2));
        }
    }

    /** Phase 2c-T: write a verbatim {@code sequences/refdiff_v2} blob
     *  inside a {@code sequences} group, bypassing the RefDiffV2
     *  codec encode. */
    static void writeBulkSequencesRefDiff(
            global.thalion.ttio.providers.StorageGroup sc, byte[] blob) {
        try (var seqGroup = sc.createGroup("sequences")) {
            global.thalion.ttio.providers.StorageDataset blobDs;
            try {
                blobDs = seqGroup.createDataset("refdiff_v2",
                    Enums.Precision.UINT8, blob.length,
                    65536, Enums.Compression.NONE, 0);
            } catch (UnsupportedOperationException e) {
                blobDs = seqGroup.createDataset("refdiff_v2",
                    Enums.Precision.UINT8, blob.length,
                    0, Enums.Compression.NONE, 0);
            }
            try (var closeMe = blobDs) {
                if (blob.length > 0) closeMe.writeAll(blob);
                closeMe.setAttribute("compression",
                    codecIdFor(Enums.Compression.REF_DIFF_V2));
            }
        }
    }

    /** Compute the canonical reference MD5 for a run as
     *  {@code md5(concat(referenceChromSeqs[k] for k in sorted(keys)))}.
     *  Mirrors the Python {@code _reference_md5_for_run} helper. */
    public static byte[] referenceMd5ForRun(WrittenGenomicRun run) {
        if (run.referenceChromSeqs() == null) {
            return new byte[0];
        }
        if (run.referenceChromSeqs() instanceof global.thalion.ttio.genomics.LazyReference lazy) {
            return lazy.setMd5();  // cached whole-FASTA digest
        }
        try {
            java.security.MessageDigest md =
                java.security.MessageDigest.getInstance("MD5");
            java.util.List<String> sortedKeys =
                new java.util.ArrayList<>(run.referenceChromSeqs().keySet());
            java.util.Collections.sort(sortedKeys);
            for (String k : sortedKeys) {
                md.update(run.referenceChromSeqs().get(k));
            }
            return md.digest();
        } catch (java.security.NoSuchAlgorithmException e) {
            throw new RuntimeException("MD5 unavailable", e);
        }
    }

    /** True iff a written genomic run's sequences default to the ref-diff path
     *  (ZLIB default codec + no explicit "sequences" override), which embeds a
     *  reference. Single source of truth for the two former inlined copies. */
    static boolean usesRefDiffDefaultPath(WrittenGenomicRun run) {
        return run.signalCompression() == Enums.Compression.ZLIB
            && !run.signalCodecOverrides().containsKey("sequences");
    }

    /** Embed each unique reference (by {@code reference_uri}) once at
     *  {@code /study/references/<uri>/}. Only runs that have
     *  {@code embedReference=true} AND a context-aware codec on
     *  {@code sequences} (or auto-default REF_DIFF) AND non-null
     *  {@code referenceChromSeqs} contribute; the dedup key is the
     *  reference URI. When the same URI carries two different MD5s
     *  across runs, raises {@link IllegalArgumentException}. */
    public static void embedReferencesForRuns(
            global.thalion.ttio.providers.StorageGroup study,
            List<WrittenGenomicRun> genomicRuns) {
        java.util.Map<String, byte[]> needsEmbedMd5 =
            new java.util.LinkedHashMap<>();
        java.util.Map<String, java.util.Map<String, byte[]>> needsEmbedSeqs =
            new java.util.LinkedHashMap<>();
        for (WrittenGenomicRun run : genomicRuns) {
            if (!run.embedReference()) continue;
            if (run.referenceChromSeqs() == null) continue;
            // Only embed when the ref-diff path will actually be taken
            // (matches the selection condition in writeGenomicRunSubtree).
            if (!usesRefDiffDefaultPath(run)) continue;

            byte[] md5 = referenceMd5ForRun(run);
            if (needsEmbedMd5.containsKey(run.referenceUri())) {
                byte[] existing = needsEmbedMd5.get(run.referenceUri());
                if (!java.util.Arrays.equals(existing, md5)) {
                    throw new IllegalArgumentException(
                        "reference_uri '" + run.referenceUri()
                        + "' carries two different MD5s across runs in "
                        + "this dataset: " + bytesToHexLocal(existing)
                        + " vs " + bytesToHexLocal(md5)
                        + " — same URI cannot map to two different "
                        + "reference contents.");
                }
                continue;
            }
            needsEmbedMd5.put(run.referenceUri(), md5);
            needsEmbedSeqs.put(run.referenceUri(),
                new java.util.LinkedHashMap<>(run.referenceChromSeqs()));
        }
        if (needsEmbedMd5.isEmpty()) return;

        global.thalion.ttio.providers.StorageGroup refsGrp;
        if (study.hasChild("references")) {
            refsGrp = study.openGroup("references");
        } else {
            refsGrp = study.createGroup("references");
        }
        try (var ignored = refsGrp) {
            for (var entry : needsEmbedMd5.entrySet()) {
                String uri = entry.getKey();
                byte[] md5 = entry.getValue();
                java.util.Map<String, byte[]> chromSeqs =
                    needsEmbedSeqs.get(uri);
                if (refsGrp.hasChild(uri)) {
                    try (var existing = refsGrp.openGroup(uri)) {
                        Object md5Attr = existing.getAttribute("md5");
                        String existingHex = md5Attr != null
                            ? md5Attr.toString() : "";
                        if (!existingHex.equals(bytesToHexLocal(md5))) {
                            throw new IllegalArgumentException(
                                "reference_uri '" + uri + "' already "
                                + "embedded with a different MD5 ("
                                + existingHex + " != "
                                + bytesToHexLocal(md5) + ")");
                        }
                    }
                    continue;
                }
                try (var refGrp = refsGrp.createGroup(uri)) {
                    refGrp.setAttribute("md5", bytesToHexLocal(md5));
                    refGrp.setAttribute("reference_uri", uri);
                    try (var chromsGrp = refGrp.createGroup("chromosomes")) {
                        java.util.List<String> sortedNames =
                            new java.util.ArrayList<>(chromSeqs.keySet());
                        java.util.Collections.sort(sortedNames);
                        for (String chromName : sortedNames) {
                            byte[] seq = chromSeqs.get(chromName);
                            try (var c = chromsGrp.createGroup(chromName)) {
                                c.setAttribute("length", (long) seq.length);
                                // data_packed when packing wins, raw
                                // data otherwise (same dispatch as
                                // ReferenceImport.writeToDataset).
                                global.thalion.ttio.genomics.PackedReference
                                    .writeChromosomeDataset(c, seq);
                            }
                        }
                    }
                }
            }
        }
    }

    static String bytesToHexLocal(byte[] buf) {
        StringBuilder sb = new StringBuilder(buf.length * 2);
        for (byte b : buf) sb.append(String.format("%02x", b & 0xFF));
        return sb.toString();
    }

    /** Write the {@code sequences} channel through the REF_DIFF codec.
     *
     *  <p>Mirrors Python's {@code _write_sequences_ref_diff}. REF_DIFF
     *  is context-aware: it needs positions, cigars, and the reference
     *  chromosome sequence in addition to the raw byte stream.
     *
     *  <p><b>v1.0 default (REF_DIFF_V2):</b> when the native JNI library
     *  is available AND the run is eligible (single-chromosome, all
     *  reads mapped, reference present), writes
     *  {@code signal_channels/sequences} as a GROUP containing a
     *  {@code refdiff_v2} child dataset ({@code @compression = 14}).
     *
     *  <p><b>v1 path (REF_DIFF, codec id 9):</b> when the native library
     *  is unavailable or eligibility checks fail — writes
     *  {@code signal_channels/sequences} as a flat uint8 dataset with
     *  {@code @compression = 9} (or BASE_PACK = 6 on fallback).
     *
     *  <p><b>Single-chromosome limitation (v1.2 first pass):</b> all
     *  reads must align to a single chromosome. Multi-chrom is M93.X.
     *
     *  <p><b>Fallback (Q5b=C):</b> when {@code referenceChromSeqs} is
     *  null (or doesn't cover the run's chromosome) OR any read has
     *  cigar="*", falls back silently to BASE_PACK
     *  ({@code @compression = 6}). */
    static void writeSequencesRefDiff(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run) {
        writeSequencesRefDiff(sc, run, null);
    }

    /** {@code precomputedMd5}, when given, replaces the per-run digest of
     *  the reference (a {@code blocks_v1} writer computes it once). */
    static void writeSequencesRefDiff(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run,
            byte[] precomputedMd5) {
        byte[] chromSeq = null;
        if (run.referenceChromSeqs() != null) {
            java.util.Set<String> uniqueChroms =
                new java.util.LinkedHashSet<>(run.chromosomes());
            if (uniqueChroms.size() > 1) {
                throw new IllegalArgumentException(
                    "REF_DIFF v1.2 first pass supports single-chromosome "
                    + "runs only; this run carries reads on chromosomes "
                    + uniqueChroms
                    + ". Multi-chromosome support is an M93.X follow-up — "
                    + "split into per-chromosome runs as a workaround.");
            }
            if (!uniqueChroms.isEmpty()) {
                String chrom = uniqueChroms.iterator().next();
                chromSeq = run.referenceChromSeqs().get(chrom);
            }
        }

        byte[] rawBytes = run.sequences();

        // v1.0 reset Phase 2b: prefer the v2 path when eligible. Unmapped
        // reads (cigar "*") are carried by the codec since v1.9 (soft-clip
        // bases plus the slice UL substream), so they no longer force
        // BASE_PACK on the whole channel.
        boolean useV2 = global.thalion.ttio.codecs.RefDiffV2.isAvailable()
            && chromSeq != null;

        if (useV2) {
            // v1.8 path: encode via RefDiffV2 and write as a GROUP with
            // a refdiff_v2 child dataset (@compression = 14).
            byte[] md5 = precomputedMd5 != null ? precomputedMd5 : referenceMd5ForRun(run);
            int n = run.readCount();
            // Build n_reads+1 offsets from run.offsets (n entries) + total.
            long[] offsets64 = run.offsets();
            long[] offsets64n1;
            if (offsets64.length == n) {
                // run.offsets has exactly n entries; append total length.
                long totalBases = n > 0
                    ? offsets64[n - 1] + run.lengths()[n - 1]
                    : 0L;
                offsets64n1 = java.util.Arrays.copyOf(offsets64, n + 1);
                offsets64n1[n] = totalBases;
            } else if (offsets64.length == n + 1) {
                offsets64n1 = offsets64;
            } else {
                throw new IllegalArgumentException(
                    "run.offsets must have n_reads or n_reads+1 entries; "
                    + "got " + offsets64.length + " for n=" + n);
            }
            String[] cigarArr = run.cigars().toArray(new String[0]);
            long[] offsetsFinal = offsets64n1;
            var refDiffCtx = global.thalion.ttio.codecs.registry.CodecContext.builder()
                .offsets(offsetsFinal)
                .positions(run.positions())
                .reference(chromSeq)
                .referenceMd5(md5)
                .referenceUri(run.referenceUri())
                .readsPerSlice(10_000)
                .sliceBytes(run.refDiffSliceBytes() > 0
                    ? run.refDiffSliceBytes() : null)
                .cigarsProvider(() -> cigarArr)
                .build();
            var layout = (global.thalion.ttio.codecs.registry.EncodedChannel.GroupLayout)
                global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                    .get(Enums.Compression.REF_DIFF_V2)
                    .encode(new global.thalion.ttio.codecs.registry.DecodedChannel.Bytes(rawBytes),
                        refDiffCtx);
            byte[] encoded = layout.children().get("refdiff_v2");
            try (var seqGroup = sc.createGroup("sequences")) {
                global.thalion.ttio.providers.StorageDataset blobDs;
                try {
                    blobDs = seqGroup.createDataset("refdiff_v2",
                        Enums.Precision.UINT8, encoded.length,
                        65536, Enums.Compression.NONE, 0);
                } catch (UnsupportedOperationException e) {
                    blobDs = seqGroup.createDataset("refdiff_v2",
                        Enums.Precision.UINT8, encoded.length,
                        0, Enums.Compression.NONE, 0);
                }
                try (var closeMe = blobDs) {
                    closeMe.writeAll(encoded);
                    closeMe.setAttribute("compression",
                        codecIdFor(Enums.Compression.REF_DIFF_V2));
                }
            }
            return;
        }

        // the v1 REF_DIFF (codec id 9) writer
        // path was removed. When v2 cannot be used (chromSeq null,
        // unmapped reads, or native lib unavailable), fall back to
        // BASE_PACK on this channel. No v1 REF_DIFF dispatch remains.
        byte[] encoded = global.thalion.ttio.codecs.BasePack.encode(rawBytes);
        int codecId = Enums.Compression.BASE_PACK.ordinal();

        global.thalion.ttio.providers.StorageDataset ds;
        try {
            ds = sc.createDataset("sequences", Enums.Precision.UINT8,
                encoded.length, 65536, Enums.Compression.NONE, 0);
        } catch (UnsupportedOperationException e) {
            ds = sc.createDataset("sequences", Enums.Precision.UINT8,
                encoded.length, 0, Enums.Compression.NONE, 0);
        }
        try (var closeMe = ds) {
            closeMe.writeAll(encoded);
            closeMe.setAttribute("compression", codecId);
        }
    }

    /** M94.Z v1.2: write the {@code qualities} channel through the
     *  FQZCOMP_NX16_Z codec.
     *
     *  <p>Mirrors Python's {@code _write_qualities_fqzcomp_nx16_z}. The
     *  codec needs per-read {@code read_lengths} and {@code revcomp_flags},
     *  derived here from {@code run.lengths} and
     *  {@code run.flags & 16} (SAM REVERSE bit). The encoded blob is
     *  written as a flat uint8 dataset with {@code @compression = 12}. */
    static void writeQualitiesFqzcompNx16Z(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run) {
        writeQualitiesFqzcompNx16Z(sc, run, -1);
    }

    /** {@code qualStrategyHint}: -1 auto (the 3-way tune), 5/6 forced
     *  V5, {@code FqzcompNx16Z.HINT_V4_AUTO} V4 with internal preset
     *  selection; the stream writer passes its per-run pin. */
    static void writeQualitiesFqzcompNx16Z(
            global.thalion.ttio.providers.StorageGroup sc,
            WrittenGenomicRun run,
            int qualStrategyHint) {
        int n = run.readCount();
        int[] readLengths = new int[n];
        for (int i = 0; i < n; i++) readLengths[i] = run.lengths()[i];
        int[] revcompFlags = new int[n];
        for (int i = 0; i < n; i++) {
            revcompFlags[i] =
                ((run.flags()[i] & SAM_REVERSE_FLAG) != 0) ? 1 : 0;
        }
        // Qualities V5 gate (spec 2.4): offer the base bytes to the
        // encoder only when the run carries a base-parallel sequences
        // channel and the caller did not opt out; V4 still wins by
        // exact size wherever sequence context does not pay.
        byte[] v5Sequences = null;
        if (!run.optDisableQualitiesV5()
                && run.sequences() != null
                && run.sequences().length == run.qualities().length) {
            v5Sequences = run.sequences();
        }
        var ctx = global.thalion.ttio.codecs.registry.CodecContext.builder()
            .readLengths(readLengths)
            .revcompFlags(revcompFlags)
            .sequences(v5Sequences)
            .qualStrategyHint(qualStrategyHint)
            .build();
        byte[] encoded = ((global.thalion.ttio.codecs.registry.EncodedChannel.DatasetBytes)
            global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                .get(Enums.Compression.FQZCOMP_NX16_Z)
                .encode(new global.thalion.ttio.codecs.registry.DecodedChannel.Bytes(run.qualities()),
                    ctx)).bytes();
        global.thalion.ttio.providers.StorageDataset ds;
        try {
            ds = sc.createDataset("qualities", Enums.Precision.UINT8,
                encoded.length, 65536, Enums.Compression.NONE, 0);
        } catch (UnsupportedOperationException e) {
            ds = sc.createDataset("qualities", Enums.Precision.UINT8,
                encoded.length, 0, Enums.Compression.NONE, 0);
        }
        try (var closeMe = ds) {
            closeMe.writeAll(encoded);
            closeMe.setAttribute("compression",
                Enums.Compression.FQZCOMP_NX16_Z.ordinal());
        }
    }

    /** write a uint8 byte channel, optionally through a TTI-O
     *  codec (rANS order-0/1, BASE_PACK, QUALITY_BINNED). When
     *  {@code codecOverride} is {@code null} the channel is written
     *  via the default HDF5-filter path (identical to M82 behaviour,
     *  no {@code @compression} attribute set). When it names a TTI-O
     *  codec, the raw bytes are encoded, written as an unfiltered
     *  uint8 dataset (— no double-compression),
     *  and the codec id is stored on the dataset's
     *  {@code @compression} attribute (uint8).
     *
     *  <p>Phase D: QUALITY_BINNED (Phase A codec id 7) added.
     *  Caller-side validation in {@link #writeGenomicRunSubtree}
     *  guarantees this branch only fires for the {@code qualities}
     *  channel ().</p> */
    static void writeByteChannelWithCodec(
            global.thalion.ttio.providers.StorageGroup sc,
            String name, byte[] data,
            Enums.Compression defaultCodec,
            Enums.Compression codecOverride) {
        if (codecOverride == null) {
            writeSignalChannel(sc, name, Enums.Precision.UINT8, data, defaultCodec);
            return;
        }
        var codec = global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY.get(codecOverride);
        if (codec == null) {
            throw new IllegalArgumentException(
                "writeByteChannelWithCodec: unsupported codec " + codecOverride);
        }
        byte[] encoded = ((global.thalion.ttio.codecs.registry.EncodedChannel.DatasetBytes)
            codec.encode(new global.thalion.ttio.codecs.registry.DecodedChannel.Bytes(data),
                global.thalion.ttio.codecs.registry.CodecContext.empty())).bytes();
        // Unfiltered uint8 dataset; codec output already entropy-coded.
        // Force a chunked layout (chunkSize > 0) so HDF5 honours our
        // explicit Compression.NONE choice rather than the legacy
        // contiguous fallback.
        global.thalion.ttio.providers.StorageDataset ds;
        try {
            ds = sc.createDataset(name, Enums.Precision.UINT8,
                encoded.length, 65536, Enums.Compression.NONE, 0);
        } catch (UnsupportedOperationException e) {
            ds = sc.createDataset(name, Enums.Precision.UINT8,
                encoded.length, 0, Enums.Compression.NONE, 0);
        }
        try (var closeMe = ds) {
            closeMe.writeAll(encoded);
            // M79 codec id (4 / 5 / 6) as a uint8 attribute on the
            // dataset itself — read path dispatches on this.
            int codecId = codecIdFor(codecOverride);
            closeMe.setAttribute("compression", codecId);
        }
    }

    /** Map a {@link Enums.Compression} enum value to its M79 codec id
     *  (the wire-format integer that travels in the
     *  {@code @compression} attribute). The enum's {@code ordinal()}
     *  already matches the M79 numbering — this helper exists so the
     *  intent is explicit at the call site. */
    static int codecIdFor(Enums.Compression codec) {
        return codec.ordinal();
    }

    /** encode a {@code List<String>} of CIGARs through
     *  the rANS codec path.
     *
     *  <p>For {@link Enums.Compression#RANS_ORDER0} and
     *  {@link Enums.Compression#RANS_ORDER1} the list is first
     *  serialised via length-prefix-concat ({@code varint(asciiLen) +
     *  asciiBytes} per CIGAR — §2.5 of the Phase C plan), then encoded
     *  through M83 rANS at the matching order. ASCII-only per the SAM
     *  spec; non-ASCII input throws {@link IllegalArgumentException}.
     *
     *  <p>the NAME_TOKENIZED branch was removed
     *  (the v1 codec is gone). Override-validation rejects
     *  NAME_TOKENIZED on cigars upfront with a clear "no longer
     *  supported" message. */
    static byte[] encodeCigars(
            List<String> cigars,
            Enums.Compression codec) {
        if (codec != Enums.Compression.RANS_ORDER0
                && codec != Enums.Compression.RANS_ORDER1) {
            // Defensive — caller-side validation rejects this first.
            throw new IllegalArgumentException(
                "encodeCigars: unsupported codec " + codec);
        }
        java.io.ByteArrayOutputStream buf = new java.io.ByteArrayOutputStream();
        for (int idx = 0; idx < cigars.size(); idx++) {
            String cig = cigars.get(idx);
            if (cig == null) {
                throw new IllegalArgumentException(
                    "signalCodecOverrides['cigars']: cigar at index "
                    + idx + " is null");
            }
            for (int j = 0; j < cig.length(); j++) {
                if (cig.charAt(j) > 0x7F) {
                    throw new IllegalArgumentException(
                        "signalCodecOverrides['cigars']: cigar at index "
                        + idx + " contains non-ASCII bytes — CIGARs "
                        + "must be 7-bit ASCII per the SAM spec");
                }
            }
            byte[] payload = cig.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
            writeUnsignedVarint(buf, payload.length);
            buf.write(payload, 0, payload.length);
        }
        int order = (codec == Enums.Compression.RANS_ORDER0) ? 0 : 1;
        return global.thalion.ttio.codecs.Rans.encode(buf.toByteArray(), order);
    }

    /** Unsigned LEB128 varint writer — low 7 bits per byte, top bit
     *  set on continuation, terminated by the first byte with the
     *  top bit clear. The cigars rANS path serialises each entry as
     *  {@code varint(asciiLen) + asciiBytes}. */
    static void writeUnsignedVarint(
            java.io.ByteArrayOutputStream out, long n) {
        if (n < 0) {
            throw new IllegalArgumentException(
                "writeUnsignedVarint: negative value " + n);
        }
        while (Long.compareUnsigned(n, 0x80L) >= 0) {
            out.write((int) ((n & 0x7FL) | 0x80L));
            n >>>= 7;
        }
        out.write((int) (n & 0x7FL));
    }

    static void writeSignalChannel(
            global.thalion.ttio.providers.StorageGroup sc,
            String name, Enums.Precision precision, Object data,
            Enums.Compression codec) {
        int len;
        if (data instanceof long[] a)        len = a.length;
        else if (data instanceof int[] a)    len = a.length;
        else if (data instanceof byte[] a)   len = a.length;
        else throw new IllegalArgumentException(
            "writeSignalChannel: unsupported data type "
            + data.getClass().getName());
        global.thalion.ttio.providers.StorageDataset ds;
        try {
            ds = sc.createDataset(name, precision, len, 65536, codec, 6);
        } catch (UnsupportedOperationException e) {
            ds = sc.createDataset(name, precision, len, 0,
                Enums.Compression.NONE, 0);
        }
        try (var closeMe = ds) { closeMe.writeAll(data); }
    }

    static void writeCompoundOneCol(
            global.thalion.ttio.providers.StorageGroup sc,
            String name,
            List<global.thalion.ttio.providers.CompoundField> fields,
            List<String> values) {
        List<Object[]> rows = new ArrayList<>(values.size());
        for (String v : values) rows.add(new Object[]{ v });
        try (var ds = sc.createCompoundDataset(name, fields, rows.size())) {
            ds.writeAll(rows);
        }
    }

    /** same as {@link #writeCompoundOneCol} but encodes
     *  values as UTF-8 byte[] for compound VL_BYTES fields (the
     *  Java-side workaround for the JHI5 VL_STRING-in-compound limit). */
    static void writeCompoundOneColBytes(
            global.thalion.ttio.providers.StorageGroup sc,
            String name,
            List<global.thalion.ttio.providers.CompoundField> fields,
            List<String> values) {
        List<Object[]> rows = new ArrayList<>(values.size());
        for (String v : values) {
            rows.add(new Object[]{
                v.getBytes(java.nio.charset.StandardCharsets.UTF_8) });
        }
        try (var ds = sc.createCompoundDataset(name, fields, rows.size())) {
            ds.writeAll(rows);
        }
    }

}
