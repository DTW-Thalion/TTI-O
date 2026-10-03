/*
 * TTIOSpectralDataset+GenomicWrite.m
 * TTI-O Objective-C Implementation
 *
 * Category: TTIOSpectralDataset (GenomicWrite)
 *
 * Genomic-modality write path for TTIOSpectralDataset, split out of the
 * 4388-LOC TTIOSpectralDataset.m god-file (OO-assessment P3.10). Pure
 * code-movement: the file-static genomic-write C helpers and the
 * genomic/MS flat-buffer class-method writers (+writeGenomicRun:,
 * +writeGenomicRunStorage:, +writeMSRunStorage:,
 * +writeMinimalGenomicViaProviderURL:, +writeMinimalToPath: family) were
 * relocated here verbatim. No public API / .tio wire / behaviour change.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#include <pthread.h>
#import "TTIOSpectralDataset.h"
#import "Dataset/TTIOSpectralDataset+AssemblyWrite.h"
#import "Assembly/TTIOWrittenAssemblyGraph.h"
#import "TTIOWrittenRun.h"
#import "TTIOIdentification.h"
#import "TTIOQuantification.h"
#import "TTIOProvenanceRecord.h"
#import "TTIOTransitionList.h"
#import "TTIOCompoundIO.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOSpectrumIndex.h"
#import "Spectra/TTIONMRSpectrum.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "Genomics/TTIOPackedReference.h"
#import "Genomics/TTIOLazyReference.h"
#import "HDF5/TTIOHDF5Dataset.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOHDF5Types.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "Protection/TTIOEncryptionManager.h"
#import "Protection/TTIOAccessPolicy.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOCompoundField.h"
#import "Providers/TTIOHDF5Provider.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOReferenceImport.h"
#import "Genomics/TTIOBulkV2Blobs.h"           // Phase 2c-T
#import "Codecs/TTIORans.h"
#import "Codecs/TTIOBasePack.h"
#import "Codecs/TTIOQuality.h"                 // M86 Phase D
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Genomics/TTIOGenomicWriteContext.h"
#import "Genomics/TTIOGenomicStreamWriter.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"             // M94.Z v1.2
#import "Codecs/TTIODeltaRans.h"                // M95 v1.2
#import "Codecs/TTIOFloatDeltaZstd.h"           // codec id 17, Phase 2 MS default
#import "Codecs/TTIOMateInfoV2.h"               // inline mate-pair codec
#import "Codecs/TTIORefDiffV2.h"               // bit-packed ref-diff v2
#import "Codecs/TTIONameTokenizerV2.h"          // v1.8 #11 ch3: adaptive name-tokenizer v2
#import "Codecs/Registry/TTIOCodecRegistry.h"   // Task 6: codec-registry routing
#import "Codecs/Registry/TTIOCodec.h"
#import "Codecs/Registry/TTIODecodedChannel.h"
#import "Codecs/Registry/TTIOEncodedChannel.h"
#import "Codecs/Registry/TTIOCodecContext.h"
#import <hdf5.h>
#include <objc/message.h>                          // Task 6: typed objc_msgSend bridge
#include <openssl/md5.h>                          // M93 v1.2 ref MD5
#import "TTIOSpectralDataset+Internal.h"
@class TTIOWrittenGenomicRun;
static BOOL _TTIO_M94_RunIsV15Candidate(TTIOWrittenGenomicRun *run);
static BOOL _TTIO_V18_UseRefDiffV2(TTIOWrittenGenomicRun *run);
// ── signal-channel codec wiring ────────────────────────────────
//
// Validation, codec dispatch (rANS / BASE_PACK), and the uint8
// @compression attribute write that the read path keys on. See
// HANDOFF.md M86 §2 + Binding Decisions §86–§89.

static NSSet *_TTIO_M86_AllowedOverrideChannels_storage = nil;
static void _TTIO_M86_AllowedOverrideChannels_init(void)
{
    // read_names joins sequences/qualities as override-eligible (only
    // NAME_TOKENIZED, §113). cigars accepts {RANS0/1, NAME_TOKENIZED}
    // (§120). mate_info_{chrom,pos,tlen} are the three per-field
    // virtual channels that trigger the mate_info schema lift
    // (§125–126); bare "mate_info" remains rejected (§143). The
    // integer channels (positions/flags/mapping_qualities) are now
    // stored only under genomic_index/ — _TTIO_M86_DroppedIntChannels
    // produces the dedicated v1.6 rejection.
    _TTIO_M86_AllowedOverrideChannels_storage = [NSSet setWithArray:@[
        @"sequences", @"qualities", @"read_names", @"cigars",
        @"mate_info_chrom", @"mate_info_pos", @"mate_info_tlen",
    ]];
}
static NSSet *_TTIO_M86_AllowedOverrideChannels(void)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, _TTIO_M86_AllowedOverrideChannels_init);
    return _TTIO_M86_AllowedOverrideChannels_storage;
}

/** Per-channel allowed-codec map (M86 Phase D §119, Phase E §113).
 *  Sequences accepts the three byte-stream codecs from Phase A
 *  (rANS-0/1, BASE_PACK). Qualities additionally accepts
 *  QUALITY_BINNED (M85 Phase A codec id 7), wired here in Phase D.
 *  read_names accepts only NAME_TOKENIZED (id 8) — the codec
 *  tokenises UTF-8 strings (digit-runs vs string-runs) and is not
 *  meaningful on the byte-stream channels (Binding Decision §113);
 *  conversely the byte-stream codecs are not valid on read_names
 *  because the source data is NSArray<NSString *>, not NSData. */
static NSDictionary<NSString *, NSSet<NSNumber *> *>
    *_TTIO_M86_AllowedOverrideCodecsByChannel_storage = nil;
static void _TTIO_M86_AllowedOverrideCodecsByChannel_init(void)
{
    // REF_DIFF v1 (id 9) removed from sequences override surface — use
    // the default (refdiff_v2) or RANS / BASE_PACK. NAME_TOKENIZED v1
    // (id 8) similarly removed from read_names; default is name_tok_v2.
    NSSet *seqAllowed = [NSSet setWithArray:@[
        @(TTIOCompressionRansOrder0),
        @(TTIOCompressionRansOrder1),
        @(TTIOCompressionBasePack),
        // M103: context-aware (the read lengths), written by
        // _TTIO_M103_WriteSequencesSeqCm*.
        @(TTIOCompressionSeqCm),
    ]];
    NSSet *qualAllowed = [NSSet setWithArray:@[
        @(TTIOCompressionRansOrder0),
        @(TTIOCompressionRansOrder1),
        @(TTIOCompressionBasePack),
        @(TTIOCompressionQualityBinned),
        @(TTIOCompressionFqzcompNx16Z),  // M94.Z v1.2
    ]];
    NSSet *nameAllowed = [NSSet set];
    // cigars: rANS pair only over length-prefix-concat CIGARs (§2.5,
    // Gotcha §139). BASE_PACK / QUALITY_BINNED are wrong-content
    // (CIGAR digits + MIDNSHP=X — neither ACGT nor Phred).
    NSSet *cigarAllowed = [NSSet setWithArray:@[
        @(TTIOCompressionRansOrder0),
        @(TTIOCompressionRansOrder1),
    ]];
    // M86 Phase B (§117): integer channels accept rANS only. The
    // others (BASE_PACK, QUALITY_BINNED, NAME_TOKENIZED) don't preserve
    // int64/uint32/uint8 values. rANS is content-agnostic over the LE
    // byte representation (§118).
    NSSet *intAllowed = [NSSet setWithArray:@[
        @(TTIOCompressionRansOrder0),
        @(TTIOCompressionRansOrder1),
        @(TTIOCompressionDeltaRansOrder0),  // delta + rANS
    ]];
    // mate_info_chrom shares cigars' allowed set (rANS pair over
    // length-prefix-concat). The two integer fields mirror the
    // existing integer channels (rANS pair only).
    NSSet *mateChromAllowed = [NSSet setWithArray:@[
        @(TTIOCompressionRansOrder0),
        @(TTIOCompressionRansOrder1),
    ]];
    _TTIO_M86_AllowedOverrideCodecsByChannel_storage = @{
        @"sequences":         seqAllowed,
        @"qualities":         qualAllowed,
        @"read_names":        nameAllowed,
        @"cigars":            cigarAllowed,
        @"mate_info_chrom":   mateChromAllowed,
        @"mate_info_pos":     intAllowed,
        @"mate_info_tlen":    intAllowed,
    };
}
static NSDictionary<NSString *, NSSet<NSNumber *> *> *_TTIO_M86_AllowedOverrideCodecsByChannel(void)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, _TTIO_M86_AllowedOverrideCodecsByChannel_init);
    return _TTIO_M86_AllowedOverrideCodecsByChannel_storage;
}

/** v1.6 (L4): per-record integer metadata channels no longer accept
 *  signal_codec_overrides — they are stored exclusively under
 *  genomic_index/ now. Validation raises a dedicated v1.6 error
 *  pointing at genomic_index/ when one of these keys is present. */
static NSSet<NSString *> *_TTIO_M86_DroppedIntChannels_storage = nil;
static void _TTIO_M86_DroppedIntChannels_init(void)
{
    _TTIO_M86_DroppedIntChannels_storage = [NSSet setWithArray:@[
        @"positions", @"flags", @"mapping_qualities",
    ]];
}
static NSSet<NSString *> *_TTIO_M86_DroppedIntChannels(void)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, _TTIO_M86_DroppedIntChannels_init);
    return _TTIO_M86_DroppedIntChannels_storage;
}

/** Validate the per-channel codec overrides BEFORE any HDF5 mutation.
 *  Raises NSInvalidArgumentException on programmer error so the file
 *  is left untouched (Binding Decision §88, HANDOFF.md M86 §3). */
static void _TTIO_M86_ValidateOverrides(NSDictionary<NSString *, NSNumber *> *overrides)
{
    if (overrides.count == 0) return;
    NSSet *allowedChans = _TTIO_M86_AllowedOverrideChannels();
    NSDictionary<NSString *, NSSet<NSNumber *> *> *allowedByChan =
        _TTIO_M86_AllowedOverrideCodecsByChannel();
    NSSet<NSString *> *droppedIntChannels = _TTIO_M86_DroppedIntChannels();
    for (NSString *chName in overrides) {
        // v1.6 (L4): per-record integer metadata channels removed from
        // the signal_channels/ override surface. They live exclusively
        // under genomic_index/ now (mirroring MS's spectrum_index/
        // pattern). Hard-error so callers with stale code learn
        // immediately. See docs/format-spec.md §4 and §10.7.
        if ([droppedIntChannels containsObject:chName]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides[\"%@\"]: removed in "
                               @"v1.6 — per-record integer metadata fields "
                               @"(positions, flags, mapping_qualities) are "
                               @"stored only under genomic_index/, not "
                               @"signal_channels/. The override no longer "
                               @"applies. See docs/format-spec.md §4 and "
                               @"§10.7.", chName];
        }
        // M86 Phase F (Binding Decision §126, Gotcha §143): the bare
        // 'mate_info' key is reserved and rejected with a discoverable
        // error pointing at the three per-field virtual channel names.
        // Producing this dedicated message before the generic
        // unknown-channel rejection makes the migration path obvious.
        if ([chName isEqualToString:@"mate_info"]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides['mate_info']: the "
                               @"bare 'mate_info' key is reserved and "
                               @"rejected — mate_info is decomposed at "
                               @"the per-field level in M86 Phase F. Use "
                               @"one or more of the three per-field keys "
                               @"instead: 'mate_info_chrom', "
                               @"'mate_info_pos', 'mate_info_tlen'. See "
                               @"docs/format-spec.md §10.9."];
        }
        if (![allowedChans containsObject:chName]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides: channel '%@' not "
                               @"supported (only sequences, qualities, "
                               @"read_names, cigars, positions, flags, "
                               @"mapping_qualities, mate_info_chrom, "
                               @"mate_info_pos, and mate_info_tlen can "
                               @"use TTIO codecs)",
                               chName];
        }
        NSNumber *codecBox = overrides[chName];
        if (![codecBox isKindOfClass:[NSNumber class]]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides['%@']: codec value "
                               @"must be an NSNumber-boxed TTIOCompression",
                               chName];
        }
        NSSet<NSNumber *> *allowed = allowedByChan[chName];
        if (![allowed containsObject:codecBox]) {
            // Phase D Binding Decision §110: explicit message for the
            // (sequences, QUALITY_BINNED) category error — names the
            // codec, the channel, and the lossy-quantisation rationale.
            TTIOCompression codec =
                (TTIOCompression)[codecBox unsignedIntegerValue];
            if (codec == TTIOCompressionQualityBinned
                && [chName isEqualToString:@"sequences"]) {
                [NSException raise:NSInvalidArgumentException
                            format:@"signalCodecOverrides['%@']: codec "
                                   @"QUALITY_BINNED is not valid on the "
                                   @"'%@' channel — quality binning is "
                                   @"lossy and only applies to Phred "
                                   @"quality scores. Applying it to ACGT "
                                   @"sequence bytes would silently destroy "
                                   @"the sequence via Phred-bin "
                                   @"quantisation. Use the 'qualities' "
                                   @"channel for QUALITY_BINNED, or "
                                   @"RansOrder0/RansOrder1/BasePack on "
                                   @"sequences.", chName, chName];
            }
            if ([chName isEqualToString:@"read_names"]) {
                // v1.0 reset: read_names is v2-only. The writer always
                // emits NAME_TOKENIZED_V2 (codec id 15) when the native
                // libttio_rans is linked; no per-call override is
                // accepted on this channel.
                [NSException raise:NSInvalidArgumentException
                            format:@"signalCodecOverrides['%@']: codec "
                                   @"%@ is not supported on the '%@' "
                                   @"channel — read_names is v2-only "
                                   @"in v1.0+ (NAME_TOKENIZED_V2, "
                                   @"codec id 15). Drop the override.",
                                   chName, codecBox, chName];
            }
            // M86 Phase C Binding Decision §120: explicit messages for
            // wrong-content codecs on the cigars channel. CIGAR strings
            // contain ASCII digits + operator letters (MIDNSHP=X), none
            // of which are ACGT bases or Phred quality values.
            if ([chName isEqualToString:@"cigars"]) {
                if (codec == TTIOCompressionBasePack) {
                    [NSException raise:NSInvalidArgumentException
                                format:@"signalCodecOverrides['%@']: codec "
                                       @"BASE_PACK is not valid on the "
                                       @"'cigars' channel — BASE_PACK "
                                       @"2-bit-packs ACGT sequence bytes "
                                       @"and would silently corrupt the "
                                       @"CIGAR strings stored on this "
                                       @"channel (CIGAR ASCII contains "
                                       @"digits and operator letters "
                                       @"MIDNSHP=X, none of which are "
                                       @"ACGT). Use RANS_ORDER0 or "
                                       @"RANS_ORDER1 on 'cigars'.",
                                       chName];
                }
                if (codec == TTIOCompressionQualityBinned) {
                    [NSException raise:NSInvalidArgumentException
                                format:@"signalCodecOverrides['%@']: codec "
                                       @"QUALITY_BINNED is not valid on "
                                       @"the 'cigars' channel — "
                                       @"QUALITY_BINNED quantises Phred "
                                       @"quality scores onto an 8-bin "
                                       @"centre table and would silently "
                                       @"destroy the CIGAR strings stored "
                                       @"on this channel. Use RANS_ORDER0 "
                                       @"or RANS_ORDER1 on 'cigars'.",
                                       chName];
                }
                [NSException raise:NSInvalidArgumentException
                            format:@"signalCodecOverrides['%@']: codec %@ "
                                   @"not supported on the '%@' channel "
                                   @"(allowed: RansOrder0, RansOrder1)",
                                   chName, codecBox, chName];
            }
            // mate_info_chrom shares cigars'
            // allowed set ({RANS_ORDER0, RANS_ORDER1}) — NAME_TOKENIZED
            // v1 (id 8) was dropped. Wrong-content rejection mirrors
            // cigars' messaging — chromosome names are short ASCII
            // strings (typically <30 distinct values), none of them
            // ACGT or Phred values.
            if ([chName isEqualToString:@"mate_info_chrom"]) {
                if (codec == TTIOCompressionBasePack) {
                    [NSException raise:NSInvalidArgumentException
                                format:@"signalCodecOverrides['%@']: codec "
                                       @"BASE_PACK is not valid on the "
                                       @"'mate_info_chrom' channel — "
                                       @"BASE_PACK 2-bit-packs ACGT "
                                       @"sequence bytes and would "
                                       @"silently corrupt the chromosome "
                                       @"names stored on this channel. "
                                       @"Use RANS_ORDER0 or RANS_ORDER1 "
                                       @"on 'mate_info_chrom'.",
                                       chName];
                }
                if (codec == TTIOCompressionQualityBinned) {
                    [NSException raise:NSInvalidArgumentException
                                format:@"signalCodecOverrides['%@']: codec "
                                       @"QUALITY_BINNED is not valid on "
                                       @"the 'mate_info_chrom' channel — "
                                       @"QUALITY_BINNED quantises Phred "
                                       @"quality scores and would "
                                       @"silently destroy the chromosome "
                                       @"names stored on this channel. "
                                       @"Use RANS_ORDER0 or RANS_ORDER1 "
                                       @"on 'mate_info_chrom'.",
                                       chName];
                }
                [NSException raise:NSInvalidArgumentException
                            format:@"signalCodecOverrides['%@']: codec %@ "
                                   @"not supported on the '%@' channel "
                                   @"(allowed: RansOrder0, RansOrder1)",
                                   chName, codecBox, chName];
            }
            NSString *allowedTail =
                [chName isEqualToString:@"qualities"] ? @", QualityBinned" : @"";
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides['%@']: codec %@ "
                               @"not supported on the '%@' channel "
                               @"(allowed: RansOrder0, RansOrder1, "
                               @"BasePack%@)",
                               chName, codecBox, chName, allowedTail];
        }
    }
}

/** Encode raw bytes through the selected M86 codec. */
static NSData *_TTIO_M86_EncodeWithCodec(NSData *raw, TTIOCompression codec)
{
    id<TTIOCodec> c = [TTIOCodecRegistry codecForId:codec];
    if (c == nil) {
        [NSException raise:NSInvalidArgumentException
                    format:@"_TTIO_M86_EncodeWithCodec: codec %lu not "
                           @"a TTIO byte-stream codec",
                           (unsigned long)codec];
        return nil;
    }
    NSError *e = nil;
    TTIOEncodedChannel *enc =
        [c encode:[[TTIODecodedBytes alloc] initWithData:raw]
           context:[TTIOCodecContext emptyContext]
             error:&e];
    return ((TTIOEncodedDatasetBytes *)enc).bytes;
}

/** Encode read-names through the NAME_TOKENIZED_V2 codec via the
 *  registry. Byte-identical to [TTIONameTokenizerV2 encodeNames:]. */
static NSData *_TTIO_M86_EncodeNamesViaRegistry(NSArray<NSString *> *names)
{
    id<TTIOCodec> c =
        [TTIOCodecRegistry codecForId:TTIOCompressionNameTokenizedV2];
    if (c == nil) return nil;
    NSError *e = nil;
    TTIOEncodedChannel *enc =
        [c encode:[[TTIODecodedStringList alloc] initWithNames:names]
           context:[TTIOCodecContext emptyContext]
             error:&e];
    return ((TTIOEncodedDatasetBytes *)enc).bytes;
}

/** Encode mate-info through the MATE_INLINE_V2 codec via the registry.
 *  Byte-identical to [TTIOMateInfoV2 encodeMateChromIds:...]. */
static NSData *_TTIO_M86_EncodeMateInfoViaRegistry(NSData *mateChromIds,
                                                   NSData *matePositions,
                                                   NSData *templateLengths,
                                                   NSData *ownChromIds,
                                                   NSData *ownPositions,
                                                   NSError **error)
{
    id<TTIOCodec> c =
        [TTIOCodecRegistry codecForId:TTIOCompressionMateInlineV2];
    if (c == nil) return nil;
    TTIOCodecContext *ctx = [TTIOCodecContext emptyContext];
    ctx.ownChromIds = ownChromIds;
    ctx.ownPositions = ownPositions;
    TTIODecodedMateInfo *mi =
        [[TTIODecodedMateInfo alloc] initWithMateChromIds:mateChromIds
                                            matePositions:matePositions
                                          templateLengths:templateLengths];
    TTIOEncodedChannel *enc = [c encode:mi context:ctx error:error];
    return ((TTIOEncodedDatasetBytes *)enc).bytes;
}

// unsigned LEB128 varint writer for the cigars rANS path.
// The serialisation contract is `varint(asciiLen) + asciiBytes` per
// CIGAR (§2.5 of the Phase C plan; mirrors NAME_TOKENIZED's verbatim
// format minus the 7-byte header). Same wire format as the codec's
// own internal varint helpers (see TTIONameTokenizer.m); reproduced
// here to avoid coupling the dataset writer to the codec module's
// private symbols.
static void _TTIO_M86_VarintWrite(NSMutableData *out, uint64_t value)
{
    uint8_t buf[10];
    size_t n = 0;
    while (value >= 0x80u) {
        buf[n++] = (uint8_t)((value & 0x7Fu) | 0x80u);
        value >>= 7;
    }
    buf[n++] = (uint8_t)(value & 0x7Fu);
    [out appendBytes:buf length:n];
}

/** encode a list of CIGAR strings via the selected codec.
 *
 *  only the rANS pair is accepted —
 *  NAME_TOKENIZED v1 (codec id 8) was removed.
 *
 *  RANS_ORDER0 / RANS_ORDER1: serialise the list as length-prefix-
 *  concat (varint(asciiLen) + asciiBytes per CIGAR — §2.5, Gotcha
 *  §139), then pass the concatenated buffer through TTIORansEncode.
 *
 *  Raises NSInvalidArgumentException if any CIGAR contains non-ASCII
 *  bytes (SAM spec is 7-bit ASCII). */
static NSData *_TTIO_M86_EncodeCigarsWithCodec(NSArray<NSString *> *cigars,
                                                TTIOCompression codec)
{
    if (codec == TTIOCompressionRansOrder0
        || codec == TTIOCompressionRansOrder1) {
        NSMutableData *buf = [NSMutableData data];
        for (NSUInteger idx = 0; idx < cigars.count; idx++) {
            NSString *cig = cigars[idx];
            const char *ascii = [cig cStringUsingEncoding:NSASCIIStringEncoding];
            if (ascii == NULL) {
                [NSException raise:NSInvalidArgumentException
                            format:@"signalCodecOverrides['cigars']: cigar "
                                   @"at index %lu contains non-ASCII bytes "
                                   @"— CIGARs must be 7-bit ASCII per the "
                                   @"SAM spec",
                                   (unsigned long)idx];
            }
            NSUInteger nBytes = strlen(ascii);
            _TTIO_M86_VarintWrite(buf, (uint64_t)nBytes);
            [buf appendBytes:ascii length:nBytes];
        }
        int order = (codec == TTIOCompressionRansOrder0) ? 0 : 1;
        return TTIORansEncode(buf, order);
    }
    [NSException raise:NSInvalidArgumentException
                format:@"_TTIO_M86_EncodeCigarsWithCodec: codec %lu not a "
                       @"valid cigars codec (only RANS_ORDER0 or "
                       @"RANS_ORDER1)",
                       (unsigned long)codec];
    return nil;
}

/** Set @compression as a uint8 attribute on an HDF5 dataset. Matches
 *  Python's ``write_int_attr(ds, "compression", n, dtype="<u1")``
 *  byte-for-byte (Binding Decision §86, HANDOFF.md M86 §5.1). */
static BOOL _TTIO_M86_WriteUInt8Attribute(hid_t did, const char *name,
                                          uint8_t value, NSError **error)
{
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2001
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"H5Screate(SCALAR) failed for @compression"}];
        return NO;
    }
    if (H5Aexists(did, name) > 0) {
        H5Adelete(did, name);
    }
    hid_t aid = H5Acreate2(did, name, H5T_NATIVE_UINT8, space,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (aid < 0) {
        H5Sclose(space);
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2002
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"H5Acreate2(@%s) failed", name]}];
        return NO;
    }
    herr_t s = H5Awrite(aid, H5T_NATIVE_UINT8, &value);
    H5Aclose(aid); H5Sclose(space);
    if (s < 0) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2003
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"H5Awrite(@%s) failed", name]}];
        return NO;
    }
    return YES;
}

/** Read @compression as a uint8 from an HDF5 dataset. Returns 0 when
 *  the attribute is absent (matches Python's read_int_attr default).
 *  ``*outExists`` (if non-NULL) signals whether the attribute was
 *  found, so callers can distinguish "absent" from "explicitly 0". */
static uint8_t _TTIO_M86_ReadUInt8Attribute(hid_t did, const char *name,
                                            BOOL *outExists)
{
    if (H5Aexists(did, name) <= 0) {
        if (outExists) *outExists = NO;
        return 0;
    }
    if (outExists) *outExists = YES;
    hid_t aid = H5Aopen(did, name, H5P_DEFAULT);
    if (aid < 0) return 0;
    uint8_t value = 0;
    H5Aread(aid, H5T_NATIVE_UINT8, &value);
    H5Aclose(aid);
    return value;
}

/** Write a uint8 byte channel either through the existing HDF5 filter
 *  (when no override) or through a TTIO codec (when overridden). For
 *  the codec path we skip the HDF5 filter entirely (Binding Decision
 *  §87 — no double-compression). The @compression attribute is set on
 *  the dataset for the read-side dispatcher. */
static BOOL _TTIO_M86_WriteByteChannel(TTIOHDF5Group *group,
                                       NSString *name,
                                       NSData *data,
                                       TTIOCompression defaultCompression,
                                       NSNumber *codecOverride,
                                       NSError **error)
{
    if (codecOverride == nil) {
        // Plain path — same behaviour as the M82 byte-channel write.
        TTIOHDF5Dataset *ds = [group createDatasetNamed:name
                                              precision:TTIOPrecisionUInt8
                                                 length:data.length
                                              chunkSize:65536
                                            compression:defaultCompression
                                       compressionLevel:6
                                                  error:error];
        if (!ds) return NO;
        return [ds writeData:data error:error];
    }

    TTIOCompression codec = (TTIOCompression)[codecOverride unsignedIntegerValue];
    NSData *encoded = _TTIO_M86_EncodeWithCodec(data, codec);
    if (!encoded) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2010
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 codec %lu encode failed for channel '%@'",
                            (unsigned long)codec, name]}];
        return NO;
    }
    // Codec-compressed datasets carry NO HDF5 filter.
    TTIOHDF5Dataset *ds = [group createDatasetNamed:name
                                          precision:TTIOPrecisionUInt8
                                             length:encoded.length
                                          chunkSize:65536
                                        compression:TTIOCompressionNone
                                   compressionLevel:0
                                              error:error];
    if (!ds) return NO;
    if (![ds writeData:encoded error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)codec, error);
}

/** Provider-path twin of _TTIO_M86_WriteByteChannel for non-HDF5
 *  backends (memory://, sqlite://). The @compression attribute uses
 *  the storage protocol's setAttributeValue:forName: which boxes as
 *  NSNumber → int64 in the HDF5 backend; non-HDF5 backends simply
 *  store the integer. The cross-language fixture matrix only covers
 *  HDF5, so the protocol path doesn't need byte-exact parity. */
static BOOL _TTIO_M86_WriteByteChannelStorage(id<TTIOStorageGroup> group,
                                              NSString *name,
                                              NSData *data,
                                              TTIOCompression defaultCompression,
                                              NSNumber *codecOverride,
                                              NSError **error)
{
    if (codecOverride == nil) {
        id<TTIOStorageDataset> ds = [group createDatasetNamed:name
                                                    precision:TTIOPrecisionUInt8
                                                       length:data.length
                                                    chunkSize:65536
                                                  compression:defaultCompression
                                             compressionLevel:6
                                                        error:error];
        if (!ds) return NO;
        return [ds writeAll:data error:error];
    }

    TTIOCompression codec = (TTIOCompression)[codecOverride unsignedIntegerValue];
    NSData *encoded = _TTIO_M86_EncodeWithCodec(data, codec);
    if (!encoded) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2011
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 codec %lu encode failed for channel '%@'",
                            (unsigned long)codec, name]}];
        return NO;
    }
    id<TTIOStorageDataset> ds = [group createDatasetNamed:name
                                                precision:TTIOPrecisionUInt8
                                                   length:encoded.length
                                                chunkSize:65536
                                              compression:TTIOCompressionNone
                                         compressionLevel:0
                                                    error:error];
    if (!ds) return NO;
    if (![ds writeAll:encoded error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)codec)
                          forName:@"compression"
                            error:error];
}

// ── integer-channel codec wiring ───────────────────────
//
// Per-channel integer dtypes for the int↔byte serialisation contract
// (Binding Decision §115). Determined by **channel name lookup**; the
// reader uses the same map to interpret the decoded byte buffer back
// to the channel's natural integer dtype, so no extra on-disk
// attribute is required beyond ``@compression``.
static TTIOPrecision _TTIO_M86_IntegerChannelPrecision(NSString *name)
{
    if ([name isEqualToString:@"positions"])         return TTIOPrecisionInt64;
    if ([name isEqualToString:@"flags"])             return TTIOPrecisionUInt32;
    if ([name isEqualToString:@"mapping_qualities"]) return TTIOPrecisionUInt8;
    return (TTIOPrecision)0;  // unreachable; validation rejects others
}

/** Serialise an integer signal-channel buffer to little-endian bytes
 *  for the rANS codec. The input is the in-memory NSData buffer the
 *  WrittenGenomicRun carries (host endianness). The output is the LE
 *  byte representation per Binding Decision §118 — non-negotiable so
 *  big-endian platforms produce identical wire bytes. We byte-swap
 *  per element on big-endian hosts; on x86/ARM this is a memcpy
 *  no-op (Gotcha §131 — uint8 is always trivially a no-op). */
static NSData *_TTIO_M86_IntChannelToLEBytes(NSString *name, NSData *data)
{
    if ([name isEqualToString:@"positions"]) {
        const int64_t *src = (const int64_t *)data.bytes;
        NSUInteger n = data.length / sizeof(int64_t);
        NSMutableData *out = [NSMutableData dataWithLength:n * sizeof(int64_t)];
        int64_t *dst = (int64_t *)out.mutableBytes;
        for (NSUInteger i = 0; i < n; i++) {
            uint64_t le = TTIO_HOST_TO_LE64((uint64_t)src[i]);
            memcpy(&dst[i], &le, sizeof(uint64_t));
        }
        return out;
    }
    if ([name isEqualToString:@"flags"]) {
        const uint32_t *src = (const uint32_t *)data.bytes;
        NSUInteger n = data.length / sizeof(uint32_t);
        NSMutableData *out = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
        uint32_t *dst = (uint32_t *)out.mutableBytes;
        for (NSUInteger i = 0; i < n; i++) {
            dst[i] = TTIO_HOST_TO_LE32(src[i]);
        }
        return out;
    }
    // mapping_qualities (uint8) — LE no-op (Gotcha §131).
    return [data copy];
}

/** element size in bytes for a named integer channel.  Used by
 *  the delta-rANS encoder which needs the element width to compute
 *  deltas across typed values rather than raw bytes. */
static uint8_t _TTIO_M95_IntChannelElementSize(NSString *name) {
    if ([name isEqualToString:@"positions"])         return 8;
    if ([name isEqualToString:@"flags"])             return 4;
    if ([name isEqualToString:@"mapping_qualities"]) return 1;
    if ([name isEqualToString:@"mate_info_pos"])     return 8;
    if ([name isEqualToString:@"mate_info_tlen"])    return 4;
    return 0;  // unreachable; validation rejects others
}

/** Write an integer signal channel either directly with the M82 typed
 *  dataset (when no override) or through the rANS/delta-rANS codec
 *  with the LE-serialisation contract (when overridden). HDF5 fast
 *  path. */
static BOOL _TTIO_M86_WriteIntChannel(TTIOHDF5Group *group,
                                      NSString *name,
                                      NSData *data,
                                      TTIOCompression defaultCompression,
                                      NSNumber *codecOverride,
                                      NSError **error)
{
    TTIOPrecision prec = _TTIO_M86_IntegerChannelPrecision(name);
    if (codecOverride == nil) {
        // M82 typed path — preserves byte parity with pre-Phase-B files.
        NSUInteger n = data.length / TTIOPrecisionElementSize(prec);
        TTIOHDF5Dataset *ds = [group createDatasetNamed:name
                                              precision:prec
                                                 length:n
                                              chunkSize:65536
                                            compression:defaultCompression
                                       compressionLevel:6
                                                  error:error];
        if (!ds) return NO;
        return [ds writeData:data error:error];
    }

    TTIOCompression codec = (TTIOCompression)[codecOverride unsignedIntegerValue];
    if (codec != TTIOCompressionRansOrder0
        && codec != TTIOCompressionRansOrder1
        && codec != TTIOCompressionDeltaRansOrder0) {
        // Defensive — _TTIO_M86_ValidateOverrides rejects this first.
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2050
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 Phase B: codec %lu is not valid on "
                            @"integer channel '%@' (only RANS_ORDER0/"
                            @"RANS_ORDER1/DELTA_RANS_ORDER0 supported)",
                            (unsigned long)codec, name]}];
        return NO;
    }
    NSData *leBytes = _TTIO_M86_IntChannelToLEBytes(name, data);
    NSData *encoded = nil;
    if (codec == TTIOCompressionDeltaRansOrder0) {
        uint8_t elemSize = _TTIO_M95_IntChannelElementSize(name);
        NSError *encErr = nil;
        encoded = TTIODeltaRansEncode(leBytes, elemSize, &encErr);
        if (!encoded) {
            if (error) *error = encErr;
            return NO;
        }
    } else {
        int order = (codec == TTIOCompressionRansOrder0) ? 0 : 1;
        encoded = TTIORansEncode(leBytes, order);
    }
    if (!encoded) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2051
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 Phase B: codec encode failed for "
                            @"integer channel '%@'", name]}];
        return NO;
    }
    // Codec-compressed datasets carry NO HDF5 filter (Binding Decision §87).
    TTIOHDF5Dataset *ds = [group createDatasetNamed:name
                                          precision:TTIOPrecisionUInt8
                                             length:encoded.length
                                          chunkSize:65536
                                        compression:TTIOCompressionNone
                                   compressionLevel:0
                                              error:error];
    if (!ds) return NO;
    if (![ds writeData:encoded error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)codec, error);
}

/** Provider-path twin of _TTIO_M86_WriteIntChannel for non-HDF5
 *  backends (memory://, sqlite://). */
static BOOL _TTIO_M86_WriteIntChannelStorage(id<TTIOStorageGroup> group,
                                             NSString *name,
                                             NSData *data,
                                             TTIOCompression defaultCompression,
                                             NSNumber *codecOverride,
                                             NSError **error)
{
    TTIOPrecision prec = _TTIO_M86_IntegerChannelPrecision(name);
    if (codecOverride == nil) {
        NSUInteger n = data.length / TTIOPrecisionElementSize(prec);
        id<TTIOStorageDataset> ds = [group createDatasetNamed:name
                                                    precision:prec
                                                       length:n
                                                    chunkSize:65536
                                                  compression:defaultCompression
                                             compressionLevel:6
                                                        error:error];
        if (!ds) return NO;
        return [ds writeAll:data error:error];
    }

    TTIOCompression codec = (TTIOCompression)[codecOverride unsignedIntegerValue];
    if (codec != TTIOCompressionRansOrder0
        && codec != TTIOCompressionRansOrder1
        && codec != TTIOCompressionDeltaRansOrder0) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2052
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 Phase B: codec %lu is not valid on "
                            @"integer channel '%@' (only RANS_ORDER0/"
                            @"RANS_ORDER1/DELTA_RANS_ORDER0 supported)",
                            (unsigned long)codec, name]}];
        return NO;
    }
    NSData *leBytes = _TTIO_M86_IntChannelToLEBytes(name, data);
    NSData *encoded = nil;
    if (codec == TTIOCompressionDeltaRansOrder0) {
        uint8_t elemSize = _TTIO_M95_IntChannelElementSize(name);
        NSError *encErr = nil;
        encoded = TTIODeltaRansEncode(leBytes, elemSize, &encErr);
        if (!encoded) {
            if (error) *error = encErr;
            return NO;
        }
    } else {
        int order = (codec == TTIOCompressionRansOrder0) ? 0 : 1;
        encoded = TTIORansEncode(leBytes, order);
    }
    if (!encoded) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2053
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"M86 Phase B: codec encode failed for "
                            @"integer channel '%@'", name]}];
        return NO;
    }
    id<TTIOStorageDataset> ds = [group createDatasetNamed:name
                                                precision:TTIOPrecisionUInt8
                                                   length:encoded.length
                                                chunkSize:65536
                                              compression:TTIOCompressionNone
                                         compressionLevel:0
                                                    error:error];
    if (!ds) return NO;
    if (![ds writeAll:encoded error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)codec)
                          forName:@"compression"
                            error:error];
}


/** returns YES when the inline_v2 path should be used.
 *  Requires native libttio_rans AND a non-empty run. (v1.0 reset:
 *  opt-out flag removed; empty runs still take the M82 compound
 *  fallback because the inline_v2 encoder requires n > 0.) */
static BOOL _TTIO_V17_UseMateInlineV2(TTIOWrittenGenomicRun *run)
{
    if (run.readCount == 0) return NO;
    return [TTIOMateInfoV2 nativeAvailable];
}

/** reject mate_info_* per-field overrides
 *  unconditionally — the per-field subgroup writer was deleted, so
 *  inline_v2 is the only mate_info layout under v1.0. Called after
 *  the standard _TTIO_M86_ValidateOverrides check so baseline
 *  unknown-channel errors are already handled. */
static NSSet<NSString *> *_TTIO_V17_MateKeys_storage = nil;
static void _TTIO_V17_MateKeys_init(void)
{
    _TTIO_V17_MateKeys_storage = [NSSet setWithArray:@[
        @"mate_info_chrom", @"mate_info_pos", @"mate_info_tlen",
    ]];
}
static void _TTIO_V17_ValidateMateInfoV2Overrides(TTIOWrittenGenomicRun *run)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, _TTIO_V17_MateKeys_init);
    for (NSString *chName in run.signalCodecOverrides) {
        if ([_TTIO_V17_MateKeys_storage containsObject:chName]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides['%@']: per-field "
                               @"mate_info_* overrides are no longer "
                               @"accepted under the v1.0 reset — the "
                               @"writer always emits the inline_v2 "
                               @"codec when the native libttio_rans "
                               @"library is linked.",
                               chName];
        }
    }
}

/** M97: QUALITY_BINNED is a fixed Illumina-8 bin table; reject it as
 *  a qualities override on long-read platforms rather than write a
 *  mistuned lossy file. */
static void _TTIO_M97_ValidateQualityBinnedPlatform(TTIOWrittenGenomicRun *run)
{
    NSNumber *qov = run.signalCodecOverrides[@"qualities"];
    if (qov != nil
        && qov.integerValue == TTIOCompressionQualityBinned
        && !TTIOQualityBinnedAllowedForPlatform(run.platform)) {
        [NSException raise:NSInvalidArgumentException
                    format:@"signalCodecOverrides['qualities']: "
                           @"QUALITY_BINNED applies the fixed Illumina-8 "
                           @"bin table, which does not fit the quality "
                           @"distribution of platform '%@' (M97).",
                           run.platform];
    }
}

// ── inline_v2 writer helpers ──────────────────────────────
//
// Build the chrom_id table (encounter-order on own chromosomes, extended
// with mate-only chroms), encode via TTIOMateInfoV2, and write:
//   signal_channels/mate_info/inline_v2   — uint8 dataset @compression=13
//   signal_channels/mate_info/chrom_names — compound (name:VL_STRING)
//
// Two variants: HDF5 fast path (TTIOHDF5Group *) and storage-protocol
// path (id<TTIOStorageGroup>).

/** Build the chrom_id map (encounter-order). Fills ownChromIds (uint16)
 *  and mateChromIds (int32: -1 for '*', own-chrom-id for '=', else
 *  chrom_id from encounter-ordered table extended with mate-only entries).
 *  Also fills chromNamesInOrder with the full chrom_id → name table. */
static BOOL _TTIO_V17_BuildChromTablesShared(TTIOWrittenGenomicRun *run,
                                              NSMutableDictionary<NSString *, NSNumber *> *shared,
                                              NSData * _Nonnull * _Nonnull ownChromIdsOut,
                                              NSData * _Nonnull * _Nonnull mateChromIdsOut,
                                              NSMutableArray<NSString *> * _Nonnull * _Nonnull chromNamesOut);

static BOOL _TTIO_V17_BuildChromTables(TTIOWrittenGenomicRun *run,
                                        NSData * _Nonnull * _Nonnull ownChromIdsOut,
                                        NSData * _Nonnull * _Nonnull mateChromIdsOut,
                                        NSMutableArray<NSString *> * _Nonnull * _Nonnull chromNamesOut)
{
    return _TTIO_V17_BuildChromTablesShared(run, nil, ownChromIdsOut, mateChromIdsOut, chromNamesOut);
}

// With a shared map (blocks_v1) ids are stable across blocks and the map
// grows in place; chromNamesOut then lists every name in id order.
static BOOL _TTIO_V17_BuildChromTablesShared(TTIOWrittenGenomicRun *run,
                                              NSMutableDictionary<NSString *, NSNumber *> *shared,
                                              NSData * _Nonnull * _Nonnull ownChromIdsOut,
                                              NSData * _Nonnull * _Nonnull mateChromIdsOut,
                                              NSMutableArray<NSString *> * _Nonnull * _Nonnull chromNamesOut)
{
    NSUInteger n = run.readCount;

    // Build encounter-ordered chrom_id table from own chromosomes.
    NSMutableDictionary<NSString *, NSNumber *> *nameToId =
        shared ?: [NSMutableDictionary dictionaryWithCapacity:32];
    NSMutableArray<NSString *> *chromNames = [NSMutableArray arrayWithCapacity:32];

    NSMutableData *ownChromIdsData =
        [NSMutableData dataWithLength:n * sizeof(uint16_t)];
    uint16_t *ownIds = (uint16_t *)ownChromIdsData.mutableBytes;

    NSArray<NSString *> *ownChroms = run.chromosomes;
    for (NSUInteger i = 0; i < n; i++) {
        NSString *name = ownChroms[i];
        NSNumber *existingId = nameToId[name];
        if (existingId == nil) {
            NSUInteger newId = nameToId.count;
            nameToId[name] = @(newId);
            ownIds[i] = (uint16_t)newId;
        } else {
            ownIds[i] = (uint16_t)[existingId unsignedIntegerValue];
        }
    }

    // Build mate_chrom_ids, extending nameToId/chromNames for mate-only chroms.
    NSMutableData *mateChromIdsData =
        [NSMutableData dataWithLength:n * sizeof(int32_t)];
    int32_t *mateIds = (int32_t *)mateChromIdsData.mutableBytes;

    NSArray<NSString *> *mateChroms = run.mateChromosomes;
    for (NSUInteger i = 0; i < n; i++) {
        NSString *name = mateChroms[i];
        if ([name isEqualToString:@"*"] || name == nil || name.length == 0) {
            mateIds[i] = -1;
        } else if ([name isEqualToString:@"="]) {
            // '=' means mate is on the same chrom as this read.
            mateIds[i] = (int32_t)ownIds[i];
        } else {
            NSNumber *existingId = nameToId[name];
            if (existingId == nil) {
                NSUInteger newId = nameToId.count;
                nameToId[name] = @(newId);
                mateIds[i] = (int32_t)newId;
            } else {
                mateIds[i] = (int32_t)[existingId unsignedIntegerValue];
            }
        }
    }
    [chromNames addObjectsFromArray:[TTIOGenomicIndex namesInIdOrder:nameToId]];

    *ownChromIdsOut  = ownChromIdsData;
    *mateChromIdsOut = mateChromIdsData;
    *chromNamesOut   = chromNames;
    return YES;
}

/** HDF5 fast path: write the inline_v2 group with the blob and chrom_names
 *  sidecar into signal_channels/mate_info/. */
// Phase 2c-T: write a verbatim mate_info/inline_v2 blob + chrom_names
// table, bypassing the v2 codec encode. Mirrors the layout produced by
// _TTIO_V17_WriteMateInfoInlineV2HDF5 but skips the encode step.
static BOOL _TTIO_PhaseT_WriteMateInfoBulkHDF5(TTIOHDF5Group *sc,
                                                NSData *blob,
                                                NSArray<NSString *> *chromNames,
                                                NSError **error)
{
    TTIOHDF5Group *mateGrp = [sc createGroupNamed:@"mate_info" error:error];
    if (!mateGrp) return NO;
    TTIOHDF5Dataset *ds = [mateGrp createDatasetNamed:@"inline_v2"
                                             precision:TTIOPrecisionUInt8
                                                length:blob.length
                                             chunkSize:65536
                                           compression:TTIOCompressionNone
                                      compressionLevel:0
                                                 error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeData:blob error:error]) return NO;
    if (!_TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                       (uint8_t)TTIOCompressionMateInlineV2,
                                       error)) return NO;
    NSArray *vlNameField = @[
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString]
    ];
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:chromNames.count];
    for (NSString *cn in chromNames) [rows addObject:@{@"name": cn}];
    if (![TTIOCompoundIO writeGeneric:rows
                             intoGroup:mateGrp datasetNamed:@"chrom_names"
                                fields:vlNameField error:error]) return NO;
    return YES;
}

// Phase 2c-T: write a verbatim read_names blob with @compression=15.
static BOOL _TTIO_PhaseT_WriteReadNamesBulkHDF5(TTIOHDF5Group *sc,
                                                  NSData *blob,
                                                  NSError **error)
{
    TTIOHDF5Dataset *ds = [sc createDatasetNamed:@"read_names"
                                       precision:TTIOPrecisionUInt8
                                          length:blob.length
                                       chunkSize:65536
                                     compression:TTIOCompressionNone
                                compressionLevel:0
                                           error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeData:blob error:error]) return NO;
    if (!_TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                       (uint8_t)TTIOCompressionNameTokenizedV2,
                                       error)) return NO;
    return YES;
}

// Phase 2c-T: write a verbatim sequences/refdiff_v2 blob (group layout).
static BOOL _TTIO_PhaseT_WriteRefDiffV2BulkHDF5(TTIOHDF5Group *sc,
                                                 NSData *blob,
                                                 NSError **error)
{
    TTIOHDF5Group *seqGrp = [sc createGroupNamed:@"sequences" error:error];
    if (!seqGrp) return NO;
    TTIOHDF5Dataset *ds = [seqGrp createDatasetNamed:@"refdiff_v2"
                                            precision:TTIOPrecisionUInt8
                                               length:blob.length
                                            chunkSize:65536
                                          compression:TTIOCompressionNone
                                     compressionLevel:0
                                                error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeData:blob error:error]) return NO;
    if (!_TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                       (uint8_t)TTIOCompressionRefDiffV2,
                                       error)) return NO;
    return YES;
}

// Phase 2c-T storage-protocol path mirrors of the HDF5-fast-path
// helpers above. Used by writeGenomicRunStorage so memory:// /
// sqlite:// / zarr:// receivers also honor bulk_v2_blobs.

static BOOL _TTIO_PhaseT_WriteMateInfoBulkStorage(
        id<TTIOStorageGroup> sc,
        NSData *blob, NSArray<NSString *> *chromNames,
        NSError **error)
{
    id<TTIOStorageGroup> mateGrp = [sc createGroupNamed:@"mate_info" error:error];
    if (!mateGrp) return NO;
    id<TTIOStorageDataset> ds = [mateGrp createDatasetNamed:@"inline_v2"
                                                  precision:TTIOPrecisionUInt8
                                                     length:blob.length
                                                  chunkSize:65536
                                                compression:TTIOCompressionNone
                                           compressionLevel:0
                                                      error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeAll:blob error:error]) return NO;
    if (![ds setAttributeValue:@((uint8_t)TTIOCompressionMateInlineV2)
                       forName:@"compression"
                         error:error]) return NO;
    NSArray *vlNameField = @[
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString]
    ];
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:chromNames.count];
    for (NSString *cn in chromNames) [rows addObject:@{@"name": cn}];
    id<TTIOStorageDataset> namesDs = [mateGrp createCompoundDatasetNamed:@"chrom_names"
                                                                     fields:vlNameField
                                                                      count:chromNames.count
                                                                      error:error];
    if (!namesDs || ![namesDs writeAll:rows error:error]) return NO;
    return YES;
}

static BOOL _TTIO_PhaseT_WriteReadNamesBulkStorage(
        id<TTIOStorageGroup> sc, NSData *blob, NSError **error)
{
    id<TTIOStorageDataset> ds = [sc createDatasetNamed:@"read_names"
                                             precision:TTIOPrecisionUInt8
                                                length:blob.length
                                             chunkSize:65536
                                           compression:TTIOCompressionNone
                                      compressionLevel:0
                                                 error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeAll:blob error:error]) return NO;
    if (![ds setAttributeValue:@((uint8_t)TTIOCompressionNameTokenizedV2)
                       forName:@"compression"
                         error:error]) return NO;
    return YES;
}

static BOOL _TTIO_PhaseT_WriteRefDiffV2BulkStorage(
        id<TTIOStorageGroup> sc, NSData *blob, NSError **error)
{
    id<TTIOStorageGroup> seqGrp = [sc createGroupNamed:@"sequences" error:error];
    if (!seqGrp) return NO;
    id<TTIOStorageDataset> ds = [seqGrp createDatasetNamed:@"refdiff_v2"
                                                 precision:TTIOPrecisionUInt8
                                                    length:blob.length
                                                 chunkSize:65536
                                               compression:TTIOCompressionNone
                                          compressionLevel:0
                                                     error:error];
    if (!ds) return NO;
    if (blob.length > 0 && ![ds writeAll:blob error:error]) return NO;
    if (![ds setAttributeValue:@((uint8_t)TTIOCompressionRefDiffV2)
                       forName:@"compression"
                         error:error]) return NO;
    return YES;
}

static BOOL _TTIO_V17_WriteMateInfoInlineV2HDF5(TTIOHDF5Group *sc,
                                                  TTIOWrittenGenomicRun *run,
                                                  NSError **error)
{
    NSData *ownChromIds = nil, *mateChromIds = nil;
    NSMutableArray<NSString *> *chromNames = nil;
    if (!_TTIO_V17_BuildChromTables(run, &ownChromIds, &mateChromIds, &chromNames)) {
        if (error) *error = [NSError errorWithDomain:@"TTIOSpectralDatasetErrorDomain"
                                                code:2100
                                            userInfo:@{NSLocalizedDescriptionKey:
                                                @"v1.7 inline_v2: chrom table build failed"}];
        return NO;
    }

    NSError *encErr = nil;
    NSData *blob = _TTIO_M86_EncodeMateInfoViaRegistry(mateChromIds,
                                                       run.matePositionsData,
                                                       run.templateLengthsData,
                                                       ownChromIds,
                                                       run.positionsData,
                                                       &encErr);
    if (!blob) {
        if (error) *error = encErr ?: [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2101
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"v1.7 inline_v2: TTIOMateInfoV2 encode failed"}];
        return NO;
    }

    TTIOHDF5Group *mateGrp = [sc createGroupNamed:@"mate_info" error:error];
    if (!mateGrp) return NO;

    // Write the inline_v2 blob as uint8 dataset with @compression = 13.
    TTIOHDF5Dataset *ds = [mateGrp createDatasetNamed:@"inline_v2"
                                             precision:TTIOPrecisionUInt8
                                                length:blob.length
                                             chunkSize:65536
                                           compression:TTIOCompressionNone
                                      compressionLevel:0
                                                 error:error];
    if (!ds) return NO;
    if (![ds writeData:blob error:error]) return NO;
    if (!_TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                       (uint8_t)TTIOCompressionMateInlineV2,
                                       error)) return NO;

    // Write the chrom_names compound sidecar (name:VL_STRING), row index = chrom_id.
    NSArray *vlNameField = @[
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString]
    ];
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:chromNames.count];
    for (NSString *cn in chromNames) [rows addObject:@{@"name": cn}];
    if (![TTIOCompoundIO writeGeneric:rows
                             intoGroup:mateGrp datasetNamed:@"chrom_names"
                                 fields:vlNameField error:error]) return NO;

    return YES;
}

/** Storage-protocol path: twin of _TTIO_V17_WriteMateInfoInlineV2HDF5. */
static BOOL _TTIO_V17_WriteMateInfoInlineV2Storage(id<TTIOStorageGroup> sc,
                                                    TTIOWrittenGenomicRun *run,
                                                    NSMutableDictionary<NSString *, NSNumber *> *shared,
                                                    NSError **error)
{
    NSData *ownChromIds = nil, *mateChromIds = nil;
    NSMutableArray<NSString *> *chromNames = nil;
    if (!_TTIO_V17_BuildChromTablesShared(run, shared, &ownChromIds, &mateChromIds, &chromNames)) {
        if (error) *error = [NSError errorWithDomain:@"TTIOSpectralDatasetErrorDomain"
                                                code:2100
                                            userInfo:@{NSLocalizedDescriptionKey:
                                                @"v1.7 inline_v2: chrom table build failed"}];
        return NO;
    }

    NSError *encErr = nil;
    NSData *blob = _TTIO_M86_EncodeMateInfoViaRegistry(mateChromIds,
                                                       run.matePositionsData,
                                                       run.templateLengthsData,
                                                       ownChromIds,
                                                       run.positionsData,
                                                       &encErr);
    if (!blob) {
        if (error) *error = encErr ?: [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2101
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"v1.7 inline_v2: TTIOMateInfoV2 encode failed"}];
        return NO;
    }

    id<TTIOStorageGroup> mateGrp = [sc createGroupNamed:@"mate_info" error:error];
    if (!mateGrp) return NO;

    // Write the inline_v2 blob.
    id<TTIOStorageDataset> ds = [mateGrp createDatasetNamed:@"inline_v2"
                                                  precision:TTIOPrecisionUInt8
                                                     length:blob.length
                                                  chunkSize:65536
                                                compression:TTIOCompressionNone
                                           compressionLevel:0
                                                      error:error];
    if (!ds) return NO;
    if (![ds writeAll:blob error:error]) return NO;
    if (![ds setAttributeValue:@((uint8_t)TTIOCompressionMateInlineV2)
                        forName:@"compression"
                          error:error]) return NO;

    // Write the chrom_names compound sidecar.
    NSArray *vlNameField = @[
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString]
    ];
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:chromNames.count];
    for (NSString *cn in chromNames) [rows addObject:@{@"name": cn}];
    id<TTIOStorageDataset> namesDs =
        [mateGrp createCompoundDatasetNamed:@"chrom_names"
                                      fields:vlNameField
                                       count:rows.count
                                       error:error];
    if (!namesDs || ![namesDs writeAll:rows error:error]) return NO;

    return YES;
}
#pragma mark - HDF5 write (flat-buffer fast path)

// ── ref-diff v2 reference-embed + write helpers ────────────────────

/** Per-run v1.5-candidacy check for the qualities auto-default gate.
 *  YES when any explicit override on the run is a v1.5 codec, OR when
 *  the ref-diff path will auto-apply to sequences (reference + signal
 *  compression "gzip" + no override). */
static BOOL _TTIO_M94_RunIsV15Candidate(TTIOWrittenGenomicRun *run)
{
    for (NSNumber *codecBox in [run.signalCodecOverrides objectEnumerator]) {
        TTIOCompression codec =
            (TTIOCompression)[codecBox unsignedIntegerValue];
        if (codec == TTIOCompressionFqzcompNx16Z) return YES;
        if (codec == TTIOCompressionDeltaRansOrder0) return YES;
    }
    // ref-diff path auto-applies to sequences → v1.5 candidate.
    if (run.signalCodecOverrides[@"sequences"] == nil
        && run.signalCompression == TTIOCompressionZlib
        && run.referenceChromSeqs != nil) {
        return YES;
    }
    return NO;
}

/** Compute the canonical reference MD5 for a single run: MD5 of the
 *  concatenation of chromosome sequences in sorted-name order. Returns
 *  empty data when ``referenceChromSeqs`` is nil. */
static NSData *_TTIO_M93_ReferenceMD5ForRun(TTIOWrittenGenomicRun *run)
{
    if (run.referenceChromSeqs == nil) return [NSData data];
    if ([run.referenceChromSeqs isKindOfClass:[TTIOLazyReference class]]) {
        return [(TTIOLazyReference *)run.referenceChromSeqs setMD5];  /* cached whole-FASTA digest */
    }
    NSArray *names = [[run.referenceChromSeqs allKeys]
        sortedArrayUsingSelector:@selector(compare:)];
    uint8_t digest[16];
    MD5_CTX c; MD5_Init(&c);
    for (NSString *name in names) {
        NSData *seq = run.referenceChromSeqs[name];
        MD5_Update(&c, seq.bytes, seq.length);
    }
    MD5_Final(digest, &c);
    return [NSData dataWithBytes:digest length:16];
}

/** Embed each unique reference (keyed by reference_uri) once at
 *  ``/study/references/<uri>/``. Per-run dedup follows Q6 = C: same
 *  URI carrying two different MD5s in one file is a hard error.
 *
 *  Embedding is pure HDF5 I/O — it does not require
 *  ``libttio_rans``. The native lib is needed only by the
 *  signal-channel REF_DIFF_V2 encode path, which is gated
 *  separately downstream (``_TTIO_V18_UseRefDiffV2``). Phase 0 Task
 *  0.11 (tio-browser): firing on ``embedReference=YES`` plus a
 *  non-nil ``referenceChromSeqs`` matches the spirit of Python's
 *  writer — embed regardless of the encode-side codec
 *  availability. */
static BOOL _TTIO_M93_EmbedReferences(TTIOHDF5Group *study,
                                       NSDictionary *genomicRuns,
                                       NSError **error)
{
    NSMutableDictionary<NSString *, NSData *> *needsEmbedMD5 =
        [NSMutableDictionary dictionary];
    NSMutableDictionary<NSString *, NSDictionary<NSString *, NSData *> *> *needsEmbedSeqs =
        [NSMutableDictionary dictionary];

    for (TTIOWrittenGenomicRun *run in [genomicRuns objectEnumerator]) {
        if (!run.embedReference) continue;
        if (run.referenceChromSeqs == nil) continue;
        // Phase 0 Task 0.11: native-lib gate removed. Writing the
        // chromosome bytes themselves is pure HDF5 I/O; the
        // REF_DIFF_V2 encode-side gate stays at its call site
        // (``_TTIO_V18_UseRefDiffV2``).

        NSData *md5 = _TTIO_M93_ReferenceMD5ForRun(run);
        NSString *uri = run.referenceUri ?: @"";
        NSData *existing = needsEmbedMD5[uri];
        if (existing) {
            if (![existing isEqualToData:md5]) {
                if (error) *error = [NSError
                    errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2200
                           userInfo:@{NSLocalizedDescriptionKey:
                               [NSString stringWithFormat:
                                    @"reference_uri '%@' carries two different "
                                    @"MD5s across runs in this dataset — same "
                                    @"URI cannot map to two different reference "
                                    @"contents.", uri]}];
                return NO;
            }
            continue;
        }
        needsEmbedMD5[uri] = md5;
        needsEmbedSeqs[uri] = run.referenceChromSeqs;
    }

    if (needsEmbedMD5.count == 0) return YES;

    TTIOHDF5Group *refsG = nil;
    if ([study hasChildNamed:@"references"]) {
        refsG = [study openGroupNamed:@"references" error:error];
    } else {
        refsG = [study createGroupNamed:@"references" error:error];
    }
    if (!refsG) return NO;

    NSArray *uris = [[needsEmbedMD5 allKeys]
        sortedArrayUsingSelector:@selector(compare:)];
    for (NSString *uri in uris) {
        if ([refsG hasChildNamed:uri]) continue;
        TTIOHDF5Group *refG = [refsG createGroupNamed:uri error:error];
        if (!refG) return NO;

        // md5 stored as hex ASCII string attribute.
        NSData *md5 = needsEmbedMD5[uri];
        const uint8_t *p = (const uint8_t *)md5.bytes;
        NSMutableString *hex = [NSMutableString stringWithCapacity:32];
        for (int i = 0; i < 16; i++) [hex appendFormat:@"%02x", p[i]];
        if (![refG setStringAttribute:@"md5" value:hex error:error]) return NO;
        if (![refG setStringAttribute:@"reference_uri"
                                 value:uri error:error]) return NO;

        TTIOHDF5Group *chromsG = [refG createGroupNamed:@"chromosomes" error:error];
        if (!chromsG) return NO;

        NSDictionary<NSString *, NSData *> *seqs = needsEmbedSeqs[uri];
        NSArray *cnames = [[seqs allKeys]
            sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *cname in cnames) {
            TTIOHDF5Group *cg = [chromsG createGroupNamed:cname error:error];
            if (!cg) return NO;
            NSData *seq = seqs[cname];
            if (![cg setIntegerAttribute:@"length"
                                    value:(int64_t)seq.length error:error]) return NO;
            /* data_packed when packing wins, raw data otherwise (same
               dispatch as TTIOReferenceImport writeToDataset). */
            NSString *dsName = nil;
            NSData *payload =
                [TTIOPackedReference payloadForSequence:seq datasetName:&dsName];
            TTIOHDF5Dataset *ds =
                [cg createDatasetNamed:dsName
                              precision:TTIOPrecisionUInt8
                                 length:payload.length
                              chunkSize:65536
                            compression:TTIOCompressionZlib
                       compressionLevel:6
                                  error:error];
            if (!ds) return NO;
            if (![ds writeData:payload error:error]) return NO;
        }
    }
    return YES;
}

/** Resolve the per-run chromosome (REF_DIFF v1.2 supports single-chrom
 *  runs only). Returns the chromosome bytes or nil with no error if the
 *  run is multi-chrom or has no covering reference, *and* sets *outChrom
 *  to the chromosome name when found. */
static NSData *_TTIO_M93_ResolveSingleChromForRun(TTIOWrittenGenomicRun *run,
                                                    NSString **outChrom,
                                                    NSError **error)
{
    if (run.referenceChromSeqs == nil) return nil;
    NSMutableSet<NSString *> *unique = [NSMutableSet set];
    for (NSString *c in run.chromosomes) {
        if (c.length > 0) [unique addObject:c];
    }
    if (unique.count == 0) return nil;
    if (unique.count > 1) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2210
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"REF_DIFF v1.2 first pass supports single-"
                            @"chromosome runs only; this run carries reads on "
                            @"%lu chromosomes. Multi-chromosome support is an "
                            @"M93.X follow-up — split into per-chromosome "
                            @"runs as a workaround.",
                            (unsigned long)unique.count]}];
        return nil;
    }
    NSString *chrom = [unique anyObject];
    NSData *seq = run.referenceChromSeqs[chrom];
    if (!seq) return nil;
    if (outChrom) *outChrom = chrom;
    return seq;
}

// v1 REF_DIFF (codec id 9) writer path removed.
// _TTIO_M93_WriteRefDiffSequences / _TTIO_M93_PerReadSequences /
// _TTIO_M93_DefaultSequencesCodec deleted alongside it. The default
// reference-aware path is now refdiff_v2 (id 14) — see the
// _TTIO_V18_* helpers below. When v2 is not eligible (no native lib,
// or no reference) the writer falls through to the
// generic byte-channel path with the run's signalCompression default
// or the explicit override (RANS / BASE_PACK).

// ── ref_diff v2 writer helpers ────────────────────────────────────
//
// When eligible (native lib + reference),
// the writer creates signal_channels/sequences as a GROUP containing a
// single child dataset "refdiff_v2" carrying the encoded blob
// @compression=14.  On opt-out or ineligibility the existing
// _TTIO_M93_WriteRefDiffSequences path is used unchanged.

/** Returns YES when the v2 path should be used for sequences.
 *  (v1.0 reset: opt-out flag removed.) */
static BOOL _TTIO_V18_UseRefDiffV2(TTIOWrittenGenomicRun *run)
{
    if (![TTIORefDiffV2 nativeAvailable]) return NO;
    if (run.referenceChromSeqs == nil) return NO;
    /* Unmapped reads (cigar "*") are carried by the codec since v1.9
     * (soft-clip bases plus the slice UL substream). */
    return YES;
}

/** Build the n_reads+1 offsets array required by TTIORefDiffV2.encode.
 *  run.offsetsData contains n_reads uint64 entries; we append
 *  offsets[n-1] + lengths[n-1] as the sentinel. */
static NSData *_TTIO_V18_BuildOffsetsPlusOne(TTIOWrittenGenomicRun *run)
{
    NSUInteger n = run.readCount;
    NSMutableData *out =
        [NSMutableData dataWithLength:(n + 1) * sizeof(uint64_t)];
    const uint64_t *src = (const uint64_t *)run.offsetsData.bytes;
    uint64_t *dst = (uint64_t *)out.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) dst[i] = src[i];
    if (n > 0) {
        const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
        dst[n] = src[n - 1] + (uint64_t)lens[n - 1];
    } else {
        dst[0] = 0;
    }
    return out;
}

/** Encode REF_DIFF_V2 sequences through the codec registry. Builds the
 *  encode context exactly as the writers sourced their args from the
 *  direct +[TTIORefDiffV2 encodeSequences:...] call: offsets (n+1),
 *  positions = run.positionsData, cigars = run.cigars (via lazy
 *  provider), reference = the resolved single-chrom sequence,
 *  referenceMd5 = the run's reference MD5, referenceUri =
 *  run.referenceUri ?: @"", readsPerSlice = 10000 (the codec's default
 *  when ctx.readsPerSlice is nil). Returns the refdiff_v2 child blob,
 *  byte-identical to the prior direct call, or nil + error on failure
 *  (caller falls back to BASE_PACK). */
static NSData *_TTIO_V18_EncodeRefDiffV2ViaRegistry(TTIOWrittenGenomicRun *run,
                                                    NSData *offsets,
                                                    NSData *reference,
                                                    NSData *referenceMd5,
                                                    NSError **error)
{
    id<TTIOCodec> c =
        [TTIOCodecRegistry codecForId:TTIOCompressionRefDiffV2];
    if (c == nil) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2120
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"REF_DIFF_V2 codec not registered"}];
        return nil;
    }
    TTIOCodecContext *ctx = [TTIOCodecContext emptyContext];
    ctx.offsets = offsets;
    ctx.positions = run.positionsData;
    NSArray<NSString *> *cigars = run.cigars;
    ctx.cigarsProvider = ^NSArray<NSString *> *{ return cigars; };
    ctx.reference = reference;
    ctx.referenceMd5 = referenceMd5;
    ctx.referenceUri = run.referenceUri ?: @"";
    // ctx.readsPerSlice left nil -> codec uses its 10000 default,
    // matching the prior direct call's readsPerSlice:10000.
    if (run.refDiffSliceBytes > 0)
        ctx.sliceBytes = @(run.refDiffSliceBytes);
    TTIOEncodedChannel *enc =
        [c encode:[[TTIODecodedBytes alloc] initWithData:run.sequencesData]
           context:ctx
             error:error];
    if (![enc isKindOfClass:[TTIOEncodedGroupLayout class]]) return nil;
    return ((TTIOEncodedGroupLayout *)enc).children[@"refdiff_v2"];
}

/** HDF5 fast path: write signal_channels/sequences as a GROUP with a
 *  refdiff_v2 child dataset @compression=14.  Falls back (with a silent
 *  BASE_PACK path) via _TTIO_M93_WriteRefDiffSequences when the C
 *  encoder returns nil. */
static BOOL _TTIO_V18_WriteRefDiffV2SequencesHDF5(TTIOHDF5Group *sc,
                                                    TTIOWrittenGenomicRun *run,
                                                    NSError **error)
{
    // v1 REF_DIFF fallback removed. When v2 is
    // ineligible (no covering reference, encode failure) we fall back
    // to BASE_PACK on the flat sequences buffer — same content the v1
    // path used as its own fallback.
    NSString *chrom = nil;
    NSData *chromSeq =
        _TTIO_M93_ResolveSingleChromForRun(run, &chrom, error);
    if (!chromSeq && error && *error) return NO;
    if (chromSeq == nil) {
        return _TTIO_M86_WriteByteChannel(sc, @"sequences",
                                           run.sequencesData,
                                           TTIOCompressionBasePack,
                                           @(TTIOCompressionBasePack),
                                           error);
    }

    NSData *offsets = _TTIO_V18_BuildOffsetsPlusOne(run);
    NSData *md5 = _TTIO_M93_ReferenceMD5ForRun(run);

    NSError *encErr = nil;
    NSData *encoded =
        _TTIO_V18_EncodeRefDiffV2ViaRegistry(run, offsets, chromSeq, md5, &encErr);
    if (!encoded) {
        // v2 encode failed — fall back to BASE_PACK so the write still
        // succeeds (no v1 REF_DIFF path under v1.0).
        return _TTIO_M86_WriteByteChannel(sc, @"sequences",
                                           run.sequencesData,
                                           TTIOCompressionBasePack,
                                           @(TTIOCompressionBasePack),
                                           error);
    }

    // Create signal_channels/sequences as a GROUP.
    TTIOHDF5Group *seqGrp = [sc createGroupNamed:@"sequences" error:error];
    if (!seqGrp) return NO;

    TTIOHDF5Dataset *ds = [seqGrp createDatasetNamed:@"refdiff_v2"
                                            precision:TTIOPrecisionUInt8
                                               length:encoded.length
                                            chunkSize:65536
                                          compression:TTIOCompressionNone
                                     compressionLevel:0
                                                error:error];
    if (!ds) return NO;
    if (![ds writeData:encoded error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)TTIOCompressionRefDiffV2,
                                         error);
}

/** Storage-protocol path: twin of _TTIO_V18_WriteRefDiffV2SequencesHDF5. */
static BOOL _TTIO_V18_WriteRefDiffV2SequencesStorage(id<TTIOStorageGroup> sc,
                                                       TTIOWrittenGenomicRun *run,
                                                       NSData *precomputedMD5,
                                                       NSError **error)
{
    NSString *chrom = nil;
    NSData *chromSeq =
        _TTIO_M93_ResolveSingleChromForRun(run, &chrom, error);
    if (!chromSeq && error && *error) return NO;
    if (chromSeq == nil) {
        // No reference: fall through to the existing BASE_PACK path via
        // _TTIO_M86_WriteByteChannelStorage with no override.
        return _TTIO_M86_WriteByteChannelStorage(sc, @"sequences",
                                                  run.sequencesData,
                                                  TTIOCompressionBasePack,
                                                  @(TTIOCompressionBasePack),
                                                  error);
    }

    NSData *offsets = _TTIO_V18_BuildOffsetsPlusOne(run);
    NSData *md5 = precomputedMD5 ?: _TTIO_M93_ReferenceMD5ForRun(run);

    NSError *encErr = nil;
    NSData *encoded =
        _TTIO_V18_EncodeRefDiffV2ViaRegistry(run, offsets, chromSeq, md5, &encErr);
    if (!encoded) {
        // Fall back to plain sequences with BASE_PACK.
        return _TTIO_M86_WriteByteChannelStorage(sc, @"sequences",
                                                  run.sequencesData,
                                                  TTIOCompressionBasePack,
                                                  @(TTIOCompressionBasePack),
                                                  error);
    }

    // Create signal_channels/sequences as a GROUP.
    id<TTIOStorageGroup> seqGrp = [sc createGroupNamed:@"sequences" error:error];
    if (!seqGrp) return NO;

    id<TTIOStorageDataset> ds =
        [seqGrp createDatasetNamed:@"refdiff_v2"
                          precision:TTIOPrecisionUInt8
                             length:encoded.length
                          chunkSize:65536
                        compression:TTIOCompressionNone
                   compressionLevel:0
                              error:error];
    if (!ds) return NO;
    if (![ds writeAll:encoded error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)TTIOCompressionRefDiffV2)
                         forName:@"compression"
                           error:error];
}

/** v1.5 default codec for qualities (M94 Q5a=B): when caller supplied
 *  no override on qualities, signal_compression is the gzip default,
 *  AND the run is a v1.5 candidate (per _TTIO_M94_RunIsV15Candidate),
 *  return FQZCOMP_NX16_Z. The v1.5-candidacy gate preserves byte-parity
 *  with M82-only writes that don't use any v1.5 codec — those keep the
 *  legacy uncompressed-qualities path. */
static NSNumber *_TTIO_M94_DefaultQualitiesCodec(TTIOWrittenGenomicRun *run)
{
    if (run.signalCodecOverrides[@"qualities"] != nil) return nil;
    if (run.signalCompression != TTIOCompressionZlib) return nil;
    if (!_TTIO_M94_RunIsV15Candidate(run)) return nil;
    return @(TTIOCompressionFqzcompNx16Z);
}


/** SEQ_CM (codec id 19, M103) over the run's sequences, given the read
 *  lengths from run.lengthsData (uint32 LE). Mirrors Python's
 *  ``_write_sequences_seq_cm``. */
static NSData *_TTIO_M103_EncodeSeqCm(TTIOWrittenGenomicRun *run, NSError **error)
{
    NSUInteger n = run.lengthsData.length / sizeof(uint32_t);
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    NSMutableArray *readLengths = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) [readLengths addObject:@(lens[i])];
    TTIOCodecContext *c = [TTIOCodecContext emptyContext];
    c.readLengths = readLengths;
    TTIOEncodedChannel *enc =
        [[TTIOCodecRegistry codecForId:TTIOCompressionSeqCm]
            encode:[[TTIODecodedBytes alloc] initWithData:run.sequencesData]
           context:c
             error:error];
    return enc ? ((TTIOEncodedDatasetBytes *)enc).bytes : nil;
}

static BOOL _TTIO_M103_WriteSequencesSeqCm(TTIOHDF5Group *sc, TTIOWrittenGenomicRun *run,
                                           NSError **error)
{
    NSData *encoded = _TTIO_M103_EncodeSeqCm(run, error);
    if (!encoded) return NO;
    TTIOHDF5Dataset *ds = [sc createDatasetNamed:@"sequences"
                                        precision:TTIOPrecisionUInt8
                                           length:encoded.length
                                        chunkSize:65536
                                      compression:TTIOCompressionNone
                                 compressionLevel:0
                                            error:error];
    if (!ds) return NO;
    if (![ds writeData:encoded error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)TTIOCompressionSeqCm, error);
}

/** Storage-protocol twin of _TTIO_M103_WriteSequencesSeqCm. */
static BOOL _TTIO_M103_WriteSequencesSeqCmStorage(id<TTIOStorageGroup> sc,
                                                  TTIOWrittenGenomicRun *run,
                                                  NSError **error)
{
    NSData *encoded = _TTIO_M103_EncodeSeqCm(run, error);
    if (!encoded) return NO;
    id<TTIOStorageDataset> ds = [sc createDatasetNamed:@"sequences"
                                              precision:TTIOPrecisionUInt8
                                                 length:encoded.length
                                              chunkSize:65536
                                            compression:TTIOCompressionNone
                                       compressionLevel:0
                                                  error:error];
    if (!ds) return NO;
    if (![ds writeAll:encoded error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)TTIOCompressionSeqCm)
                         forName:@"compression"
                           error:error];
}

/** Write the qualities channel through FQZCOMP_NX16_Z. Derives
 *  read_lengths from run.lengthsData (uint32 LE) and revcomp_flags
 *  from run.flagsData[i] & 16 (SAM REVERSE bit). Stamps the
 *  @compression attribute with codec id 12. Mirrors Python's
 *  ``_write_qualities_fqzcomp_nx16z``. */
static BOOL _TTIO_M94Z_WriteQualitiesFqzcompNx16Z(TTIOHDF5Group *sc,
                                                    TTIOWrittenGenomicRun *run,
                                                    NSError **error)
{
    NSUInteger n = run.lengthsData.length / sizeof(uint32_t);
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    const uint32_t *flgs = (const uint32_t *)run.flagsData.bytes;
    NSMutableArray *readLengths = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *revcompFlags = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        [readLengths addObject:@(lens[i])];
        uint32_t f = (run.flagsData.length >= (i + 1) * sizeof(uint32_t))
                       ? flgs[i] : 0;
        [revcompFlags addObject:(f & 16u) ? @1 : @0];
    }
    TTIOCodecContext *fqzCtx = [TTIOCodecContext emptyContext];
    fqzCtx.readLengths = readLengths;
    fqzCtx.revcompFlags = revcompFlags;
    /* Qualities V5 gate (spec 2.4): offer the base bytes to the
     * encoder only when the run carries a base-parallel sequences
     * channel and the caller did not opt out; V4 still wins by exact
     * size wherever sequence context does not pay. */
    if (!run.optDisableQualitiesV5
        && run.sequencesData.length == run.qualitiesData.length) {
        fqzCtx.sequences = run.sequencesData;
    }
    id<TTIOCodec> fqzCodec =
        [TTIOCodecRegistry codecForId:TTIOCompressionFqzcompNx16Z];
    TTIOEncodedChannel *fqzEnc =
        [fqzCodec encode:[[TTIODecodedBytes alloc] initWithData:run.qualitiesData]
                 context:fqzCtx
                   error:error];
    NSData *encoded = ((TTIOEncodedDatasetBytes *)fqzEnc).bytes;
    if (!encoded) return NO;
    TTIOHDF5Dataset *ds = [sc createDatasetNamed:@"qualities"
                                        precision:TTIOPrecisionUInt8
                                           length:encoded.length
                                        chunkSize:65536
                                      compression:TTIOCompressionNone
                                 compressionLevel:0
                                            error:error];
    if (!ds) return NO;
    if (![ds writeData:encoded error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)TTIOCompressionFqzcompNx16Z,
                                         error);
}

/** Storage-protocol twin of _TTIO_M94Z_WriteQualitiesFqzcompNx16Z.
 *  qualStrategyHint: -1 auto (the 3-way tune), 5/6 forced V5,
 *  TTIOM94ZHintV4Auto V4 with internal preset selection; the stream
 *  writer passes its per-run pin through the write context. */
static BOOL _TTIO_M94Z_WriteQualitiesFqzcompNx16ZStorage(id<TTIOStorageGroup> sc,
                                                           TTIOWrittenGenomicRun *run,
                                                           NSInteger qualStrategyHint,
                                                           NSError **error)
{
    NSUInteger n = run.lengthsData.length / sizeof(uint32_t);
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    const uint32_t *flgs = (const uint32_t *)run.flagsData.bytes;
    NSMutableArray *readLengths = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *revcompFlags = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        [readLengths addObject:@(lens[i])];
        uint32_t f = (run.flagsData.length >= (i + 1) * sizeof(uint32_t)) ? flgs[i] : 0;
        [revcompFlags addObject:(f & 16u) ? @1 : @0];
    }
    TTIOCodecContext *fqzCtx = [TTIOCodecContext emptyContext];
    fqzCtx.readLengths = readLengths;
    fqzCtx.revcompFlags = revcompFlags;
    fqzCtx.qualStrategyHint = @(qualStrategyHint);
    if (!run.optDisableQualitiesV5
        && run.sequencesData.length == run.qualitiesData.length) {
        fqzCtx.sequences = run.sequencesData;
    }
    id<TTIOCodec> fqzCodec = [TTIOCodecRegistry codecForId:TTIOCompressionFqzcompNx16Z];
    TTIOEncodedChannel *fqzEnc =
        [fqzCodec encode:[[TTIODecodedBytes alloc] initWithData:run.qualitiesData]
                 context:fqzCtx
                   error:error];
    NSData *encoded = ((TTIOEncodedDatasetBytes *)fqzEnc).bytes;
    if (!encoded) return NO;
    id<TTIOStorageDataset> ds = [sc createDatasetNamed:@"qualities"
                                              precision:TTIOPrecisionUInt8
                                                 length:encoded.length
                                              chunkSize:65536
                                            compression:TTIOCompressionNone
                                       compressionLevel:0
                                                  error:error];
    if (!ds) return NO;
    if (![ds writeAll:encoded error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)TTIOCompressionFqzcompNx16Z)
                         forName:@"compression"
                           error:error];
}

/* Write an index array as a 1-D HDF5 dataset matching what
 * TTIOSpectrumIndex -writeToGroup:error: emits (same precision,
 * chunkSize=1024, compression level 6). The format is load-bearing:
 * readers — including Java and Python — depend on exactly this
 * layout. */
static BOOL writeIndexArrayDS(TTIOHDF5Group *g, NSString *name,
                               TTIOPrecision p, NSData *data,
                               NSError **error)
{
    if (!data) return YES;
    NSUInteger n = data.length / TTIOPrecisionElementSize(p);
    TTIOHDF5Dataset *ds = [g createDatasetNamed:name
                                       precision:p
                                          length:n
                                       chunkSize:4096
                                compressionLevel:6
                                           error:error];
    if (!ds) return NO;
    return [ds writeData:data error:error];
}

// M102: a TTIOWrittenRun's pixel columns are all-or-none, each
// int32[spectrum count]. Sets *outAny when the run carries them.
static BOOL validateWrittenRunPixels(TTIOWrittenRun *run, NSString *name,
                                     BOOL *outAny, NSError **error)
{
    if (![run hasPixelCoordinates]) return YES;
    NSUInteger n = run.offsets.length / sizeof(int64_t);
    NSUInteger want = n * sizeof(int32_t);
    if (!run.pixelX || !run.pixelY || !run.pixelZ ||
        run.pixelX.length != want || run.pixelY.length != want ||
        run.pixelZ.length != want) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"run '%@': pixelX/pixelY/pixelZ must all be set, each "
            @"holding %lu int32 values", name, (unsigned long)n);
        return NO;
    }
    if (outAny) *outAny = YES;
    return YES;
}

// Task 30: provider-agnostic 1-D index array writer. Mirrors
// writeIndexArrayDS but speaks the StorageGroup protocol so it works
// for memory:// / sqlite:// / zarr:// targets. Layout matches the HDF5
// fast path (precision, length); chunkSize/compression are honoured by
// HDF5 only and harmless on other backends per Appendix B Gap 3.
//
// The ``compression`` arg is the codec the caller WOULD use against an
// HDF5 target. The ObjC Zarr provider rejects any non-None compression
// (rather than silently ignoring it like Memory/SQLite do), so the
// caller must downgrade to TTIOCompressionNone for that backend; see
// task30CompressionForProvider() below.
static BOOL writeIndexArrayStorage(id<TTIOStorageGroup> g, NSString *name,
                                    TTIOPrecision p, NSData *data,
                                    TTIOCompression compression,
                                    NSError **error)
{
    if (!data) return YES;
    NSUInteger n = data.length / TTIOPrecisionElementSize(p);
    id<TTIOStorageDataset> ds = [g createDatasetNamed:name
                                            precision:p
                                               length:n
                                            chunkSize:4096
                                          compression:compression
                                     compressionLevel:6
                                                error:error];
    if (!ds) return NO;
    return [ds writeAll:data error:error];
}

// Task 30: choose a write-time compression that the provider accepts.
// Memory/SQLite providers ignore the compression argument; HDF5 honours
// it; Zarr (ObjC) rejects anything non-None at the dataset-create step.
// This helper centralises the downgrade so per-call sites stay flat.
static TTIOCompression task30CompressionForProvider(id<TTIOStorageProvider> p)
{
    if ([p respondsToSelector:@selector(supportsCompression)]
        && [p supportsCompression]) {
        return TTIOCompressionZlib;
    }
    return TTIOCompressionNone;
}

// ── M101 SAM optional fields ──────────────────────────────────────
//
// signal_channels/tags: one flat uint8 dataset, no HDF5 filter,
// @compression = SAM_TAGS (18) as uint8. Written only when some read
// carries tags, so tag-less files are unchanged.

/** Encode the run's tags through SAM_TAGS (codec id 18).
 *
 *  MD:Z / NM:i are recomputed from the reference only when the
 *  sequences channel went through REF_DIFF_V2 (`derive`), so a reader
 *  that can decode the sequences can always decode the tags (binding
 *  decision 96). REF_DIFF_V2 runs hold one chromosome: reads on it get
 *  the reference, every other read (unmapped '*') none. */
static NSData *_TTIO_M101_EncodeTags(TTIOWrittenGenomicRun *run, BOOL derive,
                                     NSError **error)
{
    NSUInteger n = run.cigars.count;
    if (run.tags.count != n) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"WrittenGenomicRun.tags has %lu entries for %lu reads",
            (unsigned long)run.tags.count, (unsigned long)n);
        return nil;
    }
    NSArray *refs = @[];
    NSMutableData *chromIds = [NSMutableData dataWithLength:n * sizeof(uint16_t)];
    uint16_t *cid = (uint16_t *)chromIds.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) cid[i] = 0xFFFF;
    if (derive && run.referenceChromSeqs.count > 0) {
        NSMutableSet<NSString *> *names = [NSMutableSet set];
        for (NSString *c in run.chromosomes) {
            if (c.length > 0 && ![c isEqualToString:@"*"]) [names addObject:c];
        }
        NSString *only = names.count == 1 ? [names anyObject] : nil;
        NSData *ref = only ? run.referenceChromSeqs[only] : nil;
        if (ref != nil) {
            refs = @[ref];
            for (NSUInteger i = 0; i < n; i++) {
                if ([run.chromosomes[i] isEqualToString:only]) cid[i] = 0;
            }
        }
    }
    NSMutableData *offsets = [NSMutableData dataWithLength:(n + 1) * sizeof(uint64_t)];
    uint64_t *off = (uint64_t *)offsets.mutableBytes;
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    for (NSUInteger i = 0; i < n; i++) off[i + 1] = off[i] + (uint64_t)lens[i];

    id<TTIOCodec> c = [TTIOCodecRegistry codecForId:TTIOCompressionSamTags];
    TTIOCodecContext *ctx = [TTIOCodecContext emptyContext];
    ctx.sequences = refs.count ? run.sequencesData : nil;
    ctx.offsets = offsets;
    NSArray<NSString *> *cigars = run.cigars;
    ctx.cigarsProvider = ^NSArray<NSString *> *{ return cigars; };
    ctx.positions = run.positionsData;
    ctx.ownChromIds = chromIds;
    ctx.tagReferences = refs;
    TTIOEncodedChannel *enc =
        [c encode:[[TTIODecodedStringList alloc] initWithNames:run.tags]
          context:ctx error:error];
    if (![enc isKindOfClass:[TTIOEncodedDatasetBytes class]]) return nil;
    return ((TTIOEncodedDatasetBytes *)enc).bytes;
}

/** HDF5 fast path: write signal_channels/tags. */
static BOOL _TTIO_M101_WriteTagsHDF5(TTIOHDF5Group *sc, TTIOWrittenGenomicRun *run,
                                     NSError **error)
{
    BOOL derive = NO;
    if ([sc hasChildNamed:@"sequences"]) {
        H5O_info2_t info;
        memset(&info, 0, sizeof(info));
        if (H5Oget_info_by_name3([sc groupId], "sequences", &info,
                                 H5O_INFO_BASIC, H5P_DEFAULT) >= 0
            && info.type == H5O_TYPE_GROUP) {
            TTIOHDF5Group *seqGrp = [sc openGroupNamed:@"sequences" error:NULL];
            derive = [seqGrp hasChildNamed:@"refdiff_v2"];
        }
    }
    NSData *blob = _TTIO_M101_EncodeTags(run, derive, error);
    if (!blob) return NO;
    TTIOHDF5Dataset *ds = [sc createDatasetNamed:@"tags"
                                       precision:TTIOPrecisionUInt8
                                          length:blob.length
                                       chunkSize:65536
                                     compression:TTIOCompressionNone
                                compressionLevel:0
                                           error:error];
    if (!ds) return NO;
    if (![ds writeData:blob error:error]) return NO;
    return _TTIO_M86_WriteUInt8Attribute([ds datasetId], "compression",
                                         (uint8_t)TTIOCompressionSamTags, error);
}

/** Storage-protocol path: twin of _TTIO_M101_WriteTagsHDF5. */
static BOOL _TTIO_M101_WriteTagsStorage(id<TTIOStorageGroup> sc,
                                        TTIOWrittenGenomicRun *run,
                                        NSError **error)
{
    id<TTIOStorageGroup> seqGrp = [sc openGroupNamed:@"sequences" error:NULL];
    BOOL derive = seqGrp != nil && [seqGrp hasChildNamed:@"refdiff_v2"];
    NSData *blob = _TTIO_M101_EncodeTags(run, derive, error);
    if (!blob) return NO;
    id<TTIOStorageDataset> ds = [sc createDatasetNamed:@"tags"
                                             precision:TTIOPrecisionUInt8
                                                length:blob.length
                                             chunkSize:65536
                                           compression:TTIOCompressionNone
                                      compressionLevel:0
                                                 error:error];
    if (!ds) return NO;
    if (![ds writeAll:blob error:error]) return NO;
    return [ds setAttributeValue:@((uint8_t)TTIOCompressionSamTags)
                         forName:@"compression"
                           error:error];
}

@interface TTIOSpectralDataset (GenomicWriteStream)
+ (BOOL)_ttio_streamGenomicRun:(TTIOWrittenGenomicRun *)run
                          name:(NSString *)name
                         study:(id<TTIOStorageGroup>)study
                         error:(NSError **)error;
@end

@implementation TTIOSpectralDataset (GenomicWrite)

+ (void)validateGenomicCodecOverridesForRun:(TTIOWrittenGenomicRun *)run
{
    _TTIO_M86_ValidateOverrides(run.signalCodecOverrides);
    _TTIO_V17_ValidateMateInfoV2Overrides(run);
    _TTIO_M97_ValidateQualityBinnedPlatform(run);
}

+ (NSData *)referenceMD5ForRun:(TTIOWrittenGenomicRun *)run
{
    return _TTIO_M93_ReferenceMD5ForRun(run);
}

+ (BOOL)embedReferencesForRuns:(NSArray<TTIOWrittenGenomicRun *> *)runs
                       inStudy:(id<TTIOStorageGroup>)study
                         error:(NSError **)error
{
    NSMutableDictionary<NSString *, NSData *> *needsEmbedMD5 = [NSMutableDictionary dictionary];
    NSMutableDictionary<NSString *, NSDictionary<NSString *, NSData *> *> *needsEmbedSeqs =
        [NSMutableDictionary dictionary];
    for (TTIOWrittenGenomicRun *run in runs) {
        if (!run.embedReference || run.referenceChromSeqs == nil) continue;
        NSData *md5 = _TTIO_M93_ReferenceMD5ForRun(run);
        NSString *uri = run.referenceUri ?: @"";
        NSData *existing = needsEmbedMD5[uri];
        if (existing) {
            if (![existing isEqualToData:md5]) {
                if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
                    @"reference_uri '%@' carries two different MD5s across runs", uri);
                return NO;
            }
            continue;
        }
        needsEmbedMD5[uri] = md5;
        needsEmbedSeqs[uri] = run.referenceChromSeqs;
    }
    if (needsEmbedMD5.count == 0) return YES;
    id<TTIOStorageGroup> refsG = [study hasChildNamed:@"references"]
        ? [study openGroupNamed:@"references" error:error]
        : [study createGroupNamed:@"references" error:error];
    if (!refsG) return NO;
    NSArray *uris = [[needsEmbedMD5 allKeys] sortedArrayUsingSelector:@selector(compare:)];
    for (NSString *uri in uris) {
        NSData *md5 = needsEmbedMD5[uri];
        const uint8_t *pb = (const uint8_t *)md5.bytes;
        NSMutableString *hex = [NSMutableString stringWithCapacity:32];
        for (int i = 0; i < 16 && i < (int)md5.length; i++) [hex appendFormat:@"%02x", pb[i]];
        if ([refsG hasChildNamed:uri]) {
            id<TTIOStorageGroup> existingG = [refsG openGroupNamed:uri error:error];
            id have = [existingG attributeValueForName:@"md5" error:NULL];
            NSString *haveHex = [have isKindOfClass:[NSData class]]
                ? [[NSString alloc] initWithData:have encoding:NSUTF8StringEncoding]
                : [have description];
            if (haveHex.length && ![haveHex isEqualToString:hex]) {
                if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
                    @"reference_uri '%@' already embedded with a different MD5", uri);
                return NO;
            }
            continue;
        }
        id<TTIOStorageGroup> refG = [refsG createGroupNamed:uri error:error];
        if (!refG) return NO;
        if (![refG setAttributeValue:hex forName:@"md5" error:error]) return NO;
        if (![refG setAttributeValue:uri forName:@"reference_uri" error:error]) return NO;
        id<TTIOStorageGroup> chromsG = [refG createGroupNamed:@"chromosomes" error:error];
        if (!chromsG) return NO;
        NSDictionary<NSString *, NSData *> *seqs = needsEmbedSeqs[uri];
        NSArray *cnames = [[seqs allKeys] sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *cname in cnames) {
            id<TTIOStorageGroup> cg = [chromsG createGroupNamed:cname error:error];
            if (!cg) return NO;
            NSData *seq = seqs[cname];
            if (![cg setAttributeValue:@((int64_t)seq.length) forName:@"length" error:error]) return NO;
            if (![TTIOPackedReference writeChromosomeDataset:cg sequence:seq error:error]) return NO;
        }
    }
    return YES;
}

// provider-agnostic write of one /study/genomic_runs/<name>/
// subtree via the StorageGroup protocol. Used by the memory:// /
// sqlite:// / zarr:// write path. The HDF5 fast path uses
// +writeGenomicRun:toGroup:name:error: instead which goes
// HDF5-direct for byte parity.
+ (BOOL)writeGenomicRunStorage:(TTIOWrittenGenomicRun *)run
                         toGroup:(id<TTIOStorageGroup>)parent
                            name:(NSString *)name
                           error:(NSError **)error
{
    return [self writeGenomicRunStorage:run toGroup:parent name:name
                                context:[TTIOGenomicWriteContext none] error:error];
}

+ (BOOL)writeGenomicRunStorage:(TTIOWrittenGenomicRun *)run
                         toGroup:(id<TTIOStorageGroup>)parent
                            name:(NSString *)name
                         context:(TTIOGenomicWriteContext *)ctx
                           error:(NSError **)error
{
    // validate signal-channel codec overrides before any
    // mutation. Same fail-fast contract as the HDF5 fast path.
    _TTIO_M86_ValidateOverrides(run.signalCodecOverrides);
    // reject per-field mate_info_* overrides when inline_v2 active.
    _TTIO_V17_ValidateMateInfoV2Overrides(run);
    _TTIO_M97_ValidateQualityBinnedPlatform(run);

    id<TTIOStorageGroup> rg = [parent createGroupNamed:name error:error];
    if (!rg) return NO;

    // Run-level attributes via the storage protocol.
    if (![rg setAttributeValue:@(run.acquisitionMode)
                         forName:@"acquisition_mode" error:error]) return NO;
    if (![rg setAttributeValue:@"genomic_sequencing"
                         forName:@"modality" error:error]) return NO;
    if (![rg setAttributeValue:@(5)
                         forName:@"spectrum_class" error:error]) return NO;
    if (![rg setAttributeValue:run.referenceUri ?: @""
                         forName:@"reference_uri" error:error]) return NO;
    if (![rg setAttributeValue:run.platform ?: @""
                         forName:@"platform" error:error]) return NO;
    if (run.readRole.length > 0
        && ![rg setAttributeValue:run.readRole
                          forName:@"read_role" error:error]) return NO;
    if (![rg setAttributeValue:run.sampleName ?: @""
                         forName:@"sample_name" error:error]) return NO;
    if (![rg setAttributeValue:@((int64_t)run.readCount)
                         forName:@"read_count" error:error]) return NO;

    // genomic_index subgroup (already provider-agnostic).
    TTIOGenomicIndex *idx = [[TTIOGenomicIndex alloc]
        initWithOffsets:run.offsetsData
                lengths:run.lengthsData
            chromosomes:run.chromosomes
              positions:run.positionsData
       mappingQualities:run.mappingQualitiesData
                  flags:run.flagsData];
    id<TTIOStorageGroup> idxG = [rg createGroupNamed:@"genomic_index" error:error];
    if (!idxG) return NO;
    if (![idx writeToGroup:idxG nameToId:ctx.chromNameToId error:error]) return NO;

    // signal_channels subgroup.
    id<TTIOStorageGroup> sc = [rg createGroupNamed:@"signal_channels" error:error];
    if (!sc) return NO;
    TTIOCompression codec = run.signalCompression;

    // positions / flags / mapping_qualities are NOT written
    // under signal_channels/. They live exclusively under
    // genomic_index/, mirroring MS's spectrum_index/ pattern. See
    // docs/format-spec.md §4 and §10.7. Override-validation rejects
    // these channel names.
    // sequences dispatch — prefer refdiff_v2
    // group layout when eligible (no override + native lib +
    // reference + all reads mapped); otherwise fall through to the
    // M86 byte-channel writer with whatever explicit override the
    // caller supplied. The v1 REF_DIFF override is no longer
    // accepted (override-validation rejects codec id 9).
    // Phase 2c-T: when bulk_v2_blobs.refDiffBlob is set, write the
    // verbatim wire blob and skip the codec encode entirely.
    TTIOBulkV2Blobs *_bulkObjC = (TTIOBulkV2Blobs *)run.bulkV2Blobs;
    {
        NSNumber *seqOvr = run.signalCodecOverrides[@"sequences"];
        if (_bulkObjC != nil && _bulkObjC.refDiffBlob != nil) {
            if (![_bulkObjC.refDiffReferenceUri isEqualToString:run.referenceUri]) {
                if (error) *error = [NSError
                    errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2110
                           userInfo:@{NSLocalizedDescriptionKey:
                               @"BulkV2Blobs.refDiffReferenceUri does not match run.referenceUri"}];
                return NO;
            }
            if (!_TTIO_PhaseT_WriteRefDiffV2BulkStorage(sc, _bulkObjC.refDiffBlob, error)) return NO;
        } else if (seqOvr == nil && _TTIO_V18_UseRefDiffV2(run)) {
            if (!_TTIO_V18_WriteRefDiffV2SequencesStorage(sc, run, ctx.referenceMD5, error)) return NO;
        } else if (seqOvr != nil
                   && (TTIOCompression)[seqOvr unsignedIntegerValue] == TTIOCompressionSeqCm) {
            if (!_TTIO_M103_WriteSequencesSeqCmStorage(sc, run, error)) return NO;
        } else {
            if (!_TTIO_M86_WriteByteChannelStorage(sc, @"sequences",
                                                   run.sequencesData, codec,
                                                   seqOvr,
                                                   error)) return NO;
        }
    }
    NSNumber *qualOvr = run.signalCodecOverrides[@"qualities"]
                        ?: _TTIO_M94_DefaultQualitiesCodec(run);
    if (qualOvr != nil
        && (TTIOCompression)[qualOvr unsignedIntegerValue] == TTIOCompressionFqzcompNx16Z) {
        if (!_TTIO_M94Z_WriteQualitiesFqzcompNx16ZStorage(sc, run,
                                                          ctx.qualStrategyHint,
                                                          error)) return NO;
    } else if (!_TTIO_M86_WriteByteChannelStorage(sc, @"qualities",
                                                  run.qualitiesData, codec,
                                                  run.signalCodecOverrides[@"qualities"],
                                                  error)) return NO;

    // 3 compound datasets via the storage-protocol's compound API.
    NSArray *vlValueField = @[
        [TTIOCompoundField fieldWithName:@"value" kind:TTIOCompoundFieldKindVLString]
    ];

    // schema lift for cigars on the provider/storage
    // path. Same dispatch as the HDF5 fast path.
    NSNumber *cigarsOverrideS = run.signalCodecOverrides[@"cigars"];
    if (cigarsOverrideS != nil) {
        TTIOCompression cigarsCodec =
            (TTIOCompression)[cigarsOverrideS unsignedIntegerValue];
        NSData *encoded = _TTIO_M86_EncodeCigarsWithCodec(run.cigars,
                                                          cigarsCodec);
        if (!encoded) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2061
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"M86 Phase C: cigars codec %lu encode "
                                @"returned nil",
                                (unsigned long)cigarsCodec]}];
            return NO;
        }
        id<TTIOStorageDataset> cigarDs = [sc createDatasetNamed:@"cigars"
                                                      precision:TTIOPrecisionUInt8
                                                         length:encoded.length
                                                      chunkSize:65536
                                                    compression:TTIOCompressionNone
                                               compressionLevel:0
                                                          error:error];
        if (!cigarDs) return NO;
        if (![cigarDs writeAll:encoded error:error]) return NO;
        if (![cigarDs setAttributeValue:@((uint8_t)cigarsCodec)
                                forName:@"compression"
                                  error:error]) return NO;
    } else {
        NSMutableArray *cigarRows = [NSMutableArray arrayWithCapacity:run.cigars.count];
        for (NSString *c in run.cigars) [cigarRows addObject:@{@"value": c}];
        id<TTIOStorageDataset> cigarDs = [sc createCompoundDatasetNamed:@"cigars"
                                                                    fields:vlValueField
                                                                     count:run.cigars.count
                                                                     error:error];
        if (!cigarDs || ![cigarDs writeAll:cigarRows error:error]) return NO;
    }

    // read_names always written via
    // NAME_TOKENIZED_V2 (codec id 15) when libttio_rans is linked.
    // The v1 NAME_TOKENIZED override (id 8) and the M82 compound
    // fallback are gone. Empty-run short-circuit writes a zero-length
    // uint8 dataset with @compression=15 + @count=0.
    if (_bulkObjC != nil && _bulkObjC.nameTokBlob != nil) {
        // Phase 2c-T: skip codec encode.
        if (!_TTIO_PhaseT_WriteReadNamesBulkStorage(sc, _bulkObjC.nameTokBlob, error)) return NO;
    } else if (run.readCount == 0) {
        id<TTIOStorageDataset> nameDs = [sc createDatasetNamed:@"read_names"
                                                     precision:TTIOPrecisionUInt8
                                                        length:0
                                                     chunkSize:1
                                                   compression:TTIOCompressionNone
                                              compressionLevel:0
                                                         error:error];
        if (!nameDs) return NO;
        if (![nameDs setAttributeValue:@((uint8_t)TTIOCompressionNameTokenizedV2)
                               forName:@"compression"
                                 error:error]) return NO;
    } else if ([TTIONameTokenizerV2 nativeAvailable]) {
        NSData *encoded = _TTIO_M86_EncodeNamesViaRegistry(run.readNames);
        if (!encoded) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2032
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"v1.0 NAME_TOKENIZED_V2 encode of "
                           @"read_names returned nil"}];
            return NO;
        }
        id<TTIOStorageDataset> nameDs = [sc createDatasetNamed:@"read_names"
                                                     precision:TTIOPrecisionUInt8
                                                        length:encoded.length
                                                     chunkSize:65536
                                                   compression:TTIOCompressionNone
                                              compressionLevel:0
                                                         error:error];
        if (!nameDs) return NO;
        if (![nameDs writeAll:encoded error:error]) return NO;
        if (![nameDs setAttributeValue:@((uint8_t)TTIOCompressionNameTokenizedV2)
                               forName:@"compression"
                                 error:error]) return NO;
    } else {
        // native lib unavailable and run is
        // non-empty — no v1 NAME_TOKENIZED fallback exists.
        [NSException raise:NSInternalInconsistencyException
                    format:@"NAME_TOKENIZED_V2 codec requires the "
                           @"native libttio_rans library to be linked. "
                           @"Build with build.sh and ensure "
                           @"libttio_rans.so/dylib is present in "
                           @"$TTIO_NATIVE_LIB_DIR."];
    }

    // mate_info always emitted as the inline_v2
    // codec (id 13) when libttio_rans is linked. The Phase F per-
    // field subgroup writer and the M82 compound fallback are gone.
    // Empty-run short-circuit OMITS the mate_info group entirely
    // (cross-language convention shared with Python and Java; readers
    // treat absence as "no mate info").
    if (_bulkObjC != nil && _bulkObjC.mateInfoBlob != nil) {
        // Phase 2c-T: skip codec encode.
        if (!_TTIO_PhaseT_WriteMateInfoBulkStorage(sc,
                _bulkObjC.mateInfoBlob,
                _bulkObjC.mateInfoChromNames ?: @[], error)) return NO;
    } else if (run.readCount == 0) {
        // Omit the mate_info group — no children to write.
    } else if ([TTIOMateInfoV2 nativeAvailable]) {
        if (!_TTIO_V17_WriteMateInfoInlineV2Storage(sc, run, ctx.chromNameToId, error)) return NO;
    } else {
        [NSException raise:NSInternalInconsistencyException
                    format:@"mate_info inline_v2 codec requires the "
                           @"native libttio_rans library to be linked. "
                           @"Build with build.sh and ensure "
                           @"libttio_rans.so/dylib is present in "
                           @"$TTIO_NATIVE_LIB_DIR."];
    }

    // M101 — SAM optional fields. A run whose reads carry no tags
    // writes no tags channel, so tag-less files are unchanged.
    if ([run hasTags] && !_TTIO_M101_WriteTagsStorage(sc, run, error)) return NO;

    return YES;
}

// Task 30: provider-agnostic write of one /study/ms_runs/<name>/
// subtree via the StorageGroup protocol. Mirrors the per-MS-run section
// of -writeMinimalToPath: (HDF5 fast path, lines 2693-2815), but using
// only protocol-level primitives so memory:// / sqlite:// / zarr:// URLs
// work for write as they already do for read.
//
// Byte-exact parity with the HDF5 fast path is unattainable for
// non-HDF5 providers (each backend has its own physical layout), but
// the *logical* layout (groups, datasets, attribute names, attribute
// values, dataset precisions, dataset lengths) is identical so the
// existing read path (TTIOAcquisitionRun readFromStorageGroup:name:)
// reconstructs an equivalent in-memory dataset.
+ (BOOL)writeMSRunStorage:(TTIOWrittenRun *)run
                  toGroup:(id<TTIOStorageGroup>)parent
                     name:(NSString *)name
              compression:(TTIOCompression)compression
                    error:(NSError **)error
{
    id<TTIOStorageGroup> runGroup = [parent createGroupNamed:name error:error];
    if (!runGroup) return NO;

    NSUInteger spectrumCount = run.offsets.length / sizeof(int64_t);
    if (![runGroup setAttributeValue:@(run.acquisitionMode)
                              forName:@"acquisition_mode" error:error]) return NO;
    if (![runGroup setAttributeValue:@((int64_t)spectrumCount)
                              forName:@"spectrum_count" error:error]) return NO;
    if (![runGroup setAttributeValue:run.spectrumClassName ?: @""
                              forName:@"spectrum_class" error:error]) return NO;
    if (run.nucleusType.length > 0) {
        if (![runGroup setAttributeValue:run.nucleusType
                                  forName:@"nucleus_type" error:error]) return NO;
    }

    // Per-run provenance: write the JSON mirror so the storage-protocol
    // read path (which only consumes @provenance_json) finds the records.
    // The compound-dataset emission is intentionally HDF5-only — the
    // read-side fallback handles the missing /steps gracefully.
    if (run.provenanceRecords.count > 0) {
        NSMutableArray *plists =
            [NSMutableArray arrayWithCapacity:run.provenanceRecords.count];
        for (TTIOProvenanceRecord *r in run.provenanceRecords) {
            [plists addObject:[r asPlist]];
        }
        NSError *jErr = nil;
        NSData *json =
            [NSJSONSerialization dataWithJSONObject:plists
                                              options:0
                                                error:&jErr];
        if (!json) {
            if (error) *error = jErr;
            return NO;
        }
        NSString *jstr =
            [[NSString alloc] initWithData:json encoding:NSUTF8StringEncoding];
        if (![runGroup setAttributeValue:jstr
                                  forName:@"provenance_json" error:error]) return NO;
    }

    // Empty instrument_config skeleton for parity with HDF5 writeMinimal.
    id<TTIOStorageGroup> cfg =
        [runGroup createGroupNamed:@"instrument_config" error:error];
    if (!cfg) return NO;
    for (NSString *fieldName in @[@"manufacturer", @"model", @"serial_number",
                                   @"source_type", @"analyzer_type",
                                   @"detector_type"]) {
        if (![cfg setAttributeValue:@""
                            forName:fieldName error:error]) return NO;
    }

    // spectrum_index — same layout as TTIOSpectrumIndex -writeToGroup:.
    id<TTIOStorageGroup> idxG = [runGroup createGroupNamed:@"spectrum_index" error:error];
    if (!idxG) return NO;
    if (![idxG setAttributeValue:@((int64_t)spectrumCount)
                          forName:@"count" error:error]) return NO;
    // offsets is omitted on disk; readers derive it from
    // cumsum(lengths).
    if (!writeIndexArrayStorage(idxG, @"lengths",
                                 TTIOPrecisionUInt32, run.lengths,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"retention_times",
                                 TTIOPrecisionFloat64, run.retentionTimes,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"ms_levels",
                                 TTIOPrecisionInt32, run.msLevels,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"polarities",
                                 TTIOPrecisionInt32, run.polarities,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"precursor_mzs",
                                 TTIOPrecisionFloat64, run.precursorMzs,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"precursor_charges",
                                 TTIOPrecisionInt32, run.precursorCharges,
                                 compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"base_peak_intensities",
                                 TTIOPrecisionFloat64, run.basePeakIntensities,
                                 compression, error)) return NO;
    // M102 pixel columns (nil-skipped when the run has none).
    if (!writeIndexArrayStorage(idxG, @"pixel_x", TTIOPrecisionInt32,
                                 run.pixelX, compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"pixel_y", TTIOPrecisionInt32,
                                 run.pixelY, compression, error)) return NO;
    if (!writeIndexArrayStorage(idxG, @"pixel_z", TTIOPrecisionInt32,
                                 run.pixelZ, compression, error)) return NO;

    // signal_channels — pre-flattened NSData buffers, written straight
    // through. channel_names attribute is the comma-joined ordered list
    // (matches the HDF5 fast path's allKeys ordering for the parity
    // case where channels are mz + intensity).
    id<TTIOStorageGroup> channels =
        [runGroup createGroupNamed:@"signal_channels" error:error];
    if (!channels) return NO;
    NSArray *channelNames = run.channelData.allKeys;
    NSString *namesJoined = [channelNames componentsJoinedByString:@","];
    if (![channels setAttributeValue:namesJoined
                              forName:@"channel_names" error:error]) return NO;

    for (NSString *chName in channelNames) {
        NSData *buf = run.channelData[chName];
        NSUInteger total = buf.length / sizeof(double);
        NSString *dsName = [chName stringByAppendingString:@"_values"];
        id<TTIOStorageDataset> ds =
            [channels createDatasetNamed:dsName
                               precision:TTIOPrecisionFloat64
                                  length:total
                               chunkSize:65536
                             compression:compression
                        compressionLevel:6
                                   error:error];
        if (!ds) return NO;
        if (![ds writeAll:buf error:error]) return NO;
    }

    return YES;
}

// provider-agnostic minimal write — supports memory:// /
// sqlite:// / zarr:// URLs. Task 30 : MS runs are now supported
// via the StorageGroup protocol path; the HDF5 fast path in
// -writeMinimalToPath: still handles plain filesystem paths for
// byte-exact parity.
+ (BOOL)writeMinimalGenomicViaProviderURL:(NSString *)url
                                       title:(NSString *)title
                          isaInvestigationId:(NSString *)isaId
                                  msRuns:(NSDictionary *)msRuns
                                 genomicRuns:(NSDictionary *)genomicRuns
                                       error:(NSError **)error
{
    id<TTIOStorageProvider> prov =
        [[TTIOProviderRegistry sharedRegistry] openURL:url
                                                  mode:TTIOStorageOpenModeCreate
                                              provider:nil
                                                 error:error];
    if (!prov) return NO;

    @try {
        id<TTIOStorageGroup> root = [prov rootGroupWithError:error];
        if (!root) return NO;

        // Feature flags (mirror what TTIOFeatureFlags writeFormatVersion
        // does for HDF5: ttio_format_version + ttio_features attrs on
        // the root, JSON-encoded array for features). Memory provider
        // accepts NSString attribute values directly.
        NSMutableArray *features = [@[
            [TTIOFeatureFlags featureBaseV1],
            [TTIOFeatureFlags featureCompoundIdentifications],
            [TTIOFeatureFlags featureCompoundQuantifications],
            [TTIOFeatureFlags featureCompoundProvenance],
            [TTIOFeatureFlags featureCompoundPerRunProvenance],
            [TTIOFeatureFlags featureCompoundHeaders],
            [TTIOFeatureFlags featureNative2DNMR],
            [TTIOFeatureFlags featureNativeMSImageCube],
        ] mutableCopy];
        BOOL hasGenomic = genomicRuns.count > 0;
        if (hasGenomic) {
            if (![features containsObject:[TTIOFeatureFlags featureOptGenomic]]) {
                [features addObject:[TTIOFeatureFlags featureOptGenomic]];
            }
        }
        // M101: opt_sam_tags when any genomic run carries SAM optional
        // fields (a tags channel); tag-less files are unchanged.
        for (TTIOWrittenGenomicRun *gr in [genomicRuns objectEnumerator]) {
            if ([gr hasTags]) {
                if (![features containsObject:[TTIOFeatureFlags featureOptSamTags]]) {
                    [features addObject:[TTIOFeatureFlags featureOptSamTags]];
                }
                break;
            }
        }
        // M102: opt_pixel_coordinates when any MS run carries the
        // spectrum_index pixel columns.
        BOOL anyPixels = NO;
        for (NSString *rn in msRuns) {
            if (!validateWrittenRunPixels(msRuns[rn], rn, &anyPixels, error)) return NO;
        }
        if (anyPixels) {
            [features addObject:[TTIOFeatureFlags featurePixelCoordinates]];
        }
        if (![root setAttributeValue:kTTIOFormatVersion
                              forName:@"ttio_format_version" error:error]) return NO;
        NSData *featJSON = [NSJSONSerialization dataWithJSONObject:features options:0 error:NULL];
        NSString *featStr = [[NSString alloc] initWithData:featJSON encoding:NSUTF8StringEncoding];
        if (![root setAttributeValue:featStr
                              forName:@"ttio_features" error:error]) return NO;

        id<TTIOStorageGroup> study = [root createGroupNamed:@"study" error:error];
        if (!study) return NO;
        if (![study setAttributeValue:title ?: @""
                               forName:@"title" error:error]) return NO;
        if (![study setAttributeValue:isaId ?: @""
                               forName:@"isa_investigation_id" error:error]) return NO;

        // ms_runs subtree — Task 30 wires MS runs through the storage
        // protocol so memory/sqlite/zarr URLs work for write. Compression
        // is downgraded to None for backends that reject it (Zarr).
        TTIOCompression cx = task30CompressionForProvider(prov);
        id<TTIOStorageGroup> msG = [study createGroupNamed:@"ms_runs" error:error];
        if (!msG) return NO;
        NSArray *msNames = [[msRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];
        if (![msG setAttributeValue:[msNames componentsJoinedByString:@","]
                            forName:@"_run_names" error:error]) return NO;
        for (NSString *runName in msNames) {
            TTIOWrittenRun *run = msRuns[runName];
            if (![self writeMSRunStorage:run
                                  toGroup:msG
                                     name:runName
                              compression:cx
                                    error:error]) return NO;
        }
        // Empty nmr_runs for parity (readers expect the group).
        id<TTIOStorageGroup> nmrG = [study createGroupNamed:@"nmr_runs" error:error];
        if (!nmrG) return NO;
        if (![nmrG setAttributeValue:@""
                             forName:@"_run_names" error:error]) return NO;

        if (hasGenomic) {
            id<TTIOStorageGroup> gG = [study createGroupNamed:@"genomic_runs" error:error];
            if (!gG) return NO;
            NSArray *gNames = [[genomicRuns allKeys]
                sortedArrayUsingSelector:@selector(compare:)];
            if (![gG setAttributeValue:[gNames componentsJoinedByString:@","]
                               forName:@"_run_names" error:error]) return NO;
            for (NSString *gName in gNames) {
                TTIOWrittenGenomicRun *gRun = genomicRuns[gName];
                if (gRun.optLegacyWholeChannel) {
                    if (![self writeGenomicRunStorage:gRun
                                               toGroup:gG
                                                  name:gName
                                                 error:error]) return NO;
                } else {
                    if (![self _ttio_streamGenomicRun:gRun name:gName study:study error:error]) return NO;
                }
            }
        }
    }
    @finally {
        [prov close];
    }
    return YES;
}

/* Route one run through TTIOGenomicStreamWriter (blocks_v1) into the
 * /study group. `study` must already hold genomic_runs with the
 * run's name in @_run_names; the writer leaves both alone. */
+ (BOOL)_ttio_streamGenomicRun:(TTIOWrittenGenomicRun *)run
                          name:(NSString *)name
                         study:(id<TTIOStorageGroup>)study
                         error:(NSError **)error
{
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:study
                                                                             runName:name
                                                                             options:o];
    if (![w appendBatch:run error:error]) return NO;
    return [w close:error];
}

// write one /study/genomic_runs/<name>/ subtree. Mirrors the
// per-MS-run writer but for the genomic data model. Uses TTIOGenomicIndex
// for the index subgroup + TTIOCompoundIO for the 3 VL compound
// datasets (cigars, read_names, mate_info) under signal_channels/.
+ (BOOL)writeGenomicRun:(TTIOWrittenGenomicRun *)run
                 toGroup:(TTIOHDF5Group *)parent
                    name:(NSString *)name
                   error:(NSError **)error
{
    // validate signal-channel codec overrides before any HDF5
    // mutation. Raises NSInvalidArgumentException on programmer
    // error; the file is left untouched.
    _TTIO_M86_ValidateOverrides(run.signalCodecOverrides);
    // reject per-field mate_info_* overrides when inline_v2 active.
    _TTIO_V17_ValidateMateInfoV2Overrides(run);
    _TTIO_M97_ValidateQualityBinnedPlatform(run);

    TTIOHDF5Group *rg = [parent createGroupNamed:name error:error];
    if (!rg) return NO;

    // Run-level attributes.
    if (![rg setIntegerAttribute:@"acquisition_mode"
                            value:run.acquisitionMode error:error]) return NO;
    if (![rg setStringAttribute:@"modality"
                           value:@"genomic_sequencing" error:error]) return NO;
    if (![rg setIntegerAttribute:@"spectrum_class" value:5 error:error]) return NO;
    if (![rg setStringAttribute:@"reference_uri"
                           value:run.referenceUri error:error]) return NO;
    if (![rg setStringAttribute:@"platform"
                           value:run.platform error:error]) return NO;
    if (run.readRole.length > 0
        && ![rg setStringAttribute:@"read_role"
                             value:run.readRole error:error]) return NO;
    if (![rg setStringAttribute:@"sample_name"
                           value:run.sampleName error:error]) return NO;
    if (![rg setIntegerAttribute:@"read_count"
                            value:(int64_t)run.readCount error:error]) return NO;

    // genomic_index subgroup (delegates to TTIOGenomicIndex).
    TTIOGenomicIndex *idx = [[TTIOGenomicIndex alloc]
        initWithOffsets:run.offsetsData
                lengths:run.lengthsData
            chromosomes:run.chromosomes
              positions:run.positionsData
       mappingQualities:run.mappingQualitiesData
                  flags:run.flagsData];
    TTIOHDF5Group *idxG = [rg createGroupNamed:@"genomic_index" error:error];
    if (!idxG) return NO;
    // GenomicIndex.writeToGroup takes id<TTIOStorageGroup>; wrap via the
    // HDF5 provider's adapter. Use the same trick: TTIOHDF5GroupAdapter
    // is created by openProviderURL but we can construct it via
    // TTIOHDF5Provider's escape hatch. Simpler: pass the raw HDF5Group
    // via `id`-cast since GenomicIndex.writeToGroup checks
    // respondsToSelector:@selector(unwrap) and falls through to direct
    // TTIOCompoundIO + the storage protocol calls also work on
    // TTIOHDF5Group via category methods. Easiest: just reuse the
    // helper directly — TTIOGenomicIndex's internal writes use
    // createDatasetNamed:precision:length:chunkSize:compression:
    // compressionLevel:error: which TTIOHDF5Group implements with a
    // slightly different signature (no `compression` arg). Wrap via
    // adapter to bridge.
    id<TTIOStorageGroup> idxGAdapter =
        (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(idxG);
    if (!idxGAdapter) return NO;
    if (![idx writeToGroup:idxGAdapter error:error]) return NO;

    // signal_channels subgroup.
    TTIOHDF5Group *sc = [rg createGroupNamed:@"signal_channels" error:error];
    if (!sc) return NO;

    // 5 typed channels (use TTIOGenomicIndex's static writeTypedChannel
    // helper isn't accessible — inline the same pattern via the
    // adapter). These match the precision choices in the spec:
    // positions=int64, sequences=uint8, qualities=uint8, flags=uint32,
    // mapping_qualities=uint8.
    //
    // sequences and qualities go through the byte-channel codec
    // dispatcher so an override (rANS / BASE_PACK) is honoured.
    // positions / flags / mapping_qualities are NOT written
    // under signal_channels/. They live exclusively under
    // genomic_index/, mirroring MS's spectrum_index/ pattern. See
    // docs/format-spec.md §4 and §10.7.
    TTIOCompression codec = run.signalCompression;
    // sequences dispatch — prefer refdiff_v2
    // group layout when eligible (no override + native lib +
    // reference + all reads mapped); otherwise fall through to the
    // M86 byte-channel writer with whatever explicit override the
    // caller supplied. The v1 REF_DIFF override (codec id 9) is no
    // longer accepted (override-validation rejects it).
    // Phase 2c-T: when bulk_v2_blobs.refDiffBlob is set, write the
    // verbatim wire blob and skip the codec encode entirely.
    TTIOBulkV2Blobs *_bulkObjC = (TTIOBulkV2Blobs *)run.bulkV2Blobs;
    NSNumber *seqOverride = run.signalCodecOverrides[@"sequences"];
    if (_bulkObjC != nil && _bulkObjC.refDiffBlob != nil) {
        if (![_bulkObjC.refDiffReferenceUri isEqualToString:run.referenceUri]) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2110
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"BulkV2Blobs.refDiffReferenceUri does not match run.referenceUri"}];
            return NO;
        }
        if (!_TTIO_PhaseT_WriteRefDiffV2BulkHDF5(sc, _bulkObjC.refDiffBlob, error)) return NO;
    } else if (seqOverride == nil && _TTIO_V18_UseRefDiffV2(run)) {
        // group layout with refdiff_v2 child @compression=14.
        if (!_TTIO_V18_WriteRefDiffV2SequencesHDF5(sc, run, error)) return NO;
    } else if (seqOverride != nil
               && (TTIOCompression)[seqOverride unsignedIntegerValue] == TTIOCompressionSeqCm) {
        if (!_TTIO_M103_WriteSequencesSeqCm(sc, run, error)) return NO;
    } else {
        if (!_TTIO_M86_WriteByteChannel(sc, @"sequences", run.sequencesData,
                                        codec,
                                        seqOverride,
                                        error)) return NO;
    }
    // qualities (uint8) — codec-aware. M94.Z v1.2: when the override
    // (or v1.5 auto-default) selects FQZCOMP_NX16_Z, dispatch to the
    // context-aware encoder; otherwise fall through to the M86 byte-
    // channel writer.
    NSNumber *qualOverride = run.signalCodecOverrides[@"qualities"];
    if (qualOverride == nil) {
        qualOverride = _TTIO_M94_DefaultQualitiesCodec(run);
    }
    if (qualOverride != nil &&
        (TTIOCompression)[qualOverride unsignedIntegerValue]
            == TTIOCompressionFqzcompNx16Z) {
        if (!_TTIO_M94Z_WriteQualitiesFqzcompNx16Z(sc, run, error)) return NO;
    } else {
        if (!_TTIO_M86_WriteByteChannel(sc, @"qualities", run.qualitiesData,
                                        codec,
                                        qualOverride,
                                        error)) return NO;
    }

    // 3 compound datasets via TTIOCompoundIO (HDF5-direct).
    NSArray *vlValueField = @[
        [TTIOCompoundField fieldWithName:@"value" kind:TTIOCompoundFieldKindVLString]
    ];

    // schema lift for cigars. When an override is set,
    // replace the M82 compound dataset with a flat 1-D uint8 dataset
    // of the same name carrying the codec output, plus an
    // @compression attribute naming the codec id (Binding Decisions
    // §120-§122). Three codec choices are supported (rANS uses a
    // length-prefix-concat byte stream over the CIGAR list — Gotcha
    // §139 — while NAME_TOKENIZED consumes the list[str] directly).
    // No HDF5 filter applied (Binding Decision §87).
    NSNumber *cigarsOverride = run.signalCodecOverrides[@"cigars"];
    if (cigarsOverride != nil) {
        TTIOCompression cigarsCodec =
            (TTIOCompression)[cigarsOverride unsignedIntegerValue];
        NSData *encoded = _TTIO_M86_EncodeCigarsWithCodec(run.cigars,
                                                          cigarsCodec);
        if (!encoded) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2060
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"M86 Phase C: cigars codec %lu encode "
                                @"returned nil",
                                (unsigned long)cigarsCodec]}];
            return NO;
        }
        TTIOHDF5Dataset *cigarDs = [sc createDatasetNamed:@"cigars"
                                                precision:TTIOPrecisionUInt8
                                                   length:encoded.length
                                                chunkSize:65536
                                              compression:TTIOCompressionNone
                                         compressionLevel:0
                                                    error:error];
        if (!cigarDs) return NO;
        if (![cigarDs writeData:encoded error:error]) return NO;
        if (!_TTIO_M86_WriteUInt8Attribute([cigarDs datasetId], "compression",
                                           (uint8_t)cigarsCodec, error)) return NO;
    } else {
        NSMutableArray *cigarRows = [NSMutableArray arrayWithCapacity:run.cigars.count];
        for (NSString *c in run.cigars) [cigarRows addObject:@{@"value": c}];
        if (![TTIOCompoundIO writeGeneric:cigarRows
                                  intoGroup:sc datasetNamed:@"cigars"
                                      fields:vlValueField error:error]) return NO;
    }

    // read_names always written via
    // NAME_TOKENIZED_V2 (codec id 15). The v1 NAME_TOKENIZED override
    // (id 8) and the M82 compound fallback are gone. Empty-run
    // short-circuit writes a zero-length uint8 dataset with
    // @compression=15.
    if (_bulkObjC != nil && _bulkObjC.nameTokBlob != nil) {
        // Phase 2c-T: skip codec encode.
        if (!_TTIO_PhaseT_WriteReadNamesBulkHDF5(sc, _bulkObjC.nameTokBlob, error)) return NO;
    } else if (run.readCount == 0) {
        TTIOHDF5Dataset *nameDs = [sc createDatasetNamed:@"read_names"
                                               precision:TTIOPrecisionUInt8
                                                  length:0
                                               chunkSize:1
                                             compression:TTIOCompressionNone
                                        compressionLevel:0
                                                   error:error];
        if (!nameDs) return NO;
        if (!_TTIO_M86_WriteUInt8Attribute([nameDs datasetId], "compression",
                                           (uint8_t)TTIOCompressionNameTokenizedV2,
                                           error)) return NO;
    } else if ([TTIONameTokenizerV2 nativeAvailable]) {
        NSData *encoded = _TTIO_M86_EncodeNamesViaRegistry(run.readNames);
        if (!encoded) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:2032
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"v1.0 NAME_TOKENIZED_V2 encode of "
                           @"read_names returned nil"}];
            return NO;
        }
        TTIOHDF5Dataset *nameDs = [sc createDatasetNamed:@"read_names"
                                               precision:TTIOPrecisionUInt8
                                                  length:encoded.length
                                               chunkSize:65536
                                             compression:TTIOCompressionNone
                                        compressionLevel:0
                                                   error:error];
        if (!nameDs) return NO;
        if (![nameDs writeData:encoded error:error]) return NO;
        if (!_TTIO_M86_WriteUInt8Attribute([nameDs datasetId], "compression",
                                           (uint8_t)TTIOCompressionNameTokenizedV2,
                                           error)) return NO;
    } else {
        // native lib unavailable and run is
        // non-empty — no v1 NAME_TOKENIZED fallback exists.
        [NSException raise:NSInternalInconsistencyException
                    format:@"NAME_TOKENIZED_V2 codec requires the "
                           @"native libttio_rans library to be linked. "
                           @"Build with build.sh and ensure "
                           @"libttio_rans.so/dylib is present in "
                           @"$TTIO_NATIVE_LIB_DIR."];
    }

    // mate_info always emitted as the inline_v2
    // codec (id 13). The Phase F per-field subgroup writer and the
    // M82 compound fallback are gone. Empty-run short-circuit OMITS
    // the mate_info group entirely (cross-language convention shared
    // with Python and Java; readers treat absence as "no mate info").
    if (_bulkObjC != nil && _bulkObjC.mateInfoBlob != nil) {
        // Phase 2c-T: skip codec encode.
        if (!_TTIO_PhaseT_WriteMateInfoBulkHDF5(sc,
                _bulkObjC.mateInfoBlob,
                _bulkObjC.mateInfoChromNames ?: @[], error)) return NO;
    } else if (run.readCount == 0) {
        // Omit the mate_info group — no children to write.
    } else if ([TTIOMateInfoV2 nativeAvailable]) {
        if (!_TTIO_V17_WriteMateInfoInlineV2HDF5(sc, run, error)) return NO;
    } else {
        [NSException raise:NSInternalInconsistencyException
                    format:@"mate_info inline_v2 codec requires the "
                           @"native libttio_rans library to be linked. "
                           @"Build with build.sh and ensure "
                           @"libttio_rans.so/dylib is present in "
                           @"$TTIO_NATIVE_LIB_DIR."];
    }

    // M101 — SAM optional fields. A run whose reads carry no tags
    // writes no tags channel, so tag-less files are unchanged.
    if ([run hasTags] && !_TTIO_M101_WriteTagsHDF5(sc, run, error)) return NO;

    // Phase 1: per-run provenance compound at <run>/provenance/steps,
    // mirroring the AcquisitionRun MS path. Absent when the
    // WrittenGenomicRun carries no records — preserving pre-Phase-1
    // byte parity for callers that don't ship provenance.
    if (run.provenanceRecords.count > 0) {
        TTIOHDF5Group *provGroup =
            [rg createGroupNamed:@"provenance" error:error];
        if (!provGroup) return NO;
        if (![TTIOCompoundIO writeProvenance:run.provenanceRecords
                                   intoGroup:provGroup
                                datasetNamed:@"steps"
                                       error:error]) return NO;
    }

    return YES;
}

+ (BOOL)writeMinimalToPath:(NSString *)path
                      title:(NSString *)title
        isaInvestigationId:(NSString *)isaId
                    msRuns:(NSDictionary<NSString *, TTIOWrittenRun *> *)runs
            identifications:(NSArray *)identifications
            quantifications:(NSArray *)quantifications
          provenanceRecords:(NSArray *)provenance
                      error:(NSError **)error
{
    return [self writeMinimalToPath:path
                              title:title
                 isaInvestigationId:isaId
                             msRuns:runs
                         genomicRuns:nil
                     identifications:identifications
                     quantifications:quantifications
                  provenanceRecords:provenance
                              error:error];
}

+ (BOOL)writeMinimalToPath:(NSString *)path
                      title:(NSString *)title
        isaInvestigationId:(NSString *)isaId
                  mixedRuns:(NSDictionary<NSString *, id> *)mixedRuns
                genomicRuns:(NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)genomicRuns
            identifications:(NSArray *)identifications
            quantifications:(NSArray *)quantifications
          provenanceRecords:(NSArray *)provenance
                      error:(NSError **)error
{
    return [self writeMinimalToPath:path
                              title:title
                 isaInvestigationId:isaId
                          mixedRuns:mixedRuns
                        genomicRuns:genomicRuns
                    identifications:identifications
                    quantifications:quantifications
                  provenanceRecords:provenance
                           progress:nil
                              error:error];
}

+ (BOOL)writeMinimalToPath:(NSString *)path
                      title:(NSString *)title
        isaInvestigationId:(NSString *)isaId
                  mixedRuns:(NSDictionary<NSString *, id> *)mixedRuns
                genomicRuns:(NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)genomicRuns
            identifications:(NSArray *)identifications
            quantifications:(NSArray *)quantifications
          provenanceRecords:(NSArray *)provenance
                   progress:(TTIOProgressBlock)progress
                      error:(NSError **)error
{
    // Phase 2: split the mixed dict into MS-only + genomic-only maps,
    // dispatching per-value on isKindOfClass:. Pre-existing
    // genomicRuns= entries are merged in; a name appearing in BOTH
    // raises NSError with code 1100 (matches Python's ValueError on
    // collision).
    NSMutableDictionary<NSString *, TTIOWrittenRun *> *splitMS =
        [NSMutableDictionary dictionary];
    NSMutableDictionary<NSString *, TTIOWrittenGenomicRun *> *splitG =
        [NSMutableDictionary dictionaryWithDictionary:(genomicRuns ?: @{})];

    for (NSString *name in mixedRuns) {
        id value = mixedRuns[name];
        if ([value isKindOfClass:[TTIOWrittenGenomicRun class]]) {
            if (splitG[name] != nil) {
                if (error) *error = [NSError
                    errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:1100
                           userInfo:@{NSLocalizedDescriptionKey:
                               [NSString stringWithFormat:
                                    @"Phase 2 mixed runs dict: name '%@' "
                                    @"appears in both mixedRuns and "
                                    @"genomicRuns", name]}];
                return NO;
            }
            splitG[name] = (TTIOWrittenGenomicRun *)value;
        } else if ([value isKindOfClass:[TTIOWrittenRun class]]) {
            splitMS[name] = (TTIOWrittenRun *)value;
        } else {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:1101
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"Phase 2 mixed runs dict: value for '%@' "
                                @"is %@; expected TTIOWrittenRun or "
                                @"TTIOWrittenGenomicRun",
                                name, NSStringFromClass([value class])]}];
            return NO;
        }
    }

    // ──────────────────────────────────────────────────────────────────
    // Progress: compute section presence flags in §5.4 order, fire a
    // baseline (0, total), then dispatch to the inner write. After a
    // successful write, fire one (idx, total) per present section in
    // §5.4 order. Mirrors Python's _section_flags table; encryption is
    // a placeholder (writeMinimal never encrypts) so it never fires.
    // Empty sections are skipped from the count.
    // ──────────────────────────────────────────────────────────────────
    BOOL hasProvenance      = provenance.count > 0;
    BOOL hasSubjects        = NO;   // writeMinimal does not accept subjects.
    BOOL hasSamples         = NO;   // writeMinimal does not accept samples.
    BOOL hasReferences      = (splitG.count > 0);
    BOOL hasImage           = NO;   // writeMinimal does not accept image cubes.
    BOOL hasIdentifications = identifications.count > 0;
    BOOL hasQuantifications = quantifications.count > 0;
    BOOL hasRuns            = (splitMS.count > 0 || splitG.count > 0);

    NSUInteger progressTotal =
        (hasProvenance       ? 1 : 0) +
        (hasSubjects         ? 1 : 0) +
        (hasSamples          ? 1 : 0) +
        (hasReferences       ? 1 : 0) +
        (hasImage            ? 1 : 0) +
        (hasIdentifications  ? 1 : 0) +
        (hasQuantifications  ? 1 : 0) +
        (hasRuns             ? 1 : 0);

    if (progress) progress((int64_t)0, (int64_t)progressTotal);

    BOOL ok = [self writeMinimalToPath:path
                                 title:title
                    isaInvestigationId:isaId
                                msRuns:splitMS
                           genomicRuns:splitG
                       identifications:identifications
                       quantifications:quantifications
                     provenanceRecords:provenance
                                 error:error];
    if (!ok) return NO;

    // Fire section-by-section progress in §5.4 order. Because the
    // inner write is synchronous + monolithic, these fires arrive
    // clustered at end-of-write rather than interleaved. The
    // contract — "one (idx, total) per non-empty section in §5.4
    // order" — is satisfied; UI consumers see correct cumulative
    // totals even if they cluster. Threading the fires through the
    // 200+-line inner write is a follow-up if interleaved cadence
    // becomes a hard requirement.
    if (progress) {
        NSUInteger done = 0;
        if (hasProvenance)      progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasSubjects)        progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasSamples)         progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasReferences)      progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasImage)           progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasIdentifications) progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasQuantifications) progress((int64_t)(++done), (int64_t)progressTotal);
        if (hasRuns)            progress((int64_t)(++done), (int64_t)progressTotal);
    }
    return YES;
}

+ (BOOL)writeMinimalToPath:(NSString *)path
                      title:(NSString *)title
        isaInvestigationId:(NSString *)isaId
                    msRuns:(NSDictionary<NSString *, TTIOWrittenRun *> *)runs
                genomicRuns:(NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)genomicRuns
            identifications:(NSArray *)identifications
            quantifications:(NSArray *)quantifications
          provenanceRecords:(NSArray *)provenance
                      error:(NSError **)error
{
    return [self writeMinimalToPath:path
                              title:title
                 isaInvestigationId:isaId
                             msRuns:runs
                        genomicRuns:genomicRuns
                     assemblyGraphs:nil
                    identifications:identifications
                    quantifications:quantifications
                  provenanceRecords:provenance
                              error:error];
}

+ (BOOL)writeMinimalToPath:(NSString *)path
                      title:(NSString *)title
        isaInvestigationId:(NSString *)isaId
                    msRuns:(NSDictionary<NSString *, TTIOWrittenRun *> *)runs
                genomicRuns:(NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)genomicRuns
             assemblyGraphs:(NSDictionary<NSString *, TTIOWrittenAssemblyGraph *> *)assemblyGraphs
            identifications:(NSArray *)identifications
            quantifications:(NSArray *)quantifications
          provenanceRecords:(NSArray *)provenance
                      error:(NSError **)error
{
    // M82.2 + Task 30: provider-agnostic write path for non-HDF5 URLs.
    // MS runs are now wired through the StorageGroup protocol (Task 30);
    // genomic_runs were already supported via M82.2. Identifications /
    // quantifications / provenance still go HDF5-only — they require
    // compound-dataset writes through TTIOCompoundIO which is HDF5-direct
    // today (follow-up work; the JSON-mirror attribute mechanism is the
    // workaround per-run for provenance).
    if (isNonHdf5ProviderURL(path)) {
        if (assemblyGraphs.count > 0) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:1002
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"writeMinimal via provider URL does not yet "
                           @"support assembly graphs; use the HDF5 fast "
                           @"path, or writeAssemblyGraph:named:"
                           @"toStudyGroup: against the provider's study "
                           @"group directly."}];
            return NO;
        }
        if (identifications.count > 0 || quantifications.count > 0 ||
            provenance.count > 0) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOSpectralDatasetErrorDomain" code:1001
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"writeMinimal via provider URL does not yet "
                           @"support dataset-level identifications / "
                           @"quantifications / provenance (use the HDF5 "
                           @"fast path for those, or pass them as per-run "
                           @"provenance which is supported)."}];
            return NO;
        }
        return [self writeMinimalGenomicViaProviderURL:path
                                                  title:title
                                     isaInvestigationId:isaId
                                                 msRuns:runs
                                            genomicRuns:genomicRuns
                                                  error:error];
    }

    // M102: validate pixel columns before creating the file.
    BOOL anyPixels = NO;
    for (NSString *rn in runs) {
        if (!validateWrittenRunPixels(runs[rn], rn, &anyPixels, error)) return NO;
    }

    TTIOHDF5Provider *p = [[TTIOHDF5Provider alloc] init];
    if (![p openURL:path mode:TTIOStorageOpenModeCreate error:error]) return NO;
    TTIOHDF5File *f = (TTIOHDF5File *)[p nativeHandle];
    if (!f) return NO;
    TTIOHDF5Group *root = [f rootGroup];

    // Same feature-flag set as -writeToFilePath: so readers can't tell
    // the two paths apart.
    NSMutableArray *features = [@[
        [TTIOFeatureFlags featureBaseV1],
        [TTIOFeatureFlags featureCompoundIdentifications],
        [TTIOFeatureFlags featureCompoundQuantifications],
        [TTIOFeatureFlags featureCompoundProvenance],
        [TTIOFeatureFlags featureCompoundPerRunProvenance],
        [TTIOFeatureFlags featureCompoundHeaders],
        [TTIOFeatureFlags featureNative2DNMR],
        [TTIOFeatureFlags featureNativeMSImageCube],
    ] mutableCopy];

    // opt_genomic is the canonical advertisement of genomic
    // content. Add it whenever genomicRuns is non-empty, idempotent
    // if a future caller pre-populates it. Bump format_version to 1.4
    // (which implies 1.3 + 1.1 — readers gate features by flag, not
    // by version equality).
    //
    BOOL hasGenomic = genomicRuns.count > 0;
    if (hasGenomic) {
        if (![features containsObject:[TTIOFeatureFlags featureOptGenomic]]) {
            [features addObject:[TTIOFeatureFlags featureOptGenomic]];
        }
    }

    // M98: opt_assembly_graph advertises assembly-graph content the
    // same way, added only when a graph is present (the REF_DIFF
    // bump rule: absent graphs leave ms-only files byte-identical).
    BOOL hasGraphs = assemblyGraphs.count > 0;
    if (hasGraphs) {
        if (![features containsObject:[TTIOFeatureFlags featureOptAssemblyGraph]]) {
            [features addObject:[TTIOFeatureFlags featureOptAssemblyGraph]];
        }
    }
    // M101: opt_sam_tags when any genomic run carries SAM optional
    // fields (a tags channel); tag-less files are unchanged.
    for (TTIOWrittenGenomicRun *gr in [genomicRuns objectEnumerator]) {
        if ([gr hasTags]) {
            if (![features containsObject:[TTIOFeatureFlags featureOptSamTags]]) {
                [features addObject:[TTIOFeatureFlags featureOptSamTags]];
            }
            break;
        }
    }

    // M102: opt_pixel_coordinates when any MS run carries the
    // spectrum_index pixel columns; other files are unchanged.
    if (anyPixels) {
        [features addObject:[TTIOFeatureFlags featurePixelCoordinates]];
    }

    if (![TTIOFeatureFlags writeFormatVersion:kTTIOFormatVersion
                                      features:features
                                        toRoot:root
                                         error:error]) return NO;

    TTIOHDF5Group *study = [root createGroupNamed:@"study" error:error];
    if (!study) return NO;
    if (![study setStringAttribute:@"title" value:(title ?: @"") error:error]) return NO;
    if (![study setStringAttribute:@"isa_investigation_id"
                              value:(isaId ?: @"") error:error]) return NO;

    TTIOHDF5Group *msRunsGroup = [study createGroupNamed:@"ms_runs" error:error];
    if (!msRunsGroup) return NO;
    NSArray *msNames = [[runs allKeys] sortedArrayUsingSelector:@selector(compare:)];
    if (![msRunsGroup setStringAttribute:@"_run_names"
                                    value:[msNames componentsJoinedByString:@","]
                                    error:error]) return NO;

    for (NSString *runName in msNames) {
        TTIOWrittenRun *run = runs[runName];

        TTIOHDF5Group *runGroup = [msRunsGroup createGroupNamed:runName error:error];
        if (!runGroup) return NO;

        NSUInteger spectrumCount = run.offsets.length / sizeof(int64_t);
        if (![runGroup setIntegerAttribute:@"acquisition_mode"
                                     value:run.acquisitionMode error:error]) return NO;
        if (![runGroup setIntegerAttribute:@"spectrum_count"
                                     value:(int64_t)spectrumCount error:error]) return NO;
        if (![runGroup setStringAttribute:@"spectrum_class"
                                    value:run.spectrumClassName error:error]) return NO;
        if (run.nucleusType.length > 0) {
            if (![runGroup setStringAttribute:@"nucleus_type"
                                        value:run.nucleusType error:error]) return NO;
        }

        // Per-run provenance — same compound + JSON-mirror layout as
        // -[TTIOAcquisitionRun writeToGroup:name:error:]. Mirrors
        // Python's ``_write_run`` helper (spectral_dataset.py) which
        // emits ``<run>/provenance/steps`` plus a legacy
        // ``@provenance_json`` attribute. Absent when the
        // TTIOWrittenRun carries no records — preserves byte parity
        // with pre-v0.6 callers.
        if (run.provenanceRecords.count > 0) {
            TTIOHDF5Group *provGroup =
                [runGroup createGroupNamed:@"provenance" error:error];
            if (!provGroup) return NO;
            if (![TTIOCompoundIO writeProvenance:run.provenanceRecords
                                       intoGroup:provGroup
                                    datasetNamed:@"steps"
                                           error:error]) return NO;

            NSMutableArray *plists =
                [NSMutableArray arrayWithCapacity:run.provenanceRecords.count];
            for (TTIOProvenanceRecord *r in run.provenanceRecords) {
                [plists addObject:[r asPlist]];
            }
            NSError *jErr = nil;
            NSData *json =
                [NSJSONSerialization dataWithJSONObject:plists
                                                  options:0
                                                    error:&jErr];
            if (!json) {
                if (error) *error = jErr;
                return NO;
            }
            NSString *jstr =
                [[NSString alloc] initWithData:json
                                       encoding:NSUTF8StringEncoding];
            if (![runGroup setStringAttribute:@"provenance_json"
                                        value:jstr error:error]) return NO;
        }

        // instrument_config subgroup — writeMinimal callers don't ship
        // instrument metadata; emit the same empty-string skeleton that
        // Python's write_minimal does so readers don't distinguish
        // writer.
        TTIOHDF5Group *cfg =
            [runGroup createGroupNamed:@"instrument_config" error:error];
        if (!cfg) return NO;
        for (NSString *fieldName in @[@"manufacturer", @"model", @"serial_number",
                                       @"source_type", @"analyzer_type",
                                       @"detector_type"]) {
            if (![cfg setStringAttribute:fieldName value:@"" error:error]) return NO;
        }

        // spectrum_index — same layout as TTIOSpectrumIndex -writeToGroup:.
        TTIOHDF5Group *idxG = [runGroup createGroupNamed:@"spectrum_index" error:error];
        if (!idxG) return NO;
        if (![idxG setIntegerAttribute:@"count"
                                 value:(int64_t)spectrumCount error:error]) return NO;
        // offsets is omitted on disk; readers derive it
        // from cumsum(lengths).
        if (!writeIndexArrayDS(idxG, @"lengths",
                                TTIOPrecisionUInt32, run.lengths, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"retention_times",
                                TTIOPrecisionFloat64, run.retentionTimes, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"ms_levels",
                                TTIOPrecisionInt32, run.msLevels, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"polarities",
                                TTIOPrecisionInt32, run.polarities, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"precursor_mzs",
                                TTIOPrecisionFloat64, run.precursorMzs, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"precursor_charges",
                                TTIOPrecisionInt32, run.precursorCharges, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"base_peak_intensities",
                                TTIOPrecisionFloat64, run.basePeakIntensities, error)) return NO;
        // M102 pixel columns (nil-skipped when the run has none).
        if (!writeIndexArrayDS(idxG, @"pixel_x",
                                TTIOPrecisionInt32, run.pixelX, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"pixel_y",
                                TTIOPrecisionInt32, run.pixelY, error)) return NO;
        if (!writeIndexArrayDS(idxG, @"pixel_z",
                                TTIOPrecisionInt32, run.pixelZ, error)) return NO;

        // v1.1 writeMinimal intentionally SKIPS the "opt_compound_headers"
        // duplicate spectrum_index/headers compound dataset. That feature
        // (added by TTIOCompoundIO writeCompoundHeadersForIndex:) writes
        // the parallel index arrays again as a 56-byte-per-row compound,
        // uncompressed + unchunked — ~5.6 MB on 100 K spectra. The parallel
        // arrays are authoritative; the compound copy exists only for
        // h5dump readability. Python's write_minimal doesn't emit it, and
        // its absence is the single biggest file-size difference between
        // the ObjC and Python minimal paths. Callers that need h5dump-
        // friendly compound headers should use the object-mode writer.

        // signal_channels — pre-flattened NSData buffers, written
        // straight through with no per-spectrum concat.
        TTIOHDF5Group *channels =
            [runGroup createGroupNamed:@"signal_channels" error:error];
        if (!channels) return NO;
        NSArray *channelNames = run.channelData.allKeys;
        NSString *namesJoined = [channelNames componentsJoinedByString:@","];
        if (![channels setStringAttribute:@"channel_names"
                                    value:namesJoined error:error]) return NO;

        // Phase 2: "gzip" left at its default resolves to codec 17
        // on MS runs unless the caller opted out. Matches Python's
        // write_minimal and the object-mode -writeToGroup: path.
        BOOL useFloatDelta =
            [run.signalCompression isEqualToString:@"gzip"]
            && !run.optDisableFloatDelta
            && [run.spectrumClassName isEqualToString:@"TTIOMassSpectrum"];

        for (NSString *chName in channelNames) {
            NSData *buf = run.channelData[chName];
            NSUInteger total = buf.length / sizeof(double);
            NSString *dsName = [chName stringByAppendingString:@"_values"];
            if (useFloatDelta) {
                NSData *stream = [TTIOFloatDeltaZstd encodeFloat64:buf];
                if (!stream) {
                    if (error) *error = [NSError
                        errorWithDomain:@"TTIOSpectralDatasetErrorDomain"
                                   code:1102
                               userInfo:@{NSLocalizedDescriptionKey:
                            [NSString stringWithFormat:
                                @"FLOAT_DELTA_ZSTD encode failed for '%@'",
                                dsName]}];
                    return NO;
                }
                TTIOHDF5Dataset *ds =
                    [channels createDatasetNamed:dsName
                                       precision:TTIOPrecisionUInt8
                                          length:stream.length
                                       chunkSize:65536
                                     compression:TTIOCompressionNone
                                compressionLevel:0
                                           error:error];
                if (!ds) return NO;
                if (![ds writeData:stream error:error]) return NO;
                if (![ds setAttributeValue:@(TTIOCompressionFloatDeltaZstd)
                                   forName:@"compression"
                                     error:error]) return NO;
                continue;
            }
            TTIOHDF5Dataset *ds =
                [channels createDatasetNamed:dsName
                                   precision:TTIOPrecisionFloat64
                                      length:total
                                   chunkSize:65536
                                 compression:TTIOCompressionZlib
                            compressionLevel:6
                                       error:error];
            if (!ds) return NO;
            if (![ds writeData:buf error:error]) return NO;
        }
    }

    // Empty nmr_runs group for byte-parity with -writeToFilePath:.
    TTIOHDF5Group *nmrRunsGroup = [study createGroupNamed:@"nmr_runs" error:error];
    if (!nmrRunsGroup) return NO;
    if (![nmrRunsGroup setStringAttribute:@"_run_names" value:@"" error:error]) return NO;

    // genomic_runs subtree (only when non-empty — pre-M82 byte
    // parity for ms-only files).
    if (hasGenomic) {
        // M93 v1.2: embed each unique reference (by URI) once at
        // /study/references/<uri>/ BEFORE writing the genomic runs.
        // Required so the read-side resolver finds the embedded data.
        if (![self class] || !_TTIO_M93_EmbedReferences(study,
                                                        genomicRuns,
                                                        error)) return NO;

        TTIOHDF5Group *gRunsGroup = [study createGroupNamed:@"genomic_runs" error:error];
        if (!gRunsGroup) return NO;
        NSArray *gNames = [[genomicRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];
        if (![gRunsGroup setStringAttribute:@"_run_names"
                                      value:[gNames componentsJoinedByString:@","]
                                      error:error]) return NO;
        for (NSString *gName in gNames) {
            TTIOWrittenGenomicRun *gRun = genomicRuns[gName];
            if (gRun.optLegacyWholeChannel) {
                if (![self writeGenomicRun:gRun
                                    toGroup:gRunsGroup
                                       name:gName
                                      error:error]) return NO;
            } else {
                id<TTIOStorageGroup> studyAdapter =
                    (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(study);
                if (!studyAdapter) return NO;
                if (![self _ttio_streamGenomicRun:gRun name:gName study:studyAdapter error:error]) return NO;
            }
        }
    }

    // assembly_graphs subtree (M98, only when non-empty).
    if (hasGraphs) {
        id<TTIOStorageGroup> studyAdapter =
            (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(study);
        if (!studyAdapter) return NO;
        NSArray *agNames = [[assemblyGraphs allKeys]
            sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *agName in agNames) {
            if (![self writeAssemblyGraph:assemblyGraphs[agName]
                                    named:agName
                             toStudyGroup:studyAdapter
                                    error:error]) return NO;
        }
    }

    if (identifications.count > 0) {
        if (![TTIOCompoundIO writeIdentifications:identifications
                                         intoGroup:study
                                      datasetNamed:@"identifications"
                                             error:error]) return NO;
    }
    if (quantifications.count > 0) {
        if (![TTIOCompoundIO writeQuantifications:quantifications
                                         intoGroup:study
                                      datasetNamed:@"quantifications"
                                             error:error]) return NO;
    }
    if (provenance.count > 0) {
        if (![TTIOCompoundIO writeProvenance:provenance
                                    intoGroup:study
                                 datasetNamed:@"provenance"
                                        error:error]) return NO;
    }

    return [f close];
}

@end
