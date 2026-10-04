/*
 * TTIOPerAUFile.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOPerAUFile
 * Declared In:   Protection/TTIOPerAUFile.h
 *
 * File-level per-AU encryption orchestrator. Reads plaintext
 * channels + writes the segments compound layout via the storage
 * provider abstraction.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOPerAUFile.h"
#import "TTIOPerAUEncryption.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOCompoundField.h"
#import "Dataset/TTIOCompoundIO.h"
#import "Dataset/TTIOCompoundIO+Internal.h"
#import "Genomics/TTIOGenomicIndex.h"  // TTIOOffsetsFromLengths (v1.10 #10)
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "ValueClasses/TTIOEnums.h"
#import "Codecs/TTIOFloatDeltaZstd.h"  // codec id 17, Phase 2 MS default
#import "Codecs/Registry/TTIOCodecRegistry.h"  // M98 assembly sequences
#import "Genomics/TTIOBlockTable.h"            // M99 blocks_v1 walkers
#import "Genomics/TTIOBlockView.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOGenomicWriteContext.h"
#import "Genomics/TTIOGenomicStreamWriter.h"  // indexFields for the restore fallback
#import "Genomics/TTIOPackedReference.h"
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Codecs/TTIOReferenceResolver.h"  // @reference_md5s restore path
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"

#include <string.h>


static NSString *const kDomain = @"TTIOPerAUFileErrorDomain";

static NSError *makeErr(NSInteger code, NSString *fmt, ...) NS_FORMAT_FUNCTION(2, 3);
static NSError *makeErr(NSInteger code, NSString *fmt, ...)
{
    va_list args; va_start(args, fmt);
    NSString *m = [[NSString alloc] initWithFormat:fmt arguments:args];
    va_end(args);
    return [NSError errorWithDomain:kDomain code:code
                            userInfo:@{NSLocalizedDescriptionKey: m}];
}


// ---------------------------------------------------------------- helpers

// M98: read an assembly sequences channel written with an optional
// @compression codec attribute (0 or absent = raw bytes). Mirrors
// the decode in TTIOAssemblyGraph.
static NSData *decodeAssemblySequences(id<TTIOStorageDataset> ds,
                                       NSError **error)
{
    NSData *raw = (NSData *)[ds readAll:error];
    if (!raw) return nil;
    uint8_t codec = 0;
    if ([ds hasAttributeNamed:@"compression"]) {
        id v = [ds attributeValueForName:@"compression" error:NULL];
        if ([v isKindOfClass:[NSNumber class]]) {
            codec = (uint8_t)[v unsignedIntegerValue];
        }
    }
    if (codec == 0) return raw;
    id<TTIOCodec> c = [TTIOCodecRegistry codecForId:(TTIOCompression)codec];
    if (!c) {
        if (error) *error = makeErr(2,
            @"assembly sequences channel names unregistered codec %u",
            codec);
        return nil;
    }
    TTIODecodedChannel *dec =
        [c decode:[[TTIOBytesPayload alloc] initWithBytes:raw]
          context:[TTIOCodecContext emptyContext]
            error:error];
    if (![dec isKindOfClass:[TTIODecodedBytes class]]) return nil;
    return ((TTIODecodedBytes *)dec).data;
}

static NSString *readStringAttr(id<TTIOStorageGroup> g, NSString *name)
{
    if (!g || ![g hasAttributeNamed:name]) return nil;
    id v = [g attributeValueForName:name error:NULL];
    if ([v isKindOfClass:[NSString class]]) return (NSString *)v;
    if ([v isKindOfClass:[NSData class]]) {
        return [[NSString alloc] initWithData:(NSData *)v
                                      encoding:NSUTF8StringEncoding];
    }
    return nil;
}

static int64_t readIntAttr(id<TTIOStorageGroup> g, NSString *name, int64_t defaultValue)
{
    if (!g || ![g hasAttributeNamed:name]) return defaultValue;
    id v = [g attributeValueForName:name error:NULL];
    if ([v respondsToSelector:@selector(longLongValue)])
        return [v longLongValue];
    return defaultValue;
}


static NSArray<NSString *> *splitChannelNames(NSString *raw)
{
    if (!raw.length) return @[];
    NSArray *parts = [raw componentsSeparatedByString:@","];
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:parts.count];
    for (NSString *p in parts) {
        if (p.length) [out addObject:p];
    }
    return out;
}

static NSArray<NSString *> *listGroupChildren(id<TTIOStorageGroup> g)
{
    NSMutableArray *out = [NSMutableArray array];
    for (NSString *n in [g childNames]) {
        if (![n hasPrefix:@"_"]) [out addObject:n];
    }
    return out;
}


// ---------------------------------------------------------------- feature flags

// Provider-agnostic feature-flag read/write. Mirrors the Python
// _hdf5_io.read_feature_flags / write_feature_flags helpers but in
// ObjC on id<TTIOStorageGroup>.

static NSArray<NSString *> *readFeatureFlags(id<TTIOStorageGroup> root,
                                                NSString **outVersion)
{
    if (outVersion) {
        *outVersion = readStringAttr(root, @"ttio_format_version") ?: @"1.0.0";
    }
    if (![root hasAttributeNamed:@"ttio_features"]) return @[];
    NSString *json = readStringAttr(root, @"ttio_features");
    if (!json.length) return @[];
    NSData *d = [json dataUsingEncoding:NSUTF8StringEncoding];
    id parsed = [NSJSONSerialization JSONObjectWithData:d options:0 error:NULL];
    if (![parsed isKindOfClass:[NSArray class]]) return @[];
    return (NSArray *)parsed;
}

static BOOL writeFeatureFlags(id<TTIOStorageGroup> root,
                                NSString *version,
                                NSArray<NSString *> *features,
                                NSError **error)
{
    if (![root setAttributeValue:version forName:@"ttio_format_version"
                            error:error]) return NO;
    NSData *json = [NSJSONSerialization dataWithJSONObject:features
                                                    options:0 error:error];
    if (!json) return NO;
    NSString *s = [[NSString alloc] initWithData:json
                                          encoding:NSUTF8StringEncoding];
    return [root setAttributeValue:s forName:@"ttio_features" error:error];
}


// ---------------------------------------------------------------- segment I/O

static NSArray<TTIOCompoundField *> *channelSegmentsFields(void)
{
    return @[
        [TTIOCompoundField fieldWithName:@"offset" kind:TTIOCompoundFieldKindInt64],
        [TTIOCompoundField fieldWithName:@"length" kind:TTIOCompoundFieldKindUInt32],
        [TTIOCompoundField fieldWithName:@"iv" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"tag" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"ciphertext" kind:TTIOCompoundFieldKindVLBytes],
    ];
}

static NSArray<TTIOCompoundField *> *auHeaderSegmentsFields(void)
{
    return @[
        [TTIOCompoundField fieldWithName:@"iv" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"tag" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"ciphertext" kind:TTIOCompoundFieldKindVLBytes],
    ];
}

static BOOL writeChannelSegments(id<TTIOStorageGroup> parent,
                                    NSString *name,
                                    NSArray<TTIOChannelSegment *> *segments,
                                    NSError **error)
{
    if ([parent hasChildNamed:name]) {
        if (![parent deleteChildNamed:name error:error]) return NO;
    }

    // Column-oriented write: build the 5 columns in one pass over the
    // segments and hand them to TTIOCompoundIO writeColumnar, avoiding the
    // per-row NSDictionary + NSNumber boxing of the old writeGeneric path.
    //   - offset/length -> packed C arrays (int64 / uint32), zero boxing.
    //   - iv/tag/ciphertext -> NSArray<NSData*> referencing seg.* directly
    //     (no copy; writeColumnar borrows the bytes zero-copy and keeps the
    //     NSData alive through the single bulk H5Dwrite).
    NSUInteger n = segments.count;
    NSMutableData *offsetCol = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    NSMutableData *lengthCol = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    int64_t  *offsetPtr = (int64_t  *)offsetCol.mutableBytes;
    uint32_t *lengthPtr = (uint32_t *)lengthCol.mutableBytes;
    NSMutableArray<NSData *> *ivCol  = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSData *> *tagCol = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSData *> *ctCol  = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        TTIOChannelSegment *seg = segments[i];
        offsetPtr[i] = (int64_t)seg.offset;
        lengthPtr[i] = seg.length;
        [ivCol  addObject:seg.iv];
        [tagCol addObject:seg.tag];
        [ctCol  addObject:seg.ciphertext];
    }

    return [TTIOCompoundIO writeColumnar:@{
                @"offset":     offsetCol,
                @"length":     lengthCol,
                @"iv":         ivCol,
                @"tag":        tagCol,
                @"ciphertext": ctCol,
            }
                              intoGroup:parent
                           datasetNamed:name
                                 fields:channelSegmentsFields()
                                  count:n
                                  error:error];
}

// Compound-dataset reads go through TTIOCompoundIO's
// readGenericFromGroup: path because the StorageDataset protocol
// currently returns a primitive adapter on openDatasetNamed: —
// it doesn't auto-detect H5T_COMPOUND on re-open. Unwrapping the
// StorageGroup to TTIOHDF5Group is the same escape hatch used by
// TTIOSignatureManager. Write goes through the protocol; read uses
// the documented native-handle path.
static NSArray<TTIOChannelSegment *> *readChannelSegments(
    TTIOHDF5Group *hdf5Group, NSString *name, NSError **error)
{
    NSArray<NSDictionary *> *rows =
        [TTIOCompoundIO readGenericFromGroup:hdf5Group
                                  datasetNamed:name
                                        fields:channelSegmentsFields()
                                         error:error];
    if (!rows) return nil;
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:rows.count];
    for (NSDictionary *r in rows) {
        [out addObject:[[TTIOChannelSegment alloc]
            initWithOffset:[r[@"offset"] unsignedLongLongValue]
                     length:[r[@"length"] unsignedIntValue]
                         iv:r[@"iv"]
                        tag:r[@"tag"]
                 ciphertext:r[@"ciphertext"]]];
    }
    return out;
}

static BOOL writeAUHeaderSegments(id<TTIOStorageGroup> parent,
                                    NSString *name,
                                    NSArray<TTIOHeaderSegment *> *segments,
                                    NSError **error)
{
    if ([parent hasChildNamed:name]) {
        if (![parent deleteChildNamed:name error:error]) return NO;
    }
    id<TTIOStorageDataset> ds =
        [parent createCompoundDatasetNamed:name
                                      fields:auHeaderSegmentsFields()
                                       count:segments.count
                                       error:error];
    if (!ds) return NO;
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:segments.count];
    for (TTIOHeaderSegment *seg in segments) {
        [rows addObject:@{
            @"iv": seg.iv, @"tag": seg.tag, @"ciphertext": seg.ciphertext,
        }];
    }
    return [ds writeAll:rows error:error];
}

static NSArray<TTIOHeaderSegment *> *readAUHeaderSegments(
    TTIOHDF5Group *hdf5Group, NSString *name, NSError **error)
{
    NSArray<NSDictionary *> *rows =
        [TTIOCompoundIO readGenericFromGroup:hdf5Group
                                  datasetNamed:name
                                        fields:auHeaderSegmentsFields()
                                         error:error];
    if (!rows) return nil;
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:rows.count];
    for (NSDictionary *r in rows) {
        [out addObject:[[TTIOHeaderSegment alloc]
            initWithIV:r[@"iv"] tag:r[@"tag"] ciphertext:r[@"ciphertext"]]];
    }
    return out;
}


// ---------------------------------------------------------------- impl

// ---------------------------------------------------------------- M90.4 helpers

// Read /study/genomic_runs/<name>/genomic_index/chromosomes as
// NSArray<NSString *>. L1 (Task #82 Phase B.1, 2026-05-01):
// chromosomes are stored as chromosome_ids (uint16) +
// chromosome_names (compound) instead of a single VL-string compound.
static NSArray<NSString *> *readChromosomes(TTIOHDF5Group *hdf5Idx,
                                              NSError **error)
{
    id<TTIOStorageDataset> idsDs = [hdf5Idx openDatasetNamed:@"chromosome_ids" error:error];
    if (!idsDs) return nil;
    NSData *idsData = [idsDs readAll:error];
    if (!idsData) return nil;
    NSArray *fields = @[[TTIOCompoundField fieldWithName:@"name"
                                                     kind:TTIOCompoundFieldKindVLString]];
    NSArray<NSDictionary *> *nameRows =
        [TTIOCompoundIO readGenericFromGroup:hdf5Idx
                                  datasetNamed:@"chromosome_names"
                                        fields:fields
                                         error:error];
    if (!nameRows) return nil;
    NSMutableArray *nameTable = [NSMutableArray arrayWithCapacity:nameRows.count];
    for (NSDictionary *row in nameRows) {
        id v = row[@"name"];
        if ([v isKindOfClass:[NSData class]]) {
            v = [[NSString alloc] initWithData:v
                                       encoding:NSUTF8StringEncoding];
        }
        [nameTable addObject:(NSString *)v ?: @""];
    }
    const uint16_t *ids = (const uint16_t *)idsData.bytes;
    NSUInteger nIds = idsData.length / sizeof(uint16_t);
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:nIds];
    for (NSUInteger i = 0; i < nIds; i++) {
        NSUInteger idx = ids[i];
        [out addObject:idx < nameTable.count ? nameTable[idx] : @""];
    }
    return out;
}

// Per-AU dispatch for the M90.4 encrypt path. Reads on chromosomes
// in keyMap are encrypted with that key; reads on other chromosomes
// emit a clear segment (empty IV + plaintext bytes in ciphertext).
static NSArray<TTIOChannelSegment *> *encryptChannelWithDispatch(
    NSData *plaintext,
    const uint64_t *offsets, const uint32_t *lengths,
    NSArray<NSString *> *chromosomes,
    NSUInteger nReads,
    uint16_t datasetId, NSString *channelName,
    NSDictionary<NSString *, NSData *> *keyMap,
    NSError **error)
{
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:nReads];
    const uint8_t *all = (const uint8_t *)plaintext.bytes;
    for (NSUInteger i = 0; i < nReads; i++) {
        NSUInteger byteOffset = (NSUInteger)offsets[i];
        NSUInteger byteLength = (NSUInteger)lengths[i];
        NSData *chunk = [NSData dataWithBytes:all + byteOffset
                                        length:byteLength];
        NSString *chrom = (i < chromosomes.count) ? chromosomes[i] : @"";
        NSData *key = keyMap[chrom];
        if (!key) {
            // Clear segment: empty IV + tag, plaintext rides in
            // the ciphertext slot.
            [out addObject:[[TTIOChannelSegment alloc]
                initWithOffset:offsets[i]
                         length:lengths[i]
                             iv:[NSData data]
                            tag:[NSData data]
                     ciphertext:chunk]];
        } else {
            NSData *aad = [TTIOPerAUEncryption aadForChannel:channelName
                                                    datasetId:datasetId
                                                   auSequence:(uint32_t)i];
            NSData *iv = [TTIOPerAUEncryption randomIVWithError:error];
            if (!iv) return nil;
            NSData *tag = nil;
            NSData *ct = [TTIOPerAUEncryption encryptWithPlaintext:chunk
                                                                key:key
                                                                 iv:iv
                                                                aad:aad
                                                             outTag:&tag
                                                              error:error];
            if (!ct) return nil;
            [out addObject:[[TTIOChannelSegment alloc]
                initWithOffset:offsets[i]
                         length:lengths[i]
                             iv:iv tag:tag ciphertext:ct]];
        }
    }
    return out;
}

// reserved key name for encrypting genomic_index columns
// under encryptFilePathByRegion's keyMap. The presence of this entry
// is the opt-in signal for opt_encrypted_au_headers on genomic data.
static NSString *const kTTIOPerAUHeadersKeyName = @"_headers";

// serialise the chromosomes list as compact JSON
// ``["chr1","chr1","chr2"]`` — must match Python's
// json.dumps(chromosomes) which uses double-quoted strings and a
// single comma+space separator. Foundation's NSJSONSerialization
// emits without whitespace by default; we then reconstruct the
// "[\"chr1\", \"chr2\"]" shape Python uses (`json.dumps` default is
// no whitespace either, both produce ``["chr1","chr2"]``). Verified
// against test_m90_11 fixtures.
static NSData *chromosomesToJSON(NSArray<NSString *> *chromosomes)
{
    NSError *jsonErr = nil;
    NSData *d = [NSJSONSerialization dataWithJSONObject:chromosomes
                                                  options:0
                                                    error:&jsonErr];
    return d ?: [NSData data];
}

// M90.11 encrypt: replace the four plaintext genomic_index columns
// (chromosomes, positions, mapping_qualities, flags) with
// ``<col>_encrypted`` uint8 1-D blobs containing iv || tag || ct.
// AAD = "genomic_headers:<dataset_id>:<col>" (ASCII).
// offsets/lengths stay plaintext — structural framing.
static BOOL encryptGenomicIndex(id<TTIOStorageGroup> gIdx,
                                  uint16_t datasetId,
                                  NSData *key,
                                  NSArray<NSString *> *chromosomes,
                                  NSError **error)
{
    // Read the four columns into raw little-endian byte buffers.
    NSData *chrJson = chromosomesToJSON(chromosomes);

    id<TTIOStorageDataset> posDs =
        [gIdx openDatasetNamed:@"positions" error:error];
    if (!posDs) return NO;
    NSData *posBytes = [posDs readAll:error];
    if (!posBytes) return NO;

    id<TTIOStorageDataset> mqDs =
        [gIdx openDatasetNamed:@"mapping_qualities" error:error];
    if (!mqDs) return NO;
    NSData *mqBytes = [mqDs readAll:error];
    if (!mqBytes) return NO;

    id<TTIOStorageDataset> flagsDs =
        [gIdx openDatasetNamed:@"flags" error:error];
    if (!flagsDs) return NO;
    NSData *flagsBytes = [flagsDs readAll:error];
    if (!flagsBytes) return NO;

    NSDictionary<NSString *, NSData *> *columns = @{
        @"chromosomes": chrJson,
        @"positions": posBytes,
        @"mapping_qualities": mqBytes,
        @"flags": flagsBytes,
    };

    // Encrypt + write each column as ``<col>_encrypted`` uint8 blob.
    // Iterate in a fixed order so this routine is deterministic.
    NSArray<NSString *> *order = @[
        @"chromosomes", @"positions", @"mapping_qualities", @"flags",
    ];
    for (NSString *colName in order) {
        NSData *plaintext = columns[colName];
        NSString *aadStr =
            [NSString stringWithFormat:@"genomic_headers:%u:%@",
                (unsigned)datasetId, colName];
        NSData *aad = [aadStr dataUsingEncoding:NSASCIIStringEncoding];
        NSData *iv = [TTIOPerAUEncryption randomIVWithError:error];
        if (!iv) return NO;
        NSData *tag = nil;
        NSData *ct = [TTIOPerAUEncryption encryptWithPlaintext:plaintext
                                                            key:key
                                                             iv:iv
                                                            aad:aad
                                                         outTag:&tag
                                                          error:error];
        if (!ct) return NO;

        // Concat iv (12) || tag (16) || ciphertext into a single blob.
        NSMutableData *blob = [NSMutableData dataWithCapacity:
            iv.length + tag.length + ct.length];
        [blob appendData:iv];
        [blob appendData:tag];
        [blob appendData:ct];

        // Delete plaintext column (or any pre-existing _encrypted dataset).
        if ([gIdx hasChildNamed:colName]) {
            if (![gIdx deleteChildNamed:colName error:error]) return NO;
        }
        // L1 (Task #82 Phase B.1): the on-disk chromosomes column is
        // decomposed into chromosome_ids + chromosome_names — also
        // delete those when encrypting the logical "chromosomes"
        // column so plaintext doesn't linger alongside the encrypted
        // blob.
        if ([colName isEqualToString:@"chromosomes"]) {
            for (NSString *sub in @[@"chromosome_ids", @"chromosome_names"]) {
                if ([gIdx hasChildNamed:sub]) {
                    if (![gIdx deleteChildNamed:sub error:error]) return NO;
                }
            }
        }
        NSString *encName =
            [NSString stringWithFormat:@"%@_encrypted", colName];
        if ([gIdx hasChildNamed:encName]) {
            if (![gIdx deleteChildNamed:encName error:error]) return NO;
        }
        id<TTIOStorageDataset> outDs =
            [gIdx createDatasetNamed:encName
                            precision:TTIOPrecisionUInt8
                               length:blob.length
                            chunkSize:0
                          compression:TTIOCompressionNone
                     compressionLevel:0
                                error:error];
        if (!outDs) return NO;
        if (![outDs writeAll:blob error:error]) return NO;
    }
    return YES;
}

// M90.11 decrypt: read the ``<col>_encrypted`` blobs and recover the
// four plaintext columns. Returns a dictionary
// ``{"chromosomes": NSArray<NSString *>, "positions": NSData (i64
// LE), "mapping_qualities": NSData (u1), "flags": NSData (u32 LE)}``
// or nil + NSError on failure.
static NSDictionary *decryptGenomicIndex(id<TTIOStorageGroup> gIdx,
                                            uint16_t datasetId,
                                            NSData *key,
                                            NSError **error)
{
    NSArray<NSString *> *order = @[
        @"chromosomes", @"positions", @"mapping_qualities", @"flags",
    ];
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    for (NSString *colName in order) {
        NSString *encName =
            [NSString stringWithFormat:@"%@_encrypted", colName];
        if (![gIdx hasChildNamed:encName]) {
            if (error) *error = [NSError errorWithDomain:kDomain code:7
                userInfo:@{NSLocalizedDescriptionKey:
                    [NSString stringWithFormat:
                        @"genomic_index/%@ missing — file does not "
                        @"appear to carry M90.11 encrypted headers",
                        encName]}];
            return nil;
        }
        id<TTIOStorageDataset> ds =
            [gIdx openDatasetNamed:encName error:error];
        if (!ds) return nil;
        NSData *blob = [ds readAll:error];
        if (!blob) return nil;
        if (blob.length < 12 + 16) {
            if (error) *error = [NSError errorWithDomain:kDomain code:8
                userInfo:@{NSLocalizedDescriptionKey:
                    [NSString stringWithFormat:
                        @"genomic_index/%@ too short for IV+TAG", encName]}];
            return nil;
        }
        NSData *iv = [blob subdataWithRange:NSMakeRange(0, 12)];
        NSData *tag = [blob subdataWithRange:NSMakeRange(12, 16)];
        NSData *ct = [blob subdataWithRange:NSMakeRange(28, blob.length - 28)];
        NSString *aadStr =
            [NSString stringWithFormat:@"genomic_headers:%u:%@",
                (unsigned)datasetId, colName];
        NSData *aad = [aadStr dataUsingEncoding:NSASCIIStringEncoding];
        NSData *plain = [TTIOPerAUEncryption
            decryptWithCiphertext:ct key:key iv:iv tag:tag aad:aad
                            error:error];
        if (!plain) return nil;
        if ([colName isEqualToString:@"chromosomes"]) {
            NSError *jsonErr = nil;
            id parsed = [NSJSONSerialization JSONObjectWithData:plain
                                                          options:0
                                                            error:&jsonErr];
            if (![parsed isKindOfClass:[NSArray class]]) {
                if (error) *error = [NSError errorWithDomain:kDomain code:9
                    userInfo:@{NSLocalizedDescriptionKey:
                        @"chromosomes column did not decode to JSON array"}];
                return nil;
            }
            out[colName] = parsed;
        } else {
            out[colName] = plain;
        }
    }
    return out;
}

// Per-AU dispatch decrypt — branches on len(seg.iv): 0 = clear
// segment (ciphertext is plaintext), 12 = AES-256-GCM.
static NSData *decryptChannelWithDispatch(
    NSArray<TTIOChannelSegment *> *segments,
    NSArray<NSString *> *chromosomes,
    uint16_t datasetId, NSString *channelName,
    NSDictionary<NSString *, NSData *> *keyMap,
    NSError **error)
{
    NSMutableData *out = [NSMutableData data];
    for (NSUInteger i = 0; i < segments.count; i++) {
        TTIOChannelSegment *seg = segments[i];
        if (seg.iv.length == 0) {
            // Clear segment: ciphertext IS plaintext.
            [out appendData:seg.ciphertext];
            continue;
        }
        NSString *chrom = (i < chromosomes.count) ? chromosomes[i] : @"";
        NSData *key = keyMap[chrom];
        if (!key) {
            if (error) *error = [NSError errorWithDomain:kDomain code:5
                userInfo:@{NSLocalizedDescriptionKey:
                    [NSString stringWithFormat:
                        @"chromosome %@ segment %lu is encrypted but "
                        @"keyMap has no entry for %@",
                        chrom, (unsigned long)i, chrom]}];
            return nil;
        }
        NSData *aad = [TTIOPerAUEncryption aadForChannel:channelName
                                                datasetId:datasetId
                                               auSequence:(uint32_t)i];
        NSData *plain = [TTIOPerAUEncryption decryptWithCiphertext:seg.ciphertext
                                                                key:key
                                                                 iv:seg.iv
                                                                tag:seg.tag
                                                                aad:aad
                                                              error:error];
        if (!plain) return nil;
        if (plain.length != (NSUInteger)seg.length) {
            if (error) *error = [NSError errorWithDomain:kDomain code:6
                userInfo:@{NSLocalizedDescriptionKey:
                    [NSString stringWithFormat:
                        @"channel %@ segment %lu: decrypted %lu bytes, "
                        @"expected %u",
                        channelName, (unsigned long)i,
                        (unsigned long)plain.length, (unsigned)seg.length]}];
            return nil;
        }
        [out appendData:plain];
    }
    return out;
}


// ---------------------------------------------------------------- M99:
// blocks_v1 per-AU walkers. The default genomic layout stores
// codec-coded per-block blobs, so the walkers stream block by block:
// decode one block, slice its reads, encrypt one AU per read with
// GLOBAL AU numbering, append to extendable segments tables. Restore
// re-encodes each block with the stream writer's machinery; because
// writer policy the file does not persist would break that
// reproducibility, ENCRYPT re-encodes and byte-compares every block
// BEFORE deleting anything, and refuses the run when a blob is not
// reproducible.

// M101: the SAM tags are encrypted with the bases they describe (an
// MD string lists the read's mismatches). One AU per read, the
// plaintext being the read's UTF-8 tag text; "tags" is handled as
// text, not bytes.
static NSArray<NSString *> *kBlocksV1Channels(void)
{
    return @[@"sequences", @"qualities", @"tags"];
}

// Split a decrypted flat tags plaintext into one string per AU, using
// the segments' plaintext lengths.
static NSArray<NSString *> *ttioSplitTagTexts(NSData *flat,
                                              NSArray<TTIOChannelSegment *> *segs,
                                              NSError **error)
{
    NSMutableArray<NSString *> *out = [NSMutableArray arrayWithCapacity:segs.count];
    const uint8_t *p = (const uint8_t *)flat.bytes;
    NSUInteger at = 0;
    for (TTIOChannelSegment *s in segs) {
        NSUInteger ln = s.length;
        if (at + ln > flat.length) {
            if (error) *error = makeErr(5, @"per-AU tags: segment lengths exceed the plaintext");
            return nil;
        }
        NSString *t = ln ? [[NSString alloc] initWithBytes:p + at length:ln
                                                  encoding:NSUTF8StringEncoding]
                         : @"";
        if (t == nil) {
            if (error) *error = makeErr(5, @"per-AU tags: plaintext is not valid UTF-8");
            return nil;
        }
        [out addObject:t];
        at += ln;
    }
    return out;
}

static BOOL ttioIsBlocksV1(id<TTIOStorageGroup> runGroup)
{
    NSString *layout = readStringAttr(runGroup, @"layout");
    // blocks_v1_grouped (M103) is blocks_v1 plus genomic_index/input_index;
    // per-AU work walks stored rows block by block and leaves the column
    // as it is.
    return layout != nil && ([layout isEqualToString:@"blocks_v1"]
                             || [layout isEqualToString:@"blocks_v1_grouped"]);
}

// {chromosome: bytes} of an embedded reference; nil when absent.
static NSDictionary<NSString *, NSData *> *
ttioEmbeddedReferenceSeqs(id<TTIOStorageGroup> study, NSString *uri)
{
    if (uri.length == 0 || ![study hasChildNamed:@"references"]) return nil;
    id<TTIOStorageGroup> refs = [study openGroupNamed:@"references" error:NULL];
    if (!refs || ![refs hasChildNamed:uri]) return nil;
    id<TTIOStorageGroup> ref = [refs openGroupNamed:uri error:NULL];
    if (!ref || ![ref hasChildNamed:@"chromosomes"]) return nil;
    id<TTIOStorageGroup> chroms = [ref openGroupNamed:@"chromosomes" error:NULL];
    if (!chroms) return nil;
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    for (NSString *cname in [chroms childNames]) {
        id<TTIOStorageGroup> cg = [chroms openGroupNamed:cname error:NULL];
        NSData *seq = cg
            ? [TTIOPackedReference readChromosomeBytes:cg error:NULL] : nil;
        if (!seq) return nil;
        out[cname] = seq;
    }
    return out;
}

// The run-wide chromosome-id map the writer accumulated: the
// mate_info/chrom_names table dumped in row order (row index = id).
static NSMutableDictionary<NSString *, NSNumber *> *
ttioBlocksV1ChromMap(id<TTIOStorageGroup> sig)
{
    NSMutableDictionary *map = [NSMutableDictionary dictionary];
    if ([sig hasChildNamed:@"mate_info"]) {
        id<TTIOStorageGroup> mate = [sig openGroupNamed:@"mate_info" error:NULL];
        NSArray<NSString *> *names = mate
            ? [TTIOBlockView readNamesIn:mate named:@"chrom_names"] : @[];
        for (NSUInteger i = 0; i < names.count; i++) map[names[i]] = @(i);
    }
    return map;
}

// {chromosome: bytes} for re-encoding a run's REF_DIFF blocks. With
// the @reference_md5s run attribute the reference is rebuilt through
// TTIOReferenceResolver (embedded or REF_PATH); older files without
// the attribute fall back to the embedded-directory walk. Returns
// nil when neither works.
static NSDictionary<NSString *, NSData *> *
ttioBlocksV1ReferenceSeqs(id<TTIOStorageGroup> runGroup,
                          id<TTIOStorageGroup> study,
                          NSString *uri,
                          NSError **error)
{
    NSString *raw = readStringAttr(runGroup, @"reference_md5s");
    if (raw.length > 0) {
        NSDictionary *md5s = [NSJSONSerialization
            JSONObjectWithData:[raw dataUsingEncoding:NSUTF8StringEncoding]
                       options:0
                         error:error];
        if (![md5s isKindOfClass:[NSDictionary class]]) return nil;
        if (![runGroup respondsToSelector:@selector(unwrap)]) return nil;
        TTIOHDF5Group *runG = [(id)runGroup performSelector:@selector(unwrap)];
        hid_t fid = H5Iget_file_id([runG groupId]);
        if (fid < 0) return nil;
        TTIOHDF5Group *rootHDF5 = nil;
        hid_t rootId = H5Gopen2(fid, "/", H5P_DEFAULT);
        if (rootId >= 0) {
            rootHDF5 = [[TTIOHDF5Group alloc] initWithGroupId:rootId
                                                     retainer:nil];
        }
        H5Idec_ref(fid);
        if (rootHDF5 == nil) return nil;
        TTIOReferenceResolver *resolver = [[TTIOReferenceResolver alloc]
            initWithRootGroup:rootHDF5 externalReferencePath:nil];
        NSMutableDictionary *out =
            [NSMutableDictionary dictionaryWithCapacity:md5s.count];
        for (NSString *chrom in md5s) {
            NSString *hex = md5s[chrom];
            if (![hex isKindOfClass:[NSString class]]
                || hex.length != 32) return nil;
            uint8_t md5[16];
            for (NSUInteger i = 0; i < 16; i++) {
                unsigned v = 0;
                if (sscanf([[hex substringWithRange:NSMakeRange(i * 2, 2)]
                               UTF8String], "%02x", &v) != 1) return nil;
                md5[i] = (uint8_t)v;
            }
            NSData *seq = [resolver
                resolveURI:uri
               expectedMD5:[NSData dataWithBytes:md5 length:16]
                chromosome:chrom
                     error:error];
            if (!seq) return nil;
            out[chrom] = seq;
        }
        return out;
    }
    return ttioEmbeddedReferenceSeqs(study, uri);
}

// Collect reads [indexBase, indexBase+nn) from an open reader into a
// per-block TTIOWrittenGenomicRun; outSeq/outQual receive the block's
// decoded plaintext channel bytes. The encrypt walker reads block b
// through the run's own reader (indexBase = read_start[b]); the
// decrypt walker reads a materialised one-block view (indexBase = 0).
// Rows are STORED rows: a grouped run's reader presents input order
// (M103), so stored row r is read as input index inputIndex[r], which
// lands in block b and is served by the reader's cached block view.
static TTIOWrittenGenomicRun *
ttioBlocksV1BlockRun(TTIOGenomicRun *rd,
                     id<TTIOStorageGroup> runGroup,
                     id<TTIOStorageGroup> study,
                     TTIOBlockTable *table,
                     NSUInteger b,
                     unsigned long long indexBase,
                     NSData **outSeq,
                     NSData **outQual,
                     NSError **error)
{
    unsigned long long r0 = indexBase;
    NSUInteger nn = [table nReadsAt:b];
    NSMutableData *positions = [NSMutableData dataWithLength:nn * 8];
    NSMutableData *mapqs = [NSMutableData dataWithLength:nn];
    NSMutableData *flags = [NSMutableData dataWithLength:nn * 4];
    NSMutableData *offsets = [NSMutableData dataWithLength:nn * 8];
    NSMutableData *lengths = [NSMutableData dataWithLength:nn * 4];
    NSMutableData *matePos = [NSMutableData dataWithLength:nn * 8];
    NSMutableData *tlens = [NSMutableData dataWithLength:nn * 4];
    int64_t  *posP  = positions.mutableBytes;
    uint8_t  *mapqP = mapqs.mutableBytes;
    uint32_t *flagP = flags.mutableBytes;
    uint64_t *offP  = offsets.mutableBytes;
    uint32_t *lenP  = lengths.mutableBytes;
    int64_t  *mposP = matePos.mutableBytes;
    int32_t  *tlenP = tlens.mutableBytes;
    NSMutableArray *cigars = [NSMutableArray arrayWithCapacity:nn];
    NSMutableArray *readNames = [NSMutableArray arrayWithCapacity:nn];
    NSMutableArray *mateChroms = [NSMutableArray arrayWithCapacity:nn];
    NSMutableArray *chroms = [NSMutableArray arrayWithCapacity:nn];
    NSMutableArray *tags = [NSMutableArray arrayWithCapacity:nn];
    NSMutableData *seq = [NSMutableData data];
    NSMutableData *qual = [NSMutableData data];
    NSData *inputIndex = rd.inputIndex;
    const uint32_t *ii = (const uint32_t *)inputIndex.bytes;
    NSUInteger nii = inputIndex.length / sizeof(uint32_t);
    for (NSUInteger i = 0; i < nn; i++) {
        @autoreleasepool {
            NSUInteger row = (NSUInteger)(r0 + i);
            if (ii != NULL && row < nii) row = ii[row];
            TTIOAlignedRead *r = [rd readAtIndex:row error:error];
            if (!r) return nil;
            offP[i] = (uint64_t)seq.length;
            NSData *sBytes =
                [r.sequence dataUsingEncoding:NSASCIIStringEncoding];
            [seq appendData:sBytes ?: [NSData data]];
            [qual appendData:r.qualities ?: [NSData data]];
            lenP[i] = (uint32_t)(sBytes ? sBytes.length : 0);
            posP[i] = r.position;
            mapqP[i] = r.mappingQuality;
            flagP[i] = r.flags;
            mposP[i] = r.matePosition;
            tlenP[i] = r.templateLength;
            [cigars addObject:r.cigar ?: @""];
            [readNames addObject:r.readName ?: @""];
            [mateChroms addObject:r.mateChromosome ?: @""];
            [chroms addObject:r.chromosome ?: @""];
            [tags addObject:r.tags ?: @""];
        }
    }

    NSString *refUri = readStringAttr(runGroup, @"reference_uri") ?: @"";
    NSUInteger seqCodec = table.hasCodecs ? [table codecOf:@"sequences" at:b] : 0;
    NSUInteger qualCodec = table.hasCodecs ? [table codecOf:@"qualities" at:b] : 0;
    NSDictionary<NSString *, NSData *> *refSeqs = nil;
    NSMutableDictionary *overrides = [NSMutableDictionary dictionary];
    if (seqCodec == TTIOCompressionRefDiffV2) {
        refSeqs = ttioBlocksV1ReferenceSeqs(runGroup, study, refUri, error);
        if (refSeqs == nil) {
            if (error && *error == nil) *error = makeErr(4,
                @"per-AU blocks_v1: block %lu codes sequences with "
                @"REF_DIFF_V2 but reference '%@' is not embedded in "
                @"/study/references and not resolvable via REF_PATH; "
                @"restoring the blob needs the reference bytes",
                (unsigned long)b, refUri);
            return nil;
        }
    } else if (seqCodec != 0) {
        overrides[@"sequences"] = @(seqCodec);
    }
    if (qualCodec != 0) overrides[@"qualities"] = @(qualCodec);

    TTIOWrittenGenomicRun *block = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:(TTIOAcquisitionMode)readIntAttr(runGroup, @"acquisition_mode", 0)
                   referenceUri:refUri
                       platform:readStringAttr(runGroup, @"platform") ?: @""
                     sampleName:readStringAttr(runGroup, @"sample_name") ?: @""
                      positions:positions
               mappingQualities:mapqs
                          flags:flags
                      sequences:seq
                      qualities:qual
                        offsets:offsets
                        lengths:lengths
                         cigars:cigars
                      readNames:readNames
                mateChromosomes:mateChroms
                  matePositions:matePos
                templateLengths:tlens
                    chromosomes:chroms
              signalCompression:TTIOCompressionZlib
           signalCodecOverrides:overrides];
    if (refSeqs != nil) block.referenceChromSeqs = refSeqs;
    // M101: the decrypt walker overrides these with the decrypted text.
    block.tags = tags;
    NSString *role = readStringAttr(runGroup, @"read_role");
    if (role.length > 0) block.readRole = role;
    int64_t slice = readIntAttr(runGroup, @"ref_diff_slice_bytes", 0);
    if (slice > 0) block.refDiffSliceBytes = (unsigned long long)slice;
    if (readIntAttr(runGroup, @"opt_disable_qualities_v5", 0) != 0) {
        block.optDisableQualitiesV5 = YES;
    }
    if (outSeq) *outSeq = seq;
    if (outQual) *outQual = qual;
    return block;
}

// The stream writer's sticky qualities discipline: after the first
// FQZCOMP_NX16_Z block, read the winning strategy back from the
// encoded stream and pin it for the rest of the run.
static NSInteger ttioBlocksV1DeriveQualHint(TTIOBlockBlobs *blobs,
                                            NSInteger current)
{
    if (current != -1) return current;
    if ([blobs.codecs[@"qualities"] unsignedIntegerValue]
            != TTIOCompressionFqzcompNx16Z) return current;
    NSInteger strat =
        [TTIOFqzcompNx16Z strategyOfEncodedStream:blobs.blobs[@"qualities"]];
    if (strat <= 0) return current;
    return strat == 4 ? TTIOM94ZHintV4Auto : strat;
}

static BOOL ttioEncryptBlocksV1Run(id<TTIOStorageGroup> study,
                                   id<TTIOStorageGroup> runGroup,
                                   NSString *runName,
                                   uint16_t datasetId,
                                   NSData *key,
                                   NSError **error)
{
    TTIOBlockTable *table = [TTIOBlockTable readFromRunGroup:runGroup
                                                       error:error];
    if (!table) return NO;
    id<TTIOStorageGroup> sig =
        [runGroup openGroupNamed:@"signal_channels" error:error];
    if (!sig) return NO;
    NSMutableArray<NSString *> *channels = [NSMutableArray array];
    for (NSString *ch in kBlocksV1Channels()) {
        if ([sig hasChildNamed:ch]) [channels addObject:ch];
    }
    if (channels.count == 0) return YES;
    if (![sig respondsToSelector:@selector(unwrap)]) {
        if (error) *error = makeErr(3,
            @"per-AU blocks_v1 requires the HDF5 provider");
        return NO;
    }
    TTIOHDF5Group *hdf5Sig = [(id)sig performSelector:@selector(unwrap)];

    TTIOGenomicRun *rd = [TTIOGenomicRun openFromGroup:runGroup
                                                  name:runName
                                                 error:error];
    if (!rd) return NO;

    for (NSString *ch in channels) {
        NSString *segName = [NSString stringWithFormat:@"%@_segments", ch];
        if ([sig hasChildNamed:segName]
            && ![sig deleteChildNamed:segName error:error]) return NO;
        if (![TTIOCompoundIO createExtendableCompoundInGroup:hdf5Sig
                                                        name:segName
                                                      fields:channelSegmentsFields()
                                                   chunkRows:1024
                                                       error:error]) return NO;
    }

    BOOL ok = YES;
    unsigned long long tagOffset = 0;   // global plaintext offset of the tags channel
    for (NSUInteger b = 0; ok && b < table.count; b++) {
        @autoreleasepool {
            NSData *seq = nil, *qual = nil;
            TTIOWrittenGenomicRun *block = ttioBlocksV1BlockRun(
                rd, runGroup, study, table, b,
                [table readStartAt:b], &seq, &qual, error);
            if (!block) { ok = NO; break; }

            unsigned long long r0 = [table readStartAt:b];
            NSUInteger nn = [table nReadsAt:b];
            const uint64_t *localOff =
                (const uint64_t *)block.offsetsData.bytes;
            const uint32_t *localLen =
                (const uint32_t *)block.lengthsData.bytes;
            NSDictionary *plain = @{@"sequences": seq ?: [NSData data],
                                    @"qualities": qual ?: [NSData data]};
            for (NSString *ch in channels) {
                NSData *chPlain = plain[ch];
                const uint64_t *chOff = localOff;
                const uint32_t *chLen = localLen;
                unsigned long long base = [table baseStartAt:b];
                // Precise lifetime: only interior pointers are used below.
                __attribute__((objc_precise_lifetime)) NSMutableData *tagOff = nil;
                __attribute__((objc_precise_lifetime)) NSMutableData *tagLen = nil;
                if ([ch isEqualToString:@"tags"]) {
                    // One AU per read: the read's UTF-8 tag text, with
                    // the channel's own running global offset.
                    NSMutableData *flat = [NSMutableData data];
                    tagOff = [NSMutableData dataWithLength:nn * sizeof(uint64_t)];
                    tagLen = [NSMutableData dataWithLength:nn * sizeof(uint32_t)];
                    uint64_t *to = (uint64_t *)tagOff.mutableBytes;
                    uint32_t *tl = (uint32_t *)tagLen.mutableBytes;
                    for (NSUInteger i = 0; i < nn; i++) {
                        NSString *t = i < block.tags.count ? block.tags[i] : @"";
                        NSData *tb = [t dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
                        to[i] = (uint64_t)flat.length;
                        tl[i] = (uint32_t)tb.length;
                        [flat appendData:tb];
                    }
                    chPlain = flat;
                    chOff = to;
                    chLen = tl;
                    base = tagOffset;
                    tagOffset += flat.length;
                }
                NSArray<TTIOChannelSegment *> *segs =
                    [TTIOPerAUEncryption
                        encryptChannelToSegments:chPlain
                                          offsets:chOff
                                          lengths:chLen
                                         nSpectra:nn
                                  bytesPerElement:1
                                        datasetId:datasetId
                                      channelName:ch
                                              key:key
                                           auBase:(uint32_t)r0
                                       offsetBase:base
                                            error:error];
                if (!segs) { ok = NO; break; }
                NSMutableArray *rows =
                    [NSMutableArray arrayWithCapacity:segs.count];
                for (TTIOChannelSegment *s in segs) {
                    [rows addObject:@{@"offset": @((int64_t)s.offset),
                                      @"length": @(s.length),
                                      @"iv": s.iv,
                                      @"tag": s.tag,
                                      @"ciphertext": s.ciphertext}];
                }
                NSString *segName =
                    [NSString stringWithFormat:@"%@_segments", ch];
                if (![TTIOCompoundIO appendRows:rows
                                         toGroup:hdf5Sig
                                            name:segName
                                          fields:channelSegmentsFields()
                                           error:error]) { ok = NO; break; }
            }
        }
    }
    if (!ok) {
        for (NSString *ch in channels) {
            NSString *segName = [NSString stringWithFormat:@"%@_segments", ch];
            [sig deleteChildNamed:segName error:NULL];
        }
        return NO;
    }
    for (NSString *ch in channels) {
        if (![sig deleteChildNamed:ch error:error]) return NO;
        if (![sig setAttributeValue:@"aes-256-gcm"
                             forName:[NSString stringWithFormat:@"%@_algorithm", ch]
                               error:error]) return NO;
    }
    return YES;
}

// Fallback when a re-encoded blob does not land on the ranges the
// block index records: patch off/len/codec for the re-encoded
// channels and recreate blocks/index; the other columns and
// channels are carried over unchanged.
static BOOL ttioBlocksV1RewriteIndex(id<TTIOStorageGroup> runGroup,
                                     NSArray<NSString *> *channels,
                                     NSDictionary *newOff,
                                     NSDictionary *newLen,
                                     NSDictionary *newCodec,
                                     NSError **error)
{
    id<TTIOStorageGroup> blocks =
        [runGroup openGroupNamed:@"blocks" error:error];
    if (!blocks) return NO;
    id<TTIOStorageDataset> ds =
        [blocks openDatasetNamed:@"index" error:error];
    if (!ds) return NO;
    NSArray<NSDictionary *> *rows = [ds readRows:error];
    if (!rows) return NO;
    NSMutableArray *patched =
        [NSMutableArray arrayWithCapacity:rows.count];
    for (NSUInteger b = 0; b < rows.count; b++) {
        NSMutableDictionary *row = [rows[b] mutableCopy];
        for (NSString *ch in channels) {
            row[[ch stringByAppendingString:@"_off"]] = newOff[ch][b];
            row[[ch stringByAppendingString:@"_len"]] = newLen[ch][b];
            row[[ch stringByAppendingString:@"_codec"]] = newCodec[ch][b];
        }
        [patched addObject:row];
    }
    if (![blocks deleteChildNamed:@"index" error:error]) return NO;
    id<TTIOStorageDataset> out = [blocks
        createCompoundDatasetNamed:@"index"
                            fields:[TTIOGenomicStreamWriter indexFieldsForRows:patched]
                             count:0
                        extendable:YES
                         chunkRows:1024
                             error:error];
    if (!out) return NO;
    return patched.count == 0 || [out appendData:patched error:error];
}

static BOOL ttioDecryptBlocksV1RunInPlace(id<TTIOStorageGroup> study,
                                          id<TTIOStorageGroup> runGroup,
                                          uint16_t datasetId,
                                          NSData *key,
                                          NSError **error)
{
    TTIOBlockTable *table = [TTIOBlockTable readFromRunGroup:runGroup
                                                       error:error];
    if (!table) return NO;
    id<TTIOStorageGroup> sig =
        [runGroup openGroupNamed:@"signal_channels" error:error];
    if (!sig) return NO;
    NSMutableArray<NSString *> *channels = [NSMutableArray array];
    for (NSString *ch in kBlocksV1Channels()) {
        if ([sig hasChildNamed:
                [NSString stringWithFormat:@"%@_segments", ch]]) {
            [channels addObject:ch];
        }
    }
    if (channels.count == 0) return YES;
    if (![sig respondsToSelector:@selector(unwrap)]) {
        if (error) *error = makeErr(3,
            @"per-AU blocks_v1 requires the HDF5 provider");
        return NO;
    }
    TTIOHDF5Group *hdf5Sig = [(id)sig performSelector:@selector(unwrap)];

    id<TTIOStorageGroup> idxGroup =
        [runGroup openGroupNamed:@"genomic_index" error:error];
    if (!idxGroup) return NO;
    NSArray<NSString *> *chromNames =
        [TTIOBlockView readNamesIn:idxGroup named:@"chromosome_names"];
    NSArray<NSString *> *mateChromNames = @[];
    if ([sig hasChildNamed:@"mate_info"]) {
        id<TTIOStorageGroup> mate = [sig openGroupNamed:@"mate_info" error:NULL];
        if (mate) mateChromNames =
            [TTIOBlockView readNamesIn:mate named:@"chrom_names"];
    }
    NSMutableDictionary *chromMap = ttioBlocksV1ChromMap(sig);
    NSSet *skip = [NSSet setWithArray:channels];

    NSMutableDictionary *newDs = [NSMutableDictionary dictionary];
    NSMutableDictionary *written = [NSMutableDictionary dictionary];
    NSMutableDictionary *newOff = [NSMutableDictionary dictionary];
    NSMutableDictionary *newLen = [NSMutableDictionary dictionary];
    NSMutableDictionary *newCodec = [NSMutableDictionary dictionary];
    for (NSString *ch in channels) {
        written[ch] = @(0ULL);
        newOff[ch] = [NSMutableArray array];
        newLen[ch] = [NSMutableArray array];
        newCodec[ch] = [NSMutableArray array];
    }
    BOOL mismatch = NO;
    NSInteger qualHint = -1;
    NSData *refMD5 = nil;
    for (NSUInteger b = 0; b < table.count; b++) {
        @autoreleasepool {
            unsigned long long r0 = [table readStartAt:b];
            NSUInteger nn = [table nReadsAt:b];
            NSMutableDictionary *decrypted = [NSMutableDictionary dictionary];
            for (NSString *ch in channels) {
                NSString *segName =
                    [NSString stringWithFormat:@"%@_segments", ch];
                NSArray<NSDictionary *> *rows =
                    [TTIOCompoundIO readGenericFromGroup:hdf5Sig
                                             datasetNamed:segName
                                                   fields:channelSegmentsFields()
                                                   offset:(NSUInteger)r0
                                                    count:nn
                                                    error:error];
                if (!rows) return NO;
                NSMutableArray *segs =
                    [NSMutableArray arrayWithCapacity:rows.count];
                for (NSDictionary *r in rows) {
                    [segs addObject:[[TTIOChannelSegment alloc]
                        initWithOffset:[r[@"offset"] unsignedLongLongValue]
                                 length:[r[@"length"] unsignedIntValue]
                                     iv:r[@"iv"]
                                    tag:r[@"tag"]
                             ciphertext:r[@"ciphertext"]]];
                }
                NSData *plain = [TTIOPerAUEncryption
                    decryptChannelFromSegments:segs
                                bytesPerElement:1
                                      datasetId:datasetId
                                    channelName:ch
                                            key:key
                                         auBase:(uint32_t)r0
                                          error:error];
                if (!plain) return NO;
                if ([ch isEqualToString:@"tags"]) {
                    // M101: tags come back as one string per read.
                    NSArray *texts = ttioSplitTagTexts(plain, segs, error);
                    if (!texts) return NO;
                    decrypted[ch] = texts;
                } else {
                    decrypted[ch] = plain;
                }
            }

            TTIOBlockView *view =
                [TTIOBlockView materialiseBlock:b
                                          ofRun:runGroup
                                          table:table
                                     chromNames:chromNames
                                 mateChromNames:mateChromNames
                                   skipChannels:skip
                                          error:error];
            if (!view) return NO;
            id<TTIOStorageGroup> viewSig =
                [view.group openGroupNamed:@"signal_channels" error:error];
            if (!viewSig) { [view discard]; return NO; }
            for (NSString *ch in channels) {
                // The tags go straight into the block run, never into
                // the view as a codec-0 dataset.
                if ([ch isEqualToString:@"tags"]) continue;
                NSData *raw = decrypted[ch];
                id<TTIOStorageDataset> ds =
                    [viewSig createDatasetNamed:ch
                                      precision:TTIOPrecisionUInt8
                                         length:raw.length
                                      chunkSize:65536
                                    compression:TTIOCompressionNone
                               compressionLevel:0
                                          error:error];
                if (!ds || ![ds writeAll:raw error:error]
                    || ![ds setAttributeValue:@((int64_t)0)
                                       forName:@"compression"
                                         error:error]) {
                    [view discard]; return NO;
                }
            }
            TTIOGenomicRun *rd = [TTIOGenomicRun openFromGroup:view.group
                                                          name:@"block"
                                                         error:error];
            if (!rd) { [view discard]; return NO; }
            // The view is one whole block, so read [0, nn) of it.
            NSData *seq = nil, *qual = nil;
            TTIOWrittenGenomicRun *block = ttioBlocksV1BlockRun(
                rd, runGroup, study, table, b, 0, &seq, &qual, error);
            [view discard];
            if (!block) return NO;
            if (decrypted[@"tags"] != nil) block.tags = decrypted[@"tags"];
            if (refMD5 == nil && block.referenceChromSeqs != nil) {
                refMD5 = [TTIOSpectralDataset referenceMD5ForRun:block];
            }
            TTIOGenomicWriteContext *ctx = [TTIOGenomicWriteContext
                contextWithChromNameToId:chromMap referenceMD5:refMD5];
            ctx.qualStrategyHint = qualHint;
            TTIOBlockBlobs *blobs = [TTIOGenomicBlocks encodeBlock:block
                                                           context:ctx
                                                             error:error];
            if (!blobs) return NO;
            qualHint = ttioBlocksV1DeriveQualHint(blobs, qualHint);

            for (NSString *ch in channels) {
                NSData *got = blobs.blobs[ch] ?: [NSData data];
                unsigned long long off = [table offsetOf:ch at:b];
                unsigned long long ln = [table lengthOf:ch at:b];
                unsigned long long pos =
                    [written[ch] unsignedLongLongValue];
                if (got.length != (NSUInteger)ln || pos != off) {
                    mismatch = YES;
                }
                [newOff[ch] addObject:@(pos)];
                [newLen[ch] addObject:@((unsigned long long)got.length)];
                [newCodec[ch] addObject:blobs.codecs[ch] ?: @0];
                id<TTIOStorageDataset> ds = newDs[ch];
                // Like the stream writer, a channel dataset is created at
                // its first non-empty blob, so it takes that blob's codec
                // and filter (a tags channel may start after untagged
                // blocks). Empty blocks still record their range above.
                if (ds == nil && got.length == 0) continue;
                if (ds == nil) {
                    id<TTIOStorageGroup> parent;
                    NSString *dsName;
                    if ([ch isEqualToString:@"sequences"]) {
                        parent = [sig createGroupNamed:@"sequences"
                                                 error:error];
                        dsName = @"data";
                    } else {
                        parent = sig;
                        dsName = ch;
                    }
                    if (!parent) return NO;
                    NSUInteger codec =
                        [blobs.codecs[ch] unsignedIntegerValue];
                    ds = [parent createDatasetNamed:dsName
                                          precision:TTIOPrecisionUInt8
                                             length:0
                                          chunkSize:(256 << 10)
                                        compression:(codec == 0
                                            ? TTIOCompressionZlib
                                            : TTIOCompressionNone)
                                   compressionLevel:6
                                         extendable:YES
                                              error:error];
                    if (!ds) return NO;
                    if (![ds setAttributeValue:@((int64_t)codec)
                                       forName:@"compression"
                                         error:error]) return NO;
                    NSDictionary *extra = blobs.extraAttrs[ch] ?: @{};
                    for (NSString *k in [[extra allKeys]
                            sortedArrayUsingSelector:@selector(compare:)]) {
                        if (![ds setAttributeValue:extra[k]
                                           forName:k
                                             error:error]) return NO;
                    }
                    newDs[ch] = ds;
                }
                if (got.length > 0
                    && ![ds appendData:got error:error]) return NO;
                written[ch] = @(pos + got.length);
            }
        }
    }
    if (mismatch
        && !ttioBlocksV1RewriteIndex(runGroup, channels, newOff, newLen,
                                     newCodec, error)) return NO;
    for (NSString *ch in channels) {
        NSString *segName = [NSString stringWithFormat:@"%@_segments", ch];
        if (![sig deleteChildNamed:segName error:error]) return NO;
        NSString *algAttr = [NSString stringWithFormat:@"%@_algorithm", ch];
        if ([sig hasAttributeNamed:algAttr]) {
            [sig deleteAttributeNamed:algAttr error:NULL];
        }
    }
    return YES;
}


@implementation TTIOPerAUFile

+ (BOOL)encryptFilePath:(NSString *)path
                     key:(NSData *)key
         encryptHeaders:(BOOL)encryptHeaders
            providerName:(NSString *)providerName
                   error:(NSError **)error
{
    if (key.length != 32) {
        if (error) *error = makeErr(1,
            @"AES-256-GCM key must be 32 bytes, got %lu",
            (unsigned long)key.length);
        return NO;
    }

    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeReadWrite
                                                provider:providerName
                                                   error:error];
    if (!sp) return NO;
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return NO;
        NSString *version = nil;
        NSArray *featuresArr = readFeatureFlags(root, &version);
        NSMutableSet *featureSet = [NSMutableSet setWithArray:featuresArr];

        id<TTIOStorageGroup> study = [root openGroupNamed:@"study" error:error];
        if (!study) return NO;
        id<TTIOStorageGroup> msRuns = [study openGroupNamed:@"ms_runs" error:error];
        if (!msRuns) return NO;

        NSArray *runNames = listGroupChildren(msRuns);
        uint16_t datasetId = 1;
        for (NSString *runName in runNames) {
            id<TTIOStorageGroup> run = [msRuns openGroupNamed:runName error:error];
            if (!run) continue;
            id<TTIOStorageGroup> sig = [run openGroupNamed:@"signal_channels" error:error];
            id<TTIOStorageGroup> idx = [run openGroupNamed:@"spectrum_index" error:error];
            if (!sig || !idx) return NO;

            id<TTIOStorageDataset> lensDs = [idx openDatasetNamed:@"lengths" error:error];
            if (!lensDs) return NO;
            NSData *lengthsData = [lensDs readAll:error];
            if (!lengthsData) return NO;
            // offsets is omitted from disk by default;
            // synthesize from cumsum(lengths). Pre-v1.10 files have it.
            NSData *offsetsData;
            if ([idx hasChildNamed:@"offsets"]) {
                id<TTIOStorageDataset> offsDs = [idx openDatasetNamed:@"offsets" error:error];
                if (!offsDs) return NO;
                offsetsData = [offsDs readAll:error];
                if (!offsetsData) return NO;
            } else {
                offsetsData = TTIOOffsetsFromLengths(lengthsData);
            }
            NSUInteger count = lengthsData.length / 4;
            const uint64_t *offsets = (const uint64_t *)offsetsData.bytes;
            const uint32_t *lengths = (const uint32_t *)lengthsData.bytes;

            NSString *channelNamesStr = readStringAttr(sig, @"channel_names") ?: @"";
            NSArray<NSString *> *channelNames = splitChannelNames(channelNamesStr);

            for (NSString *cname in channelNames) {
                NSString *valuesName = [NSString stringWithFormat:@"%@_values", cname];
                if (![sig hasChildNamed:valuesName]) continue;
                id<TTIOStorageDataset> vDs = [sig openDatasetNamed:valuesName error:error];
                if (!vDs) return NO;
                NSData *plaintext = [vDs readAll:error];
                if (!plaintext) return NO;

                // FLOAT_DELTA_ZSTD (codec id 17, the MS default since
                // Phase 2): decode to float64 before slicing — the
                // per-AU segment contract is per-spectrum float64,
                // and decrypt writes plain float64 back.
                if ([vDs hasAttributeNamed:@"compression"]) {
                    id cval = [vDs attributeValueForName:@"compression"
                                                   error:NULL];
                    if ([cval respondsToSelector:@selector(intValue)]
                        && [cval intValue] == TTIOCompressionFloatDeltaZstd) {
                        NSData *decoded =
                            [TTIOFloatDeltaZstd decodeStream:plaintext
                                                       error:error];
                        if (!decoded) return NO;
                        plaintext = decoded;
                    }
                }

                NSArray<TTIOChannelSegment *> *segs =
                    [TTIOPerAUEncryption encryptChannelToSegments:plaintext
                                                              offsets:offsets
                                                              lengths:lengths
                                                             nSpectra:count
                                                            datasetId:datasetId
                                                          channelName:cname
                                                                  key:key
                                                                error:error];
                if (!segs) return NO;
                NSString *segName =
                    [NSString stringWithFormat:@"%@_segments", cname];
                if (!writeChannelSegments(sig, segName, segs, error)) return NO;
                if (![sig deleteChildNamed:valuesName error:error]) return NO;
                if (![sig setAttributeValue:@"aes-256-gcm"
                                     forName:[NSString stringWithFormat:@"%@_algorithm", cname]
                                       error:error]) return NO;
            }

            if (encryptHeaders) {
                int64_t acqMode = readIntAttr(run, @"acquisition_mode", 0);
                id<TTIOStorageDataset> msDs = [idx openDatasetNamed:@"ms_levels" error:error];
                id<TTIOStorageDataset> polDs = [idx openDatasetNamed:@"polarities" error:error];
                id<TTIOStorageDataset> rtDs = [idx openDatasetNamed:@"retention_times" error:error];
                id<TTIOStorageDataset> pmzDs = [idx openDatasetNamed:@"precursor_mzs" error:error];
                id<TTIOStorageDataset> pcDs = [idx openDatasetNamed:@"precursor_charges" error:error];
                id<TTIOStorageDataset> bpiDs = [idx openDatasetNamed:@"base_peak_intensities" error:error];
                if (!msDs || !polDs || !rtDs || !pmzDs || !pcDs || !bpiDs) return NO;
                NSData *msD = [msDs readAll:error];
                NSData *polD = [polDs readAll:error];
                NSData *rtD = [rtDs readAll:error];
                NSData *pmzD = [pmzDs readAll:error];
                NSData *pcD = [pcDs readAll:error];
                NSData *bpiD = [bpiDs readAll:error];
                if (!msD || !polD || !rtD || !pmzD || !pcD || !bpiD) return NO;
                const int32_t *ms = (const int32_t *)msD.bytes;
                const int32_t *pol = (const int32_t *)polD.bytes;
                const double *rt = (const double *)rtD.bytes;
                const double *pmz = (const double *)pmzD.bytes;
                const int32_t *pc = (const int32_t *)pcD.bytes;
                const double *bpi = (const double *)bpiD.bytes;

                NSMutableArray *rows = [NSMutableArray arrayWithCapacity:count];
                for (NSUInteger i = 0; i < count; i++) {
                    TTIOAUHeaderPlaintext *h = [[TTIOAUHeaderPlaintext alloc] init];
                    h.acquisitionMode = (uint8_t)(acqMode & 0xFF);
                    h.msLevel = (uint8_t)(ms[i] & 0xFF);
                    h.polarity = pol[i];
                    h.retentionTime = rt[i];
                    h.precursorMz = pmz[i];
                    h.precursorCharge = (uint8_t)(pc[i] & 0xFF);
                    h.ionMobility = 0.0;
                    h.basePeakIntensity = bpi[i];
                    [rows addObject:h];
                }
                NSArray<TTIOHeaderSegment *> *hdrSegs =
                    [TTIOPerAUEncryption encryptHeaderSegments:rows
                                                       datasetId:datasetId
                                                             key:key
                                                           error:error];
                if (!hdrSegs) return NO;
                if (!writeAUHeaderSegments(idx, @"au_header_segments",
                                             hdrSegs, error)) return NO;
                for (NSString *plainName in @[@"retention_times", @"ms_levels",
                                                @"polarities", @"precursor_mzs",
                                                @"precursor_charges",
                                                @"base_peak_intensities"]) {
                    if ([idx hasChildNamed:plainName]) {
                        if (![idx deleteChildNamed:plainName error:error]) return NO;
                    }
                }
            }
            datasetId++;
        }

        // extend encryption to genomic runs. Genomic signal
        // channels (sequences, qualities) are stored as plain uint8
        // datasets named without a "_values" suffix (different from
        // the MS layout). datasetId continues from where the MS loop
        // left off so genomic runs occupy IDs N+1..N+M (matches the
        // M89.2 transport convention and the Python reference impl).
        if ([study hasChildNamed:@"genomic_runs"]) {
            id<TTIOStorageGroup> gRuns =
                [study openGroupNamed:@"genomic_runs" error:error];
            if (!gRuns) return NO;
            NSArray *gRunNames = listGroupChildren(gRuns);
            for (NSString *gRunName in gRunNames) {
                id<TTIOStorageGroup> gRun =
                    [gRuns openGroupNamed:gRunName error:error];
                if (!gRun) continue;
                if (ttioIsBlocksV1(gRun)) {
                    if (!ttioEncryptBlocksV1Run(study, gRun, gRunName,
                                                datasetId, key, error))
                        return NO;
                    datasetId++;
                    continue;
                }
                id<TTIOStorageGroup> gSig =
                    [gRun openGroupNamed:@"signal_channels" error:error];
                if (gSig && [gSig hasChildNamed:@"tags"]) {
                    // M101: tags are encrypted with the bases; only the
                    // blocks_v1 walker carries them, and leaving them in
                    // plaintext would expose MD strings.
                    if (error) *error = makeErr(6,
                        @"per-AU encryption: genomic run '%@' uses the "
                        @"whole-channel layout and carries SAM tags; per-AU "
                        @"protection of tags needs the blocks_v1 layout. Rewrite "
                        @"the run without opt_legacy_whole_channel.", gRunName);
                    return NO;
                }
                id<TTIOStorageGroup> gIdx =
                    [gRun openGroupNamed:@"genomic_index" error:error];
                if (!gSig || !gIdx) return NO;

                id<TTIOStorageDataset> gLensDs =
                    [gIdx openDatasetNamed:@"lengths" error:error];
                if (!gLensDs) return NO;
                NSData *gLengthsData = [gLensDs readAll:error];
                if (!gLengthsData) return NO;
                // synthesize offsets from cumsum(lengths)
                // when the column is absent (default for v1.10+ files).
                NSData *gOffsetsData;
                if ([gIdx hasChildNamed:@"offsets"]) {
                    id<TTIOStorageDataset> gOffsDs =
                        [gIdx openDatasetNamed:@"offsets" error:error];
                    if (!gOffsDs) return NO;
                    gOffsetsData = [gOffsDs readAll:error];
                    if (!gOffsetsData) return NO;
                } else {
                    gOffsetsData = TTIOOffsetsFromLengths(gLengthsData);
                }
                NSUInteger gCount = gLengthsData.length / 4;
                const uint64_t *gOffsets = (const uint64_t *)gOffsetsData.bytes;
                const uint32_t *gLengths = (const uint32_t *)gLengthsData.bytes;

                for (NSString *cname in @[@"sequences", @"qualities"]) {
                    if (![gSig hasChildNamed:cname]) continue;
                    id<TTIOStorageDataset> vDs =
                        [gSig openDatasetNamed:cname error:error];
                    if (!vDs) return NO;
                    NSData *plaintext = [vDs readAll:error];
                    if (!plaintext) return NO;
                    NSArray<TTIOChannelSegment *> *segs =
                        [TTIOPerAUEncryption
                            encryptChannelToSegments:plaintext
                                              offsets:gOffsets
                                              lengths:gLengths
                                             nSpectra:gCount
                                      bytesPerElement:1
                                            datasetId:datasetId
                                          channelName:cname
                                                  key:key
                                                error:error];
                    if (!segs) return NO;
                    NSString *segName =
                        [NSString stringWithFormat:@"%@_segments", cname];
                    if (!writeChannelSegments(gSig, segName, segs, error))
                        return NO;
                    if (![gSig deleteChildNamed:cname error:error]) return NO;
                    if (![gSig setAttributeValue:@"aes-256-gcm"
                                          forName:[NSString stringWithFormat:@"%@_algorithm", cname]
                                            error:error]) return NO;
                }
                datasetId++;
            }
        }

        // M98: assembly graphs. One AU per segment record; offsets /
        // lengths come from segments/records (seq_missing rows have
        // length 0 and encrypt to empty ciphertext). The stored
        // channel is codec-encoded (@compression), so decode to raw
        // bytes before slicing; decryptFilePathInPlace writes the raw
        // channel back. datasetId continues after the genomic runs,
        // matching the Python reference impl.
        if ([study hasChildNamed:@"assembly_graphs"]) {
            id<TTIOStorageGroup> agRoot =
                [study openGroupNamed:@"assembly_graphs" error:error];
            if (!agRoot) return NO;
            for (NSString *agName in listGroupChildren(agRoot)) {
                id<TTIOStorageGroup> gGroup =
                    [agRoot openGroupNamed:agName error:error];
                id<TTIOStorageGroup> segG = nil;
                if (gGroup && [gGroup hasChildNamed:@"segments"]) {
                    segG = [gGroup openGroupNamed:@"segments" error:error];
                }
                if (!segG || ![segG hasChildNamed:@"sequences"]
                    || ![segG hasChildNamed:@"records"]) {
                    datasetId++;
                    continue;
                }
                id<TTIOStorageDataset> seqDs =
                    [segG openDatasetNamed:@"sequences" error:error];
                if (!seqDs) return NO;
                NSData *raw = decodeAssemblySequences(seqDs, error);
                if (!raw) return NO;
                id<TTIOStorageDataset> recDs =
                    [segG openDatasetNamed:@"records" error:error];
                if (!recDs) return NO;
                NSArray<NSDictionary *> *rows = [recDs readRows:error];
                if (!rows) return NO;
                NSUInteger segCount = rows.count;
                uint64_t *aOffsets = malloc(segCount * sizeof(uint64_t));
                uint32_t *aLengths = malloc(segCount * sizeof(uint32_t));
                for (NSUInteger i = 0; i < segCount; i++) {
                    aOffsets[i] =
                        [rows[i][@"seq_offset"] unsignedLongLongValue];
                    aLengths[i] =
                        (uint32_t)[rows[i][@"length"] unsignedLongLongValue];
                }
                NSArray<TTIOChannelSegment *> *segs =
                    [TTIOPerAUEncryption
                        encryptChannelToSegments:raw
                                          offsets:aOffsets
                                          lengths:aLengths
                                         nSpectra:segCount
                                  bytesPerElement:1
                                        datasetId:datasetId
                                      channelName:@"sequences"
                                              key:key
                                            error:error];
                free(aOffsets);
                free(aLengths);
                if (!segs) return NO;
                if (!writeChannelSegments(segG, @"sequences_segments",
                                          segs, error)) return NO;
                if (![segG deleteChildNamed:@"sequences" error:error])
                    return NO;
                if (![segG setAttributeValue:@"aes-256-gcm"
                                      forName:@"sequences_algorithm"
                                        error:error]) return NO;
                datasetId++;
            }
        }

        [featureSet addObject:@"opt_per_au_encryption"];
        if (encryptHeaders) [featureSet addObject:@"opt_encrypted_au_headers"];
        NSArray *sorted = [featureSet.allObjects
            sortedArrayUsingSelector:@selector(compare:)];
        if (!writeFeatureFlags(root, version, sorted, error)) return NO;
    }
    @finally {
        [sp close];
    }
    return YES;
}


+ (NSDictionary<NSString *, NSDictionary *> *)
    decryptFilePath:(NSString *)path
                key:(NSData *)key
       providerName:(NSString *)providerName
              error:(NSError **)error
{
    if (key.length != 32) {
        if (error) *error = makeErr(1,
            @"AES-256-GCM key must be 32 bytes, got %lu",
            (unsigned long)key.length);
        return nil;
    }

    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeRead
                                                provider:providerName
                                                   error:error];
    if (!sp) return nil;
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return nil;
        NSArray *features = readFeatureFlags(root, NULL);
        if (![features containsObject:@"opt_per_au_encryption"]) {
            if (error) *error = makeErr(2,
                @"%@ does not carry opt_per_au_encryption", path);
            return nil;
        }
        BOOL headersEncrypted = [features containsObject:@"opt_encrypted_au_headers"];

        // For compound-dataset reads we still need the underlying
        // TTIOHDF5Group — see readChannelSegments comment. Unwrap
        // via the provider's nativeHandle() escape hatch.
        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(3,
                @"per-AU decrypt currently requires HDF5 provider (got %@)",
                sp.providerName);
            return nil;
        }
        TTIOHDF5File *hdf5File = (TTIOHDF5File *)[sp nativeHandle];
        TTIOHDF5Group *hdf5Root = hdf5File.rootGroup;

        id<TTIOStorageGroup> study = [root openGroupNamed:@"study" error:error];
        id<TTIOStorageGroup> msRuns = [study openGroupNamed:@"ms_runs" error:error];
        if (!study || !msRuns) return nil;
        NSArray *runNames = listGroupChildren(msRuns);
        uint16_t datasetId = 1;
        for (NSString *runName in runNames) {
            id<TTIOStorageGroup> run = [msRuns openGroupNamed:runName error:error];
            id<TTIOStorageGroup> sig = [run openGroupNamed:@"signal_channels" error:error];
            id<TTIOStorageGroup> idx = [run openGroupNamed:@"spectrum_index" error:error];
            if (!run || !sig || !idx) continue;
            NSString *channelNamesStr = readStringAttr(sig, @"channel_names") ?: @"";
            NSArray<NSString *> *channelNames = splitChannelNames(channelNamesStr);

            // Raw HDF5 groups for the compound-read escape hatch.
            TTIOHDF5Group *hdf5Run =
                [[hdf5Root openGroupNamed:@"study" error:NULL]
                    openGroupNamed:@"ms_runs" error:NULL];
            TTIOHDF5Group *hdf5RunGroup =
                [hdf5Run openGroupNamed:runName error:NULL];
            TTIOHDF5Group *hdf5Sig =
                [hdf5RunGroup openGroupNamed:@"signal_channels" error:NULL];
            TTIOHDF5Group *hdf5Idx =
                [hdf5RunGroup openGroupNamed:@"spectrum_index" error:NULL];

            NSMutableDictionary *runOut = [NSMutableDictionary dictionary];
            for (NSString *cname in channelNames) {
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cname];
                if (![sig hasChildNamed:segName]) continue;
                NSArray *segs = readChannelSegments(hdf5Sig, segName, error);
                if (!segs) return nil;
                NSData *plain = [TTIOPerAUEncryption
                    decryptChannelFromSegments:segs
                                      datasetId:datasetId
                                    channelName:cname
                                            key:key
                                          error:error];
                if (!plain) return nil;
                runOut[cname] = plain;
            }
            if (headersEncrypted && [idx hasChildNamed:@"au_header_segments"]) {
                NSArray *hdrSegs = readAUHeaderSegments(hdf5Idx,
                                                          @"au_header_segments",
                                                          error);
                if (!hdrSegs) return nil;
                NSArray *rows = [TTIOPerAUEncryption
                    decryptHeaderSegments:hdrSegs
                                 datasetId:datasetId
                                       key:key
                                     error:error];
                if (!rows) return nil;
                runOut[@"__au_headers__"] = rows;
            }
            out[runName] = runOut;
            datasetId++;
        }

        // also materialise genomic_runs. datasetId continues
        // from where the MS loop left off so AAD reconstruction
        // matches the encrypt path exactly.
        if ([study hasChildNamed:@"genomic_runs"]) {
            id<TTIOStorageGroup> gRuns =
                [study openGroupNamed:@"genomic_runs" error:error];
            if (!gRuns) return nil;
            NSArray *gRunNames = listGroupChildren(gRuns);
            // Raw HDF5 access for the compound-read path.
            TTIOHDF5Group *hdf5GRuns =
                [[hdf5Root openGroupNamed:@"study" error:NULL]
                    openGroupNamed:@"genomic_runs" error:NULL];
            for (NSString *gRunName in gRunNames) {
                id<TTIOStorageGroup> gRun =
                    [gRuns openGroupNamed:gRunName error:error];
                if (!gRun) continue;
                id<TTIOStorageGroup> gSig =
                    [gRun openGroupNamed:@"signal_channels" error:error];
                if (!gSig) continue;
                TTIOHDF5Group *hdf5GRun =
                    [hdf5GRuns openGroupNamed:gRunName error:NULL];
                TTIOHDF5Group *hdf5GSig =
                    [hdf5GRun openGroupNamed:@"signal_channels" error:NULL];

                NSMutableDictionary *gRunOut = [NSMutableDictionary dictionary];
                for (NSString *cname in kBlocksV1Channels()) {
                    NSString *segName =
                        [NSString stringWithFormat:@"%@_segments", cname];
                    if (![gSig hasChildNamed:segName]) continue;
                    NSArray *segs =
                        readChannelSegments(hdf5GSig, segName, error);
                    if (!segs) return nil;
                    NSData *plain = [TTIOPerAUEncryption
                        decryptChannelFromSegments:segs
                                  bytesPerElement:1
                                          datasetId:datasetId
                                        channelName:cname
                                                key:key
                                              error:error];
                    if (!plain) return nil;
                    if ([cname isEqualToString:@"tags"]) {
                        // M101: tags come back as one string per read.
                        NSArray *texts = ttioSplitTagTexts(plain, segs, error);
                        if (!texts) return nil;
                        gRunOut[cname] = texts;
                    } else {
                        gRunOut[cname] = plain;
                    }
                }
                out[gRunName] = gRunOut;
                datasetId++;
            }
        }
    }
    @finally {
        [sp close];
    }
    return out;
}


#pragma mark - Per-AU decrypt-in-place

// Writes `data` as a new dataset of `precision` and `length` elements under
// `group`, replacing any existing child of the same name. NO-op on no data.
// Used by decryptFilePathInPlace to restore plaintext channels + headers.
static BOOL _writePlainDataset(id<TTIOStorageGroup> group, NSString *name,
                               TTIOPrecision precision, NSUInteger length,
                               NSData *data, NSError **error)
{
    if ([group hasChildNamed:name]) {
        if (![group deleteChildNamed:name error:error]) return NO;
    }
    id<TTIOStorageDataset> ds =
        [group createDatasetNamed:name
                         precision:precision
                            length:length
                         chunkSize:0
                       compression:TTIOCompressionNone
                  compressionLevel:0
                             error:error];
    if (!ds) return NO;
    return [ds writeAll:data error:error];
}

+ (BOOL)decryptFilePathInPlace:(NSString *)path
                            key:(NSData *)key
                   providerName:(NSString *)providerName
                          error:(NSError **)error
{
    if (key.length != 32) {
        if (error) *error = makeErr(1,
            @"AES-256-GCM key must be 32 bytes, got %lu",
            (unsigned long)key.length);
        return NO;
    }

    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeReadWrite
                                                provider:providerName
                                                   error:error];
    if (!sp) return NO;

    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return NO;

        NSString *version = nil;
        NSArray *features = readFeatureFlags(root, &version);
        if (![features containsObject:@"opt_per_au_encryption"]) {
            // Idempotent: file is already plaintext at the per-AU layer.
            return YES;
        }
        BOOL headersEncrypted =
            [features containsObject:@"opt_encrypted_au_headers"];

        // Need the underlying HDF5 group for compound-dataset reads, same as
        // decryptFilePath uses.
        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(3,
                @"per-AU decrypt-in-place currently requires HDF5 provider (got %@)",
                sp.providerName);
            return NO;
        }
        TTIOHDF5File *hdf5File = (TTIOHDF5File *)[sp nativeHandle];
        TTIOHDF5Group *hdf5Root = hdf5File.rootGroup;

        id<TTIOStorageGroup> study = [root openGroupNamed:@"study" error:error];
        if (!study) return NO;

        // datasetId numbering continues across MS + genomic loops, matching
        // encryptFilePath's AAD scheme (genomic IDs occupy N+1..N+M).
        uint16_t datasetId = 1;
        if ([study hasChildNamed:@"ms_runs"]) {
            id<TTIOStorageGroup> msRuns =
                [study openGroupNamed:@"ms_runs" error:error];
            if (!msRuns) return NO;
            NSArray *runNames = listGroupChildren(msRuns);
            TTIOHDF5Group *hdf5Study =
                [hdf5Root openGroupNamed:@"study" error:NULL];
            TTIOHDF5Group *hdf5MsRuns =
                [hdf5Study openGroupNamed:@"ms_runs" error:NULL];
            for (NSString *runName in runNames) {
                id<TTIOStorageGroup> run = [msRuns openGroupNamed:runName error:error];
                id<TTIOStorageGroup> sig = [run openGroupNamed:@"signal_channels" error:error];
                id<TTIOStorageGroup> idx = [run openGroupNamed:@"spectrum_index" error:error];
                if (!run || !sig || !idx) { datasetId++; continue; }

                TTIOHDF5Group *hdf5RunGroup =
                    [hdf5MsRuns openGroupNamed:runName error:NULL];
                TTIOHDF5Group *hdf5Sig =
                    [hdf5RunGroup openGroupNamed:@"signal_channels" error:NULL];
                TTIOHDF5Group *hdf5Idx =
                    [hdf5RunGroup openGroupNamed:@"spectrum_index" error:NULL];

                NSString *channelNamesStr =
                    readStringAttr(sig, @"channel_names") ?: @"";
                NSArray<NSString *> *channelNames =
                    splitChannelNames(channelNamesStr);

                for (NSString *cname in channelNames) {
                    NSString *segName =
                        [NSString stringWithFormat:@"%@_segments", cname];
                    if (![sig hasChildNamed:segName]) continue;
                    NSArray *segs =
                        readChannelSegments(hdf5Sig, segName, error);
                    if (!segs) return NO;
                    NSData *plain = [TTIOPerAUEncryption
                        decryptChannelFromSegments:segs
                                          datasetId:datasetId
                                        channelName:cname
                                                key:key
                                              error:error];
                    if (!plain) return NO;
                    NSString *valuesName =
                        [NSString stringWithFormat:@"%@_values", cname];
                    NSUInteger nValues = plain.length / 8;  // float64
                    if (!_writePlainDataset(sig, valuesName,
                                            TTIOPrecisionFloat64,
                                            nValues, plain, error)) return NO;
                    if (![sig deleteChildNamed:segName error:error]) return NO;
                    NSString *algAttr =
                        [NSString stringWithFormat:@"%@_algorithm", cname];
                    if ([sig hasAttributeNamed:algAttr]) {
                        [sig deleteAttributeNamed:algAttr error:NULL];
                    }
                }

                if (headersEncrypted && [idx hasChildNamed:@"au_header_segments"]) {
                    NSArray *hdrSegs = readAUHeaderSegments(hdf5Idx,
                                                              @"au_header_segments",
                                                              error);
                    if (!hdrSegs) return NO;
                    NSArray<TTIOAUHeaderPlaintext *> *rows =
                        [TTIOPerAUEncryption
                            decryptHeaderSegments:hdrSegs
                                         datasetId:datasetId
                                               key:key
                                             error:error];
                    if (!rows) return NO;
                    NSUInteger n = rows.count;
                    NSMutableData *msD  = [NSMutableData dataWithLength:n * 4];
                    NSMutableData *polD = [NSMutableData dataWithLength:n * 4];
                    NSMutableData *pcD  = [NSMutableData dataWithLength:n * 4];
                    NSMutableData *rtD  = [NSMutableData dataWithLength:n * 8];
                    NSMutableData *pmzD = [NSMutableData dataWithLength:n * 8];
                    NSMutableData *bpiD = [NSMutableData dataWithLength:n * 8];
                    int32_t *ms  = msD.mutableBytes,
                            *pol = polD.mutableBytes,
                            *pc  = pcD.mutableBytes;
                    double  *rt  = rtD.mutableBytes,
                            *pmz = pmzD.mutableBytes,
                            *bpi = bpiD.mutableBytes;
                    for (NSUInteger i = 0; i < n; i++) {
                        TTIOAUHeaderPlaintext *h = rows[i];
                        ms[i]  = (int32_t)h.msLevel;
                        pol[i] = (int32_t)h.polarity;
                        pc[i]  = (int32_t)h.precursorCharge;
                        rt[i]  = h.retentionTime;
                        pmz[i] = h.precursorMz;
                        bpi[i] = h.basePeakIntensity;
                    }
                    if (!_writePlainDataset(idx, @"ms_levels",
                                            TTIOPrecisionInt32, n, msD, error)
                        || !_writePlainDataset(idx, @"polarities",
                                               TTIOPrecisionInt32, n, polD, error)
                        || !_writePlainDataset(idx, @"precursor_charges",
                                               TTIOPrecisionInt32, n, pcD, error)
                        || !_writePlainDataset(idx, @"retention_times",
                                               TTIOPrecisionFloat64, n, rtD, error)
                        || !_writePlainDataset(idx, @"precursor_mzs",
                                               TTIOPrecisionFloat64, n, pmzD, error)
                        || !_writePlainDataset(idx, @"base_peak_intensities",
                                               TTIOPrecisionFloat64, n, bpiD, error)) {
                        return NO;
                    }
                    if (![idx deleteChildNamed:@"au_header_segments" error:error])
                        return NO;
                }

                datasetId++;
            }
        }

        // Genomic runs: same per-AU GCM scheme as MS, but channel data
        // (sequences / qualities) is uint8 stored under the bare channel
        // name (no _values suffix). datasetId continues from where the MS
        // loop left off, matching encryptFilePath's AAD numbering.
        if ([study hasChildNamed:@"genomic_runs"]) {
            id<TTIOStorageGroup> gRuns =
                [study openGroupNamed:@"genomic_runs" error:error];
            if (!gRuns) return NO;
            NSArray *gRunNames = listGroupChildren(gRuns);
            TTIOHDF5Group *hdf5Study =
                [hdf5Root openGroupNamed:@"study" error:NULL];
            TTIOHDF5Group *hdf5GRuns =
                [hdf5Study openGroupNamed:@"genomic_runs" error:NULL];
            for (NSString *gRunName in gRunNames) {
                id<TTIOStorageGroup> gRun =
                    [gRuns openGroupNamed:gRunName error:error];
                id<TTIOStorageGroup> gSig =
                    [gRun openGroupNamed:@"signal_channels" error:error];
                if (!gRun || !gSig) { datasetId++; continue; }
                if (ttioIsBlocksV1(gRun)) {
                    if (!ttioDecryptBlocksV1RunInPlace(study, gRun,
                                                       datasetId, key,
                                                       error))
                        return NO;
                    datasetId++;
                    continue;
                }
                TTIOHDF5Group *hdf5GRun =
                    [hdf5GRuns openGroupNamed:gRunName error:NULL];
                TTIOHDF5Group *hdf5GSig =
                    [hdf5GRun openGroupNamed:@"signal_channels" error:NULL];

                for (NSString *cname in @[@"sequences", @"qualities"]) {
                    NSString *segName =
                        [NSString stringWithFormat:@"%@_segments", cname];
                    if (![gSig hasChildNamed:segName]) continue;
                    NSArray *segs =
                        readChannelSegments(hdf5GSig, segName, error);
                    if (!segs) return NO;
                    NSData *plain = [TTIOPerAUEncryption
                        decryptChannelFromSegments:segs
                                   bytesPerElement:1
                                          datasetId:datasetId
                                        channelName:cname
                                                key:key
                                              error:error];
                    if (!plain) return NO;
                    // Write back as uint8 dataset under the bare channel name
                    // (genomic layout has no _values suffix; see
                    // encryptFilePath's genomic loop).
                    if (!_writePlainDataset(gSig, cname,
                                            TTIOPrecisionUInt8,
                                            plain.length, plain, error)) return NO;
                    if (![gSig deleteChildNamed:segName error:error]) return NO;
                    NSString *algAttr =
                        [NSString stringWithFormat:@"%@_algorithm", cname];
                    if ([gSig hasAttributeNamed:algAttr]) {
                        [gSig deleteAttributeNamed:algAttr error:NULL];
                    }
                }
                datasetId++;
            }
        }

        // M98: assembly graphs. The raw sequences bytes come back as
        // a plain uint8 dataset with no @compression; the graph
        // re-emits byte-exactly from raw bytes. datasetId numbering
        // mirrors encryptFilePath exactly.
        if ([study hasChildNamed:@"assembly_graphs"]) {
            id<TTIOStorageGroup> agRoot =
                [study openGroupNamed:@"assembly_graphs" error:error];
            if (!agRoot) return NO;
            TTIOHDF5Group *hdf5StudyAg =
                [hdf5Root openGroupNamed:@"study" error:NULL];
            TTIOHDF5Group *hdf5AgRoot =
                [hdf5StudyAg openGroupNamed:@"assembly_graphs" error:NULL];
            for (NSString *agName in listGroupChildren(agRoot)) {
                id<TTIOStorageGroup> gGroup =
                    [agRoot openGroupNamed:agName error:error];
                id<TTIOStorageGroup> segG = nil;
                if (gGroup && [gGroup hasChildNamed:@"segments"]) {
                    segG = [gGroup openGroupNamed:@"segments" error:error];
                }
                if (!segG || ![segG hasChildNamed:@"sequences_segments"]) {
                    datasetId++;
                    continue;
                }
                TTIOHDF5Group *hdf5G =
                    [hdf5AgRoot openGroupNamed:agName error:NULL];
                TTIOHDF5Group *hdf5Seg =
                    [hdf5G openGroupNamed:@"segments" error:NULL];
                NSArray *segs =
                    readChannelSegments(hdf5Seg, @"sequences_segments",
                                        error);
                if (!segs) return NO;
                NSData *plain = [TTIOPerAUEncryption
                    decryptChannelFromSegments:segs
                               bytesPerElement:1
                                     datasetId:datasetId
                                   channelName:@"sequences"
                                           key:key
                                         error:error];
                if (!plain) return NO;
                if (!_writePlainDataset(segG, @"sequences",
                                        TTIOPrecisionUInt8,
                                        plain.length, plain, error))
                    return NO;
                if (![segG deleteChildNamed:@"sequences_segments"
                                      error:error]) return NO;
                if ([segG hasAttributeNamed:@"sequences_algorithm"]) {
                    [segG deleteAttributeNamed:@"sequences_algorithm"
                                         error:NULL];
                }
                datasetId++;
            }
        }

        // Strip the per-AU feature flags + the root @encrypted attribute now
        // that all encrypted segments (MS + genomic) have been decrypted.
        NSMutableSet *featureSet = [NSMutableSet setWithArray:features];
        [featureSet removeObject:@"opt_per_au_encryption"];
        [featureSet removeObject:@"opt_encrypted_au_headers"];
        NSArray *sorted = [featureSet.allObjects
            sortedArrayUsingSelector:@selector(compare:)];
        if (!writeFeatureFlags(root, version, sorted, error)) return NO;

        if ([root hasAttributeNamed:@"encrypted"]) {
            [root deleteAttributeNamed:@"encrypted" error:NULL];
        }
    }
    @finally {
        [sp close];
    }
    return YES;
}


#pragma mark - M90.4 — region-based per-AU encryption

+ (BOOL)encryptFilePathByRegion:(NSString *)path
                          keyMap:(NSDictionary<NSString *, NSData *> *)keyMap
                    providerName:(NSString *)providerName
                           error:(NSError **)error
{
    for (NSString *chrom in keyMap) {
        NSData *k = keyMap[chrom];
        if (k.length != 32) {
            if (error) *error = makeErr(1,
                @"AES-256-GCM key for chromosome %@ must be 32 bytes, got %lu",
                chrom, (unsigned long)k.length);
            return NO;
        }
    }

    // split out the reserved "_headers" entry from the
    // chromosome-keyed entries. Headers-key encrypts the four
    // genomic_index columns (chromosomes, positions,
    // mapping_qualities, flags). Chromosome keys encrypt signal
    // channels per-AU.
    NSData *headersKey = keyMap[kTTIOPerAUHeadersKeyName];
    NSMutableDictionary<NSString *, NSData *> *chromosomeKeys =
        [keyMap mutableCopy];
    [chromosomeKeys removeObjectForKey:kTTIOPerAUHeadersKeyName];

    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeReadWrite
                                                provider:providerName
                                                   error:error];
    if (!sp) return NO;
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return NO;
        NSString *version = nil;
        NSArray *featuresArr = readFeatureFlags(root, &version);
        NSMutableSet *featureSet = [NSMutableSet setWithArray:featuresArr];

        id<TTIOStorageGroup> study =
            [root openGroupNamed:@"study" error:error];
        if (!study) return NO;
        if (![study hasChildNamed:@"genomic_runs"]) {
            // No genomic data — nothing to encrypt.
            return YES;
        }

        // Match the dataset_id_counter convention from the MS path:
        // MS runs occupy 1..N, genomic N+1..N+M. Region-only
        // encryption walks MS first to *count* runs (without
        // touching them) so genomic AAD reconstruction is correct.
        NSUInteger nMs = 0;
        if ([study hasChildNamed:@"ms_runs"]) {
            id<TTIOStorageGroup> msRuns =
                [study openGroupNamed:@"ms_runs" error:error];
            if (!msRuns) return NO;
            nMs = listGroupChildren(msRuns).count;
        }
        uint16_t datasetId = (uint16_t)(nMs + 1);

        // We need raw HDF5 access for the chromosomes compound read.
        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(3,
                @"per-AU region encrypt currently requires HDF5 provider "
                @"(got %@)", sp.providerName);
            return NO;
        }
        TTIOHDF5File *hdf5File = (TTIOHDF5File *)[sp nativeHandle];
        TTIOHDF5Group *hdf5Root = hdf5File.rootGroup;
        TTIOHDF5Group *hdf5Study =
            [hdf5Root openGroupNamed:@"study" error:NULL];
        TTIOHDF5Group *hdf5GRuns =
            [hdf5Study openGroupNamed:@"genomic_runs" error:NULL];

        id<TTIOStorageGroup> gRuns =
            [study openGroupNamed:@"genomic_runs" error:error];
        if (!gRuns) return NO;
        NSArray *gRunNames = listGroupChildren(gRuns);
        for (NSString *gRunName in gRunNames) {
            id<TTIOStorageGroup> gRun =
                [gRuns openGroupNamed:gRunName error:error];
            if (!gRun) continue;
            id<TTIOStorageGroup> gSig =
                [gRun openGroupNamed:@"signal_channels" error:error];
            if (gSig && [gSig hasChildNamed:@"tags"]) {
                // M101: region keys do not cover the SAM tags; refuse
                // rather than leave MD strings in plaintext.
                if (error) *error = makeErr(6,
                    @"per-AU region encryption: genomic run '%@' carries "
                    @"SAM tags, which region keys do not protect; use encrypt_per_au",
                    gRunName);
                return NO;
            }
            id<TTIOStorageGroup> gIdx =
                [gRun openGroupNamed:@"genomic_index" error:error];
            if (!gSig || !gIdx) return NO;

            TTIOHDF5Group *hdf5GRun =
                [hdf5GRuns openGroupNamed:gRunName error:NULL];
            TTIOHDF5Group *hdf5GIdx =
                [hdf5GRun openGroupNamed:@"genomic_index" error:NULL];

            id<TTIOStorageDataset> gLensDs =
                [gIdx openDatasetNamed:@"lengths" error:error];
            if (!gLensDs) return NO;
            NSData *gLengthsData = [gLensDs readAll:error];
            if (!gLengthsData) return NO;
            // synthesize offsets when absent.
            NSData *gOffsetsData;
            if ([gIdx hasChildNamed:@"offsets"]) {
                id<TTIOStorageDataset> gOffsDs =
                    [gIdx openDatasetNamed:@"offsets" error:error];
                if (!gOffsDs) return NO;
                gOffsetsData = [gOffsDs readAll:error];
                if (!gOffsetsData) return NO;
            } else {
                gOffsetsData = TTIOOffsetsFromLengths(gLengthsData);
            }
            NSUInteger gCount = gLengthsData.length / 4;
            const uint64_t *gOffsets = (const uint64_t *)gOffsetsData.bytes;
            const uint32_t *gLengths = (const uint32_t *)gLengthsData.bytes;

            NSArray<NSString *> *chromosomes =
                readChromosomes(hdf5GIdx, error);
            if (!chromosomes) return NO;

            // signal-channel encryption runs in two cases:
            //   (a) caller supplied chromosome keys (M90.4 path)
            //   (b) caller supplied an empty key_map (M90.4 no-op:
            //       file gets opt_per_au_encryption with all-clear
            //       segments)
            // The only path that SKIPS signal-channel encryption is
            // the headers-only case (key_map == {"_headers": K}).
            BOOL runSignalEncrypt =
                (chromosomeKeys.count > 0) || (headersKey == nil);
            if (runSignalEncrypt) {
                for (NSString *cname in @[@"sequences", @"qualities"]) {
                    if (![gSig hasChildNamed:cname]) continue;
                    id<TTIOStorageDataset> vDs =
                        [gSig openDatasetNamed:cname error:error];
                    if (!vDs) return NO;
                    NSData *plaintext = [vDs readAll:error];
                    if (!plaintext) return NO;
                    NSArray<TTIOChannelSegment *> *segs =
                        encryptChannelWithDispatch(plaintext,
                                                      gOffsets, gLengths,
                                                      chromosomes, gCount,
                                                      datasetId, cname,
                                                      chromosomeKeys, error);
                    if (!segs) return NO;
                    NSString *segName =
                        [NSString stringWithFormat:@"%@_segments", cname];
                    if (!writeChannelSegments(gSig, segName, segs, error))
                        return NO;
                    if (![gSig deleteChildNamed:cname error:error]) return NO;
                    if (![gSig setAttributeValue:@"aes-256-gcm-by-region"
                                          forName:[NSString stringWithFormat:@"%@_algorithm", cname]
                                            error:error]) return NO;
                }
            }

            // when "_headers" key is present, encrypt the
            // four genomic_index columns under it.
            if (headersKey != nil) {
                if (!encryptGenomicIndex(gIdx, datasetId, headersKey,
                                            chromosomes, error)) {
                    return NO;
                }
            }
            datasetId++;
        }

        // Feature-flag set rules:
        //  * opt_per_au_encryption — set whenever signal-channel
        //    encryption ran (chromosome keys present OR empty key_map
        //    no-op path) OR when headers_key is provided.
        //  * opt_region_keyed_encryption — only when at least one
        //    chromosome key was provided. Empty key_map leaves the
        //    file with all-clear segments — that's M90.4's no-op
        //    semantics, not a region-keyed file.
        //  * opt_encrypted_au_headers — set when _headers key was
        //    used (M90.11).
        if (chromosomeKeys.count > 0 || headersKey == nil) {
            [featureSet addObject:@"opt_per_au_encryption"];
        }
        if (chromosomeKeys.count > 0) {
            [featureSet addObject:@"opt_region_keyed_encryption"];
        }
        if (headersKey != nil) {
            [featureSet addObject:@"opt_per_au_encryption"];
            [featureSet addObject:@"opt_encrypted_au_headers"];
        }
        NSArray *sorted = [featureSet.allObjects
            sortedArrayUsingSelector:@selector(compare:)];
        if (!writeFeatureFlags(root, version, sorted, error)) return NO;
    }
    @finally {
        [sp close];
    }
    return YES;
}


+ (NSDictionary<NSString *, NSDictionary *> *)
    decryptFilePathByRegion:(NSString *)path
                      keyMap:(NSDictionary<NSString *, NSData *> *)keyMap
                providerName:(NSString *)providerName
                       error:(NSError **)error
{
    for (NSString *chrom in keyMap) {
        NSData *k = keyMap[chrom];
        if (k.length != 32) {
            if (error) *error = makeErr(1,
                @"AES-256-GCM key for chromosome %@ must be 32 bytes, got %lu",
                chrom, (unsigned long)k.length);
            return nil;
        }
    }

    // split out the reserved "_headers" entry from the
    // chromosome keys.
    NSData *headersKey = keyMap[kTTIOPerAUHeadersKeyName];
    NSMutableDictionary<NSString *, NSData *> *chromosomeKeys =
        [keyMap mutableCopy];
    [chromosomeKeys removeObjectForKey:kTTIOPerAUHeadersKeyName];

    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeRead
                                                provider:providerName
                                                   error:error];
    if (!sp) return nil;
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return nil;
        NSArray *features = readFeatureFlags(root, NULL);
        if (![features containsObject:@"opt_per_au_encryption"]) {
            if (error) *error = makeErr(2,
                @"%@ does not carry opt_per_au_encryption", path);
            return nil;
        }

        // opt_encrypted_au_headers files require the
        // "_headers" key — without it we can't even reconstruct the
        // chromosomes column needed to dispatch signal-channel
        // decryption.
        BOOL headersEncrypted =
            [features containsObject:@"opt_encrypted_au_headers"];
        if (headersEncrypted && headersKey == nil) {
            if (error) *error = makeErr(4,
                @"%@ carries opt_encrypted_au_headers; caller must "
                @"provide a '_headers' entry in keyMap to decrypt the "
                @"genomic_index columns", path);
            return nil;
        }

        id<TTIOStorageGroup> study =
            [root openGroupNamed:@"study" error:error];
        if (!study) return nil;
        if (![study hasChildNamed:@"genomic_runs"]) {
            return out;
        }

        // Walk MS runs first to keep dataset_id aligned. Region
        // decrypt does not touch MS runs; counter is purely for AAD.
        NSUInteger nMs = 0;
        if ([study hasChildNamed:@"ms_runs"]) {
            id<TTIOStorageGroup> msRuns =
                [study openGroupNamed:@"ms_runs" error:error];
            if (!msRuns) return nil;
            nMs = listGroupChildren(msRuns).count;
        }
        uint16_t datasetId = (uint16_t)(nMs + 1);

        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(3,
                @"per-AU region decrypt currently requires HDF5 provider "
                @"(got %@)", sp.providerName);
            return nil;
        }
        TTIOHDF5File *hdf5File = (TTIOHDF5File *)[sp nativeHandle];
        TTIOHDF5Group *hdf5Root = hdf5File.rootGroup;
        TTIOHDF5Group *hdf5Study =
            [hdf5Root openGroupNamed:@"study" error:NULL];
        TTIOHDF5Group *hdf5GRuns =
            [hdf5Study openGroupNamed:@"genomic_runs" error:NULL];

        id<TTIOStorageGroup> gRuns =
            [study openGroupNamed:@"genomic_runs" error:error];
        if (!gRuns) return nil;
        NSArray *gRunNames = listGroupChildren(gRuns);
        for (NSString *gRunName in gRunNames) {
            id<TTIOStorageGroup> gRun =
                [gRuns openGroupNamed:gRunName error:error];
            id<TTIOStorageGroup> gSig =
                [gRun openGroupNamed:@"signal_channels" error:error];
            id<TTIOStorageGroup> gIdx =
                [gRun openGroupNamed:@"genomic_index" error:error];
            if (!gRun || !gSig || !gIdx) continue;
            TTIOHDF5Group *hdf5GRun =
                [hdf5GRuns openGroupNamed:gRunName error:NULL];
            TTIOHDF5Group *hdf5GSig =
                [hdf5GRun openGroupNamed:@"signal_channels" error:NULL];
            TTIOHDF5Group *hdf5GIdx =
                [hdf5GRun openGroupNamed:@"genomic_index" error:NULL];

            NSMutableDictionary *gRunOut = [NSMutableDictionary dictionary];
            NSArray<NSString *> *chromosomes = nil;

            // when headers are encrypted, decrypt the four
            // genomic_index columns FIRST so the per-AU signal-channel
            // dispatch (which needs chromosomes) can proceed.
            if (headersEncrypted) {
                NSDictionary *indexPlain =
                    decryptGenomicIndex(gIdx, datasetId, headersKey, error);
                if (!indexPlain) return nil;
                chromosomes = (NSArray<NSString *> *)indexPlain[@"chromosomes"];
                gRunOut[@"__index__"] = indexPlain;
            } else {
                chromosomes = readChromosomes(hdf5GIdx, error);
                if (!chromosomes) return nil;
            }

            for (NSString *cname in @[@"sequences", @"qualities"]) {
                NSString *segName =
                    [NSString stringWithFormat:@"%@_segments", cname];
                if (![gSig hasChildNamed:segName]) continue;
                NSArray *segs = readChannelSegments(hdf5GSig, segName, error);
                if (!segs) return nil;
                NSData *plain = decryptChannelWithDispatch(segs,
                                                            chromosomes,
                                                            datasetId, cname,
                                                            chromosomeKeys, error);
                if (!plain) return nil;
                gRunOut[cname] = plain;
            }
            out[gRunName] = gRunOut;
            datasetId++;
        }
    }
    @finally {
        [sp close];
    }
    return out;
}

@end
