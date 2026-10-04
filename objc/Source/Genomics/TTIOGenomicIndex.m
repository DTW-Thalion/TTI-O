/*
 * TTIOGenomicIndex.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOGenomicIndex
 * Inherits From: NSObject
 * Declared In:   Genomics/TTIOGenomicIndex.h
 *
 * Per-read offsets, lengths, positions, mapping qualities, flags,
 * and chromosome strings for one genomic run. Loaded eagerly at
 * open time; range queries (region, unmapped, flag mask) operate
 * on the in-memory parallel arrays without touching the heavy
 * signal channels.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOGenomicIndex.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOCompoundField.h"
#import "HDF5/TTIOHDF5Types.h"  // TTIOPrecisionElementSize()
#import "HDF5/TTIOHDF5Group.h"
#import "Dataset/TTIOCompoundIO.h"

/* v1.10 #10 helper: synthesize per-record byte offsets from a uint32
 * lengths array. offsets[i] = sum(lengths[0..i]), produced as a uint64
 * NSData blob to avoid the >4 GB overflow cliff on deep WGS even when
 * the input lengths are uint32. Empty input returns 0-byte NSData. */
NSData *TTIOOffsetsFromLengths(NSData *lengths)
{
    NSUInteger n = lengths.length / sizeof(uint32_t);
    if (n == 0) {
        return [NSData data];
    }
    const uint32_t *lenBuf = (const uint32_t *)lengths.bytes;
    NSMutableData *out = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    uint64_t *outBuf = (uint64_t *)out.mutableBytes;
    outBuf[0] = 0;
    uint64_t acc = 0;
    for (NSUInteger i = 1; i < n; i++) {
        acc += (uint64_t)lenBuf[i - 1];
        outBuf[i] = acc;
    }
    return out;
}

@implementation TTIOGenomicIndex {
    NSData *_offsetsData;
    NSData *_lengthsData;
    NSArray<NSString *> *_chromosomes;
    NSData *_positionsData;
    NSData *_mappingQualitiesData;
    NSData *_flagsData;

    // PO1 (P1.2): interned chromosome table, populated by
    // +readFromGroup: from the on-disk uint16 id column (Task #82).
    // Both nil for in-memory construction — region lookup then falls
    // back to scanning the _chromosomes string array. When present,
    // -indicesForRegion: resolves the query name to its id ONCE and
    // scans a uint16 comparison instead of an O(N) per-read
    // -isEqualToString:. Mirrors Java GenomicIndex (PJ1) +
    // Python genomic_index.py:126-146.
    NSData *_chromosomeIdsData;                            // uint16 per read; nil in-memory
    NSDictionary<NSString *, NSNumber *> *_chromosomeNameToId; // nil when no ids
}

- (instancetype)initWithOffsets:(NSData *)offsets
                         lengths:(NSData *)lengths
                     chromosomes:(NSArray<NSString *> *)chromosomes
                       positions:(NSData *)positions
                mappingQualities:(NSData *)mappingQualities
                            flags:(NSData *)flags
{
    self = [super init];
    if (self) {
        _offsetsData          = [offsets copy];
        _lengthsData          = [lengths copy];
        _chromosomes          = [chromosomes copy];
        _positionsData        = [positions copy];
        _mappingQualitiesData = [mappingQualities copy];
        _flagsData            = [flags copy];
    }
    return self;
}

- (NSUInteger)count
{
    return _offsetsData.length / sizeof(uint64_t);
}

- (uint64_t)offsetAt:(NSUInteger)index
{
    return ((const uint64_t *)_offsetsData.bytes)[index];
}

- (uint32_t)lengthAt:(NSUInteger)index
{
    return ((const uint32_t *)_lengthsData.bytes)[index];
}

- (int64_t)positionAt:(NSUInteger)index
{
    return ((const int64_t *)_positionsData.bytes)[index];
}

- (uint8_t)mappingQualityAt:(NSUInteger)index
{
    return ((const uint8_t *)_mappingQualitiesData.bytes)[index];
}

- (uint32_t)flagsAt:(NSUInteger)index
{
    return ((const uint32_t *)_flagsData.bytes)[index];
}

- (NSString *)chromosomeAt:(NSUInteger)index
{
    return _chromosomes[index];
}

- (NSIndexSet *)indicesForRegion:(NSString *)chromosome
                            start:(int64_t)start
                              end:(int64_t)end
{
    NSUInteger n = self.count;
    const int64_t *positions = (const int64_t *)_positionsData.bytes;

    // PO1 fast path: on a disk-loaded index the interned uint16 ids are
    // retained alongside a name→id map, so resolve the query chromosome
    // to its id ONCE and scan a uint16 comparison per read. An absent
    // chromosome resolves to no id and yields the empty set — identical
    // to the never-matched behaviour of the string scan. Results
    // (indices + ascending order) are identical to the fallback below.
    if (_chromosomeIdsData != nil) {
        NSNumber *tid = _chromosomeNameToId[chromosome];
        if (tid == nil) return [NSIndexSet indexSet];
        uint16_t target = tid.unsignedShortValue;
        NSMutableIndexSet *result = [NSMutableIndexSet indexSet];
        const uint16_t *ids = (const uint16_t *)_chromosomeIdsData.bytes;
        for (NSUInteger i = 0; i < n; i++) {
            if (ids[i] == target
                && positions[i] >= start
                && positions[i] < end) {
                [result addIndex:i];
            }
        }
        return result;
    }

    // Fallback: in-memory index without interned ids (verbatim old loop).
    NSMutableIndexSet *result = [NSMutableIndexSet indexSet];
    for (NSUInteger i = 0; i < n; i++) {
        if ([_chromosomes[i] isEqualToString:chromosome]
            && positions[i] >= start
            && positions[i] < end) {
            [result addIndex:i];
        }
    }
    return result;
}

- (NSIndexSet *)indicesForUnmapped
{
    return [self indicesForFlag:0x4];
}

- (NSIndexSet *)indicesForFlag:(uint32_t)flagMask
{
    NSMutableIndexSet *result = [NSMutableIndexSet indexSet];
    NSUInteger n = self.count;
    const uint32_t *flags = (const uint32_t *)_flagsData.bytes;
    for (NSUInteger i = 0; i < n; i++) {
        if ((flags[i] & flagMask) != 0) {
            [result addIndex:i];
        }
    }
    return result;
}

// ── Disk I/O via the provider-agnostic StorageGroup protocol ───────

static BOOL writeTypedChannel(id<TTIOStorageGroup> g, NSString *name,
                              TTIOPrecision p, NSData *data, NSError **error)
{
    NSUInteger n = data.length / TTIOPrecisionElementSize(p);
    id<TTIOStorageDataset> ds = [g createDatasetNamed:name
                                             precision:p
                                                length:n
                                             chunkSize:65536
                                           compression:TTIOCompressionZlib
                                      compressionLevel:6
                                                 error:error];
    if (!ds) return NO;
    return [ds writeAll:data error:error];
}

static NSData *readTypedChannel(id<TTIOStorageGroup> g, NSString *name,
                                NSError **error)
{
    id<TTIOStorageDataset> ds = [g openDatasetNamed:name error:error];
    if (!ds) return nil;
    id val = [ds readAll:error];
    return [val isKindOfClass:[NSData class]] ? val : nil;
}

- (BOOL)writeToGroup:(id<TTIOStorageGroup>)group error:(NSError **)error
{
    return [self writeToGroup:group nameToId:nil error:error];
}

+ (NSArray<NSString *> *)namesInIdOrder:(NSDictionary<NSString *, NSNumber *> *)nameToId
{
    return [nameToId keysSortedByValueUsingComparator:^NSComparisonResult(NSNumber *a, NSNumber *b) {
        return [a compare:b];
    }];
}

- (BOOL)writeToGroup:(id<TTIOStorageGroup>)group
             nameToId:(NSMutableDictionary<NSString *, NSNumber *> *)sharedNameToId
                error:(NSError **)error
{
    // v1.10 #10 (offsets-cumsum): the redundant ``offsets`` column is
    // omitted on disk — readers derive it from cumsum(lengths).
    if (!writeTypedChannel(group, @"lengths",           TTIOPrecisionUInt32, _lengthsData,          error)) return NO;
    if (!writeTypedChannel(group, @"positions",         TTIOPrecisionInt64,  _positionsData,        error)) return NO;
    if (!writeTypedChannel(group, @"mapping_qualities", TTIOPrecisionUInt8,  _mappingQualitiesData, error)) return NO;
    if (!writeTypedChannel(group, @"flags",             TTIOPrecisionUInt32, _flagsData,            error)) return NO;

    // L1 (Task #82 Phase B.1, 2026-05-01): chromosomes are stored as
    // chromosome_ids (uint16) + chromosome_names (compound) instead
    // of a single VL-string compound. The old layout cost 42 MB of
    // HDF5 fractal-heap overhead per chr22 file (one heap block per
    // chunk × 432 chunks). Encounter-order id assignment —
    // first occurrence gets the next unused id.
    NSMutableDictionary<NSString *, NSNumber *> *nameToId =
        sharedNameToId ?: [NSMutableDictionary dictionary];
    NSMutableData *idsData = [NSMutableData dataWithLength:_chromosomes.count * sizeof(uint16_t)];
    uint16_t *idsBuf = (uint16_t *)idsData.mutableBytes;
    NSUInteger idx = 0;
    for (NSString *name in _chromosomes) {
        NSNumber *slot = nameToId[name];
        if (!slot) {
            if (nameToId.count > 65535) {
                if (error) *error = [NSError errorWithDomain:@"TTIOGenomicIndex"
                                                         code:1
                                                     userInfo:@{NSLocalizedDescriptionKey: @"> 65,535 unique chromosome names; uint16 chromosome_ids would overflow"}];
                return NO;
            }
            slot = @(nameToId.count);
            nameToId[name] = slot;
        }
        idsBuf[idx++] = (uint16_t)slot.unsignedShortValue;
    }
    NSArray<NSString *> *namesInOrder = [TTIOGenomicIndex namesInIdOrder:nameToId];
    if (!writeTypedChannel(group, @"chromosome_ids", TTIOPrecisionUInt16, idsData, error)) return NO;

    NSArray *nameFields = @[[TTIOCompoundField fieldWithName:@"name"
                                                        kind:TTIOCompoundFieldKindVLString]];
    NSMutableArray *nameRows = [NSMutableArray arrayWithCapacity:namesInOrder.count];
    for (NSString *n in namesInOrder) {
        [nameRows addObject:@{@"name": n}];
    }
    if ([group respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *h5 = [(id)group performSelector:@selector(unwrap)];
        return [TTIOCompoundIO writeGeneric:nameRows
                                   intoGroup:h5
                                datasetNamed:@"chromosome_names"
                                      fields:nameFields
                                       error:error];
    }
    id<TTIOStorageDataset> ds = [group createCompoundDatasetNamed:@"chromosome_names"
                                                            fields:nameFields
                                                             count:namesInOrder.count
                                                             error:error];
    if (!ds) return NO;
    return [ds writeAll:nameRows error:error];
}

+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)group error:(NSError **)error
{
    NSError *cerr = nil;
    NSData *lengths   = readTypedChannel(group, @"lengths",           &cerr);
    if (!lengths)   { if (error) *error = cerr; return nil; }
    // offsets is omitted from disk by default — synthesize
    // from cumsum(lengths). Pre-v1.10 files have it on disk.
    NSData *offsets;
    if ([group hasChildNamed:@"offsets"]) {
        offsets = readTypedChannel(group, @"offsets",           &cerr);
        if (!offsets)   { if (error) *error = cerr; return nil; }
    } else {
        offsets = TTIOOffsetsFromLengths(lengths);
    }
    NSData *positions = readTypedChannel(group, @"positions",         &cerr);
    if (!positions) { if (error) *error = cerr; return nil; }
    NSData *mapqs     = readTypedChannel(group, @"mapping_qualities", &cerr);
    if (!mapqs)     { if (error) *error = cerr; return nil; }
    NSData *flags     = readTypedChannel(group, @"flags",             &cerr);
    if (!flags)     { if (error) *error = cerr; return nil; }

    // L1 (Task #82 Phase B.1): read chromosome_ids (uint16) +
    // chromosome_names (compound) instead of a single VL-string
    // compound; materialise back to NSArray<NSString *> for callers.
    NSData *idsData = readTypedChannel(group, @"chromosome_ids", &cerr);
    if (!idsData) { if (error) *error = cerr; return nil; }
    NSArray<NSDictionary *> *nameRows = nil;
    if ([group respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *h5 = [(id)group performSelector:@selector(unwrap)];
        NSArray *fields = @[[TTIOCompoundField fieldWithName:@"name"
                                                         kind:TTIOCompoundFieldKindVLString]];
        nameRows = [TTIOCompoundIO readGenericFromGroup:h5
                                            datasetNamed:@"chromosome_names"
                                                  fields:fields
                                                   error:&cerr];
    } else {
        id<TTIOStorageDataset> nameDs = [group openDatasetNamed:@"chromosome_names" error:&cerr];
        if (nameDs) nameRows = [nameDs readAll:&cerr];
    }
    if (!nameRows) { if (error) *error = cerr; return nil; }
    NSMutableArray<NSString *> *nameTable = [NSMutableArray arrayWithCapacity:nameRows.count];
    for (NSDictionary *row in nameRows) {
        id v = row[@"name"];
        if ([v isKindOfClass:[NSData class]]) {
            v = [[NSString alloc] initWithData:v encoding:NSUTF8StringEncoding];
        }
        [nameTable addObject:(NSString *)v ?: @""];
    }
    const uint16_t *ids = (const uint16_t *)idsData.bytes;
    NSUInteger nIds = idsData.length / sizeof(uint16_t);
    NSMutableArray<NSString *> *chroms = [NSMutableArray arrayWithCapacity:nIds];
    for (NSUInteger i = 0; i < nIds; i++) {
        NSUInteger idx = ids[i];
        [chroms addObject:idx < nameTable.count ? nameTable[idx] : @""];
    }

    TTIOGenomicIndex *index = [[TTIOGenomicIndex alloc]
        initWithOffsets:offsets
                lengths:lengths
            chromosomes:chroms
              positions:positions
       mappingQualities:mapqs
                  flags:flags];

    // PO1 (P1.2): retain the per-read uint16 ids and invert the id→name
    // table (nameTable row k carries the name for id k, encounter-order
    // 1:1) into a name→uint16 map. -indicesForRegion: then resolves the
    // query chromosome once and scans uint16 comparisons. First
    // occurrence wins for any duplicate name (setObject:forKey: only on
    // absent key), matching Java/Python. Set ONLY on the disk-load path;
    // the in-memory designated initialiser leaves both ivars nil.
    NSMutableDictionary<NSString *, NSNumber *> *nameToId =
        [NSMutableDictionary dictionaryWithCapacity:nameTable.count];
    for (NSUInteger id_ = 0; id_ < nameTable.count; id_++) {
        NSString *name = nameTable[id_];
        if (nameToId[name] == nil) {
            nameToId[name] = @((uint16_t)id_);
        }
    }
    index->_chromosomeIdsData   = [idsData copy];
    index->_chromosomeNameToId  = [nameToId copy];

    return index;
}

static NSData *ttioPermute(NSData *d, NSUInteger elem, const uint32_t *p, NSUInteger n)
{
    NSMutableData *out = [NSMutableData dataWithLength:n * elem];
    const uint8_t *src = (const uint8_t *)d.bytes;
    uint8_t *dst = (uint8_t *)out.mutableBytes;
    NSUInteger have = d.length / elem;
    for (NSUInteger i = 0; i < n; i++) {
        if (p[i] < have) memcpy(dst + i * elem, src + (NSUInteger)p[i] * elem, elem);
    }
    return out;
}

- (instancetype)indexPermutedBy:(NSData *)storedOf
{
    NSUInteger n = storedOf.length / sizeof(uint32_t);
    const uint32_t *p = (const uint32_t *)storedOf.bytes;
    NSData *lengths = ttioPermute(_lengthsData, sizeof(uint32_t), p, n);
    NSMutableArray<NSString *> *chroms = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) [chroms addObject:_chromosomes[p[i]]];
    TTIOGenomicIndex *index = [[TTIOGenomicIndex alloc]
        initWithOffsets:TTIOOffsetsFromLengths(lengths)
                lengths:lengths
            chromosomes:chroms
              positions:ttioPermute(_positionsData, sizeof(int64_t), p, n)
       mappingQualities:ttioPermute(_mappingQualitiesData, sizeof(uint8_t), p, n)
                  flags:ttioPermute(_flagsData, sizeof(uint32_t), p, n)];
    if (_chromosomeIdsData != nil) {
        index->_chromosomeIdsData = ttioPermute(_chromosomeIdsData, sizeof(uint16_t), p, n);
        index->_chromosomeNameToId = _chromosomeNameToId;
    }
    return index;
}

@end
