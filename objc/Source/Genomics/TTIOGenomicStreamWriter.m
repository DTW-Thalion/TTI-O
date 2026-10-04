/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "Genomics/TTIOGenomicStreamWriter.h"
#import "Core/TTIOThreads.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOGenomicWriteContext.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Providers/TTIOCompoundField.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOHDF5Group.h"
#import "Dataset/TTIOCompoundIO.h"
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Codecs/TTIOQuality.h"
#import "Codecs/TTIOSeqGroup.h"
#import "Genomics/TTIOLazyReference.h"
#import <pthread.h>
#include <string.h>
#include <openssl/md5.h>  // reference-set digest for @reference_md5s

static NSString *const kLayout = @"blocks_v1";
/* blocks_v1 with the reads reordered by sequence and
 * genomic_index/input_index restoring input order (M103). A distinct
 * value so readers that predate it refuse the run rather than read it
 * in stored order (format-spec 10.12.7). */
static NSString *const kLayoutGrouped = @"blocks_v1_grouped";
static const NSUInteger kDefaultBlockReads = 1000000;
/* A block is the unit of encode concurrency and of residency, so a
 * big one costs both: the pool cannot start a block that is still
 * filling, and every block in flight is held whole. At 256 MiB a
 * 2x250 import of 10.6 M reads was 11 blocks, which left 7 workers
 * at 3.05x of one thread and 10.4 GiB resident. At 64 MiB the same
 * import is 43.0 s against 70.3 s and 4.50 GiB against 10.38 GiB.
 * The cost is context the codecs no longer share: 0.16% on 2x250,
 * 0.05% on 150 bp FASTQ, and 0.46% on HiFi, where reads are long
 * enough that a block holds few of them. */
static const unsigned long long kDefaultBlockBytes = 64ULL << 20;
static const NSUInteger kChannelChunk = 256 << 10;
static const NSUInteger kIndexChunkRows = 1024;
static const NSUInteger kIndexArrayChunk = 65536;

@implementation TTIOGenomicStreamWriterOptions

- (instancetype)init
{
    self = [super init];
    if (self) {
        _acquisitionMode = 0;
        _blockReads = kDefaultBlockReads;
        _blockBytes = kDefaultBlockBytes;
        _signalCodecOverrides = @{};
        _signalCompression = TTIOCompressionZlib;
        _provenanceRecords = @[];
    }
    return self;
}

+ (instancetype)defaultOptions
{
    return [[self alloc] init];
}

+ (instancetype)optionsFromRun:(TTIOWrittenGenomicRun *)run
{
    TTIOGenomicStreamWriterOptions *o = [[self alloc] init];
    o.acquisitionMode = run.acquisitionMode;
    o.referenceUri = run.referenceUri;
    o.platform = run.platform;
    o.sampleName = run.sampleName;
    o.readRole = run.readRole;
    o.refDiffSliceBytes = run.refDiffSliceBytes;
    o.referenceChromSeqs = run.referenceChromSeqs;
    o.embedReference = run.embedReference;
    o.optDisableQualitiesV5 = run.optDisableQualitiesV5;
    o.signalCodecOverrides = run.signalCodecOverrides ?: @{};
    o.signalCompression = run.signalCompression;
    o.optLegacyWholeChannel = run.optLegacyWholeChannel;
    o.groupReads = run.optGroupReads;
    o.provenanceRecords = run.provenanceRecords ?: @[];
    return o;
}

- (id)copyWithZone:(NSZone *)zone
{
    TTIOGenomicStreamWriterOptions *o = [[[self class] allocWithZone:zone] init];
    o.acquisitionMode = _acquisitionMode;
    o.referenceUri = _referenceUri;
    o.platform = _platform;
    o.sampleName = _sampleName;
    o.readRole = _readRole;
    o.refDiffSliceBytes = _refDiffSliceBytes;
    o.referenceChromSeqs = _referenceChromSeqs;
    o.embedReference = _embedReference;
    o.blockReads = _blockReads;
    o.blockBytes = _blockBytes;
    o.optDisableQualitiesV5 = _optDisableQualitiesV5;
    o.signalCodecOverrides = _signalCodecOverrides;
    o.signalCompression = _signalCompression;
    o.optLegacyWholeChannel = _optLegacyWholeChannel;
    o.groupReads = _groupReads;
    o.provenanceRecords = _provenanceRecords;
    o.threads = _threads;
    o.memoryBudgetBytes = _memoryBudgetBytes;
    return o;
}

@end

/** One block in flight: filled by the encode operation, written in
 *  sequence order by the caller's thread. */
@interface TTIOInFlightBlock : NSObject
@property (nonatomic, strong) TTIOWrittenGenomicRun *block;
@property (nonatomic) unsigned long long estimatedBytes;
@property (nonatomic, strong, nullable) TTIOBlockBlobs *blobs;
@property (nonatomic, strong, nullable) NSError *error;
@property (nonatomic) BOOL done;
@end

@implementation TTIOInFlightBlock
@end

@implementation TTIOGenomicStreamWriter {
    id<TTIOStorageGroup> _study;
    NSString *_name;
    TTIOGenomicStreamWriterOptions *_opt;
    NSMutableArray<TTIOWrittenGenomicRun *> *_pending;
    NSUInteger _pendingReads;
    unsigned long long _pendingBytes;
    NSString *_pendingChrom;
    NSMutableDictionary<NSString *, NSNumber *> *_chromMap;
    NSData *_referenceMD5;
    unsigned long long _readCount;
    unsigned long long _baseCount;
    NSUInteger _blockCount;
    /* Per-run sticky qualities strategy: block 0 auto-tunes, the winner
     * is pinned for the rest of the run. -1 = not yet pinned. */
    NSInteger _qualStrategyHint;
    BOOL _qualExhaustive;
    id<TTIOStorageGroup> _rg;
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *_channelDs;
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *_idxDs;
    id<TTIOStorageDataset> _indexDs;
    BOOL _withTags;     // the index carries the tags columns (M101)
    BOOL _embedded;
    BOOL _closed;
    NSMutableArray<TTIOWrittenGenomicRun *> *_legacyParts;
    /* M103 grouping: every batch until close, then the permutation
     * (uint32, row j = input index of the read stored at row j). */
    NSMutableArray<TTIOWrittenGenomicRun *> *_groupParts;
    NSData *_inputIndex;
    NSUInteger _threads;
    TTIOThreadPool *_pool;
    NSMutableArray<TTIOInFlightBlock *> *_inflight;
    NSCondition *_cond;
    unsigned long long _memoryBudgetBytes;
    unsigned long long _inflightBytes;
    unsigned long long _maxInFlightBytesObserved;
}

+ (NSString *)layout { return kLayout; }
+ (NSString *)groupedLayout { return kLayoutGrouped; }

+ (nullable NSString *)groupReadsRefusalForOptions:(TTIOGenomicStreamWriterOptions *)o
{
    if (!o.groupReads) return nil;
    if (o.optLegacyWholeChannel)
        return @"groupReads needs the blocks_v1 layout, not the legacy whole-channel one";
    if (o.referenceChromSeqs != nil)
        return @"groupReads applies to runs without a reference; aligned runs keep "
               @"coordinate order for REF_DIFF_V2";
    if (![TTIOSeqGroup nativeAvailable])
        return @"groupReads requires libttio_rans";
    return nil;
}
+ (NSUInteger)channelChunk { return kChannelChunk; }

static NSArray *gIndexFields = nil;
static NSArray *gIndexFieldsWithTags = nil;
static pthread_once_t gIndexFieldsOnce = PTHREAD_ONCE_INIT;

static void ttioBuildIndexFields(void)
{
    @autoreleasepool {
        NSMutableArray *f = [NSMutableArray array];
        [f addObject:[TTIOCompoundField fieldWithName:@"read_start" kind:TTIOCompoundFieldKindUInt64]];
        [f addObject:[TTIOCompoundField fieldWithName:@"n_reads" kind:TTIOCompoundFieldKindUInt32]];
        [f addObject:[TTIOCompoundField fieldWithName:@"base_start" kind:TTIOCompoundFieldKindUInt64]];
        [f addObject:[TTIOCompoundField fieldWithName:@"n_bases" kind:TTIOCompoundFieldKindUInt64]];
        for (NSString *ch in [TTIOGenomicBlocks requiredBlockChannels]) {
            [f addObject:[TTIOCompoundField fieldWithName:[ch stringByAppendingString:@"_off"]
                                                     kind:TTIOCompoundFieldKindUInt64]];
            [f addObject:[TTIOCompoundField fieldWithName:[ch stringByAppendingString:@"_len"]
                                                     kind:TTIOCompoundFieldKindUInt64]];
        }
        for (NSString *ch in [TTIOGenomicBlocks requiredBlockChannels]) {
            [f addObject:[TTIOCompoundField fieldWithName:[ch stringByAppendingString:@"_codec"]
                                                     kind:TTIOCompoundFieldKindUInt32]];
        }
        gIndexFields = [f copy];
        // M101: the tags triple trails the required columns.
        [f addObject:[TTIOCompoundField fieldWithName:@"tags_off" kind:TTIOCompoundFieldKindUInt64]];
        [f addObject:[TTIOCompoundField fieldWithName:@"tags_len" kind:TTIOCompoundFieldKindUInt64]];
        [f addObject:[TTIOCompoundField fieldWithName:@"tags_codec" kind:TTIOCompoundFieldKindUInt32]];
        gIndexFieldsWithTags = [f copy];
    }
}

+ (NSArray<TTIOCompoundField *> *)indexFields
{
    pthread_once(&gIndexFieldsOnce, ttioBuildIndexFields);
    return gIndexFields;
}

+ (NSArray<TTIOCompoundField *> *)indexFieldsWithTags
{
    pthread_once(&gIndexFieldsOnce, ttioBuildIndexFields);
    return gIndexFieldsWithTags;
}

+ (NSArray<TTIOCompoundField *> *)indexFieldsForRows:(NSArray<NSDictionary *> *)rows
{
    BOOL has = rows.count > 0 && rows[0][@"tags_off"] != nil;
    return has ? [self indexFieldsWithTags] : [self indexFields];
}

- (instancetype)initWithStudyGroup:(id<TTIOStorageGroup>)study
                           runName:(NSString *)runName
                           options:(TTIOGenomicStreamWriterOptions *)options
{
    self = [super init];
    if (self) {
        _study = study;
        _name = [runName copy];
        _opt = options ? [options copy] : [TTIOGenomicStreamWriterOptions defaultOptions];
        NSNumber *qov = _opt.signalCodecOverrides[@"qualities"];
        if (qov != nil
            && qov.integerValue == TTIOCompressionQualityBinned
            && !TTIOQualityBinnedAllowedForPlatform(_opt.platform)) {
            [NSException raise:NSInvalidArgumentException
                        format:@"signalCodecOverrides['qualities']: "
                               @"QUALITY_BINNED applies the fixed "
                               @"Illumina-8 bin table, which does not fit "
                               @"the quality distribution of platform "
                               @"'%@' (M97).", _opt.platform];
        }
        NSString *refusal = [[self class] groupReadsRefusalForOptions:_opt];
        if (refusal != nil) {
            [NSException raise:NSInvalidArgumentException format:@"%@", refusal];
        }
        if (_opt.blockReads < 1) _opt.blockReads = 1;
        if (_opt.blockBytes < 1) _opt.blockBytes = 1;
        _pending = [NSMutableArray array];
        _chromMap = [NSMutableDictionary dictionary];
        _channelDs = [NSMutableDictionary dictionary];
        _idxDs = [NSMutableDictionary dictionary];
        _legacyParts = [NSMutableArray array];
        _groupParts = [NSMutableArray array];
        _threads = [TTIOThreads resolve:_opt.threads ? @(_opt.threads) : nil];
        _pool = [TTIOThreadPool poolWithThreads:_opt.optLegacyWholeChannel ? 1 : _threads];
        _inflight = [NSMutableArray array];
        _cond = [NSCondition new];
        _memoryBudgetBytes = [TTIOThreads resolveMemoryBudget:
                                  _opt.memoryBudgetBytes ? @(_opt.memoryBudgetBytes) : nil
                                                      threads:_threads
                                                   blockBytes:_opt.blockBytes];
        _qualStrategyHint = -1;
        const char *ex = getenv("TTIO_M94Z_EXHAUSTIVE");
        _qualExhaustive = (ex != NULL && strcmp(ex, "1") == 0);
        /* A hint pinned before block 0 skips the tune and holds for the
         * run. V6 is reachable no other way, so a benchmark that wants
         * it has to say so; see TtioGenomicWriteBench. */
        const char *hint = getenv("TTIO_M94Z_HINT");
        if (!_qualExhaustive && hint && hint[0]) {
            long v = strtol(hint, NULL, 10);
            if (v > 0) _qualStrategyHint = (NSInteger)v;
        }
    }
    return self;
}

- (NSInteger)qualStrategyHint { return _qualStrategyHint; }

- (NSUInteger)threads { return _threads; }
- (unsigned long long)memoryBudgetBytes { return _memoryBudgetBytes; }
- (unsigned long long)maxInFlightBytesObserved { return _maxInFlightBytesObserved; }

/** The bytes a pending block is expected to hold while encoding: its
 *  raw channel buffers times two for codec workspace. The marginal
 *  resident cost measured on 64 MiB HiFi blocks is near raw * 1.5;
 *  the old raw * 4 halved the encode concurrency the budget admits. */
static unsigned long long ppEstimateBlockBytes(TTIOWrittenGenomicRun *b)
{
    unsigned long long raw = (unsigned long long)b.sequencesData.length
        + (unsigned long long)b.qualitiesData.length
        + (unsigned long long)b.offsetsData.length * 3ull;
    return raw * 2ull;
}

+ (void)registerBlockChromosomes:(TTIOWrittenGenomicRun *)block
                         intoMap:(NSMutableDictionary<NSString *, NSNumber *> *)map
{
    for (NSString *n in block.chromosomes) {
        if (map[n] == nil) map[n] = @(map.count);
    }
    for (NSString *n in block.mateChromosomes) {
        if (n.length > 0 && ![n isEqualToString:@"*"] && ![n isEqualToString:@"="] && map[n] == nil) {
            map[n] = @(map.count);
        }
    }
}

- (unsigned long long)readCount { return _readCount; }
- (NSUInteger)blockCount { return _blockCount; }
- (TTIOGenomicStreamWriterOptions *)options { return [_opt copy]; }

// ── run construction helpers ─────────────────────────────────────

- (TTIOWrittenGenomicRun *)_runWithMeta:(TTIOWrittenGenomicRun *)run
{
    TTIOWrittenGenomicRun *r = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:_opt.acquisitionMode
                   referenceUri:_opt.referenceUri
                       platform:_opt.platform
                     sampleName:_opt.sampleName
                      positions:run.positionsData
               mappingQualities:run.mappingQualitiesData
                          flags:run.flagsData
                      sequences:run.sequencesData
                      qualities:run.qualitiesData
                        offsets:run.offsetsData
                        lengths:run.lengthsData
                         cigars:run.cigars
                      readNames:run.readNames
                mateChromosomes:run.mateChromosomes
                  matePositions:run.matePositionsData
                templateLengths:run.templateLengthsData
                    chromosomes:run.chromosomes
              signalCompression:_opt.signalCompression
           signalCodecOverrides:_opt.signalCodecOverrides];
    r.optDisableQualitiesV5 = _opt.optDisableQualitiesV5;
    r.embedReference = _opt.embedReference;
    r.referenceChromSeqs = _opt.referenceChromSeqs;
    r.provenanceRecords = @[];
    r.optLegacyWholeChannel = NO;
    r.readRole = _opt.readRole;
    r.refDiffSliceBytes = _opt.refDiffSliceBytes;
    r.tags = run.tags;
    return r;
}

- (TTIOWrittenGenomicRun *)_singleReadRun:(TTIOAlignedRead *)r
{
    NSData *seq = r.sequence ? [r.sequence dataUsingEncoding:NSASCIIStringEncoding] : nil;
    if (!seq) seq = [NSData data];
    NSData *qual = r.qualities ?: [NSData data];
    NSString *mate = (r.mateChromosome == nil || r.mateChromosome.length == 0) ? @"*" : r.mateChromosome;
    int64_t pos = r.position;
    uint8_t mapq = r.mappingQuality;
    uint32_t flags = r.flags;
    uint64_t off = 0;
    uint32_t len = (uint32_t)seq.length;
    int64_t mpos = r.matePosition;
    int32_t tlen = r.templateLength;
    TTIOWrittenGenomicRun *run = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:_opt.acquisitionMode
                   referenceUri:_opt.referenceUri
                       platform:_opt.platform
                     sampleName:_opt.sampleName
                      positions:[NSData dataWithBytes:&pos length:sizeof(pos)]
               mappingQualities:[NSData dataWithBytes:&mapq length:1]
                          flags:[NSData dataWithBytes:&flags length:sizeof(flags)]
                      sequences:seq
                      qualities:qual
                        offsets:[NSData dataWithBytes:&off length:sizeof(off)]
                        lengths:[NSData dataWithBytes:&len length:sizeof(len)]
                         cigars:@[r.cigar ?: @"*"]
                      readNames:@[r.readName ?: @"*"]
                mateChromosomes:@[mate]
                  matePositions:[NSData dataWithBytes:&mpos length:sizeof(mpos)]
                templateLengths:[NSData dataWithBytes:&tlen length:sizeof(tlen)]
                    chromosomes:@[r.chromosome ?: @"*"]
              signalCompression:_opt.signalCompression
           signalCodecOverrides:_opt.signalCodecOverrides];
    run.optDisableQualitiesV5 = _opt.optDisableQualitiesV5;
    run.embedReference = _opt.embedReference;
    run.referenceChromSeqs = _opt.referenceChromSeqs;
    run.readRole = _opt.readRole;
    run.refDiffSliceBytes = _opt.refDiffSliceBytes;
    run.tags = @[r.tags ?: @""];
    return run;
}

// ── appending ────────────────────────────────────────────────────

- (BOOL)appendRead:(TTIOAlignedRead *)read error:(NSError **)error
{
    return [self appendBatch:[self _singleReadRun:read] error:error];
}

- (BOOL)appendBatch:(TTIOWrittenGenomicRun *)batch error:(NSError **)error
{
    if (_closed) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument, @"genomic stream writer is closed");
        return NO;
    }
    return [self _append:batch error:error];
}

- (BOOL)_append:(TTIOWrittenGenomicRun *)batch error:(NSError **)error
{
    NSUInteger n = batch.readCount;
    if (n == 0) return YES;
    if (_opt.optLegacyWholeChannel) {
        [_legacyParts addObject:batch];
        return YES;
    }
    if (_opt.groupReads && _inputIndex == nil) {
        [_groupParts addObject:batch];
        return YES;
    }
    NSArray<NSString *> *chroms = batch.chromosomes;
    const uint32_t *lens = (const uint32_t *)batch.lengthsData.bytes;
    NSUInteger start = 0;
    while (start < n) {
        NSString *chrom = chroms[start];
        NSUInteger segEnd = start + 1;
        while (segEnd < n && [chroms[segEnd] isEqualToString:chrom]) segEnd++;
        if (_pending.count > 0 && ![chrom isEqualToString:_pendingChrom]) {
            if (![self flush:error]) return NO;
        }
        _pendingChrom = chrom;
        while (start < segEnd) {
            NSUInteger roomReads = _opt.blockReads > _pendingReads ? _opt.blockReads - _pendingReads : 0;
            unsigned long long roomBytes = _opt.blockBytes > _pendingBytes ? _opt.blockBytes - _pendingBytes : 0;
            NSUInteger stop = MIN(segEnd, start + MAX(roomReads, (NSUInteger)1));
            unsigned long long cum = 0;
            NSUInteger fit = 0;
            for (NSUInteger i = start; i < stop; i++) {
                cum += lens[i];
                if (cum <= roomBytes) fit++; else break;
            }
            if (fit < stop - start) stop = start + MAX(fit, (NSUInteger)1);
            TTIOWrittenGenomicRun *part = (start == 0 && stop == n)
                ? batch : [TTIOGenomicBlocks sliceRun:batch from:start to:stop];
            [_pending addObject:part];
            _pendingReads += stop - start;
            for (NSUInteger i = start; i < stop; i++) _pendingBytes += lens[i];
            if (_pendingReads >= _opt.blockReads || _pendingBytes >= _opt.blockBytes) {
                if (![self _cutBlock:error]) return NO;
            }
            start = stop;
        }
    }
    return YES;
}

- (BOOL)flush:(NSError **)error
{
    /* The public barrier: cut the pending block and write everything in
     * flight, so an unclosed file's flushed reads are on storage. */
    if (![self _cutBlock:error]) return NO;
    return [self _drainUntil:0 error:error];
}

- (BOOL)_cutBlock:(NSError **)error
{
    if (_opt.optLegacyWholeChannel || _pending.count == 0) return YES;
    TTIOWrittenGenomicRun *block = [TTIOGenomicBlocks concatRuns:_pending];
    [_pending removeAllObjects];
    _pendingReads = 0;
    _pendingBytes = 0;
    block = [self _runWithMeta:block];
    if (_referenceMD5 == nil && _opt.referenceChromSeqs != nil) {
        _referenceMD5 = [TTIOSpectralDataset referenceMD5ForRun:block];
    }
    if (!_embedded && _opt.embedReference) {
        if (![TTIOSpectralDataset embedReferencesForRuns:@[block] inStudy:_study error:error]) return NO;
        _embedded = YES;
    }
    [TTIOSpectralDataset validateGenomicCodecOverridesForRun:block];
    [[self class] registerBlockChromosomes:block intoMap:_chromMap];
    if (!_qualExhaustive && _qualStrategyHint == -1
        && (_blockCount > 0 || _inflight.count > 0)) {
        /* Block-0 gate: the pin comes from the earliest written M94Z
         * block; hold later submissions until it is known so the
         * choice is timing-independent. */
        if (![self _drainUntil:0 error:error]) return NO;
    }
    TTIOGenomicWriteContext *ctx =
        [TTIOGenomicWriteContext contextWithChromNameToId:_chromMap referenceMD5:_referenceMD5];
    ctx.qualStrategyHint = _qualStrategyHint;
    if (_pool.queue == nil) {
        TTIOBlockBlobs *blobs = [TTIOGenomicBlocks encodeBlock:block context:ctx error:error];
        if (!blobs) return NO;
        return [self _writeEncoded:block blobs:blobs error:error];
    }
    if (![self _drainUntil:_threads error:error]) return NO;
    if (![self _drainToBudget:error]) return NO;
    /* The worker reads the map while later flushes mutate it: give each
     * block a snapshot (registration above fixed every id it needs). */
    TTIOGenomicWriteContext *bctx =
        [TTIOGenomicWriteContext contextWithChromNameToId:[_chromMap mutableCopy]
                                             referenceMD5:_referenceMD5];
    bctx.qualStrategyHint = _qualStrategyHint;
    TTIOInFlightBlock *f = [TTIOInFlightBlock new];
    f.block = block;
    f.estimatedBytes = ppEstimateBlockBytes(block);
    [_inflight addObject:f];
    _inflightBytes += f.estimatedBytes;
    if (_inflightBytes > _maxInFlightBytesObserved) _maxInFlightBytesObserved = _inflightBytes;
    /* Segment threads follow the blocks actually in flight, not the
     * pool's size. A run with fewer blocks than workers never fills the
     * pool, and sizing the segments from the worker count leaves the
     * machine idle in exactly that case: 3 blocks against 30 workers
     * asked for 2 segment threads each and used 6 cores of 32. The
     * product stays near the core count as the count settles, and the
     * pool restores the previous value when it closes. */
    [TTIOThreads applyV6SegmentThreadsForBlocksInFlight:_inflight.count];
    NSCondition *cond = _cond;
    [_pool.queue addOperationWithBlock:^{
        NSError *e = nil;
        TTIOBlockBlobs *b = nil;
        @try {
            b = [TTIOGenomicBlocks encodeBlock:block context:bctx error:&e];
        } @catch (NSException *ex) {
            /* NSOperationQueue swallows exceptions; surface them as the
             * block's error so the drain completes (and the serial path's
             * raise becomes the caller's NSError). */
            b = nil;
            e = TTIOMakeError(TTIOErrorDatasetWrite, @"%@: %@", ex.name, ex.reason);
        }
        [cond lock];
        f.blobs = b;
        f.error = e;
        f.done = YES;
        [cond broadcast];
        [cond unlock];
    }];
    return YES;
}

/** Write completed blocks in sequence order; wait on the oldest until at
 *  most blockUntil remain in flight. */
- (BOOL)_drainUntil:(NSUInteger)blockUntil error:(NSError **)error
{
    while (_inflight.count > 0) {
        TTIOInFlightBlock *f = _inflight.firstObject;
        [_cond lock];
        if (_inflight.count <= blockUntil && !f.done) {
            [_cond unlock];
            break;
        }
        while (!f.done) [_cond wait];
        [_cond unlock];
        [_inflight removeObjectAtIndex:0];
        _inflightBytes -= f.estimatedBytes <= _inflightBytes ? f.estimatedBytes : _inflightBytes;
        if (!f.blobs) {
            if (error) *error = f.error;
            return NO;
        }
        if (![self _writeEncoded:f.block blobs:f.blobs error:error]) return NO;
    }
    return YES;
}

/** Drain until the in-flight estimate fits the writer's half of the
 *  byte budget (the count window of _drainUntil: stays the upper
 *  bound; this adds the byte bound the count cannot give). */
- (BOOL)_drainToBudget:(NSError **)error
{
    unsigned long long half = _memoryBudgetBytes / 2ull;
    while (_inflight.count > 0 && _inflightBytes > half) {
        if (![self _drainUntil:(_inflight.count - 1) error:error]) return NO;
    }
    return YES;
}

- (BOOL)_writeEncoded:(TTIOWrittenGenomicRun *)block blobs:(TTIOBlockBlobs *)blobs error:(NSError **)error
{
    BOOL blockHasTags = blobs.blobs[@"tags"].length > 0;
    // The first block decides the index schema; tags arriving later
    // upgrade it once (M101).
    if (_rg == nil) _withTags = blockHasTags;
    if (![self _ensureLayout:error]) return NO;
    if (blockHasTags && !_withTags && ![self _addTagsColumns:error]) return NO;

    NSArray<TTIOCompoundField *> *fields = _withTags
        ? [[self class] indexFieldsWithTags] : [[self class] indexFields];
    NSMutableDictionary *row = [NSMutableDictionary dictionaryWithCapacity:fields.count];
    row[@"read_start"] = @(_readCount);
    row[@"n_reads"] = @((uint32_t)blobs.nReads);
    row[@"base_start"] = @(_baseCount);
    row[@"n_bases"] = @(blobs.nBases);
    for (NSString *ch in [TTIOGenomicBlocks blockChannels]) {
        if ([ch isEqualToString:@"tags"] && !_withTags) continue;
        NSData *data = blobs.blobs[ch] ?: [NSData data];
        NSNumber *codec = blobs.codecs[ch] ?: @0;
        row[[ch stringByAppendingString:@"_codec"]] = @((uint32_t)[codec unsignedIntegerValue]);
        id<TTIOStorageDataset> ds = _channelDs[ch];
        if (ds == nil) {
            if (data.length > 0) {
                ds = [self _createChannel:ch blobs:blobs error:error];
                if (!ds) return NO;
            } else {
                row[[ch stringByAppendingString:@"_off"]] = @0ULL;
                row[[ch stringByAppendingString:@"_len"]] = @0ULL;
                continue;
            }
        }
        row[[ch stringByAppendingString:@"_off"]] = @((unsigned long long)[ds length]);
        row[[ch stringByAppendingString:@"_len"]] = @((unsigned long long)data.length);
        if (data.length > 0 && ![ds appendData:data error:error]) return NO;
    }
    if (![_indexDs appendData:@[row] error:error]) return NO;
    if (![self _appendIndexArrays:block error:error]) return NO;
    _readCount += blobs.nReads;
    _baseCount += blobs.nBases;
    _blockCount++;
    if (_qualStrategyHint == -1 && !_qualExhaustive
        && [blobs.codecs[@"qualities"] unsignedIntegerValue]
             == TTIOCompressionFqzcompNx16Z) {
        NSInteger strat = [TTIOFqzcompNx16Z
            strategyOfEncodedStream:blobs.blobs[@"qualities"]];
        if (strat > 0) {
            _qualStrategyHint = strat == 4 ? TTIOM94ZHintV4Auto : strat;
        }
    }
    if (![_rg setAttributeValue:@((int64_t)_readCount) forName:@"read_count" error:error]) return NO;
    if (![_rg setAttributeValue:@((int64_t)_baseCount) forName:@"base_count" error:error]) return NO;
    return YES;
}

- (BOOL)close:(NSError **)error
{
    if (_closed) return YES;
    _closed = YES;
    if (_opt.optLegacyWholeChannel) {
        if (_legacyParts.count > 0) {
            TTIOWrittenGenomicRun *whole =
                [[[self _runWithMeta:[TTIOGenomicBlocks concatRuns:_legacyParts]]
                    copyWithProvenance:_opt.provenanceRecords]
                    copyWithOptLegacyWholeChannel:YES];
            if (_opt.embedReference) {
                if (![TTIOSpectralDataset embedReferencesForRuns:@[whole] inStudy:_study error:error]) return NO;
            }
            id<TTIOStorageGroup> g = [self _runsGroup:error];
            if (!g) return NO;
            if (![TTIOSpectralDataset writeGenomicRunStorage:whole toGroup:g name:_name
                                                     context:[TTIOGenomicWriteContext none]
                                                       error:error]) return NO;
            _readCount = whole.readCount;
        }
        [_legacyParts removeAllObjects];
        [_pool close];
        return YES;
    }
    if (_opt.groupReads) {
        /* M103: reorder the whole run by sequence, then block it as
         * usual; each written block appends its slice of the order. */
        NSArray *parts = [_groupParts copy];
        [_groupParts removeAllObjects];
        TTIOWrittenGenomicRun *whole = parts.count ? [TTIOGenomicBlocks concatRuns:parts] : nil;
        NSData *order = whole ? [TTIOGenomicBlocks groupOrderOfRun:whole error:error] : [NSData data];
        if (order == nil) { [_pool close]; return NO; }
        _inputIndex = order;
        if (whole != nil
            && ![self _append:[TTIOGenomicBlocks takeRun:whole order:order] error:error]) {
            [_pool close];
            return NO;
        }
    }
    if (![self _cutBlock:error]) { [_pool close]; return NO; }
    if (![self _drainUntil:0 error:error]) { [_pool close]; return NO; }
    [_pool close];
    if (_rg == nil && ![self _ensureLayout:error]) return NO;
    if (![self _writeCloseTables:error]) return NO;
    if (_opt.provenanceRecords.count > 0) {
        id<TTIOStorageGroup> prov = [_rg createGroupNamed:@"provenance" error:error];
        if (!prov) return NO;
        if ([prov respondsToSelector:@selector(unwrap)]) {
            TTIOHDF5Group *h5 = [(id)prov performSelector:@selector(unwrap)];
            if (h5 && ![TTIOCompoundIO writeProvenance:_opt.provenanceRecords intoGroup:h5
                                          datasetNamed:@"steps" error:error]) return NO;
        }
        NSMutableArray *plists = [NSMutableArray arrayWithCapacity:_opt.provenanceRecords.count];
        for (TTIOProvenanceRecord *r in _opt.provenanceRecords) [plists addObject:[r asPlist]];
        NSError *jErr = nil;
        NSData *json = [NSJSONSerialization dataWithJSONObject:plists options:0 error:&jErr];
        if (!json) { if (error) *error = jErr; return NO; }
        NSString *jstr = [[NSString alloc] initWithData:json encoding:NSUTF8StringEncoding];
        if (![_rg setAttributeValue:jstr forName:@"provenance_json" error:error]) return NO;
    }
    return YES;
}

// ── layout ───────────────────────────────────────────────────────

- (id<TTIOStorageGroup>)_runsGroup:(NSError **)error
{
    id<TTIOStorageGroup> g;
    if ([_study hasChildNamed:@"genomic_runs"]) {
        g = [_study openGroupNamed:@"genomic_runs" error:error];
        if (!g) return nil;
    } else {
        g = [_study createGroupNamed:@"genomic_runs" error:error];
        if (!g) return nil;
        if (![g setAttributeValue:@"" forName:@"_run_names" error:error]) return nil;
    }
    id namesAttr = [g hasAttributeNamed:@"_run_names"] ? [g attributeValueForName:@"_run_names" error:NULL] : nil;
    NSString *names = namesAttr ? [namesAttr description] : @"";
    NSMutableArray *list = [NSMutableArray array];
    for (NSString *s in [names componentsSeparatedByString:@","]) if (s.length) [list addObject:s];
    if (![list containsObject:_name]) {
        [list addObject:_name];
        if (![g setAttributeValue:[list componentsJoinedByString:@","] forName:@"_run_names" error:error]) return nil;
    }
    return g;
}

- (BOOL)_ensureLayout:(NSError **)error
{
    if (_rg != nil) return YES;
    id<TTIOStorageGroup> g = [self _runsGroup:error];
    if (!g) return NO;
    if ([g hasChildNamed:_name]) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"genomic run '%@' already exists", _name);
        return NO;
    }
    id<TTIOStorageGroup> run = [g createGroupNamed:_name error:error];
    if (!run) return NO;
    if (![run setAttributeValue:@((int64_t)_opt.acquisitionMode) forName:@"acquisition_mode" error:error]) return NO;
    if (![run setAttributeValue:@"genomic_sequencing" forName:@"modality" error:error]) return NO;
    if (![run setAttributeValue:@((int64_t)5) forName:@"spectrum_class" error:error]) return NO;
    if (![run setAttributeValue:_opt.referenceUri ?: @"" forName:@"reference_uri" error:error]) return NO;
    if (![run setAttributeValue:_opt.platform ?: @"" forName:@"platform" error:error]) return NO;
    if (_opt.readRole.length > 0
        && ![run setAttributeValue:_opt.readRole forName:@"read_role" error:error]) return NO;
    if (![run setAttributeValue:_opt.sampleName ?: @"" forName:@"sample_name" error:error]) return NO;
    if (![run setAttributeValue:@((int64_t)0) forName:@"read_count" error:error]) return NO;
    if (![run setAttributeValue:@((int64_t)0) forName:@"base_count" error:error]) return NO;
    if (![run setAttributeValue:(_opt.groupReads ? kLayoutGrouped : kLayout)
                        forName:@"layout" error:error]) return NO;
    NSString *policy = [NSString stringWithFormat:@"reads=%lu,bytes=%llu",
                        (unsigned long)_opt.blockReads, _opt.blockBytes];
    if (![run setAttributeValue:policy forName:@"block_policy" error:error]) return NO;
    /* Non-default writer policy shapes the coded blobs; a per-AU
     * restore re-encodes with it, so it is persisted on the run
     * (format-spec 9.1.1). */
    if (_opt.refDiffSliceBytes != 0
        && ![run setAttributeValue:@((int64_t)_opt.refDiffSliceBytes)
                           forName:@"ref_diff_slice_bytes" error:error]) return NO;
    if (_opt.optDisableQualitiesV5
        && ![run setAttributeValue:@((int64_t)1)
                           forName:@"opt_disable_qualities_v5" error:error]) return NO;
    if (_opt.referenceChromSeqs != nil) {
        /* The chromosome set and the reference-set digest (the md5
         * the REF_DIFF blob headers carry and the resolver verifies),
         * so a restore can rebuild the reference from
         * /study/references or REF_PATH (format-spec 9.1.1). */
        NSData *md5 = _referenceMD5;
        if (md5 == nil) {
            if ([_opt.referenceChromSeqs
                    isKindOfClass:[TTIOLazyReference class]]) {
                md5 = [(TTIOLazyReference *)_opt.referenceChromSeqs setMD5];
            } else {
                NSArray *sortedNames = [[_opt.referenceChromSeqs allKeys]
                    sortedArrayUsingSelector:@selector(compare:)];
                uint8_t digest[16];
                MD5_CTX c; MD5_Init(&c);
                for (NSString *nm in sortedNames) {
                    NSData *sq = _opt.referenceChromSeqs[nm];
                    MD5_Update(&c, sq.bytes, sq.length);
                }
                MD5_Final(digest, &c);
                md5 = [NSData dataWithBytes:digest length:16];
            }
            _referenceMD5 = md5;
        }
        NSMutableString *hex = [NSMutableString stringWithCapacity:32];
        const uint8_t *mb = md5.bytes;
        for (NSUInteger i = 0; i < md5.length; i++) {
            [hex appendFormat:@"%02x", mb[i]];
        }
        NSArray *names = [[_opt.referenceChromSeqs allKeys]
            sortedArrayUsingSelector:@selector(compare:)];
        NSMutableString *jsonStr = [NSMutableString stringWithString:@"{"];
        for (NSUInteger i = 0; i < names.count; i++) {
            if (i) [jsonStr appendString:@","];
            [jsonStr appendFormat:@"\"%@\":\"%@\"", names[i], hex];
        }
        [jsonStr appendString:@"}"];
        if (![run setAttributeValue:jsonStr
                            forName:@"reference_md5s" error:error]) return NO;
    }
    id<TTIOStorageGroup> blocks = [run createGroupNamed:@"blocks" error:error];
    if (!blocks) return NO;
    _indexDs = [blocks createCompoundDatasetNamed:@"index"
                                           fields:(_withTags ? [[self class] indexFieldsWithTags]
                                                             : [[self class] indexFields])
                                            count:0 extendable:YES chunkRows:kIndexChunkRows error:error];
    if (!_indexDs) return NO;
    id<TTIOStorageGroup> idx = [run createGroupNamed:@"genomic_index" error:error];
    if (!idx) return NO;
    NSArray *arrays = @[
        @[@"lengths", @(TTIOPrecisionUInt32)],
        @[@"positions", @(TTIOPrecisionInt64)],
        @[@"mapping_qualities", @(TTIOPrecisionUInt8)],
        @[@"flags", @(TTIOPrecisionUInt32)],
        @[@"chromosome_ids", @(TTIOPrecisionUInt16)],
    ];
    if (_opt.groupReads) {
        // M103: row j holds the input index of the read stored at row j.
        arrays = [arrays arrayByAddingObject:@[@"input_index", @(TTIOPrecisionUInt32)]];
    }
    for (NSArray *a in arrays) {
        id<TTIOStorageDataset> ds = [idx createDatasetNamed:a[0]
                                                  precision:(TTIOPrecision)[a[1] integerValue]
                                                     length:0
                                                  chunkSize:kIndexArrayChunk
                                                compression:TTIOCompressionZlib
                                           compressionLevel:6
                                                 extendable:YES
                                                      error:error];
        if (!ds) return NO;
        _idxDs[a[0]] = ds;
    }
    if (![run createGroupNamed:@"signal_channels" error:error]) return NO;
    _rg = run;
    return YES;
}

/** The first block with tags arrives after blocks without: give the
 *  index the tags columns, the earlier blocks an empty range. */
- (BOOL)_addTagsColumns:(NSError **)error
{
    id<TTIOStorageGroup> blocks = [_rg openGroupNamed:@"blocks" error:error];
    if (!blocks) return NO;
    NSMutableArray *rows = [NSMutableArray array];
    if (_blockCount > 0) {
        NSArray *old = [_indexDs readRows:error];
        if (!old) return NO;
        for (NSDictionary *r in old) {
            NSMutableDictionary *m = [r mutableCopy];
            m[@"tags_off"] = @0ULL;
            m[@"tags_len"] = @0ULL;
            m[@"tags_codec"] = @((uint32_t)0);
            [rows addObject:m];
        }
    }
    _indexDs = nil;
    if (![blocks deleteChildNamed:@"index" error:error]) return NO;
    _indexDs = [blocks createCompoundDatasetNamed:@"index" fields:[[self class] indexFieldsWithTags]
                                            count:0 extendable:YES chunkRows:kIndexChunkRows error:error];
    if (!_indexDs) return NO;
    if (rows.count > 0 && ![_indexDs appendData:rows error:error]) return NO;
    _withTags = YES;
    return YES;
}

- (id<TTIOStorageDataset>)_createChannel:(NSString *)ch
                                   blobs:(TTIOBlockBlobs *)blobs
                                   error:(NSError **)error
{
    id<TTIOStorageGroup> sc = [_rg openGroupNamed:@"signal_channels" error:error];
    if (!sc) return nil;
    id<TTIOStorageGroup> parent;
    NSString *dsName;
    if ([ch isEqualToString:@"sequences"]) {
        parent = [sc createGroupNamed:@"sequences" error:error];
        dsName = @"data";
    } else if ([ch isEqualToString:@"mate_info"]) {
        parent = [sc createGroupNamed:@"mate_info" error:error];
        dsName = @"inline_v2";
    } else {
        parent = sc;
        dsName = ch;
    }
    if (!parent) return nil;
    NSUInteger codec = [blobs.codecs[ch] unsignedIntegerValue];
    id<TTIOStorageDataset> ds = [parent createDatasetNamed:dsName
                                                 precision:TTIOPrecisionUInt8
                                                    length:0
                                                 chunkSize:kChannelChunk
                                               compression:codec == 0 ? TTIOCompressionZlib : TTIOCompressionNone
                                          compressionLevel:6
                                                extendable:YES
                                                     error:error];
    if (!ds) return nil;
    if (![ds setAttributeValue:@((int64_t)codec) forName:@"compression" error:error]) return nil;
    NSDictionary *extra = blobs.extraAttrs[ch] ?: @{};
    for (NSString *k in [[extra allKeys] sortedArrayUsingSelector:@selector(compare:)]) {
        if (![ds setAttributeValue:extra[k] forName:k error:error]) return nil;
    }
    _channelDs[ch] = ds;
    return ds;
}

- (BOOL)_appendIndexArrays:(TTIOWrittenGenomicRun *)block error:(NSError **)error
{
    NSUInteger n = block.readCount;
    NSMutableData *ids = [NSMutableData dataWithLength:n * sizeof(uint16_t)];
    uint16_t *idv = (uint16_t *)ids.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) {
        NSNumber *id_ = _chromMap[block.chromosomes[i]];
        if (id_ == nil) {
            if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
                @"chromosome '%@' missing from the shared id map", block.chromosomes[i]);
            return NO;
        }
        idv[i] = (uint16_t)[id_ unsignedIntegerValue];
    }
    if (![_idxDs[@"lengths"] appendData:block.lengthsData error:error]) return NO;
    if (![_idxDs[@"positions"] appendData:block.positionsData error:error]) return NO;
    if (![_idxDs[@"mapping_qualities"] appendData:block.mappingQualitiesData error:error]) return NO;
    if (![_idxDs[@"flags"] appendData:block.flagsData error:error]) return NO;
    if (![_idxDs[@"chromosome_ids"] appendData:ids error:error]) return NO;
    if (_opt.groupReads) {
        /* Blocks are written in order, so this block's rows start at the
         * reads written so far. */
        NSUInteger r0 = (NSUInteger)_readCount;
        if ((r0 + n) * sizeof(uint32_t) > _inputIndex.length) {
            if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
                @"grouped block past the end of input_index");
            return NO;
        }
        NSData *slice = [_inputIndex subdataWithRange:NSMakeRange(r0 * sizeof(uint32_t),
                                                                 n * sizeof(uint32_t))];
        if (![_idxDs[@"input_index"] appendData:slice error:error]) return NO;
    }
    return YES;
}

- (BOOL)_writeCloseTables:(NSError **)error
{
    NSArray *nameFields = @[[TTIOCompoundField fieldWithName:@"name" kind:TTIOCompoundFieldKindVLString]];
    NSMutableArray *rows = [NSMutableArray array];
    for (NSString *n in [TTIOGenomicIndex namesInIdOrder:_chromMap]) [rows addObject:@{@"name": n}];
    id<TTIOStorageGroup> idx = [_rg openGroupNamed:@"genomic_index" error:error];
    if (!idx) return NO;
    id<TTIOStorageDataset> ds = [idx createCompoundDatasetNamed:@"chromosome_names"
                                                          fields:nameFields count:rows.count error:error];
    if (!ds || ![ds writeAll:rows error:error]) return NO;
    id<TTIOStorageGroup> sc = [_rg openGroupNamed:@"signal_channels" error:error];
    if (!sc) return NO;
    id<TTIOStorageGroup> mate = [sc hasChildNamed:@"mate_info"]
        ? [sc openGroupNamed:@"mate_info" error:error]
        : [sc createGroupNamed:@"mate_info" error:error];
    if (!mate) return NO;
    if (![mate hasChildNamed:@"chrom_names"]) {
        id<TTIOStorageDataset> cn = [mate createCompoundDatasetNamed:@"chrom_names"
                                                              fields:nameFields count:rows.count error:error];
        if (!cn || ![cn writeAll:rows error:error]) return NO;
    }
    return YES;
}

@end
