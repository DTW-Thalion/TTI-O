/*
 * TTIOGenomicRun.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOGenomicRun
 * Inherits From: NSObject
 * Conforms To:   TTIOIndexable, TTIORun
 * Declared In:   Genomics/TTIOGenomicRun.h
 *
 * Lazy view over /study/genomic_runs/<name>/. The genomic_index/
 * subgroup is loaded eagerly; per-read materialisation of
 * TTIOAlignedRead happens on demand. Decoded read-name and CIGAR
 * lists are cached for the lifetime of the run instance.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOGenomicRun.h"
#import "Core/TTIOThreads.h"
#import "TTIOAlignedRead.h"
#import "TTIOGenomicIndex.h"
#import "TTIOBlockTable.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "TTIOBlockView.h"
#import "TTIOGenomicBlocks.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOCompoundField.h"
#import "Dataset/TTIOCompoundIO.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "HDF5/TTIOHDF5Group.h"
#import "HDF5/TTIOHDF5Dataset.h"
#import "Codecs/TTIORans.h"
#import "Codecs/TTIOBasePack.h"
#import "Codecs/TTIOQuality.h"   // M86 Phase D
// TTIONameTokenizer (v1, codec id 8) and
// TTIORefDiff (v1, codec id 9) impl files removed. Reader paths
// rejected with NSError; v2 codec headers used for the surviving
// reader dispatch.
#import "Codecs/TTIOFqzcompNx16Z.h"        // M94.Z v1.2
#import "Codecs/TTIODeltaRans.h"           // M95 v1.2
#import "Codecs/TTIOReferenceResolver.h"  // M93 v1.2
#import "Codecs/TTIOMateInfoV2.h"          // inline mate-pair codec
#import "Codecs/TTIORefDiffV2.h"          // bit-packed ref-diff v2
#import "Codecs/TTIONameTokenizerV2.h"     // v1.8 #11 ch3: adaptive name-tokenizer v2
#import "Codecs/Registry/TTIOCodecRegistry.h"  // Task 5: codec registry dispatch
#import "Codecs/Registry/TTIOCodec.h"
#import "Codecs/Registry/TTIOChannelPayload.h"
#import "Codecs/Registry/TTIODecodedChannel.h"
#import "Codecs/Registry/TTIOCodecContext.h"
#import <hdf5.h>
#include <stdlib.h>

/* How many blocks decode ahead of a serial consumer. Each one in
 * flight stays resident until it is consumed, so this is a memory
 * setting as much as a latency one, and the two do not pull the same
 * way: a consumer slower than the decoder needs almost no lookahead,
 * and paying for more of it costs throughput rather than buying any.
 * TTIO_READ_AHEAD_BLOCKS exists so the trade can be measured rather
 * than argued; see TtioGenomicReadBench. */
static NSUInteger TTIOReadAheadBlocks(void)
{
    const char *env = getenv("TTIO_READ_AHEAD_BLOCKS");
    if (env && env[0]) {
        long v = strtol(env, NULL, 10);
        if (v > 0) return (NSUInteger)v;
    }
    return 4;
}

/** One prefetched block view in flight. */
@interface TTIOInFlightView : NSObject
@property (nonatomic, strong, nullable) TTIOGenomicRun *view;
@property (nonatomic, strong, nullable) TTIOBlockView *handle;
@property (nonatomic, strong, nullable) NSError *error;
@property (nonatomic) BOOL done;
@end

@implementation TTIOInFlightView
@end

@implementation TTIOGenomicRun {
    TTIOGenomicIndex *_index;
    id<TTIOStorageGroup> _group;
    id<TTIOStorageGroup> _signalChannelsGroup;       // lazily opened, cached
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *_signalCache;
    NSMutableDictionary<NSString *, NSArray *> *_compoundCache;
    // lazy whole-channel decode cache for byte channels whose
    // @compression attribute names a TTIO codec (rANS / BASE_PACK).
    // Codec output is byte-stream non-sliceable, so the whole channel
    // is decoded once on first access and the decoded buffer is
    // sliced from memory thereafter (Binding Decision §89). Cache
    // lifetime is the TTIOGenomicRun instance — re-opening the file
    // incurs the decode cost again (Gotcha §101).
    NSMutableDictionary<NSString *, NSData *> *_decodedByteChannels;
    // lazy whole-list decode cache for read_names when
    // it's stored as a flat 1-D uint8 dataset with @compression == 8
    // (NAME_TOKENIZED). Held as NSArray<NSString *> rather than NSData
    // because the codec returns a list of names indexed by read
    // number — separate from _decodedByteChannels per Binding
    // Decision §114. Cache lifetime is the TTIOGenomicRun instance
    // (Gotcha §125 — re-opening the file incurs the decode cost
    // again; for very large runs the decoded list is materialised in
    // RAM in one shot since the codec is per-batch, Gotcha §124).
    NSArray<NSString *> *_decodedReadNames;
    // lazy whole-list decode cache for cigars when it's
    // stored as a flat 1-D uint8 dataset with @compression in
    // {RANS_ORDER0 (4), RANS_ORDER1 (5), NAME_TOKENIZED (8)}. Held as
    // NSArray<NSString *> because all three codec paths return a
    // list of CIGAR strings indexed by read number — separate from
    // _decodedReadNames per Binding Decision §123 since the two
    // channels have independent dispatch shapes (rANS uses length-
    // prefix-concat, NAME_TOKENIZED uses its own self-describing
    // wire format). Cache lifetime is the TTIOGenomicRun instance
    // (Gotcha §138 — re-opening the file incurs the decode cost
    // again).
    NSArray<NSString *> *_decodedCigars;
    // M101: decoded SAM_TAGS channel, one string per read; @[] when
    // the run has no tags channel. nil until first access.
    NSArray<NSString *> *_decodedTags;
    // lazy whole-channel decode cache for integer
    // channels (positions / flags / mapping_qualities) whose
    // @compression attribute names a TTIO rANS id. Held as NSData
    // (LE byte representation of the original integer array) keyed
    // by channel name. Separate from _decodedByteChannels per
    // Binding Decision §116 — the codec output is interpreted as
    // typed integers via channel-name dtype lookup (§115), not as
    // a uint8 byte stream. Cache lifetime is the TTIOGenomicRun
    // instance.
    // v1.6 (L4): _decodedIntChannels removed (cache for the dropped
    // intChannelArrayNamed: helper).
    // combined per-field cache for the mate_info subgroup
    // (Binding Decision §129). Single NSMutableDictionary keyed by the
    // on-disk child name (@"chrom", @"pos", @"tlen") since the three
    // fields have three different value types — chrom is
    // NSArray<NSString *>, pos is NSData carrying int64 LE bytes, tlen
    // is NSData carrying int32 LE bytes. Separate from
    // _decodedByteChannels / _decodedIntChannels / _decodedReadNames /
    // _decodedCigars per Binding Decision §129. Cache lifetime is the
    // TTIOGenomicRun instance (re-opening the file incurs the decode
    // cost again).
    NSMutableDictionary<NSString *, id> *_decodedMateInfo;
    // cached link-type query result for
    // signal_channels/mate_info. -1 = not yet probed; 0 = M82 compound
    // dataset; 1 = Phase F subgroup. Probed once via H5Oget_info_by_name
    // on first mate-field access (Binding Decision §128, Gotcha §141).
    int8_t _mateInfoLinkType;
    // cached link-type query result for
    // signal_channels/sequences. -1 = not yet probed; 0 = flat dataset
    // (v1 REF_DIFF / BASE_PACK / rANS / uncompressed); 1 = GROUP (v1.8
    // refdiff_v2 layout). Probed once on first sequences access.
    int8_t _sequencesLinkType;
    // decoded flat sequence bytes from the refdiff_v2 blob.
    // Populated on first access when _sequencesLinkType == 1.
    NSData *_decodedRefDiffV2Sequences;
    // Task 5: lazily-built, run-instance-cached codec context passed to
    // the TTIOCodecRegistry decode adapters. Mirrors the field
    // derivations the bespoke decode methods built inline.
    TTIOCodecContext *_codecCtxCache;
    // blocks_v1 (format-spec 10.12): the block table, the last
    // materialised block view, and the reference resolver shared with
    // the views. _layout is "whole" for the v1.8 whole-channel layout.
    NSString *_layout;
    TTIOBlockTable *_blockTable;
    TTIOGenomicRun *_cachedView;
    TTIOBlockView *_cachedHandle;
    NSUInteger _cachedBlock;
    TTIOReferenceResolver *_injectedResolver;
    TTIOReferenceResolver *_viewResolver;
    BOOL _viewResolverBuilt;
    NSArray<NSString *> *_chromNamesTable;
    NSArray<NSString *> *_mateChromNamesTable;
    // blocks_v1_grouped (M103): row j's input index, and its inverse
    // (input index -> stored row), uint32 each. nil for every other
    // layout. _storedIndex is the on-disk (stored-order) index, loaded
    // on first use by the bulk accessors.
    NSData *_inputIndex;
    NSData *_storedOf;
    TTIOGenomicIndex *_storedIndex;
}

@synthesize index = _index;

- (NSUInteger)readCount
{
    return _blockTable ? (NSUInteger)_blockTable.readCount : [self index].count;
}

- (TTIOGenomicIndex *)index
{
    if (_index == nil) {
        if (_storedOf != nil) {
            // M103: the index in input order, offsets recomputed.
            TTIOGenomicIndex *stored = [self _storedIndex];
            _index = stored ? [stored indexPermutedBy:_storedOf] : nil;
        } else {
            id<TTIOStorageGroup> ig = [_group openGroupNamed:@"genomic_index" error:NULL];
            if (ig) _index = [TTIOGenomicIndex readFromGroup:ig error:NULL];
        }
    }
    return _index;
}

/* The on-disk index of a grouped run, in stored order. */
- (TTIOGenomicIndex *)_storedIndex
{
    if (_storedIndex == nil) {
        id<TTIOStorageGroup> ig = [_group openGroupNamed:@"genomic_index" error:NULL];
        if (ig) _storedIndex = [TTIOGenomicIndex readFromGroup:ig error:NULL];
    }
    return _storedIndex;
}

- (NSString *)layout { return _layout ?: @"whole"; }

- (NSData *)inputIndex { return _inputIndex; }

/* The stored row of input read i (i itself for an ungrouped run, and
 * for an out-of-range i, which the block lookup then rejects). */
- (NSUInteger)_storedRow:(NSUInteger)i
{
    if (_storedOf == nil || i >= _storedOf.length / sizeof(uint32_t)) return i;
    return ((const uint32_t *)_storedOf.bytes)[i];
}

- (NSUInteger)blockCount { return _blockTable ? _blockTable.count : 1; }

- (NSArray<NSString *> *)chromosomeNames
{
    if (_chromNamesTable == nil) {
        id<TTIOStorageGroup> ig = [_group openGroupNamed:@"genomic_index" error:NULL];
        _chromNamesTable = ig ? [TTIOBlockView readNamesIn:ig named:@"chromosome_names"] : @[];
    }
    return _chromNamesTable;
}

- (void)close
{
    [self _dropCachedView];
}

- (void)dealloc
{
    [self _dropCachedView];
}

// ── blocks_v1 dispatch ─────────────────────────────────────────

/* The TTIOGenomicRun over block b, materialised on demand; the last
 * one is cached. */
/* The view handle for block b, materialised on the caller's thread
 * (storage reads). */
- (TTIOBlockView *)_materialiseHandle:(NSUInteger)b error:(NSError **)error
{
    if (_mateChromNamesTable == nil) {
        id<TTIOStorageGroup> sc = [self signalChannelsGroupWithError:NULL];
        id<TTIOStorageGroup> mate = (sc && [sc hasChildNamed:@"mate_info"])
            ? [sc openGroupNamed:@"mate_info" error:NULL] : nil;
        _mateChromNamesTable = mate ? [TTIOBlockView readNamesIn:mate named:@"chrom_names"] : @[];
    }
    return [TTIOBlockView materialiseBlock:b ofRun:_group table:_blockTable
                                chromNames:[self chromosomeNames]
                            mateChromNames:_mateChromNamesTable
                                     error:error];
}

- (TTIOGenomicRun *)_blockView:(NSUInteger)b error:(NSError **)error
{
    if (_cachedView != nil && _cachedBlock == b) return _cachedView;
    TTIOBlockView *h = [self _materialiseHandle:b error:error];
    if (!h) return nil;
    TTIOGenomicRun *sub = [TTIOGenomicRun openFromGroup:h.group name:_name
                                      referenceResolver:[self _resolverForViews] error:error];
    if (!sub) { [h discard]; return nil; }
    [self _dropCachedView];
    _cachedBlock = b;
    _cachedView = sub;
    _cachedHandle = h;
    return sub;
}

- (void)_dropCachedView
{
    _cachedView = nil;
    if (_cachedHandle) { [_cachedHandle discard]; _cachedHandle = nil; }
    _cachedBlock = NSNotFound;
}

- (TTIOReferenceResolver *)_resolverForViews
{
    if (_injectedResolver != nil) return _injectedResolver;
    // Block-decode workers call this concurrently: without the lock a
    // second caller could see the built flag before the resolver is
    // assigned and decode a REF_DIFF view with no resolver.
    @synchronized (self) {
        if (!_viewResolverBuilt) {
            _viewResolver = [self _resolverFromOwningFile];
            _viewResolverBuilt = YES;
        }
        return _viewResolver;
    }
}

/* Built exactly as the REF_DIFF_V2 decode path builds it from the run
 * group's HDF5 root; nil on non-HDF5 backends. */
- (TTIOReferenceResolver *)_resolverFromOwningFile
{
    if (![_group respondsToSelector:@selector(unwrap)]) return nil;
    TTIOHDF5Group *runG = [(id)_group performSelector:@selector(unwrap)];
    hid_t fid = H5Iget_file_id([runG groupId]);
    if (fid < 0) return nil;
    TTIOHDF5Group *rootHDF5 = nil;
    hid_t rootId = H5Gopen2(fid, "/", H5P_DEFAULT);
    if (rootId >= 0) {
        rootHDF5 = [[TTIOHDF5Group alloc] initWithGroupId:rootId retainer:nil];
    }
    H5Idec_ref(fid);
    if (rootHDF5 == nil) return nil;
    return [[TTIOReferenceResolver alloc] initWithRootGroup:rootHDF5
                                      externalReferencePath:nil];
}

- (BOOL)_blockRange:(NSUInteger)start
                 to:(NSUInteger)stop
              error:(NSError **)error
         usingBlock:(BOOL (^)(TTIOGenomicRun *view, NSUInteger blockIndex, NSUInteger localStart, NSUInteger localStop))body
{
    NSUInteger i = start;
    while (i < stop) {
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) {
            if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:0
                userInfo:@{NSLocalizedDescriptionKey: [NSString stringWithFormat:
                    @"index %lu out of range [0, %llu)", (unsigned long)i, _blockTable.readCount]}];
            return NO;
        }
        NSUInteger r0 = (NSUInteger)[_blockTable readStartAt:b];
        NSUInteger bEnd = r0 + [_blockTable nReadsAt:b];
        NSUInteger segEnd = MIN(stop, bEnd);
        TTIOGenomicRun *view = [self _blockView:b error:error];
        if (!view) return NO;
        if (!body(view, b, i - r0, segEnd - r0)) return NO;
        i = segEnd;
    }
    return YES;
}

- (BOOL)iterReadsFrom:(NSUInteger)start
                   to:(NSUInteger)stop
              threads:(NSUInteger)threads
                error:(NSError **)error
           usingBlock:(void (^)(TTIOAlignedRead *read, NSUInteger index, BOOL *stop))block
{
    NSUInteger n = [self readCount];
    NSUInteger hi = MIN(stop, n);
    NSUInteger nthreads = threads ? threads : [TTIOThreads resolve:nil];
    if (!_blockTable || nthreads <= 1 || start >= hi || _storedOf != nil) {
        return [self iterReadsFrom:start to:hi error:error usingBlock:block];
    }
    NSUInteger window = MIN(nthreads, TTIOReadAheadBlocks());
    /* The pool takes the window, not nthreads: it sizes V6's segment
     * threads from the number of workers, and the window is how many
     * blocks are ever in flight. */
    TTIOThreadPool *pool = [TTIOThreadPool poolWithThreads:window];
    NSUInteger bFirst = [_blockTable blockForRead:start];
    NSUInteger bLast = [_blockTable blockForRead:hi - 1];
    if (bFirst == NSNotFound || bLast == NSNotFound) {
        [pool close];
        return [self iterReadsFrom:start to:hi error:error usingBlock:block];
    }
    NSMutableDictionary<NSNumber *, TTIOInFlightView *> *pending = [NSMutableDictionary dictionary];
    NSCondition *cond = [NSCondition new];
    __weak typeof(self) weakSelf = self;
    void (^submit)(NSUInteger) = ^(NSUInteger b) {
        typeof(self) sself = weakSelf;
        if (!sself || b > bLast || pending[@(b)]) return;
        NSError *me = nil;
        TTIOBlockView *h = [sself _materialiseHandle:b error:&me];   /* storage reads, this thread */
        TTIOInFlightView *f = [TTIOInFlightView new];
        pending[@(b)] = f;
        if (!h) { f.error = me; f.done = YES; return; }
        [pool.queue addOperationWithBlock:^{
            NSError *ie = nil;
            TTIOGenomicRun *v = nil;
            @try {
                v = [TTIOGenomicRun openFromGroup:h.group name:sself->_name
                                referenceResolver:[sself _resolverForViews] error:&ie];
                if (v && [v readCount] > 0) [v readAtIndex:0 error:&ie];  /* warm every channel cache */
            } @catch (NSException *ex) {
                v = nil;
                ie = TTIOMakeError(TTIOErrorDatasetRead, @"%@: %@", ex.name, ex.reason);
            }
            [cond lock];
            f.view = v;
            f.handle = h;
            f.error = ie;
            f.done = YES;
            [cond broadcast];
            [cond unlock];
        }];
    };
    for (NSUInteger b = bFirst; b <= MIN(bLast, bFirst + window - 1); b++) submit(b);
    BOOL halted = NO, ok = YES;
    NSUInteger i = start;
    for (NSUInteger b = bFirst; b <= bLast && i < hi && !halted; b++) {
        TTIOInFlightView *f = pending[@(b)];
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [pending removeObjectForKey:@(b)];
        if (!f.view) {
            if (error) *error = f.error;
            ok = NO;
            [f.handle discard];
            break;
        }
        submit(b + window);
        NSUInteger r0 = (NSUInteger)[_blockTable readStartAt:b];
        NSUInteger bEnd = MIN(r0 + [_blockTable nReadsAt:b], hi);
        for (NSUInteger j = i; j < bEnd && !halted; j++) {
            NSError *re = nil;
            TTIOAlignedRead *r = [f.view readAtIndex:j - r0 error:&re];
            if (!r) {
                if (error) *error = re;
                ok = NO;
                halted = YES;
                break;
            }
            block(r, j, &halted);
        }
        [f.handle discard];
        i = bEnd;
    }
    for (TTIOInFlightView *f in pending.allValues) {
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [f.handle discard];
    }
    [pool close];
    return ok;
}

/* How many blocks may be in flight at once, given how many threads want
 * one each. A block in flight is resident for as long as the caller is
 * inside it, so this is a memory setting before it is a concurrency
 * one, and the thread count is only ever an upper bound on it.
 *
 * The cost of a block in flight is taken as eight times its sequence
 * bytes, the same figure +resolveMemoryBudget: uses on the writing
 * side, so TTIO_MEMORY_BUDGET and the half-of-physical cap mean the
 * same thing to a reader as to a writer. Measured here it is nearer
 * ten (about 100 MB resident per 10 MB block, one process per window),
 * so the estimate is a little optimistic; the thread count is the
 * binding limit on any machine with memory to spare, and the budget
 * only takes over when there is not.
 *
 * The largest block in the range sets the size rather than the mean:
 * one outsized block in flight costs its own bytes, not the average. */
- (NSUInteger)_blockWindowForThreads:(NSUInteger)nthreads
                               first:(NSUInteger)bFirst
                                last:(NSUInteger)bLast
{
    unsigned long long widest = 0;
    for (NSUInteger b = bFirst; b <= bLast; b++) {
        unsigned long long nb = (unsigned long long)[_blockTable nBasesAt:b];
        if (nb > widest) widest = nb;
    }
    if (widest == 0) return 1;
    unsigned long long budget =
        [TTIOThreads resolveMemoryBudget:nil threads:nthreads blockBytes:widest];
    unsigned long long admits = budget / (widest * 8ull);
    if (admits < 1) admits = 1;
    NSUInteger blocks = bLast - bFirst + 1;
    NSUInteger window = (NSUInteger)MIN((unsigned long long)nthreads, admits);
    /* More workers than blocks buys nothing: one block, one thread. */
    if (window > blocks) window = blocks;
    return window < 1 ? 1 : window;
}

/* Hands each decoded block to the caller on the pool rather
 * than iterating its reads on the caller's thread, because the serial
 * consumer is what bounds -iterReadsFrom:to:threads:, not the decode
 * that already runs ahead of it.
 *
 * The caller's block runs on several threads at once and in no
 * particular order, so it has to be safe to call that way; that is the
 * whole difference from -iterReadsFrom:to:threads:, whose ordering
 * guarantee is what costs the parallelism. Handles are still
 * materialised on this thread, since those are storage reads. */
- (BOOL)iterBlocksFrom:(NSUInteger)start
                    to:(NSUInteger)stop
               threads:(NSUInteger)threads
                 error:(NSError **)error
            usingBlock:(void (^)(TTIOGenomicRun *view, NSUInteger viewStart,
                                 NSUInteger firstRead, NSUInteger nReads,
                                 BOOL *stop))userBlock
{
    NSUInteger n = [self readCount];
    NSUInteger hi = MIN(stop, n);
    NSUInteger nthreads = threads ? threads : [TTIOThreads resolve:nil];
    if (start >= hi) return YES;
    /* A run with no block table is one block as far as a caller of this
     * method is concerned: hand it over whole, on this thread. */
    if (!_blockTable) {
        BOOL s = NO;
        userBlock(self, start, start, hi - start, &s);
        return YES;
    }
    NSUInteger bFirst = [_blockTable blockForRead:start];
    NSUInteger bLast = [_blockTable blockForRead:hi - 1];
    if (bFirst == NSNotFound || bLast == NSNotFound) return NO;

    /* Blocks in flight is a memory setting before it is a concurrency
     * one: each one stays resident for as long as the caller is inside
     * it. The budget resolver is the writer's, so TTIO_MEMORY_BUDGET
     * and the half-of-physical cap mean the same thing on both sides. */
    NSUInteger window = [self _blockWindowForThreads:nthreads
                                              first:bFirst
                                               last:bLast];
    /* One worker still decodes ahead. Consuming and decoding the same
     * block one after the other is what -iterReadsFrom:to:threads:
     * never does, and doing it here measured a third of that method's
     * rate: the loss is the lost overlap, not the delivery shape. So a
     * single consumer keeps a decode running one block ahead of it and
     * calls back on this thread, in order. */
    /* Below the decode-ahead depth, one consumer fed by that many
     * decoders beats the same number of threads each decoding and then
     * consuming its own block. Decoding is the slower half — one
     * decoder measured about a third of what one consumer gets through
     * — so until there are enough threads to cover it, pipelining wins
     * over dividing. Measured at 15 blocks: 1 thread 77.6 MB/s this
     * way against 23.4 the other, and 2 threads 43.6 the other way. */
    NSUInteger blocks = bLast - bFirst + 1;
    NSUInteger aheadN = MIN(blocks, MAX((NSUInteger)2, TTIOReadAheadBlocks()));
    if (window < aheadN) {
        aheadN = MIN(aheadN, MAX(window, (NSUInteger)1) * 4);
        TTIOThreadPool *ahead = [TTIOThreadPool poolWithThreads:aheadN];
        NSMutableDictionary<NSNumber *, TTIOInFlightView *> *q = [NSMutableDictionary dictionary];
        NSCondition *c = [NSCondition new];
        BOOL ok = YES, halted1 = NO;
        void (^decode)(NSUInteger) = ^(NSUInteger b) {
            if (b > bLast || q[@(b)]) return;
            NSError *me = nil;
            TTIOBlockView *h = [self _materialiseHandle:b error:&me];
            TTIOInFlightView *f = [TTIOInFlightView new];
            q[@(b)] = f;
            if (!h) { f.error = me; f.done = YES; return; }
            [ahead.queue addOperationWithBlock:^{
                NSError *ie = nil;
                TTIOGenomicRun *v = [TTIOGenomicRun openFromGroup:h.group name:self->_name
                                                referenceResolver:[self _resolverForViews]
                                                            error:&ie];
                /* openFromGroup is lazy: without this the channel
                 * decode lands on the consumer's thread at its first
                 * read, which is the whole of the overlap. */
                if (v && [v readCount] > 0) [v readAtIndex:0 error:&ie];
                [c lock];
                f.view = v; f.handle = h; f.error = ie; f.done = YES;
                [c broadcast];
                [c unlock];
            }];
        };
        for (NSUInteger b = bFirst; b <= MIN(bLast, bFirst + aheadN - 1); b++) decode(b);
        for (NSUInteger b = bFirst; b <= bLast && !halted1; b++) {
            TTIOInFlightView *f = q[@(b)];
            [c lock];
            while (!f.done) [c wait];
            [c unlock];
            [q removeObjectForKey:@(b)];
            if (!f.view) { if (error) *error = f.error; ok = NO; [f.handle discard]; break; }
            decode(b + aheadN);
            NSUInteger r0 = (NSUInteger)[_blockTable readStartAt:b];
            NSUInteger nr = (NSUInteger)[_blockTable nReadsAt:b];
            NSUInteger from = MAX(r0, start), to = MIN(r0 + nr, hi);
            if (to > from) userBlock(f.view, from - r0, from, to - from, &halted1);
            [f.handle discard];
        }
        for (TTIOInFlightView *f in q.allValues) {
            [c lock];
            while (!f.done) [c wait];
            [c unlock];
            [f.handle discard];
        }
        [ahead close];
        return ok;
    }
    TTIOThreadPool *pool = [TTIOThreadPool poolWithThreads:window];

    NSMutableDictionary<NSNumber *, TTIOInFlightView *> *pending = [NSMutableDictionary dictionary];
    NSCondition *cond = [NSCondition new];
    __block BOOL halted = NO;
    __weak typeof(self) weakSelf = self;
    void (^submit)(NSUInteger) = ^(NSUInteger b) {
        typeof(self) sself = weakSelf;
        if (!sself || b > bLast || pending[@(b)] || halted) return;
        NSError *me = nil;
        /* Storage reads, so this thread. Measured at 9 to 11 per cent
         * of the wall time, which is not what bounds this method. */
        TTIOBlockView *h = [sself _materialiseHandle:b error:&me];
        TTIOInFlightView *f = [TTIOInFlightView new];
        pending[@(b)] = f;
        if (!h) { f.error = me; f.done = YES; return; }
        NSUInteger r0 = (NSUInteger)[sself->_blockTable readStartAt:b];
        NSUInteger nr = (NSUInteger)[sself->_blockTable nReadsAt:b];
        [pool.queue addOperationWithBlock:^{
            NSError *ie = nil;
            TTIOGenomicRun *v = nil;
            @try {
                v = [TTIOGenomicRun openFromGroup:h.group name:sself->_name
                                referenceResolver:[sself _resolverForViews] error:&ie];
                if (v && !halted) {
                    BOOL s = NO;
                    NSUInteger from = MAX(r0, start);
                    NSUInteger to = MIN(r0 + nr, hi);
                    if (to > from) userBlock(v, from - r0, from, to - from, &s);
                    if (s) halted = YES;
                }
            } @catch (NSException *ex) {
                v = nil;
                ie = TTIOMakeError(TTIOErrorDatasetRead, @"%@: %@", ex.name, ex.reason);
            }
            [cond lock];
            f.view = v; f.handle = h; f.error = ie; f.done = YES;
            [cond broadcast];
            [cond unlock];
        }];
    };
    for (NSUInteger b = bFirst; b <= MIN(bLast, bFirst + window - 1); b++) submit(b);
    BOOL ok = YES;
    for (NSUInteger b = bFirst; b <= bLast && !halted; b++) {
        TTIOInFlightView *f = pending[@(b)];
        if (!f) continue;
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [pending removeObjectForKey:@(b)];
        if (!f.view) { if (error) *error = f.error; ok = NO; [f.handle discard]; break; }
        submit(b + window);
        [f.handle discard];
    }
    for (TTIOInFlightView *f in pending.allValues) {
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [f.handle discard];
    }
    [pool close];
    return ok;
}

- (BOOL)iterReadsFrom:(NSUInteger)start
                   to:(NSUInteger)stop
                error:(NSError **)error
           usingBlock:(void (^)(TTIOAlignedRead *read, NSUInteger index, BOOL *stop))block
{
    NSUInteger n = [self readCount];
    stop = MIN(stop, n);
    __block BOOL halted = NO;
    if (_storedOf != nil) return [self _iterGroupedFrom:start to:stop error:error usingBlock:block];
    if (_blockTable) {
        __block NSError *inner = nil;
        TTIOBlockTable *table = _blockTable;
        BOOL ok = [self _blockRange:start to:stop error:error
                         usingBlock:^BOOL(TTIOGenomicRun *view, NSUInteger b, NSUInteger ls, NSUInteger le) {
            NSUInteger r0 = (NSUInteger)[table readStartAt:b];
            for (NSUInteger j = ls; j < le && !halted; j++) {
                TTIOAlignedRead *r = [view readAtIndex:j error:&inner];
                if (!r) return NO;
                block(r, r0 + j, &halted);
            }
            return YES;
        }];
        if (!ok && inner && error && *error == nil) *error = inner;
        return ok;
    }
    for (NSUInteger i = start; i < stop && !halted; i++) {
        TTIOAlignedRead *r = [self readAtIndex:i error:error];
        if (!r) return NO;
        block(r, i, &halted);
    }
    return YES;
}

/* Input-order reads gathered per pass of a grouped run (M103). */
static const NSUInteger kGroupedChunk = 1u << 22;

static int _ttio_cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

/* Reads [start, stop) of a grouped run in input order. Each pass gathers
 * up to kGroupedChunk input-order reads: their stored rows are visited
 * block by block, each block decoded once per pass, and the reads are
 * delivered in input order. A pass holds its reads in memory; a range of
 * a few blocks' worth decodes only those blocks. */
- (BOOL)_iterGroupedFrom:(NSUInteger)start
                      to:(NSUInteger)stop
                   error:(NSError **)error
              usingBlock:(void (^)(TTIOAlignedRead *read, NSUInteger index, BOOL *stop))block
{
    const uint32_t *so = (const uint32_t *)_storedOf.bytes;
    BOOL halted = NO;
    for (NSUInteger c0 = start; c0 < stop && !halted; c0 += kGroupedChunk) {
        @autoreleasepool {
            NSUInteger c1 = MIN(stop, c0 + kGroupedChunk);
            NSUInteger m = c1 - c0;
            // (stored row << 32 | k), sorted: rows ascending, ties by k.
            NSMutableData *keys = [NSMutableData dataWithLength:m * sizeof(uint64_t)];
            uint64_t *kv = (uint64_t *)keys.mutableBytes;
            for (NSUInteger k = 0; k < m; k++) kv[k] = ((uint64_t)so[c0 + k] << 32) | (uint64_t)k;
            qsort(kv, m, sizeof(uint64_t), _ttio_cmp_u64);
            NSMutableArray *out = [NSMutableArray arrayWithCapacity:m];
            for (NSUInteger k = 0; k < m; k++) [out addObject:[NSNull null]];
            NSUInteger k = 0;
            while (k < m) {
                NSUInteger s = (NSUInteger)(kv[k] >> 32);
                NSUInteger b = [_blockTable blockForRead:s];
                if (b == NSNotFound) {
                    if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                        @"stored row %lu out of range", (unsigned long)s);
                    return NO;
                }
                NSUInteger r0 = (NSUInteger)[_blockTable readStartAt:b];
                NSUInteger r1 = r0 + (NSUInteger)[_blockTable nReadsAt:b];
                TTIOGenomicRun *view = [self _blockView:b error:error];
                if (!view) return NO;
                while (k < m && (NSUInteger)(kv[k] >> 32) < r1) {
                    NSUInteger row = (NSUInteger)(kv[k] >> 32);
                    TTIOAlignedRead *r = [view readAtIndex:row - r0 error:error];
                    if (!r) return NO;
                    out[(NSUInteger)(kv[k] & 0xFFFFFFFFu)] = r;
                    k++;
                }
            }
            for (NSUInteger j = 0; j < m && !halted; j++) block(out[j], c0 + j, &halted);
        }
    }
    return YES;
}

- (instancetype)initWithName:(NSString *)name
              acquisitionMode:(TTIOAcquisitionMode)mode
                     modality:(NSString *)modality
                 referenceUri:(NSString *)refUri
                     platform:(NSString *)platform
                   sampleName:(NSString *)sampleName
                        index:(TTIOGenomicIndex *)index
                        group:(id<TTIOStorageGroup>)group
{
    self = [super init];
    if (self) {
        _name             = [name copy];
        _acquisitionMode  = mode;
        _modality         = [modality copy];
        _referenceUri     = [refUri copy];
        _platform         = [platform copy];
        _sampleName       = [sampleName copy];
        _index            = index;
        _group            = group;
        _signalCache      = [NSMutableDictionary dictionary];
        _compoundCache    = [NSMutableDictionary dictionary];
        _decodedByteChannels = [NSMutableDictionary dictionary];
        _decodedMateInfo     = [NSMutableDictionary dictionary];
        _mateInfoLinkType    = -1;  // not yet probed
        _sequencesLinkType   = -1;  // not yet probed
        _cachedBlock         = NSNotFound;
        _layout              = @"whole";
    }
    return self;
}

+ (instancetype)openFromGroup:(id<TTIOStorageGroup>)runGroup
                          name:(NSString *)name
                         error:(NSError **)error
{
    return [self openFromGroup:runGroup name:name referenceResolver:nil error:error];
}

+ (instancetype)openFromGroup:(id<TTIOStorageGroup>)runGroup
                          name:(NSString *)name
             referenceResolver:(id)resolver
                         error:(NSError **)error
{
    if (!runGroup) return nil;

    NSString *layout = @"whole";
    if ([runGroup hasAttributeNamed:@"layout"]) {
        id lv = [runGroup attributeValueForName:@"layout" error:NULL];
        if (lv) layout = [lv description];
    }
    TTIOGenomicIndex *index = nil;
    TTIOBlockTable *table = nil;
    NSData *inputIndex = nil, *storedOf = nil;
    if ([layout isEqualToString:@"blocks_v1"] || [layout isEqualToString:@"blocks_v1_grouped"]) {
        table = [TTIOBlockTable readFromRunGroup:runGroup error:error];
        if (!table) return nil;
        if ([layout isEqualToString:@"blocks_v1_grouped"]) {
            /* M103: reads are stored grouped by sequence; input_index
             * restores input order (format-spec 10.12.7). */
            inputIndex = [self _readInputIndexOf:runGroup name:name
                                       readCount:(NSUInteger)table.readCount
                                        storedOf:&storedOf error:error];
            if (!inputIndex) return nil;
        }
    } else if ([layout isEqualToString:@"whole"]) {
        id<TTIOStorageGroup> idxGroup = [runGroup openGroupNamed:@"genomic_index" error:error];
        if (!idxGroup) return nil;
        index = [TTIOGenomicIndex readFromGroup:idxGroup error:error];
        if (!index) return nil;
    } else {
        if (error) *error = TTIOMakeError(TTIOErrorUnsupportedLayout,
            @"genomic run '%@': unsupported layout '%@' (this reader knows the "
            @"whole-channel layout, blocks_v1 and blocks_v1_grouped)", name, layout);
        return nil;
    }

    // Integer attribute: the storage-protocol adapter tries
    // stringAttributeNamed first which silently returns garbage bytes
    // for INT64 attrs (TTIOHDF5Group.stringAttributeNamed doesn't
    // type-check). Read directly via the underlying HDF5Group when
    // available; fall back to the protocol for non-HDF5 providers.
    int64_t modeValue = 0;
    if ([runGroup respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *h5 = [(id)runGroup performSelector:@selector(unwrap)];
        BOOL exists = NO;
        modeValue = [h5 integerAttributeNamed:@"acquisition_mode"
                                        exists:&exists error:NULL];
    } else {
        id v = [runGroup attributeValueForName:@"acquisition_mode" error:NULL];
        if ([v isKindOfClass:[NSNumber class]]) modeValue = [v longLongValue];
    }
    NSString *modality  = [runGroup attributeValueForName:@"modality"         error:error];
    NSString *refUri    = [runGroup attributeValueForName:@"reference_uri"    error:error];
    NSString *platform  = [runGroup attributeValueForName:@"platform"         error:error];
    NSString *sampleN   = [runGroup attributeValueForName:@"sample_name"      error:error];
    // @read_role is absent on pre-M97 files; never surface that as an
    // open error.
    NSString *readRole  = [runGroup attributeValueForName:@"read_role"        error:NULL];

    TTIOGenomicRun *run = [[TTIOGenomicRun alloc]
        initWithName:name
     acquisitionMode:(TTIOAcquisitionMode)modeValue
            modality:modality ?: @"genomic_sequencing"
        referenceUri:refUri ?: @""
            platform:platform ?: @""
          sampleName:sampleN ?: @""
               index:index
               group:runGroup];
    run->_layout = layout;
    run->_blockTable = table;
    run->_inputIndex = inputIndex;
    run->_storedOf = storedOf;
    run->_injectedResolver = resolver;
    run->_readRole = [readRole isKindOfClass:[NSString class]]
        ? [readRole copy] : nil;
    return run;
}

/* genomic_index/input_index of a blocks_v1_grouped run, checked to be a
 * permutation of its readCount reads; *storedOf receives the inverse. */
+ (NSData *)_readInputIndexOf:(id<TTIOStorageGroup>)runGroup
                         name:(NSString *)name
                    readCount:(NSUInteger)n
                     storedOf:(NSData **)storedOf
                        error:(NSError **)error
{
    id<TTIOStorageGroup> ig = [runGroup openGroupNamed:@"genomic_index" error:error];
    if (!ig) return nil;
    id<TTIOStorageDataset> ds = [ig hasChildNamed:@"input_index"]
        ? [ig openDatasetNamed:@"input_index" error:NULL] : nil;
    id raw = ds ? [ds readAll:NULL] : nil;
    BOOL ok = [raw isKindOfClass:[NSData class]]
        && [ds precision] == TTIOPrecisionUInt32
        && [(NSData *)raw length] == n * sizeof(uint32_t);
    NSMutableData *inv = ok ? [NSMutableData dataWithLength:n * sizeof(uint32_t)] : nil;
    if (ok) {
        const uint32_t *ii = (const uint32_t *)[(NSData *)raw bytes];
        uint32_t *so = (uint32_t *)inv.mutableBytes;
        NSMutableData *seen = [NSMutableData dataWithLength:n];
        uint8_t *sv = (uint8_t *)seen.mutableBytes;
        for (NSUInteger j = 0; j < n; j++) {
            uint32_t v = ii[j];
            if (v >= n || sv[v]) { ok = NO; break; }
            sv[v] = 1;
            so[v] = (uint32_t)j;
        }
    }
    if (!ok) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
            @"genomic run '%@': input_index is not a permutation of its %lu reads",
            name, (unsigned long)n);
        return nil;
    }
    *storedOf = inv;
    return [(NSData *)raw copy];
}

- (id<TTIOStorageGroup>)signalChannelsGroupWithError:(NSError **)error
{
    if (!_signalChannelsGroup) {
        _signalChannelsGroup = [_group openGroupNamed:@"signal_channels" error:error];
    }
    return _signalChannelsGroup;
}

- (id<TTIOStorageDataset>)signalDatasetNamed:(NSString *)name error:(NSError **)error
{
    id<TTIOStorageDataset> ds = _signalCache[name];
    if (!ds) {
        id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:error];
        if (!sig) return nil;
        ds = [sig openDatasetNamed:name error:error];
        if (ds) _signalCache[name] = ds;
    }
    return ds;
}

- (NSArray *)compoundRowsNamed:(NSString *)name
                         field:(TTIOCompoundField *)field
                         error:(NSError **)error
{
    NSArray *rows = _compoundCache[name];
    if (rows) return rows;

    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:error];
    if (!sig) return nil;

    NSArray *fields = field ? @[field] : nil;
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *h5 = [(id)sig performSelector:@selector(unwrap)];
        rows = [TTIOCompoundIO readGenericFromGroup:h5
                                        datasetNamed:name
                                              fields:fields
                                               error:error];
    } else {
        id<TTIOStorageDataset> ds = [sig openDatasetNamed:name error:error];
        if (ds) rows = [ds readAll:error];
    }
    if (rows) _compoundCache[name] = rows;
    return rows;
}

// read the @compression attribute (uint8) on an HDF5 dataset.
// Returns 0 (NONE) when the attribute is absent — equivalent to
// "uncompressed at the TTIO-codec layer". The dataset hid_t is taken
// from the underlying TTIOHDF5Dataset; non-HDF5 backends fall back to
// the storage protocol's attributeValueForName:error:.
static uint8_t _ttio_m86_read_compression_attr(hid_t did)
{
    if (H5Aexists(did, "compression") <= 0) return 0;
    hid_t aid = H5Aopen(did, "compression", H5P_DEFAULT);
    if (aid < 0) return 0;
    uint8_t value = 0;
    H5Aread(aid, H5T_NATIVE_UINT8, &value);
    H5Aclose(aid);
    return value;
}

// read the @compression attribute via the storage protocol (used
// for non-HDF5 backends). Returns 0 when absent or non-numeric.
static uint8_t _ttio_m86_read_compression_attr_protocol(id<TTIOStorageDataset> ds)
{
    if (![ds hasAttributeNamed:@"compression"]) return 0;
    NSError *e = nil;
    id v = [ds attributeValueForName:@"compression" error:&e];
    if ([v isKindOfClass:[NSNumber class]]) {
        return (uint8_t)[v unsignedIntegerValue];
    }
    return 0;
}

// probe @compression on a signal channel without decoding
// anything. Used by the transport writer to decide whether to re-
// encode each per-AU UINT8 slice with the same M86 codec on the
// wire. Returns 0 (NONE) when the attribute is absent or unreadable.
- (uint8_t)wireCompressionForChannel:(NSString *)name
{
    if (!name.length) return 0;
    if (_blockTable) {
        if (_blockTable.count == 0 || !_blockTable.hasCodecs) return 0;
        NSUInteger c = [_blockTable codecOf:name at:0];
        return (c == TTIOCompressionRansOrder0 || c == TTIOCompressionRansOrder1
                || c == TTIOCompressionBasePack) ? (uint8_t)c : 0;
    }
    id<TTIOStorageDataset> ds = [self signalDatasetNamed:name error:NULL];
    if (!ds) return 0;
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        TTIOHDF5Dataset *hds = [hg openDatasetNamed:name error:NULL];
        if (hds) return _ttio_m86_read_compression_attr([hds datasetId]);
    }
    return _ttio_m86_read_compression_attr_protocol(ds);
}

// Task 5: build (and cache) the TTIOCodecContext handed to the
// TTIOCodecRegistry decode adapters. Every field is derived EXACTLY as
// the bespoke decode methods (-byteChannelSliceNamed:, -cigarAtIndex:,
// -_decodeRefDiffV2Sequences:, -_decodeMateInfoInlineV2:) built it
// inline, so registry-routed decode stays byte-identical to the old
// switch / side paths. Cached for the lifetime of the run instance.
- (TTIOCodecContext *)_codecContext
{
    if (_codecCtxCache) return _codecCtxCache;

    TTIOCodecContext *ctx = [TTIOCodecContext emptyContext];
    NSUInteger n = [self index] ? [self index].count : 0;

    // readLengths / revcompFlags / totalBases (fqzcomp + refdiff).
    // revcomp = run.flags & 16 (mirrors -_ttio_m94z_decodeFqzcompNx16Z:
    // and -_decodeRefDiffV2Sequences:'s totalBases sum).
    NSMutableArray<NSNumber *> *readLengths =
        [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSNumber *> *revcompFlags =
        [NSMutableArray arrayWithCapacity:n];
    NSUInteger totalBases = 0;
    for (NSUInteger i = 0; i < n; i++) {
        NSUInteger len = [[self index] lengthAt:i];
        [readLengths addObject:@(len)];
        totalBases += len;
        uint32_t f = [[self index] flagsAt:i];
        [revcompFlags addObject:((f & 16u) != 0) ? @1 : @0];
    }
    ctx.readLengths  = readLengths;
    ctx.revcompFlags = revcompFlags;
    ctx.totalBases   = @(totalBases);
    ctx.readCount    = @(n);

    /* fqzcomp V5: the qualities decoder needs the decoded sequences
     * channel. Route through -byteChannelSliceNamed: so refdiff_v2
     * and plain layouts both work; the block fires only for
     * version-5 streams. __weak breaks the cycle through the cached
     * context. */
    __weak typeof(self) weakSelfV5 = self;
    NSUInteger seqTotal = totalBases;
    ctx.sequencesProvider = ^NSData * _Nullable {
        return [weakSelfV5 byteChannelSliceNamed:@"sequences"
                                        offset:0
                                         count:seqTotal
                                         error:NULL];
    };

    // positions int64-LE (mirrors -_decodeRefDiffV2Sequences:).
    NSMutableData *positions = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    int64_t *posPtr = (int64_t *)positions.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) posPtr[i] = [[self index] positionAt:i];
    ctx.positions = positions;

    // chromosomes (refdiff single-chrom constraint).
    NSMutableArray<NSString *> *chromosomes =
        [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        NSString *c = [[self index] chromosomeAt:i];
        [chromosomes addObject:c ?: @""];
    }
    ctx.chromosomes = chromosomes;

    // ownChromIds (encounter-order uint16-LE) + ownPositions (int64-LE)
    // + nRecords — MUST match -_decodeMateInfoInlineV2: EXACTLY.
    NSMutableData *ownChromIds =
        [NSMutableData dataWithLength:n * sizeof(uint16_t)];
    uint16_t *ownIdsPtr = (uint16_t *)ownChromIds.mutableBytes;
    // The ids are those of the mate_info/chrom_names sidecar when the
    // run has one (own chromosomes come first in it, so this equals the
    // encounter order for a whole-channel run and stays right for a
    // block whose first chromosome is not the run's first); encounter
    // order otherwise.
    NSMutableDictionary<NSString *, NSNumber *> *nameToId =
        [NSMutableDictionary dictionaryWithCapacity:32];
    NSArray<NSString *> *sidecar = [self readMateInfoChromNamesTable];
    for (NSUInteger j = 0; j < sidecar.count; j++) {
        if (nameToId[sidecar[j]] == nil) nameToId[sidecar[j]] = @(j);
    }
    for (NSUInteger i = 0; i < n; i++) {
        NSString *name = [[self index] chromosomeAt:i];
        NSNumber *existingId = nameToId[name];
        if (existingId == nil && sidecar.count > 0) {
            ownIdsPtr[i] = (uint16_t)0xFFFF;
        } else if (existingId == nil) {
            NSUInteger newId = nameToId.count;
            nameToId[name] = @(newId);
            ownIdsPtr[i] = (uint16_t)newId;
        } else {
            ownIdsPtr[i] = (uint16_t)[existingId unsignedIntegerValue];
        }
    }
    ctx.ownChromIds = ownChromIds;

    NSMutableData *ownPositions =
        [NSMutableData dataWithLength:n * sizeof(int64_t)];
    int64_t *ownPosPtr = (int64_t *)ownPositions.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) ownPosPtr[i] = [[self index] positionAt:i];
    ctx.ownPositions = ownPositions;
    ctx.nRecords = @(n);

    // cigarsProvider — lazy thunk that reproduces the cigar-list build
    // from -_decodeRefDiffV2Sequences: (trigger cigarAtIndex:0 to
    // populate _decodedCigars; otherwise materialise per-read).
    __weak TTIOGenomicRun *weakSelf = self;
    ctx.cigarsProvider = ^NSArray<NSString *> *(void) {
        TTIOGenomicRun *s = weakSelf;
        if (s == nil) return @[];
        NSUInteger cn = [s index] ? [s index].count : 0;
        if (cn > 0) {
            NSError *cigErr = nil;
            (void)[s cigarAtIndex:0 error:&cigErr];
        }
        if (s->_decodedCigars != nil) {
            return [NSArray arrayWithArray:s->_decodedCigars];
        }
        NSMutableArray<NSString *> *cigars =
            [NSMutableArray arrayWithCapacity:cn];
        for (NSUInteger i = 0; i < cn; i++) {
            NSError *cigErr = nil;
            NSString *cig = [s cigarAtIndex:i error:&cigErr];
            if (cig == nil) return nil;
            [cigars addObject:cig];
        }
        return cigars;
    };

    // referenceResolver — built EXACTLY as -_decodeRefDiffV2Sequences:
    // does from the HDF5 root group; nil on non-HDF5 backends / failure.
    TTIOHDF5Group *rootHDF5 = nil;
    if ([_group respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *runG = [(id)_group performSelector:@selector(unwrap)];
        hid_t fid = H5Iget_file_id([runG groupId]);
        if (fid >= 0) {
            hid_t rootId = H5Gopen2(fid, "/", H5P_DEFAULT);
            if (rootId >= 0) {
                rootHDF5 = [[TTIOHDF5Group alloc] initWithGroupId:rootId
                                                         retainer:nil];
            }
            H5Idec_ref(fid);
        }
    }
    if (_injectedResolver != nil) {
        ctx.referenceResolver = _injectedResolver;
    } else if (rootHDF5 != nil) {
        ctx.referenceResolver = [[TTIOReferenceResolver alloc]
                initWithRootGroup:rootHDF5
            externalReferencePath:nil];
    }

    // SAM_TAGS (M101): the references for MD/NM derivation, indexed
    // like ownChromIds; resolved once per decode and only for
    // REF_DIFF_V2 runs. nil = the reference could not be resolved.
    NSDictionary<NSString *, NSNumber *> *tagNameToId = [nameToId copy];
    ctx.tagReferencesProvider = ^NSArray * _Nullable {
        TTIOGenomicRun *s = weakSelf;
        return s ? [s _tagReferencesForNameToId:tagNameToId] : @[];
    };

    _codecCtxCache = ctx;
    return _codecCtxCache;
}

/* References for SAM_TAGS MD/NM derivation, indexed like ownChromIds.
 * The writer derives only on REF_DIFF_V2 runs (binding decision 96),
 * whose one chromosome resolves the same way the sequences decode
 * does (the refdiff_v2 blob header's uri and md5); any other run gets
 * none. nil when the reference cannot be resolved. */
- (NSArray *)_tagReferencesForNameToId:(NSDictionary<NSString *, NSNumber *> *)nameToId
{
    if (![self _sequencesIsRefDiffV2]) return @[];
    NSMutableSet<NSString *> *unique = [NSMutableSet set];
    NSUInteger n = [self index].count;
    for (NSUInteger i = 0; i < n; i++) {
        NSString *c = [[self index] chromosomeAt:i];
        if (c.length > 0 && ![c isEqualToString:@"*"]) [unique addObject:c];
    }
    NSMutableArray *refs = [NSMutableArray arrayWithCapacity:nameToId.count];
    for (NSUInteger j = 0; j < nameToId.count; j++) [refs addObject:[NSNull null]];
    if (unique.count != 1) return refs;
    NSString *chrom = [unique anyObject];
    NSNumber *cid = nameToId[chrom];
    if (cid == nil || cid.unsignedIntegerValue >= refs.count) return refs;

    // The outer header is 38 bytes plus a uint16-length URI:
    // [0:4] "RDF2", [20:36] md5, [36:38] uri_len LE, [38:] uri.
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    id<TTIOStorageGroup> seqGrp = [sig openGroupNamed:@"sequences" error:NULL];
    id<TTIOStorageDataset> ds = [seqGrp openDatasetNamed:@"refdiff_v2" error:NULL];
    if (!ds) return nil;
    NSUInteger want = MIN(ds.length, (NSUInteger)(38 + 0xFFFF));
    id raw = [ds readSliceAtOffset:0 count:want error:NULL];
    if (![raw isKindOfClass:[NSData class]] || [(NSData *)raw length] < 38) return nil;
    const uint8_t *hb = (const uint8_t *)[(NSData *)raw bytes];
    if (memcmp(hb, "RDF2", 4) != 0) return nil;
    uint16_t uriLen = (uint16_t)(hb[36] | (hb[37] << 8));
    if ([(NSData *)raw length] < (NSUInteger)(38 + uriLen)) return nil;
    NSString *uri = [[NSString alloc] initWithBytes:hb + 38 length:uriLen
                                           encoding:NSUTF8StringEncoding] ?: @"";
    NSData *md5 = [NSData dataWithBytes:hb + 20 length:16];
    TTIOReferenceResolver *resolver = [self _codecContext].referenceResolver;
    if (resolver == nil) return nil;
    NSData *ref = [resolver resolveURI:uri expectedMD5:md5 chromosome:chrom error:NULL];
    if (ref == nil) return nil;
    refs[cid.unsignedIntegerValue] = ref;
    return refs;
}

- (BOOL)hasTagsChannel
{
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    return sig != nil && [sig hasChildNamed:@"tags"];
}

- (NSArray<NSString *> *)allTagsWithError:(NSError **)error
{
    if (_decodedTags != nil) return _decodedTags;
    if (![self hasTagsChannel]) {
        _decodedTags = @[];
        return _decodedTags;
    }
    if (_blockTable) {
        // Not cached: the per-block views hold their own decode.
        NSMutableArray<NSString *> *all =
            [NSMutableArray arrayWithCapacity:(NSUInteger)_blockTable.readCount];
        for (NSUInteger b = 0; b < _blockTable.count; b++) {
            TTIOGenomicRun *view = [self _blockView:b error:error];
            if (!view) return nil;
            NSArray<NSString *> *part = [view allTagsWithError:error];
            if (!part) return nil;
            if (part.count == 0) {
                for (NSUInteger j = 0; j < [view readCount]; j++) [all addObject:@""];
            } else {
                [all addObjectsFromArray:part];
            }
        }
        return _storedOf ? [self _inInputOrder:all] : all;
    }
    id<TTIOStorageDataset> ds = [self signalDatasetNamed:@"tags" error:error];
    if (!ds) return nil;
    uint8_t codec_id = 0;
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        TTIOHDF5Dataset *hds = [hg openDatasetNamed:@"tags" error:NULL];
        if (hds) codec_id = _ttio_m86_read_compression_attr([hds datasetId]);
    } else {
        codec_id = _ttio_m86_read_compression_attr_protocol(ds);
    }
    if (codec_id != TTIOCompressionSamTags) {
        if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:2020
            userInfo:@{NSLocalizedDescriptionKey: [NSString stringWithFormat:
                @"signal_channel 'tags': @compression=%u is not SAM_TAGS (18), "
                @"the only codec of the tags channel", (unsigned)codec_id]}];
        return nil;
    }
    id raw = [ds readAll:error];
    if (![raw isKindOfClass:[NSData class]]) return nil;
    TTIOCodecContext *ctx = [self _codecContext];
    NSArray *refs = ctx.tagReferencesProvider ? ctx.tagReferencesProvider() : @[];
    if (refs == nil) {
        if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:2021
            userInfo:@{NSLocalizedDescriptionKey:
                @"signal_channel 'tags': the reference for MD/NM derivation "
                @"could not be resolved"}];
        return nil;
    }
    // Resolve once: hand the codec a context whose provider returns the
    // references already in hand.
    TTIOCodecContext *tagCtx = [TTIOCodecContext emptyContext];
    tagCtx.readCount = ctx.readCount;
    tagCtx.readLengths = ctx.readLengths;
    tagCtx.sequencesProvider = ctx.sequencesProvider;
    tagCtx.cigarsProvider = ctx.cigarsProvider;
    tagCtx.positions = ctx.positions;
    tagCtx.ownChromIds = ctx.ownChromIds;
    tagCtx.tagReferencesProvider = ^NSArray * _Nullable { return refs; };
    id<TTIOCodec> codec = [TTIOCodecRegistry codecForId:TTIOCompressionSamTags];
    TTIODecodedChannel *dc = [codec decode:[[TTIOBytesPayload alloc] initWithBytes:raw]
                                   context:tagCtx error:error];
    if (dc == nil) return nil;
    _decodedTags = ((TTIODecodedStringList *)dc).names;
    return _decodedTags;
}

// byte-channel slice helper.
//
// For byte channels (sequences, qualities) the read path may need to
// decode through a TTIO codec when @compression > 0. We implement the
// decode-once-then-slice tradeoff (Binding Decision §89) — the whole
// channel is decoded on first access and cached on the GenomicRun
// instance. For uncompressed channels the existing per-slice
// HDF5 hyperslab read path is preserved unchanged.
//
// for the sequences channel, probe whether signal_channels/sequences
// is a GROUP (refdiff_v2 layout) and decode via TTIORefDiffV2 when true.
/* blocks_v1: a base range [offset, offset+count) of sequences or
 * qualities, concatenated over the blocks it touches. */
- (NSData *)_blockByteChannelSliceNamed:(NSString *)name
                                 offset:(NSUInteger)offset
                                  count:(NSUInteger)count
                                  error:(NSError **)error
{
    NSMutableData *out = [NSMutableData dataWithCapacity:count];
    NSUInteger nb = _blockTable.count;
    for (NSUInteger b = 0; b < nb && count > 0; b++) {
        NSUInteger bs = (NSUInteger)[_blockTable baseStartAt:b];
        NSUInteger bn = (NSUInteger)[_blockTable nBasesAt:b];
        if (bs + bn <= offset) continue;
        if (bs >= offset + count) break;
        NSUInteger from = offset > bs ? offset - bs : 0;
        NSUInteger to = MIN(bn, offset + count - bs);
        TTIOGenomicRun *view = [self _blockView:b error:error];
        if (!view) return nil;
        NSData *part = [view byteChannelSliceNamed:name offset:from count:to - from error:error];
        if (!part) return nil;
        [out appendData:part];
    }
    return out;
}

- (NSData *)byteChannelSliceNamed:(NSString *)name
                            offset:(NSUInteger)offset
                             count:(NSUInteger)count
                             error:(NSError **)error
{
    if (_storedOf != nil) {
        // M103: base coordinates are input order (the index's offsets).
        NSData *full = [self _byteChannelFullNamed:name];
        if (offset + count > full.length) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"%@ slice [%lu, %lu) past the end (%lu)", name, (unsigned long)offset,
                (unsigned long)(offset + count), (unsigned long)full.length);
            return nil;
        }
        return [full subdataWithRange:NSMakeRange(offset, count)];
    }
    if (_blockTable) return [self _blockByteChannelSliceNamed:name offset:offset count:count error:error];
    // refdiff_v2 group layout probe for sequences channel. Routed via
    // the codec registry (REF_DIFF_V2 group-payload adapter); the
    // adapter opens the sequences GROUP / refdiff_v2 dataset from the
    // signal_channels group, parses the blob, resolves the reference
    // via context.referenceResolver, and decodes — byte-identical to
    // the old -_decodeRefDiffV2Sequences: side path. Result is cached
    // in _decodedRefDiffV2Sequences exactly as before.
    if ([name isEqualToString:@"sequences"] && [self _sequencesIsRefDiffV2]) {
        NSData *decoded = _decodedRefDiffV2Sequences;
        if (!decoded) {
            id<TTIOStorageGroup> sigGrp =
                [self signalChannelsGroupWithError:error];
            if (!sigGrp) return nil;
            id<TTIOCodec> codec =
                [TTIOCodecRegistry codecForId:TTIOCompressionRefDiffV2];
            TTIODecodedChannel *dc =
                [codec decode:[[TTIOGroupPayload alloc] initWithGroup:sigGrp]
                      context:[self _codecContext] error:error];
            if (dc == nil) return nil;
            decoded = ((TTIODecodedBytes *)dc).data;
            _decodedRefDiffV2Sequences = decoded;
        }
        NSUInteger from = MIN(offset, decoded.length);
        NSUInteger to   = MIN(from + count, decoded.length);
        return [decoded subdataWithRange:NSMakeRange(from, to - from)];
    }

    NSData *cached = _decodedByteChannels[name];
    if (cached) {
        NSUInteger from = MIN(offset, cached.length);
        NSUInteger to   = MIN(from + count, cached.length);
        return [cached subdataWithRange:NSMakeRange(from, to - from)];
    }

    id<TTIOStorageDataset> ds = [self signalDatasetNamed:name error:error];
    if (!ds) return nil;

    // Detect codec via @compression on the dataset. Two paths: HDF5
    // backend exposes the underlying TTIOHDF5Dataset whose hid_t the
    // H5A* calls need; non-HDF5 backends route through the storage
    // protocol's attributeValueForName:.
    uint8_t codec_id = 0;
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        TTIOHDF5Dataset *hds = [hg openDatasetNamed:name error:NULL];
        if (hds) {
            codec_id = _ttio_m86_read_compression_attr([hds datasetId]);
        }
    } else {
        codec_id = _ttio_m86_read_compression_attr_protocol(ds);
    }

    if (codec_id == 0) {
        // No TTIO-codec dispatch — existing hyperslab path.
        id raw = [ds readSliceAtOffset:offset count:count error:error];
        if (![raw isKindOfClass:[NSData class]]) return nil;
        return (NSData *)raw;
    }

    // Codec-compressed: read all bytes, decode, cache, slice.
    id allRaw = [ds readAll:error];
    if (![allRaw isKindOfClass:[NSData class]]) return nil;
    NSData *encoded = (NSData *)allRaw;
    NSData *decoded = nil;
    NSError *decErr = nil;

    // REF_DIFF v1 (codec id 9) is unregistered and removed — keep its
    // dedicated re-encode hint instead of the generic default-arm error
    // that the registry-nil path produces below.
    if (codec_id == 9) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOGenomicRun" code:2020
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"REF_DIFF v1 (codec id 9) is no longer "
                       @"supported in v1.0; file was written with "
                       @"an older TTI-O version. Re-encode with "
                       @"v1.0+ which uses REF_DIFF_V2 (codec id "
                       @"14)."}];
        return nil;
    }

    id<TTIOCodec> codec = [TTIOCodecRegistry codecForId:(TTIOCompression)codec_id];
    if (codec == nil) {
        // Same NSError the old default arm built.
        if (error) *error = [NSError
            errorWithDomain:@"TTIOGenomicRun" code:2020
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"signal_channel '%@': @compression=%u "
                            @"is not a supported TTIO codec id",
                            name, (unsigned)codec_id]}];
        return nil;
    }
    TTIODecodedChannel *dc =
        [codec decode:[[TTIOBytesPayload alloc] initWithBytes:encoded]
              context:[self _codecContext] error:&decErr];
    if (dc != nil) decoded = ((TTIODecodedBytes *)dc).data;

    if (!decoded) {
        if (error) *error = decErr ?: [NSError
            errorWithDomain:@"TTIOGenomicRun" code:2021
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"signal_channel '%@' codec %u decode failed",
                            name, (unsigned)codec_id]}];
        return nil;
    }
    _decodedByteChannels[name] = decoded;
    NSUInteger from = MIN(offset, decoded.length);
    NSUInteger to   = MIN(from + count, decoded.length);
    return [decoded subdataWithRange:NSMakeRange(from, to - from)];
}

// _ttio_m93_decodeRefDiff (v1 REF_DIFF reader)
// removed alongside TTIORefDiff codec impl. The byte-channel codec
// dispatcher above raises a clear NSError when @compression == 9 is
// encountered on legacy files.

// read_names dispatch helper.
//
// The read_names channel has two on-disk layouts (Binding Decisions
// §111, §112):
//
//   - **M82 compound** (no override): VL_STRING-in-compound dataset,
//     read whole-and-cache via -compoundRowsNamed:.
//   - **NAME_TOKENIZED** (override active): flat 1-D uint8 dataset
//     of the same name carrying the codec output, with
//     @compression == 8. Decoded once on first access via
//     TTIONameTokenizerDecode and cached as NSArray<NSString *>
//     on this TTIOGenomicRun instance per Binding Decision §114.
//
// Dispatch is on dataset shape — a flat uint8 dataset routes through
// the codec path; otherwise fall through to the compound path.
// All call sites that touch read_names should route through this
// helper (Gotcha §126).
- (NSString *)readNameAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) {
            if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:2040
                userInfo:@{NSLocalizedDescriptionKey: [NSString stringWithFormat:
                    @"read_names index %lu out of range [0, %llu)", (unsigned long)i, _blockTable.readCount]}];
            return nil;
        }
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view readNameAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : nil;
    }
    if (_decodedReadNames != nil) {
        if (i >= _decodedReadNames.count) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2040
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"read_names index %lu out of range "
                                @"[0, %lu)",
                                (unsigned long)i,
                                (unsigned long)_decodedReadNames.count]}];
            return nil;
        }
        return _decodedReadNames[i];
    }

    id<TTIOStorageDataset> ds = [self signalDatasetNamed:@"read_names"
                                                   error:error];
    if (!ds) return nil;

    // Shape dispatch: precision == UInt8 is a flat uint8 dataset and
    // therefore the codec path; anything else is the M82 compound.
    // The HDF5 backend's -openDatasetNamed: returns precision UInt8
    // for the schema-lifted layout (TTIOHDF5Group introspects via
    // H5Tequal(H5T_NATIVE_UINT8)); for compound datasets none of the
    // primitive H5Tequal checks match, so precision falls through
    // (here we treat anything other than UInt8 as compound).
    if ([ds precision] == TTIOPrecisionUInt8) {
        uint8_t codec_id = 0;
        id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
        if ([sig respondsToSelector:@selector(unwrap)]) {
            TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
            TTIOHDF5Dataset *hds = [hg openDatasetNamed:@"read_names"
                                                  error:NULL];
            if (hds) {
                codec_id = _ttio_m86_read_compression_attr([hds datasetId]);
            }
        } else {
            codec_id = _ttio_m86_read_compression_attr_protocol(ds);
        }
        // NAME_TOKENIZED v1 (codec id 8) reader
        // path removed — reject with a clear error so legacy files
        // surface a re-encode hint instead of silently mis-decoding.
        if (codec_id == (uint8_t)8 /* NAME_TOKENIZED v1 */) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2041
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"NAME_TOKENIZED v1 (codec id 8) is no "
                           @"longer supported in v1.0; file was "
                           @"written with an older TTI-O version. "
                           @"Re-encode with v1.0+ which uses "
                           @"NAME_TOKENIZED_V2 (codec id 15)."}];
            return nil;
        }
        if (codec_id != (uint8_t)15 /* NAME_TOKENIZED_V2 */) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2041
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"signal_channel 'read_names': "
                                @"@compression=%u is not a supported "
                                @"TTIO codec id for the read_names "
                                @"channel (only NAME_TOKENIZED_V2 = "
                                @"15 is recognised under v1.0)",
                                (unsigned)codec_id]}];
            return nil;
        }
        id allRaw = [ds readAll:error];
        if (![allRaw isKindOfClass:[NSData class]]) return nil;
        NSData *encoded = (NSData *)allRaw;
        // Empty-run short-circuit: zero-length blob → empty list.
        if (encoded.length == 0) {
            _decodedReadNames = @[];
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2040
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"read_names index %lu out of range "
                                @"[0, 0) — empty read_names blob",
                                (unsigned long)i]}];
            return nil;
        }
        // Route NAME_TOKENIZED_V2 decode through the codec registry
        // (str-list adapter). Byte-identical to the old inline
        // [TTIONameTokenizerV2 decodeData:] call.
        NSError *decErr = nil;
        id<TTIOCodec> codec =
            [TTIOCodecRegistry codecForId:TTIOCompressionNameTokenizedV2];
        TTIODecodedChannel *dc =
            [codec decode:[[TTIOBytesPayload alloc] initWithBytes:encoded]
                  context:[self _codecContext] error:&decErr];
        NSArray<NSString *> *decoded =
            (dc != nil) ? ((TTIODecodedStringList *)dc).names : nil;
        if (decoded == nil) {
            if (error) *error = decErr ?: [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2042
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"signal_channel 'read_names' "
                           @"NAME_TOKENIZED_V2 decode failed"}];
            return nil;
        }
        _decodedReadNames = [decoded copy];
        if (i >= _decodedReadNames.count) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2043
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"read_names index %lu out of range "
                                @"[0, %lu) after NAME_TOKENIZED decode",
                                (unsigned long)i,
                                (unsigned long)_decodedReadNames.count]}];
            return nil;
        }
        return _decodedReadNames[i];
    }

    // M82 read_names compound layout removed.
    // A non-UInt8 read_names dataset is from an older TTI-O version.
    if (error) *error = [NSError
        errorWithDomain:@"TTIOGenomicRun" code:2044
               userInfo:@{NSLocalizedDescriptionKey:
                   @"signal_channels/read_names is a compound layout "
                   @"(M82 / VL-string). The compound layout is no "
                   @"longer supported in v1.0; file was written with "
                   @"an older TTI-O version. Re-encode with v1.0+ "
                   @"which uses NAME_TOKENIZED_V2 (codec id 15)."}];
    return nil;
}

// unsigned LEB128 varint reader for the cigars rANS
// length-prefix-concat path. Mirrors NAME_TOKENIZED's varint_read
// (in TTIONameTokenizer.m) — reproduced here to avoid coupling the
// run reader to the codec module's private symbols. Returns 1 on
// success and advances *io_offset past the consumed bytes; returns
// 0 on truncated/oversize varints (>10 bytes / >64 bits).
static int _ttio_m86_cigars_varint_read(const uint8_t *buf, size_t buf_len,
                                         size_t *io_offset,
                                         uint64_t *out_value)
{
    uint64_t value = 0;
    int shift = 0;
    size_t pos = *io_offset;
    for (;;) {
        if (pos >= buf_len) return 0;
        const uint8_t b = buf[pos++];
        if (shift >= 64) return 0;
        value |= ((uint64_t)(b & 0x7Fu)) << shift;
        if ((b & 0x80u) == 0) {
            *io_offset = pos;
            *out_value = value;
            return 1;
        }
        shift += 7;
    }
}

// cigars dispatch helper.
//
// The cigars channel has two on-disk layouts (Binding Decisions
// §120-§123, HANDOFF M86 Phase C §2.7):
//
//   - **M82 compound** (no override): VL_STRING-in-compound dataset,
//     read whole-and-cache via -compoundRowsNamed:.
//   - **TTIO codec** (override active): flat 1-D uint8 dataset
//     of the same name carrying the codec output, with @compression
//     in {4, 5, 8}. Decoded once on first access and cached as
//     NSArray<NSString *> on this TTIOGenomicRun instance per
//     Binding Decision §123 — a separate field from
//     _decodedReadNames since the two channels have independent
//     dispatch shapes.
//
//     * @compression == 4 (RANS_ORDER0) or 5 (RANS_ORDER1): the
//       decoded byte buffer is a length-prefix-concat sequence
//       (varint(len) + bytes per CIGAR; §2.5 of the Phase C plan;
//       Gotcha §139). Walk the buffer until exhausted.
//     * @compression == 8 (NAME_TOKENIZED): pass the bytes through
//       TTIONameTokenizerDecode directly (the codec's self-describing
//       wire format records the read count internally).
//
// Dispatch is on dataset shape — a flat uint8 dataset routes through
// the codec path; otherwise fall through to the compound path (same
// pattern Phase E uses for read_names).
- (NSString *)cigarAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) {
            if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:2060
                userInfo:@{NSLocalizedDescriptionKey: [NSString stringWithFormat:
                    @"cigars index %lu out of range [0, %llu)", (unsigned long)i, _blockTable.readCount]}];
            return nil;
        }
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view cigarAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : nil;
    }
    if (_decodedCigars != nil) {
        if (i >= _decodedCigars.count) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2060
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"cigars index %lu out of range "
                                @"[0, %lu)",
                                (unsigned long)i,
                                (unsigned long)_decodedCigars.count]}];
            return nil;
        }
        return _decodedCigars[i];
    }

    id<TTIOStorageDataset> ds = [self signalDatasetNamed:@"cigars"
                                                   error:error];
    if (!ds) return nil;

    // Shape dispatch: precision == UInt8 is a flat uint8 dataset and
    // therefore the codec path; anything else (compound) falls through
    // to the M82 path. Mirrors -readNameAtIndex:'s shape check.
    if ([ds precision] == TTIOPrecisionUInt8) {
        uint8_t codec_id = 0;
        id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
        if ([sig respondsToSelector:@selector(unwrap)]) {
            TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
            TTIOHDF5Dataset *hds = [hg openDatasetNamed:@"cigars"
                                                  error:NULL];
            if (hds) {
                codec_id = _ttio_m86_read_compression_attr([hds datasetId]);
            }
        } else {
            codec_id = _ttio_m86_read_compression_attr_protocol(ds);
        }

        id allRaw = [ds readAll:error];
        if (![allRaw isKindOfClass:[NSData class]]) return nil;
        NSData *encoded = (NSData *)allRaw;

        if (codec_id == (uint8_t)4 /* RANS_ORDER0 */
            || codec_id == (uint8_t)5 /* RANS_ORDER1 */) {
            // Route ONLY the inner rANS decode through the registry; the
            // length-prefix-concat framing below stays here (the codec
            // layer is byte-stream agnostic). Byte-identical to the old
            // TTIORansDecode(encoded, ...) call.
            NSError *decErr = nil;
            id<TTIOCodec> ransCodec =
                [TTIOCodecRegistry codecForId:(TTIOCompression)codec_id];
            TTIODecodedChannel *rdc =
                [ransCodec decode:[[TTIOBytesPayload alloc] initWithBytes:encoded]
                          context:[self _codecContext] error:&decErr];
            NSData *decoded =
                (rdc != nil) ? ((TTIODecodedBytes *)rdc).data : nil;
            if (decoded == nil) {
                if (error) *error = decErr ?: [NSError
                    errorWithDomain:@"TTIOGenomicRun" code:2061
                           userInfo:@{NSLocalizedDescriptionKey:
                               @"signal_channel 'cigars' rANS decode "
                               @"failed"}];
                return nil;
            }
            // Walk length-prefix-concat: varint(len) + len bytes per
            // CIGAR, repeated until the decoded buffer is exhausted.
            const uint8_t *buf = (const uint8_t *)decoded.bytes;
            const size_t   n   = decoded.length;
            size_t off = 0;
            NSMutableArray<NSString *> *out = [NSMutableArray array];
            while (off < n) {
                uint64_t len = 0;
                if (!_ttio_m86_cigars_varint_read(buf, n, &off, &len)) {
                    if (error) *error = [NSError
                        errorWithDomain:@"TTIOGenomicRun" code:2062
                               userInfo:@{NSLocalizedDescriptionKey:
                                   @"signal_channel 'cigars' rANS "
                                   @"length-prefix-concat: truncated "
                                   @"varint length prefix"}];
                    return nil;
                }
                if (off + (size_t)len > n) {
                    if (error) *error = [NSError
                        errorWithDomain:@"TTIOGenomicRun" code:2063
                               userInfo:@{NSLocalizedDescriptionKey:
                                   [NSString stringWithFormat:
                                        @"signal_channel 'cigars' rANS "
                                        @"length-prefix-concat: entry "
                                        @"runs off end of decoded buffer "
                                        @"(offset=%zu, length=%llu, "
                                        @"buffer_size=%zu)",
                                        off,
                                        (unsigned long long)len, n]}];
                    return nil;
                }
                NSString *cig = [[NSString alloc]
                    initWithBytes:buf + off
                           length:(NSUInteger)len
                         encoding:NSASCIIStringEncoding];
                if (cig == nil) {
                    if (error) *error = [NSError
                        errorWithDomain:@"TTIOGenomicRun" code:2064
                               userInfo:@{NSLocalizedDescriptionKey:
                                   @"signal_channel 'cigars' rANS "
                                   @"length-prefix-concat: entry "
                                   @"contains non-ASCII bytes"}];
                    return nil;
                }
                [out addObject:cig];
                off += (size_t)len;
            }
            _decodedCigars = [out copy];
        } else if (codec_id == (uint8_t)8 /* NAME_TOKENIZED v1 */) {
            // NAME_TOKENIZED v1 reader removed.
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2065
                       userInfo:@{NSLocalizedDescriptionKey:
                           @"NAME_TOKENIZED v1 (codec id 8) is no "
                           @"longer supported in v1.0; cigars dataset "
                           @"was written with an older TTI-O version. "
                           @"Re-encode with v1.0+ which uses RANS "
                           @"on the cigars channel."}];
            return nil;
        } else {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2066
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"signal_channel 'cigars': "
                                @"@compression=%u is not a supported "
                                @"TTIO codec id for the cigars channel "
                                @"(only RANS_ORDER0 = 4 and "
                                @"RANS_ORDER1 = 5 are recognised "
                                @"under v1.0)",
                                (unsigned)codec_id]}];
            return nil;
        }

        if (i >= _decodedCigars.count) {
            if (error) *error = [NSError
                errorWithDomain:@"TTIOGenomicRun" code:2067
                       userInfo:@{NSLocalizedDescriptionKey:
                           [NSString stringWithFormat:
                                @"cigars index %lu out of range "
                                @"[0, %lu) after codec decode",
                                (unsigned long)i,
                                (unsigned long)_decodedCigars.count]}];
            return nil;
        }
        return _decodedCigars[i];
    }

    // Compound path (M82, no override). Materialise the whole
    // cigar list on first call and cache in _decodedCigars —
    // without this, per-record cigarAtIndex: allocates a fresh
    // NSString from NSData on every call, dominating the per-
    // record time on the genomic transport encode path (mirrors
    // Java fix / Python parity).
    TTIOCompoundField *vlValue =
        [TTIOCompoundField fieldWithName:@"value"
                                    kind:TTIOCompoundFieldKindVLString];
    NSArray *cigars = [self compoundRowsNamed:@"cigars"
                                         field:vlValue
                                         error:error];
    if (!cigars) return nil;
    if (i >= cigars.count) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOGenomicRun" code:2068
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"cigars index %lu out of range [0, %lu)",
                            (unsigned long)i,
                            (unsigned long)cigars.count]}];
        return nil;
    }
    NSMutableArray<NSString *> *out =
        [NSMutableArray arrayWithCapacity:cigars.count];
    for (NSDictionary *row in cigars) {
        id v = row[@"value"];
        if ([v isKindOfClass:[NSData class]]) {
            [out addObject:[[NSString alloc] initWithData:v
                                                  encoding:NSUTF8StringEncoding]
                          ?: @""];
        } else if ([v isKindOfClass:[NSString class]]) {
            [out addObject:(NSString *)v];
        } else {
            [out addObject:@""];
        }
    }
    _decodedCigars = [out copy];
    return _decodedCigars[i];
}

// integer-channel array reader.
//
// Returns the full integer signal-channel array (positions, flags,
// or mapping_qualities) as an NSData carrying the LE byte
// representation of the dtype implied by channel-name lookup
// (Binding Decision §115). Two paths:
//
//   - **Uncompressed (no @compression or @compression == 0):** read
//     the typed dataset directly via the storage protocol; bytes are
//     already in LE order (HDF5 stores native little-endian on
//     x86/ARM). Cache and return.
//
//   - **rANS (@compression == 4 or 5):** read the dataset whole as
//     uint8 bytes, decode through TTIORansDecode, cache and return.
//
// Per Binding Decision §119 the per-read access path
// (-readAtIndex:) does NOT consume this helper; it continues to use
// self.index.{positions,mappingQualities,flags}. This helper is
// wired for round-trip conformance and for any future reader that
// prefers signal_channels/ over genomic_index/ (Phase B is primarily
// a write-side file-size optimisation).
// v1.6 (L4): -intChannelArrayNamed:error: removed. The helper read
// positions/flags/mapping_qualities from signal_channels/ via codec
// dispatch — but those datasets no longer exist in v1.6 files. See
// docs/format-spec.md §10.7.

// HDF5 link-type query for signal_channels/mate_info.
// Per Binding Decision §128 / Gotcha §141, dispatch is on HDF5 link
// type (dataset = M82 compound; group = Phase F subgroup), NOT on
// @compression attribute presence on the bare link (the attribute
// lives on per-field child datasets within the subgroup, not on the
// bare link). Probed once on first access and cached on the run.
//
// HDF5 backend: H5Oget_info_by_name returns a struct whose `type`
// field is H5O_TYPE_GROUP or H5O_TYPE_DATASET — the cleanest signal
// for the dispatch. For non-HDF5 backends we fall back to the
// storage protocol's openGroupNamed/openDatasetNamed combination
// (one returns nil where the other doesn't).
- (BOOL)_mateInfoIsSubgroup
{
    if (_mateInfoLinkType >= 0) {
        // linkType 0 = compound, 1 = Phase-F subgroup, 2 = inline_v2 subgroup.
        return _mateInfoLinkType >= 1;
    }
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) {
        // Defensive: if signal_channels is unavailable, assume M82.
        _mateInfoLinkType = 0;
        return NO;
    }
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        if (hg) {
            H5O_info2_t info;
            herr_t s = H5Oget_info_by_name3([hg groupId], "mate_info",
                                           &info, H5O_INFO_BASIC, H5P_DEFAULT);
            if (s >= 0 && info.type == H5O_TYPE_GROUP) {
                // It is a group. Probe further for inline_v2 dataset.
                TTIOHDF5Group *mateGrp = [hg openGroupNamed:@"mate_info" error:NULL];
                if (mateGrp) {
                    H5O_info2_t dsInfo;
                    herr_t s2 = H5Oget_info_by_name3([mateGrp groupId],
                                                    "inline_v2", &dsInfo, H5O_INFO_BASIC, H5P_DEFAULT);
                    if (s2 >= 0 && dsInfo.type == H5O_TYPE_DATASET) {
                        _mateInfoLinkType = 2;  // v1.7 inline_v2
                        return YES;
                    }
                }
                _mateInfoLinkType = 1;  // Phase-F per-field subgroup
                return YES;
            }
            if (s >= 0) {
                _mateInfoLinkType = 0;  // dataset = M82 compound
                return NO;
            }
        }
        // H5Oget_info_by_name failed — assume compound (legacy default).
        _mateInfoLinkType = 0;
        return NO;
    }
    // Storage-protocol path: try openGroupNamed first.
    NSError *gErr = nil;
    id<TTIOStorageGroup> sub = [sig openGroupNamed:@"mate_info" error:&gErr];
    if (sub != nil) {
        // Probe for inline_v2 child dataset.
        NSError *dsErr = nil;
        id<TTIOStorageDataset> inlineDs =
            [sub openDatasetNamed:@"inline_v2" error:&dsErr];
        if (inlineDs != nil) {
            _mateInfoLinkType = 2;
        } else {
            _mateInfoLinkType = 1;
        }
        return YES;
    }
    _mateInfoLinkType = 0;
    return NO;
}

/** YES when signal_channels/mate_info/inline_v2 exists. */
- (BOOL)_mateInfoIsInlineV2
{
    // Force probe if not yet done.
    (void)[self _mateInfoIsSubgroup];
    return _mateInfoLinkType == 2;
}

// ── sequences GROUP probe + refdiff_v2 decoder ─────────────────────

/** probe whether signal_channels/sequences is a GROUP (refdiff_v2
 *  layout) or a flat dataset (all v1 layouts). Cached on first call. */
- (BOOL)_sequencesIsRefDiffV2
{
    if (_sequencesLinkType >= 0) return _sequencesLinkType == 1;

    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) {
        _sequencesLinkType = 0;
        return NO;
    }
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        H5O_info2_t info;
        memset(&info, 0, sizeof(info));
        herr_t s = H5Oget_info_by_name3([hg groupId], "sequences",
                                       &info, H5O_INFO_BASIC, H5P_DEFAULT);
        if (s >= 0 && info.type == H5O_TYPE_GROUP) {
            _sequencesLinkType = 1;
            return YES;
        }
        _sequencesLinkType = 0;
        return NO;
    }
    // Storage-protocol path: try openGroupNamed first.
    NSError *gErr = nil;
    id<TTIOStorageGroup> sub = [sig openGroupNamed:@"sequences" error:&gErr];
    if (sub != nil) {
        _sequencesLinkType = 1;
        return YES;
    }
    _sequencesLinkType = 0;
    return NO;
}

// Task 5: -_decodeRefDiffV2Sequences: was deleted — its body now lives
// in the REF_DIFF_V2 codec-registry adapter (_TTIORefDiffCodec in
// TTIOCodecRegistry.m). -byteChannelSliceNamed: routes the sequences
// refdiff_v2 group payload through that adapter with the context built
// by -_codecContext (positions / totalBases / chromosomes /
// cigarsProvider / referenceResolver), staying byte-identical to the
// old side path.

/** decode the inline_v2 blob; populate _decodedMateInfo with
 *  "chrom" (NSArray<NSString *>), "pos" (NSData int64), "tlen" (NSData int32).
 *  Returns NO + error on failure. Caches on success. */
- (BOOL)_decodeMateInfoInlineV2:(NSError **)error
{
    // Already cached?
    if (_decodedMateInfo[@"chrom"]) return YES;

    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:error];
    if (!sig) return NO;

    // Open the mate_info group.
    TTIOHDF5Group *mateH5 = nil;
    id<TTIOStorageGroup> mateProt = nil;
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        mateH5 = [hg openGroupNamed:@"mate_info" error:error];
        if (!mateH5) return NO;
    } else {
        mateProt = [sig openGroupNamed:@"mate_info" error:error];
        if (!mateProt) return NO;
    }

    // Read the inline_v2 blob.
    NSData *blob = nil;
    if (mateH5) {
        TTIOHDF5Dataset *ds = [mateH5 openDatasetNamed:@"inline_v2" error:error];
        if (!ds) return NO;
        id raw = [ds readDataWithError:error];
        if (![raw isKindOfClass:[NSData class]]) return NO;
        blob = (NSData *)raw;
    } else {
        id<TTIOStorageDataset> ds = [mateProt openDatasetNamed:@"inline_v2" error:error];
        if (!ds) return NO;
        id raw = [ds readAll:error];
        if (![raw isKindOfClass:[NSData class]]) return NO;
        blob = (NSData *)raw;
    }

    NSUInteger n = [self index].count;

    // own_chrom_ids (encounter-order uint16, must match writer) and
    // own_positions are now derived in -_codecContext and consumed by
    // the registry adapter below — see that method for the exact
    // encounter-order derivation this path used to inline here.

    // Decode via the codec registry (MATE_INLINE_V2 adapter). The
    // adapter consumes context.ownChromIds / ownPositions / nRecords —
    // those context fields are derived in -_codecContext using the same
    // encounter-order uint16 derivation as the ownChromIds built above,
    // so this is byte-identical to the old inline TTIOMateInfoV2 call.
    NSData *outMc = nil, *outMp = nil, *outTs = nil;
    NSError *decErr = nil;
    id<TTIOCodec> codec =
        [TTIOCodecRegistry codecForId:TTIOCompressionMateInlineV2];
    TTIODecodedChannel *dc =
        [codec decode:[[TTIOBytesPayload alloc] initWithBytes:blob]
              context:[self _codecContext] error:&decErr];
    if (dc == nil) {
        if (error) *error = decErr ?: [NSError
            errorWithDomain:@"TTIOGenomicRun" code:2090
                   userInfo:@{NSLocalizedDescriptionKey:
                       @"v1.7 inline_v2 decode failed"}];
        return NO;
    }
    TTIODecodedMateInfo *mi = (TTIODecodedMateInfo *)dc;
    outMc = mi.mateChromIds;
    outMp = mi.matePositions;
    outTs = mi.templateLengths;

    // Read the chrom_names sidecar compound to resolve mate chrom ids → names.
    NSArray *chromNameRows = nil;
    TTIOCompoundField *nameField =
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString];
    if (mateH5) {
        chromNameRows = [TTIOCompoundIO readGenericFromGroup:mateH5
                                                datasetNamed:@"chrom_names"
                                                      fields:@[nameField]
                                                       error:error];
    } else {
        id<TTIOStorageDataset> namesDs =
            [mateProt openDatasetNamed:@"chrom_names" error:error];
        if (namesDs) chromNameRows = [namesDs readAll:error];
    }
    if (!chromNameRows) return NO;

    // Build chrom_id → name table (row index = chrom_id).
    NSMutableArray<NSString *> *chromNamesById =
        [NSMutableArray arrayWithCapacity:chromNameRows.count];
    for (NSDictionary *row in chromNameRows) {
        id v = row[@"name"];
        NSString *s = [v isKindOfClass:[NSData class]]
            ? [[NSString alloc] initWithData:v encoding:NSUTF8StringEncoding]
            : (NSString *)v;
        [chromNamesById addObject:s ?: @""];
    }

    // Convert mate_chrom_ids (int32, -1=unmapped) back to chromosome name strings.
    const int32_t *mcPtr = (const int32_t *)outMc.bytes;
    NSMutableArray<NSString *> *mateChroms =
        [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        int32_t iv = mcPtr[i];
        if (iv == -1) {
            [mateChroms addObject:@"*"];
        } else if (iv >= 0 && (NSUInteger)iv < chromNamesById.count) {
            [mateChroms addObject:chromNamesById[iv]];
        } else {
            [mateChroms addObject:
                [NSString stringWithFormat:@"chr_id_%d", iv]];
        }
    }

    _decodedMateInfo[@"chrom"] = [mateChroms copy];
    _decodedMateInfo[@"pos"]   = outMp;
    _decodedMateInfo[@"tlen"]  = outTs;
    return YES;
}

// per-read mate-field accessors recognise only
// the inline_v2 layout (v1.7 codec id 13). The Phase F per-field
// subgroup (linkType 1) and M82 compound (linkType 0) layouts are
// rejected with a clear NSError directing callers at the v2 codec.
static void _ttio_v17_reject_legacy_mate_layout(NSError **error)
{
    if (error) *error = [NSError
        errorWithDomain:@"TTIOGenomicRun" code:2080
               userInfo:@{NSLocalizedDescriptionKey:
                   @"signal_channels/mate_info layout is not the "
                   @"inline_v2 codec (id 13). The Phase F per-field "
                   @"subgroup and M82 compound layouts are no longer "
                   @"supported in v1.0; file was written with an "
                   @"older TTI-O version. Re-encode with v1.0+."}];
}

- (NSString *)_mateChromAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) return nil;
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view _mateChromAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : nil;
    }
    if ([self _mateInfoIsInlineV2]) {
        if (![self _decodeMateInfoInlineV2:error]) return nil;
        NSArray<NSString *> *chroms = _decodedMateInfo[@"chrom"];
        if (!chroms || i >= chroms.count) return nil;
        return chroms[i];
    }
    _ttio_v17_reject_legacy_mate_layout(error);
    return nil;
}

- (int64_t)_matePosAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) return 0;
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view _matePosAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : 0;
    }
    if ([self _mateInfoIsInlineV2]) {
        if (![self _decodeMateInfoInlineV2:error]) return 0;
        NSData *bytes = _decodedMateInfo[@"pos"];
        if (!bytes) return 0;
        NSUInteger n = bytes.length / sizeof(int64_t);
        if (i >= n) return 0;
        int64_t v; memcpy(&v, (const int64_t *)bytes.bytes + i, sizeof(int64_t));
        return v;
    }
    _ttio_v17_reject_legacy_mate_layout(error);
    return 0;
}

- (int32_t)_mateTlenAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) return 0;
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view _mateTlenAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : 0;
    }
    if ([self _mateInfoIsInlineV2]) {
        if (![self _decodeMateInfoInlineV2:error]) return 0;
        NSData *bytes = _decodedMateInfo[@"tlen"];
        if (!bytes) return 0;
        NSUInteger n = bytes.length / sizeof(int32_t);
        if (i >= n) return 0;
        int32_t v; memcpy(&v, (const int32_t *)bytes.bytes + i, sizeof(int32_t));
        return v;
    }
    _ttio_v17_reject_legacy_mate_layout(error);
    return 0;
}

- (TTIOAlignedRead *)readAtIndex:(NSUInteger)i error:(NSError **)error
{
    if (_blockTable) {
        i = [self _storedRow:i];
        NSUInteger b = [_blockTable blockForRead:i];
        if (b == NSNotFound) {
            if (error) *error = [NSError errorWithDomain:@"TTIOGenomicRun" code:0
                userInfo:@{NSLocalizedDescriptionKey: [NSString stringWithFormat:
                    @"index %lu out of range [0, %llu)", (unsigned long)i, _blockTable.readCount]}];
            return nil;
        }
        TTIOGenomicRun *view = [self _blockView:b error:error];
        return view ? [view readAtIndex:i - (NSUInteger)[_blockTable readStartAt:b] error:error] : nil;
    }
    if (i >= [self index].count) {
        if (error) *error = [NSError
            errorWithDomain:@"TTIOGenomicRun" code:0
                   userInfo:@{NSLocalizedDescriptionKey:
                       [NSString stringWithFormat:
                            @"index %lu out of range [0, %lu)",
                            (unsigned long)i, (unsigned long)[self index].count]}];
        return nil;
    }

    uint64_t offset = [[self index] offsetAt:i];
    uint32_t length = [[self index] lengthAt:i];

    int64_t  position = [[self index] positionAt:i];
    uint8_t  mapq     = [[self index] mappingQualityAt:i];
    uint32_t flag     = [[self index] flagsAt:i];
    NSString *chrom   = [[self index] chromosomeAt:i];

    // routed through byteChannelSliceNamed: so codec-compressed
    // channels (@compression > 0) are decoded transparently before
    // slicing. Uncompressed channels go through the existing
    // hyperslab path.
    NSData *seqData = [self byteChannelSliceNamed:@"sequences"
                                            offset:offset count:length
                                             error:error];
    if (!seqData) return nil;
    NSString *sequence = [[NSString alloc] initWithData:seqData
                                               encoding:NSASCIIStringEncoding];

    NSData *qualities = [self byteChannelSliceNamed:@"qualities"
                                              offset:offset count:length
                                               error:error];
    if (!qualities) return nil;

    // route cigars through the shape-dispatching helper
    // so the schema-lifted (flat uint8 + RANS or NAME_TOKENIZED)
    // layout is decoded transparently. The compound (M82) layout
    // continues to use the existing -compoundRowsNamed: cache via
    // the helper's compound fall-through.
    NSString *cigar = [self cigarAtIndex:i error:error];
    if (!cigar && error && *error) return nil;

    // route read_names through the shape-dispatching
    // helper so the schema-lifted (flat uint8 + NAME_TOKENIZED)
    // layout is decoded transparently. The compound (M82) layout
    // continues to use the existing -compoundRowsNamed: cache.
    NSString *readName = [self readNameAtIndex:i error:error];
    if (!readName && error && *error) return nil;

    // route mate-field reads through the per-field
    // dispatch helpers. The link-type query (group vs dataset) for
    // signal_channels/mate_info is cached on the run, so the three
    // accessors are essentially free after the first call. The
    // existing M82 compound path is preserved inside the helpers
    // (see -_mateChromAtIndex: et al.).
    NSError *mErr = nil;
    NSString *mateChromosome = [self _mateChromAtIndex:i error:&mErr];
    if (!mateChromosome && mErr) {
        if (error) *error = mErr;
        return nil;
    }
    int64_t matePosition = [self _matePosAtIndex:i error:&mErr];
    int32_t templateLength = [self _mateTlenAtIndex:i error:&mErr];

    // M101: the SAM optional fields, decoded once per run / block view.
    NSArray<NSString *> *tags = [self allTagsWithError:error];
    if (!tags) return nil;

    return [[TTIOAlignedRead alloc]
        initWithReadName:readName
              chromosome:chrom
                position:position
          mappingQuality:mapq
                   cigar:cigar
                sequence:sequence
               qualities:qualities
                   flags:flag
          mateChromosome:mateChromosome
            matePosition:matePosition
          templateLength:templateLength
                    tags:(i < tags.count ? tags[i] : @"")];
}

- (NSArray<TTIOAlignedRead *> *)readsInRegion:(NSString *)chromosome
                                          start:(int64_t)start
                                            end:(int64_t)end
{
    NSIndexSet *indices = [[self index] indicesForRegion:chromosome start:start end:end];
    NSMutableArray *result = [NSMutableArray arrayWithCapacity:indices.count];
    [indices enumerateIndexesUsingBlock:^(NSUInteger idx, BOOL *stop) {
        NSError *err = nil;
        TTIOAlignedRead *r = [self readAtIndex:idx error:&err];
        if (r) [result addObject:r];
    }];
    return result;
}

#pragma mark - TTIOIndexable / TTIORun (Phase 1)

- (NSUInteger)count
{
    return [self readCount];
}

- (id)objectAtIndex:(NSUInteger)index
{
    return [self readAtIndex:index error:NULL];
}

- (NSArray<TTIOProvenanceRecord *> *)provenanceChain
{
    // Mirrors Python GenomicRun.provenance_chain() — closes the M91
    // read-side gap. Reads from <run>/provenance/steps using the same
    // compound layout as TTIOAcquisitionRun. Returns @[] for runs with
    // no provenance attached.
    if ([_group respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *runH5 = [(id)_group performSelector:@selector(unwrap)];
        if (runH5 && [runH5 hasChildNamed:@"provenance"]) {
            TTIOHDF5Group *provGroup =
                [runH5 openGroupNamed:@"provenance" error:NULL];
            if (provGroup && [provGroup hasChildNamed:@"steps"]) {
                NSArray *records =
                    [TTIOCompoundIO readProvenanceFromGroup:provGroup
                                               datasetNamed:@"steps"
                                                      error:NULL];
                if (records) return [records copy];
            }
        }
    }
    // The storage-protocol writers keep the JSON mirror only.
    if ([_group hasAttributeNamed:@"provenance_json"]) {
        id v = [_group attributeValueForName:@"provenance_json" error:NULL];
        NSString *json = [v isKindOfClass:[NSString class]] ? v : [v description];
        NSData *data = [json dataUsingEncoding:NSUTF8StringEncoding];
        NSArray *plists = data ? [NSJSONSerialization JSONObjectWithData:data options:0 error:NULL] : nil;
        if ([plists isKindOfClass:[NSArray class]]) {
            NSMutableArray *out = [NSMutableArray arrayWithCapacity:plists.count];
            for (NSDictionary *pl in plists) {
                if ([pl isKindOfClass:[NSDictionary class]]) {
                    [out addObject:[TTIOProvenanceRecord fromPlist:pl]];
                }
            }
            return out;
        }
    }
    return @[];
}

#pragma mark - Phase 2c-T verbatim v2 blob accessors

- (nullable NSData *)readMateInfoInlineV2BlobBytes
{
    if (_blockTable) {
        // A grouped run's blobs are in stored order (M103).
        if (_storedOf != nil) return nil;
        return _blockTable.count == 1 ? [[self _blockView:0 error:NULL] readMateInfoInlineV2BlobBytes] : nil;
    }
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) return nil;
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        if (![hg hasChildNamed:@"mate_info"]) return nil;
        TTIOHDF5Group *mateGrp = [hg openGroupNamed:@"mate_info" error:NULL];
        if (!mateGrp || ![mateGrp hasChildNamed:@"inline_v2"]) return nil;
        TTIOHDF5Dataset *ds = [mateGrp openDatasetNamed:@"inline_v2" error:NULL];
        if (!ds) return nil;
        id raw = [ds readDataWithError:NULL];
        return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : nil;
    }
    NSError *gErr = nil;
    id<TTIOStorageGroup> mateProt = [sig openGroupNamed:@"mate_info" error:&gErr];
    if (!mateProt) return nil;
    id<TTIOStorageDataset> ds = [mateProt openDatasetNamed:@"inline_v2" error:NULL];
    if (!ds) return nil;
    id raw = [ds readAll:NULL];
    return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : nil;
}

- (NSArray<NSString *> *)readMateInfoChromNamesTable
{
    NSMutableArray<NSString *> *out = [NSMutableArray array];
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) return out;
    NSArray *rows = nil;
    TTIOCompoundField *nameField =
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString];
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        if (![hg hasChildNamed:@"mate_info"]) return out;
        TTIOHDF5Group *mateGrp = [hg openGroupNamed:@"mate_info" error:NULL];
        if (!mateGrp || ![mateGrp hasChildNamed:@"chrom_names"]) return out;
        rows = [TTIOCompoundIO readGenericFromGroup:mateGrp
                                       datasetNamed:@"chrom_names"
                                             fields:@[nameField]
                                              error:NULL];
    } else {
        id<TTIOStorageGroup> mateProt = [sig openGroupNamed:@"mate_info" error:NULL];
        if (!mateProt) return out;
        id<TTIOStorageDataset> ds = [mateProt openDatasetNamed:@"chrom_names" error:NULL];
        if (ds) rows = [ds readAll:NULL];
    }
    for (NSDictionary *row in rows) {
        id v = row[@"name"];
        NSString *s = [v isKindOfClass:[NSData class]]
            ? [[NSString alloc] initWithData:v encoding:NSUTF8StringEncoding]
            : (NSString *)v;
        [out addObject:s ?: @""];
    }
    return out;
}

- (nullable NSData *)readNameTokV2BlobBytes
{
    if (_blockTable) {
        // A grouped run's blobs are in stored order (M103).
        if (_storedOf != nil) return nil;
        return _blockTable.count == 1 ? [[self _blockView:0 error:NULL] readNameTokV2BlobBytes] : nil;
    }
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) return nil;
    uint8_t codec = [self wireCompressionForChannel:@"read_names"];
    // 15 = TTIOCompressionNameTokenizedV2 (avoid Transport→Genomics
    // import inversion).
    if (codec != 15) return nil;
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        if (![hg hasChildNamed:@"read_names"]) return nil;
        TTIOHDF5Dataset *ds = [hg openDatasetNamed:@"read_names" error:NULL];
        if (!ds) return nil;
        id raw = [ds readDataWithError:NULL];
        return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : [NSData data];
    }
    id<TTIOStorageDataset> ds = [sig openDatasetNamed:@"read_names" error:NULL];
    if (!ds) return nil;
    id raw = [ds readAll:NULL];
    return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : [NSData data];
}

- (nullable NSData *)readRefDiffV2BlobBytes
{
    if (_blockTable) {
        // A grouped run's blobs are in stored order (M103).
        if (_storedOf != nil) return nil;
        return _blockTable.count == 1 ? [[self _blockView:0 error:NULL] readRefDiffV2BlobBytes] : nil;
    }
    id<TTIOStorageGroup> sig = [self signalChannelsGroupWithError:NULL];
    if (!sig) return nil;
    if ([sig respondsToSelector:@selector(unwrap)]) {
        TTIOHDF5Group *hg = [(id)sig performSelector:@selector(unwrap)];
        if (![hg hasChildNamed:@"sequences"]) return nil;
        H5O_info2_t info; memset(&info, 0, sizeof(info));
        if (H5Oget_info_by_name3([hg groupId], "sequences",
                                &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) return nil;
        if (info.type != H5O_TYPE_GROUP) return nil;
        TTIOHDF5Group *seqGrp = [hg openGroupNamed:@"sequences" error:NULL];
        if (!seqGrp || ![seqGrp hasChildNamed:@"refdiff_v2"]) return nil;
        TTIOHDF5Dataset *ds = [seqGrp openDatasetNamed:@"refdiff_v2" error:NULL];
        if (!ds) return nil;
        id raw = [ds readDataWithError:NULL];
        return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : nil;
    }
    NSError *gErr = nil;
    id<TTIOStorageGroup> seqGrp = [sig openGroupNamed:@"sequences" error:&gErr];
    if (!seqGrp) return nil;
    id<TTIOStorageDataset> ds = [seqGrp openDatasetNamed:@"refdiff_v2" error:NULL];
    if (!ds) return nil;
    id raw = [ds readAll:NULL];
    return [raw isKindOfClass:[NSData class]] ? (NSData *)raw : nil;
}

#pragma mark - Bulk accessors for hot serialization paths

/* A stored-order per-read list of a grouped run in input order. */
- (NSArray *)_inInputOrder:(NSArray *)stored
{
    NSUInteger n = _storedOf.length / sizeof(uint32_t);
    if (stored.count != n) return stored;
    const uint32_t *so = (const uint32_t *)_storedOf.bytes;
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) [out addObject:stored[so[i]]];
    return out;
}

- (NSUInteger)_totalBaseCount
{
    if (_blockTable) {
        NSUInteger nb = _blockTable.count;
        return nb == 0 ? 0 : (NSUInteger)([_blockTable baseStartAt:nb - 1] + [_blockTable nBasesAt:nb - 1]);
    }
    NSUInteger n = [self index].count;
    if (n == 0) return 0;
    return (NSUInteger)([[self index] offsetAt:n - 1] + [[self index] lengthAt:n - 1]);
}

- (NSData *)wholeSequencesData
{
    return [self _byteChannelFullNamed:@"sequences"];
}

- (NSData *)wholeQualitiesData
{
    return [self _byteChannelFullNamed:@"qualities"];
}

// Cache-priming whole-channel fetch. byteChannelSliceNamed: caches
// codec-compressed channels but returns the raw HDF5 buffer for
// uncompressed (codec_id == 0) without caching — fix that here so
// subsequent per-record slices hit the cache regardless of layout.
- (NSData *)_byteChannelFullNamed:(NSString *)name
{
    NSUInteger total = [self _totalBaseCount];
    if (total == 0) return [NSData data];
    NSData *cached = _decodedByteChannels[name];
    if (cached) return cached;
    NSError *err = nil;
    if (_storedOf != nil) {
        // M103: gather the stored channel into input order.
        NSData *stored = [self _blockByteChannelSliceNamed:name offset:0 count:total error:&err];
        TTIOGenomicIndex *sidx = [self _storedIndex];
        if (!stored || !sidx || stored.length < total) return [NSData data];
        NSMutableData *out = [NSMutableData dataWithLength:total];
        uint8_t *dst = (uint8_t *)out.mutableBytes;
        const uint8_t *src = (const uint8_t *)stored.bytes;
        const uint32_t *so = (const uint32_t *)_storedOf.bytes;
        NSUInteger n = _storedOf.length / sizeof(uint32_t), at = 0;
        for (NSUInteger i = 0; i < n; i++) {
            NSUInteger s = so[i];
            uint32_t len = [sidx lengthAt:s];
            uint64_t off = [sidx offsetAt:s];
            if (at + len > total || off + len > total) return [NSData data];
            memcpy(dst + at, src + off, len);
            at += len;
        }
        _decodedByteChannels[name] = out;
        return out;
    }
    NSData *full = [self byteChannelSliceNamed:name
                                         offset:0
                                          count:total
                                          error:&err];
    if (!full) return [NSData data];
    if (!_decodedByteChannels[name]) {
        _decodedByteChannels[name] = full;
    }
    return full;
}

- (NSArray<NSString *> *)allReadNames
{
    if (_blockTable) {
        NSMutableArray<NSString *> *all = [NSMutableArray arrayWithCapacity:(NSUInteger)_blockTable.readCount];
        for (NSUInteger b = 0; b < _blockTable.count; b++) {
            TTIOGenomicRun *view = [self _blockView:b error:NULL];
            if (!view) return @[];
            [all addObjectsFromArray:[view allReadNames]];
        }
        return _storedOf ? [self _inInputOrder:all] : all;
    }
    NSUInteger n = [self index].count;
    if (n == 0) return @[];
    // Touch index 0 to trigger the one-shot v2 decode + cache.
    NSError *err = nil;
    [self readNameAtIndex:0 error:&err];
    if (_decodedReadNames) return [_decodedReadNames copy];
    // Compound / uncompressed fallback path — rare under v1.0.
    NSMutableArray<NSString *> *out = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        NSString *s = [self readNameAtIndex:i error:NULL] ?: @"";
        [out addObject:s];
    }
    return [out copy];
}

@end
