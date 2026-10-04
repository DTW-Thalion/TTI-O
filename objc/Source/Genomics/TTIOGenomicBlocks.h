/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef TTIO_GENOMIC_BLOCKS_H
#define TTIO_GENOMIC_BLOCKS_H

#import <Foundation/Foundation.h>

@class TTIOWrittenGenomicRun;
@class TTIOGenomicWriteContext;

NS_ASSUME_NONNULL_BEGIN

/** One block's encoded channels (format-spec 10.12). Dictionaries are
 *  keyed by channel name; an absent channel has an empty blob and codec
 *  0. Python: <code>ttio.genomic._blocks.BlockBlobs</code>; Java:
 *  <code>GenomicBlocks.BlockBlobs</code>. */
@interface TTIOBlockBlobs : NSObject
@property (nonatomic, readonly, copy) NSDictionary<NSString *, NSData *> *blobs;
@property (nonatomic, readonly, copy) NSDictionary<NSString *, NSNumber *> *codecs;
@property (nonatomic, readonly, copy) NSDictionary<NSString *, NSDictionary<NSString *, id> *> *extraAttrs;
@property (nonatomic, readonly) NSUInteger nReads;
@property (nonatomic, readonly) uint64_t nBases;
- (instancetype)initWithBlobs:(NSDictionary<NSString *, NSData *> *)blobs
                       codecs:(NSDictionary<NSString *, NSNumber *> *)codecs
                   extraAttrs:(NSDictionary<NSString *, NSDictionary<NSString *, id> *> *)extraAttrs
                       nReads:(NSUInteger)nReads
                       nBases:(uint64_t)nBases;
@end

/** Block encoder for the blocks_v1 genomic layout. A block is a run made
 *  of a contiguous range of reads; its channel blobs come from running
 *  the storage-path whole-channel writer against a memory-provider group
 *  and harvesting each channel dataset's bytes and @compression, so a
 *  block's blob is byte-identical to what a v1.8 write of those reads
 *  alone produces. Python: <code>ttio.genomic._blocks</code>; Java:
 *  <code>GenomicBlocks</code>. */
@interface TTIOGenomicBlocks : NSObject

/** Every blob channel of a block, in block-index column order:
 *  the required channels, then the optional ones (M101: tags). */
+ (NSArray<NSString *> *)blockChannels;

/** Blob channels every blocks_v1 run indexes: sequences, qualities,
 *  read_names, cigars, mate_info. */
+ (NSArray<NSString *> *)requiredBlockChannels;

/** Channels whose index columns exist only in runs that carry them
 *  (M101: tags). Their (off, len, codec) columns follow the required
 *  channels' codec columns; a reader treats missing columns as an
 *  empty channel (format-spec 10.12.2). */
+ (NSArray<NSString *> *)optionalBlockChannels;

/** Reads [start, stop) of <code>run</code> as a run of their own,
 *  offsets rebased to 0; run-level metadata shared. */
+ (TTIOWrittenGenomicRun *)sliceRun:(TTIOWrittenGenomicRun *)run
                               from:(NSUInteger)start
                                 to:(NSUInteger)stop;

/** The inverse of <code>+sliceRun:from:to:</code> for consecutive parts. */
+ (TTIOWrittenGenomicRun *)concatRuns:(NSArray<TTIOWrittenGenomicRun *> *)parts;

/** M103: the reads of <code>run</code> in the order <code>order</code>
 *  (uint32 per read; read <code>order[j]</code> at row j), offsets
 *  recomputed. Run-level metadata is shared, provenance dropped.
 *  Python: <code>_blocks.take_run</code>. */
+ (TTIOWrittenGenomicRun *)takeRun:(TTIOWrittenGenomicRun *)run order:(NSData *)order;

/** M103: the read-grouping permutation of <code>run</code> (uint32 per
 *  read; entry j is the input index of the read stored at row j). Reads
 *  are grouped by sequence with <code>ttio_seq_group</code> within each
 *  chromosome label, labels in first-seen order, so a block never spans
 *  two labels. nil with <code>*error</code> when the native library is
 *  missing or the kernel fails. Python: <code>_blocks.group_order</code>. */
+ (nullable NSData *)groupOrderOfRun:(TTIOWrittenGenomicRun *)run error:(NSError **)error;

/** Encode one block's channels through the whole-channel writer. The
 *  forced codecs of format-spec 10.12.3 apply: cigars RANS_ORDER0,
 *  qualities FQZCOMP_NX16_Z (RANS_ORDER0 when the block holds a
 *  zero-length read), sequences RANS_ORDER1 without a reference. */
+ (nullable TTIOBlockBlobs *)encodeBlock:(TTIOWrittenGenomicRun *)block
                                 context:(TTIOGenomicWriteContext *)ctx
                                   error:(NSError **)error;

@end

NS_ASSUME_NONNULL_END

#endif
