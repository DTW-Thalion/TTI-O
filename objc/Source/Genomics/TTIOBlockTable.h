/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef TTIO_BLOCK_TABLE_H
#define TTIO_BLOCK_TABLE_H

#import <Foundation/Foundation.h>
#import "Providers/TTIOStorageProtocols.h"

NS_ASSUME_NONNULL_BEGIN

/** The decoded <code>blocks/index</code> of a blocks_v1 run
 *  (format-spec 10.12.2): one row per block. Python:
 *  <code>_block_view.BlockTable</code>; Java: <code>BlockTable</code>. */
@interface TTIOBlockTable : NSObject

/** Read <code>runGroup/blocks/index</code>. */
+ (nullable instancetype)readFromRunGroup:(id<TTIOStorageGroup>)runGroup error:(NSError **)error;

/** Number of blocks. */
@property (nonatomic, readonly) NSUInteger count;
/** Reads in the run: the last block's read_start + n_reads. */
@property (nonatomic, readonly) unsigned long long readCount;
/** YES when the index carries the <code>&lt;ch&gt;_codec</code> columns. */
@property (nonatomic, readonly) BOOL hasCodecs;
/** The channels whose columns the index carries: the required five,
 *  plus <code>tags</code> when the run has the tags triple (M101).
 *  An optional channel the index lacks reads as empty ranges. */
@property (nonatomic, readonly, copy) NSArray<NSString *> *channels;
/** YES when the index carries the tags (off, len, codec) triple. */
@property (nonatomic, readonly) BOOL hasTags;

- (unsigned long long)readStartAt:(NSUInteger)block;
- (NSUInteger)nReadsAt:(NSUInteger)block;
- (unsigned long long)baseStartAt:(NSUInteger)block;
- (unsigned long long)nBasesAt:(NSUInteger)block;
- (unsigned long long)offsetOf:(NSString *)channel at:(NSUInteger)block;
- (unsigned long long)lengthOf:(NSString *)channel at:(NSUInteger)block;
/** Codec id of a channel in a block; 0 when the table has no codec columns. */
- (NSUInteger)codecOf:(NSString *)channel at:(NSUInteger)block;

/** The block holding read <code>i</code>, or NSNotFound when out of range. */
- (NSUInteger)blockForRead:(unsigned long long)i;

@end

NS_ASSUME_NONNULL_END

#endif
