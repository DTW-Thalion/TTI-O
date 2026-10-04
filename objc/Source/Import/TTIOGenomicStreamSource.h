/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef TTIO_GENOMIC_STREAM_SOURCE_H
#define TTIO_GENOMIC_STREAM_SOURCE_H

#import <Foundation/Foundation.h>
#import "Core/TTIOProgressSink.h"
#import "Providers/TTIOStorageProtocols.h"

@class TTIOWrittenGenomicRun;

NS_ASSUME_NONNULL_BEGIN

/** Hands each batch to <code>emit</code> in order; returns NO with
 *  <code>error</code> when the source fails, and stops as soon as
 *  <code>emit</code> returns NO (the writer failed). */
typedef BOOL (^TTIOGenomicBatchProducer)(BOOL (^emit)(TTIOWrittenGenomicRun *batch, NSError **error),
                                         NSError **error);

/** A genomic run delivered as a sequence of read batches, written
 *  through TTIOGenomicStreamWriter when the importing dataset is
 *  written. Python: <code>GenomicStreamSource</code>; Java:
 *  <code>GenomicStreamSource</code>. */
@interface TTIOGenomicStreamSource : NSObject

@property (nonatomic, readonly, copy) NSString *name;
@property (nonatomic, readonly, copy) TTIOGenomicBatchProducer batches;
/** Indexed FASTA opened lazily as the writer's reference; nil for none. */
@property (nonatomic, readonly, copy, nullable) NSString *referenceFasta;
@property (nonatomic, readonly) BOOL embedReference;
/** Block policy overrides; nil keeps the writer default. */
@property (nonatomic, readonly, copy, nullable) NSNumber *blockReads;
@property (nonatomic, readonly, copy, nullable) NSNumber *blockBytes;
@property (nonatomic, readonly) BOOL optLegacyWholeChannel;
/** M103: group the reads by sequence before blocking (no reference);
 *  see <code>TTIOGenomicStreamWriterOptions.groupReads</code>. */
@property (nonatomic, readonly) BOOL groupReads;

- (instancetype)initWithName:(NSString *)name
                     batches:(TTIOGenomicBatchProducer)batches
              referenceFasta:(nullable NSString *)referenceFasta
              embedReference:(BOOL)embedReference
                  blockReads:(nullable NSNumber *)blockReads
                  blockBytes:(nullable NSNumber *)blockBytes
       optLegacyWholeChannel:(BOOL)legacy;

/** A copy with a different block policy (grouping kept). */
- (instancetype)sourceWithBlockReads:(nullable NSNumber *)blockReads
                          blockBytes:(nullable NSNumber *)blockBytes
                              legacy:(BOOL)legacy;

/** A copy with a different block policy and grouping (M103). */
- (instancetype)sourceWithBlockReads:(nullable NSNumber *)blockReads
                          blockBytes:(nullable NSNumber *)blockBytes
                              legacy:(BOOL)legacy
                          groupReads:(BOOL)groupReads;

/** Run every batch through a TTIOGenomicStreamWriter on
 *  <code>study</code>; returns the reads written, or NSNotFound with
 *  <code>error</code>. */
- (NSUInteger)writeIntoStudy:(id<TTIOStorageGroup>)study
                    progress:(nullable TTIOProgressBlock)progress
                       error:(NSError **)error;

@end

NS_ASSUME_NONNULL_END

#endif
