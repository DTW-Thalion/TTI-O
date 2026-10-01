/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef TTIO_GENOMIC_STREAM_WRITER_H
#define TTIO_GENOMIC_STREAM_WRITER_H

#import <Foundation/Foundation.h>
#import "Providers/TTIOStorageProtocols.h"
#import "ValueClasses/TTIOEnums.h"

@class TTIOWrittenGenomicRun;
@class TTIOAlignedRead;
@class TTIOProvenanceRecord;
@class TTIOCompoundField;

NS_ASSUME_NONNULL_BEGIN

/** Run-level options of a streamed genomic run. Python:
 *  <code>GenomicStreamWriter</code> keyword arguments; Java:
 *  <code>GenomicStreamWriter.Options</code>. */
@interface TTIOGenomicStreamWriterOptions : NSObject <NSCopying>
@property (nonatomic) TTIOAcquisitionMode acquisitionMode;
@property (nonatomic, copy, nullable) NSString *referenceUri;
@property (nonatomic, copy, nullable) NSString *platform;
@property (nonatomic, copy, nullable) NSString *sampleName;
/** Persisted as the <code>@read_role</code> run attribute when
 *  non-nil; see <code>TTIOWrittenGenomicRun.readRole</code>. */
@property (nonatomic, copy, nullable) NSString *readRole;
/** REF_DIFF_V2 slice byte budget; 0 = the fixed 10,000-read rule.
 *  See <code>TTIOWrittenGenomicRun.refDiffSliceBytes</code>. */
@property (nonatomic) unsigned long long refDiffSliceBytes;
@property (nonatomic, copy, nullable) NSDictionary<NSString *, NSData *> *referenceChromSeqs;
@property (nonatomic) BOOL embedReference;
/** Reads per block; default 1 000 000. */
@property (nonatomic) NSUInteger blockReads;
/** Sequence bytes per block; default 256 MiB. */
@property (nonatomic) unsigned long long blockBytes;
@property (nonatomic) BOOL optDisableQualitiesV5;
@property (nonatomic, copy) NSDictionary<NSString *, NSNumber *> *signalCodecOverrides;
@property (nonatomic) TTIOCompression signalCompression;
@property (nonatomic) BOOL optLegacyWholeChannel;
@property (nonatomic, copy) NSArray<TTIOProvenanceRecord *> *provenanceRecords;
/** Worker threads for block encode (0 = TTIO_THREADS, else cores minus 2;
 *  1 = the serial path). With more than one, completed blocks encode on a
 *  queue and are written in order by the caller's thread, at most
 *  threads + 1 in flight. The file is byte for byte the one thread's. */
@property (nonatomic) NSUInteger threads;
/** Pipeline byte budget (0 = TTIO_MEMORY_BUDGET, else
 *  max(1 GiB, min(threads x blockBytes x 16, physical memory / 2))).
 *  The writer stalls block submission while its in-flight estimate
 *  exceeds half of it. */
@property (nonatomic) unsigned long long memoryBudgetBytes;

/** Defaults: zlib, no overrides, default block policy. */
+ (instancetype)defaultOptions;
/** The run-level metadata of <code>run</code>, default block policy. */
+ (instancetype)optionsFromRun:(TTIOWrittenGenomicRun *)run;
@end

/** Writes one genomic run as blocks_v1 (format-spec 10.12) with bounded
 *  memory. Reads are buffered until a block is full (blockReads reads or
 *  blockBytes sequence bytes, whichever first; a block never spans two
 *  chromosomes), encoded through TTIOGenomicBlocks and appended to
 *  extendable per-channel datasets; blocks/index records where each
 *  block's blob lives. Python: <code>ttio.genomic.GenomicStreamWriter</code>;
 *  Java: <code>GenomicStreamWriter</code>. */
@interface TTIOGenomicStreamWriter : NSObject

/** The layout attribute value, <code>blocks_v1</code>. */
+ (NSString *)layout;
/** Chunk of the unfiltered channel datasets (256 KiB). */
+ (NSUInteger)channelChunk;
/** Block index schema, in the column order of format-spec 10.12.2:
 *  the required channels only. */
+ (NSArray<TTIOCompoundField *> *)indexFields;
/** The schema of a run that carries the tags channel (M101): the tags
 *  (off, len, codec) triple follows the required columns. */
+ (NSArray<TTIOCompoundField *> *)indexFieldsWithTags;
/** The index schema for existing index rows: with the tags triple when
 *  the first row carries <code>tags_off</code>. */
+ (NSArray<TTIOCompoundField *> *)indexFieldsForRows:(NSArray<NSDictionary *> *)rows;

/** Append reads to run <code>runName</code> of the /study group
 *  <code>study</code>; the writer creates genomic_runs when absent and
 *  maintains its @_run_names. */
- (instancetype)initWithStudyGroup:(id<TTIOStorageGroup>)study
                           runName:(NSString *)runName
                           options:(nullable TTIOGenomicStreamWriterOptions *)options;

/** Append one read. */
- (BOOL)appendRead:(TTIOAlignedRead *)read error:(NSError **)error;
/** Append the reads of <code>batch</code>; its run-level metadata is
 *  ignored, the writer's options apply. */
- (BOOL)appendBatch:(TTIOWrittenGenomicRun *)batch error:(NSError **)error;
/** Encode and write the pending reads as one block. */
- (BOOL)flush:(NSError **)error;
/** Flush, write the name tables and provenance. Idempotent. */
- (BOOL)close:(NSError **)error;

@property (nonatomic, readonly) unsigned long long readCount;
@property (nonatomic, readonly) unsigned long long memoryBudgetBytes;
/** High-water mark of the in-flight byte estimate (tests). */
@property (nonatomic, readonly) unsigned long long maxInFlightBytesObserved;
@property (nonatomic, readonly) NSUInteger blockCount;

/** The pinned per-run qualities strategy hint: -1 until block 0 is
 *  written (or always -1 with TTIO_M94Z_EXHAUSTIVE=1), then
 *  TTIOM94ZHintV4Auto, 5, or 6.
 *
 *  A positive TTIO_M94Z_HINT pins it before block 0 instead, skipping
 *  the tune. TTIOM94ZHintV6 is reachable no other way, so measuring V6
 *  through the writer needs it. TTIO_M94Z_EXHAUSTIVE=1 wins. */
@property (nonatomic, readonly) NSInteger qualStrategyHint;
@property (nonatomic, readonly) NSUInteger threads;

/** Assign ids for every chromosome name the block introduces, in the order
 *  the block encoder assigns them (own names in read order, "*" included,
 *  then mate names that are not "", "*" or "="). */
+ (void)registerBlockChromosomes:(TTIOWrittenGenomicRun *)block
                         intoMap:(NSMutableDictionary<NSString *, NSNumber *> *)map;
@property (nonatomic, readonly, copy) TTIOGenomicStreamWriterOptions *options;

@end

NS_ASSUME_NONNULL_END

#endif
