#ifndef TTIO_GENOMIC_RUN_H
#define TTIO_GENOMIC_RUN_H

#import <Foundation/Foundation.h>
#import "Protocols/TTIOIndexable.h"
#import "Protocols/TTIORun.h"
#import "ValueClasses/TTIOEnums.h"

@class TTIOAlignedRead;
@class TTIOGenomicIndex;
@class TTIOProvenanceRecord;
@class TTIOCodecContext;
@protocol TTIOStorageGroup;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Conforms To:</em> TTIOIndexable, TTIORun</p>
 * <p><em>Declared In:</em> Genomics/TTIOGenomicRun.h</p>
 *
 * <p>Lazy view over one
 * <code>/study/genomic_runs/&lt;name&gt;/</code> group.
 * Materialises <code>TTIOAlignedRead</code> objects on demand from
 * the signal channels. The <code>TTIOGenomicIndex</code> is loaded
 * eagerly at open time for cheap filtering and offset lookups; the
 * heavy signal channels (sequences, qualities, plus the inline
 * codec channels) stay lazy on disk.</p>
 *
 * <p>Genomic analogue of <code>TTIOAcquisitionRun</code>; both
 * conform to <code>TTIORun</code> so cross-modality code can
 * iterate uniformly via <code>-objectAtIndex:</code> /
 * <code>-count</code>.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.genomic_run.GenomicRun</code><br/>
 * Java: <code>global.thalion.ttio.genomics.GenomicRun</code></p>
 */
@interface TTIOGenomicRun : NSObject <TTIOIndexable, TTIORun>

/** Run identifier as stored in the .tio file (e.g.
 *  <code>@"genomic_0001"</code>). */
@property (readonly, copy) NSString *name;

/** Acquisition mode (typically
 *  <code>TTIOAcquisitionModeGenomicWGS</code> or
 *  <code>WES</code>). */
@property (readonly) TTIOAcquisitionMode acquisitionMode;

/** Omics modality identifier (typically
 *  <code>@"genomic_sequencing"</code>). */
@property (readonly, copy) NSString *modality;

/** URI of the reference genome (e.g. <code>@"GRCh38.p14"</code>). */
@property (readonly, copy) NSString *referenceUri;

/** Sequencing platform (e.g. <code>@"ILLUMINA"</code>). */
@property (readonly, copy) NSString *platform;

/** Assembly read role from the <code>@read_role</code> run attribute
 *  (M97): <code>hifi</code>, <code>ont_ul</code>, <code>hic_r1</code>,
 *  <code>hic_r2</code>, <code>parental_maternal</code>,
 *  <code>parental_paternal</code>, or <code>illumina_wgs</code>
 *  (other strings are stored unchecked). nil when the attribute is
 *  absent. Python: <code>read_role</code>; Java:
 *  <code>getReadRole()</code>. */
@property (readonly, copy, nullable) NSString *readRole;

/** Sample identifier. */
@property (readonly, copy) NSString *sampleName;

/** Per-read index. Loaded at open for the whole-channel layout, on
 *  first access under blocks_v1. */
@property (readonly, strong) TTIOGenomicIndex *index;

/** <code>blocks_v1</code> or <code>whole</code> (the v1.8 whole-channel
 *  layout). */
@property (readonly, copy) NSString *layout;

/** Number of blocks; 1 for a whole-channel run. */
@property (readonly) NSUInteger blockCount;

/** Run-level chromosome name table (genomic_index/chromosome_names),
 *  read without loading the per-read arrays. */
- (NSArray<NSString *> *)chromosomeNames;

/** Reads [start, stop) in order, holding at most one decoded block at a
 *  time under blocks_v1. The block receives each read and may set
 *  <code>*stop</code> to end early. Returns NO with <code>error</code>
 *  when a read fails to decode. */
/** As -iterReadsFrom:to:error:usingBlock: with the next `threads` blocks
 *  decoding ahead on a pool (0 = TTIO_THREADS; <= 1 is the serial walk). */
- (BOOL)iterReadsFrom:(NSUInteger)start
                   to:(NSUInteger)stop
              threads:(NSUInteger)threads
                error:(NSError **)error
           usingBlock:(void (^)(TTIOAlignedRead *read, NSUInteger index, BOOL *stop))block;

- (BOOL)iterReadsFrom:(NSUInteger)start
                   to:(NSUInteger)stop
                error:(NSError **)error
           usingBlock:(void (^)(TTIOAlignedRead *read, NSUInteger index, BOOL *stop))block;

/** One call per decoded block, on the pool.
 *
 *  <code>block</code> is called from several threads at once and in no
 *  particular order, and must be safe to call that way. That relaxed
 *  ordering is the whole difference from -iterReadsFrom:to:threads:,
 *  whose in-order delivery on the caller's thread is what bounds it.
 *
 *  <code>view</code> is a run over one block's reads and is valid only
 *  for the duration of the call. Read record <code>k</code> of the
 *  delivered range as <code>[view readAtIndex:viewStart + k]</code>;
 *  its index in the whole run is <code>firstRead + k</code>. The two
 *  differ whenever a range starts part-way into a block, so a caller
 *  that indexes the view by <code>firstRead</code> reads the wrong
 *  records. */
- (BOOL)iterBlocksFrom:(NSUInteger)start
                    to:(NSUInteger)stop
               threads:(NSUInteger)threads
                 error:(NSError **)error
            usingBlock:(void (^)(TTIOGenomicRun *view, NSUInteger viewStart,
                                 NSUInteger firstRead, NSUInteger nReads,
                                 BOOL *stop))block;

/** Drop the cached block view and its memory store. */
- (void)close;

/** @return Number of reads in the run. */
- (NSUInteger)readCount;

/** @return <code>readCount</code> (TTIOIndexable conformance). */
- (NSUInteger)count;

/**
 * Internal SPI (Task 5): builds and caches the
 * <code>TTIOCodecContext</code> handed to the codec-registry decode
 * adapters. All fields are derived from the run's index exactly as the
 * bespoke inline decode paths built them, so registry-routed decode is
 * byte-identical. Cached for the lifetime of this run instance.
 */
- (TTIOCodecContext *)_codecContext;

/**
 * @param index Zero-based read position.
 * @return The <code>TTIOAlignedRead</code> at <code>index</code>,
 *         or <code>nil</code> on error / out-of-range.
 */
- (id)objectAtIndex:(NSUInteger)index;

/**
 * @return Per-run provenance records in insertion order, read from
 *         <code>&lt;run&gt;/provenance/steps</code>. Empty array
 *         when the run carries no provenance.
 */
- (NSArray<TTIOProvenanceRecord *> *)provenanceChain;

/**
 * Materialises the read at <code>index</code>.
 *
 * @param index Zero-based read position.
 * @param error Out-parameter populated on failure.
 * @return The aligned read, or <code>nil</code> on failure.
 */
- (TTIOAlignedRead *)readAtIndex:(NSUInteger)index
                           error:(NSError **)error;

/**
 * Returns the read name at <code>index</code>. Decodes from the
 * NAME_TOKENIZED_V2 stream stored under the
 * <code>signal_channels/read_names</code> dataset; the decoded
 * list is materialised on first call and cached for the lifetime
 * of this run instance.
 *
 * @param index Zero-based read position.
 * @param error Out-parameter populated on failure.
 * @return The read name, or <code>nil</code> on failure.
 */
- (NSString *)readNameAtIndex:(NSUInteger)index
                        error:(NSError **)error;

/**
 * Returns the CIGAR string at <code>index</code>. Decodes the
 * cigars channel (rANS-O0 or rANS-O1, length-prefix-concat
 * varint+bytes); the decoded list is materialised on first call
 * and cached for the lifetime of this run instance.
 *
 * @param index Zero-based read position.
 * @param error Out-parameter populated on failure.
 * @return The CIGAR string, or <code>nil</code> on failure.
 */
- (NSString *)cigarAtIndex:(NSUInteger)index
                     error:(NSError **)error;

/**
 * @param chromosome Reference chromosome.
 * @param start      Inclusive lower bound on position.
 * @param end        Exclusive upper bound on position.
 * @return Reads on <code>chromosome</code> whose mapping position
 *         is in <code>[start, end)</code>.
 */
- (NSArray<TTIOAlignedRead *> *)readsInRegion:(NSString *)chromosome
                                        start:(int64_t)start
                                          end:(int64_t)end;

/**
 * Opens an existing
 * <code>/study/genomic_runs/&lt;name&gt;/</code> group. The caller
 * resolves the run group and passes it as <code>runGroup</code>.
 *
 * @param runGroup The run sub-group.
 * @param name     Run name.
 * @param error    Out-parameter populated on failure.
 * @return The opened run, or <code>nil</code> on failure.
 */
+ (instancetype)openFromGroup:(id<TTIOStorageGroup>)runGroup
                         name:(NSString *)name
                        error:(NSError **)error;

/** As <code>+openFromGroup:name:error:</code>; <code>resolver</code>,
 *  when given, is used for REF_DIFF decodes instead of one built from
 *  the run group's HDF5 file (block views live in a memory provider and
 *  share their parent's). */
+ (instancetype)openFromGroup:(id<TTIOStorageGroup>)runGroup
                         name:(NSString *)name
            referenceResolver:(nullable id)resolver
                        error:(NSError **)error;

/**
 * Returns the codec id (<code>TTIOCompression</code> value) declared
 * on the named signal channel via its <code>@compression</code>
 * attribute, or <code>0</code> (NONE) when the attribute is absent.
 * The transport writer probes this for sequences / qualities to
 * decide whether each per-AU slice should be re-encoded with that
 * codec on the wire.
 *
 * @param name Signal-channel dataset name.
 * @return Codec id, or <code>0</code> when no codec is declared.
 */
- (uint8_t)wireCompressionForChannel:(NSString *)name;

// ── Phase 2c-T verbatim v2 blob accessors ─────────────────────────
// Read raw on-disk codec blob bytes for the transport bulk-mode
// writer. Each returns nil when the channel is absent / not in the
// expected v2 layout (e.g. read_names with @compression != 15).

/**
 * Read the verbatim `mate_info/inline_v2` codec blob for transport
 * bulk-mode forwarding.
 *
 * Bypasses the codec-13 decode path so the transport writer can stream
 * the on-disk bytes directly without re-encoding. Returns `nil` when
 * the run does not use the inline-v2 layout.
 *
 * @return Raw inline-v2 blob bytes, or `nil` when the layout is absent.
 */
- (nullable NSData *)readMateInfoInlineV2BlobBytes;

/**
 * Read the `mate_info/chrom_names` sidecar table.
 *
 * Companion to `-readMateInfoInlineV2BlobBytes`: gives the transport
 * writer the chromosome-name dictionary it needs to reassemble mate
 * coordinates on the receiver side. Returns an empty array when the
 * table is missing.
 *
 * @return Ordered chromosome names, or `@[]` when the sidecar is absent.
 */
- (NSArray<NSString *> *)readMateInfoChromNamesTable;

/**
 * Read the verbatim `read_names` codec blob.
 *
 * Returns the on-disk bytes only when the channel is stored with
 * `@compression == NAME_TOKENIZED_V2 (15)`; returns `nil` for any
 * other layout so the transport writer can fall through to its
 * encode path.
 *
 * @return Raw NAME_TOKENIZED_V2 blob bytes, or `nil` when the channel
 *         is not stored in that layout.
 */
- (nullable NSData *)readNameTokV2BlobBytes;

/**
 * Read the verbatim `sequences/refdiff_v2` codec blob.
 *
 * Returns the on-disk bytes only when `sequences` is a group
 * containing the `refdiff_v2` dataset; returns `nil` for any other
 * layout.
 *
 * @return Raw refdiff-v2 blob bytes, or `nil` when the layout is absent.
 */
- (nullable NSData *)readRefDiffV2BlobBytes;

// ── Bulk accessors for hot serialization paths ────────────────────
//
// The per-read accessors (-objectAtIndex:, -readAtIndex:error:,
// -readNameAtIndex:error:, etc.) materialise a fresh
// TTIOAlignedRead and slice every channel on every call. For
// serialization workloads that touch every byte sequentially
// (TTIOFastqWriter, TTIOFastaWriter, TTIOTransportWriter at full-
// corpus scale) pre-fetching the whole channel once and slicing
// in-memory is dramatically faster — Python's FastqWriter saw a
// 24× speedup at 1M reads from this exact pattern.

/**
 * Return the full `signal_channels/sequences` byte buffer.
 *
 * Bulk accessor for hot serialization paths (FASTQ / FASTA / transport
 * writers). Decoded once and cached for codec-compressed channels;
 * read once for uncompressed channels. Returns an empty `NSData` for
 * zero-read runs.
 *
 * @return Concatenated sequence bytes for every read in the run, in
 *         record order.
 */
- (NSData *)wholeSequencesData;

/**
 * Return the full `signal_channels/qualities` byte buffer.
 *
 * Same caching semantics as `-wholeSequencesData`. Returns an empty
 * `NSData` for zero-read runs.
 *
 * @return Concatenated Phred+33 quality bytes for every read in the
 *         run, in record order.
 */
- (NSData *)wholeQualitiesData;

/**
 * Return the full read-names list.
 *
 * Forces the one-shot NAME_TOKENIZED_V2 decode + cache when the
 * channel is tokenised, so subsequent per-read lookups become O(1)
 * array indexing rather than re-decoding the blob.
 *
 * @return Read names in record order; empty array for zero-read runs.
 */
- (NSArray<NSString *> *)allReadNames;

// ── M101 SAM optional fields ─────────────────────────────────────

/** YES when the run carries the SAM tags channel
 *  (<code>signal_channels/tags</code>, codec SAM_TAGS) in either
 *  layout. */
- (BOOL)hasTagsChannel;

/**
 * Every read's SAM tag text (columns 12+, tab-joined), decoded once
 * and cached; <code>@[]</code> when the run has no tags channel.
 *
 * @param error Out-parameter populated on decode failure.
 * @return One string per read, <code>@[]</code>, or nil on failure.
 */
- (NSArray<NSString *> *)allTagsWithError:(NSError **)error;

@end

#endif
