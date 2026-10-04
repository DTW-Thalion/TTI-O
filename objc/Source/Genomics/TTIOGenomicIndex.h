#ifndef TTIO_GENOMIC_INDEX_H
#define TTIO_GENOMIC_INDEX_H

#import <Foundation/Foundation.h>

@protocol TTIOStorageGroup;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Genomics/TTIOGenomicIndex.h</p>
 *
 * <p>Per-read offsets, lengths, positions, mapping qualities, flags,
 * and chromosome strings for one <code>TTIOGenomicRun</code>. Held
 * in memory as <code>NSData</code> buffers backing typed C arrays;
 * exposed via typed accessors. Loaded eagerly when the run is
 * opened; persisted under
 * <code>/study/genomic_runs/&lt;name&gt;/genomic_index/</code>.</p>
 *
 * <p>Genomic analogue of <code>TTIOSpectrumIndex</code>. Range
 * queries (region, unmapped, flag mask) operate on the in-memory
 * arrays without touching the heavy signal channels.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.genomic_index.GenomicIndex</code><br/>
 * Java: <code>global.thalion.ttio.genomics.GenomicIndex</code></p>
 */
@interface TTIOGenomicIndex : NSObject

/** Number of indexed reads. */
@property (readonly) NSUInteger count;

/** @return Starting byte offset of read <code>index</code>'s
 *          sequence/quality data. */
- (uint64_t)offsetAt:(NSUInteger)index;

/** @return Number of bases in read <code>index</code>. */
- (uint32_t)lengthAt:(NSUInteger)index;

/** @return Zero-based mapping position. */
- (int64_t)positionAt:(NSUInteger)index;

/** @return SAM mapping quality (Phred-scaled). */
- (uint8_t)mappingQualityAt:(NSUInteger)index;

/** @return SAM flags bit-field. */
- (uint32_t)flagsAt:(NSUInteger)index;

/** @return Reference chromosome name. */
- (NSString *)chromosomeAt:(NSUInteger)index;

/**
 * Region query.
 *
 * @param chromosome Reference chromosome.
 * @param start      Inclusive lower bound on position.
 * @param end        Exclusive upper bound on position.
 * @return Indices on <code>chromosome</code> with
 *         <code>start &lt;= position &lt; end</code>.
 */
- (NSIndexSet *)indicesForRegion:(NSString *)chromosome
                           start:(int64_t)start
                             end:(int64_t)end;

/** @return Indices where <code>(flags &amp; 0x4) != 0</code>
 *          (unmapped reads). */
- (NSIndexSet *)indicesForUnmapped;

/**
 * @param flagMask SAM-flag bit mask to test.
 * @return Indices where
 *         <code>(flags &amp; flagMask) != 0</code>.
 */
- (NSIndexSet *)indicesForFlag:(uint32_t)flagMask;

#pragma mark - Construction

/**
 * Designated initialiser.
 *
 * @param offsets          uint64 byte offsets per read.
 * @param lengths          uint32 base counts per read.
 * @param chromosomes      Per-read chromosome names.
 * @param positions        int64 mapping positions.
 * @param mappingQualities uint8 mapping qualities.
 * @param flags            uint32 SAM flags.
 * @return An initialised index.
 */
- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                    chromosomes:(NSArray<NSString *> *)chromosomes
                      positions:(NSData *)positions
               mappingQualities:(NSData *)mappingQualities
                          flags:(NSData *)flags;

#pragma mark - Storage round-trip

/**
 * Writes the index under
 * <code>group/genomic_index/</code>.
 *
 * @param group Destination parent group.
 * @param error Out-parameter populated on failure.
 * @return <code>YES</code> on success.
 */
- (BOOL)writeToGroup:(id<TTIOStorageGroup>)group
                error:(NSError **)error;

/**
 * Writes the index; nameToId, when given, is the shared chromosome id
 * map of a blocks_v1 run: existing entries keep their ids and new names
 * are appended in place.
 */
- (BOOL)writeToGroup:(id<TTIOStorageGroup>)group
             nameToId:(NSMutableDictionary<NSString *, NSNumber *> *)nameToId
                error:(NSError **)error;

/** Names of nameToId sorted by id. */
+ (NSArray<NSString *> *)namesInIdOrder:(NSDictionary<NSString *, NSNumber *> *)nameToId;

/**
 * Reads the index from <code>group/genomic_index/</code>.
 *
 * @param group Source parent group.
 * @param error Out-parameter populated on failure.
 * @return The materialised index, or <code>nil</code> on failure.
 */
+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)group
                        error:(NSError **)error;

/** M103: this index in another read order. <code>storedOf</code> holds
 *  one uint32 per read: entry i is the row of this index that becomes
 *  row i. Offsets are recomputed from the reordered lengths; the
 *  chromosome id table, when present, is kept. Used to present a
 *  <code>blocks_v1_grouped</code> run's index in input order. */
- (instancetype)indexPermutedBy:(NSData *)storedOf;

@end

/**
 * Synthesises per-record byte offsets from a <code>uint32</code>
 * lengths array. <code>offsets[i] = sum(lengths[0..i])</code>,
 * produced as a <code>uint64</code> <code>NSData</code> blob. Empty
 * input returns a zero-byte <code>NSData</code>.
 *
 * @param lengths Per-record uint32 lengths.
 * @return uint64 cumulative offsets.
 */
extern NSData *TTIOOffsetsFromLengths(NSData *lengths);

#endif
