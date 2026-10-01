#ifndef TTIO_ALIGNED_READ_H
#define TTIO_ALIGNED_READ_H

#import <Foundation/Foundation.h>

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Conforms To:</em> NSCopying</p>
 * <p><em>Declared In:</em> Genomics/TTIOAlignedRead.h</p>
 *
 * <p>One aligned sequencing read — the genomic analogue of
 * <code>TTIOMassSpectrum</code>. Immutable value object materialised
 * by <code>TTIOGenomicRun</code> from the signal channels under
 * <code>/study/genomic_runs/&lt;name&gt;/signal_channels/</code>.
 * The class itself does not perform any storage I/O.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.aligned_read.AlignedRead</code><br/>
 * Java: <code>global.thalion.ttio.genomics.AlignedRead</code></p>
 */
@interface TTIOAlignedRead : NSObject <NSCopying>

/** Read identifier as it appears in the source SAM/BAM. */
@property (nonatomic, readonly, copy) NSString *readName;

/** Reference chromosome name; <code>@"*"</code> for unmapped reads. */
@property (nonatomic, readonly, copy) NSString *chromosome;

/** Zero-based mapping position; <code>-1</code> for unmapped. */
@property (nonatomic, readonly) int64_t position;

/** SAM mapping quality (Phred-scaled). */
@property (nonatomic, readonly) uint8_t mappingQuality;

/** CIGAR string. */
@property (nonatomic, readonly, copy) NSString *cigar;

/** Read sequence (uppercase ACGT[N] / IUPAC). */
@property (nonatomic, readonly, copy) NSString *sequence;

/** Per-base Phred quality scores. */
@property (nonatomic, readonly, copy) NSData *qualities;

/** SAM flags bit-field. */
@property (nonatomic, readonly) uint32_t flags;

/** Mate's reference chromosome. */
@property (nonatomic, readonly, copy) NSString *mateChromosome;

/** Mate's mapping position. */
@property (nonatomic, readonly) int64_t matePosition;

/** Template length (TLEN). */
@property (nonatomic, readonly) int32_t templateLength;

/** M101: the SAM optional fields (columns 12+), tab-joined as
 *  samtools prints them; <code>@""</code> when the read has none. */
@property (nonatomic, readonly, copy) NSString *tags;

/** @return <code>YES</code> if the read is mapped (SAM flag 0x4
 *          not set). */
- (BOOL)isMapped;

/** @return <code>YES</code> if the read is part of a paired
 *          alignment (SAM flag 0x1). */
- (BOOL)isPaired;

/** @return <code>YES</code> if the read is reverse-complemented
 *          (SAM flag 0x10). */
- (BOOL)isReverse;

/** @return <code>YES</code> if this is a secondary alignment
 *          (SAM flag 0x100). */
- (BOOL)isSecondary;

/** @return <code>YES</code> if this is a supplementary alignment
 *          (SAM flag 0x800). */
- (BOOL)isSupplementary;

/** @return Number of bases in the sequence. */
- (NSUInteger)readLength;

/**
 * Designated initialiser.
 *
 * @return An initialised aligned read.
 */
- (instancetype)initWithReadName:(NSString *)readName
                      chromosome:(NSString *)chromosome
                        position:(int64_t)position
                  mappingQuality:(uint8_t)mappingQuality
                           cigar:(NSString *)cigar
                        sequence:(NSString *)sequence
                       qualities:(NSData *)qualities
                           flags:(uint32_t)flags
                  mateChromosome:(NSString *)mateChromosome
                    matePosition:(int64_t)matePosition
                  templateLength:(int32_t)templateLength;

/**
 * Initialiser carrying the read's SAM optional fields (M101); the
 * designated initialiser above delegates here with <code>@""</code>.
 *
 * @return An initialised aligned read.
 */
- (instancetype)initWithReadName:(NSString *)readName
                      chromosome:(NSString *)chromosome
                        position:(int64_t)position
                  mappingQuality:(uint8_t)mappingQuality
                           cigar:(NSString *)cigar
                        sequence:(NSString *)sequence
                       qualities:(NSData *)qualities
                           flags:(uint32_t)flags
                  mateChromosome:(NSString *)mateChromosome
                    matePosition:(int64_t)matePosition
                  templateLength:(int32_t)templateLength
                            tags:(NSString *)tags;
@end

#endif
