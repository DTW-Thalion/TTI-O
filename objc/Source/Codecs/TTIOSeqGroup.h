/*
 * TTIOSeqGroup.h — read grouping for unaligned runs (M103).
 *
 * Returns the permutation that puts an unaligned run's reads from the
 * same place in the genome next to each other, so a blocks_v1 block
 * holds them together and SEQ_CM sees each read's overlap partners.
 * order[j] is the input index of the read stored at row j. The native
 * kernel breaks every tie on a total order, so the three SDKs compute
 * the same permutation. Direct link to ttio_seq_group in libttio_rans
 * (header at <ttio_rans.h>); there is no pure-ObjC fallback.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef TTIO_SEQ_GROUP_H
#define TTIO_SEQ_GROUP_H

#import <Foundation/Foundation.h>

NS_ASSUME_NONNULL_BEGIN

extern NSString *const TTIOSeqGroupErrorDomain;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Codecs/TTIOSeqGroup.h</p>
 *
 * <p>Wrapper over the shared native read-grouping kernel.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.codecs.seq_group</code> &#183;
 * Java: <code>global.thalion.ttio.codecs.SeqGroup</code></p>
 */
@interface TTIOSeqGroup : NSObject

/** @return YES iff <code>libttio_rans</code> is linked. */
+ (BOOL)nativeAvailable;

/**
 * The grouping order of reads whose bases are <code>sequences</code>
 * back to back, with the kernel's default parameters.
 *
 * @param sequences The bases of every read, concatenated.
 * @param lengths   One uint64 per read; they must sum to the base count.
 * @param names     One name per read (pairs mates by name; the kernel
 *                  strips a trailing /1 or /2), or nil to group without
 *                  them.
 * @param error     Out-error on invalid input or native failure.
 * @return uint32 per read (entry j is the input index of the read
 *         stored at row j), or nil with <code>*error</code> set.
 */
+ (nullable NSData *)groupSequences:(NSData *)sequences
                            lengths:(NSData *)lengths
                              names:(nullable NSArray<NSString *> *)names
                              error:(NSError **)error;

@end

NS_ASSUME_NONNULL_END

#endif /* TTIO_SEQ_GROUP_H */
