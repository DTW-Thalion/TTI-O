/*
 * TTIOSeqCm.h — read bases without a reference (codec id 19, M103).
 *
 * A context-mixing model (orders 11/16/24 with reverse-complement
 * training) and a binary arithmetic coder. The blob is written as
 * signal_channels/sequences with @compression = 19. Spec:
 * docs/codecs/seq_cm.md.
 *
 * Context-aware: decode takes the reads' lengths, which the run already
 * stores. Only the default parameters are used and the kernel picks the
 * model's table size from the base count, so the same reads always code
 * to the same bytes. Direct link to ttio_seq_cm_encode / _decode / _free
 * in libttio_rans (header at <ttio_rans.h>); there is no pure-ObjC
 * fallback, so without the library every call returns nil + error.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef TTIO_SEQ_CM_H
#define TTIO_SEQ_CM_H

#import <Foundation/Foundation.h>

NS_ASSUME_NONNULL_BEGIN

extern NSString *const TTIOSeqCmErrorDomain;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Codecs/TTIOSeqCm.h</p>
 *
 * <p>SEQ_CM codec wrapper over the shared native kernel, so the three
 * SDKs write the same bytes.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.codecs.seq_cm</code> &#183;
 * Java: <code>global.thalion.ttio.codecs.SeqCm</code></p>
 */
@interface TTIOSeqCm : NSObject

/** @return YES iff <code>libttio_rans</code> is linked. */
+ (BOOL)nativeAvailable;

/**
 * Encode reads' bases, back to back.
 *
 * @param sequences The bases of every read, concatenated.
 * @param lengths   One uint64 per read; they must sum to the base count.
 * @param error     Out-error on invalid input or native failure.
 * @return The SEQ_CM blob, or nil with <code>*error</code> set.
 */
+ (nullable NSData *)encodeSequences:(NSData *)sequences
                             lengths:(NSData *)lengths
                               error:(NSError **)error;

/**
 * Decode a SEQ_CM blob given the reads' lengths.
 *
 * @param blob    The SEQ_CM stream.
 * @param lengths One uint64 per read, as at encode.
 * @param error   Out-error on a corrupt blob, mismatched lengths or
 *                native failure.
 * @return The bases, or nil with <code>*error</code> set.
 */
+ (nullable NSData *)decodeData:(NSData *)blob
                        lengths:(NSData *)lengths
                          error:(NSError **)error;

@end

NS_ASSUME_NONNULL_END

#endif /* TTIO_SEQ_CM_H */
