/*
 * TTIOSamTags.h — SAM optional fields codec (codec id 18, M101).
 *
 * SAM columns 12+ (tab-joined as samtools prints them) coded as a
 * tag-line dictionary plus one column per tag key, with MD:Z and NM:i
 * recomputed from the reference when they match. The blob is written
 * as signal_channels/tags with @compression = 18. Spec:
 * docs/codecs/sam_tags.md.
 *
 * Direct link to the C library entries ttio_sam_tags_encode /
 * _decode / _free in libttio_rans (header at <ttio_rans.h>); there is
 * no pure-ObjC fallback, so without the library every call returns
 * nil + error.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef TTIO_SAM_TAGS_H
#define TTIO_SAM_TAGS_H

#import <Foundation/Foundation.h>

NS_ASSUME_NONNULL_BEGIN

extern NSString *const TTIOSamTagsErrorDomain;

/**
 * <p>MD/NM derivation context shared by encode and decode: the reads'
 * concatenated sequences with uint64 offsets (n+1), CIGAR strings,
 * int64 1-based positions (0 = unmapped) and uint16 chromosome ids
 * (0xFFFF = none), plus one reference sequence per chromosome id
 * (<code>NSNull</code> where there is none). With no non-null
 * reference, MD and NM are stored like any other tag and the other
 * fields are not read.</p>
 */
@interface TTIOSamTagsContext : NSObject
@property (nonatomic, strong, nullable) NSData *sequences;
@property (nonatomic, strong, nullable) NSData *seqOffsets;        // uint64, n+1
@property (nonatomic, copy,   nullable) NSArray<NSString *> *cigars;
@property (nonatomic, strong, nullable) NSData *positions;         // int64, n
@property (nonatomic, strong, nullable) NSData *chromIds;          // uint16, n
@property (nonatomic, copy,   nullable) NSArray *references;       // NSData or NSNull

/** YES when at least one entry of <code>references</code> is an
 *  NSData, i.e. MD/NM derivation is on. */
- (BOOL)derives;
@end

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Codecs/TTIOSamTags.h</p>
 *
 * <p>SAM_TAGS codec wrapper over the shared native kernel, so the
 * three SDKs write the same bytes.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.codecs.sam_tags</code></p>
 */
@interface TTIOSamTags : NSObject

/** @return YES iff <code>libttio_rans</code> is linked. */
+ (BOOL)nativeAvailable;

/**
 * Encode one block's per-read tag text to a SAM_TAGS blob.
 *
 * @param tags    One string per read (<code>@""</code> for none).
 * @param context Derivation context, or nil to store MD/NM verbatim.
 * @param error   Out-error on invalid input or native failure.
 * @return The blob, or nil with <code>*error</code> set.
 */
+ (nullable NSData *)encodeTags:(NSArray<NSString *> *)tags
                        context:(nullable TTIOSamTagsContext *)context
                          error:(NSError **)error;

/**
 * Decode a SAM_TAGS blob to per-read tag text (same context as
 * encode).
 *
 * @param blob    The SAM_TAGS stream.
 * @param nReads  Number of reads the blob holds.
 * @param context Derivation context, or nil.
 * @param error   Out-error on a corrupt blob or native failure.
 * @return One string per read, or nil with <code>*error</code> set.
 */
+ (nullable NSArray<NSString *> *)decodeData:(NSData *)blob
                                      nReads:(NSUInteger)nReads
                                     context:(nullable TTIOSamTagsContext *)context
                                       error:(NSError **)error;

@end

NS_ASSUME_NONNULL_END

#endif /* TTIO_SAM_TAGS_H */
