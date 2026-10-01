/*
 * TTIOSamTags.m — SAM optional fields codec (codec id 18, M101).
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "TTIOSamTags.h"

#if __has_include(<ttio_rans.h>)
#include <ttio_rans.h>
#define TTIO_HAS_NATIVE_RANS 1
#else
#define TTIO_HAS_NATIVE_RANS 0
#endif

NSString *const TTIOSamTagsErrorDomain = @"global.thalion.ttio.SamTags";

static NSError *_stgError(NSInteger code, NSString *msg) {
    return [NSError errorWithDomain:TTIOSamTagsErrorDomain code:code
                           userInfo:@{NSLocalizedDescriptionKey: msg}];
}

@implementation TTIOSamTagsContext

- (BOOL)derives {
    for (id r in _references) {
        if ([r isKindOfClass:[NSData class]]) return YES;
    }
    return NO;
}

@end

#if TTIO_HAS_NATIVE_RANS
static NSString *_stgErrorMessage(int rc) {
    switch (rc) {
        case -1: return @"invalid parameters (a NUL in the tag text, or bad offsets)";
        case -2: return @"out of memory in native code";
        case -3: return @"corrupt SAM_TAGS blob";
        default: return [NSString stringWithFormat:@"native error %d", rc];
    }
}

/* Concatenate strings as UTF-8 into buf with uint64 offsets[n+1]. */
static void _stgJoin(NSArray<NSString *> *values, NSMutableData *buf,
                     NSMutableData *offsets) {
    NSUInteger n = values.count;
    [offsets setLength:(n + 1) * sizeof(uint64_t)];
    uint64_t *off = (uint64_t *)[offsets mutableBytes];
    off[0] = 0;
    for (NSUInteger i = 0; i < n; i++) {
        NSString *s = values[i];
        NSUInteger len = [s lengthOfBytesUsingEncoding:NSUTF8StringEncoding];
        if (len) {
            NSUInteger at = buf.length;
            [buf setLength:at + len];
            memcpy((uint8_t *)[buf mutableBytes] + at, [s UTF8String], len);
        }
        off[i + 1] = off[i] + len;
    }
    // The kernel wants a non-NULL buffer even when every value is empty.
    if (buf.length == 0) [buf setLength:1];
}

/* Holds the arrays a ttio_sam_tags_ctx points into. */
typedef struct {
    ttio_sam_tags_ctx ctx;
    const uint8_t **refPtrs;
    uint64_t *refLens;
} _stgNativeCtx;

static BOOL _stgBuildCtx(_stgNativeCtx *nc, NSUInteger n,
                         TTIOSamTagsContext *c, NSMutableArray *keep,
                         NSError **error) {
    memset(nc, 0, sizeof(*nc));
    nc->ctx.n_reads = (uint64_t)n;
    if (c == nil || ![c derives]) return YES;
    if (n && (c.sequences == nil || c.seqOffsets == nil || c.cigars == nil
              || c.positions == nil || c.chromIds == nil)) {
        if (error) *error = _stgError(-1,
            @"MD/NM derivation needs sequences, cigars, positions and chrom_ids");
        return NO;
    }
    if (c.seqOffsets.length != (n + 1) * sizeof(uint64_t)
        || c.cigars.count != n
        || c.positions.length != n * sizeof(int64_t)
        || c.chromIds.length != n * sizeof(uint16_t)) {
        if (error) *error = _stgError(-1, @"SAM_TAGS context array length mismatch");
        return NO;
    }
    NSData *seq = c.sequences.length ? c.sequences
                                     : [NSMutableData dataWithLength:1];
    NSMutableData *cig = [NSMutableData data];
    NSMutableData *cigOff = [NSMutableData data];
    _stgJoin(c.cigars, cig, cigOff);
    NSData *pos = c.positions.length ? c.positions : [NSMutableData dataWithLength:8];
    NSData *chr = c.chromIds.length ? c.chromIds : [NSMutableData dataWithLength:2];
    [keep addObjectsFromArray:@[seq, cig, cigOff, pos, chr]];

    NSUInteger nRefs = c.references.count;
    NSMutableData *ptrs = [NSMutableData dataWithLength:nRefs * sizeof(uint8_t *)];
    NSMutableData *lens = [NSMutableData dataWithLength:nRefs * sizeof(uint64_t)];
    [keep addObjectsFromArray:@[ptrs, lens]];
    const uint8_t **p = (const uint8_t **)[ptrs mutableBytes];
    uint64_t *l = (uint64_t *)[lens mutableBytes];
    for (NSUInteger i = 0; i < nRefs; i++) {
        id r = c.references[i];
        if ([r isKindOfClass:[NSData class]] && [(NSData *)r length] > 0) {
            [keep addObject:r];
            p[i] = (const uint8_t *)[(NSData *)r bytes];
            l[i] = (uint64_t)[(NSData *)r length];
        } else {
            p[i] = NULL;
            l[i] = 0;
        }
    }
    nc->ctx.sequences = (const uint8_t *)seq.bytes;
    nc->ctx.seq_offsets = (const uint64_t *)c.seqOffsets.bytes;
    nc->ctx.cigars = (const uint8_t *)cig.bytes;
    nc->ctx.cigar_offsets = (const uint64_t *)cigOff.bytes;
    nc->ctx.positions = (const int64_t *)pos.bytes;
    nc->ctx.chrom_ids = (const uint16_t *)chr.bytes;
    nc->ctx.refs = p;
    nc->ctx.ref_lengths = l;
    nc->ctx.n_refs = (uint32_t)nRefs;
    return YES;
}
#endif

@implementation TTIOSamTags

+ (BOOL)nativeAvailable {
    return TTIO_HAS_NATIVE_RANS ? YES : NO;
}

+ (nullable NSData *)encodeTags:(NSArray<NSString *> *)tags
                        context:(nullable TTIOSamTagsContext *)context
                          error:(NSError **)error {
#if !TTIO_HAS_NATIVE_RANS
    if (error) *error = _stgError(-100, @"libttio_rans not linked");
    return nil;
#else
    NSUInteger n = tags.count;
    // Precise lifetime: the native call only sees interior pointers.
    __attribute__((objc_precise_lifetime)) NSMutableArray *keep = [NSMutableArray array];
    _stgNativeCtx nc;
    if (!_stgBuildCtx(&nc, n, context, keep, error)) return nil;
    __attribute__((objc_precise_lifetime)) NSMutableData *buf = [NSMutableData data];
    __attribute__((objc_precise_lifetime)) NSMutableData *off = [NSMutableData data];
    _stgJoin(tags, buf, off);
    uint8_t *out = NULL;
    size_t outLen = 0;
    int rc = ttio_sam_tags_encode(&nc.ctx, (const uint8_t *)buf.bytes,
                                  (const uint64_t *)off.bytes, &out, &outLen);
    if (rc != 0) {
        if (error) *error = _stgError(rc, _stgErrorMessage(rc));
        return nil;
    }
    NSData *blob = [NSData dataWithBytes:out length:outLen];
    ttio_sam_tags_free(out);
    return blob;
#endif
}

+ (nullable NSArray<NSString *> *)decodeData:(NSData *)blob
                                      nReads:(NSUInteger)nReads
                                     context:(nullable TTIOSamTagsContext *)context
                                       error:(NSError **)error {
#if !TTIO_HAS_NATIVE_RANS
    if (error) *error = _stgError(-100, @"libttio_rans not linked");
    return nil;
#else
    // Precise lifetime: the native call only sees interior pointers.
    __attribute__((objc_precise_lifetime)) NSMutableArray *keep = [NSMutableArray array];
    _stgNativeCtx nc;
    if (!_stgBuildCtx(&nc, nReads, context, keep, error)) return nil;
    __attribute__((objc_precise_lifetime)) NSData *src =
        blob.length ? blob : [NSMutableData dataWithLength:1];
    __attribute__((objc_precise_lifetime)) NSMutableData *offData =
        [NSMutableData dataWithLength:(nReads + 1) * sizeof(uint64_t)];
    uint64_t *off = (uint64_t *)[offData mutableBytes];
    uint8_t *out = NULL;
    int rc = ttio_sam_tags_decode(&nc.ctx, (const uint8_t *)src.bytes, blob.length,
                                  &out, off);
    if (rc != 0) {
        if (error) *error = _stgError(rc, _stgErrorMessage(rc));
        return nil;
    }
    NSMutableArray<NSString *> *result = [NSMutableArray arrayWithCapacity:nReads];
    for (NSUInteger i = 0; i < nReads; i++) {
        uint64_t a = off[i], b = off[i + 1];
        NSString *s = (b > a)
            ? [[NSString alloc] initWithBytes:out + a length:(NSUInteger)(b - a)
                                     encoding:NSUTF8StringEncoding]
            : @"";
        if (s == nil) {
            ttio_sam_tags_free(out);
            if (error) *error = _stgError(-3, @"SAM_TAGS text is not valid UTF-8");
            return nil;
        }
        [result addObject:s];
    }
    ttio_sam_tags_free(out);
    return result;
#endif
}

@end
