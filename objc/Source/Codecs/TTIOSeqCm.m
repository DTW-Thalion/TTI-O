/*
 * TTIOSeqCm.m — read bases without a reference (codec id 19, M103).
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "TTIOSeqCm.h"

#if __has_include(<ttio_rans.h>)
#include <ttio_rans.h>
#define TTIO_HAS_NATIVE_RANS 1
#else
#define TTIO_HAS_NATIVE_RANS 0
#endif

NSString *const TTIOSeqCmErrorDomain = @"global.thalion.ttio.SeqCm";

static NSError *_scmError(NSInteger code, NSString *msg) {
    return [NSError errorWithDomain:TTIOSeqCmErrorDomain code:code
                           userInfo:@{NSLocalizedDescriptionKey: msg}];
}

#if TTIO_HAS_NATIVE_RANS
static NSString *_scmErrorMessage(int rc) {
    switch (rc) {
        case -1: return @"invalid parameters (lengths that do not match the bases or the blob)";
        case -2: return @"out of memory in native code";
        case -3: return @"corrupt SEQ_CM blob";
        default: return [NSString stringWithFormat:@"native error %d", rc];
    }
}
#endif

@implementation TTIOSeqCm

+ (BOOL)nativeAvailable {
    return TTIO_HAS_NATIVE_RANS ? YES : NO;
}

+ (nullable NSData *)encodeSequences:(NSData *)sequences
                             lengths:(NSData *)lengths
                               error:(NSError **)error {
#if !TTIO_HAS_NATIVE_RANS
    (void)sequences; (void)lengths;
    if (error) *error = _scmError(-100, @"libttio_rans not linked");
    return nil;
#else
    if (lengths.length % sizeof(uint64_t)) {
        if (error) *error = _scmError(-1, @"lengths must hold uint64 values");
        return nil;
    }
    uint64_t n = lengths.length / sizeof(uint64_t);
    // The kernel wants non-NULL buffers even for an empty block.
    __attribute__((objc_precise_lifetime)) NSData *seq =
        sequences.length ? sequences : [NSMutableData dataWithLength:1];
    __attribute__((objc_precise_lifetime)) NSData *lens =
        lengths.length ? lengths : [NSMutableData dataWithLength:sizeof(uint64_t)];
    const uint64_t *lp = (const uint64_t *)lens.bytes;
    uint64_t total = 0;
    for (uint64_t i = 0; i < n; i++) total += lp[i];
    if (total != sequences.length) {
        if (error) *error = _scmError(-1, [NSString stringWithFormat:
            @"lengths sum to %llu but there are %lu bases",
            (unsigned long long)total, (unsigned long)sequences.length]);
        return nil;
    }
    uint8_t *out = NULL;
    size_t outLen = 0;
    int rc = ttio_seq_cm_encode((const uint8_t *)seq.bytes, lp, n, NULL, &out, &outLen);
    if (rc != 0) {
        if (error) *error = _scmError(rc, _scmErrorMessage(rc));
        return nil;
    }
    NSData *blob = [NSData dataWithBytes:out length:outLen];
    ttio_seq_cm_free(out);
    return blob;
#endif
}

+ (nullable NSData *)decodeData:(NSData *)blob
                        lengths:(NSData *)lengths
                          error:(NSError **)error {
#if !TTIO_HAS_NATIVE_RANS
    (void)blob; (void)lengths;
    if (error) *error = _scmError(-100, @"libttio_rans not linked");
    return nil;
#else
    if (lengths.length % sizeof(uint64_t)) {
        if (error) *error = _scmError(-1, @"lengths must hold uint64 values");
        return nil;
    }
    uint64_t n = lengths.length / sizeof(uint64_t);
    __attribute__((objc_precise_lifetime)) NSData *src =
        blob.length ? blob : [NSMutableData dataWithLength:1];
    __attribute__((objc_precise_lifetime)) NSData *lens =
        lengths.length ? lengths : [NSMutableData dataWithLength:sizeof(uint64_t)];
    uint8_t *out = NULL;
    size_t outLen = 0;
    int rc = ttio_seq_cm_decode((const uint8_t *)src.bytes, blob.length,
                                (const uint64_t *)lens.bytes, n, &out, &outLen);
    if (rc != 0) {
        if (error) *error = _scmError(rc, _scmErrorMessage(rc));
        return nil;
    }
    NSData *seq = [NSData dataWithBytes:(outLen ? out : (const uint8_t *)"") length:outLen];
    ttio_seq_cm_free(out);
    return seq;
#endif
}

@end
