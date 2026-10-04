/*
 * TTIOSeqGroup.m — read grouping for unaligned runs (M103).
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "TTIOSeqGroup.h"

#if __has_include(<ttio_rans.h>)
#include <ttio_rans.h>
#define TTIO_HAS_NATIVE_RANS 1
#else
#define TTIO_HAS_NATIVE_RANS 0
#endif

NSString *const TTIOSeqGroupErrorDomain = @"global.thalion.ttio.SeqGroup";

static NSError *_sgError(NSInteger code, NSString *msg) {
    return [NSError errorWithDomain:TTIOSeqGroupErrorDomain code:code
                           userInfo:@{NSLocalizedDescriptionKey: msg}];
}

@implementation TTIOSeqGroup

+ (BOOL)nativeAvailable {
    return TTIO_HAS_NATIVE_RANS ? YES : NO;
}

+ (nullable NSData *)groupSequences:(NSData *)sequences
                            lengths:(NSData *)lengths
                              names:(nullable NSArray<NSString *> *)names
                              error:(NSError **)error {
#if !TTIO_HAS_NATIVE_RANS
    (void)sequences; (void)lengths; (void)names;
    if (error) *error = _sgError(-100, @"read grouping requires libttio_rans");
    return nil;
#else
    if (lengths.length % sizeof(uint64_t)) {
        if (error) *error = _sgError(-1, @"lengths must hold uint64 values");
        return nil;
    }
    uint64_t n = lengths.length / sizeof(uint64_t);
    const uint64_t *lp0 = (const uint64_t *)lengths.bytes;
    uint64_t total = 0;
    for (uint64_t i = 0; i < n; i++) total += lp0[i];
    if (total != sequences.length) {
        if (error) *error = _sgError(-1, [NSString stringWithFormat:
            @"lengths sum to %llu but there are %lu bases",
            (unsigned long long)total, (unsigned long)sequences.length]);
        return nil;
    }
    if (n >= 0xFFFFFFFFull) {
        if (error) *error = _sgError(-1, @"read grouping holds at most 2^32 - 2 reads");
        return nil;
    }
    ttio_seq_group_params params;
    ttio_seq_group_default_params(&params);
    NSMutableData *nameBuf = [NSMutableData data];
    NSMutableData *nameOff = [NSMutableData dataWithLength:(NSUInteger)(n + 1) * sizeof(uint64_t)];
    if (names == nil) {
        params.flags &= (uint8_t)~TTIO_SEQ_GROUP_FLAG_MATES;
    } else {
        if (names.count != n) {
            if (error) *error = _sgError(-1, [NSString stringWithFormat:
                @"%lu names for %llu reads", (unsigned long)names.count, (unsigned long long)n]);
            return nil;
        }
        uint64_t *no = (uint64_t *)nameOff.mutableBytes;
        uint64_t acc = 0;
        for (uint64_t i = 0; i < n; i++) {
            NSData *b = [names[(NSUInteger)i] dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
            [nameBuf appendData:b];
            acc += b.length;
            no[i + 1] = acc;
        }
    }
    // The kernel wants non-NULL buffers even when empty.
    if (nameBuf.length == 0) [nameBuf setLength:1];
    __attribute__((objc_precise_lifetime)) NSData *seq =
        sequences.length ? sequences : [NSMutableData dataWithLength:1];
    __attribute__((objc_precise_lifetime)) NSData *lens =
        lengths.length ? lengths : [NSMutableData dataWithLength:sizeof(uint64_t)];
    NSMutableData *order = [NSMutableData dataWithLength:(NSUInteger)(n ? n : 1) * sizeof(uint32_t)];
    int rc = ttio_seq_group((const uint8_t *)seq.bytes, (const uint64_t *)lens.bytes, n,
                            (const uint8_t *)nameBuf.bytes, (const uint64_t *)nameOff.bytes,
                            &params, (uint32_t *)order.mutableBytes);
    if (rc == -2) {
        if (error) *error = _sgError(rc, @"out of memory in ttio_seq_group");
        return nil;
    }
    if (rc != 0) {
        if (error) *error = _sgError(rc, [NSString stringWithFormat:@"ttio_seq_group failed (rc=%d)", rc]);
        return nil;
    }
    [order setLength:(NSUInteger)n * sizeof(uint32_t)];
    return order;
#endif
}

@end
