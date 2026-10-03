// TestM103SeqCm.m — M103: SEQ_CM (codec id 19), read bases without a
// reference.
//
// The wrapper round-trips every byte value and zero-length reads, codes
// the same reads to the same bytes, rejects mismatched lengths, and the
// registry adapter takes the lengths from the codec context. The
// blocks_v1 default is covered in TestGenomicBlocks and
// TestGenomicStreamWriter.
//
// Mirrors:
//   python/tests/test_m103_seq_cm.py
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#import <Foundation/Foundation.h>
#import "Testing.h"
#import "Codecs/TTIOSeqCm.h"
#import "Codecs/Registry/TTIOCodecRegistry.h"
#import "Codecs/Registry/TTIOCodecContext.h"
#import "Codecs/Registry/TTIOChannelPayload.h"
#import "Codecs/Registry/TTIODecodedChannel.h"
#import "Codecs/Registry/TTIOEncodedChannel.h"
#import "ValueClasses/TTIOEnums.h"

/* Overlapping reads from a random genome, with N and lower case, a
 * zero-length read, and every byte value in the last read. */
static void scmReads(NSMutableData *seq, NSMutableData *lengths, unsigned seed, NSUInteger n)
{
    srandom(seed);
    static const char acgt[] = "ACGT";
    NSMutableData *genome = [NSMutableData dataWithLength:20000];
    uint8_t *g = (uint8_t *)[genome mutableBytes];
    for (NSUInteger i = 0; i < 20000; i++) g[i] = (uint8_t)acgt[random() & 3];
    for (NSUInteger i = 0; i < n; i++) {
        uint64_t L = (i == 3) ? 0 : (uint64_t)(40 + random() % 120);
        NSUInteger s = (NSUInteger)(random() % (20000 - 200));
        NSMutableData *r = [NSMutableData dataWithBytes:g + s length:(NSUInteger)L];
        uint8_t *rb = (uint8_t *)[r mutableBytes];
        if (L && i % 7 == 0) rb[L / 2] = 'N';
        if (L && i % 11 == 0) rb[0] = 'a';
        if (i == n - 1) {
            for (int b = 0; b < 256; b++) { uint8_t v = (uint8_t)b; [r appendBytes:&v length:1]; }
            L += 256;
        }
        [seq appendData:r];
        [lengths appendBytes:&L length:sizeof L];
    }
}

void testM103SeqCm(void)
{
    if (![TTIOSeqCm nativeAvailable]) {
        PASS(YES, "M103 seq_cm: libttio_rans not linked, skipped");
        return;
    }
    NSMutableData *seq = [NSMutableData data], *lengths = [NSMutableData data];
    scmReads(seq, lengths, 3, 400);
    NSError *err = nil;
    NSData *blob = [TTIOSeqCm encodeSequences:seq lengths:lengths error:&err];
    PASS(blob != nil && blob.length >= 4 && memcmp(blob.bytes, "SQC1", 4) == 0,
         "M103 seq_cm: encode (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    NSData *back = [TTIOSeqCm decodeData:blob lengths:lengths error:&err];
    PASS([back isEqualToData:seq], "M103 seq_cm: round trip, every byte value and an empty read");
    PASS([[TTIOSeqCm encodeSequences:seq lengths:lengths error:&err] isEqualToData:blob],
         "M103 seq_cm: deterministic");

    NSData *none = [NSData data];
    NSData *emptyBlob = [TTIOSeqCm encodeSequences:none lengths:none error:&err];
    PASS(emptyBlob != nil && [[TTIOSeqCm decodeData:emptyBlob lengths:none error:&err] length] == 0,
         "M103 seq_cm: empty block");

    NSData *shortLens = [lengths subdataWithRange:NSMakeRange(0, lengths.length - sizeof(uint64_t))];
    err = nil;
    PASS([TTIOSeqCm encodeSequences:seq lengths:shortLens error:&err] == nil && err != nil,
         "M103 seq_cm: encode rejects lengths that do not sum to the bases");
    err = nil;
    PASS([TTIOSeqCm decodeData:blob lengths:shortLens error:&err] == nil && err != nil,
         "M103 seq_cm: decode rejects a read count that does not match the blob");

    // Registry: lengths come from the codec context.
    NSUInteger n = lengths.length / sizeof(uint64_t);
    const uint64_t *lp = (const uint64_t *)lengths.bytes;
    NSMutableArray *rl = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) [rl addObject:@(lp[i])];
    TTIOCodecContext *ctx = [TTIOCodecContext emptyContext];
    ctx.readLengths = rl;
    id<TTIOCodec> codec = [TTIOCodecRegistry codecForId:TTIOCompressionSeqCm];
    PASS(codec != nil && [codec isContextAware], "M103 seq_cm: registered, context-aware");
    TTIOEncodedChannel *enc = [codec encode:[[TTIODecodedBytes alloc] initWithData:seq]
                                    context:ctx error:&err];
    NSData *encBytes = ((TTIOEncodedDatasetBytes *)enc).bytes;
    PASS([encBytes isEqualToData:blob], "M103 seq_cm: registry encode equals the wrapper's");
    TTIODecodedChannel *dec = [codec decode:[[TTIOBytesPayload alloc] initWithBytes:encBytes]
                                    context:ctx error:&err];
    PASS([((TTIODecodedBytes *)dec).data isEqualToData:seq], "M103 seq_cm: registry decode");
}
