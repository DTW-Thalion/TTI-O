/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOGenomicWriteContext.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"
#import "Providers/TTIOMemoryProvider.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "ValueClasses/TTIOEnums.h"
#import "Codecs/TTIOSeqCm.h"
#import <stdatomic.h>

@implementation TTIOBlockBlobs

- (instancetype)initWithBlobs:(NSDictionary<NSString *, NSData *> *)blobs
                       codecs:(NSDictionary<NSString *, NSNumber *> *)codecs
                   extraAttrs:(NSDictionary<NSString *, NSDictionary<NSString *, id> *> *)extraAttrs
                       nReads:(NSUInteger)nReads
                       nBases:(uint64_t)nBases
{
    self = [super init];
    if (self) {
        _blobs = [blobs copy];
        _codecs = [codecs copy];
        _extraAttrs = [extraAttrs copy];
        _nReads = nReads;
        _nBases = nBases;
    }
    return self;
}

@end

static NSData *ttioSubData(NSData *d, NSUInteger elem, NSUInteger start, NSUInteger stop)
{
    NSUInteger from = MIN(start * elem, d.length);
    NSUInteger to = MIN(stop * elem, d.length);
    if (to < from) to = from;
    return [d subdataWithRange:NSMakeRange(from, to - from)];
}

static id<TTIOStorageGroup> ttioTryGroup(id<TTIOStorageGroup> parent, NSString *name)
{
    if (![parent hasChildNamed:name]) return nil;
    return [parent openGroupNamed:name error:NULL];
}

@implementation TTIOGenomicBlocks

+ (NSArray<NSString *> *)requiredBlockChannels
{
    return @[@"sequences", @"qualities", @"read_names", @"cigars", @"mate_info"];
}

+ (NSArray<NSString *> *)optionalBlockChannels
{
    return @[@"tags"];
}

+ (NSArray<NSString *> *)blockChannels
{
    return [[self requiredBlockChannels] arrayByAddingObjectsFromArray:[self optionalBlockChannels]];
}

+ (TTIOWrittenGenomicRun *)sliceRun:(TTIOWrittenGenomicRun *)run
                               from:(NSUInteger)start
                                 to:(NSUInteger)stop
{
    const uint64_t *offs = (const uint64_t *)run.offsetsData.bytes;
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    NSUInteger n = stop > start ? stop - start : 0;
    uint64_t b0 = 0, b1 = 0;
    if (n > 0) {
        b0 = offs[start];
        b1 = offs[stop - 1] + lens[stop - 1];
    }
    NSMutableData *offsets = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    uint64_t *o = (uint64_t *)offsets.mutableBytes;
    for (NSUInteger i = 0; i < n; i++) o[i] = offs[start + i] - b0;
    NSRange rr = NSMakeRange(start, n);
    TTIOWrittenGenomicRun *s = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:run.acquisitionMode
                   referenceUri:run.referenceUri
                       platform:run.platform
                     sampleName:run.sampleName
                      positions:ttioSubData(run.positionsData, sizeof(int64_t), start, stop)
               mappingQualities:ttioSubData(run.mappingQualitiesData, 1, start, stop)
                          flags:ttioSubData(run.flagsData, sizeof(uint32_t), start, stop)
                      sequences:ttioSubData(run.sequencesData, 1, (NSUInteger)b0, (NSUInteger)b1)
                      qualities:ttioSubData(run.qualitiesData, 1, (NSUInteger)b0, (NSUInteger)b1)
                        offsets:offsets
                        lengths:ttioSubData(run.lengthsData, sizeof(uint32_t), start, stop)
                         cigars:[run.cigars subarrayWithRange:rr]
                      readNames:[run.readNames subarrayWithRange:rr]
                mateChromosomes:[run.mateChromosomes subarrayWithRange:rr]
                  matePositions:ttioSubData(run.matePositionsData, sizeof(int64_t), start, stop)
                templateLengths:ttioSubData(run.templateLengthsData, sizeof(int32_t), start, stop)
                    chromosomes:[run.chromosomes subarrayWithRange:rr]
              signalCompression:run.signalCompression
           signalCodecOverrides:run.signalCodecOverrides];
    s.optDisableQualitiesV5 = run.optDisableQualitiesV5;
    s.embedReference = run.embedReference;
    s.referenceChromSeqs = run.referenceChromSeqs;
    s.externalReferencePath = run.externalReferencePath;
    s.optLegacyWholeChannel = run.optLegacyWholeChannel;
    s.readRole = run.readRole;
    s.refDiffSliceBytes = run.refDiffSliceBytes;
    s.tags = run.tags ? [run.tags subarrayWithRange:rr] : nil;
    return s;
}

+ (TTIOWrittenGenomicRun *)concatRuns:(NSArray<TTIOWrittenGenomicRun *> *)parts
{
    if (parts.count == 1) return parts[0];
    TTIOWrittenGenomicRun *first = parts[0];
    NSMutableData *positions = [NSMutableData data], *mapqs = [NSMutableData data],
                  *flags = [NSMutableData data], *lengths = [NSMutableData data],
                  *matePos = [NSMutableData data], *tlens = [NSMutableData data],
                  *seqs = [NSMutableData data], *quals = [NSMutableData data];
    NSMutableArray *cigars = [NSMutableArray array], *names = [NSMutableArray array],
                   *mateChroms = [NSMutableArray array], *chroms = [NSMutableArray array];
    // M101: nil when no part carries tags; a tag-less part among tagged
    // ones contributes "" per read.
    NSMutableArray *tags = nil;
    for (TTIOWrittenGenomicRun *p in parts) {
        if (p.tags != nil) { tags = [NSMutableArray array]; break; }
    }
    for (TTIOWrittenGenomicRun *p in parts) {
        if (tags != nil) {
            if (p.tags != nil) {
                [tags addObjectsFromArray:p.tags];
            } else {
                for (NSUInteger i = 0; i < p.cigars.count; i++) [tags addObject:@""];
            }
        }
        [positions appendData:p.positionsData];
        [mapqs appendData:p.mappingQualitiesData];
        [flags appendData:p.flagsData];
        [lengths appendData:p.lengthsData];
        [matePos appendData:p.matePositionsData];
        [tlens appendData:p.templateLengthsData];
        [seqs appendData:p.sequencesData];
        [quals appendData:p.qualitiesData];
        [cigars addObjectsFromArray:p.cigars];
        [names addObjectsFromArray:p.readNames];
        [mateChroms addObjectsFromArray:p.mateChromosomes];
        [chroms addObjectsFromArray:p.chromosomes];
    }
    NSUInteger n = lengths.length / sizeof(uint32_t);
    const uint32_t *l = (const uint32_t *)lengths.bytes;
    NSMutableData *offsets = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    uint64_t *o = (uint64_t *)offsets.mutableBytes;
    uint64_t acc = 0;
    for (NSUInteger i = 0; i < n; i++) { o[i] = acc; acc += l[i]; }
    TTIOWrittenGenomicRun *c = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:first.acquisitionMode
                   referenceUri:first.referenceUri
                       platform:first.platform
                     sampleName:first.sampleName
                      positions:positions
               mappingQualities:mapqs
                          flags:flags
                      sequences:seqs
                      qualities:quals
                        offsets:offsets
                        lengths:lengths
                         cigars:cigars
                      readNames:names
                mateChromosomes:mateChroms
                  matePositions:matePos
                templateLengths:tlens
                    chromosomes:chroms
              signalCompression:first.signalCompression
           signalCodecOverrides:first.signalCodecOverrides];
    c.optDisableQualitiesV5 = first.optDisableQualitiesV5;
    c.embedReference = first.embedReference;
    c.referenceChromSeqs = first.referenceChromSeqs;
    c.externalReferencePath = first.externalReferencePath;
    c.optLegacyWholeChannel = first.optLegacyWholeChannel;
    c.readRole = first.readRole;
    c.refDiffSliceBytes = first.refDiffSliceBytes;
    c.tags = tags;
    return c;
}

+ (TTIOBlockBlobs *)encodeBlock:(TTIOWrittenGenomicRun *)block
                        context:(TTIOGenomicWriteContext *)ctx
                          error:(NSError **)error
{
    NSMutableDictionary *ov = [NSMutableDictionary dictionaryWithDictionary:block.signalCodecOverrides];
    if (ov[@"cigars"] == nil) ov[@"cigars"] = @(TTIOCompressionRansOrder0);
    if (ov[@"qualities"] == nil) {
        BOOL zero = NO;
        const uint32_t *lens = (const uint32_t *)block.lengthsData.bytes;
        NSUInteger n = block.lengthsData.length / sizeof(uint32_t);
        for (NSUInteger i = 0; i < n; i++) if (lens[i] == 0) { zero = YES; break; }
        ov[@"qualities"] = zero ? @(TTIOCompressionRansOrder0) : @(TTIOCompressionFqzcompNx16Z);
    }
    // Sequences without a reference: SEQ_CM (M103), or RANS_ORDER1 when
    // libttio_rans is not linked.
    if (ov[@"sequences"] == nil && block.referenceChromSeqs == nil) {
        ov[@"sequences"] = [TTIOSeqCm nativeAvailable] ? @(TTIOCompressionSeqCm)
                                                       : @(TTIOCompressionRansOrder1);
    }
    TTIOWrittenGenomicRun *b = [[block copyWithSignalCodecOverrides:ov] copyWithProvenance:@[]];

    static _Atomic unsigned long counter = 0;
    unsigned long seq = atomic_fetch_add(&counter, 1);
    NSString *url = [NSString stringWithFormat:@"memory://ttio-block-encode-%p-%lu", (void *)block, seq];
    [TTIOMemoryProvider discardStore:url];
    id<TTIOStorageProvider> mem = [[TTIOProviderRegistry sharedRegistry]
        openURL:url mode:TTIOStorageOpenModeCreate provider:@"memory" error:error];
    if (!mem) return nil;
    TTIOBlockBlobs *result = nil;
    @try {
        id<TTIOStorageGroup> root = [mem rootGroupWithError:error];
        if (!root) return nil;
        if (![TTIOSpectralDataset writeGenomicRunStorage:b toGroup:root name:@"b"
                                                 context:ctx error:error]) return nil;
        id<TTIOStorageGroup> sc = [[root openGroupNamed:@"b" error:error]
                                   openGroupNamed:@"signal_channels" error:error];
        if (!sc) return nil;
        NSMutableDictionary *blobs = [NSMutableDictionary dictionary];
        NSMutableDictionary *codecs = [NSMutableDictionary dictionary];
        NSMutableDictionary *extra = [NSMutableDictionary dictionary];
        for (NSString *ch in [self blockChannels]) {
            id<TTIOStorageDataset> ds = nil;
            if ([ch isEqualToString:@"sequences"]) {
                id<TTIOStorageGroup> g = ttioTryGroup(sc, @"sequences");
                if (g != nil && [g hasChildNamed:@"refdiff_v2"]) {
                    ds = [g openDatasetNamed:@"refdiff_v2" error:NULL];
                } else if (g == nil && [sc hasChildNamed:@"sequences"]) {
                    ds = [sc openDatasetNamed:@"sequences" error:NULL];
                }
            } else if ([ch isEqualToString:@"mate_info"]) {
                id<TTIOStorageGroup> g = ttioTryGroup(sc, @"mate_info");
                if (g != nil && [g hasChildNamed:@"inline_v2"]) {
                    ds = [g openDatasetNamed:@"inline_v2" error:NULL];
                }
            } else if ([sc hasChildNamed:ch] && ttioTryGroup(sc, ch) == nil) {
                ds = [sc openDatasetNamed:ch error:NULL];
            }
            if (ds == nil) {
                blobs[ch] = [NSData data];
                codecs[ch] = @0;
                extra[ch] = @{};
                continue;
            }
            id raw = [ds readAll:NULL];
            blobs[ch] = [raw isKindOfClass:[NSData class]] ? raw : [NSData data];
            NSNumber *codec = @0;
            NSMutableDictionary *attrs = [NSMutableDictionary dictionary];
            for (NSString *k in [ds attributeNames]) {
                id v = [ds attributeValueForName:k error:NULL];
                if ([k isEqualToString:@"compression"]) codec = @([v unsignedIntegerValue]);
                else if (v != nil) attrs[k] = v;
            }
            codecs[ch] = codec;
            extra[ch] = attrs;
        }
        uint64_t nBases = 0;
        const uint32_t *lens = (const uint32_t *)block.lengthsData.bytes;
        NSUInteger n = block.lengthsData.length / sizeof(uint32_t);
        for (NSUInteger i = 0; i < n; i++) nBases += lens[i];
        result = [[TTIOBlockBlobs alloc] initWithBlobs:blobs codecs:codecs extraAttrs:extra
                                                nReads:block.readCount nBases:nBases];
    } @finally {
        [mem close];
        [TTIOMemoryProvider discardStore:url];
    }
    return result;
}

@end
