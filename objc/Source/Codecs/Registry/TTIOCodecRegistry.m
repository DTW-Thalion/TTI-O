#import "Codecs/Registry/TTIOCodecRegistry.h"
#import "Codecs/TTIORans.h"
#import "Codecs/TTIOBasePack.h"
#import "Codecs/TTIOQuality.h"
#import "Codecs/TTIODeltaRans.h"
#import "Codecs/Registry/TTIOChannelPayload.h"
#import "Codecs/Registry/TTIODecodedChannel.h"
#import "Codecs/Registry/TTIOEncodedChannel.h"
#import "Codecs/Registry/TTIOCodecContext.h"
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Codecs/TTIOMateInfoV2.h"
#import "Codecs/TTIONameTokenizerV2.h"
#import "Codecs/TTIORefDiffV2.h"
#import "Codecs/TTIOSamTags.h"
#import <pthread.h>

static NSError *_TTIOCodecError(NSString *msg) {
    return [NSError errorWithDomain:@"global.thalion.ttio.CodecRegistry" code:-1
        userInfo:@{NSLocalizedDescriptionKey: msg}];
}

static NSDictionary<NSNumber *, id<TTIOCodec>> *gRegistry = nil;
static pthread_once_t gOnce = PTHREAD_ONCE_INIT;

@interface _TTIORansCodec : NSObject <TTIOCodec>
- (instancetype)initWithId:(TTIOCompression)cid order:(int)order;
@end
@implementation _TTIORansCodec {
    TTIOCompression _cid; int _order;
}
- (instancetype)initWithId:(TTIOCompression)cid order:(int)order {
    if ((self = [super init])) { _cid = cid; _order = order; }
    return self;
}
- (TTIOCompression)codecId { return _cid; }
- (BOOL)isContextAware { return NO; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = TTIORansDecode(((TTIOBytesPayload *)p).bytes, e);
    return out ? [[TTIODecodedBytes alloc] initWithData:out] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = TTIORansEncode(((TTIODecodedBytes *)v).data, _order);
    return [[TTIOEncodedDatasetBytes alloc] initWithBytes:out];
}
@end

@interface _TTIOBasePackCodec : NSObject <TTIOCodec> @end
@implementation _TTIOBasePackCodec
- (TTIOCompression)codecId { return TTIOCompressionBasePack; }
- (BOOL)isContextAware { return NO; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = TTIOBasePackDecode(((TTIOBytesPayload *)p).bytes, e);
    return out ? [[TTIODecodedBytes alloc] initWithData:out] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    return [[TTIOEncodedDatasetBytes alloc] initWithBytes:TTIOBasePackEncode(((TTIODecodedBytes *)v).data)];
}
@end

@interface _TTIOQualityCodec : NSObject <TTIOCodec> @end
@implementation _TTIOQualityCodec
- (TTIOCompression)codecId { return TTIOCompressionQualityBinned; }
- (BOOL)isContextAware { return NO; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = TTIOQualityDecode(((TTIOBytesPayload *)p).bytes, e);
    return out ? [[TTIODecodedBytes alloc] initWithData:out] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    return [[TTIOEncodedDatasetBytes alloc] initWithBytes:TTIOQualityEncode(((TTIODecodedBytes *)v).data)];
}
@end

@interface _TTIODeltaRansCodec : NSObject <TTIOCodec> @end
@implementation _TTIODeltaRansCodec
- (TTIOCompression)codecId { return TTIOCompressionDeltaRansOrder0; }
- (BOOL)isContextAware { return NO; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = TTIODeltaRansDecode(((TTIOBytesPayload *)p).bytes, e);
    return out ? [[TTIODecodedBytes alloc] initWithData:out] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    if (ctx.elementSize == nil) {
        if (e) *e = [NSError errorWithDomain:@"global.thalion.ttio.CodecRegistry" code:-1
            userInfo:@{NSLocalizedDescriptionKey: @"DELTA_RANS encode requires context.elementSize"}];
        return nil;
    }
    NSData *out = TTIODeltaRansEncode(((TTIODecodedBytes *)v).data,
        (uint8_t)ctx.elementSize.unsignedCharValue, e);
    return out ? [[TTIOEncodedDatasetBytes alloc] initWithBytes:out] : nil;
}
@end

// ── name_tok (codec 15): NOT context-aware ──────────────────────────
@interface _TTIONameTokenizedCodec : NSObject <TTIOCodec> @end
@implementation _TTIONameTokenizedCodec
- (TTIOCompression)codecId { return TTIOCompressionNameTokenizedV2; }
- (BOOL)isContextAware { return NO; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSArray<NSString *> *names =
        [TTIONameTokenizerV2 decodeData:((TTIOBytesPayload *)p).bytes error:e];
    return names ? [[TTIODecodedStringList alloc] initWithNames:names] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *out = [TTIONameTokenizerV2 encodeNames:((TTIODecodedStringList *)v).names];
    return [[TTIOEncodedDatasetBytes alloc] initWithBytes:out];
}
@end

// ── fqzcomp (codec 12): context-aware, no embedded reference ────────
@interface _TTIOFqzcompCodec : NSObject <TTIOCodec> @end
@implementation _TTIOFqzcompCodec
- (TTIOCompression)codecId { return TTIOCompressionFqzcompNx16Z; }
- (BOOL)isContextAware { return YES; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSDictionary *r = [TTIOFqzcompNx16Z decodeData:((TTIOBytesPayload *)p).bytes
                                       revcompFlags:ctx.revcompFlags
                                  sequencesProvider:ctx.sequencesProvider
                                              error:e];
    if (!r) return nil;
    return [[TTIODecodedBytes alloc] initWithData:r[@"qualities"]];
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    if (ctx.readLengths == nil || ctx.revcompFlags == nil) {
        if (e) *e = _TTIOCodecError(@"FQZCOMP_NX16_Z encode requires context.readLengths and context.revcompFlags");
        return nil;
    }
    NSInteger hint = ctx.qualStrategyHint != nil
        ? [ctx.qualStrategyHint integerValue] : -1;
    NSData *out = [TTIOFqzcompNx16Z encodeWithQualities:((TTIODecodedBytes *)v).data
                                            readLengths:ctx.readLengths
                                           revcompFlags:ctx.revcompFlags
                                              sequences:ctx.sequences
                                           strategyHint:hint
                                                  error:e];
    return out ? [[TTIOEncodedDatasetBytes alloc] initWithBytes:out] : nil;
}
@end

// ── mate_info (codec 13): context-aware, no embedded reference ──────
@interface _TTIOMateInfoCodec : NSObject <TTIOCodec> @end
@implementation _TTIOMateInfoCodec
- (TTIOCompression)codecId { return TTIOCompressionMateInlineV2; }
- (BOOL)isContextAware { return YES; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    if (ctx.ownChromIds == nil || ctx.ownPositions == nil || ctx.nRecords == nil) {
        if (e) *e = _TTIOCodecError(@"MATE_INLINE_V2 decode requires context.ownChromIds, context.ownPositions, and context.nRecords");
        return nil;
    }
    NSData *mateChromIds = nil, *matePositions = nil, *templateLengths = nil;
    BOOL ok = [TTIOMateInfoV2 decodeData:((TTIOBytesPayload *)p).bytes
                             ownChromIds:ctx.ownChromIds
                            ownPositions:ctx.ownPositions
                                nRecords:ctx.nRecords.unsignedIntegerValue
                         outMateChromIds:&mateChromIds
                        outMatePositions:&matePositions
                      outTemplateLengths:&templateLengths
                                   error:e];
    if (!ok) return nil;
    return [[TTIODecodedMateInfo alloc] initWithMateChromIds:mateChromIds
                                               matePositions:matePositions
                                             templateLengths:templateLengths];
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    if (ctx.ownChromIds == nil || ctx.ownPositions == nil) {
        if (e) *e = _TTIOCodecError(@"MATE_INLINE_V2 encode requires context.ownChromIds and context.ownPositions");
        return nil;
    }
    TTIODecodedMateInfo *mi = (TTIODecodedMateInfo *)v;
    NSData *out = [TTIOMateInfoV2 encodeMateChromIds:mi.mateChromIds
                                       matePositions:mi.matePositions
                                     templateLengths:mi.templateLengths
                                         ownChromIds:ctx.ownChromIds
                                        ownPositions:ctx.ownPositions
                                               error:e];
    return out ? [[TTIOEncodedDatasetBytes alloc] initWithBytes:out] : nil;
}
@end

// ── ref_diff (codec 14): context-aware, needs embedded reference ────
@interface _TTIORefDiffCodec : NSObject <TTIOCodec> @end
@implementation _TTIORefDiffCodec
- (TTIOCompression)codecId { return TTIOCompressionRefDiffV2; }
- (BOOL)isContextAware { return YES; }
- (BOOL)needsEmbeddedReference { return YES; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    if (![p isKindOfClass:[TTIOGroupPayload class]]) {
        if (e) *e = _TTIOCodecError(@"REF_DIFF_V2 decode requires a TTIOGroupPayload");
        return nil;
    }
    id<TTIOStorageGroup> sig = ((TTIOGroupPayload *)p).group;

    // Open the sequences GROUP and read its refdiff_v2 dataset.
    id<TTIOStorageGroup> seqGrp = [sig openGroupNamed:@"sequences" error:e];
    if (!seqGrp) return nil;
    id<TTIOStorageDataset> ds = [seqGrp openDatasetNamed:@"refdiff_v2" error:e];
    if (!ds) return nil;
    id raw = [ds readAll:e];
    if (![raw isKindOfClass:[NSData class]]) return nil;
    NSData *blob = (NSData *)raw;

    // Parse the blob header to extract reference_uri and reference_md5.
    // Header layout: [0:4] "RDF2" magic, [20:36] md5, [36:38] uri_len LE,
    // [38:38+uri_len] uri UTF-8.
    const uint8_t *blobBytes = (const uint8_t *)blob.bytes;
    NSUInteger blobLen = blob.length;
    if (blobLen < 38) {
        if (e) *e = _TTIOCodecError(@"refdiff_v2 blob too short to parse header");
        return nil;
    }
    if (memcmp(blobBytes, "RDF2", 4) != 0) {
        if (e) *e = _TTIOCodecError(@"refdiff_v2 blob magic mismatch (expected 'RDF2')");
        return nil;
    }
    NSData *blobMD5 = [NSData dataWithBytes:blobBytes + 20 length:16];
    uint16_t uriLen = 0;
    memcpy(&uriLen, blobBytes + 36, 2);  // LE uint16
    if (blobLen < (NSUInteger)(38 + uriLen)) {
        if (e) *e = _TTIOCodecError(@"refdiff_v2 blob truncated (uri)");
        return nil;
    }
    NSString *blobURI = [[NSString alloc] initWithBytes:blobBytes + 38
                                                 length:uriLen
                                               encoding:NSUTF8StringEncoding]
                        ?: @"";

    // Single-chromosome constraint (same as v1).
    if (ctx.referenceResolver == nil) {
        if (e) *e = _TTIOCodecError(@"REF_DIFF_V2 decode requires context.referenceResolver");
        return nil;
    }
    NSMutableSet<NSString *> *unique = [NSMutableSet set];
    for (NSString *c in (ctx.chromosomes ?: @[])) {
        if (c.length) [unique addObject:c];
    }
    if (unique.count != 1) {
        if (e) *e = _TTIOCodecError([NSString stringWithFormat:
            @"refdiff_v2 decode: expected single-chromosome run, got %lu chromosomes",
            (unsigned long)unique.count]);
        return nil;
    }
    NSString *chrom = [unique anyObject];
    NSData *ref = [ctx.referenceResolver resolveURI:blobURI
                                        expectedMD5:blobMD5
                                         chromosome:chrom
                                              error:e];
    if (!ref) return nil;

    NSArray<NSString *> *cigars = ctx.cigarsProvider ? ctx.cigarsProvider() : @[];

    NSData *outSeq = nil, *outOff = nil;
    BOOL ok = [TTIORefDiffV2 decodeData:blob
                             positions:ctx.positions
                          cigarStrings:cigars
                             reference:ref
                                nReads:ctx.readCount.unsignedIntegerValue
                            totalBases:ctx.totalBases.unsignedIntegerValue
                          outSequences:&outSeq
                            outOffsets:&outOff
                                 error:e];
    if (!ok) return nil;
    return [[TTIODecodedBytes alloc] initWithData:outSeq];
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSData *blob =
        [TTIORefDiffV2 encodeSequences:((TTIODecodedBytes *)v).data
                               offsets:ctx.offsets
                             positions:ctx.positions
                          cigarStrings:(ctx.cigarsProvider ? ctx.cigarsProvider() : @[])
                             reference:ctx.reference
                          referenceMd5:ctx.referenceMd5
                          referenceUri:ctx.referenceUri
                         readsPerSlice:(ctx.readsPerSlice ? ctx.readsPerSlice.unsignedIntegerValue : 10000)
                            sliceBytes:(ctx.sliceBytes ? ctx.sliceBytes.unsignedLongLongValue : 0)
                                 error:e];
    if (!blob) return nil;
    return [[TTIOEncodedGroupLayout alloc] initWithChildren:@{@"refdiff_v2": blob}
                                                      attrs:@{}];
}
@end

// ── sam_tags (codec 18, M101): context-aware, no embedded reference ─
// MD/NM derivation needs the reads' sequences, CIGARs, positions and
// chromosome ids plus the reference bases; without tag references
// none of that is read.
@interface _TTIOSamTagsCodec : NSObject <TTIOCodec> @end
@implementation _TTIOSamTagsCodec
- (TTIOCompression)codecId { return TTIOCompressionSamTags; }
- (BOOL)isContextAware { return YES; }
- (BOOL)needsEmbeddedReference { return NO; }
- (TTIODecodedChannel *)decode:(TTIOChannelPayload *)p context:(TTIOCodecContext *)ctx error:(NSError **)e {
    NSUInteger n = ctx.readCount ? ctx.readCount.unsignedIntegerValue : 0;
    TTIOSamTagsContext *tc = [[TTIOSamTagsContext alloc] init];
    tc.references = ctx.tagReferencesProvider ? ctx.tagReferencesProvider() : nil;
    if ([tc derives]) {
        NSMutableData *off = [NSMutableData dataWithLength:(n + 1) * sizeof(uint64_t)];
        uint64_t *o = (uint64_t *)[off mutableBytes];
        for (NSUInteger i = 0; i < n && i < ctx.readLengths.count; i++) {
            o[i + 1] = o[i] + (uint64_t)[ctx.readLengths[i] unsignedLongLongValue];
        }
        tc.sequences = ctx.sequencesProvider ? ctx.sequencesProvider() : nil;
        tc.seqOffsets = off;
        tc.cigars = ctx.cigarsProvider ? ctx.cigarsProvider() : nil;
        tc.positions = ctx.positions;
        tc.chromIds = ctx.ownChromIds;
    } else {
        tc = nil;
    }
    NSArray<NSString *> *tags = [TTIOSamTags decodeData:((TTIOBytesPayload *)p).bytes
                                                 nReads:n context:tc error:e];
    return tags ? [[TTIODecodedStringList alloc] initWithNames:tags] : nil;
}
- (TTIOEncodedChannel *)encode:(TTIODecodedChannel *)v context:(TTIOCodecContext *)ctx error:(NSError **)e {
    TTIOSamTagsContext *tc = [[TTIOSamTagsContext alloc] init];
    tc.references = ctx.tagReferences;
    if ([tc derives]) {
        tc.sequences = ctx.sequences;
        tc.seqOffsets = ctx.offsets;
        tc.cigars = ctx.cigarsProvider ? ctx.cigarsProvider() : nil;
        tc.positions = ctx.positions;
        tc.chromIds = ctx.ownChromIds;
    } else {
        tc = nil;
    }
    NSData *out = [TTIOSamTags encodeTags:((TTIODecodedStringList *)v).names
                                  context:tc error:e];
    return out ? [[TTIOEncodedDatasetBytes alloc] initWithBytes:out] : nil;
}
@end

static void _buildRegistry(void) {
    gRegistry = @{
        @(TTIOCompressionRansOrder0): [[_TTIORansCodec alloc] initWithId:TTIOCompressionRansOrder0 order:0],
        @(TTIOCompressionRansOrder1): [[_TTIORansCodec alloc] initWithId:TTIOCompressionRansOrder1 order:1],
        @(TTIOCompressionBasePack):   [[_TTIOBasePackCodec alloc] init],
        @(TTIOCompressionQualityBinned): [[_TTIOQualityCodec alloc] init],
        @(TTIOCompressionDeltaRansOrder0): [[_TTIODeltaRansCodec alloc] init],
        @(TTIOCompressionNameTokenizedV2): [[_TTIONameTokenizedCodec alloc] init],
        @(TTIOCompressionFqzcompNx16Z):    [[_TTIOFqzcompCodec alloc] init],
        @(TTIOCompressionMateInlineV2):    [[_TTIOMateInfoCodec alloc] init],
        @(TTIOCompressionRefDiffV2):       [[_TTIORefDiffCodec alloc] init],
        @(TTIOCompressionSamTags):         [[_TTIOSamTagsCodec alloc] init],
    };
}

@implementation TTIOCodecRegistry
+ (nullable id<TTIOCodec>)codecForId:(TTIOCompression)codecId {
    pthread_once(&gOnce, _buildRegistry);
    return gRegistry[@(codecId)];
}
@end
