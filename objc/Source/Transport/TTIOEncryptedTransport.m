/*
 * TTIOEncryptedTransport.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOEncryptedTransport
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Transport/TTIOEncryptedTransport.h
 *
 * Transport layer for per-AU encrypted .tio files. Pushes
 * ciphertext from <channel>_segments compounds onto the wire
 * unmodified (the server never decrypts in transit) and rebuilds an
 * encrypted .tio on the receiver side with the same wrapped DEK.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "TTIOEncryptedTransport.h"
#import "TTIOEncryptedTransport+Conformance.h"
#import "TTIOTransportPacket.h"
#import "TTIOTransportReader.h"
#import "TTIOAccessUnit.h"
#import "Core/TTIOPortability.h"
#import "Protection/TTIOPerAUEncryption.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOCompoundField.h"
#import "Dataset/TTIOCompoundIO.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "ValueClasses/TTIOEnums.h"
#import "Genomics/TTIOBlockTable.h"          // M99.1 blocks_v1 sidecars
#import "Genomics/TTIOBlockView.h"           // readNamesIn:named:
#import "Genomics/TTIOGenomicBlocks.h"       // blockChannels
#import "Genomics/TTIOGenomicStreamWriter.h" // indexFields

#include <string.h>

static NSString *const kDomain = @"TTIOEncryptedTransportErrorDomain";
static NSError *makeErr(NSInteger c, NSString *fmt, ...) NS_FORMAT_FUNCTION(2, 3);
static NSError *makeErr(NSInteger c, NSString *fmt, ...)
{
    va_list args; va_start(args, fmt);
    NSString *m = [[NSString alloc] initWithFormat:fmt arguments:args];
    va_end(args);
    return [NSError errorWithDomain:kDomain code:c
                            userInfo:@{NSLocalizedDescriptionKey: m}];
}

// ---------------------------------------------------------------- LE helpers

static void appendU16LE(NSMutableData *buf, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFFu), (uint8_t)((v >> 8) & 0xFFu)};
    [buf appendBytes:b length:2];
}
static void appendU32LE(NSMutableData *buf, uint32_t v) {
    uint8_t b[4];
    b[0] = (uint8_t)(v & 0xFFu);
    b[1] = (uint8_t)((v >> 8) & 0xFFu);
    b[2] = (uint8_t)((v >> 16) & 0xFFu);
    b[3] = (uint8_t)((v >> 24) & 0xFFu);
    [buf appendBytes:b length:4];
}
static void appendLEString(NSMutableData *buf, NSString *s, int width) {
    NSData *d = [(s ?: @"") dataUsingEncoding:NSUTF8StringEncoding];
    if (width == 2) appendU16LE(buf, (uint16_t)d.length);
    else            appendU32LE(buf, (uint32_t)d.length);
    [buf appendData:d];
}
static void appendU64LE(NSMutableData *buf, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (uint8_t)((v >> (8 * i)) & 0xFFu);
    [buf appendBytes:b length:8];
}
static uint64_t readU64LE(const uint8_t *b) {
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= ((uint64_t)b[i]) << (8 * i);
    return v;
}
static uint16_t readU16LE(const uint8_t *b) {
    return (uint16_t)((uint32_t)b[0] | ((uint32_t)b[1] << 8));
}
static uint32_t readU32LE(const uint8_t *b) {
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8)
         | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

// FD-1 Phase A append-only trailing block: count(u16) + per recipient
// {recipientId, kekAlgorithm, wrappedDek}. Byte-for-byte identical to the
// Python `_encode_recipient_block` / Java `encodeRecipientBlock` (the
// cross-language contract). nil/empty input -> empty data, so
// single-recipient packets stay byte-identical to transport-spec §4.4.
static NSData *encodeRecipientBlock(NSArray<NSDictionary *> *additional) {
    if (additional.count == 0) return [NSData data];
    NSMutableData *buf = [NSMutableData data];
    appendU16LE(buf, (uint16_t)additional.count);
    for (NSDictionary *r in additional) {
        appendLEString(buf, r[@"recipientId"] ?: @"", 2);
        appendLEString(buf, r[@"kekAlgorithm"] ?: @"", 2);
        NSData *wd = r[@"wrappedDek"] ?: [NSData data];
        appendU32LE(buf, (uint32_t)wd.length);
        [buf appendData:wd];
    }
    return buf;
}

// Decode the trailing block from a cursor; advances *offp. Returns []
// when fewer than 2 bytes remain.
static NSArray<NSDictionary *> *decodeRecipientBlock(const uint8_t *b,
                                                     NSUInteger len,
                                                     NSUInteger *offp) {
    NSMutableArray<NSDictionary *> *out = [NSMutableArray array];
    NSUInteger off = *offp;
    if (off + 2 > len) { *offp = off; return out; }
    uint16_t n = readU16LE(&b[off]); off += 2;
    for (uint16_t i = 0; i < n; i++) {
        if (off + 2 > len) break;
        uint16_t ridLen = readU16LE(&b[off]); off += 2;
        if (off + ridLen > len) break;
        NSString *rid = [[NSString alloc] initWithBytes:&b[off] length:ridLen
                                               encoding:NSUTF8StringEncoding];
        off += ridLen;
        if (off + 2 > len) break;
        uint16_t kaLen = readU16LE(&b[off]); off += 2;
        if (off + kaLen > len) break;
        NSString *ka = [[NSString alloc] initWithBytes:&b[off] length:kaLen
                                              encoding:NSUTF8StringEncoding];
        off += kaLen;
        if (off + 4 > len) break;
        uint32_t wl = readU32LE(&b[off]); off += 4;
        if (off + wl > len) break;
        NSData *wd = [NSData dataWithBytes:&b[off] length:wl]; off += wl;
        [out addObject:@{@"recipientId": rid, @"kekAlgorithm": ka,
                         @"wrappedDek": wd}];
    }
    *offp = off;
    return out;
}

static NSString *readStringAttr(id<TTIOStorageGroup> g, NSString *name)
{
    if (!g || ![g hasAttributeNamed:name]) return nil;
    id v = [g attributeValueForName:name error:NULL];
    if ([v isKindOfClass:[NSString class]]) return (NSString *)v;
    if ([v isKindOfClass:[NSData class]]) {
        return [[NSString alloc] initWithData:(NSData *)v
                                      encoding:NSUTF8StringEncoding];
    }
    return nil;
}

// FD-1 C-2a: append the optional trailing section (the stored recipient
// block + an optional server_kek_id) after the five transport-spec §4.4
// fields. Emitted only when a recipient block OR a server_kek_id is
// present, so single-recipient BYOK packets stay byte-identical. A
// single-recipient server-processable run emits
// additional_recipient_count = 0 followed by server_kek_id.
static void appendProtectionTrailing(NSMutableData *pm,
                                     id<TTIOStorageGroup> sig,
                                     NSString *firstChannel)
{
    NSString *recipAttr = [NSString stringWithFormat:
        @"%@_wrapped_dek_recipients", firstChannel];
    NSData *recipBlock = nil;
    if ([sig hasAttributeNamed:recipAttr]) {
        id rv = [sig attributeValueForName:recipAttr error:NULL];
        if ([rv isKindOfClass:[NSData class]]) recipBlock = (NSData *)rv;
    }
    NSString *serverKekId = readStringAttr(sig,
        [NSString stringWithFormat:@"%@_server_kek_id", firstChannel]);
    if (recipBlock.length > 0) {
        [pm appendData:recipBlock];
    } else if (serverKekId.length > 0) {
        appendU16LE(pm, 0);   // additional_recipient_count = 0
    }
    if (serverKekId.length > 0) {
        appendLEString(pm, serverKekId, 2);
    }
}
static NSArray<NSString *> *splitChannelNames(NSString *raw) {
    if (!raw.length) return @[];
    NSArray *parts = [raw componentsSeparatedByString:@","];
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:parts.count];
    for (NSString *p in parts) if (p.length) [out addObject:p];
    return out;
}

static NSArray<TTIOCompoundField *> *channelSegFields(void) {
    return @[
        [TTIOCompoundField fieldWithName:@"offset" kind:TTIOCompoundFieldKindInt64],
        [TTIOCompoundField fieldWithName:@"length" kind:TTIOCompoundFieldKindUInt32],
        [TTIOCompoundField fieldWithName:@"iv" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"tag" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"ciphertext" kind:TTIOCompoundFieldKindVLBytes],
    ];
}
static NSArray<TTIOCompoundField *> *headerSegFields(void) {
    return @[
        [TTIOCompoundField fieldWithName:@"iv" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"tag" kind:TTIOCompoundFieldKindVLBytes],
        [TTIOCompoundField fieldWithName:@"ciphertext" kind:TTIOCompoundFieldKindVLBytes],
    ];
}

static NSArray<NSString *> *readFeatures(id<TTIOStorageGroup> root) {
    if (![root hasAttributeNamed:@"ttio_features"]) return @[];
    NSString *s = readStringAttr(root, @"ttio_features");
    if (!s.length) return @[];
    id parsed = [NSJSONSerialization JSONObjectWithData:
                    [s dataUsingEncoding:NSUTF8StringEncoding]
                                                  options:0 error:NULL];
    return [parsed isKindOfClass:[NSArray class]] ? parsed : @[];
}

// ---------------------------------------------------------------- packet

static uint64_t nowNs(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static NSData *encodeHeader(TTIOTransportPacketType type, uint16_t flags,
                              uint16_t datasetId, uint32_t auSeq, uint32_t plen)
{
    TTIOTransportPacketHeader *h =
        [[TTIOTransportPacketHeader alloc] initWithPacketType:type
                                                          flags:flags
                                                      datasetId:datasetId
                                                     auSequence:auSeq
                                                  payloadLength:plen
                                                    timestampNs:nowNs()];
    return [h encode];
}

// ---------------------------------------------------------------- writer

// We need a writer that lets us set PacketFlagEncrypted (and
// PacketFlagEncryptedHeader) on outgoing AU packets. The public
// TTIOTransportWriter doesn't take a flags argument in
// -writeAccessUnit:. We work around this by writing directly to an
// NSMutableData we manage, encoding the header bytes with the
// right flags, and appending payload bytes.
//
// Rationale: the same approach Python uses in transport/encrypted.py
// where it reaches into TransportWriter._emit / _stream — TTIO's
// ObjC writer's stream sink is NSMutableData which we can append to
// directly via writeBytes on the writer's associated data buffer.
// Since the writer's internal sink isn't public, we instead build
// the stream ourselves and feed it to the writer through the public
// -writeStreamHeader.. -writeDatasetHeader..  APIs for header packets
// and a raw append for encrypted AUs via -writeBytesToInternalSink:.
//
// Implementation detail: we added a tiny internal API on
// TTIOTransportWriter, -_writeRawBytes:, to let this module emit AUs
// with arbitrary flag bits. That helper is declared in a category
// below (private to this file) so the writer's public header stays
// unchanged.

@interface TTIOTransportWriter (EncryptedTransport)
- (void)_writeRawPacketHeader:(TTIOTransportPacketType)type
                         flags:(uint16_t)flags
                     datasetId:(uint16_t)datasetId
                    auSequence:(uint32_t)auSeq
                       payload:(NSData *)payload;
@end

@implementation TTIOTransportWriter (EncryptedTransport)
// Uses the private ivars of TTIOTransportWriter via KVC / valueForKey
// since we can't import its implementation file. The public -_emit...
// methods we need are actually declared in the extension block of
// TTIOTransportWriter.m but not exposed publicly. Rather than modify
// the public header, we call the private method via performSelector
// with explicit argument marshalling.
- (void)_writeRawPacketHeader:(TTIOTransportPacketType)type
                         flags:(uint16_t)flags
                     datasetId:(uint16_t)datasetId
                    auSequence:(uint32_t)auSeq
                       payload:(NSData *)payload
{
    // Synthesise a 24-byte header with arbitrary flags, then inject the
    // raw packet (header + payload) through the writer's sink. The public
    // writer API has no flags argument, so we reach the private `_sink`
    // ivar via KVC and use its public -writeData:. This works for every
    // sink type (in-memory, file, custom). The previous code fished out
    // `fileHandle` / `dataBuffer` ivars, which were removed when the
    // writer moved to a sink model -- silently breaking every encrypted
    // write (the gnustep test runner didn't gate failures, so CI stayed
    // green).
    NSData *headerBytes = encodeHeader(type, flags, datasetId, auSeq,
                                          (uint32_t)payload.length);
    id<TTIOTransportWriterSink> sink = [self valueForKey:@"sink"];
    [sink writeData:headerBytes];
    [sink writeData:payload];
}
@end


// ---------------------------------------------------------------- reader helpers

@interface DatasetAccumulator : NSObject
@property (nonatomic, copy) NSString *name;
@property (nonatomic) uint8_t acquisitionMode;
@property (nonatomic, copy) NSString *spectrumClass;
@property (nonatomic, strong) NSMutableArray<NSString *> *channelNames;
@property (nonatomic, strong) NSMutableDictionary<NSString *, NSMutableArray *> *channelSegments;
@property (nonatomic, strong) NSMutableArray<TTIOHeaderSegment *> *headerSegments;
@property (nonatomic) BOOL usedEncryptedHeaders;
// genomic-only fields. isGenomic is set when DATASET_HEADER
// declares spectrum_class==TTIOGenomicRead. Per-AU genomic suffix
// values accumulate into the four arrays so the materialiser can
// rebuild the plaintext genomic_index columns. genomicMetadataJSON
// carries the per-run metadata (modality/platform/reference_uri/
// sample_name) from the DATASET_HEADER instrument_json slot.
@property (nonatomic) BOOL isGenomic;
@property (nonatomic, strong) NSMutableArray<NSString *> *genomicChromosomes;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *genomicPositions;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *genomicMapqs;
@property (nonatomic, strong) NSMutableArray<NSNumber *> *genomicFlags;
@property (nonatomic, copy) NSString *genomicMetadataJSON;
// M99.1 blocks_v1 sidecars; runSidecar stays nil for legacy-shaped
// genomic runs. Decoded packet dictionaries (transport-spec §4.24;
// see decodeGenomicRunSidecar / decodeBlockSidecar).
@property (nonatomic, strong) NSDictionary *runSidecar;
@property (nonatomic, strong) NSMutableArray<NSDictionary *> *blockSidecars;
@end
@implementation DatasetAccumulator
@end

@interface ProtectionMeta : NSObject
@property (nonatomic, copy) NSString *cipherSuite;
@property (nonatomic, copy) NSString *kekAlgorithm;
@property (nonatomic, strong) NSData *wrappedDek;
// FD-1 Phase A: additional DEK recipients (same DEK wrapped under other
// keys). Each element is a dict {recipientId:NSString, kekAlgorithm:
// NSString, wrappedDek:NSData}. Empty/nil for single-recipient runs.
@property (nonatomic, strong) NSArray<NSDictionary *> *additionalRecipients;
// FD-1 C-2a: the server kek_id naming the primary recipient's KEK, or nil
// (BYOK / not server-processable).
@property (nonatomic, copy) NSString *serverKekId;
@end
@implementation ProtectionMeta
@end

// ------------------------------------------- M99.1 blocks_v1 sidecars

static int64_t sidecarIntAttr(id<TTIOStorageGroup> g, NSString *name)
{
    if (![g hasAttributeNamed:name]) return 0;
    id v = [g attributeValueForName:name error:NULL];
    return [v respondsToSelector:@selector(longLongValue)]
        ? [v longLongValue] : 0;
}

// One GenomicRunSidecar, then one BlockSidecar per block, for a
// blocks_v1 per-AU-encrypted genomic run (transport-spec §4.24).
static BOOL emitBlocksV1Sidecars(TTIOTransportWriter *writer,
                                 id<TTIOStorageGroup> gRun,
                                 id<TTIOStorageGroup> gSig,
                                 uint16_t did,
                                 NSError **error)
{
    TTIOBlockTable *table = [TTIOBlockTable readFromRunGroup:gRun
                                                       error:error];
    if (!table) return NO;

    NSMutableDictionary *attrs = [NSMutableDictionary dictionary];
    int64_t slice = sidecarIntAttr(gRun, @"ref_diff_slice_bytes");
    if (slice != 0) attrs[@"ref_diff_slice_bytes"] = @(slice);
    if (sidecarIntAttr(gRun, @"opt_disable_qualities_v5") != 0) {
        attrs[@"opt_disable_qualities_v5"] = @1;
    }
    NSString *md5s = readStringAttr(gRun, @"reference_md5s");
    if (md5s.length > 0) attrs[@"reference_md5s"] = md5s;

    NSArray<NSString *> *sidecarBlobChannels =
        @[@"read_names", @"cigars", @"mate_info"];
    NSMutableArray *channels = [NSMutableArray array];
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *blobDs =
        [NSMutableDictionary dictionary];
    for (NSString *ch in sidecarBlobChannels) {
        id<TTIOStorageDataset> ds = nil;
        if ([ch isEqualToString:@"mate_info"]) {
            if ([gSig hasChildNamed:@"mate_info"]) {
                id<TTIOStorageGroup> mg =
                    [gSig openGroupNamed:@"mate_info" error:NULL];
                if (mg && [mg hasChildNamed:@"inline_v2"]) {
                    ds = [mg openDatasetNamed:@"inline_v2" error:NULL];
                }
            }
        } else if ([gSig hasChildNamed:ch]) {
            ds = [gSig openDatasetNamed:ch error:NULL];
        }
        if (ds == nil) continue;
        blobDs[ch] = ds;
        NSMutableDictionary *centry = [NSMutableDictionary dictionary];
        centry[@"name"] = ch;
        centry[@"compression"] = @([ds hasAttributeNamed:@"compression"]
            ? [[ds attributeValueForName:@"compression" error:NULL]
                  longLongValue] : 0);
        NSMutableDictionary *extra = [NSMutableDictionary dictionary];
        for (NSString *k in [ds attributeNames]) {
            if ([k isEqualToString:@"compression"]) continue;
            id v = [ds attributeValueForName:k error:NULL];
            if ([v isKindOfClass:[NSData class]]) {
                v = [[NSString alloc] initWithData:v
                                           encoding:NSUTF8StringEncoding]
                    ?: @"";
            }
            if (v != nil) extra[k] = v;
        }
        if (extra.count > 0) centry[@"extra_attrs"] = extra;
        [channels addObject:centry];
    }

    id<TTIOStorageGroup> gIdx = [gRun openGroupNamed:@"genomic_index"
                                               error:error];
    if (!gIdx) return NO;
    NSArray<NSString *> *chromNames =
        [TTIOBlockView readNamesIn:gIdx named:@"chromosome_names"];
    NSArray<NSString *> *mateNames = @[];
    if ([gSig hasChildNamed:@"mate_info"]) {
        id<TTIOStorageGroup> mg =
            [gSig openGroupNamed:@"mate_info" error:NULL];
        if (mg && [mg hasChildNamed:@"chrom_names"]) {
            mateNames = [TTIOBlockView readNamesIn:mg named:@"chrom_names"];
        }
    }
    // blocks_v1_grouped (M103, transport-spec v0.13): each BlockSidecar
    // ends with the block's rows of genomic_index/input_index.
    NSData *inputIndex = nil;
    if ([readStringAttr(gRun, @"layout") ?: @"" isEqualToString:@"blocks_v1_grouped"]) {
        id<TTIOStorageDataset> iiDs = [gIdx openDatasetNamed:@"input_index" error:error];
        id raw = iiDs ? [iiDs readAll:error] : nil;
        if (![raw isKindOfClass:[NSData class]]
            || [(NSData *)raw length] != (NSUInteger)table.readCount * sizeof(uint32_t)) {
            if (error && *error == nil) *error = makeErr(34,
                @"blocks_v1_grouped run without a usable genomic_index/input_index");
            return NO;
        }
        inputIndex = raw;
    }

    NSData *attrsJson = [NSJSONSerialization
        dataWithJSONObject:attrs options:TTIO_JSON_SORTED_KEYS error:NULL];
    NSData *channelsJson = [NSJSONSerialization
        dataWithJSONObject:channels options:TTIO_JSON_SORTED_KEYS
                     error:NULL];

    NSMutableData *payload = [NSMutableData data];
    appendLEString(payload, readStringAttr(gRun, @"layout") ?: @"", 2);
    appendLEString(payload, readStringAttr(gRun, @"block_policy") ?: @"", 2);
    appendU64LE(payload, (uint64_t)sidecarIntAttr(gRun, @"read_count"));
    appendU64LE(payload, (uint64_t)sidecarIntAttr(gRun, @"base_count"));
    appendU32LE(payload, (uint32_t)attrsJson.length);
    [payload appendData:attrsJson];
    appendU32LE(payload, (uint32_t)channelsJson.length);
    [payload appendData:channelsJson];
    appendU32LE(payload, (uint32_t)chromNames.count);
    for (NSString *n in chromNames) appendLEString(payload, n, 2);
    appendU32LE(payload, (uint32_t)mateNames.count);
    for (NSString *n in mateNames) appendLEString(payload, n, 2);
    [writer _writeRawPacketHeader:TTIOTransportPacketGenomicRunSidecar
                             flags:0
                         datasetId:did
                        auSequence:0
                           payload:payload];

    // The index's own channel set: a run without tags sends the five
    // required entries, exactly as before M101.
    NSArray<NSString *> *blockChannels = table.channels;
    for (NSUInteger b = 0; b < table.count; b++) {
        NSMutableData *bp = [NSMutableData data];
        appendU32LE(bp, (uint32_t)b);
        appendU64LE(bp, [table readStartAt:b]);
        appendU32LE(bp, (uint32_t)[table nReadsAt:b]);
        appendU64LE(bp, [table baseStartAt:b]);
        appendU64LE(bp, [table nBasesAt:b]);
        uint8_t nch = (uint8_t)blockChannels.count;
        [bp appendBytes:&nch length:1];
        for (NSString *ch in blockChannels) {
            appendLEString(bp, ch, 2);
            appendU64LE(bp, [table offsetOf:ch at:b]);
            appendU64LE(bp, [table lengthOf:ch at:b]);
            appendU32LE(bp, table.hasCodecs
                ? (uint32_t)[table codecOf:ch at:b] : 0);
        }
        uint8_t nb = (uint8_t)blobDs.count;
        [bp appendBytes:&nb length:1];
        for (NSString *ch in sidecarBlobChannels) {
            id<TTIOStorageDataset> ds = blobDs[ch];
            if (ds == nil) continue;
            unsigned long long off = [table offsetOf:ch at:b];
            unsigned long long ln = [table lengthOf:ch at:b];
            NSData *data = ln > 0
                ? [ds readSliceAtOffset:(NSUInteger)off
                                  count:(NSUInteger)ln
                                  error:error]
                : [NSData data];
            if (!data) return NO;
            appendLEString(bp, ch, 2);
            appendU32LE(bp, (uint32_t)data.length);
            [bp appendData:data];
        }
        if (inputIndex != nil) {
            // v0.13: the block's rows of input_index, n_reads x uint32 LE.
            const uint32_t *ii = (const uint32_t *)inputIndex.bytes;
            NSUInteger r0 = (NSUInteger)[table readStartAt:b];
            NSUInteger nr = (NSUInteger)[table nReadsAt:b];
            for (NSUInteger k = 0; k < nr; k++) appendU32LE(bp, ii[r0 + k]);
        }
        [writer _writeRawPacketHeader:TTIOTransportPacketBlockSidecar
                                 flags:0
                             datasetId:did
                            auSequence:(uint32_t)b
                               payload:bp];
    }
    return YES;
}

static NSDictionary *decodeGenomicRunSidecar(NSData *payload)
{
    const uint8_t *b = (const uint8_t *)payload.bytes;
    NSUInteger len = payload.length;
    NSUInteger off = 0;
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    if (off + 2 > len) return nil;
    uint16_t sl = readU16LE(&b[off]); off += 2;
    out[@"layout"] = [[NSString alloc] initWithBytes:&b[off] length:sl
                                             encoding:NSUTF8StringEncoding]
        ?: @"";
    off += sl;
    if (off + 2 > len) return nil;
    sl = readU16LE(&b[off]); off += 2;
    out[@"block_policy"] =
        [[NSString alloc] initWithBytes:&b[off] length:sl
                                encoding:NSUTF8StringEncoding] ?: @"";
    off += sl;
    if (off + 16 > len) return nil;
    out[@"read_count"] = @(readU64LE(&b[off])); off += 8;
    out[@"base_count"] = @(readU64LE(&b[off])); off += 8;
    if (off + 4 > len) return nil;
    uint32_t jl = readU32LE(&b[off]); off += 4;
    if (off + jl > len) return nil;
    id attrs = [NSJSONSerialization
        JSONObjectWithData:[payload subdataWithRange:NSMakeRange(off, jl)]
                   options:0 error:NULL];
    off += jl;
    out[@"attrs"] = [attrs isKindOfClass:[NSDictionary class]] ? attrs : @{};
    if (off + 4 > len) return nil;
    jl = readU32LE(&b[off]); off += 4;
    if (off + jl > len) return nil;
    id channels = [NSJSONSerialization
        JSONObjectWithData:[payload subdataWithRange:NSMakeRange(off, jl)]
                   options:0 error:NULL];
    off += jl;
    out[@"channels"] = [channels isKindOfClass:[NSArray class]]
        ? channels : @[];
    if (off + 4 > len) return nil;
    uint32_t n = readU32LE(&b[off]); off += 4;
    NSMutableArray *chromNames = [NSMutableArray arrayWithCapacity:n];
    for (uint32_t i = 0; i < n; i++) {
        if (off + 2 > len) return nil;
        sl = readU16LE(&b[off]); off += 2;
        if (off + sl > len) return nil;
        [chromNames addObject:
            [[NSString alloc] initWithBytes:&b[off] length:sl
                                    encoding:NSUTF8StringEncoding] ?: @""];
        off += sl;
    }
    out[@"chromosome_names"] = chromNames;
    if (off + 4 > len) return nil;
    n = readU32LE(&b[off]); off += 4;
    NSMutableArray *mateNames = [NSMutableArray arrayWithCapacity:n];
    for (uint32_t i = 0; i < n; i++) {
        if (off + 2 > len) return nil;
        sl = readU16LE(&b[off]); off += 2;
        if (off + sl > len) return nil;
        [mateNames addObject:
            [[NSString alloc] initWithBytes:&b[off] length:sl
                                    encoding:NSUTF8StringEncoding] ?: @""];
        off += sl;
    }
    out[@"mate_chrom_names"] = mateNames;
    return out;
}

static NSDictionary *decodeBlockSidecar(NSData *payload)
{
    const uint8_t *b = (const uint8_t *)payload.bytes;
    NSUInteger len = payload.length;
    NSUInteger off = 0;
    if (off + 32 > len) return nil;
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    out[@"block_index"] = @(readU32LE(&b[off])); off += 4;
    out[@"read_start"] = @(readU64LE(&b[off])); off += 8;
    out[@"n_reads"] = @(readU32LE(&b[off])); off += 4;
    out[@"base_start"] = @(readU64LE(&b[off])); off += 8;
    out[@"n_bases"] = @(readU64LE(&b[off])); off += 8;
    if (off + 1 > len) return nil;
    uint8_t nch = b[off++];
    NSMutableDictionary *channels = [NSMutableDictionary dictionary];
    for (uint8_t i = 0; i < nch; i++) {
        if (off + 2 > len) return nil;
        uint16_t sl = readU16LE(&b[off]); off += 2;
        if (off + sl + 20 > len) return nil;
        NSString *name =
            [[NSString alloc] initWithBytes:&b[off] length:sl
                                    encoding:NSUTF8StringEncoding] ?: @"";
        off += sl;
        uint64_t cOff = readU64LE(&b[off]); off += 8;
        uint64_t cLen = readU64LE(&b[off]); off += 8;
        uint32_t codec = readU32LE(&b[off]); off += 4;
        channels[name] = @[@(cOff), @(cLen), @(codec)];
    }
    out[@"channels"] = channels;
    if (off + 1 > len) return nil;
    uint8_t nb = b[off++];
    NSMutableDictionary *blobs = [NSMutableDictionary dictionary];
    for (uint8_t i = 0; i < nb; i++) {
        if (off + 2 > len) return nil;
        uint16_t sl = readU16LE(&b[off]); off += 2;
        if (off + sl > len) return nil;
        NSString *name =
            [[NSString alloc] initWithBytes:&b[off] length:sl
                                    encoding:NSUTF8StringEncoding] ?: @"";
        off += sl;
        if (off + 4 > len) return nil;
        uint32_t bl = readU32LE(&b[off]); off += 4;
        if (off + bl > len) return nil;
        blobs[name] = [payload subdataWithRange:NSMakeRange(off, bl)];
        off += bl;
    }
    out[@"blobs"] = blobs;
    // v0.13 (M103): a blocks_v1_grouped run's sidecar ends with the
    // block's input_index rows, n_reads x uint32 LE; anything else
    // trailing is malformed.
    if (off < len) {
        NSUInteger nr = [out[@"n_reads"] unsignedIntegerValue];
        if (len - off != 4 * nr) return nil;
        NSMutableData *ii = [NSMutableData dataWithLength:nr * sizeof(uint32_t)];
        uint32_t *iv = (uint32_t *)ii.mutableBytes;
        for (NSUInteger k = 0; k < nr; k++) iv[k] = readU32LE(&b[off + 4 * k]);
        out[@"input_index"] = ii;
    }
    return out;
}

// ---------------------------------------------------------------- impl

@implementation TTIOEncryptedTransport

+ (BOOL)isPerAUEncryptedAtPath:(NSString *)path
                  providerName:(NSString *)providerName
{
    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeRead
                                                provider:providerName
                                                   error:NULL];
    if (!sp) return NO;
    BOOL result = NO;
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:NULL];
        result = [readFeatures(root) containsObject:@"opt_per_au_encryption"];
    }
    @finally { [sp close]; }
    return result;
}


+ (BOOL)writeEncryptedDataset:(NSString *)ttioPath
                       writer:(TTIOTransportWriter *)writer
                 providerName:(NSString *)providerName
                        error:(NSError **)error
{
    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:ttioPath
                                                    mode:TTIOStorageOpenModeRead
                                                provider:providerName
                                                   error:error];
    if (!sp) return NO;
    BOOL ok = NO;
    @try {
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return NO;
        NSArray *features = readFeatures(root);
        if (![features containsObject:@"opt_per_au_encryption"]) {
            if (error) *error = makeErr(3,
                @"%@ does not carry opt_per_au_encryption", ttioPath);
            return NO;
        }
        BOOL headersEncrypted = [features containsObject:@"opt_encrypted_au_headers"];

        // Escape hatch for compound reads (see TTIOPerAUFile for the
        // same pattern / rationale).
        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(4,
                @"encrypted transport currently requires HDF5 provider");
            return NO;
        }
        TTIOHDF5File *hdf5File = (TTIOHDF5File *)[sp nativeHandle];
        TTIOHDF5Group *hdf5Root = hdf5File.rootGroup;

        id<TTIOStorageGroup> study = [root openGroupNamed:@"study" error:error];
        id<TTIOStorageGroup> msRuns = [study openGroupNamed:@"ms_runs" error:error];
        if (!study || !msRuns) return NO;

        NSString *title = readStringAttr(study, @"title") ?: @"";
        NSString *isa = readStringAttr(study, @"isa_investigation_id") ?: @"";

        NSMutableArray *runNames = [NSMutableArray array];
        for (NSString *n in [msRuns childNames]) {
            if (![n hasPrefix:@"_"]) [runNames addObject:n];
        }

        // also walk study/genomic_runs/. dataset_id continues
        // from MS so AAD reconstruction matches the per-AU encrypt
        // path (M90.1).
        NSMutableArray *genomicRunNames = [NSMutableArray array];
        id<TTIOStorageGroup> genomicRunsGroup = nil;
        if ([study hasChildNamed:@"genomic_runs"]) {
            genomicRunsGroup = [study openGroupNamed:@"genomic_runs" error:error];
            if (genomicRunsGroup) {
                for (NSString *n in [genomicRunsGroup childNames]) {
                    if (![n hasPrefix:@"_"]) [genomicRunNames addObject:n];
                }
            }
        }
        // blocks_v1 runs ride the same per-read AU stream plus the
        // M99.1 sidecar packets (GenomicRunSidecar + one BlockSidecar
        // per block); the StreamHeader announces them with a wire
        // feature token so receivers know the stream needs them.
        NSMutableDictionary<NSString *, NSString *> *genomicLayouts =
            [NSMutableDictionary dictionary];
        for (NSString *n in genomicRunNames) {
            id<TTIOStorageGroup> g =
                [genomicRunsGroup openGroupNamed:n error:NULL];
            id layoutVal = (g && [g hasAttributeNamed:@"layout"])
                ? [g attributeValueForName:@"layout" error:NULL] : nil;
            NSString *layout = nil;
            if ([layoutVal isKindOfClass:[NSData class]]) {
                layout = [[NSString alloc] initWithData:layoutVal
                                                encoding:NSUTF8StringEncoding];
            } else if (layoutVal != nil) {
                layout = [layoutVal description];
            }
            genomicLayouts[n] = layout ?: @"";
        }
        NSMutableArray *streamFeatures =
            [NSMutableArray arrayWithArray:features];
        if ([genomicLayouts.allValues containsObject:@"blocks_v1"]
            || [genomicLayouts.allValues containsObject:@"blocks_v1_grouped"]) {
            [streamFeatures addObject:@"transport_blocks_v1"];
        }
        // transport-spec v0.13 (M103): BlockSidecars carry input_index.
        if ([genomicLayouts.allValues containsObject:@"blocks_v1_grouped"]) {
            [streamFeatures addObject:@"transport_blocks_v1_grouped"];
        }

        // StreamHeader
        if (![writer writeStreamHeaderWithFormatVersion:@"1.2"
                                                   title:title
                                        isaInvestigation:isa
                                                features:streamFeatures
                                               nDatasets:(uint16_t)(runNames.count + genomicRunNames.count)
                                                   error:error]) return NO;

        // Per-run ProtectionMetadata + DatasetHeader
        uint16_t did = 1;
        for (NSString *runName in runNames) {
            id<TTIOStorageGroup> run = [msRuns openGroupNamed:runName error:error];
            id<TTIOStorageGroup> sig = [run openGroupNamed:@"signal_channels" error:error];
            if (!run || !sig) return NO;

            NSString *channelNamesStr = readStringAttr(sig, @"channel_names") ?: @"";
            NSArray<NSString *> *channelNames = splitChannelNames(channelNamesStr);
            NSString *firstChannel = channelNames.firstObject ?: @"intensity";

            NSString *cipherSuite = readStringAttr(sig,
                [NSString stringWithFormat:@"%@_algorithm", firstChannel]) ?: @"aes-256-gcm";
            NSString *kek = readStringAttr(sig,
                [NSString stringWithFormat:@"%@_kek_algorithm", firstChannel]) ?: @"";
            NSData *wrapped = [NSData data];
            NSString *wrappedAttr = [NSString stringWithFormat:@"%@_wrapped_dek", firstChannel];
            if ([sig hasAttributeNamed:wrappedAttr]) {
                id v = [sig attributeValueForName:wrappedAttr error:NULL];
                if ([v isKindOfClass:[NSData class]]) wrapped = (NSData *)v;
            }

            // ProtectionMetadata payload
            NSMutableData *pm = [NSMutableData data];
            appendLEString(pm, cipherSuite, 2);
            appendLEString(pm, kek, 2);
            appendU32LE(pm, (uint32_t)wrapped.length);
            [pm appendData:wrapped];
            appendLEString(pm, @"", 2);   // signature_algorithm
            appendU32LE(pm, 0);            // public_key length
            appendProtectionTrailing(pm, sig, firstChannel);
            [writer _writeRawPacketHeader:TTIOTransportPacketProtectionMetadata
                                     flags:0
                                 datasetId:did
                                auSequence:0
                                   payload:pm];

            // DatasetHeader
            TTIOHDF5Group *hdf5Sig = [[[hdf5Root openGroupNamed:@"study" error:NULL]
                                         openGroupNamed:@"ms_runs" error:NULL]
                                         openGroupNamed:runName error:NULL];
            hdf5Sig = [hdf5Sig openGroupNamed:@"signal_channels" error:NULL];
            NSString *firstSegName = [NSString stringWithFormat:@"%@_segments", firstChannel];
            NSArray *firstSegs =
                [TTIOCompoundIO readGenericFromGroup:hdf5Sig
                                          datasetNamed:firstSegName
                                                fields:channelSegFields()
                                                 error:error];
            if (!firstSegs) return NO;
            uint32_t nSpectra = (uint32_t)firstSegs.count;

            NSString *spectrumClass = readStringAttr(run, @"spectrum_class")
                                        ?: @"TTIOMassSpectrum";
            int64_t acqMode = 0;
            if ([run hasAttributeNamed:@"acquisition_mode"]) {
                id v = [run attributeValueForName:@"acquisition_mode" error:NULL];
                if ([v respondsToSelector:@selector(longLongValue)])
                    acqMode = [v longLongValue];
            }
            if (![writer writeDatasetHeaderWithDatasetId:did
                                                      name:runName
                                           acquisitionMode:(uint8_t)acqMode
                                             spectrumClass:spectrumClass
                                              channelNames:channelNames
                                            instrumentJSON:@"{}"
                                          expectedAUCount:nSpectra
                                                     error:error]) return NO;
            did++;
        }

        // per-genomic-run ProtectionMetadata + DatasetHeader.
        // dataset_id continues from MS so AAD reconstruction matches.
        for (NSString *gRunName in genomicRunNames) {
            id<TTIOStorageGroup> gRun = [genomicRunsGroup openGroupNamed:gRunName error:error];
            id<TTIOStorageGroup> gSig = [gRun openGroupNamed:@"signal_channels" error:error];
            if (!gRun || !gSig) return NO;

            // Genomic encrypts sequences + qualities (M90.1) and, when
            // the run carries them, the SAM tags (M101).
            NSMutableArray<NSString *> *gChannelNames = [NSMutableArray array];
            for (NSString *cn in @[@"sequences", @"qualities", @"tags"]) {
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cn];
                if ([gSig hasChildNamed:segName]) [gChannelNames addObject:cn];
            }
            NSString *gFirst = gChannelNames.firstObject ?: @"sequences";
            NSString *gCipherSuite = readStringAttr(gSig,
                [NSString stringWithFormat:@"%@_algorithm", gFirst]) ?: @"aes-256-gcm";
            NSString *gKek = readStringAttr(gSig,
                [NSString stringWithFormat:@"%@_kek_algorithm", gFirst]) ?: @"";
            NSData *gWrapped = [NSData data];
            NSString *gWrappedAttr = [NSString stringWithFormat:@"%@_wrapped_dek", gFirst];
            if ([gSig hasAttributeNamed:gWrappedAttr]) {
                id v = [gSig attributeValueForName:gWrappedAttr error:NULL];
                if ([v isKindOfClass:[NSData class]]) gWrapped = (NSData *)v;
            }
            NSMutableData *gPm = [NSMutableData data];
            appendLEString(gPm, gCipherSuite, 2);
            appendLEString(gPm, gKek, 2);
            appendU32LE(gPm, (uint32_t)gWrapped.length);
            [gPm appendData:gWrapped];
            appendLEString(gPm, @"", 2);
            appendU32LE(gPm, 0);
            appendProtectionTrailing(gPm, gSig, gFirst);
            [writer _writeRawPacketHeader:TTIOTransportPacketProtectionMetadata
                                     flags:0
                                 datasetId:did
                                auSequence:0
                                   payload:gPm];

            // Discover read count from sequences segments (count only).
            TTIOHDF5Group *gHdf5Sig = [[[hdf5Root openGroupNamed:@"study" error:NULL]
                                          openGroupNamed:@"genomic_runs" error:NULL]
                                          openGroupNamed:gRunName error:NULL];
            gHdf5Sig = [gHdf5Sig openGroupNamed:@"signal_channels" error:NULL];
            NSString *gFirstSegName = [NSString stringWithFormat:@"%@_segments", gFirst];
            NSArray *gFirstSegs =
                [TTIOCompoundIO readGenericFromGroup:gHdf5Sig
                                          datasetNamed:gFirstSegName
                                                fields:channelSegFields()
                                                 error:error];
            if (!gFirstSegs) return NO;
            uint32_t gNReads = (uint32_t)gFirstSegs.count;
            int64_t gAcqMode = 0;
            if ([gRun hasAttributeNamed:@"acquisition_mode"]) {
                id v = [gRun attributeValueForName:@"acquisition_mode" error:NULL];
                if ([v respondsToSelector:@selector(longLongValue)])
                    gAcqMode = [v longLongValue];
            }

            // Build genomic-run metadata JSON (modality / platform /
            // reference_uri / sample_name) per the M89.2 convention.
            NSDictionary *gMetaDict = @{
                @"modality":      readStringAttr(gRun, @"modality")      ?: @"",
                @"platform":      readStringAttr(gRun, @"platform")      ?: @"",
                /* M97: "" when the run has no @read_role attribute. */
                @"read_role":     readStringAttr(gRun, @"read_role")     ?: @"",
                @"reference_uri": readStringAttr(gRun, @"reference_uri") ?: @"",
                @"sample_name":   readStringAttr(gRun, @"sample_name")   ?: @"",
            };
            NSData *gMetaJsonData =
                [NSJSONSerialization dataWithJSONObject:gMetaDict
                                                  options:TTIO_JSON_SORTED_KEYS
                                                    error:NULL];
            NSString *gMetaJson =
                [[NSString alloc] initWithData:gMetaJsonData
                                       encoding:NSUTF8StringEncoding] ?: @"{}";

            if (![writer writeDatasetHeaderWithDatasetId:did
                                                      name:gRunName
                                           acquisitionMode:(uint8_t)gAcqMode
                                             spectrumClass:@"TTIOGenomicRead"
                                              channelNames:gChannelNames
                                            instrumentJSON:gMetaJson
                                          expectedAUCount:gNReads
                                                     error:error]) return NO;
            if (([genomicLayouts[gRunName] isEqualToString:@"blocks_v1"]
                 || [genomicLayouts[gRunName] isEqualToString:@"blocks_v1_grouped"])
                && !emitBlocksV1Sidecars(writer, gRun, gSig, did, error)) {
                return NO;
            }
            did++;
        }

        // AU emission
        did = 1;
        for (NSString *runName in runNames) {
            id<TTIOStorageGroup> run = [msRuns openGroupNamed:runName error:error];
            id<TTIOStorageGroup> sig = [run openGroupNamed:@"signal_channels" error:error];
            id<TTIOStorageGroup> idx = [run openGroupNamed:@"spectrum_index" error:error];
            if (!run || !sig || !idx) return NO;

            NSString *channelNamesStr = readStringAttr(sig, @"channel_names") ?: @"";
            NSArray<NSString *> *channelNames = splitChannelNames(channelNamesStr);
            int64_t acqMode = 0;
            if ([run hasAttributeNamed:@"acquisition_mode"]) {
                id v = [run attributeValueForName:@"acquisition_mode" error:NULL];
                if ([v respondsToSelector:@selector(longLongValue)])
                    acqMode = [v longLongValue];
            }

            // Pre-load all compound rows.
            TTIOHDF5Group *hdf5Run = [[[hdf5Root openGroupNamed:@"study" error:NULL]
                                         openGroupNamed:@"ms_runs" error:NULL]
                                         openGroupNamed:runName error:NULL];
            TTIOHDF5Group *hdf5Sig = [hdf5Run openGroupNamed:@"signal_channels" error:NULL];
            TTIOHDF5Group *hdf5Idx = [hdf5Run openGroupNamed:@"spectrum_index" error:NULL];

            NSMutableDictionary<NSString *, NSArray *> *segsByCh = [NSMutableDictionary dictionary];
            for (NSString *cname in channelNames) {
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cname];
                NSArray *rows =
                    [TTIOCompoundIO readGenericFromGroup:hdf5Sig
                                              datasetNamed:segName
                                                    fields:channelSegFields()
                                                     error:error];
                if (!rows) return NO;
                segsByCh[cname] = rows;
            }

            NSArray *headerSegRows = nil;
            if (headersEncrypted) {
                headerSegRows = [TTIOCompoundIO
                    readGenericFromGroup:hdf5Idx
                              datasetNamed:@"au_header_segments"
                                    fields:headerSegFields()
                                     error:error];
                if (!headerSegRows) return NO;
            }

            NSUInteger n = [segsByCh[channelNames.firstObject] count];
            NSString *spectrumClass = readStringAttr(run, @"spectrum_class")
                                        ?: @"TTIOMassSpectrum";
            uint8_t wireClass = 0;
            if ([spectrumClass isEqualToString:@"TTIONMRSpectrum"])        wireClass = 1;
            else if ([spectrumClass isEqualToString:@"TTIONMR2DSpectrum"]) wireClass = 2;
            else if ([spectrumClass isEqualToString:@"TTIOFreeInductionDecay"]) wireClass = 3;
            else if ([spectrumClass isEqualToString:@"TTIOMSImagePixel"]) wireClass = 4;

            for (NSUInteger i = 0; i < n; i++) {
                // Each channel's ChannelData = {name(u16+bytes), precision(u8),
                // compression(u8), n_elements(u32), data_length(u32),
                // data[IV(12)|TAG(16)|ciphertext]}.
                NSMutableArray<NSData *> *channelPayloads = [NSMutableArray array];
                for (NSString *cname in channelNames) {
                    NSDictionary *row = segsByCh[cname][i];
                    NSData *iv = row[@"iv"];
                    NSData *tag = row[@"tag"];
                    NSData *ct = row[@"ciphertext"];
                    NSMutableData *data = [NSMutableData data];
                    [data appendData:iv];
                    [data appendData:tag];
                    [data appendData:ct];
                    NSMutableData *ch = [NSMutableData data];
                    NSData *nameUtf = [cname dataUsingEncoding:NSUTF8StringEncoding];
                    appendU16LE(ch, (uint16_t)nameUtf.length);
                    [ch appendData:nameUtf];
                    uint8_t precision = TTIOPrecisionFloat64;
                    uint8_t compression = TTIOCompressionNone;
                    [ch appendBytes:&precision length:1];
                    [ch appendBytes:&compression length:1];
                    appendU32LE(ch, [row[@"length"] unsignedIntValue]);
                    appendU32LE(ch, (uint32_t)data.length);
                    [ch appendData:data];
                    [channelPayloads addObject:ch];
                }

                uint16_t flags = TTIOTransportPacketFlagEncrypted;
                NSMutableData *payload = [NSMutableData data];
                if (headersEncrypted) {
                    flags |= TTIOTransportPacketFlagEncryptedHeader;
                    uint8_t sc = wireClass;
                    uint8_t nch = (uint8_t)channelNames.count;
                    [payload appendBytes:&sc length:1];
                    [payload appendBytes:&nch length:1];
                    NSDictionary *hdrRow = headerSegRows[i];
                    [payload appendData:hdrRow[@"iv"]];
                    [payload appendData:hdrRow[@"tag"]];
                    [payload appendData:hdrRow[@"ciphertext"]];
                    for (NSData *ch in channelPayloads) [payload appendData:ch];
                } else {
                    // Plaintext filter header prefix followed by channels,
                    // matching transport-spec §4.3.2.
                    //
                    // Read per-spectrum filter-key values from the
                    // plaintext spectrum_index datasets.
                    id<TTIOStorageDataset> dsRT = [idx openDatasetNamed:@"retention_times" error:NULL];
                    id<TTIOStorageDataset> dsMS = [idx openDatasetNamed:@"ms_levels" error:NULL];
                    id<TTIOStorageDataset> dsPol = [idx openDatasetNamed:@"polarities" error:NULL];
                    id<TTIOStorageDataset> dsPMZ = [idx openDatasetNamed:@"precursor_mzs" error:NULL];
                    id<TTIOStorageDataset> dsPC = [idx openDatasetNamed:@"precursor_charges" error:NULL];
                    id<TTIOStorageDataset> dsBPI = [idx openDatasetNamed:@"base_peak_intensities" error:NULL];
                    const double *rt = (const double *)[[dsRT readAll:NULL] bytes];
                    const int32_t *ms = (const int32_t *)[[dsMS readAll:NULL] bytes];
                    const int32_t *pol = (const int32_t *)[[dsPol readAll:NULL] bytes];
                    const double *pmz = (const double *)[[dsPMZ readAll:NULL] bytes];
                    const int32_t *pc = (const int32_t *)[[dsPC readAll:NULL] bytes];
                    const double *bpi = (const double *)[[dsBPI readAll:NULL] bytes];

                    uint8_t scU8 = wireClass;
                    uint8_t acq = (uint8_t)acqMode;
                    uint8_t msLevel = (uint8_t)(ms[i] & 0xFF);
                    uint8_t polWire = 2;
                    if (pol[i] == 1) polWire = 0;
                    else if (pol[i] == -1) polWire = 1;
                    [payload appendBytes:&scU8 length:1];
                    [payload appendBytes:&acq length:1];
                    [payload appendBytes:&msLevel length:1];
                    [payload appendBytes:&polWire length:1];
                    double rtV = rt[i]; [payload appendBytes:&rtV length:8];
                    double pmzV = pmz[i]; [payload appendBytes:&pmzV length:8];
                    uint8_t pcU8 = (uint8_t)(pc[i] & 0xFF);
                    [payload appendBytes:&pcU8 length:1];
                    double ionM = 0.0; [payload appendBytes:&ionM length:8];
                    double bpiV = bpi[i]; [payload appendBytes:&bpiV length:8];
                    uint8_t nch = (uint8_t)channelNames.count;
                    [payload appendBytes:&nch length:1];
                    for (NSData *ch in channelPayloads) [payload appendData:ch];
                }

                [writer _writeRawPacketHeader:TTIOTransportPacketAccessUnit
                                         flags:flags
                                     datasetId:did
                                    auSequence:(uint32_t)i
                                       payload:payload];
            }

            if (![writer writeEndOfDatasetWithDatasetId:did
                                          finalAUSequence:(uint32_t)n
                                                     error:error]) return NO;
            did++;
        }

        // genomic AU emission. Same dataset_id space as MS so
        // AAD reconstruction (dataset_id + au_sequence) matches.
        for (NSString *gRunName in genomicRunNames) {
            id<TTIOStorageGroup> gRun = [genomicRunsGroup openGroupNamed:gRunName error:error];
            id<TTIOStorageGroup> gSig = [gRun openGroupNamed:@"signal_channels" error:error];
            id<TTIOStorageGroup> gIdx = [gRun openGroupNamed:@"genomic_index" error:error];
            if (!gRun || !gSig || !gIdx) return NO;

            NSMutableArray<NSString *> *gChannelNames = [NSMutableArray array];
            for (NSString *cn in @[@"sequences", @"qualities", @"tags"]) {
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cn];
                if ([gSig hasChildNamed:segName]) [gChannelNames addObject:cn];
            }

            // Pre-load encrypted segments per genomic channel.
            TTIOHDF5Group *gHdf5Run = [[[hdf5Root openGroupNamed:@"study" error:NULL]
                                          openGroupNamed:@"genomic_runs" error:NULL]
                                          openGroupNamed:gRunName error:NULL];
            TTIOHDF5Group *gHdf5Sig = [gHdf5Run openGroupNamed:@"signal_channels" error:NULL];
            NSMutableDictionary<NSString *, NSArray *> *gSegsByCh = [NSMutableDictionary dictionary];
            for (NSString *cn in gChannelNames) {
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cn];
                NSArray *rows =
                    [TTIOCompoundIO readGenericFromGroup:gHdf5Sig
                                              datasetNamed:segName
                                                    fields:channelSegFields()
                                                     error:error];
                if (!rows) return NO;
                gSegsByCh[cn] = rows;
            }

            // Plaintext genomic_index columns (NOT encrypted per M90.1).
            id<TTIOStorageDataset> gPosDs = [gIdx openDatasetNamed:@"positions" error:NULL];
            id<TTIOStorageDataset> gMqDs  = [gIdx openDatasetNamed:@"mapping_qualities" error:NULL];
            id<TTIOStorageDataset> gFlDs  = [gIdx openDatasetNamed:@"flags" error:NULL];
            NSData *gPosData = gPosDs ? [gPosDs readAll:NULL] : nil;
            NSData *gMqData  = gMqDs  ? [gMqDs readAll:NULL]  : nil;
            NSData *gFlData  = gFlDs  ? [gFlDs readAll:NULL]  : nil;
            const int64_t *gPosArr = (const int64_t *)gPosData.bytes;
            const uint8_t *gMqArr  = (const uint8_t *)gMqData.bytes;
            const uint32_t *gFlArr = (const uint32_t *)gFlData.bytes;
            // L1 (Task #82 Phase B.1): chromosomes are stored as
            // chromosome_ids (uint16) + chromosome_names (compound).
            TTIOHDF5Group *gHdf5Idx = [gHdf5Run openGroupNamed:@"genomic_index" error:NULL];
            id<TTIOStorageDataset> gChromIdsDs = [gHdf5Idx openDatasetNamed:@"chromosome_ids" error:NULL];
            NSData *gChromIdsData = gChromIdsDs ? [gChromIdsDs readAll:NULL] : nil;
            NSArray<TTIOCompoundField *> *chromFields = @[
                [TTIOCompoundField fieldWithName:@"name"
                                            kind:TTIOCompoundFieldKindVLString],
            ];
            NSArray *gChromNameRows =
                [TTIOCompoundIO readGenericFromGroup:gHdf5Idx
                                          datasetNamed:@"chromosome_names"
                                                fields:chromFields
                                                 error:error];
            if (!gChromIdsData || !gChromNameRows) return NO;
            NSMutableArray<NSString *> *gChromNameTable =
                [NSMutableArray arrayWithCapacity:gChromNameRows.count];
            for (NSDictionary *row in gChromNameRows) {
                id v = row[@"name"];
                if ([v isKindOfClass:[NSData class]])
                    v = [[NSString alloc] initWithData:v encoding:NSUTF8StringEncoding];
                [gChromNameTable addObject:(NSString *)v ?: @""];
            }
            const uint16_t *gChromIds = (const uint16_t *)gChromIdsData.bytes;
            NSUInteger gNChromIds = gChromIdsData.length / sizeof(uint16_t);
            NSMutableArray<NSDictionary *> *gChromRows =
                [NSMutableArray arrayWithCapacity:gNChromIds];
            for (NSUInteger i = 0; i < gNChromIds; i++) {
                NSUInteger idx = gChromIds[i];
                NSString *name = idx < gChromNameTable.count
                    ? gChromNameTable[idx] : @"";
                [gChromRows addObject:@{@"value": name}];
            }
            int64_t gAcqMode = 0;
            if ([gRun hasAttributeNamed:@"acquisition_mode"]) {
                id v = [gRun attributeValueForName:@"acquisition_mode" error:NULL];
                if ([v respondsToSelector:@selector(longLongValue)])
                    gAcqMode = [v longLongValue];
            }
            NSUInteger gN = [gSegsByCh[gChannelNames.firstObject] count];

            for (NSUInteger i = 0; i < gN; i++) {
                NSMutableArray<TTIOTransportChannelData *> *channels = [NSMutableArray array];
                for (NSString *cname in gChannelNames) {
                    NSDictionary *row = gSegsByCh[cname][i];
                    NSData *iv = row[@"iv"];
                    NSData *tag = row[@"tag"];
                    NSData *ct = row[@"ciphertext"];
                    NSMutableData *data = [NSMutableData data];
                    [data appendData:iv];
                    [data appendData:tag];
                    [data appendData:ct];
                    TTIOTransportChannelData *ch =
                        [[TTIOTransportChannelData alloc] initWithName:cname
                                                              precision:TTIOPrecisionUInt8
                                                            compression:TTIOCompressionNone
                                                              nElements:[row[@"length"] unsignedIntValue]
                                                                   data:data];
                    [channels addObject:ch];
                }
                NSString *chrom = @"";
                id chromVal = gChromRows[i][@"value"];
                if ([chromVal isKindOfClass:[NSString class]]) chrom = chromVal;
                else if ([chromVal isKindOfClass:[NSData class]]) {
                    chrom = [[NSString alloc] initWithData:chromVal
                                                   encoding:NSUTF8StringEncoding] ?: @"";
                }
                int64_t pos = gPosArr ? gPosArr[i] : 0;
                uint8_t mapq = gMqArr ? gMqArr[i] : 0;
                uint16_t flags = (uint16_t)((gFlArr ? gFlArr[i] : 0) & 0xFFFFu);
                TTIOAccessUnit *au =
                    [[TTIOAccessUnit alloc] initWithSpectrumClass:5
                                                  acquisitionMode:(uint8_t)gAcqMode
                                                          msLevel:0
                                                         polarity:2
                                                    retentionTime:0.0
                                                      precursorMz:0.0
                                                  precursorCharge:0
                                                      ionMobility:0.0
                                                basePeakIntensity:0.0
                                                         channels:channels
                                                           pixelX:0 pixelY:0 pixelZ:0
                                                       chromosome:chrom
                                                         position:pos
                                                   mappingQuality:mapq
                                                            flags:flags];
                NSData *payload = [au encode];
                [writer _writeRawPacketHeader:TTIOTransportPacketAccessUnit
                                         flags:TTIOTransportPacketFlagEncrypted
                                     datasetId:did
                                    auSequence:(uint32_t)i
                                       payload:payload];
            }
            if (![writer writeEndOfDatasetWithDatasetId:did
                                          finalAUSequence:(uint32_t)gN
                                                     error:error]) return NO;
            did++;
        }

        if (![writer writeEndOfStreamWithError:error]) return NO;
        ok = YES;
    }
    @finally {
        [sp close];
    }
    return ok;
}

// ---------------------------------------------------------------- reader

static DatasetAccumulator *makeAcc(NSString *name, uint8_t acqMode,
                                      NSString *spectrumClass,
                                      NSArray<NSString *> *channelNames)
{
    DatasetAccumulator *d = [[DatasetAccumulator alloc] init];
    d.name = name;
    d.acquisitionMode = acqMode;
    d.spectrumClass = spectrumClass;
    d.channelNames = [NSMutableArray arrayWithArray:channelNames];
    d.channelSegments = [NSMutableDictionary dictionary];
    for (NSString *c in channelNames) d.channelSegments[c] = [NSMutableArray array];
    d.headerSegments = [NSMutableArray array];
    d.usedEncryptedHeaders = NO;
    // pre-init genomic accumulators. isGenomic flips on when
    // the DATASET_HEADER spectrum_class is TTIOGenomicRead.
    d.isGenomic = [spectrumClass isEqualToString:@"TTIOGenomicRead"];
    d.genomicChromosomes = [NSMutableArray array];
    d.genomicPositions = [NSMutableArray array];
    d.genomicMapqs = [NSMutableArray array];
    d.genomicFlags = [NSMutableArray array];
    d.genomicMetadataJSON = @"";
    d.blockSidecars = [NSMutableArray array];
    return d;
}


// Parse a ProtectionMetadata payload: string cipher_suite, string
// kek_algorithm, u32 wrapped_dek_len, wrapped_dek bytes, string
// signature_algorithm, u32 public_key_len, public_key bytes.
static ProtectionMeta *parseProtection(NSData *payload)
{
    const uint8_t *b = (const uint8_t *)payload.bytes;
    NSUInteger len = payload.length;
    NSUInteger off = 0;
    ProtectionMeta *pm = [[ProtectionMeta alloc] init];
    if (off + 2 > len) return pm;
    uint16_t csLen = readU16LE(&b[off]); off += 2;
    pm.cipherSuite = [[NSString alloc] initWithBytes:&b[off] length:csLen
                                             encoding:NSUTF8StringEncoding];
    off += csLen;
    if (off + 2 > len) return pm;
    uint16_t kekLen = readU16LE(&b[off]); off += 2;
    pm.kekAlgorithm = [[NSString alloc] initWithBytes:&b[off] length:kekLen
                                              encoding:NSUTF8StringEncoding];
    off += kekLen;
    if (off + 4 > len) return pm;
    uint32_t wLen = readU32LE(&b[off]); off += 4;
    pm.wrappedDek = [NSData dataWithBytes:&b[off] length:wLen];
    off += wLen;
    // FD-1 Phase A: consume signature_algorithm + public_key, then decode
    // the trailing recipient block (absent on pre-Phase-A packets, where
    // sig/pk are empty and no bytes remain -> no additional recipients).
    pm.additionalRecipients = @[];
    if (off + 2 <= len) {
        uint16_t sigLen = readU16LE(&b[off]); off += 2;
        off += sigLen;                              // signature_algorithm
        if (off + 4 <= len) {
            uint32_t pkLen = readU32LE(&b[off]); off += 4;
            off += pkLen;                           // public_key
            pm.additionalRecipients = decodeRecipientBlock(b, len, &off);
            // FD-1 C-2a: anything after the recipient block is server_kek_id.
            if (off + 2 <= len) {
                uint16_t kidLen = readU16LE(&b[off]); off += 2;
                if (off + kidLen <= len) {
                    pm.serverKekId = [[NSString alloc]
                        initWithBytes:&b[off] length:kidLen
                             encoding:NSUTF8StringEncoding];
                    off += kidLen;
                }
            }
        }
    }
    return pm;
}


// Parse channel-data list given a payload cursor. Returns YES and
// appends one TTIOChannelSegment per channel into
// ``acc.channelSegments[cname]``. On ENCRYPTED AUs each channel's
// ``data`` is IV(12) || TAG(16) || ciphertext.
static BOOL parseEncryptedChannels(const uint8_t *buf, NSUInteger len,
                                     NSUInteger *offsetInOut,
                                     uint8_t nChannels,
                                     DatasetAccumulator *acc,
                                     NSUInteger spectrumLengthHint,
                                     NSError **error)
{
    (void)spectrumLengthHint;
    NSUInteger off = *offsetInOut;
    for (uint8_t c = 0; c < nChannels; c++) {
        if (off + 2 > len) { if (error) *error = makeErr(30, @"AU truncated: channel name len"); return NO; }
        uint16_t nameLen = readU16LE(&buf[off]); off += 2;
        if (off + nameLen > len) { if (error) *error = makeErr(30, @"AU truncated: channel name"); return NO; }
        NSString *cname = [[NSString alloc] initWithBytes:&buf[off] length:nameLen
                                                  encoding:NSUTF8StringEncoding];
        off += nameLen;
        if (off + 10 > len) { if (error) *error = makeErr(30, @"AU truncated: channel hdr"); return NO; }
        uint8_t precision = buf[off++];
        uint8_t compression = buf[off++];
        uint32_t nElements = readU32LE(&buf[off]); off += 4;
        uint32_t dataLen = readU32LE(&buf[off]); off += 4;
        if (off + dataLen > len) { if (error) *error = makeErr(30, @"AU truncated: channel data"); return NO; }
        if (dataLen < 28) { if (error) *error = makeErr(30, @"encrypted channel data shorter than IV+TAG"); return NO; }
        NSData *iv = [NSData dataWithBytes:&buf[off] length:12];
        NSData *tag = [NSData dataWithBytes:&buf[off + 12] length:16];
        NSData *ciphertext = [NSData dataWithBytes:&buf[off + 28] length:dataLen - 28];
        off += dataLen;
        (void)precision; (void)compression;

        NSMutableArray<TTIOChannelSegment *> *segs = acc.channelSegments[cname];
        if (!segs) {
            segs = [NSMutableArray array];
            acc.channelSegments[cname] = segs;
            if (![acc.channelNames containsObject:cname])
                [acc.channelNames addObject:cname];
        }
        // offset = cumulative sum of prior lengths (reconstructed).
        uint64_t prior = 0;
        for (TTIOChannelSegment *s in segs) prior += s.length;
        [segs addObject:[[TTIOChannelSegment alloc]
            initWithOffset:prior
                     length:nElements
                         iv:iv
                        tag:tag
                 ciphertext:ciphertext]];
    }
    *offsetInOut = off;
    return YES;
}


static BOOL ingestAU(TTIOTransportPacketRecord *record,
                      DatasetAccumulator *acc,
                      NSError **error)
{
    uint16_t flags = record.header.flags;
    BOOL encHeader = (flags & TTIOTransportPacketFlagEncryptedHeader) != 0;
    BOOL encChannel = (flags & TTIOTransportPacketFlagEncrypted) != 0;
    if (!encChannel) {
        if (error) *error = makeErr(30, @"ingestAU on plaintext AU");
        return NO;
    }
    acc.usedEncryptedHeaders = encHeader;
    const uint8_t *buf = (const uint8_t *)record.payload.bytes;
    NSUInteger len = record.payload.length;
    NSUInteger off = 0;

    if (encHeader) {
        // Wire: spectrum_class(u8) n_channels(u8) IV(12) TAG(16)
        // encrypted_header(36) [channels...]
        if (len < 1 + 1 + 12 + 16 + 36) {
            if (error) *error = makeErr(30, @"encrypted-header AU too short");
            return NO;
        }
        off++; // spectrum_class (already set on DatasetHeader)
        uint8_t nChannels = buf[off++];
        NSData *hdrIV = [NSData dataWithBytes:&buf[off] length:12]; off += 12;
        NSData *hdrTag = [NSData dataWithBytes:&buf[off] length:16]; off += 16;
        NSData *hdrCt = [NSData dataWithBytes:&buf[off] length:36]; off += 36;
        [acc.headerSegments addObject:[[TTIOHeaderSegment alloc]
            initWithIV:hdrIV tag:hdrTag ciphertext:hdrCt]];
        if (!parseEncryptedChannels(buf, len, &off, nChannels, acc, 0, error))
            return NO;
    } else {
        // Plaintext filter header followed by channels. Skip the
        // 38-byte header prefix — we don't need it for reconstruction
        // (the plaintext index arrays round-trip through the writer
        // via the dataset header; for encrypted-channel-only mode the
        // source file still carries plaintext spectrum_index, so the
        // MATERIALISED file will re-derive retention_times etc.)
        //
        // Actually wait — on the receiver side we DON'T have the
        // source's plaintext index arrays; they travel on the wire
        // inside the AU fixed prefix. We need to capture them.
        if (len < 38) {
            if (error) *error = makeErr(30, @"plaintext-header AU too short");
            return NO;
        }
        off++; // spectrum_class
        uint8_t acq = buf[off++];
        uint8_t msLevel = buf[off++];
        uint8_t polarityWire = buf[off++];
        double rt; memcpy(&rt, &buf[off], 8); off += 8;
        double pmz; memcpy(&pmz, &buf[off], 8); off += 8;
        uint8_t pc = buf[off++];
        double ionMob; memcpy(&ionMob, &buf[off], 8); off += 8;
        double bpi; memcpy(&bpi, &buf[off], 8); off += 8;
        uint8_t nChannels = buf[off++];

        // capture genomic suffix when this is a genomic AU
        // so the materialiser can rebuild the plaintext
        // genomic_index columns. Decode through TTIOAccessUnit
        // since it knows the chrom-len-prefixed layout.
        if (acc.isGenomic) {
            NSError *gAuErr = nil;
            TTIOAccessUnit *gAu =
                [TTIOAccessUnit decodeFromBytes:buf length:len error:&gAuErr];
            if (gAu && gAu.spectrumClass == 5) {
                [acc.genomicChromosomes addObject:(gAu.chromosome ?: @"")];
                [acc.genomicPositions addObject:@(gAu.position)];
                [acc.genomicMapqs addObject:@(gAu.mappingQuality)];
                [acc.genomicFlags addObject:@((uint32_t)gAu.flags)];
            }
        }

        // Stash plaintext filter values in a parallel array held on
        // the accumulator under keys the writer can consume.
        NSMutableArray *rts = acc.channelSegments[@"__rt__"]
                                ?: ((acc.channelSegments[@"__rt__"] = [NSMutableArray array]));
        NSMutableArray *msLevels = acc.channelSegments[@"__ms_level__"]
                                ?: ((acc.channelSegments[@"__ms_level__"] = [NSMutableArray array]));
        NSMutableArray *pols = acc.channelSegments[@"__polarity__"]
                                ?: ((acc.channelSegments[@"__polarity__"] = [NSMutableArray array]));
        NSMutableArray *pmzs = acc.channelSegments[@"__precursor_mz__"]
                                ?: ((acc.channelSegments[@"__precursor_mz__"] = [NSMutableArray array]));
        NSMutableArray *pcs = acc.channelSegments[@"__precursor_charge__"]
                                ?: ((acc.channelSegments[@"__precursor_charge__"] = [NSMutableArray array]));
        NSMutableArray *bpis = acc.channelSegments[@"__base_peak__"]
                                ?: ((acc.channelSegments[@"__base_peak__"] = [NSMutableArray array]));
        [(id)rts addObject:@(rt)];
        [(id)msLevels addObject:@(msLevel)];
        int32_t polInt = 0;
        if (polarityWire == 0) polInt = 1;
        else if (polarityWire == 1) polInt = -1;
        [(id)pols addObject:@(polInt)];
        [(id)pmzs addObject:@(pmz)];
        [(id)pcs addObject:@(pc)];
        [(id)bpis addObject:@(bpi)];
        (void)acq; (void)ionMob;

        if (!parseEncryptedChannels(buf, len, &off, nChannels, acc, 0, error))
            return NO;
    }
    return YES;
}


// Materialise one blocks_v1 per-AU-encrypted genomic run from its
// sidecar packets and AU stream, in the shape the stream writer
// creates it, so decrypt-in-place restores it (transport-spec §4.24).
static BOOL writeBlocksV1GenomicRun(id<TTIOStorageGroup> gRunsGroup,
                                    DatasetAccumulator *acc,
                                    ProtectionMeta *pm,
                                    NSError **error)
{
    NSDictionary *sc = acc.runSidecar;
    NSDictionary *gMeta = @{};
    if (acc.genomicMetadataJSON.length > 0) {
        NSData *jd =
            [acc.genomicMetadataJSON dataUsingEncoding:NSUTF8StringEncoding];
        id parsed = [NSJSONSerialization JSONObjectWithData:jd
                                                      options:0 error:NULL];
        if ([parsed isKindOfClass:[NSDictionary class]]) gMeta = parsed;
    }
    id<TTIOStorageGroup> gRun =
        [gRunsGroup createGroupNamed:acc.name error:error];
    if (!gRun) return NO;
    if (![gRun setAttributeValue:@((int64_t)acc.acquisitionMode)
                           forName:@"acquisition_mode" error:error]) return NO;
    NSString *modality = [gMeta[@"modality"] isKindOfClass:[NSString class]]
        ? gMeta[@"modality"] : @"";
    if (![gRun setAttributeValue:modality forName:@"modality" error:error])
        return NO;
    if (![gRun setAttributeValue:@((int64_t)5)
                           forName:@"spectrum_class" error:error]) return NO;
    for (NSString *fld in @[@"reference_uri", @"platform"]) {
        NSString *v = [gMeta[fld] isKindOfClass:[NSString class]]
            ? gMeta[fld] : @"";
        if (![gRun setAttributeValue:v forName:fld error:error]) return NO;
    }
    NSString *gRole = [gMeta[@"read_role"] isKindOfClass:[NSString class]]
        ? gMeta[@"read_role"] : @"";
    if (gRole.length > 0
        && ![gRun setAttributeValue:gRole forName:@"read_role" error:error])
        return NO;
    NSString *sample = [gMeta[@"sample_name"] isKindOfClass:[NSString class]]
        ? gMeta[@"sample_name"] : @"";
    if (![gRun setAttributeValue:sample forName:@"sample_name" error:error])
        return NO;
    if (![gRun setAttributeValue:sc[@"read_count"] ?: @0
                           forName:@"read_count" error:error]) return NO;
    if (![gRun setAttributeValue:sc[@"base_count"] ?: @0
                           forName:@"base_count" error:error]) return NO;
    if (![gRun setAttributeValue:sc[@"layout"] ?: @""
                           forName:@"layout" error:error]) return NO;
    if (![gRun setAttributeValue:sc[@"block_policy"] ?: @""
                           forName:@"block_policy" error:error]) return NO;
    NSDictionary *attrs = sc[@"attrs"] ?: @{};
    if ([attrs[@"ref_diff_slice_bytes"] longLongValue] != 0
        && ![gRun setAttributeValue:
                @([attrs[@"ref_diff_slice_bytes"] longLongValue])
                            forName:@"ref_diff_slice_bytes" error:error])
        return NO;
    if ([attrs[@"opt_disable_qualities_v5"] longLongValue] != 0
        && ![gRun setAttributeValue:@((int64_t)1)
                            forName:@"opt_disable_qualities_v5" error:error])
        return NO;
    if ([attrs[@"reference_md5s"] isKindOfClass:[NSString class]]
        && ![gRun setAttributeValue:attrs[@"reference_md5s"]
                            forName:@"reference_md5s" error:error])
        return NO;

    NSArray<NSDictionary *> *sidecars = [acc.blockSidecars
        sortedArrayUsingComparator:^NSComparisonResult(NSDictionary *a,
                                                       NSDictionary *b) {
            return [a[@"block_index"] compare:b[@"block_index"]];
        }];
    id<TTIOStorageGroup> blocks = [gRun createGroupNamed:@"blocks"
                                                    error:error];
    if (!blocks) return NO;
    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:sidecars.count];
    for (NSDictionary *bs in sidecars) {
        NSMutableDictionary *row = [NSMutableDictionary dictionary];
        row[@"read_start"] = bs[@"read_start"];
        row[@"n_reads"] = @([bs[@"n_reads"] unsignedIntValue]);
        row[@"base_start"] = bs[@"base_start"];
        row[@"n_bases"] = bs[@"n_bases"];
        NSDictionary *chans = bs[@"channels"];
        for (NSString *ch in chans) {
            NSArray *triple = chans[ch];
            row[[ch stringByAppendingString:@"_off"]] = triple[0];
            row[[ch stringByAppendingString:@"_len"]] = triple[1];
            row[[ch stringByAppendingString:@"_codec"]] =
                @([triple[2] unsignedIntValue]);
        }
        [rows addObject:row];
    }
    // The index carries the tags triple when the sender's did (M101).
    id<TTIOStorageDataset> idxDs = [blocks
        createCompoundDatasetNamed:@"index"
                            fields:[TTIOGenomicStreamWriter indexFieldsForRows:rows]
                             count:0
                        extendable:YES
                         chunkRows:1024
                             error:error];
    if (!idxDs) return NO;
    if (rows.count > 0 && ![idxDs appendData:rows error:error]) return NO;

    id<TTIOStorageGroup> gIdx = [gRun createGroupNamed:@"genomic_index"
                                                  error:error];
    if (!gIdx) return NO;
    NSString *firstCh = acc.channelNames.firstObject;
    NSArray<TTIOChannelSegment *> *firstSegs = acc.channelSegments[firstCh];
    NSUInteger nReads = firstSegs.count;
    NSArray<NSString *> *chromTable = sc[@"chromosome_names"] ?: @[];
    NSMutableDictionary<NSString *, NSNumber *> *nameToId =
        [NSMutableDictionary dictionary];
    for (NSUInteger i = 0; i < chromTable.count; i++) {
        nameToId[chromTable[i]] = @(i);
    }
    NSMutableData *lenData = [NSMutableData dataWithLength:nReads * 4];
    uint32_t *lenP = lenData.mutableBytes;
    for (NSUInteger i = 0; i < nReads; i++) lenP[i] = firstSegs[i].length;
    NSMutableData *posData = [NSMutableData dataWithLength:nReads * 8];
    NSMutableData *mqData = [NSMutableData dataWithLength:nReads];
    NSMutableData *flData = [NSMutableData dataWithLength:nReads * 4];
    NSMutableData *cidData = [NSMutableData dataWithLength:nReads * 2];
    int64_t *posP = posData.mutableBytes;
    uint8_t *mqP = mqData.mutableBytes;
    uint32_t *flP = flData.mutableBytes;
    uint16_t *cidP = cidData.mutableBytes;
    for (NSUInteger i = 0; i < nReads; i++) {
        posP[i] = i < acc.genomicPositions.count
            ? [acc.genomicPositions[i] longLongValue] : 0;
        mqP[i] = i < acc.genomicMapqs.count
            ? (uint8_t)[acc.genomicMapqs[i] unsignedCharValue] : 0;
        flP[i] = i < acc.genomicFlags.count
            ? (uint32_t)[acc.genomicFlags[i] unsignedIntValue] : 0;
        NSString *c = i < acc.genomicChromosomes.count
            ? acc.genomicChromosomes[i] : @"";
        NSNumber *slot = nameToId[c ?: @""];
        if (slot == nil) {
            if (error) *error = makeErr(33,
                @"AU chromosome '%@' is not in the sidecar "
                @"chromosome_names table", c);
            return NO;
        }
        cidP[i] = (uint16_t)slot.unsignedShortValue;
    }
    // Extendable, chunked and zlib-6 like the stream writer's
    // genomic_index arrays (kIndexArrayChunk).
    NSArray *idxArrays = @[
        @[@"lengths", @(TTIOPrecisionUInt32), lenData],
        @[@"positions", @(TTIOPrecisionInt64), posData],
        @[@"mapping_qualities", @(TTIOPrecisionUInt8), mqData],
        @[@"flags", @(TTIOPrecisionUInt32), flData],
        @[@"chromosome_ids", @(TTIOPrecisionUInt16), cidData],
    ];
    for (NSArray *a in idxArrays) {
        id<TTIOStorageDataset> ds =
            [gIdx createDatasetNamed:a[0]
                            precision:(TTIOPrecision)[a[1] integerValue]
                               length:0
                            chunkSize:65536
                          compression:TTIOCompressionZlib
                     compressionLevel:6
                           extendable:YES
                                error:error];
        if (!ds || ![ds appendData:a[2] error:error]) return NO;
    }
    if ([sc[@"layout"] isEqualToString:@"blocks_v1_grouped"]) {
        // M103: the permutation back to input order, from the sidecars.
        NSMutableData *ii = [NSMutableData data];
        for (NSDictionary *bs in sidecars) {
            NSData *part = bs[@"input_index"];
            if (part == nil) {
                if (error) *error = makeErr(34,
                    @"genomic run '%@': blocks_v1_grouped stream without "
                    @"input_index in its BlockSidecars", acc.name);
                return NO;
            }
            [ii appendData:part];
        }
        id<TTIOStorageDataset> iiDs =
            [gIdx createDatasetNamed:@"input_index"
                            precision:TTIOPrecisionUInt32
                               length:0
                            chunkSize:65536
                          compression:TTIOCompressionZlib
                     compressionLevel:6
                           extendable:YES
                                error:error];
        if (!iiDs || ![iiDs appendData:ii error:error]) return NO;
    }
    NSArray<TTIOCompoundField *> *nameFields = @[
        [TTIOCompoundField fieldWithName:@"name"
                                    kind:TTIOCompoundFieldKindVLString],
    ];
    id<TTIOStorageDataset> namesDs =
        [gIdx createCompoundDatasetNamed:@"chromosome_names"
                                   fields:nameFields
                                    count:chromTable.count
                                    error:error];
    if (!namesDs) return NO;
    NSMutableArray *nameRows =
        [NSMutableArray arrayWithCapacity:chromTable.count];
    for (NSString *n in chromTable) [nameRows addObject:@{ @"name": n }];
    if (![namesDs writeAll:nameRows error:error]) return NO;

    id<TTIOStorageGroup> gSig =
        [gRun createGroupNamed:@"signal_channels" error:error];
    if (!gSig) return NO;
    NSArray<TTIOCompoundField *> *chFields = channelSegFields();
    for (NSString *cname in acc.channelNames) {
        NSArray<TTIOChannelSegment *> *segs = acc.channelSegments[cname];
        NSString *segName = [NSString stringWithFormat:@"%@_segments", cname];
        id<TTIOStorageDataset> ds =
            [gSig createCompoundDatasetNamed:segName
                                       fields:chFields
                                        count:segs.count
                                        error:error];
        if (!ds) return NO;
        NSMutableArray *segRows = [NSMutableArray arrayWithCapacity:segs.count];
        for (TTIOChannelSegment *s in segs) {
            [segRows addObject:@{
                @"offset": @(s.offset), @"length": @(s.length),
                @"iv": s.iv, @"tag": s.tag, @"ciphertext": s.ciphertext,
            }];
        }
        if (![ds writeAll:segRows error:error]) return NO;
        if (![gSig setAttributeValue:(pm.cipherSuite ?: @"aes-256-gcm")
                               forName:[NSString stringWithFormat:
                                           @"%@_algorithm", cname]
                                 error:error]) return NO;
        if (pm && pm.wrappedDek.length > 0) {
            [gSig setAttributeValue:pm.wrappedDek
                              forName:[NSString stringWithFormat:
                                          @"%@_wrapped_dek", cname]
                                error:NULL];
            [gSig setAttributeValue:(pm.kekAlgorithm ?: @"")
                              forName:[NSString stringWithFormat:
                                          @"%@_kek_algorithm", cname]
                                error:NULL];
            NSData *recipBlock = encodeRecipientBlock(pm.additionalRecipients);
            if (recipBlock.length > 0) {
                [gSig setAttributeValue:recipBlock
                                forName:[NSString stringWithFormat:
                                            @"%@_wrapped_dek_recipients",
                                            cname]
                                  error:NULL];
            }
            if (pm.serverKekId.length > 0) {
                [gSig setAttributeValue:pm.serverKekId
                                forName:[NSString stringWithFormat:
                                            @"%@_server_kek_id", cname]
                                  error:NULL];
            }
        }
    }

    id<TTIOStorageGroup> mateGroup = nil;
    for (NSDictionary *centry in (NSArray *)(sc[@"channels"] ?: @[])) {
        NSString *ch = centry[@"name"];
        id<TTIOStorageGroup> parent;
        NSString *dsName;
        if ([ch isEqualToString:@"mate_info"]) {
            mateGroup = [gSig createGroupNamed:@"mate_info" error:error];
            if (!mateGroup) return NO;
            parent = mateGroup;
            dsName = @"inline_v2";
        } else {
            parent = gSig;
            dsName = ch;
        }
        int64_t codec = [centry[@"compression"] longLongValue];
        id<TTIOStorageDataset> ds =
            [parent createDatasetNamed:dsName
                             precision:TTIOPrecisionUInt8
                                length:0
                             chunkSize:[TTIOGenomicStreamWriter channelChunk]
                           compression:(codec == 0 ? TTIOCompressionZlib
                                                    : TTIOCompressionNone)
                      compressionLevel:6
                            extendable:YES
                                 error:error];
        if (!ds) return NO;
        if (![ds setAttributeValue:@(codec)
                           forName:@"compression" error:error]) return NO;
        NSDictionary *extra = centry[@"extra_attrs"] ?: @{};
        for (NSString *k in [[extra allKeys]
                sortedArrayUsingSelector:@selector(compare:)]) {
            if (![ds setAttributeValue:extra[k] forName:k error:error])
                return NO;
        }
        for (NSDictionary *bs in sidecars) {
            NSData *data = ((NSDictionary *)bs[@"blobs"])[ch];
            if (data.length > 0
                && ![ds appendData:data error:error]) return NO;
        }
    }
    // The stream writer creates mate_info and its chrom_names table
    // at close even when no mate blob was written; mirror that.
    NSArray *mateNames = sc[@"mate_chrom_names"] ?: @[];
    if (mateNames.count > 0) {
        if (mateGroup == nil) {
            mateGroup = [gSig createGroupNamed:@"mate_info" error:error];
            if (!mateGroup) return NO;
        }
        id<TTIOStorageDataset> mateDs =
            [mateGroup createCompoundDatasetNamed:@"chrom_names"
                                            fields:nameFields
                                             count:mateNames.count
                                             error:error];
        if (!mateDs) return NO;
        NSMutableArray *mateRows =
            [NSMutableArray arrayWithCapacity:mateNames.count];
        for (NSString *n in mateNames) [mateRows addObject:@{ @"name": n }];
        if (![mateDs writeAll:mateRows error:error]) return NO;
    }
    return YES;
}


static BOOL writeEncryptedFile(NSString *path,
                                  NSString *providerName,
                                  NSString *title,
                                  NSString *isa,
                                  NSArray<NSString *> *featureList,
                                  NSDictionary<NSNumber *, ProtectionMeta *> *protection,
                                  NSMutableDictionary<NSNumber *, DatasetAccumulator *> *datasets,
                                  NSError **error)
{
    id<TTIOStorageProvider> sp =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                    mode:TTIOStorageOpenModeCreate
                                                provider:providerName
                                                   error:error];
    if (!sp) return NO;
    @try {
        if (![sp.providerName isEqualToString:@"hdf5"]) {
            if (error) *error = makeErr(4,
                @"encrypted transport reader currently requires HDF5 provider");
            return NO;
        }
        id<TTIOStorageGroup> root = [sp rootGroupWithError:error];
        if (!root) return NO;

        // Feature flags.
        NSData *featuresJson =
            [NSJSONSerialization dataWithJSONObject:featureList
                                              options:0 error:error];
        if (!featuresJson) return NO;
        NSString *featuresStr = [[NSString alloc] initWithData:featuresJson
                                                       encoding:NSUTF8StringEncoding];
        if (![root setAttributeValue:@"1.1" forName:@"ttio_format_version"
                                error:error]) return NO;
        if (![root setAttributeValue:featuresStr forName:@"ttio_features"
                                error:error]) return NO;

        id<TTIOStorageGroup> study =
            [root createGroupNamed:@"study" error:error];
        if (!study) return NO;
        if (![study setAttributeValue:(title ?: @"") forName:@"title"
                                  error:error]) return NO;
        if (![study setAttributeValue:(isa ?: @"")
                                forName:@"isa_investigation_id"
                                  error:error]) return NO;

        id<TTIOStorageGroup> msRuns =
            [study createGroupNamed:@"ms_runs" error:error];
        if (!msRuns) return NO;
        // split into MS vs genomic datasets. Each gets its own
        // /study/{ms_runs,genomic_runs}/ subtree.
        NSArray *sortedDids = [datasets.allKeys sortedArrayUsingSelector:
                                @selector(compare:)];
        NSMutableArray *msDids = [NSMutableArray array];
        NSMutableArray *genomicDids = [NSMutableArray array];
        for (NSNumber *did in sortedDids) {
            if (datasets[did].isGenomic) [genomicDids addObject:did];
            else [msDids addObject:did];
        }
        NSMutableArray *runNamesList = [NSMutableArray array];
        for (NSNumber *did in msDids) [runNamesList addObject:datasets[did].name];
        if (![msRuns setAttributeValue:[runNamesList componentsJoinedByString:@","]
                                 forName:@"_run_names" error:error]) return NO;

        for (NSNumber *didKey in msDids) {
            DatasetAccumulator *acc = datasets[didKey];
            id<TTIOStorageGroup> run =
                [msRuns createGroupNamed:acc.name error:error];
            if (!run) return NO;
            if (![run setAttributeValue:@(acc.acquisitionMode)
                                   forName:@"acquisition_mode" error:error]) return NO;
            if (![run setAttributeValue:(acc.spectrumClass ?: @"TTIOMassSpectrum")
                                   forName:@"spectrum_class" error:error]) return NO;
            NSUInteger spectrumCount = acc.headerSegments.count;
            if (spectrumCount == 0) {
                spectrumCount = [acc.channelSegments[acc.channelNames.firstObject] count];
            }
            if (![run setAttributeValue:@((int64_t)spectrumCount)
                                   forName:@"spectrum_count" error:error]) return NO;

            // Empty instrument_config group (matches TTIOPerAUFile output).
            id<TTIOStorageGroup> cfg = [run createGroupNamed:@"instrument_config" error:error];
            if (!cfg) return NO;
            for (NSString *f in @[@"manufacturer", @"model", @"serial_number",
                                    @"source_type", @"analyzer_type", @"detector_type"]) {
                [cfg setAttributeValue:@"" forName:f error:NULL];
            }

            id<TTIOStorageGroup> sig = [run createGroupNamed:@"signal_channels" error:error];
            if (!sig) return NO;

            // Only the "real" channel names — drop the internal "__rt__" etc.
            NSMutableArray *realChannels = [NSMutableArray array];
            for (NSString *c in acc.channelNames) {
                if (![c hasPrefix:@"__"]) [realChannels addObject:c];
            }
            if (![sig setAttributeValue:[realChannels componentsJoinedByString:@","]
                                   forName:@"channel_names" error:error]) return NO;

            // Write each channel's segments compound via the provider
            // (createCompoundDataset + writeAll).
            ProtectionMeta *pm = protection[didKey];
            NSArray<TTIOCompoundField *> *chFields = channelSegFields();
            for (NSString *cname in realChannels) {
                NSArray<TTIOChannelSegment *> *segs = acc.channelSegments[cname];
                NSString *segName = [NSString stringWithFormat:@"%@_segments", cname];
                id<TTIOStorageDataset> ds =
                    [sig createCompoundDatasetNamed:segName
                                               fields:chFields
                                                count:segs.count
                                                error:error];
                if (!ds) return NO;
                NSMutableArray *rows = [NSMutableArray arrayWithCapacity:segs.count];
                for (TTIOChannelSegment *s in segs) {
                    [rows addObject:@{
                        @"offset": @(s.offset), @"length": @(s.length),
                        @"iv": s.iv, @"tag": s.tag, @"ciphertext": s.ciphertext,
                    }];
                }
                if (![ds writeAll:rows error:error]) return NO;
                if (![sig setAttributeValue:(pm.cipherSuite ?: @"aes-256-gcm")
                                       forName:[NSString stringWithFormat:@"%@_algorithm", cname]
                                         error:error]) return NO;
                if (pm && pm.wrappedDek.length > 0) {
                    [sig setAttributeValue:pm.wrappedDek
                                     forName:[NSString stringWithFormat:@"%@_wrapped_dek", cname]
                                       error:NULL];
                    [sig setAttributeValue:(pm.kekAlgorithm ?: @"")
                                     forName:[NSString stringWithFormat:@"%@_kek_algorithm", cname]
                                       error:NULL];
                    NSData *recipBlock = encodeRecipientBlock(pm.additionalRecipients);
                    if (recipBlock.length > 0) {
                        [sig setAttributeValue:recipBlock
                                       forName:[NSString stringWithFormat:@"%@_wrapped_dek_recipients", cname]
                                         error:NULL];
                    }
                    if (pm.serverKekId.length > 0) {
                        [sig setAttributeValue:pm.serverKekId
                                       forName:[NSString stringWithFormat:@"%@_server_kek_id", cname]
                                         error:NULL];
                    }
                }
            }

            // Spectrum index: plaintext offsets + lengths from the
            // first channel's segments; (plaintext-header mode only)
            // retention_times / ms_levels / polarities / pmzs / pcs /
            // bpis from the captured arrays; (encrypted-header mode)
            // au_header_segments compound.
            id<TTIOStorageGroup> idx = [run createGroupNamed:@"spectrum_index" error:error];
            if (!idx) return NO;
            NSArray<TTIOChannelSegment *> *firstSegs =
                acc.channelSegments[realChannels.firstObject];
            [idx setAttributeValue:@((int64_t)firstSegs.count)
                             forName:@"count" error:NULL];

            NSMutableData *offData = [NSMutableData data];
            NSMutableData *lenData = [NSMutableData data];
            for (TTIOChannelSegment *s in firstSegs) {
                uint64_t o = s.offset; [offData appendBytes:&o length:8];
                uint32_t l = s.length; [lenData appendBytes:&l length:4];
            }
            id<TTIOStorageDataset> offDs =
                [idx createDatasetNamed:@"offsets"
                                precision:TTIOPrecisionInt64
                                   length:firstSegs.count
                                chunkSize:0
                              compression:TTIOCompressionNone
                         compressionLevel:0
                                    error:error];
            if (!offDs) return NO;
            [offDs writeAll:offData error:error];
            id<TTIOStorageDataset> lenDs =
                [idx createDatasetNamed:@"lengths"
                                precision:TTIOPrecisionUInt32
                                   length:firstSegs.count
                                chunkSize:0
                              compression:TTIOCompressionNone
                         compressionLevel:0
                                    error:error];
            if (!lenDs) return NO;
            [lenDs writeAll:lenData error:error];

            if (acc.usedEncryptedHeaders) {
                // Write au_header_segments compound; omit plaintext
                // arrays since opt_encrypted_au_headers is set.
                NSArray<TTIOCompoundField *> *hdrFields = headerSegFields();
                id<TTIOStorageDataset> hdrDs =
                    [idx createCompoundDatasetNamed:@"au_header_segments"
                                               fields:hdrFields
                                                count:acc.headerSegments.count
                                                error:error];
                if (!hdrDs) return NO;
                NSMutableArray *rows = [NSMutableArray array];
                for (TTIOHeaderSegment *s in acc.headerSegments) {
                    [rows addObject:@{
                        @"iv": s.iv, @"tag": s.tag, @"ciphertext": s.ciphertext,
                    }];
                }
                if (![hdrDs writeAll:rows error:error]) return NO;
            } else {
                // Plaintext-header mode: write the six plaintext
                // index arrays from the captured per-AU filter fields.
                NSArray *rts = acc.channelSegments[@"__rt__"];
                NSArray *msLevels = acc.channelSegments[@"__ms_level__"];
                NSArray *pols = acc.channelSegments[@"__polarity__"];
                NSArray *pmzs = acc.channelSegments[@"__precursor_mz__"];
                NSArray *pcs = acc.channelSegments[@"__precursor_charge__"];
                NSArray *bpis = acc.channelSegments[@"__base_peak__"];

                NSMutableData *rtData = [NSMutableData data];
                for (NSNumber *n in rts) { double v = n.doubleValue; [rtData appendBytes:&v length:8]; }
                NSMutableData *msData = [NSMutableData data];
                for (NSNumber *n in msLevels) { int32_t v = n.intValue; [msData appendBytes:&v length:4]; }
                NSMutableData *polData = [NSMutableData data];
                for (NSNumber *n in pols) { int32_t v = n.intValue; [polData appendBytes:&v length:4]; }
                NSMutableData *pmzData = [NSMutableData data];
                for (NSNumber *n in pmzs) { double v = n.doubleValue; [pmzData appendBytes:&v length:8]; }
                NSMutableData *pcData = [NSMutableData data];
                for (NSNumber *n in pcs) { int32_t v = n.intValue; [pcData appendBytes:&v length:4]; }
                NSMutableData *bpiData = [NSMutableData data];
                for (NSNumber *n in bpis) { double v = n.doubleValue; [bpiData appendBytes:&v length:8]; }

                struct { NSString *name; TTIOPrecision prec; NSData *data; } plain[] = {
                    {@"retention_times", TTIOPrecisionFloat64, rtData},
                    {@"ms_levels", TTIOPrecisionInt32, msData},
                    {@"polarities", TTIOPrecisionInt32, polData},
                    {@"precursor_mzs", TTIOPrecisionFloat64, pmzData},
                    {@"precursor_charges", TTIOPrecisionInt32, pcData},
                    {@"base_peak_intensities", TTIOPrecisionFloat64, bpiData},
                };
                for (size_t i = 0; i < sizeof(plain) / sizeof(plain[0]); i++) {
                    id<TTIOStorageDataset> pDs =
                        [idx createDatasetNamed:plain[i].name
                                        precision:plain[i].prec
                                           length:firstSegs.count
                                        chunkSize:0
                                      compression:TTIOCompressionNone
                                 compressionLevel:0
                                            error:error];
                    if (!pDs) return NO;
                    [pDs writeAll:plain[i].data error:error];
                }
            }
        }
        // write study/genomic_runs/<name>/ for each genomic
        // dataset captured on the wire. Encrypted segments round-trip
        // verbatim; plaintext genomic_index columns are reconstructed
        // from the per-AU genomic suffix arrays captured during ingest.
        if (genomicDids.count > 0) {
            id<TTIOStorageGroup> gRunsGroup =
                [study createGroupNamed:@"genomic_runs" error:error];
            if (!gRunsGroup) return NO;
            NSMutableArray *gNamesList = [NSMutableArray array];
            for (NSNumber *did in genomicDids) [gNamesList addObject:datasets[did].name];
            if (![gRunsGroup setAttributeValue:[gNamesList componentsJoinedByString:@","]
                                          forName:@"_run_names" error:error]) return NO;

            for (NSNumber *didKey in genomicDids) {
                DatasetAccumulator *acc = datasets[didKey];
                if (acc.runSidecar != nil) {
                    if (!writeBlocksV1GenomicRun(gRunsGroup, acc,
                                                 protection[didKey],
                                                 error)) return NO;
                    continue;
                }
                id<TTIOStorageGroup> gRun =
                    [gRunsGroup createGroupNamed:acc.name error:error];
                if (!gRun) return NO;
                if (![gRun setAttributeValue:@(acc.acquisitionMode)
                                       forName:@"acquisition_mode" error:error]) return NO;
                if (![gRun setAttributeValue:@"TTIOGenomicRead"
                                       forName:@"spectrum_class" error:error]) return NO;
                NSDictionary *gMeta = @{};
                if (acc.genomicMetadataJSON.length > 0) {
                    NSData *jd = [acc.genomicMetadataJSON dataUsingEncoding:NSUTF8StringEncoding];
                    id parsed = [NSJSONSerialization JSONObjectWithData:jd
                                                                  options:0 error:NULL];
                    if ([parsed isKindOfClass:[NSDictionary class]]) gMeta = parsed;
                }
                for (NSString *fld in @[@"modality", @"platform",
                                          @"reference_uri", @"sample_name"]) {
                    NSString *v = [gMeta[fld] isKindOfClass:[NSString class]]
                        ? gMeta[fld] : @"";
                    [gRun setAttributeValue:v forName:fld error:NULL];
                }
                /* M97: absent on the container when "" on the wire. */
                NSString *gRole = [gMeta[@"read_role"] isKindOfClass:[NSString class]]
                    ? gMeta[@"read_role"] : @"";
                if (gRole.length > 0)
                    [gRun setAttributeValue:gRole forName:@"read_role" error:NULL];

                id<TTIOStorageGroup> gSig =
                    [gRun createGroupNamed:@"signal_channels" error:error];
                if (!gSig) return NO;
                NSMutableArray *gRealCh = [NSMutableArray array];
                for (NSString *c in acc.channelNames) {
                    if (![c hasPrefix:@"__"]) [gRealCh addObject:c];
                }
                if (![gSig setAttributeValue:[gRealCh componentsJoinedByString:@","]
                                       forName:@"channel_names" error:error]) return NO;

                ProtectionMeta *pm = protection[didKey];
                NSArray<TTIOCompoundField *> *chFields = channelSegFields();
                for (NSString *cname in gRealCh) {
                    NSArray<TTIOChannelSegment *> *segs = acc.channelSegments[cname];
                    NSString *segName = [NSString stringWithFormat:@"%@_segments", cname];
                    id<TTIOStorageDataset> ds =
                        [gSig createCompoundDatasetNamed:segName
                                                   fields:chFields
                                                    count:segs.count
                                                    error:error];
                    if (!ds) return NO;
                    NSMutableArray *rows = [NSMutableArray arrayWithCapacity:segs.count];
                    for (TTIOChannelSegment *s in segs) {
                        [rows addObject:@{
                            @"offset": @(s.offset), @"length": @(s.length),
                            @"iv": s.iv, @"tag": s.tag, @"ciphertext": s.ciphertext,
                        }];
                    }
                    if (![ds writeAll:rows error:error]) return NO;
                    if (![gSig setAttributeValue:(pm.cipherSuite ?: @"aes-256-gcm")
                                           forName:[NSString stringWithFormat:@"%@_algorithm", cname]
                                             error:error]) return NO;
                    if (pm && pm.wrappedDek.length > 0) {
                        [gSig setAttributeValue:pm.wrappedDek
                                          forName:[NSString stringWithFormat:@"%@_wrapped_dek", cname]
                                            error:NULL];
                        [gSig setAttributeValue:(pm.kekAlgorithm ?: @"")
                                          forName:[NSString stringWithFormat:@"%@_kek_algorithm", cname]
                                            error:NULL];
                        NSData *recipBlock = encodeRecipientBlock(pm.additionalRecipients);
                        if (recipBlock.length > 0) {
                            [gSig setAttributeValue:recipBlock
                                            forName:[NSString stringWithFormat:@"%@_wrapped_dek_recipients", cname]
                                              error:NULL];
                        }
                        if (pm.serverKekId.length > 0) {
                            [gSig setAttributeValue:pm.serverKekId
                                            forName:[NSString stringWithFormat:@"%@_server_kek_id", cname]
                                              error:NULL];
                        }
                    }
                }

                id<TTIOStorageGroup> gIdx =
                    [gRun createGroupNamed:@"genomic_index" error:error];
                if (!gIdx) return NO;
                NSArray<TTIOChannelSegment *> *gFirstSegs =
                    acc.channelSegments[gRealCh.firstObject];
                NSUInteger gNReads = gFirstSegs.count;
                [gIdx setAttributeValue:@((int64_t)gNReads)
                                  forName:@"count" error:NULL];
                NSMutableData *gOffData = [NSMutableData data];
                NSMutableData *gLenData = [NSMutableData data];
                for (TTIOChannelSegment *s in gFirstSegs) {
                    uint64_t o = s.offset; [gOffData appendBytes:&o length:8];
                    uint32_t l = s.length; [gLenData appendBytes:&l length:4];
                }
                id<TTIOStorageDataset> gOffDs =
                    [gIdx createDatasetNamed:@"offsets"
                                    precision:TTIOPrecisionUInt64
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gOffDs) return NO;
                [gOffDs writeAll:gOffData error:error];
                id<TTIOStorageDataset> gLenDs =
                    [gIdx createDatasetNamed:@"lengths"
                                    precision:TTIOPrecisionUInt32
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gLenDs) return NO;
                [gLenDs writeAll:gLenData error:error];

                NSMutableData *gPosData = [NSMutableData data];
                NSMutableData *gMqData  = [NSMutableData data];
                NSMutableData *gFlData  = [NSMutableData data];
                for (NSUInteger i = 0; i < gNReads; i++) {
                    int64_t p = i < acc.genomicPositions.count
                        ? [acc.genomicPositions[i] longLongValue] : 0;
                    [gPosData appendBytes:&p length:8];
                    uint8_t mq = i < acc.genomicMapqs.count
                        ? (uint8_t)[acc.genomicMapqs[i] unsignedCharValue] : 0;
                    [gMqData appendBytes:&mq length:1];
                    uint32_t fl = i < acc.genomicFlags.count
                        ? (uint32_t)[acc.genomicFlags[i] unsignedIntValue] : 0;
                    [gFlData appendBytes:&fl length:4];
                }
                id<TTIOStorageDataset> gPosDs =
                    [gIdx createDatasetNamed:@"positions"
                                    precision:TTIOPrecisionInt64
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gPosDs) return NO;
                [gPosDs writeAll:gPosData error:error];
                id<TTIOStorageDataset> gMqDs =
                    [gIdx createDatasetNamed:@"mapping_qualities"
                                    precision:TTIOPrecisionUInt8
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gMqDs) return NO;
                [gMqDs writeAll:gMqData error:error];
                id<TTIOStorageDataset> gFlDs =
                    [gIdx createDatasetNamed:@"flags"
                                    precision:TTIOPrecisionUInt32
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gFlDs) return NO;
                [gFlDs writeAll:gFlData error:error];

                // L1 (Task #82 Phase B.1): write chromosomes as
                // chromosome_ids (uint16) + chromosome_names (compound).
                NSMutableDictionary<NSString *, NSNumber *> *gNameToId = [NSMutableDictionary dictionary];
                NSMutableArray<NSString *> *gNamesInOrder = [NSMutableArray array];
                NSMutableData *gIdsData = [NSMutableData dataWithLength:gNReads * sizeof(uint16_t)];
                uint16_t *gIdsBuf = (uint16_t *)gIdsData.mutableBytes;
                for (NSUInteger i = 0; i < gNReads; i++) {
                    NSString *c = i < acc.genomicChromosomes.count
                        ? acc.genomicChromosomes[i] : @"";
                    if (!c) c = @"";
                    NSNumber *slot = gNameToId[c];
                    if (!slot) {
                        slot = @(gNamesInOrder.count);
                        gNameToId[c] = slot;
                        [gNamesInOrder addObject:c];
                    }
                    gIdsBuf[i] = (uint16_t)slot.unsignedShortValue;
                }
                id<TTIOStorageDataset> gCidsDs =
                    [gIdx createDatasetNamed:@"chromosome_ids"
                                    precision:TTIOPrecisionUInt16
                                       length:gNReads
                                    chunkSize:0
                                  compression:TTIOCompressionNone
                             compressionLevel:0
                                        error:error];
                if (!gCidsDs) return NO;
                if (![gCidsDs writeAll:gIdsData error:error]) return NO;
                NSArray<TTIOCompoundField *> *gNameFields = @[
                    [TTIOCompoundField fieldWithName:@"name"
                                                kind:TTIOCompoundFieldKindVLString],
                ];
                id<TTIOStorageDataset> gNamesDs =
                    [gIdx createCompoundDatasetNamed:@"chromosome_names"
                                               fields:gNameFields
                                                count:gNamesInOrder.count
                                                error:error];
                if (!gNamesDs) return NO;
                NSMutableArray *gNameRowsW = [NSMutableArray arrayWithCapacity:gNamesInOrder.count];
                for (NSString *n in gNamesInOrder) {
                    [gNameRowsW addObject:@{ @"name": n }];
                }
                if (![gNamesDs writeAll:gNameRowsW error:error]) return NO;
            }
        }
    }
    @finally {
        [sp close];
    }
    return YES;
}


+ (BOOL)readEncryptedToPath:(NSString *)outputPath
                fromStream:(NSData *)streamData
               providerName:(NSString *)providerName
                      error:(NSError **)error
{
    TTIOTransportReader *reader =
        [[TTIOTransportReader alloc] initWithData:streamData];
    NSArray<TTIOTransportPacketRecord *> *packets =
        [reader readAllPacketsWithError:error];
    if (!packets) return NO;

    NSString *title = @"";
    NSString *isa = @"";
    NSMutableArray *features = [NSMutableArray array];
    NSMutableDictionary<NSNumber *, DatasetAccumulator *> *datasets =
        [NSMutableDictionary dictionary];
    NSMutableDictionary<NSNumber *, ProtectionMeta *> *protection =
        [NSMutableDictionary dictionary];

    for (TTIOTransportPacketRecord *rec in packets) {
        TTIOTransportPacketType t = rec.header.packetType;
        const uint8_t *b = (const uint8_t *)rec.payload.bytes;
        NSUInteger len = rec.payload.length;
        NSUInteger off = 0;

        if (t == TTIOTransportPacketStreamHeader) {
            if (off + 2 > len) continue;
            uint16_t vl = readU16LE(&b[off]); off += 2 + vl;   // format_version
            if (off + 2 > len) continue;
            uint16_t tl = readU16LE(&b[off]); off += 2;
            title = [[NSString alloc] initWithBytes:&b[off] length:tl encoding:NSUTF8StringEncoding];
            off += tl;
            if (off + 2 > len) continue;
            uint16_t il = readU16LE(&b[off]); off += 2;
            isa = [[NSString alloc] initWithBytes:&b[off] length:il encoding:NSUTF8StringEncoding];
            off += il;
            if (off + 2 > len) continue;
            uint16_t nFeat = readU16LE(&b[off]); off += 2;
            for (uint16_t i = 0; i < nFeat; i++) {
                if (off + 2 > len) break;
                uint16_t fl = readU16LE(&b[off]); off += 2;
                NSString *fn = [[NSString alloc] initWithBytes:&b[off] length:fl
                                                      encoding:NSUTF8StringEncoding];
                off += fl;
                [features addObject:fn];
            }
            // n_datasets follows but we discover it from DatasetHeaders.
        } else if (t == TTIOTransportPacketProtectionMetadata) {
            ProtectionMeta *pm = parseProtection(rec.payload);
            protection[@(rec.header.datasetId)] = pm;
        } else if (t == TTIOTransportPacketDatasetHeader) {
            // dataset_id u16, name (str2), acq_mode u8, spectrum_class
            // (str2), n_channels u8, channel_names (str2 × n), instr_json
            // (str4), expected_au_count u32.
            if (off + 2 > len) continue;
            uint16_t did = readU16LE(&b[off]); off += 2;
            if (off + 2 > len) continue;
            uint16_t nl = readU16LE(&b[off]); off += 2;
            NSString *runName = [[NSString alloc] initWithBytes:&b[off] length:nl
                                                       encoding:NSUTF8StringEncoding];
            off += nl;
            if (off + 1 > len) continue;
            uint8_t acq = b[off++];
            if (off + 2 > len) continue;
            uint16_t scLen = readU16LE(&b[off]); off += 2;
            NSString *spectrumClass =
                [[NSString alloc] initWithBytes:&b[off] length:scLen
                                       encoding:NSUTF8StringEncoding];
            off += scLen;
            if (off + 1 > len) continue;
            uint8_t nch = b[off++];
            NSMutableArray *chNames = [NSMutableArray arrayWithCapacity:nch];
            for (uint8_t i = 0; i < nch; i++) {
                if (off + 2 > len) break;
                uint16_t cl = readU16LE(&b[off]); off += 2;
                [chNames addObject:[[NSString alloc] initWithBytes:&b[off]
                                                              length:cl
                                                            encoding:NSUTF8StringEncoding]];
                off += cl;
            }
            // read instrument_json — for genomic datasets it
            // carries the per-run metadata (modality / platform /
            // reference_uri / sample_name) which we re-emit on the
            // materialised genomic_run group.
            NSString *instrJson = @"";
            if (off + 4 <= len) {
                uint32_t jl = readU32LE(&b[off]); off += 4;
                if (off + jl <= len) {
                    instrJson = [[NSString alloc] initWithBytes:&b[off] length:jl
                                                         encoding:NSUTF8StringEncoding] ?: @"";
                    off += jl;
                }
            }
            // expected_au_count u32 (skipped — discovered from AUs).
            DatasetAccumulator *acc = makeAcc(runName, acq, spectrumClass, chNames);
            acc.genomicMetadataJSON = instrJson;
            datasets[@(did)] = acc;
        } else if (t == TTIOTransportPacketGenomicRunSidecar) {
            DatasetAccumulator *acc = datasets[@(rec.header.datasetId)];
            if (!acc) {
                if (error) *error = makeErr(31,
                    @"GenomicRunSidecar for unknown dataset_id %u",
                    (unsigned)rec.header.datasetId);
                return NO;
            }
            acc.runSidecar = decodeGenomicRunSidecar(rec.payload);
            if (!acc.runSidecar) {
                if (error) *error = makeErr(31,
                    @"malformed GenomicRunSidecar payload");
                return NO;
            }
        } else if (t == TTIOTransportPacketBlockSidecar) {
            DatasetAccumulator *acc = datasets[@(rec.header.datasetId)];
            if (!acc) {
                if (error) *error = makeErr(32,
                    @"BlockSidecar for unknown dataset_id %u",
                    (unsigned)rec.header.datasetId);
                return NO;
            }
            NSDictionary *bs = decodeBlockSidecar(rec.payload);
            if (!bs) {
                if (error) *error = makeErr(32,
                    @"malformed BlockSidecar payload");
                return NO;
            }
            [acc.blockSidecars addObject:bs];
        } else if (t == TTIOTransportPacketAccessUnit) {
            DatasetAccumulator *acc = datasets[@(rec.header.datasetId)];
            if (!acc) {
                if (error) *error = makeErr(30, @"AU for unknown dataset_id %u",
                                              (unsigned)rec.header.datasetId);
                return NO;
            }
            if (!ingestAU(rec, acc, error)) return NO;
        }
        // EndOfDataset / EndOfStream: skip.
    }

    NSMutableSet *featureSet = [NSMutableSet setWithArray:features];
    // The blocks_v1 wire token is transport-scoped and never a
    // container feature flag.
    [featureSet removeObject:@"transport_blocks_v1"];
    [featureSet removeObject:@"transport_blocks_v1_grouped"];
    [featureSet addObject:@"opt_per_au_encryption"];
    BOOL anyHeaderEncrypted = NO;
    for (NSNumber *k in datasets) {
        if (datasets[k].usedEncryptedHeaders) { anyHeaderEncrypted = YES; break; }
    }
    if (anyHeaderEncrypted) [featureSet addObject:@"opt_encrypted_au_headers"];
    NSArray *featureList = [featureSet.allObjects sortedArrayUsingSelector:@selector(compare:)];

    return writeEncryptedFile(outputPath, providerName, title, isa,
                                featureList, protection, datasets, error);
}

@end

// FD-1 Phase A-4: expose the recipient-block codec to the cross-language
// conformance test (see TTIOEncryptedTransport+Conformance.h).
@implementation TTIOEncryptedTransport (Conformance)

+ (NSData *)ttioConformanceEncodeRecipientBlock:(NSArray<NSDictionary *> *)additional
{
    return encodeRecipientBlock(additional);
}

+ (NSArray<NSDictionary *> *)ttioConformanceDecodeRecipientBlock:(NSData *)block
{
    NSUInteger off = 0;
    return decodeRecipientBlock((const uint8_t *)block.bytes, block.length, &off);
}

+ (NSData *)ttioConformanceEncodeTrailing:(NSArray<NSDictionary *> *)additional
                              serverKekId:(NSString *)serverKekId
{
    NSMutableData *out = [NSMutableData data];
    BOOL haveAdd = additional.count > 0;
    BOOL haveKid = serverKekId.length > 0;
    if (!haveAdd && !haveKid) return out;
    if (haveAdd) {
        [out appendData:encodeRecipientBlock(additional)];
    } else {
        appendU16LE(out, 0);   // additional_recipient_count = 0
    }
    if (haveKid) appendLEString(out, serverKekId, 2);
    return out;
}

@end
