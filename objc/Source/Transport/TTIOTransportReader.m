/*
 * TTIOTransportReader.m
 * TTI-O Objective-C Implementation
 *
 * Classes:       TTIOTransportPacketRecord, TTIOTransportReader
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Transport/TTIOTransportReader.h
 *
 * Transport stream reader. Walks a transport byte stream into
 * (header, payload) packet records and (optionally) materialises a
 * fresh .tio file from the stream. Validates magic / version /
 * CRC-32C and rejects out-of-order AU sequences.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "TTIOTransportReader.h"
#import "TTIOTransportReader+Internal.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Dataset/TTIOWrittenRun.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOBulkV2Blobs.h"
#import "Genomics/TTIOReferenceImport.h"
#import "Image/TTIOMSImage.h"
#import "Dataset/TTIOIdentification.h"
#import "Dataset/TTIOQuantification.h"
#import "Dataset/TTIOSubject.h"
#import "Dataset/TTIOSample.h"
#import "Transport/TTIOArrowIpcCodec.h"
#import "Providers/TTIOHDF5Provider.h"
#import "Providers/TTIOStorageProtocols.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import <hdf5.h>
#import "Codecs/TTIORans.h"        // rANS wire codec dispatch
#import "Codecs/TTIOBasePack.h"    // BASE_PACK wire codec dispatch
#import "Codecs/TTIOFloatDeltaZstd.h"  // codec 17 spectral AU channels
#import "ValueClasses/TTIOEnums.h"
#import <objc/runtime.h>
#import <string.h>
#import <zlib.h>
#import <zstd.h>

// ---------------------------------------------------------------- LE helpers

static inline uint16_t readU16(const uint8_t *b)
{
    return (uint16_t)((uint32_t)b[0] | ((uint32_t)b[1] << 8));
}

static inline uint32_t readU32(const uint8_t *b)
{
    return (uint32_t)b[0]
         | ((uint32_t)b[1] << 8)
         | ((uint32_t)b[2] << 16)
         | ((uint32_t)b[3] << 24);
}

static inline uint64_t readU64(const uint8_t *b)
{
    uint64_t lo = (uint64_t)readU32(b);
    uint64_t hi = (uint64_t)readU32(b + 4);
    return lo | (hi << 32);
}

static inline double readF64(const uint8_t *b)
{
    uint64_t bits = readU64(b);
    double v;
    memcpy(&v, &bits, 8);
    return v;
}

// v0.11 Task 3.3: zlib inflate for REFERENCE_CHROMOSOME encoding=1.
// expectedLen is the sequence_length field on the wire — used to
// size the destination buffer exactly so the decoder catches a size
// mismatch deterministically. Returns nil on rc != Z_OK or on a
// short / long inflation.
static NSData *zlibInflateExact(NSData *deflated,
                                 NSUInteger expectedLen,
                                 NSError **error)
{
    if (expectedLen == 0) return [NSData data];
    NSMutableData *out = [NSMutableData dataWithLength:expectedLen];
    uLongf destLen = (uLongf)expectedLen;
    int rc = uncompress((Bytef *)out.mutableBytes, &destLen,
                          (const Bytef *)deflated.bytes,
                          (uLong)deflated.length);
    if (rc != Z_OK) {
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                                 [NSString stringWithFormat:
                                     @"REFERENCE_CHROMOSOME zlib inflate failed: rc=%d",
                                     rc]}];
        return nil;
    }
    if ((NSUInteger)destLen != expectedLen) {
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                                 [NSString stringWithFormat:
                                     @"REFERENCE_CHROMOSOME zlib payload inflated "
                                     @"to %lu bytes; expected %lu",
                                     (unsigned long)destLen,
                                     (unsigned long)expectedLen]}];
        return nil;
    }
    return out;
}

// zlib inflate without a known plaintext size (IMAGE_PIXEL payloads:
// processed-mode entry counts are only knowable after inflation).
// Grows the destination geometrically up to a 128 MB cap.
static NSData *zlibInflateGrow(const uint8_t *src, NSUInteger srcLen,
                               NSError **error)
{
    NSUInteger cap = MAX((NSUInteger)4096, srcLen * 4);
    for (;;) {
        NSMutableData *out = [NSMutableData dataWithLength:cap];
        uLongf destLen = (uLongf)cap;
        int rc = uncompress((Bytef *)out.mutableBytes, &destLen,
                              (const Bytef *)src, (uLong)srcLen);
        if (rc == Z_OK) {
            [out setLength:destLen];
            return out;
        }
        if (rc == Z_BUF_ERROR && cap < ((NSUInteger)1 << 27)) {
            cap *= 2;
            continue;
        }
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                                 [NSString stringWithFormat:
                                     @"IMAGE_PIXEL zlib inflate failed: rc=%d", rc]}];
        return nil;
    }
}

static NSString *readLEString(const uint8_t *bytes, NSUInteger length,
                              NSUInteger *offset, int width)
{
    NSUInteger off = *offset;
    uint32_t strLen = 0;
    if (width == 2) {
        if (off + 2 > length) return nil;
        strLen = readU16(&bytes[off]);
        off += 2;
    } else {
        if (off + 4 > length) return nil;
        strLen = readU32(&bytes[off]);
        off += 4;
    }
    if (off + strLen > length) return nil;
    NSString *s = [[NSString alloc] initWithBytes:&bytes[off]
                                            length:strLen
                                          encoding:NSUTF8StringEncoding];
    *offset = off + strLen;
    return s ?: @"";
}

// ---------------------------------------------------------------- record

@implementation TTIOTransportPacketRecord

- (instancetype)initWithHeader:(TTIOTransportPacketHeader *)h payload:(NSData *)p
{
    if ((self = [super init])) {
        _header = h;
        _payload = [p copy];
    }
    return self;
}

@end

// ---------------------------------------------------------------- reader

@implementation TTIOTransportReader
{
    NSData *_buffer;
}

- (instancetype)initWithInputPath:(NSString *)path
{
    if ((self = [super init])) {
        _buffer = [NSData dataWithContentsOfFile:path];
    }
    return self;
}

- (instancetype)initWithData:(NSData *)data
{
    if ((self = [super init])) {
        _buffer = [data copy];
    }
    return self;
}

- (NSArray<TTIOTransportPacketRecord *> *)readAllPacketsWithError:(NSError **)error
{
    if (!_buffer) {
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorTruncated
                                             userInfo:@{NSLocalizedDescriptionKey:
                             @"empty input"}];
        return nil;
    }
    const uint8_t *bytes = (const uint8_t *)_buffer.bytes;
    NSUInteger length = _buffer.length;
    NSUInteger offset = 0;
    NSMutableArray<TTIOTransportPacketRecord *> *out = [NSMutableArray array];

    while (offset < length) {
        if (length - offset < TTIOTransportHeaderSize) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorTruncated
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"truncated packet header"}];
            return nil;
        }
        TTIOTransportPacketHeader *hdr =
            [TTIOTransportPacketHeader decodeFromBytes:&bytes[offset]
                                                 length:length - offset
                                                  error:error];
        if (!hdr) return nil;
        offset += TTIOTransportHeaderSize;

        if (length - offset < hdr.payloadLength) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorTruncated
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"truncated payload"}];
            return nil;
        }
        NSData *payload = [NSData dataWithBytes:&bytes[offset]
                                          length:hdr.payloadLength];
        offset += hdr.payloadLength;

        if (hdr.flags & TTIOTransportPacketFlagHasChecksum) {
            if (length - offset < 4) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorTruncated
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     @"truncated CRC-32C"}];
                return nil;
            }
            uint32_t expected = readU32(&bytes[offset]);
            offset += 4;
            uint32_t actual = TTIOTransportCRC32C((const uint8_t *)payload.bytes,
                                                     payload.length);
            if (expected != actual) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorChecksumFailed
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:@"CRC-32C mismatch: expected 0x%08x, got 0x%08x",
                                         expected, actual]}];
                return nil;
            }
        }

        // Forward-compat (transport-spec §6, v0.11 task 0.7): tolerate
        // unknown packet types so v0.10 readers can ingest v0.11+
        // streams. The header was length-prefixed so the payload (and
        // CRC if present) was already consumed above — just log and
        // record the packet so callers can observe it. The materialize
        // loop early-continues on the same condition.
        if (!TTIOTransportIsKnownPacketType(hdr.packetTypeByte)) {
            NSLog(@"TTIOTransportReader: skipping unknown packet type 0x%02x",
                  (unsigned)hdr.packetTypeByte);
        }

        TTIOTransportPacketRecord *rec =
            [[TTIOTransportPacketRecord alloc] initWithHeader:hdr payload:payload];
        [out addObject:rec];

        if (hdr.packetType == TTIOTransportPacketEndOfStream) break;
    }

    return out;
}

// v0.11 Task 5.3: shared HDF5 cube-group writer for raman/ir image
// embed. Mirrors the inline MS image_cube layout already in
// writeTtioToPath: but generalised over group name + axis dataset
// name + extra scalar attrs. Used by the modality dispatcher when
// materialising a stream that carries Raman or IR pixels.
//
// Returns NO + populates *error on HDF5 failure. studyGid is owned by
// the caller; this helper does not close it.
static BOOL writeImageCubeGroupAtStudy(hid_t studyGid,
                                       const char *groupName,
                                       uint32_t width, uint32_t height,
                                       uint32_t bins, NSUInteger tileSize,
                                       double pixelSizeX, double pixelSizeY,
                                       NSString *scanPattern,
                                       NSData *axis,
                                       const void *cubeBytes,
                                       NSDictionary<NSString *, NSNumber *> *doubleAttrs,
                                       NSDictionary<NSString *, NSNumber *> *intAttrs,
                                       NSDictionary<NSString *, NSString *> *stringAttrs,
                                       NSError **error)
{
    hid_t g = H5Gcreate2(studyGid, groupName,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) {
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                             [NSString stringWithFormat:
                                 @"image embed: H5Gcreate2 %s failed",
                                 groupName]}];
        return NO;
    }
    hsize_t dims[3]  = { (hsize_t)height, (hsize_t)width, (hsize_t)bins };
    hsize_t chunk[3] = { (hsize_t)MIN(tileSize, (NSUInteger)height),
                         (hsize_t)MIN(tileSize, (NSUInteger)width),
                         (hsize_t)bins };
    hid_t space = H5Screate_simple(3, dims, NULL);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 3, chunk);
    H5Pset_shuffle(plist);   /* byte-shuffle before deflate; matches the cube writers */
    H5Pset_deflate(plist, 6);
    hid_t did = H5Dcreate2(g, "intensity",
                           H5T_NATIVE_DOUBLE, space,
                           H5P_DEFAULT, plist, H5P_DEFAULT);
    if (did < 0) {
        H5Pclose(plist); H5Sclose(space); H5Gclose(g);
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                             @"image embed: H5Dcreate2 intensity failed"}];
        return NO;
    }
    herr_t st = H5Dwrite(did, H5T_NATIVE_DOUBLE,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT, cubeBytes);
    H5Dclose(did); H5Pclose(plist); H5Sclose(space);
    if (st < 0) {
        H5Gclose(g);
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                             @"image embed: H5Dwrite intensity failed"}];
        return NO;
    }
    // 1-D wavenumbers / mz_axis dataset (consumed by RamanImage /
    // IRImage / MSImage round-trip). MSImage writes "mz_axis" with
    // chunked + deflate=6; RamanImage / IRImage write "wavenumbers"
    // with default contiguous layout — match each modality's
    // production writer so the bytes round-trip identically.
    if (axis.length == bins * sizeof(double)) {
        hsize_t aDims[1] = { (hsize_t)bins };
        hid_t aSpace = H5Screate_simple(1, aDims, NULL);
        BOOL isMs = (strcmp(groupName, "image_cube") == 0);
        const char *axisName = isMs ? "mz_axis" : "wavenumbers";
        hid_t aPlist = H5P_DEFAULT;
        hid_t axisPlist = -1;
        if (isMs) {
            axisPlist = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(axisPlist, 1, aDims);
            H5Pset_shuffle(axisPlist);
            H5Pset_deflate(axisPlist, 6);
            aPlist = axisPlist;
        }
        hid_t aDid = H5Dcreate2(g, axisName,
                                 H5T_NATIVE_DOUBLE, aSpace,
                                 H5P_DEFAULT, aPlist, H5P_DEFAULT);
        if (aDid >= 0) {
            H5Dwrite(aDid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, axis.bytes);
            H5Dclose(aDid);
        }
        if (axisPlist >= 0) H5Pclose(axisPlist);
        H5Sclose(aSpace);
    }
    hid_t scalar = H5Screate(H5S_SCALAR);
    #define WRITE_INT_SHARED(name, val) do { \
        hid_t a = H5Acreate2(g, (name), H5T_NATIVE_INT64, \
                              scalar, H5P_DEFAULT, H5P_DEFAULT); \
        int64_t v = (int64_t)(val); H5Awrite(a, H5T_NATIVE_INT64, &v); H5Aclose(a); \
    } while (0)
    #define WRITE_DBL_SHARED(name, val) do { \
        hid_t a = H5Acreate2(g, (name), H5T_NATIVE_DOUBLE, \
                              scalar, H5P_DEFAULT, H5P_DEFAULT); \
        double v = (val); H5Awrite(a, H5T_NATIVE_DOUBLE, &v); H5Aclose(a); \
    } while (0)
    #define WRITE_STR_SHARED(name, val) do { \
        hid_t t = H5Tcopy(H5T_C_S1); H5Tset_size(t, H5T_VARIABLE); \
        hid_t a = H5Acreate2(g, (name), t, scalar, H5P_DEFAULT, H5P_DEFAULT); \
        const char *cs = [(val ?: @"") UTF8String]; H5Awrite(a, t, &cs); \
        H5Aclose(a); H5Tclose(t); \
    } while (0)

    WRITE_INT_SHARED("width",           (int64_t)width);
    WRITE_INT_SHARED("height",          (int64_t)height);
    WRITE_INT_SHARED("spectral_points", (int64_t)bins);
    WRITE_INT_SHARED("tile_size",       (int64_t)tileSize);
    WRITE_DBL_SHARED("pixel_size_x",    pixelSizeX);
    WRITE_DBL_SHARED("pixel_size_y",    pixelSizeY);
    for (NSString *k in doubleAttrs) {
        WRITE_DBL_SHARED([k UTF8String], [doubleAttrs[k] doubleValue]);
    }
    for (NSString *k in intAttrs) {
        WRITE_INT_SHARED([k UTF8String], [intAttrs[k] longLongValue]);
    }
    for (NSString *k in stringAttrs) {
        WRITE_STR_SHARED([k UTF8String], stringAttrs[k]);
    }
    WRITE_STR_SHARED("scan_pattern", (scanPattern ?: @""));

    #undef WRITE_INT_SHARED
    #undef WRITE_DBL_SHARED
    #undef WRITE_STR_SHARED
    H5Sclose(scalar);
    H5Gclose(g);
    return YES;
}

// ---------------------------------------------------------------- materialize

static TTIOPolarity polarityFromWire(uint8_t w)
{
    switch (w) {
        case 0: return TTIOPolarityPositive;
        case 1: return TTIOPolarityNegative;
        default: return TTIOPolarityUnknown;
    }
}

typedef struct {
    uint16_t datasetId;
    NSString *name;
    uint8_t acquisitionMode;
    NSString *spectrumClass;
    NSArray<NSString *> *channelNames;
    uint32_t expectedAUCount;
} DatasetMetaStruct;

- (BOOL)writeTtioToPath:(NSString *)outputPath error:(NSError **)error
{
    NSArray<TTIOTransportPacketRecord *> *packets =
        [self readAllPacketsWithError:error];
    if (!packets) return NO;

    NSString *title = @"";
    NSString *isa = @"";

    NSMutableDictionary<NSNumber *, NSDictionary *> *datasetMetas =
        [NSMutableDictionary dictionary];
    NSMutableDictionary<NSNumber *, NSMutableDictionary *> *runData =
        [NSMutableDictionary dictionary];
    // M89.2 / M89.4: genomic accumulators, keyed by dataset_id. Each
    // value holds the parallel arrays that ultimately feed
    // TTIOWrittenGenomicRun.
    NSMutableDictionary<NSNumber *, NSMutableDictionary *> *genomicData =
        [NSMutableDictionary dictionary];
    // Phase 2c-T: per-dataset_id verbatim v2 blob accumulators.
    NSMutableDictionary<NSNumber *, TTIOBulkV2Blobs *> *bulkBlobs =
        [NSMutableDictionary dictionary];
    NSMutableDictionary<NSNumber *, NSNumber *> *lastSeq =
        [NSMutableDictionary dictionary];
    BOOL sawStreamHeader = NO;
    BOOL bulkModeRequired = NO;

    // v0.11 Task 3.3: per-stream accumulator state for the
    // REFERENCE_GROUP_HEADER -> N x REFERENCE_CHROMOSOME ->
    // END_OF_REFERENCE_GROUP packet sequence (transport-spec
    // §4.13-§4.15). Mirrors Java
    // TransportReader.currentRefUri/currentChromNames/currentChromSeqs/
    // collectedRefs (commit 7f3dec46) and Python's reader (commit
    // 415fc24f). `collectedRefs` is consumed after the dataset is
    // materialised — see the post-loop block further down.
    NSString *currentRefUri = nil;
    NSMutableArray<NSString *> *currentChromNames = [NSMutableArray array];
    NSMutableArray<NSData *> *currentChromSeqs = [NSMutableArray array];
    NSMutableArray<TTIOReferenceImport *> *collectedRefs =
        [NSMutableArray array];

    // v0.11 Task 3.4: per-stream accumulator for the
    // ENCRYPTION_ALGORITHM (0x1B) packet. Consumed after materialise
    // by writing back as the root @encrypted attribute on the .tio
    // (mirrors how a freshly-encrypted dataset surfaces -isEncrypted /
    // -encryptedAlgorithm).
    NSString *collectedEncryptionAlgorithm = nil;

    // v0.11 Task 3.5: per-stream accumulator for the
    // DATASET_PROVENANCE (0x18) packet. Records flow through the
    // matching writeMinimalToPath: parameter so the materialised
    // .tio carries them in /study/provenance just like a source .tio.
    NSMutableArray<TTIOProvenanceRecord *> *collectedProvenance =
        [NSMutableArray array];

    // v0.11 Task 3.6 / 5.1 (Deferral 1): per-stream accumulator for the
    // IMAGE_HEADER -> N x IMAGE_PIXEL -> END_OF_IMAGE packet sequence
    // (transport-spec §4.16-§4.18). Mirrors Java
    // TransportReader.currentImageBuilder / collectedImage (commit
    // a6b1e5d9 + 1889343e) and Python's reader (commit 1f619ced +
    // 8eac605a). The MSImage is finalised on END_OF_IMAGE and embedded
    // into the materialised .tio after writeMinimalToPath: returns.
    // Both continuous-mode (is_continuous == 1) and processed-mode
    // (is_continuous == 0, sparse {channel,intensity} pairs) are
    // supported as of Task 5.1; the mode flag is cached in
    // imgIsContinuous and read by the IMAGE_PIXEL branch.
    uint32_t imgWidth = 0, imgHeight = 0, imgBins = 0;
    double imgPxX = 0.0, imgPxY = 0.0;
    NSMutableString *imgScanPattern = [NSMutableString string];
    NSMutableData *imgMzAxis = [NSMutableData data];
    NSMutableString *imgTitle = [NSMutableString string];
    NSMutableString *imgIsa = [NSMutableString string];
    NSMutableData *imgCube = [NSMutableData data];
    NSMutableData *imgSeen = [NSMutableData data];
    uint64_t imgSeenCount = 0;
    BOOL imgHeaderSeen = NO;
    BOOL imgCollected = NO;
    BOOL imgIsContinuous = YES;
    // v0.11 Stage 5.3: per-modality dispatch state. The IMAGE_HEADER
    // modality byte selects which materialiser writes the cube on
    // END_OF_IMAGE; the modality_extras tail carries the per-modality
    // metadata. imgSkipping is YES between IMAGE_HEADER (unknown
    // modality) and the matching END_OF_IMAGE so following
    // IMAGE_PIXEL packets are silently dropped (forward-compat per
    // §4.16). Per-modality snapshots are taken at END_OF_IMAGE time
    // so MS → Raman → IR on the same stream can all be materialised.
    // Java parity: TransportReader (commit f99ec47d). Python parity:
    // TransportReader (commit 6abead73).
    uint8_t imgModality = 0;
    BOOL imgRamanCollected = NO;
    BOOL imgIRCollected = NO;
    double imgRamanExcitationNm = 0.0;
    double imgRamanLaserPowerMw = 0.0;
    uint8_t imgIRModeByte = 0;
    double imgIRResolutionCmInv = 0.0;
    BOOL imgSkipping = NO;
    // Snapshots (one per modality) — set when END_OF_IMAGE fires so
    // subsequent IMAGE_HEADERs of other modalities don't clobber the
    // earlier modality's bytes.
    uint32_t msImgWidth = 0, msImgHeight = 0, msImgBins = 0;
    double msImgPxX = 0.0, msImgPxY = 0.0;
    NSData *msImgCubeSnap = nil;
    NSData *msImgAxisSnap = nil;
    NSString *msImgScanSnap = nil;
    uint32_t ramanImgWidth = 0, ramanImgHeight = 0, ramanImgBins = 0;
    double ramanImgPxX = 0.0, ramanImgPxY = 0.0;
    NSData *ramanImgCubeSnap = nil;
    NSData *ramanImgAxisSnap = nil;
    NSString *ramanImgScanSnap = nil;
    uint32_t irImgWidth = 0, irImgHeight = 0, irImgBins = 0;
    double irImgPxX = 0.0, irImgPxY = 0.0;
    NSData *irImgCubeSnap = nil;
    NSData *irImgAxisSnap = nil;
    NSString *irImgScanSnap = nil;

    // v0.11 Task 3.7: per-stream accumulators for the
    // IDENTIFICATIONS_TABLE (0x16) and QUANTIFICATIONS_TABLE (0x17)
    // packets (transport-spec §4.19 / §4.20). Multiple 0x16 / 0x17
    // packets MAY appear in a stream (spec §5.4 step 6 says "zero or
    // more"); rows accumulate in emission order. Passed into
    // writeMinimalToPath: as the identifications / quantifications
    // arguments so the on-disk study compound datasets round-trip.
    // Java parity: TransportReader.collectedIdentifications /
    // collectedQuantifications (commit a6faab16). Python parity:
    // TransportReader.materialize_to (commit 150552b6).
    NSMutableArray<TTIOIdentification *> *collectedIdentifications =
        [NSMutableArray array];
    NSMutableArray<TTIOQuantification *> *collectedQuantifications =
        [NSMutableArray array];
    // v0.11 Task 6.4 (Stage 6): per-stream accumulators for
    // SUBJECT_METADATA (0x19) + SAMPLE_METADATA (0x1A) — same shape
    // as the 0x16 / 0x17 accumulators. Materialised onto the
    // generated .tio after writeMinimalToPath: returns (the writer's
    // public API doesn't accept subjects/samples directly, mirroring
    // how Java's TransportReader.materializeTo layers them on with
    // setAttribute after the file is created).
    NSMutableArray<TTIOSubject *> *collectedSubjects =
        [NSMutableArray array];
    NSMutableArray<TTIOSample *> *collectedSamples =
        [NSMutableArray array];

    for (TTIOTransportPacketRecord *rec in packets) {
        TTIOTransportPacketHeader *h = rec.header;
        const uint8_t *bytes = (const uint8_t *)rec.payload.bytes;
        NSUInteger len = rec.payload.length;
        NSUInteger off = 0;

        // Forward-compat (transport-spec §6, v0.11 task 0.7): silently
        // skip packets whose wire type byte was not a known
        // TTIOTransportPacketType. The per-packet loop already logged
        // + length-prefix-consumed the bytes; we just ignore the
        // record here so unknown packets don't trigger
        // MissingStreamHeader / UnexpectedPayload errors. Java parity:
        // `if (h.packetType == null) continue;`.
        if (!TTIOTransportIsKnownPacketType(h.packetTypeByte)) continue;

        if (h.packetType == TTIOTransportPacketStreamHeader) {
            if (sawStreamHeader) continue;
            sawStreamHeader = YES;
            NSString *formatVersion = readLEString(bytes, len, &off, 2); (void)formatVersion;
            title = readLEString(bytes, len, &off, 2) ?: @"";
            isa = readLEString(bytes, len, &off, 2) ?: @"";
            if (off + 2 > len) break;
            uint16_t nFeatures = readU16(&bytes[off]); off += 2;
            for (uint16_t i = 0; i < nFeatures; i++) {
                NSString *f = readLEString(bytes, len, &off, 2);
                if ([f isEqualToString:TTIOTransportBulkModeV2BlobsFeature]) {
                    bulkModeRequired = YES;
                }
            }
            // n_datasets (not needed on the read side)
            continue;
        }
        if (!sawStreamHeader) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorMissingStreamHeader
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"first packet must be StreamHeader"}];
            return NO;
        }

        if (h.packetType == TTIOTransportPacketDatasetHeader) {
            if (off + 2 > len) continue;
            uint16_t did = readU16(&bytes[off]); off += 2;
            NSString *name = readLEString(bytes, len, &off, 2);
            if (off + 1 > len) continue;
            uint8_t acqMode = bytes[off]; off += 1;
            NSString *spectrumClass = readLEString(bytes, len, &off, 2);
            if (off + 1 > len) continue;
            uint8_t nch = bytes[off]; off += 1;
            NSMutableArray<NSString *> *chNames = [NSMutableArray array];
            for (uint8_t i = 0; i < nch; i++) {
                NSString *c = readLEString(bytes, len, &off, 2);
                if (c) [chNames addObject:c];
            }
            (void)readLEString(bytes, len, &off, 4);  // instrument_json
            // expected_au_count
            uint32_t expected = 0;
            if (off + 4 <= len) { expected = readU32(&bytes[off]); off += 4; }

            // Re-read instrument_json (M89.2: genomic dataset header
            // carries reference_uri / platform / sample_name / modality
            // in this slot; we already advanced past it above with a
            // discarded value, so re-extract by walking the payload
            // again before the n_channels byte). Cheaper to keep a
            // local copy from the first pass.
            //
            // (We deliberately skip rewinding — readLEString returned
            // the JSON string we discarded with `(void)`. To minimise
            // churn we recompute by re-reading from a fresh offset.)
            NSUInteger jsonOff = 0;
            jsonOff += 2;  // dataset_id
            // skip name
            (void)readLEString(bytes, len, &jsonOff, 2);
            jsonOff += 1;  // acq_mode
            (void)readLEString(bytes, len, &jsonOff, 2);  // spectrum_class
            jsonOff += 1;  // n_channels
            for (uint8_t i = 0; i < nch; i++) {
                (void)readLEString(bytes, len, &jsonOff, 2);
            }
            NSString *instrumentJSON = readLEString(bytes, len, &jsonOff, 4) ?: @"";

            datasetMetas[@(did)] = @{
                @"name": name ?: @"",
                @"acquisitionMode": @(acqMode),
                @"spectrumClass": spectrumClass ?: @"TTIOMassSpectrum",
                @"channelNames": [chNames copy],
                @"expectedAUCount": @(expected),
                @"instrumentJSON": instrumentJSON,
            };

            // route genomic datasets to a parallel accumulator.
            if ([spectrumClass isEqualToString:@"TTIOGenomicRead"]) {
                NSMutableDictionary *gd = [NSMutableDictionary dictionary];
                gd[@"runningOffset"] = @(0);
                gd[@"chromosomes"] = [NSMutableArray array];
                gd[@"positions"] = [NSMutableArray array];
                gd[@"mappingQualities"] = [NSMutableArray array];
                gd[@"flags"] = [NSMutableArray array];
                gd[@"sequences"] = [NSMutableData data];
                gd[@"qualities"] = [NSMutableData data];
                gd[@"offsets"] = [NSMutableArray array];
                gd[@"lengths"] = [NSMutableArray array];
                // M90.9 compound-field accumulators.
                gd[@"cigars"] = [NSMutableArray array];
                gd[@"readNames"] = [NSMutableArray array];
                gd[@"mateChromosomes"] = [NSMutableArray array];
                // M101: per-read SAM tag text; anyTags marks a stream
                // that carried the tags channel at all.
                gd[@"tags"] = [NSMutableArray array];
                gd[@"anyTags"] = @NO;
                gd[@"matePositions"] = [NSMutableArray array];
                gd[@"templateLengths"] = [NSMutableArray array];
                genomicData[@(did)] = gd;
                continue;
            }

            NSMutableDictionary *rd = [NSMutableDictionary dictionary];
            rd[@"runningOffset"] = @(0);
            rd[@"offsets"] = [NSMutableArray array];
            rd[@"lengths"] = [NSMutableArray array];
            rd[@"retentionTimes"] = [NSMutableArray array];
            rd[@"msLevels"] = [NSMutableArray array];
            rd[@"polarities"] = [NSMutableArray array];
            rd[@"precursorMzs"] = [NSMutableArray array];
            rd[@"precursorCharges"] = [NSMutableArray array];
            rd[@"basePeakIntensities"] = [NSMutableArray array];
            NSMutableDictionary *chans = [NSMutableDictionary dictionary];
            for (NSString *c in chNames) chans[c] = [NSMutableData data];
            rd[@"channels"] = chans;
            runData[@(did)] = rd;
            continue;
        }

        if (h.packetType == TTIOTransportPacketAccessUnit) {
            NSNumber *didKey = @(h.datasetId);
            NSDictionary *meta = datasetMetas[didKey];
            if (!meta) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:@"AccessUnit before DatasetHeader for id %u",
                                         (unsigned)h.datasetId]}];
                return NO;
            }
            NSNumber *prev = lastSeq[didKey];
            if (prev && h.auSequence <= prev.unsignedIntValue) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorNonMonotonicAU
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:@"non-monotonic au_sequence in dataset %u",
                                         (unsigned)h.datasetId]}];
                return NO;
            }
            lastSeq[didKey] = @(h.auSequence);

            NSError *auErr = nil;
            TTIOAccessUnit *au =
                [TTIOAccessUnit decodeFromBytes:bytes length:len error:&auErr];
            if (!au) {
                if (error) *error = auErr;
                return NO;
            }

            // route to genomic accumulator if this dataset is
            // a TTIOGenomicRead stream.
            NSMutableDictionary *gd = genomicData[didKey];
            if (gd) {
                if (au.spectrumClass != 5) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:@"genomic accumulator received spectrum_class %u",
                                             (unsigned)au.spectrumClass]}];
                    return NO;
                }
                [(NSMutableArray *)gd[@"chromosomes"] addObject:(au.chromosome ?: @"")];
                [(NSMutableArray *)gd[@"positions"] addObject:@(au.position)];
                [(NSMutableArray *)gd[@"mappingQualities"] addObject:@(au.mappingQuality)];
                [(NSMutableArray *)gd[@"flags"] addObject:@((uint32_t)au.flags)];
                // AU mate extension — pulled directly off the
                // decoded AU, defaults to -1 / 0 for M89.1 fixtures.
                [(NSMutableArray *)gd[@"matePositions"] addObject:@(au.matePosition)];
                [(NSMutableArray *)gd[@"templateLengths"] addObject:@(au.templateLength)];
                NSMutableData *seqSink = gd[@"sequences"];
                NSMutableData *qualSink = gd[@"qualities"];
                NSUInteger length = 0;
                // M90.9 compound-field defaults — empty when the AU
                // omits the channel (M89.2-era stream).
                NSString *cigarStr = @"";
                NSString *readNameStr = @"";
                NSString *mateChrStr = @"";
                NSString *tagsStr = @"";
                for (TTIOTransportChannelData *ch in au.channels) {
                    if (ch.precision != TTIOPrecisionUInt8) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:@"genomic channel precision %u not supported (UINT8 only)",
                                                 (unsigned)ch.precision]}];
                        return NO;
                    }
                    // dispatch on the wire compression byte.
                    // NONE → identity; RANS_ORDER0/1 → TTIORansDecode;
                    // BASE_PACK → TTIOBasePackDecode. Other codecs
                    // unsupported on the genomic transport path.
                    NSData *decoded = ch.data;
                    if (ch.compression != TTIOCompressionNone) {
                        NSError *decErr = nil;
                        NSData *out = nil;
                        if (ch.compression == TTIOCompressionRansOrder0
                            || ch.compression == TTIOCompressionRansOrder1) {
                            out = TTIORansDecode(ch.data, &decErr);
                        } else if (ch.compression == TTIOCompressionBasePack) {
                            out = TTIOBasePackDecode(ch.data, &decErr);
                        } else {
                            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                     code:TTIOTransportErrorUnexpectedPayload
                                                                 userInfo:@{NSLocalizedDescriptionKey:
                                                 [NSString stringWithFormat:@"genomic channel compression %u unsupported on transport (M90.10)",
                                                     (unsigned)ch.compression]}];
                            return NO;
                        }
                        if (!out) {
                            if (error) *error = decErr ?: [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                                code:TTIOTransportErrorUnexpectedPayload
                                                                            userInfo:@{NSLocalizedDescriptionKey:
                                                  [NSString stringWithFormat:@"genomic channel '%@' codec decode failed", ch.name]}];
                            return NO;
                        }
                        decoded = out;
                    }
                    if ([ch.name isEqualToString:@"sequences"]) {
                        [seqSink appendData:decoded];
                        length = decoded.length;
                    } else if ([ch.name isEqualToString:@"qualities"]) {
                        [qualSink appendData:decoded];
                        if (length == 0) length = decoded.length;
                    } else if ([ch.name isEqualToString:@"cigar"]) {
                        cigarStr = [[NSString alloc] initWithData:decoded
                                                          encoding:NSUTF8StringEncoding] ?: @"";
                    } else if ([ch.name isEqualToString:@"read_name"]) {
                        readNameStr = [[NSString alloc] initWithData:decoded
                                                              encoding:NSUTF8StringEncoding] ?: @"";
                    } else if ([ch.name isEqualToString:@"mate_chromosome"]) {
                        mateChrStr = [[NSString alloc] initWithData:decoded
                                                              encoding:NSUTF8StringEncoding] ?: @"";
                    } else if ([ch.name isEqualToString:@"tags"]) {
                        tagsStr = [[NSString alloc] initWithData:decoded
                                                         encoding:NSUTF8StringEncoding] ?: @"";
                        gd[@"anyTags"] = @YES;
                    }
                }
                [(NSMutableArray *)gd[@"tags"] addObject:tagsStr];
                [(NSMutableArray *)gd[@"cigars"] addObject:cigarStr];
                [(NSMutableArray *)gd[@"readNames"] addObject:readNameStr];
                [(NSMutableArray *)gd[@"mateChromosomes"] addObject:mateChrStr];
                uint64_t curOffset = ((NSNumber *)gd[@"runningOffset"]).unsignedLongLongValue;
                [(NSMutableArray *)gd[@"offsets"] addObject:@(curOffset)];
                [(NSMutableArray *)gd[@"lengths"] addObject:@((uint32_t)length)];
                gd[@"runningOffset"] = @(curOffset + length);
                continue;
            }

            NSMutableDictionary *rd = runData[didKey];
            NSMutableDictionary<NSString *, NSMutableData *> *chans = rd[@"channels"];
            NSUInteger spectrumLength = 0;
            for (TTIOTransportChannelData *ch in au.channels) {
                if (ch.precision != TTIOPrecisionFloat64) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         @"reader supports FLOAT64 precision only"}];
                    return NO;
                }
                NSData *decoded = ch.data;
                if (ch.compression == TTIOCompressionZlib) {
                    // Allocate a generous output buffer. For
                    // float64 payloads the decompressed size is
                    // ch.nElements * 8 exactly.
                    NSMutableData *out = [NSMutableData dataWithLength:(NSUInteger)ch.nElements * 8];
                    uLongf destLen = out.length;
                    int rc = uncompress((Bytef *)out.mutableBytes, &destLen,
                                         (const Bytef *)ch.data.bytes, (uLong)ch.data.length);
                    if (rc != Z_OK) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:@"zlib inflate failed: rc=%d",
                                                 rc]}];
                        return NO;
                    }
                    [out setLength:destLen];
                    decoded = out;
                } else if (ch.compression == TTIOCompressionZstd) {
                    // Exact-size output: the channel header carries
                    // the plaintext element count.
                    NSMutableData *out = [NSMutableData dataWithLength:(NSUInteger)ch.nElements * 8];
                    size_t n = ZSTD_decompress(out.mutableBytes, out.length,
                                               ch.data.bytes, ch.data.length);
                    if (ZSTD_isError(n) || n != out.length) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:@"zstd inflate failed or short: %s",
                                                 ZSTD_isError(n) ? ZSTD_getErrorName(n) : "size mismatch"]}];
                        return NO;
                    }
                    decoded = out;
                } else if (ch.compression == TTIOCompressionFloatDeltaZstd) {
                    NSError *fdzErr = nil;
                    NSData *out = [TTIOFloatDeltaZstd decodeStream:ch.data error:&fdzErr];
                    if (!out) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:@"FLOAT_DELTA_ZSTD decode failed: %@",
                                                 fdzErr.localizedDescription ?: @"unknown"]}];
                        return NO;
                    }
                    if (out.length != (NSUInteger)ch.nElements * 8) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:@"FLOAT_DELTA_ZSTD channel decoded "
                                                 @"%lu values but the channel header declares n_elements=%u",
                                                 (unsigned long)(out.length / 8), (unsigned)ch.nElements]}];
                        return NO;
                    }
                    decoded = out;
                } else if (ch.compression != TTIOCompressionNone) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:@"unsupported compression on reader: %u",
                                             (unsigned)ch.compression]}];
                    return NO;
                }
                NSMutableData *sink = chans[ch.name];
                if (!sink) {
                    sink = [NSMutableData data];
                    chans[ch.name] = sink;
                }
                [sink appendData:decoded];
                NSUInteger n = decoded.length / 8;
                if (spectrumLength != 0 && spectrumLength != n) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         @"channels have mismatched lengths in AU"}];
                    return NO;
                }
                spectrumLength = n;
            }

            uint64_t curOffset = ((NSNumber *)rd[@"runningOffset"]).unsignedLongLongValue;
            [(NSMutableArray *)rd[@"offsets"] addObject:@(curOffset)];
            [(NSMutableArray *)rd[@"lengths"] addObject:@((uint32_t)spectrumLength)];
            rd[@"runningOffset"] = @(curOffset + spectrumLength);
            [(NSMutableArray *)rd[@"retentionTimes"] addObject:@(au.retentionTime)];
            [(NSMutableArray *)rd[@"msLevels"] addObject:@((int32_t)au.msLevel)];
            [(NSMutableArray *)rd[@"polarities"]
                addObject:@((int32_t)polarityFromWire(au.polarity))];
            [(NSMutableArray *)rd[@"precursorMzs"] addObject:@(au.precursorMz)];
            [(NSMutableArray *)rd[@"precursorCharges"]
                addObject:@((int32_t)au.precursorCharge)];
            [(NSMutableArray *)rd[@"basePeakIntensities"]
                addObject:@(au.basePeakIntensity)];
            continue;
        }

        if (h.packetType == TTIOTransportPacketBlobV2MateInfo) {
            if (off + 5 > len) continue;
            uint16_t did = readU16(&bytes[off]); off += 2;
            uint8_t codecId = bytes[off]; off += 1;
            (void)codecId;
            uint16_t nNames = readU16(&bytes[off]); off += 2;
            NSMutableArray<NSString *> *names = [NSMutableArray array];
            for (uint16_t i = 0; i < nNames; i++) {
                NSString *n = readLEString(bytes, len, &off, 2);
                if (n) [names addObject:n];
            }
            if (off + 4 > len) continue;
            uint32_t blobLen = readU32(&bytes[off]); off += 4;
            if (off + blobLen > len) continue;
            NSData *blob = [NSData dataWithBytes:&bytes[off] length:blobLen];
            TTIOBulkV2Blobs *slot = bulkBlobs[@(did)];
            if (!slot) { slot = [TTIOBulkV2Blobs new]; bulkBlobs[@(did)] = slot; }
            slot.mateInfoBlob = blob;
            slot.mateInfoChromNames = names;
            continue;
        }
        if (h.packetType == TTIOTransportPacketBlobV2RefDiff) {
            if (off + 3 > len) continue;
            uint16_t did = readU16(&bytes[off]); off += 2;
            (void)bytes[off]; off += 1;  // codec_id
            NSString *refUri = readLEString(bytes, len, &off, 2);
            if (off + 4 > len) continue;
            uint32_t blobLen = readU32(&bytes[off]); off += 4;
            if (off + blobLen > len) continue;
            NSData *blob = [NSData dataWithBytes:&bytes[off] length:blobLen];
            TTIOBulkV2Blobs *slot = bulkBlobs[@(did)];
            if (!slot) { slot = [TTIOBulkV2Blobs new]; bulkBlobs[@(did)] = slot; }
            slot.refDiffBlob = blob;
            slot.refDiffReferenceUri = refUri ?: @"";
            continue;
        }
        if (h.packetType == TTIOTransportPacketBlobV2NameTok) {
            if (off + 7 > len) continue;
            uint16_t did = readU16(&bytes[off]); off += 2;
            (void)bytes[off]; off += 1;  // codec_id
            uint32_t blobLen = readU32(&bytes[off]); off += 4;
            if (off + blobLen > len) continue;
            NSData *blob = [NSData dataWithBytes:&bytes[off] length:blobLen];
            TTIOBulkV2Blobs *slot = bulkBlobs[@(did)];
            if (!slot) { slot = [TTIOBulkV2Blobs new]; bulkBlobs[@(did)] = slot; }
            slot.nameTokBlob = blob;
            continue;
        }
        // v0.11 Task 3.3: REFERENCE_GROUP_HEADER (0x10) — prime
        // per-group accumulator. chromosome_count / total_bases /
        // md5_hex are parsed for buffer-advance only; the actual
        // values come from the per-chromosome accumulator
        // (TTIOReferenceImport recomputes MD5 from chromosome bytes).
        if (h.packetType == TTIOTransportPacketReferenceGroupHeader) {
            if (off + 2 > len) continue;
            uint16_t uriLen = readU16(&bytes[off]); off += 2;
            if (off + uriLen > len) continue;
            currentRefUri = [[NSString alloc] initWithBytes:&bytes[off]
                                                       length:uriLen
                                                     encoding:NSUTF8StringEncoding]
                ?: @"";
            off += uriLen;
            // chromosome_count (uint32) — advance only.
            if (off + 4 > len) continue;
            off += 4;
            // total_bases (uint64) — advance only.
            if (off + 8 > len) continue;
            off += 8;
            // md5_hex[32] — advance only.
            if (off + 32 > len) continue;
            off += 32;
            [currentChromNames removeAllObjects];
            [currentChromSeqs removeAllObjects];
            continue;
        }
        if (h.packetType == TTIOTransportPacketReferenceChromosome) {
            if (off + 2 > len) continue;
            uint16_t nameLen = readU16(&bytes[off]); off += 2;
            if (off + nameLen > len) continue;
            NSString *name = [[NSString alloc] initWithBytes:&bytes[off]
                                                       length:nameLen
                                                     encoding:NSUTF8StringEncoding]
                ?: @"";
            off += nameLen;
            if (off + 8 > len) continue;
            uint64_t seqLen = readU64(&bytes[off]); off += 8;
            if (off + 1 > len) continue;
            uint8_t encoding = bytes[off]; off += 1;
            if (off + 4 > len) continue;
            uint32_t dataLen = readU32(&bytes[off]); off += 4;
            if (off + dataLen > len) continue;
            NSData *seqBytes = nil;
            if (encoding == 0) {
                seqBytes = [NSData dataWithBytes:&bytes[off] length:dataLen];
                if ((uint64_t)seqBytes.length != seqLen) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:
                                             @"REFERENCE_CHROMOSOME raw payload length %u "
                                             @"!= sequence_length %llu",
                                             (unsigned)dataLen,
                                             (unsigned long long)seqLen]}];
                    return NO;
                }
            } else if (encoding == 1) {
                NSData *deflated = [NSData dataWithBytesNoCopy:(void *)&bytes[off]
                                                          length:dataLen
                                                    freeWhenDone:NO];
                NSError *infErr = nil;
                seqBytes = zlibInflateExact(deflated,
                                               (NSUInteger)seqLen,
                                               &infErr);
                if (!seqBytes) {
                    if (error) *error = infErr;
                    return NO;
                }
            } else {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"unknown REFERENCE_CHROMOSOME encoding: %u",
                                         (unsigned)encoding]}];
                return NO;
            }
            [currentChromNames addObject:name];
            [currentChromSeqs addObject:seqBytes];
            continue;
        }
        if (h.packetType == TTIOTransportPacketEndOfReferenceGroup) {
            if (currentRefUri == nil) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     @"END_OF_REFERENCE_GROUP without prior "
                                     @"REFERENCE_GROUP_HEADER"}];
                return NO;
            }
            TTIOReferenceImport *ref =
                [[TTIOReferenceImport alloc]
                    initWithUri:currentRefUri
                    chromosomes:[currentChromNames copy]
                      sequences:[currentChromSeqs copy]];
            [collectedRefs addObject:ref];
            currentRefUri = nil;
            [currentChromNames removeAllObjects];
            [currentChromSeqs removeAllObjects];
            continue;
        }
        // v0.11 Task 3.4: ENCRYPTION_ALGORITHM (0x1B). Wire layout per
        // transport-spec §4.23: `uint16 algorithm_length + UTF-8
        // bytes[algorithm_length]`. Stashed for post-materialise
        // attachment as the root @encrypted attribute. Java parity:
        // TransportReader.materializeTo (commit 530a5833) which calls
        // markRootEncryptedWithEncryptionAlgorithm. Python parity:
        // TransportReader.materialize_to (commit bf38bdc9) which sets
        // the @encrypted root attr through the provider.
        if (h.packetType == TTIOTransportPacketEncryptionAlgorithm) {
            if (off + 2 > len) continue;
            uint16_t algoLen = readU16(&bytes[off]); off += 2;
            if (off + algoLen > len) continue;
            collectedEncryptionAlgorithm =
                [[NSString alloc] initWithBytes:&bytes[off]
                                          length:algoLen
                                        encoding:NSUTF8StringEncoding]
                ?: @"";
            off += algoLen;
            continue;
        }
        // v0.11 Task 3.5: DATASET_PROVENANCE (0x18). Wire layout per
        // transport-spec §4.21: `uint32 record_count` then N records,
        // each `int64 timestamp_unix + uint16 software_length + UTF-8
        // bytes + uint16 parameters_length + UTF-8 JSON +
        // uint16 input_refs_length + UTF-8 CSV +
        // uint16 output_refs_length + UTF-8 CSV`. The parameters JSON
        // is decoded into an NSDictionary for the TTIOProvenanceRecord
        // initialiser; CSV refs split on `,` with empty-string -> [].
        // Java parity: TransportReader.materializeTo (commit
        // 563e09c3). Python parity: TransportReader.materialize_to
        // (commit 434d45a6).
        if (h.packetType == TTIOTransportPacketDatasetProvenance) {
            if (off + 4 > len) continue;
            uint32_t recordCount = readU32(&bytes[off]); off += 4;
            for (uint32_t recI = 0; recI < recordCount; recI++) {
                if (off + 8 > len) break;
                uint64_t tsBits = readU64(&bytes[off]); off += 8;
                int64_t timestampUnix = (int64_t)tsBits;
                NSString *software = readLEString(bytes, len, &off, 2) ?: @"";
                NSString *paramsJson = readLEString(bytes, len, &off, 2) ?: @"{}";
                NSString *inputsCsv = readLEString(bytes, len, &off, 2) ?: @"";
                NSString *outputsCsv = readLEString(bytes, len, &off, 2) ?: @"";
                // Decode parameters_json -> NSDictionary. Empty dict
                // if parse fails (treat as `{}`).
                NSDictionary *paramsDict = @{};
                NSData *pjData =
                    [paramsJson dataUsingEncoding:NSUTF8StringEncoding];
                if (pjData.length > 0) {
                    id parsed = [NSJSONSerialization JSONObjectWithData:pjData
                                                                 options:0
                                                                   error:NULL];
                    if ([parsed isKindOfClass:[NSDictionary class]]) {
                        paramsDict = parsed;
                    }
                }
                NSArray<NSString *> *inputRefs = inputsCsv.length > 0
                    ? [inputsCsv componentsSeparatedByString:@","]
                    : @[];
                NSArray<NSString *> *outputRefs = outputsCsv.length > 0
                    ? [outputsCsv componentsSeparatedByString:@","]
                    : @[];
                TTIOProvenanceRecord *rec =
                    [[TTIOProvenanceRecord alloc]
                        initWithInputRefs:inputRefs
                                 software:software
                               parameters:paramsDict
                               outputRefs:outputRefs
                            timestampUnix:timestampUnix];
                [collectedProvenance addObject:rec];
            }
            continue;
        }
        // v0.11 Task 3.6: IMAGE_HEADER (0x13). Wire layout per
        // transport-spec §4.16: u8 modality + u32 width + u32 height
        // + u32 spectrum_bins + f64 pixel_size_x + f64 pixel_size_y
        // + u8 scan_pattern + u8 axis_kind + u32 axis_length +
        // N x f64 axis + u8 is_continuous + u16 + UTF-8 title +
        // u16 + UTF-8 isa_id. Java parity: TransportReader.startImage
        // (commit a6b1e5d9). Python parity: TransportReader._start_image
        // (commit 1f619ced).
        if (h.packetType == TTIOTransportPacketImageHeader) {
            if (off + 1 > len) continue;
            uint8_t modality = bytes[off]; off += 1;
            if (off + 4 > len) continue;
            uint32_t hdrW = readU32(&bytes[off]); off += 4;
            if (off + 4 > len) continue;
            uint32_t hdrH = readU32(&bytes[off]); off += 4;
            if (off + 4 > len) continue;
            uint32_t hdrBins = readU32(&bytes[off]); off += 4;
            if (off + 8 > len) continue;
            double pxX = readF64(&bytes[off]); off += 8;
            if (off + 8 > len) continue;
            double pxY = readF64(&bytes[off]); off += 8;
            if (off + 1 > len) continue;
            uint8_t scanPat = bytes[off]; off += 1;
            if (off + 1 > len) continue;
            uint8_t axisKind = bytes[off]; off += 1;
            (void)axisKind;  // accepted but not currently surfaced
            if (off + 4 > len) continue;
            uint32_t axisLen = readU32(&bytes[off]); off += 4;
            if (off + 8 * axisLen > len) continue;
            [imgMzAxis setLength:0];
            for (uint32_t i = 0; i < axisLen; i++) {
                double v = readF64(&bytes[off]); off += 8;
                [imgMzAxis appendBytes:&v length:sizeof(double)];
            }
            if (off + 1 > len) continue;
            uint8_t isCont = bytes[off]; off += 1;
            NSString *title = readLEString(bytes, len, &off, 2) ?: @"";
            NSString *isa = readLEString(bytes, len, &off, 2) ?: @"";

            // v0.11 Stage 5.3: parse modality_extras tail. Older
            // (pre-5.3) streams emit no such slot, so probe the
            // remaining length: an unfilled slot stays empty rather
            // than fault the stream.
            NSData *imgExtras = [NSData data];
            if (off + 2 <= len) {
                uint16_t extrasLen = readU16(&bytes[off]); off += 2;
                if (off + extrasLen > len) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:
                                             @"IMAGE_HEADER: modality_extras_length "
                                             @"%u exceeds remaining payload",
                                             (unsigned)extrasLen]}];
                    return NO;
                }
                imgExtras = [NSData dataWithBytes:&bytes[off]
                                            length:extrasLen];
                off += extrasLen;
            }

            if (isCont != 0 && isCont != 1) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"IMAGE_HEADER: is_continuous must be "
                                         @"0 or 1; got %u", (unsigned)isCont]}];
                return NO;
            }

            // v0.11 Stage 5.3: modality dispatch — 0=MS, 1=Raman,
            // 2=IR. Unknown modalities are logged + skipped (forward
            // compat per §4.16); the self-describing extras_len has
            // already advanced the buffer past the IMAGE_HEADER.
            imgSkipping = NO;
            if (modality == 0) {
                // MS — no extras consumed.
            } else if (modality == 1) {
                if (imgExtras.length != 16) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:
                                             @"IMAGE_HEADER (modality=1, Raman) "
                                             @"expects 16-byte modality_extras "
                                             @"(excitation + laser_power); got %lu",
                                             (unsigned long)imgExtras.length]}];
                    return NO;
                }
                const uint8_t *eb = imgExtras.bytes;
                imgRamanExcitationNm = readF64(&eb[0]);
                imgRamanLaserPowerMw = readF64(&eb[8]);
            } else if (modality == 2) {
                if (imgExtras.length != 9) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:
                                             @"IMAGE_HEADER (modality=2, IR) "
                                             @"expects 9-byte modality_extras "
                                             @"(ir_mode + resolution); got %lu",
                                             (unsigned long)imgExtras.length]}];
                    return NO;
                }
                const uint8_t *eb = imgExtras.bytes;
                imgIRModeByte = eb[0];
                imgIRResolutionCmInv = readF64(&eb[1]);
            } else {
                NSLog(@"[TTIOTransportReader] IMAGE_HEADER: unknown "
                      @"modality=%u; skipping image block (extras_len=%lu, "
                      @"width=%u, height=%u)",
                      (unsigned)modality, (unsigned long)imgExtras.length,
                      (unsigned)hdrW, (unsigned)hdrH);
                imgSkipping = YES;
                imgHeaderSeen = YES;  // need a matching END_OF_IMAGE to clear
                continue;
            }
            imgModality = modality;
            imgWidth = hdrW;
            imgHeight = hdrH;
            imgBins = hdrBins;
            imgPxX = pxX;
            imgPxY = pxY;
            [imgScanPattern setString:@""];
            switch (scanPat) {
                case 0:  [imgScanPattern setString:@"raster"];  break;
                case 1:  [imgScanPattern setString:@"meander"]; break;
                case 2:  [imgScanPattern setString:@"random"];  break;
                default: [imgScanPattern setString:@"raster"];  break;
            }
            [imgTitle setString:title];
            [imgIsa setString:isa];
            [imgCube setLength:(NSUInteger)hdrW * hdrH * hdrBins * sizeof(double)];
            memset(imgCube.mutableBytes, 0, imgCube.length);
            [imgSeen setLength:(NSUInteger)hdrW * hdrH];
            memset(imgSeen.mutableBytes, 0, imgSeen.length);
            imgSeenCount = 0;
            imgHeaderSeen = YES;
            imgIsContinuous = (isCont == 1);
            continue;
        }
        // v0.11 Task 3.6: IMAGE_PIXEL (0x14). Wire layout per
        // transport-spec §4.17: u32 x + u32 y + u8 precision +
        // u8 compression + u32 payload_length + intensities[..].
        // Continuous-mode only — precision MUST be 0 (float32) or
        // 1 (float64); compression MUST be 0 (NONE).
        if (h.packetType == TTIOTransportPacketImagePixel) {
            // v0.11 Stage 5.3: unknown-modality stream — silently drop
            // the pixel until the matching END_OF_IMAGE.
            if (imgSkipping) continue;
            if (!imgHeaderSeen) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     @"IMAGE_PIXEL received before IMAGE_HEADER"}];
                return NO;
            }
            if (off + 4 > len) continue;
            uint32_t px = readU32(&bytes[off]); off += 4;
            if (off + 4 > len) continue;
            uint32_t py = readU32(&bytes[off]); off += 4;
            if (off + 1 > len) continue;
            uint8_t precision = bytes[off]; off += 1;
            if (off + 1 > len) continue;
            uint8_t compression = bytes[off]; off += 1;
            if (off + 4 > len) continue;
            uint32_t payloadLen = readU32(&bytes[off]); off += 4;
            if (off + payloadLen > len) continue;
            /* §4.17 pixel enum: 0=none, 1=zstd, 2=zlib. Writers emit
               0 today; the inflate paths make the enum real. The
               decoded buffer replaces (pix, pixLen) for the rest of
               this branch. */
            const uint8_t *pix = &bytes[off];
            NSUInteger pixLen = payloadLen;
            NSData *pixInflated = nil;
            if (compression == 1) {
                unsigned long long plain =
                    ZSTD_getFrameContentSize(pix, pixLen);
                if (plain == ZSTD_CONTENTSIZE_ERROR
                    || plain == ZSTD_CONTENTSIZE_UNKNOWN
                    || plain > (1ull << 27)) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         @"IMAGE_PIXEL zstd frame without a sane content size"}];
                    return NO;
                }
                NSMutableData *out = [NSMutableData dataWithLength:(NSUInteger)plain];
                size_t n = ZSTD_decompress(out.mutableBytes, out.length, pix, pixLen);
                if (ZSTD_isError(n) || n != out.length) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         @"IMAGE_PIXEL zstd inflate failed"}];
                    return NO;
                }
                pixInflated = out;
            } else if (compression == 2) {
                pixInflated = zlibInflateGrow(pix, pixLen, error);
                if (!pixInflated) return NO;
            } else if (compression != 0) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"IMAGE_PIXEL compression=%u not "
                                         @"supported (0=none, 1=zstd, 2=zlib)",
                                         (unsigned)compression]}];
                return NO;
            }
            if (pixInflated) {
                pix = pixInflated.bytes;
                pixLen = pixInflated.length;
            }
            if (px >= imgWidth || py >= imgHeight) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"IMAGE_PIXEL coordinates out of bounds: "
                                         @"x=%u, y=%u (width=%u, height=%u)",
                                         (unsigned)px, (unsigned)py,
                                         (unsigned)imgWidth, (unsigned)imgHeight]}];
                return NO;
            }
            NSUInteger pixelIdx = (NSUInteger)py * imgWidth + px;
            uint8_t *seenBytes = imgSeen.mutableBytes;
            if (seenBytes[pixelIdx] != 0) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"duplicate IMAGE_PIXEL at (x=%u, y=%u)",
                                         (unsigned)px, (unsigned)py]}];
                return NO;
            }
            NSUInteger base = pixelIdx * imgBins;
            double *cubeWrite = (double *)imgCube.mutableBytes;
            if (precision != 0 && precision != 1) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"IMAGE_PIXEL precision=%u not supported "
                                         @"(expected 0=float32 or 1=float64)",
                                         (unsigned)precision]}];
                return NO;
            }
            if (imgIsContinuous) {
                // Continuous mode: dense intensity vector
                // (spec §4.17). One f32/f64 per channel; total
                // bytes == bins * sizeof(precision).
                if (precision == 1) {
                    // FLOAT64
                    NSUInteger n = pixLen / 8;
                    if (n != imgBins) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:
                                                 @"IMAGE_PIXEL intensity count %lu does not "
                                                 @"match IMAGE_HEADER.spectrum_bins=%u",
                                                 (unsigned long)n, (unsigned)imgBins]}];
                        return NO;
                    }
                    for (NSUInteger k = 0; k < n; k++) {
                        cubeWrite[base + k] = readF64(&pix[8 * k]);
                    }
                } else {
                    // FLOAT32 — widen into the float64 cube.
                    NSUInteger n = pixLen / 4;
                    if (n != imgBins) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:
                                                 @"IMAGE_PIXEL intensity count %lu does not "
                                                 @"match IMAGE_HEADER.spectrum_bins=%u",
                                                 (unsigned long)n, (unsigned)imgBins]}];
                        return NO;
                    }
                    for (NSUInteger k = 0; k < n; k++) {
                        uint32_t bits = readU32(&pix[4 * k]);
                        float f;
                        memcpy(&f, &bits, 4);
                        cubeWrite[base + k] = (double)f;
                    }
                }
            } else {
                // v0.11 Task 5.1: processed-mode payload —
                //   u32 nonzero_count
                //   nonzero_count × { u32 channel + fXX intensity }.
                // Cube was zero-initialised at IMAGE_HEADER time, so
                // unmentioned channels stay 0.0.
                if (pixLen < 4) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         @"IMAGE_PIXEL (processed): payload < 4 bytes"}];
                    return NO;
                }
                uint32_t nonzero = readU32(pix);
                NSUInteger entrySize = (precision == 1) ? (4 + 8) : (4 + 4);
                if (pixLen != 4 + (NSUInteger)nonzero * entrySize) {
                    if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                             code:TTIOTransportErrorUnexpectedPayload
                                                         userInfo:@{NSLocalizedDescriptionKey:
                                         [NSString stringWithFormat:
                                             @"IMAGE_PIXEL (processed): payload_length "
                                             @"%u does not match 4 + nonzero_count(%u)*"
                                             @"entry_size(%lu)",
                                             (unsigned)payloadLen,
                                             (unsigned)nonzero,
                                             (unsigned long)entrySize]}];
                    return NO;
                }
                NSUInteger entryOff = 4;
                for (uint32_t e = 0; e < nonzero; e++) {
                    uint32_t ch = readU32(&pix[entryOff]);
                    entryOff += 4;
                    if (ch >= imgBins) {
                        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                                 code:TTIOTransportErrorUnexpectedPayload
                                                             userInfo:@{NSLocalizedDescriptionKey:
                                             [NSString stringWithFormat:
                                                 @"IMAGE_PIXEL (processed): channel_index "
                                                 @"%u out of range [0, %u) at pixel "
                                                 @"(x=%u, y=%u)",
                                                 (unsigned)ch, (unsigned)imgBins,
                                                 (unsigned)px, (unsigned)py]}];
                        return NO;
                    }
                    double v;
                    if (precision == 1) {
                        v = readF64(&pix[entryOff]);
                        entryOff += 8;
                    } else {
                        uint32_t bits = readU32(&pix[entryOff]);
                        float f;
                        memcpy(&f, &bits, 4);
                        v = (double)f;
                        entryOff += 4;
                    }
                    cubeWrite[base + ch] = v;
                }
            }
            seenBytes[pixelIdx] = 1;
            imgSeenCount++;
            continue;
        }
        // v0.11 Task 3.6: END_OF_IMAGE (0x15). Wire layout per
        // transport-spec §4.18: u32 pixel_count_seen. Verifies the
        // declared count matches the per-pixel ingest count and stages
        // the built MSImage for write-out after writeMinimalToPath:.
        if (h.packetType == TTIOTransportPacketEndOfImage) {
            // v0.11 Stage 5.3: drain an unknown-modality block.
            if (imgSkipping) {
                imgSkipping = NO;
                imgHeaderSeen = NO;
                continue;
            }
            if (!imgHeaderSeen) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     @"END_OF_IMAGE without prior IMAGE_HEADER"}];
                return NO;
            }
            if (off + 4 > len) continue;
            uint64_t declared = (uint64_t)readU32(&bytes[off]); off += 4;
            uint64_t expected = (uint64_t)imgWidth * (uint64_t)imgHeight;
            if (declared != imgSeenCount) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"END_OF_IMAGE pixel_count_seen mismatch: "
                                         @"declared=%llu, actual=%llu (width*height=%llu)",
                                         (unsigned long long)declared,
                                         (unsigned long long)imgSeenCount,
                                         (unsigned long long)expected]}];
                return NO;
            }
            if (imgSeenCount != expected) {
                if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                         code:TTIOTransportErrorUnexpectedPayload
                                                     userInfo:@{NSLocalizedDescriptionKey:
                                     [NSString stringWithFormat:
                                         @"END_OF_IMAGE pixel count %llu does not "
                                         @"equal width*height=%llu",
                                         (unsigned long long)imgSeenCount,
                                         (unsigned long long)expected]}];
                return NO;
            }
            // Per-modality collection flags + snapshots. The
            // materialiser branches on these after writeMinimalToPath:
            // returns. Snapshots use [imgCube copy] so an immediate
            // following IMAGE_HEADER (of a different modality) can
            // reset the shared imgCube buffer without invalidating
            // the earlier modality's data. ARC keeps the strong
            // references alive via the assignments below.
            NSData *cubeSnap = [imgCube copy];
            NSData *axisSnap = [imgMzAxis copy];
            NSString *scanSnap = [imgScanPattern copy];
            if (imgModality == 0) {
                imgCollected = YES;
                msImgWidth = imgWidth; msImgHeight = imgHeight; msImgBins = imgBins;
                msImgPxX = imgPxX; msImgPxY = imgPxY;
                msImgCubeSnap = cubeSnap;
                msImgAxisSnap = axisSnap;
                msImgScanSnap = scanSnap;
            } else if (imgModality == 1) {
                imgRamanCollected = YES;
                ramanImgWidth = imgWidth; ramanImgHeight = imgHeight; ramanImgBins = imgBins;
                ramanImgPxX = imgPxX; ramanImgPxY = imgPxY;
                ramanImgCubeSnap = cubeSnap;
                ramanImgAxisSnap = axisSnap;
                ramanImgScanSnap = scanSnap;
            } else if (imgModality == 2) {
                imgIRCollected = YES;
                irImgWidth = imgWidth; irImgHeight = imgHeight; irImgBins = imgBins;
                irImgPxX = imgPxX; irImgPxY = imgPxY;
                irImgCubeSnap = cubeSnap;
                irImgAxisSnap = axisSnap;
                irImgScanSnap = scanSnap;
            }
            imgHeaderSeen = NO;
            continue;
        }
        // v0.11 Task 3.7: IDENTIFICATIONS_TABLE (0x16). Wire layout per
        // transport-spec §4.19: `uint32 arrow_ipc_length + bytes
        // arrow_ipc[length]` — a single self-describing Apache Arrow
        // IPC stream. Decoded via TTIOArrowIpcCodec; rows accumulate
        // across multiple 0x16 packets (spec §5.4 step 6 "zero or
        // more"). Java parity: TransportReader.decodeIdentificationsTable
        // (commit a6faab16). Python parity:
        // TransportReader.materialize_to (commit 150552b6).
        if (h.packetType == TTIOTransportPacketIdentificationsTable) {
            if (off + 4 > len) continue;
            uint32_t ipcLen = readU32(&bytes[off]); off += 4;
            if (off + ipcLen > len) continue;
            NSData *ipcBytes = [NSData dataWithBytes:&bytes[off]
                                              length:ipcLen];
            off += ipcLen;
            NSArray<TTIOIdentification *> *decoded =
                [TTIOArrowIpcCodec decodeIdentifications:ipcBytes];
            if (decoded.count > 0) {
                [collectedIdentifications addObjectsFromArray:decoded];
            }
            continue;
        }
        // v0.11 Task 3.7: QUANTIFICATIONS_TABLE (0x17). Wire layout per
        // transport-spec §4.20 — identical shape to §4.19.
        if (h.packetType == TTIOTransportPacketQuantificationsTable) {
            if (off + 4 > len) continue;
            uint32_t ipcLen = readU32(&bytes[off]); off += 4;
            if (off + ipcLen > len) continue;
            NSData *ipcBytes = [NSData dataWithBytes:&bytes[off]
                                              length:ipcLen];
            off += ipcLen;
            NSArray<TTIOQuantification *> *decoded =
                [TTIOArrowIpcCodec decodeQuantifications:ipcBytes];
            if (decoded.count > 0) {
                [collectedQuantifications addObjectsFromArray:decoded];
            }
            continue;
        }
        // v0.11 Task 6.4 (Stage 6): SUBJECT_METADATA (0x19) +
        // SAMPLE_METADATA (0x1A). Wire layout per transport-spec
        // §4.22 — `uint32 arrow_ipc_length + bytes
        // arrow_ipc[length]`. Rows accumulate across multiple packets
        // per spec §5.4 step 5 "zero or more". Java parity:
        // TransportReader.decodeSubjectMetadata /
        // decodeSampleMetadata (commit dd211600). Python parity:
        // TransportReader.materialize_to (commit 00c7e1b7).
        if (h.packetType == TTIOTransportPacketSubjectMetadata) {
            if (off + 4 > len) continue;
            uint32_t ipcLen = readU32(&bytes[off]); off += 4;
            if (off + ipcLen > len) continue;
            NSData *ipcBytes = [NSData dataWithBytes:&bytes[off]
                                              length:ipcLen];
            off += ipcLen;
            NSArray<TTIOSubject *> *decoded =
                [TTIOArrowIpcCodec decodeSubjects:ipcBytes];
            if (decoded.count > 0) {
                [collectedSubjects addObjectsFromArray:decoded];
            }
            continue;
        }
        if (h.packetType == TTIOTransportPacketSampleMetadata) {
            if (off + 4 > len) continue;
            uint32_t ipcLen = readU32(&bytes[off]); off += 4;
            if (off + ipcLen > len) continue;
            NSData *ipcBytes = [NSData dataWithBytes:&bytes[off]
                                              length:ipcLen];
            off += ipcLen;
            NSArray<TTIOSample *> *decoded =
                [TTIOArrowIpcCodec decodeSamples:ipcBytes];
            if (decoded.count > 0) {
                [collectedSamples addObjectsFromArray:decoded];
            }
            continue;
        }
        if (h.packetType == TTIOTransportPacketEndOfDataset) continue;
        if (h.packetType == TTIOTransportPacketEndOfStream) break;
        // Annotation/Provenance/Chromatogram/Protection — skipped in M67.
    }

    // Phase 2c-T: a stream that declared bulk_mode_v2_blobs but
    // shipped zero blob packets is malformed.
    if (bulkModeRequired && bulkBlobs.count == 0) {
        if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                 code:TTIOTransportErrorUnexpectedPayload
                                             userInfo:@{NSLocalizedDescriptionKey:
                                 @"StreamHeader declared bulk_mode_v2_blobs "
                                 @"but no BlobV2* packets were received"}];
        return NO;
    }

    // Build TTIOWrittenRun objects.
    NSMutableDictionary<NSString *, TTIOWrittenRun *> *runs =
        [NSMutableDictionary dictionary];
    for (NSNumber *didKey in datasetMetas) {
        NSDictionary *meta = datasetMetas[didKey];
        // Skip genomic datasets — built separately below.
        if (genomicData[didKey]) continue;
        NSDictionary *rd = runData[didKey];
        NSMutableDictionary *channelDataOut = [NSMutableDictionary dictionary];
        for (NSString *c in (NSArray *)meta[@"channelNames"]) {
            NSMutableData *src = rd[@"channels"][c];
            channelDataOut[c] = src ? [src copy] : [NSData data];
        }

        NSArray *offArr = rd[@"offsets"];
        NSMutableData *offsetsData = [NSMutableData dataWithCapacity:offArr.count * 8];
        for (NSNumber *n in offArr) {
            uint64_t v = n.unsignedLongLongValue;
            [offsetsData appendBytes:&v length:8];
        }
        NSArray *lenArr = rd[@"lengths"];
        NSMutableData *lengthsData = [NSMutableData dataWithCapacity:lenArr.count * 4];
        for (NSNumber *n in lenArr) {
            uint32_t v = n.unsignedIntValue;
            [lengthsData appendBytes:&v length:4];
        }
        NSArray *rtArr = rd[@"retentionTimes"];
        NSMutableData *rtData = [NSMutableData dataWithCapacity:rtArr.count * 8];
        for (NSNumber *n in rtArr) {
            double v = n.doubleValue;
            [rtData appendBytes:&v length:8];
        }
        NSArray *msArr = rd[@"msLevels"];
        NSMutableData *msData = [NSMutableData dataWithCapacity:msArr.count * 4];
        for (NSNumber *n in msArr) {
            int32_t v = n.intValue;
            [msData appendBytes:&v length:4];
        }
        NSArray *polArr = rd[@"polarities"];
        NSMutableData *polData = [NSMutableData dataWithCapacity:polArr.count * 4];
        for (NSNumber *n in polArr) {
            int32_t v = n.intValue;
            [polData appendBytes:&v length:4];
        }
        NSArray *pmzArr = rd[@"precursorMzs"];
        NSMutableData *pmzData = [NSMutableData dataWithCapacity:pmzArr.count * 8];
        for (NSNumber *n in pmzArr) {
            double v = n.doubleValue;
            [pmzData appendBytes:&v length:8];
        }
        NSArray *pcArr = rd[@"precursorCharges"];
        NSMutableData *pcData = [NSMutableData dataWithCapacity:pcArr.count * 4];
        for (NSNumber *n in pcArr) {
            int32_t v = n.intValue;
            [pcData appendBytes:&v length:4];
        }
        NSArray *bpiArr = rd[@"basePeakIntensities"];
        NSMutableData *bpiData = [NSMutableData dataWithCapacity:bpiArr.count * 8];
        for (NSNumber *n in bpiArr) {
            double v = n.doubleValue;
            [bpiData appendBytes:&v length:8];
        }

        TTIOWrittenRun *wr =
            [[TTIOWrittenRun alloc]
                initWithSpectrumClassName:(NSString *)meta[@"spectrumClass"]
                          acquisitionMode:((NSNumber *)meta[@"acquisitionMode"]).longLongValue
                              channelData:channelDataOut
                                  offsets:offsetsData
                                  lengths:lengthsData
                           retentionTimes:rtData
                                 msLevels:msData
                               polarities:polData
                             precursorMzs:pmzData
                         precursorCharges:pcData
                      basePeakIntensities:bpiData];
        runs[(NSString *)meta[@"name"]] = wr;
    }

    // build TTIOWrittenGenomicRun objects for each genomic
    // dataset_id. These travel through the extended writeMinimalToPath
    // overload alongside any MS runs.
    NSMutableDictionary<NSString *, TTIOWrittenGenomicRun *> *genomicRuns =
        [NSMutableDictionary dictionary];
    for (NSNumber *didKey in genomicData) {
        NSDictionary *meta = datasetMetas[didKey];
        NSMutableDictionary *gd = genomicData[didKey];
        NSString *instrumentJSON = meta[@"instrumentJSON"] ?: @"";
        NSString *referenceUri = @"";
        NSString *platform = @"";
        NSString *sampleName = @"";
        NSString *readRole = @"";
        if (instrumentJSON.length > 0) {
            NSData *jdata = [instrumentJSON dataUsingEncoding:NSUTF8StringEncoding];
            id parsed = [NSJSONSerialization JSONObjectWithData:jdata options:0 error:NULL];
            if ([parsed isKindOfClass:[NSDictionary class]]) {
                NSDictionary *jd = parsed;
                referenceUri = [jd[@"reference_uri"] isKindOfClass:[NSString class]]
                    ? jd[@"reference_uri"] : @"";
                platform     = [jd[@"platform"] isKindOfClass:[NSString class]]
                    ? jd[@"platform"]     : @"";
                sampleName   = [jd[@"sample_name"] isKindOfClass:[NSString class]]
                    ? jd[@"sample_name"]   : @"";
                /* M97: "" on the wire means no role. */
                readRole     = [jd[@"read_role"] isKindOfClass:[NSString class]]
                    ? jd[@"read_role"]    : @"";
            }
        }

        NSArray *posArr = gd[@"positions"];
        NSMutableData *positionsData = [NSMutableData dataWithCapacity:posArr.count * 8];
        for (NSNumber *n in posArr) {
            int64_t v = n.longLongValue;
            [positionsData appendBytes:&v length:8];
        }
        NSArray *mqArr = gd[@"mappingQualities"];
        NSMutableData *mqData = [NSMutableData dataWithCapacity:mqArr.count];
        for (NSNumber *n in mqArr) {
            uint8_t v = (uint8_t)n.unsignedCharValue;
            [mqData appendBytes:&v length:1];
        }
        NSArray *flagsArr = gd[@"flags"];
        NSMutableData *flagsData = [NSMutableData dataWithCapacity:flagsArr.count * 4];
        for (NSNumber *n in flagsArr) {
            uint32_t v = (uint32_t)n.unsignedIntValue;
            [flagsData appendBytes:&v length:4];
        }
        NSArray *offArr = gd[@"offsets"];
        NSMutableData *offsetsData = [NSMutableData dataWithCapacity:offArr.count * 8];
        for (NSNumber *n in offArr) {
            uint64_t v = n.unsignedLongLongValue;
            [offsetsData appendBytes:&v length:8];
        }
        NSArray *lenArr = gd[@"lengths"];
        NSMutableData *lengthsData = [NSMutableData dataWithCapacity:lenArr.count * 4];
        for (NSNumber *n in lenArr) {
            uint32_t v = n.unsignedIntValue;
            [lengthsData appendBytes:&v length:4];
        }

        // compound fields ride on the wire as 3 string channels
        // + a 12-byte mate extension on the AU genomic suffix. The
        // accumulator captured them per-AU; materialise into the
        // run-level shapes the WrittenGenomicRun expects. M89.1-only
        // streams default to "" / -1 / 0 because the AU decoder
        // returns those defaults when the extension is absent.
        NSUInteger n = posArr.count;
        NSArray *cigarsCollected = gd[@"cigars"] ?: @[];
        NSArray *readNamesCollected = gd[@"readNames"] ?: @[];
        NSArray *mateChromsCollected = gd[@"mateChromosomes"] ?: @[];
        NSArray *matePositionsCollected = gd[@"matePositions"] ?: @[];
        NSArray *templateLengthsCollected = gd[@"templateLengths"] ?: @[];
        NSMutableArray *cigars = [NSMutableArray arrayWithCapacity:n];
        NSMutableArray *readNames = [NSMutableArray arrayWithCapacity:n];
        NSMutableArray *mateChroms = [NSMutableArray arrayWithCapacity:n];
        for (NSUInteger i = 0; i < n; i++) {
            [cigars addObject:(i < cigarsCollected.count
                                ? cigarsCollected[i] : @"")];
            [readNames addObject:(i < readNamesCollected.count
                                   ? readNamesCollected[i] : @"")];
            [mateChroms addObject:(i < mateChromsCollected.count
                                    ? mateChromsCollected[i] : @"")];
        }
        NSMutableData *matePosData = [NSMutableData dataWithLength:n * sizeof(int64_t)];
        int64_t *matePosBuf = (int64_t *)matePosData.mutableBytes;
        for (NSUInteger i = 0; i < n; i++) {
            matePosBuf[i] = i < matePositionsCollected.count
                ? [(NSNumber *)matePositionsCollected[i] longLongValue]
                : -1;
        }
        NSMutableData *tlenData = [NSMutableData dataWithLength:n * sizeof(int32_t)];
        int32_t *tlenBuf = (int32_t *)tlenData.mutableBytes;
        for (NSUInteger i = 0; i < n; i++) {
            tlenBuf[i] = i < templateLengthsCollected.count
                ? (int32_t)[(NSNumber *)templateLengthsCollected[i] intValue]
                : 0;
        }

        TTIOWrittenGenomicRun *wgr = [[TTIOWrittenGenomicRun alloc]
            initWithAcquisitionMode:(TTIOAcquisitionMode)((NSNumber *)meta[@"acquisitionMode"]).unsignedIntegerValue
                       referenceUri:referenceUri
                           platform:platform
                         sampleName:sampleName
                          positions:positionsData
                   mappingQualities:mqData
                              flags:flagsData
                          sequences:[gd[@"sequences"] copy]
                          qualities:[gd[@"qualities"] copy]
                            offsets:offsetsData
                            lengths:lengthsData
                             cigars:cigars
                          readNames:readNames
                    mateChromosomes:mateChroms
                      matePositions:matePosData
                    templateLengths:tlenData
                        chromosomes:[gd[@"chromosomes"] copy]
                  signalCompression:TTIOCompressionNone];
        if (readRole.length > 0) wgr.readRole = readRole;
        if ([gd[@"anyTags"] boolValue]) {
            NSArray *tagsCollected = gd[@"tags"] ?: @[];
            NSMutableArray *tags = [NSMutableArray arrayWithCapacity:n];
            for (NSUInteger i = 0; i < n; i++) {
                [tags addObject:(i < tagsCollected.count ? tagsCollected[i] : @"")];
            }
            wgr.tags = tags;
        }
        // Phase 2c-T: attach verbatim blobs collected for this dataset_id.
        TTIOBulkV2Blobs *slot = bulkBlobs[didKey];
        if (slot) wgr.bulkV2Blobs = slot;
        genomicRuns[(NSString *)meta[@"name"]] = wgr;
    }

    BOOL wrote = [TTIOSpectralDataset writeMinimalToPath:outputPath
                                                    title:title
                                       isaInvestigationId:isa
                                                   msRuns:runs
                                              genomicRuns:(genomicRuns.count ? genomicRuns : nil)
                                          identifications:(collectedIdentifications.count
                                                            ? collectedIdentifications : nil)
                                          quantifications:(collectedQuantifications.count
                                                            ? collectedQuantifications : nil)
                                        provenanceRecords:(collectedProvenance.count
                                                            ? collectedProvenance : nil)
                                                    error:error];
    if (!wrote) return NO;

    // v0.11 Task 6.4 (Stage 6): layer any decoded Subjects + Samples
    // onto /study/subjects/<external_id>/ + /study/samples/<sample_id>/
    // per-row HDF5 groups. writeMinimalToPath doesn't accept Stage 6
    // accessors yet, so we reopen the file RW and write the typed
    // attributes directly — mirrors the pattern Java's
    // TransportReader.materializeTo established (commit dd211600,
    // "wrote subjects via setAttribute after the create call"). Soft-FK
    // mismatches surface via NSLog (validation runs upstream, doesn't
    // fail materialise — spec §4.4).
    if (collectedSubjects.count > 0 || collectedSamples.count > 0) {
        @try {
            [TTIOSpectralDataset validateSubjects:collectedSubjects
                                            samples:collectedSamples];
        } @catch (NSException *exc) {
            // Duplicate IDs from a malformed stream — surface as error.
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 [NSString stringWithFormat:
                                     @"SUBJECT/SAMPLE_METADATA embed: %@",
                                     exc.reason]}];
            return NO;
        }
        TTIOHDF5File *f =
            [TTIOHDF5File openAtPath:outputPath error:error];
        if (f == nil) return NO;
        TTIOHDF5Group *rootGroup = [f rootGroup];
        TTIOHDF5Group *studyGroup =
            [rootGroup openGroupNamed:@"study" error:error];
        if (studyGroup == nil) { [f close]; return NO; }
        if (collectedSubjects.count > 0) {
            TTIOHDF5Group *sg =
                [studyGroup createGroupNamed:@"subjects" error:error];
            if (sg == nil) { [f close]; return NO; }
            for (TTIOSubject *s in collectedSubjects) {
                TTIOHDF5Group *row =
                    [sg createGroupNamed:s.externalId error:error];
                if (row == nil) { [f close]; return NO; }
                if (![row setStringAttribute:@"external_id"
                                        value:s.externalId
                                        error:error]) { [f close]; return NO; }
                if (s.project.length > 0) {
                    if (![row setStringAttribute:@"project" value:s.project
                                            error:error]) { [f close]; return NO; }
                }
                if (s.sex.length > 0) {
                    if (![row setStringAttribute:@"sex" value:s.sex
                                            error:error]) { [f close]; return NO; }
                }
                if (![row setIntegerAttribute:@"birth_year"
                                          value:s.birthYear
                                          error:error]) { [f close]; return NO; }
                if (![row setStringAttribute:@"attributes_json"
                                        value:[s attributesJson]
                                        error:error]) { [f close]; return NO; }
            }
        }
        if (collectedSamples.count > 0) {
            TTIOHDF5Group *smg =
                [studyGroup createGroupNamed:@"samples" error:error];
            if (smg == nil) { [f close]; return NO; }
            for (TTIOSample *s in collectedSamples) {
                TTIOHDF5Group *row =
                    [smg createGroupNamed:s.sampleId error:error];
                if (row == nil) { [f close]; return NO; }
                if (![row setStringAttribute:@"sample_id"
                                        value:s.sampleId
                                        error:error]) { [f close]; return NO; }
                if (s.subjectExternalId.length > 0) {
                    if (![row setStringAttribute:@"subject_external_id"
                                            value:s.subjectExternalId
                                            error:error]) { [f close]; return NO; }
                }
                if (s.sampleKind.length > 0) {
                    if (![row setStringAttribute:@"sample_kind"
                                            value:s.sampleKind
                                            error:error]) { [f close]; return NO; }
                }
                if (![row setIntegerAttribute:@"collected_at"
                                          value:s.collectedAt
                                          error:error]) { [f close]; return NO; }
                if (![row setStringAttribute:@"attributes_json"
                                        value:[s attributesJson]
                                        error:error]) { [f close]; return NO; }
            }
        }
        if (![f close]) return NO;
    }

    // v0.11 Task 3.6: embed any image cube decoded from the stream's
    // IMAGE_* packets. The TTI-O on-disk layout puts the image_cube
    // group directly under /study/ (see TTIOMSImage -writeToFilePath:);
    // here we replicate the same H5 layout inline so the production
    // MSImage code remains untouched. Java parity:
    // TransportReader.materializeTo embeds collectedImage via
    // MSImage.writeTo(studyGroup) (commit a6b1e5d9). Python parity:
    // TransportReader.materialize_to writes the cube via the storage
    // provider (commit 1f619ced).
    if (imgCollected) {
        hid_t fid = H5Fopen([outputPath fileSystemRepresentation],
                              H5F_ACC_RDWR, H5P_DEFAULT);
        if (fid < 0) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"image embed: H5Fopen RDWR failed"}];
            return NO;
        }
        hid_t studyGid = H5Gopen2(fid, "study", H5P_DEFAULT);
        if (studyGid < 0) {
            H5Fclose(fid);
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"image embed: H5Gopen2 /study failed"}];
            return NO;
        }
        BOOL msOk = writeImageCubeGroupAtStudy(
            studyGid, "image_cube",
            msImgWidth, msImgHeight, msImgBins, /*tileSize=*/32,
            msImgPxX, msImgPxY, msImgScanSnap,
            msImgAxisSnap, msImgCubeSnap.bytes,
            @{}, @{}, @{},
            error);
        H5Gclose(studyGid);
        H5Fclose(fid);
        if (!msOk) return NO;
    }

    // v0.11 Stage 5.3 (Deferral 1): embed Raman + IR cubes alongside
    // the MS image_cube when the stream carried IMAGE blocks with
    // modality 1 / 2. Each modality lives in its own
    // /study/{raman,ir}_image_cube/ subgroup so all three coexist on
    // the same materialised .tio. Java parity:
    // TransportReader.materializeTo (commit f99ec47d). Python parity:
    // TransportReader.materialize_to (commit 6abead73).
    if (imgRamanCollected || imgIRCollected) {
        hid_t fid = H5Fopen([outputPath fileSystemRepresentation],
                              H5F_ACC_RDWR, H5P_DEFAULT);
        if (fid < 0) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"raman/ir image embed: H5Fopen RDWR failed"}];
            return NO;
        }
        hid_t studyGid = H5Gopen2(fid, "study", H5P_DEFAULT);
        if (studyGid < 0) {
            H5Fclose(fid);
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"raman/ir image embed: H5Gopen2 /study failed"}];
            return NO;
        }
        BOOL ok = YES;
        if (imgRamanCollected) {
            // Raman modality_extras carries excitation_wavelength_nm
            // + laser_power_mw — both materialise as group-level
            // double attrs on /study/raman_image_cube/, matching
            // TTIORamanImage -writeToFilePath:.
            ok = writeImageCubeGroupAtStudy(
                studyGid, "raman_image_cube",
                ramanImgWidth, ramanImgHeight, ramanImgBins, /*tileSize=*/32,
                ramanImgPxX, ramanImgPxY, ramanImgScanSnap,
                ramanImgAxisSnap, ramanImgCubeSnap.bytes,
                @{ @"excitation_wavelength_nm": @(imgRamanExcitationNm),
                   @"laser_power_mw":           @(imgRamanLaserPowerMw) },
                @{},
                @{},
                error);
        }
        if (ok && imgIRCollected) {
            // IR modality_extras → resolution_cm_inv (double attr)
            // + ir_mode (i64 enum: 0=transmittance, 1=absorbance —
            // matches Python's int(IRMode) convention and the typed
            // wire-form now used by TTIOIRImage -writeToFilePath:
            // post-Stage-5.6 cross-language parity fix).
            ok = writeImageCubeGroupAtStudy(
                studyGid, "ir_image_cube",
                irImgWidth, irImgHeight, irImgBins, /*tileSize=*/32,
                irImgPxX, irImgPxY, irImgScanSnap,
                irImgAxisSnap, irImgCubeSnap.bytes,
                @{ @"resolution_cm_inv": @(imgIRResolutionCmInv) },
                @{ @"ir_mode":           @((imgIRModeByte == 1) ? 1LL : 0LL) },
                @{},
                error);
        }
        H5Gclose(studyGid);
        H5Fclose(fid);
        if (!ok) return NO;
    }

    // v0.11 Task 3.4: persist the collected encryption algorithm as the
    // root @encrypted HDF5 attribute. Mirrors how a freshly-encrypted
    // .tio surfaces -isEncrypted / -encryptedAlgorithm — readers consult
    // the same root attribute on open. Java parity:
    // TransportReader.materializeTo (commit 530a5833) which uses
    // SpectralDataset.markRootEncryptedWithEncryptionAlgorithm. Python
    // parity: TransportReader.materialize_to (commit bf38bdc9).
    if (collectedEncryptionAlgorithm.length > 0) {
        TTIOHDF5File *f =
            [TTIOHDF5File openAtPath:outputPath error:error];
        if (!f) return NO;
        TTIOHDF5Group *root = [f rootGroup];
        if (![root setStringAttribute:@"encrypted"
                                  value:collectedEncryptionAlgorithm
                                  error:error]) {
            [f close];
            return NO;
        }
        if (![f close]) return NO;
    }

    // v0.11 Task 3.3: embed any reference groups collected from the
    // stream's REFERENCE_* packets. TTIOReferenceImport.writeToDataset
    // requires an open writable provider; reopen the just-written .tio
    // through TTIOHDF5Provider in ReadWrite mode and attach it to a
    // stub TTIOSpectralDataset via object_setIvar so the public
    // `provider` getter surfaces it (same pattern as
    // TTIOReferenceImportWriteToDatasetTests.m). Java parity:
    // TransportReader.materializeTo embeds collectedRefs before
    // returning. Python parity: TransportReader.materialize_to opens
    // SpectralDataset(writable=True), embeds, then reopens read-only.
    if (collectedRefs.count > 0) {
        TTIOHDF5Provider *p = [[TTIOHDF5Provider alloc] init];
        if (![p openURL:outputPath
                   mode:TTIOStorageOpenModeReadWrite
                  error:error]) {
            return NO;
        }
        TTIOSpectralDataset *embedDs =
            [[TTIOSpectralDataset alloc] initWithTitle:title
                                    isaInvestigationId:isa
                                                msRuns:@{}
                                               nmrRuns:@{}
                                       identifications:@[]
                                       quantifications:@[]
                                     provenanceRecords:@[]
                                           transitions:nil];
        Ivar provIvar =
            class_getInstanceVariable([TTIOSpectralDataset class],
                                        "_provider");
        if (provIvar == NULL) {
            if (error) *error = [NSError errorWithDomain:TTIOTransportErrorDomain
                                                     code:TTIOTransportErrorUnexpectedPayload
                                                 userInfo:@{NSLocalizedDescriptionKey:
                                 @"TTIOSpectralDataset _provider ivar not found "
                                 @"— cannot attach writable provider for "
                                 @"reference embed"}];
            [p close];
            return NO;
        }
        object_setIvar(embedDs, provIvar, p);
        for (TTIOReferenceImport *ref in collectedRefs) {
            NSError *embedErr = nil;
            if (![ref writeToDataset:embedDs
                            overwrite:YES
                                error:&embedErr]) {
                if (error) *error = embedErr;
                [embedDs closeFile];
                return NO;
            }
        }
        [embedDs closeFile];
    }
    return YES;
}

@end

// ---------------------------------------------------------- internal

@implementation TTIOTransportReader (Internal)

- (NSArray<TTIOTransportPacketRecord *> *)recordsForTest
{
    NSError *err = nil;
    return [self readAllPacketsWithError:&err];
}

@end
