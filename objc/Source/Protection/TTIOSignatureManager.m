/*
 * TTIOSignatureManager.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOSignatureManager
 * Declared In:   Protection/TTIOSignatureManager.h
 *
 * HMAC-SHA256 + post-quantum signing/verification helpers.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOSignatureManager.h"
#import "TTIOCipherSuite.h"
#import "TTIOPostQuantumCrypto.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOFeatureFlags.h"
#import <hdf5.h>
#import <openssl/hmac.h>
#import <openssl/evp.h>
#import <string.h>

@implementation TTIOSignatureManager

#pragma mark - Primitive

+ (NSData *)hmacSHA256OfData:(NSData *)data withKey:(NSData *)key
{
    unsigned char mac[EVP_MAX_MD_SIZE];
    unsigned int macLen = 0;
    HMAC(EVP_sha256(),
         key.bytes, (int)key.length,
         data.bytes, data.length,
         mac, &macLen);
    return [NSData dataWithBytes:mac length:macLen];
}

#pragma mark - Helpers

static hid_t openDatasetByPath(hid_t fid, NSString *path)
{
    // H5Dopen2 accepts absolute paths "/a/b/c".
    return H5Dopen2(fid, [path UTF8String], H5P_DEFAULT);
}

static NSData *readDatasetAsBytes(hid_t did)
{
    hid_t typeId = H5Dget_type(did);
    size_t typeSize = H5Tget_size(typeId);

    hid_t space = H5Dget_space(did);
    int rank = H5Sget_simple_extent_ndims(space);
    hsize_t dims[16] = {0};
    H5Sget_simple_extent_dims(space, dims, NULL);

    NSUInteger total = 1;
    for (int i = 0; i < rank; i++) total *= (NSUInteger)dims[i];
    NSUInteger totalBytes = total * typeSize;

    NSMutableData *buf = [NSMutableData dataWithLength:totalBytes];
    H5Dread(did, typeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.mutableBytes);

    H5Sclose(space);
    H5Tclose(typeId);
    return buf;
}

#pragma mark - M18 canonical byte-order helpers

/** Map an atomic file type to its little-endian equivalent memory type.
 *  Returns H5I_INVALID_HID for unsupported classes so the caller can
 *  fall back to a native read. */
static hid_t canonicalLEForAtomic(hid_t fileType)
{
    H5T_class_t cls = H5Tget_class(fileType);
    size_t size = H5Tget_size(fileType);
    if (cls == H5T_FLOAT) {
        if (size == 4) return H5T_IEEE_F32LE;
        if (size == 8) return H5T_IEEE_F64LE;
    } else if (cls == H5T_INTEGER) {
        H5T_sign_t sign = H5Tget_sign(fileType);
        if (sign == H5T_SGN_2) {
            if (size == 1) return H5T_STD_I8LE;
            if (size == 2) return H5T_STD_I16LE;
            if (size == 4) return H5T_STD_I32LE;
            if (size == 8) return H5T_STD_I64LE;
        } else {
            if (size == 1) return H5T_STD_U8LE;
            if (size == 2) return H5T_STD_U16LE;
            if (size == 4) return H5T_STD_U32LE;
            if (size == 8) return H5T_STD_U64LE;
        }
    }
    return H5I_INVALID_HID;
}

static void appendUInt32LE(NSMutableData *out, uint32_t v)
{
    uint8_t b[4] = {
        (uint8_t)(v & 0xff),
        (uint8_t)((v >> 8) & 0xff),
        (uint8_t)((v >> 16) & 0xff),
        (uint8_t)((v >> 24) & 0xff),
    };
    [out appendBytes:b length:4];
}

/** Compound canonical byte stream: for each record in declaration
 *  order, walk members. Numeric members are read via LE memory types
 *  (so the on-disk buffer is already canonical) and appended as-is.
 *  Variable-length string members are emitted as
 *  ``u32_le(byte_length) || bytes`` so padding and pointer layouts
 *  cannot influence the result. */
static NSData *canonicalBytesForCompoundDataset(hid_t did, hid_t fileType,
                                                 hid_t space, hsize_t total)
{
    int nmembers = H5Tget_nmembers(fileType);
    if (nmembers <= 0) {
        NSMutableData *buf = [NSMutableData dataWithLength:0];
        H5Sclose(space);
        H5Tclose(fileType);
        return buf;
    }

    // Build a packed memory type: numerics -> LE equivalent, VL strings
    // stay as VL strings (pointer in the native buffer). Fixed-size
    // strings are emitted via their native type class.
    size_t packedSize = 0;
    size_t *memOffset  = (size_t *)calloc(nmembers, sizeof(size_t));
    hid_t  *memType    = (hid_t  *)calloc(nmembers, sizeof(hid_t));
    BOOL   *memIsVL    = (BOOL   *)calloc(nmembers, sizeof(BOOL));
    BOOL   *memOwnType = (BOOL   *)calloc(nmembers, sizeof(BOOL));

    for (int i = 0; i < nmembers; i++) {
        hid_t native = H5Tget_member_type(fileType, i);
        H5T_class_t cls = H5Tget_class(native);
        if (cls == H5T_STRING && H5Tis_variable_str(native) > 0) {
            // VL string: native memory type holds a char*.
            memOffset[i] = packedSize;
            memType[i]   = native;
            memIsVL[i]   = YES;
            memOwnType[i] = YES;  // close after use
            packedSize += sizeof(char *);
        } else if (cls == H5T_FLOAT || cls == H5T_INTEGER) {
            hid_t le = canonicalLEForAtomic(native);
            if (le == H5I_INVALID_HID) le = native;
            memOffset[i] = packedSize;
            memType[i]   = le;
            memIsVL[i]   = NO;
            memOwnType[i] = NO;  // built-in IDs, don't close
            packedSize += H5Tget_size(le);
            H5Tclose(native);
        } else {
            // Fall back to native member type as-is (fixed strings,
            // nested compounds, ...). Rare in practice for TTIO.
            memOffset[i] = packedSize;
            memType[i]   = native;
            memIsVL[i]   = NO;
            memOwnType[i] = YES;
            packedSize += H5Tget_size(native);
        }
    }

    hid_t packedType = H5Tcreate(H5T_COMPOUND, packedSize);
    for (int i = 0; i < nmembers; i++) {
        char *mname = H5Tget_member_name(fileType, i);
        H5Tinsert(packedType, mname, memOffset[i], memType[i]);
        H5free_memory(mname);
    }

    void *rawBuf = calloc((size_t)total, packedSize);
    H5Dread(did, packedType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawBuf);

    NSMutableData *out = [NSMutableData data];
    for (hsize_t r = 0; r < total; r++) {
        const uint8_t *record = (const uint8_t *)rawBuf + r * packedSize;
        for (int i = 0; i < nmembers; i++) {
            const uint8_t *fieldPtr = record + memOffset[i];
            if (memIsVL[i]) {
                const char *cstr = *(const char *const *)fieldPtr;
                uint32_t blen = cstr ? (uint32_t)strlen(cstr) : 0;
                appendUInt32LE(out, blen);
                if (cstr && blen > 0) [out appendBytes:cstr length:blen];
            } else {
                size_t fsize = H5Tget_size(memType[i]);
                [out appendBytes:fieldPtr length:fsize];
            }
        }
    }

    H5Dvlen_reclaim(packedType, space, H5P_DEFAULT, rawBuf);
    free(rawBuf);
    H5Tclose(packedType);
    for (int i = 0; i < nmembers; i++) {
        if (memOwnType[i]) H5Tclose(memType[i]);
    }
    free(memOffset);
    free(memType);
    free(memIsVL);
    free(memOwnType);
    H5Sclose(space);
    H5Tclose(fileType);
    return out;
}

/** v2 canonical byte reader.
 *
 *  Atomic numeric datasets are read via the host's little-endian mem
 *  type so the resulting buffer is canonical on any architecture.
 *  Compound datasets dispatch to :func:`canonicalBytesForCompoundDataset`.
 *  Any other class (fixed-size strings, enums, ...) falls through to
 *  the native-bytes form, matching the v1 behaviour.
 */
static NSData *readDatasetCanonical(hid_t did)
{
    hid_t fileType = H5Dget_type(did);
    hid_t space    = H5Dget_space(did);
    int rank = H5Sget_simple_extent_ndims(space);
    hsize_t dims[16] = {0};
    H5Sget_simple_extent_dims(space, dims, NULL);
    hsize_t total = 1;
    for (int i = 0; i < rank; i++) total *= dims[i];

    H5T_class_t cls = H5Tget_class(fileType);

    if (cls == H5T_FLOAT || cls == H5T_INTEGER) {
        hid_t memType = canonicalLEForAtomic(fileType);
        if (memType == H5I_INVALID_HID) memType = fileType;
        size_t elemSize = H5Tget_size(memType);
        NSMutableData *buf = [NSMutableData dataWithLength:(NSUInteger)total * elemSize];
        H5Dread(did, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.mutableBytes);
        H5Sclose(space);
        H5Tclose(fileType);
        return buf;
    }

    if (cls == H5T_COMPOUND) {
        return canonicalBytesForCompoundDataset(did, fileType, space, total);
    }

    // Unsupported class: fall back to native bytes (v1 path behaviour).
    size_t size = H5Tget_size(fileType);
    NSMutableData *buf = [NSMutableData dataWithLength:(NSUInteger)total * size];
    H5Dread(did, fileType, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.mutableBytes);
    H5Sclose(space);
    H5Tclose(fileType);
    return buf;
}

static NSString *const kTTIOSignatureV2Prefix = @"v2:";
static NSString *const kTTIOSignatureV3Prefix = @"v3:";  // v0.8 M49.1

static BOOL writeStringAttribute(hid_t locId, const char *name, NSString *value)
{
    if (H5Aexists(locId, name) > 0) H5Adelete(locId, name);
    hid_t strType = H5Tcopy(H5T_C_S1);
    H5Tset_size(strType, H5T_VARIABLE);
    hid_t scalar = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(locId, name, strType, scalar, H5P_DEFAULT, H5P_DEFAULT);
    const char *cs = [value UTF8String];
    herr_t rc = H5Awrite(attr, strType, &cs);
    H5Aclose(attr);
    H5Sclose(scalar);
    H5Tclose(strType);
    return rc >= 0;
}

static NSString *readStringAttribute(hid_t locId, const char *name)
{
    if (H5Aexists(locId, name) <= 0) return nil;
    hid_t attr = H5Aopen(locId, name, H5P_DEFAULT);
    hid_t t = H5Aget_type(attr);
    NSString *out = @"";
    if (H5Tis_variable_str(t) > 0) {
        char *cs = NULL;
        H5Aread(attr, t, &cs);
        if (cs) {
            out = [NSString stringWithUTF8String:cs];
            H5free_memory(cs);
        }
    } else {
        // Fixed-size string: TTIOHDF5Group.setStringAttribute writes
        // these by default.
        size_t size = H5Tget_size(t);
        char *buf = (char *)calloc(size + 1, 1);
        if (H5Aread(attr, t, buf) >= 0) {
            out = [NSString stringWithUTF8String:buf];
        }
        free(buf);
    }
    H5Tclose(t);
    H5Aclose(attr);
    return out;
}

static BOOL ensureSignatureFeatureFlags(TTIOHDF5Group *root, NSError **error)
{
    NSArray *features = [TTIOFeatureFlags featuresForRoot:root];
    BOOL hasDig   = [features containsObject:[TTIOFeatureFlags featureDigitalSignatures]];
    BOOL hasCanon = [features containsObject:[TTIOFeatureFlags featureCanonicalSignatures]];
    if (hasDig && hasCanon) return YES;
    NSMutableArray *updated = [features mutableCopy];
    if (!hasDig)   [updated addObject:[TTIOFeatureFlags featureDigitalSignatures]];
    if (!hasCanon) [updated addObject:[TTIOFeatureFlags featureCanonicalSignatures]];
    NSString *version = [TTIOFeatureFlags formatVersionForRoot:root] ?: @"1.1";
    return [TTIOFeatureFlags writeFormatVersion:version
                                       features:updated
                                         toRoot:root
                                          error:error];
}

#pragma mark - Dataset signing

+ (BOOL)signDataset:(NSString *)datasetPath
             inFile:(NSString *)filePath
            withKey:(NSData *)hmacKey
              error:(NSError **)error
{
    TTIOHDF5File *file = [TTIOHDF5File openAtPath:filePath error:error];
    if (!file) return NO;

    hid_t did = openDatasetByPath([file rootGroup].groupId, datasetPath);
    if (did < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"cannot open dataset %@ for signing", datasetPath);
        [file close];
        return NO;
    }

    // canonical little-endian serialization before hashing so
    // the MAC is portable across host endianness. The stored string
    // carries a ``v2:`` prefix to distinguish it from v0.2 native
    // signatures. Verifiers still accept unprefixed v1 signatures.
    NSData *bytes = readDatasetCanonical(did);
    NSData *mac = [self hmacSHA256OfData:bytes withKey:hmacKey];
    NSString *b64 = [kTTIOSignatureV2Prefix stringByAppendingString:
                     [mac base64EncodedStringWithOptions:0]];

    BOOL ok = writeStringAttribute(did, "ttio_signature", b64);
    H5Dclose(did);
    if (!ok) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeWrite,
            @"failed to write @ttio_signature on %@", datasetPath);
        [file close];
        return NO;
    }

    // Ensure the root feature flags include opt_digital_signatures and
    // opt_canonical_signatures.
    ensureSignatureFeatureFlags([file rootGroup], error);
    return [file close];
}

+ (BOOL)verifyDataset:(NSString *)datasetPath
               inFile:(NSString *)filePath
              withKey:(NSData *)hmacKey
                error:(NSError **)error
{
    TTIOHDF5File *file = [TTIOHDF5File openReadOnlyAtPath:filePath error:error];
    if (!file) return NO;

    hid_t did = openDatasetByPath([file rootGroup].groupId, datasetPath);
    if (did < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"cannot open dataset %@ for verification", datasetPath);
        [file close];
        return NO;
    }

    NSString *storedB64 = readStringAttribute(did, "ttio_signature");
    if (!storedB64) {
        H5Dclose(did);
        [file close];
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"dataset %@ has no @ttio_signature attribute", datasetPath);
        return NO;
    }

    BOOL isV2 = [storedB64 hasPrefix:kTTIOSignatureV2Prefix];
    NSString *payloadB64 = isV2
        ? [storedB64 substringFromIndex:kTTIOSignatureV2Prefix.length]
        : storedB64;
    NSData *storedMac = [[NSData alloc] initWithBase64EncodedString:payloadB64
                                                              options:0];

    NSData *bytes = isV2 ? readDatasetCanonical(did) : readDatasetAsBytes(did);
    NSData *computedMac = [self hmacSHA256OfData:bytes withKey:hmacKey];
    H5Dclose(did);
    [file close];

    if (!storedMac || ![storedMac isEqualToData:computedMac]) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"signature mismatch on %@ (tampered or wrong key)", datasetPath);
        return NO;
    }
    return YES;
}

#pragma mark - Algorithm-dispatched sign / verify (v0.8 M49.1)

// Add opt_pqc_preview to the root feature list (idempotent). Used
// whenever a v3: signature is written or read.
static BOOL markPQCPreviewFeature(TTIOHDF5Group *root, NSError **error)
{
    NSArray<NSString *> *features = [TTIOFeatureFlags featuresForRoot:root];
    NSString *flag = [TTIOFeatureFlags featurePQCPreview];
    if ([features containsObject:flag]) return YES;
    NSMutableArray *updated = [features mutableCopy] ?: [NSMutableArray array];
    [updated addObject:flag];
    NSString *version = [TTIOFeatureFlags formatVersionForRoot:root] ?: @"1.2";
    return [TTIOFeatureFlags writeFormatVersion:version
                                        features:updated
                                          toRoot:root
                                           error:error];
}

+ (BOOL)signDataset:(NSString *)datasetPath
             inFile:(NSString *)filePath
            withKey:(NSData *)key
           algorithm:(NSString *)algorithm
              error:(NSError **)error
{
    if ([algorithm isEqualToString:@"hmac-sha256"]) {
        return [self signDataset:datasetPath
                          inFile:filePath
                         withKey:key
                           error:error];
    }
    if (![algorithm isEqualToString:@"ml-dsa-87"]) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"signDataset: algorithm %@ not supported", algorithm);
        return NO;
    }
    if (![TTIOCipherSuite validatePrivateKey:key
                                    algorithm:@"ml-dsa-87"
                                        error:error]) {
        return NO;
    }

    TTIOHDF5File *file = [TTIOHDF5File openAtPath:filePath error:error];
    if (!file) return NO;

    hid_t did = openDatasetByPath([file rootGroup].groupId, datasetPath);
    if (did < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"cannot open dataset %@ for signing", datasetPath);
        [file close];
        return NO;
    }

    NSData *canonical = readDatasetCanonical(did);
    NSData *sig = [TTIOPostQuantumCrypto sigSignWithPrivateKey:key
                                                         message:canonical
                                                           error:error];
    if (!sig) {
        H5Dclose(did);
        [file close];
        return NO;
    }
    NSString *stored = [kTTIOSignatureV3Prefix stringByAppendingString:
        [sig base64EncodedStringWithOptions:0]];

    BOOL ok = writeStringAttribute(did, "ttio_signature", stored);
    H5Dclose(did);
    if (!ok) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeWrite,
            @"failed to write @ttio_signature on %@", datasetPath);
        [file close];
        return NO;
    }

    // Feature flag updates: v3 implies both opt_digital_signatures and
    // opt_pqc_preview. The canonical-sig flag is informative (v3 is
    // canonical by construction), so we keep it on for consistency.
    ensureSignatureFeatureFlags([file rootGroup], error);
    markPQCPreviewFeature([file rootGroup], error);
    return [file close];
}

+ (BOOL)verifyDataset:(NSString *)datasetPath
               inFile:(NSString *)filePath
              withKey:(NSData *)key
            algorithm:(NSString *)algorithm
                error:(NSError **)error
{
    TTIOHDF5File *file = [TTIOHDF5File openReadOnlyAtPath:filePath error:error];
    if (!file) return NO;

    hid_t did = openDatasetByPath([file rootGroup].groupId, datasetPath);
    if (did < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"cannot open dataset %@ for verification", datasetPath);
        [file close];
        return NO;
    }

    NSString *stored = readStringAttribute(did, "ttio_signature");
    if (!stored) {
        H5Dclose(did);
        [file close];
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"dataset %@ has no @ttio_signature attribute", datasetPath);
        return NO;
    }

    BOOL storedIsV3 = [stored hasPrefix:kTTIOSignatureV3Prefix];

    // Algorithm/prefix must match — don't silently accept a crossed check.
    if (storedIsV3 && ![algorithm isEqualToString:@"ml-dsa-87"]) {
        H5Dclose(did);
        [file close];
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"stored signature is v3 (ml-dsa-87) but caller passed "
            @"algorithm=%@", algorithm);
        return NO;
    }
    if ([algorithm isEqualToString:@"ml-dsa-87"] && !storedIsV3) {
        H5Dclose(did);
        [file close];
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"stored signature is not v3 (ml-dsa-87) — pass "
            @"algorithm=hmac-sha256 to verify legacy signatures");
        return NO;
    }

    if (storedIsV3) {
        if (![TTIOCipherSuite validatePublicKey:key
                                       algorithm:@"ml-dsa-87"
                                           error:error]) {
            H5Dclose(did);
            [file close];
            return NO;
        }
        NSString *payloadB64 = [stored substringFromIndex:
            kTTIOSignatureV3Prefix.length];
        NSData *sig = [[NSData alloc] initWithBase64EncodedString:payloadB64
                                                            options:0];
        if (!sig) {
            H5Dclose(did);
            [file close];
            if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
                @"v3: base64 payload malformed on %@", datasetPath);
            return NO;
        }
        NSData *canonical = readDatasetCanonical(did);
        H5Dclose(did);
        [file close];
        NSError *pqcErr = nil;
        BOOL ok = [TTIOPostQuantumCrypto sigVerifyWithPublicKey:key
                                                         message:canonical
                                                       signature:sig
                                                           error:&pqcErr];
        if (!ok) {
            if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
                @"v3 signature failed verification on %@", datasetPath);
            return NO;
        }
        return YES;
    }

    // HMAC (v1/v2) path — delegate to the legacy entry point.
    H5Dclose(did);
    [file close];
    return [self verifyDataset:datasetPath
                        inFile:filePath
                       withKey:key
                         error:error];
}

#pragma mark - M90.2 Genomic run-level convenience

// Channels signed by signGenomicRun: / verifyGenomicRun:. M90.15:
// chromosomes (a VL_BYTES row compound) is now included alongside the
// atomic columns. readDatasetCanonical handles compound datasets via
// canonicalBytesForCompoundDataset, which already serialises
// VL_STRING fields as u32_le(length) || utf-8_bytes — so signing the
// chromosomes compound parallels Python's M90.15 path.
// M101: the SAM tags channel (a flat uint8 SAM_TAGS stream) is signed
// with the bases it describes, when present.
static NSString *const kSignalChannelNames[] = {
    @"sequences", @"qualities", @"tags",
};
// L1 (Task #82 Phase B.1, 2026-05-01): chromosomes are decomposed
// into chromosome_ids (uint16) + chromosome_names (compound) — both
// new columns are signed in place of the old single column.
static NSString *const kIndexColumnNames[] = {
    @"offsets", @"lengths", @"positions", @"mapping_qualities", @"flags",
    @"chromosome_ids", @"chromosome_names",
};

// Test whether ``parentPath/childName`` exists in the file. Used to
// silently skip datasets that aren't present (e.g. encrypted files
// where signal channels become *_segments compounds).
static BOOL fileHasObjectAtPath(NSString *filePath, NSString *objectPath)
{
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:filePath error:NULL];
    if (!f) return NO;
    htri_t exists = H5Lexists([f rootGroup].groupId,
                                [objectPath UTF8String], H5P_DEFAULT);
    [f close];
    return exists > 0;
}

// The signal-channel datasets of a run, as sub-paths under the run
// group. A channel stored as a group (blocks_v1 sequences/data, the
// v1.8 refdiff_v2 group) contributes each child dataset; a plain
// dataset contributes itself. blocks/index follows when present.
static NSArray<NSString *> *genomicSignedSubPaths(NSString *filePath, NSString *runPath)
{
    NSMutableArray<NSString *> *out = [NSMutableArray array];
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:filePath error:NULL];
    if (!f) return out;
    TTIOHDF5Group *rootGroup = [f rootGroup];
    hid_t root = rootGroup.groupId;
    for (size_t i = 0; i < sizeof(kSignalChannelNames) / sizeof(kSignalChannelNames[0]); i++) {
        NSString *sub = [NSString stringWithFormat:@"signal_channels/%@", kSignalChannelNames[i]];
        NSString *full = [NSString stringWithFormat:@"%@/%@", runPath, sub];
        if (H5Lexists(root, [full UTF8String], H5P_DEFAULT) <= 0) continue;
        H5O_info2_t info; memset(&info, 0, sizeof(info));
        if (H5Oget_info_by_name3(root, [full UTF8String], &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) continue;
        if (info.type != H5O_TYPE_GROUP) { [out addObject:sub]; continue; }
        hid_t g = H5Gopen2(root, [full UTF8String], H5P_DEFAULT);
        if (g < 0) continue;
        H5G_info_t ginfo; memset(&ginfo, 0, sizeof(ginfo));
        if (H5Gget_info(g, &ginfo) >= 0) {
            for (hsize_t k = 0; k < ginfo.nlinks; k++) {
                char nameBuf[512];
                ssize_t len = H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, k,
                                                 nameBuf, sizeof(nameBuf), H5P_DEFAULT);
                if (len <= 0) continue;
                H5O_info2_t cinfo; memset(&cinfo, 0, sizeof(cinfo));
                if (H5Oget_info_by_name3(g, nameBuf, &cinfo, H5O_INFO_BASIC, H5P_DEFAULT) < 0) continue;
                if (cinfo.type == H5O_TYPE_DATASET) {
                    [out addObject:[NSString stringWithFormat:@"%@/%s", sub, nameBuf]];
                }
            }
        }
        H5Gclose(g);
    }
    NSString *blocks = [NSString stringWithFormat:@"%@/blocks", runPath];
    if (H5Lexists(root, [blocks UTF8String], H5P_DEFAULT) > 0
        && H5Lexists(root, [[blocks stringByAppendingString:@"/index"] UTF8String], H5P_DEFAULT) > 0) {
        [out addObject:@"blocks/index"];
    }
    (void)rootGroup;
    [f close];
    return out;
}

+ (NSDictionary<NSString *, NSString *> *)
    signGenomicRun:(NSString *)runName
            inFile:(NSString *)filePath
           withKey:(NSData *)hmacKey
             error:(NSError **)error
{
    NSMutableDictionary<NSString *, NSString *> *out =
        [NSMutableDictionary dictionary];
    NSString *runPath =
        [NSString stringWithFormat:@"/study/genomic_runs/%@", runName];

    // ── signal_channels/{sequences, qualities} datasets, blocks/index
    for (NSString *subPath in genomicSignedSubPaths(filePath, runPath)) {
        NSString *fullPath =
            [NSString stringWithFormat:@"%@/%@", runPath, subPath];
        if (!fileHasObjectAtPath(filePath, fullPath)) continue;
        NSError *signErr = nil;
        if (![self signDataset:fullPath
                        inFile:filePath
                       withKey:hmacKey
                         error:&signErr]) {
            if (error) *error = signErr;
            return nil;
        }
        // Read back the stored attribute so the caller has the
        // same return shape as Python's sign_genomic_run.
        TTIOHDF5File *file =
            [TTIOHDF5File openReadOnlyAtPath:filePath error:NULL];
        if (file) {
            hid_t did = openDatasetByPath([file rootGroup].groupId, fullPath);
            if (did >= 0) {
                NSString *stored =
                    readStringAttribute(did, "ttio_signature");
                H5Dclose(did);
                if (stored) out[subPath] = stored;
            }
            [file close];
        }
    }

    // ── genomic_index/{offsets, lengths, positions,
    //                    mapping_qualities, flags}
    for (size_t i = 0; i < sizeof(kIndexColumnNames) / sizeof(kIndexColumnNames[0]); i++) {
        NSString *cname = kIndexColumnNames[i];
        NSString *subPath =
            [NSString stringWithFormat:@"genomic_index/%@", cname];
        NSString *fullPath =
            [NSString stringWithFormat:@"%@/%@", runPath, subPath];
        if (!fileHasObjectAtPath(filePath, fullPath)) continue;
        NSError *signErr = nil;
        if (![self signDataset:fullPath
                        inFile:filePath
                       withKey:hmacKey
                         error:&signErr]) {
            if (error) *error = signErr;
            return nil;
        }
        TTIOHDF5File *file =
            [TTIOHDF5File openReadOnlyAtPath:filePath error:NULL];
        if (file) {
            hid_t did = openDatasetByPath([file rootGroup].groupId, fullPath);
            if (did >= 0) {
                NSString *stored =
                    readStringAttribute(did, "ttio_signature");
                H5Dclose(did);
                if (stored) out[subPath] = stored;
            }
            [file close];
        }
    }
    return out;
}

+ (BOOL)verifyGenomicRun:(NSString *)runName
                  inFile:(NSString *)filePath
                 withKey:(NSData *)hmacKey
                   error:(NSError **)error
{
    NSString *runPath =
        [NSString stringWithFormat:@"/study/genomic_runs/%@", runName];
    for (NSString *subPath in genomicSignedSubPaths(filePath, runPath)) {
        NSString *fullPath =
            [NSString stringWithFormat:@"%@/%@", runPath, subPath];
        if (!fileHasObjectAtPath(filePath, fullPath)) continue;
        NSError *vErr = nil;
        if (![self verifyDataset:fullPath
                          inFile:filePath
                         withKey:hmacKey
                           error:&vErr]) {
            if (error) *error = vErr;
            return NO;
        }
    }
    for (size_t i = 0; i < sizeof(kIndexColumnNames) / sizeof(kIndexColumnNames[0]); i++) {
        NSString *cname = kIndexColumnNames[i];
        NSString *fullPath =
            [NSString stringWithFormat:@"%@/genomic_index/%@",
                runPath, cname];
        if (!fileHasObjectAtPath(filePath, fullPath)) continue;
        NSError *vErr = nil;
        if (![self verifyDataset:fullPath
                          inFile:filePath
                         withKey:hmacKey
                           error:&vErr]) {
            if (error) *error = vErr;
            return NO;
        }
    }
    return YES;
}


#pragma mark - M98 assembly graphs

// The signed set is the segments/sequences byte channel — signatures
// cover it the way they cover genomic-run channels. The compound
// tables (records, links, paths, extras, line_index) are outside the
// signed set the way the pre-M90.15 chromosomes compound was.

+ (NSDictionary<NSString *, NSString *> *)
    signAssemblyGraph:(NSString *)graphName
               inFile:(NSString *)filePath
              withKey:(NSData *)hmacKey
                error:(NSError **)error
{
    NSMutableDictionary<NSString *, NSString *> *out =
        [NSMutableDictionary dictionary];
    NSString *subPath = @"segments/sequences";
    NSString *fullPath = [NSString stringWithFormat:
        @"/study/assembly_graphs/%@/%@", graphName, subPath];
    if (!fileHasObjectAtPath(filePath, fullPath)) return out;
    NSError *signErr = nil;
    if (![self signDataset:fullPath
                    inFile:filePath
                   withKey:hmacKey
                     error:&signErr]) {
        if (error) *error = signErr;
        return nil;
    }
    TTIOHDF5File *file =
        [TTIOHDF5File openReadOnlyAtPath:filePath error:NULL];
    if (file) {
        hid_t did = openDatasetByPath([file rootGroup].groupId, fullPath);
        if (did >= 0) {
            NSString *stored = readStringAttribute(did, "ttio_signature");
            H5Dclose(did);
            if (stored) out[subPath] = stored;
        }
        [file close];
    }
    return out;
}

+ (BOOL)verifyAssemblyGraph:(NSString *)graphName
                     inFile:(NSString *)filePath
                    withKey:(NSData *)hmacKey
                      error:(NSError **)error
{
    NSString *fullPath = [NSString stringWithFormat:
        @"/study/assembly_graphs/%@/segments/sequences", graphName];
    if (!fileHasObjectAtPath(filePath, fullPath)) return YES;
    return [self verifyDataset:fullPath
                        inFile:filePath
                       withKey:hmacKey
                         error:error];
}


#pragma mark - Provenance signing

+ (BOOL)signProvenanceInRun:(NSString *)runPath
                     inFile:(NSString *)filePath
                    withKey:(NSData *)hmacKey
                      error:(NSError **)error
{
    TTIOHDF5File *file = [TTIOHDF5File openAtPath:filePath error:error];
    if (!file) return NO;

    hid_t gid = H5Gopen2([file rootGroup].groupId, [runPath UTF8String], H5P_DEFAULT);
    if (gid < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorGroupOpen,
            @"cannot open run %@ for provenance signing", runPath);
        [file close];
        return NO;
    }

    NSString *json = readStringAttribute(gid, "provenance_json");
    if (!json) {
        H5Gclose(gid);
        [file close];
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"run %@ has no @provenance_json to sign", runPath);
        return NO;
    }
    NSData *jsonBytes = [json dataUsingEncoding:NSUTF8StringEncoding];
    NSData *mac = [self hmacSHA256OfData:jsonBytes withKey:hmacKey];
    NSString *b64 = [mac base64EncodedStringWithOptions:0];

    BOOL ok = writeStringAttribute(gid, "provenance_signature", b64);
    H5Gclose(gid);
    if (!ok) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeWrite,
            @"failed to write @provenance_signature on %@", runPath);
        [file close];
        return NO;
    }

    ensureSignatureFeatureFlags([file rootGroup], error);
    return [file close];
}

+ (BOOL)verifyProvenanceInRun:(NSString *)runPath
                       inFile:(NSString *)filePath
                      withKey:(NSData *)hmacKey
                        error:(NSError **)error
{
    TTIOHDF5File *file = [TTIOHDF5File openReadOnlyAtPath:filePath error:error];
    if (!file) return NO;

    hid_t gid = H5Gopen2([file rootGroup].groupId, [runPath UTF8String], H5P_DEFAULT);
    if (gid < 0) {
        if (error) *error = TTIOMakeError(TTIOErrorGroupOpen,
            @"cannot open run %@ for provenance verification", runPath);
        [file close];
        return NO;
    }

    NSString *json = readStringAttribute(gid, "provenance_json");
    NSString *storedB64 = readStringAttribute(gid, "provenance_signature");
    H5Gclose(gid);
    [file close];

    if (!json || !storedB64) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"run %@ missing provenance_json or provenance_signature", runPath);
        return NO;
    }
    NSData *jsonBytes = [json dataUsingEncoding:NSUTF8StringEncoding];
    NSData *computedMac = [self hmacSHA256OfData:jsonBytes withKey:hmacKey];
    NSData *storedMac = [[NSData alloc] initWithBase64EncodedString:storedB64
                                                              options:0];
    if (!storedMac || ![storedMac isEqualToData:computedMac]) {
        if (error) *error = TTIOMakeError(TTIOErrorAttributeRead,
            @"provenance signature mismatch on %@", runPath);
        return NO;
    }
    return YES;
}

@end
