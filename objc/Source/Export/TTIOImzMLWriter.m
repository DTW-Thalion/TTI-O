/*
 * TTIOImzMLWriter.m
 * TTI-O Objective-C Implementation
 *
 * Classes:       TTIOImzMLWriteResult, TTIOImzMLWriter
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Export/TTIOImzMLWriter.h
 *
 * imzML + .ibd exporter. Reverses TTIOImzMLReader: emits the canonical
 * paired .imzML / .ibd files for both continuous and processed mode
 * mass-spectrometry imaging data. Output passes pyimzml's ImzMLParser
 * and round-trips through TTIOImzMLReader bit-identically.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#import "TTIOImzMLWriter.h"
#import "Import/TTIOImzMLReader.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Spectra/TTIOSpectrum.h"
#import "Core/TTIOSignalArray.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import <openssl/sha.h>

// Mirrors Java ImzMLWriter.PROGRESS_INTERVAL_PIXELS (100).
const NSUInteger TTIOImzMLWriterProgressIntervalPixels = 100;


@implementation TTIOImzMLWriteResult {
    NSString *_imzmlPath;
    NSString *_ibdPath;
    NSString *_uuidHex;
    NSString *_mode;
    NSUInteger _nPixels;
}
- (instancetype)initInternal:(NSString *)imzml
                         ibd:(NSString *)ibd
                        uuid:(NSString *)uuid
                        mode:(NSString *)mode
                     nPixels:(NSUInteger)n
{
    self = [super init];
    if (self) {
        _imzmlPath = [imzml copy];
        _ibdPath = [ibd copy];
        _uuidHex = [uuid copy];
        _mode = [mode copy];
        _nPixels = n;
    }
    return self;
}
- (NSString *)imzmlPath { return _imzmlPath; }
- (NSString *)ibdPath   { return _ibdPath; }
- (NSString *)uuidHex   { return _uuidHex; }
- (NSString *)mode      { return _mode; }
- (NSUInteger)nPixels   { return _nPixels; }
@end


// -------------------------------------------------------------- helpers ---

static NSString *NormaliseUUID(NSString *raw) {
    // Strip braces / dashes / whitespace then lowercase. Doing the
    // lowercase via NSString's UTF8 lowercaseString gives correct
    // results for ASCII hex characters; the tolower(unichar) trick
    // only works for the ASCII subset and round-trips weirdly when
    // %c is fed a wider unichar value.
    NSCharacterSet *strip = [NSCharacterSet characterSetWithCharactersInString:@"{}- \t\n\r"];
    NSString *kept = [[raw componentsSeparatedByCharactersInSet:strip]
                       componentsJoinedByString:@""];
    return [kept lowercaseString];
}

static NSString *RandomUUIDHex(void) {
    NSUUID *u = [NSUUID UUID];
    uuid_t bytes;
    [u getUUIDBytes:bytes];
    NSMutableString *hex = [NSMutableString stringWithCapacity:32];
    for (int i = 0; i < 16; i++) [hex appendFormat:@"%02x", bytes[i]];
    return hex;
}

static NSData *HexToBytes(NSString *hex) {
    NSMutableData *data = [NSMutableData dataWithCapacity:hex.length / 2];
    for (NSUInteger i = 0; i + 1 < hex.length; i += 2) {
        NSString *byteStr = [hex substringWithRange:NSMakeRange(i, 2)];
        unsigned v = 0;
        [[NSScanner scannerWithString:byteStr] scanHexInt:&v];
        uint8_t b = (uint8_t)v;
        [data appendBytes:&b length:1];
    }
    return data;
}

static NSString *SHA1Hex(NSData *data) {
    uint8_t digest[SHA_DIGEST_LENGTH];
    SHA1(data.bytes, data.length, digest);
    NSMutableString *hex = [NSMutableString stringWithCapacity:SHA_DIGEST_LENGTH * 2];
    for (int i = 0; i < SHA_DIGEST_LENGTH; i++) [hex appendFormat:@"%02x", digest[i]];
    return hex;
}

static NSString *XmlEscape(NSString *s) {
    NSMutableString *out = [NSMutableString stringWithCapacity:s.length];
    for (NSUInteger i = 0; i < s.length; i++) {
        unichar c = [s characterAtIndex:i];
        switch (c) {
            case '&':  [out appendString:@"&amp;"];  break;
            case '<':  [out appendString:@"&lt;"];   break;
            case '>':  [out appendString:@"&gt;"];   break;
            case '"':  [out appendString:@"&quot;"]; break;
            case '\'': [out appendString:@"&apos;"]; break;
            default:   [out appendFormat:@"%C", c];
        }
    }
    return out;
}

static NSError *MakeError(NSInteger code, NSString *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    NSString *msg = [[NSString alloc] initWithFormat:fmt arguments:ap];
    va_end(ap);
    return [NSError errorWithDomain:@"TTIOImzMLWriter" code:code
                             userInfo:@{NSLocalizedDescriptionKey: msg}];
}


@implementation TTIOImzMLWriter

+ (nullable TTIOImzMLWriteResult *)writePixels:(NSArray<TTIOImzMLPixelSpectrum *> *)pixels
                                   toImzMLPath:(NSString *)imzmlPath
                                       ibdPath:(nullable NSString *)ibdPath
                                          mode:(NSString *)mode
                                     gridMaxX:(NSInteger)gridMaxX
                                     gridMaxY:(NSInteger)gridMaxY
                                     gridMaxZ:(NSInteger)gridMaxZ
                                   pixelSizeX:(double)pixelSizeX
                                   pixelSizeY:(double)pixelSizeY
                                   scanPattern:(NSString *)scanPattern
                                       uuidHex:(nullable NSString *)uuidHex
                                         error:(NSError **)error
{
    return [self writePixels:pixels
                 toImzMLPath:imzmlPath
                     ibdPath:ibdPath
                        mode:mode
                    gridMaxX:gridMaxX
                    gridMaxY:gridMaxY
                    gridMaxZ:gridMaxZ
                  pixelSizeX:pixelSizeX
                  pixelSizeY:pixelSizeY
                 scanPattern:scanPattern
                     uuidHex:uuidHex
                    progress:nil
                       error:error];
}

+ (nullable TTIOImzMLWriteResult *)writePixels:(NSArray<TTIOImzMLPixelSpectrum *> *)pixels
                                   toImzMLPath:(NSString *)imzmlPath
                                       ibdPath:(nullable NSString *)ibdPath
                                          mode:(NSString *)mode
                                     gridMaxX:(NSInteger)gridMaxX
                                     gridMaxY:(NSInteger)gridMaxY
                                     gridMaxZ:(NSInteger)gridMaxZ
                                   pixelSizeX:(double)pixelSizeX
                                   pixelSizeY:(double)pixelSizeY
                                   scanPattern:(NSString *)scanPattern
                                       uuidHex:(nullable NSString *)uuidHex
                                       progress:(TTIOProgressBlock)progress
                                         error:(NSError **)error
{
    if (progress == nil) progress = TTIOProgressDiscard();
    if (![mode isEqualToString:@"continuous"] &&
        ![mode isEqualToString:@"processed"]) {
        if (error) *error = MakeError(1, @"mode must be 'continuous' or 'processed', got '%@'", mode);
        return nil;
    }
    if (pixels.count == 0) {
        if (error) *error = MakeError(2, @"at least one pixel spectrum is required");
        return nil;
    }

    NSString *ibd = ibdPath;
    if (!ibd) {
        ibd = [imzmlPath stringByDeletingPathExtension];
        ibd = [ibd stringByAppendingPathExtension:@"ibd"];
    }

    NSString *uuid = uuidHex ? NormaliseUUID(uuidHex) : RandomUUIDHex();
    if (uuid.length != 32) {
        if (error) *error = MakeError(3, @"uuidHex must be 32 hex chars after normalisation, got %lu",
                                       (unsigned long)uuid.length);
        return nil;
    }

    // Derive grid extents when not supplied.
    if (gridMaxX == 0) {
        NSInteger m = 0;
        for (TTIOImzMLPixelSpectrum *p in pixels) if (p.x > m) m = p.x;
        gridMaxX = m;
    }
    if (gridMaxY == 0) {
        NSInteger m = 0;
        for (TTIOImzMLPixelSpectrum *p in pixels) if (p.y > m) m = p.y;
        gridMaxY = m;
    }
    if (gridMaxZ == 0) gridMaxZ = 1;

    // ── .ibd assembly ──────────────────────────────────────────
    NSMutableData *ibdBuf = [NSMutableData dataWithCapacity:1024];
    [ibdBuf appendData:HexToBytes(uuid)];    // 16-byte header
    NSUInteger cursor = 16;

    NSMutableArray *offsets = [NSMutableArray arrayWithCapacity:pixels.count];

    if ([mode isEqualToString:@"continuous"]) {
        NSData *sharedMz = pixels[0].mzArray;
        [ibdBuf appendData:sharedMz];
        NSUInteger mzOffset = cursor;
        NSUInteger mzLen = pixels[0].mzCount;
        cursor += sharedMz.length;

        for (NSUInteger i = 0; i < pixels.count; i++) {
            TTIOImzMLPixelSpectrum *p = pixels[i];
            if (p.mzArray.length != sharedMz.length
                || memcmp(p.mzArray.bytes, sharedMz.bytes, sharedMz.length) != 0) {
                if (error) *error = MakeError(4,
                    @"continuous-mode imzML requires all pixels to share the same m/z axis; "
                    @"pixel %lu (x=%ld, y=%ld) differs",
                    (unsigned long)i, (long)p.x, (long)p.y);
                return nil;
            }
            NSData *inten = p.intensityArray;
            [ibdBuf appendData:inten];
            NSUInteger intOffset = cursor;
            NSUInteger intLen = inten.length / sizeof(double);
            cursor += inten.length;
            [offsets addObject:@[@(mzOffset), @(mzLen), @(intOffset), @(intLen)]];
        }
    } else {
        for (TTIOImzMLPixelSpectrum *p in pixels) {
            if (p.mzArray.length != p.intensityArray.length) {
                if (error) *error = MakeError(5,
                    @"processed-mode pixel (x=%ld, y=%ld): mz and intensity arrays"
                    @" must be the same length", (long)p.x, (long)p.y);
                return nil;
            }
            NSUInteger mzLen = p.mzCount;
            [ibdBuf appendData:p.mzArray];
            NSUInteger mzOffset = cursor;
            cursor += p.mzArray.length;

            [ibdBuf appendData:p.intensityArray];
            NSUInteger intOffset = cursor;
            NSUInteger intLen = p.intensityArray.length / sizeof(double);
            cursor += p.intensityArray.length;
            [offsets addObject:@[@(mzOffset), @(mzLen), @(intOffset), @(intLen)]];
        }
    }

    NSError *ioErr = nil;
    BOOL ok = [ibdBuf writeToFile:ibd options:NSDataWritingAtomic error:&ioErr];
    if (!ok) { if (error) *error = ioErr; return nil; }

    NSString *ibdSha1 = SHA1Hex(ibdBuf);

    // ── .imzML XML ─────────────────────────────────────────────
    NSMutableString *xml = [NSMutableString stringWithCapacity:8192];
    NSString *modeAcc = [mode isEqualToString:@"continuous"] ? @"IMS:1000030" : @"IMS:1000031";
    NSString *modeName = [mode isEqualToString:@"continuous"] ? @"continuous" : @"processed";

    [xml appendString:@"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"];
    [xml appendString:@"<mzML xmlns=\"http://psi.hupo.org/ms/mzml\""
                      @" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
                      @" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml"
                      @" http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\""
                      @" version=\"1.1\">\n"];
    [xml appendString:@"  <cvList count=\"3\">\n"];
    [xml appendString:@"    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\""
                      @" version=\"4.1.0\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n"];
    [xml appendString:@"    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"2020-03-10\""
                      @" URI=\"http://ontologies.berkeleybop.org/uo.obo\"/>\n"];
    [xml appendString:@"    <cv id=\"IMS\" fullName=\"Mass Spectrometry Imaging Ontology\""
                      @" version=\"1.1.0\" URI=\"https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo\"/>\n"];
    [xml appendString:@"  </cvList>\n"];
    [xml appendString:@"  <fileDescription>\n"];
    [xml appendString:@"    <fileContent>\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n"];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000080\""
                      @" name=\"universally unique identifier\" value=\"%@\"/>\n", uuid];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000091\""
                      @" name=\"ibd SHA-1\" value=\"%@\"/>\n", ibdSha1];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"%@\" name=\"%@\" value=\"\"/>\n",
                      modeAcc, modeName];
    [xml appendString:@"    </fileContent>\n"];
    [xml appendString:@"  </fileDescription>\n"];

    // referenceableParamGroups for pyimzml + MSIqr compatibility.
    [xml appendString:@"  <referenceableParamGroupList count=\"2\">\n"];
    [xml appendString:@"    <referenceableParamGroup id=\"mzArray\">\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\""
                      @" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"];
    [xml appendString:@"    </referenceableParamGroup>\n"];
    [xml appendString:@"    <referenceableParamGroup id=\"intensityArray\">\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\""
                      @" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n"];
    [xml appendString:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"];
    [xml appendString:@"    </referenceableParamGroup>\n"];
    [xml appendString:@"  </referenceableParamGroupList>\n"];

    [xml appendString:@"  <softwareList count=\"1\">\n"];
    [xml appendString:@"    <software id=\"ttio\" version=\"0.9.0\">\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000799\""
                      @" name=\"custom unreleased software tool\" value=\"ttio\"/>\n"];
    [xml appendString:@"    </software>\n"];
    [xml appendString:@"  </softwareList>\n"];

    [xml appendString:@"  <scanSettingsList count=\"1\">\n"];
    [xml appendString:@"    <scanSettings id=\"scansettings1\">\n"];
    [xml appendFormat:@"      <userParam name=\"scan pattern\" value=\"%@\" type=\"xsd:string\"/>\n",
                      XmlEscape(scanPattern)];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000040\""
                      @" name=\"linescan sequence\" value=\"%@\"/>\n",
                      XmlEscape(scanPattern)];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000042\""
                      @" name=\"max count of pixels x\" value=\"%ld\"/>\n", (long)gridMaxX];
    [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000043\""
                      @" name=\"max count of pixels y\" value=\"%ld\"/>\n", (long)gridMaxY];
    if (pixelSizeX > 0.0) {
        [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000046\""
                          @" name=\"pixel size (x)\" value=\"%g\""
                          @" unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n",
                          pixelSizeX];
    }
    if (pixelSizeY > 0.0) {
        [xml appendFormat:@"      <cvParam cvRef=\"IMS\" accession=\"IMS:1000047\""
                          @" name=\"pixel size y\" value=\"%g\""
                          @" unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n",
                          pixelSizeY];
    }
    [xml appendString:@"    </scanSettings>\n"];
    [xml appendString:@"  </scanSettingsList>\n"];

    [xml appendString:@"  <instrumentConfigurationList count=\"1\">\n"];
    [xml appendString:@"    <instrumentConfiguration id=\"IC1\">\n"];
    [xml appendString:@"      <cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n"];
    [xml appendString:@"    </instrumentConfiguration>\n"];
    [xml appendString:@"  </instrumentConfigurationList>\n"];

    [xml appendString:@"  <dataProcessingList count=\"1\">\n"];
    [xml appendString:@"    <dataProcessing id=\"dp_export\">\n"];
    [xml appendString:@"      <processingMethod order=\"0\" softwareRef=\"ttio\">\n"];
    [xml appendString:@"        <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\"/>\n"];
    [xml appendString:@"      </processingMethod>\n"];
    [xml appendString:@"    </dataProcessing>\n"];
    [xml appendString:@"  </dataProcessingList>\n"];

    [xml appendString:@"  <run id=\"ttio_imzml_export\" defaultInstrumentConfigurationRef=\"IC1\">\n"];
    [xml appendFormat:@"    <spectrumList count=\"%lu\" defaultDataProcessingRef=\"dp_export\">\n",
                      (unsigned long)pixels.count];

    for (NSUInteger i = 0; i < pixels.count; i++) {
        TTIOImzMLPixelSpectrum *p = pixels[i];
        NSArray *o = offsets[i];
        NSUInteger mzOff = [o[0] unsignedIntegerValue];
        NSUInteger mzLen = [o[1] unsignedIntegerValue];
        NSUInteger inOff = [o[2] unsignedIntegerValue];
        NSUInteger inLen = [o[3] unsignedIntegerValue];
        NSUInteger mzEnc = mzLen * 8, inEnc = inLen * 8;

        [xml appendFormat:@"      <spectrum id=\"Scan=%lu\" index=\"%lu\" defaultArrayLength=\"0\">\n",
                          (unsigned long)(i + 1), (unsigned long)i];
        [xml appendString:@"        <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n"];
        [xml appendString:@"        <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>\n"];
        [xml appendString:@"        <scanList count=\"1\">\n"];
        [xml appendString:@"          <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n"];
        [xml appendString:@"          <scan instrumentConfigurationRef=\"IC1\">\n"];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000050\""
                          @" name=\"position x\" value=\"%ld\"/>\n", (long)p.x];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000051\""
                          @" name=\"position y\" value=\"%ld\"/>\n", (long)p.y];
        if (p.z != 1) {
            [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000052\""
                              @" name=\"position z\" value=\"%ld\"/>\n", (long)p.z];
        }
        [xml appendString:@"          </scan>\n"];
        [xml appendString:@"        </scanList>\n"];
        [xml appendString:@"        <binaryDataArrayList count=\"2\">\n"];

        // m/z array
        [xml appendString:@"          <binaryDataArray encodedLength=\"0\">\n"];
        [xml appendString:@"            <referenceableParamGroupRef ref=\"mzArray\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\""
                          @" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000102\""
                          @" name=\"external offset\" value=\"%lu\"/>\n", (unsigned long)mzOff];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000103\""
                          @" name=\"external array length\" value=\"%lu\"/>\n", (unsigned long)mzLen];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000104\""
                          @" name=\"external encoded length\" value=\"%lu\"/>\n", (unsigned long)mzEnc];
        [xml appendString:@"            <binary/>\n"];
        [xml appendString:@"          </binaryDataArray>\n"];

        // intensity array
        [xml appendString:@"          <binaryDataArray encodedLength=\"0\">\n"];
        [xml appendString:@"            <referenceableParamGroupRef ref=\"intensityArray\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\""
                          @" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n"];
        [xml appendString:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000102\""
                          @" name=\"external offset\" value=\"%lu\"/>\n", (unsigned long)inOff];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000103\""
                          @" name=\"external array length\" value=\"%lu\"/>\n", (unsigned long)inLen];
        [xml appendFormat:@"            <cvParam cvRef=\"IMS\" accession=\"IMS:1000104\""
                          @" name=\"external encoded length\" value=\"%lu\"/>\n", (unsigned long)inEnc];
        [xml appendString:@"            <binary/>\n"];
        [xml appendString:@"          </binaryDataArray>\n"];
        [xml appendString:@"        </binaryDataArrayList>\n"];
        [xml appendString:@"      </spectrum>\n"];

        // Per-N progress fire. Total = pixels.count.
        if (((i + 1) % TTIOImzMLWriterProgressIntervalPixels) == 0) {
            progress((int64_t)(i + 1), (int64_t)pixels.count);
        }
    }
    // Final fire.
    progress((int64_t)pixels.count, (int64_t)pixels.count);

    [xml appendString:@"    </spectrumList>\n"];
    [xml appendString:@"  </run>\n"];
    [xml appendString:@"</mzML>\n"];

    NSData *xmlBytes = [xml dataUsingEncoding:NSUTF8StringEncoding];
    ok = [xmlBytes writeToFile:imzmlPath options:NSDataWritingAtomic error:&ioErr];
    if (!ok) { if (error) *error = ioErr; return nil; }

    return [[TTIOImzMLWriteResult alloc] initInternal:imzmlPath ibd:ibd
                                                  uuid:uuid mode:mode
                                               nPixels:pixels.count];
}

#pragma mark - Pixel runs (M102)

/* NSNumber / numeric NSString parameter values. */
static NSInteger ParamInteger(id v)
{
    if ([v isKindOfClass:[NSNumber class]] || [v isKindOfClass:[NSString class]])
        return [v integerValue];
    return 0;
}

static double ParamDouble(id v)
{
    if ([v isKindOfClass:[NSNumber class]] || [v isKindOfClass:[NSString class]])
        return [v doubleValue];
    return 0.0;
}

static NSString *ParamString(id v)
{
    if ([v isKindOfClass:[NSString class]] && [(NSString *)v length]) return v;
    return nil;
}

+ (nullable TTIOImzMLWriteResult *)writeRun:(TTIOAcquisitionRun *)run
                          datasetProvenance:(nullable NSArray<TTIOProvenanceRecord *> *)datasetProvenance
                                toImzMLPath:(NSString *)imzmlPath
                                    ibdPath:(nullable NSString *)ibdPath
                                      error:(NSError **)error
{
    NSData *coords = [TTIOImzMLReader pixelCoordinatesForRun:run
                                           datasetProvenance:datasetProvenance];
    if (!coords) {
        if (error) *error = MakeError(5,
            @"run '%@' has no pixel coordinates to export as imzML",
            run.name ?: @"");
        return nil;
    }

    // Run records win over dataset records, later records over
    // earlier (Python write_from_run).
    NSMutableArray<TTIOProvenanceRecord *> *records = [NSMutableArray array];
    [records addObjectsFromArray:datasetProvenance ?: @[]];
    [records addObjectsFromArray:[run provenanceChain] ?: @[]];
    NSMutableDictionary<NSString *, id> *params = [NSMutableDictionary dictionary];
    for (TTIOProvenanceRecord *r in [records reverseObjectEnumerator]) {
        [r.parameters enumerateKeysAndObjectsUsingBlock:^(NSString *k, id v, BOOL *stop) {
            (void)stop;
            if ([k hasPrefix:@"imzml_"] && params[k] == nil) params[k] = v;
        }];
    }

    const int32_t *xyz = (const int32_t *)coords.bytes;
    NSUInteger n = coords.length / (3 * sizeof(int32_t));
    NSMutableArray<TTIOImzMLPixelSpectrum *> *pixels =
        [NSMutableArray arrayWithCapacity:n];
    // The pool drains per pixel; a failure's NSError is kept in a
    // strong local so it outlives the pool.
    NSError *loopErr = nil;
    BOOL loopOk = YES;
    for (NSUInteger i = 0; i < n && loopOk; i++) {
        @autoreleasepool {
            NSError *e = nil;
            TTIOSpectrum *spec = [run spectrumAtIndex:i error:&e];
            TTIOSignalArray *mzA = spec.signalArrays[@"mz"];
            TTIOSignalArray *inA = spec.signalArrays[@"intensity"];
            TTIOImzMLPixelSpectrum *px = nil;
            if (spec && (!mzA || !inA)) {
                e = MakeError(6, @"run '%@' spectrum %lu lacks mz/intensity arrays",
                              run.name ?: @"", (unsigned long)i);
            } else if (spec) {
                px = [[TTIOImzMLPixelSpectrum alloc] initWithX:xyz[i * 3]
                                                             y:xyz[i * 3 + 1]
                                                             z:xyz[i * 3 + 2]
                                                       mzArray:[[mzA float64Buffer] copy]
                                                intensityArray:[[inA float64Buffer] copy]
                                                         error:&e];
            }
            if (px) {
                [pixels addObject:px];
            } else {
                loopErr = e;
                loopOk = NO;
            }
        }
    }
    if (!loopOk) {
        if (error) *error = loopErr;
        return nil;
    }

    NSString *mode = ParamString(params[@"imzml_mode"]);
    if (![mode isEqualToString:@"continuous"] && ![mode isEqualToString:@"processed"]) {
        BOOL shared = YES;
        NSData *first = pixels.firstObject.mzArray;
        for (TTIOImzMLPixelSpectrum *p in pixels) {
            if (![p.mzArray isEqualToData:first]) { shared = NO; break; }
        }
        mode = shared ? @"continuous" : @"processed";
    }

    return [self writePixels:pixels
                 toImzMLPath:imzmlPath
                     ibdPath:ibdPath
                        mode:mode
                    gridMaxX:ParamInteger(params[@"imzml_grid_max_x"])
                    gridMaxY:ParamInteger(params[@"imzml_grid_max_y"])
                    gridMaxZ:ParamInteger(params[@"imzml_grid_max_z"])
                  pixelSizeX:ParamDouble(params[@"imzml_pixel_size_x"])
                  pixelSizeY:ParamDouble(params[@"imzml_pixel_size_y"])
                 scanPattern:ParamString(params[@"imzml_scan_pattern"]) ?: @"flyback"
                     uuidHex:ParamString(params[@"imzml_uuid_hex"])
                    progress:nil
                       error:error];
}

+ (nullable TTIOImzMLWriteResult *)writeFromImport:(TTIOImzMLImport *)import
                                         toImzMLPath:(NSString *)imzmlPath
                                             ibdPath:(nullable NSString *)ibdPath
                                               error:(NSError **)error
{
    return [self writeFromImport:import
                     toImzMLPath:imzmlPath
                         ibdPath:ibdPath
                        progress:nil
                           error:error];
}

+ (nullable TTIOImzMLWriteResult *)writeFromImport:(TTIOImzMLImport *)import
                                         toImzMLPath:(NSString *)imzmlPath
                                             ibdPath:(nullable NSString *)ibdPath
                                            progress:(TTIOProgressBlock)progress
                                               error:(NSError **)error
{
    return [self writePixels:import.spectra
                 toImzMLPath:imzmlPath
                     ibdPath:ibdPath
                        mode:import.mode
                    gridMaxX:import.gridMaxX
                    gridMaxY:import.gridMaxY
                    gridMaxZ:import.gridMaxZ
                  pixelSizeX:import.pixelSizeX
                  pixelSizeY:import.pixelSizeY
                 scanPattern:import.scanPattern
                     uuidHex:import.uuidHex
                    progress:progress
                       error:error];
}

@end
