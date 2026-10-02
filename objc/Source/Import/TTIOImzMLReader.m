/*
 * TTIOImzMLReader.m
 * TTI-O Objective-C Implementation
 *
 * Classes:       TTIOImzMLPixelSpectrum, TTIOImzMLImport,
 *                TTIOImzMLReader
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Import/TTIOImzMLReader.h
 *
 * imzML + .ibd importer for mass-spectrometry imaging data. Parses
 * the .imzML metadata (continuous or processed mode), validates the
 * 16-byte UUID prefix in .ibd, and reads each pixel's m/z +
 * intensity arrays at the offsets named in the metadata.
 *
 * Licensed under the Apache License, Version 2.0.
 * See LICENSE-IMPORT-EXPORT in the repository root.
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#import "TTIOImzMLReader.h"
#import "Import/TTIOXMLStreamParser.h"
#import "Import/TTIOImportedDataset.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Dataset/TTIOWrittenRun.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOSpectrumIndex.h"
#import "ValueClasses/TTIOEnums.h"
#import <time.h>

#import <unistd.h>

NSString *const TTIOImzMLReaderErrorDomain = @"TTIOImzMLReaderErrorDomain";
NSString *const TTIOImzMLPixelRunName = @"imzml_pixels";
NSString *const TTIOImzMLLegacyCoordinatesParameter = @"imzml_pixel_coordinates_csv";

// Mirrors Java ImzMLReader.PROGRESS_INTERVAL_PIXELS (100).
const NSUInteger TTIOImzMLReaderProgressIntervalPixels = 100;

#pragma mark - Pixel spectrum

@implementation TTIOImzMLPixelSpectrum
{
    NSData *_mzArray;
    NSData *_intensityArray;
}

- (instancetype)initWithX:(NSInteger)x
                        y:(NSInteger)y
                        z:(NSInteger)z
                       mz:(NSData *)mz
                intensity:(NSData *)intensity
{
    if ((self = [super init])) {
        _x = x;
        _y = y;
        _z = z;
        _mzArray = [mz copy];
        _intensityArray = [intensity copy];
    }
    return self;
}

- (NSData *)mzArray { return _mzArray; }
- (NSData *)intensityArray { return _intensityArray; }
- (NSUInteger)mzCount { return _mzArray.length / sizeof(double); }

- (nullable instancetype)initWithX:(NSInteger)x
                                  y:(NSInteger)y
                                  z:(NSInteger)z
                            mzArray:(NSData *)mzArray
                     intensityArray:(NSData *)intensityArray
                              error:(NSError **)error
{
    if (!mzArray || mzArray.length == 0 || mzArray.length % 8 != 0) {
        if (error) *error = [NSError errorWithDomain:TTIOImzMLReaderErrorDomain
            code:TTIOImzMLReaderErrorMissingMetadata userInfo:@{
            NSLocalizedDescriptionKey: @"mzArray must be non-empty and a multiple of 8 bytes"
        }];
        return nil;
    }
    if (!intensityArray || intensityArray.length == 0
        || intensityArray.length != mzArray.length) {
        if (error) *error = [NSError errorWithDomain:TTIOImzMLReaderErrorDomain
            code:TTIOImzMLReaderErrorMissingMetadata userInfo:@{
            NSLocalizedDescriptionKey: @"intensityArray must be non-empty and match mzArray length"
        }];
        return nil;
    }
    return [self initWithX:x y:y z:z mz:mzArray intensity:intensityArray];
}

@end

#pragma mark - Import value object

@implementation TTIOImzMLImport

- (instancetype)initWithMode:(NSString *)mode
                     uuidHex:(NSString *)uuidHex
                    gridMaxX:(NSInteger)gridMaxX
                    gridMaxY:(NSInteger)gridMaxY
                    gridMaxZ:(NSInteger)gridMaxZ
                  pixelSizeX:(double)pixelSizeX
                  pixelSizeY:(double)pixelSizeY
                 scanPattern:(NSString *)scanPattern
                     spectra:(NSArray<TTIOImzMLPixelSpectrum *> *)spectra
                 sourceImzML:(NSString *)sourceImzML
                   sourceIbd:(NSString *)sourceIbd
{
    if ((self = [super init])) {
        _mode = [mode copy];
        _uuidHex = [uuidHex copy];
        _gridMaxX = gridMaxX;
        _gridMaxY = gridMaxY;
        _gridMaxZ = gridMaxZ;
        _pixelSizeX = pixelSizeX;
        _pixelSizeY = pixelSizeY;
        _scanPattern = [scanPattern copy] ?: @"";
        _spectra = [spectra copy];
        _sourceImzML = [sourceImzML copy] ?: @"";
        _sourceIbd = [sourceIbd copy] ?: @"";
    }
    return self;
}

@end

#pragma mark - Reader implementation

@interface TTIOImzMLReaderState : NSObject
@property (nonatomic, copy) NSString *mode;
@property (nonatomic, copy) NSString *uuidHex;
@property (nonatomic) NSInteger gridMaxX;
@property (nonatomic) NSInteger gridMaxY;
@property (nonatomic) NSInteger gridMaxZ;
@property (nonatomic) double pixelSizeX;
@property (nonatomic) double pixelSizeY;
@property (nonatomic, copy) NSString *scanPattern;
@property (nonatomic, strong) NSMutableArray *stubs; // array of NSMutableDictionary
@end

@implementation TTIOImzMLReaderState
- (instancetype)init {
    if ((self = [super init])) {
        _mode = @"";
        _uuidHex = @"";
        _gridMaxZ = 1;
        _scanPattern = @"";
        _stubs = [NSMutableArray array];
    }
    return self;
}
@end

@interface TTIOImzMLReader () <NSXMLParserDelegate>
@end

@implementation TTIOImzMLReader
{
    TTIOImzMLReaderState *_state;
    NSMutableDictionary *_currentStub;
    BOOL _inSpectrum;
    BOOL _inBinaryArray;
    BOOL _inScan;
    NSString *_currentArrayKind;     // "mz" / "intensity" / @""
    NSString *_currentArrayPrecision; // "32" / "64"
}

#pragma mark - CV term constants

// imzML storage mode accessions. Only the IMS-namespaced forms are
// real: MS:1000030 = "vendor processing software", MS:1000031 =
// "instrument model" — completely unrelated terms that would
// false-positive on every well-formed mzML/imzML file otherwise.
static NSString *const kCVContinuous30 = @"IMS:1000030";
static NSString *const kCVProcessed31  = @"IMS:1000031";
// Canonical IMS accessions (real-world imzML 1.1, pyimzML test corpus,
// TTIO writer output v0.9+).
static NSString *const kCVUUID         = @"IMS:1000080";
static NSString *const kCVMaxX         = @"IMS:1000042";
static NSString *const kCVMaxY         = @"IMS:1000043";
// Legacy TTIO synthetic-fixture accessions (pre-v0.9): kept for
// backward-compat with any old test files. Importer handler probes
// both canonical + legacy.
static NSString *const kCVUUIDLegacy   = @"IMS:1000042";
static NSString *const kCVMaxXLegacy   = @"IMS:1000003";
static NSString *const kCVMaxYLegacy   = @"IMS:1000004";
static NSString *const kCVMaxZ         = @"IMS:1000005";
static NSString *const kCVPixelSizeX   = @"IMS:1000046";
static NSString *const kCVPixelSizeY   = @"IMS:1000047";
static NSString *const kCVScanPattern1 = @"IMS:1000040";
static NSString *const kCVScanPattern2 = @"IMS:1000048";
static NSString *const kCVPositionX    = @"IMS:1000050";
static NSString *const kCVPositionY    = @"IMS:1000051";
static NSString *const kCVPositionZ    = @"IMS:1000052";
static NSString *const kCVMzArray      = @"MS:1000514";
static NSString *const kCVIntensity    = @"MS:1000515";
static NSString *const kCV64Bit        = @"MS:1000523";
static NSString *const kCV32Bit        = @"MS:1000521";
static NSString *const kCVExtOffset    = @"IMS:1000102";
static NSString *const kCVExtLength    = @"IMS:1000103";
static NSString *const kCVExtEncoded   = @"IMS:1000104";

#pragma mark - Helpers

static NSString *normaliseUUID(NSString *value) {
    NSMutableString *m = [[value lowercaseString] mutableCopy];
    [m replaceOccurrencesOfString:@"{" withString:@"" options:0 range:NSMakeRange(0, m.length)];
    [m replaceOccurrencesOfString:@"}" withString:@"" options:0 range:NSMakeRange(0, m.length)];
    [m replaceOccurrencesOfString:@"-" withString:@"" options:0 range:NSMakeRange(0, m.length)];
    return [m stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
}

+ (NSError *)errorWithCode:(TTIOImzMLReaderErrorCode)code message:(NSString *)message {
    return [NSError errorWithDomain:TTIOImzMLReaderErrorDomain
                               code:code
                           userInfo:@{NSLocalizedDescriptionKey: message}];
}

#pragma mark - Public API

+ (nullable TTIOImzMLImport *)readFromImzMLPath:(NSString *)imzmlPath
                                         ibdPath:(nullable NSString *)ibdPath
                                           error:(NSError **)error
{
    return [self readFromImzMLPath:imzmlPath ibdPath:ibdPath progress:nil error:error];
}

+ (nullable TTIOImzMLImport *)readFromImzMLPath:(NSString *)imzmlPath
                                         ibdPath:(nullable NSString *)ibdPath
                                        progress:(nullable TTIOProgressBlock)progress
                                           error:(NSError **)error
{
    NSFileManager *fm = [NSFileManager defaultManager];
    if (![fm fileExistsAtPath:imzmlPath]) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingFile
                                        message:[NSString stringWithFormat:@"imzML metadata not found: %@", imzmlPath]];
        return nil;
    }
    NSString *resolvedIbd = ibdPath ?: [[imzmlPath stringByDeletingPathExtension] stringByAppendingPathExtension:@"ibd"];
    if (![fm fileExistsAtPath:resolvedIbd]) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingFile
                                        message:[NSString stringWithFormat:@"imzML binary not found: %@", resolvedIbd]];
        return nil;
    }

    TTIOImzMLReader *reader = [[TTIOImzMLReader alloc] init];
    reader->_state = [[TTIOImzMLReaderState alloc] init];
    reader->_currentArrayKind = @"";
    reader->_currentArrayPrecision = @"64";

    NSError *xmlError = nil;
    if (![TTIOXMLStreamParser parseFileAtPath:imzmlPath
                                     delegate:reader
                                        error:&xmlError]) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorParseFailed
                                        message:[NSString stringWithFormat:@"cannot parse imzML: %@",
                                                 xmlError.localizedDescription ?: @"(unknown)"]];
        return nil;
    }

    if (reader->_state.mode.length == 0) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingMetadata
                                        message:@"no continuous/processed mode CV term found"];
        return nil;
    }
    if (reader->_state.uuidHex.length == 0) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingMetadata
                                        message:@"missing IMS:1000042 universally unique identifier"];
        return nil;
    }
    if (reader->_state.stubs.count == 0) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingMetadata
                                        message:@"no <spectrum> elements parsed"];
        return nil;
    }

    NSDictionary *ibdAttrs = [fm attributesOfItemAtPath:resolvedIbd error:NULL];
    if (!ibdAttrs) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingFile
                                        message:[NSString stringWithFormat:@"cannot read .ibd: %@", resolvedIbd]];
        return nil;
    }
    unsigned long long ibdSize = [ibdAttrs fileSize];
    if (ibdSize < 16) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorBinaryShorterThanUUID
                                        message:[NSString stringWithFormat:@"%@ shorter than the 16-byte UUID header", resolvedIbd]];
        return nil;
    }
    NSFileHandle *ibd = [NSFileHandle fileHandleForReadingAtPath:resolvedIbd];
    if (!ibd) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingFile
                                        message:[NSString stringWithFormat:@"cannot read .ibd: %@", resolvedIbd]];
        return nil;
    }
    NSData *uuidBytes = [ibd readDataOfLength:16];
    if (uuidBytes.length != 16) {
        [ibd closeFile];
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorBinaryShorterThanUUID
                                        message:[NSString stringWithFormat:@"%@ shorter than the 16-byte UUID header", resolvedIbd]];
        return nil;
    }
    NSMutableString *ibdUUIDHex = [NSMutableString stringWithCapacity:32];
    const unsigned char *bytes = uuidBytes.bytes;
    for (NSUInteger i = 0; i < 16; i++) {
        [ibdUUIDHex appendFormat:@"%02x", bytes[i]];
    }
    if (![ibdUUIDHex isEqualToString:reader->_state.uuidHex]) {
        [ibd closeFile];
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorUUIDMismatch
                                        message:[NSString stringWithFormat:@"UUID mismatch: imzML declares %@ but .ibd header is %@", reader->_state.uuidHex, ibdUUIDHex]];
        return nil;
    }

    NSError *materialiseError = nil;
    NSArray<TTIOImzMLPixelSpectrum *> *pixels =
        [reader materialiseSpectraWithIBD:ibd.fileDescriptor
                                   ibdSize:ibdSize
                                   ibdPath:resolvedIbd
                                  progress:progress
                                     error:&materialiseError];
    [ibd closeFile];
    if (!pixels) {
        if (error) *error = materialiseError;
        return nil;
    }

    return [[TTIOImzMLImport alloc] initWithMode:reader->_state.mode
                                          uuidHex:reader->_state.uuidHex
                                         gridMaxX:reader->_state.gridMaxX
                                         gridMaxY:reader->_state.gridMaxY
                                         gridMaxZ:reader->_state.gridMaxZ
                                       pixelSizeX:reader->_state.pixelSizeX
                                       pixelSizeY:reader->_state.pixelSizeY
                                      scanPattern:reader->_state.scanPattern
                                          spectra:pixels
                                      sourceImzML:imzmlPath
                                        sourceIbd:resolvedIbd];
}

#pragma mark - NSXMLParserDelegate

- (void)parser:(NSXMLParser *)parser
  didStartElement:(NSString *)elementName
    namespaceURI:(NSString *)namespaceURI
   qualifiedName:(NSString *)qName
      attributes:(NSDictionary<NSString *, NSString *> *)attrs
{
    if ([elementName isEqualToString:@"spectrum"]) {
        _currentStub = [NSMutableDictionary dictionary];
        _currentStub[@"x"] = @0; _currentStub[@"y"] = @0; _currentStub[@"z"] = @1;
        _currentStub[@"mz_offset"] = @(-1);  _currentStub[@"mz_length"] = @0;  _currentStub[@"mz_precision"] = @"64";
        _currentStub[@"int_offset"] = @(-1); _currentStub[@"int_length"] = @0; _currentStub[@"int_precision"] = @"64";
        _inSpectrum = YES;
    } else if ([elementName isEqualToString:@"binaryDataArray"]) {
        _inBinaryArray = YES;
        _currentArrayKind = @"";
    } else if ([elementName isEqualToString:@"scan"]) {
        _inScan = YES;
    } else if ([elementName isEqualToString:@"cvParam"]) {
        [self handleCVParam:attrs];
    }
}

- (void)parser:(NSXMLParser *)parser
   didEndElement:(NSString *)elementName
    namespaceURI:(NSString *)namespaceURI
   qualifiedName:(NSString *)qName
{
    if ([elementName isEqualToString:@"spectrum"]) {
        if (_currentStub) [_state.stubs addObject:_currentStub];
        _currentStub = nil;
        _inSpectrum = NO;
    } else if ([elementName isEqualToString:@"binaryDataArray"]) {
        _inBinaryArray = NO;
        _currentArrayKind = @"";
    } else if ([elementName isEqualToString:@"scan"]) {
        _inScan = NO;
    }
}

- (void)handleCVParam:(NSDictionary<NSString *, NSString *> *)attrs {
    NSString *acc = attrs[@"accession"] ?: @"";
    NSString *value = attrs[@"value"] ?: @"";

    if ([acc isEqualToString:kCVContinuous30]) {
        _state.mode = @"continuous";
    } else if ([acc isEqualToString:kCVProcessed31]) {
        _state.mode = @"processed";
    } else if ([acc isEqualToString:kCVUUID] && value.length > 0) {
        _state.uuidHex = normaliseUUID(value);
    } else if ([acc isEqualToString:kCVMaxX] && value.length > 0) {
        _state.gridMaxX = value.integerValue;
    } else if ([acc isEqualToString:kCVMaxY] && value.length > 0) {
        _state.gridMaxY = value.integerValue;
    } else if ([acc isEqualToString:kCVUUIDLegacy] && value.length > 0
               && _state.uuidHex.length == 0) {
        // Legacy TTIO pre-v0.9 fallback: only consume IMS:1000042 as
        // UUID when IMS:1000080 hasn't appeared yet AND the value
        // normalises to a 32-hex-char UUID. Real imzML uses
        // IMS:1000042 for "max count of pixels x" with an integer
        // value, never a UUID.
        NSString *cand = normaliseUUID(value);
        if (cand.length == 32) _state.uuidHex = cand;
    } else if ([acc isEqualToString:kCVMaxXLegacy] && value.length > 0) {
        _state.gridMaxX = value.integerValue;
    } else if ([acc isEqualToString:kCVMaxYLegacy] && value.length > 0) {
        _state.gridMaxY = value.integerValue;
    } else if ([acc isEqualToString:kCVMaxZ] && value.length > 0) {
        _state.gridMaxZ = value.integerValue;
    } else if ([acc isEqualToString:kCVPixelSizeX] && value.length > 0) {
        _state.pixelSizeX = value.doubleValue;
    } else if ([acc isEqualToString:kCVPixelSizeY] && value.length > 0) {
        _state.pixelSizeY = value.doubleValue;
    } else if (([acc isEqualToString:kCVScanPattern1] || [acc isEqualToString:kCVScanPattern2]) && value.length > 0) {
        if (_state.scanPattern.length == 0) _state.scanPattern = value;
    } else if (_inScan && _currentStub) {
        if ([acc isEqualToString:kCVPositionX] && value.length > 0) {
            _currentStub[@"x"] = @(value.integerValue);
        } else if ([acc isEqualToString:kCVPositionY] && value.length > 0) {
            _currentStub[@"y"] = @(value.integerValue);
        } else if ([acc isEqualToString:kCVPositionZ] && value.length > 0) {
            _currentStub[@"z"] = @(value.integerValue);
        }
    } else if (_inBinaryArray && _currentStub) {
        if ([acc isEqualToString:kCVMzArray]) {
            _currentArrayKind = @"mz";
        } else if ([acc isEqualToString:kCVIntensity]) {
            _currentArrayKind = @"intensity";
        } else if ([acc isEqualToString:kCV64Bit]) {
            if ([_currentArrayKind isEqualToString:@"mz"]) _currentStub[@"mz_precision"] = @"64";
            else if ([_currentArrayKind isEqualToString:@"intensity"]) _currentStub[@"int_precision"] = @"64";
        } else if ([acc isEqualToString:kCV32Bit]) {
            if ([_currentArrayKind isEqualToString:@"mz"]) _currentStub[@"mz_precision"] = @"32";
            else if ([_currentArrayKind isEqualToString:@"intensity"]) _currentStub[@"int_precision"] = @"32";
        } else if ([acc isEqualToString:kCVExtOffset] && value.length > 0) {
            if ([_currentArrayKind isEqualToString:@"mz"]) _currentStub[@"mz_offset"] = @(value.longLongValue);
            else if ([_currentArrayKind isEqualToString:@"intensity"]) _currentStub[@"int_offset"] = @(value.longLongValue);
        } else if ([acc isEqualToString:kCVExtLength] && value.length > 0) {
            if ([_currentArrayKind isEqualToString:@"mz"]) _currentStub[@"mz_length"] = @(value.longLongValue);
            else if ([_currentArrayKind isEqualToString:@"intensity"]) _currentStub[@"int_length"] = @(value.longLongValue);
        }
    }
}

#pragma mark - Binary materialisation

- (NSArray<TTIOImzMLPixelSpectrum *> *)materialiseSpectraWithIBD:(int)ibd
                                                          ibdSize:(unsigned long long)ibdSize
                                                          ibdPath:(NSString *)ibdPath
                                                         progress:(TTIOProgressBlock)progress
                                                            error:(NSError **)error
{
    NSMutableArray<TTIOImzMLPixelSpectrum *> *pixels = [NSMutableArray array];
    NSData *sharedMz = nil;
    BOOL continuous = [_state.mode isEqualToString:@"continuous"];
    NSUInteger total = _state.stubs.count;

    for (NSDictionary *stub in _state.stubs) {
        NSData *mzData = [self readArrayFromIBD:ibd
                                          offset:[stub[@"mz_offset"] longLongValue]
                                          length:[stub[@"mz_length"] longLongValue]
                                       precision:stub[@"mz_precision"]
                                         ibdSize:ibdSize
                                         ibdPath:ibdPath
                                           label:@"m/z"
                                           error:error];
        if (!mzData) return nil;
        NSData *intData = [self readArrayFromIBD:ibd
                                           offset:[stub[@"int_offset"] longLongValue]
                                           length:[stub[@"int_length"] longLongValue]
                                        precision:stub[@"int_precision"]
                                          ibdSize:ibdSize
                                          ibdPath:ibdPath
                                            label:@"intensity"
                                            error:error];
        if (!intData) return nil;

        NSData *effectiveMz;
        if (continuous) {
            if (!sharedMz) sharedMz = mzData;
            effectiveMz = sharedMz;
        } else {
            effectiveMz = mzData;
        }
        if (effectiveMz.length / sizeof(double) != intData.length / sizeof(double)) {
            if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                    message:[NSString stringWithFormat:@"%@: mz/intensity size mismatch", ibdPath]];
            return nil;
        }
        TTIOImzMLPixelSpectrum *pixel = [[TTIOImzMLPixelSpectrum alloc] initWithX:[stub[@"x"] integerValue]
                                                                                y:[stub[@"y"] integerValue]
                                                                                z:[stub[@"z"] integerValue]
                                                                               mz:effectiveMz
                                                                        intensity:intData];
        [pixels addObject:pixel];

        // Per-N progress fire INSIDE the materialisation loop. imzML
        // total pixel count IS known (from <spectrum> stubs), so emit
        // a real total rather than -1.
        if (progress &&
            (pixels.count % TTIOImzMLReaderProgressIntervalPixels) == 0) {
            progress((int64_t)pixels.count, (int64_t)total);
        }
    }
    // Final fire.
    if (progress) progress((int64_t)pixels.count, (int64_t)total);
    return pixels;
}

- (NSData *)readArrayFromIBD:(int)ibd
                       offset:(long long)offset
                       length:(long long)length
                    precision:(NSString *)precision
                      ibdSize:(unsigned long long)ibdSize
                      ibdPath:(NSString *)ibdPath
                        label:(NSString *)label
                        error:(NSError **)error
{
    if (offset < 0 || length < 0) {
        if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                message:[NSString stringWithFormat:@"%@: negative offset/length for %@ array", ibdPath, label]];
        return nil;
    }
    if (length == 0) return [NSData data];
    NSUInteger bytesPer = [precision isEqualToString:@"64"] ? 8 : 4;
    NSUInteger nbytes = (NSUInteger)length * bytesPer;
    if ((unsigned long long)offset + (unsigned long long)nbytes > ibdSize) {
        if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                message:[NSString stringWithFormat:@"%@: %@ array reads past end of file (offset=%lld, bytes=%lu, size=%llu)",
                                                         ibdPath, label, offset, (unsigned long)nbytes, ibdSize]];
        return nil;
    }
    void *raw = malloc(nbytes);
    if (!raw) {
        if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                message:[NSString stringWithFormat:@"%@: out of memory for a %lu-byte %@ array",
                                                         ibdPath, (unsigned long)nbytes, label]];
        return nil;
    }
    ssize_t got = pread(ibd, raw, nbytes, (off_t)offset);
    if (got != (ssize_t)nbytes) {
        free(raw);
        if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                message:[NSString stringWithFormat:@"%@: %@ array is short (offset=%lld, wanted=%lu, got=%ld)",
                                                         ibdPath, label, offset, (unsigned long)nbytes, (long)got]];
        return nil;
    }
    /* -dataWithBytesNoCopy: takes the malloc'd buffer and yields an
     * immutable NSData, so -copy on it is identity. The pixel
     * initialiser copies what it is handed, and continuous mode
     * depends on every pixel aliasing one m/z object. */
    if (bytesPer == 8) {
        return [NSData dataWithBytesNoCopy:raw length:nbytes freeWhenDone:YES];
    }
    // 32-bit -> promote to 64-bit float NSData.
    NSUInteger n = (NSUInteger)length;
    double *dst = malloc(n * sizeof(double));
    if (!dst) {
        free(raw);
        if (error) *error = [[self class] errorWithCode:TTIOImzMLReaderErrorOffsetOverflow
                                                message:[NSString stringWithFormat:@"%@: out of memory promoting a %lu-value %@ array",
                                                         ibdPath, (unsigned long)n, label]];
        return nil;
    }
    const float *src = (const float *)raw;
    for (NSUInteger i = 0; i < n; i++) dst[i] = (double)src[i];
    free(raw);
    return [NSData dataWithBytesNoCopy:dst length:n * sizeof(double) freeWhenDone:YES];
}

#pragma mark - Pixel runs (M102)

/* Parses "x,y,z;x,y,z;..." into int32_t[count * 3]; nil unless the
 * value is a string listing exactly `count` integer triples. */
static NSData *parseLegacyCoordinates(id value, NSUInteger count)
{
    if (![value isKindOfClass:[NSString class]] || [(NSString *)value length] == 0) return nil;
    NSArray<NSString *> *triples = [(NSString *)value componentsSeparatedByString:@";"];
    if (triples.count != count) return nil;
    NSMutableData *out = [NSMutableData dataWithLength:count * 3 * sizeof(int32_t)];
    int32_t *p = (int32_t *)out.mutableBytes;
    NSCharacterSet *ws = [NSCharacterSet whitespaceAndNewlineCharacterSet];
    for (NSUInteger i = 0; i < count; i++) {
        NSArray<NSString *> *parts = [triples[i] componentsSeparatedByString:@","];
        if (parts.count != 3) return nil;
        for (NSUInteger k = 0; k < 3; k++) {
            NSString *t = [parts[k] stringByTrimmingCharactersInSet:ws];
            NSScanner *sc = [NSScanner scannerWithString:t];
            long long v = 0;
            if (t.length == 0 || ![sc scanLongLong:&v] || ![sc isAtEnd]) return nil;
            if (v < INT32_MIN || v > INT32_MAX) return nil;
            p[i * 3 + k] = (int32_t)v;
        }
    }
    return out;
}

+ (nullable NSData *)pixelCoordinatesForRun:(TTIOAcquisitionRun *)run
                          datasetProvenance:(nullable NSArray<TTIOProvenanceRecord *> *)datasetProvenance
{
    TTIOSpectrumIndex *idx = run.spectrumIndex;
    NSUInteger n = idx.count;
    if (idx.hasPixelCoordinates) {
        NSMutableData *out = [NSMutableData dataWithLength:n * 3 * sizeof(int32_t)];
        int32_t *p = (int32_t *)out.mutableBytes;
        for (NSUInteger i = 0; i < n; i++) {
            p[i * 3]     = [idx pixelXAt:i];
            p[i * 3 + 1] = [idx pixelYAt:i];
            p[i * 3 + 2] = [idx pixelZAt:i];
        }
        return out;
    }
    NSMutableArray<TTIOProvenanceRecord *> *records = [NSMutableArray array];
    [records addObjectsFromArray:[run provenanceChain] ?: @[]];
    [records addObjectsFromArray:datasetProvenance ?: @[]];
    for (TTIOProvenanceRecord *r in records) {
        NSData *c = parseLegacyCoordinates(r.parameters[TTIOImzMLLegacyCoordinatesParameter], n);
        if (c) return c;
    }
    return nil;
}

+ (TTIOProvenanceRecord *)provenanceRecordForImport:(TTIOImzMLImport *)import
{
    NSDictionary *params = @{
        @"imzml_mode":         import.mode ?: @"",
        @"imzml_uuid_hex":     import.uuidHex ?: @"",
        @"imzml_grid_max_x":   @((long long)import.gridMaxX),
        @"imzml_grid_max_y":   @((long long)import.gridMaxY),
        @"imzml_grid_max_z":   @((long long)import.gridMaxZ),
        @"imzml_pixel_size_x": @(import.pixelSizeX),
        @"imzml_pixel_size_y": @(import.pixelSizeY),
        @"imzml_scan_pattern": import.scanPattern ?: @"",
    };
    return [[TTIOProvenanceRecord alloc]
        initWithInputRefs:@[import.sourceImzML ?: @"", import.sourceIbd ?: @""]
                 software:@"ttio imzml importer v0.9"
               parameters:params
               outputRefs:@[]
            timestampUnix:(int64_t)time(NULL)];
}

+ (nullable TTIOWrittenRun *)pixelRunFromImport:(TTIOImzMLImport *)import
                                          error:(NSError **)error
{
    NSArray<TTIOImzMLPixelSpectrum *> *spectra = import.spectra;
    NSUInteger n = spectra.count;
    if (n == 0) {
        if (error) *error = [self errorWithCode:TTIOImzMLReaderErrorMissingMetadata
            message:[NSString stringWithFormat:@"%@: no spectra parsed",
                     import.sourceImzML]];
        return nil;
    }
    NSUInteger total = 0;
    for (TTIOImzMLPixelSpectrum *s in spectra) total += s.mzCount;

    NSMutableData *mz  = [NSMutableData dataWithCapacity:total * sizeof(double)];
    NSMutableData *it  = [NSMutableData dataWithCapacity:total * sizeof(double)];
    NSMutableData *off = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    NSMutableData *len = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *rt  = [NSMutableData dataWithLength:n * sizeof(double)];
    NSMutableData *ml  = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *pol = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *pmz = [NSMutableData dataWithLength:n * sizeof(double)];
    NSMutableData *pc  = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *bp  = [NSMutableData dataWithLength:n * sizeof(double)];
    NSMutableData *px  = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *py  = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *pz  = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    int64_t  *offP = off.mutableBytes;
    uint32_t *lenP = len.mutableBytes;
    int32_t  *mlP  = ml.mutableBytes;
    int32_t  *polP = pol.mutableBytes;
    double   *bpP  = bp.mutableBytes;
    int32_t  *pxP  = px.mutableBytes;
    int32_t  *pyP  = py.mutableBytes;
    int32_t  *pzP  = pz.mutableBytes;

    int64_t cursor = 0;
    for (NSUInteger i = 0; i < n; i++) {
        TTIOImzMLPixelSpectrum *s = spectra[i];
        NSUInteger m = s.mzCount;
        [mz appendBytes:s.mzArray.bytes length:m * sizeof(double)];
        [it appendBytes:s.intensityArray.bytes length:m * sizeof(double)];
        offP[i] = cursor;
        lenP[i] = (uint32_t)m;
        mlP[i]  = 1;
        polP[i] = (int32_t)TTIOPolarityUnknown;
        const double *ip = (const double *)s.intensityArray.bytes;
        double maxI = 0.0;
        for (NSUInteger j = 0; j < m; j++) {
            if (j == 0 || ip[j] > maxI) maxI = ip[j];
        }
        bpP[i] = maxI;
        pxP[i] = (int32_t)s.x;
        pyP[i] = (int32_t)s.y;
        pzP[i] = (int32_t)s.z;
        cursor += (int64_t)m;
    }

    TTIOWrittenRun *run = [[TTIOWrittenRun alloc]
        initWithSpectrumClassName:@"TTIOMassSpectrum"
                  acquisitionMode:0
                      channelData:@{@"mz": mz, @"intensity": it}
                          offsets:off
                          lengths:len
                   retentionTimes:rt
                         msLevels:ml
                       polarities:pol
                     precursorMzs:pmz
                 precursorCharges:pc
              basePeakIntensities:bp];
    run.pixelX = px;
    run.pixelY = py;
    run.pixelZ = pz;
    run.provenanceRecords = @[[self provenanceRecordForImport:import]];
    return run;
}

+ (nullable TTIOImportedDataset *)importedDatasetFromImport:(TTIOImzMLImport *)import
                                                      title:(nullable NSString *)title
                                                      error:(NSError **)error
{
    TTIOWrittenRun *run = [self pixelRunFromImport:import error:error];
    if (!run) return nil;
    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    d.title = title.length
        ? title
        : [NSString stringWithFormat:@"imzML import: %@",
           [import.sourceImzML lastPathComponent]];
    d.msRuns[TTIOImzMLPixelRunName] = run;
    [d.provenanceRecords addObjectsFromArray:run.provenanceRecords];
    return d;
}

@end
