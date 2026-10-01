/*
 * TTIOSpectrumIndex.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOSpectrumIndex
 * Inherits From: NSObject
 * Declared In:   Run/TTIOSpectrumIndex.h
 *
 * Per-spectrum offsets, lengths, and queryable metadata for one
 * acquisition run. Range queries (RT, ms_level, polarity) operate
 * on the in-memory parallel arrays without touching the signal
 * channels — the compressed-domain query property of the
 * access-unit storage model.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOSpectrumIndex.h"
#import "ValueClasses/TTIOValueRange.h"
#import "ValueClasses/TTIOIsolationWindow.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOHDF5Types.h"
#import "Genomics/TTIOGenomicIndex.h"  // TTIOOffsetsFromLengths (v1.10 #10)

@implementation TTIOSpectrumIndex
{
    NSData *_offsets;          // uint64_t[count]
    NSData *_lengths;          // uint32_t[count]
    NSData *_retentionTimes;   // double[count]
    NSData *_msLevels;         // uint32_t[count] — stored as int32 for HDF5
    NSData *_polarities;       // int32_t[count]
    NSData *_precursorMzs;     // double[count]
    NSData *_precursorCharges; // int32_t[count]
    NSData *_basePeakIntensities; // double[count]
    // nil (legacy) or all-four populated.
    NSData *_activationMethods;     // int32_t[count]
    NSData *_isolationTargetMzs;    // double[count]
    NSData *_isolationLowerOffsets; // double[count]
    NSData *_isolationUpperOffsets; // double[count]
    // nil for legacy files; otherwise int32_t[count] with 0 = profile,
    // 1 = centroided. Mirrors mzML CV terms MS:1000127 / MS:1000128.
    NSData *_centroideds;
    // M102: nil, or all three int32_t[count] imaging-grid positions.
    NSData *_pixelX;
    NSData *_pixelY;
    NSData *_pixelZ;
    NSUInteger _count;
}

- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
             basePeakIntensities:(NSData *)basePeakIntensities
{
    return [self initWithOffsets:offsets
                         lengths:lengths
                  retentionTimes:retentionTimes
                        msLevels:msLevels
                      polarities:polarities
                    precursorMzs:precursorMzs
                precursorCharges:precursorCharges
             basePeakIntensities:basePeakIntensities
               activationMethods:nil
              isolationTargetMzs:nil
           isolationLowerOffsets:nil
           isolationUpperOffsets:nil];
}

- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
             basePeakIntensities:(NSData *)basePeakIntensities
               activationMethods:(NSData *)activationMethods
             isolationTargetMzs:(NSData *)isolationTargetMzs
          isolationLowerOffsets:(NSData *)isolationLowerOffsets
          isolationUpperOffsets:(NSData *)isolationUpperOffsets
{
    return [self initWithOffsets:offsets
                         lengths:lengths
                  retentionTimes:retentionTimes
                        msLevels:msLevels
                      polarities:polarities
                    precursorMzs:precursorMzs
                precursorCharges:precursorCharges
             basePeakIntensities:basePeakIntensities
               activationMethods:activationMethods
              isolationTargetMzs:isolationTargetMzs
           isolationLowerOffsets:isolationLowerOffsets
           isolationUpperOffsets:isolationUpperOffsets
                      centroideds:nil];
}

- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
             basePeakIntensities:(NSData *)basePeakIntensities
               activationMethods:(NSData *)activationMethods
             isolationTargetMzs:(NSData *)isolationTargetMzs
          isolationLowerOffsets:(NSData *)isolationLowerOffsets
          isolationUpperOffsets:(NSData *)isolationUpperOffsets
                      centroideds:(NSData *)centroideds
{
    BOOL anyNil = !activationMethods || !isolationTargetMzs
               || !isolationLowerOffsets || !isolationUpperOffsets;
    BOOL allNil = !activationMethods && !isolationTargetMzs
               && !isolationLowerOffsets && !isolationUpperOffsets;
    NSAssert(allNil || !anyNil,
        @"TTIOSpectrumIndex: M74 columns must be all-nil or all-non-nil");
    self = [super init];
    if (self) {
        _offsets             = [offsets copy];
        _lengths             = [lengths copy];
        _retentionTimes      = [retentionTimes copy];
        _msLevels            = [msLevels copy];
        _polarities          = [polarities copy];
        _precursorMzs        = [precursorMzs copy];
        _precursorCharges    = [precursorCharges copy];
        _basePeakIntensities = [basePeakIntensities copy];
        _activationMethods     = [activationMethods copy];
        _isolationTargetMzs    = [isolationTargetMzs copy];
        _isolationLowerOffsets = [isolationLowerOffsets copy];
        _isolationUpperOffsets = [isolationUpperOffsets copy];
        _centroideds         = [centroideds copy];
        _count               = offsets.length / sizeof(uint64_t);
    }
    return self;
}

- (NSUInteger)count { return _count; }

- (uint64_t)offsetAt:(NSUInteger)i { return ((const uint64_t *)_offsets.bytes)[i]; }
- (uint32_t)lengthAt:(NSUInteger)i { return ((const uint32_t *)_lengths.bytes)[i]; }
- (double)retentionTimeAt:(NSUInteger)i { return ((const double *)_retentionTimes.bytes)[i]; }
- (uint8_t)msLevelAt:(NSUInteger)i      { return (uint8_t)((const int32_t *)_msLevels.bytes)[i]; }
- (TTIOPolarity)polarityAt:(NSUInteger)i { return (TTIOPolarity)((const int32_t *)_polarities.bytes)[i]; }
- (double)precursorMzAt:(NSUInteger)i   { return ((const double *)_precursorMzs.bytes)[i]; }
- (uint8_t)precursorChargeAt:(NSUInteger)i { return (uint8_t)((const int32_t *)_precursorCharges.bytes)[i]; }
- (double)basePeakIntensityAt:(NSUInteger)i { return ((const double *)_basePeakIntensities.bytes)[i]; }

- (BOOL)hasActivationDetail { return _activationMethods != nil; }

- (TTIOActivationMethod)activationMethodAt:(NSUInteger)i
{
    if (!_activationMethods) return TTIOActivationMethodNone;
    return (TTIOActivationMethod)((const int32_t *)_activationMethods.bytes)[i];
}

- (TTIOIsolationWindow *)isolationWindowAt:(NSUInteger)i
{
    if (!_isolationTargetMzs) return nil;
    double t = ((const double *)_isolationTargetMzs.bytes)[i];
    double lo = ((const double *)_isolationLowerOffsets.bytes)[i];
    double hi = ((const double *)_isolationUpperOffsets.bytes)[i];
    if (t == 0.0 && lo == 0.0 && hi == 0.0) return nil;
    return [TTIOIsolationWindow windowWithTargetMz:t
                                        lowerOffset:lo
                                        upperOffset:hi];
}

- (BOOL)hasCentroided { return _centroideds != nil; }

- (BOOL)centroidedAt:(NSUInteger)i
{
    if (!_centroideds) return NO;
    return ((const int32_t *)_centroideds.bytes)[i] != 0;
}

#pragma mark - Pixel coordinates (M102)

- (NSData *)pixelX { return _pixelX; }
- (NSData *)pixelY { return _pixelY; }
- (NSData *)pixelZ { return _pixelZ; }
- (BOOL)hasPixelCoordinates { return _pixelX != nil; }

- (int32_t)pixelXAt:(NSUInteger)i
{
    return _pixelX ? ((const int32_t *)_pixelX.bytes)[i] : 0;
}

- (int32_t)pixelYAt:(NSUInteger)i
{
    return _pixelY ? ((const int32_t *)_pixelY.bytes)[i] : 0;
}

- (int32_t)pixelZAt:(NSUInteger)i
{
    return _pixelZ ? ((const int32_t *)_pixelZ.bytes)[i] : 0;
}

/* All three columns or none, each int32_t[count]. */
static BOOL validatePixelColumns(NSData *x, NSData *y, NSData *z,
                                 NSUInteger count, NSError **error)
{
    if (!x && !y && !z) return YES;
    if (!x || !y || !z) {
        if (error) *error = TTIOMakeError(TTIOErrorUnsupportedLayout,
            @"spectrum_index: pixel_x/pixel_y/pixel_z must be present "
            @"together or not at all");
        return NO;
    }
    NSUInteger want = count * sizeof(int32_t);
    if (x.length != want || y.length != want || z.length != want) {
        if (error) *error = TTIOMakeError(TTIOErrorUnsupportedLayout,
            @"spectrum_index: pixel columns must hold %lu int32 values "
            @"(got %lu/%lu/%lu bytes)", (unsigned long)count,
            (unsigned long)x.length, (unsigned long)y.length,
            (unsigned long)z.length);
        return NO;
    }
    return YES;
}

- (instancetype)indexWithPixelX:(NSData *)pixelX
                         pixelY:(NSData *)pixelY
                         pixelZ:(NSData *)pixelZ
                          error:(NSError **)error
{
    if (!validatePixelColumns(pixelX, pixelY, pixelZ, _count, error)) return nil;
    TTIOSpectrumIndex *out =
        [[[self class] alloc] initWithOffsets:_offsets
                                      lengths:_lengths
                               retentionTimes:_retentionTimes
                                     msLevels:_msLevels
                                   polarities:_polarities
                                 precursorMzs:_precursorMzs
                             precursorCharges:_precursorCharges
                          basePeakIntensities:_basePeakIntensities
                            activationMethods:_activationMethods
                           isolationTargetMzs:_isolationTargetMzs
                        isolationLowerOffsets:_isolationLowerOffsets
                        isolationUpperOffsets:_isolationUpperOffsets
                                  centroideds:_centroideds];
    out->_pixelX = [pixelX copy];
    out->_pixelY = [pixelY copy];
    out->_pixelZ = [pixelZ copy];
    return out;
}

- (NSIndexSet *)indicesInRetentionTimeRange:(TTIOValueRange *)range
{
    const double *rts = _retentionTimes.bytes;
    NSMutableIndexSet *out = [NSMutableIndexSet indexSet];
    for (NSUInteger i = 0; i < _count; i++) {
        if ([range containsValue:rts[i]]) [out addIndex:i];
    }
    return out;
}

- (NSIndexSet *)indicesForMsLevel:(uint8_t)msLevel
{
    const int32_t *ml = _msLevels.bytes;
    NSMutableIndexSet *out = [NSMutableIndexSet indexSet];
    for (NSUInteger i = 0; i < _count; i++) {
        if (ml[i] == msLevel) [out addIndex:i];
    }
    return out;
}

#pragma mark - Storage round-trip (provider-agnostic)

static BOOL writeArray(id<TTIOStorageGroup> g, NSString *name, TTIOPrecision p,
                       NSData *data, NSError **error)
{
    NSUInteger n = data.length / TTIOPrecisionElementSize(p);
    id<TTIOStorageDataset> ds = [g createDatasetNamed:name
                                            precision:p
                                               length:n
                                            chunkSize:4096
                                          compression:TTIOCompressionZlib
                                     compressionLevel:6
                                                error:error];
    if (!ds) return NO;
    return [ds writeAll:data error:error];
}

static NSData *readArray(id<TTIOStorageGroup> g, NSString *name, NSError **error)
{
    id<TTIOStorageDataset> ds = [g openDatasetNamed:name error:error];
    if (!ds) return nil;
    id val = [ds readAll:error];
    return [val isKindOfClass:[NSData class]] ? val : nil;
}

- (BOOL)writeToGroup:(id<TTIOStorageGroup>)parent error:(NSError **)error
{
    id<TTIOStorageGroup> g = [parent createGroupNamed:@"spectrum_index" error:error];
    if (!g) return NO;
    if (![g setAttributeValue:@((int64_t)_count) forName:@"count" error:error]) return NO;
    // offsets is omitted on disk; readers compute it from
    // cumsum(lengths).
    if (!writeArray(g, @"lengths",          TTIOPrecisionUInt32,  _lengths,          error)) return NO;
    if (!writeArray(g, @"retention_times",  TTIOPrecisionFloat64, _retentionTimes,   error)) return NO;
    if (!writeArray(g, @"ms_levels",        TTIOPrecisionInt32,   _msLevels,         error)) return NO;
    if (!writeArray(g, @"polarities",       TTIOPrecisionInt32,   _polarities,       error)) return NO;
    if (!writeArray(g, @"precursor_mzs",    TTIOPrecisionFloat64, _precursorMzs,     error)) return NO;
    if (!writeArray(g, @"precursor_charges", TTIOPrecisionInt32,  _precursorCharges, error)) return NO;
    if (!writeArray(g, @"base_peak_intensities", TTIOPrecisionFloat64, _basePeakIntensities, error)) return NO;
    // M74 schema-gating: emit the four optional columns only when the
    // index was built with them. The designated initializer enforces
    // all-or-nothing, so probing one column is sufficient.
    if (_activationMethods) {
        if (!writeArray(g, @"activation_methods",      TTIOPrecisionInt32,   _activationMethods,     error)) return NO;
        if (!writeArray(g, @"isolation_target_mzs",    TTIOPrecisionFloat64, _isolationTargetMzs,    error)) return NO;
        if (!writeArray(g, @"isolation_lower_offsets", TTIOPrecisionFloat64, _isolationLowerOffsets, error)) return NO;
        if (!writeArray(g, @"isolation_upper_offsets", TTIOPrecisionFloat64, _isolationUpperOffsets, error)) return NO;
    }
    // Optional centroided column — independent of M74 gating.
    if (_centroideds) {
        if (!writeArray(g, @"centroideds", TTIOPrecisionInt32, _centroideds, error)) return NO;
    }
    // M102 pixel columns (format-spec �4b): emitted only when present,
    // so indexes without them stay byte-identical.
    if (_pixelX) {
        if (!writeArray(g, @"pixel_x", TTIOPrecisionInt32, _pixelX, error)) return NO;
        if (!writeArray(g, @"pixel_y", TTIOPrecisionInt32, _pixelY, error)) return NO;
        if (!writeArray(g, @"pixel_z", TTIOPrecisionInt32, _pixelZ, error)) return NO;
    }
    return YES;
}

+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)parent error:(NSError **)error
{
    id<TTIOStorageGroup> g = [parent openGroupNamed:@"spectrum_index" error:error];
    if (!g) return nil;
    NSData *lengths = readArray(g, @"lengths", error);
    if (!lengths) return nil;
    // offsets is omitted from disk by default; synthesize
    // from cumsum(lengths). Pre-v1.10 files have it on disk.
    NSData *offsets;
    if ([g hasChildNamed:@"offsets"]) {
        offsets = readArray(g, @"offsets", error);
        if (!offsets) return nil;
    } else {
        offsets = TTIOOffsetsFromLengths(lengths);
    }
    NSData *rts     = readArray(g, @"retention_times", error);
    NSData *ml      = readArray(g, @"ms_levels", error);
    NSData *pol     = readArray(g, @"polarities", error);
    NSData *pmz     = readArray(g, @"precursor_mzs", error);
    NSData *pc      = readArray(g, @"precursor_charges", error);
    NSData *bp      = readArray(g, @"base_peak_intensities", error);
    if (!rts || !ml || !pol || !pmz || !pc || !bp) return nil;
    // M74 schema-gating: probe for the four optional columns.
    NSData *am = nil, *itm = nil, *ilo = nil, *iup = nil;
    if ([g hasChildNamed:@"activation_methods"]) {
        am  = readArray(g, @"activation_methods",      error); if (!am)  return nil;
        itm = readArray(g, @"isolation_target_mzs",    error); if (!itm) return nil;
        ilo = readArray(g, @"isolation_lower_offsets", error); if (!ilo) return nil;
        iup = readArray(g, @"isolation_upper_offsets", error); if (!iup) return nil;
    }
    // Optional centroided column — independent of M74 gating.
    NSData *cent = nil;
    if ([g hasChildNamed:@"centroideds"]) {
        cent = readArray(g, @"centroideds", error);
        if (!cent) return nil;
    }
    // M102 pixel columns: all three or none; partial presence is a
    // malformed file.
    BOOL hasPX = [g hasChildNamed:@"pixel_x"];
    BOOL hasPY = [g hasChildNamed:@"pixel_y"];
    BOOL hasPZ = [g hasChildNamed:@"pixel_z"];
    NSData *px = nil, *py = nil, *pz = nil;
    if (hasPX || hasPY || hasPZ) {
        if (!(hasPX && hasPY && hasPZ)) {
            if (error) *error = TTIOMakeError(TTIOErrorUnsupportedLayout,
                @"spectrum_index is malformed: partial pixel_x/pixel_y/"
                @"pixel_z columns present");
            return nil;
        }
        px = readArray(g, @"pixel_x", error); if (!px) return nil;
        py = readArray(g, @"pixel_y", error); if (!py) return nil;
        pz = readArray(g, @"pixel_z", error); if (!pz) return nil;
    }
    TTIOSpectrumIndex *idx = [[self alloc] initWithOffsets:offsets
                                 lengths:lengths
                          retentionTimes:rts
                                msLevels:ml
                              polarities:pol
                            precursorMzs:pmz
                        precursorCharges:pc
                     basePeakIntensities:bp
                        activationMethods:am
                      isolationTargetMzs:itm
                   isolationLowerOffsets:ilo
                   isolationUpperOffsets:iup
                              centroideds:cent];
    if (!px) return idx;
    return [idx indexWithPixelX:px pixelY:py pixelZ:pz error:error];
}

+ (instancetype)readFromStorageGroup:(id)parent error:(NSError **)error
{
    // Legacy alias — readFromGroup: now takes id<TTIOStorageGroup>
    // directly. Kept for source-compatibility with v0.9 callers.
    id<TTIOStorageGroup> par = (id<TTIOStorageGroup>)parent;
    if (![par hasChildNamed:@"spectrum_index"]) return nil;
    return [self readFromGroup:par error:error];
}

@end
