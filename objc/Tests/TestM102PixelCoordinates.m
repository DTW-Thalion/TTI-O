// TestM102PixelCoordinates.m — M102: imzML pixel coordinates as
// spectrum_index columns (format-spec §4b, opt_pixel_coordinates).
//
//   - a processed-mode imzML with 15,000 pixels imports, stores its
//     positions as pixel_x/y/z columns (the pre-M102 provenance CSV
//     overflowed the 64 KB @provenance_json attribute past ~9,600
//     pixels on a 100-wide grid), and exports back to imzML intact;
//   - a pre-M102 file whose positions live in the
//     imzml_pixel_coordinates_csv parameter still reads and exports;
//   - spectrum_index all-or-none validation, partial columns on disk,
//     absent columns.
//
// Mirrors:
//   python/tests/test_m102_pixel_coordinates.py
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#import <Foundation/Foundation.h>
#import "Testing.h"
#import "Core/TTIOSignalArray.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Dataset/TTIOWrittenRun.h"
#import "Export/TTIOExporterRegistry.h"
#import "Export/TTIOImzMLWriter.h"
#import "Export/TTIOWriterAdapters.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "Import/TTIOImportedDataset.h"
#import "Import/TTIOImporterRegistry.h"
#import "Import/TTIOImzMLReader.h"
#import "Import/TTIOReaderAdapters.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOInstrumentConfig.h"
#import "Run/TTIOSpectrumIndex.h"
#import "Spectra/TTIOMassSpectrum.h"
#import "ValueClasses/TTIOEncodingSpec.h"
#import "ValueClasses/TTIOEnums.h"
#include <unistd.h>

static const NSInteger kM102W = 150;
static const NSInteger kM102H = 100;
static NSString *const kM102Uuid = @"0123456789abcdef0123456789abcdef";

static NSString *m102Dir(NSString *tag)
{
    NSString *d = [NSTemporaryDirectory() stringByAppendingPathComponent:
        [NSString stringWithFormat:@"ttio-m102-%@-%d", tag, (int)getpid()]];
    [[NSFileManager defaultManager] removeItemAtPath:d error:NULL];
    [[NSFileManager defaultManager] createDirectoryAtPath:d
                              withIntermediateDirectories:YES attributes:nil error:NULL];
    return d;
}

static void m102Rm(NSString *path)
{
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
}

static BOOL m102HasFeature(NSString *path, NSString *flag)
{
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
        if (!f) return NO;
        BOOL has = [TTIOFeatureFlags root:f.rootGroup supportsFeature:flag];
        [f close];
        return has;
    }
}

static BOOL m102AnyRecordHas(NSArray<TTIOProvenanceRecord *> *records, NSString *key)
{
    for (TTIOProvenanceRecord *r in records) if (r.parameters[key] != nil) return YES;
    return NO;
}

/* Source pixel i of the synthetic processed-mode grid: 3-5 peaks with
 * an m/z axis that differs per pixel. */
static TTIOImzMLPixelSpectrum *m102SourcePixel(NSInteger x, NSInteger y)
{
    NSUInteger k = 3 + (NSUInteger)((x + y) % 3);
    double mz[5], in[5];
    for (NSUInteger j = 0; j < k; j++) {
        mz[j] = 100.0 + (double)j * 50.0 + (double)x * 0.01 + (double)y * 0.0001;
        in[j] = (double)((x * 7 + y * 13 + (NSInteger)j * 31) % 1000) + 1.0;
    }
    return [[TTIOImzMLPixelSpectrum alloc]
        initWithX:x y:y z:1
          mzArray:[NSData dataWithBytes:mz length:k * sizeof(double)]
   intensityArray:[NSData dataWithBytes:in length:k * sizeof(double)]
            error:NULL];
}

static TTIOSignalArray *m102F64(NSData *buf)
{
    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
    return [[TTIOSignalArray alloc] initWithBuffer:buf
                                            length:buf.length / sizeof(double)
                                          encoding:enc
                                              axis:nil];
}

static TTIOAcquisitionRun *m102InMemoryRun(NSArray<TTIOImzMLPixelSpectrum *> *pixels)
{
    NSMutableArray *spectra = [NSMutableArray array];
    for (NSUInteger i = 0; i < pixels.count; i++) {
        TTIOMassSpectrum *s = [[TTIOMassSpectrum alloc]
            initWithMzArray:m102F64(pixels[i].mzArray)
             intensityArray:m102F64(pixels[i].intensityArray)
                    msLevel:1
                   polarity:TTIOPolarityUnknown
                 scanWindow:nil
              indexPosition:i
            scanTimeSeconds:0
                precursorMz:0
            precursorCharge:0
                      error:NULL];
        [spectra addObject:s];
    }
    TTIOInstrumentConfig *cfg =
        [[TTIOInstrumentConfig alloc] initWithManufacturer:@"" model:@""
                                              serialNumber:@"" sourceType:@""
                                              analyzerType:@"" detectorType:@""];
    return [[TTIOAcquisitionRun alloc] initWithSpectra:spectra
                                       acquisitionMode:TTIOAcquisitionModeMS1DDA
                                      instrumentConfig:cfg];
}

static BOOL m102SameArrays(TTIOMassSpectrum *s, TTIOImzMLPixelSpectrum *p)
{
    return [[s.mzArray float64Buffer] isEqualToData:p.mzArray]
        && [[s.intensityArray float64Buffer] isEqualToData:p.intensityArray];
}

/* Re-read an exported imzML and compare every pixel with `src`. */
static BOOL m102ExportMatches(TTIOImzMLImport *imp,
                              NSArray<TTIOImzMLPixelSpectrum *> *src)
{
    if (imp.spectra.count != src.count) return NO;
    for (NSUInteger i = 0; i < src.count; i++) {
        TTIOImzMLPixelSpectrum *a = imp.spectra[i], *b = src[i];
        if (a.x != b.x || a.y != b.y || a.z != b.z) return NO;
        if (![a.mzArray isEqualToData:b.mzArray]) return NO;
        if (![a.intensityArray isEqualToData:b.intensityArray]) return NO;
    }
    return YES;
}

#pragma mark - Large processed-mode import

static void testLargeProcessedImport(void)
{
    NSString *dir = m102Dir(@"large");
    NSString *imzml = [dir stringByAppendingPathComponent:@"grid.imzML"];
    NSString *tio = [dir stringByAppendingPathComponent:@"grid.tio"];
    NSString *out = [dir stringByAppendingPathComponent:@"back.imzML"];
    NSString *outIbd = [dir stringByAppendingPathComponent:@"back.ibd"];

    NSMutableArray<TTIOImzMLPixelSpectrum *> *src = [NSMutableArray array];
    for (NSInteger y = 1; y <= kM102H; y++)
        for (NSInteger x = 1; x <= kM102W; x++)
            [src addObject:m102SourcePixel(x, y)];
    NSUInteger n = src.count;
    PASS(n == 15000, "M102 large: 150 x 100 = 15,000 source pixels");

    NSError *err = nil;
    TTIOImzMLWriteResult *w = [TTIOImzMLWriter
        writePixels:src toImzMLPath:imzml ibdPath:nil mode:@"processed"
           gridMaxX:kM102W gridMaxY:kM102H gridMaxZ:1
         pixelSizeX:10.0 pixelSizeY:20.0 scanPattern:@"meandering"
            uuidHex:kM102Uuid error:&err];
    PASS(w != nil, "M102 large: processed-mode imzML fixture written");

    // The adapter now imports processed mode as the pixel run.
    err = nil;
    TTIOImportedDataset *draft = [[TTIOImzMLReaderAdapter new]
        readInputs:@[imzml] options:@{} progress:nil error:&err];
    PASS(draft != nil, "M102 large: adapter accepts processed mode");
    PASS(draft.msImage == nil && draft.writeDelegate == nil,
         "M102 large: processed draft carries no image cube");
    PASS(draft.msRuns[@"imzml_pixels"] != nil,
         "M102 large: draft carries the imzml_pixels run");
    PASS([draft.title isEqualToString:@"imzML import: grid.imzML"],
         "M102 large: draft title names the imzML file");
    PASS(draft.provenanceRecords.count == 1,
         "M102 large: draft carries one dataset-level record");

    err = nil;
    BOOL ok = [TTIOImporterRegistry encodeFormat:@"imzml" inputs:@[imzml]
                                          output:tio options:@{} error:&err];
    PASS(ok, "M102 large: registry imports 15,000 pixels to .tio");
    if (!ok) NSLog(@"M102 large import error: %@", err);

    PASS(m102HasFeature(tio, [TTIOFeatureFlags featurePixelCoordinates]),
         "M102 large: opt_pixel_coordinates set");

    err = nil;
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:tio error:&err];
    PASS(ds != nil, "M102 large: .tio reopens");
    TTIOAcquisitionRun *run = ds.msRuns[@"imzml_pixels"];
    TTIOSpectrumIndex *idx = run.spectrumIndex;
    PASS(idx.hasPixelCoordinates && idx.count == n,
         "M102 large: spectrum index carries pixel columns for every pixel");
    BOOL same = idx.hasPixelCoordinates && idx.count == n;
    for (NSUInteger i = 0; same && i < n; i++) {
        same = [idx pixelXAt:i] == src[i].x && [idx pixelYAt:i] == src[i].y
            && [idx pixelZAt:i] == src[i].z;
    }
    PASS(same, "M102 large: pixel columns equal the source coordinates");

    NSData *helper = [TTIOImzMLReader pixelCoordinatesForRun:run
                                           datasetProvenance:ds.provenanceRecords];
    const int32_t *h = helper.bytes;
    PASS(helper.length == n * 3 * sizeof(int32_t)
         && h[0] == 1 && h[1] == 1 && h[2] == 1
         && h[(n - 1) * 3] == kM102W && h[(n - 1) * 3 + 1] == kM102H,
         "M102 large: helper returns the interleaved columns");

    PASS(!m102AnyRecordHas([run provenanceChain], TTIOImzMLLegacyCoordinatesParameter)
         && !m102AnyRecordHas(ds.provenanceRecords, TTIOImzMLLegacyCoordinatesParameter),
         "M102 large: no record carries imzml_pixel_coordinates_csv");
    TTIOProvenanceRecord *rec = [run provenanceChain].firstObject;
    PASS([rec.parameters[@"imzml_mode"] isEqual:@"processed"]
         && [rec.parameters[@"imzml_grid_max_x"] integerValue] == kM102W
         && [rec.parameters[@"imzml_uuid_hex"] isEqual:kM102Uuid]
         && [rec.software isEqualToString:@"ttio imzml importer v0.9"],
         "M102 large: run-level record keeps the imzml_* scalars");
    PASS(ds.provenanceRecords.count == 1
         && [ds.provenanceRecords[0].parameters[@"imzml_scan_pattern"] isEqual:@"meandering"],
         "M102 large: dataset-level record keeps the imzml_* scalars");

    TTIOMassSpectrum *s0 = [run spectrumAtIndex:0 error:NULL];
    TTIOMassSpectrum *sl = [run spectrumAtIndex:n - 1 error:NULL];
    PASS(m102SameArrays(s0, src[0]) && m102SameArrays(sl, src[n - 1]),
         "M102 large: first/last pixel arrays survive the import");
    PASS([idx msLevelAt:7] == 1 && [idx basePeakIntensityAt:0] > 0,
         "M102 large: ms level 1, base peak set");
    [ds closeFile];

    err = nil;
    ok = [TTIOExporterRegistry exportFormat:@"imzml" tioPath:tio layer:nil
                                     output:out options:@{} error:&err];
    PASS(ok, "M102 large: writer adapter exports the pixel run as imzML");
    if (!ok) NSLog(@"M102 large export error: %@", err);

    err = nil;
    TTIOImzMLImport *back = [TTIOImzMLReader readFromImzMLPath:out ibdPath:outIbd error:&err];
    PASS(back != nil && [back.mode isEqualToString:@"processed"],
         "M102 large: exported imzML re-reads in processed mode");
    PASS([back.uuidHex isEqualToString:kM102Uuid] && back.gridMaxX == kM102W
         && back.gridMaxY == kM102H && back.pixelSizeX == 10.0
         && back.pixelSizeY == 20.0,
         "M102 large: UUID, grid and pixel size come from the parameters");
    PASS(m102ExportMatches(back, src),
         "M102 large: exported coordinates and arrays equal the source");

    m102Rm(dir);
}

#pragma mark - Pre-M102 provenance CSV

static void testLegacyCsv(void)
{
    NSString *dir = m102Dir(@"legacy");
    NSString *tio = [dir stringByAppendingPathComponent:@"legacy.tio"];
    NSString *out = [dir stringByAppendingPathComponent:@"legacy.imzML"];
    NSString *outIbd = [dir stringByAppendingPathComponent:@"legacy.ibd"];

    NSArray<TTIOImzMLPixelSpectrum *> *src = @[
        m102SourcePixel(1, 1), m102SourcePixel(2, 1),
        m102SourcePixel(1, 2), m102SourcePixel(2, 2)];
    NSString *csv = @"1,1,1;2,1,1;1,2,1;2,2,1";
    NSDictionary *params = @{
        TTIOImzMLLegacyCoordinatesParameter: csv,
        @"imzml_mode": @"processed",
        @"imzml_uuid_hex": kM102Uuid,
        @"imzml_grid_max_x": @2, @"imzml_grid_max_y": @2, @"imzml_grid_max_z": @1,
        @"imzml_pixel_size_x": @5.0, @"imzml_pixel_size_y": @5.0,
        @"imzml_scan_pattern": @"flyback",
    };
    TTIOProvenanceRecord *rec = [[TTIOProvenanceRecord alloc]
        initWithInputRefs:@[@"legacy.imzML", @"legacy.ibd"]
                 software:@"ttio imzml importer v0.9"
               parameters:params outputRefs:@[] timestampUnix:1700000000];

    TTIOAcquisitionRun *run = m102InMemoryRun(src);
    [run addProcessingStep:rec];
    TTIOSpectralDataset *ds = [[TTIOSpectralDataset alloc]
        initWithTitle:@"legacy" isaInvestigationId:@""
               msRuns:@{@"imzml_pixels": run} nmrRuns:@{}
      identifications:@[] quantifications:@[] provenanceRecords:@[rec]
          transitions:nil];
    NSError *err = nil;
    PASS([ds writeToFilePath:tio error:&err], "M102 legacy: pre-M102-shaped .tio written");
    PASS(!m102HasFeature(tio, [TTIOFeatureFlags featurePixelCoordinates]),
         "M102 legacy: no pixel columns, no opt_pixel_coordinates");

    TTIOSpectralDataset *re = [TTIOSpectralDataset readFromFilePath:tio error:&err];
    TTIOAcquisitionRun *rr = re.msRuns[@"imzml_pixels"];
    PASS(rr != nil && !rr.spectrumIndex.hasPixelCoordinates
         && rr.spectrumIndex.pixelX == nil,
         "M102 legacy: reopened run has no pixel columns");
    NSData *c = [TTIOImzMLReader pixelCoordinatesForRun:rr datasetProvenance:@[]];
    const int32_t *p = c.bytes;
    PASS(c.length == 12 * sizeof(int32_t)
         && p[0] == 1 && p[1] == 1 && p[2] == 1
         && p[3] == 2 && p[4] == 1 && p[6] == 1 && p[7] == 2
         && p[9] == 2 && p[10] == 2 && p[11] == 1,
         "M102 legacy: helper parses the run-level CSV");
    [re closeFile];

    err = nil;
    BOOL ok = [TTIOExporterRegistry exportFormat:@"imzml" tioPath:tio layer:nil
                                          output:out options:@{} error:&err];
    PASS(ok, "M102 legacy: writer adapter exports the CSV-positioned run");
    if (!ok) NSLog(@"M102 legacy export error: %@", err);
    TTIOImzMLImport *back = [TTIOImzMLReader readFromImzMLPath:out ibdPath:outIbd error:&err];
    PASS(back != nil && [back.mode isEqualToString:@"processed"]
         && [back.uuidHex isEqualToString:kM102Uuid] && back.gridMaxX == 2,
         "M102 legacy: mode, UUID and grid come from the parameters");
    PASS(m102ExportMatches(back, src),
         "M102 legacy: exported coordinates and arrays equal the source");

    // Dataset-level fallback, and a CSV with the wrong triple count.
    TTIOAcquisitionRun *bare = m102InMemoryRun(src);
    NSData *dc = [TTIOImzMLReader pixelCoordinatesForRun:bare datasetProvenance:@[rec]];
    PASS(dc.length == 12 * sizeof(int32_t),
         "M102 legacy: helper falls back to dataset-level records");
    PASS([TTIOImzMLReader pixelCoordinatesForRun:bare datasetProvenance:@[]] == nil,
         "M102 legacy: no columns and no CSV -> nil");
    TTIOProvenanceRecord *shortRec = [[TTIOProvenanceRecord alloc]
        initWithInputRefs:@[] software:@"x"
               parameters:@{TTIOImzMLLegacyCoordinatesParameter: @"1,1,1;2,1,1;1,2,1"}
               outputRefs:@[] timestampUnix:0];
    TTIOProvenanceRecord *badRec = [[TTIOProvenanceRecord alloc]
        initWithInputRefs:@[] software:@"x"
               parameters:@{TTIOImzMLLegacyCoordinatesParameter: @"1,1,1;2,1;1,2,1;2,2,1"}
               outputRefs:@[] timestampUnix:0];
    NSArray *badRecs = @[shortRec, badRec];
    PASS([TTIOImzMLReader pixelCoordinatesForRun:bare datasetProvenance:badRecs] == nil,
         "M102 legacy: CSV without one triple per spectrum is ignored");

    // Continuous-mode inference when imzml_mode is absent.
    NSArray<TTIOImzMLPixelSpectrum *> *shared = @[
        [[TTIOImzMLPixelSpectrum alloc] initWithX:1 y:1 z:1
            mzArray:src[0].mzArray intensityArray:src[0].intensityArray error:NULL],
        [[TTIOImzMLPixelSpectrum alloc] initWithX:2 y:1 z:1
            mzArray:src[0].mzArray intensityArray:src[0].intensityArray error:NULL]];
    TTIOAcquisitionRun *cont = m102InMemoryRun(shared);
    NSMutableData *xs = [NSMutableData data], *ys = [NSMutableData data],
                  *zs = [NSMutableData data];
    for (TTIOImzMLPixelSpectrum *s in shared) {
        int32_t v[3] = {(int32_t)s.x, (int32_t)s.y, (int32_t)s.z};
        [xs appendBytes:&v[0] length:4];
        [ys appendBytes:&v[1] length:4];
        [zs appendBytes:&v[2] length:4];
    }
    PASS([cont setPixelX:xs pixelY:ys pixelZ:zs error:&err],
         "M102 legacy: in-memory run accepts pixel columns");
    NSString *cOut = [dir stringByAppendingPathComponent:@"cont.imzML"];
    TTIOImzMLWriteResult *cw = [TTIOImzMLWriter writeRun:cont datasetProvenance:nil
                                             toImzMLPath:cOut ibdPath:nil error:&err];
    PASS(cw != nil && [cw.mode isEqualToString:@"continuous"],
         "M102 legacy: shared m/z axis without imzml_mode -> continuous");

    m102Rm(dir);
}

#pragma mark - Spectrum index columns

static void testSpectrumIndexColumns(void)
{
    NSString *dir = m102Dir(@"index");
    NSString *tio = [dir stringByAppendingPathComponent:@"index.tio"];

    NSArray<TTIOImzMLPixelSpectrum *> *src = @[
        m102SourcePixel(3, 4), m102SourcePixel(5, 6), m102SourcePixel(7, 8)];
    TTIOAcquisitionRun *run = m102InMemoryRun(src);
    TTIOSpectrumIndex *idx = run.spectrumIndex;
    PASS(!idx.hasPixelCoordinates && idx.pixelX == nil && [idx pixelXAt:0] == 0,
         "M102 index: columns absent by default");

    int32_t xv[3] = {3, 5, 7}, yv[3] = {4, 6, 8}, zv[3] = {1, 1, 1};
    NSData *xs = [NSData dataWithBytes:xv length:sizeof(xv)];
    NSData *ys = [NSData dataWithBytes:yv length:sizeof(yv)];
    NSData *zs = [NSData dataWithBytes:zv length:sizeof(zv)];
    NSError *err = nil;
    PASS([idx indexWithPixelX:xs pixelY:nil pixelZ:nil error:&err] == nil && err != nil,
         "M102 index: partial columns rejected");
    err = nil;
    PASS([idx indexWithPixelX:xs pixelY:ys pixelZ:[zs subdataWithRange:NSMakeRange(0, 8)]
                        error:&err] == nil && err != nil,
         "M102 index: wrong column length rejected");
    TTIOSpectrumIndex *with = [idx indexWithPixelX:xs pixelY:ys pixelZ:zs error:&err];
    PASS(with.hasPixelCoordinates && [with pixelXAt:1] == 5 && [with pixelYAt:2] == 8
         && [with pixelZAt:0] == 1 && with.count == 3
         && [with lengthAt:1] == [idx lengthAt:1] && !idx.hasPixelCoordinates,
         "M102 index: copy carries the columns, original unchanged");
    PASS(![[with indexWithPixelX:nil pixelY:nil pixelZ:nil error:NULL] hasPixelCoordinates],
         "M102 index: three nils drop the columns");

    // Object-mode writer: columns + flag.
    PASS([run setPixelX:xs pixelY:ys pixelZ:zs error:&err],
         "M102 index: run accepts pixel columns");
    TTIOSpectralDataset *ds = [[TTIOSpectralDataset alloc]
        initWithTitle:@"idx" isaInvestigationId:@"" msRuns:@{@"r": run} nmrRuns:@{}
      identifications:@[] quantifications:@[] provenanceRecords:@[] transitions:nil];
    PASS([ds writeToFilePath:tio error:&err], "M102 index: dataset with columns written");
    PASS(m102HasFeature(tio, [TTIOFeatureFlags featurePixelCoordinates]),
         "M102 index: object-mode writer sets opt_pixel_coordinates");
    TTIOSpectralDataset *re = [TTIOSpectralDataset readFromFilePath:tio error:&err];
    TTIOSpectrumIndex *ri = re.msRuns[@"r"].spectrumIndex;
    PASS([ri.pixelX isEqualToData:xs] && [ri.pixelY isEqualToData:ys]
         && [ri.pixelZ isEqualToData:zs],
         "M102 index: columns round-trip through the object-mode writer");
    [re closeFile];

    // writeMinimal rejects a run with only some columns.
    PASS([TTIOImzMLReader pixelRunFromImport:[[TTIOImzMLImport alloc] init]
                                       error:NULL] == nil,
         "M102 index: empty import builds no pixel run");
    double one = 1.0;
    NSData *d1 = [NSData dataWithBytes:&one length:sizeof(one)];
    int64_t off0 = 0; uint32_t len1 = 1; int32_t i1 = 1;
    TTIOWrittenRun *wr = [[TTIOWrittenRun alloc]
        initWithSpectrumClassName:@"TTIOMassSpectrum" acquisitionMode:0
                      channelData:@{@"mz": d1, @"intensity": d1}
                          offsets:[NSData dataWithBytes:&off0 length:8]
                          lengths:[NSData dataWithBytes:&len1 length:4]
                   retentionTimes:d1
                         msLevels:[NSData dataWithBytes:&i1 length:4]
                       polarities:[NSData dataWithBytes:&i1 length:4]
                     precursorMzs:d1
                 precursorCharges:[NSData dataWithBytes:&i1 length:4]
              basePeakIntensities:d1];
    wr.pixelX = [NSData dataWithBytes:&i1 length:4];
    NSString *partial = [dir stringByAppendingPathComponent:@"partial.tio"];
    err = nil;
    PASS(![TTIOSpectralDataset writeMinimalToPath:partial title:@"p"
                               isaInvestigationId:@"" msRuns:@{@"r": wr}
                                  identifications:nil quantifications:nil
                                provenanceRecords:nil error:&err] && err != nil,
         "M102 index: writeMinimal rejects a run with only pixelX");
    wr.pixelY = wr.pixelX; wr.pixelZ = wr.pixelX;
    err = nil;
    PASS([TTIOSpectralDataset writeMinimalToPath:partial title:@"p"
                              isaInvestigationId:@"" msRuns:@{@"r": wr}
                                 identifications:nil quantifications:nil
                               provenanceRecords:nil error:&err]
         && m102HasFeature(partial, [TTIOFeatureFlags featurePixelCoordinates]),
         "M102 index: writeMinimal writes complete columns and the flag");
    wr.pixelX = nil; wr.pixelY = nil; wr.pixelZ = nil;
    PASS([TTIOSpectralDataset writeMinimalToPath:partial title:@"p"
                              isaInvestigationId:@"" msRuns:@{@"r": wr}
                                 identifications:nil quantifications:nil
                               provenanceRecords:nil error:&err]
         && !m102HasFeature(partial, [TTIOFeatureFlags featurePixelCoordinates]),
         "M102 index: writeMinimal without columns adds no flag");

    // Partial columns on disk -> malformed.
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openAtPath:tio error:&err];
        TTIOHDF5Group *g = [[[[f.rootGroup openGroupNamed:@"study" error:NULL]
                                openGroupNamed:@"ms_runs" error:NULL]
                               openGroupNamed:@"r" error:NULL]
                              openGroupNamed:@"spectrum_index" error:NULL];
        PASS([g deleteChildNamed:@"pixel_z" error:&err],
             "M102 index: pixel_z removed to simulate a partial file");
        g = nil;
        [f close];
    }
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:tio error:&err];
        TTIOHDF5Group *rg = [[[f.rootGroup openGroupNamed:@"study" error:NULL]
                                 openGroupNamed:@"ms_runs" error:NULL]
                                openGroupNamed:@"r" error:NULL];
        err = nil;
        TTIOSpectrumIndex *bad = [TTIOSpectrumIndex readFromGroup:(id)rg error:&err];
        PASS(bad == nil && err != nil,
             "M102 index: partial pixel columns on disk -> error");
        rg = nil;
        [f close];
    }
    err = nil;
    TTIOSpectralDataset *badDs = [TTIOSpectralDataset readFromFilePath:tio error:&err];
    PASS(badDs == nil && err != nil,
         "M102 index: dataset open fails on a partial pixel column set");

    m102Rm(dir);
}

void testM102PixelCoordinates(void)
{
    @autoreleasepool { testSpectrumIndexColumns(); }
    @autoreleasepool { testLegacyCsv(); }
    @autoreleasepool { testLargeProcessedImport(); }
}
