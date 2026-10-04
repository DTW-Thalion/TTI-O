/*
 * TTIOReaderAdapters.m
 * TTI-O Objective-C Implementation
 *
 * Licensed under the Apache License, Version 2.0.
 * See LICENSE-IMPORT-EXPORT in the repository root.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * OT5: the 11 per-format TTIOReader adapters. Each adapter normalizes one
 * importer into a TTIOImportedDataset draft, mirroring the merged Java
 * importers.readers.*Adapter classes + the Python readers behavior. The write
 * happens at -[TTIOImportedDataset writeToPath:error:], never in -readInputs:.
 */
#import "Import/TTIOReaderAdapters.h"
#import "Import/TTIOImportedDataset.h"

#import "Import/TTIOMzMLReader.h"
#import "Import/TTIONmrMLReader.h"
#import "Import/TTIOThermoRawReader.h"
#import "Import/TTIOWatersMassLynxReader.h"
#import "Import/TTIOMzTabReader.h"
#import "Import/TTIOImzMLReader.h"
#import "Import/TTIOJcampDxReader.h"
#import "Import/TTIOBrukerTDFReader.h"
#import "Import/TTIOBamReader.h"
#import "Import/TTIOSamReader.h"
#import "Import/TTIOCramReader.h"
#import "Import/TTIOGenomicStreamSource.h"
#import "Import/TTIOSpectralStreamSource.h"

#import "Dataset/TTIOSpectralDataset.h"
#import "Image/TTIOMSImage.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOInstrumentConfig.h"
#import "Spectra/TTIOSpectrum.h"
#import "Spectra/TTIORamanSpectrum.h"
#import "Spectra/TTIOIRSpectrum.h"
#import "Spectra/TTIOUVVisSpectrum.h"
#import "Core/TTIOSignalArray.h"
#import "ValueClasses/TTIOEnums.h"

static NSString *const kTTIOReaderAdaptersErrorDomain =
    @"global.thalion.ttio.ReaderAdapters";

static NSError *_ofrError(NSInteger code, NSString *msg)
{
    return [NSError errorWithDomain:kTTIOReaderAdaptersErrorDomain
                               code:code
                           userInfo:@{NSLocalizedDescriptionKey: msg}];
}

/* opts string helper: returns the string for `key`, else `dflt`. */
static NSString *_ofrOptString(NSDictionary<NSString *, id> *opts,
                               NSString *key, NSString *dflt)
{
    id v = opts[key];
    return [v isKindOfClass:[NSString class]] && [(NSString *)v length]
        ? (NSString *)v : dflt;
}

// The streaming knobs an importer accepts in its opts (the --extra k=v
// pairs of TtioEncode): reference, embed_reference, batch_reads,
// batch_spectra, block_reads, block_bytes, legacy_whole_channel.
static BOOL _ofrOptFlag(NSDictionary<NSString *, id> *opts, NSString *key)
{
    id v = opts[key];
    if ([v isKindOfClass:[NSNumber class]]) return [v boolValue];
    if ([v isKindOfClass:[NSString class]]) {
        NSString *l = [(NSString *)v lowercaseString];
        return [l isEqualToString:@"1"] || [l isEqualToString:@"true"] || [l isEqualToString:@"yes"];
    }
    return NO;
}

static NSNumber *_ofrOptNumber(NSDictionary<NSString *, id> *opts, NSString *key)
{
    id v = opts[key];
    if ([v isKindOfClass:[NSNumber class]]) return v;
    if ([v isKindOfClass:[NSString class]] && [(NSString *)v length]) {
        return @([(NSString *)v longLongValue]);
    }
    return nil;
}

static NSUInteger _ofrBatchReads(NSDictionary<NSString *, id> *opts)
{
    NSNumber *n = _ofrOptNumber(opts, @"batch_reads");
    return n && [n unsignedIntegerValue] > 0 ? [n unsignedIntegerValue] : TTIOBamReaderDefaultBatchReads;
}

static NSUInteger _ofrBatchSpectra(NSDictionary<NSString *, id> *opts)
{
    NSNumber *n = _ofrOptNumber(opts, @"batch_spectra");
    return n && [n unsignedIntegerValue] > 0 ? [n unsignedIntegerValue] : [TTIOMzMLReader defaultBatchSpectra];
}

#pragma mark - mzML (write-through delegate)

@implementation TTIOMzMLReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"mzML import: no input path");
        return nil;
    }
    NSString *stem = [[inputs[0] lastPathComponent] stringByDeletingPathExtension];
    NSString *name = _ofrOptString(options, @"name", stem ?: @"run");
    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    d.title = stem ?: @"";
    d.spectralStreams[name] = [TTIOMzMLReader streamFromPath:inputs[0] runName:name
                                                batchSpectra:_ofrBatchSpectra(options)
                                                    progress:progress];
    return d;
}
@end

#pragma mark - nmrML (write-through delegate)

@implementation TTIONmrMLReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    (void)options;
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"nmrML import: no input path");
        return nil;
    }
    TTIOSpectralDataset *parsed = progress
        ? [TTIONmrMLReader readFromFilePath:inputs[0] progress:progress error:error]
        : [TTIONmrMLReader readFromFilePath:inputs[0] error:error];
    if (!parsed) return nil;
    return [TTIOImportedDataset datasetWithWriteDelegate:
        ^BOOL(NSString *out, NSError *_Nullable *_Nullable e) {
            return [parsed writeToFilePath:out error:e];
        }];
}
@end

#pragma mark - Thermo .raw (write-through delegate; SDK stub surfaces via reader)

@implementation TTIOThermoRawReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    (void)options; (void)progress;  // ThermoRawReader has no progress overload.
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"Thermo .raw import: no input path");
        return nil;
    }
    // The current ObjC ThermoRawReader is a stub: returns nil + an
    // "SDK dependency missing" NSError. The adapter surfaces that verbatim.
    TTIOSpectralDataset *parsed =
        [TTIOThermoRawReader readFromFilePath:inputs[0] error:error];
    if (!parsed) return nil;
    return [TTIOImportedDataset datasetWithWriteDelegate:
        ^BOOL(NSString *out, NSError *_Nullable *_Nullable e) {
            return [parsed writeToFilePath:out error:e];
        }];
}
@end

#pragma mark - Waters MassLynx .raw (write-through delegate)

@implementation TTIOWatersMassLynxReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    (void)options; (void)progress;  // WatersMassLynxReader has no progress overload.
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"Waters import: no input directory");
        return nil;
    }
    TTIOSpectralDataset *parsed =
        [TTIOWatersMassLynxReader readFromDirectoryPath:inputs[0] error:error];
    if (!parsed) return nil;
    return [TTIOImportedDataset datasetWithWriteDelegate:
        ^BOOL(NSString *out, NSError *_Nullable *_Nullable e) {
            return [parsed writeToFilePath:out error:e];
        }];
}
@end

#pragma mark - mzTab (idents + quants into the draft; non-delegate path)

@implementation TTIOMzTabReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    (void)options;
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"mzTab import: no input path");
        return nil;
    }
    TTIOMzTabImport *im = progress
        ? [TTIOMzTabReader readFromFilePath:inputs[0] progress:progress error:error]
        : [TTIOMzTabReader readFromFilePath:inputs[0] error:error];
    if (!im) return nil;

    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    if (im.title.length) d.title = im.title;
    [d.identifications addObjectsFromArray:im.identifications];
    [d.quantifications addObjectsFromArray:im.quantifications];
    return d;
}
@end

#pragma mark - imzML (image cube + write-through delegate; NEW capability)

@implementation TTIOImzMLReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"imzML import: no input path");
        return nil;
    }
    NSString *imzml = inputs[0];

    // .ibd location: opts["ibd"] if present, else inputs[1], else sibling (nil).
    NSString *ibd = nil;
    id ibdOpt = options[@"ibd"];
    if ([ibdOpt isKindOfClass:[NSString class]] && [(NSString *)ibdOpt length]) {
        ibd = (NSString *)ibdOpt;
    } else if (inputs.count > 1 && inputs[1] != nil &&
               [inputs[1] isKindOfClass:[NSString class]]) {
        ibd = inputs[1];
    }

    TTIOImzMLImport *imp = progress
        ? [TTIOImzMLReader readFromImzMLPath:imzml ibdPath:ibd progress:progress error:error]
        : [TTIOImzMLReader readFromImzMLPath:imzml ibdPath:ibd error:error];
    if (!imp) return nil;

    if (imp.spectra.count == 0) {
        if (error) *error = _ofrError(2,
            [NSString stringWithFormat:@"imzML import: no pixels parsed from %@",
             imzml]);
        return nil;
    }
    if (![imp.mode isEqualToString:@"continuous"]) {
        // M102: processed mode (per-pixel m/z axes) imports as the
        // imzml_pixels MS run, positions in the spectrum_index pixel
        // columns — the shape Python's to_imported_dataset produces.
        return [TTIOImzMLReader importedDatasetFromImport:imp
                                                    title:nil
                                                    error:error];
    }

    NSInteger width  = imp.gridMaxX;
    NSInteger height = imp.gridMaxY;
    TTIOImzMLPixelSpectrum *first = imp.spectra[0];
    NSUInteger sp = first.mzCount;

    // mzAxis from spectrum[0] (float64).
    NSData *mzAxis = first.mzArray;

    // Flat row-major cube[(row * width + col) * sp + s], 1-indexed pixels.
    NSUInteger cubeCount = (NSUInteger)width * (NSUInteger)height * sp;
    NSMutableData *cubeData =
        [NSMutableData dataWithLength:cubeCount * sizeof(double)];
    double *cube = (double *)cubeData.mutableBytes;
    for (TTIOImzMLPixelSpectrum *pix in imp.spectra) {
        NSInteger col = pix.x - 1;  // imzML is 1-indexed
        NSInteger row = pix.y - 1;
        if (row < 0 || row >= height || col < 0 || col >= width) continue;
        const double *pi = (const double *)pix.intensityArray.bytes;
        NSUInteger piCount = pix.intensityArray.length / sizeof(double);
        NSUInteger base = ((NSUInteger)(row * width + col)) * sp;
        NSUInteger n = MIN(piCount, sp);
        for (NSUInteger i = 0; i < n; i++) cube[base + i] = pi[i];
    }

    TTIOMSImage *image =
        [[TTIOMSImage alloc] initWithTitle:@""
                        isaInvestigationId:@""
                           identifications:@[]
                           quantifications:@[]
                         provenanceRecords:@[]
                                     width:(NSUInteger)width
                                    height:(NSUInteger)height
                            spectralPoints:sp
                                  tileSize:0
                                pixelSizeX:imp.pixelSizeX
                                pixelSizeY:imp.pixelSizeY
                               scanPattern:imp.scanPattern
                                      cube:cubeData
                                    mzAxis:mzAxis];

    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    d.msImage = image;
    d.writeDelegate = ^BOOL(NSString *out, NSError *_Nullable *_Nullable e) {
        return [image writeToFilePath:out error:e];
    };
    return d;
}
@end

#pragma mark - JCAMP-DX (single-spectrum run; write-through delegate)

@implementation TTIOJcampDxReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"JCAMP-DX import: no input path");
        return nil;
    }
    NSString *src = inputs[0];
    TTIOSpectrum *spectrum = progress
        ? [TTIOJcampDxReader readSpectrumFromPath:src progress:progress error:error]
        : [TTIOJcampDxReader readSpectrumFromPath:src error:error];
    if (!spectrum) return nil;

    // Mode from the concrete spectrum subclass (Raman / IR / UV-Vis,
    // defaulting to Raman) — mirrors JcampDxReaderAdapter (Java).
    TTIOAcquisitionMode mode;
    if ([spectrum isKindOfClass:[TTIORamanSpectrum class]]) {
        mode = TTIOAcquisitionModeRaman;
    } else if ([spectrum isKindOfClass:[TTIOIRSpectrum class]]) {
        mode = TTIOAcquisitionModeIR;
    } else if ([spectrum isKindOfClass:[TTIOUVVisSpectrum class]]) {
        mode = TTIOAcquisitionModeUVVis;
    } else {
        mode = TTIOAcquisitionModeRaman;
    }

    NSString *runName = _ofrOptString(options, @"name", @"spectrum_0001");
    TTIOInstrumentConfig *cfg =
        [[TTIOInstrumentConfig alloc] initWithManufacturer:@"" model:@""
                                              serialNumber:@"" sourceType:@""
                                              analyzerType:@"" detectorType:@""];
    TTIOAcquisitionRun *run =
        [[TTIOAcquisitionRun alloc] initWithSpectra:@[spectrum]
                                    acquisitionMode:mode
                                   instrumentConfig:cfg];

    NSString *stem = [[src lastPathComponent] stringByDeletingPathExtension];
    TTIOSpectralDataset *dataset =
        [[TTIOSpectralDataset alloc] initWithTitle:(stem ?: @"")
                                isaInvestigationId:@""
                                            msRuns:@{runName: run}
                                           nmrRuns:@{}
                                   identifications:@[]
                                   quantifications:@[]
                                 provenanceRecords:@[]
                                       transitions:nil];

    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    d.title = (stem ?: @"");
    d.writeDelegate = ^BOOL(NSString *out, NSError *_Nullable *_Nullable e) {
        return [dataset writeToFilePath:out error:e];
    };
    return d;
}
@end

#pragma mark - Bruker timsTOF .d (OT4 write-through draft)

@implementation TTIOBrukerReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    (void)options; (void)progress;
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"Bruker import: no input directory");
        return nil;
    }
    // OT4: returns a draft whose writeDelegate runs the Python subprocess at
    // -writeToPath: time; up-front SQLite metadata validation fails fast here.
    return [TTIOBrukerTDFReader readDatasetFromPath:inputs[0] error:error];
}
@end

#pragma mark - Genomic adapters (BAM / SAM / CRAM) — run into draft.genomicRuns

/* Shared genomic-run build: name (default "genomic_0001") / region / sample. */
static TTIOImportedDataset *_ofrBuildGenomicDraft(TTIOBamReader *reader,
                                                  NSDictionary<NSString *, id> *opts,
                                                  TTIOProgressBlock progress,
                                                  NSError *_Nullable *_Nullable error)
{
    (void)error;
    NSString *name   = _ofrOptString(opts, @"name", @"genomic_0001");
    NSString *region = _ofrOptString(opts, @"region", nil);
    NSString *sample = _ofrOptString(opts, @"sample", nil);

    TTIOGenomicStreamSource *src =
        [[reader streamWithName:name region:region sampleName:sample
                 referenceFasta:_ofrOptString(opts, @"reference", nil)
                 embedReference:_ofrOptFlag(opts, @"embed_reference")
                     batchReads:_ofrBatchReads(opts)
                       progress:progress]
            sourceWithBlockReads:_ofrOptNumber(opts, @"block_reads")
                      blockBytes:_ofrOptNumber(opts, @"block_bytes")
                          legacy:_ofrOptFlag(opts, @"legacy_whole_channel")
                      // M103; the cross-language CLIs pass "1"/"0".
                      groupReads:_ofrOptFlag(opts, @"group_reads")];
    TTIOImportedDataset *d = [[TTIOImportedDataset alloc] init];
    d.genomicStreams[name] = src;
    return d;
}

@implementation TTIOBamReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"BAM import: no input path");
        return nil;
    }
    TTIOBamReader *reader = [[TTIOBamReader alloc] initWithPath:inputs[0]];
    return _ofrBuildGenomicDraft(reader, options, progress, error);
}
@end

@implementation TTIOSamReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"SAM import: no input path");
        return nil;
    }
    TTIOSamReader *reader = [[TTIOSamReader alloc] initWithPath:inputs[0]];
    return _ofrBuildGenomicDraft(reader, options, progress, error);
}
@end

@implementation TTIOCramReaderAdapter
- (nullable TTIOImportedDataset *)readInputs:(NSArray<NSString *> *)inputs
                                     options:(NSDictionary<NSString *, id> *)options
                                    progress:(nullable TTIOProgressBlock)progress
                                       error:(NSError *_Nullable *_Nullable)error
{
    if (inputs.count == 0) {
        if (error) *error = _ofrError(1, @"CRAM import: no input path");
        return nil;
    }
    NSString *reference = _ofrOptString(options, @"reference", nil);
    if (!reference) {
        if (error) *error = _ofrError(4,
            @"CRAM import requires a reference FASTA (opts \"reference\")");
        return nil;
    }
    TTIOCramReader *reader =
        [[TTIOCramReader alloc] initWithPath:inputs[0] referenceFasta:reference];
    return _ofrBuildGenomicDraft(reader, options, progress, error);
}
@end
