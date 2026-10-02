/*
 * TTIOWriterAdapters.m
 * TTI-O Objective-C Implementation
 *
 * Licensed under the Apache License, Version 2.0.
 * See LICENSE-IMPORT-EXPORT in the repository root.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * OT6: the 8 per-format TTIOWriter adapters. Each adapter serializes one
 * layer of an opened TTIOSpectralDataset to an output path, mirroring the
 * merged Java exporters.writers.*Adapter classes + the Python
 * ttio.exporters.writers behavior. Run selection is shared via
 * TTIORunSelection; the read-side -> write-side genomic conversion for
 * BAM / CRAM goes through +[TTIORunSelection writtenFromGenomicRun:].
 * Error-message strings are kept byte-identical to Python.
 */
#import "Export/TTIOWriterAdapters.h"
#import "Export/TTIORunSelection.h"

#import "Export/TTIOMzMLWriter.h"
#import "Export/TTIOMzTabWriter.h"
#import "Export/TTIONmrMLWriter.h"
#import "Export/TTIOImzMLWriter.h"
#import "Export/TTIOJcampDxWriter.h"
#import "Export/TTIOJcampDxEncoding.h"
#import "Export/TTIOISAExporter.h"
#import "Export/TTIOBamWriter.h"
#import "Export/TTIOCramWriter.h"

#import "Dataset/TTIOSpectralDataset.h"
#import "Image/TTIOImage.h"
#import "Image/TTIOMSImage.h"
#import "Image/TTIOPixelSpectrum.h"
#import "Import/TTIOImzMLReader.h"  /* TTIOImzMLPixelSpectrum */
#import "Run/TTIOAcquisitionRun.h"
#import "Spectra/TTIONMRSpectrum.h"
#import "Spectra/TTIORamanSpectrum.h"
#import "Spectra/TTIOIRSpectrum.h"
#import "Spectra/TTIOUVVisSpectrum.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOWrittenGenomicRun.h"

static NSString *const kTTIOWriterAdaptersErrorDomain =
    @"global.thalion.ttio.WriterAdapters";

/* Build an NSError whose localizedDescription is `msg` — mirroring how
 * Python raises ValueError with that exact text. */
static NSError *_owaError(NSString *msg)
{
    return [NSError errorWithDomain:kTTIOWriterAdaptersErrorDomain
                               code:1
                           userInfo:@{NSLocalizedDescriptionKey: msg}];
}

/* Python's `layer or '(only)'!r` repr: '(only)' when nil/empty, else the
 * single-quoted layer name. */
static NSString *_owaLayerRepr(NSString *layer)
{
    if (layer.length == 0) return @"'(only)'";
    return [NSString stringWithFormat:@"'%@'", layer];
}

/* opts[key] as a non-empty string, else nil. */
static NSString *_owaOptString(NSDictionary<NSString *, id> *opts, NSString *key)
{
    id v = opts[key];
    return ([v isKindOfClass:[NSString class]] && [(NSString *)v length])
        ? (NSString *)v : nil;
}

/* Replace the output file's extension with `ext` (Python
 * Path.with_suffix). `ext` includes the leading dot. */
static NSString *_owaWithSuffix(NSString *output, NSString *ext)
{
    return [[output stringByDeletingPathExtension] stringByAppendingString:ext];
}


#pragma mark - mzML (whole dataset, zlib on)

@implementation TTIOMzMLWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    (void)layer; (void)options;
    // Java GUI: MzMLWriter.write(run, targetPath, /*zlib=*/true, sink).
    return [TTIOMzMLWriter writeDataset:dataset
                                toPath:output
                       zlibCompression:YES
                                 error:error];
}
@end


#pragma mark - mzTab (idents + quants; dialect default "1.0")

@implementation TTIOMzTabWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    (void)layer;
    NSString *version = _owaOptString(options, @"dialect") ?: @"1.0";
    // Java GUI: MzTabWriter.write(targetPath, idents, quants, List.of(),
    //              dialect, title, "", sink).
    TTIOMzTabWriteResult *r =
        [TTIOMzTabWriter writeToPath:output
                    identifications:dataset.identifications
                    quantifications:dataset.quantifications
                            version:version
                              title:dataset.title
                        description:@""
                              error:error];
    return r != nil;
}
@end


#pragma mark - nmrML (NMR run, first spectrum; Python guard)

@implementation TTIONmrMLWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    (void)options;
    TTIOAcquisitionRun *run =
        [TTIORunSelection nmrRunIn:dataset layer:layer error:error];
    if (!run) return NO;

    NSArray *spectra = [run spectra];
    if (spectra.count == 0) {
        if (error) *error = _owaError([NSString stringWithFormat:
            @"run %@ has no spectra", _owaLayerRepr(layer)]);
        return NO;
    }
    id first = spectra.firstObject;
    if (![first isKindOfClass:[TTIONMRSpectrum class]]) {
        if (error) *error = _owaError([NSString stringWithFormat:
            @"run %@ is %@, not an NMR spectrum; pass --layer to select an "
            @"NMR run", _owaLayerRepr(layer),
            NSStringFromClass([first class])]);
        return NO;
    }
    // Java GUI: NmrMLWriter.write(run, targetPath, sink). ObjC nmrML is
    // single-spectrum: export the run's first spectrum.
    return [TTIONmrMLWriter writeSpectrum:(TTIONMRSpectrum *)first
                                     fid:nil
                           sweepWidthPPM:0
                                  toPath:output
                                   error:error];
}
@end


#pragma mark - imzML (ds.msImage; null-guarded; sibling .ibd)

@implementation TTIOImzMLWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    // The ObjC -msImage accessor never returns nil for a valid .tio: when
    // the file has no /study/image_cube it yields a degenerate empty image
    // (width/height/spectralPoints == 0). Treat that as "no image" so the
    // guard matches Java's `ds.image() == null` and Python's `img is None`.
    TTIOMSImage *img = (TTIOMSImage *)[dataset imageForKind:TTIOImageKindMS];
    // Python: ibd = Path(output).with_suffix(".ibd").
    NSString *ibd = _owaWithSuffix(output, @".ibd");
    if (img == nil ||
        img.width == 0 || img.height == 0 || img.spectralPoints == 0) {
        // M102: no cube — an imported pixel run (spectrum_index pixel
        // columns, or the pre-M102 provenance CSV) exports with its own
        // positions. First matching MS run in name order; `layer`
        // restricts the search to that run.
        NSArray<TTIOProvenanceRecord *> *dsProv = dataset.provenanceRecords;
        NSArray *names = [[dataset.msRuns allKeys]
            sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *name in names) {
            if (layer.length && ![name isEqualToString:layer]) continue;
            TTIOAcquisitionRun *run = dataset.msRuns[name];
            if ([TTIOImzMLReader pixelCoordinatesForRun:run
                                      datasetProvenance:dsProv] == nil) continue;
            return [TTIOImzMLWriter writeRun:run
                           datasetProvenance:dsProv
                                 toImzMLPath:output
                                     ibdPath:ibd
                                       error:error] != nil;
        }
        if (error) *error =
            _owaError(@"dataset has no MS image to export as imzML");
        return NO;
    }
    NSString *mode = _owaOptString(options, @"mode") ?: @"continuous";

    // Project the image's pixel spectra (TTIOPixelSpectrum) into the
    // writer's TTIOImzMLPixelSpectrum form. Mirrors Java's
    // ImzMLWriter.write(img.toPixelSpectra(), ...).
    NSArray<TTIOPixelSpectrum *> *src = [img pixelSpectra];
    NSMutableArray<TTIOImzMLPixelSpectrum *> *pixels =
        [NSMutableArray arrayWithCapacity:src.count];
    for (TTIOPixelSpectrum *p in src) {
        NSError *pe = nil;
        TTIOImzMLPixelSpectrum *px =
            [[TTIOImzMLPixelSpectrum alloc] initWithX:(NSInteger)p.x
                                                    y:(NSInteger)p.y
                                                    z:(NSInteger)p.z
                                              mzArray:p.mz
                                       intensityArray:p.intensity
                                                error:&pe];
        if (!px) { if (error) *error = pe; return NO; }
        [pixels addObject:px];
    }

    NSString *scanPattern = img.scanPattern.length ? img.scanPattern : @"flyback";
    // Java GUI: ImzMLWriter.write(pixels, targetPath, ibd, mode, width,
    //              height, 1, pixelSizeX, pixelSizeY, scanPattern, null, sink).
    TTIOImzMLWriteResult *r =
        [TTIOImzMLWriter writePixels:pixels
                        toImzMLPath:output
                            ibdPath:ibd
                               mode:mode
                           gridMaxX:(NSInteger)img.width
                           gridMaxY:(NSInteger)img.height
                           gridMaxZ:1
                         pixelSizeX:img.pixelSizeX
                         pixelSizeY:img.pixelSizeY
                        scanPattern:scanPattern
                            uuidHex:nil
                              error:error];
    return r != nil;
}
@end


#pragma mark - JCAMP-DX (analytical run, first spectrum; encoding dispatch)

/* Parse the encoding opt (default affn). Mirrors Java
 * JcampDxEncoding.fromString. */
static TTIOJcampDxEncoding _owaJcampEncoding(NSString *s)
{
    NSString *low = [(s ?: @"affn") lowercaseString];
    if ([low isEqualToString:@"pac"]) return TTIOJcampDxEncodingPAC;
    if ([low isEqualToString:@"sqz"]) return TTIOJcampDxEncodingSQZ;
    if ([low isEqualToString:@"dif"]) return TTIOJcampDxEncodingDIF;
    return TTIOJcampDxEncodingAFFN;
}

@implementation TTIOJcampDxWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    TTIOJcampDxEncoding enc =
        _owaJcampEncoding(_owaOptString(options, @"encoding"));
    TTIOAcquisitionRun *run =
        [TTIORunSelection analyticalRunIn:dataset layer:layer error:error];
    if (!run) return NO;

    NSArray *spectra = [run spectra];
    if (spectra.count == 0) {
        if (error) *error = _owaError([NSString stringWithFormat:
            @"run %@ has no spectra", _owaLayerRepr(layer)]);
        return NO;
    }
    id first = spectra.firstObject;
    NSString *title = dataset.title ?: @"";
    if ([first isKindOfClass:[TTIOIRSpectrum class]]) {
        return [TTIOJcampDxWriter writeIRSpectrum:(TTIOIRSpectrum *)first
                                          toPath:output
                                           title:title
                                        encoding:enc
                                           error:error];
    } else if ([first isKindOfClass:[TTIORamanSpectrum class]]) {
        return [TTIOJcampDxWriter writeRamanSpectrum:(TTIORamanSpectrum *)first
                                             toPath:output
                                              title:title
                                           encoding:enc
                                              error:error];
    } else if ([first isKindOfClass:[TTIOUVVisSpectrum class]]) {
        return [TTIOJcampDxWriter writeUVVisSpectrum:(TTIOUVVisSpectrum *)first
                                             toPath:output
                                              title:title
                                           encoding:enc
                                              error:error];
    }
    if (error) *error = _owaError([NSString stringWithFormat:
        @"run %@ is %@, not a vibrational (IR/Raman/UV-Vis) spectrum",
        _owaLayerRepr(layer), NSStringFromClass([first class])]);
    return NO;
}
@end


#pragma mark - ISA (.json => ISA-JSON; else ISA-Tab bundle directory)

@implementation TTIOIsaWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    (void)layer; (void)options;
    NSString *lower = [output.lastPathComponent lowercaseString];
    if ([lower hasSuffix:@".json"]) {
        // ISA-JSON: serialise the bundle's investigation.json to the file.
        NSDictionary<NSString *, NSData *> *bundle =
            [TTIOISAExporter bundleForDataset:dataset error:error];
        if (!bundle) return NO;
        NSData *json = bundle[@"investigation.json"];
        if (!json) {
            if (error) *error =
                _owaError(@"ISA bundle has no investigation.json");
            return NO;
        }
        return [json writeToFile:output options:NSDataWritingAtomic error:error];
    }
    // Directory-style ISA-Tab bundle.
    return [TTIOISAExporter writeBundleForDataset:dataset
                                     toDirectory:output
                                           error:error];
}
@end


#pragma mark - BAM (genomic run -> written; sorted; provenance)

@implementation TTIOBamWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    (void)options;
    TTIOGenomicRun *run =
        [TTIORunSelection genomicRunIn:dataset layer:layer error:error];
    if (!run) return NO;
    // Java GUI: new BamWriter(targetPath).write(run, provenance, /*sort=*/true, sink).
    TTIOBamWriter *w = [[TTIOBamWriter alloc] initWithPath:output];
    return [w writeReadSideRun:run
             provenanceRecords:dataset.provenanceRecords
                          sort:YES
                         error:error];
}
@end


#pragma mark - CRAM (requires reference; reference-compressed)

@implementation TTIOCramWriterAdapter
- (BOOL)writeDataset:(TTIOSpectralDataset *)dataset
               layer:(nullable NSString *)layer
            toOutput:(NSString *)output
             options:(NSDictionary<NSString *, id> *)options
               error:(NSError *_Nullable *_Nullable)error
{
    NSString *reference = _owaOptString(options, @"reference");
    if (reference == nil) {
        if (error) *error = _owaError(
            @"CRAM export is reference-compressed; pass the reference FASTA "
            @"via --extra --reference <path>");
        return NO;
    }
    TTIOGenomicRun *run =
        [TTIORunSelection genomicRunIn:dataset layer:layer error:error];
    if (!run) return NO;
    // Java GUI: new CramWriter(targetPath, reference)
    //              .write(run, provenance, /*sort=*/true, sink).
    TTIOCramWriter *w =
        [[TTIOCramWriter alloc] initWithPath:output referenceFasta:reference];
    return [w writeReadSideRun:run
             provenanceRecords:dataset.provenanceRecords
                          sort:YES
                         error:error];
}
@end
