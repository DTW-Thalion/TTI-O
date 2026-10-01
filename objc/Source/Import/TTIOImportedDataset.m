/*
 * TTIOImportedDataset.m
 * TTI-O Objective-C Implementation
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "Import/TTIOImportedDataset.h"
#import "Import/TTIOGenomicStreamSource.h"
#import "Import/TTIOSpectralStreamSource.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "HDF5/TTIOHDF5Group.h"
#import "HDF5/TTIOHDF5Errors.h"

@implementation TTIOImportedDataset

- (instancetype)init
{
    self = [super init];
    if (self) {
        _title = @"";
        _isaInvestigationId = @"";
        _msRuns = [[NSMutableDictionary alloc] init];
        _genomicRuns = [[NSMutableDictionary alloc] init];
        _genomicStreams = [[NSMutableDictionary alloc] init];
        _spectralStreams = [[NSMutableDictionary alloc] init];
        _identifications = [[NSMutableArray alloc] init];
        _quantifications = [[NSMutableArray alloc] init];
        _provenanceRecords = [[NSMutableArray alloc] init];
    }
    return self;
}

+ (instancetype)datasetWithWriteDelegate:(TTIOImportedDatasetWriteDelegate)delegate
{
    TTIOImportedDataset *ds = [[self alloc] init];
    ds.writeDelegate = delegate;
    return ds;
}

- (BOOL)writeToPath:(NSString *)path error:(NSError *_Nullable *_Nullable)error
{
    return [self writeToPath:path progress:nil error:error];
}

- (BOOL)writeToPath:(NSString *)path
           progress:(TTIOProgressBlock)progress
              error:(NSError *_Nullable *_Nullable)error
{
    // Write-through delegate (subprocess / image importers) wins.
    if (self.writeDelegate) {
        return self.writeDelegate(path, error);
    }
    if (![self _writeStaticToPath:path error:error]) return NO;
    if (self.genomicStreams.count == 0 && self.spectralStreams.count == 0) return YES;
    return [self _writeStreamsToPath:path progress:progress error:error];
}

/* The streams are appended through a read-write reopen of the file
 * the static write produced, the genomic feature flag added when a
 * genomic stream is the first genomic content. */
- (BOOL)_writeStreamsToPath:(NSString *)path
                   progress:(TTIOProgressBlock)progress
                      error:(NSError *_Nullable *_Nullable)error
{
    id<TTIOStorageProvider> p = [[TTIOProviderRegistry sharedRegistry]
        openURL:path mode:TTIOStorageOpenModeReadWrite provider:@"hdf5" error:error];
    if (!p) return NO;
    BOOL ok = YES;
    @try {
        id<TTIOStorageGroup> root = [p rootGroupWithError:error];
        id<TTIOStorageGroup> study = root ? [root openGroupNamed:@"study" error:error] : nil;
        if (!study) { ok = NO; return NO; }
        if (self.genomicStreams.count > 0 && [root respondsToSelector:@selector(unwrap)]) {
            TTIOHDF5Group *h5root = [(id)root performSelector:@selector(unwrap)];
            NSString *flag = [TTIOFeatureFlags featureOptGenomic];
            if (h5root && ![TTIOFeatureFlags root:h5root supportsFeature:flag]) {
                NSMutableArray *features =
                    [[TTIOFeatureFlags featuresForRoot:h5root] mutableCopy] ?: [NSMutableArray array];
                [features addObject:flag];
                NSString *version = [TTIOFeatureFlags formatVersionForRoot:h5root];
                if (![TTIOFeatureFlags writeFormatVersion:version features:features
                                                    toRoot:h5root error:error]) { ok = NO; return NO; }
            }
        }
        NSArray *gNames = [[self.genomicStreams allKeys] sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *name in gNames) {
            if ([self.genomicStreams[name] writeIntoStudy:study progress:progress error:error] == NSNotFound) {
                ok = NO; return NO;
            }
        }
        // M101: opt_sam_tags once a streamed genomic run turned out to
        // carry SAM tags; the static write could not know.
        if (gNames.count > 0 && [root respondsToSelector:@selector(unwrap)]) {
            TTIOHDF5Group *h5root = [(id)root performSelector:@selector(unwrap)];
            id<TTIOStorageGroup> runs = [study hasChildNamed:@"genomic_runs"]
                ? [study openGroupNamed:@"genomic_runs" error:NULL] : nil;
            BOOL anyTags = NO;
            for (NSString *key in gNames) {
                NSString *runName = self.genomicStreams[key].name ?: key;
                if (![runs hasChildNamed:runName]) continue;
                id<TTIOStorageGroup> sc = [[runs openGroupNamed:runName error:NULL]
                                           openGroupNamed:@"signal_channels" error:NULL];
                if ([sc hasChildNamed:@"tags"]) { anyTags = YES; break; }
            }
            NSString *flag = [TTIOFeatureFlags featureOptSamTags];
            if (anyTags && h5root && ![TTIOFeatureFlags root:h5root supportsFeature:flag]) {
                NSMutableArray *features =
                    [[TTIOFeatureFlags featuresForRoot:h5root] mutableCopy] ?: [NSMutableArray array];
                [features addObject:flag];
                NSString *version = [TTIOFeatureFlags formatVersionForRoot:h5root];
                if (![TTIOFeatureFlags writeFormatVersion:version features:features
                                                    toRoot:h5root error:error]) { ok = NO; return NO; }
            }
        }
        NSArray *sNames = [[self.spectralStreams allKeys] sortedArrayUsingSelector:@selector(compare:)];
        for (NSString *name in sNames) {
            if ([self.spectralStreams[name] writeIntoStudy:study progress:progress error:error] == NSNotFound) {
                ok = NO; return NO;
            }
        }
    } @finally {
        [p close];
    }
    return ok;
}

- (BOOL)_writeStaticToPath:(NSString *)path error:(NSError *_Nullable *_Nullable)error
{

    // In-memory draft -> canonical mixed-dictionary write API. Images are
    // NOT written here (writeMinimal has no image parameter); they are
    // persisted by the imzML adapter's write-through delegate (OT5).
    NSString *title = (self.title.length ? self.title : @"imported");
    NSDictionary<NSString *, id> *genomic =
        (self.genomicRuns.count ? self.genomicRuns : nil);
    NSArray *idents = (self.identifications.count ? self.identifications : nil);
    NSArray *quants = (self.quantifications.count ? self.quantifications : nil);
    NSArray *prov   = (self.provenanceRecords.count ? self.provenanceRecords : nil);

    return [TTIOSpectralDataset writeMinimalToPath:path
                                             title:title
                                isaInvestigationId:self.isaInvestigationId
                                         mixedRuns:self.msRuns
                                       genomicRuns:genomic
                                   identifications:idents
                                   quantifications:quants
                                 provenanceRecords:prov
                                             error:error];
}

@end
