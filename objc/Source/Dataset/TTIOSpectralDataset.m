/*
 * TTIOSpectralDataset.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOSpectralDataset
 * Inherits From: NSObject
 * Conforms To:   TTIOEncryptable
 * Declared In:   Dataset/TTIOSpectralDataset.h
 *
 * Root container for a .tio file. Owns the top-level study/ group
 * plus per-modality run dictionaries (MS / NMR / genomic). Provides
 * the +writeMinimalToPath: flat-buffer write paths and the
 * +decryptInPlaceAtPath: in-place decryption helper.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#include <pthread.h>
#import "TTIOSpectralDataset.h"
#import "TTIOWrittenRun.h"
#import "TTIOIdentification.h"
#import "TTIOQuantification.h"
#import "TTIOProvenanceRecord.h"
#import "TTIOTransitionList.h"
#import "TTIOCompoundIO.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOSpectrumIndex.h"
#import "Spectra/TTIONMRSpectrum.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "HDF5/TTIOHDF5Dataset.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOHDF5Types.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "Protection/TTIOEncryptionManager.h"
#import "Protection/TTIOAccessPolicy.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOCompoundField.h"
#import "Providers/TTIOHDF5Provider.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOReferenceImport.h"
#import "Assembly/TTIOAssemblyGraph.h"          // M98
#import "Genomics/TTIOBulkV2Blobs.h"           // Phase 2c-T
#import "Codecs/TTIORans.h"
#import "Codecs/TTIOBasePack.h"
#import "Codecs/TTIOQuality.h"                 // M86 Phase D
#import "Codecs/TTIOFqzcompNx16Z.h"             // M94.Z v1.2
#import "Codecs/TTIODeltaRans.h"                // M95 v1.2
#import "Codecs/TTIOMateInfoV2.h"               // inline mate-pair codec
#import "Codecs/TTIORefDiffV2.h"               // bit-packed ref-diff v2
#import "Codecs/TTIONameTokenizerV2.h"          // v1.8 #11 ch3: adaptive name-tokenizer v2
#import "Codecs/Registry/TTIOCodecRegistry.h"   // Task 6: codec-registry routing
#import "Codecs/Registry/TTIOCodec.h"
#import "Codecs/Registry/TTIODecodedChannel.h"
#import "Codecs/Registry/TTIOEncodedChannel.h"
#import "Codecs/Registry/TTIOCodecContext.h"
#import <hdf5.h>
#include <objc/message.h>                          // Task 6: typed objc_msgSend bridge
#include <openssl/md5.h>                          // M93 v1.2 ref MD5
#import "TTIOSpectralDataset+Internal.h"          // P3.10: shared with +GenomicWrite category

// Bridge to the dynamically-resolved TTIOHDF5GroupAdapter. The adapter
// lives in a separate compilation unit and is looked up by name to
// avoid a hard build dependency. `initWithGroup:` is an init-family
// selector, so a direct -performSelector: trips ARC's hard
// "performSelector names a selector which retains the object" error
// once the selector is in scope. We dispatch through a typed
// objc_msgSend pointer instead: identical runtime behaviour (the
// alloc/init pair is unchanged) but ARC tracks ownership normally.
// Non-static: shared with the +GenomicWrite category via +Internal.h.
id _TTIO_MakeHDF5GroupAdapter(id group)
{
    Class cls = NSClassFromString(@"TTIOHDF5GroupAdapter");
    if (cls == Nil) return nil;
    id obj = [cls alloc];
    SEL sel = NSSelectorFromString(@"initWithGroup:");
    id (*msg)(id, SEL, id) = (id (*)(id, SEL, id))objc_msgSend;
    return msg(obj, sel, group);
}

// M98: enumerate /study/assembly_graphs/ (storage-protocol group)
// into name -> TTIOAssemblyGraph. Empty dict when the subtree is
// absent or malformed; individual unreadable graphs are skipped the
// way unreadable genomic runs are.
static NSDictionary<NSString *, TTIOAssemblyGraph *> *
_ttio_loadAssemblyGraphs(id<TTIOStorageGroup> study)
{
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    if (![study hasChildNamed:@"assembly_graphs"]) return out;
    id<TTIOStorageGroup> ag = [study openGroupNamed:@"assembly_graphs"
                                              error:NULL];
    if (!ag) return out;
    id namesObj = [ag attributeValueForName:@"_graph_names" error:NULL];
    if (![namesObj isKindOfClass:[NSString class]]) return out;
    for (NSString *gn in [(NSString *)namesObj
             componentsSeparatedByString:@","]) {
        NSString *trimmed = [gn stringByTrimmingCharactersInSet:
            [NSCharacterSet whitespaceCharacterSet]];
        if (trimmed.length == 0 || ![ag hasChildNamed:trimmed]) continue;
        id<TTIOStorageGroup> graphG = [ag openGroupNamed:trimmed
                                                   error:NULL];
        if (!graphG) continue;
        TTIOAssemblyGraph *graph = [TTIOAssemblyGraph
            openFromGroup:graphG name:trimmed error:NULL];
        if (graph) out[trimmed] = graph;
    }
    return out;
}

// Internal SPI surfaced by TTIOAcquisitionRun for the dataset-level
// decrypt lifecycle. Not part of the public header.
@interface TTIOAcquisitionRun (TTIOSpectralDatasetInternal)
- (NSData *)decryptedChannelNamed:(NSString *)chName;
- (BOOL)reattachSignalHandlesFromGroup:(id<TTIOStorageGroup>)channels error:(NSError **)error;
@end

// v1.0 single format-version stamp. Readers gate optional features
// by the feature-flag list (opt_*), not by version equality.
// Non-static: shared with the +GenomicWrite category via +Internal.h.
NSString *const kTTIOFormatVersion = @"1.0";

/** v0.12 M74 Slice E: scan the ms_runs dict for any run whose
 *  spectrum_index carries the four optional activation/isolation
 *  columns. When present, the writer upgrades the feature flag list
 *  with opt_ms2_activation_detail and bumps the on-disk format version
 *  to 1.3. Returns NO when every run has the legacy layout. */
static BOOL datasetRunsHaveActivationDetail(NSDictionary *msRuns)
{
    for (TTIOAcquisitionRun *run in [msRuns objectEnumerator]) {
        if (run.spectrumIndex.hasActivationDetail) return YES;
    }
    return NO;
}

/** M102: YES when any MS run's spectrum_index carries the pixel_x /
 *  pixel_y / pixel_z columns; the writer then adds
 *  opt_pixel_coordinates. */
static BOOL datasetRunsHavePixelCoordinates(NSDictionary *msRuns)
{
    for (TTIOAcquisitionRun *run in [msRuns objectEnumerator]) {
        if (run.spectrumIndex.hasPixelCoordinates) return YES;
    }
    return NO;
}


@implementation TTIOSpectralDataset
{
    TTIOHDF5File     *_file;       // retained while alive for lazy reads
    NSString         *_filePath;
    TTIOAccessPolicy *_accessPolicy;
    NSString         *_encryptedAlgorithm;  // empty string when not encrypted
    id<TTIOStorageProvider> _provider;  // owns _file
}

@synthesize filePath = _filePath;
@synthesize provider = _provider;
@synthesize encryptedAlgorithm = _encryptedAlgorithm;
@synthesize genomicRuns = _genomicRuns;
@synthesize assemblyGraphs = _assemblyGraphs;
@synthesize references = _references;

- (BOOL)isEncrypted
{
    return _encryptedAlgorithm.length > 0;
}

- (instancetype)initWithTitle:(NSString *)title
           isaInvestigationId:(NSString *)isaId
                       msRuns:(NSDictionary *)msRuns
                      nmrRuns:(NSDictionary *)nmrRuns
              identifications:(NSArray *)identifications
              quantifications:(NSArray *)quantifications
            provenanceRecords:(NSArray *)provenance
                  transitions:(TTIOTransitionList *)transitions
{
    self = [super init];
    if (self) {
        _title              = [title copy];
        _isaInvestigationId = [isaId copy];
        _msRuns             = [msRuns copy] ?: @{};
        _nmrRuns            = [nmrRuns copy] ?: @{};
        _genomicRuns        = @{};   // populated by +readFromFilePath: when present
        _assemblyGraphs     = @{};   // populated by +readFromFilePath: when present
        _references         = @{};   // populated by +readFromFilePath: when present
        _identifications    = [identifications copy] ?: @[];
        _quantifications    = [quantifications copy] ?: @[];
        _provenanceRecords  = [provenance copy] ?: @[];
        _transitions        = transitions;
        _encryptedAlgorithm = @"";
    }
    return self;
}

- (void)dealloc
{
    [self closeFile];
}

#pragma mark - Access policy JSON helpers

static NSString *encodeAccessPolicy(TTIOAccessPolicy *p)
{
    if (!p || !p.policy) return nil;
    NSData *d = [NSJSONSerialization dataWithJSONObject:p.policy options:0 error:NULL];
    if (!d) return nil;
    return [[NSString alloc] initWithData:d encoding:NSUTF8StringEncoding];
}

static TTIOAccessPolicy *decodeAccessPolicy(NSString *json)
{
    if (json.length == 0) return nil;
    NSData *d = [json dataUsingEncoding:NSUTF8StringEncoding];
    NSDictionary *parsed = [NSJSONSerialization JSONObjectWithData:d options:0 error:NULL];
    if (![parsed isKindOfClass:[NSDictionary class]]) return nil;
    return [[TTIOAccessPolicy alloc] initWithPolicy:parsed];
}

#pragma mark - JSON-plist helpers (v0.1 fallback only)

static NSArray *decodePlistArray(NSString *json, Class cls, NSError **error)
{
    if (!json) return @[];
    NSData *data = [json dataUsingEncoding:NSUTF8StringEncoding];
    NSArray *plists = [NSJSONSerialization JSONObjectWithData:data options:0 error:error];
    if (!plists) return nil;
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:plists.count];
    for (NSDictionary *p in plists) [out addObject:[cls fromPlist:p]];
    return out;
}

#pragma mark - URL scheme routing (v0.9 M64.5)

// Non-static: shared with the +GenomicWrite category via +Internal.h.
BOOL isNonHdf5ProviderURL(NSString *url) {
    if (url.length == 0) return NO;
    return [url hasPrefix:@"memory://"]
        || [url hasPrefix:@"sqlite://"]
        || [url hasPrefix:@"zarr://"];
}

static NSError *makeProviderWriteNotImplementedError(NSString *url) {
    NSString *msg = [NSString stringWithFormat:
        @"ObjC SpectralDataset *write* via URL '%@' not implemented "
        @"in v0.9 (read is supported via +readViaProviderURL:). "
        @"Produce non-HDF5 .tio files through Python / Java which "
        @"have the full write-side caller refactor.",
        url];
    return [NSError errorWithDomain:@"TTIOSpectralDatasetErrorDomain"
                                code:999
                            userInfo:@{NSLocalizedDescriptionKey: msg}];
}

#pragma mark - HDF5 write

- (BOOL)writeToFilePath:(NSString *)path error:(NSError **)error
{
    // Task 31: narrow the non-HDF5 rejection. MS-only datasets without
    // Image-subclass hooks can write to memory/sqlite/zarr URLs through
    // the protocol path. NMR runs and Image subclasses use HDF5-direct
    // features (TTIONMR2DSpectrum H5DSset_scale dimension scales,
    // MSImage 3-D cube via H5Pset_chunk) that don't have protocol
    // equivalents — those still require an HDF5 backing file.
    if (isNonHdf5ProviderURL(path)) {
        BOOL hasNmrRuns      = (_nmrRuns.count > 0);
        BOOL isImageSubclass = ([self class] != [TTIOSpectralDataset class]);
        if (hasNmrRuns || isImageSubclass) {
            if (error) *error = makeProviderWriteNotImplementedError(path);
            return NO;
        }
    }

    // Open the appropriate provider based on URL scheme. Plain paths
    // and file:// URLs go to TTIOHDF5Provider; memory:// / sqlite:// /
    // zarr:// to their respective providers (registered via +load).
    id<TTIOStorageProvider> provider =
        [[TTIOProviderRegistry sharedRegistry] openURL:path
                                                  mode:TTIOStorageOpenModeCreate
                                              provider:nil
                                                 error:error];
    if (!provider) return NO;

    // For HDF5, use the raw TTIOHDF5Group (which conforms to
    // <TTIOStorageGroup> via Task 31 bridge methods). This keeps Image
    // subclass and 2D-NMR HDF5-direct features (H5DSset_scale, native
    // 3D cubes) reachable without wrapper-unwrap dance. Non-HDF5
    // providers continue to expose their group via rootGroupWithError:.
    id<TTIOStorageGroup> root = nil;
    if ([provider isKindOfClass:[TTIOHDF5Provider class]]) {
        TTIOHDF5File *f = (TTIOHDF5File *)[provider nativeHandle];
        if (f) root = [f rootGroup];
    } else {
        root = [provider rootGroupWithError:error];
    }
    if (!root) { [provider close]; return NO; }

    // Emit v0.2 format + feature flags. The per-run compound provenance
    // flag (M17) is emitted unconditionally: every v0.3 writer produces
    // compound-form per-run provenance when any run carries records, and
    // the flag advertises that capability to future readers even when the
    // current in-memory dataset happens to have no provenance to persist.
    NSMutableArray *features = [@[
        [TTIOFeatureFlags featureBaseV1],
        [TTIOFeatureFlags featureCompoundIdentifications],
        [TTIOFeatureFlags featureCompoundQuantifications],
        [TTIOFeatureFlags featureCompoundProvenance],
        [TTIOFeatureFlags featureCompoundPerRunProvenance],
        [TTIOFeatureFlags featureCompoundHeaders],
        [TTIOFeatureFlags featureNative2DNMR],
        [TTIOFeatureFlags featureNativeMSImageCube],
    ] mutableCopy];
    BOOL anyM74 = datasetRunsHaveActivationDetail(_msRuns);
    if (anyM74) {
        [features addObject:[TTIOFeatureFlags featureMS2ActivationDetail]];
    }
    if (datasetRunsHavePixelCoordinates(_msRuns)) {
        [features addObject:[TTIOFeatureFlags featurePixelCoordinates]];
    }
    if (![root setAttributeValue:kTTIOFormatVersion
                         forName:@"ttio_format_version" error:error]) {
        [provider close]; return NO;
    }
    NSData *featJSON =
        [NSJSONSerialization dataWithJSONObject:features options:0 error:NULL];
    NSString *featStr =
        [[NSString alloc] initWithData:featJSON encoding:NSUTF8StringEncoding];
    if (![root setAttributeValue:featStr
                         forName:@"ttio_features" error:error]) {
        [provider close]; return NO;
    }

    // Access policy, if set.
    NSString *apJson = encodeAccessPolicy(_accessPolicy);
    if (apJson) {
        if (![root setAttributeValue:apJson
                             forName:@"access_policy_json" error:error]) {
            [provider close]; return NO;
        }
    }

    id<TTIOStorageGroup> study = [root createGroupNamed:@"study" error:error];
    if (!study) { [provider close]; return NO; }
    if (![study setAttributeValue:(_title ?: @"")
                          forName:@"title" error:error]) {
        [provider close]; return NO;
    }
    if (![study setAttributeValue:(_isaInvestigationId ?: @"")
                          forName:@"isa_investigation_id" error:error]) {
        [provider close]; return NO;
    }

    // MS runs
    id<TTIOStorageGroup> msRunsGroup = [study createGroupNamed:@"ms_runs" error:error];
    if (!msRunsGroup) { [provider close]; return NO; }
    NSArray *msNames = [[_msRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];
    if (![msRunsGroup setAttributeValue:[msNames componentsJoinedByString:@","]
                                forName:@"_run_names" error:error]) {
        [provider close]; return NO;
    }
    for (NSString *runName in msNames) {
        TTIOAcquisitionRun *run = _msRuns[runName];
        if (![run writeToGroup:msRunsGroup name:runName error:error]) {
            [provider close]; return NO;
        }

        // Write compound headers alongside the parallel index datasets.
        // HDF5-only feature (h5dump readability via VL compound type);
        // skip silently for non-HDF5 providers.
        id<TTIOStorageGroup> runG = [msRunsGroup openGroupNamed:runName error:NULL];
        id<TTIOStorageGroup> idxG = [runG openGroupNamed:@"spectrum_index" error:NULL];
        if (idxG && [idxG isKindOfClass:[TTIOHDF5Group class]]) {
            [TTIOCompoundIO writeCompoundHeadersForIndex:run.spectrumIndex
                                                intoGroup:(TTIOHDF5Group *)idxG
                                                    error:NULL];
        }
    }

    // NMR runs (legacy nmrRuns dict, kept for backward compat). The
    // writer above guards non-HDF5 URLs from reaching here when nmrRuns
    // is non-empty, so this path is effectively HDF5-only. We still use
    // protocol methods for byte-exact symmetry with the MS path.
    id<TTIOStorageGroup> nmrRunsGroup = [study createGroupNamed:@"nmr_runs" error:error];
    if (!nmrRunsGroup) { [provider close]; return NO; }
    NSArray *nmrNames = [[_nmrRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];
    if (![nmrRunsGroup setAttributeValue:[nmrNames componentsJoinedByString:@","]
                                 forName:@"_run_names" error:error]) {
        [provider close]; return NO;
    }
    for (NSString *runName in nmrNames) {
        id<TTIOStorageGroup> nmrRun = [nmrRunsGroup createGroupNamed:runName error:error];
        if (!nmrRun) { [provider close]; return NO; }
        NSArray<TTIONMRSpectrum *> *spectra = _nmrRuns[runName];
        if (![nmrRun setAttributeValue:@((int64_t)spectra.count)
                               forName:@"count" error:error]) {
            [provider close]; return NO;
        }
        for (NSUInteger i = 0; i < spectra.count; i++) {
            NSString *name = [NSString stringWithFormat:@"spec_%06lu", (unsigned long)i];
            if (![spectra[i] writeToGroup:nmrRun name:name error:error]) {
                [provider close]; return NO;
            }
        }
    }

    // Compound identifications / quantifications / provenance.
    // Task 31: TTIOCompoundIO routes HDF5 through its fast path and
    // non-HDF5 through createCompoundDatasetNamed: + writeAll:.
    if (_identifications.count > 0) {
        if (![TTIOCompoundIO writeIdentifications:_identifications
                                         intoGroup:study
                                      datasetNamed:@"identifications"
                                             error:error]) {
            [provider close]; return NO;
        }
    }
    if (_quantifications.count > 0) {
        if (![TTIOCompoundIO writeQuantifications:_quantifications
                                         intoGroup:study
                                      datasetNamed:@"quantifications"
                                             error:error]) {
            [provider close]; return NO;
        }
    }
    if (_provenanceRecords.count > 0) {
        if (![TTIOCompoundIO writeProvenance:_provenanceRecords
                                    intoGroup:study
                                 datasetNamed:@"provenance"
                                        error:error]) {
            [provider close]; return NO;
        }
    }

    if (_transitions) {
        NSData *tdata = [NSJSONSerialization dataWithJSONObject:[_transitions asPlist]
                                                        options:0
                                                          error:error];
        if (!tdata) { [provider close]; return NO; }
        NSString *tjson = [[NSString alloc] initWithData:tdata encoding:NSUTF8StringEncoding];
        if (![study setAttributeValue:tjson
                              forName:@"transitions_json" error:error]) {
            [provider close]; return NO;
        }
    }

    _filePath = [path copy];
    [provider close];
    return YES;
}



#pragma mark - HDF5 read

+ (instancetype)readViaProviderURL:(NSString *)url
                              mode:(TTIOStorageOpenMode)mode
                             error:(NSError **)error
{
    // v0.9 M64.5-objc-java: read a non-HDF5 .tio by routing through
    // the provider registry. Metadata (idents/quants/prov) comes from
    // the JSON mirror attributes; runs are reconstructed via
    // +[TTIOAcquisitionRun readFromStorageGroup:].
    id<TTIOStorageProvider> prov = [[TTIOProviderRegistry sharedRegistry]
        openURL:url mode:mode provider:nil error:error];
    if (!prov) return nil;
    id<TTIOStorageGroup> root = [prov rootGroupWithError:error];
    if (!root) return nil;

    NSString *title = @"", *isaId = @"";
    NSMutableDictionary *msRuns = [NSMutableDictionary dictionary];
    NSMutableDictionary *genomicRunsMap = [NSMutableDictionary dictionary];
    NSMutableDictionary *assemblyGraphsMap = [NSMutableDictionary dictionary];
    NSMutableDictionary<NSString *, TTIOReferenceImport *> *refsMap =
        [NSMutableDictionary dictionary];
    NSArray *idents = @[], *quants = @[], *provRecs = @[];

    if ([root hasChildNamed:@"study"]) {
        id<TTIOStorageGroup> study = [root openGroupNamed:@"study" error:error];
        if (!study) return nil;

        id titleObj = [study attributeValueForName:@"title" error:NULL];
        if ([titleObj isKindOfClass:[NSString class]]) title = titleObj;
        id isaObj = [study attributeValueForName:@"isa_investigation_id" error:NULL];
        if ([isaObj isKindOfClass:[NSString class]]) isaId = isaObj;

        if ([study hasChildNamed:@"ms_runs"]) {
            id<TTIOStorageGroup> ms = [study openGroupNamed:@"ms_runs" error:NULL];
            id namesObj = [ms attributeValueForName:@"_run_names" error:NULL];
            if ([namesObj isKindOfClass:[NSString class]]) {
                for (NSString *rn in [(NSString *)namesObj componentsSeparatedByString:@","]) {
                    if (rn.length == 0) continue;
                    TTIOAcquisitionRun *run = [TTIOAcquisitionRun readFromStorageGroup:ms
                                                                                   name:rn
                                                                                  error:NULL];
                    if (run) msRuns[rn] = run;
                }
            }
        }

        // provider-agnostic genomic_runs read.
        if ([study hasChildNamed:@"genomic_runs"]) {
            id<TTIOStorageGroup> gG = [study openGroupNamed:@"genomic_runs" error:NULL];
            id namesObj = [gG attributeValueForName:@"_run_names" error:NULL];
            if ([namesObj isKindOfClass:[NSString class]]) {
                for (NSString *rn in [(NSString *)namesObj componentsSeparatedByString:@","]) {
                    NSString *trimmed = [rn stringByTrimmingCharactersInSet:
                        [NSCharacterSet whitespaceCharacterSet]];
                    if (trimmed.length == 0 || ![gG hasChildNamed:trimmed]) continue;
                    id<TTIOStorageGroup> runG = [gG openGroupNamed:trimmed error:NULL];
                    TTIOGenomicRun *gr = [TTIOGenomicRun openFromGroup:runG
                                                                    name:trimmed
                                                                   error:NULL];
                    if (gr) genomicRunsMap[trimmed] = gr;
                }
            }
        }

        // provider-agnostic assembly_graphs read (M98).
        [assemblyGraphsMap addEntriesFromDictionary:
            _ttio_loadAssemblyGraphs(study)];

        // /study/references/<uri>/ — embedded references read-back
        // (1.1.0 — Phase 0 tio-browser). Empty dict when absent.
        if ([study hasChildNamed:@"references"]) {
            id<TTIOStorageGroup> refsG =
                [study openGroupNamed:@"references" error:NULL];
            if (refsG != nil) {
                for (NSString *uri in [refsG childNames]) {
                    id<TTIOStorageGroup> oneRefG =
                        [refsG openGroupNamed:uri error:NULL];
                    if (oneRefG == nil) continue;
                    TTIOReferenceImport *r =
                        [TTIOReferenceImport readFromGroup:oneRefG];
                    if (r != nil) refsMap[uri] = r;
                }
            }
        }

        id iObj = [study attributeValueForName:@"identifications_json" error:NULL];
        if ([iObj isKindOfClass:[NSString class]]) {
            NSArray *plists = [NSJSONSerialization
                JSONObjectWithData:[(NSString *)iObj dataUsingEncoding:NSUTF8StringEncoding]
                           options:0 error:NULL];
            NSMutableArray *arr = [NSMutableArray array];
            for (NSDictionary *d in plists) {
                id rec = [TTIOIdentification fromPlist:d];
                if (rec) [arr addObject:rec];
            }
            idents = arr;
        }
        id qObj = [study attributeValueForName:@"quantifications_json" error:NULL];
        if ([qObj isKindOfClass:[NSString class]]) {
            NSArray *plists = [NSJSONSerialization
                JSONObjectWithData:[(NSString *)qObj dataUsingEncoding:NSUTF8StringEncoding]
                           options:0 error:NULL];
            NSMutableArray *arr = [NSMutableArray array];
            for (NSDictionary *d in plists) {
                id rec = [TTIOQuantification fromPlist:d];
                if (rec) [arr addObject:rec];
            }
            quants = arr;
        }
        id pObj = [study attributeValueForName:@"provenance_json" error:NULL];
        if ([pObj isKindOfClass:[NSString class]]) {
            NSArray *plists = [NSJSONSerialization
                JSONObjectWithData:[(NSString *)pObj dataUsingEncoding:NSUTF8StringEncoding]
                           options:0 error:NULL];
            NSMutableArray *arr = [NSMutableArray array];
            for (NSDictionary *d in plists) {
                id rec = [TTIOProvenanceRecord fromPlist:d];
                if (rec) [arr addObject:rec];
            }
            provRecs = arr;
        }
    }

    TTIOSpectralDataset *ds = [[self alloc] initWithTitle:title
                                        isaInvestigationId:isaId
                                                    msRuns:msRuns
                                                   nmrRuns:@{}
                                           identifications:idents
                                           quantifications:quants
                                         provenanceRecords:provRecs
                                               transitions:nil];
    ds->_filePath    = [url copy];
    ds->_genomicRuns = [genomicRunsMap copy];
    ds->_assemblyGraphs = [assemblyGraphsMap copy];
    ds->_references  = [refsMap copy];
    // Surface the root `encrypted` attr for provider-backed reads too.
    id encObj = [root attributeValueForName:@"encrypted" error:NULL];
    if ([encObj isKindOfClass:[NSString class]]) {
        ds->_encryptedAlgorithm = [(NSString *)encObj copy];
    } else {
        ds->_encryptedAlgorithm = @"";
    }
    // _file / _provider stay nil — the provider instance was transient;
    // close() is a no-op for provider-backed datasets in v0.9.
    return ds;
}

+ (instancetype)readFromFilePath:(NSString *)path error:(NSError **)error
{
    return [self readFromFilePath:path writable:NO error:error];
}

+ (instancetype)readFromFilePath:(NSString *)path
                        writable:(BOOL)writable
                           error:(NSError **)error
{
    TTIOStorageOpenMode mode = writable ? TTIOStorageOpenModeReadWrite
                                        : TTIOStorageOpenModeRead;
    if (isNonHdf5ProviderURL(path)) {
        return [self readViaProviderURL:path mode:mode error:error];
    }
    // route through TTIOHDF5Provider; the native handle is the
    // TTIOHDF5File previously obtained directly.
    TTIOHDF5Provider *p = [[TTIOHDF5Provider alloc] init];
    if (![p openURL:path mode:mode error:error]) return nil;
    TTIOHDF5File *f = (TTIOHDF5File *)[p nativeHandle];
    if (!f) return nil;
    TTIOHDF5Group *root = [f rootGroup];

    BOOL isV1 = [TTIOFeatureFlags isLegacyV1File:root];

    TTIOHDF5Group *study = [root openGroupNamed:@"study" error:error];
    if (!study) return nil;

    NSString *title  = [study stringAttributeNamed:@"title" error:error];
    NSString *isaId  = [study stringAttributeNamed:@"isa_investigation_id" error:error];

    // MS runs (v0.2: any TTIOSpectrum subclass)
    NSMutableDictionary *msRuns = [NSMutableDictionary dictionary];
    if ([study hasChildNamed:@"ms_runs"]) {
        TTIOHDF5Group *msg = [study openGroupNamed:@"ms_runs" error:error];
        NSString *namesStr = [msg stringAttributeNamed:@"_run_names" error:error];
        for (NSString *rname in [namesStr componentsSeparatedByString:@","]) {
            if (rname.length == 0) continue;
            TTIOAcquisitionRun *run = [TTIOAcquisitionRun readFromGroup:msg name:rname error:error];
            if (!run) return nil;
            [run setPersistenceFilePath:path runName:rname];
            msRuns[rname] = run;
        }
    }

    // genomic_runs subtree (absent on pre-M82 files → empty dict).
    NSMutableDictionary *genomicRuns = [NSMutableDictionary dictionary];
    if ([study hasChildNamed:@"genomic_runs"]) {
        TTIOHDF5Group *gg = [study openGroupNamed:@"genomic_runs" error:error];
        NSString *gNames = [gg stringAttributeNamed:@"_run_names" error:NULL] ?: @"";
        for (NSString *rname in [gNames componentsSeparatedByString:@","]) {
            NSString *trimmed =
                [rname stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceCharacterSet]];
            if (trimmed.length == 0 || ![gg hasChildNamed:trimmed]) continue;
            TTIOHDF5Group *runG = [gg openGroupNamed:trimmed error:error];
            if (!runG) return nil;
            id<TTIOStorageGroup> runAdapter =
                (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(runG);
            TTIOGenomicRun *gr =
                [TTIOGenomicRun openFromGroup:runAdapter name:trimmed error:error];
            if (!gr) return nil;
            genomicRuns[trimmed] = gr;
        }
    }

    // assembly_graphs subtree (M98) — read through the storage
    // adapter so the same loader serves the provider-backed path.
    NSDictionary *assemblyGraphs = @{};
    if ([study hasChildNamed:@"assembly_graphs"]) {
        id<TTIOStorageGroup> studyAdapter =
            (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(study);
        if (studyAdapter) {
            assemblyGraphs = _ttio_loadAssemblyGraphs(studyAdapter);
        }
    }

    // NMR runs (legacy)
    NSMutableDictionary *nmrRuns = [NSMutableDictionary dictionary];
    if ([study hasChildNamed:@"nmr_runs"]) {
        TTIOHDF5Group *ng = [study openGroupNamed:@"nmr_runs" error:error];
        NSString *namesStr = [ng stringAttributeNamed:@"_run_names" error:error];
        for (NSString *rname in [namesStr componentsSeparatedByString:@","]) {
            if (rname.length == 0) continue;
            TTIOHDF5Group *runG = [ng openGroupNamed:rname error:error];
            BOOL exists = NO;
            NSUInteger n = (NSUInteger)[runG integerAttributeNamed:@"count"
                                                            exists:&exists error:error];
            NSMutableArray *spectra = [NSMutableArray arrayWithCapacity:n];
            for (NSUInteger i = 0; i < n; i++) {
                NSString *sname = [NSString stringWithFormat:@"spec_%06lu", (unsigned long)i];
                TTIONMRSpectrum *s = [TTIONMRSpectrum readFromGroup:runG name:sname error:error];
                if (!s) return nil;
                [spectra addObject:s];
            }
            nmrRuns[rname] = spectra;
        }
    }

    // Identifications, quantifications, provenance: compound if present
    // (v0.2 feature flags), JSON fallback otherwise .
    NSArray *idents = @[];
    NSArray *quants = @[];
    NSArray *prov   = @[];

    if (!isV1 &&
        [TTIOFeatureFlags root:root supportsFeature:[TTIOFeatureFlags featureCompoundIdentifications]] &&
        [study hasChildNamed:@"identifications"]) {
        idents = [TTIOCompoundIO readIdentificationsFromGroup:study
                                                 datasetNamed:@"identifications"
                                                        error:NULL] ?: @[];
    } else if ([study hasAttributeNamed:@"identifications_json"]) {
        idents = decodePlistArray([study stringAttributeNamed:@"identifications_json" error:NULL],
                                  [TTIOIdentification class], NULL) ?: @[];
    }

    if (!isV1 &&
        [TTIOFeatureFlags root:root supportsFeature:[TTIOFeatureFlags featureCompoundQuantifications]] &&
        [study hasChildNamed:@"quantifications"]) {
        quants = [TTIOCompoundIO readQuantificationsFromGroup:study
                                                 datasetNamed:@"quantifications"
                                                        error:NULL] ?: @[];
    } else if ([study hasAttributeNamed:@"quantifications_json"]) {
        quants = decodePlistArray([study stringAttributeNamed:@"quantifications_json" error:NULL],
                                  [TTIOQuantification class], NULL) ?: @[];
    }

    if (!isV1 &&
        [TTIOFeatureFlags root:root supportsFeature:[TTIOFeatureFlags featureCompoundProvenance]] &&
        [study hasChildNamed:@"provenance"]) {
        prov = [TTIOCompoundIO readProvenanceFromGroup:study
                                          datasetNamed:@"provenance"
                                                 error:NULL] ?: @[];
    } else if ([study hasAttributeNamed:@"provenance_json"]) {
        prov = decodePlistArray([study stringAttributeNamed:@"provenance_json" error:NULL],
                                [TTIOProvenanceRecord class], NULL) ?: @[];
    }

    TTIOTransitionList *trans = nil;
    if ([study hasAttributeNamed:@"transitions_json"]) {
        NSString *tjson = [study stringAttributeNamed:@"transitions_json" error:error];
        NSData *tdata = [tjson dataUsingEncoding:NSUTF8StringEncoding];
        NSDictionary *plist = [NSJSONSerialization JSONObjectWithData:tdata options:0 error:error];
        if (plist) trans = [TTIOTransitionList fromPlist:plist];
    }

    TTIOSpectralDataset *ds = [[self alloc] initWithTitle:title
                                        isaInvestigationId:isaId
                                                    msRuns:msRuns
                                                   nmrRuns:nmrRuns
                                           identifications:idents
                                           quantifications:quants
                                         provenanceRecords:prov
                                               transitions:trans];
    ds->_file        = f;
    ds->_provider    = p;
    ds->_filePath    = [path copy];
    ds->_genomicRuns = [genomicRuns copy];
    ds->_assemblyGraphs = [assemblyGraphs copy];

    // /study/references/<uri>/ — embedded references read-back
    // (1.1.0 — Phase 0 tio-browser). Empty dict when absent.
    NSMutableDictionary<NSString *, TTIOReferenceImport *> *refs =
        [NSMutableDictionary dictionary];
    if ([study hasChildNamed:@"references"]) {
        TTIOHDF5Group *refsG = [study openGroupNamed:@"references" error:NULL];
        if (refsG != nil) {
            for (NSString *uri in [refsG childNames]) {
                TTIOHDF5Group *oneRefG =
                    [refsG openGroupNamed:uri error:NULL];
                if (oneRefG == nil) continue;
                id<TTIOStorageGroup> oneRefAdapter =
                    (id<TTIOStorageGroup>)_TTIO_MakeHDF5GroupAdapter(oneRefG);
                TTIOReferenceImport *r =
                    [TTIOReferenceImport readFromGroup:oneRefAdapter];
                if (r != nil) refs[uri] = r;
            }
        }
    }
    ds->_references = [refs copy];

    if ([root hasAttributeNamed:@"access_policy_json"]) {
        ds->_accessPolicy = decodeAccessPolicy(
            [root stringAttributeNamed:@"access_policy_json" error:NULL]);
    }

    // Surface the root `encrypted` attribute (written by
    // -markRootEncryptedWithError:) so -isEncrypted / -encryptedAlgorithm
    // round-trip across close/reopen. Absent → empty string.
    if ([root hasAttributeNamed:@"encrypted"]) {
        NSString *alg = [root stringAttributeNamed:@"encrypted" error:NULL];
        ds->_encryptedAlgorithm = [(alg ?: @"") copy];
    } else {
        ds->_encryptedAlgorithm = @"";
    }

    return ds;
}

- (BOOL)closeFile
{
    // Cascade: runs hold open HDF5 group/dataset handles that would
    // otherwise keep the file alive even after [_file close].
    for (NSString *runName in _msRuns) {
        [[_msRuns objectForKey:runName] releaseHDF5Handles];
    }
    if (_file) {
        BOOL ok = [_file close];
        _file = nil;
        [_provider close];
        _provider = nil;
        return ok;
    }
    return YES;
}

- (NSArray<TTIOProvenanceRecord *> *)provenanceRecordsForInputRef:(NSString *)ref
{
    NSMutableArray *out = [NSMutableArray array];
    for (TTIOProvenanceRecord *r in _provenanceRecords) {
        if ([r containsInputRef:ref]) [out addObject:r];
    }
    return out;
}

#pragma mark - Phase 1 / Phase 2: modality-agnostic run accessors

- (NSDictionary<NSString *, id<TTIORun>> *)runs
{
    NSMutableDictionary *merged = [NSMutableDictionary dictionary];
    // MS first; preserve existing iteration order for caller stability.
    for (NSString *k in _msRuns) {
        merged[k] = _msRuns[k];
    }
    // Genomic runs second; do not overwrite an existing key (parity
    // with Python's dict.setdefault path).
    for (NSString *k in _genomicRuns) {
        if (!merged[k]) merged[k] = _genomicRuns[k];
    }
    return [merged copy];
}

- (NSDictionary<NSString *, id<TTIORun>> *)allRunsUnified
{
    return [self runs];
}

- (NSDictionary<NSString *, id<TTIORun>> *)runsForSample:(NSString *)sampleURI
{
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    if (sampleURI.length == 0) return @{};
    NSDictionary<NSString *, id<TTIORun>> *all = [self runs];
    for (NSString *name in all) {
        id<TTIORun> run = all[name];
        NSArray<TTIOProvenanceRecord *> *chain = nil;
        NS_DURING
            chain = [run provenanceChain];
        NS_HANDLER
            chain = nil;
        NS_ENDHANDLER
        for (TTIOProvenanceRecord *prov in chain) {
            if ([prov.inputRefs containsObject:sampleURI]) {
                out[name] = run;
                break;
            }
        }
    }
    return [out copy];
}

- (NSDictionary<NSString *, id<TTIORun>> *)runsOfModality:(Class)runClass
{
    NSMutableDictionary *out = [NSMutableDictionary dictionary];
    if (runClass == Nil) return @{};
    NSDictionary<NSString *, id<TTIORun>> *all = [self runs];
    for (NSString *name in all) {
        id<TTIORun> run = all[name];
        if ([(NSObject *)run isKindOfClass:runClass]) {
            out[name] = run;
        }
    }
    return [out copy];
}

#pragma mark - TTIOEncryptable

- (BOOL)encryptWithKey:(NSData *)key
                 level:(TTIOEncryptionLevel)level
                 error:(NSError **)error
{
    if (!_filePath) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"TTIOSpectralDataset: cannot encrypt before the dataset has been persisted");
        return NO;
    }

    // Release our handle so the encryption manager can reopen rw.
    [self closeFile];

    // 1. Encrypt each MS run's intensity channel by delegating to the
    //    run's own protocol method (persistence context was set during
    //    the load or by the caller after initial write).
    for (NSString *runName in _msRuns) {
        TTIOAcquisitionRun *run = _msRuns[runName];
        // Use the full HDF5 path since runs live under /study/ms_runs/
        // when persisted by TTIOSpectralDataset. H5Gopen2 accepts slash-
        // separated paths, so the encryption manager can locate the run.
        NSString *fullPath = [NSString stringWithFormat:@"/study/ms_runs/%@", runName];
        [run setPersistenceFilePath:_filePath runName:fullPath];
        if (![run encryptWithKey:key level:level error:error]) return NO;
    }

    // 2. Seal compound identifications + quantifications into encrypted
    //    byte blobs under /study/, dropping the plaintext compound
    //    datasets.
    if (![self sealCompoundDatasetsWithKey:key error:error]) return NO;

    // 3. Mark the root + persist access policy.
    if (![self markRootEncryptedWithError:error]) return NO;

    // 4. Mirror the on-disk attr in memory so -isEncrypted /
    //    -encryptedAlgorithm return the new state without a reopen.
    _encryptedAlgorithm = @"aes-256-gcm";

    return YES;
}

/**
 * Decrypt every MS run's intensity channel into an in-memory overlay.
 *
 * <b>Read-only / asymmetric</b>: the on-disk file is NOT modified,
 * the root <code>@encrypted</code> attribute is left in place, and
 * <code>-isEncrypted</code> continues to return YES on this instance
 * and on any reopen. This is the asymmetric counterpart to
 * <code>-encryptWithKey:level:error:</code> (which IS persistent +
 * flag-flipping) by design: in-memory rehydration lets a process
 * read encrypted data without rewriting the file.
 *
 * To fully reverse encryption on disk and clear the
 * <code>@encrypted</code> attribute, use the class method
 * <code>+decryptInPlaceAtPath:withKey:error:</code> after closing
 * any open instance.
 *
 * <p>Cross-language equivalents: Java
 * <code>SpectralDataset.decryptWithKey(byte[])</code> (same), Python
 * <code>SpectralDataset.decrypt_with_key</code> (same).</p>
 */
- (BOOL)decryptWithKey:(NSData *)key error:(NSError **)error
{
    if (!_filePath) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"TTIOSpectralDataset: no file path to decrypt against");
        return NO;
    }
    [self closeFile];

    for (NSString *runName in _msRuns) {
        TTIOAcquisitionRun *run = _msRuns[runName];
        NSString *fullPath = [NSString stringWithFormat:@"/study/ms_runs/%@", runName];
        [run setPersistenceFilePath:_filePath runName:fullPath];
        if (![run decryptWithKey:key error:error]) return NO;
    }

    if (![self unsealCompoundDatasetsWithKey:key error:error]) return NO;

    // M5-handoff: reopen the file read-only and reattach each run's
    // signal-channel handles so -spectrumAtIndex: can serve both the
    // decrypted intensity channel (from the run's in-memory cache) and
    // any unencrypted channels (mz, chemical_shift) from disk. The
    // on-disk file still carries the `encrypted` attribute and the
    // ciphertext datasets — decryption does not modify the file.
    return [self reopenAfterDecryptWithError:error];
}

- (BOOL)reopenAfterDecryptWithError:(NSError **)error
{
    TTIOHDF5Provider *p = [[TTIOHDF5Provider alloc] init];
    if (![p openURL:_filePath mode:TTIOStorageOpenModeRead error:error]) return NO;
    TTIOHDF5File *f = (TTIOHDF5File *)[p nativeHandle];
    if (!f) return NO;
    TTIOHDF5Group *root = [f rootGroup];
    if (!root) { [p close]; return NO; }
    TTIOHDF5Group *study = [root openGroupNamed:@"study" error:error];
    if (!study) { [p close]; return NO; }
    TTIOHDF5Group *msRunsG = nil;
    if ([study hasChildNamed:@"ms_runs"]) {
        msRunsG = [study openGroupNamed:@"ms_runs" error:error];
        if (!msRunsG) { [p close]; return NO; }
    }
    for (NSString *runName in _msRuns) {
        if (!msRunsG) break;
        TTIOHDF5Group *runG = [msRunsG openGroupNamed:runName error:NULL];
        if (!runG) continue;
        TTIOHDF5Group *channels = [runG openGroupNamed:@"signal_channels" error:NULL];
        if (!channels) continue;
        TTIOAcquisitionRun *run = _msRuns[runName];
        (void)[run reattachSignalHandlesFromGroup:channels error:NULL];
    }
    _file     = f;
    _provider = p;
    return YES;
}

- (TTIOAccessPolicy *)accessPolicy { return _accessPolicy; }
- (void)setAccessPolicy:(TTIOAccessPolicy *)policy { _accessPolicy = policy; }

#pragma mark - Compound dataset sealing (encryption of /study compound datasets)

- (BOOL)sealCompoundDatasetsWithKey:(NSData *)key error:(NSError **)error
{
    TTIOHDF5File *f = [TTIOHDF5File openAtPath:_filePath error:error];
    if (!f) return NO;
    TTIOHDF5Group *root = [f rootGroup];
    TTIOHDF5Group *study = [root openGroupNamed:@"study" error:error];
    if (!study) { [f close]; return NO; }

    if ([study hasChildNamed:@"identifications"]) {
        NSArray *idents = [TTIOCompoundIO readIdentificationsFromGroup:study
                                                          datasetNamed:@"identifications"
                                                                 error:error];
        if (!idents) { [f close]; return NO; }
        NSMutableArray *plists = [NSMutableArray array];
        for (TTIOIdentification *i in idents) [plists addObject:[i asPlist]];
        NSData *json = [NSJSONSerialization dataWithJSONObject:plists options:0 error:error];
        if (!json) { [f close]; return NO; }
        H5Ldelete(study.groupId, "identifications", H5P_DEFAULT);
        // also strip the JSON attribute mirror so sealed files are
        // not readable without decryption.
        if ([study hasAttributeNamed:@"identifications_json"])
            H5Adelete(study.groupId, "identifications_json");
        if (![self writeSealedBlob:json name:@"identifications_sealed"
                           inGroup:study key:key error:error]) { [f close]; return NO; }
    }

    if ([study hasChildNamed:@"quantifications"]) {
        NSArray *quants = [TTIOCompoundIO readQuantificationsFromGroup:study
                                                          datasetNamed:@"quantifications"
                                                                 error:error];
        if (!quants) { [f close]; return NO; }
        NSMutableArray *plists = [NSMutableArray array];
        for (TTIOQuantification *q in quants) [plists addObject:[q asPlist]];
        NSData *json = [NSJSONSerialization dataWithJSONObject:plists options:0 error:error];
        if (!json) { [f close]; return NO; }
        H5Ldelete(study.groupId, "quantifications", H5P_DEFAULT);
        if ([study hasAttributeNamed:@"quantifications_json"])
            H5Adelete(study.groupId, "quantifications_json");
        if (![self writeSealedBlob:json name:@"quantifications_sealed"
                           inGroup:study key:key error:error]) { [f close]; return NO; }
    }

    return [f close];
}

- (BOOL)unsealCompoundDatasetsWithKey:(NSData *)key error:(NSError **)error
{
    TTIOHDF5File *f = [TTIOHDF5File openAtPath:_filePath error:error];
    if (!f) return NO;
    TTIOHDF5Group *root = [f rootGroup];
    TTIOHDF5Group *study = [root openGroupNamed:@"study" error:error];
    if (!study) { [f close]; return NO; }

    if ([study hasChildNamed:@"identifications_sealed"]) {
        NSData *json = [self readSealedBlob:@"identifications_sealed"
                                    inGroup:study key:key error:error];
        if (!json) { [f close]; return NO; }
        NSArray *plists = [NSJSONSerialization JSONObjectWithData:json options:0 error:error];
        NSMutableArray *idents = [NSMutableArray array];
        for (NSDictionary *p in plists) [idents addObject:[TTIOIdentification fromPlist:p]];
        H5Ldelete(study.groupId, "identifications_sealed", H5P_DEFAULT);
        H5Ldelete(study.groupId, "identifications_sealed_iv", H5P_DEFAULT);
        H5Ldelete(study.groupId, "identifications_sealed_tag", H5P_DEFAULT);
        H5Ldelete(study.groupId, "identifications_sealed_bytes", H5P_DEFAULT);
        if (![TTIOCompoundIO writeIdentifications:idents
                                        intoGroup:study
                                     datasetNamed:@"identifications"
                                            error:error]) { [f close]; return NO; }
        _identifications = [idents copy];
    }

    if ([study hasChildNamed:@"quantifications_sealed"]) {
        NSData *json = [self readSealedBlob:@"quantifications_sealed"
                                    inGroup:study key:key error:error];
        if (!json) { [f close]; return NO; }
        NSArray *plists = [NSJSONSerialization JSONObjectWithData:json options:0 error:error];
        NSMutableArray *quants = [NSMutableArray array];
        for (NSDictionary *p in plists) [quants addObject:[TTIOQuantification fromPlist:p]];
        H5Ldelete(study.groupId, "quantifications_sealed", H5P_DEFAULT);
        H5Ldelete(study.groupId, "quantifications_sealed_iv", H5P_DEFAULT);
        H5Ldelete(study.groupId, "quantifications_sealed_tag", H5P_DEFAULT);
        H5Ldelete(study.groupId, "quantifications_sealed_bytes", H5P_DEFAULT);
        if (![TTIOCompoundIO writeQuantifications:quants
                                        intoGroup:study
                                     datasetNamed:@"quantifications"
                                            error:error]) { [f close]; return NO; }
        _quantifications = [quants copy];
    }

    return [f close];
}

- (BOOL)writeSealedBlob:(NSData *)plaintext
                   name:(NSString *)name
                inGroup:(TTIOHDF5Group *)group
                    key:(NSData *)key
                  error:(NSError **)error
{
    NSData *iv = nil, *tag = nil;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    NSData *cipher = [TTIOEncryptionManager encryptData:plaintext
                                                withKey:key
                                                     iv:&iv
                                                authTag:&tag
                                                  error:error];
#pragma clang diagnostic pop
    if (!cipher) return NO;

    NSMutableData *padded = [NSMutableData dataWithData:cipher];
    while (padded.length % 4 != 0) {
        uint8_t zero = 0;
        [padded appendBytes:&zero length:1];
    }
    TTIOHDF5Dataset *ds = [group createDatasetNamed:name
                                           precision:TTIOPrecisionInt32
                                              length:padded.length / 4
                                           chunkSize:0
                                    compressionLevel:0
                                               error:error];
    if (!ds) return NO;
    if (![ds writeData:padded error:error]) return NO;

    NSString *ivName  = [name stringByAppendingString:@"_iv"];
    NSString *tagName = [name stringByAppendingString:@"_tag"];
    NSString *lenName = [name stringByAppendingString:@"_bytes"];

    TTIOHDF5Dataset *ivDs = [group createDatasetNamed:ivName
                                             precision:TTIOPrecisionInt32
                                                length:3   // 12 bytes
                                             chunkSize:0
                                      compressionLevel:0
                                                 error:error];
    if (![ivDs writeData:iv error:error]) return NO;

    TTIOHDF5Dataset *tagDs = [group createDatasetNamed:tagName
                                              precision:TTIOPrecisionInt32
                                                 length:4   // 16 bytes
                                              chunkSize:0
                                       compressionLevel:0
                                                  error:error];
    if (![tagDs writeData:tag error:error]) return NO;

    // Store original cipher length (before padding) as 1-element dataset
    uint32_t lenBytes = (uint32_t)cipher.length;
    TTIOHDF5Dataset *lenDs = [group createDatasetNamed:lenName
                                              precision:TTIOPrecisionUInt32
                                                 length:1
                                              chunkSize:0
                                       compressionLevel:0
                                                  error:error];
    return [lenDs writeData:[NSData dataWithBytes:&lenBytes length:sizeof(lenBytes)] error:error];
}

- (NSData *)readSealedBlob:(NSString *)name
                   inGroup:(TTIOHDF5Group *)group
                       key:(NSData *)key
                     error:(NSError **)error
{
    TTIOHDF5Dataset *ds = [group openDatasetNamed:name error:error];
    if (!ds) return nil;
    NSData *padded = [ds readDataWithError:error];
    if (!padded) return nil;

    NSString *lenName = [name stringByAppendingString:@"_bytes"];
    TTIOHDF5Dataset *lenDs = [group openDatasetNamed:lenName error:error];
    NSData *lenData = [lenDs readDataWithError:error];
    uint32_t cipherLen = ((const uint32_t *)lenData.bytes)[0];

    NSData *cipher = [padded subdataWithRange:NSMakeRange(0, cipherLen)];

    NSString *ivName  = [name stringByAppendingString:@"_iv"];
    NSString *tagName = [name stringByAppendingString:@"_tag"];
    NSData *iv  = [[group openDatasetNamed:ivName  error:error] readDataWithError:error];
    NSData *tag = [[group openDatasetNamed:tagName error:error] readDataWithError:error];

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    return [TTIOEncryptionManager decryptData:cipher
                                      withKey:key
                                           iv:iv
                                      authTag:tag
                                        error:error];
#pragma clang diagnostic pop
}

- (BOOL)markRootEncryptedWithError:(NSError **)error
{
    TTIOHDF5File *f = [TTIOHDF5File openAtPath:_filePath error:error];
    if (!f) return NO;
    TTIOHDF5Group *root = [f rootGroup];
    if (![root setStringAttribute:@"encrypted"
                             value:@"aes-256-gcm"
                             error:error]) { [f close]; return NO; }

    NSString *ap = encodeAccessPolicy(_accessPolicy);
    if (ap) {
        if (![root setStringAttribute:@"access_policy_json"
                                value:ap
                                error:error]) { [f close]; return NO; }
    }
    return [f close];
}

#pragma mark - v1.1.1: persist-to-disk decrypt

+ (BOOL)decryptInPlaceAtPath:(NSString *)path
                     withKey:(NSData *)key
                       error:(NSError **)error
{
    if (key.length != 32) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"AES-256-GCM requires a 32-byte key, got %lu",
            (unsigned long)key.length);
        return NO;
    }

    // 1. Enumerate MS run names while the file is closed to readers.
    NSMutableArray<NSString *> *runNames = [NSMutableArray array];
    {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:error];
        if (!f) return NO;
        TTIOHDF5Group *root = [f rootGroup];
        if ([root hasChildNamed:@"study"]) {
            TTIOHDF5Group *study = [root openGroupNamed:@"study" error:error];
            if (!study) { [f close]; return NO; }
            if ([study hasChildNamed:@"ms_runs"]) {
                TTIOHDF5Group *msRunsG =
                    [study openGroupNamed:@"ms_runs" error:error];
                if (!msRunsG) { [f close]; return NO; }
                for (NSString *name in [msRunsG childNames]) {
                    [runNames addObject:name];
                }
            }
        }
        if (![f close]) return NO;
    }

    // 2. Decrypt each run's intensity channel in place. The encryption
    //    manager opens/closes the file for each call, mirroring the
    //    encrypt side's per-run lifecycle.
    for (NSString *name in runNames) {
        NSString *fullPath =
            [NSString stringWithFormat:@"/study/ms_runs/%@", name];
        if (![TTIOEncryptionManager
                decryptIntensityChannelInRunInPlace:fullPath
                                         atFilePath:path
                                            withKey:key
                                              error:error]) {
            return NO;
        }
    }

    // 3. Clear the root @encrypted attribute so a reopen sees the
    //    file as unprotected.
    TTIOHDF5File *fw = [TTIOHDF5File openAtPath:path error:error];
    if (!fw) return NO;
    TTIOHDF5Group *root = [fw rootGroup];
    if ([root hasAttributeNamed:@"encrypted"]) {
        if (![root deleteAttributeNamed:@"encrypted" error:error]) {
            [fw close];
            return NO;
        }
    }
    return [fw close];
}

@end

#import "../Image/TTIOMSImage.h"
#import "../Image/TTIORamanImage.h"
#import "../Image/TTIOIRImage.h"
#import "TTIOSubject.h"
#import "TTIOSample.h"
#import "../HDF5/TTIOHDF5File.h"
#import "../HDF5/TTIOHDF5Group.h"
#import <objc/runtime.h>

// ─── Stage 6 (transport-spec v0.11, Deferral 2): Subjects + Samples ──
//
// /study/subjects/<external_id>/ and /study/samples/<sample_id>/ per-row
// HDF5 group readers, mirroring Java SpectralDataset.readSubjects /
// readSamples (commit dd39f4e6) and Python SpectralDataset.subjects /
// samples (commit 721ad21c). Lazy + cached via objc_setAssociatedObject
// in the same pattern as -msImage / -ramanImage / -irImage above.

static NSDictionary<NSString *, NSString *> *
TTIOParseAttributesJsonString(NSString *blob)
{
    if (blob == nil || blob.length == 0 || [blob isEqualToString:@"{}"]) {
        return @{};
    }
    NSData *data = [blob dataUsingEncoding:NSUTF8StringEncoding];
    if (data == nil) return @{};
    NSError *err = nil;
    id parsed = [NSJSONSerialization JSONObjectWithData:data
                                                options:0
                                                  error:&err];
    if (![parsed isKindOfClass:[NSDictionary class]]) return @{};
    NSMutableDictionary<NSString *, NSString *> *out =
        [NSMutableDictionary dictionaryWithCapacity:[parsed count]];
    for (id key in (NSDictionary *)parsed) {
        if (![key isKindOfClass:[NSString class]]) continue;
        id v = ((NSDictionary *)parsed)[key];
        if ([v isKindOfClass:[NSString class]]) {
            out[(NSString *)key] = (NSString *)v;
        } else if (v != [NSNull null]) {
            out[(NSString *)key] = [v description];
        }
    }
    return [out copy];
}

static NSArray<TTIOSubject *> *
TTIOReadSubjectsFromFile(NSString *path)
{
    if (path == nil) return @[];
    NSError *err = nil;
    // Read-only — the lazy accessor coexists with the
    // readFromFilePath: provider that also holds the file RO. HDF5
    // allows multiple RO handles per file but rejects RDWR while RO
    // is held.
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:&err];
    if (f == nil) return @[];
    TTIOHDF5Group *root = [f rootGroup];
    TTIOHDF5Group *study = nil;
    if ([root hasChildNamed:@"study"]) {
        study = [root openGroupNamed:@"study" error:NULL];
    }
    if (study == nil || ![study hasChildNamed:@"subjects"]) {
        [f close];
        return @[];
    }
    TTIOHDF5Group *subjects = [study openGroupNamed:@"subjects" error:NULL];
    if (subjects == nil) { [f close]; return @[]; }
    NSMutableArray<TTIOSubject *> *out = [NSMutableArray array];
    for (NSString *name in [subjects childNames]) {
        TTIOHDF5Group *row = [subjects openGroupNamed:name error:NULL];
        if (row == nil) continue;
        NSString *externalId = name;
        if ([row hasAttributeNamed:@"external_id"]) {
            NSString *v = [row stringAttributeNamed:@"external_id" error:NULL];
            if (v != nil) externalId = v;
        }
        NSString *project = @"";
        if ([row hasAttributeNamed:@"project"]) {
            NSString *v = [row stringAttributeNamed:@"project" error:NULL];
            if (v != nil) project = v;
        }
        NSString *sex = @"";
        if ([row hasAttributeNamed:@"sex"]) {
            NSString *v = [row stringAttributeNamed:@"sex" error:NULL];
            if (v != nil) sex = v;
        }
        int64_t birthYear = 0;
        if ([row hasAttributeNamed:@"birth_year"]) {
            BOOL exists = NO;
            int64_t v = [row integerAttributeNamed:@"birth_year"
                                            exists:&exists error:NULL];
            if (exists) birthYear = v;
        }
        NSDictionary<NSString *, NSString *> *attrs = @{};
        if ([row hasAttributeNamed:@"attributes_json"]) {
            NSString *v = [row stringAttributeNamed:@"attributes_json" error:NULL];
            attrs = TTIOParseAttributesJsonString(v);
        }
        @try {
            TTIOSubject *s = [[TTIOSubject alloc]
                initWithExternalId:externalId
                            project:project
                                sex:sex
                          birthYear:birthYear
                         attributes:attrs];
            [out addObject:s];
        } @catch (NSException *exc) {
            // Skip rows that fail validation (e.g. empty externalId on
            // a corrupted file) — readers MAY tolerate these per design
            // spec §4.4 forward-compat note.
            NSLog(@"TTIOSpectralDataset.subjects: skipping invalid row "
                  @"'%@': %@", name, exc.reason);
        }
    }
    [f close];
    return [out copy];
}

static NSArray<TTIOSample *> *
TTIOReadSamplesFromFile(NSString *path)
{
    if (path == nil) return @[];
    NSError *err = nil;
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:&err];
    if (f == nil) return @[];
    TTIOHDF5Group *root = [f rootGroup];
    TTIOHDF5Group *study = nil;
    if ([root hasChildNamed:@"study"]) {
        study = [root openGroupNamed:@"study" error:NULL];
    }
    if (study == nil || ![study hasChildNamed:@"samples"]) {
        [f close];
        return @[];
    }
    TTIOHDF5Group *samples = [study openGroupNamed:@"samples" error:NULL];
    if (samples == nil) { [f close]; return @[]; }
    NSMutableArray<TTIOSample *> *out = [NSMutableArray array];
    for (NSString *name in [samples childNames]) {
        TTIOHDF5Group *row = [samples openGroupNamed:name error:NULL];
        if (row == nil) continue;
        NSString *sampleId = name;
        if ([row hasAttributeNamed:@"sample_id"]) {
            NSString *v = [row stringAttributeNamed:@"sample_id" error:NULL];
            if (v != nil) sampleId = v;
        }
        NSString *subjectExternalId = @"";
        if ([row hasAttributeNamed:@"subject_external_id"]) {
            NSString *v = [row stringAttributeNamed:@"subject_external_id" error:NULL];
            if (v != nil) subjectExternalId = v;
        }
        NSString *sampleKind = @"";
        if ([row hasAttributeNamed:@"sample_kind"]) {
            NSString *v = [row stringAttributeNamed:@"sample_kind" error:NULL];
            if (v != nil) sampleKind = v;
        }
        int64_t collectedAt = 0;
        if ([row hasAttributeNamed:@"collected_at"]) {
            BOOL exists = NO;
            int64_t v = [row integerAttributeNamed:@"collected_at"
                                            exists:&exists error:NULL];
            if (exists) collectedAt = v;
        }
        NSDictionary<NSString *, NSString *> *attrs = @{};
        if ([row hasAttributeNamed:@"attributes_json"]) {
            NSString *v = [row stringAttributeNamed:@"attributes_json" error:NULL];
            attrs = TTIOParseAttributesJsonString(v);
        }
        @try {
            TTIOSample *s = [[TTIOSample alloc]
                initWithSampleId:sampleId
               subjectExternalId:subjectExternalId
                      sampleKind:sampleKind
                     collectedAt:collectedAt
                      attributes:attrs];
            [out addObject:s];
        } @catch (NSException *exc) {
            NSLog(@"TTIOSpectralDataset.samples: skipping invalid row "
                  @"'%@': %@", name, exc.reason);
        }
    }
    [f close];
    return [out copy];
}

@implementation TTIOSpectralDataset (SubjectsSamples)

- (NSArray<TTIOSubject *> *)subjects
{
    // Lazy + cached, mirroring -msImage / -ramanImage / -irImage above
    // (Stage 5.2 commit 50ef8bc3).
    static const void * const kSubjectsCacheKey = &kSubjectsCacheKey;
    NSArray<TTIOSubject *> *cached =
        objc_getAssociatedObject(self, kSubjectsCacheKey);
    if (cached != nil) return cached;
    NSArray<TTIOSubject *> *out = TTIOReadSubjectsFromFile([self filePath]);
    if (out == nil) out = @[];
    objc_setAssociatedObject(self, kSubjectsCacheKey, out,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC);
    return out;
}

- (NSArray<TTIOSample *> *)samples
{
    static const void * const kSamplesCacheKey = &kSamplesCacheKey;
    NSArray<TTIOSample *> *cached =
        objc_getAssociatedObject(self, kSamplesCacheKey);
    if (cached != nil) return cached;
    NSArray<TTIOSample *> *out = TTIOReadSamplesFromFile([self filePath]);
    if (out == nil) out = @[];
    objc_setAssociatedObject(self, kSamplesCacheKey, out,
                              OBJC_ASSOCIATION_RETAIN_NONATOMIC);
    return out;
}

+ (void)validateSubjects:(NSArray<TTIOSubject *> *)subjects
                  samples:(NSArray<TTIOSample *> *)samples
{
    // Design spec §4.4: duplicate IDs raise; soft-FK mismatch is a
    // WARNING. Mirrors Java SpectralDataset.validateSubjectsAndSamples
    // (commit dd39f4e6) and Python SpectralDataset._validate_stage6
    // (commit 721ad21c).
    NSMutableSet<NSString *> *seenSubjects = [NSMutableSet set];
    for (TTIOSubject *s in subjects) {
        if ([seenSubjects containsObject:s.externalId]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"duplicate Subject.externalId: %@",
                                s.externalId];
        }
        [seenSubjects addObject:s.externalId];
    }
    NSMutableSet<NSString *> *seenSamples = [NSMutableSet set];
    for (TTIOSample *s in samples) {
        if ([seenSamples containsObject:s.sampleId]) {
            [NSException raise:NSInvalidArgumentException
                        format:@"duplicate Sample.sampleId: %@", s.sampleId];
        }
        [seenSamples addObject:s.sampleId];
    }
    for (TTIOSample *s in samples) {
        NSString *fk = s.subjectExternalId;
        if (fk == nil || fk.length == 0) continue;
        if (![seenSubjects containsObject:fk]) {
            NSLog(@"WARNING: Sample '%@' references unknown "
                  @"Subject.externalId '%@' — soft-FK mismatch, "
                  @"writing anyway (spec §4.4).",
                  s.sampleId, fk);
        }
    }
}

@end

// Private lazy/cached per-modality materialisers backing -imageForKind:.
// Each reads + caches its image from the dataset file on first access
// (associated-object cache, preserved verbatim from the former public
// -msImage / -ramanImage / -irImage property getters).
@interface TTIOSpectralDataset (ImagePrivate)
- (nullable TTIOMSImage *)_lazyMSImage;
- (nullable TTIORamanImage *)_lazyRamanImage;
- (nullable TTIOIRImage *)_lazyIRImage;
@end

@implementation TTIOSpectralDataset (Image)

- (nullable TTIOImage *)imageForKind:(TTIOImageKind)kind
{
    switch (kind) {
        case TTIOImageKindMS:    return [self _lazyMSImage];
        case TTIOImageKindRaman: return [self _lazyRamanImage];
        case TTIOImageKindIR:    return [self _lazyIRImage];
    }
    return nil;
}

- (NSArray<TTIOImage *> *)images
{
    NSMutableArray<TTIOImage *> *out = [NSMutableArray array];
    for (TTIOImageKind k = TTIOImageKindMS; k <= TTIOImageKindIR; k++) {
        TTIOImage *img = [self imageForKind:k];
        if (img) [out addObject:img];
    }
    return out;
}

- (TTIOMSImage *)_lazyMSImage
{
    // TTIOMSImage no longer inherits from TTIOSpectralDataset; always
    // materialise from file path via the image class factory.
    static const void * const kMsImageCacheKey = &kMsImageCacheKey;
    TTIOMSImage *cached = objc_getAssociatedObject(self, kMsImageCacheKey);
    if (cached != nil) {
        return cached;
    }
    NSString *path = [self filePath];
    if (path == nil) return nil;
    NSError *err = nil;
    TTIOMSImage *img = [TTIOMSImage readFromFilePath:path error:&err];
    if (img != nil) {
        objc_setAssociatedObject(self, kMsImageCacheKey, img,
                                  OBJC_ASSOCIATION_RETAIN_NONATOMIC);
    }
    return img;
}

- (TTIORamanImage *)_lazyRamanImage
{
    // v0.11 Stage 5.2: lazy /study/raman_image_cube accessor mirroring
    // -msImage. Java parity: SpectralDataset.ramanImage(). Python parity:
    // SpectralDataset.raman_image (lazy property).
    static const void * const kRamanImageCacheKey = &kRamanImageCacheKey;
    TTIORamanImage *cached = objc_getAssociatedObject(self, kRamanImageCacheKey);
    if (cached != nil) return cached;
    NSString *path = [self filePath];
    if (path == nil) return nil;
    NSError *err = nil;
    TTIORamanImage *img = [TTIORamanImage readFromFilePath:path error:&err];
    if (img != nil) {
        objc_setAssociatedObject(self, kRamanImageCacheKey, img,
                                  OBJC_ASSOCIATION_RETAIN_NONATOMIC);
    }
    return img;
}

- (TTIOIRImage *)_lazyIRImage
{
    // v0.11 Stage 5.2: lazy /study/ir_image_cube accessor mirroring
    // -msImage / -ramanImage. Java parity: SpectralDataset.irImage()
    // (commit 97fb065e). Python parity: SpectralDataset.ir_image
    // (commit 8b57baa7).
    static const void * const kIRImageCacheKey = &kIRImageCacheKey;
    TTIOIRImage *cached = objc_getAssociatedObject(self, kIRImageCacheKey);
    if (cached != nil) return cached;
    NSString *path = [self filePath];
    if (path == nil) return nil;
    NSError *err = nil;
    TTIOIRImage *img = [TTIOIRImage readFromFilePath:path error:&err];
    if (img != nil) {
        objc_setAssociatedObject(self, kIRImageCacheKey, img,
                                  OBJC_ASSOCIATION_RETAIN_NONATOMIC);
    }
    return img;
}

@end
