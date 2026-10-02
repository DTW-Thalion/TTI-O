/*
 * TTIOAcquisitionRun.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOAcquisitionRun
 * Inherits From: NSObject
 * Conforms To:   TTIOIndexable, TTIOStreamable, TTIOProvenanceable,
 *                TTIOEncryptable, TTIORun
 * Declared In:   Run/TTIOAcquisitionRun.h
 *
 * Ordered run of spectra sharing instrument configuration and
 * acquisition mode. Provider-agnostic write/read; cooperates with
 * TTIOSpectralDataset for in-place encryption / decryption.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOAcquisitionRun.h"
#import "Core/TTIOThreads.h"
#import "TTIOInstrumentConfig.h"
#import "TTIOSpectrumIndex.h"
#import "Genomics/TTIOGenomicIndex.h"  // TTIOOffsetsFromLengths (v1.10 #10)
#import "Spectra/TTIOSpectrum.h"
#import "Spectra/TTIOMassSpectrum.h"
#import "Spectra/TTIONMRSpectrum.h"
#import "Spectra/TTIOIRSpectrum.h"
#import "Spectra/TTIORamanSpectrum.h"
#import "Spectra/TTIOUVVisSpectrum.h"
#import "Spectra/TTIOChromatogram.h"
#import "Codecs/TTIOFloatDeltaZstd.h"
#import "Run/TTIOSpectralBlockIndex.h"
#import "Run/TTIOSpectralUnitPlan.h"
#import "Core/TTIOSignalArray.h"
#import "ValueClasses/TTIOEnums.h"
#import "ValueClasses/TTIOEncodingSpec.h"
#import "ValueClasses/TTIOValueRange.h"
#import "ValueClasses/TTIOIsolationWindow.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Dataset/TTIOCompoundIO.h"
#import "Core/TTIONumpress.h"
#import "Protection/TTIOEncryptionManager.h"
#import "Protection/TTIOAccessPolicy.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "HDF5/TTIOHDF5Dataset.h"
#import "HDF5/TTIOHDF5Errors.h"
#import "HDF5/TTIOHDF5Types.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOHDF5Provider.h"
#import <pthread.h>

// Shared immutable flyweight for the standard float64 / zlib /
// little-endian channel encoding used by every spectrum read in
// -spectrumAtIndex:error:. The spec is value-identical on every call
// and TTIOEncodingSpec is immutable (all properties readonly, no
// mutators), so a single process-wide instance is safe to hand to
// every SignalArray — eliminating one allocation per channel per
// spectrum on the hot read path. Built once via pthread_once (the
// codebase's established once-init idiom on GNUstep/Linux).
static TTIOEncodingSpec *gStdChannelEncoding = nil;
static pthread_once_t   gStdChannelEncodingOnce = PTHREAD_ONCE_INIT;
static void _buildStdChannelEncoding(void)
{
    gStdChannelEncoding =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
}

/** One FLOAT_DELTA_ZSTD block decoding in flight. */
@interface TTIOInFlightFdzDecode : NSObject
@property (nonatomic, strong, nullable) NSData *values;
@property (nonatomic) BOOL done;
@end

@implementation TTIOInFlightFdzDecode
@end

@interface TTIOInFlightUnit : NSObject
@property (nonatomic, strong, nullable) NSError *error;
@property (nonatomic) BOOL done;
@end
@implementation TTIOInFlightUnit
@end

@interface TTIOAcquisitionRun ()
/* Defined below, used by the two view builders above it. */
- (instancetype)_initAsUnitViewOfRun:(TTIOAcquisitionRun *)parent
                                unit:(TTIOSpectralUnit)u
                             columns:(NSDictionary<NSString *, NSData *> *)cols;
@end

@implementation TTIOAcquisitionRun
{
    NSUInteger _iterThreads;
    // Phase 1: Run protocol name. Set by readFromGroup:name: /
    // readFromStorageGroup:name: at load time, by
    // setPersistenceFilePath:runName: post-load (kept in sync with
    // _persistenceRunName below), or remains @"" for in-memory runs
    // not yet persisted.
    NSString                    *_name;
    NSArray                     *_inMemorySpectra;       // nil when read-from-disk
    // storage-protocol iVars. Populated by
    // +readFromGroup:name:error: via TTIOHDF5Provider's adapter
    // factory; the hot spectrum-read path goes through the protocol
    // (readSliceAtOffset:count:error:) so a future non-HDF5 provider
    // can host a run without per-class migration.
    id<TTIOStorageGroup>         _storageSignalGroup;    // nil when in-memory
    // blocks/index of a codec-17 run. nil for every run written
    // before the group existed and for runs whose channels fell
    // out of step, which is not an error.
    TTIOSpectralBlockIndex      *_spectralBlockIndex;
    // Unit-view state. A view handed to an -iterBlocksFrom: visitor
    // holds its unit's decoded channel columns and builds a spectrum
    // only when one is asked for. Materialising a whole unit up front
    // cost 14 per cent against the ordered reader at one thread.
    NSDictionary<NSString *, NSData *> *_unitColumns;   // nil unless a view
    NSUInteger                   _unitFirstSpectrum;
    NSUInteger                   _unitCount;
    unsigned long long           _unitValueStart;
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *_storageDatasets;
    NSArray<NSString *>         *_channelNames;          // ordered list
    NSUInteger                   _streamPosition;

    NSMutableArray<TTIOProvenanceRecord *> *_provenance;
    TTIOAccessPolicy            *_accessPolicy;

    // eagerly decoded Numpress-delta channels, keyed by channel
    // name. When a channel is present here, spectrumAtIndex: slices
    // into this float64 buffer instead of reading the HDF5 dataset,
    // because Numpress decoding needs the running sum prefix.
    NSMutableDictionary<NSString *, NSData *> *_numpressChannels;

    // M5-handoff: in-memory plaintext channels populated by
    // -decryptWithKey:error:. Keyed by channel name. When present,
    // spectrumAtIndex: slices into this float64 buffer so spectra are
    // readable through the normal API after decrypt without modifying
    // the on-disk file (mirrors the Python rehydrate-in-memory
    // semantics in TTIOAcquisitionRun).
    NSMutableDictionary<NSString *, NSData *> *_decryptedChannels;

    // Perf (Fix #1): lazy per-channel full-column cache for read-only
    // disk-backed runs. spectrumAtIndex: previously issued one HDF5
    // hyperslab read per channel per spectrum (~200k round-trips for a
    // 100k-AU 2-channel run). On first access of a channel we read the
    // WHOLE float64 column once via the storage protocol's readAll: and
    // retain it here, then slice element-wise — byte-identical to the
    // prior per-slice read (readAll and readSliceAtOffset:count: return
    // the same packed-LE bytes; M43 cross-backend identity covers this).
    // Java loads all channels at open; this loads only accessed channels,
    // lazily. Invalidated by -releaseHDF5Handles / -reattachSignalHandles
    // so it can never serve stale bytes after the dataset handles change.
    NSMutableDictionary<NSString *, NSData *> *_cachedFullChannels;

    // FLOAT_DELTA_ZSTD channels are left encoded on disk at open; a
    // range read decodes the blocks it covers through the block table,
    // keeping the last block per channel.
    NSMutableSet<NSString *> *_fdzChannels;
    NSMutableDictionary<NSString *, TTIOFDZBlockTable *> *_fdzTables;
    NSMutableDictionary<NSString *, NSData *> *_fdzBlockCache;
    NSMutableDictionary<NSString *, NSNumber *> *_fdzBlockCacheIndex;

    // Persistence context attached post-load for protocol encryption
    NSString *_persistenceFilePath;
    NSString *_persistenceRunName;

    // chromatogram traces carried with this run.
    NSArray<TTIOChromatogram *> *_chromatograms;

    // Vibrational-spectrum run metadata (IR / Raman / UV-Vis), parity
    // with Python/Java. Captured from the first spectrum on the
    // in-memory path and from run attributes by the readers; consumed by
    // spectrumAtIndex: to rebuild the right subclass. UV-Vis solvent
    // reuses _solvent.
    TTIOIRMode  _irMode;
    double      _irResolutionCmInv;
    NSUInteger  _irNumberOfScans;
    double      _ramanExcitationWavelengthNm;
    double      _ramanLaserPowerMw;
    double      _ramanIntegrationTimeSec;
    double      _uvvisPathLengthCm;
}

@synthesize chromatograms = _chromatograms;
@synthesize name = _name;

#pragma mark - Construction

- (instancetype)initWithSpectra:(NSArray *)spectra
                acquisitionMode:(TTIOAcquisitionMode)mode
               instrumentConfig:(TTIOInstrumentConfig *)config
{
    self = [super init];
    if (self) {
        _name             = @"";
        _inMemorySpectra  = [spectra copy];
        _acquisitionMode  = mode;
        _instrumentConfig = config;
        _streamPosition   = 0;
        _provenance       = [NSMutableArray array];
        _signalCompression = TTIOCompressionZlib;  // M21 default

        if (spectra.count > 0) {
            TTIOSpectrum *first = spectra[0];
            _spectrumClassName = NSStringFromClass([first class]);
            _channelNames = [[first.signalArrays allKeys]
                sortedArrayUsingSelector:@selector(compare:)];

            if ([first isKindOfClass:[TTIONMRSpectrum class]]) {
                TTIONMRSpectrum *n = (TTIONMRSpectrum *)first;
                _nucleusType = [n.nucleusType copy];
                _spectrometerFrequencyMHz = n.spectrometerFrequencyMHz;
            } else if ([first isKindOfClass:[TTIOIRSpectrum class]]) {
                TTIOIRSpectrum *ir = (TTIOIRSpectrum *)first;
                _irMode = ir.mode;
                _irResolutionCmInv = ir.resolutionCmInv;
                _irNumberOfScans = ir.numberOfScans;
            } else if ([first isKindOfClass:[TTIORamanSpectrum class]]) {
                TTIORamanSpectrum *r = (TTIORamanSpectrum *)first;
                _ramanExcitationWavelengthNm = r.excitationWavelengthNm;
                _ramanLaserPowerMw = r.laserPowerMw;
                _ramanIntegrationTimeSec = r.integrationTimeSec;
            } else if ([first isKindOfClass:[TTIOUVVisSpectrum class]]) {
                TTIOUVVisSpectrum *u = (TTIOUVVisSpectrum *)first;
                _uvvisPathLengthCm = u.pathLengthCm;
                _solvent = [u.solvent copy];
            }
        } else {
            _spectrumClassName = @"TTIOMassSpectrum";
            _channelNames = @[@"mz", @"intensity"];
        }

        _spectrumIndex = [self buildIndexFromSpectra:spectra];
        _chromatograms = @[];
        _modality = @"mass_spectrometry";
        // Default solvent to empty unless a UV-Vis first spectrum set it.
        if (!_solvent) _solvent = @"";
    }
    return self;
}

- (instancetype)initWithSpectra:(NSArray *)spectra
                  chromatograms:(NSArray<TTIOChromatogram *> *)chromatograms
                acquisitionMode:(TTIOAcquisitionMode)mode
               instrumentConfig:(TTIOInstrumentConfig *)config
{
    self = [self initWithSpectra:spectra acquisitionMode:mode instrumentConfig:config];
    if (self) {
        _chromatograms = chromatograms ? [chromatograms copy] : @[];
    }
    return self;
}

#pragma mark - Index construction

- (TTIOSpectrumIndex *)buildIndexFromSpectra:(NSArray *)spectra
{
    NSUInteger n = spectra.count;
    NSMutableData *offsets = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    NSMutableData *lengths = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *rts     = [NSMutableData dataWithLength:n * sizeof(double)];
    NSMutableData *ml      = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *pol     = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *pmz     = [NSMutableData dataWithLength:n * sizeof(double)];
    NSMutableData *pc      = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    NSMutableData *bp      = [NSMutableData dataWithLength:n * sizeof(double)];

    uint64_t *off = offsets.mutableBytes;
    uint32_t *len = lengths.mutableBytes;
    double   *rt  = rts.mutableBytes;
    int32_t  *mlp = ml.mutableBytes;
    int32_t  *plp = pol.mutableBytes;
    double   *pmp = pmz.mutableBytes;
    int32_t  *pcp = pc.mutableBytes;
    double   *bpp = bp.mutableBytes;

    // scan once to see if any MS spectrum carries activation/isolation
    // detail. If so, build four parallel optional columns; otherwise pass nil
    // so the index reflects the legacy (no opt_ms2_activation_detail) layout.
    BOOL anyM74 = NO;
    for (NSUInteger i = 0; i < n; i++) {
        TTIOSpectrum *s = spectra[i];
        if (![s isKindOfClass:[TTIOMassSpectrum class]]) continue;
        TTIOMassSpectrum *ms = (TTIOMassSpectrum *)s;
        if (ms.activationMethod != TTIOActivationMethodNone ||
            ms.isolationWindow != nil) { anyM74 = YES; break; }
    }
    NSMutableData *actM = nil;
    NSMutableData *isoT = nil;
    NSMutableData *isoL = nil;
    NSMutableData *isoU = nil;
    int32_t *actMp = NULL;
    double  *isoTp = NULL;
    double  *isoLp = NULL;
    double  *isoUp = NULL;
    if (anyM74) {
        actM = [NSMutableData dataWithLength:n * sizeof(int32_t)];
        isoT = [NSMutableData dataWithLength:n * sizeof(double)];
        isoL = [NSMutableData dataWithLength:n * sizeof(double)];
        isoU = [NSMutableData dataWithLength:n * sizeof(double)];
        actMp = actM.mutableBytes;
        isoTp = isoT.mutableBytes;
        isoLp = isoL.mutableBytes;
        isoUp = isoU.mutableBytes;
    }

    NSString *firstChannel = _channelNames.firstObject;
    uint64_t cursor = 0;
    for (NSUInteger i = 0; i < n; i++) {
        TTIOSpectrum *s = spectra[i];
        TTIOSignalArray *primary = s.signalArrays[firstChannel];
        off[i] = cursor;
        len[i] = (uint32_t)primary.length;
        rt[i]  = s.scanTimeSeconds;
        pmp[i] = s.precursorMz;
        pcp[i] = (int32_t)s.precursorCharge;

        if ([s isKindOfClass:[TTIOMassSpectrum class]]) {
            TTIOMassSpectrum *ms = (TTIOMassSpectrum *)s;
            mlp[i] = (int32_t)ms.msLevel;
            plp[i] = (int32_t)ms.polarity;

            double maxI = 0;
            TTIOSignalArray *inA = ms.intensityArray;
            /* -float64Buffer hands back a fresh conversion buffer for
             * any precision other than float64, and nothing else holds
             * it: keep it alive across the scan. */
            NSData *intBuf = [inA float64Buffer];
            const double *intP = intBuf.bytes;
            NSUInteger m = inA.length;
            for (NSUInteger j = 0; j < m; j++) if (intP[j] > maxI) maxI = intP[j];
            bpp[i] = maxI;

            if (anyM74) {
                actMp[i] = (int32_t)ms.activationMethod;
                TTIOIsolationWindow *iw = ms.isolationWindow;
                if (iw) {
                    isoTp[i] = iw.targetMz;
                    isoLp[i] = iw.lowerOffset;
                    isoUp[i] = iw.upperOffset;
                } else {
                    isoTp[i] = 0.0;
                    isoLp[i] = 0.0;
                    isoUp[i] = 0.0;
                }
            }
        } else {
            // NMR or other non-MS spectra: sentinel values.
            mlp[i] = 0;
            plp[i] = (int32_t)TTIOPolarityUnknown;

            double maxI = 0;
            TTIOSignalArray *inA = s.signalArrays[@"intensity"];
            if (inA) {
                NSData *intBuf = [inA float64Buffer];
                const double *intP = intBuf.bytes;
                NSUInteger m = inA.length;
                for (NSUInteger j = 0; j < m; j++) if (intP[j] > maxI) maxI = intP[j];
            }
            bpp[i] = maxI;

            if (anyM74) {
                actMp[i] = (int32_t)TTIOActivationMethodNone;
                isoTp[i] = 0.0;
                isoLp[i] = 0.0;
                isoUp[i] = 0.0;
            }
        }

        cursor += primary.length;
    }
    return [[TTIOSpectrumIndex alloc] initWithOffsets:offsets
                                              lengths:lengths
                                       retentionTimes:rts
                                             msLevels:ml
                                           polarities:pol
                                         precursorMzs:pmz
                                     precursorCharges:pc
                                  basePeakIntensities:bp
                                    activationMethods:actM
                                   isolationTargetMzs:isoT
                                isolationLowerOffsets:isoL
                                isolationUpperOffsets:isoU];
}

#pragma mark - Pixel coordinates (M102)

- (BOOL)setPixelX:(NSData *)pixelX
           pixelY:(NSData *)pixelY
           pixelZ:(NSData *)pixelZ
            error:(NSError **)error
{
    TTIOSpectrumIndex *idx = [_spectrumIndex indexWithPixelX:pixelX
                                                      pixelY:pixelY
                                                      pixelZ:pixelZ
                                                       error:error];
    if (!idx) return NO;
    _spectrumIndex = idx;
    return YES;
}

#pragma mark - HDF5 write

- (BOOL)writeToGroup:(id<TTIOStorageGroup>)parent name:(NSString *)name error:(NSError **)error
{
    NSParameterAssert(_inMemorySpectra != nil);  // disk-backed runs are read-only

    id<TTIOStorageGroup> runGroup = [parent createGroupNamed:name error:error];
    if (!runGroup) return NO;

    if (![runGroup setAttributeValue:@((int64_t)_acquisitionMode)
                             forName:@"acquisition_mode" error:error]) return NO;
    if (![runGroup setAttributeValue:@((int64_t)_inMemorySpectra.count)
                             forName:@"spectrum_count" error:error]) return NO;
    if (![runGroup setAttributeValue:_spectrumClassName
                             forName:@"spectrum_class" error:error]) return NO;

    if (_nucleusType) {
        if (![runGroup setAttributeValue:_nucleusType
                                 forName:@"nucleus_type" error:error]) return NO;
        id<TTIOStorageDataset> fd = [runGroup createDatasetNamed:@"_spectrometer_freq_mhz"
                                                       precision:TTIOPrecisionFloat64
                                                          length:1
                                                       chunkSize:0
                                                     compression:TTIOCompressionZlib
                                                compressionLevel:0
                                                           error:error];
        if (!fd) return NO;
        double f[1] = { _spectrometerFrequencyMHz };
        if (![fd writeAll:[NSData dataWithBytes:f length:sizeof(f)] error:error]) return NO;
    }

    if (_solvent && _solvent.length > 0) {
        if (![runGroup setAttributeValue:_solvent
                                 forName:@"solvent" error:error]) return NO;
    }

    // Vibrational-spectrum run metadata (parity with Python/Java
    // _write_run). Emitted only for the matching class so MS/NMR runs
    // stay byte-identical; ir_mode is always written for IR (0 =
    // transmittance is meaningful), the float/scan fields only when set.
    // _spectrumClassName remains the persisted source of truth; the enum
    // is an in-code dispatch key only (P3.8).
    TTIOSpectrumKind k = TTIOSpectrumKindFromPersisted(_spectrumClassName);
    if (k == TTIOSpectrumKindIR) {
        if (![runGroup setAttributeValue:@((int64_t)_irMode)
                                 forName:@"ir_mode" error:error]) return NO;
        if (_irResolutionCmInv != 0.0 &&
            ![runGroup setAttributeValue:@(_irResolutionCmInv)
                                 forName:@"ir_resolution_cm_inv" error:error]) return NO;
        if (_irNumberOfScans != 0 &&
            ![runGroup setAttributeValue:@((int64_t)_irNumberOfScans)
                                 forName:@"ir_number_of_scans" error:error]) return NO;
    } else if (k == TTIOSpectrumKindRaman) {
        if (_ramanExcitationWavelengthNm != 0.0 &&
            ![runGroup setAttributeValue:@(_ramanExcitationWavelengthNm)
                                 forName:@"raman_excitation_wavelength_nm" error:error]) return NO;
        if (_ramanLaserPowerMw != 0.0 &&
            ![runGroup setAttributeValue:@(_ramanLaserPowerMw)
                                 forName:@"raman_laser_power_mw" error:error]) return NO;
        if (_ramanIntegrationTimeSec != 0.0 &&
            ![runGroup setAttributeValue:@(_ramanIntegrationTimeSec)
                                 forName:@"raman_integration_time_sec" error:error]) return NO;
    } else if (k == TTIOSpectrumKindUVVis) {
        if (_uvvisPathLengthCm != 0.0 &&
            ![runGroup setAttributeValue:@(_uvvisPathLengthCm)
                                 forName:@"uvvis_path_length_cm" error:error]) return NO;
    }

    // Per-run provenance.
    //
    // v0.3 writes the records as a compound HDF5 dataset at
    //     /study/ms_runs/<run>/provenance/steps
    // using the same compound type as the dataset-level `/study/provenance`
    // (see TTIOCompoundIO). The `compound_per_run_provenance` feature flag
    // on the root group advertises this layout.
    //
    // For backward compatibility with v0.2 readers (including the in-tree
    // signature manager which still operates on the JSON blob), the writer
    // keeps `@provenance_json` as a legacy mirror. M18 will replace the
    // mirror with a canonical-byte-order signature path that covers the
    // compound dataset directly; until then the mirror is intentional.
    if (_provenance.count > 0) {
        if (![[self class] writeProvenance:_provenance toRunGroup:runGroup error:error]) return NO;
    }

    if (![_instrumentConfig writeToGroup:runGroup error:error]) return NO;
    if (![_spectrumIndex    writeToGroup:runGroup error:error]) return NO;

    id<TTIOStorageGroup> channels = [runGroup createGroupNamed:@"signal_channels" error:error];
    if (!channels) return NO;

    NSString *namesJoined = [_channelNames componentsJoinedByString:@","];
    if (![channels setAttributeValue:namesJoined
                             forName:@"channel_names" error:error]) return NO;

    NSUInteger total = 0;
    for (TTIOSpectrum *s in _inMemorySpectra) {
        total += [s.signalArrays[_channelNames.firstObject] length];
    }

    // Phase 2: Zlib left at its default resolves to codec 17 on MS
    // runs unless the caller opted out.
    TTIOCompression effectiveCompression = _signalCompression;
    if (_signalCompression == TTIOCompressionZlib
        && !_optDisableFloatDelta
        && [_spectrumClassName isEqualToString:@"TTIOMassSpectrum"]) {
        effectiveCompression = TTIOCompressionFloatDeltaZstd;
    }

    for (NSString *chName in _channelNames) {
        // Concat per-spectrum channel buffers into one flat NSData.
        // NSMutableData would zero-fill the backing store on
        // -dataWithLength: before we memcpy over it; allocate a bare
        // C buffer instead and hand it to NSData via -dataWithBytesNoCopy:
        // so the compressor sees a single contiguous region without the
        // zero-fill / NSMutableData bookkeeping tax. ~3× faster concat
        // on 100K-spectrum runs.
        NSUInteger totalBytes = total * sizeof(double);
        void *raw = malloc(totalBytes);
        if (!raw) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetCreate,
                @"out of memory concatenating signal channel '%@'", chName);
            return NO;
        }
        NSUInteger cursor = 0;
        for (TTIOSpectrum *s in _inMemorySpectra) {
            TTIOSignalArray *arr = s.signalArrays[chName];
            NSUInteger n = arr.length;
            memcpy((uint8_t *)raw + cursor * sizeof(double),
                   [arr float64Buffer].bytes, n * sizeof(double));
            cursor += n;
        }
        NSData *all = [NSData dataWithBytesNoCopy:raw length:totalBytes freeWhenDone:YES];
        NSString *dsName = [chName stringByAppendingString:@"_values"];

        if (_signalCompression == TTIOCompressionNumpressDelta) {
            // Fixed-point + first-difference transform. The dataset
            // stores int64 deltas; the reader detects the
            // ``@numpress_fixed_point`` attribute and reverses.
            const double *src = (const double *)all.bytes;
            double minV = src[0], maxV = src[0];
            for (NSUInteger k = 1; k < total; k++) {
                if (src[k] < minV) minV = src[k];
                if (src[k] > maxV) maxV = src[k];
            }
            int64_t scale = [TTIONumpress scaleForValueRangeMin:minV max:maxV];
            NSMutableData *deltas = [NSMutableData dataWithLength:total * sizeof(int64_t)];
            if (![TTIONumpress encodeFloat64:src
                                        count:total
                                        scale:scale
                                    outDeltas:(int64_t *)deltas.mutableBytes]) {
                if (error) *error = TTIOMakeError(TTIOErrorDatasetCreate,
                    @"numpress encode failed for '%@'", dsName);
                return NO;
            }
            id<TTIOStorageDataset> ds =
                [channels createDatasetNamed:dsName
                                   precision:TTIOPrecisionInt64
                                      length:total
                                   chunkSize:65536
                                 compression:TTIOCompressionZlib
                            compressionLevel:6
                                       error:error];
            if (!ds) return NO;
            if (![ds writeAll:deltas error:error]) return NO;
            NSString *spAttr =
                [NSString stringWithFormat:@"%@_numpress_fixed_point", chName];
            if (![channels setAttributeValue:@(scale)
                                     forName:spAttr error:error]) return NO;
        } else if (effectiveCompression == TTIOCompressionFloatDeltaZstd) {
            /* Codec id 17: the dataset bytes ARE the FDZ1 stream;
               @compression on the dataset is the dispatch signal and
               no HDF5 filter is applied. */
            NSData *stream = [TTIOFloatDeltaZstd encodeFloat64:all];
            if (!stream) {
                if (error) *error = TTIOMakeError(TTIOErrorDatasetCreate,
                    @"FLOAT_DELTA_ZSTD encode failed for '%@'", dsName);
                return NO;
            }
            id<TTIOStorageDataset> ds =
                [channels createDatasetNamed:dsName
                                   precision:TTIOPrecisionUInt8
                                      length:stream.length
                                   chunkSize:65536
                                 compression:TTIOCompressionNone
                            compressionLevel:0
                                       error:error];
            if (!ds) return NO;
            if (![ds writeAll:stream error:error]) return NO;
            if (![ds setAttributeValue:@(TTIOCompressionFloatDeltaZstd)
                               forName:@"compression" error:error]) return NO;
        } else {
            id<TTIOStorageDataset> ds =
                [channels createDatasetNamed:dsName
                                   precision:TTIOPrecisionFloat64
                                      length:total
                                   chunkSize:65536
                                 compression:_signalCompression
                            compressionLevel:6
                                       error:error];
            if (!ds) return NO;
            if (![ds writeAll:all error:error]) return NO;
        }
    }

    // chromatograms under <run>/chromatograms/
    if (_chromatograms.count > 0) {
        if (![[self class] writeChromatograms:_chromatograms toRunGroup:runGroup error:error]) return NO;
    }

    return YES;
}

+ (BOOL)writeProvenance:(NSArray<TTIOProvenanceRecord *> *)records
             toRunGroup:(id<TTIOStorageGroup>)runGroup
                  error:(NSError **)error
{
    id<TTIOStorageGroup> provGroup =
        [runGroup createGroupNamed:@"provenance" error:error];
    if (!provGroup) return NO;
    if (![TTIOCompoundIO writeProvenance:records
                               intoGroup:provGroup
                            datasetNamed:@"steps"
                                   error:error]) return NO;

    NSMutableArray *plists = [NSMutableArray arrayWithCapacity:records.count];
    for (TTIOProvenanceRecord *r in records) [plists addObject:[r asPlist]];
    NSError *jErr = nil;
    NSData *json = [NSJSONSerialization dataWithJSONObject:plists options:0 error:&jErr];
    if (!json) {
        if (error) *error = jErr;
        return NO;
    }
    NSString *jstr = [[NSString alloc] initWithData:json encoding:NSUTF8StringEncoding];
    return [runGroup setAttributeValue:jstr forName:@"provenance_json" error:error];
}

// M24 helper — lays out /chromatograms/ with concatenated time/intensity
// datasets and a chromatogram_index/ subgroup of parallel metadata.
// chromatogram_index/offsets is omitted on disk; readers
// compute it from cumsum(lengths).
+ (BOOL)writeChromatograms:(NSArray<TTIOChromatogram *> *)_chromatograms
                toRunGroup:(id<TTIOStorageGroup>)runGroup
                     error:(NSError **)error
{
    NSUInteger nChroms = _chromatograms.count;
    id<TTIOStorageGroup> chromGroup =
        [runGroup createGroupNamed:@"chromatograms" error:error];
    if (!chromGroup) return NO;
    if (![chromGroup setAttributeValue:@((int64_t)nChroms)
                               forName:@"count" error:error]) return NO;

    NSUInteger totalPoints = 0;
    for (TTIOChromatogram *c in _chromatograms) totalPoints += c.timeArray.length;

    NSMutableData *timeAll = [NSMutableData dataWithLength:totalPoints * sizeof(double)];
    NSMutableData *intAll  = [NSMutableData dataWithLength:totalPoints * sizeof(double)];

    uint32_t *lengths      = calloc(nChroms, sizeof(uint32_t));
    int32_t  *types        = calloc(nChroms, sizeof(int32_t));
    double   *targetMzs    = calloc(nChroms, sizeof(double));
    double   *precursorMzs = calloc(nChroms, sizeof(double));
    double   *productMzs   = calloc(nChroms, sizeof(double));

    NSUInteger cursor = 0;
    for (NSUInteger i = 0; i < nChroms; i++) {
        TTIOChromatogram *c = _chromatograms[i];
        NSUInteger n = c.timeArray.length;
        memcpy((uint8_t *)timeAll.mutableBytes + cursor * sizeof(double),
               c.timeArray.buffer.bytes, n * sizeof(double));
        memcpy((uint8_t *)intAll.mutableBytes + cursor * sizeof(double),
               c.intensityArray.buffer.bytes, n * sizeof(double));
        lengths[i]      = (uint32_t)n;
        types[i]        = (int32_t)c.type;
        targetMzs[i]    = c.targetMz;
        precursorMzs[i] = c.precursorProductMz;
        productMzs[i]   = c.productMz;
        cursor += n;
    }

    BOOL ok = YES;

    #define WRITE_DS(_grp, _dname, _prec, _nelem, _data) do { \
        id<TTIOStorageDataset> _ds = [(_grp) createDatasetNamed:(_dname) \
                                                      precision:(_prec) \
                                                         length:(_nelem) \
                                                      chunkSize:0 \
                                                    compression:TTIOCompressionZlib \
                                               compressionLevel:0 \
                                                          error:error]; \
        if (!_ds) { ok = NO; break; } \
        if (![_ds writeAll:(_data) error:error]) { ok = NO; break; } \
    } while (0)

    do {
        WRITE_DS(chromGroup, @"time_values",      TTIOPrecisionFloat64, totalPoints, timeAll);
        WRITE_DS(chromGroup, @"intensity_values", TTIOPrecisionFloat64, totalPoints, intAll);

        id<TTIOStorageGroup> idx =
            [chromGroup createGroupNamed:@"chromatogram_index" error:error];
        if (!idx) { ok = NO; break; }

        WRITE_DS(idx, @"lengths",      TTIOPrecisionUInt32,  nChroms,
                 [NSData dataWithBytesNoCopy:lengths      length:nChroms*sizeof(uint32_t) freeWhenDone:NO]);
        WRITE_DS(idx, @"types",        TTIOPrecisionInt32,   nChroms,
                 [NSData dataWithBytesNoCopy:types        length:nChroms*sizeof(int32_t)  freeWhenDone:NO]);
        WRITE_DS(idx, @"target_mzs",   TTIOPrecisionFloat64, nChroms,
                 [NSData dataWithBytesNoCopy:targetMzs    length:nChroms*sizeof(double)   freeWhenDone:NO]);
        WRITE_DS(idx, @"precursor_mzs",TTIOPrecisionFloat64, nChroms,
                 [NSData dataWithBytesNoCopy:precursorMzs length:nChroms*sizeof(double)   freeWhenDone:NO]);
        WRITE_DS(idx, @"product_mzs",  TTIOPrecisionFloat64, nChroms,
                 [NSData dataWithBytesNoCopy:productMzs   length:nChroms*sizeof(double)   freeWhenDone:NO]);
    } while (0);

    #undef WRITE_DS

    free(lengths); free(types);
    free(targetMzs); free(precursorMzs); free(productMzs);
    return ok;
}

#pragma mark - HDF5 read

/** Restore vibrational-spectrum run metadata (IR / Raman / UV-Vis) from
 *  run-group attributes into the ivars, keyed by the spectrum_class.
 *  Shared by both readers; parity with Python AcquisitionRun.open and
 *  Java AcquisitionRun.readFrom. No-op for MS / NMR runs. */
+ (void)loadVibrationalMetadataInto:(TTIOAcquisitionRun *)run
                               from:(id<TTIOStorageGroup>)runGroup
                          className:(NSString *)className
{
    // className is the persisted spectrum_class string (source of truth);
    // the enum is an in-code dispatch key only (P3.8).
    TTIOSpectrumKind k = TTIOSpectrumKindFromPersisted(className);
    if (k == TTIOSpectrumKindIR) {
        id m = [runGroup attributeValueForName:@"ir_mode" error:NULL];
        run->_irMode = (TTIOIRMode)
            ([m respondsToSelector:@selector(longLongValue)] ? [m longLongValue] : 0);
        run->_irResolutionCmInv = [self doubleAttr:runGroup name:@"ir_resolution_cm_inv"];
        id sc = [runGroup attributeValueForName:@"ir_number_of_scans" error:NULL];
        run->_irNumberOfScans = (NSUInteger)
            ([sc respondsToSelector:@selector(longLongValue)] ? [sc longLongValue] : 0);
    } else if (k == TTIOSpectrumKindRaman) {
        run->_ramanExcitationWavelengthNm =
            [self doubleAttr:runGroup name:@"raman_excitation_wavelength_nm"];
        run->_ramanLaserPowerMw = [self doubleAttr:runGroup name:@"raman_laser_power_mw"];
        run->_ramanIntegrationTimeSec =
            [self doubleAttr:runGroup name:@"raman_integration_time_sec"];
    } else if (k == TTIOSpectrumKindUVVis) {
        run->_uvvisPathLengthCm = [self doubleAttr:runGroup name:@"uvvis_path_length_cm"];
    }
}

+ (double)doubleAttr:(id<TTIOStorageGroup>)runGroup name:(NSString *)name
{
    id v = [runGroup attributeValueForName:name error:NULL];
    return [v respondsToSelector:@selector(doubleValue)] ? [v doubleValue] : 0.0;
}

+ (instancetype)readFromStorageGroup:(id)parent
                                 name:(NSString *)name
                                error:(NSError **)error
{
    id<TTIOStorageGroup> par = (id<TTIOStorageGroup>)parent;
    if (![par hasChildNamed:name]) return nil;
    id<TTIOStorageGroup> runGroup = [par openGroupNamed:name error:error];
    if (!runGroup) return nil;

    id modeObj = [runGroup attributeValueForName:@"acquisition_mode" error:NULL];
    TTIOAcquisitionMode mode = (TTIOAcquisitionMode)
        ([modeObj respondsToSelector:@selector(longLongValue)]
            ? [modeObj longLongValue] : 0);
    id classObj = [runGroup attributeValueForName:@"spectrum_class" error:NULL];
    NSString *className = [classObj isKindOfClass:[NSString class]]
        ? (NSString *)classObj : @"TTIOMassSpectrum";
    id nucObj = [runGroup attributeValueForName:@"nucleus_type" error:NULL];
    NSString *nucleus = [nucObj isKindOfClass:[NSString class]] ? nucObj : nil;
    // @modality fallback to "mass_spectrometry" so pre-v0.11
    // runs read back as mass-spec.
    id modObj = [runGroup attributeValueForName:@"modality" error:NULL];
    NSString *modality = ([modObj isKindOfClass:[NSString class]]
                          && [(NSString *)modObj length] > 0)
        ? (NSString *)modObj : @"mass_spectrometry";

    TTIOSpectrumIndex *idx = [TTIOSpectrumIndex readFromStorageGroup:runGroup error:error];
    if (!idx) return nil;

    // signal_channels: read channel_names attr; full signal read is
    // deferred to HDF5-only paths for v0.9.
    NSArray<NSString *> *channelNames = @[];
    if ([runGroup hasChildNamed:@"signal_channels"]) {
        id<TTIOStorageGroup> sc = [runGroup openGroupNamed:@"signal_channels" error:NULL];
        id names = [sc attributeValueForName:@"channel_names" error:NULL];
        if ([names isKindOfClass:[NSString class]]) {
            channelNames = [(NSString *)names componentsSeparatedByString:@","];
        }
    }

    // Provenance via the JSON mirror (compound-dataset decode is HDF5-only).
    NSMutableArray<TTIOProvenanceRecord *> *provenance = [NSMutableArray array];
    id provObj = [runGroup attributeValueForName:@"provenance_json" error:NULL];
    if ([provObj isKindOfClass:[NSString class]] && [(NSString *)provObj length] > 0) {
        NSData *jdata = [(NSString *)provObj dataUsingEncoding:NSUTF8StringEncoding];
        NSArray *plists = [NSJSONSerialization JSONObjectWithData:jdata options:0 error:NULL];
        for (NSDictionary *p in plists) {
            TTIOProvenanceRecord *r = [TTIOProvenanceRecord fromPlist:p];
            if (r) [provenance addObject:r];
        }
    }

    // Default InstrumentConfig; non-HDF5 writers don't persist it today.
    TTIOInstrumentConfig *cfg = [[TTIOInstrumentConfig alloc] initWithManufacturer:@""
                                                                             model:@""
                                                                      serialNumber:@""
                                                                        sourceType:@""
                                                                      analyzerType:@""
                                                                      detectorType:@""];

    TTIOAcquisitionRun *run = [[self alloc] init];
    run->_name                 = [name copy] ?: @"";
    run->_acquisitionMode      = mode;
    run->_instrumentConfig     = cfg;
    run->_spectrumIndex        = idx;
    run->_storageSignalGroup   = nil;
    run->_storageDatasets      = nil;
    run->_channelNames         = [channelNames copy];
    run->_spectrumClassName    = [className copy];
    run->_nucleusType          = [nucleus copy];
    run->_spectrometerFrequencyMHz = 0.0;
    run->_inMemorySpectra      = nil;
    run->_streamPosition       = 0;
    run->_provenance           = provenance;
    run->_numpressChannels     = nil;
    run->_signalCompression    = TTIOCompressionNone;
    run->_chromatograms        = @[];
    run->_modality             = [modality copy];
    // solvent (UV-Vis / NMR label) — read it so UV-Vis runs round-trip
    // through this lighter protocol reader too.
    id solvObj = [runGroup attributeValueForName:@"solvent" error:NULL];
    run->_solvent = ([solvObj isKindOfClass:[NSString class]]
                     && [(NSString *)solvObj length] > 0)
        ? [(NSString *)solvObj copy] : @"";
    [self loadVibrationalMetadataInto:run from:runGroup className:className];
    return run;
}

+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)parent name:(NSString *)name error:(NSError **)error
{
    id<TTIOStorageGroup> runGroup = [parent openGroupNamed:name error:error];
    if (!runGroup) return nil;

    NSNumber *modeNum = [runGroup attributeValueForName:@"acquisition_mode" error:error];
    TTIOAcquisitionMode mode = modeNum ? (TTIOAcquisitionMode)[modeNum longLongValue] : 0;

    TTIOInstrumentConfig *cfg = [TTIOInstrumentConfig readFromGroup:runGroup error:error];
    if (!cfg) return nil;

    TTIOSpectrumIndex *idx = [TTIOSpectrumIndex readFromGroup:runGroup error:error];
    if (!idx) return nil;

    // v0.2 additions; v0.1 fallback if missing.
    NSString *className = @"TTIOMassSpectrum";
    if ([runGroup hasAttributeNamed:@"spectrum_class"]) {
        NSString *cn = [runGroup attributeValueForName:@"spectrum_class" error:NULL];
        if (cn.length > 0) className = cn;
    }

    // @modality with pre-v0.11 mass-spec fallback.
    NSString *modality = @"mass_spectrometry";
    if ([runGroup hasAttributeNamed:@"modality"]) {
        NSString *m = [runGroup attributeValueForName:@"modality" error:NULL];
        if (m.length > 0) modality = m;
    }

    NSString *nucleus = nil;
    double freqMHz = 0.0;
    if ([runGroup hasAttributeNamed:@"nucleus_type"]) {
        nucleus = [runGroup attributeValueForName:@"nucleus_type" error:NULL];
        if ([runGroup hasChildNamed:@"_spectrometer_freq_mhz"]) {
            id<TTIOStorageDataset> fd = [runGroup openDatasetNamed:@"_spectrometer_freq_mhz" error:NULL];
            NSData *fdata = [fd readAll:NULL];
            if (fdata.length >= sizeof(double)) {
                freqMHz = ((const double *)fdata.bytes)[0];
            }
        }
    }

    NSString *solvent = @"";
    if ([runGroup hasAttributeNamed:@"solvent"]) {
        NSString *s = [runGroup attributeValueForName:@"solvent" error:NULL];
        if (s.length > 0) solvent = s;
    }

    // Per-run provenance: prefer the v0.3 compound layout at
    // runGroup/provenance/steps; fall back to the v0.2 @provenance_json
    // attribute if the compound subgroup is absent. Pre-v0.2 files had
    // neither form, in which case `provenance` remains an empty array.
    NSMutableArray<TTIOProvenanceRecord *> *provenance = [NSMutableArray array];
    if ([runGroup hasChildNamed:@"provenance"]) {
        id<TTIOStorageGroup> provGroup = [runGroup openGroupNamed:@"provenance" error:NULL];
        if (provGroup && [provGroup hasChildNamed:@"steps"]) {
            NSArray *compound =
                [TTIOCompoundIO readProvenanceFromGroup:provGroup
                                           datasetNamed:@"steps"
                                                  error:NULL];
            if (compound) [provenance addObjectsFromArray:compound];
        }
    }
    if (provenance.count == 0 && [runGroup hasAttributeNamed:@"provenance_json"]) {
        NSString *jstr = [runGroup attributeValueForName:@"provenance_json" error:NULL];
        NSData *jdata = [jstr dataUsingEncoding:NSUTF8StringEncoding];
        NSArray *plists = [NSJSONSerialization JSONObjectWithData:jdata
                                                           options:0
                                                             error:NULL];
        for (NSDictionary *p in plists) {
            TTIOProvenanceRecord *r = [TTIOProvenanceRecord fromPlist:p];
            if (r) [provenance addObject:r];
        }
    }

    id<TTIOStorageGroup> channels = [runGroup openGroupNamed:@"signal_channels" error:error];
    if (!channels) return nil;

    NSArray<NSString *> *channelNames = nil;
    if ([channels hasAttributeNamed:@"channel_names"]) {
        NSString *joined = [channels attributeValueForName:@"channel_names" error:NULL];
        channelNames = [joined componentsSeparatedByString:@","];
    } else {
        // v0.1 fallback
        channelNames = @[@"mz", @"intensity"];
    }

    // v0.7 M44 / Task 31: channelDatasets is a protocol-valued
    // dictionary so the hot-path read routes through TTIOStorageDataset.
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *channelDatasets =
        [NSMutableDictionary dictionaryWithCapacity:channelNames.count];
    NSMutableDictionary<NSString *, NSData *> *numpressChannels =
        [NSMutableDictionary dictionary];
    NSMutableSet<NSString *> *fdzChannels = [NSMutableSet set];
    TTIOCompression runCompression = TTIOCompressionZlib;
    for (NSString *chName in channelNames) {
        NSString *dsName = [chName stringByAppendingString:@"_values"];
        if (![channels hasChildNamed:dsName]) {
            // Channel is absent — most likely the file is encrypted
            // and this channel lives as `<name>_values_encrypted`. Keep
            // metadata load going; spectrumAtIndex: will error cleanly
            // if anyone later asks for data from this channel.
            continue;
        }

        // detect Numpress-delta encoding via the per-channel
        // ``@<chName>_numpress_fixed_point`` attribute.
        NSString *scaleAttr = [NSString stringWithFormat:@"%@_numpress_fixed_point", chName];
        if ([channels hasAttributeNamed:scaleAttr]) {
            NSNumber *scaleNum =
                [channels attributeValueForName:scaleAttr error:NULL];
            int64_t scale = scaleNum ? [scaleNum longLongValue] : 0;
            id<TTIOStorageDataset> ds =
                [channels openDatasetNamed:dsName error:error];
            if (!ds) return nil;
            NSData *raw = [ds readAll:error];
            if (!raw) return nil;
            NSUInteger nElems = raw.length / sizeof(int64_t);
            NSMutableData *decoded =
                [NSMutableData dataWithLength:nElems * sizeof(double)];
            if (![TTIONumpress decodeInt64:(const int64_t *)raw.bytes
                                       count:nElems
                                       scale:scale
                                  outValues:(double *)decoded.mutableBytes]) {
                if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
                    @"numpress decode failed for '%@'", chName);
                return nil;
            }
            numpressChannels[chName] = decoded;
            runCompression = TTIOCompressionNumpressDelta;
            continue;
        }

        id<TTIOStorageDataset> ds = [channels openDatasetNamed:dsName error:error];
        if (!ds) return nil;

        /* FLOAT_DELTA_ZSTD (codec id 17): the dataset is a flat uint8
           FDZ1 stream with @compression = 17. Decode once into the
           same eager buffer the numpress path uses. */
        NSNumber *codecNum = [ds attributeValueForName:@"compression" error:NULL];
        NSUInteger codecId = codecNum ? [codecNum unsignedIntegerValue] : 0;
        if (codecId == TTIOCompressionFloatDeltaZstd) {
            /* Left encoded: -channelRange:start:count:error: decodes
               block-wise through the FDZ1 block table. */
            [fdzChannels addObject:chName];
            channelDatasets[chName] = ds;
            runCompression = TTIOCompressionFloatDeltaZstd;
            continue;
        }
        if (codecId != 0) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
                @"signal channel '%@': @compression=%lu is not a spectral "
                @"channel codec (FLOAT_DELTA_ZSTD=17 is the only one wired)",
                chName, (unsigned long)codecId);
            return nil;
        }
        channelDatasets[chName] = ds;
    }

    TTIOAcquisitionRun *run = [[self alloc] init];
    run->_name                 = [name copy] ?: @"";
    run->_acquisitionMode      = mode;
    run->_instrumentConfig     = cfg;
    run->_spectrumIndex        = idx;
    run->_storageSignalGroup   = channels;
    run->_spectralBlockIndex   = [TTIOSpectralBlockIndex readFromRunGroup:runGroup
                                                                   error:NULL];
    run->_storageDatasets      = channelDatasets;
    run->_channelNames         = [channelNames copy];
    run->_spectrumClassName    = [className copy];
    run->_nucleusType          = [nucleus copy];
    run->_spectrometerFrequencyMHz = freqMHz;
    run->_modality             = [modality copy];
    run->_solvent              = [solvent copy];
    run->_inMemorySpectra      = nil;
    run->_streamPosition       = 0;
    run->_provenance           = provenance;
    run->_numpressChannels     = numpressChannels.count > 0 ? numpressChannels : nil;
    run->_fdzChannels          = fdzChannels.count > 0 ? fdzChannels : nil;
    run->_signalCompression    = runCompression;

    // read chromatograms if present. Absence means v0.3 file → empty list.
    run->_chromatograms = [self readChromatogramsFromRunGroup:runGroup];
    [self loadVibrationalMetadataInto:run from:runGroup className:className];
    return run;
}

+ (NSArray<TTIOChromatogram *> *)readChromatogramsFromRunGroup:(id<TTIOStorageGroup>)runGroup
{
    if (![runGroup hasChildNamed:@"chromatograms"]) return @[];
    id<TTIOStorageGroup> chromGroup = [runGroup openGroupNamed:@"chromatograms" error:NULL];
    if (!chromGroup) return @[];

    NSNumber *cntNum = [chromGroup attributeValueForName:@"count" error:NULL];
    if (!cntNum) return @[];
    int64_t count = [cntNum longLongValue];
    if (count <= 0) return @[];

    id<TTIOStorageDataset> timeDs = [chromGroup openDatasetNamed:@"time_values" error:NULL];
    id<TTIOStorageDataset> intDs  = [chromGroup openDatasetNamed:@"intensity_values" error:NULL];
    if (!timeDs || !intDs) return @[];
    NSData *timeAll = [timeDs readAll:NULL];
    NSData *intAll  = [intDs  readAll:NULL];
    if (!timeAll || !intAll) return @[];

    id<TTIOStorageGroup> idxGroup = [chromGroup openGroupNamed:@"chromatogram_index" error:NULL];
    if (!idxGroup) return @[];

    NSData *lengthsData   = [[idxGroup openDatasetNamed:@"lengths" error:NULL] readAll:NULL];
    if (!lengthsData) return @[];
    // offsets is omitted from disk by default; synthesize.
    NSData *offsetsData;
    if ([idxGroup hasChildNamed:@"offsets"]) {
        offsetsData = [[idxGroup openDatasetNamed:@"offsets" error:NULL] readAll:NULL];
    } else {
        offsetsData = TTIOOffsetsFromLengths(lengthsData);
    }
    NSData *typesData     = [[idxGroup openDatasetNamed:@"types" error:NULL] readAll:NULL];
    NSData *targetData    = [[idxGroup openDatasetNamed:@"target_mzs" error:NULL] readAll:NULL];
    NSData *precursorData = [[idxGroup openDatasetNamed:@"precursor_mzs" error:NULL] readAll:NULL];
    NSData *productData   = [[idxGroup openDatasetNamed:@"product_mzs" error:NULL] readAll:NULL];
    if (!offsetsData || !typesData ||
        !targetData || !precursorData || !productData) return @[];

    const int64_t  *offsets      = offsetsData.bytes;
    const uint32_t *lengths      = lengthsData.bytes;
    const int32_t  *types        = typesData.bytes;
    const double   *targetMzs    = targetData.bytes;
    const double   *precursorMzs = precursorData.bytes;
    const double   *productMzs   = productData.bytes;

    pthread_once(&gStdChannelEncodingOnce, _buildStdChannelEncoding);
    TTIOEncodingSpec *enc = gStdChannelEncoding;

    NSMutableArray<TTIOChromatogram *> *out =
        [NSMutableArray arrayWithCapacity:(NSUInteger)count];
    for (int64_t i = 0; i < count; i++) {
        NSUInteger off = (NSUInteger)offsets[i];
        NSUInteger len = (NSUInteger)lengths[i];
        NSData *tSlice = [NSData dataWithBytes:(const uint8_t *)timeAll.bytes + off*sizeof(double)
                                         length:len*sizeof(double)];
        NSData *iSlice = [NSData dataWithBytes:(const uint8_t *)intAll.bytes  + off*sizeof(double)
                                         length:len*sizeof(double)];
        TTIOSignalArray *tArr = [[TTIOSignalArray alloc] initWithBuffer:tSlice
                                                                  length:len
                                                                encoding:enc
                                                                    axis:nil];
        TTIOSignalArray *iArr = [[TTIOSignalArray alloc] initWithBuffer:iSlice
                                                                  length:len
                                                                encoding:enc
                                                                    axis:nil];
        TTIOChromatogram *c =
            [[TTIOChromatogram alloc] initWithTimeArray:tArr
                                          intensityArray:iArr
                                                    type:(TTIOChromatogramType)types[i]
                                                targetMz:targetMzs[i]
                                             precursorMz:precursorMzs[i]
                                               productMz:productMzs[i]
                                                   error:NULL];
        if (c) [out addObject:c];
    }
    return [out copy];
}

#pragma mark - Random access

// The whole column of a channel, decoded, from the memory caches or
// the dataset (readAll, or the full FDZ1 stream decoded), and kept in
// _cachedFullChannels. nil with *error on failure.
- (NSData *)_fullChannel:(NSString *)chName error:(NSError **)error
{
    NSData *plaintext = _decryptedChannels[chName];
    NSData *decoded = plaintext ?: _numpressChannels[chName];
    if (decoded) return decoded;
    NSData *full = _cachedFullChannels[chName];
    if (full) return full;
    @synchronized (self) {
        full = _cachedFullChannels[chName];
        if (!full) {
            id<TTIOStorageDataset> ds = _storageDatasets[chName];
            if (!ds) {
                if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
                    @"signal channel '%@' has no open dataset", chName);
                return nil;
            }
            full = [ds readAll:error];
            if (full && [_fdzChannels containsObject:chName]) {
                full = [TTIOFloatDeltaZstd decodeStream:full error:error];
            }
            if (full) {
                if (!_cachedFullChannels) _cachedFullChannels = [NSMutableDictionary dictionary];
                _cachedFullChannels[chName] = full;
            }
        }
    }
    return full;
}

/** As -_fdzRange:start:count:error: with block decodes on a pool: block
 *  bodies are read on the caller's thread, decoded concurrently, and
 *  appended in block order. */
- (NSData *)_fdzRange:(NSString *)chName
                start:(NSUInteger)start
                count:(NSUInteger)count
              threads:(NSUInteger)nthreads
                error:(NSError **)error
{
    id<TTIOStorageDataset> ds = _storageDatasets[chName];
    if (!ds) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"signal channel '%@' has no open dataset", chName);
        return nil;
    }
    TTIOFDZByteRangeReader reader = ^NSData *(NSUInteger offset, NSUInteger n) {
        id raw = [ds readSliceAtOffset:offset count:n error:NULL];
        return [raw isKindOfClass:[NSData class]] ? raw : nil;
    };
    TTIOFDZBlockTable *table;
    @synchronized (self) {
        if (!_fdzTables) _fdzTables = [NSMutableDictionary dictionary];
        table = _fdzTables[chName];
        if (!table) {
            table = [TTIOFloatDeltaZstd readBlockTableWithReader:reader error:error];
            if (!table) return nil;
            _fdzTables[chName] = table;
        }
    }
    if ((uint64_t)start + count > table.nValues) {
        if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
            @"channel '%@' range [%lu, %lu) beyond %llu values", chName,
            (unsigned long)start, (unsigned long)(start + count), (unsigned long long)table.nValues);
        return nil;
    }
    NSUInteger bs = table.blockSize;
    NSUInteger k0 = start / bs;
    NSUInteger k1 = (start + count - 1) / bs;
    if (nthreads <= 1 || count == 0 || k1 == k0) {
        return [self _fdzRange:chName start:start count:count error:error];
    }
    TTIOThreadPool *pool = [TTIOThreadPool poolWithThreads:nthreads];
    if (pool.queue == nil) {
        return [self _fdzRange:chName start:start count:count error:error];
    }
    NSMutableDictionary<NSNumber *, TTIOInFlightFdzDecode *> *pending = [NSMutableDictionary dictionary];
    NSCondition *cond = [NSCondition new];
    /* The one-block cache still serves a block the previous range call
     * already decoded. */
    NSData *cachedVals = nil;
    NSUInteger cachedK = NSNotFound;
    @synchronized (self) {
        NSNumber *ck = _fdzBlockCacheIndex[chName];
        if (ck != nil && _fdzBlockCache[chName] != nil
            && ck.unsignedIntegerValue >= k0 && ck.unsignedIntegerValue <= k1) {
            cachedK = ck.unsignedIntegerValue;
            cachedVals = _fdzBlockCache[chName];
        }
    }
    void (^submit)(NSUInteger) = ^(NSUInteger k) {
        if (k > k1 || pending[@(k)]) return;
        if (k == cachedK) {
            TTIOInFlightFdzDecode *f = [TTIOInFlightFdzDecode new];
            f.values = cachedVals;
            f.done = YES;
            pending[@(k)] = f;
            return;
        }
        TTIOInFlightFdzDecode *f = [TTIOInFlightFdzDecode new];
        pending[@(k)] = f;
        NSData *body = reader((NSUInteger)[table offsetAt:k], [table lengthAt:k]);   /* storage read, this thread */
        if (body.length != [table lengthAt:k]) { f.done = YES; return; }
        NSUInteger base = (NSUInteger)[table offsetAt:k];
        TTIOFDZByteRangeReader memReader = ^NSData *(NSUInteger offset, NSUInteger n) {
            if (offset < base || offset + n > base + body.length) return nil;
            return [body subdataWithRange:NSMakeRange(offset - base, n)];
        };
        [pool.queue addOperationWithBlock:^{
            NSData *v = nil;
            @try {
                v = [TTIOFloatDeltaZstd decodeBlock:k table:table reader:memReader error:NULL];
            } @catch (NSException *ex) {
                v = nil;
            }
            [cond lock];
            f.values = v;
            f.done = YES;
            [cond broadcast];
            [cond unlock];
        }];
    };
    for (NSUInteger k = k0; k <= MIN(k1, k0 + nthreads - 1); k++) submit(k);
    NSMutableData *out = [NSMutableData dataWithCapacity:count * sizeof(double)];
    NSUInteger pos = start, end = start + count;
    BOOL ok = YES;
    for (NSUInteger k = k0; k <= k1; k++) {
        TTIOInFlightFdzDecode *f = pending[@(k)];
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [pending removeObjectForKey:@(k)];
        if (!f.values) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"FLOAT_DELTA_ZSTD block %lu of '%@' failed to decode",
                (unsigned long)k, chName);
            ok = NO;
            break;
        }
        submit(k + nthreads);
        NSUInteger blockStart = k * bs;
        NSUInteger from = pos - blockStart;
        NSUInteger to = MIN(end - blockStart, f.values.length / sizeof(double));
        [out appendBytes:(const uint8_t *)f.values.bytes + from * sizeof(double)
                  length:(to - from) * sizeof(double)];
        pos = blockStart + to;
        if (k == k1) {
            @synchronized (self) {
                if (!_fdzBlockCache) _fdzBlockCache = [NSMutableDictionary dictionary];
                if (!_fdzBlockCacheIndex) _fdzBlockCacheIndex = [NSMutableDictionary dictionary];
                _fdzBlockCache[chName] = f.values;
                _fdzBlockCacheIndex[chName] = @(k);
            }
        }
    }
    for (TTIOInFlightFdzDecode *f in pending.allValues) {
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
    }
    [pool close];
    return ok ? out : nil;
}

- (NSData *)_fdzRange:(NSString *)chName
                start:(NSUInteger)start
                count:(NSUInteger)count
                error:(NSError **)error
{
    id<TTIOStorageDataset> ds = _storageDatasets[chName];
    if (!ds) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"signal channel '%@' has no open dataset", chName);
        return nil;
    }
    TTIOFDZByteRangeReader reader = ^NSData *(NSUInteger offset, NSUInteger n) {
        id raw = [ds readSliceAtOffset:offset count:n error:NULL];
        return [raw isKindOfClass:[NSData class]] ? raw : nil;
    };
    TTIOFDZBlockTable *table;
    @synchronized (self) {
        if (!_fdzTables) _fdzTables = [NSMutableDictionary dictionary];
        table = _fdzTables[chName];
        if (!table) {
            table = [TTIOFloatDeltaZstd readBlockTableWithReader:reader error:error];
            if (!table) return nil;
            _fdzTables[chName] = table;
        }
    }
    if ((uint64_t)start + count > table.nValues) {
        if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
            @"channel '%@' range [%lu, %lu) beyond %llu values", chName,
            (unsigned long)start, (unsigned long)(start + count), (unsigned long long)table.nValues);
        return nil;
    }
    NSMutableData *out = [NSMutableData dataWithCapacity:count * sizeof(double)];
    NSUInteger bs = table.blockSize;
    NSUInteger pos = start, end = start + count;
    while (pos < end) {
        NSUInteger k = pos / bs;
        NSData *block = nil;
        @synchronized (self) {
            if ([_fdzBlockCacheIndex[chName] unsignedIntegerValue] == k && _fdzBlockCache[chName]) {
                block = _fdzBlockCache[chName];
            }
        }
        if (!block) {
            block = [TTIOFloatDeltaZstd decodeBlock:k table:table reader:reader error:error];
            if (!block) return nil;
            @synchronized (self) {
                if (!_fdzBlockCache) _fdzBlockCache = [NSMutableDictionary dictionary];
                if (!_fdzBlockCacheIndex) _fdzBlockCacheIndex = [NSMutableDictionary dictionary];
                _fdzBlockCache[chName] = block;
                _fdzBlockCacheIndex[chName] = @(k);
            }
        }
        NSUInteger blockStart = k * bs;
        NSUInteger from = pos - blockStart;
        NSUInteger to = MIN(end - blockStart, block.length / sizeof(double));
        [out appendBytes:(const uint8_t *)block.bytes + from * sizeof(double) length:(to - from) * sizeof(double)];
        pos = blockStart + to;
    }
    return out;
}

- (NSData *)channelRange:(NSString *)chName
                  offset:(NSUInteger)start
                   count:(NSUInteger)count
                   error:(NSError **)error
{
    NSData *plaintext = _decryptedChannels[chName];
    NSData *decoded = plaintext ?: _numpressChannels[chName] ?: _cachedFullChannels[chName];
    if (decoded) {
        if ((start + count) * sizeof(double) > decoded.length) {
            if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
                @"channel '%@' range [%lu, %lu) beyond %lu values", chName,
                (unsigned long)start, (unsigned long)(start + count),
                (unsigned long)(decoded.length / sizeof(double)));
            return nil;
        }
        return [NSData dataWithBytes:(const uint8_t *)decoded.bytes + start * sizeof(double)
                              length:count * sizeof(double)];
    }
    if ([_fdzChannels containsObject:chName]) {
        return [self _fdzRange:chName start:start count:count error:error];
    }
    id<TTIOStorageDataset> ds = _storageDatasets[chName];
    if (!ds) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetOpen,
            @"signal channel '%@' has no open dataset", chName);
        return nil;
    }
    if (count == 0) return [NSData data];
    id raw = [ds readSliceAtOffset:start count:count error:error];
    if (![raw isKindOfClass:[NSData class]]) return nil;
    return [NSData dataWithData:raw];
}

- (NSData *)channelRange:(NSString *)channelName
                   start:(NSUInteger)start
                   count:(NSUInteger)count
                   error:(NSError **)error
{
    if (_inMemorySpectra) {
        NSMutableData *out = [NSMutableData dataWithCapacity:count * sizeof(double)];
        NSUInteger cursor = 0;
        for (TTIOSpectrum *sp in _inMemorySpectra) {
            TTIOSignalArray *a = sp.signalArrays[channelName];
            NSUInteger n = a.length;
            NSUInteger from = start > cursor ? start - cursor : 0;
            NSUInteger to = MIN(n, start + count > cursor ? start + count - cursor : 0);
            if (to > from) {
                [out appendBytes:(const uint8_t *)[a float64Buffer].bytes + from * sizeof(double)
                          length:(to - from) * sizeof(double)];
            }
            cursor += n;
            if (cursor >= start + count) break;
        }
        return out;
    }
    return [self channelRange:channelName offset:start count:count error:error];
}

- (NSData *)channelRange:(NSString *)channelName
                   start:(NSUInteger)start
                   count:(NSUInteger)count
                 threads:(NSUInteger)threads
                   error:(NSError **)error
{
    NSUInteger nthreads = threads ? threads : [TTIOThreads resolve:nil];
    if (nthreads > 1 && !_inMemorySpectra && !_decryptedChannels[channelName]
        && !_numpressChannels[channelName] && !_cachedFullChannels[channelName]
        && [_fdzChannels containsObject:channelName]) {
        return [self _fdzRange:channelName start:start count:count threads:nthreads error:error];
    }
    return [self channelRange:channelName start:start count:count error:error];
}

#pragma mark - Parallel block consumer

/* Units covering [from,to). A spectrum belongs to the block holding its
 * first value; see TTIOSpectralUnitPlan. */
- (void)_testDropBlockIndex { _spectralBlockIndex = nil; }

- (NSArray<NSValue *> *)_unitsFrom:(NSUInteger)from to:(NSUInteger)to
{
    NSUInteger n = [self count];
    NSUInteger hi = MIN(to, n);
    if (from >= hi) return @[];

    NSMutableData *offD = [NSMutableData dataWithLength:hi * sizeof(unsigned long long)];
    NSMutableData *lenD = [NSMutableData dataWithLength:hi * sizeof(unsigned int)];
    unsigned long long *off = offD.mutableBytes;
    unsigned int *len = lenD.mutableBytes;
    for (NSUInteger i = 0; i < hi; i++) {
        off[i] = (unsigned long long)[_spectrumIndex offsetAt:i];
        len[i] = (unsigned int)[_spectrumIndex lengthAt:i];
    }

    NSMutableData *bsD = nil;
    NSUInteger nBlocks = 0;

    TTIOSpectralBlockIndex *bi = _spectralBlockIndex;
    if (bi && bi.count > 0) {
        /* Tier 1: blocks/index. One compound read, no header walking. */
        nBlocks = bi.count;
        bsD = [NSMutableData dataWithLength:nBlocks * sizeof(unsigned long long)];
        unsigned long long *bs = bsD.mutableBytes;
        for (NSUInteger b = 0; b < nBlocks; b++) bs[b] = [bi valueStartAt:b];
    } else {
        /* Tier 2: the FDZ1 header walk. blocks/index is written only
         * when every channel cut at the same value boundaries, so with
         * it absent that has to be re-checked rather than assumed: the
         * tables must agree on block size and count, or the boundaries
         * are not shared and there is no single unit definition. */
        NSDictionary<NSString *, TTIOFDZBlockTable *> *tables =
            [self _fdzTablesForAllChannels];
        TTIOFDZBlockTable *first = nil;
        BOOL agree = (tables.count > 0);
        for (NSString *ch in tables) {
            TTIOFDZBlockTable *t = tables[ch];
            if (!first) { first = t; continue; }
            if (t.blockSize != first.blockSize || t.nBlocks != first.nBlocks) {
                agree = NO;
                break;
            }
        }
        if (agree && first && first.nBlocks > 0) {
            nBlocks = first.nBlocks;
            bsD = [NSMutableData dataWithLength:nBlocks * sizeof(unsigned long long)];
            unsigned long long *bs = bsD.mutableBytes;
            for (NSUInteger b = 0; b < nBlocks; b++) {
                bs[b] = (unsigned long long)b * first.blockSize;
            }
        }
    }

    if (nBlocks == 0) {
        /* Tier 3: no FDZ1 stream to honour, which covers numpress,
         * decrypted and cached channels. Cut fixed spectrum-count
         * batches instead. The units are arbitrary but still whole
         * spectra, so the visitor contract does not change. */
        const NSUInteger batch = 512;
        NSMutableArray<NSValue *> *out = [NSMutableArray array];
        for (NSUInteger i = from; i < hi; i += batch) {
            NSUInteger last = MIN(i + batch, hi) - 1;
            TTIOSpectralUnit u;
            u.block = out.count;
            u.firstSpectrum = i;
            u.nSpectra = last - i + 1;
            u.valueStart = off[i];
            u.valueEnd = off[last] + (unsigned long long)len[last];
            [out addObject:[NSValue valueWithBytes:&u objCType:@encode(TTIOSpectralUnit)]];
        }
        return out;
    }

    return [TTIOSpectralUnitPlan unitsForSpectraFrom:from to:hi
                                             offsets:off lengths:len
                                    blockValueStarts:bsD.mutableBytes count:nBlocks];
}

/* A view over one unit. Holds the unit's decoded columns and nothing
 * else; -spectrumAtIndex: slices one spectrum out of them on demand, so
 * a caller that reads part of a unit pays for that part. */
- (instancetype)_initAsUnitViewOfRun:(TTIOAcquisitionRun *)parent
                                unit:(TTIOSpectralUnit)u
                             columns:(NSDictionary<NSString *, NSData *> *)cols
{
    if ((self = [super init])) {
        _name                     = parent->_name;
        _acquisitionMode          = parent->_acquisitionMode;
        _instrumentConfig         = parent->_instrumentConfig;
        _spectrumIndex            = parent->_spectrumIndex;
        _channelNames             = parent->_channelNames;
        _spectrumClassName        = parent->_spectrumClassName;
        _modality                 = parent->_modality;
        _nucleusType              = parent->_nucleusType;
        _spectrometerFrequencyMHz = parent->_spectrometerFrequencyMHz;
        _solvent                  = parent->_solvent;
        _unitColumns              = cols;
        _unitFirstSpectrum        = u.firstSpectrum;
        _unitCount                = u.nSpectra;
        _unitValueStart           = u.valueStart;
    }
    return self;
}

/* View-local index in, run-global metadata out: spectrum k of the view
 * is run spectrum _unitFirstSpectrum + k, and its msLevel, polarity and
 * retention time come from the run's spectrum index at that global
 * position. */
- (NSData *)unitColumnForChannel:(NSString *)channel
                      valueStart:(unsigned long long *)outValueStart
{
    if (!_unitColumns) return nil;
    NSData *col = _unitColumns[channel];
    if (!col) return nil;
    if (outValueStart) *outValueStart = _unitValueStart;
    return col;
}

- (id)_unitSpectrumAtIndex:(NSUInteger)index error:(NSError **)error
{
    if (index >= _unitCount) {
        if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
            @"index %lu beyond unit spectrum count %lu",
            (unsigned long)index, (unsigned long)_unitCount);
        return nil;
    }
    NSUInteger i = _unitFirstSpectrum + index;
    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
    NSUInteger off = (NSUInteger)((unsigned long long)[_spectrumIndex offsetAt:i]
                                  - _unitValueStart);
    NSUInteger len = [_spectrumIndex lengthAt:i];
    NSMutableDictionary<NSString *, TTIOSignalArray *> *arrays =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
    for (NSString *ch in _channelNames) {
        NSData *col = _unitColumns[ch];
        if (!col || (off + len) * sizeof(double) > col.length) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"unit column '%@' is short for spectrum %lu", ch, (unsigned long)i);
            return nil;
        }
        NSData *d = [NSData dataWithBytes:(const uint8_t *)col.bytes + off * sizeof(double)
                                   length:len * sizeof(double)];
        arrays[ch] = [[TTIOSignalArray alloc] initWithOwnedBuffer:d length:len
                                                         encoding:enc axis:nil];
    }
    return [self _spectrumAtIndex:i channels:arrays error:error];
}

/* A sibling run holding just this unit's spectra.
 *
 * Each channel is read once for the unit's whole value extent and then
 * sliced, the way -iterSpectraWithBatch: reads a batch. Calling
 * -spectrumAtIndex: per spectrum instead costs one storage read per
 * spectrum per channel and measured 148 MB/s against the ordered
 * reader's 213 on the Exploris corpus. */
- (TTIOAcquisitionRun *)_viewForUnit:(TTIOSpectralUnit)u error:(NSError **)error
{
    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
    NSUInteger start = (NSUInteger)u.valueStart;
    NSUInteger count = (NSUInteger)(u.valueEnd - u.valueStart);

    NSMutableDictionary<NSString *, NSData *> *cols =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
    for (NSString *ch in _channelNames) {
        NSData *d = [self channelRange:ch offset:start count:count error:error];
        if (!d) return nil;
        cols[ch] = d;
    }

    return [[TTIOAcquisitionRun alloc] _initAsUnitViewOfRun:self
                                                       unit:u
                                                    columns:cols];
}

/* Units in flight is a memory setting before it is a concurrency one:
 * each stays resident for as long as the caller is inside it. The
 * resolver is the writer's, so TTIO_MEMORY_BUDGET means the same thing
 * on both sides.
 *
 * Sized from the WIDEST unit in the range, not the mean. Sizing from
 * the mean makes the budget non-binding exactly when it matters. */
- (NSUInteger)_unitWindowForThreads:(NSUInteger)nthreads
                              units:(NSArray<NSValue *> *)units
{
    if (units.count == 0) return 1;

    unsigned long long widest = 0;
    for (NSValue *v in units) {
        TTIOSpectralUnit u; [v getValue:&u];
        unsigned long long vals = u.valueEnd - u.valueStart;
        if (vals > widest) widest = vals;
    }
    if (widest == 0) return 1;

    /* A unit in flight costs its decoded values across every channel. */
    NSUInteger nch = MAX((NSUInteger)1, _channelNames.count);
    unsigned long long unitBytes = widest * 8ull * (unsigned long long)nch;

    unsigned long long budget = [TTIOThreads resolveMemoryBudget:nil
                                                         threads:nthreads
                                                      blockBytes:unitBytes];
    unsigned long long admits = budget / unitBytes;
    if (admits < 1) admits = 1;

    NSUInteger window = (NSUInteger)MIN((unsigned long long)nthreads, admits);
    /* More workers than units buys nothing: one unit, one thread. */
    if (window > units.count) window = units.count;
    return window < 1 ? 1 : window;
}

/* Per-channel FDZ1 block tables, built here because building one is a
 * storage read. Returns nil when any channel has no FDZ1 stream, which
 * sends the caller to the serial path. */
- (NSDictionary<NSString *, TTIOFDZBlockTable *> *)_fdzTablesForAllChannels
{
    if (_inMemorySpectra || !_storageDatasets) return nil;
    NSMutableDictionary<NSString *, TTIOFDZBlockTable *> *out =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
    for (NSString *ch in _channelNames) {
        if (![_fdzChannels containsObject:ch]) return nil;
        if (_decryptedChannels[ch] || _numpressChannels[ch] || _cachedFullChannels[ch]) {
            return nil;
        }
        id<TTIOStorageDataset> ds = _storageDatasets[ch];
        if (!ds) return nil;
        TTIOFDZByteRangeReader reader = ^NSData *(NSUInteger offset, NSUInteger n) {
            id raw = [ds readSliceAtOffset:offset count:n error:NULL];
            return [raw isKindOfClass:[NSData class]] ? raw : nil;
        };
        TTIOFDZBlockTable *t;
        @synchronized (self) {
            if (!_fdzTables) _fdzTables = [NSMutableDictionary dictionary];
            t = _fdzTables[ch];
            if (!t) {
                t = [TTIOFloatDeltaZstd readBlockTableWithReader:reader error:NULL];
                if (t) _fdzTables[ch] = t;
            }
        }
        if (!t) return nil;
        out[ch] = t;
    }
    return out;
}

/* Raw bytes of every FDZ1 block the unit's value extent touches, keyed
 * by channel then by the block's byte offset.
 *
 * ⛔ Storage reads. HDF5 is not thread-safe, so this runs on the calling
 * thread and only the decode that consumes the result goes to the pool.
 * A unit's extent can run past its own block's end, so this reads every
 * block the extent spans, not just the owning one. */
- (NSDictionary *)_rawBlocksForUnit:(TTIOSpectralUnit)u
                             tables:(NSDictionary<NSString *, TTIOFDZBlockTable *> *)tables
{
    NSMutableDictionary *byChannel = [NSMutableDictionary dictionaryWithCapacity:tables.count];
    for (NSString *ch in _channelNames) {
        TTIOFDZBlockTable *t = tables[ch];
        id<TTIOStorageDataset> ds = _storageDatasets[ch];
        if (!t || !ds) return nil;
        NSUInteger bs = t.blockSize;
        if (bs == 0) return nil;
        NSUInteger k0 = (NSUInteger)(u.valueStart / bs);
        NSUInteger k1 = (NSUInteger)((u.valueEnd - 1) / bs);
        NSMutableDictionary<NSNumber *, NSData *> *blocks = [NSMutableDictionary dictionary];
        for (NSUInteger k = k0; k <= k1; k++) {
            NSUInteger off = (NSUInteger)[t offsetAt:k];
            NSUInteger len = (NSUInteger)[t lengthAt:k];
            id raw = [ds readSliceAtOffset:off count:len error:NULL];
            if (![raw isKindOfClass:[NSData class]]) return nil;
            blocks[@(off)] = raw;
        }
        byChannel[ch] = blocks;
    }
    return byChannel;
}

/* Decode + materialise. Runs on the pool: every byte it needs is already
 * in `raw`, so it performs no storage read. */
- (TTIOAcquisitionRun *)_viewForUnit:(TTIOSpectralUnit)u
                                 raw:(NSDictionary *)raw
                              tables:(NSDictionary<NSString *, TTIOFDZBlockTable *> *)tables
                               error:(NSError **)error
{
    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
    NSUInteger nVals = (NSUInteger)(u.valueEnd - u.valueStart);
    NSMutableDictionary<NSString *, NSData *> *cols =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];

    for (NSString *ch in _channelNames) {
        TTIOFDZBlockTable *t = tables[ch];
        NSDictionary<NSNumber *, NSData *> *blocks = raw[ch];
        if (!t || !blocks) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"channel '%@' has no pre-read bytes for the unit", ch);
            return nil;
        }
        /* Serves any sub-range of a block already in memory. */
        TTIOFDZByteRangeReader memReader = ^NSData *(NSUInteger offset, NSUInteger n) {
            for (NSNumber *baseNum in blocks) {
                NSUInteger base = baseNum.unsignedIntegerValue;
                NSData *d = blocks[baseNum];
                if (offset >= base && offset + n <= base + d.length) {
                    return [d subdataWithRange:NSMakeRange(offset - base, n)];
                }
            }
            return nil;
        };
        NSUInteger bs = t.blockSize;
        NSMutableData *acc = [NSMutableData dataWithCapacity:nVals * sizeof(double)];
        unsigned long long pos = u.valueStart;
        while (pos < u.valueEnd) {
            NSUInteger k = (NSUInteger)(pos / bs);
            NSData *blk = [TTIOFloatDeltaZstd decodeBlock:k table:t reader:memReader error:error];
            if (!blk) return nil;
            unsigned long long blkStart = (unsigned long long)k * bs;
            NSUInteger inBlk = (NSUInteger)(pos - blkStart);
            NSUInteger avail = (NSUInteger)(blk.length / sizeof(double)) - inBlk;
            NSUInteger want = (NSUInteger)MIN((unsigned long long)avail, u.valueEnd - pos);
            [acc appendBytes:(const uint8_t *)blk.bytes + inBlk * sizeof(double)
                      length:want * sizeof(double)];
            pos += want;
        }
        cols[ch] = acc;
    }

    return [[TTIOAcquisitionRun alloc] _initAsUnitViewOfRun:self
                                                       unit:u
                                                    columns:cols];
}

- (BOOL)iterBlocksFrom:(NSUInteger)from
                    to:(NSUInteger)to
               threads:(NSUInteger)threads
                 error:(NSError **)error
            usingBlock:(void (^)(TTIOAcquisitionRun *view, NSUInteger viewStart,
                                 NSUInteger firstSpectrum, NSUInteger nSpectra,
                                 BOOL *stop))userBlock
{
    NSArray<NSValue *> *units = [self _unitsFrom:from to:to];
    if (units.count == 0) return YES;
    NSUInteger nthreads = threads ? threads : [TTIOThreads resolve:nil];

    NSDictionary<NSString *, TTIOFDZBlockTable *> *tables =
        (nthreads > 1 && units.count > 1) ? [self _fdzTablesForAllChannels] : nil;

    /* One unit, one thread, or a channel the pool path cannot serve:
     * deliver on this thread. */
    if (!tables) {
        BOOL halted = NO;
        for (NSValue *v in units) {
            if (halted) break;
            TTIOSpectralUnit u; [v getValue:&u];
            NSError *ie = nil;
            TTIOAcquisitionRun *view = [self _viewForUnit:u error:&ie];
            if (!view) { if (error) *error = ie; return NO; }
            userBlock(view, 0, u.firstSpectrum, u.nSpectra, &halted);
        }
        return YES;
    }

    NSUInteger window = [self _unitWindowForThreads:nthreads units:units];
    TTIOThreadPool *pool = [TTIOThreadPool poolWithThreads:window];
    if (pool.queue == nil) {
        return [self iterBlocksFrom:from to:to threads:1 error:error usingBlock:userBlock];
    }

    NSMutableDictionary<NSNumber *, TTIOInFlightUnit *> *pending =
        [NSMutableDictionary dictionary];
    NSCondition *cond = [NSCondition new];
    __block BOOL halted = NO;
    BOOL ok = YES;
    __weak typeof(self) weakSelf = self;

    void (^submit)(NSUInteger) = ^(NSUInteger idx) {
        typeof(self) sself = weakSelf;
        if (!sself || idx >= units.count || pending[@(idx)] || halted) return;
        TTIOSpectralUnit u; [units[idx] getValue:&u];
        TTIOInFlightUnit *f = [TTIOInFlightUnit new];
        pending[@(idx)] = f;
        /* Storage reads, this thread. */
        NSDictionary *raw = [sself _rawBlocksForUnit:u tables:tables];
        if (!raw) {
            f.error = TTIOMakeError(TTIOErrorDatasetRead,
                @"unit %lu could not be read", (unsigned long)idx);
            f.done = YES;
            return;
        }
        [pool.queue addOperationWithBlock:^{
            NSError *ie = nil;
            @try {
                TTIOAcquisitionRun *view = [sself _viewForUnit:u raw:raw tables:tables
                                                          error:&ie];
                if (view && !halted) {
                    BOOL st = NO;
                    userBlock(view, 0, u.firstSpectrum, u.nSpectra, &st);
                    if (st) halted = YES;
                } else if (!view && !ie) {
                    ie = TTIOMakeError(TTIOErrorDatasetRead,
                        @"unit %lu did not materialise", (unsigned long)idx);
                }
            } @catch (NSException *ex) {
                ie = TTIOMakeError(TTIOErrorDatasetRead, @"%@: %@", ex.name, ex.reason);
            }
            [cond lock];
            f.error = ie; f.done = YES;
            [cond broadcast];
            [cond unlock];
        }];
    };

    for (NSUInteger i = 0; i < MIN(units.count, window); i++) submit(i);
    for (NSUInteger i = 0; i < units.count && !halted; i++) {
        TTIOInFlightUnit *f = pending[@(i)];
        if (!f) continue;
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
        [pending removeObjectForKey:@(i)];
        if (f.error) { if (error) *error = f.error; ok = NO; break; }
        submit(i + window);
    }
    for (TTIOInFlightUnit *f in pending.allValues) {
        [cond lock];
        while (!f.done) [cond wait];
        [cond unlock];
    }
    [pool close];
    return ok;
}

- (BOOL)iterSpectraWithBatch:(NSUInteger)batch
                     threads:(NSUInteger)threads
                       error:(NSError **)error
                  usingBlock:(void (^)(id spectrum, NSUInteger index, BOOL *stop))block
{
    NSUInteger keep = _iterThreads;
    _iterThreads = threads ? threads : [TTIOThreads resolve:nil];
    BOOL ok = [self iterSpectraWithBatch:batch error:error usingBlock:block];
    _iterThreads = keep;
    return ok;
}

- (BOOL)iterSpectraWithBatch:(NSUInteger)batch
                       error:(NSError **)error
                  usingBlock:(void (^)(id spectrum, NSUInteger index, BOOL *stop))block
{
    NSUInteger n = [self count];
    if (batch < 1) batch = 1;
    BOOL halted = NO;
    if (_inMemorySpectra) {
        for (NSUInteger i = 0; i < n && !halted; i++) block(_inMemorySpectra[i], i, &halted);
        return YES;
    }
    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];
    for (NSUInteger b0 = 0; b0 < n && !halted; b0 += batch) {
        NSUInteger b1 = MIN(n, b0 + batch);
        NSUInteger start = (NSUInteger)[_spectrumIndex offsetAt:b0];
        NSUInteger end = (NSUInteger)([_spectrumIndex offsetAt:b1 - 1] + [_spectrumIndex lengthAt:b1 - 1]);
        NSMutableDictionary<NSString *, NSData *> *cols = [NSMutableDictionary dictionary];
        for (NSString *c in _channelNames) {
            NSData *d = _iterThreads > 1
                ? [self channelRange:c start:start count:end - start threads:_iterThreads error:error]
                : [self channelRange:c offset:start count:end - start error:error];
            if (!d) return NO;
            cols[c] = d;
        }
        for (NSUInteger i = b0; i < b1 && !halted; i++) {
            NSUInteger off = (NSUInteger)[_spectrumIndex offsetAt:i] - start;
            NSUInteger len = [_spectrumIndex lengthAt:i];
            NSMutableDictionary *arrays = [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
            for (NSString *c in _channelNames) {
                NSData *d = [NSData dataWithBytes:(const uint8_t *)cols[c].bytes + off * sizeof(double)
                                           length:len * sizeof(double)];
                arrays[c] = [[TTIOSignalArray alloc] initWithOwnedBuffer:d length:len encoding:enc axis:nil];
            }
            id sp = [self _spectrumAtIndex:i channels:arrays error:error];
            if (!sp) return NO;
            block(sp, i, &halted);
        }
    }
    return YES;
}

- (id)spectrumAtIndex:(NSUInteger)index error:(NSError **)error
{
    if (_unitColumns) return [self _unitSpectrumAtIndex:index error:error];

    if (_inMemorySpectra) {
        if (index >= _inMemorySpectra.count) {
            if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
                @"index %lu beyond spectrum count %lu",
                (unsigned long)index, (unsigned long)_inMemorySpectra.count);
            return nil;
        }
        return _inMemorySpectra[index];
    }

    if (index >= _spectrumIndex.count) {
        if (error) *error = TTIOMakeError(TTIOErrorOutOfRange,
            @"index %lu beyond spectrum count %lu",
            (unsigned long)index, (unsigned long)_spectrumIndex.count);
        return nil;
    }

    uint64_t off = [_spectrumIndex offsetAt:index];
    uint32_t len = [_spectrumIndex lengthAt:index];

    TTIOEncodingSpec *enc =
        [TTIOEncodingSpec specWithPrecision:TTIOPrecisionFloat64
                       compressionAlgorithm:TTIOCompressionZlib
                                  byteOrder:TTIOByteOrderLittleEndian];

    NSMutableDictionary<NSString *, TTIOSignalArray *> *channels =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
    for (NSString *chName in _channelNames) {
        NSData *d = [self channelRange:chName offset:(NSUInteger)off count:len error:error];
        if (!d) return nil;
        TTIOSignalArray *sa = [[TTIOSignalArray alloc] initWithOwnedBuffer:d
                                                                    length:len
                                                                  encoding:enc
                                                                      axis:nil];
        channels[chName] = sa;
    }
    return [self _spectrumAtIndex:index channels:channels error:error];
}

// Build the spectrum object of the persisted class from its channel
// arrays and the index row.
- (id)_spectrumAtIndex:(NSUInteger)index
              channels:(NSDictionary<NSString *, TTIOSignalArray *> *)channels
                 error:(NSError **)error
{
    // _spectrumClassName remains the persisted source of truth; the enum
    // is an in-code dispatch key only (P3.8). An unrecognised class falls
    // through to the unknown-class error below, exactly as before.
    TTIOSpectrumKind k = TTIOSpectrumKindFromPersisted(_spectrumClassName);
    if (k == TTIOSpectrumKindMass) {
        TTIOMassSpectrum *ms = [[TTIOMassSpectrum alloc]
                initWithMzArray:channels[@"mz"]
                 intensityArray:channels[@"intensity"]
                        msLevel:[_spectrumIndex msLevelAt:index]
                       polarity:[_spectrumIndex polarityAt:index]
                     scanWindow:nil
                  indexPosition:index
                scanTimeSeconds:[_spectrumIndex retentionTimeAt:index]
                    precursorMz:[_spectrumIndex precursorMzAt:index]
                precursorCharge:[_spectrumIndex precursorChargeAt:index]
                          error:error];
        if (ms) ms.isCentroided = [_spectrumIndex centroidedAt:index];
        return ms;
    }

    if (k == TTIOSpectrumKindNMR) {
        return [[TTIONMRSpectrum alloc]
                initWithChemicalShiftArray:channels[@"chemical_shift"]
                            intensityArray:channels[@"intensity"]
                               nucleusType:_nucleusType
                  spectrometerFrequencyMHz:_spectrometerFrequencyMHz
                             indexPosition:index
                           scanTimeSeconds:[_spectrumIndex retentionTimeAt:index]
                                     error:error];
    }

    // Vibrational types dispatch on the stored spectrum_class (parity
    // with Python / Java); channels are wavenumber/intensity (IR/Raman)
    // or wavelength/absorbance (UV-Vis).
    if (k == TTIOSpectrumKindIR) {
        return [[TTIOIRSpectrum alloc]
                initWithWavenumberArray:channels[@"wavenumber"]
                         intensityArray:channels[@"intensity"]
                                   mode:_irMode
                        resolutionCmInv:_irResolutionCmInv
                          numberOfScans:_irNumberOfScans
                          indexPosition:index
                        scanTimeSeconds:[_spectrumIndex retentionTimeAt:index]
                                  error:error];
    }

    if (k == TTIOSpectrumKindRaman) {
        return [[TTIORamanSpectrum alloc]
                initWithWavenumberArray:channels[@"wavenumber"]
                         intensityArray:channels[@"intensity"]
                 excitationWavelengthNm:_ramanExcitationWavelengthNm
                           laserPowerMw:_ramanLaserPowerMw
                     integrationTimeSec:_ramanIntegrationTimeSec
                          indexPosition:index
                        scanTimeSeconds:[_spectrumIndex retentionTimeAt:index]
                                  error:error];
    }

    if (k == TTIOSpectrumKindUVVis) {
        return [[TTIOUVVisSpectrum alloc]
                initWithWavelengthArray:channels[@"wavelength"]
                        absorbanceArray:channels[@"absorbance"]
                           pathLengthCm:_uvvisPathLengthCm
                                solvent:_solvent
                          indexPosition:index
                        scanTimeSeconds:[_spectrumIndex retentionTimeAt:index]
                                  error:error];
    }

    if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
        @"unknown spectrum_class %@ in acquisition run", _spectrumClassName);
    return nil;
}

- (NSArray<NSNumber *> *)indicesInRetentionTimeRange:(TTIOValueRange *)range
{
    NSIndexSet *set = [_spectrumIndex indicesInRetentionTimeRange:range];
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:set.count];
    NSUInteger idx = [set firstIndex];
    while (idx != NSNotFound) {
        [out addObject:@(idx)];
        idx = [set indexGreaterThanIndex:idx];
    }
    return out;
}

#pragma mark - TTIOIndexable

- (id)objectAtIndex:(NSUInteger)index
{
    return [self spectrumAtIndex:index error:NULL];
}

- (NSUInteger)count
{
    if (_unitColumns) return _unitCount;
    return _inMemorySpectra ? _inMemorySpectra.count : _spectrumIndex.count;
}

- (NSArray *)objectsInRange:(NSRange)range
{
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:range.length];
    for (NSUInteger i = range.location; i < NSMaxRange(range); i++) {
        id obj = [self objectAtIndex:i];
        if (obj) [out addObject:obj];
    }
    return out;
}

- (NSArray *)spectra
{
    NSUInteger n = [self count];
    if (!_inMemorySpectra) {
        // Whole columns once, then per-spectrum slices from memory.
        for (NSString *c in _channelNames) (void)[self _fullChannel:c error:NULL];
    }
    NSMutableArray *out = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        id obj = [self spectrumAtIndex:i error:NULL];
        if (obj) [out addObject:obj];
    }
    return [out copy];
}

#pragma mark - TTIOStreamable

- (id)nextObject
{
    if (![self hasMore]) return nil;
    id obj = [self objectAtIndex:_streamPosition];
    _streamPosition++;
    return obj;
}

- (BOOL)hasMore               { return _streamPosition < [self count]; }
- (NSUInteger)currentPosition { return _streamPosition; }
- (BOOL)seekToPosition:(NSUInteger)position
{
    if (position > [self count]) return NO;
    _streamPosition = position;
    return YES;
}
- (void)reset                 { _streamPosition = 0; }

#pragma mark - Persistence context

- (void)setPersistenceFilePath:(NSString *)path runName:(NSString *)runName
{
    _persistenceFilePath = [path copy];
    _persistenceRunName  = [runName copy];
    // Phase 1: keep the public Run-protocol name in sync with the
    // persistence-context run name when the latter is supplied.
    if (runName.length > 0) {
        _name = [runName copy];
    }
}

- (void)releaseHDF5Handles
{
    _storageDatasets     = nil;
    _storageSignalGroup  = nil;
    // Fix #1: drop the whole-column cache when the backing handles go
    // away so it can never serve bytes from a closed/replaced dataset.
    @synchronized (self) {
        _cachedFullChannels = nil;
        _fdzTables = nil;
        _fdzBlockCache = nil;
        _fdzBlockCacheIndex = nil;
    }
}

#pragma mark - TTIOProvenanceable

- (void)addProcessingStep:(TTIOProvenanceRecord *)step
{
    if (!_provenance) _provenance = [NSMutableArray array];
    if (step) [_provenance addObject:step];
}

- (NSArray<TTIOProvenanceRecord *> *)provenanceChain
{
    return _provenance ? [_provenance copy] : @[];
}

- (NSArray<NSString *> *)inputEntities
{
    NSMutableSet *set = [NSMutableSet set];
    for (TTIOProvenanceRecord *r in _provenance) [set addObjectsFromArray:r.inputRefs];
    return [set allObjects];
}

- (NSArray<NSString *> *)outputEntities
{
    NSMutableSet *set = [NSMutableSet set];
    for (TTIOProvenanceRecord *r in _provenance) [set addObjectsFromArray:r.outputRefs];
    return [set allObjects];
}

#pragma mark - TTIOEncryptable

- (BOOL)encryptWithKey:(NSData *)key
                 level:(TTIOEncryptionLevel)level
                 error:(NSError **)error
{
    (void)level;
    if (!_persistenceFilePath || !_persistenceRunName) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"TTIOAcquisitionRun: cannot encrypt in-memory run; persist via "
            @"TTIOSpectralDataset first so the run has a file context");
        return NO;
    }
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    return [TTIOEncryptionManager encryptIntensityChannelInRun:_persistenceRunName
                                                    atFilePath:_persistenceFilePath
                                                       withKey:key
                                                         error:error];
#pragma clang diagnostic pop
}

- (BOOL)decryptWithKey:(NSData *)key error:(NSError **)error
{
    if (!_persistenceFilePath || !_persistenceRunName) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"TTIOAcquisitionRun: no persistence context for decrypt");
        return NO;
    }
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    NSData *plain = [TTIOEncryptionManager
                      decryptIntensityChannelInRun:_persistenceRunName
                                        atFilePath:_persistenceFilePath
                                           withKey:key
                                             error:error];
#pragma clang diagnostic pop
    if (!plain) return NO;

    // M5-handoff: cache the concatenated plaintext so
    // -spectrumAtIndex: can slice it directly. The on-disk file stays
    // encrypted — only the open handle sees plaintext.
    if (!_decryptedChannels) {
        _decryptedChannels = [NSMutableDictionary dictionary];
    }
    _decryptedChannels[@"intensity"] = plain;
    return YES;
}

/** Expose the decrypted plaintext for channel ``chName`` if
 *  -decryptWithKey:error: has populated the cache. Returns nil
 *  otherwise. Consumed by TTIOSpectralDataset so the dataset-level
 *  -decryptWithKey:error: can return a {runName: plaintext} NSDictionary
 *  matching the Python surface. Internal API. */
- (NSData *)decryptedChannelNamed:(NSString *)chName
{
    return _decryptedChannels[chName];
}

/** Reattach storage handles after a dataset-level decrypt that had to
 *  close the file for compound-dataset unsealing. Accepts the fresh
 *  signal_channels ``TTIOHDF5Group`` from the reopened file and
 *  rebuilds ``_storageSignalGroup`` / ``_storageDatasets`` so
 *  ``spectrumAtIndex:`` can once again read unencrypted channels
 *  (mz, chemical_shift, ...) from disk. The decrypted intensity
 *  channel continues to serve from the in-memory cache. Internal API
 *  — called only by TTIOSpectralDataset. */
- (BOOL)reattachSignalHandlesFromGroup:(id<TTIOStorageGroup>)channels error:(NSError **)error
{
    if (!channels) return NO;
    NSMutableDictionary<NSString *, id<TTIOStorageDataset>> *datasets =
        [NSMutableDictionary dictionaryWithCapacity:_channelNames.count];
    for (NSString *chName in _channelNames) {
        NSString *dsName = [chName stringByAppendingString:@"_values"];
        if (![channels hasChildNamed:dsName]) continue;  // encrypted / absent
        id<TTIOStorageDataset> ds = [channels openDatasetNamed:dsName error:error];
        if (!ds) return NO;
        datasets[chName] = ds;
    }
    _storageSignalGroup = channels;
    _storageDatasets    = datasets;
    // Fix #1: the datasets were just rebuilt from a reopened file; any
    // previously cached columns belonged to the old handles. Invalidate
    // so the next read repopulates from the fresh datasets.
    @synchronized (self) {
        _cachedFullChannels = nil;
        _fdzTables = nil;
        _fdzBlockCache = nil;
        _fdzBlockCacheIndex = nil;
    }
    return YES;
}

- (TTIOAccessPolicy *)accessPolicy         { return _accessPolicy; }
- (void)setAccessPolicy:(TTIOAccessPolicy *)policy { _accessPolicy = policy; }

@end
