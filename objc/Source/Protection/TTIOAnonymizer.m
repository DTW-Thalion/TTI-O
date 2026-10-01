/*
 * TTIOAnonymizer.m
 * TTI-O Objective-C Implementation
 *
 * Classes:       TTIOAnonymizationPolicy, TTIOAnonymizationResult,
 *                TTIOAnonymizer
 * Declared In:   Protection/TTIOAnonymizer.h
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOAnonymizer.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Dataset/TTIOIdentification.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Dataset/TTIOWrittenRun.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOInstrumentConfig.h"
#import "Run/TTIOSpectrumIndex.h"
#import "Spectra/TTIOMassSpectrum.h"
#import "Spectra/TTIONMRSpectrum.h"
#import "Spectra/TTIOChromatogram.h"
#import "Core/TTIOSignalArray.h"
#import "ValueClasses/TTIOEncodingSpec.h"
#import "ValueClasses/TTIOEnums.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOWrittenGenomicRun.h"

// Forward-declare the private helper for building transformed
// genomic runs (without writing). The unified write path threads
// the resulting dict straight into ``writeMinimalToPath:msRuns:
// genomicRuns:...`` alongside the MS dict — no per-modality
// fallback writers.
@interface TTIOAnonymizer ()
+ (NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)
    _buildTransformedGenomicRuns:(TTIOSpectralDataset *)source
                            policy:(TTIOAnonymizationPolicy *)policy
                            result:(TTIOAnonymizationResult *)result
                       appliedList:(NSMutableArray<NSString *> *)applied
                             error:(NSError **)error;
@end

#pragma mark - Policy

@implementation TTIOAnonymizationPolicy

- (instancetype)init
{
    self = [super init];
    if (self) {
        _rareMetaboliteThreshold = 0.05;
        _coarsenMzDecimals = -1;
        _coarsenChemicalShiftDecimals = -1;
        // M90.3 defaults
        _randomiseQualitiesConstant = 30;
    }
    return self;
}

@end

#pragma mark - Result

@implementation TTIOAnonymizationResult
@end

#pragma mark - Helpers

static BOOL isSAAV(TTIOIdentification *ident)
{
    NSString *upper = [ident.chemicalEntity uppercaseString];
    return ([upper rangeOfString:@"SAAV"].location != NSNotFound ||
            [upper rangeOfString:@"VARIANT"].location != NSNotFound);
}

static double *copyDoubleArray(TTIOSignalArray *arr)
{
    NSUInteger n = arr.length;
    double *out = (double *)malloc(n * sizeof(double));
    memcpy(out, arr.buffer.bytes, n * sizeof(double));
    return out;
}

static void roundArray(double *arr, NSUInteger n, NSInteger decimals)
{
    double scale = pow(10.0, (double)decimals);
    for (NSUInteger i = 0; i < n; i++) {
        arr[i] = round(arr[i] * scale) / scale;
    }
}

// M90.13 — CIGAR reference span. Returns the number of reference
// bases consumed by the alignment described by ``cigar``. Ops that
// consume reference: M, D, N, =, X. Ops that don't: I, S, H, P. A
// return of 0 means "unknown" (empty / "*" / non-parseable) and the
// caller falls back to the M90.3 position-only check.
static int64_t cigarRefSpan(NSString *cigar)
{
    if (cigar == nil || cigar.length == 0) return 0;
    if ([cigar isEqualToString:@"*"]) return 0;
    int64_t total = 0;
    int64_t accum = 0;
    BOOL haveDigits = NO;
    NSUInteger n = cigar.length;
    for (NSUInteger i = 0; i < n; i++) {
        unichar ch = [cigar characterAtIndex:i];
        if (ch >= '0' && ch <= '9') {
            accum = accum * 10 + (int64_t)(ch - '0');
            haveDigits = YES;
        } else {
            if (!haveDigits) return 0;  // op without count = malformed
            if (ch == 'M' || ch == 'D' || ch == 'N'
                || ch == '=' || ch == 'X') {
                total += accum;
            } else if (ch == 'I' || ch == 'S' || ch == 'H' || ch == 'P') {
                // Don't consume reference bases.
            } else {
                return 0;  // unknown op — bail out
            }
            accum = 0;
            haveDigits = NO;
        }
    }
    if (haveDigits) return 0;  // trailing digits with no op = malformed
    return total;
}

// M90.14 — xoshiro256** PRNG. Self-contained; seeded by splitmix64
// from a single 64-bit seed. Reproducible within ObjC: same seed →
// same byte sequence. NOT byte-equal to numpy's PCG64 (cross-
// language byte parity is not a goal of M90.14 — that's a follow-up).
typedef struct {
    uint64_t s[4];
} TTIORngState;

static uint64_t splitmix64Next(uint64_t *state)
{
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void ttioRngSeed(TTIORngState *r, uint64_t seed)
{
    uint64_t sm = seed;
    r->s[0] = splitmix64Next(&sm);
    r->s[1] = splitmix64Next(&sm);
    r->s[2] = splitmix64Next(&sm);
    r->s[3] = splitmix64Next(&sm);
}

static uint64_t ttioRngNextU64(TTIORngState *r)
{
    // xoshiro256**
    uint64_t result = r->s[1] * 5;
    result = ((result << 7) | (result >> (64 - 7))) * 9;
    uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = (r->s[3] << 45) | (r->s[3] >> (64 - 45));
    return result;
}

// Draw a uniform integer in [0, bound) via Lemire's debiased method.
// Bound must be > 0.
static uint32_t ttioRngNextBoundedU32(TTIORngState *r, uint32_t bound)
{
    uint64_t x = (uint32_t)(ttioRngNextU64(r) >> 32);
    uint64_t m = x * (uint64_t)bound;
    uint32_t l = (uint32_t)m;
    if (l < bound) {
        uint32_t t = (uint32_t)(-(int32_t)bound) % bound;
        while (l < t) {
            x = (uint32_t)(ttioRngNextU64(r) >> 32);
            m = x * (uint64_t)bound;
            l = (uint32_t)m;
        }
    }
    return (uint32_t)(m >> 32);
}

#pragma mark - Anonymizer

@implementation TTIOAnonymizer

+ (TTIOAnonymizationResult *)anonymizeDataset:(TTIOSpectralDataset *)source
                                   outputPath:(NSString *)outputPath
                                       policy:(TTIOAnonymizationPolicy *)policy
                                        error:(NSError **)error
{
    TTIOAnonymizationResult *result = [[TTIOAnonymizationResult alloc] init];
    NSMutableArray<NSString *> *appliedPolicies = [NSMutableArray array];

    NSArray<TTIOIdentification *> *identifications = source.identifications;

    // Build SAAV index set
    NSMutableSet<NSString *> *saavKeys = [NSMutableSet set];
    if (policy.redactSAAVSpectra) {
        for (TTIOIdentification *ident in identifications) {
            if (isSAAV(ident)) {
                NSString *key = [NSString stringWithFormat:@"%@:%lu",
                                 ident.runName, (unsigned long)ident.spectrumIndex];
                [saavKeys addObject:key];
            }
        }
    }

    // Rare metabolite lookup
    NSDictionary<NSString *, NSNumber *> *prevalence = policy.rareMetaboliteTable;

    // Process each run. Build TTIOWrittenRun objects in-flight rather
    // than TTIOAcquisitionRun + TTIOMassSpectrum object graphs: that
    // lets the mixed MS+genomic write path below feed a single
    // ``writeMinimalToPath:msRuns:genomicRuns:...`` call instead of
    // the legacy ``writeToFilePath:`` + open-RW + append-genomic
    // dance. Mirrors Python's ``anonymization.py`` flow which builds
    // ``new_runs: dict[str, WrittenRun]`` and ``new_genomic_runs``
    // up front and hands them to ``SpectralDataset.write_minimal`` in
    // one shot.
    NSMutableDictionary<NSString *, TTIOWrittenRun *> *newRuns =
        [NSMutableDictionary dictionary];
    NSArray<NSString *> *runNames =
        [[source.msRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];

    for (NSString *runName in runNames) {
        TTIOAcquisitionRun *run = source.msRuns[runName];
        NSUInteger nSpectra = run.spectrumIndex.count;

        // Concatenated channel buffers + parallel index columns.
        NSMutableData *mzBufConcat = [NSMutableData data];
        NSMutableData *intBufConcat = [NSMutableData data];
        NSMutableData *offsets = [NSMutableData data];
        NSMutableData *lengths = [NSMutableData data];
        NSMutableData *rts = [NSMutableData data];
        NSMutableData *mlsBuf = [NSMutableData data];
        NSMutableData *polsBuf = [NSMutableData data];
        NSMutableData *pmzsBuf = [NSMutableData data];
        NSMutableData *pcsBuf = [NSMutableData data];
        NSMutableData *bpisBuf = [NSMutableData data];
        // M102: imaging pixels keep their grid positions.
        TTIOSpectrumIndex *srcIdx = run.spectrumIndex;
        BOOL keepPixels = srcIdx.hasPixelCoordinates;
        NSMutableData *pxBuf = keepPixels ? [NSMutableData data] : nil;
        NSMutableData *pyBuf = keepPixels ? [NSMutableData data] : nil;
        NSMutableData *pzBuf = keepPixels ? [NSMutableData data] : nil;
        int64_t cursor = 0;
        NSUInteger keptCount = 0;

        for (NSUInteger i = 0; i < nSpectra; i++) {
            // SAAV redaction
            if (policy.redactSAAVSpectra) {
                NSString *key = [NSString stringWithFormat:@"%@:%lu",
                                 runName, (unsigned long)i];
                if ([saavKeys containsObject:key]) {
                    result.spectraRedacted++;
                    continue;
                }
            }

            NSError *specErr = nil;
            id specObj = [run spectrumAtIndex:i error:&specErr];
            if (!specObj) continue;

            // The MS-only object path: NMR / non-mass spectra are
            // silently skipped here, matching the pre-refactor
            // behaviour. Extending the anonymiser to NMR is a separate
            // workstream tracked by the python reference impl's
            // chemical-shift coarsening tests.
            if (![specObj isKindOfClass:[TTIOMassSpectrum class]]) continue;

            TTIOMassSpectrum *ms = (TTIOMassSpectrum *)specObj;
            NSUInteger n = ms.mzArray.length;
            double *mzBuf = copyDoubleArray(ms.mzArray);
            double *intBuf = copyDoubleArray(ms.intensityArray);

            if (policy.coarsenMzDecimals >= 0) {
                roundArray(mzBuf, n, policy.coarsenMzDecimals);
                result.mzValuesCoarsened += n;
            }

            if (policy.maskIntensityBelowQuantile > 0.0) {
                double sorted[n];
                memcpy(sorted, intBuf, n * sizeof(double));
                for (NSUInteger a = 0; a < n; a++)
                    for (NSUInteger b = a + 1; b < n; b++)
                        if (sorted[a] > sorted[b]) {
                            double tmp = sorted[a];
                            sorted[a] = sorted[b];
                            sorted[b] = tmp;
                        }
                NSUInteger qIdx =
                    (NSUInteger)(policy.maskIntensityBelowQuantile * (double)(n - 1));
                double threshold = sorted[qIdx];
                for (NSUInteger j = 0; j < n; j++) {
                    if (intBuf[j] < threshold) {
                        intBuf[j] = 0.0;
                        result.intensitiesZeroed++;
                    }
                }
            }

            if (policy.maskRareMetabolites && prevalence) {
                for (TTIOIdentification *ident in identifications) {
                    if (![ident.runName isEqualToString:runName]) continue;
                    if (ident.spectrumIndex != i) continue;
                    NSNumber *prev = prevalence[ident.chemicalEntity];
                    if (prev && prev.doubleValue < policy.rareMetaboliteThreshold) {
                        memset(intBuf, 0, n * sizeof(double));
                        result.metabolitesMasked++;
                        break;
                    }
                }
            }

            // Recompute base peak intensity from the (possibly
            // masked) intensity buffer so anonymised files don't
            // carry stale base-peak metadata that contradicts the
            // signal channel. Mirrors Python which copies it from
            // the source index column unchanged; using the
            // post-mask max keeps the field consistent under
            // intensity-zero masking.
            double bpi = 0.0;
            for (NSUInteger j = 0; j < n; j++) {
                if (intBuf[j] > bpi) bpi = intBuf[j];
            }

            // Append this spectrum's transformed buffers to the
            // run-level concatenated channel data and parallel
            // per-spectrum metadata columns.
            [mzBufConcat appendBytes:mzBuf length:n * sizeof(double)];
            [intBufConcat appendBytes:intBuf length:n * sizeof(double)];
            free(mzBuf); free(intBuf);

            int64_t off = cursor;
            uint32_t len32 = (uint32_t)n;
            double rt = ms.scanTimeSeconds;
            int32_t mlv = (int32_t)ms.msLevel;
            int32_t pol = (int32_t)ms.polarity;
            double pmz = ms.precursorMz;
            int32_t pc = (int32_t)ms.precursorCharge;

            [offsets appendBytes:&off length:sizeof(off)];
            [lengths appendBytes:&len32 length:sizeof(len32)];
            [rts appendBytes:&rt length:sizeof(rt)];
            [mlsBuf appendBytes:&mlv length:sizeof(mlv)];
            [polsBuf appendBytes:&pol length:sizeof(pol)];
            [pmzsBuf appendBytes:&pmz length:sizeof(pmz)];
            [pcsBuf appendBytes:&pc length:sizeof(pc)];
            [bpisBuf appendBytes:&bpi length:sizeof(bpi)];
            if (keepPixels) {
                int32_t px = [srcIdx pixelXAt:i];
                int32_t py = [srcIdx pixelYAt:i];
                int32_t pz = [srcIdx pixelZAt:i];
                [pxBuf appendBytes:&px length:sizeof(px)];
                [pyBuf appendBytes:&py length:sizeof(py)];
                [pzBuf appendBytes:&pz length:sizeof(pz)];
            }

            cursor += (int64_t)n;
            keptCount++;
        }

        // Use the source run's spectrum class so NMR-shape runs
        // (although currently skipped above) at least carry the
        // right class string when downstream NMR support lands.
        NSString *spectrumClass = run.spectrumClassName ?: @"TTIOMassSpectrum";

        TTIOWrittenRun *newRun = [[TTIOWrittenRun alloc]
            initWithSpectrumClassName:spectrumClass
                      acquisitionMode:(int64_t)run.acquisitionMode
                          channelData:@{@"mz": mzBufConcat,
                                        @"intensity": intBufConcat}
                              offsets:offsets
                              lengths:lengths
                       retentionTimes:rts
                             msLevels:mlsBuf
                           polarities:polsBuf
                         precursorMzs:pmzsBuf
                     precursorCharges:pcsBuf
                  basePeakIntensities:bpisBuf];
        if (run.nucleusType.length > 0) {
            newRun.nucleusType = run.nucleusType;
        }
        if (keepPixels) {
            newRun.pixelX = pxBuf;
            newRun.pixelY = pyBuf;
            newRun.pixelZ = pzBuf;
        }
        // Carry per-run provenance forward so the anonymised file's
        // MS run keeps its history (now that TTIOWrittenRun has the
        // ``provenanceRecords`` field — the deferred parity gap was
        // closed in the same workstream as this refactor).
        NSArray<TTIOProvenanceRecord *> *srcChain = [run provenanceChain];
        if (srcChain.count > 0) {
            newRun.provenanceRecords = srcChain;
        }
        newRuns[runName] = newRun;
        (void)keptCount;
    }

    // Track which policies fired
    if (policy.redactSAAVSpectra && result.spectraRedacted > 0)
        [appliedPolicies addObject:@"redact_saav_spectra"];
    if (policy.coarsenMzDecimals >= 0 && result.mzValuesCoarsened > 0)
        [appliedPolicies addObject:@"coarsen_mz_decimals"];
    if (policy.maskIntensityBelowQuantile > 0.0 && result.intensitiesZeroed > 0)
        [appliedPolicies addObject:@"mask_intensity_below_quantile"];
    if (policy.maskRareMetabolites && result.metabolitesMasked > 0)
        [appliedPolicies addObject:@"mask_rare_metabolites"];

    NSString *title = source.title;
    if (policy.stripMetadataFields) {
        title = @"";
        result.metadataFieldsStripped = 1;
        [appliedPolicies addObject:@"strip_metadata_fields"];
    }

    result.policiesApplied = appliedPolicies;

    // Build provenance
    NSDictionary *params = @{
        @"policies": appliedPolicies,
        @"spectraRedacted": @(result.spectraRedacted),
        @"intensitiesZeroed": @(result.intensitiesZeroed),
        @"mzValuesCoarsened": @(result.mzValuesCoarsened),
        @"metabolitesMasked": @(result.metabolitesMasked),
        @"metadataFieldsStripped": @(result.metadataFieldsStripped),
    };
    TTIOProvenanceRecord *prov =
        [[TTIOProvenanceRecord alloc]
            initWithInputRefs:@[source.filePath ?: @""]
                     software:@"ttio anonymizer v0.4"
                   parameters:params
                   outputRefs:@[outputPath]
                timestampUnix:(int64_t)time(NULL)];

    // Single unified write path: build the transformed genomic runs
    // (if any) and hand both ms+genomic dicts to
    // ``writeMinimalToPath:msRuns:genomicRuns:...`` in one call.
    // Replaces the pre-refactor split between a "genomic-only" fast
    // path through writeMinimalToPath: and a "mixed MS+genomic"
    // legacy fallback that wrote MS via writeToFilePath: then
    // re-opened the file RW to append /study/genomic_runs/. Mirrors
    // Python's ``anonymization.anonymize()`` flow which always feeds
    // ``SpectralDataset.write_minimal(..., runs=new_runs,
    // genomic_runs=new_genomic_runs)`` regardless of which
    // modalities are populated.
    NSMutableArray<NSString *> *applied =
        [NSMutableArray arrayWithArray:result.policiesApplied ?: @[]];
    NSDictionary<NSString *, TTIOWrittenGenomicRun *> *transformedGenomic =
        @{};
    if (source.genomicRuns.count > 0) {
        transformedGenomic =
            [self _buildTransformedGenomicRuns:source
                                          policy:policy
                                          result:result
                                     appliedList:applied
                                           error:error];
        if (transformedGenomic == nil) return nil;
        result.policiesApplied = applied;
    }

    if (![TTIOSpectralDataset
            writeMinimalToPath:outputPath
                          title:title
            isaInvestigationId:source.isaInvestigationId
                        msRuns:newRuns
                    genomicRuns:transformedGenomic
                identifications:source.identifications
                quantifications:source.quantifications
              provenanceRecords:@[prov]
                          error:error]) {
        return nil;
    }
    return result;
}


#pragma mark - M90.3 Genomic policies

+ (NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)
    _buildTransformedGenomicRuns:(TTIOSpectralDataset *)source
                            policy:(TTIOAnonymizationPolicy *)policy
                            result:(TTIOAnonymizationResult *)result
                       appliedList:(NSMutableArray<NSString *> *)applied
                             error:(NSError **)error
{
    NSDictionary *built =
        [self _buildOrApplyGenomicPolicies:source
                                       policy:policy
                                       result:result
                                  appliedList:applied
                                        error:error];
    return built;
}

+ (NSDictionary<NSString *, TTIOWrittenGenomicRun *> *)
    _buildOrApplyGenomicPolicies:(TTIOSpectralDataset *)source
                            policy:(TTIOAnonymizationPolicy *)policy
                            result:(TTIOAnonymizationResult *)result
                       appliedList:(NSMutableArray<NSString *> *)applied
                             error:(NSError **)error
{
    NSArray<NSString *> *runNames =
        [[source.genomicRuns allKeys] sortedArrayUsingSelector:@selector(compare:)];

    NSMutableDictionary<NSString *, TTIOWrittenGenomicRun *> *transformed =
        [NSMutableDictionary dictionary];

    for (NSString *runName in runNames) {
        TTIOGenomicRun *gr = source.genomicRuns[runName];
        NSUInteger n = gr.readCount;

        // Materialise per-read fields by iterating the lazy run.
        // O(N) reads — anonymizer is one-shot offline so the cost
        // is acceptable.
        NSMutableArray<NSString *> *readNames =
            [NSMutableArray arrayWithCapacity:n];
        NSMutableArray<NSString *> *cigars =
            [NSMutableArray arrayWithCapacity:n];
        NSMutableArray<NSData *> *sequencesList =
            [NSMutableArray arrayWithCapacity:n];
        NSMutableArray<NSMutableData *> *qualitiesList =
            [NSMutableArray arrayWithCapacity:n];
        NSMutableArray<NSString *> *mateChromosomes =
            [NSMutableArray arrayWithCapacity:n];
        NSMutableData *matePositionsData =
            [NSMutableData dataWithLength:n * sizeof(int64_t)];
        int64_t *matePositions = (int64_t *)matePositionsData.mutableBytes;
        NSMutableData *templateLengthsData =
            [NSMutableData dataWithLength:n * sizeof(int32_t)];
        int32_t *templateLengths =
            (int32_t *)templateLengthsData.mutableBytes;

        for (NSUInteger i = 0; i < n; i++) {
            NSError *readErr = nil;
            TTIOAlignedRead *r = [gr readAtIndex:i error:&readErr];
            if (!r) {
                if (error) *error = readErr;
                return NO;
            }
            [readNames addObject:r.readName ?: @""];
            [cigars    addObject:r.cigar ?: @""];
            NSData *seqBytes =
                [r.sequence dataUsingEncoding:NSASCIIStringEncoding] ?: [NSData data];
            [sequencesList addObject:seqBytes];
            NSMutableData *q = [r.qualities mutableCopy] ?: [NSMutableData data];
            [qualitiesList addObject:q];
            [mateChromosomes addObject:r.mateChromosome ?: @""];
            matePositions[i] = r.matePosition;
            templateLengths[i] = r.templateLength;
        }

        // ── strip_read_names ─────────────────────────────────────
        if (policy.stripReadNames) {
            for (NSUInteger i = 0; i < n; i++) {
                [readNames replaceObjectAtIndex:i withObject:@""];
            }
            result.readNamesStripped += n;
            if (![applied containsObject:@"strip_read_names"]) {
                [applied addObject:@"strip_read_names"];
            }
        }

        // ── randomise_qualities ──────────────────────────────────
        if (policy.randomiseQualities) {
            if (policy.randomiseQualitiesSeed != nil) {
                // seeded RNG path. Range [0, 93] matches the
                // SAM spec valid Phred range. Reproducible within
                // ObjC; cross-language byte equality is NOT a goal.
                TTIORngState rng;
                uint64_t seed =
                    (uint64_t)[policy.randomiseQualitiesSeed unsignedLongLongValue];
                ttioRngSeed(&rng, seed);
                for (NSUInteger i = 0; i < n; i++) {
                    NSMutableData *q = qualitiesList[i];
                    uint8_t *bytes = (uint8_t *)q.mutableBytes;
                    NSUInteger qlen = q.length;
                    for (NSUInteger j = 0; j < qlen; j++) {
                        bytes[j] = (uint8_t)ttioRngNextBoundedU32(&rng, 94);
                    }
                }
            } else {
                uint8_t k = policy.randomiseQualitiesConstant;
                for (NSUInteger i = 0; i < n; i++) {
                    NSMutableData *q = qualitiesList[i];
                    uint8_t *bytes = (uint8_t *)q.mutableBytes;
                    memset(bytes, k, q.length);
                }
            }
            result.qualitiesRandomised += n;
            if (![applied containsObject:@"randomise_qualities"]) {
                [applied addObject:@"randomise_qualities"];
            }
        }

        // ── mask_regions (M90.13: SAM-overlap by CIGAR walk) ─────
        if (policy.maskRegions.count > 0) {
            NSMutableIndexSet *alreadyMasked = [NSMutableIndexSet indexSet];
            for (NSArray *region in policy.maskRegions) {
                if (region.count != 3) continue;
                NSString *rChr = region[0];
                int64_t rStart = [region[1] longLongValue];
                int64_t rEnd   = [region[2] longLongValue];
                for (NSUInteger i = 0; i < n; i++) {
                    if ([alreadyMasked containsIndex:i]) continue;
                    NSString *chromI = [gr.index chromosomeAt:i];
                    if (![chromI isEqualToString:rChr]) continue;
                    int64_t posI = [gr.index positionAt:i];
                    int64_t span = cigarRefSpan(cigars[i]);
                    BOOL hit = NO;
                    if (span > 0) {
                        // SAM-overlap: read covers [pos, pos+span-1].
                        int64_t readEnd = posI + span - 1;
                        if (!(readEnd < rStart || posI > rEnd)) {
                            hit = YES;
                        }
                    } else {
                        // Empty / unparseable CIGAR — M90.3 fallback.
                        if (posI >= rStart && posI <= rEnd) {
                            hit = YES;
                        }
                    }
                    if (hit) {
                        NSMutableData *zeroSeq = [NSMutableData
                            dataWithLength:[sequencesList[i] length]];
                        [sequencesList replaceObjectAtIndex:i
                                                  withObject:zeroSeq];
                        NSMutableData *q = qualitiesList[i];
                        memset(q.mutableBytes, 0, q.length);
                        result.readsInMaskedRegion += 1;
                        [alreadyMasked addIndex:i];
                    }
                }
            }
            if (![applied containsObject:@"mask_regions"]) {
                [applied addObject:@"mask_regions"];
            }
        }

        // ── Re-pack into the flat WrittenGenomicRun layout. ──────
        NSMutableData *lengthsData =
            [NSMutableData dataWithLength:n * sizeof(uint32_t)];
        uint32_t *lengths = (uint32_t *)lengthsData.mutableBytes;
        NSMutableData *offsetsData =
            [NSMutableData dataWithLength:n * sizeof(uint64_t)];
        uint64_t *offsets = (uint64_t *)offsetsData.mutableBytes;
        uint64_t running = 0;
        NSMutableData *sequencesFlat = [NSMutableData data];
        NSMutableData *qualitiesFlat = [NSMutableData data];
        for (NSUInteger i = 0; i < n; i++) {
            NSData *s = sequencesList[i];
            offsets[i] = running;
            lengths[i] = (uint32_t)s.length;
            running += s.length;
            [sequencesFlat appendData:s];
            [qualitiesFlat appendData:qualitiesList[i]];
        }

        // Repack positions / mapping_qualities / flags / chromosomes
        // from the source index (no per-read transform on these — only
        // sequences/qualities/read_names/cigars are touched).
        NSMutableData *positionsData =
            [NSMutableData dataWithLength:n * sizeof(int64_t)];
        int64_t *positions = (int64_t *)positionsData.mutableBytes;
        NSMutableData *mapqsData =
            [NSMutableData dataWithLength:n * sizeof(uint8_t)];
        uint8_t *mapqs = (uint8_t *)mapqsData.mutableBytes;
        NSMutableData *flagsData =
            [NSMutableData dataWithLength:n * sizeof(uint32_t)];
        uint32_t *flags = (uint32_t *)flagsData.mutableBytes;
        NSMutableArray<NSString *> *chromosomes =
            [NSMutableArray arrayWithCapacity:n];
        for (NSUInteger i = 0; i < n; i++) {
            positions[i] = [gr.index positionAt:i];
            mapqs[i]     = [gr.index mappingQualityAt:i];
            flags[i]     = [gr.index flagsAt:i];
            [chromosomes addObject:[gr.index chromosomeAt:i] ?: @""];
        }

        TTIOWrittenGenomicRun *written =
            [[TTIOWrittenGenomicRun alloc]
                initWithAcquisitionMode:gr.acquisitionMode
                           referenceUri:gr.referenceUri ?: @""
                               platform:gr.platform ?: @""
                             sampleName:gr.sampleName ?: @""
                              positions:positionsData
                       mappingQualities:mapqsData
                                  flags:flagsData
                              sequences:sequencesFlat
                              qualities:qualitiesFlat
                                offsets:offsetsData
                                lengths:lengthsData
                                 cigars:cigars
                              readNames:readNames
                        mateChromosomes:mateChromosomes
                          matePositions:matePositionsData
                        templateLengths:templateLengthsData
                            chromosomes:chromosomes
                     signalCompression:TTIOCompressionZlib];
        transformed[runName] = written;
    }

    return transformed;
}

@end
