/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "Export/TTIORunSelection.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "ValueClasses/TTIOEnums.h"

static NSString *const kTTIORunSelectionErrorDomain = @"TTIORunSelection";

/** Python's <code>spectrum_class == "TTIONMRSpectrum"</code> discriminant. */
static NSString *const kTTIONMRSpectrumClass = @"TTIONMRSpectrum";

@implementation TTIORunSelection

/* Build an NSError whose localizedDescription is `message` — mirroring
 * how Python raises KeyError/ValueError with that exact text. */
static NSError *rsError(NSString *message)
{
    return [NSError errorWithDomain:kTTIORunSelectionErrorDomain
                               code:1
                           userInfo:@{NSLocalizedDescriptionKey: message}];
}

/* ", ".join(sorted(runs)) — sorted, comma-space-joined run names. */
static NSString *rsSortedNames(NSDictionary *runs)
{
    NSArray *sorted =
        [runs.allKeys sortedArrayUsingSelector:@selector(compare:)];
    return [sorted componentsJoinedByString:@", "];
}

+ (TTIOAcquisitionRun *)analyticalRunIn:(TTIOSpectralDataset *)ds
                                  layer:(NSString *)layer
                                  error:(NSError **)error
{
    NSDictionary<NSString *, TTIOAcquisitionRun *> *runs = ds.msRuns;
    if (runs.count == 0) {
        if (error) *error = rsError(@"no analytical runs in dataset");
        return nil;
    }
    if (layer.length > 0) {
        TTIOAcquisitionRun *r = runs[layer];
        if (r == nil) {
            if (error) *error = rsError([NSString stringWithFormat:
                @"run '%@' not found; have: %@", layer, rsSortedNames(runs)]);
            return nil;
        }
        return r;
    }
    if (runs.count == 1) {
        return runs.allValues.firstObject;
    }
    if (error) *error = rsError(@"multiple runs present; pass --layer <name>");
    return nil;
}

+ (TTIOAcquisitionRun *)nmrRunIn:(TTIOSpectralDataset *)ds
                           layer:(NSString *)layer
                           error:(NSError **)error
{
    NSDictionary<NSString *, TTIOAcquisitionRun *> *runs = ds.msRuns;
    if (runs.count == 0) {
        if (error) *error = rsError(@"no analytical runs in dataset");
        return nil;
    }
    if (layer.length > 0) {
        TTIOAcquisitionRun *r = runs[layer];
        if (r == nil) {
            if (error) *error = rsError([NSString stringWithFormat:
                @"run '%@' not found; have: %@", layer, rsSortedNames(runs)]);
            return nil;
        }
        return r;
    }
    NSMutableArray<TTIOAcquisitionRun *> *nmr = [NSMutableArray array];
    for (TTIOAcquisitionRun *r in runs.allValues) {
        if ([r.spectrumClassName isEqualToString:kTTIONMRSpectrumClass]) {
            [nmr addObject:r];
        }
    }
    if (nmr.count == 1) {
        return nmr.firstObject;
    }
    if (nmr.count > 1) {
        if (error) *error =
            rsError(@"multiple NMR runs present; pass --layer <name>");
        return nil;
    }
    if (runs.count == 1) {
        return runs.allValues.firstObject;
    }
    if (error) *error = rsError(@"multiple runs present; pass --layer <name>");
    return nil;
}

+ (TTIOGenomicRun *)genomicRunIn:(TTIOSpectralDataset *)ds
                           layer:(NSString *)layer
                           error:(NSError **)error
{
    NSDictionary<NSString *, TTIOGenomicRun *> *runs = ds.genomicRuns;
    if (runs.count == 0) {
        if (error) *error = rsError(@"no genomic runs in dataset");
        return nil;
    }
    if (layer.length > 0) {
        TTIOGenomicRun *r = runs[layer];
        if (r == nil) {
            if (error) *error = rsError([NSString stringWithFormat:
                @"genomic run '%@' not found; have: %@",
                layer, rsSortedNames(runs)]);
            return nil;
        }
        return r;
    }
    if (runs.count == 1) {
        return runs.allValues.firstObject;
    }
    if (error) *error = rsError([NSString stringWithFormat:
        @"multiple genomic runs present; pass --layer <name>: %@",
        rsSortedNames(runs)]);
    return nil;
}

+ (TTIOWrittenGenomicRun *)writtenFromGenomicRun:(TTIOGenomicRun *)readSideRun
{
    NSUInteger n = readSideRun.readCount;
    TTIOGenomicIndex *idx = readSideRun.index;

    NSMutableData *positions  = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    NSMutableData *mapqs      = [NSMutableData dataWithLength:n * sizeof(uint8_t)];
    NSMutableData *flags      = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *offsets    = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    NSMutableData *lengths    = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *matePos    = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    NSMutableData *tlens      = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    int64_t  *posPtr  = (int64_t  *)positions.mutableBytes;
    uint8_t  *mapqPtr = (uint8_t  *)mapqs.mutableBytes;
    uint32_t *flagPtr = (uint32_t *)flags.mutableBytes;
    uint64_t *offPtr  = (uint64_t *)offsets.mutableBytes;
    uint32_t *lenPtr  = (uint32_t *)lengths.mutableBytes;
    int64_t  *mposPtr = (int64_t  *)matePos.mutableBytes;
    int32_t  *tlenPtr = (int32_t  *)tlens.mutableBytes;

    NSMutableArray<NSString *> *chromosomes    = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSString *> *readNames      = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSString *> *cigars         = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray<NSString *> *mateChromosomes = [NSMutableArray arrayWithCapacity:n];
    // M101: the SAM optional fields, carried when the run stores them.
    BOOL hasTags = [readSideRun hasTagsChannel];
    NSMutableArray<NSString *> *tags = hasTags ? [NSMutableArray arrayWithCapacity:n] : nil;

    NSMutableData *sequences = [NSMutableData data];
    NSMutableData *qualities = [NSMutableData data];

    for (NSUInteger i = 0; i < n; i++) {
        posPtr[i]  = [idx positionAt:i];
        mapqPtr[i] = [idx mappingQualityAt:i];
        flagPtr[i] = [idx flagsAt:i];
        [chromosomes addObject:([idx chromosomeAt:i] ?: @"*")];
    }
    // Reads in input order, like the index above; the sequential walk
    // decodes each block once, also for a grouped run (M103).
    __block uint64_t running = 0;
    [readSideRun iterReadsFrom:0 to:n error:NULL
                    usingBlock:^(TTIOAlignedRead *read, NSUInteger i, BOOL *stop) {
        (void)stop;
        NSString *seq  = read.sequence ?: @"";
        NSData   *qual = read.qualities ?: [NSData data];
        NSData   *seqBytes = [seq dataUsingEncoding:NSASCIIStringEncoding]
                             ?: [NSData data];

        // offsets/lengths slice the concatenated sequence/quality buffers.
        offPtr[i] = running;
        lenPtr[i] = (uint32_t)seqBytes.length;
        running  += (uint64_t)seqBytes.length;
        [sequences appendData:seqBytes];
        [qualities appendData:qual];

        [readNames addObject:(read.readName ?: @"")];
        [cigars addObject:(read.cigar ?: @"*")];
        [mateChromosomes addObject:(read.mateChromosome ?: @"*")];
        [tags addObject:(read.tags ?: @"")];
        mposPtr[i] = read.matePosition;
        tlenPtr[i] = read.templateLength;
    }];
    // A read that fails to decode ends the walk; the rest stay empty,
    // as a nil read did before.
    for (NSUInteger i = readNames.count; i < n; i++) {
        offPtr[i] = running;
        lenPtr[i] = 0;
        [readNames addObject:@""];
        [cigars addObject:@"*"];
        [mateChromosomes addObject:@"*"];
        [tags addObject:@""];
    }

    TTIOAcquisitionMode mode = readSideRun.acquisitionMode;
    if (mode == 0) mode = TTIOAcquisitionModeGenomicWGS;

    TTIOWrittenGenomicRun *w = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:mode
                   referenceUri:(readSideRun.referenceUri ?: @"")
                       platform:(readSideRun.platform ?: @"")
                     sampleName:(readSideRun.sampleName ?: @"")
                      positions:positions
               mappingQualities:mapqs
                          flags:flags
                      sequences:sequences
                      qualities:qualities
                        offsets:offsets
                        lengths:lengths
                         cigars:cigars
                      readNames:readNames
                mateChromosomes:mateChromosomes
                  matePositions:matePos
                templateLengths:tlens
                    chromosomes:chromosomes
              signalCompression:TTIOCompressionNone];
    w.tags = tags;
    return w;
}

@end
