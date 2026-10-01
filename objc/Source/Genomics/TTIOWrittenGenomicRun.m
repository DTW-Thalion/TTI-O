/*
 * TTIOWrittenGenomicRun.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOWrittenGenomicRun
 * Inherits From: NSObject
 * Declared In:   Genomics/TTIOWrittenGenomicRun.h
 *
 * Write-side container for a single genomic run. Pure data class —
 * accessors plus the two designated initialisers. The writer
 * (TTIOSpectralDataset) consumes these to produce the on-disk
 * channel layout described in docs/format-spec.md §10.4-§10.10.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOWrittenGenomicRun.h"

@implementation TTIOWrittenGenomicRun

- (instancetype)initWithAcquisitionMode:(TTIOAcquisitionMode)mode
                            referenceUri:(NSString *)referenceUri
                                platform:(NSString *)platform
                              sampleName:(NSString *)sampleName
                                positions:(NSData *)positions
                         mappingQualities:(NSData *)mappingQualities
                                    flags:(NSData *)flags
                                sequences:(NSData *)sequences
                                qualities:(NSData *)qualities
                                  offsets:(NSData *)offsets
                                  lengths:(NSData *)lengths
                                   cigars:(NSArray<NSString *> *)cigars
                                readNames:(NSArray<NSString *> *)readNames
                          mateChromosomes:(NSArray<NSString *> *)mateChromosomes
                            matePositions:(NSData *)matePositions
                          templateLengths:(NSData *)templateLengths
                              chromosomes:(NSArray<NSString *> *)chromosomes
                       signalCompression:(TTIOCompression)signalCompression
{
    return [self initWithAcquisitionMode:mode
                            referenceUri:referenceUri
                                platform:platform
                              sampleName:sampleName
                                positions:positions
                         mappingQualities:mappingQualities
                                    flags:flags
                                sequences:sequences
                                qualities:qualities
                                  offsets:offsets
                                  lengths:lengths
                                   cigars:cigars
                                readNames:readNames
                          mateChromosomes:mateChromosomes
                            matePositions:matePositions
                          templateLengths:templateLengths
                              chromosomes:chromosomes
                       signalCompression:signalCompression
                     signalCodecOverrides:@{}];
}

- (instancetype)initWithAcquisitionMode:(TTIOAcquisitionMode)mode
                            referenceUri:(NSString *)referenceUri
                                platform:(NSString *)platform
                              sampleName:(NSString *)sampleName
                                positions:(NSData *)positions
                         mappingQualities:(NSData *)mappingQualities
                                    flags:(NSData *)flags
                                sequences:(NSData *)sequences
                                qualities:(NSData *)qualities
                                  offsets:(NSData *)offsets
                                  lengths:(NSData *)lengths
                                   cigars:(NSArray<NSString *> *)cigars
                                readNames:(NSArray<NSString *> *)readNames
                          mateChromosomes:(NSArray<NSString *> *)mateChromosomes
                            matePositions:(NSData *)matePositions
                          templateLengths:(NSData *)templateLengths
                              chromosomes:(NSArray<NSString *> *)chromosomes
                       signalCompression:(TTIOCompression)signalCompression
                     signalCodecOverrides:(NSDictionary<NSString *, NSNumber *> *)signalCodecOverrides
{
    self = [super init];
    if (self) {
        _acquisitionMode      = mode;
        _referenceUri         = [referenceUri copy];
        _platform             = [platform copy];
        _sampleName           = [sampleName copy];
        _positionsData        = [positions copy];
        _mappingQualitiesData = [mappingQualities copy];
        _flagsData            = [flags copy];
        _sequencesData        = [sequences copy];
        _qualitiesData        = [qualities copy];
        _offsetsData          = [offsets copy];
        _lengthsData          = [lengths copy];
        _cigars               = [cigars copy];
        _readNames            = [readNames copy];
        _mateChromosomes      = [mateChromosomes copy];
        _matePositionsData    = [matePositions copy];
        _templateLengthsData  = [templateLengths copy];
        _chromosomes          = [chromosomes copy];
        _signalCompression    = signalCompression;
        _signalCodecOverrides = signalCodecOverrides
            ? [signalCodecOverrides copy]
            : @{};
        _provenanceRecords    = @[];
        // L3 (Task #82 Phase B.1, 2026-05-01): embedReference now
        // defaults to NO so chr22-style benchmarks don't carry the
        // ~10 MB embedded reference blob by default. CRAM 3.1's
        // default is also external-reference; users who want
        // self-contained .tio files set embedReference = YES
        // explicitly.
        _embedReference        = NO;
        _referenceChromSeqs    = nil;
        _externalReferencePath = nil;
    }
    return self;
}

- (NSUInteger)readCount
{
    return _offsetsData.length / sizeof(uint64_t);
}

- (BOOL)hasTags
{
    for (NSString *t in _tags) {
        if (t.length > 0) return YES;
    }
    return NO;
}

- (instancetype)copyWithSignalCodecOverrides:(NSDictionary<NSString *, NSNumber *> *)overrides
{
    TTIOWrittenGenomicRun *c = [[[self class] alloc]
        initWithAcquisitionMode:_acquisitionMode
                   referenceUri:_referenceUri
                       platform:_platform
                     sampleName:_sampleName
                      positions:_positionsData
               mappingQualities:_mappingQualitiesData
                          flags:_flagsData
                      sequences:_sequencesData
                      qualities:_qualitiesData
                        offsets:_offsetsData
                        lengths:_lengthsData
                         cigars:_cigars
                      readNames:_readNames
                mateChromosomes:_mateChromosomes
                  matePositions:_matePositionsData
                templateLengths:_templateLengthsData
                    chromosomes:_chromosomes
              signalCompression:_signalCompression
           signalCodecOverrides:overrides];
    c.optDisableQualitiesV5 = _optDisableQualitiesV5;
    c.provenanceRecords = _provenanceRecords;
    c.embedReference = _embedReference;
    c.referenceChromSeqs = _referenceChromSeqs;
    c.externalReferencePath = _externalReferencePath;
    c.bulkV2Blobs = _bulkV2Blobs;
    c.optLegacyWholeChannel = _optLegacyWholeChannel;
    c.readRole = _readRole;
    c.refDiffSliceBytes = _refDiffSliceBytes;
    c.tags = _tags;
    return c;
}

- (instancetype)copyWithProvenance:(NSArray<TTIOProvenanceRecord *> *)records
{
    TTIOWrittenGenomicRun *c = [self copyWithSignalCodecOverrides:_signalCodecOverrides];
    c.provenanceRecords = records ?: @[];
    return c;
}

- (instancetype)copyWithOptLegacyWholeChannel:(BOOL)legacy
{
    TTIOWrittenGenomicRun *c = [self copyWithSignalCodecOverrides:_signalCodecOverrides];
    c.optLegacyWholeChannel = legacy;
    return c;
}

@end
