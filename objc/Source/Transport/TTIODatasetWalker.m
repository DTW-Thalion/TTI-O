/*
 * TTIODatasetWalker.m
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * Filtered dataset walker — produces transport-stream events
 * (StreamHeader → DatasetHeaders → AccessUnits → EndOfDataset →
 * EndOfStream) and dispatches them to a visitor. Owners decide what
 * to do with each event: encode as transport packets, build stats
 * records, count matches, etc.
 *
 * The walker is intentionally independent of TTIOTransportWriter —
 * downstream code can drive it without producing any wire bytes.
 */
#import "TTIODatasetWalker.h"

#import "Core/TTIOPortability.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Run/TTIOAcquisitionRun.h"
#import "Run/TTIOInstrumentConfig.h"
#import "Run/TTIOSpectrumIndex.h"
#import "Spectra/TTIOSpectrum.h"
#import "Spectra/TTIOMassSpectrum.h"
#import "Core/TTIOSignalArray.h"
#import "ValueClasses/TTIOEnums.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOAlignedRead.h"
#import "TTIOAUFilter.h"
#import "TTIOAccessUnit.h"
#import "Codecs/TTIORans.h"
#import "Codecs/TTIOBasePack.h"
// #140: v0.11 prelude accessors and image-cube types reachable via
// TTIOSpectralDataset's read-only properties.
#import "Genomics/TTIOReferenceImport.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "Dataset/TTIOSubject.h"
#import "Dataset/TTIOSample.h"
#import "Image/TTIOImage.h"
#import "Image/TTIOMSImage.h"
#import "Image/TTIORamanImage.h"
#import "Image/TTIOIRImage.h"
#import "Dataset/TTIOIdentification.h"
#import "Dataset/TTIOQuantification.h"

// ---------------------------------------------------------------- helpers

static uint8_t wireFromSpectrumClassName(NSString *name)
{
    if ([name isEqualToString:@"TTIONMRSpectrum"])      return 1;
    if ([name isEqualToString:@"TTIONMR2D"])            return 2;
    if ([name isEqualToString:@"TTIOFID"])              return 3;
    if ([name isEqualToString:@"TTIOMSImagePixel"])     return 4;
    if ([name isEqualToString:@"TTIOGenomicRead"])      return 5;
    return 0;  // TTIOMassSpectrum / default
}

static uint8_t wireFromPolarity(TTIOPolarity p)
{
    switch (p) {
        case TTIOPolarityPositive: return 0;
        case TTIOPolarityNegative: return 1;
        case TTIOPolarityUnknown:
        default:                   return 2;
    }
}

static NSString *instrumentConfigJSON(TTIOInstrumentConfig *cfg)
{
    if (!cfg) return @"{}";
    // Mirrors the field set produced by
    // TTIOTransportWriter.m's instrumentConfigJSON so the
    // visitor-driven walk emits byte-identical dataset headers.
    NSDictionary *d = @{
        @"analyzer_type": cfg.analyzerType ?: @"",
        @"detector_type": cfg.detectorType ?: @"",
        @"manufacturer":  cfg.manufacturer ?: @"",
        @"model":         cfg.model        ?: @"",
        @"serial_number": cfg.serialNumber ?: @"",
        @"source_type":   cfg.sourceType   ?: @"",
    };
    NSData *data = [NSJSONSerialization dataWithJSONObject:d
                                                    options:TTIO_JSON_SORTED_KEYS
                                                      error:NULL];
    return data
        ? [[NSString alloc] initWithData:data encoding:NSUTF8StringEncoding]
        : @"{}";
}

static NSString *genomicRunMetadataJSON(TTIOGenomicRun *run)
{
    if (!run) return @"{}";
    NSDictionary *d = @{
        @"modality":      run.modality      ?: @"",
        @"platform":      run.platform      ?: @"",
        /* M97: always present, "" when the run has no read role, so
         * the JSON stays byte-identical across the 3 SDK writers. */
        @"read_role":     run.readRole      ?: @"",
        @"reference_uri": run.referenceUri  ?: @"",
        @"sample_name":   run.sampleName    ?: @"",
    };
    NSData *data = [NSJSONSerialization dataWithJSONObject:d
                                                    options:TTIO_JSON_SORTED_KEYS
                                                      error:NULL];
    return data
        ? [[NSString alloc] initWithData:data encoding:NSUTF8StringEncoding]
        : @"{}";
}

// #140: mirror TTIOTransportWriter.m's applyWireCodecGenomic so the
// walker's genomic AU emission produces wire-equivalent channel bytes
// to the writer's writeDataset: / writeGenomicRun: path.
static NSData *walkerApplyWireCodecGenomic(NSData *plaintext, uint8_t codec)
{
    if (codec == TTIOCompressionNone) return plaintext;
    switch (codec) {
        case TTIOCompressionRansOrder0:
            return TTIORansEncode(plaintext, 0);
        case TTIOCompressionRansOrder1:
            return TTIORansEncode(plaintext, 1);
        case TTIOCompressionBasePack:
            return TTIOBasePackEncode(plaintext);
        default:
            return plaintext;
    }
}

static TTIOAccessUnit *accessUnitFromGenomicRead(TTIOGenomicRun *run,
                                                  NSUInteger readIndex,
                                                  TTIOAlignedRead *r,
                                                  NSData *seqData,
                                                  NSData *qualData,
                                                  NSString *nameStr,
                                                  uint8_t seqCodec,
                                                  uint8_t qualCodec)
{
    uint32_t seqLen = (uint32_t)seqData.length;
    uint32_t qualLen = (uint32_t)qualData.length;
    NSData *seqPayload = walkerApplyWireCodecGenomic(seqData, seqCodec);
    NSData *qualPayload = walkerApplyWireCodecGenomic(qualData, qualCodec);
    TTIOTransportChannelData *seqCh =
        [[TTIOTransportChannelData alloc] initWithName:@"sequences"
                                              precision:TTIOPrecisionUInt8
                                            compression:seqCodec
                                              nElements:seqLen
                                                   data:seqPayload];
    TTIOTransportChannelData *qualCh =
        [[TTIOTransportChannelData alloc] initWithName:@"qualities"
                                              precision:TTIOPrecisionUInt8
                                            compression:qualCodec
                                              nElements:qualLen
                                                   data:qualPayload];
    NSData *cigarData =
        [(r.cigar ?: @"") dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
    NSData *nameData =
        [(nameStr ?: @"") dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
    NSData *mateChrData =
        [(r.mateChromosome ?: @"") dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
    TTIOTransportChannelData *cigarCh =
        [[TTIOTransportChannelData alloc] initWithName:@"cigar"
                                              precision:TTIOPrecisionUInt8
                                            compression:TTIOCompressionNone
                                              nElements:(uint32_t)cigarData.length
                                                   data:cigarData];
    TTIOTransportChannelData *nameCh =
        [[TTIOTransportChannelData alloc] initWithName:@"read_name"
                                              precision:TTIOPrecisionUInt8
                                            compression:TTIOCompressionNone
                                              nElements:(uint32_t)nameData.length
                                                   data:nameData];
    TTIOTransportChannelData *mateChrCh =
        [[TTIOTransportChannelData alloc] initWithName:@"mate_chromosome"
                                              precision:TTIOPrecisionUInt8
                                            compression:TTIOCompressionNone
                                              nElements:(uint32_t)mateChrData.length
                                                   data:mateChrData];
    // M101: a run carrying SAM tags adds a UTF-8 "tags" channel.
    NSMutableArray *auChannels = [@[seqCh, qualCh, cigarCh, nameCh, mateChrCh] mutableCopy];
    if ([run hasTagsChannel]) {
        NSData *tagData =
            [(r.tags ?: @"") dataUsingEncoding:NSUTF8StringEncoding] ?: [NSData data];
        [auChannels addObject:[[TTIOTransportChannelData alloc] initWithName:@"tags"
                                                                   precision:TTIOPrecisionUInt8
                                                                 compression:TTIOCompressionNone
                                                                   nElements:(uint32_t)tagData.length
                                                                        data:tagData]];
    }
    // Mirror the writer: prefer the genomic index for chromosome /
    // position / mappingQuality / flags when available; otherwise fall
    // back to the AlignedRead fields.
    NSString *chrom = r.chromosome;
    int64_t pos = r.position;
    uint8_t mapq = r.mappingQuality;
    uint16_t flags = (uint16_t)(r.flags & 0xFFFFu);
    TTIOGenomicIndex *idx = run.index;
    if (idx && readIndex < idx.count) {
        chrom = [idx chromosomeAt:readIndex] ?: chrom;
        pos = [idx positionAt:readIndex];
        mapq = [idx mappingQualityAt:readIndex];
        flags = (uint16_t)([idx flagsAt:readIndex] & 0xFFFFu);
    }
    return [[TTIOAccessUnit alloc]
        initWithSpectrumClass:5
              acquisitionMode:(uint8_t)run.acquisitionMode
                      msLevel:0
                     polarity:2
                retentionTime:0.0
                  precursorMz:0.0
              precursorCharge:0
                  ionMobility:0.0
            basePeakIntensity:0.0
                     channels:auChannels
                       pixelX:0 pixelY:0 pixelZ:0
                   chromosome:(chrom ?: @"")
                     position:pos
               mappingQuality:mapq
                        flags:flags
                 matePosition:r.matePosition
               templateLength:r.templateLength];
}

static TTIOAccessUnit *accessUnitFromSpectrum(TTIOSpectrum *spectrum,
                                                TTIOAcquisitionRun *run,
                                                NSArray<NSString *> *channelNames)
{
    uint8_t wireClass = wireFromSpectrumClassName(run.spectrumClassName);
    uint8_t msLevel = 0;
    uint8_t polarityWire = 2;
    if ([spectrum isKindOfClass:[TTIOMassSpectrum class]]) {
        TTIOMassSpectrum *ms = (TTIOMassSpectrum *)spectrum;
        msLevel = (uint8_t)MIN((NSUInteger)255, ms.msLevel);
        polarityWire = wireFromPolarity(ms.polarity);
    }
    double bpi = 0.0;
    if (run.spectrumIndex && spectrum.indexPosition < run.spectrumIndex.count) {
        bpi = [run.spectrumIndex basePeakIntensityAt:spectrum.indexPosition];
    }
    NSMutableArray<TTIOTransportChannelData *> *channels = [NSMutableArray array];
    for (NSString *cname in channelNames) {
        TTIOSignalArray *sa = spectrum.signalArrays[cname];
        if (!sa) continue;
        NSData *raw = sa.buffer ?: [NSData data];
        uint32_t nElements = (uint32_t)(raw.length / 8);
        TTIOTransportChannelData *ch =
            [[TTIOTransportChannelData alloc]
                initWithName:cname
                   precision:TTIOPrecisionFloat64
                 compression:TTIOCompressionNone
                   nElements:nElements
                        data:raw];
        [channels addObject:ch];
    }
    return [[TTIOAccessUnit alloc]
        initWithSpectrumClass:wireClass
                acquisitionMode:(uint8_t)run.acquisitionMode
                        msLevel:msLevel
                       polarity:polarityWire
                  retentionTime:spectrum.scanTimeSeconds
                    precursorMz:spectrum.precursorMz
                precursorCharge:(uint8_t)MIN((NSUInteger)255, spectrum.precursorCharge)
                    ionMobility:0.0
              basePeakIntensity:bpi
                       channels:channels
                         pixelX:0 pixelY:0 pixelZ:0];
}

// ---------------------------------------------------------------- walker

@implementation TTIODatasetWalker

- (BOOL)walkDataset:(TTIOSpectralDataset *)dataset
              filter:(TTIOAUFilter *)filter
             visitor:(id<TTIOTransportEventVisitor>)visitor
               error:(NSError **)error
{
    (void)error;  // No structural errors today; reserved for future.
    if (!dataset || !visitor) return NO;

    NSArray<NSString *> *msNames = [dataset.msRuns.allKeys
        sortedArrayUsingSelector:@selector(compare:)];
    NSArray<NSString *> *gNames = [dataset.genomicRuns.allKeys
        sortedArrayUsingSelector:@selector(compare:)];
    NSUInteger nDatasets = msNames.count + gNames.count;

    // ── 1. StreamHeader ────────────────────────────────────────────
    if ([visitor respondsToSelector:
            @selector(walker:visitStreamHeaderWithFormatVersion:title:isaInvestigation:features:nDatasets:)]) {
        [visitor walker:self
visitStreamHeaderWithFormatVersion:@"1.2"
                                title:dataset.title ?: @""
                     isaInvestigation:dataset.isaInvestigationId ?: @""
                             features:@[]
                           nDatasets:(uint16_t)nDatasets];
    }

    // ── 1.5. v0.11 prelude (§5.4) ─────────────────────────────────
    // Sub-sections emitted in spec order, mirroring the gating logic
    // in TTIOTransportWriter writeDataset: (the "v0.11 §5.4 prelude"
    // block). Filed as #140 — walker previously emitted only MS AUs
    // and dropped every v0.11 accessor on the workbench daemon's
    // download path.
    //
    //   §5.4.1 ENCRYPTION_ALGORITHM
    //   §5.4.2 DATASET_PROVENANCE
    //   §5.4.3 SUBJECT_METADATA → SAMPLE_METADATA
    //   §5.4.4 reference groups (one call per import, sorted-uri order)
    //   §5.4.5 image cubes (MS → Raman → IR)
    //   §5.4.6 IDENTIFICATIONS_TABLE → QUANTIFICATIONS_TABLE
    if (dataset.isEncrypted && dataset.encryptedAlgorithm.length > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitEncryptionAlgorithm:)]) {
            [visitor walker:self
visitEncryptionAlgorithm:dataset.encryptedAlgorithm];
        }
    }
    NSArray<TTIOProvenanceRecord *> *provRecords =
        dataset.provenanceRecords ?: @[];
    if (provRecords.count > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitDatasetProvenance:)]) {
            [visitor walker:self visitDatasetProvenance:provRecords];
        }
    }
    NSArray<TTIOSubject *> *subjects = dataset.subjects ?: @[];
    if (subjects.count > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitSubjectMetadata:)]) {
            [visitor walker:self visitSubjectMetadata:subjects];
        }
    }
    NSArray<TTIOSample *> *samples = dataset.samples ?: @[];
    if (samples.count > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitSampleMetadata:)]) {
            [visitor walker:self visitSampleMetadata:samples];
        }
    }
    NSDictionary<NSString *, TTIOReferenceImport *> *refs =
        dataset.references ?: @{};
    if (refs.count > 0) {
        // Sort by URI so cross-call output is reproducible — matches
        // the writer's writeDataset: emission order.
        NSArray<NSString *> *refUris =
            [refs.allKeys sortedArrayUsingSelector:@selector(compare:)];
        if ([visitor respondsToSelector:
                @selector(walker:visitReferenceGroup:)]) {
            for (NSString *uri in refUris) {
                TTIOReferenceImport *ref = refs[uri];
                if (!ref) continue;
                [visitor walker:self visitReferenceGroup:ref];
            }
        }
    }
    // -msImage / -ramanImage / -irImage may return non-nil placeholders
    // (width=0, height=0) when their cube group is absent; gate on a
    // dimension check to match the writer.
    TTIOMSImage *msImg = (TTIOMSImage *)[dataset imageForKind:TTIOImageKindMS];
    if (msImg && msImg.width > 0 && msImg.height > 0) {
        if ([visitor respondsToSelector:@selector(walker:visitImage:)]) {
            [visitor walker:self visitImage:msImg];
        }
    }
    TTIORamanImage *ramanImg = (TTIORamanImage *)[dataset imageForKind:TTIOImageKindRaman];
    if (ramanImg && ramanImg.width > 0 && ramanImg.height > 0) {
        if ([visitor respondsToSelector:@selector(walker:visitRamanImage:)]) {
            [visitor walker:self visitRamanImage:ramanImg];
        }
    }
    TTIOIRImage *irImg = (TTIOIRImage *)[dataset imageForKind:TTIOImageKindIR];
    if (irImg && irImg.width > 0 && irImg.height > 0) {
        if ([visitor respondsToSelector:@selector(walker:visitIRImage:)]) {
            [visitor walker:self visitIRImage:irImg];
        }
    }
    NSArray<TTIOIdentification *> *idents =
        dataset.identifications ?: @[];
    if (idents.count > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitIdentificationsTable:)]) {
            [visitor walker:self visitIdentificationsTable:idents];
        }
    }
    NSArray<TTIOQuantification *> *quants =
        dataset.quantifications ?: @[];
    if (quants.count > 0) {
        if ([visitor respondsToSelector:
                @selector(walker:visitQuantificationsTable:)]) {
            [visitor walker:self visitQuantificationsTable:quants];
        }
    }

    // ── 2. DatasetHeaders (MS runs first, then genomic) ────────────
    uint16_t did = 1;
    NSNumber *filterDid = filter.datasetId;
    for (NSString *name in msNames) {
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOAcquisitionRun *run = dataset.msRuns[name];
        NSArray<NSString *> *channelNames =
            [run valueForKey:@"channelNames"] ?: @[@"mz", @"intensity"];
        if ([visitor respondsToSelector:
                @selector(walker:visitDatasetHeaderWithDatasetId:name:acquisitionMode:spectrumClass:channelNames:instrumentJSON:expectedAUCount:)]) {
            [visitor walker:self
visitDatasetHeaderWithDatasetId:did
                                  name:name
                       acquisitionMode:(uint8_t)run.acquisitionMode
                         spectrumClass:run.spectrumClassName ?: @"TTIOMassSpectrum"
                          channelNames:channelNames
                        instrumentJSON:instrumentConfigJSON(run.instrumentConfig)
                      expectedAUCount:(uint32_t)[run count]];
        }
        did++;
    }
    NSArray<NSString *> *gChannelNames = @[@"sequences", @"qualities",
                                            @"cigar", @"read_name",
                                            @"mate_chromosome"];
    for (NSString *name in gNames) {
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOGenomicRun *grun = dataset.genomicRuns[name];
        if ([visitor respondsToSelector:
                @selector(walker:visitDatasetHeaderWithDatasetId:name:acquisitionMode:spectrumClass:channelNames:instrumentJSON:expectedAUCount:)]) {
            [visitor walker:self
visitDatasetHeaderWithDatasetId:did
                                  name:name
                       acquisitionMode:(uint8_t)grun.acquisitionMode
                         spectrumClass:@"TTIOGenomicRead"
                          channelNames:([grun hasTagsChannel]
                                            ? [gChannelNames arrayByAddingObject:@"tags"]
                                            : gChannelNames)
                        instrumentJSON:genomicRunMetadataJSON(grun)
                      expectedAUCount:(uint32_t)grun.readCount];
        }
        did++;
    }

    // ── 3. AccessUnits ────────────────────────────────────────────
    uint32_t emitted = 0;
    uint32_t maxAU = filter.maxAU ? filter.maxAU.unsignedIntValue : UINT32_MAX;
    did = 1;
    BOOL hasAccessUnitVisitor = [visitor respondsToSelector:
        @selector(walker:visitAccessUnit:datasetId:auSequence:)];

    for (NSString *name in msNames) {
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOAcquisitionRun *run = dataset.msRuns[name];
        NSArray<NSString *> *channelNames =
            [run valueForKey:@"channelNames"] ?: @[@"mz", @"intensity"];
        NSUInteger count = [run count];
        for (NSUInteger i = 0; i < count && emitted < maxAU; i++) {
            TTIOSpectrum *sp = [run objectAtIndex:i];
            TTIOAccessUnit *au =
                accessUnitFromSpectrum(sp, run, channelNames);
            if (filter && ![filter matches:au datasetId:did]) continue;
            if (hasAccessUnitVisitor) {
                [visitor walker:self
             visitAccessUnit:au
                   datasetId:did
                  auSequence:(uint32_t)i];
            }
            emitted++;
        }
        did++;
        if (emitted >= maxAU) break;
    }
    // #140: Genomic AU emission — mirrors writeDataset:'s inline
    // genomic block (TTIOTransportWriter.m). Each read becomes one
    // AccessUnit with 5 channels (sequences, qualities, cigar,
    // read_name, mate_chromosome) plus the genomic-suffix fields
    // (chromosome / position / mappingQuality / flags / matePosition
    // / templateLength). auSequence is per-dataset (#139 ingest now
    // enforces per-dataset monotonicity).
    did = (uint16_t)(msNames.count + 1);
    for (NSString *name in gNames) {
        if (emitted >= maxAU) break;
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOGenomicRun *grun = dataset.genomicRuns[name];
        NSUInteger nReads = grun.readCount;
        TTIOGenomicIndex *gIdx = grun.index;
        uint8_t seqCodec = [grun wireCompressionForChannel:@"sequences"];
        uint8_t qualCodec = [grun wireCompressionForChannel:@"qualities"];
        // Bulk-fetch byte channels once per run (mirrors writer's
        // optimisation — per-record cost dominated by NSData slicing).
        NSData *seqAll = (nReads > 0)
            ? [grun wholeSequencesData] : [NSData data];
        NSData *qualAll = (nReads > 0)
            ? [grun wholeQualitiesData] : [NSData data];
        NSArray<NSString *> *namesAll = [grun allReadNames];
        const uint8_t *seqBytes  = seqAll.bytes;
        const uint8_t *qualBytes = qualAll.bytes;
        NSUInteger qualLenTotal = qualAll.length;
        for (NSUInteger i = 0; i < nReads && emitted < maxAU; i++) {
            uint64_t offset = gIdx ? [gIdx offsetAt:i] : 0;
            uint32_t length = gIdx ? [gIdx lengthAt:i] : 0;
            NSData *seqData = (length > 0)
                ? [NSData dataWithBytes:seqBytes + offset length:length]
                : [NSData data];
            NSData *qualData;
            if (qualLenTotal >= offset + length && length > 0) {
                qualData = [NSData dataWithBytes:qualBytes + offset
                                          length:length];
            } else {
                qualData = [NSData data];
            }
            NSError *readErr = nil;
            TTIOAlignedRead *r = [grun readAtIndex:i error:&readErr];
            if (!r) continue;
            NSString *nameStr = (i < namesAll.count)
                ? namesAll[i] : (r.readName ?: @"");
            TTIOAccessUnit *au =
                accessUnitFromGenomicRead(grun, i, r, seqData, qualData,
                                          nameStr, seqCodec, qualCodec);
            if (filter && ![filter matches:au datasetId:did]) continue;
            if (hasAccessUnitVisitor) {
                [visitor walker:self
                visitAccessUnit:au
                      datasetId:did
                     auSequence:(uint32_t)i];
            }
            emitted++;
        }
        did++;
    }

    // ── 4. EndOfDataset per dataset ────────────────────────────────
    did = 1;
    BOOL hasEodVisitor = [visitor respondsToSelector:
        @selector(walker:visitEndOfDatasetWithDatasetId:finalAUSequence:)];
    for (NSString *name in msNames) {
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOAcquisitionRun *run = dataset.msRuns[name];
        if (hasEodVisitor) {
            [visitor walker:self
visitEndOfDatasetWithDatasetId:did
              finalAUSequence:(uint32_t)[run count]];
        }
        did++;
    }
    for (NSString *name in gNames) {
        if (filterDid && did != filterDid.unsignedIntValue) {
            did++;
            continue;
        }
        TTIOGenomicRun *grun = dataset.genomicRuns[name];
        if (hasEodVisitor) {
            [visitor walker:self
visitEndOfDatasetWithDatasetId:did
              finalAUSequence:(uint32_t)grun.readCount];
        }
        did++;
    }

    // ── 5. EndOfStream ────────────────────────────────────────────
    if ([visitor respondsToSelector:@selector(walkerVisitEndOfStream:)]) {
        [visitor walkerVisitEndOfStream:self];
    }

    return YES;
}

@end
