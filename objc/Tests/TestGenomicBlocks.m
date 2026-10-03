/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * The storage-path genomic writer over a memory-provider root, and the
 * block encoder built on it (format-spec 10.12).
 */
#import <Foundation/Foundation.h>
#import "Testing.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicWriteContext.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Import/TTIOBamReader.h"
#import "Providers/TTIOMemoryProvider.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "ValueClasses/TTIOEnums.h"
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Codecs/TTIOSeqCm.h"
#include <unistd.h>

static NSString *gbBamPath(void)
{
    NSString *rel = @"Tests/Fixtures/genomic/m87_test.bam";
    if ([[NSFileManager defaultManager] fileExistsAtPath:rel]) return rel;
    NSString *abs = @"/home/toddw/TTI-O/objc/Tests/Fixtures/genomic/m87_test.bam";
    return [[NSFileManager defaultManager] fileExistsAtPath:abs] ? abs : nil;
}

TTIOWrittenGenomicRun *gbM87Run(NSString *region)
{
    NSString *path = gbBamPath();
    if (!path) return nil;
    NSError *err = nil;
    TTIOBamReader *r = [[TTIOBamReader alloc] initWithPath:path];
    return [r toGenomicRunWithName:@"genomic_0001" region:region sampleName:nil error:&err];
}

id<TTIOStorageGroup> gbMemRoot(NSString *url)
{
    [TTIOMemoryProvider discardStore:url];
    NSError *err = nil;
    id<TTIOStorageProvider> p = [[TTIOProviderRegistry sharedRegistry]
        openURL:url mode:TTIOStorageOpenModeCreate provider:@"memory" error:&err];
    return [p rootGroupWithError:&err];
}

static void gbStoragePathRoundTrip(void)
{
    // chr1 only: the unmapped zero-length read r005 trips the fqzcomp
    // zero-length decode bug, which is out of scope here.
    TTIOWrittenGenomicRun *run = gbM87Run(@"chr1");
    if (!run) { PASS(YES, "genomic blocks: m87 BAM unavailable, skipped"); return; }
    NSString *url = [NSString stringWithFormat:@"memory://gb-rt-%d", (int)getpid()];
    id<TTIOStorageGroup> root = gbMemRoot(url);
    NSError *err = nil;
    // Force FQZCOMP on qualities: the storage path must take the same
    // codec the HDF5 path takes.
    TTIOWrittenGenomicRun *fq = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:run.acquisitionMode referenceUri:run.referenceUri
                       platform:run.platform sampleName:run.sampleName
                      positions:run.positionsData mappingQualities:run.mappingQualitiesData
                          flags:run.flagsData sequences:run.sequencesData
                      qualities:run.qualitiesData offsets:run.offsetsData
                        lengths:run.lengthsData cigars:run.cigars readNames:run.readNames
                mateChromosomes:run.mateChromosomes matePositions:run.matePositionsData
                templateLengths:run.templateLengthsData chromosomes:run.chromosomes
              signalCompression:run.signalCompression
           signalCodecOverrides:@{@"qualities": @(TTIOCompressionFqzcompNx16Z),
                                  @"cigars": @(TTIOCompressionRansOrder0)}];
    BOOL ok = [TTIOSpectralDataset writeGenomicRunStorage:fq toGroup:root name:@"r"
                                                  context:[TTIOGenomicWriteContext none]
                                                    error:&err];
    PASS(ok, "genomic blocks: storage-path write with FQZCOMP qualities (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    id<TTIOStorageGroup> rg = [root openGroupNamed:@"r" error:&err];
    id<TTIOStorageGroup> sc = [rg openGroupNamed:@"signal_channels" error:&err];
    id<TTIOStorageDataset> q = [sc openDatasetNamed:@"qualities" error:&err];
    id comp = [q attributeValueForName:@"compression" error:NULL];
    PASS([comp integerValue] == TTIOCompressionFqzcompNx16Z,
         "genomic blocks: qualities @compression is FQZCOMP (%ld)", (long)[comp integerValue]);
    TTIOGenomicRun *g = [TTIOGenomicRun openFromGroup:rg name:@"r" error:&err];
    PASS(g != nil, "genomic blocks: memory run opens");
    if (!g) return;
    PASS(g.readCount == run.readCount, "genomic blocks: readCount %lu", (unsigned long)g.readCount);
    BOOL allEqual = YES;
    for (NSUInteger i = 0; i < run.readCount && allEqual; i++) {
        err = nil;
        TTIOAlignedRead *a = [g readAtIndex:i error:&err];
        if (!a) { allEqual = NO; NSLog(@"readAtIndex %lu failed: %@", (unsigned long)i, err); break; }
        const uint64_t *offs = (const uint64_t *)run.offsetsData.bytes;
        const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
        NSData *seq = [run.sequencesData subdataWithRange:NSMakeRange((NSUInteger)offs[i], lens[i])];
        NSData *qual = [run.qualitiesData subdataWithRange:NSMakeRange((NSUInteger)offs[i], lens[i])];
        NSString *seqStr = [[NSString alloc] initWithData:seq encoding:NSASCIIStringEncoding];
        if (![a.readName isEqualToString:run.readNames[i]] || ![a.cigar isEqualToString:run.cigars[i]]
            || ![a.chromosome isEqualToString:run.chromosomes[i]]
            || ![a.sequence isEqualToString:seqStr] || ![a.qualities isEqualToData:qual]) {
            allEqual = NO;
            NSLog(@"read %lu differs: %@ vs %@ / %@ vs %@", (unsigned long)i, a.readName,
                  run.readNames[i], a.sequence, seqStr);
        }
    }
    PASS(allEqual, "genomic blocks: every read reads back through the storage path");
    [TTIOMemoryProvider discardStore:url];
}

static void gbSharedChromMap(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    NSString *url = [NSString stringWithFormat:@"memory://gb-shared-%d", (int)getpid()];
    id<TTIOStorageGroup> root = gbMemRoot(url);
    NSError *err = nil;
    NSMutableDictionary *shared = [NSMutableDictionary dictionaryWithObject:@0 forKey:@"chrZ"];
    TTIOGenomicWriteContext *ctx = [TTIOGenomicWriteContext contextWithChromNameToId:shared referenceMD5:nil];
    BOOL ok = [TTIOSpectralDataset writeGenomicRunStorage:run toGroup:root name:@"a" context:ctx error:&err];
    PASS(ok, "genomic blocks: write with a shared chromosome map");
    NSString *first = run.chromosomes[0];
    PASS([shared[first] unsignedIntegerValue] >= 1, "genomic blocks: pre-seeded id 0 survives");
    id<TTIOStorageGroup> idx = [[root openGroupNamed:@"a" error:&err] openGroupNamed:@"genomic_index" error:&err];
    NSData *ids = [[idx openDatasetNamed:@"chromosome_ids" error:&err] readAll:&err];
    const uint16_t *idv = (const uint16_t *)ids.bytes;
    PASS(ids.length >= 2 && idv[0] == [shared[first] unsignedShortValue],
         "genomic blocks: chromosome_ids use the shared map");
    NSArray *rows = [[idx openDatasetNamed:@"chromosome_names" error:&err] readRows:&err];
    PASS(rows.count > 0 && [[rows[0][@"name"] description] isEqualToString:@"chrZ"],
         "genomic blocks: name table in id order, chrZ first");
    [TTIOMemoryProvider discardStore:url];
}

static void gbSliceConcat(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    NSUInteger n = run.readCount;
    TTIOWrittenGenomicRun *a = [TTIOGenomicBlocks sliceRun:run from:0 to:4];
    TTIOWrittenGenomicRun *b = [TTIOGenomicBlocks sliceRun:run from:4 to:n];
    PASS(a.readCount == 4 && b.readCount == n - 4, "genomic blocks: slice read counts");
    const uint64_t *bo = (const uint64_t *)b.offsetsData.bytes;
    PASS(b.readCount > 0 && bo[0] == 0, "genomic blocks: slice offsets rebased to 0");
    PASS(a.sequencesData.length + b.sequencesData.length == run.sequencesData.length,
         "genomic blocks: slice bases partition the run");
    TTIOWrittenGenomicRun *back = [TTIOGenomicBlocks concatRuns:@[a, b]];
    PASS([back.offsetsData isEqualToData:run.offsetsData]
         && [back.lengthsData isEqualToData:run.lengthsData]
         && [back.sequencesData isEqualToData:run.sequencesData]
         && [back.qualitiesData isEqualToData:run.qualitiesData]
         && [back.positionsData isEqualToData:run.positionsData]
         && [back.flagsData isEqualToData:run.flagsData]
         && [back.mappingQualitiesData isEqualToData:run.mappingQualitiesData]
         && [back.matePositionsData isEqualToData:run.matePositionsData]
         && [back.templateLengthsData isEqualToData:run.templateLengthsData]
         && [back.cigars isEqualToArray:run.cigars]
         && [back.readNames isEqualToArray:run.readNames]
         && [back.mateChromosomes isEqualToArray:run.mateChromosomes]
         && [back.chromosomes isEqualToArray:run.chromosomes],
         "genomic blocks: concat is the inverse of slice");
    TTIOWrittenGenomicRun *empty = [TTIOGenomicBlocks sliceRun:run from:3 to:3];
    PASS(empty.readCount == 0 && empty.sequencesData.length == 0, "genomic blocks: empty slice");
}

static void gbEncodeBlock(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(@"chr1");
    if (!run) return;
    NSError *err = nil;
    TTIOBlockBlobs *bb = [TTIOGenomicBlocks encodeBlock:run context:[TTIOGenomicWriteContext none] error:&err];
    PASS(bb != nil, "genomic blocks: encodeBlock (%s)", [[err localizedDescription] UTF8String] ?: "");
    if (!bb) return;
    PASS(bb.nReads == run.readCount && bb.nBases == run.sequencesData.length,
         "genomic blocks: block counts %lu reads %llu bases",
         (unsigned long)bb.nReads, (unsigned long long)bb.nBases);
    // No reference: SEQ_CM when libttio_rans is linked (M103), else RANS_ORDER1.
    TTIOCompression expSeq = [TTIOSeqCm nativeAvailable] ? TTIOCompressionSeqCm
                                                         : TTIOCompressionRansOrder1;
    PASS([bb.codecs[@"qualities"] integerValue] == TTIOCompressionFqzcompNx16Z
         && [bb.codecs[@"cigars"] integerValue] == TTIOCompressionRansOrder0
         && [bb.codecs[@"sequences"] integerValue] == (NSInteger)expSeq,
         "genomic blocks: forced codecs (q=%ld c=%ld s=%ld)",
         (long)[bb.codecs[@"qualities"] integerValue], (long)[bb.codecs[@"cigars"] integerValue],
         (long)[bb.codecs[@"sequences"] integerValue]);
    PASS([bb.blobs[@"read_names"] length] > 0 && [bb.blobs[@"mate_info"] length] > 0,
         "genomic blocks: read_names and mate_info blobs present");
    // Byte parity with a storage-path write of the same reads under the
    // same overrides.
    NSString *url = [NSString stringWithFormat:@"memory://gb-enc-%d", (int)getpid()];
    id<TTIOStorageGroup> root = gbMemRoot(url);
    TTIOWrittenGenomicRun *same = [run copyWithSignalCodecOverrides:@{
        @"qualities": @(TTIOCompressionFqzcompNx16Z),
        @"cigars": @(TTIOCompressionRansOrder0),
        @"sequences": @(expSeq)}];
    BOOL ok = [TTIOSpectralDataset writeGenomicRunStorage:same toGroup:root name:@"r"
                                                  context:[TTIOGenomicWriteContext none] error:&err];
    PASS(ok, "genomic blocks: reference storage-path write");
    id<TTIOStorageGroup> sc = [[root openGroupNamed:@"r" error:&err] openGroupNamed:@"signal_channels" error:&err];
    BOOL parity = YES;
    for (NSString *ch in @[@"sequences", @"qualities", @"read_names", @"cigars"]) {
        NSData *ref = [[sc openDatasetNamed:ch error:&err] readAll:&err];
        if (![ref isEqualToData:bb.blobs[ch]]) { parity = NO; NSLog(@"channel %@ differs", ch); }
    }
    NSData *mate = [[[sc openGroupNamed:@"mate_info" error:&err] openDatasetNamed:@"inline_v2" error:&err] readAll:&err];
    if (![mate isEqualToData:bb.blobs[@"mate_info"]]) parity = NO;
    PASS(parity, "genomic blocks: block blobs equal a whole-run write of the same reads");
    [TTIOMemoryProvider discardStore:url];

    // A zero-length read forces RANS_ORDER0 on qualities.
    TTIOWrittenGenomicRun *whole = gbM87Run(nil);
    TTIOBlockBlobs *zb = [TTIOGenomicBlocks encodeBlock:whole context:[TTIOGenomicWriteContext none] error:&err];
    PASS(zb != nil && [zb.codecs[@"qualities"] integerValue] == TTIOCompressionRansOrder0,
         "genomic blocks: zero-length read selects RANS_ORDER0 qualities");
}

static void gbRunCopies(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    run.optDisableQualitiesV5 = YES;
    run.embedReference = YES;
    TTIOWrittenGenomicRun *c = [run copyWithOptLegacyWholeChannel:YES];
    PASS(c.optLegacyWholeChannel && !run.optLegacyWholeChannel && c.optDisableQualitiesV5
         && c.embedReference && c.readCount == run.readCount,
         "genomic blocks: copyWithOptLegacyWholeChannel carries the options");
    TTIOWrittenGenomicRun *p = [run copyWithProvenance:@[]];
    PASS(p.provenanceRecords.count == 0 && [p.signalCodecOverrides isEqualToDictionary:run.signalCodecOverrides],
         "genomic blocks: copyWithProvenance keeps the overrides");
}

static void gbQualStrategyHint(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(@"chr1");
    if (!run) return;
    NSError *err = nil;
    TTIOBlockBlobs *def = [TTIOGenomicBlocks encodeBlock:run
        context:[TTIOGenomicWriteContext none] error:&err];
    PASS(def != nil
         && [TTIOFqzcompNx16Z strategyOfEncodedStream:def.blobs[@"qualities"]] == 4,
         "genomic blocks: default qualities strategy is V4 on the small block");
    TTIOGenomicWriteContext *c5 = [TTIOGenomicWriteContext none];
    c5.qualStrategyHint = 5;
    TTIOBlockBlobs *forced = [TTIOGenomicBlocks encodeBlock:run
        context:c5 error:&err];
    PASS(forced != nil
         && [TTIOFqzcompNx16Z strategyOfEncodedStream:forced.blobs[@"qualities"]] == 5,
         "genomic blocks: context hint 5 forces a V5-S5 qualities stream");
    TTIOGenomicWriteContext *cm1 = [TTIOGenomicWriteContext none];
    cm1.qualStrategyHint = -1;
    TTIOBlockBlobs *expl = [TTIOGenomicBlocks encodeBlock:run
        context:cm1 error:&err];
    PASS(expl != nil
         && [expl.blobs[@"qualities"] isEqualToData:def.blobs[@"qualities"]],
         "genomic blocks: explicit hint -1 is the default path byte for byte");
}

void testGenomicBlocks(void)
{
    @autoreleasepool {
        gbStoragePathRoundTrip();
        gbSharedChromMap();
        gbSliceConcat();
        gbEncodeBlock();
        gbQualStrategyHint();
        gbRunCopies();
    }
}
