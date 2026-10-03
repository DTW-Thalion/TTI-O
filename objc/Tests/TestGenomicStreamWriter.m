/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 *
 * TTIOGenomicStreamWriter: the blocks_v1 layout (format-spec 10.12),
 * chromosome cuts, the legacy flag, single-read appends and the default
 * flip in writeMinimalToPath.
 */
#import <Foundation/Foundation.h>
#import "Testing.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"
#import "Genomics/TTIOGenomicStreamWriter.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOLazyReference.h"
#import "Codecs/TTIOFqzcompNx16Z.h"
#import "Providers/TTIOMemoryProvider.h"
#import "Providers/TTIOProviderRegistry.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Providers/TTIOCompoundField.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "ValueClasses/TTIOEnums.h"
#import "Codecs/TTIOSeqCm.h"
#include <unistd.h>

extern TTIOWrittenGenomicRun *gbM87Run(NSString *region);
extern id<TTIOStorageGroup> gbMemRoot(NSString *url);

static id<TTIOStorageGroup> gswStudy(NSString *url)
{
    id<TTIOStorageGroup> root = gbMemRoot(url);
    return [root createGroupNamed:@"study" error:NULL];
}

static NSArray *gswIndexRows(id<TTIOStorageGroup> study, NSString *name)
{
    NSError *err = nil;
    id<TTIOStorageGroup> rg = [[study openGroupNamed:@"genomic_runs" error:&err] openGroupNamed:name error:&err];
    id<TTIOStorageDataset> idx = [[rg openGroupNamed:@"blocks" error:&err] openDatasetNamed:@"index" error:&err];
    return [idx readRows:&err];
}

static void gswLayout(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) { PASS(YES, "stream writer: m87 BAM unavailable, skipped"); return; }
    NSString *url = [NSString stringWithFormat:@"memory://gsw-layout-%d", (int)getpid()];
    id<TTIOStorageGroup> study = gswStudy(url);
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = 3;
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:study runName:@"g" options:o];
    PASS([w appendBatch:run error:&err], "stream writer: appendBatch (%s)", [[err localizedDescription] UTF8String] ?: "");
    PASS([w close:&err], "stream writer: close (%s)", [[err localizedDescription] UTF8String] ?: "");
    PASS(w.readCount == 10, "stream writer: readCount %llu", w.readCount);
    // chr1 5 reads -> 3+2, chr2 2 -> 2, * 3 -> 3
    PASS(w.blockCount == 4, "stream writer: chromosome cuts give 4 blocks (%lu)", (unsigned long)w.blockCount);
    id<TTIOStorageGroup> runs = [study openGroupNamed:@"genomic_runs" error:&err];
    PASS([[[runs attributeValueForName:@"_run_names" error:NULL] description] isEqualToString:@"g"],
         "stream writer: _run_names maintained");
    id<TTIOStorageGroup> rg = [runs openGroupNamed:@"g" error:&err];
    PASS([[[rg attributeValueForName:@"layout" error:NULL] description] isEqualToString:@"blocks_v1"],
         "stream writer: @layout blocks_v1");
    PASS([[[rg attributeValueForName:@"block_policy" error:NULL] description]
             isEqualToString:@"reads=3,bytes=67108864"], "stream writer: @block_policy");
    PASS([[rg attributeValueForName:@"read_count" error:NULL] longLongValue] == 10
         && [[rg attributeValueForName:@"base_count" error:NULL] longLongValue] == 720,
         "stream writer: @read_count/@base_count");
    NSArray *rows = gswIndexRows(study, @"g");
    PASS(rows.count == 4, "stream writer: 4 index rows (%lu)", (unsigned long)rows.count);
    PASS([TTIOGenomicStreamWriter indexFields].count == 19 && [rows[0] count] == 19,
         "stream writer: 19 index columns");
    NSDictionary *r0 = rows[0], *r1 = rows.count > 1 ? rows[1] : nil;
    PASS([r0[@"read_start"] unsignedLongLongValue] == 0 && [r0[@"n_reads"] unsignedIntValue] == 3
         && [r1[@"read_start"] unsignedLongLongValue] == 3 && [r1[@"n_reads"] unsignedIntValue] == 2
         && [r1[@"base_start"] unsignedLongLongValue] == 300,
         "stream writer: read_start/base_start chain");
    unsigned long long off1 = [r1[@"qualities_off"] unsignedLongLongValue];
    PASS(off1 == [r0[@"qualities_len"] unsignedLongLongValue],
         "stream writer: block 1 qualities_off follows block 0 qualities_len");
    PASS([r0[@"qualities_codec"] unsignedIntValue] == TTIOCompressionFqzcompNx16Z
         && [r0[@"cigars_codec"] unsignedIntValue] == TTIOCompressionRansOrder0
         && [r0[@"sequences_codec"] unsignedIntValue]
                == ([TTIOSeqCm nativeAvailable] ? TTIOCompressionSeqCm : TTIOCompressionRansOrder1),
         "stream writer: forced codecs in the index");
    NSDictionary *r3 = rows.count > 3 ? rows[3] : nil;
    PASS([r3[@"qualities_codec"] unsignedIntValue] == TTIOCompressionRansOrder0,
         "stream writer: unmapped block with a zero-length read uses RANS_ORDER0 qualities");
    id<TTIOStorageGroup> sc = [rg openGroupNamed:@"signal_channels" error:&err];
    id<TTIOStorageDataset> seq = [[sc openGroupNamed:@"sequences" error:&err] openDatasetNamed:@"data" error:&err];
    PASS(seq != nil && [seq isExtendable], "stream writer: sequences/data extendable");
    unsigned long long seqTotal = 0;
    for (NSDictionary *r in rows) seqTotal += [r[@"sequences_len"] unsignedLongLongValue];
    PASS([seq length] == seqTotal, "stream writer: sequences/data length is the sum of the block lens");
    id<TTIOStorageDataset> q = [sc openDatasetNamed:@"qualities" error:&err];
    PASS([[q attributeValueForName:@"compression" error:NULL] integerValue] == TTIOCompressionFqzcompNx16Z,
         "stream writer: qualities @compression is the first block's codec");
    id<TTIOStorageGroup> idx = [rg openGroupNamed:@"genomic_index" error:&err];
    PASS([[[idx openDatasetNamed:@"lengths" error:&err] readAll:&err] isEqualToData:run.lengthsData]
         && [[[idx openDatasetNamed:@"flags" error:&err] readAll:&err] isEqualToData:run.flagsData]
         && [[[idx openDatasetNamed:@"positions" error:&err] readAll:&err] isEqualToData:run.positionsData],
         "stream writer: index arrays equal the run's");
    PASS(![idx hasChildNamed:@"offsets"], "stream writer: no offsets on disk");
    NSArray *names = [[idx openDatasetNamed:@"chromosome_names" error:&err] readRows:&err];
    PASS(names.count == 3 && [[names[0][@"name"] description] isEqualToString:@"chr1"],
         "stream writer: chromosome_names in first-seen order (%lu)", (unsigned long)names.count);
    NSArray *mateNames = [[[sc openGroupNamed:@"mate_info" error:&err] openDatasetNamed:@"chrom_names" error:&err] readRows:&err];
    PASS(mateNames.count == names.count, "stream writer: mate_info/chrom_names shares the map");
    [TTIOMemoryProvider discardStore:url];
}

static void gswLegacy(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    NSString *url = [NSString stringWithFormat:@"memory://gsw-legacy-%d", (int)getpid()];
    id<TTIOStorageGroup> study = gswStudy(url);
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.optLegacyWholeChannel = YES;
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:study runName:@"g" options:o];
    PASS([w appendBatch:run error:&err] && [w close:&err], "stream writer: legacy write");
    id<TTIOStorageGroup> rg = [[study openGroupNamed:@"genomic_runs" error:&err] openGroupNamed:@"g" error:&err];
    PASS(rg != nil && ![rg hasAttributeNamed:@"layout"] && ![rg hasChildNamed:@"blocks"],
         "stream writer: legacy flag gives the whole-channel layout");
    PASS([[[rg openGroupNamed:@"signal_channels" error:&err] openDatasetNamed:@"cigars" error:&err] compoundFields].count == 1,
         "stream writer: legacy layout keeps the compound cigars");
    [TTIOMemoryProvider discardStore:url];
}

static void gswSingleReadEqualsBatch(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    NSError *err = nil;
    NSString *urlA = [NSString stringWithFormat:@"memory://gsw-a-%d", (int)getpid()];
    NSString *urlB = [NSString stringWithFormat:@"memory://gsw-b-%d", (int)getpid()];
    id<TTIOStorageGroup> sa = gswStudy(urlA), sb = gswStudy(urlB);
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = 4;
    TTIOGenomicStreamWriter *wa = [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:sa runName:@"g" options:o];
    TTIOGenomicStreamWriter *wb = [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:sb runName:@"g" options:o];
    BOOL ok = [wa appendBatch:run error:&err] && [wa close:&err];
    const uint64_t *offs = (const uint64_t *)run.offsetsData.bytes;
    const uint32_t *lens = (const uint32_t *)run.lengthsData.bytes;
    const int64_t *pos = (const int64_t *)run.positionsData.bytes;
    const uint8_t *mapq = (const uint8_t *)run.mappingQualitiesData.bytes;
    const uint32_t *flags = (const uint32_t *)run.flagsData.bytes;
    const int64_t *mpos = (const int64_t *)run.matePositionsData.bytes;
    const int32_t *tlen = (const int32_t *)run.templateLengthsData.bytes;
    for (NSUInteger i = 0; i < run.readCount && ok; i++) {
        NSData *seq = [run.sequencesData subdataWithRange:NSMakeRange((NSUInteger)offs[i], lens[i])];
        NSData *qual = [run.qualitiesData subdataWithRange:NSMakeRange((NSUInteger)offs[i], lens[i])];
        TTIOAlignedRead *r = [[TTIOAlignedRead alloc]
            initWithReadName:run.readNames[i] chromosome:run.chromosomes[i] position:pos[i]
              mappingQuality:mapq[i] cigar:run.cigars[i]
                    sequence:[[NSString alloc] initWithData:seq encoding:NSASCIIStringEncoding]
                   qualities:qual flags:flags[i] mateChromosome:run.mateChromosomes[i]
                matePosition:mpos[i] templateLength:tlen[i]];
        ok = [wb appendRead:r error:&err];
    }
    ok = ok && [wb close:&err];
    PASS(ok, "stream writer: per-read appends (%s)", [[err localizedDescription] UTF8String] ?: "");
    NSArray *ra = gswIndexRows(sa, @"g"), *rb = gswIndexRows(sb, @"g");
    PASS(ra.count == rb.count && ra.count == 4, "stream writer: same block count (%lu vs %lu)",
         (unsigned long)ra.count, (unsigned long)rb.count);
    BOOL same = ra.count == rb.count;
    for (NSUInteger i = 0; i < ra.count && same; i++) {
        for (TTIOCompoundField *f in [TTIOGenomicStreamWriter indexFields]) {
            if (![[ra[i][f.name] description] isEqualToString:[rb[i][f.name] description]]) { same = NO; break; }
        }
    }
    PASS(same, "stream writer: identical block index");
    id<TTIOStorageGroup> sca = [[[sa openGroupNamed:@"genomic_runs" error:&err] openGroupNamed:@"g" error:&err] openGroupNamed:@"signal_channels" error:&err];
    id<TTIOStorageGroup> scb = [[[sb openGroupNamed:@"genomic_runs" error:&err] openGroupNamed:@"g" error:&err] openGroupNamed:@"signal_channels" error:&err];
    BOOL bytesSame = YES;
    for (NSString *ch in @[@"qualities", @"read_names", @"cigars"]) {
        NSData *a = [[sca openDatasetNamed:ch error:&err] readAll:&err];
        NSData *b = [[scb openDatasetNamed:ch error:&err] readAll:&err];
        if (![a isEqualToData:b]) bytesSame = NO;
    }
    NSData *qa = [[[sca openGroupNamed:@"sequences" error:&err] openDatasetNamed:@"data" error:&err] readAll:&err];
    NSData *qb = [[[scb openGroupNamed:@"sequences" error:&err] openDatasetNamed:@"data" error:&err] readAll:&err];
    if (![qa isEqualToData:qb]) bytesSame = NO;
    PASS(bytesSame, "stream writer: identical channel bytes");
    [TTIOMemoryProvider discardStore:urlA];
    [TTIOMemoryProvider discardStore:urlB];
}

static void gswWriteMinimalDefault(void)
{
    TTIOWrittenGenomicRun *run = gbM87Run(nil);
    if (!run) return;
    NSError *err = nil;
    NSString *path = [NSString stringWithFormat:@"/tmp/gsw-default-%d.tio", (int)getpid()];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
    BOOL ok = [TTIOSpectralDataset writeMinimalToPath:path title:@"t" isaInvestigationId:@"i"
                                              msRuns:@{} genomicRuns:@{@"genomic_0001": run}
                                     identifications:nil quantifications:nil provenanceRecords:nil
                                               error:&err];
    PASS(ok, "stream writer: writeMinimalToPath default (%s)", [[err localizedDescription] UTF8String] ?: "");
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:&err];
    TTIOHDF5Group *rg = [[[f rootGroup] openGroupNamed:@"study" error:&err] openGroupNamed:@"genomic_runs" error:&err];
    rg = [rg openGroupNamed:@"genomic_0001" error:&err];
    PASS(rg != nil && [[rg stringAttributeNamed:@"layout" error:NULL] isEqualToString:@"blocks_v1"]
         && [rg hasChildNamed:@"blocks"],
         "stream writer: writeMinimalToPath writes blocks_v1 by default");
    [f close];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];

    TTIOWrittenGenomicRun *legacy = [run copyWithOptLegacyWholeChannel:YES];
    ok = [TTIOSpectralDataset writeMinimalToPath:path title:@"t" isaInvestigationId:@"i"
                                          msRuns:@{} genomicRuns:@{@"genomic_0001": legacy}
                                 identifications:nil quantifications:nil provenanceRecords:nil
                                           error:&err];
    f = [TTIOHDF5File openReadOnlyAtPath:path error:&err];
    rg = [[[[f rootGroup] openGroupNamed:@"study" error:&err] openGroupNamed:@"genomic_runs" error:&err]
          openGroupNamed:@"genomic_0001" error:&err];
    PASS(ok && rg != nil && ![rg hasAttributeNamed:@"layout"] && ![rg hasChildNamed:@"blocks"],
         "stream writer: optLegacyWholeChannel keeps the v1.8 layout in writeMinimalToPath");
    [f close];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];

    // Storage-path writeMinimal (memory://) takes the same default.
    NSString *url = [NSString stringWithFormat:@"memory://gsw-wm-%d", (int)getpid()];
    [TTIOMemoryProvider discardStore:url];
    ok = [TTIOSpectralDataset writeMinimalToPath:url title:@"t" isaInvestigationId:@"i"
                                          msRuns:@{} genomicRuns:@{@"genomic_0001": run}
                                 identifications:nil quantifications:nil provenanceRecords:nil
                                           error:&err];
    PASS(ok, "stream writer: memory:// writeMinimal (%s)", [[err localizedDescription] UTF8String] ?: "");
    id<TTIOStorageProvider> p = [[TTIOProviderRegistry sharedRegistry] openURL:url mode:TTIOStorageOpenModeRead provider:@"memory" error:&err];
    id<TTIOStorageGroup> mrg = [[[[p rootGroupWithError:&err] openGroupNamed:@"study" error:&err]
                                 openGroupNamed:@"genomic_runs" error:&err] openGroupNamed:@"genomic_0001" error:&err];
    PASS([[[mrg attributeValueForName:@"layout" error:NULL] description] isEqualToString:@"blocks_v1"],
         "stream writer: memory:// writeMinimal writes blocks_v1");
    [TTIOMemoryProvider discardStore:url];
}

static void gswLazyReference(void)
{
    NSError *err = nil;
    NSString *fa = [NSString stringWithFormat:@"/tmp/gsw-ref-%d.fa", (int)getpid()];
    NSString *text = @">chrA desc\nACGTACGTAC\nGTAC\n>chrB\nNNNNACGT\nACG\n>chrE\n";
    [text writeToFile:fa atomically:YES encoding:NSASCIIStringEncoding error:NULL];
    [[NSFileManager defaultManager] removeItemAtPath:[fa stringByAppendingString:@".fai"] error:NULL];
    TTIOLazyReference *ref = [[TTIOLazyReference alloc] initWithFastaPath:fa error:&err];
    PASS(ref != nil, "lazy reference: opens without a .fai (%s)", [[err localizedDescription] UTF8String] ?: "");
    if (!ref) return;
    PASS([[NSFileManager defaultManager] fileExistsAtPath:[fa stringByAppendingString:@".fai"]],
         "lazy reference: writes the .fai");
    NSArray *expectNames = @[@"chrA", @"chrB", @"chrE"];
    PASS(ref.count == 3 && [ref.chromosomeNames isEqualToArray:expectNames],
         "lazy reference: three entries in FASTA order");
    PASS([ref lengthOf:@"chrA"] == 14 && [ref lengthOf:@"chrB"] == 11 && [ref lengthOf:@"chrE"] == 0,
         "lazy reference: lengths from the index");
    NSData *a = ref[@"chrA"];
    PASS([a isEqualToData:[@"ACGTACGTACGTAC" dataUsingEncoding:NSASCIIStringEncoding]],
         "lazy reference: chrA bytes without newlines");
    PASS([ref[@"chrB"] isEqualToData:[@"NNNNACGTACG" dataUsingEncoding:NSASCIIStringEncoding]],
         "lazy reference: chrB across uneven lines");
    PASS(ref[@"chrE"].length == 0 && ref[@"nope"] == nil, "lazy reference: empty and unknown");
    PASS([[ref copy] isKindOfClass:[TTIOLazyReference class]], "lazy reference: copy is the same object");
    NSMutableArray *keys = [NSMutableArray array];
    for (NSString *k in ref) [keys addObject:k];
    PASS(keys.count == 3, "lazy reference: fast enumeration over keys");
    // A second open reads the .fai written by the first.
    TTIOLazyReference *ref2 = [[TTIOLazyReference alloc] initWithFastaPath:fa error:&err];
    PASS(ref2 != nil && [ref2[@"chrA"] isEqualToData:a], "lazy reference: reopens from the .fai");
    [[NSFileManager defaultManager] removeItemAtPath:fa error:NULL];
    [[NSFileManager defaultManager] removeItemAtPath:[fa stringByAppendingString:@".fai"] error:NULL];
}

void testGenomicStreamWriter(void)
{
    @autoreleasepool {
        gswLayout();
        gswLegacy();
        gswSingleReadEqualsBatch();
        gswWriteMinimalDefault();
        gswLazyReference();
    }
}

// ── block-parallel writer ─────────────────────────────────────────────

/* Two chromosomes, placed-unmapped reads (every 97th, cigar "*"),
 * cross-chromosome mates (every 13th), "=" mates (every 3rd); mirrors the
 * Python and Java tests. */
static int gswCmpInt64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static TTIOWrittenGenomicRun *gswBigSyntheticRun(NSUInteger n, unsigned seed)
{
    const NSUInteger L = 100;
    srand(seed);
    const char *alphabet = "ACGT";
    NSMutableData *ref1 = [NSMutableData dataWithLength:400000];
    NSMutableData *ref2 = [NSMutableData dataWithLength:400000];
    uint8_t *r1 = ref1.mutableBytes, *r2 = ref2.mutableBytes;
    for (NSUInteger i = 0; i < 400000; i++) { r1[i] = alphabet[rand() % 4]; r2[i] = alphabet[rand() % 4]; }
    NSDictionary *refs = @{@"chr1": ref1, @"chr2": ref2};
    NSUInteger half = n / 2;
    NSMutableData *posD = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    int64_t *pos = posD.mutableBytes;
    for (NSUInteger i = 0; i < half; i++) pos[i] = 1 + rand() % 399000;
    for (NSUInteger i = half; i < n; i++) pos[i] = 1 + rand() % 399000;
    qsort(pos, half, sizeof(int64_t), gswCmpInt64);
    qsort(pos + half, n - half, sizeof(int64_t), gswCmpInt64);
    NSMutableData *seqD = [NSMutableData dataWithLength:n * L];
    NSMutableData *qualD = [NSMutableData dataWithLength:n * L];
    NSMutableData *mqD = [NSMutableData dataWithLength:n];
    NSMutableData *flD = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *offD = [NSMutableData dataWithLength:n * sizeof(uint64_t)];
    NSMutableData *lenD = [NSMutableData dataWithLength:n * sizeof(uint32_t)];
    NSMutableData *mpD = [NSMutableData dataWithLength:n * sizeof(int64_t)];
    NSMutableData *tlD = [NSMutableData dataWithLength:n * sizeof(int32_t)];
    uint8_t *seq = seqD.mutableBytes, *qual = qualD.mutableBytes, *mq = mqD.mutableBytes;
    uint32_t *fl = flD.mutableBytes, *len = lenD.mutableBytes;
    uint64_t *off = offD.mutableBytes;
    int64_t *mp = mpD.mutableBytes;
    NSMutableArray *chroms = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *cigars = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *names = [NSMutableArray arrayWithCapacity:n];
    NSMutableArray *mates = [NSMutableArray arrayWithCapacity:n];
    for (NSUInteger i = 0; i < n; i++) {
        NSString *c = i < half ? @"chr1" : @"chr2";
        [chroms addObject:c];
        const uint8_t *ref = [refs[c] bytes];
        memcpy(seq + i * L, ref + pos[i] - 1, L);
        for (int k = 0; k < 3; k++) seq[i * L + rand() % L] = alphabet[rand() % 4];
        for (NSUInteger k = 0; k < L; k++) qual[i * L + k] = (uint8_t)(2 + rand() % 38);
        off[i] = i * L; len[i] = (uint32_t)L; mq[i] = 60;
        [names addObject:[NSString stringWithFormat:@"r%06lu", (unsigned long)i]];
        fl[i] = 0x3; mp[i] = -1;
        if (i % 97 == 0) { [cigars addObject:@"*"]; fl[i] = 0x5; }
        else [cigars addObject:[NSString stringWithFormat:@"%luM", (unsigned long)L]];
        if (i % 13 == 0) { [mates addObject:[c isEqualToString:@"chr1"] ? @"chr2" : @"chr1"]; mp[i] = pos[(i * 7) % n]; }
        else if (i % 3 == 0) { [mates addObject:@"="]; mp[i] = pos[i] + 200; }
        else [mates addObject:@""];
    }
    TTIOWrittenGenomicRun *run = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:TTIOAcquisitionModeGenomicWGS
                   referenceUri:@"synthetic" platform:@"ILLUMINA" sampleName:@"s"
                      positions:posD mappingQualities:mqD flags:flD
                      sequences:seqD qualities:qualD offsets:offD lengths:lenD
                         cigars:cigars readNames:names mateChromosomes:mates
                  matePositions:mpD templateLengths:tlD chromosomes:chroms
              signalCompression:TTIOCompressionZlib signalCodecOverrides:@{}];
    run.referenceChromSeqs = refs;
    run.embedReference = YES;
    return run;
}

static BOOL gswWriteThreadsToStudy(id<TTIOStorageGroup> study, TTIOWrittenGenomicRun *run,
                                   NSUInteger threads, NSUInteger blockReads)
{
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = blockReads;
    o.threads = threads;
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"g" options:o];
    NSUInteger n = run.readCount;
    for (NSUInteger a = 0; a < n; a += 7001) {
        TTIOWrittenGenomicRun *part = [TTIOGenomicBlocks sliceRun:run from:a to:MIN(n, a + 7001)];
        if (![w appendBatch:part error:&err]) { PASS(NO, "bp: appendBatch (%s)", [[err localizedDescription] UTF8String] ?: ""); return NO; }
    }
    PASS(w.threads == threads, "bp: writer threads resolved to %lu", (unsigned long)threads);
    if (![w close:&err]) { PASS(NO, "bp: close (%s)", [[err localizedDescription] UTF8String] ?: ""); return NO; }
    return YES;
}

static id<TTIOStorageGroup> gswWriteThreads(NSString *url, TTIOWrittenGenomicRun *run,
                                            NSUInteger threads, NSUInteger blockReads)
{
    id<TTIOStorageGroup> study = gswStudy(url);
    return gswWriteThreadsToStudy(study, run, threads, blockReads) ? study : nil;
}

static NSString *gswCanon(id v)
{
    if (v == nil) return @"nil";
    if ([v isKindOfClass:[NSData class]]) {
        NSData *d = (NSData *)v;
        uint8_t md[16] = {0};
        // cheap digest: fold length + bytes
        unsigned long h = 1469598103u ^ d.length;
        const uint8_t *b = d.bytes;
        for (NSUInteger i = 0; i < d.length; i++) h = (h * 1099511628211ul) ^ b[i];
        (void)md;
        return [NSString stringWithFormat:@"data:%lu:%lx", (unsigned long)d.length, h];
    }
    if ([v isKindOfClass:[NSArray class]]) {
        NSMutableString *s = [NSMutableString stringWithString:@"["];
        for (id o in (NSArray *)v) { [s appendString:gswCanon(o)]; [s appendString:@","]; }
        [s appendString:@"]"];
        return s;
    }
    if ([v isKindOfClass:[NSDictionary class]]) {
        NSMutableString *s = [NSMutableString stringWithString:@"{"];
        for (NSString *k in [[(NSDictionary *)v allKeys] sortedArrayUsingSelector:@selector(compare:)]) {
            [s appendFormat:@"%@=%@;", k, gswCanon(((NSDictionary *)v)[k])];
        }
        [s appendString:@"}"];
        return s;
    }
    return [v description];
}

static void gswCollect(id<TTIOStorageGroup> g, NSString *prefix, NSMutableDictionary *out)
{
    NSMutableString *attrs = [NSMutableString string];
    for (NSString *a in [[g attributeNames] sortedArrayUsingSelector:@selector(compare:)]) {
        [attrs appendFormat:@"%@=%@;", a, gswCanon([g attributeValueForName:a error:NULL])];
    }
    out[prefix] = attrs;
    for (NSString *c in [g childNames]) {
        NSError *e = nil;
        id<TTIOStorageGroup> sub = [g openGroupNamed:c error:&e];
        if (sub) {
            gswCollect(sub, [NSString stringWithFormat:@"%@/%@", prefix, c], out);
        } else {
            id<TTIOStorageDataset> ds = [g openDatasetNamed:c error:&e];
            if (!ds) continue;
            NSMutableString *da = [NSMutableString string];
            for (NSString *a in [[ds attributeNames] sortedArrayUsingSelector:@selector(compare:)]) {
                [da appendFormat:@"%@=%@;", a, gswCanon([ds attributeValueForName:a error:NULL])];
            }
            out[[NSString stringWithFormat:@"%@/%@", prefix, c]] =
                [NSString stringWithFormat:@"%@|%@", da, gswCanon([ds readAll:NULL])];
        }
    }
}

static void gswThreadedIdentical(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(60000, 7);
    id<TTIOStorageGroup> a = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-bp-a-%d", (int)getpid()], run, 1, 20000);
    id<TTIOStorageGroup> b = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-bp-b-%d", (int)getpid()], run, 6, 20000);
    if (!a || !b) return;
    NSMutableDictionary *ma = [NSMutableDictionary dictionary], *mb = [NSMutableDictionary dictionary];
    gswCollect(a, @"", ma);
    gswCollect(b, @"", mb);
    PASS([[NSSet setWithArray:ma.allKeys] isEqualToSet:[NSSet setWithArray:mb.allKeys]],
         "bp: same object set (%lu vs %lu)", (unsigned long)ma.count, (unsigned long)mb.count);
    BOOL same = YES;
    NSString *bad = nil;
    for (NSString *k in ma) {
        if (![ma[k] isEqualToString:mb[k] ?: @""]) { same = NO; bad = k; break; }
    }
    PASS(same, "bp: threads=1 and threads=6 files identical%s%s",
         bad ? " first diff at " : "", bad ? [bad UTF8String] : "");
}

static void gswRegisterOrder(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(200, 3);
    NSMutableDictionary *m = [NSMutableDictionary dictionary];
    [TTIOGenomicStreamWriter registerBlockChromosomes:run intoMap:m];
    PASS([m[@"chr1"] intValue] == 0, "bp: own names first (chr1 -> 0)");
    NSUInteger size = m.count;
    [TTIOGenomicStreamWriter registerBlockChromosomes:run intoMap:m];
    PASS(m.count == size, "bp: registration is idempotent");
}

/* One iterReads pass as name|sequence|cigar|mateChromosome lines. */
static NSArray<NSString *> *gswIterStrings(TTIOGenomicRun *g, NSUInteger start, NSUInteger stop,
                                           NSUInteger threads, BOOL *okOut)
{
    NSMutableArray *out = [NSMutableArray array];
    NSError *err = nil;
    BOOL ok = [g iterReadsFrom:start to:stop threads:threads error:&err
                    usingBlock:^(TTIOAlignedRead *r, NSUInteger index, BOOL *stop2) {
        [out addObject:[NSString stringWithFormat:@"%lu:%@|%@|%@|%@", (unsigned long)index,
                        r.readName, r.sequence, r.cigar, r.mateChromosome ?: @""]];
    }];
    if (!ok) PASS(NO, "bp: iterReads threads=%lu (%s)", (unsigned long)threads,
                  [[err localizedDescription] UTF8String] ?: "");
    *okOut = ok;
    return out;
}

static void gswIterReadsThreaded(void)
{
    /* A real HDF5 file: the default REF_DIFF reference resolver is built
     * from the run group's owning file (nil on the memory provider). */
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(30000, 11);
    NSString *path = [NSString stringWithFormat:@"/tmp/gsw-bp-it-%d.tio", (int)getpid()];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
    NSError *err = nil;
    id<TTIOStorageProvider> pw = [[TTIOProviderRegistry sharedRegistry]
        openURL:path mode:TTIOStorageOpenModeCreate provider:nil error:&err];
    id<TTIOStorageGroup> ws = [[pw rootGroupWithError:&err] createGroupNamed:@"study" error:&err];
    if (!gswWriteThreadsToStudy(ws, run, 1, 5000)) { [pw close]; return; }
    [pw close];
    id<TTIOStorageProvider> pr = [[TTIOProviderRegistry sharedRegistry]
        openURL:path mode:TTIOStorageOpenModeRead provider:nil error:&err];
    id<TTIOStorageGroup> study = [[pr rootGroupWithError:&err] openGroupNamed:@"study" error:&err];
    id<TTIOStorageGroup> rg = [[study openGroupNamed:@"genomic_runs" error:&err]
                               openGroupNamed:@"g" error:&err];
    TTIOGenomicRun *g = [TTIOGenomicRun openFromGroup:rg name:@"g" error:&err];
    PASS(g != nil && [g readCount] == 30000, "bp: reader run open (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    if (!g) return;
    BOOL ok1 = NO, ok4 = NO;
    NSArray *serial = gswIterStrings(g, 0, 30000, 1, &ok1);
    NSArray *threaded = gswIterStrings(g, 0, 30000, 4, &ok4);
    PASS(ok1 && ok4 && [serial isEqualToArray:threaded],
         "bp: iterReads threads=4 equals serial (%lu reads)", (unsigned long)threaded.count);
    BOOL okA = NO, okB = NO;
    NSArray *subS = gswIterStrings(g, 12345, 17890, 1, &okA);
    NSArray *subT = gswIterStrings(g, 12345, 17890, 3, &okB);
    PASS(okA && okB && subS.count == 17890 - 12345 && [subS isEqualToArray:subT],
         "bp: iterReads sub-range threads=3 equals serial (%lu reads)", (unsigned long)subT.count);
    [pr close];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
}

/* -iterBlocksFrom:to:threads: is weaker than -iterReadsFrom:to:threads:
 * by exactly one thing, that blocks arrive in no order and on several
 * threads. So what has to hold is that every read is delivered exactly
 * once whatever the thread count, and that the ranges tile the span.
 * Python: test_genomic_for_each_block.py; Java: GenomicRunBlocksTest. */
static void gswIterBlocks(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(30000, 12);
    NSString *path = [NSString stringWithFormat:@"/tmp/gsw-bp-blk-%d.tio", (int)getpid()];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
    NSError *err = nil;
    id<TTIOStorageProvider> pw = [[TTIOProviderRegistry sharedRegistry]
        openURL:path mode:TTIOStorageOpenModeCreate provider:nil error:&err];
    id<TTIOStorageGroup> ws = [[pw rootGroupWithError:&err] createGroupNamed:@"study" error:&err];
    if (!gswWriteThreadsToStudy(ws, run, 1, 5000)) { [pw close]; return; }
    [pw close];
    id<TTIOStorageProvider> pr = [[TTIOProviderRegistry sharedRegistry]
        openURL:path mode:TTIOStorageOpenModeRead provider:nil error:&err];
    id<TTIOStorageGroup> study = [[pr rootGroupWithError:&err] openGroupNamed:@"study" error:&err];
    id<TTIOStorageGroup> rg = [[study openGroupNamed:@"genomic_runs" error:&err]
                               openGroupNamed:@"g" error:&err];
    TTIOGenomicRun *g = [TTIOGenomicRun openFromGroup:rg name:@"g" error:&err];
    if (!g) { PASS(NO, "bp blocks: reader run open"); [pr close]; return; }

    NSUInteger threadCounts[] = { 1, 2, 3, 8 };
    for (int ti = 0; ti < 4; ti++) {
        NSLock *lock = [NSLock new];
        NSMutableIndexSet *seen = [NSMutableIndexSet indexSet];
        __block NSUInteger delivered = 0, covered = 0;
        __block BOOL contiguous = YES;
        NSMutableArray *ranges = [NSMutableArray array];
        BOOL ok = [g iterBlocksFrom:0 to:30000 threads:threadCounts[ti] error:&err
                         usingBlock:^(TTIOGenomicRun *v, NSUInteger vs, NSUInteger f0,
                                      NSUInteger nr, BOOL *st) {
            (void)st;
            NSMutableIndexSet *local = [NSMutableIndexSet indexSet];
            for (NSUInteger k = 0; k < nr; k++) {
                TTIOAlignedRead *r = [v readAtIndex:vs + k error:NULL];
                if (r) [local addIndex:f0 + k];
            }
            [lock lock];
            [seen addIndexes:local];
            delivered += local.count;
            [ranges addObject:@[@(f0), @(nr)]];
            [lock unlock];
        }];
        /* delivered counts calls, seen counts distinct indices: a read
         * handed over twice would separate them. */
        PASS(ok && delivered == 30000 && seen.count == 30000,
             "bp blocks: threads=%lu delivers every read exactly once (%lu calls, %lu distinct)",
             (unsigned long)threadCounts[ti], (unsigned long)delivered,
             (unsigned long)seen.count);
        [ranges sortUsingComparator:^NSComparisonResult(NSArray *a, NSArray *b) {
            return [a[0] compare:b[0]];
        }];
        for (NSArray *r in ranges) {
            if ([r[0] unsignedIntegerValue] != covered) contiguous = NO;
            covered += [r[1] unsignedIntegerValue];
        }
        PASS(contiguous && covered == 30000,
             "bp blocks: threads=%lu ranges tile the span (%lu covered)",
             (unsigned long)threadCounts[ti], (unsigned long)covered);
    }

    /* A range starting part-way into a block is where viewStart and
     * firstRead part company. Checking the reported indices alone
     * cannot see that; the names have to be compared. */
    __block NSUInteger subCount = 0;
    __block BOOL inRange = YES;
    NSMutableDictionary<NSNumber *, NSString *> *subNames = [NSMutableDictionary dictionary];
    NSLock *sublock = [NSLock new];
    BOOL okSub = [g iterBlocksFrom:12345 to:17890 threads:4 error:&err
                        usingBlock:^(TTIOGenomicRun *v, NSUInteger vs, NSUInteger f0,
                                     NSUInteger nr, BOOL *st) {
        (void)st;
        NSMutableDictionary *local = [NSMutableDictionary dictionary];
        for (NSUInteger k = 0; k < nr; k++) {
            TTIOAlignedRead *r = [v readAtIndex:vs + k error:NULL];
            if (r) local[@(f0 + k)] = r.readName;
        }
        [sublock lock];
        if (f0 < 12345 || f0 + nr > 17890) inRange = NO;
        subCount += nr;
        [subNames addEntriesFromDictionary:local];
        [sublock unlock];
    }];
    PASS(okSub && inRange && subCount == 17890 - 12345,
         "bp blocks: a sub-range delivers only that range (%lu reads)",
         (unsigned long)subCount);
    NSUInteger mismatched = 0;
    for (NSUInteger i = 12345; i < 17890; i++) {
        TTIOAlignedRead *want = [g readAtIndex:i error:NULL];
        if (![subNames[@(i)] isEqualToString:want.readName]) mismatched++;
    }
    PASS(mismatched == 0,
         "bp blocks: a sub-range starting mid-block returns the right records "
         "(%lu wrong)", (unsigned long)mismatched);

    __block NSUInteger emptyCalls = 0;
    [g iterBlocksFrom:200 to:200 threads:4 error:&err
           usingBlock:^(TTIOGenomicRun *v, NSUInteger vs, NSUInteger f0,
                                     NSUInteger nr, BOOL *st) {
        (void)v; (void)vs; (void)f0; (void)nr; (void)st;
        emptyCalls++;
    }];
    PASS(emptyCalls == 0, "bp blocks: an empty range calls nothing");

    /* The window is a memory setting before a concurrency one: a budget
     * that admits one block at a time must still deliver every read. */
    setenv("TTIO_MEMORY_BUDGET", "1048576", 1);
    __block NSUInteger tightCount = 0;
    NSLock *tlock = [NSLock new];
    BOOL okTight = [g iterBlocksFrom:0 to:30000 threads:8 error:&err
                          usingBlock:^(TTIOGenomicRun *v, NSUInteger vs, NSUInteger f0,
                                     NSUInteger nr, BOOL *st) {
        (void)v; (void)vs; (void)f0; (void)st;
        [tlock lock]; tightCount += nr; [tlock unlock];
    }];
    unsetenv("TTIO_MEMORY_BUDGET");
    PASS(okTight && tightCount == 30000,
         "bp blocks: a tight memory budget still delivers every read (%lu)",
         (unsigned long)tightCount);

    [pr close];
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
}

static id<TTIOStorageGroup> gswWriteBudget(NSString *url, TTIOWrittenGenomicRun *run,
                                           NSUInteger threads, unsigned long long budget,
                                           unsigned long long *maxObservedOut)
{
    id<TTIOStorageGroup> study = gswStudy(url);
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = 2000;
    o.threads = threads;
    o.memoryBudgetBytes = budget;
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"g" options:o];
    NSUInteger n = run.readCount;
    for (NSUInteger a = 0; a < n; a += 7001) {
        TTIOWrittenGenomicRun *part = [TTIOGenomicBlocks sliceRun:run from:a to:MIN(n, a + 7001)];
        if (![w appendBatch:part error:&err]) { PASS(NO, "pp: appendBatch (%s)", [[err localizedDescription] UTF8String] ?: ""); return nil; }
    }
    if (![w close:&err]) { PASS(NO, "pp: close (%s)", [[err localizedDescription] UTF8String] ?: ""); return nil; }
    if (maxObservedOut) *maxObservedOut = w.maxInFlightBytesObserved;
    return study;
}

static void gswBudgetBounded(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(40000, 5);
    unsigned long long budget = 4ull << 20;   // 4 MiB, ~1 block in flight
    unsigned long long maxObs = 0;
    id<TTIOStorageGroup> a = gswWriteBudget(
        [NSString stringWithFormat:@"memory://gsw-pp-a-%d", (int)getpid()], run, 6, budget, &maxObs);
    unsigned long long maxSerial = ~0ull;
    id<TTIOStorageGroup> b = gswWriteBudget(
        [NSString stringWithFormat:@"memory://gsw-pp-b-%d", (int)getpid()], run, 1, budget, &maxSerial);
    if (!a || !b) return;
    PASS(maxObs > 0 && maxObs <= budget,
         "pp: in-flight bytes bounded by the budget (%llu <= %llu)", maxObs, budget);
    NSMutableDictionary *ma = [NSMutableDictionary dictionary], *mb = [NSMutableDictionary dictionary];
    gswCollect(a, @"", ma);
    gswCollect(b, @"", mb);
    BOOL same = [ma isEqualToDictionary:mb];
    PASS(same, "pp: budget-stalled file identical to serial (%lu objects)", (unsigned long)ma.count);
}

static void gswStickyMatchesExhaustive(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(40000, 21);
    id<TTIOStorageGroup> a = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-st-a-%d", (int)getpid()], run, 6, 20000);
    setenv("TTIO_M94Z_EXHAUSTIVE", "1", 1);
    id<TTIOStorageGroup> b = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-st-b-%d", (int)getpid()], run, 6, 20000);
    unsetenv("TTIO_M94Z_EXHAUSTIVE");
    if (!a || !b) return;
    NSMutableDictionary *ma = [NSMutableDictionary dictionary], *mb = [NSMutableDictionary dictionary];
    gswCollect(a, @"", ma);
    gswCollect(b, @"", mb);
    PASS([ma isEqualToDictionary:mb],
         "sticky: pinned file identical to exhaustive (%lu objects)",
         (unsigned long)ma.count);
}

static void gswStickyDeterministic(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(40000, 22);
    id<TTIOStorageGroup> a = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-sd-a-%d", (int)getpid()], run, 6, 20000);
    id<TTIOStorageGroup> b = gswWriteThreads([NSString stringWithFormat:@"memory://gsw-sd-b-%d", (int)getpid()], run, 6, 20000);
    if (!a || !b) return;
    NSMutableDictionary *ma = [NSMutableDictionary dictionary], *mb = [NSMutableDictionary dictionary];
    gswCollect(a, @"", ma);
    gswCollect(b, @"", mb);
    PASS([ma isEqualToDictionary:mb],
         "sticky: repeated runs identical (%lu objects)",
         (unsigned long)ma.count);
}

static void gswStickyPinSet(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(40000, 23);
    id<TTIOStorageGroup> study = gswStudy([NSString stringWithFormat:@"memory://gsw-sp-%d", (int)getpid()]);
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = 20000;
    o.threads = 2;
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"g" options:o];
    BOOL ok = [w appendBatch:run error:&err] && [w close:&err];
    PASS(ok && w.qualStrategyHint != -1,
         "sticky: pin set after the first block (hint %ld)",
         (long)w.qualStrategyHint);
}

/* TTIO_M94Z_HINT pins the strategy before block 0 rather than after
 * it. V6 is reachable no other way, which is what the knob is for. */
static void gswHintEnvPins(void)
{
    TTIOWrittenGenomicRun *run = gswBigSyntheticRun(40000, 24);
    id<TTIOStorageGroup> study = gswStudy(
        [NSString stringWithFormat:@"memory://gsw-hint-%d", (int)getpid()]);
    NSError *err = nil;
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
    o.blockReads = 20000;
    o.threads = 2;

    setenv("TTIO_M94Z_HINT", "8", 1);
    TTIOGenomicStreamWriter *w = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"g" options:o];
    NSInteger atInit = w.qualStrategyHint;
    BOOL ok = [w appendBatch:run error:&err] && [w close:&err];
    unsetenv("TTIO_M94Z_HINT");
    PASS(atInit == TTIOM94ZHintV6,
         "hint env: pinned at init, before any block (hint %ld)", (long)atInit);
    PASS(ok && w.qualStrategyHint == TTIOM94ZHintV6,
         "hint env: the pin survives the run (%s)",
         [[err localizedDescription] UTF8String] ?: "");

    /* Junk and zero leave the tune alone, and the exhaustive flag wins
     * over the knob rather than fighting it. */
    setenv("TTIO_M94Z_HINT", "junk", 1);
    TTIOGenomicStreamWriter *wj = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"gj" options:o];
    PASS(wj.qualStrategyHint == -1, "hint env: junk falls through to the tune");
    setenv("TTIO_M94Z_HINT", "0", 1);
    TTIOGenomicStreamWriter *wz = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"gz" options:o];
    PASS(wz.qualStrategyHint == -1, "hint env: zero falls through to the tune");
    setenv("TTIO_M94Z_HINT", "8", 1);
    setenv("TTIO_M94Z_EXHAUSTIVE", "1", 1);
    TTIOGenomicStreamWriter *we = [[TTIOGenomicStreamWriter alloc]
        initWithStudyGroup:study runName:@"ge" options:o];
    PASS(we.qualStrategyHint == -1, "hint env: TTIO_M94Z_EXHAUSTIVE wins over the hint");
    unsetenv("TTIO_M94Z_EXHAUSTIVE");
    unsetenv("TTIO_M94Z_HINT");
}

void testGenomicStreamWriterThreads(void);
void testGenomicStreamWriterThreads(void)
{
    gswRegisterOrder();
    gswThreadedIdentical();
    gswIterReadsThreaded();
    gswBudgetBounded();
    gswStickyMatchesExhaustive();
    gswStickyDeterministic();
    gswStickyPinSet();
    gswHintEnvPins();
    gswIterBlocks();
}
