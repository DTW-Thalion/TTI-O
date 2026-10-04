// TestM103ReadGrouping.m — M103: reads grouped by sequence before
// blocking (layout blocks_v1_grouped).
//
// The writer reorders an unaligned run's reads with the native kernel,
// cuts blocks over the reordered run and stores genomic_index/input_index;
// every read accessor presents input order, block iteration stays in
// stored order, and the grouping survives per-AU encryption, the
// encrypted transport (BlockSidecar input_index, transport-spec v0.13)
// and signatures. Plaintext transport delivers input order and the
// receiver writes it ungrouped.
//
// Mirrors:
//   python/tests/test_m103_read_grouping.py
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#import <Foundation/Foundation.h>
#import "Testing.h"
#import "Codecs/TTIOSeqCm.h"
#import "Codecs/TTIOSeqGroup.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Dataset/TTIOSpectralDataset+GenomicWrite.h"
#import "Genomics/TTIOGenomicStreamWriter.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOGenomicIndex.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Protection/TTIOPerAUFile.h"
#import "Protection/TTIOSignatureManager.h"
#import "Transport/TTIOEncryptedTransport.h"
#import "Transport/TTIOTransportWriter.h"
#import "Transport/TTIOTransportReader.h"
#import "Providers/TTIOHDF5Provider.h"
#import "Providers/TTIOStorageProtocols.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "ValueClasses/TTIOEnums.h"
#include <hdf5.h>
#include <unistd.h>

static NSString *grpTmp(NSString *tag)
{
    return [NSTemporaryDirectory() stringByAppendingPathComponent:
        [NSString stringWithFormat:@"ttio-m103g-%@-%d.tio", tag, (int)getpid()]];
}

static void grpRm(NSString *path)
{
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
}

static NSData *grpKey(void)
{
    NSMutableData *k = [NSMutableData dataWithLength:32];
    uint8_t *p = k.mutableBytes;
    for (int i = 0; i < 32; i++) p[i] = (uint8_t)i;
    return k;
}

static uint32_t grpNext(uint64_t *s)
{
    *s = *s * 6364136223846793005ull + 1442695040888963407ull;
    return (uint32_t)(*s >> 33);
}

static NSData *grpRevComp(const uint8_t *src, NSUInteger L)
{
    NSMutableData *d = [NSMutableData dataWithLength:L];
    uint8_t *o = d.mutableBytes;
    for (NSUInteger i = 0; i < L; i++) {
        uint8_t c = src[L - 1 - i];
        o[i] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    }
    return d;
}

/* Paired reads from both strands of a random genome, shuffled, with a
 * few N bases and a zero-length read: about 12x coverage. Positions are
 * distinct per read (its input index), to check order. */
static TTIOWrittenGenomicRun *grpRun(uint64_t seed, NSUInteger n)
{
    const NSUInteger L = 120, G = 30000;
    uint64_t s = seed;
    NSMutableData *genome = [NSMutableData dataWithLength:G];
    uint8_t *g = genome.mutableBytes;
    for (NSUInteger i = 0; i < G; i++) g[i] = (uint8_t)"ACGT"[grpNext(&s) & 3];
    NSMutableArray<NSData *> *seqs = [NSMutableArray array];
    NSMutableArray<NSString *> *names = [NSMutableArray array];
    for (NSUInteger i = 0; i < n / 2; i++) {
        NSUInteger p = grpNext(&s) % (G - 500);
        NSMutableData *r1 = [NSMutableData dataWithBytes:g + p length:L];
        NSData *r2 = grpRevComp(g + p + 300, L);
        if (i % 37 == 0) ((uint8_t *)r1.mutableBytes)[50] = 'N';
        [seqs addObject:r1];
        [seqs addObject:r2];
        [names addObject:[NSString stringWithFormat:@"frag%lu/1", (unsigned long)i]];
        [names addObject:[NSString stringWithFormat:@"frag%lu/2", (unsigned long)i]];
    }
    seqs[5] = [NSData data];
    NSUInteger m = seqs.count;
    for (NSUInteger i = m - 1; i > 0; i--) {   // Fisher-Yates
        NSUInteger j = grpNext(&s) % (i + 1);
        [seqs exchangeObjectAtIndex:i withObjectAtIndex:j];
        [names exchangeObjectAtIndex:i withObjectAtIndex:j];
    }
    NSMutableData *seq = [NSMutableData data];
    NSMutableData *lengths = [NSMutableData dataWithLength:m * 4];
    NSMutableData *offsets = [NSMutableData dataWithLength:m * 8];
    NSMutableData *positions = [NSMutableData dataWithLength:m * 8];
    NSMutableData *mapqs = [NSMutableData dataWithLength:m];
    NSMutableData *flags = [NSMutableData dataWithLength:m * 4];
    NSMutableData *matePos = [NSMutableData dataWithLength:m * 8];
    NSMutableData *tlens = [NSMutableData dataWithLength:m * 4];
    uint32_t *lp = lengths.mutableBytes;
    uint64_t *op = offsets.mutableBytes;
    int64_t *pp = positions.mutableBytes, *mpp = matePos.mutableBytes;
    uint8_t *mq = mapqs.mutableBytes;
    uint32_t *fp = flags.mutableBytes;
    NSMutableArray *stars = [NSMutableArray arrayWithCapacity:m];
    for (NSUInteger i = 0; i < m; i++) {
        op[i] = seq.length;
        lp[i] = (uint32_t)seqs[i].length;
        [seq appendData:seqs[i]];
        pp[i] = (int64_t)i;
        mq[i] = (uint8_t)(i % 61);
        fp[i] = 0x4;
        mpp[i] = -1;
        [stars addObject:@"*"];
    }
    NSMutableData *qual = [NSMutableData dataWithLength:seq.length];
    uint8_t *qp = qual.mutableBytes;
    for (NSUInteger i = 0; i < qual.length; i++) qp[i] = (uint8_t)(33 + grpNext(&s) % 40);
    return [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:(TTIOAcquisitionMode)7
                   referenceUri:@""
                       platform:@"ILLUMINA"
                     sampleName:@"M103G"
                      positions:positions
               mappingQualities:mapqs
                          flags:flags
                      sequences:seq
                      qualities:qual
                        offsets:offsets
                        lengths:lengths
                         cigars:stars
                      readNames:names
                mateChromosomes:stars
                  matePositions:matePos
                templateLengths:tlens
                    chromosomes:stars
              signalCompression:TTIOCompressionZlib];
}

static NSString *grpSeqAt(TTIOWrittenGenomicRun *run, NSUInteger i)
{
    const uint64_t *op = run.offsetsData.bytes;
    const uint32_t *lp = run.lengthsData.bytes;
    return [[NSString alloc] initWithBytes:(const uint8_t *)run.sequencesData.bytes + op[i]
                                    length:lp[i] encoding:NSASCIIStringEncoding];
}

static BOOL grpWriteImpl(NSString *path, TTIOWrittenGenomicRun *run, BOOL group,
                     NSUInteger blockReads, NSError **error)
{
    grpRm(path);
    if (![TTIOSpectralDataset writeMinimalToPath:path title:@"m103" isaInvestigationId:@"i"
                                          msRuns:@{} genomicRuns:@{} identifications:nil
                                 quantifications:nil provenanceRecords:nil error:error]) return NO;
    TTIOHDF5File *f = [TTIOHDF5File openAtPath:path error:error];
    if (!f) return NO;
    BOOL ok = YES;
    NSError *e = nil;
    @autoreleasepool {   // the writer and its groups go before the file closes
        TTIOHDF5Group *study = [f.rootGroup openGroupNamed:@"study" error:&e];
        TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions optionsFromRun:run];
        o.blockReads = blockReads;
        o.groupReads = group;
        TTIOGenomicStreamWriter *w = study
            ? [[TTIOGenomicStreamWriter alloc] initWithStudyGroup:[TTIOHDF5Provider adapterForGroup:study]
                                                          runName:@"g" options:o]
            : nil;
        ok = w != nil;
        NSUInteger n = run.readCount;
        for (NSUInteger st = 0; ok && st < n; st += 500) {
            ok = [w appendBatch:[TTIOGenomicBlocks sliceRun:run from:st to:MIN(st + 500, n)]
                          error:&e];
        }
        if (ok) ok = [w close:&e];
        [e retain];        // Tests build without ARC
        [w release];
    }
    [f close];
    if (!ok && error) *error = [e autorelease];
    return ok;
}

/* The run group of "g" through the provider adapter; *file to close. */
static id<TTIOStorageGroup> grpRunGroup(NSString *path, TTIOHDF5File **file)
{
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
    *file = [f retain];   // outlives the caller's pool; the caller releases it
    TTIOHDF5Group *rg = [[[f.rootGroup openGroupNamed:@"study" error:NULL]
        openGroupNamed:@"genomic_runs" error:NULL] openGroupNamed:@"g" error:NULL];
    return rg ? [TTIOHDF5Provider adapterForGroup:rg] : nil;
}

static NSString *grpLayoutAttrImpl(NSString *path)
{
    TTIOHDF5File *f = nil;
    NSString *s = nil;
    @autoreleasepool {
        id<TTIOStorageGroup> rg = grpRunGroup(path, &f);
        id v = [rg attributeValueForName:@"layout" error:NULL];
        s = [v isKindOfClass:[NSData class]]
            ? [[NSString alloc] initWithData:v encoding:NSUTF8StringEncoding] : [[v description] copy];
    }
    [f close];
    [f release];
    return [s autorelease];
}

static NSData *grpReadImpl(NSString *path, NSString *sub)
{
    TTIOHDF5File *f = nil;
    NSData *out = nil;
    @autoreleasepool {
        id<TTIOStorageGroup> g = grpRunGroup(path, &f);
        NSArray *parts = [sub componentsSeparatedByString:@"/"];
        for (NSUInteger i = 0; g && i + 1 < parts.count; i++) g = [g openGroupNamed:parts[i] error:NULL];
        id<TTIOStorageDataset> ds = [g openDatasetNamed:parts.lastObject error:NULL];
        id v = ds ? ([sub isEqualToString:@"blocks/index"] ? [ds readCanonicalBytes:NULL]
                                                           : [ds readAll:NULL]) : nil;
        if ([v isKindOfClass:[NSData class]]) out = [[NSData alloc] initWithData:v];
    }
    [f close];
    [f release];
    return [out autorelease];
}

static NSArray<NSString *> *grpNamesImpl(NSString *path)
{
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:NULL];
    TTIOGenomicRun *g = ds.genomicRuns[@"g"] ?: ds.genomicRuns.allValues.firstObject;
    NSMutableArray *out = [NSMutableArray array];
    [g iterReadsFrom:0 to:[g readCount] error:NULL
          usingBlock:^(TTIOAlignedRead *r, NSUInteger i, BOOL *stop) {
        (void)i; (void)stop;
        [out addObject:r.readName ?: @""];
    }];
    [g close];
    return out;
}

/* The helpers above run in their own pool: HDF5 objects (groups, the
 * stream writer) keep a file open until released, and a later open of
 * the same file read-write or a copy of it needs them gone. */
static BOOL grpWrite(NSString *path, TTIOWrittenGenomicRun *run, BOOL group,
                     NSUInteger blockReads, NSError **error)
{
    BOOL ok;
    NSError *e = nil;
    @autoreleasepool { ok = grpWriteImpl(path, run, group, blockReads, &e); [e retain]; }
    if (!ok && error) *error = e;
    [e autorelease];
    return ok;
}

static NSString *grpLayoutAttr(NSString *path)
{
    NSString *s;
    @autoreleasepool { s = [grpLayoutAttrImpl(path) retain]; }
    return [s autorelease];
}

static NSData *grpRead(NSString *path, NSString *sub)
{
    NSData *d;
    @autoreleasepool { d = [grpReadImpl(path, sub) retain]; }
    return [d autorelease];
}

static NSArray<NSString *> *grpNames(NSString *path)
{
    NSArray *a;
    @autoreleasepool { a = [grpNamesImpl(path) retain]; }
    return [a autorelease];
}

/* Rewrite genomic_index/input_index in place through HDF5. */
static void grpEditInputIndex(NSString *path, void (^edit)(uint32_t *ii, NSUInteger n))
{
    TTIOHDF5File *f = [TTIOHDF5File openAtPath:path error:NULL];
    hid_t did = H5Dopen2(f.rootGroup.groupId, "/study/genomic_runs/g/genomic_index/input_index",
                         H5P_DEFAULT);
    if (did >= 0) {
        hid_t sp = H5Dget_space(did);
        hssize_t n = H5Sget_simple_extent_npoints(sp);
        H5Sclose(sp);
        uint32_t *buf = calloc((size_t)(n > 0 ? n : 1), sizeof(uint32_t));
        H5Dread(did, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        edit(buf, (NSUInteger)n);
        H5Dwrite(did, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        free(buf);
        H5Dclose(did);
    }
    [f close];
}

static void grpLayoutAndInputOrder(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 3000);
    NSUInteger n = run.readCount;
    NSString *path = grpTmp(@"grp");
    NSError *err = nil;
    PASS(grpWrite(path, run, YES, 700, &err), "M103 grouping: grouped run written (%s)",
         [[err localizedDescription] UTF8String] ?: "ok");
    PASS([grpLayoutAttr(path) isEqualToString:@"blocks_v1_grouped"], "M103 grouping: @layout blocks_v1_grouped");

    NSData *ii = grpRead(path, @"genomic_index/input_index");
    PASS(ii.length == n * 4, "M103 grouping: input_index has one uint32 per read");
    const uint32_t *iv = ii.bytes;
    NSMutableData *seen = [NSMutableData dataWithLength:n];
    BOOL perm = ii.length == n * 4, identity = YES;
    for (NSUInteger j = 0; perm && j < n; j++) {
        if (iv[j] >= n || ((uint8_t *)seen.mutableBytes)[iv[j]]) { perm = NO; break; }
        ((uint8_t *)seen.mutableBytes)[iv[j]] = 1;
        if (iv[j] != j) identity = NO;
    }
    PASS(perm, "M103 grouping: input_index is a permutation");
    PASS(perm && !identity, "M103 grouping: the reads were reordered");
    NSData *storedPos = grpRead(path, @"genomic_index/positions");
    BOOL moved = storedPos.length == n * 8;
    for (NSUInteger j = 0; moved && j < n; j++) moved = ((const int64_t *)storedPos.bytes)[j] == (int64_t)iv[j];
    PASS(moved, "M103 grouping: every index column moved with its read");
    TTIOHDF5File *f = nil;
    id<TTIOStorageGroup> rg = grpRunGroup(path, &f);
    NSArray *rows = [[[rg openGroupNamed:@"blocks" error:NULL] openDatasetNamed:@"index" error:NULL]
                     readRows:NULL];
    PASS(rows.count > 1, "M103 grouping: more than one block (%lu)", (unsigned long)rows.count);
    [f close];

    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"g"];
    PASS(g != nil && [g.layout isEqualToString:@"blocks_v1_grouped"] && [g readCount] == n,
         "M103 grouping: run opens grouped with every read");
    PASS(g.inputIndex != nil && [g.inputIndex isEqualToData:ii], "M103 grouping: inputIndex exposed");

    __block BOOL orderOk = YES;
    __block NSUInteger got = 0;
    [g iterReadsFrom:0 to:n error:&err usingBlock:^(TTIOAlignedRead *r, NSUInteger i, BOOL *stop) {
        (void)stop;
        if (i != got || ![r.readName isEqualToString:run.readNames[i]]
            || ![r.sequence ?: @"" isEqualToString:grpSeqAt(run, i)] || r.position != (int64_t)i) orderOk = NO;
        got++;
    }];
    PASS(orderOk && got == n, "M103 grouping: iterReads delivers input order");

    __block BOOL thrOk = YES;
    __block NSUInteger thrGot = 0;
    [g iterReadsFrom:0 to:n threads:4 error:&err usingBlock:^(TTIOAlignedRead *r, NSUInteger i, BOOL *stop) {
        (void)stop;
        if (i != thrGot || ![r.readName isEqualToString:run.readNames[i]]) thrOk = NO;
        thrGot++;
    }];
    PASS(thrOk && thrGot == n, "M103 grouping: threaded iterReads delivers input order");

    BOOL atOk = YES;
    NSUInteger picks[6] = {0, 1, 5, n / 2, n - 1, 2999};
    for (int k = 0; k < 6 && atOk; k++) {
        TTIOAlignedRead *r = [g readAtIndex:picks[k] error:&err];
        atOk = r && [r.readName isEqualToString:run.readNames[picks[k]]]
            && [r.sequence ?: @"" isEqualToString:grpSeqAt(run, picks[k])]
            && [[g readNameAtIndex:picks[k] error:NULL] isEqualToString:run.readNames[picks[k]]]
            && [(TTIOAlignedRead *)[g objectAtIndex:picks[k]] position] == (int64_t)picks[k];
    }
    PASS(atOk, "M103 grouping: readAtIndex / readNameAtIndex / objectAtIndex in input order");

    TTIOGenomicIndex *idx = g.index;
    BOOL idxOk = idx.count == n;
    const uint32_t *lens = run.lengthsData.bytes;
    const uint64_t *offs = run.offsetsData.bytes;
    for (NSUInteger i = 0; idxOk && i < n; i++) {
        idxOk = [idx positionAt:i] == (int64_t)i && [idx lengthAt:i] == lens[i]
            && [idx offsetAt:i] == offs[i] && [idx mappingQualityAt:i] == (uint8_t)(i % 61);
    }
    PASS(idxOk, "M103 grouping: index in input order, offsets recomputed");

    NSMutableArray *part = [NSMutableArray array];
    [g iterReadsFrom:1000 to:1100 error:&err usingBlock:^(TTIOAlignedRead *r, NSUInteger i, BOOL *stop) {
        (void)i; (void)stop;
        [part addObject:r.readName];
    }];
    PASS([part isEqualToArray:[run.readNames subarrayWithRange:NSMakeRange(1000, 100)]],
         "M103 grouping: ranged iteration in input order");

    PASS([[g allReadNames] isEqualToArray:run.readNames], "M103 grouping: allReadNames in input order");
    PASS([[g wholeSequencesData] isEqualToData:run.sequencesData]
         && [[g wholeQualitiesData] isEqualToData:run.qualitiesData],
         "M103 grouping: whole sequences / qualities in input order");

    // Block iteration is in stored order; inputIndex maps its rows.
    NSMutableDictionary<NSNumber *, NSString *> *visit = [NSMutableDictionary dictionary];
    NSData *map = g.inputIndex;
    [g iterBlocksFrom:0 to:n threads:1 error:&err
           usingBlock:^(TTIOGenomicRun *view, NSUInteger viewStart, NSUInteger firstRead,
                        NSUInteger nReads, BOOL *stop) {
        (void)stop;
        for (NSUInteger k = 0; k < nReads; k++) {
            NSString *nm = [view readNameAtIndex:viewStart + k error:NULL] ?: @"";
            @synchronized (visit) {
                visit[@(((const uint32_t *)map.bytes)[firstRead + k])] = nm;
            }
        }
    }];
    BOOL blocksOk = visit.count == n;
    for (NSUInteger i = 0; blocksOk && i < n; i++) blocksOk = [visit[@(i)] isEqualToString:run.readNames[i]];
    PASS(blocksOk, "M103 grouping: block iteration stored order, mapped through inputIndex");
    [g close];
    grpRm(path);
}

static void grpShrinks(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 6000);
    NSString *a = grpTmp(@"plain"), *b = grpTmp(@"grp6");
    NSError *err = nil;
    BOOL ok = grpWrite(a, run, NO, 1000, &err) && grpWrite(b, run, YES, 1000, &err);
    NSUInteger sa = grpRead(a, @"signal_channels/sequences/data").length;
    NSUInteger sb = grpRead(b, @"signal_channels/sequences/data").length;
    PASS(ok && sa > 0 && sb < sa * 0.8, "M103 grouping: grouped sequences smaller (%lu vs %lu)",
         (unsigned long)sb, (unsigned long)sa);
    grpRm(a);
    grpRm(b);
}

static void grpDefaultWrite(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 800);
    run.optGroupReads = YES;
    TTIOWrittenGenomicRun *copy = [run copyWithProvenance:@[]];
    PASS(copy.optGroupReads, "M103 grouping: copies keep optGroupReads");
    PASS([TTIOGenomicStreamWriterOptions optionsFromRun:run].groupReads
         && [[[TTIOGenomicStreamWriterOptions optionsFromRun:run] copy] groupReads],
         "M103 grouping: options carry groupReads");
    NSString *path = grpTmp(@"wm");
    grpRm(path);
    NSError *err = nil;
    BOOL ok = [TTIOSpectralDataset writeMinimalToPath:path title:@"m103" isaInvestigationId:@"i"
                                               msRuns:@{} genomicRuns:@{@"g": run}
                                      identifications:nil quantifications:nil provenanceRecords:nil
                                                error:&err];
    PASS(ok && [grpLayoutAttr(path) isEqualToString:@"blocks_v1_grouped"],
         "M103 grouping: default dataset write honours optGroupReads (%s)",
         [[err localizedDescription] UTF8String] ?: "ok");
    PASS([grpNames(path) isEqualToArray:run.readNames], "M103 grouping: default write reads back in input order");
    grpRm(path);
}

static void grpRefusals(void)
{
    TTIOGenomicStreamWriterOptions *o = [TTIOGenomicStreamWriterOptions defaultOptions];
    o.groupReads = YES;
    PASS([TTIOGenomicStreamWriter groupReadsRefusalForOptions:o] == nil, "M103 grouping: plain options may group");
    o.referenceChromSeqs = @{@"chr1": [@"ACGT" dataUsingEncoding:NSASCIIStringEncoding]};
    id<TTIOStorageGroup> noStudy = nil;   // refused before the study is used
    BOOL raised = NO;
    @try {
        (void)[[TTIOGenomicStreamWriter alloc] initWithStudyGroup:noStudy
                                                          runName:@"a" options:o];
    } @catch (NSException *e) {
        raised = [e.name isEqualToString:NSInvalidArgumentException];
    }
    PASS(raised, "M103 grouping: refused with a reference");
    o.referenceChromSeqs = nil;
    o.optLegacyWholeChannel = YES;
    raised = NO;
    @try {
        (void)[[TTIOGenomicStreamWriter alloc] initWithStudyGroup:noStudy
                                                          runName:@"b" options:o];
    } @catch (NSException *e) {
        raised = [e.name isEqualToString:NSInvalidArgumentException];
    }
    PASS(raised, "M103 grouping: refused with the legacy layout");

    // The dataset write turns the refusal into an error.
    TTIOWrittenGenomicRun *run = grpRun(3, 200);
    run.optGroupReads = YES;
    run.referenceChromSeqs = @{@"*": [@"ACGT" dataUsingEncoding:NSASCIIStringEncoding]};
    NSString *path = grpTmp(@"refuse");
    grpRm(path);
    NSError *err = nil;
    BOOL ok = [TTIOSpectralDataset writeMinimalToPath:path title:@"m103" isaInvestigationId:@"i"
                                               msRuns:@{} genomicRuns:@{@"g": run}
                                      identifications:nil quantifications:nil provenanceRecords:nil
                                                error:&err];
    PASS(!ok && err != nil, "M103 grouping: dataset write with a reference fails with an error");
    grpRm(path);
}

static void grpPerAU(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 1500);
    NSString *path = grpTmp(@"perau"), *pristine = grpTmp(@"pristine");
    NSError *err = nil;
    PASS(grpWrite(path, run, YES, 400, &err), "M103 grouping per-AU: written");
    grpRm(pristine);
    [[NSFileManager defaultManager] copyItemAtPath:path toPath:pristine error:NULL];
    PASS([TTIOPerAUFile encryptFilePath:path key:grpKey() encryptHeaders:NO providerName:nil error:&err],
         "M103 grouping per-AU: encrypt (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    PASS([TTIOPerAUFile decryptFilePathInPlace:path key:grpKey() providerName:nil error:&err],
         "M103 grouping per-AU: decrypt in place (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    for (NSString *sub in @[@"blocks/index", @"genomic_index/input_index",
                            @"signal_channels/sequences/data", @"signal_channels/qualities"]) {
        NSData *a = grpRead(path, sub), *b = grpRead(pristine, sub);
        PASS(a != nil && [a isEqualToData:b], "M103 grouping per-AU: %s byte-identical", sub.UTF8String);
    }
    PASS([grpNames(path) isEqualToArray:run.readNames], "M103 grouping per-AU: restored reads in input order");
    grpRm(path);
    grpRm(pristine);
}

static void grpEncryptedTransport(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 1500);
    NSString *path = grpTmp(@"et"), *stream = grpTmp(@"et-stream"), *got = grpTmp(@"et-got");
    NSError *err = nil;
    PASS(grpWrite(path, run, YES, 400, &err), "M103 grouping transport: written");
    NSData *wantII = grpRead(path, @"genomic_index/input_index");
    PASS([TTIOPerAUFile encryptFilePath:path key:grpKey() encryptHeaders:NO providerName:nil error:&err],
         "M103 grouping transport: encrypt");
    grpRm(stream);
    TTIOTransportWriter *tw = [[TTIOTransportWriter alloc] initWithOutputPath:stream];
    BOOL ok = [TTIOEncryptedTransport writeEncryptedDataset:path writer:tw providerName:nil error:&err];
    [tw close];
    PASS(ok, "M103 grouping transport: send (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    NSData *bytes = [NSData dataWithContentsOfFile:stream];
    NSRange tok = [bytes rangeOfData:[@"transport_blocks_v1_grouped" dataUsingEncoding:NSUTF8StringEncoding]
                             options:0 range:NSMakeRange(0, bytes.length)];
    PASS(tok.location != NSNotFound, "M103 grouping transport: StreamHeader carries the grouped token");
    grpRm(got);
    PASS([TTIOEncryptedTransport readEncryptedToPath:got fromStream:bytes providerName:nil error:&err],
         "M103 grouping transport: receive (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    PASS([grpLayoutAttr(got) isEqualToString:@"blocks_v1_grouped"], "M103 grouping transport: layout kept");
    PASS([grpRead(got, @"genomic_index/input_index") isEqualToData:wantII],
         "M103 grouping transport: input_index carried");
    NSString *feats = nil;
    @autoreleasepool {   // the root group must go before decrypt reopens the file
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:got error:NULL];
        feats = [[f.rootGroup stringAttributeNamed:@"ttio_features" error:NULL] copy];
        [f close];
    }
    [feats autorelease];
    PASS(feats != nil && [feats rangeOfString:@"transport_blocks_v1"].location == NSNotFound,
         "M103 grouping transport: wire tokens dropped from the file features");
    PASS([TTIOPerAUFile decryptFilePathInPlace:got key:grpKey() providerName:nil error:&err],
         "M103 grouping transport: decrypt received (%s)", [[err localizedDescription] UTF8String] ?: "ok");
    PASS([grpNames(got) isEqualToArray:run.readNames], "M103 grouping transport: reads in input order");
    grpRm(path);
    grpRm(stream);
    grpRm(got);
}

static void grpPlaintextTransport(void)
{
    TTIOWrittenGenomicRun *run = grpRun(7, 1200);
    NSString *path = grpTmp(@"pt"), *out = grpTmp(@"pt-out");
    NSError *err = nil;
    PASS(grpWrite(path, run, YES, 300, &err), "M103 grouping plaintext: written");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:&err];
    NSMutableData *buf = [NSMutableData data];
    TTIOTransportWriter *tw = [[TTIOTransportWriter alloc] initWithMutableData:buf];
    PASS([tw writeDataset:ds error:&err], "M103 grouping plaintext: send (%s)",
         [[err localizedDescription] UTF8String] ?: "ok");
    [tw close];
    grpRm(out);
    TTIOTransportReader *tr = [[TTIOTransportReader alloc] initWithData:buf];
    PASS([tr writeTtioToPath:out error:&err], "M103 grouping plaintext: receive (%s)",
         [[err localizedDescription] UTF8String] ?: "ok");
    TTIOSpectralDataset *back = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = back.genomicRuns.allValues.firstObject;
    PASS(g != nil && ![g.layout isEqualToString:@"blocks_v1_grouped"] && g.inputIndex == nil,
         "M103 grouping plaintext: receiver writes an ungrouped run");
    NSMutableArray *names = [NSMutableArray array];
    [g iterReadsFrom:0 to:[g readCount] error:NULL usingBlock:^(TTIOAlignedRead *r, NSUInteger i, BOOL *stop) {
        (void)i; (void)stop;
        [names addObject:r.readName ?: @""];
    }];
    PASS([names isEqualToArray:run.readNames], "M103 grouping plaintext: input order delivered");
    [g close];
    grpRm(path);
    grpRm(out);
}

static void grpSignatures(void)
{
    NSString *path = grpTmp(@"sig");
    NSError *err = nil;
    PASS(grpWrite(path, grpRun(7, 600), YES, 200, &err), "M103 grouping signatures: written");
    NSMutableData *key = [NSMutableData dataWithLength:32];
    memset(key.mutableBytes, 'k', 32);
    NSDictionary *sigs = [TTIOSignatureManager signGenomicRun:@"g" inFile:path withKey:key error:&err];
    PASS(sigs[@"genomic_index/input_index"] != nil, "M103 grouping signatures: input_index signed");
    PASS([TTIOSignatureManager verifyGenomicRun:@"g" inFile:path withKey:key error:&err],
         "M103 grouping signatures: verifies");
    grpEditInputIndex(path, ^(uint32_t *ii, NSUInteger n) {
        if (n > 1) { uint32_t t = ii[0]; ii[0] = ii[1]; ii[1] = t; }
    });
    PASS(![TTIOSignatureManager verifyGenomicRun:@"g" inFile:path withKey:key error:NULL],
         "M103 grouping signatures: a swapped input_index fails verification");
    grpRm(path);
}

static void grpCorrupt(void)
{
    NSString *path = grpTmp(@"corrupt");
    NSError *err = nil;
    PASS(grpWrite(path, grpRun(7, 600), YES, 200, &err), "M103 grouping corrupt: written");
    grpEditInputIndex(path, ^(uint32_t *ii, NSUInteger n) { if (n > 1) ii[0] = ii[1]; });
    TTIOHDF5File *f = nil;
    id<TTIOStorageGroup> rg = grpRunGroup(path, &f);
    err = nil;
    TTIOGenomicRun *g = [TTIOGenomicRun openFromGroup:rg name:@"g" error:&err];
    PASS(g == nil && [[err localizedDescription] rangeOfString:@"permutation"].location != NSNotFound,
         "M103 grouping corrupt: a non-permutation input_index is refused (%s)",
         [[err localizedDescription] UTF8String] ?: "no error");
    [f close];
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:NULL];
    PASS(ds.genomicRuns[@"g"] == nil, "M103 grouping corrupt: the dataset does not present the run");
    grpRm(path);
}

void testM103ReadGrouping(void)
{
    if (![TTIOSeqCm nativeAvailable] || ![TTIOSeqGroup nativeAvailable]) {
        PASS(YES, "M103 read grouping: libttio_rans not linked, skipped");
        return;
    }
    @autoreleasepool {
        grpLayoutAndInputOrder();
        grpShrinks();
        grpDefaultWrite();
        grpRefusals();
        grpPerAU();
        grpEncryptedTransport();
        grpPlaintextTransport();
        grpSignatures();
        grpCorrupt();
    }
}
