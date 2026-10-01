// TestM101SamTags.m — M101: SAM optional tags through import, storage,
// read and export.
//
// The SAM fixture is built in the test: random reference, reads with
// mismatches, deletions, insertions, soft clips and an unmapped read,
// and MD/NM written by `samtools calmd` as the independent ground
// truth.
//
// Mirrors:
//   python/tests/test_m101_sam_tags.py
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#import <Foundation/Foundation.h>
#import <hdf5.h>
#import "Testing.h"
#import "Codecs/TTIOSamTags.h"
#import "Dataset/TTIOSpectralDataset.h"
#import "Export/TTIOBamWriter.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "HDF5/TTIOFeatureFlags.h"
#import "HDF5/TTIOHDF5File.h"
#import "HDF5/TTIOHDF5Group.h"
#import "Import/TTIOBamReader.h"
#import "Import/TTIOImporterRegistry.h"
#import "Protection/TTIOPerAUFile.h"
#import "Protection/TTIOSignatureManager.h"
#import "Providers/TTIOHDF5Provider.h"
#import "Providers/TTIOStorageProtocols.h"
#import "Transport/TTIOAccessUnit.h"
#import "Transport/TTIOEncryptedTransport.h"
#import "Transport/TTIOTransportPacket.h"
#import "Transport/TTIOTransportReader.h"
#import "Transport/TTIOTransportWriter.h"
#import "ValueClasses/TTIOEnums.h"
#include <stdlib.h>
#include <unistd.h>

static const NSUInteger kM101RefLen = 20000;

static NSString *m101Dir(NSString *tag)
{
    NSString *d = [NSTemporaryDirectory() stringByAppendingPathComponent:
        [NSString stringWithFormat:@"ttio-m101-%@-%d", tag, (int)getpid()]];
    [[NSFileManager defaultManager] removeItemAtPath:d error:NULL];
    [[NSFileManager defaultManager] createDirectoryAtPath:d
                              withIntermediateDirectories:YES attributes:nil error:NULL];
    return d;
}

static void m101Rm(NSString *path)
{
    [[NSFileManager defaultManager] removeItemAtPath:path error:NULL];
}

static NSString *m101Samtools(void)
{
    NSString *path = [[NSProcessInfo processInfo] environment][@"PATH"];
    for (NSString *dir in [path componentsSeparatedByString:@":"]) {
        if (dir.length == 0) continue;
        NSString *full = [dir stringByAppendingPathComponent:@"samtools"];
        if ([[NSFileManager defaultManager] isExecutableFileAtPath:full]) return full;
    }
    return nil;
}

/* Run samtools with `args`; stdout goes to `outPath` (or is returned
 * when outPath is nil). nil on a non-zero exit. */
static NSData *m101RunSamtools(NSArray<NSString *> *args, NSString *outPath)
{
    NSTask *t = [[NSTask alloc] init];
    t.launchPath = m101Samtools();
    t.arguments = args;
    NSPipe *pipe = nil;
    NSFileHandle *fh = nil;
    if (outPath) {
        [[NSFileManager defaultManager] createFileAtPath:outPath contents:nil attributes:nil];
        fh = [NSFileHandle fileHandleForWritingAtPath:outPath];
        t.standardOutput = fh;
    } else {
        pipe = [NSPipe pipe];
        t.standardOutput = pipe;
    }
    t.standardError = [NSFileHandle fileHandleWithNullDevice];
    @try { [t launch]; } @catch (NSException *e) { return nil; }
    NSData *out = pipe ? [[pipe fileHandleForReading] readDataToEndOfFile] : [NSData data];
    [t waitUntilExit];
    [fh closeFile];
    return t.terminationStatus == 0 ? out : nil;
}

/* Deterministic generator for the fixture (the values need not match
 * the Python fixture's; calmd is the ground truth either way). */
static uint64_t gM101Rng = 11;
static uint32_t m101Rand(uint32_t n)
{
    gM101Rng = gM101Rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (uint32_t)((gM101Rng >> 33) % n);
}

/* Write ref.fa and in.sam (coordinate-sorted, one chromosome) into
 * `dir`. The first `taglessPrefix` reads carry no optional fields. */
static BOOL m101MakeSam(NSString *dir, NSUInteger n, BOOL tags,
                        NSUInteger taglessPrefix,
                        NSString **outSam, NSString **outRef)
{
    gM101Rng = 11;
    const char *acgt = "ACGT";
    NSMutableString *ref = [NSMutableString stringWithCapacity:kM101RefLen];
    for (NSUInteger i = 0; i < kM101RefLen; i++) {
        [ref appendFormat:@"%c", acgt[m101Rand(4)]];
    }
    NSString *refPath = [dir stringByAppendingPathComponent:@"ref.fa"];
    NSMutableString *fa = [NSMutableString stringWithString:@">chr1\n"];
    for (NSUInteger i = 0; i < kM101RefLen; i += 60) {
        [fa appendString:[ref substringWithRange:NSMakeRange(i, MIN((NSUInteger)60, kM101RefLen - i))]];
        [fa appendString:@"\n"];
    }
    [fa writeToFile:refPath atomically:YES encoding:NSASCIIStringEncoding error:NULL];
    if (!m101RunSamtools(@[@"faidx", refPath], nil)) return NO;

    NSMutableArray<NSString *> *lines = [NSMutableArray array];
    [lines addObject:@"@HD\tVN:1.6\tSO:coordinate"];
    [lines addObject:[NSString stringWithFormat:@"@SQ\tSN:chr1\tLN:%lu", (unsigned long)kM101RefLen]];
    [lines addObject:@"@RG\tID:g1\tSM:S1\tPL:ILLUMINA"];
    const char *qa = "#+5:?FI";
    for (NSUInteger k = 0; k <= n; k++) {
        NSMutableArray<NSString *> *rec = [NSMutableArray array];
        if (k == n) {
            [rec addObjectsFromArray:@[@"u0", @"4", @"*", @"0", @"0", @"*", @"*", @"0", @"0",
                                       @"ACGTACGTAC", @"IIIIIIIIII"]];
        } else {
            NSUInteger i = k;
            NSUInteger pos = 1 + (i * (kM101RefLen - 400)) / n;
            NSUInteger kind = i % 6;
            NSMutableString *q = [NSMutableString string];
            for (int j = 0; j < 100; j++) [q appendFormat:@"%c", qa[m101Rand(7)]];
            NSUInteger p0 = pos - 1;
            NSMutableString *seq = [[ref substringWithRange:NSMakeRange(p0, 100)] mutableCopy];
            NSString *cigar = @"100M";
            if (kind == 1) {                                   // mismatches
                NSMutableSet *used = [NSMutableSet set];
                while (used.count < 3) {
                    NSUInteger at = m101Rand(100);
                    if ([used containsObject:@(at)]) continue;
                    [used addObject:@(at)];
                    unichar c = [seq characterAtIndex:at];
                    NSString *sub = c == 'A' ? @"C" : c == 'C' ? @"G" : c == 'G' ? @"T" : @"A";
                    [seq replaceCharactersInRange:NSMakeRange(at, 1) withString:sub];
                }
            } else if (kind == 2) {                            // deletion
                seq = [[[ref substringWithRange:NSMakeRange(p0, 40)]
                        stringByAppendingString:[ref substringWithRange:NSMakeRange(p0 + 45, 60)]] mutableCopy];
                cigar = @"40M5D60M";
            } else if (kind == 3) {                            // insertion
                seq = [[NSString stringWithFormat:@"%@TTT%@",
                        [ref substringWithRange:NSMakeRange(p0, 50)],
                        [ref substringWithRange:NSMakeRange(p0 + 50, 47)]] mutableCopy];
                cigar = @"50M3I47M";
            } else if (kind == 4) {                            // soft clip
                seq = [[@"GGGGG" stringByAppendingString:[ref substringWithRange:NSMakeRange(p0, 95)]] mutableCopy];
                cigar = @"5S95M";
            }
            BOOL first = (i % 2 == 0);
            [rec addObjectsFromArray:@[
                [NSString stringWithFormat:@"r%05lu", (unsigned long)i],
                first ? @"99" : @"147", @"chr1",
                [NSString stringWithFormat:@"%lu", (unsigned long)pos], @"60", cigar, @"=",
                [NSString stringWithFormat:@"%lu", (unsigned long)(pos + 150)],
                first ? @"250" : @"-250", seq, q]];
        }
        if (tags && k >= taglessPrefix) {
            uint32_t a = m101Rand(300);
            [rec addObject:[NSString stringWithFormat:@"AS:i:%u", a]];
            [rec addObject:[NSString stringWithFormat:@"UQ:i:%u", a]];
            [rec addObject:@"RG:Z:g1"];
            [rec addObject:[NSString stringWithFormat:@"XS:i:%d", (int)m101Rand(10) - 5]];
            if (k % 7 == 0) {
                [rec addObject:[NSString stringWithFormat:@"XA:Z:chr1,+%u,100M,%u;",
                                1 + m101Rand(9998), m101Rand(4)]];
            }
            if (k % 11 == 0) [rec addObject:@"BC:B:c,-3,0,7"];
            if (k % 13 == 0) [rec addObject:@"XF:f:0.25"];
        }
        [lines addObject:[rec componentsJoinedByString:@"\t"]];
    }
    NSString *raw = [dir stringByAppendingPathComponent:@"raw.sam"];
    [[[lines componentsJoinedByString:@"\n"] stringByAppendingString:@"\n"]
        writeToFile:raw atomically:YES encoding:NSASCIIStringEncoding error:NULL];
    NSString *sam = [dir stringByAppendingPathComponent:@"in.sam"];
    if (tags) {
        // calmd appends NM/MD where they are missing (all tagged reads).
        if (!m101RunSamtools(@[@"calmd", raw, refPath], sam)) return NO;
        if (taglessPrefix > 0) {
            // Strip calmd's tags from the tag-less prefix again.
            NSString *text = [NSString stringWithContentsOfFile:sam encoding:NSASCIIStringEncoding error:NULL];
            NSMutableArray *out = [NSMutableArray array];
            for (NSString *line in [text componentsSeparatedByString:@"\n"]) {
                if (line.length == 0) continue;
                NSArray *c = [line componentsSeparatedByString:@"\t"];
                if (![line hasPrefix:@"@"] && [c[0] hasPrefix:@"r"]
                    && (NSUInteger)[[c[0] substringFromIndex:1] integerValue] < taglessPrefix) {
                    c = [c subarrayWithRange:NSMakeRange(0, 11)];
                }
                [out addObject:[c componentsJoinedByString:@"\t"]];
            }
            [[[out componentsJoinedByString:@"\n"] stringByAppendingString:@"\n"]
                writeToFile:sam atomically:YES encoding:NSASCIIStringEncoding error:NULL];
        }
    } else {
        [[NSFileManager defaultManager] copyItemAtPath:raw toPath:sam error:NULL];
    }
    if (outSam) *outSam = sam;
    if (outRef) *outRef = refPath;
    return YES;
}

/* Sorted SAM records of `samtools view <path>`, RNEXT "=" expanded to
 * RNAME; `withTags` NO keeps columns 1-11 only. */
static NSArray<NSString *> *m101SamLines(NSString *path, BOOL withTags)
{
    NSData *d = m101RunSamtools(@[@"view", path], nil);
    if (!d) return nil;
    NSString *text = [[NSString alloc] initWithData:d encoding:NSUTF8StringEncoding];
    NSMutableArray *out = [NSMutableArray array];
    for (NSString *line in [text componentsSeparatedByString:@"\n"]) {
        if (line.length == 0) continue;
        NSMutableArray *c = [[line componentsSeparatedByString:@"\t"] mutableCopy];
        if ([c[6] isEqualToString:@"="]) c[6] = c[2];
        if (!withTags && c.count > 11) [c removeObjectsInRange:NSMakeRange(11, c.count - 11)];
        [out addObject:[c componentsJoinedByString:@"\t"]];
    }
    return [out sortedArrayUsingSelector:@selector(compare:)];
}

/* The same records computed from a stored run's reads. */
static NSArray<NSString *> *m101RunLines(TTIOGenomicRun *run, BOOL withTags)
{
    NSMutableArray *out = [NSMutableArray array];
    NSError *err = nil;
    for (NSUInteger i = 0; i < run.readCount; i++) {
        TTIOAlignedRead *r = [run readAtIndex:i error:&err];
        if (!r) return nil;
        NSString *seq = r.sequence.length ? r.sequence : @"*";
        NSString *qual = @"*";
        const uint8_t *qb = r.qualities.bytes;
        BOOL allFF = YES;
        for (NSUInteger k = 0; k < r.qualities.length; k++) if (qb[k] != 0xFF) { allFF = NO; break; }
        if (r.qualities.length > 0 && !allFF) {
            qual = [[NSString alloc] initWithData:r.qualities encoding:NSISOLatin1StringEncoding];
        }
        NSMutableString *line = [NSMutableString stringWithFormat:
            @"%@\t%u\t%@\t%lld\t%u\t%@\t%@\t%lld\t%d\t%@\t%@",
            r.readName.length ? r.readName : @"*", (unsigned)r.flags,
            r.chromosome.length ? r.chromosome : @"*", (long long)r.position,
            (unsigned)r.mappingQuality, r.cigar.length ? r.cigar : @"*",
            r.mateChromosome.length ? r.mateChromosome : @"*",
            (long long)r.matePosition, (int)r.templateLength, seq, qual];
        if (withTags && r.tags.length) [line appendFormat:@"\t%@", r.tags];
        [out addObject:line];
    }
    return [out sortedArrayUsingSelector:@selector(compare:)];
}

static BOOL m101Import(NSString *sam, NSString *ref, NSString *out,
                       NSString *blockReads, NSError **err)
{
    NSMutableDictionary *opts = [NSMutableDictionary dictionary];
    if (ref) {
        opts[@"reference"] = ref;
        opts[@"embed_reference"] = @"true";
    }
    if (blockReads) opts[@"block_reads"] = blockReads;
    m101Rm(out);
    return [TTIOImporterRegistry encodeFormat:@"sam" inputs:@[sam] output:out
                                      options:opts error:err];
}

static id<TTIOStorageGroup> m101RunGroup(TTIOHDF5File *f, NSString *name)
{
    TTIOHDF5Group *rg = [[[f.rootGroup openGroupNamed:@"study" error:NULL]
        openGroupNamed:@"genomic_runs" error:NULL]
        openGroupNamed:name error:NULL];
    return rg ? [TTIOHDF5Provider adapterForGroup:rg] : nil;
}

// The readers below drain their own pool so no HDF5 handle outlives
// the call (a later read-write open would fail). The test driver is
// MRC, so results are retained across the drain.
static BOOL m101HasFeature(NSString *path, NSString *flag)
{
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
        if (!f) return NO;
        BOOL has = [TTIOFeatureFlags root:f.rootGroup supportsFeature:flag];
        [f close];
        return has;
    }
}

/* blocks/index rows, signal_channels/tags bytes (nil when absent). */
static NSArray *m101IndexRows(NSString *path, NSString *run)
{
    NSArray *rows = nil;
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
        if (!f) return nil;
        rows = [[[m101RunGroup(f, run) openGroupNamed:@"blocks" error:NULL]
            openDatasetNamed:@"index" error:NULL] readRows:NULL];
        [rows retain];
        [f close];
    }
    return [rows autorelease];
}

static NSData *m101TagsBytes(NSString *path, NSString *run, uint8_t *codec)
{
    NSData *d = nil;
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
        if (!f) return nil;
        id<TTIOStorageGroup> sig = [m101RunGroup(f, run) openGroupNamed:@"signal_channels" error:NULL];
        if ([sig hasChildNamed:@"tags"]) {
            id<TTIOStorageDataset> ds = [sig openDatasetNamed:@"tags" error:NULL];
            d = [ds readAll:NULL];
            if (codec) *codec = (uint8_t)[[ds attributeValueForName:@"compression" error:NULL] unsignedIntegerValue];
        }
        [d retain];
        [f close];
    }
    return [d autorelease];
}

/* index canonical bytes, tags, sequences/data, qualities. */
static NSArray<NSData *> *m101ChannelState(NSString *path)
{
    NSData *idx = nil, *tags = nil, *seq = nil, *qual = nil;
    @autoreleasepool {
        TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
        if (!f) return nil;
        id<TTIOStorageGroup> a = m101RunGroup(f, @"genomic_0001");
        idx = [[[a openGroupNamed:@"blocks" error:NULL]
            openDatasetNamed:@"index" error:NULL] readCanonicalBytes:NULL];
        id<TTIOStorageGroup> sig = [a openGroupNamed:@"signal_channels" error:NULL];
        tags = [[sig openDatasetNamed:@"tags" error:NULL] readAll:NULL];
        seq = [[[sig openGroupNamed:@"sequences" error:NULL]
            openDatasetNamed:@"data" error:NULL] readAll:NULL];
        qual = [[sig openDatasetNamed:@"qualities" error:NULL] readAll:NULL];
        [idx retain]; [tags retain]; [seq retain]; [qual retain];
        [f close];
    }
    [idx autorelease]; [tags autorelease]; [seq autorelease]; [qual autorelease];
    if (!idx || !tags || !seq || !qual) return nil;
    return @[idx, tags, seq, qual];
}

static NSData *m101Key(void)
{
    NSMutableData *k = [NSMutableData dataWithLength:32];
    memset(k.mutableBytes, 0x65, 32);
    return k;
}

// ── codec ────────────────────────────────────────────────────────────

static void m101CodecRoundTrip(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 120, YES, 0, &sam, &ref), "M101 codec: fixture written");
    NSError *err = nil;
    TTIOWrittenGenomicRun *run = [[[TTIOBamReader alloc] initWithPath:sam]
        toGenomicRunWithName:@"g" region:nil sampleName:nil error:&err];
    PASS(run != nil && run.tags.count == 121 && [run.tags[0] rangeOfString:@"MD:Z:"].location != NSNotFound,
         "M101 codec: importer keeps columns 12+ (%s)", [[err localizedDescription] UTF8String] ?: "");

    // Reference bases of chr1, from the FASTA the fixture wrote.
    NSString *fa = [NSString stringWithContentsOfFile:ref encoding:NSASCIIStringEncoding error:NULL];
    NSMutableString *bases = [NSMutableString string];
    for (NSString *l in [fa componentsSeparatedByString:@"\n"]) {
        if (l.length && ![l hasPrefix:@">"]) [bases appendString:l];
    }
    NSData *refBytes = [bases dataUsingEncoding:NSASCIIStringEncoding];

    NSUInteger n = run.readCount;
    NSMutableData *off = [NSMutableData dataWithLength:(n + 1) * 8];
    NSMutableData *cid = [NSMutableData dataWithLength:n * 2];
    uint64_t *o = off.mutableBytes;
    uint16_t *c = cid.mutableBytes;
    const uint32_t *lens = run.lengthsData.bytes;
    for (NSUInteger i = 0; i < n; i++) {
        o[i + 1] = o[i] + lens[i];
        c[i] = [run.chromosomes[i] isEqualToString:@"chr1"] ? 0 : 0xFFFF;
    }
    TTIOSamTagsContext *ctx = [[TTIOSamTagsContext alloc] init];
    ctx.sequences = run.sequencesData;
    ctx.seqOffsets = off;
    ctx.cigars = run.cigars;
    ctx.positions = run.positionsData;
    ctx.chromIds = cid;
    ctx.references = @[refBytes];

    NSData *derived = [TTIOSamTags encodeTags:run.tags context:ctx error:&err];
    NSData *plain = [TTIOSamTags encodeTags:run.tags context:nil error:&err];
    PASS(derived.length > 0 && plain.length > 0 && derived.length < plain.length,
         "M101 codec: MD/NM derivation shrinks the blob (%lu < %lu)",
         (unsigned long)derived.length, (unsigned long)plain.length);
    PASS(derived.length >= 4 && memcmp(derived.bytes, "STG1", 4) == 0, "M101 codec: STG1 magic");
    NSArray *back = [TTIOSamTags decodeData:derived nReads:n context:ctx error:&err];
    PASS([back isEqualToArray:run.tags], "M101 codec: derived round trip");
    back = [TTIOSamTags decodeData:plain nReads:n context:nil error:&err];
    PASS([back isEqualToArray:run.tags], "M101 codec: verbatim round trip");

    // Non-canonical text is kept verbatim; empty input is fine.
    NSArray *odd = @[@"", @"NM:i:00\tMD:Z:5", @"XX:Z:café", @"not a tag", @"AS:i:1\t"];
    NSData *ob = [TTIOSamTags encodeTags:odd context:nil error:&err];
    PASS([[TTIOSamTags decodeData:ob nReads:odd.count context:nil error:&err] isEqualToArray:odd],
         "M101 codec: verbatim storage of non-canonical text");
    NSData *eb = [TTIOSamTags encodeTags:@[] context:nil error:&err];
    PASS(eb != nil && [[TTIOSamTags decodeData:eb nReads:0 context:nil error:&err] count] == 0,
         "M101 codec: empty input");
    err = nil;
    NSString *withNul = [NSString stringWithFormat:@"AS:i:1%C", (unichar)0];
    PASS([TTIOSamTags encodeTags:@[withNul] context:nil error:&err] == nil && err != nil,
         "M101 codec: NUL in the tag text is refused");
}

// ── import / read / export ───────────────────────────────────────────

static void m101ImportWithReference(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 0, &sam, &ref), "M101 import: fixture written");
    NSString *out = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, out, @"200", &err), "M101 import: encode with reference (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    NSArray *want = m101SamLines(sam, YES);
    PASS(g != nil && want.count == 601 && [m101RunLines(g, YES) isEqualToArray:want],
         "M101 import: every record (tags included) reads back");
    [g close];
    [ds closeFile];
    PASS(m101HasFeature(out, @"opt_sam_tags"), "M101 import: opt_sam_tags set");
    NSArray *rows = m101IndexRows(out, @"genomic_0001");
    BOOL codecOk = rows.count == 4 && rows[0][@"tags_off"] != nil;
    for (NSDictionary *r in rows) {
        if ([r[@"tags_len"] unsignedLongLongValue] > 0
            && [r[@"tags_codec"] unsignedIntValue] != TTIOCompressionSamTags) codecOk = NO;
    }
    PASS(codecOk, "M101 import: index carries the tags triple, codec 18");
    NSData *stored = m101TagsBytes(out, @"genomic_0001", NULL);
    NSUInteger text = 0;
    for (NSString *line in want) {
        NSArray *c = [line componentsSeparatedByString:@"\t"];
        if (c.count > 11) text += [[[c subarrayWithRange:NSMakeRange(11, c.count - 11)]
                                     componentsJoinedByString:@"\t"] length];
    }
    PASS(stored.length > 0 && stored.length * 8 < text,
         "M101 import: tags stored well under the text (%lu vs %lu)",
         (unsigned long)stored.length, (unsigned long)text);

    // Export: the optional fields follow column 11.
    NSString *bam = [dir stringByAppendingPathComponent:@"o.bam"];
    ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    g = ds.genomicRuns[@"genomic_0001"];
    TTIOBamWriter *w = [[TTIOBamWriter alloc] initWithPath:bam];
    PASS([w writeReadSideRun:g provenanceRecords:@[] sort:NO error:&err],
         "M101 export: BAM written (%s)", [[err localizedDescription] UTF8String] ?: "");
    [g close];
    [ds closeFile];
    PASS([m101SamLines(bam, YES) isEqualToArray:want], "M101 export: records equal the source SAM");
}

static void m101ImportWithoutReference(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 0, &sam, &ref), "M101 no-ref: fixture written");
    NSString *out = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, nil, out, @"250", &err), "M101 no-ref: encode (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    PASS(g != nil && [m101RunLines(g, YES) isEqualToArray:m101SamLines(sam, YES)],
         "M101 no-ref: every record reads back");
    [g close];
}

static void m101Tagless(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, NO, 0, &sam, &ref), "M101 tagless: fixture written");
    NSString *out = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, out, @"200", &err), "M101 tagless: encode (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    PASS(m101TagsBytes(out, @"genomic_0001", NULL) == nil, "M101 tagless: no tags dataset");
    PASS(!m101HasFeature(out, @"opt_sam_tags"), "M101 tagless: no opt_sam_tags");
    NSArray *rows = m101IndexRows(out, @"genomic_0001");
    PASS(rows.count > 0 && rows[0][@"tags_off"] == nil && rows[0][@"mate_info_codec"] != nil,
         "M101 tagless: index without tags columns");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    PASS(g != nil && ![g hasTagsChannel]
         && [m101RunLines(g, NO) isEqualToArray:m101SamLines(sam, NO)],
         "M101 tagless: columns 1-11 read back");
    BOOL empty = YES;
    for (NSUInteger i = 0; i < g.readCount && empty; i++) {
        empty = [[g readAtIndex:i error:NULL].tags isEqualToString:@""];
    }
    PASS(empty, "M101 tagless: every read's tags are empty");
    [g close];

    // Plaintext transport: no tags channel in the AUs.
    NSMutableData *buf = [NSMutableData data];
    ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOTransportWriter *tw = [[TTIOTransportWriter alloc] initWithMutableData:buf];
    PASS([tw writeDataset:ds error:&err], "M101 tagless: transport stream written");
    [tw close];
    BOOL noTags = YES;
    NSUInteger nAU = 0;
    for (TTIOTransportPacketRecord *rec in [[[TTIOTransportReader alloc] initWithData:buf]
                                            readAllPacketsWithError:&err]) {
        if (rec.header.packetType != TTIOTransportPacketAccessUnit) continue;
        TTIOAccessUnit *au = [TTIOAccessUnit decodeFromBytes:rec.payload.bytes
                                                       length:rec.payload.length error:&err];
        nAU++;
        for (TTIOTransportChannelData *c in au.channels) {
            if ([c.name isEqualToString:@"tags"]) noTags = NO;
        }
    }
    PASS(nAU == 601 && noTags, "M101 tagless: plaintext stream has no tags channel");
}

static void m101LateTags(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 250, &sam, &ref), "M101 late tags: fixture written");
    NSString *out = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, out, @"100", &err), "M101 late tags: encode (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    NSArray *rows = m101IndexRows(out, @"genomic_0001");
    PASS(rows.count == 7 && rows[0][@"tags_len"] != nil
         && [rows[0][@"tags_len"] unsignedLongLongValue] == 0
         && [rows[1][@"tags_len"] unsignedLongLongValue] == 0
         && [rows.lastObject[@"tags_len"] unsignedLongLongValue] > 0,
         "M101 late tags: index upgraded once, early blocks empty");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    PASS(g != nil && [m101RunLines(g, YES) isEqualToArray:m101SamLines(sam, YES)],
         "M101 late tags: every record reads back");
    [g close];
}

// ── protection ───────────────────────────────────────────────────────

static void m101PerAU(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 0, &sam, &ref), "M101 per-AU: fixture written");
    NSString *path = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, path, @"200", &err), "M101 per-AU: encode (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    NSArray<NSData *> *before = m101ChannelState(path);
    PASS(before != nil, "M101 per-AU: pre-encrypt snapshot");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    NSMutableArray *want = [NSMutableArray array];
    for (NSUInteger i = 0; i < g.readCount; i++) [want addObject:[g readAtIndex:i error:NULL].tags ?: @""];
    [g close];
    [ds closeFile];
    ds = nil;

    PASS([TTIOPerAUFile encryptFilePath:path key:m101Key() encryptHeaders:NO
                           providerName:nil error:&err],
         "M101 per-AU: encrypt (%s)", [[err localizedDescription] UTF8String] ?: "");
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
    id<TTIOStorageGroup> sig = [m101RunGroup(f, @"genomic_0001") openGroupNamed:@"signal_channels" error:NULL];
    PASS(![sig hasChildNamed:@"tags"] && [sig hasChildNamed:@"tags_segments"],
         "M101 per-AU: tags replaced by tags_segments");
    [f close];
    NSDictionary *plain = [TTIOPerAUFile decryptFilePath:path key:m101Key() providerName:nil error:&err];
    PASS([plain[@"genomic_0001"][@"tags"] isEqualToArray:want],
         "M101 per-AU: read-only decrypt returns the tags");

    PASS([TTIOPerAUFile decryptFilePathInPlace:path key:m101Key() providerName:nil error:&err],
         "M101 per-AU: decrypt in place (%s)", [[err localizedDescription] UTF8String] ?: "");
    NSArray<NSData *> *after = m101ChannelState(path);
    PASS(after != nil && [after isEqualToArray:before],
         "M101 per-AU: index, tags, sequences, qualities byte-identical");
    ds = [TTIOSpectralDataset readFromFilePath:path error:&err];
    g = ds.genomicRuns[@"genomic_0001"];
    PASS(g != nil && [m101RunLines(g, YES) isEqualToArray:m101SamLines(sam, YES)],
         "M101 per-AU: restored records equal the source");
    [g close];
    [ds closeFile];
    ds = nil;

    // Encrypted transport: the receiver resolves the REF_DIFF reference
    // (and the MD/NM derivation) through REF_PATH.
    PASS(m101Import(sam, ref, path, @"200", &err), "M101 enc transport: encode");
    setenv("REF_PATH", ref.fileSystemRepresentation, 1);
    PASS([TTIOPerAUFile encryptFilePath:path key:m101Key() encryptHeaders:NO
                           providerName:nil error:&err], "M101 enc transport: encrypt");
    NSString *stream = [dir stringByAppendingPathComponent:@"e.tis"];
    m101Rm(stream);
    TTIOTransportWriter *tw = [[TTIOTransportWriter alloc] initWithOutputPath:stream];
    BOOL ok = [TTIOEncryptedTransport writeEncryptedDataset:path writer:tw
                                                providerName:nil error:&err];
    [tw close];
    PASS(ok, "M101 enc transport: send (%s)", [[err localizedDescription] UTF8String] ?: "");
    NSString *received = [dir stringByAppendingPathComponent:@"r.tio"];
    m101Rm(received);
    PASS([TTIOEncryptedTransport readEncryptedToPath:received
                                          fromStream:[NSData dataWithContentsOfFile:stream]
                                        providerName:nil error:&err],
         "M101 enc transport: receive (%s)", [[err localizedDescription] UTF8String] ?: "");
    PASS([TTIOPerAUFile decryptFilePathInPlace:received key:m101Key() providerName:nil error:&err],
         "M101 enc transport: decrypt received (%s)", [[err localizedDescription] UTF8String] ?: "");
    NSArray<NSData *> *recv = m101ChannelState(received);
    PASS(recv != nil && [recv[0] isEqualToData:before[0]] && [recv[1] isEqualToData:before[1]],
         "M101 enc transport: index and tags byte-identical");
    unsetenv("REF_PATH");
}

/* Number of HDF5 filters on a dataset; -1 when it cannot be opened. */
static int m101FilterCount(NSString *path, const char *dsPath)
{
    hid_t fid = H5Fopen(path.fileSystemRepresentation, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) return -1;
    int n = -1;
    hid_t did = H5Dopen2(fid, dsPath, H5P_DEFAULT);
    if (did >= 0) {
        hid_t dcpl = H5Dget_create_plist(did);
        n = H5Pget_nfilters(dcpl);
        H5Pclose(dcpl);
        H5Dclose(did);
    }
    H5Fclose(fid);
    return n;
}

/* The restore creates a channel dataset at its first non-empty blob:
 * tags first appearing after untagged blocks come back with the stream
 * writer's codec 18 and no filter, byte-identical. */
static void m101PerAULateTags(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 250, &sam, &ref), "M101 per-AU late tags: fixture written");
    NSString *path = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, path, @"100", &err), "M101 per-AU late tags: encode (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    NSArray *rows = m101IndexRows(path, @"genomic_0001");
    PASS(rows.count > 2 && [rows[0][@"tags_len"] unsignedLongLongValue] == 0,
         "M101 per-AU late tags: block 0 has no tags");
    NSArray<NSData *> *before = m101ChannelState(path);
    const char *tagsPath = "/study/genomic_runs/genomic_0001/signal_channels/tags";
    PASS(before != nil && m101FilterCount(path, tagsPath) == 0,
         "M101 per-AU late tags: written tags channel has no filter");
    PASS([TTIOPerAUFile encryptFilePath:path key:m101Key() encryptHeaders:NO
                           providerName:nil error:&err],
         "M101 per-AU late tags: encrypt (%s)", [[err localizedDescription] UTF8String] ?: "");
    PASS([TTIOPerAUFile decryptFilePathInPlace:path key:m101Key() providerName:nil error:&err],
         "M101 per-AU late tags: decrypt in place (%s)", [[err localizedDescription] UTF8String] ?: "");
    NSArray<NSData *> *after = m101ChannelState(path);
    PASS(after != nil && [after[0] isEqualToData:before[0]],
         "M101 per-AU late tags: blocks/index byte-identical");
    PASS(after != nil && [after[1] isEqualToData:before[1]]
         && [after[2] isEqualToData:before[2]] && [after[3] isEqualToData:before[3]],
         "M101 per-AU late tags: tags, sequences, qualities byte-identical");
    uint8_t codec = 0;
    PASS(m101TagsBytes(path, @"genomic_0001", &codec) != nil && codec == TTIOCompressionSamTags,
         "M101 per-AU late tags: restored tags @compression = 18");
    PASS(m101FilterCount(path, tagsPath) == 0,
         "M101 per-AU late tags: restored tags dataset has no HDF5 filter");
}

static void m101WholeChannel(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 120, YES, 0, &sam, &ref), "M101 whole: fixture written");
    NSError *err = nil;
    TTIOWrittenGenomicRun *run = [[[[TTIOBamReader alloc] initWithPath:sam]
        toGenomicRunWithName:@"g" region:nil sampleName:nil error:&err]
        copyWithOptLegacyWholeChannel:YES];
    NSString *path = [dir stringByAppendingPathComponent:@"o.tio"];
    m101Rm(path);
    PASS([TTIOSpectralDataset writeMinimalToPath:path title:@"t" isaInvestigationId:@"i"
                                          msRuns:@{} genomicRuns:@{@"g": run}
                                 identifications:nil quantifications:nil
                               provenanceRecords:nil error:&err],
         "M101 whole: written (%s)", [[err localizedDescription] UTF8String] ?: "");
    uint8_t codec = 0;
    PASS(m101TagsBytes(path, @"g", &codec) != nil && codec == TTIOCompressionSamTags,
         "M101 whole: tags dataset @compression = 18");
    PASS(m101HasFeature(path, @"opt_sam_tags"), "M101 whole: opt_sam_tags set");
    @autoreleasepool {
        // Scoped so every HDF5 handle of the read is gone before the
        // read-write open below.
        TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:path error:&err];
        TTIOGenomicRun *g = ds.genomicRuns[@"g"];
        PASS(g != nil && ![g.layout isEqualToString:@"blocks_v1"]
             && [m101RunLines(g, YES) isEqualToArray:m101SamLines(sam, YES)],
             "M101 whole: every record reads back");
        [g close];
        [ds closeFile];
    }

    err = nil;
    PASS(![TTIOPerAUFile encryptFilePath:path key:m101Key() encryptHeaders:NO
                            providerName:nil error:&err]
         && [[err localizedDescription] rangeOfString:@"SAM tags"].location != NSNotFound,
         "M101 whole: per-AU encryption refuses a whole-channel tagged run");
    TTIOHDF5File *f = [TTIOHDF5File openReadOnlyAtPath:path error:NULL];
    id<TTIOStorageGroup> sig = [m101RunGroup(f, @"g") openGroupNamed:@"signal_channels" error:NULL];
    PASS([sig hasChildNamed:@"tags"] && [sig hasChildNamed:@"sequences"],
         "M101 whole: refused run left in plaintext");
    [f close];
}

static void m101Signatures(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 120, YES, 0, &sam, &ref), "M101 signatures: fixture written");
    NSString *path = [dir stringByAppendingPathComponent:@"o.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, path, nil, &err), "M101 signatures: encode");
    NSMutableData *key = [NSMutableData dataWithLength:32];
    memset(key.mutableBytes, 'k', 32);
    NSDictionary *sigs = [TTIOSignatureManager signGenomicRun:@"genomic_0001" inFile:path
                                                      withKey:key error:&err];
    PASS(sigs[@"signal_channels/tags"] != nil, "M101 signatures: tags dataset signed");
    PASS([TTIOSignatureManager verifyGenomicRun:@"genomic_0001" inFile:path withKey:key error:&err],
         "M101 signatures: verify");

    hid_t fid = H5Fopen(path.fileSystemRepresentation, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t did = H5Dopen2(fid, "/study/genomic_runs/genomic_0001/signal_channels/tags", H5P_DEFAULT);
    hid_t sp = H5Dget_space(did);
    hssize_t n = H5Sget_simple_extent_npoints(sp);
    NSMutableData *buf = [NSMutableData dataWithLength:(NSUInteger)n];
    H5Dread(did, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.mutableBytes);
    ((uint8_t *)buf.mutableBytes)[n - 1] ^= 1;
    H5Dwrite(did, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.bytes);
    H5Sclose(sp); H5Dclose(did); H5Fclose(fid);
    err = nil;
    PASS(![TTIOSignatureManager verifyGenomicRun:@"genomic_0001" inFile:path withKey:key error:&err],
         "M101 signatures: a flipped tags byte fails verification");
}

static void m101PlaintextTransport(NSString *dir)
{
    NSString *sam = nil, *ref = nil;
    PASS(m101MakeSam(dir, 600, YES, 0, &sam, &ref), "M101 transport: fixture written");
    NSString *src = [dir stringByAppendingPathComponent:@"src.tio"];
    NSError *err = nil;
    PASS(m101Import(sam, ref, src, @"200", &err), "M101 transport: encode");
    TTIOSpectralDataset *ds = [TTIOSpectralDataset readFromFilePath:src error:&err];
    NSMutableData *buf = [NSMutableData data];
    TTIOTransportWriter *tw = [[TTIOTransportWriter alloc] initWithMutableData:buf];
    PASS([tw writeDataset:ds error:&err], "M101 transport: stream written (%s)",
         [[err localizedDescription] UTF8String] ?: "");
    [tw close];
    [ds closeFile];
    ds = nil;
    NSString *out = [dir stringByAppendingPathComponent:@"out.tio"];
    m101Rm(out);
    PASS([[[TTIOTransportReader alloc] initWithData:buf] writeTtioToPath:out error:&err],
         "M101 transport: materialised (%s)", [[err localizedDescription] UTF8String] ?: "");
    ds = [TTIOSpectralDataset readFromFilePath:out error:&err];
    TTIOGenomicRun *g = ds.genomicRuns[@"genomic_0001"];
    PASS(g != nil && [m101RunLines(g, YES) isEqualToArray:m101SamLines(sam, YES)],
         "M101 transport: every record (tags included) survives");
    [g close];
    PASS(m101HasFeature(out, @"opt_sam_tags"), "M101 transport: opt_sam_tags on the received file");
}

void testM101SamTags(void);
void testM101SamTags(void)
{
    if (m101Samtools() == nil || ![TTIOSamTags nativeAvailable]) {
        PASS(YES, "M101: samtools or libttio_rans unavailable, skipped");
        return;
    }
    NSString *dir = m101Dir(@"codec");
    m101CodecRoundTrip(dir);
    m101Rm(dir);
    dir = m101Dir(@"import");
    m101ImportWithReference(dir);
    m101Rm(dir);
    dir = m101Dir(@"noref");
    m101ImportWithoutReference(dir);
    m101Rm(dir);
    dir = m101Dir(@"tagless");
    m101Tagless(dir);
    m101Rm(dir);
    dir = m101Dir(@"late");
    m101LateTags(dir);
    m101Rm(dir);
    dir = m101Dir(@"perau");
    m101PerAU(dir);
    m101Rm(dir);
    dir = m101Dir(@"perau-late");
    m101PerAULateTags(dir);
    m101Rm(dir);
    dir = m101Dir(@"whole");
    m101WholeChannel(dir);
    m101Rm(dir);
    dir = m101Dir(@"sign");
    m101Signatures(dir);
    m101Rm(dir);
    dir = m101Dir(@"transport");
    m101PlaintextTransport(dir);
    m101Rm(dir);
}
