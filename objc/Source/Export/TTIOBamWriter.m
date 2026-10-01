/*
 * TTIOBamWriter.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOBamWriter
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Export/TTIOBamWriter.h
 *
 * BAM exporter. Formats a TTIOWrittenGenomicRun as SAM text and pipes
 * it via stdin to samtools view -bS - (optionally piped through
 * samtools sort -O bam). Subprocess-only; no htslib linkage.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#import "TTIOBamWriter.h"
#import "Core/TTIOThreads.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Genomics/TTIOGenomicRun.h"
#import "Genomics/TTIOAlignedRead.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "HDF5/TTIOHDF5Errors.h"

// ── Default @SQ length when the writer doesn't know the true reference
// length. SAM requires LN: on every @SQ; we pick INT32_MAX so the
// emitted header is valid for any plausible coordinate. Cross-language
// fixed value (Python's _DEFAULT_SQ_LENGTH).
static const int64_t kTTIODefaultSqLength = 2147483647;

// Mirrors Java BamWriter.PROGRESS_INTERVAL_READS (1000).
const NSUInteger TTIOBamWriterProgressIntervalReads = 1000;

static NSString *const kTTIOWriterSamtoolsHelp =
    @"samtools is required by TTIOBamWriter but was not found on PATH. "
    @"Install it via your platform's package manager:\n"
    @"  Debian/Ubuntu: apt install samtools\n"
    @"  macOS:         brew install samtools\n"
    @"  Conda:         conda install -c bioconda samtools\n"
    @"Then re-run.";

static NSString *bamWriterFindOnPath(NSString *exe)
{
    NSString *path = [[NSProcessInfo processInfo] environment][@"PATH"];
    if (path.length == 0) return nil;
    NSArray *parts = [path componentsSeparatedByString:@":"];
    NSFileManager *fm = [NSFileManager defaultManager];
    for (NSString *dir in parts) {
        if (dir.length == 0) continue;
        NSString *full = [dir stringByAppendingPathComponent:exe];
        if ([fm isExecutableFileAtPath:full]) return full;
    }
    return nil;
}

static BOOL bamWriterSamtoolsAvailable(NSString **outBinary, NSError **error)
{
    NSString *bin = bamWriterFindOnPath(@"samtools");
    if (!bin) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen, @"%@",
                                          kTTIOWriterSamtoolsHelp);
        return NO;
    }
    NSTask *task = [[NSTask alloc] init];
    task.launchPath = bin;
    task.arguments = @[@"--version"];
    NSPipe *out = [NSPipe pipe];
    task.standardOutput = out;
    task.standardError  = out;
    @try {
        [task launch];
    } @catch (NSException *exc) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"%@\n(invocation failed: %@)",
            kTTIOWriterSamtoolsHelp, exc.reason ?: @"unknown");
        return NO;
    }
    [task waitUntilExit];
    [[out fileHandleForReading] readDataToEndOfFile];
    if (task.terminationStatus != 0) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"%@\n(samtools --version exited %d)",
            kTTIOWriterSamtoolsHelp, task.terminationStatus);
        return NO;
    }
    if (outBinary) *outBinary = bin;
    return YES;
}

@implementation TTIOBamWriter

- (instancetype)initWithPath:(NSString *)path
{
    self = [super init];
    if (self) {
        _path = [path copy];
    }
    return self;
}

// ── Header assembly ──────────────────────────────────────────────────
// Mirrors Python BamWriter._build_header in bam.py.
- (NSString *)buildHeaderForRun:(TTIOWrittenGenomicRun *)run
              provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                           sort:(BOOL)sort
{
    return [self buildHeaderForChromosomes:run.chromosomes sampleName:run.sampleName
                                  platform:run.platform provenanceRecords:provenance sort:sort];
}

- (NSString *)buildHeaderForChromosomes:(NSArray<NSString *> *)chromosomes
                             sampleName:(NSString *)sampleName
                               platform:(NSString *)platform
                      provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                                   sort:(BOOL)sort
{
    NSMutableArray<NSString *> *lines = [NSMutableArray array];
    NSString *so = sort ? @"coordinate" : @"unsorted";
    [lines addObject:[NSString stringWithFormat:@"@HD\tVN:1.6\tSO:%@", so]];

    // @SQ — first-seen order, dropping "*" and empty strings.
    NSMutableSet<NSString *> *seen = [NSMutableSet set];
    for (NSString *chrom in chromosomes) {
        if (chrom.length == 0) continue;
        if ([chrom isEqualToString:@"*"]) continue;
        if ([seen containsObject:chrom]) continue;
        [seen addObject:chrom];
        [lines addObject:[NSString stringWithFormat:@"@SQ\tSN:%@\tLN:%lld",
                                                    chrom, (long long)kTTIODefaultSqLength]];
    }

    // @RG — single line if either sample_name or platform is set.
    if (sampleName.length > 0 || platform.length > 0) {
        NSMutableArray<NSString *> *parts = [NSMutableArray array];
        [parts addObject:@"@RG"];
        [parts addObject:@"ID:rg1"];
        if (sampleName.length > 0) {
            [parts addObject:[NSString stringWithFormat:@"SM:%@", sampleName]];
        }
        if (platform.length > 0) {
            [parts addObject:[NSString stringWithFormat:@"PL:%@", platform]];
        }
        [lines addObject:[parts componentsJoinedByString:@"\t"]];
    }

    // @PG — one line per provenance record, with .1/.2 collision suffix
    // matching Python's behaviour.
    NSMutableSet<NSString *> *usedIds = [NSMutableSet set];
    NSUInteger idx = 0;
    for (TTIOProvenanceRecord *prov in provenance) {
        NSString *baseId = (prov.software.length > 0)
            ? prov.software
            : [NSString stringWithFormat:@"pg%lu", (unsigned long)idx];
        NSString *pgId = baseId;
        NSUInteger n = 1;
        while ([usedIds containsObject:pgId]) {
            pgId = [NSString stringWithFormat:@"%@.%lu", baseId, (unsigned long)n];
            n++;
        }
        [usedIds addObject:pgId];

        NSMutableArray<NSString *> *parts = [NSMutableArray array];
        [parts addObject:@"@PG"];
        [parts addObject:[NSString stringWithFormat:@"ID:%@", pgId]];
        [parts addObject:[NSString stringWithFormat:@"PN:%@", prov.software ?: @""]];
        id clVal = prov.parameters ? prov.parameters[@"CL"] : nil;
        if ([clVal isKindOfClass:[NSString class]] && [(NSString *)clVal length] > 0) {
            [parts addObject:[NSString stringWithFormat:@"CL:%@", (NSString *)clVal]];
        }
        [lines addObject:[parts componentsJoinedByString:@"\t"]];
        idx++;
    }

    return [[lines componentsJoinedByString:@"\n"] stringByAppendingString:@"\n"];
}

// ── Per-read alignment line assembly ─────────────────────────────────
// Mirrors Python BamWriter._iter_alignment_lines in bam.py.
- (NSString *)buildAlignmentLinesForRun:(TTIOWrittenGenomicRun *)run
{
    return [self buildAlignmentLinesForRun:run progress:nil];
}

- (NSString *)buildAlignmentLinesForRun:(TTIOWrittenGenomicRun *)run
                                progress:(TTIOProgressBlock)progress
{
    if (progress == nil) progress = TTIOProgressDiscard();
    NSUInteger n = run.readNames.count;
    NSMutableString *out = [NSMutableString string];

    const int64_t  *pos      = (const int64_t  *)run.positionsData.bytes;
    const uint8_t  *mapq     = (const uint8_t  *)run.mappingQualitiesData.bytes;
    const uint32_t *flags    = (const uint32_t *)run.flagsData.bytes;
    const uint64_t *offsets  = (const uint64_t *)run.offsetsData.bytes;
    const uint32_t *lengths  = (const uint32_t *)run.lengthsData.bytes;
    const int64_t  *matePos  = (const int64_t  *)run.matePositionsData.bytes;
    const int32_t  *tlens    = (const int32_t  *)run.templateLengthsData.bytes;
    const uint8_t  *seqBuf   = (const uint8_t  *)run.sequencesData.bytes;
    const uint8_t  *qualBuf  = (const uint8_t  *)run.qualitiesData.bytes;
    NSUInteger seqBufLen     = run.sequencesData.length;
    NSUInteger qualBufLen    = run.qualitiesData.length;

    for (NSUInteger i = 0; i < n; i++) {
        NSString *qname = run.readNames[i];
        if (qname.length == 0) qname = @"*";

        uint32_t flag = flags ? flags[i] : 0;
        NSString *rname = (i < run.chromosomes.count) ? run.chromosomes[i] : @"*";
        if (rname.length == 0) rname = @"*";
        int64_t  p     = pos  ? pos[i]  : 0;
        uint8_t  mq    = mapq ? mapq[i] : 0;
        NSString *cigar = (i < run.cigars.count) ? run.cigars[i] : @"*";
        if (cigar.length == 0) cigar = @"*";

        NSString *mateChrom = (i < run.mateChromosomes.count)
            ? run.mateChromosomes[i] : @"*";
        if (mateChrom.length == 0) mateChrom = @"*";

        // RNEXT collapse §136
        NSString *rnext;
        if ([mateChrom isEqualToString:rname] && ![rname isEqualToString:@"*"]) {
            rnext = @"=";
        } else {
            rnext = mateChrom;
        }

        // PNEXT mapping §138: mate_position < 0 → 0
        int64_t mp = matePos ? matePos[i] : 0;
        int64_t pnext = (mp < 0) ? 0 : mp;

        int32_t tlen = tlens ? tlens[i] : 0;

        uint64_t off = offsets ? offsets[i] : 0;
        uint32_t len = lengths ? lengths[i] : 0;

        NSString *seqStr;
        NSString *qualStr;
        if (len == 0) {
            seqStr  = @"*";
            qualStr = @"*";
        } else {
            // SEQ from sequences buffer (ASCII bases)
            if (off + len > seqBufLen) {
                // Defensive — should never happen on a sane WGR.
                seqStr = @"*";
            } else {
                seqStr = [[NSString alloc]
                          initWithBytes:seqBuf + off
                                 length:len
                               encoding:NSASCIIStringEncoding] ?: @"*";
            }

            // QUAL from qualities buffer. M87 fills with 0xFF when source
            // SAM had QUAL '*' but a non-empty SEQ; map back to '*' on
            // write so the round trip canonicalises to the source.
            if (qualBuf && off + len <= qualBufLen) {
                BOOL allFF = YES;
                for (uint32_t k = 0; k < len; k++) {
                    if (qualBuf[off + k] != 0xFF) { allFF = NO; break; }
                }
                if (allFF) {
                    qualStr = @"*";
                } else {
                    // Bytes already stored as ASCII Phred+33 — emit verbatim
                    // (Latin-1 keeps the bytes intact for any value > 0x7F
                    // that might sneak through).
                    qualStr = [[NSString alloc]
                               initWithBytes:qualBuf + off
                                      length:len
                                    encoding:NSISOLatin1StringEncoding] ?: @"*";
                }
            } else {
                qualStr = @"*";
            }
        }

        // M101: the optional fields follow column 11 verbatim.
        NSString *tags = (i < run.tags.count) ? run.tags[i] : @"";
        [out appendFormat:@"%@\t%u\t%@\t%lld\t%u\t%@\t%@\t%lld\t%d\t%@\t%@%@%@\n",
                          qname, (unsigned)flag, rname, (long long)p,
                          (unsigned)mq, cigar, rnext, (long long)pnext,
                          (int)tlen, seqStr, qualStr,
                          tags.length ? @"\t" : @"", tags];

        // Per-N progress fire. Total = n.
        if (((i + 1) % TTIOBamWriterProgressIntervalReads) == 0) {
            progress((int64_t)(i + 1), (int64_t)n);
        }
    }
    // Final fire.
    progress((int64_t)n, (int64_t)n);
    return out;
}

// One SAM alignment line from a decoded read; the same field rules as
// the run-side builder above.
- (NSString *)alignmentLineForRead:(TTIOAlignedRead *)r
{
    NSString *qname = r.readName.length ? r.readName : @"*";
    NSString *rname = r.chromosome.length ? r.chromosome : @"*";
    NSString *cigar = r.cigar.length ? r.cigar : @"*";
    NSString *mateChrom = r.mateChromosome.length ? r.mateChromosome : @"*";
    NSString *rnext = ([mateChrom isEqualToString:rname] && ![rname isEqualToString:@"*"]) ? @"=" : mateChrom;
    int64_t pnext = r.matePosition < 0 ? 0 : r.matePosition;
    NSString *seqStr, *qualStr;
    NSUInteger len = r.sequence.length;
    if (len == 0) {
        seqStr = @"*"; qualStr = @"*";
    } else {
        seqStr = r.sequence;
        NSData *q = r.qualities;
        BOOL allFF = q.length > 0;
        const uint8_t *qb = q.bytes;
        for (NSUInteger k = 0; k < q.length; k++) if (qb[k] != 0xFF) { allFF = NO; break; }
        qualStr = (q.length == 0 || allFF) ? @"*"
            : ([[NSString alloc] initWithData:q encoding:NSISOLatin1StringEncoding] ?: @"*");
    }
    // M101: the optional fields follow column 11 verbatim.
    NSString *tags = r.tags ?: @"";
    return [NSString stringWithFormat:@"%@\t%u\t%@\t%lld\t%u\t%@\t%@\t%lld\t%d\t%@\t%@%@%@\n",
            qname, (unsigned)r.flags, rname, (long long)r.position, (unsigned)r.mappingQuality,
            cigar, rnext, (long long)pnext, (int)r.templateLength, seqStr, qualStr,
            tags.length ? @"\t" : @"", tags];
}

// ── samtools command builder (overridable). Returns an array of arg
// arrays. Single element for one-stage; two for view + sort.
- (NSArray<NSArray<NSString *> *> *)samtoolsCommandsForSort:(BOOL)sort
{
    if (sort) {
        return @[
            @[@"view", @"-bS", @"-"],
            @[@"sort", @"-O", @"bam", @"-o", self.path, @"-"],
        ];
    } else {
        return @[
            @[@"view", @"-bS", @"-o", self.path, @"-"],
        ];
    }
}

// ── samtools subprocess invocation. Single- or two-stage pipeline. ──
- (BOOL)runSamtoolsPipeline:(NSArray<NSArray<NSString *> *> *)commands
                samtoolsBin:(NSString *)samtoolsBin
                   samBytes:(NSData *)samBytes
                      error:(NSError **)error
{
    return [self runSamtoolsPipeline:commands samtoolsBin:samtoolsBin
                              feeder:^BOOL(NSFileHandle *stdinHandle, NSError **e) {
        (void)e;
        [stdinHandle writeData:samBytes];
        return YES;
    } error:error];
}

// The feeder writes the SAM text into samtools' stdin as it goes and
// returns NO (with *error) to abort; the handle is closed here.
- (BOOL)runSamtoolsPipeline:(NSArray<NSArray<NSString *> *> *)commands
                samtoolsBin:(NSString *)samtoolsBin
                     feeder:(BOOL (^)(NSFileHandle *stdinHandle, NSError **error))feeder
                      error:(NSError **)error
{
    if (commands.count == 1) {
        NSTask *task = [[NSTask alloc] init];
        task.launchPath = samtoolsBin;
        task.arguments = commands[0];
        NSPipe *inPipe  = [NSPipe pipe];
        NSPipe *outPipe = [NSPipe pipe];
        NSPipe *errPipe = [NSPipe pipe];
        task.standardInput  = inPipe;
        task.standardOutput = outPipe;
        task.standardError  = errPipe;
        @try {
            [task launch];
        } @catch (NSException *exc) {
            if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
                @"failed to launch samtools: %@", exc.reason ?: @"unknown");
            return NO;
        }
        BOOL fed = YES;
        @try {
            fed = feeder([inPipe fileHandleForWriting], error);
            [[inPipe fileHandleForWriting] closeFile];
        } @catch (NSException *exc) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetWrite,
                @"failed to pipe SAM text to samtools: %@",
                exc.reason ?: @"unknown");
            fed = NO;
        }
        [[outPipe fileHandleForReading] readDataToEndOfFile];
        NSData *errData = [[errPipe fileHandleForReading] readDataToEndOfFile];
        [task waitUntilExit];
        if (!fed) return NO;
        if (task.terminationStatus != 0) {
            NSString *errText = [[NSString alloc] initWithData:errData
                                                       encoding:NSUTF8StringEncoding] ?: @"";
            errText = [errText stringByTrimmingCharactersInSet:
                [NSCharacterSet whitespaceAndNewlineCharacterSet]];
            if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
                @"samtools exited %d: %@", task.terminationStatus,
                errText.length ? errText : @"(no stderr)");
            return NO;
        }
        return YES;
    }

    // Two-stage: stage[0].stdout -> stage[1].stdin.
    NSTask *first  = [[NSTask alloc] init];
    NSTask *second = [[NSTask alloc] init];
    first.launchPath  = samtoolsBin;
    second.launchPath = samtoolsBin;
    first.arguments  = commands[0];
    second.arguments = commands[1];

    NSPipe *firstIn  = [NSPipe pipe];
    NSPipe *between  = [NSPipe pipe];
    NSPipe *firstErr = [NSPipe pipe];
    NSPipe *secondOut = [NSPipe pipe];
    NSPipe *secondErr = [NSPipe pipe];

    first.standardInput   = firstIn;
    first.standardOutput  = between;
    first.standardError   = firstErr;
    second.standardInput  = between;
    second.standardOutput = secondOut;
    second.standardError  = secondErr;

    @try {
        [first launch];
        [second launch];
    } @catch (NSException *exc) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"failed to launch samtools: %@", exc.reason ?: @"unknown");
        return NO;
    }

    BOOL fed = YES;
    @try {
        fed = feeder([firstIn fileHandleForWriting], error);
        [[firstIn fileHandleForWriting] closeFile];
    } @catch (NSException *exc) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetWrite,
            @"failed to pipe SAM text to samtools: %@",
            exc.reason ?: @"unknown");
        fed = NO;
    }

    [first waitUntilExit];
    // Give samtools sort a chance to finish reading the pipe.
    [[between fileHandleForWriting] closeFile];
    NSData *firstErrData  = [[firstErr fileHandleForReading] readDataToEndOfFile];
    [[secondOut fileHandleForReading] readDataToEndOfFile];
    NSData *secondErrData = [[secondErr fileHandleForReading] readDataToEndOfFile];
    [second waitUntilExit];

    if (!fed) return NO;
    if (first.terminationStatus != 0) {
        NSString *errText = [[NSString alloc] initWithData:firstErrData
                                                   encoding:NSUTF8StringEncoding] ?: @"";
        errText = [errText stringByTrimmingCharactersInSet:
            [NSCharacterSet whitespaceAndNewlineCharacterSet]];
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"samtools (stage 1) exited %d: %@",
            first.terminationStatus,
            errText.length ? errText : @"(no stderr)");
        return NO;
    }
    if (second.terminationStatus != 0) {
        NSString *errText = [[NSString alloc] initWithData:secondErrData
                                                   encoding:NSUTF8StringEncoding] ?: @"";
        errText = [errText stringByTrimmingCharactersInSet:
            [NSCharacterSet whitespaceAndNewlineCharacterSet]];
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"samtools (stage 2) exited %d: %@",
            second.terminationStatus,
            errText.length ? errText : @"(no stderr)");
        return NO;
    }
    return YES;
}

- (BOOL)writeRun:(TTIOWrittenGenomicRun *)run
   provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                sort:(BOOL)sort
               error:(NSError **)error
{
    return [self writeRun:run
        provenanceRecords:provenance
                     sort:sort
                 progress:nil
                    error:error];
}

- (BOOL)writeReadSideRun:(TTIOGenomicRun *)run
       provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                    sort:(BOOL)sort
                progress:(TTIOProgressBlock)progress
                   error:(NSError **)error
{
    NSString *samtoolsBin = nil;
    if (!bamWriterSamtoolsAvailable(&samtoolsBin, error)) {
        return NO;
    }
    if (!run) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"writeReadSideRun: run must not be nil");
        return NO;
    }
    TTIOProgressBlock sink = progress ?: TTIOProgressDiscard();
    NSArray<TTIOProvenanceRecord *> *provs = provenance ?: @[];
    NSString *header = [self buildHeaderForChromosomes:[run chromosomeNames]
                                            sampleName:run.sampleName
                                              platform:run.platform
                                     provenanceRecords:provs
                                                  sort:sort];
    NSUInteger total = run.readCount;
    NSArray<NSArray<NSString *> *> *commands = [self samtoolsCommandsForSort:sort];
    return [self runSamtoolsPipeline:commands samtoolsBin:samtoolsBin
                              feeder:^BOOL(NSFileHandle *stdinHandle, NSError **e) {
        [stdinHandle writeData:[header dataUsingEncoding:NSASCIIStringEncoding] ?: [NSData data]];
        NSMutableData *chunk = [NSMutableData dataWithCapacity:1 << 20];
        __block NSUInteger done = 0;
        BOOL ok = [run iterReadsFrom:0 to:total threads:[TTIOThreads resolve:nil] error:e
                          usingBlock:^(TTIOAlignedRead *r, NSUInteger index, BOOL *stop) {
            (void)index; (void)stop;
            NSString *line = [self alignmentLineForRead:r];
            [chunk appendData:[line dataUsingEncoding:NSISOLatin1StringEncoding] ?: [NSData data]];
            done++;
            if (chunk.length >= (1 << 20)) {
                [stdinHandle writeData:chunk];
                [chunk setLength:0];
            }
            if ((done % TTIOBamWriterProgressIntervalReads) == 0 && done < total) {
                sink((int64_t)done, (int64_t)total);
            }
        }];
        if (chunk.length > 0) [stdinHandle writeData:chunk];
        if (ok) sink((int64_t)total, (int64_t)total);
        return ok;
    } error:error];
}

- (BOOL)writeReadSideRun:(TTIOGenomicRun *)run
       provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                    sort:(BOOL)sort
                   error:(NSError **)error
{
    return [self writeReadSideRun:run provenanceRecords:provenance sort:sort progress:nil error:error];
}

- (BOOL)writeRun:(TTIOWrittenGenomicRun *)run
   provenanceRecords:(NSArray<TTIOProvenanceRecord *> *)provenance
                sort:(BOOL)sort
            progress:(TTIOProgressBlock)progress
               error:(NSError **)error
{
    NSString *samtoolsBin = nil;
    if (!bamWriterSamtoolsAvailable(&samtoolsBin, error)) {
        return NO;
    }
    if (!run) {
        if (error) *error = TTIOMakeError(TTIOErrorInvalidArgument,
            @"writeRun: run must not be nil");
        return NO;
    }

    NSArray<TTIOProvenanceRecord *> *provs = provenance ?: @[];

    NSString *header = [self buildHeaderForRun:run
                             provenanceRecords:provs
                                          sort:sort];
    NSString *aligns = [self buildAlignmentLinesForRun:run progress:progress];
    NSString *sam = [header stringByAppendingString:aligns];
    NSData *samBytes = [sam dataUsingEncoding:NSASCIIStringEncoding] ?: [NSData data];

    NSArray<NSArray<NSString *> *> *commands = [self samtoolsCommandsForSort:sort];

    return [self runSamtoolsPipeline:commands
                          samtoolsBin:samtoolsBin
                             samBytes:samBytes
                                error:error];
}

@end
