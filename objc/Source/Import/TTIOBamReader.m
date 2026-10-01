/*
 * TTIOBamReader.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOBamReader
 * Inherits From: NSObject
 * Conforms To:   NSObject (NSObject)
 * Declared In:   Import/TTIOBamReader.h
 *
 * SAM/BAM importer. Wraps the user-installed samtools binary as an
 * NSTask subprocess to materialise SAM/BAM files into
 * TTIOWrittenGenomicRun instances; no htslib source is linked.
 * Probes samtools availability lazily on the first read call.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#import "TTIOBamReader.h"
#import "Import/TTIOOrderedBatchAssembler.h"
#import "Core/TTIOThreads.h"
#import "Import/TTIOGenomicStreamSource.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#import "Dataset/TTIOProvenanceRecord.h"
#import "ValueClasses/TTIOEnums.h"
#import "HDF5/TTIOHDF5Errors.h"
#import <sys/stat.h>

// Mirrors Java BamReader.PROGRESS_INTERVAL_READS (1000).
const NSUInteger TTIOBamReaderProgressIntervalReads = 1000;

// ── Install help text per Binding Decision §135.
// Must contain the substrings "apt", "brew", "conda" so cross-language
// tests can grep for at least one.
static NSString *const kTTIOSamtoolsInstallHelp =
    @"samtools is required by TTIOBamReader but was not found on PATH. "
    @"Install it via your platform's package manager:\n"
    @"  Debian/Ubuntu: apt install samtools\n"
    @"  macOS:         brew install samtools\n"
    @"  Conda:         conda install -c bioconda samtools\n"
    @"Then re-run.";

// ── PATH lookup ───────────────────────────────────────────────────────

static NSString *bamFindOnPath(NSString *exe)
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

// Returns YES iff `samtools --version` exits 0.
static BOOL bamSamtoolsAvailable(NSString **outBinary, NSError **error)
{
    NSString *bin = bamFindOnPath(@"samtools");
    if (!bin) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen, @"%@",
                                          kTTIOSamtoolsInstallHelp);
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
            kTTIOSamtoolsInstallHelp, exc.reason ?: @"unknown");
        return NO;
    }
    [task waitUntilExit];
    [[out fileHandleForReading] readDataToEndOfFile];
    if (task.terminationStatus != 0) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"%@\n(samtools --version exited %d)",
            kTTIOSamtoolsInstallHelp, task.terminationStatus);
        return NO;
    }
    if (outBinary) *outBinary = bin;
    return YES;
}

// Split a UTF8 line into tab fields, with `maxFields` analogous to
// Python's `str.split("\t", maxsplit)`. The last element absorbs any
// remaining tabs verbatim.
static NSArray<NSString *> *bamSplitTabsLimited(NSString *line, NSUInteger maxFields)
{
    NSMutableArray<NSString *> *fields = [NSMutableArray array];
    NSUInteger len = line.length;
    NSUInteger start = 0;
    for (NSUInteger i = 0; i < len; i++) {
        unichar c = [line characterAtIndex:i];
        if (c != '\t') continue;
        if (fields.count + 1 >= maxFields) {
            // We've collected (maxFields - 1) fields already; the
            // remainder of the line (including any further tabs) is
            // the final field. Mirrors Python's split(..., maxsplit).
            break;
        }
        [fields addObject:[line substringWithRange:NSMakeRange(start, i - start)]];
        start = i + 1;
    }
    [fields addObject:[line substringFromIndex:start]];
    return fields;
}

// Parse a SAM header line (after the @TAG token) into a {KEY: VALUE}
// dict. The header lines have the form `@TAG\tKEY:VALUE\tKEY:VALUE...`
static NSDictionary<NSString *, NSString *> *bamParseHeaderFields(NSString *line)
{
    NSMutableDictionary *fields = [NSMutableDictionary dictionary];
    NSArray *tokens = [line componentsSeparatedByString:@"\t"];
    for (NSUInteger i = 1; i < tokens.count; i++) {
        NSString *tok = tokens[i];
        NSRange colon = [tok rangeOfString:@":"];
        if (colon.location == NSNotFound) continue;
        NSString *key = [tok substringToIndex:colon.location];
        NSString *val = [tok substringFromIndex:colon.location + 1];
        fields[key] = val;
    }
    return fields;
}

// ─────────────────────────────────────────────────────────────────────

const NSUInteger TTIOBamReaderDefaultBatchReads = 100000;

/* Consumes samtools SAM text one line at a time: header lines set the
 * run-level metadata and provenance, alignment lines accumulate until
 * -drainBatchWithSampleName: hands them out as a run of their own. */
@interface TTIOSamBatchAccumulator : NSObject
@property (nonatomic, readonly) NSMutableArray<NSString *> *sqNames;
@property (nonatomic, readonly) NSMutableArray<TTIOProvenanceRecord *> *provenance;
@property (nonatomic, copy) NSString *rgSample;
@property (nonatomic, copy) NSString *rgPlatform;
@property (nonatomic) int64_t fileMtime;
@property (nonatomic, readonly) NSUInteger readCount;
- (BOOL)consumeLine:(NSString *)line lineNo:(NSUInteger)lineNo error:(NSError **)error;
- (BOOL)consumeLineBytes:(const uint8_t *)line length:(NSUInteger)len
                  lineNo:(NSUInteger)lineNo error:(NSError **)error;
- (TTIOWrittenGenomicRun *)drainBatchWithSampleName:(NSString *)sampleName;
@end

@implementation TTIOSamBatchAccumulator {
    NSMutableArray<NSString *> *_readNames;
    NSMutableArray<NSString *> *_chromosomes;
    NSMutableArray<NSString *> *_cigars;
    NSMutableArray<NSString *> *_mateChromosomes;
    NSMutableArray<NSString *> *_tags;     // M101: columns 12+, "" for none
    NSMutableData *_positionsData;
    NSMutableData *_mappingQualitiesData;
    NSMutableData *_flagsData;
    NSMutableData *_offsetsData;
    NSMutableData *_lengthsData;
    NSMutableData *_matePositionsData;
    NSMutableData *_templateLengthsData;
    NSMutableData *_sequencesData;
    NSMutableData *_qualitiesData;
    uint64_t _runningOffset;
}

- (instancetype)init
{
    self = [super init];
    if (self) {
        _sqNames = [NSMutableArray array];
        _provenance = [NSMutableArray array];
        _rgSample = @"";
        _rgPlatform = @"";
        [self _reset];
    }
    return self;
}

- (void)_reset
{
    _readNames = [NSMutableArray array];
    _chromosomes = [NSMutableArray array];
    _cigars = [NSMutableArray array];
    _mateChromosomes = [NSMutableArray array];
    _tags = [NSMutableArray array];
    _positionsData = [NSMutableData data];
    _mappingQualitiesData = [NSMutableData data];
    _flagsData = [NSMutableData data];
    _offsetsData = [NSMutableData data];
    _lengthsData = [NSMutableData data];
    _matePositionsData = [NSMutableData data];
    _templateLengthsData = [NSMutableData data];
    _sequencesData = [NSMutableData data];
    _qualitiesData = [NSMutableData data];
    _runningOffset = 0;
}

- (NSUInteger)readCount { return _readNames.count; }

/* SAM integers, read straight from the record. NSString's
 * longLongValue takes the leading integer and yields 0 for a field
 * with no digits; this matches that on the shapes SAM permits. */
static int64_t bamFieldInt(const uint8_t *p, NSUInteger n)
{
    NSUInteger i = 0;
    int neg = 0;
    if (i < n && (p[i] == '-' || p[i] == '+')) { neg = (p[i] == '-'); i++; }
    int64_t v = 0;
    for (; i < n && p[i] >= '0' && p[i] <= '9'; i++) v = v * 10 + (p[i] - '0');
    return neg ? -v : v;
}

/* One alignment record, parsed over the bytes samtools wrote.
 *
 * The string path built an NSString for the record, split it into
 * twelve more with characterAtIndex: per character, and copied SEQ and
 * QUAL into their own NSData before appending them. That is about
 * sixteen objects and a message send per byte for every read. This
 * walks the record once, keeps only the fields that are stored as
 * text (M101: the optional tags too), and appends SEQ and QUAL
 * straight from the slice. */
- (BOOL)consumeLineBytes:(const uint8_t *)line length:(NSUInteger)len
                  lineNo:(NSUInteger)lineNo error:(NSError **)error
{
    if (len == 0) return YES;
    if (line[0] == '@') {
        /* Headers are a few hundred lines in a file of millions. */
        NSString *h = [[NSString alloc] initWithBytes:line length:len
                                             encoding:NSUTF8StringEncoding];
        if (h == nil) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"SAM header not valid UTF-8 at line %lu", (unsigned long)lineNo);
            return NO;
        }
        return [self consumeLine:h lineNo:lineNo error:error];
    }

    /* Fields 1-11, trailing optional tags left in the twelfth, exactly
     * as bamSplitTabsLimited(line, 12) cut them. */
    const uint8_t *fp[12];
    NSUInteger fn[12];
    NSUInteger nf = 0, start = 0;
    for (NSUInteger i = 0; i < len; i++) {
        if (line[i] != '\t') continue;
        if (nf + 1 >= 12) break;
        fp[nf] = line + start; fn[nf] = i - start; nf++;
        start = i + 1;
    }
    fp[nf] = line + start; fn[nf] = len - start; nf++;
    if (nf < 11) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
            @"Malformed SAM alignment at line %lu: expected >=11 "
            @"tab-separated fields, got %lu",
            (unsigned long)lineNo, (unsigned long)nf);
        return NO;
    }

    NSString *qname = [[NSString alloc] initWithBytes:fp[0] length:fn[0]
                                             encoding:NSUTF8StringEncoding];
    NSString *rname = [[NSString alloc] initWithBytes:fp[2] length:fn[2]
                                             encoding:NSUTF8StringEncoding];
    NSString *cigar = [[NSString alloc] initWithBytes:fp[5] length:fn[5]
                                             encoding:NSUTF8StringEncoding];
    /* M101: the optional fields, kept as their tab-joined text. */
    NSString *tags = @"";
    if (nf > 11 && fn[11] > 0) {
        tags = [[NSString alloc] initWithBytes:fp[11] length:fn[11]
                                      encoding:NSUTF8StringEncoding];
    }
    if (qname == nil || rname == nil || cigar == nil || tags == nil) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
            @"samtools output not valid UTF-8 for %@ (line %lu)",
            @"this file", (unsigned long)lineNo);
        return NO;
    }
    /* Binding Decision 131: RNEXT "=" expands to RNAME. Sharing the
     * object is also one allocation fewer on the common case. */
    NSString *expandedRnext;
    if (fn[6] == 1 && fp[6][0] == '=') {
        expandedRnext = rname;
    } else {
        expandedRnext = [[NSString alloc] initWithBytes:fp[6] length:fn[6]
                                               encoding:NSUTF8StringEncoding];
        if (expandedRnext == nil) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"samtools output not valid UTF-8 for %@ (line %lu)",
                @"this file", (unsigned long)lineNo);
            return NO;
        }
    }

    uint32_t flag  = (uint32_t)bamFieldInt(fp[1], fn[1]);
    int64_t  pos   = bamFieldInt(fp[3], fn[3]);
    uint8_t  mapq  = (uint8_t)bamFieldInt(fp[4], fn[4]);
    int64_t  pnext = bamFieldInt(fp[7], fn[7]);
    int32_t  tlen  = (int32_t)bamFieldInt(fp[8], fn[8]);

    BOOL seqStar  = (fn[9] == 1 && fp[9][0] == '*');
    BOOL qualStar = (fn[10] == 1 && fp[10][0] == '*');
    NSUInteger sLen = seqStar ? 0 : fn[9];
    NSUInteger qLen = qualStar ? sLen : fn[10];
    /* SEQ/QUAL length-mismatch defence, as the string path had it. */
    if (qLen != sLen) {
        if (seqStar) {
            qLen = 0;
        } else if (!qualStar) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"SEQ/QUAL length mismatch at line %lu: SEQ=%lu QUAL=%lu",
                (unsigned long)lineNo, (unsigned long)sLen, (unsigned long)fn[10]);
            return NO;
        }
    }

    [_readNames addObject:qname];
    [_chromosomes addObject:rname];
    [_cigars addObject:cigar];
    [_mateChromosomes addObject:expandedRnext];
    [_tags addObject:tags];
    [_positionsData appendBytes:&pos length:sizeof(int64_t)];
    [_flagsData appendBytes:&flag length:sizeof(uint32_t)];
    [_mappingQualitiesData appendBytes:&mapq length:sizeof(uint8_t)];
    [_matePositionsData appendBytes:&pnext length:sizeof(int64_t)];
    [_templateLengthsData appendBytes:&tlen length:sizeof(int32_t)];

    uint64_t offset = _runningOffset;
    uint32_t length = (uint32_t)sLen;
    [_offsetsData appendBytes:&offset length:sizeof(uint64_t)];
    [_lengthsData appendBytes:&length length:sizeof(uint32_t)];
    if (sLen > 0) [_sequencesData appendBytes:fp[9] length:sLen];
    if (qualStar) {
        /* QUAL "*" with SEQ present fills 0xFF, without allocating. */
        NSUInteger before = _qualitiesData.length;
        [_qualitiesData increaseLengthBy:qLen];
        if (qLen > 0) memset((uint8_t *)_qualitiesData.mutableBytes + before, 0xFF, qLen);
    } else if (qLen > 0) {
        [_qualitiesData appendBytes:fp[10] length:qLen];
    }
    _runningOffset += sLen;
    return YES;
}

- (BOOL)consumeLine:(NSString *)line lineNo:(NSUInteger)lineNo error:(NSError **)error
{
    if (line.length == 0) return YES;
    if ([line hasPrefix:@"@"]) {
        if ([line hasPrefix:@"@SQ"]) {
            NSDictionary *f = bamParseHeaderFields(line);
            NSString *sn = f[@"SN"];
            if (sn.length > 0) [_sqNames addObject:sn];
        } else if ([line hasPrefix:@"@RG"]) {
            NSDictionary *f = bamParseHeaderFields(line);
            NSString *sm = f[@"SM"];
            NSString *pl = f[@"PL"];
            // First @RG wins per Binding Decision §133.
            if (_rgSample.length == 0 && sm.length > 0) _rgSample = [sm copy];
            if (_rgPlatform.length == 0 && pl.length > 0) _rgPlatform = [pl copy];
        } else if ([line hasPrefix:@"@PG"]) {
            NSDictionary *f = bamParseHeaderFields(line);
            NSString *program = f[@"PN"] ?: @"";
            NSString *commandLine = f[@"CL"];
            NSMutableDictionary *params = [NSMutableDictionary dictionary];
            if (commandLine.length > 0) params[@"CL"] = commandLine;
            for (NSString *k in @[@"ID", @"VN", @"PP"]) {
                NSString *v = f[k];
                if (v.length > 0) params[k] = v;
            }
            TTIOProvenanceRecord *pg = [[TTIOProvenanceRecord alloc]
                initWithInputRefs:@[]
                         software:program
                       parameters:params
                       outputRefs:@[]
                    timestampUnix:_fileMtime];
            [_provenance addObject:pg];
        }
        // @HD, @CO: noted but not mapped to TTI-O fields in v0.
        return YES;
    }

    // Alignment record: fields 1-11, then the optional tags (M101)
    // kept as their tab-joined text in cols[11].
    NSArray<NSString *> *cols = bamSplitTabsLimited(line, 12);
    if (cols.count < 11) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
            @"Malformed SAM alignment at line %lu: expected >=11 "
            @"tab-separated fields, got %lu",
            (unsigned long)lineNo, (unsigned long)cols.count);
        return NO;
    }

    NSString *qname  = cols[0];
    NSString *flagS  = cols[1];
    NSString *rname  = cols[2];
    NSString *posS   = cols[3];
    NSString *mapqS  = cols[4];
    NSString *cigar  = cols[5];
    NSString *rnext  = cols[6];
    NSString *pnextS = cols[7];
    NSString *tlenS  = cols[8];
    NSString *seq    = cols[9];
    NSString *qual   = cols[10];

    uint32_t flag = (uint32_t)[flagS longLongValue];
    int64_t  pos  = (int64_t)[posS  longLongValue];
    uint8_t  mapq = (uint8_t)[mapqS intValue];
    int64_t  pnext = (int64_t)[pnextS longLongValue];
    int32_t  tlen  = (int32_t)[tlenS  intValue];

    // Binding Decision §131: RNEXT "=" expands to RNAME so
    // downstream consumers don't need to remember the convention.
    NSString *expandedRnext = [rnext isEqualToString:@"="] ? rname : rnext;

    [_readNames addObject:qname];
    [_chromosomes addObject:rname];
    [_cigars addObject:cigar];
    [_mateChromosomes addObject:expandedRnext];
    [_tags addObject:(cols.count > 11 ? cols[11] : @"")];
    [_positionsData appendBytes:&pos length:sizeof(int64_t)];
    [_flagsData appendBytes:&flag length:sizeof(uint32_t)];
    [_mappingQualitiesData appendBytes:&mapq length:sizeof(uint8_t)];
    [_matePositionsData appendBytes:&pnext length:sizeof(int64_t)];
    [_templateLengthsData appendBytes:&tlen length:sizeof(int32_t)];

    // SEQ / QUAL handling per HANDOFF §2.5 + §11 §153.
    // SEQ "*" means absent → contributes 0 bytes.
    // QUAL "*" with SEQ present → fill with 0xFF * len(seq).
    // QUAL "*" with SEQ "*" → 0 bytes.
    NSData *seqBytes;
    NSData *qualBytes;
    if ([seq isEqualToString:@"*"]) {
        seqBytes = [NSData data];
    } else {
        seqBytes = [seq dataUsingEncoding:NSASCIIStringEncoding] ?: [NSData data];
    }
    if ([qual isEqualToString:@"*"]) {
        if ([seq isEqualToString:@"*"]) {
            qualBytes = [NSData data];
        } else {
            NSMutableData *q = [NSMutableData dataWithLength:seqBytes.length];
            memset(q.mutableBytes, 0xFF, q.length);
            qualBytes = q;
        }
    } else {
        qualBytes = [qual dataUsingEncoding:NSASCIIStringEncoding] ?: [NSData data];
    }

    // SEQ/QUAL length-mismatch defence (samtools normally
    // validates, but be defensive).
    if (qualBytes.length != seqBytes.length) {
        if ([seq isEqualToString:@"*"]) {
            qualBytes = [NSData data];
        } else if (![qual isEqualToString:@"*"]) {
            if (error) *error = TTIOMakeError(TTIOErrorDatasetRead,
                @"SEQ/QUAL length mismatch at line %lu: "
                @"SEQ=%lu QUAL=%lu",
                (unsigned long)lineNo,
                (unsigned long)seqBytes.length,
                (unsigned long)qualBytes.length);
            return NO;
        }
    }

    uint64_t offset = _runningOffset;
    uint32_t length = (uint32_t)seqBytes.length;
    [_offsetsData appendBytes:&offset length:sizeof(uint64_t)];
    [_lengthsData appendBytes:&length length:sizeof(uint32_t)];
    [_sequencesData appendData:seqBytes];
    [_qualitiesData appendData:qualBytes];
    _runningOffset += length;
    return YES;
}

- (TTIOWrittenGenomicRun *)drainBatchWithSampleName:(NSString *)sampleName
{
    // Effective sample_name: caller override per Binding Decision §133.
    NSString *effSample = (sampleName != nil) ? sampleName : _rgSample;
    NSString *referenceUri = (_sqNames.count > 0) ? _sqNames[0] : @"";
    TTIOWrittenGenomicRun *run = [[TTIOWrittenGenomicRun alloc]
        initWithAcquisitionMode:TTIOAcquisitionModeGenomicWGS
                   referenceUri:referenceUri
                       platform:_rgPlatform
                     sampleName:effSample
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
              signalCompression:TTIOCompressionZlib];
    run.tags = _tags;
    [self _reset];
    return run;
}

@end

@interface TTIOBamReader ()
@property (nonatomic, readwrite, copy) NSArray<TTIOProvenanceRecord *> *provenanceRecords;
@end

@implementation TTIOBamReader

- (instancetype)initWithPath:(NSString *)path
{
    self = [super init];
    if (self) {
        _path = [path copy];
        _provenanceRecords = @[];
    }
    return self;
}

- (nullable TTIOWrittenGenomicRun *)toGenomicRunWithName:(nullable NSString *)name
                                                   region:(nullable NSString *)region
                                               sampleName:(nullable NSString *)sampleName
                                                    error:(NSError **)error
{
    return [self toGenomicRunWithName:name
                                region:region
                            sampleName:sampleName
                              progress:nil
                                 error:error];
}

- (nullable TTIOWrittenGenomicRun *)toGenomicRunWithName:(nullable NSString *)name
                                                   region:(nullable NSString *)region
                                               sampleName:(nullable NSString *)sampleName
                                                 progress:(nullable TTIOProgressBlock)progress
                                                    error:(NSError **)error
{
    (void)name;  // used only for caller-side bookkeeping; not stored on WGR.
    // One batch of unbounded size: the whole file as a single run.
    __block TTIOWrittenGenomicRun *whole = nil;
    BOOL ok = [self iterBatchesWithRegion:region
                               sampleName:sampleName
                               batchReads:NSUIntegerMax
                                 progress:progress
                                    error:error
                               usingBlock:^BOOL(TTIOWrittenGenomicRun *batch, NSError **e) {
        (void)e;
        whole = batch;
        return YES;
    }];
    if (!ok) return nil;
    return whole;
}

- (NSArray<NSString *> *)samtoolsArgumentsForRegion:(NSString *)region error:(NSError **)error
{
    (void)error;
    /* BGZF decode is the subprocess's own bottleneck and is serial
     * without -@. The flag counts threads additional to the main one,
     * and the gain flattens early because this process still has to
     * parse everything the subprocess emits, so ask for a few rather
     * than for the pool's worth. */
    NSUInteger decode = [TTIOThreads resolveImportThreads];
    if (decode > 4) decode = 4;
    NSMutableArray<NSString *> *args = [NSMutableArray arrayWithObjects:
        @"view", @"-h", nil];
    if (decode > 1) {
        [args addObject:@"-@"];
        [args addObject:[NSString stringWithFormat:@"%lu", (unsigned long)(decode - 1)]];
    }
    [args addObject:_path];
    if (region.length > 0) [args addObject:region];
    return args;
}

- (BOOL)iterBatchesWithRegion:(NSString *)region
                   sampleName:(NSString *)sampleName
                   batchReads:(NSUInteger)batchReads
                     progress:(TTIOProgressBlock)progress
                        error:(NSError **)error
                   usingBlock:(BOOL (^)(TTIOWrittenGenomicRun *batch, NSError **error))block
{
    if (progress == nil) progress = TTIOProgressDiscard();
    if (batchReads < 1) batchReads = 1;

    NSString *samtoolsBin = nil;
    if (!bamSamtoolsAvailable(&samtoolsBin, error)) {
        return NO;
    }

    NSFileManager *fm = [NSFileManager defaultManager];
    if (![fm fileExistsAtPath:_path]) {
        if (error) *error = TTIOMakeError(TTIOErrorFileNotFound,
            @"BAM/SAM file not found: %@", _path);
        return NO;
    }
    NSArray<NSString *> *args = [self samtoolsArgumentsForRegion:region error:error];
    if (!args) return NO;

    NSTask *task = [[NSTask alloc] init];
    task.launchPath = samtoolsBin;
    task.arguments = args;
    NSPipe *outPipe = [NSPipe pipe];
    NSPipe *errPipe = [NSPipe pipe];
    task.standardOutput = outPipe;
    task.standardError  = errPipe;

    @try {
        [task launch];
    } @catch (NSException *exc) {
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"failed to launch samtools: %@ (%@)",
            samtoolsBin, exc.reason ?: @"unknown");
        return NO;
    }

    // stderr is drained on a background thread so a chatty samtools
    // never blocks on a full pipe while stdout is being consumed.
    NSFileHandle *errHandle = [errPipe fileHandleForReading];
    NSMutableData *errData = [NSMutableData data];
    NSThread *errThread = [[NSThread alloc] initWithBlock:^{
        NSData *d;
        while ((d = [errHandle availableData]).length > 0) [errData appendData:d];
    }];
    [errThread start];

    TTIOSamBatchAccumulator *acc = [[TTIOSamBatchAccumulator alloc] init];
    {
        struct stat st;
        acc.fileMtime = (stat([_path fileSystemRepresentation], &st) == 0)
            ? (int64_t)st.st_mtime : (int64_t)time(NULL);
    }

    NSFileHandle *outHandle = [outPipe fileHandleForReading];
    NSMutableData *carry = [NSMutableData data];
    NSUInteger lineNo = 0;
    unsigned long long total = 0;
    BOOL ok = YES;
    BOOL anyBatch = NO;
    __block NSError *innerErr = nil;

    BOOL (^emitBatch)(void) = ^BOOL(void) {
        self.provenanceRecords = acc.provenance;
        TTIOWrittenGenomicRun *batch = [acc drainBatchWithSampleName:sampleName];
        NSError *e = nil;
        if (!block(batch, &e)) {
            innerErr = e;
            return NO;
        }
        return YES;
    };

    /* Drain per chunk: samtools output temporaries (and each emitted
     * batch's) otherwise live until return, so memory grows with the
     * file instead of the batch. The error is carried across the pool
     * in a strong local.
     *
     * With threads, header lines feed the shared accumulator on this
     * thread; record lines gather into slices that pool workers parse
     * with their own accumulators (seeded with the header state), and
     * the ordered assembler returns their batches in file order. The
     * caller is both submitter and consumer, so submission never
     * blocks and the window is bounded by pulling before submitting,
     * the same shape as the FASTQ pipeline producer. */
    NSUInteger bamThreads = [TTIOThreads resolveImportThreads];
    TTIOThreadPool *bamPool = bamThreads > 1 ? [TTIOThreadPool poolWithThreads:bamThreads] : nil;
    TTIOOrderedBatchAssembler *bamAsm = bamPool.queue
        ? [[TTIOOrderedBatchAssembler alloc] initWithPool:bamPool] : nil;
    NSUInteger bamWindow = bamThreads + 2;
    __block NSUInteger bamSeq = 0;
    NSUInteger bamPulled = 0;
    BOOL headerDone = NO;
    /* The slice is the bytes samtools wrote, newline terminated;
     * the worker cuts it. Building an NSString per record here only
     * to split it again on the worker cost two passes and an object
     * per read on the thread that has to keep the pipe drained. */
    NSMutableData *sliceBytes = [NSMutableData data];
    __block NSUInteger sliceCount = 0;
    __block unsigned long long totalParallel = 0;

    BOOL (^pullOneBam)(BOOL) = ^BOOL(BOOL emitIt) {
        BOOL done = NO;
        NSError *pe = nil;
        TTIOWrittenGenomicRun *batch = [bamAsm nextBatchWithError:&pe done:&done];
        if (done) return YES;
        if (batch == nil) {
            if (emitIt) innerErr = pe ?: TTIOMakeError(TTIOErrorDatasetRead, @"SAM slice parse failed");
            return NO;
        }
        if (!emitIt) return YES;   /* failure drain: discard, never re-emit */
        self.provenanceRecords = acc.provenance;
        NSError *e = nil;
        if (!block(batch, &e)) { innerErr = e; return NO; }
        return YES;
    };

    NSString *sampleCopy = sampleName;
    BOOL (^submitSlice)(void) = ^BOOL {
        if (sliceCount == 0) return YES;
        NSData *lines = [sliceBytes copy];
        [sliceBytes setLength:0];
        sliceCount = 0;
        NSString *rgS = acc.rgSample, *rgP = acc.rgPlatform;
        int64_t mtime = acc.fileMtime;
        NSArray *sq = [acc.sqNames copy];
        [bamAsm submitSlot:bamSeq producer:^TTIOWrittenGenomicRun *(NSError * __autoreleasing *pe) {
            @autoreleasepool {
                TTIOSamBatchAccumulator *wacc = [[TTIOSamBatchAccumulator alloc] init];
                wacc.rgSample = rgS;
                wacc.rgPlatform = rgP;
                wacc.fileMtime = mtime;
                [wacc.sqNames addObjectsFromArray:sq];
                const uint8_t *lb = lines.bytes;
                NSUInteger ll = lines.length, ls = 0, ln = 0;
                for (NSUInteger i = 0; i < ll; i++) {
                    if (lb[i] != '\n') continue;
                    ln++;
                    NSError *ce = nil;
                    if (![wacc consumeLineBytes:lb + ls length:i - ls
                                         lineNo:ln error:&ce]) {
                        if (pe) *pe = ce;
                        return nil;
                    }
                    ls = i + 1;
                }
                return [wacc drainBatchWithSampleName:sampleCopy];
            }
        }];
        return YES;
    };

    NSError *loopErr = nil;
    BOOL reachedEof = NO;
    while (ok && !reachedEof) {
        @autoreleasepool {
        NSData *chunk = [outHandle availableData];
        BOOL eof = chunk.length == 0;
        if (!eof) [carry appendData:chunk];
        const uint8_t *b = carry.bytes;
        NSUInteger n = carry.length, start = 0;
        for (NSUInteger i = 0; i < n && ok; i++) {
            if (b[i] != '\n') continue;
            NSUInteger end = i;
            if (end > start && b[end - 1] == '\r') end--;
            lineNo++;
            BOOL isHeader = end > start && b[start] == '@';
            if (bamAsm != nil && !isHeader) headerDone = YES;
            if (bamAsm != nil && headerDone) {
                [sliceBytes appendBytes:b + start length:end - start];
                [sliceBytes appendBytes:"\n" length:1];
                sliceCount++;
                start = i + 1;
                totalParallel++;
                if ((totalParallel % TTIOBamReaderProgressIntervalReads) == 0) {
                    progress((int64_t)totalParallel, (int64_t)-1);
                }
                if (sliceCount >= batchReads) {
                    anyBatch = YES;
                    if (!submitSlice()) { ok = NO; break; }
                    bamSeq++;
                    while (ok && bamSeq - bamPulled >= bamWindow) {
                        if (!pullOneBam(YES)) { ok = NO; break; }
                        bamPulled++;
                    }
                }
                continue;
            }
            NSString *line = [[NSString alloc] initWithBytes:b + start length:end - start
                                                    encoding:NSUTF8StringEncoding];
            start = i + 1;
            if (!line) {
                loopErr = TTIOMakeError(TTIOErrorDatasetRead,
                    @"samtools output not valid UTF-8 for %@ (line %lu)", _path, (unsigned long)lineNo);
                ok = NO; break;
            }
            NSUInteger before = acc.readCount;
            NSError *ce = nil;
            if (![acc consumeLine:line lineNo:lineNo error:&ce]) { loopErr = ce; ok = NO; break; }
            if (acc.readCount > before) {
                total++;
                if ((total % TTIOBamReaderProgressIntervalReads) == 0) {
                    progress((int64_t)total, (int64_t)-1);
                }
                if (acc.readCount >= batchReads) {
                    anyBatch = YES;
                    if (!emitBatch()) { ok = NO; break; }
                }
            }
        }
        if (start > 0) [carry replaceBytesInRange:NSMakeRange(0, start) withBytes:NULL length:0];
        if (eof) {
            if (ok && carry.length > 0) {
                NSString *line = [[NSString alloc] initWithData:carry encoding:NSUTF8StringEncoding];
                lineNo++;
                if (bamAsm != nil && headerDone && line.length > 0) {
                    [sliceBytes appendData:carry];
                    [sliceBytes appendBytes:"\n" length:1];
                    sliceCount++;
                    totalParallel++;
                } else {
                    NSUInteger before = acc.readCount;
                    NSError *ce = nil;
                    if (line && ![acc consumeLine:line lineNo:lineNo error:&ce]) { loopErr = ce; ok = NO; }
                    else if (acc.readCount > before) total++;
                }
            }
            reachedEof = YES;
        }
        }
    }
    if (bamAsm != nil) {
        if (ok && sliceCount > 0) {
            anyBatch = YES;
            if (!submitSlice()) ok = NO;
            else bamSeq++;
        }
        [bamAsm finishAfterSlots:bamSeq];
        while (bamPulled < bamSeq) {
            @autoreleasepool {
                if (!pullOneBam(ok)) { ok = NO; }
                bamPulled++;
            }
        }
        [bamPool close];
        total += totalParallel;
    }
    if (!ok && loopErr && error && *error == nil) *error = loopErr;
    if (ok && (acc.readCount > 0 || !anyBatch)) {
        // The remainder, or the empty run of a read-less file.
        if (!emitBatch()) ok = NO;
    }
    [outHandle readDataToEndOfFile];
    [task waitUntilExit];
    while (![errThread isFinished]) [NSThread sleepForTimeInterval:0.001];

    if (ok && task.terminationStatus != 0) {
        NSString *errText = [[NSString alloc] initWithData:errData
                                                  encoding:NSUTF8StringEncoding] ?: @"";
        errText = [errText stringByTrimmingCharactersInSet:
            [NSCharacterSet whitespaceAndNewlineCharacterSet]];
        if (error) *error = TTIOMakeError(TTIOErrorFileOpen,
            @"samtools view exited %d for %@: %@",
            task.terminationStatus, _path,
            errText.length ? errText : @"(no stderr)");
        return NO;
    }
    if (!ok) {
        if (innerErr && error && *error == nil) *error = innerErr;
        return NO;
    }
    self.provenanceRecords = acc.provenance;
    progress((int64_t)total, (int64_t)total);
    return YES;
}

- (TTIOGenomicStreamSource *)streamWithName:(NSString *)name
                                     region:(NSString *)region
                                 sampleName:(NSString *)sampleName
                             referenceFasta:(NSString *)referenceFasta
                             embedReference:(BOOL)embedReference
                                 batchReads:(NSUInteger)batchReads
                                   progress:(TTIOProgressBlock)progress
{
    // The source keeps the reader alive; the reader holds no reference
    // back, so the capture is not a cycle.
    TTIOGenomicBatchProducer producer = ^BOOL(BOOL (^emit)(TTIOWrittenGenomicRun *, NSError **), NSError **error) {
        return [self iterBatchesWithRegion:region sampleName:sampleName batchReads:batchReads
                                  progress:progress error:error usingBlock:emit];
    };
    return [[TTIOGenomicStreamSource alloc] initWithName:name ?: @"genomic_0001"
                                                 batches:producer
                                          referenceFasta:referenceFasta
                                          embedReference:embedReference
                                              blockReads:nil blockBytes:nil
                                   optLegacyWholeChannel:NO];
}

@end
