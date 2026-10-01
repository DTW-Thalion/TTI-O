/*
 * TtioBamDump — canonical-JSON dump of a SAM/BAM/CRAM file for the
 * M87/M88.1 cross-language conformance harness.
 *
 * Usage:
 *     TtioBamDump <bam_or_sam_or_cram_path> [--name NAME] [--reference FASTA]
 *
 * When the input path ends in `.cram` (case-insensitive) the tool
 * dispatches to TTIOCramReader, which requires `--reference <fasta>`.
 * For BAM/SAM input `--reference` is accepted but unused.
 *
 * Reads the file via TTIOBamReader and emits a canonical JSON document
 * on stdout matching the schema documented in HANDOFF.md §7. The same
 * shape is produced by the Python `python -m ttio.importers.bam_dump`
 * CLI and the Java `BamDump` CLI; cross-language tests diff the three
 * outputs to verify field-equality decoding.
 *
 * The JSON keys are sorted alphabetically and the document is indented
 * two spaces with one element per line for arrays — byte-identical to
 * Python's `json.dumps(payload, sort_keys=True, indent=2)` followed by
 * a trailing newline. NSJSONSerialization on GNUstep doesn't expose
 * NSJSONWritingSortedKeys reliably, so we serialise canonically by
 * hand.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#import <Foundation/Foundation.h>
#import "Import/TTIOBamReader.h"
#import "Import/TTIOCramReader.h"
#import "Genomics/TTIOWrittenGenomicRun.h"
#include <openssl/md5.h>

// ── Helpers for canonical JSON emission ─────────────────────────────

static NSString *bdHexFromMD5(NSData *bytes)
{
    unsigned char digest[MD5_DIGEST_LENGTH];
    MD5(bytes.bytes, bytes.length, digest);
    NSMutableString *hex = [NSMutableString stringWithCapacity:32];
    for (int i = 0; i < MD5_DIGEST_LENGTH; i++) {
        [hex appendFormat:@"%02x", digest[i]];
    }
    return hex;
}

// Escape a string per RFC 8259 (the subset Python's json.dumps emits
// for ASCII non-control input). Python's default `ensure_ascii=True`
// escapes any non-ASCII as \uXXXX — for the M87 fixture all values
// are pure ASCII so we never hit that branch, but handle it for
// future safety.
static void bdAppendJsonString(NSMutableString *out, NSString *s)
{
    [out appendString:@"\""];
    NSUInteger n = s.length;
    for (NSUInteger i = 0; i < n; i++) {
        unichar c = [s characterAtIndex:i];
        switch (c) {
            case '\"': [out appendString:@"\\\""]; break;
            case '\\': [out appendString:@"\\\\"]; break;
            case '\b': [out appendString:@"\\b"]; break;
            case '\f': [out appendString:@"\\f"]; break;
            case '\n': [out appendString:@"\\n"]; break;
            case '\r': [out appendString:@"\\r"]; break;
            case '\t': [out appendString:@"\\t"]; break;
            default:
                if (c < 0x20) {
                    [out appendFormat:@"\\u%04x", (unsigned)c];
                } else if (c < 0x7F) {
                    [out appendFormat:@"%c", (char)c];
                } else {
                    // Python ensure_ascii=True default → \uXXXX.
                    [out appendFormat:@"\\u%04x", (unsigned)c];
                }
        }
    }
    [out appendString:@"\""];
}

static NSString *bdIndent(NSUInteger level)
{
    NSMutableString *s = [NSMutableString string];
    for (NSUInteger i = 0; i < level * 2; i++) [s appendString:@" "];
    return s;
}

// Render a string-array as a Python-style multi-line block:
//
//   [
//     "a",
//     "b"
//   ]
//
// Empty arrays render as `[]`.
static void bdAppendStringArray(NSMutableString *out, NSArray<NSString *> *arr,
                                  NSUInteger indentLevel)
{
    if (arr.count == 0) { [out appendString:@"[]"]; return; }
    [out appendString:@"[\n"];
    NSString *itemIndent = bdIndent(indentLevel + 1);
    for (NSUInteger i = 0; i < arr.count; i++) {
        [out appendString:itemIndent];
        bdAppendJsonString(out, arr[i]);
        if (i + 1 < arr.count) [out appendString:@","];
        [out appendString:@"\n"];
    }
    [out appendString:bdIndent(indentLevel)];
    [out appendString:@"]"];
}

// Same shape, but for int64 arrays.
static void bdAppendInt64Array(NSMutableString *out, const int64_t *vals,
                                 NSUInteger n, NSUInteger indentLevel)
{
    if (n == 0) { [out appendString:@"[]"]; return; }
    [out appendString:@"[\n"];
    NSString *itemIndent = bdIndent(indentLevel + 1);
    for (NSUInteger i = 0; i < n; i++) {
        [out appendString:itemIndent];
        [out appendFormat:@"%lld", (long long)vals[i]];
        if (i + 1 < n) [out appendString:@","];
        [out appendString:@"\n"];
    }
    [out appendString:bdIndent(indentLevel)];
    [out appendString:@"]"];
}

static void bdAppendUInt32Array(NSMutableString *out, const uint32_t *vals,
                                  NSUInteger n, NSUInteger indentLevel)
{
    if (n == 0) { [out appendString:@"[]"]; return; }
    [out appendString:@"[\n"];
    NSString *itemIndent = bdIndent(indentLevel + 1);
    for (NSUInteger i = 0; i < n; i++) {
        [out appendString:itemIndent];
        [out appendFormat:@"%u", (unsigned)vals[i]];
        if (i + 1 < n) [out appendString:@","];
        [out appendString:@"\n"];
    }
    [out appendString:bdIndent(indentLevel)];
    [out appendString:@"]"];
}

static void bdAppendUInt8Array(NSMutableString *out, const uint8_t *vals,
                                 NSUInteger n, NSUInteger indentLevel)
{
    if (n == 0) { [out appendString:@"[]"]; return; }
    [out appendString:@"[\n"];
    NSString *itemIndent = bdIndent(indentLevel + 1);
    for (NSUInteger i = 0; i < n; i++) {
        [out appendString:itemIndent];
        [out appendFormat:@"%u", (unsigned)vals[i]];
        if (i + 1 < n) [out appendString:@","];
        [out appendString:@"\n"];
    }
    [out appendString:bdIndent(indentLevel)];
    [out appendString:@"]"];
}

static void bdAppendInt32Array(NSMutableString *out, const int32_t *vals,
                                 NSUInteger n, NSUInteger indentLevel)
{
    if (n == 0) { [out appendString:@"[]"]; return; }
    [out appendString:@"[\n"];
    NSString *itemIndent = bdIndent(indentLevel + 1);
    for (NSUInteger i = 0; i < n; i++) {
        [out appendString:itemIndent];
        [out appendFormat:@"%d", (int)vals[i]];
        if (i + 1 < n) [out appendString:@","];
        [out appendString:@"\n"];
    }
    [out appendString:bdIndent(indentLevel)];
    [out appendString:@"]"];
}

// Render a top-level dict with sorted string keys + 2-space indent.
// Values are dispatched via the @class of the entry.
static void bdAppendKeyString(NSMutableString *out, NSString *key,
                                NSString *val, NSUInteger lvl, BOOL trailingComma)
{
    [out appendString:bdIndent(lvl)];
    bdAppendJsonString(out, key);
    [out appendString:@": "];
    bdAppendJsonString(out, val);
    [out appendString:trailingComma ? @",\n" : @"\n"];
}

static void bdAppendKeyInt(NSMutableString *out, NSString *key,
                             long long val, NSUInteger lvl, BOOL trailingComma)
{
    [out appendString:bdIndent(lvl)];
    bdAppendJsonString(out, key);
    [out appendString:@": "];
    [out appendFormat:@"%lld", val];
    [out appendString:trailingComma ? @",\n" : @"\n"];
}

// ─────────────────────────────────────────────────────────────────────

int main(int argc, const char *argv[])
{
    @autoreleasepool {
        if (argc < 2) {
            fprintf(stderr, "usage: %s <bam_or_sam_or_cram_path> [--name NAME] [--reference FASTA]\n", argv[0]);
            return 2;
        }
        NSString *path = [NSString stringWithUTF8String:argv[1]];
        NSString *runName = @"genomic_0001";
        NSString *reference = nil;
        for (int i = 2; i < argc; i++) {
            NSString *arg = [NSString stringWithUTF8String:argv[i]];
            if ([arg isEqualToString:@"--name"] && i + 1 < argc) {
                runName = [NSString stringWithUTF8String:argv[++i]];
            } else if ([arg isEqualToString:@"--reference"] && i + 1 < argc) {
                reference = [NSString stringWithUTF8String:argv[++i]];
            }
        }

        // Dispatch on extension: .cram → TTIOCramReader (requires --reference);
        // anything else → TTIOBamReader (BAM/SAM). For BAM/SAM, --reference
        // is accepted but ignored (defensive — keeps the CLI tolerant).
        BOOL isCram = [path.lowercaseString hasSuffix:@".cram"];
        if (isCram && reference == nil) {
            fprintf(stderr,
                    "TtioBamDump: --reference <fasta> is required for .cram input\n");
            return 2;
        }

        TTIOBamReader *reader = isCram
            ? [[TTIOCramReader alloc] initWithPath:path referenceFasta:reference]
            : [[TTIOBamReader alloc] initWithPath:path];
        NSError *err = nil;
        TTIOWrittenGenomicRun *run = [reader toGenomicRunWithName:runName
                                                            region:nil
                                                        sampleName:nil
                                                             error:&err];
        if (!run) {
            fprintf(stderr, "TtioBamDump: read failed: %s\n",
                    err.localizedDescription.UTF8String ?: "(unknown)");
            return 1;
        }

        // ── Derive scalar/array views over the parallel channels.
        NSUInteger n = run.readNames.count;
        const int64_t  *positions = run.positionsData.bytes;
        const uint32_t *flags     = run.flagsData.bytes;
        const uint8_t  *mapqs     = run.mappingQualitiesData.bytes;
        const int64_t  *matePos   = run.matePositionsData.bytes;
        const int32_t  *tlens     = run.templateLengthsData.bytes;

        NSString *seqMd5  = bdHexFromMD5(run.sequencesData);
        NSString *qualMd5 = bdHexFromMD5(run.qualitiesData);

        // ── Emit canonical JSON. Keys MUST be alphabetically sorted
        //    to match Python's `json.dumps(..., sort_keys=True)`.
        // Sorted order:
        //   chromosomes, cigars, flags, mapping_qualities,
        //   mate_chromosomes, mate_positions, name, platform,
        //   positions, provenance_count, qualities_md5, read_count,
        //   read_names, reference_uri, sample_name, sequences_md5,
        //   tags, template_lengths
        NSMutableString *out = [NSMutableString string];
        [out appendString:@"{\n"];

        // chromosomes
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"chromosomes");
        [out appendString:@": "];
        bdAppendStringArray(out, run.chromosomes, 1);
        [out appendString:@",\n"];

        // cigars
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"cigars");
        [out appendString:@": "];
        bdAppendStringArray(out, run.cigars, 1);
        [out appendString:@",\n"];

        // flags
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"flags");
        [out appendString:@": "];
        bdAppendUInt32Array(out, flags, n, 1);
        [out appendString:@",\n"];

        // mapping_qualities
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"mapping_qualities");
        [out appendString:@": "];
        bdAppendUInt8Array(out, mapqs, n, 1);
        [out appendString:@",\n"];

        // mate_chromosomes
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"mate_chromosomes");
        [out appendString:@": "];
        bdAppendStringArray(out, run.mateChromosomes, 1);
        [out appendString:@",\n"];

        // mate_positions
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"mate_positions");
        [out appendString:@": "];
        bdAppendInt64Array(out, matePos, n, 1);
        [out appendString:@",\n"];

        // name
        bdAppendKeyString(out, @"name", runName, 1, YES);

        // platform
        bdAppendKeyString(out, @"platform", run.platform ?: @"", 1, YES);

        // positions
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"positions");
        [out appendString:@": "];
        bdAppendInt64Array(out, positions, n, 1);
        [out appendString:@",\n"];

        // provenance_count
        bdAppendKeyInt(out, @"provenance_count",
                        (long long)reader.provenanceRecords.count, 1, YES);

        // qualities_md5
        bdAppendKeyString(out, @"qualities_md5", qualMd5, 1, YES);

        // read_count
        bdAppendKeyInt(out, @"read_count", (long long)n, 1, YES);

        // read_names
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"read_names");
        [out appendString:@": "];
        bdAppendStringArray(out, run.readNames, 1);
        [out appendString:@",\n"];

        // reference_uri
        bdAppendKeyString(out, @"reference_uri", run.referenceUri ?: @"", 1, YES);

        // sample_name
        bdAppendKeyString(out, @"sample_name", run.sampleName ?: @"", 1, YES);

        // sequences_md5
        bdAppendKeyString(out, @"sequences_md5", seqMd5, 1, YES);

        // tags — M101: the SAM optional fields of each read, as
        // samtools prints them ("" when a read has none).
        NSMutableArray<NSString *> *tags = [NSMutableArray arrayWithCapacity:n];
        for (NSUInteger i = 0; i < n; i++) {
            [tags addObject:(i < run.tags.count ? run.tags[i] : @"")];
        }
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"tags");
        [out appendString:@": "];
        bdAppendStringArray(out, tags, 1);
        [out appendString:@",\n"];

        // template_lengths (last — no trailing comma after value array)
        [out appendString:bdIndent(1)];
        bdAppendJsonString(out, @"template_lengths");
        [out appendString:@": "];
        bdAppendInt32Array(out, tlens, n, 1);
        [out appendString:@"\n"];

        [out appendString:@"}\n"];

        fputs(out.UTF8String, stdout);
    }
    return 0;
}
