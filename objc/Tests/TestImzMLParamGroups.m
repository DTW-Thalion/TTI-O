/*
 * TestImzMLParamGroups — imzML arrays declared through
 * referenceableParamGroups.
 *
 * Many writers (the HR2MSI PXD001283 files among them) state each
 * binaryDataArray's kind, precision and compression in a
 * <referenceableParamGroup> and reference it with
 * <referenceableParamGroupRef>, keeping only the external
 * offset/length cvParams inline. The reader must resolve the reference
 * to know which array an offset belongs to. The fixture declares m/z
 * 64-bit and intensity 32-bit only in the groups, so the precision is
 * checked as well as the array kind.
 *
 * Cross-language equivalents:
 *   python/tests/test_imzml_param_groups.py
 *   java/.../importers/ImzMLParamGroupsTest.java
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#import <Foundation/Foundation.h>
#import "Testing.h"
#import <unistd.h>

#import "Import/TTIOImzMLReader.h"

static NSString *const kPGUUIDHex = @"0123456789abcdef0123456789abcdef";
static const int kPGPixels = 3;

static double pgMzAt(int pixel, int j) { return 100.0 + 50.0 * j + pixel; }
static float pgIntensityAt(int pixel, int j) { return 10.0f * (pixel + 1) + j; }

static NSString *pgArray(NSString *ref, NSUInteger offset, int length, int encoded)
{
    return [NSString stringWithFormat:
        @"<binaryDataArray encodedLength=\"0\">"
        @"<referenceableParamGroupRef ref=\"%@\"/>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000102\" name=\"external offset\" value=\"%lu\"/>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000103\" name=\"external array length\" value=\"%d\"/>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000104\" name=\"external encoded length\" value=\"%d\"/>"
        @"<binary/></binaryDataArray>",
        ref, (unsigned long)offset, length, encoded];
}

static BOOL pgWriteFixture(NSString *imzmlPath, NSString *ibdPath)
{
    NSMutableData *ibd = [NSMutableData data];
    for (NSUInteger i = 0; i < 16; i++) {
        unsigned int byte = 0;
        sscanf([[kPGUUIDHex substringWithRange:NSMakeRange(i * 2, 2)] UTF8String], "%2x", &byte);
        uint8_t b = (uint8_t)byte;
        [ibd appendBytes:&b length:1];
    }

    NSMutableString *spectra = [NSMutableString string];
    for (int p = 0; p < kPGPixels; p++) {
        int n = 3 + p;
        NSUInteger mzOffset = ibd.length;
        for (int j = 0; j < n; j++) {
            double v = pgMzAt(p, j);  // host is little-endian (x86_64 / arm64)
            [ibd appendBytes:&v length:sizeof v];
        }
        NSUInteger intOffset = ibd.length;
        for (int j = 0; j < n; j++) {
            float v = pgIntensityAt(p, j);
            [ibd appendBytes:&v length:sizeof v];
        }
        [spectra appendFormat:
            @"<spectrum id=\"s%d\" index=\"%d\" defaultArrayLength=\"0\">"
            @"<scanList count=\"1\"><scan>"
            @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000050\" name=\"position x\" value=\"%d\"/>"
            @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>"
            @"</scan></scanList>"
            @"<binaryDataArrayList count=\"2\">%@%@</binaryDataArrayList></spectrum>",
            p, p, p + 1,
            pgArray(@"mzArray", mzOffset, n, n * 8),
            pgArray(@"intensityArray", intOffset, n, n * 4)];
    }

    NSString *xml = [NSString stringWithFormat:
        @"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        @"<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1\">"
        @"<fileDescription><fileContent>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000031\" name=\"processed\"/>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000080\" name=\"universally unique identifier\" value=\"{%@}\"/>"
        @"</fileContent></fileDescription>"
        @"<referenceableParamGroupList count=\"2\">"
        @"<referenceableParamGroup id=\"mzArray\">"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>"
        @"</referenceableParamGroup>"
        @"<referenceableParamGroup id=\"intensityArray\">"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\"/>"
        @"<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>"
        @"</referenceableParamGroup>"
        @"</referenceableParamGroupList>"
        @"<scanSettingsList count=\"1\"><scanSettings id=\"ss\">"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000042\" name=\"max count of pixels x\" value=\"%d\"/>"
        @"<cvParam cvRef=\"IMS\" accession=\"IMS:1000043\" name=\"max count of pixels y\" value=\"1\"/>"
        @"</scanSettings></scanSettingsList>"
        @"<run id=\"r\"><spectrumList count=\"%d\">%@</spectrumList></run></mzML>",
        kPGUUIDHex, kPGPixels, kPGPixels, spectra];

    return [xml writeToFile:imzmlPath atomically:YES encoding:NSUTF8StringEncoding error:NULL]
        && [ibd writeToFile:ibdPath atomically:YES];
}

void testImzMLParamGroups(void)
{
    @autoreleasepool {
        NSString *imzml = [NSString stringWithFormat:@"/tmp/ttio_imzml_pg_%d.imzML", (int)getpid()];
        NSString *ibd = [[imzml stringByDeletingPathExtension] stringByAppendingPathExtension:@"ibd"];
        PASS(pgWriteFixture(imzml, ibd), "param-group imzML fixture written");

        NSError *err = nil;
        TTIOImzMLImport *imp = [TTIOImzMLReader readFromImzMLPath:imzml ibdPath:nil error:&err];
        PASS(imp != nil, "param-group imzML parses (%s)",
             err ? [[err localizedDescription] UTF8String] : "ok");
        PASS([imp.mode isEqualToString:@"processed"], "processed mode");
        PASS([imp.uuidHex isEqualToString:kPGUUIDHex], "UUID read");
        PASS(imp.spectra.count == (NSUInteger)kPGPixels, "every pixel read");

        BOOL allMatch = imp.spectra.count == (NSUInteger)kPGPixels;
        for (int p = 0; allMatch && p < kPGPixels; p++) {
            TTIOImzMLPixelSpectrum *s = imp.spectra[p];
            int n = 3 + p;
            if (s.x != p + 1 || s.y != 1 || s.mzCount != (NSUInteger)n
                || s.intensityArray.length != n * sizeof(double)) {
                allMatch = NO;
                break;
            }
            const double *mz = s.mzArray.bytes;
            const double *in = s.intensityArray.bytes;
            for (int j = 0; j < n; j++) {
                if (mz[j] != pgMzAt(p, j) || in[j] != (double)pgIntensityAt(p, j)) {
                    allMatch = NO;
                    break;
                }
            }
        }
        PASS(allMatch, "m/z (64-bit) and intensity (32-bit) arrays resolved via param groups");

        unlink([imzml fileSystemRepresentation]);
        unlink([ibd fileSystemRepresentation]);
    }
}
