/*
 * TTIOFeatureFlags.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOFeatureFlags
 * Inherits From: NSObject
 * Declared In:   HDF5/TTIOFeatureFlags.h
 *
 * Versioning + feature-flag utility for the TTI-O file format.
 * Reads and writes the @ttio_format_version + @ttio_features
 * preamble on the root group; detects legacy files (those without
 * the preamble) for backward-compatible decode dispatch.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOFeatureFlags.h"
#import "TTIOHDF5Group.h"

@implementation TTIOFeatureFlags

+ (NSString *)featureBaseV1                  { return @"base_v1"; }
+ (NSString *)featureCompoundIdentifications { return @"compound_identifications"; }
+ (NSString *)featureCompoundQuantifications { return @"compound_quantifications"; }
+ (NSString *)featureCompoundProvenance      { return @"compound_provenance"; }
+ (NSString *)featureCompoundHeaders         { return @"opt_compound_headers"; }
+ (NSString *)featureDatasetEncryption       { return @"opt_dataset_encryption"; }
+ (NSString *)featureNative2DNMR             { return @"opt_native_2d_nmr"; }
+ (NSString *)featureNativeMSImageCube       { return @"opt_native_msimage_cube"; }
+ (NSString *)featureDigitalSignatures       { return @"opt_digital_signatures"; }
+ (NSString *)featureCompoundPerRunProvenance { return @"compound_per_run_provenance"; }
+ (NSString *)featureCanonicalSignatures     { return @"opt_canonical_signatures"; }
+ (NSString *)featureKeyRotation             { return @"opt_key_rotation"; }
+ (NSString *)featureAnonymized              { return @"opt_anonymized"; }
+ (NSString *)featurePQCPreview              { return @"opt_pqc_preview"; }
+ (NSString *)featureMS2ActivationDetail     { return @"opt_ms2_activation_detail"; }
+ (NSString *)featureOptGenomic              { return @"opt_genomic"; }
+ (NSString *)featureOptAssemblyGraph        { return @"opt_assembly_graph"; }
+ (NSString *)featureOptSamTags              { return @"opt_sam_tags"; }
+ (NSString *)featureNoSignalIntDups         { return @"opt_no_signal_int_dups"; }

+ (NSString *)formatVersionForRoot:(TTIOHDF5Group *)root
{
    if ([root hasAttributeNamed:@"ttio_format_version"]) {
        return [root stringAttributeNamed:@"ttio_format_version" error:NULL];
    }
    if ([root hasAttributeNamed:@"ttio_version"]) {
        return [root stringAttributeNamed:@"ttio_version" error:NULL];
    }
    return nil;
}

+ (NSArray<NSString *> *)featuresForRoot:(TTIOHDF5Group *)root
{
    if (![root hasAttributeNamed:@"ttio_features"]) return @[];
    NSString *json = [root stringAttributeNamed:@"ttio_features" error:NULL];
    if (json.length == 0) return @[];
    NSData *data = [json dataUsingEncoding:NSUTF8StringEncoding];
    NSError *err = nil;
    id parsed = [NSJSONSerialization JSONObjectWithData:data options:0 error:&err];
    if (![parsed isKindOfClass:[NSArray class]]) return @[];
    return (NSArray *)parsed;
}

+ (BOOL)root:(TTIOHDF5Group *)root supportsFeature:(NSString *)feature
{
    NSArray *features = [self featuresForRoot:root];
    return [features containsObject:feature];
}

+ (BOOL)isLegacyV1File:(TTIOHDF5Group *)root
{
    return ![root hasAttributeNamed:@"ttio_features"];
}

+ (BOOL)writeFormatVersion:(NSString *)version
                  features:(NSArray<NSString *> *)features
                    toRoot:(TTIOHDF5Group *)root
                     error:(NSError **)error
{
    if (![root setStringAttribute:@"ttio_format_version"
                            value:version
                            error:error]) return NO;

    NSError *jErr = nil;
    NSData *data = [NSJSONSerialization dataWithJSONObject:features
                                                    options:0
                                                      error:&jErr];
    if (!data) {
        if (error) *error = jErr;
        return NO;
    }
    NSString *json = [[NSString alloc] initWithData:data encoding:NSUTF8StringEncoding];
    return [root setStringAttribute:@"ttio_features" value:json error:error];
}

@end
