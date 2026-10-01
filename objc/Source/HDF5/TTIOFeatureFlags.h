#ifndef TTIO_FEATURE_FLAGS_H
#define TTIO_FEATURE_FLAGS_H

#import <Foundation/Foundation.h>

@class TTIOHDF5File;
@class TTIOHDF5Group;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> HDF5/TTIOFeatureFlags.h</p>
 *
 * <p>Versioning and feature-flag utility for the TTI-O file format.
 * Every TTI-O container carries two root-group attributes:</p>
 *
 * <pre>
 *     @ttio_format_version = "1.0"   (major.minor string)
 *     @ttio_features       = JSON array of feature strings
 * </pre>
 *
 * <p>Feature strings without an <code>opt_</code> prefix are
 * <strong>required</strong>: a reader MUST refuse a file that names
 * a required feature it does not support. Strings with an
 * <code>opt_</code> prefix are informational and MAY be ignored by
 * readers that do not implement the optional capability.</p>
 *
 * <p>Pre-format-version files (those without
 * <code>@ttio_format_version</code> and <code>@ttio_features</code>)
 * are detected by <code>+isLegacyV1File:</code>; the reader falls
 * back to the legacy JSON-attribute paths in that case.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.feature_flags</code><br/>
 * Java: <code>global.thalion.ttio.FeatureFlags</code></p>
 */
@interface TTIOFeatureFlags : NSObject

#pragma mark - Standard feature strings

/** @return <code>@"base_v1"</code> */
+ (NSString *)featureBaseV1;

/** @return <code>@"compound_identifications"</code> */
+ (NSString *)featureCompoundIdentifications;

/** @return <code>@"compound_quantifications"</code> */
+ (NSString *)featureCompoundQuantifications;

/** @return <code>@"compound_provenance"</code> */
+ (NSString *)featureCompoundProvenance;

/** @return <code>@"opt_compound_headers"</code> */
+ (NSString *)featureCompoundHeaders;

/** @return <code>@"opt_dataset_encryption"</code> */
+ (NSString *)featureDatasetEncryption;

/** @return <code>@"opt_native_2d_nmr"</code> */
+ (NSString *)featureNative2DNMR;

/** @return <code>@"opt_native_msimage_cube"</code> */
+ (NSString *)featureNativeMSImageCube;

/** @return <code>@"opt_digital_signatures"</code> */
+ (NSString *)featureDigitalSignatures;

/** @return <code>@"compound_per_run_provenance"</code> */
+ (NSString *)featureCompoundPerRunProvenance;

/** @return <code>@"opt_canonical_signatures"</code> */
+ (NSString *)featureCanonicalSignatures;

/** @return <code>@"opt_key_rotation"</code> */
+ (NSString *)featureKeyRotation;

/** @return <code>@"opt_anonymized"</code> */
+ (NSString *)featureAnonymized;

/** @return <code>@"opt_pqc_preview"</code> — post-quantum signature
 *          / encapsulation envelope preview. */
+ (NSString *)featurePQCPreview;

/** @return <code>@"opt_ms2_activation_detail"</code> — opt-in for
 *          per-spectrum activation method + energy in
 *          <code>spectrum_index</code>. */
+ (NSString *)featureMS2ActivationDetail;

/** @return <code>@"opt_genomic"</code> — file declares genomic
 *          alignment-run content under
 *          <code>/study/genomic_runs/</code>. */
+ (NSString *)featureOptGenomic;

/** @return <code>@"opt_assembly_graph"</code> — file declares
 *          assembly-graph content under
 *          <code>/study/assembly_graphs/</code> (M98,
 *          format-spec 11a). */
+ (NSString *)featureOptAssemblyGraph;

/** @return <code>@"opt_sam_tags"</code> — a genomic run carries SAM
 *          optional fields as a <code>signal_channels/tags</code>
 *          channel (M101, codec SAM_TAGS). */
+ (NSString *)featureOptSamTags;

/** @return <code>@"opt_no_signal_int_dups"</code> — file omits the
 *          legacy integer-channel duplicates from
 *          <code>signal_channels/</code>. */
+ (NSString *)featureNoSignalIntDups;

#pragma mark - Read

/**
 * @param root Root group of the file.
 * @return The <code>@ttio_format_version</code> string, or
 *         <code>nil</code> if absent (legacy file).
 */
+ (NSString *)formatVersionForRoot:(TTIOHDF5Group *)root;

/**
 * @param root Root group of the file.
 * @return The parsed <code>@ttio_features</code> JSON array, or an
 *         empty array if absent.
 */
+ (NSArray<NSString *> *)featuresForRoot:(TTIOHDF5Group *)root;

/**
 * @param root    Root group of the file.
 * @param feature Feature string to test for.
 * @return <code>YES</code> if the feature is listed in
 *         <code>@ttio_features</code>.
 */
+ (BOOL)root:(TTIOHDF5Group *)root supportsFeature:(NSString *)feature;

/**
 * Detects a legacy file that pre-dates the
 * <code>@ttio_format_version</code> + <code>@ttio_features</code>
 * preamble. Such files carry a single <code>@ttio_version</code>
 * attribute instead.
 *
 * @param root Root group of the file.
 * @return <code>YES</code> if neither
 *         <code>@ttio_format_version</code> nor
 *         <code>@ttio_features</code> is present.
 */
+ (BOOL)isLegacyV1File:(TTIOHDF5Group *)root;

#pragma mark - Write

/**
 * Writes the format version + feature list to the root group.
 *
 * @param version  Format-version string (e.g. <code>@"1.0"</code>).
 * @param features Feature string list.
 * @param root     Root group of the file.
 * @param error    Out-parameter populated on failure.
 * @return <code>YES</code> on success.
 */
+ (BOOL)writeFormatVersion:(NSString *)version
                  features:(NSArray<NSString *> *)features
                    toRoot:(TTIOHDF5Group *)root
                     error:(NSError **)error;

@end

#endif /* TTIO_FEATURE_FLAGS_H */
