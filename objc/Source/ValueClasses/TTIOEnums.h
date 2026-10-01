#ifndef TTIO_ENUMS_H
#define TTIO_ENUMS_H

#import <Foundation/Foundation.h>

/*
 * TTIOEnums.h
 *
 * Enum declarations shared across the TTI-O Objective-C
 * implementation. Each NS_ENUM declaration order is authoritative;
 * Python's IntEnum values and Java's enum declaration order must
 * match. The HDF5 type mapping (H5T_NATIVE_DOUBLE &c.) for
 * TTIOPrecision lives in the HDF5 provider layer so non-HDF5
 * providers (SQLite, in-memory, Zarr) can use these enums without
 * pulling the HDF5 library onto the link graph.
 *
 * Cross-language equivalents:
 *   Python: ttio.enums
 *   Java:   global.thalion.ttio.Enums
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */

/**
 * Sampling regime for an axis: uniformly-spaced (constant step) or
 * non-uniformly spaced (centroided peak lists, irregular timepoints).
 */
typedef NS_ENUM(NSUInteger, TTIOSamplingMode) {
    TTIOSamplingModeUniform = 0,
    TTIOSamplingModeNonUniform
};

/**
 * Numeric precision of a signal buffer.
 */
typedef NS_ENUM(NSUInteger, TTIOPrecision) {
    TTIOPrecisionFloat32 = 0,
    TTIOPrecisionFloat64,
    TTIOPrecisionInt32,
    TTIOPrecisionInt64,
    TTIOPrecisionUInt32,
    TTIOPrecisionComplex128,
    TTIOPrecisionUInt8,        // Genomic quality scores + packed bases.
    TTIOPrecisionUInt16 = 7,   // Genomic chromosome IDs.
    TTIOPrecisionUInt64 = 9    // Genomic index offsets (8 reserved for INT8).
};

/**
 * Compression / codec applied to a signal channel.
 *
 * Codec ids 0-3 map onto HDF5 filter pipeline stages. Ids 4+ are
 * dedicated per-channel codecs signalled via the per-dataset
 * <code>@compression</code> uint8 attribute. Slots 8 (NAME_TOKENIZED
 * v1), 9 (REF_DIFF v1), and 10 (FQZCOMP_NX16 v1) are reserved on
 * the wire for backward-compatibility ordinal stability; v1.0
 * readers reject them with a migration error.
 */
typedef NS_ENUM(NSUInteger, TTIOCompression) {
    TTIOCompressionNone = 0,
    TTIOCompressionZlib,
    TTIOCompressionLZ4,
    TTIOCompressionNumpressDelta,         // Fixed-point + first-difference (lossy MS m/z).
    TTIOCompressionRansOrder0,            // rANS order-0 entropy coder.
    TTIOCompressionRansOrder1,            // rANS order-1 entropy coder.
    TTIOCompressionBasePack,              // 2-bit ACGT packed bases + sidecar mask.
    TTIOCompressionQualityBinned,         // Illumina-8 quality binning (lossy).
    // Slots 8, 9, 10 are reserved on the wire (removed v1 codecs).
    TTIOCompressionDeltaRansOrder0 = 11,  // Delta + rANS for sortable integer channels.
    TTIOCompressionFqzcompNx16Z = 12,     // CRAM-mimic rANS-Nx16 quality codec (V4 wire format).
    TTIOCompressionMateInlineV2 = 13,     // CRAM-style inline mate-pair codec.
    TTIOCompressionRefDiffV2 = 14,        // Bit-packed reference-diff v2 (substream layout).
    TTIOCompressionNameTokenizedV2 = 15,  // CRAM-style adaptive name-tokenizer v2.
    TTIOCompressionZstd = 16,             // Zstandard (RFC 8878); wire-only today.
    TTIOCompressionFloatDeltaZstd = 17,   // Lossless float64 channels: none/delta + transpose + zstd.
    TTIOCompressionSamTags = 18           // M101: SAM optional fields (docs/codecs/sam_tags.md).
};

/**
 * Byte order for multi-byte numeric values inside a signal buffer.
 * TTI-O canonicalises to little-endian on disk; big-endian is
 * supported only for legacy ingest paths.
 */
typedef NS_ENUM(NSUInteger, TTIOByteOrder) {
    TTIOByteOrderLittleEndian = 0,
    TTIOByteOrderBigEndian
};

/**
 * Polarity of a mass-spectrometry scan. Stored as int32 in the
 * <code>polarities</code> column of <code>spectrum_index</code>;
 * <code>Unknown</code> is the sentinel for non-MS modalities and
 * for scans whose polarity was not reported.
 */
typedef NS_ENUM(NSInteger, TTIOPolarity) {
    TTIOPolarityUnknown  =  0,
    TTIOPolarityPositive =  1,
    TTIOPolarityNegative = -1
};

/**
 * Type of chromatogram trace.
 */
typedef NS_ENUM(NSUInteger, TTIOChromatogramType) {
    TTIOChromatogramTypeTIC = 0,  // Total ion current.
    TTIOChromatogramTypeXIC,      // Extracted ion chromatogram.
    TTIOChromatogramTypeSRM       // Selected reaction monitoring transition.
};

/**
 * Acquisition mode of a run, identifying the instrument or protocol
 * context. Stored on each run group's <code>@acquisition_mode</code>
 * attribute.
 */
typedef NS_ENUM(NSUInteger, TTIOAcquisitionMode) {
    TTIOAcquisitionModeMS1DDA = 0,
    TTIOAcquisitionModeMS2DDA,
    TTIOAcquisitionModeDIA,
    TTIOAcquisitionModeSRM,
    TTIOAcquisitionMode1DNMR,
    TTIOAcquisitionMode2DNMR,
    TTIOAcquisitionModeImaging,
    TTIOAcquisitionModeGenomicWGS,    // Whole-genome sequencing.
    TTIOAcquisitionModeGenomicWES,    // Whole-exome sequencing.
    TTIOAcquisitionModeRaman,         // Raman spectroscopy.
    TTIOAcquisitionModeIR,            // Infrared spectroscopy.
    TTIOAcquisitionModeUVVis          // UV-Visible spectroscopy.
};

/**
 * Granularity at which dataset-level encryption is applied. See the
 * <code>TTIOEncryptable</code> protocol.
 */
typedef NS_ENUM(NSUInteger, TTIOEncryptionLevel) {
    TTIOEncryptionLevelNone = 0,
    TTIOEncryptionLevelDatasetGroup,
    TTIOEncryptionLevelDataset,
    TTIOEncryptionLevelDescriptorStream,
    TTIOEncryptionLevelAccessUnit
};

/**
 * Mid-IR measurement mode. <code>absorbance = -log10(transmittance)</code>;
 * the two are convertible so the exporter / importer preserves
 * whichever the source specified.
 */
typedef NS_ENUM(NSUInteger, TTIOIRMode) {
    TTIOIRModeTransmittance = 0,
    TTIOIRModeAbsorbance
};

/**
 * MS/MS precursor activation (dissociation) method. Stored as int32
 * in the optional <code>activation_methods</code> column of
 * <code>spectrum_index</code> (gated by feature flag
 * <code>opt_ms2_activation_detail</code>). <code>None</code> is the
 * sentinel for MS1 scans and for MS2+ scans whose activation method
 * was not reported by the source instrument.
 *
 * Cross-language equivalents: Python
 * <code>ttio.enums.ActivationMethod</code>, Java
 * <code>global.thalion.ttio.Enums.ActivationMethod</code>.
 */
typedef NS_ENUM(NSInteger, TTIOActivationMethod) {
    TTIOActivationMethodNone  = 0,
    TTIOActivationMethodCID   = 1,
    TTIOActivationMethodHCD   = 2,
    TTIOActivationMethodETD   = 3,
    TTIOActivationMethodUVPD  = 4,
    TTIOActivationMethodECD   = 5,
    TTIOActivationMethodEThcD = 6
};

/**
 * Discriminator for a spectrum run's concrete type. Derived from the
 * persisted <code>@spectrum_class</code> string, which remains the
 * on-disk source of truth; this enum is an in-code dispatch key only
 * (P3.8) and is never written to the <code>.tio</code> / wire format.
 *
 * Use <code>TTIOSpectrumKindFromPersisted()</code> to map a stored
 * string to its member (nil/empty falls back to
 * <code>TTIOSpectrumKindMass</code> for v0.1 files; an unrecognised
 * string yields <code>TTIOSpectrumKindUnknown</code>) and
 * <code>TTIOSpectrumKindPersisted()</code> for the reverse.
 *
 * Cross-language equivalents: Python
 * <code>ttio.enums.SpectrumKind</code>, Java
 * <code>global.thalion.ttio.Enums.SpectrumKind</code>.
 */
typedef NS_ENUM(NSInteger, TTIOSpectrumKind) {
    TTIOSpectrumKindMass = 0,
    TTIOSpectrumKindNMR,
    TTIOSpectrumKindNMR2D,
    TTIOSpectrumKindIR,
    TTIOSpectrumKindRaman,
    TTIOSpectrumKindUVVis,
    TTIOSpectrumKindFreeInductionDecay,
    TTIOSpectrumKindMSImagePixel,
    TTIOSpectrumKindUnknown
};

/**
 * Maps a persisted <code>@spectrum_class</code> string to its
 * <code>TTIOSpectrumKind</code>. A nil or empty string returns
 * <code>TTIOSpectrumKindMass</code> (v0.1 fallback); a string that
 * matches no known member returns <code>TTIOSpectrumKindUnknown</code>.
 */
TTIOSpectrumKind TTIOSpectrumKindFromPersisted(NSString * _Nullable s);

/**
 * Maps a <code>TTIOSpectrumKind</code> back to its canonical persisted
 * <code>@spectrum_class</code> string. <code>TTIOSpectrumKindUnknown</code>
 * maps to the empty string.
 */
NSString *TTIOSpectrumKindPersisted(TTIOSpectrumKind k);

#endif /* TTIO_ENUMS_H */
