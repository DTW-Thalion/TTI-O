/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

/**
 * Enumerations mirroring the ObjC {@code TTIOEnums.h} integer values.
 *
 * <p>Declaration order in each {@code enum} matches the ObjC
 * {@code NS_ENUM} integer values so that {@link Enum#ordinal()} is a
 * direct pass-through of the on-disk attribute. No translation table
 * can go stale.</p>
 *
 * <p><b>API status:</b> Stable.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOEnums.h}, Python {@code ttio.enums}.</p>
 *
 *
 */
public final class Enums {

    private Enums() {}

    /** Numeric precision of a signal buffer.
     *
     * <p><b>:</b> this enum used to hold a
     * {@code HDF5Constants.H5T_NATIVE_*} id per constant, which made
     * the HDF5 JNI wrapper a load-time dependency of every consumer,
     * including non-HDF5 providers such as SQLite. The HDF5 type ids
     * have moved to {@code Hdf5Group.hdf5TypeFor(Precision)} so this
     * enum no longer pulls in the native library.</p>
     */
    public enum Precision {
        /** 32-bit IEEE 754 single-precision float (4 bytes). */
        FLOAT32(4),
        /** 64-bit IEEE 754 double-precision float (8 bytes). */
        FLOAT64(8),
        /** Signed 32-bit integer (4 bytes). */
        INT32(4),
        /** Signed 64-bit integer (8 bytes). */
        INT64(8),
        /** Unsigned 32-bit integer stored as signed {@code int} (4 bytes). */
        UINT32(4),
        /** 128-bit complex: two {@code double} values per element (16 bytes). */
        COMPLEX128(16),
        /** Unsigned 8-bit integer (1 byte).. */
        UINT8(1),
        /**
         * Unsigned 16-bit integer (2 bytes). v1.2.0 L1 (Task #82
         * Phase B.1, 2026-05-01): used by
         * {@code genomic_index/chromosome_ids} after the
         * VL-string-compound chromosomes column was decomposed into a
         * {@code (uint16 ids, compound names)} pair (recovered 42 MB
         * of HDF5 fractal-heap overhead per chr22 .tio file). Ordinal
         * 7 matches Python {@code Precision.UINT16 = 7}.
         */
        UINT16(2),
        /**
         * Reserved for cross-language INT8 parity (ordinal slot 8).
         * @deprecated never use directly; reserved for future extension.
         */
        @Deprecated _RESERVED_INT8(1),
        /**
         * Unsigned 64-bit integer (8 bytes).: genomic index
         * offsets. Cross-language ordinal {@code = 9} matches Python's
         * {@code Precision.UINT64} and ObjC's {@code TTIOPrecisionUInt64}.
         * Wire format uses HDF5 native types.
         */
        UINT64(8);

        private final int elementSize;

        Precision(int elementSize) {
            this.elementSize = elementSize;
        }

        /** @return size in bytes of a single element at this precision. */
        public int elementSize() { return elementSize; }
    }

    /** Compression algorithm applied to a signal buffer. */
    public enum Compression {
        /** No compression; raw binary layout. */
        NONE,
        /** Deflate (zlib/gzip-compatible) compression. */
        ZLIB,
        /** LZ4 block compression. */
        LZ4,
        /** NumPRESS delta-integer compression for ordered m/z or retention-time arrays. */
        NUMPRESS_DELTA,
        /** rANS order-0 entropy coder (genomic CRAM-style). reservation; codec lands in M75. */
        RANS_ORDER0,
        /** rANS order-1 entropy coder (context-aware). reservation; codec lands in M75. */
        RANS_ORDER1,
        /** 2-bit packed nucleotide bases (A/C/G/T). reservation; codec lands in M75. */
        BASE_PACK,
        /** Quality-score binning (Illumina-style). reservation; codec lands in M75. */
        QUALITY_BINNED,
        /** @deprecated Codec id 8 (NAME_TOKENIZED v1) removed in the
         *  v1.0 reset. Slot retained as a placeholder so subsequent
         *  ordinals match the on-disk wire format. */
        @Deprecated
        _RESERVED_8,
        /** @deprecated Codec id 9 (REF_DIFF v1) removed in the v1.0
         *  reset. Slot retained as a placeholder so subsequent ordinals
         *  match the on-disk wire format. */
        @Deprecated
        _RESERVED_9,
        /** @deprecated Codec id 10 reserved — never shipped. */
        @Deprecated
        _RESERVED_10,
        /** Delta + zigzag + varint + rANS order-0 for sorted integer channels
         *  (v1.2, codec id 11). */
        DELTA_RANS_ORDER0,
        /** CRAM-mimic rANS-Nx16 quality codec (M94.Z v1.2, codec id 12). */
        FQZCOMP_NX16_Z,
        /**
         * CRAM-style inline mate-pair encoding (mate_info v2, codec id 13).
         * v1.0 default ON for the {@code signal_channels/mate_info/inline_v2}
         * blob when the native JNI library is available.
         * Cross-language ordinal {@code = 13} matches Python
         * {@code Compression.MATE_INLINE_V2} and ObjC
         * {@code TTIOCompressionMateInlineV2}.
         */
        MATE_INLINE_V2,
        /**
         * Bit-packed REF_DIFF v2 codec (sequences v2, codec id 14).
         * v1.0 default ON for the
         * {@code signal_channels/sequences/refdiff_v2} blob when the
         * native JNI library is available and the run is eligible
         * (single-chromosome, all reads mapped, reference present).
         * Cross-language ordinal {@code = 14} matches Python
         * {@code Compression.REF_DIFF_V2} and ObjC
         * {@code TTIOCompressionRefDiffV2}.
         */
        REF_DIFF_V2,
        /**
         * CRAM-style adaptive name-tokenizer v2 codec (codec id 15).
         * v1.0 default ON for the
         * {@code signal_channels/read_names} blob when the native JNI
         * library is available.
         * Cross-language ordinal {@code = 15} matches Python
         * {@code Compression.NAME_TOKENIZED_V2} and ObjC
         * {@code TTIOCompressionNameTokenizedV2}.
         */
        NAME_TOKENIZED_V2,
        /**
         * Zstandard (RFC 8878), codec id 16. Wire-only today: an
         * opt-in codec for spectral access-unit channels on the
         * transport stream; no on-disk {@code @compression} dispatch
         * is wired for it yet. Cross-language ordinal {@code = 16}
         * matches Python {@code Compression.ZSTD} and ObjC
         * {@code TTIOCompressionZstd}.
         */
        ZSTD,
        /**
         * Lossless float64 channel codec (codec id 17): per-block
         * none/delta on the uint64 bit view + byte-plane transpose +
         * zstd.
         * Cross-language ordinal {@code = 17} matches Python
         * {@code Compression.FLOAT_DELTA_ZSTD} and ObjC
         * {@code TTIOCompressionFloatDeltaZstd}.
         */
        FLOAT_DELTA_ZSTD,
        /**
         * SAM optional fields (codec id 18, M101): a tag-line
         * dictionary plus one column per tag key, MD:Z / NM:i
         * recomputed from the reference. The genomic {@code tags}
         * channel's only codec (docs/codecs/sam_tags.md).
         * Cross-language ordinal {@code = 18} matches Python
         * {@code Compression.SAM_TAGS} and ObjC
         * {@code TTIOCompressionSamTags}.
         */
        SAM_TAGS
    }

    /** Ion polarity for mass spectrometry. */
    public enum Polarity {
        /** Polarity not specified or not applicable. */
        UNKNOWN(0),
        /** Positive-ion mode. */
        POSITIVE(1),
        /** Negative-ion mode. */
        NEGATIVE(-1);

        private final int value;
        Polarity(int value) { this.value = value; }

        /** @return the integer on-disk representation of this polarity. */
        public int intValue() { return value; }

        /**
         * Resolve an integer on-disk value to a {@link Polarity} constant.
         *
         * @param v the integer value read from the file
         * @return the matching constant, or {@link #UNKNOWN} if unrecognised
         */
        public static Polarity fromInt(int v) {
            return switch (v) {
                case 1 -> POSITIVE;
                case -1 -> NEGATIVE;
                default -> UNKNOWN;
            };
        }
    }

    /** Axis sampling regularity. */
    public enum SamplingMode {
        /** Evenly-spaced samples (fixed step size). */
        UNIFORM,
        /** Arbitrarily-spaced samples (coordinate array required). */
        NON_UNIFORM
    }

    /** High-level acquisition scheme for a run. */
    public enum AcquisitionMode {
        /** Data-dependent acquisition, MS1 survey scan. */
        MS1_DDA,
        /** Data-dependent acquisition, MS2 fragmentation scan. */
        MS2_DDA,
        /** Data-independent acquisition. */
        DIA,
        /** Selected reaction monitoring. */
        SRM,
        /** One-dimensional NMR experiment. */
        NMR_1D,
        /** Two-dimensional NMR experiment. */
        NMR_2D,
        /** Mass spectrometry imaging. */
        IMAGING,
        /** Whole-genome sequencing.. */
        GENOMIC_WGS,
        /** Whole-exome sequencing.. */
        GENOMIC_WES,
        /** Raman spectroscopy. */
        RAMAN,
        /** Infrared spectroscopy. */
        IR,
        /** UV-Visible spectroscopy. */
        UV_VIS
    }

    /** Chromatogram kind. */
    public enum ChromatogramType {
        /** Total ion chromatogram. */
        TIC,
        /** Extracted ion chromatogram. */
        XIC,
        /** Selected reaction monitoring chromatogram. */
        SRM
    }

    /** Byte order of a signal buffer on disk. */
    public enum ByteOrder {
        /** Least-significant byte first (x86 native). */
        LITTLE_ENDIAN,
        /** Most-significant byte first (network order). */
        BIG_ENDIAN
    }

    /** Modality of a spatial image cube (MS / Raman / IR). */
    public enum ImageKind {
        /** Imaging mass spectrometry. */
        MS,
        /** Raman hyperspectral imaging. */
        RAMAN,
        /** Infrared hyperspectral imaging. */
        IR
    }

    /** Physical interpretation of an image's spectral axis. */
    public enum SpectralAxisKind {
        /** Mass-to-charge ratio (m/z). */
        MZ,
        /** Wavenumber in {@code cm^-1}. */
        WAVENUMBER
    }

    /** Infrared y-axis interpretation (transmittance vs. absorbance). */
    public enum IRMode {
        /** y-values are transmittance (fraction of incident light). */
        TRANSMITTANCE,
        /** y-values are absorbance (log10 of reciprocal transmittance). */
        ABSORBANCE
    }

    /** Multi-level protection granularity (MPEG-G style). */
    public enum EncryptionLevel {
        /** No encryption applied. */
        NONE,
        /** Protection spans a complete dataset group. */
        DATASET_GROUP,
        /** Protection spans a single dataset. */
        DATASET,
        /** Protection spans a descriptor stream within a dataset. */
        DESCRIPTOR_STREAM,
        /** Protection spans a single access unit (finest granularity). */
        ACCESS_UNIT
    }

    /**
     * MS/MS precursor activation (dissociation) method.
     *
     * <p>Stored as {@code int32} in the optional {@code activation_methods}
     * column of {@code spectrum_index} (gated by feature flag
     * {@code opt_ms2_activation_detail}). {@link #NONE} is the sentinel
     * for MS1 scans and for MS2+ scans whose activation method was not
     * reported by the source instrument.</p>
     */
    public enum ActivationMethod {
        /** No activation (MS1 or not reported). */
        NONE(0),
        /** Collision-induced dissociation. */
        CID(1),
        /** Higher-energy collisional dissociation. */
        HCD(2),
        /** Electron-transfer dissociation. */
        ETD(3),
        /** Ultraviolet photodissociation. */
        UVPD(4),
        /** Electron-capture dissociation. */
        ECD(5),
        /** Electron-transfer and higher-energy collision dissociation. */
        EThcD(6);

        private final int value;
        ActivationMethod(int value) { this.value = value; }

        /** @return the integer on-disk representation of this activation method. */
        public int intValue() { return value; }

        /**
         * Resolve an integer on-disk value to an {@link ActivationMethod}
         * constant.
         *
         * @param v integer value read from the file
         * @return the matching constant, or {@link #NONE} if unrecognised
         */
        public static ActivationMethod fromInt(int v) {
            return switch (v) {
                case 1 -> CID;
                case 2 -> HCD;
                case 3 -> ETD;
                case 4 -> UVPD;
                case 5 -> ECD;
                case 6 -> EThcD;
                default -> NONE;
            };
        }
    }

    /**
     * Kind of spectrum, derived from the persisted {@code @spectrum_class}
     * attribute string.
     *
     * <p>The persisted string remains the source of truth and is written
     * verbatim; this enum exists only to replace stringly-typed dispatch
     * with typed switches. Each member carries its canonical persisted
     * string. {@link #fromPersisted(String)} maps an absent
     * ({@code null}/empty) class to {@link #MASS} (the v0.1 fallback) and
     * any unrecognized string to {@link #UNKNOWN}.</p>
     *
     * <p><b>Cross-language equivalents:</b> Python {@code SpectrumKind}.</p>
     */
    public enum SpectrumKind {
        /** {@code TTIOMassSpectrum} — also the absent/default fallback. */
        MASS("TTIOMassSpectrum"),
        /** {@code TTIONMRSpectrum}. */
        NMR("TTIONMRSpectrum"),
        /** {@code TTIONMR2DSpectrum}. */
        NMR_2D("TTIONMR2DSpectrum"),
        /** {@code TTIOIRSpectrum}. */
        IR("TTIOIRSpectrum"),
        /** {@code TTIORamanSpectrum}. */
        RAMAN("TTIORamanSpectrum"),
        /** {@code TTIOUVVisSpectrum}. */
        UVVIS("TTIOUVVisSpectrum"),
        /** {@code TTIOFreeInductionDecay}. */
        FREE_INDUCTION_DECAY("TTIOFreeInductionDecay"),
        /** {@code TTIOMSImagePixel}. */
        MS_IMAGE_PIXEL("TTIOMSImagePixel"),
        /** Unrecognized persisted class string (forward-compatibility). */
        UNKNOWN(null);

        private final String persisted;

        SpectrumKind(String persisted) {
            this.persisted = persisted;
        }

        /**
         * The canonical persisted {@code @spectrum_class} string for this
         * kind, or {@code null} for {@link #UNKNOWN}.
         */
        public String persisted() {
            return persisted;
        }

        /**
         * Map a persisted {@code @spectrum_class} string to its kind.
         *
         * @param s the stored class string; {@code null} or empty maps to
         *          {@link #MASS} (v0.1 fallback)
         * @return the matching member, or {@link #UNKNOWN} if unrecognized
         */
        public static SpectrumKind fromPersisted(String s) {
            if (s == null || s.isEmpty()) {
                return MASS;
            }
            for (SpectrumKind k : values()) {
                if (s.equals(k.persisted)) {
                    return k;
                }
            }
            return UNKNOWN;
        }
    }
}
