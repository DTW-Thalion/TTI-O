"""Enumerations mirroring the Objective-C ``TTIOEnums.h`` integer values.

All enums are ``IntEnum`` subclasses because every enum in the format
is persisted on disk as an integer attribute (for example
``@acquisition_mode = 0``). Keeping the Python-side values identical
to the Objective-C ``NS_ENUM`` values makes the HDF5 layer a direct
pass-through and removes any translation table to go stale.

Cross-language equivalents
--------------------------
Objective-C: ``TTIOEnums.h`` ·
Java: ``global.thalion.ttio.Enums``
"""
from __future__ import annotations

from enum import Enum, IntEnum


class SamplingMode(IntEnum):
    """Axis sampling regularity.

    Cross-language: ObjC ``TTIOSamplingMode`` · Java
    ``Enums.SamplingMode``.
    """

    UNIFORM = 0
    NON_UNIFORM = 1


class Precision(IntEnum):
    """Numeric precision of a signal buffer.

    Cross-language: ObjC ``TTIOPrecision`` · Java
    ``Enums.Precision``.

    : ``Precision`` is a pure enum. The HDF5 type
    mapping (``H5T_NATIVE_DOUBLE`` &c.) lives in
    ``ttio.providers.hdf5`` so non-HDF5 providers (SQLite, Memory,
    future Zarr) can import ``Precision`` without pulling HDF5
    dependencies onto their classpath. See v0.6.1 Appendix B gap 7
    in ``docs/api-review-v0.6.md`` for the motivating bug.
    """

    FLOAT32 = 0
    FLOAT64 = 1
    INT32 = 2
    INT64 = 3
    UINT32 = 4
    COMPLEX128 = 5
    UINT8 = 6           # genomic quality scores + packed bases
    UINT16 = 7          # v1.2.0 L1 (Task #82): genomic_index/chromosome_ids
    UINT64 = 9          # genomic index offsets (8 reserved for INT8)

    def numpy_dtype(self) -> str:
        """Return the little-endian NumPy dtype string for this precision."""
        return {
            Precision.FLOAT32: "<f4",
            Precision.FLOAT64: "<f8",
            Precision.INT32: "<i4",
            Precision.INT64: "<i8",
            Precision.UINT32: "<u4",
            Precision.COMPLEX128: "<c16",
            Precision.UINT8: "u1",
            Precision.UINT16: "<u2",
            Precision.UINT64: "<u8",
        }[self]


class Compression(IntEnum):
    """Compression algorithm applied to a signal buffer.

    Cross-language: ObjC ``TTIOCompression`` · Java
    ``Enums.Compression``.
    """

    NONE = 0
    ZLIB = 1
    LZ4 = 2
    NUMPRESS_DELTA = 3
    # genomic codecs (clean-room implementations land in M75+).
    RANS_ORDER0 = 4
    RANS_ORDER1 = 5
    BASE_PACK = 6
    QUALITY_BINNED = 7
    # Slots 8 (NAME_TOKENIZED v1) and 9 (REF_DIFF v1) removed in the
    # v1.0 reset (Phase 2d). v2 successors live at slots 14/15.
    # Slot 10 was reserved for an experimental codec that never shipped.
    # delta + zigzag + varint + rANS order-0 — sorted integer channels.
    DELTA_RANS_ORDER0 = 11
    #.Z: CRAM-mimic rANS-Nx16 quality codec. Static-per-block
    # frequency tables, L=2^15, B=16, N=4, bit-pack 15-bit context.
    # Sole quality codec for the qualities channel; wire magic ``M94Z``.
    FQZCOMP_NX16_Z = 12
    # CRAM-style inline mate-pair encoding.
    # Single blob at signal_channels/mate_info/inline_v2 with
    # @compression = 13. Default mate_info codec from v1.7 onward
    # (no opt-out — v1.0 reset removed the user-facing flag).
    MATE_INLINE_V2 = 13
    # CRAM-style bit-packed sequence diff codec.
    # Single blob at signal_channels/sequences/refdiff_v2 with
    # @compression = 14. Default sequences codec from v1.8 onward
    # (no opt-out — v1.0 reset removed the user-facing flag).
    REF_DIFF_V2 = 14
    # v1.8 #11 ch3: CRAM-style adaptive name-tokenizer codec.
    # Single blob at signal_channels/read_names/name_tok_v2 with
    # @compression = 15. Default read_names codec from v1.8 onward
    # (no opt-out — v1.0 reset removed the user-facing flag).
    NAME_TOKENIZED_V2 = 15
    # Zstandard (RFC 8878). Wire-only today: an opt-in codec for
    # spectral access-unit channels on the transport stream
    # (TransportWriter compression_codec="zstd"); no on-disk
    # @compression dispatch is wired for it yet.
    ZSTD = 16
    # Lossless float64 channel codec: per-block none/delta on the
    # u64 bit view + byte-plane transpose + zstd.
    # On disk: the MS float64 channel default (opt_disable_float_delta
    # to opt out). On the transport wire: the default spectral AU
    # channel codec when TransportWriter(use_compression=True), one
    # FDZ1 stream per channel.
    FLOAT_DELTA_ZSTD = 17
    # M101: SAM optional fields as a tag-line dictionary plus one
    # column per tag key, MD/NM recomputed from the reference. The
    # genomic ``tags`` channel's only codec (docs/codecs/sam_tags.md).
    SAM_TAGS = 18


class ByteOrder(IntEnum):
    """Byte order of a signal buffer on disk.

    Cross-language: ObjC ``TTIOByteOrder`` · Java
    ``Enums.ByteOrder``.
    """

    LITTLE_ENDIAN = 0
    BIG_ENDIAN = 1


class Polarity(IntEnum):
    """Ion polarity for mass spectrometry.

    Cross-language: ObjC ``TTIOPolarity`` · Java
    ``Enums.Polarity``.
    """

    UNKNOWN = 0
    POSITIVE = 1
    NEGATIVE = -1


class ChromatogramType(IntEnum):
    """Chromatogram kind.

    Cross-language: ObjC ``TTIOChromatogramType`` · Java
    ``Enums.ChromatogramType``.
    """

    TIC = 0
    XIC = 1
    SRM = 2


class AcquisitionMode(IntEnum):
    """High-level acquisition scheme for a run.

    Cross-language: ObjC ``TTIOAcquisitionMode`` · Java
    ``Enums.AcquisitionMode``.
    """

    MS1_DDA = 0
    MS2_DDA = 1
    DIA = 2
    SRM = 3
    NMR_1D = 4
    NMR_2D = 5
    IMAGING = 6
    GENOMIC_WGS = 7   # whole-genome sequencing
    GENOMIC_WES = 8   # whole-exome sequencing
    RAMAN = 9         # Raman spectroscopy
    IR = 10           # Infrared spectroscopy
    UV_VIS = 11       # UV-Visible spectroscopy


class IRMode(IntEnum):
    """Infrared y-axis interpretation.

    Cross-language: ObjC ``TTIOIRMode`` · Java
    ``Enums.IRMode``.
    """

    TRANSMITTANCE = 0
    ABSORBANCE = 1


class ImageKind(IntEnum):
    """Discriminator for the shared :class:`ttio.image.Image` base.

    Identifies which concrete imaging modality an :class:`~ttio.image.Image`
    instance is, without inspecting its distinct fields. Purely an in-memory
    aid; not written to the ``.tio`` wire format (each image subclass keeps
    its own on-disk group).
    """

    MS = 0
    RAMAN = 1
    IR = 2


class SpectralAxisKind(IntEnum):
    """Interpretation of an :class:`~ttio.image.Image`'s generic
    ``spectral_axis``.

    ``MZ`` for mass-spectrometry m/z axes, ``WAVENUMBER`` for Raman/IR
    wavenumber axes. In-memory aid only; not part of the wire format.
    """

    MZ = 0
    WAVENUMBER = 1


class SpectrumKind(Enum):
    """Discriminator for a spectrum run's concrete type, derived from the
    persisted ``@spectrum_class`` string (which stays the on-disk source of
    truth — this enum is an in-code dispatch key, P3.8)."""
    MASS = "TTIOMassSpectrum"
    NMR = "TTIONMRSpectrum"
    NMR_2D = "TTIONMR2DSpectrum"
    IR = "TTIOIRSpectrum"
    RAMAN = "TTIORamanSpectrum"
    UVVIS = "TTIOUVVisSpectrum"
    FREE_INDUCTION_DECAY = "TTIOFreeInductionDecay"
    MS_IMAGE_PIXEL = "TTIOMSImagePixel"
    UNKNOWN = ""

    @property
    def persisted(self) -> str:
        return self.value

    @classmethod
    def from_persisted(cls, s):
        if not s:
            return cls.MASS  # v0.1 fallback
        for k in cls:
            if k is not cls.UNKNOWN and k.value == s:
                return k
        return cls.UNKNOWN


class EncryptionLevel(IntEnum):
    """Multi-level protection granularity (MPEG-G style).

    Cross-language: ObjC ``TTIOEncryptionLevel`` · Java
    ``Enums.EncryptionLevel``.
    """

    NONE = 0
    DATASET_GROUP = 1
    DATASET = 2
    DESCRIPTOR_STREAM = 3
    ACCESS_UNIT = 4


class ActivationMethod(IntEnum):
    """MS/MS precursor activation (dissociation) method.

    Stored as ``int32`` in the optional ``activation_methods`` column of
    ``spectrum_index`` (see ``opt_ms2_activation_detail``). ``NONE`` is
    the sentinel for MS1 scans and for MS2+ scans whose activation method
    was not reported by the source instrument.

    Cross-language: ObjC ``TTIOActivationMethod`` · Java
    ``Enums.ActivationMethod``.
    """

    NONE = 0
    CID = 1
    HCD = 2
    ETD = 3
    UVPD = 4
    ECD = 5
    EThcD = 6
