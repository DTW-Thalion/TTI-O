"""Packet encoding primitives for the TTI-O transport format.

All integers and floats are little-endian on the wire. Strings are
UTF-8 with explicit ``uint16`` or ``uint32`` length prefixes and are
NOT NUL-terminated.
"""
from __future__ import annotations

import struct
import time
from dataclasses import dataclass, field
from enum import IntEnum, IntFlag

HEADER_MAGIC = b"TI"
VERSION = 1
HEADER_SIZE = 24

# Little-endian, 2-byte magic, uint8 version, uint8 packet_type,
# uint16 flags, uint16 dataset_id, uint32 au_sequence, uint32
# payload_length, uint64 timestamp_ns.
_HEADER_FMT = "<2sBBHHIIQ"
_HEADER_STRUCT = struct.Struct(_HEADER_FMT)


class PacketType(IntEnum):
    """Transport packet types. See ``docs/transport-spec.md`` §3.2."""

    STREAM_HEADER = 0x01
    DATASET_HEADER = 0x02
    ACCESS_UNIT = 0x03
    PROTECTION_METADATA = 0x04
    ANNOTATION = 0x05
    PROVENANCE = 0x06
    CHROMATOGRAM = 0x07
    END_OF_DATASET = 0x08
    # Bulk-mode v2-blob carriage (Phase 2c-T, transport-spec §4.10-§4.12).
    # Emitted only when the StreamHeader features list contains
    # "bulk_mode_v2_blobs". One packet per (dataset_id, codec_id) at
    # most, between the matching DatasetHeader and EndOfDataset.
    BLOB_V2_MATE_INFO = 0x09
    BLOB_V2_REF_DIFF = 0x0A
    BLOB_V2_NAME_TOK = 0x0B
    # ---- v0.11 (transport-spec-complete-coverage 2026-05-25) ----
    # See transport-spec §4.13-§4.23. Emitted only when the
    # StreamHeader features list contains "transport_v0_11".
    # Java parity: global.thalion.ttio.transport.PacketType.
    REFERENCE_GROUP_HEADER = 0x10
    REFERENCE_CHROMOSOME   = 0x11
    END_OF_REFERENCE_GROUP = 0x12
    IMAGE_HEADER           = 0x13
    IMAGE_PIXEL            = 0x14
    END_OF_IMAGE           = 0x15
    IDENTIFICATIONS_TABLE  = 0x16
    QUANTIFICATIONS_TABLE  = 0x17
    DATASET_PROVENANCE     = 0x18
    SUBJECT_METADATA       = 0x19
    SAMPLE_METADATA        = 0x1A
    ENCRYPTION_ALGORITHM   = 0x1B
    # ---- M99.1 blocks_v1 per-AU carriage (transport-spec §4.24) ----
    # Emitted only for genomic runs with layout blocks_v1 in an
    # encrypted stream, announced by the StreamHeader feature token
    # "transport_blocks_v1". One GenomicRunSidecar per run after its
    # DatasetHeader, then one BlockSidecar per block before the AUs.
    GENOMIC_RUN_SIDECAR    = 0x1C
    BLOCK_SIDECAR          = 0x1D
    END_OF_STREAM = 0xFF


_KNOWN_PACKET_TYPE_BYTES = frozenset(int(pt) for pt in PacketType)


def is_known_packet_type(type_byte: int) -> bool:
    """Return True iff ``type_byte`` is a defined :class:`PacketType`.

    Used by the reader's forward-compat path (v0.11 task 0.5 / Java
    parity ``PacketType.fromWireOrNull``) to decide whether to skip an
    unrecognised packet rather than fail. The header has already been
    length-prefix-decoded; only the dispatch arm differs.
    """
    return (int(type_byte) & 0xFF) in _KNOWN_PACKET_TYPE_BYTES


# Bulk-mode feature flag — appears in StreamHeader.features when v2
# codec blobs ride on the wire. Receivers without bulk-mode support
# MUST refuse the stream (no opt_ prefix → required).
BULK_MODE_V2_BLOBS_FEATURE = "bulk_mode_v2_blobs"

# v0.11 feature flag — appears in StreamHeader.features when any of
# the 0x10-0x1B packet types ride on the wire. Receivers without
# v0.11 support MUST refuse the stream (no opt_ prefix → required).
TRANSPORT_V0_11_FEATURE = "transport_v0_11"

# M99.1 feature token — appears in StreamHeader.features when any
# genomic run in the stream uses the blocks_v1 layout, so its
# GenomicRunSidecar (0x1C) and BlockSidecar (0x1D) packets ride on
# the wire. Wire-scoped: receivers strip it before writing container
# feature flags. Receivers without sidecar support cannot rebuild a
# restorable container from such a stream (no opt_ prefix → required).
TRANSPORT_BLOCKS_V1_FEATURE = "transport_blocks_v1"
#: Required StreamHeader token (transport-spec v0.13) when a run in the
#: stream has the blocks_v1_grouped layout: its BlockSidecars end with the
#: block's slice of genomic_index/input_index (M103). Wire-scoped, like
#: transport_blocks_v1.
TRANSPORT_BLOCKS_V1_GROUPED_FEATURE = "transport_blocks_v1_grouped"


class PacketFlag(IntFlag):
    ENCRYPTED = 0x0001
    COMPRESSED = 0x0002
    HAS_CHECKSUM = 0x0004
    # set in addition to ENCRYPTED when the AU's semantic header
    # fields are also AES-GCM encrypted. See transport-spec §4.3.3.
    # Readers MUST reject ENCRYPTED_HEADER without ENCRYPTED.
    ENCRYPTED_HEADER = 0x0008


# -------------------------------------------------------------- CRC-32C

_CRC32C_POLY_REFLECTED = 0x82F63B78


def _build_crc32c_table() -> tuple[int, ...]:
    table: list[int] = [0] * 256
    for b in range(256):
        crc = b
        for _ in range(8):
            crc = (crc >> 1) ^ (_CRC32C_POLY_REFLECTED if crc & 1 else 0)
        table[b] = crc
    return tuple(table)


_CRC32C_TABLE = _build_crc32c_table()


def crc32c(data: bytes) -> int:
    """CRC-32C (Castagnoli, reflected) of ``data``.

    Matches ``google-crc32c`` and ``java.util.zip.CRC32C`` output.
    Used when :attr:`PacketFlag.HAS_CHECKSUM` is set.
    """
    crc = 0xFFFFFFFF
    table = _CRC32C_TABLE
    for byte in data:
        crc = (crc >> 8) ^ table[(crc ^ byte) & 0xFF]
    return crc ^ 0xFFFFFFFF


# -------------------------------------------------------------- PacketHeader


@dataclass(frozen=True, slots=True)
class PacketHeader:
    """24-byte transport packet header.

    :attr:`packet_type` is stored as a raw uint8 wire byte (not the
    :class:`PacketType` enum) so headers decoded from forward-compat
    streams whose type byte is outside the enum still round-trip
    cleanly. See :attr:`packet_type_byte` and the v0.11 skip-unknown
    contract (``transport-spec`` §6 / task 0.5).
    """

    packet_type: int
    flags: int = 0
    dataset_id: int = 0
    au_sequence: int = 0
    payload_length: int = 0
    timestamp_ns: int = 0

    @property
    def packet_type_byte(self) -> int:
        """Raw uint8 wire byte for the packet type. Always populated,
        regardless of whether the byte names a known :class:`PacketType`.
        Java parity: :meth:`PacketHeader.packetTypeByte`.
        """
        return int(self.packet_type) & 0xFF

    def to_bytes(self) -> bytes:
        return _HEADER_STRUCT.pack(
            HEADER_MAGIC,
            VERSION,
            int(self.packet_type) & 0xFF,
            int(self.flags) & 0xFFFF,
            int(self.dataset_id) & 0xFFFF,
            int(self.au_sequence) & 0xFFFFFFFF,
            int(self.payload_length) & 0xFFFFFFFF,
            int(self.timestamp_ns) & 0xFFFFFFFFFFFFFFFF,
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> "PacketHeader":
        if len(data) < HEADER_SIZE:
            raise ValueError(
                f"packet header needs {HEADER_SIZE} bytes, got {len(data)}"
            )
        magic, version, ptype, flags, did, aus, plen, ts = _HEADER_STRUCT.unpack_from(
            data, 0
        )
        if magic != HEADER_MAGIC:
            raise ValueError(f"invalid packet magic: {magic!r}")
        if version != VERSION:
            raise ValueError(f"unsupported transport version: {version}")
        return cls(ptype, flags, did, aus, plen, ts)


# -------------------------------------------------------------- ChannelData


_CHANNEL_NAMELEN_STRUCT = struct.Struct("<H")
_CHANNEL_SUFFIX_STRUCT = struct.Struct("<BBII")


@dataclass(slots=True)
class ChannelData:
    """One signal channel inside an :class:`AccessUnit`.

    ``precision`` and ``compression`` are the wire encoding — the
    :attr:`data` bytes are already encoded. The codec does not
    transcode channels; round-tripping preserves the source encoding.
    """

    name: str
    precision: int  # matches ``Precision`` enum
    compression: int  # matches ``Compression`` enum
    n_elements: int
    data: bytes

    def to_bytes(self) -> bytes:
        name_bytes = self.name.encode("utf-8")
        return b"".join((
            _CHANNEL_NAMELEN_STRUCT.pack(len(name_bytes)),
            name_bytes,
            _CHANNEL_SUFFIX_STRUCT.pack(
                int(self.precision) & 0xFF,
                int(self.compression) & 0xFF,
                int(self.n_elements) & 0xFFFFFFFF,
                len(self.data) & 0xFFFFFFFF,
            ),
            self.data,
        ))

    @classmethod
    def from_buffer(cls, buf: bytes, offset: int) -> tuple["ChannelData", int]:
        (name_len,) = _CHANNEL_NAMELEN_STRUCT.unpack_from(buf, offset)
        offset += 2
        name = bytes(buf[offset:offset + name_len]).decode("utf-8")
        offset += name_len
        precision, compression, n_elements, data_length = _CHANNEL_SUFFIX_STRUCT.unpack_from(
            buf, offset
        )
        offset += 10
        data = bytes(buf[offset:offset + data_length])
        offset += data_length
        return (
            cls(name=name, precision=precision, compression=compression,
                n_elements=n_elements, data=data),
            offset,
        )


# -------------------------------------------------------------- AccessUnit

# AU fixed-prefix layout (§4.3): spectrum_class(u8) + acquisition_mode(u8) +
# ms_level(u8) + polarity(u8) + retention_time(f64) + precursor_mz(f64) +
# precursor_charge(u8) + ion_mobility(f64) + base_peak_intensity(f64) +
# n_channels(u8). Total 38 bytes.
_AU_PREFIX_STRUCT = struct.Struct("<BBBBddBddB")
_AU_PIXEL_STRUCT = struct.Struct("<III")
# M89.1 GenomicRead suffix (§4.3.4): position(i64) + mapq(u8) + flags(u16),
# preceded by a uint16-length-prefixed UTF-8 chromosome string. The
# chromosome is variable-length so it can't fit in the fixed prefix —
# this is a deliberate tradeoff for filter-key flexibility (decoy
# contigs, T2T-style names) over byte parsimony.
_AU_GENOMIC_FIXED_STRUCT = struct.Struct("<qBH")
# M90.9 mate extension: appended after the M89.1 fixed suffix when
# present. mate_position(i64) + template_length(i32). Decoder treats
# the extension as optional — payloads ending right after flags
# default these to -1 / 0 (preserves M89.1 wire compatibility).
_AU_GENOMIC_MATE_STRUCT = struct.Struct("<qi")

# Named spectrum_class values. reserved the GenomicRead value;
# the genomic-specific AU payload extension (chromosome, position, mapq,
# flags) ships in M89.1.
SPECTRUM_CLASS_MASS_SPECTRUM = 0
SPECTRUM_CLASS_NMR_SPECTRUM = 1
SPECTRUM_CLASS_NMR_2D = 2
SPECTRUM_CLASS_FID = 3
SPECTRUM_CLASS_MS_IMAGE_PIXEL = 4
SPECTRUM_CLASS_GENOMIC_READ = 5


@dataclass(slots=True)
class AccessUnit:
    """One spectrum as a transport-layer value.

    The fixed-prefix fields are the **filter keys** a server uses for
    selective access. Channel bytes travel in :attr:`channels`.
    """

    spectrum_class: int  # 0=MassSpectrum, 1=NMRSpectrum, 2=NMR2D,
    #                     3=FID, 4=MSImagePixel, 5=GenomicRead ()
    acquisition_mode: int
    ms_level: int
    polarity: int  # wire: 0=positive, 1=negative, 2=unknown
    retention_time: float
    precursor_mz: float
    precursor_charge: int
    ion_mobility: float
    base_peak_intensity: float
    channels: list[ChannelData] = field(default_factory=list)

    # MSImagePixel extension (written only when ``spectrum_class == 4``).
    pixel_x: int = 0
    pixel_y: int = 0
    pixel_z: int = 0

    # GenomicRead extension (written only when ``spectrum_class == 5``).
    # See M89.1 / transport-spec §4.3.4. position is signed i64 to match
    # the BAM convention of -1 for unmapped reads (and to allow
    # T2T-style concatenated assemblies > 4 Gbp). flags is u16 (BAM/SAM
    # convention); mapping_quality is u8 (BAM range 0-255).
    chromosome: str = ""
    position: int = 0
    mapping_quality: int = 0
    flags: int = 0
    # mate extension fields. Optional on the wire — when
    # absent (file or empty AU) they default to BAM unmapped
    # sentinels (-1 mate_position, 0 template_length).
    mate_position: int = -1
    template_length: int = 0

    def to_bytes(self) -> bytes:
        prefix = _AU_PREFIX_STRUCT.pack(
            int(self.spectrum_class) & 0xFF,
            int(self.acquisition_mode) & 0xFF,
            int(self.ms_level) & 0xFF,
            int(self.polarity) & 0xFF,
            float(self.retention_time),
            float(self.precursor_mz),
            int(self.precursor_charge) & 0xFF,
            float(self.ion_mobility),
            float(self.base_peak_intensity),
            len(self.channels) & 0xFF,
        )
        body = prefix + b"".join(ch.to_bytes() for ch in self.channels)
        if self.spectrum_class == 4:
            body += _AU_PIXEL_STRUCT.pack(
                int(self.pixel_x) & 0xFFFFFFFF,
                int(self.pixel_y) & 0xFFFFFFFF,
                int(self.pixel_z) & 0xFFFFFFFF,
            )
        elif self.spectrum_class == 5:
            body += pack_string(self.chromosome, width=2)
            body += _AU_GENOMIC_FIXED_STRUCT.pack(
                int(self.position),
                int(self.mapping_quality) & 0xFF,
                int(self.flags) & 0xFFFF,
            )
            # M90.9 mate extension — always emitted by Python writers
            # post-M90.9. Decoders fall back to defaults when missing
            # so M89.1 fixtures still decode.
            body += _AU_GENOMIC_MATE_STRUCT.pack(
                int(self.mate_position),
                int(self.template_length),
            )
        return body

    @classmethod
    def from_bytes(cls, data: bytes) -> "AccessUnit":
        if len(data) < 38:
            raise ValueError(f"access unit payload too short: {len(data)}")
        (
            spectrum_class, acquisition_mode, ms_level, polarity,
            retention_time, precursor_mz,
            precursor_charge,
            ion_mobility, base_peak_intensity,
            n_channels,
        ) = _AU_PREFIX_STRUCT.unpack_from(data, 0)
        offset = 38
        channels: list[ChannelData] = []
        for _ in range(n_channels):
            ch, offset = ChannelData.from_buffer(data, offset)
            channels.append(ch)
        pixel_x = pixel_y = pixel_z = 0
        chromosome = ""
        position = mapping_quality = flags = 0
        # M90.9 mate extension defaults — match the dataclass
        # defaults so M89.1 AUs decode unchanged.
        mate_position = -1
        template_length = 0
        if spectrum_class == 4:
            if len(data) - offset < 12:
                raise ValueError("MSImagePixel AU missing pixel coordinates")
            pixel_x, pixel_y, pixel_z = _AU_PIXEL_STRUCT.unpack_from(data, offset)
            offset += 12
        elif spectrum_class == 5:
            try:
                chromosome, offset = unpack_string(data, offset, width=2)
            except struct.error as exc:
                raise ValueError(
                    "GenomicRead AU missing chromosome length prefix"
                ) from exc
            if len(data) - offset < _AU_GENOMIC_FIXED_STRUCT.size:
                raise ValueError(
                    "GenomicRead AU missing position/mapq/flags suffix"
                )
            position, mapping_quality, flags = _AU_GENOMIC_FIXED_STRUCT.unpack_from(
                data, offset
            )
            offset += _AU_GENOMIC_FIXED_STRUCT.size
            # M90.9 mate extension — optional. M89.1 payloads end
            # right after flags; M90.9+ payloads carry 12 more bytes.
            if len(data) - offset >= _AU_GENOMIC_MATE_STRUCT.size:
                mate_position, template_length = (
                    _AU_GENOMIC_MATE_STRUCT.unpack_from(data, offset)
                )
                offset += _AU_GENOMIC_MATE_STRUCT.size
        return cls(
            spectrum_class=spectrum_class,
            acquisition_mode=acquisition_mode,
            ms_level=ms_level,
            polarity=polarity,
            retention_time=retention_time,
            precursor_mz=precursor_mz,
            precursor_charge=precursor_charge,
            ion_mobility=ion_mobility,
            base_peak_intensity=base_peak_intensity,
            channels=channels,
            pixel_x=pixel_x,
            pixel_y=pixel_y,
            pixel_z=pixel_z,
            chromosome=chromosome,
            position=position,
            mapping_quality=mapping_quality,
            flags=flags,
            mate_position=mate_position,
            template_length=template_length,
        )


# -------------------------------------------------------------- helpers


def now_ns() -> int:
    return time.time_ns()


def pack_string(s: str, *, width: int = 2) -> bytes:
    """Encode a UTF-8 string with a length prefix (``width=2`` → uint16,
    ``width=4`` → uint32). Matches the ``{uint<N> len, bytes[len]}``
    wire convention used throughout the transport format.
    """
    encoded = s.encode("utf-8")
    if width == 2:
        if len(encoded) > 0xFFFF:
            raise ValueError(f"string too long for uint16 prefix: {len(encoded)}")
        return struct.pack("<H", len(encoded)) + encoded
    if width == 4:
        return struct.pack("<I", len(encoded) & 0xFFFFFFFF) + encoded
    raise ValueError(f"unsupported prefix width: {width}")


def unpack_string(buf: bytes, offset: int, *, width: int = 2) -> tuple[str, int]:
    """Inverse of :func:`pack_string`."""
    if width == 2:
        (length,) = struct.unpack_from("<H", buf, offset)
        offset += 2
    elif width == 4:
        (length,) = struct.unpack_from("<I", buf, offset)
        offset += 4
    else:
        raise ValueError(f"unsupported prefix width: {width}")
    value = bytes(buf[offset:offset + length]).decode("utf-8")
    return value, offset + length


# ----------------------------------------------- Bulk-mode blob payloads
#
# See ``docs/transport-spec.md`` §4.10-§4.12. Each packet ships one
# v2 codec blob verbatim so the receiver writes it byte-for-byte to
# the matching HDF5 path, bypassing the v2 codec encode step.

CODEC_ID_MATE_INLINE_V2 = 13
CODEC_ID_REF_DIFF_V2 = 14
CODEC_ID_NAME_TOKENIZED_V2 = 15


def pack_blob_mate_info(
    *, dataset_id: int, chrom_names: list[str], blob: bytes
) -> bytes:
    """Encode a :data:`PacketType.BLOB_V2_MATE_INFO` payload.

    Mirrors transport-spec §4.10. The chromosome name table is the
    one written alongside ``signal_channels/mate_info/inline_v2`` and
    is required to resolve mate chromosome ids back to strings.
    """
    chrom_table = b"".join(pack_string(n, width=2) for n in chrom_names)
    return b"".join((
        struct.pack("<HBH", dataset_id & 0xFFFF,
                    CODEC_ID_MATE_INLINE_V2, len(chrom_names) & 0xFFFF),
        chrom_table,
        struct.pack("<I", len(blob) & 0xFFFFFFFF),
        bytes(blob),
    ))


def unpack_blob_mate_info(payload: bytes) -> tuple[int, list[str], bytes]:
    """Inverse of :func:`pack_blob_mate_info`. Returns
    ``(dataset_id, chrom_names, blob)``.
    """
    if len(payload) < 5:
        raise ValueError(f"BlobV2MateInfo payload too short: {len(payload)}")
    dataset_id, codec_id, n_names = struct.unpack_from("<HBH", payload, 0)
    if codec_id != CODEC_ID_MATE_INLINE_V2:
        raise ValueError(
            f"BlobV2MateInfo codec_id {codec_id} != {CODEC_ID_MATE_INLINE_V2}"
        )
    offset = 5
    chrom_names: list[str] = []
    for _ in range(n_names):
        n, offset = unpack_string(payload, offset, width=2)
        chrom_names.append(n)
    if len(payload) - offset < 4:
        raise ValueError("BlobV2MateInfo missing blob_length")
    (blob_length,) = struct.unpack_from("<I", payload, offset)
    offset += 4
    if len(payload) - offset != blob_length:
        raise ValueError(
            f"BlobV2MateInfo trailing bytes mismatch: "
            f"declared {blob_length}, actual {len(payload) - offset}"
        )
    blob = bytes(payload[offset:offset + blob_length])
    return dataset_id, chrom_names, blob


def pack_blob_ref_diff(
    *, dataset_id: int, reference_uri: str, blob: bytes
) -> bytes:
    """Encode a :data:`PacketType.BLOB_V2_REF_DIFF` payload (§4.11)."""
    return b"".join((
        struct.pack("<HB", dataset_id & 0xFFFF, CODEC_ID_REF_DIFF_V2),
        pack_string(reference_uri, width=2),
        struct.pack("<I", len(blob) & 0xFFFFFFFF),
        bytes(blob),
    ))


def unpack_blob_ref_diff(payload: bytes) -> tuple[int, str, bytes]:
    """Inverse of :func:`pack_blob_ref_diff`. Returns
    ``(dataset_id, reference_uri, blob)``.
    """
    if len(payload) < 3:
        raise ValueError(f"BlobV2RefDiff payload too short: {len(payload)}")
    dataset_id, codec_id = struct.unpack_from("<HB", payload, 0)
    if codec_id != CODEC_ID_REF_DIFF_V2:
        raise ValueError(
            f"BlobV2RefDiff codec_id {codec_id} != {CODEC_ID_REF_DIFF_V2}"
        )
    offset = 3
    reference_uri, offset = unpack_string(payload, offset, width=2)
    if len(payload) - offset < 4:
        raise ValueError("BlobV2RefDiff missing blob_length")
    (blob_length,) = struct.unpack_from("<I", payload, offset)
    offset += 4
    if len(payload) - offset != blob_length:
        raise ValueError(
            f"BlobV2RefDiff trailing bytes mismatch: "
            f"declared {blob_length}, actual {len(payload) - offset}"
        )
    return dataset_id, reference_uri, bytes(payload[offset:offset + blob_length])


def pack_blob_name_tok(*, dataset_id: int, blob: bytes) -> bytes:
    """Encode a :data:`PacketType.BLOB_V2_NAME_TOK` payload (§4.12)."""
    return b"".join((
        struct.pack("<HB", dataset_id & 0xFFFF, CODEC_ID_NAME_TOKENIZED_V2),
        struct.pack("<I", len(blob) & 0xFFFFFFFF),
        bytes(blob),
    ))


def unpack_blob_name_tok(payload: bytes) -> tuple[int, bytes]:
    """Inverse of :func:`pack_blob_name_tok`. Returns
    ``(dataset_id, blob)``.
    """
    if len(payload) < 7:
        raise ValueError(f"BlobV2NameTok payload too short: {len(payload)}")
    dataset_id, codec_id = struct.unpack_from("<HB", payload, 0)
    if codec_id != CODEC_ID_NAME_TOKENIZED_V2:
        raise ValueError(
            f"BlobV2NameTok codec_id {codec_id} != {CODEC_ID_NAME_TOKENIZED_V2}"
        )
    (blob_length,) = struct.unpack_from("<I", payload, 3)
    offset = 7
    if len(payload) - offset != blob_length:
        raise ValueError(
            f"BlobV2NameTok trailing bytes mismatch: "
            f"declared {blob_length}, actual {len(payload) - offset}"
        )
    return dataset_id, bytes(payload[offset:offset + blob_length])
