"""Transport-stream writer: :class:`TransportWriter` + writer helpers.

Pure code-movement split of ``codec.py`` (OO-assessment P3.10). The
writer walks a :class:`SpectralDataset` and emits the full packet
sequence specified in ``docs/transport-spec.md``. Bodies are verbatim
from the pre-split ``codec.py``.
"""
from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path
from typing import BinaryIO

import numpy as np

from ..acquisition_run import AcquisitionRun
from ..enums import Compression, Polarity, Precision
from ..mass_spectrum import MassSpectrum
from ..spectral_dataset import SpectralDataset
from ..spectrum import Spectrum
from .._hdf5_io import read_int_attr as io_attr_int  # M90.10 wire codec probe
from .packets import (
    BULK_MODE_V2_BLOBS_FEATURE,
    CODEC_ID_NAME_TOKENIZED_V2,
    HEADER_MAGIC,
    TRANSPORT_V0_11_FEATURE,
    VERSION,
    AccessUnit,
    ChannelData,
    PacketFlag,
    PacketHeader,
    PacketType,
    _AU_PIXEL_STRUCT,
    _AU_PREFIX_STRUCT,
    _CHANNEL_NAMELEN_STRUCT,
    _CHANNEL_SUFFIX_STRUCT,
    _HEADER_STRUCT,
    crc32c,
    now_ns,
    pack_blob_mate_info,
    pack_blob_name_tok,
    pack_blob_ref_diff,
    pack_string,
)
from ._common import (
    _CHECKSUM_STRUCT,
    _POLARITY_TO_WIRE,
    _SPECTRUM_CLASS_TO_WIRE,
    _iter_genomic_run_access_units,
    _read_mate_chrom_names_table,
    genomic_channel_names,
)

#: Writer-side names for the spectral AU channel codecs, in the order
#: a reader gained them. Wire ids: zlib 1, zstd 16, float_delta_zstd 17.
_WIRE_CODECS = ("float_delta_zstd", "zstd", "zlib")


def _wire_channel_encoder(codec: str):
    """Return ``(compression id, encode)`` for a wire codec name, where
    ``encode`` maps a little-endian float64 array to the channel's
    ``data`` bytes."""
    if codec == "float_delta_zstd":
        from ..codecs import float_delta_zstd
        return int(Compression.FLOAT_DELTA_ZSTD), float_delta_zstd.encode
    if codec == "zstd":
        import zstandard
        zc = zstandard.ZstdCompressor(level=3)
        return int(Compression.ZSTD), lambda arr: zc.compress(arr.tobytes())
    if codec == "zlib":
        return int(Compression.ZLIB), lambda arr: zlib.compress(arr.tobytes())
    raise ValueError(f"unsupported wire codec {codec!r}")


def _bulk_carriable(run) -> bool:
    """Whether a genomic run's channel blobs can be carried verbatim:
    every whole-channel run, and a blocks_v1 run with exactly one block
    (its blobs are the whole-channel blobs). Multi-block runs go per-AU,
    and so does a grouped run (M103): its blobs are in stored order while
    its access units are in input order."""
    layout = getattr(run, "layout", "whole")
    if layout == "blocks_v1_grouped":
        return False
    return layout != "blocks_v1" or run.block_count == 1


class TransportWriter:
    """Serialize a :class:`SpectralDataset` as a transport byte stream."""

    def __init__(
        self,
        output: BinaryIO | str | Path,
        *,
        use_checksum: bool = False,
        use_compression: bool = False,
        use_bulk_mode: bool = False,
        compression_codec: str = "float_delta_zstd",
    ):
        """Construct a writer targeting ``output``.

        Parameters
        ----------
        output : BinaryIO, str, or pathlib.Path
            Destination. A path is opened in ``"wb"`` mode and
            closed by :meth:`close`; a file-like object is borrowed
            (not closed by the writer).
        use_checksum : bool, optional
            When ``True``, every emitted packet carries a trailing
            CRC-32C of its payload and the ``HAS_CHECKSUM`` flag is
            set on the header. Default ``False``.
        use_compression : bool, optional
            When ``True``, signal channels in spectral access units
            are compressed before emission with ``compression_codec``.
            The reader detects the per-channel compression byte and
            decompresses transparently. Default ``False``.
        compression_codec : str, optional
            ``"float_delta_zstd"`` (default, wire id 17: the codec-17
            FDZ1 stream per channel, the smallest of the three on
            per-AU float64 payloads), ``"zstd"`` (wire id 16, level 3)
            or ``"zlib"`` (wire id 1). Only consulted when
            ``use_compression`` is ``True``. Readers older than a
            codec's addition reject its id with an
            unsupported-compression error (id 16: pre-1.8.0 readers;
            id 17: 1.8.0 and older), so name ``"zlib"`` explicitly
            until a deployment's readers are current.
        use_bulk_mode : bool, optional
            When ``True``, the writer probes each genomic run for
            v2 codec blobs on disk and emits ``BlobV2*`` packets
            carrying them verbatim, preserving SAM mate sentinels
            byte-for-byte. Adds the
            :data:`BULK_MODE_V2_BLOBS_FEATURE` to the
            ``StreamHeader`` feature list when any blob is emitted.
            Default ``False``.
        """
        self._owns_stream = isinstance(output, (str, Path))
        if self._owns_stream:
            self._stream: BinaryIO = open(output, "wb")  # noqa: SIM115
        else:
            self._stream = output  # type: ignore[assignment]
        self._use_checksum = use_checksum
        self._use_compression = use_compression
        if compression_codec not in _WIRE_CODECS:
            raise ValueError(
                f"unsupported compression_codec {compression_codec!r}; "
                "expected 'float_delta_zstd', 'zstd' or 'zlib'"
            )
        self._compression_codec = compression_codec
        # Phase 2c-T: when True, the writer probes each genomic run for
        # v2 blobs on disk and emits BlobV2* packets carrying them
        # verbatim. The receiver writes the blobs back without going
        # through the v2 codec encode step, preserving SAM mate
        # sentinels (`=`, `""`) byte-for-byte. When the source has no
        # v2 blobs, bulk mode is silently a no-op for that channel.
        self._use_bulk_mode = use_bulk_mode
        self._stream_header_written = False

    @property
    def use_compression(self) -> bool:
        """Whether per-channel zlib compression is enabled for signal data."""
        return self._use_compression

    def __enter__(self) -> "TransportWriter":
        """Return ``self`` so the writer can be used as a context manager."""
        return self

    def __exit__(self, *exc: object) -> None:
        """Close the writer on context exit (delegates to :meth:`close`)."""
        self.close()

    def close(self) -> None:
        """Close the underlying stream if the writer opened it.

        No-op when the caller passed an externally-managed file-like
        object at construction. Safe to call more than once.
        """
        if self._owns_stream and not self._stream.closed:
            self._stream.close()

    def _emit(
        self,
        packet_type: PacketType,
        payload: bytes,
        *,
        dataset_id: int = 0,
        au_sequence: int = 0,
    ) -> None:
        # Inlined from PacketHeader.to_bytes to skip dataclass construct
        # + method dispatch in the per-spectrum hot path. Same layout.
        flags = int(PacketFlag.HAS_CHECKSUM) if self._use_checksum else 0
        header = _HEADER_STRUCT.pack(
            HEADER_MAGIC,
            VERSION,
            int(packet_type) & 0xFF,
            flags & 0xFFFF,
            dataset_id & 0xFFFF,
            au_sequence & 0xFFFFFFFF,
            len(payload) & 0xFFFFFFFF,
            now_ns() & 0xFFFFFFFFFFFFFFFF,
        )
        if self._use_checksum:
            self._stream.write(header + payload + _CHECKSUM_STRUCT.pack(crc32c(payload)))
        else:
            self._stream.write(header + payload)

    def write_stream_header(
        self,
        *,
        format_version: str,
        title: str,
        isa_investigation: str,
        features: list[str],
        n_datasets: int,
    ) -> None:
        """Emit the leading :data:`PacketType.STREAM_HEADER` packet.

        Must be the first packet on the wire (transport-spec §5.4).
        The reader uses :paramref:`features` to enable opt-in wire
        modes such as ``BULK_MODE_V2_BLOBS_FEATURE``.

        Parameters
        ----------
        format_version : str
            Container format version string (e.g. ``"1.2"``).
        title : str
            Free-form container title.
        isa_investigation : str
            ISA-Tab investigation identifier (may be empty).
        features : list of str
            Feature flag strings declaring optional wire-format
            extensions present in the stream.
        n_datasets : int
            Number of dataset blocks (spectral + genomic) the
            stream contains.
        """
        payload = (
            pack_string(format_version, width=2)
            + pack_string(title, width=2)
            + pack_string(isa_investigation, width=2)
            + struct.pack("<H", len(features) & 0xFFFF)
            + b"".join(pack_string(f, width=2) for f in features)
            + struct.pack("<H", n_datasets & 0xFFFF)
        )
        self._emit(PacketType.STREAM_HEADER, payload)
        self._stream_header_written = True

    def write_dataset_header(
        self,
        *,
        dataset_id: int,
        name: str,
        acquisition_mode: int,
        spectrum_class: str,
        channel_names: list[str],
        instrument_json: str,
        expected_au_count: int = 0,
    ) -> None:
        """Emit a :data:`PacketType.DATASET_HEADER` packet.

        One :data:`DATASET_HEADER` precedes every dataset's access
        units. The packet declares the dataset's identity, schema,
        and instrument metadata so the reader can allocate
        per-dataset buffers before AU ingest.

        Parameters
        ----------
        dataset_id : int
            1-based dataset identifier within the stream
            (``uint16``).
        name : str
            Dataset name (run name) as exposed by the source
            container.
        acquisition_mode : int
            Wire encoding of :class:`ttio.AcquisitionMode`.
        spectrum_class : str
            ObjC class name for the spectrum type (e.g.
            ``"TTIOMassSpectrum"``).
        channel_names : list of str
            Ordered signal-channel names.
        instrument_json : str
            JSON-encoded instrument configuration (mass-spectrum
            runs) or genomic-run metadata.
        expected_au_count : int, optional
            Total AU count for the dataset (``uint32``). ``0``
            (default) when unknown.
        """
        payload = (
            struct.pack("<H", dataset_id & 0xFFFF)
            + pack_string(name)
            + struct.pack("<B", int(acquisition_mode) & 0xFF)
            + pack_string(spectrum_class)
            + struct.pack("<B", len(channel_names) & 0xFF)
            + b"".join(pack_string(c) for c in channel_names)
            + pack_string(instrument_json, width=4)
            + struct.pack("<I", expected_au_count & 0xFFFFFFFF)
        )
        self._emit(PacketType.DATASET_HEADER, payload, dataset_id=dataset_id)

    def write_access_unit(
        self,
        *,
        dataset_id: int,
        au_sequence: int,
        au: AccessUnit,
    ) -> None:
        """Emit one :data:`PacketType.ACCESS_UNIT` packet.

        The :class:`AccessUnit` is serialised via
        :meth:`AccessUnit.to_bytes`. The hot per-spectrum path is in
        :meth:`_emit_run_access_units`; this method is the
        general-purpose entry point used by genomic runs and by
        callers driving emission manually.

        Parameters
        ----------
        dataset_id : int
            Owning dataset id (matches the prior
            ``DATASET_HEADER``).
        au_sequence : int
            0-based monotonically increasing AU index within the
            dataset.
        au : AccessUnit
            The unit to emit.
        """
        self._emit(
            PacketType.ACCESS_UNIT,
            au.to_bytes(),
            dataset_id=dataset_id,
            au_sequence=au_sequence,
        )

    def write_blob_v2_mate_info(
        self, *, dataset_id: int, chrom_names: list[str], blob: bytes
    ) -> None:
        """Emit a :data:`PacketType.BLOB_V2_MATE_INFO` packet (§4.10)."""
        self._emit(
            PacketType.BLOB_V2_MATE_INFO,
            pack_blob_mate_info(
                dataset_id=dataset_id, chrom_names=chrom_names, blob=blob,
            ),
            dataset_id=dataset_id,
        )

    def write_blob_v2_ref_diff(
        self, *, dataset_id: int, reference_uri: str, blob: bytes
    ) -> None:
        """Emit a :data:`PacketType.BLOB_V2_REF_DIFF` packet (§4.11)."""
        self._emit(
            PacketType.BLOB_V2_REF_DIFF,
            pack_blob_ref_diff(
                dataset_id=dataset_id,
                reference_uri=reference_uri,
                blob=blob,
            ),
            dataset_id=dataset_id,
        )

    def write_blob_v2_name_tok(
        self, *, dataset_id: int, blob: bytes
    ) -> None:
        """Emit a :data:`PacketType.BLOB_V2_NAME_TOK` packet (§4.12)."""
        self._emit(
            PacketType.BLOB_V2_NAME_TOK,
            pack_blob_name_tok(dataset_id=dataset_id, blob=blob),
            dataset_id=dataset_id,
        )

    # ----------------------------------------------- v0.11 §4.21

    def write_dataset_provenance(self, records) -> None:
        """v0.11 Task 2.5: emit a
        :data:`PacketType.DATASET_PROVENANCE` (0x18) packet carrying
        the dataset-level provenance chain (format-spec §6.3).

        One packet carries all records. Wire layout per transport-spec
        §4.21::

            record_count:        uint32
            # repeated record_count times:
            timestamp_unix:      int64
            software_length:     uint16, software bytes[..]   (UTF-8)
            parameters_length:   uint16, parameters_json[..]  (UTF-8 JSON)
            input_refs_length:   uint16, input_refs_csv[..]   (UTF-8 CSV)
            output_refs_length:  uint16, output_refs_csv[..]  (UTF-8 CSV)

        All multi-byte integers LITTLE-ENDIAN per spec §1.7. The
        input_refs / output_refs lists ride as comma-joined UTF-8 — a
        single empty string for an empty list (no separators).

        Distinct from the per-run ``Provenance`` (0x06) packet which
        carries one JSON record per packet.

        Java parity: :meth:`TransportWriter.writeDatasetProvenance`
        (commit ``563e09c3``). Per-field empty handling matches Java:
        empty parameters render as ``"{}"`` (not omitted); empty refs
        render as ``""`` (the empty-join result).

        :param records: iterable of :class:`ttio.provenance.ProvenanceRecord`.
            ``None`` raises ``ValueError``; an empty list is a no-op
            (no packet emitted) per spec §5.4 "zero or more".
        """
        if records is None:
            raise ValueError(
                "write_dataset_provenance: records must not be None"
            )
        records = list(records)
        if not records:
            return
        # Pre-compute UTF-8 byte arrays so we can size the buffer
        # exactly. Mirrors the StreamHeader/DatasetHeader emit pattern.
        software_bytes: list[bytes] = []
        params_bytes: list[bytes] = []
        inputs_bytes: list[bytes] = []
        outputs_bytes: list[bytes] = []
        total = 4  # record_count
        for r in records:
            sb = (r.software or "").encode("utf-8")
            pb = _provenance_params_json(r.parameters).encode("utf-8")
            ib = _provenance_csv_join(r.input_refs).encode("utf-8")
            ob = _provenance_csv_join(r.output_refs).encode("utf-8")
            for b in (sb, pb, ib, ob):
                if len(b) > 0xFFFF:
                    raise ValueError(
                        f"DATASET_PROVENANCE: per-field length "
                        f"{len(b)} exceeds uint16 max"
                    )
            software_bytes.append(sb)
            params_bytes.append(pb)
            inputs_bytes.append(ib)
            outputs_bytes.append(ob)
            total += (
                8                       # timestamp_unix
                + 2 + len(sb)
                + 2 + len(pb)
                + 2 + len(ib)
                + 2 + len(ob)
            )
        chunks: list[bytes] = [struct.pack("<I", len(records) & 0xFFFFFFFF)]
        for i, r in enumerate(records):
            chunks.append(struct.pack("<q", int(r.timestamp_unix)))
            chunks.append(struct.pack("<H", len(software_bytes[i]) & 0xFFFF))
            chunks.append(software_bytes[i])
            chunks.append(struct.pack("<H", len(params_bytes[i]) & 0xFFFF))
            chunks.append(params_bytes[i])
            chunks.append(struct.pack("<H", len(inputs_bytes[i]) & 0xFFFF))
            chunks.append(inputs_bytes[i])
            chunks.append(struct.pack("<H", len(outputs_bytes[i]) & 0xFFFF))
            chunks.append(outputs_bytes[i])
        payload = b"".join(chunks)
        assert len(payload) == total, (
            f"DATASET_PROVENANCE size mismatch: predicted {total}, "
            f"actual {len(payload)}"
        )
        self._emit(PacketType.DATASET_PROVENANCE, payload)

    # ----------------------------------------------- v0.11 §4.23

    def write_encryption_algorithm(self, algorithm: str) -> None:
        """v0.11 Task 2.4: emit an
        :data:`PacketType.ENCRYPTION_ALGORITHM` (0x1B) packet carrying
        the dataset-level ``@encrypted`` algorithm name
        (e.g. ``"aes-256-gcm"``).

        Wire layout per transport-spec §4.23::

            algorithm_length:  uint16
            algorithm_utf8:    bytes[algorithm_length]

        All multi-byte integers are LITTLE-ENDIAN (spec §1.7). This
        packet only conveys the algorithm-name string; per-AU key
        material continues to ride on ``ProtectionMetadata`` (0x04).

        Java parity: :meth:`TransportWriter.writeEncryptionAlgorithm`
        (commit ``530a5833``).
        """
        if algorithm is None:
            raise ValueError(
                "write_encryption_algorithm: algorithm must not be None"
            )
        algo_bytes = algorithm.encode("utf-8")
        if len(algo_bytes) > 0xFFFF:
            raise ValueError(
                f"ENCRYPTION_ALGORITHM: algorithm name {len(algo_bytes)} "
                "bytes exceeds uint16 max"
            )
        payload = (
            struct.pack("<H", len(algo_bytes) & 0xFFFF)
            + algo_bytes
        )
        self._emit(PacketType.ENCRYPTION_ALGORITHM, payload)

    # ----------------------------------------------- v0.11 §4.13-§4.15

    #: Threshold below which a chromosome rides as raw UINT8
    #: (encoding=0). Mirrors transport-spec §4.14: ZLIB framing costs
    #: dominate short sequences, so the writer skips compression below
    #: 4 KiB and lets the reader handle both encodings. Java parity:
    #: ``TransportWriter.REFERENCE_CHROMOSOME_ZLIB_THRESHOLD``.
    REFERENCE_CHROMOSOME_ZLIB_THRESHOLD = 4096

    def write_reference_group(self, ref) -> None:
        """v0.11 Stage 1: emit a :class:`ReferenceImport` as the packet
        sequence
        ``REFERENCE_GROUP_HEADER (0x10) -> N x REFERENCE_CHROMOSOME (0x11)
        -> END_OF_REFERENCE_GROUP (0x12)``.

        Wire layout matches transport-spec §4.13-§4.15. All multi-byte
        integers are LITTLE-ENDIAN (spec §1.7). The chromosome index
        rides in the packet header's ``au_sequence`` field (0-based).
        The MD5 hex string from ``ReferenceImport.md5.hex()`` is
        emitted verbatim as 32 ASCII bytes.

        The encoding byte on each chromosome record is 0 (uncompressed
        UINT8) when the raw sequence is shorter than
        :attr:`REFERENCE_CHROMOSOME_ZLIB_THRESHOLD`, otherwise 1
        (zlib via :func:`zlib.compress` with default settings).

        Reader-side materialisation is added by Task 2.3; this method
        only emits the wire bytes. Java parity:
        :meth:`TransportWriter.writeReferenceGroup`
        (commit ``622aa8bd``).

        :param ref: :class:`ttio.genomic.reference_import.ReferenceImport`
            to emit. ``ref.md5`` must be the 16-byte digest the
            constructor populates by default.
        """
        chrom_names = ref.chromosomes
        seqs = ref.sequences
        chrom_count = len(chrom_names)
        total_bases = ref.total_bases
        md5_hex = ref.md5.hex()
        if len(md5_hex) != 32:
            raise ValueError(
                f"ReferenceImport.md5.hex() must be 32 hex chars, got "
                f"{len(md5_hex)}"
            )

        # -- REFERENCE_GROUP_HEADER (0x10) ------------------------------
        uri_bytes = ref.uri.encode("utf-8")
        md5_hex_bytes = md5_hex.encode("ascii")
        header_payload = b"".join((
            struct.pack("<H", len(uri_bytes) & 0xFFFF),
            uri_bytes,
            struct.pack("<I", chrom_count & 0xFFFFFFFF),
            struct.pack("<Q", total_bases & 0xFFFFFFFFFFFFFFFF),
            md5_hex_bytes,
        ))
        self._emit(PacketType.REFERENCE_GROUP_HEADER, header_payload)

        # -- REFERENCE_CHROMOSOME (0x11) — one per contig ---------------
        for i, name in enumerate(chrom_names):
            seq = seqs[i]
            name_bytes = name.encode("utf-8")
            if len(seq) < self.REFERENCE_CHROMOSOME_ZLIB_THRESHOLD:
                encoding = 0
                payload_bytes = seq
            else:
                encoding = 1
                payload_bytes = zlib.compress(seq)
            chrom_payload = b"".join((
                struct.pack("<H", len(name_bytes) & 0xFFFF),
                name_bytes,
                struct.pack("<Q", len(seq) & 0xFFFFFFFFFFFFFFFF),
                struct.pack("<B", encoding & 0xFF),
                struct.pack("<I", len(payload_bytes) & 0xFFFFFFFF),
                bytes(payload_bytes),
            ))
            self._emit(
                PacketType.REFERENCE_CHROMOSOME,
                chrom_payload,
                au_sequence=i,
            )

        # -- END_OF_REFERENCE_GROUP (0x12) ------------------------------
        self._emit(
            PacketType.END_OF_REFERENCE_GROUP,
            struct.pack("<I", chrom_count & 0xFFFFFFFF),
        )

    # ------------------------------------------------ v0.11 §4.16-§4.18
    def _emit_image_header(
        self,
        *,
        modality: int,
        width: int,
        height: int,
        bins: int,
        pixel_size_x: float,
        pixel_size_y: float,
        scan_pattern_byte: int,
        axis_kind: int,
        axis: np.ndarray,
        is_continuous: int,
        title: str,
        isa_id: str,
        extras: bytes,
    ) -> None:
        """v0.11 Task 5.3: shared IMAGE_HEADER (0x13) packing routine.

        Mirrors Java's ``emitImageHeader`` helper so the common
        header shape stays byte-stable across modalities (MS / Raman
        / IR) and the ``modality_extras`` slot is appended once per
        call. The continuous-mode bit and the modality-specific tail
        come from the caller.
        """
        axis = np.asarray(axis, dtype=np.float64)
        title_bytes = (title or "").encode("utf-8")
        isa_bytes = (isa_id or "").encode("utf-8")
        if len(title_bytes) > 0xFFFF:
            raise ValueError(
                f"IMAGE_HEADER: title {len(title_bytes)} bytes exceeds "
                f"uint16 max"
            )
        if len(isa_bytes) > 0xFFFF:
            raise ValueError(
                f"IMAGE_HEADER: isa_id {len(isa_bytes)} bytes exceeds "
                f"uint16 max"
            )
        if len(extras) > 0xFFFF:
            raise ValueError(
                f"IMAGE_HEADER: modality_extras {len(extras)} bytes "
                f"exceeds uint16 max"
            )
        header_payload = b"".join((
            struct.pack("<B", modality & 0xFF),
            struct.pack("<I", width & 0xFFFFFFFF),
            struct.pack("<I", height & 0xFFFFFFFF),
            struct.pack("<I", bins & 0xFFFFFFFF),
            struct.pack("<d", float(pixel_size_x)),
            struct.pack("<d", float(pixel_size_y)),
            struct.pack("<B", scan_pattern_byte & 0xFF),
            struct.pack("<B", axis_kind & 0xFF),
            struct.pack("<I", int(axis.size) & 0xFFFFFFFF),
            axis.astype("<f8", copy=False).tobytes(),
            struct.pack("<B", is_continuous & 0xFF),
            struct.pack("<H", len(title_bytes) & 0xFFFF),
            title_bytes,
            struct.pack("<H", len(isa_bytes) & 0xFFFF),
            isa_bytes,
            struct.pack("<H", len(extras) & 0xFFFF),
            bytes(extras),
        ))
        self._emit(PacketType.IMAGE_HEADER, header_payload)

    def write_image(self, image) -> None:
        """v0.11 Task 2.6: emit an :class:`~ttio.MSImage` as the packet
        sequence
        ``IMAGE_HEADER (0x13) -> N x IMAGE_PIXEL (0x14)
        -> END_OF_IMAGE (0x15)``, where ``N = width * height``.

        Wire layout matches transport-spec §4.16-§4.18 byte-for-byte
        with Java's :meth:`TransportWriter.writeImage` (commit
        ``a6b1e5d9``). All multi-byte integers are LITTLE-ENDIAN
        (spec §1.7).

        Continuous-mode only — the shared m/z axis rides on the
        IMAGE_HEADER and every IMAGE_PIXEL carries only its
        intensities (FLOAT64, uncompressed). For sparse cubes the
        opt-in :meth:`write_image_processed` sibling emits the same
        packet sequence with ``is_continuous=0`` and per-pixel
        ``(channel_index, intensity)`` pairs (spec §4.17).

        :param image: :class:`ttio.MSImage` to emit. ``image`` must
            not be ``None`` (caller's responsibility to gate on
            :meth:`SpectralDataset.image_for_kind` is not None).
        """
        if image is None:
            raise ValueError("write_image: image must not be None")

        width = int(image.width)
        height = int(image.height)
        bins = int(image.spectral_points)
        axis = np.asarray(image.mz_axis, dtype=np.float64)
        scan_pattern_byte = _scan_pattern_to_byte(image.scan_pattern)

        # -- IMAGE_HEADER (0x13) ----------------------------------------
        # Modality=0 (MS) has no modality-specific extras (spec §4.16).
        self._emit_image_header(
            modality=0,
            width=width,
            height=height,
            bins=bins,
            pixel_size_x=float(image.pixel_size_x),
            pixel_size_y=float(image.pixel_size_y),
            scan_pattern_byte=scan_pattern_byte,
            axis_kind=0,  # mz
            axis=axis,
            is_continuous=1,
            title=image.title or "",
            isa_id=image.isa_investigation_id or "",
            extras=b"",
        )

        # -- IMAGE_PIXEL (0x14) — one per pixel -------------------------
        # Continuous-mode payload: x(u32) + y(u32) + precision(u8) +
        # compression(u8) + payload_length(u32) + intensities[..].
        # Always FLOAT64 (precision=1), uncompressed (compression=0).
        precision = 1  # FLOAT64
        compression = 0  # NONE
        payload_len = 8 * bins
        cube = np.ascontiguousarray(image.intensity, dtype=np.float64)
        # cube shape: (height, width, spectral_points). Iterate raster:
        # row-major over (y, x) — matches Java's y-outer / x-inner loop.
        pixel_index = 0
        for y in range(height):
            for x in range(width):
                intensities_bytes = cube[y, x].tobytes()
                pixel_payload = b"".join((
                    struct.pack("<I", x & 0xFFFFFFFF),
                    struct.pack("<I", y & 0xFFFFFFFF),
                    struct.pack("<B", precision & 0xFF),
                    struct.pack("<B", compression & 0xFF),
                    struct.pack("<I", payload_len & 0xFFFFFFFF),
                    intensities_bytes,
                ))
                self._emit(
                    PacketType.IMAGE_PIXEL,
                    pixel_payload,
                    au_sequence=pixel_index,
                )
                pixel_index += 1

        # -- END_OF_IMAGE (0x15) ----------------------------------------
        # pixel_count_seen: uint32 per spec §4.18.
        self._emit(
            PacketType.END_OF_IMAGE,
            struct.pack("<I", pixel_index & 0xFFFFFFFF),
        )

    def write_image_processed(self, image) -> None:
        """v0.11 Task 5.1 (Deferral 1): emit an :class:`~ttio.MSImage`
        as the packet sequence
        ``IMAGE_HEADER (0x13) -> N x IMAGE_PIXEL (0x14) -> END_OF_IMAGE
        (0x15)`` in **processed mode** (sparse), where each pixel
        carries only its nonzero ``(channel_index, intensity)`` pairs
        indexed into the shared ``mz_axis``. The dense cube is
        reconstructed by the reader.

        Wire layout per transport-spec §4.17 (LITTLE-ENDIAN). The
        IMAGE_HEADER is identical to :meth:`write_image` except for
        ``is_continuous=0``; each IMAGE_PIXEL payload is::

          x(u32) + y(u32) + precision(u8) + compression(u8)
            + payload_length(u32)
            + payload_bytes = nonzero_count(u32)
                + nonzero_count × { channel_index(u32) + intensity(f64) }

        Nonzero is defined strictly as ``v != 0.0``; NaN is preserved
        verbatim (NaN compares unequal to 0.0). The MSImage data model
        stays dense; processed mode is purely a wire optimisation for
        sparse cubes.

        This is an opt-in sibling of :meth:`write_image`. Callers pick
        continuous vs processed mode explicitly today; an automatic
        heuristic (emit whichever is smaller) is a follow-up. Java
        parity: :meth:`TransportWriter.writeImageProcessed` (commit
        ``1889343e``).

        :param image: :class:`ttio.MSImage` to emit. ``image`` must
            not be ``None``.
        """
        if image is None:
            raise ValueError("write_image_processed: image must not be None")

        width = int(image.width)
        height = int(image.height)
        bins = int(image.spectral_points)
        axis = np.asarray(image.mz_axis, dtype=np.float64)
        scan_pattern_byte = _scan_pattern_to_byte(image.scan_pattern)

        # -- IMAGE_HEADER (0x13) — shared layout, is_continuous=0 -------
        # Modality=0 (MS) has no modality-specific extras (spec §4.16).
        self._emit_image_header(
            modality=0,
            width=width,
            height=height,
            bins=bins,
            pixel_size_x=float(image.pixel_size_x),
            pixel_size_y=float(image.pixel_size_y),
            scan_pattern_byte=scan_pattern_byte,
            axis_kind=0,  # mz
            axis=axis,
            is_continuous=0,  # processed
            title=image.title or "",
            isa_id=image.isa_investigation_id or "",
            extras=b"",
        )

        # -- IMAGE_PIXEL (0x14) — sparse per spec §4.17 -----------------
        # Always FLOAT64 (precision=1) uncompressed (compression=0) —
        # mirrors write_image above so the wire round-trip stays
        # byte-exact with the cube.
        precision = 1  # FLOAT64
        compression = 0  # NONE
        cube = np.ascontiguousarray(image.intensity, dtype=np.float64)
        pixel_index = 0
        for y in range(height):
            for x in range(width):
                spec = cube[y, x]
                # Nonzero mask (v != 0.0); NaN preserved verbatim.
                nz_mask = spec != 0.0
                nz_indices = np.nonzero(nz_mask)[0]
                nz_count = int(nz_indices.size)
                payload_len = 4 + nz_count * (4 + 8)
                # Pack sparse entries: nonzero_count(u32) + repeating
                # (channel_index u32 + intensity f64). Match Java's
                # ascending-channel iteration order (numpy.nonzero
                # returns sorted indices for rank-1 input).
                if nz_count > 0:
                    fmt = "<" + ("Id" * nz_count)
                    pairs: list = []
                    for ch in nz_indices:
                        pairs.append(int(ch))
                        pairs.append(float(spec[int(ch)]))
                    entries = struct.pack(fmt, *pairs)
                else:
                    entries = b""
                pixel_payload = b"".join((
                    struct.pack("<I", x & 0xFFFFFFFF),
                    struct.pack("<I", y & 0xFFFFFFFF),
                    struct.pack("<B", precision & 0xFF),
                    struct.pack("<B", compression & 0xFF),
                    struct.pack("<I", payload_len & 0xFFFFFFFF),
                    struct.pack("<I", nz_count & 0xFFFFFFFF),
                    entries,
                ))
                self._emit(
                    PacketType.IMAGE_PIXEL,
                    pixel_payload,
                    au_sequence=pixel_index,
                )
                pixel_index += 1

        # -- END_OF_IMAGE (0x15) ----------------------------------------
        self._emit(
            PacketType.END_OF_IMAGE,
            struct.pack("<I", pixel_index & 0xFFFFFFFF),
        )

    def write_raman_image(self, image) -> None:
        """v0.11 Task 5.3 (Deferral 1): emit a :class:`~ttio.RamanImage`
        as the packet sequence
        ``IMAGE_HEADER (0x13) -> N x IMAGE_PIXEL (0x14) -> END_OF_IMAGE
        (0x15)`` with ``modality=1``.

        Wire layout per transport-spec §4.16. The shared axis on the
        IMAGE_HEADER carries the Raman wavenumbers vector
        (``axis_kind = 1 = wavenumber``). The ``modality_extras`` slot
        at the tail of the IMAGE_HEADER carries the Raman-specific
        fields::

            excitation_wavelength_nm:  float64
            laser_power_mw:            float64

        (16 bytes total.)

        Each pixel rides as a continuous-mode IMAGE_PIXEL whose
        ``payload_bytes`` is a dense vector of ``spectrum_bins``
        FLOAT64 intensities at the shared wavenumber axis. Java
        parity: :meth:`TransportWriter.writeRamanImage` (commit
        ``f99ec47d``).
        """
        if image is None:
            raise ValueError("write_raman_image: image must not be None")

        width = int(image.width)
        height = int(image.height)
        bins = int(image.spectral_points)
        axis = np.asarray(image.wavenumbers, dtype=np.float64)
        scan_pattern_byte = _scan_pattern_to_byte(image.scan_pattern)
        # Raman modality_extras: 8B excitation_wavelength_nm + 8B laser_power_mw.
        extras = struct.pack(
            "<dd",
            float(image.excitation_wavelength_nm),
            float(image.laser_power_mw),
        )

        self._emit_image_header(
            modality=1,
            width=width,
            height=height,
            bins=bins,
            pixel_size_x=float(image.pixel_size_x),
            pixel_size_y=float(image.pixel_size_y),
            scan_pattern_byte=scan_pattern_byte,
            axis_kind=1,  # wavenumber
            axis=axis,
            is_continuous=1,
            title=image.title or "",
            isa_id=image.isa_investigation_id or "",
            extras=extras,
        )

        precision = 1  # FLOAT64
        compression = 0  # NONE
        payload_len = 8 * bins
        cube = np.ascontiguousarray(image.intensity, dtype=np.float64)
        pixel_index = 0
        for y in range(height):
            for x in range(width):
                intensities_bytes = cube[y, x].tobytes()
                pixel_payload = b"".join((
                    struct.pack("<I", x & 0xFFFFFFFF),
                    struct.pack("<I", y & 0xFFFFFFFF),
                    struct.pack("<B", precision & 0xFF),
                    struct.pack("<B", compression & 0xFF),
                    struct.pack("<I", payload_len & 0xFFFFFFFF),
                    intensities_bytes,
                ))
                self._emit(
                    PacketType.IMAGE_PIXEL,
                    pixel_payload,
                    au_sequence=pixel_index,
                )
                pixel_index += 1

        self._emit(
            PacketType.END_OF_IMAGE,
            struct.pack("<I", pixel_index & 0xFFFFFFFF),
        )

    def write_ir_image(self, image) -> None:
        """v0.11 Task 5.3 (Deferral 1): emit an :class:`~ttio.IRImage`
        as the packet sequence
        ``IMAGE_HEADER (0x13) -> N x IMAGE_PIXEL (0x14) -> END_OF_IMAGE
        (0x15)`` with ``modality=2``.

        Wire layout per transport-spec §4.16. The shared axis on the
        IMAGE_HEADER carries the IR wavenumbers vector
        (``axis_kind = 1 = wavenumber``). The ``modality_extras`` slot
        at the tail of the IMAGE_HEADER carries the IR-specific
        fields::

            ir_mode:            uint8   # 0=transmittance, 1=absorbance
            resolution_cm_inv:  float64

        (9 bytes total.) Java parity:
        :meth:`TransportWriter.writeIRImage` (commit ``f99ec47d``).
        """
        if image is None:
            raise ValueError("write_ir_image: image must not be None")

        from ..enums import IRMode

        width = int(image.width)
        height = int(image.height)
        bins = int(image.spectral_points)
        axis = np.asarray(image.wavenumbers, dtype=np.float64)
        scan_pattern_byte = _scan_pattern_to_byte(image.scan_pattern)
        # IR modality_extras: u8 ir_mode + f64 resolution_cm_inv = 9B.
        ir_mode_byte = 1 if image.mode == IRMode.ABSORBANCE else 0
        extras = struct.pack(
            "<Bd",
            ir_mode_byte & 0xFF,
            float(image.resolution_cm_inv),
        )

        self._emit_image_header(
            modality=2,
            width=width,
            height=height,
            bins=bins,
            pixel_size_x=float(image.pixel_size_x),
            pixel_size_y=float(image.pixel_size_y),
            scan_pattern_byte=scan_pattern_byte,
            axis_kind=1,  # wavenumber
            axis=axis,
            is_continuous=1,
            title=image.title or "",
            isa_id=image.isa_investigation_id or "",
            extras=extras,
        )

        precision = 1  # FLOAT64
        compression = 0  # NONE
        payload_len = 8 * bins
        cube = np.ascontiguousarray(image.intensity, dtype=np.float64)
        pixel_index = 0
        for y in range(height):
            for x in range(width):
                intensities_bytes = cube[y, x].tobytes()
                pixel_payload = b"".join((
                    struct.pack("<I", x & 0xFFFFFFFF),
                    struct.pack("<I", y & 0xFFFFFFFF),
                    struct.pack("<B", precision & 0xFF),
                    struct.pack("<B", compression & 0xFF),
                    struct.pack("<I", payload_len & 0xFFFFFFFF),
                    intensities_bytes,
                ))
                self._emit(
                    PacketType.IMAGE_PIXEL,
                    pixel_payload,
                    au_sequence=pixel_index,
                )
                pixel_index += 1

        self._emit(
            PacketType.END_OF_IMAGE,
            struct.pack("<I", pixel_index & 0xFFFFFFFF),
        )

    # ----------------------------------------------- v0.11 §4.19 / §4.20

    def write_identifications_table(self, rows) -> None:
        """v0.11 Task 2.7: emit an
        :data:`PacketType.IDENTIFICATIONS_TABLE` (0x16) packet carrying
        the full identifications table as a single length-prefixed
        Apache Arrow IPC stream. Wire layout per transport-spec §4.19::

            arrow_ipc_length:    uint32
            arrow_ipc:           bytes[arrow_ipc_length]   # self-describing IPC

        All multi-byte integers LITTLE-ENDIAN per spec §1.7. The Arrow
        IPC stream carries its own schema, row count, and null bitmaps,
        so no per-row TLV envelope is needed. Empty lists are a no-op
        (no packet emitted) per spec §5.4 step 6 ("zero or more").

        Java parity: :meth:`TransportWriter.writeIdentifications`
        (commit ``a6faab16``).

        :param rows: iterable of
            :class:`ttio.identification.Identification`. ``None``
            raises ``ValueError``.
        """
        if rows is None:
            raise ValueError(
                "write_identifications_table: rows must not be None"
            )
        rows = list(rows)
        if not rows:
            return
        from .arrow_ipc import encode_identifications
        ipc = encode_identifications(rows)
        payload = struct.pack("<I", len(ipc) & 0xFFFFFFFF) + ipc
        self._emit(PacketType.IDENTIFICATIONS_TABLE, payload)

    def write_quantifications_table(self, rows) -> None:
        """v0.11 Task 2.7: emit a
        :data:`PacketType.QUANTIFICATIONS_TABLE` (0x17) packet carrying
        the full quantifications table as a single length-prefixed
        Apache Arrow IPC stream. Wire layout per transport-spec §4.20 —
        identical shape to §4.19 but with a distinct packet type so
        receivers can dispatch without parsing the IPC payload first.

        All multi-byte integers LITTLE-ENDIAN per spec §1.7. Empty
        lists are a no-op (spec §5.4 step 6).

        Java parity: :meth:`TransportWriter.writeQuantifications`
        (commit ``a6faab16``).

        :param rows: iterable of
            :class:`ttio.quantification.Quantification`. ``None``
            raises ``ValueError``.
        """
        if rows is None:
            raise ValueError(
                "write_quantifications_table: rows must not be None"
            )
        rows = list(rows)
        if not rows:
            return
        from .arrow_ipc import encode_quantifications
        ipc = encode_quantifications(rows)
        payload = struct.pack("<I", len(ipc) & 0xFFFFFFFF) + ipc
        self._emit(PacketType.QUANTIFICATIONS_TABLE, payload)

    # ----------------------------------------------- v0.11 §4.22 (Stage 6)

    def write_subject_metadata(self, rows) -> None:
        """Stage 6 / Task 6.3: emit a
        :data:`PacketType.SUBJECT_METADATA` (0x19) packet carrying the
        full Subject table as a single length-prefixed Apache Arrow
        IPC stream. Wire layout per transport-spec §4.22::

            arrow_ipc_length:    uint32
            arrow_ipc:           bytes[arrow_ipc_length]   # self-describing IPC

        All multi-byte integers LITTLE-ENDIAN per spec §1.7. The Arrow
        IPC stream carries its own schema, row count, and null bitmaps,
        so no per-row TLV envelope is needed. Empty lists are a no-op
        (spec §5.4 step 5: "zero or more").

        Java parity: :meth:`TransportWriter.writeSubjectMetadata`
        (commit ``dd211600``).

        :param rows: iterable of :class:`ttio.subject.Subject`.
            ``None`` raises ``ValueError``.
        """
        if rows is None:
            raise ValueError(
                "write_subject_metadata: rows must not be None"
            )
        rows = list(rows)
        if not rows:
            return
        from .arrow_ipc import encode_subjects
        ipc = encode_subjects(rows)
        payload = struct.pack("<I", len(ipc) & 0xFFFFFFFF) + ipc
        self._emit(PacketType.SUBJECT_METADATA, payload)

    def write_sample_metadata(self, rows) -> None:
        """Stage 6 / Task 6.3: emit a
        :data:`PacketType.SAMPLE_METADATA` (0x1A) packet carrying the
        full Sample table as a single length-prefixed Apache Arrow IPC
        stream. Wire layout per transport-spec §4.22 — identical shape
        to the SUBJECT_METADATA framing with a distinct packet type so
        receivers can dispatch without parsing the IPC payload first.

        All multi-byte integers LITTLE-ENDIAN per spec §1.7. Empty
        lists are a no-op (spec §5.4 step 5).

        Java parity: :meth:`TransportWriter.writeSampleMetadata`
        (commit ``dd211600``).

        :param rows: iterable of :class:`ttio.sample.Sample`. ``None``
            raises ``ValueError``.
        """
        if rows is None:
            raise ValueError(
                "write_sample_metadata: rows must not be None"
            )
        rows = list(rows)
        if not rows:
            return
        from .arrow_ipc import encode_samples
        ipc = encode_samples(rows)
        payload = struct.pack("<I", len(ipc) & 0xFFFFFFFF) + ipc
        self._emit(PacketType.SAMPLE_METADATA, payload)

    def write_end_of_dataset(
        self, *, dataset_id: int, final_au_sequence: int
    ) -> None:
        """Emit a :data:`PacketType.END_OF_DATASET` sentinel packet.

        Terminates a dataset's access-unit run. The reader uses
        :paramref:`final_au_sequence` to verify it observed every
        expected AU.

        Parameters
        ----------
        dataset_id : int
            Owning dataset id.
        final_au_sequence : int
            One past the last ``au_sequence`` emitted for the
            dataset (i.e. the dataset's AU count).
        """
        payload = struct.pack(
            "<HI", dataset_id & 0xFFFF, final_au_sequence & 0xFFFFFFFF
        )
        self._emit(PacketType.END_OF_DATASET, payload, dataset_id=dataset_id)

    def write_end_of_stream(self) -> None:
        """Emit a :data:`PacketType.END_OF_STREAM` sentinel packet.

        Marks the end of the transport stream. The reader stops
        ingesting after this packet.
        """
        self._emit(PacketType.END_OF_STREAM, b"")

    def write_dataset(self, dataset: SpectralDataset) -> None:
        """Walk ``dataset`` and emit the complete transport packet sequence.

        The emission order follows the v0.11 prelude rules:
        ``StreamHeader`` then the optional v0.11 prelude
        (encryption algorithm, dataset provenance, subjects,
        samples, reference groups, image cubes, identifications,
        quantifications) then per-dataset
        ``DatasetHeader`` / access-unit run /
        ``EndOfDataset``, finishing with ``EndOfStream``.
        Spectral runs are emitted with dataset ids ``1..N`` and
        genomic runs with ``N+1..N+M``.

        Parameters
        ----------
        dataset : SpectralDataset
            Source container to serialise. Borrowed; not closed.
        """
        runs = list(dataset.all_runs.items())
        genomic_runs = list(getattr(dataset, "genomic_runs", {}).items())
        features = list(dataset.feature_flags.features)
        # Phase 2c-T: declare bulk-mode in the StreamHeader features
        # list when the writer is bulk-enabled AND there is at least
        # one genomic run that has v2 blobs to carry. Receivers see
        # this flag and dispatch to the bulk path.
        bulk_active = (
            self._use_bulk_mode
            and any(_bulk_carriable(g) for _, g in genomic_runs)
            and BULK_MODE_V2_BLOBS_FEATURE not in features
        )
        if bulk_active:
            features = features + [BULK_MODE_V2_BLOBS_FEATURE]

        # v0.11 Task 2.4/2.5/2.6/2.7: detect v0.11 content (references +
        # encryption algorithm + dataset provenance + image cube +
        # identifications/quantifications tables today; subjects and
        # samples land at the same prelude insertion point per §5.4
        # ordering in a follow-up task). Java parity:
        # TransportWriter.writeDataset (commits 530a5833 + 563e09c3
        # + a6b1e5d9 + a6faab16 + dc0de926).
        refs = getattr(dataset, "references", {}) or {}
        # ``provenance`` is a method on SpectralDataset (not a
        # property) — call it eagerly so the empty-vs-nonempty check
        # is cheap and the actual records are reused below.
        try:
            dataset_provenance = list(dataset.provenance())
        except Exception:  # pragma: no cover - defensive
            dataset_provenance = []
        # Images are read via image_for_kind() — read once so the cache
        # hit is shared between the flag-detect and the emit branch.
        from ..enums import ImageKind
        dataset_image = dataset.image_for_kind(ImageKind.MS)
        dataset_raman_image = dataset.image_for_kind(ImageKind.RAMAN)
        dataset_ir_image = dataset.image_for_kind(ImageKind.IR)
        # identifications + quantifications are methods on
        # SpectralDataset (matching ``provenance()``). Empty lists do
        # not emit a packet (spec §5.4 step 6 says "zero or more")
        # and do not trigger the v0.11 feature flag.
        try:
            dataset_identifications = list(dataset.identifications())
        except Exception:  # pragma: no cover - defensive
            dataset_identifications = []
        try:
            dataset_quantifications = list(dataset.quantifications())
        except Exception:  # pragma: no cover - defensive
            dataset_quantifications = []
        # Stage 6 / Task 6.3 (transport-spec v0.11 §4.22): subjects +
        # samples are lazy properties on SpectralDataset (not methods).
        # Read once so the cache hit is shared between the flag-detect
        # and the emit branch. Defensive try/except follows the same
        # pattern as the identifications/quantifications block above.
        try:
            dataset_subjects = list(getattr(dataset, "subjects", []) or [])
        except Exception:  # pragma: no cover - defensive
            dataset_subjects = []
        try:
            dataset_samples = list(getattr(dataset, "samples", []) or [])
        except Exception:  # pragma: no cover - defensive
            dataset_samples = []
        has_v011_content = (
            len(refs) > 0
            or bool(getattr(dataset, "is_encrypted", False))
            or len(dataset_provenance) > 0
            or dataset_image is not None
            or dataset_raman_image is not None
            or dataset_ir_image is not None
            or len(dataset_identifications) > 0
            or len(dataset_quantifications) > 0
            or len(dataset_subjects) > 0
            or len(dataset_samples) > 0
        )
        if has_v011_content and TRANSPORT_V0_11_FEATURE not in features:
            features = features + [TRANSPORT_V0_11_FEATURE]

        self.write_stream_header(
            format_version="1.2",
            title=dataset.title or "",
            isa_investigation=dataset.isa_investigation_id or "",
            features=features,
            n_datasets=len(runs) + len(genomic_runs),
        )

        # v0.11 Task 2.4/2.5/2.6/2.7/6.3: v0.11 prelude — per §5.4
        # ordering, v0.11 sections come BEFORE the v0.10 dataset/run
        # sections, and the sub-sections appear in this order:
        #   §5.4.1 ENCRYPTION_ALGORITHM
        #   §5.4.2 DATASET_PROVENANCE
        #   §5.4.3 SUBJECT_METADATA / SAMPLE_METADATA  (subjects first)
        #   §5.4.4 reference groups
        #   §5.4.5 image cubes
        #   §5.4.6 IDENTIFICATIONS_TABLE / QUANTIFICATIONS_TABLE
        if has_v011_content:
            if getattr(dataset, "is_encrypted", False):
                algo = getattr(dataset, "encrypted_algorithm", "") or ""
                if algo:
                    self.write_encryption_algorithm(algo)
            if dataset_provenance:
                self.write_dataset_provenance(dataset_provenance)
            # §5.4.3 Stage 6 / Task 6.3: SUBJECT_METADATA (0x19) emits
            # before SAMPLE_METADATA (0x1A) so a downstream reader sees
            # subjects ahead of any samples that soft-FK into them.
            # Empty lists emit NO packet (spec §5.4 step 5: "zero or
            # more"). Java parity: TransportWriter.writeDataset (commit
            # dd211600).
            if dataset_subjects:
                self.write_subject_metadata(dataset_subjects)
            if dataset_samples:
                self.write_sample_metadata(dataset_samples)
            # §5.4.4: reference groups, one packet sequence per
            # ReferenceImport. Empty refs dict emits nothing.
            for ref in refs.values():
                self.write_reference_group(ref)
            # §5.4.5 image cubes: MS → Raman → IR (deterministic
            # emission order when more than one modality is populated
            # on the same dataset). Java parity: TransportWriter.write-
            # Dataset (commit f99ec47d).
            if dataset_image is not None:
                self.write_image(dataset_image)
            if dataset_raman_image is not None:
                self.write_raman_image(dataset_raman_image)
            if dataset_ir_image is not None:
                self.write_ir_image(dataset_ir_image)
            # §5.4.6: identifications first, then quantifications.
            # Each empty list is a no-op (spec §5.4 step 6).
            if dataset_identifications:
                self.write_identifications_table(dataset_identifications)
            if dataset_quantifications:
                self.write_quantifications_table(dataset_quantifications)
        for i, (name, run) in enumerate(runs, start=1):
            self.write_dataset_header(
                dataset_id=i,
                name=name,
                acquisition_mode=int(run.acquisition_mode),
                spectrum_class=run.spectrum_class,
                channel_names=list(run.channel_names),
                instrument_json=_instrument_config_json(run),
                expected_au_count=len(run),
            )
        for j, (name, grun) in enumerate(genomic_runs, start=len(runs) + 1):
            self.write_dataset_header(
                dataset_id=j,
                name=name,
                acquisition_mode=int(grun.acquisition_mode),
                spectrum_class="TTIOGenomicRead",
                channel_names=genomic_channel_names(grun),
                instrument_json=_genomic_run_metadata_json(grun),
                expected_au_count=len(grun),
            )
        for i, (name, run) in enumerate(runs, start=1):
            self._emit_run_access_units(dataset_id=i, run=run)
            self.write_end_of_dataset(dataset_id=i, final_au_sequence=len(run))
        for j, (name, grun) in enumerate(genomic_runs, start=len(runs) + 1):
            if self._use_bulk_mode:
                self._emit_genomic_run_v2_blobs(dataset_id=j, run=grun)
            self._emit_genomic_run_access_units(dataset_id=j, run=grun)
            self.write_end_of_dataset(dataset_id=j, final_au_sequence=len(grun))
        self.write_end_of_stream()

    def write_genomic_run(
        self, *, dataset_id: int, name: str, run
    ) -> None:
        """Write a single GenomicRun as a stream segment.

        Used by callers that drive emission manually (multiplexed
        streams, M89.4). The dataset header + AUs + end-of-dataset
        are emitted; the caller is responsible for stream framing.
        """
        self.write_dataset_header(
            dataset_id=dataset_id,
            name=name,
            acquisition_mode=int(run.acquisition_mode),
            spectrum_class="TTIOGenomicRead",
            channel_names=["sequences", "qualities"],
            instrument_json=_genomic_run_metadata_json(run),
            expected_au_count=len(run),
        )
        self._emit_genomic_run_access_units(dataset_id=dataset_id, run=run)
        self.write_end_of_dataset(
            dataset_id=dataset_id, final_au_sequence=len(run)
        )

    def _emit_genomic_run_v2_blobs(self, *, dataset_id: int, run) -> None:
        """Phase 2c-T: probe ``run`` for v2 codec blobs and emit the
        matching ``BlobV2*`` packets.

        Each blob is independently optional — a run with no
        ``mate_info/inline_v2`` group (rare; happens when
        ``len(mate_chromosomes) == 0`` at write time) emits no
        ``BLOB_V2_MATE_INFO`` packet, and so on. The receiver fills
        the missing channels from the per-AU stream just as in
        per-AU mode.

        A blocks_v1 run (format-spec 10.12) with more than one
        block has no single per-channel blob to carry verbatim; such
        runs are sent per-AU. A one-block run's blobs are exactly the
        whole-channel blobs and are carried as before.
        """
        if not _bulk_carriable(run):
            return
        sig = run.group.open_group("signal_channels")

        # mate_info/inline_v2 + chrom_names table
        try:
            mate_grp = sig.open_group("mate_info")
            try:
                mate_ds = mate_grp.open_dataset("inline_v2")
                mate_blob = bytes(
                    mate_ds.read(offset=0, count=int(mate_ds.length))
                )
                # chrom_names compound dataset → list of strings.
                chrom_names = _read_mate_chrom_names_table(mate_grp)
                self.write_blob_v2_mate_info(
                    dataset_id=dataset_id,
                    chrom_names=chrom_names,
                    blob=mate_blob,
                )
            except KeyError:
                pass
        except KeyError:
            pass

        # read_names — flat uint8 dataset with @compression=15.
        try:
            rn_ds = sig.open_dataset("read_names")
            codec_id = io_attr_int(rn_ds, "compression", default=0) or 0
            if codec_id == CODEC_ID_NAME_TOKENIZED_V2 and int(rn_ds.length) > 0:
                rn_blob = bytes(
                    rn_ds.read(offset=0, count=int(rn_ds.length))
                )
                self.write_blob_v2_name_tok(
                    dataset_id=dataset_id, blob=rn_blob,
                )
        except KeyError:
            pass

        # sequences — group with refdiff_v2 child dataset (v1.8 layout).
        try:
            seq_grp = sig.open_group("sequences")
            try:
                rd_ds = seq_grp.open_dataset("refdiff_v2")
                rd_blob = bytes(
                    rd_ds.read(offset=0, count=int(rd_ds.length))
                )
                self.write_blob_v2_ref_diff(
                    dataset_id=dataset_id,
                    reference_uri=run.reference_uri or "",
                    blob=rd_blob,
                )
            except KeyError:
                pass
        except KeyError:
            # sequences is a flat dataset (no ref-diff) — nothing to ship.
            pass

    def _emit_genomic_run_access_units(self, *, dataset_id: int, run) -> None:
        """Emit one ACCESS_UNIT packet per AlignedRead in ``run``.

        M89.2: per-read fixed fields go into the AU's genomic suffix
        (chromosome / position / mapping_quality / flags). The
        variable-length sequences and qualities arrays ride as two
        UINT8 channels with the per-read slice as data.

        M90.9: compound fields now also round-trip on the wire.
        cigar, read_name, mate_chromosome ride as additional UINT8
        string channels (one per AU). mate_position + template_length
        live in the M90.9 mate extension at the end of the AU genomic
        suffix.

        M90.10: when the source channel carries an ``@compression``
        attribute naming an M86 codec (RANS_ORDER0/1, BASE_PACK), the
        writer re-encodes each per-AU slice with the same codec on
        the wire. The wire ChannelData.compression byte tells the
        reader which decoder to dispatch.

        Body now delegates to :func:`_iter_genomic_run_access_units`
        so the walker (#141) can yield the same AUs without
        re-implementing the construction.
        """
        for i, au in _iter_genomic_run_access_units(run):
            self.write_access_unit(
                dataset_id=dataset_id, au_sequence=i, au=au
            )

    def _emit_run_access_units(
        self, *, dataset_id: int, run: AcquisitionRun
    ) -> None:
        """Hot path: emit AccessUnit packets for every spectrum in ``run``.

        Bulk-reads each channel dataset once up-front and slices per AU.
        Skips per-spectrum ``_materialize_spectrum`` (which was ~60% of
        encode walltime through h5py hyperslab reads) and dataclass
        constructions.
        """
        # Pre-compute everything stable across spectra once per run.
        channel_names = list(run.channel_names)
        channel_name_prefixes = [
            _CHANNEL_NAMELEN_STRUCT.pack(len(nb)) + nb
            for nb in (cn.encode("utf-8") for cn in channel_names)
        ]
        wire_class = _SPECTRUM_CLASS_TO_WIRE.get(run.spectrum_class, 0) & 0xFF
        acq_mode = int(run.acquisition_mode) & 0xFF
        is_ms_class = run.spectrum_class == "TTIOMassSpectrum"
        is_pixel_class = wire_class == 4
        use_compression = self._use_compression
        if use_compression:
            compression_enum, compress = _wire_channel_encoder(self._compression_codec)
            compression_enum &= 0xFF
        else:
            compression_enum = int(Compression.NONE) & 0xFF
            compress = None
        precision_enum = int(Precision.FLOAT64) & 0xFF
        unknown_polarity_wire = _POLARITY_TO_WIRE[Polarity.UNKNOWN]

        # Channel data is read in windows of WINDOW spectra through
        # AcquisitionRun.channel_range (per-block codec 17 decode,
        # hyperslab reads for plain channels), so a whole channel is
        # never held at once. Channels missing from the run are skipped.
        index = run.index
        present: list[int] = []
        for ci, cname in enumerate(channel_names):
            try:
                run._signal_dataset(cname)
            except KeyError:
                if run._numpress_channels.get(cname) is None and cname not in run._decrypted_channels:
                    continue
            present.append(ci)
        WINDOW = 4096

        # Index columns are numpy arrays already — slice per-i.
        offsets = index.offsets
        lengths = index.lengths
        rts = index.retention_times
        pmzs = index.precursor_mzs
        pcs = index.precursor_charges
        bpis = index.base_peak_intensities
        ms_levels = index.ms_levels
        polarities_wire = (
            np.array(
                [_POLARITY_TO_WIRE.get(Polarity(int(p)), 2) for p in index.polarities],
                dtype="<i4",
            )
            if is_ms_class
            else None
        )

        # Hoist method lookups out of the loop.
        stream_write = self._stream.write
        header_pack = _HEADER_STRUCT.pack
        au_prefix_pack = _AU_PREFIX_STRUCT.pack
        channel_suffix_pack = _CHANNEL_SUFFIX_STRUCT.pack
        pixel_pack = _AU_PIXEL_STRUCT.pack
        crc32c_ = crc32c
        checksum_pack = _CHECKSUM_STRUCT.pack
        use_checksum = self._use_checksum
        flags = int(PacketFlag.HAS_CHECKSUM) if use_checksum else 0
        ac_type = int(PacketType.ACCESS_UNIT) & 0xFF
        now_ns_ = now_ns
        did = dataset_id & 0xFFFF

        n_spectra = len(run)
        window_base = 0
        window_end = 0
        window_arrays: list[tuple[int, np.ndarray]] = []
        for j in range(n_spectra):
            start = int(offsets[j])
            length = int(lengths[j])
            stop = start + length
            if j >= window_end:
                window_end = min(n_spectra, j + WINDOW)
                window_base = start
                total = int(offsets[window_end - 1]) + int(lengths[window_end - 1]) - window_base
                window_arrays = [
                    (ci, np.ascontiguousarray(
                        run.channel_range(channel_names[ci], window_base, total), dtype="<f8"))
                    for ci in present]

            # Channel data collection from the current window.
            channel_chunks: list[bytes] = []
            n_channels = 0
            for ci, win_arr in window_arrays:
                seg = win_arr[start - window_base:stop - window_base]
                payload_bytes = compress(seg) if use_compression else seg.tobytes()
                channel_chunks.append(channel_name_prefixes[ci])
                channel_chunks.append(channel_suffix_pack(
                    precision_enum,
                    compression_enum,
                    length & 0xFFFFFFFF,
                    len(payload_bytes) & 0xFFFFFFFF,
                ))
                channel_chunks.append(payload_bytes)
                n_channels += 1

            if is_ms_class:
                ms_level = int(ms_levels[j])
                polarity_wire = int(polarities_wire[j]) if polarities_wire is not None else unknown_polarity_wire
            else:
                ms_level = 0
                polarity_wire = unknown_polarity_wire

            au_prefix = au_prefix_pack(
                wire_class,
                acq_mode,
                ms_level & 0xFF,
                polarity_wire & 0xFF,
                float(rts[j]),
                float(pmzs[j]),
                int(pcs[j]) & 0xFF,
                0.0,
                float(bpis[j]),
                n_channels & 0xFF,
            )
            payload_parts = [au_prefix, *channel_chunks]
            if is_pixel_class:
                payload_parts.append(pixel_pack(0, 0, 0))
            payload = b"".join(payload_parts)

            header = header_pack(
                HEADER_MAGIC,
                VERSION,
                ac_type,
                flags & 0xFFFF,
                did,
                j & 0xFFFFFFFF,
                len(payload) & 0xFFFFFFFF,
                now_ns_() & 0xFFFFFFFFFFFFFFFF,
            )
            if use_checksum:
                stream_write(header + payload + checksum_pack(crc32c_(payload)))
            else:
                stream_write(header + payload)


def _scan_pattern_to_byte(scan_pattern: str | None) -> int:
    """v0.11 Task 2.6: map an :class:`~ttio.MSImage` scan-pattern
    string to the wire byte per transport-spec §4.16
    (``0=raster/flyback, 1=meander, 2=random``).

    The on-disk format uses ``"raster"`` as the default name for the
    flyback pattern. Unknown values map to 0 (raster) defensively.
    Java parity: :meth:`TransportWriter.scanPatternToByte`.
    """
    if not scan_pattern:
        return 0
    table = {"raster": 0, "flyback": 0, "meander": 1, "random": 2}
    return table.get(scan_pattern, 0)


def _provenance_params_json(params) -> str:
    """v0.11 Task 2.5: serialise the ``parameters`` dict of a
    :class:`~ttio.provenance.ProvenanceRecord` to the canonical wire
    JSON form.

    Format: ``{"k":"v","k2":"v2"}`` with keys sort_keys-ordered and
    no whitespace between separators. Empty dict renders as ``"{}"``.
    Java parity: :meth:`ProvenanceRecord.parametersJson` produces the
    same braces / quote / no-whitespace shape; Python sorts keys for
    a deterministic on-wire ordering (Java preserves insertion order
    via ``Map.copyOf``, which for the small Maps used in practice is
    LinkedHashMap-equivalent).
    """
    if not params:
        return "{}"
    # Coerce values to str so a dict[str, Any] (per the Python
    # ProvenanceRecord type annotation) still produces the
    # string-valued shape Java emits. Sort keys for stability.
    coerced = {str(k): str(v) for k, v in params.items()}
    return json.dumps(coerced, sort_keys=True, separators=(",", ":"))


def _provenance_csv_join(refs) -> str:
    """v0.11 Task 2.5: comma-join a list of refs into the wire form.

    Empty list → empty string. No quoting/escaping is performed — per
    spec §4.21, refs are URIs that have been URL-encoded so they
    cannot themselves contain commas. Java parity:
    :meth:`TransportWriter.csvJoin`.
    """
    if not refs:
        return ""
    return ",".join(str(r) for r in refs)


def _provenance_csv_split(csv: str) -> list[str]:
    """Reverse of :func:`_provenance_csv_join`. Empty string → empty
    list. Java parity: :meth:`TransportReader.parseCsv`."""
    if not csv:
        return []
    # No quoting / escaping in v0.11 — URIs are URL-encoded so they
    # cannot themselves contain commas. Plain split mirrors Java's
    # ``csv.split(",", -1)``.
    return csv.split(",")


def _provenance_params_parse(json_blob: str) -> dict:
    """Reverse of :func:`_provenance_params_json`. Empty / ``"{}"``
    blob → empty dict. Tolerant of parse failure (returns empty dict).
    Java parity: :meth:`TransportReader.parseParametersJson`."""
    if not json_blob or json_blob == "{}":
        return {}
    try:
        parsed = json.loads(json_blob)
        if isinstance(parsed, dict):
            return {str(k): str(v) for k, v in parsed.items()}
    except (json.JSONDecodeError, TypeError):
        pass
    return {}


def _instrument_config_json(run: AcquisitionRun) -> str:
    cfg = run.instrument_config
    return json.dumps({
        "manufacturer": cfg.manufacturer,
        "model": cfg.model,
        "serial_number": cfg.serial_number,
        "source_type": cfg.source_type,
        "analyzer_type": cfg.analyzer_type,
        "detector_type": cfg.detector_type,
    }, sort_keys=True)


def _genomic_run_metadata_json(run) -> str:
    """Serialise per-genomic-run metadata for the dataset header.

    Reuses the instrument_json slot in the DATASET_HEADER packet —
    GenomicRun has its own metadata fields (reference_uri, platform,
    sample_name) instead of an InstrumentConfig. M89.2.
    """
    return json.dumps({
        "reference_uri": getattr(run, "reference_uri", "") or "",
        "platform": getattr(run, "platform", "") or "",
        "sample_name": getattr(run, "sample_name", "") or "",
        "modality": getattr(run, "modality", "") or "",
        # M97: always present, "" when the run has no read role, so the
        # JSON stays byte-identical across the 3 SDK writers.
        "read_role": getattr(run, "read_role", "") or "",
    }, sort_keys=True)


def _spectrum_to_access_unit(
    spectrum: Spectrum,
    run: AcquisitionRun,
    *,
    use_compression: bool = False,
    compression_codec: str = "float_delta_zstd",
) -> AccessUnit:
    wire_class = _SPECTRUM_CLASS_TO_WIRE.get(run.spectrum_class, 0)
    ms_level = 0
    polarity_wire = _POLARITY_TO_WIRE[Polarity.UNKNOWN]
    if isinstance(spectrum, MassSpectrum):
        ms_level = spectrum.ms_level
        polarity_wire = _POLARITY_TO_WIRE.get(spectrum.polarity, 2)

    bpi = float(run.index.base_peak_intensity_at(spectrum.index_position))

    channels: list[ChannelData] = []
    for cname in run.channel_names:
        if not spectrum.has_signal_array(cname):
            continue
        sa = spectrum.signal_array(cname)
        arr = np.asarray(sa.data).astype("<f8", copy=False)
        if use_compression:
            compression, encode = _wire_channel_encoder(compression_codec)
            payload = encode(arr)
        else:
            payload = arr.tobytes()
            compression = int(Compression.NONE)
        channels.append(ChannelData(
            name=cname,
            precision=int(Precision.FLOAT64),
            compression=compression,
            n_elements=int(arr.size),
            data=payload,
        ))

    return AccessUnit(
        spectrum_class=wire_class,
        acquisition_mode=int(run.acquisition_mode),
        ms_level=ms_level,
        polarity=polarity_wire,
        retention_time=float(spectrum.scan_time_seconds),
        precursor_mz=float(spectrum.precursor_mz),
        precursor_charge=int(spectrum.precursor_charge),
        ion_mobility=0.0,
        base_peak_intensity=bpi,
        channels=channels,
    )
