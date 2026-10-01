"""Transport-stream reader: :class:`TransportReader` + reader helpers.

Pure code-movement split of ``codec.py`` (OO-assessment P3.10). The
reader ingests a packet stream and materializes it back into a ``.tio``
file. Bodies are verbatim from the pre-split ``codec.py``.
"""
from __future__ import annotations

import json
import logging
import struct
import zlib
from pathlib import Path
from typing import BinaryIO, Iterator

import numpy as np

from ..enums import Compression, Polarity, Precision
from ..spectral_dataset import SpectralDataset, WrittenRun
from .packets import (
    BULK_MODE_V2_BLOBS_FEATURE,
    HEADER_SIZE,
    AccessUnit,
    ChannelData,
    PacketFlag,
    PacketHeader,
    PacketType,
    _AU_PREFIX_STRUCT,
    _CHANNEL_NAMELEN_STRUCT,
    _CHANNEL_SUFFIX_STRUCT,
    crc32c,
    is_known_packet_type,
    unpack_blob_mate_info,
    unpack_blob_name_tok,
    unpack_blob_ref_diff,
    unpack_string,
)
from ._common import (
    _WIRE_TO_POLARITY,
    _decode_wire_codec,
)
from ._writer import (
    _provenance_csv_split,
    _provenance_params_parse,
)

# Logger name pinned to the historical ``ttio.transport.codec`` so that
# callers / tests that filter on that name keep capturing reader debug
# records after the P3.10 module split.
_LOG = logging.getLogger("ttio.transport.codec")


def _scan_pattern_from_byte(b: int) -> str:
    """Inverse of :func:`_scan_pattern_to_byte`. Java parity:
    :meth:`TransportReader.scanPatternFromByte`.
    """
    return {0: "raster", 1: "meander", 2: "random"}.get(b, "raster")


class TransportReader:
    """Deserialize a transport byte stream.

    The low-level API :meth:`iter_packets` yields ``(header, payload)``
    pairs. The high-level :meth:`read_to_dataset` materializes the
    stream into a new ``.tio`` file.
    """

    def __init__(self, source: BinaryIO | str | Path):
        """Construct a reader sourcing from a path or file-like object.

        Parameters
        ----------
        source : BinaryIO, str, or pathlib.Path
            A path is opened in ``"rb"`` mode and closed by
            :meth:`close`; a file-like object is borrowed (not
            closed by the reader).
        """
        self._owns_stream = isinstance(source, (str, Path))
        if self._owns_stream:
            self._stream: BinaryIO = open(source, "rb")  # noqa: SIM115
        else:
            self._stream = source  # type: ignore[assignment]

    def __enter__(self) -> "TransportReader":
        """Return ``self`` so the reader can be used as a context manager."""
        return self

    def __exit__(self, *exc: object) -> None:
        """Close the reader on context exit (delegates to :meth:`close`)."""
        self.close()

    def close(self) -> None:
        """Close the underlying stream if the reader opened it.

        No-op when the caller passed an externally-managed file-like
        object at construction. Safe to call more than once.
        """
        if self._owns_stream and not self._stream.closed:
            self._stream.close()

    def iter_packets(self) -> Iterator[tuple[PacketHeader, bytes]]:
        """Yield ``(header, payload)`` for every packet in the stream.

        Reads sequentially, verifying any CRC-32C trailer when the
        ``HAS_CHECKSUM`` flag is set. Unknown packet types are
        logged at DEBUG level and still yielded (forward-compat per
        transport-spec §6). Iteration terminates after an
        ``END_OF_STREAM`` packet or when the stream is exhausted.

        Yields
        ------
        tuple of (PacketHeader, bytes)

        Raises
        ------
        ValueError
            On truncated header, truncated payload, missing CRC, or
            CRC mismatch.
        """
        while True:
            header_bytes = self._stream.read(HEADER_SIZE)
            if not header_bytes:
                return
            if len(header_bytes) < HEADER_SIZE:
                raise ValueError(
                    f"truncated header: {len(header_bytes)}/{HEADER_SIZE} bytes"
                )
            header = PacketHeader.from_bytes(header_bytes)
            payload = self._stream.read(header.payload_length)
            if len(payload) != header.payload_length:
                raise ValueError(
                    f"truncated payload: {len(payload)}/{header.payload_length}"
                )
            if header.flags & int(PacketFlag.HAS_CHECKSUM):
                crc_bytes = self._stream.read(4)
                if len(crc_bytes) != 4:
                    raise ValueError("truncated CRC-32C")
                (expected_crc,) = struct.unpack("<I", crc_bytes)
                actual_crc = crc32c(payload)
                if expected_crc != actual_crc:
                    raise ValueError(
                        f"CRC-32C mismatch on packet type "
                        f"0x{header.packet_type_byte:02x}: "
                        f"expected 0x{expected_crc:08x}, got 0x{actual_crc:08x}"
                    )
            # Forward-compat (transport-spec §6 / v0.11 task 0.5):
            # tolerate unknown packet types so v0.10 readers can ingest
            # v0.11+ streams. The header was length-prefixed so the
            # payload (and CRC if present) was already consumed above
            # — just log and yield it so the caller can see it.
            if not is_known_packet_type(header.packet_type):
                _LOG.debug(
                    "skipping unknown packet type 0x%02x",
                    header.packet_type_byte,
                )
            yield header, payload
            if header.packet_type == int(PacketType.END_OF_STREAM):
                return

    def records_for_test(self) -> list["PacketRecord"]:
        """Materialise :meth:`iter_packets` into a list of
        :class:`PacketRecord` values. Test-only inspection hook
        mirroring Java's ``recordsForTest`` (forward-compat skip-unknown
        path, v0.11 task 0.5). Underscore semantics: this is not
        part of the stable transport-reader surface — production
        consumers should drive :meth:`iter_packets` directly.
        """
        # Local import to avoid pulling the ingest module into the
        # codec import cycle in the hot path.
        from .ingest import PacketRecord
        return [PacketRecord(header=h, payload=p)
                for h, p in self.iter_packets()]

    def read_to_dataset(
        self,
        *,
        output_path: str | Path,
        provider: str = "hdf5",
    ) -> SpectralDataset:
        """Materialise the stream into a ``.tio`` file.

        Drains :meth:`iter_packets`, accumulates per-dataset state
        (spectral runs, genomic runs, references, images,
        identifications, quantifications, subjects, samples,
        provenance, encryption algorithm), then materialises the
        result via :meth:`SpectralDataset.write_minimal`. The
        returned dataset is freshly opened on the written file.

        Parameters
        ----------
        output_path : str or pathlib.Path
            Destination ``.tio`` path. Created if absent; truncated
            if present.
        provider : str, optional
            Storage provider name. Default ``"hdf5"``.

        Returns
        -------
        SpectralDataset
            Reopened on the freshly written container.

        Raises
        ------
        ValueError
            On duplicate ``StreamHeader``, missing ``StreamHeader``,
            non-monotonic ``au_sequence``, malformed prelude
            packets, or a ``bulk_v2_blobs`` feature declaration
            without matching ``BlobV2*`` packets.
        """
        stream_meta: dict = {}
        dataset_metas: dict[int, dict] = {}
        run_data: dict[int, dict] = {}
        genomic_data: dict[int, dict] = {}
        # Phase 2c-T: per-genomic-dataset_id buffers for verbatim
        # v2 blobs. Populated when BlobV2* packets arrive, drained
        # into the WrittenGenomicRun before write_minimal.
        bulk_blobs: dict[int, dict] = {}
        last_seq: dict[int, int] = {}
        saw_stream_header = False
        bulk_mode_required = False
        # v0.11 Stage 1 / Task 2.3 — per-stream accumulator state for
        # the REFERENCE_GROUP_HEADER -> N x REFERENCE_CHROMOSOME ->
        # END_OF_REFERENCE_GROUP packet sequence. Java parity:
        # TransportReader.currentRefUri / currentChromNames /
        # currentChromSeqs / collectedRefs.
        current_ref_uri: str | None = None
        current_chrom_names: list[str] = []
        current_chrom_seqs: list[bytes] = []
        collected_refs: list = []  # list[ReferenceImport]
        # v0.11 Task 2.4: dataset-level @encrypted algorithm string
        # carried by ENCRYPTION_ALGORITHM (0x1B) packets. ``None`` when
        # no such packet appears in the stream. Multiple 0x1B packets
        # are tolerated — last-write-wins (spec §5.4 says "zero or
        # more"; in practice the writer emits exactly one). Java
        # parity: TransportReader.collectedEncryptionAlgorithm (commit
        # 530a5833).
        collected_encryption_algorithm: str | None = None
        # v0.11 Task 2.5: dataset-level provenance chain decoded from
        # DATASET_PROVENANCE (0x18) packets. Multiple 0x18 packets MAY
        # appear in a stream (spec §5.4 "zero or more"); each carries
        # its own record_count + records and they accumulate in
        # emission order. Passed into ``write_minimal`` as the
        # ``provenance`` kwarg so the on-disk
        # ``/study/provenance_json`` attribute round-trips. Java
        # parity: TransportReader.collectedProvenance (commit
        # 563e09c3).
        collected_provenance: list = []  # list[ProvenanceRecord]
        # v0.11 Task 2.6: per-stream image-cube accumulator state for
        # the IMAGE_HEADER (0x13) -> N x IMAGE_PIXEL (0x14) ->
        # END_OF_IMAGE (0x15) packet sequence. ``current_image_builder``
        # is non-None between IMAGE_HEADER and END_OF_IMAGE. Java
        # parity: TransportReader.currentImageBuilder /
        # collectedImage (commit a6b1e5d9).
        current_image_builder: dict | None = None
        collected_image = None  # type: MSImage | None
        # v0.11 Task 5.3 (Deferral 1): per-modality image accumulators.
        # The IMAGE_HEADER's modality byte selects which collected_*
        # slot is filled at END_OF_IMAGE time. ``current_image_skipping``
        # is True between IMAGE_HEADER (unknown modality) and the
        # matching END_OF_IMAGE so following IMAGE_PIXEL packets are
        # silently dropped (forward-compat per §4.16).
        collected_raman_image = None  # type: RamanImage | None
        collected_ir_image = None  # type: IRImage | None
        current_image_skipping = False
        # v0.11 Task 2.7: identification / quantification rows decoded
        # from IDENTIFICATIONS_TABLE (0x16) / QUANTIFICATIONS_TABLE
        # (0x17) packets. Multiple 0x16 / 0x17 packets MAY appear in a
        # stream (spec §5.4 step 6 says "zero or more"); rows accumulate
        # in emission order. Passed into ``write_minimal`` as the
        # ``identifications`` / ``quantifications`` kwargs so the
        # on-disk study compound datasets round-trip. Java parity:
        # TransportReader.collectedIdentifications /
        # collectedQuantifications (commit a6faab16).
        collected_identifications: list = []  # list[Identification]
        collected_quantifications: list = []  # list[Quantification]
        # Stage 6 / Task 6.3 (transport-spec v0.11 §4.22): subject +
        # sample rows decoded from SUBJECT_METADATA (0x19) /
        # SAMPLE_METADATA (0x1A) packets. Multiple 0x19 / 0x1A packets
        # MAY appear in a stream (spec §5.4 step 5: "zero or more");
        # rows accumulate in emission order. Passed into ``write_minimal``
        # as the ``subjects`` / ``samples`` kwargs so the on-disk
        # per-row groups round-trip via SpectralDataset's lazy properties.
        # Java parity: TransportReader.collectedSubjects /
        # collectedSamples (commit dd211600).
        collected_subjects: list = []  # list[Subject]
        collected_samples: list = []  # list[Sample]

        for header, payload in self.iter_packets():
            ptype = header.packet_type
            if ptype == int(PacketType.STREAM_HEADER):
                if saw_stream_header:
                    raise ValueError("duplicate StreamHeader")
                stream_meta = _decode_stream_header(payload)
                saw_stream_header = True
                bulk_mode_required = (
                    BULK_MODE_V2_BLOBS_FEATURE
                    in stream_meta.get("features", [])
                )
                continue
            if not saw_stream_header:
                raise ValueError(
                    f"first packet must be StreamHeader, got type 0x{ptype:02x}"
                )
            if ptype == int(PacketType.DATASET_HEADER):
                meta = _decode_dataset_header(payload)
                did = meta["dataset_id"]
                dataset_metas[did] = meta
                # genomic datasets get a parallel accumulator.
                if meta["spectrum_class"] == "TTIOGenomicRead":
                    genomic_data[did] = _new_genomic_accumulator()
                else:
                    run_data[did] = {
                        "channels": {c: [] for c in meta["channel_names"]},
                        "offsets": [],
                        "lengths": [],
                        "retention_times": [],
                        "ms_levels": [],
                        "polarities": [],
                        "precursor_mzs": [],
                        "precursor_charges": [],
                        "base_peak_intensities": [],
                        "running_offset": 0,
                    }
            elif ptype == int(PacketType.ACCESS_UNIT):
                did = header.dataset_id
                if did not in dataset_metas:
                    raise ValueError(
                        f"AccessUnit before DatasetHeader for id {did}"
                    )
                prev = last_seq.get(did, -1)
                if header.au_sequence <= prev:
                    raise ValueError(
                        f"non-monotonic au_sequence in dataset {did}: "
                        f"prev={prev}, got={header.au_sequence}"
                    )
                last_seq[did] = header.au_sequence
                if did in genomic_data:
                    _ingest_genomic_access_unit_bytes(
                        genomic_data[did], payload
                    )
                else:
                    _ingest_access_unit_bytes(run_data[did], payload)
            elif ptype == int(PacketType.BLOB_V2_MATE_INFO):
                # Phase 2c-T: verbatim mate_info blob for one
                # genomic dataset_id. At most one per dataset.
                ds_id, chrom_names, blob = unpack_blob_mate_info(payload)
                if ds_id != header.dataset_id:
                    raise ValueError(
                        f"BlobV2MateInfo dataset_id {ds_id} != "
                        f"header.dataset_id {header.dataset_id}"
                    )
                slot = bulk_blobs.setdefault(ds_id, {})
                if "mate_info" in slot:
                    raise ValueError(
                        f"duplicate BlobV2MateInfo for dataset_id {ds_id}"
                    )
                slot["mate_info"] = (chrom_names, blob)
            elif ptype == int(PacketType.BLOB_V2_REF_DIFF):
                ds_id, ref_uri, blob = unpack_blob_ref_diff(payload)
                if ds_id != header.dataset_id:
                    raise ValueError(
                        f"BlobV2RefDiff dataset_id {ds_id} != "
                        f"header.dataset_id {header.dataset_id}"
                    )
                slot = bulk_blobs.setdefault(ds_id, {})
                if "ref_diff" in slot:
                    raise ValueError(
                        f"duplicate BlobV2RefDiff for dataset_id {ds_id}"
                    )
                slot["ref_diff"] = (ref_uri, blob)
            elif ptype == int(PacketType.BLOB_V2_NAME_TOK):
                ds_id, blob = unpack_blob_name_tok(payload)
                if ds_id != header.dataset_id:
                    raise ValueError(
                        f"BlobV2NameTok dataset_id {ds_id} != "
                        f"header.dataset_id {header.dataset_id}"
                    )
                slot = bulk_blobs.setdefault(ds_id, {})
                if "name_tok" in slot:
                    raise ValueError(
                        f"duplicate BlobV2NameTok for dataset_id {ds_id}"
                    )
                slot["name_tok"] = blob
            elif ptype == int(PacketType.ENCRYPTION_ALGORITHM):
                # v0.11 Task 2.4: decode and stash the algorithm
                # string for application at materialize time. Multiple
                # 0x1B packets are tolerated (last-write-wins). Java
                # parity: TransportReader.decodeEncryptionAlgorithm.
                (algo_len,) = struct.unpack_from("<H", payload, 0)
                collected_encryption_algorithm = payload[
                    2:2 + algo_len
                ].decode("utf-8")
            elif ptype == int(PacketType.DATASET_PROVENANCE):
                # v0.11 Task 2.5: decode the per-packet record_count +
                # records and append to ``collected_provenance``.
                # Java parity: TransportReader.decodeDatasetProvenance.
                from ..provenance import ProvenanceRecord
                pl = payload
                off = 0
                (record_count,) = struct.unpack_from("<I", pl, off)
                off += 4
                for _ in range(record_count):
                    (timestamp,) = struct.unpack_from("<q", pl, off)
                    off += 8
                    (sw_len,) = struct.unpack_from("<H", pl, off)
                    off += 2
                    software = pl[off:off + sw_len].decode("utf-8")
                    off += sw_len
                    (pj_len,) = struct.unpack_from("<H", pl, off)
                    off += 2
                    params_json = pl[off:off + pj_len].decode("utf-8")
                    off += pj_len
                    (in_len,) = struct.unpack_from("<H", pl, off)
                    off += 2
                    inputs_csv = pl[off:off + in_len].decode("utf-8")
                    off += in_len
                    (out_len,) = struct.unpack_from("<H", pl, off)
                    off += 2
                    outputs_csv = pl[off:off + out_len].decode("utf-8")
                    off += out_len
                    collected_provenance.append(ProvenanceRecord(
                        timestamp_unix=int(timestamp),
                        software=software,
                        parameters=_provenance_params_parse(params_json),
                        input_refs=_provenance_csv_split(inputs_csv),
                        output_refs=_provenance_csv_split(outputs_csv),
                    ))
            elif ptype == int(PacketType.REFERENCE_GROUP_HEADER):
                # v0.11 Stage 1 / Task 2.3: decode header and prime the
                # per-group accumulator. chromosome_count / total_bases
                # / md5_hex are parsed for buffer-position advance only —
                # the actual values come from the per-chromosome
                # accumulator (ReferenceImport.__post_init__ recomputes
                # MD5 from the sequences).
                pl = payload
                off = 0
                (uri_len,) = struct.unpack_from("<H", pl, off)
                off += 2
                current_ref_uri = pl[off:off + uri_len].decode("utf-8")
                off += uri_len
                # chromosome_count (uint32) — advance only.
                off += 4
                # total_bases (uint64) — advance only.
                off += 8
                # md5_hex[32] — advance only.
                off += 32
                current_chrom_names = []
                current_chrom_seqs = []
            elif ptype == int(PacketType.REFERENCE_CHROMOSOME):
                # v0.11 Stage 1 / Task 2.3: decode chromosome + append.
                pl = payload
                off = 0
                (name_len,) = struct.unpack_from("<H", pl, off)
                off += 2
                name = pl[off:off + name_len].decode("utf-8")
                off += name_len
                (seq_len,) = struct.unpack_from("<Q", pl, off)
                off += 8
                encoding = pl[off]
                off += 1
                (data_len,) = struct.unpack_from("<I", pl, off)
                off += 4
                raw = bytes(pl[off:off + data_len])
                off += data_len
                if encoding == 0:
                    seq_bytes = raw
                elif encoding == 1:
                    seq_bytes = zlib.decompress(raw)
                    if len(seq_bytes) != seq_len:
                        raise ValueError(
                            f"REFERENCE_CHROMOSOME zlib payload inflated to "
                            f"{len(seq_bytes)} bytes; expected {seq_len}"
                        )
                else:
                    raise ValueError(
                        f"unknown REFERENCE_CHROMOSOME encoding: {encoding}"
                    )
                current_chrom_names.append(name)
                current_chrom_seqs.append(seq_bytes)
            elif ptype == int(PacketType.END_OF_REFERENCE_GROUP):
                # v0.11 Stage 1 / Task 2.3: close out the current group.
                if current_ref_uri is None:
                    raise ValueError(
                        "END_OF_REFERENCE_GROUP without prior "
                        "REFERENCE_GROUP_HEADER"
                    )
                # Local import to avoid pulling the genomic module into
                # the codec import cycle when no references are present.
                from ..genomic.reference_import import ReferenceImport
                collected_refs.append(ReferenceImport(
                    uri=current_ref_uri,
                    chromosomes=list(current_chrom_names),
                    sequences=list(current_chrom_seqs),
                ))
                current_ref_uri = None
                current_chrom_names = []
                current_chrom_seqs = []
            elif ptype == int(PacketType.IMAGE_HEADER):
                # v0.11 Task 2.6 / 5.1 (Deferral 1): decode an
                # IMAGE_HEADER (0x13) payload and prime the per-image
                # accumulator. Wire layout matches transport-spec
                # §4.16. Both continuous-mode (is_continuous == 1) and
                # processed-mode (is_continuous == 0, sparse
                # {channel,intensity} pairs indexed into the shared
                # axis) are supported; the mode flag is cached on the
                # builder and read by the IMAGE_PIXEL branch. Only
                # modality == 0 (MS) materialises today;
                # Raman/IR/UV-Vis modalities are handled by the
                # Task 5.3 modality-dispatch follow-up.
                pl = payload
                off = 0
                modality = pl[off]; off += 1
                (img_width,) = struct.unpack_from("<I", pl, off); off += 4
                (img_height,) = struct.unpack_from("<I", pl, off); off += 4
                (img_bins,) = struct.unpack_from("<I", pl, off); off += 4
                (img_px_x,) = struct.unpack_from("<d", pl, off); off += 8
                (img_px_y,) = struct.unpack_from("<d", pl, off); off += 8
                img_scan_byte = pl[off]; off += 1
                img_axis_kind = pl[off]; off += 1
                (img_axis_len,) = struct.unpack_from("<I", pl, off); off += 4
                img_axis = np.frombuffer(
                    pl, dtype="<f8", count=img_axis_len, offset=off
                ).copy()
                off += 8 * img_axis_len
                img_continuous = pl[off]; off += 1
                (title_len,) = struct.unpack_from("<H", pl, off); off += 2
                img_title = pl[off:off + title_len].decode("utf-8")
                off += title_len
                (isa_len,) = struct.unpack_from("<H", pl, off); off += 2
                img_isa = pl[off:off + isa_len].decode("utf-8")
                off += isa_len
                # v0.11 Task 5.3: optional modality_extras slot at the
                # tail. v0.10 streams emitted no such field, so older
                # fixtures may end at `isa_id`. Probe length to stay
                # backwards-compatible.
                if off + 2 <= len(pl):
                    (extras_len,) = struct.unpack_from("<H", pl, off)
                    off += 2
                    img_extras = bytes(pl[off:off + extras_len])
                    off += extras_len
                else:
                    img_extras = b""
                if img_continuous not in (0, 1):
                    raise ValueError(
                        f"IMAGE_HEADER: is_continuous must be 0 or 1; "
                        f"got {img_continuous}"
                    )
                # v0.11 Task 5.3: modality dispatch. modality 0=MS, 1=Raman,
                # 2=IR. Unknown modalities are logged + skipped (forward
                # compat per §4.16); the self-describing extras_len has
                # already advanced the buffer past the IMAGE_HEADER.
                if modality == 0:
                    builder_modality = 0
                    builder_extras = {}
                elif modality == 1:
                    if len(img_extras) != 16:
                        raise ValueError(
                            f"IMAGE_HEADER (modality=1, Raman) expects "
                            f"16-byte modality_extras (excitation + "
                            f"laser_power); got {len(img_extras)}"
                        )
                    exc, laser = struct.unpack("<dd", img_extras)
                    builder_modality = 1
                    builder_extras = {
                        "excitation_wavelength_nm": float(exc),
                        "laser_power_mw": float(laser),
                    }
                elif modality == 2:
                    if len(img_extras) != 9:
                        raise ValueError(
                            f"IMAGE_HEADER (modality=2, IR) expects "
                            f"9-byte modality_extras (ir_mode + "
                            f"resolution); got {len(img_extras)}"
                        )
                    ir_mode_byte, resolution = struct.unpack(
                        "<Bd", img_extras
                    )
                    builder_modality = 2
                    builder_extras = {
                        "ir_mode_byte": int(ir_mode_byte),
                        "resolution_cm_inv": float(resolution),
                    }
                else:
                    _LOG.warning(
                        "IMAGE_HEADER: unknown modality=%d; skipping "
                        "image block (extras_len=%d, width=%d, height=%d)",
                        modality, len(img_extras), int(img_width),
                        int(img_height),
                    )
                    current_image_skipping = True
                    current_image_builder = None
                    continue

                current_image_builder = {
                    "modality": builder_modality,
                    "extras": builder_extras,
                    "width": int(img_width),
                    "height": int(img_height),
                    "spectral_points": int(img_bins),
                    "pixel_size_x": float(img_px_x),
                    "pixel_size_y": float(img_px_y),
                    "scan_pattern": _scan_pattern_from_byte(img_scan_byte),
                    "axis_kind": int(img_axis_kind),
                    "axis": img_axis,
                    "title": img_title,
                    "isa_investigation_id": img_isa,
                    "is_continuous": bool(img_continuous),
                    "cube": np.zeros(
                        (int(img_height), int(img_width), int(img_bins)),
                        dtype=np.float64,
                    ),
                    "seen": np.zeros(
                        int(img_height) * int(img_width), dtype=bool
                    ),
                    "seen_count": 0,
                }
            elif ptype == int(PacketType.IMAGE_PIXEL):
                # v0.11 Task 2.6 / 5.1 (Deferral 1): decode an
                # IMAGE_PIXEL (0x14) payload per §4.17 and stash the
                # intensities at the pixel's (y, x) slot. The wire
                # shape inside ``payload_bytes`` branches on the
                # cached ``is_continuous`` from the IMAGE_HEADER:
                #
                #   * continuous: dense ``spectrum_bins`` intensities
                #   * processed:  ``u32 nonzero_count`` + that many
                #                 ``u32 channel_index + fXX intensity``
                #                 pairs; unmentioned channels stay 0.0.
                #
                # Java parity: TransportReader.appendPixel.
                if current_image_skipping:
                    # v0.11 Task 5.3: unknown-modality stream — silently
                    # drop the pixel until the matching END_OF_IMAGE.
                    continue
                if current_image_builder is None:
                    raise ValueError(
                        "IMAGE_PIXEL received before IMAGE_HEADER"
                    )
                pl = payload
                off = 0
                (px_x,) = struct.unpack_from("<I", pl, off); off += 4
                (px_y,) = struct.unpack_from("<I", pl, off); off += 4
                precision = pl[off]; off += 1
                compression = pl[off]; off += 1
                (payload_len,) = struct.unpack_from("<I", pl, off); off += 4
                raw = bytes(pl[off:off + payload_len])
                # §4.17 pixel enum: 0=none, 1=zstd, 2=zlib. Writers
                # emit 0 today; the inflate paths make the enum real.
                if compression == 1:
                    raw = _zstd_decompress(raw, 1 << 27)
                elif compression == 2:
                    raw = zlib.decompress(raw)
                elif compression != 0:
                    raise ValueError(
                        f"IMAGE_PIXEL compression={compression} not "
                        "supported (0=none, 1=zstd, 2=zlib)"
                    )
                if precision not in (0, 1):
                    raise ValueError(
                        f"IMAGE_PIXEL precision={precision} not supported "
                        "(expected 0=float32 or 1=float64)"
                    )
                w_img = current_image_builder["width"]
                h_img = current_image_builder["height"]
                sp_img = current_image_builder["spectral_points"]
                is_cont = current_image_builder["is_continuous"]
                if px_x >= w_img or px_y >= h_img:
                    raise ValueError(
                        f"IMAGE_PIXEL coordinates out of bounds: x={px_x}, "
                        f"y={px_y} (width={w_img}, height={h_img})"
                    )
                if is_cont:
                    if precision == 1:
                        intensities = np.frombuffer(
                            raw, dtype="<f8"
                        ).astype(np.float64, copy=True)
                    else:
                        intensities = np.frombuffer(
                            raw, dtype="<f4"
                        ).astype(np.float64, copy=True)
                    if intensities.size != sp_img:
                        raise ValueError(
                            f"IMAGE_PIXEL intensity count {intensities.size} "
                            f"does not match IMAGE_HEADER.spectrum_bins={sp_img}"
                        )
                else:
                    # Processed-mode: u32 nonzero_count + entries.
                    intensities = np.zeros(sp_img, dtype=np.float64)
                    if payload_len < 4:
                        raise ValueError(
                            "IMAGE_PIXEL (processed) payload too short to "
                            "carry nonzero_count"
                        )
                    (nonzero_count,) = struct.unpack_from("<I", raw, 0)
                    entry_off = 4
                    val_size = 8 if precision == 1 else 4
                    val_fmt = "<d" if precision == 1 else "<f"
                    expected_payload_len = 4 + nonzero_count * (4 + val_size)
                    if payload_len != expected_payload_len:
                        raise ValueError(
                            f"IMAGE_PIXEL (processed) payload_length "
                            f"{payload_len} does not match nonzero_count="
                            f"{nonzero_count} (expected "
                            f"{expected_payload_len})"
                        )
                    for _ in range(nonzero_count):
                        (ch,) = struct.unpack_from(
                            "<I", raw, entry_off
                        )
                        entry_off += 4
                        (v,) = struct.unpack_from(
                            val_fmt, raw, entry_off
                        )
                        entry_off += val_size
                        if ch >= sp_img:
                            raise ValueError(
                                f"IMAGE_PIXEL (processed) channel_index "
                                f"{ch} out of range [0, {sp_img}) at "
                                f"pixel (x={px_x}, y={px_y})"
                            )
                        intensities[int(ch)] = float(v)
                pixel_idx = int(px_y) * w_img + int(px_x)
                if current_image_builder["seen"][pixel_idx]:
                    raise ValueError(
                        f"duplicate IMAGE_PIXEL at (x={px_x}, y={px_y})"
                    )
                current_image_builder["seen"][pixel_idx] = True
                current_image_builder["seen_count"] += 1
                current_image_builder["cube"][int(px_y), int(px_x), :] = (
                    intensities
                )
            elif ptype == int(PacketType.END_OF_IMAGE):
                # v0.11 Task 2.6 / 5.3 (Deferral 1): close out the
                # current image cube on END_OF_IMAGE (0x15). Verifies
                # pixel_count_seen against the per-pixel ingest count
                # + width*height. Dispatches to the per-modality
                # collected_* slot (MS / Raman / IR). Java parity:
                # TransportReader.finishImage.
                if current_image_skipping:
                    # v0.11 Task 5.3: drain the END_OF_IMAGE of a
                    # skipped (unknown-modality) block. The
                    # pixel_count_seen field is still consumed for
                    # stream hygiene but not validated against any
                    # per-pixel count (we never accumulated one).
                    current_image_skipping = False
                    continue
                if current_image_builder is None:
                    raise ValueError(
                        "END_OF_IMAGE without prior IMAGE_HEADER"
                    )
                (declared,) = struct.unpack_from("<I", payload, 0)
                actual = current_image_builder["seen_count"]
                w_img = current_image_builder["width"]
                h_img = current_image_builder["height"]
                if declared != actual:
                    raise ValueError(
                        f"END_OF_IMAGE pixel_count_seen mismatch: "
                        f"declared={declared}, actual={actual} "
                        f"(width*height={w_img * h_img})"
                    )
                if actual != w_img * h_img:
                    raise ValueError(
                        f"END_OF_IMAGE pixel count {actual} does not "
                        f"equal width*height={w_img * h_img}"
                    )
                builder_modality = current_image_builder["modality"]
                builder_extras = current_image_builder["extras"]
                if builder_modality == 0:
                    from ..ms_image import MSImage
                    collected_image = MSImage(
                        width=w_img,
                        height=h_img,
                        spectral_points=current_image_builder["spectral_points"],
                        intensity=current_image_builder["cube"],
                        mz_axis=current_image_builder["axis"],
                        pixel_size_x=current_image_builder["pixel_size_x"],
                        pixel_size_y=current_image_builder["pixel_size_y"],
                        scan_pattern=current_image_builder["scan_pattern"],
                        title=current_image_builder["title"],
                        isa_investigation_id=current_image_builder[
                            "isa_investigation_id"
                        ],
                    )
                elif builder_modality == 1:
                    from ..raman_image import RamanImage
                    collected_raman_image = RamanImage(
                        width=w_img,
                        height=h_img,
                        spectral_points=current_image_builder["spectral_points"],
                        intensity=current_image_builder["cube"],
                        wavenumbers=current_image_builder["axis"],
                        pixel_size_x=current_image_builder["pixel_size_x"],
                        pixel_size_y=current_image_builder["pixel_size_y"],
                        scan_pattern=current_image_builder["scan_pattern"],
                        excitation_wavelength_nm=builder_extras[
                            "excitation_wavelength_nm"
                        ],
                        laser_power_mw=builder_extras["laser_power_mw"],
                        title=current_image_builder["title"],
                        isa_investigation_id=current_image_builder[
                            "isa_investigation_id"
                        ],
                    )
                elif builder_modality == 2:
                    from ..enums import IRMode
                    from ..ir_image import IRImage
                    ir_mode = (IRMode.ABSORBANCE
                               if builder_extras["ir_mode_byte"] == 1
                               else IRMode.TRANSMITTANCE)
                    collected_ir_image = IRImage(
                        width=w_img,
                        height=h_img,
                        spectral_points=current_image_builder["spectral_points"],
                        intensity=current_image_builder["cube"],
                        wavenumbers=current_image_builder["axis"],
                        pixel_size_x=current_image_builder["pixel_size_x"],
                        pixel_size_y=current_image_builder["pixel_size_y"],
                        scan_pattern=current_image_builder["scan_pattern"],
                        mode=ir_mode,
                        resolution_cm_inv=builder_extras["resolution_cm_inv"],
                        title=current_image_builder["title"],
                        isa_investigation_id=current_image_builder[
                            "isa_investigation_id"
                        ],
                    )
                current_image_builder = None
            elif ptype == int(PacketType.IDENTIFICATIONS_TABLE):
                # v0.11 Task 2.7: decode the uint32-length-prefixed
                # Arrow IPC payload and append the resulting rows to
                # ``collected_identifications``. Multiple 0x16 packets
                # accumulate in emission order. Java parity:
                # TransportReader.decodeIdentificationsTable.
                from .arrow_ipc import decode_identifications
                (ipc_len,) = struct.unpack_from("<I", payload, 0)
                ipc_bytes = bytes(payload[4:4 + ipc_len])
                collected_identifications.extend(
                    decode_identifications(ipc_bytes)
                )
            elif ptype == int(PacketType.QUANTIFICATIONS_TABLE):
                # v0.11 Task 2.7: decode the uint32-length-prefixed
                # Arrow IPC payload and append the resulting rows to
                # ``collected_quantifications``. Java parity:
                # TransportReader.decodeQuantificationsTable.
                from .arrow_ipc import decode_quantifications
                (ipc_len,) = struct.unpack_from("<I", payload, 0)
                ipc_bytes = bytes(payload[4:4 + ipc_len])
                collected_quantifications.extend(
                    decode_quantifications(ipc_bytes)
                )
            elif ptype == int(PacketType.SUBJECT_METADATA):
                # Stage 6 / Task 6.3: decode the uint32-length-prefixed
                # Arrow IPC payload and append the resulting rows to
                # ``collected_subjects``. Java parity:
                # TransportReader.decodeSubjectMetadata.
                from .arrow_ipc import decode_subjects
                (ipc_len,) = struct.unpack_from("<I", payload, 0)
                ipc_bytes = bytes(payload[4:4 + ipc_len])
                collected_subjects.extend(decode_subjects(ipc_bytes))
            elif ptype == int(PacketType.SAMPLE_METADATA):
                # Stage 6 / Task 6.3: decode the uint32-length-prefixed
                # Arrow IPC payload and append the resulting rows to
                # ``collected_samples``. Java parity:
                # TransportReader.decodeSampleMetadata.
                from .arrow_ipc import decode_samples
                (ipc_len,) = struct.unpack_from("<I", payload, 0)
                ipc_bytes = bytes(payload[4:4 + ipc_len])
                collected_samples.extend(decode_samples(ipc_bytes))
            elif ptype == int(PacketType.END_OF_DATASET):
                continue
            elif ptype == int(PacketType.END_OF_STREAM):
                break
            else:
                # Annotation / Provenance / Chromatogram /
                # ProtectionMetadata — recognized but not yet
                # materialized (M70/M71 scope).
                continue

        # Phase 2c-T: fail closed when bulk-mode was declared but no
        # blobs arrived for any genomic dataset. The feature flag is
        # required, so a bulk stream that ships zero blobs is malformed.
        if bulk_mode_required and not bulk_blobs:
            raise ValueError(
                f"StreamHeader declared {BULK_MODE_V2_BLOBS_FEATURE!r} "
                "but no BlobV2* packets were received"
            )

        runs: dict[str, WrittenRun] = {}
        for did, meta in dataset_metas.items():
            if did in genomic_data:
                continue
            rd = run_data[did]
            channel_data = {
                c: (np.concatenate(rd["channels"][c])
                    if rd["channels"][c]
                    else np.array([], dtype="<f8"))
                for c in meta["channel_names"]
            }
            runs[meta["name"]] = WrittenRun(
                spectrum_class=meta["spectrum_class"],
                acquisition_mode=meta["acquisition_mode"],
                channel_data=channel_data,
                offsets=np.array(rd["offsets"], dtype="<u8"),
                lengths=np.array(rd["lengths"], dtype="<u4"),
                retention_times=np.array(rd["retention_times"], dtype="<f8"),
                ms_levels=np.array(rd["ms_levels"], dtype="<i4"),
                polarities=np.array(rd["polarities"], dtype="<i4"),
                precursor_mzs=np.array(rd["precursor_mzs"], dtype="<f8"),
                precursor_charges=np.array(rd["precursor_charges"], dtype="<i4"),
                base_peak_intensities=np.array(
                    rd["base_peak_intensities"], dtype="<f8"
                ),
                signal_compression="gzip",
            )

        # build WrittenGenomicRun for each genomic dataset.
        from ..written_genomic_run import BulkV2Blobs, WrittenGenomicRun
        genomic_runs: dict[str, WrittenGenomicRun] = {}
        for did, gd in genomic_data.items():
            meta = dataset_metas[did]
            n = len(gd["chromosomes"])
            instrument_meta = json.loads(meta.get("instrument_json") or "{}")
            # Phase 2c-T: attach any verbatim blobs collected for this
            # dataset_id so write_minimal skips the v2 codec encode
            # for those channels.
            slot = bulk_blobs.get(did, {})
            bulk_obj: BulkV2Blobs | None = None
            if slot:
                mate = slot.get("mate_info")
                rdif = slot.get("ref_diff")
                ntok = slot.get("name_tok")
                bulk_obj = BulkV2Blobs(
                    mate_info_blob=mate[1] if mate else None,
                    mate_info_chrom_names=list(mate[0]) if mate else None,
                    ref_diff_blob=rdif[1] if rdif else None,
                    ref_diff_reference_uri=rdif[0] if rdif else None,
                    name_tok_blob=ntok if ntok is not None else None,
                )
            genomic_runs[meta["name"]] = WrittenGenomicRun(
                acquisition_mode=meta["acquisition_mode"],
                reference_uri=instrument_meta.get("reference_uri", ""),
                platform=instrument_meta.get("platform", ""),
                sample_name=instrument_meta.get("sample_name", ""),
                # M97: "" on the wire means no role (attribute absent).
                read_role=instrument_meta.get("read_role") or None,
                positions=np.array(gd["positions"], dtype=np.int64),
                mapping_qualities=np.array(gd["mapping_qualities"], dtype=np.uint8),
                flags=np.array(gd["flags"], dtype=np.uint32),
                sequences=(np.concatenate(gd["sequences_chunks"])
                           if gd["sequences_chunks"]
                           else np.array([], dtype=np.uint8)),
                qualities=(np.concatenate(gd["qualities_chunks"])
                           if gd["qualities_chunks"]
                           else np.array([], dtype=np.uint8)),
                offsets=np.array(gd["offsets"], dtype=np.uint64),
                lengths=np.array(gd["lengths"], dtype=np.uint32),
                # compound fields now round-trip on the wire.
                # When the source is an M89.2-era stream the per-AU
                # decoders default the missing strings to "" and the
                # mate scalars to -1 / 0 (preserved by the AU
                # decoder + accumulator paths).
                cigars=list(gd["cigars"]) if gd["cigars"]
                        else ["" for _ in range(n)],
                read_names=list(gd["read_names"]) if gd["read_names"]
                            else ["" for _ in range(n)],
                mate_chromosomes=list(gd["mate_chromosomes"])
                                  if gd["mate_chromosomes"]
                                  else ["" for _ in range(n)],
                mate_positions=(np.array(gd["mate_positions"], dtype=np.int64)
                                if gd["mate_positions"]
                                else np.full(n, -1, dtype=np.int64)),
                template_lengths=(np.array(gd["template_lengths"], dtype=np.int32)
                                   if gd["template_lengths"]
                                   else np.zeros(n, dtype=np.int32)),
                chromosomes=list(gd["chromosomes"]),
                tags=list(gd["tags"]) if gd.get("any_tags") else None,
                bulk_v2_blobs=bulk_obj,
            )

        path = SpectralDataset.write_minimal(
            output_path,
            title=stream_meta.get("title", ""),
            isa_investigation_id=stream_meta.get("isa_investigation", ""),
            runs=runs,
            genomic_runs=genomic_runs or None,
            features=list(stream_meta.get("features", [])) or None,
            # v0.11 Task 2.5: surface decoded DATASET_PROVENANCE
            # records on the materialised .tio. ``None`` when no 0x18
            # packets appeared, matching write_minimal's signature.
            # Java parity: TransportReader.materializeTo passes
            # collectedProvenance into SpectralDataset.create.
            provenance=collected_provenance or None,
            # v0.11 Task 2.7: surface decoded IDENTIFICATIONS_TABLE /
            # QUANTIFICATIONS_TABLE rows on the materialised .tio.
            # ``None`` when the wire had no 0x16 / 0x17 packets, so
            # the matching study compound datasets are omitted.
            identifications=collected_identifications or None,
            quantifications=collected_quantifications or None,
            # Stage 6 / Task 6.3: surface decoded SUBJECT_METADATA /
            # SAMPLE_METADATA rows on the materialised .tio as per-row
            # groups under /study/subjects/ + /study/samples/. ``None``
            # when the wire had no 0x19 / 0x1A packets so the matching
            # groups are absent (pre-Stage-6 compat). Java parity:
            # TransportReader.materializeTo layering loop.
            subjects=collected_subjects or None,
            samples=collected_samples or None,
            provider=provider,
        )
        # v0.11 Stage 1 / Task 2.3: embed any reference groups decoded
        # from the stream's REFERENCE_* packets. ReferenceImport.write_-
        # to_dataset(ds) requires an open writable HDF5 provider, so
        # reopen the just-written .tio in r+ mode, embed, then close
        # before the final read-only open returned to the caller. Java
        # parity: TransportReader.materializeTo (commit 7f3dec46).
        if collected_refs:
            with SpectralDataset.open(path, writable=True) as ds_w:
                for ref in collected_refs:
                    ref.write_to_dataset(ds_w)
        # v0.11 Task 2.6 / 5.3 (Deferral 1): embed any image cubes
        # decoded from the stream's IMAGE_* packets. Each modality
        # lives in its own /study/{image_cube,raman_image_cube,
        # ir_image_cube} subgroup, so all three can coexist on the
        # same dataset. SpectralDataset caches the per-modality image
        # at construction time, so the close+reopen at the end forces
        # a fresh read of each cube on the returned handle. Java
        # parity: TransportReader.materializeTo (commit f99ec47d).
        if (collected_image is not None
                or collected_raman_image is not None
                or collected_ir_image is not None):
            with SpectralDataset.open(path, writable=True) as ds_w:
                study_grp = ds_w.provider.root_group().open_group("study")
                if collected_image is not None:
                    collected_image.write_to(study_grp)
                if collected_raman_image is not None:
                    collected_raman_image.write_to(study_grp)
                if collected_ir_image is not None:
                    collected_ir_image.write_to(study_grp)
        # v0.11 Task 2.4: persist the dataset-level @encrypted root
        # attribute so the materialised file reports is_encrypted ==
        # True on reopen. SpectralDataset caches encrypted_algorithm
        # in an instance field at construction time, so we must close
        # + reopen to surface the value on the returned dataset. Java
        # parity: TransportReader.materializeTo (commit 530a5833).
        if collected_encryption_algorithm is not None:
            with SpectralDataset.open(path, writable=True) as ds_w:
                ds_w.provider.root_group().set_attribute(
                    "encrypted", collected_encryption_algorithm
                )
        return SpectralDataset.open(path)


def _new_genomic_accumulator() -> dict:
    """Per-dataset accumulator for genomic AUs. See M89.2 + M90.9."""
    return {
        "chromosomes": [],
        "positions": [],
        "mapping_qualities": [],
        "flags": [],
        "sequences_chunks": [],
        "qualities_chunks": [],
        "offsets": [],
        "lengths": [],
        "running_offset": 0,
        # compound-field accumulators.
        "cigars": [],
        "read_names": [],
        "mate_chromosomes": [],
        "mate_positions": [],
        "template_lengths": [],
        # M101: per-read SAM tag text; any_tags marks a stream that
        # carried the tags channel at all.
        "tags": [],
        "any_tags": False,
    }


def _ingest_genomic_access_unit_bytes(gd: dict, payload: bytes) -> None:
    """Parse a genomic AU payload (spectrum_class==5) into ``gd``.

    Mirrors :func:`_ingest_access_unit_bytes` but extracts the genomic
    suffix (chromosome / position / mapq / flags) and accumulates
    sequences + qualities as concatenated uint8 buffers. M89.2.
    """
    au = AccessUnit.from_bytes(payload)
    if au.spectrum_class != 5:
        raise ValueError(
            f"genomic accumulator received spectrum_class {au.spectrum_class}"
        )
    gd["chromosomes"].append(au.chromosome)
    gd["positions"].append(int(au.position))
    gd["mapping_qualities"].append(int(au.mapping_quality))
    gd["flags"].append(int(au.flags) & 0xFFFFFFFF)
    # mate extension fields ride on the AU genomic suffix.
    gd["mate_positions"].append(int(au.mate_position))
    gd["template_lengths"].append(int(au.template_length))
    length = 0
    # compound-string channels default to "" if absent (an
    # M89.2-era AU). Channel-name dispatch covers both layouts.
    cigar_str = ""
    name_str = ""
    mate_chr_str = ""
    tags_str = ""
    for ch in au.channels:
        if ch.precision != int(Precision.UINT8):
            raise NotImplementedError(
                f"genomic channel precision {ch.precision} not yet supported "
                "(UINT8 only in M89.2)"
            )
        # dispatch on wire compression byte (NONE / RANS_*
        # / BASE_PACK). See _decode_wire_codec.
        decoded = _decode_wire_codec(bytes(ch.data), int(ch.compression))
        if ch.name == "sequences":
            arr = np.frombuffer(decoded, dtype=np.uint8).copy()
            gd["sequences_chunks"].append(arr)
            length = len(arr)
        elif ch.name == "qualities":
            arr = np.frombuffer(decoded, dtype=np.uint8).copy()
            gd["qualities_chunks"].append(arr)
            if length == 0:
                length = len(arr)
        elif ch.name == "cigar":
            cigar_str = decoded.decode("utf-8")
        elif ch.name == "read_name":
            name_str = decoded.decode("utf-8")
        elif ch.name == "mate_chromosome":
            mate_chr_str = decoded.decode("utf-8")
        elif ch.name == "tags":
            tags_str = decoded.decode("utf-8")
            gd["any_tags"] = True
    gd["tags"].append(tags_str)
    gd["cigars"].append(cigar_str)
    gd["read_names"].append(name_str)
    gd["mate_chromosomes"].append(mate_chr_str)
    gd["offsets"].append(gd["running_offset"])
    gd["lengths"].append(length)
    gd["running_offset"] += length


def _decode_stream_header(payload: bytes) -> dict:
    offset = 0
    format_version, offset = unpack_string(payload, offset, width=2)
    title, offset = unpack_string(payload, offset, width=2)
    isa_investigation, offset = unpack_string(payload, offset, width=2)
    (n_features,) = struct.unpack_from("<H", payload, offset)
    offset += 2
    features: list[str] = []
    for _ in range(n_features):
        f, offset = unpack_string(payload, offset, width=2)
        features.append(f)
    (n_datasets,) = struct.unpack_from("<H", payload, offset)
    offset += 2
    return {
        "format_version": format_version,
        "title": title,
        "isa_investigation": isa_investigation,
        "features": features,
        "n_datasets": n_datasets,
    }


def _decode_dataset_header(payload: bytes) -> dict:
    offset = 0
    (dataset_id,) = struct.unpack_from("<H", payload, offset)
    offset += 2
    name, offset = unpack_string(payload, offset, width=2)
    (acquisition_mode,) = struct.unpack_from("<B", payload, offset)
    offset += 1
    spectrum_class, offset = unpack_string(payload, offset, width=2)
    (n_channels,) = struct.unpack_from("<B", payload, offset)
    offset += 1
    channel_names: list[str] = []
    for _ in range(n_channels):
        c, offset = unpack_string(payload, offset, width=2)
        channel_names.append(c)
    instrument_json, offset = unpack_string(payload, offset, width=4)
    (expected_au_count,) = struct.unpack_from("<I", payload, offset)
    offset += 4
    return {
        "dataset_id": dataset_id,
        "name": name,
        "acquisition_mode": acquisition_mode,
        "spectrum_class": spectrum_class,
        "channel_names": channel_names,
        "instrument_json": instrument_json,
        "expected_au_count": expected_au_count,
    }


def _ingest_access_unit_bytes(rd: dict, payload: bytes) -> None:
    """Parse an AU payload directly into ``rd``, skipping dataclass construction.

    Equivalent to ``_ingest_access_unit(rd, AccessUnit.from_bytes(payload))`` but
    avoids creating 1 AccessUnit + N ChannelData dataclasses per AU.
    """
    if len(payload) < 38:
        raise ValueError(f"access unit payload too short: {len(payload)}")
    (
        _spectrum_class, _acq_mode, ms_level, polarity_wire,
        retention_time, precursor_mz,
        precursor_charge,
        _ion_mobility, base_peak_intensity,
        n_channels,
    ) = _AU_PREFIX_STRUCT.unpack_from(payload, 0)

    channel_map = rd["channels"]
    offset = 38
    length = 0
    seen: dict[str, bool] = {}
    for _ in range(n_channels):
        (name_len,) = _CHANNEL_NAMELEN_STRUCT.unpack_from(payload, offset)
        offset += 2
        name = bytes(payload[offset:offset + name_len]).decode("utf-8")
        offset += name_len
        precision, compression, n_elements, data_length = _CHANNEL_SUFFIX_STRUCT.unpack_from(
            payload, offset
        )
        offset += 10
        data = bytes(payload[offset:offset + data_length])
        offset += data_length
        if precision != _FLOAT64_WIRE:
            raise NotImplementedError(
                f"precision {precision} not yet supported (FLOAT64 only)"
            )
        arr = _decode_wire_channel(compression, data, n_elements)
        seen[name] = True
        if length == 0:
            length = len(arr)
        elif len(arr) != length:
            raise ValueError(
                f"channels in one AU have mismatched lengths: {length} vs {len(arr)}"
            )
        if name in channel_map:
            channel_map[name].append(arr)

    rd["offsets"].append(rd["running_offset"])
    rd["lengths"].append(length)
    rd["running_offset"] += length
    for cname, buckets in channel_map.items():
        if cname not in seen:
            buckets.append(np.zeros(length, dtype="<f8"))

    rd["retention_times"].append(retention_time)
    rd["ms_levels"].append(ms_level)
    rd["polarities"].append(int(_WIRE_TO_POLARITY.get(polarity_wire, Polarity.UNKNOWN)))
    rd["precursor_mzs"].append(precursor_mz)
    rd["precursor_charges"].append(precursor_charge)
    rd["base_peak_intensities"].append(base_peak_intensity)


_FLOAT64_WIRE = int(Precision.FLOAT64)
_COMPRESSION_NONE_WIRE = int(Compression.NONE)
_COMPRESSION_ZLIB_WIRE = int(Compression.ZLIB)
_COMPRESSION_ZSTD_WIRE = int(Compression.ZSTD)
_COMPRESSION_FDZ_WIRE = int(Compression.FLOAT_DELTA_ZSTD)


def _zstd_decompress(data: bytes, max_output_size: int) -> bytes:
    """One-shot zstd decompress. ``max_output_size`` backstops frames
    without an embedded content size (the channel header supplies the
    exact plaintext size)."""
    import zstandard

    return zstandard.ZstdDecompressor().decompress(
        data, max_output_size=max_output_size)


def _decode_wire_channel(compression: int, data: bytes, n_elements: int) -> np.ndarray:
    """Decode one spectral AU channel's ``data`` bytes to a fresh
    little-endian float64 array, dispatching on the wire compression
    byte (spec §4.3.1)."""
    if compression == _COMPRESSION_NONE_WIRE:
        raw = data
    elif compression == _COMPRESSION_ZLIB_WIRE:
        raw = zlib.decompress(data)
    elif compression == _COMPRESSION_ZSTD_WIRE:
        raw = _zstd_decompress(data, n_elements * 8)
    elif compression == _COMPRESSION_FDZ_WIRE:
        from ..codecs import float_delta_zstd
        arr = float_delta_zstd.decode(data)
        if len(arr) != n_elements:
            raise ValueError(
                f"FLOAT_DELTA_ZSTD channel decoded {len(arr)} values but the "
                f"channel header declares n_elements={n_elements}"
            )
        return arr
    else:
        raise NotImplementedError(
            f"compression {compression} not yet supported "
            "(current codec: NONE, ZLIB, ZSTD, FLOAT_DELTA_ZSTD)"
        )
    return np.frombuffer(raw, dtype="<f8").copy()


def _ingest_access_unit(rd: dict, au: AccessUnit) -> None:
    arr_by_name: dict[str, np.ndarray] = {}
    for ch in au.channels:
        if ch.precision != int(Precision.FLOAT64):
            raise NotImplementedError(
                f"precision {ch.precision} not yet supported (FLOAT64 only)"
            )
        arr_by_name[ch.name] = _decode_wire_channel(ch.compression, ch.data, ch.n_elements)

    lengths = {len(a) for a in arr_by_name.values()}
    if len(lengths) > 1:
        raise ValueError(
            f"channels in one AU have mismatched lengths: {sorted(lengths)}"
        )
    length = next(iter(lengths)) if lengths else 0

    rd["offsets"].append(rd["running_offset"])
    rd["lengths"].append(length)
    rd["running_offset"] += length
    for cname in rd["channels"]:
        arr = arr_by_name.get(cname)
        if arr is None:
            arr = np.zeros(length, dtype="<f8")
        rd["channels"][cname].append(arr)

    rd["retention_times"].append(au.retention_time)
    rd["ms_levels"].append(au.ms_level)
    rd["polarities"].append(int(_WIRE_TO_POLARITY.get(au.polarity, Polarity.UNKNOWN)))
    rd["precursor_mzs"].append(au.precursor_mz)
    rd["precursor_charges"].append(au.precursor_charge)
    rd["base_peak_intensities"].append(au.base_peak_intensity)
