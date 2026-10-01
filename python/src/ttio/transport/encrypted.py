"""v1.0 encrypted transport paths.

Bypasses the high-level :class:`Spectrum` / :class:`SignalArray`
decode machinery (which requires the key to materialise) and reads
HDF5 segment compounds directly, so encrypted bytes travel through
transport without being decrypted in transit.

Cross-language equivalents: ObjC
``TTIOEncryptedTransport`` · Java
``global.thalion.ttio.transport.EncryptedTransport``.
"""
from __future__ import annotations

import json
import struct
from pathlib import Path
from typing import BinaryIO

import numpy as np

from .. import _hdf5_io as io
from ..enums import Precision
from ..feature_flags import (
    OPT_ENCRYPTED_AU_HEADERS,
    OPT_PER_AU_ENCRYPTION,
)
from .codec import (
    TransportWriter,
    _SPECTRUM_CLASS_TO_WIRE,
    _instrument_config_json,  # reused via wrapper below
    unpack_string,
)
from .packets import (
    AccessUnit,
    ChannelData,
    PacketFlag,
    PacketHeader,
    PacketType,
    TRANSPORT_BLOCKS_V1_FEATURE,
    now_ns,
    pack_string,
)


def is_per_au_encrypted(path: str | Path,
                          *,
                          provider: str | None = None) -> bool:
    """Return True if the file on disk carries ``opt_per_au_encryption``.

    Routes through :func:`ttio.providers.open_provider` so any
    backend works uniformly."""
    from ..providers.registry import open_provider
    sp = open_provider(str(path), provider=provider, mode="r")
    try:
        _, features = io.read_feature_flags(sp.root_group())
    finally:
        sp.close()
    return OPT_PER_AU_ENCRYPTION in features


def write_encrypted_dataset(
    writer: TransportWriter,
    ttio_path: str | Path,
    *,
    provider: str | None = None,
) -> None:
    """Emit a full transport stream from a per-AU-encrypted .tio
    file. Bypasses :class:`SpectralDataset.open` (which would refuse
    the file because the plaintext ``<channel>_values`` are absent)
    and walks the source via the StorageProvider abstraction so every
    backend that supports VL_BYTES compound reads works.

    ProtectionMetadata is emitted after ``StreamHeader``; encrypted
    AUs carry the flag bits defined in ``docs/transport-spec.md``
    §3.1.1.
    """
    from ..providers.registry import open_provider

    sp = open_provider(str(ttio_path), provider=provider, mode="r")
    try:
        root = sp.root_group()
        _, features = io.read_feature_flags(root)
        if OPT_PER_AU_ENCRYPTION not in features:
            raise ValueError(
                f"{ttio_path!r} does not carry opt_per_au_encryption"
            )
        headers_encrypted = OPT_ENCRYPTED_AU_HEADERS in features

        study = root.open_group("study")
        title = io.read_string_attr(study, "title") or ""
        isa = io.read_string_attr(study, "isa_investigation_id") or ""
        ms_runs = study.open_group("ms_runs")
        run_items = [(n, ms_runs.open_group(n))
                      for n in ms_runs.child_names()
                      if not n.startswith("_") and ms_runs.has_child(n)]
        # also walk genomic_runs after MS. dataset_id_counter
        # continues from MS so AAD reconstruction matches the
        # per-AU encrypt path (M90.1).
        genomic_layouts: dict[str, str] = {}
        if study.has_child("genomic_runs"):
            g_runs_group = study.open_group("genomic_runs")
            genomic_run_items = [(n, g_runs_group.open_group(n))
                                  for n in g_runs_group.child_names()
                                  if not n.startswith("_")
                                  and g_runs_group.has_child(n)]
            for run_name, run_group in genomic_run_items:
                layout = ""
                if run_group.has_attribute("layout"):
                    raw = run_group.get_attribute("layout")
                    layout = (bytes(raw).decode("utf-8")
                              if isinstance(raw, (bytes, bytearray))
                              else str(raw))
                genomic_layouts[run_name] = layout
        else:
            genomic_run_items = []

        # blocks_v1 runs ride the same per-read AU stream plus the
        # M99.1 sidecar packets (GenomicRunSidecar + one BlockSidecar
        # per block); the StreamHeader announces them with a wire
        # feature token so receivers know the stream needs them.
        stream_features = list(features)
        if any(v == "blocks_v1" for v in genomic_layouts.values()):
            stream_features.append(TRANSPORT_BLOCKS_V1_FEATURE)

        writer.write_stream_header(
            format_version="1.2",
            title=title,
            isa_investigation=isa,
            features=stream_features,
            n_datasets=len(run_items) + len(genomic_run_items),
        )

        # One ProtectionMetadata per dataset (per run). The wrapped
        # DEK is read from the channel's @<channel>_wrapped_dek attr
        # if present; otherwise an empty byte string is emitted. The
        # receiver is responsible for KEK unwrap via out-of-band
        # key management.
        for dataset_id, (run_name, run_group) in enumerate(run_items, start=1):
            sig = run_group.open_group("signal_channels")
            # Probe the first channel for the algorithm / wrapped DEK
            # metadata. All channels in a run share the same DEK in
            # the v1.0 design.
            channel_names = [
                c for c in (io.read_string_attr(sig, "channel_names") or "").split(",")
                if c
            ]
            first_channel = channel_names[0] if channel_names else "intensity"
            cipher_suite = (io.read_string_attr(sig, f"{first_channel}_algorithm")
                              or "aes-256-gcm")
            kek_algorithm = (io.read_string_attr(sig, f"{first_channel}_kek_algorithm")
                               or "")
            wrapped_dek_attr = f"{first_channel}_wrapped_dek"
            if sig.has_attribute(wrapped_dek_attr):
                wrapped_dek = bytes(sig.get_attribute(wrapped_dek_attr))
            else:
                wrapped_dek = b""
            _emit_protection_metadata(
                writer,
                dataset_id=dataset_id,
                cipher_suite=cipher_suite,
                kek_algorithm=kek_algorithm,
                wrapped_dek=wrapped_dek,
                signature_algorithm="",
                public_key=b"",
                additional_recipients=_read_additional_recipients(
                    sig, first_channel),
                server_kek_id=_read_server_kek_id(sig, first_channel),
            )

            spectrum_class = (io.read_string_attr(run_group, "spectrum_class")
                              or "TTIOMassSpectrum")
            acquisition_mode = io.read_int_attr(run_group, "acquisition_mode",
                                                  default=0) or 0
            # count discovered from any channel's segment dataset length
            first_segs = io.read_channel_segments(sig, f"{first_channel}_segments")
            n_spectra = len(first_segs)

            writer.write_dataset_header(
                dataset_id=dataset_id,
                name=run_name,
                acquisition_mode=int(acquisition_mode),
                spectrum_class=spectrum_class,
                channel_names=list(channel_names),
                instrument_json="{}",   # skipped for brevity; non-PHI
                expected_au_count=n_spectra,
            )

        # Emit AUs for each run
        for dataset_id, (run_name, run_group) in enumerate(run_items, start=1):
            sig = run_group.open_group("signal_channels")
            idx = run_group.open_group("spectrum_index")
            channel_names = [
                c for c in (io.read_string_attr(sig, "channel_names") or "").split(",")
                if c
            ]
            # Pre-load all segments for each channel + header segments.
            channel_segments_by_name = {
                c: io.read_channel_segments(sig, f"{c}_segments")
                for c in channel_names
            }
            if headers_encrypted:
                header_segs = io.read_au_header_segments(idx, "au_header_segments")
            else:
                header_segs = None
                # Read plaintext index arrays via the provider.
                def _read_or_none(name):
                    return (idx.open_dataset(name).read()
                             if idx.has_child(name) else None)
                rts = _read_or_none("retention_times")
                ms_levels = _read_or_none("ms_levels")
                polarities = _read_or_none("polarities")
                precursor_mzs = _read_or_none("precursor_mzs")
                precursor_charges = _read_or_none("precursor_charges")
                base_peak = _read_or_none("base_peak_intensities")

            n = len(next(iter(channel_segments_by_name.values())))
            wire_class = _SPECTRUM_CLASS_TO_WIRE.get(spectrum_class, 0)

            for i in range(n):
                # Build encrypted ChannelData list.
                channels = []
                for cname in channel_names:
                    seg = channel_segments_by_name[cname][i]
                    data = seg.iv + seg.tag + seg.ciphertext
                    channels.append(ChannelData(
                        name=cname,
                        precision=1,     # float64
                        compression=0,   # NONE (inner plaintext would be f64 raw)
                        n_elements=seg.length,
                        data=data,
                    ))

                if headers_encrypted:
                    # Wire: spectrum_class(u8) n_channels(u8) IV(12) TAG(16)
                    # ct(36) [channels] [pixel optional].
                    hdr_seg = header_segs[i]
                    payload = (
                        struct.pack("<BB", int(wire_class) & 0xFF, len(channels) & 0xFF)
                        + hdr_seg.iv + hdr_seg.tag + hdr_seg.ciphertext
                    )
                    for ch in channels:
                        payload += ch.to_bytes()
                    _emit_raw_au(
                        writer,
                        dataset_id=dataset_id,
                        au_sequence=i,
                        payload=payload,
                        flags=int(PacketFlag.ENCRYPTED) | int(PacketFlag.ENCRYPTED_HEADER),
                    )
                else:
                    # Plaintext filter header, encrypted channels.
                    au = AccessUnit(
                        spectrum_class=wire_class,
                        acquisition_mode=int(acquisition_mode),
                        ms_level=int(ms_levels[i]) if ms_levels is not None else 0,
                        polarity=_wire_polarity(int(polarities[i])) if polarities is not None else 2,
                        retention_time=float(rts[i]) if rts is not None else 0.0,
                        precursor_mz=float(precursor_mzs[i]) if precursor_mzs is not None else 0.0,
                        precursor_charge=int(precursor_charges[i]) if precursor_charges is not None else 0,
                        ion_mobility=0.0,
                        base_peak_intensity=float(base_peak[i]) if base_peak is not None else 0.0,
                        channels=channels,
                    )
                    _emit_raw_au(
                        writer,
                        dataset_id=dataset_id,
                        au_sequence=i,
                        payload=au.to_bytes(),
                        flags=int(PacketFlag.ENCRYPTED),
                    )

            writer.write_end_of_dataset(
                dataset_id=dataset_id,
                final_au_sequence=n,
            )

        # emit genomic_runs after MS. Same dataset_id space
        # (continues from MS) so AAD reconstruction stays symmetric.
        for genomic_offset, (g_run_name, g_run_group) in enumerate(
            genomic_run_items, start=1,
        ):
            g_dataset_id = len(run_items) + genomic_offset
            g_sig = g_run_group.open_group("signal_channels")
            g_idx = g_run_group.open_group("genomic_index")
            # Genomic encrypts sequences + qualities (M90.1) and, when
            # the run carries them, the SAM tags (M101).
            g_channel_names = [c for c in ("sequences", "qualities", "tags")
                                if g_sig.has_child(f"{c}_segments")]
            g_first = g_channel_names[0] if g_channel_names else "sequences"
            cipher_suite = (io.read_string_attr(g_sig, f"{g_first}_algorithm")
                              or "aes-256-gcm")
            kek_algorithm = (io.read_string_attr(g_sig,
                                                   f"{g_first}_kek_algorithm")
                               or "")
            wrapped_dek_attr = f"{g_first}_wrapped_dek"
            if g_sig.has_attribute(wrapped_dek_attr):
                wrapped_dek = bytes(g_sig.get_attribute(wrapped_dek_attr))
            else:
                wrapped_dek = b""
            _emit_protection_metadata(
                writer,
                dataset_id=g_dataset_id,
                cipher_suite=cipher_suite,
                kek_algorithm=kek_algorithm,
                wrapped_dek=wrapped_dek,
                signature_algorithm="",
                public_key=b"",
                additional_recipients=_read_additional_recipients(
                    g_sig, g_first),
                server_kek_id=_read_server_kek_id(g_sig, g_first),
            )
            g_acquisition_mode = io.read_int_attr(g_run_group,
                                                    "acquisition_mode",
                                                    default=0) or 0
            # Genomic-run metadata JSON (convention).
            g_metadata_json = json.dumps({
                "modality": io.read_string_attr(g_run_group,
                                                  "modality") or "",
                "platform": io.read_string_attr(g_run_group,
                                                  "platform") or "",
                "reference_uri": io.read_string_attr(g_run_group,
                                                       "reference_uri") or "",
                "sample_name": io.read_string_attr(g_run_group,
                                                     "sample_name") or "",
                # M97: "" when the run has no @read_role attribute.
                "read_role": io.read_string_attr(g_run_group,
                                                   "read_role") or "",
            }, sort_keys=True)
            # Read the plaintext genomic_index columns (these are NOT
            # encrypted by M90.1 — only signal channels are).
            g_chromosomes = _read_chromosomes_compound(g_idx)
            g_positions = np.asarray(
                g_idx.open_dataset("positions").read(), dtype=np.int64,
            )
            g_mapqs = np.asarray(
                g_idx.open_dataset("mapping_qualities").read(),
                dtype=np.uint8,
            )
            g_flags = np.asarray(
                g_idx.open_dataset("flags").read(), dtype=np.uint32,
            )
            g_segments_by_name = {
                c: io.read_channel_segments(g_sig, f"{c}_segments")
                for c in g_channel_names
            }
            n_reads = (len(next(iter(g_segments_by_name.values())))
                        if g_segments_by_name else 0)

            writer.write_dataset_header(
                dataset_id=g_dataset_id,
                name=g_run_name,
                acquisition_mode=int(g_acquisition_mode),
                spectrum_class="TTIOGenomicRead",
                channel_names=list(g_channel_names),
                instrument_json=g_metadata_json,
                expected_au_count=n_reads,
            )

            if genomic_layouts.get(g_run_name) == "blocks_v1":
                _emit_blocks_v1_sidecars(writer, g_run_group, g_sig,
                                          dataset_id=g_dataset_id)

            for i in range(n_reads):
                # Build encrypted ChannelData list (UINT8 for genomic).
                g_channels = []
                for cname in g_channel_names:
                    seg = g_segments_by_name[cname][i]
                    data = bytes(seg.iv) + bytes(seg.tag) + bytes(seg.ciphertext)
                    g_channels.append(ChannelData(
                        name=cname,
                        precision=int(Precision.UINT8) & 0xFF,
                        compression=0,  # NONE (inner plaintext is raw uint8)
                        n_elements=int(seg.length),
                        data=data,
                    ))
                au = AccessUnit(
                    spectrum_class=5,
                    acquisition_mode=int(g_acquisition_mode),
                    ms_level=0,
                    polarity=2,
                    retention_time=0.0,
                    precursor_mz=0.0,
                    precursor_charge=0,
                    ion_mobility=0.0,
                    base_peak_intensity=0.0,
                    channels=g_channels,
                    chromosome=g_chromosomes[i],
                    position=int(g_positions[i]),
                    mapping_quality=int(g_mapqs[i]),
                    flags=int(g_flags[i]) & 0xFFFF,
                )
                _emit_raw_au(
                    writer,
                    dataset_id=g_dataset_id,
                    au_sequence=i,
                    payload=au.to_bytes(),
                    flags=int(PacketFlag.ENCRYPTED),
                )
            writer.write_end_of_dataset(
                dataset_id=g_dataset_id,
                final_au_sequence=n_reads,
            )

        writer.write_end_of_stream()
    finally:
        sp.close()


#: blocks_v1 channels whose coded blobs stay plaintext under per-AU
#: protection and ride in BlockSidecar packets as verbatim slices.
_BLOCKS_V1_SIDECAR_BLOBS = ("read_names", "cigars", "mate_info")


def _json_safe_attr(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray) and v.size == 1:
        return _json_safe_attr(v.item())
    return v


def _pack_json(obj) -> bytes:
    data = json.dumps(obj, sort_keys=True).encode("utf-8")
    return struct.pack("<I", len(data)) + data


def _unpack_json(payload: bytes, off: int):
    (n,) = struct.unpack_from("<I", payload, off)
    off += 4
    return json.loads(payload[off:off + n].decode("utf-8")), off + n


def _compound_names(group, name: str) -> list[str]:
    out = []
    for row in io.read_compound_dataset(group, name):
        v = row["name"]
        out.append(v.decode("utf-8") if isinstance(v, bytes) else v)
    return out


def _emit_blocks_v1_sidecars(writer, run_group, sig, *,
                              dataset_id: int) -> None:
    """One GenomicRunSidecar, then one BlockSidecar per block, for a
    blocks_v1 per-AU-encrypted genomic run (transport-spec §4.24).

    The run sidecar carries the run scalars (layout, block_policy,
    read_count, base_count), the optional restore attrs
    (ref_diff_slice_bytes, opt_disable_qualities_v5, reference_md5s)
    the receiver writes back verbatim, the sidecar channels' dataset
    attrs, and the chromosome_names and mate chrom_names tables in
    row order. Each block sidecar carries that block's index row and
    its verbatim slices of the plaintext channel blobs."""
    from ..genomic._block_view import BlockTable

    table = BlockTable.read(run_group)
    # The index's own channel set: a run without tags sends the five
    # required entries, exactly as before M101.
    BLOCK_CHANNELS = table.channels
    attrs: dict = {}
    for name in ("ref_diff_slice_bytes", "opt_disable_qualities_v5"):
        v = io.read_int_attr(run_group, name, default=0)
        if v:
            attrs[name] = int(v)
    md5s = io.read_string_attr(run_group, "reference_md5s")
    if md5s:
        attrs["reference_md5s"] = md5s

    channels = []
    blob_ds = {}
    for ch in _BLOCKS_V1_SIDECAR_BLOBS:
        ds = None
        if ch == "mate_info":
            if sig.has_child("mate_info"):
                g = sig.open_group("mate_info")
                if g.has_child("inline_v2"):
                    ds = g.open_dataset("inline_v2")
        elif sig.has_child(ch):
            ds = sig.open_dataset(ch)
        if ds is None:
            continue
        blob_ds[ch] = ds
        centry = {"name": ch,
                  "compression": (int(ds.get_attribute("compression"))
                                   if ds.has_attribute("compression")
                                   else 0)}
        extra = {k: _json_safe_attr(ds.get_attribute(k))
                 for k in ds.attribute_names() if k != "compression"}
        if extra:
            centry["extra_attrs"] = extra
        channels.append(centry)

    idx = run_group.open_group("genomic_index")
    chrom_names = _compound_names(idx, "chromosome_names")
    mate_names: list[str] = []
    if sig.has_child("mate_info"):
        mg = sig.open_group("mate_info")
        if mg.has_child("chrom_names"):
            mate_names = _compound_names(mg, "chrom_names")

    payload = (
        pack_string(io.read_string_attr(run_group, "layout") or "",
                    width=2)
        + pack_string(io.read_string_attr(run_group, "block_policy")
                      or "", width=2)
        + struct.pack("<QQ",
                       int(io.read_int_attr(run_group, "read_count",
                                              default=0) or 0),
                       int(io.read_int_attr(run_group, "base_count",
                                              default=0) or 0))
        + _pack_json(attrs)
        + _pack_json(channels)
        + struct.pack("<I", len(chrom_names))
        + b"".join(pack_string(n, width=2) for n in chrom_names)
        + struct.pack("<I", len(mate_names))
        + b"".join(pack_string(n, width=2) for n in mate_names)
    )
    writer._emit(PacketType.GENOMIC_RUN_SIDECAR, payload,
                  dataset_id=dataset_id)

    for b in range(table.count):
        parts = [struct.pack("<IQIQQ", b,
                              int(table.read_start[b]),
                              int(table.n_reads[b]),
                              int(table.base_start[b]),
                              int(table.n_bases[b])),
                 struct.pack("<B", len(BLOCK_CHANNELS))]
        for ch in BLOCK_CHANNELS:
            codec = int(table.codecs[ch][b]) if table.codecs else 0
            parts.append(pack_string(ch, width=2)
                         + struct.pack("<QQI",
                                        int(table.ranges[ch][0][b]),
                                        int(table.ranges[ch][1][b]),
                                        codec))
        parts.append(struct.pack("<B", len(blob_ds)))
        for ch, ds in blob_ds.items():
            off = int(table.ranges[ch][0][b])
            ln = int(table.ranges[ch][1][b])
            data = (bytes(np.asarray(ds.read(off, ln),
                                      dtype=np.uint8).tobytes())
                    if ln else b"")
            parts.append(pack_string(ch, width=2)
                         + struct.pack("<I", len(data)) + data)
        writer._emit(PacketType.BLOCK_SIDECAR, b"".join(parts),
                      dataset_id=dataset_id, au_sequence=b)


def _decode_genomic_run_sidecar(payload: bytes) -> dict:
    off = 0
    layout, off = unpack_string(payload, off, width=2)
    block_policy, off = unpack_string(payload, off, width=2)
    read_count, base_count = struct.unpack_from("<QQ", payload, off)
    off += 16
    attrs, off = _unpack_json(payload, off)
    channels, off = _unpack_json(payload, off)
    (n,) = struct.unpack_from("<I", payload, off); off += 4
    chrom_names = []
    for _ in range(n):
        s, off = unpack_string(payload, off, width=2)
        chrom_names.append(s)
    (n,) = struct.unpack_from("<I", payload, off); off += 4
    mate_chrom_names = []
    for _ in range(n):
        s, off = unpack_string(payload, off, width=2)
        mate_chrom_names.append(s)
    return {"layout": layout, "block_policy": block_policy,
            "read_count": int(read_count), "base_count": int(base_count),
            "attrs": attrs, "channels": channels,
            "chromosome_names": chrom_names,
            "mate_chrom_names": mate_chrom_names}


def _decode_block_sidecar(payload: bytes) -> dict:
    b, read_start, n_reads, base_start, n_bases = struct.unpack_from(
        "<IQIQQ", payload, 0)
    off = 32
    (nch,) = struct.unpack_from("<B", payload, off); off += 1
    channels: dict = {}
    for _ in range(nch):
        name, off = unpack_string(payload, off, width=2)
        c_off, c_len, codec = struct.unpack_from("<QQI", payload, off)
        off += 20
        channels[name] = (int(c_off), int(c_len), int(codec))
    (nb,) = struct.unpack_from("<B", payload, off); off += 1
    blobs: dict = {}
    for _ in range(nb):
        name, off = unpack_string(payload, off, width=2)
        (ln,) = struct.unpack_from("<I", payload, off); off += 4
        blobs[name] = bytes(payload[off:off + ln]); off += ln
    return {"block_index": int(b), "read_start": int(read_start),
            "n_reads": int(n_reads), "base_start": int(base_start),
            "n_bases": int(n_bases), "channels": channels,
            "blobs": blobs}


def _read_chromosomes_compound(idx_group) -> list[str]:
    """Read the genomic_index chromosome columns → list[str].

    L1 (Task #82 Phase B.1, 2026-05-01): chromosomes are now stored as
    a uint16 id column + compound name lookup table. Materialize back
    to ``list[str]`` for callers that still want the per-read view.
    """
    import numpy as np
    ids_ds = idx_group.open_dataset("chromosome_ids")
    ids = np.asarray(ids_ds.read(), dtype=np.uint16)
    name_rows = io.read_compound_dataset(idx_group, "chromosome_names")
    name_table: list[str] = []
    for row in name_rows:
        v = row["name"]
        name_table.append(v.decode("utf-8") if isinstance(v, bytes) else v)
    return [name_table[i] for i in ids.tolist()]


def _wire_polarity(raw: int) -> int:
    if raw == 1: return 0      # POSITIVE
    if raw == -1: return 1     # NEGATIVE
    return 2                   # UNKNOWN


def _storage_polarity(wire: int) -> int:
    """Inverse of :func:`_wire_polarity`: map the on-the-wire polarity
    encoding back to the stored ``spectrum_index/polarities`` form."""
    if wire == 0: return 1     # POSITIVE
    if wire == 1: return -1    # NEGATIVE
    return 0                   # UNKNOWN


def _encode_recipient_block(
    additional_recipients: "list[tuple[str, str, bytes]]",
) -> bytes:
    """Encode the FD-1 Phase A append-only trailing block:
    ``additional_recipient_count`` (u16) + per recipient
    ``{recipient_id, kek_algorithm, wrapped_dek}``. Empty input → empty
    bytes (the single-recipient case emits NO trailing block, keeping the
    packet byte-identical to transport-spec §4.4)."""
    if not additional_recipients:
        return b""
    out = struct.pack("<H", len(additional_recipients))
    for recipient_id, kek_algorithm, wrapped_dek in additional_recipients:
        out += (pack_string(recipient_id, width=2)
                + pack_string(kek_algorithm, width=2)
                + struct.pack("<I", len(wrapped_dek)) + wrapped_dek)
    return out


def _decode_recipient_block(payload: bytes, off: int):
    """Inverse of :func:`_encode_recipient_block`. Returns
    ``(list[(recipient_id, kek_algorithm, wrapped_dek)], new_off)``."""
    (n,) = struct.unpack_from("<H", payload, off)
    off += 2
    out = []
    for _ in range(n):
        recipient_id, off = unpack_string(payload, off, width=2)
        kek_algorithm, off = unpack_string(payload, off, width=2)
        (wl,) = struct.unpack_from("<I", payload, off)
        off += 4
        wrapped = bytes(payload[off:off + wl])
        off += wl
        out.append((recipient_id, kek_algorithm, wrapped))
    return out, off


def _read_additional_recipients(sig, first_channel: str):
    """Read additional (non-primary) DEK recipients stamped on a
    ``signal_channels`` group as the ``<channel>_wrapped_dek_recipients``
    attribute (an encoded recipient block). Returns ``[]`` when absent —
    i.e. single-recipient runs, which is the common case."""
    attr = f"{first_channel}_wrapped_dek_recipients"
    if not sig.has_attribute(attr):
        return []
    blob = bytes(sig.get_attribute(attr))
    if not blob:
        return []
    recips, _ = _decode_recipient_block(blob, 0)
    return recips


def _read_server_kek_id(sig, first_channel: str):
    """Read the FD-1 C-2a server kek_id stamped as
    ``<channel>_server_kek_id``. Returns ``None`` when absent (BYOK /
    not server-processable)."""
    return io.read_string_attr(sig, f"{first_channel}_server_kek_id") or None


def _stamp_additional_recipients_attr(sig, channel: str, pm: dict) -> None:
    """Persist a decoded packet's *additional* recipients (everything past
    the primary) onto ``<channel>_wrapped_dek_recipients`` as a uint8 array,
    and the C-2a ``server_kek_id`` onto ``<channel>_server_kek_id``. No-op
    for single-recipient BYOK packets, so MS/NMR/BYOK files keep only the
    existing ``<channel>_wrapped_dek`` attribute."""
    recipients = pm.get("recipients") or []
    additional = recipients[1:]
    if additional:
        sig.set_attribute(
            f"{channel}_wrapped_dek_recipients",
            np.frombuffer(_encode_recipient_block(additional), dtype=np.uint8))
    server_kek_id = pm.get("server_kek_id")
    if server_kek_id:
        sig.set_attribute(f"{channel}_server_kek_id", server_kek_id)


def _emit_protection_metadata(
    writer: TransportWriter,
    *,
    dataset_id: int,
    cipher_suite: str,
    kek_algorithm: str,
    wrapped_dek: bytes,
    signature_algorithm: str,
    public_key: bytes,
    additional_recipients: "list[tuple[str, str, bytes]]" = (),
    server_kek_id: "str | None" = None,
) -> None:
    payload = (
        pack_string(cipher_suite, width=2)
        + pack_string(kek_algorithm, width=2)
        + struct.pack("<I", len(wrapped_dek))
        + wrapped_dek
        + pack_string(signature_algorithm, width=2)
        + struct.pack("<I", len(public_key))
        + public_key
        + _encode_protection_trailing(additional_recipients, server_kek_id)
    )
    writer._emit(PacketType.PROTECTION_METADATA, payload,
                  dataset_id=dataset_id)


def _encode_protection_trailing(
    additional_recipients: "list[tuple[str, str, bytes]]",
    server_kek_id: "str | None",
) -> bytes:
    """FD-1 C-2a trailing section after the five transport-spec §4.4 fields.
    Emitted ONLY when there are additional recipients OR a ``server_kek_id``
    (so pure BYOK / single-recipient packets stay byte-identical to §4.4):

        additional_recipient_count u16
        <count> recipient entries
        [ server_kek_id  u16 len + UTF-8 ]   # iff present

    A single-recipient server-processable container emits ``count = 0``
    followed by ``server_kek_id`` (Phase A readers tolerate ``count = 0``)."""
    if not additional_recipients and not server_kek_id:
        return b""
    if additional_recipients:
        out = _encode_recipient_block(list(additional_recipients))
    else:
        out = struct.pack("<H", 0)
    if server_kek_id:
        out += pack_string(server_kek_id, width=2)
    return out


def _emit_raw_au(
    writer: TransportWriter,
    *,
    dataset_id: int,
    au_sequence: int,
    payload: bytes,
    flags: int,
) -> None:
    """Bypass the public write_access_unit so custom flag bits (including
    ENCRYPTED_HEADER) make it onto the packet header."""
    header = PacketHeader(
        packet_type=int(PacketType.ACCESS_UNIT),
        flags=flags,
        dataset_id=dataset_id,
        au_sequence=au_sequence,
        payload_length=len(payload),
        timestamp_ns=now_ns(),
    )
    writer._stream.write(header.to_bytes())
    writer._stream.write(payload)
    if writer._use_checksum:
        from .packets import crc32c
        writer._stream.write(struct.pack("<I", crc32c(payload)))


# ---------------------------------------------------------------------
# Reader side: encrypted stream → new .tio file
# ---------------------------------------------------------------------


def read_encrypted_to_file(
    stream_source,
    output_path: str | Path,
    *,
    provider: str | None = None,
) -> dict:
    """Materialise an encrypted transport stream into a new .tio
    file preserving the encrypted ChannelData bytes verbatim. Routes
    through :func:`ttio.providers.open_provider` so the output
    can live in any backend that supports VL_BYTES compound writes.

    ``stream_source`` may be a ``BinaryIO`` or a path to a ``.tis``
    file. The output file is written with ``opt_per_au_encryption``
    (and ``opt_encrypted_au_headers`` when the source stream used
    encrypted headers). The receiver does NOT decrypt in transit —
    the emitted file carries the same ciphertext the sender stored.

    Returns a metadata dict summarising the stream
    (``{"title", "runs": {run_name: {"n_spectra", "channels"}}}``)
    mainly for testing / introspection.
    """
    from .codec import TransportReader
    from ..encryption_per_au import ChannelSegment, HeaderSegment
    from ..providers.registry import open_provider

    # Accumulate stream → in-memory structure, then emit .tio at the end.
    stream_meta: dict = {}
    datasets: dict[int, dict] = {}  # dataset_id -> {name, channel_names, ...}
    protection: dict[int, dict] = {}

    reader = TransportReader(stream_source)
    for header, payload in reader.iter_packets():
        ptype = header.packet_type
        if ptype == int(PacketType.STREAM_HEADER):
            from .codec import _decode_stream_header
            stream_meta = _decode_stream_header(payload)
        elif ptype == int(PacketType.DATASET_HEADER):
            from .codec import _decode_dataset_header
            meta = _decode_dataset_header(payload)
            did = meta["dataset_id"]
            datasets[did] = {
                "meta": meta,
                "channel_segments": {c: [] for c in meta["channel_names"]},
                "header_segments": [],
                "used_encrypted_headers": False,
                # MS per-spectrum metadata carried in the plaintext AU
                # filter header (encrypt_headers=False path). Accumulated by
                # _ingest_encrypted_au and written back to spectrum_index so
                # the round-trip preserves it -- mirrors the Java/ObjC
                # readers (TTI-O#199).
                "plain_rts": [],
                "plain_ms_levels": [],
                "plain_polarities": [],
                "plain_precursor_mzs": [],
                "plain_precursor_charges": [],
                "plain_base_peaks": [],
                # genomic-only accumulator, populated by
                # _ingest_encrypted_au when spectrum_class == 5.
                "is_genomic": meta["spectrum_class"] == "TTIOGenomicRead",
                "genomic_chromosomes": [],
                "genomic_positions": [],
                "genomic_mapqs": [],
                "genomic_flags": [],
                # M99.1 blocks_v1 sidecars; run_sidecar stays None for
                # legacy-shaped genomic runs.
                "run_sidecar": None,
                "block_sidecars": [],
            }
        elif ptype == int(PacketType.PROTECTION_METADATA):
            pm = _decode_protection_metadata(payload)
            protection[header.dataset_id] = pm
        elif ptype == int(PacketType.GENOMIC_RUN_SIDECAR):
            did = header.dataset_id
            if did not in datasets:
                raise ValueError(
                    f"GenomicRunSidecar for unknown dataset_id {did}")
            datasets[did]["run_sidecar"] = _decode_genomic_run_sidecar(
                payload)
        elif ptype == int(PacketType.BLOCK_SIDECAR):
            did = header.dataset_id
            if did not in datasets:
                raise ValueError(
                    f"BlockSidecar for unknown dataset_id {did}")
            datasets[did]["block_sidecars"].append(
                _decode_block_sidecar(payload))
        elif ptype == int(PacketType.ACCESS_UNIT):
            did = header.dataset_id
            if did not in datasets:
                raise ValueError(f"AU for unknown dataset_id {did}")
            _ingest_encrypted_au(
                datasets[did],
                header=header,
                payload=payload,
                dataset_id=did,
                au_sequence=header.au_sequence,
            )
        elif ptype == int(PacketType.END_OF_STREAM):
            break
        # EndOfDataset / Annotation / Provenance / Chromatogram: skip for now.

    # Emit the .tio with encrypted compounds. The blocks_v1 wire
    # token is transport-scoped and never a container feature flag.
    features = set(stream_meta.get("features", []))
    features.discard(TRANSPORT_BLOCKS_V1_FEATURE)
    features.add(OPT_PER_AU_ENCRYPTION)
    any_encrypted_headers = any(d["used_encrypted_headers"] for d in datasets.values())
    if any_encrypted_headers:
        features.add(OPT_ENCRYPTED_AU_HEADERS)

    sp_out = open_provider(str(output_path), provider=provider, mode="w")
    try:
        root = sp_out.root_group()
        io.write_feature_flags(root, "1.1", sorted(features))
        study = root.create_group("study")
        io.write_fixed_string_attr(study, "title", stream_meta.get("title", ""))
        io.write_fixed_string_attr(study, "isa_investigation_id",
                                     stream_meta.get("isa_investigation", ""))
        # split datasets into MS vs genomic for output.
        ms_datasets = {did: d for did, d in datasets.items()
                        if not d.get("is_genomic")}
        genomic_datasets = {did: d for did, d in datasets.items()
                             if d.get("is_genomic")}

        ms_runs = study.create_group("ms_runs")
        ms_names = ",".join(d["meta"]["name"]
                              for _, d in sorted(ms_datasets.items()))
        io.write_fixed_string_attr(ms_runs, "_run_names", ms_names)

        for did, d in sorted(ms_datasets.items()):
            meta = d["meta"]
            run_group = ms_runs.create_group(meta["name"])
            io.write_int_attr(run_group, "acquisition_mode",
                                meta["acquisition_mode"])
            io.write_int_attr(run_group, "spectrum_count",
                                len(d["header_segments"])
                                if d["used_encrypted_headers"]
                                else len(next(iter(d["channel_segments"].values()))))
            io.write_fixed_string_attr(run_group, "spectrum_class",
                                         meta["spectrum_class"])
            cfg = run_group.create_group("instrument_config")
            for fname in ("manufacturer", "model", "serial_number",
                            "source_type", "analyzer_type", "detector_type"):
                io.write_fixed_string_attr(cfg, fname, "")

            sig = run_group.create_group("signal_channels")
            io.write_fixed_string_attr(sig, "channel_names",
                                         ",".join(meta["channel_names"]))
            for cname in meta["channel_names"]:
                segs = d["channel_segments"][cname]
                io.write_channel_segments(sig, f"{cname}_segments", segs)
                sig.set_attribute(f"{cname}_algorithm", "aes-256-gcm")
                pm = protection.get(did)
                if pm and pm.get("wrapped_dek"):
                    # Store the wrapped DEK as a byte array, not a string:
                    # the v1.2 / ML-KEM blob contains embedded NULs that a
                    # VLEN-string attribute rejects.
                    sig.set_attribute(
                        f"{cname}_wrapped_dek",
                        np.frombuffer(pm["wrapped_dek"], dtype=np.uint8))
                    sig.set_attribute(f"{cname}_kek_algorithm",
                                         pm["kek_algorithm"])
                    _stamp_additional_recipients_attr(sig, cname, pm)

            idx = run_group.create_group("spectrum_index")
            first_segs = next(iter(d["channel_segments"].values()))
            io.write_int_attr(idx, "count", len(first_segs))
            offsets_arr = np.array([s.offset for s in first_segs], dtype="<u8")
            lengths_arr = np.array([s.length for s in first_segs], dtype="<u4")
            ds_off = idx.create_dataset("offsets", Precision.INT64,
                                           len(first_segs))
            ds_off.write(offsets_arr)
            ds_len = idx.create_dataset("lengths", Precision.UINT32,
                                           len(first_segs))
            ds_len.write(lengths_arr)
            if d["used_encrypted_headers"]:
                io.write_au_header_segments(idx, "au_header_segments",
                                              d["header_segments"])
            elif d["plain_rts"]:
                # Restore the per-spectrum metadata carried in the plaintext
                # AU filter header (TTI-O#199), mirroring the Java/ObjC
                # readers so the encrypted round-trip is lossless.
                n_idx = len(first_segs)
                _idx_cols = (
                    ("retention_times", Precision.FLOAT64, d["plain_rts"], "<f8"),
                    ("ms_levels", Precision.INT32, d["plain_ms_levels"], "<i4"),
                    ("polarities", Precision.INT32, d["plain_polarities"], "<i4"),
                    ("precursor_mzs", Precision.FLOAT64,
                     d["plain_precursor_mzs"], "<f8"),
                    ("precursor_charges", Precision.INT32,
                     d["plain_precursor_charges"], "<i4"),
                    ("base_peak_intensities", Precision.FLOAT64,
                     d["plain_base_peaks"], "<f8"),
                )
                for _name, _prec, _vals, _dt in _idx_cols:
                    _ds = idx.create_dataset(_name, _prec, n_idx)
                    _ds.write(np.array(_vals, dtype=_dt))

        # write genomic_runs/ if any genomic datasets came
        # through the stream.
        if genomic_datasets:
            g_runs_group = study.create_group("genomic_runs")
            g_names = ",".join(d["meta"]["name"]
                                 for _, d in sorted(genomic_datasets.items()))
            io.write_fixed_string_attr(g_runs_group, "_run_names", g_names)
            for did, d in sorted(genomic_datasets.items()):
                if d.get("run_sidecar"):
                    _write_blocks_v1_genomic_run(
                        g_runs_group, d, protection.get(did))
                    continue
                meta = d["meta"]
                g_metadata = json.loads(meta.get("instrument_json") or "{}")
                g_run_group = g_runs_group.create_group(meta["name"])
                io.write_int_attr(g_run_group, "acquisition_mode",
                                    meta["acquisition_mode"])
                io.write_fixed_string_attr(g_run_group, "spectrum_class",
                                             meta["spectrum_class"])
                io.write_fixed_string_attr(g_run_group, "modality",
                                             g_metadata.get("modality", ""))
                io.write_fixed_string_attr(g_run_group, "platform",
                                             g_metadata.get("platform", ""))
                io.write_fixed_string_attr(g_run_group, "reference_uri",
                                             g_metadata.get("reference_uri", ""))
                io.write_fixed_string_attr(g_run_group, "sample_name",
                                             g_metadata.get("sample_name", ""))
                # M97: absent on the container when "" on the wire.
                if g_metadata.get("read_role"):
                    io.write_fixed_string_attr(g_run_group, "read_role",
                                                 g_metadata["read_role"])

                g_sig = g_run_group.create_group("signal_channels")
                io.write_fixed_string_attr(g_sig, "channel_names",
                                             ",".join(meta["channel_names"]))
                for cname in meta["channel_names"]:
                    segs = d["channel_segments"][cname]
                    io.write_channel_segments(g_sig, f"{cname}_segments", segs)
                    g_sig.set_attribute(f"{cname}_algorithm", "aes-256-gcm")
                    pm = protection.get(did)
                    if pm and pm.get("wrapped_dek"):
                        g_sig.set_attribute(
                            f"{cname}_wrapped_dek",
                            np.frombuffer(pm["wrapped_dek"], dtype=np.uint8))
                        g_sig.set_attribute(f"{cname}_kek_algorithm",
                                             pm["kek_algorithm"])
                        _stamp_additional_recipients_attr(g_sig, cname, pm)

                g_idx = g_run_group.create_group("genomic_index")
                first_segs = next(iter(d["channel_segments"].values()))
                n_reads = len(first_segs)
                io.write_int_attr(g_idx, "count", n_reads)
                offsets_arr = np.array(
                    [s.offset for s in first_segs], dtype="<u8",
                )
                lengths_arr = np.array(
                    [s.length for s in first_segs], dtype="<u4",
                )
                ds_off = g_idx.create_dataset(
                    "offsets", Precision.UINT64, n_reads,
                )
                ds_off.write(offsets_arr)
                ds_len = g_idx.create_dataset(
                    "lengths", Precision.UINT32, n_reads,
                )
                ds_len.write(lengths_arr)
                # Per-read genomic suffix columns from the AU stream.
                ds_pos = g_idx.create_dataset(
                    "positions", Precision.INT64, n_reads,
                )
                ds_pos.write(np.array(d["genomic_positions"], dtype=np.int64))
                ds_mq = g_idx.create_dataset(
                    "mapping_qualities", Precision.UINT8, n_reads,
                )
                ds_mq.write(np.array(d["genomic_mapqs"], dtype=np.uint8))
                ds_fl = g_idx.create_dataset(
                    "flags", Precision.UINT32, n_reads,
                )
                ds_fl.write(np.array(d["genomic_flags"], dtype=np.uint32))
                # L1 (Task #82 Phase B.1): write chromosomes as
                # uint16 id column + compound name lookup table —
                # the M82-era VL-string compound cost 42 MB of
                # fractal-heap overhead per chr22 file.
                _chroms = d["genomic_chromosomes"]
                _name_to_id: dict[str, int] = {}
                _names: list[str] = []
                _ids = np.empty(len(_chroms), dtype=np.uint16)
                for _i, _name in enumerate(_chroms):
                    _slot = _name_to_id.get(_name)
                    if _slot is None:
                        _slot = len(_names)
                        _name_to_id[_name] = _slot
                        _names.append(_name)
                    _ids[_i] = _slot
                ds_cids = g_idx.create_dataset(
                    "chromosome_ids", Precision.UINT16, n_reads,
                )
                ds_cids.write(_ids)
                io.write_compound_dataset(
                    g_idx,
                    "chromosome_names",
                    [{"name": n} for n in _names],
                    [("name", io.vl_str())],
                )
    finally:
        sp_out.close()

    return {
        "title": stream_meta.get("title", ""),
        "runs": {d["meta"]["name"]: {
            "n_spectra": len(d["header_segments"])
                           if d["used_encrypted_headers"]
                           else len(next(iter(d["channel_segments"].values()))),
            "channels": list(d["meta"]["channel_names"]),
            "encrypted_headers": d["used_encrypted_headers"],
        } for d in datasets.values()},
    }


def _decode_protection_metadata(payload: bytes) -> dict:
    off = 0
    cipher_suite, off = unpack_string(payload, off, width=2)
    kek_algorithm, off = unpack_string(payload, off, width=2)
    (wrapped_len,) = struct.unpack_from("<I", payload, off); off += 4
    wrapped_dek = bytes(payload[off:off + wrapped_len]); off += wrapped_len
    signature_algorithm, off = unpack_string(payload, off, width=2)
    (pk_len,) = struct.unpack_from("<I", payload, off); off += 4
    public_key = bytes(payload[off:off + pk_len]); off += pk_len
    # FD-1 Phase A: the primary recipient is the in-band wrapped DEK;
    # additional recipients (if any) follow in the trailing block. Older
    # single-recipient packets have no trailing bytes -> one recipient.
    recipients = [("", kek_algorithm, wrapped_dek)]
    server_kek_id = None
    if off < len(payload):
        extra, off = _decode_recipient_block(payload, off)
        recipients.extend(extra)
        # FD-1 C-2a: anything after the recipient block is the optional
        # server_kek_id (names the primary recipient's KEK).
        if off < len(payload):
            server_kek_id, off = unpack_string(payload, off, width=2)
    return {
        "cipher_suite": cipher_suite,
        "kek_algorithm": kek_algorithm,
        "wrapped_dek": wrapped_dek,
        "signature_algorithm": signature_algorithm,
        "public_key": public_key,
        "recipients": recipients,
        "server_kek_id": server_kek_id,
    }


def _write_blocks_v1_genomic_run(g_runs_group, d: dict,
                                   pm: dict | None) -> None:
    """Materialise one blocks_v1 per-AU-encrypted genomic run from
    its sidecar packets and AU stream, in the shape the stream
    writer creates it, so decrypt-in-place restores it."""
    from ..enums import Compression
    from ..genomic.stream_writer import (CHANNEL_CHUNK, _INDEX_ARRAYS,
                                           index_fields_for)

    meta = d["meta"]
    sc = d["run_sidecar"]
    g_metadata = json.loads(meta.get("instrument_json") or "{}")
    rg = g_runs_group.create_group(meta["name"])
    io.write_int_attr(rg, "acquisition_mode", meta["acquisition_mode"])
    io.write_fixed_string_attr(rg, "modality",
                                 g_metadata.get("modality", ""))
    io.write_int_attr(rg, "spectrum_class", 5)
    io.write_fixed_string_attr(rg, "reference_uri",
                                 g_metadata.get("reference_uri", ""))
    io.write_fixed_string_attr(rg, "platform",
                                 g_metadata.get("platform", ""))
    if g_metadata.get("read_role"):
        io.write_fixed_string_attr(rg, "read_role",
                                     g_metadata["read_role"])
    io.write_fixed_string_attr(rg, "sample_name",
                                 g_metadata.get("sample_name", ""))
    io.write_int_attr(rg, "read_count", sc["read_count"])
    io.write_int_attr(rg, "base_count", sc["base_count"])
    io.write_fixed_string_attr(rg, "layout", sc["layout"])
    io.write_fixed_string_attr(rg, "block_policy", sc["block_policy"])
    attrs = sc.get("attrs") or {}
    if attrs.get("ref_diff_slice_bytes"):
        io.write_int_attr(rg, "ref_diff_slice_bytes",
                            int(attrs["ref_diff_slice_bytes"]),
                            dtype="<u8")
    if attrs.get("opt_disable_qualities_v5"):
        io.write_int_attr(rg, "opt_disable_qualities_v5", 1)
    if attrs.get("reference_md5s"):
        io.write_fixed_string_attr(rg, "reference_md5s",
                                     attrs["reference_md5s"])

    sidecars = sorted(d["block_sidecars"],
                       key=lambda x: x["block_index"])
    blocks = rg.create_group("blocks")
    rows = []
    for bs in sidecars:
        row = {"read_start": bs["read_start"],
               "n_reads": bs["n_reads"],
               "base_start": bs["base_start"],
               "n_bases": bs["n_bases"]}
        for ch, (c_off, c_len, codec) in bs["channels"].items():
            row[f"{ch}_off"] = c_off
            row[f"{ch}_len"] = c_len
            row[f"{ch}_codec"] = codec
        rows.append(row)
    # The index carries the tags triple when the sender's did (M101).
    idx_ds = blocks.create_compound_dataset(
        "index", index_fields_for(rows), 0, extendable=True, chunk_rows=1024)
    if rows:
        idx_ds.append(rows)

    idx_group = rg.create_group("genomic_index")
    first_segs = next(iter(d["channel_segments"].values()))
    chrom_table = sc["chromosome_names"]
    name_to_id = {n: i for i, n in enumerate(chrom_table)}
    arrays = {
        "lengths": np.array([s.length for s in first_segs],
                             dtype=np.uint32),
        "positions": np.array(d["genomic_positions"], dtype=np.int64),
        "mapping_qualities": np.array(d["genomic_mapqs"],
                                       dtype=np.uint8),
        "flags": np.array(d["genomic_flags"], dtype=np.uint32),
        "chromosome_ids": np.array(
            [name_to_id[c] for c in d["genomic_chromosomes"]],
            dtype=np.uint16),
    }
    for name, prec, dt in _INDEX_ARRAYS:
        ds = idx_group.create_dataset(
            name, prec, 0, chunk_size=io.DEFAULT_SIGNAL_CHUNK,
            compression=Compression.ZLIB, compression_level=6,
            extendable=True)
        ds.append(arrays[name].astype(dt))
    io.write_compound_dataset(
        idx_group, "chromosome_names",
        [{"name": n} for n in chrom_table], [("name", io.vl_str())])

    sig = rg.create_group("signal_channels")
    for cname in meta["channel_names"]:
        segs = d["channel_segments"][cname]
        io.write_channel_segments(sig, f"{cname}_segments", segs)
        sig.set_attribute(f"{cname}_algorithm", "aes-256-gcm")
        if pm and pm.get("wrapped_dek"):
            sig.set_attribute(
                f"{cname}_wrapped_dek",
                np.frombuffer(pm["wrapped_dek"], dtype=np.uint8))
            sig.set_attribute(f"{cname}_kek_algorithm",
                                 pm["kek_algorithm"])
            _stamp_additional_recipients_attr(sig, cname, pm)

    mate_group = None
    for centry in sc.get("channels") or []:
        ch = centry["name"]
        if ch == "mate_info":
            mate_group = sig.create_group("mate_info")
            parent, ds_name = mate_group, "inline_v2"
        else:
            parent, ds_name = sig, ch
        codec = int(centry.get("compression", 0))
        ds = parent.create_dataset(
            ds_name, Precision.UINT8, 0, chunk_size=CHANNEL_CHUNK,
            compression=(Compression.ZLIB if codec == 0
                         else Compression.NONE),
            compression_level=6, extendable=True)
        io.write_int_attr(ds, "compression", codec, dtype="<u1")
        for k, v in (centry.get("extra_attrs") or {}).items():
            ds.set_attribute(k, v)
        for bs in sidecars:
            data = bs["blobs"].get(ch)
            if data:
                ds.append(np.frombuffer(data, dtype=np.uint8))
    # The stream writer creates mate_info and its chrom_names table
    # at close even when no mate blob was written; mirror that.
    if sc.get("mate_chrom_names"):
        if mate_group is None:
            mate_group = sig.create_group("mate_info")
        io.write_compound_dataset(
            mate_group, "chrom_names",
            [{"name": n} for n in sc["mate_chrom_names"]],
            [("name", io.vl_str())])


def _ingest_encrypted_au(d: dict, *, header, payload: bytes,
                           dataset_id: int, au_sequence: int) -> None:
    from ..encryption_per_au import ChannelSegment, HeaderSegment

    flags = header.flags
    encrypted_header = bool(flags & int(PacketFlag.ENCRYPTED_HEADER))
    encrypted_channel = bool(flags & int(PacketFlag.ENCRYPTED))
    if not encrypted_channel:
        raise ValueError("encrypted-transport reader saw plaintext AU")

    d["used_encrypted_headers"] = encrypted_header

    if encrypted_header:
        # Wire: spectrum_class(u8) n_channels(u8) IV(12) TAG(16)
        # ciphertext(36) [channels].
        off = 0
        spectrum_class = payload[off]; off += 1
        n_channels = payload[off]; off += 1
        hdr_iv = bytes(payload[off:off + 12]); off += 12
        hdr_tag = bytes(payload[off:off + 16]); off += 16
        hdr_ct = bytes(payload[off:off + 36]); off += 36
        d["header_segments"].append(HeaderSegment(
            iv=hdr_iv, tag=hdr_tag, ciphertext=hdr_ct,
        ))
        remaining = payload[off:]
    else:
        # Plaintext filter header: the encrypted-channel-only AU
        # variant. We don't decode the filter header; we do need to
        # skip it to reach the channels. Delegate to AccessUnit.
        au = AccessUnit.from_bytes(payload)
        # capture genomic suffix when this is a genomic AU.
        if d.get("is_genomic"):
            d["genomic_chromosomes"].append(au.chromosome)
            d["genomic_positions"].append(int(au.position))
            d["genomic_mapqs"].append(int(au.mapping_quality))
            d["genomic_flags"].append(int(au.flags) & 0xFFFFFFFF)
        else:
            # MS per-spectrum metadata rides in the plaintext filter
            # header; accumulate it so spectrum_index round-trips
            # (TTI-O#199). Reverse the wire polarity encoding back to
            # storage form -- the inverse of _wire_polarity.
            d["plain_rts"].append(float(au.retention_time))
            d["plain_ms_levels"].append(int(au.ms_level))
            d["plain_polarities"].append(_storage_polarity(int(au.polarity)))
            d["plain_precursor_mzs"].append(float(au.precursor_mz))
            d["plain_precursor_charges"].append(int(au.precursor_charge))
            d["plain_base_peaks"].append(float(au.base_peak_intensity))
        # AccessUnit.from_bytes already parsed channels; retrieve
        # them directly.
        for cname in list(d["channel_segments"].keys()):
            ch = next((c for c in au.channels if c.name == cname), None)
            if ch is None:
                continue
            if len(ch.data) < 28:
                raise ValueError(
                    f"encrypted channel {cname!r} data shorter than IV+TAG"
                )
            iv = ch.data[:12]
            tag = ch.data[12:28]
            ciphertext = ch.data[28:]
            # Offset is derived from the sum of previous lengths. We
            # don't have it in the payload, so reconstruct.
            prior = sum(s.length for s in d["channel_segments"][cname])
            d["channel_segments"][cname].append(ChannelSegment(
                offset=prior,
                length=ch.n_elements,
                iv=iv, tag=tag, ciphertext=ciphertext,
            ))
        return

    # After reading the encrypted header (encrypted_header branch),
    # walk the channel entries from ``remaining``.
    buf = remaining
    off = 0
    for cname in list(d["channel_segments"].keys())[:n_channels]:
        (name_len,) = struct.unpack_from("<H", buf, off); off += 2
        name = bytes(buf[off:off + name_len]).decode("utf-8"); off += name_len
        precision = buf[off]; off += 1
        compression = buf[off]; off += 1
        (n_elements,) = struct.unpack_from("<I", buf, off); off += 4
        (data_len,) = struct.unpack_from("<I", buf, off); off += 4
        data = bytes(buf[off:off + data_len]); off += data_len
        if len(data) < 28:
            raise ValueError(
                f"encrypted channel {name!r} data shorter than IV+TAG"
            )
        iv = data[:12]
        tag = data[12:28]
        ciphertext = data[28:]
        prior = sum(s.length for s in d["channel_segments"][name])
        d["channel_segments"][name].append(ChannelSegment(
            offset=prior,
            length=n_elements,
            iv=iv, tag=tag, ciphertext=ciphertext,
        ))


def stamp_transport_wrapped_dek(path, wrapped_dek, kek_algorithm, *,
                                  additional_recipients=(),
                                  server_kek_id=None, provider=None):
    """Stamp a wrapped DEK + KEK algorithm onto every run's
    ``signal_channels`` so :func:`write_encrypted_dataset` emits it in the
    ProtectionMetadata packet.

    Used for envelope / PQC per-AU upload, where the DEK is wrapped under
    the recipient's KEK (``aes-256-gcm``) or KEM public key
    (``ml-kem-1024``). The DEK itself is shared across all channels in a
    run (v1.0 design), so the wrapper is stamped on each run's first
    channel -- the same attribute :func:`write_encrypted_dataset` reads.

    ``additional_recipients`` (FD-1 Phase A) is an optional list of
    ``(recipient_id, kek_algorithm, wrapped_dek)`` for the *same* DEK
    wrapped under other recipients' keys; it is stamped as the
    ``<channel>_wrapped_dek_recipients`` attribute and emitted as the
    packet's trailing block. Empty (the default) keeps single-recipient
    runs byte-identical to pre-Phase-A.

    ``server_kek_id`` (FD-1 C-2a) is an optional opaque label naming the KEK
    under which the primary ``wrapped_dek`` is wrapped; it is stamped as
    ``<channel>_server_kek_id`` and emitted in the packet so the daemon can
    decide server-processability. ``None`` (the default) marks the container
    BYOK / not server-processable.
    """
    from ..providers.registry import open_provider

    wd = np.frombuffer(wrapped_dek, dtype=np.uint8)
    extra_blob = (np.frombuffer(
        _encode_recipient_block(list(additional_recipients)), dtype=np.uint8)
        if additional_recipients else None)
    sp = open_provider(path, provider=provider, mode="a")
    try:
        root = sp.root_group()
        if not root.has_child("study"):
            return
        study = root.open_group("study")
        for parent in ("ms_runs", "genomic_runs"):
            if not study.has_child(parent):
                continue
            g = study.open_group(parent)
            for n in g.child_names():
                if n.startswith("_") or not g.has_child(n):
                    continue
                run = g.open_group(n)
                if not run.has_child("signal_channels"):
                    continue
                sig = run.open_group("signal_channels")
                names = [c for c in (io.read_string_attr(sig, "channel_names")
                                     or "").split(",") if c]
                if not names:
                    continue
                sig.set_attribute(f"{names[0]}_wrapped_dek", wd)
                sig.set_attribute(f"{names[0]}_kek_algorithm", kek_algorithm)
                if extra_blob is not None:
                    sig.set_attribute(
                        f"{names[0]}_wrapped_dek_recipients", extra_blob)
                if server_kek_id:
                    sig.set_attribute(
                        f"{names[0]}_server_kek_id", server_kek_id)
    finally:
        sp.close()


def read_transport_wrapped_dek(path, *, provider=None):
    """Return ``(wrapped_dek_bytes, kek_algorithm)`` stamped on the first
    run's ``signal_channels``, or ``(b"", "")`` if none. The receiver
    unwraps the DEK with its KEK / KEM private key, then calls
    :func:`ttio.encryption_per_au.decrypt_per_au`."""
    from ..providers.registry import open_provider

    sp = open_provider(path, provider=provider, mode="r")
    try:
        root = sp.root_group()
        if not root.has_child("study"):
            return b"", ""
        study = root.open_group("study")
        for parent in ("ms_runs", "genomic_runs"):
            if not study.has_child(parent):
                continue
            g = study.open_group(parent)
            for n in g.child_names():
                if n.startswith("_") or not g.has_child(n):
                    continue
                run = g.open_group(n)
                if not run.has_child("signal_channels"):
                    continue
                sig = run.open_group("signal_channels")
                names = [c for c in (io.read_string_attr(sig, "channel_names")
                                     or "").split(",") if c]
                if not names:
                    continue
                attr = f"{names[0]}_wrapped_dek"
                if sig.has_attribute(attr):
                    return (bytes(sig.get_attribute(attr)),
                            io.read_string_attr(sig, f"{names[0]}_kek_algorithm")
                            or "")
        return b"", ""
    finally:
        sp.close()


def read_transport_recipients(path, *, provider=None):
    """Return the full recipient list stamped on the first run's
    ``signal_channels`` as ``[(recipient_id, kek_algorithm, wrapped_dek)]``
    (FD-1 Phase A). The primary recipient is index 0 with id ``""``;
    additional recipients follow. Empty list if no wrapped DEK is present.

    The single-recipient counterpart :func:`read_transport_wrapped_dek`
    stays the convenience accessor for the BYOK / envelope / PQC client
    paths that hold one key."""
    from ..providers.registry import open_provider

    sp = open_provider(path, provider=provider, mode="r")
    try:
        root = sp.root_group()
        if not root.has_child("study"):
            return []
        study = root.open_group("study")
        for parent in ("ms_runs", "genomic_runs"):
            if not study.has_child(parent):
                continue
            g = study.open_group(parent)
            for n in g.child_names():
                if n.startswith("_") or not g.has_child(n):
                    continue
                run = g.open_group(n)
                if not run.has_child("signal_channels"):
                    continue
                sig = run.open_group("signal_channels")
                names = [c for c in (io.read_string_attr(sig, "channel_names")
                                     or "").split(",") if c]
                if not names or not sig.has_attribute(f"{names[0]}_wrapped_dek"):
                    continue
                primary = (
                    "",
                    io.read_string_attr(sig, f"{names[0]}_kek_algorithm") or "",
                    bytes(sig.get_attribute(f"{names[0]}_wrapped_dek")),
                )
                return [primary] + _read_additional_recipients(sig, names[0])
        return []
    finally:
        sp.close()


def read_transport_server_kek_id(path, *, provider=None):
    """Return the FD-1 C-2a ``server_kek_id`` stamped on the first run's
    ``signal_channels`` (the KEK the daemon resolves to process the
    container server-side), or ``None`` if absent (BYOK / not
    server-processable)."""
    from ..providers.registry import open_provider

    sp = open_provider(path, provider=provider, mode="r")
    try:
        root = sp.root_group()
        if not root.has_child("study"):
            return None
        study = root.open_group("study")
        for parent in ("ms_runs", "genomic_runs"):
            if not study.has_child(parent):
                continue
            g = study.open_group(parent)
            for n in g.child_names():
                if n.startswith("_") or not g.has_child(n):
                    continue
                run = g.open_group(n)
                if not run.has_child("signal_channels"):
                    continue
                sig = run.open_group("signal_channels")
                names = [c for c in (io.read_string_attr(sig, "channel_names")
                                     or "").split(",") if c]
                if not names:
                    continue
                return _read_server_kek_id(sig, names[0])
        return None
    finally:
        sp.close()
