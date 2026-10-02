"""imzML + .ibd importer ().

imzML is the dominant interchange format for mass-spectrometry imaging
data. The format is a pair of files:

* ``.imzML`` — XML metadata that mirrors mzML's structure but adds a
  ``<scanSettings>`` block carrying the pixel grid (max count of x, y,
  z) and per-spectrum coordinates. Each ``<spectrum>`` element points
  at byte offsets in a sibling binary file.
* ``.ibd``   — concatenated binary mass / intensity arrays plus a
  16-byte UUID header that must match the UUID encoded in the .imzML.

Two storage modes are supported by the spec:

* **Continuous** — a single shared m/z axis stored once at the start of
  the .ibd; per-pixel intensity arrays follow.
* **Processed** — every pixel carries its own m/z + intensity arrays.

This importer covers both. It produces one spectrum per pixel; the
spatial grid (x, y, optionally z) and the ``IMS:1000030`` continuous /
``IMS:1000031`` processed designation are preserved as a per-run
:class:`ProvenanceRecord` parameter dict, from which an MSImage cube
writer can allocate a 3-D ``[height, width, spectral_points]`` dataset.

Cross-language equivalents
--------------------------
Objective-C: ``TTIOImzMLReader`` (Import/)
Java:        ``ImzMLReader.java`` (importers/)

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import iterparse

import numpy as np

from ..io.progress import ProgressSinkLike, _fire
from ..provenance import ProvenanceRecord
from ..spectral_dataset import SpectralDataset, WrittenRun


#: Mirror Java's ``ImzMLReader.PROGRESS_INTERVAL_SPECTRA``.
PROGRESS_INTERVAL_SPECTRA = 100


class ImzMLParseError(ValueError):
    """Raised when the imzML XML is structurally invalid for our needs."""


class ImzMLBinaryError(ValueError):
    """Raised when the .ibd binary disagrees with the .imzML metadata.

    Covers UUID mismatches and offset / length values that would read
    past the end of the .ibd -49).
    """


# imzML storage mode accessions. Only the IMS-namespaced forms are
# used — the MS-namespaced forms at the same numeric suffix are
# UNRELATED terms (MS:1000030 = "vendor processing software",
# MS:1000031 = "instrument model") and must not be matched here, or
# they collide with every MS:1000031 cvParam emitted elsewhere in the
# file.
_CONTINUOUS_ACCESSIONS = frozenset({"IMS:1000030"})
_PROCESSED_ACCESSIONS = frozenset({"IMS:1000031"})


@dataclass(slots=True)
class ImzMLPixelSpectrum:
    """One pixel's parsed spectrum + spatial coordinates."""
    x: int
    y: int
    z: int
    mz: np.ndarray   # 1-D float64
    intensity: np.ndarray  # 1-D float64 (same length as mz)


@dataclass(slots=True)
class ImzMLImport:
    """Result of parsing an imzML + .ibd pair.

    Attributes
    ----------
    mode
        Either ``"continuous"`` or ``"processed"``.
    uuid_hex
        16-byte UUID encoded as 32 lowercase hex characters. Cross-
        validated between the .imzML and the .ibd header — a mismatch
        is a hard error, not a warning.
    grid_max_x / grid_max_y / grid_max_z
        Pixel grid extents declared in the .imzML. Used by downstream
        MSImage writers to allocate the cube.
    pixel_size_x / pixel_size_y
        Spatial pixel size in micrometres (or whatever unit the
        producer claimed; we pass the value through unchanged).
    scan_pattern
        Free-text label, e.g. ``"flyback"`` / ``"meandering"``.
    spectra
        One :class:`ImzMLPixelSpectrum` per ``<spectrum>`` element.
    source_imzml / source_ibd
        Resolved Path objects for traceability.
    """
    mode: str
    uuid_hex: str
    grid_max_x: int
    grid_max_y: int
    grid_max_z: int
    pixel_size_x: float
    pixel_size_y: float
    scan_pattern: str
    spectra: list[ImzMLPixelSpectrum] = field(default_factory=list)
    source_imzml: str = ""
    source_ibd: str = ""

    # ---------------------------------------------------------------- #
    # Persistence — one spectrum per pixel, spatial metadata via
    # provenance.
    # ---------------------------------------------------------------- #

    def to_imported_dataset(
        self,
        *,
        title: str | None = None,
        isa_investigation_id: str = "",
        path: str | Path | None = None,
    ) -> "ImportedDataset":
        """Build the normalized :class:`ImportedDataset` draft.

        ``path`` is recorded in the provenance ``output_refs`` exactly as
        the legacy :meth:`to_ttio` did; pass the eventual destination so
        the round-trip provenance is byte-identical.
        """
        n = len(self.spectra)
        if n == 0:
            raise ImzMLParseError(f"{self.source_imzml}: no spectra parsed")

        lengths = np.array([s.mz.size for s in self.spectra], dtype=np.uint32)
        offsets = np.zeros(n, dtype=np.uint64)
        if n > 0:
            offsets[1:] = np.cumsum(lengths[:-1], dtype=np.uint64)
        total = int(lengths.sum())
        mz_buf = np.empty(total, dtype=np.float64)
        it_buf = np.empty(total, dtype=np.float64)
        pos = 0
        for s, length in zip(self.spectra, lengths):
            ln = int(length)
            mz_buf[pos:pos + ln] = s.mz
            it_buf[pos:pos + ln] = s.intensity
            pos += ln

        base_peaks = np.array(
            [float(np.max(s.intensity)) if s.intensity.size else 0.0 for s in self.spectra],
            dtype=np.float64,
        )

        # Per-spectrum spatial coordinates fit in retention_times only
        # awkwardly. Encode them in run-level provenance instead so the
        # spectrum_index keeps its conventional shape; downstream tools
        # can recover the grid from the provenance record.
        coords = [(s.x, s.y, s.z) for s in self.spectra]
        prov = ProvenanceRecord(
            timestamp_unix=int(time.time()),
            software="ttio imzml importer v0.9",
            parameters={
                "imzml_mode": self.mode,
                "imzml_uuid_hex": self.uuid_hex,
                "imzml_grid_max_x": int(self.grid_max_x),
                "imzml_grid_max_y": int(self.grid_max_y),
                "imzml_grid_max_z": int(self.grid_max_z),
                "imzml_pixel_size_x": float(self.pixel_size_x),
                "imzml_pixel_size_y": float(self.pixel_size_y),
                "imzml_scan_pattern": self.scan_pattern,
                "imzml_pixel_coordinates_csv": ";".join(f"{x},{y},{z}" for x, y, z in coords),
            },
            input_refs=[self.source_imzml, self.source_ibd],
            output_refs=[str(path)],
        )

        run = WrittenRun(
            spectrum_class="TTIOMassSpectrum",
            acquisition_mode=0,
            channel_data={"mz": mz_buf, "intensity": it_buf},
            offsets=offsets,
            lengths=lengths,
            retention_times=np.zeros(n, dtype=np.float64),
            ms_levels=np.ones(n, dtype=np.int32),
            polarities=np.zeros(n, dtype=np.int32),
            precursor_mzs=np.zeros(n, dtype=np.float64),
            precursor_charges=np.zeros(n, dtype=np.int32),
            base_peak_intensities=base_peaks,
            provenance_records=[prov],
        )
        from .imported_dataset import ImportedDataset
        return ImportedDataset(
            title=title or f"imzML import: {Path(self.source_imzml).name}",
            isa_investigation_id=isa_investigation_id,
            runs={"imzml_pixels": run},
            provenance=[prov],
        )

    def to_ttio(
        self,
        path: str | Path,
        *,
        title: str | None = None,
        isa_investigation_id: str = "",
    ) -> Path:
        """Write the imported pixels to an .tio as a single MS run."""
        return self.to_imported_dataset(
            title=title,
            isa_investigation_id=isa_investigation_id,
            path=path,
        ).write(path)


# --------------------------------------------------------------------------- #
# Parsing — XML with iterparse, .ibd as random-access reads.
# --------------------------------------------------------------------------- #

def _local(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


@dataclass(slots=True)
class _SpectrumStub:
    """In-progress spectrum collected while iterating the imzML tree."""
    x: int = 0
    y: int = 0
    z: int = 1
    mz_offset: int = -1
    mz_length: int = 0
    mz_encoded_length: int = 0
    int_offset: int = -1
    int_length: int = 0
    int_encoded_length: int = 0
    mz_precision: str = "64"
    int_precision: str = "64"


def read(
    imzml_path: str | Path,
    ibd_path: str | Path | None = None,
    *,
    progress: "ProgressSinkLike | None" = None,
) -> ImzMLImport:
    """Parse an imzML + .ibd pair and return an :class:`ImzMLImport`.

    Parameters
    ----------
    imzml_path
        Path to the ``.imzML`` XML metadata file.
    ibd_path
        Optional explicit path to the binary file. When ``None``, the
        sibling ``<stem>.ibd`` is used. UUID matching is still
        enforced regardless of how the .ibd was located.
    progress
        Optional :class:`~ttio.io.progress.ProgressSink` (or bare
        callable). Fires every :data:`PROGRESS_INTERVAL_SPECTRA` pixels
        during .ibd materialisation with ``total = -1`` and once more
        at the end with ``(n_pixels, n_pixels)``. The hook is inside
        the materialise loop, not at the boundary.
    """
    imzml = Path(imzml_path)
    if not imzml.is_file():
        raise FileNotFoundError(f"imzML metadata not found: {imzml}")
    ibd = Path(ibd_path) if ibd_path is not None else imzml.with_suffix(".ibd")
    if not ibd.is_file():
        raise FileNotFoundError(f"imzML binary not found: {ibd}")

    state_uuid = ""
    state_mode = ""
    grid_max = [0, 0, 1]
    pixel_size = [0.0, 0.0]
    scan_pattern = ""
    spectra_stubs: list[_SpectrumStub] = []
    current: _SpectrumStub | None = None
    in_spectrum = False
    in_binary_array = False
    in_position = False
    in_scan_settings = False  # disambiguate IMS:1000042 by context
    array_kind = ""  # "mz" or "intensity"
    # <referenceableParamGroup id=...> cvParams, applied wherever a
    # <referenceableParamGroupRef ref=...> appears. Writers commonly
    # declare the array kind / precision / compression of every
    # binaryDataArray this way (e.g. the HR2MSI PXD001283 files).
    param_groups: dict[str, list[tuple[str, str]]] = {}
    in_group: str | None = None

    def handle_cv(accession: str, value: str) -> None:
        nonlocal state_mode, state_uuid, scan_pattern, array_kind
        if accession in _CONTINUOUS_ACCESSIONS:
            state_mode = "continuous"
        elif accession in _PROCESSED_ACCESSIONS:
            state_mode = "processed"
        elif accession == "IMS:1000080" and value:  # canonical UUID accession
            state_uuid = _normalise_uuid(value)
        elif accession == "IMS:1000042" and in_scan_settings and value:
            # canonical max count of pixels x (inside scanSettings)
            grid_max[0] = int(value)
        elif accession == "IMS:1000043" and value:  # canonical max count of pixels y
            grid_max[1] = int(value)
        elif accession == "IMS:1000003" and value:  # legacy UUID-adjacent accession
            # Seen in the wild + our pre-0.9 synthetic fixtures as
            # a placeholder for max count of pixels x. Accept for
            # round-trip with older TTIO outputs.
            grid_max[0] = int(value)
        elif accession == "IMS:1000004" and value:  # legacy max count of pixels y
            grid_max[1] = int(value)
        elif accession == "IMS:1000005" and value:  # legacy max count of pixels z
            grid_max[2] = int(value)
        elif accession == "IMS:1000042" and value and state_uuid == "":
            # Legacy TTIO pre-0.9 synthetic fixtures used IMS:1000042
            # for UUID (outside scanSettings). Accept as fallback when
            # the canonical IMS:1000080 hasn't been seen yet and the
            # value looks like a UUID hex string (dashes + braces OK).
            candidate = _normalise_uuid(value)
            if len(candidate) == 32:
                state_uuid = candidate
        elif accession == "IMS:1000046" and value:  # pixel size x
            pixel_size[0] = float(value)
        elif accession == "IMS:1000047" and value:  # pixel size y
            pixel_size[1] = float(value)
        elif accession in {"IMS:1000040", "IMS:1000048"} and value:  # scan pattern / type
            scan_pattern = scan_pattern or value
        elif in_position and current is not None:
            if accession == "IMS:1000050" and value:  # position x
                current.x = int(value)
            elif accession == "IMS:1000051" and value:  # position y
                current.y = int(value)
            elif accession == "IMS:1000052" and value:  # position z
                current.z = int(value)
        elif in_binary_array and current is not None:
            if accession == "MS:1000514":  # m/z array
                array_kind = "mz"
            elif accession == "MS:1000515":  # intensity array
                array_kind = "intensity"
            elif accession == "MS:1000523":  # 64-bit float
                if array_kind == "mz":
                    current.mz_precision = "64"
                elif array_kind == "intensity":
                    current.int_precision = "64"
            elif accession == "MS:1000521":  # 32-bit float
                if array_kind == "mz":
                    current.mz_precision = "32"
                elif array_kind == "intensity":
                    current.int_precision = "32"
            elif accession == "IMS:1000102" and value:  # external offset
                if array_kind == "mz":
                    current.mz_offset = int(value)
                elif array_kind == "intensity":
                    current.int_offset = int(value)
            elif accession == "IMS:1000103" and value:  # external array length
                if array_kind == "mz":
                    current.mz_length = int(value)
                elif array_kind == "intensity":
                    current.int_length = int(value)
            elif accession == "IMS:1000104" and value:  # external encoded length
                if array_kind == "mz":
                    current.mz_encoded_length = int(value)
                elif array_kind == "intensity":
                    current.int_encoded_length = int(value)

    for event, elem in iterparse(str(imzml), events=("start", "end")):
        tag = _local(elem.tag)
        if event == "start":
            if tag == "spectrum":
                current = _SpectrumStub()
                in_spectrum = True
            elif tag == "binaryDataArray":
                in_binary_array = True
                array_kind = ""
            elif tag == "scan":
                in_position = True
            elif tag == "scanSettings":
                in_scan_settings = True
            elif tag == "referenceableParamGroup":
                in_group = elem.attrib.get("id", "")
                param_groups[in_group] = []
            continue

        # event == "end"
        if tag == "spectrum":
            if current is not None:
                spectra_stubs.append(current)
            current = None
            in_spectrum = False
        elif tag == "scanSettings":
            in_scan_settings = False
        elif tag == "binaryDataArray":
            in_binary_array = False
            array_kind = ""
        elif tag == "scan":
            in_position = False
        elif tag == "referenceableParamGroup":
            in_group = None
        elif tag == "referenceableParamGroupRef":
            for accession, value in param_groups.get(elem.attrib.get("ref", ""), ()):
                handle_cv(accession, value)
        elif tag == "cvParam":
            accession = elem.attrib.get("accession", "")
            value = elem.attrib.get("value", "")
            if in_group is not None:
                param_groups[in_group].append((accession, value))
            handle_cv(accession, value)

        elem.clear()

    if not state_mode:
        raise ImzMLParseError(f"{imzml}: no continuous/processed mode CV term found")
    if not state_uuid:
        raise ImzMLParseError(f"{imzml}: missing IMS:1000042 universally unique identifier")
    if not spectra_stubs:
        raise ImzMLParseError(f"{imzml}: no <spectrum> elements parsed")

    ibd_uuid = _read_ibd_uuid(ibd)
    if ibd_uuid != state_uuid:
        raise ImzMLBinaryError(
            f"UUID mismatch: imzML declares {state_uuid} but .ibd header is {ibd_uuid}"
        )

    ibd_size = ibd.stat().st_size
    pixels = _materialise_spectra(
        spectra_stubs, ibd, ibd_size, mode=state_mode, progress=progress,
    )

    return ImzMLImport(
        mode=state_mode,
        uuid_hex=state_uuid,
        grid_max_x=grid_max[0],
        grid_max_y=grid_max[1],
        grid_max_z=grid_max[2],
        pixel_size_x=pixel_size[0],
        pixel_size_y=pixel_size[1],
        scan_pattern=scan_pattern,
        spectra=pixels,
        source_imzml=str(imzml),
        source_ibd=str(ibd),
    )


def _normalise_uuid(value: str) -> str:
    """Strip ``{}``, dashes, whitespace and lowercase the hex string."""
    return value.replace("{", "").replace("}", "").replace("-", "").strip().lower()


def _read_ibd_uuid(ibd: Path) -> str:
    """Read the 16-byte UUID at the head of the .ibd file."""
    with ibd.open("rb") as fh:
        head = fh.read(16)
    if len(head) < 16:
        raise ImzMLBinaryError(f"{ibd}: shorter than the 16-byte UUID header")
    return head.hex()


def _materialise_spectra(
    stubs: list[_SpectrumStub],
    ibd: Path,
    ibd_size: int,
    *,
    mode: str,
    progress: ProgressSinkLike | None = None,
) -> list[ImzMLPixelSpectrum]:
    """Read every pixel's mz + intensity arrays from the .ibd.

    Fires ``progress`` every :data:`PROGRESS_INTERVAL_SPECTRA` pixels
    with ``total = -1`` and once more at the end with the true total.
    """
    pixels: list[ImzMLPixelSpectrum] = []
    shared_mz: np.ndarray | None = None

    with ibd.open("rb") as fh:
        for idx, stub in enumerate(stubs, 1):
            mz_array = _read_external_array(
                fh, stub.mz_offset, stub.mz_length, stub.mz_precision,
                ibd_size, ibd, label="m/z",
            )
            int_array = _read_external_array(
                fh, stub.int_offset, stub.int_length, stub.int_precision,
                ibd_size, ibd, label="intensity",
            )
            if mz_array.size != int_array.size:
                raise ImzMLBinaryError(
                    f"{ibd}: pixel ({stub.x},{stub.y}) mz/intensity size mismatch "
                    f"({mz_array.size} vs {int_array.size})"
                )
            if mode == "continuous":
                if shared_mz is None:
                    shared_mz = mz_array
                effective_mz = shared_mz
            else:
                effective_mz = mz_array
            pixels.append(ImzMLPixelSpectrum(
                x=stub.x, y=stub.y, z=stub.z,
                mz=effective_mz, intensity=int_array,
            ))
            if idx % PROGRESS_INTERVAL_SPECTRA == 0:
                _fire(progress, idx, -1)
    # Final fire — total is now known.
    _fire(progress, len(pixels), len(pixels))
    return pixels


def _read_external_array(
    fh,
    offset: int,
    length: int,
    precision: str,
    ibd_size: int,
    ibd: Path,
    *,
    label: str,
) -> np.ndarray:
    if offset < 0 or length < 0:
        raise ImzMLBinaryError(f"{ibd}: negative offset/length for {label} array")
    if length == 0:
        return np.empty(0, dtype=np.float64)
    bytes_per = 8 if precision == "64" else 4
    nbytes = length * bytes_per
    if offset + nbytes > ibd_size:
        raise ImzMLBinaryError(
            f"{ibd}: {label} array reads past end of file "
            f"(offset={offset}, bytes={nbytes}, size={ibd_size})"
        )
    fh.seek(offset)
    raw = fh.read(nbytes)
    if len(raw) != nbytes:
        raise ImzMLBinaryError(
            f"{ibd}: short read on {label} array at offset {offset}"
        )
    dtype = "<f8" if precision == "64" else "<f4"
    arr = np.frombuffer(raw, dtype=dtype).astype(np.float64, copy=True)
    return arr


__all__ = [
    "ImzMLImport",
    "ImzMLPixelSpectrum",
    "ImzMLParseError",
    "ImzMLBinaryError",
    "read",
]
