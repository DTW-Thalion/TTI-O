"""imzML + .ibd exporter (v0.9+).

Reverses the M59 importer: takes an :class:`ttio.importers.imzml
.ImzMLImport` (or an equivalent list of pixel spectra plus grid metadata,
or an imported pixel run read back from a ``.tio``, :func:`write_from_run`)
and emits a paired ``.imzML`` / ``.ibd`` on disk. Both the continuous
mode (one shared m/z axis at the head of the .ibd, per-pixel intensity
arrays following) and the processed mode (per-pixel m/z + intensity)
are supported.

The emitted XML uses the canonical imzML IMS accessions that match what
real-world files (e.g. the pyimzML test corpus, imzML 1.1 spec) use:
``IMS:1000080`` for the universally unique identifier and
``IMS:1000042`` / ``IMS:1000043`` for the max pixel counts. The
importer accepts both the canonical form and the legacy form that
earlier TTIO synthetic tests emitted.

SPDX-License-Identifier: Apache-2.0

Cross-language equivalents
--------------------------
Objective-C: ``TTIOImzMLWriter`` · Java:
``global.thalion.ttio.exporters.ImzMLWriter``

API status: Provisional (v0.9+ — will stabilise after the first round
of production feedback on the imaging writer).
"""
from __future__ import annotations

import hashlib
import struct
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from ..importers.imzml import (
    ImzMLImport,
    ImzMLPixelSpectrum,
    run_pixel_coordinates,
)
from ..io.progress import ProgressSinkLike, _fire


__all__ = ["write", "write_from_import", "write_from_run", "WriteResult",
           "PROGRESS_INTERVAL_SPECTRA"]


#: Mirror Java's ``ImzMLWriter.PROGRESS_INTERVAL_SPECTRA``.
PROGRESS_INTERVAL_SPECTRA = 100


@dataclass(slots=True)
class WriteResult:
    """Return value of :func:`write` — resolved paths + the UUID that
    ties the ``.imzML`` and ``.ibd`` together."""
    imzml_path: Path
    ibd_path: Path
    uuid_hex: str
    mode: str
    n_pixels: int


def write(
    pixels: Iterable[ImzMLPixelSpectrum],
    imzml_path: str | Path,
    ibd_path: str | Path | None = None,
    *,
    mode: str = "continuous",
    grid_max_x: int = 0,
    grid_max_y: int = 0,
    grid_max_z: int = 1,
    pixel_size_x: float = 0.0,
    pixel_size_y: float = 0.0,
    scan_pattern: str = "flyback",
    uuid_hex: str | None = None,
    progress: ProgressSinkLike | None = None,
) -> WriteResult:
    """Write ``pixels`` as an imzML + .ibd pair.

    Parameters
    ----------
    pixels
        Iterable of :class:`ImzMLPixelSpectrum`. In ``mode="continuous"``
        the first pixel's ``mz`` array is written once at the head of
        the .ibd and every subsequent pixel is checked to share that
        exact axis (raises if not). In ``mode="processed"`` every
        pixel's ``mz`` and ``intensity`` arrays are stored in full.
    imzml_path
        Path to the ``.imzML`` output. The sibling ``.ibd`` is derived
        by swapping the extension unless ``ibd_path`` is given.
    mode
        ``"continuous"`` or ``"processed"``. Default is continuous.
    grid_max_x, grid_max_y, grid_max_z
        Pixel grid extents. When ``0``, the writer derives them from
        the actual pixel coordinates.
    pixel_size_x, pixel_size_y
        Spatial pixel sizes in micrometres; emitted as-is via
        ``IMS:1000046`` / ``IMS:1000047`` cvParams.
    scan_pattern
        Free-text scan pattern label (e.g. ``"flyback"``,
        ``"meandering"``). Default ``"flyback"``.
    uuid_hex
        Optional explicit UUID (32 lowercase hex characters or a
        standard UUID string). When ``None``, a random UUID4 is
        generated.
    """
    if mode not in ("continuous", "processed"):
        raise ValueError(
            f"imzML writer mode must be 'continuous' or 'processed', got {mode!r}"
        )

    imzml = Path(imzml_path)
    ibd = Path(ibd_path) if ibd_path else imzml.with_suffix(".ibd")
    pixel_list = list(pixels)
    if not pixel_list:
        raise ValueError("imzML writer: at least one pixel spectrum is required")

    if uuid_hex is None:
        uuid_hex = uuid.uuid4().hex
    else:
        uuid_hex = _normalise_uuid(uuid_hex)
        if len(uuid_hex) != 32:
            raise ValueError(f"uuid_hex must be 32 hex chars, got {uuid_hex!r}")

    # Derive grid extents from pixel coordinates when not specified.
    if grid_max_x == 0:
        grid_max_x = max(p.x for p in pixel_list)
    if grid_max_y == 0:
        grid_max_y = max(p.y for p in pixel_list)
    if grid_max_z == 0:
        grid_max_z = max((p.z for p in pixel_list), default=1)

    # ── .ibd layout ────────────────────────────────────────────────────
    # Byte 0-15: UUID (raw 16 bytes, not hex)
    # Continuous mode: shared m/z axis, then intensity[0], intensity[1], ...
    # Processed mode: mz[0], intensity[0], mz[1], intensity[1], ...
    # m/z and intensity are both little-endian float64.
    ibd_parts: list[bytes] = [bytes.fromhex(uuid_hex)]
    offsets: list[tuple[int, int, int, int]] = []
    #   per pixel: (mz_offset, mz_length, int_offset, int_length)

    cursor = 16  # bytes

    if mode == "continuous":
        # Validate shared axis, write once.
        shared_mz = np.ascontiguousarray(pixel_list[0].mz, dtype="<f8")
        shared_mz_bytes = shared_mz.tobytes()
        mz_offset = cursor
        mz_len = int(shared_mz.size)
        ibd_parts.append(shared_mz_bytes)
        cursor += len(shared_mz_bytes)
        n_pixels = len(pixel_list)
        for i, pix in enumerate(pixel_list):
            mz_i = np.ascontiguousarray(pix.mz, dtype="<f8")
            if mz_i.size != mz_len or not np.array_equal(mz_i, shared_mz):
                raise ValueError(
                    f"continuous-mode imzML requires all pixels to share the "
                    f"same m/z axis; pixel {i} (x={pix.x}, y={pix.y}) differs"
                )
            inten = np.ascontiguousarray(pix.intensity, dtype="<f8")
            inten_bytes = inten.tobytes()
            int_offset = cursor
            ibd_parts.append(inten_bytes)
            cursor += len(inten_bytes)
            offsets.append((mz_offset, mz_len, int_offset, int(inten.size)))
            if (i + 1) % PROGRESS_INTERVAL_SPECTRA == 0:
                _fire(progress, i + 1, n_pixels)
    else:  # processed
        n_pixels = len(pixel_list)
        for i, pix in enumerate(pixel_list):
            mz_i = np.ascontiguousarray(pix.mz, dtype="<f8")
            inten = np.ascontiguousarray(pix.intensity, dtype="<f8")
            if mz_i.size != inten.size:
                raise ValueError(
                    f"processed-mode pixel at (x={pix.x}, y={pix.y}): mz ({mz_i.size})"
                    f" and intensity ({inten.size}) arrays must be the same length"
                )
            mz_bytes = mz_i.tobytes()
            mz_offset = cursor
            ibd_parts.append(mz_bytes)
            cursor += len(mz_bytes)
            inten_bytes = inten.tobytes()
            int_offset = cursor
            ibd_parts.append(inten_bytes)
            cursor += len(inten_bytes)
            offsets.append((mz_offset, int(mz_i.size),
                            int_offset, int(inten.size)))
            if (i + 1) % PROGRESS_INTERVAL_SPECTRA == 0:
                _fire(progress, i + 1, n_pixels)

    ibd_bytes = b"".join(ibd_parts)
    ibd.write_bytes(ibd_bytes)
    # IMS:1000091 "ibd SHA-1" is a common audit cvParam; compute it so
    # downstream validators see a match.
    ibd_sha1 = hashlib.sha1(ibd_bytes).hexdigest()

    # ── .imzML XML ────────────────────────────────────────────────────
    xml = _build_imzml_xml(
        uuid_hex=uuid_hex,
        ibd_sha1=ibd_sha1,
        mode=mode,
        grid_max_x=grid_max_x,
        grid_max_y=grid_max_y,
        grid_max_z=grid_max_z,
        pixel_size_x=pixel_size_x,
        pixel_size_y=pixel_size_y,
        scan_pattern=scan_pattern,
        pixels=pixel_list,
        offsets=offsets,
    )
    imzml.write_text(xml, encoding="utf-8")

    _fire(progress, len(pixel_list), len(pixel_list))
    return WriteResult(
        imzml_path=imzml,
        ibd_path=ibd,
        uuid_hex=uuid_hex,
        mode=mode,
        n_pixels=len(pixel_list),
    )


def write_from_import(
    import_result: ImzMLImport,
    imzml_path: str | Path,
    ibd_path: str | Path | None = None,
    *,
    progress: ProgressSinkLike | None = None,
) -> WriteResult:
    """Re-emit an :class:`ImzMLImport` to disk.

    Convenience shortcut for round-trip tests and imzML-to-imzML
    rewriters (e.g. re-baselining a corrupt UUID).
    """
    return write(
        import_result.spectra,
        imzml_path=imzml_path,
        ibd_path=ibd_path,
        mode=import_result.mode,
        grid_max_x=import_result.grid_max_x,
        grid_max_y=import_result.grid_max_y,
        grid_max_z=import_result.grid_max_z,
        pixel_size_x=import_result.pixel_size_x,
        pixel_size_y=import_result.pixel_size_y,
        scan_pattern=import_result.scan_pattern,
        uuid_hex=import_result.uuid_hex,
        progress=progress,
    )


def write_from_run(
    run,
    imzml_path: str | Path,
    ibd_path: str | Path | None = None,
    *,
    dataset_provenance=(),
    progress: ProgressSinkLike | None = None,
) -> WriteResult:
    """Write an imaging pixel run read from a ``.tio`` back to imzML.

    ``run`` is an :class:`~ttio.acquisition_run.AcquisitionRun` whose
    spectra are pixels: the imzML importer's ``imzml_pixels`` run, or
    any run with the M102 ``spectrum_index`` pixel columns. Positions
    come from :func:`ttio.importers.imzml.run_pixel_coordinates`, so a
    pre-M102 run that kept them in the ``imzml_pixel_coordinates_csv``
    provenance parameter also works. The storage mode, UUID, grid
    extents, pixel size and scan pattern come from the importer's
    ``imzml_*`` provenance parameters (run records first, then
    ``dataset_provenance``) when present. Without ``imzml_mode`` the
    mode is continuous when every pixel shares one m/z axis and
    processed otherwise.

    Raises ``ValueError`` when the run carries no pixel positions.
    """
    coords = run_pixel_coordinates(run, dataset_provenance)
    if coords is None:
        raise ValueError(
            f"run {run.name!r} has no pixel coordinates to export as imzML")

    # Run records win over dataset records, later records over earlier.
    params: dict = {}
    for record in reversed([*dataset_provenance, *run.provenance()]):
        for key, value in record.parameters.items():
            if key.startswith("imzml_"):
                params.setdefault(key, value)

    pixels: list[ImzMLPixelSpectrum] = []
    for i, (x, y, z) in enumerate(coords):
        arrays = run[i].signal_arrays
        pixels.append(ImzMLPixelSpectrum(
            x=int(x), y=int(y), z=int(z),
            mz=np.asarray(arrays["mz"].data, dtype=np.float64),
            intensity=np.asarray(arrays["intensity"].data, dtype=np.float64),
        ))

    mode = params.get("imzml_mode")
    if mode not in ("continuous", "processed"):
        first = pixels[0].mz if pixels else np.zeros(0)
        shared = all(np.array_equal(p.mz, first) for p in pixels)
        mode = "continuous" if shared else "processed"

    return write(
        pixels,
        imzml_path=imzml_path,
        ibd_path=ibd_path,
        mode=mode,
        grid_max_x=int(params.get("imzml_grid_max_x", 0) or 0),
        grid_max_y=int(params.get("imzml_grid_max_y", 0) or 0),
        grid_max_z=int(params.get("imzml_grid_max_z", 0) or 0),
        pixel_size_x=float(params.get("imzml_pixel_size_x", 0.0) or 0.0),
        pixel_size_y=float(params.get("imzml_pixel_size_y", 0.0) or 0.0),
        scan_pattern=str(params.get("imzml_scan_pattern") or "flyback"),
        uuid_hex=params.get("imzml_uuid_hex") or None,
        progress=progress,
    )


# --------------------------------------------------------------------------- #
# Helpers.
# --------------------------------------------------------------------------- #

def _normalise_uuid(value: str) -> str:
    """Strip ``{}``, dashes, whitespace and lowercase the hex string."""
    return value.replace("{", "").replace("}", "").replace("-", "").strip().lower()


def _build_imzml_xml(
    *,
    uuid_hex: str,
    ibd_sha1: str,
    mode: str,
    grid_max_x: int,
    grid_max_y: int,
    grid_max_z: int,
    pixel_size_x: float,
    pixel_size_y: float,
    scan_pattern: str,
    pixels: list[ImzMLPixelSpectrum],
    offsets: list[tuple[int, int, int, int]],
) -> str:
    mode_acc = "IMS:1000030" if mode == "continuous" else "IMS:1000031"
    mode_name = "continuous" if mode == "continuous" else "processed"

    parts: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>\n',
        '<mzML xmlns="http://psi.hupo.org/ms/mzml"'
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
        ' xsi:schemaLocation="http://psi.hupo.org/ms/mzml'
        ' http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"'
        ' version="1.1">\n',
        '  <cvList count="3">\n',
        '    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"'
        ' version="4.1.0"'
        ' URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>\n',
        '    <cv id="UO" fullName="Unit Ontology" version="2020-03-10"'
        ' URI="http://ontologies.berkeleybop.org/uo.obo"/>\n',
        '    <cv id="IMS" fullName="Mass Spectrometry Imaging Ontology"'
        ' version="1.1.0"'
        ' URI="https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo"/>\n',
        '  </cvList>\n',
        '  <fileDescription>\n',
        '    <fileContent>\n',
        '      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>\n',
        f'      <cvParam cvRef="IMS" accession="IMS:1000080"'
        f' name="universally unique identifier" value="{uuid_hex}"/>\n',
        f'      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1"'
        f' value="{ibd_sha1}"/>\n',
        f'      <cvParam cvRef="IMS" accession="{mode_acc}" name="{mode_name}" value=""/>\n',
        '    </fileContent>\n',
        '  </fileDescription>\n',
        # Real-world imzML tools (pyimzml, MSIqr, etc.) resolve every
        # <binaryDataArray>'s meaning via <referenceableParamGroupRef>
        # and fail hard when the ref target is missing. Emit the two
        # group declarations here; each per-spectrum binaryDataArray
        # below also references one of these IDs AND emits inline
        # cvParams so simple readers (including our own importer) work
        # without a ref-resolution pass.
        '  <referenceableParamGroupList count="2">\n',
        '    <referenceableParamGroup id="mzArray">\n',
        '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>\n',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>\n',
        '      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"'
        ' unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n',
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>\n',
        '    </referenceableParamGroup>\n',
        '    <referenceableParamGroup id="intensityArray">\n',
        '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>\n',
        '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>\n',
        '      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"'
        ' unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>\n',
        '      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>\n',
        '    </referenceableParamGroup>\n',
        '  </referenceableParamGroupList>\n',
        '  <softwareList count="1">\n',
        '    <software id="ttio" version="0.9.0">\n',
        '      <cvParam cvRef="MS" accession="MS:1000799"'
        ' name="custom unreleased software tool" value="ttio"/>\n',
        '    </software>\n',
        '  </softwareList>\n',
        '  <scanSettingsList count="1">\n',
        '    <scanSettings id="scansettings1">\n',
        # Scan pattern — IMS:1000040 is what our importer looks for,
        # even though pyimzml flags its canonical name as "linescan
        # sequence". Emit both forms (cvParam + userParam) so every
        # reader finds the value without triggering ontology warnings.
        f'      <userParam name="scan pattern"'
        f' value="{_xml_escape(scan_pattern)}" type="xsd:string"/>\n',
        f'      <cvParam cvRef="IMS" accession="IMS:1000040"'
        f' name="linescan sequence" value="{_xml_escape(scan_pattern)}"/>\n',
        f'      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixels x"'
        f' value="{grid_max_x}"/>\n',
        f'      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixels y"'
        f' value="{grid_max_y}"/>\n',
    ]
    if pixel_size_x > 0.0:
        parts.append(
            f'      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size (x)"'
            f' value="{pixel_size_x:g}"'
            f' unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>\n'
        )
    if pixel_size_y > 0.0:
        parts.append(
            f'      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y"'
            f' value="{pixel_size_y:g}"'
            f' unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>\n'
        )
    parts += [
        '    </scanSettings>\n',
        '  </scanSettingsList>\n',
        '  <instrumentConfigurationList count="1">\n',
        '    <instrumentConfiguration id="IC1">\n',
        '      <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value=""/>\n',
        '    </instrumentConfiguration>\n',
        '  </instrumentConfigurationList>\n',
        '  <dataProcessingList count="1">\n',
        '    <dataProcessing id="dp_export">\n',
        '      <processingMethod order="0" softwareRef="ttio">\n',
        '        <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML"/>\n',
        '      </processingMethod>\n',
        '    </dataProcessing>\n',
        '  </dataProcessingList>\n',
        '  <run id="ttio_imzml_export" defaultInstrumentConfigurationRef="IC1">\n',
        f'    <spectrumList count="{len(pixels)}" defaultDataProcessingRef="dp_export">\n',
    ]

    # Per-spectrum (pixel) blocks.
    for i, pix in enumerate(pixels):
        mz_offset, mz_length, int_offset, int_length = offsets[i]
        mz_enc_len = mz_length * 8
        int_enc_len = int_length * 8
        parts.append(
            f'      <spectrum id="Scan={i + 1}" index="{i}" defaultArrayLength="0">\n'
        )
        parts.append(
            '        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>\n'
        )
        parts.append(
            '        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>\n'
        )
        parts.append('        <scanList count="1">\n')
        parts.append(
            '          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>\n'
        )
        parts.append('          <scan instrumentConfigurationRef="IC1">\n')
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000050"'
            f' name="position x" value="{int(pix.x)}"/>\n'
        )
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000051"'
            f' name="position y" value="{int(pix.y)}"/>\n'
        )
        if int(pix.z) != 1:
            parts.append(
                f'            <cvParam cvRef="IMS" accession="IMS:1000052"'
                f' name="position z" value="{int(pix.z)}"/>\n'
            )
        parts.append('          </scan>\n')
        parts.append('        </scanList>\n')
        parts.append('        <binaryDataArrayList count="2">\n')
        # m/z array — the referenceableParamGroupRef lets strict tools
        # (pyimzml, MSIqr) identify which binaryDataArray is which; the
        # inline cvParams let simple readers identify it without
        # resolving the ref.
        parts.append('          <binaryDataArray encodedLength="0">\n')
        parts.append('            <referenceableParamGroupRef ref="mzArray"/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000523"'
                     ' name="64-bit float" value=""/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000576"'
                     ' name="no compression" value=""/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000514"'
                     ' name="m/z array" value=""'
                     ' unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n')
        parts.append('            <cvParam cvRef="IMS" accession="IMS:1000101"'
                     ' name="external data" value="true"/>\n')
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000102"'
            f' name="external offset" value="{mz_offset}"/>\n'
        )
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000103"'
            f' name="external array length" value="{mz_length}"/>\n'
        )
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000104"'
            f' name="external encoded length" value="{mz_enc_len}"/>\n'
        )
        parts.append('            <binary/>\n')
        parts.append('          </binaryDataArray>\n')
        # intensity array
        parts.append('          <binaryDataArray encodedLength="0">\n')
        parts.append('            <referenceableParamGroupRef ref="intensityArray"/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000523"'
                     ' name="64-bit float" value=""/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000576"'
                     ' name="no compression" value=""/>\n')
        parts.append('            <cvParam cvRef="MS" accession="MS:1000515"'
                     ' name="intensity array" value=""'
                     ' unitCvRef="MS" unitAccession="MS:1000131"'
                     ' unitName="number of detector counts"/>\n')
        parts.append('            <cvParam cvRef="IMS" accession="IMS:1000101"'
                     ' name="external data" value="true"/>\n')
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000102"'
            f' name="external offset" value="{int_offset}"/>\n'
        )
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000103"'
            f' name="external array length" value="{int_length}"/>\n'
        )
        parts.append(
            f'            <cvParam cvRef="IMS" accession="IMS:1000104"'
            f' name="external encoded length" value="{int_enc_len}"/>\n'
        )
        parts.append('            <binary/>\n')
        parts.append('          </binaryDataArray>\n')
        parts.append('        </binaryDataArrayList>\n')
        parts.append('      </spectrum>\n')

    parts += [
        '    </spectrumList>\n',
        '  </run>\n',
        '</mzML>\n',
    ]
    return "".join(parts)


def _xml_escape(s: str) -> str:
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;").replace("'", "&apos;"))
