"""imzML + .ibd importer integration tests (v0.9 M59).

Covers HANDOFF M59 acceptance:

* Continuous-mode reference imports correctly (single shared m/z axis,
  varying intensity per pixel).
* Processed-mode reference imports correctly (per-pixel m/z arrays).
* Spatial coordinates map onto the pixel grid.
* m/z and intensity bytes round-trip through ``ImzMLImport`` to a
  readable .tio via :meth:`SpectralDataset.open`.
* Hard error on UUID mismatch (gotcha 49).
* Hard error on offsets that read past .ibd EOF (gotcha 48).

Synthetic-fixture only — pyimzML is not required and the tests stay
hermetic (no network). The pinned imzml.org reference files are
TBD entries in ``tests/fixtures/download.py``; their wiring lands
when the URLs are confirmed.
"""
from __future__ import annotations

import random
import struct
import uuid
from pathlib import Path

import numpy as np
import pytest

from ttio import SpectralDataset
from ttio.importers.imzml import (
    ImzMLBinaryError,
    ImzMLImport,
    ImzMLParseError,
    read as imzml_read,
)


# --------------------------------------------------------------------------- #
# Synthetic .imzML + .ibd builder.
# --------------------------------------------------------------------------- #

def _build_pair(
    tmp_path: Path,
    *,
    mode: str,
    grid_x: int,
    grid_y: int,
    n_peaks: int,
    pixel_size: tuple[float, float] = (10.0, 10.0),
    rng: np.random.Generator | None = None,
    uuid_override_ibd: bytes | None = None,
    truncate_ibd_to: int | None = None,
) -> tuple[Path, Path, dict]:
    """Build a deterministic imzML + .ibd pair.

    Returns the imzML path, .ibd path, and a metadata dict the test
    can use for assertions (uuid_hex, expected mz arrays per pixel,
    expected intensities, grid extents).
    """
    if rng is None:
        rng = np.random.default_rng(20260418)
    # Pull 128 random bits via stdlib so we don't fight numpy's uint64
    # bound. The seed is derived from rng so the fixture stays
    # deterministic per test.
    seed = int(rng.integers(0, 2**32 - 1, dtype=np.uint32))
    file_uuid = uuid.UUID(int=random.Random(seed).getrandbits(128))
    uuid_hex = file_uuid.hex
    uuid_bytes = file_uuid.bytes

    n_pixels = grid_x * grid_y

    # Build .ibd payload.
    payload = bytearray()
    payload += uuid_override_ibd if uuid_override_ibd is not None else uuid_bytes

    expected_mz: list[np.ndarray] = []
    expected_intensity: list[np.ndarray] = []
    mz_offsets: list[int] = []
    int_offsets: list[int] = []

    # Continuous: one shared m/z block referenced by every pixel.
    shared_mz = np.linspace(100.0, 1000.0, n_peaks)
    if mode == "continuous":
        shared_offset = len(payload)
        payload += shared_mz.astype("<f8").tobytes()

    for pixel in range(n_pixels):
        if mode == "continuous":
            mz_offset = shared_offset
            mz = shared_mz
        else:
            mz_offset = len(payload)
            mz = np.linspace(100.0 + pixel, 1000.0 + pixel, n_peaks)
            payload += mz.astype("<f8").tobytes()
        mz_offsets.append(mz_offset)
        expected_mz.append(mz)

        int_offset = len(payload)
        intensity = rng.uniform(0.0, 1e6, size=n_peaks)
        payload += intensity.astype("<f8").tobytes()
        int_offsets.append(int_offset)
        expected_intensity.append(intensity)

    # Build .imzML XML.
    mode_acc = "IMS:1000030" if mode == "continuous" else "IMS:1000031"
    spectrum_blocks: list[str] = []
    for pixel in range(n_pixels):
        x = (pixel % grid_x) + 1
        y = (pixel // grid_x) + 1
        spectrum_blocks.append(f"""    <spectrum index="{pixel}" id="px={pixel}">
      <scanList count="1">
        <scan>
          <cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="{x}"/>
          <cvParam cvRef="IMS" accession="IMS:1000051" name="position y" value="{y}"/>
        </scan>
      </scanList>
      <binaryDataArrayList count="2">
        <binaryDataArray encodedLength="{n_peaks * 8}">
          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
          <cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>
          <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="{mz_offsets[pixel]}"/>
          <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="{n_peaks}"/>
          <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="{n_peaks * 8}"/>
        </binaryDataArray>
        <binaryDataArray encodedLength="{n_peaks * 8}">
          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
          <cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>
          <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="{int_offsets[pixel]}"/>
          <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="{n_peaks}"/>
          <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="{n_peaks * 8}"/>
        </binaryDataArray>
      </binaryDataArrayList>
    </spectrum>""")

    imzml_text = f"""<?xml version="1.0" encoding="UTF-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <fileDescription>
    <fileContent>
      <cvParam cvRef="IMS" accession="IMS:1000042" name="universally unique identifier" value="{{{file_uuid}}}"/>
      <cvParam cvRef="IMS" accession="{mode_acc}" name="{mode} mode"/>
    </fileContent>
  </fileDescription>
  <scanSettingsList count="1">
    <scanSettings id="scanset1">
      <cvParam cvRef="IMS" accession="IMS:1000003" name="max count of pixels x" value="{grid_x}"/>
      <cvParam cvRef="IMS" accession="IMS:1000004" name="max count of pixels y" value="{grid_y}"/>
      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size x" value="{pixel_size[0]}"/>
      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" value="{pixel_size[1]}"/>
      <cvParam cvRef="IMS" accession="IMS:1000040" name="scan pattern" value="flyback"/>
    </scanSettings>
  </scanSettingsList>
  <run id="ims_run">
    <spectrumList count="{n_pixels}">
{chr(10).join(spectrum_blocks)}
    </spectrumList>
  </run>
</mzML>
"""
    imzml_path = tmp_path / f"synth_{mode}.imzML"
    ibd_path = tmp_path / f"synth_{mode}.ibd"
    imzml_path.write_text(imzml_text)
    if truncate_ibd_to is not None:
        ibd_path.write_bytes(bytes(payload[:truncate_ibd_to]))
    else:
        ibd_path.write_bytes(bytes(payload))

    return imzml_path, ibd_path, {
        "uuid_hex": uuid_hex,
        "expected_mz": expected_mz,
        "expected_intensity": expected_intensity,
        "grid_x": grid_x,
        "grid_y": grid_y,
        "n_peaks": n_peaks,
    }


# --------------------------------------------------------------------------- #
# Continuous-mode happy path.
# --------------------------------------------------------------------------- #

def test_continuous_mode_pixel_count_and_grid(tmp_path: Path) -> None:
    imzml, _ibd, meta = _build_pair(tmp_path, mode="continuous", grid_x=3, grid_y=2, n_peaks=8)
    result = imzml_read(imzml)
    assert result.mode == "continuous"
    assert result.uuid_hex == meta["uuid_hex"]
    assert len(result.spectra) == meta["grid_x"] * meta["grid_y"]
    assert result.grid_max_x == meta["grid_x"]
    assert result.grid_max_y == meta["grid_y"]
    assert result.scan_pattern == "flyback"
    # First pixel sits at logical (1,1); last at (grid_x, grid_y).
    assert (result.spectra[0].x, result.spectra[0].y) == (1, 1)
    assert (result.spectra[-1].x, result.spectra[-1].y) == (meta["grid_x"], meta["grid_y"])


def test_continuous_mode_shares_mz_axis(tmp_path: Path) -> None:
    """All pixels in continuous mode reference the same m/z block."""
    imzml, _ibd, meta = _build_pair(tmp_path, mode="continuous", grid_x=2, grid_y=2, n_peaks=16)
    result = imzml_read(imzml)
    expected_mz = meta["expected_mz"][0]
    for spec in result.spectra:
        np.testing.assert_array_equal(spec.mz, expected_mz)
    # Every pixel should literally share the array (continuous mode
    # contract — proves we aren't materialising N copies).
    first_mz = result.spectra[0].mz
    assert all(spec.mz is first_mz for spec in result.spectra)


def test_continuous_mode_intensities_per_pixel(tmp_path: Path) -> None:
    imzml, _ibd, meta = _build_pair(tmp_path, mode="continuous", grid_x=2, grid_y=2, n_peaks=8)
    result = imzml_read(imzml)
    for i, spec in enumerate(result.spectra):
        np.testing.assert_array_equal(spec.intensity, meta["expected_intensity"][i])


# --------------------------------------------------------------------------- #
# Processed-mode happy path.
# --------------------------------------------------------------------------- #

def test_processed_mode_per_pixel_mz(tmp_path: Path) -> None:
    imzml, _ibd, meta = _build_pair(tmp_path, mode="processed", grid_x=2, grid_y=3, n_peaks=12)
    result = imzml_read(imzml)
    assert result.mode == "processed"
    assert len(result.spectra) == 6
    for i, spec in enumerate(result.spectra):
        np.testing.assert_array_equal(spec.mz, meta["expected_mz"][i])
        np.testing.assert_array_equal(spec.intensity, meta["expected_intensity"][i])
    # Sanity: each pixel got its own m/z array, not the first.
    assert result.spectra[1].mz is not result.spectra[0].mz


# --------------------------------------------------------------------------- #
# Error contracts.
# --------------------------------------------------------------------------- #

def test_uuid_mismatch_raises(tmp_path: Path) -> None:
    """When the .ibd UUID header disagrees with the .imzML CV term the
    importer fails fast (HANDOFF gotcha 49)."""
    bad_uuid_bytes = uuid.UUID("ffffffff-ffff-ffff-ffff-ffffffffffff").bytes
    imzml, _ibd, _meta = _build_pair(
        tmp_path, mode="continuous", grid_x=1, grid_y=1, n_peaks=4,
        uuid_override_ibd=bad_uuid_bytes,
    )
    with pytest.raises(ImzMLBinaryError, match="UUID mismatch"):
        imzml_read(imzml)


def test_offset_overflow_raises(tmp_path: Path) -> None:
    """When external_offset + bytes exceed the .ibd size the importer
    refuses the read (HANDOFF gotcha 48)."""
    imzml, _ibd, _meta = _build_pair(
        tmp_path, mode="continuous", grid_x=1, grid_y=1, n_peaks=8,
        truncate_ibd_to=20,  # keeps the UUID, kills the binary tail
    )
    with pytest.raises(ImzMLBinaryError, match="reads past end"):
        imzml_read(imzml)


def test_missing_imzml_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        imzml_read(tmp_path / "absent.imzML")


def test_missing_ibd_raises(tmp_path: Path) -> None:
    bare = tmp_path / "bare.imzML"
    bare.write_text("<?xml version='1.0'?><mzML/>")
    with pytest.raises(FileNotFoundError):
        imzml_read(bare)


def test_no_mode_cv_raises(tmp_path: Path) -> None:
    bare = tmp_path / "no_mode.imzML"
    ibd = tmp_path / "no_mode.ibd"
    bare.write_text(
        "<?xml version='1.0'?><mzML><fileDescription><fileContent>"
        "<cvParam accession='IMS:1000042' value='{00000000-0000-0000-0000-000000000000}'/>"
        "</fileContent></fileDescription></mzML>"
    )
    ibd.write_bytes(b"\x00" * 16)
    with pytest.raises(ImzMLParseError, match="no continuous/processed mode"):
        imzml_read(bare)


# --------------------------------------------------------------------------- #
# Round-trip into .tio and back through SpectralDataset.
# --------------------------------------------------------------------------- #

def test_to_ttio_round_trip(tmp_path: Path) -> None:
    imzml, _ibd, meta = _build_pair(
        tmp_path, mode="processed", grid_x=2, grid_y=2, n_peaks=10,
    )
    parsed = imzml_read(imzml)
    out = tmp_path / "imaging.tio"
    parsed.to_ttio(out, title="synth IMS")
    with SpectralDataset.open(out) as ds:
        run = ds.ms_runs["imzml_pixels"]
        assert len(run) == 4
        assert ds.title == "synth IMS"

        # First pixel data round-trips bit-for-bit.
        spec = run[0]
        np.testing.assert_array_equal(
            spec.signal_arrays["mz"].data, meta["expected_mz"][0]
        )
        np.testing.assert_array_equal(
            spec.signal_arrays["intensity"].data, meta["expected_intensity"][0]
        )

        # Provenance preserves the imzML grid metadata.
        prov = ds.provenance()
        assert len(prov) >= 1
        params = prov[0].parameters
        assert params["imzml_mode"] == "processed"
        assert params["imzml_grid_max_x"] == 2
        assert params["imzml_grid_max_y"] == 2
        assert params["imzml_uuid_hex"] == meta["uuid_hex"]
        # M102: positions are spectrum_index columns in scan order,
        # not a provenance parameter.
        assert "imzml_pixel_coordinates_csv" not in params
        assert ds.feature_flags.has("opt_pixel_coordinates")
        assert run.index.pixel_coordinates_at(0) == (1, 1, 1)
        assert run.index.pixel_coordinates_at(3) == (2, 2, 1)


def test_continuous_round_trip_records_mode(tmp_path: Path) -> None:
    imzml, _ibd, _meta = _build_pair(
        tmp_path, mode="continuous", grid_x=2, grid_y=1, n_peaks=4,
    )
    parsed = imzml_read(imzml)
    out = tmp_path / "cont.tio"
    parsed.to_ttio(out)
    with SpectralDataset.open(out) as ds:
        prov = ds.provenance()
        assert prov[0].parameters["imzml_mode"] == "continuous"
