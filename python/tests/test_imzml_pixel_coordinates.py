"""M102: imzML pixel positions as spectrum_index columns.

Before M102 the importer put every pixel's position in one provenance
parameter, ``imzml_pixel_coordinates_csv``. Both provenance levels are
mirrored into a fixed-length ``@provenance_json`` attribute, which HDF5
caps at 64 KB, so an import failed with "object header message is too
large" once the CSV passed 64 KB: from about 9,000 pixels on a 100-pixel
wide grid (PRIDE PXD001283 has 34,840 pixels, a 306 KB CSV).
The positions now live in ``spectrum_index/pixel_x|pixel_y|pixel_z``
(format-spec §4b, ``opt_pixel_coordinates``); older files that only
carry the CSV still read and export.
"""
from __future__ import annotations

import time
from pathlib import Path

import h5py
import numpy as np
import pytest

from ttio import SpectralDataset, WrittenRun
from ttio.exporters import imzml as imzml_writer
from ttio.exporters import registry as export_registry
from ttio.importers import imzml as imzml_reader
from ttio.importers import registry as import_registry
from ttio.importers.imzml import ImzMLPixelSpectrum, run_pixel_coordinates
from ttio.provenance import ProvenanceRecord

LEGACY = imzml_reader.LEGACY_COORDINATES_PARAMETER


def _pixels(width: int, height: int, seed: int = 7) -> list[ImzMLPixelSpectrum]:
    """Processed-mode pixels in flyback order, 3-5 peaks each with
    per-pixel m/z values."""
    rng = np.random.default_rng(seed)
    out = []
    for y in range(1, height + 1):
        for x in range(1, width + 1):
            n = int(rng.integers(3, 6))
            mz = np.sort(rng.uniform(100.0, 1000.0, n))
            out.append(ImzMLPixelSpectrum(
                x=x, y=y, z=1, mz=mz, intensity=rng.uniform(1.0, 1e5, n)))
    return out


def _assert_no_legacy_parameter(ds: SpectralDataset, run) -> None:
    for record in [*ds.provenance(), *run.provenance()]:
        assert LEGACY not in record.parameters


def test_large_pixel_count_imports_exports_and_round_trips(tmp_path: Path) -> None:
    # 150 x 100 = 15,000 pixels: a 123 KB CSV. Pre-M102 this import
    # raised OSError (measured: 8,000 pixels / 62 KB still wrote,
    # 9,600 / 77 KB failed).
    src = _pixels(150, 100)
    src_imzml = tmp_path / "big.imzML"
    imzml_writer.write(src, src_imzml, mode="processed",
                       pixel_size_x=20.0, pixel_size_y=20.0)

    tio = tmp_path / "big.tio"
    import_registry.encode("imzml", [str(src_imzml)], str(tio))

    expected = np.array([(p.x, p.y, p.z) for p in src], dtype=np.int32)
    with SpectralDataset.open(tio) as ds:
        assert ds.feature_flags.has("opt_pixel_coordinates")
        run = ds.ms_runs["imzml_pixels"]
        assert len(run) == 15000
        assert run.index.has_pixel_coordinates
        np.testing.assert_array_equal(run.index.pixel_x, expected[:, 0])
        np.testing.assert_array_equal(run.index.pixel_y, expected[:, 1])
        np.testing.assert_array_equal(run.index.pixel_z, expected[:, 2])
        np.testing.assert_array_equal(run_pixel_coordinates(run), expected)
        _assert_no_legacy_parameter(ds, run)
        params = run.provenance()[0].parameters
        assert params["imzml_mode"] == "processed"
        assert params["imzml_grid_max_x"] == 150
        assert params["imzml_grid_max_y"] == 100

    out = tmp_path / "back.imzML"
    export_registry.export("imzml", str(tio), None, str(out))
    back = imzml_reader.read(out)
    assert back.mode == "processed"
    assert (back.grid_max_x, back.grid_max_y) == (150, 100)
    assert back.pixel_size_x == 20.0
    assert len(back.spectra) == 15000
    for a, b in zip(src, back.spectra):
        assert (a.x, a.y, a.z) == (b.x, b.y, b.z)
        np.testing.assert_array_equal(a.mz, b.mz)
        np.testing.assert_array_equal(a.intensity, b.intensity)


def test_continuous_import_keeps_mode_and_positions(tmp_path: Path) -> None:
    rng = np.random.default_rng(3)
    axis = np.linspace(100.0, 200.0, 8)
    src = [ImzMLPixelSpectrum(x=x, y=y, z=1, mz=axis,
                              intensity=rng.uniform(0, 10, 8))
           for y in (1, 2) for x in (1, 2, 3)]
    src_imzml = tmp_path / "cont.imzML"
    written = imzml_writer.write(src, src_imzml, mode="continuous")

    tio = tmp_path / "cont.tio"
    imzml_reader.read(src_imzml).to_ttio(tio)
    out = tmp_path / "back.imzML"
    with SpectralDataset.open(tio) as ds:
        run = ds.ms_runs["imzml_pixels"]
        imzml_writer.write_from_run(run, out,
                                    dataset_provenance=ds.provenance())
    back = imzml_reader.read(out)
    assert back.mode == "continuous"
    assert back.uuid_hex == written.uuid_hex
    assert [(p.x, p.y) for p in back.spectra] == [(p.x, p.y) for p in src]


def _legacy_tio(path: Path, coords: list[tuple[int, int, int]], *,
                run_level: bool, csv: str | None = None) -> Path:
    """A pre-M102-shaped file: no pixel columns, positions in the CSV."""
    n = len(coords)
    mz = np.tile(np.array([100.0, 150.0, 200.0]), n)
    intensity = np.arange(3 * n, dtype=np.float64) + 1.0
    params = {
        "imzml_mode": "continuous",
        "imzml_uuid_hex": "0123456789abcdef0123456789abcdef",
        "imzml_grid_max_x": 2,
        "imzml_grid_max_y": 2,
        "imzml_grid_max_z": 1,
        "imzml_pixel_size_x": 5.0,
        "imzml_pixel_size_y": 5.0,
        "imzml_scan_pattern": "flyback",
        LEGACY: csv if csv is not None else ";".join(
            f"{x},{y},{z}" for x, y, z in coords),
    }
    prov = ProvenanceRecord(timestamp_unix=int(time.time()),
                            software="ttio imzml importer v0.9",
                            parameters=params)
    run = WrittenRun(
        spectrum_class="TTIOMassSpectrum",
        acquisition_mode=0,
        channel_data={"mz": mz, "intensity": intensity},
        offsets=np.arange(n, dtype=np.uint64) * 3,
        lengths=np.full(n, 3, dtype=np.uint32),
        retention_times=np.zeros(n),
        ms_levels=np.ones(n, dtype=np.int32),
        polarities=np.zeros(n, dtype=np.int32),
        precursor_mzs=np.zeros(n),
        precursor_charges=np.zeros(n, dtype=np.int32),
        base_peak_intensities=intensity.reshape(n, 3).max(axis=1),
        provenance_records=[prov] if run_level else [],
    )
    return SpectralDataset.write_minimal(
        path, title="legacy", isa_investigation_id="",
        runs={"imzml_pixels": run}, provenance=[prov])


@pytest.mark.parametrize("run_level", [True, False])
def test_legacy_csv_still_reads_and_exports(tmp_path: Path, run_level: bool) -> None:
    coords = [(1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1)]
    tio = _legacy_tio(tmp_path / "legacy.tio", coords, run_level=run_level)
    with SpectralDataset.open(tio) as ds:
        assert not ds.feature_flags.has("opt_pixel_coordinates")
        run = ds.ms_runs["imzml_pixels"]
        assert not run.index.has_pixel_coordinates
        assert run.index.pixel_coordinates_at(0) is None
        got = run_pixel_coordinates(run, ds.provenance())
        np.testing.assert_array_equal(got, np.array(coords, dtype=np.int32))
        if not run_level:
            # The dataset-level record is only consulted when passed in.
            assert run_pixel_coordinates(run) is None

    out = tmp_path / "legacy.imzML"
    export_registry.export("imzml", str(tio), None, str(out))
    back = imzml_reader.read(out)
    assert [(p.x, p.y, p.z) for p in back.spectra] == coords
    assert back.uuid_hex == "0123456789abcdef0123456789abcdef"
    assert back.mode == "continuous"


def test_legacy_csv_must_match_the_spectrum_count(tmp_path: Path) -> None:
    tio = _legacy_tio(tmp_path / "bad.tio", [(1, 1, 1), (2, 1, 1)],
                      run_level=True, csv="1,1,1")
    with SpectralDataset.open(tio) as ds:
        assert run_pixel_coordinates(ds.ms_runs["imzml_pixels"],
                                     ds.provenance()) is None
    parse = imzml_reader._parse_legacy_coordinates
    assert parse("1,1;2,1", 2) is None
    assert parse("a,b,c;1,1,1", 2) is None
    np.testing.assert_array_equal(parse("1,2,3;4,5,6", 2), [[1, 2, 3], [4, 5, 6]])


def test_partial_pixel_columns_rejected(tmp_path: Path) -> None:
    src = _pixels(3, 2)
    src_imzml = tmp_path / "s.imzML"
    imzml_writer.write(src, src_imzml, mode="processed")
    draft = imzml_reader.read(src_imzml).to_imported_dataset()
    draft.runs["imzml_pixels"].pixel_z = None
    with pytest.raises(ValueError, match="pixel_x/pixel_y/pixel_z"):
        draft.write(tmp_path / "partial.tio")

    tio = imzml_reader.read(src_imzml).to_ttio(tmp_path / "ok.tio")
    with h5py.File(tio, "r+") as f:
        del f["study/ms_runs/imzml_pixels/spectrum_index/pixel_z"]
    with pytest.raises(ValueError, match="partial pixel"):
        with SpectralDataset.open(tio) as ds:
            ds.ms_runs["imzml_pixels"].index  # noqa: B018


def test_runs_without_pixels_carry_no_flag(tmp_path: Path) -> None:
    tio = _legacy_tio(tmp_path / "plain.tio", [(1, 1, 1)], run_level=False)
    with h5py.File(tio, "r") as f:
        assert "pixel_x" not in f["study/ms_runs/imzml_pixels/spectrum_index"]
    with SpectralDataset.open(tio) as ds:
        assert not ds.feature_flags.has("opt_pixel_coordinates")
