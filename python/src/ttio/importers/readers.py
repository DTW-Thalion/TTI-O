"""Per-format :class:`Reader` implementations. Each parses its source
into an :class:`ImportedDataset`; lazy module imports keep optional
dependencies off the import path.

The importer registry (:mod:`ttio.importers.registry`) holds one of these
per format and dispatches ``encode`` through ``Reader.read(...).write(...)``.
"""
from __future__ import annotations

from pathlib import Path

from .imported_dataset import ImportedDataset


class _ImportResultReader:
    """mzML / mzTab / nmrML / Thermo / Waters: module ``read()`` returns a
    result object exposing ``to_imported_dataset()``.

    ``_supports_progress`` is ``False`` for modules whose ``read`` does
    not accept a ``progress`` kwarg (Thermo, Waters); passing it would
    raise ``TypeError``, so the flag matches the old generic adapter,
    which only ever threaded ``progress`` for readers that accept it.
    """
    _module = ""
    _supports_progress = True

    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        import importlib
        mod = importlib.import_module(f"ttio.importers.{self._module}")
        kwargs = {}
        if self._supports_progress and progress is not None:
            kwargs["progress"] = progress
        return mod.read(inputs[0], **kwargs).to_imported_dataset()


class MzMLReader(_ImportResultReader):
    _module = "mzml"

    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from .mzml import MzMLStream
        stream = MzMLStream(inputs[0], batch_spectra=int(opts.get("batch_spectra", 4096)),
                            progress=progress)
        return ImportedDataset(title=Path(inputs[0]).stem,
                               spectral_streams={"run_0001": stream.stream_source()})


class MzTabReader(_ImportResultReader):
    _module = "mztab"


class NmrMLReader(_ImportResultReader):
    _module = "nmrml"


class ThermoRawReader(_ImportResultReader):
    _module = "thermo_raw"
    _supports_progress = False

    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from . import thermo_raw
        src = thermo_raw.stream_source(
            inputs[0], batch_spectra=int(opts.get("batch_spectra", 4096)), progress=progress)
        return ImportedDataset(title=Path(inputs[0]).stem, spectral_streams={"run_0001": src})


class WatersMassLynxReader(_ImportResultReader):
    _module = "waters_masslynx"
    _supports_progress = False

    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from . import waters_masslynx
        src = waters_masslynx.stream_source(
            inputs[0], batch_spectra=int(opts.get("batch_spectra", 4096)), progress=progress)
        stem = Path(inputs[0]).name
        stem = stem[:-4] if stem.lower().endswith(".raw") else stem
        return ImportedDataset(title=stem, spectral_streams={"run_0001": src})


class ImzMLReader:
    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from . import imzml
        ibd = opts.get("ibd")
        if ibd is None and len(inputs) > 1:
            ibd = inputs[1]
        kwargs = {"progress": progress} if progress is not None else {}
        return imzml.read(inputs[0], ibd_path=ibd, **kwargs).to_imported_dataset()


class BrukerReader:
    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from . import bruker_tdf
        stream = bruker_tdf.BrukerTDFStream(
            inputs[0], batch_frames=int(opts.get("batch_frames", 256)),
            ms2=bool(opts.get("ms2", False)))
        return ImportedDataset(title=Path(inputs[0]).stem,
                               spectral_streams=stream.stream_sources())


class JcampDxReader:
    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from pathlib import Path

        from . import jcamp_dx
        kwargs = {"progress": progress} if progress is not None else {}
        spectrum = jcamp_dx.read_spectrum(inputs[0], **kwargs)
        run = jcamp_dx.build_written_run(spectrum)
        return ImportedDataset(title=Path(inputs[0]).stem,
                               isa_investigation_id="",
                               runs={"run_0001": run})


class _GenomicReader:
    _attr = ""

    def read(self, inputs, opts, progress=None) -> ImportedDataset:
        from . import bam, cram, sam
        cls = {"BamReader": bam.BamReader, "SamReader": sam.SamReader,
               "CramReader": cram.CramReader}[self._attr]
        name = opts.get("name", "genomic_0001")
        kwargs = {"progress": progress} if progress is not None else {}
        reference = opts.get("reference")
        reader = cls(inputs[0], reference) if self._attr == "CramReader" else cls(inputs[0])
        src = reader.stream_source(
            name=name, sample_name=opts.get("sample"), reference_fasta=reference,
            embed_reference=bool(opts.get("embed_reference", False)),
            batch_reads=int(opts.get("batch_reads", 100_000)), **kwargs)
        src.block_reads = opts.get("block_reads")
        src.block_bytes = opts.get("block_bytes")
        src.opt_legacy_whole_channel = bool(opts.get("legacy_whole_channel", False))
        # M103; the cross-language CLIs pass "1"/"0" strings.
        src.group_reads = str(opts.get("group_reads", "0")).strip().lower() in ("1", "true", "yes")
        return ImportedDataset(genomic_streams={name: src})


class BamReader(_GenomicReader):
    _attr = "BamReader"


class SamReader(_GenomicReader):
    _attr = "SamReader"


class CramReader(_GenomicReader):
    _attr = "CramReader"
