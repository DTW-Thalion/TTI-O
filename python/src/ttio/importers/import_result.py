"""``ImportResult`` — lightweight in-memory container produced by importers.

The main ``SpectralDataset`` class wraps an open HDF5 file, which would be
awkward to construct from an importer that has no backing file yet. Instead,
importers produce an :class:`ImportResult` that can be inspected in memory
and then flushed to a real ``.tio`` file with :meth:`ImportResult.to_ttio`.

Notes
-----
API status: Stable (Python-idiomatic helper; ObjC and Java use
inline constructs).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from ..identification import Identification
from ..provenance import ProvenanceRecord
from ..quantification import Quantification
from ..spectral_dataset import SpectralDataset, WrittenRun
from .._threads import resolve_import_threads


@dataclass(slots=True)
class ImportedSpectrum:
    """One decoded spectrum from an import-time XML parse."""

    mz_or_chemical_shift: np.ndarray
    intensity: np.ndarray
    retention_time: float = 0.0
    ms_level: int = 1
    polarity: int = 0  # Polarity enum value
    precursor_mz: float = 0.0
    precursor_charge: int = 0
    # MS/MS activation + isolation. `activation_method` is an
    # ActivationMethod IntEnum value (NONE=0 for MS1 or unreported).
    # `isolation_*` are raw mzML fields; zero sentinel means MS1.
    activation_method: int = 0
    isolation_target_mz: float = 0.0
    isolation_lower_offset: float = 0.0
    isolation_upper_offset: float = 0.0


@dataclass(slots=True)
class ImportedChromatogram:
    """One decoded chromatogram trace from mzML (M24)."""

    retention_times: np.ndarray
    intensities: np.ndarray
    chromatogram_type: int = 0  # ChromatogramType IntEnum value
    target_mz: float = 0.0
    precursor_mz: float = 0.0
    product_mz: float = 0.0


@dataclass(slots=True)
class ImportResult:
    """Container returned by the mzML / nmrML importers."""

    title: str = ""
    isa_investigation_id: str = ""
    ms_spectra: list[ImportedSpectrum] = field(default_factory=list)
    nmr_spectra: list[ImportedSpectrum] = field(default_factory=list)
    chromatograms: list[ImportedChromatogram] = field(default_factory=list)
    nucleus_type: str = ""
    identifications: list[Identification] = field(default_factory=list)
    quantifications: list[Quantification] = field(default_factory=list)
    provenance: list[ProvenanceRecord] = field(default_factory=list)
    source_file: str = ""
    # nmrML acquisition-parameter fields. Populated by the nmrML
    # reader when present in the source file; zero / empty when the
    # file does not declare them. Cross-language convention: these
    # belong on the import result rather than on each ImportedSpectrum
    # because they are run-level metadata, not per-spectrum.
    spectrometer_frequency_mhz: float = 0.0
    number_of_scans: int = 0
    # nmrML free induction decay arrays. When the source file carries
    # ``<fidData byteFormat="complex128">`` the reader decodes the
    # interleaved real / imaginary float64 stream and exposes the two
    # halves here. Zero-length arrays when no FID is present.
    fid_real: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    fid_imag: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )

    def __iter__(self) -> Iterator[ImportedSpectrum]:
        yield from self.ms_spectra
        yield from self.nmr_spectra

    @property
    def spectrum_count(self) -> int:
        return len(self.ms_spectra) + len(self.nmr_spectra)

    def build_runs(self) -> dict[str, WrittenRun]:
        """Convert the parsed spectra into ``WrittenRun`` buffers ready for
        :func:`SpectralDataset.write_minimal`.
        """
        from ..chromatogram import Chromatogram
        from ..enums import ChromatogramType
        from ..signal_array import SignalArray

        chrom_objs = [
            Chromatogram(
                signal_arrays={
                    "time": SignalArray(data=c.retention_times),
                    "intensity": SignalArray(data=c.intensities),
                },
                axes=[],
                chromatogram_type=ChromatogramType(c.chromatogram_type),
                target_mz=c.target_mz,
                precursor_mz=c.precursor_mz,
                product_mz=c.product_mz,
            )
            for c in self.chromatograms
        ]

        runs: dict[str, WrittenRun] = {}
        if self.ms_spectra:
            runs["run_0001"] = _pack_run(
                self.ms_spectra, spectrum_class="TTIOMassSpectrum",
                acquisition_mode=0, channel_x="mz",
                chromatograms=chrom_objs,
            )
        elif chrom_objs:
            # Chromatograms with no spectra — synthesize an empty run to
            # carry them so the /chromatograms/ group has somewhere to live.
            runs["run_0001"] = _empty_run_with_chromatograms(chrom_objs)
        if self.nmr_spectra:
            runs["nmr_run"] = _pack_run(
                self.nmr_spectra, spectrum_class="TTIONMRSpectrum",
                acquisition_mode=4, channel_x="chemical_shift",
                nucleus_type=self.nucleus_type,
            )
        return runs

    def to_imported_dataset(self) -> "ImportedDataset":
        from .imported_dataset import ImportedDataset
        return ImportedDataset(
            title=self.title or "imported",
            isa_investigation_id=self.isa_investigation_id,
            runs=self.build_runs(),
            identifications=list(self.identifications),
            quantifications=list(self.quantifications),
            provenance=list(self.provenance),
        )

    def to_ttio(
        self,
        path: str | Path,
        features: list[str] | None = None,
        *,
        provider: str = "hdf5",
    ) -> Path:
        """Persist the parsed result as a ``.tio`` container.

: ``provider`` selects the storage backend
        (``"hdf5"``, ``"memory"``, ``"sqlite"``, ``"zarr"``). Passed
        through to :meth:`SpectralDataset.write_minimal` unchanged.
        """
        return self.to_imported_dataset().write(path, features=features,
                                                provider=provider)


def _pack_run(
    spectra: list[ImportedSpectrum],
    *,
    spectrum_class: str,
    acquisition_mode: int,
    channel_x: str,
    nucleus_type: str = "",
    chromatograms: list | None = None,
) -> WrittenRun:
    n = len(spectra)
    lengths = np.array([s.mz_or_chemical_shift.shape[0] for s in spectra], dtype=np.uint32)
    offsets = np.zeros(n, dtype=np.uint64)
    if n > 0:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.uint64)

    total = int(lengths.sum())
    x_buf = np.empty(total, dtype=np.float64)
    i_buf = np.empty(total, dtype=np.float64)
    pos = 0
    for s, length in zip(spectra, lengths):
        ln = int(length)
        x_buf[pos:pos + ln] = s.mz_or_chemical_shift
        i_buf[pos:pos + ln] = s.intensity
        pos += ln

    def _col(attr: str, dtype: type) -> np.ndarray:
        return np.array([getattr(s, attr) for s in spectra], dtype=dtype)

    base_peaks = np.array(
        [float(np.max(s.intensity)) if s.intensity.size else 0.0 for s in spectra],
        dtype=np.float64,
    )

    # M74 schema-gating: emit the four optional spectrum_index columns when
    # at least one spectrum carries activation metadata or a non-zero
    # isolation offset. Writer will store zero sentinels for the rest
    # (MS1 / unreported), matching the "all-or-nothing" on-disk schema.
    any_m74 = any(
        s.activation_method != 0
        or s.isolation_target_mz != 0.0
        or s.isolation_lower_offset != 0.0
        or s.isolation_upper_offset != 0.0
        for s in spectra
    )
    if any_m74:
        activation_methods = _col("activation_method", np.int32)
        isolation_target_mzs = _col("isolation_target_mz", np.float64)
        isolation_lower_offsets = _col("isolation_lower_offset", np.float64)
        isolation_upper_offsets = _col("isolation_upper_offset", np.float64)
    else:
        activation_methods = None
        isolation_target_mzs = None
        isolation_lower_offsets = None
        isolation_upper_offsets = None

    return WrittenRun(
        spectrum_class=spectrum_class,
        acquisition_mode=acquisition_mode,
        channel_data={channel_x: x_buf, "intensity": i_buf},
        offsets=offsets,
        lengths=lengths,
        retention_times=_col("retention_time", np.float64),
        ms_levels=_col("ms_level", np.int32),
        polarities=_col("polarity", np.int32),
        precursor_mzs=_col("precursor_mz", np.float64),
        precursor_charges=_col("precursor_charge", np.int32),
        base_peak_intensities=base_peaks,
        activation_methods=activation_methods,
        isolation_target_mzs=isolation_target_mzs,
        isolation_lower_offsets=isolation_lower_offsets,
        isolation_upper_offsets=isolation_upper_offsets,
        nucleus_type=nucleus_type,
        chromatograms=list(chromatograms or []),
    )


def _empty_run_with_chromatograms(chromatograms: list) -> WrittenRun:
    """Build an empty MS run that carries only chromatograms."""
    z8 = np.zeros(0, dtype=np.uint64)
    z4 = np.zeros(0, dtype=np.uint32)
    f8 = np.zeros(0, dtype=np.float64)
    i4 = np.zeros(0, dtype=np.int32)
    return WrittenRun(
        spectrum_class="TTIOMassSpectrum",
        acquisition_mode=0,
        channel_data={"mz": f8, "intensity": f8},
        offsets=z8, lengths=z4,
        retention_times=f8, ms_levels=i4, polarities=i4,
        precursor_mzs=f8, precursor_charges=i4, base_peak_intensities=f8,
        chromatograms=list(chromatograms),
    )


@dataclass
class GenomicStreamSource:
    """A genomic run delivered as consecutive WrittenGenomicRun batches
    (see BamReader.iter_batches, FastqReader.iter_batches). Written by
    ImportedDataset.write through GenomicStreamWriter with bounded
    memory. ``reference_fasta`` enables REF_DIFF_V2 via LazyReference."""

    name: str
    iter_batches: object                     # () -> Iterator[WrittenGenomicRun]
    reference_fasta: object | None = None    # Path or None
    embed_reference: bool = False
    block_reads: int | None = None           # None -> writer defaults
    block_bytes: int | None = None
    opt_legacy_whole_channel: bool = False
    # M103: group the reads by sequence before blocking (no reference).
    group_reads: bool = False
    threads: int | None = None                # None -> TTIO_THREADS
    # M97 — @read_role run attribute; falls back to the first batch's
    # read_role when None.
    read_role: str | None = None
    # M97 — REF_DIFF_V2 slice byte budget; 0 = the fixed-count rule.
    ref_diff_slice_bytes: int = 0

    def write_into(self, study_group, *, progress=None) -> int:
        from ..genomic.stream_writer import GenomicStreamWriter
        from ..genomic.lazy_reference import LazyReference
        ref = LazyReference(self.reference_fasta) if self.reference_fasta else None
        writer = None
        n = 0
        for batch in self.iter_batches():
            if writer is None:
                kw = {}
                if self.block_reads is not None:
                    kw["block_reads"] = self.block_reads
                if self.block_bytes is not None:
                    kw["block_bytes"] = self.block_bytes
                writer = GenomicStreamWriter(
                    study_group, self.name,
                    acquisition_mode=batch.acquisition_mode,
                    reference_uri=batch.reference_uri, platform=batch.platform,
                    sample_name=batch.sample_name,
                    reference_chrom_seqs=ref if ref is not None else batch.reference_chrom_seqs,
                    embed_reference=self.embed_reference or batch.embed_reference,
                    opt_disable_qualities_v5=batch.opt_disable_qualities_v5,
                    signal_codec_overrides=batch.signal_codec_overrides,
                    signal_compression=batch.signal_compression,
                    opt_legacy_whole_channel=self.opt_legacy_whole_channel,
                    group_reads=self.group_reads,
                    provenance_records=batch.provenance_records,
                    read_role=self.read_role or batch.read_role,
                    ref_diff_slice_bytes=(self.ref_diff_slice_bytes
                                          or batch.ref_diff_slice_bytes),
                    # Same rule as the producer side, so the two
                    # halves of the pipeline budget are sized off one
                    # count.
                    threads=(self.threads if self.threads is not None
                             else resolve_import_threads()), **kw)
            writer.append_batch(batch)
            n += int(len(batch.lengths))
        if writer is not None:
            writer.close()
        return n


@dataclass
class SpectralStreamSource:
    """A spectral run delivered as consecutive WrittenRun batches
    (see mzml.iter_batches). Written through SpectralStreamWriter."""

    name: str
    iter_batches: object                     # () -> Iterator[WrittenRun]
    instrument_config: object | None = None
    batch_spectra: int = 4096
    chromatograms_after: object | None = None   # () -> list[Chromatogram], read after the batches

    def write_into(self, study_group, *, progress=None) -> int:
        from ..spectral_stream_writer import SpectralStreamWriter
        writer = None
        n = 0
        for batch in self.iter_batches():
            if writer is None:
                writer = SpectralStreamWriter(
                    study_group, self.name, spectrum_class=batch.spectrum_class,
                    acquisition_mode=batch.acquisition_mode,
                    channel_names=list(batch.channel_data.keys()),
                    instrument_config=self.instrument_config,
                    batch_spectra=self.batch_spectra,
                    opt_disable_float_delta=batch.opt_disable_float_delta,
                    signal_compression=batch.signal_compression,
                    nucleus_type=batch.nucleus_type, solvent=getattr(batch, "solvent", ""),
                    provenance_records=batch.provenance_records)
            writer.append_batch(batch)
            n += int(batch.offsets.shape[0])
        if writer is not None:
            if self.chromatograms_after is not None:
                writer.set_chromatograms(self.chromatograms_after())
            writer.close()
        return n
