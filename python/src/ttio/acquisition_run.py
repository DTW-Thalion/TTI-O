"""``AcquisitionRun`` — lazy view over one ``/study/ms_runs/<name>/`` group."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator

import json

import numpy as np

from . import _hdf5_io as io
from .access_policy import AccessPolicy
from .axis_descriptor import AxisDescriptor
from .chromatogram import Chromatogram
from .enums import (
    AcquisitionMode,
    ActivationMethod,
    ChromatogramType,
    EncryptionLevel,
    IRMode,
    Polarity,
    SpectrumKind,
)
from .instrument_config import InstrumentConfig
from .ir_spectrum import IRSpectrum
from .isolation_window import IsolationWindow
from .mass_spectrum import MassSpectrum
from .nmr_spectrum import NMRSpectrum
from .raman_spectrum import RamanSpectrum
from .uv_vis_spectrum import UVVisSpectrum
from .protocols import Indexable, Streamable, Provenanceable
from .provenance import ProvenanceRecord
from .providers.base import StorageDataset, StorageGroup
from .signal_array import SignalArray
from .spectrum import Spectrum
from .value_range import ValueRange


def _wrap_hdf5_group(obj: object) -> StorageGroup:
    """helper: adapt an h5py.Group (legacy caller) to a
    :class:`StorageGroup` view. Pass-through for StorageGroup inputs."""
    if isinstance(obj, StorageGroup):
        return obj
    # Assume h5py.Group.
    from .providers.hdf5 import _Group as _Hdf5Group
    return _Hdf5Group(obj)


# Channel -> default axis metadata for the two spectrum classes we currently
# materialize lazily. Writers may store additional channels; reading an unknown
# channel falls back to a generic "amplitude" axis.
_CHANNEL_AXIS: dict[str, AxisDescriptor] = {
    "mz": AxisDescriptor(name="mz", unit="m/z"),
    "intensity": AxisDescriptor(name="intensity", unit="counts"),
    "chemical_shift": AxisDescriptor(name="chemical_shift", unit="ppm"),
}


@dataclass(slots=True)
class SpectrumIndex:
    """Parallel per-spectrum arrays loaded eagerly at run open time.

    The arrays map 1:1 to the datasets under
    ``/study/ms_runs/<name>/spectrum_index/`` described in §4 of
    ``docs/format-spec.md``. They are small (length = spectrum_count) and
    cheap to hold in memory; signal channels remain lazy.

    Parameters
    ----------
    offsets : numpy.ndarray
        Starting element index of each spectrum in the mz_values channel.
    lengths : numpy.ndarray
        Number of elements (peaks) in each spectrum.
    retention_times : numpy.ndarray
        Retention time in seconds for each spectrum.
    ms_levels : numpy.ndarray
        MS level (1, 2, …) for each spectrum.
    polarities : numpy.ndarray
        Polarity (1=positive, -1=negative, 0=unknown) for each spectrum.
    precursor_mzs : numpy.ndarray
        Precursor m/z for each spectrum (0.0 for MS1).
    precursor_charges : numpy.ndarray
        Precursor charge state for each spectrum (0 for MS1).
    base_peak_intensities : numpy.ndarray
        Base-peak intensity for each spectrum.
    activation_methods : numpy.ndarray | None
        (optional) Activation method int value per spectrum; see
        :class:`~ttio.enums.ActivationMethod`. Present only when the
        file was written with the ``opt_ms2_activation_detail`` feature
        flag; absent columns indicate MS1-only or legacy files.
    isolation_target_mzs : numpy.ndarray | None
        (optional) Target m/z of the isolation window per spectrum.
        Zero when no isolation applied.
    isolation_lower_offsets : numpy.ndarray | None
        (optional) Lower offset of the isolation window per
        spectrum. Zero when no isolation applied.
    isolation_upper_offsets : numpy.ndarray | None
        (optional) Upper offset of the isolation window per
        spectrum. Zero when no isolation applied.
    pixel_x, pixel_y, pixel_z : numpy.ndarray | None
        (optional) Imaging-grid position of each spectrum (int32),
        present together on runs written with the
        ``opt_pixel_coordinates`` feature flag; ``None`` otherwise.

    Notes
    -----
    API status: Stable.

    Cross-language equivalents
    --------------------------
    Objective-C: ``TTIOSpectrumIndex`` ·
    Java: ``global.thalion.ttio.SpectrumIndex``.
    """

    offsets: np.ndarray
    lengths: np.ndarray
    retention_times: np.ndarray
    ms_levels: np.ndarray
    polarities: np.ndarray
    precursor_mzs: np.ndarray
    precursor_charges: np.ndarray
    base_peak_intensities: np.ndarray
    activation_methods: np.ndarray | None = None
    isolation_target_mzs: np.ndarray | None = None
    isolation_lower_offsets: np.ndarray | None = None
    isolation_upper_offsets: np.ndarray | None = None
    # Optional per-spectrum centroided flag. ``None`` for files written
    # before the centroided column was added; otherwise ``0 = profile,
    # 1 = centroided`` per spectrum. Mirrors mzML CV terms MS:1000127
    # (centroid) and MS:1000128 (profile).
    centroideds: np.ndarray | None = None
    # M102 imaging-grid position per spectrum (format-spec §4b). All
    # three or none.
    pixel_x: np.ndarray | None = None
    pixel_y: np.ndarray | None = None
    pixel_z: np.ndarray | None = None

    @property
    def count(self) -> int:
        """Number of spectra indexed.

        Returns
        -------
        int
            Length of the parallel ``offsets`` array; equals the number
            of spectra resolvable through :meth:`offset_at` /
            :meth:`length_at`.
        """
        return int(self.offsets.shape[0])

    # ------------------------------------------------------------------ #
    # Element-at accessors                                                 #
    # ------------------------------------------------------------------ #

    def offset_at(self, index: int) -> int:
        """Return element offset of spectrum ``index`` in mz_values."""
        return int(self.offsets[index])

    def length_at(self, index: int) -> int:
        """Return element count (peaks) of spectrum ``index``."""
        return int(self.lengths[index])

    def retention_time_at(self, index: int) -> float:
        """Return retention time in seconds of spectrum ``index``."""
        return float(self.retention_times[index])

    def ms_level_at(self, index: int) -> int:
        """Return MS level of spectrum ``index``."""
        return int(self.ms_levels[index])

    def polarity_at(self, index: int) -> Polarity:
        """Return :class:`~ttio.enums.Polarity` of spectrum ``index``."""
        return Polarity(int(self.polarities[index]))

    def precursor_mz_at(self, index: int) -> float:
        """Return precursor m/z of spectrum ``index`` (0.0 for MS1)."""
        return float(self.precursor_mzs[index])

    def precursor_charge_at(self, index: int) -> int:
        """Return precursor charge state of spectrum ``index`` (0 for MS1)."""
        return int(self.precursor_charges[index])

    def base_peak_intensity_at(self, index: int) -> float:
        """Return base-peak intensity of spectrum ``index``."""
        return float(self.base_peak_intensities[index])

    def activation_method_at(self, index: int) -> ActivationMethod:
        """Return :class:`~ttio.enums.ActivationMethod` of spectrum
        ``index``. Returns ``ActivationMethod.NONE`` when the M74
        column is absent (legacy file or MS1-only run)."""
        if self.activation_methods is None:
            return ActivationMethod.NONE
        return ActivationMethod(int(self.activation_methods[index]))

    def isolation_window_at(self, index: int) -> IsolationWindow | None:
        """Return :class:`~ttio.isolation_window.IsolationWindow` of
        spectrum ``index``, or ``None`` when the M74 columns are absent
        or the stored target+offsets are all zero (MS1 sentinel)."""
        if (self.isolation_target_mzs is None
                or self.isolation_lower_offsets is None
                or self.isolation_upper_offsets is None):
            return None
        target = float(self.isolation_target_mzs[index])
        lower = float(self.isolation_lower_offsets[index])
        upper = float(self.isolation_upper_offsets[index])
        if target == 0.0 and lower == 0.0 and upper == 0.0:
            return None
        return IsolationWindow(target_mz=target,
                                lower_offset=lower,
                                upper_offset=upper)

    def centroided_at(self, index: int) -> bool:
        """Return whether spectrum ``index`` is centroided (mzML
        MS:1000127). Returns ``False`` when the column is absent —
        callers needing to distinguish "unknown" from "profile" should
        check ``self.centroideds is None``.

        Cross-language equivalents
        --------------------------
        Java: ``SpectrumIndex.centroidedAt(i)`` · Objective-C:
        ``-[TTIOSpectrumIndex centroidedAt:]``.
        """
        if self.centroideds is None:
            return False
        return bool(self.centroideds[index])

    @property
    def has_pixel_coordinates(self) -> bool:
        """Whether the index carries the ``pixel_x/y/z`` columns."""
        return self.pixel_x is not None

    def pixel_coordinates_at(self, index: int) -> tuple[int, int, int] | None:
        """Return the ``(x, y, z)`` imaging-grid position of spectrum
        ``index``, or ``None`` when the run has no pixel columns.

        Cross-language equivalents
        --------------------------
        Java: ``SpectrumIndex.pixelCoordinatesAt(i)`` · Objective-C:
        ``-[TTIOSpectrumIndex pixelXAt:]`` / ``pixelYAt:`` / ``pixelZAt:``.
        """
        if self.pixel_x is None or self.pixel_y is None or self.pixel_z is None:
            return None
        return (int(self.pixel_x[index]), int(self.pixel_y[index]),
                int(self.pixel_z[index]))

    # ------------------------------------------------------------------ #
    # Range queries                                                        #
    # ------------------------------------------------------------------ #

    def indices_in_retention_time_range(self, value_range: ValueRange) -> list[int]:
        """Return indices whose retention time lies in
        ``[value_range.minimum, value_range.maximum]`` (inclusive)."""
        rt = self.retention_times
        mask = (rt >= value_range.minimum) & (rt <= value_range.maximum)
        return np.where(mask)[0].tolist()

    def indices_for_ms_level(self, ms_level: int) -> list[int]:
        """Return indices whose ``ms_level`` equals ``ms_level``."""
        return np.where(self.ms_levels == ms_level)[0].tolist()

    @classmethod
    def read(cls, idx_group: StorageGroup | "h5py.Group") -> "SpectrumIndex":
        """Load all index columns.

        Accepts either a :class:`StorageGroup` (protocol path)
        or an :class:`h5py.Group` (legacy). The legacy h5py path uses
        direct ``group[name]`` indexing so existing callers that
        dereferenced child groups continue to work.
        """
        if isinstance(idx_group, StorageGroup):
            # Protocol path: read columns via StorageDataset.read().
            def col(name: str, dtype: str) -> np.ndarray:
                ds = idx_group.open_dataset(name)
                arr = ds.read()
                return np.asarray(arr).astype(dtype, copy=False)

            def present(name: str) -> bool:
                return idx_group.has_child(name)
        else:
            # Legacy h5py.Group path.
            def col(name: str, dtype: str) -> np.ndarray:
                return idx_group[name][()].astype(dtype, copy=False)

            def present(name: str) -> bool:
                return name in idx_group

        # M74 schema-gating: load the four optional parallel columns only
        # when the writer emitted them (gated by opt_ms2_activation_detail).
        activation_methods = (
            col("activation_methods", "<i4")
            if present("activation_methods") else None
        )
        isolation_target_mzs = (
            col("isolation_target_mzs", "<f8")
            if present("isolation_target_mzs") else None
        )
        isolation_lower_offsets = (
            col("isolation_lower_offsets", "<f8")
            if present("isolation_lower_offsets") else None
        )
        isolation_upper_offsets = (
            col("isolation_upper_offsets", "<f8")
            if present("isolation_upper_offsets") else None
        )

        # Optional centroided column — independent of M74 gating.
        centroideds = col("centroideds", "<i4") if present("centroideds") else None

        # M102 pixel columns: all three or none.
        pixel_present = [present(n) for n in ("pixel_x", "pixel_y", "pixel_z")]
        if any(pixel_present) and not all(pixel_present):
            raise ValueError(
                "spectrum_index is malformed: partial pixel_x/pixel_y/pixel_z "
                "columns present")
        if all(pixel_present):
            pixel_x = col("pixel_x", "<i4")
            pixel_y = col("pixel_y", "<i4")
            pixel_z = col("pixel_z", "<i4")
        else:
            pixel_x = pixel_y = pixel_z = None

        # offsets is never on disk — synthesize from cumsum(lengths).
        from .genomic_index import _offsets_from_lengths
        lengths = col("lengths", "<u4")
        offsets = _offsets_from_lengths(lengths)
        return cls(
            offsets=offsets,
            lengths=lengths,
            retention_times=col("retention_times", "<f8"),
            ms_levels=col("ms_levels", "<i4"),
            polarities=col("polarities", "<i4"),
            precursor_mzs=col("precursor_mzs", "<f8"),
            precursor_charges=col("precursor_charges", "<i4"),
            base_peak_intensities=col("base_peak_intensities", "<f8"),
            activation_methods=activation_methods,
            isolation_target_mzs=isolation_target_mzs,
            isolation_lower_offsets=isolation_lower_offsets,
            isolation_upper_offsets=isolation_upper_offsets,
            centroideds=centroideds,
            pixel_x=pixel_x,
            pixel_y=pixel_y,
            pixel_z=pixel_z,
        )


@dataclass(slots=True)
class AcquisitionRun:
    """Lazy view over one acquisition run inside an ``.tio`` file.

    Conforms to :class:`~ttio.protocols.Indexable`,
    :class:`~ttio.protocols.Streamable`, and
    :class:`~ttio.protocols.Provenanceable`.
    :class:`~ttio.protocols.Encryptable` conformance is delivered
    in slice 41.5 when the encryption manager subsystem lands.

    Spectrum access is zero-copy-aware: the spectrum index is
    pre-loaded into numpy arrays at open time but signal channels are
    sliced on demand, so random access to spectrum ``i`` touches only
    the dataset bytes it needs.

    Notes
    -----
    API status: Stable (Encryptable surface pending).

    Cross-language equivalents
    --------------------------
    Objective-C: ``TTIOAcquisitionRun`` · Java:
    ``global.thalion.ttio.AcquisitionRun``.
    """

    name: str
    group: StorageGroup  # was h5py.Group; wrap legacy callers via _wrap_hdf5_group.
    spectrum_class: str
    acquisition_mode: AcquisitionMode
    index: SpectrumIndex
    channel_names: tuple[str, ...]
    instrument_config: InstrumentConfig
    nucleus_type: str = ""
    provenance_json: str = ""
    # top-level modality string. Read-side only — write-side
    # support arrives with GenomicRun in M74. Defaults to mass-spec for
    # backward-compat with files written by v0.10 and earlier (which
    # didn't carry the attribute at all).
    modality: str = "mass_spectrometry"
    # Optional NMR solvent label (e.g. "CDCl3", "DMSO-d6"). Empty when
    # not specified or when the run is not NMR. Stored as the
    # ``@solvent`` string attribute on the run group. Reused as the
    # UV-Vis solvent label.
    solvent: str = ""
    # Vibrational-spectrum run metadata (IR / Raman / UV-Vis), read from
    # scalar run-group attributes and consumed by _materialize_spectrum.
    # Defaults match the spectrum subclasses for files that don't carry
    # the attributes (i.e. all MS/NMR files).
    ir_mode: int = int(IRMode.TRANSMITTANCE)
    ir_resolution_cm_inv: float = 0.0
    ir_number_of_scans: int = 0
    raman_excitation_wavelength_nm: float = 0.0
    raman_laser_power_mw: float = 0.0
    raman_integration_time_sec: float = 0.0
    uvvis_path_length_cm: float = 0.0
    # chromatogram traces. Empty list on v0.3 files (group absent).
    chromatograms: list[Chromatogram] = field(default_factory=list)
    # signal cache holds protocol datasets, not h5py.Dataset.
    _signal_cache: dict[str, StorageDataset] = field(default_factory=dict, repr=False)
    # eagerly decoded Numpress-delta channels, keyed by channel
    # name. When present, :meth:`_materialize_spectrum` slices from
    # this float64 buffer instead of hitting the HDF5 dataset, because
    # Numpress decoding needs the running-sum prefix of the run.
    _numpress_channels: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    # FLOAT_DELTA_ZSTD (codec 17) channels are decoded per block, not
    # whole-channel: the FDZ1 block table per channel, a one-block
    # cache (block number, values) per channel, and a decode counter.
    _fdz_channels: set = field(default_factory=set, repr=False)
    _fdz_tables: dict = field(default_factory=dict, repr=False)
    _fdz_cache: dict = field(default_factory=dict, repr=False)
    _fdz_blocks_decoded: dict = field(default_factory=dict, repr=False)
    # M5-handoff: decrypted in-memory channels populated by
    # :meth:`decrypt_with_key`, keyed by channel name. When present,
    # :meth:`_materialize_spectrum` slices from this buffer so spectra
    # are readable through the normal API after decrypt while the
    # on-disk file stays encrypted. Mirrors ObjC rehydration.
    _decrypted_channels: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    # Lazy full-column cache for the plain (non-decrypted, non-numpress)
    # read path. On first access of a channel, :meth:`_materialize_spectrum`
    # bulk-reads the WHOLE channel column once and stores it here, then
    # slices in memory for every spectrum — mirroring the ObjC/Java fast
    # path and the existing _decrypted/_numpress caches. The underlying
    # run is read-only after open, so no invalidation is needed.
    _full_channel_cache: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    # Streamable cursor and Provenanceable cache.
    _cursor: int = field(default=0, repr=False)
    _provenance_cache: list[ProvenanceRecord] | None = field(default=None, repr=False)
    # Encryptable conformance.
    _access_policy: AccessPolicy | None = field(default=None, repr=False)
    # persistence context — set by SpectralDataset.open so that
    # encrypt_with_key / decrypt_with_key can delegate to the encryption module.
    _persistence_file_path: str | None = field(default=None, repr=False)
    _persistence_run_name: str | None = field(default=None, repr=False)
    # When True (local files), :meth:`_materialize_spectrum` bulk-reads
    # each channel column once and slices in memory (the fast path). When
    # False (remote fsspec-backed runs), it keeps the per-spectrum
    # hyperslab read so a single random access does not pull the whole
    # column over the network — preserving lazy HTTP range-request
    # behavior. Set by SpectralDataset._from_provider.
    _bulk_read: bool = field(default=True, repr=False)

    @property
    def kind(self) -> SpectrumKind:
        """The :class:`~ttio.enums.SpectrumKind` derived from the persisted
        ``@spectrum_class`` string.

        In-code dispatch key only (P3.8). The ``spectrum_class`` field stays
        the on-disk source of truth and is never normalized through this
        enum; an unrecognized class maps to
        :attr:`~ttio.enums.SpectrumKind.UNKNOWN` and dispatches to the
        :class:`~ttio.mass_spectrum.MassSpectrum` default, preserving the
        v0.1 fallback behavior.
        """
        return SpectrumKind.from_persisted(self.spectrum_class)

    @classmethod
    def open(
        cls,
        group: StorageGroup | "h5py.Group",
        name: str,
        *,
        bulk_read: bool = True,
    ) -> "AcquisitionRun":
        """Open a run from a storage-provider group .

        Accepts either a :class:`StorageGroup` (protocol path) or an
        :class:`h5py.Group` (legacy path — wrapped transparently via
        :class:`~ttio.providers.hdf5._Group`). made every
        cold-path read route through the StorageGroup protocol so
        Memory / SQLite / Zarr backends work without touching h5py.

        ``bulk_read`` (default ``True``) controls the per-spectrum read
        strategy: local files bulk-read each channel column once and
        slice in memory; pass ``False`` for remote fsspec-backed runs so
        a single random access only fetches its own slice. The materialized
        spectra are byte-identical either way.
        """
        sgroup = _wrap_hdf5_group(group)

        mode_raw = io.read_int_attr(sgroup, "acquisition_mode", default=0) or 0
        spectrum_class = io.read_string_attr(
            sgroup, "spectrum_class", default="TTIOMassSpectrum"
        ) or "TTIOMassSpectrum"
        nucleus = io.read_string_attr(sgroup, "nucleus_type", default="") or ""
        prov = io.read_string_attr(sgroup, "provenance_json", default="") or ""
        modality = io.read_string_attr(
            sgroup, "modality", default="mass_spectrometry"
        ) or "mass_spectrometry"
        solvent = io.read_string_attr(sgroup, "solvent", default="") or ""

        # Vibrational run metadata (absent on MS/NMR files -> defaults).
        ir_mode = io.read_int_attr(
            sgroup, "ir_mode", default=int(IRMode.TRANSMITTANCE)
        )
        ir_resolution_cm_inv = io.read_float_attr(
            sgroup, "ir_resolution_cm_inv", default=0.0) or 0.0
        ir_number_of_scans = io.read_int_attr(
            sgroup, "ir_number_of_scans", default=0) or 0
        raman_excitation_wavelength_nm = io.read_float_attr(
            sgroup, "raman_excitation_wavelength_nm", default=0.0) or 0.0
        raman_laser_power_mw = io.read_float_attr(
            sgroup, "raman_laser_power_mw", default=0.0) or 0.0
        raman_integration_time_sec = io.read_float_attr(
            sgroup, "raman_integration_time_sec", default=0.0) or 0.0
        uvvis_path_length_cm = io.read_float_attr(
            sgroup, "uvvis_path_length_cm", default=0.0) or 0.0

        idx = SpectrumIndex.read(sgroup.open_group("spectrum_index"))

        sig_group = sgroup.open_group("signal_channels")
        channel_names_raw = io.read_string_attr(
            sig_group, "channel_names", default="mz,intensity"
        ) or "mz,intensity"
        channel_names = tuple(c for c in channel_names_raw.split(",") if c)

        if sgroup.has_child("instrument_config"):
            cfg_group = sgroup.open_group("instrument_config")
            config = InstrumentConfig(
                manufacturer=io.read_string_attr(cfg_group, "manufacturer", "") or "",
                model=io.read_string_attr(cfg_group, "model", "") or "",
                serial_number=io.read_string_attr(cfg_group, "serial_number", "") or "",
                source_type=io.read_string_attr(cfg_group, "source_type", "") or "",
                analyzer_type=io.read_string_attr(cfg_group, "analyzer_type", "") or "",
                detector_type=io.read_string_attr(cfg_group, "detector_type", "") or "",
            )
        else:
            config = InstrumentConfig()

        # detect Numpress-delta channels via the
        # ``<chName>_numpress_fixed_point`` attribute on the
        # signal_channels group, and eagerly decode them here so
        # :meth:`_materialize_spectrum` can just slice a float64 buffer.
        numpress_channels: dict[str, np.ndarray] = {}
        fdz_channels: set = set()
        for chName in channel_names:
            scale_attr = f"{chName}_numpress_fixed_point"
            if sig_group.has_attribute(scale_attr):
                from ._numpress import decode as _np_decode
                ds_name = f"{chName}_values"
                if not sig_group.has_child(ds_name):
                    continue
                raw = sig_group.open_dataset(ds_name).read()
                scale = int(sig_group.get_attribute(scale_attr))
                numpress_channels[chName] = _np_decode(raw, scale)
                continue
            # FLOAT_DELTA_ZSTD (codec id 17): the channel dataset is a
            # flat uint8 FDZ1 stream with @compression = 17 on the
            # dataset. Decode-once-cache into the same eager buffer the
            # numpress path uses; every downstream consumer (spectrum
            # materialisation, transport writer) already prefers it.
            ds_name = f"{chName}_values"
            if not sig_group.has_child(ds_name):
                continue
            ds = sig_group.open_dataset(ds_name)
            from . import _hdf5_io as _io
            codec_id = _io.read_int_attr(ds, "compression", default=0) or 0
            if codec_id == 17:
                fdz_channels.add(chName)
            elif codec_id != 0:
                raise ValueError(
                    f"signal channel {chName!r}: @compression={codec_id} "
                    "is not a spectral channel codec (FLOAT_DELTA_ZSTD=17 "
                    "is the only one wired)"
                )

        return cls(
            name=name,
            group=sgroup,  # StorageGroup protocol value.
            spectrum_class=spectrum_class,
            acquisition_mode=AcquisitionMode(mode_raw),
            index=idx,
            channel_names=channel_names,
            instrument_config=config,
            nucleus_type=nucleus,
            provenance_json=prov,
            modality=modality,
            solvent=solvent,
            ir_mode=int(ir_mode if ir_mode is not None else IRMode.TRANSMITTANCE),
            ir_resolution_cm_inv=ir_resolution_cm_inv,
            ir_number_of_scans=ir_number_of_scans,
            raman_excitation_wavelength_nm=raman_excitation_wavelength_nm,
            raman_laser_power_mw=raman_laser_power_mw,
            raman_integration_time_sec=raman_integration_time_sec,
            uvvis_path_length_cm=uvvis_path_length_cm,
            chromatograms=_read_chromatograms(sgroup),
            _numpress_channels=numpress_channels,
            _fdz_channels=fdz_channels,
            _bulk_read=bulk_read,
        )

    # ----------------------------------------------------- spectrum access

    def __len__(self) -> int:
        return self.index.count

    def __iter__(self) -> Iterator[Spectrum]:
        return self.iter_spectra()

    def channel_range(self, channel: str, start: int, count: int, *,
                      threads: int | None = None) -> np.ndarray:
        """Values ``[start, start + count)`` of a signal channel as
        float64, touching only what is needed: decrypted or numpress
        buffers are sliced, codec 17 channels decode the blocks the
        range covers (one-block cache per channel), plain channels are
        read as a hyperslab (or sliced from the full-column cache when
        it is already loaded)."""
        if count <= 0:
            return np.zeros(0, dtype=np.float64)
        decrypted = self._decrypted_channels.get(channel)
        if decrypted is not None:
            return decrypted[start:start + count]
        decoded = self._numpress_channels.get(channel)
        if decoded is not None:
            return decoded[start:start + count]
        if channel in self._fdz_channels:
            return self._fdz_range(channel, start, count, threads)
        full = self._full_channel_cache.get(channel)
        if full is not None:
            return full[start:start + count]
        ds = self._signal_dataset(channel)
        return np.asarray(ds.read(offset=start, count=count))

    def _fdz_table(self, channel: str):
        table = self._fdz_tables.get(channel)
        if table is None:
            from .codecs import float_delta_zstd as _fdz
            ds = self._signal_dataset(channel)
            table = _fdz.read_block_table(
                lambda off, n: np.asarray(ds.read(offset=off, count=n), dtype=np.uint8).tobytes())
            self._fdz_tables[channel] = table
        return table

    def _fdz_block(self, channel: str, k: int) -> np.ndarray:
        cached = self._fdz_cache.get(channel)
        if cached is not None and cached[0] == k:
            return cached[1]
        from .codecs import float_delta_zstd as _fdz
        ds = self._signal_dataset(channel)
        table = self._fdz_table(channel)
        vals = _fdz.decode_block(
            lambda off, n: np.asarray(ds.read(offset=off, count=n), dtype=np.uint8).tobytes(),
            table, k)
        self._fdz_cache[channel] = (k, vals)
        self._fdz_blocks_decoded[channel] = self._fdz_blocks_decoded.get(channel, 0) + 1
        return vals

    def _fdz_range(self, channel: str, start: int, count: int,
                   threads: int | None = None) -> np.ndarray:
        if count <= 0:
            return np.zeros(0, dtype=np.float64)
        table = self._fdz_table(channel)
        bs = table.block_size
        k0, k1 = start // bs, (start + count - 1) // bs
        if k0 == k1:
            blk = self._fdz_block(channel, k0)
            o = start - k0 * bs
            return blk[o:o + count]
        from ._threads import resolve_threads, pool_context
        ks = list(range(k0, k1 + 1))
        nthreads = resolve_threads(threads)
        if nthreads > 1:
            # storage reads on this thread, decodes on the pool
            from .codecs import float_delta_zstd as _fdz
            import concurrent.futures as _cf
            ds = self._signal_dataset(channel)
            cached = self._fdz_cache.get(channel)
            todo = [k for k in ks if cached is None or cached[0] != k]
            raw = {k: np.asarray(ds.read(offset=int(table.offsets[k]), count=int(table.lengths[k])),
                                 dtype=np.uint8).tobytes() for k in todo}
            with pool_context(nthreads), _cf.ThreadPoolExecutor(
                    max_workers=nthreads, thread_name_prefix="ttio-fdz-decode") as pool:
                futs = {k: pool.submit(_fdz.decode_block_bytes, int(table.transforms[k]), raw[k],
                                       table.block_values(k)) for k in todo}
                blocks = {k: futs[k].result() for k in todo}
            if cached is not None and cached[0] in ks:
                blocks[cached[0]] = cached[1]
            self._fdz_cache[channel] = (k1, blocks[k1])
            self._fdz_blocks_decoded[channel] = self._fdz_blocks_decoded.get(channel, 0) + len(todo)
        else:
            blocks = {k: self._fdz_block(channel, k) for k in ks}
        parts = []
        for k in ks:
            blk = blocks[k]
            lo = start - k * bs if k == k0 else 0
            hi = (start + count) - k * bs if k == k1 else len(blk)
            parts.append(blk[lo:hi])
        return np.concatenate(parts)

    def iter_spectra(self, batch: int = 4096, *, threads: int | None = None) -> Iterator[Spectrum]:
        """Yield every spectrum in order, reading channel data in
        windows of ``batch`` spectra (bounded memory; codec 17 blocks
        decode once per window, on a pool of ``threads`` when a window
        spans several blocks)."""
        n = len(self)
        offsets, lengths = self.index.offsets, self.index.lengths
        for j0 in range(0, n, max(1, batch)):
            j1 = min(n, j0 + batch)
            base = int(offsets[j0])
            total = int(offsets[j1 - 1]) + int(lengths[j1 - 1]) - base
            arrays = {}
            for c in self.channel_names:
                try:
                    arrays[c] = self.channel_range(c, base, total, threads=threads)
                except KeyError:
                    continue
            for i in range(j0, j1):
                o = int(offsets[i]) - base
                ln = int(lengths[i])
                yield self._build_spectrum(i, {c: a[o:o + ln] for c, a in arrays.items()})

    # ── parallel block consumer ──────────────────────────────────────

    def _spectral_units(self, lo: int, hi: int) -> list[tuple]:
        """Scheduling units covering spectra ``[lo, hi)``.

        Spectral blocks are cut on value boundaries while a spectrum is
        located by value offset, so block edges fall inside spectra. A
        spectrum belongs to the block holding its **first** value, which
        partitions the spectra exactly once. A unit's value extent then
        runs to the end of its last spectrum and may lie past the owning
        block's end, and a block holding no spectrum start owns nothing
        and is skipped rather than emitted empty.

        Units come from the FDZ1 block tables when every channel has one
        and they agree, and from fixed spectrum counts otherwise, which
        covers numpress, decrypted and plain channels. Objective-C reads
        ``blocks/index`` as a shortcut for the first case; the block
        starts it records are the same ``b * block_size``, so the units
        and everything observable about them are identical.

        Returns ``(block, first_spectrum, n_spectra, value_start,
        value_end)`` tuples.
        """
        offsets, lengths = self.index.offsets, self.index.lengths
        units: list[tuple] = []
        if lo >= hi:
            return units

        block_size = n_blocks = 0
        names = list(self.channel_names)
        if names and all(c in self._fdz_channels for c in names):
            tables = [self._fdz_table(c) for c in names]
            first = tables[0]
            if all(t.block_size == first.block_size and t.n_blocks == first.n_blocks
                   for t in tables):
                block_size, n_blocks = int(first.block_size), int(first.n_blocks)

        if n_blocks:
            b, i = 0, lo
            while i < hi:
                while b + 1 < n_blocks and (b + 1) * block_size <= int(offsets[i]):
                    b += 1
                end = ((b + 1) * block_size) if b + 1 < n_blocks else (1 << 62)
                start_i = i
                while i < hi and int(offsets[i]) < end:
                    i += 1
                units.append((b, start_i, i - start_i, int(offsets[start_i]),
                              int(offsets[i - 1]) + int(lengths[i - 1])))
            return units

        batch = 512
        for i in range(lo, hi, batch):
            last = min(i + batch, hi) - 1
            units.append((len(units), i, last - i + 1, int(offsets[i]),
                          int(offsets[last]) + int(lengths[last])))
        return units

    def for_each_block(self, fn, start: int = 0, stop: int | None = None, *,
                       threads: int | None = None) -> None:
        """Call ``fn(view, view_start, first_spectrum, n_spectra)`` once
        per scheduling unit, on the pool.

        Read spectrum ``k`` of the delivered range as
        ``view[view_start + k]``; its index in the whole run is
        ``first_spectrum + k``. The two differ whenever a range starts
        part-way into a block, so a caller that indexes the view by
        ``first_spectrum`` reads the wrong spectra.

        ``fn`` is called from several threads at once and in no
        particular order and must be safe that way. ``view`` is a run
        over one unit's spectra, indexed from 0, valid only for the
        duration of the call, and it builds a spectrum only when one is
        asked for. That relaxed ordering is the whole difference from
        :meth:`iter_spectra`, whose in-order delivery on the calling
        thread is what bounds it.

        ⚠️ **In CPython this is close to flat.** Measured on an
        Orbitrap Exploris run, 342546 spectra over 16 units, MB/s:
        :meth:`iter_spectra` 97.7, this method 99.6 at one worker, 102.5
        at four, 100.4 at sixteen. So it buys about 5 per cent and does
        not degrade.

        That is *not* what the genomic equivalent measures, and the
        difference is the callback's work rather than the method: there
        the consumer builds read objects in pure Python and is
        GIL-bound, going 6.2 MB/s serial, 6.8 at four workers and 4.5 at
        eight and above. Here the consumer's work is numpy reductions
        over channel arrays, and numpy releases the GIL, so the
        consumers overlap instead of serialising. A callback that stays
        in pure Python should expect the genomic shape.

        The compiled SDKs gain far more: Objective-C measures 213.3 MB/s
        ordered against 780.5 at 16 threads on the same run.

        An unset ``threads`` still resolves to at most
        ``_read_ahead_blocks()`` workers rather than to ``TTIO_THREADS``.
        With unit counts in the tens there is nothing to gain past a few
        workers, so the cap costs nothing and protects a pure-Python
        callback from the genomic outcome. An explicit ``threads`` is
        honoured as given.

        Java: ``AcquisitionRun.iterBlocks``; Objective-C:
        ``-iterBlocksFrom:to:threads:error:usingBlock:``.
        """
        n = len(self)
        if stop is None or stop > n:
            stop = n
        lo = max(start + n if start < 0 else start, 0)
        if lo >= stop:
            return

        units = self._spectral_units(lo, stop)
        if not units:
            return

        from ._threads import resolve_threads, pool_context
        from .genomic_run import _read_ahead_blocks
        nthreads = (resolve_threads(threads) if threads is not None
                    else min(resolve_threads(None), _read_ahead_blocks()))
        if nthreads <= 1 or len(units) == 1:
            for u in units:
                view = _SpectralUnitView(self, u)
                fn(view, 0, u[1], u[2])
            return

        import concurrent.futures as _cf
        window = min(nthreads, len(units))
        with pool_context(window), _cf.ThreadPoolExecutor(
                max_workers=window, thread_name_prefix="ttio-sblock") as pool:
            def one(u: tuple) -> None:
                fn(_SpectralUnitView(self, u), 0, u[1], u[2])
            for fut in [pool.submit(one, u) for u in units]:
                fut.result()


    def spectra(self) -> list[Spectrum]:
        """Return all spectra as a list.

        Convenience for parity with the Java ``AcquisitionRun.spectra()``
        and Objective-C ``-[TTIOAcquisitionRun spectra]`` accessors.
        Equivalent to ``list(run)`` but provided as a named method so
        cross-language code reads the same.
        """
        return list(self)

    def provenance(self) -> list[ProvenanceRecord]:
        """Per-run provenance records.

        Prefers the v0.3 compound layout at ``<run>/provenance/steps`` and
        falls back to the v0.2 ``@provenance_json`` attribute. Pre-v0.2
        files (no per-run provenance of any kind) return an empty list.
        """
        # cold-path — navigate the compound-provenance subgroup through
        # the StorageGroup protocol. read_compound_dataset (via
        # _decode_provenance_compound) unwraps to the native h5py handle
        # for HDF5 backends and falls back to read_rows() otherwise.
        if self.group.has_child("provenance"):
            prov = self.group.open_group("provenance")
            if prov.has_child("steps"):
                return _decode_provenance_compound(prov, "steps")
        if self.provenance_json:
            return _decode_provenance_json(self.provenance_json)
        return []

    def __getitem__(self, i: int) -> Spectrum:
        if i < 0:
            i += len(self)
        if not 0 <= i < len(self):
            raise IndexError(f"spectrum index {i} out of range [0, {len(self)})")
        return self._materialize_spectrum(i)

    # ---- Indexable conformance ----

    def object_at_index(self, index: int) -> Spectrum:
        """Return the spectrum at ``index``. Negative indices are supported."""
        return self[index]

    def count(self) -> int:
        """Return the total number of spectra."""
        return len(self)

    def object_for_key(self, key: object) -> Spectrum:
        """Not supported — AcquisitionRun uses integer indexing only."""
        raise NotImplementedError("AcquisitionRun does not support key-based access")

    def objects_in_range(self, start: int, stop: int) -> list[Spectrum]:
        """Return spectra in the half-open slice ``[start, stop)``."""
        return [self[i] for i in range(start, stop)]

    # ---- Streamable conformance ----

    def next_object(self) -> Spectrum:
        """Return the next spectrum and advance the cursor."""
        if self._cursor >= len(self):
            raise StopIteration
        s = self[self._cursor]
        self._cursor += 1
        return s

    def has_more(self) -> bool:
        """Return ``True`` if ``next_object`` can be called."""
        return self._cursor < len(self)

    def current_position(self) -> int:
        """0-based position of the next spectrum to be yielded."""
        return self._cursor

    def seek_to_position(self, position: int) -> bool:
        """Reposition the cursor. Returns ``True`` on success."""
        if not 0 <= position <= len(self):
            return False
        self._cursor = position
        return True

    def reset(self) -> None:
        """Reposition the cursor to 0."""
        self._cursor = 0

    # ---- Provenanceable conformance ----

    def add_processing_step(self, step: ProvenanceRecord) -> None:
        """Append a processing step to this run's provenance chain."""
        if self._provenance_cache is None:
            self._provenance_cache = self.provenance()
        self._provenance_cache.append(step)

    def provenance_chain(self) -> list[ProvenanceRecord]:
        """Return this run's provenance records in insertion order."""
        if self._provenance_cache is not None:
            return list(self._provenance_cache)
        return self.provenance()

    def input_entities(self) -> list[str]:
        """Distinct input entity identifiers referenced by the chain."""
        seen: list[str] = []
        for r in self.provenance_chain():
            for e in r.input_refs:
                if e not in seen:
                    seen.append(e)
        return seen

    def output_entities(self) -> list[str]:
        """Distinct output entity identifiers referenced by the chain."""
        seen: list[str] = []
        for r in self.provenance_chain():
            for e in r.output_refs:
                if e not in seen:
                    seen.append(e)
        return seen

    # ---- Encryptable conformance ----

    def _set_persistence_context(self, file_path: str, run_name: str) -> None:
        """Attach file + run path so ``encrypt_with_key`` can delegate.

        Internal API — called by SpectralDataset._from_provider after
        loading each run.
        """
        object.__setattr__(self, "_persistence_file_path", file_path)
        object.__setattr__(self, "_persistence_run_name", run_name)

    def encrypt_with_key(self, key: bytes, level: EncryptionLevel) -> None:
        """Encrypt this run's intensity channel in place.

        Operates through the already-open HDF5 group so no second file
        handle is required — the file must be open for writing (``"r+"``
        or ``"w"``). Matches ObjC
        ``-[TTIOAcquisitionRun encryptWithKey:level:error:]`` semantics.

        Requires a persistence context — call only after opening via
        :meth:`SpectralDataset.open`.
        """
        if not self._persistence_file_path or not self._persistence_run_name:
            raise RuntimeError(
                "AcquisitionRun.encrypt_with_key requires a persistence "
                "context; call via a run obtained from SpectralDataset.open"
            )
        from .encryption import encrypt_intensity_channel_in_group
        # phase B: the in-place encrypter now accepts either
        # h5py.Group or StorageGroup, so route through the protocol —
        # works on Memory / SQLite / Zarr backends without unwrapping
        # to a raw h5py handle.
        sig_group = self.group.open_group("signal_channels")
        encrypt_intensity_channel_in_group(sig_group, key)

    def decrypt_with_key(self, key: bytes) -> bytes:
        """Decrypt this run's intensity channel.

        Returns the plaintext bytes. The on-disk file is NOT modified.
        Operates through the already-open HDF5 group. Matches ObjC
        ``-[TTIOAcquisitionRun decryptWithKey:]`` semantics (NSData → bytes).

        Side effect: the decrypted ndarray is cached on the run so
        subsequent :meth:`object_at_index` calls can return
        :class:`SignalArray` views without re-decrypting. See the M5
        handoff note — spectra become readable through the normal API
        path after decrypt.

        Requires a persistence context — call only after opening via
        :meth:`SpectralDataset.open`.
        """
        if not self._persistence_file_path or not self._persistence_run_name:
            raise RuntimeError(
                "AcquisitionRun.decrypt_with_key requires a persistence "
                "context; call via a run obtained from SpectralDataset.open"
            )
        from .encryption import read_encrypted_channel
        # phase B: read_encrypted_channel accepts both h5py
        # and StorageGroup; no native unwrap needed.
        sig_group = self.group.open_group("signal_channels")
        arr = read_encrypted_channel(sig_group, "intensity", key, dtype="<f8")
        self._decrypted_channels["intensity"] = arr
        return arr.tobytes()

    def access_policy(self) -> AccessPolicy | None:
        """Return the current access policy, or ``None`` if not set."""
        return self._access_policy

    def set_access_policy(self, policy: AccessPolicy | None) -> None:
        """Replace the current access policy."""
        object.__setattr__(self, "_access_policy", policy)

    def _signal_dataset(self, channel: str) -> StorageDataset:
        """Resolve ``channel`` to its ``<channel>_values`` StorageDataset
        (). Cached across calls; missing channels raise
        :class:`KeyError` so :meth:`_materialize_spectrum` can skip."""
        ds = self._signal_cache.get(channel)
        if ds is not None:
            return ds
        name = f"{channel}_values"
        sig_group = self.group.open_group("signal_channels")
        if not sig_group.has_child(name):
            raise KeyError(
                f"signal channel {channel!r} missing under run {self.name!r}")
        ds = sig_group.open_dataset(name)
        self._signal_cache[channel] = ds
        return ds

    def _materialize_spectrum(self, i: int) -> Spectrum:
        offset = int(self.index.offsets[i])
        length = int(self.index.lengths[i])

        arrays: dict[str, np.ndarray] = {}
        for c in self.channel_names:
            decrypted = self._decrypted_channels.get(c)
            if decrypted is not None:
                arr = decrypted[offset:offset + length]
            elif (decoded := self._numpress_channels.get(c)) is not None:
                arr = decoded[offset:offset + length]
            elif c in self._fdz_channels:
                arr = self._fdz_range(c, offset, length)
            else:
                full = self._full_channel_cache.get(c)
                if full is not None:
                    # Cache HIT: slice the already-loaded full column.
                    arr = full[offset:offset + length]
                elif self._bulk_read:
                    # Local fast path: bulk-read the WHOLE channel column
                    # once through the storage protocol (offset=0,
                    # count=-1 -> full read), cache it, then slice in
                    # memory for every spectrum. This mirrors the ObjC/Java
                    # fast path and the _decrypted/_numpress caches above.
                    # The per-spectrum slice ``full[offset:offset+length]``
                    # is byte-identical to the old per-call
                    # ``ds.read(offset=offset, count=length)`` (same dtype,
                    # same element range). Works uniformly across
                    # Hdf5/Memory/Sqlite providers — the cross-backend
                    # byte-identity tests in M43 guarantee equivalence.
                    try:
                        ds = self._signal_dataset(c)
                    except KeyError:
                        continue
                    full = np.asarray(ds.read())
                    self._full_channel_cache[c] = full
                    arr = full[offset:offset + length]
                else:
                    # Remote (fsspec) path: keep the per-spectrum hyperslab
                    # read so a single random access only pulls its own
                    # slice over the network, not the whole column. The
                    # result is identical to the bulk-then-slice path.
                    try:
                        ds = self._signal_dataset(c)
                    except KeyError:
                        continue
                    arr = np.asarray(ds.read(offset=offset, count=length))
            arrays[c] = arr
        return self._build_spectrum(i, arrays)

    def _build_spectrum(self, i: int, arrays: dict) -> Spectrum:
        signal_arrays: dict[str, SignalArray] = {}
        for c, arr in arrays.items():
            axis = _CHANNEL_AXIS.get(c, AxisDescriptor(name=c, unit=""))
            signal_arrays[c] = SignalArray.from_numpy(arr, axis=axis)

        polarity = Polarity(int(self.index.polarities[i]))
        base_kwargs = dict(
            signal_arrays=signal_arrays,
            index_position=i,
            scan_time_seconds=float(self.index.retention_times[i]),
            precursor_mz=float(self.index.precursor_mzs[i]),
            precursor_charge=int(self.index.precursor_charges[i]),
        )

        kind = self.kind  # one lookup per spectrum (hot path)
        if kind is SpectrumKind.NMR:
            return NMRSpectrum(nucleus_type=self.nucleus_type, **base_kwargs)
        if kind is SpectrumKind.IR:
            return IRSpectrum(
                mode=IRMode(self.ir_mode),
                resolution_cm_inv=self.ir_resolution_cm_inv,
                number_of_scans=self.ir_number_of_scans,
                **base_kwargs,
            )
        if kind is SpectrumKind.RAMAN:
            return RamanSpectrum(
                excitation_wavelength_nm=self.raman_excitation_wavelength_nm,
                laser_power_mw=self.raman_laser_power_mw,
                integration_time_sec=self.raman_integration_time_sec,
                **base_kwargs,
            )
        if kind is SpectrumKind.UVVIS:
            return UVVisSpectrum(
                path_length_cm=self.uvvis_path_length_cm,
                solvent=self.solvent,
                **base_kwargs,
            )
        return MassSpectrum(
            ms_level=int(self.index.ms_levels[i]),
            polarity=polarity,
            is_centroided=self.index.centroided_at(i),
            **base_kwargs,
        )


# ------------------------------------------------ provenance decoders ---


def _decode_provenance_compound(
    prov_group: "StorageGroup", dataset_name: str
) -> list[ProvenanceRecord]:
    records = io.read_compound_dataset(prov_group, dataset_name)
    out: list[ProvenanceRecord] = []
    for r in records:
        out.append(ProvenanceRecord(
            timestamp_unix=int(r.get("timestamp_unix", 0)),
            software=str(r.get("software", "")),
            parameters=_safe_json_dict(r.get("parameters_json", "{}")),
            input_refs=_safe_json_list(r.get("input_refs_json", "[]")),
            output_refs=_safe_json_list(r.get("output_refs_json", "[]")),
        ))
    return out


def _decode_provenance_json(blob: str) -> list[ProvenanceRecord]:
    try:
        data = json.loads(blob) if blob else []
    except json.JSONDecodeError:
        return []
    if not isinstance(data, list):
        return []
    out: list[ProvenanceRecord] = []
    for r in data:
        if not isinstance(r, dict):
            continue
        out.append(ProvenanceRecord(
            timestamp_unix=int(r.get("timestampUnix") or r.get("timestamp_unix") or 0),
            software=str(r.get("software", "")),
            parameters=r.get("parameters", {}) if isinstance(r.get("parameters"), dict) else {},
            input_refs=[str(x) for x in (r.get("inputRefs") or r.get("input_refs") or [])],
            output_refs=[str(x) for x in (r.get("outputRefs") or r.get("output_refs") or [])],
        ))
    return out


def _safe_json_list(value: str | list) -> list[str]:
    if isinstance(value, list):
        return [str(x) for x in value]
    try:
        parsed = json.loads(value) if value else []
    except json.JSONDecodeError:
        return []
    return [str(x) for x in parsed] if isinstance(parsed, list) else []


def _safe_json_dict(value: str | dict) -> dict[str, object]:
    if isinstance(value, dict):
        return value
    try:
        parsed = json.loads(value) if value else {}
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


# ----------------------------------------------------- M24 chromatograms ---


def _read_chromatograms(run_group) -> list[Chromatogram]:
    """Read ``<run>/chromatograms/`` into a list of :class:`Chromatogram`.

    Accepts either an ``h5py.Group`` or a :class:`StorageGroup`
    (). Returns an empty list when the group is absent
    (v0.3 backward compat).
    """
    sgroup = _wrap_hdf5_group(run_group)
    if not sgroup.has_child("chromatograms"):
        return []
    g = sgroup.open_group("chromatograms")
    count = int(io.read_int_attr(g, "count", default=0) or 0)
    if count <= 0:
        return []
    time_all = np.asarray(g.open_dataset("time_values").read())
    int_all  = np.asarray(g.open_dataset("intensity_values").read())
    idx = g.open_group("chromatogram_index")
    lengths       = np.asarray(idx.open_dataset("lengths").read())
    # offsets always synthesized from cumsum(lengths).
    from .genomic_index import _offsets_from_lengths
    offsets = _offsets_from_lengths(lengths.astype(np.uint64, copy=False))
    types         = np.asarray(idx.open_dataset("types").read())
    target_mzs    = np.asarray(idx.open_dataset("target_mzs").read())
    precursor_mzs = np.asarray(idx.open_dataset("precursor_mzs").read())
    product_mzs   = np.asarray(idx.open_dataset("product_mzs").read())

    from .signal_array import SignalArray

    out: list[Chromatogram] = []
    for i in range(count):
        off = int(offsets[i])
        n   = int(lengths[i])
        out.append(Chromatogram(
            signal_arrays={
                "time": SignalArray(data=np.asarray(time_all[off:off+n], dtype="<f8").copy()),
                "intensity": SignalArray(data=np.asarray(int_all[off:off+n], dtype="<f8").copy()),
            },
            axes=[],
            chromatogram_type=ChromatogramType(int(types[i])),
            target_mz=float(target_mzs[i]),
            precursor_mz=float(precursor_mzs[i]),
            product_mz=float(product_mzs[i]),
        ))
    return out


def write_chromatograms_to_run_group(
    run_group, chromatograms: list[Chromatogram],
) -> None:
    """Write ``chromatograms`` under ``<run>/chromatograms/``.

    Does nothing when the list is empty. The ``run_group`` argument
    accepts either an ``h5py.Group`` or a
    :class:`~ttio.providers.base.StorageGroup`; non-HDF5 providers
    route through :meth:`StorageGroup.create_dataset` so the same
    logical data lands in Memory / SQLite / Zarr backends.

    ``chromatogram_index/offsets`` is never written — readers derive
    it from ``cumsum(lengths)``.
    """
    if not chromatograms:
        return
    native = getattr(run_group, "_grp", None)
    if isinstance(run_group, StorageGroup) and native is None:
        # Non-HDF5 provider — write through the protocol.
        from . import _hdf5_io as _io
        from .enums import Precision
        g = run_group.create_group("chromatograms")
        _io.write_int_attr(g, "count", len(chromatograms))
        total = sum(len(c.time_array) for c in chromatograms)
        time_all = np.empty(total, dtype="<f8")
        int_all  = np.empty(total, dtype="<f8")
        offsets  = np.empty(len(chromatograms), dtype="<i8")
        lengths  = np.empty(len(chromatograms), dtype="<u4")
        types    = np.empty(len(chromatograms), dtype="<i4")
        targets  = np.empty(len(chromatograms), dtype="<f8")
        precs    = np.empty(len(chromatograms), dtype="<f8")
        prods    = np.empty(len(chromatograms), dtype="<f8")
        cursor = 0
        for i, c in enumerate(chromatograms):
            n = len(c.time_array)
            time_all[cursor:cursor+n] = c.time_array.data.astype("<f8", copy=False)
            int_all [cursor:cursor+n] = c.intensity_array.data.astype("<f8", copy=False)
            offsets[i] = cursor
            lengths[i] = n
            types[i]   = int(c.chromatogram_type)
            targets[i] = c.target_mz
            precs[i]   = c.precursor_mz
            prods[i]   = c.product_mz
            cursor += n
        g.create_dataset("time_values", Precision.FLOAT64, total).write(time_all)
        g.create_dataset("intensity_values", Precision.FLOAT64, total).write(int_all)
        idx = g.create_group("chromatogram_index")
        n_chrom = len(chromatograms)
        idx.create_dataset("lengths",       Precision.UINT32,  n_chrom).write(lengths)
        idx.create_dataset("types",         Precision.INT32,   n_chrom).write(types)
        idx.create_dataset("target_mzs",    Precision.FLOAT64, n_chrom).write(targets)
        idx.create_dataset("precursor_mzs", Precision.FLOAT64, n_chrom).write(precs)
        idx.create_dataset("product_mzs",   Precision.FLOAT64, n_chrom).write(prods)
        return

    # HDF5 fast path (preserves byte parity with pre-M64.5 files).
    target = native if native is not None else run_group
    g = target.create_group("chromatograms")
    g.attrs["count"] = np.int64(len(chromatograms))

    total = sum(len(c.time_array) for c in chromatograms)
    time_all = np.empty(total, dtype="<f8")
    int_all  = np.empty(total, dtype="<f8")
    offsets  = np.empty(len(chromatograms), dtype="<i8")
    lengths  = np.empty(len(chromatograms), dtype="<u4")
    types    = np.empty(len(chromatograms), dtype="<i4")
    targets  = np.empty(len(chromatograms), dtype="<f8")
    precs    = np.empty(len(chromatograms), dtype="<f8")
    prods    = np.empty(len(chromatograms), dtype="<f8")

    cursor = 0
    for i, c in enumerate(chromatograms):
        n = len(c.time_array)
        time_all[cursor:cursor+n] = c.time_array.data.astype("<f8", copy=False)
        int_all [cursor:cursor+n] = c.intensity_array.data.astype("<f8", copy=False)
        offsets[i] = cursor
        lengths[i] = n
        types[i]   = int(c.chromatogram_type)
        targets[i] = c.target_mz
        precs[i]   = c.precursor_mz
        prods[i]   = c.product_mz
        cursor += n

    g.create_dataset("time_values",      data=time_all)
    g.create_dataset("intensity_values", data=int_all)
    idx = g.create_group("chromatogram_index")
    idx.create_dataset("lengths",       data=lengths)
    idx.create_dataset("types",         data=types)
    idx.create_dataset("target_mzs",    data=targets)
    idx.create_dataset("precursor_mzs", data=precs)
    idx.create_dataset("product_mzs",   data=prods)

class _SpectralUnitView:
    """A view over one scheduling unit.

    Holds the unit's channel columns and builds a spectrum only when one
    is indexed, so a caller that reads part of a unit pays for that part.
    """

    __slots__ = ("_run", "_first", "_n", "_value_start", "_columns")

    def __init__(self, run: "AcquisitionRun", unit: tuple):
        _, first, n, value_start, value_end = unit
        self._run, self._first, self._n = run, first, n
        self._value_start = value_start
        count = value_end - value_start
        self._columns = {}
        for c in run.channel_names:
            try:
                self._columns[c] = run.channel_range(c, value_start, count)
            except KeyError:
                continue

    def __len__(self) -> int:
        return self._n

    @property
    def value_start(self) -> int:
        """Run-global value index :meth:`column` starts at."""
        return self._value_start

    def column(self, channel: str):
        """The unit's decoded values for one channel, without building a
        single spectrum; ``None`` when the view has no such channel.

        Building spectra dominates a bulk read, and in CPython it
        dominates it completely. Measured on an Orbitrap Exploris run,
        342546 spectra over 16 units, through ``for_each_block``, MB/s:

        =========  =========  =========  ======
        threads    spectra    columns    gain
        =========  =========  =========  ======
        1          98.9       1045.5     10.6x
        4          102.4      1523.2     14.9x
        16         100.6      1614.9     16.1x
        =========  =========  =========  ======

        Two things follow. Python's columns very nearly match
        Objective-C's 1637.8, so the SDK's usual tenfold deficit is
        per-spectrum object construction, not the data path. And unlike
        the spectrum path, which is flat, the column path **scales with
        workers**, because the numpy reduction releases the GIL. A
        caller reading columns can raise ``threads`` past the capped
        default and gain from it.

        Spectrum ``k`` of the unit covers
        ``run.index.offsets[first_spectrum + k] - view.value_start``
        onward, for ``run.index.lengths[first_spectrum + k]`` values.

        ⚠️ The array is the view's own buffer, not a copy. Do not
        mutate it, and do not keep it past the visitor call.

        Java: ``SpectralUnitView.column``; Objective-C:
        ``-unitColumnForChannel:valueStart:``.
        """
        return self._columns.get(channel)

    def __getitem__(self, k: int) -> Spectrum:
        if k < 0:
            k += self._n
        if k < 0 or k >= self._n:
            raise IndexError(f"index {k} beyond unit spectrum count {self._n}")
        i = self._first + k
        o = int(self._run.index.offsets[i]) - self._value_start
        ln = int(self._run.index.lengths[i])
        return self._run._build_spectrum(
            i, {c: a[o:o + ln] for c, a in self._columns.items()})

    def __iter__(self):
        for k in range(self._n):
            yield self[k]
