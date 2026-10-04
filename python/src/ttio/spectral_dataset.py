"""``SpectralDataset`` — root ``.tio`` reader / writer façade.

This class is the main public entry point. It delegates storage to a
:class:`~ttio.providers.base.StorageProvider`, provides mapping-style access to runs, and exposes the
feature-flag / identifications / quantifications / provenance metadata.

Only reading and minimal writing are implemented in M16. Full writing
support (feature flags, compound datasets, signal channels for new runs)
uses the same helpers in :mod:`_hdf5_io` and is fleshed out alongside the
mzML/nmrML importers in M16.7.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType, TracebackType
from typing import Any, Iterator, Mapping

import h5py
import numpy as np

from . import _hdf5_io as io
from ._rwlock import RWLock
from .access_policy import AccessPolicy
from .io.progress import ProgressSinkLike, _fire
from .acquisition_run import AcquisitionRun
from .assembly import (  # M98
    AssemblyGraph,
    WrittenAssemblyGraph,
    write_assembly_graph,
)
from .genomic.reference_import import ReferenceImport  # tio-browser Phase 0
from .genomic_run import GenomicRun  # M82
from .enums import EncryptionLevel, SpectrumKind
from .feature_flags import FeatureFlags
from .ir_image import IRImage
from .ms_image import MSImage
from .raman_image import RamanImage
from .identification import Identification
from .providers import StorageProvider, open_provider
from .providers.base import StorageGroup
from .providers.hdf5 import Hdf5Provider
from .provenance import ProvenanceRecord
from .quantification import Quantification
from .sample import Sample  # Stage 6 (transport-spec v0.11, Deferral 2)
from .subject import Subject  # Stage 6 (transport-spec v0.11, Deferral 2)
from .written_genomic_run import WrittenGenomicRun  # M82

# P3.10: genomic-write + metadata-IO helper subsystems extracted into
# private submodules (pure code movement; no API/wire/behaviour change).
from . import _dataset_write_genomic as _gw


def _write_genomic_run_default(study, g_group, name: str, run: WrittenGenomicRun) -> None:
    """Write one genomic run: blocks_v1 through the stream writer unless
    the run asks for the v1.8 whole-channel layout."""
    if getattr(run, "opt_legacy_whole_channel", False):
        _gw._write_genomic_run(g_group, name, run)
        return
    from .genomic.stream_writer import GenomicStreamWriter
    from .providers.hdf5 import _Group as _H5Group
    study_sg = _H5Group(study) if isinstance(study, h5py.Group) else study
    with GenomicStreamWriter(
            study_sg, name,
            acquisition_mode=run.acquisition_mode, reference_uri=run.reference_uri,
            platform=run.platform, sample_name=run.sample_name,
            reference_chrom_seqs=run.reference_chrom_seqs,
            embed_reference=run.embed_reference,
            opt_disable_qualities_v5=run.opt_disable_qualities_v5,
            signal_codec_overrides=run.signal_codec_overrides,
            signal_compression=run.signal_compression,
            provenance_records=run.provenance_records,
            read_role=run.read_role,
            ref_diff_slice_bytes=run.ref_diff_slice_bytes,
            group_reads=run.opt_group_reads) as w:
        w.append_batch(run)
from . import _dataset_write_metadata as _mw
# Back-compat re-export: tests/test_references_accessor.py imports this
# private from the old path.
from ._dataset_write_genomic import _embed_references_for_runs  # noqa: F401

# M23 sentinel: returned by ``read_lock``/``write_lock`` when ``thread_safe``
# is False so call sites can use ``with ds.read_lock(): ...`` unconditionally.
class _NullGuard:
    def __enter__(self) -> "_NullGuard":
        return self

    def __exit__(self, *exc: object) -> None:
        return None


_NULL_GUARD = _NullGuard()


def _split_run_names(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    return tuple(name for name in value.split(",") if name)


@dataclass(slots=True)
class SpectralDataset:
    """Root container for an TTI-O ``.tio`` file.

    Owns a top-level ``study/`` group plus zero or more named MS
    acquisition runs, zero or more named NMR-spectrum collections,
    the dataset-wide identifications, quantifications, provenance
    records, and an optional transition list.

    Persistence uses a
    :class:`~ttio.providers.base.StorageProvider` (HDF5 by
    default); callers may supply another via the ``provider`` kwarg.

    Notes
    -----
    API status: Stable. ``Encryptable`` conformance is delivered in
    slice 41.5 when the encryption manager lands in Python.

    Cross-language equivalents
    --------------------------
    Objective-C: ``TTIOSpectralDataset`` · Java:
    ``global.thalion.ttio.SpectralDataset``.
    """

    path: Path
    feature_flags: FeatureFlags
    title: str
    isa_investigation_id: str
    ms_runs: dict[str, AcquisitionRun] = field(default_factory=dict)
    nmr_runs: dict[str, AcquisitionRun] = field(default_factory=dict)
    genomic_runs: dict[str, GenomicRun] = field(default_factory=dict)  # M82
    # Assembly graphs at /study/assembly_graphs/<name>/, populated on
    # open() (M98, format-spec §11a). Empty for files without the
    # opt_assembly_graph feature.
    assembly_graphs: dict[str, AssemblyGraph] = field(default_factory=dict)
    # Embedded references at /study/references/<uri>/, populated on
    # open() from the canonical embed layout (see
    # :func:`_embed_references_for_runs` for the writer side). Empty
    # for files without an embedded reference, including those whose
    # genomic runs use ``embed_reference=False``. Mirrors Java's
    # ``SpectralDataset.references()`` (added Phase 0 of tio-browser).
    _references: dict[str, "ReferenceImport"] = field(default_factory=dict)
    encrypted_algorithm: str = ""
    _closed: bool = False
    _remote_fileobj: Any = None  # fsspec file-like kept alive when remote
    _lock: RWLock | None = None  # set when opened with thread_safe=True
    provider: StorageProvider | None = None  # owning storage provider
    # ``provider`` is the backend abstraction introduced in M39. All code
    # reaches storage through ``provider.root_group()`` and the StorageGroup
    # protocol; the raw-handle escape hatch ``native_handle()`` is deprecated
    # (P3.9) and no longer used by mainline code.
    # Encryptable conformance.
    _access_policy: AccessPolicy | None = field(default=None, repr=False)
    # Lazy image cache — populated on first access of the `image` property.
    _image_cache_loaded: bool = field(default=False, repr=False)
    _image_cache: "MSImage | None" = field(default=None, repr=False)
    # Lazy raman_image cache — populated on first access of the `raman_image` property.
    _raman_image_cache_loaded: bool = field(default=False, repr=False)
    _raman_image_cache: "RamanImage | None" = field(default=None, repr=False)
    # Lazy ir_image cache — populated on first access of the `ir_image` property.
    _ir_image_cache_loaded: bool = field(default=False, repr=False)
    _ir_image_cache: "IRImage | None" = field(default=None, repr=False)
    # Stage 6 (transport-spec v0.11, Deferral 2): lazy subject + sample
    # caches. Populated on first access of the `subjects` / `samples`
    # property by enumerating `/study/subjects/` / `/study/samples/`.
    _subjects_cache_loaded: bool = field(default=False, repr=False)
    _subjects_cache: list[Subject] = field(default_factory=list, repr=False)
    _samples_cache_loaded: bool = field(default=False, repr=False)
    _samples_cache: list[Sample] = field(default_factory=list, repr=False)

    # ------------------------------------------------------------- lifecycle

    @classmethod
    def open(
        cls,
        path: str | Path,
        *,
        provider: str | None = None,
        thread_safe: bool = False,
        writable: bool = False,
        **fsspec_kwargs: Any,
    ) -> "SpectralDataset":
        """Open a ``.tio`` dataset from a local path or cloud URL.

        URLs with a scheme recognised by :data:`ttio.remote.REMOTE_SCHEMES`
        (``s3://``, ``http(s)://``, ``gs://``, ``az://``, ``file://``) are
        routed through fsspec and read lazily — only the HDF5 metadata and
        any actively touched chunks are fetched. Extra keyword arguments are
        forwarded to :func:`fsspec.open` and are typically used for
        cloud-backend options (``anon=True``, ``key=...``, ...).

        Parameters
        ----------
        writable:
            If ``True``, open the file in read-write mode (``"r+"``) so
            that in-place operations such as
            :meth:`AcquisitionRun.encrypt_with_key` can write back to the
            same file handle. Ignored for remote URLs (which are always
            read-only). Default: ``False``.

        M39: a :class:`~ttio.providers.Hdf5Provider` is constructed for
        the target and exposed as ``dataset.provider``. The dataset is
        provider-backed; reach storage through ``dataset.provider``.
        """
        from .remote import is_remote_url, open_remote_file

        if is_remote_url(path):
            fileobj = open_remote_file(str(path), **fsspec_kwargs)
            try:
                f = h5py.File(fileobj, "r")
            except Exception:
                fileobj.close()
                raise
            try:
                provider = Hdf5Provider._from_open_h5py(f)
                return cls._from_provider(Path(str(path)), provider,
                                          thread_safe=thread_safe,
                                          _remote_fileobj=fileobj)
            except Exception:
                f.close()
                fileobj.close()
                raise

        # URL scheme detection routes non-HDF5 providers
        # (``memory://...``, ``sqlite://...``, ``dir://...``) through
        # the storage protocol. Bare paths still open via HDF5 for
        # byte-parity with pre-M64.5 files.
        path_str = str(path)
        mode = "r+" if writable else "r"
        if "://" in path_str and not path_str.startswith("file://"):
            provider = open_provider(path_str, mode=mode)
            try:
                return cls._from_provider(Path(path_str), provider,
                                           thread_safe=thread_safe)
            except Exception:
                provider.close()
                raise

        p = Path(path)
        # if an explicit provider name is given for a bare path, route
        # through open_provider so Memory / SQLite / Zarr backends work.
        if provider is not None and provider not in ("hdf5", "h5", "h5py"):
            sp = open_provider(str(path), provider=provider, mode=mode)
            try:
                return cls._from_provider(p, sp, thread_safe=thread_safe)
            except Exception:
                sp.close()
                raise
        provider = Hdf5Provider.open(str(p), mode=mode)
        try:
            return cls._from_provider(p, provider, thread_safe=thread_safe)
        except Exception:
            provider.close()
            raise

    @classmethod
    def _from_provider(
        cls,
        path: Path,
        provider: StorageProvider,
        *,
        thread_safe: bool = False,
        _remote_fileobj: Any = None,
    ) -> "SpectralDataset":
        """Open-side constructor for non-HDF5 providers.

        Reads everything through the :class:`StorageGroup` protocol so
        Memory / SQLite / Zarr backends work without touching h5py.
        """
        root = provider.root_group()
        version, features = io.read_feature_flags(root)
        flags = FeatureFlags.from_iterable(version, features)
        encrypted = io.read_string_attr(root, "encrypted", default="") or ""
        if not root.has_child("study"):
            raise ValueError(f"{path}: missing /study group; not an .tio file")
        study = root.open_group("study")
        title = io.read_string_attr(study, "title", default="") or ""
        isa = io.read_string_attr(study, "isa_investigation_id", default="") or ""

        # Remote (fsspec) datasets keep the per-spectrum hyperslab read so
        # a single random access only fetches its own slice over the
        # network; local datasets bulk-read each channel column once. The
        # materialized spectra are byte-identical either way.
        bulk_read = _remote_fileobj is None

        ms_runs: dict[str, AcquisitionRun] = {}
        if study.has_child("ms_runs"):
            ms_group = study.open_group("ms_runs")
            names = _split_run_names(
                io.read_string_attr(ms_group, "_run_names", default="") or ""
            )
            for name in names:
                if ms_group.has_child(name):
                    run = AcquisitionRun.open(
                        ms_group.open_group(name), name, bulk_read=bulk_read)
                    run._set_persistence_context(str(path), name)
                    ms_runs[name] = run

        nmr_runs: dict[str, AcquisitionRun] = {}
        if study.has_child("nmr_runs"):
            nmr_group = study.open_group("nmr_runs")
            names = _split_run_names(
                io.read_string_attr(nmr_group, "_run_names", default="") or ""
            )
            for name in names:
                if nmr_group.has_child(name):
                    run = AcquisitionRun.open(
                        nmr_group.open_group(name), name, bulk_read=bulk_read)
                    run._set_persistence_context(str(path), name)
                    nmr_runs[name] = run

        genomic_runs_map: dict[str, GenomicRun] = {}  # M82
        if study.has_child("genomic_runs"):
            g_group = study.open_group("genomic_runs")
            # /study/references threaded into each run for REF_DIFF decode
            # via the StorageGroup protocol (P3.9). None when absent.
            refs_group = (
                study.open_group("references")
                if study.has_child("references")
                else None
            )
            names = _split_run_names(
                io.read_string_attr(g_group, "_run_names", default="") or ""
            )
            for name in names:
                if g_group.has_child(name):
                    genomic_runs_map[name] = GenomicRun.open(
                        g_group.open_group(name), name,
                        references_group=refs_group,
                        bulk_read=bulk_read,
                    )

        # M98: provider-agnostic assembly_graphs read. Unreadable
        # graphs are skipped the way unreadable genomic runs are.
        assembly_graphs_map: dict[str, AssemblyGraph] = {}
        if study.has_child("assembly_graphs"):
            ag_group = study.open_group("assembly_graphs")
            ag_names = io.read_string_attr(
                ag_group, "_graph_names", default="") or ""
            for gname in ag_names.split(","):
                gname = gname.strip()
                if not gname or not ag_group.has_child(gname):
                    continue
                try:
                    assembly_graphs_map[gname] = AssemblyGraph.open(
                        ag_group.open_group(gname), gname)
                except ValueError:
                    continue

        references_map = _gw._load_references_provider(study)

        return cls(
            path=path,
            feature_flags=flags,
            title=title,
            isa_investigation_id=isa,
            ms_runs=ms_runs,
            nmr_runs=nmr_runs,
            genomic_runs=genomic_runs_map,  # M82
            assembly_graphs=assembly_graphs_map,  # M98
            _references=references_map,
            encrypted_algorithm=encrypted,
            _remote_fileobj=_remote_fileobj,
            _lock=(RWLock() if thread_safe else None),
            provider=provider,
        )

    # ----------------------------------------------------- thread safety (M23)

    @property
    def is_thread_safe(self) -> bool:
        """True iff this dataset was opened with ``thread_safe=True``."""
        return self._lock is not None

    def read_lock(self) -> Any:
        """Context manager acquiring the shared read lock.

        A no-op when ``thread_safe`` was not set at open time, so call sites
        can use ``with ds.read_lock(): ...`` unconditionally.
        """
        return self._lock.read() if self._lock is not None else _NULL_GUARD

    def write_lock(self) -> Any:
        """Context manager acquiring the exclusive write lock (no-op when
        ``thread_safe`` was not set at open time)."""
        return self._lock.write() if self._lock is not None else _NULL_GUARD

    def close(self) -> None:
        """Release the underlying HDF5 file, provider, and remote handles.

        Idempotent: calling :meth:`close` on an already-closed dataset is
        a no-op. Invoked automatically by the context-manager protocol
        (``with SpectralDataset.open(...) as ds:``) and from ``__del__``
        as a best-effort safety net.
        """
        with self.write_lock():
            if not self._closed:
                # Close via the provider when we have one — it owns the
                # h5py.File and any fsspec file-like. Fall back to direct
                # close for legacy instances constructed without M39.
                if self.provider is not None:
                    self.provider.close()
                if self._remote_fileobj is not None:
                    try:
                        self._remote_fileobj.close()
                    except Exception:
                        pass
                    self._remote_fileobj = None
                self._closed = True

    def __enter__(self) -> "SpectralDataset":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    # --------------------------------------------------------------- queries

    @property
    def is_encrypted(self) -> bool:
        """Whether the dataset's signal payload is encrypted.

        Returns
        -------
        bool
            True when an ``@encrypted_algorithm`` attribute was recorded
            on the root group at write time (regardless of which
            channels carry ciphertext). False for plaintext datasets.
        """
        return bool(self.encrypted_algorithm)

    @property
    def references(self) -> Mapping[str, ReferenceImport]:
        """Embedded references at ``/study/references/<uri>/``, keyed
        by URI string.

        Populated on :meth:`open` from the canonical embed layout (see
        :func:`_embed_references_for_runs`). Empty for files without
        embedded references — including those whose genomic runs use
        ``embed_reference=False`` (the v1.2+ default), which keep the
        ``reference_uri`` attribute but rely on the
        :class:`ReferenceResolver` for byte resolution.

        The return is a read-only view: mutations on the returned
        mapping do not affect the dataset. Mirrors Java's
        ``SpectralDataset.references()`` (added Phase 0 of
        tio-browser).

        :since: 1.1.0
        """
        # Deliberate tightening over `all_runs` etc.: /study/references/
        # is a write-once-on-open snapshot, so we expose it as a true
        # read-only view (mutation raises TypeError) rather than a
        # cheap dict copy. Added in 1.1.0.
        return MappingProxyType(self._references)

    def _lazy_ms_image(self) -> "MSImage | None":
        """Lazy-read the embedded MSImage at /study/image_cube.

        Reads the cube on first access, caches the result. Returns
        ``None`` when no image group exists.
        """
        if not self._image_cache_loaded:
            self._image_cache_loaded = True
            if self.provider is None:
                self._image_cache = None
            else:
                root = self.provider.root_group()
                if root.has_child("study"):
                    study = root.open_group("study")
                    self._image_cache = MSImage.read_from(study)
                # else: _image_cache stays None (already initialized)
        return self._image_cache

    def _lazy_raman_image(self) -> "RamanImage | None":
        """Lazy-read the embedded RamanImage at /study/raman_image_cube.

        Reads the cube on first access, caches the result. Returns
        ``None`` when no raman image group exists.
        """
        if not self._raman_image_cache_loaded:
            self._raman_image_cache_loaded = True
            if self.provider is None:
                self._raman_image_cache = None
            else:
                root = self.provider.root_group()
                if root.has_child("study"):
                    study = root.open_group("study")
                    self._raman_image_cache = RamanImage.read_from(study)
                # else: _raman_image_cache stays None (already initialized)
        return self._raman_image_cache

    def _lazy_ir_image(self) -> "IRImage | None":
        """Lazy-read the embedded IRImage at /study/ir_image_cube.

        Reads the cube on first access, caches the result. Returns
        ``None`` when no IR image group exists. Stage 5.2
        (transport-spec v0.11, Deferral 1).
        """
        if not self._ir_image_cache_loaded:
            self._ir_image_cache_loaded = True
            if self.provider is None:
                self._ir_image_cache = None
            else:
                root = self.provider.root_group()
                if root.has_child("study"):
                    study = root.open_group("study")
                    self._ir_image_cache = IRImage.read_from(study)
                # else: _ir_image_cache stays None (already initialized)
        return self._ir_image_cache

    def image_for_kind(self, kind: "ImageKind") -> "Image | None":
        """The embedded :class:`~ttio.image.Image` for the given
        :class:`~ttio.enums.ImageKind`, or ``None`` when that modality
        is absent.

        Lazy per-kind: reads the cube on first access and caches the
        result. Replaces the former typed ``image`` / ``raman_image``
        / ``ir_image`` accessors with a uniform lookup.

        :since: 1.2.0
        """
        from .enums import ImageKind
        if kind == ImageKind.MS:
            return self._lazy_ms_image()
        if kind == ImageKind.RAMAN:
            return self._lazy_raman_image()
        if kind == ImageKind.IR:
            return self._lazy_ir_image()
        raise ValueError(f"unknown ImageKind: {kind!r}")

    @property
    def images(self) -> "dict[ImageKind, Image]":
        """Embedded images keyed by :class:`~ttio.enums.ImageKind`,
        containing only the modalities present on this dataset.

        Each entry is materialised lazily via :meth:`image_for_kind`
        (per-kind cache preserved). Absent modalities are omitted.

        :since: 1.2.0
        """
        from .enums import ImageKind
        out: "dict[ImageKind, Image]" = {}
        for k in (ImageKind.MS, ImageKind.RAMAN, ImageKind.IR):
            img = self.image_for_kind(k)
            if img is not None:
                out[k] = img
        return out

    @property
    def subjects(self) -> list[Subject]:
        """List of :class:`~ttio.subject.Subject` rows persisted under
        ``/study/subjects/`` on this dataset, in on-disk iteration
        order. Empty list when no Subjects were written (which is
        most pre-Stage-6 files). Lazily populated on first access.

        Stage 6 (transport-spec v0.11, Deferral 2). Mirrors Java
        :meth:`global.thalion.ttio.SpectralDataset.subjects`.

        :since: 1.4.0
        """
        if not self._subjects_cache_loaded:
            self._subjects_cache_loaded = True
            self._subjects_cache = _mw._read_subjects(self.provider)
        return self._subjects_cache

    @property
    def samples(self) -> list[Sample]:
        """List of :class:`~ttio.sample.Sample` rows persisted under
        ``/study/samples/`` on this dataset, in on-disk iteration
        order. Empty list when no Samples were written. Lazily
        populated on first access.

        Stage 6 (transport-spec v0.11, Deferral 2). Mirrors Java
        :meth:`global.thalion.ttio.SpectralDataset.samples`.

        :since: 1.4.0
        """
        if not self._samples_cache_loaded:
            self._samples_cache_loaded = True
            self._samples_cache = _mw._read_samples(self.provider)
        return self._samples_cache

    @property
    def all_runs(self) -> Mapping[str, AcquisitionRun]:
        """Union of MS and NMR runs, keyed by run name.

        Note: This DOES NOT include ``genomic_runs`` — those have a
        different return type (:class:`GenomicRun` rather than
        :class:`AcquisitionRun`). For modality-agnostic iteration
        across all run types, use :meth:`all_runs_unified` (Phase 1).
        """
        merged: dict[str, AcquisitionRun] = dict(self.ms_runs)
        for k, v in self.nmr_runs.items():
            merged.setdefault(k, v)
        return merged

    # ── Phase 2 (post-M91) — canonical unified runs accessor ────────

    @property
    def runs(self) -> Mapping[str, "Run"]:
        """Canonical mapping over every run in the file (MS + NMR +
        genomic), keyed by run name.

        Values conform to the :class:`ttio.protocols.run.Run`
        Protocol so callers can iterate uniformly without knowing
        the underlying modality:

            for name, run in ds.runs.items():
                print(f"{name}: {len(run)} measurements")

        Use :meth:`runs_of_modality` to narrow by class, or
        :meth:`runs_for_sample` to filter by provenance sample URI.

        Phase 2 promotes this to the canonical access pattern.
        Backward-compat: the legacy ``ms_runs`` / ``nmr_runs`` /
        ``genomic_runs`` dicts and the MS+NMR-only ``all_runs``
        property continue to work; new code should prefer ``runs``.
        """
        merged: dict[str, Any] = dict(self.ms_runs)
        for k, v in self.nmr_runs.items():
            merged.setdefault(k, v)
        for k, v in self.genomic_runs.items():
            merged.setdefault(k, v)
        return merged

    @property
    def all_runs_unified(self) -> Mapping[str, "Run"]:
        """Deprecated alias for :attr:`runs`. Kept for the brief
        Phase 1 → Phase 2 transition; remove in v1.0."""
        return self.runs

    def runs_for_sample(self, sample_uri: str) -> Mapping[str, "Run"]:
        """Return every run associated with ``sample_uri``.

        A run is considered associated when its
        :meth:`ttio.protocols.run.Run.provenance_chain` carries
        ``sample_uri`` in any record's ``input_refs``. Walks all
        modalities (MS, NMR, genomic) uniformly via the Run
        Protocol — closes the M91 cross-modality query gap that
        previously had to fork on access pattern.

        Returns a dict keyed by run name; empty when no run
        matches.
        """
        out: dict[str, Any] = {}
        for name, run in self.runs.items():
            try:
                chain = run.provenance_chain()
            except Exception:
                continue
            for prov in chain:
                if sample_uri in prov.input_refs:
                    out[name] = run
                    break
        return out

    def runs_of_modality(self, run_type: type) -> Mapping[str, "Run"]:
        """Return every run whose value is an instance of ``run_type``.

        Pass :class:`AcquisitionRun` to get the union of MS + NMR
        runs (any spectrum-class subtype); pass :class:`GenomicRun`
        to get genomic only. The return is a thin filter over
        :attr:`all_runs_unified`.
        """
        return {
            name: run
            for name, run in self.runs.items()
            if isinstance(run, run_type)
        }

    @property
    def study_group(self) -> StorageGroup:
        """The ``/study`` :class:`StorageGroup` (open with
        ``writable=True`` to add runs through a stream writer)."""
        return self.provider.root_group().open_group("study")

    def _study_target(self) -> Any:
        """Return the IO target representing ``/study``.

        Returns the :class:`StorageGroup` wrapping ``/study``. The
        helpers in :mod:`_hdf5_io` accept either form.
        """
        return self.provider.root_group().open_group("study")

    def _study_has_child(self, name: str) -> bool:
        return self.provider.root_group().open_group("study").has_child(name)

    def identifications(self) -> list[Identification]:
        """Return the dataset's identification records.

        Reads from the compound ``/study/identifications`` dataset when
        present, falling back to the legacy ``@identifications_json``
        attribute on ``/study`` for older files.

        Returns
        -------
        list[Identification]
            Identification records in stored order. Empty when no
            identifications group or JSON attribute is present.
        """
        with self.read_lock():
            study = self._study_target()
            if self._study_has_child("identifications"):
                return [
                    Identification(
                        run_name=r["run_name"],
                        spectrum_index=int(r["spectrum_index"]),
                        chemical_entity=r["chemical_entity"],
                        confidence_score=float(r["confidence_score"]),
                        evidence_chain=_mw._maybe_json_list(r.get("evidence_chain_json", "[]")),
                    )
                    for r in io.read_compound_dataset(study, "identifications")
                ]
            blob = io.read_string_attr(study, "identifications_json", default="")
            return _mw._decode_identifications_json(blob) if blob else []

    def quantifications(self) -> list[Quantification]:
        """Return the dataset's quantification records.

        Reads from the compound ``/study/quantifications`` dataset when
        present (with optional per-row units from the
        ``@quantification_units`` JSON sidecar attribute), falling back
        to the legacy ``@quantifications_json`` attribute for older
        files.

        Returns
        -------
        list[Quantification]
            Quantification records in stored order. Empty when no
            quantifications group or JSON attribute is present.
        """
        with self.read_lock():
            study = self._study_target()
            if self._study_has_child("quantifications"):
                # Optional sidecar `@quantification_units` JSON-array
                # attribute carries one unit per row. Absent on legacy
                # files; readers default missing rows to empty string.
                units_blob = io.read_string_attr(
                    study, "quantification_units", default="")
                units = json.loads(units_blob) if units_blob else []
                rows = list(io.read_compound_dataset(study, "quantifications"))
                out: list[Quantification] = []
                for i, r in enumerate(rows):
                    unit = units[i] if i < len(units) else ""
                    out.append(Quantification(
                        chemical_entity=r["chemical_entity"],
                        sample_ref=r["sample_ref"],
                        abundance=float(r["abundance"]),
                        normalization_method=r.get("normalization_method", ""),
                        unit=unit if isinstance(unit, str) else "",
                    ))
                return out
            blob = io.read_string_attr(study, "quantifications_json", default="")
            return _mw._decode_quantifications_json(blob) if blob else []

    def provenance(self) -> list[ProvenanceRecord]:
        """Return the dataset-level provenance chain.

        Reads from the compound ``/study/provenance`` dataset when
        present, falling back to the legacy ``@provenance_json``
        attribute on ``/study`` for older files. Per-run provenance is
        exposed separately by :class:`AcquisitionRun.provenance_records`.

        Returns
        -------
        list[ProvenanceRecord]
            Provenance records in stored order. Empty when neither the
            compound dataset nor the JSON attribute is present.
        """
        with self.read_lock():
            study = self._study_target()
            if self._study_has_child("provenance"):
                out: list[ProvenanceRecord] = []
                for r in io.read_compound_dataset(study, "provenance"):
                    out.append(
                        ProvenanceRecord(
                            timestamp_unix=int(r["timestamp_unix"]),
                            software=r["software"],
                            parameters=_mw._maybe_json_dict(r.get("parameters_json", "{}")),
                            input_refs=_mw._maybe_json_list(r.get("input_refs_json", "[]")),
                            output_refs=_mw._maybe_json_list(r.get("output_refs_json", "[]")),
                        )
                    )
                return out
            blob = io.read_string_attr(study, "provenance_json", default="")
            return _mw._decode_provenance_json(blob) if blob else []

    # ---- Encryptable conformance ----

    def encrypt_with_key(
        self, key: bytes, level: EncryptionLevel | None = None
    ) -> None:
        """Encrypt protectable content at the given granularity.

        For DATASET-level encryption (ObjC ``TTIOEncryptionLevelDataset``),
        encrypts every MS run's intensity channel in place. For finer-grained
        levels, callers should use ``run.encrypt_with_key`` directly.

        Matches ObjC ``-[TTIOSpectralDataset encryptWithKey:level:error:]``:
        after run-level encryption, marks the root with
        ``encrypted="aes-256-gcm"`` so :attr:`is_encrypted` and
        :attr:`encrypted_algorithm` round-trip across close/reopen.
        """
        for run in self.ms_runs.values():
            run.encrypt_with_key(key, level)
        self._mark_root_encrypted()

    def _mark_root_encrypted(self) -> None:
        """Persist ``encrypted=<algorithm>`` on the root, mirroring ObjC's
        ``-[TTIOSpectralDataset markRootEncryptedWithError:]``.

        Updates the in-memory ``encrypted_algorithm`` field as well so
        :attr:`is_encrypted` becomes True without requiring a reopen.
        """
        from .encryption import DEFAULT_ENCRYPTION_ALGORITHM

        self.provider.root_group().set_attribute("encrypted", DEFAULT_ENCRYPTION_ALGORITHM)
        object.__setattr__(self, "encrypted_algorithm", DEFAULT_ENCRYPTION_ALGORITHM)

    def decrypt_with_key(self, key: bytes) -> dict[str, bytes]:
        """Decrypt every MS run's intensity channel into an in-memory overlay.

        Returns a mapping of ``{run_name: plaintext_bytes}``.

        **Read-only / asymmetric**: the on-disk file is NOT modified,
        the root ``@encrypted`` attribute is left in place, and
        :attr:`is_encrypted` continues to return ``True`` on this
        instance and on any reopen. This is the asymmetric counterpart
        to :meth:`encrypt_with_key` (which IS persistent +
        flag-flipping) by design: in-memory rehydration lets a process
        read encrypted data without rewriting the file.

        To fully reverse encryption on disk and clear the
        ``@encrypted`` attribute, use the classmethod
        :meth:`decrypt_in_place` — close any open instance first.

        Side effect: each run's decrypted channel is cached in memory
        so ``run.object_at_index(i).intensity_array`` works without
        re-decrypting (Option 1 of the MCP-Server M5 handoff; mirrors
        ObjC ``-[TTIOSpectralDataset decryptWithKey:]`` rehydration
        semantics).

        Cross-language equivalents
        --------------------------
        Java: ``SpectralDataset.decryptWithKey(byte[])`` (same
        in-memory-only semantics) · Objective-C:
        ``-[TTIOSpectralDataset decryptWithKey:error:]`` (same).
        """
        return {name: run.decrypt_with_key(key) for name, run in self.ms_runs.items()}

    @classmethod
    def decrypt_in_place(cls, path: str | Path, key: bytes) -> None:
        """Strip AES-256-GCM encryption from a ``.tio`` file on disk.

        For every MS run with an encrypted intensity channel, replaces
        ``intensity_values_encrypted`` + IV + tag with a plaintext
        ``intensity_values`` dataset, then clears the root ``@encrypted``
        attribute. After this call the file is byte-compatible with the
        pre-encryption state and :attr:`is_encrypted` is ``False`` when
        reopened.

        Symmetric with :meth:`encrypt_with_key`: that method leaves the
        root ``@encrypted`` attribute set, this one removes it.

        The file must not be held open by another writer.

        Raises ``FileNotFoundError`` if the file does not exist,
        ``ValueError`` if ``key`` is not 32 bytes or any channel's tag
        does not verify.
        """
        from .encryption import (
            AES_KEY_LEN,
            decrypt_intensity_channel_in_run_in_place,
        )

        if len(key) != AES_KEY_LEN:
            raise ValueError(
                f"AES-256-GCM key must be {AES_KEY_LEN} bytes, got {len(key)}"
            )

        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

        with h5py.File(str(p), "r") as f:
            ms_group = f.get("study/ms_runs")
            if ms_group is None:
                run_names: list[str] = []
            else:
                names_attr = io.read_string_attr(
                    ms_group, "_run_names", default=""
                ) or ""
                run_names = [n for n in names_attr.split(",") if n]

        for run_name in run_names:
            decrypt_intensity_channel_in_run_in_place(str(p), run_name, key)

        with h5py.File(str(p), "r+") as f:
            if "encrypted" in f.attrs:
                del f.attrs["encrypted"]

    def access_policy(self) -> AccessPolicy | None:
        """Return the current access policy, or ``None`` if not set."""
        return self._access_policy

    def set_access_policy(self, policy: AccessPolicy | None) -> None:
        """Replace the current access policy."""
        object.__setattr__(self, "_access_policy", policy)

    # ---------------------------------------------------------------- writer

    @classmethod
    def write_minimal(
        cls,
        path: str | Path,
        *,
        title: str,
        isa_investigation_id: str,
        runs: Mapping[str, "WrittenRun"],
        genomic_runs: Mapping[str, WrittenGenomicRun] | None = None,  # M82
        assembly_graphs: Mapping[str, WrittenAssemblyGraph] | None = None,  # M98
        identifications: list[Identification] | None = None,
        quantifications: list[Quantification] | None = None,
        provenance: list[ProvenanceRecord] | None = None,
        features: list[str] | None = None,
        provider: str | StorageProvider = "hdf5",
        image: "MSImage | None" = None,
        raman_image: "RamanImage | None" = None,
        ir_image: "IRImage | None" = None,
        subjects: list[Subject] | None = None,
        samples: list[Sample] | None = None,
        progress: "ProgressSinkLike | None" = None,
    ) -> Path:
        """Write a minimal v1.1 ``.tio`` file from in-memory data.

        Parameters
        ----------
        provider
            : which storage backend to write through. The
            string ``"hdf5"`` (default) keeps byte-for-byte parity with
            pre-M64.5 files. Other values dispatch through
            :func:`open_provider` — ``"memory"``, ``"sqlite"``,
            ``"zarr"``, or any registered backend. A pre-opened
            :class:`StorageProvider` may also be passed; the caller
            owns its lifecycle in that case.
        """
        p = Path(path)
        # Stage 6 (transport-spec v0.11, Deferral 2): validate
        # Subject + Sample lists upfront so the writer fails fast on
        # duplicate IDs before any backend mutation. Soft-FK warnings
        # are emitted by the same call.
        subjects_list = list(subjects) if subjects else []
        samples_list = list(samples) if samples else []
        if subjects_list or samples_list:
            _mw._validate_subjects_and_samples(subjects_list, samples_list)
        feature_list = features or [
            "base_v1",
            "compound_identifications",
            "compound_quantifications",
            "compound_provenance",
            "compound_per_run_provenance",
            "opt_compound_headers",
        ]

        # Phase 2: ``runs`` may be a MIXED dict carrying both WrittenRun
        # (MS / NMR) and WrittenGenomicRun entries. Split them into the
        # legacy two-kwarg internal layout BEFORE any MS-only
        # introspection (e.g. activation_methods below). Callers using
        # the pre-Phase-2 form (separate ``runs=`` and ``genomic_runs=``
        # kwargs) are unaffected.
        if any(isinstance(v, WrittenGenomicRun) for v in runs.values()):
            split_ms: dict[str, WrittenRun] = {}
            split_g: dict[str, WrittenGenomicRun] = dict(genomic_runs or {})
            for name, value in runs.items():
                if isinstance(value, WrittenGenomicRun):
                    if name in split_g:
                        raise ValueError(
                            f"Phase 2 mixed runs dict: name {name!r} "
                            f"appears in both runs= and genomic_runs="
                        )
                    split_g[name] = value
                else:
                    split_ms[name] = value
            runs = split_ms
            genomic_runs = split_g

        # v1.0 single format-version stamp. Readers gate optional
        # features by the feature-flag list (opt_*), not by version
        # equality, so per-feature version bumps are unnecessary.
        any_m74 = any(
            run.activation_methods is not None for run in runs.values()
        )
        if features is None and any_m74:
            feature_list = feature_list + ["opt_ms2_activation_detail"]
        # M102: opt_pixel_coordinates when any run carries pixel columns.
        if (any(run.pixel_x is not None for run in runs.values())
                and "opt_pixel_coordinates" not in feature_list):
            feature_list = feature_list + ["opt_pixel_coordinates"]
        has_genomic = bool(genomic_runs)
        if has_genomic and "opt_genomic" not in feature_list:
            feature_list = feature_list + ["opt_genomic"]
        # M98: opt_assembly_graph advertises assembly-graph content
        # the same way, added only when a graph is present so ms-only
        # files stay byte-identical.
        has_graphs = bool(assembly_graphs)
        if has_graphs and "opt_assembly_graph" not in feature_list:
            feature_list = feature_list + ["opt_assembly_graph"]
        # M101: opt_sam_tags when any genomic run carries SAM optional
        # fields (a tags channel); tag-less files are unchanged.
        if (any(r.tags and any(r.tags) for r in (genomic_runs or {}).values())
                and "opt_sam_tags" not in feature_list):
            feature_list = feature_list + ["opt_sam_tags"]
        format_version = "1.0"

        # ------------------------------------------------------------------
        # Progress: emit one tick per non-empty section in §5.4 order:
        #   encryption -> provenance -> subjects -> samples ->
        #   references -> image -> identifications ->
        #   quantifications -> runs (MS + NMR + Genomic combined).
        # We start with a baseline (0, total) so consumers can show
        # an immediate determinate bar.
        # ------------------------------------------------------------------
        _section_flags: list[tuple[str, bool]] = [
            ("encryption", False),  # write_minimal never encrypts; placeholder
            ("provenance", bool(provenance)),
            ("subjects", bool(subjects_list)),
            ("samples", bool(samples_list)),
            # references are embedded as part of genomic runs; gate on those
            ("references", has_genomic),
            ("image", image is not None or raman_image is not None
                or ir_image is not None),
            ("identifications", bool(identifications)),
            ("quantifications", bool(quantifications)),
            ("runs", bool(runs) or bool(genomic_runs)),
        ]
        _progress_total = sum(1 for _, present in _section_flags if present)
        _progress_done = 0
        _fire(progress, _progress_done, _progress_total)

        def _section_done(name: str) -> None:
            nonlocal _progress_done
            _progress_done += 1
            _fire(progress, _progress_done, _progress_total)

        # HDF5 fast path keeps the legacy byte layout (fixed-length
        # string attrs, padded compound types) so existing tests and
        # cross-language readers continue to round-trip bit-for-bit.
        if isinstance(provider, str) and provider in ("hdf5", "h5", "h5py"):
            with h5py.File(p, "w") as f:
                io.write_feature_flags(f, format_version, feature_list)
                study = f.create_group("study")
                io.write_fixed_string_attr(study, "title", title)
                io.write_fixed_string_attr(study, "isa_investigation_id", isa_investigation_id)

                # §5.4 ordering: provenance first, then subjects/samples,
                # then references, then image, then identifications /
                # quantifications, then runs (last because they tend to
                # be the largest payload).
                if provenance:
                    _mw._write_provenance(study, provenance)
                    _section_done("provenance")
                if subjects_list:
                    _mw._write_subjects_h5(study, subjects_list)
                    _section_done("subjects")
                if samples_list:
                    _mw._write_samples_h5(study, samples_list)
                    _section_done("samples")

                ms_group = study.create_group("ms_runs")
                io.write_fixed_string_attr(ms_group, "_run_names", ",".join(runs.keys()))
                for rname, run in runs.items():
                    _write_run(ms_group, rname, run)

                nmr_group = study.create_group("nmr_runs")
                io.write_fixed_string_attr(nmr_group, "_run_names", "")

                if has_genomic:
                    # M93 v1.2: embed referenced chromosome sequences at
                    # /study/references/<uri>/ before writing genomic runs,
                    # so the writer's REF_DIFF dispatch can resolve the
                    # md5 attribute back from disk if needed.
                    _gw._embed_references_for_runs(study, genomic_runs)
                    _section_done("references")
                    g_group = study.create_group("genomic_runs")
                    io.write_fixed_string_attr(
                        g_group, "_run_names", ",".join(genomic_runs.keys())
                    )
                    for gname, grun in genomic_runs.items():
                        _write_genomic_run_default(study, g_group, gname, grun)

                if has_graphs:
                    # M98: assembly_graphs subtree, written through the
                    # StorageGroup adapter so the layout is identical
                    # to the provider path below.
                    from ttio.providers.hdf5 import _Group as _Hdf5Group
                    study_sg = _Hdf5Group(study)
                    for ag_name in sorted(assembly_graphs):
                        write_assembly_graph(
                            assembly_graphs[ag_name], ag_name, study_sg)

                if (image is not None or raman_image is not None
                        or ir_image is not None):
                    # Wrap the raw h5py.Group in the package-private adapter
                    # so write_to (which expects a StorageGroup) works
                    # uniformly across both the fast-path and protocol-path branches.
                    from ttio.providers.hdf5 import _Group as _Hdf5Group
                    if image is not None:
                        image.write_to(_Hdf5Group(study))
                    if raman_image is not None:
                        raman_image.write_to(_Hdf5Group(study))
                    if ir_image is not None:
                        ir_image.write_to(_Hdf5Group(study))
                    _section_done("image")

                if identifications:
                    _mw._write_identifications(study, identifications)
                    _section_done("identifications")
                if quantifications:
                    _mw._write_quantifications(study, quantifications)
                    _section_done("quantifications")
                if runs or genomic_runs:
                    _section_done("runs")
            return p

        # Provider-driven write path — Memory / SQLite / Zarr / future.
        # Use the raw ``path`` string rather than ``Path(path)`` because
        # ``Path("memory://x")`` collapses ``//`` → ``/``, breaking the
        # MemoryProvider URL convention.
        owns_provider = False
        if isinstance(provider, str):
            url = str(path)
            sp = open_provider(url, provider=provider, mode="w")
            owns_provider = True
        else:
            sp = provider
        try:
            root = sp.root_group()
            io.write_feature_flags(root, format_version, feature_list)
            study = root.create_group("study")
            io.write_fixed_string_attr(study, "title", title)
            io.write_fixed_string_attr(study, "isa_investigation_id", isa_investigation_id)

            # §5.4 ordering: provenance first, then subjects/samples,
            # then references, then image, then identifications /
            # quantifications, then runs.
            if provenance:
                _mw._write_provenance(study, provenance)
                _section_done("provenance")
            if subjects_list:
                _mw._write_subjects_provider(study, subjects_list)
                _section_done("subjects")
            if samples_list:
                _mw._write_samples_provider(study, samples_list)
                _section_done("samples")

            ms_group = study.create_group("ms_runs")
            io.write_fixed_string_attr(ms_group, "_run_names", ",".join(runs.keys()))
            for rname, run in runs.items():
                _write_run(ms_group, rname, run)

            nmr_group = study.create_group("nmr_runs")
            io.write_fixed_string_attr(nmr_group, "_run_names", "")

            if has_genomic:
                # M93 v1.2: embed referenced chromosome sequences before
                # writing the runs (provider path mirror of the HDF5 path).
                _gw._embed_references_for_runs(study, genomic_runs)
                _section_done("references")
                g_group = study.create_group("genomic_runs")
                io.write_fixed_string_attr(
                    g_group, "_run_names", ",".join(genomic_runs.keys())
                )
                for gname, grun in genomic_runs.items():
                    _write_genomic_run_default(study, g_group, gname, grun)

            if has_graphs:
                # M98: study is already a StorageGroup here.
                for ag_name in sorted(assembly_graphs):
                    write_assembly_graph(
                        assembly_graphs[ag_name], ag_name, study)

            if image is not None:
                image.write_to(study)   # study is already a StorageGroup here
            if raman_image is not None:
                raman_image.write_to(study)   # study is already a StorageGroup here
            if ir_image is not None:
                ir_image.write_to(study)   # study is already a StorageGroup here
            if (image is not None or raman_image is not None
                    or ir_image is not None):
                _section_done("image")

            if identifications:
                _mw._write_identifications(study, identifications)
                _section_done("identifications")
            if quantifications:
                _mw._write_quantifications(study, quantifications)
                _section_done("quantifications")
            if runs or genomic_runs:
                _section_done("runs")
        finally:
            if owns_provider:
                sp.close()
        return p


# ------------------------------------------------------------ writer helpers

@dataclass(slots=True)
class WrittenRun:
    """Simple container passed to :meth:`SpectralDataset.write_minimal`."""

    spectrum_class: str  # "TTIOMassSpectrum" or "TTIONMRSpectrum"
    acquisition_mode: int
    channel_data: dict[str, np.ndarray]  # concatenated signal buffers
    offsets: np.ndarray
    lengths: np.ndarray
    retention_times: np.ndarray
    ms_levels: np.ndarray
    polarities: np.ndarray
    precursor_mzs: np.ndarray
    precursor_charges: np.ndarray
    base_peak_intensities: np.ndarray
    # optional parallel per-spectrum arrays. Writer emits the
    # four datasets only when all four are non-None (schema-gating per
    # the opt_ms2_activation_detail feature flag).
    activation_methods: np.ndarray | None = None
    isolation_target_mzs: np.ndarray | None = None
    isolation_lower_offsets: np.ndarray | None = None
    isolation_upper_offsets: np.ndarray | None = None
    # Optional per-spectrum centroided flag (0 = profile, 1 = centroided).
    # Independent of M74 gating; written to ``spectrum_index/centroideds``.
    centroideds: np.ndarray | None = None
    # M102 imaging-grid position per spectrum, written to
    # ``spectrum_index/pixel_x|pixel_y|pixel_z`` (int32) and flagged
    # ``opt_pixel_coordinates``. All three or none.
    pixel_x: np.ndarray | None = None
    pixel_y: np.ndarray | None = None
    pixel_z: np.ndarray | None = None
    nucleus_type: str = ""
    # Optional NMR solvent label (e.g. "CDCl3", "DMSO-d6"). Empty when
    # not specified or when the run is not NMR. Written as the
    # ``@solvent`` string attribute on the run group. Also reused as the
    # UV-Vis solvent label.
    solvent: str = ""
    # Optional vibrational-spectrum run metadata (IR / Raman / UV-Vis).
    # Written as scalar run-group attributes only when set, so MS/NMR
    # files stay byte-identical. Read back in AcquisitionRun.open and
    # used by _materialize_spectrum to reconstruct the subclass.
    ir_mode: int | None = None              # IRMode value (0/1)
    ir_resolution_cm_inv: float = 0.0
    ir_number_of_scans: int = 0
    raman_excitation_wavelength_nm: float = 0.0
    raman_laser_power_mw: float = 0.0
    raman_integration_time_sec: float = 0.0
    uvvis_path_length_cm: float = 0.0
    provenance_records: list[ProvenanceRecord] = field(default_factory=list)
    # signal compression codec. Valid values are the strings
    # recognised by :func:`ttio._hdf5_io.write_signal_channel` plus
    # the TTIO-level ``"numpress_delta"`` codec, which transforms the
    # float64 buffer into an int64 first-difference array and stores
    # the fixed-point scaling factor on the signal_channels group.
    # "float_delta_zstd" stores each channel as a lossless FDZ1 codec
    # stream (codec id 17, @compression on the dataset).
    # Phase 2 of that spec: "gzip" on a TTIOMassSpectrum run resolves
    # to float_delta_zstd unless opt_disable_float_delta is set;
    # non-MS runs keep the chunked-zlib layout.
    signal_compression: str = "gzip"
    # Opt-out for the MS float_delta_zstd default, same pattern as
    # WrittenGenomicRun.opt_disable_inline_mate_info_v2. Java/ObjC:
    # optDisableFloatDelta.
    opt_disable_float_delta: bool = False
    # optional chromatogram traces for this run. Empty list
    # results in no /chromatograms/ group, preserving byte parity with
    # v0.3 files written by callers that don't supply chromatograms.
    chromatograms: list = field(default_factory=list)  # list[Chromatogram]


def _write_run(parent: h5py.Group, name: str, run: WrittenRun) -> None:
    g = parent.create_group(name)
    io.write_int_attr(g, "acquisition_mode", run.acquisition_mode)
    io.write_int_attr(g, "spectrum_count", int(run.offsets.shape[0]))
    io.write_fixed_string_attr(g, "spectrum_class", run.spectrum_class)
    if run.nucleus_type:
        io.write_fixed_string_attr(g, "nucleus_type", run.nucleus_type)

    if getattr(run, "solvent", ""):
        io.write_fixed_string_attr(g, "solvent", run.solvent)

    # Vibrational-spectrum run metadata (only emitted when set, so MS/NMR
    # files keep byte parity). Mirrored in Java/ObjC readers.
    if run.ir_mode is not None:
        io.write_int_attr(g, "ir_mode", int(run.ir_mode))
    if run.ir_resolution_cm_inv:
        io.write_float_attr(g, "ir_resolution_cm_inv", run.ir_resolution_cm_inv)
    if run.ir_number_of_scans:
        io.write_int_attr(g, "ir_number_of_scans", int(run.ir_number_of_scans))
    if run.raman_excitation_wavelength_nm:
        io.write_float_attr(g, "raman_excitation_wavelength_nm",
                            run.raman_excitation_wavelength_nm)
    if run.raman_laser_power_mw:
        io.write_float_attr(g, "raman_laser_power_mw", run.raman_laser_power_mw)
    if run.raman_integration_time_sec:
        io.write_float_attr(g, "raman_integration_time_sec",
                            run.raman_integration_time_sec)
    if run.uvvis_path_length_cm:
        io.write_float_attr(g, "uvvis_path_length_cm", run.uvvis_path_length_cm)

    if run.provenance_records:
        prov = g.create_group("provenance")
        _mw._write_provenance(prov, run.provenance_records, dataset_name="steps")
        # Legacy @provenance_json mirror so ObjC signature manager and
        # v0.2 readers keep working — matches TTIOAcquisitionRun.m.
        legacy = json.dumps([
            {
                "inputRefs": r.input_refs,
                "software": r.software,
                "parameters": r.parameters,
                "outputRefs": r.output_refs,
                "timestampUnix": int(r.timestamp_unix),
            } for r in run.provenance_records
        ])
        io.write_fixed_string_attr(g, "provenance_json", legacy)

    cfg = g.create_group("instrument_config")
    for field_name in ("manufacturer", "model", "serial_number",
                       "source_type", "analyzer_type", "detector_type"):
        io.write_fixed_string_attr(cfg, field_name, "")

    idx = g.create_group("spectrum_index")
    io.write_int_attr(idx, "count", int(run.offsets.shape[0]))
    # offsets is never written — readers derive from cumsum(lengths).
    columns: list[tuple[str, np.ndarray, str]] = [
        ("lengths", run.lengths, "<u4"),
        ("retention_times", run.retention_times, "<f8"),
        ("ms_levels", run.ms_levels, "<i4"),
        ("polarities", run.polarities, "<i4"),
        ("precursor_mzs", run.precursor_mzs, "<f8"),
        ("precursor_charges", run.precursor_charges, "<i4"),
        ("base_peak_intensities", run.base_peak_intensities, "<f8"),
    ]
    # M74 schema-gating: only emit the four optional columns when all
    # four are supplied. The opt_ms2_activation_detail feature flag is
    # the author-level gate; this block translates that gate into
    # physical column presence on disk.
    m74_cols = (run.activation_methods, run.isolation_target_mzs,
                run.isolation_lower_offsets, run.isolation_upper_offsets)
    if all(c is not None for c in m74_cols):
        columns += [
            ("activation_methods", run.activation_methods, "<i4"),
            ("isolation_target_mzs", run.isolation_target_mzs, "<f8"),
            ("isolation_lower_offsets", run.isolation_lower_offsets, "<f8"),
            ("isolation_upper_offsets", run.isolation_upper_offsets, "<f8"),
        ]
    elif any(c is not None for c in m74_cols):
        raise ValueError(
            "WrittenRun M74 columns must be either all-None or all-set; "
            "partial population is not a valid schema state."
        )
    # Optional centroided column — independent of M74 gating.
    if run.centroideds is not None:
        columns.append(("centroideds", run.centroideds, "<i4"))
    # M102 pixel columns — all three or none.
    pixel_cols = (run.pixel_x, run.pixel_y, run.pixel_z)
    if all(c is not None for c in pixel_cols):
        columns += [
            ("pixel_x", run.pixel_x, "<i4"),
            ("pixel_y", run.pixel_y, "<i4"),
            ("pixel_z", run.pixel_z, "<i4"),
        ]
    elif any(c is not None for c in pixel_cols):
        raise ValueError(
            "WrittenRun pixel_x/pixel_y/pixel_z must be either all-None or "
            "all-set; partial population is not a valid schema state."
        )
    for dname, data, dtype in columns:
        io.write_signal_channel(idx, dname, data.astype(dtype, copy=False),
                                chunk_size=io.DEFAULT_INDEX_CHUNK)

    sig = g.create_group("signal_channels")
    io.write_fixed_string_attr(sig, "channel_names", ",".join(run.channel_data.keys()))
    codec = run.signal_compression
    if (codec == "gzip"
            and not run.opt_disable_float_delta
            and run.spectrum_class == SpectrumKind.MASS.value):
        codec = "float_delta_zstd"
    for cname, buffer in run.channel_data.items():
        if codec == "numpress_delta":
            from ._numpress import encode as _np_encode
            deltas, scale = _np_encode(buffer.astype(np.float64, copy=False))
            ds_name = f"{cname}_values"
            io.write_signal_channel(
                sig, ds_name, deltas, compression="gzip",
            )
            # Per-channel fixed-point attribute, matching the ObjC
            # writer's ``@<chName>_numpress_fixed_point``.
            io.write_int_attr(sig, f"{cname}_numpress_fixed_point", int(scale))
        elif codec == "float_delta_zstd":
            # Codec id 17: the dataset bytes ARE the FDZ1 stream;
            # @compression on the dataset is the dispatch signal and
            # no HDF5 filter is applied (§10.5 discipline).
            from .codecs import float_delta_zstd as _fdz
            from .enums import Compression as _Compression
            stream = _fdz.encode(
                np.ascontiguousarray(buffer, dtype=np.float64))
            io.write_codec_stream_channel(
                sig, f"{cname}_values", stream,
                compression_id=int(_Compression.FLOAT_DELTA_ZSTD))
        else:
            io.write_signal_channel(
                sig, f"{cname}_values",
                buffer.astype("<f8", copy=False),
                compression=codec,
            )

    # chromatograms
    if run.chromatograms:
        from .acquisition_run import write_chromatograms_to_run_group
        write_chromatograms_to_run_group(g, run.chromatograms)


