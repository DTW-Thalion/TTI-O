"""Dataset walker — produces transport-stream events from a
:class:`SpectralDataset` filtered through an :class:`AUFilter`.

This is the Python equivalent of the ObjC ``TTIODatasetWalker`` and
Java ``DatasetWalker``. Where ObjC/Java use a visitor protocol, the
Python form is a generator that yields :class:`WalkerEvent`
instances — more idiomatic and trivially composable with ``async``
producers via ``async for`` adapters.

Used by:
  * :func:`ttio.transport.server._emit_stream` — encodes each event
    as a transport packet.
  * Workbench Server S3 (out-of-tree) — drives binary, stats-only,
    and stats-with-payload download modes from the same iteration.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Optional, TYPE_CHECKING

from ._common import genomic_channel_names
from .codec import _iter_genomic_run_access_units, _spectrum_to_access_unit
from .filters import AUFilter
from .packets import AccessUnit
from ..spectral_dataset import SpectralDataset

if TYPE_CHECKING:  # pragma: no cover - import-cycle avoidance
    from ..genomic.reference_import import ReferenceImport
    from ..identification import Identification
    from ..ir_image import IRImage
    from ..ms_image import MSImage
    from ..provenance import ProvenanceRecord
    from ..quantification import Quantification
    from ..raman_image import RamanImage
    from ..sample import Sample
    from ..subject import Subject


# ── Event types ────────────────────────────────────────────────────


class WalkerEvent:
    """Marker base class for walker events. Use ``isinstance`` /
    pattern matching in consumers."""


@dataclass(frozen=True, slots=True)
class StreamHeaderEvent(WalkerEvent):
    format_version: str
    title: str
    isa_investigation: str
    features: list[str]
    n_datasets: int


@dataclass(frozen=True, slots=True)
class DatasetHeaderEvent(WalkerEvent):
    dataset_id: int
    name: str
    acquisition_mode: int
    spectrum_class: str
    channel_names: list[str]
    instrument_json: str
    expected_au_count: int


@dataclass(frozen=True, slots=True)
class AccessUnitEvent(WalkerEvent):
    au: AccessUnit
    dataset_id: int
    au_sequence: int


@dataclass(frozen=True, slots=True)
class EndOfDatasetEvent(WalkerEvent):
    dataset_id: int
    final_au_sequence: int


@dataclass(frozen=True, slots=True)
class EndOfStreamEvent(WalkerEvent):
    pass


# ── v0.11 §5.4 prelude events (#141) ───────────────────────────────


@dataclass(frozen=True, slots=True)
class EncryptionAlgorithmEvent(WalkerEvent):
    """§5.4.1 — dataset-level ``@encrypted`` algorithm name."""
    algorithm: str


@dataclass(frozen=True, slots=True)
class DatasetProvenanceEvent(WalkerEvent):
    """§5.4.2 — dataset-level provenance chain."""
    records: list  # list[ProvenanceRecord]


@dataclass(frozen=True, slots=True)
class SubjectMetadataEvent(WalkerEvent):
    """§5.4.3 — Subject rows."""
    rows: list  # list[Subject]


@dataclass(frozen=True, slots=True)
class SampleMetadataEvent(WalkerEvent):
    """§5.4.3 — Sample rows."""
    rows: list  # list[Sample]


@dataclass(frozen=True, slots=True)
class ReferenceGroupEvent(WalkerEvent):
    """§5.4.4 — one embedded reference (single ReferenceImport)."""
    reference: object  # ReferenceImport


@dataclass(frozen=True, slots=True)
class ImageEvent(WalkerEvent):
    """§5.4.5 — embedded MSImage cube."""
    image: object  # MSImage


@dataclass(frozen=True, slots=True)
class RamanImageEvent(WalkerEvent):
    """§5.4.5 — embedded RamanImage cube."""
    image: object  # RamanImage


@dataclass(frozen=True, slots=True)
class IRImageEvent(WalkerEvent):
    """§5.4.5 — embedded IRImage cube."""
    image: object  # IRImage


@dataclass(frozen=True, slots=True)
class IdentificationsTableEvent(WalkerEvent):
    """§5.4.6 — Identification rows."""
    rows: list  # list[Identification]


@dataclass(frozen=True, slots=True)
class QuantificationsTableEvent(WalkerEvent):
    """§5.4.6 — Quantification rows."""
    rows: list  # list[Quantification]


# ── Walker ─────────────────────────────────────────────────────────


def walk_dataset(
    dataset: SpectralDataset,
    filter: Optional[AUFilter] = None,
) -> Iterator[WalkerEvent]:
    """Generator that yields the canonical transport-stream event
    sequence for the given dataset:

      1. :class:`StreamHeaderEvent` once.
      1a. v0.11 §5.4 prelude events (when populated):
            :class:`EncryptionAlgorithmEvent`,
            :class:`DatasetProvenanceEvent`,
            :class:`SubjectMetadataEvent`,
            :class:`SampleMetadataEvent`,
            :class:`ReferenceGroupEvent` (one per reference, sorted
            by URI key),
            :class:`ImageEvent`,
            :class:`RamanImageEvent`,
            :class:`IRImageEvent`,
            :class:`IdentificationsTableEvent`,
            :class:`QuantificationsTableEvent`.
      2. :class:`DatasetHeaderEvent` per matched dataset.
      3. :class:`AccessUnitEvent` per matched AU (spectral + genomic).
      4. :class:`EndOfDatasetEvent` per matched dataset.
      5. :class:`EndOfStreamEvent` once.

    Iteration order matches :meth:`TransportWriter.write_dataset` —
    the byte form emitted by encoding the events one-by-one is
    identical to the unfiltered writer output (mod filter cuts).

    :param dataset: A :class:`SpectralDataset` to walk.
    :param filter:  Optional :class:`AUFilter`. ``None`` emits every
                    AU in dataset iteration order.
    """
    from .codec import _instrument_config_json, _genomic_run_metadata_json

    flt = filter or AUFilter()
    ms_runs = list(dataset.all_runs.items())
    genomic_runs = list(getattr(dataset, "genomic_runs", {}).items())
    features = list(dataset.feature_flags.features)

    # 1. StreamHeader
    yield StreamHeaderEvent(
        format_version="1.2",
        title=dataset.title or "",
        isa_investigation=dataset.isa_investigation_id or "",
        features=features,
        n_datasets=len(ms_runs) + len(genomic_runs),
    )

    # 1a. v0.11 §5.4 prelude — match TransportWriter.write_dataset
    # ordering verbatim. Each gate uses the same "populated?" check the
    # writer uses, so the on-wire packet sequence is byte-identical.
    # See codec.py::TransportWriter.write_dataset (lines ~1180-1225).
    refs = getattr(dataset, "references", {}) or {}
    try:
        dataset_provenance = list(dataset.provenance())
    except Exception:  # pragma: no cover - defensive
        dataset_provenance = []
    from ..enums import ImageKind
    dataset_image = dataset.image_for_kind(ImageKind.MS)
    dataset_raman_image = dataset.image_for_kind(ImageKind.RAMAN)
    dataset_ir_image = dataset.image_for_kind(ImageKind.IR)
    try:
        dataset_identifications = list(dataset.identifications())
    except Exception:  # pragma: no cover - defensive
        dataset_identifications = []
    try:
        dataset_quantifications = list(dataset.quantifications())
    except Exception:  # pragma: no cover - defensive
        dataset_quantifications = []
    try:
        dataset_subjects = list(getattr(dataset, "subjects", []) or [])
    except Exception:  # pragma: no cover - defensive
        dataset_subjects = []
    try:
        dataset_samples = list(getattr(dataset, "samples", []) or [])
    except Exception:  # pragma: no cover - defensive
        dataset_samples = []

    # §5.4.1 ENCRYPTION_ALGORITHM
    if getattr(dataset, "is_encrypted", False):
        algo = getattr(dataset, "encrypted_algorithm", "") or ""
        if algo:
            yield EncryptionAlgorithmEvent(algorithm=algo)
    # §5.4.2 DATASET_PROVENANCE
    if dataset_provenance:
        yield DatasetProvenanceEvent(records=dataset_provenance)
    # §5.4.3 SUBJECT_METADATA → SAMPLE_METADATA (subjects first so a
    # reader sees the soft-FK target ahead of any referencing sample).
    if dataset_subjects:
        yield SubjectMetadataEvent(rows=dataset_subjects)
    if dataset_samples:
        yield SampleMetadataEvent(rows=dataset_samples)
    # §5.4.4 reference groups — sorted by URI key (matches ObjC walker
    # which sorts the NSDictionary keys for deterministic order).
    for uri in sorted(refs.keys()):
        yield ReferenceGroupEvent(reference=refs[uri])
    # §5.4.5 image cubes — MS → Raman → IR.
    if dataset_image is not None:
        yield ImageEvent(image=dataset_image)
    if dataset_raman_image is not None:
        yield RamanImageEvent(image=dataset_raman_image)
    if dataset_ir_image is not None:
        yield IRImageEvent(image=dataset_ir_image)
    # §5.4.6 IDENTIFICATIONS_TABLE → QUANTIFICATIONS_TABLE
    if dataset_identifications:
        yield IdentificationsTableEvent(rows=dataset_identifications)
    if dataset_quantifications:
        yield QuantificationsTableEvent(rows=dataset_quantifications)

    # 2. DatasetHeaders — spectral runs first (dataset_ids 1..N),
    # then genomic runs (dataset_ids N+1..N+M). Matches
    # TransportWriter.write_dataset ordering.
    for i, (name, run) in enumerate(ms_runs, start=1):
        if flt.dataset_id is not None and i != flt.dataset_id:
            continue
        yield DatasetHeaderEvent(
            dataset_id=i,
            name=name,
            acquisition_mode=int(run.acquisition_mode),
            spectrum_class=run.spectrum_class,
            channel_names=list(run.channel_names),
            instrument_json=_instrument_config_json(run),
            expected_au_count=len(run),
        )
    for j, (name, grun) in enumerate(genomic_runs, start=len(ms_runs) + 1):
        if flt.dataset_id is not None and j != flt.dataset_id:
            continue
        yield DatasetHeaderEvent(
            dataset_id=j,
            name=name,
            acquisition_mode=int(grun.acquisition_mode),
            spectrum_class="TTIOGenomicRead",
            channel_names=genomic_channel_names(grun),
            instrument_json=_genomic_run_metadata_json(grun),
            expected_au_count=len(grun),
        )

    # 3. AccessUnits + 4. EndOfDataset, interleaved per dataset to
    # match TransportWriter.write_dataset: AUs(run1) → EOD(run1) →
    # AUs(run2) → EOD(run2) → … → AUs(genomic_runK) → EOD(genomic_runK).
    # Spectral runs first (dataset_ids 1..N), then genomic runs
    # (dataset_ids N+1..N+M).
    emitted = 0
    max_au = flt.max_au
    for i, (_name, run) in enumerate(ms_runs, start=1):
        if flt.dataset_id is not None and i != flt.dataset_id:
            continue
        for j, spectrum in enumerate(run):
            if max_au is not None and emitted >= max_au:
                break
            au = _spectrum_to_access_unit(spectrum, run)
            if not flt.matches(au, i):
                continue
            yield AccessUnitEvent(au=au, dataset_id=i, au_sequence=j)
            emitted += 1
        yield EndOfDatasetEvent(dataset_id=i, final_au_sequence=len(run))
        if max_au is not None and emitted >= max_au:
            break
    for j, (_name, grun) in enumerate(genomic_runs, start=len(ms_runs) + 1):
        if flt.dataset_id is not None and j != flt.dataset_id:
            continue
        for au_seq, au in _iter_genomic_run_access_units(grun):
            if max_au is not None and emitted >= max_au:
                break
            if not flt.matches(au, j):
                continue
            yield AccessUnitEvent(au=au, dataset_id=j, au_sequence=au_seq)
            emitted += 1
        yield EndOfDatasetEvent(dataset_id=j, final_au_sequence=len(grun))
        if max_au is not None and emitted >= max_au:
            break

    # 5. EndOfStream
    yield EndOfStreamEvent()
