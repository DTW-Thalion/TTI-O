"""``ImportedDataset`` — the normalized in-memory draft every importer
produces. The single call site of :meth:`SpectralDataset.write_minimal`,
collapsing the per-format adapter normalization the registry used to do.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from ..identification import Identification
from ..provenance import ProvenanceRecord
from ..quantification import Quantification
from ..spectral_dataset import SpectralDataset, WrittenRun

if TYPE_CHECKING:
    from ..written_genomic_run import WrittenGenomicRun


def _flag_sam_tags(ds: SpectralDataset, streams: dict) -> None:
    """Add ``opt_sam_tags`` once a streamed genomic run turned out to
    carry SAM tags (M101); the static write could not know."""
    from .. import _hdf5_io as io
    runs = ds.study_group.open_group("genomic_runs")
    names = [getattr(src, "name", key) for key, src in streams.items()]
    if not any(runs.has_child(name)
               and runs.open_group(name).open_group("signal_channels").has_child("tags")
               for name in names):
        return
    root = ds.provider.root_group()
    version, features = io.read_feature_flags(root)
    if "opt_sam_tags" not in features:
        io.write_feature_flags(root, version, list(features) + ["opt_sam_tags"])


@dataclass(slots=True)
class ImportedDataset:
    """In-memory bundle of built runs + dataset metadata, ready to write."""

    title: str = ""
    isa_investigation_id: str = ""
    runs: dict[str, WrittenRun] = field(default_factory=dict)
    genomic_runs: dict[str, "WrittenGenomicRun"] = field(default_factory=dict)
    identifications: list[Identification] = field(default_factory=list)
    quantifications: list[Quantification] = field(default_factory=list)
    provenance: list[ProvenanceRecord] = field(default_factory=list)
    image: object | None = None
    raman_image: object | None = None
    ir_image: object | None = None
    subjects: list = field(default_factory=list)
    samples: list = field(default_factory=list)
    # Streaming sources (GenomicStreamSource / SpectralStreamSource):
    # written after write_minimal through the stream writers, so a run
    # of any size never sits in memory whole.
    genomic_streams: dict = field(default_factory=dict)
    spectral_streams: dict = field(default_factory=dict)

    def write(self, path: str | Path, *, features: list[str] | None = None,
              provider: str = "hdf5", progress=None) -> Path:
        out = self._write_static(path, features=features, provider=provider, progress=progress)
        if self.genomic_streams or self.spectral_streams:
            with SpectralDataset.open(out, writable=True, provider=provider) as ds:
                for src in self.genomic_streams.values():
                    src.write_into(ds.study_group, progress=progress)
                for src in self.spectral_streams.values():
                    src.write_into(ds.study_group, progress=progress)
                if self.genomic_streams:
                    _flag_sam_tags(ds, self.genomic_streams)
        return out

    def _write_static(self, path: str | Path, *, features: list[str] | None = None,
                      provider: str = "hdf5", progress=None) -> Path:
        return SpectralDataset.write_minimal(
            path,
            title=self.title or "imported",
            isa_investigation_id=self.isa_investigation_id,
            runs=self.runs,
            genomic_runs=self.genomic_runs or None,
            identifications=self.identifications or None,
            quantifications=self.quantifications or None,
            provenance=self.provenance or None,
            features=features,
            provider=provider,
            image=self.image,
            raman_image=self.raman_image,
            ir_image=self.ir_image,
            subjects=self.subjects or None,
            samples=self.samples or None,
            progress=progress,
        )
