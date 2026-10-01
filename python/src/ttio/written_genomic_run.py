"""WrittenGenomicRun — write-side container for a single genomic run.

Passed to :meth:`ttio.spectral_dataset.SpectralDataset.write_minimal`
via the ``genomic_runs`` parameter. Genomic analogue of
:class:`ttio.spectral_dataset.WrittenRun`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .enums import Compression
from .provenance import ProvenanceRecord


@dataclass(slots=True)
class BulkV2Blobs:
    """Verbatim v2 codec blobs for direct on-disk write.

    Set on a :class:`WrittenGenomicRun` to bypass the v2 codec encode
    step in ``write_minimal`` and write the blob bytes directly to
    the matching HDF5 paths. Used by the transport bulk-mode receiver
    (see ``docs/transport-spec.md`` §6.4).

    Each field is independently optional. When ``mate_info_blob`` is
    set the writer also requires ``mate_info_chrom_names``.
    ``ref_diff_blob`` requires ``ref_diff_reference_uri`` to validate
    against the run's ``reference_uri`` attribute.

    Cross-language: ObjC ``TTIOBulkV2Blobs`` · Java ``BulkV2Blobs``.
    """

    mate_info_blob: bytes | None = None
    mate_info_chrom_names: list[str] | None = None
    name_tok_blob: bytes | None = None
    ref_diff_blob: bytes | None = None
    ref_diff_reference_uri: str | None = None

    def __post_init__(self) -> None:
        if self.mate_info_blob is not None and self.mate_info_chrom_names is None:
            raise ValueError(
                "BulkV2Blobs.mate_info_blob requires mate_info_chrom_names"
            )
        if self.ref_diff_blob is not None and self.ref_diff_reference_uri is None:
            raise ValueError(
                "BulkV2Blobs.ref_diff_blob requires ref_diff_reference_uri"
            )


@dataclass(slots=True)
class WrittenGenomicRun:
    """Data container for writing a genomic run via SpectralDataset."""

    acquisition_mode: int             # AcquisitionMode.GENOMIC_WGS or _WES (.value)
    reference_uri: str                # e.g., "GRCh38.p14"
    platform: str                     # e.g., "ILLUMINA"
    sample_name: str                  # e.g., "NA12878"

    # Per-read parallel arrays (all length == read_count)
    positions: np.ndarray             # int64
    mapping_qualities: np.ndarray     # uint8
    flags: np.ndarray                 # uint32

    # Concatenated signal data
    sequences: np.ndarray             # uint8 — one ASCII byte per base (M82)
    qualities: np.ndarray             # uint8 — Phred scores, concatenated

    # Per-read offsets into sequences/qualities
    offsets: np.ndarray               # uint64
    lengths: np.ndarray               # uint32

    # Per-read variable-length fields
    cigars: list[str]                 # one CIGAR string per read
    read_names: list[str]             # one read name per read

    # Mate info (per-read)
    mate_chromosomes: list[str]
    mate_positions: np.ndarray        # int64 (-1 if unpaired)
    template_lengths: np.ndarray      # int32 (0 if unpaired)

    # Chromosomes (per-read, for the index)
    chromosomes: list[str]

    # Optional
    provenance_records: list[ProvenanceRecord] = field(default_factory=list)
    signal_compression: str = "gzip"  # "gzip" → ZLIB; "none" → NONE

    # Removes the V5 sequence-context strategies from the qualities
    # auto-tune set (float-delta and mate_info opt-out pattern; spec
    # 2.4). Java/ObjC: optDisableQualitiesV5.
    opt_disable_qualities_v5: bool = False

    # blocks_v1 (streaming): a preassigned chromosome name -> id map
    # shared across the blocks of one run. When set, genomic_index and
    # mate_info use these ids and extend the map in place for names not
    # yet present; when None (whole-run write) ids are assigned in
    # first-seen order from this run alone.
    chrom_name_to_id: dict[str, int] | None = None

    # Write this run in the v1.8 whole-channel layout instead of the
    # blocks_v1 default (readers up to v1.8.0 need it). Memory is then
    # unbounded in the run size.
    opt_legacy_whole_channel: bool = False

    # Precomputed reference MD5 (the sorted-name concatenation digest,
    # see genomic.reference_import.compute_reference_md5). The stream
    # writer computes it once per run so a large reference is not
    # digested per block; None means compute from reference_chrom_seqs.
    reference_md5: bytes | None = None

    # per-channel codec opt-in. Maps channel name to a TTI-O
    # internal codec id. Only "sequences" and "qualities" are
    # accepted; only RANS_ORDER0, RANS_ORDER1, BASE_PACK are
    # accepted as codec values. Channels not in this dict use the
    # existing signal_compression string path.
    signal_codec_overrides: dict[str, Compression] = field(default_factory=dict)

    # M93 v1.2 — reference embed for the REF_DIFF_V2 codec on the
    # ``sequences`` channel. When ``embed_reference=True`` AND
    # ``reference_chrom_seqs`` is provided AND the REF_DIFF_V2
    # default applies (signal_compression="gzip", reference present),
    # the writer embeds the chromosome sequences at
    # ``/study/references/<reference_uri>/`` in the output file.
    # When ``embed_reference=False`` (the default since L3, Task #82
    # Phase B.1, 2026-05-01), the writer records ``reference_uri``
    # and ``reference_md5`` only and expects the reader to resolve
    # via REF_PATH or an explicit external path.
    #
    # The default flipped to ``False`` to match CRAM 3.1's default
    # (external reference) and to drop the ~10 MB chr22 reference
    # blob from the v1.2.0 chr22 benchmark; users who want
    # self-contained files set ``embed_reference=True`` explicitly.
    #
    # v1.0 reset (Phase 2c): the v1 REF_DIFF (codec id 9) writer was
    # removed; sequences default to REF_DIFF_V2 (codec id 14).
    embed_reference: bool = False

    # Mapping ``chromosome_name → uppercase ACGTN bytes``, supplied at
    # write time for any chromosome that has at least one read aligned
    # to it. Required when ``embed_reference=True`` so the v1.8
    # REF_DIFF_V2 writer can resolve the per-read diff against the
    # per-chromosome reference; if absent, the writer falls back to
    # BASE_PACK on the sequences channel (Q5b = C).
    reference_chrom_seqs: dict[str, bytes] | None = None

    # External reference path stamped into the file's metadata for
    # decoder fallback when the embedded reference is absent. The
    # writer never reads this path; it is metadata only.
    external_reference_path: Path | None = None

    # Phase 2c-T : verbatim v2 codec blobs from the transport
    # bulk-mode receiver. When set, ``write_minimal`` writes the
    # blob bytes directly to the matching HDF5 paths and SKIPS the
    # corresponding v2 codec encode step. This is the only
    # mechanism that preserves ``mate_chromosome`` SAM sentinels
    # (``=``, ``""``) byte-for-byte across transport.
    bulk_v2_blobs: BulkV2Blobs | None = None

    # M97 — assembly read role, persisted as the ``@read_role`` UTF-8
    # run attribute when non-None. Recognised values: ``hifi``,
    # ``ont_ul``, ``hic_r1``, ``hic_r2``, ``parental_maternal``,
    # ``parental_paternal``, ``illumina_wgs``; other strings are
    # stored unchecked. Java/ObjC: readRole.
    read_role: str | None = None

    # M101 — SAM optional fields per read (columns 12+, tab-joined as
    # samtools prints them; "" for a read without tags). None, or a list
    # whose entries are all "", writes no ``tags`` channel. Stored with
    # SAM_TAGS (codec id 18); MD:Z / NM:i are recomputed from the
    # reference when the sequences channel is REF_DIFF_V2. Java/ObjC:
    # tags.
    tags: list[str] | None = None

    # M97 — REF_DIFF_V2 slice byte budget: a slice closes before the
    # read that would push it past this many bases (the 10,000-read
    # cap still applies). 0 keeps the fixed-count rule. Writer policy
    # only — wire format and decoder unchanged. Java/ObjC:
    # refDiffSliceBytes.
    ref_diff_slice_bytes: int = 0
