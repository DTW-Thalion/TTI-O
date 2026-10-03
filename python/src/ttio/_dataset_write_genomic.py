"""Genomic-run write helpers extracted from spectral_dataset.py (P3.10).

Pure code-movement extraction (P3.10 PR-1): the genomic-run write
subsystem (reference embed, REF_DIFF_V2 / FQZCOMP_NX16_Z / NAME_TOKENIZED_V2 /
MATE_INLINE_V2 dispatch, bulk-verbatim writers). No behaviour change.
"""
from __future__ import annotations

from typing import Mapping

import h5py
import numpy as np

from . import _hdf5_io as io
from .genomic import packed_reference
from .genomic.reference_import import ReferenceImport
from .written_genomic_run import WrittenGenomicRun
from ._dataset_write_metadata import _write_provenance


def _any_v1_5_codec(
    genomic_runs: "Mapping[str, WrittenGenomicRun] | None",
) -> bool:
    """Return True if any run carries a v1.5 codec (REF_DIFF_V2 or FQZCOMP_NX16_Z).

    Used to gate the format-version bump from 1.4 → 1.5: only files that
    actually exercise an M93+ codec get the new version string, so
    M82-only writes preserve byte-parity with existing fixtures.

    v1.0 reset (Phase 2c): the v1 REF_DIFF (codec id 9) implementation was
    removed; only REF_DIFF_V2 (codec id 14) and FQZCOMP_NX16_Z (codec id 12)
    register as v1.5 codecs for the format-version gate. FQZCOMP_NX16_Z
    carries its sibling-channel metadata (``read_lengths`` / ``revcomp_flags``)
    inside the codec wire format, so it does not need an embedded
    reference (``needs_embedded_reference`` is False) — but it IS a v1.5
    codec for the gate.
    """
    if not genomic_runs:
        return False
    from .enums import Compression as _Compression
    _V1_5_CODECS = frozenset({
        _Compression.REF_DIFF_V2,      # v1.8 #11
        _Compression.FQZCOMP_NX16_Z,   # M94.Z (CRAM-mimic rANS-Nx16)
    })
    for run in genomic_runs.values():
        for codec in run.signal_codec_overrides.values():
            try:
                ce = _Compression(codec)
            except ValueError:
                continue
            if ce in _V1_5_CODECS:
                return True
    return False


# Back-compat alias — pre-M94 callers used this name.
_any_context_aware_codec = _any_v1_5_codec


def _reference_md5_for_run(run: WrittenGenomicRun) -> bytes:
    """MD5 of concatenated chromosome sequences, in sorted-name order.

    Mirrors the on-disk ``@md5`` attribute computation. Returns an empty
    digest when ``reference_chrom_seqs`` is absent.
    """
    import hashlib
    if run.reference_md5 is not None:
        return run.reference_md5
    if run.reference_chrom_seqs is None:
        return b""
    set_md5 = getattr(run.reference_chrom_seqs, "set_md5", None)
    if callable(set_md5):  # LazyReference: cached whole-FASTA digest
        return set_md5()
    md5 = hashlib.md5()
    for chrom_name in sorted(run.reference_chrom_seqs):
        md5.update(run.reference_chrom_seqs[chrom_name])
    return md5.digest()


def _load_references_provider(study) -> dict[str, ReferenceImport]:
    """Read ``/study/references/<uri>/`` via the StorageGroup protocol.

    Inverse-side helper for :func:`_embed_references_for_runs`. Returns
    an empty dict when the ``references`` sub-group is absent, mirroring
    Java's ``Map.of()`` empty-map contract.
    """
    if not study.has_child("references"):
        return {}
    refs_grp = study.open_group("references")
    out: dict[str, ReferenceImport] = {}
    try:
        for uri in refs_grp.child_names():
            sub = refs_grp.open_group(uri)
            try:
                out[uri] = ReferenceImport.read_from_group(sub)
            finally:
                sub.close()
    finally:
        refs_grp.close()
    return out


def _embed_references_for_runs(
    study, genomic_runs: "Mapping[str, WrittenGenomicRun]",
) -> None:
    """Embed each unique reference (by ``reference_uri``) once at
    ``/study/references/<reference_uri>/``.

    Only runs that have ``embed_reference=True`` AND a context-aware
    codec override on ``sequences`` AND non-None ``reference_chrom_seqs``
    contribute; the dedup key is ``reference_uri``. When the same URI
    carries two different MD5s across runs, raises :class:`ValueError`
    (Q6 = C, single source of truth per file).

    Accepts either a raw ``h5py.Group`` (HDF5 fast path) or a
    :class:`StorageGroup` (provider path). Internally normalises to
    ``StorageGroup``.
    """
    from .codecs._registry import CODEC_REGISTRY
    from .enums import Compression as _Compression
    from .providers.hdf5 import _Group as _H5Group

    from .codecs import ref_diff_v2 as _rdv2_meta

    needs_embed: dict[str, tuple[bytes, dict[str, bytes]]] = {}
    for run in genomic_runs.values():
        if not run.embed_reference:
            continue
        if run.reference_chrom_seqs is None:
            continue
        # Embed if a context-aware codec override is set on this run,
        # OR if the v1.8 REF_DIFF_V2 default path will be used (when the
        # native lib is available).
        _has_context_aware_override = any(
            getattr(CODEC_REGISTRY.get(_Compression(c)), "needs_embedded_reference", False)
            for c in run.signal_codec_overrides.values()
            if _is_valid_compression(c)
        )
        _uses_ref_diff_v2_default = _rdv2_meta.HAVE_NATIVE_LIB
        if not (_has_context_aware_override or _uses_ref_diff_v2_default):
            continue
        md5 = _reference_md5_for_run(run)
        if run.reference_uri in needs_embed:
            existing_md5, _ = needs_embed[run.reference_uri]
            if existing_md5 != md5:
                raise ValueError(
                    f"reference_uri {run.reference_uri!r} carries two "
                    "different MD5s across runs in this dataset: "
                    f"{existing_md5.hex()} vs {md5.hex()} — same URI "
                    "cannot map to two different reference contents."
                )
            continue
        needs_embed[run.reference_uri] = (md5, dict(run.reference_chrom_seqs))

    if not needs_embed:
        return

    # Normalise study to a StorageGroup so create_group / create_dataset
    # have a single API surface.
    if isinstance(study, h5py.Group):
        study_sg = _H5Group(study)
    else:
        study_sg = study

    if study_sg.has_child("references"):
        refs_grp = study_sg.open_group("references")
    else:
        refs_grp = study_sg.create_group("references")

    for uri, (md5, chrom_seqs) in needs_embed.items():
        if refs_grp.has_child(uri):
            existing = refs_grp.open_group(uri)
            existing_md5_hex = io.read_string_attr(existing, "md5") or ""
            if existing_md5_hex != md5.hex():
                raise ValueError(
                    f"reference_uri {uri!r} already embedded with a "
                    f"different MD5 ({existing_md5_hex!r} != "
                    f"{md5.hex()!r}); same URI cannot map to two "
                    "different reference contents in one file."
                )
            continue
        ref_grp = refs_grp.create_group(uri)
        io.write_fixed_string_attr(ref_grp, "md5", md5.hex())
        io.write_fixed_string_attr(ref_grp, "reference_uri", uri)
        chroms_grp = ref_grp.create_group("chromosomes")
        for chrom_name in sorted(chrom_seqs):
            seq = chrom_seqs[chrom_name]
            c = chroms_grp.create_group(chrom_name)
            io.write_int_attr(c, "length", len(seq))
            packed_reference.write_chromosome_dataset(c, seq)


def _is_valid_compression(value: object) -> bool:
    from .enums import Compression as _Compression
    try:
        _Compression(value)
        return True
    except ValueError:
        return False


def _write_sequences_ref_diff_v2(sc, run: WrittenGenomicRun) -> None:
    """Write the ``sequences`` channel through the REF_DIFF_V2 codec.

    v1.0 reset (Phase 2c): the v1 REF_DIFF (codec id 9) writer was
    removed. The v2 path (codec id 14) is now the only reference-diff
    sequences writer.

    Eligibility: requires libttio_rans loadable, a single-chromosome
    run and a reference present. When any precondition fails, falls
    back to BASE_PACK on a flat dataset (Q5b = C) — same fallback
    semantics as the original v1.5 REF_DIFF path. The fallback uses
    the canonical, codec-free sequences dataset layout. Unmapped reads
    (``cigar == "*"``, for example a read placed on its mate's
    chromosome) are carried by the codec itself since v1.9: their bases
    go to the soft-clip substream and their lengths to the slice's UL
    substream (docs/codecs/ref_diff_v2.md section 4.4).

    **Single-chromosome limitation (v1.8 first pass):** REF_DIFF_V2
    requires all reads aligned to a single chromosome. Multi-chromosome
    runs raise :class:`ValueError`.
    """
    from .codecs import ref_diff_v2 as _rdv2
    from .codecs._registry import CODEC_REGISTRY
    from .codecs._context import CodecContext, DecodedChannel
    from .codecs.base_pack import encode as _base_pack_encode
    from .enums import Compression as _Compression, Precision as _Precision

    # Resolve the reference sequence for this run.
    chrom_seq: bytes | None = None
    if run.reference_chrom_seqs is not None:
        unique_chroms = set(run.chromosomes)
        if len(unique_chroms) == 0:
            chrom_seq = None
        elif len(unique_chroms) > 1:
            raise ValueError(
                "REF_DIFF_V2 v1.8 supports single-chromosome runs only; "
                f"this run carries reads on chromosomes {sorted(unique_chroms)}. "
                "Multi-chromosome support is a follow-up — split into "
                "per-chromosome runs as a workaround."
            )
        else:
            chrom = next(iter(unique_chroms))
            chrom_seq = run.reference_chrom_seqs.get(chrom)

    raw_bytes = bytes(run.sequences.tobytes())

    use_v2 = (
        _rdv2.HAVE_NATIVE_LIB
        and chrom_seq is not None
    )

    if use_v2:
        # v1.8 path: encode via ref_diff_v2 and write as a GROUP with
        # a refdiff_v2 child dataset (@compression = 14).
        positions = np.asarray(run.positions, dtype=np.int64)
        n = len(run.cigars)
        # Build n_reads+1 offsets array from run.offsets (n entries)
        # and run.lengths (n entries): append the total base count.
        offsets_arr = np.asarray(run.offsets, dtype=np.uint64)
        if offsets_arr.shape[0] == n:
            total_len = int(offsets_arr[-1]) + int(run.lengths[-1]) if n > 0 else 0
            offsets_arr = np.append(offsets_arr, np.uint64(total_len))
        elif offsets_arr.shape[0] != n + 1:
            raise ValueError(
                f"run.offsets must have n_reads or n_reads+1 entries; "
                f"got {offsets_arr.shape[0]} for n={n}"
            )

        md5 = _reference_md5_for_run(run)
        sequences_arr = np.asarray(run.sequences, dtype=np.uint8)
        # Codec-registry refactor (Task 5c): route REF_DIFF_V2 encode
        # through the registry. Inputs that get written into the blob
        # header (offsets/reference/md5/uri/cigars/positions) ride on the
        # encode-only CodecContext fields. The codec returns a GROUP
        # layout ({"refdiff_v2": blob}); we materialise the `sequences`
        # group + child + @compression below — byte-identical to the
        # prior direct ref_diff_v2.encode call.
        encoded_channel = CODEC_REGISTRY[_Compression.REF_DIFF_V2].encode(
            DecodedChannel.of_bytes(bytes(sequences_arr.tobytes())),
            CodecContext(
                offsets=offsets_arr,
                positions=positions,
                cigar_strings=list(run.cigars),
                reference=chrom_seq,
                reference_md5=md5,
                reference_uri=run.reference_uri,
                reads_per_slice=10_000,
                slice_bytes=int(run.ref_diff_slice_bytes or 0),
            ),
        )
        seq_group = sc.create_group("sequences")
        for child_name, blob in encoded_channel.group_children.items():
            arr = np.frombuffer(blob, dtype=np.uint8)
            ds = seq_group.create_dataset(
                child_name,
                _Precision.UINT8,
                length=int(arr.shape[0]),
                chunk_size=io.DEFAULT_SIGNAL_CHUNK,
                compression=_Compression.NONE,
            )
            ds.write(arr)
            io.write_int_attr(
                ds, "compression", int(_Compression.REF_DIFF_V2), dtype="<u1")
        return

    # Fallback: flat dataset with BASE_PACK (Q5b = C).
    encoded = _base_pack_encode(raw_bytes)
    codec_id = int(_Compression.BASE_PACK)

    arr = np.frombuffer(encoded, dtype=np.uint8)
    ds = sc.create_dataset(
        "sequences",
        _Precision.UINT8,
        length=int(arr.shape[0]),
        chunk_size=io.DEFAULT_SIGNAL_CHUNK,
        compression=_Compression.NONE,
    )
    ds.write(arr)
    io.write_int_attr(ds, "compression", codec_id, dtype="<u1")


# M94.Z — FQZCOMP_NX16_Z quality codec. Same dispatch pattern as
# _write_sequences_ref_diff but for the qualities channel.
SAM_REVERSE_FLAG = 16


def _write_sequences_seq_cm(sc, run: WrittenGenomicRun) -> None:
    """Write the ``sequences`` channel through SEQ_CM (codec id 19, M103).

    Context-aware like FQZCOMP_NX16_Z: the encoder takes the read
    lengths, which the run stores in genomic_index, so the registry is
    given a CodecContext carrying them.
    """
    from .codecs._registry import CODEC_REGISTRY
    from .codecs._context import CodecContext, DecodedChannel
    from .enums import Compression as _Compression, Precision as _Precision

    encoded = CODEC_REGISTRY[_Compression.SEQ_CM].encode(
        DecodedChannel.of_bytes(bytes(np.asarray(run.sequences, dtype=np.uint8).tobytes())),
        CodecContext(read_lengths=np.asarray(run.lengths, dtype=np.uint64)),
    ).dataset_bytes
    arr = np.frombuffer(encoded, dtype=np.uint8)
    ds = sc.create_dataset(
        "sequences",
        _Precision.UINT8,
        length=int(arr.shape[0]),
        chunk_size=io.DEFAULT_SIGNAL_CHUNK,
        compression=_Compression.NONE,
    )
    ds.write(arr)
    io.write_int_attr(ds, "compression", int(_Compression.SEQ_CM), dtype="<u1")


def _write_qualities_fqzcomp_nx16_z(sc, run: WrittenGenomicRun,
                                    qual_strategy_hint: int = -1) -> None:
    """Write the ``qualities`` channel through the FQZCOMP_NX16_Z codec.

    M94.Z is the CRAM-mimic rANS-Nx16 variant — parallel to v1, same
    sibling-channel inputs (read_lengths + revcomp_flags) but a different
    on-wire format (magic ``M94Z`` instead of ``FQZN``). Codec id 12.

    ``qual_strategy_hint``: -1 auto (default, today's 3-way tune), 5/6
    forced V5, 7 V4 with internal preset selection — the streaming
    writers pass their per-run pin here; -1 keeps every existing call
    site byte-identical.
    """
    from .codecs._registry import CODEC_REGISTRY
    from .codecs._context import CodecContext, DecodedChannel
    from .enums import Compression as _Compression, Precision as _Precision

    qualities = bytes(run.qualities.tobytes())
    read_lengths = [int(x) for x in run.lengths]
    revcomp_flags = [
        1 if (int(f) & SAM_REVERSE_FLAG) else 0 for f in run.flags
    ]

    # Codec-registry refactor (Task 5b): route FQZCOMP_NX16_Z through the
    # registry. Its encode adapter requires read_lengths + revcomp_flags
    # in the CodecContext — sourced exactly as the prior direct call did
    # (index lengths + the SAM reverse flag bit). Bytes are identical.
    # Qualities V5 gate (spec 2.4): offer the sequence bytes to the
    # encoder only when the run carries a base-parallel sequences
    # channel and the caller did not opt out; the encoder still emits
    # V4 whenever a V4 strategy wins by exact size.
    v5_sequences = None
    if not getattr(run, "opt_disable_qualities_v5", False):
        seq_arr = np.asarray(run.sequences, dtype=np.uint8)
        if seq_arr.shape[0] == len(qualities):
            v5_sequences = bytes(seq_arr.tobytes())
    encoded = CODEC_REGISTRY[_Compression.FQZCOMP_NX16_Z].encode(
        DecodedChannel.of_bytes(qualities),
        CodecContext(
            read_lengths=np.asarray(read_lengths),
            revcomp_flags=np.asarray(revcomp_flags),
            sequences=v5_sequences,
            qual_strategy_hint=qual_strategy_hint,
        ),
    ).dataset_bytes

    arr = np.frombuffer(encoded, dtype=np.uint8)
    ds = sc.create_dataset(
        "qualities",
        _Precision.UINT8,
        length=int(arr.shape[0]),
        chunk_size=io.DEFAULT_SIGNAL_CHUNK,
        compression=_Compression.NONE,
    )
    ds.write(arr)
    io.write_int_attr(
        ds, "compression",
        int(_Compression.FQZCOMP_NX16_Z), dtype="<u1",
    )


def _write_genomic_run(parent, name: str, run: WrittenGenomicRun,
                       qual_strategy_hint: int = -1) -> None:
    """Write one /study/genomic_runs/<name>/ subtree.

    Mirrors :func:`_write_run` but for the genomic data model. Uses the
    M82 signal-channel helpers from ``_hdf5_io`` and the existing
    compound-dataset writer for variable-length per-read fields.

    ``parent`` may be either a raw h5py group (HDF5 fast path) or a
    :class:`~ttio.providers.base.StorageGroup` (provider path). Raw
    h5py groups are wrapped in :class:`~ttio.providers.hdf5._Group` so
    the signal-channel helpers (which expect the StorageGroup API) work
    identically on both paths.
    """
    from ttio.genomic_index import GenomicIndex
    from ttio.providers.hdf5 import _Group as _H5Group

    # Normalise: the signal-channel helpers call StorageGroup.create_dataset
    # with the (name, Precision, length=N, ...) signature, which differs
    # from h5py's positional API.  Wrap any bare h5py group so both paths
    # use the same StorageGroup interface.
    if isinstance(parent, h5py.Group):
        parent = _H5Group(parent)

    # validate any per-channel codec overrides before we touch
    # the file. The override surface covers the four byte/string
    # channels (sequences, qualities, read_names, cigars). Anything
    # outside the per-channel whitelist is a caller error and must
    # surface immediately ().
    from .enums import Compression as _Compression
    _ALLOWED_OVERRIDE_CODECS_BY_CHANNEL = {
        "sequences": frozenset({
            _Compression.RANS_ORDER0,
            _Compression.RANS_ORDER1,
            _Compression.BASE_PACK,
            # M103: context-aware (the read lengths), written by
            # _write_sequences_seq_cm.
            _Compression.SEQ_CM,
        }),
        "qualities": frozenset({
            _Compression.RANS_ORDER0,
            _Compression.RANS_ORDER1,
            _Compression.BASE_PACK,
            _Compression.QUALITY_BINNED,
            # M94.Z v1.2: CRAM-mimic rANS-Nx16 quality codec. Carries its
            # own read_lengths + needs revcomp_flags from run.flags & 16.
            _Compression.FQZCOMP_NX16_Z,
        }),
        # v1.0 reset: read_names is auto-encoded with NAME_TOKENIZED_V2
        # (codec id 15); no explicit override is supported.
        "read_names": frozenset(),
        # cigars accepts the rANS pair on a length-prefix-
        # concat byte stream of the CIGAR strings (varint(len) + bytes
        # per CIGAR). BASE_PACK and QUALITY_BINNED are wrong-content
        # (CIGARs are not ACGT bytes nor Phred values) and are
        # explicitly rejected with named error messages below.
        "cigars": frozenset({
            _Compression.RANS_ORDER0,
            _Compression.RANS_ORDER1,
        }),
    }
    # per-record integer metadata channels removed from the
    # signal_channels/ override surface. They live exclusively in
    # genomic_index/ now (see comment above). Hard-error so callers
    # with stale code learn immediately.
    _DROPPED_INT_CHANNELS = frozenset(
        {"positions", "flags", "mapping_qualities"}
    )
    for ch_name, codec in run.signal_codec_overrides.items():
        if ch_name in _DROPPED_INT_CHANNELS:
            raise ValueError(
                f"signal_codec_overrides[{ch_name!r}]: removed in v1.6 — "
                f"per-record integer metadata fields ({sorted(_DROPPED_INT_CHANNELS)!r}) "
                f"are stored only under genomic_index/, not "
                f"signal_channels/. The override no longer applies. "
                f"See docs/format-spec.md §4 and §10.7."
            )
        # per-field mate_info_* overrides are disallowed —
        # the inline_v2 codec encodes all three fields together.
        if ch_name in ("mate_info_chrom", "mate_info_pos", "mate_info_tlen"):
            raise ValueError(
                f"signal_codec_overrides[{ch_name!r}]: per-field "
                "mate_info_* overrides are disallowed — the v1.7+ "
                "inline_v2 codec encodes all three mate fields into "
                "a single blob with no per-field codec choice."
            )
        # M86 Phase F / Gotcha §143: the bare
        # "mate_info" key is reserved and rejected with a message
        # pointing at the three per-field names. Without the
        # explicit reject, the bare key would fall through to the
        # generic "channel not supported" branch and the caller
        # would not learn about the per-field surface.
        if ch_name == "mate_info":
            raise ValueError(
                "signal_codec_overrides['mate_info']: the bare "
                "'mate_info' key is reserved and rejected — "
                "mate_info is decomposed at the per-field level in "
                "M86 Phase F. Use one or more of the three per-"
                "field virtual channel names instead: "
                "'mate_info_chrom', 'mate_info_pos', "
                "'mate_info_tlen'. See docs/format-spec.md §10.9."
            )
        if ch_name not in _ALLOWED_OVERRIDE_CODECS_BY_CHANNEL:
            raise ValueError(
                f"signal_codec_overrides: channel '{ch_name}' not supported "
                f"(only sequences, qualities, read_names, and cigars "
                f"can use TTIO codecs)"
            )
        try:
            codec_enum = _Compression(codec)
        except ValueError as exc:
            raise ValueError(
                f"signal_codec_overrides['{ch_name}']: codec {codec!r} "
                "is not a valid Compression value"
            ) from exc
        allowed = _ALLOWED_OVERRIDE_CODECS_BY_CHANNEL[ch_name]
        if codec_enum not in allowed:
            # Phase D : explicit message for the
            # (sequences, QUALITY_BINNED) category error — naming the
            # codec, the channel, and the lossy-quantisation rationale.
            if (
                codec_enum == _Compression.QUALITY_BINNED
                and ch_name == "sequences"
            ):
                raise ValueError(
                    f"signal_codec_overrides['{ch_name}']: codec "
                    f"QUALITY_BINNED is not valid on the '{ch_name}' "
                    "channel — quality binning is lossy and only "
                    "applies to Phred quality scores. Applying it to "
                    "ACGT sequence bytes would silently destroy the "
                    "sequence via Phred-bin quantisation. Use the "
                    "'qualities' channel for QUALITY_BINNED, or "
                    "RANS_ORDER0/RANS_ORDER1/BASE_PACK/SEQ_CM on sequences."
                )
            # Phase C Binding Decisions §120, §121: explicit messages
            # for the wrong-content codecs on the cigars channel. The
            # cigars channel holds variable-length ASCII CIGAR strings
            # — neither ACGT bytes (BASE_PACK) nor Phred quality
            # values (QUALITY_BINNED) match. The error names the
            # codec, the channel, and the wrong-content rationale.
            if ch_name == "cigars":
                if codec_enum == _Compression.BASE_PACK:
                    raise ValueError(
                        f"signal_codec_overrides['{ch_name}']: codec "
                        f"BASE_PACK is not valid on the '{ch_name}' "
                        "channel — BASE_PACK 2-bit-packs ACGT sequence "
                        "bytes and would silently corrupt the structured "
                        "ASCII strings stored on this channel. Use "
                        f"RANS_ORDER0 or RANS_ORDER1 on '{ch_name}'."
                    )
                if codec_enum == _Compression.QUALITY_BINNED:
                    raise ValueError(
                        f"signal_codec_overrides['{ch_name}']: codec "
                        f"QUALITY_BINNED is not valid on the '{ch_name}' "
                        "channel — QUALITY_BINNED quantises Phred "
                        "quality scores onto an 8-bin centre table and "
                        "would silently destroy the structured ASCII "
                        "strings stored on this channel. Use "
                        f"RANS_ORDER0 or RANS_ORDER1 on '{ch_name}'."
                    )
            raise ValueError(
                f"signal_codec_overrides['{ch_name}']: codec {codec!r} "
                f"not supported on the '{ch_name}' channel "
                f"(allowed: {sorted(c.name for c in allowed)})"
            )
        # M97: QUALITY_BINNED is a fixed Illumina-8 bin table; reject
        # it on long-read platforms rather than write a mistuned lossy
        # file.
        from .codecs import quality as _quality
        if (
            ch_name == "qualities"
            and codec_enum == _Compression.QUALITY_BINNED
            and not _quality.binned_allowed_for_platform(run.platform)
        ):
            raise ValueError(
                "signal_codec_overrides['qualities']: QUALITY_BINNED "
                "applies the fixed Illumina-8 bin table, which does "
                "not fit the quality distribution of platform "
                f"{run.platform!r} (M97)."
            )

    rg = parent.create_group(name)

    # Run-level attributes (mirrors _write_run pattern).
    io.write_int_attr(rg, "acquisition_mode", run.acquisition_mode)
    io.write_fixed_string_attr(rg, "modality", "genomic_sequencing")
    io.write_int_attr(rg, "spectrum_class", 5)
    io.write_fixed_string_attr(rg, "reference_uri", run.reference_uri)
    io.write_fixed_string_attr(rg, "platform", run.platform)
    if run.read_role:
        io.write_fixed_string_attr(rg, "read_role", run.read_role)
    io.write_fixed_string_attr(rg, "sample_name", run.sample_name)
    io.write_int_attr(rg, "read_count", int(run.offsets.shape[0]))

    # Genomic index (parallel arrays, including chromosomes as compound).
    idx = GenomicIndex(
        offsets=run.offsets,
        lengths=run.lengths,
        chromosomes=run.chromosomes,
        positions=run.positions,
        mapping_qualities=run.mapping_qualities,
        flags=run.flags,
    )
    idx_group = rg.create_group("genomic_index")
    idx.write(idx_group, name_to_id=run.chrom_name_to_id)

    # Signal channels — these honour run.signal_compression by default;
    # M86 lets per-channel overrides route sequences/qualities through
    # the rANS / BASE_PACK codecs instead.
    sc = rg.create_group("signal_channels")

    # M93 v1.2: REF_DIFF is a context-aware codec — encoding requires
    # positions, cigars, and the reference sequence in addition to the
    # raw byte stream. Dispatch on a special branch when the override
    # selects it; everything else falls through to the existing helper.
    _seq_codec = run.signal_codec_overrides.get("sequences")
    # v1.0 reset (Phase 2c): the v1 REF_DIFF (codec id 9) writer was
    # removed. The default codec lookup now resolves to REF_DIFF_V2
    # (codec id 14) when caller has not selected a per-channel codec
    # AND signal_compression is the "auto-pick best lossless" gzip
    # default AND a reference is available. The
    # _write_sequences_ref_diff_v2 helper handles BASE_PACK fallback
    # (Q5b=C) when the v2 native lib is unavailable / single-chrom
    # check fails / unmapped reads are present.
    if (
        _seq_codec is None
        and run.signal_compression == "gzip"
        and run.reference_chrom_seqs is not None
    ):
        from .genomic._default_codecs import default_codec_for
        _default = default_codec_for("sequences")
        if _default is not None:
            _seq_codec = _default

    # M94.Z v1.2: FQZCOMP_NX16_Z is a v1.5 quality codec — carries
    # read_lengths + revcomp_flags inside the codec wire format. Apply
    # auto-default (Q5a=B): when signal_compression="gzip" AND empty
    # qualities override AND the run is ALREADY a v1.5 candidate (i.e.
    # at least one v1.5 codec is active on this run, whether through an
    # explicit override or the REF_DIFF_V2 auto-default we just resolved
    # for sequences), use FQZCOMP_NX16_Z.
    #
    # The "v1.5 candidate" gate preserves byte-parity for pure-M82
    # baseline writes (no reference, no v1.5 overrides) — those keep
    # the legacy uncompressed/zlib qualities path so existing M82+
    # fixtures remain stable.
    _qual_codec = run.signal_codec_overrides.get("qualities")
    _is_v1_5_candidate = False
    if (
        _qual_codec is None
        and run.signal_compression == "gzip"
    ):
        # Detect v1.5 candidacy: any explicit override is a v1.5 codec,
        # or the sequences channel is going through REF_DIFF_V2 (resolved
        # above into _seq_codec).
        if (_seq_codec is not None
                and _is_valid_compression(_seq_codec)
                and _Compression(_seq_codec) == _Compression.REF_DIFF_V2):
            _is_v1_5_candidate = True
        else:
            for _ovr in run.signal_codec_overrides.values():
                if _is_valid_compression(_ovr):
                    _ce = _Compression(_ovr)
                    if _ce in (
                        _Compression.REF_DIFF_V2,
                        _Compression.FQZCOMP_NX16_Z,
                        _Compression.DELTA_RANS_ORDER0,
                    ):
                        _is_v1_5_candidate = True
                        break
        if _is_v1_5_candidate:
            from .genomic._default_codecs import default_codec_for
            _default = default_codec_for("qualities")
            if _default is not None:
                _qual_codec = _default

    # positions / flags / mapping_qualities are NOT written
    # under signal_channels/. They live exclusively in genomic_index/,
    # mirroring MS's spectrum_index/ pattern (per-record metadata =
    # index; signal_channels = bulk data). See docs/format-spec.md
    # §4 and §10.7. Override-validation rejects these channel names.
    # Old files (v1.5 and earlier) may carry these under
    # signal_channels/ — readers ignore them; the genomic_index/ copy
    # is canonical.
    if (
        run.bulk_v2_blobs is not None
        and run.bulk_v2_blobs.ref_diff_blob is not None
    ):
        # Phase 2c-T verbatim path: skip codec encode and write the
        # wire blob bytes directly.
        if run.bulk_v2_blobs.ref_diff_reference_uri != run.reference_uri:
            raise ValueError(
                f"BulkV2Blobs.ref_diff_reference_uri "
                f"{run.bulk_v2_blobs.ref_diff_reference_uri!r} != "
                f"run.reference_uri {run.reference_uri!r}"
            )
        _write_sequences_ref_diff_bulk_verbatim(
            sc, run.bulk_v2_blobs.ref_diff_blob,
        )
    elif (
        _seq_codec is not None
        and _is_valid_compression(_seq_codec)
        and _Compression(_seq_codec) == _Compression.REF_DIFF_V2
    ):
        _write_sequences_ref_diff_v2(sc, run)
    elif (
        _seq_codec is not None
        and _is_valid_compression(_seq_codec)
        and _Compression(_seq_codec) == _Compression.SEQ_CM
    ):
        _write_sequences_seq_cm(sc, run)
    else:
        io._write_byte_channel_with_codec(
            sc, "sequences", run.sequences, run.signal_compression,
            _seq_codec,
        )
    if (
        _qual_codec is not None
        and _is_valid_compression(_qual_codec)
        and _Compression(_qual_codec) == _Compression.FQZCOMP_NX16_Z
    ):
        _write_qualities_fqzcomp_nx16_z(sc, run, qual_strategy_hint)
    else:
        io._write_byte_channel_with_codec(
            sc, "qualities", run.qualities, run.signal_compression,
            _qual_codec,
        )
    # Variable-length per-read string fields — cigars and read_names are
    # 7-bit ASCII; vl_str() (ASCII encoding) matches the ObjC reader.
    # schema lift for cigars. When an override is
    # present the writer replaces the M82 compound dataset with a
    # flat 1-D uint8 dataset of the same name carrying the codec
    # output, plus an @compression attribute (Binding Decisions
    # §120-§122). The rANS pair operate on a length-prefix-concat
    # byte stream over the CIGAR list (varint(len) + bytes per
    # CIGAR).
    #
    # v1.0 reset (Phase 2c): the v1 NAME_TOKENIZED (codec id 8) writer
    # branch was removed.
    if "cigars" in run.signal_codec_overrides:
        from .enums import Precision as _Precision
        cigars_codec = _Compression(
            run.signal_codec_overrides["cigars"]
        )
        if cigars_codec in (
            _Compression.RANS_ORDER0,
            _Compression.RANS_ORDER1,
        ):
            from .codecs.rans import encode as _rans_enc
            from .codecs._varint import varint_encode as _ve
            # Validate ASCII early so non-ASCII surfaces a clear
            # error before we touch the file (§2.5 contract).
            buf = bytearray()
            for idx, cig in enumerate(run.cigars):
                try:
                    payload = cig.encode("ascii")
                except UnicodeEncodeError as exc:
                    raise ValueError(
                        f"signal_codec_overrides['cigars']: cigar "
                        f"at index {idx} contains non-ASCII bytes "
                        "— CIGARs must be 7-bit ASCII per the SAM "
                        "spec"
                    ) from exc
                buf.extend(_ve(len(payload)))
                buf.extend(payload)
            order = (
                0
                if cigars_codec == _Compression.RANS_ORDER0
                else 1
            )
            encoded = _rans_enc(bytes(buf), order=order)
        else:  # pragma: no cover — validation above rejects this
            raise ValueError(
                f"signal_codec_overrides['cigars']: codec "
                f"{cigars_codec!r} is not supported"
            )
        arr = np.frombuffer(encoded, dtype=np.uint8)
        ds = sc.create_dataset(
            "cigars",
            _Precision.UINT8,
            length=int(arr.shape[0]),
            chunk_size=io.DEFAULT_SIGNAL_CHUNK,
            compression=_Compression.NONE,
        )
        ds.write(arr)
        io.write_int_attr(
            ds, "compression",
            int(cigars_codec), dtype="<u1",
        )
    else:
        io.write_compound_dataset(
            sc,
            "cigars",
            [{"value": c} for c in run.cigars],
            [("value", io.vl_str())],
        )
    # v1.0 reset (Phase 2c) — read_names is always written via the
    # NAME_TOKENIZED_V2 codec (codec id 15) when there is at least one
    # read. The v1 NAME_TOKENIZED writer (codec id 8) and the M82
    # compound fallback were both removed; per-channel overrides for
    # read_names are no longer accepted (validation rejects them
    # earlier in this function).
    #
    # Empty-run case: the writer short-circuits and writes a placeholder
    # NAME_TOKENIZED_V2-tagged empty dataset so readers can detect the
    # layout uniformly. When there are reads but the native library is
    # unavailable, raise a clear RuntimeError pointing at the install
    # path.
    from .codecs import name_tokenizer_v2 as _nt2
    from .enums import Precision as _Precision
    if (
        run.bulk_v2_blobs is not None
        and run.bulk_v2_blobs.name_tok_blob is not None
    ):
        # Phase 2c-T verbatim path.
        _write_read_names_bulk_verbatim(sc, run.bulk_v2_blobs.name_tok_blob)
    elif len(run.read_names) == 0:
        # Short-circuit: write a zero-length uint8 dataset tagged with
        # the v2 codec id so readers dispatch to the v2 path uniformly.
        ds = sc.create_dataset(
            "read_names",
            _Precision.UINT8,
            length=0,
            chunk_size=io.DEFAULT_SIGNAL_CHUNK,
            compression=_Compression.NONE,
        )
        io.write_int_attr(
            ds, "compression",
            int(_Compression.NAME_TOKENIZED_V2), dtype="<u1",
        )
    else:
        if not _nt2.HAVE_NATIVE_LIB:
            raise RuntimeError(
                "NAME_TOKENIZED_V2 codec requires the native libttio_rans "
                "library. Install via 'pip install ttio[native]' or build "
                "from source with --with-native (set TTIO_RANS_LIB_PATH if "
                "the library is at a non-default location)."
            )
        # Codec-registry refactor (Task 5b): NAME_TOKENIZED_V2 needs no
        # run context — encode through the registry with an empty
        # CodecContext. Bytes are identical to the prior direct call.
        from .codecs._registry import CODEC_REGISTRY
        from .codecs._context import CodecContext, DecodedChannel
        encoded = CODEC_REGISTRY[_Compression.NAME_TOKENIZED_V2].encode(
            DecodedChannel.of_str_list(list(run.read_names)),
            CodecContext.empty(),
        ).dataset_bytes
        arr = np.frombuffer(encoded, dtype=np.uint8)
        ds = sc.create_dataset(
            "read_names",
            _Precision.UINT8,
            length=int(arr.shape[0]),
            chunk_size=io.DEFAULT_SIGNAL_CHUNK,
            compression=_Compression.NONE,
        )
        ds.write(arr)
        io.write_int_attr(
            ds, "compression",
            int(_Compression.NAME_TOKENIZED_V2), dtype="<u1",
        )
    # v1.0 reset (Phase 2c) — mate_info is always written via the
    # inline_v2 codec (codec id 13) under
    # signal_channels/mate_info/inline_v2. The v1 per-field subgroup
    # writer (Phase F) and the M82 compound fallback were both
    # removed; per-field mate_info_* overrides are rejected earlier
    # in this function.
    #
    # Empty-run case: if there are no reads, no mate_info group is
    # emitted (the reader treats absence as "no mates"). When there
    # are reads but the native library is unavailable, raise a clear
    # RuntimeError pointing at the install path.
    from .codecs import mate_info_v2 as _miv2
    if (
        run.bulk_v2_blobs is not None
        and run.bulk_v2_blobs.mate_info_blob is not None
    ):
        # Phase 2c-T verbatim path: skip codec encode entirely and write
        # the wire blob bytes + chrom_names table.
        _write_mate_info_bulk_verbatim(
            sc,
            run.bulk_v2_blobs.mate_info_blob,
            run.bulk_v2_blobs.mate_info_chrom_names or [],
        )
    elif len(run.mate_chromosomes) > 0:
        if not _miv2.HAVE_NATIVE_LIB:
            raise RuntimeError(
                "MATE_INLINE_V2 codec requires the native libttio_rans "
                "library. Install via 'pip install ttio[native]' or build "
                "from source with --with-native (set TTIO_RANS_LIB_PATH if "
                "the library is at a non-default location)."
            )
        _write_mate_info_inline_v2(sc, run)

    # M101 — SAM optional fields. A run whose reads carry no tags
    # writes no tags channel, so tag-less files are unchanged.
    if run.tags is not None and any(run.tags):
        _write_tags(sc, run)

    # Per-run provenance — same pattern as _write_run.
    if run.provenance_records:
        prov = rg.create_group("provenance")
        _write_provenance(prov, run.provenance_records, dataset_name="steps")


def _write_tags(sc, run: WrittenGenomicRun) -> None:
    """Write ``signal_channels/tags`` through SAM_TAGS (codec id 18).

    MD:Z / NM:i are recomputed from the reference only when the
    sequences channel went through REF_DIFF_V2, so a reader that can
    decode the sequences can always decode the tags (binding
    decision §96). REF_DIFF_V2 runs hold one chromosome: reads on it
    get the reference, every other read (unmapped '*') none.
    """
    from .codecs._context import CodecContext, DecodedChannel
    from .codecs._registry import CODEC_REGISTRY
    from .enums import Compression as _Compression, Precision as _Precision

    n = len(run.cigars)
    if len(run.tags) != n:
        raise ValueError(f"WrittenGenomicRun.tags has {len(run.tags)} entries for {n} reads")
    try:
        derive = sc.open_group("sequences").has_child("refdiff_v2")
    except Exception:
        derive = False
    refs: list = []
    chrom_ids = np.full(n, 0xFFFF, dtype=np.uint16)
    if derive and run.reference_chrom_seqs:
        names = sorted(set(run.chromosomes) - {"*", ""})
        if len(names) == 1 and run.reference_chrom_seqs.get(names[0]) is not None:
            refs = [run.reference_chrom_seqs[names[0]]]
            chrom_ids[np.asarray(run.chromosomes, dtype=object) == names[0]] = 0
    offsets = np.zeros(n + 1, dtype=np.uint64)
    if n:
        np.cumsum(np.asarray(run.lengths, dtype=np.uint64), out=offsets[1:])
    ctx = CodecContext(
        sequences=bytes(np.asarray(run.sequences, dtype=np.uint8).tobytes()) if refs else None,
        offsets=offsets,
        cigar_strings=list(run.cigars),
        positions=np.asarray(run.positions, dtype=np.int64),
        own_chrom_ids=chrom_ids,
        tag_references=refs,
    )
    blob = CODEC_REGISTRY[_Compression.SAM_TAGS].encode(
        DecodedChannel.of_str_list(list(run.tags)), ctx).dataset_bytes
    arr = np.frombuffer(blob, dtype=np.uint8)
    ds = sc.create_dataset("tags", _Precision.UINT8, length=int(arr.shape[0]),
                           chunk_size=io.DEFAULT_SIGNAL_CHUNK, compression=_Compression.NONE)
    ds.write(arr)
    io.write_int_attr(ds, "compression", int(_Compression.SAM_TAGS), dtype="<u1")


def _build_chrom_id_table(chromosomes: list[str]) -> "tuple[np.ndarray, dict[str, int]]":
    """Encounter-order chrom_id assignment matching the L1 contract.

    Returns (uint16 array of chrom_ids per record, dict name -> id).
    Uses 0xFFFF for unmapped records ('*' or empty string).
    """
    name_to_id: dict[str, int] = {}
    ids = np.empty(len(chromosomes), dtype=np.uint16)
    for i, name in enumerate(chromosomes):
        if name == "*" or not name:
            ids[i] = 0xFFFF
            continue
        if name not in name_to_id:
            name_to_id[name] = len(name_to_id)
        ids[i] = name_to_id[name]
    return ids, name_to_id


def _chrom_ids_with_map(chromosomes: list[str],
                        name_to_id: "dict[str, int]") -> "np.ndarray":
    """Encounter-order ids from (and extending) a shared map."""
    ids = np.empty(len(chromosomes), dtype=np.uint16)
    for i, name in enumerate(chromosomes):
        if name == "*" or not name:
            ids[i] = 0xFFFF
            continue
        if name not in name_to_id:
            if len(name_to_id) > 65535:
                raise ValueError(
                    "genomic_index: > 65,535 unique chromosome names; "
                    "uint16 chromosome_ids would overflow.")
            name_to_id[name] = len(name_to_id)
        ids[i] = name_to_id[name]
    return ids


def _resolve_mate_chrom_ids(
    mate_chromosomes: list[str],
    own_chrom_ids: "np.ndarray",
    name_to_id: "dict[str, int]",
    *,
    extend_in_place: bool = False,
) -> "np.ndarray":
    """Map mate chromosome names to int32 ids; -1 for '*'.

    Uses the same encounter-order dict as own_chrom_ids; extends the
    dict if a mate references a chrom that never appears as own
    (rare cross-chrom case). The '=' SAM shortcut is resolved to the
    record's own chrom_id. name_to_id is copied and not mutated unless
    ``extend_in_place`` is set (shared map across blocks).
    """
    n = len(mate_chromosomes)
    out = np.empty(n, dtype=np.int32)
    local_map = name_to_id if extend_in_place else dict(name_to_id)
    for i, name in enumerate(mate_chromosomes):
        if name == "*" or not name:
            out[i] = -1
        elif name == "=":
            own = own_chrom_ids[i]
            out[i] = -1 if own == 0xFFFF else int(own)
        else:
            if name not in local_map:
                local_map[name] = len(local_map)
            out[i] = local_map[name]
    return out


def _write_bulk_v2_blob(
    parent, *, dataset_name: str, blob: bytes, codec_id: int
) -> None:
    """Write a verbatim v2 codec blob as a 1-D uint8 dataset.

    Used by the transport bulk-mode receiver (Phase 2c-T) to bypass
    the v2 codec encode step and write the wire bytes directly. The
    ``@compression`` attribute is stamped with ``codec_id`` so the
    file remains self-describing — readers dispatch on the same
    compression byte they would for a freshly-encoded blob.
    """
    from .enums import Compression as _Compression, Precision as _Precision

    arr = np.frombuffer(bytes(blob), dtype=np.uint8)
    ds = parent.create_dataset(
        dataset_name,
        _Precision.UINT8,
        length=int(arr.shape[0]),
        chunk_size=io.DEFAULT_SIGNAL_CHUNK,
        compression=_Compression.NONE,
    )
    if arr.shape[0] > 0:
        ds.write(arr)
    io.write_int_attr(ds, "compression", int(codec_id), dtype="<u1")


def _write_mate_info_bulk_verbatim(
    sc, blob: bytes, chrom_names: list[str]
) -> None:
    """Phase 2c-T: write a verbatim ``inline_v2`` blob + chrom_names
    table. Mirrors the layout produced by
    :func:`_write_mate_info_inline_v2` but skips the codec encode."""
    from .enums import Compression as _Compression
    mate_group = sc.create_group("mate_info")
    _write_bulk_v2_blob(
        mate_group,
        dataset_name="inline_v2",
        blob=blob,
        codec_id=int(_Compression.MATE_INLINE_V2),
    )
    io.write_compound_dataset(
        mate_group,
        "chrom_names",
        [{"name": n} for n in chrom_names],
        [("name", io.vl_str())],
    )


def _write_sequences_ref_diff_bulk_verbatim(sc, blob: bytes) -> None:
    """Phase 2c-T: write a verbatim ``refdiff_v2`` blob under
    ``signal_channels/sequences/refdiff_v2``. Mirrors the v1.8 v2
    layout produced by :func:`_write_sequences_ref_diff_v2`."""
    from .enums import Compression as _Compression
    seq_group = sc.create_group("sequences")
    _write_bulk_v2_blob(
        seq_group,
        dataset_name="refdiff_v2",
        blob=blob,
        codec_id=int(_Compression.REF_DIFF_V2),
    )


def _write_read_names_bulk_verbatim(sc, blob: bytes) -> None:
    """Phase 2c-T: write a verbatim ``name_tok_v2`` blob to
    ``signal_channels/read_names``. Mirrors the v1.8 v2 layout
    produced by the inline NAME_TOKENIZED_V2 path."""
    from .enums import Compression as _Compression
    _write_bulk_v2_blob(
        sc,
        dataset_name="read_names",
        blob=blob,
        codec_id=int(_Compression.NAME_TOKENIZED_V2),
    )


def _write_mate_info_inline_v2(sc, run: "WrittenGenomicRun") -> None:
    """v1.7+ inline_v2 writer per spec §4.

    Encodes the full mate triple via libttio_rans (the cross-language
    byte-exact codec from T11) and writes the result as a single
    uint8 blob at signal_channels/mate_info/inline_v2.

    Also writes signal_channels/mate_info/chrom_names — a compound
    dataset mapping chrom_id (uint16) → name (VL_STRING). This covers
    mate-only chromosomes (e.g. a cross-chromosome mate on a chrom that
    no own-read lands on) which are absent from genomic_index/chromosome_names.
    The reader uses this table to resolve mate_chrom_ids returned by
    ttio_mate_info_v2_decode back to string names.
    """
    from .codecs._registry import CODEC_REGISTRY
    from .codecs._context import CodecContext, DecodedChannel
    from .enums import Compression as _Compression, Precision as _Precision

    if run.chrom_name_to_id is not None:
        # Shared map (blocks_v1): ids are stable across blocks and the
        # map grows in place.
        name_to_id = run.chrom_name_to_id
        own_chrom_ids = _chrom_ids_with_map(run.chromosomes, name_to_id)
        mate_chrom_ids = _resolve_mate_chrom_ids(
            run.mate_chromosomes, own_chrom_ids, name_to_id,
            extend_in_place=True)
        full_name_to_id = name_to_id
    else:
        own_chrom_ids, name_to_id = _build_chrom_id_table(run.chromosomes)
        mate_chrom_ids = _resolve_mate_chrom_ids(
            run.mate_chromosomes, own_chrom_ids, name_to_id)

        # After _resolve_mate_chrom_ids, name_to_id may have been extended
        # for mate-only chroms. Reconstruct the full ordered list from the
        # (possibly extended) local map used by _resolve_mate_chrom_ids.
        # Since _resolve_mate_chrom_ids uses a copy, we rebuild from scratch.
        full_name_to_id = dict(name_to_id)
        for name in run.mate_chromosomes:
            if name and name not in ("*", "=") and name not in full_name_to_id:
                full_name_to_id[name] = len(full_name_to_id)
    chrom_names_in_order = sorted(full_name_to_id.keys(),
                                  key=lambda n: full_name_to_id[n])

    # Codec-registry refactor (Task 5b): route MATE_INLINE_V2 through the
    # registry. Its encode adapter consumes the mate triple via the
    # DecodedChannel.of_mate_info value and own_chrom_ids/own_positions
    # from the CodecContext — sourced exactly as the prior direct call.
    # Bytes are identical.
    own_positions = np.asarray(run.positions, dtype=np.int64)
    encoded = CODEC_REGISTRY[_Compression.MATE_INLINE_V2].encode(
        DecodedChannel.of_mate_info({
            "mate_chrom_ids": mate_chrom_ids,
            "mate_positions": np.asarray(run.mate_positions, dtype=np.int64),
            "template_lengths": np.asarray(run.template_lengths, dtype=np.int32),
        }),
        CodecContext(own_chrom_ids=own_chrom_ids, own_positions=own_positions),
    ).dataset_bytes
    arr = np.frombuffer(encoded, dtype=np.uint8)

    mate_group = sc.create_group("mate_info")
    ds = mate_group.create_dataset(
        "inline_v2",
        _Precision.UINT8,
        length=int(arr.shape[0]),
        chunk_size=io.DEFAULT_SIGNAL_CHUNK,
        compression=_Compression.NONE,
    )
    ds.write(arr)
    io.write_int_attr(ds, "compression",
                      int(_Compression.MATE_INLINE_V2), dtype="<u1")

    # Write the full chrom_id → name lookup table (encounter-order, id = row index).
    io.write_compound_dataset(
        mate_group,
        "chrom_names",
        [{"name": n} for n in chrom_names_in_order],
        [("name", io.vl_str())],
    )
