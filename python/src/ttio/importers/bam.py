"""SAM/BAM importer — M87.

Wraps the user-installed ``samtools`` binary as a subprocess to read
SAM and BAM (Sequence Alignment/Map) files into
:class:`~ttio.written_genomic_run.WrittenGenomicRun` instances. No
htslib source is linked or consulted; SAM/BAM format parsing is from
the public SAMv1 specification (https://samtools.github.io/hts-specs).

The subprocess approach mirrors :mod:`ttio.importers.thermo_raw`
(M38) and :mod:`ttio.importers.bruker_tdf` (M53). ``samtools`` is a
runtime dependency only — ``import ttio.importers.bam`` succeeds on
systems without samtools; only :meth:`BamReader.to_genomic_run`
requires the binary on PATH ().

samtools auto-detects SAM vs BAM format from magic bytes; one parser
handles both. The companion :class:`~ttio.importers.sam.SamReader`
exists as a discoverable convenience alias.

Cross-language equivalents
--------------------------
Objective-C: ``TTIOBamReader`` · Java:
``global.thalion.ttio.importers.BamReader``.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable

import numpy as np

from ..enums import AcquisitionMode
from ..io.progress import ProgressSinkLike, _fire
from ..provenance import ProvenanceRecord
from ..written_genomic_run import WrittenGenomicRun


__all__ = ["BamReader", "SamtoolsNotFoundError", "PROGRESS_INTERVAL_READS"]


#: Mirror Java's ``BamReader.PROGRESS_INTERVAL_READS``.
PROGRESS_INTERVAL_READS = 1000


_INSTALL_HELP = (
    "samtools is required by ttio.importers.bam but was not found on "
    "PATH. Install it via your platform's package manager:\n"
    "  Debian/Ubuntu: apt install samtools\n"
    "  macOS:         brew install samtools\n"
    "  Conda:         conda install -c bioconda samtools\n"
    "Then re-run."
)


class SamtoolsNotFoundError(RuntimeError):
    """Raised at first use when ``samtools`` is not available on PATH.

    The class is a subclass of :class:`RuntimeError` so callers can
    catch it loosely; the message includes platform-appropriate
    install guidance (apt / brew / conda).
    """


def _samtools_on_path() -> bool:
    """Return True iff ``samtools`` is resolvable via :func:`shutil.which`."""
    return shutil.which("samtools") is not None


def _check_samtools() -> None:
    """Raise :class:`SamtoolsNotFoundError` if samtools is missing.

    Performs the PATH check via :func:`shutil.which` and additionally
    invokes ``samtools --version`` to verify the binary is callable.
    Per this happens at first use, NOT at module
    import time.
    """
    if not _samtools_on_path():
        raise SamtoolsNotFoundError(_INSTALL_HELP)
    try:
        # capture binary; samtools --version prints copyright bytes
        # that aren't always strict UTF-8.
        proc = subprocess.run(
            ["samtools", "--version"],
            capture_output=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise SamtoolsNotFoundError(
            f"{_INSTALL_HELP}\n(invocation failed: {exc})"
        ) from exc
    if proc.returncode != 0:
        stderr_text = (proc.stderr or b"").decode("utf-8", errors="replace")
        stdout_text = (proc.stdout or b"").decode("utf-8", errors="replace")
        raise SamtoolsNotFoundError(
            f"{_INSTALL_HELP}\n"
            f"(samtools --version exited {proc.returncode}: "
            f"{(stderr_text or stdout_text).strip()[:200]})"
        )


class BamReader:
    """Read a SAM or BAM file via the ``samtools view`` subprocess.

    Parameters
    ----------
    path : str or :class:`pathlib.Path`
        Filesystem path to a SAM or BAM file. samtools auto-detects
        the format from magic bytes.

    Notes
    -----
    The ``samtools`` binary is a runtime dependency, not a build
    dependency. Construction succeeds without samtools on PATH;
    :meth:`to_genomic_run` raises :class:`SamtoolsNotFoundError` when
    samtools cannot be located at first use ().
    """

    def __init__(self, path: str | os.PathLike[str]):
        """Configure the reader; does not invoke samtools.

        Parameters
        ----------
        path : str or os.PathLike
            Filesystem path to a SAM or BAM file. The path is stored
            verbatim — existence is not checked until
            :meth:`to_genomic_run` runs.
        """
        self._path = Path(path)

    @property
    def path(self) -> Path:
        """Return the SAM/BAM input path as a :class:`pathlib.Path`.

        Returns
        -------
        pathlib.Path
            The path supplied at construction, unchanged.
        """
        return self._path

    def _view_cmd(self, region: str | None) -> list[str]:
        """The ``samtools view -h`` command for this input; CRAM adds
        ``--reference``."""
        cmd = ["samtools", "view", "-h", str(self._path)]
        if region is not None:
            cmd.append(region)
        return cmd

    def to_genomic_run(
        self,
        name: str = "genomic_0001",
        region: str | None = None,
        sample_name: str | None = None,
        *,
        progress: ProgressSinkLike | None = None,
    ) -> WrittenGenomicRun:
        """Read the BAM/SAM and return a :class:`WrittenGenomicRun`.

        Holds the whole run in memory; for large inputs use
        :meth:`iter_batches` with a
        :class:`~ttio.genomic.GenomicStreamWriter` (what ``ttio encode``
        does).

        Parameters
        ----------
        name : str
            The genomic-run name (becomes the subgroup name under
            ``/study/genomic_runs/<name>/``). Default
            ``"genomic_0001"``.
        region : str or None
            Optional region filter passed verbatim to
            ``samtools view`` (e.g. ``"chr1:1000-2000"`` or ``"*"``
            for unmapped reads).
        sample_name : str or None
            Optional override for the run's ``sample_name``. If
            ``None``, derived from the first ``@RG SM:`` tag in the
            header (or the empty string if no @RG present).

        Raises
        ------
        SamtoolsNotFoundError
            If ``samtools`` is not on PATH at first call.
        FileNotFoundError
            If the input path does not exist.
        RuntimeError
            If ``samtools view`` exits non-zero (stderr included in
            message) or if a SAM line is malformed.
        """
        from ..genomic._blocks import concat_runs
        batches = list(self.iter_batches(region=region, sample_name=sample_name,
                                         batch_reads=1 << 62, progress=progress))
        if not batches:
            return self._empty_run(sample_name)
        return concat_runs(batches)

    def stream_source(self, *, name: str = "genomic_0001", region: str | None = None,
                      sample_name: str | None = None, batch_reads: int = 100_000,
                      progress: ProgressSinkLike | None = None,
                      reference_fasta: str | os.PathLike[str] | None = None,
                      embed_reference: bool = False):
        """A :class:`~ttio.importers.import_result.GenomicStreamSource`
        that feeds :class:`~ttio.genomic.GenomicStreamWriter` batch by
        batch. ``reference_fasta`` enables REF_DIFF_V2 through a
        :class:`~ttio.genomic.lazy_reference.LazyReference`."""
        from .import_result import GenomicStreamSource
        return GenomicStreamSource(
            name=name,
            iter_batches=lambda: self.iter_batches(region=region, sample_name=sample_name,
                                                   batch_reads=batch_reads, progress=progress),
            reference_fasta=Path(reference_fasta) if reference_fasta else None,
            embed_reference=embed_reference,
        )

    def _empty_run(self, sample_name: str | None) -> WrittenGenomicRun:
        z = np.zeros(0, dtype=np.uint8)
        return WrittenGenomicRun(
            acquisition_mode=int(AcquisitionMode.GENOMIC_WGS), reference_uri="",
            platform="", sample_name=sample_name or "",
            positions=np.zeros(0, dtype=np.int64), mapping_qualities=z,
            flags=np.zeros(0, dtype=np.uint32), sequences=z, qualities=z,
            offsets=np.zeros(0, dtype=np.uint64), lengths=np.zeros(0, dtype=np.uint32),
            cigars=[], read_names=[], mate_chromosomes=[],
            mate_positions=np.zeros(0, dtype=np.int64),
            template_lengths=np.zeros(0, dtype=np.int32), chromosomes=[],
        )

    def iter_batches(
        self,
        *,
        region: str | None = None,
        sample_name: str | None = None,
        batch_reads: int = 100_000,
        progress: ProgressSinkLike | None = None,
    ):
        """Yield the input as consecutive :class:`WrittenGenomicRun`
        batches of at most ``batch_reads`` reads, parsing the
        ``samtools view`` pipe line by line. Every batch carries the
        run-level metadata (reference_uri from the first ``@SQ``,
        platform and sample from ``@RG``); the header ``@PG``
        provenance rides on the first batch only.
        """
        _check_samtools()
        if not self._path.exists():
            raise FileNotFoundError(f"BAM/SAM file not found: {self._path}")
        if batch_reads < 1:
            raise ValueError("batch_reads must be >= 1")

        proc = subprocess.Popen(
            self._view_cmd(region),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Header state
        sq_names: list[str] = []
        rg_sample: str = ""
        rg_platform: str = ""
        provenance: list[ProvenanceRecord] = []
        provenance_sent = False

        # Per-read accumulators (reset per batch)
        acc = _BatchAccumulator()
        total_reads = 0

        # Provenance timestamp comes from the file mtime per
        try:
            file_mtime = int(self._path.stat().st_mtime)
        except OSError:
            file_mtime = int(time.time())

        def _emit() -> WrittenGenomicRun:
            nonlocal provenance_sent
            effective_sample = sample_name if sample_name is not None else rg_sample
            reference_uri = sq_names[0] if sq_names else ""
            prov = [] if provenance_sent else list(provenance)
            provenance_sent = True
            return acc.to_run(
                acquisition_mode=int(AcquisitionMode.GENOMIC_WGS),
                reference_uri=reference_uri, platform=rg_platform,
                sample_name=effective_sample, provenance=prov)

        try:
            assert proc.stdout is not None
            for line_no, raw_line in enumerate(proc.stdout, 1):
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("@"):
                    self._parse_header_line(
                        line,
                        sq_names=sq_names,
                        provenance=provenance,
                        rg_state=[rg_sample, rg_platform],
                        file_mtime=file_mtime,
                    )
                    # rg_state mutation happens via local function; use
                    # a different style to avoid the closure trick:
                    if line.startswith("@RG") and not rg_sample:
                        sm, pl = self._parse_rg_fields(line)
                        if sm and not rg_sample:
                            rg_sample = sm
                        if pl and not rg_platform:
                            rg_platform = pl
                    continue

                # Alignment record: fields 1-11, then the optional tags
                # (M101) kept as their tab-joined text in cols[11].
                cols = line.split("\t", 11)
                if len(cols) < 11:
                    raise RuntimeError(
                        f"Malformed SAM alignment at line {line_no}: "
                        f"expected >=11 tab-separated fields, got {len(cols)}"
                        f" — {line[:120]}"
                    )
                qname, flag_s, rname, pos_s, mapq_s, cigar, \
                    rnext, pnext_s, tlen_s, seq, qual = cols[:11]

                try:
                    flag = int(flag_s)
                    pos = int(pos_s)
                    mapq = int(mapq_s)
                    pnext = int(pnext_s)
                    tlen = int(tlen_s)
                except ValueError as exc:
                    raise RuntimeError(
                        f"Malformed SAM numeric field at line {line_no}: "
                        f"{exc} — {line[:120]}"
                    ) from exc

                # RNEXT special handling — :
                # "=" expands to RNAME so downstream consumers don't
                # need to remember the convention.
                if rnext == "=":
                    rnext = rname

                # SEQ / QUAL: "*" means absent — contributes 0 bytes.
                # Per Gotcha §153 cigars[i] keeps "*" literally; SEQ/
                # QUAL are reduced to empty bytes in the buffer (the
                # offsets/lengths pair carries the "absent" signal).
                if seq == "*":
                    seq_bytes = b""
                else:
                    seq_bytes = seq.encode("ascii")
                if qual == "*":
                    qual_bytes = b"" if seq == "*" else b"\xff" * len(seq_bytes)
                else:
                    qual_bytes = qual.encode("ascii")

                # SAM spec: SEQ and QUAL must be the same length when
                # both present. We don't try to "fix" inputs; we
                # truncate qual to seq length if mismatched (samtools
                # already validated on the wire side).
                if len(qual_bytes) != len(seq_bytes):
                    if seq == "*":
                        # SEQ absent but qual present — discard qual.
                        qual_bytes = b""
                    elif qual == "*":
                        # Already handled above (filled to seq length).
                        pass
                    else:
                        raise RuntimeError(
                            f"SEQ/QUAL length mismatch at line {line_no}: "
                            f"SEQ={len(seq_bytes)} QUAL={len(qual_bytes)}"
                        )

                acc.add(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen,
                        seq_bytes, qual_bytes, cols[11] if len(cols) > 11 else "")
                total_reads += 1

                if total_reads % PROGRESS_INTERVAL_READS == 0:
                    # samtools subprocess doesn't pre-count, so total
                    # stays -1 until the final fire below.
                    _fire(progress, total_reads, -1)
                if acc.n >= batch_reads:
                    yield _emit()
                    acc = _BatchAccumulator()

            proc.wait()
            if proc.returncode != 0:
                stderr_text = (proc.stderr.read()
                               if proc.stderr else "") or ""
                raise RuntimeError(
                    f"samtools view exited {proc.returncode} for "
                    f"{self._path}: {stderr_text.strip()[:500]}"
                )
        finally:
            try:
                if proc.stdout is not None:
                    proc.stdout.close()
            except Exception:
                pass
            try:
                if proc.stderr is not None:
                    proc.stderr.close()
            except Exception:
                pass
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()

        # Final progress fire: total is now known.
        _fire(progress, total_reads, total_reads)
        if acc.n > 0 or total_reads == 0:
            yield _emit()

    # ------------------------------------------------------------------
    # Header-line parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_header_fields(line: str) -> dict[str, str]:
        """Split a SAM header line into a {KEY: VALUE} dict.

        Skips the leading ``@TAG`` token. Tolerates fields that lack
        a colon by silently dropping them (samtools never emits
        these but be defensive).
        """
        fields: dict[str, str] = {}
        for token in line.split("\t")[1:]:
            if ":" not in token:
                continue
            k, _, v = token.partition(":")
            fields[k] = v
        return fields

    @classmethod
    def _parse_rg_fields(cls, line: str) -> tuple[str, str]:
        """Return (sample, platform) from an ``@RG`` header line."""
        fields = cls._parse_header_fields(line)
        return fields.get("SM", ""), fields.get("PL", "")

    @classmethod
    def _parse_header_line(
        cls,
        line: str,
        *,
        sq_names: list[str],
        provenance: list[ProvenanceRecord],
        rg_state: list[str],   # unused; kept for legacy signature
        file_mtime: int,
    ) -> None:
        """Dispatch a header line to the appropriate accumulator.

        Only @SQ and @PG are accumulated into structured state here;
        @RG is handled inline in :meth:`to_genomic_run` so the
        first-wins rule () is obvious at the
        callsite. @HD and @CO are read but not mapped to TTI-O
        fields in v0
        """
        if line.startswith("@SQ"):
            fields = cls._parse_header_fields(line)
            sn = fields.get("SN")
            if sn:
                sq_names.append(sn)
        elif line.startswith("@PG"):
            fields = cls._parse_header_fields(line)
            program = fields.get("PN", "")
            command_line = fields.get("CL", "")
            params: dict[str, object] = {}
            if command_line:
                params["CL"] = command_line
            for k in ("ID", "VN", "PP"):
                if k in fields:
                    params[k] = fields[k]
            provenance.append(
                ProvenanceRecord(
                    timestamp_unix=file_mtime,
                    software=program,
                    parameters=params,
                )
            )
        # @HD, @CO, @RG: handled elsewhere or ignored in v0.


class _BatchAccumulator:
    """Per-read accumulators for one batch of SAM records."""

    __slots__ = ("read_names", "chromosomes", "positions", "mapqs", "flags", "cigars",
                 "mate_chromosomes", "mate_positions", "template_lengths", "offsets",
                 "lengths", "seq_chunks", "qual_chunks", "tags", "running", "n")

    def __init__(self) -> None:
        self.read_names: list[str] = []
        self.chromosomes: list[str] = []
        self.positions: list[int] = []
        self.mapqs: list[int] = []
        self.flags: list[int] = []
        self.cigars: list[str] = []
        self.mate_chromosomes: list[str] = []
        self.mate_positions: list[int] = []
        self.template_lengths: list[int] = []
        self.offsets: list[int] = []
        self.lengths: list[int] = []
        self.seq_chunks: list[bytes] = []
        self.qual_chunks: list[bytes] = []
        self.tags: list[str] = []
        self.running = 0
        self.n = 0

    def add(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen,
            seq_bytes, qual_bytes, tags: str = "") -> None:
        self.tags.append(tags)
        self.read_names.append(qname)
        self.flags.append(flag)
        self.chromosomes.append(rname)
        self.positions.append(pos)
        self.mapqs.append(mapq)
        self.cigars.append(cigar)
        self.mate_chromosomes.append(rnext)
        self.mate_positions.append(pnext)
        self.template_lengths.append(tlen)
        length = len(seq_bytes)
        self.offsets.append(self.running)
        self.lengths.append(length)
        self.seq_chunks.append(seq_bytes)
        self.qual_chunks.append(qual_bytes)
        self.running += length
        self.n += 1

    def to_run(self, *, acquisition_mode, reference_uri, platform, sample_name,
               provenance) -> WrittenGenomicRun:
        return WrittenGenomicRun(
            acquisition_mode=acquisition_mode,
            reference_uri=reference_uri,
            platform=platform,
            sample_name=sample_name,
            positions=np.asarray(self.positions, dtype=np.int64),
            mapping_qualities=np.asarray(self.mapqs, dtype=np.uint8),
            flags=np.asarray(self.flags, dtype=np.uint32),
            sequences=np.frombuffer(b"".join(self.seq_chunks), dtype=np.uint8).copy(),
            qualities=np.frombuffer(b"".join(self.qual_chunks), dtype=np.uint8).copy(),
            offsets=np.asarray(self.offsets, dtype=np.uint64),
            lengths=np.asarray(self.lengths, dtype=np.uint32),
            cigars=self.cigars,
            read_names=self.read_names,
            mate_chromosomes=self.mate_chromosomes,
            mate_positions=np.asarray(self.mate_positions, dtype=np.int64),
            template_lengths=np.asarray(self.template_lengths, dtype=np.int32),
            chromosomes=self.chromosomes,
            tags=self.tags,
            provenance_records=list(provenance),
        )
