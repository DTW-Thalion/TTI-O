"""BAM exporter — M88.

Writes a :class:`~ttio.written_genomic_run.WrittenGenomicRun` to BAM
by formatting the in-memory parallel-array representation as SAM
text and piping that text via stdin to the user-installed
``samtools`` binary (``samtools view -bS -``, optionally piped
through ``samtools sort -O bam``). Subprocess-only — no htslib
linkage; SAM line layout is from the public SAMv1 spec.

Quality byte encoding
---------------------
M87's :class:`~ttio.importers.bam.BamReader` stores SAM's QUAL field
bytes verbatim into ``WrittenGenomicRun.qualities`` — i.e. the
buffer holds **ASCII Phred+33** characters (so a Phred-40 score is
stored as the byte value 73, the ASCII code for ``'I'``). This
writer mirrors that convention: each ``qualities[i]`` byte is
written directly as the SAM QUAL character with no arithmetic
adjustment. The pair is therefore lossless byte-for-byte across the
M87 read → M88 write round trip.

Cross-language note: ObjC and Java implementations must adopt the
same convention (store QUAL bytes verbatim on read; emit them
verbatim on write) so that conformance dumps match.

Cross-language equivalents
--------------------------
Objective-C: ``TTIOBamWriter`` · Java:
``global.thalion.ttio.exporters.BamWriter``.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Iterable

from ..importers.bam import _check_samtools
from ..io.progress import ProgressSinkLike, _fire
from ..provenance import ProvenanceRecord
from ..written_genomic_run import WrittenGenomicRun


__all__ = ["BamWriter", "PROGRESS_INTERVAL_READS"]


#: Mirror Java's ``BamWriter.PROGRESS_INTERVAL_READS``.
PROGRESS_INTERVAL_READS = 1000


# Default @SQ length when the writer doesn't know the true reference
# length. SAM requires LN: on every @SQ; we pick INT32_MAX so the
# emitted header is valid for any plausible coordinate. samtools'
# downstream consumers (IGV, GATK) accept this fallback. The same
# value should be used by the ObjC and Java writers for cross-
# language byte-equality on the unsorted code path.
_DEFAULT_SQ_LENGTH = 2147483647


class BamWriter:
    """Write a :class:`~ttio.written_genomic_run.WrittenGenomicRun` to BAM.

    Parameters
    ----------
    path : str or :class:`pathlib.Path`
        Output BAM file path. The ``.bam`` extension is honoured by
        samtools' file-format auto-detection (Gotcha §165).

    Notes
    -----
    The ``samtools`` binary is a runtime dependency. Construction
    succeeds without samtools on PATH; :meth:`write` raises
    :class:`~ttio.importers.bam.SamtoolsNotFoundError` when samtools
    is missing at first call (from M87).
    """

    def __init__(self, path: str | os.PathLike[str]):
        """Configure the writer with an output BAM path.

        Parameters
        ----------
        path : str or os.PathLike
            Destination BAM file path. The path is stored verbatim;
            the file is not opened or created until :meth:`write` runs.
        """
        self._path = Path(path)

    @property
    def path(self) -> Path:
        """Return the configured output path as a :class:`pathlib.Path`.

        Returns
        -------
        pathlib.Path
            The destination path supplied at construction.
        """
        return self._path

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write(
        self,
        run: WrittenGenomicRun,
        provenance_records: list[ProvenanceRecord] | None = None,
        sort: bool = True,
        *,
        progress: "ProgressSinkLike | None" = None,
    ) -> None:
        """Serialise ``run`` to the configured output path.

        Parameters
        ----------
        run : WrittenGenomicRun
            The genomic-run container to write.
        provenance_records : list[ProvenanceRecord] or None
            Optional provenance records to inject as ``@PG`` header
            lines. If ``None``, falls back to ``run.provenance_records``
            so the most common Python-side call is one-arg. Java and
            ObjC pass this explicitly because their
            ``WrittenGenomicRun`` analogues don't carry provenance.
        sort : bool, default True
            When ``True`` (the default ),
            pipes the SAM text through ``samtools sort -O bam`` so
            the output BAM is coordinate-sorted (the precondition
            most BAM consumers expect — IGV, GATK, ``samtools
            index``). When ``False``, output is written in the input
            ``run``'s read order and the ``@HD SO:`` tag is set to
            ``unsorted``.
        """
        _check_samtools()
        if provenance_records is None:
            provenance_records = list(_provenance_of(run))
        # SAM text is streamed into the samtools pipe line by line, so
        # a run of any size (a GenomicRun walked through iter_reads, one
        # decoded block at a time) exports with bounded memory.
        chunks = self._iter_sam_text(run, provenance_records, sort=sort, progress=progress)
        self._invoke_samtools(chunks, sort=sort)

    # ------------------------------------------------------------------
    # SAM text assembly
    # ------------------------------------------------------------------

    def _build_sam_text(
        self,
        run: WrittenGenomicRun,
        provenance_records: list[ProvenanceRecord],
        *,
        sort: bool,
        progress: ProgressSinkLike | None = None,
    ) -> str:
        """The full SAM text (header + alignment lines) as one string;
        :meth:`write` streams :meth:`_iter_sam_text` instead."""
        return "".join(self._iter_sam_text(run, provenance_records, sort=sort, progress=progress))

    def _iter_sam_text(
        self,
        run,
        provenance_records: list[ProvenanceRecord],
        *,
        sort: bool,
        progress: ProgressSinkLike | None = None,
    ):
        """Yield the header, then one SAM line per read.

        Fires ``progress`` every :data:`PROGRESS_INTERVAL_READS` reads
        with ``total = len(run)`` and once more at the end with
        ``(total, total)``.
        """
        yield self._build_header(run, provenance_records, sort=sort)
        total = _read_count_of(run)
        idx = 0
        for idx, line in enumerate(self._iter_alignment_lines(run), 1):
            yield line
            if idx % PROGRESS_INTERVAL_READS == 0:
                _fire(progress, idx, total)
        _fire(progress, total, total)

    @staticmethod
    def _build_header(
        run: WrittenGenomicRun,
        provenance_records: list[ProvenanceRecord],
        *,
        sort: bool,
    ) -> str:
        """Emit the @HD / @SQ / @RG / @PG header block."""
        lines: list[str] = []

        so = "coordinate" if sort else "unsorted"
        lines.append(f"@HD\tVN:1.6\tSO:{so}")

        # @SQ — one per unique chromosome (excluding "*" which is
        # the SAM unmapped sentinel and not a real reference). Emit
        # in first-seen order so writer output is deterministic.
        seen: set[str] = set()
        for chrom in _chromosome_names_of(run):
            if not chrom or chrom == "*" or chrom in seen:
                continue
            seen.add(chrom)
            lines.append(f"@SQ\tSN:{chrom}\tLN:{_DEFAULT_SQ_LENGTH}")

        # @RG — single line if either sample_name or platform is set.
        if run.sample_name or run.platform:
            rg_parts = ["@RG", "ID:rg1"]
            if run.sample_name:
                rg_parts.append(f"SM:{run.sample_name}")
            if run.platform:
                rg_parts.append(f"PL:{run.platform}")
            lines.append("\t".join(rg_parts))

        # @PG — one line per provenance record. SAM requires ID;
        # synthesize "pg<idx>" if the record's software field is
        # blank or collides.
        used_ids: set[str] = set()
        for idx, prov in enumerate(provenance_records):
            base_id = prov.software or f"pg{idx}"
            pg_id = base_id
            n = 1
            while pg_id in used_ids:
                pg_id = f"{base_id}.{n}"
                n += 1
            used_ids.add(pg_id)
            pg_parts = [
                "@PG",
                f"ID:{pg_id}",
                f"PN:{prov.software}",
            ]
            cl = prov.parameters.get("CL") if prov.parameters else None
            if cl:
                pg_parts.append(f"CL:{cl}")
            lines.append("\t".join(pg_parts))

        return "\n".join(lines) + "\n"

    @staticmethod
    def _iter_alignment_lines(run: WrittenGenomicRun) -> Iterable[str]:
        """Yield one SAM alignment text line per read in ``run``.

        Field handling per :
        - QNAME / RNAME / CIGAR: ``"*"`` sentinel preserved.
        - FLAG / MAPQ / TLEN: decimal ints (signed for TLEN,
          unsigned otherwise).
        - POS / PNEXT: decimal ints; ``mate_position == -1`` is
          mapped to SAM's ``0`` .
        - RNEXT: collapsed to ``=`` when equal to RNAME per Binding
          Decision §136 (writer-side reverse of M87's expansion).
        - SEQ / QUAL: ASCII bytes from the concatenated
          sequences/qualities buffers, sliced by
          ``offsets[i]:offsets[i]+lengths[i]``. Empty slice -> ``"*"``.
        """
        if not isinstance(run, WrittenGenomicRun):
            yield from BamWriter._iter_alignment_lines_lazy(run)
            return
        seq_buf = bytes(run.sequences)
        qual_buf = bytes(run.qualities)
        tags = run.tags

        n = len(run.read_names)
        for i in range(n):
            qname = run.read_names[i] or "*"
            flag = int(run.flags[i])
            rname = run.chromosomes[i] or "*"
            pos = int(run.positions[i])
            mapq = int(run.mapping_qualities[i])
            cigar = run.cigars[i] or "*"

            # RNEXT collapse (§136).
            mate_chrom = run.mate_chromosomes[i] or "*"
            if mate_chrom == rname and rname != "*":
                rnext = "="
            else:
                rnext = mate_chrom

            # PNEXT mapping (§138).
            mate_pos = int(run.mate_positions[i])
            pnext = 0 if mate_pos < 0 else mate_pos

            tlen = int(run.template_lengths[i])

            offset = int(run.offsets[i])
            length = int(run.lengths[i])
            if length == 0:
                seq = "*"
                qual = "*"
            else:
                seq_bytes = seq_buf[offset:offset + length]
                qual_bytes = qual_buf[offset:offset + length]
                seq = seq_bytes.decode("ascii")
                # M87's reader produces an all-0xff buffer when the
                # source SAM had QUAL '*' but a non-empty SEQ. Map
                # that back to SAM's '*' on write so the round trip
                # canonicalises to the source convention.
                if qual_bytes and all(b == 0xff for b in qual_bytes):
                    qual = "*"
                else:
                    # qual stored as ASCII Phred+33 already (see
                    # module docstring). Just decode latin-1 to keep
                    # the bytes round-tripping when any value > 127
                    # ever sneaks through; in practice samtools
                    # rejects QUAL > '~' (0x7e).
                    qual = qual_bytes.decode("latin-1")

            # M101: the optional fields follow column 11 verbatim.
            tail = f"\t{tags[i]}" if tags and tags[i] else ""
            yield (
                f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t"
                f"{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}{tail}\n"
            )

    # ------------------------------------------------------------------
    # samtools subprocess invocation
    # ------------------------------------------------------------------

    def _invoke_samtools(self, sam_text, *, sort: bool) -> None:
        """Pipe SAM text (a string or an iterable of string chunks)
        through samtools to produce the BAM file.

        Subclasses (CramWriter) override this to inject reference and
        format flags.
        """
        cmd_view, cmd_sort = self._build_samtools_commands(sort=sort)

        if cmd_sort is None:
            # Single-stage: samtools view -bS -o <path> -
            self._run_pipeline([cmd_view], sam_text)
        else:
            # Two-stage: view -bS -  | sort -O bam -o <path>
            self._run_pipeline([cmd_view, cmd_sort], sam_text)

    def _build_samtools_commands(
        self, *, sort: bool,
    ) -> tuple[list[str], list[str] | None]:
        """Return (view-cmd, sort-cmd-or-None) for the BAM pipeline.

        Subclasses override to swap in CRAM flags.
        """
        if sort:
            view = ["samtools", "view", "-bS", "-"]
            sort_cmd = ["samtools", "sort", "-O", "bam", "-o",
                        str(self._path), "-"]
            return view, sort_cmd
        else:
            view = ["samtools", "view", "-bS", "-o",
                    str(self._path), "-"]
            return view, None

    @staticmethod
    def _run_pipeline(commands: list[list[str]], stdin_text) -> None:
        """Run a 1- or 2-stage samtools pipeline; raise on non-zero exit.
        ``stdin_text`` is a string or an iterable of string chunks that
        are written to the first stage as they come."""
        chunks = [stdin_text] if isinstance(stdin_text, str) else stdin_text
        if len(commands) == 1:
            proc = subprocess.Popen(
                commands[0],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            try:
                assert proc.stdin is not None
                for c in chunks:
                    proc.stdin.write(c.encode("ascii"))
            finally:
                proc.stdin.close()
            err = proc.stderr.read() if proc.stderr else b""
            proc.wait()
            for h in (proc.stdout, proc.stderr):
                if h is not None:
                    h.close()
            if proc.returncode != 0:
                stderr = (err or b"").decode("utf-8", errors="replace")
                raise RuntimeError(
                    f"samtools exited {proc.returncode}: "
                    f"{stderr.strip()[:500]}"
                )
            return

        # Two-stage pipeline: stage[0].stdout -> stage[1].stdin.
        first = subprocess.Popen(
            commands[0],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        second = subprocess.Popen(
            commands[1],
            stdin=first.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Allow first to receive SIGPIPE if second exits.
        if first.stdout is not None:
            first.stdout.close()

        try:
            assert first.stdin is not None
            try:
                for c in chunks:
                    first.stdin.write(c.encode("ascii"))
            finally:
                first.stdin.close()
            second_out, second_err = second.communicate()
            first.wait()
        except subprocess.TimeoutExpired:
            first.kill()
            second.kill()
            raise

        first_err = first.stderr.read() if first.stderr else b""
        if first.stderr is not None:
            first.stderr.close()

        if first.returncode != 0:
            err = first_err.decode("utf-8", errors="replace")
            raise RuntimeError(
                f"samtools (stage 1, {commands[0][:3]}) exited "
                f"{first.returncode}: {err.strip()[:500]}"
            )
        if second.returncode != 0:
            err = second_err.decode("utf-8", errors="replace")
            raise RuntimeError(
                f"samtools (stage 2, {commands[1][:3]}) exited "
                f"{second.returncode}: {err.strip()[:500]}"
            )


def _provenance_of(run) -> list[ProvenanceRecord]:
    if isinstance(run, WrittenGenomicRun):
        return list(run.provenance_records)
    try:
        return list(run.provenance_chain())
    except Exception:
        return []


def _read_count_of(run) -> int:
    if isinstance(run, WrittenGenomicRun):
        return len(run.read_names)
    return len(run)


def _chromosome_names_of(run):
    """Own chromosome names in first-seen order. For a GenomicRun the
    run-level genomic_index/chromosome_names table is used (no per-read
    array is loaded)."""
    if isinstance(run, WrittenGenomicRun):
        return list(run.chromosomes)
    try:
        from .. import _hdf5_io as io
        rows = io.read_compound_dataset(run.group.open_group("genomic_index"), "chromosome_names")
        return [(r["name"].decode() if isinstance(r["name"], bytes) else r["name"]) for r in rows]
    except Exception:
        return list(run.index.chromosome_names or [])


def _lazy_lines(run):
    for r in run.iter_reads():
        qname = r.read_name or "*"
        rname = r.chromosome or "*"
        mate_chrom = r.mate_chromosome or "*"
        rnext = "=" if (mate_chrom == rname and rname != "*") else mate_chrom
        mate_pos = int(r.mate_position)
        pnext = 0 if mate_pos < 0 else mate_pos
        seq = r.sequence or "*"
        q = bytes(r.qualities)
        if not r.sequence:
            seq, qual = "*", "*"
        elif q and all(b == 0xff for b in q):
            qual = "*"
        else:
            qual = q.decode("latin-1")
        tail = f"\t{r.tags}" if r.tags else ""
        yield (
            f"{qname}\t{int(r.flags)}\t{rname}\t{int(r.position)}\t{int(r.mapping_quality)}\t"
            f"{r.cigar or '*'}\t{rnext}\t{pnext}\t{int(r.template_length)}\t{seq}\t{qual}{tail}\n"
        )


BamWriter._iter_alignment_lines_lazy = staticmethod(_lazy_lines)
