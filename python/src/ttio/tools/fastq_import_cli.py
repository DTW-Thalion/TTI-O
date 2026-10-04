"""FASTQ import CLI — converts a FASTQ file into an unaligned
genomic run inside a ``.tio`` container.

Usage
-----
::

    python -m ttio.tools.fastq_import_cli \\
        --fastq reads.fq.gz --out study.tio --sample S1

Exit codes
----------
- ``0`` — wrote output successfully.
- ``1`` — argument error.
- ``2`` — read / write failure.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..importers.fastq import FastqParseError, FastqReader
from ..spectral_dataset import SpectralDataset


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ttio.tools.fastq_import_cli",
        description="Import FASTQ into a .tio container.",
    )
    p.add_argument("--fastq", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--name", default="genomic_0001")
    p.add_argument("--sample", default="")
    p.add_argument("--platform", default="")
    p.add_argument("--read-role", default=None,
                   help="@read_role run attribute (hifi, ont_ul, hic_r1, "
                        "hic_r2, parental_maternal, parental_paternal, "
                        "illumina_wgs)")
    p.add_argument("--phred", type=int, choices=(33, 64), default=None,
                   help="force Phred offset (default: auto-detect)")
    p.add_argument("--block-reads", type=int, default=None,
                   help="reads per blocks_v1 block (default 1000000)")
    p.add_argument("--block-bytes", type=int, default=None,
                   help="sequence bytes per blocks_v1 block (default 256 MiB)")
    p.add_argument("--legacy-whole-channel", action="store_true",
                   help="write the v1.8 whole-channel layout (memory-unbounded)")
    p.add_argument("--group-reads", action="store_true",
                   help="reorder the reads by sequence before blocking and store the "
                        "permutation (blocks_v1_grouped; memory-unbounded)")
    return p


def main(argv: list[str] | None = None) -> int:
    """Convert a FASTQ file into an unaligned genomic run.

    Writes the resulting run under ``/study/genomic_runs/<name>/``
    inside a freshly created ``.tio`` container. The Phred offset is
    auto-detected unless ``--phred`` is supplied.

    Parameters
    ----------
    argv : list[str], optional
        Argument vector. Defaults to ``sys.argv[1:]`` when ``None``.

    Returns
    -------
    int
        ``0`` on success, ``2`` on FASTQ parse / I/O failure.
        Argparse exits with ``2`` on usage errors.
    """
    args = _parser().parse_args(argv)

    try:
        reader = FastqReader(args.fastq, force_phred=args.phred)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    try:
        src = reader.stream_source(name=args.name, sample_name=args.sample,
                                   platform=args.platform)
        src.block_reads = args.block_reads
        src.block_bytes = args.block_bytes
        src.opt_legacy_whole_channel = args.legacy_whole_channel
        src.group_reads = args.group_reads
        src.read_role = args.read_role
        SpectralDataset.write_minimal(
            args.out,
            title="",
            isa_investigation_id="",
            runs={},
        )
        with SpectralDataset.open(args.out, writable=True) as ds:
            n = src.write_into(ds.study_group)
        print(
            f"wrote unaligned run {args.name!r} "
            f"({n} reads, "
            f"Phred+{reader.detected_phred_offset}) to {args.out}"
        )
        return 0
    except FastqParseError as e:
        print(f"FASTQ parse error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
