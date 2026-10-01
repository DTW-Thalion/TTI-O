"""bam_dump — canonical-JSON dump of a SAM/BAM/CRAM file for the
M87 / M88.1 cross-language conformance harness.

Usage::

    # BAM/SAM (M87):
    python -m ttio.importers.bam_dump <path>

    # CRAM (M88.1):
    python -m ttio.importers.bam_dump <path.cram> --reference <fa>

Reads the file via :class:`~ttio.importers.bam.BamReader` for SAM/BAM
or :class:`~ttio.importers.cram.CramReader` for CRAM (auto-dispatched
on the `.cram` extension) and emits a canonical JSON document on
stdout matching the schema documented in HANDOFF.md M87 §7. The
same shape is produced by the ObjC ``TtioBamDump`` and Java
``BamDump`` CLIs.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from typing import Any

from .bam import BamReader
from .cram import CramReader


__all__ = ["dump", "main"]


def dump(
    path: str,
    name: str = "genomic_0001",
    reference: str | None = None,
) -> dict[str, Any]:
    """Read ``path`` and return the canonical-JSON-shaped dict.

    If ``path`` ends in ``.cram`` (case-insensitive), a
    :class:`CramReader` is used; ``reference`` must then be provided.
    Otherwise a :class:`BamReader` handles the file (samtools
    auto-detects SAM vs BAM); ``reference`` is accepted but unused.

    Returned keys are unchanged from M87 — see this module's
    docstring for the full schema.
    """
    if path.lower().endswith(".cram"):
        if reference is None:
            raise ValueError(
                "--reference <fasta> is required for .cram input"
            )
        reader = CramReader(path, reference)
    else:
        reader = BamReader(path)

    run = reader.to_genomic_run(name=name)

    seq_md5 = hashlib.md5(bytes(run.sequences)).hexdigest()
    qual_md5 = hashlib.md5(bytes(run.qualities)).hexdigest()

    return {
        "name": name,
        "read_count": len(run.read_names),
        "sample_name": run.sample_name,
        "platform": run.platform,
        "reference_uri": run.reference_uri,
        "read_names": list(run.read_names),
        "positions": [int(x) for x in run.positions],
        "chromosomes": list(run.chromosomes),
        "flags": [int(x) for x in run.flags],
        "mapping_qualities": [int(x) for x in run.mapping_qualities],
        "cigars": list(run.cigars),
        "mate_chromosomes": list(run.mate_chromosomes),
        "mate_positions": [int(x) for x in run.mate_positions],
        "template_lengths": [int(x) for x in run.template_lengths],
        "sequences_md5": seq_md5,
        "qualities_md5": qual_md5,
        "provenance_count": len(run.provenance_records),
        # M101: the SAM optional fields of each read, as samtools
        # prints them ("" when a read has none).
        "tags": list(run.tags) if run.tags is not None else [""] * len(run.read_names),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m ttio.importers.bam_dump",
        description=(
            "Emit canonical M87/M88.1 JSON for a SAM, BAM, or CRAM file."
        ),
    )
    parser.add_argument(
        "path",
        help="Path to a SAM, BAM, or CRAM file (.cram dispatches to CramReader).",
    )
    parser.add_argument(
        "--reference", default=None,
        help="Path to reference FASTA — required for .cram input, ignored otherwise.",
    )
    parser.add_argument(
        "--name", default="genomic_0001",
        help="Genomic-run name to embed in the JSON (default: genomic_0001).",
    )
    args = parser.parse_args(argv)

    try:
        payload = dump(args.path, name=args.name, reference=args.reference)
    except ValueError as exc:
        parser.error(str(exc))  # exits 2 with message on stderr

    sys.stdout.write(json.dumps(payload, sort_keys=True, indent=2))
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
