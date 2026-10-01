"""Order-independent digests over the information TTI-O keeps for a
genomic run (SAM columns 1-11, RNEXT expanded), used to check that an
import or export preserved every read.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import hashlib
import subprocess


def _md5_lines(lines) -> str:
    h = hashlib.md5()
    for l in sorted(lines):
        h.update(l.encode()); h.update(b"\n")
    return h.hexdigest()


def sam11_md5(path, *, with_tags: bool = False) -> str:
    """md5 over sorted SAM columns 1-11 of ``samtools view <path>``,
    RNEXT ``=`` expanded to RNAME. ``with_tags`` (M101) adds the
    optional fields, so the digest covers the whole record."""
    p = subprocess.run(["samtools", "view", str(path)], capture_output=True, text=True, check=True)
    lines = []
    for line in p.stdout.splitlines():
        c = line.split("\t", 11)
        if c[6] == "=":
            c[6] = c[2]
        if not with_tags:
            c = c[:11]
        lines.append("\t".join(c))
    return _md5_lines(lines)


def sam_full_md5(path) -> str:
    """md5 over sorted whole SAM records (tags included)."""
    return sam11_md5(path, with_tags=True)


def genomic_run_sam_full_md5(run) -> str:
    """The whole-record digest computed from a GenomicRun's reads."""
    return genomic_run_sam11_md5(run, with_tags=True)


def genomic_run_sam11_md5(run, *, with_tags: bool = False) -> str:
    """The same digest computed from a GenomicRun's reads."""
    lines = []
    for r in run.iter_reads():
        seq = r.sequence or "*"
        q = bytes(r.qualities)
        if not q:
            qual = "*"
        elif all(b == 0xFF for b in q):
            qual = "*"
        else:
            qual = q.decode("latin-1")
        rnext = r.mate_chromosome or "*"
        cols = [
            r.read_name or "*", str(int(r.flags)), r.chromosome or "*", str(int(r.position)),
            str(int(r.mapping_quality)), r.cigar or "*", rnext, str(int(r.mate_position)),
            str(int(r.template_length)), seq, qual]
        if with_tags and r.tags:
            cols.append(r.tags)
        lines.append("\t".join(cols))
    return _md5_lines(lines)


def fastq_md5(path) -> str:
    """md5 over sorted (name, seq, qual) triples of a FASTQ file (gz or plain)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    triples = []
    with opener(path, "rt") as f:
        while True:
            name = f.readline()
            if not name:
                break
            seq = f.readline().rstrip("\n"); f.readline(); qual = f.readline().rstrip("\n")
            triples.append(name[1:].split()[0] + "\t" + seq + "\t" + qual)
    return _md5_lines(triples)


def genomic_run_fastq_md5(run) -> str:
    lines = []
    for r in run.iter_reads():
        lines.append(r.read_name + "\t" + r.sequence + "\t" + bytes(r.qualities).decode("latin-1"))
    return _md5_lines(lines)
