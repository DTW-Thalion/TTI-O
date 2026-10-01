"""SAM_TAGS codec (codec id 18, M101) through the ctypes wrapper.

Round trips per-read tag text with and without MD/NM derivation, and
checks the blob is smaller when MD/NM are recomputed from the reference.
"""
from __future__ import annotations

import numpy as np
import pytest

from ttio.codecs import sam_tags

if not sam_tags.HAVE_NATIVE_LIB:
    pytest.skip(
        "requires native libttio_rans.so via TTIO_RANS_LIB_PATH",
        allow_module_level=True,
    )

REF = b"ACGTACGTACGTTTGGCCAANNNNACGTACGTACGTACGTACGTTTTTTT"


def _ctx(seqs, cigars, positions, chrom_ids, refs):
    off = np.zeros(len(seqs) + 1, dtype=np.uint64)
    np.cumsum([len(s) for s in seqs], out=off[1:])
    return dict(sequences=b"".join(seqs), seq_offsets=off, cigars=cigars,
                positions=np.asarray(positions, dtype=np.int64),
                chrom_ids=np.asarray(chrom_ids, dtype=np.uint16), references=refs)


def _reads():
    seqs = [b"ACGTACGTAC", b"ACCTACGTAC", b"ACGGTAC", b"ACGT", b"", b"ACGTACGTAC"]
    cigars = ["10M", "10M", "3M3D4M", "4M", "*", "10M"]
    tags = [
        "NM:i:0\tMD:Z:10\tAS:i:30\tUQ:i:30\tRG:Z:rg1",
        "NM:i:1\tMD:Z:2G7\tAS:i:47\tUQ:i:47\tRG:Z:rg1",
        "NM:i:3\tMD:Z:3^TAC4\tXA:Z:chr1,+100,10M,0;",
        "",
        "RG:Z:rg1",
        "MD:Z:9\tNM:i:7\tBC:B:C,1,2,3\tXF:f:0.5\tAS:i:007",
    ]
    return seqs, cigars, tags


def test_round_trip_with_reference():
    seqs, cigars, tags = _reads()
    ctx = _ctx(seqs, cigars, [1, 1, 1, 1, 0, 1], [0, 0, 0, 0, 0xFFFF, 0], [REF])
    blob = sam_tags.encode(tags, **ctx)
    assert blob[:4] == b"STG1"
    assert sam_tags.decode(blob, len(tags), **ctx) == tags
    assert sam_tags.encode(tags, **ctx) == blob  # deterministic


def test_round_trip_without_reference():
    _, _, tags = _reads()
    blob = sam_tags.encode(tags)
    assert sam_tags.decode(blob, len(tags)) == tags


def test_derivation_saves_bytes():
    rng = np.random.default_rng(7)
    seqs, cigars, tags, pos = [], [], [], []
    for _ in range(2000):
        p = int(rng.integers(1, len(REF) - 10))
        s = bytearray(REF[p - 1:p + 9])
        k = int(rng.integers(0, 10))
        s[k] = ord("A") if s[k] != ord("A") else ord("C")
        seqs.append(bytes(s))
        cigars.append("10M")
        pos.append(p)
        tags.append("")
    ctx = _ctx(seqs, cigars, pos, [0] * len(seqs), [REF])
    # Let the codec tell us the true MD/NM: encode text that matches calmd.
    for i, (s, p) in enumerate(zip(seqs, pos)):
        ref = REF[p - 1:p + 9]
        md, run, nm = [], 0, 0
        for a, b in zip(s, ref):
            if a == b and b != ord("N"):
                run += 1
            else:
                md.append(f"{run}{chr(b)}")
                run, nm = 0, nm + 1
        tags[i] = f"NM:i:{nm}\tMD:Z:{''.join(md)}{run}"
    with_ref = sam_tags.encode(tags, **ctx)
    without = sam_tags.encode(tags)
    assert sam_tags.decode(with_ref, len(tags), **ctx) == tags
    assert len(with_ref) * 10 < len(without)


def test_empty_block():
    assert sam_tags.decode(sam_tags.encode([]), 0) == []


def test_rejects_nul_and_corruption():
    with pytest.raises(ValueError):
        sam_tags.encode(["XA:Z:a\0b"])
    blob = sam_tags.encode(["RG:Z:rg1", "RG:Z:rg2"])
    with pytest.raises(ValueError):
        sam_tags.decode(blob[:-1], 2)
    with pytest.raises(ValueError):
        sam_tags.decode(blob, 3)
