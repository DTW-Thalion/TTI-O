#!/usr/bin/env python3
"""M101 Phase 0: a column codec for SAM optional tags, proved and sized.

Reads alignment records as SAM text (``samtools view``), cuts them into
blocks the way the blocks_v1 writer does, and runs every block through a
prototype of the SAM_TAGS column codec:

* each read's tags become one entry of a per-block tag-line dictionary
  (the ordered ``KEY:TYPE`` sequence) plus one value per tag;
* values go to one column per (key, type, kind); integer columns try a
  zigzag-varint and a fixed-width byte-plane transform, text columns try
  NUL-joined bytes and NAME_TOKENIZED_V2, and each column keeps the
  smallest of raw, rANS order-0 and rANS order-1;
* MD:Z and NM:i that equal the values recomputed from the read, its
  CIGAR and the reference are marked ``derived`` in the tag line and
  store nothing.

Every block is decoded again and compared with the input tag text, so a
run that finishes is a byte-exact round trip. The report gives bytes per
column against the plain alternative (the tab-joined tag text of each
read, length-prefixed, through rANS order-1) and the input text size.

Usage (WSL, from the repo root, ttio installed and TTIO_RANS_LIB_PATH set):

    python tools/prototypes/m101_sam_tags/tags_codec_proof.py \
        data/genomic/hg002_2x250/hg002_2x250.chr22.bam \
        --reference data/genomic/hg002_2x250/GRCh38_analysis_set.chr22.fa
"""

from __future__ import annotations

import argparse
import collections
import io
import json
import multiprocessing as mp
import struct
import subprocess
import sys
import time

import numpy as np

from ttio.codecs import name_tokenizer_v2, rans

MAGIC = b"STG0"  # prototype framing; the spec assigns the real one

# Column kinds. A tag-line entry is (key, sam_type, kind).
KIND_INT = 0      # canonical decimal integer (type i)
KIND_TEXT = 1     # verbatim value text after "KEY:T:"
KIND_DERIVED = 2  # MD:Z / NM:i equal to the recomputed value; nothing stored
KIND_DUP = 3      # integer equal to entry <arg> earlier in the same read; nothing stored

# Substream methods.
M_RAW, M_RANS0, M_RANS1, M_NAMETOK = 0, 1, 2, 3
# Integer transforms.
T_VARINT, T_PLANES, T_CONST = 0, 1, 2
# Text transforms.
TT_JOINED, TT_NAMETOK, TT_CONST = 0, 1, 2


# --------------------------------------------------------------------------
# varints

def put_uvarint(out: bytearray, v: int) -> None:
    while v >= 0x80:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)


def get_uvarint(buf: memoryview, pos: int) -> tuple[int, int]:
    shift = v = 0
    while True:
        b = buf[pos]
        pos += 1
        v |= (b & 0x7F) << shift
        if b < 0x80:
            return v, pos
        shift += 7


def zigzag(v: int) -> int:
    return (v << 1) ^ (v >> 63)


def unzigzag(u: int) -> int:
    return (u >> 1) ^ -(u & 1)


# --------------------------------------------------------------------------
# substreams: pick the smallest of raw / rANS-O0 / rANS-O1

def pack_bytes(data: bytes, allow=(M_RAW, M_RANS0, M_RANS1)) -> bytes:
    best = None
    for m in allow:
        if m == M_RAW:
            body = data
        elif m == M_RANS0:
            body = rans.encode(data, 0)
        else:
            body = rans.encode(data, 1)
        if best is None or len(body) < len(best[1]):
            best = (m, body)
    out = bytearray([best[0]])
    put_uvarint(out, len(best[1]))
    return bytes(out) + best[1]


def unpack_bytes(buf: memoryview, pos: int) -> tuple[bytes, int, int]:
    m = buf[pos]
    n, pos = get_uvarint(buf, pos + 1)
    body = bytes(buf[pos:pos + n])
    pos += n
    if m == M_RAW:
        return body, pos, m
    return rans.decode(body), pos, m


# --------------------------------------------------------------------------
# columns

def encode_int_column(vals: list[int]) -> bytes:
    if all(v == vals[0] for v in vals):
        out = bytearray([T_CONST])
        put_uvarint(out, zigzag(vals[0]))
        return bytes(out)
    cands = []
    vb = bytearray()
    for v in vals:
        put_uvarint(vb, zigzag(v))
    cands.append(bytes([T_VARINT]) + pack_bytes(bytes(vb)))
    # Byte planes, each its own substream with its own method: the low
    # plane carries the information, the high planes are mostly zero.
    zz = np.fromiter((zigzag(v) for v in vals), dtype=np.uint64, count=len(vals))
    mx = int(zz.max())
    w = 1 if mx < 1 << 8 else 2 if mx < 1 << 16 else 4 if mx < 1 << 32 else 8
    planes = zz.astype(f"<u{w}").view(np.uint8).reshape(-1, w).T
    body = bytearray([T_PLANES, w])
    for k in range(w):
        body += pack_bytes(planes[k].tobytes())
    cands.append(bytes(body))
    return min(cands, key=len)


def decode_int_column(buf: memoryview, pos: int, n: int) -> tuple[list[int], int]:
    t = buf[pos]
    pos += 1
    if t == T_CONST:
        u, pos = get_uvarint(buf, pos)
        return [unzigzag(u)] * n, pos
    if t == T_VARINT:
        raw, pos, _ = unpack_bytes(buf, pos)
        mv, p, out = memoryview(raw), 0, []
        for _ in range(n):
            u, p = get_uvarint(mv, p)
            out.append(unzigzag(u))
        return out, pos
    w = buf[pos]
    pos += 1
    planes = []
    for _ in range(w):
        raw, pos, _ = unpack_bytes(buf, pos)
        planes.append(np.frombuffer(raw, dtype=np.uint8))
    zz = np.stack(planes).T.copy().view(f"<u{w}").ravel()
    return [unzigzag(int(u)) for u in zz], pos


def encode_text_column(vals: list[str]) -> bytes:
    if all(v == vals[0] for v in vals):
        b = vals[0].encode("ascii")
        out = bytearray([TT_CONST])
        put_uvarint(out, len(b))
        return bytes(out) + b
    joined = "\0".join(vals).encode("ascii")
    cands = [bytes([TT_JOINED]) + pack_bytes(joined)]
    if all(v and "\0" not in v for v in vals):
        nt = name_tokenizer_v2.encode(vals)
        hdr = bytearray([TT_NAMETOK])
        put_uvarint(hdr, len(nt))
        cands.append(bytes(hdr) + nt)
    return min(cands, key=len)


def decode_text_column(buf: memoryview, pos: int, n: int) -> tuple[list[str], int]:
    t = buf[pos]
    pos += 1
    if t == TT_CONST:
        ln, pos = get_uvarint(buf, pos)
        return [bytes(buf[pos:pos + ln]).decode("ascii")] * n, pos + ln
    if t == TT_JOINED:
        raw, pos, _ = unpack_bytes(buf, pos)
        return raw.decode("ascii").split("\0"), pos
    ln, pos = get_uvarint(buf, pos)
    vals = name_tokenizer_v2.decode(bytes(buf[pos:pos + ln]))
    return vals, pos + ln


# --------------------------------------------------------------------------
# MD / NM recomputation (samtools calmd rules: an N on either side is a
# mismatch; '=' in SEQ is a match; soft clips and insertions are skipped
# for MD; NM = mismatches + inserted + deleted bases)

def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    ops, n = [], 0
    for ch in cigar:
        if "0" <= ch <= "9":
            n = n * 10 + ord(ch) - 48
        else:
            ops.append((n, ch))
            n = 0
    return ops


def calc_md_nm(seq: bytes, cigar: str, pos1: int, ref: bytes) -> tuple[str, int] | None:
    """Return (MD, NM) for one read, or None if not computable."""
    if cigar == "*" or seq == b"*":
        return None
    rpos = pos1 - 1
    qpos = 0
    md: list[str] = []
    run = 0
    nm = 0
    for n, op in parse_cigar(cigar):
        if op in "M=X":
            if rpos + n > len(ref):
                return None
            q = seq[qpos:qpos + n]
            r = ref[rpos:rpos + n]
            for i in range(n):
                qb = q[i]
                rb = r[i]
                if qb == 61:  # '='
                    run += 1
                    continue
                if (qb | 0x20) == (rb | 0x20) and (rb | 0x20) != 0x6E and (qb | 0x20) != 0x6E:
                    run += 1
                else:
                    md.append(str(run))
                    md.append(chr(rb).upper())
                    run = 0
                    nm += 1
            qpos += n
            rpos += n
        elif op == "I":
            qpos += n
            nm += n
        elif op == "S":
            qpos += n
        elif op == "D":
            if rpos + n > len(ref):
                return None
            md.append(str(run))
            md.append("^" + ref[rpos:rpos + n].decode("ascii").upper())
            run = 0
            rpos += n
            nm += n
        elif op == "N":
            rpos += n
        elif op in "HP":
            pass
        else:
            return None
    md.append(str(run))
    return "".join(md), nm


# --------------------------------------------------------------------------
# block codec

def split_tags(text: str) -> list[tuple[str, str, str]]:
    if not text:
        return []
    out = []
    for f in text.split("\t"):
        # KEY:T:VALUE ; VALUE may contain ':'
        out.append((f[0:2], f[3], f[5:]))
    return out


def encode_block(tag_texts, seqs, cigars, pos1s, ref, derive: bool):
    """Return (blob, per-column byte sizes, derivation stats)."""
    line_ids: dict[tuple, int] = {}
    lines: list[tuple] = []
    tl: list[int] = []
    cols: dict[tuple, list] = collections.OrderedDict()
    stats = collections.Counter()
    for i, text in enumerate(tag_texts):
        entries = split_tags(text)
        mdnm = None
        line = []
        ints: list[int | None] = []  # integer value of each entry so far, for DUP
        for key, typ, val in entries:
            kind = KIND_TEXT
            arg = 0
            iv = None
            if typ == "i":
                try:
                    iv = int(val)
                    if str(iv) == val:
                        kind = KIND_INT
                    else:
                        iv = None
                except ValueError:
                    pass
            if derive and ((key == "MD" and typ == "Z") or (key == "NM" and typ == "i")):
                if mdnm is None:
                    mdnm = calc_md_nm(seqs[i], cigars[i], pos1s[i], ref) or ("", -1)
                want = mdnm[0] if key == "MD" else str(mdnm[1])
                if val == want:
                    kind = KIND_DERIVED
                    stats[key + "_derived"] += 1
                else:
                    stats[key + "_verbatim"] += 1
            if kind == KIND_INT and iv in ints:
                kind = KIND_DUP
                arg = ints.index(iv)
                stats[key + "_dup"] += 1
            ints.append(iv if kind in (KIND_INT, KIND_DUP) or (kind == KIND_DERIVED and key == "NM") else None)
            line.append((key, typ, kind, arg))
            if kind == KIND_INT:
                cols.setdefault((key, typ, kind), []).append(int(val))
            elif kind == KIND_TEXT:
                cols.setdefault((key, typ, kind), []).append(val)
        lt = tuple(line)
        li = line_ids.get(lt)
        if li is None:
            li = line_ids[lt] = len(lines)
            lines.append(lt)
        tl.append(li)

    out = bytearray(MAGIC)
    put_uvarint(out, len(tag_texts))
    put_uvarint(out, len(lines))
    for lt in lines:
        put_uvarint(out, len(lt))
        for key, typ, kind, arg in lt:
            out += key.encode("ascii") + typ.encode("ascii") + bytes([kind, arg])
    sizes = collections.OrderedDict()
    n0 = len(out)
    tlb = bytearray()
    for v in tl:
        put_uvarint(tlb, v)
    out += pack_bytes(bytes(tlb))
    sizes["(tag-line ids)"] = len(out) - n0
    sizes["(dictionary+header)"] = n0
    put_uvarint(out, len(cols))
    for (key, typ, kind), vals in cols.items():
        n0 = len(out)
        out += key.encode("ascii") + typ.encode("ascii") + bytes([kind])
        put_uvarint(out, len(vals))
        if kind == KIND_INT:
            out += encode_int_column(vals)
        else:
            out += encode_text_column(vals)
        sizes[f"{key}:{typ}"] = sizes.get(f"{key}:{typ}", 0) + len(out) - n0
    return bytes(out), sizes, stats


def decode_block(blob: bytes, seqs, cigars, pos1s, ref) -> list[str]:
    buf = memoryview(blob)
    assert bytes(buf[:4]) == MAGIC
    pos = 4
    n, pos = get_uvarint(buf, pos)
    nl, pos = get_uvarint(buf, pos)
    lines = []
    for _ in range(nl):
        k, pos = get_uvarint(buf, pos)
        lt = []
        for _ in range(k):
            lt.append((bytes(buf[pos:pos + 2]).decode(), chr(buf[pos + 2]), buf[pos + 3], buf[pos + 4]))
            pos += 5
        lines.append(lt)
    tlraw, pos, _ = unpack_bytes(buf, pos)
    tmv, p, tl = memoryview(tlraw), 0, []
    for _ in range(n):
        v, p = get_uvarint(tmv, p)
        tl.append(v)
    ncols, pos = get_uvarint(buf, pos)
    cols = {}
    for _ in range(ncols):
        key = bytes(buf[pos:pos + 2]).decode()
        typ = chr(buf[pos + 2])
        kind = buf[pos + 3]
        pos += 4
        cnt, pos = get_uvarint(buf, pos)
        if kind == KIND_INT:
            vals, pos = decode_int_column(buf, pos, cnt)
        else:
            vals, pos = decode_text_column(buf, pos, cnt)
        cols[(key, typ, kind)] = collections.deque(vals)
    out = []
    for i in range(n):
        fields = []
        vals: list[str] = []
        mdnm = None
        for key, typ, kind, arg in lines[tl[i]]:
            if kind == KIND_DERIVED:
                if mdnm is None:
                    mdnm = calc_md_nm(seqs[i], cigars[i], pos1s[i], ref)
                val = mdnm[0] if key == "MD" else str(mdnm[1])
            elif kind == KIND_DUP:
                val = vals[arg]
            else:
                val = str(cols[(key, typ, kind)].popleft())
            vals.append(val)
            fields.append(f"{key}:{typ}:{val}")
        out.append("\t".join(fields))
    return out


def plain_blob(tag_texts) -> bytes:
    raw = bytearray()
    for t in tag_texts:
        b = t.encode("ascii")
        put_uvarint(raw, len(b))
        raw += b
    return rans.encode(bytes(raw), 1)


# --------------------------------------------------------------------------
# driver

_REF: bytes = b""


def _init_worker(ref_path: str | None):
    global _REF
    _REF = load_ref(ref_path) if ref_path else b""


def _run_block(args):
    idx, tag_texts, seqs, cigars, pos1s, derive = args
    t0 = time.perf_counter()
    blob, sizes, stats = encode_block(tag_texts, seqs, cigars, pos1s, _REF, derive)
    t1 = time.perf_counter()
    back = decode_block(blob, seqs, cigars, pos1s, _REF)
    t2 = time.perf_counter()
    if back != list(tag_texts):
        bad = next(i for i, (a, b) in enumerate(zip(back, tag_texts)) if a != b)
        raise SystemExit(f"block {idx}: round trip differs at read {bad}:\n  in  {tag_texts[bad]!r}\n  out {back[bad]!r}")
    plain = plain_blob(tag_texts)
    return idx, len(tag_texts), len(blob), len(plain), sizes, stats, t1 - t0, t2 - t1


def load_ref(path: str) -> bytes:
    seq = io.BytesIO()
    with open(path, "rb") as f:
        for line in f:
            if not line.startswith(b">"):
                seq.write(line.strip())
    return seq.getvalue()


def iter_blocks(bam: str, region: str | None, block_reads: int, block_bytes: int, max_reads: int | None):
    cmd = ["samtools", "view", "-@4", bam] + ([region] if region else [])
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    tags, seqs, cigars, pos1s = [], [], [], []
    nbytes = total = 0
    try:
        for line in proc.stdout:
            cols = line.rstrip("\n").split("\t", 11)
            tags.append(cols[11] if len(cols) > 11 else "")
            seqs.append(cols[9].encode("ascii"))
            cigars.append(cols[5])
            pos1s.append(int(cols[3]))
            nbytes += 2 * len(cols[9])  # sequences + qualities, as the block cut counts
            total += 1
            if len(tags) >= block_reads or nbytes >= block_bytes:
                yield tags, seqs, cigars, pos1s
                tags, seqs, cigars, pos1s, nbytes = [], [], [], [], 0
            if max_reads and total >= max_reads:
                break
        if tags:
            yield tags, seqs, cigars, pos1s
    finally:
        proc.kill()
        proc.wait()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("--reference", help="FASTA of the single chromosome in the slice (enables MD/NM derivation)")
    ap.add_argument("--region")
    ap.add_argument("--block-reads", type=int, default=1_000_000)
    ap.add_argument("--block-bytes", type=int, default=64 << 20)
    ap.add_argument("--max-reads", type=int)
    ap.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 2))
    ap.add_argument("--json", help="write the report as JSON here")
    a = ap.parse_args()

    derive = a.reference is not None
    text_bytes = 0
    reads = blob_total = plain_total = 0
    sizes = collections.Counter()
    stats = collections.Counter()
    enc_s = dec_s = 0.0
    nblocks = 0
    t0 = time.time()

    def jobs():
        nonlocal text_bytes
        for i, (tg, sq, cg, ps) in enumerate(iter_blocks(a.bam, a.region, a.block_reads, a.block_bytes, a.max_reads)):
            text_bytes += sum(len(t) + 1 for t in tg)
            yield i, tg, sq, cg, ps, derive

    with mp.Pool(a.workers, initializer=_init_worker, initargs=(a.reference,)) as pool:
        for idx, n, nb, npl, sz, st, te, td in pool.imap_unordered(_run_block, jobs()):
            nblocks += 1
            reads += n
            blob_total += nb
            plain_total += npl
            sizes.update(sz)
            stats.update(st)
            enc_s += te
            dec_s += td
            print(f"  block {idx}: {n} reads, column {nb} B, plain {npl} B", file=sys.stderr)

    wall = time.time() - t0
    report = {
        "bam": a.bam,
        "reads": reads,
        "blocks": nblocks,
        "derive_md_nm": derive,
        "input_tag_text_bytes": text_bytes,
        "plain_rans1_bytes": plain_total,
        "column_codec_bytes": blob_total,
        "column_vs_plain": round(blob_total / plain_total, 4) if plain_total else None,
        "column_bytes_per_read": round(blob_total / reads, 3) if reads else None,
        "per_column_bytes": dict(sorted(sizes.items(), key=lambda kv: -kv[1])),
        "derivation": dict(stats),
        "encode_cpu_s": round(enc_s, 1),
        "decode_cpu_s": round(dec_s, 1),
        "wall_s": round(wall, 1),
        "round_trip": "byte-exact",
    }
    print(json.dumps(report, indent=2))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
