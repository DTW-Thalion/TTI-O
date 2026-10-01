"""Python ctypes wrapper for SAM_TAGS (codec id 18, M101).

SAM optional fields (columns 12+, tab-joined as samtools prints them)
coded as a tag-line dictionary plus one column per tag key, with MD:Z
and NM:i recomputed from the reference when they match. The blob is
written as ``signal_channels/tags`` with ``@compression = 18``. Spec:
``docs/codecs/sam_tags.md``.

Encode and decode take the same context: the reads' sequences, CIGARs,
1-based positions and chromosome ids, and one reference sequence per
chromosome id (``None`` where there is none). With no references, MD
and NM are stored like any other tag. Requires TTIO_RANS_LIB_PATH to
point at a built libttio_rans.
"""
from __future__ import annotations

import ctypes
from collections.abc import Sequence

import numpy as np

from .fqzcomp_nx16_z import _HAVE_NATIVE_LIB, _native_lib

HAVE_NATIVE_LIB: bool = _HAVE_NATIVE_LIB

ERR_PARAM = -1
ERR_ALLOC = -2
ERR_CORRUPT = -3

_ERR_MESSAGES = {
    ERR_PARAM: "invalid parameters (a NUL in the tag text, or bad offsets)",
    ERR_ALLOC: "out of memory in native code",
    ERR_CORRUPT: "corrupt SAM_TAGS blob",
}

_u8p = ctypes.POINTER(ctypes.c_uint8)
_u64p = ctypes.POINTER(ctypes.c_uint64)


class _Ctx(ctypes.Structure):
    _fields_ = [
        ("n_reads", ctypes.c_uint64),
        ("sequences", _u8p),
        ("seq_offsets", _u64p),
        ("cigars", _u8p),
        ("cigar_offsets", _u64p),
        ("positions", ctypes.POINTER(ctypes.c_int64)),
        ("chrom_ids", ctypes.POINTER(ctypes.c_uint16)),
        ("refs", ctypes.POINTER(_u8p)),
        ("ref_lengths", _u64p),
        ("n_refs", ctypes.c_uint32),
    ]


if HAVE_NATIVE_LIB:
    _lib = _native_lib
    _lib.ttio_sam_tags_encode.argtypes = [
        ctypes.POINTER(_Ctx), _u8p, _u64p,
        ctypes.POINTER(_u8p), ctypes.POINTER(ctypes.c_size_t),
    ]
    _lib.ttio_sam_tags_encode.restype = ctypes.c_int
    _lib.ttio_sam_tags_decode.argtypes = [
        ctypes.POINTER(_Ctx), _u8p, ctypes.c_size_t,
        ctypes.POINTER(_u8p), _u64p,
    ]
    _lib.ttio_sam_tags_decode.restype = ctypes.c_int
    _lib.ttio_sam_tags_free.argtypes = [ctypes.c_void_p]
    _lib.ttio_sam_tags_free.restype = None


def _require_lib(what: str) -> None:
    if not HAVE_NATIVE_LIB:
        raise RuntimeError(f"sam_tags.{what} requires libttio_rans (set TTIO_RANS_LIB_PATH)")


def _joined(values: Sequence[str] | Sequence[bytes]) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate str/bytes values into (uint8 buffer, uint64 offsets[n+1])."""
    enc = [v.encode("utf-8") if isinstance(v, str) else bytes(v) for v in values]
    off = np.zeros(len(enc) + 1, dtype=np.uint64)
    if enc:
        np.cumsum([len(b) for b in enc], out=off[1:])
    buf = np.frombuffer(b"".join(enc) or b"\0", dtype=np.uint8)
    return buf, off


class _Context:
    """Owns the arrays a ttio_sam_tags_ctx points into."""

    def __init__(self, n_reads: int, sequences, seq_offsets, cigars, positions,
                 chrom_ids, references) -> None:
        self.ctx = _Ctx()
        self.ctx.n_reads = n_reads
        refs = list(references or [])
        if not any(r is not None for r in refs):
            self.ctx.n_refs = 0
            return
        if n_reads and (sequences is None or cigars is None or positions is None or chrom_ids is None):
            raise ValueError("MD/NM derivation needs sequences, cigars, positions and chrom_ids")
        self._seq = np.ascontiguousarray(np.frombuffer(bytes(sequences), dtype=np.uint8)
                                         if isinstance(sequences, (bytes, bytearray, memoryview))
                                         else sequences, dtype=np.uint8)
        if self._seq.size == 0:
            self._seq = np.zeros(1, dtype=np.uint8)
        self._so = np.ascontiguousarray(seq_offsets, dtype=np.uint64)
        self._cig, self._co = _joined(cigars)
        self._pos = np.ascontiguousarray(positions, dtype=np.int64)
        self._chr = np.ascontiguousarray(chrom_ids, dtype=np.uint16)
        for name, arr, want in (("seq_offsets", self._so, n_reads + 1),
                                ("cigars", self._co, n_reads + 1),
                                ("positions", self._pos, n_reads),
                                ("chrom_ids", self._chr, n_reads)):
            if arr.shape[0] != want:
                raise ValueError(f"{name} must have {want} entries, has {arr.shape[0]}")
        self._refs = [None if r is None else np.frombuffer(bytes(r), dtype=np.uint8) for r in refs]
        self._ref_ptrs = (_u8p * len(refs))(*[
            ctypes.cast(None, _u8p) if r is None or r.size == 0 else r.ctypes.data_as(_u8p)
            for r in self._refs])
        self._ref_lens = np.array([0 if r is None else r.size for r in self._refs], dtype=np.uint64)
        c = self.ctx
        c.sequences = self._seq.ctypes.data_as(_u8p)
        c.seq_offsets = self._so.ctypes.data_as(_u64p)
        c.cigars = self._cig.ctypes.data_as(_u8p)
        c.cigar_offsets = self._co.ctypes.data_as(_u64p)
        c.positions = self._pos.ctypes.data_as(ctypes.POINTER(ctypes.c_int64))
        c.chrom_ids = self._chr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        c.refs = self._ref_ptrs
        c.ref_lengths = self._ref_lens.ctypes.data_as(_u64p)
        c.n_refs = len(refs)


def encode(tags: Sequence[str], *, sequences=None, seq_offsets=None, cigars=None,
           positions=None, chrom_ids=None, references=None) -> bytes:
    """Encode one block's per-read tag text to a SAM_TAGS blob.

    ``references[c]`` is the bases of chromosome id ``c`` (or ``None``);
    pass no references to store MD/NM verbatim.
    """
    _require_lib("encode")
    n = len(tags)
    ctx = _Context(n, sequences, seq_offsets, cigars, positions, chrom_ids, references)
    buf, off = _joined(tags)
    out = _u8p()
    out_len = ctypes.c_size_t(0)
    rc = _lib.ttio_sam_tags_encode(ctypes.byref(ctx.ctx), buf.ctypes.data_as(_u8p),
                                   off.ctypes.data_as(_u64p), ctypes.byref(out),
                                   ctypes.byref(out_len))
    if rc != 0:
        raise ValueError(_ERR_MESSAGES.get(rc, f"native error {rc}"))
    try:
        return ctypes.string_at(out, out_len.value)
    finally:
        _lib.ttio_sam_tags_free(out)


def decode(blob: bytes, n_reads: int, *, sequences=None, seq_offsets=None, cigars=None,
           positions=None, chrom_ids=None, references=None) -> list[str]:
    """Decode a SAM_TAGS blob to per-read tag text (same context as encode)."""
    _require_lib("decode")
    ctx = _Context(n_reads, sequences, seq_offsets, cigars, positions, chrom_ids, references)
    src = np.frombuffer(bytes(blob) or b"\0", dtype=np.uint8)
    out = _u8p()
    off = np.zeros(n_reads + 1, dtype=np.uint64)
    rc = _lib.ttio_sam_tags_decode(ctypes.byref(ctx.ctx), src.ctypes.data_as(_u8p), len(blob),
                                   ctypes.byref(out), off.ctypes.data_as(_u64p))
    if rc != 0:
        raise ValueError(_ERR_MESSAGES.get(rc, f"native error {rc}"))
    try:
        text = ctypes.string_at(out, int(off[-1])).decode("utf-8")
    finally:
        _lib.ttio_sam_tags_free(out)
    o = off.tolist()
    # Offsets are byte offsets; split on the byte form when non-ASCII.
    if text.isascii():
        return [text[o[i]:o[i + 1]] for i in range(n_reads)]
    raw = text.encode("utf-8")
    return [raw[o[i]:o[i + 1]].decode("utf-8") for i in range(n_reads)]
