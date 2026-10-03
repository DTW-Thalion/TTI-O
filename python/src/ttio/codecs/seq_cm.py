"""Python ctypes wrapper for SEQ_CM (codec id 19, M103).

Read bases coded without a reference by a context-mixing model (orders
11/16/24 with reverse-complement training) and a binary arithmetic coder.
The blob is written as ``signal_channels/sequences`` with
``@compression = 19``. Spec: ``docs/codecs/seq_cm.md``.

The codec is context-aware: decode needs the reads' lengths, which the run
already stores. Only the default parameters are used, and the encoder picks
the model's table size from the base count, so the same reads always code
to the same bytes (a per-AU decrypt restore re-encodes and must match).
Requires TTIO_RANS_LIB_PATH to point at a built libttio_rans.
"""
from __future__ import annotations

import ctypes

import numpy as np

from .fqzcomp_nx16_z import _HAVE_NATIVE_LIB, _native_lib

HAVE_NATIVE_LIB: bool = _HAVE_NATIVE_LIB

ERR_PARAM = -1
ERR_ALLOC = -2
ERR_CORRUPT = -3

_ERR_MESSAGES = {
    ERR_PARAM: "invalid parameters (lengths that do not match the bases or the blob)",
    ERR_ALLOC: "out of memory in native code",
    ERR_CORRUPT: "corrupt SEQ_CM blob",
}

_u8p = ctypes.POINTER(ctypes.c_uint8)
_u64p = ctypes.POINTER(ctypes.c_uint64)


class _Params(ctypes.Structure):
    _fields_ = [
        ("n_orders", ctypes.c_uint8),
        ("orders", ctypes.c_uint8 * 3),
        ("table_bits", ctypes.c_uint8),
        ("flags", ctypes.c_uint8),
        ("limit", ctypes.c_uint16),
        ("lr", ctypes.c_uint16),
    ]


if HAVE_NATIVE_LIB:
    _lib = _native_lib
    _lib.ttio_seq_cm_encode.argtypes = [
        _u8p, _u64p, ctypes.c_uint64, ctypes.POINTER(_Params),
        ctypes.POINTER(_u8p), ctypes.POINTER(ctypes.c_size_t),
    ]
    _lib.ttio_seq_cm_encode.restype = ctypes.c_int
    _lib.ttio_seq_cm_decode.argtypes = [
        _u8p, ctypes.c_size_t, _u64p, ctypes.c_uint64,
        ctypes.POINTER(_u8p), ctypes.POINTER(ctypes.c_size_t),
    ]
    _lib.ttio_seq_cm_decode.restype = ctypes.c_int
    _lib.ttio_seq_cm_free.argtypes = [ctypes.c_void_p]
    _lib.ttio_seq_cm_free.restype = None


def _require_lib(what: str) -> None:
    if not HAVE_NATIVE_LIB:
        raise RuntimeError(f"seq_cm.{what} requires libttio_rans (set TTIO_RANS_LIB_PATH)")


def _lengths(lengths) -> np.ndarray:
    arr = np.ascontiguousarray(lengths, dtype=np.uint64)
    if arr.ndim != 1:
        raise ValueError("lengths must be one-dimensional")
    return arr


def encode(seq: bytes, lengths) -> bytes:
    """Encode the bases of ``len(lengths)`` reads, back to back in ``seq``."""
    _require_lib("encode")
    lens = _lengths(lengths)
    if int(lens.sum()) != len(seq):
        raise ValueError(f"lengths sum to {int(lens.sum())} but there are {len(seq)} bases")
    src = np.frombuffer(bytes(seq) or b"\0", dtype=np.uint8)
    lp = lens if lens.size else np.zeros(1, dtype=np.uint64)
    out = _u8p()
    out_len = ctypes.c_size_t(0)
    rc = _lib.ttio_seq_cm_encode(src.ctypes.data_as(_u8p), lp.ctypes.data_as(_u64p),
                                 lens.size, None, ctypes.byref(out), ctypes.byref(out_len))
    if rc != 0:
        raise ValueError(_ERR_MESSAGES.get(rc, f"native error {rc}"))
    try:
        return ctypes.string_at(out, out_len.value)
    finally:
        _lib.ttio_seq_cm_free(out)


def decode(blob: bytes, lengths) -> bytes:
    """Decode a SEQ_CM blob to the reads' bases, given their lengths."""
    _require_lib("decode")
    lens = _lengths(lengths)
    src = np.frombuffer(bytes(blob) or b"\0", dtype=np.uint8)
    lp = lens if lens.size else np.zeros(1, dtype=np.uint64)
    out = _u8p()
    out_len = ctypes.c_size_t(0)
    rc = _lib.ttio_seq_cm_decode(src.ctypes.data_as(_u8p), len(blob), lp.ctypes.data_as(_u64p),
                                 lens.size, ctypes.byref(out), ctypes.byref(out_len))
    if rc != 0:
        raise ValueError(_ERR_MESSAGES.get(rc, f"native error {rc}"))
    try:
        return ctypes.string_at(out, out_len.value)
    finally:
        _lib.ttio_seq_cm_free(out)
