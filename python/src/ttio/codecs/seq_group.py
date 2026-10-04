"""Python ctypes wrapper for read grouping (M103, ``ttio_seq_group``).

Returns the permutation that puts an unaligned run's reads from the same
place in the genome next to each other, so a blocks_v1 block holds them
together and SEQ_CM sees each read's overlap partners. ``order[j]`` is the
input index of the read stored at row ``j``. The native kernel breaks every
tie on a total order, so the three SDKs compute the same permutation.
Requires TTIO_RANS_LIB_PATH to point at a built libttio_rans.
"""
from __future__ import annotations

import ctypes
from collections.abc import Sequence

import numpy as np

from .fqzcomp_nx16_z import _HAVE_NATIVE_LIB, _native_lib

HAVE_NATIVE_LIB: bool = _HAVE_NATIVE_LIB

FLAG_FILL = 0x01
FLAG_MATES = 0x02

_u8p = ctypes.POINTER(ctypes.c_uint8)
_u64p = ctypes.POINTER(ctypes.c_uint64)
_u32p = ctypes.POINTER(ctypes.c_uint32)


class _Params(ctypes.Structure):
    _fields_ = [
        ("k", ctypes.c_uint8),
        ("window", ctypes.c_uint8),
        ("min_votes", ctypes.c_uint8),
        ("flags", ctypes.c_uint8),
        ("max_occ", ctypes.c_uint32),
    ]


if HAVE_NATIVE_LIB:
    _lib = _native_lib
    _lib.ttio_seq_group_default_params.argtypes = [ctypes.POINTER(_Params)]
    _lib.ttio_seq_group_default_params.restype = None
    _lib.ttio_seq_group.argtypes = [
        _u8p, _u64p, ctypes.c_uint64, _u8p, _u64p, ctypes.POINTER(_Params), _u32p,
    ]
    _lib.ttio_seq_group.restype = ctypes.c_int


def group(seq: bytes, lengths, names: Sequence[str] | None) -> np.ndarray:
    """The grouping order of reads whose bases are ``seq`` back to back.

    ``names`` (one per read) pair mates by name; ``None`` groups without
    them. Returns a uint32 array: entry ``j`` is the input index of the
    read stored at row ``j``.
    """
    if not HAVE_NATIVE_LIB:
        raise RuntimeError("seq_group.group requires libttio_rans (set TTIO_RANS_LIB_PATH)")
    lens = np.ascontiguousarray(lengths, dtype=np.uint64)
    n = int(lens.size)
    if int(lens.sum()) != len(seq):
        raise ValueError(f"lengths sum to {int(lens.sum())} but there are {len(seq)} bases")
    if n >= 0xFFFFFFFF:
        raise ValueError("read grouping holds at most 2^32 - 2 reads")
    params = _Params()
    _lib.ttio_seq_group_default_params(ctypes.byref(params))
    name_buf = np.zeros(1, dtype=np.uint8)
    name_off = np.zeros(n + 1, dtype=np.uint64)
    if names is None:
        params.flags &= ~FLAG_MATES & 0xFF
    else:
        if len(names) != n:
            raise ValueError(f"{len(names)} names for {n} reads")
        enc = [s.encode("utf-8") if isinstance(s, str) else bytes(s) for s in names]
        if n:
            np.cumsum([len(b) for b in enc], out=name_off[1:])
        name_buf = np.frombuffer(b"".join(enc) or b"\0", dtype=np.uint8)
    order = np.zeros(max(n, 1), dtype=np.uint32)
    src = np.frombuffer(bytes(seq) or b"\0", dtype=np.uint8)
    lp = lens if n else np.zeros(1, dtype=np.uint64)
    rc = _lib.ttio_seq_group(src.ctypes.data_as(_u8p), lp.ctypes.data_as(_u64p), n,
                             name_buf.ctypes.data_as(_u8p), name_off.ctypes.data_as(_u64p),
                             ctypes.byref(params), order.ctypes.data_as(_u32p))
    if rc == -2:
        raise MemoryError("out of memory in ttio_seq_group")
    if rc != 0:
        raise ValueError(f"ttio_seq_group failed (rc={rc})")
    return order[:n]
