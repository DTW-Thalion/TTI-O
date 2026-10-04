"""HMAC-SHA256 digital signatures matching the ObjC reference implementation.

v0.2 (``v1`` in this file) signed over the raw bytes returned by
``H5Dread`` in the dataset's native memory type; on little-endian hosts
that happens to be the same as the canonical form, but signatures would
not validate across host endianness. v0.3 introduces a ``v2`` canonical
form that normalizes atomic numeric datasets to little-endian before
hashing and emits compound records with a fixed per-field layout
(little-endian numerics plus ``u32_le(length) || bytes`` for VL strings).

Stored signatures carry a ``v2:`` prefix; verifiers accept unprefixed
``v1`` signatures for backward compatibility and silently route them
through the native-bytes path.

The format is fully interoperable with the Objective-C reference
implementation — see ``objc/Source/Protection/TTIOSignatureManager.m``.

Cross-language equivalents
--------------------------
Objective-C: ``TTIOSignatureManager`` · Java:
``global.thalion.ttio.protection.SignatureManager``.

API status: Stable.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import warnings
from typing import Any

import h5py
import numpy as np

SIGNATURE_ATTR = "ttio_signature"
PROVENANCE_SIGNATURE_ATTR = "provenance_signature"
SIGNATURE_V2_PREFIX = "v2:"
SIGNATURE_V3_PREFIX = "v3:"  # ML-DSA-87,

def hmac_sha256(data: bytes, key: bytes) -> bytes:
    """Return the raw 32-byte HMAC-SHA256 MAC."""
    return hmac.new(key, data, hashlib.sha256).digest()


def hmac_sha256_b64(data: bytes, key: bytes) -> str:
    """Return the base64-encoded MAC as produced by the ObjC writer."""
    return base64.b64encode(hmac_sha256(data, key)).decode("ascii")


# ---------------------------------------------------- dataset signatures ---


def _dataset_native_bytes(dataset: h5py.Dataset) -> bytes:
    """Return the raw dataset bytes in native type order (v1 path).

    Matches ``H5Dread`` with the file's native memory type, which is
    what the v0.2 ObjC signer hashed over. h5py's ``dataset[()]`` already
    performs the read; we just take the underlying buffer. This path is
    retained purely for backward compatibility with v0.2 files.
    """
    arr = dataset[()]
    if isinstance(arr, np.ndarray):
        return arr.tobytes()
    return np.asarray(arr).tobytes()


def _dataset_canonical_bytes(dataset: h5py.Dataset) -> bytes:
    """Return the canonical little-endian byte stream for ``dataset`` (v2).

    This helper delegates to the storage-provider protocol's
    :meth:`StorageDataset.read_canonical_bytes`, so the byte stream is
    guaranteed bit-identical across HDF5, Memory, and SQLite backends —
    a file signed through any provider verifies through any other. The
    HDF5 backend's zero-copy speed comes from its :meth:`_Dataset.read`
    override, which the shared canonicaliser consumes.
    """
    from .providers.hdf5 import _Dataset as _Hdf5Dataset
    return _Hdf5Dataset(dataset).read_canonical_bytes()


def sign_dataset(
    dataset: Any,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> str:
    """Sign ``dataset`` with a canonical signature.

    For ``algorithm="hmac-sha256"`` (default) ``key`` is a shared secret
    and the output is ``"v2:" + base64(hmac)``. For
    ``algorithm="ml-dsa-87"`` () ``key`` is the 4896-byte
    ML-DSA-87 signing private key and the output is
    ``"v3:" + base64(signature)``. Use :func:`verify_dataset` to
    validate; it dispatches on the stored prefix.

    phase B: ``dataset`` may be either an ``h5py.Dataset``
    (legacy fast path) or a :class:`StorageDataset` from any provider.
    Non-h5py inputs delegate to :func:`sign_storage_dataset` so the
    same signature ends up on Memory / SQLite / Zarr backends.

    : ML-DSA-87 requires the ``[pqc]`` optional extra (Python /
    ObjC backend is ``liboqs``; Java uses Bouncy Castle — see
    :file:`docs/pqc.md`).
    """
    from .providers.base import StorageDataset
    if isinstance(dataset, StorageDataset):
        return sign_storage_dataset(dataset, key, algorithm=algorithm)
    # P3.9: a raw ``h5py.Dataset`` is wrapped in the HDF5 provider's
    # StorageDataset adapter and signed through the provider-agnostic
    # path, so the live signing path no longer touches raw h5py here.
    # The wrapper's ``set_attribute`` writes a byte-identical
    # ``@ttio_signature`` attribute (proven by the cross-path tests).
    from .providers.hdf5 import _Dataset as _Hdf5Dataset
    return sign_storage_dataset(_Hdf5Dataset(dataset), key, algorithm=algorithm)


def verify_dataset(
    dataset: Any,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> bool:
    """Verify the stored ``@ttio_signature`` against ``key``.

    Accepts v0.2 unprefixed HMAC (native-bytes path), v0.3 ``v2:``
    canonical HMAC, and v0.8 ``v3:`` ML-DSA-87 signatures. Uses
    timing-safe comparison for HMAC; ML-DSA verification itself runs
    in constant time via liboqs.

    phase B: ``dataset`` may be either an ``h5py.Dataset``
    or a :class:`StorageDataset`; non-h5py inputs delegate to
    :func:`verify_storage_dataset`.

    The ``algorithm`` keyword tells the verifier what key shape to
    expect. If the on-disk prefix does not match, raises
    :class:`~ttio.cipher_suite.UnsupportedAlgorithmError` so callers
    don't silently pass verification of a file encrypted with a
    different algorithm. For ML-DSA-87, ``key`` is the 2592-byte
    verification public key.
    """
    from .providers.base import StorageDataset
    if isinstance(dataset, StorageDataset):
        return verify_storage_dataset(dataset, key, algorithm=algorithm)
    from . import cipher_suite
    from .providers.hdf5 import _Dataset as _Hdf5Dataset
    # P3.9: wrap the raw h5py.Dataset once and read the stored signature
    # through the provider adapter, mirroring verify_storage_dataset.
    wrapped = _Hdf5Dataset(dataset)
    if not wrapped.has_attribute(SIGNATURE_ATTR):
        return False
    stored = wrapped.get_attribute(SIGNATURE_ATTR)
    if stored is None:
        return False
    if isinstance(stored, bytes):
        stored = stored.decode("utf-8", errors="replace")
    stored = str(stored)

    # v3 (ml-dsa-87) and v2 (hmac) prefixes route entirely through the
    # provider-agnostic verifier — including its algorithm-mismatch
    # guards, which are identical to the previous inline ones.
    if stored.startswith(SIGNATURE_V3_PREFIX) or stored.startswith(
        SIGNATURE_V2_PREFIX
    ):
        return verify_storage_dataset(wrapped, key, algorithm=algorithm)

    # Reject caller passing "ml-dsa-87" against a non-v3 stored blob —
    # saves a confusing empty-verify return. (Same guard as the v2/v3
    # path above, applied here before the v1 native-bytes fallback.)
    if algorithm == "ml-dsa-87":
        raise cipher_suite.UnsupportedAlgorithmError(
            "stored signature is not v3 (ml-dsa-87) — pass "
            "algorithm='hmac-sha256' to verify legacy signatures"
        )

    # --- SOLE remaining raw-h5py island (v1 legacy) -------------------
    # v0.2 legacy path — unprefixed base64 HMAC over the dataset's
    # *native* byte layout (not the canonical stream), which requires
    # _dataset_native_bytes()'s direct h5py read. This is the only
    # signature code path that still touches raw h5py and is scheduled
    # for separate removal at v1.0 per docs/api-stability-v0.8.md §6.
    cipher_suite.validate_key(algorithm, key)
    warnings.warn(
        "verifying an unprefixed v1 HMAC signature — the v0.2 "
        "native-byte fallback is scheduled for removal at v1.0. "
        "Re-sign with algorithm='hmac-sha256' (emits v2: prefix) "
        "or 'ml-dsa-87' (v3:) to stay compatible.",
        DeprecationWarning,
        stacklevel=2,
    )
    expected = hmac_sha256_b64(_dataset_native_bytes(dataset), key)
    return hmac.compare_digest(stored, expected)


def verify_provenance(run_group: h5py.Group, key: bytes) -> bool:
    """Verify the ``@provenance_signature`` over ``@provenance_json``."""
    stored = _read_vl_string_attr(run_group, PROVENANCE_SIGNATURE_ATTR)
    if stored is None:
        return False
    prov_json = _read_vl_string_attr(run_group, "provenance_json") or ""
    expected = hmac_sha256_b64(prov_json.encode("utf-8"), key)
    return hmac.compare_digest(stored, expected)


# run-level convenience helpers for genomic_runs. These walk
# every signal channel (sequences, qualities) and every genomic_index
# column (offsets, lengths, positions, mapping_qualities, flags) and
# sign each dataset individually with sign_dataset(). The chromosomes
# compound is intentionally skipped — sign_dataset signs atomic
# datasets, and the compound row format requires a separate code path
# that the existing canonical-bytes API doesn't yet cover for VL string
# rows.

# M101: the SAM tags channel (a flat uint8 SAM_TAGS stream) is signed
# with the bases it describes, when present.
_GENOMIC_SIGNAL_CHANNELS = ("sequences", "qualities", "tags")
# L1 (Task #82 Phase B.1, 2026-05-01): the M82-era VL-string
# `chromosomes` compound was replaced with a `chromosome_ids` (uint16)
# + `chromosome_names` (VL-compound) pair, eliminating 42 MB of HDF5
# fractal-heap overhead per chr22 file. Both new columns are signed
# in place of the old single column.
_GENOMIC_INDEX_COLUMNS = (
    "offsets", "lengths", "positions", "mapping_qualities", "flags",
    "chromosome_ids", "chromosome_names",
    # blocks_v1_grouped only (M103): the permutation back to input order.
    "input_index",
)


def _signal_datasets(sig_group: h5py.Group, cname: str) -> list[tuple[str, h5py.Dataset]]:
    """The dataset(s) behind a genomic signal channel: the flat dataset,
    or the datasets inside the channel's group (sequences/refdiff_v2
    for a whole-channel run, sequences/data for blocks_v1)."""
    obj = sig_group[cname]
    if isinstance(obj, h5py.Dataset):
        return [(cname, obj)]
    return [(f"{cname}/{k}", v) for k, v in obj.items() if isinstance(v, h5py.Dataset)]


def sign_genomic_run(
    run_group: h5py.Group,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> dict[str, str]:
    """Sign every signal channel + index column under one genomic
    run with one call.

    Returns a dict mapping ``"<sub>/<dataset>"`` (e.g. ``"signal_channels/sequences"``,
    ``"genomic_index/positions"``) to the prefixed signature stored on
    that dataset's ``@ttio_signature`` attribute. Datasets that don't
    exist are silently skipped (e.g. encrypted files have segments
    instead of plaintext signal channels).

    Algorithm dispatches identically to :func:`sign_dataset` —
    ``"hmac-sha256"`` (HMAC + shared key) or ``"ml-dsa-87"`` (PQC
    private key).
    """
    out: dict[str, str] = {}
    if "signal_channels" in run_group:
        sig_group = run_group["signal_channels"]
        for cname in _GENOMIC_SIGNAL_CHANNELS:
            if cname in sig_group:
                for path, ds in _signal_datasets(sig_group, cname):
                    out[f"signal_channels/{path}"] = sign_dataset(
                        ds, key, algorithm=algorithm,
                    )
    if "genomic_index" in run_group:
        idx_group = run_group["genomic_index"]
        for cname in _GENOMIC_INDEX_COLUMNS:
            if cname in idx_group:
                out[f"genomic_index/{cname}"] = sign_dataset(
                    idx_group[cname], key, algorithm=algorithm,
                )
    # blocks_v1 (format-spec 10.12.6): the block index is part of the
    # signed set, after genomic_index and before signal_channels in the
    # documented order (each dataset carries its own signature, so the
    # order only fixes the dict layout).
    if "blocks" in run_group and "index" in run_group["blocks"]:
        out["blocks/index"] = sign_dataset(
            run_group["blocks"]["index"], key, algorithm=algorithm,
        )
    return out


def verify_genomic_run(
    run_group: h5py.Group,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> bool:
    """Verify every dataset that :func:`sign_genomic_run` would sign.

    Returns ``True`` iff every present, signed dataset verifies under
    the key. A dataset that has no ``@ttio_signature`` (i.e. wasn't
    signed) returns ``False`` from :func:`verify_dataset` and so is
    treated as a verification failure — that's intentional, since a
    partial-signature run is not a fully-signed run. Datasets that
    don't exist on disk are skipped.
    """
    if "signal_channels" in run_group:
        sig_group = run_group["signal_channels"]
        for cname in _GENOMIC_SIGNAL_CHANNELS:
            if cname in sig_group:
                for _path, ds in _signal_datasets(sig_group, cname):
                    if not verify_dataset(ds, key, algorithm=algorithm):
                        return False
    if "genomic_index" in run_group:
        idx_group = run_group["genomic_index"]
        for cname in _GENOMIC_INDEX_COLUMNS:
            if cname in idx_group and not verify_dataset(
                idx_group[cname], key, algorithm=algorithm,
            ):
                return False
    if "blocks" in run_group and "index" in run_group["blocks"]:
        if not verify_dataset(run_group["blocks"]["index"], key, algorithm=algorithm):
            return False
    return True


# M98: assembly graphs. The signed set is the segments/sequences byte
# channel — signatures cover it the way they cover genomic-run
# channels. The compounds (records, links, paths, extras, line_index)
# are skipped the way the genomic chromosomes compound was: VL-string
# compound rows are outside the canonical-bytes API.


def sign_assembly_graph(
    graph_group: h5py.Group,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> dict[str, str]:
    """Sign the byte channels of one assembly graph
    (``/study/assembly_graphs/<name>/``).

    Returns a dict mapping ``"segments/sequences"`` to the prefixed
    signature stored on the dataset's ``@ttio_signature`` attribute.
    Absent datasets are skipped (a graph whose segments all carry
    ``*`` has no sequences channel).
    """
    out: dict[str, str] = {}
    if "segments" in graph_group:
        seg_group = graph_group["segments"]
        if "sequences" in seg_group:
            out["segments/sequences"] = sign_dataset(
                seg_group["sequences"], key, algorithm=algorithm,
            )
    return out


def verify_assembly_graph(
    graph_group: h5py.Group,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> bool:
    """Verify every dataset :func:`sign_assembly_graph` would sign.

    Returns ``True`` iff every present, signed dataset verifies under
    the key; an unsigned present dataset fails, matching
    :func:`verify_genomic_run`. Absent datasets are skipped.
    """
    if "segments" in graph_group:
        seg_group = graph_group["segments"]
        if "sequences" in seg_group and not verify_dataset(
            seg_group["sequences"], key, algorithm=algorithm,
        ):
            return False
    return True


# ── Provider-agnostic sign / verify () ─────────────────────────


def sign_storage_dataset(
    dataset: Any,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> str:
    """Sign any :class:`~ttio.providers.base.StorageDataset` with the
    named algorithm.

    Mirrors :func:`sign_dataset` (h5py-native) but works across every
    provider — HDF5, Memory, SQLite, Zarr. The canonical byte stream
    comes from :meth:`StorageDataset.read_canonical_bytes` so the
    signature is identical regardless of backend.

    Used by the M54.1 cross-language conformance matrix. The dataset
    must support ``set_attribute("ttio_signature", str)`` — every
    shipping provider does.

    @since 0.8 ()
    """
    from . import cipher_suite
    canonical = dataset.read_canonical_bytes()
    if algorithm == "hmac-sha256":
        cipher_suite.validate_key(algorithm, key)
        mac_b64 = hmac_sha256_b64(canonical, key)
        prefixed = SIGNATURE_V2_PREFIX + mac_b64
    elif algorithm == "ml-dsa-87":
        cipher_suite.validate_private_key(algorithm, key)
        from . import pqc
        sig = pqc.sig_sign(key, canonical)
        prefixed = SIGNATURE_V3_PREFIX + base64.b64encode(sig).decode("ascii")
    else:
        raise cipher_suite.UnsupportedAlgorithmError(
            f"{algorithm}: signature path not yet implemented"
        )
    dataset.set_attribute(SIGNATURE_ATTR, prefixed)
    return prefixed


def verify_storage_dataset(
    dataset: Any,
    key: bytes,
    *,
    algorithm: str = "hmac-sha256",
) -> bool:
    """Verify any :class:`~ttio.providers.base.StorageDataset`'s
    stored ``@ttio_signature``.

    Prefix/algorithm must match — a ``v3:`` attribute with
    ``algorithm="hmac-sha256"`` raises
    :class:`~ttio.cipher_suite.UnsupportedAlgorithmError` so silent
    acceptance of a wrong-scheme file is impossible.

    @since 0.8 ()
    """
    from . import cipher_suite
    if not dataset.has_attribute(SIGNATURE_ATTR):
        return False
    stored = dataset.get_attribute(SIGNATURE_ATTR)
    if stored is None:
        return False
    if isinstance(stored, bytes):
        stored = stored.decode("utf-8", errors="replace")
    stored = str(stored)

    canonical = dataset.read_canonical_bytes()

    if stored.startswith(SIGNATURE_V3_PREFIX):
        if algorithm != "ml-dsa-87":
            raise cipher_suite.UnsupportedAlgorithmError(
                f"stored signature is v3 (ml-dsa-87) but caller "
                f"passed algorithm={algorithm!r}"
            )
        cipher_suite.validate_public_key(algorithm, key)
        from . import pqc
        sig = base64.b64decode(stored[len(SIGNATURE_V3_PREFIX):])
        return pqc.sig_verify(key, canonical, sig)

    if algorithm == "ml-dsa-87":
        raise cipher_suite.UnsupportedAlgorithmError(
            "stored signature is not v3 (ml-dsa-87) — pass "
            "algorithm='hmac-sha256' to verify legacy signatures"
        )
    cipher_suite.validate_key(algorithm, key)
    if stored.startswith(SIGNATURE_V2_PREFIX):
        payload = stored[len(SIGNATURE_V2_PREFIX):]
    else:
        payload = stored
    expected = hmac_sha256_b64(canonical, key)
    return hmac.compare_digest(payload, expected)


# --------------------------------------------------- VL string attr helpers ---


def _write_vl_string_attr(obj: Any, name: str, value: str) -> None:
    """Write a variable-length UTF-8 string attribute (for parity with the
    ObjC signature writer which uses ``H5T_VARIABLE``).
    """
    nbytes = name.encode("utf-8")
    if h5py.h5a.exists(obj.id, nbytes):
        h5py.h5a.delete(obj.id, nbytes)
    tid = h5py.h5t.C_S1.copy()
    tid.set_size(h5py.h5t.VARIABLE)
    tid.set_strpad(h5py.h5t.STR_NULLTERM)
    tid.set_cset(h5py.h5t.CSET_UTF8)
    space = h5py.h5s.create(h5py.h5s.SCALAR)
    aid = h5py.h5a.create(obj.id, nbytes, tid, space)
    aid.write(np.array([value.encode("utf-8")], dtype=h5py.string_dtype()))
    aid.close()


def _read_vl_string_attr(obj: Any, name: str) -> str | None:
    if name not in obj.attrs:
        return None
    raw = obj.attrs[name]
    if isinstance(raw, bytes):
        return raw.decode("utf-8")
    if isinstance(raw, np.bytes_):
        return raw.tobytes().decode("utf-8")
    if isinstance(raw, str):
        return raw
    return str(raw)
