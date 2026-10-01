"""M101: SAM optional tags across the three SDKs.

* import  — each SDK imports the same tagged SAM; the tags blob, the
  block index and the sequences blob are byte-identical to Python's,
  and the records round-trip (tags included).
* transport — the 3x3 plaintext transport matrix carries the tags.
* per-AU — the 3x3 encryptor x decryptor matrix restores the tags
  channel byte-identically, on a run whose tags first appear in a later
  block.
* bam_dump — the three dump CLIs agree on the ``tags`` key.

Fixtures come from ``tests/test_m101_sam_tags._make_sam`` (MD/NM by
``samtools calmd``). Cells whose SDK binary is absent skip, the M89
convention.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("h5py")
import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from test_m89_cross_language import (  # type: ignore[import-not-found]  # noqa: E402
    _DECODERS, _ENCODERS, _MATRIX, _resolve_java_tool, _resolve_objc_tool)
from test_m99_blocks_v1_matrix import _PERAU_CLIS  # type: ignore[import-not-found]  # noqa: E402
from test_m101_sam_tags import _make_sam  # type: ignore[import-not-found]  # noqa: E402
from _digests import genomic_run_sam_full_md5, sam_full_md5  # noqa: E402

from ttio.codecs import sam_tags  # noqa: E402
from ttio.importers import registry as im  # noqa: E402
from ttio.spectral_dataset import SpectralDataset  # noqa: E402

pytestmark = [
    pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools"),
    pytest.mark.skipif(not sam_tags.HAVE_NATIVE_LIB, reason="libttio_rans"),
]

RUN = "genomic_0001"
SDKS = ("python", "objc", "java")
# variant -> (_make_sam kwargs, use the reference)
VARIANTS = {
    "refdiff": ({}, True),
    "no_reference": ({}, False),
    "late_tags": ({"tagless_prefix": 250}, True),
}


def _import(sdk: str, sam: Path, ref: Path | None, out: Path) -> None:
    extras = {"block_reads": "100"}
    if ref is not None:
        extras.update(reference=str(ref), embed_reference="1")
    if sdk == "python":
        im.encode("sam", [str(sam)], str(out), block_reads=100,
                  **({"reference": str(ref), "embed_reference": True} if ref else {}))
        return
    if sdk == "objc":
        tool = _resolve_objc_tool("TtioEncode")
        if tool is None:
            pytest.skip("ObjC TtioEncode not built")
        argv, env = [str(tool[0])], tool[1]
    else:
        tool = _resolve_java_tool("global.thalion.ttio.tools.EncodeCli")
        if tool is None:
            pytest.skip("Java classpath not available")
        argv, env = tool
    cmd = [*argv, "--format", "sam", "--input", str(sam), "--output", str(out)]
    for k, v in extras.items():
        cmd += ["--extra", f"{k}={v}"]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    assert proc.returncode == 0, f"{sdk} import failed: {proc.stderr[-800:]}"


def _state(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        rg = f[f"study/genomic_runs/{RUN}"]
        sc = rg["signal_channels"]
        tags = sc.get("tags")
        feats = f.attrs["ttio_features"]
        feats = feats.decode() if isinstance(feats, bytes) else feats
        return {
            "index": rg["blocks/index"][...].tobytes(),
            "tags": None if tags is None else bytes(tags[...]),
            "tags_compression": None if tags is None else int(tags.attrs["compression"]),
            "sequences": bytes(sc["sequences/data"][...]),
            "opt_sam_tags": "opt_sam_tags" in json.loads(feats),
        }


@pytest.mark.parametrize("variant", list(VARIANTS))
@pytest.mark.parametrize("sdk", SDKS)
def test_m101_import_byte_identical(sdk, variant, tmp_path):
    kwargs, use_ref = VARIANTS[variant]
    sam, ref = _make_sam(tmp_path, **kwargs)
    ref = ref if use_ref else None
    want = tmp_path / "python.tio"
    _import("python", sam, ref, want)
    got = tmp_path / f"{sdk}.tio"
    _import(sdk, sam, ref, got)
    w, g = _state(want), _state(got)
    assert g["tags_compression"] == 18 and g["opt_sam_tags"]
    for key in ("index", "tags", "sequences"):
        assert g[key] == w[key], f"{sdk} {variant}: {key} differs from Python"
    with SpectralDataset.open(str(got)) as ds:
        assert genomic_run_sam_full_md5(ds.genomic_runs[RUN]) == sam_full_md5(sam)


@pytest.mark.parametrize("writer,reader", _MATRIX,
                         ids=[f"{w}-encode_{r}-decode" for w, r in _MATRIX])
def test_m101_transport_3x3(writer, reader, tmp_path):
    sam, ref = _make_sam(tmp_path)
    src = tmp_path / "src.tio"
    _import("python", sam, ref, src)
    tis = tmp_path / "cell.tis"
    out = tmp_path / "cell.tio"
    _ENCODERS[writer](src, tis)
    _DECODERS[reader](tis, out)
    with SpectralDataset.open(str(out)) as ds:
        assert genomic_run_sam_full_md5(ds.genomic_runs[RUN]) == sam_full_md5(sam)


_PERAU_MATRIX = [(a, b) for a in SDKS for b in SDKS]


@pytest.mark.parametrize("encryptor,decryptor", _PERAU_MATRIX,
                         ids=[f"{a}-encrypt_{b}-decrypt" for a, b in _PERAU_MATRIX])
def test_m101_per_au_3x3(encryptor, decryptor, tmp_path):
    sam, ref = _make_sam(tmp_path, tagless_prefix=250)
    plain = tmp_path / "plain.tio"
    _import("python", sam, ref, plain)
    before = _state(plain)
    key = tmp_path / "key.bin"
    key.write_bytes(bytes(range(32)))

    enc = tmp_path / "enc.tio"
    _PERAU_CLIS[encryptor](["encrypt", str(plain), str(enc), str(key)])
    with h5py.File(enc, "r") as f:
        sc = f[f"study/genomic_runs/{RUN}/signal_channels"]
        assert "tags" not in sc and "tags_segments" in sc

    work = tmp_path / "work.tio"
    shutil.copyfile(enc, work)
    _PERAU_CLIS[decryptor](["decrypt-in-place", str(work), str(key)])
    after = _state(work)
    for key_ in ("index", "tags", "tags_compression", "sequences"):
        assert after[key_] == before[key_], f"{encryptor}->{decryptor}: {key_}"
    with h5py.File(work, "r") as f:
        assert f[f"study/genomic_runs/{RUN}/signal_channels/tags"].compression is None


def _bam_dump(sdk: str, sam: Path) -> dict:
    if sdk == "python":
        from ttio.importers.bam_dump import dump
        return dump(str(sam))
    if sdk == "objc":
        tool = _resolve_objc_tool("TtioBamDump")
        if tool is None:
            pytest.skip("ObjC TtioBamDump not built")
        argv, env = [str(tool[0])], tool[1]
    else:
        tool = _resolve_java_tool("global.thalion.ttio.importers.BamDump")
        if tool is None:
            pytest.skip("Java classpath not available")
        argv, env = tool
    proc = subprocess.run([*argv, str(sam)], capture_output=True, text=True, env=env, timeout=300)
    assert proc.returncode == 0, proc.stderr[-800:]
    return json.loads(proc.stdout)


@pytest.mark.parametrize("sdk", ["objc", "java"])
def test_m101_bam_dump_tags_match_python(sdk, tmp_path):
    sam, _ = _make_sam(tmp_path, n=200)
    bam = tmp_path / "in.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    for path in (sam, bam):
        want, got = _bam_dump("python", path), _bam_dump(sdk, path)
        assert got["tags"] == want["tags"], f"{sdk} {path.suffix}"
        assert got["read_names"] == want["read_names"]
