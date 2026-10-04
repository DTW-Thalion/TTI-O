"""M103: grouped runs (layout blocks_v1_grouped) across the three SDKs.

* import — each SDK imports the same unaligned paired SAM with
  ``group_reads=1``; layout, block index, every blob and
  ``genomic_index/input_index`` are byte-identical to Python's, and the
  reads read back in SAM order.
* transport — the 3x3 plaintext transport matrix over a grouped file
  delivers the reads in input order (the readers' input-order iteration).
* per-AU — the 3x3 encryptor x decryptor matrix restores the grouped run
  byte-identically, ``input_index`` included.
* encrypted transport — the 3x3 sender x receiver matrix carries
  ``input_index`` in the BlockSidecars (transport-spec v0.13).

Cells whose SDK binary is absent skip, the M89 convention.
"""
from __future__ import annotations

import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("h5py")
import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from test_m89_cross_language import (  # type: ignore[import-not-found]  # noqa: E402
    _DECODERS, _ENCODERS, _MATRIX, _resolve_java_tool, _resolve_objc_tool)
from test_m99_blocks_v1_matrix import _PERAU_CLIS  # type: ignore[import-not-found]  # noqa: E402

from ttio.codecs import seq_group  # noqa: E402
from ttio.importers import registry as im  # noqa: E402
from ttio.spectral_dataset import SpectralDataset  # noqa: E402

pytestmark = pytest.mark.skipif(not seq_group.HAVE_NATIVE_LIB, reason="libttio_rans")

RUN = "genomic_0001"
SDKS = ("python", "objc", "java")


def _comp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _make_sam(tmp: Path, pairs: int = 900, L: int = 100, G: int = 30_000) -> tuple[Path, list[str]]:
    """Unaligned read pairs from both strands of a random genome, pairs in
    random order, mates adjacent; returns the SAM and its read names in
    record order."""
    rnd = random.Random(103)
    genome = "".join(rnd.choice("ACGT") for _ in range(G))
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    names = []
    order = list(range(pairs))
    rnd.shuffle(order)
    for k in order:
        p = rnd.randrange(0, G - 400)
        r1 = genome[p:p + L]
        r2 = _comp(genome[p + 250:p + 250 + L])
        if k % 41 == 0:
            r1 = r1[:30] + "N" + r1[31:]
        q = "".join(chr(33 + rnd.randrange(2, 40)) for _ in range(L))
        name = f"frag{k}"
        lines.append(f"{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1}\t{q}")
        lines.append(f"{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2}\t{q}")
        names += [name, name]
    sam = tmp / "in.sam"
    sam.write_text("\n".join(lines) + "\n")
    return sam, names


def _import(sdk: str, sam: Path, out: Path) -> None:
    if sdk == "python":
        im.encode("sam", [str(sam)], str(out), block_reads=100, group_reads="1")
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
    cmd = [*argv, "--format", "sam", "--input", str(sam), "--output", str(out),
           "--extra", "block_reads=100", "--extra", "group_reads=1"]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    assert proc.returncode == 0, f"{sdk} import failed: {proc.stderr[-800:]}"


def _state(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        rg = f[f"study/genomic_runs/{RUN}"]
        sc = rg["signal_channels"]
        layout = rg.attrs["layout"]
        out = {
            "layout": layout.decode() if isinstance(layout, bytes) else str(layout),
            "index": rg["blocks/index"][...].tobytes(),
            "input_index": bytes(rg["genomic_index/input_index"][...].astype("<u4").tobytes()),
        }
        for name in ("sequences/data", "qualities", "read_names", "cigars", "mate_info/inline_v2"):
            out[name] = bytes(sc[name][...]) if name in sc else None
        return out


def _names(path: Path) -> list[str]:
    with SpectralDataset.open(str(path)) as ds:
        g = next(iter(ds.genomic_runs.values()))
        return [r.read_name for r in g.iter_reads()]


@pytest.mark.parametrize("sdk", SDKS)
def test_m103_grouped_import_byte_identical(sdk, tmp_path):
    sam, names = _make_sam(tmp_path)
    want = tmp_path / "python.tio"
    _import("python", sam, want)
    w = _state(want)
    assert w["layout"] == "blocks_v1_grouped"
    ii = np.frombuffer(w["input_index"], dtype="<u4")
    assert sorted(ii.tolist()) == list(range(len(names))) and ii.tolist() != list(range(len(names)))
    got = tmp_path / f"{sdk}.tio"
    _import(sdk, sam, got)
    g = _state(got)
    for key in w:
        assert g[key] == w[key], f"{sdk}: {key} differs from Python"
    assert _names(got) == names


@pytest.mark.parametrize("writer,reader", _MATRIX,
                         ids=[f"{w}-encode_{r}-decode" for w, r in _MATRIX])
def test_m103_grouped_transport_3x3(writer, reader, tmp_path):
    sam, names = _make_sam(tmp_path)
    src = tmp_path / "src.tio"
    _import("python", sam, src)
    tis = tmp_path / "cell.tis"
    out = tmp_path / "cell.tio"
    _ENCODERS[writer](src, tis)
    _DECODERS[reader](tis, out)
    assert _names(out) == names, f"{writer}->{reader}: reads not in input order"


_PAIRS = [(a, b) for a in SDKS for b in SDKS]


@pytest.mark.parametrize("encryptor,decryptor", _PAIRS,
                         ids=[f"{a}-encrypt_{b}-decrypt" for a, b in _PAIRS])
def test_m103_grouped_per_au_3x3(encryptor, decryptor, tmp_path):
    sam, names = _make_sam(tmp_path)
    plain = tmp_path / "plain.tio"
    _import("python", sam, plain)
    before = _state(plain)
    key = tmp_path / "key.bin"
    key.write_bytes(bytes(range(32)))
    enc = tmp_path / "enc.tio"
    _PERAU_CLIS[encryptor](["encrypt", str(plain), str(enc), str(key)])
    work = tmp_path / "work.tio"
    shutil.copyfile(enc, work)
    _PERAU_CLIS[decryptor](["decrypt-in-place", str(work), str(key)])
    after = _state(work)
    for k in before:
        assert after[k] == before[k], f"{encryptor}->{decryptor}: {k}"
    assert _names(work) == names


@pytest.mark.parametrize("sender,receiver", _PAIRS,
                         ids=[f"{a}-send_{b}-recv" for a, b in _PAIRS])
def test_m103_grouped_encrypted_transport_3x3(sender, receiver, tmp_path):
    sam, names = _make_sam(tmp_path)
    plain = tmp_path / "plain.tio"
    _import("python", sam, plain)
    before = _state(plain)
    key = tmp_path / "key.bin"
    key.write_bytes(bytes(range(32)))
    enc = tmp_path / "enc.tio"
    _PERAU_CLIS["python"](["encrypt", str(plain), str(enc), str(key)])
    stream = tmp_path / "s.tis"
    received = tmp_path / "received.tio"
    _PERAU_CLIS[sender](["send", str(enc), str(stream)])
    _PERAU_CLIS[receiver](["recv", str(stream), str(received)])
    _PERAU_CLIS["python"](["decrypt-in-place", str(received), str(key)])
    after = _state(received)
    for k in before:
        assert after[k] == before[k], f"{sender}->{receiver}: {k}"
    assert _names(received) == names
