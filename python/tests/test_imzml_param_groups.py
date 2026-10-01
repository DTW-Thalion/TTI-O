"""imzML arrays declared through referenceableParamGroups.

Many writers (the HR2MSI PXD001283 files among them) state each
binaryDataArray's kind, precision and compression in a
<referenceableParamGroup> and reference it with
<referenceableParamGroupRef>, keeping only the external offset/length
cvParams inline. The reader must resolve the reference to know which
array an offset belongs to.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ttio.exporters import imzml as imzml_writer
from ttio.importers.imzml import ImzMLPixelSpectrum, read

_KIND = ("MS:1000514", "MS:1000515", "MS:1000521", "MS:1000523", "MS:1000574", "MS:1000576")


def _write_inline(path: Path) -> list[ImzMLPixelSpectrum]:
    rng = np.random.default_rng(5)
    pixels = []
    for i in range(12):
        n = 5 + i
        pixels.append(ImzMLPixelSpectrum(
            x=i % 4 + 1, y=i // 4 + 1, z=1,
            mz=np.sort(rng.uniform(100, 900, n)),
            intensity=rng.uniform(0, 1e4, n)))
    imzml_writer.write(pixels, path, mode="processed", grid_max_x=4, grid_max_y=3,
                       pixel_size_x=10.0, pixel_size_y=10.0)
    return pixels


def _to_param_groups(src: Path, dst: Path) -> None:
    doc = src.read_text(encoding="utf-8")
    groups: dict[str, str] = {}

    def move(m: re.Match) -> str:
        block = m.group(0)
        name = "mzArray" if "MS:1000514" in block else "intensityArray"
        params = [p for p in re.findall(r"\s*<cvParam [^>]*/>", block)
                  if any(a in p for a in _KIND)]
        groups.setdefault(name, "".join(params))
        for p in params:
            block = block.replace(p, "", 1)
        return re.sub(r"(<binaryDataArray\b[^>]*>)",
                      lambda b: b.group(1) + f'<referenceableParamGroupRef ref="{name}"/>',
                      block, count=1)

    doc = re.sub(r"<binaryDataArray\b.*?</binaryDataArray>", move, doc, flags=re.S)
    assert "MS:1000514" not in doc.split("<run")[1], "array kind still inline"
    group_list = (f'<referenceableParamGroupList count="{len(groups)}">'
                  + "".join(f'<referenceableParamGroup id="{k}">{v}</referenceableParamGroup>'
                            for k, v in groups.items())
                  + "</referenceableParamGroupList>")
    doc = re.sub(r"(</fileDescription>)", lambda m: m.group(1) + group_list, doc, count=1)
    dst.write_text(doc, encoding="utf-8")
    dst.with_suffix(".ibd").write_bytes(src.with_suffix(".ibd").read_bytes())


def test_arrays_declared_by_param_group_ref(tmp_path: Path) -> None:
    inline = tmp_path / "inline.imzML"
    pixels = _write_inline(inline)
    grouped = tmp_path / "grouped.imzML"
    _to_param_groups(inline, grouped)

    got = read(grouped).spectra
    want = read(inline).spectra
    assert len(got) == len(want) == len(pixels)
    for g, w, p in zip(got, want, pixels):
        assert (g.x, g.y) == (w.x, w.y)
        np.testing.assert_array_equal(g.mz, w.mz)
        np.testing.assert_array_equal(g.intensity, w.intensity)
        np.testing.assert_allclose(g.mz, p.mz)
