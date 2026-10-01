"""M101: SAM optional tags through import, storage, read and export.

The SAM fixture is built in the test: random reference, reads with
mismatches, deletions, insertions, soft clips and unmapped reads, and
MD/NM written by ``samtools calmd`` as the independent ground truth.
"""
from __future__ import annotations

import random
import shutil
import subprocess
from pathlib import Path

import pytest

from _digests import genomic_run_sam_full_md5, sam_full_md5, sam11_md5, genomic_run_sam11_md5
from ttio import _hdf5_io as io
from ttio.codecs import sam_tags
from ttio.enums import Compression
from ttio.exporters import registry as ex
from ttio.importers import registry as im
from ttio.spectral_dataset import SpectralDataset

needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools")
pytestmark = [needs_samtools,
              pytest.mark.skipif(not sam_tags.HAVE_NATIVE_LIB, reason="libttio_rans")]

REF_LEN = 20_000


def _make_sam(tmp: Path, n: int = 600, *, tags: bool = True, tagless_prefix: int = 0,
              seed: int = 11) -> tuple[Path, Path]:
    """Write ref.fa and in.sam (coordinate-sorted, one chromosome).
    The first ``tagless_prefix`` reads carry no optional fields."""
    rng = random.Random(seed)
    ref_seq = "".join(rng.choice("ACGT") for _ in range(REF_LEN))
    ref = tmp / "ref.fa"
    ref.write_text(">chr1\n" + "\n".join(ref_seq[i:i + 60] for i in range(0, REF_LEN, 60)) + "\n")
    subprocess.run(["samtools", "faidx", str(ref)], check=True)

    recs = []
    for i in range(n):
        pos = 1 + (i * (REF_LEN - 400)) // n
        kind = i % 6
        q = "".join(rng.choice("#+5:?FI") for _ in range(100))
        seq = list(ref_seq[pos - 1:pos - 1 + 100])
        cigar = "100M"
        if kind == 1:                                   # mismatches
            for k in rng.sample(range(100), 3):
                seq[k] = {"A": "C", "C": "G", "G": "T", "T": "A"}[seq[k]]
        elif kind == 2:                                 # deletion
            seq = list(ref_seq[pos - 1:pos + 39] + ref_seq[pos + 44:pos + 104])
            cigar = "40M5D60M"
        elif kind == 3:                                 # insertion
            seq = list(ref_seq[pos - 1:pos + 49] + "TTT" + ref_seq[pos + 49:pos + 96])
            cigar = "50M3I47M"
        elif kind == 4:                                 # soft clip
            seq = list("GGGGG" + ref_seq[pos - 1:pos + 94])
            cigar = "5S95M"
        flag = 99 if i % 2 == 0 else 147
        recs.append([f"r{i:05d}", str(flag), "chr1", str(pos), "60", cigar, "=",
                     str(pos + 150), "250" if flag == 99 else "-250", "".join(seq), q])
    recs.append(["u0", "4", "*", "0", "0", "*", "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII"])

    header = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:chr1\tLN:{REF_LEN}",
              "@RG\tID:g1\tSM:S1\tPL:ILLUMINA"]
    raw = tmp / "raw.sam"
    lines = []
    for k, r in enumerate(recs):
        extra = []
        if tags and k >= tagless_prefix:
            a = rng.randrange(0, 300)
            extra = [f"AS:i:{a}", f"UQ:i:{a}", "RG:Z:g1", f"XS:i:{rng.randrange(-5, 5)}"]
            if k % 7 == 0:
                extra.append(f"XA:Z:chr1,+{rng.randrange(1, 9999)},100M,{rng.randrange(0, 4)};")
            if k % 11 == 0:
                extra.append("BC:B:c,-3,0,7")
            if k % 13 == 0:
                extra.append("XF:f:0.25")
        lines.append("\t".join(r + extra))
    raw.write_text("\n".join(header + lines) + "\n")
    sam = tmp / "in.sam"
    if tags:
        # calmd appends NM/MD where they are missing (all tagged reads).
        with open(sam, "w") as f:
            subprocess.run(["samtools", "calmd", str(raw), str(ref)], stdout=f,
                           stderr=subprocess.DEVNULL, check=True)
        if tagless_prefix:
            # Strip calmd's tags from the tag-less prefix again.
            out = []
            for line in sam.read_text().splitlines():
                c = line.split("\t")
                if not line.startswith("@") and c[0].startswith("r") and int(c[0][1:]) < tagless_prefix:
                    c = c[:11]
                out.append("\t".join(c))
            sam.write_text("\n".join(out) + "\n")
    else:
        shutil.copy(raw, sam)
    return sam, ref


def _tags_ds(g):
    sc = g.group.open_group("signal_channels")
    return sc.open_dataset("tags") if sc.has_child("tags") else None


def test_import_with_reference_derives_md_nm(tmp_path):
    sam, ref = _make_sam(tmp_path)
    out = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(out), reference=str(ref), embed_reference=True,
              block_reads=200)
    with SpectralDataset.open(str(out)) as ds:
        g = ds.genomic_runs["genomic_0001"]
        assert genomic_run_sam_full_md5(g) == sam_full_md5(sam)
        assert "opt_sam_tags" in ds.feature_flags.features
        rows = g.group.open_group("blocks").open_dataset("index").read_rows()
        assert "tags_off" in rows[0]
        assert all(int(r["tags_codec"]) == int(Compression.SAM_TAGS) for r in rows if int(r["tags_len"]))
        stored = int(_tags_ds(g).length)
    # Recomputing MD/NM keeps the channel well under the tag text.
    text = sum(len(line.split("\t", 11)[11]) for line in sam.read_text().splitlines()
               if not line.startswith("@") and line.count("\t") > 10)
    assert stored * 8 < text


def test_import_without_reference(tmp_path):
    sam, _ = _make_sam(tmp_path)
    out = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(out), block_reads=250)
    with SpectralDataset.open(str(out)) as ds:
        assert genomic_run_sam_full_md5(ds.genomic_runs["genomic_0001"]) == sam_full_md5(sam)


def test_export_writes_the_tags(tmp_path):
    sam, ref = _make_sam(tmp_path)
    tio = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(tio), reference=str(ref), embed_reference=True,
              block_reads=200)
    out = tmp_path / "o.bam"
    ex.export("bam", str(tio), "genomic_0001", str(out))
    assert sam_full_md5(out) == sam_full_md5(sam)


def test_tagless_input_writes_no_tags_channel(tmp_path):
    sam, ref = _make_sam(tmp_path, tags=False)
    out = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(out), reference=str(ref), embed_reference=True,
              block_reads=200)
    with SpectralDataset.open(str(out)) as ds:
        g = ds.genomic_runs["genomic_0001"]
        assert _tags_ds(g) is None
        assert "opt_sam_tags" not in ds.feature_flags.features
        rows = g.group.open_group("blocks").open_dataset("index").read_rows()
        assert "tags_off" not in rows[0]
        assert genomic_run_sam11_md5(g) == sam11_md5(sam)
        assert all(r.tags == "" for r in g.iter_reads())


def test_tags_first_appearing_in_a_later_block(tmp_path):
    """The index gains its tags columns when the first tagged block
    arrives; the earlier blocks keep an empty tags range."""
    sam, ref = _make_sam(tmp_path, tagless_prefix=250)
    out = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(out), reference=str(ref), embed_reference=True,
              block_reads=100)
    with SpectralDataset.open(str(out)) as ds:
        g = ds.genomic_runs["genomic_0001"]
        rows = g.group.open_group("blocks").open_dataset("index").read_rows()
        assert [int(r["tags_len"]) == 0 for r in rows][:2] == [True, True]
        assert int(rows[-1]["tags_len"]) > 0
        assert genomic_run_sam_full_md5(g) == sam_full_md5(sam)


KEY = b"\x65" * 32


def _channel_state(path):
    import h5py
    with h5py.File(path, "r") as f:
        rg = f["study/genomic_runs/genomic_0001"]
        sc = rg["signal_channels"]
        return (rg["blocks/index"][...].tobytes(), bytes(sc["tags"][...]),
                bytes(sc["sequences/data"][...]), bytes(sc["qualities"][...]))


def test_per_au_encrypts_tags_and_restores_byte_identical(tmp_path):
    import h5py
    from ttio.encryption_per_au import decrypt_per_au, decrypt_per_au_in_place, encrypt_per_au
    sam, ref = _make_sam(tmp_path)
    path = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(path), reference=str(ref), embed_reference=True,
              block_reads=200)
    before = _channel_state(path)
    with SpectralDataset.open(str(path)) as ds:
        want = [r.tags for r in ds.genomic_runs["genomic_0001"].iter_reads()]

    encrypt_per_au(str(path), KEY)
    with h5py.File(path, "r") as f:
        sc = f["study/genomic_runs/genomic_0001/signal_channels"]
        assert "tags" not in sc and "tags_segments" in sc
    assert decrypt_per_au(str(path), KEY)["genomic_0001"]["tags"] == want

    decrypt_per_au_in_place(str(path), KEY)
    assert _channel_state(path) == before
    with SpectralDataset.open(str(path)) as ds:
        assert genomic_run_sam_full_md5(ds.genomic_runs["genomic_0001"]) == sam_full_md5(sam)


def test_per_au_refuses_whole_channel_tags(tmp_path):
    import dataclasses
    from ttio.encryption_per_au import encrypt_per_au
    from ttio.importers.bam import BamReader
    sam, _ = _make_sam(tmp_path, n=60)
    run = dataclasses.replace(BamReader(sam).to_genomic_run(), opt_legacy_whole_channel=True)
    path = tmp_path / "o.tio"
    SpectralDataset.write_minimal(str(path), title="t", isa_investigation_id="i",
                                  runs={}, genomic_runs={"g": run})
    with pytest.raises(ValueError, match="SAM tags"):
        encrypt_per_au(str(path), KEY)


def test_signatures_cover_the_tags(tmp_path):
    import h5py
    from ttio.signatures import sign_genomic_run, verify_genomic_run
    sam, ref = _make_sam(tmp_path, n=120)
    path = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(path), reference=str(ref), embed_reference=True)
    key = b"k" * 32
    with h5py.File(path, "r+") as f:
        rg = f["study/genomic_runs/genomic_0001"]
        assert "signal_channels/tags" in sign_genomic_run(rg, key)
        assert verify_genomic_run(rg, key)
        data = rg["signal_channels/tags"][...]
        data[-1] ^= 1
        rg["signal_channels/tags"][...] = data
        assert not verify_genomic_run(rg, key)


def test_plaintext_transport_carries_tags(tmp_path):
    from ttio.transport.codec import file_to_transport, transport_to_file
    sam, ref = _make_sam(tmp_path)
    src = tmp_path / "src.tio"
    im.encode("sam", [str(sam)], str(src), reference=str(ref), embed_reference=True,
              block_reads=200)
    stream = tmp_path / "s.tis"
    file_to_transport(src, stream)
    out = tmp_path / "out.tio"
    transport_to_file(stream, out)
    with SpectralDataset.open(str(out)) as ds:
        assert genomic_run_sam_full_md5(ds.genomic_runs["genomic_0001"]) == sam_full_md5(sam)
        assert "opt_sam_tags" in ds.feature_flags.features


def test_tagless_plaintext_stream_is_unchanged(tmp_path):
    """No tags channel in the DatasetHeader or the AUs of a tag-less run."""
    from ttio.transport.codec import TransportReader, file_to_transport
    from ttio.transport.packets import AccessUnit, PacketType
    sam, _ = _make_sam(tmp_path, n=60, tags=False)
    src = tmp_path / "src.tio"
    im.encode("sam", [str(sam)], str(src))
    stream = tmp_path / "s.tis"
    file_to_transport(src, stream)
    with TransportReader(stream) as tr:
        for hdr, payload in tr.iter_packets():
            if hdr.packet_type == PacketType.ACCESS_UNIT:
                assert "tags" not in {c.name for c in AccessUnit.from_bytes(payload).channels}


def test_encrypted_transport_restores_tags(tmp_path, monkeypatch):
    """The stream does not carry embedded reference bytes, so the
    receiver resolves the REF_DIFF reference (and with it the MD/NM
    derivation) through REF_PATH."""
    from ttio.encryption_per_au import decrypt_per_au_in_place, encrypt_per_au
    from ttio.transport.codec import TransportWriter
    from ttio.transport.encrypted import read_encrypted_to_file, write_encrypted_dataset
    sam, ref = _make_sam(tmp_path)
    monkeypatch.setenv("REF_PATH", str(ref))
    path = tmp_path / "o.tio"
    im.encode("sam", [str(sam)], str(path), reference=str(ref), embed_reference=True,
              block_reads=200)
    before = _channel_state(path)
    encrypt_per_au(str(path), KEY)
    stream = tmp_path / "e.tis"
    with TransportWriter(str(stream)) as writer:
        write_encrypted_dataset(writer, str(path))
    received = tmp_path / "r.tio"
    read_encrypted_to_file(str(stream), str(received))
    decrypt_per_au_in_place(str(received), KEY)
    after = _channel_state(received)
    assert after[0] == before[0] and after[1] == before[1]


def test_whole_channel_layout(tmp_path):
    from ttio.importers.bam import BamReader
    sam, ref = _make_sam(tmp_path, n=120)
    run = BamReader(sam).to_genomic_run()
    import dataclasses
    run = dataclasses.replace(run, opt_legacy_whole_channel=True)
    out = tmp_path / "o.tio"
    SpectralDataset.write_minimal(str(out), title="t", isa_investigation_id="i",
                                  runs={}, genomic_runs={"g": run})
    with SpectralDataset.open(str(out)) as ds:
        g = ds.genomic_runs["g"]
        assert g.layout != "blocks_v1"
        ds_tags = _tags_ds(g)
        assert io.read_int_attr(ds_tags, "compression") == int(Compression.SAM_TAGS)
        assert "opt_sam_tags" in ds.feature_flags.features
        assert genomic_run_sam_full_md5(g) == sam_full_md5(sam)
