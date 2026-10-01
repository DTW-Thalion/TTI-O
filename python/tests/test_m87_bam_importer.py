"""M87 SAM/BAM importer acceptance tests.

The whole file is skipped when ``samtools`` is not on PATH so CI
runners without the binary stay green (HANDOFF.md Gotcha §156).
:class:`~ttio.importers.bam.BamReader` itself remains importable
without samtools per Binding Decision §135 — that property is
exercised by ``test_samtools_missing_error``.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


def _samtools_on_path() -> bool:
    return shutil.which("samtools") is not None


pytestmark = pytest.mark.skipif(
    not _samtools_on_path(),
    reason="samtools not installed; M87 importer tests require it on PATH",
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "genomic"
BAM_PATH = FIXTURE_DIR / "m87_test.bam"
SAM_PATH = FIXTURE_DIR / "m87_test.sam"


# NOTE: the fixture is committed in coordinate-sorted order (chr1 by
# pos ascending, then chr2 by pos ascending, then unmapped reads).
# This is required so `samtools index` succeeds — region filtering
# (tests #12, #13) needs an indexed BAM. Read indexes within these
# arrays therefore mirror the on-disk record order, NOT the r0NN
# numeric labels.
EXPECTED_READ_NAMES = ["r000", "r001", "r002", "r008", "r009",
                       "r003", "r004", "r005", "r006", "r007"]
EXPECTED_POSITIONS = [1000, 1100, 2000, 3000, 4000, 5000, 5100, 0, 0, 0]
EXPECTED_CHROMOSOMES = ["chr1", "chr1", "chr1", "chr1", "chr1",
                        "chr2", "chr2", "*", "*", "*"]
EXPECTED_FLAGS = [99, 147, 0, 16, 0, 99, 147, 4, 77, 141]
EXPECTED_MAPQ = [60, 60, 30, 30, 30, 60, 60, 0, 0, 0]
EXPECTED_CIGARS = ["100M", "100M", "50M50S", "100M", "100M",
                   "100M", "100M", "*", "*", "*"]
EXPECTED_MATE_CHROMS = ["chr1", "chr1", "*", "*", "*",
                        "chr2", "chr2", "*", "*", "*"]
EXPECTED_MATE_POS = [1100, 1000, 0, 0, 0, 5100, 5000, 0, 0, 0]
EXPECTED_TLEN = [200, -200, 0, 0, 0, 200, -200, 0, 0, 0]


# 1
def test_samtools_available():
    """``samtools --version`` succeeds when the binary is on PATH."""
    proc = subprocess.run(
        ["samtools", "--version"], capture_output=True, timeout=10,
    )
    assert proc.returncode == 0, "samtools --version failed"


# 2
def test_read_full_bam():
    """Read m87_test.bam → exactly 10 reads, names r000..r009 in order."""
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert len(run.read_names) == 10
    assert run.read_names == EXPECTED_READ_NAMES


# 3
def test_read_positions():
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert list(run.positions) == EXPECTED_POSITIONS


# 4
def test_read_chromosomes():
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert run.chromosomes == EXPECTED_CHROMOSOMES


# 5
def test_read_flags():
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert list(run.flags) == EXPECTED_FLAGS


# 6
def test_read_mapping_qualities():
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert list(run.mapping_qualities) == EXPECTED_MAPQ


# 7
def test_read_cigars():
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert run.cigars == EXPECTED_CIGARS


# 8
def test_read_sequences_concat():
    """Concatenated SEQ buffer is 720 bytes (per fixture properties).

    "*" reads contribute 0 bytes; everything else contributes its
    SEQ length (100, 100, 100, 100, 100, 0, 10, 10, 100, 100).
    """
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert len(run.sequences) == 720
    # Verify per-read offsets/lengths reconstruct the buffer. Lengths
    # follow the coordinate-sorted on-disk order (r000, r001, r002,
    # r008, r009, r003, r004, r005, r006, r007).
    expected_lengths = [100, 100, 100, 100, 100, 100, 100, 0, 10, 10]
    assert list(run.lengths) == expected_lengths
    assert list(run.offsets) == [
        sum(expected_lengths[:i]) for i in range(10)
    ]


# 9
def test_read_mate_info():
    """Mate fields with RNEXT '=' expanded to RNAME (Binding Decision §131)."""
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert run.mate_chromosomes == EXPECTED_MATE_CHROMS
    assert list(run.mate_positions) == EXPECTED_MATE_POS
    assert list(run.template_lengths) == EXPECTED_TLEN


# 10
def test_read_metadata_from_header():
    """Sample/platform from first @RG; reference_uri from first @SQ."""
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert run.sample_name == "M87_TEST_SAMPLE"
    assert run.platform == "ILLUMINA"
    # First @SQ wins per HANDOFF §2.4.
    assert run.reference_uri == "chr1"


# 11
def test_round_trip_through_writer(tmp_path):
    """BAM → WrittenGenomicRun → .tio → GenomicRun → AlignedRead.

    Iterates reads through the M82 read-side API and confirms the
    parallel-array → AlignedRead materialisation matches the
    importer's view of the source BAM.
    """
    from ttio.importers.bam import BamReader
    from ttio.spectral_dataset import SpectralDataset

    written = BamReader(BAM_PATH).to_genomic_run(name="genomic_0001")

    out = tmp_path / "m87_round_trip.tio"
    SpectralDataset.write_minimal(
        out,
        title="M87 round-trip",
        isa_investigation_id="ISA-M87",
        runs={},
        genomic_runs={"genomic_0001": written},
    )

    ds = SpectralDataset.open(out)
    try:
        gr = ds.genomic_runs["genomic_0001"]
        assert len(gr) == 10
        for i, expected_name in enumerate(EXPECTED_READ_NAMES):
            r = gr[i]
            assert r.read_name == expected_name, f"read {i}"
            assert r.position == EXPECTED_POSITIONS[i], f"read {i}"
            assert r.chromosome == EXPECTED_CHROMOSOMES[i], f"read {i}"
            assert r.cigar == EXPECTED_CIGARS[i], f"read {i}"
            assert r.flags == EXPECTED_FLAGS[i], f"read {i}"
            assert r.mapping_quality == EXPECTED_MAPQ[i], f"read {i}"
            assert r.mate_chromosome == EXPECTED_MATE_CHROMS[i], f"read {i}"
            assert r.mate_position == EXPECTED_MATE_POS[i], f"read {i}"
            assert r.template_length == EXPECTED_TLEN[i], f"read {i}"
        # Spot-check a SEQ round-trip on the mapped paired read r000
        # (index 0 in the coordinate-sorted fixture).
        r0 = gr[0]
        assert r0.sequence == "ACGT" * 25
        assert r0.qualities == b"I" * 100
        # And the empty SEQ for the wholly-unmapped r005 (index 7
        # in coordinate-sorted order).
        r_unmapped = gr[7]
        assert r_unmapped.read_name == "r005"
        assert r_unmapped.sequence == ""
        assert r_unmapped.qualities == b""
    finally:
        ds.close()


# 12
def test_region_filter():
    """region='chr2:5000-5200' → only the two chr2 reads come back."""
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run(region="chr2:5000-5200")
    # r003 + r004 are the chr2 reads in this window.
    assert run.read_names == ["r003", "r004"]
    assert run.chromosomes == ["chr2", "chr2"]


# 13
def test_region_unmapped():
    """region='*' → only unmapped reads (those with no chromosome).

    samtools' '*' selector returns reads whose RNAME is '*' AND whose
    POS is 0; r005, r006, r007 from the fixture all qualify.
    """
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run(region="*")
    assert sorted(run.read_names) == ["r005", "r006", "r007"]
    for chrom in run.chromosomes:
        assert chrom == "*"


# 14
def test_provenance_from_pg():
    """The @PG bwa entry becomes a ProvenanceRecord with the right CL."""
    from ttio.importers.bam import BamReader

    run = BamReader(BAM_PATH).to_genomic_run()
    assert len(run.provenance_records) >= 1
    bwa_records = [p for p in run.provenance_records if p.software == "bwa"]
    assert len(bwa_records) == 1
    bwa = bwa_records[0]
    assert "bwa mem ref.fa reads.fq" in bwa.parameters.get("CL", "")


# 15
def test_sam_input():
    """SamReader on m87_test.sam matches BamReader on m87_test.bam.

    samtools auto-detects format from magic bytes; the two readers
    should produce equal field-by-field outputs (modulo the @PG
    chain, since `samtools view -bS` injects a @PG of its own when
    making the BAM).
    """
    from ttio.importers.bam import BamReader
    from ttio.importers.sam import SamReader

    sam_run = SamReader(SAM_PATH).to_genomic_run()
    bam_run = BamReader(BAM_PATH).to_genomic_run()

    assert sam_run.read_names == bam_run.read_names
    assert list(sam_run.positions) == list(bam_run.positions)
    assert sam_run.chromosomes == bam_run.chromosomes
    assert list(sam_run.flags) == list(bam_run.flags)
    assert list(sam_run.mapping_qualities) == list(bam_run.mapping_qualities)
    assert sam_run.cigars == bam_run.cigars
    assert sam_run.mate_chromosomes == bam_run.mate_chromosomes
    assert list(sam_run.mate_positions) == list(bam_run.mate_positions)
    assert list(sam_run.template_lengths) == list(bam_run.template_lengths)
    assert bytes(sam_run.sequences) == bytes(bam_run.sequences)
    assert bytes(sam_run.qualities) == bytes(bam_run.qualities)
    assert sam_run.sample_name == bam_run.sample_name
    assert sam_run.platform == bam_run.platform
    assert sam_run.reference_uri == bam_run.reference_uri


# 16
def test_samtools_missing_error(monkeypatch):
    """samtools-not-on-PATH must raise with apt/brew/conda guidance.

    Per Binding Decision §135 this is a runtime error at first use,
    NOT an import error. The error text must include actionable
    install guidance for the major OSes.
    """
    from ttio.importers import bam as bam_mod

    monkeypatch.setattr(bam_mod.shutil, "which", lambda _name: None)

    with pytest.raises(bam_mod.SamtoolsNotFoundError) as excinfo:
        bam_mod.BamReader(BAM_PATH).to_genomic_run()

    msg = str(excinfo.value)
    # At least one of the install hints must be present.
    assert any(token in msg for token in ("apt", "brew", "conda"))


# Bonus: the bam_dump CLI emits the canonical-JSON shape from §7.
def test_bam_dump_canonical_json_shape():
    """python -m ttio.importers.bam_dump <bam> emits the §7 schema."""
    from ttio.importers.bam_dump import dump

    payload = dump(str(BAM_PATH))
    expected_keys = {
        "name", "read_count", "sample_name", "platform", "reference_uri",
        "read_names", "positions", "chromosomes", "flags",
        "mapping_qualities", "cigars", "mate_chromosomes", "mate_positions",
        "template_lengths", "sequences_md5", "qualities_md5",
        "provenance_count", "tags",
    }
    assert set(payload.keys()) == expected_keys
    assert payload["read_count"] == 10
    assert payload["sample_name"] == "M87_TEST_SAMPLE"
    assert payload["platform"] == "ILLUMINA"
    # MD5 hex digest of the 720-byte sequences buffer (deterministic).
    assert len(payload["sequences_md5"]) == 32
    assert len(payload["qualities_md5"]) == 32
    # Sorted-keys + indent=2 round-trips to the same dict.
    s = json.dumps(payload, sort_keys=True, indent=2)
    assert json.loads(s) == payload
