# M101 re-measure against CRAM 3.1 — chr22 slices (2026-10-01)

WORKPLAN item 2 asked for the sizes on record to be measured again after
the REF_DIFF_V2 unmapped-read fix. This run also includes the SAM
optional tags that M101 keeps; the sizes on record had dropped them.

**Code:** branch `m101-sam-optional-tags`, Python importer
(`registry.encode("bam", ..., reference=...)`), blocks_v1 defaults,
external reference (no embedding). Every file round-trips: the full SAM
records, tags included, read back equal to `samtools view` of the
source.

**CRAM:** samtools/htslib 1.22.1,
`samtools view -T <ref> -O cram,version=3.1,<profile>`. The reference
covers every `@SQ` of the BAM header, so CRAM codes against the true
reference and drops MD/NM. That condition matters; see *Correction*
below.

**Machine:** WSL2 Ubuntu 26.04 on this workstation (8 threads, 16 GB).
This is a different machine from the one the earlier numbers came from,
so the timings below are not comparable with them.

## NA12878 WES chr22

GIAB `Garvan_NA12878_HG001_HiSeq_Exome` (bwa 0.7.5a, hg19 names), chr22
slice: 992,974 reads (992,135 mapped + 839 unmapped placed). This is the
slice on record. Reference: b37/hs37d5 chromosome 22 renamed `chr22`;
the other header contigs are N at their header lengths, and no read in
the slice uses them.

| Format | Bytes | vs BAM |
|---|---:|---:|
| BAM | 72,803,408 | 1.000 |
| **TTI-O (M101)** | **44,713,903** | **0.614** |
| CRAM 3.1 normal | 35,801,443 | 0.492 |
| CRAM 3.1 small | 33,322,851 | 0.458 |
| CRAM 3.1 archive | 32,141,988 | 0.441 |

The size on record was 64.1 MB. The REF_DIFF_V2 unmapped-read fix takes
TTI-O to 44.7 MB with the tags included.

| Channel | TTI-O | CRAM 3.1 small |
|---|---:|---:|
| qualities | 26,641,548 | 26,589,391 (QS) |
| read names | 6,699,048 | 2,986,714 (RN) |
| tags | 3,218,274 | 1,259,902 (aux + TL) |
| — of which `XA:Z` | 2,189,033¹ | 932,297 |
| — tag lines | 180,210¹ | 28,016 |
| — MD/NM | ~0 (recomputed) | ~0 (recomputed) |

¹ Phase 0 prototype breakdown of the same reads; the C kernel's total is
within 1% of it.

## HG002 2x250 chr22

GIAB HG002 NIST Illumina 2x250, novoalign against the GRCh38 analysis
set, chr22 slice: 10,633,980 reads. Reference for both formats:
`GRCh38_full_analysis_set_plus_decoy_hla.fa` (TTI-O resolves the chr22
sequence from it).

| Format | Bytes | vs BAM |
|---|---:|---:|
| BAM | 1,635,690,815 | 1.000 |
| **TTI-O (M101)** | **1,020,546,498** | **0.624** |
| CRAM 3.1 normal | 900,680,668 | 0.551 |
| CRAM 3.1 small | 834,585,157 | 0.510 |
| CRAM 3.1 archive | 820,548,458 | 0.502 |

| Channel | TTI-O | CRAM 3.1 small |
|---|---:|---:|
| qualities | 726,874,512 | 712,236,645 (QS) |
| read names | 152,350,325 | 35,233,679 (RN) |
| tags | 20,547,597 | 27,304,364 (aux + TL) |

Read names are 117 MB of the 186 MB gap to CRAM 3.1 small. Qualities
are within 2%, and the tags are 25% smaller than CRAM's.

Wall time on this machine: TTI-O import 98.7 s (peak RSS 4.6 GB); CRAM
3.1 normal 30.2 s, small 48.4 s, archive 73.4 s.

## Correction to the M101 Phase 0 CRAM figure

The first CRAM baseline in Phase 0 (`tools/prototypes/m101_sam_tags/`)
passed a reference that held only chr22. With contigs of the header
missing from the reference, samtools 1.22 switches to `embed_ref=2`: it
embeds a reference built from the reads. Against that reference it
cannot recompute MD/NM, so it stored them verbatim (29.1 + 3.2 MB). The
whole file also grew from 834.6 MB to 1,247.3 MB. The M101 tag
comparison therefore read 59.6 MB for CRAM where the true figure is
27.3 MB, and claimed a 2.9× advantage instead of 25%. The Phase 0
README, `docs/codecs/sam_tags.md`, CHANGELOG, WORKPLAN §97 and HANDOFF
carry the corrected figures.

Whole-genome HG002 is pending.
