# M101 Phase 0: SAM optional tags as a column codec

`tags_codec_proof.py` prototypes the SAM_TAGS column codec. It runs every
block of a BAM through the codec, decodes the block again and compares it
with the input tag text, and it sizes the result against a plain text
channel. CRAM 3.1 is measured separately on the same reads with `samtools`.

## What the prototype codes

Each read's optional fields (SAM columns 12 and up, tab-joined) are cut
into entries `KEY:TYPE:VALUE`. A block stores:

- **a tag-line dictionary.** A tag line is the ordered list of
  `(key, type, kind, arg)` entries of one read. Each read stores one
  index into the dictionary.
- **one column per `(key, type, kind)`**, holding the values in read order.

The entry kinds are:

| Kind | Meaning | Stored |
|---|---|---|
| INT | type `i`, canonical decimal text | integer column |
| TEXT | anything else (`A`, `Z`, `H`, `f`, `B`, or a non-canonical `i`) | verbatim value text |
| DERIVED | `MD:Z` or `NM:i` equal to the value recomputed from the read's sequence, CIGAR and position against the reference | nothing |
| DUP | integer equal to entry `arg` earlier in the same read | nothing |

A non-canonical `i` value is one such as `+5` or `007`. It stays TEXT, so
the decoded text is byte-identical to the input.

Integer columns choose the smallest of three forms:

- constant;
- zigzag varints;
- fixed-width zigzag byte planes, each plane its own substream.

Text columns choose the smallest of three forms:

- constant;
- NUL-joined bytes;
- NAME_TOKENIZED_V2.

Every substream picks the smallest of raw, rANS order-0 and rANS order-1.

MD and NM are recomputed by the samtools `calmd` rules:

- an N on either side is a mismatch;
- `=` in SEQ is a match;
- MD skips soft clips and insertions;
- NM is mismatches + inserted + deleted bases.

## Results, 2026-10-01

**Data:** GIAB HG002 NIST Illumina 2x250, novoalign against the GRCh38
analysis set, chr22 slice (`samtools view <bam> chr22`):

- 10,633,980 reads.
- BAM 1,635,690,815 bytes.
- Reference: chr22 of `GRCh38_full_analysis_set_plus_decoy_hla.fa` (M5 `ac37ec46683600f808cdd41eac1d55cd`).

**Tags present:** PG:Z, AS:i, UQ:i, NM:i, MD:Z, PQ:i, SM:i, AM:i, plus
ZS:Z and NH:i on a few reads.

**Blocks:** blocks_v1 defaults (1M reads or 64 MiB of sequence plus
qualities), 79 blocks.

**Round trip:** byte-exact on every read.

| | Tag bytes | Bytes/read |
|---|---:|---:|
| SAM tag text (with separators) | 776,524,017 | 73.02 |
| Text channel: varint length + text, rANS order-1 | 181,623,170 | 17.08 |
| CRAM 3.1 `small`, aux series + TL | 59,616,292 | 5.61 |
| **Column codec, MD/NM derived** | **20,532,130** | **1.93** |

Per tag, prototype against CRAM 3.1 `small` (from `samtools cram-size -v`):

| Tag | Prototype | CRAM 3.1 small | Note |
|---|---:|---:|---|
| MD:Z | 0 | 29,079,158 | derived on 10,457,612 / 10,457,612 reads |
| NM:i | 0 | 3,232,952 | derived on every read |
| UQ:i | 0 | 5,937,126 | DUP of AS on every read |
| PQ:i | 8,486,542 | 9,062,003 | DUP on 1,871,844 reads |
| AS:i | 6,538,917 | 5,939,784 | |
| SM:i | 2,126,247 | 2,142,639 | |
| AM:i | 1,039,709 | 3,177,497 | DUP of SM on 8,264,860 reads |
| tag lines | 2,287,480 | 984,859 | dictionary + per-read ids |
| PG, ZS, NH | 53,156 | 60,274 | |

**CRAM did not drop MD/NM here.** It stored both verbatim, although every
value matches the recomputation.

**Without a reference** (first 400,000 reads), MD/NM must be stored:

- The codec takes 8.15 B/read, of which MD:Z is 5.1 B/read.
- CRAM's bzip2 MD column is about 2.7 B/read.
- An MD-specific transform is a follow-up for that case.

**Speed:** the pure-Python prototype runs at 30k reads/s encode and
39k reads/s decode per core. That is not a target; the codec belongs in
the C kernel.

## Reproduce

Run under WSL, with ttio installed in a venv and `TTIO_RANS_LIB_PATH`
pointing at a Linux build of `native/`:

```
samtools view -b -o hg002_2x250.chr22.bam \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam chr22
samtools faidx https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa chr22 > GRCh38_analysis_set.chr22.fa
python tools/prototypes/m101_sam_tags/tags_codec_proof.py hg002_2x250.chr22.bam \
  --reference GRCh38_analysis_set.chr22.fa --json proof_full.json
samtools view -C -T GRCh38_analysis_set.chr22.fa --output-fmt-option version=3.1 \
  --output-fmt-option small -o chr22.3.1.small.cram hg002_2x250.chr22.bam
samtools cram-size -v chr22.3.1.small.cram
```
