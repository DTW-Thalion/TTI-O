# M103 Phase 0: a reference-free model for bases

WORKPLAN item 5. Unaligned reads code `sequences` with rANS order-1 at
about 1.94 bits per base, which loses to xz on FASTQ. This prototype sizes
a context-model replacement before any codec is built.

`seq_model_proof.c` reads FASTQ on stdin and reports the ideal code length
(sum of -log2 p) of the A/C/G/T bases under an adaptive model; an arithmetic
coder reaches that to within a fraction of a percent. Nothing is written.

## Model

- Each base is two binary decisions (high bit, then low bit), so three
  counters per context.
- One table per order k: a 16-bit probability with an adaptive rate
  (1/2, 1/3, ... saturating at 1/61.5). Contexts of 2k bits index the
  table directly when they fit, otherwise through a Fibonacci hash, with
  no collision check.
- With two or three orders, a logistic mixer combines the stretched
  predictions; weights are selected by node and by the highest count among
  the orders.
- `--rc`: after a read is coded, the model also trains on its reverse
  complement, so a read from the other strand finds its context. The
  decoder can do the same; it costs time, not bytes.
- Non-ACGT bases reset the context and are left to an exception list.

## Results

bits/base and the ideal size of the bases alone. TTI-O today is rANS
order-1 at about 1.94 bits/base.

### NA12878 WES chr22 (992,974 reads, 95.0 M bases), 2^24 tables

`run_wes.sh`. Shuffled approximates sequencer order; sorted is coordinate
order from the BAM.

| Model | Shuffled | Sorted |
|---|---|---|
| TTI-O today (rANS O1) | ~1.94 (23.0 MB) | |
| xz -9 (whole sequence stream) | 0.70 (8.3 MB) | |
| order 12 | 0.798 | 0.781 |
| order 16 | 0.736 | 0.729 |
| order 20 | 0.797 | 0.790 |
| order 24 | 0.858 | 0.852 |
| --rc order 16 | 0.687 | 0.675 |
| --rc 12 + 20 | 0.571 | 0.552 |
| --rc 11 + 16 + 24 | **0.558 (6.63 MB)** | 0.544 (6.46 MB) |
| --rc 12 + 24 | 0.588 | 0.567 |

### HG002 2x250 chr22 (10.6 M reads, 2.64 G bases, ~65x)

`run_hg002.sh`. Reads in samtools-collate order (hashed by name).

| Model | Tables | bits/base | Bases |
|---|---|---|---|
| TTI-O today (rANS O1) | | ~1.94 | ~640 MB |
| --rc 11 + 16 + 24 | 2^24 | 1.122 | 371 MB |
| --rc 12 + 20 | 2^26 | 0.818 | 270 MB |
| --rc 11 + 16 + 24 | 2^26 | 0.656 | 217 MB |
| --rc 11 + 16 + 22 | 2^28 | 0.453 | 150 MB |
| --rc 11 + 16 + 24 | 2^28 | **0.452** | **149 MB** |

### SEQ_CM kernel in 64 MiB blocks (the blocks_v1 writer's)

`native/tools/seq_cm_fastq` codes each block on its own (model restarted)
and checks the round trip. Default parameters (orders 11/16/24, RC,
automatic table_bits = 24).

| Data | Read order | Blocks | bits/base | Bases | Encode / decode |
|---|---|---|---|---|---|
| WES chr22 | coordinate | 2 | 0.496 | 5.90 MB | 8.6 / 8.3 MB/s |
| WES chr22 | shuffled | 2 | 0.624 | 7.42 MB | 4.8 / 4.6 MB/s |
| HG002 chr22 | coordinate | 40 | **0.348** | 114.9 MB | 7.4 / 6.8 MB/s |
| HG002 chr22 | name hash (collate) | 40 | 1.453 | 480.4 MB | 3.3 / 3.2 MB/s |

Per-block coding is order-sensitive: a block of coordinate-ordered reads
holds one region at full depth, so it beats the whole-run model (0.35
against 0.45); a block of shuffled reads holds ~1.7x of the whole
chromosome, so the long contexts rarely see a read's overlap partner.
Shuffled or sequencer-order FASTQ needs either larger blocks for this
channel or reads grouped by sequence before blocking (with the input
order stored).

### HG002 HiSeq WGS slice (22.5 M reads, 3.33 G bases, ~1x, sequencer order)

`run_wgs.sh`, whole-run prototype. xz -9 on the sequence lines: 769.4 MB
(1.85 bits/base).

| Model | Tables | bits/base |
|---|---|---|
| --rc 11 + 16 + 24 | 2^24 | 1.675 |
| --rc 11 + 16 + 24 | 2^26 | 1.645 |
| SEQ_CM kernel, 64 MiB blocks (50), round trip checked | 2^24 | 1.737 (723.6 MB; 4.3 / 4.1 MB/s) |

At ~1x the reads barely overlap, so little is left to model.

## Findings

1. A mixed order-11/16/24 model with reverse-complement training reaches
   0.56 bits/base on the exome slice and 0.45 on the 65x chr22 slice:
   3.5x and 4.3x smaller than today, and below xz on the exome.
2. With one model over the whole run, read order barely matters
   (shuffled vs sorted differ by 2-3%). Coded per 64 MiB block, as the
   writer does, it matters a great deal: HG002 chr22 is 0.35 bits/base
   in coordinate order and 1.45 in name-hash order (table above).
3. Model memory is the main lever on high-coverage data. On HG002 the
   hashed long orders saturate: 2^24 to 2^26 entries per order cuts 41%,
   2^26 to 2^28 another 31%. At 2^28 (4-byte counters, about 3.2 GB per
   hashed order, 6.5 GB in all) the curve is still falling. The table size
   must be a codec parameter written to the stream, since the decoder
   needs identical tables, and a writer has to choose it from the run size.
4. Order 22 vs 24 makes no difference once the tables are large; the third
   order is worth 0.16 bits/base on HG002 against two orders.

## Open for the codec design

- Counter packing: 12-bit probability + 4-bit state in 2 bytes would
  halve memory at the same table size; check the cost in bits.
- Collision handling (check bits per slot) against plain hashing.
- Table size policy: fixed tiers (2^24 / 2^26 / 2^28) chosen from the
  base count, recorded in the stream header.
- Block policy for unaligned runs (see the per-block table above).
- Arithmetic coder: the kernel uses its own carry-less binary coder
  (docs/codecs/seq_cm.md section 4) rather than `rc_cram`, whose
  multi-symbol interface costs a division per decision.
