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
| --rc 11 + 16 + 24 | 2^28 | 1.590 |
| --rc 12 + 20 | 2^28 | 1.620 |
| --rc 16 | 2^24 | 1.770 |
| SEQ_CM kernel, 64 MiB blocks (50), round trip checked | 2^24 | 1.737 (723.6 MB; 4.3 / 4.1 MB/s) |

At ~1x the reads barely overlap, so little is left to model.

## Grouping reads before blocking

Decision (2026-10-02): unaligned runs reorder their reads so that reads
from the same place in the genome share a block, and store the
permutation that restores the input order. The alternative, larger blocks
for this channel, was rejected.

`native/tools/seq_cm_fastq.c` reorders the whole run before cutting
64 MiB blocks, codes every block with SEQ_CM (orders 11/16/24, RC,
automatic table_bits) and checks every round trip. The permutation is
reported at its ideal size, log2(n!) bits, and as fixed-width indices;
"total" below uses the ideal size. `run_grouping.sh` reproduces the runs.

- `--group`: sort reads by their canonical minimizer (k = 20, least hash
  over both strands), then by offset from it. `--pairs` keeps mates
  together and permutes pairs.
- `--chain`: overlap layout. Each read is indexed by its (w, k = 20)
  window minimizers on both strands, w = ceil(median read length / 8)
  clamped to 8..32 (32 for 250 bp reads, 13 for 100 bp), so a read carries
  about 16 whatever its length. From a seed, walk to the unplaced
  read with the least positive offset, tracking strand and offset, then
  the same leftwards. A neighbour counts only when 2 minimizers agree on
  its offset (`--min-votes`), and minimizers shared by over 256 reads are
  ignored (`--max-occ`). The next chain is seeded from the most recently
  placed read that still has an unplaced neighbour (a stack: depth-first
  with backtracking). A read that chains to nothing is moved to the end of
  the run.

### HG002 2x250 chr22, name-hash (collate) order, 40 blocks

| Order | Bases | + permutation | Total |
|---|---|---|---|
| input order | 1.453 | | 1.453 |
| `--group --pairs` | 0.925 | 0.042 | 0.967 |
| `--group` | 0.568 | 0.088 | 0.656 |
| `--chain`, first version (1 vote, next seed from the last 64 reads) | 0.511 | 0.088 | 0.599 |
| `--chain` (w 32) | **0.407** | 0.088 | **0.495** |
| `--chain --group-w 16` | 0.426 | 0.088 | 0.514 |
| coordinate order (needs the alignment) | 0.348 | | 0.348 |

Fixed-width indices cost 0.096 bits/base here (0.503 total). Grouping
by pairs loses: mates come from opposite ends of a fragment, so keying a
pair on one mate leaves the other among unrelated reads.

`--chain` (w 32): 381,995 chains, 307,158 reads moved to the end (2.9%),
58,597 seeds from input order; peak RSS about 7.3 GB; coding 8.3 / 7.8
MB/s. Indexing takes 33.5 s and grouping 93.7 s in all, down from 60 s and
237.5 s when the index was qsorted and every step sorted its candidates
(now a stable radix sort and a vote hash table; the layout is identical).
A denser window does not help 250 bp reads: w 16 is 0.426 bits/base,
with 11.7 GB peak RSS and 165 s of grouping.

### Where the rest of the gap is

`--pos` takes each read's mapping position (used for this diagnostic
only) and reports how local an order is. "Busiest bins" is the share of a
block's reads in its 20 busiest 100 kb bins.

| Order | Consecutive reads within 1 kb | 100 kb bins per block | Busiest bins |
|---|---|---|---|
| coordinate | 100% | 10.6 | 100% |
| input order (mates adjacent) | 49.9% | 383 | |
| `--group` | 76.9% | 383 | 9.8% |
| `--chain`, 1 vote | 80.4% | 383 | |
| `--chain` (w 32) | 93.1% | 355 | 26.7% |

The chains are locally sound but the blocks are still far from one region
at full depth. On simulated error-free reads (3 Mbp, 30x, half reverse
strand) the same code makes 1 seed from input order and every block is
local, so the remaining breaks come from repeats and sequencing errors,
not from the layout itself. Raising `--max-occ` to 1024 cut the
input-order seeds from 111k to 67k but made blocks less local (busiest
bins 24.2% against 26.1%), because chains then cross between repeat
copies.

### WES chr22 (~100 bp reads), 2 blocks

| Order | Bases | + permutation | Total |
|---|---|---|---|
| shuffled | 0.624 | | 0.624 |
| `--group` | 0.556 | 0.193 | 0.749 |
| `--group --pairs` | 0.598 | 0.091 | 0.689 |
| `--chain`, first version (1 vote, next seed from the last 64 reads) | 0.514 | 0.193 | 0.707 |
| `--chain --group-w 32` | 0.526 | 0.193 | 0.719 |
| `--chain` (w 13 from read length) | **0.503** | 0.193 | 0.696 |
| `--chain --min-votes 1` (w 13) | 0.500 | 0.193 | 0.693 |
| coordinate order | 0.496 | | 0.496 |

At exome depth the permutation costs more than grouping saves, so a
writer must choose per run whether to group. At w = 32 a ~100 bp read
carries about 5 minimizers, 2 agreeing ones are rare, and 21% of reads
end up moved to the end; at w = 13 that is 13%, and the bases come
within 0.007 bits/base of coordinate order. A window of 8 or 16 gives
0.501 and 0.505; 3 votes 0.506. Grouping takes 3.6 s.

### Closing the HG002 gap: what was tried

Each option is added on top of `--chain`; HG002 numbers are the bases.

| Option | HG002 bases | HG002 total | WES shuffled bases |
|---|---|---|---|
| `--chain` | 0.4069 | 0.4949 | 0.5027 |
| `--fill`: after each chain read, its unplaced neighbours on even one minimizer within a read length, by offset | 0.3972 | 0.4853 | **0.4967** |
| `--fill --mates`: seed the next chain from the unplaced mate (paired by read name) of a recent read | 0.3960 | 0.4841 | 0.4967 |
| `--fill --mates --scaffold`: union-find clustering of chains along overlap and mate links (>= 4), clusters capped at 16 MiB, emitted whole | **0.3928** | **0.4808** | 0.4967 |
| coordinate order | 0.348 | | 0.496 |

WES is closed: `--fill` reaches coordinate order. HG002 is not. Two
other attempts lost and stay only as options: `--layout` (a coordinate
per read from a breadth-first search over confident overlaps, then a
sort) is exact on simulated reads and 0.5004 on WES, but on HG002 a few
overlaps through repeats put whole regions into one frame and only 2.2%
of consecutive reads land within 1 kb; joining chains into paths (at most
two neighbours each) instead of clusters gave 0.3944. A first `--mates`
assumed reads 2i and 2i+1 were mates; in this FASTQ that breaks at pair
733, so mates are paired by name.

Diagnostics (`--oracle-sort`, `--dump-order` and `frag.py`-style
analysis of the dumped order against mapping positions) say where the
rest is:

- Order within a block does not matter: sorting each block of the
  `--fill --mates` order by true position codes 0.3958 against 0.3960.
  Which reads share a block is everything.
- In coordinate order a 100 kb bin's reads span 1.11 of the 40 blocks;
  in every grouped order here, 39.5-40. The grouped order is 1.59 M runs
  of genome-contiguous reads (median 1 read, 99th percentile 127 reads,
  none longer than 11.6 kb of genome), and consecutive runs are a median
  6.8 Mbp apart. A 100 kb bin's ~28,000 reads sit in 1,000-1,900
  separate stretches of the output.
- So the local order is good and the long-range order is random. What
  is missing is which ~10 kb run follows which along the chromosome,
  across repeats longer than a read: scaffolding. Overlaps cannot bridge
  those repeats, and mates (fragments of ~500 bp) bridge only the short
  ones.

Grouping with `--scaffold` takes 238 s, most of it the pass that counts
links between chains.

### Open

- The HG002 gap (0.393 against 0.348 for the bases): long-range order
  across repeats, i.e. assembly-grade scaffolding of the ~10 kb runs.
- Grouping speed: 93.7 s for 10.6 M reads with `--chain`, 238 s with
  `--scaffold`, single-threaded; indexing, the walk and the link count can
  all be split across threads.
- Reads moved to the end and the remaining input-order seeds.
- When a writer groups: a rule from coverage and read count.
- The permutation's wire format, and random access by input index (a
  fixed-width column keeps a lookup at one block).

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
- Block policy for unaligned runs: decided, grouped with a stored
  permutation (see "Grouping reads before blocking"); the writer's rule
  for when to group is still open.
- Arithmetic coder: the kernel uses its own carry-less binary coder
  (docs/codecs/seq_cm.md section 4) rather than `rc_cram`, whose
  multi-symbol interface costs a division per decision.
