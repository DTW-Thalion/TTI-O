# SEQ_CM codec (codec id 19)

> **Status:** M103. Reference kernel in C (`native/src/seq_cm.{c,h}`),
> shared by the three SDKs (Python `ttio.codecs.seq_cm`, Java
> `codecs.SeqCm` over JNI, ObjC `TTIOSeqCm`); the `blocks_v1` writers
> code `sequences` with it when a block has no reference. Writers use
> the default parameters with `table_bits` chosen from the base count
> (§2.2), so the same reads always code to the same bytes and a per-AU
> decrypt restore re-encodes byte-identically. This document is binding
> for the wire format.

The codec stores read bases without a reference. It replaces rANS
order-1 for the `sequences` channel of runs with no `REF_DIFF_V2`
layout (FASTQ imports, unaligned reads, runs without a reference),
where rANS order-1 reaches only about 1.94 bits per base.

Each A/C/G/T base is predicted from the bases before it in the same read
by a context-mixing model: up to three context lengths (default 11, 16
and 24 bases), each with its own table of adaptive counters, combined by
a logistic mixer. A long context recognises sequence that an earlier read
covered, so overlap between reads, which rANS order-1 cannot see, turns
into compression. After a read is coded, the model optionally also
trains on its reverse complement, so a read from the other strand finds
its context. Decode repeats every step, so the model is never stored.

All arithmetic is integer, so the stream is the same on every platform.

Phase 0 (`tools/prototypes/m103_seq_model/`): 0.56 bits per base on the
NA12878 WES chr22 reads and 0.45 on HG002 2x250 chr22 with a model
trained over the whole run.

The codec is **context-aware**: decode takes the reads' lengths, which
the run already stores, as side input.

---

## 1. Input

The encoder input is `n_reads` reads, their bases concatenated into one
byte string `S` of `n_bases` bytes, and the read lengths `L[0..n_reads)`.
Every byte value is allowed. The bytes `A C G T` (0x41 0x43 0x47 0x54)
are **bases**; every other byte (`N`, IUPAC codes, lower case, `NUL`) is
an **exception**.

## 2. Blob layout

All integers are little-endian.

| Offset | Size | Field |
|---:|---:|---|
| 0 | 4 | magic `SQC1` |
| 4 | 1 | version, `0x01` |
| 5 | 1 | flags: bit 0 `RC` (reverse-complement training); other bits 0 |
| 6 | 1 | `n_orders`, 1 to 3 |
| 7 | 3 | `orders[3]`: context lengths in bases, strictly rising, each 1 to 31; unused entries 0 |
| 10 | 1 | `table_bits`, 10 to 30 |
| 11 | 1 | reserved, 0 |
| 12 | 2 | `limit`, counter adaptation limit, 1 to 1023 |
| 14 | 2 | `lr`, mixer learning rate, 1 to 4096 |
| 16 | 8 | `n_reads` |
| 24 | 8 | `n_bases` |
| 32 | 8 | `n_runs`, exception runs |
| 40 | 8 | `exc_len`, bytes in the exception section |
| 48 | `exc_len` | exception section (§2.1) |
| 48 + `exc_len` | rest | arithmetic-coder stream (§4) |

A decoder rejects the blob when the magic, version, a reserved bit or
byte, or any parameter is out of range, or when `n_reads` and `n_bases`
disagree with the lengths it was given.

### 2.1 Exception runs

A run is a maximal stretch of equal exception bytes in `S` (it may cross
read boundaries). Each run is three fields, in order of position:

- `gap` (varint): its start minus the end of the previous run (0 for the
  first run measured from position 0);
- `length` (varint, at least 1);
- `byte` (1 byte, never a base).

Varints are unsigned LEB128. Runs must lie within `n_bases`, and the
section must hold exactly `n_runs` runs and end at `exc_len`.

### 2.2 Choice of `table_bits`

The writer chooses `table_bits`; the decoder only reads it. With
`table_bits` 0 at the API, the reference encoder takes
`ceil(log2(n_bases)) + 1`, clamped to 16 to 24. The model holds
`16 << table_bits` bytes per hashed order (§3.1), on encode and on
decode.

## 3. Model

### 3.1 Counters and contexts

The context of a base is `h`, the previous bases of the same read as
2-bit codes (`A`=0, `C`=1, `G`=2, `T`=3), most recent in the low bits.
`h` is 0 at the start of every read and after every exception byte.

Every order `k` has a table of **buckets** of 16 counters. A bucket is
chosen by the context's older `k - 1` bases,
`key = (h >> 2) & (4^(k-1) - 1)`, and holds the four contexts that differ
only in the newest base `h & 3`: that context's three counters sit at
`bucket * 16 + (h & 3) * 4`, and the fourth slot of each group is unused.
The bucket index is

- direct, when `2k <= table_bits`: `key` (the table has `4^(k-1)`
  buckets);
- hashed otherwise: `(key * 0x9E3779B97F4A7C15) >> (64 - (table_bits - 2))`,
  a 64-bit unsigned multiply keeping the low 64 bits (the table has
  `2^(table_bits - 2)` buckets).

A table of `table_bits` holds `16 << table_bits` bytes. Because the next
base's bucket depends only on bases already known, a coder can fetch it
while the current base is coded; the layout is part of the format only
through which counters contexts share.

A base `c` is coded as two binary decisions: `hi = c >> 1` at node 0,
then `lo = c & 1` at node `1 + hi`, each predicted by the context's
counter for that node.

A counter holds a probability `p` (of the bit being 1, in 1/65536, start
32768) and a count `n` (start 0). After the bit is known, with
`r = 65536 / (n + 2)` (integer division):

- bit 1: `p += ((65535 - p) * r) >> 16`;
- bit 0: `p -= (p * r) >> 16`;
- then `n += 1` while `n < limit`.

### 3.2 Logistic domain

`squash(d)` maps `d` in [-2047, 2047] (values outside are clamped) to a
12-bit probability:

```
T = {1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,1546,2047,2549,
     2994,3348,3607,3785,3901,3975,4022,4050,4068,4079,4085,4089,4092,
     4093,4094}
x = d + 2048;  i = x >> 7;  w = x & 127
squash(d) = (T[i] * (128 - w) + T[i + 1] * w + 64) >> 7
```

`stretch(p)` for `p` in [0, 4095] is the smallest `d` in [-2047, 2047]
with `squash(d) >= p`, or 2047 when there is none.

### 3.3 Mixer

Inputs: `s_i = stretch(p_i >> 4)` for each order `i`, plus a bias input
`s_n = 64`. The weight set is chosen by node and confidence:
`set = node * 64 + min(max_i n_i, 63)`. There are 192 sets of
`n_orders + 1` weights (16.16 fixed point), each order's weight starting
at `65536 / n_orders` (integer division) and the bias weight at 0.

```
dot = sum_i w[set][i] * s_i                (64-bit)
p12 = clamp(squash(floor(dot / 65536)), 1, 4095)
```

The coder codes the bit with `P(1) = p12 * 16 / 65536` (§4). After the
bit is known:

```
err = ((bit << 12) - p12) * lr
w[set][i] += floor(err * s_i / 65536)      for every input i
```

then every order's counter updates (§3.1). Divisions marked `floor`
round towards minus infinity.

### 3.4 Coding order and reverse-complement training

Reads are coded in order. Within a read, each base is coded at its
context, then `h = (h << 2) | c`; an exception byte codes nothing and
sets `h = 0`.

With the `RC` flag, once a read is coded the model trains on its reverse
complement: the read's bytes from last to first, each base replaced by
`3 - c`, exceptions resetting `h`. Training performs the same predict
and update steps as coding (mixer included) but codes nothing.

## 4. Arithmetic coder

A binary carry-less coder with 32-bit bounds `x1 = 0`, `x2 = 0xFFFFFFFF`.
To code `bit` with `P(1) = q / 65536` (`q = p12 * 16`):

```
xmid = x1 + (((x2 - x1) * q) >> 16)        (64-bit product)
bit ? x2 = xmid : x1 = xmid + 1
while ((x1 ^ x2) & 0xFF000000) == 0:
    emit x2 >> 24;  x1 <<= 8;  x2 = (x2 << 8) | 0xFF
```

At the end the encoder emits the four bytes of `x1`, high first. The
decoder starts with `x` = the first four stream bytes (high first),
decides `bit = (x <= xmid)`, updates `x1` and `x2` the same way and
shifts in one byte per emitted byte; bytes past the end of the stream
read as 0.

## 5. Open for M103

- Read order: blocks_v1 codes each block on its own, so the model
  restarts every block and the result depends on which reads share a
  block. Unaligned runs will be able to reorder their reads before
  blocking and store the permutation (opt-in; M103 Phase 0 README,
  "Grouping reads before blocking").
- Throughput: about 8 MB/s per core on grouped or coordinate-ordered
  reads, about 3 MB/s on shuffled reads, in the reference kernel.
- Model memory: `16 << table_bits` bytes per hashed order on encode and
  decode (about 0.5 GB at the automatic maximum of 24), per block coded
  at once.
