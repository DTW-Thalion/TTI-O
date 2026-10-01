# SAM_TAGS codec (codec id 18)

> **Status:** current (M101). Default and only codec for the genomic
> `tags` channel. Reference kernel in C (`native/src/sam_tags.{c,h}`);
> language wrappers in Python (ctypes — `ttio.codecs.sam_tags`), Java
> (JNI — `SamTags.java`) and Objective-C (direct link —
> `TTIOSamTags.{h,m}`). All three call the same C entry points, so a
> block encodes to the same bytes whichever SDK writes it.

The codec stores the SAM optional fields of a run of reads: everything
after SAM column 11, tab-joined, exactly as `samtools view` prints it
(`NM:i:1\tMD:Z:130C117\tRG:Z:rg1`). It is lossless on that text: decode
returns, for every read, the exact bytes the encoder was given.

It is a CRAM-style column codec. Each read's tags become one entry in a
per-blob **tag-line dictionary** plus one value per tag, and the values
of one tag key go to one **column**. Two kinds of value store nothing:

- `MD:Z` and `NM:i` equal to the values recomputed from the read's
  sequence, CIGAR and position against the reference (`DERIVED`);
- an integer equal to an integer earlier in the same read (`DUP`), such
  as novoalign's `UQ` that repeats `AS`.

The codec is **context-aware**, like REF_DIFF_V2 and MATE_INLINE_V2:
encode and decode take the reads' sequences, CIGARs, positions and
chromosome ids, and the reference bases, as side input. The blob does
not store them.

Phase 0 (`tools/prototypes/m101_sam_tags/`): 1.93 bytes per read on the
GIAB HG002 2x250 chr22 slice against 2.57 for CRAM 3.1 `small`; on
the bwa-aligned NA12878 WES chr22 slice CRAM is ahead (1.27 against
3.24 B/read) because of the long `XA:Z` strings
(`docs/benchmarks/2026-10-01-m101-cram-remeasure.md`).

---

## 1. Tag text

### 1.1 Parsing

The encoder input for read *i* is a byte string `T_i` (empty when the
read has no optional fields). `T_i` must not contain `NUL`. It is split
on `TAB` into fields. A field is **well-formed** when

- it is at least 5 bytes long,
- bytes 2 and 4 are `:`,
- byte 0 is in `[A-Za-z]` and byte 1 in `[A-Za-z0-9]` (the tag key),
- byte 3 is one of `A i f Z H B` (the type).

The value is everything after byte 4 (it may contain `:`, and may be
empty). If `T_i` is non-empty and any field is not well-formed (an empty
field from a doubled or trailing `TAB` included), the read is coded
**verbatim** (§1.4).

### 1.2 Entry kinds

Each field of a read becomes one entry `(key, type, kind, arg)`:

| Kind | Id | Condition | Stored |
|---|---:|---|---|
| `INT` | 0 | type `i` and the value is canonical: `0` or `-?[1-9][0-9]*`, within int64 | integer column `(key, i, INT)` |
| `TEXT` | 1 | any other well-formed field | text column `(key, type, TEXT)` |
| `DERIVED` | 2 | §1.3 | nothing |
| `DUP` | 3 | an `INT` field whose value equals the integer value of an earlier entry `j` of the same read, `j < 256`; `arg = j`, the smallest such `j` | nothing |
| `VERBATIM` | 4 | the whole read, §1.4 | text column `(0x00 0x00, 0x00, VERBATIM)` |

`arg` is 0 for every kind but `DUP`. An entry **has an integer value**
when it is `INT`, `DUP`, or a `DERIVED` `NM`. `DERIVED` is decided before
`DUP`.

### 1.3 MD and NM recomputation

A read is **derivable** when the encoder (and the decoder) is given its
reference, its POS is ≥ 1, its CIGAR is not `*` and parses as
`([0-9]+[MIDNSHP=X])+`, its SEQ is not empty and not `*`, the CIGAR's
query length (the sum of `M I S = X` lengths) equals the SEQ length, and
every reference position the CIGAR touches lies inside the reference.

For a derivable read the codec computes `(MD, NM)` with the samtools
`calmd` rules. Let `up(b)` map `a..z` to `A..Z` and leave other bytes
alone. Walk the CIGAR with `r = POS - 1` and `q = 0`:

- `M`, `=`, `X` of length *n*: for each of the *n* positions, the base
  **matches** when `SEQ[q] == '='`, or when `up(SEQ[q]) == up(REF[r])`
  and neither is `N`. A match adds 1 to the current run; a mismatch
  appends the run's decimal digits and `up(REF[r])` to MD, resets the
  run to 0, and adds 1 to NM. Advance `q` and `r`.
- `I` of length *n*: NM += *n*; `q` += *n*.
- `D` of length *n*: append the run's digits, `^`, and
  `up(REF[r .. r+n-1])` to MD; the run becomes 0; NM += *n*; `r` += *n*.
- `N` of length *n*: `r` += *n*.
- `S` of length *n*: `q` += *n*.
- `H`, `P`: nothing.

At the end, append the run's digits. (So MD always starts and ends with
digits, and two adjacent mismatches are separated by `0`.)

A field `MD:Z:v` is `DERIVED` when the read is derivable and `v` equals
the recomputed MD. A field `NM:i:v` is `DERIVED` when the read is
derivable and `v` is the decimal text of the recomputed NM. Anything
else falls to `INT` / `TEXT` / `DUP` as usual: a value the aligner
computed differently simply costs its bytes.

### 1.4 Verbatim reads

A verbatim read's tag line has a single entry `(0x00 0x00, 0x00,
VERBATIM, 0)` and its whole `T_i` goes to the verbatim text column.

---

## 2. Columns

Columns are keyed by `(key, type, kind)` and numbered in order of first
appearance. Each holds its values in read order.

### 2.1 Substreams

A substream is a byte string `S` framed as

```
mode     u8         0 = raw, 1 = rANS order-0 with a sparse table
raw_len  uvarint    len(S), at most 2^32 - 1
body_len uvarint
body     body_len bytes
```

Mode 0's body is `S`. Mode 1's body is the RANS_ORDER0 (codec 4) encoding
of `S` with its fixed 1024-byte frequency table listed sparsely:

```
n_syms   uvarint                      1..256
entries  n_syms × { sym u8, freq uvarint }   ascending sym, freq 1..4096
payload  the rest of the body         codec 4's payload (final state + renorm bytes)
```

A decoder rebuilds the codec-4 blob (order byte 0, `orig_len = raw_len`,
`payload_len`, the 256 frequencies with zeros for absent symbols, the
payload) and decodes it with RANS_ORDER0. The encoder uses mode 1 when
`S` is non-empty and the mode-1 body is shorter than `S`; otherwise raw.

The sparse table matters because blobs are per block: a codec-4 stream
spends 1024 bytes on its table whatever its length. rANS order-1 is not
used: in Phase 0 it saved 3.3% of the tag bytes; adding it is a codec
version bump.

### 2.2 Integer columns

Values are int64, written through zigzag `u = (v << 1) ^ (v >> 63)`.

```
transform u8
  2 CONST   : uvarint u                         (every value equal)
  0 VARINT  : substream of uvarint(u) per value
  1 PLANES  : width u8 ∈ {1,2,4,8}, then `width` substreams;
              substream k holds byte k (little-endian) of each u
```

`width` is the smallest of 1, 2, 4, 8 bytes that holds the largest `u`.
The encoder emits CONST when all values are equal; otherwise it builds
VARINT and PLANES and keeps the shorter, VARINT on a tie.

### 2.3 Text columns

```
transform u8
  2 CONST   : uvarint len, bytes                (every value equal)
  0 JOINED  : substream of the values joined by NUL (no trailing NUL)
  1 NAMETOK : uvarint len, a NAME_TOKENIZED_V2 blob (codec 15) of the values
```

NAMETOK is a candidate only when every value is non-empty 7-bit ASCII
and the NAME_TOKENIZED_V2 encoder accepts the list. The encoder emits
CONST when all values are equal; otherwise it keeps the shortest of the
candidates, JOINED on a tie.

---

## 3. Wire format

All integers are little-endian; `uvarint` is LEB128.

```
magic       4 bytes   "STG1"
version     u8        1
n_reads     uvarint
n_lines     uvarint
lines       n_lines × { n_entries uvarint,
                        n_entries × { key[2], type u8, kind u8, arg u8 } }
line_ids    substream of n_reads uvarints, each < n_lines
n_columns   uvarint
columns     n_columns × { key[2], type u8, kind u8, n_values uvarint,
                          integer column (kind INT) | text column (TEXT, VERBATIM) }
```

Lines are numbered in order of first appearance; a read with no tags
uses a line with zero entries. Column `n_values` equals the number of
entries of that `(key, type, kind)` across all reads. The encoder is
deterministic: the same reads and context produce the same bytes.

### 3.1 Decode

For each read, take its line; for each entry emit `key:type:value`,
joined by `TAB`:

- `INT`, `TEXT`, `VERBATIM` take the next value of their column (a
  verbatim read emits the value alone);
- `DERIVED` recomputes MD or NM by §1.3 (the decoder needs the
  sequences, CIGARs, positions and reference, so the `sequences` channel
  is decoded first);
- `DUP` repeats the value text of entry `arg`.

The decoder rejects (`TTIO_RANS_ERR_CORRUPT`) a bad magic or version, a
`line_id ≥ n_lines`, a `DUP` whose `arg` is not an earlier entry with an
integer value, a `DERIVED` entry on a read that is not derivable, a
column that runs out or keeps values left over, and trailing bytes.

---

## 4. Public API

### C

```c
#include "ttio_rans.h"

typedef struct {
    uint64_t        n_reads;
    const uint8_t  *sequences;      /* concatenated SEQ bytes          */
    const uint64_t *seq_offsets;    /* n_reads + 1                     */
    const uint8_t  *cigars;         /* concatenated CIGAR text         */
    const uint64_t *cigar_offsets;  /* n_reads + 1                     */
    const int64_t  *positions;      /* 1-based POS, 0 = unmapped       */
    const uint16_t *chrom_ids;      /* 0xFFFF = none                   */
    const uint8_t * const *refs;    /* n_refs chromosome sequences, or NULL entries */
    const uint64_t *ref_lengths;
    uint32_t        n_refs;         /* 0 disables MD/NM derivation     */
} ttio_sam_tags_ctx;

int ttio_sam_tags_encode(const ttio_sam_tags_ctx *ctx,
                         const uint8_t *tags, const uint64_t *tag_offsets,
                         uint8_t **out, size_t *out_len);
int ttio_sam_tags_decode(const ttio_sam_tags_ctx *ctx,
                         const uint8_t *encoded, size_t encoded_len,
                         uint8_t **out_tags, uint64_t *out_tag_offsets);
void ttio_sam_tags_free(void *p);
```

`out` and `out_tags` are allocated by the library and released with
`ttio_sam_tags_free`; `out_tag_offsets` is caller-allocated with
`n_reads + 1` entries.

---

## 5. Channel routing

The codec writes the genomic `tags` channel
(`signal_channels/tags`, `@compression = 18`), one blob per blocks_v1
block, and one blob for the whole run under the legacy whole-channel
layout. See `docs/genomic-runs.md` and format-spec §10.12.
