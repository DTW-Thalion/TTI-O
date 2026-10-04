# Genomic Sequencing Data Model

TTI-O models aligned sequencing reads as a **run-and-element hierarchy**
that runs in parallel with the spectrum-based modalities (mass
spectrometry, NMR, vibrational, UV/Vis). Where the MS side has
`MassSpectrum` elements grouped under an `AcquisitionRun` and indexed by
a `SpectrumIndex`, the genomic side has `AlignedRead` elements grouped
under a `GenomicRun` and indexed by a `GenomicIndex`. The two hierarchies
share the same storage primitives (`signal_channels/`, an eagerly-loaded
parallel-array index, per-run provenance), the same `StorageProvider`
protocol (see [`providers.md`](providers.md)), and the same on-disk
container (`/study/...`, see [`format-spec.md`](format-spec.md)).

This document is the data-model walkthrough for that genomic hierarchy:
the value/container classes, the on-disk layout under
`/study/genomic_runs/<name>/`, the selective-access query APIs, the codec
applied to each channel, and the BAM/CRAM/SAM import path. It is written
against the Python reference implementation
(`python/src/ttio/genomic_run.py`, `genomic_index.py`, `aligned_read.py`,
`written_genomic_run.py`, `_dataset_write_genomic.py`, and the
`genomic/` package — `stream_writer.py`, `_blocks.py`,
`reference_resolver.py`); the Objective-C and Java SDKs expose the same
classes and produce byte-exact files.

## 1. Overview — the run-and-element hierarchy

A genomic dataset is held by a `SpectralDataset` (the same top-level
container used by every modality). Its `genomic_runs` attribute is a
`dict[str, GenomicRun]` mapping run name → lazy run view, populated from
`/study/genomic_runs/` at open time. Each `GenomicRun` is a lazy,
sequence-protocol view over one run: `len(gr)` returns the read count and
`gr[i]` materialises the i-th read as an `AlignedRead` value object. The
parallel `GenomicIndex` (per-read scalar arrays) is loaded eagerly so
that counting and region/flag filtering are cheap; the heavy per-base
channels (`sequences`, `qualities`) and the variable-length channels
(`cigars`, `read_names`, `mate_info`) stay lazy on disk and are decoded
on first access.

Modality analogues (for readers familiar with the MS side):

| Genomic class        | MS analogue        | Role                          |
|----------------------|--------------------|-------------------------------|
| `GenomicRun`         | `AcquisitionRun`   | Lazy element-collection view  |
| `AlignedRead`        | `MassSpectrum`     | Per-element frozen value      |
| `GenomicIndex`       | `SpectrumIndex`    | Eager parallel scalar arrays  |
| `WrittenGenomicRun`  | `WrittenRun`       | Write-side data container     |

The model is deliberately shaped after **SAM/BAM**: an `AlignedRead`
carries the SAM core fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
PNEXT, TLEN, SEQ, QUAL), and the importers (§6) map a SAM record's
columns 1–11 directly onto those attributes.

## 2. Value and container classes

### 2.1 `AlignedRead` — the per-read value object

`ttio.aligned_read.AlignedRead` is a frozen, slotted dataclass — one
aligned sequencing read, materialised on demand by
`GenomicRun.__getitem__`. Fields (all populated from the index + signal
channels):

| Attribute          | Type    | Source                                    |
|--------------------|---------|-------------------------------------------|
| `read_name`        | `str`   | `read_names` channel (NAME_TOKENIZED_V2)  |
| `chromosome`       | `str`   | `genomic_index/` chromosome table         |
| `position`         | `int`   | `genomic_index/positions` (0-based)       |
| `mapping_quality`  | `int`   | `genomic_index/mapping_qualities`         |
| `cigar`            | `str`   | `cigars` channel                          |
| `sequence`         | `str`   | `sequences` channel slice, ASCII-decoded  |
| `qualities`        | `bytes` | `qualities` channel slice (Phred bytes)   |
| `flags`            | `int`   | `genomic_index/flags` (SAM flags)         |
| `mate_chromosome`  | `str`   | `mate_info` (MATE_INLINE_V2)              |
| `mate_position`    | `int`   | `mate_info` (MATE_INLINE_V2)              |
| `template_length`  | `int`   | `mate_info` (MATE_INLINE_V2)              |

It exposes SAM-flag convenience properties that decode `self.flags`:
`is_mapped` (flag `0x4` cleared), `is_paired` (`0x1`), `is_reverse`
(`0x10`), `is_secondary` (`0x100`), `is_supplementary` (`0x800`), plus
`read_length` (`len(self.sequence)`; `0` when the source record stored
`SEQ=*`).

Note on README naming: the README summarises mate-pair info as one
field. In code it is three discrete attributes — `mate_chromosome`,
`mate_position`, `template_length` — decoded together from the
`mate_info` channel.

### 2.2 `GenomicRun` — the lazy run view

`ttio.genomic_run.GenomicRun` is a slotted dataclass holding the run-level
attributes plus a loaded `GenomicIndex` and an open `StorageGroup` handle.
Run-level attributes (eagerly read from the run group at open time):
`name`, `acquisition_mode` (`AcquisitionMode`), `modality`,
`reference_uri`, `platform`, `sample_name`, and `channel_names` (the list
of `signal_channels/` child names, kept for introspection — it is not
consulted by `__getitem__`).

Public surface:

- `len(gr)` → read count (delegates to `index.count`).
- `gr[i]` / iteration → materialise the i-th `AlignedRead`. Negative
  indices wrap; out-of-range raises `IndexError`.
- `gr.layout` → `"blocks_v1"` or `"whole"` (the legacy whole-channel
  layout); `gr.block_count` → number of blocks (`1` for whole-channel
  files).
- `gr.iter_reads(start=0, stop=None, *, threads=None)` → yield reads
  `[start, stop)` in order. Under `blocks_v1` the next few blocks
  decode ahead on a thread pool; `threads=1` keeps the serial
  one-block-at-a-time path.
- `gr.for_each_block(fn, start=0, stop=None, *, threads=None)` → call
  `fn(view, view_start, first_read, n_reads)` once per decoded block,
  from several threads at once and in no particular order. `view` is a
  run over one block's reads, indexed from 0 and valid only for the
  duration of the call. Java: `GenomicRun.iterBlocks`; Objective-C:
  `-iterBlocksFrom:to:threads:error:usingBlock:`.
- `gr.reads_in_region(chromosome, start, end)` → `list[AlignedRead]` for
  reads whose **mapping start position** falls in `[start, end)` (see
  §4 for the overlap caveat).
- `gr.provenance_chain()` → per-run provenance records (from
  `<run>/provenance/steps`), `[]` when absent.
- `GenomicRun.open(group, name, *, references_group=None, bulk_read=True)`
  → classmethod factory. The caller resolves the child group first;
  `open` eagerly loads the `GenomicIndex` and run attributes and lists
  the signal-channel names, leaving the channel datasets closed until
  first access. `references_group` is the `/study/references` group
  threaded through to REF_DIFF_V2 decode (§5); `bulk_read=False` keeps
  per-record hyperslab reads for remote fsspec-backed files instead of
  whole-channel bulk reads.

Internally, `GenomicRun` decodes each channel lazily and caches the
decoded result on the instance: byte channels are decoded whole on first
slice (codec output is not sliceable), and `cigars` / `read_names` /
`mate_info` are decoded once into per-read lists/arrays. A run-derived
`CodecContext` (read lengths, reverse-complement flags from `flags & 16`,
positions, chromosomes, encounter-order chrom ids, and the
`ReferenceResolver`) is built once and shared across all registry decode
calls.

### 2.3 `GenomicIndex` — parallel per-read scalars

`ttio.genomic_index.GenomicIndex` is a slotted dataclass of parallel
arrays, all length == read count, loaded eagerly when a `GenomicRun`
opens. It is the genomic analogue of `SpectrumIndex` and the backing
store for all selective-access queries.

| Field                | Type / dtype          | Meaning                              |
|----------------------|-----------------------|--------------------------------------|
| `offsets`            | `np.ndarray` uint64   | byte offset of each read into `sequences`/`qualities` |
| `lengths`            | `np.ndarray` uint32   | read length in bases                 |
| `chromosomes`        | `list[str]`           | reference name per read (`"chr1"`, `"*"` …) |
| `positions`          | `np.ndarray` int64    | 0-based mapping position             |
| `mapping_qualities`  | `np.ndarray` uint8    | Phred-scaled MAPQ                    |
| `flags`              | `np.ndarray` uint32   | SAM flags (uint32 to allow extended bits) |
| `chromosome_ids`     | `np.ndarray` uint16 \| None | interned chrom id per read (disk-loaded) |
| `chromosome_names`   | `list[str]` \| None   | unique names in id order             |

`offsets` is **not stored on disk** in v1.0+ files: `GenomicIndex.read`
derives it as `cumsum(lengths)` (always uint64 to avoid the >4 GB uint32
overflow cliff on deep WGS). `count` is a property returning
`offsets.shape[0]`.

On disk the chromosomes are stored as a uint16 `chromosome_ids` column
plus a compound `chromosome_names` lookup table (encounter-order id
assignment, first occurrence of a name gets the next id; >65 535 unique
names overflows uint16 and raises). When loaded from disk, the interned
`chromosome_ids` / `chromosome_names` enable a vectorised id comparison
in `indices_for_region` instead of an O(N) Python name scan;
in-memory-constructed indexes (no interned table) fall back to scanning
the names list.

Query methods are covered in §4. The disk round-trip is
`GenomicIndex.read(idx_group)` / `GenomicIndex.write(idx_group)`.

### 2.4 `WrittenGenomicRun` — the write-side container

`ttio.written_genomic_run.WrittenGenomicRun` is the data container the
caller hands to `SpectralDataset.write_minimal(..., genomic_runs={name:
wgr})` (or via the unified `runs=` dict — `write_minimal` splits
`WrittenGenomicRun` entries out automatically). It is the explicit,
columnar write-side form: every per-read field is supplied as a numpy
array or list.

Required fields: `acquisition_mode` (the `AcquisitionMode.GENOMIC_WGS`
= 7 or `GENOMIC_WES` = 8 integer value), `reference_uri`, `platform`,
`sample_name`; the per-read arrays `positions` (int64),
`mapping_qualities` (uint8), `flags` (uint32); the concatenated
per-base channels `sequences` (uint8, one ASCII byte per base) and
`qualities` (uint8 Phred); the per-read `offsets` (uint64) and `lengths`
(uint32); the variable-length lists `cigars` and `read_names`; the
mate-pair fields `mate_chromosomes` (list[str]), `mate_positions` (int64,
`-1` if unpaired), `template_lengths` (int32, `0` if unpaired); and
`chromosomes` (list[str], for the index).

Optional fields: `provenance_records`, `signal_compression` (`"gzip"` →
ZLIB, the default; `"none"` → NONE), `signal_codec_overrides`
(`dict[str, Compression]` per-channel codec opt-in for `sequences`,
`qualities`, `cigars`), `embed_reference` (default `False`),
`reference_chrom_seqs` (`chrom_name → uppercase ACGTN bytes`, required
for the REF_DIFF_V2 path), `external_reference_path`,
`opt_disable_qualities_v5` (keep the qualities codec on the V4
strategy set), `opt_legacy_whole_channel` (write the legacy
whole-channel layout instead of `blocks_v1`, §3), `bulk_v2_blobs`
(a `BulkV2Blobs` carrying verbatim v2 codec blobs from the transport
bulk-mode receiver, which bypass the codec encode step),
`read_role` (M97: persisted as the `@read_role` run attribute when
set — recognised values `hifi`, `ont_ul`, `hic_r1`, `hic_r2`,
`parental_maternal`, `parental_paternal`, `illumina_wgs`; other
strings are stored unchecked), and `ref_diff_slice_bytes` (M97: the
REF_DIFF_V2 slice byte budget, `docs/codecs/ref_diff_v2.md` §3;
0 keeps the fixed 10,000-read rule).

Since M97 a `QUALITY_BINNED` override on the `qualities` channel is
rejected when `platform` names a long-read instrument (case-
insensitive `hifi`, `pacbio`, `nanopore`, or the whole token `ont`)
— the fixed Illumina-8 bin table does not fit those quality
distributions. The predicate is
`ttio.codecs.quality.binned_allowed_for_platform` /
`Quality.binnedAllowedForPlatform` /
`TTIOQualityBinnedAllowedForPlatform`.

### 2.5 `GenomicStreamWriter` — the streaming write path

`ttio.genomic.stream_writer.GenomicStreamWriter` appends reads to one
genomic run of an open-for-write dataset without holding the run in
memory: `append(read)` / `append_batch(written_genomic_run)` accumulate
a pending block, `flush()` closes it out, and `close()` finalises the
run. Blocks are cut at `block_reads` / `block_bytes` bounds and at
every chromosome change. With `threads` > 1 (the `threads` argument,
else `TTIO_THREADS`, else cores minus 2) completed blocks are
encoded on a pool and written in order by the caller's thread, bounded
by `memory_budget_bytes`; the file is byte-for-byte the one thread's.
`SpectralDataset.write_minimal` routes every genomic run through this
writer, so the `blocks_v1` layout (§3) is the default for batch and
streaming writes alike; `opt_legacy_whole_channel=True` restores the
whole-channel layout.

## 3. On-disk storage layout

A genomic run is written under `/study/genomic_runs/<name>/`. The run
group carries attributes `acquisition_mode` (int), `modality`
(`"genomic_sequencing"`), `spectrum_class` (`5`), `reference_uri`,
`platform`, `sample_name`, `read_count`, and `base_count`, plus
`read_role` (M97) when the writer was given one — absent otherwise,
and readers surface absence as `None`/`null`/`nil`. The
`genomic_runs` parent group carries a `_run_names` attribute
(comma-joined run names) used to enumerate runs at open time.

Two layouts exist, discriminated by the run group's `@layout`
attribute (`format-spec.md` §10.12):

- **`blocks_v1`** — the default for every write path
  (`SpectralDataset.write_minimal` and `GenomicStreamWriter` alike).
  The run is a sequence of independently coded blocks; a block index
  records where each block's blob sits in each channel dataset.
- **`blocks_v1_grouped`** (M103) — `blocks_v1` with an unaligned
  run's reads grouped by sequence before blocking and
  `genomic_index/input_index` holding the permutation. Opt-in
  (`group_reads` on the stream writer and importers, `opt_group_reads`
  on a written run). Readers present the run in input order; block
  iteration stays in stored order and maps rows through
  `run.input_index`. See format-spec §10.12.7.
- **Whole-channel** — the legacy layout (no `@layout` attribute); one
  codec blob per channel covering the whole run. Written only when
  `opt_legacy_whole_channel=True`; readers keep full support.

### 3.1 `blocks_v1` (default)

```
/study/genomic_runs/<name>/          @layout="blocks_v1", @block_policy,
│                                    @read_count, @base_count
├── blocks/
│   └── index              compound  (one row per block: read_start, n_reads,
│                                     base_start, n_bases, per-channel byte
│                                     ranges, per-channel codec ids)
├── genomic_index/                   (group — eager, parallel scalars;
│   │                                 same element types as §3.2, extendable)
│   ├── lengths            uint32  1-D
│   ├── positions          int64   1-D
│   ├── mapping_qualities  uint8   1-D
│   ├── flags              uint32  1-D
│   ├── chromosome_ids     uint16  1-D
│   └── chromosome_names   compound[(name, VL_STRING)]
├── signal_channels/                 (blocks' blobs back to back per channel)
│   ├── sequences/         (GROUP, whatever the codec)
│   │   └── data           uint8  1-D
│   ├── qualities          uint8  1-D
│   ├── read_names         uint8  1-D
│   ├── cigars             uint8  1-D
│   └── mate_info/         (GROUP)
│       ├── inline_v2      uint8  1-D
│       └── chrom_names    compound[(name, VL_STRING)]    (chrom_id → name)
└── provenance/                      (group — optional)
    └── steps              compound
```

Key facts for `blocks_v1`, verified against the writer:

- **Every blob channel is codec-coded.** Per-block codec ids live in
  the block index; the `@compression` attribute on a channel dataset
  carries the first block's codec and is informative only. Defaults
  when the caller sets no override: `sequences` → REF_DIFF_V2 (14)
  with a reference, SEQ_CM (19) without (RANS_ORDER1 (5) when the
  native library is not loaded); `qualities` →
  FQZCOMP_NX16_Z (12), unconditionally; `read_names` →
  NAME_TOKENIZED_V2 (15); `cigars` → RANS_ORDER0 (4); `mate_info` →
  MATE_INLINE_V2 (13).
- **A block never spans two chromosomes** — the writer flushes the
  pending block at every chromosome change, so a coordinate-sorted
  whole-genome BAM streams through REF_DIFF_V2 as long
  same-chromosome stretches. Unmapped reads inside a mapped block
  (FLAG `0x4`, CIGAR `*`, placed on the mate's contig) stay in the
  REF_DIFF_V2 blob via its UL substream; fully unmapped stretches
  form their own blocks and fall back to BASE_PACK.
- **Blocks are cut** at `block_reads` / `block_bytes` bounds (and at
  chromosome changes); block boundaries do not change bytes within a
  block, so a block's blob is byte-identical to what a whole-channel
  writer would produce for a run holding that block's reads alone.
- **Partial files are readable**: a writer that stops before `close`
  leaves the block index and datasets agreeing up to the last flushed
  block; readers use the index row count, not `@read_count`, when the
  two disagree.

### 3.2 Legacy whole-channel layout

```
/study/genomic_runs/<name>/
├── genomic_index/                       (group — eager, parallel scalars)
│   ├── lengths            uint32  1-D    (per read; offsets derived at read time)
│   ├── positions          int64   1-D    (0-based mapping position)
│   ├── mapping_qualities  uint8   1-D    (Phred MAPQ)
│   ├── flags              uint32  1-D    (SAM flags)
│   ├── chromosome_ids     uint16  1-D    (interned chrom id per read)
│   └── chromosome_names   compound[(name, VL_STRING)]  (unique names, id order)
├── signal_channels/                     (group — lazy, bulk + VL data)
│   ├── sequences/         (GROUP)        (REF_DIFF_V2 layout; v1 was a flat dataset)
│   │   └── refdiff_v2     uint8  1-D     @compression = 14
│   ├── qualities          uint8  1-D     @compression = 12 (FQZCOMP_NX16_Z)
│   ├── read_names         uint8  1-D     @compression = 15 (NAME_TOKENIZED_V2)
│   ├── cigars             compound[(value, VL_STRING)]   (or uint8 + rANS override)
│   └── mate_info/         (GROUP)
│       ├── inline_v2      uint8  1-D     @compression = 13 (MATE_INLINE_V2)
│       └── chrom_names    compound[(name, VL_STRING)]    (chrom_id → name)
└── provenance/                          (group — optional)
    └── steps              compound       (per-run provenance records)
```

Key layout facts, verified against the writer:

- **Integer per-record fields live only under `genomic_index/`** —
  `positions`, `flags`, `mapping_qualities` are *not* duplicated under
  `signal_channels/` (a v1.6 change; see `format-spec.md` §10.7). The
  writer rejects `signal_codec_overrides` for those names, and the
  README's "per-read parallel arrays (positions, flags,
  mapping_qualities)" describes the `genomic_index/` contents, not
  `signal_channels/`.
- **`sequences` is a GROUP, not a dataset**, in the REF_DIFF_V2
  layout: it contains a single `refdiff_v2` child dataset tagged
  `@compression = 14`. Readers dispatch on the link type (group → v2
  path; dataset → fallback path). When the run has no reference, the
  writer falls back to a flat `sequences` dataset encoded with
  BASE_PACK. Unmapped reads (CIGAR `*`) no longer force that
  fallback — REF_DIFF_V2 carries them in its UL substream
  (`codecs/ref_diff_v2.md` §4.4). Multi-chromosome runs raise on this
  whole-channel path; write them through the default `blocks_v1`
  layout (§3.1) instead, which cuts a block at every chromosome
  change.
- **`mate_info` is a GROUP** containing the `inline_v2` blob
  (`@compression = 13`) plus a `chrom_names` compound sidecar that maps
  chrom id → name. The sidecar covers mate-only chromosomes that no own
  read aligns to (and would therefore be absent from
  `genomic_index/chromosome_names`).
- **`read_names` is a flat uint8 dataset** carrying the
  NAME_TOKENIZED_V2 stream (`@compression = 15`). The legacy M82
  compound `{value: VL_STRING}` layout is rejected by v1.0 readers.
- **`cigars`** defaults to a compound `{value: VL_STRING}` dataset; with
  a `signal_codec_overrides["cigars"]` of `RANS_ORDER0`/`RANS_ORDER1`
  it becomes a flat uint8 dataset holding a length-prefix-concat byte
  stream (`varint(len) + ASCII bytes` per CIGAR) under the rANS codec.
- The `@compression` attribute is a single `uint8` (the M79 codec id).
  Channels with a codec carry no HDF5 filter (the codec output is
  high-entropy); uncompressed channels carry the modality's default
  ZLIB filter and no `@compression` attribute.

The byte channels are addressed by the `genomic_index` `offsets`/`lengths`
pair: read `i`'s sequence and qualities are the slices
`[offsets[i], offsets[i] + lengths[i])` of the (decoded) `sequences` and
`qualities` channels.

### Embedded references

When `embed_reference=True` and a context-aware sequences codec applies,
`_embed_references_for_runs` writes each unique reference (deduped by
`reference_uri`) once at `/study/references/<reference_uri>/`, with an
`md5` attribute and a `chromosomes/<chrom>/data` uint8 dataset per
covered chromosome. The default is `embed_reference=False` (external
reference, matching CRAM 3.1): the file records `reference_uri` and
`reference_md5` only, and the reader resolves the bytes externally
(§5). The same URI mapping to two different MD5s is a hard error.

## 4. Query APIs (region / unmapped / flag)

Selective access is served by the eagerly-loaded `GenomicIndex` without
materialising any `AlignedRead`. All three methods return plain
`list[int]` read indices:

- `GenomicIndex.indices_for_region(chromosome, start, end)` — indices on
  `chromosome` with `start <= positions[i] < end`. Disk-loaded indexes
  resolve the name to its interned uint16 id once, then compare the id
  column vectorised; in-memory indexes scan the names list. Returns `[]`
  when the chromosome name is unknown.
- `GenomicIndex.indices_for_unmapped()` — indices where SAM flag `0x4`
  (unmapped) is set: `np.where(flags & 0x4)`.
- `GenomicIndex.indices_for_flag(flag_mask)` — indices where
  `(flags & flag_mask) != 0`, for arbitrary SAM flag-bit queries.

`GenomicRun.reads_in_region(chromosome, start, end)` is the
read-materialising wrapper over `indices_for_region`. **Overlap
semantics:** filtering is by **mapping start position only**, not by the
read's end coordinate — a read whose start lies outside the window but
whose alignment extends into it is *not* returned. Full SAM-style
interval-overlap is a possible future enhancement.

## 5. Codecs per channel

Each `signal_channels/` channel has a default codec; `sequences`,
`qualities`, and `cigars` additionally accept a per-channel override.
Codecs are dispatched through the shared codec registry and are
byte-exact across the Python / Objective-C / Java SDKs.

| Channel       | Default codec / id                                            | Override surface                          | Codec doc |
|---------------|---------------------------------------------------------------|-------------------------------------------|-----------|
| `sequences`   | REF_DIFF_V2 (14) with a reference; SEQ_CM (19) without under `blocks_v1` (RANS_ORDER1 (5) without the native library), BASE_PACK (6) under whole-channel | RANS_ORDER0 (4) / RANS_ORDER1 (5) / BASE_PACK (6) / SEQ_CM (19) | [`codecs/ref_diff_v2.md`](codecs/ref_diff_v2.md), [`codecs/seq_cm.md`](codecs/seq_cm.md) |
| `qualities`   | FQZCOMP_NX16_Z (12)             | RANS_ORDER0 / RANS_ORDER1 / BASE_PACK / QUALITY_BINNED (7) / FQZCOMP_NX16_Z | [`codecs/fqzcomp_nx16_z.md`](codecs/fqzcomp_nx16_z.md) |
| `read_names`  | NAME_TOKENIZED_V2 (15)          | none (auto-only)                          | [`codecs/name_tokenizer_v2.md`](codecs/name_tokenizer_v2.md) |
| `cigars`      | RANS_ORDER0 (4) under `blocks_v1`; compound VL_STRING (uncoded) under whole-channel | RANS_ORDER0 / RANS_ORDER1                 | [`codecs/rans.md`](codecs/rans.md) |
| `mate_info`   | MATE_INLINE_V2 (13)             | none (auto-only)                          | [`codecs/mate_info_v2.md`](codecs/mate_info_v2.md) |

Notes:

- **REF_DIFF_V2** (sequences) is context-aware: encoding needs the
  per-read positions, CIGARs, and the per-chromosome reference
  sequence. One blob codes one chromosome; under `blocks_v1` the
  writer cuts a block at every chromosome change, so multi-chromosome
  runs stream through it, while the whole-channel path raises on a
  multi-chromosome run. Unmapped reads (CIGAR `*`) are carried by the
  codec itself — bases as soft clip, lengths in the per-slice UL
  substream — rather than forcing a fallback. A missing reference
  falls back to RANS_ORDER1 (`blocks_v1`) or BASE_PACK
  (whole-channel).
- **FQZCOMP_NX16_Z** (qualities) carries its sibling-channel inputs
  (read lengths, reverse-complement flags from `flags & 16`) inside
  its wire format. Under `blocks_v1` it is the unconditional
  qualities default; on the legacy whole-channel path it auto-applies
  only when the run is already a "v1.5 candidate" (i.e. the sequences
  channel is going through REF_DIFF_V2), preserving byte-parity for
  plain reference-less M82 writes. When the run carries a
  base-parallel `sequences` channel the encoder also tries the
  sequence-context V5 strategies and keeps the smallest stream;
  `opt_disable_qualities_v5` pins the V4 strategy set
  (`codecs/fqzcomp_nx16_z.md`).
- **Override validation is content-aware**: applying BASE_PACK or
  QUALITY_BINNED to `cigars` (ASCII strings, not ACGT bytes or Phred
  values) is rejected with a named error, as is QUALITY_BINNED on
  `sequences`. Per-field `mate_info_chrom/pos/tlen` overrides and the
  bare `mate_info` key are rejected — the inline codec encodes all three
  mate fields together.

The rANS, BASE_PACK, and QUALITY_BINNED codecs are documented at
[`codecs/rans.md`](codecs/rans.md),
[`codecs/base_pack.md`](codecs/base_pack.md), and
[`codecs/quality.md`](codecs/quality.md); the integer
DELTA_RANS_ORDER0 codec (used for sorted integer channels) is at
[`codecs/delta_rans.md`](codecs/delta_rans.md). The on-disk
`@compression` scheme and per-channel wire formats are specified in
`format-spec.md` §§10.5–10.11; the `blocks_v1` block layout in §10.12.

## 6. Import path — BAM / CRAM / SAM → model

The importers wrap the user-installed **`samtools`** binary as a
subprocess (no htslib is linked); SAM/BAM/CRAM parsing follows the public
SAMv1 specification.

- `ttio.importers.bam.BamReader(path)` reads SAM or BAM (samtools
  auto-detects the format from magic bytes). `to_genomic_run(name=...,
  region=..., sample_name=...)` runs `samtools view -h <path> [region]`
  and returns a `WrittenGenomicRun`. `ttio.importers.sam.SamReader` is a
  discoverable convenience alias for the same class.
- `ttio.importers.cram.CramReader(path, reference_fasta)` subclasses
  `BamReader` and reuses its SAM-text parsing; the only difference is
  that its `samtools view` invocation adds `--reference <fasta>` so the
  reference-compressed CRAM bytes are reconstituted. The reference FASTA
  is a required positional argument.
- `samtools` is a **runtime** dependency only: importing the module
  succeeds without it; `to_genomic_run` raises `SamtoolsNotFoundError`
  (with apt/brew/conda install guidance) at first use if the binary is
  missing.

`to_genomic_run` materialises the whole input in memory. For inputs
that should not be held whole, the readers also stream:
`BamReader.iter_batches(batch_reads=100_000, ...)` yields the input as
consecutive `WrittenGenomicRun` batches parsed from the `samtools
view` pipe line by line, and `BamReader.stream_source(...)` wraps that
iterator as a `GenomicStreamSource` that feeds a `GenomicStreamWriter`
(§2.5) batch by batch — `reference_fasta=` enables REF_DIFF_V2 through
a `LazyReference` without loading the FASTA whole. The FASTQ and
FASTA importers expose the same `iter_batches` / `stream_source`
surface (the FASTA batches carry the unaligned-import sentinels).

Column mapping (SAM fields 1–11 → model, then the optional tags):

| SAM field | Model destination                       |
|-----------|-----------------------------------------|
| QNAME     | `read_names[i]`                         |
| FLAG      | `flags[i]`                              |
| RNAME     | `chromosomes[i]`                        |
| POS       | `positions[i]`                          |
| MAPQ      | `mapping_qualities[i]`                  |
| CIGAR     | `cigars[i]` (kept literally, incl. `*`) |
| RNEXT     | `mate_chromosomes[i]` (`=` expands to RNAME) |
| PNEXT     | `mate_positions[i]`                     |
| TLEN      | `template_lengths[i]`                   |
| SEQ       | concatenated into `sequences`; `*` → 0 bytes |
| QUAL      | concatenated into `qualities`; `*` → `0xFF`×len(SEQ), or 0 bytes if SEQ also `*` |
| 12+       | `tags[i]`: the optional fields, tab-joined as samtools prints them (`""` when none); stored as the SAM_TAGS `tags` channel (M101, format-spec §10.13) |

Per-read `offsets`/`lengths` are accumulated from the SEQ byte length as
records stream in (a `SEQ=*` record contributes a zero-length slice — the
offsets/lengths pair carries the "absent" signal). Header handling:
`@SQ SN:` names populate the reference candidates (first `@SQ` wins for
`reference_uri`); the first `@RG SM:`/`PL:` tags seed `sample_name` and
`platform`; `@PG` lines become provenance records (timestamped from the
file mtime). All imports are tagged `AcquisitionMode.GENOMIC_WGS`.

### Reference resolution

On the read side, REF_DIFF_V2 decode resolves a chromosome's reference
bytes through `ttio.genomic.reference_resolver.ReferenceResolver`, with a
fixed lookup chain (hard error on miss — no partial decode):

1. **Embedded** — `/study/references/<uri>/chromosomes/<chrom>/data` in
   the open `.tio` file (the `md5` attribute is verified against the
   codec's expected MD5).
2. **External FASTA** — an explicit `external_reference_path`, else the
   `REF_PATH` environment variable; the extracted chromosome's MD5 is
   verified.
3. **`RefMissingError`** — when neither resolves, or on any MD5 mismatch.

The resolver navigates `/study/references` through the `StorageGroup`
protocol; the `GenomicRun` threads that group in at `open` time
(`references_group=`), so embedded-reference decode works across every
storage provider, not just raw HDF5.

## See also

- [`format-spec.md`](format-spec.md) §§10.5–10.12 — on-disk `@compression`
  scheme and per-channel wire formats.
- [`providers.md`](providers.md) — the `StorageProvider` / `StorageGroup`
  / `StorageDataset` protocol the model is stored through.
- [`codecs/ref_diff_v2.md`](codecs/ref_diff_v2.md),
  [`codecs/mate_info_v2.md`](codecs/mate_info_v2.md),
  [`codecs/fqzcomp_nx16_z.md`](codecs/fqzcomp_nx16_z.md),
  [`codecs/name_tokenizer_v2.md`](codecs/name_tokenizer_v2.md) — the
  genomic codecs.
