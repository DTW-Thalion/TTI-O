# CHANGELOG

All notable changes to the TTI-O multi-omics data standard reference
implementation.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning follows [Semantic Versioning](https://semver.org/); the
public API is stable from onward.

---

## [Unreleased]

### Added
- **M103 SEQ_CM, read bases without a reference.** Runs with no
  reference (FASTQ imports, unaligned reads) now code their bases with
  SEQ_CM (codec id 19) under `blocks_v1`, instead of rANS order-1: a
  context-mixing model over 11-, 16- and 24-base contexts with
  reverse-complement training and a binary arithmetic coder, in one
  native kernel shared by the three SDKs (Python `ttio.codecs.seq_cm`,
  Java `codecs.SeqCm`, ObjC `TTIOSeqCm`). On the M103 Phase 0 data the
  bases take 0.35-0.62 bits each against about 1.94 before. Writers
  without the native library keep RANS_ORDER1, recorded per block;
  SEQ_CM is also accepted as a `sequences` override on the
  whole-channel layout. Readers that predate it reject codec id 19 by
  name (binding decision §104; `docs/codecs/seq_cm.md`).
- **M102 imzML pixel coordinates.** Large imaging imports no longer
  fail. The Python imzML importer stored every pixel's position in one
  provenance parameter, `imzml_pixel_coordinates_csv`, which both
  provenance levels mirror into a fixed-length `@provenance_json`
  attribute. HDF5 caps an attribute at 64 KB, so an import failed with
  "object header message is too large" from about 9,000 pixels (PRIDE
  PXD001283, 34,840 pixels, could not be imported). Positions are now
  optional `spectrum_index/pixel_x|pixel_y|pixel_z` int32 columns,
  flagged `opt_pixel_coordinates` (format-spec §4b), read and written
  by all three SDKs (`SpectrumIndex.pixel_x` / `pixelX()` / `-pixelX`,
  `pixel_coordinates_at`). The imzML exporters write an imported pixel
  run back to imzML when the dataset has no MS image cube, keeping its
  mode, UUID, grid and pixel size, and older files that only carry the
  CSV still read and export. The ObjC and Java importers now accept
  processed-mode imzML, as the same `imzml_pixels` run (continuous
  mode still becomes an MSImage cube there). Files without pixel runs
  are unchanged. Not yet carried by `.tis` transport. Binding
  decisions §101–§103.

- **M101 SAM optional tags.** The BAM, SAM and CRAM importers keep the
  optional fields (SAM columns 12 and up) instead of discarding them,
  so a `.tio` made from a BAM is a lossless copy of its records. Each
  read's tags are its tab-joined tag text as `samtools view` prints
  it (`WrittenGenomicRun.tags`, `AlignedRead.tags`), stored as the
  genomic `tags` channel with the new SAM_TAGS codec (id 18), and the
  exporters write them back after column 11. SAM_TAGS is a CRAM-style
  column codec in the shared native library: a tag-line dictionary per
  blob, one column per tag key, MD:Z and NM:i recomputed from the
  reference when the block's sequences are REF_DIFF_V2, and integers
  equal to an earlier tag in the read (novoalign's UQ repeats AS)
  stored as back-references. On the GIAB HG002 2x250 chr22 slice
  (10.6 M reads) the tags take 20.5 MB, 1.93 bytes per read, against
  27.3 MB for CRAM 3.1 `small` and 776.5 MB of tag text, round-trip
  byte-exact; on the bwa-aligned NA12878 WES chr22 slice CRAM is ahead
  (1.26 MB against 3.22 MB), mostly on `XA:Z` strings
  (`docs/benchmarks/2026-10-01-m101-cram-remeasure.md`). Files without tags are unchanged: no `tags` dataset, no
  `blocks/index` columns, no `opt_sam_tags` flag. Per-AU encryption
  encrypts the tags with the sequences (an MD string lists the read's
  variants), signatures cover them, plaintext and encrypted transport
  carry them, and `bam_dump` gains a `tags` key. Format-spec §10.13,
  `docs/codecs/sam_tags.md`, binding decisions §96–§100.

- **M100 writable dataset open.** Objective-C and Java gain the
  public writable reopen Python already has
  (`SpectralDataset.open(path, writable=True)`):
  `+[TTIOSpectralDataset readFromFilePath:writable:error:]` and
  `SpectralDataset.open(pathOrUrl, writable)`. A dataset opened
  writable surfaces a read-write provider, so in-place writers
  such as `ReferenceImport.writeToDataset` work on a reopened
  container instead of failing at `H5Gcreate2` for
  `/study/references/`. Provider-URL opens route the same mode
  through the provider registry; the file must exist in either
  mode. The references cross-language writers now drive the
  writable reopen in all three languages, so the 9-cell
  writer x reader matrix exercises the new entry points
  end-to-end.

- **M99.1 `blocks_v1` restore hardening and transport carriage.**
  The stream writers persist the writer policy that shapes the coded
  blobs (`@ref_diff_slice_bytes`, `@opt_disable_qualities_v5`) and,
  for REF_DIFF_V2 runs, the `@reference_md5s` chromosome map
  carrying the reference-set digest, so a per-AU restore re-encodes
  under the original policy and rebuilds the reference from
  `/study/references/` or a `REF_PATH` FASTA through the reference
  resolver (format-spec §9.1.1). The M99 encrypt-time re-encode and
  byte-compare gate is removed in all 3 SDKs, and a run written with
  `embed_reference=False` now encrypts and restores. When a
  re-encoded blob still lands on different ranges — an older file
  without the attrs — restore rewrites the block index to the ranges
  actually written instead of refusing, keeping the file consistent
  and readable. Transport-spec v0.12 adds the GenomicRunSidecar
  (0x1C) and BlockSidecar (0x1D) packets and the required
  `transport_blocks_v1` feature token, so the encrypted stream
  carries `blocks_v1` runs and the M99 sender refusal is removed; a
  received container decrypts in place byte-identically (binding
  decisions 93–95 revised). The cross-language matrix grows a 3×3
  send × receive plane: 27/27 cells byte-identical. The Java
  transport receiver now always writes `study/ms_runs` like the
  other receivers, and the Python per-AU decrypt no longer silently
  no-ops on a container without it.

- **M99 per-AU protection over `blocks_v1`.** Per-AU encryption
  raised on any genomic run in the default `blocks_v1` layout — the
  walkers only knew the legacy v1.8 whole-channel shape, so nothing
  the current genomic write path produces could be encrypted. The
  walkers in all 3 SDKs now stream `blocks_v1` runs block by block
  (format-spec §9.1.1): per-block codec decode, one AES-GCM AU per
  read with global AU numbering, extendable segments tables appended
  per block, so peak memory follows the block policy (64 MB default)
  rather than the channel size — a run with 100 MB of decoded
  channels encrypts with no measurable RSS growth over baseline.
  Decrypt-in-place restores the original layout byte-identically
  (per-block re-encode under the writer's sticky qualities
  discipline). The 3×3 encryptor × decryptor cross-language matrix
  is byte-identical in every cell, on a cross-chromosome-mates
  fixture and an embedded-reference REF_DIFF fixture. At M99 the
  restore contract was enforced by an encrypt-time re-encode and
  byte-compare gate, an embedded-reference requirement and a
  sender-side transport refusal; M99.1 (above) replaces all three.

- **M98 AssemblyGraph.** One new container primitive: a GFA 1.x
  assembly graph at `/study/assembly_graphs/<name>/` (format-spec
  §11a), parsed into `segments`/`links`/`paths`/`extras` tables plus
  a `line_index` that replays the original line order, so re-emission
  is byte-exact against the source file — proven on hifiasm output
  and a synthetic covering every 1.x line type. `GfaReader` /
  `GfaWriter` in all 3 SDKs (Apache-2.0 import/export set); segment
  sequences store through the existing codec stack (BASE_PACK for
  ACGTN alphabets, RANS_ORDER1 otherwise, `@compression` on the
  dataset). `write_minimal` / `create` / `writeMinimalToPath:` take
  an `assembly_graphs` mapping and stamp the `opt_assembly_graph`
  feature flag only when a graph is present; datasets expose an
  `assembly_graphs` accessor on open. The `GfaDump` CLIs
  (`python -m ttio.importers.gfa_dump`, `TtioGfaDump`, Java
  `GfaDump`) emit one canonical JSON shape and carry `--write-tio` /
  `--emit-gfa` conformance modes; the 3×3 writer × emitter container
  matrix is byte-exact in every cell. Per-AU encryption covers the
  sequences channel with one AES-GCM operation per segment record
  and decrypts in place in all 3 languages; `v2:`/`v3:` signatures
  cover it the way they cover genomic-run channels
  (`sign_assembly_graph` / `signAssemblyGraph`).

- **Per-AU CLI `decrypt-in-place`.** The three PerAU conformance CLIs
  (`python -m ttio.tools.per_au_cli`, `TtioPerAU`, Java `PerAUCli`)
  gain a `decrypt-in-place <file> <key-file>` subcommand over the
  existing in-place APIs. The M98 matrix gains 3×3 encryptor ×
  decryptor cells proving a graph container encrypted by one SDK
  decrypts in place under any other and re-emits GFA byte-exactly,
  with the per-AU feature flags stripped.

- **M97 long-read profile.** Genomic runs gain the `@read_role` UTF-8
  attribute (`hifi`, `ont_ul`, `hic_r1`, `hic_r2`, `parental_maternal`,
  `parental_paternal`, `illumina_wgs`; other strings stored unchecked),
  set via `WrittenGenomicRun.read_role` / `readRole` and the stream
  writers, read back via `GenomicRun.read_role` / `getReadRole()` /
  `-readRole` (absent → `None`/`null`/`nil`); the FASTQ/FASTA import
  CLIs take `--read-role`. REF_DIFF_V2 accepts a `slice_bytes` byte
  budget (`WrittenGenomicRun.ref_diff_slice_bytes` /
  `refDiffSliceBytes`): a slice closes before the read that would push
  it past that many bases, so a 10,000-read ONT-UL run no longer
  produces a gigabyte-scale decode unit. Writer policy only — the wire
  format and decoder are unchanged, 0 keeps the fixed-count rule, and
  a cross-language fixture pins the non-uniform boundaries byte-exact
  in all 3 SDKs. A `QUALITY_BINNED` qualities override is now rejected
  when `@platform` names a HiFi/PacBio/ONT/Nanopore instrument — the
  fixed Illumina-8 bin table does not fit those quality distributions
  (`quality.binned_allowed_for_platform` and equivalents). The `.tis`
  genomic run-metadata JSON carries a `read_role` key (plain and
  encrypted transport; always present, `""` means no role, receivers
  leave the container attribute absent on `""`), so the attribute
  survives transport round-trips. HiFi- and ONT-UL-shaped fixtures
  run the 3x3 transport matrix, which asserts the attribute in every
  cell.

- **The FASTA importer streams.** The unaligned path gains the batch
  surface the FASTQ and BAM readers already have: Python
  `FastaReader.iter_batches` / `stream_source`, Java `iterBatches` /
  `stream`, Objective-C `+iterBatchesFromPath:…usingBlock:` /
  `+streamFromPath:…`. Batches cut at whichever of the read-count /
  byte limits is hit first (100 000 reads / 64 MiB by default), carry
  the SAM-unmapped sentinels of the eager read, and concatenate to
  exactly its run — including CRLF line ends, blank lines, multi-line
  records and a final line without a newline, which the new tests pin
  against the eager walk in all 3 SDKs.
  `python -m ttio.tools.fasta_import_cli unaligned` now writes
  through a `GenomicStreamWriter` and takes `--block-reads`,
  `--block-bytes` and `--legacy-whole-channel`, so a FASTA of any
  size imports with bounded memory.

- **`-[TTIOAcquisitionRun unitColumnForChannel:valueStart:]`, the unit's
  decoded values without building a spectrum.** Only a view handed to an
  `-iterBlocksFrom:to:threads:error:usingBlock:` visitor has columns;
  every other run returns nil.

  Constructing spectra is 73.7% of the ordered reader's wall time on
  `5_TAP_Pdr1_R376W_PTM` (342546 spectra, 256 MB): decoding and scanning
  the same bytes runs at 820.7 MB/s against 215.5 through
  `-iterSpectraWithBatch:`, single-threaded, minima of 5 interleaved
  passes with a contended control row. Broken down per step, the two
  `NSData` copies are 23% of construction, the two `TTIOSignalArray`
  wrappers 30%, the channel dictionary 11% and the `TTIOMassSpectrum`
  itself 36%, so two thirds of the cost is the objects rather than the
  data, and no amount of copy elimination reaches it.

  Through `iterBlocks`, reading columns instead of spectra measures
  207.4 to 1050.9 MB/s at one thread, 439.0 to 1188.3 at four and 771.4
  to 1637.8 at sixteen. One thread reading columns beats sixteen threads
  building spectra. A caller doing bulk arithmetic over a channel should
  read the column; a caller that needs per-spectrum metadata still wants
  the spectra.

- **`-[TTIOAcquisitionRun iterBlocksFrom:to:threads:error:usingBlock:]`,
  a parallel block consumer for spectral runs.** The visitor is called
  once per scheduling unit, on a thread pool, out of order; that relaxed
  ordering is the whole difference from
  `-iterSpectraWithBatch:threads:`, whose ordering guarantee is what
  costs the parallelism. Spectrum `k` of a call is
  `[view spectrumAtIndex:viewStart + k]` and its run index is
  `firstSpectrum + k`.

  Spectral blocks are cut on value boundaries while a spectrum is
  located by value offset, so block edges fall inside spectra. A
  spectrum belongs to the block holding its first value, which
  partitions the spectra exactly once; a unit's extent can then run past
  its own block's end, and a block holding no spectrum start emits
  nothing rather than an empty unit. Units are planned from
  `blocks/index` when present, from the FDZ1 block tables when not, and
  from fixed spectrum counts for a channel with no FDZ1 stream, so the
  method works on every MS run rather than only on one written since
  `blocks/index` existed.

  HDF5 is not thread-safe, so storage reads stay on the calling thread
  and only decode and spectrum construction go to the pool. Memory sizes
  the in-flight window, charged from the widest unit in the range across
  every channel, and `TTIO_MEMORY_BUDGET` binds it. Setting `*stop`
  stops scheduling; units already in flight still run.

  On `5_TAP_Pdr1_R376W_PTM` (Orbitrap Exploris, 342546 spectra, 256.0 MB
  of float64), minima of 5 interleaved passes: 213.3 MB/s for the
  ordered reader against 443.5 at 4 threads and 780.5 at 16, so 3.66x.
  Parallelism is capped by the unit count, 16 here, so 30 threads is
  within noise of 16. At one thread the new method is 5.6% slower than
  the ordered reader; the two cross at two threads. Objective-C only for
  now.

- **MS runs carry a `blocks/index`, the table genomic runs already
  had.** One compound row per FLOAT_DELTA_ZSTD block: the value range
  it covers, and for every signal channel the byte offset, length and
  codec of that block. The extent includes the block's 5-byte header,
  so the bytes it names are a self-describing block. Without it a
  consumer had to walk each channel's FDZ1 stream reading block
  headers to learn the same offsets; with it, one compound read plans
  a range read or a parallel decode, which is what the genomic side
  has been able to do since `blocks_v1`. Written by all three SDKs and
  read by `TTIOSpectralBlockIndex`.

  The channel set of a spectral run is not fixed, so a reader recovers
  it from the compound type rather than a constant: every
  `<name>_off` column names one channel. The group is written only for
  runs stored with codec 17 and only when every channel cut its blocks
  at the same value boundaries; a run that fell out of step gets no
  table rather than a wrong one. It is optional on read, so MS runs
  written before it exists open unchanged.

### Fixed
- **ObjC tools read freed memory at exit and occasionally died with
  SIGSEGV.** libobjc2 keeps its weak-reference table in a C++ static
  whose destructor is registered on the first weak store, which in a
  TTI-O tool happens inside `main()`; `exit()` therefore freed the
  table before `_dl_fini` ran gnustep-base's clean-up, and
  `+[NSUserDefaults atExit]` then released a weakly referenced object
  into the freed table. Valgrind saw it on every exit of
  `TtioTransportEncode` and `TtioDumpIdentifications`; it crashed about
  one CI run in 30. A libTTIO constructor (`Core/TTIOExitOrder.m`) now
  builds the table before `main()`, so it is torn down from libobjc's
  own finaliser after gnustep-base's, and CI runs the encoder under
  valgrind to keep it that way.

- **ObjC and Java compound readers built ASCII VL-string memory types.**
  HDF5 has no ASCII↔UTF-8 string conversion path, so reading an
  h5py-written compound (UTF-8 members) failed on a cold
  conversion-path cache; a warm cache aliases the pair and masked the
  bug in every earlier cross-language read. First cold reader was the
  per-AU graph walker: ObjC errored, Java swallowed the exception and
  encrypted zero AUs. Both readers now mirror the file member's
  character set (`TTIOCompoundIO` generic read, `Hdf5CompoundIO`
  VL-string pass); the encrypted-exchange matrix asserts one AU per
  segment record so an empty walk can no longer pass.

- **`TTIOSpectralStreamWriter` capped encode work in flight by count, not
  by memory.** `-_emitFdzBlock:` drained only until `threads` blocks
  remained in flight for a channel. An FDZ1 block is 2^20 float64 values
  and the method takes a copy of each, so at 30 threads across 2 channels
  that admitted about 480 MB of copies with no memory bound. The genomic
  writer has charged its in-flight work against
  `+[TTIOThreads resolveMemoryBudget:threads:blockBytes:]` since the
  writer-ceiling work; the spectral writer never did.

  Bytes in flight are now summed across every channel, because memory is
  global while the existing drain is per channel, and
  `-_drainFdzToBudget:` drains the fullest channel first until they fit
  the writer's half of the budget. The count window stays as the upper
  bound. A new `memoryBudgetBytes` option takes an explicit ceiling and
  0 resolves, matching the genomic writer, so `TTIO_MEMORY_BUDGET` means
  one thing across the library.

  On 20000 synthetic spectra the high-water mark falls from 18628608
  bytes with the resolved budget to 8388608, exactly one block, with the
  budget set to 1. Output is byte-identical either way. The writer also
  reports `maxInFlightBytesObserved`: identical output shows the writer
  is correct under a tight budget, not that the budget was consulted,
  and the first version of this change was a silent no-op that passed
  every other assertion because `-copyWithZone:` on the options class
  dropped the new property.

- **`base_peak_intensities` was computed from a freed buffer.**
  `-[TTIOSignalArray float64Buffer]` returns the array's own buffer
  when the precision is already float64 and allocates a conversion
  buffer for every other precision. Four callers scanned the result
  through a raw pointer that nothing kept alive, so ARC released the
  temporary and the scan read freed memory whenever the array was not
  float64. Real mzML stores intensity as float32 and took that path on
  every spectrum: on the 1.73 GB Orbitrap import the column came out
  with a maximum of 6.859e+78 against a true 3.111e+09.

  The damage was confined to `spectrum_index/base_peak_intensities`, a
  queryable column, so the base-peak filter returned wrong results;
  channel data goes through a different call site and is unaffected.
  It was intermittent, because freed memory usually still holds the
  old bytes. Two sites were in `TTIOWrittenSpectralBatch` and two in
  the eager writer in `TTIOAcquisitionRun`. Repeated imports of one
  file now produce byte-identical output, so the run to run variation
  in .tio size was this bug rather than the threaded writer.

- **Three tests scanned a buffer nothing kept alive.** An audit of the
  Objective-C tree for the shape above — a raw pointer taken from a
  method's return value and read after the statement that produced it
  — found three sites, all in the test suite. Two are in
  `TestSpectralStreamWriter`, which reads `-float64Buffer` inside the
  `iterSpectra` blocks; one is in `TestM90GenomicProtection`, which
  reads the ASCII bytes of a masked read's sequence, where
  `-dataUsingEncoding:` allocates on every call and the pointer was
  therefore always dangling. Each now holds the `NSData` in a local for
  the life of the scan. The masking check also derives its loop bound
  from the data it scans and requires the expected length, rather than
  reading a fixed eight bytes from a buffer it never measured.

  The library sources are clean. Every remaining pointer of this shape
  either reads through a property backed by an instance variable or an
  element of a collection the caller holds, or is used inside the
  expression that produced it, which ARC keeps alive until the full
  expression ends.

### Changed
- **`TTIOImzMLReader` reads the `.ibd` at the offsets the metadata
  names instead of loading the whole file.** The reader called
  `+[NSData dataWithContentsOfFile:]` on the binary and then cut each
  pixel's arrays out of it with `-subdataWithRange:`, so peak memory
  held the entire `.ibd` on top of the materialised pixels. Imaging
  `.ibd` files run to many GB. Each array is now read with `pread`
  into a buffer handed to `+[NSData dataWithBytesNoCopy:length:freeWhenDone:]`,
  which is what the Java and Python importers have always done:
  `RandomAccessFile.seek` and `fh.seek` respectively. This closes a
  gap where Objective-C was the only SDK that loaded the file whole.

  On a synthetic continuous-mode image of 20000 pixels by 25000 points,
  a 4.0 GB `.ibd`, peak resident set falls from 7869984 KB to 3967776 KB
  and wall time from 4.39 s to 2.37 s, minima of 3 runs. The remaining
  term is the array of materialised pixels, which the return type of
  `+readFromImzMLPath:ibdPath:error:` requires; the file itself no
  longer contributes. The pixel checksum is identical before and after.

  `-[NSFileHandle readDataOfLength:]` was not used for the read because
  GNUstep returns `NSMutableDataMalloc` from it, and `-copy` on mutable
  data allocates. The pixel initialiser copies what it is handed and
  continuous mode requires every pixel to alias one m/z object, so that
  route silently broke the aliasing. `+dataWithBytesNoCopy:` yields
  immutable data, on which `-copy` is identity.

- **FLOAT_DELTA_ZSTD (codec id 17) chooses the byte-plane transpose per
  block instead of always applying it.** Bit 1 of the per-block
  transform byte now means the values enter the zstd frame as plain
  little-endian uint64; the encoder compares all four combinations of
  that bit with the delta bit and keeps the smallest, as it already did
  for the delta alone. The transpose pays on intensity arrays and costs
  on m/z: on `5_TAP_Pdr1_R376W_PTM.mzML` (Orbitrap Exploris, 1.73 GB,
  342546 spectra) the signal channels go from 72509598 to 53494074
  bytes, 26.2% smaller, and the whole .tio by about 23%. Intensity
  blocks still choose the transpose; m/z blocks no longer do. The two
  extra zstd passes per block cost about 1.5% of import wall time,
  39.74 s against 40.32 s on the minima of 3 and 4 runs.

  A stream carrying the new bit is rejected by readers built before
  this change, which report the unknown transform and stop. Streams
  already written use transform 0x00 or 0x01 and decode unchanged.

- **mzML is parsed incrementally.** `TTIOMzMLReader` read the whole
  document with `+[NSData dataWithContentsOfFile:]` and parsed it with
  `NSXMLParser`, which does not set `XML_PARSE_HUGE`. A 1.73 GB
  Orbitrap Exploris run with 342546 spectra and uncompressed arrays
  failed with libxml2's "Huge input lookup" after 339968 spectra, exit
  2, 21.4 GB peak RSS, leaving a partial .tio. `TTIOXMLStreamParser`
  drives libxml2's `xmlTextReader` and sends the same
  `NSXMLParserDelegate` callbacks, so the reader's element handling is
  unchanged: that import completes at 255 MiB peak RSS in 39.8 s, the
  minimum of 3 runs, and a 0.94 GB timsTOF PASEF run takes 25.3 s at
  249 MiB. The Python and Java readers already parsed from the file.


### Added
- **M94.Z V6, a segmented qualities variant.** V4 and V5 code a block
  as one adaptive chain, so a 64 MiB block cannot be split across
  compute units, and the qualities codec is about 73% of encode CPU.
  V6 splits a block into independent segments at read boundaries, each
  with its own model and coder, so one block spreads across cores. It
  is reached only by strategy hint 8 (`HINT_V6` in Python and Java,
  `TTIOM94ZHintV6` in Objective-C); auto-tune never selects it and V5
  remains the default, so nothing already written changes. Output bytes are
  identical for any thread count, and a cross-SDK golden fixture pins
  them across Python, Java and Objective-C. V6 decodes without the
  sequences channel.

  It trades ratio for that parallelism: measured against V4 at the
  default 1 Mi segment, +3.75% on low-coverage chr22, +4.31% on 2x250
  chr22, +5.42% on HiFi and +5.89% on NovaSeq WGS. The cost is
  per-segment model warm-up and falls as segments grow.

  The default segment size is 1 Mi, not the 256 Ki first shipped in this
  branch. 256 Ki was sized against the number of chains a GPU could run
  at once, and that backend was removed; swept at the shipped context,
  1 Mi is better on ratio for all four corpus classes and on encode
  throughput, and still leaves 64 segments in a 64 MiB block. Segment
  size is written into the body header, so streams already written keep
  their own and decode unchanged. Full table and method in
  `docs/codecs/m94z_v6.md`.

- **Two environment variables that make the writer's thread split
  measurable.** `TTIO_V6_SEGMENT_THREADS` overrides
  `clamp(cores / workers, 2, 8)`, and reaches outside that clamp on
  purpose: a sweep has to be able to ask for one segment thread and for
  more than eight. `TTIO_M94Z_HINT` pins the qualities strategy before
  block 0 rather than after it, which is the only way to reach V6
  through a stream writer, since auto-tune never selects it;
  `TTIO_M94Z_EXHAUSTIVE` still wins. Both are read in all three SDKs,
  Java also taking `-Dttio.v6.segmentThreads` and `-Dttio.m94z.hint`. A
  positive integer wins and anything else falls through to the existing
  rule, matching `TTIO_READ_AHEAD_BLOCKS`. Defaults are unchanged.

### Fixed
- **Exports stream their output.** The FASTQ and FASTA exporters in all
  three SDKs built the whole output in memory before one write, so
  exporting a run cost its own size again on top of the decode (the
  suite's 20 GB HiFi decode-verify hit a 20 GB cap this way). They now
  stream records to the file (a running byte position keeps the .fai
  offsets); bytes are unchanged. The threaded sequential readers also
  held one decoded block per worker; the decode-ahead window is now at
  most 4 blocks, which a serial consumer cannot outrun anyway.

### Changed
- **The default thread count is now `cpu_count - 2`, was `cpu_count - 8`.**
  `TTIO_THREADS` unset or 0 resolves to this in all three SDKs. Measured
  on a 32-thread machine encoding a corpus with more blocks than cores,
  throughput kept climbing to roughly one writer per core: at a fixed
  two segment threads, 1462 MB/s at 32 writers against 1128 at the 24
  the old default chose. The floor of
  1 stays, so a one or two core machine still resolves to a single
  worker rather than none. Callers that set `TTIO_THREADS` or pass
  `threads=` are unaffected.

  The memory budget derives from the thread count when it is not given
  explicitly (`threads * block_bytes * 16`, capped at half of physical
  memory), so a machine that was hitting the thread-derived figure
  rather than the physical cap now allows a proportionally larger
  budget.

- **Qualities strategy is decided once per run.** The M94.Z auto-tune
  encoded every block three times (V4, V5-S5, V5-S6; smallest wins).
  The streaming writers now tune block 0 only, read the winner from
  the encoded stream and pin it for the rest of the run; later blocks
  encode once. The pin is taken from block 0 by index (a block that
  arrives at the encoder before the verdict waits for it), so repeated
  runs and all three SDKs stay byte-identical for identical input.
   (Java also )
  restores the every-block tune. Wire format and decoders unchanged;
  on a corpus whose winner drifts mid-run the pinned stream is still
  valid, merely up to about a percent larger on the affected blocks.
  New kernel entry points: strategy hint 7 (V4 with its internal
  preset selection) and  (the wire
  sniffer), exposed in all three SDKs.
- **In-flight block estimate halved.** The writers budgeted raw*4
  bytes per in-flight block; the marginal resident cost measured on
  64 MiB HiFi blocks is near raw*1.5, so the old estimate halved the
  encode concurrency the memory budget admits. Now raw*2. The 50 GB
  HiFi import: 5m 40s -> 5m 11s (153 MB/s) at a 20.1 GB peak.
- **The Python import producer parallelises.** The third phase of
  the parallel producer: shard ranges for plain FASTQ, sliced records
  for gzip, numpy-vectorised boundary scanning and slice parsing
  (clean four-line fast path; a tolerant serial fallback keeps stray
  blank lines and zero-length reads correct), batches re-chunked by
  the serial cut rule so the parallel producer is batch-for-batch
  identical to the serial reader. The 50 GB HiFi corpus at 24
  threads: 11m 19s (74 MB/s) with 32 MiB batches, 64 MiB
  blocks and an 8 GB budget; the container is identical to the ObjC
  import of the same file in every genomic dataset and attribute.
  The throughput needs the `_rans` Cython accelerator (wheels build
  it; an editable install builds it with `cythonize -i`) — without it
  the pure-Python order-1 rANS holds the GIL and the pool cannot
  help. Producer read chunks clamp to the byte budget so small
  machines are not swamped structurally.
- **The Java import producer parallelises.** The same architecture as
  the ObjC phase on the SDK's executors: shard ranges for plain files,
  sliced records for gzip, ordered future/queue assembly, byte-sized
  batches and the shared byte budget. The 50 GB HiFi corpus at 24
  threads: 8m 13s (97 MB/s) with 32 MiB batches, 64 MiB blocks, an 8 GB
  budget and a 14 GiB heap, within the phase's 1.5x-of-ObjC bar; codec
  bytes are identical to the ObjC import of the same file (the only
  container difference is the pre-existing string-attribute type
  class). On the JVM the codec's JNI workspace and input copies live
  off heap, so size `TTIO_MEMORY_BUDGET` against the heap rather than
  trusting the physical-memory default. BAM keeps its serial htsjdk
  reader: its decode dominates and is not sliceable without
  index-chunking, unlike the ObjC samtools pipe.
- **The ObjC import producer parallelises.** FASTQ imports split into
  per-thread shard ranges on record boundaries (plain files) or slice
  records off a single scan (gzip, pipes), parse on the shared pool,
  and reach the writer through an ordered assembler, so the output is
  byte for byte the serial import's; SAM field parsing moves to the
  pool the same way. Batches cut by bytes (`batchBytes`, 64 MiB
  default) and the whole pipeline stalls against one byte budget
  (`TTIO_MEMORY_BUDGET`, default sixteen blockBytes per thread clamped
  to half the physical memory), replacing the block-count window. The
  50 GB HiFi corpus at 24 threads: 13m 12s before this work, 5m 40s
  after (140 MB/s, 17.7 GB peak) with 64 MiB blocks and a 20 GB
  budget; the serial path is unchanged and 64 MiB blocks compress the
  corpus marginally better than the 256 MiB default.
- **Block-parallel encode and decode.** `GenomicStreamWriter`, `SpectralStreamWriter`
  and the sequential readers (`iter_reads`, the exporters, `iter_spectra`,
  `channel_range`) in Python, Java and ObjC run codec work for several blocks at
  once: writers do the ordered work (chromosome ids, reference md5, sequence
  number) on the caller's thread, encode blocks on a pool of `TTIO_THREADS`
  workers (unset or 0 = cores minus 8, 1 = the serial path) with at most
  `threads + 1` blocks in flight, and append them in order from the caller's
  thread, so the file is byte for byte the serial writer's; readers decode the
  next blocks ahead (at most 4 in flight). `threads=` on the writers and readers and
  `--threads` on the CLIs override the environment. While a pool exists the
  FQZCOMP auto-tune runs its candidates in sequence
  (`ttio_m94z_set_autotune_threads`). Measured on the chr22
  slices (two to three blocks, so the window barely fills): low-coverage encode
  30.7 s to 25.8 s and export 40.2 s to 30.2 s at 24 threads; the WES slice,
  dominated by one block, stays put (encode 24.1 s serial, 27.7 s threaded;
  export 13.5 s both). Outputs are byte-identical (cmp) in every pairing.
- **FQZCOMP auto-tune runs its three candidate encodes concurrently.** The V5
  auto-tune encodes every block three times (V4, then the S5 and S6
  sequence-context models) and keeps the smallest; a profile of the NA12878 chr22
  encode put 58 percent of the wall time there. The three encodes now run on
  three threads with the sequential tie order kept, so the stream is byte for
  byte the one the sequential loop chose. WES chr22 11-column encode 35.5 s to
  25.1 s, peak RSS +430 MB. `TTIO_M94Z_SEQUENTIAL=1` restores the one-thread loop.
- **The reference-set md5 of an external FASTA is cached beside it.** Writing or
  reading a run against an external reference computed the digest of the whole
  FASTA on every run (hs37d5: about 7 s, 15 percent of the chr22 encode and of
  the decode). `LazyReference` in Python, Java and ObjC keeps it in
  `<fasta>.ttio-md5` as `<hex> <size> <mtime_s>`, reuses it while size and mtime
  match, and skips the write when the location is read-only (the digest is still
  cached in the process); the writers take it from the LazyReference.

### Fixed
- **REF_DIFF_V2 carries unmapped reads.** A read whose CIGAR is `*` (an unmapped
  read placed on its mate's contig, FLAG 0x4) used to switch the whole channel,
  and under `blocks_v1` the whole block, to BASE_PACK. On paired-end data every
  block has such reads: the NA12878 WES chr22 slice has 839 in 992,974 and its
  single block coded 24 MB of sequences as BASE_PACK, so TTI-O (64.1 MB) came
  out next to BAM (66.4 MB) where CRAM 3.1 small is 32.0 MB. The shared C
  kernel now stores an unmapped read's bases as soft clip and its length in a
  new per-slice UL substream (the sub-header word that was reserved and 0; a
  slice without unmapped reads still writes 0, so every earlier blob decodes
  unchanged), and the Python, Java and ObjC writers no longer fall back on such
  runs. docs/codecs/ref_diff_v2.md sections 4.3 and 4.4.
- **External reference resolution (`REF_PATH`) validates a multi-chromosome
  FASTA.** Writers record the reference-set md5 (every chromosome, alphabetic
  order, case preserved; format-spec 10.10), but the external-FASTA branch of
  the resolvers compared it with the md5 of one upper-cased chromosome, so a run
  written against any FASTA with more than one contig failed to decode from
  `REF_PATH` (`MD5 mismatch for external reference`) unless the reference was
  embedded. Python `ReferenceResolver`, Java `codecs.ReferenceResolver` and
  ObjC `TTIOReferenceResolver` now check the single-chromosome digests first
  (raw and upper-cased) and then the reference-set digest, read the chromosome
  through the `.fai`-indexed `LazyReference` (one index and one whole-FASTA
  digest per path per process instead of a full-file scan per block), and all
  return the chromosome upper-cased. `LazyReference` gains `set_md5()` /
  `setMd5()` / `-setMD5`; the Python one builds the `.fai` itself when
  samtools is absent.

### Changed
- **Objective-C: streaming import and export; genomic runs written as `blocks_v1`.**
  The ObjC reference implementation reads the `blocks_v1` layout (format-spec
  10.12) and writes it by default: `TTIOGenomicRun` opens both layouts and
  decodes one block at a time (`iterReadsFrom:to:usingBlock:`, `layout`,
  `blockCount`, `chromosomeNames`, `close`), `+writeMinimalToPath:…genomicRuns:`
  and the storage-protocol writer route every `TTIOWrittenGenomicRun` through
  the new `TTIOGenomicStreamWriter` unless `optLegacyWholeChannel` is set. MS
  runs get `TTIOSpectralStreamWriter` and a lazy `TTIOAcquisitionRun`
  (`channelRange:start:count:error:`, `iterSpectraWithBatch:usingBlock:`; a
  codec-17 channel decodes only the FDZ1 blocks a range covers, so opening a
  dataset no longer decodes every channel). Importers stream: `TTIOBamReader
  iterBatchesWithRegion:…`/`streamWithName:…` (SAM/CRAM inherit),
  `TTIOFastqReader iterBatchesFromPath:…`/`streamFromPath:…`, `TTIOMzMLReader
  streamFromPath:…` (parser thread and bounded queue); `TTIOImportedDataset`
  carries `genomicStreams`/`spectralStreams` and `TTIOImporterRegistry
  encodeFormat:` routes BAM/SAM/CRAM and mzML through them (`TtioEncode --extra
  block_reads=`, `block_bytes=`, `legacy_whole_channel=1`, `reference=<fasta>`,
  `embed_reference=1`, `batch_reads=`, `batch_spectra=`); `TTIOBamWriter
  writeReadSideRun:…`, `TTIOFastqWriter writeReadSideRun:…` and
  `TTIOMzMLWriter writeDataset:toPath:…` write read by read / spectrum by
  spectrum. Storage providers gain extendable datasets
  (`createDatasetNamed:…extendable:`, `appendData:`, `writeSlice:atOffset:`)
  and `TTIOCompoundFieldKindUInt64`; `TTIOSignatureManager` signs
  `blocks/index` and the datasets inside a channel group; `TTIOLazyReference`
  (indexed FASTA, `.fai` written in-process when absent) supplies REF_DIFF_V2
  references chromosome by chromosome. `TTIOSignalArray float64Buffer`
  converts float32 arrays, so an mzML with 32-bit binary arrays now writes the
  values it holds. `TtioWriteGenomicFixture --blocks` writes a `blocks_v1`
  file; an ObjC test decodes the Python-written golden fixture and a Python
  test reads an ObjC-written `blocks_v1` file.
- **Java: streaming import and export; genomic runs written as `blocks_v1`.**
  The Java SDK reads the `blocks_v1` layout (format-spec 10.12) and writes it by
  default: `GenomicRun` opens both layouts and decodes one block at a time
  (`iterReads`, `layout()`, `blockCount()`, `chromosomeNames()`),
  `SpectralDataset.create` and every genomic write go through the new
  `genomics.GenomicStreamWriter` unless `WrittenGenomicRun.optLegacyWholeChannel`
  is set. MS runs get `SpectralStreamWriter` and a lazy `AcquisitionRun`
  (`channelRange`, `iterSpectra`; `channels()` keeps its meaning but decodes on
  first call, so `SpectralDataset.open` no longer decodes every channel).
  Importers stream: `BamReader.iterBatches`/`stream`, `FastqReader.iterBatches`/
  `stream`, `MzMLReader.stream`; `ImportedDataset` carries `genomicStreams` /
  `spectralStreams` and `ImporterRegistry.encode` routes BAM/SAM/CRAM and mzML
  through them (`--extra block_reads=`, `block_bytes=`, `legacy_whole_channel=1`,
  `reference=<fasta>`, `embed_reference=1`); `BamWriter`, `FastqWriter` and
  `MzMLWriter` write read by read / spectrum by spectrum. Storage providers gain
  extendable datasets (`createDataset(..., extendable)`, `append`, `writeSlice`)
  and `CompoundField.Kind.UINT64`; `SignatureManager` signs `blocks/index` and
  the datasets inside a channel group; the transport writer sends multi-block
  runs per AU. `LazyReference` (htsjdk indexed FASTA) supplies REF_DIFF_V2
  references chromosome by chromosome. tio-browser streams FASTQ imports. A Java
  test decodes the Python-written golden fixture and a Python test reads a
  Java-written `blocks_v1` file; the cross-language fixtures stay on the
  whole-channel layout until ObjC reads blocks.
- **Streaming import and export; genomic runs stored as `blocks_v1`.** Every
  Python importer (BAM/SAM/CRAM, FASTQ, mzML, Thermo RAW, Waters .raw, Bruker .d)
  and exporter (SAM/BAM, FASTQ, mzML) now streams: a run of any size is written
  and read with bounded memory (one block or batch plus, for REF_DIFF_V2, the
  reference chromosomes touched). Genomic runs are written as independently coded
  blocks (format-spec 10.12: `@layout="blocks_v1"`, a `blocks/index` compound with
  per-block byte ranges and codec ids, per-channel extendable datasets, blocks
  cut at chromosome boundaries so a coordinate-sorted whole-genome BAM streams
  through REF_DIFF_V2); every blob channel is codec-coded (cigars RANS_ORDER0,
  qualities FQZCOMP_NX16_Z, reference-less sequences RANS_ORDER1). New API:
  `ttio.genomic.GenomicStreamWriter`, `ttio.spectral_stream_writer.
  SpectralStreamWriter`, `GenomicRun.iter_reads`, `AcquisitionRun.iter_spectra` /
  `channel_range`, `BamReader.iter_batches`, `FastqReader.iter_batches`,
  `MzMLStream`, `LazyReference`; storage providers gain extendable datasets
  (`create_dataset(..., extendable=True)`, `append`, `write_slice`) and the
  `UINT64` compound field kind. `ttio encode` gains `--block-reads`,
  `--block-bytes`, `--legacy-whole-channel`, `--reference`, `--embed-reference`.
  Measured on the NA12878 chr22 low-coverage BAM: 73.5 MB as blocks_v1 vs
  113.4 MB whole-channel, peak RSS 1.8 GB vs 3.3 GB. **Reader compatibility:**
  genomic files written with default settings are not readable by releases up to
  and including v1.8.0; `opt_legacy_whole_channel=True` (CLI
  `--legacy-whole-channel`) writes the v1.8 layout. Java and ObjC read support
  for `blocks_v1` follows in the next two PRs; until then the cross-language
  fixtures use the whole-channel layout. Per-AU and region encryption of genomic
  channels operate on the whole-channel layout only.
- **FLOAT_DELTA_ZSTD (compression id 17) on the transport wire, and the default
  wire codec.** Spectral AU channels can carry one FDZ1 stream per channel, and it
  is now what `use_compression=True` / `setUseCompression(true)` / `useCompression =
  YES` emit when no codec is named; `zstd` (id 16) and `zlib` (id 1) stay
  selectable through the same selector, and all three readers decode id 17 (a
  payload whose value count differs from the channel header is rejected). Measured
  on 4,341 real MS2 AUs (PXD000001, mean 175 points): −6.0% bytes vs zstd-3 and
  −14.4% vs zlib-6. Per-run zstd dictionaries were measured in the same pass and
  not adopted: −3.5% on raw float64, and the Java reader's pure-JVM zstd
  (aircompressor) does not support dictionaries. **Reader compatibility:**
  compressed streams written with default settings are not readable by releases up
  to and including v1.8.0; name `zlib` explicitly until a deployment's readers are
  current. The three transport-encode CLIs gain `--compress
  float_delta_zstd|zstd|zlib`.

## [1.8.0] - 2026-08-16

The compression release. Spectral float64 channels gain byte-shuffle and the
FLOAT_DELTA_ZSTD codec (id 17, the MS default), the qualities channel gains the V5
sequence-context flavor, embedded references pack to 2 bits per base, and the AU
transport wire gains zstd (id 16) — measured together on the audit corpora these
take the spectral side -53% lossless and chr22 qualities -22.7%. Also: a
cross-language `dek_wrapped` storage-corruption fix, a CI-hang fix, PyPI packaging
(workflow present, not yet published), the R7/R8 coverage work, and the conformance
matrix actually running in CI. **Reader compatibility:** `dek_wrapped` changes on
disk for Java and ObjC (both now write the spec layout Python always wrote; every
reader still accepts the old one), and **files written with default settings —
codec-17 MS channels and V5 qualities where it wins — are not readable by releases
up to and including v1.7.1**; `opt_disable_float_delta` and
`opt_disable_qualities_v5` restore the old layouts per run.

### Fixed
- **Cross-language `dek_wrapped` storage corruption.** `/protection/key_info/dek_wrapped`
  was written as spec uint8 exact-length by Python but as int32 zero-padded to 4 bytes by
  Java and ObjC. Java threw `ClassCastException` on a Python-written file, ObjC truncated
  71- and 1639-byte blobs to a hard-coded 60 (silent decryption failure), and Python
  value-cast the int32 files. Java and ObjC now write spec uint8; all three readers
  precision-dispatch, so pre-fix files still read. A new `conformance/key_rotation/` suite
  covers the AES-GCM 3×3 language matrix and the PQC ML-KEM Python↔ObjC pair. (#269)
- **`TransportClient` WebSocket thread leak.** `fetchPackets` closed the client only on the
  success path, so a timed-out download leaked non-daemon Java-WebSocket threads and wedged
  the forked surefire JVM until GitHub killed the tio-browser job at its 6h limit. The close
  now runs in a `finally`, `TransportServer` selector threads are daemons, and tio-browser's
  surefire runs headless JavaFX. Leaked non-daemon threads 2 → 0. (#270)

### Changed
- **MS float64 channels default to FLOAT_DELTA_ZSTD (codec id 17).** Phase 2 of the
  float-delta spec: leaving `signal_compression` at its default on a `TTIOMassSpectrum`
  run now writes each channel as a flat uint8 FDZ1 stream instead of the chunked
  shuffle+zlib layout (−33% on the PXD000001 MS channels). Opt out per run with
  `WrittenRun.opt_disable_float_delta` (Python) / `setOptDisableFloatDelta` (Java) /
  `optDisableFloatDelta` (ObjC); NMR, vibrational and genomic channels are unchanged.
  Consequences of the whole-stream layout: readers older than this release cannot open
  default-written MS files; remote range-read consumers download the full channel at
  open, so writers targeting that use the opt-out; channel and per-AU encryption decode
  to float64 before encrypting, so an encrypt/decrypt cycle rewrites the channel in the
  plain float64 layout, as it already did for numpress. The Java PQC signer gained the
  uint8 branch its canonical reader was missing, and the ObjC HDF5 provider's dataset
  adapter now routes attribute reads to the real implementation instead of a stub that
  reported every dataset attribute as absent.

### Added
- **Qualities V5: sequence-motif context for FQZCOMP_NX16_Z (codec id 12).** The
  encoder auto-tunes across the V4 presets plus 2 sequence-context strategies
  (S5 q6/p7/s5, S6 q8/p4/s6, 18 context bits) and keeps the smallest stream by
  exact size, so files where sequence context loses stay byte-identical V4 and
  readable by existing releases. Measured on the bake-off corpora: chr22 NA12878
  -22.7% B/q (0.3600 to 0.2782), PacBio HiFi -4.4%; WES capture and HG002 2x250
  keep V4. A version-5 stream decodes against the run's decoded sequences
  channel (readers order sequences before qualities; a V5 stream without
  sequences fails with a distinct error, and readers older than this release
  reject the version byte cleanly). The writer tries sequence strategies only
  for runs with a base-parallel sequences channel, skips them below 1 MiB of
  qualities, and `opt_disable_qualities_v5` / `optDisableQualitiesV5` removes
  them per run. One native implementation serves all three languages; a shared
  golden fixture pins the decode side and a cross-language file-level edge
  covers Python-written V5 files opened by Java and ObjC. Bake-off:
  `docs/benchmarks/2026-08-16-qualities-v5-bakeoff.md`.
- **FLOAT_DELTA_ZSTD (codec id 17): lossless float64 channel codec.** Per block:
  none/delta on the uint64 bit view (picked by exact size comparison), byte-plane
  transpose, one zstd frame. Opt-in via `signal_compression="float_delta_zstd"`
  (Python) / `setSignalCompression` (Java) / `signalCompression` (ObjC); the channel
  dataset becomes a flat uint8 FDZ1 stream with `@compression = 17`, decoded once at
  open. Measured on PXD000001: the four MS channels drop 193.9 → 130.8 MB at the
  level-9 default vs the shipping shuffle+gzip pipeline, bit-exact, ~100 MB/s encode.
  Decoders in all three languages; a shared golden fixture pins the decode side
  (encoders may differ byte-wise per the spec's Option B). Shipped opt-in; the
  default flipped for MS runs in the same release (see Changed).
- **HDF5 byte-shuffle ahead of the channel compressor.** All three writers now set the
  core HDF5 shuffle filter before deflate/LZ4 on chunked numeric datasets with multi-byte
  elements (signal channels, index arrays, image cubes, 2-D NMR matrices). Measured on
  the PXD000001 Orbitrap corpus: profile m/z 174.0 → 137.1 MB (−21%), all four signal
  channels together −16%, at lower encode cost. The filter is self-describing, so every existing
  reader — including pre-change releases and plain h5py/h5dump — reads the new files
  unchanged, and old files stay valid. (#280)
- **Packed embedded references.** Reference chromosomes under
  `/study/references/…/chromosomes/…/` now store a 2-bit ACGT body plus a run mask
  (`data_packed`) instead of raw bytes when packing wins — chr22 drops 9.71 → 8.24 MB
  (−15%) with encode 13× faster. Sequences that would not benefit (soft-masked,
  IUPAC-dense) keep the raw `data` layout, chosen deterministically from content in all
  three languages, and readers fall back to `data` transparently, so old files read
  unchanged. Pre-change readers fail with a missing-`data` error on packed chromosomes
  instead of misreading bytes. Cross-language byte-exact via a shared golden stream. (#281)
- **zstd on the transport wire (compression id 16).** Spectral AU channels can now be
  zstd-compressed (opt-in per writer; zlib stays the default) and all three readers
  inflate id 16; the ImagePixel packet's long-specced `1=zstd / 2=zlib` ids are now
  actually decoded instead of rejected. Measured on PXD000001 MS2 payloads: −9% bytes
  at 3.2× the encode speed of zlib-6 (zstd-3). Deps: `zstandard` (Python),
  `io.airlift:aircompressor` (Java, pure JVM), system `libzstd` (ObjC).
- **`ttio` builds as a PyPI package.** A scikit-build-core backend, an sdist that vendors
  the `native/` sources, and cibuildwheel manylinux/macOS/Windows wheels that bundle
  `libttio_rans` into `ttio/.libs` so the native codecs resolve with no environment
  variable. Publishes to TestPyPI on a `ttio-v*` tag. (#271)

### Internal
- **NUMPRESS_DELTA demoted from the size-saving tier in the format spec.** The 2026-08
  compression audit measured it at 119.4 MB (lossy) on PXD000001 MS1 profile m/z where
  a modern lossless numerical codec stores the identical channel in 61.9 MB. The spec
  now states its role plainly: mzML / msNumpress interchange parity only. Docs change;
  no code or format behaviour touched.
- **Documentation ground-truth audit.** All 219 tracked `.md` files checked against the
  code: `@mpgo_*`→`@ttio_*` rebrand leftovers, the transport header magic, the version/JDK/
  HDF5 numbers and the mate-info MF taxonomy corrected, and the missing `ref_diff_v2`,
  `mate_info_v2` and genomic-runs docs written. (#268)
- **Coverage R7 and R8**: a live-daemon workbench-client coverage gate, plus native C and
  Cython codec coverage visibility. (#272, #273)
- **Dead `_fqzcomp_nx16_z` Cython extension removed** — the wrapper is V4-native-only and
  had zero call sites. (#274)
- **Perf baseline re-captured** at 2026-06-09 numbers; every metric was inside the ±15%
  gate before the refresh. (#267)
- **tio-browser CI hang diagnostics**: a job timeout and a thread dump on hang. (#275)

PRs: #267–#287.

## [1.7.1] - 2026-06-08

Performance campaign (post-1.7.0): repaired the cross-SDK performance suite and landed
behavior-identical perf wins across all three SDKs. **No `.tio` on-disk, transport wire, or
breaking public-API change** — every change verified byte-identical with cross-language
conformance preserved.

### Performance
- **Objective-C transport & per-AU encryption.** `.mots` transport encode ~976→629ms (a
  whole-channel read cache in `spectrumAtIndex:` eliminates ~200k per-AU HDF5 hyperslab reads);
  per-AU `encryption.encrypt` ~825→576ms via zero-copy variable-length writes in the shared
  compound-write path (`writeGeneric` — benefits **all** compound writers), `EVP_CIPHER_CTX`
  reuse, and per-AU autorelease pools. (#256, #257, #258, #260)
- **Python `delta_rans` codec.** decode ~12× faster (271→22ms, now at Java/ObjC parity), encode
  ~2× — a new optional Cython varint/int64 accelerator with a pure-Python fallback;
  byte-identical `DRA0` output. (#262)
- **Python bulk reads.** `streaming.read` ~2.2×, `genomic.read` ~3.7×, genomic random-access p99
  ~100× faster — channel columns are now read once and sliced in memory instead of per-spectrum
  HDF5 hyperslab reads. Remote (fsspec/HTTP) datasets retain lazy range requests. (#263)

### Fixed
- **Java `.tio` storage bloat.** Small spectral `.tio` files no longer carry ~8MB of unused HDF5
  metadata-block space (~8MB → ~36KB on disk). The large meta-block is now used only for
  genome-reference writes that need it. FAPL-only change — files stay valid HDF5 with
  byte-identical content. (#251, #261)

### Internal
- **Cross-SDK performance suite repaired + hardened** (`tools/perf/`). It had been silently
  green while executing nothing; it is now a trustworthy **manual** regression gate — min-of-N
  timing, an absolute floor + two-tier per-metric thresholds, real-format import / PQC /
  `mate_info_v2` benches, honest cross-SDK size metrics, and a cross-SDK parity checker. (#248,
  #249, #250, #252, #253, #254, #255, #259)
- **Cross-SDK coverage campaign R1–R6.** In-process CLI coverage (Java/Python) + happy-path
  ObjC NSTask runs; dead fqzcomp V1/V2/V3 removed across all 3 SDKs (live V4 path tested);
  coverage gates ratcheted (Java line + Python total to 0.86, ObjC gate enforced in CI,
  per-unit floors). Docs refreshed. (#239, #240–#247, #264, #265)

PRs: #239–#265.

## [1.7.0] - 2026-06-05

Object-oriented design sweep: completes the entire 2026-06-02 OO design-assessment
backlog across the Python, Java, and Objective-C SDKs, plus the associated
performance work. Pure-refactor / additive release — **no `.tio` on-disk, transport
wire, or breaking public-API change**; cross-language conformance preserved throughout.

Highlights:
- **P2.6** importer/exporter `Reader`/`Writer` registry parity across all 3 SDKs +
  tio-browser delegation.
- **P1.2–P1.4** perf: Java/ObjC GenomicIndex region/flag vectorization; Python
  DELTA_RANS numpy vectorization; GenomicRun signal-channel handle caching.
- **P2.5** shared `Image` base + `ImageKind`/generic spectral axis (MS/Raman/IR).
- **P2.7** Java `SqliteProvider` compound-JSON reader → Jackson.
- **P3.8** `SpectrumKind` enum + factory dispatch (replaces stringly-typed
  `spectrum_class` dispatch; the persisted string stays the source of truth).
- **P3.9** retire the raw-`h5py` leak (Python): `SpectralDataset.file` removed,
  signatures + provenance + genomic reference resolution routed through the
  StorageProvider protocol, `native_handle()` **deprecated**.
- **P3.10** split the `SpectralDataset` / `transport/codec.py` god-files into focused
  units (public APIs byte-identical).
- **P3.11** encapsulation parity: read-only zero-copy `SignalArray.data` view +
  `Run` protocol promoted to Stable (Python); `SignalArray.asX()` defensive copies
  (Java).

PRs: #213–#237.

Java `AcquisitionRun`'s stringly-typed `spectrum_class` dispatch (the
materialize switch and the read/write `equals` chains) now routes through a new
`Enums.SpectrumKind` enum (`fromPersisted(String)` / `persisted()`); the
persisted `@spectrum_class` attribute string remains the source of truth and is
written verbatim. `RunSelection`'s NMR discriminant was converted to the same
enum. No `.tio` wire, transport, or API-shape change. (OO-assessment P3.8.)

### Changed — `AcquisitionRun` spectrum dispatch goes through a `SpectrumKind` enum (Python)

The stringly-typed `spectrum_class` dispatch in `AcquisitionRun._materialize_spectrum`
now routes through a new `SpectrumKind` enum (`ttio.enums.SpectrumKind`) via a derived
`AcquisitionRun.kind` property, instead of comparing the raw `@spectrum_class` string in
each branch. The persisted `@spectrum_class` string is unchanged — it stays the on-disk
source of truth and is written verbatim; the enum is an in-code dispatch key only.
Unrecognized classes map to `SpectrumKind.UNKNOWN` and dispatch to the `MassSpectrum`
default, preserving the v0.1 fallback. No `.tio` wire, transport, or API-shape change.
(OO-assessment P3.8.)

### Changed — read-only `SignalArray.data` view + `Run` protocol promoted to Stable (Python)

`SignalArray.data` is now stored as a zero-copy, read-only numpy view
(`flags.writeable=False`) so callers can no longer mutate the value object in
place; an in-place write such as `sa.data[i] = x` now raises `ValueError`. The
freeze is applied to an internal view, never to the caller's source array. The
`Run` protocol is promoted from Provisional to Stable. No `.tio` wire or
API-shape change — `data` stays an `np.ndarray` field. (OO-assessment P3.11.)

### Changed — `SignalArray` typed accessors return defensive copies (Java)

Java `SignalArray.asDoubles()`/`asFloats()`/`asInts()`/`asLongs()` now return a
defensive copy of the backing array (they no longer leak the internal array by
reference, so a caller can no longer corrupt the `SignalArray`'s state).
`buffer()` remains the documented raw accessor for callers needing zero-copy
access. No public API shape or `.tio` wire change. (OO-assessment P3.11.)

### Changed — `TTIOSpectralDataset.m` god-file split: genomic-write category (ObjC)

The genomic-modality write path is extracted out of the 4388-LOC
`TTIOSpectralDataset.m` into a new Objective-C category file,
`TTIOSpectralDataset+GenomicWrite.m`, with a new internal header
`TTIOSpectralDataset+Internal.h` sharing the handful of file-static helpers
(`kTTIOFormatVersion`, `_TTIO_MakeHDF5GroupAdapter`, `isNonHdf5ProviderURL`,
and the little-endian serialisation macros) used by both translation units.
The core `.m` (init / read / encrypt / decrypt / `writeToFilePath:` instance
methods) stays in place. Pure internal restructure — no public API, `.tio`
wire, or behaviour change. (OO-assessment P3.10.)

### Changed — `SpectralDataset` god-file split: genomic-write + metadata-IO helpers (Java)

The Java `SpectralDataset` static write/IO machinery is extracted into two new
package-private helper classes, `SpectralDatasetGenomicWriter` (genomic-run
write statics) and `SpectralDatasetMetadataIO` (metadata + subjects/samples
statics). The public `SpectralDataset` class (instance accessors, public
`open`/`create`/`createWithImages`, encrypt/decrypt/close, and create
orchestration) stays in place and calls into the helpers. Pure internal
restructure — no public API, `.tio` wire, or behaviour change. (OO-assessment
P3.10.)

### Changed — `spectral_dataset.py` god-file split: genomic-write + metadata-IO submodules (Python)

The genomic-run write helpers and the metadata + subjects/samples IO helpers
are extracted out of `spectral_dataset.py` into two new private submodules
(`_dataset_write_genomic.py`, `_dataset_write_metadata.py`). The
`SpectralDataset` class and its write orchestration stay in place and call
into the submodules. Pure internal restructure — no public API, `.tio` wire,
or behaviour change. (OO-assessment P3.10.)

### Changed — `transport/codec.py` god-file split: writer / reader / common modules (Python)

`transport/codec.py` is split into three new private modules: `_writer.py`
(`TransportWriter` + writer-only helpers), `_reader.py` (`TransportReader` +
reader/ingest helpers), and `_common.py` (the wire-mapping constants and the
shared genomic/codec helpers used by both). `codec.py` is now a thin
re-export facade that keeps every historical
`from ttio.transport.codec import …` path working. Pure internal
restructure — no public API, `.tio` wire, or behaviour change.
(OO-assessment P3.10.)

### Deprecated — `StorageProvider.native_handle()` (Python)

`StorageProvider.native_handle()` — the raw-backend escape hatch — is now
deprecated and emits a `DeprecationWarning`. Reach storage through
`root_group()` / the `StorageGroup` protocol instead. The method and its return
values are unchanged for now; this mirrors the Java SDK's
`@Deprecated(forRemoval=true)` to keep the three SDKs in parity, with hard
removal deferred to a future coordinated major. All mainline TTI-O code already
routes through the protocol, so the warning never fires on our own code.
Completes OO-assessment P3.9.

### Changed — genomic `ReferenceResolver` resolves embedded references via the StorageGroup protocol (Python)

The genomic REF_DIFF `ReferenceResolver` now navigates the
`/study/references` subtree through the `StorageGroup` protocol
(`has_child`/`open_group`/`open_dataset` + `StorageDataset.read`) instead of a
raw `h5py.File` handle; `SpectralDataset` threads the references group into each
`GenomicRun` at open time. The `_native_h5py` escape shim — its last caller —
is removed. No REF_DIFF wire/decode change: encoded→decoded sequences stay
byte-identical and the external-FASTA / `REF_PATH` / Q5c hard-error fallbacks
are preserved. (OO-assessment P3.9.)

### Changed — per-run provenance cold path reads via the StorageGroup protocol (Python)

The per-run provenance cold path (`AcquisitionRun.provenance` and
`GenomicRun.provenance_chain`) now reads the compound `provenance/steps` dataset
through the `StorageGroup` protocol (`has_child`/`open_group` navigation +
`read_compound_dataset`) instead of casting to a raw h5py handle via
`_native_h5py`. Decoded provenance is identical and on-disk `.tio` bytes are
unchanged. (OO-assessment P3.9.)

### Changed — signature sign/verify route raw-h5py through the protocol (Python)

`sign_dataset`/`verify_dataset` now wrap raw `h5py.Dataset` inputs in the HDF5
provider's `StorageDataset` adapter and route through the provider-agnostic
`sign_storage_dataset`/`verify_storage_dataset`, so the live signature path no
longer touches raw h5py directly. The deprecated v1 unprefixed native-bytes
verify is the sole remaining raw-h5py path (scheduled for separate removal at
v1.0). No signature/`.tio` wire change — signatures stay byte-identical and
cross-language conformance is preserved. (OO-assessment P3.9.)

### Changed — `SpectralDataset.file` raw-h5py handle removed (Python)

The public `SpectralDataset.file` attribute (a raw `h5py.File` handle) is
removed; all storage is now reached through `dataset.provider`, the
protocol-based `StorageProvider`. The remaining internal consumers (reference
embed `ReferenceImport.write_to_dataset`, the FASTA-export CLI helper, and the
subject/sample readers) were migrated to navigate via `provider.root_group()`
and the `StorageGroup` protocol. On-disk `.tio` bytes are unchanged. (OO-assessment
P3.9.)

### Changed — Java SqliteProvider JSON reader uses Jackson (Java)

The Java `SqliteProvider`'s hand-rolled compound-JSON reader (a brittle
string-splitter vulnerable to whitespace, key reordering, and values containing
JSON tokens — the class of bug behind #205) is replaced by Jackson
(`jackson-databind`, declared as a direct dependency; already on the classpath
via Arrow). The byte-canonical JSON *serializer* is unchanged, so cross-language
compound byte-parity is preserved; a cross-language compound round-trip
conformance test is added. No `.tio`/SQLite schema change. (OO-assessment P2.7.)

### Changed — Shared Image base + uniform image collection (Java + tio-browser)

`MSImage`/`RamanImage`/`IRImage` now share an abstract `Image` base (common
geometry, intensity cube, metadata) with an `ImageKind` discriminator and a
generic `spectralAxis()`. `SpectralDataset`'s `image()`/`ramanImage()`/
`irImage()` are replaced by `imageForKind(ImageKind)` + `images()`; the transport
writer/walker, exporter adapter, CLI, and tio-browser export eligibility
migrated. No `.tio` wire/format or transport-protocol change; image-class
constructor signatures unchanged. Second of the 3-SDK P2.5 ports (Python #219;
ObjC follows). (OO-assessment P2.5.)

### Changed — Shared Image base + uniform image collection (Objective-C)

`TTIOMSImage`/`TTIORamanImage`/`TTIOIRImage` now share a `TTIOImage` base (common
geometry, intensity cube, metadata) with a `TTIOImageKind` discriminator and a
generic `spectralAxis`. `TTIOSpectralDataset`'s `msImage`/`ramanImage`/`irImage`
are replaced by `-imageForKind:` + `-images`; the transport writer/walker,
exporter adapter, and CLI migrated. No `.tio` wire/format or transport-protocol
change; initializer signatures unchanged. Completes the 3-SDK P2.5 shared-image
abstraction (Python #219, Java #220). (OO-assessment P2.5.)

### Changed — Shared Image base + uniform image collection (Python)

`MSImage`/`RamanImage`/`IRImage` now share an `Image` base (common geometry,
intensity cube, dataset-level metadata) with an `ImageKind` discriminator and a
generic `spectral_axis`. `SpectralDataset`'s `image`/`raman_image`/`ir_image`
accessors are replaced by `images` (collection of present kinds) +
`image_for_kind(kind)`; every consumer (exporters, transport walker/writer, CLI)
migrated. No `.tio` wire/format or transport-protocol change; each image keeps
its own on-disk group + I/O. First of the 3-SDK P2.5 ports (Java + tio-browser
and ObjC follow). (OO-assessment P2.5.)

### Performance — Vectorized genomic region/flag queries (Java + Objective-C)

`GenomicIndex.indicesForRegion`/`indicesForFlag` (Java) and `TTIOGenomicIndex`
(Objective-C) now scan the interned `chromosome_ids` (uint16) with integer
comparisons instead of a per-read string compare, resolving the query
chromosome to its id once — finishing the parity with the Python
`indices_for_region` vectorization (PR #202). Identical returned indices; no
wire/format change. (OO-assessment P1.2.)

### Performance — Vectorized DELTA_RANS + cached signal-channels handle (Python)

`DELTA_RANS` encode/decode now compute delta + zigzag via numpy for
element_size 1 and 4 (only the variable-length varint stream stays serial;
element_size 8 stays scalar for byte-exact safety) — ~2.3× faster on a 1M-element
int32 channel. `GenomicRun` caches the `signal_channels` group handle instead of
re-opening it on every channel access. Byte-identical output; no wire/format
change. (OO-assessment P1.3 + P1.4.)

### Changed — Importer/exporter dispatch unified behind Reader/Writer interfaces (Python)

Python importers now implement a uniform `Reader` protocol returning an
`ImportedDataset` draft (the single `SpectralDataset.write_minimal` call site),
and exporters a uniform `Writer` protocol over an opened dataset. The
importer/exporter registries dispatch through these interfaces instead of
per-format adapter callables; the run-selection helpers are shared. No `.tio`
wire/on-disk change; supported formats, aliases, and the `ttio encode`/`export`
CLI are unchanged. First of the 3-SDK P2.6 parity ports (Java + ObjC follow).

### Changed — Importer/exporter dispatch unified behind Reader/Writer interfaces (Java SDK)

The Java SDK gains uniform `Reader`/`Writer` interfaces, an `ImportedDataset`
draft (the single dataset-write call site, supporting a write-through delegate
for the subprocess-backed Bruker importer), shared `RunSelection` export helpers,
per-format Reader/Writer adapters, `ImporterRegistry`/`ExporterRegistry`
mirroring the Python registries (11 import / 8 export formats; `fasta`/`fastq`
remain CLI-delegated), and unified `encode`/`export` CLI tools. The dataset
write path (`SpectralDataset.createMixed`) learned MS/Raman/IR image writing
(reusing `MSImage.writeTo`; image-free output is byte-identical). No `.tio` wire
change; format keys/aliases match Python (asserted by a cross-language parity
test). JCAMP-DX / NMR export are scoped to the selected run's first spectrum
(matching the Python contract). The tio-browser GUI dispatch migration is the
next PR. Second of the 3-SDK P2.6 importer/exporter ports (Python shipped in
#213; ObjC follows).

### Changed — tio-browser dispatches imports/exports via the SDK registry (Java)

The desktop GUI's `ImportTask`/`ExportTask` now dispatch the registry-covered
formats (11 import / 8 export) through the Java SDK `ImporterRegistry`/
`ExporterRegistry` and shared `RunSelection`, replacing the per-format
`importX`/`exportX` bodies and the duplicated `toWritten`/`pickRun`/
`pickGenomicRun`; the GUI format registries now source extensions and
required-tool from the SDK registry. `fasta`/`fastq` (and the FASTA
reference/reads + FASTQ export rows) remain GUI-local (CLI-delegated in the
SDK). Two-phase import progress is preserved; registry-dispatched **exports**
no longer emit per-section writer-phase progress (the bar advances 50%→100% on
completion) — fasta/fastq exports retain full progress. The SDK registry's
genomic `requiredTool` was corrected to null (Java uses bundled htsjdk, not
samtools). No `.tio` wire change. Completes the P2.6 importer/exporter parity
for Java (SDK in #214); ObjC port follows.

### Changed — Importer/exporter dispatch unified behind Reader/Writer protocols (Objective-C)

The ObjC SDK gains `TTIOReader`/`TTIOWriter` protocols, a `TTIOImportedDataset`
draft (single write site, with a write-through delegate for the subprocess-backed
Bruker importer, for imzML image output via `TTIOMSImage`, and for the
dataset-returning XML importers), `TTIORunSelection` helpers (incl. the read-side
genomic→written conversion), per-format Reader/Writer adapters,
`TTIOImporterRegistry`/`TTIOExporterRegistry` mirroring the Python registries
(11 import / 8 export; `fasta`/`fastq` CLI-delegated; genomic `requiredTool`
="samtools" per ObjC's samtools-subprocess reality), a cross-language registry
parity test, and `TtioEncode`/`TtioExport` CLI tools. imzML import now produces
a `.tio` (previously parse-only). No `.tio` wire change. Completes the 3-SDK
P2.6 importer/exporter parity (Python #213, Java #214/#215).

## [1.6.5] - 2026-06-03

Refactor-only release completing the **3-SDK codec-registry parity**. Genomic
codec dispatch in all three reference SDKs — Python, Java, and Objective-C — now
routes through a single `Codec` registry keyed by the `Compression` id, fronted
by a uniform codec interface, a `CodecContext` value object, and closed
`DecodedChannel`/`EncodedChannel`/`ChannelPayload` unions, replacing the former
decode ladders, encode switches, and per-codec side-paths. Adding a genomic
codec is now one registry entry per language. Codecs expose two distinct flags,
`isContextAware` (needs run context) and `needsEmbeddedReference` (the
reference-embed predicate, REF_DIFF_V2 only) — the latter is parity metadata and
does not alter the embed decision. No wire/on-disk format change; every
byte-equality and cross-language fixture is unchanged.

PRs: #209 (Python registry) · #210 (Java port) · #211 (Objective-C port).

### Changed — Codec dispatch unified behind a registry (Objective-C)

ObjC's genomic codec dispatch (the decode switch + five bespoke
ref_diff/fqzcomp/name_tok/mate_info/cigars side-paths and the two near-identical
encode bodies) now routes through a single `TTIOCodec` registry keyed by
`TTIOCompression` (`Codecs/Registry/`), fronted by a uniform `TTIOCodec`
protocol, a `TTIOCodecContext` value object, and abstract-class-cluster unions
`TTIODecodedChannel`/`TTIOEncodedChannel`/`TTIOChannelPayload`. Codecs expose
`isContextAware` and `needsEmbeddedReference` (REF_DIFF_V2 only). No
wire/on-disk format change; all byte-equality and cross-language fixtures
unchanged. Completes the 3-SDK codec-registry parity (Python #209, Java #210).

### Changed — Codec dispatch unified behind a registry (Java)

Java's genomic codec dispatch (the decode ladder + four bespoke
ref_diff/fqzcomp/name_tok/mate_info side-paths and the encode switch + writer
methods) is replaced by a single `Codec` registry keyed by `Compression`
(`codecs/registry/`), fronted by a uniform `Codec` interface, a `CodecContext`
value object, and sealed `DecodedChannel`/`EncodedChannel`/`ChannelPayload`
unions. Codecs expose `isContextAware` and `needsEmbeddedReference` (REF_DIFF_V2
only); the duplicated `useRefDiffPath` embed predicate is consolidated into one
helper. `DELTA_RANS_ORDER0` is now registered (previously unwired). No
wire/on-disk format change; all byte-equality and cross-language fixtures
unchanged. Mirrors the Python registry (PR #209); ObjC port is the remaining
parity follow-on.

### Changed — Codec dispatch unified behind a registry (Python)

Python's genomic codec dispatch (the decode/encode `if`/`elif` ladders, the
four bespoke ref_diff/fqzcomp/name_tok/mate_info side-paths, and the
`_codec_meta` context-aware set) is replaced by a single `Codec` registry keyed
by `Compression` id, fronted by a uniform `Codec` interface, a `CodecContext`
value object, and closed `DecodedChannel`/`EncodedChannel` unions
(`ttio/codecs/_registry.py`, `_context.py`). Adding a codec is now one registry
entry. Codecs expose `is_context_aware` (needs run context) and
`needs_embedded_reference` (the reference-embed predicate, formerly
`_codec_meta._CONTEXT_AWARE`). No wire/on-disk format change; all byte-equality
and cross-language fixtures unchanged. Java and ObjC parity ports are tracked as
follow-on work.

## [1.6.4] - 2026-06-02

Patch release bundling two cross-language SDK fixes. (1) The Java
`SqliteProvider` could not read `.tio.sqlite` compound datasets written by
the Python provider — a whitespace-intolerant JSON parser broke the
Python→Java direction; fixed, with a new cross-language round-trip test.
(2) Reference-embed progress reached cross-language parity: Python
`ReferenceImport.write_to_dataset` and ObjC `TTIOReferenceImport` gained
per-contig `progress` support matching Java, the three stale Java importer
parity TODOs were removed, and the previously-untested ObjC `TTIOBamReader`
progress sink got a test. No wire/on-disk format change.

PRs: #205 (Java SQLite cross-language read) · #206 (reference-embed progress parity).

### Added — Reference-embed progress parity (Python + ObjC)

`ReferenceImport.write_to_dataset` (Python) gains a keyword-only
`progress` parameter, and `TTIOReferenceImport` (ObjC) gains a
`-writeToDataset:overwrite:progress:error:` overload, both firing
per-contig progress (`(0, N)` then `(i+1, N)` ending at `(N, N)`) to
match Java's `ReferenceImport.writeToDataset(..., ProgressSink)`
(`@since 1.3.0`). Closes the last cross-language gap in the importer
progress-sink surface; the three stale "future parity PR" TODOs in the
Java importers are removed (the reader sinks they described already
shipped in all three SDKs). Added an ObjC test for the existing
`TTIOBamReader` progress sink (was untested). No wire/on-disk format
change. New tests: `test_reference_import_progress.py` (Python),
`writeToDataset…progress:` + BAM-reader cases (ObjC).

### Fixed — Java SQLite provider couldn't read Python-written compound datasets

The Java `SqliteProvider` failed to open a `.tio.sqlite` compound dataset
written by the Python `ttio.providers.sqlite` provider: its hand-rolled JSON
parser for compound field descriptors (`extractJsonString`) searched for the
literal token `"key":"` with no whitespace, but Python's `json.dumps` emits a
space after the colon (`"key": "value"`). The lookup returned `null`, and
reading the field kind threw `NullPointerException`. Java read its own
(whitespace-free) output fine and Python reads Java's via standard
`json.loads`, so only the Python→Java direction was broken — a silent
cross-language compatibility gap, since no test exercised it. `extractJsonString`
now tolerates optional whitespace around the colon (spec-valid JSON is
whitespace-insensitive). New regression test
`SqliteProviderTest::crossLanguagePythonWrittenFileReadback` writes a fixture
via the Python provider (`python/tests/fixtures/make_sqlite_fixture.py`) and
reads it back through the Java provider; it is skipped when the Python `ttio`
package isn't importable and runs in the cross-language parity CI job.

## [1.6.3] - 2026-06-02

Patch release bundling the post-v1.6.2 Python fixes and FD-1 client SDK
additions. Three threads: (1) the encrypted-reader per-spectrum metadata
data-loss fix and 3-language parity restoration (#199); (2) the
wrap-for-server / unwrap-for-server client SDK helpers + `ServerRecipient`
and `download_via_server`, completing the HSM server-side key-custody
round-trip from the client side (FD-1-PF-6, FD-1-PF-7; paired with
`tti-workbench-server` #77/#78); (3) a ~17× speedup of `GenomicIndex`
region lookups via the interned chromosome id column. No on-disk format
change; the Java and ObjC reference implementations are unchanged.

PRs: #200 (#199) · #197 (FD-1-PF-6) · #198 (FD-1-PF-7) · #202 (perf).

### Changed — Performance: vectorized `GenomicIndex.indices_for_region` (Python)

`GenomicIndex.indices_for_region` (Python) scanned the per-read
`chromosomes` Python list on every call — an O(N) name comparison that
dominated region queries on deep WGS indexes. The interned `chromosome_ids`
(uint16) + `chromosome_names` lookup table is now retained on the dataclass
when an index is loaded from disk via `GenomicIndex.read`, and the region
lookup resolves the chromosome name to its id once (O(unique)) then compares
the id column with a vectorized NumPy operation. In-memory construction
(without an id table) keeps the original list-scan fallback, so behaviour is
unchanged for every existing caller. Measured ~17× faster on a 2 M-read
index; results are byte-identical to the previous path. No on-disk format
change. New fields `chromosome_ids` / `chromosome_names` (both default
`None`) are additive. Regression tests:
`test_m82_genomic_run.py::test_genomic_index_read_retains_interned_chromosomes`
and `::test_genomic_index_indices_for_region_disk_loaded`.

### Fixed — TTI-O#199: Python encrypted reader dropped per-spectrum metadata

`read_encrypted_to_file` (Python) discarded the MS per-spectrum metadata
columns — `retention_times`, `ms_levels`, `polarities`, `precursor_mzs`,
`precursor_charges`, `base_peak_intensities` — when reconstructing a `.tio`
from an encrypted transport stream produced with the default
`encrypt_headers=False`. The data was present on the wire (in the plaintext
AU filter header) but `_ingest_encrypted_au` decoded and then dropped it, and
the `spectrum_index` reconstruction wrote only `offsets`/`lengths`. The Java
and ObjC readers already reconstruct these columns, so this was a Python-only
data-loss bug **and** a 3-language parity violation: every encrypted download
path (`download_via_server`, `download_decrypted`, `download_decrypted_envelope`,
`download_decrypted_pqc`, `download_decrypted_multi`) silently lost per-spectrum
metadata. The reader now accumulates the AU header fields (reversing the wire
polarity encoding) and writes them back, mirroring the Java/ObjC readers.
Regression test: `test_encrypted_transport.py::TestEncryptedChannelRoundTrip::
test_full_roundtrip_preserves_spectrum_metadata`.

### Added — FD-1-PF-6: wrap-for-server client SDK + `ServerRecipient`

Python SDK half of FD-1-PF-6. Complements `tti-workbench-server` #77
(FD-1-PF-4), which ships the daemon-side `/v1/key-custody/wrap-for-server`
endpoint. Together they allow the client to ask the daemon to wrap a DEK under
an HSM-resident PKCS#11 key without the KEK ever leaving the server.

New public surface (`python/src/ttio/workbench/client.py`):

- **`ServerRecipient(recipient_id, kek_id, algorithm="aes-256-gcm")`** —
  frozen dataclass mirroring `EnvelopeRecipient`, but the caller holds no KEK
  bytes; the daemon wraps on its behalf.
- **`WorkbenchClient.wrap_for_server(*, dek, kek_id) -> bytes`** — thin async
  method that POSTs `{"dek": <b64>, "kek_id": <str>}` to
  `/v1/key-custody/wrap-for-server` and returns the wrapped blob. Raises
  `ValueError` for non-32-byte DEK; raises `WorkbenchHttpError` on 401/403/404/409.
- **`upload_encrypted_multi`** now accepts `list[EnvelopeRecipient |
  ServerRecipient]`. When a `ServerRecipient` is present it calls
  `wrap_for_server` instead of `_wrap_dek`, and auto-stamps `server_kek_id`
  from the first `ServerRecipient`'s `kek_id`. Passing both explicitly with
  inconsistent values raises `ValueError`. Pure-`EnvelopeRecipient` callers are
  unaffected.

Both `EnvelopeRecipient` and `ServerRecipient` are re-exported from
`ttio.workbench`. 11 new unit tests in
`python/tests/workbench/test_fd1_pf6_wrap_for_server.py` (daemon mocked via
`patch("ttio.workbench._http.http_json")`).

### Added -- FD-1-PF-7: unwrap-for-server SDK helper + download_via_server (Python)

Python SDK complement to `tti-workbench-server` #81 (FD-1-PF-7 unwrap
endpoint). Closes the FD-1 client round-trip started by #80 (PF-6 SDK upload)
and #77 (PF-4 daemon wrap): callers can now download and decrypt containers
that were uploaded via `ServerRecipient` without ever holding KEK bytes.

New public surface (`python/src/ttio/workbench/client.py`):

- **`WorkbenchClient.unwrap_for_server(*, wrapped_dek, kek_id) -> bytes`** --
  inverse of `wrap_for_server`. POSTs the wrapped blob to
  `/v1/key-custody/unwrap-for-server`; daemon returns the 32-byte DEK.
  Raises `WorkbenchHttpError` on 401/403 (auth), 404 (unknown `kek_id`), or
  422 (AEAD authentication failure / tampered blob). Raises `ValueError` if
  the daemon returns a non-32-byte DEK (defensive guard).
- **`WorkbenchClient.download_via_server(*, container_uri, out_tio_path,
  filters=None, max_au=0)`** -- full server-mediated download. Downloads the
  `.tis`, materialises the encrypted `.tio` to disk, reads `wrapped_dek` and
  `server_kek_id` from `ProtectionMetadata`, calls `unwrap_for_server` to
  recover the DEK, then performs per-AU AES-256-GCM decryption client-side.
  Raises `ValueError` with explicit guidance ("use download_decrypted_multi
  instead") when the container has no server-recipient.

6 new unit tests extend `python/tests/workbench/test_fd1_pf6_wrap_for_server.py`
(daemon mocked via `patch("ttio.workbench._http.http_json")`); total in that
file is now 17.  The happy-path test exercises the full download round-trip by
uploading with a parallel `EnvelopeRecipient` to recover the actual DEK for
mock injection.


## [1.6.2] - 2026-05-28

Documentation-only patch. Five-pass per-instance-method documentation completion across all three SDKs. Class-level OpenStep-style headers shipped in May 2026; this release fills the deferred per-method `/** */` blocks (ObjC), NumPy docstrings (Python), and Javadoc (Java) so the generated `autogsdoc` / Sphinx / Javadoc output now provides robust SDK reference for every public surface.

~110 files touched, ~5892 lines of structured documentation added. No logic changes; tests green throughout.

Coverage by pass:

| Pass | PR | Surface |
|---|---|---|
| 1 | #191 | Python transport + workbench client + providers/hdf5; ObjC Transport core + Providers |
| 2 | #192 | Workbench subsurface (jobs/sessions/containers/pipeline/cohort/auth/proxy); ObjC Transport Packet/AccessUnit + Import readers |
| 3 | #193 | Importers + exporters; ObjC Protection, Genomics, Run, Codecs, CVTermMapper, ArrowIpcCodec |
| 4 | #194 | Python core (spectral_dataset, signal_array, aligned_read); ObjC Dataset; **Java top-level 187 methods** (earlier audit measurement was broken) |
| 5 | #195 | Java fringe subpackages (protocols, io, hdf5, analysis, providers); Python CLI tools |

Drive-by milestone-reference scrub on touched files (M-numbers, `@since`, version-history prose) per DOC-AUDIT.md rules. Restored one `v1.0` reference in `cmd_provenance` runtime error message — that version reference is part of the user-facing contract a test exercises.

## [1.6.1] - 2026-05-27

Patch release. Closes the v0.11 silent-drop pattern in the two non-daemon transport servers, surfaced by writing end-to-end round-trip tests against the new live coverage from v1.6.0.

- **#144 — Python `transport.serving` framing.** `_emit_stream` packed every multi-packet writer call (`write_reference_group` = HEADER + N × CHROMOSOME + FOOTER, `write_image` = HEADER + N × PIXEL + FOOTER, etc.) into a single WebSocket frame, but the client parses one packet per frame. Every packet after the first of each multi-packet emitter was silently dropped. Fix: per-packet reframing on the server side.
- **#145 — Java `TransportServer` walker delegation.** The Java reference server hand-rolled its emission loop and walked only `msRuns()`, so every v0.11 prelude accessor + every genomic AU was silently dropped on the wire — same shape as the now-fixed ObjC daemon walker (#140) and Python framing bug (#144), in a third location. Fix: delegate to `DatasetWalker` via a `WriterDispatchVisitor` with the same per-packet reframing as the Python fix.

Both fixes are validated by new end-to-end round-trip tests (`TestV011RoundTrip` in Python, `v011FullAccessorRoundTrip` in Java) mirroring the workbench-live coverage from v1.6.0.

PRs: #188 #144 · #189 #145.

Paired server CI tightening (in tti-workbench-server #61): server CI now runs the TTI-O Python live-smoke against the built daemon, so a server PR that breaks the WS visitor / handshake / framing fails its own CI instead of waiting for a TTI-O PR to surface it.

### Fixed — Java `TransportServer` walks the dataset (#145)

`TransportServer.streamDataset` hand-rolled its emission loop and walked only `dataset.msRuns()`, so every v0.11 prelude accessor (references, subjects, samples, identifications, quantifications, image cubes, dataset_provenance, encryption_algorithm) plus all genomic AUs were silently dropped on the wire — same root cause as the now-fixed ObjC daemon walker (#140) and the Python framing bug (#144), in a third location.

The server now delegates to `DatasetWalker` via a `WriterDispatchVisitor` that pipes every event through a per-call `TransportWriter` sinked at a `ByteArrayOutputStream`, then splits the writer's concatenated output (`writeReferenceGroup` etc. emit multiple packets per call) back into individual packets via `PacketHeader.decode` and sends each as its own WebSocket binary frame — same per-packet reframing the Python fix introduced.

New `TransportServerTest.v011FullAccessorRoundTrip` builds `FixtureBuilder.buildEverything`, serves it, downloads via `TransportClient.streamToFile`, and asserts every `AccessorSpec` round-trips byte-equivalent. Mirrors the Python `TestV011RoundTrip` (#144) and the ObjC daemon's `test_v011_full_accessor_round_trip` (#140).

Java `AccessorSpec.RAMAN_IMAGE` / `IR_IMAGE` comparators now treat "both null" as trivially equal — `buildEverything` populates only the MS image cube, and the strict-non-null check was tripping immediately on otherwise-correct round-trips. Matches the Python comparator polish from the v1.6.0 release.

### Fixed — Python transport server frames v0.11 prelude packets individually (#144)

`_emit_stream` invoked `TransportWriter.write_reference_group` / `write_image` / `write_subject_metadata` / etc. into a single `BytesIO`, then sent the concatenated bytes as one WebSocket binary frame. Each of those writer methods emits **multiple** packets (e.g. `REFERENCE_GROUP_HEADER + N × REFERENCE_CHROMOSOME + END_OF_REFERENCE_GROUP`). `TransportClient._split_packet` only parses one packet per frame, so every chromosome / pixel / subject row after the first packet of each multi-packet emitter was silently dropped on download.

Fix: split the writer's buffer back into per-packet slices on the server side and send each as its own binary frame.

Added `tests/test_transport_server.py::TestV011RoundTrip` — builds `v0_11_fixtures.build_everything`, serves via `serving()`, downloads via `TransportClient.stream_to_file`, asserts every accessor round-trips byte-equivalent through `ACCESSOR_SPECS`. Mirrors the workbench-live coverage for the Python reference server path.

## [1.6.0] - 2026-05-27

Patch release closing the v0.11 transport round-trip end-to-end through the workbench daemon. v1.5.0 shipped the v0.11 spec + per-SDK encode/decode parity; v1.6.0 wires the missing pieces so a multi-accessor `.tio` uploaded to the daemon and downloaded back is byte-equivalent across every populated accessor.

Theme: three layered bugs surfaced once a single live-smoke test (`test_v011_full_accessor_round_trip`) exercised the end-to-end path for the first time.

1. **AU sequence ingest (#139)** — TransportIngest enforced stream-wide monotonicity, but walkers emit per-dataset sequences. Fix: per-`datasetId` tracking in all three ingesters. No wire-format change.
2. **Daemon download walker (#140)** — `TTIODatasetWalker` and the WS download visitor only spoke MS access units; every v0.11 prelude packet + every genomic AU was silently dropped on the way back to the client. Fix: 10 new visitor methods + §5.4 prelude emission + per-read genomic AU iteration. Paired daemon-side fix shipped in tti-workbench-server #59.
3. **Python + Java walker parity (#141)** — same bug shape as #140 in the two non-daemon walkers, so `ttio.transport.serving` (Python's reference WS server) had the same silent-drop. Fix: 10 new event types / visitor selectors in each SDK + shared `_iter_genomic_run_access_units` helper for byte-form equivalence with the writer.

tio-browser polish: PR #183 adds a `Size` column to the Transfers workspace and a live MiB-counter on download progress (#135, #136).

PRs in 1.6.0 (in order): #181 v0.11 round-trip live-smoke test · #183 tio-browser Download progress + Size column · #184 #139 · #185 #140 · #186 #141.

### Fixed — Python + Java walker v0.11 parity (#141)

Follow-up to #140 (ObjC walker). The Python `walk_dataset` generator and Java `DatasetWalker` previously emitted only MS access units, so every workbench-daemon download routed through the Python/Java reference servers silently dropped v0.11 accessor content — same root cause as #140 but in the other two SDKs.

Python: `walker.py` gains `EncryptionAlgorithmEvent`, `DatasetProvenanceEvent`, `SubjectMetadataEvent`, `SampleMetadataEvent`, `ReferenceGroupEvent`, `ImageEvent`, `RamanImageEvent`, `IRImageEvent`, `IdentificationsTableEvent`, `QuantificationsTableEvent`. The v0.11 §5.4 prelude (encryption → provenance → subjects → samples → references → image → identifications → quantifications) is yielded between `StreamHeaderEvent` and the first `DatasetHeaderEvent`, gated on each accessor being populated. Genomic AUs now iterate via a shared `_iter_genomic_run_access_units(run)` helper extracted from `_emit_genomic_run_access_units`. `server._emit_stream` dispatches each new event through the matching `TransportWriter.write_*` method.

Java: `AccessUnitVisitor` gains 10 optional default-method callbacks (`visitEncryptionAlgorithm`, `visitDatasetProvenance`, `visitSubjectMetadata`, `visitSampleMetadata`, `visitReferenceGroup`, `visitImage`, `visitRamanImage`, `visitIRImage`, `visitIdentificationsTable`, `visitQuantificationsTable`). `DatasetWalker.walk` emits the same §5.4 prelude block plus per-read genomic AUs via a shared `TransportWriter.genomicRunAccessUnits(run)` helper. AU + EOD events now interleave per-dataset (matching `TransportWriter.writeDataset`).

### Fixed — multi-accessor AU sequence ingest (#139)

Walker emits AU sequences that reset per dataset (`for j, spectrum in enumerate(run)`), but TransportIngest in all three SDKs enforced a single stream-wide monotonicity counter. v0.11 multi-accessor `.tio` uploads (the first to exercise more than one dataset over the workbench daemon) hit `AU sequence regressed: got 0, last seen N` mid-upload on the second accessor.

Per-dataset tracking in the ingester. `PacketHeader.datasetId` is already on the wire, so no format change — each dataset's monotonicity is checked against its own last-seen sequence (`Map<datasetId, lastSeq>` / `dict[int, int]` / `NSMutableDictionary`). Within-dataset regression still fails loudly.

### Fixed — daemon download emits v0.11 accessors (#140)

Pairs with #139 to close the full v0.11 round-trip end-to-end. `TTIODatasetWalker` previously walked only MS access units; the workbench daemon's download path re-encodes through that walker, so every reference, subject, sample, identifications, quantifications, image, encryption-algorithm, and `dataset_provenance` row was silently dropped on the way back to the client. Genomic AU emission was also gated off behind an "intentionally not yet wired" comment.

Walker now emits the v0.11 §5.4 prelude (encryption → provenance → subjects → samples → references → image → identifications → quantifications) between the StreamHeader and the first DatasetHeader, and iterates each genomic run's reads as five-channel AccessUnits with the GenomicRead suffix (chromosome / position / mappingQuality / flags / matePosition / templateLength). `TTIOTransportEventVisitor` gained the matching optional `walker:visit*:` methods. The workbench daemon's `TTIOWBDownloadWriterVisitor` + `TTIOWBDownloadComboVisitor` implement each new event by dispatching to the corresponding `TTIOTransportWriter` method, with structured logging on per-event failure so a single bad accessor doesn't abort the whole download. `test_v011_full_accessor_round_trip` is now an active gate (xfail-strict removed).

## [1.5.0] - 2026-05-27

Major SDK release. Three themes:

1. **Transport-spec v0.11** — complete `.tio` accessor coverage across Java/Python/ObjC (was: silent-drop on every accessor besides MS/genomic runs). 12 new packet types (0x10–0x1B), Subject + Sample as first-class TTI-O entities, IRImage promoted to first-class on SpectralDataset.
2. **Progress reporting end-to-end** — `ProgressSink` callbacks throughout the SDK reader + writer chain in all three languages. tio-browser's ImportTask + ExportTask now show quantitative percent + ETA instead of an indeterminate spinner. PhaseProgress splits read (0..50%) and write (50..100%) cleanly.
3. **Streaming upload** — `WorkbenchTransportClient.upload(Path)` + `UploadClient.upload_path` replace the `Files.readAllBytes` slurp; peak heap is now O(chunkSize) instead of O(payload size). 100 GB+ uploads no longer OOM client heap.

PRs merged into 1.5.0 (in order): #172 transport v0.11 · #173 TransferRow refresh fix · #174 Stage B (per-read reader progress) · #175 Stage C (per-spectrum reader progress) · #176 Stage D (writer progress + per-section) · #177 Stage E (two-phase ImportTask/ExportTask + heartbeat removal) · #178 streaming upload Java · #179 Python parity · #180 ObjC parity.

7 stage tags pinned on the v0.11 milestone history (`stage-0-transport-v0-11-foundation` through `stage-6-transport-v0-11-subjects-samples`).

### Added — Transport-spec v0.11: complete `.tio` coverage across all 3 SDKs (2026-05-26)

The transport (`.tis`) protocol now round-trips every first-class
`SpectralDataset` accessor in Java / Python / Objective-C. Resolves
the silent-drop bug where `writeDataset(...)` on a reference-bearing
`.tio` was producing a 180-byte `.tis` (only `StreamHeader` +
`EndOfStream`) because the writer iterated `msRuns` + `genomicRuns`
only.

**Wire format — 12 new packet types (0x10–0x1B)** documented in
`docs/transport-spec.md` §4.13–§4.23 and `docs/format-spec.md` §11
(new):

- `REFERENCE_GROUP_HEADER` 0x10, `REFERENCE_CHROMOSOME` 0x11,
  `END_OF_REFERENCE_GROUP` 0x12 — per-chromosome with encoding=0
  raw / encoding=1 zlib at a 4 KiB threshold.
- `IMAGE_HEADER` 0x13, `IMAGE_PIXEL` 0x14, `END_OF_IMAGE` 0x15 —
  modality dispatch byte (0=MS, 1=Raman, 2=IR) with a self-describing
  `u16 modality_extras_length + bytes` slot. MS supports continuous
  (dense) and processed (sparse `{channel,intensity}` pairs indexed
  into the shared `mzAxis`).
- `IDENTIFICATIONS_TABLE` 0x16, `QUANTIFICATIONS_TABLE` 0x17,
  `SUBJECT_METADATA` 0x19, `SAMPLE_METADATA` 0x1A — Arrow IPC stream
  with a `uint32 LE length` prefix.
- `DATASET_PROVENANCE` 0x18 — single packet with `uint32 record_count`
  + N inline records.
- `ENCRYPTION_ALGORITHM` 0x1B — `uint16 length + UTF-8 algorithm`.

Ordering follows transport-spec §5.4: encryption → provenance →
subjects → samples → references → image → identifications →
quantifications, BEFORE the v0.10 dataset/run sections.

**Backward compat.** `transport_v0_11` StreamHeader feature flag is
set only when v0.11 content is present. v0.10 streams are byte-
equivalent to pre-Stage-0 output. Readers ingest unknown packet types
via a skip-unknown forward-compat path (length-prefixed wire frames,
log + skip on unrecognized type byte).

**New first-class entities and accessors:**

- `Subject` / `Sample` — tight core (`external_id`/`project`/`sex`/
  `birth_year` and `sample_id`/`subject_external_id`/`sample_kind`/
  `collected_at`) plus an open `attributes: Map<String,String>` slot.
  Stored as per-row HDF5 groups under `/study/subjects/<external_id>/`
  and `/study/samples/<sample_id>/`. Validation: duplicate / empty /
  slash-bearing IDs raise; soft-FK mismatch (Sample referencing an
  unknown Subject) logs a WARNING. `AcquisitionRun.sampleName` stays
  the canonical run→sample link string (no breaking change). New
  `SpectralDataset.subjects()` / `.samples()` accessors (Java),
  `subjects` / `samples` properties (Python), `-subjects` / `-samples`
  (ObjC). New 9-arg `SpectralDataset.create(..., subjects, samples)`
  overload (Java); `subjects=` + `samples=` kwargs on Python's
  `write_minimal`.
- `IRImage` — promoted to first-class on SpectralDataset in all 3
  SDKs (`.irImage()` / `.ir_image` / `-irImage`) alongside the
  existing `.image()` (MSImage) and `.ramanImage()` accessors.
  IR-specific HDF5 attrs are now typed int64 (`ir_mode`) + float64
  (`resolution_cm_inv`, `pixel_size_*`) in all 3 SDKs for cross-lang
  byte parity; Java + ObjC readers retain backward compat for legacy
  VL-string attribute form.

**Conformance — three nets prevent future silent-drop regressions:**

- `AccessorMatrixConformanceTest` parameterised over 13 accessors
  in each SDK (REFERENCES, MS_RUNS, GENOMIC_RUNS, IMAGE, MS_IMAGE_
  PROCESSED, RAMAN_IMAGE, IR_IMAGE, IDENTIFICATIONS, QUANTIFICATIONS,
  DATASET_PROVENANCE, ENCRYPTION_ALGORITHM, SUBJECTS, SAMPLES). Each
  builds a single-accessor fixture, round-trips through writer +
  reader, asserts byte- or content-equivalence.
- `CoverageGapWatchdogTest` runs the `everything.tio` fixture
  (every accessor populated) through `writeDataset`; if `.tis < 1%`
  of `.tio` size, the watchdog fires. Caught the silent-drop bug
  immediately when run against pre-fix code.
- Cross-language conformance: `Tests/cross_lang/transport_v0_11/
  accessor_matrix_xlang.sh` drives every encoder × every decoder
  pair (9 directional pairs × 13 accessors = 117 cells, plus 3
  bytes-equal sanity tests = 120 cells). 111 pass / 9 skipped on
  envs without `libttio_rans` (only GENOMIC_RUNS gated).

**Cross-language byte equivalence:**

- Hand-written packets (REFERENCES, IMAGE, ENCRYPTION_ALGORITHM,
  DATASET_PROVENANCE) are byte-identical across Java / Python / ObjC
  when given identical inputs. Java's `ProvenanceRecord.parametersJson`
  was patched to sort keys for parity with Python's `sort_keys=True`
  and ObjC's `NSJSONWritingSortedKeys`.
- Arrow IPC packets (0x16, 0x17, 0x19, 0x1A) are logically
  equivalent — Arrow Java 16 / pyarrow 16 / libarrow-C++ 24 each
  emit different flatbuffer envelopes, but row content cross-decodes
  identically. Schema column names match exactly across SDKs.

**Tooling and build dependencies:**

- `pyarrow >= 16` added to `python/pyproject.toml` (required for
  Arrow IPC).
- Apache Arrow Java 16 added to `java/pom.xml`. Surefire `argLine`
  gains `--add-opens=java.base/java.nio=ALL-UNNAMED` for `arrow-
  memory-core`'s reflective `Buffer.address` lookup on JDK 17+.
- libarrow-C++ 24 wired into `objc/Source/GNUmakefile.preamble` via
  `pkg-config arrow` with `-std=c++20`. New `TTIOArrowIpcBridge.mm`
  (Objective-C++) keeps libarrow symbols quarantined to a single
  translation unit; pure-ObjC callers see `TTIOArrowIpcCodec` only.

**Tags:** `stage-0-transport-v0-11-foundation` through
`stage-6-transport-v0-11-subjects-samples` are pushed for incremental
review. 72 commits ahead of `main`.

### Added -- Genomic-run support in per-AU decrypt-in-place (all 3 languages) (2026-05-23)

Extends `decryptFilePathInPlace` / `decrypt_per_au_in_place` /
`decryptFileInPlace` (ObjC / Python / Java) to also walk
`study/genomic_runs/<name>/signal_channels/`, mirroring the encrypt-side
genomic loop. Per-AU GCM-decrypts each `<sequences|qualities>_segments`
compound (uint8) and writes the plaintext back as the bare
`<sequences|qualities>` dataset (no `_values` suffix, matching the
genomic layout). `dataset_id` continues from the MS loop, so the AAD
numbering matches the encrypt path exactly.

With genomic now handled, all three languages strip
`opt_per_au_encryption` / `opt_encrypted_au_headers` / root `@encrypted`
**unconditionally** after decryption (previously Python + Java gated the
strip on "no genomic_runs" -- that gate is removed here and ObjC's
already-unconditional strip is now correct because genomic IS handled).

Tests:
- ObjC: `testM90DecryptFilePathInPlaceGenomic` -- byte-equal restoration
  of sequences + qualities, segments + algorithm attrs removed, feature
  flag + @encrypted cleared.
- Python: `TestGenomicDecryptInPlace.test_genomic_round_trip` -- same.
- Java: `decryptInPlaceGenomicRoundTrip` in `M90GenomicProtectionTest`
  -- same (via re-encrypt + decryptFile sandwich to assert byte-equal).

### Added -- Python `decrypt_per_au_in_place` + Java `PerAUFile.decryptFileInPlace` (2026-05-23)

Cross-language mirrors of the ObjC `decryptFilePathInPlace:` API
added in the prior commit, closing the parity gap noted in that PR
(both languages previously had only the read-only `decrypt_per_au` /
`decryptFile` and the legacy single-dataset `decryptInPlace`, which
silently no-ops on per-AU files).

- **Python:** `ttio.encryption_per_au.decrypt_per_au_in_place(path, key)`.
  Same contract as ObjC: per-AU GCM decrypt → writes `<channel>_values`
  back, removes `<channel>_segments` and the `<channel>_algorithm`
  attribute, restores the 6 plaintext index columns when
  `opt_encrypted_au_headers` was set, then (when the file has no genomic
  runs) strips `opt_per_au_encryption` / `opt_encrypted_au_headers`
  feature flags and the root `@encrypted` attribute. Idempotent on
  plaintext files.
- **Java:** `global.thalion.ttio.protection.PerAUFile.decryptFileInPlace(path, key, providerName)`.
  Same contract.
- **Genomic-run / mixed-MS+genomic refinement:** all three languages
  now leave the feature flags intact when `genomic_runs` is present
  (those segments aren't yet handled by the MS-only decrypt-in-place
  path -- the genomic follow-up will lift this gate). The ObjC impl
  shipped in the prior commit will be updated to match in the genomic
  decrypt-in-place PR.

Covered by `TestDecryptPerAuInPlace` (Python) and `decryptInPlace*`
tests in `PerAUFileTest` (Java): channels-only round-trip with
byte-equal intensity recovery + feature-flag clearing, headers-mode
restoration of all 6 plaintext index columns, and idempotency on a
plaintext file.

### Added -- ObjC `+[TTIOPerAUFile decryptFilePathInPlace:withKey:providerName:error:]` (2026-05-23)

Persist-to-disk decrypt counterpart to `+encryptFilePath:` for the per-AU
compound layout. Mirrors the legacy
`+[TTIOSpectralDataset decryptInPlaceAtPath:withKey:error:]` (which only
handles the older `intensity_values_encrypted` single-dataset format and
is a silent idempotent no-op on per-AU files).

For each MS run with `<channel>_segments` under `signal_channels`,
decrypts each spectrum row with the per-AU GCM scheme
(AAD = `dataset_id || au_sequence || channel_name`), writes the
concatenated plaintext back as `<channel>_values` (float64), deletes
`<channel>_segments`, and removes the `<channel>_algorithm` attribute.
When the file carries `opt_encrypted_au_headers`, the six plaintext index
datasets (`retention_times`, `ms_levels`, `polarities`, `precursor_mzs`,
`precursor_charges`, `base_peak_intensities`) are restored from
`au_header_segments` and the encrypted compound is removed. Strips
`opt_per_au_encryption` / `opt_encrypted_au_headers` feature flags and
the root `@encrypted` attribute. Idempotent on already-plaintext files.

Unblocks the TTI-O Workbench server FD-1 D-1 / D-2b pipeline path: D-1's
`decryptStoredContainerAtPath:` previously called the legacy decrypt and
was a silent no-op on per-AU containers; the server pipeline ran on
still-encrypted data and Phase E's "round-trip" was a false positive
(documented in tti-workbench-server #41). Genomic-run per-AU
decrypt-in-place is a follow-up; Python / Java equivalents are also
gaps to mirror.

Covered by `testPerAUFileDecryptInPlace` (channels-only round-trip with
byte-equal intensity values + feature-flag clearing, headers-mode
restoration of all 6 plaintext index columns, idempotency on a
plaintext file).

### Added -- FD-1 Phase C-2a-4: server_kek_id cross-language conformance (2026-05-22)

Pins `server_kek_id` byte-parity across Python / Java / ObjC in the shared
`conformance/multi_recipient/` golden vectors.

- Two new golden vectors (`prot_server_kek_id_single`,
  `prot_server_kek_id_multi`) and a `trailing_hex` field (the full trailing
  section after the five §4.4 fields = recipient block + optional
  `server_kek_id`); `gen_vectors.py` produces them from the Python encoder.
- Python/Java pin `body_hex` (now including `server_kek_id`) and assert the
  decoded `server_kek_id`; ObjC pins `trailing_hex` via a new
  `ttioConformanceEncodeTrailing:serverKekId:`. All three assert the same
  committed hex, so byte-parity is transitive.

### Added -- FD-1 Phase C-2a: server_kek_id in ProtectionMetadata (ObjC) (2026-05-22)

ObjC mirror of the C-2a `server_kek_id` field, completing the three-language
trio (byte-compatible with Python/Java).

- `TTIOEncryptedTransport`: a shared `appendProtectionTrailing` helper at
  both (MS + genomic) emit sites appends the stored recipient block + the
  optional `<channel>_server_kek_id` (count=0 + kid for single-recipient
  server-processable; nothing for BYOK → byte-identical). `parseProtection`
  reads `server_kek_id` after the recipient block into
  `ProtectionMeta.serverKekId`; both materialize paths re-stamp it.
- Covered by `TestC2aServerKekId` (round-trip through write→read, and BYOK →
  no attribute).

### Added -- FD-1 Phase C-2a: server_kek_id in ProtectionMetadata (Java) (2026-05-22)

Java mirror of the C-2a Python `server_kek_id` field (byte-compatible).

- `EncryptedTransport.encodeProtection` gains a `serverKekId` overload that
  appends the field after the Phase A recipient block (trailing section
  emitted iff additional recipients OR server_kek_id; single-recipient
  server-processable emits `count=0 + server_kek_id`); the 4-arg overload
  delegates with null. `parseProtection` / `ProtectionMeta` expose
  `serverKekId`.
- `stampTransportWrappedDek` gains a 6-arg overload that stamps
  `<channel>_server_kek_id`; the materialize paths re-stamp it;
  `readTransportServerKekId` reads it back. `WorkbenchClient.uploadEncryptedMulti`
  grows a `serverKekId` overload.
- Covered by `ServerKekIdProtectionTest` (byte-identity when absent,
  round-trip with/without additional recipients, storage carriage, BYOK →
  null); existing multi-recipient + conformance tests unaffected.

### Added -- FD-1 Phase C-2a: server_kek_id in ProtectionMetadata (Python) (2026-05-22)

Append-only `server_kek_id` field in the `ProtectionMetadata` packet, so the
workbench daemon can record a container's server-resolvable `kek_id` at
upload and decide server-processability (FD-1 Phase C-2).

- `_emit_protection_metadata(..., server_kek_id=None)` appends the field
  after the Phase A recipient block; the trailing section is emitted iff
  there are additional recipients OR a `server_kek_id`, so pure
  BYOK/single-recipient packets stay byte-identical to transport-spec §4.4.
- `_decode_protection_metadata` returns `server_kek_id` (None when absent).
- `stamp_transport_wrapped_dek(..., server_kek_id=None)` stamps
  `<channel>_server_kek_id`; write/read carry it; `read_transport_server_kek_id`
  reads it back. `WorkbenchClient.upload_encrypted_multi` grows a
  `server_kek_id` arg.
- Covered by `tests/test_fd1_c2a_server_kek_id.py` (byte-identity when
  absent, round-trip with/without additional recipients, storage carriage,
  BYOK → None). Java/ObjC mirrors + conformance follow in C-2a sub-steps.

### Fixed -- ObjC encrypted-transport round-trip + CI test gating (2026-05-22)

The ObjC `make check` runner exited 0 even when tests failed, so CI was
green while several encrypted-transport tests were actually broken. Made
the runner gate, then fixed the real defects it exposed:

- **Test gating:** `objc/build.sh check` now fails the build if the
  gnustep Testing framework reports any failed test/set, so local runs and
  CI catch failures instead of passing silently.
- **Raw packet write:** `TTIOEncryptedTransport._writeRawPacketHeader:`
  reached `TTIOTransportWriter`'s removed `fileHandle`/`dataBuffer` ivars
  via KVC (a casualty of the sink refactor), throwing on every encrypted
  write. Now writes through the writer's `_sink` (works for all sink
  types).
- **Binary attributes:** the HDF5 provider rejected `NSData` attributes, so
  `<channel>_wrapped_dek` / `_wrapped_dek_recipients` were never persisted.
  `TTIOHDF5Group` now stores binary attributes as HDF5 `OPAQUE`
  (`setDataAttribute:` / `dataAttributeNamed:`), and the provider adapter
  delegates reads/writes to the group's type-dispatching path.
- **Double attributes:** the provider truncated every `NSNumber` through
  `longLongValue`, so floating-point attributes (e.g. Raman
  `integration_time_sec`) lost their fractional part. Double/float
  `NSNumber`s are now stored as `H5T_NATIVE_DOUBLE` and read back as
  doubles.
- Fixes `TestFD1MultiRecipient` (incl. a `stream:`→`fromStream:` selector
  typo) and `TestJcampVibrationalRoundTrip` (Raman); the full ObjC suite is
  green under the new gate.
- Environment-dependent tests now skip instead of failing the gate:
  `TestTransportClient` skips (not fails) when no Python transport server
  is available (the standalone ObjC CI job has no venv), and the
  `TestM94ZFqzcompPerf` throughput floor is advisory unless
  `TTIO_PERF_STRICT` is set (it flaked on shared CI runners ~16 MB/s).

### Added -- FD-1 Phase C-0: standalone ObjC key-wrap primitive (2026-05-22)

Precursor for the server key-custody seam (FD-1 Phase C). ObjC lacked a
fileless key-wrap primitive that Java (`EncryptionManager.wrapKey`) and
Python (`key_rotation._wrap_dek`) already expose.

- `+[TTIOKeyRotationManager wrapKey:withKEK:algorithm:error:]` and
  `+[... unwrapKey:withKEK:algorithm:error:]` — wrap/unwrap a 32-byte DEK
  under a KEK (`aes-256-gcm` symmetric key or `ml-kem-1024` public/private
  key) to/from the canonical v1.2 wrapped-key blob, with no HDF5 file.
  They reuse the existing self-contained blob logic, so output is
  byte-identical to the file-bound `-wrapDEK:` path and to the Java/Python
  primitives. The Phase C-1 software-KMS stub calls these to (un)wrap a
  DEK under a tenant KEK.
- Covered by `TestC0StandaloneKeyWrap` (aes-256-gcm round-trip, wrong-KEK
  auth failure, unsupported-algorithm rejection).

### Added -- FD-1 Phase B-2: multi-recipient envelope client API (Java) (2026-05-22)

Java mirror of the B-1 client API on `WorkbenchClient`, byte-compatible
with the Python side (both stamp the Phase A multi-recipient packet).

- `WorkbenchClient.EnvelopeRecipient(recipientId, key, algorithm)` — one
  recipient's wrapping key (`aes-256-gcm` symmetric KEK or `ml-kem-1024`
  public key); a 2-arg constructor defaults to `aes-256-gcm`.
- `WorkbenchClient.uploadEncryptedMulti(project, containerUri, tioPath,
  recipients, preview, encryptHeaders)` — one DEK wrapped once per
  recipient via `EncryptionManager.wrapKey`; `recipients.get(0)` is the
  packet primary (wire id `""`), the rest go in the Phase A trailing
  block. Preview-gated iff any recipient uses `ml-kem-1024`.
- `WorkbenchClient.downloadDecryptedMulti(containerUri, key, outTioPath,
  recipientId, preview)` — selects the recipient entry the caller holds a
  key for (`""` = primary) and unwraps it.
- Covered by `MultiRecipientClientTest` (daemon-free wrap → stamp → read →
  unwrap-per-recipient → decrypt) and the live smoke
  `WorkbenchLiveTest.multiRecipientUploadRoundTrip` (server-KEK +
  researcher-ML-KEM end-to-end through the actual client methods).

### Added -- FD-1 Phase B-1: multi-recipient envelope client API (Python) (2026-05-22)

Client-side API to encrypt a `.tio` for **multiple recipients** — the FD-1
output shape (wrap the per-run DEK for both a server KEK and the
researcher's key). Pure client-side; reuses the BYOK / envelope / PQC wrap
primitives.

- `EnvelopeRecipient(recipient_id, key, algorithm)` — one recipient's
  wrapping key (`aes-256-gcm` symmetric KEK or `ml-kem-1024` public key);
  re-exported from `ttio.workbench`.
- `WorkbenchClient.upload_encrypted_multi(..., recipients=[...])` —
  generates one DEK, `encrypt_per_au`, wraps it once per recipient, and
  stamps `recipients[0]` as the packet primary (wire id `""`) with the
  rest in the Phase A trailing block. Preview-gated iff any recipient uses
  `ml-kem-1024`.
- `WorkbenchClient.download_decrypted_multi(..., recipient_id="")` —
  selects the recipient entry the caller holds a key for and unwraps it
  (`""` = primary). The daemon never holds a key.
- Covered by `tests/workbench/test_multi_recipient_client.py` (daemon-free
  in-memory data plane), incl. the server-KEK + researcher-ML-KEM shape.

### Added -- FD-1 Phase A-4: multi-recipient cross-language conformance (2026-05-22)

Completes Phase A: golden byte vectors pin the multi-recipient
`ProtectionMetadata` wire format byte-identically across Python / Java /
ObjC.

- `conformance/multi_recipient/vectors.json` — five `prot_*` golden vectors
  (spec §6) generated from the Python reference encoder by `gen_vectors.py`;
  each carries the full protection-metadata `body_hex`, the trailing
  `recipient_block_hex`, and a frozen pre-Phase-A primary-recovery
  expectation.
- All three suites assert against the *same* committed hex, so byte-parity
  is transitive (Python == golden ∧ Java == golden ∧ ObjC == golden ⇒ all
  equal): Python `tests/conformance/test_multi_recipient_xlang.py`, Java
  `MultiRecipientXLangTest`, ObjC `TestMultiRecipientXLang`.
- The recipient-block codec was exposed for testing without a public API
  change: package-private in Java; an internal `(Conformance)` category in
  ObjC. Single-recipient vectors confirm no trailing block is emitted.

### Added -- FD-1 Phase A-3: multi-recipient ProtectionMetadata (ObjC) (2026-05-21)

ObjC mirror completing the Phase A wire-format trio across Python / Java /
ObjC.

- `TTIOEncryptedTransport`: `encodeRecipientBlock` / `decodeRecipientBlock`
  (byte-identical to the Python/Java block); `parseProtection` consumes
  `signature_algorithm` + `public_key` then decodes the trailing block
  into `ProtectionMeta.additionalRecipients`; the MS + genomic write paths
  append the stored block to the packet, and both materialize paths
  persist it as `<channel>_wrapped_dek_recipients`. Single-recipient
  packets stay byte-identical.
- Covered by `TestFD1MultiRecipient`. Cross-language conformance vectors
  (A-4) follow to pin byte-parity across all three.

### Added -- FD-1 Phase A-2: multi-recipient ProtectionMetadata (Java) (2026-05-21)

Java mirror of the A-1 multi-recipient `ProtectionMetadata` wire format,
using the same append-only layout (single-recipient byte-identical).

- `EncryptedTransport`: public `Recipient` record; `encodeProtection`
  appends the trailing recipient block; `parseProtection` decodes it into
  `ProtectionMeta.additionalRecipients`; the write + materialize paths
  carry a base64 `<channel>_wrapped_dek_recipients` storage attribute
  (per-language on-disk encoding need not match across languages — only
  the transport packet bytes do).
- `stampTransportWrappedDek` gains a multi-recipient overload;
  `readTransportRecipients` returns the full list (primary + additional);
  `readTransportWrappedDek` stays the single-recipient accessor.
- Covered by `MultiRecipientProtectionTest`. ObjC (A-3) and cross-language
  conformance vectors (A-4) follow.

### Added -- FD-1 Phase A-1: multi-recipient ProtectionMetadata (Python) (2026-05-21)

First implementation step of the FD-1 server-side-compute groundwork: the transport
`ProtectionMetadata` packet can now carry the per-run DEK wrapped for
**multiple recipients** (e.g. a server KEK + a researcher key), the
prerequisite for encrypted pipelines whose output is decryptable both
server-side and client-side.

- Append-only wire layout: the existing five §4.4 fields are the
  **primary** recipient; additional recipients follow in an optional
  trailing block emitted only when present, so single-recipient packets
  (BYOK / envelope / PQC) stay **byte-identical** to the prior format and
  current readers parse them unchanged.
- `_emit_protection_metadata` gains `additional_recipients`;
  `_decode_protection_metadata` returns a `recipients` list (primary at
  index 0). Storage carries the extra recipients on a
  `<channel>_wrapped_dek_recipients` run attribute via
  `write_encrypted_dataset` / `read_encrypted_to_file`.
- `stamp_transport_wrapped_dek(..., additional_recipients=...)` +
  new `read_transport_recipients(path)`; `read_transport_wrapped_dek`
  remains the single-recipient accessor.
- Java + ObjC parity for the packet (Phase A-2/A-3) and cross-language
  conformance vectors (A-4) follow.

### Added -- vibrational .tio materialization, ObjC parity (2026-05-21)

ObjC mirror completing the §3.1 cross-language trio. `TTIOAcquisitionRun`
now materializes `TTIOIRSpectrum` / `TTIORamanSpectrum` /
`TTIOUVVisSpectrum` and reads/writes the per-class run-group attributes
matching the Python and Java contract, so a vibrational `.tio` written
by any language reads everywhere.

- `spectrumAtIndex:` dispatches on `spectrum_class` to build the right
  vibrational subclass; the in-memory initializer captures the metadata
  from the first spectrum; `writeToGroup:name:error:` emits `ir_mode` /
  `ir_resolution_cm_inv` / `ir_number_of_scans`, `raman_*`,
  `uvvis_path_length_cm` (UV-Vis solvent via `@solvent`); both readers
  restore them. The lighter `readFromStorageGroup:` reader also now reads
  `@solvent` so UV-Vis round-trips through it too.
- `TestJcampVibrationalRoundTrip`: write → read → materialize for each
  vibrational type. §3.1 (JCAMP-DX vibrational round-trip) is now fully
  resolved across Python + Java + ObjC.

### Added -- vibrational .tio materialization, Java parity (2026-05-21)

Java mirror of the §3.1 vibrational round-trip (PR after the Python
landing): `AcquisitionRun` now materializes `IRSpectrum` /
`RamanSpectrum` / `UVVisSpectrum` and persists/reads their metadata as
run-group attributes, matching the Python attribute contract so a
Python-written vibrational `.tio` reads in Java (and vice versa).

- `AcquisitionRun.setIRMetadata` / `setRamanMetadata` / `setUVVisMetadata`
  tag a run with a vibrational `@spectrum_class` + its scalar metadata;
  `objectAtIndex` dispatches on that class to build the right subclass
  (channels `wavenumber`/`intensity` or `wavelength`/`absorbance`).
- `writeTo` emits the same attributes Python writes (`ir_mode` /
  `ir_resolution_cm_inv` / `ir_number_of_scans`, `raman_*`,
  `uvvis_path_length_cm`; UV-Vis solvent via the existing `@solvent`),
  and `readFrom` restores them from `@spectrum_class` + the attrs.
  MS/NMR runs are unaffected (override unset).
- `JcampVibrationalRoundTripTest`: write → read → materialize for each
  type plus an MS-unaffected guard. ObjC parity follows.

### Added -- JCAMP-DX vibrational .tio round-trip (Python) (2026-05-21)

Closes parity-audit v1.0 §3.1 on the Python side: the `.tio` container
now round-trips the vibrational spectrum types (IR / Raman / UV-Vis),
unblocking `ttio encode --format jcamp-dx` and `ttio export --format
jcamp-dx`.

- `AcquisitionRun._materialize_spectrum` reconstructs `IRSpectrum`,
  `RamanSpectrum`, and `UVVisSpectrum` (previously only `MassSpectrum` /
  `NMRSpectrum`).
- Per-class metadata persists as scalar run-group attributes
  (`ir_mode` / `ir_resolution_cm_inv` / `ir_number_of_scans`,
  `raman_excitation_wavelength_nm` / `_laser_power_mw` /
  `_integration_time_sec`, `uvvis_path_length_cm`; UV-Vis solvent reuses
  the existing `solvent` attribute), via new `WrittenRun` fields written
  in `_write_run` and read in `AcquisitionRun.open`. MS/NMR files stay
  byte-identical (attributes emitted only when set). Added
  `_hdf5_io.write_float_attr` / `read_float_attr`.
- `jcamp-dx` removed from `DEFERRED_PYTHON` in both the importer and
  exporter registries; `importers.jcamp_dx.build_written_run` bridges a
  parsed vibrational `Spectrum` to a one-spectrum `WrittenRun`. CLI
  unsupported-format messages no longer special-case JCAMP-DX.
- Cross-language Java + ObjC parity for the vibrational `.tio`
  materialization follows in separate PRs.

### Removed -- workbench blob-level encryption path (2026-05-21)

Removed the daemon-incompatible W6.2/W6.3 blob upload path, superseded by
the per-AU encrypted upload (BYOK / envelope / PQC). It had no remaining
caller after the rework.

- Python: deleted `ttio/workbench/encryption.py` (`seal` / `open_sealed`
  / `ProtectionMode` / `ProtectionMetadata`) and the blob helpers in
  `ttio/workbench/pqc.py` (`seal_pqc` / `open_pqc` / `verify_pqc` /
  `PqcSealed` / `sig_keygen`); dropped the `WorkbenchClient.upload_protected`
  / `download_and_open` stubs. `ttio.workbench.pqc` keeps `kem_keygen`,
  `ML_KEM_1024`, and the `opt_pqc_preview` gate used by the per-AU path.
- Java: deleted the `global.thalion.ttio.workbench.encryption` package
  (`WorkbenchEncryptor` / `ProtectionMode` / `ProtectionMetadata`) and the
  blob methods on `WorkbenchPqc` (`sealPqc` / `openPqc` / `verifyPqc` /
  `PqcSealed` / `sigKeygen`); dropped `WorkbenchClient.uploadProtected`
  / `downloadAndOpen`. The transport-layer `ProtectionMetadata` (the real
  cross-language wire anchor) is unaffected.
- Trimmed the blob unit tests; the preview gate + keypair generator are
  still covered (unit + live).

### Changed -- workbench encrypted-upload follow-ups: docs + live coverage (2026-05-21)

Post-rework cleanup of the per-AU encrypted-upload feature:

- Fixed the `workbench_quickstart.ipynb` tutorial, which still pointed
  at the removed `upload_protected` / `download_and_open`; it now shows
  `upload_encrypted` / `download_decrypted` (+ envelope / PQC variants).
- Added live round-trips closing coverage gaps: `encrypt_headers=True`
  (BYOK, Python + Java — exercises the encrypted-AU-headers path) and a
  genomic-runs container (Python; the client path is content-agnostic
  and the daemon stores the `.tis` opaquely).
- Hardened the flaky `WorkbenchLiveTest.sessionCreateListTerminate`:
  wait out the `starting` state before terminating and tolerate a 409
  (already terminal/terminating) on DELETE.

### Added -- workbench client per-AU encrypted upload/download, envelope variant (2026-05-21)

Phase 4 (final) of the per-AU encrypted-upload rework: the **envelope**
variant, which wraps the per-run DEK under a symmetric AES-256-GCM
key-encryption key (KEK) — the GA, non-preview analogue of the PQC
ML-KEM path (Python + Java lockstep). This completes the BYOK + envelope
+ PQC trio, each with a live-daemon round-trip.

- Python `WorkbenchClient.upload_encrypted_envelope(*, project,
  container_uri, tio_path, kek, encrypt_headers=False)` /
  `download_decrypted_envelope(*, container_uri, kek, out_tio_path)`.
  Java `WorkbenchClient.uploadEncryptedEnvelope(...)` /
  `downloadDecryptedEnvelope(...)`. The fresh per-run DEK is wrapped with
  a 32-byte symmetric KEK (`aes-256-gcm`) and carried in the
  `ProtectionMetadata`; the daemon never holds the KEK. Recover with the
  same KEK; not preview-gated (unlike the PQC variant).
- Live round-trips added for all three variants in `workbench-live`
  (`test_per_au_encrypted_envelope_upload_round_trip` /
  `perAuEncryptedEnvelopeUploadRoundTrip`): encrypt → upload → download →
  decrypt → channel data matches; the wrong KEK fails to decrypt.

### Added -- workbench client per-AU encrypted upload/download, PQC variant (2026-05-21)

Phase 3 of the per-AU encrypted-upload rework: the post-quantum
(ML-KEM-1024) variant, gated behind `opt_pqc_preview` (Python + Java
lockstep, Decision 2). Unlike BYOK (Phases 1-2, caller-held key), the
per-run DEK is randomly generated, ML-KEM-wrapped into the
`ProtectionMetadata` packet, and recoverable only with the recipient's
ML-KEM private key — the daemon never holds a key.

- Python `WorkbenchClient.upload_encrypted_pqc(*, project,
  container_uri, tio_path, recipient_public_key, preview,
  encrypt_headers=False)` / `download_decrypted_pqc(*, container_uri,
  recipient_private_key, out_tio_path, preview)`. Java
  `WorkbenchClient.uploadEncryptedPqc(...)` /
  `downloadDecryptedPqc(...)`. All four refuse unless `preview` is set
  (`PQCPreviewDisabledError` / `PqcPreviewDisabledException`), mirroring
  the server's `opt_pqc_preview`.
- New transport helpers stamp/read the wrapped DEK on a run's
  `signal_channels` so `write_encrypted_dataset` carries it and the
  receiver can unwrap it: Python
  `transport.encrypted.stamp_transport_wrapped_dek` /
  `read_transport_wrapped_dek`; Java
  `EncryptedTransport.stampTransportWrappedDek` /
  `readTransportWrappedDek`. The wrapped DEK is stored as a `uint8`
  attribute array (not a VLEN string) so the v1.2 / ML-KEM blob's
  embedded NULs survive the round-trip.

Validated end-to-end against a live daemon (Python
`test_per_au_encrypted_pqc_upload_round_trip`, Java
`perAuEncryptedPqcUploadRoundTrip`: ML-KEM-wrapped encrypt → upload →
download → unwrap → decrypt → channel data matches; the un-previewed
call refuses and the wrong private key fails to decrypt).

### Added -- workbench client per-AU encrypted upload/download (Java) (2026-05-21)

Phase 2 of the per-AU encrypted-upload rework: the Java mirror of the
Python Phase 1 client (lockstep, Decision 2).

- `WorkbenchClient.uploadEncrypted(project, containerUri, tioPath, key,
  encryptHeaders)` — encrypts a *copy* of the plaintext `.tio` per-AU
  (`PerAUFile.encryptFile`) into a valid `.tis`
  (`EncryptedTransport.writeEncryptedDataset` via a
  `ByteArrayOutputStream` `TransportWriter`) and uploads it. The source
  is not mutated; the daemon stores/serves it opaque (server #31).
- `WorkbenchClient.downloadDecrypted(containerUri, key, outTioPath)` —
  downloads, `readEncryptedToPath` materialises the still-encrypted
  `.tio`, and returns `PerAUFile.decryptFile(...)` channels per run.
- The blob `uploadProtected` / `downloadAndOpen` (W6.2) now throw
  `UnsupportedOperationException` pointing to the per-AU methods.

Validated end-to-end against a live daemon
(`WorkbenchLiveTest.perAuEncryptedUploadRoundTrip`: encrypt → upload →
download → decrypt → channel bytes match the plaintext source).

### Added -- workbench client per-AU encrypted upload/download (Python) (2026-05-21)

Phase 1 of the per-AU encrypted-upload rework (replaces the
daemon-incompatible W6.2 blob BYOK).

- `WorkbenchClient.upload_encrypted(*, project, container_uri,
  tio_path, key, encrypt_headers=False)` — encrypts a *copy* of the
  plaintext `.tio` per-AU (AES-256-GCM, channel payloads + optional AU
  headers) into a valid `.tis` carrying a `ProtectionMetadata` packet,
  and uploads it. The daemon stores/serves it opaque (server #31); it
  never sees plaintext or holds a key.
- `WorkbenchClient.download_decrypted(*, container_uri, key,
  out_tio_path)` — downloads the encrypted container, materialises the
  still-encrypted `.tio`, and returns the decrypted channels
  (`{run: {channel: ndarray}}`).
- The blob-level `upload_protected` / `download_and_open` (W6.2) now
  raise `NotImplementedError` with a pointer to the per-AU methods:
  sealing a whole payload into one opaque ciphertext blob is rejected
  by the daemon (it validates uploads as transport streams).

Validated end-to-end against a live daemon
(`test_per_au_encrypted_upload_round_trip`: encrypt → upload → download
→ decrypt → channel values match). Java mirror + PQC variant follow.

### Fixed -- workbench upload ack-drain on websockets >= 14 (2026-05-20)

`UploadClient._drain_acks` read `self._ws.messages`, a deque the
`websockets` >= 14 asyncio `ClientConnection` no longer exposes, so
**every** live upload crashed with `AttributeError`. No live test had
exercised `upload_bytes` against a real daemon, so this went uncaught.
Rewritten to a cancel-safe non-blocking `recv()` drain (acks are still
fully processed by `_wait_for_done`).

Added a valid-`.tis` upload -> ingest -> download -> decode round-trip
to the `workbench-live` smoke (the first upload/download e2e there;
9/9 against a local daemon). This also documents -- via
`docs/parity-audit-v1.0.md` 3.2 -- that **W6.2 blob-level BYOK is not
daemon-compatible**: the daemon validates uploads as transport streams
(rejects opaque ciphertext) and re-encodes on download, so encrypted
upload must use per-AU encryption yielding a valid `.tis`. That rework
is the real W6.2 follow-up.

### Added -- W6.6: SDK reference docs + quickstart + finalisation (2026-05-20)

Closes **W6** (the final workbench-client milestone) — W1–W6 ship a
complete client SDK (Python + Java), the tio-browser GUI, the CLI,
and now reference docs.

- **Sphinx** (Python) already autoapi-covers the whole `ttio`
  package; added a `workbench` landing link + a **quickstart
  notebook** (`python/docs/tutorials/workbench_quickstart.ipynb`)
  rendered via `myst-nb` (execution off — the live flow is verified
  by the `workbench-live` smoke, not the docs build). The notebook
  walks the spec §8.3 path: connect → encode → upload → query →
  submit a pipeline → download, plus the federation no-op.
- **Javadoc** for the Java SDK (the `maven-javadoc-plugin` was
  configured but never run in CI).
- New **`docs` CI workflow** builds both (Sphinx + Javadoc) on PRs
  touching the SDK or docs, so doc-breaking changes are caught.

### Added -- W6.5: federation client (graceful no-op vs v1.0) (2026-05-20)

Client surface for the v1.1+ federation endpoint (spec §12.3) that
degrades gracefully against a v1.0 single-node server.

- New `ttio.workbench.federation.FederationClient` (Python) /
  `global.thalion.ttio.workbench.federation.FederationClient` (Java),
  reached via `client.federation()`.
- `peers()` hits `GET /v1/federation/peers`; on **404** (the v1.0
  server doesn't expose it) it returns an **empty list** instead of
  raising, so callers never special-case version detection. Other
  non-2xx statuses still raise. `is_federated()` / `isFederated()`
  is a convenience over an empty peer list.

Tested both languages incl. the 404→empty contract (Java via a
local `HttpServer` stub).

### Added -- W6.4b: ttio export --format expansion (2026-05-20)

Export-side mirror of W6.4a: `ttio export --format` graduates from
the `fastq | fasta` stub to the dataset-level exporters.

- New `ttio.exporters.registry` (mirrors `ttio.importers.registry`).
  Each format wraps its exporter in a uniform
  `adapter(tio_path, layer, output, **opts)`:
  - mzML / mzTab / ISA-Tab/JSON (whole-dataset writers),
  - BAM / CRAM (genomic-run writers; samtools at runtime; CRAM needs
    `--extra --reference <fasta>`).
- `fasta` / `fastq` keep their dedicated CLIs.
- Unknown formats fail with rc 3 and a clear message.

Parity: a test pins the GUI export writers against the CLI's. The
documented gaps are **nmrML / JCAMP-DX / imzML** -- those export
from per-spectrum / per-pixel objects and the Python side has no
`.tio`-layer→object extraction helper yet (GUI/Java only).

### Added -- ExportPanel progress bar (2026-05-20)

The tio-browser "Export container" dialog now shows a progress bar
during the open→export run instead of only a static status label,
mirroring the W6.1b EncodingPanel. The bar is indeterminate (the
open/export tasks don't report granular progress) but distinguishes
working from hung. Closes the last surface noted in #113 (a local
`.tio`→file export, not a client-server interaction).

### Added -- W6.3: PQC client (ML-KEM-1024 + ML-DSA-87, opt_pqc_preview) (2026-05-20)

Post-quantum payload protection on the workbench client, gated
behind the `opt_pqc_preview` flag (spec Decision 9), built on the
W6.2 envelope path + core `ttio.pqc` / `PostQuantumCrypto`.

- New `ttio.workbench.pqc` (Python) /
  `global.thalion.ttio.workbench.pqc` (Java):
  - `seal_pqc` / `open_pqc` -- envelope whose DEK is wrapped under an
    **ML-KEM-1024** encapsulation key; the recipient decapsulates
    with the matching private key.
  - Optional **ML-DSA-87** signing of the sealed ciphertext; the
    signer key + algorithm land in `ProtectionMetadata`, the
    detached signature rides alongside. `verify_pqc` checks it.
  - `kem_keygen` / `sig_keygen` passthroughs.
  - Every entry point refuses unless `preview=True`
    (`PQCPreviewDisabledError` / `PqcPreviewDisabledException`),
    mirroring the server's `opt_pqc_preview` gating.
- Fix: the W6.2 Java `WorkbenchEncryptor.openSealed` now dispatches
  unwrap on `kekAlgorithm` (the 2-arg `unwrapKey` only handled
  AES-256-GCM), so ML-KEM envelopes unwrap. Mirrors the Python
  `open_sealed`, which already passed the algorithm through.

Cross-language: the PQC-envelope `ProtectionMetadata` shape
(`kek_algorithm="ml-kem-1024"` / `signature_algorithm="ml-dsa-87"`)
is byte-identical across languages (matching anchor literals).
Round-trips are unit-level; live-daemon variant deferred. Tests:
`test_pqc.py` (6, crypto cases skip without liboqs),
`WorkbenchPqcTest` (6).

### Added -- W6.2: BYOK + envelope encryption client (Python + Java) (2026-05-20)

Client-side payload protection for encrypted workbench uploads
(spec UC-03.2/3), a thin wrapper over the existing core crypto
(`ttio.encryption` / `EncryptionManager` + the v1.2 DEK wrap).

- New `ttio.workbench.encryption` (Python) /
  `global.thalion.ttio.workbench.encryption` (Java):
  - `ProtectionMode` (BYOK / ENVELOPE).
  - `seal(payload, ...)` / `open_sealed(...)` -- **BYOK** seals under
    a researcher-supplied 32-byte DEK (key never leaves the client;
    empty `wrapped_dek`); **ENVELOPE** generates a fresh DEK, seals,
    and wraps the DEK under a KEK via the shared v1.2 wrap. Framing
    is `iv(12) || tag(16) || ciphertext` in both languages.
  - `ProtectionMetadata` with a canonical `to_json` / `from_json`
    (sorted keys, compact separators, standard base64) -- the
    cross-language anchor.
- `WorkbenchClient.upload_protected` / `download_and_open` (Python)
  and `uploadProtected` / `downloadAndOpen` (Java) seal-then-upload
  and download-then-decrypt, returning/consuming the
  `ProtectionMetadata`.

Cross-language: the BYOK `ProtectionMetadata` JSON is byte-identical
across Python and Java (anchored by matching test literals).
Round-trips are unit-level (BYOK + envelope); the live-daemon
variant is deferred. Tests: `test_encryption.py` (9),
`WorkbenchEncryptionTest` (8).

### Added -- W6.1b: tio-browser determinate transfer progress (2026-05-20)

Second slice of W6.1: wires the W6.1a `TransferProgress` callback
through the GUI so a running upload shows a real percentage
instead of an indeterminate spinner (the gap that made the
`whale_sequences` upload look hung).

- `WorkbenchClient` gains progress-bearing `upload` / `download`
  facade overloads delegating to the transport client.
- `TransferManager` passes a coalesced progress callback (at most
  one pending FX update, so a fast transfer can't flood the event
  loop) that drives each `Transfer`'s byte count + a
  `"Uploading... NN%"` / bytes-so-far message.
- `TransferQueueView` progress column is now a determinate
  fraction (`bytesTransferred / sizeBytes`) for uploads; downloads
  stream without a known total and stay indeterminate.
- `EncodingPanel` shows a progress bar bound to the encode task's
  message during the local encode phase (indeterminate -- the
  importer doesn't yet report granular progress, see #114), then
  hands off to the determinate Transfers queue for the upload.

### Added -- W6.1a: transport progress callback (Python + Java) (2026-05-20)

First slice of W6.1 (progress feedback). Adds a progress
callback to the workbench transport client so uploads/downloads
can drive a determinate progress bar instead of an indeterminate
spinner -- the root SDK gap behind issue #113.

Java (`global.thalion.ttio.workbench.transport`):
- New `TransferProgress` functional interface --
  `onProgress(bytesDone, bytesTotal)`, with an `UNKNOWN_TOTAL`
  (-1) sentinel for streamed downloads. Throwing callbacks are
  swallowed so they can't abort a transfer.
- `WorkbenchTransportClient.upload(..., TransferProgress)` and
  `.download(..., TransferProgress)` overloads. Upload reports
  `(bytesSent, payload.length)` per chunk (determinate); download
  reports `(bytesReceived, UNKNOWN_TOTAL)` per binary frame.

Python (`ttio.workbench.transport`):
- `UploadClient.upload_bytes(..., progress=)` and
  `DownloadClient.download(..., progress=)` accept a
  `Callable[[int, int], None]` with the same `(done, total)`
  contract; `_report_progress` helper swallows callback
  exceptions. Cross-language equivalent of the Java interface.

Tests: `TransferProgressTest` (Java, 3) pins the sentinel +
contract + that the overloads exist; `test_transport_progress.py`
(Python, 6) pins `_report_progress` + the `progress=` kwargs.
The end-to-end "callback fires with rising bytes" assertion lands
with W6.1b (GUI wiring) / a live-upload test.

### Changed -- W6 plan: progress feedback promoted to W6.1 (2026-05-20)

Reprioritised docs/workbench-client/W6-plan.md so **progress
feedback for all client-server interactions** leads the milestone
(was unscheduled; now W6.1, ahead of encryption/PQC). Encryption
-> W6.2, PQC -> W6.3, formats -> W6.4, federation -> W6.5, docs
-> W6.6. Driven by a live cross-environment test where a
multi-contig FASTA encode+upload showed no feedback for minutes
(working, not hung). Requirement: every client-server op (and the
encode/decode bracketing it) shows phase (reading/writing,
encoding/decoding, uploading/downloading) + % by source size;
no bare spinners. Tracked in issues #113 (requirement) and #114
(the ReferenceImport per-attribute H5Acreate slowness that
exposed it).

### Added -- W6.0: kickoff plan (SDK polish + formats + PQC/BYOK + federation) (2026-05-20)

Kickoff for the final workplan milestone. docs/workbench-client/
W6-plan.md lays out six sub-phases: W6.0 kickoff (this), W6.1
BYOK + envelope encryption client, W6.2 PQC client
(opt_pqc_preview), W6.3 format expansion (wire the spec §4 codecs
into `ttio encode --format` for CLI/SDK/GUI parity), W6.4
federation client (graceful no-op vs v1.0), W6.5 SDK reference
docs + tutorial + finalisation. Plan-only PR; sub-phases follow.

Grounding: the core crypto (ttio/pqc.py, ttio/encryption.py,
ttio/transport/encrypted.py) and the 13-format codec set already
exist at the library level; W6 exposes them through the workbench
client surface (no ttio/workbench/{encryption,pqc}.py yet) and
finalises docs. Python + Java lockstep per Decision 2; Rust SDK
deferred to v1.2.

### Changed -- Un-xfail cohort live test after server fix (2026-05-19)

tti-workbench-server PR #29 registered TTIOWBCohortsHandler (the
cohort REST plane was 404 in v1.0 -- the gap the live smoke
exposed). With the server fix on main, the live smoke's cohort
test (test_cohort_preview_count_round_trips) is un-xfailed and
now passes normally: the workbench-live workflow checks out the
server at main. Live smoke is now 8 passed (was 7 passed, 1
xfailed). Doc updated in docs/workbench-client/live-daemon-smoke.md.

### Added -- Live-daemon end-to-end smoke (closes the W1-W5 live-acceptance deferral) (2026-05-19)

Wires a real tti-workbench-server daemon into CI and drives the
actual ttio.workbench.* client SDK against it -- the live half
of the W1-W5 acceptance gates that every prior sub-phase
deferred.

Pieces:
- python/tests/integration/test_workbench_live.py -- env-gated
  (TTIO_WORKBENCH_URL + TTIO_WORKBENCH_STAGING) test driving
  connect (BootstrapAdminAuth) -> containers.list ->
  pipelines.register/list/get -> jobs.submit + poll + events
  (SSE) + cancel -> sessions.create/list/terminate. SKIPS when
  the env vars are unset, so the normal unit CI is untouched.
- scripts/workbench-live-smoke.sh -- local runner; boots the
  daemon with a temp SQLite config, seeds the admin project,
  runs pytest, tears down. Validated locally: 7 passed, 1
  xfailed.
- .github/workflows/workbench-live.yml -- CI that builds the
  GNUstep/ObjC toolchain + libTTIO (from the PR checkout) + the
  pinned tti-workbench-server, boots the daemon, runs the smoke.
  Runs on workbench-client-path PRs + manual dispatch (not a
  blanket per-PR gate; cold GNUstep build is ~10-15 min). Needs
  the TTIO_LIBRARY_CHECKOUT_TOKEN secret for cross-repo checkout.

Bugs the smoke caught:
- W5.2 containers.py keyword-arg bug (FIXED in this change).
  ContainersClient called http_json() and WorkbenchHttpError()
  with positional args for keyword-only params (scheme / token /
  body / status). The W5.2 unit tests only covered the pure
  dataclasses (HTTP methods are coverage-excluded), so this was
  invisible until the live round-trip. All five methods fixed.
- Server-side gap (tti-workbench-server repo follow-up):
  TTIOWBCohortsHandler is implemented but never registered in
  Source/Core/TTIOWBServer.m, so /v1/cohorts/{query,preview-count}
  return 404 on the v1.0 daemon. The cohort live test is marked
  xfail(strict=False) and flips to XPASS once the server wires
  the handler. The client SDK request is correct (raises a clean
  WorkbenchHttpError(404)).

Full detail in docs/workbench-client/live-daemon-smoke.md.

### Added -- W5.7: tio-browser Encoding + Export panels + 1.5.0 (closes W5) (2026-05-19)

Seventh and final W5 sub-phase. Adds the last two spec section
8.1 GUI components, bumps tio-browser to 1.5.0, and lands the
GUI-assembly end-to-end smoke. Completes the WC Desktop GUI
track: all nine spec section 8.1 components now have a working
tio-browser panel.

tio-browser new files in workbench/:
- EncodingPanel: modal encode + upload coordinator. Picks a
  source file, detects format via FormatSniffer, encodes to a
  temp .tio via the existing Phase-8 ImportTask, then enqueues
  an upload through the W5.3 TransferManager under a derived
  container URI. Static helpers deriveContainerUri (project +
  filename to uri:tio:...), deriveTempTio, isValidProject.
- ExportPanel: modal client-side export. Opens a local .tio
  (typically just downloaded via the W5.3 download dialog) and
  runs the existing Phase-9 ExportTask to a target format.
  Static helpers extensionFor (format to conventional
  extension), deriveExportTarget, isValidTioPath. v1.0 is
  client-side export; server-side export (export pipeline on the
  daemon) is a follow-up.

MainWindow gains two new menu items under Workbench: Encode +
upload and Export container.

Version: tio-browser 1.4.1 -> 1.5.0 (workbench-aware GUI). The
java/python ttio library stays at 1.3.0 (the W5.0 SDK bump);
tio-browser versions independently per the carry-forward rule.

Tests:
- EncodingPanelTest 10 tests: deriveContainerUri (simple name,
  path + extension strip, lowercase + hyphenate, no-project,
  empty-base fallback, dot-prefixed name, Windows backslash
  path); deriveTempTio (.tio suffix + name); isValidProject.
- ExportPanelTest 7 tests: extensionFor (known formats +
  unknown fallback); deriveExportTarget (extension swap,
  directory preserved, no-directory, null-source fallback);
  isValidTioPath.
- WorkbenchMenuSmokeTest 2 TestFX tests (end-to-end GUI
  assembly): the Workbench menu exposes every W5.1-W5.7 action;
  every Workbench action has an onAction handler. This is the
  GUI-assembly half of the W5 acceptance gate; the live-daemon
  round-trip (login to browse to upload to submit to download)
  remains a shared cross-W follow-up needing the workbench-server
  Docker image in CI.

W5 COMPLETE. Spec section 8.1 component coverage:
- Connection Manager (W5.1)
- Container Browser (W5.2)
- Upload/Download Manager + Selective Access Panel (W5.3)
- Cohort Query Builder (W5.4)
- Pipeline Launcher + Job Monitor (W5.5)
- Interactive Session Launcher (W5.6)
- Encoding Panel + Export Panel (W5.7)

Deferred across W5 (tracked in docs/workbench-client/W5.7-progress.md):
live-daemon round-trip smoke; schema-driven pipeline form;
save-as-cohort (needs server v1.1); embedded Jupyter WebView;
server-side export pipeline; tree-style cohort editor.

### Added -- W5.6: tio-browser Interactive Session Launcher (2026-05-19)

Sixth W5 sub-phase. JavaFX surfaces wrap the W4 SessionsClient
SDK -- no new wire surface, no new SDK code.

Decision (W5-plan open-question 2): interactive attach happens
through the operator's own WS-capable client (CLI
`ttio sessions attach` or a terminal), not an embedded JavaFX
WebView. The session list copies the `wss://` attach URL to the
clipboard and shows it in a copyable dialog. Embedded Jupyter
WebView (spec 7.4 step 3 option a) is a v1.1 enhancement.

tio-browser new files in workbench/:
- SessionLauncher: modal create form (project / engine pin /
  image / command / bind-mounts / env). Builds a
  SessionsClient.CreateRequest and calls create. Static parsers
  parseCommand, parseBindMounts, parseEnv are the testable
  boundary.
- SessionList: non-modal TableView Session with session id,
  status, project, engine, host-port columns; refresh,
  copy-attach-URL, terminate controls. Static attachUrl Session,
  WorkbenchClient returns the wss URL for running sessions, null
  otherwise.

MainWindow gains two new menu items under Workbench: Launch
session and Sessions after Jobs.

Tests:
- SessionLauncherTest 15 tests: parseCommand whitespace split +
  blank/null empty; parseBindMounts basic host:container, drops
  :mode suffix, skips blank lines, rejects missing/leading/
  trailing colon; parseEnv KEY=VALUE, value-with-equals
  preserved, skips blanks, rejects no-equals/leading-equals;
  isValidProject blank rejection.
- SessionListTest 4 tests: attachUrl builds the WS proxy URL for
  running sessions (wss + ws scheme variants); null for
  non-running / null args.

Cross-language scope: GUI-only Java; no new SDK code. The W4
Session record + SessionProxy URL builder already have
cross-language byte-equivalence pinned in W4 tests.

Deferred follow-ups: live-daemon round-trip smoke (shared);
embedded Jupyter WebView (spec 7.4 step 3 option a); auto-refresh
/ live session status.

### Added -- W5.5: tio-browser Pipeline Launcher + Job Monitor (2026-05-19)

Fifth W5 sub-phase. JavaFX surfaces wrap the W3 PipelinesClient
and JobsClient SDKs -- no new wire surface, no new SDK code.

tio-browser new files in workbench/:
- PipelineLauncher: modal pipeline picker. ChoiceBox loaded from
  PipelinesClient.list; pipeline metadata + schema preview pane;
  raw JSON textareas for inputs and params (v1.0 -- schema-driven
  form generation is v1.1); submit button calls
  JobsClient.submit and shows an info Alert with the resulting
  job id. Static isValidJsonObject(String) is the testable
  JSON-object validator.
- JobMonitor: non-modal job dashboard. TableView Job with job
  id, pipeline, status, project, queued, started, completed
  columns; status-filter ChoiceBox driving JobsClient.list;
  refresh, cancel-selected, tail-events controls. Static helpers
  formatTimestamp Long (ISO-8601 UTC) and filterValue String
  (resolves "(all)" to null) testable without FX.
- JobEventsView: non-modal SSE tail viewer for a single job.
  Worker thread calls JobsClient.events jobId, consumer; the
  consumer marshals each frame to the FX thread via
  Platform.runLater. Auto-closes the stream on terminal-state
  event (Completed / Failed / Cancelled). Static formatFrame
  JobEvent and isTerminalEvent JobEvent testable without FX.

MainWindow gains two new menu items under Workbench:
Launch pipeline and Jobs after Cohort query.

v1.0 scope: raw JSON textareas for pipeline inputs / params.
Schema-driven form rendering is a v1.1 enhancement; the SDK
already accepts Map String, Object so only the form pair
changes.

Tests:
- PipelineLauncherTest 9 tests: isValidJsonObject accepts
  populated objects, rejects blank / null / array / scalar /
  malformed; renderSchemaPreview surfaces identifier, version,
  project, owner, engine-pin, schemas; empty schemas render as
  none; engine-pin omitted when null/empty; stable output.
- JobMonitorTest 6 tests: formatTimestamp UTC ISO-8601 for
  positive seconds, blank for null / 0 / negative; filterValue
  resolves "(all)" / "" / null to null, passes through known
  statuses.
- JobEventsViewTest 7 tests: formatFrame event-name and
  key=value pairs, empty data map, null event name; isTerminalEvent
  detects completed / failed / cancelled, false for queued /
  starting / running, false for non-state events, missing
  status, null.

Cross-language scope: GUI-only Java; no new SDK code. The W3
Pipeline / Job / JobEvent records already have cross-language
byte-equivalence pinned in W3 tests.

Deferred follow-ups: live-daemon round-trip smoke (shared);
schema-driven inputs form; richer job filters (project, owner,
since); re-submit job.

### Added -- W5.4: tio-browser Cohort Query Builder (2026-05-19)

Fourth W5 sub-phase. Adds a visual predicate-composition window
to tio-browser that drives the W3 `CohortQuery` SDK; no new
wire surface (the v1.0 server's POST /v1/cohorts/query +
/v1/cohorts/preview-count endpoints already cover the
GUI-supported flows).

tio-browser new files in workbench/:
- CohortLeafRow -- mutable row for the leaf-predicate table
  with JavaFX-property-backed fields (kind / field / op /
  rawValue). toPredicate builds the matching CohortPredicate
  leaf; coerceValue raw, op converts the raw text input to a
  typed value (int / double / bool / string / comma-separated
  list for in).
- CohortQueryBuilder -- window with composite-root choice (AND
  / OR / NOT), select-kind choice (containers / subjects /
  samples), TableView CohortLeafRow with editable cells, Run
  and Preview Count buttons, result TableView. Static
  buildPredicate(composite, rows) is the testable boundary
  between form and SDK.

MainWindow integration: new Workbench menu item Cohort query
after Transfers.

Client-side rule enforcement (mirrors v1.0 server rules):
- phenotype rejected under OR / NOT (server enforces this
  too; client-side catch gives a clearer error than a 400).
- NOT requires exactly one leaf.

v1.0 scope: flat leaf list under a single composite root.
Nested composite trees + drag-drop reorder are a v1.1
enhancement. The W3 SDK already supports nested composites;
the GUI just does not surface them in v1.0.

Tests:
- tio-browser CohortLeafRowTest 13 tests: coerceValue parsers
  (int / double / bool / string / list for in / blank for
  exists), toPredicate builds the right leaf subclass per
  kind, blank field rejection, Kind.fromLabel round-trip.
- tio-browser CohortQueryBuilderTest 11 tests: AND with one
  leaf collapses to leaf, AND/OR composite construction,
  phenotype rejected under OR/NOT, NOT with non-1 leaves
  rejected, empty leaf list rejected, unknown composite
  rejected.

Cross-language scope: GUI-only; no new SDK code. The W3
cohort-predicate AST is already cross-language byte-equivalent.

Deferred follow-ups: live-daemon round-trip smoke (shared);
save-as-cohort needs server v1.1 POST /v1/cohorts; tree-style
editor with nested composites + drag-drop; field-level
autocomplete.

### Added -- W5.3: Transfer Manager + Selective Access Panel + filter builder (2026-05-19)

Third W5 sub-phase. Wires the W1 `WorkbenchTransportClient`
into a JavaFX queue UI, adds a visual filter form for selective
downloads, and pins a typed-builder filter API on both SDKs.

No new server wire surface; W5.3 builds on top of the existing
W1 `/transport` WS handshake. The download-filter allowlist
(`ms_level`, `polarity`, `retention_time_{min,max}`,
`precursor_mz_{min,max}`, `precursor_charge`, `max_au`) is
unchanged.

Java SDK new file `workbench/transport/SelectiveAccessFilter.java`:
- Fluent typed setters per allowed filter key.
- Per-key range checks throw `IllegalArgumentException`.
- Cross-key `validate()` checks `rt_max >= rt_min` and
  `mz_max >= mz_min`, throws `IllegalStateException`.
- `build()` returns a `LinkedHashMap<String, Object>` ready
  for `WorkbenchTransportClient.download`.

Python SDK new file `workbench/transport/selective_access.py`:
- Single-class mirror of the Java builder. Same method names
  (snake_case), same exception semantics (`ValueError` per-key,
  `RuntimeError` cross-key).

tio-browser new package `browser/workbench/`:
- `TransferKind` / `TransferState` enums + `Transfer` mutable
  entry with JavaFX properties.
- `TransferManager` -- process-wide singleton; daemon-thread
  executor; FX-thread-safe state mutation; `enqueueUpload`,
  `enqueueDownload`.
- `SelectiveAccessPanel` -- GridPane filter form with one
  input per allowed filter key; `buildFilter()` delegates to
  the SDK `SelectiveAccessFilter` and runs `validate()`.
- `UploadStartDialog` -- modal source file + project + URI
  picker; container-URI validator.
- `DownloadStartDialog` -- modal URI + destination + embedded
  selective-access panel.
- `TransferQueueView` -- non-modal TableView<Transfer> with
  indeterminate ProgressBar while RUNNING; bound to the
  manager's observable list.

`MainWindow` gains three new menu items under `Workbench`:
`Upload to workbench...`, `Download from workbench...`,
`Transfers...`.

Tests:
- Java `SelectiveAccessFilterTest` (18 tests) covers accept /
  reject / cross-key validation + cross-language anchor.
- Python `test_selective_access.py` (20 tests) mirrors the
  Java suite.
- tio-browser: `SelectiveAccessPanelTest`, `TransferStateTest`,
  `TransferTest`, `UploadDownloadDialogTest` (17 pure-unit
  tests total) cover parsers, state-machine predicates, queue
  defensive copy, and URI/project validators.

Cross-language byte-equivalence: canonical filter dict pinned
in both test suites is the fifth independent cross-language
anchor (after W1 handshake, W3 cohort predicate, W4 attach
handshake, W5.2 container list-page).

Coverage: no new excludes; `SelectiveAccessFilter` /
`selective_access.py` are pure data and stay measured.

Deferred follow-ups: live-daemon round-trip smoke (shared with
W1/W3/W4/W5.1/W5.2); true progress percentage (needs a
progress-callback API on `WorkbenchTransportClient`); pause /
cancel (needs a cancellation primitive on the W1 client);
chromosome / position genomic filters (not in v1.0 server
allowlist).

### Added -- W5.2: Container Browser + /v1/containers SDK (Python + Java) (2026-05-19)

Second W5 sub-phase. Adds the `/v1/containers` REST surface to
both SDKs (Decision-2 lockstep) and a JavaFX Container Browser
window to tio-browser.

Pre-implementation: surveyed `tti-workbench-server/Source/HTTP/
handlers/TTIOWBContainersHandler.{h,m}` for the wire contract.
Five endpoints: list / get / layers / manifest / delete; opaque
base64url cursor pagination; container shape
`{uri, project, owner, encrypted, storage_path, created_at,
updated_at}`; gates via `containers.read.any_project` /
`containers.delete.{any,own_uploads}` server-side; existence
never leaked (404 vs 403).

Java SDK (new package `global.thalion.ttio.workbench.containers`):
- `Container` / `ContainerDetail` records (list + detail shapes).
- `ContainerListPage` with `nextCursor` + `hasMore()`.
- `ContainerLayer` record.
- `ContainerManifest` outer record + nested `MsRunSummary` /
  `NmrRunSummary` / `GenomicRunSummary` records.
- `ContainersClient`: `list`, `get`, `layers`, `manifest`,
  `delete`.
- `WorkbenchClient.containers()` factory.

Python SDK (`ttio/workbench/containers.py`):
- Frozen-dataclass mirrors of every Java record.
- `ContainersClient` with identical method shape.
- `WorkbenchClient.containers()` factory.

tio-browser (`workbench/ContainerBrowser.java`):
- Modal-but-non-modal window. TableView<Container> with
  sortable URI / project / owner / encrypted / created /
  updated columns.
- Filter row (project / owner / limit) + Refresh button +
  Load-more button driving cursor pagination.
- Manifest pane in a SplitPane; selecting a row fetches and
  renders the manifest as plain-text summary.
- Opened via new `MainWindow` menu item `Workbench -> Browse
  containers...`; gated on `ConnectionManager.isConnected()`.

Tests:
- `ContainersTest` (Java, 13 tests) covers all records +
  parsing edge cases + the cross-language anchor.
- `test_containers.py` (Python, 14 tests) mirrors the Java
  suite.
- `ContainerBrowserTest` (tio-browser, 11 pure-unit tests)
  covers static helpers (`parseLimit`, `formatTimestamp`,
  `renderManifest`).
- `WorkbenchClientTest.w3W4W5SubClientsAreLive` /
  `test_client.test_w3_w4_w5_sub_clients_are_live` extended
  to assert `containers()` returns non-null.

Cross-language byte-equivalence: the GET /v1/containers
list-page parser is the fourth independent cross-language anchor
(after W1 handshake, W3 cohort predicate, W4 attach handshake).
Same literal JSON pinned in both test suites.

Coverage adjustments:
- `java/pom.xml` JaCoCo excludes extended to
  `ContainersClient*` (HTTP-method wrappers need a live daemon).
  Records stay measured.
- `python/pyproject.toml` `[tool.coverage.run].omit` extended to
  `*/workbench/containers.py` (matches jobs.py / pipeline.py).

Deferred follow-ups: live-daemon round-trip smoke (shared with
W1/W3/W4/W5.1); unified DataSourceTree merging local + remote
sources (per W5-plan); GUI surface for the layer breakdown
(client method present, no UI yet).

### Added -- W5.1: tio-browser Connection Manager (2026-05-19)

First substantive W5 sub-phase. Adds the JavaFX foundation that
every other W5 panel will hang off: an observable workbench
connection holder, a modal login dialog, and a status-bar
indicator. Wires the W1+W3+W4 `global.thalion.ttio.workbench.
WorkbenchClient` SDK into the GUI for the first time.

New package `tio-browser/src/main/java/global/thalion/ttio/
browser/workbench/`:
- `ConnectionState` enum (DISCONNECTED / CONNECTING / CONNECTED
  / FAILED).
- `ConnectionListener` functional interface.
- `ConnectionManager` -- process-wide singleton holding the
  `WorkbenchClient` instance + state machine. Thread-safe
  listeners; `connect(url, auth)` and `disconnect()` drive
  state transitions; throwing listeners do not block siblings.
- `LoginDialog` -- modal `Stage`-based form with server URL +
  username + password + TOTP. Static validators
  (`isValidUrl`, `isValidTotp`, `isValidUsername`,
  `isValidPassword`). `Bindings.createBooleanBinding`
  disables the Connect button until valid; worker `Task<Session>`
  calls `ConnectionManager.connect()` off the FX thread.
- `StatusIndicator` -- small HBox (coloured Circle + Label +
  Tooltip) for the status bar; subscribes to `ConnectionManager`
  and renders the four colour states.

`MainWindow` integration:
- New `Workbench` menu between `Transport` and `Tools` with
  `Connect...`, `Disconnect`, `Status...` items.
- Status bar gains the right-aligned indicator.
- `showWorkbenchStatus()` modal Alert summarises endpoint /
  user / provider / projects / capability count / session id;
  surfaces the last failure message when disconnected.
- `dispose()` detaches the indicator listener.

Tests:
- `ConnectionManagerTest` (8 tests, pure unit) covers state
  transitions, listener dispatch, idempotent disconnect,
  throwing-listener tolerance, post-failure reconnect.
- `LoginDialogTest` (9 tests) covers the static validators
  (URL forms, TOTP 6-digit rule, blank-rejection).
- `StatusIndicatorSmokeTest` (3 TestFX tests) verifies the
  three colour transitions render correctly and the tooltip
  carries the failure message.

CI: new `tio-browser-test` job in `ci.yml` mirrors the
release-shaded-jar workflow's linux-x64 leg (JDK 25 + HDF5 +
native libttio_rans build + `mvn install` of `java/` to local
M2 + `mvn -P linux-x64 test` of `tio-browser/`). Every PR is
now compile- and TestFX-verified at PR time, not just on tag
push.

Deferred to W5.1 follow-up: live-daemon round-trip (shared
deferral with W1/W3/W4).

### Changed -- W5.0: TTI-O Java SDK 1.3.0 (kickoff) (2026-05-19)

Fifth workbench-client milestone kickoff. Bumps `java/pom.xml`
1.2.0 -> 1.3.0 (line marker for the workbench-client Java
surface added in W1+W3+W4) and `tio-browser/pom.xml`'s
`<ttio.version>` in lockstep so the JavaFX GUI can consume the
new `global.thalion.ttio.workbench.*` classes. Corrects the
workplan's W5 cross-repo language: tio-browser is a sibling
subdirectory of `java/`, not a separate GitHub repository.

W5 phasing recorded in `docs/workbench-client/W5-plan.md` --
eight sub-phases (W5.0 kickoff -> W5.7 encoding+export+smoke)
delivering spec section 8.1's nine GUI components into
tio-browser. Subsequent W5.x PRs add the panels themselves;
this PR is admin-only (version bumps + plan + workplan
correction).

Library SemVer bump: minor. No breaking changes to existing
TTI-O format APIs; the workbench-client surface is purely
additive on top of v1.2.0.

### Added -- W4: interactive sessions client (Python + Java) (2026-05-19)

Fourth workbench-client milestone. Wraps the workbench server's
`/v1/sessions` REST surface + the `ttio-session-proxy` WS attach
helper (spec UC-11). Python + Java in lockstep per Decision-2.

Pre-implementation: deep-surveyed the v1.0.0 server contract --
`Documentation/session-protocol.md` cross-referenced with
`Source/HTTP/handlers/TTIOWBSessionsHandler.m`,
`Source/Sessions/{TTIOWBSessionRegistry,TTIOWBSessionLifecycle,
TTIOWBSessionProxy}.m`. Survey findings recorded in
`docs/workbench-client/W4-progress.md`. v1.0 deferrals
respected: idle-timeout sweep, host-port allocator, ring-buffer
backpressure are all server-side -- client just observes.

Python (`python/src/ttio/workbench/`):
- `sessions.py` -- `Session` dataclass (5 status enum, runtime
  fields, `is_terminal` / `is_attachable` properties),
  `SessionsClient` (create / list / get / terminate),
  `validate_bind_mounts()` client-side validator mirroring the
  server's rules (absolute, no `..`, project-scope check when
  `container_storage_root` is known).
- `session_proxy.py` -- `build_attach_handshake()` and
  `session_proxy_url()` pure helpers, `SessionProxyAttach` async
  context manager that opens the `ttio-session-proxy` WS, sends
  the JSON attach frame, pumps bytes bidirectionally between
  caller-supplied byte streams (stdin/stdout or in-memory).
- `client.py` -- `WorkbenchClient.sessions()` /
  `.session_create()` / `.session_proxy()` promoted from W2-era
  stub to live methods. No remaining `NotImplementedError`
  paths -- the v1.0 client SDK is feature-complete for the
  spec section 8.3 sample.

Java (`java/src/main/java/global/thalion/ttio/workbench/`):
- `sessions/Session.java` -- mirror record.
- `sessions/BindMountValidator.java` -- same validation rules.
- `sessions/SessionsClient.java` -- REST surface with
  fluent `CreateRequest` builder.
- `sessions/SessionProxy.java` -- pure builders + URL constructor.
- `sessions/SessionProxyAttach.java` -- callback-driven WS
  attach (built on `org.java_websocket`), pumps `InputStream`
  <-> `OutputStream` until close.
- `WorkbenchClient.sessions()` / `.sessionProxy()` live.

CLI (`python/src/ttio/tools/workbench_cli.py`): `ttio sessions`
promoted from W3-era stub to live verb subcommand:
  - `ttio sessions create --engine X --project Y [--image Z]
                         [--command CMD...] [--env K=V]
                         [--bind-mount HOST:CONT]`
  - `ttio sessions ls [--status X] [--limit N]`
  - `ttio sessions status <id>`
  - `ttio sessions attach <id> [--path /]` -- proxies stdin /
    stdout against the engine subprocess.
  - `ttio sessions terminate <id>` -- DELETE /v1/sessions/{id}.

Tests:
- **Python:** 156 workbench tests pass locally (+25 over W3):
  `test_sessions.py` (25 tests) covering Session parsing across
  all statuses, bind-mount validator (5 failure paths + happy +
  noop), attach handshake builder, URL constructor, status-set
  pinning, cross-language anchor literal.
- **Java:** `SessionsTest` mirrors the Python suite including
  the cross-language attach-handshake JSON literal.

Cross-language anchor: identical attach-handshake JSON literal
pinned in both Python + Java session test suites.

JaCoCo + Python coverage excludes extended to the new
daemon-required classes (`SessionsClient`, `SessionProxyAttach`
on Java; `session_proxy.py` on Python). The pure-data `Session`
record + `BindMountValidator` + attach builders stay measured.

### Added -- W3: cohort + pipeline + job client (Python + Java) (2026-05-19)

Third workbench-client milestone. Wraps the workbench server's
REST + SSE surface for the cohort / pipeline / job plane (spec
UCs 6, 7, 9, 10). Python + Java in lockstep per the workplan
Decision-2 parity rule.

Pre-implementation: deep-surveyed the v1.0.0 server wire
contract for the W3 endpoints (`/v1/cohorts/{query,preview-count}`,
`/v1/pipelines{,/{id}}`, `/v1/jobs{,/{id},/{id}/events}`).
Survey findings recorded in
`docs/workbench-client/W3-progress.md`. v1.0 deferrals
documented:
  - No `GET /v1/cohorts` / `POST /v1/cohorts` (saved cohorts);
    queries are ephemeral.
  - No `has_layer(...)` assay-availability filters in the cohort
    AST.
  - No `Last-Event-Id` SSE resumption; reconnect for full replay.
  - No `GET /v1/containers/{uri}/provenance`; `provenance_edges`
    lives in the DB but isn't surfaced over HTTP.

Python (`python/src/ttio/workbench/`):
- `cohort.py` -- predicate AST (4 leaf kinds + 3 composites).
  `OR` / `NOT` reject nested phenotype leaves per the server's
  column-join semantics. Operator overloading (`&`, `|`, `~`)
  for fluent composition. `CohortQuery` builder, `CohortResult`
  parser.
- `pipeline.py` -- `Pipeline` dataclass + `PipelinesClient`.
- `jobs.py` -- `Job` dataclass with `is_terminal`; `JobsClient`
  (submit / list / get / cancel); `JobEvent` + async-iterator
  `events(job_id)` SSE parser. `build_cohort_input()` builds
  the Decision-4 envelope.
- `_http.py` -- internal REST helper (urllib-based; zero new
  runtime deps).
- `client.py` -- `WorkbenchClient.query()` /
  `.preview_count()` / `.pipelines()` / `.submit_pipeline()` /
  `.jobs()` promoted from W2-era stubs to live implementations.
  Only `session_create()` remains a W4 stub.

Java (`java/src/main/java/global/thalion/ttio/workbench/`):
- `cohort/` -- `CohortPredicate` abstract + 7 subclasses with the
  same validation rules as Python. `CohortQuery` builder,
  `CohortResult` record.
- `pipeline/` -- `Pipeline` record + `PipelinesClient`.
- `jobs/` -- `Job` + `JobEvent` records + `JobsClient` with
  callback-driven SSE (`events(jobId, Consumer<JobEvent>)`).
- `WorkbenchHttp` -- internal REST helper over
  `java.net.http.HttpClient`; emits compact JSON byte-matching
  Python.
- `WorkbenchClient.query()` / `.previewCount()` / `.pipelines()` /
  `.jobs()` live.

CLI: W2-era placeholder subcommands promoted to live impls --
`ttio query` / `submit` / `jobs {ls,status,cancel,events}` /
`pipelines {ls,get,register}` / `cohorts` (alias for `query`).
`ttio provenance` surfaces the v1.0 "endpoint not exposed"
deferral. `ttio sessions` still W4-stubbed.

Tests:
- **Python:** 131 workbench tests pass locally (+52 over W2):
  `test_cohort.py` (33), `test_jobs.py` (13), `test_pipelines.py`
  (4), plus updates to W2 tests for W3-promoted methods.
- **Java:** `CohortPredicateTest` (predicate AST + CohortQuery),
  `JobAndPipelineTest` (Job / Pipeline / JobEvent / clients).
- Cross-language anchor: identical predicate JSON literal pinned
  in both suites.

JaCoCo + Python coverage excludes extended to the new
daemon-required classes (`WorkbenchHttp`, `PipelinesClient`,
`JobsClient` in Java; `_http.py`, `pipeline.py`, `jobs.py` in
Python). Pure-data records + predicate AST + parsers stay
measured.

### Added -- W2: `ttio` CLI umbrella + Python SDK foundation (2026-05-19)

Second milestone of the
[Workbench Client workplan](docs/workbench-client-workplan.md).
Lands the spec section 8.2 CLI surface verbatim and the section
8.3 SDK shape (`ttio.connect(...)`, auth providers,
`WorkbenchClient`). All Python; the Java equivalent for W2 lives
inside W5 when tio-browser bumps its TTI-O dep.

SDK foundation (`python/src/ttio/workbench/`):
- `auth_providers.py` -- `AuthProvider` ABC plus four concrete
  providers: `PasswordTotpAuth` (interactive creds via W1
  `login_password`), `BearerAuth` (caller already holds a
  bearer; synthesises a `Session` without round-trip),
  `BootstrapAdminAuth` (reads `<staging_root>/bootstrap-
  credentials.json`; smoke / dev path), and `OIDCAuth` (v1.1
  stub that raises a clear NotImplementedError).
- `client.py` -- `connect(url, auth=...)` factory + the
  `WorkbenchClient` class. Resolves WSS/WS/HTTPS/HTTP URLs,
  authenticates through the provider, exposes
  `upload_client(...)` / `download_client(...)` builders + the
  `upload_bytes(...)` / `download_bytes(...)` async
  convenience methods. W3 surfaces (`query`, `submit_pipeline`,
  `jobs`) and W4 surfaces (`session_create`) registered as
  methods that raise NotImplementedError pointing to the
  milestone.
- `cohort.py`, `pipeline.py`, `jobs.py`, `sessions.py` --
  namespace stubs so the eventual W3/W4 implementations have
  the import path reserved and IDE-completion works today.
- Top-level `ttio` re-exports `connect`, `WorkbenchClient`,
  `Session`, and all four `*Auth` providers so the spec section
  8.3 sample (`ttio.connect(..., auth=ttio.OIDCAuth())`) works
  without operators digging into sub-modules.
- `parse_filter_kv()` helper turns repeated `--filter k=v`
  arguments into a dict with numeric coercion.

CLI (`python/src/ttio/tools/workbench_cli.py`):
- `ttio` umbrella console-script with subcommands matching
  spec section 8.2:
    - `login` -- resolve credentials, print the auth JSON.
    - `upload` -- WS upload of a local `.tio`.
    - `download` -- WS download to a local `.tio` with optional
      selective-access filters (`--filter k=v`, repeatable).
    - `stream` -- WS download saved as raw `.tis`.
    - `inspect` -- stats-only WS read; prints the per-AU
      summary frames.
    - `encode` / `export` -- dispatch into the existing
      `ttio.tools.{fastq,fasta}_{import,export}_cli`.
    - `query`, `submit`, `jobs`, `cohorts` -- W3 placeholders
      (exit 2 with milestone deferral message).
    - `sessions` -- W4 placeholder.
- Auth-mode resolver enforces exactly one of `--token+--owner`,
  `--staging-root`, or `--username+--password+--totp`.
- `pyproject.toml`: `[project.scripts] ttio = ...` so
  `pip install -e .` exposes `ttio` on PATH.

Tests (`python/tests/workbench/`):
- `test_client.py` -- 21 tests covering top-level re-exports,
  URL parsing, all four auth providers, the `connect()`
  factory, W3/W4 placeholder method dispatch, and the
  filter-key parser.
- `test_cli.py` -- 21 tests covering every spec section 8.2
  verb's `--help`, the auth-mode resolver's three failure
  modes, the W3/W4 placeholder exit codes, the encode-format
  pointer to W6.

Java (`java/src/main/java/global/thalion/ttio/workbench/`) --
**added in the same PR after a cross-language parity review:**
- `auth/AuthProvider.java` interface + four concrete providers
  (`PasswordTotpAuth`, `BearerAuth`, `BootstrapAdminAuth`,
  `OIDCAuth` stub).
- `WorkbenchClient.java` top-level entry: `connect(url, auth)`
  factory, `transportClient()` builder, `upload()` / `download()`
  convenience methods, W3 / W4 placeholder methods, and
  `parseUrl()` byte-matching the Python URL parser.
- JaCoCo excludes extended to the new daemon-required classes
  (`BootstrapAdminAuth`, `PasswordTotpAuth`).
- 23 new Java tests covering URL parser, all four auth
  providers, `connect()`, `reauth()`, `close()`, W3 / W4
  placeholders.

Workplan amendment in this PR
(`docs/workbench-client-workplan.md` Decision 2): **Python +
Java SDK ship in lockstep at every milestone.** The `ttio` CLI
stays Python-only by design (Decision 1: console-script). ObjC
stays server-runtime (not extended for client purposes).

Coverage: 79 Python workbench tests (37 W1 + 42 W2) plus 83
Java workbench tests (60 W1 + 23 W2), all pass locally. The W1
cross-language byte-equivalence anchor (handshake JSON literals
in both test suites) carries forward into W2.

### Added -- W1: Workbench client (Python + Java) (2026-05-19)

First milestone of the
[Workbench Client workplan](docs/workbench-client-workplan.md).
Ships the workbench-aware transport client (Python + Java) that
speaks `tti-workbench-server` v1.0.0's auth-bearing handshake.
Replaces the existing reference-protocol-only clients
(`ttio.transport.client` / `global.thalion.ttio.transport.TransportClient`)
which target the Python reference server.

Python (`python/src/ttio/workbench/`):
- `auth` -- RFC 6238 TOTP (HMAC-SHA1, 30s, 6 digits), `Session`
  dataclass, `login_password(host, port, user, pass, totp)` POSTing
  to `/v1/auth/login`. Typed exceptions for 401 / 423 / 429 paths;
  `Retry-After` surfaced on rate-limit.
- `transport.handshake` -- pure JSON builders for upload + download
  first-frames. Client-side filter-key validation matching the
  daemon's accept list.
- `transport.upload.UploadClient` -- async context manager that
  drives one upload over `ws://host:port/transport` with the
  `ttio-transport` subprotocol. Per-AU acks, EndOfStream
  acknowledgment, resumable-upload support via `ResumeState`.
- `transport.download.DownloadClient` -- analogous async download
  with selective-access filtering, three output modes (binary /
  stats-only / stats-with-payload), stats-frame collection.

Java (`java/src/main/java/global/thalion/ttio/workbench/`):
- `WorkbenchJson` -- minimal compact JSON encoder + parser scoped
  to the handshake / ack frames. No Jackson / Gson dependency
  (matches the existing hand-rolled pattern in `BamDump` +
  `ProvenanceJsonParse`).
- `auth.Totp`, `auth.Session`, `auth.Login` -- Java mirror of the
  Python auth module, using `java.net.http.HttpClient` for the
  REST POST.
- `transport.WorkbenchHandshake` -- pure JSON builders + parser.
- `transport.WorkbenchTransportClient` -- end-to-end upload +
  download built on `org.java_websocket.WebSocketClient` (the
  existing TTI-O Java WS dep), with builder construction
  (`WorkbenchTransportClient.forSession(host, port, session)`).
- `transport.ResumeState` -- resume bookkeeping record.
- `transport.WorkbenchTransportException` -- base + Handshake +
  Upload + Download subclasses carrying WS close code + reason.

Tests:
- Python (`python/tests/workbench/`): 37 tests across `test_auth.py`,
  `test_handshake.py`, `test_cross_language.py`. All pass.
- Java (`java/src/test/java/global/thalion/ttio/workbench/`):
  `TotpTest`, `WorkbenchHandshakeTest`. Both suites pin against
  the same RFC 6238 TOTP vectors + the same handshake JSON
  literals as the Python tests -- this is the cross-language
  byte-equivalence anchor (a Python or Java drift will fail both
  sides).
- Daemon round-trip integration tests deferred to a W1 follow-up:
  they need a running `tti-workbench-server` binary in CI, which
  requires building the binary on the runner or vendoring a
  prebuilt -- out of scope for this PR.

Per the [W1 progress doc](docs/workbench-client/W1-progress.md),
W2 (`ttio` CLI umbrella + Python SDK foundation) follows.

## [1.4.1] - 2026-05-11

This release was re-tagged in flight: the initial v1.4.1 build shipped
class files compiled with `--enable-preview` (Java 21 preview FFM API),
which a stock JDK launcher refused to load with
"Preview features are not enabled for `Hdf5CompoundIO$FieldKind`
(class file version 65.65535)". The retag drops `--enable-preview`
and targets Java 22 (where FFM is a stable API). Tag and release
assets at https://github.com/DTW-Thalion/TTI-O/releases/tag/v1.4.1
were replaced.

### Fixed

- **`tio-browser` Windows: opens `.tio` files with genomic data.** v1.4.0
  launched cleanly on Windows (after the MinGW DLL closure was bundled),
  but it could not open `.tio` files containing BAM/CRAM-style genomic
  data because the Java implementation shelled out to the
  `samtools` CLI binary — not installed on a fresh Windows machine.
  Replaced the samtools subprocess across `BamReader`, `BamWriter`,
  `CramReader`, and `CramWriter` with [htsjdk](https://github.com/samtools/htsjdk)
  4.1.3 — the pure-Java SAM/BAM/CRAM library used by GATK, Picard, and
  IGV. No external binary required at runtime; tio-browser now works
  end-to-end on Windows with only a JDK 22+ installed.
- **`tio-browser`: opens `.tio` files on a fresh JDK install.** The
  first v1.4.1 build hit a class-loader rejection on
  `Hdf5CompoundIO$FieldKind` because the compile used `--enable-preview`
  (Java 21 preview FFM API). Bumped the compile target from Java
  21+preview to Java 22-stable; FFM is a stable API in JDK 22 and the
  API surface is identical, so no source changes were needed.

### Changed

- **Minimum JDK bumped from 17 to 22.** The FFM API used by
  `Hdf5CompoundIO` / `VlBytesFFM` is stable in JDK 22+; an older JDK
  cannot load the compiled classes (class file version 66).
  Adoptium/Temurin ships JDK 22 binaries for all 3 supported
  platforms.

### Internal

- `htsjdk` added as a Maven runtime dependency (`com.github.samtools:htsjdk:4.1.3`).
- `BamReader.SamtoolsNotFoundException` retained as a no-throw alias
  for source compat with callers and tests that catch it.
- `BamReader.isSamtoolsAvailable()` returns `true` unconditionally
  (htsjdk is always available as a Maven dep).
- CRAM reference handling: custom `InMemoryFastaReferenceSource`
  bypasses htsjdk's stock `ReferenceSource` strict length/MD5
  validation, matching samtools' lenient default behaviour. Existing
  samtools-produced CRAMs (with placeholder `@SQ LN`) decode without
  modification.
- `BamReader` static initializer sets
  `samjdk.cram.use_alignment_md5_check=false` for the same reason.
- `WrittenGenomicRun.qualities` byte semantics preserved: ASCII Phred+33
  on the cross-language wire (M87/M88 convention); htsjdk's raw-Phred
  byte arrays are converted with ±33 in the reader/writer.
- 0-byte BAM/SAM rejection: explicit `Files.size(path) == 0` check
  before opening (htsjdk would otherwise treat as a 0-record BAM).
- `enable-preview` compiler flag removed from `maven-compiler-plugin`
  configuration in `java/pom.xml`.
- `enable-preview` removed from `surefire-plugin` `<argLine>` in both
  `java/pom.xml` and `tio-browser/pom.xml`.
- CI workflows (`ci.yml`, `release-shaded-jar.yml`): `setup-java`
  `java-version: '21'` → `'22'` (4 invocations).
- `tio-browser/README.md` + `docs/tio-browser.md`: "JDK 17+" → "JDK 22+".

## [1.4.0] - 2026-05-09

### Changed

- **`tio-browser` distribution model: per-platform JARs.** Replaces the
  prior universal shaded JAR with three platform-specific JARs
  (`tio-browser-1.4.0-linux-x64.jar`, `tio-browser-1.4.0-mac-aarch64.jar`,
  `tio-browser-1.4.0-win-x64.jar`). Each JAR is ~31 MB instead of the
  universal-with-HDF5 alternative (~64 MB), carries only its own
  platform's natives, and **bundles HDF5 1.14 + the LZ4 filter plugin**
  so `java -jar tio-browser-1.4.0-<your-os>.jar` on a fresh machine
  works out of the box with only a JDK 17+ — no system HDF5 install
  required.
- New `Hdf5NativeLoader` extracts the bundled HDF5 native libs to a
  per-JVM temp dir at `App.start()`, calls `System.load` in dependency
  order (hdf5 → hdf5_hl → hdf5_java), and registers the LZ4 plugin
  search path via `H5.H5PLappend`. Idempotent; throws
  `Hdf5NativeLoadException` on hard failures (modal Alert + exit, with
  a headless detector that suppresses the exit during TestFX runs).
- Wrong-JAR-for-OS detection: running `tio-browser-1.4.0-linux-x64.jar`
  on a Mac shows a clear modal Alert with the correct download name.
- `release-shaded-jar.yml` workflow restructured: each platform's
  build job now produces its own complete shaded JAR end-to-end (no
  separate assembly job). Workflow grants `contents: write`
  permission so the auto-publish step works without the v1.3.0 manual
  workaround.
- `jarhdf5` switched from `<scope>system</scope>` to a vendored
  `org.hdfgroup:jarhdf5:1.14.6` Maven dep at `tio-browser/local-repo/`,
  so the JHI5 classes ship in each per-platform shaded JAR.

## [1.3.0] - 2026-05-09

### Added

- **`tio-browser` desktop GUI (Phases 0–13 + native bundling)** —
  JavaFX desktop application for inspecting `.tio` datasets, peer to
  the Java / Python / ObjC reference implementations. Cross-platform
  shaded jar bundles `libttio_rans_jni` for **Linux x86_64**, **macOS
  Apple Silicon (arm64)**, and **Windows x86_64**; end users can run
  `java -jar tio-browser-<ver>-shaded.jar` without any toolchain
  setup beyond a JDK 17+ runtime.

  - **Phase 8** — Import wizard: 13-format dispatch (mzML, ImzML,
    nmrML, JCAMP-DX, BAM/SAM/CRAM, FASTA, FASTQ, mzTab, Thermo,
    Waters, Bruker), drag-and-drop with format auto-detection.
  - **Phase 9** — Export dialog: 11 export formats with eligibility-
    based greying when the open `.tio` doesn't contain the required
    run kind.
  - **Phase 10/11** — Transport: download `.tis` streams from
    `http(s)`/`ws(s)` URLs into a local `.tio`; upload a local `.tio`
    as a `.tis` byte stream to the same URL families. Server-side
    filters (run kind, dataset-id list, RT range) let clients fetch
    subsets without downloading the whole file.
  - **Phase 12** — Diagnostics dialog (Tools → Diagnostics): probes
    HDF5 JNI, `samtools`, `ThermoRawFileParser`, `masslynxraw`, and
    the Bruker Python helper. Greys out Import/Export format rows
    whose binary isn't available, with tooltips listing the missing
    dep. Re-probe button picks up newly-installed binaries without
    restarting the app via a listener bus.
  - **Phase 13** — Distribution: `--open <path>` CLI flag for the
    shaded JAR opens a dataset at launch; `mvn -P native-package
    package` runs `jpackage` to produce platform-native installers
    (`.deb`/`.rpm`/`.dmg`/`.msi`) at `target/installer/`. `NativeLibraryLoader` resolves the
  library via `System.loadLibrary` → bundled-resource extract →
  graceful degradation (Intel Mac and other unbundled platforms get
  a placeholder in the genomic Read Inspector; all non-genomic
  features keep working). Built and validated via the
  `release-shaded-jar.yml` GitHub Actions matrix workflow
  (`ubuntu-22.04` + `macos-14` arm64 + `windows-2022` MinGW UCRT64).
  See [`tio-browser/README.md`](tio-browser/README.md) for the
  install matrix, build-from-source instructions, release path, and
  rationale for arm64-only macOS / MinGW-w64 Windows.

- `Quantification.unit` field (Java / Python / ObjC) — optional
  per-quantification unit string (e.g. `"fmol"`, `"ng/mL"`). Stored
  as a JSON-array sidecar attribute `@quantification_units` on
  `/study/quantifications` for backward-compat (legacy datasets read
  back with empty units). Surfaced as a column in `tio-browser`'s
  Quantifications tab.

- `MassSpectrum.isCentroided()` / `is_centroided` / `isCentroided`
  (Java / Python / ObjC) — per-spectrum centroid-vs-profile
  classification, stored as a parallel-array column
  `spectrum_index/centroideds` (int32, 0=profile, 1=centroided).
  Wire-format additive optional column; legacy files read as
  `false` for all rows. Used by `tio-browser`'s spectrum plot to
  auto-select stem rendering for centroided MS.

- `AcquisitionRun.spectra() : List<Spectrum>` (Java) /
  `acquisition_run.spectra` property (Python) /
  `-[TTIOAcquisitionRun spectra]` (ObjC) — modality-uniform spectrum
  enumeration, replacing the `for(i)+objectAtIndex+instanceof`
  pattern. Includes new `AcquisitionMode` constants
  `RAMAN` / `IR` / `UV_VIS` (ordinals 9 / 10 / 11) and an
  `AcquisitionRun.solvent` attribute (default `""`) for NMR runs.

- `ReferenceImport.writeToDataset` (Java),
  `-[TTIOReferenceImport writeToDataset:overwrite:error:]` (ObjC) —
  public counterpart to Python's `ReferenceImport.write_to_dataset`,
  closing the cross-language API parity gap on the reference-write
  path. Each language now has matching `readFromGroup` /
  `writeToDataset` symmetry. All three writers produce a
  byte-identical `/study/references/<uri>/` subtree (same `@md5`,
  `@reference_uri`, per-chromosome `@length`, and ZLIB-compressed
  UINT8 `data`); verified by structural comparison across Python,
  Java, and ObjC fixtures.

---

## [1.2.0] — 2026-05-08

### Added
- **`MSImage.mz_axis`**: shared m/z spectral axis on `MSImage` across
  Java, Python, and ObjC. Persisted as a 1-D FLOAT64 dataset under
  `/study/image_cube/mz_axis`. Required for imzML export of an
  `MSImage`-bearing `.tio`.
- **`MSImage.toPixelSpectra()` / `to_pixel_spectra()` / `-pixelSpectra`**:
  continuous-mode projection of the cube into per-pixel `(mz, intensity)`
  records suitable for `ImzMLWriter.write`.
- **`SpectralDataset.image()` / `.image` / `-msImage`**: accessor on the
  open dataset returning the materialised `MSImage` if `/study/image_cube`
  is present. Pattern mirrors the 1.1.0 `references()` accessor.
- **Python `SpectralDataset.write_minimal(image=...)`**: high-level
  kwarg writes the image cube alongside runs.
- **Python `MSImage.write_to(study_group)` / `MSImage.read_from(study_group)`**:
  standalone storage methods mirroring Java's `writeTo` / `readFrom`.
- **ObjC `TTIOPixelSpectrum`**: new value class for `-pixelSpectra` output.

### Backwards compatibility
- v1.1.x `.tio` files without `mz_axis` read as empty axis; the imzML
  exporter raises a clear error pointing at re-import.
- v1.1.x readers transparently skip the `mz_axis` dataset when reading
  v1.2.0-written files.

### Cross-language conformance
- New `python/tests/conformance/test_msimage_xlang.py` asserts byte-equal
  `mz_axis` payloads when written by Python and read back by Java + ObjC.
- Bug fix: `Hdf5Group.openDataset` (Java) now correctly handles N-D datasets,
  fixing a silent length-truncation that previously made Python-written
  `.tio` files unreadable from Java.

## [1.1.0] — 2026-05-06

Pure additive release. No wire-format change: `.tio` files written
by 1.0.0 are read identically by 1.1.0 and vice versa.

### Added

- `SpectralDataset.references()` (Java),
  `SpectralDataset.references` property (Python), and
  `[TTIOSpectralDataset references]` (ObjC) — enumerates embedded
  references at `/study/references/<reference_uri>/` for opened
  datasets, keyed by reference URI. Datasets written without
  embedded references (writer flag `embedReference = false`) return
  an empty map regardless of whether individual genomic runs carry a
  `referenceUri`. Cross-language parity verified by
  `python/tests/conformance/test_references_xlang.py` (9 directed
  pairs).
- `ReferenceImport.readFromGroup` (Java) /
  `ReferenceImport.read_from_group` (Python) /
  `+[TTIOReferenceImport readFromGroup:]` (ObjC) — factory that
  materialises a `ReferenceImport` from its on-disk group.

### Fixed

- ObjC writer's reference-embed path no longer requires `libttio_rans`
  to be available. Embedding `/study/references/<uri>/...` is pure HDF5
  I/O and now fires whenever `embedReference=YES` on a
  `TTIOWrittenGenomicRun`, matching Python's behavior. Signal-channel
  encoding via REF_DIFF_V2 still gates on the native lib (unchanged).
  Resolves the writer-gate asymmetry finding from Phase 0 Task 0.6.

### Changed

- Unified `@md5` attribute computation on
  `/study/references/<uri>/` to a single seq-only form across all
  writers and helpers. Previously, REF_DIFF_V2 auto-embed used
  `MD5(seq_a || seq_b || ...)` (sorted by name) while FASTA-import
  writers and the public canonical helpers
  (`compute_reference_md5` / `ReferenceImport.computeMd5` /
  `+[TTIOReferenceImport computeMd5WithChromosomes:sequences:]`)
  used `MD5(name_a || 0x0A || seq_a || 0x0A || ...)`. All three
  paths now agree on the seq-only form, which was already the
  authoritative on-disk digest. Existing v1.0.0 files written via
  REF_DIFF_V2 auto-embed are unchanged (their on-disk `@md5` was
  already seq-only). Existing v1.0.0 files written via FASTA-import
  retain their on-disk `@md5` verbatim through the v1.1.0
  `read_from_group` / `readFromGroup` path; only the auto-recompute
  fallback (when `md5=None` / `md5=null` / `md5:nil` is passed to
  the constructor) now produces a seq-only digest. Resolves the
  three-form `@md5` finding from Phase 0 Task 0.6.

### Notes

- ObjC has no canonical library-version constant; the version bump
  applies only to the Java pom and Python `__version__` /
  `pyproject.toml` metadata.
- The on-disk reference layout itself is unchanged from 1.0.0 —
  v1.1.0 only fills in the previously-missing read path. See
  `docs/format-spec.md` §10.10 (subsection "Reading embedded
  references") for the exact byte layout and the `@md5` form note.

---

## [post-v1.0.0 perf + parity tweaks — included in 1.1.0]

All correctness-neutral (same wire bytes, same on-disk container).
Headline numbers + reproducer instructions consolidated in
`docs/benchmarks/2026-05-05-v1.0-comprehensive-perf-report.md`
§11.

### Performance

- **24× FASTQ re-export speedup** (Python). The hot loop in
  `FastqWriter.write(GenomicRun, ...)` materialised one
  `AlignedRead` per record, which decoded cigar + mate triple
  for every read — fields FASTQ does not need. Now pre-fetches
  the whole `sequences` + `qualities` byte buffers + read-names
  list once, slices in-memory. 11K reads/s → 265K reads/s on
  1M reads × 100bp (commit `ae9441d`).
- **Java + ObjC FASTQ writer bulk-fetch parity** mirrors the
  same fix. Java now sustains ~750K reads/s on the same 1M
  workload; ObjC ~635K reads/s (commit `0f99852`). New
  microbenches: `FastqBulkBenchTest` (Java, opt-in via
  `-DTTIO_FASTQ_BENCH=1`) + `objc/Tools/obj/TtioFastqBench`.
- **Java transport genomic encode 33 → 235K reads/s (+612%)** at
  100K reads × 100bp (commits `758b340` + `701f310`). Two
  fixes: (1) memoised `GenomicRun.isMateInfoInlineV2()` after
  noticing the probe reopened the `mate_info` HDF5 group 3× per
  record (300K group opens at 100K reads, ~2.2s of pure
  framework overhead — ObjC already cached this via
  `_mateInfoLinkType`), and (2) eager-cached the M82 compound
  fall-through of `cigarAt` + bypassed per-record `AlignedRead`
  materialisation in `TransportWriter`. At 1M reads: 328K rps.
  Java is now ~40% faster than ObjC on this workload. New
  microbench: `TransportEncodeBenchTest`.
- **ObjC transport genomic encode 70 → 164K reads/s (+134%)**
  (commit `d5e2e25`). Same per-record `[grun readAtIndex:i]` →
  `dataUsingEncoding:` re-encode roundtrip Java had. Now bulk-
  fetches `wholeSequencesData` / `wholeQualitiesData` /
  `allReadNames` once. New microbench: `TtioTransportEncodeBench`.
- **Python + ObjC `cigar_at` eager-cache (parity with Java)**
  (commit `7ac32e4`). M82 compound fall-through no longer
  re-decodes per call. Python: 10 → 36K rps at 100K reads
  transport encode (+260%). ObjC: net-neutral on the test
  fixture but protective for any NSData-typed compound returns.
- **Byte-channel cache audit** (`GenomicRun.byteChannelSlice` /
  `-byteChannelSliceNamed:`). The codec-compressed path cached;
  the uncompressed path returned the raw HDF5 buffer per call,
  so `sequencesFull` / `-wholeSequencesData` warmups were
  silently a no-op for files written with `signal_compression =
  NONE`. Fixed in both Java + ObjC (commit `221611c`).

### Cross-language parity gates

- **nmrML 3-way probe parity** — `NmrMLProbe.java`,
  `TtioNmrMLProbe.m`, and a Python harness drive all three
  readers against synthetic + `bmse000325` inputs and assert
  bit-exact JSON for `numberOfScans` /
  `spectrometerFrequencyMHz` / `fidReal` / `fidImag` (commit
  `8c6b8b0`).
- **mzML 3-way probe parity** — `MzMLProbe.java`,
  `TtioMzMLProbe.m`, and `test_mzml_cross_lang_parity.py` cover
  synthetic + `tiny.pwiz.1.1.mzML` fixtures with full mz +
  intensity arrays plus precursor / polarity / RT scalars.

### mzML reader bug-fixes

- **Java**: `endElement("binary")` no longer skips spectra with
  empty `<binary></binary>` arrays (PSI-MS reference fixture
  intentionally tests this via its "spectrum with no data"
  userParam). 4-spectrum `tiny.pwiz` now reports 4 spectra (was
  3) (commit `01ca2b4`).
- **Python + ObjC**: `<referenceableParamGroupRef>` is now
  resolved. Both readers buffer cvParams under each
  `<referenceableParamGroup id="...">` and replay them when a
  spectrum / chromatogram cites the group via
  `<referenceableParamGroupRef ref="...">`. Polarity (and any
  other CV param) on referenced groups now reaches the
  per-spectrum surface (same commit).

### Documentation

- `docs/cross-language-matrix.md` gains entries for the two new
  probe-style parity tests + CLI inventory rows for the four
  new probe binaries.
- `docs/benchmarks/2026-05-05-v1.0-comprehensive-perf-report.md`
  §11 documents the post-v1.0.0 perf tweaks with reproducer
  invocations.

---

## [v1.0.0] — 2026-05-04 — first stable release

This is the first stable release of TTI-O. The format string is
`ttio_format_version = "1.0"`; container ABI, codec wire formats,
encryption envelope, and digital-signature canonicalisations are
contractually frozen at this point. Pre-v1.0 development was never
publicly released; that history lives in `git log`.

### Format

- HDF5-backed `.tio` container; opaque `study/` group with per-modality
  child groups (`ms_runs/`, `genomic_runs/`, `chromatograms/`,
  `nmr_runs/`, `image_cubes/`, …).
- Deterministic write order; the Python, Java, and Objective-C
  reference implementations all produce byte-identical output for the
  same input. Cross-language byte-equality is part of the contract,
  not a coincidence — see `pytest -m integration`.
- Feature-flag preamble (`ttio_features` JSON array attribute) for
  forward-compatible optional capabilities. ISA-Tab investigation
  linkage on every container.

### Codecs

| Id | Symbol               | Description                                          | Channels                                              |
|---:|----------------------|------------------------------------------------------|-------------------------------------------------------|
| 0  | NONE                 | Passthrough                                          | any                                                   |
| 1  | ZLIB                 | HDF5 deflate filter (level 6 default)                | any                                                   |
| 2  | LZ4                  | HDF5 filter id 32004 (~35× faster write than zlib)   | any                                                   |
| 3  | NUMPRESS_DELTA       | Numpress + delta encode (sub-ppm lossy)              | numeric MS m/z channels                               |
| 4  | RANS_ORDER0          | rANS order-0 entropy coder                           | sequences / qualities / cigars / integers             |
| 5  | RANS_ORDER1          | rANS order-1 entropy coder                           | sequences / qualities / cigars / integers             |
| 6  | BASE_PACK            | 2-bit ACGT pack with sidecar mask for IUPAC bases    | sequences                                             |
| 7  | QUALITY_BINNED       | Illumina-8 binning (lossy, CRUMBLE-derived)          | qualities                                             |
| 11 | DELTA_RANS_ORDER0    | Delta + rANS-O0                                      | sortable integer channels                             |
| 12 | FQZCOMP_NX16_Z       | CRAM-mimic adaptive quality (V4 only, magic `M94Z`)  | qualities                                             |
| 13 | MATE_INLINE_V2       | Inlined mate_info v2 (single channel)                | mate_info compound                                    |
| 14 | REF_DIFF_V2          | Reference-diff v2 (slice-based, embedded reference)  | sequences                                             |
| 15 | NAME_TOKENIZED_V2    | 8-substream multi-token columnar codec               | read_names                                            |

Ids 8, 9, 10 are reserved on the wire (Java enum ordinal stability)
but carry no live codec. Reader paths reject them with migration
errors. Codec wire formats are documented in `docs/codecs/*.md`;
per-channel pipeline wiring is documented in `docs/format-spec.md`
§10.4–§10.10.

### Modalities

- Mass spectrometry: LC-MS, MS-image cubes, ion mobility, profile +
  centroid spectra.
- Nuclear magnetic resonance: 1-D and native 2-D (HSQC, COSY, NOESY).
- Vibrational imaging: Raman, IR.
- UV-Vis spectra.
- Two-dimensional correlation spectroscopy (2DCOS).
- Chromatograms.
- Genomic alignment runs: full BAM/CRAM importer parity, per-record
  metadata, codec-aware channel wiring.

### Format I/O — FASTA / FASTQ

- **FASTA importer**: `FastaReader` reads reference genomes (for
  embedding at `/study/references/<uri>/`, paired with BAM/CRAM
  input) or unaligned reads (panels, target lists, quality-stripped
  reads → `WrittenGenomicRun` with SAM unmapped sentinels). gzip
  auto-detected via magic bytes regardless of extension.
- **FASTA exporter**: `FastaWriter` writes a `ReferenceImport` or a
  `WrittenGenomicRun` to FASTA with configurable line wrap
  (default 60 chars) and a samtools-compatible `.fai` index
  emitted alongside.
- **FASTQ importer**: `FastqReader` parses 4-line records into
  unaligned `WrittenGenomicRun` instances. Phred offset is auto-
  detected (`33` modern Illumina / Sanger vs `64` legacy
  Illumina); detected source recorded for round-trip planning.
  Internal storage normalises to Phred+33.
- **FASTQ exporter**: `FastqWriter` writes a run to FASTQ with
  Phred+33 default and Phred+64 selectable. The `0xFF` "qualities
  unknown" sentinel is mapped to Phred 0 (`!`) on output so the
  result is always parseable.
- **Cross-language byte equality**: Python, Java, and ObjC produce
  byte-identical FASTA + FASTQ output for the same input — proven
  by the `test_fasta_fastq_cross_language.py` 3-way harness.
- **CLIs**: Python `python -m ttio.tools.{fasta,fastq}_{import,export}_cli`,
  Java `FastaRoundTrip` / `FastqRoundTrip`, ObjC `TtioFastaRoundTrip`
  / `TtioFastqRoundTrip`.

### Encryption + signing

- **Per-AU encryption** (AES-256-GCM) on signal-channel datasets and
  compound-metadata payloads. Versioned wrapped-key blob carries DEK
  rotation history; envelope decryption supported via local key, KMS,
  or user-supplied callback.
- **Digital signatures**: HMAC-SHA256 (canonical) plus post-quantum
  ML-DSA via liboqs. Signatures verify identically across all three
  reference implementations.

### Language bindings

- **Python** (`pip install ttio`): full read/write/encryption/sign;
  ctypes wrapper for the native rANS / v2-codec library.
- **Java** (Maven Central `global.thalion:ttio`): full parity; JNI
  wrapper for the same native library.
- **Objective-C** (GNUstep): full parity; native library linked
  directly. `objc/Tools/MakeFixtures` produces the canonical
  cross-language reference fixtures.

### Cross-language guarantee

Byte-equal output for shared codec paths under the test corpora in
`data/genomic/` (NA12878 chr22, NA12878 WES, HG002 Illumina 2×250,
HG002 PacBio HiFi subset). Verified on every commit via
`pytest -m integration`; SHA-256 hashes match Python ↔ Java ↔
Objective-C.

### Native library

`libttio_rans` (CMake / clang) ships the v2 codec kernels (rANS,
ref_diff_v2, mate_info_v2, name_tokenized_v2, fqzcomp_nx16_z V4).
Mandatory at runtime for genomic-run write/read on all three
language bindings (`TTIO_RANS_LIB_PATH` env var, or `libttio_rans.so` /
`.dylib` / `.jni` on the loader search path).

### Transport — genomic bulk mode (Phase 2c-T)

- New packet types `BlobV2MateInfo` (0x09), `BlobV2RefDiff` (0x0A),
  `BlobV2NameTok` (0x0B) carry the verbatim v2 codec blobs
  (`mate_info/inline_v2`, `sequences/refdiff_v2`,
  `read_names/name_tok_v2`) on the wire. See
  `docs/transport-spec.md` §3.2 / §4.10–4.12 / §5.8 / §6.4.
- Stream-level feature flag `bulk_mode_v2_blobs` (required, no
  `opt_` prefix). Receivers that cannot honor verbatim blob
  injection refuse the stream.
- CLIs accept `--bulk` on encode: Python
  `python -m ttio.tools.transport_encode_cli --bulk`, Java
  `TransportEncodeCli --bulk`, ObjC `TtioTransportEncode --bulk`.
- Cross-language byte-identity verified by the 9-cell
  Python/Java/ObjC matrix in
  `python/tests/validation/test_phase_2c_t_bulk_mode.py`.
- Storage-provider parity: HDF5, memory://, sqlite://, and zarr://
  write paths all honor `bulk_v2_blobs` and write the verbatim
  blob bytes — see
  `python/tests/validation/test_phase_2c_t_storage_providers.py`.
- Measured speedup: receiver-side decode runs **1.36×–1.43× faster**
  in bulk mode (10K and 50K-read fixtures); encode is near-parity
  (≤3% delta), wire size grows ~4% from the additional blob
  packets. See
  `docs/benchmarks/2026-05-05-phase-2c-T-bulk-mode.md`.

### Cross-language nmrML reader parity

The Python `ImportResult` gained four nmrML acquisition-parameter
fields on 2026-05-05 (`spectrometer_frequency_mhz`,
`number_of_scans`, `fid_real`, `fid_imag`); Java and ObjC sibling
readers now surface the same fields:

- **Java** `NmrMLReader.NmrMLResult` exposes
  `spectrometerFrequencyMHz()`, `numberOfScans()`, `fidReal()`,
  `fidImag()`. The parser also accepts
  `<irradiationFrequency value="...">` directly inside
  `<acquisitionParameterSet>` (matches Python; previously
  required `<directDimensionParameterSet>`).
- **ObjC** `TTIONmrMLReader` now also exposes deinterleaved
  `fidReal` / `fidImag` `NSData` properties alongside the
  pre-existing `spectrometerFrequencyMHz` / `numberOfScans`.

### Native codec fixes

- **NAME_TOKENIZED_V2** (codec id 15): decoder MATCH path no
  longer rejects valid encoded blobs whose pool entry has a
  different total token count than the block's column shape.
  Surfaced by the production-corpus decode benchmark on real
  Illumina BAMs with mixed flowcell prefixes (e.g. `H2YHMBCXX`
  tokenises to 3 tokens vs `H2YT5BCXX`'s 5 tokens because of
  the internal digit). Permanent regression guard at
  `python/tests/test_name_tokenizer_v2_native.py::test_mixed_flowcell_token_count_regression`;
  minimal failing fixture preserved at
  `python/tests/fixtures/codecs/name_tok_v2_corrupt_94.txt`.
- **Bruker TDF importer**: `frame2retention_time` returned a 1-D
  ndarray under opentimspy ≥ 1.2 even for scalar input; the
  per-frame list comprehension produced a `(n, 1)` 2-D
  `retention_times` buffer that the writer rejected. Now passes
  the whole `frame_ids` array at once.
- **Zarr 3.x empty-chunk**: `ZarrProvider.create_dataset` now
  clamps chunk dims to ≥ 1 so empty datasets (length == 0) build
  cleanly under zarr-python 3.x.

### Format support — Bruker .tsf

`ttio.importers.bruker_tdf.read_metadata()` now recognises both
`analysis.tdf` (TIMS) and `analysis.tsf` (non-TIMS Bruker QTOF /
MALDI) `.d` directories. The SQLite metadata schema is shared,
so frame counts / retention times / instrument-vendor strings
parse identically. Full per-frame `read()` is TDF-only in v1.0
(opentimspy is TDF-only); calling it on a `.tsf` directory raises
`BrukerTDFUnavailableError` with a pointer at the
`msconvert` → mzML workaround.

### Performance — benchmark suites + microbench tooling

Three new perf harnesses for release-to-release tracking:

- **`python/tests/stress/test_fasta_fastq_benchmark.py`** —
  five scenarios per fixture size (FASTQ export / import,
  FASTA export / import, FASTQ → `.tio` → FASTQ round-trip)
  at 1K, 10K, plus opt-in 100K and 1M reads via
  `TTIO_INCLUDE_LONG_TAIL=1`.
- **`python/tests/stress/test_production_corpus_benchmark.py`** —
  BAM → `.tio` → decode-all-reads cycle against the real corpora
  under `data/genomic/` (synthetic, na12878 chr22, na12878 WES,
  hg002 Illumina subset, hg002 PacBio). The full 1.6 GB
  hg002 chr22 BAM is opt-in via `TTIO_INCLUDE_FULL_CORPUS=1`.
- **`global.thalion.ttio.tools.Benchmark` (Java) +
  `TtioBenchmark` (ObjC)** — pair with the Python harnesses to
  give cross-language perf-tracking parity. Both emit the same
  JSON schema so a single `jq` over the three result files
  diffs cleanly.

Headline numbers consolidated in
`docs/benchmarks/2026-05-05-v1.0-comprehensive-perf-report.md`.

### Build + dev workflow

- **`scripts/dev-setup.sh`** — one-shot Python developer setup:
  builds `libttio_rans.so` + JNI wrapper, installs the package
  with the broadest test extras, prints required env-var
  exports. PEP 668 (Ubuntu 24.04+ externally-managed) friendly.
- **`scripts/fetch-vendor-fixtures.sh`** — downloads + sha256-
  verifies the public Thermo `small.RAW` (MIT, ~1.5 MB) and
  Bruker `diaPASEF.d` (Apache-2.0, ~1 MB) fixtures from upstream
  repos. Manifests pinned at `data/vendor/{thermo,bruker}/*.sha256`.
- **CI** (`.github/workflows/ci.yml`):
  - All Python jobs (`python-test`, `python-validation`,
    `python-stress`) now build the native rANS library + JNI
    wrapper before running pytest. Previously skipped ~90 v2-codec
    tests silently because the runtime library was never present.
  - New `python-vendor-fixtures` job exercises the Bruker TDF +
    Thermo `.raw` integration paths on every push/PR (mono +
    ThermoRawFileParser v1.4.5 + sha256-pinned fixtures).
- **`CONTRIBUTING.md`** added — top-level entry point with
  quick-start, repo layout, per-language test commands, optional
  fixture flow, code style notes.
