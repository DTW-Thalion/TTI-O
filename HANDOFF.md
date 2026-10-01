# HANDOFF — M101 SAM optional tags

**As of 2026-10-01.** M101 closes the first gap on WORKPLAN's
"Take these first" list. The BAM, SAM and CRAM importers kept SAM
fields 1–11 and discarded the optional tags in all three SDKs, so a
`.tio` made from a BAM was not a lossless copy of it. The tags now ride
a genomic `tags` channel coded with SAM_TAGS (codec id 18), a CRAM-style
column codec in the shared native library:

- a tag-line dictionary per blob, and one column per tag key;
- MD:Z / NM:i recomputed from the reference;
- repeated integers stored as back-references.

Phase 0 proved it on the GIAB HG002 2x250 chr22 slice before any SDK
code: 1.93 B/read against 2.57 for CRAM 3.1 `small`, byte-exact over
10.6 M reads. (The Phase 0 CRAM baseline first read 5.61 B/read: its
chr22-only reference made samtools embed a reference built from the
reads and keep MD/NM. Corrected 2026-10-01. On bwa-aligned NA12878 WES
CRAM is ahead on tags, 1.27 against 3.24 B/read, because of `XA:Z`.) Format-spec §10.13, `docs/codecs/sam_tags.md`, binding
decisions §96–§100. Branch: `m101-sam-optional-tags`.

| Task | Scope | Status | Spec proof |
|---|---|---|---|
| **Phase 0** | Prototype column codec (`tools/prototypes/m101_sam_tags/`) over the HG002 chr22 slice; CRAM 3.1 normal/small baselines via `samtools cram-size`; entropy check of the integer columns; O0-only vs O0+O1 | ✅ 2026-10-01 | required before the codec design froze |
| **Native kernel** | `native/src/sam_tags.{c,h}`, `ttio_sam_tags_encode/_decode/_free` in `ttio_rans.h`. rANS-O0 substreams with a sparse codec-4 table. `test_sam_tags` passes under ASan/UBSan: calmd cases, canonical ints, every-prefix and trailing-byte rejection, determinism | ✅ | wire spec written first |
| **A Python** | `codecs/sam_tags.py` + registry adapter, `WrittenGenomicRun.tags`, `AlignedRead.tags`, BAM/SAM/CRAM import, whole-run writer, blocks_v1 writer, index triple and late upgrade, readers, exporters, `opt_sam_tags`, per-AU walkers, signatures, plaintext and encrypted transport, `bam_dump` `tags` key | ✅ `test_m101_sam_tags.py` 12 + `test_sam_tags_native.py` 5 | — |
| **B ObjC** | Mirror A: `TTIOSamTags` over the kernel, `TTIOBamReader` keeps `fp[11]`, writer/reader/exporter/per-AU/transport/signatures, `TtioBamDump` `tags` | ✅ `TestM101SamTags` 76 checks | — |
| **C Java** | Mirror A: JNI entry points, `BamReader.addRecord` formats `getAttributes()` as samtools text (float and `B` arrays included), `WrittenGenomicRun.tags` record component, writer/reader/`BamWriter`/per-AU/transport/signatures, `BamDump` `tags` | ✅ `M101SamTagsTest` 17 tests | `SamTagText` matches `samtools view` on every tag type and a seeded float sweep (`%g`, half-even on the exact binary value) |
| **D Conformance** | 3×3 writer × reader matrix over a tagged fixture (REF_DIFF with derivation, no reference, tags first appearing in a later block), `bam_dump` cross-language JSON over a tagged SAM fixture, per-AU encryptor × decryptor cells | ✅ `test_m101_sam_tags_matrix.py` 29/29 | — |
| **E Release + re-measure** | WORKPLAN item 2: cut the release and re-measure NA12878 WES chr22, HG002 chr22 and whole HG002 against CRAM 3.1 | chr22 slices ✅ (`docs/benchmarks/2026-10-01-m101-cram-remeasure.md`: WES 44.7 MB vs CRAM small 33.3 MB, HG002 chr22 1,020.5 MB vs 834.6 MB); whole HG002 and the release ⏳ | CRAM reference must cover every header contig |
| **F Docs** | format-spec §10.4 row 18, §10.12.2 triple, §10.12.6, §10.13; codec doc; feature-flags; transport-spec 4.3.1 / 4.24; genomic-runs; vendor-formats; binding decisions §96–§100; CHANGELOG; this file | ✅ | — |

**Fidelity contract:** decode returns each read's tag text exactly as
the importer saw it. Two limits on that:

- **BAM integer widths are not kept.** `samtools view` prints every
  integer type as `i`, so the width is chosen again when a BAM is
  written.
- **MD/NM are recomputed only in REF_DIFF_V2 blocks.** There the reader
  already needs the reference. Every other block stores them like any
  other tag.

**Known limits, deliberate:**

- **Whole-channel runs with tags refuse per-AU and region encryption**
  rather than leave MD strings in plaintext; rewrite them as blocks_v1.
- **No MD transform without a reference.** Without a reference, MD costs
  5.1 B/read, above CRAM's bzip2 column; an MD-specific transform is a
  follow-up.
- **Java CRAM tags are tag-sorted.** htsjdk's CRAM records have neither
  aux bytes nor text, so Java takes their tags from the tag-sorted
  attribute map; Java BAM export writes `H` tags as `Z`.
- **Encrypted transport carries no reference bytes.** It does not carry
  embedded reference bytes, so a received tagged REF_DIFF container
  restores through `REF_PATH` (the M99.1 limit).

**Suite state at handoff (M101):**

- **Python:** 2312 pass, 0 M101 failures.
- **ObjC:** 5374 pass, 0 fail.
- **Java:** 1664 tests, 0 failures, 18 skipped.
- **Cross-language matrices** (M82–M101, run from a space-free mirror):
  343 pass.

**Local environment:** the Python suite runs under WSL Ubuntu, with
venv `~/ttio-venv` and a Linux build of `native/`
(`TTIO_RANS_LIB_PATH`). The Windows Python build lacks zlib and the
MinGW runtime. `scripts/setup-objc-wsl.sh` installs the ObjC
toolchain the way CI does (needs sudo).

---

## When to overwrite this file

This `HANDOFF.md` is replaced *per active milestone* — the git
history (`git log -- HANDOFF.md`) shows that pattern (M81 →
M82 → … → M99 → this). When the next multi-language milestone kicks
off, overwrite this file with the milestone's plan + task table;
otherwise small post-v1.0 follow-ups go to PRs + CHANGELOG only.

For ongoing work not coordinated through HANDOFF, see:

- `CHANGELOG.md` § `[Unreleased]` — what's landed since the last
  tag.
- `WORKPLAN.md` — milestone history + binding decisions (§96–§100
  are M101's).
- `tti-workbench-server` repository — daemon-side workstreams.
