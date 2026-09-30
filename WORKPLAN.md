# TTI-O Workplan

Milestones with concrete deliverables and acceptance criteria.
Historical milestones (M1–M8) use the original MPEG-O class names
as a record of what was built; current milestones use TTI-O names.

---

## Take these first when work resumes (2026-09-30)

These gaps decide whether a .tio can stand in for BAM or CRAM 3.1 on
aligned reads, so they come before new features.

1. **Carry SAM optional tags.** The BAM and SAM importers keep fields
   1-11 and discard the rest (Gotcha §152) in Python, Java and ObjC,
   so a .tio made from a BAM is not a lossless copy of it. On the
   HG002 2x250 GRCh38 BAM the tags are about 11.3% of each record:
   PG 15.0 bytes, MD 13.6, PQ 7.8, AS and UQ 7.6, SM 7.4, AM 7.2,
   NM 6.9. This needs a format-spec channel for them, a codec choice
   (MD and NM can be recomputed from the reference and the alignment,
   as CRAM does), import and export in Python, Java and ObjC, and a
   conformance round trip.
2. **Release and re-measure REF_DIFF_V2 with unmapped reads.** The fix
   is in `[Unreleased]`: the kernel carries a CIGAR `*` read in the UL
   substream and the writers no longer fall back to BASE_PACK. The
   sizes on record predate it: on NA12878 WES chr22 TTI-O was 64.1 MB
   against BAM 66.4 MB and CRAM 3.1 small 32.0 MB, and on HG002 2x250
   chr22 1,000 MB with GRCh38 against CRAM 3.1 small 825 MB. After the
   release, measure both slices and the whole HG002 BAM against CRAM
   3.1 again. The comment above the BASE_PACK fallback in
   `SpectralDatasetGenomicWriter.java` still names unmapped reads as a
   cause of it, and they no longer are.

---

> **Status (2026-05-24).** v1.0.0 shipped 2026-05-04. The
> three-language reference implementation (ObjC normative,
> Python `ttio`, Java `global.thalion.ttio`) is feature-frozen
> for the v1.0 surface; post-v1.0 work appears in CHANGELOG
> `[Unreleased]` and ships as small per-PR follow-ups rather
> than coordinated multi-language milestones. Highlights since
> v1.0.0:
>
> - **TTI-O Workbench Client (W1–W6)** shipped through PR #95
>   (auth + transport) → PR #126 (W6.6 SDK reference docs +
>   quickstart). Surface: `ttio` Python CLI umbrella + Python SDK
>   (`python/src/ttio/workbench/`) + Java SDK
>   (`java/.../workbench/`) + JavaFX desktop GUI `tio-browser`
>   panels (Connection Manager / Container Browser / Transfer
>   Manager + Selective Access / Cohort Query Builder / Pipeline
>   Launcher + Job Monitor / Interactive Session Launcher /
>   Encoding + Export). `tio-browser` rebumped to **v1.5.0** at
>   W5.7 (#108). PQC + BYOK + federation client capabilities + W6
>   live-acceptance smoke green.
> - **FD-1 multi-recipient encryption (Phases A / B / C-0 /
>   C-2a)** shipped 2026-05-21 → 2026-05-22 across all 3
>   languages with cross-language conformance vectors at
>   `conformance/multi_recipient/`. Append-only on §4.4 →
>   single-recipient packets stay byte-identical, so existing
>   BYOK / envelope / PQC fixtures are unchanged. Server-side
>   FD-1 D+ continues in the `tti-workbench-server` repo.
> - **Per-AU decrypt-in-place** parity APIs in all 3 languages
>   (PRs #161 + #162 + **#163**, 2026-05-23) — closes the
>   long-standing "silent no-op on per-AU files" gap on the
>   legacy `decryptInPlace` path. PR #163 (`5462489d`) brought
>   genomic-run signal-channel coverage in line with MS, so all
>   three languages now walk `study/genomic_runs/<name>/signal_channels/`
>   in addition to `ms_runs/<name>/signal_channels/` and
>   unconditionally clear `opt_per_au_encryption` /
>   `opt_encrypted_au_headers` / root `@encrypted` after a
>   successful decrypt. Unblocks `tti-workbench-server` FD-1 D-1
>   pipeline step.
> - **Vibrational JCAMP-DX `.tio` materialization parity**
>   (parity-audit 3.1) shipped through PRs #144–#146 (Python +
>   Java + ObjC).
> - **CI** bumped to JDK 25 + GitHub Actions Node 24 runtimes
>   (PR #92 + #142). Minimum supported runtime remains
>   JDK 17+ for the library, JDK 22+ for `tio-browser` (JavaFX).
>
> For the active milestone-style handoff (or the absence
> thereof), see [`HANDOFF.md`](HANDOFF.md).

---

> **Status (2026-05-02).** Phase 6 (M89 transport extension, M90
> encryption/anonymisation, M91 multi-omics integration) **shipped**
> alongside V- and C-series debt repayment. Phase 8 (post-M91
> abstraction polish — `Run` protocol, modality-agnostic
> `runs` / `runsForSample` / `runsOfModality`, mixed-dict write API,
> per-run provenance compound dual-write/dual-read) **shipped**.
> **Phase 9 (M93/M94/M94.Z/M95 codec parity) shipped 2026-04-30.**
> **Phase 10 — FQZCOMP acceleration + ObjC M44 catch-up shipped
> 2026-05-01:** removed legacy FQZCOMP_NX16 (M94 v1, codec id 10);
> built `libttio_rans` C library at `native/` (AVX2/SSE4.1/scalar
> SIMD dispatch, pthread thread pool, V2 multi-block wire format,
> ~510 MiB/s encode / ~605 MiB/s decode); wired V2 dispatch through
> Python ctypes / Java JNI / ObjC direct linkage with `prefer_native`
> opt-in; refactored ObjC writer chain (TTIOAcquisitionRun,
> SpectrumIndex, Spectrum + 7 subclasses, InstrumentConfig,
> SignalArray, CompoundIO, SpectralDataset.writeToFilePath:) to
> accept `id<TTIOStorageGroup>` — closing the 2026-04 gap with
> Python's M44 and Java's M44 migrations. ObjC MS-only datasets
> now write through memory:// / sqlite:// / zarr:// URLs via either
> the class-method (Task 30) or instance-method (Task 31) entry
> point. Test counts: Python 540/540, Java 845/845, ObjC 3256/2
> (the 2 are pre-existing TestMilestone29 Thermo mock-binary tests,
> unrelated). NMR runs and Image-subclass datasets continue to
> require HDF5 (H5DSset_scale dimension scales / native 3D cubes
> have no protocol equivalents — same scope as Python and Java).
> M93 REF_DIFF (codec id 9), M94 v1 FQZCOMP_NX16 (codec id 10),
> M94.Z FQZCOMP_NX16_Z (codec id 12, CRAM-mimic), and M95
> DELTA_RANS_ORDER0 (codec id 11, delta + zigzag + varint + rANS
> order-0 for sorted-ascending integer channels) all byte-exact
> across Python / ObjC / Java. Cython acceleration of REF_DIFF
> (44× pack / 35× unpack) and NAME_TOKENIZED (11% / 8% chr22 drop)
> ships alongside. Pipeline defects fixed
> (`_mate_info_is_subgroup` caching + per-read codec-import hoist).
> chr22 encode 18 min → 27.91 s (38.7×); decode 24.6 min → 21.76 s
> (67.8×); now 9.2× / 13.4× off CRAM 3.1, down from 355× / 1162×
> pre-session. Codec compute is no longer the bottleneck; remaining
> gap is HDF5 framework + multi-omics infrastructure, accepted for
> v1.2.0 scope per user direction. All codec-tier items complete;
> M92 release prep follows.
>
> **Phase 11 — cross-language perf sweep shipped 2026-05-01.** Five
> focused wins, no wire-format breaks. Tasks #78–#81 + #83 from the
> v1.2.0 codec parity workplan all closed in a single session:
> * **Task #83 — Python rANS Cython accelerator.** New
>   `python/src/ttio/codecs/_rans/_rans.pyx` (byte-exact mirror of
>   `rans.py`) drops rans_o0 17–29× and cascades 6.8× into REF_DIFF and
>   2× into DELTA_RANS — both internally call `rans.encode/decode`.
>   Python rans_o0 now runs faster than ObjC's hand-rolled C.
> * **Task #81 — Native M94.Z decode entry.** New
>   `ttio_rans_decode_block_m94z` C entry point bakes the M94.Z context
>   formula directly into native code. C decode kernel ~107 MiB/s vs
>   Cython's ~96 MiB/s; full V2 wrapper still has ~110 ms metadata-setup
>   overhead in Python (follow-up). Java JNI + ObjC linkage to the new
>   entry deferred.
> * **Task #80 — Python lazy zarr import.** `discover_providers()`
>   eager `from .zarr import ZarrProvider` was loading zarr (~135 ms)
>   on every non-HDF5 first call, inflating ms.memory write to 167 ms.
>   Made lazy via `_try_register_zarr()`. Result: ms.memory write
>   167 ms → 5.5 ms (30×).
> * **Task #79 — Java jcamp reader.** Replaced
>   `HashSet<Character>.contains` autobox-per-char in
>   `JcampDxDecode.hasCompression` with `boolean[128]`; pre-allocated
>   `double[]` from NPOINTS hint and replaced
>   `ArrayList<Double> + line.split("\\s+") + Double.parseDouble`
>   with a manual whitespace tokenizer in `JcampDxReader`. jcamp:
>   280 ms → 134 ms (2.1×).
> * **Task #78 — Java FQZCOMP_NX16_Z encode.** Per-stream chunk buffers
>   as `short[]` (one store per renorm replaces two byte stores);
>   pre-padded qualities to drop the `(i < n)` branch from the hot
>   loop. Warm encode: ~25 MB/s → ~34 MB/s (36% faster).
>
> Test counts after Phase 11: Python 1800/1800, Java 879/879, native
> ctest 5/5. chr22 ratio against CRAM 3.1 unchanged at 1.965× — the
> perf wins are at the codec layer; the residual gap is still HDF5
> multi-omics framing (Task #82, deferred as multi-week scope).
>
> **Phase B.2 — L2 adaptive M94.Z V3 (Range Coder) infrastructure
> shipped 2026-05-02; encode-perf follow-up landed same day.**
> Initial V3 wrapper at HEAD `72eb845`; numpy-vectorized context
> derivation + first-encounter pass at the next commit drops V3
> encode wall **124.79 s → 25.83 s (4.8×)** on chr22, beating V1's
> ~28 s. Output byte-exact, 1811 Python tests green. Range-Coder
> pivot during
> Phase 2 (rANS state-range constraint blocked variable-T close to
> 2¹⁶; see memory `rans_nx16_variable_t_invariant`); native C kernel
> + V3 wire format + Python ctypes wrapper all green (8/8 native
> ctests, 553/553 Python M94.Z tests). **Phase 4 chr22 hard gate:
> FAILED** — measured 113.33 MB / 1.316× CRAM (target ≤ 99 MB / 1.15×),
> over by 14.32 MB. Qualities 69.34 MB at 0.393 B/qual (vs V1: 69.73 MB
> at 0.395 B/qual) — adaptive moves the needle by ~0.4 MB / 0.002 B/qual.
> Per-symbol adaptive freqs alone, on top of the existing context
> formula (`prev_q × pos_bucket × revcomp`, sloc=14), do not move the
> needle on chr22 — block sizes already let static freqs converge to
> the empirical distribution.
> Per spec §10, this **blocks Phase 5 (Java JNI) and Phase 6 (ObjC)**
> wrappers; the V3 infrastructure stays Python-only until a model
> beats V1 substantially. The Range-Coder kernel itself is reusable.
> Next step is **Task #84** (below): brainstorm a richer context
> model before any further codec work. See
> `docs/benchmarks/2026-05-01-chr22-byte-breakdown.md` §8.
>
> * **Task #84 — Richer-context M94.Z (Stage 1 + Stage 2 done
>   2026-05-02; Stage 3 done 2026-05-03; V4 SHIPPED across all 3 langs).**
>
>   **Stage 3 (Java JNI + ObjC V4 dispatch) SHIPPED 2026-05-03.** V4
>   reaches feature parity across all 3 reference implementations:
>   Python (ctypes), Java (JNI through `libttio_rans_jni`), and
>   Objective-C (direct link). All three call the same deterministic
>   `ttio_m94z_v4_encode` C entry, so the resulting M94Z V4 streams are
>   byte-identical end-to-end across all 4 benchmark corpora — proven
>   by the cross-language Python ↔ Java ↔ ObjC integration test
>   `python/tests/integration/test_m94z_v4_cross_language.py` (4
>   parametrized cases, all PASS). V4 is now the default emit format in
>   each language when its respective native library is loaded; V1/V2
>   read-compat preserved as the no-native fallback.
>   Follow-ups (out of scope, separate tasks): Java/ObjC V4 perf
>   benchmarks, V4 empty-input handling (V3 fallback at the Python
>   layer was deferred from Stage 2 Task 12).
>
>
>   **Stage 2 (V4 CRAM 3.1 fqzcomp port) SHIPPED 2026-05-02.** The
>   native C port of `htscodecs/fqzcomp_qual.c` lives in
>   `native/src/{rc_cram,fqzcomp_qual,m94z_v4_wire}.{c,h}`; htscodecs
>   SHA `7dd27f4` is vendored at `tools/perf/htscodecs/` for test-time
>   byte-equality only (no runtime link). The Python ctypes wrapper at
>   `python/src/ttio/codecs/fqzcomp_nx16_z.py` exposes
>   `_encode_v4_native` / `_decode_v4_via_native`; **V4 is the default
>   encoded format when `_HAVE_NATIVE_LIB`**, V3 stays as the no-native
>   fallback and read-compat path. V4 byte-equality with htscodecs is
>   guaranteed across all 4 benchmark corpora. Per-corpus B/qual:
>   - chr22 NA12878 100bp WGS: **0.365 B/qual**
>   - NA12878 WES (capture, ~95 bp): **0.273 B/qual**
>   - HG002 Illumina 2×250: **0.259 B/qual**
>   - HG002 PacBio HiFi (~14K reads, 18.5 kb mean): **0.398 B/qual
>     (4% better than V3 — V4 actually beats V3 on PacBio HiFi
>     despite Phase 0's auto-tune-doesn't-save-PacBio finding;
>     fqzcomp_qual's per-block strategy selection captures structure
>     the bit-pack candidates missed)**.
>   chr22 V4 encode wall is **3.4× faster** than V3 end-to-end. Test
>   counts: 11/11 native ctests, 48/48 V4-related Python tests, full
>   1800/1800 baseline. Java JNI + ObjC wrappers are **Stage 3
>   follow-up** — Python is the load-bearing V4 path; the C plumbing
>   is proven byte-exact and ready for binding without further model
>   work. Results at
>   `docs/benchmarks/2026-05-02-m94z-v4-stage2-results.md`; V4 wire
>   format documented in `docs/codecs/fqzcomp_nx16_z.md` §2.
>
>   **Stage 1 (cross-corpus prototype, 2026-05-02 — historical).**
>   Stage 1 prototype harness at `tools/perf/m94z_v4_prototype/`
>   measured 5 candidate context-model designs on three Illumina
>   corpora: chr22 NA12878 (100 bp WGS), NA12878 WES (capture, ~95 bp),
>   HG002 Illumina 2×250 (250 bp WGS, 1 M-read subset of 10.6 M-read
>   chr22). All candidates miss the 1.15× CRAM gate on chr22; the
>   gate is chr22-v1.2.0-specific and not directly applicable to other
>   corpora. **Cross-corpus signal:**
>   - **c3** (length-heavy: 8 bit full-Phred prev_q + 4 bit pos +
>     4 bit length + 1 bit revcomp, sloc=17) is the cross-corpus
>     winner — outright on WES and HG002 2×250, within 0.7% of c2
>     on chr22.
>   - **c2** (equal-precision 4+4+4 prev_q, drop length, sloc=17)
>     wins only on chr22 — its uniform-100bp Illumina assumption
>     does not generalize to length-variable or longer-read corpora.
>   - **c4** (SplitMix64 hash, CRAM-exact) is *worse than V3* on
>     every corpus (8-21% worse) — the hash-escalation path is
>     conclusively refuted across Illumina diversity.
>   - **HG002 PacBio HiFi (~14K reads, 18.5 kb mean): c0 (V3
>     baseline) wins, every richer-context candidate is *worse*.**
>     Sourced via raw GIAB `*.fastq.gz` → `samtools import` (public
>     PacBio BAMs strip QUAL; see `feedback_pacbio_hifi_qual_stripped`).
>     PacBio HiFi qualities cluster at Q60+ with a peak at Q93;
>     narrow distribution defeats richer context conditioning.
>     B/qual ~0.42 may be near-optimal for this data shape.
>
>   **Re-charter pending — design space changed.** No single static
>   bit-pack design wins across platforms: c3 dominates Illumina;
>   c0 dominates PacBio HiFi. Stage 2 must adopt an adaptive
>   mechanism (per-platform header flag, per-block adaptive sloc,
>   or two-pass encoder picking smaller of c3/c0) rather than lock
>   in a single bit budget. The 1.15× CRAM gate is chr22-v1.2.0-
>   specific; multi-platform robustness is a separate acceptance
>   criterion. Results:
>   `docs/benchmarks/2026-05-02-m94z-v4-multi-corpus.md` (cross-corpus
>   summary, 4 corpora), per-corpus docs
>   `2026-05-02-m94z-v4-{candidates, na12878_wes_chr22,
>   hg002_illumina_2x250_chr22, hg002_pacbio_hifi}.md`.
>
>   **Stage 2 (V4 CRAM 3.1 fqzcomp port) — see top-level status above
>   for shipped outcomes; this block kept as historical context.** All
>   15 tasks landed 2026-05-02 (Phase 0 sanity → Phase 1 RC primitives
>   → Phase 2 model → Phase 3 strategies/auto-tune → Phase 4 wire
>   format + Python ctypes → Phase 5 cross-corpus + docs). Phase 0
>   verdict (auto-tune does NOT save PacBio at the C level) was
>   **superseded by the end-to-end measurement** — V4 beats V3 by 4%
>   on PacBio HiFi via the auto-tuner's per-block strategy selection.
>
>
> **Phase 7 — M92 release prep (v1.2.0) follows the codec trio.** Active
> sibling workplans (separate CHANGELOG sections under
> `[Unreleased]`, do not interleave with M-series):
> * [`docs/verification-workplan.md`](docs/verification-workplan.md) — **V-series** (V1-V9 + P1-P4 perf follow-ups; mostly complete as of 2026-04-27).
> * [`docs/coverage-workplan.md`](docs/coverage-workplan.md) — **C-series** (C1-C8, coverage debt repayment; targets Python ≥92%, Java ≥88%, ObjC ≥85% line coverage).
>
> **#11 codec engineering — Channel 1 (mate_info) SHIPPED v1.7
> 2026-05-03.** 14-task plan executed; mate_info v2 codec (MATE_INLINE_V2,
> codec id 13) ships as default across Python/Java/ObjC. Key outcomes:
> * **C kernel:** `native/src/mate_info_v2.{c,h}` + rANS-O0 helpers
>   in `native/src/rans_o0.c`. Entry points
>   `ttio_mate_info_v2_encode` / `ttio_mate_info_v2_decode`.
> * **4-substream wire format:** MF / NS / NP / TS (34-byte container
>   header + auto-pick rANS-O0).
> * **Cross-language byte-exact gate:** 12/12 PASS (4 corpora ×
>   Python ↔ Java ↔ ObjC) in
>   `python/tests/integration/test_mate_info_v2_cross_language.py`.
> * **chr22 compression gate (T15):** 47.922 MB savings (hard gate 5 MB,
>   PASS). Test at
>   `python/tests/integration/test_mate_info_v2_compression_gate.py`.
> * **chrom_names sidecar:** `signal_channels/mate_info/chrom_names`
>   compound dataset (encounter-order, own first) preserves mate-only
>   chroms absent from `genomic_index/chromosome_names`.
> * **opt-out flag:** `WrittenGenomicRun.opt_disable_inline_mate_info_v2`
>   preserves v1 round-trip for v1.6 readers.
> * Last code commit: `b68e132`. Task 15 (docs/gate/changelog): this
>   commit.
>
> **#11 Channel 2 (REF_DIFF v2) SHIPPED v1.8 2026-05-03.** 14-task plan
> executed; REF_DIFF v2 codec (REF_DIFF_V2, codec id 14) ships as default
> sequences codec across Python/Java/ObjC. Key outcomes:
> * **C kernel:** `native/src/ref_diff_v2.{c,h}` + rANS-O0 helpers.
>   Entry points `ttio_ref_diff_v2_encode` / `ttio_ref_diff_v2_decode`.
> * **5-substream wire format:** FLAG / BS / IN / SC / ESC (38 + uri_len
>   byte container header + per-slice bodies).
> * **Cross-language byte-exact gate:** 3/3 PASS + 1 SKIP (hg002_pacbio
>   BAM has SEQ=`*`) in
>   `python/tests/integration/test_ref_diff_v2_cross_language.py`.
> * **chr22 compression gate (T15):** 4.314 MB savings (hard gate 2 MB,
>   PASS). Test at
>   `python/tests/integration/test_ref_diff_v2_compression_gate.py`.
> * **opt-out flag:** `WrittenGenomicRun.opt_disable_ref_diff_v2`
>   preserves v1 round-trip for v1.7 readers.
> * Last code commit: `f4f0c38`. Task 15 (docs/gate/changelog): this
>   commit.
>
> **#11 Channel 3 (NameTokenized v2) SHIPPED v1.9 2026-05-04.** 16-task
> plan executed; NAME_TOKENIZED v2 codec (NAME_TOKENIZED_V2, codec id 15)
> ships as default `read_names` codec across Python/Java/ObjC. Key
> outcomes:
> * **C kernel:** `native/src/name_tok_v2.{c,h}` (~1370 lines).
>   Entry points `ttio_name_tok_v2_encode` / `_decode`. Multi-substream
>   + DUP-pool (N=8) + PREFIX-MATCH + per-block reset (4096 reads).
> * **8-substream wire format:** FLAG / POOL_IDX / MATCH_K / COL_TYPES
>   / NUM_DELTA / DICT_CODE / DICT_LIT / VERB_LIT, each auto-picked
>   rANS-O0 vs raw.
> * **Phase 0 prototype** at `tools/perf/name_tok_v2_prototype/`
>   validated the algorithm + chr22 ≥ 3 MB savings gate before any C
>   work, per `feedback_phase_0_spec_proof`.
> * **Cross-language byte-exact gate:** 12/12 PASS (4 corpora × 3
>   languages) in
>   `python/tests/integration/test_name_tok_v2_cross_language.py`.
> * **chr22 compression gate (T15):** 67.76 MB savings vs pre-v1.9
>   default (gate floor 3 MB, PASS). Test at
>   `python/tests/integration/test_name_tok_v2_compression_gate.py`.
> * **Three-way writer dispatch:** v2 default, v1 via override, M82
>   compound via opt-out for v1.x reader compat.
> * **opt-out flag:** `WrittenGenomicRun.opt_disable_name_tokenized_v2`.
>
> **#11 channel work COMPLETE.** Channels 1+2+3 all shipped (v1.7
> mate_info v2 → v1.8 REF_DIFF v2 → v1.9 NAME_TOKENIZED v2).
>
> **#10 offsets-cumsum SHIPPED v1.10 2026-05-04.** Drops the redundant
> `offsets` column from `genomic_index/`, `spectrum_index/`, and
> `chromatogram_index/` — `offsets[i]` is mathematically derivable from
> `cumsum(lengths[0..i])`, so storing both was pure on-disk redundancy.
> Readers synthesize via uint64 accumulator (overflow-safe on deep WGS).
> **chr22 savings: 3.50 MB.** Opt-out flag
> `opt_keep_offsets_columns` (Python) / `optKeepOffsetsColumns`
> (Java/ObjC) preserves byte-equivalent backward compat with pre-v1.10
> readers.
>
> **#13 V5 multi-stream rANS** (3-6× speedup for qualities, wire-format
> break, needs Phase 0 spec-proof per `feedback_phase_0_spec_proof`)
> remains as a separate future cycle.

---

## Phase 1 — Core Data Model (v0.1.0-alpha, complete)

### Milestone 1 — Foundation (Protocols + Value Classes) ✓
### Milestone 2 — SignalArray + HDF5 Wrapper ✓
### Milestone 3 — Spectrum + Concrete Spectrum Classes ✓
### Milestone 4 — AcquisitionRun + SpectrumIndex ✓
### Milestone 5 — SpectralDataset + Identification + Quantification + Provenance ✓
### Milestone 6 — MSImage (Spatial Extension) ✓
### Milestone 7 — Protection / Encryption ✓
### Milestone 8 — Query API + Streaming ✓

> All eight milestones shipped with v0.1.0-alpha. 379 tests on
> `ubuntu-latest`. See `docs/version-history.md` for per-milestone
> detail.

---

## Phase 2 — Cross-Language Parity + Production Features (v0.2–v1.0, complete)

### M9–M78 (v0.2.0 through v1.0.0) ✓

72+ milestones delivered across three languages (Objective-C, Python,
Java). Key capabilities shipped:

- Four storage/transport providers (HDF5, Memory, SQLite, Zarr v3)
- Seven importers, seven exporters
- .tis streaming transport protocol with WebSocket server/client
- Per-AU encryption (AES-256-GCM, ML-KEM-1024 key wrap, ML-DSA-87 signatures)
- Spectral anonymisation (proteomics, metabolomics, NMR)
- LZ4 + Numpress compression codecs
- MS2 activation/isolation data model (M74)
- ISA-Tab and ISA-JSON support
- Cross-language 3×3 conformance matrix

> See `CHANGELOG.md` for per-milestone detail. v1.0.0 tagged
> 2026-04-23.

---

## Phase 3 — TTI-O Rebrand + Genomic Foundation (v0.11, in progress)

### M79 — Modality Abstraction + Genomic Enumerations ✓

- `Precision.UINT8`, `Compression` values 4–8 (RANS_ORDER0 through
  NAME_TOKENIZED), `AcquisitionMode.GENOMIC_WGS/WES`, transport
  `spectrum_class=5`, `@modality` attribute, `opt_genomic` feature
  flag. All three languages.

### M80 — TTI-O Rebrand (Clean Sweep) ✓

- Repository-wide rename: MPEG-O → TTI-O, MPGO → TTIO, `.mpgo` →
  `.tio`, `.mots` → `.tis`, transport magic `"MO"` → `"TI"`. No
  backward compatibility.

### M81 — Java Reverse-DNS Correction ✓

- `com.dtwthalion.ttio` → `global.thalion.ttio`. Maven groupId
  corrected to `global.thalion`.

### M82 — GenomicRun + AlignedRead + Signal Channel Layout ✓

Shipped as five sub-milestones (M82.1–M82.5):

- M82.1: Python reference — `AlignedRead`, `GenomicIndex`,
  `GenomicRun`, `WrittenGenomicRun`. `/study/genomic_runs/` group
  with signal channels + genomic index. `SpectralDataset` extended.
- M82.2: ObjC normative — `TTIOAlignedRead`, `TTIOGenomicRun`,
  `TTIOGenomicIndex`, `TTIOWrittenGenomicRun`.
- M82.3: Java — `AlignedRead` record, `GenomicRun`, `GenomicIndex`,
  `WrittenGenomicRun`. `Precision.UINT64` added.
- M82.4: Cross-language wire parity — Java VL_STRING-in-compound
  fix. Full 3×3 conformance matrix. `TtioWriteGenomicFixture` +
  `TtioVerify` genomic extensions.
- M82.5: Documentation — `docs/M82.md`, `ARCHITECTURE.md` genomic
  section, format-spec §10 codec gap correction.

> Python 854 tests, ObjC 1935 assertions, Java 402 tests. Genomic
> data stored with zlib compression; CRAM-derived codecs ship in
> M83–M86.

---

## Phase 4 — Genomic Compression Codecs (v0.11, in progress)

**Architectural premise:** All compression codecs are clean-room
implementations from published academic literature. No htslib source
code is consulted. htslib (via samtools) is used only as a runtime
bridge for BAM/CRAM file import and export. Sources: rANS from Duda
2014 (public domain); base packing from standard 2-bit nucleotide
encoding; quality quantisation from published Illumina binning tables;
name tokenisation from Bonfield, Bioinformatics 2022.

### M83 — rANS Entropy Codec + MPGO Remnant Cleanup (in progress)

- Clean-room rANS order-0 and order-1 encoder/decoder across all
  three languages.
- Self-contained wire format (header + frequency table + payload).
- Canonical test vectors with cross-language byte-exact conformance.
- `MPGOInstrumentConfig.h` renamed to `TTIOInstrumentConfig.h`;
  zero MPGO references remaining in ObjC source.

### M84 — Base-Packing + Quality Quantiser Codecs

- `ttio.codecs.base_pack`: 2-bit ACGT packing (4 bases per byte),
  N positions in side-channel. 4:1 raw compression before entropy
  coding.
- `ttio.codecs.quality`: Phred score quantisation with schemes
  `"illumina-8"` (8 bins), `"crumble-40"` (near-lossless),
  `"lossless"` (passthrough). Per-run configurable.
- All three languages. Cross-language byte-exact conformance.

### M85 — Quality Quantiser + Name Tokeniser Codecs

**Status: Phase A and Phase B both shipped (2026-04-26).**

The original M84 sketch in this WORKPLAN bundled "Base-Packing +
Quality Quantiser Codecs"; the M84 milestone shipped only
BASE_PACK on 2026-04-26, and quality_binned slipped to M85. M85
is now restructured into two phases: Phase A is the catch-up
quality codec, Phase B is the original name_tokenizer scope.

#### Phase A — quality_binned codec (SHIPPED 2026-04-26)

- [x] `ttio.codecs.quality` — fixed Illumina-8 / CRUMBLE-derived
      8-bin Phred quantisation table, 4-bit-packed bin indices,
      big-endian within byte. Lossy by construction;
      `decode(encode(x)) == bin_centre[bin_of[x]]`.
- [x] All three languages, cross-language byte-exact fixtures
      (`quality_{a,b,c,d}.bin`).
- [x] Phred 41+ saturates to bin 7 / centre 40 (documented).
- [x] Wire format: `[1 B version=0x00][1 B scheme_id=0x00][4 B
      BE orig_len][packed body of ceil(orig/2) bytes]`.
- [x] M79 codec id `7` slot now has a working encoder + decoder.
- [x] `docs/codecs/quality.md` codec spec.
- [x] `docs/format-spec.md` §10.4 quality-binned row flipped.

Commits: `9cfb08b` (Python), `5ad8952` (Java), `1449502` (ObjC),
plus M85 Phase A docs landing alongside.

Throughput on the M85 Phase A reference host (4 MiB random Phred
mod 41): Python 61 / 471 MB/s, Java 2001 / 425 MB/s, ObjC 3203 /
2196 MB/s.

#### Phase B — name_tokenizer codec (SHIPPED 2026-04-26)

- [x] `ttio.codecs.name_tokenizer` — lean two-token-type
      columnar codec (numeric digit-runs without leading zeros +
      string non-digit-runs absorbing leading-zero digit-runs),
      per-column type detection (columnar mode vs verbatim
      fallback), delta-encoded numeric columns,
      inline-dictionary-encoded string columns.
- [x] All three languages, cross-language byte-exact fixtures
      (`name_tok_{a,b,c,d}.bin`).
- [x] M79 codec id `8` slot now has a working encoder + decoder.
- [x] `docs/codecs/name_tokenizer.md` codec spec.
- [x] `docs/format-spec.md` §10.4 name-tokenized row flipped.
- Compression target was originally ≥ 20:1 (CRAM 3.1 / Bonfield
  2022 style); Phase B's lean implementation achieves ~3–7:1.
  Reaching 20:1 requires the full Bonfield-style encoder (eight
  token types — DIGIT0, MATCH, DUP, ALPHA-vs-CHAR distinction,
  per-token-type encoding variants — multi-thousand lines per
  language with substantial cross-language byte-exact
  conformance work). Tracked as a future optimisation milestone.

Commits: `cf665e7` (Python), `6a7e22a` (Java), `a66e750` (ObjC),
plus M85 Phase B docs landing alongside.

Throughput on the M85 Phase B reference host (100k Illumina-style
names, 2.18 MB raw): Python 5.4 / 16.5 MB/s, Java ~14 / ~63
MB/s, ObjC 37.0 / 310.7 MB/s.

#### Forward reference

A future M86 phase will wire QUALITY_BINNED (id 7) into the
`qualities` channel and NAME_TOKENIZED (id 8) into the
`read_names` channel of `signal_channels/`. The wiring
infrastructure (per-channel `signal_codec_overrides` dict,
`@compression` attribute, lazy-decode cache) is already in
place from M86 Phase A for byte channels; the `read_names`
channel will additionally need lifting from VL_STRING-in-
compound storage to a flat dataset that can carry the
`@compression` attribute.

#### Codec stack status (post-M85 Phase B)

All five M79 codec slots (4–8) now have working encoders +
decoders across Python, ObjC, and Java with cross-language
byte-exact conformance fixtures. The genomic codec library is
conceptually complete; remaining work is pipeline-wiring
(future M86 phases) and optional optimisation milestones for
higher compression ratios.

### M86 — Codec Integration into Signal Channel Pipeline

**Status: Phases A, B, C, D, E, and F all shipped (2026-04-26). The genomic codec pipeline-wiring is complete for ALL M82 channels.**

#### Phase A — byte channels (SHIPPED 2026-04-26)

- [x] Wire rANS order-0/1 and BASE_PACK into the genomic signal
      channel write/read paths for the **byte channels only** —
      `sequences` and `qualities`.
- [x] `@compression` attribute on coded channels (uint8 holding
      the M79 codec id; provider-agnostic data transform, not
      HDF5 filter; lives on the dataset, not the parent group).
      Note: original WORKPLAN sketched `@_ttio_codec`; final
      shipped name is `@compression` to match the format-spec
      convention.
- [x] `WrittenGenomicRun.signal_codec_overrides` (per-channel
      opt-in dict) on Python, ObjC, Java.
- [x] Validation rejects overrides for non-byte channels and
      non-codec values (sequences/qualities + RANS_ORDER0/1 +
      BASE_PACK only).
- [x] Codec-compressed datasets carry no HDF5 filter (no
      double-compression).
- [x] Lazy whole-channel decode + per-instance cache on the open
      `GenomicRun` object.
- [x] All three languages + cross-language `.tio` fixture
      conformance: 3 fixtures (one per codec) generated by
      Python; ObjC and Java each read all three byte-exact.
- [x] BASE_PACK on pure-ACGT sequences hits ~19% of uncompressed
      size (target was < 30%).

Commits: `31b0fa48` (Python), `0450913` (Java), `d2befc20` (ObjC),
plus M86 docs landing alongside.

#### Phase D — QUALITY_BINNED wiring (SHIPPED 2026-04-26)

- [x] Extend the M86 Phase A dispatch to accept
      `Compression.QUALITY_BINNED` (M79 codec id `7`) on the
      qualities byte channel of genomic runs.
- [x] Per-channel allowed-codec map: sequences accepts
      RANS_ORDER0/1 + BASE_PACK; qualities accepts those plus
      QUALITY_BINNED. Validation rejects `(sequences,
      QUALITY_BINNED)` with a clear lossy-quantisation error
      message (Binding Decision §108) — applying Phred-bin
      quantisation to ACGT bytes would silently destroy the
      sequence data.
- [x] Reuses M85 Phase A `Quality.encode/decode` codec and the
      M86 Phase A wiring infrastructure (per-channel
      `signal_codec_overrides` dict, `@compression` attribute,
      lazy-decode cache); zero new infrastructure.
- [x] All three languages, cross-language `.tio` fixture
      conformance: 1 new fixture
      (`m86_codec_quality_binned.tio`, 48 432 bytes) generated
      by Python; ObjC and Java both read it byte-exact (BASE_PACK
      on sequences + QUALITY_BINNED on qualities, bin-centre
      values for byte-exact lossy round-trip).
- [x] QUALITY_BINNED on a 100k-byte qualities channel hits ~50%
      of the HDF5-chunked-ZLIB baseline.

Commits: `425393a` (Python), `5d09294` (Java), `e7a85c6` (ObjC),
plus M86 Phase D docs landing alongside.

#### Phase B — integer channels (SHIPPED 2026-04-26)

- [x] Wire `RANS_ORDER0` and `RANS_ORDER1` (M79 codec ids 4 + 5)
      into the three integer channels of `signal_channels/`:
      `positions` (int64), `flags` (uint32),
      `mapping_qualities` (uint8).
- [x] Defines the int↔byte serialisation contract: integer
      arrays serialise to **little-endian** bytes per element
      before encoding; reader looks up original dtype by
      channel name (no on-disk dtype attribute, per Binding
      Decision §115).
- [x] Per-channel allowed-codec map gains positions/flags/
      mapping_qualities entries (rANS only). Validation
      rejects BASE_PACK / QUALITY_BINNED / NAME_TOKENIZED on
      integer channels (wrong-content codecs) per §117.
- [x] New per-instance `_decoded_int_channels` cache (separate
      from byte-channel and read-names caches per Binding
      Decision §116).
- [x] All three languages, cross-language `.tio` fixture
      conformance: 1 new fixture
      (`m86_codec_integer_channels.tio`, 60 720 bytes) with
      all three integer channels under rANS — Python, ObjC,
      Java all decode byte-exact.
- [x] Compression on clustered-positions int64 (10000-read /
      100-locus pattern): ~18% of raw LE bytes, byte-identical
      across all three languages (M83 rANS conformance).

**Important caveat (Binding Decision §119):** The current M82
read path for per-read integer fields uses `genomic_index/`,
not `signal_channels/`. Phase B compression is therefore
primarily a **write-side file-size optimisation**; it does not
currently affect read performance through `__getitem__` /
`alignedReadAt(int)`. The new `_int_channel_array(name)`
helper IS callable for direct access; future readers that
prefer `signal_channels/` (streaming, M89 transport) will
benefit transparently.

Commits: `798ec88` (Python), `46076b1` (Java), `6213711`
(ObjC), plus M86 Phase B docs landing alongside.

#### Phase C — VL_STRING / compound channels (PARTIAL: cigars SHIPPED 2026-04-26; mate_info DEFERRED)

**cigars (SHIPPED):**

- [x] Schema-lift cigars from M82 compound to flat 1-D uint8
      via `@compression` attribute, mirroring Phase E for
      read_names.
- [x] cigars accepts THREE codec choices:
      `Compression.RANS_ORDER0`, `Compression.RANS_ORDER1`,
      `Compression.NAME_TOKENIZED`.
- [x] rANS path uses length-prefix-concat serialisation
      (`varint(asciiLen) + asciiBytes` per CIGAR), then
      `Rans.encode(bytes, order)`. Reader reverses.
- [x] NAME_TOKENIZED path calls `NameTokenizer.encode(cigars)`
      directly (codec's own self-describing wire format).
- [x] Validation rejects `(cigars, BASE_PACK)` and
      `(cigars, QUALITY_BINNED)` with clear messages.
- [x] `_decoded_cigars: list[str]` lazy cache (separate from
      `_decoded_read_names` per Binding Decision §123).
- [x] All three languages, two cross-language fixtures (one
      per codec path) — both read byte-exact across all three
      implementations.
- [x] Codec selection guidance documented in `docs/codecs/`
      and `docs/format-spec.md` §10.8: rANS is the
      recommended default for real WGS data; NAME_TOKENIZED
      is the niche choice for known-uniform CIGARs.
- [x] Empirical compression on 1000-read mixed CIGARs (80%
      "100M" + 10% "99M1D" + 10% "50M50S"), byte-identical
      across all three languages via M83 + M85B conformance:
      RANS_ORDER1 = 1111 bytes (~17× smaller than M82
      compound baseline ~18-29 KB); NAME_TOKENIZED falls back
      to verbatim mode and produces 5307 bytes.

Commits: `61fb946` (Python), `133762e` (Java), `f5ff91f`
(ObjC), plus M86 Phase C docs landing alongside.

**mate_info (SHIPPED 2026-04-26 — see Phase F below).**

The mate_info portion was originally deferred from Phase C
because the 3-field compound (chrom VL_STRING + pos int64 +
tlen int32) needed either per-field schema decomposition or
per-compound-field codec dispatch. Phase F shipped the
former approach: subgroup with three per-field flat datasets,
each independently codec-compressible.

#### Phase F — mate_info per-field decomposition (SHIPPED 2026-04-26)

- [x] Schema-lift mate_info from M82 compound to subgroup
      `signal_channels/mate_info/` containing three child
      datasets (`chrom`, `pos`, `tlen`) when any per-field
      override is set.
- [x] Three per-field virtual channel names exposed via
      `signal_codec_overrides`: `mate_info_chrom`,
      `mate_info_pos`, `mate_info_tlen`. The bare key
      `mate_info` is rejected with a clear error pointing at
      the per-field names (Binding Decision §126 / Gotcha
      §143).
- [x] Per-field codec applicability: chrom takes the cigars-
      style codec set ({RANS_ORDER0/1, NAME_TOKENIZED}); pos
      and tlen take the integer-channel set ({RANS_ORDER0/1}).
- [x] Partial overrides allowed: any one per-field override
      triggers the subgroup; un-overridden fields use natural
      dtype with HDF5 ZLIB inside the subgroup (Binding
      Decision §127).
- [x] Read-side dispatch on HDF5 link type for the bare
      `mate_info` link (group = Phase F; dataset = M82) per
      Binding Decision §128. First M86 phase introducing
      this dispatch axis. Per-language link-type query API:
      Python h5py exception-based; ObjC `H5Oget_info_by_name`;
      Java provider-abstract via `openGroup` exception (or
      H5.H5Oget_info_by_name).
- [x] New `_decoded_mate_info` combined dict cache keyed by
      field name (separate from existing five caches per
      Binding Decision §129).
- [x] All three languages, one cross-language fixture
      (`m86_codec_mate_info_full.tio`, 60 757 bytes) — both
      ObjC and Java decode all three mate fields byte-exact
      against the Python input over 100 reads with realistic
      mate distributions (90% chr1 / 10% other; monotonic
      positions for paired mates; tlen clustered around
      350bp).
- [x] Reuses Phase C cigars helpers for the chrom rANS path
      (length-prefix-concat) and Phase B integer-channel
      helpers for pos/tlen (LE byte serialisation).

Commits: `20ca474` (Python), `b0b3926` (Java), `4d28629`
(ObjC), plus M86 Phase F docs landing alongside.

#### Codec stack pipeline-wiring status (post-M86 Phase F)

The genomic codec pipeline-wiring is now **complete for ALL
M82 channels**. Every channel under `signal_channels/` has at
least one accepted codec; every M79 codec slot (4–8) is wired
into its applicable channels with cross-language byte-exact
conformance. The M82-era genomic codec story is functionally
done.

#### Phase E — NAME_TOKENIZED wiring + read_names schema lift (SHIPPED 2026-04-26)

- [x] When `signal_codec_overrides["read_names"]` is set, the
      writer replaces the M82 compound `read_names` dataset
      with a flat 1-D uint8 dataset of the same name containing
      the codec output, and sets `@compression == 8` on it.
- [x] The reader dispatches on dataset shape (compound → M82
      path; 1-D uint8 → codec dispatch via the lazy-decode
      list cache).
- [x] Per-channel allowed-codec map gains `read_names:
      {NAME_TOKENIZED}`. NAME_TOKENIZED is rejected on
      `sequences` and `qualities` (would mis-tokenise binary
      byte streams) per Binding Decision §113.
- [x] New per-instance cache `_decoded_read_names: list[str]`
      separate from the byte-channel cache per Binding
      Decision §114.
- [x] All three languages, cross-language `.tio` fixture
      conformance: 1 new fixture (`m86_codec_name_tokenized.tio`,
      48 432 bytes) generated by Python; ObjC and Java both
      read it byte-exact.
- [x] All read-side call sites audited (only `__getitem__`
      directly touches `read_names`; region queries inherit
      via `self[i]`).
- [x] Compression on 1000 structured Illumina names: ~50%
      (Python/Java via H5 storage-size baseline) to ~19%
      (ObjC via file-size delta baseline; both methodologies
      undercount differently due to VL_STRING heap accounting).

Commits: `d12d135` (Python), `6624406` (Java), `bc6158c`
(ObjC), plus M86 Phase E docs landing alongside.

The schema-lift pattern is documented in `docs/format-spec.md`
§10.6 alongside the M86 Phase A `@compression` attribute scheme
in §10.5.

#### Codec stack pipeline-wiring status (post-M86 Phase E)

All five M79 codec slots (4–8) are now wired into the genomic
signal-channel pipeline. The genomic codec stack is conceptually
complete for the byte and string channels; remaining work
(integer channels — Phase B; cigars/mate_info — Phase C) is
deferred future scope.

#### Original wishlist preserved for future scope

- Default pipeline: positions → delta + rANS order-1; sequences →
  base_pack + rANS order-0; qualities → quality_quantise + rANS
  order-0; flags → rANS order-0; read_names → name_tokenizer;
  cigars → RLE + rANS order-0.
- Compression ratio benchmarking vs. uncompressed (target ≥ 5:1).

---

## Phase 5 — Import/Export + Interoperability

### M87 — SAM/BAM Importer (SHIPPED 2026-04-26)

- [x] `ttio.importers.bam`: reads BAM via `samtools view -h`
      subprocess (no htslib link). `BamReader.to_genomic_run()`
      with optional region filter passed verbatim to samtools.
- [x] SAM header parsing: @SQ → `reference_uri`, @RG (first wins
      per Binding Decision §133) → `sample_name` + `platform`,
      @PG → `ProvenanceRecord` list (samtools also injects @PG
      records on view-bS / view-h, so the M87 fixture's
      provenance_count is 3 not 1).
- [x] `ttio.importers.sam`: thin subclass; samtools auto-detects
      SAM vs BAM format from magic bytes.
- [x] Per-language `bam_dump` CLI emitting canonical JSON
      (sorted keys, 2-space indent, MD5 fingerprints for the
      sequences/qualities buffers) — 1341 bytes byte-identical
      across Python, ObjC, and Java.
- [x] Cross-language harness
      (`python/tests/integration/test_m87_cross_language.py`)
      runs all three CLIs on the canonical fixture and asserts
      byte-equality.
- [x] All three languages, 16 test cases per language
      (samtools-availability, full-BAM read, per-field
      verification, mate-info, header parsing, region filter,
      SAM input wrapper, samtools-missing error message,
      round-trip-through-writer BAM → `.tio` → GenomicRun →
      AlignedRead).
- [x] samtools-not-on-PATH detected at first call, NOT at import
      time per Binding Decision §135. Error message includes
      apt/brew/conda install hints.
- [x] RNEXT `=` expanded to RNAME per Binding Decision §131.
      1-based positions preserved per §132.

Commits: `9504d1a3` (Python), `126410f` (Java), `b0ec762` (ObjC),
plus M87 docs landing alongside.

Test deltas: Python 898 → 915 (+17 incl. JSON shape check), ObjC
2477 → 2532 (+55 assertions across 16 test methods), Java 512 →
529 (+17 incl. JSON shape check).

### M88 — CRAM Importer + BAM/CRAM Exporter (SHIPPED 2026-04-26)

- `ttio.importers.cram` (and `TTIOCramReader` / Java `CramReader`):
  reads CRAM via `samtools view -h --reference <fasta>`. Subclasses
  the M87 BAM reader; reuses the SAM-text parsing path. Reference
  FASTA enforced at construction (Python `TypeError` / ObjC
  `NS_UNAVAILABLE` / Java null-rejection). Binding Decision §139.
- `ttio.exporters.bam` (and `TTIOBamWriter` / Java `BamWriter`):
  composes SAM text from `WrittenGenomicRun`, pipes through
  `samtools view -b -o <out.bam>`. All 11 SAM columns always emitted
  with sentinel fills; RNEXT collapses to `=` when chrom matches a
  mapped read (§136); negative `mate_position` maps to 0 on emit
  (§138); QUAL ASCII Phred+33 verbatim.
- `ttio.exporters.cram` (and `TTIOCramWriter` / Java `CramWriter`):
  two-stage pipeline `samtools view -C --reference <fa> | samtools
  sort -O cram --reference <fa> -o <out.cram>` because CRAM slices
  require sorted input. Subprocess chaining: ObjC `NSTask`+`NSPipe`,
  Java `ProcessBuilder` + transferTo pump thread, Python `Popen`
  with `stdin=prev.stdout`.
- Lossless round-trips verified per language: BAM↔BAM (per-field),
  CRAM↔CRAM (sequence buffer + names), BAM↔CRAM, CRAM↔BAM.
- Cross-language harness `test_m88_cross_language.py` re-uses M87
  `bam_dump` CLIs against the new M88 BAM fixture; byte-identical
  canonical JSON across Python / ObjC / Java. CRAM cross-language
  parity verified implicitly via shared canonical CRAM fixture in
  per-language unit suites.
- Test counts: Python 915 → 932 (+17), ObjC 2532 → 2595 (+63 across
  14 test methods), Java 529 → 543 (+14).

### M88.1 — `bam_dump --reference` flag for CRAM cross-language conformance (SHIPPED 2026-04-26)

- Closes the implicit-parity gap from M88. Extends the existing
  M87 `bam_dump` CLI in each language (Python `bam_dump`, ObjC
  `TtioBamDump`, Java `BamDump`) with an optional
  `--reference <fasta>` flag and case-insensitive `.cram`
  extension dispatch to the M88 `CramReader` / `TTIOCramReader`
  / Java `CramReader`.
- User chose Option A (extend existing CLI) over Option B
  (parallel `cram_dump` per language) — single CLI surface,
  single harness file, mirrors the
  `CramReader extends BamReader` inheritance.
- For BAM/SAM paths, `--reference` is accepted but unused
  (defensive). For `.cram` paths without `--reference`, exits 2
  with a stderr error message.
- Cross-language harness `test_m88_cross_language.py` now has
  six tests (3 BAM + 3 CRAM); CRAM canonical JSON for the M88
  fixture is **914 bytes** with md5
  `2be5c5bccc95635240aa60337406cb35` byte-identical across
  Python / ObjC / Java.
- Test counts: Python 1045 → 1047 (+2), ObjC 2595 → 2597, Java
  543 → 543 (no delta; coverage is in the cross-language
  harness).

---

## Phase 6 — Framework Integration

### M89 — Transport Layer Extension ✓ (shipped 2026-04-27)

- `.tis` GenomicRead AU payload: chromosome + position + mapq + flags
  prefix (replaces zeroed spectral fields from M79). `spectrum_class
  == 5` discriminator.
- `TransportWriter.write_genomic_run()` / `TransportReader`
  materialisation in all three languages.
- `AUFilter` extended with chromosome + position_range predicates.
- Multiplexed streams: MS + genomic runs interleaved in one `.tis`.
- Per-AU encryption on genomic AUs verified end-to-end.
- 3×3 cross-language transport matrix green.

### M90 — Encryption, Signatures, and Anonymisation for Genomic Data ✓ (shipped 2026-04-28)

Shipped as 15 sub-mileposts (M90.1–M90.15):

- M90.1–M90.6 — Per-AU AES-256-GCM on genomic signal channels with
  AAD = `dataset_id || au_sequence || channel_name`. ML-DSA-87
  signatures on genomic datasets.
- M90.7 — Java VL_STRING attribute writer/reader (
  `H5Awrite_VLStrings` / `H5Aread_VLStrings`). ObjC attribute reader
  follow-up landed in `a3495d4`.
- M90.8/9/10 — AU compound-field round-trip; UINT8 wire compression
  via M86 codecs; M90.10 cross-language parity in Java + ObjC.
- M90.11 — Encrypted genomic AU headers with per-region key map
  (reserved `_headers` key).
- M90.12 — UINT8-aware MPAD format ("MPA1" magic + per-entry dtype
  byte). Fixes the float64 cast bug.
- M90.13 — Region masking by SAM overlap.
- M90.14 — Seeded-RNG random quality scores (anonymiser).
- M90.15 — Sign chromosomes VL compound (Python +
  cross-language follow-up).
- Genomic anonymiser strips read names, randomises quality, masks
  HLA-style regions. Java + ObjC final parity in
  `f1728dc` / `cb728f7`.

### M91 — Multi-Omics Integration Test ✓ (shipped 2026-04-28)

- Single `.tio` carrying WGS + proteomics MS + NMR metabolomics with
  shared provenance and a unified encryption envelope.
- Cross-modality query (`runs_for_sample("sample://NA12878")` returns
  all three modalities).
- `.tis` transport multiplexing verified end-to-end.
- All three languages. Python ref impl in `9038f76`.

---

## Phase 8 — Post-M91 Abstraction Polish ✓ (shipped 2026-04-28)

OO design pass on the modality abstraction surface, driven by the
M91 cross-modality findings. Closes the gap where MS and genomic
runs each had their own accessor surface despite both being
indexable, streamable, and provenanceable.

### Phase 1 — Run protocol + modality-agnostic helpers ✓

- `Run` protocol (Python `runtime_checkable Protocol`, ObjC
  `@protocol TTIORun`, Java `interface Run`) — unified surface for
  `name`, `acquisition_mode`, `__len__` / `count` / `numberOfRuns`,
  `__getitem__` / `get` / `objectAtIndex:`, `provenance_chain` /
  `provenanceChain`. Both `AcquisitionRun` and `GenomicRun` conform.
- `dataset.runs_for_sample(uri)` / `runsForSample:` and
  `runs_of_modality(cls)` / `runsOfModality:` modality-agnostic
  accessors in all three languages.
- `GenomicRun.provenance_chain()` / `provenanceChain` exposed (closes
  the M91 read-side gap where genomic runs had no provenance API).

### Phase 2 — Mixed-dict write + per-run provenance dual-write ✓

- `SpectralDataset.write_minimal(runs={...})` accepts a mixed dict of
  MS + genomic runs and dispatches by isinstance / class. Same
  surface in Java (mixed `Map<String, Object>` overload) and ObjC
  (`mixedRuns:` parameter).
- Per-run provenance now writes the canonical compound dataset
  `<run>/provenance/steps` on the HDF5 fast path in all three
  languages, plus the legacy `@provenance_json` mirror for non-HDF5
  providers and pre-Phase-2 readers. Reader prefers compound, falls
  back to JSON.
- Anonymiser (ObjC) refactored onto the unified
  `writeMinimalToPath:` path; −270 / +150 lines.
- M51 cross-language byte-parity harness extended with an
  `ms_per_run_provenance` section so the Python / Java / ObjC
  dumpers' per-run output is byte-identical (`6200b4f`).

Commits: `145485c`, `772eb00`, `6992ae9`, `7a2ffef`, `54ef6f1`,
`6ceba4a`, `6200b4f`. All three test suites green
(Python 1324 / Java 755 / ObjC 3070).

---

## Phase 9 — Codec parity (v1.2.0 dependency for M92)

Three new codecs across three languages, all cross-language byte-exact
via canonical conformance fixtures. Closes the chr22 compression gap
identified in M92's smoke benchmark (TTI-O at 2.5× CRAM 3.1, target
1.15×).

### M93 — REF_DIFF reference-based sequence-diff codec ✓ (Python landed 2026-04-28; ObjC + Java in progress)

- Codec id `9`. Replaces `BASE_PACK` as the default for
  `signal_channels/sequences` when a reference is available.
- Context-aware per-channel codec — receives `(positions, cigars,
  reference_resolver)` from sibling channels at write/read time.
- Slice-based wire format (10 K reads/slice, CRAM-aligned) for
  random-access decode.
- Embedded reference at `/study/references/<reference_uri>/` with
  auto-deduplication across runs sharing a URI.
- Format-version bumps `1.4 → 1.5` only when REF_DIFF is actually
  used (M82-only writes stay at `1.4` for byte-parity).

### M94 — FQZCOMP_NX16 lossless quality codec ✓ (Python + ObjC + Java landed 2026-04-29)

- Codec id `10`. Replaces `RANS_ORDER0` as the default for
  `signal_channels/qualities` when on the v1.5 path.
- fqzcomp-Nx16 (Bonfield 2022 / CRAM 3.1 default) — interleaved 4-way
  rANS for SIMD parallelism, context model on `(prev_q[0..2],
  position_bucket, revcomp_flag, length_bucket)` hashed via
  SplitMix64 to 4096 contexts. Adaptive `+16` LR with halve-with-
  floor-1 renormalisation at 4096 max-count boundary.
- Wire-format header **54 + L bytes** (field-by-field sum); body
  has a 16-byte substream-length prefix before round-robin bytes.
- Auto-default on `qualities` gated on **v1.5 candidacy** to
  preserve M82 byte-parity (binding decision §80h).
- Python implementation links to a Cython C extension; ObjC + Java
  implement natively. **8 canonical fixtures byte-exact across all
  three languages.**
- Test counts post-M94: Python 313 / ObjC 3204 / Java 813 (all
  baselines unchanged + ~146 new M94 tests).
- Spec: §3 M94. Codec spec: `docs/codecs/fqzcomp_nx16.md`.

### M95 — DELTA_RANS_ORDER0 + integer-channel auto-defaults ✓ (Python + ObjC + Java landed 2026-04-30)

- Codec id `11`. Delta + zigzag + unsigned LEB128 varint + rANS
  order-0 wrapper for sorted-ascending integer channels.
- Wire format: 8-byte header (`DRA0` magic, version 1,
  `element_size` uint8, 2 reserved bytes) + rANS order-0 body.
- Auto-default integer channel compression on v1.5 genomic runs:
  `positions → DELTA_RANS_ORDER0`, `flags / mapping_qualities /
  template_lengths / mate_info_pos / mate_info_tlen → RANS_ORDER0`,
  `mate_info_chrom → NAME_TOKENIZED`.
- 4 canonical cross-language fixtures (`delta_rans_{a,b,c,d}.bin`)
  byte-exact across Python / ObjC / Java.
- Python reference: `python/src/ttio/codecs/delta_rans.py`.
  ObjC: `objc/Source/Codecs/TTIODeltaRans.{h,m}`.
  Java: `java/src/main/java/global/thalion/ttio/codecs/DeltaRans.java`.
- Codec spec: `docs/codecs/delta_rans.md`.

### M94.X — FQZCOMP_NX16 variable-total rANS (ABANDONED 2026-04-29)

Path 2 (variable-total rANS, eliminating the per-symbol M-normalisation
step) was approved 2026-04-29 but ultimately abandoned in favour of
**M94.Z** (CRAM-mimic 16-bit renormalisation with mathematically
guaranteed byte-pairing). The variable-total approach required a
wire-format break across three languages plus full fixture
re-canonicalisation; the CRAM-mimic codec achieves the same
algorithmic speedup with simpler invariants and ships as the
production v1.2.0 quality codec while M94 v1 (codec id 10) is
retained for backward compatibility.

### M94.Z — FQZCOMP_NX16_Z CRAM-mimic quality codec ✓ (Python + ObjC + Java landed 2026-04-29)

- Codec id `12` (registered alongside the v1 `FQZCOMP_NX16 = 10`,
  which is retained for backward-compat fixture readability).
- **Replaces M94.X as the v1.2.0 quality codec.** Mirrors the CRAM
  3.1 fqzcomp encoder's 16-bit renormalisation strategy:
  per-symbol normalisation is eliminated by maintaining state
  invariants across adaptive count updates, with rounding handled
  via a deterministic 16-bit lower-bound renorm guaranteed to
  produce byte-paired encoder/decoder output by construction
  (no insertion-sort tie-break, no fixed-`M` boundary).
- **Cross-language perf** (100K reads × 100bp Illumina synthetic,
  end-to-end API including wire-format pack/unpack):
  ObjC 52 MB/s encode / 30 MB/s decode; Python (Cython) 50 / 17
  MB/s (kernel alone ~95 / 50, wrapper adds ~50% encode / ~66%
  decode overhead from zlib freq tables + struct.pack); Java 33 / 14
  MB/s.
- 7 canonical `m94z_{a,b,c,d,f,g,h}.bin` fixtures byte-exact across
  Python / ObjC / Java.
- Wired into the M86 pipeline (`Compression.FQZCOMP_NX16_Z = 12`
  on `signal_channels/qualities`); benchmark adapter
  (`tools/benchmarks/formats.py`) updated to use the new codec id
  for quality compression.
- Codec spec: `docs/codecs/fqzcomp_nx16_z.md`.

### M93/M94 — Cython acceleration ✓ (Python only, 2026-04-29)

Pure-Python codec hot paths replaced with thin Cython kernels under
`python/src/ttio/codecs/_<codec>/_<codec>.pyx`, building on the
M94.Z infrastructure pattern. The `.pyx` source is committed; the
`.c` transpilation output and compiled `.so` are gitignored and
regenerated by `python setup.py build_ext`.

| Codec | Pack speedup | Unpack speedup | Notes |
|---|---|---|---|
| REF_DIFF (M93) | **44×** | **35×** | `_ref_diff.pyx` |
| NAME_TOKENIZED (M86 Phase E) | 11% chr22 drop | 8% chr22 drop | `_name_tokenizer.pyx` |
| FQZCOMP_NX16 (M94 v1) | (built-in) | (built-in) | `_fqzcomp_nx16.pyx` |
| FQZCOMP_NX16_Z (M94.Z) | 50 MB/s (e2e) | 17 MB/s (e2e) | `_fqzcomp_nx16_z.pyx`; kernel ~95/50, wrapper overhead ~50%/66% |

ObjC and Java carry their native fast paths separately (no Cython
analogue). Pure-Python `_*_py` fallbacks remain available for
environments that haven't built the C extensions; gated on
`_HAVE_C_EXTENSION`.

### Pipeline defects — `_mate_info_is_subgroup` caching + codec-import hoist ✓ (2026-04-29)

Two cumulative fixes in `python/src/ttio/genomic_run.py` discovered
during M94.Z chr22 profiling:

- **`_mate_info_is_subgroup()` caching.** Was burning ~50% of decode
  wall on redundant HDF5 link probes (5.3M calls all returning the
  same answer). Now memoised on `_mate_info_subgroup_cached` per
  run instance.
- **Codec import hoist.** Per-read `importlib._handle_fromlist`
  was eating ~6% of decode wall on attribute lookups; codec module
  imports (`Compression`, `Precision`, `_hdf5_io`, `codecs.rans`,
  `codecs.name_tokenizer`) hoisted to module load.

### Cumulative chr22 perf (post-M94.Z + Cython + pipeline fixes)

| Metric | Pre-session | Post-session | Speedup | vs CRAM 3.1 |
|---|---|---|---|---|
| Encode | 18 min | 27.91 s | 38.7× | 9.2× off (was 355×) |
| Decode | 24.6 min | 21.76 s | 67.8× | 13.4× off (was 1162×) |

Codec compute is now ~4% of TTI-O wall time; the remaining gap to
CRAM 3.1 is HDF5 framework overhead + multi-omics infrastructure
(metadata round-trip, provenance compounds, modality dispatch),
which is expected per user direction — TTI-O is a multi-modal
container and not aiming for byte-tight parity with a
genomics-only single-modality format.

#### v1.2.0 acceptance gate status

**Acceptance gate (carried forward from M92).** TTI-O lossless within
1.15× of CRAM 3.1 on the chr22 lean fixture
(`data/genomic/na12878/na12878.chr22.lean.bam`, 1.78M reads,
aux-stripped). Status: **9.2× / 13.4× off**, well off the 1.15×
target but down from the 355× / 1162× originally observed pre-codec
work. Codec compute is no longer the bottleneck; remaining gap is
HDF5 + multi-omics framework overhead, which user has signed off as
acceptable for v1.2.0 scope. Verification test:
`python/tests/integration/test_m93_compression_gate.py`.

---

## Phase 10 — Performance + Scale (v1.3+)

### M96 — FQZCOMP_NX16 slice-level parallelism (Path 4b, v1.3)

Independent of M94.X (algorithmic single-thread speedup). M96
adds CRAM-style slice-level parallelism on top of M94.X to scale
FQZCOMP_NX16 encode/decode to N cores for full-WGS workloads.

#### Motivation

M94.X (variable-total rANS) lifts single-thread encode from
~0.2 MB/s to ~30+ MB/s on each language. That's sufficient for
chr22-scale acceptance gates and most real-world reference-genome-
scale workloads. For full-WGS scale (50-500 GB BAM inputs), where
per-format encode wall time matters at the seconds-vs-minutes
level, slice-level parallelism multiplies single-thread throughput
by core count.

#### Design

- **Slicing strategy.** Mirror M93 REF_DIFF: 10 K reads/slice,
  CRAM-aligned. Per-slice independence means each slice carries
  its own freq-table state (no cross-slice adaptive context).
- **Wire format.** Add a slice index to the FQZCOMP_NX16 header
  identical in shape to REF_DIFF's slice index (per-slice byte
  offset + length + first/last position). Body becomes a
  concatenation of independent per-slice encoded blobs.
- **Single-slice (`num_slices == 1`)** preserves M94.X output
  byte-exact; users who don't enable parallelism see no change.
- **Parallelism dispatch.** Each language uses its idiomatic
  thread-pool: Python `concurrent.futures.ThreadPoolExecutor` (the
  Cython extension releases GIL during encode); ObjC `NSOperationQueue`
  with `maxConcurrentOperationCount = NSProcessInfo.processorCount`;
  Java `java.util.concurrent.ForkJoinPool.commonPool()`.

#### Trade-off vs M94.X scope

M96 requires a SECOND wire-format change to FQZCOMP_NX16 after
M94.X's first one. We deliberately separate the two milestones:

- **M94.X** ships in v1.2.0 (release blocker) and resolves the
  catastrophic single-thread speed problem.
- **M96** ships in v1.3 once user demand for full-WGS scale
  emerges. Requires re-canonicalisation of fixtures a second time;
  v1.2 users keep working files via single-slice fallback.

Combining both into one milestone would conflate two distinct
problems (algorithm vs orchestration) and double the cross-
language byte-exact verification matrix. Keeping them sequential
is cheaper.

#### Acceptance gates

- 8 canonical fixtures byte-exact across Python / ObjC / Java
  with `num_slices == N` for N ∈ {1, 4, 16}.
- chr22 mapped-only encode wall time scales sub-linearly with
  core count (4 cores → ≥3× speedup over M94.X single-thread).
- Single-slice output is byte-identical to M94.X output (no
  regression for users not opting in to parallelism).
- Decoder is parallel-or-serial agnostic (decoder picks a
  parallel strategy based on its own runtime, not the encoder's).

#### Scope notes

- **Decoder also parallelisable** but each slice's decode is
  small enough that thread-startup overhead dominates for chr22-
  scale. v1.3 ships parallel encode; parallel decode is a v1.3.x
  follow-up if profiling shows it matters.
- **Same parallelism model could extend to M93 REF_DIFF** which
  already has slice structure — M93's encoder is currently
  single-thread despite the slice layout. Bundle into M96 as a
  unified "parallel slice encoders for genomic codecs" milestone.

---

## Phase 7 — Release

### M92 — Benchmarking, Documentation, and v1.2.0 Tag

- Compression benchmarking report: TTI-O genomic vs. BAM, CRAM 3.1,
  and MPEG-G (Genie) on NA12878 WGS (downsampled), ERR194147 WES,
  and a synthetic mixed-chromosome dataset.
- Documentation refresh: README, ARCHITECTURE, migration guide.
- v1.2.0 CHANGELOG entry. (Tag chosen to follow the existing
  v1.0.0 / v1.1.0 / v1.1.1 line; the prior workplan said "v0.11.0"
  but that tag already exists from the pre-stable line and would
  go backwards.)
- Tag gated on user sign-off.

**Acceptance Criteria**

- [ ] Compression within 15% of CRAM 3.1 on all benchmarks (lossless).
- [ ] All three language test suites pass.
- [ ] Cross-language 3×3 conformance matrix green for .tio + .tis
      with genomic data.
- [ ] Multi-omics integration test (M91) green.

---

## Binding Decisions (Genomic Integration)

| # | Decision | Rationale |
|---|---|---|
| 60 | `modality` is a UTF-8 string attribute, not an enum integer. | Extensible for future modalities without enum-value coordination. |
| 61 | Codecs are provider-agnostic data transforms with `@_ttio_codec` attribute, not HDF5 filter plugins. | Works identically on Memory, SQLite, Zarr, and HDF5. HDF5 filter registration deferred to v1.1+. |
| 62 | Base-packing uses 2-bit big-endian within each byte (MSB-first). N positions in separate side-channel. | Bit-compatible with CRAM convention; simplifies cross-validation. |
| 63 | Quality quantisation scheme is a per-run attribute (`@quality_quantisation`). | Different runs in the same file may use different quantisation. |
| 64 | `samtools` is a runtime dependency only for BAM/CRAM import/export. TTI-O's compression pipeline has zero htslib dependency. | htslib's value is as a BAM/CRAM format parser, not as a compression library. |
| 65 | GenomicRead AU payload uses a genomic-specific prefix (chromosome + position + mapq + flags). | Genomic filter keys are fundamentally different from spectral filter keys. |
| 66 | All compression codecs are clean-room implementations from published literature. No htslib source consulted. | Independence and credibility: TTI-O's pipeline is an original implementation of public-domain algorithms. |
| 67 | No backward compatibility with `.mpgo` / `mpeg_o_*` attributes. Clean break. | No external users. Dual-read logic adds complexity with zero benefit. |
| 68 | Transport magic: `"TI"` (Thalion Initiative). | Two-byte mnemonic, no known collision. |
| 69 | `"MPAD"` debug dump magic not renamed. | Internal diagnostic format, not public wire spec. |
| 70 | Genomic runs under `/study/genomic_runs/`, separate from ms_runs. | Clean dispatch; ms_runs iterators never see genomic groups. |
| 71 | M82 stores one ASCII byte per base. Base-packing deferred to codec milestone. | Separation of concerns: validate data model first. |
| 72 | Flags stored as UINT32 (not UINT16). | Future-proofing for extended flag bits. |
| 73 | GenomicIndex.chromosomes is `list[str]`, stored as compound VL_BYTES. | Variable-length strings; established pattern. |
| 74 | Data model classes ship in all three languages in the same milestone. | Language parity is a core principle. |
| 75 | rANS state: 64-bit, L=2^23, b=2^8. | Standard Duda 2014 parameterisation. |
| 76 | Frequency table normalisation total M=4096. | Power-of-two for fast modulo. |
| 77 | Wire format: big-endian, self-contained header. | Embeddable in HDF5 datasets or .tis packets without external metadata. |
| 78 | Frequency rounding: descending count, stable by symbol value. | Deterministic cross-language byte-exact conformance. |
| 79 | MPGO remnant cleanup mandatory in M83. | Zero tolerance for stale pre-rebrand identifiers. |
| 80 | M94.Z V2 wire format (version byte = 2) carries the libttio_rans byte layout `[4×states LE][4×lane_sizes LE][per-lane data]`; V1 (version byte = 1) remains canonical. | New V2 streams encode at native rANS speed without breaking V1 readers. |
| 81 | `prefer_native` is opt-in for V2 writes (parameter or `TTIO_M94Z_USE_NATIVE` env var); V1 is the default. | Existing call sites are unaffected; new builds opt in explicitly. |
| 82 | ObjC writer chain takes `id<TTIOStorageGroup>` (Python/Java M44 parity, Task 31). NMR runs and Image-subclass datasets stay HDF5-only. | The HDF5-direct features (H5DSset_scale, native 3D cubes) have no protocol equivalents — same scope as Python and Java. |
| 83 | Provider compression-unsupported behaviour is "accept argument and silently ignore" per `<TTIOStorageProvider>` protocol. | Memory / SQLite / Zarr accept `TTIOCompressionZlib` from upstream callers without erroring; HDF5 honors it. |
| 84 | `@read_role` (M97) is a free UTF-8 run attribute with the recognised vocabulary `hifi`, `ont_ul`, `hic_r1`, `hic_r2`, `parental_maternal`, `parental_paternal`, `illumina_wgs`; writers store other strings unchecked, readers surface absence as `None`/`null`/`nil`. | Same rationale as §60: new roles must not need a coordinated 3-language enum release. |
| 85 | REF_DIFF_V2 `slice_bytes` (M97) is writer policy only: a slice closes before the read that would push it past the budget, `reads_per_slice` still caps the count, every slice keeps at least one read, and 0 reproduces the fixed-count output byte for byte. | Phase 0 proof 2026-08-24: the decoder walks the slice index with no boundary assumption; byte-budget blobs decode byte-exact with the unmodified decoder. |
| 86 | A QUALITY_BINNED `qualities` override is rejected (M97) when `@platform` matches `hifi` / `pacbio` / `nanopore` as case-insensitive substrings or `ont` as a whole token. | The fixed Illumina-8 bin table mistunes long-read quality distributions; the token rule keeps IONTORRENT allowed. |
| 87 | `ttio_ref_diff_v2_max_encoded_size2` bounds the budgeted slice count by `4 * total_bases / slice_bytes + 2`, capped at `n_reads`. | A slice under half the budget is always followed by one at least that big, so of any two consecutive slices one holds at least `slice_bytes / 2` bases. |
| 88 | An assembly graph (M98) stores at `/study/assembly_graphs/<name>/` as four compound tables (`segments/records`, `links`, `paths`, `extras`) plus a concatenated `segments/sequences` uint8 channel and a `line_index` compound (`line_type` 0=S 1=L 2=P 3=extra, `row`) that replays the original line order on emission. | The GFA file is the interchange artefact; storing the file's line order makes re-emission byte-exact without a verbatim text blob, and the tables stay queryable. |
| 89 | Empty assembly-graph tables are absent (M98); readers treat a missing table as empty, and the `@final_newline` attribute is the structural marker of a graph group. | 0-row non-extendable compounds do not round-trip on the memory provider; the attribute marker keeps open-validation independent of table presence. |
| 90 | The assembly `segments/sequences` channel (M98) stores through BASE_PACK when every byte is in `ACGTNacgtn`, RANS_ORDER1 otherwise, uncompressed when empty, with the codec id in the dataset's `@compression` attribute. | Reuses the §10.5 signal-channel convention and the existing codec stack; the byte-level rule is cheap to evaluate and identical in 3 SDKs. |
| 91 | GFA parsing (M98) is structural, not semantic: `S` needs 3 fields, `L` 6, `P` 4; every other line goes verbatim into `extras`; `tags` is the tab-joined verbatim remainder; `*` stores no bytes and sets `seq_missing`; `gfa_version` comes from the first extras `H` line's `VN:Z:` value. | hifiasm `A` lines, comments and future GFA record types survive unparsed; nothing the parser does not understand can break the byte-exact round-trip. |
| 92 | Per-AU encryption over an assembly graph (M98) is one AES-GCM operation per segment record, offsets/lengths from `segments/records`, dataset_id continuing after the genomic runs; the channel is codec-decoded before slicing and decrypt-in-place restores the raw uint8 channel with no `@compression`. Signatures cover `segments/sequences` only. | Segment records are the graph's access units the way reads are a run's; writing raw bytes back keeps re-emission byte-exact while treating the codec as an at-rest optimisation, the M90.12 MS precedent. |
| 93 | Per-AU encryption of a `blocks_v1` genomic run (M99) walks block by block: per-block codec decode from the block index, one AU per read with global AU numbering and global plaintext offsets, extendable segments tables appended per block. Decrypt-in-place re-encodes each block through the block writer under the sticky qualities discipline (block 0 auto-tunes, the winner read back from the encoded stream pins the rest). The block index is rewritten only when a re-encoded blob does not land on its recorded ranges (M99.1); byte-identity holds whenever the persisted writer policy attrs are honoured. | The block policy already bounds a block's decoded size, so streaming falls out of the layout; the sticky pin is recoverable from the first re-encoded blob, so no restore metadata is stored; and the index rewrite turns a policy mismatch into a consistency event instead of a refusal. |
| 94 | The blocks_v1 stream writers persist non-default writer policy on the run group (`@ref_diff_slice_bytes`, `@opt_disable_qualities_v5`) and the reference set of a REF_DIFF_V2 run (`@reference_md5s`, each chromosome mapped to the hex reference-set digest); restore re-encodes under those attrs, resolving the reference from `/study/references/` or `REF_PATH` through the reference resolver (M99.1, replacing the M99 encrypt-time re-encode and byte-compare gate). When a re-encoded blob still differs, as for an older file without the attrs, restore rewrites the block index instead of refusing. | Persisting the policy makes restore reproducible by construction instead of verified per encrypt; the reference-set digest is the one the blob headers carry and the resolver verifies, so a REF_PATH restore needs no embedded copy; the fallback keeps pre-M99.1 files restorable, readable and consistent at the cost of byte-identity. |
| 95 | The encrypted transport stream carries `blocks_v1` genomic runs via the transport-spec v0.12 sidecar packets — one GenomicRunSidecar (0x1C) per run with the run scalars, restore attrs and name tables, one BlockSidecar (0x1D) per block with the index row and the verbatim plaintext blob slices — announced by the required `transport_blocks_v1` StreamHeader feature token (M99.1, replacing the M99 sender refusal). Pre-v0.12 receivers skip the sidecar packets and rebuild an unrestorable legacy-shaped container, so senders must not target them; the stream does not carry embedded reference bytes. | The AU stream already carries everything per-read; the sidecars carry exactly the run- and block-scoped remainder, so a received container decrypts in place byte-identically, and the feature token is the wire-level detection handle the skip-unknown-packets rule cannot provide. |

---

## Deferred to v1.1+

| Item | Description |
|---|---|
| HDF5 filter registration | Register rANS/base_pack/quality/name_tokenizer as official HDF5 filter IDs with The HDF Group. |
| Population-level deduplication | Cross-sample content-defined chunking and dedup at the genomic region level. |
| Long-read support | PacBio HiFi / Oxford Nanopore: longer sequences, higher error rates, different codec tuning. **Profile piece pulled forward to M97** (2026-08-24): `@read_role`, the QUALITY_BINNED platform guard, and the REF_DIFF_V2 `slice_bytes` budget; codec tuning itself stays deferred (the M94.Z auto-tuner already picks V4 on HiFi and V5 on ONT per block). |
| VCF/gVCF importer | Import variant call files as VariantAnnotation records. |
| htsget-style REST API | HTTP range-request protocol for serving genomic regions from .tio files. |
| Methylation / epigenomics | Extended base alphabet (5mC, 5hmC) and modified-base probability arrays. |
| CRAM native interop | Direct htslib linking for high-throughput BAM/CRAM conversion (bypass samtools subprocess). |
| Transcriptomics modality | RNA-seq: splice-junction awareness, expression quantification, gene-level aggregation. |
| ParquetProvider | Apache Parquet as a fifth storage provider. |
| FIPS compliance mode | FIPS 140-3 validated crypto backend. |
| DBMS transport | Database-backed streaming for enterprise deployments. |
| ~~rANS-Nx16 / fqzcomp~~ | **Pulled forward to v1.2 / M94** (2026-04-28). Required to close the M92 chr22 compression gap. |
| M40 PyPI + Maven Central publishing | Package registry publication (internal-only until further notice). |
