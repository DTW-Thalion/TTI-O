# Vendor Format Integration Guide

This document covers the status and integration patterns for vendor-specific
mass spectrometry and NMR data formats in TTI-O.

## Thermo .raw (v0.6+ — delegation to ThermoRawFileParser)

**Status:** Implemented via delegation. No proprietary code ships
with TTI-O; the converter is resolved at runtime.

Thermo Scientific `.raw` files are the most common proprietary MS
format in proteomics and metabolomics. The format has no open spec;
reading it requires Thermo's closed-source `RawFileReader` .NET
assembly. The CompOmics team's [ThermoRawFileParser][trfp] wraps
that assembly in a CLI that emits mzML / MGF / parquet. TTI-O
delegates to that CLI and parses the emitted mzML with its existing
reader.

[trfp]: https://github.com/compomics/ThermoRawFileParser

### Installation

| Method | Command | Notes |
|---|---|---|
| **Conda/Mamba** | `conda install -c bioconda thermorawfileparser` | Ships a self-contained build; no Mono needed. |
| **BioContainers (Docker)** | `docker pull quay.io/biocontainers/thermorawfileparser:<tag>` | Bind-mount your `.raw` into the container. |
| **Release tarball** | Download `ThermoRawFileParser.zip` from the [releases page][trfp-rel] | Requires Mono ≥ 6.x for the `.exe` build; the self-contained Linux build needs no extra runtime. |
| **Nix** | `nix-shell -p thermorawfileparser` | |

[trfp-rel]: https://github.com/compomics/ThermoRawFileParser/releases

### Binary resolution (all three languages)

1. Explicit argument (`thermorawfileparser=` in Python;
   `ThermoRawReader.read(path, binPath)` in Java; env var for ObjC).
2. `THERMORAWFILEPARSER` environment variable (absolute path).
3. `ThermoRawFileParser` on `PATH` (.NET 8 self-contained build).
4. `ThermoRawFileParser.exe` on `PATH` — invoked via `mono` if present.

If none resolves, the importer raises a clear error that references
this document. The input file is validated **before** the binary
lookup, so missing-file errors surface even when the parser isn't
installed.

### Command invoked

```
<binary> -i <path/to/sample.raw> -o <tmpdir> -f 2
```

`-f 2` selects mzML output. The importer reads
`<tmpdir>/<stem>.mzML` — or the first `*.mzML` in the temp dir if
naming differs — via the existing mzML parser, then deletes the temp
directory.

### What is recovered

Everything the mzML bridge preserves: scan headers (RT, MS level,
polarity, precursor m/z + charge), profile or centroid m/z + intensity
arrays, TIC/BPC chromatograms. Thermo-specific extended method
metadata that the bridge drops is not carried into TTI-O by this
path — if that data is load-bearing, you will need a direct reader.

## Bruker TDF (v0.8 M53 — shipped)

**Status:** Implemented. Metadata reads natively in all three
languages; binary frame decompression uses `opentimspy` (Python) or
delegates to the Python helper via subprocess (Java, Objective-C).

Bruker timsTOF `.d` directories hold two files:

* `analysis.tdf` — a plain SQLite database with metadata tables
  (`Frames`, `GlobalMetadata`, `Properties`, `Precursors`, ...).
* `analysis.tdf_bin` or `analysis.tdf_raw` — a binary blob with
  ZSTD-compressed frame data plus a scan-to-ion index.

The SQLite half is openly readable from the standard library in
every language. The binary half is read via the open-source
[`opentimspy`][opentimspy] Python package (wraps `libtimsdata.so` from
the paired [`opentims-bruker-bridge`][bridge] wheel). **No proprietary
Bruker SDK is involved** — the bridge wheel embeds the open reference
implementation of the documented frame format.

[opentimspy]: https://pypi.org/project/opentimspy/
[bridge]: https://pypi.org/project/opentims-bruker-bridge/

### Installation

```bash
pip install 'ttio[bruker]'
```

This pulls `opentimspy` + `opentims-bruker-bridge`. The bridge wheel
ships prebuilt `libtimsdata.so` (Linux), `libtimsdata.dylib` (macOS),
and `libtimsdata.dll` (Windows) — no extra toolchain required.

### Binary path resolution (Java + Objective-C)

The Java `BrukerTDFReader` and Objective-C `TTIOBrukerTDFReader`
read SQLite metadata natively but subprocess to Python for binary
frame data. Python interpreter lookup order:

1. `TTIO_PYTHON` environment variable (absolute path).
2. `python3` on `PATH`.
3. `python` on `PATH`.

The chosen interpreter must have `ttio[bruker]` installed.

### What is recovered

Per-peak arrays for every MS1 frame (MS2 optional via `--ms2` /
`ms2=True`):

* **m/z** — calibrated mass-to-charge, float64 Da.
* **intensity** — raw peak intensity, float64.
* **inv_ion_mobility** — inverse reduced ion mobility (1/K₀),
  float64 Vs/cm². This is the **third signal channel** added in
  v0.8 M53 — timsTOF frames are 2-D acquisitions and the ion-mobility
  axis must round-trip per-peak, not per-spectrum.

Frame-level metadata (retention time, MS level) and instrument
config (vendor, model, acquisition software) are populated from the
`Frames`, `GlobalMetadata`, and `Properties` tables.

### Round-trip verification

The Python test `tests/test_bruker_tdf.py::test_real_tdf_round_trip`
takes an optional `TTIO_BRUKER_TDF_FIXTURE` environment variable
pointing at a real Bruker `.d` directory. When set, it round-trips
through `read()` and asserts that:

* Frame count matches `analysis.tdf`'s `Frames` table.
* The written `.tio` has `mz`, `intensity`, and `inv_ion_mobility`
  signal channels with identical shapes.
* m/z and intensity match the opentimspy reference extraction.

The test is skipped when the environment variable is not set
(i.e., on CI and when no real fixture is available).

### Command invoked (Java + ObjC)

```
<python> -m ttio.importers.bruker_tdf_cli \
    --input <path/to/run.d> \
    --output <path/to/target.tio> [--title "..."] [--ms2]
```

Exit codes: 0 = success, 2 = bad args, 3 = `opentimspy` missing, 4 = I/O.

### Scope

v0.8 M53 ships metadata + binary-data extraction. MS2 precursor
threading into the compound schema and native C ports of the frame
decoder (removing the Python helper dependency for Java and ObjC)
are v0.9 concerns.

## Waters MassLynx (v0.9 M63 — delegation pattern)

**Status:** Implemented via delegation, in all three languages. No
proprietary code ships with TTI-O; the converter is resolved at
runtime exactly like the Thermo integration.

Waters MassLynx raw data lives in `.raw` **directories** (same
extension as Thermo but a directory, not a file — a common source of
confusion). The format is undocumented; access requires Waters'
`MassLynxRaw` SDK or a community CLI wrapper. TTI-O delegates to
a user-installed `masslynxraw` command (or `MassLynxRaw.exe` on
Windows / via Mono), which emits mzML; TTI-O then parses the mzML
with its existing reader.

### Installation

Waters tooling is less standardised than Thermo's; the most common
deployment patterns are:

| Method | Command / path | Notes |
|---|---|---|
| **Community `masslynxraw` wrapper** | `pip install masslynxraw` then `masslynxraw -i <raw-dir> -o <out>` | Python wrapper; invokes the SDK via cffi. Requires the Waters SDK installed separately. |
| **Waters Connect API** | Vendor-provided MSI/DMG | Ships a `MassLynxRaw.exe` or equivalent CLI. |
| **In-house scripts** | Any argv-compatible `-i <raw-dir> -o <out-dir>` CLI | Pass an explicit `converter=` argument; TTI-O invokes it verbatim. |

### Binary resolution order

Same pattern as Thermo:

1. Explicit `converter=` argument (Python / Java) or
   `converter:` kwarg (ObjC).
2. `MASSLYNXRAW` environment variable.
3. `masslynxraw` on `PATH` (native).
4. `MassLynxRaw.exe` on `PATH` — invoked via `mono` on non-Windows.

When resolution fails, all three implementations raise a
`FileNotFoundError` / `IOException` / `NSError` with a pointer to
this document.

### API

| Language | Entry point |
|---|---|
| Python | `ttio.importers.waters_masslynx.read(raw_dir, *, converter=None)` |
| ObjC | `+[TTIOWatersMassLynxReader readFromDirectoryPath:converter:error:]` |
| Java | `WatersMassLynxReader.read(String dirPath, String converter)` |

### CLI contract

TTI-O invokes the converter with:

```
<converter> -i <input-raw-dir> -o <output-dir>
```

and expects an `.tio`-ready mzML in the output directory, named
`<basename>.mzML` where `<basename>` is the `.raw` directory name with
its trailing `.raw` stripped. If the converter uses different flags,
wrap it with a small shim that accepts the TTI-O contract.

### Testing without Waters tooling

The test suites ship a POSIX shell **mock converter** that emits a
fixed-size stub mzML. This proves the delegation pipeline
(resolver → subprocess → output parse) end-to-end without installing
the vendor SDK:

* Python: `tests/integration/test_waters_masslynx.py::test_mock_converter_roundtrip`
* ObjC: `Tests/TestWatersMassLynxReader.m` (the "mock converter" set)
* Java: `WatersMassLynxReaderTest.mockConverter_roundTrip`

Set `TTIO_MASSLYNX_FIXTURE` to a real Waters `.raw` directory and
ensure `masslynxraw` is on PATH to exercise the vendor-tooling path
in nightly CI.

## JCAMP-DX 5.01 (v0.11 M73 — vibrational spectra; v0.11.1 M73.1 adds compression + UV-Vis)

**Status:** Native reader and writer in all three languages.
No external dependency; the format is plain ASCII with a documented
IUPAC spec.

JCAMP-DX is the de-facto interchange format for 1-D vibrational
spectra (Raman, IR, UV-Vis, NMR). M73 targets the AFFN
`##XYDATA=(X++(Y..Y))` dialect at spec level 5.01 — the most
interoperable subset. Writers in all three languages emit AFFN only
(bit-accurate round-trips over the byte savings). Readers in v0.11.1
additionally decode the §5.9 compressed dialects (PAC / SQZ / DIF /
DUP); see "Compressed-form reader" below.

### API

| Language | Reader | Writer |
|---|---|---|
| Python | `ttio.importers.jcamp_dx.read_spectrum(path)` | `ttio.exporters.jcamp_dx.write_raman_spectrum` / `write_ir_spectrum` / `write_uv_vis_spectrum` |
| ObjC | `+[TTIOJcampDxReader readSpectrumFromPath:error:]` | `+[TTIOJcampDxWriter writeRamanSpectrum:...]` / `writeIRSpectrum:` / `writeUVVisSpectrum:` |
| Java | `JcampDxReader.readSpectrum(Path)` | `JcampDxWriter.writeRamanSpectrum` / `writeIRSpectrum` / `writeUVVisSpectrum` |

The reader dispatches on `##DATA TYPE=`:
`RAMAN SPECTRUM` → `RamanSpectrum`;
`INFRARED ABSORBANCE` / `INFRARED TRANSMITTANCE` → `IRSpectrum`
with the matching `TTIOIRMode`; `INFRARED SPECTRUM` falls back to
`##YUNITS=` for mode detection (`ABSORBANCE` substring →
absorbance, otherwise transmittance);
`UV/VIS SPECTRUM` / `UV-VIS SPECTRUM` / `UV/VISIBLE SPECTRUM` →
`UVVisSpectrum` (v0.11.1). Any other `DATA TYPE` is rejected.

### Compressed-form reader (v0.11.1)

The reader auto-detects SQZ / DIF / DUP / PAC bodies by scanning
for compression sentinel chars (`@A-IJ-Ra-ij-rSTUVWXYZs`), excluding
`e`/`E` so AFFN scientific notation (`1.2e-3`) doesn't false-trigger.
On a compressed body, it delegates to a per-language decoder:

| Language | Decoder |
|---|---|
| Python | `ttio.importers._jcamp_decode.decode_xydata` |
| ObjC | `+[TTIOJcampDxDecode decodeLines:firstx:deltax:xfactor:yfactor:outXs:outYs:error:]` |
| Java | `global.thalion.ttio.importers.JcampDxDecode.decode` |

Each decoder implements the full SQZ alphabet (`@`, `A-I`, `a-i`),
DIF alphabet (`%`, `J-R`, `j-r`), DUP alphabet (`S-Z`, `s`), plus
the DIF Y-check convention (the repeated leading Y of the next
line, within 1e-9 of the previous line's last Y, is dropped) and
X-reconstruction from `##FIRSTX=` / `##LASTX=` / `##NPOINTS=`.

### Wire-level determinism

All three writers emit LDRs in a fixed order with `%.10g`
floating-point formatting. Given the same logical spectrum, the
three writers produce byte-identical output. This is verified by
`python/tests/integration/test_raman_ir_cross_language.py::test_jcamp_layout_is_deterministic`,
which locks the LDR prefix so a drift between language
implementations is caught in code review.

The cross-language conformance harness additionally feeds a
Python-generated `.jdx` to small subprocess drivers built on top
of the ObjC (`objc/Tools/TtioJcampDxDump`) and Java (compiled
ad-hoc into `/tmp/mpgo_m73_driver/`) readers, then compares the
parsed `(wavenumber, intensity)` arrays bit-for-bit. The tests
skip on dev boxes where the ObjC/Java sides aren't built and run
in full in CI.

### What is NOT covered

* 2-D JCAMP-DX (`##NTUPLES=`, `##PAGE=`) for imaging / 2-D NMR —
  ASCII cubes are impractical at 10–100 MB per map. Raman and IR
  cubes are stored in native HDF5 groups (`docs/format-spec.md`
  §7a) instead.
* Compressed **writers** — all three languages emit AFFN only. A
  compressed input can be round-tripped through the reader and
  re-emitted as AFFN without loss.
* Mass-spectrum `##DATA TYPE=` variants (`MASS SPECTRUM`) — the
  reader rejects them explicitly. Mass spectra round-trip through
  mzML.

### M73.1 additions (v0.11.1)

* PAC / SQZ / DIF / DUP compressed `##XYDATA=` reader in all three
  languages.
* `UVVisSpectrum` class and `UV/VIS SPECTRUM` / `UV-VIS SPECTRUM` /
  `UV/VISIBLE SPECTRUM` reader dispatch.
* `UVVisSpectrum` JCAMP-DX writer emitting `##DATA TYPE=UV/VIS
  SPECTRUM` with `##XUNITS=NANOMETERS`, `##YUNITS=ABSORBANCE`, and
  `##$PATH LENGTH CM` / `##$SOLVENT` custom LDRs.

## SAM/BAM (M87, post-v1.1.1 — delegation to samtools)

**Status:** Implemented via subprocess delegation. No `htslib`
links into TTI-O; the converter is resolved at runtime.

The Sequence Alignment/Map (SAM) text format and its binary
counterpart (BAM) are the de facto exchange format for aligned
sequencing reads. M87 ingests both via the canonical
[`samtools`][samtools] CLI and converts them to
`WrittenGenomicRun` instances ready for the M82-era write path.

[samtools]: http://www.htslib.org/

### Installation

| Method | Command | Notes |
|---|---|---|
| **apt (Debian/Ubuntu)** | `sudo apt install samtools` | Universe repo; tracks current htslib. |
| **brew (macOS)** | `brew install samtools` | Bottled; depends on `htslib`. |
| **conda/mamba** | `conda install -c bioconda samtools` | Self-contained. |
| **From source** | Build htslib + samtools from [GitHub][samtools-src] | For non-x86 architectures. |

[samtools-src]: https://github.com/samtools/samtools

The reference development environment for M87 used `samtools 1.19.2` (Ubuntu 24.04 universe build). Other 1.x versions are expected to work — TTI-O parses the SAM text emitted by `samtools view -h` and does not consume any version-specific output.

### Binary resolution (all three languages)

The reader checks `samtools` on `PATH` at first use of `to_genomic_run()` (NOT at import time per Binding Decision §135 in the M87 plan). If missing, raises a clear error including install guidance for apt / brew / conda. The error fires only on the first import operation; the `BamReader` class itself is loadable on systems without `samtools` (e.g. for documentation generation, type checking).

### Command invoked

```
samtools view -h <path> [region]
```

stdout is consumed line-by-line by the language-native subprocess wrapper:

* Python: `subprocess.Popen(..., stdout=PIPE, stderr=PIPE, text=True)`.
* ObjC: `NSTask` with `setLaunchPath:` / `setStandardOutput:`.
* Java: `ProcessBuilder([...]).start()`.

### Field mapping (SAM column → TTI-O field)

| Col | SAM Name | TTI-O destination                                          |
|-----|----------|------------------------------------------------------------|
| 1   | QNAME    | `read_names[i]`                                            |
| 2   | FLAG     | `flags[i]` (uint32)                                        |
| 3   | RNAME    | `chromosomes[i]` (or `"*"` for unmapped)                   |
| 4   | POS      | `positions[i]` (int64; SAM 1-based; 0 = unmapped)          |
| 5   | MAPQ     | `mapping_qualities[i]` (uint8)                             |
| 6   | CIGAR    | `cigars[i]` (literal `"*"` preserved if absent)            |
| 7   | RNEXT    | `mate_chromosomes[i]` (`"="` expanded to RNAME)            |
| 8   | PNEXT    | `mate_positions[i]` (int64)                                |
| 9   | TLEN     | `template_lengths[i]` (int32; signed)                      |
| 10  | SEQ      | concatenated into `sequences` byte array (`"*"` → 0 bytes) |
| 11  | QUAL     | concatenated into `qualities` byte array                   |

Columns 1–11 are parsed into the model; the optional tag fields (12+) are kept as each read's tab-joined tag text (`tags[i]`, M101) and stored with the SAM_TAGS codec. The `sequences` and `qualities` byte arrays are concatenated across all reads with `offsets[i]` / `lengths[i]` parallel arrays giving each read's slice.

### Header line handling

* `@SQ` (sequence dictionary) — `SN:` populates a chromosome list; `reference_uri` is set to the first `SN:` value.
* `@RG` (read group) — `SM:` → `sample_name`, `PL:` → `platform`. First `@RG` line wins for multi-`@RG` BAMs.
* `@PG` (program/tool) — each entry becomes a `ProvenanceRecord` accessible via `reader.lastProvenance()` (Java) / `reader.provenanceRecords` (ObjC) / set on the returned `WrittenGenomicRun.provenance_records` (Python). Note that `samtools view -bS` and `samtools view -h` themselves inject `@PG` records on the BAM/SAM stream — the `provenance_count` for the M87 fixture is 3 (1 user-supplied `bwa` plus 2 samtools-injected entries).
* `@HD` (header version) and `@CO` (comments) are read but not mapped to TTI-O fields.

### Region filter

The `region` parameter is passed verbatim to `samtools view` as a positional argument:

```python
reader = BamReader("sample.bam")
chr1_only = reader.to_genomic_run(name="chr1", region="chr1")
window   = reader.to_genomic_run(name="window", region="chr1:1000-2000")
unmapped = reader.to_genomic_run(name="unmapped", region="*")
```

The BAM must be indexed (a `.bai` companion file alongside the `.bam`) for region queries to work. Build with `samtools index <file.bam>` if the index is missing.

### Round trip via the M82 write path

```python
from ttio.spectral_dataset import SpectralDataset
from ttio.importers.bam import BamReader

reader = BamReader("alignments.bam")
run = reader.to_genomic_run(name="sample1", sample_name="NA12878")

SpectralDataset.write_minimal(
    "out.tio",
    title="Sample 1",
    isa_investigation_id="ISA-001",
    genomic_runs={run.name: run},
)
```

The resulting `.tio` accepts any of the M83–M86 codec choices on the genomic signal channels via `WrittenGenomicRun.signal_codec_overrides`.

### Cross-language conformance

A small canonical fixture (`m87_test.sam` + `m87_test.bam` + `m87_test.bam.bai`, 10 reads on two chromosomes with mixed mapped/unmapped/clipped) is committed under `python/tests/fixtures/genomic/`. Each language ships a `bam_dump` CLI emitting a fixed canonical-JSON shape (sorted keys, 2-space indent, MD5 fingerprints for the byte buffers). The M87 cross-language harness (`python/tests/integration/test_m87_cross_language.py`) runs all three CLIs on the same BAM and asserts byte-identical output. Java and ObjC currently match Python byte-exact at 1341 bytes per the canonical-JSON serialisation.

### What is NOT covered

* **CRAM input** — requires a reference FASTA. M88 (separate milestone) handles CRAM import.
* **BAM/SAM writing** — TTI-O is the read direction in M87. M88 covers BAM/CRAM writers.
* **BAM integer tag widths** — the optional fields are kept as the SAM text samtools prints (M101), where every integer type (`c C s S i I`) reads as `i`; the width is chosen again when a BAM is written.
* **Multi-`@RG` aggregation** — only the first `@RG` is parsed. Caller can override `sample_name=` if needed.
* **htslib direct linking** — subprocess via `samtools` is the Phase 5 design choice. A future optimisation milestone could add htslib-Java / pysam fast paths if the subprocess startup overhead (~50 ms per import) becomes a bottleneck.
* **Whole-run materialisation** — `to_genomic_run` builds the full `WrittenGenomicRun` in memory. Inputs that should not be held whole go through the streaming surface instead: `BamReader.iter_batches` / `stream_source` feed a `GenomicStreamWriter` batch by batch (see `docs/genomic-runs.md` §2.5, §6).

## CRAM (M88, post-M87 — reference-aware delegation to samtools)

**Status:** Implemented via subprocess delegation. CRAM extends the
M87 SAM/BAM read path with a mandatory reference FASTA argument
(Binding Decision §139); no `htslib` is linked into TTI-O.

The CRAM (Compressed Reference-oriented Alignment Map) format is
the de facto reference-compressed successor to BAM. M88 ingests it
via the same `samtools` CLI as M87, with `--reference <fasta>`
injected into the `samtools view` command line. The reader subclass
(`CramReader` / `TTIOCramReader`) extends the M87 BAM reader and
reuses every line of its SAM-text parsing path.

### Command invoked

```
samtools view -h --reference <fasta> <cram> [region]
```

Reference passed verbatim through the language-native subprocess
wrapper. CRAM cannot be decoded without a matching reference; the
reader rejects construction without one (Python: `TypeError`; ObjC:
the single-arg `initWithPath:` selector is `NS_UNAVAILABLE`; Java:
the inherited single-arg ctor delegates to `BamReader` rather than
`CramReader`, and the two-arg ctor null-rejects the FASTA path).

### Constructor signature (all three languages)

| Language | Signature |
|---|---|
| **Python** | `CramReader(path: str \| Path, reference_fasta: str \| Path)` |
| **Objective-C** | `-[TTIOCramReader initWithPath:referenceFasta:]` |
| **Java** | `new CramReader(Path path, Path referenceFasta)` |

### Reference index (`.fai`) handling

`samtools` autogenerates an `<fasta>.fai` index alongside the FASTA
on first read. The M88 conformance fixture commits the reference
FASTA but **NOT** the autogenerated index — let samtools regenerate
it locally. This avoids spurious diffs across machines (the index
encodes byte offsets that change with line-ending normalisation).

### Field mapping

Identical to SAM/BAM (M87). The reader produces a
`WrittenGenomicRun` with the same shape as `BamReader`. CRAM-specific
features (reference-relative encoding, slice-level metadata) are
opaque; samtools materialises them as plain SAM text before our
parser sees them.

### CLI dispatch via `bam_dump --reference` (M88.1)

The M87 `bam_dump` / `TtioBamDump` / `BamDump` CLI in each language
auto-dispatches to the CRAM reader when the input path ends in
`.cram` (case-insensitive). The reference FASTA is supplied via a
`--reference <fasta>` flag; omitting it on a `.cram` path is a
hard error (exit 2) with a stderr message identifying the missing
flag. For BAM/SAM paths the flag is accepted but unused (defensive
— scripts can always pass it). The canonical-JSON output schema is
unchanged from M87.

## SAM/BAM/CRAM Export (M88, post-M87)

**Status:** Implemented in all three languages. Subprocess
delegation to `samtools view -b` / `samtools view -C` / `samtools
sort`; no `htslib` in TTI-O.

M88 closes the read/write loop: a `WrittenGenomicRun` produced by
the M82 write path or the M87 BAM importer can be re-emitted as
SAM, BAM, or CRAM. The writer composes SAM text in memory from
the `WrittenGenomicRun` field arrays, then pipes it through one or
two `samtools` invocations to produce the final binary.

### Pipeline

* **BAM**: `(SAM text) | samtools view -b -o <out.bam>`
* **CRAM**: `(SAM text) | samtools view -C --reference <fasta> | samtools sort -O cram --reference <fasta> -o <out.cram>`

The CRAM path uses two `samtools` invocations because `samtools
view -C` produces unsorted CRAM and CRAM slices require sorted
input for efficient reference-relative encoding. The intermediate
pipe is in-process: ObjC chains two `NSTask`s with an `NSPipe`,
Java uses two `ProcessBuilder` processes with a pump thread
(`Process.getInputStream().transferTo(...)`), Python uses two
`subprocess.Popen` instances chained via `stdin=prev.stdout`.

### Writer constructor signatures

| Format | Python | Objective-C | Java |
|---|---|---|---|
| **BAM** | `BamWriter(path)` | `-[TTIOBamWriter initWithPath:]` | `new BamWriter(Path path)` |
| **CRAM** | `CramWriter(path, reference_fasta)` | `-[TTIOCramWriter initWithPath:referenceFasta:]` | `new CramWriter(Path path, Path referenceFasta)` |

### SAM text emission rules (Binding Decisions §136 + §138)

* **All 11 SAM columns are always emitted.** Sentinel values fill
  unset fields: `*` for QNAME / RNAME / CIGAR / MRNM / SEQ / QUAL,
  `0` for POS / MAPQ / MPOS / ISIZE.
* **RNEXT collapse** (§136): mate chromosome is written as `=` only
  when the read is mapped (`chromosome != "*"`) AND the mate
  chromosome equals the read chromosome. Unmapped reads always emit
  `*` for RNEXT regardless of the stored mate-chromosome value.
* **PNEXT mapping** (§138): negative `mate_position` values map to
  `0` on emit (SAM spec disallows negative POS / PNEXT).
* **QUAL pass-through**: ASCII Phred+33 verbatim. No arithmetic at
  read or write — the bytes that came in via samtools are the same
  bytes that go back out.
* **Read group**: a single `@RG ID:rg1` line is emitted, populated
  with the run's `sample_name` (`SM:`) and `platform` (`PL:`).
* **Unknown chromosome lengths**: `@SQ LN:` fallback is `2147483647`
  (SAM `LN` is a 31-bit signed integer; this is the max value, used
  when the input genomic run does not carry per-chromosome lengths).

### Round-trip semantics

| Round trip | Lossless? | Caveats |
|---|---|---|
| **BAM → `.tio` → BAM** | Yes (per-field) | All 11 SAM columns preserved; tag fields (NM, MD, etc.) discarded — they were never read in M87. |
| **CRAM → `.tio` → CRAM** | Sequence buffer + names | Quality bytes preserved verbatim; per-field positions equal. CRAM's reference-relative encoding round-trips through the SAM text intermediate without loss. |
| **BAM → `.tio` → CRAM** | Yes | Same as BAM, with reference-relative re-compression on emit. |
| **CRAM → `.tio` → BAM** | Yes | Same as CRAM read; output is plain BAM. |

The M88 conformance suite (`test_m88_cram_bam_round_trip`, all
three languages) verifies each of these matrix entries against the
canonical M88 fixture (5 reads, 2 chromosomes, multi-RG, `M88_TEST_SAMPLE`).

### Cross-language conformance

The M88 cross-language harness
(`python/tests/integration/test_m88_cross_language.py`) drives both
read paths against the M88 fixtures and asserts byte-identical
canonical JSON across Python / ObjC / Java:

* **BAM** (3 tests) — the existing M87 `bam_dump` CLIs run against
  `m88_test.bam`; output is 1341-ish bytes of canonical JSON
  (sorted keys, 2-space indent, MD5 fingerprints for the byte
  buffers).
* **CRAM** (3 tests, added in M88.1) — the same `bam_dump` /
  `TtioBamDump` / `BamDump` CLIs are re-used with a new
  `--reference <fasta>` flag; the path's `.cram` extension
  auto-dispatches to `CramReader` / `TTIOCramReader` / Java
  `CramReader`. Output for the M88 CRAM fixture is exactly 914
  bytes (md5 `2be5c5bccc95635240aa60337406cb35`) and is
  byte-identical across all three languages. The CRAM fixture's
  `provenance_count` differs from the BAM fixture's by design —
  samtools injects different `@PG` records on each format's read
  path, and that's correct.

### What is NOT covered

* **Optional SAM tag fields** — still ignored in v0 (writers cannot
  emit tags they never read).
* **Multi-`@RG` aggregation on emit** — a single `@RG ID:rg1` is
  always emitted; if the input had multiple read groups they were
  collapsed at import time.
* **CRAM 4.0** — samtools versions <1.20 emit CRAM 3.x. The format
  emitted depends on the installed samtools version; TTI-O does
  not pin or override it.

---

## FASTA (reference + unaligned reads)

The FASTA importer / exporter pair is **native to TTI-O** — no
`samtools` subprocess, no third-party FASTA library. Pure-text
parsers in each language read FASTA and FASTQ records directly.
gzip-compressed input is auto-detected via the `1f 8b` magic bytes
regardless of file extension.

### Two import modes

The same `FastaReader` services two distinct workflows:

| Mode | Method | Result | Use case |
|---|---|---|---|
| **Reference** | `read_reference()` / `+readReferenceFromPath:uri:error:` / `readReference()` | `ReferenceImport` (chromosomes + content MD5) | Pair with BAM/CRAM input. Embed under `/study/references/<uri>/` so the `.tio` is self-contained — REF_DIFF_V2 decode no longer needs an external `.fa` on disk. |
| **Unaligned** | `read_unaligned()` / `+readUnalignedFromPath:...:error:` / `readUnaligned(...)` | `WrittenGenomicRun` with SAM unmapped sentinels | Amplicon panels, target lists, quality-stripped reads. `flags=4` (unmapped), `chrom="*"`, `pos=0`, `mapq=255`, `cigar="*"`, qualities filled with the `0xFF` "qualities unknown" sentinel matching `BamReader`'s convention. |

### Reference-FASTA embedding layout

```
/study/references/<uri>/
  attr "reference_uri"   = <URI string>
  attr "md5"             = <32-char lowercase hex string>
  chromosomes/
    <chrom_name>/
      attr "length"      = <int64 sequence length in bases>
      data               = uint8 dataset, zlib-compressed, uppercase ACGTN
    ...
```

`ReferenceImport.write_to_dataset()` and the REF_DIFF_V2 auto-embed
path emit the **same** layout — `format-spec.md` §10.10 documents
the canonical shape and is the single source of truth. (Pre-v1.1.1
`write_to_dataset()` emitted a flat `<group>/<chromName>` variant
with `@total_bases` and case-preserving lowercase data; that variant
was unified in v1.1.1 — see CHANGELOG and `format-spec.md` for the
historical note.)

The MD5 attribute is computed by sorting chromosomes by name (so
the digest is order-invariant), then digesting the concatenated
sequence bytes: `MD5(seq_for_chr_a || seq_for_chr_b || ...)` where
`||` denotes byte concatenation. Unified in v1.1.0 (previously
some writers used the `(name + 0x0A + seq + 0x0A)` form; see the
migration guide and CHANGELOG for details). Cross-language
byte-equal: Python, Java, and ObjC all agree.

The `ReferenceResolver` lookup chain is:
1. Embedded reference at `/study/references/<uri>/` in the open
   `.tio`.
2. External FASTA via the `REF_PATH` env var (or an explicit path
   passed to the resolver constructor).
3. `RefMissingError` — hard-fail rather than partial decode.

### Exporter — line wrap, `.fai`, gzip

`FastaWriter` writes either a `ReferenceImport` (chromosomes) or a
`WrittenGenomicRun` / read-side `GenomicRun` (one record per read).
Output guarantees:

* **Header**: exactly `>name\n` — no description preserved.
* **Line wrap**: configurable, default 60 chars (NCBI convention).
  Last line of each record may be shorter; the index records the
  canonical wrap width.
* **Line endings**: LF only.
* **`.fai` index**: emitted alongside non-gzip output by default.
  Format: `<name>\t<length>\t<offset>\t<linebases>\t<linewidth>\n`
  matching samtools' `faidx` convention. Skipped silently for gzip
  output (samtools requires bgzip for indexed gzip-FASTA, which
  TTI-O does not yet emit).
* **gzip**: enabled when the destination ends in `.gz` or
  `gzipOutput=True`.

### Cross-language byte-equality

For uncompressed output, all three languages produce **byte-
identical** FASTA bytes for the same input — proven by
`python/tests/integration/test_fasta_fastq_cross_language.py` (the
test invokes Python directly, then subprocesses the Java
`FastaRoundTrip` and ObjC `TtioFastaRoundTrip` CLIs and diffs the
three outputs). gzip wrappers are deflate-non-deterministic across
implementations; the inner uncompressed payload is byte-equal.

### Round-trip through `.tio`

The full chain `FASTA → WrittenGenomicRun → SpectralDataset.write_minimal
→ SpectralDataset.open → FastaWriter.write_run → FASTA` produces
byte-identical output, verified by
`test_fasta_fastq_tio_roundtrip.py` (Python),
`FastaFastqTioRoundTripTest` (Java), and
`TestFastaFastqTioRoundTrip` (ObjC).

### CLI

| Language | Tool | Notes |
|---|---|---|
| Python | `python -m ttio.tools.fasta_import_cli {reference,unaligned}` <br> `python -m ttio.tools.fasta_export_cli {reference,run}` | Production CLIs. |
| Java | `global.thalion.ttio.tools.FastaRoundTrip <in.fa> <out.fa> [line_width]` | Round-trip CLI used by the conformance harness. |
| ObjC | `TtioFastaRoundTrip <in.fa> <out.fa> [line_width]` | Same. |

---

## FASTQ (unaligned reads)

The FASTQ importer / exporter pair complements FASTA's unaligned-
mode read path with quality scores preserved verbatim.

### Phred auto-detect

`FastqReader` inspects the qualities byte range to choose the
encoding offset:

| Observed range | Decision |
|---|---|
| any byte `< 59` | **Phred+33** (Phred+64 starts at byte 64; modern Illumina / Sanger range). |
| every byte in `[64, 104]` | **Phred+64** (legacy Illumina / Solexa pre-1.8). |
| otherwise | **Phred+33** (default). |

The detected source offset is exposed on the reader
(`detected_phred_offset` / `detectedPhredOffset()` /
`outDetected:` out-parameter). Override the heuristic with
`force_phred=33|64`. Internal storage normalises to Phred+33 ASCII
verbatim — the byte you store is the byte the writer emits, no
arithmetic mid-pipeline.

### SAM-unmapped sentinels

Each FASTQ record becomes one read in a `WrittenGenomicRun` with
`flags=4`, `chromosome="*"`, `position=0`, `mapping_quality=255`,
`cigar="*"`. Sequence and qualities are populated from the FASTQ
record verbatim. Mate fields default to BAM unmapped sentinels
(`mate_position=-1`, `template_length=0`).

### Exporter — Phred selection, gzip

`FastqWriter.write` defaults to Phred+33 output. Pass
`phred_offset=64` for legacy Illumina output (each byte rewritten
as `b + 31` on the way out). The internal `0xFF` "qualities
unknown" sentinel is mapped to Phred 0 (`!`) on output so the
result is always a parseable FASTQ.

gzip is enabled when the destination ends in `.gz`.

### Cross-language byte-equality

Same guarantee as FASTA: uncompressed output is byte-equal across
Python / Java / ObjC for the same input. Verified by the same
3-way harness.

### CLI

| Language | Tool |
|---|---|
| Python | `python -m ttio.tools.fastq_import_cli` <br> `python -m ttio.tools.fastq_export_cli` |
| Java | `global.thalion.ttio.tools.FastqRoundTrip <in.fq> <out.fq>` |
| ObjC | `TtioFastqRoundTrip <in.fq> <out.fq>` |

### What is NOT covered

* **Phred+64 round-trip preservation** — internal storage is
  Phred+33 only. A Phred+64 input is converted on read; on export,
  the user explicitly opts back into Phred+64 with the
  `phred_offset` argument if needed.
* **Multi-line FASTQ** — most FASTQ aligners produce one-line
  sequence + one-line qualities. Multi-line records (rarely seen
  in modern data) are not supported.
* **bgzip-FASTQ** — gzip output is plain gzip, not bgzip; samtools
  random-access via `.fqi` will not work.
* **Read pairing in a single FASTQ** — paired-end FASTQ is two
  separate files (`R1.fq` + `R2.fq`); the importer treats each as
  an independent unaligned run.
