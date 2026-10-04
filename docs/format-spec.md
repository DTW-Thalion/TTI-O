# TTI-O `.tio` File Format Specification — v0.10.0

This document specifies the on-disk layout of an `.tio` file as written
by libTTIO v0.10.0. It is detailed enough for a reader implemented in a
different language (Python, Java, Rust, Go) to open, validate, and fully
decode a file without consulting the reference source. Three
interoperable implementations (ObjC, Python, Java) read and write this
format; SQLite and Zarr ship as alternative chunked-array container
backends (see `docs/providers.md`).

Each point release is a strict superset of the previous release on
disk:

- **v0.3** — compound per-run provenance, `v2:` canonical signatures,
  LZ4 / Numpress-delta compression codecs.
- **v0.4** — envelope encryption + key rotation, spectral anonymization,
  nmrML writer, chromatogram API.
- **v0.7** — `ttio_format_version` bumps from `"1.1"` to `"1.2"`; the
  versioned wrapped-key blob (§10b) replaces the fixed 60-byte v1.1
  blob; `read_canonical_bytes` becomes the byte-level contract for
  signatures and encryption (§10c).
- **v0.8** — post-quantum crypto preview (`opt_pqc_preview`:
  ML-KEM-1024 for KEM, ML-DSA-87 for signatures); see `docs/pqc.md`.
- **v0.9** — provider abstraction hardening: SQLite and Zarr v3
  backends in all three languages.
- **v0.10** — streaming transport layer (`.tis`; see
  `docs/transport-spec.md`) and v1.0 per-Access-Unit encryption
  (`opt_per_au_encryption`, optional `opt_encrypted_au_headers`).
  Adds the `VL_BYTES` compound field kind and the
  `<channel>_segments` / `spectrum_index/au_header_segments`
  compound layouts (§9.1). The file format version remains `"1.2"`;
  per-AU encryption is additive and feature-flagged.

Every feature added after v0.2 is gated by a feature flag (see
`docs/feature-flags.md`) so readers can detect capability support at
open time.

A conforming `.tio` file is a plain HDF5 file (format 1.x) with the
group, dataset, and attribute hierarchy described below. Any program
that can read HDF5 can introspect an `.tio` file with `h5dump`.

---

## 1. Versioning

Every v1.0+ file carries two attributes on the root group `/`:

| Attribute                | Type              | Value                              |
|--------------------------|-------------------|------------------------------------|
| `ttio_format_version`  | fixed-len string  | `"1.0"` (current; all v1.0.0 writers stamp this unconditionally) |
| `ttio_features`        | fixed-len string  | JSON array of feature strings      |

**Version semantics:**

- **Major** (`1.x -> 2.0`): backward-incompatible layout changes.
  Readers that only support 1.x must refuse to open 2.0 files.
- **Minor** (`1.0 -> 1.1`): backward-compatible additions. Existing
  readers may ignore new attributes/datasets they do not understand.

v1.0 is the first stable release. Pre-v1.0 development files
(historical `"0.x"` stamps) were never publicly released; v1.0+
readers are not expected to ingest them.

**Feature flags** are documented in
[`feature-flags.md`](feature-flags.md). Features without an `opt_`
prefix are **required**: per the spec, a reader must refuse to
open a file that lists a required feature it does not support.
`opt_`-prefixed features are informational; readers may ignore
them. *Implementation note:* the current v1.0 reference readers
parse the feature list but do not yet enforce the
required-feature-refusal gate — that enforcement is tracked as a
v1.x follow-up.

---

## 2. Top-level layout

```
/                                       (root group)
├── @ttio_format_version              ("1.0")
├── @ttio_features                    (JSON array)
├── @encrypted                          ("aes-256-gcm") — optional
├── @access_policy_json                 (JSON) — optional
├── study/                              (group, required)
│   ├── @title                          (string)
│   ├── @isa_investigation_id           (string)
│   ├── @transitions_json               (string) — optional
│   ├── ms_runs/                        (group, always present)
│   │   ├── @_run_names                 (comma-separated run names)
│   │   └── <run_name>/                 (one group per run)
│   ├── nmr_runs/                       (group, always present)
│   │   ├── @_run_names                 (comma-separated run names)
│   │   └── <run_name>/                 (legacy NMR runs)
│   ├── genomic_runs/                   (group) — optional, M82+ genomic data
│   │   ├── @_run_names                 (comma-separated run names)
│   │   └── <run_name>/                 (per-genomic-run group; see §10)
│   ├── references/                     (group) — optional, M93+ reference embedding
│   │   └── <reference_uri>/            (per-reference group; FASTA + chromosomes)
│   ├── subjects/                       (group) — optional, v0.11 cohort metadata (see §11)
│   │   └── <external_id>/              (per-Subject group)
│   ├── samples/                        (group) — optional, v0.11 cohort metadata (see §11)
│   │   └── <sample_id>/                (per-Sample group)
│   ├── identifications                 (compound dataset) — optional
│   ├── quantifications                 (compound dataset) — optional
│   ├── provenance                      (compound dataset) — optional
│   ├── identifications_sealed          (encrypted byte blob) — optional
│   ├── identifications_sealed_iv       (3 × int32) — optional
│   ├── identifications_sealed_tag      (4 × int32) — optional
│   ├── identifications_sealed_bytes    (1 × uint32) — optional
│   ├── quantifications_sealed          (encrypted byte blob) — optional
│   ├── quantifications_sealed_iv       (3 × int32) — optional
│   ├── quantifications_sealed_tag      (4 × int32) — optional
│   ├── quantifications_sealed_bytes    (1 × uint32) — optional
│   └── image_cube/                     (group) — optional, MSImage only
└── protection/                         (group) — optional, M70+ envelope encryption
    ├── key_info/                       (group; wrapped-DEK blob + KEK metadata; see §5b)
    └── access_policies/                (compound dataset) — optional
```

Runs under `/study/ms_runs/` may contain any spectrum subclass, not
just mass spectra, despite the legacy name. NMR runs written via
`TTIOSpectralDataset.msRuns` (the post-M10 idiom) land here. The
`/study/nmr_runs/` dict is retained for backward compatibility with
files written by pre-M10 code.

---

## 3. Acquisition run layout

Each `/study/ms_runs/<name>/` group contains:

```
<run_name>/
├── @acquisition_mode                   (int64; TTIOAcquisitionMode enum)
├── @spectrum_count                     (int64)
├── @spectrum_class                     (string)  e.g. "TTIOMassSpectrum"
├── @nucleus_type                       (string) — optional, NMR only
├── @provenance_json                    (string) — optional
├── @provenance_signature               (base64 HMAC) — optional
├── _spectrometer_freq_mhz              (1 × float64) — optional, NMR only
├── instrument_config/                  (group with string attrs)
├── spectrum_index/                     (group, described in §4)
├── signal_channels/                    (group, described in §5)
└── chromatograms/                      (group, optional, M24 v0.4, described in §5a)
```

**Spectrum classes** currently recognized:

| String                  | Meaning                                           |
|-------------------------|---------------------------------------------------|
| `TTIOMassSpectrum`      | MS spectra; required channels `mz` + `intensity` |
| `TTIONMRSpectrum`       | 1-D NMR; required channels `chemical_shift` + `intensity` |

Absence of `@spectrum_class` triggers v0.1 fallback: the reader
assumes `TTIOMassSpectrum` and hardcoded `mz_values`/`intensity_values`
channel names.

**Instrument config** holds only string attributes
(`manufacturer`, `model`, `serial_number`, `source_type`,
`analyzer_type`, `detector_type`), any of which may be empty.

---

## 3a. Run modality (M79, v0.11)

Each run group MAY carry an optional `@modality` UTF-8 string
attribute identifying the omics modality the run represents. The
attribute is purely informational at v0.11 — readers continue to
dispatch on `@spectrum_class` for record decoding — but it scopes
which downstream metadata + analytics apply to the run.

| `@modality`           | Meaning                                                                       |
|-----------------------|-------------------------------------------------------------------------------|
| `mass_spectrometry`   | Default. Mass-spec, NMR, vibrational, UV-Vis runs (every v0.10 / pre-v0.11 file). |
| `genomic_sequencing`  | Genomic short-read / long-read runs. Reserved for the v0.11 genomic milestones (M74+). |

Absence of `@modality` MUST be interpreted as `mass_spectrometry`
so v0.10 files load unchanged. Future modality strings (proteomic,
metabolomic, …) MAY be added without a format-version bump because
unrecognised values surface as the literal string at the API layer.

---

## 4. `spectrum_index/`

Parallel 1-D datasets, one per field. All datasets have length
`spectrum_count`.

| Dataset               | Type            | Semantics                              |
|-----------------------|-----------------|----------------------------------------|
| `offsets`             | uint64[N]       | **OPT-OUT in v1.10+** — see §4a; computed from `cumsum(lengths)` on read |
| `lengths`             | uint32[N]       | Number of elements per spectrum        |
| `retention_times`     | float64[N]      | Scan time in seconds                   |
| `ms_levels`           | int32[N]        | MS level (1, 2, …); 0 for non-MS       |
| `polarities`          | int32[N]        | 1 = +, -1 = -, 0 = unknown             |
| `precursor_mzs`       | float64[N]      | Precursor m/z; 0 for non-tandem        |
| `precursor_charges`   | int32[N]        | Precursor charge state; 0 if unknown   |
| `base_peak_intensities`| float64[N]     | Max intensity per spectrum             |
| `headers`             | compound[N]     | Optional `opt_compound_headers` view   |

### 4a. `offsets` redundancy elimination (v1.10 #10)

`offsets[i]` is mathematically `sum(lengths[0..i])`, so storing
both `offsets` and `lengths` was pure on-disk redundancy. v1.10+
writers omit the `offsets` column from `spectrum_index/`,
`chromatogram_index/`, and `genomic_index/` by default. Readers
synthesize `offsets` from `cumsum(lengths)` at load time using a
uint64 accumulator (overflow-safe on >4 GB genomic runs even when
`lengths` is uint32).

The omission is unconditional in v1.10+ writers — there is no
opt-out to restore the redundant column. Pre-v1.10 readers cannot
open v1.10-default files.

All parallel datasets are chunked (chunk = 1024) and `zlib -6`
compressed, with the HDF5 byte-shuffle filter ahead of the
compressor on multi-byte element types (current writers; the
filter is part of core HDF5 and self-describing, so any reader
decodes it transparently and files written without it remain
valid). `@count` (int64) on the `spectrum_index/` group mirrors
`spectrum_count` for quick access.

### 4b. Pixel coordinates (M102, `opt_pixel_coordinates`)

A run whose spectra are the pixels of an imaging acquisition (the
imzML importer's `imzml_pixels` run) carries each spectrum's spatial
position as three more parallel columns:

| Dataset    | Type       | Semantics                                   |
|------------|------------|---------------------------------------------|
| `pixel_x`  | int32[N]   | Pixel x position, as the source numbered it |
| `pixel_y`  | int32[N]   | Pixel y position                            |
| `pixel_z`  | int32[N]   | Pixel z position; `1` for a 2-D acquisition |

The three columns are present together or not at all; a reader that
finds only some of them treats the file as malformed. They use the
chunking and compression of the other `spectrum_index/` columns.
imzML numbers pixels from 1 (`IMS:1000050` / `IMS:1000051` /
`IMS:1000052`), and the importers store those numbers unchanged, so
an export writes back the coordinates it read. Values are
non-negative; the transport MSImagePixel extension carries the same
numbers as uint32 (transport-spec §4.3.1).

Writers add the `opt_pixel_coordinates` feature flag when any run
carries the columns, so files without them stay byte-identical to
pre-M102 output.

**Older files.** Before M102 the Python imzML importer kept the
coordinates in a provenance parameter, `imzml_pixel_coordinates_csv`,
on both the run-level and the dataset-level record: one `x,y,z`
triple per spectrum in spectrum order, triples joined by `;`. That
parameter is no longer written. Both provenance levels are mirrored
into a fixed-length `@provenance_json` attribute, and HDF5 limits an
attribute to 64 KB, so an import failed once the parameter passed
that size: from about 9,000 pixels on a 100-pixel-wide grid. A reader that needs coordinates for a run without
the columns looks for the parameter in the run's provenance records,
then in the dataset-level records, and uses it only when it lists
exactly one triple per spectrum.

The other imzML scalars stay in the run- and dataset-level provenance
parameters as before: `imzml_mode`, `imzml_uuid_hex`,
`imzml_grid_max_x` / `_y` / `_z`, `imzml_pixel_size_x` / `_y` and
`imzml_scan_pattern`. An imzML exporter that writes a pixel run reads
the mode, UUID, grid, pixel size and scan pattern from them when they
are present.

The **optional compound `headers` dataset** packs all of the above
into one rank-1 dataset of compound records for external tooling
readability. Layout:

```
struct {
    uint64   offset;
    uint32   length;
    float64  retention_time;
    uint8    ms_level;
    int8     polarity;
    float64  precursor_mz;
    int32    precursor_charge;
    float64  base_peak_intensity;
}
```

---

## 5. `signal_channels/`

Name-driven channel storage.

| Attribute / dataset        | Type              | Notes                            |
|----------------------------|-------------------|----------------------------------|
| `@channel_names`           | fixed-len string  | Comma-separated channel names    |
| `<channel>_values`         | float64[N_total]  | Concatenation of all spectra     |
| `<channel>_values_encrypted` | int32[M]        | Present if channel is encrypted  |
| `<channel>_iv`             | int32[3]          | 12-byte IV (if encrypted)        |
| `<channel>_tag`            | int32[4]          | 16-byte GCM auth tag             |
| `@<channel>_ciphertext_bytes` | int64          | Exact ciphertext size            |
| `@<channel>_original_count`| int64             | Original element count           |
| `@<channel>_algorithm`     | string            | e.g. `"AES-256-GCM"`             |
| `@ttio_signature`          | string (base64)   | HMAC-SHA256 if signed            |

For an MS run the channels are `{mz, intensity}` producing
`mz_values` and `intensity_values`. For an NMR run they are
`{chemical_shift, intensity}` producing `chemical_shift_values` and
`intensity_values`. Channel data is chunked (chunk = 16384) with
`zlib -6` behind the HDF5 byte-shuffle filter (see §4).

Each channel's concatenated dataset is indexed by `offsets[i]` and
`lengths[i]` from `spectrum_index/`: spectrum `i`'s contents for
channel `c` is `<c>_values[offsets[i] .. offsets[i] + lengths[i])`.

---

## 5a. Chromatograms (M24, v0.4)

The `chromatograms/` group is **optional** — v0.3 files lack it and
readers return an empty list. When present:

```
chromatograms/
├── @count                              (int64, number of chromatograms)
├── time_values                         (float64[total_points], concatenated)
├── intensity_values                    (float64[total_points], concatenated)
└── chromatogram_index/                 (subgroup, parallel metadata arrays)
    ├── offsets                          (int64[count])
    ├── lengths                          (uint32[count])
    ├── types                            (int32[count]; TTIOChromatogramType enum)
    ├── target_mzs                       (float64[count]; XIC target m/z, 0.0 otherwise)
    ├── precursor_mzs                    (float64[count]; SRM precursor, 0.0 otherwise)
    └── product_mzs                      (float64[count]; SRM product, 0.0 otherwise)
```

Chromatogram `i`'s time/intensity slice is
`time_values[offsets[i] .. offsets[i] + lengths[i])`.

**TTIOChromatogramType enum**: 0 = TIC, 1 = XIC, 2 = SRM.

---

## 5b. Envelope encryption key info (M25, v0.4; v1.2 blob in v0.7)

The `opt_key_rotation` feature flag gates the `/protection/key_info/`
group. When present:

```
/protection/key_info/
├── @kek_id                             (string, caller-supplied KEK identifier)
├── @kek_algorithm                      (string, "aes-256-gcm" default; identifies the KEK cipher)
├── @wrapped_at                         (string, ISO-8601 timestamp)
├── @key_history_json                   (string, JSON array of prior entries)
├── @dek_wrapped_bytes                  (int64, actual blob length — v0.7+, see below)
└── dek_wrapped                         (uint8[N]; layout depends on the blob version)
```

### 5b.1 v1.1 wrapped-key layout (pre-v0.7 writers)

Fixed 60 bytes, AES-256-GCM-only:

```
offset  len  field
0       32   AES-256-GCM ciphertext (wrapped 32-byte DEK)
32      12   IV
44      16   auth tag
```

v0.7+ readers accept this layout indefinitely (binding decision 38);
the dispatch rule is **"if blob length == 60 and magic bytes are not
`'M','W'`, treat as v1.1"**.

### 5b.2 v1.2 versioned wrapped-key blob (v0.7, `wrapped_key_v2` flag)

When the `wrapped_key_v2` feature flag is present, `dek_wrapped` is a
variable-length blob:

```
offset  len  field
0       2    magic         = 0x4D 0x57  ('M','W' — TTIO Wrap)
2       1    version       = 0x02
3       2    algorithm_id  (big-endian)
               0x0000 = AES-256-GCM
               0x0001 = ML-KEM-1024 (v0.8 M49 — active; FIPS 203)
               0x0002 = reserved
5       4    ciphertext_len (big-endian, u32)
9       2    metadata_len   (big-endian, u16)
11      M    metadata       (algorithm-specific — see below)
11+M    C    ciphertext     (algorithm-specific — see below)
```

Algorithm-specific metadata/ciphertext layouts:

| `algorithm_id` | metadata                                        | ciphertext          | total blob |
|----------------|-------------------------------------------------|---------------------|-----------:|
| `0x0000` AES-GCM  | `iv(12) \|\| tag(16)` = 28                    | wrapped DEK = 32    | 71 bytes   |
| `0x0001` ML-KEM-1024 (v0.8) | `kem_ct(1568) \|\| aes_iv(12) \|\| aes_tag(16)` = 1596 | AES-GCM-wrapped DEK = 32 | 1639 bytes |

For ML-KEM-1024 the outer envelope contains a classical KEM + AEAD
construction: `ML_KEM.encapsulate(recipient_pk) → (kem_ct,
shared_secret)`, then `shared_secret` (32 bytes) is used as the AES-
256-GCM key to wrap the DEK. Decryption reverses the chain; AES-GCM
authenticates end-to-end (ML-KEM decapsulation on its own is
unauthenticated). See `docs/pqc.md` for the full story including the
language-specific library choices (Python/ObjC use liboqs;
Java uses Bouncy Castle).

Total length = `11 + metadata_len + ciphertext_len`. For AES-256-GCM
this equals 11 + 28 + 32 = **71 bytes**; for ML-KEM-1024 it equals
11 + 1596 + 32 = **1639 bytes**. The `@dek_wrapped_bytes`
attribute records the exact length so readers can avoid relying on
dataset-size probes through storage adapters that pad to a fixed
width.

Writers default to v1.2 when the feature flag is set; readers fall
back to v1.1 for any blob that doesn't start with the `'M','W'` magic.
`docs/feature-flags.md §v0.7` has the flag definition; `CipherSuite`
(M48) is the runtime dispatch catalog that maps `algorithm_id` to a
concrete cipher.

The DEK wraps signal data via AES-256-GCM (same as
`opt_dataset_encryption`). The KEK wraps the DEK. Rotation re-wraps
only the DEK — signal datasets are not touched, so rotation cost is
O(1) in file size.

---

## 6. Compound metadata datasets

All compound datasets live directly under `/study/`. Fields use HDF5
variable-length C strings (`H5Tvar_str`) where marked **VL**.

**Supported compound field kinds** (provider capability floor, all
three languages expose the same set):

| Kind        | HDF5 mapping                                      | Added   |
| ----------- | ------------------------------------------------- | ------- |
| `UINT32`    | `H5T_NATIVE_UINT32`                               | v0.3    |
| `INT64`     | `H5T_NATIVE_INT64`                                | v0.3    |
| `FLOAT64`   | `H5T_NATIVE_DOUBLE`                               | v0.3    |
| `VL_STRING` | `H5Tcopy(H5T_C_S1) + H5Tset_size(H5T_VARIABLE)`   | v0.3    |
| `VL_BYTES`  | `H5Tvlen_create(H5T_NATIVE_UCHAR)` (hvl_t slot)   | v0.10   |

`VL_BYTES` carries the {IV, tag, ciphertext} triplet of per-AU
encryption (§9.1). Providers that can't serialise `hvl_t` inside a
compound (SQLite + Zarr as of v0.10) raise
`NotImplementedError`/`UnsupportedOperationException` at the
`create_compound_dataset` boundary.

### 6.1 `identifications`

```
struct {
    VL string  run_name;
    uint32     spectrum_index;
    VL string  chemical_entity;
    float64    confidence_score;
    VL string  evidence_chain_json;  (JSON array of ref strings)
}
```

Feature flag: `compound_identifications`.

### 6.2 `quantifications`

```
struct {
    VL string  chemical_entity;
    VL string  sample_ref;
    float64    abundance;
    VL string  normalization_method;  (empty string = null)
}
```

Feature flag: `compound_quantifications`.

### 6.3 `provenance` (dataset-level)

```
struct {
    int64      timestamp_unix;
    VL string  software;
    VL string  parameters_json;      (JSON dict)
    VL string  input_refs_json;      (JSON array of refs)
    VL string  output_refs_json;     (JSON array of refs)
}
```

Feature flag: `compound_provenance`.

### 6.4 Per-run provenance

v0.3 (M17) migrates per-run provenance from the v0.2 `@provenance_json`
string attribute to a native compound HDF5 dataset at
`/study/ms_runs/<run>/provenance/steps`. The dataset uses the same
5-field compound type as the dataset-level `/study/provenance`
described in §6.3 above.

The v0.3 writer emits **both** forms during the transition window: the
compound dataset is the primary record, and the `@provenance_json`
legacy mirror is kept in place so the v0.2 signature manager (which
hashes the UTF-8 bytes of the JSON attribute) keeps working. A future
release will drop the legacy mirror once canonical-byte-order
signatures (§10 below) are wired to the compound dataset directly.

The v0.3 reader prefers the compound subgroup when present and falls
back to the `@provenance_json` attribute only when the subgroup is
absent; a run with neither form decodes to an empty provenance chain.

Feature flag: `compound_per_run_provenance`.

---

## 7. Image cube (MSImage)

`/study/image_cube/` is present when the dataset is an `TTIOMSImage`.

```
image_cube/
├── @width               (int64)
├── @height              (int64)
├── @spectral_points     (int64)
├── @tile_size           (int64)
├── @pixel_size_x        (float64)
├── @pixel_size_y        (float64)
├── @scan_pattern        (VL string)
└── intensity            (float64[H][W][SP])
```

The `intensity` dataset is rank-3, chunked with shape
`(tile_size, tile_size, spectral_points)`, and `zlib -6` compressed
behind the byte-shuffle filter (see §4).
This chunking ensures that reading a `(tileSize × tileSize)` tile
hits exactly one chunk.

Files written by v0.1 code have the equivalent layout at the root
group (`/image_cube/`) instead of `/study/image_cube/`. v0.2 readers
auto-detect both locations.

Feature flag (opt_): `opt_native_msimage_cube`.

---

## 7a. Vibrational imaging cubes (M73, v0.11)

`/study/raman_image_cube/` is present when the dataset carries an
`TTIORamanImage`, and `/study/ir_image_cube/` when it carries an
`TTIOIRImage`. Both groups share the MSImage cube layout but add
the modality-specific scalars and a shared wavenumber axis.

```
raman_image_cube/ (or ir_image_cube/)
├── @width                     (int64)
├── @height                    (int64)
├── @spectral_points           (int64)
├── @tile_size                 (int64)
├── @pixel_size_x              (float64)
├── @pixel_size_y              (float64)
├── @scan_pattern              (VL string)
├── @excitation_wavelength_nm  (float64)   ← raman_image_cube only
├── @laser_power_mw            (float64)   ← raman_image_cube only
├── @ir_mode                   (int64 enum) ← ir_image_cube only  (0 = TRANSMITTANCE, 1 = ABSORBANCE)
├── @resolution_cm_inv         (float64)   ← ir_image_cube only
├── intensity                  (float64[H][W][SP])
└── wavenumbers                (float64[SP])
```

`intensity` is chunked at `(tile_size, tile_size, spectral_points)`
with `zlib -6`, matching the MSImage convention so reading a
`(tileSize × tileSize)` tile hits exactly one chunk. `wavenumbers`
is a rank-1 companion that names the spectral axis and is
identical across all pixels — store once, not per-pixel.

The two groups are mutually exclusive per study; a file is either
a Raman map or an IR map, not both.

---

## 7b. UV-Vis spectra (M73.1, v0.11.1)

`TTIOUVVisSpectrum` is a plain `TTIOSpectrum` subclass — no new
group layout is required. It rides the generic `Spectrum`
persistence path, keyed by the following named signal arrays:

```
spec_NNNNNN/
├── @ttio_class = "TTIOUVVisSpectrum"
├── @path_length_cm   (float64, optional)
├── @solvent          (VL string, optional)
└── arrays/
    ├── wavelength    (float64[N], nm)
    └── absorbance    (float64[N])
```

No feature flag — pre-v0.11.1 readers open the group as a generic
`Spectrum` and retain the two channels unchanged.

---

## 7c. Two-dimensional correlation spectra (M73.1, v0.11.1)

`TTIOTwoDimensionalCorrelationSpectrum` is an `TTIOSpectrum`
subclass that carries a 1-D variable axis and two rank-2
correlation matrices of equal size.

```
spec_NNNNNN/
├── @ttio_class = "TTIOTwoDimensionalCorrelationSpectrum"
├── arrays/
│   ├── variable_axis  (float64[N])
│   ├── synchronous    (float64[N][N], row-major; in-phase, symmetric)
│   └── asynchronous   (float64[N][N], row-major; quadrature, antisymmetric)
```

Both matrices share the single variable axis — `nu_1 == nu_2`, so
no separate F1/F2 dimension scales are attached (differs from
§8's native 2-D NMR layout). Construction validates squareness
(`synchronous.shape == asynchronous.shape == (N, N)`).

Feature flag (opt_): `opt_native_2d_cos`. Pre-v0.11.1 readers
without the flag see the three arrays as opaque channels on a
generic `Spectrum` and round-trip them unchanged.

---

## 8. Native 2-D NMR

An `TTIONMR2DSpectrum` group (written via the generic
`TTIOSpectrum` persistence path under
`/study/nmr_runs/<run>/spec_NNNNNN/`) carries:

```
spec_NNNNNN/
├── @ttio_class = "TTIONMR2DSpectrum"
├── @matrix_width, @matrix_height
├── @nucleus_f1, @nucleus_f2
├── arrays/
│   └── intensity_matrix (1-D flattened; v0.1 fallback)
├── intensity_matrix_2d  (float64[H][W])           ← v0.2 addition
├── f1_scale             (float64[H])  H5T_SCALE
└── f2_scale             (float64[W])  H5T_SCALE
```

`intensity_matrix_2d` is the native rank-2 dataset, chunked at
`(min(128, H), min(128, W))` with `zlib -6`. `f1_scale` and
`f2_scale` are attached via `H5DSattach_scale` to dimensions 0 and 1
respectively so `h5dump` renders them as dimension labels.

Feature flag (opt_): `opt_native_2d_nmr`.

v0.2 readers prefer `intensity_matrix_2d` when present and fall back
to the flattened `intensity_matrix` in `arrays/` otherwise.

---

## 9. Encryption

`TTIOSpectralDataset.encryptWithKey:level:error:` performs two
operations in a single call:

1. **Per-run intensity channel encryption.** For each run under
   `/study/ms_runs/`, the plaintext `<channel>_values` dataset is
   read, encrypted with AES-256-GCM via OpenSSL, and replaced by
   `<channel>_values_encrypted` (int32-padded byte blob) alongside
   `<channel>_iv`, `<channel>_tag`, and the sizing attributes listed
   in §5.

2. **Compound dataset sealing.** If `/study/identifications` or
   `/study/quantifications` exist, each is read back into memory,
   serialized to JSON, encrypted with AES-256-GCM, and written as
   `<name>_sealed` (int32-padded bytes) plus `_iv`, `_tag`, and
   `_bytes` sibling datasets. The original compound dataset is
   deleted.

After encryption the root group carries:

- `@encrypted = "aes-256-gcm"`
- `@access_policy_json` (JSON representation of `TTIOAccessPolicy`)

`decryptWithKey:` reverses both operations.

Feature flag (opt_): `opt_dataset_encryption`. Required features may
also list `compound_identifications`/`compound_quantifications` even
after encryption, since those compound datasets will be restored on
decrypt.

### 9.1 Per-AU encryption (v1.0, `opt_per_au_encryption`)

v0.x encryption (above) encrypts an entire channel as a single
AES-GCM operation, which is incompatible with the per-Access-Unit
streaming transport introduced in v0.10. v1.0 adds a
per-spectrum encryption mode gated on the
`opt_per_au_encryption` feature flag.

New on-disk layout, one compound dataset per channel under
`/study/ms_runs/<run>/signal_channels/` — and, for genomic data,
under `/study/genomic_runs/<run>/signal_channels/`:

```
<channel>_segments          HDF5 compound, spectrum_count rows
    offset     INT64         index into the plaintext element stream
    length     INT64         plaintext element count
    iv         VL_BYTES      per-row AES-256-GCM IV (12 bytes)
    tag        VL_BYTES      per-row GCM tag (16 bytes)
    ciphertext VL_BYTES      ciphertext bytes (length-prefixed)
@<channel>_algorithm        "aes-256-gcm"
@<channel>_wrapped_dek      VL uint8
@<channel>_kek_algorithm    "rsa-oaep-sha256" | "ml-kem-1024" | …
@<channel>_wrapped_dek_recipients  (optional) base64 multi-recipient block; see §4.4 of `transport-spec.md` and `transport-encryption-design.md` §4.1.1
@<channel>_server_kek_id    (optional) UTF-8 server-resolvable KEK id; see `transport-encryption-design.md` §4.1.2
```

The IV / tag / ciphertext fields are HDF5 `VL_BYTES` (variable-
length byte arrays) for all three reference implementations. IV
is always 12 bytes, tag is always 16 bytes, and ciphertext varies
per row.

Each row's AES-GCM operation uses authenticated data
`dataset_id (u16 LE) || au_sequence (u32 LE) || channel_name_utf8`
so ciphertext cannot be replayed against a different spectrum
or channel. (`au_sequence` is the per-AU index; older drafts of
this doc called it `row_index` — the runtime AAD bytes are
identical, only the doc nomenclature is normalised.)

For a mixed-modality file (MS runs + genomic runs), the writer
processes MS first and then genomic, and `dataset_id` continues
incrementing through the genomic loop so a re-encrypt path that
walks runs in the same order produces byte-identical ciphertext.

When `opt_encrypted_au_headers` is also set, the plaintext
`spectrum_index/*` arrays from §4 are **omitted** and replaced
with:

```
spectrum_index/au_header_segments   HDF5 compound, one row per spectrum
    iv         VL_BYTES     12 bytes
    tag        VL_BYTES     16 bytes
    ciphertext VL_BYTES     36 bytes plaintext content, AES-GCM encrypted:
        acquisition_mode(u8) || ms_level(u8) || polarity(u8)
        || retention_time(f64) || precursor_mz(f64)
        || precursor_charge(u8) || ion_mobility(f64)
        || base_peak_intensity(f64)
@wrapped_dek                same DEK as the channel segments
```

AAD for the header row = `dataset_id (u16 LE) || au_sequence (u32 LE) || b"header"`.

The `decryptInPlace` family (`+[TTIOPerAUFile decryptFilePathInPlace:...]`,
`decrypt_per_au_in_place`, `PerAUFile.decryptFileInPlace`) reverses
the above: each row is GCM-decrypted, the plaintext is concatenated
back into the bare `<channel>_values` (MS) or `<channel>` (genomic)
dataset, `<channel>_segments` and `@<channel>_algorithm` are
removed, and on full success `opt_per_au_encryption` /
`opt_encrypted_au_headers` / root `@encrypted` are stripped from
the feature flags + root attributes.

No v0.x back-compat: files using the v0.x channel-grained
encryption layout cannot be opened by v1.0 encryption paths.
The legacy decryption code remains for migration; see the
`--transcode` flow described in
`docs/transport-encryption-design.md` §6.

#### 9.1.1 Genomic runs under `blocks_v1` (M99)

For a genomic run in the `blocks_v1` layout (§10.12), the walkers
stream **block by block** instead of materialising the channel:
each block's `sequences` / `qualities` blob is decoded with its
per-block codec from the block index, sliced per read using the
`genomic_index/lengths` slice, and encrypted as one AU per read
with **global** AU numbering (`au_sequence` = the read's index-row
ordinal; the stored segment `offset` is the read's global plaintext
offset). The segments compounds are extendable and appended one
block at a time, so peak memory follows the block policy, not the
channel size. The on-disk `<channel>_segments` layout, the AAD, and
the `dataset_id` ordering are identical to the whole-channel form
above — a reader of the segments tables cannot tell which walker
produced them.

Decrypt-in-place restores the **original layout**: each block's AUs
are decrypted, the block is re-encoded through the block writer
(replicating the stream writer's sticky qualities strategy: block 0
auto-tunes, the winner read back from the encoded stream pins the
rest), and the blobs are appended into recreated channel datasets.

So that the re-encode reproduces the original blobs, the stream
writers persist the writer policy that shapes them as run-group
attributes (M99.1), and the walkers set them on every reconstructed
block:

| Attribute | Type | Written when |
|---|---|---|
| `@ref_diff_slice_bytes` | uint64 | The REF_DIFF slice byte budget is non-zero. |
| `@opt_disable_qualities_v5` | int `1` | The qualities-V5 opt-out is set. |
| `@reference_md5s` | string, JSON object | The run's sequences code through `REF_DIFF_V2`. Maps each chromosome name of the writer's reference set to the hex **reference-set** md5 — the digest of every chromosome's bytes concatenated in sorted-name order, the same value the blob headers carry and the on-disk reference `@md5` records. |

A `REF_DIFF_V2` run's restore rebuilds `reference_chrom_seqs`
through the reference resolver from `@reference_md5s`: the embedded
`/study/references/` copy when present, else a `REF_PATH` FASTA
verified against the digests. Files without the attribute keep the
embedded-directory walk; when neither source resolves, the walkers
refuse with an error naming both.

When every re-encoded blob lands on the ranges the block index
records — the case whenever the persisted policy is honoured — the
index is untouched and the restore is **byte-identical**. When any
blob differs (an older file without the policy attrs), restore does
not refuse: it appends the blobs it produced and rewrites
`blocks/index` with the offsets, lengths and codec ids actually
written, keeping the file consistent and readable at the cost of
byte-identity.

The encrypted transport stream carries `blocks_v1` genomic runs via
the transport-spec v0.12 sidecar packets (`GenomicRunSidecar` 0x1C,
`BlockSidecar` 0x1D, feature token `transport_blocks_v1` — see
`transport-spec.md` §4.24): the AU stream carries the encrypted
sequences and qualities per read, the sidecars carry the run
scalars, the restore attributes, both chromosome name tables, the
block index rows and the verbatim plaintext
`read_names`/`cigars`/`mate_info` blob slices, so decrypt-in-place
on the received container restores it byte-identically. The stream
does not carry embedded reference bytes; a received `REF_DIFF_V2`
run restores through `@reference_md5s` and `REF_PATH`.

---

## 10. Digital signatures

`TTIOSignatureManager.signDataset:inFile:withKey:error:` computes an
HMAC-SHA256 over a canonical byte stream derived from the target
dataset and stores a **prefixed** base64 MAC in `@ttio_signature`.
Provenance chain signing stores its MAC in `@provenance_signature`
on the run group, computed over the UTF-8 bytes of
`@provenance_json` (legacy v0.2 path, kept for v0.2 compatibility).

### 10.1 v2 canonical signatures (M18, v0.3 default)

Stored as `"v2:" + base64(mac)`. The canonical byte stream is:

- **Atomic numeric datasets** (float / int / uint, 1–8 bytes) — read
  via an explicit little-endian HDF5 memory type
  (`H5T_IEEE_F64LE`, `H5T_STD_U32LE`, ...). The resulting byte buffer
  is canonical on any host architecture.
- **Compound datasets** — each record is walked in declaration order.
  Numeric members are read via a packed memory type that maps each
  atomic member to its LE equivalent. Variable-length string members
  are emitted as `u32_le(byte_length) || utf8_bytes`, so struct
  padding and pointer layouts cannot influence the hash.
- **Other classes** (fixed strings, enums, nested compounds) — fall
  back to the native-bytes path, matching v0.2 behaviour.

The v0.3 signer adds both `opt_digital_signatures` and
`opt_canonical_signatures` to the root feature list on first sign.

### 10.2 v1 native-byte signatures (v0.2 compatibility)

Signatures without the `v2:` prefix are treated as v0.2 native-byte
HMACs and verified by hashing `H5Dread` output in the dataset's
native type. This is how the v0.2 `signed.tio` reference fixture
continues to verify under v0.3 readers.

### 10.2b v3 post-quantum signatures (v0.8 M49)

ML-DSA-87 (FIPS 204) signatures are stored as:

```
"v3:" + base64(ml_dsa_87_signature_bytes)
```

The signature covers the same canonical little-endian byte stream
as v2, so a file can carry either flavor without format changes
beyond the prefix. A reader that sees a `v3:` prefix but does not
support PQC raises `UnsupportedAlgorithmError` — it does not
silently pass verification. The presence of any v3 signature on a
file causes `opt_pqc_preview` to be added to the root feature
list; see `docs/pqc.md` and `docs/feature-flags.md §v0.8` for the
full story.

### 10.3 Cross-language parity

The TTIO Objective-C and Python implementations produce byte-identical
`v2:` MACs for the same input (see the `TtioSign` CLI test harness
under `objc/Tools/` and `python/tests/test_canonical_signatures.py`).

Feature flags: `opt_digital_signatures` (first sign), plus
`opt_canonical_signatures` when any `v2:` signature is present.

---

## 10c. Byte-level protocol contract (M43, v0.7)

All cryptographic paths — signatures and dataset / envelope encryption
— consume their input through the `StorageDataset.read_canonical_bytes`
method defined by the protocol abstraction
(`TTIOStorageDataset`, `global.thalion.ttio.providers.StorageDataset`,
`ttio.providers.base.StorageDataset`). The canonical stream is:

- **Primitive numeric datasets** — little-endian packed values.
- **Compound datasets** — rows in storage order; fields in declaration
  order. Variable-length strings encoded as `u32_le(length) ||
  utf-8_bytes`. Numeric fields little-endian.

On big-endian hosts the conversion is an explicit byteswap (HDF5's
automatic type conversion is **not** relied on). Every provider that
ships with v0.7 (`Hdf5Provider`, `MemoryProvider`, `SqliteProvider`,
`ZarrProvider`) emits bit-equal bytes for the same logical data;
cross-backend round-trip tests in
`python/tests/test_canonical_bytes_cross_backend.py` and
`python/tests/test_zarr_provider.py::test_compound_canonical_bytes_matches_hdf5`
lock that guarantee.

Binding decision 37: a signed or encrypted dataset verifies
identically regardless of which provider wrote it.

---

## 10.4 Compression codecs

Signal-channel datasets carry their compression codec via either the
HDF5 filter pipeline (codec ids 1–3) or a dedicated per-channel
`@compression` uint8 attribute (codec ids 4+):

| Id | Codec                  | Transport                                                                                                                                                                     |
|---:|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 0  | NONE                   | Passthrough.                                                                                                                                                                  |
| 1  | ZLIB                   | `H5P_DEFLATE` filter at level 6 (default for non-genomic channels). Lossless. Readable by any HDF5 library without extra plugins.                                              |
| 2  | LZ4                    | HDF5 filter id **32004**. Requires the LZ4 filter plugin (`libh5lz4.so`) to be loadable at runtime via `HDF5_PLUGIN_PATH`. Lossless. ~35× faster write / ~2× faster read than zlib, at ~20% larger files on random data. |
| 3  | NUMPRESS_DELTA         | Per-channel transform implemented inside TTIO, **not** an HDF5 filter. The dataset stores an `int64` array of first differences of a fixed-point quantised signal. The signal_channels group carries `@<channel>_numpress_fixed_point` (int64) giving the scaling factor. Readers detect the codec via that attribute. Lossy, sub-ppm relative error for typical mass-spectrometry m/z. Clean-room implementation from Teleman et al., *MCP* 13(6), 2014. **Not a size optimisation** — on the PXD000001 Orbitrap corpus, MS1 profile m/z stored as numpress deltas + zlib is 119.4 MB where a modern lossless numerical codec stores the identical channel bit-exact in 61.9 MB (2026-08 compression audit). Keep it for mzML / msNumpress interchange parity; do not choose it to shrink files. |
| 4  | RANS_ORDER0            | Range Asymmetric Numeral Systems entropy coder, order-0 (per-byte) frequency model. Clean-room from Duda 2014, no htslib source consulted. Wire format + deterministic frequency-table normalisation specified in `docs/codecs/rans.md`. |
| 5  | RANS_ORDER1            | Order-1 (preceding-byte context) rANS variant. Same wire format as order-0; per-context frequency tables (256 of them) with run-length-encoded sparse rows. See `docs/codecs/rans.md`. |
| 6  | BASE_PACK              | 2-bit ACGT packed-base codec for genomic sequences. Lossless on the full byte alphabet via a sparse position+byte sidecar mask: bases that are uppercase `{A,C,G,T}` pack into 2-bit slots (4 bases per output byte, big-endian within byte); everything else (`N`, IUPAC ambiguity codes, lowercase soft-masking, gaps) is recorded in the mask alongside its input position so the decoder restores it byte-for-byte. See `docs/codecs/base_pack.md`. |
| 7  | QUALITY_BINNED         | Illumina-8 / CRUMBLE-derived 8-bin Phred quantisation. 4-bit-packed bin indices, big-endian within byte. Lossy (decode returns the bin centre, not the original Phred). See `docs/codecs/quality.md`. |
| 8  | _RESERVED_8            | Reserved on the wire (Java enum ordinal stability). Was the v1 NAME_TOKENIZED codec; superseded by id 15. Reader paths reject id 8 with a v1.0 migration error. |
| 9  | _RESERVED_9            | Reserved on the wire. Was the v1 REF_DIFF codec; superseded by id 14. |
| 10 | _RESERVED_10           | Reserved on the wire. Was the FQZCOMP_NX16 v1 codec; superseded by id 12. |
| 11 | DELTA_RANS_ORDER0      | Running delta + zigzag-varint + rANS-O0 wrapper for sorted-ascending integer arrays. Read-compatibility codec: no current write path selects it (per-read integer fields moved to `genomic_index/` with the zlib filter in v1.6); readers decode it from files written under the earlier integer-channel defaults, where `positions` carried it by default. See `docs/codecs/delta_rans.md`. |
| 12 | FQZCOMP_NX16_Z         | CRAM-mimic lossless quality codec. Magic `M94Z`, V4 wire format. Static-per-block frequency tables (zlib-deflated in the header so the decoder skips the build pass), 16-bit renormalisation (`B = 16`, `b = 2^16`), `T = 4096` fixed (12-bit shift); `T \| b·L = 2^31` exactly, byte-pairing mathematically guaranteed. Bit-pack context model: 12 bits `prev_q` + 2 bits position bucket + 1 bit revcomp, masked to `2^14` slots. 4-way interleaved rANS states. **Default for the `qualities` channel.** Since qualities V5, the encoder also tries the S5/S6 sequence-context strategies when the run carries a base-parallel `sequences` channel and keeps the smallest stream by exact size; a winning sequence strategy is emitted as a `M94Z` version-5 stream whose body is the explicit 8-byte parameter block plus a range-coded stream, and whose decoder requires the decoded sequences as side input. Files where V4 wins stay byte-identical version-4 streams. Writers opt out per run with `opt_disable_qualities_v5` / `optDisableQualitiesV5`. Clean-room implementation of CRAM 3.1 `rANS-Nx16` discipline; no htslib / tools-Java source consulted. See `docs/codecs/fqzcomp_nx16_z.md`. |
| 13 | MATE_INLINE_V2         | Inlined per-record mate_info (chrom + pos + tlen) as a single channel. Replaces the M82 compound + per-field subgroup decomposition. **Default for the `mate_info` channel.** See `docs/codecs/mate_info_v2.md`. |
| 14 | REF_DIFF_V2            | **Context-aware** reference-based sequence-diff codec. Encoder/decoder consume sibling channels (`positions`, `cigars`) and an external reference resolver alongside the channel bytes. Slice-based wire format with embedded reference at `/study/references/<reference_uri>/`. **Default for the `sequences` channel** when a reference is available; falls back to BASE_PACK silently when not. See `docs/codecs/ref_diff_v2.md`. |
| 15 | NAME_TOKENIZED_V2      | 8-substream multi-token columnar codec for read names. Substreams: FLAG / POOL_IDX / MATCH_K / COL_TYPES / NUM_DELTA / DICT_CODE / DICT_LIT / VERB_LIT, each auto-picked between rANS-O0 and raw passthrough. Per-block reset every 4096 reads; magic `NTK2`. **Default for the `read_names` channel.** See `docs/codecs/name_tokenizer_v2.md`. |
| 16 | ZSTD                   | Zstandard (RFC 8878). Wire-only: an opt-in codec for spectral access-unit channels on the transport stream (`transport-spec.md` §4.3). No on-disk `@compression` dispatch. |
| 17 | FLOAT_DELTA_ZSTD       | Lossless float64 channel codec: per block of 2^20 values, none/delta on the uint64 bit view (chosen by exact size comparison), byte-plane transpose, one zstd frame. Magic `FDZ1`. Values round-trip bit-exactly (NaN payloads, signed zeros, Inf). The default for float64 channels of `TTIOMassSpectrum` runs (writers opt out via `opt_disable_float_delta` / `optDisableFloatDelta`); other spectral classes opt in via `signal_compression="float_delta_zstd"`. The dataset becomes a flat uint8 stream with `@compression = 17` and no HDF5 filter. Encoders MAY differ byte-wise across languages (zstd builds differ); decoders MUST accept any conforming stream — a shared golden fixture pins the decode side. See `docs/codecs/float_delta_zstd.md`. |

| 18 | SAM_TAGS               | SAM optional fields (columns 12+) of a run of reads (M101): a per-blob tag-line dictionary, one column per (tag key, type), MD:Z and NM:i recomputed from the reference when they match, integers equal to an earlier tag in the read stored as back-references, verbatim storage for non-canonical text. Lossless on the samtools tag text. Magic `STG1`. Context-aware (sequences, CIGARs, positions, chromosome ids, reference). One native kernel shared by the three SDKs. See §10.13 and `docs/codecs/sam_tags.md`. |
| 19 | SEQ_CM                 | Read bases without a reference (M103): up to three base-context lengths (default 11, 16, 24) in adaptive-counter tables, a logistic mixer, training on each read's reverse complement, and a carry-less binary arithmetic coder; bytes other than A/C/G/T round-trip as exception runs. All arithmetic is integer, so every platform codes the same bytes. Magic `SQC1`. Context-aware: decode takes the reads' lengths. **Default for the `sequences` channel under `blocks_v1` when there is no reference.** One native kernel shared by the three SDKs. See `docs/codecs/seq_cm.md`. |

Ids `0`–`3` ride the HDF5 filter pipeline; ids `4`+ are signalled via
the per-channel `@compression` attribute (see §10.5). Reserved ids
`8` / `9` / `10` retain their wire-format slots so cross-language
readers can defensively reject removed-codec files with a clear v1.0
migration error.

### Codec-channel applicability

- Ids `4`, `5` (rANS order-0/1) apply to every byte and integer
  channel via the LE byte-serialisation contract in §10.7.
- Id `6` (BASE_PACK) applies to the `sequences` channel.
- Id `7` (QUALITY_BINNED) applies to the `qualities` channel only —
  validation rejects QUALITY_BINNED on `sequences` because Phred-bin
  quantisation would silently destroy ACGT data.
- Id `11` (DELTA_RANS_ORDER0) applies to sortable integer channels
  (`positions` is the canonical case).
- Id `12` (FQZCOMP_NX16_Z) applies to the `qualities` channel and
  is the v1.0 default. Its M94.Z stream is version 4, or version 5
  (sequence-context body) when the writer had base-parallel
  sequences and the sequence strategy won by exact size; version-5
  streams decode against the run's decoded sequences channel.
- Id `13` (MATE_INLINE_V2) applies to the `mate_info` channel and
  is the v1.0 default.
- Id `14` (REF_DIFF_V2) applies to the `sequences` channel of
  reference-aligned `WrittenGenomicRun` instances and is the v1.0
  default when a reference is available.
- Id `15` (NAME_TOKENIZED_V2) applies to the `read_names` channel
  and is the v1.0 default.
- Id `17` (FLOAT_DELTA_ZSTD) applies to spectral float64 signal
  channels (`mz`, `intensity`, `chemical_shift`, FID channels) and
  is the default for `TTIOMassSpectrum` runs (Phase 2 of the codec
  spec; `opt_disable_float_delta` preserves the previous layout).
  Non-MS spectral channels keep the HDF5 shuffle + zlib filter
  pipeline as their default and opt in explicitly.
- Id `18` (SAM_TAGS) applies to the genomic `tags` channel only, and
  is its only codec (§10.13).
- Id `19` (SEQ_CM) applies to the `sequences` channel and is the
  `blocks_v1` default when a block has no reference; writers accept it
  as a `sequences` override on either layout.

See §10.5 for the `@compression` attribute scheme, §10.6 for the
`read_names` channel format, §10.7 for the integer-channel
serialisation contract, §10.8 for the cigars channel, §10.9 for
the mate_info channel, §10.9b for the MATE_INLINE_V2 wire format,
§10.10b for the NAME_TOKENIZED_V2 wire format.

The genomic codec pipeline-wiring covers ALL M82 channels — every
channel under `signal_channels/` has at least one applicable codec
with cross-language byte-exact conformance.

## 10.5 `@compression` attribute on signal-channel datasets (M86)

Genomic signal-channel datasets (`signal_channels/sequences` and
`signal_channels/qualities` under `/study/genomic_runs/<name>/`)
that use a TTI-O internal compression codec (rANS order-0, rANS
order-1, BASE_PACK) carry a `@compression` attribute holding the
M79 codec id. The attribute type is `H5T_NATIVE_UINT8` (one byte).
The dataset bytes ARE the self-contained codec stream specified in
`docs/codecs/rans.md` (ids `4`, `5`) or `docs/codecs/base_pack.md`
(id `6`). **No HDF5 filter is applied to such datasets** — the
codec output is high-entropy and would not benefit from deflate.

Absence of the attribute, or value `0` (`Compression.NONE`), means
the dataset is stored as-is and any HDF5 filter applies (typically
zlib level 6, the default for genomic byte channels). The attribute
is written ONLY when an override is in effect; uncompressed channels
have no `@compression` attribute at all.

Pre-M86 readers that ignore `@compression` will silently
misinterpret a v0.12-encoded channel: the read path slices into a
non-sliceable codec stream and returns garbage for any read whose
offset/length walks past the encoded payload boundary. The
attribute is the canonical signal for codec dispatch — any
TTI-O-conformant reader from M86 onwards must check it before
slicing.

M86 wires this attribute scheme for the **byte channels only**:
`sequences` and `qualities`. Integer channels (`positions`,
`flags`, `mapping_qualities`) and VL_STRING channels (`cigars`,
`read_names`, `mate_info`) do not yet support TTI-O codecs; they
ignore `@compression` if set and stay on HDF5-filter ZLIB. Lifting
that restriction (integer-channel codecs, plus M85's
name-tokenizer for read_names) is a future milestone.

### Read-side dispatch (informative)

When opening a `sequences` or `qualities` dataset, an M86 reader:

1. Checks for the `@compression` attribute. If absent or `0`, uses
   the existing slice-based read path (no change from M82).
2. If `4`, `5`, or `6`, reads ALL dataset bytes, decodes the whole
   stream through the corresponding `decode()` function, and
   caches the decoded buffer on the open `GenomicRun` instance.
   Subsequent per-read access slices the cached buffer in memory.

This decode-once-cache strategy is the natural shape for the
M83/M84 codecs, which produce non-sliceable byte streams. The
memory cost is one decoded channel per open run instance —
acceptable for typical sequencing workloads.

## 10.6 `read_names` channel layout

`signal_channels/read_names` is a flat 1-D `UINT8` dataset with
`@compression = 15` (NAME_TOKENIZED_V2). Dataset bytes are the v2
stream specified in §10.6b. **No HDF5 filter applied.**

A v1.0 reader rejects any other layout (including the M82
compound `{value: VL_STRING}` and the v1 NAME_TOKENIZED
`@compression = 8` flat-byte layout) with a v1.0 migration error.

Under `blocks_v1` (section 10.12) this dataset holds one such blob per block, back to back, addressed through `blocks/index`.

## 10.6b `NAME_TOKENIZED_V2` wire format (codec id 15)

NAME_TOKENIZED_V2 is a multi-substream + DUP-pool + PREFIX-MATCH
codec for read names. Wire magic `NTK2`, container version
`0x01`. Default for the `read_names` channel. §4 of the design spec
has the authoritative byte-level layout.

Summary:

```
Container header (12 + 4·n_blocks bytes):
  4 bytes "NTK2" magic
  1 byte version = 0x01
  1 byte flags (bit 0 = empty stream)
  4 bytes n_reads u32 LE
  2 bytes n_blocks u16 LE (≤ 65535)
  n_blocks × 4 bytes block_offset[i] u32 LE
    (offset of block i body relative to start of first block)

Per block (≤ 4096 reads):
  4 bytes block_n_reads u32 LE
  4 bytes block_body_len u32 LE
  block body: 8 substreams in fixed order, each prefixed by:
    4 bytes substream_body_len u32 LE
    1 byte mode (0x00 raw, 0x01 rANS-O0)
    body bytes
```

Substream order: FLAG (2-bit per read), POOL_IDX (3-bit per
DUP/MATCH row), MATCH_K (varint per MATCH row), COL_TYPES
(per-block column-type bitmap), NUM_DELTA (numeric column
deltas), DICT_CODE (string column codes), DICT_LIT (string
column literals), VERB_LIT (verbatim escape rows). All
emissions row-major within NUM_DELTA / DICT_CODE.

Reader auto-picks the per-substream mode (0x00 raw vs 0x01
rANS-O0) and dispatches accordingly. Block boundaries are
independently decodable; the block-offset table enables O(1)
seek to any block's body.

**Pre-v1.9 reader behaviour:** Files with `@compression = 15`
are unreadable by pre-v1.9 readers — they raise "unknown codec
id" at read time. v1.x → v1.9 forward-compat is read-only via
the opt-out flag (writer downgrades to M82 compound).

**Pre-M86-Phase-E reader behaviour:** A v0.12 file with the
override is **unreadable** by pre-M86 readers — they expect
the M82 compound layout and will silently misinterpret the
flat-uint8 dataset as a corrupt compound. Discipline matches
M80 / M82 / M86 Phase A (write-forward, no back-compat shim).
Files written without the override remain identical to M82
output and are read identically by all reader versions.

The other VL_STRING channels (`cigars`, `mate_info`) do NOT
currently support a codec override; they remain in compound
storage. `cigars` would want an RLE-then-rANS pipeline (no
codec match in M79); `mate_info` is an integer-tuple compound
with no codec match.

## 10.7 Integer-channel serialisation contract

The per-record integer channels (`positions` int64, `flags` uint32,
`mapping_qualities` uint8) are stored exclusively under
`genomic_index/`, mirroring the MS modality's `spectrum_index/`
pattern: per-record metadata is eagerly loaded; `signal_channels/`
holds bulk per-base / variable-length data only.

`signal_codec_overrides[positions|flags|mapping_qualities]` is
rejected at write time with a `ValueError` (Python) /
`IllegalArgumentException` (Java) / `NSInvalidArgumentException`
(Objective-C) — these channels have no `signal_channels/` presence
to override.

When integer arrays appear inside other codec wire formats (e.g.
the MATE_INLINE_V2 NS / NP / TS substreams; the cigars channel's
length-prefix-concat stream), they are serialised in **little-endian
byte order**:

- Python: numpy dtype strings `<i8`, `<u4`, `<u1`.
- Objective-C: `OSSwapHostToLittleInt64` / `htole32` etc.
- Java: `ByteBuffer.LITTLE_ENDIAN` with `putLong` / `putInt` / `put`.

LE is fixed and non-negotiable across all three implementations;
codec wire formats that embed integer values rely on this contract.

## 10.8 `cigars` channel codec wiring

The `cigars` channel under `signal_channels/` accepts the rANS
entropy coders (codec ids `4` = RANS_ORDER0 and `5` = RANS_ORDER1)
via `@compression`. The default is `RANS_ORDER1`.

### 10.8.1 On-disk schema

Flat 1-D `UINT8` dataset of length = encoded byte count, with
`@compression` set to `4` or `5`. **No HDF5 filter** is applied.
A `@compression` value outside `{4, 5}` is malformed.

Under `blocks_v1` (section 10.12) this dataset holds one such blob per block, back to back, addressed through `blocks/index`.

### 10.8.2 The rANS-on-cigars serialisation contract

The encoder serialises the `list[str]` of CIGARs to a flat byte
stream using length-prefix concatenation:

```
For each cigar in cigars:
    emit varint(len(cigar.encode('ascii')))
    emit cigar.encode('ascii')
```

Varints are unsigned LEB128 (low 7 bits + continuation flag). The
serialised buffer is then encoded via `Rans.encode(buf, order)`.
The decoder reverses: `Rans.decode(stream)` → byte buffer → walk
varint-length-prefix entries until exhausted.

CIGAR strings are 7-bit ASCII per the SAM spec; encoders reject
non-ASCII input.

### 10.8.3 Empirical compression

Measured byte-identical across Python / ObjC / Java on a 1000-read
mixed-CIGAR corpus:

- M82 compound (no override): ~18–29 KB depending on HDF5 filter.
- `RANS_ORDER1`: **1111 bytes** (~17× smaller than the M82 baseline).

rANS exploits byte-level repetition over the limited CIGAR
alphabet (digits 0–9 + the operator letters `MIDNSHP=X`), so the
codec is robust across uniform-CIGAR and mixed-CIGAR inputs.

## 10.9 `mate_info` channel layout

The `mate_info` channel is encoded via MATE_INLINE_V2 (codec id
13) — see §10.9b for the wire format.

The per-field per-record `signal_codec_overrides[mate_info_chrom |
mate_info_pos | mate_info_tlen]` keys are rejected at write time
with a `ValueError` pointing at this section: under v1.0 the v1
M86 Phase F per-field subgroup decomposition is gone, and the bare
key `"mate_info"` is also rejected (use the inline v2 codec
default, or write a custom decode path against the on-disk inline
blob).

## 10.9b mate_info v2 inline codec (codec id 13)

Encodes the full mate triple (mate_chrom_id, mate_pos, tlen) as a
single CRAM-style inline blob exploiting SAM mate-pair invariants.
Saves ~47.9 MB on chr22 vs the M82 compound baseline (see
`docs/benchmarks/2026-05-03-mate-info-v2-results.md`).

### 10.9b.1 On-disk schema

Two sibling datasets under `signal_channels/mate_info/`:

```
signal_channels/mate_info/
├── inline_v2     uint8 1-D blob, @compression = 13 (MATE_INLINE_V2)
└── chrom_names   compound[(name, VL_STRING)], one row per chrom_id
```

The `inline_v2` blob carries a 34-byte container header + 4 substreams
(MF / NS / NP / TS); full wire format in the design spec §4.

The `chrom_names` sidecar is a compound dataset that maps chrom_ids
(row index) to chromosome names. **This is necessary because mate
chromosomes can reference chroms that no own-read uses** (e.g. a
properly aligned read on chr22 with a cross-chrom mate on chr11
where no other read aligns to chr11). The L1 `genomic_index/chromosome_names`
table only covers chroms that appear as own_chrom; mate-only chroms
would be lost without `chrom_names`.

The chrom_id assignment is encounter-order over `(own_chromosomes ∪
mate_chromosomes)`, with own chroms enumerated first. The `'='`
SAM shortcut is canonicalised at write time to the record's own
chrom_id; `'*'` maps to -1.

The 4-substream container wire format:

| Substream | Content | Encoding |
|-----------|---------|----------|
| MF | Per-record mate-flag (0=SAME_CHROM, 1=CROSS_CHROM, 2=NO_MATE, 3=RESERVED — never emitted) | raw-pack or rANS-O0 (auto-pick) |
| NS | Chrom_id for CROSS_CHROM records (0 elsewhere) | varint + rANS-O0 auto-pick |
| NP | Mate_pos: zigzag-varint delta vs read position for SAME_CHROM, absolute for CROSS_CHROM; 0 for others | varint + rANS-O0 auto-pick |
| TS | Template_length (zigzag-varint); 0 for NO_MATE | varint + rANS-O0 auto-pick |

Container header (34 bytes): 4-byte magic `b"MIv2"` + 1-byte version
`\x01` + 1-byte flags `\x00` + `n_records` (uint64 LE) + `num_cross`
(uint32 LE) + 4 × uint32 LE substream byte lengths (MF, NS, NP, TS).
See [docs/codecs/mate_info_v2.md](codecs/mate_info_v2.md) for the full codec.

Under `blocks_v1` (section 10.12) this dataset holds one such blob per block, back to back, addressed through `blocks/index`.

### 10.9b.2 Reader-side dependency

Decoding `inline_v2` requires `genomic_index/positions` and
`genomic_index/chromosome_ids` to be loaded first. The decoder needs
own_pos and own_chrom_id per record to reconstruct mate_pos for
SAME_CHROM records (which use a delta encoding) and to validate the
MF taxonomy.

Readers must enforce this read order; the v1.7+ Python/Java/ObjC
implementations do so transparently.

### 10.9b.3 Backward compatibility

A v1.6 reader on a v1.7 file fails with "unknown compression id 13"
when it encounters the `inline_v2` dataset. The user must upgrade
the reader OR write the source file with
`opt_disable_inline_mate_info_v2 = True` to keep the v1 layout
from §10.9.

A v1.7 reader on a v1.6 file finds no `inline_v2` dataset and falls
through to the v1 layout transparently — the reader dispatches on
whether the subgroup contains `inline_v2` before checking for the
v1 per-field children.

### 10.9b.4 signal_codec_overrides interaction

Setting `signal_codec_overrides[mate_info_chrom / mate_info_pos /
mate_info_tlen]` when `opt_disable_inline_mate_info_v2 == False`
raises a write-time error pointing at the opt-out flag. The v1
per-field codec dispatch from §10.9 is only available under the
opt-out path.

### 10.9b.5 Cross-language byte-exact

The encode/decode primitives are implemented as a shared C kernel
in libttio_rans (`ttio_mate_info_v2_encode` / `ttio_mate_info_v2_decode`
entry points in `native/src/mate_info_v2.{c,h}`). All three language
implementations (Python ctypes, Java JNI, ObjC direct link) call
the same C functions, so the encoded byte stream is byte-exact
identical regardless of which language wrote the file. Verified at
test time by `python/tests/integration/test_mate_info_v2_cross_language.py`
(4 corpora × 3 languages = 12 byte-exact assertions, all PASS).

### Precision additions (M79, v0.11)

`TTIOPrecision` gains `UINT8` (id `6`) for byte-typed datasets —
genomic packed-base buffers, quality-score arrays, and any future
per-element symbol stream that does not need wider integers. The
existing storage providers (HDF5, Memory, SQLite, Zarr) honour
`UINT8` byte-exactly; canonical bytes for a `UINT8` dataset are the
raw payload (endian-neutral).

### Numpress-delta algorithm

1. Compute scale `S = floor((2^62 - 1) / max|v|)`; degenerate ranges
   default to `S = 1`.
2. Quantise: `q[i] = llround(v[i] * S)` (IEEE-754 round-to-even).
3. Emit `deltas[0] = q[0]`, `deltas[i] = q[i] - q[i-1]` for `i ≥ 1`.
4. Store `deltas` as the `<channel>_values` int64 HDF5 dataset with
   zlib on top.

Decoding is the exact inverse: cumsum the int64 array, cast to
double, divide by the scale. The TTIO ObjC and Python encoders agree
byte-for-byte on any input (see `test_numpress_scale_matches_objc_formula`).

The codec's role is mzML / msNumpress interchange parity, not file
size: the 2026-08 compression audit measured numpress + zlib at
119.4 MB on PXD000001 MS1 profile m/z against 61.9 MB for the same
channel stored losslessly by a modern numerical codec. Accepting
sub-ppm loss to end up with a file nearly twice the lossless size is
the wrong trade everywhere except when byte-level msNumpress
compatibility is itself the requirement.

## 10.10 `sequences` channel reference storage

REF_DIFF_V2 (codec id 14) is **context-aware** — encoder and
decoder consume sibling channels (`positions`, `cigars`) and an
external reference resolver alongside the channel bytes. This
section documents the on-disk layout for the embedded reference;
the codec wire format is in §10.10b.

### Reference storage group `/study/references/<reference_uri>/`

When `WrittenGenomicRun.embed_reference == True` (the default) and
the run uses REF_DIFF_V2 on its `sequences` channel, the writer
embeds the covered chromosome sequences at:

```
/study/references/<reference_uri>/
  @md5             : fixed-length string — 32-char hex of md5(concat(
                     sorted_chromosomes_in_uri_order))
  @reference_uri   : fixed-length string — the URI itself
  chromosomes/
    <chrom_name>/
      @length      : int64 — chromosome length in bases
      data_packed  : 1-D uint8 dataset holding the packed stream
                     (zlib-compressed) — written when packing wins
      data         : 1-D uint8 dataset of raw sequence bytes
                     (zlib-compressed) — the fallback layout, and the
                     only one pre-change readers understand
```

Exactly one of `data_packed` / `data` is present per chromosome.
The packed stream (big-endian, version byte `0x01`) is a 2-bit ACGT
body plus a run mask: `version(u8) || original_length(u32) ||
run_count(u32) || run_count × (position(u32) || length(u32)) ||
run_bytes || packed_body`. Runs are maximal stretches of non-ACGT
bytes (N blocks, IUPAC codes, soft-masked lowercase), recorded with
their original bytes so decode is byte-exact; the body packs the
remaining ACGT bytes 4-per-byte, first base in the two
highest-order bits. Writers pack only when the sequence is ≥ 50%
uppercase ACGT AND the packed stream is smaller than the raw bytes
— the decision is deterministic on content, so all three reference
implementations choose the same layout. On chr22 the packed layout
lands at 8.24 MB against 9.71 MB for raw + zlib. Readers probe
`data_packed` first and fall back to `data`; pre-change readers
fail with a missing-`data` error on packed chromosomes rather than
misreading bytes.

**Auto-deduplication:** when multiple runs in the same `.tio` file
share a `reference_uri`, the writer embeds the reference at this
path **once**. Subsequent runs that name the same URI link by
reference; their `@reference_uri` attribute matches the embedded
group's URI. A second run with the same URI but a different MD5
raises `ValueError` at write time — the URI–MD5 binding is
one-to-one within a file.

### Single canonical layout for all writers

Both the REF_DIFF_V2 auto-embed path (triggered when
`WrittenGenomicRun.embed_reference == True` and a run uses
REF_DIFF_V2) and the explicit `ReferenceImport.write_to_dataset()`
path used by the `FastaReader` / `FastaWriter` pair (Python
`ttio.importers.fasta` / `ttio.exporters.fasta`, Java
`global.thalion.ttio.importers.FastaReader` /
`...exporters.FastaWriter`, ObjC `TTIOFastaReader` /
`TTIOFastaWriter`) emit the **same** group shape documented above:
3-level `<uri>/chromosomes/<chrom_name>/data`, with
`@reference_uri` and `@md5` on the URI group and `@length` on each
per-chromosome group. The `@md5` attribute is computed by both
paths using the same single canonical seq-only formula (see "MD5
computation" below); byte-identical digests are guaranteed for the
same content.

#### Historical note (pre-v1.1.1)

Prior to v1.1.1, `ReferenceImport.write_to_dataset()` emitted a
distinct **flat-layout variant** with attributes `@md5` and
`@total_bases` on the URI group and flat `<group>/<chromName>`
datasets (no nested `chromosomes/` sub-group, no `@reference_uri`,
no per-chromosome `@length`). The chromosome `data` datasets were
gzip-compressed and **case-preserving** (lowercase soft-masking
bytes survived the round-trip), in contrast with the auto-embed
path's uppercase-normalised ACGTN bytes. The flat layout was
unified with the canonical 3-level form in v1.1.1 (commits
`ab70a27` + `586d6bd`), making `write_to_dataset()` round-trippable
through `SpectralDataset.references()`.

Files written with the pre-1.1.1 flat layout retain their on-disk
`@md5` verbatim (the canonical seq-only formula was already in
effect from v1.1.0), and `read_from_group`'s fallback path reads
their chromosome bytes directly. However, those files are **not**
exposed through `SpectralDataset.references()` — readers needing
the new accessor must rewrite the file via the v1.1.1+ writer. The
`@total_bases` attribute and case-preserving lowercase round-trip
are no longer emitted by any current writer.

### Reading embedded references

A `.tio` may embed one or more references at
`/study/references/<reference_uri>/`. Each child group is one
reference; under it sits a `chromosomes/` sub-group whose children
are per-chromosome sub-groups. Each per-chromosome sub-group
contains a `data` UINT8 dataset holding the raw sequence bytes
(uppercase ACGTN; current writers normalise case for codec
determinism — see "Historical note" above for pre-v1.1.1 layout).

URI-group attributes:

| Attribute        | Type   | Description                                                |
|------------------|--------|------------------------------------------------------------|
| `reference_uri`  | string | The reference URI (matches the group's leaf-name segment). |
| `md5`            | string | 32-character lowercase hex digest of the sequences.        |

Per-chromosome group attribute `length` carries the sequence length
in bases (informational; equals the `data` dataset's length and may
be ignored by readers).

**MD5 computation.** The `@md5` attribute on each
`/study/references/<uri>/` group is the MD5 of the **concatenated
sequence bytes** in **alphabetic order of chromosome name** —
i.e. `MD5(seq_for_chr_a || seq_for_chr_b || ...)` where `||`
denotes byte concatenation, with no inter-chromosome framing.
The same formula is used by all writers and by the
`compute_reference_md5` (Python) / `ReferenceImport.computeMd5`
(Java) / `+[TTIOReferenceImport computeMd5WithChromosomes:sequences:]`
(ObjC) helpers. Unified in v1.1.0 (previously: the REF_DIFF_V2
auto-embed path used this form while
`ReferenceImport.write_to_dataset` and the standalone helpers used
a name-framed `name + 0x0A + seq + 0x0A` form; see CHANGELOG for
details).

**External reference (`REF_PATH`).** When a run's reference is not
embedded, the readers resolve it from the FASTA named by the
`REF_PATH` environment variable (or an explicit resolver argument),
through its `.fai` index. The recorded 16-byte `reference_md5` is
checked against, in order: the md5 of the requested chromosome's
case-preserved bytes, of its upper-cased bytes (the pre-1.9 check,
which only matched a single-contig FASTA), and the reference-set md5
of the whole FASTA computed with the formula above (every
chromosome, alphabetic order, case preserved). A run written against
a multi-chromosome FASTA therefore validates against that same
FASTA. Since v1.9 all three implementations apply the same order and
return the chromosome upper-cased. The whole-FASTA digest is cached in
`<fasta>.ttio-md5` as `<hex> <size> <mtime_s>` by whichever
implementation computes it first, and reused by all of them while the
FASTA's size and modification time match; a location that cannot be
written is simply recomputed per process.

**Read-side accessor (v1.1.0+).** A freshly-opened dataset exposes
embedded references through the `references()` accessor
(`references` property in Python,
`-[TTIOSpectralDataset references]` in ObjC), keyed by reference
URI. Datasets written without embedded references (writer flag
`embedReference = false`) return an empty map regardless of whether
individual genomic runs carry a `referenceUri`.

**Write-side helper (v1.1.0+).** Outside the auto-embed path,
`ReferenceImport.write_to_dataset` (Python) /
`ReferenceImport.writeToDataset` (Java) /
`-[TTIOReferenceImport writeToDataset:overwrite:error:]` (ObjC)
embed a `ReferenceImport` directly at
`/study/references/<uri>/` on an open writable dataset. The
on-disk layout is byte-identical to the auto-embed path (same
`@md5` formula, same per-chromosome `chromosomes/<name>/data` +
`@length` shape). All three implementations honour an `overwrite`
parameter (default `False`/`NO`) — a write that collides with an
existing URI raises (`FileExistsError` /
`IllegalStateException` / `NSError` code 2201) unless overwrite
is explicitly requested.

The `@md5` attribute is preserved verbatim through `references()` —
i.e. the returned `ReferenceImport.md5` matches the on-disk
attribute byte-for-byte when present and well-formed. Missing or
malformed values fall back to recomputation via the canonical helper
(`compute_reference_md5` etc.), which means a read-back digest may
differ from the writer-stamped one in that specific edge case.

### Default codec selection

When a run provides `reference_chrom_seqs` (or otherwise resolves a
reference via the supplied resolver) AND `signal_codec_overrides` is
empty for `sequences`, the writer auto-applies REF_DIFF_V2 (§10.10b).
Without a reference, the channel falls through to BASE_PACK.

## 10.10b REF_DIFF_V2 — bit-packed sequence diff codec (codec id 14)

**Default in v1.8+.** Encodes the per-base diff stream as a 5-substream
blob (FLAG / BS / IN / SC / ESC) instead of v1's single rANS-encoded
bitstream. Saves ~4.5 MB on chr22 vs v1 by 2-bit-packing substitution
and IN/SC bases (was 8-bit literals in v1).

### 10.10b.1 On-disk schema

`signal_channels/sequences` is a GROUP (in v1 it's a flat dataset)
containing one child:

```
signal_channels/sequences/
└── refdiff_v2     uint8 1-D blob, @compression = 14 (REF_DIFF_V2)
```

The `refdiff_v2` blob carries an outer container header (38 + uri_len
bytes, magic "RDF2"), a slice index (32 bytes/entry, byte-compatible
with v1), and per-slice bodies each with a 24-byte sub-header + 5
rANS-O0-encoded substreams. Full wire format in the design spec §4.

Under `blocks_v1` (section 10.12) this dataset holds one such blob per block, back to back, addressed through `blocks/index`.

### 10.10b.2 Substream taxonomy

| Substream | Content | Encoding |
|-----------|---------|----------|
| FLAG | 1 byte per M/=/X cigar base; 0=match, 1=substitution | rANS-O0 |
| BS | 2-bit ACGT code per FLAG=1 base (substitution) | 2-bit packed → rANS-O0 |
| IN | 2-bit ACGT code per cigar I (insertion) base | 2-bit packed → rANS-O0 |
| SC | 2-bit ACGT code per cigar S (soft-clip) base | 2-bit packed → rANS-O0 |
| ESC | (stream_id, varint base_index, literal byte) per non-ACGT | rANS-O0 |

### 10.10b.3 Reader/writer dispatch

The reader probes the link type of `signal_channels/sequences`:
- GROUP → v2 path (`refdiff_v2` child + this codec id 14)
- DATASET → v1 path (codec id 9 / 6 / etc on the @compression attribute)

This means v1 → v2 transition produces a structural change at the
`sequences` link. v1.7 readers fail with "unknown link type" or
"sequences not a dataset" depending on the implementation.

### 10.10b.4 Backward compatibility

v1.7 reader on v1.8 file: encounters `sequences` as a GROUP, fails to
find a flat dataset there. Upgrade reader OR write the source file with
`opt_disable_ref_diff_v2 = True` to keep v1 layout.

v1.8 reader on v1.7 file: probes `open_group("sequences")` → KeyError →
falls through to v1 dispatch. Transparent.

### 10.10b.5 Cross-language byte-exact

Encode/decode primitives are a shared C kernel in libttio_rans
(`ttio_ref_diff_v2_encode/decode`), called by Python ctypes, Java JNI,
and ObjC direct linkage. The encoded byte stream is byte-exact identical
regardless of which language wrote the file. Verified by
`python/tests/integration/test_ref_diff_v2_cross_language.py` (3 PASS
+ 1 SKIP — hg002_pacbio's BAM has SEQ=`*`).

## 10.11 FQZCOMP_NX16_Z codec — CRAM-mimic quality codec

FQZCOMP_NX16_Z (codec id `12`) is the v1.0 default codec for the
`qualities` channel. Magic `M94Z`, V4 wire format. Clean-room
implementation of CRAM 3.1's `rANS-Nx16` discipline; see
`docs/codecs/fqzcomp_nx16_z.md` for the full algorithm and wire
format.

### On-disk schema

The `signal_channels/qualities` dataset under
`/study/genomic_runs/<run>/` carries `@compression == 12` and a flat
1-D `uint8` body (the M94.Z V4 codec stream — header + body as
specified in `docs/codecs/fqzcomp_nx16_z.md` §2). The dataset has
**no HDF5 filter applied** — codec output is high-entropy and
double-compression is a CPU loss.

Under `blocks_v1` (section 10.12) this dataset holds one such blob per block, back to back, addressed through `blocks/index`.

### Default codec selection

When a run has empty `signal_codec_overrides["qualities"]`, the
writer auto-applies FQZCOMP_NX16_Z. Override values
RANS_ORDER0 / RANS_ORDER1 / QUALITY_BINNED are accepted on the
`qualities` channel; BASE_PACK and the v2 codecs (13/14/15) are
rejected as wrong-content.

## 10.12 `blocks_v1` genomic block layout (v1.9)

From v1.9 every genomic writer emits a run as a sequence of
independently coded blocks. A block is a contiguous range of reads;
each blob channel of a block is coded with the same codec and wire
format as sections 10.6 to 10.11 define for a whole run, so a block's
blob is byte-identical to what a v1.8 writer would produce for a run
consisting of that block's reads alone. The only new structure is the
block index. `opt_legacy_whole_channel=True` on a run restores the
v1.8 whole-channel layout.

### 10.12.1 Run-level

| Attribute | Type | Meaning |
|---|---|---|
| `@layout` | fixed string `"blocks_v1"`, or `"blocks_v1_grouped"` (§10.12.7) | Selects this layout. Absent means the v1.8 whole-channel layout. A reader that does not know the value MUST fail with an unsupported-layout error. |
| `@block_policy` | fixed string | Writer policy, informative, e.g. `reads=1000000,bytes=67108864`. |
| `@read_count`, `@base_count` | int64 | Totals over the blocks written so far; updated at every flush and at close. |
| `@ref_diff_slice_bytes`, `@opt_disable_qualities_v5`, `@reference_md5s` | see §9.1.1 | Persisted writer policy and reference set for the per-AU restore re-encode (M99.1). Present only when non-default / applicable. |

### 10.12.2 `blocks/index`

Compound dataset, one row per block, extendable, no filter. Field
order is part of the contract:

| Field | Type | Meaning |
|---|---|---|
| `read_start` | uint64 | index of the block's first read |
| `n_reads` | uint32 | reads in the block |
| `base_start` | uint64 | offset of the block's first base in the concatenated base space |
| `n_bases` | uint64 | bases in the block |
| `sequences_off`, `sequences_len` | uint64 | byte range of the block's blob inside `sequences/data` |
| `qualities_off`, `qualities_len` | uint64 | same for `qualities` |
| `read_names_off`, `read_names_len` | uint64 | same for `read_names` |
| `cigars_off`, `cigars_len` | uint64 | same for `cigars` |
| `mate_info_off`, `mate_info_len` | uint64 | same for `mate_info/inline_v2` |
| `sequences_codec`, `qualities_codec`, `read_names_codec`, `cigars_codec`, `mate_info_codec` | uint32 | the codec id (section 10.4) of that block's blob, one column per channel in the same order |

A run that carries the `tags` channel (§10.13) has three more columns
after `mate_info_codec`: `tags_off`, `tags_len` (uint64) and
`tags_codec` (uint32). A run without tags has no such columns, so its
index is unchanged; readers treat absent tag columns as an empty tags
channel in every block. A writer whose first tagged block arrives after
untagged blocks rewrites the index with the three columns, the earlier
rows at `tags_len = 0`.

A channel a run does not carry has `_len = 0` in every row. The codec
columns exist because a block's codec can differ from its
neighbours': a mapped block codes sequences with REF_DIFF_V2 while
the unmapped block of the same run falls back to BASE_PACK. The
`@compression` attribute on a channel dataset carries the first
block's codec and is informative only.

### 10.12.3 Channel datasets

Each blob channel is one extendable 1-D `uint8` dataset holding the
blocks' blobs back to back: `sequences/data` (always a group under
this layout, whatever the codec), `qualities`, `read_names`,
`cigars`, `mate_info/inline_v2`, and `tags` when the run carries SAM
tags (§10.13). Codec output is unfiltered; a channel
whose codec is 0 keeps the zlib filter. Chunk size is 256 KiB.

A block never spans two chromosomes: the writer flushes the pending
block when the chromosome changes (REF_DIFF_V2 codes one chromosome
per blob; unmapped `*` reads form their own blocks). A
coordinate-sorted whole-genome BAM therefore streams through
REF_DIFF_V2 as long same-chromosome stretches.
Blocks are independent, so writers may encode
several at once; the file does not record thread counts and is
identical whatever the count.

Every blob channel is codec-coded under this layout. Defaults when the
caller sets no override: `cigars` RANS_ORDER0 (id 4, section 10.8;
the v1.8 compound VL-string default has no blob form), `qualities`
FQZCOMP_NX16_Z (id 12; v1.8 only used it when another v1.5 codec was
active on the run and otherwise left zlib-filtered raw bytes), and
`sequences` REF_DIFF_V2 (id 14) with a reference, or without one
SEQ_CM (id 19), falling back to RANS_ORDER1 (id 5) when the writer has
no native library; the block's `sequences_codec` column records which. Unmapped reads inside a mapped block (FLAG 0x4,
CIGAR `*`, placed on the mate's contig) stay in the REF_DIFF_V2 blob:
the codec carries their bases as soft clip and their lengths in the
slice's UL substream (docs/codecs/ref_diff_v2.md section 4.4). A block
that contains a zero-length read (SEQ `*`, e.g. a secondary
alignment) codes qualities with RANS_ORDER0 because the FQZCOMP_NX16_Z
kernel does not decode a stream in which a zero-length read precedes
another read; the codec column records it.

`mate_info/chrom_names`, `genomic_index/chromosome_names` and the
reference tables stay run-level; chromosome ids are assigned once per
run in first-seen order across blocks and both name tables are
written at close from that map. Readers derive a read's own
chromosome id for the mate decoder from the `mate_info/chrom_names`
row index (encounter order restarts per block, so it must not be
rebuilt from the block's own names).

Codec consequences the codecs already tolerate: FQZCOMP_NX16_Z
auto-tunes per block, REF_DIFF_V2 carries its slice index per block,
NAME_TOKENIZED_V2 restarts its tokenizer per block. Measured on the
NA12878 chr22 low-coverage BAM (151 MB) with the hs37 chr22 reference:
blocks_v1 73.5 MB, whole-channel layout 113.4 MB (the whole-channel
writer falls back to BASE_PACK on the whole channel when any read is
unmapped); NA12878 WES chr22 (72.8 MB): 64.1 MB vs 87.5 MB. Peak RSS
of the streaming import at the default 1 M-read blocks: 1.8 GB and
1.6 GB.

### 10.12.4 `genomic_index/`

`lengths`, `positions`, `mapping_qualities`, `flags`,
`chromosome_ids` keep their v1.8 element types and are extendable
chunked datasets. `read_start`/`base_start` in the block index give
`run[i]` its block with one binary search. A `blocks_v1_grouped` run
adds `input_index` (§10.12.7).

### 10.12.5 Close and partial files

A writer that stops before close leaves a file whose block index and
datasets agree up to the last flushed block; readers ignore bytes past
the last indexed block and use the index row count, not
`@read_count`, when the two disagree.

### 10.12.6 Signatures

`sign_genomic_run` / `verify_genomic_run` cover the same datasets as
for the whole-channel layout (the datasets inside a channel group,
`sequences/data`, included), the `tags` channel when present (M101),
plus `blocks/index` (canonical compound bytes) and, for a
`blocks_v1_grouped` run, `genomic_index/input_index`. Per-AU encryption
walks this layout block by block (§9.1.1, M99) and encrypts
`sequences`, `qualities` and, when present, `tags`. Region encryption
operates on the whole-channel layout only and refuses a run that
carries tags.

Cross-language: Java `GenomicStreamWriter` / ObjC
`TTIOGenomicStreamWriter` and their readers implement this section;
the golden fixture is `python/tests/fixtures/genomic/blocks_v1_golden.tio`.

### 10.12.7 `blocks_v1_grouped` — reads grouped by sequence (M103)

An unaligned run (no reference) can be written with its reads
reordered so that reads from the same place in the genome share a
block, which is what SEQ_CM (codec id 19) needs to see each read's
overlap partners (`docs/codecs/seq_cm.md` §6). Writers do this only on
request (`group_reads` / `opt_group_reads`); a run with a reference
keeps its order for REF_DIFF_V2 and is refused.

The run is a `blocks_v1` run in every respect except:

- `@layout` is `"blocks_v1_grouped"`. Readers that predate it reject
  the run as an unknown layout instead of reading it in stored order.
- Every channel and every `genomic_index` column holds the reads in
  stored (grouped) order, and `genomic_index/input_index` (uint32,
  extendable, chunked and compressed like the other index columns)
  holds, at row `j`, the input index of the read stored at row `j`.
  It is a permutation of `0 .. read_count - 1`; a reader MUST refuse a
  run whose `input_index` has the wrong length or is not a permutation.
- The order is `ttio_seq_group` (`native/include/ttio_rans.h`) at its
  default parameters, applied within each chromosome label (labels in
  first-seen order), so a block still never spans two labels. The
  kernel is shared by the three SDKs and breaks every tie on a total
  order, so they write identical files.

Readers present the run in input order: `run[i]` is the `i`-th read as
imported (one lookup through the inverse permutation, one block
decode), the index arrays are permuted into input order (`offsets`
recomputed from the input-order lengths), sequential iteration gathers
input-order chunks block by block, and exporters therefore restore the
original order. Block-level iteration stays in stored order and maps
its rows through `input_index`. Per-AU encryption and restore walk the
stored rows block by block and leave `input_index` as it is; plaintext
transport sends reads in input order (never in bulk mode) and the
receiver writes an ordinary run; the encrypted transport carries
`input_index` in each BlockSidecar (transport-spec v0.13, §4.24).

## 10.13 `tags` channel — SAM optional fields (M101)

A genomic run may carry the SAM optional fields of its reads: for each
read, SAM columns 12 and up joined by TAB, exactly as `samtools view`
prints them (`""` for a read without tags). They are stored as
`signal_channels/tags`, a flat `uint8` dataset with
`@compression = 18` (SAM_TAGS, `docs/codecs/sam_tags.md`), one blob
per block under `blocks_v1` and one blob for the run under the
whole-channel layout. A run none of whose reads carries a tag has no
`tags` dataset, and a file with at least one such dataset lists
`opt_sam_tags` in `@ttio_features`. Readers without M101 ignore the
dataset and lose the tags; they read everything else unchanged.

Fidelity is the SAM text. BAM stores integer tags with a width (`c`,
`C`, `s`, `S`, `i`, `I`) that `samtools view` prints as `i`; the width
is not kept, and an exporter that writes BAM lets samtools choose it
again. Float (`f`) and array (`B`) values keep the text samtools
printed.

MD:Z and NM:i are recomputed from the reference (SAM_TAGS `DERIVED`
entries) only in a blob whose block codes `sequences` with REF_DIFF_V2,
so a reader that can decode the sequences can decode the tags. The
reference is the one the REF_DIFF_V2 blob header names, resolved the
same way (embedded under `/study/references/`, or `REF_PATH`). Any
other blob stores MD and NM like other tags.

Per-AU encryption (§9.1.1) encrypts the tags with the sequences, one
AU per read whose plaintext is the read's tag text; the encrypted
container holds `tags_segments` in place of `tags`. A whole-channel
run with tags refuses per-AU and region encryption rather than leave
MD strings in plaintext.

## 11. Subjects + Samples (v0.11)

`/study/subjects/` and `/study/samples/` carry per-dataset cohort
metadata: who the data came from (`Subject`) and what was collected
from them (`Sample`). Each is a parent group containing one HDF5
**sub-group per row**, keyed by the row's primary key. The
attribute schema mirrors the data model in the v1.4 subjects/samples
design.

The transport-stream encoding of these two tables is documented in
[`transport-spec.md`](transport-spec.md) §4.22 (`SUBJECT_METADATA`
0x19 + `SAMPLE_METADATA` 0x1A, length-prefixed Arrow IPC).

### 11.1 Layout

```
/study/subjects/
└── <external_id>/                       (HDF5 group, one per Subject)
    # Attributes (all optional except external_id):
    @external_id        : VL string  — primary key; matches the group leaf name
    @project            : VL string  — empty string when absent in source
    @sex                : VL string  — empty string when absent in source
    @birth_year         : int64      — sentinel 0 = unknown
    @attributes_json    : VL string  — sort_keys JSON object; "{}" when empty

/study/samples/
└── <sample_id>/                         (HDF5 group, one per Sample)
    # Attributes (all optional except sample_id):
    @sample_id            : VL string  — primary key; matches the group leaf name
    @subject_external_id  : VL string  — soft FK; empty string when absent in source
    @sample_kind          : VL string  — empty string when absent in source
    @collected_at         : int64      — unix seconds since epoch; sentinel 0 = unknown
    @attributes_json      : VL string  — sort_keys JSON object; "{}" when empty
```

Per-row groups (rather than a single compound dataset) are used
because HDF5 compound datasets are awkward to extend (adding a
column rewrites the table) and do not accept variable-length string
fields cleanly. Per-row groups inspect cleanly with `h5dump` and let
`attributes_json` carry the open-extension slot without
table-schema gymnastics. The cost is metadata overhead for large
cohorts (≈ hundreds of bytes per row), acceptable for v1.4.

### 11.2 Attribute schemas

**Subject attributes** (on `/study/subjects/<external_id>/`):

| Attribute          | Type      | Required | Notes                                                                                |
|--------------------|-----------|----------|--------------------------------------------------------------------------------------|
| `external_id`      | VL string | **yes**  | Unique within the dataset. Matches the group leaf name verbatim. Primary key.        |
| `project`          | VL string | no       | Free string, often a study acronym. Empty string `""` when absent in source.         |
| `sex`              | VL string | no       | Free string (e.g. `"M"`, `"F"`, `"NA"`). Empty string `""` when absent in source.    |
| `birth_year`       | int64     | no       | 4-digit YYYY. Sentinel `0` denotes unknown.                                          |
| `attributes_json`  | VL string | no       | Open-extension slot: JSON object, `sort_keys=true`, separators `","`/`":"`. `"{}"` when empty (never absent). See §11.4. |

**Sample attributes** (on `/study/samples/<sample_id>/`):

| Attribute              | Type      | Required | Notes                                                                                |
|------------------------|-----------|----------|--------------------------------------------------------------------------------------|
| `sample_id`            | VL string | **yes**  | Unique within the dataset. Matches the group leaf name verbatim. Primary key.        |
| `subject_external_id`  | VL string | no       | Soft foreign key to `Subject.external_id` in this dataset (see §11.3). Empty string `""` when absent in source. |
| `sample_kind`          | VL string | no       | Free string. Empty string `""` when absent in source.                                |
| `collected_at`         | int64     | no       | Unix seconds since epoch. Sentinel `0` denotes unknown.                              |
| `attributes_json`      | VL string | no       | Open-extension slot — same convention as on Subject. `"{}"` when empty.              |

### 11.3 ID safety, validation, and soft-FK

**ID safety.** Because `external_id` and `sample_id` are used
verbatim as HDF5 group names, they MUST be valid HDF5 group names:
non-empty and containing no `/` characters (the HDF5 path separator).
Writers SHALL validate these constraints and raise on violation.
Readers SHOULD tolerate corruption — an unreadable row name surfaces
as a warning and the row is skipped, mirroring the §11.4 forward
compatibility posture.

**Soft-FK.** `Sample.subject_external_id` is a *soft* reference:
- Present and matching a Subject in this dataset → consistent.
- Present but no matching Subject in this dataset → log a WARNING
  (NOT an error). This allows datasets to ship Samples without
  full Subject metadata.
- Absent (empty string) → fine. Anonymous samples are valid.

**Duplicates.** Duplicate `external_id` or `sample_id` on write
SHALL raise. The HDF5 group-name uniqueness check enforces this at
the on-disk layer.

### 11.4 Cardinality and the empty case

- Zero or more Subjects per dataset.
- Zero or more Samples per dataset.

If a dataset has **no Subjects**, `/study/subjects/` is **absent
entirely** (NOT an empty group). Likewise for `/study/samples/`.
Readers MUST treat an absent parent group as zero rows and surface
the corresponding accessor as an empty list — never an error.

### 11.5 Relationship to `AcquisitionRun.sampleName`

`AcquisitionRun.sampleName` (a free string on each acquisition run;
see §3) remains the canonical run→sample link. When both Sample
rows in `/study/samples/` and `AcquisitionRun.sampleName` are
present, applications SHOULD treat `sampleName` as a foreign key
into `Sample.sample_id`: a run's `sampleName` matches one Sample's
`sample_id`. No automatic enrichment is performed — Sample lookup
is the caller's job.

Anonymous runs (no matching Sample row) remain valid. Mismatched
`sampleName` values surface as application-level data quality
warnings, NOT reader errors. Datasets without explicit Subject
or Sample groups are unchanged in interpretation: the run→sample
link stays a free string.

### 11.6 `attributes_json` byte form

`attributes_json` MUST be encoded as a JSON object with:

- `sort_keys=true` — keys appear in code-point sorted order.
- No whitespace between tokens; separators `","` (item) and `":"`
  (key-value).
- An empty map encodes to the two-byte string `"{}"` (NEVER absent,
  NEVER literal HDF5 null).

This matches Python `json.dumps(d, sort_keys=True,
separators=(",", ":"))`, the Java reference writer's
`TreeMap`-walk encoder, and the Objective-C reference writer's
`NSJSONWritingSortedKeys` option. Byte equality across the three
SDKs is required for cross-language conformance and stable
content hashes.

### 11.7 Backward compatibility

`/study/subjects/` and `/study/samples/` were introduced in v0.11.
Pre-v0.11 readers naturally see neither group (the writer simply
does not emit them when the caller passes no subjects/samples) and
materialise empty cohort accessors. v0.11+ readers opening a
pre-v0.11 file MUST surface empty `subjects()` / `samples()`
lists and MUST NOT raise. No feature flag gates the layout — the
group's presence is itself the signal — though the transport
layer's `transport_v0_11` flag covers the wire-format counterpart
(`SUBJECT_METADATA` 0x19 + `SAMPLE_METADATA` 0x1A; see
transport-spec §4.22).

---

## 11a. Assembly graphs (M98)

One assembly graph stores a GFA 1.x file at
`/study/assembly_graphs/<name>/` so that re-emission is byte-exact
against the parsed input. The `assembly_graphs/` group carries a
`@_graph_names` comma-list attribute (the `_run_names` convention);
files with at least one graph set the `opt_assembly_graph` feature
flag, and graph-less files carry neither the flag nor the subtree.

Per-graph attributes:

| Attribute | Type | Semantics |
|---|---|---|
| `@gfa_version` | UTF-8 string | The `VN:Z:` value of the first header line routed to extras, else `"1.0"`. |
| `@producer` | UTF-8 string | Reserved; `""` from the GFA readers. |
| `@final_newline` | integer | 1 when the source file ended in LF, else 0. Structural marker: a group without this attribute is not an M98 graph. |

Per-graph children, all 1-D:

| Child | Kind | Fields / element |
|---|---|---|
| `segments/records` | compound | `name` VL_STRING, `length` UINT64, `seq_offset` UINT64, `seq_missing` UINT32, `tags` VL_STRING |
| `segments/sequences` | uint8 | Concatenation of every present segment sequence; `seq_offset`/`length` index into the decoded bytes. Absent when every segment carries `*`. |
| `links` | compound | `from`, `from_orient`, `to`, `to_orient`, `overlap`, `tags` — all VL_STRING |
| `paths` | compound | `name`, `segment_list`, `overlaps`, `tags` — all VL_STRING |
| `extras` | compound | `line` VL_STRING — the verbatim source line |
| `line_index` | compound | `line_type` UINT32 (0 = S, 1 = L, 2 = P, 3 = extra), `row` UINT64 |

Empty tables are ABSENT: a 0-row non-extendable compound does not
round-trip on every provider, and readers in all 3 SDKs treat a
missing table as empty.

**Parser structural rules** (identical in the 3 SDKs): split on LF
after final-newline detection; fields split on TAB with trailing
empty fields preserved. `S` needs at least 3 fields, `L` 6, `P` 4;
every other line (`H`, `C`, comments, hifiasm `A` lines, short
S/L/P) goes verbatim into `extras`. `tags` is the tab-joined
verbatim remainder, `""` when none. A `*` sequence stores no bytes
and sets `seq_missing`. An empty file has 0 lines; a file holding
one LF is one empty extras line. Emission replays `line_index`,
writes `*` for a missing sequence, appends `tags` with a TAB only
when non-empty, and restores the final newline from the attribute.

**Sequences codec**: `segments/sequences` is stored through
BASE_PACK (id 6) when every byte is in `ACGTNacgtn`, RANS_ORDER1
(id 5) otherwise, uncompressed when empty; the codec id sits in the
dataset's `@compression` attribute per §10.5 (0 or absent = raw).

**Protection**: per-AU encryption (§9.1) covers
`segments/sequences` with one AES-256-GCM operation per segment
record — offsets and lengths from `segments/records`, AAD channel
name `sequences`, dataset_id continuing after the genomic runs. The
channel is decoded to raw bytes before slicing; decrypt-in-place
restores the raw uint8 dataset with no `@compression`. `v2:`/`v3:`
signatures (§10) cover `segments/sequences` the way they cover
genomic-run channels; the compound tables are outside the signed
set.

---

## 12. Backward compatibility

A v0.2 reader recognizes a v0.1 file by the absence of
`@ttio_features`. The fallback paths:

| v0.2 location                              | v0.1 fallback                              |
|--------------------------------------------|--------------------------------------------|
| `/study/identifications` compound dataset  | `/study/@identifications_json` string attr |
| `/study/quantifications` compound dataset  | `/study/@quantifications_json` string attr |
| `/study/provenance` compound dataset       | `/study/@provenance_json` string attr      |
| `/study/image_cube/` group                 | `/image_cube/` at root                     |
| `spectrum_index/headers` compound dataset  | parallel 1-D datasets only (still present) |
| `intensity_matrix_2d` rank-2 dataset       | `arrays/intensity_matrix` flattened 1-D    |

All v0.1 files written by libTTIO v0.1.0-alpha are readable by v0.2.0
code without modification.

### 12.1 Compound JSON mirror (v0.6 / M37)

All three writers (ObjC, Python, Java) emit the `*_json` string
attribute <em>alongside</em> the `/study/{identifications,quantifications,provenance}`
compound dataset. The attribute carries the same records as the
compound, encoded as a JSON array of objects:

```
[{"run_name": "...", "spectrum_index": 0, "chemical_entity": "...",
  "confidence_score": 0.95, "evidence_chain": ["..."]}, ...]
```

The mirror exists because the HDF5 Java binding JHI5 1.10.x (the
version on current apt/homebrew HDF5 packages) cannot marshal
variable-length-string fields out of a compound dataset: the JNI
rejects any H5Dread whose mem-type contains `H5T_STRING` or
`H5T_VLEN`. Java therefore prefers the JSON attribute on read and
only falls back to primitive-field projection of the compound when
the mirror is absent, which yields empty strings for every VL field.

The mirror is <em>not</em> emitted for the per-run `<run>/provenance/steps`
compound dataset (§6.4) — per-run provenance was added after the
v0.2 attribute fallback and Java's reader does not descend into
run-level compound metadata.

Sealing (§10): encryption moves the compound dataset into
`*_sealed` blobs; the sealing code deletes the matching `*_json`
attribute so sealed files stay opaque without the key.

A future release will remove the mirror once Java's HDF5 binding
gains compound-with-VL read support, at which point every reader
will go through the compound dataset directly.

### 12.2 Per-AU encrypted layout (v0.10+)

A v0.10 reader recognises a per-AU-encrypted file by
`opt_per_au_encryption` in `@ttio_features`. The channel layout
under `signal_channels/` flips from plaintext `<channel>_values`
to the `<channel>_segments` compound (§9.1) with VL_BYTES members
for `iv` / `tag` / `ciphertext`. Pre-v0.10 readers without
VL_BYTES compound support will refuse to open the file (the flag
is non-optional when set).

When `opt_encrypted_au_headers` is also present, the six plaintext
index arrays (retention_times, ms_levels, polarities,
precursor_mzs, precursor_charges, base_peak_intensities) are
absent; the semantic header travels in the
`spectrum_index/au_header_segments` compound instead. `offsets`
and `lengths` remain plaintext because they frame the compound
rows.

Migration: `python -m ttio.tools.per_au_cli transcode` rewrites
a plaintext or previously-encrypted file through the v0.10 path,
with optional `--rekey`. v0.x `opt_dataset_encryption` files must
be decrypted via the v0.x `SpectralDataset.decrypt()` API first —
channel-level AES-GCM cannot be converted in place (the plaintext
must be materialised to re-slice per-spectrum).

---

## 13. Example `h5dump` output (minimal MS file)

```
HDF5 "minimal_ms.tio" {
GROUP "/" {
   ATTRIBUTE "ttio_format_version" { DATATYPE H5T_STRING ... DATA { "1.0" } }
   ATTRIBUTE "ttio_features" { DATA { "[\"base_v1\",\"compound_identifications\",...]" } }
   GROUP "study" {
      ATTRIBUTE "title" { ... }
      GROUP "ms_runs" {
         ATTRIBUTE "_run_names" { DATA { "run_0001" } }
         GROUP "run_0001" {
            ATTRIBUTE "acquisition_mode"  { DATA { 0 } }
            ATTRIBUTE "spectrum_count"    { DATA { 10 } }
            ATTRIBUTE "spectrum_class"    { DATA { "TTIOMassSpectrum" } }
            GROUP "instrument_config" { ... }
            GROUP "spectrum_index" {
               DATASET "offsets"               { DATATYPE H5T_STD_U64LE ... }
               DATASET "lengths"               { DATATYPE H5T_STD_U32LE ... }
               DATASET "retention_times"       { DATATYPE H5T_IEEE_F64LE ... }
               DATASET "ms_levels"             { DATATYPE H5T_STD_I32LE ... }
               DATASET "polarities"            { DATATYPE H5T_STD_I32LE ... }
               DATASET "precursor_mzs"         { DATATYPE H5T_IEEE_F64LE ... }
               DATASET "precursor_charges"     { DATATYPE H5T_STD_I32LE ... }
               DATASET "base_peak_intensities" { DATATYPE H5T_IEEE_F64LE ... }
               DATASET "headers"               { DATATYPE H5T_COMPOUND { ... } }
            }
            GROUP "signal_channels" {
               ATTRIBUTE "channel_names" { DATA { "mz,intensity" } }
               DATASET "mz_values"        { DATATYPE H5T_IEEE_F64LE ... }
               DATASET "intensity_values" { DATATYPE H5T_IEEE_F64LE ... }
            }
         }
      }
      GROUP "nmr_runs" { ATTRIBUTE "_run_names" { DATA { "" } } }
   }
}
}
```

---

## 14. Conformance checklist for a new reader

1. Open the file with any HDF5 library.
2. Read `@ttio_format_version`. If absent, treat as v0.1.
3. Read `@ttio_features` as a JSON array. Refuse the file if any
   non-`opt_`-prefixed feature is unknown.
4. Open `/study/`. Read `@title`, `@isa_investigation_id`.
5. Enumerate runs under `/study/ms_runs/` via the `_run_names`
   attribute.
6. For each run, read `@spectrum_class` (default
   `TTIOMassSpectrum`). Read `spectrum_index/` parallel datasets.
7. Read `signal_channels/@channel_names` (default `"mz,intensity"`).
8. For random access to spectrum `i`: slice each
   `<channel>_values[offsets[i] : offsets[i] + lengths[i]]` and
   reconstruct the spectrum from the resulting buffers plus the
   index metadata.
9. For identifications/quantifications/provenance, prefer the
   compound datasets if the feature flags are set; fall back to the
   JSON attributes otherwise.
10. If `@encrypted` is set and no key is available, report the file
    as read-only-with-restrictions and skip the encrypted channels.

A conforming reader need not implement writing. The spec is
deliberately write-once: the Objective-C reference implementation
handles the complex sealing and compound marshalling paths, and
third-party readers are expected to consume files written by it.
