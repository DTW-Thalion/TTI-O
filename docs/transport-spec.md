# TTI-O Transport Format Specification — v0.11.0

This document specifies the TTI-O streaming transport format: a
self-describing binary framing protocol that carries the same logical
data as the `.tio` file format over sequential byte streams and
network connections. It is the MPEG-G ISO/IEC 23092-1 §7 ("Transport
Format") equivalent for TTI-O.

The file format (see [`format-spec.md`](format-spec.md)) is the
authoritative at-rest representation. The transport format is a
streaming wire protocol. The two are bidirectionally convertible with
no loss of information on signal data, annotations, provenance, or
protection metadata.

Transport streams use the `.tis` file extension (TTI-O Transport
Stream) when serialized to a file. This parallels MPEG-G's `.mgg`
(file) vs `.mggts` (transport) convention.

The codec, the WebSocket client / server, the acquisition simulator,
bidirectional conformance, selective access, ProtectionMetadata, and
the v1.0 per-Access-Unit encryption paths referenced throughout this
document all shipped in the v0.10.0 release across the Python, ObjC,
and Java reference implementations. The inline "v1.0" labels on
encrypted variants denote the encryption feature level, not a future
version; these variants are live and exercised by the cross-language
conformance harness in `tests/integration/test_per_au_cross_language.py`.
See `docs/transport-encryption-design.md` for the normative design
notes and `docs/format-spec.md` §9.1 for the on-disk layout.

## 1. Design Goals

1. **Bit-identical signal round-trip.** File → transport → file
   produces the same m/z, intensity, chromatogram, and spectral
   values to within float64 epsilon; encoded/compressed channel
   bytes are preserved verbatim when the receiver's container can
   accept them.
2. **Streamable.** Packets are self-delimiting and parseable in a
   single forward pass with fixed-size header reads. No backward
   seeks are required to materialize a dataset.
3. **Selective access.** Access Unit headers expose the filter keys
   (RT, MS level, polarity, precursor m/z) needed for server-side
   filtering without decoding the signal payload.
4. **Multiplexing.** AUs from different runs within one dataset
   group may be interleaved on the wire; the receiver demultiplexes
   by `dataset_id`.
5. **Provider-agnostic.** Transport writers iterate the source
   provider's `SpectrumIndex`; transport readers write through the
   target provider's `StorageGroup` API. The transport layer never
   touches HDF5 / SQLite / Zarr directly.
6. **Encryption-preserving.** Datasets encrypted at rest stream as
   ciphertext AUs with a `ProtectionMetadata` prologue. The receiver
   reconstructs the file-format protection model without decrypting
   in transit.
7. **Little-endian, explicit.** All multi-byte integers and IEEE-754
   floats are little-endian on the wire, independent of host byte
   order.

## 2. Versioning

Every transport stream carries `version = 0x01` in every packet
header. Version 1 is defined by this document.

- **Major** (`1 → 2`): incompatible framing changes. Version-1 readers
  MUST refuse to parse version-2 streams.
- **Stream-level feature flags** travel in the `StreamHeader`
  `features` list and inherit semantics from
  [`feature-flags.md`](feature-flags.md): flags without an `opt_`
  prefix are required; `opt_`-prefixed flags are informational.

### 2.1 Spec version history

| Version | Date       | Changes |
|---------|------------|---------|
| v0.12   | 2026-08-25 | blocks_v1 per-AU carriage: `GENOMIC_RUN_SIDECAR` (0x1C) + `BLOCK_SIDECAR` (0x1D) packets and the required `transport_blocks_v1` feature token (§4.24). Senders no longer refuse blocks_v1 containers. |
| v0.11   | 2026-05-25 | Complete `.tio` coverage: references, image cubes, identifications, quantifications, dataset-level provenance, subjects/samples metadata, encryption-algorithm name. New packet types 0x10–0x1B. `SUBJECT_METADATA` (0x19) + `SAMPLE_METADATA` (0x1A) wire format finalised 2026-05-26 (§4.22). Backward-compatible: v0.10 readers skip unknown packet types via length-prefixed wire frames; readers that need v0.11 semantics check the `transport_v0_11` feature flag in StreamHeader. |

## 3. Wire Format

### 3.1 Packet Header

Every packet begins with a 24-byte fixed header. All fields are
little-endian.

| Offset | Field            | Type       | Value / Meaning                                    |
|-------:|------------------|------------|----------------------------------------------------|
|      0 | `magic`          | `bytes[2]` | ASCII `"TI"` (0x54 0x49)                           |
|      2 | `version`        | `uint8`    | `0x01` in this spec                                |
|      3 | `packet_type`    | `uint8`    | See §3.2                                           |
|      4 | `flags`          | `uint16`   | Bit flags — see §3.1.1                             |
|      6 | `dataset_id`     | `uint16`   | Identifies the AcquisitionRun (0 for stream-scope packets) |
|      8 | `au_sequence`    | `uint32`   | Monotonic per `dataset_id`; 0 for non-AU packets   |
|     12 | `payload_length` | `uint32`   | Bytes of payload following the header              |
|     16 | `timestamp_ns`   | `uint64`   | Nanosecond Unix timestamp at packet emission       |

`HEADER_SIZE = 24`. Readers MUST validate `magic == "TI"` and reject
the stream otherwise.

### 3.1.1 Flag bits

| Bit | Name                | Meaning                                                     |
|----:|---------------------|-------------------------------------------------------------|
|   0 | `ENCRYPTED`         | Channel data in this AU is AES-GCM encrypted (v0.10+)        |
|   1 | `COMPRESSED`        | Reserved (packet-level compression; unused today)           |
|   2 | `HAS_CHECKSUM`      | Payload is followed by a 4-byte CRC-32C (§3.3)              |
|   3 | `ENCRYPTED_HEADER`  | AU semantic header is also AES-GCM encrypted (v0.10+)        |

Readers MUST reject `ENCRYPTED_HEADER` without `ENCRYPTED` —
encrypting the filter header while leaving channel data
plaintext is not a meaningful mode.

### 3.2 Packet Types

| Code | Name                  | Purpose                                                  |
|-----:|-----------------------|----------------------------------------------------------|
| 0x01 | `StreamHeader`        | Format version, dataset group metadata, feature flags    |
| 0x02 | `DatasetHeader`       | One AcquisitionRun's identity and channel layout         |
| 0x03 | `AccessUnit`          | One spectrum's signal + filter metadata                  |
| 0x04 | `ProtectionMetadata`  | Cipher suite, wrapped DEK, signature public key          |
| 0x05 | `Annotation`          | One identification or quantification record              |
| 0x06 | `Provenance`          | One processing-step record                               |
| 0x07 | `Chromatogram`        | A batch of chromatogram data points                      |
| 0x08 | `EndOfDataset`        | Terminates a specific `dataset_id`                       |
| 0x09 | `BlobV2MateInfo`      | Verbatim `mate_info/inline_v2` blob (bulk mode, §4.10)   |
| 0x0A | `BlobV2RefDiff`       | Verbatim `sequences/refdiff_v2` blob (bulk mode, §4.11)  |
| 0x0B | `BlobV2NameTok`       | Verbatim `read_names/name_tok_v2` blob (bulk mode, §4.12)|
| 0x10 | `ReferenceGroupHeader`   | Reference-import header: URI + chromosome count + md5 (v0.11, §4.13) |
| 0x11 | `ReferenceChromosome`    | One contig within a reference group (v0.11, §4.14)       |
| 0x12 | `EndOfReferenceGroup`    | Terminates a reference group (v0.11, §4.15)              |
| 0x13 | `ImageHeader`            | Imaging cube grid + axis metadata (v0.11, §4.16)         |
| 0x14 | `ImagePixel`             | One pixel of an imaging cube (v0.11, §4.17)              |
| 0x15 | `EndOfImage`             | Terminates an imaging cube (v0.11, §4.18)                |
| 0x16 | `IdentificationsTable`   | Full identifications table as an Arrow IPC stream (v0.11, §4.19) |
| 0x17 | `QuantificationsTable`   | Full quantifications table as an Arrow IPC stream (v0.11, §4.20) |
| 0x18 | `DatasetProvenance`      | Dataset-level provenance chain (v0.11, §4.21)            |
| 0x19 | `SubjectMetadata`        | Subjects table as an Arrow IPC stream (v0.11, §4.22)     |
| 0x1A | `SampleMetadata`         | Samples table as an Arrow IPC stream (v0.11, §4.22)      |
| 0x1B | `EncryptionAlgorithm`    | Dataset-level `@encrypted` algorithm name (v0.11, §4.23) |
| 0xFF | `EndOfStream`         | Terminates the entire transport stream                   |

### 3.3 Checksum

When `flags & 0x04` (`has_checksum`), the payload bytes are followed
by a 4-byte CRC-32C (Castagnoli polynomial 0x1EDC6F41, initial
0xFFFFFFFF, final XOR 0xFFFFFFFF) of the payload. The checksum is NOT
counted in `payload_length`. Receivers MUST verify the checksum when
the flag is set and MAY reject the packet on mismatch.

CRC-32C was chosen for hardware acceleration availability
(`_mm_crc32_u64` on x86-64, PMULL on ARMv8, `java.util.zip.CRC32C`
on JDK 9+, `google-crc32c` on Python).

## 4. Payload Formats

Integer sizes match the field type; strings are UTF-8 and NOT
NUL-terminated.

### 4.1 StreamHeader (`0x01`)

Exactly one StreamHeader MUST appear as the first packet of every
transport stream.

```
format_version_len:  uint16
format_version:      bytes[format_version_len]   # e.g. "1.2"
title_len:           uint16
title:               bytes[title_len]
isa_id_len:          uint16
isa_investigation:   bytes[isa_id_len]
n_features:          uint16
features:            repeated { uint16 len, bytes[len] }
n_datasets:          uint16                       # number of DatasetHeaders to follow
```

`format_version` SHOULD match the `ttio_format_version` root
attribute of the source file (e.g. `"1.2"`) so receivers can enable
version-appropriate parsing of nested JSON payloads.

### 4.2 DatasetHeader (`0x02`)

One DatasetHeader per AcquisitionRun *or* GenomicRun. MUST precede
any AccessUnit carrying the matching `dataset_id`.

```
dataset_id:          uint16                       # matches header.dataset_id
name_len:            uint16
name:                bytes[name_len]              # run name, e.g. "run_0001"
acquisition_mode:    uint8                        # TTIOAcquisitionMode enum
spectrum_class_len:  uint16
spectrum_class:      bytes[spectrum_class_len]    # e.g. "TTIOMassSpectrum",
                                                  # "TTIOGenomicRead" (M89.2)
n_channels:          uint8
channel_names:       repeated { uint16 len, bytes[len] }
instrument_json_len: uint32
instrument_json:     bytes[instrument_json_len]   # InstrumentConfig OR
                                                  # genomic-run metadata JSON
expected_au_count:   uint32                       # 0 if unknown / real-time
```

`spectrum_class` uses the ObjC class name as the canonical wire
token; Python and Java map it to their native class via
`class_hierarchy.md`.

For genomic runs (`spectrum_class == "TTIOGenomicRead"`), the
`instrument_json` slot carries
``{"reference_uri", "platform", "sample_name", "modality",
"read_role"}`` instead of an `InstrumentConfig` (M89.2; `read_role`
added in M97 — always present, `""` means the run has no
`@read_role` attribute, receivers leave the container attribute
absent on `""`, and pre-M97 streams simply lack the key). Channel
names for genomic runs are `["sequences", "qualities"]` in M89.2;
compound channels (cigars, read_names, mate_*) are not yet carried
on the wire.

`expected_au_count = 0` is the real-time acquisition signal — the
receiver MUST allocate growable structures and MUST NOT pre-size
indices.

### 4.3 AccessUnit (`0x03`)

One AU carries exactly one spectrum (or one pixel of an MS image).
The payload takes one of three forms selected by the packet-header
flags:

- **Plaintext** (neither `ENCRYPTED` nor `ENCRYPTED_HEADER`) —
  §4.3.1.
- **Encrypted channel data** (`ENCRYPTED` set) — §4.3.2.
- **Fully encrypted** (`ENCRYPTED | ENCRYPTED_HEADER` set) —
  §4.3.3.

#### 4.3.1 Plaintext AU

```
spectrum_class:      uint8        # 0=MassSpectrum, 1=NMRSpectrum,
                                  # 2=NMR2D, 3=FID, 4=MSImagePixel,
                                  # 5=GenomicRead (M89.1)
acquisition_mode:    uint8        # TTIOAcquisitionMode
ms_level:            uint8        # 1, 2, ... (0 for NMR)
polarity:            uint8        # 0=positive, 1=negative, 2=unknown
retention_time:      float64      # seconds
precursor_mz:        float64      # 0.0 if MS1 or NMR
precursor_charge:    uint8        # 0 if MS1 or NMR
ion_mobility:        float64      # 0.0 if not applicable
base_peak_intensity: float64
n_channels:          uint8

# Repeated n_channels times:
channel_name_len:    uint16
channel_name:        bytes[channel_name_len]      # e.g. "mz", "intensity"
precision:           uint8        # matches Precision enum: 0=float32,
                                  # 1=float64, 2=int32, 3=int64,
                                  # 4=uint32, 5=complex128
compression:         uint8        # M79 Compression enum id: 0=none,
                                  # 1=zlib, 2=lz4, 3=numpress_delta,
                                  # 16=zstd (RFC 8878, single frame),
                                  # 17=float_delta_zstd (one FDZ1 stream)
n_elements:          uint32
data_length:         uint32       # compressed byte length
data:                bytes[data_length]

# MSImagePixel extension (spectrum_class == 4):
pixel_x:             uint32
pixel_y:             uint32
pixel_z:             uint32

# GenomicRead extension (spectrum_class == 5, M89.1):
chromosome_len:      uint16
chromosome:          bytes[chromosome_len]    # UTF-8, e.g. "chr1", "*", "chr22_KI270739v1_random"
position:            int64           # 0-based; -1 for unmapped (BAM convention)
mapping_quality:     uint8           # 0-255 (BAM convention)
flags:               uint16          # SAM/BAM bit flags
# M90.9 mate extension (optional — M89.1 payloads end after flags):
mate_position:       int64           # 0-based mate position; -1 if unpaired/unmapped
template_length:     int32           # observed template length (TLEN); 0 if unpaired
```

In addition to the M89.1 sequences/qualities UINT8 channels, M90.9
genomic AUs may carry per-AU UINT8 string channels named ``cigar``,
``read_name``, and ``mate_chromosome`` carrying the SAM-equivalent
per-read fields. Decoders dispatch on channel name; missing channels
default to the empty string. M101 adds a fourth string channel,
``tags``: the read's SAM optional fields, tab-joined as samtools
prints them. A sender lists ``tags`` in the DatasetHeader and carries
it in every AU only for a run that has a tags channel, so streams of
tag-less runs are unchanged. Under per-AU encryption the ``tags``
channel is encrypted with ``sequences`` and ``qualities``. The mate extension on the suffix is
optional — M89.1 fixtures decode unchanged with mate_position=-1 and
template_length=0.

The header fields (`retention_time` through `base_peak_intensity`)
constitute the **filter keys**: a server can evaluate any AU filter
(§7) against these bytes without touching `data`. For
GenomicRead AUs (`spectrum_class == 5`) the genomic suffix
(`chromosome`, `position`) is *also* a filter key — see §7 for
the genomic predicate set added in M89.3.

Channel payloads are conveyed in their native container encoding.
When `compression = 0` (`none`) the receiver decodes raw IEEE-754
values of the declared precision. When nonzero, the receiver MUST
apply the matching decoder (`zlib` / `lz4` / `numpress_delta` /
`zstd` / `float_delta_zstd`). `complex128` packs Re/Im as two
consecutive float64s.
For `zstd` (id 16) the payload is one standard RFC 8878 frame; the
plaintext size is `n_elements` times the element width, so readers
decode into an exact buffer.
For `float_delta_zstd` (id 17) the payload is one self-contained
FDZ1 stream exactly as the codec-17 on-disk format defines it
(magic, header, per-block transform byte + zstd frame of the
transposed byte planes). Its `n_values` MUST equal the channel's
`n_elements`; readers reject a mismatch. Encoders MAY differ
byte-wise across languages; decoders MUST accept any conforming
stream. On per-AU float64 payloads it is the smallest of the three
codecs (measured on real MS2 spectra: -6% vs zstd, -14% vs zlib).
Writers choose the codec through a selector (`compression_codec` on
the Python writer, `setCompressionCodec` on Java, `compressionCodec`
on ObjC); `float_delta_zstd` is the default when compression is
enabled, `zstd` (level 3) and `zlib` remain selectable. Readers older
than a codec's addition reject its id with their unsupported-
compression error (id 16: releases before 1.8.0; id 17: 1.8.0 and
older), so name `zlib` explicitly until a deployment's readers are
current.

#### 4.3.2 Encrypted-channel AU (`ENCRYPTED` set, v0.10+)

Same layout as §4.3.1, except each channel's `data` field carries
`[IV (12)] [TAG (16)] [ciphertext(plaintext-of-length n_elements ×
precision_size)]` instead of the raw plaintext bytes. `data_length
= 28 + ciphertext_bytes`.

The filter header (`retention_time` through `base_peak_intensity`)
remains plaintext; server-side filtering works keyless.

Each channel's AES-GCM operation uses **authenticated data**
`dataset_id (u16 LE) || au_sequence (u32 LE) || channel_name_utf8`.

MSImagePixel fields (`pixel_x / pixel_y / pixel_z`) stay
plaintext. GenomicRead fields (`chromosome / position /
mapping_quality / flags`) likewise stay plaintext under
`ENCRYPTED` only — server-side filtering on chromosome and
position remains keyless. (Encryption of the genomic suffix
under `ENCRYPTED_HEADER` lands in M89.5.)

#### 4.3.3 Fully-encrypted AU (`ENCRYPTED | ENCRYPTED_HEADER` set, v0.10+)

```
spectrum_class:      uint8               # plaintext — needed for dispatch
n_channels:          uint8               # plaintext — needed for parsing
IV_header:           bytes[12]
TAG_header:          bytes[16]
encrypted_semantic_header: bytes[36]
    # AES-GCM plaintext (36 bytes: 1+1+1+8+8+1+8+8) =
    #   acquisition_mode(u8) || ms_level(u8) || polarity(u8)
    #   || retention_time(f64) || precursor_mz(f64)
    #   || precursor_charge(u8) || ion_mobility(f64)
    #   || base_peak_intensity(f64)
    # AAD = dataset_id || au_sequence || "header"

# Repeated n_channels times (§4.3.2 layout):
(channel framing plaintext; data = IV || TAG || ciphertext)

# MSImagePixel extension (spectrum_class == 4):
IV_pixel:            bytes[12]
TAG_pixel:           bytes[16]
encrypted_pixel_xyz: bytes[12]
    # AES-GCM plaintext = pixel_x(u32) || pixel_y(u32) || pixel_z(u32)
    # AAD = dataset_id || au_sequence || "pixel"
```

Only `spectrum_class` and `n_channels` stay plaintext — the
minimum a reader needs to parse without the key. Server-side
filtering is **disabled** for streams carrying fully-encrypted
AUs; clients pull the whole stream and filter after decrypt.

#### 4.3.4 AAD summary

| Envelope                   | AAD                                                          |
|----------------------------|--------------------------------------------------------------|
| channel (§4.3.2, §4.3.3)   | `dataset_id || au_sequence || channel_name`                  |
| semantic header (§4.3.3)   | `dataset_id || au_sequence || "header"` (literal 6 bytes)    |
| pixel xyz (§4.3.3)         | `dataset_id || au_sequence || "pixel"` (literal 5 bytes)     |

Integrity check failures (bad tag) are fatal; the receiver
rejects the stream.

### 4.4 ProtectionMetadata (`0x04`)

Emitted once per stream before any encrypted AUs. Mirrors the
file-format KEK/DEK envelope defined in [`pqc.md`](pqc.md) §2 and
[`format-spec.md`](format-spec.md) §10b.

```
cipher_suite_len:    uint16
cipher_suite:        bytes[cipher_suite_len]      # e.g. "aes-256-gcm"
kek_algorithm_len:   uint16
kek_algorithm:       bytes[kek_algorithm_len]     # e.g. "ml-kem-1024"
wrapped_dek_len:     uint32
wrapped_dek:         bytes[wrapped_dek_len]
signature_algo_len:  uint16
signature_algorithm: bytes[signature_algo_len]    # e.g. "ml-dsa-87"
public_key_len:      uint32
public_key:          bytes[public_key_len]
```

Encrypted AUs carry the ciphertext of the channel payload; the
receiver materializes them into the file-format encryption layout
without decrypting in transit. Key management stays at the file
format layer.

### 4.5 Annotation (`0x05`)

```
record_kind:         uint8        # 0=identification, 1=quantification
record_json_len:     uint32
record_json:         bytes[record_json_len]       # JSON per format-spec.md §6
```

Annotations MAY be emitted at any point in the stream; they attach
to the most recent AU with the same `dataset_id` matching the
embedded `spectrum_index` JSON field (if present) or to the dataset
as a whole (if absent).

### 4.6 Provenance (`0x06`)

```
record_json_len:     uint32
record_json:         bytes[record_json_len]       # one ProvenanceRecord as JSON
```

Schema is defined in [`format-spec.md`](format-spec.md) §7. Multiple
Provenance packets MAY appear per dataset; they accumulate in
emission order.

### 4.7 Chromatogram (`0x07`)

```
chrom_id_len:        uint16
chrom_id:            bytes[chrom_id_len]          # e.g. "TIC", "BPC"
n_points:            uint32
time_precision:      uint8                        # 0=float32, 1=float64
intensity_precision: uint8
time_data_len:       uint32
time_data:           bytes[time_data_len]         # n_points * precision bytes
intensity_data_len:  uint32
intensity_data:      bytes[intensity_data_len]
```

No compression field — chromatograms are emitted raw. Use the
`flags & 0x02` (compressed) bit on the packet header if stream-level
compression of the whole payload is wanted.

### 4.8 EndOfDataset (`0x08`)

```
dataset_id:          uint16                       # matches header.dataset_id
final_au_sequence:   uint32                       # last au_sequence seen for this dataset
```

Allows the receiver to assert completeness. Emitted once per
DatasetHeader.

### 4.9 EndOfStream (`0xFF`)

Empty payload (`payload_length = 0`). MUST be the final packet.

### 4.10 BlobV2MateInfo (`0x09`) — bulk mode

Carries the `<run>/signal_channels/mate_info/inline_v2` codec blob
verbatim, plus the adjacent `chrom_names` table needed to reconstruct
mate references. Emitted only in bulk mode (`bulk_mode_v2_blobs`
feature flag, §6.4); per-AU streams MUST NOT emit this packet type.

```
dataset_id:           uint16            # matches header.dataset_id
codec_id:             uint8             # 13 = MATE_INLINE_V2
n_chrom_names:        uint16
chrom_names:          repeated { uint16 len, bytes[len] }
blob_length:          uint32
blob:                 bytes[blob_length]
```

When the receiver consumes this packet, it writes the blob bytes
verbatim to `signal_channels/mate_info/inline_v2` in the target
file (with `@compression = 13`) and the `chrom_names` table to the
adjacent compound dataset, BYPASSING the `MATE_INLINE_V2` codec
encode step. This is the only mechanism that preserves the SAM
sentinels `=` and `""` byte-for-byte across transport.

### 4.11 BlobV2RefDiff (`0x0A`) — bulk mode

Carries the `<run>/signal_channels/sequences/refdiff_v2` codec blob
verbatim. Emitted only when the source `.tio` has ref-diff-encoded
sequences (the codec is engaged when a reference is registered and
the alignment ratio passes the codec's heuristic).

```
dataset_id:           uint16            # matches header.dataset_id
codec_id:             uint8             # 14 = REF_DIFF_V2
reference_uri_len:    uint16
reference_uri:        bytes[len]
blob_length:          uint32
blob:                 bytes[blob_length]
```

The `reference_uri` MUST equal the run's `reference_uri` attribute.
A receiver decoding the blob MUST have the matching reference embedded
under `/study/references/<reference_uri>/`; otherwise it MUST reject
the stream.

### 4.12 BlobV2NameTok (`0x0B`) — bulk mode

Carries the `<run>/signal_channels/read_names` blob verbatim
(`@compression = 15`, NAME_TOKENIZED_V2). Emitted only in bulk mode.

```
dataset_id:           uint16            # matches header.dataset_id
codec_id:             uint8             # 15 = NAME_TOKENIZED_V2
blob_length:          uint32
blob:                 bytes[blob_length]
```

### 4.13 ReferenceGroupHeader (`0x10`) — v0.11

Declares a `ReferenceImport` in the dataset: an
`(reference_uri, chromosome_count, total_bases, md5_hex)` tuple.
Subsequent `ReferenceChromosome` packets carry the actual contig data,
terminated by `EndOfReferenceGroup`. One header per reference import.

```
uri_length:          uint16
uri_utf8:            bytes[uri_length]
chromosome_count:    uint32
total_bases:         uint64
md5_hex:             bytes[32]          # ASCII hex; matches format-spec §10.10 @md5 attr
```

### 4.14 ReferenceChromosome (`0x11`) — v0.11

One per contig within a reference group. Order MUST match the
sorted-name order used on disk (format-spec §10.10). The `encoding`
byte mirrors the format-spec's Perf-A decision (skip ZLIB below
4 KB) and lets the writer choose per chromosome; readers MUST
handle both.

```
name_length:         uint16
name_utf8:           bytes[name_length]
length:              uint64             # bases (also encoded as the length attribute on disk)
encoding:            uint8              # 0 = uncompressed UINT8, 1 = ZLIB-compressed
payload_length:      uint32
payload:             bytes[payload_length]   # raw bases or zlib stream
```

### 4.15 EndOfReferenceGroup (`0x12`) — v0.11

Terminator for a reference group; payload is a single
`chromosome_count_seen` field so the reader can assert against the
matching `ReferenceGroupHeader`.

```
chromosome_count_seen: uint32
```

### 4.16 ImageHeader (`0x13`) — v0.11

One per `MSImage` / vibrational / UV-Vis imaging cube (format-spec
§7, §7a, §7b). Declares grid + axis metadata. The `is_continuous`
flag selects which `ImagePixel` payload shape readers expect
(continuous-mode pixels carry intensities only; processed-mode
pixels add a per-pixel axis).

```
modality:            uint8              # 0=MS, 1=Raman, 2=IR, 3=UV-Vis
                                        # (matches AcquisitionMode ordinals for imaging modes)
width:               uint32
height:              uint32
spectrum_bins:       uint32             # per-pixel intensity samples
pixel_size_x:        float64
pixel_size_y:        float64
scan_pattern:        uint8              # 0=flyback, 1=meander, 2=random
axis_kind:           uint8              # 0=mz, 1=wavenumber, 2=wavelength, 3=ppm
axis_length:         uint32
axis:                float64[axis_length]   # shared axis (continuous mode)
                                            # or zeroes (processed mode)
is_continuous:       uint8              # 0/1; if 0, each pixel carries its own axis
title_length:        uint16
title_utf8:          bytes[title_length]
isa_id_length:       uint16
isa_id_utf8:         bytes[isa_id_length]
modality_extras_length: uint16          # length of modality-specific tail (may be 0)
modality_extras:     bytes[modality_extras_length]
                                        # modality-specific fields; layout
                                        # determined by `modality` (see below).
                                        # Readers that don't recognise the
                                        # modality MUST skip these bytes and
                                        # the following IMAGE_PIXEL+EOI block
                                        # rather than abort the stream.
```

**Modality-specific extras** (v0.11):

* `modality == 0` (MS): `modality_extras` is empty
  (`modality_extras_length == 0`).
* `modality == 1` (Raman):
  ```
  excitation_wavelength_nm: float64
  laser_power_mw:           float64
  ```
  (`modality_extras_length == 16`).
* `modality == 2` (IR):
  ```
  ir_mode:           uint8       # 0=transmittance, 1=absorbance
  resolution_cm_inv: float64
  ```
  (`modality_extras_length == 9`).
* `modality == 3` (UV-Vis): not defined in v0.11; readers MUST skip
  the IMAGE block. Future modalities populate this slot without
  breaking existing readers.

### 4.17 ImagePixel (`0x14`) — v0.11

One per pixel. The packet wraps a wire layout shared by both
emission modes (continuous and processed); the mode is fixed for the
entire image by `ImageHeader.is_continuous` and dictates how readers
parse the `payload_bytes` field below.

The fixed outer header (identical in both modes):
```
x:                   uint32
y:                   uint32
precision:           uint8              # 0=float32, 1=float64
compression:         uint8              # 0=none, 1=zstd, 2=zlib
payload_length:      uint32             # length of payload_bytes
payload_bytes:       bytes[payload_length]
```

**Continuous mode** (`ImageHeader.is_continuous == 1`):
`payload_bytes` is a dense array of `spectrum_bins` intensities at
the declared `precision` (FLOAT32 or FLOAT64), in `mzAxis` order.
No per-pixel axis — every pixel reuses the shared axis on the
`ImageHeader`.

**Processed mode** (`ImageHeader.is_continuous == 0`): each pixel
carries a sparse list of nonzero entries indexed into the shared
`mzAxis`. The MSImage data model remains a dense cube; readers
expand the sparse payload into the cube at decode time (zero-filling
unmentioned channels). `payload_bytes` carries:
```
nonzero_count:       uint32
entries:             nonzero_count × {
                       channel_index:  uint32  # 0 ≤ index < spectrum_bins
                       intensity:      bytes[precision_size]
                                               # FLOAT32 if precision=0
                                               # FLOAT64 if precision=1
                     }
```
NaN intensities, if produced, are emitted verbatim (no filtering)
so round-trips preserve them. Writers SHOULD pick processed mode
when it produces a smaller payload than continuous mode (a typical
heuristic: average nonzero count per pixel below
`spectrum_bins / 4`). Readers MUST tolerate both modes invisibly —
materialised MSImages look identical regardless of the mode chosen
on the wire.

### 4.18 EndOfImage (`0x15`) — v0.11

Terminator for an imaging cube; payload is a single
`pixel_count_seen` field for cross-check against `width * height`
declared in the matching `ImageHeader`.

```
pixel_count_seen:    uint32
```

### 4.19 IdentificationsTable (`0x16`) — v0.11

Single packet carrying the full `identifications` compound dataset
as a length-prefixed Apache Arrow IPC stream. The Arrow IPC stream
is the canonical Arrow inter-process format (magic `ARROW1` +
schema flatbuffer + record-batch flatbuffers + EOS); it carries its
own schema, dictionary encoding, and null bitmaps, so no TLV
envelope on top is needed. Writers MUST emit one IPC stream per
packet covering all rows of that table.

```
arrow_ipc_length:    uint32
arrow_ipc:           bytes[arrow_ipc_length]   # self-describing Arrow IPC stream
                                               # (schema message + record-batch messages)
```

### 4.20 QuantificationsTable (`0x17`) — v0.11

Identical wire shape to `IdentificationsTable` (Arrow IPC stream);
distinct packet type so receivers can dispatch without parsing the
payload first.

```
arrow_ipc_length:    uint32
arrow_ipc:           bytes[arrow_ipc_length]   # self-describing Arrow IPC stream
```

### 4.21 DatasetProvenance (`0x18`) — v0.11

Distinct from the per-run `Provenance` packet (0x06). Carries the
dataset-level provenance chain (format-spec §6.3). Per-record
structure mirrors the v0.10 per-run `Provenance` record layout
verbatim.

```
record_count:        uint32

# Per record (repeated record_count times):
timestamp_unix:      int64
software_length:     uint16
software:            bytes[software_length]
parameters_length:   uint16
parameters_json:     bytes[parameters_length]
input_refs_length:   uint16
input_refs_csv:      bytes[input_refs_length]
output_refs_length:  uint16
output_refs_csv:     bytes[output_refs_length]
```

### 4.22 SubjectMetadata (`0x19`) and SampleMetadata (`0x1A`) — v0.11

Two distinct packet types carrying the dataset's Subject and Sample
tables (see [`format-spec.md`](format-spec.md) §11 for the on-disk
HDF5 layout). Each packet wraps a single length-prefixed Apache
Arrow IPC stream — same wire shape as `IdentificationsTable` (§4.19)
and `QuantificationsTable` (§4.20). Writers MUST emit at most one
packet of each type per stream, covering all rows of that table.
Distinct packet types let receivers dispatch without parsing the
payload first; receivers ingesting these MUST write the rows into
the on-disk `/study/subjects/<external_id>/` and
`/study/samples/<sample_id>/` per-row groups.

```
arrow_ipc_length:    uint32
arrow_ipc:           bytes[arrow_ipc_length]   # self-describing Arrow IPC stream
                                               # (schema message + record-batch messages)
```

**`SUBJECT_METADATA (0x19)` Arrow schema:**

| Column            | Arrow type | Nullable | Notes                                                                                |
|-------------------|-----------:|---------:|--------------------------------------------------------------------------------------|
| `external_id`     | utf8       | no       | Primary key. Matches the on-disk group leaf name (`/study/subjects/<external_id>/`). |
| `project`         | utf8       | yes      | Empty-string in source → Arrow null on the wire.                                     |
| `sex`             | utf8       | yes      | Empty-string in source → Arrow null on the wire.                                     |
| `birth_year`      | int32      | yes      | Sentinel `0` in source → Arrow null on the wire. Widened from on-disk int64 to int32 for Arrow column-width consistency with `Identification.spectrum_index`; 4-digit years fit comfortably. |
| `attributes_json` | utf8       | yes      | Open-extension slot. Schema-nullable BUT always emitted with `"{}"` for empty (never literal Arrow null). Encoding matches the on-disk attribute (sort_keys, separators `","`/`":"`; see format-spec §11). |

**`SAMPLE_METADATA (0x1A)` Arrow schema:**

| Column                | Arrow type | Nullable | Notes                                                                                          |
|-----------------------|-----------:|---------:|------------------------------------------------------------------------------------------------|
| `sample_id`           | utf8       | no       | Primary key. Matches the on-disk group leaf name (`/study/samples/<sample_id>/`).              |
| `subject_external_id` | utf8       | yes      | Soft foreign key to `SUBJECT_METADATA.external_id`. Empty-string in source → Arrow null on the wire. Mismatch logs WARNING but is not an error (format-spec §11.3). |
| `sample_kind`         | utf8       | yes      | Empty-string in source → Arrow null on the wire.                                               |
| `collected_at`        | int64      | yes      | Unix seconds since epoch; sentinel `0` in source → Arrow null on the wire. Wider than `birth_year` intentionally — covers the full unix-seconds range. |
| `attributes_json`     | utf8       | yes      | Open-extension slot. Same `"{}"` convention as in `SUBJECT_METADATA`.                          |

**Null-handling convention** (applies symmetrically to encode + decode):

| Source representation                            | On-the-wire (Arrow)  | Materialised row                       |
|--------------------------------------------------|----------------------|----------------------------------------|
| Optional string column: empty-string (`""`)      | Arrow null           | empty-string (`""`)                    |
| Optional int column: sentinel `0`                | Arrow null           | sentinel `0`                           |
| `attributes_json`: always present                | utf8 (never null)    | `"{}"` if empty, sort_keys JSON otherwise |

Writers MUST emit Arrow null (not empty-string / not `0`) for absent
optional values so Python / Objective-C readers correctly produce
`None` / `nil` (or the documented sentinel in materialised rows).
The `attributes_json` column never carries literal Arrow null —
empty extension maps surface as the two-byte string `"{}"`.

**Cardinality.** Emit zero packets when the corresponding table is
empty (no subjects ⇒ no `SUBJECT_METADATA` packet; no samples ⇒ no
`SAMPLE_METADATA` packet). An empty table is never encoded as a
zero-row Arrow batch.

**Cross-language byte equivalence.** Logical equivalence, not byte
equality — same rule as `IdentificationsTable` / `QuantificationsTable`.
Arrow Java, pyarrow, and libarrow-C++ each emit different flatbuffer
envelopes; the row content cross-decodes correctly across all three
SDKs. The leading `uint32 arrow_ipc_length` prefix is byte-identical
across implementations.

**Ordering.** Per §5.4, when `transport_v0_11` is set, both
`SUBJECT_METADATA` and `SAMPLE_METADATA` appear in section 3 (after
`ENCRYPTION_ALGORITHM` + `DATASET_PROVENANCE`, before reference
groups). `SUBJECT_METADATA` MUST precede `SAMPLE_METADATA` so a
streaming receiver can resolve `subject_external_id` forward references
during materialise (a Sample whose `subject_external_id` has not yet
been seen still decodes — the FK is soft per format-spec §11.3 — but
emitting subjects first keeps the soft-FK warning surface clean).

The on-disk relationship to `AcquisitionRun.sampleName` (the canonical
run→sample link string) is unchanged by transport: see format-spec §11.

### 4.23 EncryptionAlgorithm (`0x1B`) — v0.11

Carries the dataset-level `@encrypted` algorithm name when present
so `isEncrypted()` round-trips correctly. Per-AU
`ProtectionMetadata` (0x04) continues to carry the per-key material;
this packet only conveys the algorithm-name string.

```
algorithm_length:    uint16
algorithm_utf8:      bytes[algorithm_length]
```

### 4.24 GenomicRunSidecar (`0x1C`) and BlockSidecar (`0x1D`) — v0.12

Carriage for per-AU-encrypted genomic runs in the `blocks_v1`
layout (format-spec §10.12). The per-read AU stream carries the
encrypted sequences and qualities plus the genomic filter suffix; the
sidecar packets carry everything else the receiver needs to rebuild a
container that `decrypt-in-place` restores byte-identically: the run
scalars and restore-policy attributes, the block index rows, and the
verbatim slices of the plaintext channel blobs (`read_names`,
`cigars`, `mate_info`).

When any genomic run in the stream uses `blocks_v1`, the sender adds
the required feature token `transport_blocks_v1` to
`StreamHeader.features`. The token is wire-scoped: receivers MUST NOT
write it into the output container's feature flags. Pre-v0.12
receivers skip the sidecar packets per §3.3 and rebuild a
legacy-shaped container that cannot be restored, so senders MUST NOT
target them with `blocks_v1` streams; the feature token is the
detection handle.

Exactly one GenomicRunSidecar follows the run's DatasetHeader, before
any of its BlockSidecars or AccessUnits:

```
layout:              uint16 len + UTF-8       # "blocks_v1"
block_policy:        uint16 len + UTF-8       # e.g. "reads=1000000,bytes=67108864"
read_count:          uint64
base_count:          uint64
attrs_json_len:      uint32
attrs_json:          bytes[attrs_json_len]    # JSON object; only keys present on
                                              # the run group are emitted:
                                              # ref_diff_slice_bytes (int),
                                              # opt_disable_qualities_v5 (int 1),
                                              # reference_md5s (string, verbatim)
channels_json_len:   uint32
channels_json:       bytes[channels_json_len] # JSON array, one object per plaintext
                                              # sidecar channel present in the run:
                                              # {"name", "compression" (dataset attr),
                                              #  "extra_attrs" {k: v} when present}
n_chrom_names:       uint32
chrom_names:         repeated uint16 len + UTF-8   # genomic_index/chromosome_names
                                                   # rows, in table order
n_mate_chrom_names:  uint32
mate_chrom_names:    repeated uint16 len + UTF-8   # mate_info/chrom_names rows
```

One BlockSidecar per block, in block order, `au_sequence` = block
index:

```
block_index:         uint32
read_start:          uint64
n_reads:             uint32
base_start:          uint64
n_bases:             uint64
n_channels:          uint8                    # the index's channels: the 5
                                              # required, plus "tags" when the
                                              # run's index has the M101 triple
# repeated n_channels times:
channel_name:        uint16 len + UTF-8
off:                 uint64                   # blocks/index <ch>_off
len:                 uint64                   # blocks/index <ch>_len
codec:               uint32                   # blocks/index <ch>_codec
n_blobs:             uint8                    # plaintext sidecar channels only
# repeated n_blobs times:
channel_name:        uint16 len + UTF-8
blob_length:         uint32
blob:                bytes[blob_length]       # the block's verbatim blob slice
```

The receiver writes the run attributes (including the verbatim
`reference_md5s` string and the policy attributes), recreates
`blocks/index` from the sidecar rows (with the `tags` triple when the
rows carry a `tags` entry, M101; the tags blob itself is encrypted and
travels as the per-AU `tags` channel, not as a sidecar blob), appends the blob slices into
the plaintext channel datasets with the carried dataset attributes,
builds the `genomic_index` arrays from the AU stream (lengths from
the segment table, positions, mapping qualities, flags and
chromosome ids from the filter suffixes, ids resolved against the
carried `chromosome_names` order), and writes both name tables from
the sidecar rather than from first-seen AU order, so the restored
tables match the sender's container.

The stream does not carry the embedded reference bytes. A received
REF_DIFF run restores through its `reference_md5s` attribute and
`REF_PATH` (format-spec §9.1.1); receivers that need the embedded
copy must obtain it out of band.

## 5. Ordering Rules

1. **StreamHeader first.** The first packet of every stream MUST be
   of type `0x01`.
2. **DatasetHeader before its AUs.** For every `dataset_id > 0`,
   the `DatasetHeader` packet MUST precede any `AccessUnit`,
   `Annotation`, `Provenance`, `Chromatogram`, or `EndOfDataset`
   carrying that `dataset_id`.
3. **ProtectionMetadata before encrypted AUs.** When any AU will
   carry `flags & 0x01` (`ENCRYPTED`) or `flags & 0x08`
   (`ENCRYPTED_HEADER`), a `ProtectionMetadata` packet for that
   `dataset_id` MUST precede it. `ENCRYPTED_HEADER` without
   `ENCRYPTED` is illegal (§3.1.1).
4. **Monotonic AU sequence per dataset.** Within one `dataset_id`,
   `au_sequence` MUST strictly increase. Gaps are permitted
   (e.g. after server-side filtering) but reordering is not.
5. **Interleaving across datasets is allowed.** AUs from different
   `dataset_id`s MAY be multiplexed. Receivers MUST buffer or route
   per-dataset without assuming global ordering.
6. **EndOfDataset per dataset.** Exactly one `EndOfDataset` per
   `DatasetHeader`, emitted after the last AU (or immediately after
   the header if the dataset is empty).
7. **EndOfStream last.** Exactly one `EndOfStream` as the final
   packet.
8. **Blob packets bracket their dataset.** `BlobV2MateInfo`,
   `BlobV2RefDiff`, and `BlobV2NameTok` MUST appear after the
   matching `DatasetHeader` and before the matching `EndOfDataset`
   for the same `dataset_id`. They MAY be interleaved freely with
   `AccessUnit` packets. At most one of each blob type per
   `dataset_id`; receivers MUST reject duplicates.

Receivers that encounter a violation SHOULD reject the stream and
SHOULD surface the violated rule to the application.

### 5.4 v0.11 ordering (after StreamHeader, before any DatasetHeader)

When `transport_v0_11` is in the StreamHeader features list, the
following sections appear in this order, each section optional:

1. Zero or more `ENCRYPTION_ALGORITHM` (0x1B) packets.
2. Zero or more `DATASET_PROVENANCE` (0x18) packets.
3. Zero or more `SUBJECT_METADATA` (0x19) and `SAMPLE_METADATA` (0x1A) packets.
4. Zero or more reference groups:
   `REFERENCE_GROUP_HEADER` → N × `REFERENCE_CHROMOSOME` → `END_OF_REFERENCE_GROUP`.
5. Zero or more image cubes:
   `IMAGE_HEADER` → N × `IMAGE_PIXEL` → `END_OF_IMAGE`.
6. Zero or more `IDENTIFICATIONS_TABLE` and `QUANTIFICATIONS_TABLE` packets.

The pre-existing v0.10 run sequences (DatasetHeader + AccessUnit + ...
+ EndOfDataset) follow. `EndOfStream` (0xFF) terminates the stream
regardless of section count.

## 6. Bidirectional Conversion

### 6.1 File → Transport

The writer iterates the source provider's `SpectrumIndex` in
ascending `(retention_time, au_sequence)` order for each
AcquisitionRun. For each spectrum it:

1. Reads the channel bytes via
   `provider.read_canonical_bytes(channel_path)` (see
   [`providers.md`](providers.md)) to preserve compressed/encrypted
   layouts byte-for-byte.
2. Copies the filter-key fields (RT, MS level, polarity, precursor
   m/z, ion mobility, base peak intensity) from the index without
   decoding signal data.
3. Emits an `AccessUnit` packet.

The writer MUST NOT re-compress or re-encode data that is already
in the target wire encoding. If the source uses a codec absent from
§4.3 (`compression` field), the writer either transcodes or
rejects the dataset; transcoding is reported through the
`ProvenanceRecord` chain.

### 6.2 Transport → File

The reader opens (or creates) a target provider, writes a
`DatasetHeader`-matching `AcquisitionRun` via the provider's
`StorageGroup` API, and appends each AU's channel bytes via
`storage_group.write_canonical_bytes`. Filter-key fields go into
the `SpectrumIndex`. Annotations, provenance, and chromatograms
route to the matching container groups.

The reader MUST NOT decode channel payloads unless the target
container demands a different encoding (e.g. an encrypted source
stream written into an unencrypted file on a different key).

### 6.3 Signal Round-Trip

File → transport → file round-trip is **bit-identical on signal
data** when source and target use the same storage backend and
codec. When backends differ (e.g. HDF5 → Zarr), signal values
match within float64 epsilon but byte layouts may differ.

Annotations, provenance, feature flags, and instrument config are
preserved verbatim. `ProtectionMetadata` is preserved for encrypted
datasets.

For **genomic** runs, both per-AU and bulk modes deliver the same
user-visible round-trip (the v2 mate codec normalizes
`mate_chromosome` sentinels `=` and `""` to resolved chrom names
and `*` respectively at write time, so neither mode "preserves"
the SAM-level sentinels — the source `.tio` no longer has them).
**Bulk mode** (§6.4) skips the receiver's v2 encode pass and
ships the v2 codec blobs verbatim, which preserves blob
byte-identity across deterministic codecs and avoids the encode
cost on the receiver.

### 6.4 Bulk Mode (genomic v2-blob carriage)

When the source `.tio` has v2-encoded genomic blobs on disk
(`mate_info/inline_v2`, optionally `sequences/refdiff_v2`,
`read_names` with `@compression = 15`), a writer MAY enable bulk
mode by:

1. Adding the `bulk_mode_v2_blobs` feature flag to the
   `StreamHeader.features` list. Flag has no `opt_` prefix → it is
   **required** per `feature-flags.md`. Receivers that cannot honor
   verbatim blob injection MUST refuse the stream.
2. Emitting one `BlobV2MateInfo`, one `BlobV2RefDiff` (if
   applicable), and one `BlobV2NameTok` per genomic dataset_id,
   between the `DatasetHeader` and `EndOfDataset` for that dataset
   (§5 rule 8).

Per-AU `AccessUnit` packets are still emitted in bulk mode. Their
chromosomes / positions / mapq / flags / lengths feed the
reconstructed `genomic_index/`; the per-AU mate_chromosome and
read_name string channels are **discarded** by the receiver in
favor of the verbatim blobs.

Bulk mode is post-acquisition only — live-streamed genomic AUs
have no v2 blobs to carry until block-flush. A writer SHOULD
degrade to per-AU mode (with a stderr warning) when bulk-mode is
requested but the source has no v2 blobs.

## 7. Selective Access

Servers supporting query-filtered streaming (see M68/M71) evaluate
client filters against AU-header fields BEFORE emitting the AU.
Supported filters:

| Filter                          | AU field                       | Comparison                    |
|---------------------------------|--------------------------------|-------------------------------|
| `rt_min` / `rt_max`             | `retention_time`               | closed interval               |
| `ms_level`                      | `ms_level`                     | equality                      |
| `precursor_mz_*`                | `precursor_mz`                 | closed interval               |
| `polarity`                      | `polarity`                     | equality                      |
| `dataset_id`                    | packet-header `dataset_id`     | equality                      |
| `max_au`                        | count                          | cap across all datasets       |
| `chromosome` (M89.3)            | genomic-suffix `chromosome`    | exact string equality         |
| `position_min` / `position_max` (M89.3) | genomic-suffix `position` | closed interval (inclusive)   |

In multiplexed streams (MS + genomic in one `.tis`, M89.4),
genomic predicates filter the AU types cleanly: `chromosome` and
`position_*` predicates exclude every non-genomic AU
(`spectrum_class != 5`), so a region query such as
`{"chromosome": "chr1", "position_min": 0, "position_max": 1e6}`
yields the genomic AUs in that window and skips MS entirely. To
combine modalities in one query, send two separate filtered
queries and union the resulting streams.

`StreamHeader`, `DatasetHeader`, `ProtectionMetadata`,
`EndOfDataset`, and `EndOfStream` packets are ALWAYS emitted
regardless of AU filters so the receiver has a complete container
skeleton.

Filtered streams are valid transport streams — a client MAY tee a
filtered stream to disk as a smaller `.tis` file, or materialize
it as a subset `.tio`.

## 8. Gotchas

(These entries supplement [`format-spec.md`](format-spec.md) §14.)

1. **Endianness.** All multi-byte fields are little-endian. Readers
   on big-endian platforms MUST byte-swap. Python: `struct.pack('<...')`.
   Java: `ByteBuffer.order(ByteOrder.LITTLE_ENDIAN)`. ObjC: explicit
   `OSSwapHostToLittleInt*` / explicit shift-and-OR.
2. **WebSocket frame splitting.** The transport format is
   frame-agnostic; a single AU MAY span multiple WebSocket
   continuation frames or, conversely, multiple small packets MAY
   share one frame. Receivers MUST buffer at the byte level, not
   the frame level.
3. **CRC-32C vs CRC-32.** The checksum is CRC-32C (Castagnoli), not
   the more common ISO 3309 CRC-32 used by zlib/gzip. Don't confuse
   the two.
4. **Multiplexed AUs.** Multi-run streams interleave AUs across
   `dataset_id`s. Receivers must route per `dataset_id` and must
   NOT assume AUs arrive in packed runs.
5. **Real-time `expected_au_count = 0`.** Live acquisitions signal
   unknown length by setting `expected_au_count = 0`. The receiver
   cannot pre-size indices and must grow them dynamically.
6. **Empty datasets are legal.** A DatasetHeader MAY be followed
   immediately by an EndOfDataset with no intervening AUs. This is
   the in-stream representation of an AcquisitionRun that was
   created but never had spectra written.
7. **ProtectionMetadata per dataset.** One dataset may be encrypted
   and another in the same stream unencrypted. The presence or
   absence of a ProtectionMetadata for a given `dataset_id` is
   authoritative; AU flags MUST match.

## 9. Conformance

A conforming transport implementation:

- Validates `magic`, `version`, and packet-type values on read;
  rejects malformed packets.
- Emits and validates CRC-32C when `has_checksum` is set.
- Honors ordering rules §5; rejects streams that violate them.
- Round-trips file → transport → file with bit-identical signal
  data on the same provider and codec.
- Round-trips file → transport → file with epsilon-equal signal
  data across providers.
- Preserves `ProtectionMetadata`, `Annotation`, `Provenance`, and
  `Chromatogram` records.

## 10. References

- [`format-spec.md`](format-spec.md) — file format specification
- [`providers.md`](providers.md) — storage provider protocol
- [`pqc.md`](pqc.md) — post-quantum crypto envelope model
- [`feature-flags.md`](feature-flags.md) — feature flag semantics
- [`class-hierarchy.md`](class-hierarchy.md) — spectrum class names
- ISO/IEC 23092-1 §7 — MPEG-G Transport Format (conceptual model)
- RFC 6455 — The WebSocket Protocol (M68 network layer)
