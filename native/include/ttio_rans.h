#ifndef TTIO_RANS_H
#define TTIO_RANS_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TTIO_RANS_L       (1u << 15)
#define TTIO_RANS_B_BITS  16
#define TTIO_RANS_B       (1u << 16)
#define TTIO_RANS_B_MASK  (TTIO_RANS_B - 1)
#define TTIO_RANS_T       (1u << 12)
#define TTIO_RANS_T_BITS  12
#define TTIO_RANS_T_MASK  (TTIO_RANS_T - 1)
#define TTIO_RANS_STREAMS 4
#define TTIO_RANS_X_MAX_PREFACTOR  ((TTIO_RANS_L >> TTIO_RANS_T_BITS) << TTIO_RANS_B_BITS)

#define TTIO_RANS_OK           0
#define TTIO_RANS_ERR_PARAM   -1
#define TTIO_RANS_ERR_ALLOC   -2
#define TTIO_RANS_ERR_CORRUPT -3
#define TTIO_RANS_ERR_RESERVED_MF        -4   /* mate_info_v2: MF value 3 seen */
#define TTIO_RANS_ERR_NS_LENGTH_MISMATCH -5   /* mate_info_v2: NS varint count != NUM_CROSS */
#define TTIO_RANS_ERR_ESC_LENGTH_MISMATCH -6  /* ref_diff_v2: ESC count mismatch */
#define TTIO_RANS_ERR_RESERVED_ESC_STREAM -7  /* ref_diff_v2: ESC stream_id >= 3 */
#define TTIO_RANS_ERR_NTV2_BAD_FLAG       -8  /* name_tok_v2: invalid 2-bit FLAG */
#define TTIO_RANS_ERR_NTV2_POOL_OOB       -9  /* name_tok_v2: pool_idx out of range */
#define TTIO_RANS_ERR_NTV2_BAD_K          -10 /* name_tok_v2: K=0 or K>=n_cols */
#define TTIO_RANS_ERR_NTV2_DICT_OVERFLOW  -11 /* name_tok_v2: dict code > dict size */
#define TTIO_RANS_ERR_NTV2_BAD_VERSION    -12 /* name_tok_v2: bad container version */
#define TTIO_RANS_ERR_NTV2_BAD_MAGIC      -13 /* name_tok_v2: magic != "NTK2" */

typedef struct ttio_rans_pool ttio_rans_pool;

int ttio_rans_encode_block(
    const uint8_t  *symbols,
    const uint16_t *contexts,
    size_t          n_symbols,
    uint16_t        n_contexts,
    const uint32_t (*freq)[256],
    uint8_t        *out,
    size_t         *out_len
);

int ttio_rans_decode_block(
    const uint8_t  *compressed,
    size_t          comp_len,
    const uint16_t *contexts,
    uint16_t        n_contexts,
    const uint32_t (*freq)[256],
    const uint32_t (*cum)[256],
    const uint8_t  (*dtab)[TTIO_RANS_T],
    uint8_t        *symbols,
    size_t          n_symbols
);

/*
 * Caller-provided context resolver.
 *
 * Called before decoding each symbol.  Receives:
 *   user_data : opaque pointer passed at decode time
 *   i         : current symbol index (0-based, in range [0, n_symbols))
 *   prev_sym  : the symbol just decoded (or 0 for i==0)
 * Returns the context ID for position i.
 *
 * Must be deterministic and side-effect free except for the caller's
 * own bookkeeping.  Must not return a context >= n_contexts; doing so
 * causes ttio_rans_decode_block_streaming to return TTIO_RANS_ERR_PARAM.
 */
typedef uint16_t (*ttio_rans_context_resolver)(
    void    *user_data,
    size_t   i,
    uint8_t  prev_sym
);

/*
 * Decode a block with on-the-fly context derivation.
 *
 * Same compressed-byte layout as ttio_rans_decode_block.  Calls
 * `resolver(user_data, i, prev_sym)` to obtain the context for each
 * position before decoding it.  Intended for codecs whose context
 * depends on previously decoded symbols (e.g. M94.Z order-1 cascades),
 * where the contexts[] array is unavailable up front.
 *
 * Note: the streaming decoder is scalar-only — it is bottlenecked by
 * the per-symbol callback, so SIMD acceleration would not help.
 */
int ttio_rans_decode_block_streaming(
    const uint8_t              *compressed,
    size_t                      comp_len,
    uint16_t                    n_contexts,
    const uint32_t            (*freq)[256],
    const uint32_t            (*cum)[256],
    const uint8_t             (*dtab)[TTIO_RANS_T],
    uint8_t                    *symbols,
    size_t                      n_symbols,
    ttio_rans_context_resolver  resolver,
    void                       *user_data
);

/*
 * M94.Z context-derivation parameters — mirrors the Python
 * ``ContextParams`` and Java ``ContextParams`` records used by the
 * pure-language reference implementations.  Defaults are
 * qbits=12, pbits=2, sloc=14 (see DEFAULT_QBITS / DEFAULT_PBITS /
 * DEFAULT_SLOC in fqzcomp_nx16_z.py and FqzcompNx16Z.java).
 */
typedef struct {
    uint32_t qbits;
    uint32_t pbits;
    uint32_t sloc;
} ttio_m94z_params;

/*
 * Decode a V2 block whose contexts follow the M94.Z scheme, with the
 * context derivation done inline in C.
 *
 * Replaces the per-symbol cross-language callback path of
 * `ttio_rans_decode_block_streaming` for the M94.Z codec — the JNI/
 * ctypes/objc round-trip per symbol made that approach slower than
 * the pure-language decoder.  By baking the (prev_q ring + position
 * bucket + revcomp) → context formula directly into C, the entire
 * decode loop runs without leaving native code.
 *
 * Inputs:
 *   compressed / comp_len  — same V2 byte layout as
 *                            `ttio_rans_decode_block`
 *   n_contexts             — number of dense contexts in freq/cum/dtab
 *   freq / cum / dtab      — per-DENSE-context frequency tables
 *   params                 — qbits / pbits / sloc (CRAM-Nx16 discipline)
 *   ctx_remap              — optional sparse→dense map of length
 *                            `1u << params->sloc`.  NULL means identity
 *                            (each sparse ctx == its dense index).  Any
 *                            sparse ctx not present in the active set
 *                            should map to `pad_ctx_dense` (typically 0).
 *   read_lengths           — uint32 lengths of each read, total ==
 *                            n_symbols
 *   n_reads                — number of entries in read_lengths /
 *                            revcomp_flags
 *   revcomp_flags          — 0/1 reverse-complement flag per read
 *   pad_ctx_dense          — dense ctx ID assigned to padding positions
 *                            (i >= n_symbols) and any sparse->dense miss
 *
 * Outputs:
 *   symbols  — decoded bytes, length n_symbols
 *
 * Returns TTIO_RANS_OK on success.
 */
int ttio_rans_decode_block_m94z(
    const uint8_t            *compressed,
    size_t                    comp_len,
    uint16_t                  n_contexts,
    const uint32_t          (*freq)[256],
    const uint32_t          (*cum)[256],
    const uint8_t           (*dtab)[TTIO_RANS_T],
    const ttio_m94z_params   *params,
    const uint16_t           *ctx_remap,
    const uint32_t           *read_lengths,
    size_t                    n_reads,
    const uint8_t            *revcomp_flags,
    uint16_t                  pad_ctx_dense,
    uint8_t                  *symbols,
    size_t                    n_symbols
);

int ttio_rans_build_decode_table(
    uint16_t        n_contexts,
    const uint32_t (*freq)[256],
    const uint32_t (*cum)[256],
    uint8_t        (*dtab)[TTIO_RANS_T]
);

/*
 * L2 (Task #82 Phase B.2, 2026-05-01): adaptive M94.Z (CRAM-mimic).
 *
 * Per-symbol adaptive freq updates: count[sym] += STEP, halve when
 * T > T_max - STEP. T_max = 65519, STEP = 16. Encoder maintains
 * the freq tables internally; decoder rebuilds them via the same
 * update rules. Wire format omits the freq-tables sidecar.
 */

#define TTIO_RANS_ADAPTIVE_STEP     16u
#define TTIO_RANS_ADAPTIVE_T_MAX    65519u

int ttio_rans_encode_block_adaptive(
    const uint8_t  *symbols,         /* n_symbols quality bytes      */
    const uint16_t *contexts,        /* n_symbols dense ctx indices  */
    size_t          n_symbols,
    uint16_t        n_contexts,      /* dense context count          */
    uint16_t        max_sym,         /* active range [0, max_sym)    */
    uint8_t        *out,             /* output buffer                */
    size_t         *out_len);        /* in: cap; out: actual         */

int ttio_rans_decode_block_adaptive_m94z(
    const uint8_t           *compressed,
    size_t                   comp_len,
    uint16_t                 n_contexts,
    uint16_t                 max_sym,
    const ttio_m94z_params  *params,
    const uint16_t          *ctx_remap,        /* len 1<<sloc        */
    const uint32_t          *read_lengths,
    size_t                   n_reads,
    const uint8_t           *revcomp_flags,
    uint16_t                 pad_ctx_dense,
    uint8_t                 *symbols,
    size_t                   n_symbols);

ttio_rans_pool *ttio_rans_pool_create(int n_threads);

int ttio_rans_encode_mt(
    ttio_rans_pool *pool,
    const uint8_t  *symbols,
    const uint16_t *contexts,
    size_t          n_symbols,
    uint16_t        n_contexts,
    size_t          reads_per_block,
    const size_t   *read_lengths,
    size_t          n_reads,
    uint8_t        *out,
    size_t         *out_len
);

int ttio_rans_decode_mt(
    ttio_rans_pool *pool,
    const uint8_t  *compressed,
    size_t          comp_len,
    uint8_t        *symbols,
    size_t         *n_symbols
);

void ttio_rans_pool_destroy(ttio_rans_pool *pool);

/* Diagnostic: name of the kernel selected at library-load time by cpuid
 * dispatch.  Returns a pointer to a static string — one of
 * "scalar", "sse4.1", "avx2".  Never returns NULL. */
const char *ttio_rans_kernel_name(void);

/* M94.Z V4: CRAM 3.1 fqzcomp port. See native/src/m94z_v4_wire.h
 * + native/src/fqzcomp_qual.h for details. The V4 outer wire format
 * wraps a CRAM-byte-compatible fqzcomp body with an M94.Z header
 * (magic "M94Z", version=4) so codec layers can dispatch on version.
 *
 * Encode:
 *   qual_in        — n_qualities bytes (Phred-33 ASCII)
 *   read_lengths   — n_reads uint32 (sum must equal n_qualities)
 *   flags          — n_reads bytes (bit 4 = SAM_REVERSE)
 *   strategy_hint  — -1 = auto-tune, 0..4 = preset
 *   pad_count      — 0..3 (V3 pad-count convention, packed in flags)
 *   out, *out_len  — caller-owned buffer + capacity-in/length-out
 *
 * Decode:
 *   read_lengths   — caller-allocated, n_reads entries; populated
 *                    from the V4 header's deflated RLT
 *   flags          — n_reads bytes (must match those used at encode)
 *   out_qual       — caller-owned, n_qualities bytes
 */
int ttio_m94z_v4_encode(
    const uint8_t  *qual_in,
    size_t          n_qualities,
    const uint32_t *read_lengths,
    size_t          n_reads,
    const uint8_t  *flags,
    int             strategy_hint,
    uint8_t         pad_count,
    uint8_t        *out,
    size_t         *out_len);

int ttio_m94z_v4_decode(
    const uint8_t  *in,
    size_t          in_len,
    uint32_t       *read_lengths,
    size_t          n_reads,
    const uint8_t  *flags,
    uint8_t        *out_qual,
    size_t          n_qualities);

/* Qualities V5: sequence-context strategies. The
 * umbrella auto-tunes across the V4 presets plus S5/S6 and keeps the
 * smallest stream by exact size; ties go to V4. S5/S6 are tried only
 * when seq_in is non-NULL and n_qualities is at least the floor
 * below; forced hints 5/6 bypass the floor but require seq_in.
 *
 * Encode:
 *   seq_in         — n_qualities base bytes parallel to qual_in, or
 *                    NULL for V4-only behaviour
 *   strategy_hint  — -1 auto, 0..4 V4 preset, 5..6 forced sequence
 *                    strategy, 7 V4 with internal preset selection
 *                    (the no-seq auto path, sequences ignored)
 * Decode dispatches on the version byte: a version-5 stream requires
 * seq_in (the decoded sequences channel) and fails with a distinct
 * error without it; version-4 streams ignore seq_in. */
#define TTIO_M94Z_V5_VERSION 5
#define TTIO_M94Z_V5_MIN_QUALITIES (1u << 20)
#define TTIO_M94Z_HINT_V4_AUTO 7

/* 8 forces M94.Z V6, the segmented adaptive variant. V6 needs no
 * sequences, on either side. Auto-tune never selects it: V6 does not
 * beat V4 or V5 on size, so it must not enter the size race. */
#define TTIO_M94Z_HINT_V6 8

/* Name of the engine that encodes M94.Z V6 blocks: "cpu", or
 * "vulkan:<device>" when a GPU engine has been asked for and came up.
 * Never NULL. */
const char *ttio_engine_active_name(void);

/* 1 when a GPU engine has been asked for and came up, else 0. */
int ttio_engine_gpu_available(void);

/* Strategy of an encoded M94.Z stream: 4 = V4, 5/6 = V5 S5/S6,
 * 8 = V6. <0: -1 args, -2 magic/version, -3 truncated or unknown id. */
int ttio_m94z_qual_stream_strategy(const uint8_t *in, size_t in_len);

/* Threads the FQZCOMP auto-tune uses for its V4/S5/S6 candidate encodes
 * (default 3; <= 1 runs them in sequence). Process-global; the initial
 * value is 1 when TTIO_M94Z_SEQUENTIAL=1 is set. A caller that already
 * runs blocks on a pool sets 1 for the life of the pool. */
void ttio_m94z_set_autotune_threads(int n);

/* Segments of one V6 block encoded at once.
 *
 * This is a different question from the auto-tune candidate count, and
 * it wants a different answer. Auto-tune races V4 against two V5
 * strategies, so a caller that already runs blocks on a pool sets it to
 * 1 to avoid three candidate encodes per worker. V6 has no candidates:
 * the same setting would simply switch off intra-block parallelism and
 * leave the segments to encode one after another.
 *
 * 0, the default, means follow the auto-tune count, which is the
 * behaviour callers had before this existed. */
void ttio_m94z_set_v6_threads(int n);
int  ttio_m94z_get_v6_threads(void);

/* Width of M94.Z V6's sequence-context field, for streams this
 * process writes. 0, the default, is the context V6 shipped with and
 * needs no sequences to decode. 255 asks the encoder to choose per
 * block. Any other value is used as given, and encoding without
 * sequences then fails rather than quietly dropping the field.
 *
 * The width travels in the stream, so decoding never consults this. */
void ttio_m94z_set_v6_sbits(int n);
int  ttio_m94z_get_v6_sbits(void);
int  ttio_m94z_get_autotune_threads(void);

int ttio_m94z_qual_encode(
    const uint8_t  *qual_in,
    size_t          n_qualities,
    const uint32_t *read_lengths,
    size_t          n_reads,
    const uint8_t  *flags,
    const uint8_t  *seq_in,
    int             strategy_hint,
    uint8_t         pad_count,
    uint8_t        *out,
    size_t         *out_len);

int ttio_m94z_qual_decode(
    const uint8_t  *in,
    size_t          in_len,
    uint32_t       *read_lengths,
    size_t          n_reads,
    const uint8_t  *flags,
    const uint8_t  *seq_in,
    uint8_t        *out_qual,
    size_t          n_qualities);

/* ──────────────────────────────────────────────────────────────────────
 * Plain rANS-O0 — byte-exact port of python/src/ttio/codecs/rans.py
 * (order=0 path) and Java/ObjC equivalents. Used by mate_info v2
 * substreams (NS / NP / TS / MF auto-pick), REF_DIFF v2, and
 * NameTokenized v2.
 *
 * Wire format: 9-byte header + 1024-byte freq table + payload.
 *   [order=0x00][orig_len u32 BE][payload_len u32 BE]
 *   [256 × u32 BE freq table]
 *   [4 × u8 final state BE][renorm bytes — read forward at decode]
 *
 * Algorithm: M=4096, L=2^23, b=256, state width 64 bits.
 * ────────────────────────────────────────────────────────────────────── */
size_t ttio_rans_o0_max_encoded_size(size_t in_len);

int ttio_rans_o0_encode(
    const uint8_t *in,
    size_t         in_len,
    uint8_t       *out,
    size_t        *out_len);   /* in: capacity, out: actual size */

int ttio_rans_o0_decode(
    const uint8_t *in,
    size_t         in_len,
    uint8_t       *out,
    size_t         out_capacity,
    size_t        *out_len);   /* writes original-length on success */

/* ──────────────────────────────────────────────────────────────────────
 * mate_info v2 — CRAM-style inline mate-pair encoding.
 *
 * Encode produces a self-contained uint8 blob written as
 * signal_channels/mate_info/inline_v2 with @compression = 13.
 *
 * Encode inputs:
 *   mate_chrom_ids[N]    — int32, -1 if RNEXT='*'
 *   mate_positions[N]    — int64, 0-based POS
 *   template_lengths[N]  — int32, signed tlen
 *   own_chrom_ids[N]     — uint16 from genomic_index/chromosome_ids;
 *                          0xFFFF treated as -1
 *   own_positions[N]     — int64 from genomic_index/positions
 *
 * Decode requires (own_chrom_ids, own_positions, n_records) to
 * reconstruct mate_pos for SAME_CHROM records.
 *
 * Returns 0 on success; negative TTIO_RANS_ERR_* on framing,
 * length, or reserved-value violations.
 * ────────────────────────────────────────────────────────────────────── */
size_t ttio_mate_info_v2_max_encoded_size(uint64_t n_records);

int ttio_mate_info_v2_encode(
    const int32_t  *mate_chrom_ids,
    const int64_t  *mate_positions,
    const int32_t  *template_lengths,
    const uint16_t *own_chrom_ids,
    const int64_t  *own_positions,
    uint64_t        n_records,
    uint8_t        *out,
    size_t         *out_len);  /* in: capacity, out: actual size */

int ttio_mate_info_v2_decode(
    const uint8_t  *encoded,
    size_t          encoded_size,
    const uint16_t *own_chrom_ids,
    const int64_t  *own_positions,
    uint64_t        n_records,
    int32_t        *out_mate_chrom_ids,
    int64_t        *out_mate_positions,
    int32_t        *out_template_lengths);

/* ──────────────────────────────────────────────────────────────────────
 * REF_DIFF v2 — CRAM-style bit-packed sequence diff codec (codec id 14).
 *
 * Encoded blob written as signal_channels/sequences/refdiff_v2 with
 * @compression = 14. Outer container preserves v1 slice index;
 * each slice body is a 24-byte sub-header + 5 rANS-O0-encoded
 * substreams (FLAG / BS / IN / SC / ESC).
 *
 * Returns 0 on success; negative TTIO_RANS_ERR_* on framing,
 * length, or reserved-value violations.
 * ────────────────────────────────────────────────────────────────────── */
typedef struct {
    const uint8_t  *sequences;        /* concatenated ACGTN bytes */
    const uint64_t *offsets;          /* n_reads + 1 entries */
    const int64_t  *positions;        /* 1-based reference position */
    const char    **cigar_strings;    /* per-read CIGAR */
    uint64_t        n_reads;
    const uint8_t  *reference;        /* reference chrom bytes */
    uint64_t        reference_length;
    uint64_t        reads_per_slice;  /* default 10000 */
    const uint8_t  *reference_md5;   /* 16 bytes */
    const char     *reference_uri;   /* UTF-8 nul-terminated */
    uint64_t        slice_bytes;      /* 0 = the reads_per_slice rule.
                                       * > 0: a slice closes before the
                                       * read that would push it past
                                       * this many bases (each slice
                                       * keeps >= 1 read; reads_per_slice
                                       * still caps the read count).
                                       * Writer policy only — wire format
                                       * and decoder are unchanged. */
} ttio_ref_diff_v2_input;

size_t ttio_ref_diff_v2_max_encoded_size(uint64_t n_reads, uint64_t total_bases);

/* Capacity bound honouring a byte-budget slice policy (slice_bytes > 0
 * can produce more slices than the reads_per_slice rule would).
 * slice_bytes == 0 matches ttio_ref_diff_v2_max_encoded_size. */
size_t ttio_ref_diff_v2_max_encoded_size2(uint64_t n_reads, uint64_t total_bases,
                                          uint64_t slice_bytes);

int ttio_ref_diff_v2_encode(
    const ttio_ref_diff_v2_input *in,
    uint8_t *out,
    size_t  *out_len);

int ttio_ref_diff_v2_decode(
    const uint8_t  *encoded,
    size_t          encoded_size,
    const int64_t  *positions,
    const char    **cigar_strings,
    uint64_t        n_reads,
    const uint8_t  *reference,
    uint64_t        reference_length,
    uint8_t        *out_sequences,
    uint64_t       *out_offsets);

/* ──────────────────────────────────────────────────────────────────────
 * NAME_TOKENIZED v2 — multi-substream + DUP-pool + PREFIX-MATCH codec
 * (codec id 15).
 *
 * Encoded blob written to read_names HDF5 dataset with @compression = 15.
 * Wire magic "NTK2", version 0x01.
 *
 * Returns 0 on success; negative TTIO_RANS_ERR_* on framing or value
 * violations.
 * ────────────────────────────────────────────────────────────────────── */

size_t ttio_name_tok_v2_max_encoded_size(uint64_t n_reads, uint64_t total_name_bytes);

/* Encodes n_reads names. names[i] is a NUL-terminated 7-bit-ASCII string.
 * Caller allocates `out` of at least max_encoded_size bytes; *out_len
 * is set to the actual encoded size. */
int ttio_name_tok_v2_encode(
    const char * const *names,
    uint64_t            n_reads,
    uint8_t            *out,
    size_t             *out_len);

/* Decodes a v2 stream. *out_names is malloc'd as an array of n_reads
 * c-string pointers (each entry malloc'd separately); caller frees
 * each entry plus the array. *out_n_reads is set. */
int ttio_name_tok_v2_decode(
    const uint8_t  *encoded,
    size_t          encoded_size,
    char         ***out_names,
    uint64_t       *out_n_reads);

/* ──────────────────────────────────────────────────────────────────────
 * SAM_TAGS — SAM optional fields as a tag-line dictionary plus one
 * column per (key, type, kind), MD/NM recomputed from the reference
 * (codec id 18, M101). Spec: docs/codecs/sam_tags.md.
 *
 * The blob is written as signal_channels/tags with @compression = 18.
 * Wire magic "STG1", version 0x01.
 *
 * Tag text per read is SAM columns 12+ tab-joined, as samtools prints
 * them; tags[tag_offsets[i] .. tag_offsets[i+1]) is read i's text.
 * Encode and decode take the same context; the reads' sequences,
 * CIGARs, positions and chromosome ids, and the reference bases, are
 * needed only for MD/NM derivation (n_refs == 0 disables it, and the
 * other context arrays may then be NULL).
 *
 * Returns 0 on success; TTIO_RANS_ERR_PARAM on bad input (a NUL in the
 * tag text included), TTIO_RANS_ERR_ALLOC, or TTIO_RANS_ERR_CORRUPT on
 * a blob that violates the spec.
 * ────────────────────────────────────────────────────────────────────── */
typedef struct {
    uint64_t        n_reads;
    const uint8_t  *sequences;      /* concatenated SEQ bytes          */
    const uint64_t *seq_offsets;    /* n_reads + 1                     */
    const uint8_t  *cigars;         /* concatenated CIGAR text         */
    const uint64_t *cigar_offsets;  /* n_reads + 1                     */
    const int64_t  *positions;      /* 1-based POS, 0 = unmapped       */
    const uint16_t *chrom_ids;      /* 0xFFFF = none                   */
    const uint8_t * const *refs;    /* n_refs chromosome sequences; NULL entries allowed */
    const uint64_t *ref_lengths;
    uint32_t        n_refs;         /* 0 disables MD/NM derivation     */
} ttio_sam_tags_ctx;

/* *out is allocated by the library; release it with ttio_sam_tags_free. */
int ttio_sam_tags_encode(
    const ttio_sam_tags_ctx *ctx,
    const uint8_t  *tags,
    const uint64_t *tag_offsets,     /* n_reads + 1 */
    uint8_t       **out,
    size_t         *out_len);

/* *out_tags is allocated by the library (release with ttio_sam_tags_free);
 * out_tag_offsets is caller-allocated with n_reads + 1 entries. */
int ttio_sam_tags_decode(
    const ttio_sam_tags_ctx *ctx,
    const uint8_t  *encoded,
    size_t          encoded_len,
    uint8_t       **out_tags,
    uint64_t       *out_tag_offsets);

void ttio_sam_tags_free(void *p);

/* ──────────────────────────────────────────────────────────────────────
 * SEQ_CM — read bases coded with a reference-free context-mixing model
 * (up to three base-context lengths, a logistic mixer and an optional
 * reverse-complement training pass) and a binary arithmetic coder
 * (codec id 19, M103). Spec: docs/codecs/seq_cm.md.
 *
 * The blob is written as signal_channels/sequences with @compression = 19.
 * Wire magic "SQC1", version 0x01.
 *
 * seq holds the reads' bases back to back; lengths[i] is read i's length.
 * Decode takes the same lengths. Bytes other than A/C/G/T (N, IUPAC,
 * lower case) round-trip as exception runs. params == NULL takes the
 * defaults of ttio_seq_cm_default_params; table_bits == 0 lets the
 * encoder choose from the base count. The model holds 16 << table_bits
 * bytes per hashed order, on encode and on decode.
 *
 * Returns 0 on success; TTIO_RANS_ERR_PARAM on bad input (lengths that
 * do not match the blob included), TTIO_RANS_ERR_ALLOC, or
 * TTIO_RANS_ERR_CORRUPT on a blob that violates the spec.
 * ────────────────────────────────────────────────────────────────────── */
#define TTIO_SEQ_CM_FLAG_RC 0x01    /* train on each read's reverse complement */

typedef struct {
    uint8_t  n_orders;              /* 1 .. 3                              */
    uint8_t  orders[3];             /* context lengths in bases, rising    */
    uint8_t  table_bits;            /* 10 .. 30, or 0 for automatic        */
    uint8_t  flags;                 /* TTIO_SEQ_CM_FLAG_*                  */
    uint16_t limit;                 /* counter adaptation limit, 1 .. 1023 */
    uint16_t lr;                    /* mixer learning rate, 1 .. 4096      */
} ttio_seq_cm_params;

void ttio_seq_cm_default_params(ttio_seq_cm_params *p);

/* *out is allocated by the library; release it with ttio_seq_cm_free. */
int ttio_seq_cm_encode(
    const uint8_t  *seq,
    const uint64_t *lengths,
    uint64_t        n_reads,
    const ttio_seq_cm_params *params,
    uint8_t       **out,
    size_t         *out_len);

/* *out_seq is allocated by the library (release with ttio_seq_cm_free). */
int ttio_seq_cm_decode(
    const uint8_t  *encoded,
    size_t          encoded_len,
    const uint64_t *lengths,
    uint64_t        n_reads,
    uint8_t       **out_seq,
    size_t         *out_len);

void ttio_seq_cm_free(void *p);

#ifdef __cplusplus
}
#endif

#endif /* TTIO_RANS_H */
