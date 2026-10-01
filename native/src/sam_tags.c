/*
 * sam_tags.c — SAM_TAGS codec (codec id 18): SAM optional fields as a
 * tag-line dictionary plus one column per (key, type, kind), with MD/NM
 * recomputed from the reference and repeated integers stored as
 * back-references. Spec: docs/codecs/sam_tags.md.
 */
#include "sam_tags.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ──────── growable byte buffer ──────── */

typedef struct {
    uint8_t *p;
    size_t len, cap;
} stg_buf;

static int buf_reserve(stg_buf *b, size_t extra) {
    if (b->len + extra <= b->cap) return 0;
    size_t cap = b->cap ? b->cap : 256;
    while (cap < b->len + extra) cap *= 2;
    uint8_t *np = (uint8_t *)realloc(b->p, cap);
    if (!np) return -1;
    b->p = np;
    b->cap = cap;
    return 0;
}

static int buf_put(stg_buf *b, const void *src, size_t n) {
    if (n == 0) return 0;
    if (buf_reserve(b, n)) return -1;
    memcpy(b->p + b->len, src, n);
    b->len += n;
    return 0;
}

static int buf_u8(stg_buf *b, uint8_t v) { return buf_put(b, &v, 1); }

static int buf_uvarint(stg_buf *b, uint64_t v) {
    uint8_t tmp[10];
    size_t n = 0;
    while (v >= 0x80) {
        tmp[n++] = (uint8_t)(v | 0x80);
        v >>= 7;
    }
    tmp[n++] = (uint8_t)v;
    return buf_put(b, tmp, n);
}

static void buf_free(stg_buf *b) {
    free(b->p);
    b->p = NULL;
    b->len = b->cap = 0;
}

/* ──────── reader over an input span ──────── */

typedef struct {
    const uint8_t *p;
    size_t len, pos;
} stg_rd;

static int rd_u8(stg_rd *r, uint8_t *v) {
    if (r->pos >= r->len) return -1;
    *v = r->p[r->pos++];
    return 0;
}

static int rd_uvarint(stg_rd *r, uint64_t *v) {
    uint64_t x = 0;
    for (int shift = 0; shift < 64; shift += 7) {
        uint8_t b;
        if (rd_u8(r, &b)) return -1;
        x |= (uint64_t)(b & 0x7F) << shift;
        if (b < 0x80) {
            *v = x;
            return 0;
        }
    }
    return -1;
}

static int rd_bytes(stg_rd *r, size_t n, const uint8_t **out) {
    if (n > r->len - r->pos) return -1;
    *out = r->p + r->pos;
    r->pos += n;
    return 0;
}

static uint64_t zz_enc(int64_t v) { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
static int64_t zz_dec(uint64_t u) { return (int64_t)(u >> 1) ^ -(int64_t)(u & 1); }

/* ──────── substreams (spec §2.1) ──────── */

/* RANS_ORDER0 (codec 4) blob layout: order u8, orig_len u32 BE,
 * payload_len u32 BE, 256 × u32 BE frequencies, payload. */
#define O0_HDR 9
#define O0_TABLE 1024

static uint32_t be32(const uint8_t *p) {
    return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | p[3];
}
static void put_be32(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)(v >> 24); p[1] = (uint8_t)(v >> 16); p[2] = (uint8_t)(v >> 8); p[3] = (uint8_t)v;
}

/* Mode 1 body (spec §2.1): the codec-4 stream with its frequency table
 * listed sparsely — uvarint n_syms, n_syms × { sym u8, freq uvarint }
 * in ascending symbol order — followed by the codec-4 payload. */
static int put_substream(stg_buf *out, const uint8_t *s, size_t n) {
    stg_buf body = {0};
    int use_rans = 0;
    if (n > 0 && n <= UINT32_MAX) {
        size_t cap = ttio_rans_o0_max_encoded_size(n);
        uint8_t *rb = (uint8_t *)malloc(cap);
        if (!rb) return -1;
        size_t rl = cap;
        if (ttio_rans_o0_encode(s, n, rb, &rl) == 0 && rl >= O0_HDR + O0_TABLE) {
            const uint8_t *tab = rb + O0_HDR;
            uint32_t payload_len = be32(rb + 5);
            unsigned n_syms = 0;
            for (int k = 0; k < 256; k++) n_syms += be32(tab + 4 * k) != 0;
            int rc = buf_uvarint(&body, n_syms);
            for (int k = 0; k < 256 && !rc; k++) {
                uint32_t f = be32(tab + 4 * k);
                if (f) rc = buf_u8(&body, (uint8_t)k) || buf_uvarint(&body, f);
            }
            if (!rc) rc = buf_put(&body, rb + O0_HDR + O0_TABLE, payload_len);
            if (rc) { free(rb); buf_free(&body); return -1; }
            use_rans = body.len < n;
        }
        free(rb);
    }
    int rc;
    if (use_rans)
        rc = buf_u8(out, 1) || buf_uvarint(out, n) || buf_uvarint(out, body.len) ||
             buf_put(out, body.p, body.len);
    else
        rc = buf_u8(out, 0) || buf_uvarint(out, n) || buf_uvarint(out, n) || buf_put(out, s, n);
    buf_free(&body);
    return rc ? -1 : 0;
}

/* Rebuilds the codec-4 blob from a mode-1 body and decodes it. */
static int rans_compact_decode(const uint8_t *body, size_t body_len, uint8_t *dst, size_t raw_len) {
    if (raw_len == 0 || raw_len > UINT32_MAX) return TTIO_RANS_ERR_CORRUPT;
    stg_rd r = {body, body_len, 0};
    uint64_t n_syms;
    if (rd_uvarint(&r, &n_syms) || n_syms == 0 || n_syms > 256) return TTIO_RANS_ERR_CORRUPT;
    uint8_t *blob = (uint8_t *)calloc(1, O0_HDR + O0_TABLE + body_len);
    if (!blob) return TTIO_RANS_ERR_ALLOC;
    int prev = -1;
    for (uint64_t k = 0; k < n_syms; k++) {
        uint8_t sym;
        uint64_t f;
        if (rd_u8(&r, &sym) || (int)sym <= prev || rd_uvarint(&r, &f) || f == 0 || f > 4096) {
            free(blob);
            return TTIO_RANS_ERR_CORRUPT;
        }
        put_be32(blob + O0_HDR + 4 * sym, (uint32_t)f);
        prev = sym;
    }
    size_t payload_len = body_len - r.pos;
    blob[0] = 0;
    put_be32(blob + 1, (uint32_t)raw_len);
    put_be32(blob + 5, (uint32_t)payload_len);
    memcpy(blob + O0_HDR + O0_TABLE, body + r.pos, payload_len);
    size_t got = 0;
    int rc = ttio_rans_o0_decode(blob, O0_HDR + O0_TABLE + payload_len, dst, raw_len, &got);
    free(blob);
    return (rc == 0 && got == raw_len) ? 0 : TTIO_RANS_ERR_CORRUPT;
}

/* Decodes one substream into a malloc'd buffer (*out, *out_len). */
static int get_substream(stg_rd *r, uint8_t **out, size_t *out_len) {
    uint8_t mode;
    uint64_t raw_len, body_len;
    const uint8_t *body;
    if (rd_u8(r, &mode) || rd_uvarint(r, &raw_len) || rd_uvarint(r, &body_len) ||
        rd_bytes(r, (size_t)body_len, &body))
        return TTIO_RANS_ERR_CORRUPT;
    if (raw_len > UINT32_MAX) return TTIO_RANS_ERR_CORRUPT;
    uint8_t *buf = (uint8_t *)malloc(raw_len ? (size_t)raw_len : 1);
    if (!buf) return TTIO_RANS_ERR_ALLOC;
    if (mode == 0) {
        if (body_len != raw_len) { free(buf); return TTIO_RANS_ERR_CORRUPT; }
        memcpy(buf, body, (size_t)raw_len);
    } else if (mode == 1) {
        int rc = rans_compact_decode(body, (size_t)body_len, buf, (size_t)raw_len);
        if (rc) {
            free(buf);
            return rc;
        }
    } else {
        free(buf);
        return TTIO_RANS_ERR_CORRUPT;
    }
    *out = buf;
    *out_len = (size_t)raw_len;
    return 0;
}

/* ──────── canonical integers, MD/NM ──────── */

int stg_parse_canonical_int(const uint8_t *s, size_t n, int64_t *out) {
    size_t i = 0;
    int neg = 0;
    if (n == 0) return 0;
    if (s[0] == '-') {
        neg = 1;
        i = 1;
        if (n == 1) return 0;
    }
    if (s[i] == '0') {
        if (n == 1) { *out = 0; return 1; }
        return 0; /* leading zero, or "-0" */
    }
    uint64_t mag = 0;
    for (; i < n; i++) {
        if (s[i] < '0' || s[i] > '9') return 0;
        uint64_t d = (uint64_t)(s[i] - '0');
        if (mag > (UINT64_MAX - d) / 10) return 0;
        mag = mag * 10 + d;
    }
    if (neg) {
        if (mag > (uint64_t)INT64_MAX + 1) return 0;
        *out = (int64_t)(0 - mag);
    } else {
        if (mag > (uint64_t)INT64_MAX) return 0;
        *out = (int64_t)mag;
    }
    return 1;
}

static uint8_t up(uint8_t b) { return (b >= 'a' && b <= 'z') ? (uint8_t)(b - 32) : b; }

static int put_run(stg_buf *b, uint64_t run) {
    char tmp[24];
    int n = snprintf(tmp, sizeof tmp, "%" PRIu64, run);
    return buf_put(b, tmp, (size_t)n);
}

int stg_calc_md_nm(const ttio_sam_tags_ctx *ctx, uint64_t i,
                   char **md, size_t *md_len, int64_t *nm) {
    if (ctx->n_refs == 0 || !ctx->chrom_ids || !ctx->positions || !ctx->sequences ||
        !ctx->seq_offsets || !ctx->cigars || !ctx->cigar_offsets)
        return 0;
    uint16_t chrom = ctx->chrom_ids[i];
    if (chrom == 0xFFFF || chrom >= ctx->n_refs || !ctx->refs[chrom]) return 0;
    const uint8_t *ref = ctx->refs[chrom];
    uint64_t ref_len = ctx->ref_lengths[chrom];
    int64_t pos = ctx->positions[i];
    if (pos < 1) return 0;
    const uint8_t *seq = ctx->sequences + ctx->seq_offsets[i];
    uint64_t seq_len = ctx->seq_offsets[i + 1] - ctx->seq_offsets[i];
    if (seq_len == 0 || (seq_len == 1 && seq[0] == '*')) return 0;
    const uint8_t *cig = ctx->cigars + ctx->cigar_offsets[i];
    uint64_t cig_len = ctx->cigar_offsets[i + 1] - ctx->cigar_offsets[i];
    if (cig_len == 0 || (cig_len == 1 && cig[0] == '*')) return 0;

    /* Pass 1: validate the CIGAR, the query length and the reference span. */
    uint64_t qlen = 0, rlen = 0, num = 0;
    int digits = 0;
    for (uint64_t k = 0; k < cig_len; k++) {
        uint8_t c = cig[k];
        if (c >= '0' && c <= '9') {
            if (++digits > 10) return 0;
            num = num * 10 + (uint64_t)(c - '0');
            continue;
        }
        if (digits == 0) return 0;
        switch (c) {
        case 'M': case '=': case 'X': qlen += num; rlen += num; break;
        case 'I': case 'S': qlen += num; break;
        case 'D': case 'N': rlen += num; break;
        case 'H': case 'P': break;
        default: return 0;
        }
        num = 0;
        digits = 0;
    }
    if (digits != 0) return 0;
    if (qlen != seq_len) return 0;
    if ((uint64_t)(pos - 1) > ref_len || rlen > ref_len - (uint64_t)(pos - 1)) return 0;

    /* Pass 2: build MD and NM. */
    stg_buf b = {0};
    uint64_t r = (uint64_t)(pos - 1), q = 0, run = 0;
    int64_t edits = 0;
    num = 0;
    for (uint64_t k = 0; k < cig_len; k++) {
        uint8_t c = cig[k];
        if (c >= '0' && c <= '9') {
            num = num * 10 + (uint64_t)(c - '0');
            continue;
        }
        switch (c) {
        case 'M': case '=': case 'X':
            for (uint64_t j = 0; j < num; j++, q++, r++) {
                uint8_t qb = seq[q], rb = up(ref[r]);
                if (qb == '=' || (up(qb) == rb && rb != 'N' && up(qb) != 'N')) {
                    run++;
                } else {
                    if (put_run(&b, run) || buf_u8(&b, rb)) goto oom;
                    run = 0;
                    edits++;
                }
            }
            break;
        case 'I':
            q += num;
            edits += (int64_t)num;
            break;
        case 'S':
            q += num;
            break;
        case 'D':
            if (put_run(&b, run) || buf_u8(&b, '^')) goto oom;
            for (uint64_t j = 0; j < num; j++, r++)
                if (buf_u8(&b, up(ref[r]))) goto oom;
            run = 0;
            edits += (int64_t)num;
            break;
        case 'N':
            r += num;
            break;
        default: /* H, P */
            break;
        }
        num = 0;
    }
    if (put_run(&b, run) || buf_u8(&b, 0)) goto oom;
    *md = (char *)b.p;
    *md_len = b.len - 1;
    *nm = edits;
    return 1;
oom:
    buf_free(&b);
    return -1;
}

/* ──────── encoder state ──────── */

/* Direct index of the (key, type, kind) column for kinds INT and TEXT:
 * key0 ∈ [A-Za-z] (52), key1 ∈ [A-Za-z0-9] (62), type ∈ AifZHB (6),
 * kind ∈ {INT, TEXT} (2). */
#define STG_SLOTS (52 * 62 * 6 * 2)

static int key0_idx(uint8_t c) {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return 26 + c - 'a';
    return -1;
}
static int key1_idx(uint8_t c) {
    if (c >= '0' && c <= '9') return 52 + c - '0';
    return key0_idx(c);
}
static int type_idx(uint8_t c) {
    switch (c) {
    case 'A': return 0; case 'i': return 1; case 'f': return 2;
    case 'Z': return 3; case 'H': return 4; case 'B': return 5;
    default: return -1;
    }
}
static int slot_of(const uint8_t key[2], uint8_t type, uint8_t kind) {
    int a = key0_idx(key[0]), b = key1_idx(key[1]), t = type_idx(type);
    if (a < 0 || b < 0 || t < 0 || kind > STG_KIND_TEXT) return -1;
    return ((a * 62 + b) * 6 + t) * 2 + kind;
}

typedef struct {
    uint8_t key[2], type, kind;
    uint64_t n;
    /* INT */
    int64_t *ints;
    size_t ints_cap;
    /* TEXT / VERBATIM: values joined by NUL, plus each value's start */
    stg_buf text;
    uint64_t *starts;
    size_t starts_cap;
} stg_col;

typedef struct {
    stg_col *cols;
    size_t n_cols, cap_cols;
    int32_t *slot; /* STG_SLOTS entries, -1 = none */
    int32_t verbatim_col;
    /* line dictionary: serialized entries, 5 bytes each */
    stg_buf line_bytes;
    uint64_t *line_off; /* start in line_bytes */
    uint32_t *line_n;   /* entries */
    size_t n_lines, cap_lines;
    uint32_t *hash;     /* open addressing, value = line index + 1 */
    size_t hash_cap;
} stg_enc;

static void enc_free(stg_enc *e) {
    for (size_t c = 0; c < e->n_cols; c++) {
        free(e->cols[c].ints);
        buf_free(&e->cols[c].text);
        free(e->cols[c].starts);
    }
    free(e->cols);
    free(e->slot);
    buf_free(&e->line_bytes);
    free(e->line_off);
    free(e->line_n);
    free(e->hash);
}

static int enc_col(stg_enc *e, const uint8_t key[2], uint8_t type, uint8_t kind) {
    int s = -1;
    if (kind == STG_KIND_VERBATIM) {
        if (e->verbatim_col >= 0) return e->verbatim_col;
    } else {
        s = slot_of(key, type, kind);
        if (s < 0) return -1;
        if (e->slot[s] >= 0) return e->slot[s];
    }
    if (e->n_cols == e->cap_cols) {
        size_t cap = e->cap_cols ? e->cap_cols * 2 : 16;
        stg_col *nc = (stg_col *)realloc(e->cols, cap * sizeof *nc);
        if (!nc) return -1;
        e->cols = nc;
        e->cap_cols = cap;
    }
    stg_col *c = &e->cols[e->n_cols];
    memset(c, 0, sizeof *c);
    c->key[0] = key[0];
    c->key[1] = key[1];
    c->type = type;
    c->kind = kind;
    int idx = (int)e->n_cols++;
    if (kind == STG_KIND_VERBATIM) e->verbatim_col = idx;
    else e->slot[s] = idx;
    return idx;
}

static int col_add_int(stg_col *c, int64_t v) {
    if (c->n == c->ints_cap) {
        size_t cap = c->ints_cap ? c->ints_cap * 2 : 1024;
        int64_t *p = (int64_t *)realloc(c->ints, cap * sizeof *p);
        if (!p) return -1;
        c->ints = p;
        c->ints_cap = cap;
    }
    c->ints[c->n++] = v;
    return 0;
}

static int col_add_text(stg_col *c, const uint8_t *s, size_t n) {
    if (c->n == c->starts_cap) {
        size_t cap = c->starts_cap ? c->starts_cap * 2 : 1024;
        uint64_t *p = (uint64_t *)realloc(c->starts, cap * sizeof *p);
        if (!p) return -1;
        c->starts = p;
        c->starts_cap = cap;
    }
    c->starts[c->n++] = c->text.len;
    if (buf_put(&c->text, s, n) || buf_u8(&c->text, 0)) return -1;
    return 0;
}

static uint32_t fnv1a(const uint8_t *p, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; i++) h = (h ^ p[i]) * 16777619u;
    return h;
}

static int hash_grow(stg_enc *e) {
    size_t cap = e->hash_cap ? e->hash_cap * 2 : 1024;
    uint32_t *h = (uint32_t *)calloc(cap, sizeof *h);
    if (!h) return -1;
    for (size_t l = 0; l < e->n_lines; l++) {
        uint32_t k = fnv1a(e->line_bytes.p + e->line_off[l], (size_t)e->line_n[l] * 5);
        size_t at = k & (cap - 1);
        while (h[at]) at = (at + 1) & (cap - 1);
        h[at] = (uint32_t)l + 1;
    }
    free(e->hash);
    e->hash = h;
    e->hash_cap = cap;
    return 0;
}

/* Returns the id of the line `ent[0 .. n*5)`, adding it if new; -1 on OOM. */
static int64_t enc_line(stg_enc *e, const uint8_t *ent, uint32_t n) {
    if ((e->n_lines + 1) * 2 > e->hash_cap && hash_grow(e)) return -1;
    size_t nb = (size_t)n * 5;
    uint32_t k = fnv1a(ent, nb);
    size_t at = k & (e->hash_cap - 1);
    while (e->hash[at]) {
        size_t l = e->hash[at] - 1;
        if (e->line_n[l] == n && memcmp(e->line_bytes.p + e->line_off[l], ent, nb) == 0)
            return (int64_t)l;
        at = (at + 1) & (e->hash_cap - 1);
    }
    if (e->n_lines == e->cap_lines) {
        size_t cap = e->cap_lines ? e->cap_lines * 2 : 256;
        uint64_t *o = (uint64_t *)realloc(e->line_off, cap * sizeof *o);
        if (!o) return -1;
        e->line_off = o;
        uint32_t *c = (uint32_t *)realloc(e->line_n, cap * sizeof *c);
        if (!c) return -1;
        e->line_n = c;
        e->cap_lines = cap;
    }
    e->line_off[e->n_lines] = e->line_bytes.len;
    e->line_n[e->n_lines] = n;
    if (buf_put(&e->line_bytes, ent, nb)) return -1;
    e->hash[at] = (uint32_t)e->n_lines + 1;
    return (int64_t)e->n_lines++;
}

/* ──────── column payloads (spec §2.2, §2.3) ──────── */

static int put_int_column(stg_buf *out, const stg_col *c) {
    int all_eq = 1;
    for (uint64_t k = 1; k < c->n && all_eq; k++) all_eq = c->ints[k] == c->ints[0];
    if (all_eq) return buf_u8(out, 2) || buf_uvarint(out, zz_enc(c->ints[0]));

    stg_buf raw = {0}, a = {0}, b = {0};
    int rc = -1;
    uint64_t mx = 0;
    for (uint64_t k = 0; k < c->n; k++) {
        uint64_t u = zz_enc(c->ints[k]);
        if (u > mx) mx = u;
        if (buf_uvarint(&raw, u)) goto done;
    }
    if (buf_u8(&a, 0) || put_substream(&a, raw.p, raw.len)) goto done;

    uint8_t w = mx < (1ull << 8) ? 1 : mx < (1ull << 16) ? 2 : mx < (1ull << 32) ? 4 : 8;
    if (buf_u8(&b, 1) || buf_u8(&b, w)) goto done;
    raw.len = 0;
    if (buf_reserve(&raw, (size_t)c->n)) goto done;
    for (uint8_t p = 0; p < w; p++) {
        for (uint64_t k = 0; k < c->n; k++)
            raw.p[k] = (uint8_t)(zz_enc(c->ints[k]) >> (8 * p));
        if (put_substream(&b, raw.p, (size_t)c->n)) goto done;
    }
    rc = (b.len < a.len) ? buf_put(out, b.p, b.len) : buf_put(out, a.p, a.len);
done:
    buf_free(&raw);
    buf_free(&a);
    buf_free(&b);
    return rc ? -1 : 0;
}

static int put_text_column(stg_buf *out, const stg_col *c) {
    const uint8_t *t = c->text.p;
    size_t len0 = (size_t)(c->n > 1 ? c->starts[1] : c->text.len) - 1;
    int all_eq = 1;
    for (uint64_t k = 1; k < c->n && all_eq; k++) {
        size_t lk = (size_t)((k + 1 < c->n ? c->starts[k + 1] : c->text.len) - c->starts[k]) - 1;
        all_eq = lk == len0 && memcmp(t + c->starts[k], t, len0) == 0;
    }
    if (all_eq) return buf_u8(out, 2) || buf_uvarint(out, len0) || buf_put(out, t, len0);

    stg_buf a = {0};
    int rc = -1;
    /* JOINED: the stored text minus its final NUL. */
    if (buf_u8(&a, 0) || put_substream(&a, t, c->text.len - 1)) goto done;

    int eligible = 1;
    for (uint64_t k = 0; k < c->n && eligible; k++)
        if (t[c->starts[k]] == 0) eligible = 0;
    for (size_t k = 0; k < c->text.len && eligible; k++)
        if (t[k] >= 0x80) eligible = 0;
    if (eligible) {
        const char **names = (const char **)malloc((size_t)c->n * sizeof *names);
        if (!names) goto done;
        for (uint64_t k = 0; k < c->n; k++) names[k] = (const char *)(t + c->starts[k]);
        size_t cap = ttio_name_tok_v2_max_encoded_size(c->n, c->text.len);
        uint8_t *nt = (uint8_t *)malloc(cap);
        size_t nl = cap;
        if (nt && ttio_name_tok_v2_encode((const char *const *)names, c->n, nt, &nl) == 0) {
            stg_buf b = {0};
            if (buf_u8(&b, 1) || buf_uvarint(&b, nl) || buf_put(&b, nt, nl)) {
                buf_free(&b);
                free(nt);
                free(names);
                goto done;
            }
            if (b.len < a.len) {
                buf_free(&a);
                a = b;
            } else {
                buf_free(&b);
            }
        }
        free(nt);
        free(names);
    }
    rc = buf_put(out, a.p, a.len);
done:
    buf_free(&a);
    return rc ? -1 : 0;
}

/* ──────── encode ──────── */

static int is_key0(uint8_t c) { return key0_idx(c) >= 0; }
static int is_key1(uint8_t c) { return key1_idx(c) >= 0; }

int ttio_sam_tags_encode(const ttio_sam_tags_ctx *ctx,
                         const uint8_t *tags, const uint64_t *tag_offsets,
                         uint8_t **out, size_t *out_len) {
    if (!ctx || !tag_offsets || !out || !out_len || (ctx->n_reads && !tags)) return TTIO_RANS_ERR_PARAM;
    if (ctx->n_refs && (!ctx->refs || !ctx->ref_lengths)) return TTIO_RANS_ERR_PARAM;
    const uint64_t n = ctx->n_reads;
    for (uint64_t i = 0; i < n; i++)
        if (tag_offsets[i + 1] < tag_offsets[i]) return TTIO_RANS_ERR_PARAM;

    stg_enc e;
    memset(&e, 0, sizeof e);
    e.verbatim_col = -1;
    stg_buf ids = {0}, ent = {0}, o = {0};
    int64_t *ivals = NULL;
    uint8_t *has_int = NULL;
    size_t ent_cap = 0;
    int rc = TTIO_RANS_ERR_ALLOC;

    e.slot = (int32_t *)malloc(STG_SLOTS * sizeof *e.slot);
    if (!e.slot) goto out;
    for (size_t s = 0; s < STG_SLOTS; s++) e.slot[s] = -1;
    if (hash_grow(&e)) goto out;

    for (uint64_t i = 0; i < n; i++) {
        const uint8_t *t = tags + tag_offsets[i];
        size_t tl = (size_t)(tag_offsets[i + 1] - tag_offsets[i]);
        if (memchr(t, 0, tl)) { rc = TTIO_RANS_ERR_PARAM; goto out; }

        /* Count fields and check that all are well-formed (§1.1). */
        size_t nf = 0;
        int verbatim = 0;
        if (tl > 0) {
            size_t a = 0;
            for (;;) {
                const uint8_t *tab = (const uint8_t *)memchr(t + a, '\t', tl - a);
                size_t b = tab ? (size_t)(tab - t) : tl;
                size_t fl = b - a;
                const uint8_t *f = t + a;
                if (fl < 5 || f[2] != ':' || f[4] != ':' || !is_key0(f[0]) || !is_key1(f[1]) ||
                    type_idx(f[3]) < 0)
                    verbatim = 1;
                nf++;
                if (!tab) break;
                a = b + 1;
            }
        }
        if (nf > ent_cap) {
            size_t cap = nf * 2;
            int64_t *iv = (int64_t *)realloc(ivals, cap * sizeof *iv);
            if (!iv) goto out;
            ivals = iv;
            uint8_t *hi = (uint8_t *)realloc(has_int, cap);
            if (!hi) goto out;
            has_int = hi;
            ent_cap = cap;
        }
        ent.len = 0;
        uint32_t n_ent = 0;
        if (verbatim) {
            static const uint8_t vk[2] = {0, 0};
            int c = enc_col(&e, vk, 0, STG_KIND_VERBATIM);
            if (c < 0 || col_add_text(&e.cols[c], t, tl)) goto out;
            uint8_t row[5] = {0, 0, 0, STG_KIND_VERBATIM, 0};
            if (buf_put(&ent, row, 5)) goto out;
            n_ent = 1;
        } else if (tl > 0) {
            char *md = NULL;
            size_t md_len = 0;
            int64_t nm = 0;
            int derivable = -2; /* not yet computed */
            size_t a = 0;
            for (;;) {
                const uint8_t *tab = (const uint8_t *)memchr(t + a, '\t', tl - a);
                size_t b = tab ? (size_t)(tab - t) : tl;
                const uint8_t *f = t + a, *v = f + 5;
                size_t vl = b - a - 5;
                uint8_t kind = STG_KIND_TEXT, arg = 0;
                int64_t iv = 0;
                int is_int = f[3] == 'i' && stg_parse_canonical_int(v, vl, &iv);
                if (is_int) kind = STG_KIND_INT;

                int md_tag = f[0] == 'M' && f[1] == 'D' && f[3] == 'Z';
                int nm_tag = f[0] == 'N' && f[1] == 'M' && f[3] == 'i';
                if (md_tag || nm_tag) {
                    if (derivable == -2) {
                        derivable = stg_calc_md_nm(ctx, i, &md, &md_len, &nm);
                        if (derivable < 0) goto out;
                    }
                    if (derivable == 1) {
                        if (md_tag && vl == md_len && memcmp(v, md, vl) == 0) kind = STG_KIND_DERIVED;
                        if (nm_tag && is_int && iv == nm) kind = STG_KIND_DERIVED;
                    }
                }
                if (kind == STG_KIND_INT) {
                    for (uint32_t j = 0; j < n_ent && j < 256; j++) {
                        if (has_int[j] && ivals[j] == iv) {
                            kind = STG_KIND_DUP;
                            arg = (uint8_t)j;
                            break;
                        }
                    }
                }
                has_int[n_ent] = kind == STG_KIND_INT || kind == STG_KIND_DUP ||
                                 (kind == STG_KIND_DERIVED && nm_tag);
                ivals[n_ent] = (kind == STG_KIND_DERIVED && nm_tag) ? nm : iv;

                if (kind == STG_KIND_INT || kind == STG_KIND_TEXT) {
                    int c = enc_col(&e, f, f[3], kind);
                    if (c < 0) { free(md); goto out; }
                    if (kind == STG_KIND_INT ? col_add_int(&e.cols[c], iv)
                                             : col_add_text(&e.cols[c], v, vl)) {
                        free(md);
                        goto out;
                    }
                }
                uint8_t row[5] = {f[0], f[1], f[3], kind, arg};
                if (buf_put(&ent, row, 5)) { free(md); goto out; }
                n_ent++;
                if (!tab) break;
                a = b + 1;
            }
            free(md);
        }
        int64_t lid = enc_line(&e, ent.p, n_ent);
        if (lid < 0 || buf_uvarint(&ids, (uint64_t)lid)) goto out;
    }

    /* Assemble (spec §3). */
    if (buf_put(&o, STG_MAGIC, 4) || buf_u8(&o, STG_VERSION) || buf_uvarint(&o, n) ||
        buf_uvarint(&o, e.n_lines))
        goto out;
    for (size_t l = 0; l < e.n_lines; l++)
        if (buf_uvarint(&o, e.line_n[l]) ||
            buf_put(&o, e.line_bytes.p + e.line_off[l], (size_t)e.line_n[l] * 5))
            goto out;
    if (put_substream(&o, ids.p, ids.len) || buf_uvarint(&o, e.n_cols)) goto out;
    for (size_t c = 0; c < e.n_cols; c++) {
        const stg_col *col = &e.cols[c];
        if (buf_put(&o, col->key, 2) || buf_u8(&o, col->type) || buf_u8(&o, col->kind) ||
            buf_uvarint(&o, col->n))
            goto out;
        if (col->kind == STG_KIND_INT ? put_int_column(&o, col) : put_text_column(&o, col)) goto out;
    }
    *out = o.p;
    *out_len = o.len;
    o.p = NULL;
    rc = 0;
out:
    enc_free(&e);
    buf_free(&ids);
    buf_free(&ent);
    buf_free(&o);
    free(ivals);
    free(has_int);
    return rc;
}

/* ──────── decode ──────── */

typedef struct {
    uint8_t key[2], type, kind;
    uint64_t n, next;
    int64_t *ints;              /* INT */
    uint8_t *text;              /* JOINED / CONST backing store */
    const uint8_t **vals;       /* TEXT: value starts */
    size_t *lens;
    char **names;               /* NAMETOK result, owned */
} stg_dcol;

static void dcol_free(stg_dcol *c) {
    free(c->ints);
    free(c->text);
    free(c->vals);
    free(c->lens);
    if (c->names) {
        for (uint64_t k = 0; k < c->n; k++) free(c->names[k]);
        free(c->names);
    }
}

static int read_int_column(stg_rd *r, stg_dcol *c) {
    uint8_t t;
    if (rd_u8(r, &t)) return TTIO_RANS_ERR_CORRUPT;
    c->ints = (int64_t *)malloc((size_t)(c->n ? c->n : 1) * sizeof *c->ints);
    if (!c->ints) return TTIO_RANS_ERR_ALLOC;
    if (t == 2) {
        uint64_t u;
        if (rd_uvarint(r, &u)) return TTIO_RANS_ERR_CORRUPT;
        for (uint64_t k = 0; k < c->n; k++) c->ints[k] = zz_dec(u);
        return 0;
    }
    if (t == 0) {
        uint8_t *raw;
        size_t rl;
        int rc = get_substream(r, &raw, &rl);
        if (rc) return rc;
        stg_rd s = {raw, rl, 0};
        for (uint64_t k = 0; k < c->n; k++) {
            uint64_t u;
            if (rd_uvarint(&s, &u)) { free(raw); return TTIO_RANS_ERR_CORRUPT; }
            c->ints[k] = zz_dec(u);
        }
        rc = s.pos == rl ? 0 : TTIO_RANS_ERR_CORRUPT;
        free(raw);
        return rc;
    }
    if (t == 1) {
        uint8_t w;
        if (rd_u8(r, &w) || (w != 1 && w != 2 && w != 4 && w != 8)) return TTIO_RANS_ERR_CORRUPT;
        uint64_t *u = (uint64_t *)calloc((size_t)(c->n ? c->n : 1), sizeof *u);
        if (!u) return TTIO_RANS_ERR_ALLOC;
        for (uint8_t p = 0; p < w; p++) {
            uint8_t *raw;
            size_t rl;
            int rc = get_substream(r, &raw, &rl);
            if (rc || rl != c->n) {
                if (!rc) free(raw);
                free(u);
                return rc ? rc : TTIO_RANS_ERR_CORRUPT;
            }
            for (uint64_t k = 0; k < c->n; k++) u[k] |= (uint64_t)raw[k] << (8 * p);
            free(raw);
        }
        for (uint64_t k = 0; k < c->n; k++) c->ints[k] = zz_dec(u[k]);
        free(u);
        return 0;
    }
    return TTIO_RANS_ERR_CORRUPT;
}

static int read_text_column(stg_rd *r, stg_dcol *c) {
    uint8_t t;
    if (rd_u8(r, &t)) return TTIO_RANS_ERR_CORRUPT;
    size_t n = (size_t)(c->n ? c->n : 1);
    c->vals = (const uint8_t **)malloc(n * sizeof *c->vals);
    c->lens = (size_t *)malloc(n * sizeof *c->lens);
    if (!c->vals || !c->lens) return TTIO_RANS_ERR_ALLOC;
    if (t == 2) {
        uint64_t len;
        const uint8_t *p;
        if (rd_uvarint(r, &len) || rd_bytes(r, (size_t)len, &p)) return TTIO_RANS_ERR_CORRUPT;
        for (uint64_t k = 0; k < c->n; k++) {
            c->vals[k] = p;
            c->lens[k] = (size_t)len;
        }
        return 0;
    }
    if (t == 0) {
        size_t rl;
        int rc = get_substream(r, &c->text, &rl);
        if (rc) return rc;
        uint64_t k = 0;
        size_t a = 0;
        for (size_t j = 0; j <= rl; j++) {
            if (j == rl || c->text[j] == 0) {
                if (k >= c->n) return TTIO_RANS_ERR_CORRUPT;
                c->vals[k] = c->text + a;
                c->lens[k] = j - a;
                k++;
                a = j + 1;
            }
        }
        return k == c->n ? 0 : TTIO_RANS_ERR_CORRUPT;
    }
    if (t == 1) {
        uint64_t len;
        const uint8_t *p;
        uint64_t got = 0;
        if (rd_uvarint(r, &len) || rd_bytes(r, (size_t)len, &p)) return TTIO_RANS_ERR_CORRUPT;
        if (ttio_name_tok_v2_decode(p, (size_t)len, &c->names, &got) != 0) {
            c->names = NULL;
            return TTIO_RANS_ERR_CORRUPT;
        }
        if (got != c->n) {
            for (uint64_t k = 0; k < got; k++) free(c->names[k]);
            free(c->names);
            c->names = NULL;
            return TTIO_RANS_ERR_CORRUPT;
        }
        for (uint64_t k = 0; k < c->n; k++) {
            c->vals[k] = (const uint8_t *)c->names[k];
            c->lens[k] = strlen(c->names[k]);
        }
        return 0;
    }
    return TTIO_RANS_ERR_CORRUPT;
}

int ttio_sam_tags_decode(const ttio_sam_tags_ctx *ctx,
                         const uint8_t *encoded, size_t encoded_len,
                         uint8_t **out_tags, uint64_t *out_tag_offsets) {
    if (!ctx || !encoded || !out_tags || !out_tag_offsets) return TTIO_RANS_ERR_PARAM;
    if (ctx->n_refs && (!ctx->refs || !ctx->ref_lengths)) return TTIO_RANS_ERR_PARAM;
    stg_rd r = {encoded, encoded_len, 0};
    const uint8_t *magic;
    uint8_t ver;
    uint64_t n, n_lines, n_cols;
    int rc = TTIO_RANS_ERR_CORRUPT;
    const uint8_t **line_ent = NULL;
    uint64_t *line_n = NULL;
    uint8_t *ids_raw = NULL;
    uint64_t *ids = NULL;
    stg_dcol *cols = NULL;
    int32_t *slot = NULL;
    int32_t verbatim_col = -1;
    stg_buf o = {0};
    uint64_t *spans = NULL; /* per entry: start, len in o */
    size_t spans_cap = 0;
    char *md = NULL;

    if (rd_bytes(&r, 4, &magic) || memcmp(magic, STG_MAGIC, 4) != 0 || rd_u8(&r, &ver) ||
        ver != STG_VERSION || rd_uvarint(&r, &n) || n != ctx->n_reads || rd_uvarint(&r, &n_lines) ||
        n_lines > encoded_len)
        goto out;

    line_ent = (const uint8_t **)malloc((size_t)(n_lines ? n_lines : 1) * sizeof *line_ent);
    line_n = (uint64_t *)malloc((size_t)(n_lines ? n_lines : 1) * sizeof *line_n);
    if (!line_ent || !line_n) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
    for (uint64_t l = 0; l < n_lines; l++) {
        if (rd_uvarint(&r, &line_n[l]) || line_n[l] > encoded_len ||
            rd_bytes(&r, (size_t)line_n[l] * 5, &line_ent[l]))
            goto out;
        /* Validate the line's entries (spec §3.1). */
        const uint8_t *e = line_ent[l];
        for (uint64_t j = 0; j < line_n[l]; j++) {
            const uint8_t *x = e + j * 5;
            uint8_t kind = x[3], arg = x[4];
            switch (kind) {
            case STG_KIND_INT:
                if (x[2] != 'i' || slot_of(x, x[2], kind) < 0 || arg) goto out;
                break;
            case STG_KIND_TEXT:
                if (slot_of(x, x[2], kind) < 0 || arg) goto out;
                break;
            case STG_KIND_DERIVED:
                if (!((x[0] == 'M' && x[1] == 'D' && x[2] == 'Z') ||
                      (x[0] == 'N' && x[1] == 'M' && x[2] == 'i')) || arg)
                    goto out;
                break;
            case STG_KIND_DUP: {
                if (x[2] != 'i' || slot_of(x, 'i', STG_KIND_INT) < 0 || arg >= j) goto out;
                const uint8_t *y = e + (size_t)arg * 5;
                int has = y[3] == STG_KIND_INT || y[3] == STG_KIND_DUP ||
                          (y[3] == STG_KIND_DERIVED && y[0] == 'N');
                if (!has) goto out;
                break;
            }
            case STG_KIND_VERBATIM:
                if (line_n[l] != 1 || x[0] || x[1] || x[2] || arg) goto out;
                break;
            default:
                goto out;
            }
        }
    }

    size_t ids_len;
    if ((rc = get_substream(&r, &ids_raw, &ids_len)) != 0) goto out;
    rc = TTIO_RANS_ERR_CORRUPT;
    ids = (uint64_t *)malloc((size_t)(n ? n : 1) * sizeof *ids);
    if (!ids) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
    uint64_t total_entries = 0; /* bounds every column's value count */
    {
        stg_rd s = {ids_raw, ids_len, 0};
        for (uint64_t i = 0; i < n; i++) {
            if (rd_uvarint(&s, &ids[i]) || ids[i] >= n_lines) goto out;
            total_entries += line_n[ids[i]];
        }
        if (s.pos != ids_len) goto out;
    }

    if (rd_uvarint(&r, &n_cols) || n_cols > encoded_len) goto out;
    cols = (stg_dcol *)calloc((size_t)(n_cols ? n_cols : 1), sizeof *cols);
    slot = (int32_t *)malloc(STG_SLOTS * sizeof *slot);
    if (!cols || !slot) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
    for (size_t s = 0; s < STG_SLOTS; s++) slot[s] = -1;
    for (uint64_t c = 0; c < n_cols; c++) {
        stg_dcol *dc = &cols[c];
        const uint8_t *k;
        if (rd_bytes(&r, 2, &k) || rd_u8(&r, &dc->type) || rd_u8(&r, &dc->kind) ||
            rd_uvarint(&r, &dc->n) || dc->n > total_entries)
            goto out;
        dc->key[0] = k[0];
        dc->key[1] = k[1];
        if (dc->kind == STG_KIND_VERBATIM) {
            if (k[0] || k[1] || dc->type || verbatim_col >= 0) goto out;
            verbatim_col = (int32_t)c;
        } else {
            int s = slot_of(dc->key, dc->type, dc->kind);
            if (s < 0 || slot[s] >= 0 || (dc->kind == STG_KIND_INT && dc->type != 'i')) goto out;
            slot[s] = (int32_t)c;
        }
        rc = dc->kind == STG_KIND_INT ? read_int_column(&r, dc) : read_text_column(&r, dc);
        if (rc) goto out;
        rc = TTIO_RANS_ERR_CORRUPT;
    }
    if (r.pos != r.len) goto out;

    /* Rebuild each read's tag text. */
    out_tag_offsets[0] = 0;
    for (uint64_t i = 0; i < n; i++) {
        const uint8_t *e = line_ent[ids[i]];
        uint64_t ne = line_n[ids[i]];
        if (ne * 2 > spans_cap) {
            size_t cap = (size_t)ne * 4;
            uint64_t *sp = (uint64_t *)realloc(spans, cap * sizeof *sp);
            if (!sp) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
            spans = sp;
            spans_cap = cap;
        }
        int derivable = -2;
        size_t md_len = 0;
        int64_t nm = 0;
        for (uint64_t j = 0; j < ne; j++) {
            const uint8_t *x = e + j * 5;
            uint8_t kind = x[3];
            if (kind == STG_KIND_VERBATIM) {
                stg_dcol *dc = &cols[verbatim_col < 0 ? 0 : verbatim_col];
                if (verbatim_col < 0 || dc->next >= dc->n) goto out;
                if (buf_put(&o, dc->vals[dc->next], dc->lens[dc->next])) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
                dc->next++;
                continue;
            }
            if (j > 0 && buf_u8(&o, '\t')) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
            uint8_t pre[5] = {x[0], x[1], ':', x[2], ':'};
            if (buf_put(&o, pre, 5)) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
            uint64_t start = o.len;
            int ok = 0;
            if (kind == STG_KIND_INT || kind == STG_KIND_TEXT) {
                int s = slot_of(x, x[2], kind);
                if (slot[s] < 0) goto out;
                stg_dcol *dc = &cols[slot[s]];
                if (dc->next >= dc->n) goto out;
                if (kind == STG_KIND_INT) {
                    char tmp[24];
                    int tn = snprintf(tmp, sizeof tmp, "%" PRId64, dc->ints[dc->next]);
                    ok = buf_put(&o, tmp, (size_t)tn) == 0;
                } else {
                    ok = buf_put(&o, dc->vals[dc->next], dc->lens[dc->next]) == 0;
                }
                dc->next++;
            } else if (kind == STG_KIND_DERIVED) {
                if (derivable == -2) {
                    derivable = stg_calc_md_nm(ctx, i, &md, &md_len, &nm);
                    if (derivable < 0) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
                }
                if (derivable != 1) goto out;
                if (x[0] == 'M') {
                    ok = buf_put(&o, md, md_len) == 0;
                } else {
                    char tmp[24];
                    int tn = snprintf(tmp, sizeof tmp, "%" PRId64, nm);
                    ok = buf_put(&o, tmp, (size_t)tn) == 0;
                }
            } else { /* DUP */
                uint64_t a = spans[x[4] * 2], l = spans[x[4] * 2 + 1];
                if (buf_reserve(&o, (size_t)l)) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
                memcpy(o.p + o.len, o.p + a, (size_t)l);
                o.len += (size_t)l;
                ok = 1;
            }
            if (!ok) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
            spans[j * 2] = start;
            spans[j * 2 + 1] = o.len - start;
        }
        free(md);
        md = NULL;
        out_tag_offsets[i + 1] = o.len;
    }
    for (uint64_t c = 0; c < n_cols; c++)
        if (cols[c].next != cols[c].n) goto out;

    if (!o.p && buf_reserve(&o, 1)) { rc = TTIO_RANS_ERR_ALLOC; goto out; }
    *out_tags = o.p;
    o.p = NULL;
    rc = 0;
out:
    free(md);
    free(line_ent);
    free(line_n);
    free(ids_raw);
    free(ids);
    if (cols)
        for (uint64_t c = 0; c < n_cols; c++) dcol_free(&cols[c]);
    free(cols);
    free(slot);
    free(spans);
    buf_free(&o);
    return rc;
}

void ttio_sam_tags_free(void *p) { free(p); }
