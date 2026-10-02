/*
 * seq_cm.c — SEQ_CM (codec id 19, M103): read bases coded with a
 * reference-free context-mixing model and a binary arithmetic coder.
 *
 * Each A/C/G/T base is two binary decisions (high bit, then low bit given
 * the high bit). Every decision is predicted by one adaptive counter per
 * model order, indexed by the previous k bases of the read, and the
 * predictions are combined by a logistic mixer. With the RC flag the model
 * also trains on each read's reverse complement once the read is coded, so
 * a read from the other strand finds its context; the decoder does the
 * same. Bytes other than A/C/G/T are stored as exception runs and reset
 * the context. Everything is integer arithmetic, so the stream is the same
 * on every platform. Spec: docs/codecs/seq_cm.md.
 *
 * Copyright (c) 2026 Thalion Global. All rights reserved.
 */
#include "seq_cm.h"

#include <stdlib.h>
#include <string.h>

#define SQC_MAX_ORDERS 3
#define SQC_MIX_SETS (3 * 64)       /* node x confidence bucket */
#define SQC_BIAS_INPUT 64           /* stretch units: 0.25 */

/* ── Integer logistic functions (spec §3.2) ──────────────────────────── */

static const int SQUASH_T[33] = {
    1, 2, 3, 6, 10, 16, 27, 45, 73, 120, 194, 310, 488, 747, 1101,
    1546, 2047, 2549, 2994, 3348, 3607, 3785, 3901, 3975, 4022,
    4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094
};

int sqc_squash(int d)
{
    if (d > 2047) d = 2047;
    if (d < -2047) d = -2047;
    int x = d + 2048;               /* 1 .. 4095, never negative */
    int i = x >> 7, w = x & 127;
    return (SQUASH_T[i] * (128 - w) + SQUASH_T[i + 1] * w + 64) >> 7;
}

static int16_t STRETCH_T[4096];
static int stretch_ready;

static void stretch_init(void)
{
    /* Inverse of squash: the smallest d whose squash reaches p. Built
     * once; every caller builds the same table, so a race only writes the
     * same values twice. */
    int pi = 0;
    for (int d = -2047; d <= 2047; d++) {
        int v = sqc_squash(d);
        for (int p = pi; p <= v; p++) STRETCH_T[p] = (int16_t)d;
        pi = v + 1;
    }
    for (int p = pi; p < 4096; p++) STRETCH_T[p] = 2047;
    stretch_ready = 1;
}

int sqc_stretch(int p12)
{
    if (!stretch_ready) stretch_init();
    if (p12 < 0) p12 = 0;
    if (p12 > 4095) p12 = 4095;
    return STRETCH_T[p12];
}

/* Floor division by 2^s for signed values, without relying on the
 * implementation-defined right shift of a negative number. */
static inline int64_t asr64(int64_t x, int s)
{
    return x < 0 ? ~((~x) >> s) : x >> s;
}

uint8_t sqc_auto_table_bits(uint64_t n_bases)
{
    int b = 0;
    while (b < 63 && ((uint64_t)1 << b) < n_bases) b++;
    b += 1;
    if (b < 16) b = 16;
    if (b > 24) b = 24;
    return (uint8_t)b;
}

/* ── Byte buffer ─────────────────────────────────────────────────────── */

typedef struct {
    uint8_t *p;
    size_t len, cap;
    int err;
} buf_t;

static void buf_put(buf_t *b, uint8_t v)
{
    if (b->len == b->cap) {
        size_t nc = b->cap ? b->cap * 2 : 4096;
        uint8_t *np = (uint8_t *)realloc(b->p, nc);
        if (!np) { b->err = 1; return; }
        b->p = np;
        b->cap = nc;
    }
    b->p[b->len++] = v;
}

static void buf_varint(buf_t *b, uint64_t v)
{
    while (v >= 0x80) { buf_put(b, (uint8_t)(v | 0x80)); v >>= 7; }
    buf_put(b, (uint8_t)v);
}

static int read_varint(const uint8_t *p, size_t n, size_t *pos, uint64_t *out)
{
    uint64_t v = 0;
    for (int shift = 0; shift < 64; shift += 7) {
        if (*pos >= n) return -1;
        uint8_t c = p[(*pos)++];
        v |= (uint64_t)(c & 0x7F) << shift;
        if (!(c & 0x80)) { *out = v; return 0; }
    }
    return -1;
}

static void put_u16(uint8_t *p, uint16_t v) { p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8); }
static void put_u64(uint8_t *p, uint64_t v) { for (int i = 0; i < 8; i++) p[i] = (uint8_t)(v >> (8 * i)); }
static uint16_t get_u16(const uint8_t *p) { return (uint16_t)(p[0] | (p[1] << 8)); }
static uint64_t get_u64(const uint8_t *p)
{
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8 * i);
    return v;
}

/* ── Binary arithmetic coder (spec §4) ───────────────────────────────── */

typedef struct {
    uint32_t x1, x2;
    buf_t *out;
} bac_enc;

typedef struct {
    uint32_t x1, x2, x;
    const uint8_t *in;
    size_t pos, len;
} bac_dec;

/* p16 = P(bit = 1) in 1/65536, 1 .. 65535. */
static inline void bac_encode(bac_enc *e, int bit, uint32_t p16)
{
    uint32_t xmid = e->x1 + (uint32_t)(((uint64_t)(e->x2 - e->x1) * p16) >> 16);
    if (bit) e->x2 = xmid; else e->x1 = xmid + 1;
    while (((e->x1 ^ e->x2) & 0xFF000000u) == 0) {
        buf_put(e->out, (uint8_t)(e->x2 >> 24));
        e->x1 <<= 8;
        e->x2 = (e->x2 << 8) | 0xFF;
    }
}

static void bac_flush(bac_enc *e)
{
    for (int i = 0; i < 4; i++) { buf_put(e->out, (uint8_t)(e->x1 >> 24)); e->x1 <<= 8; }
}

static inline uint8_t bac_byte(bac_dec *d)
{
    return d->pos < d->len ? d->in[d->pos++] : 0;
}

static void bac_dec_init(bac_dec *d, const uint8_t *in, size_t len)
{
    d->x1 = 0;
    d->x2 = 0xFFFFFFFFu;
    d->in = in;
    d->len = len;
    d->pos = 0;
    d->x = 0;
    for (int i = 0; i < 4; i++) d->x = (d->x << 8) | bac_byte(d);
}

static inline int bac_decode(bac_dec *d, uint32_t p16)
{
    uint32_t xmid = d->x1 + (uint32_t)(((uint64_t)(d->x2 - d->x1) * p16) >> 16);
    int bit = d->x <= xmid;
    if (bit) d->x2 = xmid; else d->x1 = xmid + 1;
    while (((d->x1 ^ d->x2) & 0xFF000000u) == 0) {
        d->x1 <<= 8;
        d->x2 = (d->x2 << 8) | 0xFF;
        d->x = (d->x << 8) | bac_byte(d);
    }
    return bit;
}

/* ── Model (spec §3) ─────────────────────────────────────────────────── */

/* A counter is one uint32: probability P(bit = 1) in the high 16 bits,
 * stored XOR 0x8000 so that a zeroed table means p = 1/2, and the
 * observation count in the low 16 bits. */
typedef struct {
    int k;
    int bits;
    int direct;
    uint64_t mask;
    uint32_t *t;            /* [1 << bits][3] */
} order_t;

typedef struct {
    int n_orders;
    order_t o[SQC_MAX_ORDERS];
    int32_t w[SQC_MIX_SETS][SQC_MAX_ORDERS + 1];
    uint32_t limit;
    uint32_t lr;
    uint32_t *rate_tab;     /* [limit + 1]: 65536 / (n + 2) */
    /* per-decision scratch */
    uint32_t *slot[SQC_MAX_ORDERS];
    int st[SQC_MAX_ORDERS + 1];
    int set;
    int p12;
} model_t;

static void model_free(model_t *m)
{
    for (int i = 0; i < m->n_orders; i++) free(m->o[i].t);
    free(m->rate_tab);
}

static int model_init(model_t *m, int n_orders, const uint8_t *orders,
                      int table_bits, uint32_t limit, uint32_t lr)
{
    memset(m, 0, sizeof *m);
    if (!stretch_ready) stretch_init();
    m->n_orders = n_orders;
    m->limit = limit;
    m->lr = lr;
    m->rate_tab = (uint32_t *)malloc(((size_t)limit + 1) * sizeof(uint32_t));
    if (!m->rate_tab) return TTIO_RANS_ERR_ALLOC;
    for (uint32_t n = 0; n <= limit; n++) m->rate_tab[n] = 65536u / (n + 2);
    for (int i = 0; i < n_orders; i++) {
        order_t *o = &m->o[i];
        o->k = orders[i];
        o->mask = o->k >= 32 ? ~(uint64_t)0 : (((uint64_t)1 << (2 * o->k)) - 1);
        o->direct = 2 * o->k <= table_bits;
        o->bits = o->direct ? 2 * o->k : table_bits;
        o->t = (uint32_t *)calloc((size_t)3 << o->bits, sizeof(uint32_t));
        if (!o->t) { m->n_orders = i; model_free(m); return TTIO_RANS_ERR_ALLOC; }
    }
    m->n_orders = n_orders;
    for (int s = 0; s < SQC_MIX_SETS; s++)
        for (int j = 0; j <= SQC_MAX_ORDERS; j++)
            m->w[s][j] = j < n_orders ? (int32_t)(65536 / n_orders) : 0;
    return 0;
}

static inline uint64_t ctx_row(const order_t *o, uint64_t hist)
{
    uint64_t h = hist & o->mask;
    if (o->direct) return h;
    return (h * 0x9E3779B97F4A7C15ull) >> (64 - o->bits);
}

/* Predict one decision: fills the scratch and returns p16 for the coder. */
static inline uint32_t model_predict(model_t *m, const uint64_t *rows, int node)
{
    uint32_t conf = 0;
    int64_t dot = 0;
    for (int i = 0; i < m->n_orders; i++) {
        uint32_t *c = &m->o[i].t[rows[i] * 3 + (uint64_t)node];
        m->slot[i] = c;
        uint32_t p = (*c >> 16) ^ 0x8000u;
        uint32_t n = *c & 0xFFFFu;
        if (n > conf) conf = n;
        m->st[i] = sqc_stretch((int)(p >> 4));
    }
    m->st[m->n_orders] = SQC_BIAS_INPUT;
    m->set = node * 64 + (int)(conf > 63 ? 63 : conf);
    const int32_t *w = m->w[m->set];
    for (int i = 0; i <= m->n_orders; i++) dot += (int64_t)w[i] * m->st[i];
    int p12 = sqc_squash((int)asr64(dot, 16));
    if (p12 < 1) p12 = 1;
    if (p12 > 4095) p12 = 4095;
    m->p12 = p12;
    return (uint32_t)p12 << 4;
}

static inline void model_update(model_t *m, int bit)
{
    int64_t err = (int64_t)((bit << 12) - m->p12) * m->lr;
    int32_t *w = m->w[m->set];
    for (int i = 0; i <= m->n_orders; i++)
        w[i] += (int32_t)asr64(err * m->st[i], 16);
    for (int i = 0; i < m->n_orders; i++) {
        uint32_t c = *m->slot[i];
        uint32_t p = (c >> 16) ^ 0x8000u;
        uint32_t n = c & 0xFFFFu;
        uint32_t r = m->rate_tab[n];
        if (bit) p += ((65535u - p) * r) >> 16;
        else     p -= (p * r) >> 16;
        if (n < m->limit) n++;
        *m->slot[i] = ((p ^ 0x8000u) << 16) | n;
    }
}

static inline void rows_of(const model_t *m, uint64_t hist, uint64_t *rows)
{
    for (int i = 0; i < m->n_orders; i++) rows[i] = ctx_row(&m->o[i], hist);
}

/* Train on one base without coding it. */
static inline void model_train_base(model_t *m, uint64_t hist, int base)
{
    uint64_t rows[SQC_MAX_ORDERS];
    rows_of(m, hist, rows);
    int hi = (base >> 1) & 1, lo = base & 1;
    model_predict(m, rows, 0);
    model_update(m, hi);
    model_predict(m, rows, 1 + hi);
    model_update(m, lo);
}

static inline int base_code(uint8_t c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:  return -1;
    }
}

static const uint8_t BASE_CHAR[4] = { 'A', 'C', 'G', 'T' };

/* Reverse-complement training on read bytes s[0..len) (spec §3.4). */
static void model_train_rc(model_t *m, const uint8_t *s, uint64_t len)
{
    uint64_t hist = 0;
    for (uint64_t i = len; i-- > 0;) {
        int c = base_code(s[i]);
        if (c < 0) { hist = 0; continue; }
        c = 3 - c;
        model_train_base(m, hist, c);
        hist = (hist << 2) | (uint64_t)c;
    }
}

/* ── Parameters ──────────────────────────────────────────────────────── */

static int check_params(int n_orders, const uint8_t *orders, int table_bits,
                        uint32_t limit, uint32_t lr)
{
    if (n_orders < 1 || n_orders > SQC_MAX_ORDERS) return -1;
    for (int i = 0; i < n_orders; i++) {
        if (orders[i] < 1 || orders[i] > 31) return -1;
        if (i && orders[i] <= orders[i - 1]) return -1;
    }
    if (table_bits < SQC_MIN_TABLE_BITS || table_bits > SQC_MAX_TABLE_BITS) return -1;
    if (limit < 1 || limit > 1023) return -1;
    if (lr < 1 || lr > 4096) return -1;
    return 0;
}

void ttio_seq_cm_default_params(ttio_seq_cm_params *p)
{
    memset(p, 0, sizeof *p);
    p->n_orders = 3;
    p->orders[0] = 11;
    p->orders[1] = 16;
    p->orders[2] = 24;
    p->table_bits = 0;          /* chosen from the base count */
    p->flags = TTIO_SEQ_CM_FLAG_RC;
    p->limit = 60;
    p->lr = 82;
}

/* ── Encode ──────────────────────────────────────────────────────────── */

int ttio_seq_cm_encode(const uint8_t *seq, const uint64_t *lengths,
                       uint64_t n_reads, const ttio_seq_cm_params *params,
                       uint8_t **out, size_t *out_len)
{
    if (!out || !out_len || (n_reads && !lengths)) return TTIO_RANS_ERR_PARAM;
    *out = NULL;
    *out_len = 0;
    ttio_seq_cm_params p;
    if (params) p = *params; else ttio_seq_cm_default_params(&p);
    if (p.flags & ~(uint8_t)TTIO_SEQ_CM_FLAG_RC) return TTIO_RANS_ERR_PARAM;

    uint64_t n_bases = 0;
    for (uint64_t i = 0; i < n_reads; i++) {
        if (lengths[i] > UINT64_MAX - n_bases) return TTIO_RANS_ERR_PARAM;
        n_bases += lengths[i];
    }
    if (n_bases && !seq) return TTIO_RANS_ERR_PARAM;
    if (p.table_bits == 0) p.table_bits = sqc_auto_table_bits(n_bases);
    if (check_params(p.n_orders, p.orders, p.table_bits, p.limit, p.lr))
        return TTIO_RANS_ERR_PARAM;

    /* Exception runs: (gap from the previous run's end, run length, byte). */
    buf_t exc = {0};
    uint64_t n_runs = 0, prev_end = 0;
    for (uint64_t i = 0; i < n_bases;) {
        if (base_code(seq[i]) >= 0) { i++; continue; }
        uint64_t j = i + 1;
        while (j < n_bases && seq[j] == seq[i]) j++;
        buf_varint(&exc, i - prev_end);
        buf_varint(&exc, j - i);
        buf_put(&exc, seq[i]);
        n_runs++;
        prev_end = j;
        i = j;
    }
    if (exc.err) { free(exc.p); return TTIO_RANS_ERR_ALLOC; }

    buf_t ob = {0};
    uint64_t exc_len = exc.len;
    for (int i = 0; i < SQC_HEADER_LEN; i++) buf_put(&ob, 0);
    for (size_t i = 0; i < exc.len; i++) buf_put(&ob, exc.p[i]);
    free(exc.p);
    if (ob.err) { free(ob.p); return TTIO_RANS_ERR_ALLOC; }

    model_t m;
    int rc = model_init(&m, p.n_orders, p.orders, p.table_bits, p.limit, p.lr);
    if (rc) { free(ob.p); return rc; }

    bac_enc e = { 0, 0xFFFFFFFFu, &ob };
    uint64_t off = 0;
    for (uint64_t r = 0; r < n_reads; r++) {
        const uint8_t *s = seq + off;
        uint64_t hist = 0, rows[SQC_MAX_ORDERS];
        for (uint64_t i = 0; i < lengths[r]; i++) {
            int c = base_code(s[i]);
            if (c < 0) { hist = 0; continue; }
            int hi = (c >> 1) & 1, lo = c & 1;
            rows_of(&m, hist, rows);
            bac_encode(&e, hi, model_predict(&m, rows, 0));
            model_update(&m, hi);
            bac_encode(&e, lo, model_predict(&m, rows, 1 + hi));
            model_update(&m, lo);
            hist = (hist << 2) | (uint64_t)c;
        }
        if (p.flags & TTIO_SEQ_CM_FLAG_RC) model_train_rc(&m, s, lengths[r]);
        off += lengths[r];
    }
    bac_flush(&e);
    model_free(&m);
    if (ob.err) { free(ob.p); return TTIO_RANS_ERR_ALLOC; }

    uint8_t *h = ob.p;
    memcpy(h, SQC_MAGIC, 4);
    h[4] = SQC_VERSION;
    h[5] = p.flags;
    h[6] = p.n_orders;
    for (int i = 0; i < 3; i++) h[7 + i] = i < p.n_orders ? p.orders[i] : 0;
    h[10] = p.table_bits;
    h[11] = 0;
    put_u16(h + 12, p.limit);
    put_u16(h + 14, p.lr);
    put_u64(h + 16, n_reads);
    put_u64(h + 24, n_bases);
    put_u64(h + 32, n_runs);
    put_u64(h + 40, exc_len);
    *out = ob.p;
    *out_len = ob.len;
    return 0;
}

/* ── Decode ──────────────────────────────────────────────────────────── */

int ttio_seq_cm_decode(const uint8_t *in, size_t in_len,
                       const uint64_t *lengths, uint64_t n_reads,
                       uint8_t **out_seq, size_t *out_len)
{
    if (!in || !out_seq || !out_len || (n_reads && !lengths)) return TTIO_RANS_ERR_PARAM;
    *out_seq = NULL;
    *out_len = 0;
    if (in_len < SQC_HEADER_LEN || memcmp(in, SQC_MAGIC, 4) || in[4] != SQC_VERSION)
        return TTIO_RANS_ERR_CORRUPT;
    uint8_t flags = in[5];
    int n_orders = in[6];
    uint8_t orders[3] = { in[7], in[8], in[9] };
    int table_bits = in[10];
    uint32_t limit = get_u16(in + 12), lr = get_u16(in + 14);
    uint64_t h_reads = get_u64(in + 16), n_bases = get_u64(in + 24);
    uint64_t n_runs = get_u64(in + 32), exc_len = get_u64(in + 40);
    if (flags & ~(uint8_t)TTIO_SEQ_CM_FLAG_RC || in[11]) return TTIO_RANS_ERR_CORRUPT;
    if (check_params(n_orders, orders, table_bits, limit, lr)) return TTIO_RANS_ERR_CORRUPT;
    for (int i = n_orders; i < 3; i++) if (orders[i]) return TTIO_RANS_ERR_CORRUPT;
    if (h_reads != n_reads) return TTIO_RANS_ERR_PARAM;
    uint64_t total = 0;
    for (uint64_t i = 0; i < n_reads; i++) {
        if (lengths[i] > UINT64_MAX - total) return TTIO_RANS_ERR_PARAM;
        total += lengths[i];
    }
    if (total != n_bases) return TTIO_RANS_ERR_PARAM;
    if (exc_len > in_len - SQC_HEADER_LEN) return TTIO_RANS_ERR_CORRUPT;
    if (n_bases > (uint64_t)SIZE_MAX - 1) return TTIO_RANS_ERR_ALLOC;

    uint8_t *s = (uint8_t *)malloc(n_bases ? (size_t)n_bases : 1);
    if (!s) return TTIO_RANS_ERR_ALLOC;
    /* Exception runs are read one at a time as the walk reaches them. */
    const uint8_t *ex = in + SQC_HEADER_LEN;
    size_t epos = 0;
    uint64_t run_start = UINT64_MAX, run_end = 0, prev_end = 0, runs_left = n_runs;
    uint8_t run_byte = 0;

#define NEXT_RUN()                                                          \
    do {                                                                    \
        if (runs_left) {                                                    \
            uint64_t gap, rl;                                               \
            if (read_varint(ex, (size_t)exc_len, &epos, &gap)               \
                || read_varint(ex, (size_t)exc_len, &epos, &rl)             \
                || epos >= exc_len || rl == 0                               \
                || gap > n_bases - prev_end || rl > n_bases - prev_end - gap) \
                { free(s); return TTIO_RANS_ERR_CORRUPT; }                  \
            run_byte = ex[epos++];                                          \
            if (base_code(run_byte) >= 0) { free(s); return TTIO_RANS_ERR_CORRUPT; } \
            run_start = prev_end + gap;                                     \
            run_end = run_start + rl;                                       \
            prev_end = run_end;                                             \
            runs_left--;                                                    \
        } else {                                                            \
            run_start = UINT64_MAX;                                         \
            run_end = UINT64_MAX;                                           \
        }                                                                   \
    } while (0)

    NEXT_RUN();

    model_t m;
    int rc = model_init(&m, n_orders, orders, table_bits, limit, lr);
    if (rc) { free(s); return rc; }
    bac_dec d;
    bac_dec_init(&d, in + SQC_HEADER_LEN + exc_len, (size_t)(in_len - SQC_HEADER_LEN - exc_len));

    uint64_t off = 0;
    for (uint64_t r = 0; r < n_reads; r++) {
        uint64_t hist = 0, rows[SQC_MAX_ORDERS];
        for (uint64_t i = 0; i < lengths[r]; i++) {
            uint64_t g = off + i;
            if (g >= run_start) {
                s[g] = run_byte;
                hist = 0;
                if (g + 1 == run_end) NEXT_RUN();
                continue;
            }
            rows_of(&m, hist, rows);
            int hi = bac_decode(&d, model_predict(&m, rows, 0));
            model_update(&m, hi);
            int lo = bac_decode(&d, model_predict(&m, rows, 1 + hi));
            model_update(&m, lo);
            int c = (hi << 1) | lo;
            s[g] = BASE_CHAR[c];
            hist = (hist << 2) | (uint64_t)c;
        }
        if (flags & TTIO_SEQ_CM_FLAG_RC) model_train_rc(&m, s + off, lengths[r]);
        off += lengths[r];
    }
#undef NEXT_RUN
    model_free(&m);
    if (runs_left || epos != exc_len) { free(s); return TTIO_RANS_ERR_CORRUPT; }
    *out_seq = s;
    *out_len = (size_t)n_bases;
    return 0;
}

void ttio_seq_cm_free(void *p) { free(p); }
