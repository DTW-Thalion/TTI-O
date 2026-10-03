/*
 * M103 Phase 0: how many bits per base can a reference-free context model
 * reach on read bases?
 *
 * Reads FASTQ on stdin, codes every read's A/C/G/T bases with an adaptive
 * model and reports the ideal code length (sum of -log2 p), which an
 * arithmetic coder reaches to within a fraction of a percent. Nothing is
 * written; this sizes the model before any codec is built.
 *
 * Model: each base is two binary decisions (high bit, then low bit given
 * the high bit). Each decision is predicted by one 16-bit probability
 * counter per (order, hashed context of the previous k bases, node), with
 * an adaptive learning rate. With several orders, the predictions are
 * combined by a logistic mixer (weights chosen by node and the most
 * confident order's count). Optionally, after a read is coded, the model
 * is also trained on the read's reverse complement, so that a read from the
 * other strand finds its context (the decoder can do the same).
 *
 * Bases other than ACGT (N, IUPAC) are not modelled: they reset the
 * context and are counted separately; the codec stores them as exceptions.
 *
 * --block-bases N resets the whole model (tables and mixer) at the first
 * read boundary after every N bases, as a codec that codes each blocks_v1
 * block on its own must (format-spec 10.12; blocks hold 64 MiB of bases).
 *
 * usage: seq_model_proof [--rc] [--table-bits B] [--block-bases N] k1 [k2 [k3]]
 *        < reads.fastq
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ORDERS 3

typedef struct {
    uint16_t p;  /* P(bit = 1) in 1/65536 */
    uint8_t n;   /* observations, saturating: sets the learning rate */
} counter;

typedef struct {
    int k;
    uint64_t mask;      /* 2k bits of history */
    int bits;           /* table index bits */
    counter *t;         /* [1 << bits][3] */
} order_model;

static int g_limit = 60;

   /* counter rate saturates at 1/(limit+1.5) */

static inline int squash(int d) { /* stretch domain (x/256) -> 12-bit prob */
    if (d > 2047) d = 2047;
    if (d < -2047) d = -2047;
    double x = d / 256.0;
    return (int)(4096.0 / (1.0 + exp(-x)));
}

static float stretch_tab[4096];
static inline float stretch(int p12) { return stretch_tab[p12 < 1 ? 1 : p12 > 4095 ? 4095 : p12]; }

static inline uint32_t ctx_index(const order_model *m, uint64_t hist) {
    uint64_t h = hist & m->mask;
    if (2 * m->k <= m->bits) return (uint32_t)h;
    h *= 0x9E3779B97F4A7C15ull;     /* Fibonacci hash */
    return (uint32_t)(h >> (64 - m->bits));
}

static inline void update(counter *c, int bit) {
    int target = bit ? 65535 : 0;
    int rate_div = c->n + 2;        /* 1/2, 1/3, ... then fixed */
    c->p += (target - (int)c->p) / rate_div;
    if (c->n < g_limit) c->n++;
}

typedef struct {
    int n_orders;
    order_model o[MAX_ORDERS];
    float w[3 * 64][MAX_ORDERS + 1];  /* mixer weights per (node, confidence bucket) */
    double bits;
} model;

static void model_reset(model *M) {
    for (int i = 0; i < M->n_orders; i++) {
        order_model *m = &M->o[i];
        for (size_t j = 0; j < ((size_t)3 << m->bits); j++) { m->t[j].p = 32768; m->t[j].n = 0; }
    }
    for (int i = 0; i < 3 * 64; i++)
        for (int j = 0; j <= MAX_ORDERS; j++) M->w[i][j] = (j < M->n_orders) ? 1.0f / M->n_orders : 0;
}

/* Code (or just train on) one base; returns nothing, accumulates bits. */
static void code_base(model *M, uint64_t hist, int base, int count_bits) {
    int node = 0;
    for (int b = 1; b >= 0; b--) {
        int bit = (base >> b) & 1;
        counter *c[MAX_ORDERS];
        float in[MAX_ORDERS + 1];
        int conf = 0;
        for (int i = 0; i < M->n_orders; i++) {
            order_model *m = &M->o[i];
            c[i] = &m->t[(size_t)ctx_index(m, hist) * 3 + node];
            in[i] = stretch(c[i]->p >> 4);
            if (c[i]->n > conf) conf = c[i]->n;
        }
        in[M->n_orders] = 0.25f;   /* bias */
        double p;
        float *w = M->w[node * 64 + (conf > 63 ? 63 : conf)];
        if (M->n_orders == 1) {
            p = (c[0]->p + 0.5) / 65536.0;
        } else {
            float dot = 0;
            for (int i = 0; i <= M->n_orders; i++) dot += w[i] * in[i];
            int p12 = squash((int)(dot * 256.0f));
            p = p12 / 4096.0;
            if (p < 1.0 / 4096) p = 1.0 / 4096;
            if (p > 4095.0 / 4096) p = 4095.0 / 4096;
            float err = (float)((bit ? 1.0 : 0.0) - p) * 0.02f;   /* learning rate */
            for (int i = 0; i <= M->n_orders; i++) w[i] += err * in[i];
        }
        if (count_bits) M->bits += -log2(bit ? p : 1.0 - p);
        for (int i = 0; i < M->n_orders; i++) update(c[i], bit);
        node = 1 + bit;
    }
}

static int code_of(char ch) {
    switch (ch) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

static void code_read(model *M, const char *s, size_t len, int count_bits,
                      uint64_t *n_bases, uint64_t *n_other) {
    uint64_t hist = 0;
    for (size_t i = 0; i < len; i++) {
        int c = code_of(s[i]);
        if (c < 0) { hist = 0; if (count_bits) (*n_other)++; continue; }
        code_base(M, hist, c, count_bits);
        hist = (hist << 2) | (uint64_t)c;
        if (count_bits) (*n_bases)++;
    }
}

int main(int argc, char **argv) {
    int rc = 0, table_bits = 24, ai = 1;
    uint64_t block_bases = 0, since_reset = 0, n_resets = 0;
    while (ai < argc && argv[ai][0] == '-') {
        if (!strcmp(argv[ai], "--rc")) rc = 1;
        else if (!strcmp(argv[ai], "--table-bits") && ai + 1 < argc) table_bits = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--block-bases") && ai + 1 < argc) block_bases = strtoull(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "--limit") && ai + 1 < argc) g_limit = atoi(argv[++ai]);
        ai++;
    }
    model M;
    memset(&M, 0, sizeof M);
    for (; ai < argc && M.n_orders < MAX_ORDERS; ai++) {
        order_model *m = &M.o[M.n_orders++];
        m->k = atoi(argv[ai]);
        m->mask = (m->k >= 32) ? ~0ull : ((1ull << (2 * m->k)) - 1);
        m->bits = (2 * m->k < table_bits) ? 2 * m->k : table_bits;
        m->t = calloc((size_t)3 << m->bits, sizeof(counter));
        if (!m->t) { fprintf(stderr, "out of memory\n"); return 1; }
        for (size_t j = 0; j < ((size_t)3 << m->bits); j++) m->t[j].p = 32768;
    }
    if (!M.n_orders) { fprintf(stderr, "usage: %s [--rc] [--table-bits B] k1 [k2 [k3]]\n", argv[0]); return 2; }
    for (int i = 0; i < 4096; i++) {
        double p = (i + 0.5) / 4096.0;
        stretch_tab[i] = (float)log(p / (1 - p));
    }
    for (int i = 0; i < 3 * 64; i++)
        for (int j = 0; j <= MAX_ORDERS; j++) M.w[i][j] = (j < M.n_orders) ? 1.0f / M.n_orders : 0;

    char *line = NULL, *rcbuf = NULL;
    size_t cap = 0, rccap = 0;
    ssize_t len;
    uint64_t n_bases = 0, n_other = 0, n_reads = 0, ln = 0;
    while ((len = getline(&line, &cap, stdin)) > 0) {
        if ((ln++ & 3) != 1) continue;          /* FASTQ: sequence is line 2 of 4 */
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) len--;
        if (block_bases && since_reset >= block_bases) { model_reset(&M); since_reset = 0; n_resets++; }
        since_reset += (uint64_t)len;
        code_read(&M, line, (size_t)len, 1, &n_bases, &n_other);
        if (rc) {
            if ((size_t)len > rccap) { rccap = (size_t)len * 2; rcbuf = realloc(rcbuf, rccap); }
            for (ssize_t i = 0; i < len; i++) {
                char ch = line[len - 1 - i];
                rcbuf[i] = ch == 'A' ? 'T' : ch == 'C' ? 'G' : ch == 'G' ? 'C' : ch == 'T' ? 'A' : 'N';
            }
            code_read(&M, rcbuf, (size_t)len, 0, NULL, NULL);
        }
        n_reads++;
    }
    printf("orders");
    for (int i = 0; i < M.n_orders; i++) printf(" %d", M.o[i].k);
    if (block_bases) printf(" block_bases %llu (%llu resets)", (unsigned long long)block_bases, (unsigned long long)n_resets);
    printf("%s table_bits %d: reads %llu bases %llu other %llu  bits/base %.4f  bytes %.0f\n",
           rc ? " +rc" : "", table_bits, (unsigned long long)n_reads, (unsigned long long)n_bases,
           (unsigned long long)n_other, n_bases ? M.bits / n_bases : 0.0, M.bits / 8);
    return 0;
}
