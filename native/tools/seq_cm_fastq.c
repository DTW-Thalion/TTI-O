/*
 * seq_cm_fastq — code a FASTQ file's bases with SEQ_CM (codec id 19) the
 * way the blocks_v1 writer does: the reads are cut into blocks of about
 * --block-bases bases (default 64 MiB, the writer's block size), each block
 * coded on its own, then decoded and compared.
 *
 * --group reorders the whole run before blocking: reads are sorted by their
 * canonical minimizer (the k-mer of least hash over both strands, k =
 * --group-k) and then by the read's offset from it, so reads that overlap
 * tend to land in the same block. The permutation that restores the input
 * order is reported as its ideal size (log2 n!) and as fixed-width indices
 * (n * ceil(log2 n) bits). --pairs keeps reads 2i and 2i+1 together and
 * permutes pairs, keyed by the mate with the smaller minimizer.
 *
 * --chain lays reads out by overlap instead (see chain_order). --as-is
 * keeps the input order through the same path. --pos takes one mapping
 * position per read (input order, 0 for unmapped) and reports how local the
 * resulting order is; --no-code skips coding for that diagnostic.
 *
 * usage: seq_cm_fastq [--block-bases N] [--table-bits B] [--no-rc]
 *                     [--group [--group-k K] [--pairs] | --chain [--group-k K]
 *                      [--group-w W] [--max-occ N] [--min-votes V] | --as-is]
 *                     [--pos POSITIONS] [--no-code]
 *                     [k1 [k2 [k3]]] < reads.fastq
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ttio_rans.h"

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

typedef struct {
    ttio_seq_cm_params p;
    unsigned long long total_bases, total_bytes, blocks;
    double t_enc, t_dec;
} coder;

static int code_block(coder *c, const uint8_t *seq, const uint64_t *len, size_t n, size_t used)
{
    uint8_t *enc = NULL, *dec = NULL;
    size_t enc_len = 0, dec_len = 0;
    double t0 = now();
    int rc = ttio_seq_cm_encode(seq, len, n, &c->p, &enc, &enc_len);
    double t1 = now();
    if (!rc) rc = ttio_seq_cm_decode(enc, enc_len, len, n, &dec, &dec_len);
    double t2 = now();
    if (rc || dec_len != used || memcmp(dec, seq, used)) {
        fprintf(stderr, "block %llu: round trip failed (rc=%d)\n", c->blocks, rc);
        return 1;
    }
    c->t_enc += t1 - t0;
    c->t_dec += t2 - t1;
    c->total_bases += used;
    c->total_bytes += enc_len;
    c->blocks++;
    ttio_seq_cm_free(enc);
    ttio_seq_cm_free(dec);
    return 0;
}

/* ---- grouping ---------------------------------------------------------- */

typedef struct {
    uint64_t hash;  /* canonical minimizer hash; UINT64_MAX when the read has no k-mer */
    int64_t pos;    /* read start relative to the minimizer, on the minimizer's strand */
    uint64_t idx;   /* read (or pair) index in input order */
} group_key;

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ull;
    x ^= x >> 27; x *= 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

static int base_code(uint8_t b)
{
    switch (b) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
    }
}

static void read_key(const uint8_t *s, uint64_t l, int k, group_key *out)
{
    const uint64_t mask = (k == 32) ? ~0ull : ((1ull << (2 * k)) - 1);
    uint64_t fw = 0, rv = 0;
    int valid = 0;
    out->hash = UINT64_MAX;
    out->pos = 0;
    for (uint64_t i = 0; i < l; i++) {
        int c = base_code(s[i]);
        if (c < 0) { valid = 0; fw = rv = 0; continue; }
        fw = ((fw << 2) | (uint64_t)c) & mask;
        rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
        if (++valid < k) continue;
        uint64_t start = i + 1 - (uint64_t)k;
        int forward = fw <= rv;
        uint64_t h = mix64(forward ? fw : rv);
        if (h < out->hash) {
            out->hash = h;
            /* Forward: the read starts `start` bases before the k-mer.
             * Reverse: on the k-mer's strand, the read starts
             * l - start - k bases before it. */
            out->pos = forward ? -(int64_t)start : -(int64_t)(l - start - (uint64_t)k);
        }
    }
}

static int key_cmp(const void *a, const void *b)
{
    const group_key *x = a, *y = b;
    if (x->hash != y->hash) return x->hash < y->hash ? -1 : 1;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return x->idx < y->idx ? -1 : (x->idx > y->idx);
}

/* ---- overlap chaining (--chain) ---------------------------------------- *
 *
 * Every read is indexed by its (w, k) window minimizers over both strands.
 * Reads are then laid out greedily: from a seed, walk left to the unvisited
 * read that starts the least distance before it, then from the seed walk
 * right to the one that starts the least distance after it (most overlap),
 * tracking each read's strand and start in the chain's frame. Minimizers
 * shared by more than --max-occ reads (repeats) are ignored; how the next
 * seed is chosen is described at chain_order. */

typedef struct {
    uint64_t h;
    uint32_t r;
    uint16_t off;   /* k-mer start in the read, read's own orientation */
    uint8_t s;      /* 1: the canonical k-mer is the read's forward k-mer */
    uint8_t pad;
} mz_ent;

typedef struct {
    uint32_t r;
    uint16_t off;
    uint8_t s;
    uint8_t pad;
} mz_hit;

typedef struct {
    uint32_t b;     /* bucket */
    uint16_t off;
    uint8_t s;
    uint8_t pad;
} mz_own;

static int ent_cmp(const void *a, const void *b)
{
    const mz_ent *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->r != y->r) return x->r < y->r ? -1 : 1;
    return (int)x->off - (int)y->off;
}

/* Window minimizers of one read; returns the number written to out. */
static size_t read_minimizers(const uint8_t *s, uint64_t l, int k, int w, uint32_t r,
                              mz_ent *out, uint64_t *hb, uint8_t *sb, uint32_t *dq)
{
    if (l > 65535) l = 65535;
    if (l < (uint64_t)k) return 0;
    const uint64_t mask = (k == 32) ? ~0ull : ((1ull << (2 * k)) - 1);
    size_t nk = (size_t)(l - (uint64_t)k + 1);
    uint64_t fw = 0, rv = 0;
    int valid = 0;
    for (uint64_t i = 0; i < l; i++) {
        int c = base_code(s[i]);
        if (c < 0) { valid = 0; fw = rv = 0; }
        else {
            fw = ((fw << 2) | (uint64_t)c) & mask;
            rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
            valid++;
        }
        if (i + 1 >= (uint64_t)k) {
            size_t p = (size_t)(i + 1 - (uint64_t)k);
            if (valid >= k) {
                sb[p] = fw <= rv;
                hb[p] = mix64(sb[p] ? fw : rv);
                if (hb[p] == UINT64_MAX) hb[p]--;
            } else {
                hb[p] = UINT64_MAX;
                sb[p] = 0;
            }
        }
    }
    size_t win = (size_t)w < nk ? (size_t)w : nk, head = 0, tail = 0, n = 0, last = (size_t)-1;
    for (size_t p = 0; p < nk; p++) {
        while (tail > head && hb[dq[tail - 1]] > hb[p]) tail--;
        dq[tail++] = (uint32_t)p;
        if (dq[head] + win <= p) head++;
        if (p + 1 >= win) {
            size_t m = dq[head];
            if (m != last && hb[m] != UINT64_MAX) {
                out[n].h = hb[m];
                out[n].r = r;
                out[n].off = (uint16_t)m;
                out[n].s = sb[m];
                out[n].pad = 0;
                n++;
                last = m;
            }
        }
    }
    return n;
}

typedef struct {
    uint32_t r;
    uint8_t o;
    int64_t P;
} cand;

typedef struct {
    const uint64_t *len;
    const uint32_t *bstart;     /* bucket b: hits[bstart[b] .. bstart[b+1]) */
    const mz_hit *hits;
    const uint64_t *ostart;     /* read r: own[ostart[r] .. ostart[r+1]) */
    const mz_own *own;
    uint64_t *visited;          /* bitset */
    int k;
    uint32_t max_occ;
    unsigned min_votes;         /* shared minimizers that must agree on the offset */
    cand *scratch;
    size_t scratch_cap;
} chain_ctx;

#define VISITED(x, r) (((x)->visited[(r) >> 6] >> ((r) & 63)) & 1)
#define SET_VISITED(x, r) ((x)->visited[(r) >> 6] |= 1ull << ((r) & 63))

static int cand_cmp(const void *a, const void *b)
{
    const cand *x = a, *y = b;
    if (x->r != y->r) return x->r < y->r ? -1 : 1;
    if (x->o != y->o) return (int)x->o - (int)y->o;
    return x->P < y->P ? -1 : (x->P > y->P);
}

enum { NB_LEFT, NB_RIGHT };

/* A neighbour of read cur (orientation o, chain start P). NB_RIGHT: the
 * unvisited read with the least shift >= 0; NB_LEFT: the unvisited read
 * with the greatest shift < 0. A neighbour counts only when at least
 * `votes` of cur's minimizers place it at the same start and strand, which
 * rejects overlaps through a single repeated k-mer. */
static int best_neighbour(const chain_ctx *x, uint32_t cur, int o, int64_t P, int mode, unsigned votes,
                          uint32_t *nr, int *no, int64_t *nP)
{
    size_t nc = 0;
    const int64_t k = x->k, lc = (int64_t)x->len[cur];
    for (uint64_t j = x->ostart[cur]; j < x->ostart[(uint64_t)cur + 1]; j++) {
        const mz_own *e = &x->own[j];
        uint32_t lo = x->bstart[e->b], hi = x->bstart[e->b + 1];
        if (hi - lo < 2 || hi - lo > x->max_occ) continue;
        int64_t kpos = o ? P + e->off : P + lc - e->off - k;
        int c = (o == e->s);
        for (uint32_t i = lo; i < hi && nc < x->scratch_cap; i++) {
            const mz_hit *h = &x->hits[i];
            if (h->r == cur || VISITED(x, h->r)) continue;
            int on = c ? h->s : !h->s;
            x->scratch[nc].r = h->r;
            x->scratch[nc].o = (uint8_t)on;
            x->scratch[nc].P = on ? kpos - h->off : kpos - ((int64_t)x->len[h->r] - h->off - k);
            nc++;
        }
    }
    if (votes > 1) qsort(x->scratch, nc, sizeof *x->scratch, cand_cmp);
    int found = 0;
    int64_t best = mode == NB_RIGHT ? INT64_MAX : INT64_MIN;
    for (size_t i = 0; i < nc;) {
        size_t g = i + 1;
        if (votes > 1)
            while (g < nc && x->scratch[g].r == x->scratch[i].r && x->scratch[g].o == x->scratch[i].o &&
                   x->scratch[g].P == x->scratch[i].P) g++;
        const cand *cd = &x->scratch[i];
        if (g - i >= votes) {
            int64_t sh = cd->P - P;
            int better;
            if (mode == NB_RIGHT) better = sh >= 0 && (sh < best || (sh == best && cd->r < *nr));
            else better = sh < 0 && (sh > best || (sh == best && cd->r < *nr));
            if (better) {
                best = sh;
                *nr = cd->r;
                *no = cd->o;
                *nP = cd->P;
                found = 1;
            }
        }
        i = g;
    }
    return found;
}

static unsigned g_min_votes = 2;
static long long g_backtrack_seeds, g_fallbacks, g_singletons;

/* Lay the reads out by overlap; fills order[0..n) and returns the number
 * of chains, or -1 on allocation failure.
 *
 * A new chain is seeded from the most recently placed read that still has
 * an unplaced neighbour (a stack of placed reads, popped once a read has
 * none left: placing never gives a read new unplaced neighbours, so each
 * read is examined about once). That is a depth-first walk of the overlap
 * graph that backtracks into the nearest unfinished region; the next read
 * in input order seeds only once every placed read's neighbourhood is
 * exhausted, i.e. a region not connected to anything placed. */
static long long chain_order(const uint8_t *seq, const uint64_t *off, const uint64_t *len, size_t n,
                             int k, int w, uint32_t max_occ, uint64_t *order, double *t_index)
{
    double t0 = now();
    size_t ecap = n * 8 + 1024, ne = 0;
    mz_ent *ents = malloc(ecap * sizeof *ents);
    uint64_t maxl = 0;
    for (size_t r = 0; r < n; r++) if (len[r] > maxl) maxl = len[r];
    if (maxl > 65535) maxl = 65535;
    uint64_t *hb = malloc((maxl + 1) * sizeof *hb);
    uint8_t *sb = malloc(maxl + 1);
    uint32_t *dq = malloc((maxl + 1) * sizeof *dq);
    mz_ent *tmp = malloc((maxl + 1) * sizeof *tmp);
    if (!ents || !hb || !sb || !dq || !tmp) return -1;
    for (size_t r = 0; r < n; r++) {
        size_t m = read_minimizers(seq + off[r], len[r], k, w, (uint32_t)r, tmp, hb, sb, dq);
        if (ne + m > ecap) {
            ecap = (ecap + m) * 3 / 2;
            ents = realloc(ents, ecap * sizeof *ents);
            if (!ents) return -1;
        }
        memcpy(ents + ne, tmp, m * sizeof *tmp);
        ne += m;
    }
    free(hb); free(sb); free(dq); free(tmp);
    qsort(ents, ne, sizeof *ents, ent_cmp);

    /* Buckets, per-bucket hits and per-read own minimizers. */
    size_t nb = 0;
    for (size_t i = 0; i < ne; i++) if (i == 0 || ents[i].h != ents[i - 1].h) nb++;
    uint32_t *bstart = malloc((nb + 1) * sizeof *bstart);
    mz_hit *hits = malloc((ne ? ne : 1) * sizeof *hits);
    uint64_t *ostart = calloc(n + 1, sizeof *ostart);
    mz_own *own = malloc((ne ? ne : 1) * sizeof *own);
    if (!bstart || !hits || !ostart || !own) return -1;
    for (size_t i = 0; i < ne; i++) ostart[ents[i].r + 1]++;
    for (size_t r = 0; r < n; r++) ostart[r + 1] += ostart[r];
    uint64_t *fill = malloc((n + 1) * sizeof *fill);
    if (!fill) return -1;
    memcpy(fill, ostart, (n + 1) * sizeof *fill);
    size_t b = (size_t)-1;
    for (size_t i = 0; i < ne; i++) {
        if (i == 0 || ents[i].h != ents[i - 1].h) bstart[++b] = (uint32_t)i;
        hits[i].r = ents[i].r;
        hits[i].off = ents[i].off;
        hits[i].s = ents[i].s;
        hits[i].pad = 0;
        mz_own *o = &own[fill[ents[i].r]++];
        o->b = (uint32_t)b;
        o->off = ents[i].off;
        o->s = ents[i].s;
        o->pad = 0;
    }
    bstart[nb] = (uint32_t)ne;
    free(fill);
    free(ents);
    *t_index = now() - t0;

    chain_ctx x = { len, bstart, hits, ostart, own, calloc((n + 63) / 64, sizeof(uint64_t)), k, max_occ,
                    g_min_votes, malloc((1u << 16) * sizeof(cand)), 1u << 16 };
    uint64_t *left = malloc(n * sizeof *left);
    uint64_t *deferred = malloc(n * sizeof *deferred);
    uint32_t *stack = malloc(n * sizeof *stack);   /* placed reads, most recent on top */
    size_t nd = 0, ns = 0;
    if (!x.visited || !left || !x.scratch || !deferred || !stack) return -1;
    size_t no = 0, next_input = 0;
    long long chains = 0;
    for (;;) {
        size_t seed = (size_t)-1;
        while (ns && seed == (size_t)-1) {
            uint32_t r = stack[ns - 1], nr2 = 0;
            int on2 = 1;
            int64_t nP2 = 0;
            if (best_neighbour(&x, r, 1, 0, NB_RIGHT, 1, &nr2, &on2, &nP2) ||
                best_neighbour(&x, r, 1, 0, NB_LEFT, 1, &nr2, &on2, &nP2)) {
                seed = nr2;
                g_backtrack_seeds++;
            } else {
                ns--;
            }
        }
        if (seed == (size_t)-1) {
            while (next_input < n && VISITED(&x, next_input)) next_input++;
            if (next_input == n) break;
            seed = next_input;
            g_fallbacks++;
        }
        chains++;
        SET_VISITED(&x, seed);
        /* Left of the seed first (stored in reverse), then the seed, then right. */
        size_t nl = 0;
        uint32_t cur = (uint32_t)seed, nr = 0;
        int o = 1, on = 1;
        int64_t P = 0, nP = 0;
        while (best_neighbour(&x, cur, o, P, NB_LEFT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            left[nl++] = nr;
            cur = nr; o = on; P = nP;
        }
        size_t chain_start = no;
        while (nl) order[no++] = left[--nl];
        order[no++] = seed;
        cur = (uint32_t)seed; o = 1; P = 0;
        while (best_neighbour(&x, cur, o, P, NB_RIGHT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            order[no++] = nr;
            cur = nr; o = on; P = nP;
        }
        if (no - chain_start == 1) {
            /* A read with no confident neighbour goes to the end of the
             * run, so it does not scatter into the local blocks. */
            g_singletons++;
            deferred[nd++] = order[--no];
        } else {
            for (size_t i = chain_start; i < no; i++) stack[ns++] = (uint32_t)order[i];
        }
    }
    memcpy(order + no, deferred, nd * sizeof *order);
    no += nd;
    free(stack);
    free(deferred);
    free(left);
    free(x.scratch);
    free(x.visited);
    free(bstart);
    free(hits);
    free(ostart);
    free(own);
    return chains;
}

/* Diagnostic for --pos: how local the grouped order is. Reads the mapping
 * position of each read (input order, 0 = unmapped) and reports the share
 * of consecutive mapped reads within 1 kb, the mean number of distinct
 * 100 kb bins each block touches, and the share of each block's mapped
 * reads that fall in its 20 busiest bins (coordinate order: ~10 bins and
 * 100%). The bin count is swamped by a few stray reads; the busiest-bins
 * share is the one to read. */
static void report_locality(const char *path, const uint64_t *order, size_t n_units, size_t unit,
                            const uint64_t *len, size_t n, unsigned long long block_bases)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return; }
    int64_t *pos = malloc(n * sizeof *pos);
    size_t m = 0;
    long long v;
    while (m < n && fscanf(f, "%lld", &v) == 1) pos[m++] = v;
    fclose(f);
    if (m != n) { fprintf(stderr, "%s: %zu positions for %zu reads\n", path, m, n); free(pos); return; }
    enum { NBINS = 4096, TOP = 20 };
    unsigned char seen[NBINS];
    static unsigned cnt[NBINS];
    memset(seen, 0, sizeof seen);
    memset(cnt, 0, sizeof cnt);
    unsigned long long in_top = 0, mapped = 0, b_mapped = 0;
    unsigned long long pairs = 0, near = 0, blocks = 0, bins_total = 0, bins = 0, b_used = 0;
    int64_t prev = 0;
    for (size_t u = 0; u < n_units; u++) {
        for (size_t k = 0; k < unit; k++) {
            size_t r = order[u] * unit + k;
            int64_t p = pos[r];
            if (p > 0) {
                if (prev > 0) { pairs++; if (llabs(p - prev) < 1000) near++; }
                size_t bin = (size_t)(p / 100000);
                if (bin >= NBINS) bin = NBINS - 1;
                if (!seen[bin]) { seen[bin] = 1; bins++; }
                cnt[bin]++;
                b_mapped++;
            }
            prev = p;
            b_used += len[r];
        }
        if (b_used >= block_bases || u + 1 == n_units) {
            blocks++; bins_total += bins; bins = 0; b_used = 0;
            for (int t = 0; t < TOP; t++) {
                size_t bi = 0;
                for (size_t q = 1; q < NBINS; q++) if (cnt[q] > cnt[bi]) bi = q;
                in_top += cnt[bi];
                cnt[bi] = 0;
            }
            mapped += b_mapped;
            b_mapped = 0;
            memset(seen, 0, sizeof seen);
            memset(cnt, 0, sizeof cnt);
        }
    }
    printf("locality: consecutive mapped reads within 1 kb %.2f%%, mean 100 kb bins per block %.1f, "
           "reads in each block's %d busiest bins %.2f%% (%llu blocks)\n",
           pairs ? 100.0 * (double)near / (double)pairs : 0.0, blocks ? (double)bins_total / (double)blocks : 0.0,
           (int)TOP, mapped ? 100.0 * (double)in_top / (double)mapped : 0.0, blocks);
    free(pos);
}

int main(int argc, char **argv)
{
    coder c;
    memset(&c, 0, sizeof c);
    ttio_seq_cm_default_params(&c.p);
    unsigned long long block_bases = 64ull << 20;
    int ai = 1, n_orders = 0, group = 0, chain = 0, as_is = 0, no_code = 0, pairs = 0, group_k = 20, group_w = 32;
    unsigned max_occ = 256;
    const char *pos_path = NULL;  /* diagnostic: one mapping position per read, input order */
    while (ai < argc) {
        if (!strcmp(argv[ai], "--block-bases") && ai + 1 < argc) block_bases = strtoull(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "--table-bits") && ai + 1 < argc) c.p.table_bits = (uint8_t)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--no-rc")) c.p.flags = 0;
        else if (!strcmp(argv[ai], "--group")) group = 1;
        else if (!strcmp(argv[ai], "--chain")) group = chain = 1;
        else if (!strcmp(argv[ai], "--as-is")) group = as_is = 1;
        else if (!strcmp(argv[ai], "--no-code")) no_code = 1;
        else if (!strcmp(argv[ai], "--pairs")) pairs = 1;
        else if (!strcmp(argv[ai], "--group-k") && ai + 1 < argc) group_k = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--group-w") && ai + 1 < argc) group_w = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--max-occ") && ai + 1 < argc) max_occ = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--min-votes") && ai + 1 < argc) g_min_votes = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--pos") && ai + 1 < argc) pos_path = argv[++ai];
        else if (n_orders < 3) c.p.orders[n_orders++] = (uint8_t)atoi(argv[ai]);
        ai++;
    }
    if (n_orders) c.p.n_orders = (uint8_t)n_orders;
    if (group_k < 8 || group_k > 32) { fprintf(stderr, "--group-k must be 8 to 32\n"); return 2; }
    if (group_w < 1) { fprintf(stderr, "--group-w must be at least 1\n"); return 2; }
    if (chain && pairs) { fprintf(stderr, "--chain orders single reads; --pairs is not supported with it\n"); return 2; }
    long long n_chains = 0;
    double t_index = 0;

    size_t cap = 1 << 26, used = 0, rcap = 1 << 20, n = 0;
    uint8_t *seq = malloc(cap);
    uint64_t *len = malloc(rcap * sizeof *len);
    char *line = NULL;
    size_t lcap = 0;
    ssize_t l;
    unsigned long long ln = 0;
    int done = 0;
    double perm_ideal_bits = 0, perm_fixed_bits = 0, t_group = 0;
    while (!done) {
        l = getline(&line, &lcap, stdin);
        if (l > 0 && (ln++ & 3) == 1) {
            while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) l--;
            while (used + (size_t)l > cap) { cap *= 2; seq = realloc(seq, cap); }
            if (n == rcap) { rcap *= 2; len = realloc(len, rcap * sizeof *len); }
            memcpy(seq + used, line, (size_t)l);
            used += (size_t)l;
            len[n++] = (uint64_t)l;
        }
        if (l <= 0) done = 1;
        if (!group && (used >= block_bases || done) && n) {
            if (code_block(&c, seq, len, n, used)) return 1;
            used = 0;
            n = 0;
        }
    }

    if (group && n) {
        double t0 = now();
        if (pairs && (n & 1)) { fprintf(stderr, "--pairs needs an even read count\n"); return 2; }
        size_t unit = pairs ? 2 : 1, n_units = n / unit;
        uint64_t *off = malloc((n + 1) * sizeof *off);
        off[0] = 0;
        for (size_t i = 0; i < n; i++) off[i + 1] = off[i] + len[i];
        uint64_t *order = malloc(n_units * sizeof *order);
        if (as_is) {
            for (size_t u = 0; u < n_units; u++) order[u] = u;
        } else if (chain) {
            n_chains = chain_order(seq, off, len, n, group_k, group_w, max_occ, order, &t_index);
            if (n_chains < 0) { fprintf(stderr, "out of memory while chaining\n"); return 1; }
            printf("chain k %d w %d max-occ %u min-votes %u: %lld chains (%.1f reads each; "
                   "%lld single reads, %lld backtrack seeds, %lld input-order seeds), index %.1f s\n",
                   group_k, group_w, max_occ, g_min_votes, n_chains,
                   (double)n / (double)(n_chains ? n_chains : 1), g_singletons, g_backtrack_seeds, g_fallbacks, t_index);
        } else {
            group_key *keys = malloc(n_units * sizeof *keys);
            for (size_t u = 0; u < n_units; u++) {
                group_key best;
                read_key(seq + off[u * unit], len[u * unit], group_k, &best);
                for (size_t m = 1; m < unit; m++) {
                    group_key k2;
                    read_key(seq + off[u * unit + m], len[u * unit + m], group_k, &k2);
                    if (k2.hash < best.hash) best = k2;
                }
                best.idx = u;
                keys[u] = best;
            }
            qsort(keys, n_units, sizeof *keys, key_cmp);
            for (size_t u = 0; u < n_units; u++) order[u] = keys[u].idx;
            free(keys);
        }
        if (pos_path) report_locality(pos_path, order, n_units, unit, len, n, block_bases);
        uint8_t *gseq = malloc(used ? used : 1);
        uint64_t *glen = malloc(n * sizeof *glen);
        size_t gused = 0, gn = 0;
        for (size_t u = 0; u < n_units; u++) {
            for (size_t m = 0; m < unit; m++) {
                size_t r = order[u] * unit + m;
                memcpy(gseq + gused, seq + off[r], len[r]);
                gused += len[r];
                glen[gn++] = len[r];
            }
        }
        if (gused != used || gn != n) { fprintf(stderr, "grouping lost reads (%zu/%zu bases)\n", gused, used); return 1; }
        free(order);
        free(off);
        free(seq);
        free(len);
        t_group = now() - t0;
        perm_ideal_bits = lgamma((double)n_units + 1.0) / log(2.0);
        perm_fixed_bits = (double)n_units * ceil(log2((double)(n_units > 1 ? n_units : 2)));

        size_t b_start = 0, r_start = 0, b_used = 0;
        for (size_t r = 0; r < gn; r++) {
            b_used += glen[r];
            /* Cut on a unit boundary so pairs never split across blocks. */
            if ((b_used >= block_bases && (r + 1) % unit == 0) || r + 1 == gn) {
                if (!no_code && code_block(&c, gseq + b_start, glen + r_start, r + 1 - r_start, b_used)) return 1;
                b_start += b_used;
                r_start = r + 1;
                b_used = 0;
            }
        }
        free(gseq);
        free(glen);
    }

    printf("orders");
    for (int i = 0; i < c.p.n_orders; i++) printf(" %u", c.p.orders[i]);
    printf("%s table_bits %u block_bases %llu: blocks %llu bases %llu bytes %llu  bits/base %.4f  "
           "encode %.1f MB/s decode %.1f MB/s\n",
           c.p.flags ? " +rc" : "", c.p.table_bits, block_bases, c.blocks, c.total_bases, c.total_bytes,
           c.total_bases ? 8.0 * (double)c.total_bytes / (double)c.total_bases : 0.0,
           c.total_bases / 1e6 / c.t_enc, c.total_bases / 1e6 / c.t_dec);
    if (group && c.total_bases) {
        double b = (double)c.total_bases;
        printf("group k %d%s: grouping %.1f s  permutation ideal %.0f bytes (%.4f bits/base), "
               "fixed-width %.0f bytes  total with ideal permutation %.4f bits/base, "
               "with fixed-width %.4f bits/base\n",
               group_k, pairs ? " pairs" : "", t_group,
               perm_ideal_bits / 8, perm_ideal_bits / b, perm_fixed_bits / 8,
               (8.0 * (double)c.total_bytes + perm_ideal_bits) / b,
               (8.0 * (double)c.total_bytes + perm_fixed_bits) / b);
    }
    return 0;
}
