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
 *                      [--group-w W (0: from read length)] [--max-occ N]
 *                      [--min-votes V] [--seed-votes V] [--fill] [--mates]
 *                      [--scaffold [--scaffold-min N] [--cluster-cap BASES]]
 *                      | --layout ... | --as-is]
 *                     [--lib]   (ttio_seq_group from the library, defaults)
 *                     [--pos POSITIONS] [--dump-order FILE] [--no-code]
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

/* Sort by hash with an LSD radix sort, 16 bits a pass. Stable, and the
 * entries arrive in (read, offset) order, so equal hashes keep that order.
 * Returns 0, or -1 when the buffer cannot be allocated. */
static int sort_ents(mz_ent *a, size_t n)
{
    mz_ent *t = malloc((n ? n : 1) * sizeof *t);
    size_t *cnt = malloc(((size_t)1 << 16) * sizeof *cnt);
    if (!t || !cnt) { free(t); free(cnt); return -1; }
    for (int shift = 0; shift < 64; shift += 16) {
        memset(cnt, 0, ((size_t)1 << 16) * sizeof *cnt);
        for (size_t i = 0; i < n; i++) cnt[(a[i].h >> shift) & 0xFFFF]++;
        size_t sum = 0;
        for (size_t d = 0; d < ((size_t)1 << 16); d++) { size_t c = cnt[d]; cnt[d] = sum; sum += c; }
        for (size_t i = 0; i < n; i++) t[cnt[(a[i].h >> shift) & 0xFFFF]++] = a[i];
        mz_ent *sw = a; a = t; t = sw;
    }
    /* Four passes: the sorted data is back in the caller's buffer. */
    free(t);
    free(cnt);
    return 0;
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
    int64_t P;      /* start in the chain's frame */
    uint32_t r;     /* UINT32_MAX: empty slot */
    uint8_t o;
    uint16_t votes;
} cand;

enum { VOTE_BITS = 17 };    /* vote table slots: 2^17, kept under half full */

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
    cand *vt;                   /* vote table, open addressing on (r, o, P) */
    uint32_t *used;             /* its occupied slots, in insertion order */
} chain_ctx;

#define VISITED(x, r) (((x)->visited[(r) >> 6] >> ((r) & 63)) & 1)
#define SET_VISITED(x, r) ((x)->visited[(r) >> 6] |= 1ull << ((r) & 63))

enum { NB_LEFT, NB_RIGHT };

/* Count, for every unvisited read sharing a minimizer with read cur
 * (orientation o, start P in the frame), how many of cur's minimizers place
 * it at each start and strand. Fills x->vt and lists the slots in x->used;
 * returns their number; with `all`, placed reads count too. clear_tally
 * empties the table again. */
static size_t tally(const chain_ctx *x, uint32_t cur, int o, int64_t P, int all)
{
    size_t nu = 0;
    const size_t mask = ((size_t)1 << VOTE_BITS) - 1, cap = (size_t)1 << (VOTE_BITS - 1);
    const int64_t k = x->k, lc = (int64_t)x->len[cur];
    for (uint64_t j = x->ostart[cur]; j < x->ostart[(uint64_t)cur + 1]; j++) {
        const mz_own *e = &x->own[j];
        uint32_t lo = x->bstart[e->b], hi = x->bstart[e->b + 1];
        if (hi - lo < 2 || hi - lo > x->max_occ) continue;
        int64_t kpos = o ? P + e->off : P + lc - e->off - k;
        int c = (o == e->s);
        for (uint32_t i = lo; i < hi && nu < cap; i++) {
            const mz_hit *h = &x->hits[i];
            if (h->r == cur || (!all && VISITED(x, h->r))) continue;
            uint8_t on = (uint8_t)(c ? h->s : !h->s);
            int64_t Pn = on ? kpos - h->off : kpos - ((int64_t)x->len[h->r] - h->off - k);
            size_t slot = mix64(((uint64_t)h->r << 1 | on) ^ ((uint64_t)Pn * 0x9E3779B97F4A7C15ull)) & mask;
            for (;;) {
                cand *v = &x->vt[slot];
                if (v->r == UINT32_MAX) {
                    v->r = h->r; v->o = on; v->P = Pn; v->votes = 1;
                    x->used[nu++] = (uint32_t)slot;
                    break;
                }
                if (v->r == h->r && v->o == on && v->P == Pn) {
                    if (v->votes < UINT16_MAX) v->votes++;
                    break;
                }
                slot = (slot + 1) & mask;
            }
        }
    }
    return nu;
}

static void clear_tally(const chain_ctx *x, size_t nu)
{
    for (size_t i = 0; i < nu; i++) x->vt[x->used[i]].r = UINT32_MAX;
}

/* A neighbour of read cur (orientation o, chain start P). NB_RIGHT: the
 * unvisited read with the least shift >= 0; NB_LEFT: the unvisited read
 * with the greatest shift < 0. A neighbour counts only when at least
 * `votes` of cur's minimizers place it at the same start and strand, which
 * rejects overlaps through a single repeated k-mer. */
static int best_neighbour(const chain_ctx *x, uint32_t cur, int o, int64_t P, int mode, unsigned votes,
                          uint32_t *nr, int *no, int64_t *nP)
{
    size_t nu = tally(x, cur, o, P, 0);
    int found = 0;
    int64_t best = mode == NB_RIGHT ? INT64_MAX : INT64_MIN;
    for (size_t i = 0; i < nu; i++) {
        cand *cd = &x->vt[x->used[i]];
        if (cd->votes >= votes) {
            int64_t sh = cd->P - P;
            int better;
            /* Ties go to the lower read, then orientation 0, so the layout
             * does not depend on the vote table's slot order. */
            int tie = sh == best && (cd->r < *nr || (cd->r == *nr && cd->o < *no));
            if (mode == NB_RIGHT) better = sh >= 0 && (sh < best || tie);
            else better = sh < 0 && (sh > best || tie);
            if (better) {
                best = sh;
                *nr = cd->r;
                *no = cd->o;
                *nP = cd->P;
                found = 1;
            }
        }
    }
    clear_tally(x, nu);
    return found;
}

static unsigned g_min_votes = 2;
static unsigned g_seed_votes = 1;   /* votes a backtrack seed needs (--seed-votes) */
static int g_fill;                  /* --fill: place weakly linked reads beside the chain */
static long long g_filled;
/* --mates: reads with the same name (any /1 or /2 suffix dropped) are
 * mates; g_mate[r] is r's mate or UINT32_MAX. When a chain ends, the next seed is
 * the unplaced mate of one of the last MATE_LOOKBACK placed reads, newest
 * first: a mate lies one fragment further along the genome, across repeats
 * that no unique overlap spans. */
static int g_mates;
static const uint32_t *g_mate;
static int g_scaffold;              /* --scaffold: join chains into paths before emitting */
static unsigned g_scaffold_min = 4; /* links two chains need to be joined (--scaffold-min) */
static uint64_t g_cluster_cap = 16ull << 20;  /* bases per cluster (--cluster-cap) */
static long long g_scaffold_joins, g_scaffold_paths;
static long long g_mate_seeds;
enum { MATE_LOOKBACK = 4096 };
static long long g_backtrack_seeds, g_fallbacks, g_singletons;
/* Diagnostic, with --pos: each read's mapping position, and per seed kind
 * (0 backtrack, 1 input order, 2 mate) the chains of 2+ reads it started and how
 * many of those began more than 10 kb from the read placed before them. */
static const int64_t *g_pos;
static long long g_kind_chains[3], g_kind_far[3], g_kind_mapped[3];

typedef struct {
    int64_t c;      /* coordinate in the component's frame */
    uint32_t r;
} placed;

static int placed_cmp(const void *a, const void *b)
{
    const placed *x = a, *y = b;
    if (x->c != y->c) return x->c < y->c ? -1 : 1;
    return x->r < y->r ? -1 : (x->r > y->r);
}

static int hit_cmp(const void *a, const void *b)
{
    const cand *x = a, *y = b;
    if (x->r != y->r) return x->r < y->r ? -1 : 1;
    if (x->votes != y->votes) return x->votes > y->votes ? -1 : 1;
    if (x->o != y->o) return (int)x->o - (int)y->o;
    return x->P < y->P ? -1 : (x->P > y->P);
}

static int shift_cmp(const void *a, const void *b)
{
    const cand *x = a, *y = b;
    if (x->P != y->P) return x->P < y->P ? -1 : 1;
    return x->r < y->r ? -1 : (x->r > y->r);
}

/* --layout: give every read a coordinate instead of emitting reads in walk
 * order. A breadth-first search from a seed follows every confident overlap
 * (min_votes agreeing minimizers), carrying each read's strand and start
 * from the read that reached it; the component is then sorted by start. A
 * read the walk would skip lands where its coordinate puts it, so no holes
 * are left to backtrack into. Components of one read go to the end. */
static long long layout_order(chain_ctx *x, size_t n, uint64_t *order)
{
    uint32_t *queue = malloc((n ? n : 1) * sizeof *queue);
    int64_t *crd = malloc((n ? n : 1) * sizeof *crd);
    uint8_t *ori = malloc(n ? n : 1);
    placed *pl = malloc((n ? n : 1) * sizeof *pl);
    uint64_t *deferred = malloc((n ? n : 1) * sizeof *deferred);
    cand *hits = malloc(((size_t)1 << (VOTE_BITS - 1)) * sizeof *hits);
    if (!queue || !crd || !ori || !pl || !deferred || !hits) return -1;
    size_t no = 0, nd = 0;
    long long comps = 0;
    for (size_t seed = 0; seed < n; seed++) {
        if (VISITED(x, seed)) continue;
        comps++;
        SET_VISITED(x, seed);
        size_t qh = 0, qt = 0;
        queue[qt++] = (uint32_t)seed;
        crd[seed] = 0;
        ori[seed] = 1;
        while (qh < qt) {
            uint32_t cur = queue[qh++];
            size_t nu = tally(x, cur, ori[cur], crd[cur], 0), nh = 0;
            for (size_t i = 0; i < nu; i++) {
                const cand *cd = &x->vt[x->used[i]];
                if (cd->votes >= x->min_votes) hits[nh++] = *cd;
            }
            clear_tally(x, nu);
            /* Best-supported placement per read, deterministically. */
            qsort(hits, nh, sizeof *hits, hit_cmp);
            for (size_t i = 0; i < nh; i++) {
                if (i && hits[i].r == hits[i - 1].r) continue;
                uint32_t r = hits[i].r;
                SET_VISITED(x, r);
                crd[r] = hits[i].P;
                ori[r] = hits[i].o;
                queue[qt++] = r;
            }
        }
        if (qt == 1) {
            g_singletons++;
            deferred[nd++] = seed;
            continue;
        }
        g_kind_chains[1]++;
        for (size_t i = 0; i < qt; i++) { pl[i].c = crd[queue[i]]; pl[i].r = queue[i]; }
        qsort(pl, qt, sizeof *pl, placed_cmp);
        for (size_t i = 0; i < qt; i++) order[no++] = pl[i].r;
    }
    memcpy(order + no, deferred, nd * sizeof *order);
    free(queue); free(crd); free(ori); free(pl); free(deferred); free(hits);
    return comps;
}

static int u64_cmp(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y);
}

typedef struct { uint32_t a, b, w; } seg_edge;

static int edge_cmp(const void *p, const void *q)
{
    const seg_edge *x = p, *y = q;
    if (x->w != y->w) return x->w > y->w ? -1 : 1;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    return x->b < y->b ? -1 : (x->b > y->b);
}

static uint32_t uf_find(uint32_t *par, uint32_t a)
{
    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
    return a;
}

/* --scaffold: the chains are local, but their sequence wanders over the
 * whole chromosome, so a region's reads end up spread over every block. At
 * high coverage a region is covered by many parallel chains, not a line of
 * them, so the chains are clustered rather than joined into paths. Count
 * links between chains -- a confident overlap (min_votes) between reads of
 * two chains, or a mate pair across them, worth 2 -- and merge clusters
 * along the strongest links (union-find, by descending weight). Links
 * through a repeat are few and links through unique overlap many, so
 * --scaffold-min drops the weak ones, and --cluster-cap bounds a cluster's
 * bases. Each cluster is emitted whole. seg[0..nseg] holds each chain's
 * first position in order; co/cP each read's strand and start in its chain. */
static int scaffold_order(const chain_ctx *x, uint64_t *order, size_t n, const uint64_t *seg, size_t nseg,
                          const uint8_t *co, const int64_t *cP)
{
    uint32_t *cid = malloc((n ? n : 1) * sizeof *cid);
    size_t lcap = n * 4 + 1024, nl = 0;
    uint64_t *links = malloc(lcap * sizeof *links);
    if (!cid || !links) return -1;
    for (size_t sg = 0; sg < nseg; sg++)
        for (uint64_t i = seg[sg]; i < seg[sg + 1]; i++) cid[order[i]] = (uint32_t)sg;
    for (size_t r = 0; r < n; r++) {
        size_t nu = tally(x, (uint32_t)r, co[r], cP[r], 1);
        for (size_t j = 0; j < nu; j++) {
            const cand *cd = &x->vt[x->used[j]];
            if (cd->votes < x->min_votes || cd->r <= r || cid[cd->r] == cid[r]) continue;
            uint32_t a = cid[r], b = cid[cd->r];
            if (a > b) { uint32_t t = a; a = b; b = t; }
            if (nl == lcap) { lcap = lcap * 3 / 2; links = realloc(links, lcap * sizeof *links); if (!links) return -1; }
            links[nl++] = (uint64_t)a << 32 | b;
        }
        clear_tally(x, nu);
        if (g_mate && g_mate[r] != UINT32_MAX && g_mate[r] > r && cid[g_mate[r]] != cid[r]) {
            uint32_t a = cid[r], b = cid[g_mate[r]];
            if (a > b) { uint32_t t = a; a = b; b = t; }
            for (int t = 0; t < 2; t++) {
                if (nl == lcap) { lcap = lcap * 3 / 2; links = realloc(links, lcap * sizeof *links); if (!links) return -1; }
                links[nl++] = (uint64_t)a << 32 | b;
            }
        }
    }
    qsort(links, nl, sizeof *links, u64_cmp);
    size_t ne = 0;
    seg_edge *edges = malloc((nl ? nl : 1) * sizeof *edges);
    if (!edges) return -1;
    for (size_t i = 0; i < nl;) {
        size_t j = i + 1;
        while (j < nl && links[j] == links[i]) j++;
        if (j - i >= g_scaffold_min)
            edges[ne++] = (seg_edge){ (uint32_t)(links[i] >> 32), (uint32_t)links[i], (uint32_t)(j - i) };
        i = j;
    }
    free(links);
    qsort(edges, ne, sizeof *edges, edge_cmp);
    /* Merge along the strongest links first; a merge that would take a
     * cluster past g_cluster_cap bases is refused, so links through repeats
     * cannot snowball a cluster across the chromosome. */
    uint32_t *par = malloc((nseg ? nseg : 1) * sizeof *par);
    uint64_t *bases = malloc((nseg ? nseg : 1) * sizeof *bases);
    if (!par || !bases) return -1;
    for (size_t sg = 0; sg < nseg; sg++) {
        par[sg] = (uint32_t)sg;
        bases[sg] = 0;
        for (uint64_t i = seg[sg]; i < seg[sg + 1]; i++) bases[sg] += x->len[order[i]];
    }
    for (size_t e = 0; e < ne; e++) {
        uint32_t ra = uf_find(par, edges[e].a), rb = uf_find(par, edges[e].b);
        if (ra == rb || bases[ra] + bases[rb] > g_cluster_cap) continue;
        par[ra] = rb;
        bases[rb] += bases[ra];
        g_scaffold_joins++;
    }
    free(edges);
    /* Emit clusters whole, in the order their first chain was laid out;
     * within a cluster, chains keep their order. */
    uint32_t *root = malloc((nseg ? nseg : 1) * sizeof *root);
    uint32_t *head = malloc((nseg ? nseg : 1) * sizeof *head), *next = malloc((nseg ? nseg : 1) * sizeof *next);
    uint32_t *tail = malloc((nseg ? nseg : 1) * sizeof *tail);
    uint64_t *out = malloc((n ? n : 1) * sizeof *out);
    uint8_t *done = calloc(nseg ? nseg : 1, 1);
    if (!root || !head || !next || !tail || !out || !done) return -1;
    for (size_t sg = 0; sg < nseg; sg++) { head[sg] = tail[sg] = next[sg] = UINT32_MAX; }
    for (size_t sg = 0; sg < nseg; sg++) {
        uint32_t r = root[sg] = uf_find(par, (uint32_t)sg);
        if (head[r] == UINT32_MAX) head[r] = (uint32_t)sg; else next[tail[r]] = (uint32_t)sg;
        tail[r] = (uint32_t)sg;
    }
    size_t no = 0;
    for (size_t sg = 0; sg < nseg; sg++) {
        uint32_t r = root[sg];
        if (done[r]) continue;
        done[r] = 1;
        g_scaffold_paths++;
        for (uint32_t c = head[r]; c != UINT32_MAX; c = next[c])
            for (uint64_t i = seg[c]; i < seg[c + 1]; i++) out[no++] = order[i];
    }
    if (no != n) { fprintf(stderr, "clustering lost reads (%zu of %zu)\n", no, n); return -1; }
    memcpy(order, out, n * sizeof *order);
    free(out); free(done); free(root); free(head); free(next); free(tail); free(par); free(bases); free(cid);
    return 0;
}

/* Lay the reads out by overlap; fills order[0..n) and returns the number
 * of chains (components with --layout), or -1 on allocation failure.
 *
 * With --mates, a new chain is seeded first from the unplaced mate of one
 * of the last MATE_LOOKBACK placed reads. Otherwise, or failing that, from
 * the most recently placed read that still has an unplaced neighbour (a
 * stack of placed reads, popped once a read has none left: placing never
 * gives a read new unplaced neighbours, so each read is examined about
 * once). That is a depth-first walk of the overlap graph that backtracks
 * into the nearest unfinished region; the next read in input order seeds
 * only once every placed read's neighbourhood is exhausted, i.e. a region
 * not connected to anything placed. --fill and --scaffold act on the
 * chains as they are made and once all are made. */
static long long chain_order(const uint8_t *seq, const uint64_t *off, const uint64_t *len, size_t n,
                             int k, int w, uint32_t max_occ, int layout, uint64_t *order, double *t_index)
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
    if (sort_ents(ents, ne)) return -1;

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
                    g_min_votes, malloc(((size_t)1 << VOTE_BITS) * sizeof(cand)),
                    malloc(((size_t)1 << (VOTE_BITS - 1)) * sizeof(uint32_t)) };
    uint64_t *left = malloc(n * sizeof *left);
    uint64_t *deferred = malloc(n * sizeof *deferred);
    uint32_t *stack = malloc(n * sizeof *stack);   /* placed reads, most recent on top */
    /* Each placed read's strand and start in its chain's frame, the chain as
     * walked, and the weak neighbours --fill places beside a read. */
    uint8_t *co = malloc(n ? n : 1);
    int64_t *cP = malloc((n ? n : 1) * sizeof *cP);
    uint64_t *walked = malloc((n ? n : 1) * sizeof *walked);
    cand *weak = malloc(((size_t)1 << (VOTE_BITS - 1)) * sizeof *weak);
    uint64_t *seg = malloc((n + 1) * sizeof *seg);     /* each chain's first position in order */
    size_t nd = 0, ns = 0, nseg = 0;
    if (!x.visited || !left || !x.vt || !x.used || !deferred || !stack || !co || !cP || !walked || !weak || !seg) return -1;
    for (size_t i = 0; i < ((size_t)1 << VOTE_BITS); i++) x.vt[i].r = UINT32_MAX;
    if (layout) {
        long long comps = layout_order(&x, n, order);
        free(stack); free(deferred); free(left);
        free(x.vt); free(x.used); free(x.visited);
        free(bstart); free(hits); free(ostart); free(own);
        return comps;
    }
    size_t no = 0, next_input = 0;
    long long chains = 0;
    for (;;) {
        size_t seed = (size_t)-1;
        int kind = 0;
        if (g_mates)
            for (size_t back = 0; back < MATE_LOOKBACK && back < no && seed == (size_t)-1; back++) {
                uint32_t m = g_mate[order[no - 1 - back]];
                if (m != UINT32_MAX && !VISITED(&x, m)) { seed = m; kind = 2; g_mate_seeds++; }
            }
        while (ns && seed == (size_t)-1) {
            uint32_t r = stack[ns - 1], nr2 = 0;
            int on2 = 1;
            int64_t nP2 = 0;
            if (best_neighbour(&x, r, 1, 0, NB_RIGHT, g_seed_votes, &nr2, &on2, &nP2) ||
                best_neighbour(&x, r, 1, 0, NB_LEFT, g_seed_votes, &nr2, &on2, &nP2)) {
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
            kind = 1;
            g_fallbacks++;
        }
        chains++;
        SET_VISITED(&x, seed);
        /* Left of the seed first (stored in reverse), then the seed, then right. */
        size_t nl = 0;
        uint32_t cur = (uint32_t)seed, nr = 0;
        int o = 1, on = 1;
        int64_t P = 0, nP = 0;
        co[seed] = 1;
        cP[seed] = 0;
        while (best_neighbour(&x, cur, o, P, NB_LEFT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            left[nl++] = nr;
            co[nr] = (uint8_t)on; cP[nr] = nP;
            cur = nr; o = on; P = nP;
        }
        size_t chain_start = no;
        while (nl) order[no++] = left[--nl];
        order[no++] = seed;
        cur = (uint32_t)seed; o = 1; P = 0;
        while (best_neighbour(&x, cur, o, P, NB_RIGHT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            order[no++] = nr;
            co[nr] = (uint8_t)on; cP[nr] = nP;
            cur = nr; o = on; P = nP;
        }
        if (g_fill) {
            /* After each chain read, its unplaced neighbours on even one
             * minimizer that start within a read length of it, by start: a
             * read with errors that no step reached sits beside the reads
             * it overlaps instead of becoming a hole or a single. */
            size_t nw = no - chain_start;
            memcpy(walked, order + chain_start, nw * sizeof *walked);
            no = chain_start;
            for (size_t i = 0; i < nw; i++) {
                uint32_t r = (uint32_t)walked[i];
                order[no++] = r;
                size_t nu = tally(&x, r, co[r], cP[r], 0), nk = 0;
                for (size_t j = 0; j < nu; j++) {
                    const cand *cd = &x.vt[x.used[j]];
                    if (llabs(cd->P - cP[r]) <= (int64_t)len[r]) weak[nk++] = *cd;
                }
                clear_tally(&x, nu);
                qsort(weak, nk, sizeof *weak, hit_cmp);
                size_t nk2 = 0;
                for (size_t j = 0; j < nk; j++)
                    if (!j || weak[j].r != weak[j - 1].r) weak[nk2++] = weak[j];
                qsort(weak, nk2, sizeof *weak, shift_cmp);
                for (size_t j = 0; j < nk2; j++) {
                    SET_VISITED(&x, weak[j].r);
                    co[weak[j].r] = weak[j].o; cP[weak[j].r] = weak[j].P;
                    order[no++] = weak[j].r;
                    g_filled++;
                }
            }
        }
        if (no - chain_start == 1) {
            /* A read with no confident neighbour goes to the end of the
             * run, so it does not scatter into the local blocks. */
            g_singletons++;
            deferred[nd++] = order[--no];
        } else {
            for (size_t i = chain_start; i < no; i++) stack[ns++] = (uint32_t)order[i];
            seg[nseg++] = chain_start;
            g_kind_chains[kind]++;
            if (g_pos && chain_start > 0) {
                int64_t a = g_pos[order[chain_start - 1]], b = g_pos[order[chain_start]];
                if (a > 0 && b > 0) {
                    g_kind_mapped[kind]++;
                    if (llabs(a - b) > 10000) g_kind_far[kind]++;
                }
            }
        }
    }
    memcpy(order + no, deferred, nd * sizeof *order);
    for (size_t i = 0; i < nd; i++) seg[nseg++] = no + i;
    no += nd;
    seg[nseg] = no;
    if (g_scaffold && scaffold_order(&x, order, n, seg, nseg, co, cP)) return -1;
    free(seg);
    free(co); free(cP); free(walked); free(weak);
    free(stack);
    free(deferred);
    free(left);
    free(x.vt);
    free(x.used);
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
    /* How many blocks each 100 kb bin's reads are spread over, weighted by
     * reads: coordinate order keeps a bin in one or two blocks. */
    {
        uint16_t *last = calloc(NBINS, sizeof *last);
        unsigned *nblk = calloc(NBINS, sizeof *nblk), *nrd = calloc(NBINS, sizeof *nrd);
        unsigned long long bu = 0, blk = 1, wsum = 0, wn = 0;
        if (last && nblk && nrd) {
            for (size_t u = 0; u < n_units; u++) {
                for (size_t k = 0; k < unit; k++) {
                    size_t r = order[u] * unit + k;
                    if (pos[r] > 0) {
                        size_t bin = (size_t)(pos[r] / 100000);
                        if (bin >= NBINS) bin = NBINS - 1;
                        if (last[bin] != blk) { last[bin] = (uint16_t)blk; nblk[bin]++; }
                        nrd[bin]++;
                    }
                    bu += len[r];
                }
                if (bu >= block_bases) { blk++; bu = 0; }
            }
            for (size_t b = 0; b < NBINS; b++) { wsum += (unsigned long long)nblk[b] * nrd[b]; wn += nrd[b]; }
            printf("locality: a 100 kb bin's reads span %.2f blocks (read-weighted mean)\n",
                   wn ? (double)wsum / (double)wn : 0.0);
        }
        free(last); free(nblk); free(nrd);
    }
    printf("locality: consecutive mapped reads within 1 kb %.2f%%, mean 100 kb bins per block %.1f, "
           "reads in each block's %d busiest bins %.2f%% (%llu blocks)\n",
           pairs ? 100.0 * (double)near / (double)pairs : 0.0, blocks ? (double)bins_total / (double)blocks : 0.0,
           (int)TOP, mapped ? 100.0 * (double)in_top / (double)mapped : 0.0, blocks);
    free(pos);
}

/* Hash of a FASTQ name line: up to the first whitespace, without a
 * trailing /1 or /2. */
static uint64_t name_hash(const char *line, ssize_t l)
{
    ssize_t e = 1;
    while (e < l && line[e] != ' ' && line[e] != '\t' && line[e] != '\n' && line[e] != '\r') e++;
    if (e >= 3 && line[e - 2] == '/' && (line[e - 1] == '1' || line[e - 1] == '2')) e -= 2;
    uint64_t h = 0xcbf29ce484222325ull;
    for (ssize_t i = 1; i < e; i++) { h ^= (uint8_t)line[i]; h *= 0x100000001b3ull; }
    return h;
}

typedef struct { uint64_t h; uint32_t r; } name_ent;

static int name_cmp(const void *a, const void *b)
{
    const name_ent *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    return x->r < y->r ? -1 : (x->r > y->r);
}

/* Mates from name hashes: two reads sharing a name are mates. */
static uint32_t *pair_by_name(const uint64_t *names, size_t n)
{
    name_ent *e = malloc((n ? n : 1) * sizeof *e);
    uint32_t *mate = malloc((n ? n : 1) * sizeof *mate);
    if (!e || !mate) { free(e); free(mate); return NULL; }
    for (size_t i = 0; i < n; i++) { e[i].h = names[i]; e[i].r = (uint32_t)i; mate[i] = UINT32_MAX; }
    qsort(e, n, sizeof *e, name_cmp);
    for (size_t i = 0; i < n;) {
        size_t j = i + 1;
        while (j < n && e[j].h == e[i].h) j++;
        if (j - i == 2) { mate[e[i].r] = e[i + 1].r; mate[e[i + 1].r] = e[i].r; }
        i = j;
    }
    free(e);
    return mate;
}

/* One mapping position per read for the --pos diagnostics, or NULL. */
static int64_t *load_positions(const char *path, size_t n)
{
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    int64_t *pos = malloc((n ? n : 1) * sizeof *pos);
    size_t m = 0;
    long long v;
    while (pos && m < n && fscanf(f, "%lld", &v) == 1) pos[m++] = v;
    fclose(f);
    if (m != n) { free(pos); return NULL; }
    return pos;
}

/* The minimizer window for --group-w 0: ceil(median read length / 8),
 * clamped to 8..32, so a read carries ~16 minimizers whatever its length
 * (32 for 250 bp reads, 13 for 100 bp). Two of them must agree on an
 * offset (--min-votes), which a read with only ~5 rarely manages. */
static int auto_window(const uint64_t *len, size_t n)
{
    enum { MAXL = 65536 };
    size_t *hist = calloc(MAXL, sizeof *hist);
    if (!hist || !n) { free(hist); return 32; }
    for (size_t i = 0; i < n; i++) hist[len[i] < MAXL ? len[i] : MAXL - 1]++;
    size_t acc = 0, med = 0;
    while (med < MAXL && (acc += hist[med]) < (n + 1) / 2) med++;
    free(hist);
    int w = (int)((med + 7) / 8);
    return w < 8 ? 8 : w > 32 ? 32 : w;
}

int main(int argc, char **argv)
{
    coder c;
    memset(&c, 0, sizeof c);
    ttio_seq_cm_default_params(&c.p);
    unsigned long long block_bases = 64ull << 20;
    int ai = 1, n_orders = 0, group = 0, chain = 0, lib = 0, layout = 0, oracle_sort = 0, as_is = 0, no_code = 0, pairs = 0, group_k = 20, group_w = 0;
    unsigned max_occ = 256;
    const char *pos_path = NULL;  /* diagnostic: one mapping position per read, input order */
    const char *dump_path = NULL; /* diagnostic: write the grouped order, one input index per line */
    while (ai < argc) {
        if (!strcmp(argv[ai], "--block-bases") && ai + 1 < argc) block_bases = strtoull(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "--table-bits") && ai + 1 < argc) c.p.table_bits = (uint8_t)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--no-rc")) c.p.flags = 0;
        else if (!strcmp(argv[ai], "--group")) group = 1;
        else if (!strcmp(argv[ai], "--chain")) group = chain = 1;
        else if (!strcmp(argv[ai], "--lib")) group = lib = 1;
        else if (!strcmp(argv[ai], "--layout")) group = chain = layout = 1;
        else if (!strcmp(argv[ai], "--as-is")) group = as_is = 1;
        else if (!strcmp(argv[ai], "--no-code")) no_code = 1;
        else if (!strcmp(argv[ai], "--pairs")) pairs = 1;
        else if (!strcmp(argv[ai], "--group-k") && ai + 1 < argc) group_k = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--group-w") && ai + 1 < argc) group_w = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--max-occ") && ai + 1 < argc) max_occ = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--min-votes") && ai + 1 < argc) g_min_votes = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--seed-votes") && ai + 1 < argc) g_seed_votes = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--fill")) g_fill = 1;
        else if (!strcmp(argv[ai], "--mates")) g_mates = 1;
        else if (!strcmp(argv[ai], "--oracle-sort")) oracle_sort = 1;
        else if (!strcmp(argv[ai], "--scaffold")) g_scaffold = 1;
        else if (!strcmp(argv[ai], "--scaffold-min") && ai + 1 < argc) g_scaffold_min = (unsigned)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--cluster-cap") && ai + 1 < argc) g_cluster_cap = strtoull(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "--pos") && ai + 1 < argc) pos_path = argv[++ai];
        else if (!strcmp(argv[ai], "--dump-order") && ai + 1 < argc) dump_path = argv[++ai];
        else if (n_orders < 3) c.p.orders[n_orders++] = (uint8_t)atoi(argv[ai]);
        ai++;
    }
    if (n_orders) c.p.n_orders = (uint8_t)n_orders;
    if (group_k < 8 || group_k > 32) { fprintf(stderr, "--group-k must be 8 to 32\n"); return 2; }
    if (group_w < 0) { fprintf(stderr, "--group-w must be 0 (from read length) or more\n"); return 2; }
    if (chain && pairs) { fprintf(stderr, "--chain orders single reads; --pairs is not supported with it\n"); return 2; }
    long long n_chains = 0;
    double t_index = 0;

    size_t cap = 1 << 26, used = 0, rcap = 1 << 20, n = 0;
    uint8_t *seq = malloc(cap);
    uint64_t *len = malloc(rcap * sizeof *len);
    uint64_t *names = g_mates ? malloc(rcap * sizeof *names) : NULL;
    /* --lib: the raw names (no '@', no newline) for ttio_seq_group. */
    size_t rncap = lib ? (size_t)1 << 24 : 0, rnused = 0;
    char *rawnames = lib ? malloc(rncap) : NULL;
    uint64_t *rnoff = lib ? malloc(rcap * sizeof *rnoff) : NULL;
    if (rnoff) rnoff[0] = 0;
    char *line = NULL;
    size_t lcap = 0;
    ssize_t l;
    unsigned long long ln = 0;
    int done = 0;
    double perm_ideal_bits = 0, perm_fixed_bits = 0, t_group = 0;
    while (!done) {
        l = getline(&line, &lcap, stdin);
        if (l > 0 && (ln & 3) == 0 && names) names[n] = name_hash(line, l);
        if (l > 0 && (ln & 3) == 0 && rawnames) {
            ssize_t e = l;
            while (e > 1 && (line[e - 1] == '\n' || line[e - 1] == '\r')) e--;
            while (rnused + (size_t)e > rncap) { rncap *= 2; rawnames = realloc(rawnames, rncap); }
            memcpy(rawnames + rnused, line + 1, (size_t)(e - 1));
            rnused += (size_t)(e - 1);
        }
        if (l > 0 && (ln++ & 3) == 1) {
            while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) l--;
            while (used + (size_t)l > cap) { cap *= 2; seq = realloc(seq, cap); }
            if (n + 1 >= rcap) {
                rcap *= 2;
                len = realloc(len, rcap * sizeof *len);
                if (names) names = realloc(names, rcap * sizeof *names);
                if (rnoff) rnoff = realloc(rnoff, rcap * sizeof *rnoff);
            }
            memcpy(seq + used, line, (size_t)l);
            used += (size_t)l;
            len[n++] = (uint64_t)l;
            if (rnoff) rnoff[n] = rnused;
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
        } else if (lib) {
            /* The library's grouping (native/src/seq_group.c), defaults:
             * the same as --chain --fill --mates. */
            uint32_t *o32 = malloc((n ? n : 1) * sizeof *o32);
            int rc = ttio_seq_group(seq, len, n, (const uint8_t *)rawnames, rnoff, NULL, o32);
            if (rc) { fprintf(stderr, "ttio_seq_group failed: %d\n", rc); return 1; }
            for (size_t u = 0; u < n; u++) order[u] = o32[u];
            free(o32);
            printf("lib: ttio_seq_group defaults\n");
        } else if (chain) {
            if (group_w == 0) group_w = auto_window(len, n);
            if (pos_path) g_pos = load_positions(pos_path, n);
            if (g_mates) {
                uint32_t *mate = pair_by_name(names, n);
                if (!mate) { fprintf(stderr, "out of memory pairing mates\n"); return 1; }
                g_mate = mate;
            }
            n_chains = chain_order(seq, off, len, n, group_k, group_w, max_occ, layout, order, &t_index);
            if (n_chains < 0) { fprintf(stderr, "out of memory while chaining\n"); return 1; }
            printf("chain k %d w %d max-occ %u min-votes %u: %lld chains (%.1f reads each; "
                   "%lld single reads, %lld backtrack seeds, %lld input-order seeds), index %.1f s\n",
                   group_k, group_w, max_occ, g_min_votes, n_chains,
                   (double)n / (double)(n_chains ? n_chains : 1), g_singletons, g_backtrack_seeds, g_fallbacks, t_index);
            if (g_fill) printf("filled %lld weakly linked reads beside their chains\n", g_filled);
            if (g_mates) printf("mate seeds %lld\n", g_mate_seeds);
            if (g_scaffold) printf("scaffold: %lld merges, %lld clusters\n", g_scaffold_joins, g_scaffold_paths);
            if (g_pos)
                printf("chain starts (2+ reads): backtrack %lld, %lld of %lld mapped more than 10 kb away; "
                       "input order %lld, %lld of %lld; mate %lld, %lld of %lld\n",
                       g_kind_chains[0], g_kind_far[0], g_kind_mapped[0],
                       g_kind_chains[1], g_kind_far[1], g_kind_mapped[1],
                       g_kind_chains[2], g_kind_far[2], g_kind_mapped[2]);
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
        if (dump_path) {
            FILE *df = fopen(dump_path, "w");
            if (df) {
                for (size_t u = 0; u < n_units; u++) fprintf(df, "%llu\n", (unsigned long long)order[u]);
                fclose(df);
            }
        }
        if (oracle_sort && g_pos && unit == 1) {
            /* Diagnostic only: same reads per block, each block sorted by
             * true position, to separate block membership from the order
             * within a block. */
            placed *pl = malloc((n ? n : 1) * sizeof *pl);
            size_t bs = 0, bu = 0;
            for (size_t u = 0; u <= n_units; u++) {
                if (u == n_units || bu >= block_bases) {
                    for (size_t i = bs; i < u; i++) { pl[i - bs].c = g_pos[order[i]]; pl[i - bs].r = (uint32_t)order[i]; }
                    qsort(pl, u - bs, sizeof *pl, placed_cmp);
                    for (size_t i = bs; i < u; i++) order[i] = pl[i - bs].r;
                    bs = u;
                    bu = 0;
                }
                if (u < n_units) bu += len[order[u]];
            }
            free(pl);
            printf("oracle: each block sorted by mapping position\n");
        }
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
