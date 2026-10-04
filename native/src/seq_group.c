/*
 * seq_group.c — group an unaligned run's reads by sequence before blocking
 * (M103). Returns a permutation that puts reads from the same place in the
 * genome next to each other, so a blocks_v1 block holds them together and
 * SEQ_CM (codec id 19) sees each read's overlap partners. The writer stores
 * the permutation so readers restore input order.
 *
 * Method (M103 Phase 0 README, "Grouping reads before blocking"):
 *  1. Index every read by its (w, k) window minimizers over both strands;
 *     minimizers shared by more than max_occ reads (repeats) are ignored.
 *  2. Chain: from a seed, walk left to the unplaced read that starts the
 *     least distance before the current one, then from the seed right to
 *     the one that starts the least distance after it, tracking each read's
 *     strand and start in the chain's frame. A neighbour counts only when
 *     min_votes of the current read's minimizers agree on its start.
 *  3. FILL: after each chain read, its unplaced neighbours on even one
 *     minimizer that start within a read length of it, by start.
 *  4. Seeds: with MATES, the unplaced mate (same name, /1 or /2 dropped) of
 *     one of the last MATE_LOOKBACK placed reads; then the most recently
 *     placed read that still has an unplaced neighbour (a stack, popped once
 *     a read has none left); then the next read in input order.
 *  5. A read that chains to nothing goes to the end of the run.
 *
 * Every tie is broken by a total order, so the permutation is the same on
 * every platform and in every SDK.
 *
 * Copyright (c) 2026 Thalion Global. All rights reserved.
 */
#include "ttio_rans.h"

#include <stdlib.h>
#include <string.h>

enum { MATE_LOOKBACK = 4096, VOTE_BITS = 17 };

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

/* ── Minimizer index ─────────────────────────────────────────────────── */

typedef struct {
    uint64_t h;
    uint32_t r;
    uint16_t off;   /* k-mer start in the read, read's own orientation */
    uint8_t s;      /* 1: the canonical k-mer is the read's forward k-mer */
    uint8_t pad;
} mz_ent;

typedef struct { uint32_t r; uint16_t off; uint8_t s; uint8_t pad; } mz_hit;
typedef struct { uint32_t b; uint16_t off; uint8_t s; uint8_t pad; } mz_own;

/* Stable LSD radix sort on the hash, 16 bits a pass. Entries arrive in
 * (read, offset) order, which equal hashes keep. */
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
    free(t);        /* four passes: the result is back in the caller's buffer */
    free(cnt);
    return 0;
}

/* Window minimizers of one read into out; returns how many. */
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

/* ── Neighbour search ────────────────────────────────────────────────── */

typedef struct {
    int64_t P;      /* start in the chain's frame */
    uint32_t r;     /* UINT32_MAX: empty slot */
    uint8_t o;
    uint16_t votes;
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
    unsigned min_votes;
    cand *vt;                   /* vote table, open addressing on (r, o, P) */
    uint32_t *used;             /* its occupied slots, in insertion order */
} ctx_t;

#define VISITED(x, r) (((x)->visited[(r) >> 6] >> ((r) & 63)) & 1)
#define SET_VISITED(x, r) ((x)->visited[(r) >> 6] |= 1ull << ((r) & 63))

/* For every unvisited read sharing a minimizer with read cur (orientation
 * o, start P), count how many of cur's minimizers place it at each start
 * and strand; returns the number of distinct placements in x->used. */
static size_t tally(const ctx_t *x, uint32_t cur, int o, int64_t P)
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
            if (h->r == cur || VISITED(x, h->r)) continue;
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

static void clear_tally(const ctx_t *x, size_t nu)
{
    for (size_t i = 0; i < nu; i++) x->vt[x->used[i]].r = UINT32_MAX;
}

enum { NB_LEFT, NB_RIGHT };

/* NB_RIGHT: the unvisited read with the least shift >= 0; NB_LEFT: the one
 * with the greatest shift < 0; at least `votes` agreeing minimizers. Ties go
 * to the lower read, then orientation 0. */
static int best_neighbour(const ctx_t *x, uint32_t cur, int o, int64_t P, int mode, unsigned votes,
                          uint32_t *nr, int *no, int64_t *nP)
{
    size_t nu = tally(x, cur, o, P);
    int found = 0;
    int64_t best = mode == NB_RIGHT ? INT64_MAX : INT64_MIN;
    for (size_t i = 0; i < nu; i++) {
        const cand *cd = &x->vt[x->used[i]];
        if (cd->votes < votes) continue;
        int64_t sh = cd->P - P;
        int tie = sh == best && (cd->r < *nr || (cd->r == *nr && cd->o < *no));
        int better = mode == NB_RIGHT ? sh >= 0 && (sh < best || tie)
                                      : sh < 0 && (sh > best || tie);
        if (better) {
            best = sh;
            *nr = cd->r;
            *no = cd->o;
            *nP = cd->P;
            found = 1;
        }
    }
    clear_tally(x, nu);
    return found;
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

/* ── Mates by name ───────────────────────────────────────────────────── */

typedef struct { uint64_t h; uint32_t r; } name_ent;

static int name_cmp(const void *a, const void *b)
{
    const name_ent *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    return x->r < y->r ? -1 : (x->r > y->r);
}

/* A name up to its first whitespace, without a trailing /1 or /2; two reads
 * whose names hash equal (and no third) are mates. */
static int pair_by_name(const uint8_t *names, const uint64_t *off, size_t n, uint32_t *mate)
{
    name_ent *e = malloc((n ? n : 1) * sizeof *e);
    if (!e) return -1;
    for (size_t i = 0; i < n; i++) {
        const uint8_t *s = names + off[i];
        uint64_t l = off[i + 1] - off[i], end = 0;
        while (end < l && s[end] != ' ' && s[end] != '\t') end++;
        if (end >= 2 && s[end - 2] == '/' && (s[end - 1] == '1' || s[end - 1] == '2')) end -= 2;
        uint64_t h = 0xcbf29ce484222325ull;
        for (uint64_t j = 0; j < end; j++) { h ^= s[j]; h *= 0x100000001b3ull; }
        e[i].h = h;
        e[i].r = (uint32_t)i;
        mate[i] = UINT32_MAX;
    }
    qsort(e, n, sizeof *e, name_cmp);
    for (size_t i = 0; i < n;) {
        size_t j = i + 1;
        while (j < n && e[j].h == e[i].h) j++;
        if (j - i == 2) { mate[e[i].r] = e[i + 1].r; mate[e[i + 1].r] = e[i].r; }
        i = j;
    }
    free(e);
    return 0;
}

/* ── Public API ──────────────────────────────────────────────────────── */

void ttio_seq_group_default_params(ttio_seq_group_params *p)
{
    if (!p) return;
    p->k = 20;
    p->window = 0;
    p->min_votes = 2;
    p->flags = TTIO_SEQ_GROUP_FLAG_FILL | TTIO_SEQ_GROUP_FLAG_MATES;
    p->max_occ = 256;
}

/* ceil(median read length / 8), clamped to 8..32: ~16 minimizers a read. */
static int auto_window(const uint64_t *len, size_t n)
{
    enum { MAXL = 65536 };
    size_t *hist = calloc(MAXL, sizeof *hist);
    if (!hist) return -1;
    for (size_t i = 0; i < n; i++) hist[len[i] < MAXL ? len[i] : MAXL - 1]++;
    size_t acc = 0, med = 0;
    while (med < MAXL && (acc += hist[med]) < (n + 1) / 2) med++;
    free(hist);
    int w = (int)((med + 7) / 8);
    return w < 8 ? 8 : w > 32 ? 32 : w;
}

int ttio_seq_group(const uint8_t *seq, const uint64_t *lengths, uint64_t n_reads,
                   const uint8_t *names, const uint64_t *name_offsets,
                   const ttio_seq_group_params *params, uint32_t *order)
{
    if (!order || (n_reads && !lengths) || n_reads >= UINT32_MAX) return TTIO_RANS_ERR_PARAM;
    ttio_seq_group_params p;
    if (params) p = *params; else ttio_seq_group_default_params(&p);
    if (p.k < 8 || p.k > 32 || p.min_votes < 1 || p.max_occ < 2
        || (p.flags & ~(uint8_t)(TTIO_SEQ_GROUP_FLAG_FILL | TTIO_SEQ_GROUP_FLAG_MATES)))
        return TTIO_RANS_ERR_PARAM;
    if ((p.flags & TTIO_SEQ_GROUP_FLAG_MATES) && n_reads && (!names || !name_offsets))
        return TTIO_RANS_ERR_PARAM;
    const size_t n = (size_t)n_reads;
    if (n == 0) return 0;
    uint64_t total = 0;
    for (size_t r = 0; r < n; r++) {
        if (lengths[r] > UINT64_MAX - total) return TTIO_RANS_ERR_PARAM;
        total += lengths[r];
    }
    if (total && !seq) return TTIO_RANS_ERR_PARAM;
    int w = p.window ? p.window : auto_window(lengths, n);
    if (w < 0) return TTIO_RANS_ERR_ALLOC;
    const int k = p.k;

    int rc = TTIO_RANS_ERR_ALLOC;
    uint64_t *off = NULL, *ostart = NULL, *fill = NULL;
    mz_ent *ents = NULL, *tmp = NULL;
    uint64_t *hb = NULL;
    uint8_t *sb = NULL;
    uint32_t *dq = NULL, *bstart = NULL, *mate = NULL;
    mz_hit *hits = NULL;
    mz_own *own = NULL;
    ctx_t x = {0};
    uint32_t *left = NULL, *deferred = NULL, *stack = NULL, *walked = NULL;
    uint8_t *co = NULL;
    int64_t *cP = NULL;
    cand *weak = NULL;

    off = malloc((n + 1) * sizeof *off);
    if (!off) goto done;
    off[0] = 0;
    for (size_t r = 0; r < n; r++) off[r + 1] = off[r] + lengths[r];

    /* 1. Minimizer index. */
    uint64_t maxl = 0;
    for (size_t r = 0; r < n; r++) if (lengths[r] > maxl) maxl = lengths[r];
    if (maxl > 65535) maxl = 65535;
    size_t ecap = n * 8 + 1024, ne = 0;
    ents = malloc(ecap * sizeof *ents);
    hb = malloc((maxl + 1) * sizeof *hb);
    sb = malloc(maxl + 1);
    dq = malloc((maxl + 1) * sizeof *dq);
    tmp = malloc((maxl + 1) * sizeof *tmp);
    if (!ents || !hb || !sb || !dq || !tmp) goto done;
    for (size_t r = 0; r < n; r++) {
        size_t m = read_minimizers(seq + off[r], lengths[r], k, w, (uint32_t)r, tmp, hb, sb, dq);
        if (ne + m > ecap) {
            size_t ncap = (ecap + m) * 3 / 2;
            mz_ent *g = realloc(ents, ncap * sizeof *ents);
            if (!g) goto done;
            ents = g;
            ecap = ncap;
        }
        memcpy(ents + ne, tmp, m * sizeof *tmp);
        ne += m;
    }
    free(hb); hb = NULL;
    free(sb); sb = NULL;
    free(dq); dq = NULL;
    free(tmp); tmp = NULL;
    if (ne >= UINT32_MAX) { rc = TTIO_RANS_ERR_PARAM; goto done; }
    if (sort_ents(ents, ne)) goto done;

    size_t nb = 0;
    for (size_t i = 0; i < ne; i++) if (i == 0 || ents[i].h != ents[i - 1].h) nb++;
    bstart = malloc((nb + 1) * sizeof *bstart);
    hits = malloc((ne ? ne : 1) * sizeof *hits);
    ostart = calloc(n + 1, sizeof *ostart);
    own = malloc((ne ? ne : 1) * sizeof *own);
    fill = malloc((n + 1) * sizeof *fill);
    if (!bstart || !hits || !ostart || !own || !fill) goto done;
    for (size_t i = 0; i < ne; i++) ostart[ents[i].r + 1]++;
    for (size_t r = 0; r < n; r++) ostart[r + 1] += ostart[r];
    memcpy(fill, ostart, (n + 1) * sizeof *fill);
    size_t b = (size_t)-1;
    for (size_t i = 0; i < ne; i++) {
        if (i == 0 || ents[i].h != ents[i - 1].h) bstart[++b] = (uint32_t)i;
        hits[i] = (mz_hit){ ents[i].r, ents[i].off, ents[i].s, 0 };
        own[fill[ents[i].r]++] = (mz_own){ (uint32_t)b, ents[i].off, ents[i].s, 0 };
    }
    bstart[nb] = (uint32_t)ne;
    free(fill); fill = NULL;
    free(ents); ents = NULL;

    if (p.flags & TTIO_SEQ_GROUP_FLAG_MATES) {
        mate = malloc(n * sizeof *mate);
        if (!mate || pair_by_name(names, name_offsets, n, mate)) goto done;
    }

    /* 2-5. Chains. */
    x.len = lengths;
    x.bstart = bstart;
    x.hits = hits;
    x.ostart = ostart;
    x.own = own;
    x.visited = calloc((n + 63) / 64, sizeof(uint64_t));
    x.k = k;
    x.max_occ = p.max_occ;
    x.min_votes = p.min_votes;
    x.vt = malloc(((size_t)1 << VOTE_BITS) * sizeof(cand));
    x.used = malloc(((size_t)1 << (VOTE_BITS - 1)) * sizeof(uint32_t));
    left = malloc(n * sizeof *left);
    deferred = malloc(n * sizeof *deferred);
    stack = malloc(n * sizeof *stack);
    walked = malloc(n * sizeof *walked);
    co = malloc(n);
    cP = malloc(n * sizeof *cP);
    weak = malloc(((size_t)1 << (VOTE_BITS - 1)) * sizeof *weak);
    if (!x.visited || !x.vt || !x.used || !left || !deferred || !stack || !walked || !co || !cP || !weak)
        goto done;
    for (size_t i = 0; i < ((size_t)1 << VOTE_BITS); i++) x.vt[i].r = UINT32_MAX;

    size_t no = 0, nd = 0, ns = 0, next_input = 0;
    for (;;) {
        size_t seed = (size_t)-1;
        if (mate)
            for (size_t back = 0; back < MATE_LOOKBACK && back < no && seed == (size_t)-1; back++) {
                uint32_t m = mate[order[no - 1 - back]];
                if (m != UINT32_MAX && !VISITED(&x, m)) seed = m;
            }
        while (ns && seed == (size_t)-1) {
            uint32_t r = stack[ns - 1], nr2 = 0;
            int on2 = 1;
            int64_t nP2 = 0;
            if (best_neighbour(&x, r, 1, 0, NB_RIGHT, 1, &nr2, &on2, &nP2) ||
                best_neighbour(&x, r, 1, 0, NB_LEFT, 1, &nr2, &on2, &nP2))
                seed = nr2;
            else
                ns--;
        }
        if (seed == (size_t)-1) {
            while (next_input < n && VISITED(&x, next_input)) next_input++;
            if (next_input == n) break;
            seed = next_input;
        }
        SET_VISITED(&x, seed);
        co[seed] = 1;
        cP[seed] = 0;
        size_t nl = 0;
        uint32_t cur = (uint32_t)seed, nr = 0;
        int o = 1, on = 1;
        int64_t P = 0, nP = 0;
        while (best_neighbour(&x, cur, o, P, NB_LEFT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            left[nl++] = nr;
            co[nr] = (uint8_t)on; cP[nr] = nP;
            cur = nr; o = on; P = nP;
        }
        size_t chain_start = no;
        while (nl) order[no++] = left[--nl];
        order[no++] = (uint32_t)seed;
        cur = (uint32_t)seed; o = 1; P = 0;
        while (best_neighbour(&x, cur, o, P, NB_RIGHT, x.min_votes, &nr, &on, &nP)) {
            SET_VISITED(&x, nr);
            order[no++] = nr;
            co[nr] = (uint8_t)on; cP[nr] = nP;
            cur = nr; o = on; P = nP;
        }
        if (p.flags & TTIO_SEQ_GROUP_FLAG_FILL) {
            size_t nw = no - chain_start;
            memcpy(walked, order + chain_start, nw * sizeof *walked);
            no = chain_start;
            for (size_t i = 0; i < nw; i++) {
                uint32_t r = walked[i];
                order[no++] = r;
                size_t nu = tally(&x, r, co[r], cP[r]), nk = 0;
                for (size_t j = 0; j < nu; j++) {
                    const cand *cd = &x.vt[x.used[j]];
                    int64_t d = cd->P - cP[r];
                    if (d <= (int64_t)lengths[r] && -d <= (int64_t)lengths[r]) weak[nk++] = *cd;
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
                }
            }
        }
        if (no - chain_start == 1) {
            deferred[nd++] = order[--no];
        } else {
            for (size_t i = chain_start; i < no; i++) stack[ns++] = order[i];
        }
    }
    memcpy(order + no, deferred, nd * sizeof *order);
    rc = (no + nd == n) ? 0 : TTIO_RANS_ERR_PARAM;

done:
    free(off); free(ents); free(tmp); free(hb); free(sb); free(dq);
    free(bstart); free(hits); free(ostart); free(own); free(fill); free(mate);
    free(x.visited); free(x.vt); free(x.used);
    free(left); free(deferred); free(stack); free(walked); free(co); free(cP); free(weak);
    return rc;
}
