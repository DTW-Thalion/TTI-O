/*
 * test_seq_group.c — read grouping for unaligned runs (M103).
 *
 * The result is a permutation, the same on every call, and on shuffled
 * error-free reads from both strands of a random genome it puts reads from
 * the same place next to each other. Bad parameters are rejected.
 */
#include "ttio_rans.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
#define CHECK(c, msg) do { if (!(c)) { fprintf(stderr, "FAIL: %s\n", msg); failures++; } } while (0)

static uint64_t rng_state = 7;
static uint64_t rnd(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}

int main(void)
{
    uint32_t dummy = 0;
    CHECK(ttio_seq_group(NULL, NULL, 0, NULL, NULL, NULL, &dummy) == 0, "empty run");
    uint64_t one = 4;
    CHECK(ttio_seq_group((const uint8_t *)"ACGT", &one, 1, NULL, NULL, NULL, &dummy)
          == TTIO_RANS_ERR_PARAM, "MATES without names is rejected");
    ttio_seq_group_params bad;
    ttio_seq_group_default_params(&bad);
    bad.k = 40;
    CHECK(ttio_seq_group((const uint8_t *)"ACGT", &one, 1, NULL, NULL, &bad, &dummy)
          == TTIO_RANS_ERR_PARAM, "k out of range is rejected");

    /* 200 kb genome, 20x of 150 bp reads, half reverse-complemented, shuffled,
     * paired names, a few N bases and zero-length reads. */
    const size_t G = 200000, L = 150, n = G * 20 / L;
    uint8_t *genome = malloc(G);
    for (size_t i = 0; i < G; i++) genome[i] = "ACGT"[rnd() & 3];
    uint64_t *lengths = malloc(n * sizeof *lengths), *pos = malloc(n * sizeof *pos);
    uint8_t *seq = malloc(n * L);
    char *names = malloc(n * 24);
    uint64_t *noff = malloc((n + 1) * sizeof *noff);
    size_t used = 0;
    noff[0] = 0;
    for (size_t i = 0; i < n; i++) {
        pos[i] = (size_t)(rnd() % (G - L));
        lengths[i] = (i % 997 == 5) ? 0 : L;
        int nl = snprintf(names + noff[i], 24, "r%zu/%d", i / 2, (int)(i % 2) + 1);
        noff[i + 1] = noff[i] + (uint64_t)nl;
    }
    /* Each read comes from one strand. */
    for (size_t i = 0; i < n; i++) {
        int rc = (int)(rnd() & 1);
        for (size_t j = 0; j < lengths[i]; j++) {
            uint8_t b = rc ? genome[pos[i] + L - 1 - j] : genome[pos[i] + j];
            if (rc) b = b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
            seq[used + j] = (i % 101 == 0 && j == 70) ? 'N' : b;
        }
        used += lengths[i];
    }

    uint32_t *order = malloc(n * sizeof *order), *again = malloc(n * sizeof *again);
    int rc = ttio_seq_group(seq, lengths, n, (const uint8_t *)names, noff, NULL, order);
    CHECK(rc == 0, "group");
    uint8_t *seen = calloc(n, 1);
    int perm = 1;
    for (size_t j = 0; j < n; j++) {
        if (order[j] >= n || seen[order[j]]) { perm = 0; break; }
        seen[order[j]] = 1;
    }
    CHECK(perm, "order is a permutation");
    CHECK(ttio_seq_group(seq, lengths, n, (const uint8_t *)names, noff, NULL, again) == 0
          && memcmp(order, again, n * sizeof *order) == 0, "deterministic");

    size_t near = 0, pairs = 0;
    for (size_t j = 1; j < n; j++) {
        if (!lengths[order[j]] || !lengths[order[j - 1]]) continue;
        pairs++;
        int64_t d = (int64_t)pos[order[j]] - (int64_t)pos[order[j - 1]];
        if (d < 1000 && d > -1000) near++;
    }
    CHECK(pairs && near * 100 >= pairs * 95, "grouped reads are local (>= 95% within 1 kb)");
    printf("seq_group: %zu reads, %.1f%% of consecutive reads within 1 kb\n",
           n, pairs ? 100.0 * (double)near / (double)pairs : 0.0);

    /* Without mates or fill it still returns a permutation. */
    ttio_seq_group_params plain;
    ttio_seq_group_default_params(&plain);
    plain.flags = 0;
    CHECK(ttio_seq_group(seq, lengths, n, NULL, NULL, &plain, again) == 0, "no mates, no fill");
    memset(seen, 0, n);
    perm = 1;
    for (size_t j = 0; j < n; j++) {
        if (again[j] >= n || seen[again[j]]) { perm = 0; break; }
        seen[again[j]] = 1;
    }
    CHECK(perm, "plain order is a permutation");

    free(genome); free(lengths); free(pos); free(seq); free(names); free(noff);
    free(order); free(again); free(seen);
    if (failures) { fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    printf("test_seq_group: all passed\n");
    return 0;
}
