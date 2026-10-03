/* Round-trip, parameter and corruption tests for the SEQ_CM codec
 * (codec id 19, M103). Spec: docs/codecs/seq_cm.md. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ttio_rans.h"
#include "seq_cm.h"

/* Always evaluated, unlike assert under NDEBUG. */
#define CHECK(cond)                                                    \
    do {                                                               \
        if (!(cond)) {                                                 \
            fprintf(stderr, "%s:%d: check failed: %s\n",               \
                    __FILE__, __LINE__, #cond);                        \
            exit(1);                                                   \
        }                                                              \
    } while (0)

static uint64_t rng_state = 0x243F6A8885A308D3ull;
static uint32_t rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (uint32_t)rng_state;
}

static int roundtrip(const uint8_t *seq, const uint64_t *len, uint64_t n,
                     const ttio_seq_cm_params *p, size_t *enc_size)
{
    uint8_t *enc = NULL, *dec = NULL;
    size_t enc_len = 0, dec_len = 0;
    int rc = ttio_seq_cm_encode(seq, len, n, p, &enc, &enc_len);
    if (rc) { fprintf(stderr, "encode rc=%d\n", rc); return 0; }
    rc = ttio_seq_cm_decode(enc, enc_len, len, n, &dec, &dec_len);
    if (rc) { fprintf(stderr, "decode rc=%d\n", rc); ttio_seq_cm_free(enc); return 0; }
    uint64_t total = 0;
    for (uint64_t i = 0; i < n; i++) total += len[i];
    int ok = dec_len == total && (total == 0 || memcmp(dec, seq, (size_t)total) == 0);
    if (!ok) fprintf(stderr, "mismatch (len %zu vs %llu)\n", dec_len, (unsigned long long)total);
    /* Encoding is deterministic. */
    uint8_t *enc2 = NULL;
    size_t enc2_len = 0;
    if (ttio_seq_cm_encode(seq, len, n, p, &enc2, &enc2_len) || enc2_len != enc_len
        || memcmp(enc, enc2, enc_len)) { fprintf(stderr, "non-deterministic\n"); ok = 0; }
    if (enc_size) *enc_size = enc_len;
    ttio_seq_cm_free(enc);
    ttio_seq_cm_free(enc2);
    ttio_seq_cm_free(dec);
    return ok;
}

/* Reads sampled from a random genome, from both strands. */
static uint8_t *make_reads(uint64_t genome_len, uint64_t n, uint64_t rlen,
                           uint64_t *len, double err)
{
    uint8_t *g = (uint8_t *)malloc(genome_len);
    for (uint64_t i = 0; i < genome_len; i++) g[i] = "ACGT"[rnd() & 3];
    uint8_t *s = (uint8_t *)malloc(n * rlen);
    for (uint64_t r = 0; r < n; r++) {
        uint64_t at = rnd() % (genome_len - rlen);
        int rev = rnd() & 1;
        for (uint64_t i = 0; i < rlen; i++) {
            uint8_t b = rev ? g[at + rlen - 1 - i] : g[at + i];
            if (rev) b = b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
            if ((rnd() % 1000000) < err * 1e6) b = "ACGT"[rnd() & 3];
            s[r * rlen + i] = b;
        }
        len[r] = rlen;
    }
    free(g);
    return s;
}

static void test_logistic(void)
{
    int prev = -1;
    for (int d = -2047; d <= 2047; d++) {
        int p = sqc_squash(d);
        CHECK(p >= prev && p >= 0 && p <= 4095);
        prev = p;
    }
    CHECK(sqc_squash(0) == 2047);
    for (int p = 1; p < 4095; p++) {
        int d = sqc_stretch(p);
        CHECK(sqc_squash(d) >= p);
        CHECK(d == -2047 || sqc_squash(d - 1) < p);
    }
    CHECK(sqc_auto_table_bits(0) == 16);
    CHECK(sqc_auto_table_bits(1000) == 16);
    CHECK(sqc_auto_table_bits(1u << 20) == 21);
    CHECK(sqc_auto_table_bits((uint64_t)1 << 40) == 24);
    printf("  logistic and auto table bits ok\n");
}

static void test_edges(void)
{
    ttio_seq_cm_params p;
    ttio_seq_cm_default_params(&p);
    /* No reads. */
    CHECK(roundtrip(NULL, NULL, 0, &p, NULL));
    /* Zero-length reads among others. */
    {
        const char *s = "ACGTNNacgtRYACGT";
        uint64_t len[5] = { 0, 4, 0, 12, 0 };
        CHECK(roundtrip((const uint8_t *)s, len, 5, &p, NULL));
    }
    /* All exceptions, runs spanning read boundaries, a NUL byte. */
    {
        uint8_t s[40];
        memset(s, 'N', sizeof s);
        s[10] = 0;
        s[11] = 'n';
        s[25] = 'A';
        uint64_t len[4] = { 7, 13, 10, 10 };
        CHECK(roundtrip(s, len, 4, &p, NULL));
    }
    /* One base. */
    {
        uint64_t len[1] = { 1 };
        CHECK(roundtrip((const uint8_t *)"G", len, 1, &p, NULL));
    }
    printf("  edge cases ok\n");
}

static void test_params(void)
{
    uint64_t n = 3000, rlen = 150;
    uint64_t *len = (uint64_t *)malloc(n * sizeof *len);
    uint8_t *s = make_reads(20000, n, rlen, len, 0.002);
    for (int i = 0; i < 37; i++) s[(rnd() % (n * rlen))] = 'N';

    struct { uint8_t no, o0, o1, o2, tb, fl; } cases[] = {
        { 1, 12, 0, 0, 0, 0 },
        { 1, 16, 0, 0, 12, TTIO_SEQ_CM_FLAG_RC },   /* hashed, tiny table */
        { 2, 12, 20, 0, 18, TTIO_SEQ_CM_FLAG_RC },
        { 3, 11, 16, 24, 0, TTIO_SEQ_CM_FLAG_RC },
        { 3, 11, 16, 24, 10, 0 },
        { 3, 1, 2, 31, 14, TTIO_SEQ_CM_FLAG_RC },
    };
    size_t sz_rc = 0, sz_norc = 0;
    for (size_t c = 0; c < sizeof cases / sizeof cases[0]; c++) {
        ttio_seq_cm_params p;
        ttio_seq_cm_default_params(&p);
        p.n_orders = cases[c].no;
        p.orders[0] = cases[c].o0;
        p.orders[1] = cases[c].o1;
        p.orders[2] = cases[c].o2;
        p.table_bits = cases[c].tb;
        p.flags = cases[c].fl;
        size_t sz = 0;
        CHECK(roundtrip(s, len, n, &p, &sz));
        printf("    orders %u/%u/%u bits %u flags %u: %zu bytes (%.3f bits/base)\n",
               cases[c].o0, cases[c].o1, cases[c].o2, cases[c].tb, cases[c].fl, sz,
               8.0 * (double)sz / (double)(n * rlen));
        if (c == 3) sz_rc = sz;
    }
    /* Reverse-complement training must help on two-strand reads. */
    {
        ttio_seq_cm_params p;
        ttio_seq_cm_default_params(&p);
        p.flags = 0;
        CHECK(roundtrip(s, len, n, &p, &sz_norc));
        CHECK(sz_rc < sz_norc);
        /* 20 kb genome at ~22x coverage: far below 2 bits per base. */
        CHECK(8.0 * (double)sz_rc / (double)(n * rlen) < 0.5);
    }
    free(s);
    free(len);
    printf("  parameter sweep ok\n");
}

static void test_random_bases(void)
{
    /* Incompressible bases cost about 2 bits each, not much more. */
    uint64_t n = 2000, rlen = 100;
    uint64_t *len = (uint64_t *)malloc(n * sizeof *len);
    uint8_t *s = (uint8_t *)malloc(n * rlen);
    for (uint64_t i = 0; i < n * rlen; i++) s[i] = "ACGT"[rnd() & 3];
    for (uint64_t i = 0; i < n; i++) len[i] = rlen;
    size_t sz = 0;
    CHECK(roundtrip(s, len, n, NULL, &sz));
    double bpb = 8.0 * (double)sz / (double)(n * rlen);
    printf("  random bases: %.3f bits/base\n", bpb);
    CHECK(bpb < 2.1);
    free(s);
    free(len);
}

static void test_bad_input(void)
{
    const char *s = "ACGTACGTNNACGT";
    uint64_t len[2] = { 6, 8 };
    uint8_t *enc = NULL, *dec = NULL;
    size_t enc_len = 0, dec_len = 0;
    ttio_seq_cm_params p;
    ttio_seq_cm_default_params(&p);

    p.n_orders = 0;
    CHECK(ttio_seq_cm_encode((const uint8_t *)s, len, 2, &p, &enc, &enc_len) == TTIO_RANS_ERR_PARAM);
    ttio_seq_cm_default_params(&p);
    p.orders[1] = p.orders[0];
    CHECK(ttio_seq_cm_encode((const uint8_t *)s, len, 2, &p, &enc, &enc_len) == TTIO_RANS_ERR_PARAM);
    ttio_seq_cm_default_params(&p);
    p.table_bits = 31;
    CHECK(ttio_seq_cm_encode((const uint8_t *)s, len, 2, &p, &enc, &enc_len) == TTIO_RANS_ERR_PARAM);
    ttio_seq_cm_default_params(&p);
    p.flags = 0x80;
    CHECK(ttio_seq_cm_encode((const uint8_t *)s, len, 2, &p, &enc, &enc_len) == TTIO_RANS_ERR_PARAM);

    CHECK(ttio_seq_cm_encode((const uint8_t *)s, len, 2, NULL, &enc, &enc_len) == 0);
    /* Lengths that disagree with the blob. */
    uint64_t bad_len[2] = { 7, 8 };
    CHECK(ttio_seq_cm_decode(enc, enc_len, bad_len, 2, &dec, &dec_len) == TTIO_RANS_ERR_PARAM);
    CHECK(ttio_seq_cm_decode(enc, enc_len, len, 1, &dec, &dec_len) == TTIO_RANS_ERR_PARAM);
    /* Header damage. */
    uint8_t *cp = (uint8_t *)malloc(enc_len);
    memcpy(cp, enc, enc_len);
    cp[0] = 'X';
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    memcpy(cp, enc, enc_len);
    cp[4] = 2;
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    memcpy(cp, enc, enc_len);
    cp[10] = 40;
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    memcpy(cp, enc, enc_len);
    cp[40] = 0xFF;          /* exception section longer than the blob */
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    memcpy(cp, enc, enc_len);
    cp[32] = 2;             /* claims a second exception run */
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    memcpy(cp, enc, enc_len);
    cp[SQC_HEADER_LEN + 2] = 'C';   /* exception byte that is a base */
    CHECK(ttio_seq_cm_decode(cp, enc_len, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    CHECK(ttio_seq_cm_decode(enc, SQC_HEADER_LEN - 1, len, 2, &dec, &dec_len) == TTIO_RANS_ERR_CORRUPT);
    /* A truncated coder stream decodes without overrunning (bytes past
     * the end read as zero); the bases may differ, the call must not crash. */
    if (ttio_seq_cm_decode(enc, enc_len - 3, len, 2, &dec, &dec_len) == 0) ttio_seq_cm_free(dec);
    free(cp);
    ttio_seq_cm_free(enc);
    printf("  bad input ok\n");
}

int main(void)
{
    printf("test_seq_cm\n");
    test_logistic();
    test_edges();
    test_params();
    test_random_bases();
    test_bad_input();
    printf("all seq_cm tests passed\n");
    return 0;
}
