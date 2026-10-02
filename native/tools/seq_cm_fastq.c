/*
 * seq_cm_fastq — code a FASTQ file's bases with SEQ_CM (codec id 19) the
 * way the blocks_v1 writer does: the reads are cut into blocks of about
 * --block-bases bases (default 64 MiB, the writer's block size), each block
 * coded on its own, then decoded and compared.
 *
 * usage: seq_cm_fastq [--block-bases N] [--table-bits B] [--no-rc]
 *                     [k1 [k2 [k3]]] < reads.fastq
 */
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

int main(int argc, char **argv)
{
    ttio_seq_cm_params p;
    ttio_seq_cm_default_params(&p);
    unsigned long long block_bases = 64ull << 20;
    int ai = 1, n_orders = 0;
    while (ai < argc) {
        if (!strcmp(argv[ai], "--block-bases") && ai + 1 < argc) block_bases = strtoull(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "--table-bits") && ai + 1 < argc) p.table_bits = (uint8_t)atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--no-rc")) p.flags = 0;
        else if (n_orders < 3) p.orders[n_orders++] = (uint8_t)atoi(argv[ai]);
        ai++;
    }
    if (n_orders) p.n_orders = (uint8_t)n_orders;

    size_t cap = 1 << 26, used = 0, rcap = 1 << 20, n = 0;
    uint8_t *seq = malloc(cap);
    uint64_t *len = malloc(rcap * sizeof *len);
    char *line = NULL;
    size_t lcap = 0;
    ssize_t l;
    unsigned long long ln = 0, total_bases = 0, total_bytes = 0, blocks = 0;
    double t_enc = 0, t_dec = 0;
    int done = 0;
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
        if ((used >= block_bases || done) && n) {
            uint8_t *enc = NULL, *dec = NULL;
            size_t enc_len = 0, dec_len = 0;
            double t0 = now();
            int rc = ttio_seq_cm_encode(seq, len, n, &p, &enc, &enc_len);
            double t1 = now();
            if (!rc) rc = ttio_seq_cm_decode(enc, enc_len, len, n, &dec, &dec_len);
            double t2 = now();
            if (rc || dec_len != used || memcmp(dec, seq, used)) {
                fprintf(stderr, "block %llu: round trip failed (rc=%d)\n", blocks, rc);
                return 1;
            }
            t_enc += t1 - t0;
            t_dec += t2 - t1;
            total_bases += used;
            total_bytes += enc_len;
            blocks++;
            ttio_seq_cm_free(enc);
            ttio_seq_cm_free(dec);
            used = 0;
            n = 0;
        }
    }
    printf("orders");
    for (int i = 0; i < p.n_orders; i++) printf(" %u", p.orders[i]);
    printf("%s table_bits %u block_bases %llu: blocks %llu bases %llu bytes %llu  bits/base %.4f  "
           "encode %.1f MB/s decode %.1f MB/s\n",
           p.flags ? " +rc" : "", p.table_bits, block_bases, blocks, total_bases, total_bytes,
           total_bases ? 8.0 * (double)total_bytes / (double)total_bases : 0.0,
           total_bases / 1e6 / t_enc, total_bases / 1e6 / t_dec);
    return 0;
}
