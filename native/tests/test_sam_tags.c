/* Round-trip, derivation and corruption tests for the SAM_TAGS codec
 * (codec id 18, M101). Spec: docs/codecs/sam_tags.md. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ttio_rans.h"
#include "sam_tags.h"

/* One test run: up to 64 reads against one reference. */
typedef struct {
    uint64_t n;
    const char *seq[64], *cigar[64], *tags[64];
    int64_t pos[64];
    uint16_t chrom[64];
} run_t;

static const char *REF =
    "ACGTACGTAC" "GTTTGGCCAA" "NNNNACGTAC" "GTACGTACGT" "ACGTTTTTTT";

typedef struct {
    uint8_t *seqs, *cigs, *tags;
    uint64_t *so, *co, *to;
    ttio_sam_tags_ctx ctx;
    const uint8_t *refs[1];
    uint64_t ref_len[1];
} built_t;

static uint8_t *concat(const char *const *v, uint64_t n, uint64_t *off) {
    size_t total = 0;
    for (uint64_t i = 0; i < n; i++) total += strlen(v[i]);
    uint8_t *b = (uint8_t *)malloc(total + 1);
    off[0] = 0;
    for (uint64_t i = 0; i < n; i++) {
        size_t l = strlen(v[i]);
        memcpy(b + off[i], v[i], l);
        off[i + 1] = off[i] + l;
    }
    return b;
}

static void build(const run_t *r, int with_ref, built_t *b) {
    memset(b, 0, sizeof *b);
    b->so = calloc(r->n + 1, sizeof(uint64_t));
    b->co = calloc(r->n + 1, sizeof(uint64_t));
    b->to = calloc(r->n + 1, sizeof(uint64_t));
    b->seqs = concat(r->seq, r->n, b->so);
    b->cigs = concat(r->cigar, r->n, b->co);
    b->tags = concat(r->tags, r->n, b->to);
    b->refs[0] = (const uint8_t *)REF;
    b->ref_len[0] = strlen(REF);
    b->ctx.n_reads = r->n;
    b->ctx.sequences = b->seqs;
    b->ctx.seq_offsets = b->so;
    b->ctx.cigars = b->cigs;
    b->ctx.cigar_offsets = b->co;
    b->ctx.positions = r->pos;
    b->ctx.chrom_ids = r->chrom;
    b->ctx.refs = b->refs;
    b->ctx.ref_lengths = b->ref_len;
    b->ctx.n_refs = with_ref ? 1 : 0;
}

static void unbuild(built_t *b) {
    free(b->seqs); free(b->cigs); free(b->tags);
    free(b->so); free(b->co); free(b->to);
}

/* Encodes, checks determinism and the round trip; returns the blob size. */
static size_t round_trip(const run_t *r, int with_ref) {
    built_t b;
    build(r, with_ref, &b);
    uint8_t *enc = NULL, *enc2 = NULL;
    size_t el = 0, el2 = 0;
    int rc = ttio_sam_tags_encode(&b.ctx, b.tags, b.to, &enc, &el);
    assert(rc == 0);
    assert(memcmp(enc, "STG1", 4) == 0);
    rc = ttio_sam_tags_encode(&b.ctx, b.tags, b.to, &enc2, &el2);
    assert(rc == 0 && el2 == el && memcmp(enc, enc2, el) == 0);

    uint8_t *dec = NULL;
    uint64_t *doff = calloc(r->n + 1, sizeof(uint64_t));
    rc = ttio_sam_tags_decode(&b.ctx, enc, el, &dec, doff);
    assert(rc == 0);
    for (uint64_t i = 0; i < r->n; i++) {
        size_t l = (size_t)(doff[i + 1] - doff[i]);
        if (l != strlen(r->tags[i]) || memcmp(dec + doff[i], r->tags[i], l) != 0) {
            fprintf(stderr, "read %llu: want '%s' got '%.*s'\n", (unsigned long long)i,
                    r->tags[i], (int)l, (const char *)(dec + doff[i]));
            assert(0);
        }
    }

    /* Every strict prefix is rejected, and so is a trailing byte. */
    for (size_t cut = 0; cut < el; cut++) {
        uint8_t *d2 = NULL;
        int rc2 = ttio_sam_tags_decode(&b.ctx, enc, cut, &d2, doff);
        assert(rc2 == TTIO_RANS_ERR_CORRUPT);
    }
    uint8_t *longer = malloc(el + 1);
    memcpy(longer, enc, el);
    longer[el] = 0;
    uint8_t *d3 = NULL;
    assert(ttio_sam_tags_decode(&b.ctx, longer, el + 1, &d3, doff) == TTIO_RANS_ERR_CORRUPT);
    free(longer);

    ttio_sam_tags_free(enc);
    ttio_sam_tags_free(enc2);
    ttio_sam_tags_free(dec);
    free(doff);
    unbuild(&b);
    return el;
}

static void expect_md(const char *seq, const char *cigar, int64_t pos, const char *md, int64_t nm) {
    run_t r = {.n = 1, .seq = {seq}, .cigar = {cigar}, .tags = {""}, .pos = {pos}, .chrom = {0}};
    built_t b;
    build(&r, 1, &b);
    char *got = NULL;
    size_t gl = 0;
    int64_t gnm = -1;
    int ok = stg_calc_md_nm(&b.ctx, 0, &got, &gl, &gnm);
    if (md == NULL) {
        assert(ok == 0);
    } else {
        if (ok != 1 || strcmp(got, md) != 0 || gnm != nm) {
            fprintf(stderr, "%s %s @%lld: want %s/%lld got %s/%lld\n", seq, cigar, (long long)pos, md,
                    (long long)nm, ok == 1 ? got : "(none)", (long long)gnm);
            assert(0);
        }
        assert(gl == strlen(md));
    }
    free(got);
    unbuild(&b);
}

static void test_md_nm(void) {
    /* REF[0..9] = ACGTACGTAC */
    expect_md("ACGTACGTAC", "10M", 1, "10", 0);
    expect_md("ACCTACGTAC", "10M", 1, "2G7", 1);
    expect_md("TCGTACGTAA", "10M", 1, "0A8C0", 2);
    expect_md("acgtacgtac", "10M", 1, "10", 0);           /* case-insensitive */
    expect_md("AC=TACGTAC", "10M", 1, "10", 0);           /* '=' matches */
    expect_md("ACNTACGTAC", "10M", 1, "2G7", 1);          /* read N mismatches */
    expect_md("ACGGTAC", "3M3D4M", 1, "3^TAC4", 3);       /* ref ACG [TAC] GTAC */
    expect_md("ACGGTCC", "3M3D4M", 1, "3^TAC2A1", 4);
    expect_md("ACGTTTACG", "4M2I3M", 1, "7", 2);
    expect_md("GGACGTA", "2S5M", 1, "5", 0);              /* soft clip */
    expect_md("ACCGT", "2M3N3M", 1, "5", 0);              /* ref AC [GTA] CGT */
    expect_md("ACGTACGTAC", "10M", 21, "0N0N0N0N6", 4);   /* REF[20..23] = NNNN */
    expect_md("ACGTACGTAC", "1M1D", 1, NULL, 0);          /* query length 1 != 10 */
    expect_md("A", "1M1D", 1, "1^C0", 1);
    expect_md("ACGTACGTAC", "5H10M", 1, "10", 0);
    expect_md("ACGT", "4M", 0, NULL, 0);                  /* unmapped */
    expect_md("ACGT", "*", 1, NULL, 0);
    expect_md("*", "1M", 1, NULL, 0);
    expect_md("ACGT", "5M", 1, NULL, 0);                  /* query length != SEQ */
    expect_md("ACGTACGTAC", "10M", 45, NULL, 0);          /* runs off the reference */
    expect_md("ACGT", "4Q", 1, NULL, 0);
    printf("test_md_nm: PASS\n");
}

static void test_canonical_int(void) {
    int64_t v;
    assert(stg_parse_canonical_int((const uint8_t *)"0", 1, &v) && v == 0);
    assert(stg_parse_canonical_int((const uint8_t *)"-17", 3, &v) && v == -17);
    assert(stg_parse_canonical_int((const uint8_t *)"9223372036854775807", 19, &v) && v == INT64_MAX);
    assert(stg_parse_canonical_int((const uint8_t *)"-9223372036854775808", 20, &v) && v == INT64_MIN);
    assert(!stg_parse_canonical_int((const uint8_t *)"9223372036854775808", 19, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"-0", 2, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"007", 3, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"+5", 2, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"", 0, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"-", 1, &v));
    assert(!stg_parse_canonical_int((const uint8_t *)"1e3", 3, &v));
    printf("test_canonical_int: PASS\n");
}

static void test_empty_run(void) {
    run_t r = {.n = 0};
    round_trip(&r, 0);
    round_trip(&r, 1);
    printf("test_empty_run: PASS\n");
}

static void test_mixed(void) {
    run_t r = {.n = 12};
    const char *tags[12] = {
        "NM:i:0\tMD:Z:10\tAS:i:30\tUQ:i:30\tRG:Z:rg1",
        "NM:i:1\tMD:Z:2G7\tAS:i:47\tUQ:i:47\tRG:Z:rg1",
        "",                                              /* no tags */
        "MD:Z:9\tNM:i:0\tRG:Z:rg2",                      /* MD differs: stored */
        "XA:Z:chr1,+100,10M,0;\tAS:i:007\tXS:i:-0",      /* non-canonical ints */
        "BC:B:C,1,2,3\tXF:f:0.5\tXH:H:1AE301\tXC:A:c",
        "NM:i:3\tMD:Z:3^TAC4\tAS:i:9223372036854775807\tXN:i:-9223372036854775808",
        "bad\tNM:i:0",                                    /* malformed: verbatim */
        "RG:Z:rg1\t",                                     /* trailing tab: verbatim */
        "ZZ:Z:",                                          /* empty Z value */
        "NM:i:0\tMD:Z:10\tAS:i:0\tXS:i:0\tYS:i:0",        /* NM derived; DUP of NM */
        "MD:Z:10\tNM:i:0",
    };
    const char *seqs[12] = {"ACGTACGTAC", "ACCTACGTAC", "ACGTACGTAC", "ACGTACGTAC", "ACGT",
                            "ACGT", "ACGGTAC", "ACGT", "ACGT", "ACGT", "ACGTACGTAC", "ACGTACGTAC"};
    const char *cigs[12] = {"10M", "10M", "10M", "10M", "4M", "4M", "3M3D4M", "4M", "4M", "4M",
                            "10M", "10M"};
    for (int i = 0; i < 12; i++) {
        r.tags[i] = tags[i];
        r.seq[i] = seqs[i];
        r.cigar[i] = cigs[i];
        r.pos[i] = 1;
        r.chrom[i] = 0;
    }
    size_t with = round_trip(&r, 1);
    size_t without = round_trip(&r, 0);
    assert(with < without); /* MD/NM derivation saves bytes */
    printf("test_mixed: PASS (%zu bytes with reference, %zu without)\n", with, without);
}

static void test_columns(void) {
    /* 64 reads: constant text, varying ints across widths, token-rich text. */
    static char tagbuf[64][128];
    run_t r = {.n = 64};
    for (int i = 0; i < 64; i++) {
        snprintf(tagbuf[i], sizeof tagbuf[i],
                 "PG:Z:novoalign\tAS:i:%d\tPQ:i:%d\tXL:i:%lld\tCB:Z:AAACCTG%04dTCAG-1",
                 (i * 30) % 255, i * 1000, (long long)i * 5000000000LL, i % 7);
        r.tags[i] = tagbuf[i];
        r.seq[i] = "ACGT";
        r.cigar[i] = "4M";
        r.pos[i] = 1;
        r.chrom[i] = 0xFFFF;
    }
    round_trip(&r, 1);
    printf("test_columns: PASS\n");
}

static void test_bad_input(void) {
    uint64_t off[2] = {0, 3};
    const uint8_t t[3] = {'A', 0, 'B'};
    ttio_sam_tags_ctx ctx = {.n_reads = 1};
    uint8_t *out = NULL;
    size_t ol = 0;
    assert(ttio_sam_tags_encode(&ctx, t, off, &out, &ol) == TTIO_RANS_ERR_PARAM);
    assert(ttio_sam_tags_encode(NULL, t, off, &out, &ol) == TTIO_RANS_ERR_PARAM);

    /* A DERIVED entry decoded without the reference is corrupt. */
    run_t r = {.n = 1, .seq = {"ACGTACGTAC"}, .cigar = {"10M"}, .tags = {"MD:Z:10"}, .pos = {1},
               .chrom = {0}};
    built_t b;
    build(&r, 1, &b);
    assert(ttio_sam_tags_encode(&b.ctx, b.tags, b.to, &out, &ol) == 0);
    b.ctx.n_refs = 0;
    uint8_t *dec = NULL;
    uint64_t doff[2];
    assert(ttio_sam_tags_decode(&b.ctx, out, ol, &dec, doff) == TTIO_RANS_ERR_CORRUPT);
    /* n_reads must match the blob. */
    b.ctx.n_refs = 1;
    b.ctx.n_reads = 2;
    assert(ttio_sam_tags_decode(&b.ctx, out, ol, &dec, doff) == TTIO_RANS_ERR_CORRUPT);
    b.ctx.n_reads = 1;
    out[0] = 'X';
    assert(ttio_sam_tags_decode(&b.ctx, out, ol, &dec, doff) == TTIO_RANS_ERR_CORRUPT);
    ttio_sam_tags_free(out);
    unbuild(&b);
    printf("test_bad_input: PASS\n");
}

int main(void) {
    test_canonical_int();
    test_md_nm();
    test_empty_run();
    test_mixed();
    test_columns();
    test_bad_input();
    printf("ALL sam_tags tests PASS\n");
    return 0;
}
