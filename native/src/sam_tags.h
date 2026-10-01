/*
 * sam_tags.h — internal helpers of the SAM_TAGS codec (codec id 18).
 *
 * The public entry points (ttio_sam_tags_encode / _decode / _free) and
 * the context struct live in include/ttio_rans.h. These are exposed for
 * the native unit tests. Spec: docs/codecs/sam_tags.md.
 */
#ifndef TTIO_SAM_TAGS_H
#define TTIO_SAM_TAGS_H

#include <stddef.h>
#include <stdint.h>

#include "ttio_rans.h"

#define STG_MAGIC "STG1"
#define STG_VERSION 1

enum {
    STG_KIND_INT = 0,
    STG_KIND_TEXT = 1,
    STG_KIND_DERIVED = 2,
    STG_KIND_DUP = 3,
    STG_KIND_VERBATIM = 4
};

/* Recompute MD and NM for read `i` (spec §1.3). On success writes the
 * NUL-terminated MD into a malloc'd *md (caller frees), its length into
 * *md_len, NM into *nm, and returns 1. Returns 0 when the read is not
 * derivable, -1 on allocation failure. */
int stg_calc_md_nm(const ttio_sam_tags_ctx *ctx, uint64_t i,
                   char **md, size_t *md_len, int64_t *nm);

/* 1 when s[0..n) is a canonical decimal int64 (spec §1.2), with the
 * value in *out. */
int stg_parse_canonical_int(const uint8_t *s, size_t n, int64_t *out);

#endif /* TTIO_SAM_TAGS_H */
