/*
 * seq_cm.h — internal helpers of the SEQ_CM codec (codec id 19, M103):
 * a reference-free context-mixing model for read bases.
 *
 * The public entry points (ttio_seq_cm_encode / _decode / _free) and the
 * parameter struct live in include/ttio_rans.h. These are exposed for the
 * native unit tests. Spec: docs/codecs/seq_cm.md.
 */
#ifndef TTIO_SEQ_CM_H
#define TTIO_SEQ_CM_H

#include <stddef.h>
#include <stdint.h>

#include "ttio_rans.h"

#define SQC_MAGIC "SQC1"
#define SQC_VERSION 1
#define SQC_HEADER_LEN 48

#define SQC_MIN_TABLE_BITS 10
#define SQC_MAX_TABLE_BITS 30

/* Logistic helpers in the 12-bit domain (spec §3.2): squash maps a
 * stretch value in [-2047, 2047] to a probability in [0, 4095]; stretch
 * is its inverse. Integer only, so every platform codes the same bytes. */
int sqc_squash(int d);
int sqc_stretch(int p12);

/* The table_bits an encoder picks for n_bases when the caller passes 0
 * (spec §2.3). */
uint8_t sqc_auto_table_bits(uint64_t n_bases);

#endif /* TTIO_SEQ_CM_H */
