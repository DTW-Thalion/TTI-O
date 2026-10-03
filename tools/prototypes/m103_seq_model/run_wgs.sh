#!/usr/bin/env bash
# M103 Phase 0 on a sequencer-order whole-genome FASTQ: GIAB HG002 HiSeq 2500
# 2x148, flowcell 140528_D00360_0018_AH8VC6ADXX, sample 2A1, lanes 1-2
# (8 files, ~1x genome). Files are read in sequencer order: each chunk's R1
# then R2. xz -9 on the bare sequence lines is the general-purpose baseline.
set -uo pipefail
D="${1:-/mnt/c/Users/toddw/Documents/Business/Thalion Initiative/BAFAR Source/TTI-O/data/genomic/hg002_hiseq_wgs}"
cd ~/m103
reads() {
  for l in L001 L002; do for c in 001 002; do for r in R1 R2; do
    zcat "$D/2A1_CGATGT_${l}_${r}_${c}.fastq.gz"
  done; done; done
}
printf '%-32s ' "xz -9 sequence lines"
reads | awk 'NR % 4 == 2' | xz -9 -T4 | wc -c
for args in "--rc 11 16 24" "--rc --table-bits 26 11 16 24" "--rc --table-bits 28 11 16 24" \
            "--rc --table-bits 28 12 20" "--rc 16"; do
  printf '%-32s ' "$args"
  reads | ./seqp $args | sed 's/.*reads/reads/'
done
