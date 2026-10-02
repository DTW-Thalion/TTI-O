#!/usr/bin/env bash
# M103 Phase 0 on HG002 2x250 chr22 (10.6 M reads, 2.64 G bases). Reads in
# samtools-collate order (grouped by a hash of the read name, so effectively
# shuffled with mates adjacent). Needs ~seqp built by run_wes.sh; 2^28 tables
# need about 6.5 GB of RAM.
set -uo pipefail
cd ~/m103
for args in "--rc 11 16 24" "--rc --table-bits 26 11 16 24" "--rc --table-bits 26 12 20" \
            "--rc --table-bits 28 11 16 24" "--rc --table-bits 28 11 16 22"; do
  printf '%-32s ' "$args"
  samtools collate -O -u -@4 ~/m101/hg002_2x250.chr22.bam /tmp/coll 2>/dev/null \
    | samtools fastq -n - 2>/dev/null | ./seqp $args | sed 's/.*bits/bits/'
done
