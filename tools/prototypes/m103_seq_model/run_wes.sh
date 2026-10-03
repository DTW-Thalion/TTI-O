#!/usr/bin/env bash
# M103 Phase 0 sweep on the WES chr22 reads (sorted and shuffled).
set -uo pipefail
SRC="$(dirname "$0")/seq_model_proof.c"
W=~/m103; mkdir -p $W; cd $W
gcc -O3 -march=native -o seqp "$SRC" -lm || exit 1
[ -f wes.sorted.fq ] || samtools fastq -n ~/m101/wes/na12878_wes.chr22.bam > wes.sorted.fq 2>/dev/null
[ -f wes.shuf.fq ] || paste - - - - < wes.sorted.fq | shuf --random-source=<(yes) | tr '\t' '\n' > wes.shuf.fq
for f in wes.shuf.fq wes.sorted.fq; do
  echo "== $f"
  for args in "12" "16" "20" "24" "--rc 12" "--rc 16" "--rc 20" "--rc 24" "--rc 12 20" "--rc 11 16 24" "--rc 12 24"; do
    printf '%-16s ' "$args"; ./seqp $args < $f | sed 's/.*bits/bits/'
  done
done
