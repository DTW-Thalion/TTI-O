#!/usr/bin/env bash
# M103: reads grouped by sequence before blocking (README, "Grouping reads
# before blocking"). Builds native/tools/seq_cm_fastq.c against the kernel in
# ~/m103/build and runs it on the HG002 2x250 chr22 reads in name-hash
# (collate) order and the WES chr22 slice. hg002.collp.pos holds each read's
# mapping position in the same order, for the locality diagnostic only.
set -euo pipefail
N="$(cd "$(dirname "$0")/../../../native" && pwd)"
cd ~/m103
gcc -O3 -march=native -Wall -Wextra -o seqcm_group "$N/tools/seq_cm_fastq.c" -I"$N/include" \
    -Lbuild -lttio_rans -Wl,-rpath,"$HOME/m103/build" -lm
if [ ! -s hg002.collp.fq ]; then
  samtools collate -u -@4 -o /tmp/hgcoll.bam ~/m101/hg002_2x250.chr22.bam /tmp/collp
  samtools fastq -n /tmp/hgcoll.bam > hg002.collp.fq
  samtools view -F 0x900 /tmp/hgcoll.bam | cut -f4 > hg002.collp.pos
  rm -f /tmp/hgcoll.bam
fi
if [ ! -s hg002.sorted.fq ]; then
  mkdir -p tmp
  paste - - - - < hg002.collp.fq | paste - hg002.collp.pos | sort -t$'\t' -k5,5n -S4G -T tmp > tmp/hgs.tsv
  cut -f5 tmp/hgs.tsv > hg002.sorted.pos
  cut -f1-4 tmp/hgs.tsv | tr '\t' '\n' > hg002.sorted.fq
  rm tmp/hgs.tsv
fi
echo "== HG002 chr22, coordinate order";  ./seqcm_group --as-is --pos hg002.sorted.pos < hg002.sorted.fq
echo "== HG002 chr22, input order";       ./seqcm_group --as-is --pos hg002.collp.pos < hg002.collp.fq
echo "== HG002 chr22, single k-mer";      ./seqcm_group --group --pos hg002.collp.pos < hg002.collp.fq
echo "== HG002 chr22, single k-mer pairs"; ./seqcm_group --group --pairs < hg002.collp.fq
echo "== HG002 chr22, chained";           ./seqcm_group --chain --pos hg002.collp.pos < hg002.collp.fq
echo "== HG002 chr22, chained, w 16";    ./seqcm_group --chain --group-w 16 --pos hg002.collp.pos < hg002.collp.fq
for opt in "--fill" "--fill --mates" "--fill --mates --scaffold"; do
  echo "== HG002 chr22, chained $opt"; ./seqcm_group --chain $opt --pos hg002.collp.pos < hg002.collp.fq
done
echo "== HG002 chr22, --fill --mates, blocks sorted by position (diagnostic)"
./seqcm_group --chain --fill --mates --oracle-sort --pos hg002.collp.pos < hg002.collp.fq
echo "== HG002 chr22, --layout (locality only)"; ./seqcm_group --layout --no-code --pos hg002.collp.pos < hg002.collp.fq
for f in wes.shuf.fq wes.sorted.fq; do
  echo "== WES chr22 $f"
  ./seqcm_group < $f
  ./seqcm_group --group < $f
  ./seqcm_group --group --pairs < $f
  ./seqcm_group --chain --group-w 32 < $f
  ./seqcm_group --chain < $f
  ./seqcm_group --chain --min-votes 1 < $f
  ./seqcm_group --chain --fill < $f
  ./seqcm_group --chain --fill --mates --scaffold < $f
  ./seqcm_group --layout < $f
done
