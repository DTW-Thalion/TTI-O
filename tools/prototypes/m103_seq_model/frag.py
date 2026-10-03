"""How a region's reads are fragmented in a grouped order.

usage: frag.py ORDER POS
ORDER: one input read index per line (grouped order); POS: mapping position
per input read (0 = unmapped). Reports, for a few 100 kb bins, the runs of
consecutive output positions holding that bin's reads, and overall run
statistics: a run is a maximal stretch of output order whose reads all map
within 100 kb of each other.
"""
import sys
from array import array

order = array('I', (int(l) for l in open(sys.argv[1])))
pos = array('q', (int(l) for l in open(sys.argv[2])))
n = len(order)
print('reads', n)

# Runs: break the output order wherever consecutive mapped reads are > 2 kb apart.
runs = []  # (start index in output, length, min pos, max pos)
cur_start, cur_min, cur_max, prev = 0, None, None, None
for i, r in enumerate(order):
    p = pos[r]
    if p <= 0:
        continue
    if prev is not None and abs(p - prev) > 2000:
        runs.append((cur_start, i - cur_start, cur_min, cur_max))
        cur_start, cur_min, cur_max = i, p, p
    else:
        if cur_min is None:
            cur_min = cur_max = p
        cur_min, cur_max = min(cur_min, p), max(cur_max, p)
    prev = p
runs.append((cur_start, n - cur_start, cur_min, cur_max))
lens = sorted(r[1] for r in runs)
spans = sorted((r[3] - r[2]) for r in runs if r[2] is not None)
def q(v, f):
    return v[min(len(v) - 1, int(f * len(v)))]
print('runs (breaks at > 2 kb jumps):', len(runs))
print('run length in reads: median', q(lens, 0.5), 'p90', q(lens, 0.9), 'p99', q(lens, 0.99), 'max', lens[-1])
print('run genome span bp: median', q(spans, 0.5), 'p90', q(spans, 0.9), 'p99', q(spans, 0.99), 'max', spans[-1])
tot = sum(lens)
big = sum(l for l in lens if l >= 1000)
print('share of reads in runs of >= 1000 reads: %.1f%%' % (100.0 * big / tot))
big = sum(l for l in lens if l >= 100)
print('share of reads in runs of >= 100 reads: %.1f%%' % (100.0 * big / tot))

# Jump distances between runs.
jumps = []
for a, b in zip(runs, runs[1:]):
    if a[3] is not None and b[2] is not None:
        jumps.append(abs(b[2] - a[3]))
jumps.sort()
print('jump between runs bp: median', q(jumps, 0.5), 'p10', q(jumps, 0.1), 'p90', q(jumps, 0.9))

# A few bins: where their reads sit in the output.
where = {}
for i, r in enumerate(order):
    p = pos[r]
    if p > 0:
        b = p // 100000
        where.setdefault(b, []).append(i)
for b in sorted(where)[100::100][:4]:
    idx = where[b]
    segs, s0, last = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - last > 50:
            segs.append((s0, last - s0 + 1))
            s0 = i
        last = i
    segs.append((s0, last - s0 + 1))
    segs.sort(key=lambda t: -t[1])
    print('bin %d (%d reads): %d stretches in output; largest %s' % (b, len(idx), len(segs), [t[1] for t in segs[:8]]))
