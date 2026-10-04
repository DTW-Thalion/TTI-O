"""Genomic read benchmark: does the decode-ahead window bind?

Usage
-----
::

    python -m ttio.tools.genomic_read_bench make FILE N_READS READ_LEN BLOCK_READS
    python -m ttio.tools.genomic_read_bench read FILE [WINDOWS] [RUN]

``make`` writes a blocks_v1 run; ``read`` iterates it end to end once
per decode-ahead window in ``WINDOWS`` (default ``1,2,4,8,16``) and
prints one ``[py-bench]`` line each.

The point of the sweep is what it does *not* show. Decoding a block is
parallel and fast; consuming one is a serial Python loop building an
AlignedRead per record. If the window is what limits a reader, raising
it raises the rate; if the consumer is, the rate is flat and the window
is holding memory for nothing. ``v6_acceptance`` cannot answer this
because it consumes at memory speed, which no caller does.

Rows are best-of-N over interleaved rounds, because a single run on a
loaded machine can be a third slow. A control row repeats one window
under a second name: when it disagrees with its twin by as much as the
windows disagree with each other, the sweep says nothing.

SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np


def _build_run(n_reads: int, read_len: int):
    """A run whose qualities are worth compressing: per-read quality
    means drift and decline along the read, as real ones do, so the
    codec models something rather than a flat line."""
    from ..written_genomic_run import WrittenGenomicRun

    rng = np.random.default_rng(11)
    total = n_reads * read_len
    ramp = np.linspace(38.0, 26.0, read_len, dtype=np.float32)
    per_read = rng.normal(0.0, 2.5, n_reads).astype(np.float32)
    q = ramp[None, :] + per_read[:, None] + rng.normal(0.0, 2.0, (n_reads, read_len)).astype(np.float32)
    quals = np.clip(q, 2, 41).astype(np.uint8).reshape(-1)

    seqs = np.frombuffer(b"ACGT", dtype=np.uint8)[rng.integers(0, 4, total)].astype(np.uint8)
    lengths = np.full(n_reads, read_len, dtype=np.uint32)
    offsets = (np.arange(n_reads, dtype=np.uint64) * np.uint64(read_len))
    positions = (np.arange(n_reads, dtype=np.int64) * 10 + 1)

    return WrittenGenomicRun(
        acquisition_mode=7,
        reference_uri="synthetic.ref",
        platform="ILLUMINA",
        sample_name="bench",
        positions=positions,
        mapping_qualities=np.full(n_reads, 60, dtype=np.uint8),
        flags=np.zeros(n_reads, dtype=np.uint32),
        sequences=seqs,
        qualities=quals,
        offsets=offsets,
        lengths=lengths,
        cigars=[f"{read_len}M"] * n_reads,
        read_names=[f"r{i}" for i in range(n_reads)],
        mate_chromosomes=["*"] * n_reads,
        mate_positions=np.full(n_reads, -1, dtype=np.int64),
        template_lengths=np.zeros(n_reads, dtype=np.int32),
        chromosomes=["chr1"] * n_reads,
    )


def _make(path: Path, n_reads: int, read_len: int, block_reads: int) -> int:
    from ..genomic import GenomicStreamWriter
    from ..spectral_dataset import SpectralDataset

    run = _build_run(n_reads, read_len)
    SpectralDataset.write_minimal(str(path), title="bench",
                                  isa_investigation_id="i", runs={})
    t0 = time.monotonic()
    ds = SpectralDataset.open(str(path), writable=True)
    with ds, GenomicStreamWriter(ds.study_group, "run",
                                 acquisition_mode=run.acquisition_mode,
                                 reference_uri=run.reference_uri,
                                 platform=run.platform,
                                 sample_name=run.sample_name,
                                 block_reads=block_reads) as w:
        w.append_batch(run)
    dt = time.monotonic() - t0
    print(f"[py-bench] wrote reads={n_reads} read_len={read_len} "
          f"block_reads={block_reads} seconds={dt:.1f} "
          f"file_bytes={path.stat().st_size}")
    return 0


def _iterate(path: Path, run_name: str | None, window: int) -> tuple[int, int, float]:
    """One full pass. Touches every quality byte, so the decode cannot
    be skipped and the consumer cost is the real one."""
    from .. import genomic_run as _gr
    from ..spectral_dataset import SpectralDataset

    saved = _gr._READ_AHEAD_BLOCKS
    _gr._READ_AHEAD_BLOCKS = window
    try:
        with SpectralDataset.open(str(path)) as ds:
            name = run_name or next(iter(ds.genomic_runs))
            g = ds.genomic_runs[name]
            t0 = time.monotonic()
            n = 0
            nq = 0
            for r in g.iter_reads():
                nq += len(r.qualities)
                n += 1
            dt = time.monotonic() - t0
        return n, nq, dt
    finally:
        _gr._READ_AHEAD_BLOCKS = saved


def _read(path: Path, windows: list[int], run_name: str | None,
          rounds: int = 3) -> int:
    from .. import genomic_run as _gr
    from .._threads import resolve_threads
    from ..spectral_dataset import SpectralDataset

    with SpectralDataset.open(str(path)) as ds:
        name = run_name or next(iter(ds.genomic_runs))
        g = ds.genomic_runs[name]
        layout, blocks, n_reads = g.layout, g.block_count, len(g)
    if layout not in ("blocks_v1", "blocks_v1_grouped"):
        print(f"{path} run {name} is {layout}, not blocks_v1; the window "
              f"does not apply", file=sys.stderr)
        return 1

    print(f"[py-bench] file={path.name} run={name} layout={layout} "
          f"blocks={blocks} reads={n_reads} threads={resolve_threads(None)} "
          f"default_window={_gr._READ_AHEAD_BLOCKS} rounds={rounds}")

    # The control repeats the first window, so a flat sweep can be told
    # apart from a machine too noisy to show anything.
    plan = [(str(w), w) for w in windows] + [("control", windows[0])]
    best: dict[str, tuple[int, int, float]] = {}
    for _ in range(rounds):
        for label, w in plan:
            n, nq, dt = _iterate(path, name, w)
            prev = best.get(label)
            if prev is None or dt < prev[2]:
                best[label] = (n, nq, dt)

    for label, w in plan:
        n, nq, dt = best[label]
        print(f"[py-bench] window={label:<8} reads={n} seconds={dt:.2f} "
              f"reads_per_s={n / dt:,.0f} qual_mb_per_s={nq / dt / 1e6:.1f}")

    first, ctl = best[str(windows[0])][2], best["control"][2]
    drift = abs(ctl / first - 1) * 100
    slowest = max(best[str(w)][2] for w in windows)
    fastest = min(best[str(w)][2] for w in windows)
    spread = (slowest / fastest - 1) * 100
    print(f"[py-bench] control drift={drift:.1f}% window spread={spread:.1f}%")
    if spread <= drift * 1.5:
        print("[py-bench] verdict: the window does NOT bind — the spread "
              "across windows is within the drift between two runs of the "
              "same window")
    else:
        print("[py-bench] verdict: the window BINDS — widening it changes "
              "the rate by more than the machine's own drift")
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__, file=sys.stderr)
        return 1
    if argv[0] == "make":
        if len(argv) != 5:
            print("usage: genomic_read_bench make FILE N_READS READ_LEN "
                  "BLOCK_READS", file=sys.stderr)
            return 1
        return _make(Path(argv[1]), int(argv[2]), int(argv[3]), int(argv[4]))
    if argv[0] == "read":
        if not 2 <= len(argv) <= 4:
            print("usage: genomic_read_bench read FILE [WINDOWS] [RUN]",
                  file=sys.stderr)
            return 1
        windows = ([int(x) for x in argv[2].split(",")]
                   if len(argv) > 2 else [1, 2, 4, 8, 16])
        return _read(Path(argv[1]), windows, argv[3] if len(argv) > 3 else None)
    print(__doc__, file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
