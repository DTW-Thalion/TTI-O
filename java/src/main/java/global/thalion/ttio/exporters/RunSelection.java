/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.exporters;

import global.thalion.ttio.AcquisitionRun;
import global.thalion.ttio.Enums;
import global.thalion.ttio.SpectralDataset;
import global.thalion.ttio.genomics.GenomicIndex;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.WrittenGenomicRun;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

/**
 * Run-selection helpers shared by the export registry and the per-format
 * {@code Writer} adapters. Each picks the right {@link AcquisitionRun}
 * (or {@link GenomicRun}) from an opened {@link SpectralDataset} given an
 * optional {@code layer} name.
 *
 * <p>This is the Java port of the Python reference
 * {@code ttio.exporters._select} ({@code analytical_run} /
 * {@code nmr_run} / {@code genomic_run}). The error-message strings are
 * kept byte-identical to Python so cross-language error parity holds.</p>
 *
 * <p><b>Structural difference from Python:</b> Python merges two run maps
 * ({@code ds.ms_runs} and {@code ds.nmr_runs}). The Java
 * {@link SpectralDataset} exposes a single analytical-run map via
 * {@link SpectralDataset#msRuns()}; there is no separate
 * {@code nmrRuns()} accessor. NMR runs live inside {@code msRuns()} and
 * are distinguished by
 * {@code Enums.SpectrumKind.fromPersisted(run.spectrumClassName()) == NMR}
 * (the same discriminant Python applies via {@code spectrum_class}).
 * Selection behaviour is therefore equivalent to Python's intent.</p>
 *
 * <p>Exception mapping: Python's {@code KeyError} (not-found / ambiguous)
 * maps to {@link IllegalArgumentException} here — matching how
 * {@code exporters} / {@code ExportTask} already signal missing runs.</p>
 *
 * @since 1.7.0
 */
public final class RunSelection {

    private RunSelection() { }

    /**
     * Select an analytical run (any spectrum class — MS or NMR) by
     * {@code layer} name, or the single run when unambiguous. Mirrors
     * Python's {@code analytical_run}.
     *
     * @throws IllegalArgumentException if there are no analytical runs,
     *     the named {@code layer} is absent, or {@code layer} is null and
     *     multiple runs are present (ambiguous).
     */
    public static AcquisitionRun analyticalRun(SpectralDataset ds, String layer) {
        Map<String, AcquisitionRun> runs = ds.msRuns();
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no analytical runs in dataset");
        }
        if (layer != null && !layer.isEmpty()) {
            AcquisitionRun r = runs.get(layer);
            if (r == null) {
                throw new IllegalArgumentException(
                    "run '" + layer + "' not found; have: " + sortedNames(runs));
            }
            return r;
        }
        if (runs.size() == 1) {
            return runs.values().iterator().next();
        }
        throw new IllegalArgumentException(
            "multiple runs present; pass --layer <name>");
    }

    /**
     * Select an NMR run, preferring the NMR-classed run, falling back to
     * the sole analytical run. Mirrors Python's {@code nmr_run}.
     *
     * @throws IllegalArgumentException if there are no analytical runs,
     *     the named {@code layer} is absent, or the choice is ambiguous.
     */
    public static AcquisitionRun nmrRun(SpectralDataset ds, String layer) {
        Map<String, AcquisitionRun> runs = ds.msRuns();
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no analytical runs in dataset");
        }
        if (layer != null && !layer.isEmpty()) {
            AcquisitionRun r = runs.get(layer);
            if (r == null) {
                throw new IllegalArgumentException(
                    "run '" + layer + "' not found; have: " + sortedNames(runs));
            }
            return r;
        }
        List<AcquisitionRun> nmr = new ArrayList<>();
        for (AcquisitionRun r : runs.values()) {
            if (Enums.SpectrumKind.fromPersisted(r.spectrumClassName())
                    == Enums.SpectrumKind.NMR) {
                nmr.add(r);
            }
        }
        if (nmr.size() == 1) {
            return nmr.get(0);
        }
        if (nmr.size() > 1) {
            throw new IllegalArgumentException(
                "multiple NMR runs present; pass --layer <name>");
        }
        if (runs.size() == 1) {
            return runs.values().iterator().next();
        }
        throw new IllegalArgumentException(
            "multiple runs present; pass --layer <name>");
    }

    /**
     * Select a genomic run by {@code layer} name, or the single run when
     * unambiguous. Mirrors Python's {@code genomic_run}.
     *
     * @throws IllegalArgumentException if there are no genomic runs, the
     *     named {@code layer} is absent, or {@code layer} is null and
     *     multiple genomic runs are present (ambiguous).
     */
    public static GenomicRun genomicRun(SpectralDataset ds, String layer) {
        Map<String, GenomicRun> runs = ds.genomicRuns();
        if (runs.isEmpty()) {
            throw new IllegalArgumentException("no genomic runs in dataset");
        }
        if (layer != null && !layer.isEmpty()) {
            GenomicRun r = runs.get(layer);
            if (r == null) {
                throw new IllegalArgumentException(
                    "genomic run '" + layer + "' not found; have: "
                    + sortedNames(runs));
            }
            return r;
        }
        if (runs.size() == 1) {
            return runs.values().iterator().next();
        }
        throw new IllegalArgumentException(
            "multiple genomic runs present; pass --layer <name>: "
            + sortedNames(runs));
    }

    /**
     * Materialise a read-side {@link GenomicRun} into a write-side
     * {@link WrittenGenomicRun} for BAM / CRAM export.
     *
     * <p>This is the single shared copy of the conversion the tio-browser
     * GUI performs in {@code ExportTask.toWritten}; the {@code BamWriterAdapter}
     * and {@code CramWriterAdapter} both call it so the materialisation logic
     * lives in exactly one place (and PR-J2's GUI can adopt it later).</p>
     */
    public static WrittenGenomicRun toWritten(GenomicRun run) {
        int n = run.readCount();
        GenomicIndex idx = run.index();
        long[] positions = new long[n];
        byte[] mapqs     = new byte[n];
        int[]  flags     = new int[n];
        long[] offsets   = new long[n];
        int[]  lengths   = new int[n];
        List<String> chromosomes = new ArrayList<>(n);
        List<String> readNames   = new ArrayList<>(n);
        List<String> cigars      = new ArrayList<>(n);
        List<String> mateChroms  = new ArrayList<>(n);
        long[] matePos   = new long[n];
        int[]  tlens     = new int[n];
        // A grouped run (M103) presents input order everywhere, but its
        // per-read accessors hop between stored blocks: walk its
        // input-order iterator instead, which decodes each block once.
        java.util.Iterator<global.thalion.ttio.genomics.AlignedRead> grouped =
            "blocks_v1_grouped".equals(run.layout()) ? run.iterReads() : null;
        for (int i = 0; i < n; i++) {
            positions[i] = idx.positionAt(i);
            mapqs[i]     = (byte) idx.mappingQualityAt(i);
            flags[i]     = idx.flagsAt(i);
            offsets[i]   = idx.offsetAt(i);
            lengths[i]   = idx.lengthAt(i);
            chromosomes.add(idx.chromosomeAt(i));
            if (grouped != null) {
                global.thalion.ttio.genomics.AlignedRead r = grouped.next();
                readNames.add(r.readName());
                cigars.add(r.cigar());
                mateChroms.add(r.mateChromosome());
                matePos[i] = r.matePosition();
                tlens[i]   = r.templateLength();
                continue;
            }
            readNames.add(run.readNameAt(i));
            cigars.add(run.cigarAt(i));
            mateChroms.add(run.mateChromAt(i));
            matePos[i]   = run.matePosAt(i);
            tlens[i]     = run.mateTlenAt(i);
        }
        byte[] seqs  = n > 0 ? run.sequencesFull() : new byte[0];
        byte[] quals = n > 0 ? run.qualitiesFull() : new byte[0];
        return new WrittenGenomicRun(
            run.acquisitionMode() != null
                ? run.acquisitionMode() : Enums.AcquisitionMode.GENOMIC_WGS,
            run.referenceUri() != null ? run.referenceUri() : "",
            run.platform() != null ? run.platform() : "",
            run.sampleName() != null ? run.sampleName() : "",
            positions, mapqs, flags,
            seqs, quals,
            offsets, lengths,
            cigars, readNames, mateChroms, matePos, tlens, chromosomes,
            Enums.Compression.NONE
        );
    }

    /** {@code ", ".join(sorted(runs))} — sorted, comma-space-joined keys. */
    private static String sortedNames(Map<String, ?> runs) {
        return String.join(", ", new TreeSet<>(runs.keySet()));
    }
}
