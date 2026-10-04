/* TTI-O Java Implementation / Copyright (c) 2026 The Thalion Initiative / SPDX-License-Identifier: LGPL-3.0-or-later */
package global.thalion.ttio.protection;

import global.thalion.ttio.*;
import global.thalion.ttio.Enums.*;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Applies caller-selected policies to a dataset and writes a new
 * {@code .tio} file. The original is never modified.
 *
 * <p>Supported policies (matching ObjC {@code TTIOAnonymizationPolicy}):
 * redact SAAV spectra, mask low-intensity samples, mask rare
 * metabolites, coarsen m/z decimals, coarsen chemical-shift
 * decimals, strip metadata fields.</p>
 *
 * <p>Output carries the {@code opt_anonymized} feature flag and a
 * {@link global.thalion.ttio.ProvenanceRecord} documenting which
 * policies ran.</p>
 *
 * <p><b>API status:</b> Stable.</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOAnonymizer}, Python {@code ttio.anonymization}.</p>
 *
 *
 */
public class Anonymizer {

    /** a single masked region (chromosome + half-open
     *  interval). Reads whose mapping position falls inside any
     *  registered region have their sequence + qualities zeroed, but
     *  their index entries are preserved so downstream readers see
     *  the same read count and per-read offsets. */
    public record MaskRegion(String chromosome, long start, long end) {}

    public record AnonymizationPolicy(
        boolean redactSaavSpectra,
        double maskIntensityBelowQuantile,  // 0.0 = disabled
        boolean maskRareMetabolites,
        double rareMetaboliteThreshold,     // default 0.05
        int coarsenMzDecimals,              // -1 = disabled
        int coarsenChemicalShiftDecimals,   // -1 = disabled
        boolean stripMetadata,
        // M90.3 genomic policies. None / false / null = no-op.
        boolean stripReadNames,
        boolean randomiseQualities,
        int randomiseQualitiesConstant,    // default 30
        java.util.List<MaskRegion> maskRegions,  // null = no masking
        /** when non-null, qualities are replaced with
         *  deterministic random Phred bytes in [0, 93] from a seeded
         *  RNG. When null, the M90.3 constant path is used. */
        Long randomiseQualitiesSeed
    ) {
        /** Backward-compatible 7-arg constructor for the pre-M90.3
         *  surface. Defaults the four genomic fields to no-op
         *  (genomic policies disabled, randomise constant 30, no
         *  mask regions). Exists so the M90 work is purely additive
         *  for callers that don't touch genomic content. */
        public AnonymizationPolicy(boolean redactSaavSpectra,
                                     double maskIntensityBelowQuantile,
                                     boolean maskRareMetabolites,
                                     double rareMetaboliteThreshold,
                                     int coarsenMzDecimals,
                                     int coarsenChemicalShiftDecimals,
                                     boolean stripMetadata) {
            this(redactSaavSpectra, maskIntensityBelowQuantile,
                 maskRareMetabolites, rareMetaboliteThreshold,
                 coarsenMzDecimals, coarsenChemicalShiftDecimals,
                 stripMetadata,
                 false, false, 30, null, null);
        }

        /** Backward-compatible 11-arg constructor for the M90.3
         *  surface. Defaults {@code randomiseQualitiesSeed} to
         *  {@code null} (constant-replacement path). M90.14 callers
         *  use the canonical 12-arg constructor instead. */
        public AnonymizationPolicy(boolean redactSaavSpectra,
                                     double maskIntensityBelowQuantile,
                                     boolean maskRareMetabolites,
                                     double rareMetaboliteThreshold,
                                     int coarsenMzDecimals,
                                     int coarsenChemicalShiftDecimals,
                                     boolean stripMetadata,
                                     boolean stripReadNames,
                                     boolean randomiseQualities,
                                     int randomiseQualitiesConstant,
                                     java.util.List<MaskRegion> maskRegions) {
            this(redactSaavSpectra, maskIntensityBelowQuantile,
                 maskRareMetabolites, rareMetaboliteThreshold,
                 coarsenMzDecimals, coarsenChemicalShiftDecimals,
                 stripMetadata,
                 stripReadNames, randomiseQualities,
                 randomiseQualitiesConstant, maskRegions, null);
        }

        public static AnonymizationPolicy defaults() {
            return new AnonymizationPolicy(
                true, 0.0, false, 0.05, -1, -1, true);
        }
    }

    public record AnonymizationResult(
        SpectralDataset dataset,
        int spectraRedacted,
        int intensitiesZeroed,
        int mzValuesCoarsened,
        int chemicalShiftValuesCoarsened,
        int metabolitesMasked,
        int metadataFieldsStripped,
        // M90.3 genomic counters. Zero = either no genomic content
        // or the corresponding policy was disabled.
        int readNamesStripped,
        int qualitiesRandomised,
        int readsInMaskedRegion
    ) {
        /** Backward-compatible 7-arg constructor for the pre-M90.3
         *  surface. Defaults the three genomic counters to 0. */
        public AnonymizationResult(SpectralDataset dataset,
                                     int spectraRedacted,
                                     int intensitiesZeroed,
                                     int mzValuesCoarsened,
                                     int chemicalShiftValuesCoarsened,
                                     int metabolitesMasked,
                                     int metadataFieldsStripped) {
            this(dataset, spectraRedacted, intensitiesZeroed,
                 mzValuesCoarsened, chemicalShiftValuesCoarsened,
                 metabolitesMasked, metadataFieldsStripped,
                 0, 0, 0);
        }
    }

    /** Apply anonymization policies to a dataset, writing a new anonymized file. */
    public static AnonymizationResult anonymize(
            SpectralDataset source, String outputPath, AnonymizationPolicy policy) {
        return anonymize(source, outputPath, policy, loadDefaultPrevalence());
    }

    public static AnonymizationResult anonymize(
            SpectralDataset source, String outputPath,
            AnonymizationPolicy policy, Map<String, Double> prevalenceTable) {

        int spectraRedacted = 0;
        int intensitiesZeroed = 0;
        int mzValuesCoarsened = 0;
        int csValuesCoarsened = 0;
        int metabolitesMasked = 0;
        int metadataFieldsStripped = 0;

        // Build set of SAAV spectrum indices to redact
        Set<String> saavRunSpecs = new HashSet<>();
        if (policy.redactSaavSpectra()) {
            for (Identification id : source.identifications()) {
                if (isSaav(id.chemicalEntity())) {
                    saavRunSpecs.add(id.runName() + ":" + id.spectrumIndex());
                }
            }
            spectraRedacted = saavRunSpecs.size();
        }

        // Process runs
        List<AcquisitionRun> anonymizedRuns = new ArrayList<>();
        for (var entry : source.msRuns().entrySet()) {
            String runName = entry.getKey();
            AcquisitionRun run = entry.getValue();
            SpectrumIndex idx = run.spectrumIndex();

            // Copy channels
            Map<String, double[]> newChannels = new LinkedHashMap<>();
            for (var chEntry : run.channels().entrySet()) {
                newChannels.put(chEntry.getKey(), chEntry.getValue().clone());
            }

            // Per-spectrum processing
            for (int i = 0; i < idx.count(); i++) {
                String specKey = runName + ":" + i;
                long offset = idx.offsetAt(i);
                int length = idx.lengthAt(i);

                // SAAV redaction: zero out both channels
                if (policy.redactSaavSpectra() && saavRunSpecs.contains(specKey)) {
                    for (double[] ch : newChannels.values()) {
                        Arrays.fill(ch, (int) offset, (int) offset + length, 0.0);
                    }
                    continue; // skip further processing for redacted spectra
                }

                // Rare metabolite masking
                if (policy.maskRareMetabolites()) {
                    for (Identification id : source.identifications()) {
                        if (id.runName().equals(runName) && id.spectrumIndex() == i) {
                            Double prev = prevalenceTable.get(id.chemicalEntity());
                            if (prev != null && prev < policy.rareMetaboliteThreshold()) {
                                double[] intensity = newChannels.get("intensity");
                                if (intensity != null) {
                                    Arrays.fill(intensity, (int) offset,
                                            (int) offset + length, 0.0);
                                    metabolitesMasked++;
                                }
                            }
                        }
                    }
                }
            }

            // Intensity quantile masking (global across run)
            if (policy.maskIntensityBelowQuantile() > 0) {
                double[] intensity = newChannels.get("intensity");
                if (intensity != null) {
                    double threshold = quantileThreshold(intensity,
                            policy.maskIntensityBelowQuantile());
                    for (int i = 0; i < intensity.length; i++) {
                        if (intensity[i] < threshold && intensity[i] != 0) {
                            intensity[i] = 0.0;
                            intensitiesZeroed++;
                        }
                    }
                }
            }

            // m/z coarsening
            if (policy.coarsenMzDecimals() >= 0) {
                double[] mz = newChannels.get("mz");
                if (mz != null) {
                    double factor = Math.pow(10, policy.coarsenMzDecimals());
                    for (int i = 0; i < mz.length; i++) {
                        double rounded = Math.round(mz[i] * factor) / factor;
                        if (rounded != mz[i]) mzValuesCoarsened++;
                        mz[i] = rounded;
                    }
                }
            }

            // Chemical shift coarsening
            if (policy.coarsenChemicalShiftDecimals() >= 0) {
                double[] cs = newChannels.get("chemical_shift");
                if (cs != null) {
                    double factor = Math.pow(10, policy.coarsenChemicalShiftDecimals());
                    for (int i = 0; i < cs.length; i++) {
                        double rounded = Math.round(cs[i] * factor) / factor;
                        if (rounded != cs[i]) csValuesCoarsened++;
                        cs[i] = rounded;
                    }
                }
            }

            anonymizedRuns.add(new AcquisitionRun(runName, run.acquisitionMode(),
                    idx, run.instrumentConfig(), newChannels,
                    run.chromatograms(), run.provenanceRecords(),
                    run.nucleusType(), run.spectrometerFrequencyMHz()));
        }

        // Metadata stripping
        String title = source.title();
        if (policy.stripMetadata()) {
            title = "";
            metadataFieldsStripped = 1;
        }

        // walk genomic_runs and apply genomic policies. Returns
        // an empty list when the source carries no genomic runs (and
        // create()'s genomic-runs branch is then a no-op). The three
        // genomic counters in the result mirror the three policies.
        int[] genomicCounters = new int[3];  // [stripped, randomised, masked]
        List<WrittenGenomicRun> anonymizedGenomicRuns =
            applyGenomicPolicies(source, policy, genomicCounters);

        // Build provenance record
        Map<String, String> params = new LinkedHashMap<>();
        params.put("spectra_redacted", String.valueOf(spectraRedacted));
        params.put("intensities_zeroed", String.valueOf(intensitiesZeroed));
        params.put("mz_values_coarsened", String.valueOf(mzValuesCoarsened));
        params.put("chemical_shift_values_coarsened", String.valueOf(csValuesCoarsened));
        params.put("metabolites_masked", String.valueOf(metabolitesMasked));
        params.put("metadata_fields_stripped", String.valueOf(metadataFieldsStripped));
        params.put("read_names_stripped", String.valueOf(genomicCounters[0]));
        params.put("qualities_randomised", String.valueOf(genomicCounters[1]));
        params.put("reads_in_masked_region", String.valueOf(genomicCounters[2]));

        ProvenanceRecord anonProv = ProvenanceRecord.of(
                "ttio anonymizer v0.4", params, List.of(), List.of());

        List<ProvenanceRecord> prov = new ArrayList<>(source.provenanceRecords());
        prov.add(anonProv);

        FeatureFlags flags = source.featureFlags().with(FeatureFlags.OPT_ANONYMIZED);

        // Filter identifications (remove SAAV if redacted)
        List<Identification> idents = source.identifications();
        if (policy.redactSaavSpectra()) {
            idents = idents.stream()
                    .filter(id -> !isSaav(id.chemicalEntity()))
                    .toList();
        }

        SpectralDataset result = SpectralDataset.create(outputPath, title,
                source.isaInvestigationId(), anonymizedRuns,
                anonymizedGenomicRuns,
                idents, source.quantifications(), prov, flags);

        return new AnonymizationResult(result, spectraRedacted, intensitiesZeroed,
                mzValuesCoarsened, csValuesCoarsened, metabolitesMasked,
                metadataFieldsStripped,
                genomicCounters[0], genomicCounters[1], genomicCounters[2]);
    }

    // ─────────────────────────────────────── M90.3 genomic policy walker

    /** Walk {@code source.genomicRuns()} and produce a list of
     *  {@link WrittenGenomicRun} copies with {@code policy}-driven
     *  transformations applied (no source mutation).
     *
     *  <p>Returns an empty list if the source carries no genomic
     *  runs OR if no genomic policy is set — the create() path
     *  treats both cases as "no genomic content". The three counters
     *  in {@code countersOut} are populated:
     *  <ul>
     *    <li>{@code [0]} read_names stripped (one per read).</li>
     *    <li>{@code [1]} reads whose qualities were randomised.</li>
     *    <li>{@code [2]} reads whose mapping position fell inside any
     *        registered mask region.</li>
     *  </ul>
     *
     *  <p>When the source has genomic runs but no genomic policy is
     *  set the runs are still copied verbatim — that's required for
     *  parity with the Python reference (a no-op anonymize on a
     *  genomic-bearing file must preserve the genomic content). */
    private static List<WrittenGenomicRun> applyGenomicPolicies(
            SpectralDataset source, AnonymizationPolicy policy,
            int[] countersOut) {
        Map<String, GenomicRun> grs = source.genomicRuns();
        if (grs == null || grs.isEmpty()) {
            return List.of();
        }
        List<WrittenGenomicRun> out = new ArrayList<>(grs.size());
        // Sort by run name so the output is deterministic.
        List<String> names = new ArrayList<>(grs.keySet());
        java.util.Collections.sort(names);
        for (String runName : names) {
            GenomicRun gr = grs.get(runName);
            int n = gr.readCount();

            // Materialise the per-read fields by iterating the lazy
            // GenomicRun. This is O(N reads) but the anonymizer is a
            // one-shot offline tool so the overhead is acceptable.
            List<String> readNames = new ArrayList<>(n);
            List<String> cigars = new ArrayList<>(n);
            byte[][] sequences = new byte[n][];
            byte[][] qualities = new byte[n][];
            List<String> mateChromosomes = new ArrayList<>(n);
            long[] matePositions = new long[n];
            int[] templateLengths = new int[n];
            // Input order, like gr.index() below; the iterator decodes each
            // block once (a grouped run's reads are stored out of input
            // order, M103).
            java.util.Iterator<AlignedRead> reads = gr.iterReads();
            for (int i = 0; i < n; i++) {
                AlignedRead r = reads.next();
                readNames.add(r.readName());
                cigars.add(r.cigar());
                sequences[i] = r.sequence().getBytes(StandardCharsets.US_ASCII);
                qualities[i] = r.qualities().clone();
                mateChromosomes.add(r.mateChromosome());
                matePositions[i] = r.matePosition();
                templateLengths[i] = r.templateLength();
            }

            // ── strip_read_names ────────────────────────────────────
            if (policy.stripReadNames()) {
                for (int i = 0; i < n; i++) readNames.set(i, "");
                countersOut[0] += n;
            }

            // ── randomise_qualities ─────────────────────────────────
            if (policy.randomiseQualities()) {
                if (policy.randomiseQualitiesSeed() != null) {
                    // seeded random Phred per byte. Range
                    // [0, 93] matches the SAM spec valid Phred range
                    // (0 = lowest, 93 = highest representable in
                    // Illumina-style ASCII offset 33 + 60). Cross-
                    // language byte-equality with numpy's PCG64 is
                    // explicitly NOT a goal — Python's docstring
                    // calls out that as a follow-up.
                    java.util.Random rng = new java.util.Random(
                        policy.randomiseQualitiesSeed());
                    for (int i = 0; i < n; i++) {
                        byte[] q = new byte[qualities[i].length];
                        for (int j = 0; j < q.length; j++) {
                            q[j] = (byte) (rng.nextInt(94) & 0xFF);
                        }
                        qualities[i] = q;
                    }
                } else {
                    byte constByte = (byte) (policy.randomiseQualitiesConstant() & 0xFF);
                    for (int i = 0; i < n; i++) {
                        byte[] q = new byte[qualities[i].length];
                        Arrays.fill(q, constByte);
                        qualities[i] = q;
                    }
                }
                countersOut[1] += n;
            }

            // ── mask_regions ────────────────────────────────────────
            if (policy.maskRegions() != null && !policy.maskRegions().isEmpty()) {
                List<String> chroms = new ArrayList<>(n);
                long[] positions = new long[n];
                for (int i = 0; i < n; i++) {
                    chroms.add(gr.index().chromosomeAt(i));
                    positions[i] = gr.index().positionAt(i);
                }
                // SAM-overlap semantics. Walk the CIGAR to
                // compute each read's reference end coordinate; a
                // read overlaps a region iff [pos, pos+span-1]
                // intersects [region_start, region_end] on inclusive
                // endpoints. Falls back to position-only check when
                // CIGAR is empty / "*" / non-parseable (M90.3
                // backward-compat).
                java.util.BitSet alreadyMasked = new java.util.BitSet(n);
                for (MaskRegion region : policy.maskRegions()) {
                    String chrName = region.chromosome();
                    long start = region.start();
                    long end = region.end();
                    for (int i = 0; i < n; i++) {
                        if (alreadyMasked.get(i)) continue;
                        if (!chroms.get(i).equals(chrName)) continue;
                        long pos = positions[i];
                        long span = cigarRefSpan(cigars.get(i));
                        boolean overlaps;
                        if (span > 0) {
                            long readEnd = pos + span - 1;  // inclusive
                            overlaps = !(readEnd < start || pos > end);
                        } else {
                            // Empty / unparseable CIGAR — fall back to
                            // position-only check (behaviour).
                            overlaps = (pos >= start && pos <= end);
                        }
                        if (overlaps) {
                            Arrays.fill(sequences[i], (byte) 0);
                            Arrays.fill(qualities[i], (byte) 0);
                            countersOut[2] += 1;
                            alreadyMasked.set(i);
                        }
                    }
                }
            }

            // Re-pack into the flat WrittenGenomicRun layout.
            int[] lengths = new int[n];
            long[] offsets = new long[n];
            long running = 0;
            int totalSeq = 0, totalQual = 0;
            for (int i = 0; i < n; i++) {
                lengths[i] = sequences[i].length;
                offsets[i] = running;
                running += sequences[i].length;
                totalSeq += sequences[i].length;
                totalQual += qualities[i].length;
            }
            byte[] sequencesFlat = new byte[totalSeq];
            byte[] qualitiesFlat = new byte[totalQual];
            int sCursor = 0, qCursor = 0;
            for (int i = 0; i < n; i++) {
                System.arraycopy(sequences[i], 0, sequencesFlat, sCursor, sequences[i].length);
                sCursor += sequences[i].length;
                System.arraycopy(qualities[i], 0, qualitiesFlat, qCursor, qualities[i].length);
                qCursor += qualities[i].length;
            }

            // Pull the per-read integer fields directly from the index;
            // these are unchanged by anonymization (the index entries
            // are preserved for masked reads — only the sequence /
            // qualities bytes are zeroed).
            long[] positionsOut = new long[n];
            byte[] mappingQualities = new byte[n];
            int[] flagsArr = new int[n];
            List<String> chromosomesList = new ArrayList<>(n);
            for (int i = 0; i < n; i++) {
                positionsOut[i] = gr.index().positionAt(i);
                mappingQualities[i] = (byte) gr.index().mappingQualityAt(i);
                flagsArr[i] = gr.index().flagsAt(i);
                chromosomesList.add(gr.index().chromosomeAt(i));
            }

            out.add(new WrittenGenomicRun(
                gr.acquisitionMode(),
                gr.referenceUri(),
                gr.platform(),
                gr.sampleName(),
                positionsOut,
                mappingQualities,
                flagsArr,
                sequencesFlat,
                qualitiesFlat,
                offsets,
                lengths,
                cigars,
                readNames,
                mateChromosomes,
                matePositions,
                templateLengths,
                chromosomesList,
                Compression.ZLIB
            ));
        }
        return out;
    }

    // CIGAR reference span (number of bases consumed on the
    // reference). Ops that consume reference: M, D, N, =, X. Ops that
    // do NOT: I, S, H, P. Returns 0 for empty / "*" / non-parseable
    // CIGAR — caller falls back to position-only masking (M90.3
    // backward compatibility).
    static long cigarRefSpan(String cigar) {
        if (cigar == null || cigar.isEmpty() || cigar.equals("*")) {
            return 0L;
        }
        long total = 0L;
        long acc = 0L;
        boolean haveDigit = false;
        for (int i = 0; i < cigar.length(); i++) {
            char ch = cigar.charAt(i);
            if (ch >= '0' && ch <= '9') {
                acc = acc * 10 + (ch - '0');
                haveDigit = true;
            } else {
                if (!haveDigit) return 0L;  // malformed
                switch (ch) {
                    case 'M': case 'D': case 'N': case '=': case 'X':
                        total += acc;
                        break;
                    case 'I': case 'S': case 'H': case 'P':
                        // Don't consume reference bases.
                        break;
                    default:
                        // Unknown op — bail out and fall back to pos-only.
                        return 0L;
                }
                acc = 0L;
                haveDigit = false;
            }
        }
        if (haveDigit) return 0L;  // trailing digits with no op
        return total;
    }

    static boolean isSaav(String chemicalEntity) {
        if (chemicalEntity == null) return false;
        String upper = chemicalEntity.toUpperCase();
        return upper.contains("SAAV") || upper.contains("VARIANT");
    }

    static double quantileThreshold(double[] data, double quantile) {
        double[] sorted = data.clone();
        Arrays.sort(sorted);
        int idx = (int) (quantile * (sorted.length - 1));
        return sorted[Math.min(idx, sorted.length - 1)];
    }

    static Map<String, Double> loadDefaultPrevalence() {
        Map<String, Double> table = new LinkedHashMap<>();
        try (InputStream is = Anonymizer.class.getResourceAsStream(
                "/data/metabolite_prevalence.json")) {
            if (is == null) return table;
            String json = new String(is.readAllBytes(), StandardCharsets.UTF_8);
            // Simple JSON parse for {"key": value, ...}
            for (String line : json.split("\n")) {
                line = line.strip();
                if (line.startsWith("\"CHEBI:")) {
                    int colonIdx = line.indexOf(':');
                    // Find the key
                    int keyEnd = line.indexOf('"', 1);
                    String key = line.substring(1, keyEnd);
                    // Find the value
                    int valStart = line.indexOf(':', keyEnd) + 1;
                    String valStr = line.substring(valStart).replaceAll("[,\\s]", "");
                    try {
                        table.put(key, Double.parseDouble(valStr));
                    } catch (NumberFormatException ignored) {}
                }
            }
        } catch (IOException ignored) {}
        return table;
    }
}
