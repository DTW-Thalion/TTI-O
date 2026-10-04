/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.Enums;
import global.thalion.ttio.SpectralDatasetGenomicWriter;
import global.thalion.ttio.providers.MemoryProvider;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;
import global.thalion.ttio.providers.StorageProvider;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Block encoder for the {@code blocks_v1} genomic layout (format-spec
 * 10.12). A block is a run made of a contiguous range of reads; its
 * channel blobs come from running the whole-channel writer against a
 * memory-provider group and harvesting each channel dataset's bytes and
 * {@code @compression}, so a block's blob is byte-identical to what a
 * v1.8 write of those reads alone produces. Python:
 * {@code ttio.genomic._blocks}.
 */
public final class GenomicBlocks {

    /** Blob channels every blocks_v1 run indexes, in block-index column
     *  order. */
    public static final List<String> REQUIRED_BLOCK_CHANNELS =
        List.of("sequences", "qualities", "read_names", "cigars", "mate_info");
    /** Channels whose index columns exist only in runs that carry them
     *  (M101: the SAM tags). Their (off, len, codec) columns follow the
     *  required channels' codec columns; a reader treats missing columns
     *  as an empty channel (format-spec 10.12.2). */
    public static final List<String> OPTIONAL_BLOCK_CHANNELS = List.of("tags");
    /** Every blob channel of a block. */
    public static final List<String> BLOCK_CHANNELS =
        List.of("sequences", "qualities", "read_names", "cigars", "mate_info", "tags");

    private GenomicBlocks() {}

    /** One block's encoded channels.
     *
     *  @param blobs      channel to blob bytes (empty array when absent)
     *  @param codecs     channel to codec id (0 when absent)
     *  @param extraAttrs channel to the dataset attributes other than
     *                    {@code compression}
     *  @param nReads     reads in the block
     *  @param nBases     bases in the block */
    public record BlockBlobs(Map<String, byte[]> blobs, Map<String, Integer> codecs,
                             Map<String, Map<String, Object>> extraAttrs,
                             int nReads, long nBases) {}

    /** Reads {@code [start, stop)} of {@code run} as a run of their own,
     *  offsets rebased to 0; run-level metadata shared. */
    public static WrittenGenomicRun sliceRun(WrittenGenomicRun run, int start, int stop) {
        long b0, b1;
        if (stop <= start) {
            b0 = b1 = 0;
        } else {
            b0 = run.offsets()[start];
            b1 = run.offsets()[stop - 1] + run.lengths()[stop - 1];
        }
        int n = Math.max(stop - start, 0);
        long[] offsets = new long[n];
        for (int i = 0; i < n; i++) offsets[i] = run.offsets()[start + i] - b0;
        return new WrittenGenomicRun(
            run.acquisitionMode(), run.referenceUri(), run.platform(), run.sampleName(),
            Arrays.copyOfRange(run.positions(), start, stop),
            Arrays.copyOfRange(run.mappingQualities(), start, stop),
            Arrays.copyOfRange(run.flags(), start, stop),
            Arrays.copyOfRange(run.sequences(), (int) b0, (int) b1),
            Arrays.copyOfRange(run.qualities(), (int) b0, (int) b1),
            offsets,
            Arrays.copyOfRange(run.lengths(), start, stop),
            new ArrayList<>(run.cigars().subList(start, stop)),
            new ArrayList<>(run.readNames().subList(start, stop)),
            new ArrayList<>(run.mateChromosomes().subList(start, stop)),
            Arrays.copyOfRange(run.matePositions(), start, stop),
            Arrays.copyOfRange(run.templateLengths(), start, stop),
            new ArrayList<>(run.chromosomes().subList(start, stop)),
            run.signalCompression(), run.signalCodecOverrides(), List.of(),
            run.embedReference(), run.referenceChromSeqs(), run.externalReferencePath(),
            null, run.optDisableQualitiesV5(), run.optLegacyWholeChannel(),
            run.readRole(), run.refDiffSliceBytes(),
            run.tags() == null ? null : new ArrayList<>(run.tags().subList(start, stop)));
    }

    /** The inverse of {@link #sliceRun} for consecutive parts. */
    public static WrittenGenomicRun concatRuns(List<WrittenGenomicRun> parts) {
        if (parts.size() == 1) return parts.get(0);
        WrittenGenomicRun first = parts.get(0);
        int n = 0;
        long bases = 0;
        for (WrittenGenomicRun p : parts) { n += p.readCount(); bases += p.sequences().length; }
        long[] positions = new long[n];
        byte[] mapqs = new byte[n];
        int[] flags = new int[n];
        int[] lengths = new int[n];
        long[] matePos = new long[n];
        int[] tlens = new int[n];
        byte[] seqs = new byte[(int) bases];
        byte[] quals = new byte[(int) bases];
        List<String> cigars = new ArrayList<>(n), names = new ArrayList<>(n),
                     mateChroms = new ArrayList<>(n), chroms = new ArrayList<>(n);
        boolean anyTags = false;
        for (WrittenGenomicRun p : parts) anyTags |= p.tags() != null;
        List<String> tags = anyTags ? new ArrayList<>(n) : null;
        int r = 0, b = 0;
        for (WrittenGenomicRun p : parts) {
            if (tags != null) {
                if (p.tags() != null) tags.addAll(p.tags());
                else for (int i = 0; i < p.readCount(); i++) tags.add("");
            }
            int pn = p.readCount();
            System.arraycopy(p.positions(), 0, positions, r, pn);
            System.arraycopy(p.mappingQualities(), 0, mapqs, r, pn);
            System.arraycopy(p.flags(), 0, flags, r, pn);
            System.arraycopy(p.lengths(), 0, lengths, r, pn);
            System.arraycopy(p.matePositions(), 0, matePos, r, pn);
            System.arraycopy(p.templateLengths(), 0, tlens, r, pn);
            System.arraycopy(p.sequences(), 0, seqs, b, p.sequences().length);
            System.arraycopy(p.qualities(), 0, quals, b, p.qualities().length);
            cigars.addAll(p.cigars());
            names.addAll(p.readNames());
            mateChroms.addAll(p.mateChromosomes());
            chroms.addAll(p.chromosomes());
            r += pn;
            b += p.sequences().length;
        }
        return new WrittenGenomicRun(
            first.acquisitionMode(), first.referenceUri(), first.platform(), first.sampleName(),
            positions, mapqs, flags, seqs, quals,
            GenomicIndex.offsetsFromLengths(lengths), lengths,
            cigars, names, mateChroms, matePos, tlens, chroms,
            first.signalCompression(), first.signalCodecOverrides(), List.of(),
            first.embedReference(), first.referenceChromSeqs(), first.externalReferencePath(),
            null, first.optDisableQualitiesV5(), first.optLegacyWholeChannel(),
            first.readRole(), first.refDiffSliceBytes(), tags);
    }

    /** The reads of {@code run} in the order {@code order} (read
     *  {@code order[j]} at row {@code j}), offsets recomputed (M103).
     *  Run-level metadata is shared; provenance is dropped. */
    public static WrittenGenomicRun takeRun(WrittenGenomicRun run, int[] order) {
        int n = order.length;
        long[] positions = new long[n];
        byte[] mapqs = new byte[n];
        int[] flags = new int[n];
        int[] lengths = new int[n];
        long[] matePos = new long[n];
        int[] tlens = new int[n];
        List<String> cigars = new ArrayList<>(n), names = new ArrayList<>(n),
                     mateChroms = new ArrayList<>(n), chroms = new ArrayList<>(n);
        List<String> tags = run.tags() == null ? null : new ArrayList<>(n);
        long total = 0;
        for (int j = 0; j < n; j++) {
            int i = order[j];
            positions[j] = run.positions()[i];
            mapqs[j] = run.mappingQualities()[i];
            flags[j] = run.flags()[i];
            lengths[j] = run.lengths()[i];
            matePos[j] = run.matePositions()[i];
            tlens[j] = run.templateLengths()[i];
            cigars.add(run.cigars().get(i));
            names.add(run.readNames().get(i));
            mateChroms.add(run.mateChromosomes().get(i));
            chroms.add(run.chromosomes().get(i));
            if (tags != null) tags.add(run.tags().get(i));
            total += lengths[j];
        }
        boolean withQuals = run.qualities().length > 0;
        byte[] seqs = new byte[(int) total];
        byte[] quals = withQuals ? new byte[(int) total] : new byte[0];
        int o = 0;
        for (int j = 0; j < n; j++) {
            int src = (int) run.offsets()[order[j]];
            int len = lengths[j];
            System.arraycopy(run.sequences(), src, seqs, o, len);
            if (withQuals) System.arraycopy(run.qualities(), src, quals, o, len);
            o += len;
        }
        return new WrittenGenomicRun(
            run.acquisitionMode(), run.referenceUri(), run.platform(), run.sampleName(),
            positions, mapqs, flags, seqs, quals,
            GenomicIndex.offsetsFromLengths(lengths), lengths,
            cigars, names, mateChroms, matePos, tlens, chroms,
            run.signalCompression(), run.signalCodecOverrides(), List.of(),
            run.embedReference(), run.referenceChromSeqs(), run.externalReferencePath(),
            null, run.optDisableQualitiesV5(), run.optLegacyWholeChannel(),
            run.readRole(), run.refDiffSliceBytes(), tags, run.optGroupReads());
    }

    /** The read-grouping permutation of {@code run} (M103): reads are
     *  grouped by sequence within each chromosome label, labels in
     *  first-seen order, so a block still never spans two labels. Entry
     *  {@code j} is the input index of the read stored at row {@code j}. */
    public static int[] groupOrder(WrittenGenomicRun run) {
        int n = run.readCount();
        if (n == 0) return new int[0];
        Map<String, List<Integer>> labels = new LinkedHashMap<>();
        List<String> chroms = run.chromosomes();
        for (int i = 0; i < n; i++) labels.computeIfAbsent(chroms.get(i), k -> new ArrayList<>()).add(i);
        int[] out = new int[n];
        int w = 0;
        for (List<Integer> members : labels.values()) {
            int[] m = new int[members.size()];
            for (int k = 0; k < m.length; k++) m[k] = members.get(k);
            WrittenGenomicRun sub = labels.size() == 1 ? run : takeRun(run, m);
            int[] local = global.thalion.ttio.codecs.SeqGroup.group(
                sub.sequences(), sub.lengths(), sub.readNames());
            for (int k : local) out[w++] = m[k];
        }
        return out;
    }

    /** Encode one block's channels through the whole-channel writer. The
     *  forced codecs of format-spec 10.12.3 apply: cigars RANS_ORDER0,
     *  qualities FQZCOMP_NX16_Z (RANS_ORDER0 when the block holds a
     *  zero-length read), sequences RANS_ORDER1 without a reference. */
    public static BlockBlobs encodeBlock(WrittenGenomicRun block, GenomicWriteContext ctx) {
        Map<String, Enums.Compression> ov = new LinkedHashMap<>(block.signalCodecOverrides());
        ov.putIfAbsent("cigars", Enums.Compression.RANS_ORDER0);
        if (!ov.containsKey("qualities")) {
            boolean zero = false;
            for (int l : block.lengths()) if (l == 0) { zero = true; break; }
            ov.put("qualities", zero ? Enums.Compression.RANS_ORDER0
                                     : Enums.Compression.FQZCOMP_NX16_Z);
        }
        // Sequences without a reference: SEQ_CM (M103), or RANS_ORDER1 when
        // libttio_rans_jni is not loaded.
        if (!ov.containsKey("sequences") && block.referenceChromSeqs() == null) {
            ov.put("sequences", global.thalion.ttio.codecs.SeqCm.isAvailable()
                ? Enums.Compression.SEQ_CM : Enums.Compression.RANS_ORDER1);
        }
        WrittenGenomicRun b = block.withSignalCodecOverrides(ov).withProvenance(List.of());
        String url = "memory://ttio-block-encode-" + System.identityHashCode(block)
                   + "-" + System.nanoTime();
        try (StorageProvider mem = new MemoryProvider().open(url, StorageProvider.Mode.CREATE)) {
            StorageGroup root = mem.rootGroup();
            SpectralDatasetGenomicWriter.writeGenomicRunSubtree(root, "b", b, ctx);
            StorageGroup sc = root.openGroup("b").openGroup("signal_channels");
            Map<String, byte[]> blobs = new LinkedHashMap<>();
            Map<String, Integer> codecs = new LinkedHashMap<>();
            Map<String, Map<String, Object>> extra = new LinkedHashMap<>();
            for (String ch : BLOCK_CHANNELS) {
                StorageDataset ds = null;
                if (ch.equals("sequences")) {
                    StorageGroup g = tryGroup(sc, "sequences");
                    if (g != null && g.hasChild("refdiff_v2")) {
                        ds = g.openDataset("refdiff_v2");
                    } else if (g == null && sc.hasChild("sequences")) {
                        ds = sc.openDataset("sequences");
                    }
                } else if (ch.equals("mate_info")) {
                    StorageGroup g = tryGroup(sc, "mate_info");
                    if (g != null && g.hasChild("inline_v2")) ds = g.openDataset("inline_v2");
                } else if (sc.hasChild(ch) && tryGroup(sc, ch) == null) {
                    ds = sc.openDataset(ch);
                }
                if (ds == null) {
                    blobs.put(ch, new byte[0]);
                    codecs.put(ch, 0);
                    extra.put(ch, Map.of());
                    continue;
                }
                Object raw = ds.readAll();
                blobs.put(ch, raw == null ? new byte[0] : (byte[]) raw);
                int codec = 0;
                Map<String, Object> attrs = new LinkedHashMap<>();
                for (String k : ds.attributeNames()) {
                    Object v = ds.getAttribute(k);
                    if (k.equals("compression")) {
                        codec = ((Number) v).intValue();
                    } else {
                        attrs.put(k, v);
                    }
                }
                codecs.put(ch, codec);
                extra.put(ch, attrs);
            }
            long nBases = 0;
            for (int l : block.lengths()) nBases += l;
            return new BlockBlobs(blobs, codecs, extra, block.readCount(), nBases);
        } finally {
            MemoryProvider.discardStore(url);
        }
    }

    /** {@code parent}'s child {@code name} as a group, or {@code null}
     *  when it is a dataset or absent. */
    public static StorageGroup tryGroup(StorageGroup parent, String name) {
        if (!parent.hasChild(name)) return null;
        try {
            return parent.openGroup(name);
        } catch (RuntimeException e) {
            return null;
        }
    }
}
