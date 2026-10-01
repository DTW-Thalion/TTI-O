/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** The {@code blocks/index} rows of a {@code blocks_v1} run as column
 *  arrays (format-spec 10.12.2). Public accessors serve the per-AU
 *  blocks_v1 walkers in the protection package (M99). */
public final class BlockTable {

    final long[] readStart;
    final int[] nReads;
    final long[] baseStart;
    final long[] nBases;
    /** channel -> byte offsets of each block's blob */
    final Map<String, long[]> off = new LinkedHashMap<>();
    /** channel -> byte lengths of each block's blob */
    final Map<String, long[]> len = new LinkedHashMap<>();
    /** channel -> codec id per block; {@code null} when the file has no
     *  codec columns */
    final Map<String, int[]> codec;
    /** The channels whose columns the index carries: the required five,
     *  plus {@code tags} when the run has the tags triple (M101). */
    List<String> channels = GenomicBlocks.REQUIRED_BLOCK_CHANNELS;

    private BlockTable(int n, boolean hasCodecs) {
        readStart = new long[n];
        nReads = new int[n];
        baseStart = new long[n];
        nBases = new long[n];
        codec = hasCodecs ? new LinkedHashMap<>() : null;
    }

    public long readStartAt(int b) { return readStart[b]; }
    public int nReadsAt(int b) { return nReads[b]; }
    public long baseStartAt(int b) { return baseStart[b]; }
    public long nBasesAt(int b) { return nBases[b]; }
    public long offsetOf(String channel, int b) { return off.get(channel)[b]; }
    public long lengthOf(String channel, int b) { return len.get(channel)[b]; }
    public boolean hasCodecs() { return codec != null; }
    /** The index's channel columns, in column order (M101). */
    public List<String> channels() { return channels; }
    /** True when the index carries the tags (off, len, codec) triple. */
    public boolean hasTagsColumns() { return channels.contains("tags"); }
    /** Codec id of a channel in a block; 0 without codec columns. */
    public int codecOf(String channel, int b) {
        return codec != null ? codec.get(channel)[b] : 0;
    }

    public static BlockTable read(StorageGroup runGroup) {
        List<Map<String, Object>> rows;
        try (StorageGroup blocks = runGroup.openGroup("blocks");
             StorageDataset ds = blocks.openDataset("index")) {
            rows = ds.readRows();
        }
        boolean hasCodecs = !rows.isEmpty()
            && rows.get(0).containsKey(GenomicBlocks.BLOCK_CHANNELS.get(0) + "_codec");
        BlockTable t = new BlockTable(rows.size(), hasCodecs);
        // An optional channel's columns are absent from runs that do
        // not carry it (M101: tags): every block then has an empty range.
        List<String> present = new java.util.ArrayList<>(GenomicBlocks.REQUIRED_BLOCK_CHANNELS);
        for (String ch : GenomicBlocks.OPTIONAL_BLOCK_CHANNELS) {
            if (!rows.isEmpty() && rows.get(0).containsKey(ch + "_off")) present.add(ch);
        }
        t.channels = List.copyOf(present);
        for (String ch : GenomicBlocks.BLOCK_CHANNELS) {
            t.off.put(ch, new long[rows.size()]);
            t.len.put(ch, new long[rows.size()]);
            if (hasCodecs) t.codec.put(ch, new int[rows.size()]);
        }
        for (int i = 0; i < rows.size(); i++) {
            Map<String, Object> r = rows.get(i);
            t.readStart[i] = num(r, "read_start");
            t.nReads[i] = (int) num(r, "n_reads");
            t.baseStart[i] = num(r, "base_start");
            t.nBases[i] = num(r, "n_bases");
            for (String ch : t.channels) {
                t.off.get(ch)[i] = num(r, ch + "_off");
                t.len.get(ch)[i] = num(r, ch + "_len");
                if (hasCodecs) t.codec.get(ch)[i] = (int) num(r, ch + "_codec");
            }
        }
        return t;
    }

    private static long num(Map<String, Object> row, String key) {
        Object v = row.get(key);
        if (v == null) throw new IllegalStateException("blocks/index: missing column " + key);
        return ((Number) v).longValue();
    }

    public int count() { return readStart.length; }

    long readCount() {
        int n = count();
        return n == 0 ? 0 : readStart[n - 1] + nReads[n - 1];
    }

    /** Index of the block holding read {@code i}. */
    int blockFor(long i) {
        int b = Arrays.binarySearch(readStart, i);
        if (b < 0) b = -b - 2;
        if (b < 0 || i >= readStart[b] + nReads[b]) {
            throw new IndexOutOfBoundsException(
                "read index " + i + " out of range [0, " + readCount() + ")");
        }
        return b;
    }
}
