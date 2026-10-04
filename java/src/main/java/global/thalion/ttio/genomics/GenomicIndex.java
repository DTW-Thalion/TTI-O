/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.Enums.Precision;
import global.thalion.ttio.providers.CompoundField;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Per-read offsets, lengths, positions, mapping qualities, flags, and
 * chromosome strings for one {@link GenomicRun}. Held in memory as
 * parallel arrays; loaded eagerly when the run is opened.
 *
 * <p>Genomic analogue of {@link global.thalion.ttio.SpectrumIndex}.</p>
 *
 * <p><b>API status:</b> Stable ().</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOGenomicIndex}, Python {@code ttio.genomic_index.GenomicIndex}.</p>
 */
public final class GenomicIndex {

    static final int CHUNK_SIZE = 65536;
    static final int COMPRESSION_LEVEL = 6;

    private final long[]   offsets;          // uint64 — byte offset into sequence channel
    private final int[]    lengths;          // uint32 — read length in bases
    private final List<String> chromosomes;  // one per read
    private final long[]   positions;        // int64 — 0-based mapping position
    private final byte[]   mappingQualities; // uint8
    private final int[]    flags;            // uint32

    // PJ1: interned chromosome table, populated by readFrom() from the
    // on-disk uint16 id column (Task #82). Both null for in-memory
    // construction — region lookup then falls back to scanning the
    // chromosomes string list. When present, indicesForRegion resolves
    // the query name to its id ONCE and scans an int comparison instead
    // of an O(N) per-read String.equals. Mirrors Python
    // genomic_index.py:126-146.
    private final short[] chromosomeIds;                          // uint16 per read; null in-memory
    private final java.util.Map<String, Integer> chromosomeNameToId; // null when no ids

    public GenomicIndex(long[] offsets, int[] lengths,
                         List<String> chromosomes, long[] positions,
                         byte[] mappingQualities, int[] flags) {
        this(offsets, lengths, chromosomes, positions, mappingQualities,
             flags, null, null);
    }

    /** PJ1 internal constructor that also retains the interned
     *  chromosome ids (one uint16 per read) and a name→id map, enabling
     *  the vectorized region scan. Both id arguments are {@code null}
     *  for in-memory construction (string-fallback path). */
    private GenomicIndex(long[] offsets, int[] lengths,
                         List<String> chromosomes, long[] positions,
                         byte[] mappingQualities, int[] flags,
                         short[] chromosomeIds,
                         java.util.Map<String, Integer> chromosomeNameToId) {
        Objects.requireNonNull(offsets);
        Objects.requireNonNull(lengths);
        Objects.requireNonNull(chromosomes);
        Objects.requireNonNull(positions);
        Objects.requireNonNull(mappingQualities);
        Objects.requireNonNull(flags);
        if (lengths.length != offsets.length
                || chromosomes.size() != offsets.length
                || positions.length != offsets.length
                || mappingQualities.length != offsets.length
                || flags.length != offsets.length) {
            throw new IllegalArgumentException(
                "GenomicIndex column lengths must match");
        }
        this.offsets = offsets;
        this.lengths = lengths;
        this.chromosomes = List.copyOf(chromosomes);
        this.positions = positions;
        this.mappingQualities = mappingQualities;
        this.flags = flags;
        this.chromosomeIds = chromosomeIds;
        this.chromosomeNameToId = chromosomeNameToId;
    }

    /** This index with row {@code i} taken from row {@code order[i]},
     *  offsets recomputed from the reordered lengths (M103: a
     *  {@code blocks_v1_grouped} run presents its index in input order).
     *  The interned chromosome ids move with their rows. */
    public GenomicIndex permuted(int[] order) {
        int n = order.length;
        int[] lens = new int[n];
        long[] pos = new long[n];
        byte[] mapq = new byte[n];
        int[] fl = new int[n];
        List<String> chroms = new ArrayList<>(n);
        short[] ids = chromosomeIds == null ? null : new short[n];
        for (int i = 0; i < n; i++) {
            int s = order[i];
            lens[i] = lengths[s];
            pos[i] = positions[s];
            mapq[i] = mappingQualities[s];
            fl[i] = flags[s];
            chroms.add(chromosomes.get(s));
            if (ids != null) ids[i] = chromosomeIds[s];
        }
        return new GenomicIndex(offsetsFromLengths(lens), lens, chroms, pos, mapq, fl,
                                ids, chromosomeNameToId);
    }

    /** Number of reads. */
    public int count() { return offsets.length; }

    /** Byte offset of read {@code i} in the sequences/qualities channels. */
    public long offsetAt(int i) { return offsets[i]; }
    /** Length of read {@code i} in bases. */
    public int lengthAt(int i) { return lengths[i]; }
    /** 0-based mapping position of read {@code i}. */
    public long positionAt(int i) { return positions[i]; }
    /** Phred-scaled mapping quality of read {@code i}, unsigned (0–255). */
    public int mappingQualityAt(int i) { return mappingQualities[i] & 0xFF; }
    /** SAM flags of read {@code i}. */
    public int flagsAt(int i) { return flags[i]; }
    /** Reference sequence name of read {@code i}. */
    public String chromosomeAt(int i) { return chromosomes.get(i); }

    /** Read indices on {@code chromosome} with {@code start <= position < end}.
     *
     *  <p>PJ1: on disk-loaded indices (interned {@code chromosome_ids}
     *  present) the query chromosome is resolved to its uint16 id once,
     *  then the per-read id column is scanned with an int comparison —
     *  avoiding an O(N) {@code String.equals} per read. An absent
     *  chromosome resolves to no id and yields the empty list, matching
     *  the old never-matched behavior. Results are byte-identical to the
     *  string-scan fallback (same indices, same ascending order).</p>
     */
    public List<Integer> indicesForRegion(String chromosome, long start, long end) {
        if (chromosomeIds != null) {
            Integer tid = chromosomeNameToId.get(chromosome);
            if (tid == null) return List.of();
            int target = tid;
            int[] buf = new int[8];
            int n = 0;
            for (int i = 0; i < count(); i++) {
                if ((chromosomeIds[i] & 0xFFFF) == target
                        && positions[i] >= start && positions[i] < end) {
                    if (n == buf.length) buf = java.util.Arrays.copyOf(buf, n * 2);
                    buf[n++] = i;
                }
            }
            List<Integer> out = new ArrayList<>(n);
            for (int k = 0; k < n; k++) out.add(buf[k]);
            return out;
        }
        // Fallback: in-memory index without interned ids.
        List<Integer> out = new ArrayList<>();
        for (int i = 0; i < count(); i++) {
            if (chromosomes.get(i).equals(chromosome)
                    && positions[i] >= start && positions[i] < end) {
                out.add(i);
            }
        }
        return out;
    }

    /** Read indices where {@code (flags & 0x4) != 0}. */
    public List<Integer> indicesForUnmapped() { return indicesForFlag(0x4); }

    /** Read indices where {@code (flags & flagMask) != 0}. */
    public List<Integer> indicesForFlag(int flagMask) {
        int[] buf = new int[8];
        int n = 0;
        for (int i = 0; i < count(); i++) {
            if ((flags[i] & flagMask) != 0) {
                if (n == buf.length) buf = java.util.Arrays.copyOf(buf, n * 2);
                buf[n++] = i;
            }
        }
        List<Integer> out = new ArrayList<>(n);
        for (int k = 0; k < n; k++) out.add(buf[k]);
        return out;
    }

    // ── Disk I/O via the StorageGroup protocol ─────────────────────

    /** Write this index into {@code idxGroup} (typically created via
     *  {@code parent.createGroup("genomic_index")}). The mathematically
     *  redundant {@code offsets} column is omitted; readers synthesize
     *  it from {@code cumsum(lengths)}.
     */
    public void writeTo(StorageGroup idxGroup) {
        writeTo(idxGroup, null);
    }

    /** Write the index; {@code nameToId}, when given, is the shared
     *  chromosome id map of a {@code blocks_v1} run: existing entries
     *  keep their ids and new names are appended in place. */
    public void writeTo(StorageGroup idxGroup, java.util.Map<String, Integer> nameToId) {
        writeInts (idxGroup, "lengths",          Precision.UINT32, lengths);
        writeLongs(idxGroup, "positions",        Precision.INT64,  positions);
        writeBytes(idxGroup, "mapping_qualities", Precision.UINT8, mappingQualities);
        writeInts (idxGroup, "flags",            Precision.UINT32, flags);

        // L1 (Task #82 Phase B.1, 2026-05-01): chromosomes are stored
        // as `chromosome_ids` (uint16) + `chromosome_names` (compound)
        // instead of a single VL-string compound. The old layout cost
        // 42 MB of HDF5 fractal-heap overhead per chr22 .tio file
        // (one heap block per chunk × 432 chunks) just to repeat one
        // byte-string 1.77M times. Encounter-order id assignment —
        // first occurrence gets the next unused id; cross-language
        // byte-exact contract.
        if (nameToId == null) nameToId = new java.util.LinkedHashMap<>();
        short[] ids = new short[chromosomes.size()];
        for (int i = 0; i < chromosomes.size(); i++) {
            String name = chromosomes.get(i);
            Integer slot = nameToId.get(name);
            if (slot == null) {
                if (nameToId.size() > 65535) {
                    throw new IllegalStateException(
                        "genomic_index: > 65,535 unique chromosome names; "
                        + "uint16 chromosome_ids would overflow.");
                }
                slot = nameToId.size();
                nameToId.put(name, slot);
            }
            ids[i] = slot.shortValue();
        }
        // Write chromosome_ids as uint16 (Java has no unsigned short,
        // so we pass the ids array via the StorageDataset writeAll
        // path which interprets the raw 16-bit pattern).
        StorageDataset cids;
        try {
            cids = idxGroup.createDataset(
                "chromosome_ids", Precision.UINT16, ids.length,
                CHUNK_SIZE, Compression.ZLIB, COMPRESSION_LEVEL);
        } catch (UnsupportedOperationException e) {
            cids = idxGroup.createDataset(
                "chromosome_ids", Precision.UINT16, ids.length,
                0, Compression.NONE, 0);
        }
        try (StorageDataset closeMe = cids) {
            closeMe.writeAll(ids);
        }
        // Write chromosome_names as compound[(name, VL_str)].
        List<CompoundField> nameFields = List.of(
            new CompoundField("name", CompoundField.Kind.VL_STRING));
        List<Object[]> nameRows = new ArrayList<>(nameToId.size());
        for (String n : namesInIdOrder(nameToId)) nameRows.add(new Object[]{ n });
        try (StorageDataset ds = idxGroup.createCompoundDataset(
                "chromosome_names", nameFields, nameRows.size())) {
            ds.writeAll(nameRows);
        }
    }

    /** Names of {@code nameToId} sorted by id. */
    public static List<String> namesInIdOrder(java.util.Map<String, Integer> nameToId) {
        List<java.util.Map.Entry<String, Integer>> entries = new ArrayList<>(nameToId.entrySet());
        entries.sort(java.util.Map.Entry.comparingByValue());
        List<String> out = new ArrayList<>(entries.size());
        for (var e : entries) out.add(e.getKey());
        return out;
    }

    /** Read a {@link GenomicIndex} from an existing {@code genomic_index/}
     *  group (typically opened via {@code parent.openGroup("genomic_index")}).
     *
     *  <p>v1.10 #10 (2026-05-04): {@code offsets} is no longer stored
     *  on disk by default — it's mathematically derivable from
     *  {@code cumsum(lengths)}. Reader handles both layouts: pre-v1.10
     *  files have the column on disk (read directly); v1.10+ files
     *  synthesize it from lengths. Always uint64.</p>
     */
    @SuppressWarnings("unchecked")
    public static GenomicIndex readFrom(StorageGroup idxGroup) {
        int[]   lengths   = readInts (idxGroup, "lengths");
        long[]  offsets   = idxGroup.hasChild("offsets")
            ? readLongs(idxGroup, "offsets")
            : offsetsFromLengths(lengths);
        long[]  positions = readLongs(idxGroup, "positions");
        byte[]  mapqs     = readBytes(idxGroup, "mapping_qualities");
        int[]   flags     = readInts (idxGroup, "flags");

        // L1: read chromosome_ids (uint16) + chromosome_names
        // (compound) and materialise back to a List<String> so the
        // GenomicIndex API surface stays unchanged for callers.
        short[] ids;
        try (StorageDataset ds = idxGroup.openDataset("chromosome_ids")) {
            ids = (short[]) ds.readAll();
        }
        List<Object[]> nameRows;
        try (StorageDataset ds = idxGroup.openDataset("chromosome_names")) {
            nameRows = (List<Object[]>) ds.readAll();
        }
        List<String> nameTable = new ArrayList<>(nameRows.size());
        for (Object[] row : nameRows) {
            Object v = row[0];
            if (v instanceof byte[] b) {
                nameTable.add(new String(b, java.nio.charset.StandardCharsets.UTF_8));
            } else {
                nameTable.add(v == null ? "" : v.toString());
            }
        }
        List<String> chroms = new ArrayList<>(ids.length);
        for (short id : ids) {
            int idx = Short.toUnsignedInt(id);
            chroms.add(idx < nameTable.size() ? nameTable.get(idx) : "");
        }

        // PJ1: invert the id→name table (chromosome_names is written in
        // encounter-order id sequence — row k carries the name for id k,
        // a 1:1 correspondence) into a name→id map so indicesForRegion
        // can resolve the query chromosome to its uint16 id once and scan
        // the retained `ids` column with int comparisons.
        java.util.Map<String, Integer> nameToId =
            new java.util.HashMap<>(nameTable.size() * 2);
        for (int id = 0; id < nameTable.size(); id++) {
            nameToId.putIfAbsent(nameTable.get(id), id);
        }

        return new GenomicIndex(offsets, lengths, chroms, positions, mapqs,
                                flags, ids, nameToId);
    }

    // ── Typed-channel helpers (chunked + zlib, falling back to raw) ─

    private static void writeLongs(StorageGroup g, String name, Precision p, long[] data) {
        StorageDataset ds;
        try {
            ds = g.createDataset(name, p, data.length, CHUNK_SIZE,
                    Compression.ZLIB, COMPRESSION_LEVEL);
        } catch (UnsupportedOperationException e) {
            ds = g.createDataset(name, p, data.length, 0, Compression.NONE, 0);
        }
        try (StorageDataset closeMe = ds) { closeMe.writeAll(data); }
    }

    private static void writeInts(StorageGroup g, String name, Precision p, int[] data) {
        StorageDataset ds;
        try {
            ds = g.createDataset(name, p, data.length, CHUNK_SIZE,
                    Compression.ZLIB, COMPRESSION_LEVEL);
        } catch (UnsupportedOperationException e) {
            ds = g.createDataset(name, p, data.length, 0, Compression.NONE, 0);
        }
        try (StorageDataset closeMe = ds) { closeMe.writeAll(data); }
    }

    private static void writeBytes(StorageGroup g, String name, Precision p, byte[] data) {
        StorageDataset ds;
        try {
            ds = g.createDataset(name, p, data.length, CHUNK_SIZE,
                    Compression.ZLIB, COMPRESSION_LEVEL);
        } catch (UnsupportedOperationException e) {
            ds = g.createDataset(name, p, data.length, 0, Compression.NONE, 0);
        }
        try (StorageDataset closeMe = ds) { closeMe.writeAll(data); }
    }

    private static long[] readLongs(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (long[]) ds.readAll();
        }
    }

    private static int[] readInts(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (int[]) ds.readAll();
        }
    }

    private static byte[] readBytes(StorageGroup g, String name) {
        try (StorageDataset ds = g.openDataset(name)) {
            return (byte[]) ds.readAll();
        }
    }

    /**
     * v1.10 #10 helper: synthesize per-record byte offsets from a
     * lengths array. {@code offsets[i] = sum(lengths[0..i])}, returned
     * as a {@code long[]} (i.e. uint64) to avoid the >4 GB overflow
     * cliff on deep WGS even though the input lengths are uint32.
     *
     * <p>Empty input returns {@code new long[0]}.</p>
     */
    public static long[] offsetsFromLengths(int[] lengths) {
        long[] out = new long[lengths.length];
        if (lengths.length == 0) return out;
        out[0] = 0L;
        long acc = 0L;
        for (int i = 1; i < lengths.length; i++) {
            // Mask up the int to its uint32 value before adding to long.
            acc += ((long) lengths[i - 1]) & 0xFFFFFFFFL;
            out[i] = acc;
        }
        return out;
    }
}
