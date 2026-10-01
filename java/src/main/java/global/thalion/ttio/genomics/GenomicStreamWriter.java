/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.Enums.Precision;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.SpectralDatasetGenomicWriter;
import global.thalion.ttio.providers.CompoundField;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Writes one genomic run as {@code blocks_v1} (format-spec 10.12) with
 * bounded memory. Reads are buffered until a block is full
 * ({@code blockReads} reads or {@code blockBytes} sequence bytes,
 * whichever first; a block never spans two chromosomes), encoded through
 * {@link GenomicBlocks#encodeBlock}, and appended to extendable
 * per-channel datasets; {@code blocks/index} records where each block's
 * blob lives. Python: {@code ttio.genomic.GenomicStreamWriter}.
 */
public final class GenomicStreamWriter implements AutoCloseable {

    public static final String LAYOUT = "blocks_v1";
    public static final int DEFAULT_BLOCK_READS = 1_000_000;
    /* A block is the unit of encode concurrency and of residency, so a
 * big one costs both: the pool cannot start a block that is still
 * filling, and every block in flight is held whole. At 256 MiB a
 * 2x250 import of 10.6 M reads was 11 blocks, which left 7 workers
 * at 3.05x of one thread and 10.4 GiB resident. At 64 MiB the same
 * import is 43.0 s against 70.3 s and 4.50 GiB against 10.38 GiB.
 * The cost is context the codecs no longer share: 0.16% on 2x250,
 * 0.05% on 150 bp FASTQ, and 0.46% on HiFi, where reads are long
 * enough that a block holds few of them. */
    public static final long DEFAULT_BLOCK_BYTES = 64L << 20;
    /** Chunk of the unfiltered channel datasets. */
    public static final int CHANNEL_CHUNK = 256 << 10;

    /** Block index schema, in the column order of format-spec 10.12.2. */
    public static final List<CompoundField> INDEX_FIELDS;
    static {
        List<CompoundField> f = new ArrayList<>();
        f.add(new CompoundField("read_start", CompoundField.Kind.UINT64));
        f.add(new CompoundField("n_reads", CompoundField.Kind.UINT32));
        f.add(new CompoundField("base_start", CompoundField.Kind.UINT64));
        f.add(new CompoundField("n_bases", CompoundField.Kind.UINT64));
        for (String ch : GenomicBlocks.REQUIRED_BLOCK_CHANNELS) {
            f.add(new CompoundField(ch + "_off", CompoundField.Kind.UINT64));
            f.add(new CompoundField(ch + "_len", CompoundField.Kind.UINT64));
        }
        for (String ch : GenomicBlocks.REQUIRED_BLOCK_CHANNELS) {
            f.add(new CompoundField(ch + "_codec", CompoundField.Kind.UINT32));
        }
        INDEX_FIELDS = Collections.unmodifiableList(f);
        List<CompoundField> t = new ArrayList<>(f);
        t.add(new CompoundField("tags_off", CompoundField.Kind.UINT64));
        t.add(new CompoundField("tags_len", CompoundField.Kind.UINT64));
        t.add(new CompoundField("tags_codec", CompoundField.Kind.UINT32));
        INDEX_FIELDS_WITH_TAGS = Collections.unmodifiableList(t);
    }

    /** The schema of a run that carries the tags channel (M101): the
     *  tags (off, len, codec) triple follows the required columns. */
    public static final List<CompoundField> INDEX_FIELDS_WITH_TAGS;

    /** The index schema for a run with or without the tags triple. */
    public static List<CompoundField> indexFields(boolean withTags) {
        return withTags ? INDEX_FIELDS_WITH_TAGS : INDEX_FIELDS;
    }

    /** Column of {@code ch}'s offset in a positional index row; its
     *  length follows, and {@link #codecColumn} gives its codec. */
    public static int offsetColumn(String ch) {
        int ci = GenomicBlocks.REQUIRED_BLOCK_CHANNELS.indexOf(ch);
        if (ci >= 0) return 4 + 2 * ci;
        if (ch.equals("tags")) return INDEX_FIELDS.size();
        throw new IllegalArgumentException("not a block channel: " + ch);
    }

    /** Column of {@code ch}'s codec in a positional index row. */
    public static int codecColumn(String ch) {
        int ci = GenomicBlocks.REQUIRED_BLOCK_CHANNELS.indexOf(ch);
        if (ci >= 0) return 4 + 2 * GenomicBlocks.REQUIRED_BLOCK_CHANNELS.size() + ci;
        if (ch.equals("tags")) return INDEX_FIELDS.size() + 2;
        throw new IllegalArgumentException("not a block channel: " + ch);
    }

    /** Per-read index arrays: name, precision. */
    private static final String[][] INDEX_ARRAYS = {
        {"lengths", "UINT32"}, {"positions", "INT64"}, {"mapping_qualities", "UINT8"},
        {"flags", "UINT32"}, {"chromosome_ids", "UINT16"},
    };

    /** Run-level options of a streamed run. */
    public record Options(AcquisitionMode acquisitionMode, String referenceUri, String platform,
                          String sampleName, Map<String, byte[]> referenceChromSeqs,
                          boolean embedReference, int blockReads, long blockBytes,
                          boolean optDisableQualitiesV5,
                          Map<String, Compression> signalCodecOverrides,
                          Compression signalCompression, boolean optLegacyWholeChannel,
                          List<ProvenanceRecord> provenanceRecords,
                          /** M97 — persisted as the {@code @read_role} run
                           *  attribute when non-null; see
                           *  {@link WrittenGenomicRun#readRole()}. */
                          String readRole,
                          /** M97 — REF_DIFF_V2 slice byte budget; 0 = the
                           *  fixed 10,000-read rule. */
                          long refDiffSliceBytes) {
        public Options {
            if (blockReads < 1 || blockBytes < 1) {
                throw new IllegalArgumentException("blockReads and blockBytes must be >= 1");
            }
            signalCodecOverrides = signalCodecOverrides == null ? Map.of() : Map.copyOf(signalCodecOverrides);
            provenanceRecords = provenanceRecords == null ? List.of() : List.copyOf(provenanceRecords);
            if (signalCompression == null) signalCompression = Compression.ZLIB;
        }

        /** Pre-M97 signature (13 components); no read role, the fixed
         *  REF_DIFF_V2 slice rule. */
        public Options(AcquisitionMode acquisitionMode, String referenceUri, String platform,
                       String sampleName, Map<String, byte[]> referenceChromSeqs,
                       boolean embedReference, int blockReads, long blockBytes,
                       boolean optDisableQualitiesV5,
                       Map<String, Compression> signalCodecOverrides,
                       Compression signalCompression, boolean optLegacyWholeChannel,
                       List<ProvenanceRecord> provenanceRecords) {
            this(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, blockReads, blockBytes,
                optDisableQualitiesV5, signalCodecOverrides, signalCompression,
                optLegacyWholeChannel, provenanceRecords, null, 0L);
        }

        /** The run-level metadata of {@code run}, default block policy. */
        public static Options fromRun(WrittenGenomicRun run) {
            return new Options(run.acquisitionMode(), run.referenceUri(), run.platform(),
                run.sampleName(), run.referenceChromSeqs(), run.embedReference(),
                DEFAULT_BLOCK_READS, DEFAULT_BLOCK_BYTES, run.optDisableQualitiesV5(),
                run.signalCodecOverrides(), run.signalCompression(),
                run.optLegacyWholeChannel(), run.provenanceRecords(),
                run.readRole(), run.refDiffSliceBytes());
        }

        public Options withBlockPolicy(int reads, long bytes) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, reads, bytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, optLegacyWholeChannel, provenanceRecords,
                readRole, refDiffSliceBytes);
        }

        public Options withLegacy(boolean legacy) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, blockReads, blockBytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, legacy, provenanceRecords,
                readRole, refDiffSliceBytes);
        }

        public Options withReference(Map<String, byte[]> reference, boolean embed) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                reference, embed, blockReads, blockBytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, optLegacyWholeChannel, provenanceRecords,
                readRole, refDiffSliceBytes);
        }

        public Options withProvenance(List<ProvenanceRecord> records) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, blockReads, blockBytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, optLegacyWholeChannel, records,
                readRole, refDiffSliceBytes);
        }

        /** Same options with the given {@code @read_role} value (M97). */
        public Options withReadRole(String role) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, blockReads, blockBytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, optLegacyWholeChannel, provenanceRecords,
                role, refDiffSliceBytes);
        }

        /** Same options with the given REF_DIFF_V2 slice byte budget (M97). */
        public Options withRefDiffSliceBytes(long sliceBytes) {
            return new Options(acquisitionMode, referenceUri, platform, sampleName,
                referenceChromSeqs, embedReference, blockReads, blockBytes, optDisableQualitiesV5,
                signalCodecOverrides, signalCompression, optLegacyWholeChannel, provenanceRecords,
                readRole, sliceBytes);
        }
    }

    private final StorageGroup study;
    private final String name;
    private final Options opt;
    private final List<WrittenGenomicRun> pending = new ArrayList<>();
    private int pendingReads;
    private long pendingBytes;
    private String pendingChrom;
    private final Map<String, Integer> chromMap = new LinkedHashMap<>();
    private byte[] referenceMd5;
    private long readCount;
    private long baseCount;
    private int blockCount;
    private StorageGroup rg;
    private final Map<String, StorageDataset> channelDs = new LinkedHashMap<>();
    private final Map<String, StorageDataset> idxDs = new LinkedHashMap<>();
    private StorageDataset indexDs;
    /** The index carries the tags columns (M101). */
    private boolean withTags;
    private boolean embedded;
    private boolean closed;
    private final List<WrittenGenomicRun> legacyParts = new ArrayList<>();
    private final int threads;
    private final global.thalion.ttio.Threads.PoolScope scope;
    private final java.util.ArrayDeque<InFlight> inflight = new java.util.ArrayDeque<>();
    private final long memoryBudgetBytes;
    private long inflightBytes;
    private long maxInFlightBytesObserved;
    private record InFlight(WrittenGenomicRun block,
                            java.util.concurrent.Future<GenomicBlocks.BlockBlobs> blobs,
                            long estimatedBytes) {}

    /** Per-run sticky qualities strategy: block 0 auto-tunes, the winner
     *  is pinned for the rest of the run. -1 = not yet pinned. */
    private int qualStrategyHint = pinnedHint();
    private final boolean qualExhaustive = exhaustiveTune();

    private static boolean exhaustiveTune() {
        String raw = System.getProperty("ttio.m94z.exhaustive");
        if (raw == null || raw.isBlank()) raw = System.getenv("TTIO_M94Z_EXHAUSTIVE");
        return "1".equals(raw);
    }

    /** A positive {@code -Dttio.m94z.hint} or {@code TTIO_M94Z_HINT}
     *  pins the strategy before block 0 instead, skipping the tune. V6
     *  is reachable no other way, so measuring it through the writer
     *  needs it. The exhaustive flag wins. */
    private static int pinnedHint() {
        if (exhaustiveTune()) return -1;
        String raw = System.getProperty("ttio.m94z.hint");
        if (raw == null || raw.isBlank()) raw = System.getenv("TTIO_M94Z_HINT");
        if (raw != null && !raw.isBlank()) {
            try {
                int v = Integer.parseInt(raw.trim());
                if (v > 0) return v;
            } catch (NumberFormatException ignored) {
                // fall through to the tune
            }
        }
        return -1;
    }

    int qualStrategyHintForTests() { return qualStrategyHint; }

    /** Append reads to run {@code runName} of the {@code /study} group
     *  {@code studyGroup} (the writer creates {@code genomic_runs} when
     *  absent and maintains {@code @_run_names}). */
    public GenomicStreamWriter(StorageGroup studyGroup, String runName, Options options) {
        this(studyGroup, runName, options, global.thalion.ttio.Threads.resolve(null));
    }

    /** With {@code threads} > 1 completed blocks are encoded on a pool and
     *  written in order by the caller's thread; at most {@code threads + 1}
     *  blocks are pending or encoding, so memory is about that many block
     *  working sets. The file is byte for byte the one thread's. */
    public GenomicStreamWriter(StorageGroup studyGroup, String runName, Options options, int threads) {
        this(studyGroup, runName, options, threads, 0L);
    }

    /** As above with an explicit pipeline byte budget (0 = the
     *  {@link global.thalion.ttio.Threads#resolveMemoryBudget} default).
     *  The writer stalls block submission while its in-flight estimate
     *  exceeds half of it; the count window stays the upper bound. */
    public GenomicStreamWriter(StorageGroup studyGroup, String runName, Options options,
                               int threads, long memoryBudgetBytes) {
        this.study = studyGroup;
        this.name = runName;
        // M97: QUALITY_BINNED is a fixed Illumina-8 bin table; reject
        // it on long-read platforms before anything is written.
        if (options.signalCodecOverrides().get("qualities")
                == Compression.QUALITY_BINNED
                && !global.thalion.ttio.codecs.Quality
                    .binnedAllowedForPlatform(options.platform())) {
            throw new IllegalArgumentException(
                "signalCodecOverrides['qualities']: QUALITY_BINNED "
                + "applies the fixed Illumina-8 bin table, which does "
                + "not fit the quality distribution of platform '"
                + options.platform() + "' (M97).");
        }
        this.opt = options;
        this.threads = Math.max(1, threads);
        this.scope = global.thalion.ttio.Threads.pool(this.threads);
        this.memoryBudgetBytes = global.thalion.ttio.Threads.resolveMemoryBudget(
            memoryBudgetBytes > 0 ? memoryBudgetBytes : null, this.threads, options.blockBytes());
    }

    public long readCount() { return readCount; }
    public int threads() { return threads; }
    public long memoryBudgetBytes() { return memoryBudgetBytes; }
    /** High-water mark of the in-flight byte estimate (tests). */
    public long maxInFlightBytesObserved() { return maxInFlightBytesObserved; }

    private static long estimateBlockBytes(WrittenGenomicRun b) {
        long raw = (long) b.sequences().length + b.qualities().length + (long) b.offsets().length * 24L;
        return raw * 2L;
    }

    /** Assign ids for every chromosome name the block introduces, in the
     *  order the block encoder assigns them (own names in read order, "*"
     *  included, then mate names that are not "", "*" or "="), so the
     *  encoder only reads the map and blocks can encode concurrently. */
    static void registerBlockChromosomes(WrittenGenomicRun block, Map<String, Integer> map) {
        for (String n : block.chromosomes()) map.putIfAbsent(n, map.size());
        for (String n : block.mateChromosomes()) {
            if (n != null && !n.isEmpty() && !"*".equals(n) && !"=".equals(n)) map.putIfAbsent(n, map.size());
        }
    }
    public int blockCount() { return blockCount; }

    /** Append one read. */
    public void append(AlignedRead read) {
        appendBatch(singleReadRun(read));
    }

    /** Append the reads of {@code batch}; its run-level metadata is
     *  ignored, the writer's options apply. */
    public void appendBatch(WrittenGenomicRun batch) {
        if (closed) throw new IllegalStateException("writer is closed");
        int n = batch.readCount();
        if (n == 0) return;
        if (opt.optLegacyWholeChannel()) {
            legacyParts.add(batch);
            return;
        }
        List<String> chroms = batch.chromosomes();
        int start = 0;
        while (start < n) {
            String chrom = chroms.get(start);
            int segEnd = start + 1;
            while (segEnd < n && chroms.get(segEnd).equals(chrom)) segEnd++;
            if (!pending.isEmpty() && !chrom.equals(pendingChrom)) flush();
            pendingChrom = chrom;
            while (start < segEnd) {
                int roomReads = opt.blockReads() - pendingReads;
                long roomBytes = opt.blockBytes() - pendingBytes;
                int stop = Math.min(segEnd, start + Math.max(roomReads, 1));
                long cum = 0;
                int fit = 0;
                for (int i = start; i < stop; i++) {
                    cum += batch.lengths()[i];
                    if (cum <= roomBytes) fit++; else break;
                }
                if (fit < stop - start) stop = start + Math.max(fit, 1);
                WrittenGenomicRun part = (start == 0 && stop == n) ? batch
                    : GenomicBlocks.sliceRun(batch, start, stop);
                pending.add(part);
                pendingReads += stop - start;
                for (int l : part.lengths()) pendingBytes += l;
                if (pendingReads >= opt.blockReads() || pendingBytes >= opt.blockBytes()) flush();
                start = stop;
            }
        }
    }

    /** Encode and write the pending reads as one block. */
    public void flush() {
        if (opt.optLegacyWholeChannel() || pending.isEmpty()) return;
        WrittenGenomicRun block = GenomicBlocks.concatRuns(pending);
        pending.clear();
        pendingReads = 0;
        pendingBytes = 0;
        if (referenceMd5 == null && opt.referenceChromSeqs() != null) {
            referenceMd5 = SpectralDatasetGenomicWriter.referenceMd5ForRun(applyMeta(block));
        }
        block = applyMeta(block);
        if (!embedded && opt.embedReference()) {
            SpectralDatasetGenomicWriter.embedReferencesForRuns(study, List.of(block));
            embedded = true;
        }
        registerBlockChromosomes(block, chromMap);
        if (!qualExhaustive && qualStrategyHint == -1
                && (blockCount > 0 || !inflight.isEmpty())) {
            // Block-0 gate: the pin comes from the earliest written M94Z
            // block; hold later submissions until it is known so the
            // choice is timing-independent.
            drain(0);
        }
        GenomicWriteContext ctx =
            new GenomicWriteContext(chromMap, referenceMd5, qualStrategyHint);
        if (scope.executor() == null) {
            writeEncoded(block, GenomicBlocks.encodeBlock(block, ctx));
            return;
        }
        drain(threads);   // window: threads in the pool plus this one
        // The worker reads the map while later flushes mutate it: give each
        // block a snapshot (registration above fixed every id it needs).
        GenomicWriteContext bctx =
            new GenomicWriteContext(new java.util.LinkedHashMap<>(chromMap),
                                    referenceMd5, qualStrategyHint);
        drainToBudget();
        WrittenGenomicRun fb = block;
        long est = estimateBlockBytes(block);
        // Segment threads follow the blocks actually in flight, not the
        // pool's size; see Threads.resolveV6SegmentThreads.
        global.thalion.ttio.Threads.applyV6SegmentThreadsForBlocksInFlight(
            inflight.size() + 1);
        inflight.add(new InFlight(block, scope.executor().submit(() -> GenomicBlocks.encodeBlock(fb, bctx)), est));
        inflightBytes += est;
        if (inflightBytes > maxInFlightBytesObserved) maxInFlightBytesObserved = inflightBytes;
    }

    /** Write completed blocks in sequence order; wait on the oldest until
     *  at most {@code blockUntil} remain in flight. */
    /** Drain until the in-flight estimate fits the writer's half of the
     *  byte budget. */
    private void drainToBudget() {
        long half = memoryBudgetBytes / 2L;
        while (!inflight.isEmpty() && inflightBytes > half) {
            drain(inflight.size() - 1);
        }
    }

    private void drain(int blockUntil) {
        while (!inflight.isEmpty() && (inflight.size() > blockUntil || inflight.peekFirst().blobs().isDone())) {
            InFlight f = inflight.pollFirst();
            inflightBytes -= Math.min(f.estimatedBytes(), inflightBytes);
            try {
                writeEncoded(f.block(), f.blobs().get());
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException(e);
            } catch (java.util.concurrent.ExecutionException e) {
                Throwable c = e.getCause();
                throw c instanceof RuntimeException r ? r : new IllegalStateException(c);
            }
        }
    }

    private void writeEncoded(WrittenGenomicRun block, GenomicBlocks.BlockBlobs blobs) {
        ensureLayout(blobs);
        byte[] tagBlob = blobs.blobs().get("tags");
        if (tagBlob != null && tagBlob.length > 0 && !withTags) addTagsColumns();
        // Rows are filled by position: the required channels keep
        // their columns and the tags triple trails them (M101).
        Object[] row = new Object[indexFields(withTags).size()];
        row[0] = readCount;
        row[1] = blobs.nReads();
        row[2] = baseCount;
        row[3] = blobs.nBases();
        for (String ch : GenomicBlocks.BLOCK_CHANNELS) {
            if (ch.equals("tags") && !withTags) continue;
            int col = offsetColumn(ch);
            byte[] data = blobs.blobs().get(ch);
            if (data == null) data = new byte[0];
            int codec = blobs.codecs().getOrDefault(ch, 0);
            row[codecColumn(ch)] = codec;
            StorageDataset ds = channelDs.get(ch);
            if (ds == null) {
                if (data.length > 0) {
                    ds = createChannel(ch, blobs);
                } else {
                    row[col] = 0L;
                    row[col + 1] = 0L;
                    continue;
                }
            }
            row[col] = ds.length();
            row[col + 1] = (long) data.length;
            if (data.length > 0) ds.append(data);
        }
        indexDs.append(Collections.singletonList(row));
        appendIndexArrays(block);
        readCount += blobs.nReads();
        baseCount += blobs.nBases();
        blockCount++;
        Integer qc = blobs.codecs().get("qualities");
        if (qualStrategyHint == -1 && !qualExhaustive && qc != null
                && qc == global.thalion.ttio.Enums.Compression.FQZCOMP_NX16_Z.ordinal()) {
            int strat = global.thalion.ttio.codecs.FqzcompNx16Z
                .streamStrategy(blobs.blobs().get("qualities"));
            qualStrategyHint = strat == 4
                ? global.thalion.ttio.codecs.FqzcompNx16Z.HINT_V4_AUTO : strat;
        }
        rg.setAttribute("read_count", readCount);
        rg.setAttribute("base_count", baseCount);
    }

    @Override
    public void close() {
        if (closed) return;
        closed = true;
        if (opt.optLegacyWholeChannel()) {
            if (!legacyParts.isEmpty()) {
                WrittenGenomicRun whole = applyMeta(GenomicBlocks.concatRuns(legacyParts))
                    .withProvenance(opt.provenanceRecords())
                    .withOptLegacyWholeChannel(true);
                if (opt.embedReference()) {
                    SpectralDatasetGenomicWriter.embedReferencesForRuns(study, List.of(whole));
                }
                SpectralDatasetGenomicWriter.writeGenomicRunSubtree(runsGroup(), name, whole, GenomicWriteContext.none());
                readCount = whole.readCount();
            }
            legacyParts.clear();
            scope.close();
            return;
        }
        try {
            flush();
            drain(0);
            if (rg == null) ensureLayout(null);
            writeCloseTables();
            if (!opt.provenanceRecords().isEmpty()) {
                SpectralDatasetGenomicWriter.writeRunProvenance(rg, opt.provenanceRecords());
            }
        } finally {
            scope.close();
        }
    }

    // ------------------------------------------------------------------

    private StorageGroup runsGroup() {
        StorageGroup g;
        if (study.hasChild("genomic_runs")) {
            g = study.openGroup("genomic_runs");
        } else {
            g = study.createGroup("genomic_runs");
            g.setAttribute("_run_names", "");
        }
        Object namesAttr = g.hasAttribute("_run_names") ? g.getAttribute("_run_names") : null;
        String names = namesAttr == null ? "" : namesAttr.toString();
        List<String> list = new ArrayList<>();
        for (String s : names.split(",")) if (!s.isEmpty()) list.add(s);
        if (!list.contains(name)) {
            list.add(name);
            g.setAttribute("_run_names", String.join(",", list));
        }
        return g;
    }

    /** MD5 of the chromosome sequences concatenated in sorted-name
     *  order: the on-disk {@code @md5} / blob-header digest. */
    private static byte[] referenceSetMd5(Map<String, byte[]> seqs) {
        try {
            java.security.MessageDigest md =
                java.security.MessageDigest.getInstance("MD5");
            List<String> names = new ArrayList<>(seqs.keySet());
            Collections.sort(names);
            for (String n : names) {
                md.update(seqs.get(n));
            }
            return md.digest();
        } catch (java.security.NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    private void ensureLayout(GenomicBlocks.BlockBlobs first) {
        if (rg != null) return;
        StorageGroup g = runsGroup();
        if (g.hasChild(name)) {
            throw new IllegalArgumentException("genomic run '" + name + "' already exists");
        }
        StorageGroup run = g.createGroup(name);
        run.setAttribute("acquisition_mode", (long) opt.acquisitionMode().ordinal());
        run.setAttribute("modality", "genomic_sequencing");
        run.setAttribute("spectrum_class", 5L);
        run.setAttribute("reference_uri", opt.referenceUri());
        run.setAttribute("platform", opt.platform());
        if (opt.readRole() != null && !opt.readRole().isEmpty())
            run.setAttribute("read_role", opt.readRole());
        run.setAttribute("sample_name", opt.sampleName());
        run.setAttribute("read_count", 0L);
        run.setAttribute("base_count", 0L);
        run.setAttribute("layout", LAYOUT);
        run.setAttribute("block_policy", "reads=" + opt.blockReads() + ",bytes=" + opt.blockBytes());
        // Non-default writer policy shapes the coded blobs; a per-AU
        // restore re-encodes with it, so it is persisted on the run
        // (format-spec 9.1.1).
        if (opt.refDiffSliceBytes() != 0) {
            run.setAttribute("ref_diff_slice_bytes", opt.refDiffSliceBytes());
        }
        if (opt.optDisableQualitiesV5()) {
            run.setAttribute("opt_disable_qualities_v5", 1L);
        }
        if (opt.referenceChromSeqs() != null) {
            // The chromosome set and the reference-set digest (the
            // md5 the REF_DIFF blob headers carry and the resolver
            // verifies), so a restore can rebuild the reference from
            // /study/references or REF_PATH (format-spec 9.1.1).
            Map<String, byte[]> seqs = opt.referenceChromSeqs();
            if (referenceMd5 == null) {
                referenceMd5 = referenceSetMd5(seqs);
            }
            StringBuilder hex = new StringBuilder(32);
            for (byte b : referenceMd5) {
                hex.append(String.format("%02x", b));
            }
            List<String> names = new ArrayList<>(seqs.keySet());
            Collections.sort(names);
            StringBuilder json = new StringBuilder("{");
            for (int i = 0; i < names.size(); i++) {
                if (i > 0) json.append(',');
                json.append('"').append(names.get(i)).append("\":\"")
                    .append(hex).append('"');
            }
            json.append('}');
            run.setAttribute("reference_md5s", json.toString());
        }
        StorageGroup blocks = run.createGroup("blocks");
        byte[] firstTags = first == null ? null : first.blobs().get("tags");
        withTags = firstTags != null && firstTags.length > 0;
        indexDs = blocks.createCompoundDataset("index", indexFields(withTags), 0, true, 1024);
        StorageGroup idx = run.createGroup("genomic_index");
        for (String[] a : INDEX_ARRAYS) {
            idxDs.put(a[0], idx.createDataset(a[0], Precision.valueOf(a[1]), 0,
                GenomicIndex.CHUNK_SIZE, Compression.ZLIB, GenomicIndex.COMPRESSION_LEVEL, true));
        }
        run.createGroup("signal_channels");
        rg = run;
    }

    /** The first block with tags arrives after blocks without: give the
     *  index the tags columns, the earlier blocks an empty range (M101). */
    @SuppressWarnings("unchecked")
    private void addTagsColumns() {
        StorageGroup blocks = rg.openGroup("blocks");
        List<Object[]> rows = new ArrayList<>();
        if (blockCount > 0) {
            for (Object[] r : (List<Object[]>) indexDs.readAll()) {
                Object[] wide = java.util.Arrays.copyOf(r, INDEX_FIELDS_WITH_TAGS.size());
                wide[offsetColumn("tags")] = 0L;
                wide[offsetColumn("tags") + 1] = 0L;
                wide[codecColumn("tags")] = 0;
                rows.add(wide);
            }
        }
        indexDs.close();
        blocks.deleteChild("index");
        indexDs = blocks.createCompoundDataset("index", INDEX_FIELDS_WITH_TAGS, 0, true, 1024);
        if (!rows.isEmpty()) indexDs.append(rows);
        withTags = true;
    }

    private StorageDataset createChannel(String ch, GenomicBlocks.BlockBlobs blobs) {
        StorageGroup sc = rg.openGroup("signal_channels");
        StorageGroup parent;
        String dsName;
        if (ch.equals("sequences")) {
            parent = sc.createGroup("sequences");
            dsName = "data";
        } else if (ch.equals("mate_info")) {
            parent = sc.createGroup("mate_info");
            dsName = "inline_v2";
        } else {
            parent = sc;
            dsName = ch;
        }
        int codec = blobs.codecs().get(ch);
        StorageDataset ds = parent.createDataset(dsName, Precision.UINT8, 0, CHANNEL_CHUNK,
            codec == 0 ? Compression.ZLIB : Compression.NONE, 6, true);
        ds.setAttribute("compression", codec);
        for (Map.Entry<String, Object> e : blobs.extraAttrs().get(ch).entrySet()) {
            ds.setAttribute(e.getKey(), e.getValue());
        }
        channelDs.put(ch, ds);
        return ds;
    }

    private void appendIndexArrays(WrittenGenomicRun block) {
        int n = block.readCount();
        short[] ids = new short[n];
        for (int i = 0; i < n; i++) {
            Integer id = chromMap.get(block.chromosomes().get(i));
            if (id == null) throw new IllegalStateException(
                "chromosome '" + block.chromosomes().get(i) + "' missing from the shared id map");
            ids[i] = id.shortValue();
        }
        idxDs.get("lengths").append(block.lengths());
        idxDs.get("positions").append(block.positions());
        idxDs.get("mapping_qualities").append(block.mappingQualities());
        idxDs.get("flags").append(block.flags());
        idxDs.get("chromosome_ids").append(ids);
    }

    private void writeCloseTables() {
        List<CompoundField> nameFields = List.of(
            new CompoundField("name", CompoundField.Kind.VL_STRING));
        List<Object[]> rows = new ArrayList<>();
        for (String n : GenomicIndex.namesInIdOrder(chromMap)) rows.add(new Object[]{ n });
        StorageGroup idx = rg.openGroup("genomic_index");
        try (StorageDataset ds = idx.createCompoundDataset("chromosome_names", nameFields, rows.size())) {
            ds.writeAll(rows);
        }
        StorageGroup sc = rg.openGroup("signal_channels");
        StorageGroup mate = sc.hasChild("mate_info") ? sc.openGroup("mate_info") : sc.createGroup("mate_info");
        if (!mate.hasChild("chrom_names")) {
            try (StorageDataset ds = mate.createCompoundDataset("chrom_names", nameFields, rows.size())) {
                ds.writeAll(rows);
            }
        }
    }

    private WrittenGenomicRun applyMeta(WrittenGenomicRun run) {
        return new WrittenGenomicRun(
            opt.acquisitionMode(), opt.referenceUri(), opt.platform(), opt.sampleName(),
            run.positions(), run.mappingQualities(), run.flags(), run.sequences(), run.qualities(),
            run.offsets(), run.lengths(), run.cigars(), run.readNames(), run.mateChromosomes(),
            run.matePositions(), run.templateLengths(), run.chromosomes(),
            opt.signalCompression(), opt.signalCodecOverrides(), List.of(),
            opt.embedReference(), opt.referenceChromSeqs(), null, null,
            opt.optDisableQualitiesV5(), false,
            opt.readRole(), opt.refDiffSliceBytes(), run.tags());
    }

    private WrittenGenomicRun singleReadRun(AlignedRead r) {
        byte[] seq = r.sequence() == null ? new byte[0] : r.sequence().getBytes(StandardCharsets.US_ASCII);
        byte[] qual = r.qualities() == null ? new byte[0] : r.qualities();
        String mate = (r.mateChromosome() == null || r.mateChromosome().isEmpty()) ? "*" : r.mateChromosome();
        return new WrittenGenomicRun(
            opt.acquisitionMode(), opt.referenceUri(), opt.platform(), opt.sampleName(),
            new long[]{ r.position() }, new byte[]{ (byte) r.mappingQuality() }, new int[]{ r.flags() },
            seq, qual, new long[]{ 0L }, new int[]{ seq.length },
            List.of(r.cigar() == null ? "*" : r.cigar()), List.of(r.readName() == null ? "*" : r.readName()),
            List.of(mate), new long[]{ r.matePosition() }, new int[]{ r.templateLength() },
            List.of(r.chromosome() == null ? "*" : r.chromosome()),
            opt.signalCompression(), opt.signalCodecOverrides(), List.of(),
            opt.embedReference(), opt.referenceChromSeqs(), null, null,
            opt.optDisableQualitiesV5(), false,
            opt.readRole(), opt.refDiffSliceBytes(),
            List.of(r.tags() == null ? "" : r.tags()));
    }
}
