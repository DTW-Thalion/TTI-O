/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.providers.StorageDataset;
import global.thalion.ttio.providers.StorageGroup;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

/**
 * Lazy view over one {@code /study/genomic_runs/<name>/} group.
 *
 * <p>Materialises {@link AlignedRead} objects on demand from the
 * signal channels. The {@link GenomicIndex} is loaded eagerly at
 * open time for cheap filtering and offset lookups; the heavy signal
 * channels (sequences, qualities, plus 3 compounds) stay lazy.</p>
 *
 * <p>Genomic analogue of
 * {@link global.thalion.ttio.AcquisitionRun}.</p>
 *
 * <p><b>API status:</b> Stable ().</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOGenomicRun}, Python {@code ttio.genomic_run.GenomicRun}.</p>
 */
public class GenomicRun
        implements global.thalion.ttio.protocols.Indexable<AlignedRead>,
                   global.thalion.ttio.protocols.Streamable<AlignedRead>,
                   global.thalion.ttio.protocols.Run,
                   AutoCloseable {

    private final String name;
    private final AcquisitionMode acquisitionMode;
    private final String modality;
    private final String referenceUri;
    private final String platform;
    private final String sampleName;
    private GenomicIndex index;
    private final StorageGroup runGroup;
    // blocks_v1 (format-spec 10.12): the block table, the last
    // materialised block view and the run-level name tables. Null /
    // "whole" for the v1.8 whole-channel layout.
    private final String layout;
    private final BlockTable blockTable;
    private int cachedBlock = -1;
    private GenomicRun cachedView;
    private BlockView.Handle cachedHandle;
    private List<String> chromNamesTable;
    private List<String> mateChromNamesTable;
    private global.thalion.ttio.codecs.ReferenceResolver injectedResolver;
    private global.thalion.ttio.codecs.ReferenceResolver viewResolver;
    private boolean viewResolverBuilt;
    // M97: the @read_role run attribute; null when absent.
    private String readRole;
    // Phase 1 (post-M91): per-run provenance, cached at open time so
    // provenanceChain() is a pure accessor. Eager because the on-disk
    // form is a small JSON attribute on the run group; lazy decode
    // would buy nothing and would complicate the Run protocol surface.
    private final List<ProvenanceRecord> provenanceRecords;

    // ── Per-record HDF5-probe cache invariants ────────────────────────
    //
    // Every cache below memoises a per-record-path operation that
    // would otherwise hit HDF5 (open group / open dataset / read
    // attribute) on every call. objectAtIndex(i) drives one call to
    // each of: byteChannelSlice("sequences"), byteChannelSlice(
    // "qualities"), cigarAt(i), readNameAt(i), mateChromAt(i),
    // matePosAt(i), mateTlenAt(i). Any uncached probe in that chain
    // multiplies per-record cost by HDF5 round-trip overhead, which
    // dominates the genomic transport encode hot path at scale.
    //
    //   * sequencesIsV2Cached    — covers isSequencesRefDiffV2
    //   * mateInfoInlineV2Cached — covers isMateInfoInlineV2
    //   * decodedByteChannels    — covers byteChannelSlice (compressed +
    //                              uncompressed via byteChannelFull)
    //   * decodedReadNames       — covers readNameAt
    //   * decodedCigars          — covers cigarAt
    //   * decodedMateV2          — covers mate{Chrom,Pos,Tlen}At
    //   * compoundCache          — covers compoundRows
    //
    // signalChannelCompressionCode opens a dataset per call but is
    // only invoked twice per run (sequences + qualities at the top
    // of TransportWriter.emitGenomicRunAccessUnits), not per record.
    private StorageGroup signalChannels;                       // lazy
    private StorageDataset sequencesDs;                        // lazy
    private StorageDataset qualitiesDs;                        // lazy
    private final Map<String, List<Object[]>> compoundCache = new HashMap<>();
    // lazy whole-channel decode cache for byte channels whose
    // @compression attribute names a TTI-O codec (rANS / BASE_PACK).
    // Codec output is byte-stream non-sliceable, so the whole channel
    // is decoded once on first access and the decoded buffer is sliced
    // from memory thereafter (). Cache lifetime is
    // this GenomicRun instance — re-opening the file incurs the decode
    // cost again (Gotcha §101). Mutable HashMap, not Map.of(), so the
    // dispatch helper can populate it.
    private final Map<String, byte[]> decodedByteChannels = new HashMap<>();
    // Task 4 (codec registry): cached run-derived CodecContext. Built
    // once on first decode and reused for every codec.decode() call on
    // this GenomicRun instance. Mirrors the per-run lifetime of the
    // other decode caches above.
    private global.thalion.ttio.codecs.registry.CodecContext codecCtxCache;
    // lazy decode cache for the read_names channel when it
    // carries a NAME_TOKENIZED codec override. Held as a List<String>
    // (not byte[]) — different value type and
    // semantics from decodedByteChannels (which holds raw byte buffers
    // sliced by per-read offset/length). The whole list is materialised
    // on first access regardless of the access pattern.
    private List<String> decodedReadNames = null;
    // lazy decode cache for the cigars channel when it
    // carries an RANS_ORDER0 / RANS_ORDER1 / NAME_TOKENIZED codec
    // override. Held as a List<String> mirroring decodedReadNames per
    // — separate cache from decodedReadNames
    // (Option A from §2.3, lower-risk than a generalised dict).
    private List<String> decodedCigars = null;
    /** M101: decoded SAM_TAGS channel, one string per read; empty when
     *  the run has no tags channel. Null until first access. */
    private List<String> decodedTags = null;
    // lazy decode cache for integer channels. Per Binding
    // Decision §116 this is a separate cache from decodedByteChannels
    // (byte[]) and decodedReadNames (List<String>) because the value
    // v1.6 (L4): decodedIntChannels removed. The cache supported the
    // intChannelArray helper which read positions/flags/mapping_qualities
    // from signal_channels/ via codec dispatch — but those datasets no
    // longer exist in v1.6 files (they live exclusively in
    // genomic_index/). See docs/format-spec.md §10.7.
    // Task 13 (mate_info v2): lazy decoded triple from inline_v2 blob.
    // Null until first access to a mate field on a v2-layout file.
    private global.thalion.ttio.codecs.MateInfoV2.Triple decodedMateV2 = null;
    /** Memoised result of {@link #isMateInfoInlineV2}. The probe opens
     *  the {@code mate_info} HDF5 group; without caching, calling it
     *  three times per record from {@link #objectAtIndex} blew up the
     *  TransportWriter genomic encode path (300K HDF5 group opens
     *  per 100K reads, ~2.2s of pure framework overhead). */
    private Boolean mateInfoInlineV2Cached = null;
    // Resolved chrom name table for the v2 path: index → chrom name.
    // Task 13 (ref_diff v2): lazy decoded flat byte stream from the
    // signal_channels/sequences/refdiff_v2 blob. Null until first access
    // on a v1.8-layout file. Separate from decodedByteChannels["sequences"]
    // because the source is a group child, not the sequences dataset.
    private byte[] decodedRefDiffV2 = null;
    // Tri-state cache for isSequencesRefDiffV2(): null = not yet probed,
    // TRUE/FALSE = result. Avoids repeated group-open on every byteChannelSlice call.
    private Boolean sequencesIsV2Cached = null;
    private List<String> mateV2ChromNames = null;

    private int cursor = 0;  // Streamable

    private GenomicRun(String name, AcquisitionMode acquisitionMode,
                       String modality, String referenceUri,
                       String platform, String sampleName,
                       GenomicIndex index, StorageGroup runGroup,
                       List<ProvenanceRecord> provenanceRecords,
                       String layout, BlockTable blockTable) {
        this.name = name;
        this.acquisitionMode = acquisitionMode;
        this.modality = modality;
        this.referenceUri = referenceUri;
        this.platform = platform;
        this.sampleName = sampleName;
        this.index = index;
        this.runGroup = runGroup;
        this.provenanceRecords = provenanceRecords != null
            ? List.copyOf(provenanceRecords) : List.of();
        this.layout = layout;
        this.blockTable = blockTable;
    }

    public String name()                       { return name; }
    public AcquisitionMode acquisitionMode()   { return acquisitionMode; }
    public String modality()                   { return modality; }
    public String referenceUri()               { return referenceUri; }
    public String platform()                   { return platform; }
    public String sampleName()                 { return sampleName; }

    /** The {@code @read_role} run attribute (M97): {@code hifi},
     *  {@code ont_ul}, {@code hic_r1}, {@code hic_r2},
     *  {@code parental_maternal}, {@code parental_paternal}, or
     *  {@code illumina_wgs} (other strings are stored unchecked);
     *  null when the attribute is absent. Python: {@code read_role};
     *  ObjC: {@code readRole}. */
    public String getReadRole()                { return readRole; }

    /** The per-read index. Under {@code blocks_v1} it is loaded from
     *  {@code genomic_index/} on first call. */
    public GenomicIndex index() {
        if (index == null) {
            try (StorageGroup ig = runGroup.openGroup("genomic_index")) {
                index = GenomicIndex.readFrom(ig);
            }
        }
        return index;
    }

    public int readCount() {
        return blockTable != null ? (int) blockTable.readCount() : index().count();
    }

    /** {@code "blocks_v1"} or {@code "whole"} (the v1.8 whole-channel
     *  layout). */
    public String layout() { return layout; }

    /** Number of blocks; 1 for a whole-channel run. */
    public int blockCount() { return blockTable != null ? blockTable.count() : 1; }

    /** Run-level chromosome name table ({@code genomic_index/chromosome_names}),
     *  read without loading the per-read arrays. */
    public List<String> chromosomeNames() {
        if (chromNamesTable == null) {
            try (StorageGroup ig = runGroup.openGroup("genomic_index")) {
                chromNamesTable = BlockView.readNames(ig, "chromosome_names");
            }
        }
        return chromNamesTable;
    }

    /** Open an existing genomic_runs/&lt;name&gt;/ group. The caller
     *  resolves the run group and passes it as {@code runGroup}. */
    public static GenomicRun readFrom(StorageGroup runGroup, String name) {
        return readFrom(runGroup, name, null);
    }

    /** {@code resolver}, when given, is used for REF_DIFF decodes instead
     *  of one built from the run group's HDF5 file (block views live in
     *  a memory provider and share their parent's). */
    static GenomicRun readFrom(StorageGroup runGroup, String name,
                               global.thalion.ttio.codecs.ReferenceResolver resolver) {
        String layout = stringAttr(runGroup, "layout", "whole");
        GenomicIndex idx = null;
        BlockTable table = null;
        if ("blocks_v1".equals(layout)) {
            table = BlockTable.read(runGroup);
        } else if ("whole".equals(layout)) {
            try (StorageGroup ig = runGroup.openGroup("genomic_index")) {
                idx = GenomicIndex.readFrom(ig);
            }
        } else {
            throw new IllegalStateException(
                "genomic run '" + name + "': unsupported layout '" + layout
                + "' (this reader knows the whole-channel layout and blocks_v1)");
        }
        Object modeObj = runGroup.getAttribute("acquisition_mode");
        AcquisitionMode mode = AcquisitionMode.values()[
            modeObj == null ? 0 : ((Number) modeObj).intValue()];
        String modality   = stringAttr(runGroup, "modality",       "genomic_sequencing");
        String refUri     = stringAttr(runGroup, "reference_uri",  "");
        String platform   = stringAttr(runGroup, "platform",       "");
        String sampleName = stringAttr(runGroup, "sample_name",    "");
        List<ProvenanceRecord> prov = readPerRunProvenance(runGroup);
        GenomicRun run = new GenomicRun(name, mode, modality, refUri, platform,
                                        sampleName, idx, runGroup, prov, layout, table);
        run.injectedResolver = resolver;
        // @read_role is absent on pre-M97 files; null, not "".
        String role = stringAttr(runGroup, "read_role", null);
        run.readRole = (role == null || role.isEmpty()) ? null : role;
        return run;
    }

    // ── blocks_v1 dispatch ─────────────────────────────────────────

    /** The {@link GenomicRun} over block {@code b}, materialised on
     *  demand; the last one is cached. */
    /** The view handle for block {@code b}, materialised on the caller's
     *  thread (storage reads). */
    private BlockView.Handle materialiseHandle(int b) {
        if (mateChromNamesTable == null) {
            try (StorageGroup sc = runGroup.openGroup("signal_channels")) {
                mateChromNamesTable = sc.hasChild("mate_info")
                    ? BlockView.readNames(sc.openGroup("mate_info"), "chrom_names")
                    : List.of();
            }
        }
        return BlockView.materialise(runGroup, blockTable, b,
                chromosomeNames(), mateChromNamesTable);
    }

    private GenomicRun blockView(int b) {
        if (cachedView != null && cachedBlock == b) return cachedView;
        BlockView.Handle h = materialiseHandle(b);
        GenomicRun sub = GenomicRun.readFrom(h.group(), name, resolverForViews());
        dropCachedView();
        cachedBlock = b;
        cachedView = sub;
        cachedHandle = h;
        return sub;
    }

    private void dropCachedView() {
        if (cachedView != null) cachedView.close();
        if (cachedHandle != null) cachedHandle.discard();
        cachedView = null;
        cachedHandle = null;
        cachedBlock = -1;
    }

    private global.thalion.ttio.codecs.ReferenceResolver resolverForViews() {
        if (injectedResolver != null) return injectedResolver;
        if (!viewResolverBuilt) {
            viewResolverBuilt = true;
            try {
                global.thalion.ttio.hdf5.Hdf5Group h5g = global.thalion.ttio.providers
                    .Hdf5Provider.tryUnwrapHdf5Group(runGroup);
                viewResolver = h5g == null ? null
                    : new global.thalion.ttio.codecs.ReferenceResolver(h5g.owningFile());
            } catch (RuntimeException e) {
                viewResolver = null;
            }
        }
        return viewResolver;
    }

    /** Reads {@code [start, stop)} in order, holding at most one decoded
     *  block at a time under {@code blocks_v1}. */
    public java.util.Iterator<AlignedRead> iterReads(int start, int stop) {
        int n = readCount();
        int lo = Math.max(start, 0), hi = Math.min(stop, n);
        return new java.util.Iterator<>() {
            int i = lo;
            @Override public boolean hasNext() { return i < hi; }
            @Override public AlignedRead next() {
                if (i >= hi) throw new NoSuchElementException();
                return objectAtIndex(i++);
            }
        };
    }

    /** Every read in order; see {@link #iterReads(int, int)}. */
    public java.util.Iterator<AlignedRead> iterReads() { return iterReads(0, readCount()); }

    private record InFlightView(GenomicRun view, BlockView.Handle handle) {}

    /** Reads {@code [start, stop)} in order; under {@code blocks_v1} the
     *  next blocks decode ahead on a pool, at most
     *  {@link #READ_AHEAD_BLOCKS} of them decoded at once. {@code threads}
     *  {@code <= 1} is the serial iterator. */
    /** Decoded blocks a threaded sequential walk holds at most (the
     *  decode-ahead window; memory is about this many block working sets). */
    static final int READ_AHEAD_BLOCKS = 4;

    /** The window, with {@code TTIO_READ_AHEAD_BLOCKS} allowed to
     *  override it. Each block in flight stays resident until it is
     *  consumed, so this is a memory setting as much as a latency one,
     *  and the two do not pull the same way. The override exists so the
     *  trade can be measured rather than argued. */
    static int readAheadBlocks() {
        String env = System.getenv("TTIO_READ_AHEAD_BLOCKS");
        if (env != null && !env.isEmpty()) {
            try {
                int v = Integer.parseInt(env.trim());
                if (v > 0) return v;
            } catch (NumberFormatException ignored) {
                // fall through to the default
            }
        }
        return READ_AHEAD_BLOCKS;
    }

    public java.util.Iterator<AlignedRead> iterReads(int start, int stop, int threads) {
        int n = readCount();
        int lo = Math.max(start, 0), hi = Math.min(stop, n);
        int nthreads = Math.max(1, threads);
        if (blockTable == null || nthreads <= 1 || lo >= hi) return iterReads(lo, hi);
        final int bFirst = blockTable.blockFor(lo), bLast = blockTable.blockFor(hi - 1);
        // The serial consumer only needs enough decode-ahead to never
        // stall; more in flight is pure memory.
        final int window = Math.min(nthreads, readAheadBlocks());
        // The scope takes the window, not nthreads: it sizes V6's
        // segment threads from the number of workers, and the window is
        // how many blocks are ever in flight.
        final global.thalion.ttio.Threads.PoolScope scope = global.thalion.ttio.Threads.pool(window);
        final java.util.Map<Integer, java.util.concurrent.Future<InFlightView>> pending = new java.util.HashMap<>();
        final java.util.function.IntConsumer submit = b -> {
            if (b <= bLast && !pending.containsKey(b)) {
                BlockView.Handle h = materialiseHandle(b);       // storage reads, this thread
                pending.put(b, scope.executor().submit(() -> {
                    GenomicRun v = GenomicRun.readFrom(h.group(), name, resolverForViews());
                    if (v.readCount() > 0) v.readAt(0);           // warm every channel cache
                    return new InFlightView(v, h);
                }));
            }
        };
        for (int b = bFirst; b <= Math.min(bLast, bFirst + window - 1); b++) submit.accept(b);
        return new java.util.Iterator<>() {
            int i = lo, b = bFirst;
            InFlightView cur;
            int r0, bEnd;
            @Override public boolean hasNext() {
                if (i < hi) return true;
                release();
                scope.close();
                for (var f : pending.values()) {
                    try { f.get().handle().discard(); } catch (Exception ignored) { }
                }
                pending.clear();
                return false;
            }
            private void release() {
                if (cur != null) { cur.view().close(); cur.handle().discard(); cur = null; }
            }
            @Override public AlignedRead next() {
                if (i >= hi) throw new NoSuchElementException();
                if (cur == null || i >= bEnd) {
                    release();
                    try {
                        cur = pending.remove(b).get();
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        scope.close();
                        throw new IllegalStateException(e);
                    } catch (java.util.concurrent.ExecutionException e) {
                        scope.close();
                        Throwable c = e.getCause();
                        throw c instanceof RuntimeException r ? r : new IllegalStateException(c);
                    }
                    submit.accept(b + window);
                    r0 = (int) blockTable.readStart[b];
                    bEnd = Math.min(r0 + blockTable.nReads[b], hi);
                    b++;
                }
                return cur.view().objectAtIndex(i++ - r0);
            }
        };
    }

    /** One decoded block, handed over whole. */
    @FunctionalInterface
    public interface BlockVisitor {
        /** {@code view} is a run over one block's reads and is valid
         *  only for the duration of the call. Read record {@code k} of
         *  the delivered range as {@code view.readAt(viewStart + k)};
         *  its index in the whole run is {@code firstRead + k}. The two
         *  differ whenever a range starts part-way into a block, so a
         *  caller that indexes the view by {@code firstRead} reads the
         *  wrong records. */
        void visit(GenomicRun view, int viewStart, int firstRead, int nReads);
    }

    /** How many blocks may be in flight at once, given how many threads
     *  want one each. A block in flight is resident for as long as the
     *  caller is inside it, so this is a memory setting before it is a
     *  concurrency one, and the thread count is only an upper bound.
     *
     *  <p>A block in flight is charged 8 times its sequence bytes, the
     *  figure {@link global.thalion.ttio.Threads#resolveMemoryBudget}
     *  uses on the writing side, so {@code TTIO_MEMORY_BUDGET} and the
     *  half-of-physical cap mean the same thing to a reader as to a
     *  writer. The widest block in the range sets the size, not the
     *  mean.
     *
     *  <p>Objective-C: {@code -_blockWindowForThreads:first:last:}. */
    private int blockWindow(int nthreads, int bFirst, int bLast) {
        long widest = 0;
        for (int b = bFirst; b <= bLast; b++) widest = Math.max(widest, blockTable.nBases[b]);
        if (widest == 0) return 1;
        long budget = global.thalion.ttio.Threads.resolveMemoryBudget(null, nthreads, widest);
        long admits = Math.max(1, budget / (widest * 8L));
        int blocks = bLast - bFirst + 1;
        return (int) Math.max(1, Math.min(Math.min(nthreads, admits), blocks));
    }

    /** One call per decoded block, on the pool.
     *
     *  <p>{@code fn} is called from several threads at once and in no
     *  particular order, and must be safe to call that way. That
     *  relaxed ordering is the whole difference from
     *  {@link #iterReads(int, int, int)}, whose in-order delivery on the
     *  caller's thread is what bounds it: the decoder outruns the
     *  consumer, so dividing the consumer is what buys the throughput.
     *
     *  <p>Python: {@code GenomicRun.for_each_block}; Objective-C:
     *  {@code -iterBlocksFrom:to:threads:error:usingBlock:}. */
    public void iterBlocks(int start, int stop, int threads, BlockVisitor fn) {
        int n = readCount();
        final int lo = Math.max(start, 0), hi = Math.min(stop, n);
        if (lo >= hi) return;
        if (blockTable == null) { fn.visit(this, lo, lo, hi - lo); return; }
        int nthreads = Math.max(1, threads == 0
            ? global.thalion.ttio.Threads.resolve(null) : threads);
        final int bFirst = blockTable.blockFor(lo), bLast = blockTable.blockFor(hi - 1);
        final int window = blockWindow(nthreads, bFirst, bLast);
        final int blocks = bLast - bFirst + 1;
        final int aheadN = Math.min(blocks, Math.max(2, readAheadBlocks()));

        /* Below the decode-ahead depth, one consumer fed by that many
         * decoders beats the same number of threads each decoding and
         * then consuming its own block: decoding is the slower half, so
         * until there are threads enough to cover it, pipelining wins
         * over dividing. */
        if (window < aheadN) {
            iterBlocksPipelined(lo, hi, bFirst, bLast,
                                Math.min(aheadN, Math.max(window, 1) * 4), fn);
            return;
        }
        final global.thalion.ttio.Threads.PoolScope scope =
            global.thalion.ttio.Threads.pool(window);
        final java.util.Map<Integer, java.util.concurrent.Future<?>> pending =
            new java.util.HashMap<>();
        try {
            final java.util.function.IntConsumer submit = b -> {
                if (b > bLast || pending.containsKey(b)) return;
                BlockView.Handle h = materialiseHandle(b);       // storage reads, this thread
                final int r0 = (int) blockTable.readStart[b];
                final int nr = blockTable.nReads[b];
                pending.put(b, scope.executor().submit(() -> {
                    try {
                        GenomicRun v = GenomicRun.readFrom(h.group(), name, resolverForViews());
                        int from = Math.max(r0, lo), to = Math.min(r0 + nr, hi);
                        if (to > from) fn.visit(v, from - r0, from, to - from);
                        v.close();
                    } finally {
                        h.discard();
                    }
                    return null;
                }));
            };
            for (int b = bFirst; b <= Math.min(bLast, bFirst + window - 1); b++) submit.accept(b);
            for (int b = bFirst; b <= bLast; b++) {
                var f = pending.remove(b);
                if (f == null) continue;
                try {
                    f.get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new IllegalStateException(e);
                } catch (java.util.concurrent.ExecutionException e) {
                    Throwable c = e.getCause();
                    throw c instanceof RuntimeException r ? r : new IllegalStateException(c);
                }
                submit.accept(b + window);
            }
        } finally {
            for (var f : pending.values()) { try { f.get(); } catch (Exception ignored) { } }
            scope.close();
        }
    }

    /** One consumer on this thread, {@code aheadN} decoders behind it.
     *  {@code readFrom} is lazy, so the warming read has to happen on
     *  the pool: without it the channel decode lands on the consumer's
     *  thread at its first read and the overlap is lost entirely. */
    private void iterBlocksPipelined(int lo, int hi, int bFirst, int bLast,
                                     int aheadN, BlockVisitor fn) {
        final global.thalion.ttio.Threads.PoolScope scope =
            global.thalion.ttio.Threads.pool(Math.max(2, aheadN));
        final java.util.Map<Integer, java.util.concurrent.Future<InFlightView>> pending =
            new java.util.HashMap<>();
        try {
            final java.util.function.IntConsumer submit = b -> {
                if (b > bLast || pending.containsKey(b)) return;
                BlockView.Handle h = materialiseHandle(b);
                pending.put(b, scope.executor().submit(() -> {
                    GenomicRun v = GenomicRun.readFrom(h.group(), name, resolverForViews());
                    if (v.readCount() > 0) v.readAt(0);          // decode, off the consumer
                    return new InFlightView(v, h);
                }));
            };
            for (int b = bFirst; b <= Math.min(bLast, bFirst + aheadN - 1); b++) submit.accept(b);
            for (int b = bFirst; b <= bLast; b++) {
                InFlightView cur;
                try {
                    cur = pending.remove(b).get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new IllegalStateException(e);
                } catch (java.util.concurrent.ExecutionException e) {
                    Throwable c = e.getCause();
                    throw c instanceof RuntimeException r ? r : new IllegalStateException(c);
                }
                submit.accept(b + aheadN);
                int r0 = (int) blockTable.readStart[b];
                int from = Math.max(r0, lo), to = Math.min(r0 + blockTable.nReads[b], hi);
                try {
                    if (to > from) fn.visit(cur.view(), from - r0, from, to - from);
                } finally {
                    cur.view().close();
                    cur.handle().discard();
                }
            }
        } finally {
            for (var f : pending.values()) {
                try { var v = f.get(); v.view().close(); v.handle().discard(); }
                catch (Exception ignored) { }
            }
            scope.close();
        }
    }

    /** Phase 2 (post-M91): read per-run provenance. Prefers the
     *  canonical compound dataset {@code provenance/steps} (matches
     *  Python's writer and Java's HDF5 fast path), falling back to
     *  the {@code provenance_json} attribute for non-HDF5 providers
     *  (memory/sqlite/zarr) and legacy Java-written files. */
    private static List<ProvenanceRecord> readPerRunProvenance(StorageGroup runGroup) {
        if (runGroup.hasChild("provenance")) {
            try (StorageGroup prov = runGroup.openGroup("provenance")) {
                global.thalion.ttio.hdf5.Hdf5Group h5 =
                    global.thalion.ttio.providers.Hdf5Provider
                        .tryUnwrapHdf5Group(prov);
                if (h5 != null && h5.hasChild("steps")) {
                    List<Object[]> rows =
                        global.thalion.ttio.hdf5.Hdf5CompoundIO
                            .readCompoundFull(h5, "steps",
                                global.thalion.ttio.hdf5.Hdf5CompoundIO
                                    .provenanceSchema());
                    List<ProvenanceRecord> out = new ArrayList<>(rows.size());
                    for (Object[] r : rows) {
                        out.add(new ProvenanceRecord(
                            ((Number) r[0]).longValue(),
                            (String) r[1],
                            global.thalion.ttio.MiniJson.parseStringMap(
                                (String) r[2]),
                            global.thalion.ttio.MiniJson.parseArrayOfStrings(
                                (String) r[3]),
                            global.thalion.ttio.MiniJson.parseArrayOfStrings(
                                (String) r[4])));
                    }
                    return out;
                }
            }
        }
        if (!runGroup.hasAttribute("provenance_json")) {
            return List.of();
        }
        Object v = runGroup.getAttribute("provenance_json");
        if (v == null) return List.of();
        String json = v instanceof String s ? s
                    : v instanceof byte[] b ? new String(b, StandardCharsets.UTF_8)
                    : v.toString();
        return global.thalion.ttio.ProvenanceJsonParse.parseArray(json);
    }

    /** probe the {@code @compression} attribute on a
     *  signal_channels child dataset. Returns the codec id (an
     *  {@link global.thalion.ttio.Enums.Compression} ordinal), or 0
     *  ({@code NONE}) when the attribute is absent or the channel
     *  doesn't exist. Used by
     *  {@link global.thalion.ttio.transport.TransportWriter} to mirror
     *  the file's per-channel codec choice on the wire. */
    public int signalChannelCompressionCode(String channelName) {
        if (blockTable != null) {
            if (blockTable.count() == 0 || blockTable.codec == null
                    || !blockTable.codec.containsKey(channelName)) return 0;
            int c = blockTable.codec.get(channelName)[0];
            return (c == global.thalion.ttio.Enums.Compression.RANS_ORDER0.ordinal()
                 || c == global.thalion.ttio.Enums.Compression.RANS_ORDER1.ordinal()
                 || c == global.thalion.ttio.Enums.Compression.BASE_PACK.ordinal()) ? c : 0;
        }
        ensureSignalChannels();
        if (!signalChannels.hasChild(channelName)) return 0;
        try (StorageDataset ds = signalChannels.openDataset(channelName)) {
            Object v = ds.getAttribute("compression");
            if (v instanceof Number n) return n.intValue();
            return 0;
        } catch (Exception e) {
            return 0;
        }
    }

    /** Materialise read at index {@code i}. {@link Indexable} requires
     *  this signature. The shorthand {@code readAt} is provided as a
     *  domain-natural alias. */
    @Override
    public AlignedRead objectAtIndex(int i) {
        if (i < 0 || i >= readCount()) {
            throw new IndexOutOfBoundsException(
                "read index " + i + " out of range [0, " + readCount() + ")");
        }
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).objectAtIndex(i - (int) blockTable.readStart[b]);
        }
        long offset = index().offsetAt(i);
        int  length = index().lengthAt(i);

        ensureSignalChannels();
        // routed through byteChannelSlice so that channels written
        // with a TTIO codec override (@compression > 0) are decoded
        // transparently before slicing.
        byte[] seqBytes = byteChannelSlice("sequences", offset, length);
        String sequence = new String(seqBytes, StandardCharsets.US_ASCII);
        byte[] qualities = byteChannelSlice("qualities", offset, length);

        // read_names dispatched separately via
        // readNameAt() (the dataset shape varies — compound vs flat
        // uint8 codec). M86 Phase C: cigars likewise via cigarAt().
        // mate fields dispatched via three per-field
        // accessors (mateChromAt / matePosAt / mateTlenAt) since the
        // mate_info link can be either an M82 compound dataset OR a
        // Phase F subgroup containing three child datasets (Binding
        // Decision §128, link-type dispatch).
        String cigar     = cigarAt(i);
        String readName  = readNameAt(i);
        String mateChrom = mateChromAt(i);
        long   matePos   = matePosAt(i);
        int    tlen      = mateTlenAt(i);

        List<String> tags = allTags();

        return new AlignedRead(
            readName,
            index().chromosomeAt(i),
            index().positionAt(i),
            index().mappingQualityAt(i),
            cigar,
            sequence,
            qualities,
            index().flagsAt(i),
            mateChrom,
            matePos,
            tlen,
            tags.isEmpty() ? "" : tags.get(i));
    }

    /** True when the run carries the SAM tags channel (M101). */
    public boolean hasTagsChannel() {
        if (hasTagsCached == null) {
            boolean has = false;
            if (runGroup.hasChild("signal_channels")) {
                try (StorageGroup sig = runGroup.openGroup("signal_channels")) {
                    has = sig.hasChild("tags");
                }
            }
            hasTagsCached = has;
        }
        return hasTagsCached;
    }
    private Boolean hasTagsCached;

    /** Read {@code i}'s SAM tag text (M101); {@code ""} when the run
     *  has no tags channel. */
    public String tagsAt(int i) {
        if (i < 0 || i >= readCount()) {
            throw new IndexOutOfBoundsException(
                "read index " + i + " out of range [0, " + readCount() + ")");
        }
        if (blockTable != null) {
            if (!hasTagsChannel()) return "";
            int b = blockTable.blockFor(i);
            return blockView(b).tagsAt(i - (int) blockTable.readStart[b]);
        }
        List<String> tags = allTags();
        return tags.isEmpty() ? "" : tags.get(i);
    }

    /** Every read's SAM tag text (M101), decoded once; empty when the
     *  run has no {@code tags} channel. */
    private List<String> allTags() {
        if (decodedTags != null) return decodedTags;
        ensureSignalChannels();
        if (!signalChannels.hasChild("tags")) {
            decodedTags = List.of();
            return decodedTags;
        }
        byte[] blob;
        long codecId;
        try (StorageDataset ds = signalChannels.openDataset("tags")) {
            Object v = ds.getAttribute("compression");
            codecId = v instanceof Number num ? num.longValue() : 0L;
            blob = (byte[]) ds.readSlice(0L, ds.shape()[0]);
        }
        if (codecId != global.thalion.ttio.Enums.Compression.SAM_TAGS.ordinal()) {
            throw new IllegalStateException(
                "signal_channel 'tags': @compression=" + codecId
                + " is not SAM_TAGS (18), the only codec of the tags channel");
        }
        decodedTags = ((global.thalion.ttio.codecs.registry.DecodedChannel.StrList)
            global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                .get(global.thalion.ttio.Enums.Compression.SAM_TAGS)
                .decode(new global.thalion.ttio.codecs.registry.ChannelPayload.BytesPayload(blob),
                    codecContext())).names();
        return decodedTags;
    }

    /** References for SAM_TAGS MD/NM derivation, indexed like
     *  {@code ownChromIds}. The writer derives only on REF_DIFF_V2 runs
     *  (binding decision §96), whose one chromosome resolves the same
     *  way the sequences decode does; any other run gets none. */
    private List<byte[]> tagReferences(List<String> chromTable,
            global.thalion.ttio.codecs.ReferenceResolver resolver) {
        if (!isSequencesRefDiffV2() || chromTable == null) return List.of();
        byte[] head;
        try (StorageGroup seqGrp = signalChannels.openGroup("sequences");
             StorageDataset ds = seqGrp.openDataset("refdiff_v2")) {
            // The outer header is 38 bytes plus a uint16-length URI.
            long len = Math.min(ds.shape()[0], 38L + 0xFFFF);
            head = (byte[]) ds.readSlice(0L, len);
        }
        RefDiffHeader h = RefDiffHeader.of(head);
        java.util.TreeSet<String> names = new java.util.TreeSet<>();
        for (int i = 0; i < index().count(); i++) names.add(index().chromosomeAt(i));
        names.remove("*");
        names.remove("");
        // Same name -> id map as the ownChromIds derivation.
        java.util.Map<String, Integer> nameToId = new java.util.LinkedHashMap<>();
        for (int j = 0; j < chromTable.size(); j++) nameToId.put(chromTable.get(j), j);
        List<byte[]> refs = new ArrayList<>(java.util.Collections.nCopies(chromTable.size(), (byte[]) null));
        if (names.size() == 1 && nameToId.containsKey(names.first())) {
            if (resolver == null) {
                throw new IllegalStateException(
                    "SAM_TAGS MD/NM decode needs a reference resolver");
            }
            refs.set(nameToId.get(names.first()),
                resolver.resolve(h.uri(), h.md5(), names.first()));
        }
        return refs;
    }

    private record RefDiffHeader(String uri, byte[] md5) {
        static RefDiffHeader of(byte[] head) {
            var bh = global.thalion.ttio.codecs.RefDiffV2.parseBlobHeader(head);
            return new RefDiffHeader(bh.referenceUri(), bh.referenceMd5());
        }
    }

    /** Domain-natural alias for {@link #objectAtIndex(int)}. */
    public AlignedRead readAt(int i) { return objectAtIndex(i); }

    /** Reads on {@code chromosome} whose mapping position is in
     *  {@code [start, end)}. */
    public List<AlignedRead> readsInRegion(String chromosome,
                                            long start, long end) {
        List<Integer> indices = index().indicesForRegion(chromosome, start, end);
        List<AlignedRead> out = new ArrayList<>(indices.size());
        for (int i : indices) out.add(objectAtIndex(i));
        return out;
    }

    // ── Indexable<AlignedRead> ─────────────────────────────────────

    @Override public int count() { return readCount(); }

    // ── Run conformance (Phase 1, post-M91) ────────────────────────

    /** Phase 1: modality-agnostic accessor required by
     *  {@link global.thalion.ttio.protocols.Run}. Delegates to
     *  {@link #objectAtIndex(int)}; the typed return is widened to
     *  {@code Object} so callers iterating uniformly over
     *  AcquisitionRun + GenomicRun see a single signature. */
    @Override
    public Object get(int index) { return objectAtIndex(index); }

    /** Phase 1 (post-M91): per-run provenance chain in insertion
     *  order. Closes the M91 read-side gap — until Phase 1 the lazy
     *  {@code GenomicRun} view didn't expose provenance, so cross-
     *  modality queries had to fall back to the {@code @sample_name}
     *  attribute. Returns an empty list for runs that carry no
     *  provenance.
     *
     *  <p>Source of record: the {@code provenance_json} attribute on
     *  the {@code /study/genomic_runs/<name>/} group, written by
     *  {@link global.thalion.ttio.SpectralDataset#writeGenomicRunSubtree}
     *  (Phase 1) — same on-disk pattern as
     *  {@link global.thalion.ttio.AcquisitionRun#writeProvenance}.</p> */
    @Override
    public List<ProvenanceRecord> provenanceChain() {
        return provenanceRecords;
    }

    // ── Streamable<AlignedRead> ────────────────────────────────────

    @Override public boolean hasMore() { return cursor < readCount(); }
    @Override public AlignedRead nextObject() {
        if (!hasMore()) throw new NoSuchElementException();
        return objectAtIndex(cursor++);
    }
    @Override public int currentPosition() { return cursor; }
    @Override public boolean seekToPosition(int position) {
        if (position < 0 || position > readCount()) return false;
        cursor = position;
        return true;
    }
    @Override public void reset() { cursor = 0; }

    @Override
    public void close() {
        dropCachedView();
        if (sequencesDs != null) { sequencesDs.close(); sequencesDs = null; }
        if (qualitiesDs != null) { qualitiesDs.close(); qualitiesDs = null; }
        if (signalChannels != null) { signalChannels.close(); signalChannels = null; }
        sequencesIsV2Cached = null;
        decodedRefDiffV2 = null;
    }

    // ── Internal helpers ───────────────────────────────────────────

    private void ensureSignalChannels() {
        if (signalChannels == null) {
            signalChannels = runGroup.openGroup("signal_channels");
            // v1.8 (ref_diff v2): sequences may be a GROUP rather than a
            // dataset. Probe and leave sequencesDs null when it's a group;
            // the byteChannelSlice dispatch will route via isSequencesRefDiffV2().
            if (!isSequencesRefDiffV2()) {
                sequencesDs = signalChannels.openDataset("sequences");
            }
            qualitiesDs = signalChannels.openDataset("qualities");
        }
    }

    /** Task 13 (ref_diff v2): return {@code true} iff
     *  {@code signal_channels/sequences} is a GROUP containing a
     *  {@code refdiff_v2} child dataset (v1.8 layout).
     *
     *  <p>Uses a try-openGroup pattern (): an
     *  exception from {@code openGroup} means it's a dataset, not a
     *  group. Result is cached in {@link #sequencesIsV2Cached}. */
    private boolean isSequencesRefDiffV2() {
        if (sequencesIsV2Cached != null) return sequencesIsV2Cached;
        // Ensure signal_channels is open before probing.
        if (signalChannels == null) {
            signalChannels = runGroup.openGroup("signal_channels");
        }
        if (!signalChannels.hasChild("sequences")) {
            sequencesIsV2Cached = false;
            return false;
        }
        try (StorageGroup seqGrp = signalChannels.openGroup("sequences")) {
            // It's a group — check for the refdiff_v2 child dataset.
            sequencesIsV2Cached = seqGrp.hasChild("refdiff_v2");
        } catch (Exception e) {
            // sequences is a dataset (v1) — not v2.
            sequencesIsV2Cached = false;
        }
        return sequencesIsV2Cached;
    }

    /** slice {@code count} bytes starting at {@code offset} from a
     *  uint8 byte channel. For codec-compressed channels
     *  ({@code @compression > 0}) the whole channel is decoded once on
     *  first access, the decoded buffer is cached on this
     *  {@link GenomicRun} instance, and subsequent slices come from the
     *  cached array. For uncompressed channels (no attribute or value
     *  0) the existing per-slice {@link StorageDataset#readSlice} path
     *  is used unchanged. */
    private byte[] byteChannelSlice(String name, long offset, int count) {
        byte[] cached = decodedByteChannels.get(name);
        if (cached != null) {
            byte[] out = new byte[count];
            System.arraycopy(cached, (int) offset, out, 0, count);
            return out;
        }
        // v1.8 probe: for sequences, check for the group layout first.
        if ("sequences".equals(name) && isSequencesRefDiffV2()) {
            // Task 4: route REF_DIFF_V2 through the codec registry. The
            // adapter opens the refdiff_v2 child dataset, parses the blob
            // header, resolves the reference via CodecContext, and calls
            // the same RefDiffV2.decode() the old side-path used.
            ensureSignalChannels();
            byte[] decoded;
            try (StorageGroup sigGrp = signalChannels.openGroup("sequences")) {
                decoded = ((global.thalion.ttio.codecs.registry.DecodedChannel.Bytes)
                    global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                        .get(global.thalion.ttio.Enums.Compression.REF_DIFF_V2)
                        .decode(new global.thalion.ttio.codecs.registry.ChannelPayload
                            .GroupPayload(sigGrp), codecContext())).data();
            }
            decodedRefDiffV2 = decoded;
            decodedByteChannels.put(name, decoded);
            byte[] out = new byte[count];
            System.arraycopy(decoded, (int) offset, out, 0, count);
            return out;
        }
        StorageDataset ds = "sequences".equals(name) ? sequencesDs
                          : "qualities".equals(name) ? qualitiesDs
                          : signalChannels.openDataset(name);
        Object codecAttr = ds.getAttribute("compression");
        long codecId = (codecAttr instanceof Number n) ? n.longValue() : 0L;
        if (codecId == 0L) {
            return (byte[]) ds.readSlice(offset, count);
        }
        // Compressed: read the whole channel, decode once, cache.
        long total = ds.shape()[0];
        byte[] all = (byte[]) ds.readSlice(0L, total);
        // Task 4: route through the codec registry. The registry's
        // RANS/BASE_PACK/QUALITY/FQZCOMP decode adapters wrap the same
        // static codec calls verbatim — byte-identical to the old ladder
        // (FQZCOMP pulls revcomp_flags from CodecContext, derived from
        // flags & 16 exactly as before).
        byte[] decoded;
        global.thalion.ttio.Enums.Compression comp = null;
        var vals = global.thalion.ttio.Enums.Compression.values();
        if (codecId >= 0 && codecId < vals.length) comp = vals[(int) codecId];
        var codec = comp == null ? null
            : global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY.get(comp);
        if (codec == null) {
            throw new IllegalStateException(
                "signal_channel '" + name + "': @compression=" + codecId
                + " is not a supported TTIO codec id");
        }
        decoded = ((global.thalion.ttio.codecs.registry.DecodedChannel.Bytes)
            codec.decode(new global.thalion.ttio.codecs.registry.ChannelPayload.BytesPayload(all),
                codecContext())).data();
        decodedByteChannels.put(name, decoded);
        byte[] out = new byte[count];
        System.arraycopy(decoded, (int) offset, out, 0, count);
        return out;
    }

    /** Task 4 (codec registry): build (and cache) the run-derived
     *  {@link global.thalion.ttio.codecs.registry.CodecContext} consumed
     *  by the registry decode adapters. Field derivations are byte-for-byte
     *  copies of the inline logic the pre-registry decode ladder used:
     *  <ul>
     *    <li>{@code revcompFlags} ← {@code (flags & 16) != 0} (FQZCOMP path).</li>
     *    <li>{@code positions}/{@code totalBases}/{@code chromosomes}/
     *        {@code readCount} ← the index (REF_DIFF_V2 path).</li>
     *    <li>{@code ownChromIds} ← the {@code mate_info/chrom_names} sidecar
     *        + index chromosomes, EXACTLY as {@link #_decodeMateV2} derives
     *        them (encounter-order LinkedHashMap, {@code 0xFFFF} sentinel for
     *        chroms absent from the sidecar). When the run has no inline_v2
     *        mate layout the sidecar is absent and {@code ownChromIds} is left
     *        null — the MATE codec is then never invoked.</li>
     *    <li>{@code referenceResolver} ← constructed exactly as the old
     *        decodeRefDiffV2Sequences side-path did (unwrap HDF5 run group →
     *        owning file → new ReferenceResolver). Left null for non-HDF5
     *        backends; the REF_DIFF codec then raises a clear error.</li>
     *  </ul> */
    private global.thalion.ttio.codecs.registry.CodecContext codecContext() {
        if (codecCtxCache != null) return codecCtxCache;
        int n = index().count();
        int[] readLengths = new int[n];
        int[] revcomp = new int[n];
        long[] positions = new long[n];
        long totalBases = 0L;
        for (int i = 0; i < n; i++) {
            readLengths[i] = index().lengthAt(i);
            revcomp[i] = ((index().flagsAt(i) & 16) != 0) ? 1 : 0;
            positions[i] = index().positionAt(i);
            totalBases += index().lengthAt(i);
        }
        String[] chromosomes = new String[n];
        for (int i = 0; i < n; i++) chromosomes[i] = index().chromosomeAt(i);

        // own_chrom_ids: rebuild encounter-order id-per-read from the
        // mate_info/chrom_names sidecar, byte-identical to _decodeMateV2.
        // Only meaningful when an inline_v2 mate layout exists.
        short[] ownChromIds = null;
        List<String> chromTable = null;
        if (isMateInfoInlineV2()) {
            chromTable = readMateInfoChromNamesTable();
            java.util.LinkedHashMap<String, Integer> nameToId =
                new java.util.LinkedHashMap<>();
            for (int j = 0; j < chromTable.size(); j++) {
                nameToId.put(chromTable.get(j), j);
            }
            ownChromIds = new short[n];
            for (int i = 0; i < n; i++) {
                String chr = index().chromosomeAt(i);
                Integer id = nameToId.get(chr);
                ownChromIds[i] = (id == null) ? (short) 0xFFFF
                               : id.shortValue();
            }
        }

        global.thalion.ttio.codecs.ReferenceResolver resolver = resolverForViews();
        final List<String> tagChromTable = chromTable;

        codecCtxCache = global.thalion.ttio.codecs.registry.CodecContext.builder()
            .readLengths(readLengths).revcompFlags(revcomp).readCount(n)
            .positions(positions).totalBases(totalBases).chromosomes(chromosomes)
            .ownChromIds(ownChromIds).ownPositions(positions).nRecords(n)
            .cigarsProvider(() -> allCigars().toArray(new String[0]))
            .referenceResolver(resolver)
            .sequencesProvider(() -> byteChannelFull("sequences"))
            .tagReferencesProvider(() -> tagReferences(tagChromTable, resolver))
            .build();
        return codecCtxCache;
    }

    // decodeRefDiffSequences removed — the v1
    // REF_DIFF reader is no longer supported. Files written with the
    // v1 codec (@compression == 9) raise IllegalStateException at
    // byteChannelSlice (see codec dispatch above).

    // Task 4 (codec registry): decodeRefDiffV2Sequences() removed — the
    // REF_DIFF_V2 sequences side-path in byteChannelSlice now routes
    // through CodecRegistry's RefDiffCodec adapter (which performs the
    // same blob-header parse + ReferenceResolver chain + RefDiffV2.decode),
    // with the run-derived inputs supplied by codecContext().

    /** Return the full list of CIGAR strings for this run. Honours the
     *  M86 Phase C codec dispatch on the cigars channel (RANS /
     *  NAME_TOKENIZED override → uint8 dataset; no override → M82
     *  compound dataset). Caches the result on
     *  {@link #decodedCigars}. */
    private java.util.List<String> allCigars() {
        if (decodedCigars != null) return decodedCigars;
        java.util.List<String> out = new java.util.ArrayList<>(index().count());
        for (int i = 0; i < index().count(); i++) {
            out.add(cigarAt(i));
        }
        // cigarAt() populates decodedCigars when the codec path is hit;
        // the compound-path doesn't set the cache, so set it explicitly
        // here to avoid re-walking on subsequent calls.
        if (decodedCigars == null) decodedCigars = out;
        return decodedCigars;
    }

    /** return the read name at index {@code i}.
     *
     *  <p>Only the NAME_TOKENIZED_V2 (codec id 15) layout is supported
     *  in v1.0+. Legacy v1 layouts raise {@code IllegalStateException}:
     *  <ul>
     *    <li>flat uint8 + {@code @compression == NAME_TOKENIZED (8)} →
     *        v1 codec rejected, see message.</li>
     *    <li>VL-string compound dataset (layout) → also rejected;
     *        the v1.0 writer produces flat uint8 v2 only.</li>
     *  </ul>
     *
     *  <p>If {@code readCount == 0} the writer emits an empty group
     *  (no child datasets); this method short-circuits there. */
    public String readNameAt(int i) {
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).readNameAt(i - (int) blockTable.readStart[b]);
        }
        if (index().count() == 0) {
            // Defensive: read at index 0 on an empty run is an
            // out-of-range error caught upstream; return-empty-string
            // here keeps the codepath safe.
            return "";
        }
        List<String> cached = decodedReadNames;
        if (cached != null) {
            return cached.get(i);
        }
        ensureSignalChannels();
        try (StorageDataset ds = signalChannels.openDataset("read_names")) {
            global.thalion.ttio.Enums.Precision p = ds.precision();
            if (p == global.thalion.ttio.Enums.Precision.UINT8) {
                Object codecAttr = ds.getAttribute("compression");
                long codecId = (codecAttr instanceof Number n)
                    ? n.longValue() : 0L;
                if (codecId == global.thalion.ttio.Enums.Compression
                        .NAME_TOKENIZED_V2.ordinal()) {
                    // v1.8 #11 ch3: name_tok_v2 codec output (NTK2 magic).
                    // Task 4: route through the codec registry (StrList
                    // variant) — the adapter wraps NameTokenizerV2.decode.
                    long total = ds.shape()[0];
                    byte[] all = (byte[]) ds.readSlice(0L, total);
                    decodedReadNames = ((global.thalion.ttio.codecs.registry
                            .DecodedChannel.StrList)
                        global.thalion.ttio.codecs.registry.CodecRegistry
                            .CODEC_REGISTRY.get(global.thalion.ttio.Enums
                                .Compression.NAME_TOKENIZED_V2)
                            .decode(new global.thalion.ttio.codecs.registry
                                .ChannelPayload.BytesPayload(all),
                                codecContext())).names();
                    return decodedReadNames.get(i);
                }
                throw new IllegalStateException(
                    "signal_channel 'read_names': @compression="
                    + codecId + " is not a supported TTIO codec id "
                    + "for the read_names channel (only "
                    + "NAME_TOKENIZED_V2 = 15 is recognised in v1.0+)");
            }
        }
        // Compound (VL_STRING) path was removed in Phase 2c — the
        // v1.0+ writer always emits a flat uint8 dataset (or an empty
        // group for readCount == 0). Files with the M82 compound were
        // produced by older writers; reject with a clear message.
        throw new IllegalStateException(
            "signal_channels/read_names is a compound (VL_STRING) "
            + "dataset — that legacy M82 layout was removed in "
            + "Phase 2c (v1.0 reset). Re-encode the file with v1.0+ "
            + "which writes read_names as NAME_TOKENIZED_V2 (codec "
            + "id 15) on a flat uint8 dataset.");
    }

    /** return the cigar string at index {@code i},
     *  dispatching on dataset shape (Binding Decisions §120-§123).
     *
     *  <p>Two on-disk layouts:
     *  <ul>
     *    <li><b>M82 compound</b> (no override): VL_STRING-in-compound
     *        dataset, read whole-and-cache via {@link #compoundRows}.</li>
     *    <li><b>rANS codec</b> (override active): flat 1-D uint8
     *        dataset. The whole channel is read, decoded once, and
     *        cached as a {@code List<String>} on this instance per
     *        . The decoded byte buffer is a
     *        length-prefix-concat sequence ({@code varint(len) + bytes}
     *        per CIGAR); walk the buffer to reconstruct the list.</li>
     *  </ul>
     */
    public String cigarAt(int i) {
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).cigarAt(i - (int) blockTable.readStart[b]);
        }
        List<String> cached = decodedCigars;
        if (cached != null) {
            return cached.get(i);
        }
        ensureSignalChannels();
        try (StorageDataset ds = signalChannels.openDataset("cigars")) {
            global.thalion.ttio.Enums.Precision p = ds.precision();
            if (p == global.thalion.ttio.Enums.Precision.UINT8) {
                Object codecAttr = ds.getAttribute("compression");
                long codecId = (codecAttr instanceof Number n)
                    ? n.longValue() : 0L;
                long total = ds.shape()[0];
                byte[] all = (byte[]) ds.readSlice(0L, total);
                if (codecId == global.thalion.ttio.Enums.Compression
                        .RANS_ORDER0.ordinal()
                    || codecId == global.thalion.ttio.Enums.Compression
                        .RANS_ORDER1.ordinal()) {
                    byte[] decoded = global.thalion.ttio.codecs
                        .Rans.decode(all);
                    decodedCigars = decodeLengthPrefixConcat(decoded);
                    return decodedCigars.get(i);
                }
                throw new IllegalStateException(
                    "signal_channel 'cigars': @compression="
                    + codecId + " is not a supported TTIO codec id "
                    + "for the cigars channel (only RANS_ORDER0 = 4 "
                    + "and RANS_ORDER1 = 5 are recognised in v1.0+)");
            }
        }
        // Compound path (M82, no override). Materialise the whole
        // list on first call and cache in decodedCigars — without
        // this, per-record cigarAt(i) allocates a fresh String via
        // stringFromCompound per call (UTF-8 decode + char[] alloc),
        // which dominates the genomic transport encode hot path
        // when the source carries no codec override (audit
        // 2026-05-06).
        List<Object[]> rows = compoundRows("cigars");
        int rn = rows.size();
        java.util.List<String> out = new java.util.ArrayList<>(rn);
        for (Object[] row : rows) out.add(stringFromCompound(row[0]));
        decodedCigars = out;
        return out.get(i);
    }

    /** Walk a length-prefix-concat byte buffer back into a list of
     *  ASCII strings. Each entry is {@code varint(len) + len bytes} of
     *  ASCII payload; iteration stops when the buffer is exhausted. */
    private static List<String> decodeLengthPrefixConcat(byte[] buf) {
        List<String> out = new ArrayList<>();
        int offset = 0;
        int n = buf.length;
        long[] tmp = new long[1];
        while (offset < n) {
            offset = readUnsignedVarint(buf, offset, tmp);
            long lengthL = tmp[0];
            if (lengthL < 0 || lengthL > Integer.MAX_VALUE) {
                throw new IllegalArgumentException(
                    "cigars rANS stream: length-prefix-concat entry "
                    + "length " + lengthL + " out of int range");
            }
            int length = (int) lengthL;
            if (offset + length > n) {
                throw new IllegalArgumentException(
                    "cigars rANS stream: length-prefix-concat entry "
                    + "runs off end of decoded buffer (offset="
                    + offset + ", length=" + length
                    + ", buffer_size=" + n + ")");
            }
            for (int k = 0; k < length; k++) {
                int b = Byte.toUnsignedInt(buf[offset + k]);
                if (b > 0x7F) {
                    throw new IllegalArgumentException(
                        "cigars rANS stream: entry contains "
                        + "non-ASCII bytes");
                }
            }
            out.add(new String(buf, offset, length,
                StandardCharsets.US_ASCII));
            offset += length;
        }
        return out;
    }

    /** Unsigned LEB128 varint reader matching the writer in
     *  {@link global.thalion.ttio.SpectralDataset}. Returns the new
     *  offset; the decoded value is stored in {@code out[0]}. */
    private static int readUnsignedVarint(byte[] buf, int offset, long[] out) {
        long value = 0;
        int shift = 0;
        int pos = offset;
        int n = buf.length;
        while (true) {
            if (pos >= n) {
                throw new IllegalArgumentException(
                    "cigars rANS stream: varint runs off end of buffer "
                    + "at offset " + offset);
            }
            int b = Byte.toUnsignedInt(buf[pos]);
            pos++;
            value |= ((long) (b & 0x7F)) << shift;
            if ((b & 0x80) == 0) {
                out[0] = value;
                return pos;
            }
            shift += 7;
            if (shift > 63) {
                throw new IllegalArgumentException(
                    "cigars rANS stream: varint overflow at offset "
                    + offset);
            }
        }
    }

    // ── mate_info per-field dispatch ──────────────────

    /** Task 13 (mate_info v2): true iff
     *  {@code signal_channels/mate_info/inline_v2} exists (v1.7 layout).
     *  Called first in the mate accessor dispatch chain; when true,
     *  {@link #_decodeMateV2()} is used instead of the Phase F subgroup
     *  or M82 compound paths. */
    private boolean isMateInfoInlineV2() {
        if (mateInfoInlineV2Cached != null) return mateInfoInlineV2Cached;
        ensureSignalChannels();
        if (!signalChannels.hasChild("mate_info")) {
            mateInfoInlineV2Cached = Boolean.FALSE;
            return false;
        }
        try (StorageGroup mateGrp = signalChannels.openGroup("mate_info")) {
            mateInfoInlineV2Cached = mateGrp.hasChild("inline_v2");
            return mateInfoInlineV2Cached;
        } catch (Exception e) {
            mateInfoInlineV2Cached = Boolean.FALSE;
            return false;
        }
    }

    /** Task 13: lazily decode the inline_v2 blob + chrom_names sidecar.
     *  Caches the result in {@link #decodedMateV2} and
     *  {@link #mateV2ChromNames}. */
    @SuppressWarnings("unchecked")
    private void _decodeMateV2() {
        if (decodedMateV2 != null) return;
        ensureSignalChannels();
        try (StorageGroup mateGrp = signalChannels.openGroup("mate_info")) {
            // Read the blob.
            byte[] blob;
            try (StorageDataset blobDs = mateGrp.openDataset("inline_v2")) {
                long total = blobDs.shape()[0];
                blob = (byte[]) blobDs.readSlice(0L, total);
            }
            // Read chrom_names sidecar.
            List<Object[]> nameRows;
            try (StorageDataset nameDs = mateGrp.openDataset("chrom_names")) {
                nameRows = (List<Object[]>) nameDs.readAll();
            }
            List<String> chromTable = new ArrayList<>(nameRows.size());
            for (Object[] row : nameRows) {
                Object v = row[0];
                if (v instanceof byte[] b) {
                    chromTable.add(new String(b, StandardCharsets.UTF_8));
                } else {
                    chromTable.add(v == null ? "" : v.toString());
                }
            }
            mateV2ChromNames = chromTable;
            // Task 4: route the inline_v2 decode through the codec
            // registry. The own_chrom_ids / own_positions / nRecords the
            // MATE codec needs come from codecContext(), which derives
            // own_chrom_ids from this same chrom_names sidecar in the same
            // encounter-order LinkedHashMap (0xFFFF sentinel) — byte-
            // identical to the old inline derivation. We rebuild a
            // MateInfoV2.Triple from the registry result to preserve the
            // decodedMateV2 field type + downstream array accessors.
            var mi = (global.thalion.ttio.codecs.registry.DecodedChannel.MateInfo)
                global.thalion.ttio.codecs.registry.CodecRegistry.CODEC_REGISTRY
                    .get(global.thalion.ttio.Enums.Compression.MATE_INLINE_V2)
                    .decode(new global.thalion.ttio.codecs.registry.ChannelPayload
                        .BytesPayload(blob), codecContext());
            decodedMateV2 = new global.thalion.ttio.codecs.MateInfoV2.Triple(
                mi.mateChromIds(), mi.matePositions(), mi.templateLengths());
        }
    }

    // isMateInfoSubgroup removed — the M86
    // Phase F per-field subgroup reader is no longer reached. The
    // mate accessors below short-circuit to inline_v2 or throw.

    /** return the mate chromosome at index
     *  {@code i}. Only the inline_v2 blob layout is supported now;
     *  the M86 Phase F per-field subgroup and the M82 compound layout
     *  raise {@code IllegalStateException}. */
    public String mateChromAt(int i) {
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).mateChromAt(i - (int) blockTable.readStart[b]);
        }
        if (isMateInfoInlineV2()) {
            _decodeMateV2();
            int mateChromId = decodedMateV2.mateChromIds[i];
            if (mateChromId == -1) return "*";
            if (mateV2ChromNames != null && mateChromId < mateV2ChromNames.size()) {
                return mateV2ChromNames.get(mateChromId);
            }
            return "*";  // defensive fallback
        }
        throw mateInfoLegacyLayoutError();
    }

    /** return the mate position at index
     *  {@code i}. Inline_v2 only — see {@link #mateChromAt}. */
    public long matePosAt(int i) {
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).matePosAt(i - (int) blockTable.readStart[b]);
        }
        if (isMateInfoInlineV2()) {
            _decodeMateV2();
            return decodedMateV2.matePositions[i];
        }
        throw mateInfoLegacyLayoutError();
    }

    /** return the template length at index
     *  {@code i}. Inline_v2 only — see {@link #mateChromAt}. */
    public int mateTlenAt(int i) {
        if (blockTable != null) {
            int b = blockTable.blockFor(i);
            return blockView(b).mateTlenAt(i - (int) blockTable.readStart[b]);
        }
        if (isMateInfoInlineV2()) {
            _decodeMateV2();
            return decodedMateV2.templateLengths[i];
        }
        throw mateInfoLegacyLayoutError();
    }

    /** Common error for legacy mate_info layouts (Phase F
     *  per-field subgroup or M82 compound). Both were removed in
     *  Phase 2c; only the v2 inline_v2 blob is read. */
    private static IllegalStateException mateInfoLegacyLayoutError() {
        return new IllegalStateException(
            "signal_channels/mate_info legacy layout (Phase F "
            + "per-field subgroup or M82 compound dataset) is no "
            + "longer supported in v1.0; this file was written with "
            + "an older TTI-O version. Re-encode with v1.0+ which "
            + "uses MATE_INLINE_V2 (codec id 13) at "
            + "signal_channels/mate_info/inline_v2.");
    }

    // decodeMateChrom + decodeMateIntField +
    // decodeLengthPrefixConcatMate removed — the M86 Phase F per-
    // field subgroup readers are gone. Only the v2 inline_v2 blob
    // path survives via _decodeMateV2.

    @SuppressWarnings("unchecked")
    private List<Object[]> compoundRows(String name) {
        List<Object[]> rows = compoundCache.get(name);
        if (rows == null) {
            ensureSignalChannels();
            try (StorageDataset ds = signalChannels.openDataset(name)) {
                rows = (List<Object[]>) ds.readAll();
            }
            compoundCache.put(name, rows);
        }
        return rows;
    }

    private static String stringFromCompound(Object v) {
        if (v == null) return "";
        if (v instanceof byte[] b) return new String(b, StandardCharsets.UTF_8);
        return (String) v;
    }

    private static String stringAttr(StorageGroup g, String name, String fallback) {
        try {
            Object v = g.getAttribute(name);
            if (v instanceof String s) return s;
            if (v instanceof byte[] b) return new String(b, StandardCharsets.UTF_8);
            return v == null ? fallback : v.toString();
        } catch (Exception e) {
            return fallback;
        }
    }

    // ---------------------------------------------------------- Phase 2c-T
    //
    // Verbatim v2 blob accessors used by TransportWriter in bulk mode.
    // Each returns the raw on-disk codec blob bytes (or null when the
    // channel is absent from this run's signal_channels/ group).

    /** Phase 2c-T: read the verbatim {@code mate_info/inline_v2} blob.
     *  Returns null when the run has no inline_v2 layout (empty
     *  mate_chromosomes at write time). */
    @SuppressWarnings("unchecked")
    public byte[] readMateInfoInlineV2BlobBytes() {
        if (blockTable != null) {
            return blockTable.count() == 1 ? blockView(0).readMateInfoInlineV2BlobBytes() : null;
        }
        ensureSignalChannels();
        if (!signalChannels.hasChild("mate_info")) return null;
        try (StorageGroup mateGrp = signalChannels.openGroup("mate_info")) {
            if (!mateGrp.hasChild("inline_v2")) return null;
            try (StorageDataset ds = mateGrp.openDataset("inline_v2")) {
                long total = ds.shape()[0];
                return (byte[]) ds.readSlice(0L, total);
            }
        }
    }

    /** Phase 2c-T: read the {@code mate_info/chrom_names} sidecar
     *  table. Returns an empty list when the table is missing. */
    @SuppressWarnings("unchecked")
    public List<String> readMateInfoChromNamesTable() {
        if (blockTable != null) {
            try (StorageGroup sc = runGroup.openGroup("signal_channels")) {
                return sc.hasChild("mate_info")
                    ? BlockView.readNames(sc.openGroup("mate_info"), "chrom_names")
                    : new ArrayList<>();
            }
        }
        ensureSignalChannels();
        List<String> out = new ArrayList<>();
        if (!signalChannels.hasChild("mate_info")) return out;
        try (StorageGroup mateGrp = signalChannels.openGroup("mate_info")) {
            if (!mateGrp.hasChild("chrom_names")) return out;
            try (StorageDataset nameDs = mateGrp.openDataset("chrom_names")) {
                List<Object[]> rows = (List<Object[]>) nameDs.readAll();
                for (Object[] row : rows) {
                    Object v = row[0];
                    if (v instanceof byte[] b) {
                        out.add(new String(b, StandardCharsets.UTF_8));
                    } else {
                        out.add(v == null ? "" : v.toString());
                    }
                }
            }
        }
        return out;
    }

    /** Phase 2c-T: read the verbatim {@code read_names} blob when
     *  {@code @compression == NAME_TOKENIZED_V2 (15)}. Returns null
     *  when read_names is absent or carries a different codec. */
    public byte[] readNameTokV2BlobBytes() {
        if (blockTable != null) {
            return blockTable.count() == 1 ? blockView(0).readNameTokV2BlobBytes() : null;
        }
        ensureSignalChannels();
        if (!signalChannels.hasChild("read_names")) return null;
        try (StorageDataset ds = signalChannels.openDataset("read_names")) {
            int codec = signalChannelCompressionCode("read_names");
            if (codec != 15) return null;  // not NAME_TOKENIZED_V2
            long total = ds.shape()[0];
            if (total <= 0) return new byte[0];
            return (byte[]) ds.readSlice(0L, total);
        }
    }

    /** Phase 2c-T: read the verbatim
     *  {@code sequences/refdiff_v2} blob when sequences is the v1.8
     *  group layout. Returns null otherwise. */
    public byte[] readRefDiffV2BlobBytes() {
        if (blockTable != null) {
            return blockTable.count() == 1 ? blockView(0).readRefDiffV2BlobBytes() : null;
        }
        ensureSignalChannels();
        if (!signalChannels.hasChild("sequences")) return null;
        // sequences may be a flat dataset or a group containing refdiff_v2.
        try (StorageGroup seqGrp = signalChannels.openGroup("sequences")) {
            if (!seqGrp.hasChild("refdiff_v2")) return null;
            try (StorageDataset ds = seqGrp.openDataset("refdiff_v2")) {
                long total = ds.shape()[0];
                if (total <= 0) return new byte[0];
                return (byte[]) ds.readSlice(0L, total);
            }
        } catch (Exception e) {
            // sequences is a flat dataset, not a group.
            return null;
        }
    }

    // ── Bulk accessors for hot serialization paths ────────────────
    //
    // The per-read accessors (objectAtIndex / readAt /
    // sequence-on-AlignedRead / etc.) are convenient but materialise
    // a fresh AlignedRead and slice every channel on every call.
    // For serialization workloads that touch every byte sequentially
    // (FastqWriter, FastaWriter, TransportWriter at full-corpus
    // scale), pre-fetching the whole channel once and slicing
    // in-memory is dramatically faster — Python's FastqWriter saw a
    // 24× speedup at 1M reads from this exact pattern.

    /** Return the full ``signal_channels/sequences`` byte array.
     *  Side-effect: populates the per-channel cache so that
     *  subsequent per-record {@link #byteChannelSlice} calls slice
     *  the warm buffer instead of issuing fresh HDF5 reads. This
     *  matters for uncompressed channels, which the codec path
     *  does not cache automatically. */
    public byte[] sequencesFull() {
        if (blockTable != null) return concatBlocks(GenomicRun::sequencesFull);
        ensureSignalChannels();
        return byteChannelFull("sequences");
    }

    private byte[] concatBlocks(java.util.function.Function<GenomicRun, byte[]> f) {
        java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
        for (int b = 0; b < blockTable.count(); b++) out.writeBytes(f.apply(blockView(b)));
        return out.toByteArray();
    }

    /** Return the full ``signal_channels/qualities`` byte array.
     *  Same warm-cache semantics as {@link #sequencesFull}. */
    public byte[] qualitiesFull() {
        if (blockTable != null) return concatBlocks(GenomicRun::qualitiesFull);
        ensureSignalChannels();
        return byteChannelFull("qualities");
    }

    /** Cache-priming whole-channel fetch. For codec-compressed
     *  channels {@link #byteChannelSlice} already caches; for
     *  uncompressed channels we explicitly read the full extent and
     *  put it in {@link #decodedByteChannels} so subsequent slices
     *  are O(arraycopy). */
    private byte[] byteChannelFull(String name) {
        byte[] cached = decodedByteChannels.get(name);
        if (cached != null) {
            return cached.clone();
        }
        long total = totalBaseCount();
        byte[] full = byteChannelSlice(name, 0L, (int) total);
        // byteChannelSlice's compressed-codec path already populated
        // the cache. The uncompressed path (codecId == 0) returns the
        // raw HDF5 buffer without caching — fix that here so the
        // warmup is actually effective for those channels too.
        if (!decodedByteChannels.containsKey(name)) {
            decodedByteChannels.put(name, full);
        }
        return full;
    }

    /** Return the full read-names list, forcing the one-shot
     *  NAME_TOKENIZED_V2 decode + cache. Mirrors the Python
     *  ``GenomicRun._read_name_at`` cache priming idiom. */
    public List<String> readNamesAll() {
        if (blockTable != null) {
            List<String> all = new ArrayList<>(readCount());
            for (int b = 0; b < blockTable.count(); b++) all.addAll(blockView(b).readNamesAll());
            return all;
        }
        int n = index().count();
        if (n == 0) return java.util.Collections.emptyList();
        // Touch index 0 to trigger the one-shot decode for v2 layouts.
        readNameAt(0);
        if (decodedReadNames != null) return decodedReadNames;
        // Compound / uncompressed fallback path.
        List<String> out = new ArrayList<>(n);
        for (int i = 0; i < n; i++) out.add(readNameAt(i));
        return out;
    }

    private long totalBaseCount() {
        int n = index().count();
        if (n == 0) return 0L;
        return index().offsetAt(n - 1) + index().lengthAt(n - 1);
    }
}
