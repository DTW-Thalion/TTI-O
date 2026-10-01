/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.codecs;

/**
 * JNI bridge to libttio_rans native library.
 *
 * <p>Thread-safe. Methods marshal Java arrays to C, call the rANS kernels
 * (which dispatch to AVX2/SSE4.1/scalar at runtime), and return results.
 *
 * <p>Falls back gracefully: if the native library cannot be loaded
 * ({@link UnsatisfiedLinkError} on first call to {@code System.loadLibrary}),
 * {@link #isAvailable()} returns false and callers should use the pure-Java
 * implementation in {@link FqzcompNx16Z}.
 *
 * <p>This class is the Java analogue of the Python ctypes wrapper in
 * {@code python/src/ttio/codecs/fqzcomp_nx16_z.py}. The JNI C glue lives at
 * {@code native/src/ttio_rans_jni.c}; build it with
 * {@code cmake -DTTIO_RANS_BUILD_JNI=ON}.
 */
public final class TtioRansNative {
    private static final boolean LOADED;

    static {
        boolean ok = false;
        try {
            System.loadLibrary("ttio_rans_jni");
            ok = true;
        } catch (UnsatisfiedLinkError e) {
            // Library not on java.library.path; fall back to pure Java.
            ok = false;
        }
        LOADED = ok;
    }

    private TtioRansNative() {}

    /**
     * @return {@code true} if the native library was loaded successfully and
     *         the {@code encodeBlock}/{@code decodeBlock}/{@code kernelName}
     *         native methods are callable.
     */
    public static boolean isAvailable() { return LOADED; }

    /**
     * Encode a block of symbols via libttio_rans.
     *
     * @param symbols    raw symbols (0..255 each, treated as unsigned)
     * @param contexts   parallel context vector (one uint16 per symbol)
     * @param nContexts  total number of contexts (= freq.length)
     * @param freq       freq table, shape [nContexts][256]; entries must sum
     *                   to T = 4096 per context with no zero-frequency symbols
     *                   actually used
     * @param out        output buffer (caller-allocated, sized worst-case)
     * @param outLen     in/out single-element array: input is buffer capacity,
     *                   output is bytes actually written
     * @return TTIO_RANS_OK (0) or negative error code (TTIO_RANS_ERR_PARAM=-1,
     *         ERR_ALLOC=-2, ERR_CORRUPT=-3)
     */
    public static native int encodeBlock(
        byte[] symbols, short[] contexts, int nContexts,
        int[][] freq, byte[] out, int[] outLen);

    /**
     * Decode a block of symbols via libttio_rans.
     *
     * @param compressed encoded bytes (exactly the slice produced by encode;
     *                   no trailing slack)
     * @param contexts   parallel context vector (one uint16 per output symbol)
     * @param nContexts  total number of contexts
     * @param freq       freq table, shape [nContexts][256]
     * @param cum        cumulative table, shape [nContexts][256], where
     *                   cum[c][s] = sum_{s'<s} freq[c][s']
     * @param symbols    output buffer, length nSymbols, written in place
     * @param nSymbols   number of symbols to decode
     * @return TTIO_RANS_OK (0) or negative error code
     */
    public static native int decodeBlock(
        byte[] compressed, short[] contexts, int nContexts,
        int[][] freq, int[][] cum,
        byte[] symbols, int nSymbols);

    /**
     * @return Selected SIMD kernel name: "scalar", "sse4.1", or "avx2".
     *         Determined at library-load time by CPUID dispatch in
     *         {@code dispatch.c}.
     */
    public static native String kernelName();

    /** Threads the FQZCOMP auto-tune uses for its candidate encodes
     *  (default 3; {@code n <= 1} runs them in sequence). */
    public static native void setAutotuneThreads(int n);

    /** Threads M94.Z V6 uses for the segments of one block. 0 defers to
     *  {@link #setAutotuneThreads}, which is what V6 did before it had a
     *  count of its own — and which meant a writer pool, correctly
     *  standing the auto-tune down to 1, also switched V6's segments to
     *  running in sequence. */
    public static native void setV6Threads(int n);

    public static native void setV6Sbits(int n);

    public static native int getV6Sbits();

    /** @return the V6 segment thread count; 0 means defer to the
     *  auto-tune knob. */
    public static native int getV6Threads();
    public static native int getAutotuneThreads();

    /** Strategy of an encoded M94.Z stream: 4 = V4, 5/6 = V5 S5/S6;
     *  negative on a stream that is not M94.Z qualities. */
    public static native int qualStreamStrategy(byte[] stream);

    /**
     * Per-symbol context resolver used by {@link #decodeBlockStreaming}.
     *
     * <p>Implementations must derive a context id from the current symbol
     * index {@code i} and the previously decoded symbol {@code prevSym}
     * (0 when {@code i == 0}). Must be deterministic. The returned int
     * is treated as an unsigned 16-bit context id; it must be less than
     * the {@code nContexts} passed to {@link #decodeBlockStreaming}, or
     * the native call returns {@code TTIO_RANS_ERR_PARAM}.
     *
     * <p><b>Performance reality (Task 26c)</b>: per-symbol JNI dispatch
     * has a JVM context-switch cost that is much higher than ctypes
     * callback dispatch in Python (which itself was a wash with
     * pure-Python in Task 26b). For realistic block sizes the streaming
     * path is expected to be SLOWER than the pure-Java V2 decoder. The
     * binding is shipped as infrastructure; codec callers gate it
     * behind an opt-in path with a safe pure-Java fallback.
     */
    @FunctionalInterface
    public interface ContextResolver {
        /**
         * @param i        current symbol index (0-based, in [0, nSymbols))
         * @param prevSym  previously decoded symbol (or 0 for i==0); the
         *                 byte is passed as an unsigned int via JNI
         * @return         context id, treated as unsigned 16-bit
         */
        int resolve(long i, int prevSym);
    }

    /**
     * Decode a block via libttio_rans's streaming context API.
     *
     * <p>The native library calls back into {@code resolver.resolve(i, prevSym)}
     * before decoding each non-padding symbol, so codecs whose context
     * derives from previously-decoded symbols (e.g. M94.Z order-1
     * cascades) can drive the C decode kernel without materialising the
     * contexts vector up front. Padding positions use {@code ctx=0}
     * (matching the encoder) and skip the callback.
     *
     * <p>See {@link ContextResolver} for the per-call cost reality. Callers
     * should treat this as an opt-in path with a safe pure-Java fallback
     * on {@code rc != 0} or any thrown exception.
     *
     * @param compressed encoded body bytes (V2 native body — same layout
     *                   as {@link #decodeBlock} input)
     * @param nContexts  number of contexts (= freq.length = cum.length)
     * @param freq       freq table, shape [nContexts][256]
     * @param cum        cumulative table, shape [nContexts][256]
     * @param symbols    output buffer, length >= nSymbols
     * @param nSymbols   number of symbols to decode
     * @param resolver   per-symbol context derivation callback
     * @return TTIO_RANS_OK (0) or negative error code
     */
    public static native int decodeBlockStreaming(
        byte[] compressed, int nContexts,
        int[][] freq, int[][] cum,
        byte[] symbols, int nSymbols,
        ContextResolver resolver);

    /**
     * Encode flat qualities to an M94.Z V4 stream via the native library.
     *
     * <p>Mirrors {@code ttio_m94z_v4_encode} in {@code native/include/ttio_rans.h}.
     * Auto-tunes (CRAM 3.1) when {@code strategyHint == -1}; uses the
     * specified preset 0..3 otherwise.
     *
     * @param qualities flat Phred bytes (length == sum of readLengths)
     * @param readLengths per-read quality counts
     * @param flags per-read SAM flags (Phase 2 strategy 1 ignores these; auto-tune
     *              uses bit 4 = SAM_REVERSE_FLAG)
     * @param strategyHint -1 = auto-tune (default); 0..3 = explicit preset
     * @param padCount 0..3 (low-2 bits of the V4 flags byte)
     * @return encoded V4 stream
     * @throws RuntimeException if the native call returns a non-zero rc
     */
    public static byte[] encodeV4(byte[] qualities, int[] readLengths, int[] flags,
                                   int strategyHint, int padCount) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return encodeV4Native(qualities, readLengths, flags, strategyHint, padCount);
    }

    /**
     * Decode an M94.Z V4 stream via the native library.
     *
     * @param encoded V4 stream (must start with "M94Z" + version 4)
     * @param numReads expected read count (decoder pre-allocates lengths array)
     * @param numQualities expected quality count
     * @param flags per-read SAM flags
     * @return [qualities[], readLengths[]] as a 2-element Object[]; element 0 is
     *         the byte[] of qualities, element 1 is the int[] of recovered
     *         read lengths
     */
    public static Object[] decodeV4(byte[] encoded, int numReads, int numQualities,
                                     int[] flags) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return decodeV4Native(encoded, numReads, numQualities, flags);
    }

    private static native byte[] encodeV4Native(byte[] qualities, int[] readLengths,
                                                  int[] flags, int strategyHint,
                                                  int padCount);

    private static native Object[] decodeV4Native(byte[] encoded, int numReads,
                                                    int numQualities, int[] flags);

    /**
     * Encode flat qualities through the M94.Z umbrella (V4 presets plus
     * the S5/S6 sequence-context strategies when {@code sequences} is
     * non-null; smallest stream wins by exact size).
     */
    public static byte[] encodeQual(byte[] qualities, int[] readLengths,
                                    int[] flags, byte[] sequences,
                                    int strategyHint, int padCount) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return encodeQualNative(qualities, readLengths, flags, sequences,
                                strategyHint, padCount);
    }

    /**
     * Decode an M94.Z V4 or V5 stream. {@code sequences} is required for
     * V5 streams (the decoded sequences channel) and ignored for V4.
     */
    public static Object[] decodeQual(byte[] encoded, int numReads,
                                      int numQualities, int[] flags,
                                      byte[] sequences) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return decodeQualNative(encoded, numReads, numQualities, flags,
                                sequences);
    }

    private static native byte[] encodeQualNative(byte[] qualities,
                                                    int[] readLengths,
                                                    int[] flags,
                                                    byte[] sequences,
                                                    int strategyHint,
                                                    int padCount);

    private static native Object[] decodeQualNative(byte[] encoded,
                                                      int numReads,
                                                      int numQualities,
                                                      int[] flags,
                                                      byte[] sequences);

    /**
     * Encode a mate triple via libttio_rans (mate_info v2).
     *
     * @param mateChromIds    int[N] — -1 if RNEXT='*', else id >= 0
     * @param matePositions   long[N] — 0-based POS
     * @param templateLengths int[N] — signed tlen
     * @param ownChromIds     short[N] — uint16 from L1 chrom_ids;
     *                        (short)0xFFFF treated as unmapped sentinel
     * @param ownPositions    long[N] — own POS
     * @return encoded inline_v2 blob ready for HDF5 write
     * @throws RuntimeException on native error
     */
    public static byte[] encodeMateInfoV2(
            int[] mateChromIds, long[] matePositions, int[] templateLengths,
            short[] ownChromIds, long[] ownPositions) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return encodeMateInfoV2Native(mateChromIds, matePositions, templateLengths,
                                       ownChromIds, ownPositions);
    }

    /**
     * Decode an inline_v2 blob via libttio_rans.
     *
     * @return Object[3]: int[] mateChromIds, long[] matePositions, int[] templateLengths
     * @throws RuntimeException on native error
     */
    public static Object[] decodeMateInfoV2(
            byte[] encoded, short[] ownChromIds, long[] ownPositions, int nRecords) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return decodeMateInfoV2Native(encoded, ownChromIds, ownPositions, nRecords);
    }

    private static native byte[] encodeMateInfoV2Native(
        int[] mateChromIds, long[] matePositions, int[] templateLengths,
        short[] ownChromIds, long[] ownPositions);

    private static native Object[] decodeMateInfoV2Native(
        byte[] encoded, short[] ownChromIds, long[] ownPositions, int nRecords);

    /**
     * Encode a slice of reads via libttio_rans (ref_diff v2).
     */
    public static byte[] encodeRefDiffV2(
            byte[]   sequences,
            long[]   offsets,
            long[]   positions,
            String[] cigarStrings,
            byte[]   reference,
            byte[]   referenceMd5,
            String   referenceUri,
            int      readsPerSlice) {
        return encodeRefDiffV2(sequences, offsets, positions, cigarStrings,
                               reference, referenceMd5, referenceUri,
                               readsPerSlice, 0L);
    }

    /**
     * Encode with a byte budget on the slice partition (M97): with
     * {@code sliceBytes} > 0 a slice closes before the read that would
     * push it past that many bases; 0 keeps the fixed-count rule.
     * Writer policy only — the wire format and decoder are unchanged.
     */
    public static byte[] encodeRefDiffV2(
            byte[]   sequences,
            long[]   offsets,
            long[]   positions,
            String[] cigarStrings,
            byte[]   reference,
            byte[]   referenceMd5,
            String   referenceUri,
            int      readsPerSlice,
            long     sliceBytes) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return encodeRefDiffV2Native(sequences, offsets, positions,
                                      cigarStrings, reference, referenceMd5,
                                      referenceUri, readsPerSlice, sliceBytes);
    }

    /**
     * Decode a refdiff_v2 blob.
     * @return Object[2]: byte[] sequences, long[] offsets
     */
    public static Object[] decodeRefDiffV2(
            byte[]   encoded,
            long[]   positions,
            String[] cigarStrings,
            byte[]   reference,
            int      nReads,
            long     totalBases) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return decodeRefDiffV2Native(encoded, positions, cigarStrings,
                                      reference, nReads, totalBases);
    }

    private static native byte[] encodeRefDiffV2Native(
        byte[] sequences, long[] offsets, long[] positions,
        String[] cigarStrings, byte[] reference, byte[] referenceMd5,
        String referenceUri, int readsPerSlice, long sliceBytes);

    private static native Object[] decodeRefDiffV2Native(
        byte[] encoded, long[] positions, String[] cigarStrings,
        byte[] reference, int nReads, long totalBases);

    /**
     * Encode a list of read names to a NAME_TOKENIZED v2 (codec id 15) blob.
     *
     * @param names ASCII read names; UTF-8 marshalling matches the C entry's
     *              {@code const char * const *} contract
     * @return encoded blob (starts with magic {@code "NTK2"})
     * @throws RuntimeException on native error (rc != 0)
     */
    public static native byte[] encodeNameTokV2Native(String[] names);

    /**
     * Decode a NAME_TOKENIZED v2 (codec id 15) blob to read names.
     *
     * @param blob encoded blob produced by {@link #encodeNameTokV2Native}
     * @return recovered read names in original order
     * @throws RuntimeException on native error (rc != 0)
     */
    public static native String[] decodeNameTokV2Native(byte[] blob);

    /**
     * Encode per-read SAM tag text to a SAM_TAGS (codec id 18, M101)
     * blob. Text arrays are UTF-8 bytes plus {@code n + 1} offsets;
     * {@code refs} null (or all-null) disables MD/NM derivation and the
     * other context arrays may then be null. See {@link SamTags}.
     */
    public static byte[] encodeSamTags(
            byte[] tags, long[] tagOffsets,
            byte[] sequences, long[] seqOffsets, byte[] cigars, long[] cigarOffsets,
            long[] positions, short[] chromIds, byte[][] refs) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return encodeSamTagsNative(tags, tagOffsets, sequences, seqOffsets,
                                   cigars, cigarOffsets, positions, chromIds, refs);
    }

    /**
     * Decode a SAM_TAGS blob (same context as {@link #encodeSamTags}).
     * @return Object[2]: byte[] UTF-8 text, long[] offsets (n + 1)
     */
    public static Object[] decodeSamTags(
            byte[] encoded, int nReads,
            byte[] sequences, long[] seqOffsets, byte[] cigars, long[] cigarOffsets,
            long[] positions, short[] chromIds, byte[][] refs) {
        if (!LOADED) throw new IllegalStateException("libttio_rans_jni not loaded");
        return decodeSamTagsNative(encoded, nReads, sequences, seqOffsets,
                                   cigars, cigarOffsets, positions, chromIds, refs);
    }

    private static native byte[] encodeSamTagsNative(
        byte[] tags, long[] tagOffsets,
        byte[] sequences, long[] seqOffsets, byte[] cigars, long[] cigarOffsets,
        long[] positions, short[] chromIds, byte[][] refs);

    private static native Object[] decodeSamTagsNative(
        byte[] encoded, int nReads,
        byte[] sequences, long[] seqOffsets, byte[] cigars, long[] cigarOffsets,
        long[] positions, short[] chromIds, byte[][] refs);
}
