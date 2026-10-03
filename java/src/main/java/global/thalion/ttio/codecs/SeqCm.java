/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.codecs;

/**
 * SEQ_CM — read bases without a reference (codec id 19, M103).
 *
 * <p>A context-mixing model (orders 11/16/24 with reverse-complement
 * training) and a binary arithmetic coder. The blob is written as
 * {@code signal_channels/sequences} with {@code @compression = 19}.
 * Spec: {@code docs/codecs/seq_cm.md}.</p>
 *
 * <p>Context-aware: decode takes the reads' lengths, which the run
 * already stores. Only the default parameters are used and the kernel
 * picks the model's table size from the base count, so the same reads
 * always code to the same bytes. Delegates to
 * {@code ttio_seq_cm_encode/_decode} through {@link TtioRansNative}.</p>
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.codecs.seq_cm}, Objective-C {@code TTIOSeqCm}.</p>
 */
public final class SeqCm {

    private SeqCm() {}

    /** @return {@code true} iff the native JNI library loaded. */
    public static boolean isAvailable() { return TtioRansNative.isAvailable(); }

    /** Encode reads' bases, back to back, given one length per read. */
    public static byte[] encode(byte[] sequences, int[] lengths) {
        return TtioRansNative.encodeSeqCm(sequences, widen(lengths));
    }

    /** Decode a SEQ_CM blob given the reads' lengths. */
    public static byte[] decode(byte[] blob, int[] lengths) {
        return TtioRansNative.decodeSeqCm(blob, widen(lengths));
    }

    private static long[] widen(int[] lengths) {
        long[] out = new long[lengths.length];
        for (int i = 0; i < lengths.length; i++) {
            if (lengths[i] < 0) throw new IllegalArgumentException("negative read length at " + i);
            out[i] = lengths[i];
        }
        return out;
    }
}
