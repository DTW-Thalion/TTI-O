/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.codecs;

import java.nio.charset.StandardCharsets;
import java.util.List;

/**
 * Read grouping (M103, {@code ttio_seq_group}).
 *
 * <p>Returns the permutation that puts an unaligned run's reads from the
 * same place in the genome next to each other, so a blocks_v1 block holds
 * them together and SEQ_CM sees each read's overlap partners.
 * {@code order[j]} is the input index of the read stored at row {@code j}.
 * The native kernel breaks every tie on a total order, so the three SDKs
 * compute the same permutation. Delegates through {@link TtioRansNative}.</p>
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.codecs.seq_group}, Objective-C {@code TTIOSeqGroup}.</p>
 */
public final class SeqGroup {

    private SeqGroup() {}

    /** @return {@code true} iff the native JNI library loaded. */
    public static boolean isAvailable() { return TtioRansNative.isAvailable(); }

    /**
     * The grouping order of reads whose bases are {@code sequences} back to
     * back. {@code names} (one per read, passed as stored: the kernel strips
     * a trailing {@code /1} or {@code /2}) pair mates by name; {@code null}
     * groups without them.
     */
    public static int[] group(byte[] sequences, int[] lengths, List<String> names) {
        int n = lengths.length;
        long[] lens = new long[n];
        long total = 0;
        for (int i = 0; i < n; i++) {
            if (lengths[i] < 0) throw new IllegalArgumentException("negative read length at " + i);
            lens[i] = lengths[i];
            total += lengths[i];
        }
        if (total != sequences.length) {
            throw new IllegalArgumentException("lengths sum to " + total + " but there are "
                + sequences.length + " bases");
        }
        if (names == null) return TtioRansNative.groupReads(sequences, lens, null, null);
        if (names.size() != n) {
            throw new IllegalArgumentException(names.size() + " names for " + n + " reads");
        }
        byte[][] enc = new byte[n][];
        long[] off = new long[n + 1];
        for (int i = 0; i < n; i++) {
            String s = names.get(i);
            enc[i] = s == null ? new byte[0] : s.getBytes(StandardCharsets.UTF_8);
            off[i + 1] = off[i] + enc[i].length;
        }
        if (off[n] > Integer.MAX_VALUE - 8) {
            throw new IllegalArgumentException("read names too large for one grouping call");
        }
        byte[] buf = new byte[(int) off[n]];
        for (int i = 0; i < n; i++) System.arraycopy(enc[i], 0, buf, (int) off[i], enc[i].length);
        return TtioRansNative.groupReads(sequences, lens, buf, off);
    }
}
