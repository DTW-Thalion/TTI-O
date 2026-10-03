/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio;

import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.codecs.SeqCm;
import global.thalion.ttio.codecs.registry.ChannelPayload;
import global.thalion.ttio.codecs.registry.Codec;
import global.thalion.ttio.codecs.registry.CodecContext;
import global.thalion.ttio.codecs.registry.CodecRegistry;
import global.thalion.ttio.codecs.registry.DecodedChannel;
import global.thalion.ttio.codecs.registry.EncodedChannel;

import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.Random;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * M103: SEQ_CM (codec id 19), read bases without a reference.
 *
 * <p>The wrapper round-trips every byte value and zero-length reads,
 * codes the same reads to the same bytes, rejects mismatched lengths,
 * and the registry adapter takes the lengths from the codec context.
 * The blocks_v1 default is covered in {@code GenomicBlocksTest}.</p>
 *
 * <p>Mirrors {@code python/tests/test_m103_seq_cm.py} and
 * {@code objc/Tests/TestM103SeqCm.m}.</p>
 */
final class M103SeqCmTest {

    private byte[] seq;
    private int[] lengths;

    @BeforeEach
    void reads() {
        Assumptions.assumeTrue(SeqCm.isAvailable(), "libttio_rans_jni not loaded");
        Random rng = new Random(3);
        byte[] genome = new byte[20000];
        byte[] acgt = {'A', 'C', 'G', 'T'};
        for (int i = 0; i < genome.length; i++) genome[i] = acgt[rng.nextInt(4)];
        int n = 400;
        lengths = new int[n];
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (int i = 0; i < n; i++) {
            int len = i == 3 ? 0 : 40 + rng.nextInt(120);
            int s = rng.nextInt(genome.length - 200);
            byte[] r = Arrays.copyOfRange(genome, s, s + len);
            if (len > 0 && i % 7 == 0) r[len / 2] = 'N';
            if (len > 0 && i % 11 == 0) r[0] = 'a';
            out.write(r, 0, r.length);
            if (i == n - 1) {
                for (int b = 0; b < 256; b++) out.write(b);
                len += 256;
            }
            lengths[i] = len;
        }
        seq = out.toByteArray();
    }

    @Test
    void roundTripEveryByteAndEmptyRead() {
        byte[] blob = SeqCm.encode(seq, lengths);
        assertArrayEquals(new byte[] {'S', 'Q', 'C', '1'}, Arrays.copyOf(blob, 4));
        assertArrayEquals(seq, SeqCm.decode(blob, lengths));
        assertArrayEquals(blob, SeqCm.encode(seq, lengths), "deterministic");
    }

    @Test
    void emptyBlock() {
        byte[] blob = SeqCm.encode(new byte[0], new int[0]);
        assertEquals(0, SeqCm.decode(blob, new int[0]).length);
    }

    @Test
    void lengthsMustMatch() {
        int[] shortLens = Arrays.copyOf(lengths, lengths.length - 1);
        assertThrows(IllegalArgumentException.class, () -> SeqCm.encode(seq, shortLens));
        byte[] blob = SeqCm.encode(seq, lengths);
        assertThrows(IllegalArgumentException.class, () -> SeqCm.decode(blob, shortLens));
    }

    @Test
    void registryTakesLengthsFromContext() {
        Codec c = CodecRegistry.CODEC_REGISTRY.get(Compression.SEQ_CM);
        assertNotNull(c);
        assertTrue(c.isContextAware());
        CodecContext ctx = CodecContext.builder().readLengths(lengths).build();
        byte[] enc = ((EncodedChannel.DatasetBytes) c.encode(new DecodedChannel.Bytes(seq), ctx)).bytes();
        assertArrayEquals(SeqCm.encode(seq, lengths), enc);
        byte[] dec = ((DecodedChannel.Bytes) c.decode(new ChannelPayload.BytesPayload(enc), ctx)).data();
        assertArrayEquals(seq, dec);
    }
}
