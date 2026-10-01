/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.importers;

import global.thalion.ttio.Enums.AcquisitionMode;
import global.thalion.ttio.Enums.Compression;
import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.genomics.WrittenGenomicRun;

import java.util.ArrayList;
import java.util.List;

/** Per-read accumulator of the BAM/SAM/CRAM importers: collects records
 *  and packs them into a {@link WrittenGenomicRun} (whole file or one
 *  batch of a stream). */
final class BatchAccumulator {

    private final List<String> readNames = new ArrayList<>();
    private final List<String> chromosomes = new ArrayList<>();
    private final List<Long> positions = new ArrayList<>();
    private final List<Integer> mappingQualities = new ArrayList<>();
    private final List<Integer> flags = new ArrayList<>();
    private final List<String> cigars = new ArrayList<>();
    private final List<String> mateChromosomes = new ArrayList<>();
    private final List<Long> matePositions = new ArrayList<>();
    private final List<Integer> templateLengths = new ArrayList<>();
    private final List<byte[]> seqChunks = new ArrayList<>();
    private final List<byte[]> qualChunks = new ArrayList<>();
    /** M101: per-read SAM tag text ("" when a read has none). */
    private final List<String> tags = new ArrayList<>();
    private long totalBases;

    int size() { return readNames.size(); }

    void add(String qname, int flag, String rname, long pos, int mapq, String cigar,
             String rnext, long pnext, int tlen, byte[] seq, byte[] qual) {
        add(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, "");
    }

    void add(String qname, int flag, String rname, long pos, int mapq, String cigar,
             String rnext, long pnext, int tlen, byte[] seq, byte[] qual, String tagText) {
        tags.add(tagText == null ? "" : tagText);
        readNames.add(qname);
        flags.add(flag);
        chromosomes.add(rname);
        positions.add(pos);
        mappingQualities.add(mapq);
        cigars.add(cigar);
        mateChromosomes.add(rnext);
        matePositions.add(pnext);
        templateLengths.add(tlen);
        seqChunks.add(seq);
        qualChunks.add(qual);
        totalBases += seq.length;
    }

    void clear() {
        readNames.clear(); chromosomes.clear(); positions.clear(); mappingQualities.clear();
        flags.clear(); cigars.clear(); mateChromosomes.clear(); matePositions.clear();
        templateLengths.clear(); seqChunks.clear(); qualChunks.clear(); tags.clear();
        totalBases = 0;
    }

    WrittenGenomicRun toRun(AcquisitionMode mode, String referenceUri, String platform,
                            String sampleName, List<ProvenanceRecord> provenance) {
        int n = readNames.size();
        long[] pos = new long[n];
        byte[] mapq = new byte[n];
        int[] flg = new int[n];
        long[] offsets = new long[n];
        int[] lengths = new int[n];
        long[] matePos = new long[n];
        int[] tlens = new int[n];
        int totalQual = 0;
        for (int i = 0; i < n; i++) {
            pos[i] = positions.get(i);
            mapq[i] = (byte) (mappingQualities.get(i) & 0xFF);
            flg[i] = flags.get(i);
            lengths[i] = seqChunks.get(i).length;
            matePos[i] = matePositions.get(i);
            tlens[i] = templateLengths.get(i);
            totalQual += qualChunks.get(i).length;
        }
        byte[] sequences = new byte[(int) totalBases];
        byte[] qualities = new byte[totalQual];
        int so = 0, qo = 0;
        for (int i = 0; i < n; i++) {
            offsets[i] = so;
            byte[] s = seqChunks.get(i), q = qualChunks.get(i);
            System.arraycopy(s, 0, sequences, so, s.length);
            so += s.length;
            System.arraycopy(q, 0, qualities, qo, q.length);
            qo += q.length;
        }
        return new WrittenGenomicRun(mode, referenceUri, platform, sampleName,
            pos, mapq, flg, sequences, qualities, offsets, lengths,
            new ArrayList<>(cigars), new ArrayList<>(readNames), new ArrayList<>(mateChromosomes),
            matePos, tlens, new ArrayList<>(chromosomes), Compression.ZLIB, java.util.Map.of(),
            provenance == null ? List.of() : List.copyOf(provenance),
            false, null, null, null, false, false, null, 0L, new ArrayList<>(tags));
    }
}
