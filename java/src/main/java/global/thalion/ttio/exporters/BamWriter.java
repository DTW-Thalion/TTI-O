/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.exporters;

import global.thalion.ttio.ProvenanceRecord;
import global.thalion.ttio.genomics.AlignedRead;
import global.thalion.ttio.genomics.GenomicRun;
import global.thalion.ttio.genomics.WrittenGenomicRun;
import global.thalion.ttio.importers.BamReader;
import global.thalion.ttio.io.ProgressSink;

import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.cram.ref.ReferenceSource;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/**
 * BAM exporter — M88.
 *
 * <p>v1.5.0 swapped the previous {@code samtools} subprocess approach for
 * <a href="https://github.com/samtools/htsjdk">htsjdk</a>'s pure-Java
 * SAM/BAM/CRAM writers (used by GATK, Picard, IGV). No external binary
 * required; tio-browser works on Windows without a system samtools
 * install. SAM line semantics still match the SAMv1 specification.</p>
 *
 * <h2>Quality byte encoding</h2>
 * <p>M87's {@link BamReader} stores SAM's QUAL field bytes verbatim
 * into {@link WrittenGenomicRun#qualities()} — i.e. the buffer holds
 * <b>ASCII Phred+33</b> characters (so a Phred-40 score is stored as
 * the byte value 73, the ASCII code for {@code 'I'}). This writer
 * mirrors that convention: each {@code qualities[i]} byte is
 * interpreted as ASCII Phred+33 on write, then converted to raw Phred
 * for htsjdk's {@code SAMRecord.setBaseQualities(byte[])} (which expects
 * raw bytes). The all-{@code 0xFF} sentinel produced by the reader when
 * source SAM had {@code QUAL '*'} maps back to "no qualities" via
 * {@link SAMRecord#NULL_QUALS}.</p>
 *
 * <p><b>Cross-language equivalents:</b> Python
 * {@code ttio.exporters.bam.BamWriter},
 * Objective-C {@code TTIOBamWriter}. Cross-language byte-equality on
 * BAM output is verified by the {@code cross-compat} CI job.</p>
 *
 * (M88, v1.5.0 htsjdk swap)
 */
public class BamWriter {

    /**
     * Default {@code @SQ LN:} length when the writer doesn't know the
     * true reference length. SAM requires {@code LN:} on every
     * {@code @SQ}; we pick INT32_MAX so the emitted header is valid
     * for any plausible coordinate. Matches the Python reference's
     * {@code _DEFAULT_SQ_LENGTH} and the ObjC writer's constant for
     * cross-language byte-equality.
     */
    protected static final long DEFAULT_SQ_LENGTH = 2147483647L;

    /** Stage D: emit-every-N cadence for {@link ProgressSink}
     *  callbacks during BAM/CRAM record serialisation. Mirrors
     *  {@code BamReader.PROGRESS_INTERVAL_READS = 1000}. */
    public static final int PROGRESS_INTERVAL_READS = 1000;

    private final Path path;

    /**
     * Construct a {@code BamWriter}.
     *
     * @param path output BAM file path. htsjdk uses the file extension
     *             ({@code .bam}, {@code .sam}, {@code .cram}) to pick
     *             the output format.
     */
    public BamWriter(Path path) {
        this.path = Objects.requireNonNull(path);
    }

    /** @return the output path passed at construction time. */
    public Path path() { return path; }

    // ------------------------------------------------------------------
    // Public API
    // ------------------------------------------------------------------

    /**
     * Serialise {@code run} to the configured output path.
     *
     * @param run        the genomic-run container to write.
     * @param provenance optional provenance records to inject as
     *                   {@code @PG} header lines. Pass an empty list
     *                   for none.
     * @param sort       when {@code true}, set {@code @HD SO:coordinate}
     *                   so htsjdk's writer sorts records on close
     *                   (precondition most downstream tools expect —
     *                   IGV, GATK, {@code samtools index}). When
     *                   {@code false}, output is unsorted.
     * @throws IOException on write failures.
     */
    public void write(WrittenGenomicRun run, List<ProvenanceRecord> provenance,
                      boolean sort) throws IOException {
        write(run, provenance, sort, ProgressSink.discard());
    }

    /**
     * Stage D overload of
     * {@link #write(WrittenGenomicRun, List, boolean)} that fires
     * {@code progress.onProgress(readsDone, totalReads)} every
     * {@link #PROGRESS_INTERVAL_READS} records during the alignment
     * write loop and a final {@code onProgress(total, total)} once
     * the writer closes.
     *
     * <p>{@code total} is known up front (= {@code run.readNames().size()})
     * so every mid-emit callback already passes the determinate total —
     * listeners can drive a determinate progress bar from the first
     * fire.</p>
     *
     * @since 1.5.0
     */
    public void write(WrittenGenomicRun run, List<ProvenanceRecord> provenance,
                      boolean sort, ProgressSink progress) throws IOException {
        Objects.requireNonNull(run, "run");
        if (provenance == null) provenance = List.of();
        if (progress == null) progress = ProgressSink.discard();
        SAMFileHeader header = buildHeader(run, provenance, sort);
        try (SAMFileWriter writer = makeWriter(header, sort)) {
            int n = run.readNames().size();
            long total = n;
            for (int i = 0; i < n; i++) {
                SAMRecord rec = buildSamRecord(run, header, i);
                writer.addAlignment(rec);
                long done = (long) (i + 1);
                if (done % PROGRESS_INTERVAL_READS == 0 && done < total) {
                    progress.onProgress(done, total);
                }
            }
            // Final fire always — even for empty / sub-cadence runs.
            progress.onProgress(total, total);
        }
    }

    /** Stream a stored {@link GenomicRun} out read by read through
     *  {@link GenomicRun#iterReads()}: one decoded block resident for a
     *  {@code blocks_v1} run. The header comes from the run-level
     *  chromosome table and attributes. */
    public void write(GenomicRun run, List<ProvenanceRecord> provenance,
                      boolean sort, ProgressSink progress) throws IOException {
        Objects.requireNonNull(run, "run");
        if (provenance == null) provenance = List.of();
        if (progress == null) progress = ProgressSink.discard();
        SAMFileHeader header = buildHeader(run.chromosomeNames(), run.sampleName(),
            run.platform(), provenance, sort);
        long total = run.readCount();
        try (SAMFileWriter writer = makeWriter(header, sort)) {
            java.util.Iterator<AlignedRead> it = run.iterReads(0, run.readCount(), global.thalion.ttio.Threads.resolve(null));
            long done = 0;
            while (it.hasNext()) {
                writer.addAlignment(buildSamRecord(it.next(), header));
                done++;
                if (done % PROGRESS_INTERVAL_READS == 0 && done < total) {
                    progress.onProgress(done, total);
                }
            }
            progress.onProgress(total, total);
        }
    }

    /** {@link #write(GenomicRun, List, boolean, ProgressSink)} without
     *  progress. */
    public void write(GenomicRun run, List<ProvenanceRecord> provenance, boolean sort)
            throws IOException {
        write(run, provenance, sort, ProgressSink.discard());
    }

    // ------------------------------------------------------------------
    // Header
    // ------------------------------------------------------------------

    /**
     * Build a SAMFileHeader from {@code run} + provenance. Public-by-
     * package-default so {@link CramWriter} can read it back if needed
     * and so tests can introspect the generated header.
     */
    SAMFileHeader buildHeader(WrittenGenomicRun run,
                              List<ProvenanceRecord> provenance,
                              boolean sort) {
        return buildHeader(run.chromosomes(), run.sampleName(), run.platform(),
            provenance, sort);
    }

    /** {@code chromosomes} may list every read's chromosome or the
     *  run-level table; {@code @SQ} lines are emitted once per unique
     *  name in first-seen order, {@code *} excluded. */
    SAMFileHeader buildHeader(List<String> chromosomes, String sample, String platform,
                              List<ProvenanceRecord> provenance, boolean sort) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(sort
            ? SAMFileHeader.SortOrder.coordinate
            : SAMFileHeader.SortOrder.unsorted);

        // @SQ — one per unique chromosome in first-seen order
        // (excluding "*" SAM unmapped sentinel).
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        Set<String> seen = new LinkedHashSet<>();
        for (String chrom : chromosomes) {
            if (chrom == null || chrom.isEmpty() || "*".equals(chrom)) continue;
            if (!seen.add(chrom)) continue;
            dict.addSequence(new SAMSequenceRecord(chrom, (int) DEFAULT_SQ_LENGTH));
        }
        header.setSequenceDictionary(dict);

        // @RG — single line if either sample or platform is set.
        boolean hasSample = sample != null && !sample.isEmpty();
        boolean hasPlatform = platform != null && !platform.isEmpty();
        if (hasSample || hasPlatform) {
            SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
            if (hasSample) rg.setSample(sample);
            if (hasPlatform) rg.setPlatform(platform);
            header.addReadGroup(rg);
        }

        // @PG — one per provenance record. SAM requires unique ID;
        // synthesise "pg<idx>" if software is blank, suffix duplicates
        // with .1/.2/... matching the cross-language convention.
        Set<String> usedIds = new HashSet<>();
        for (int idx = 0; idx < provenance.size(); idx++) {
            ProvenanceRecord prov = provenance.get(idx);
            String software = prov.software();
            String baseId = (software == null || software.isEmpty())
                ? ("pg" + idx) : software;
            String pgId = baseId;
            int n = 1;
            while (usedIds.contains(pgId)) {
                pgId = baseId + "." + n;
                n++;
            }
            usedIds.add(pgId);
            SAMProgramRecord pg = new SAMProgramRecord(pgId);
            pg.setProgramName(software == null ? "" : software);
            String cl = prov.parameters() != null
                ? prov.parameters().get("CL") : null;
            if (cl != null && !cl.isEmpty()) {
                pg.setCommandLine(cl);
            }
            header.addProgramRecord(pg);
        }

        return header;
    }

    // ------------------------------------------------------------------
    // Records
    // ------------------------------------------------------------------

    /** One {@link AlignedRead} as a SAM record (same field rules as the
     *  array-based builder). */
    private SAMRecord buildSamRecord(AlignedRead r, SAMFileHeader header) {
        OrderedTagRecord rec = new OrderedTagRecord(header);
        rec.setTagText(r.tags());
        String qname = r.readName();
        rec.setReadName((qname == null || qname.isEmpty()) ? "*" : qname);
        rec.setFlags(r.flags());
        String rname = r.chromosome();
        rec.setReferenceName((rname == null || rname.isEmpty()) ? "*" : rname);
        rec.setAlignmentStart(r.position() > Integer.MAX_VALUE ? 0 : (int) r.position());
        rec.setMappingQuality(r.mappingQuality() & 0xFF);
        String cigar = r.cigar();
        rec.setCigarString((cigar == null || cigar.isEmpty()) ? "*" : cigar);
        String mate = r.mateChromosome();
        rec.setMateReferenceName((mate == null || mate.isEmpty()) ? "*" : mate);
        rec.setMateAlignmentStart(r.matePosition() < 0 ? 0 : (int) r.matePosition());
        rec.setInferredInsertSize(r.templateLength());
        String seq = r.sequence();
        if (seq == null || seq.isEmpty()) {
            rec.setReadBases(SAMRecord.NULL_SEQUENCE);
            rec.setBaseQualities(SAMRecord.NULL_QUALS);
        } else {
            rec.setReadBases(seq.getBytes(StandardCharsets.US_ASCII));
            byte[] q = r.qualities();
            boolean allFF = q != null && q.length > 0;
            if (q != null) for (byte b : q) if ((b & 0xFF) != 0xFF) { allFF = false; break; }
            if (q == null || q.length == 0 || allFF) {
                rec.setBaseQualities(SAMRecord.NULL_QUALS);
            } else {
                byte[] raw = new byte[q.length];
                for (int k = 0; k < q.length; k++) raw[k] = (byte) ((q[k] & 0xFF) - 33);
                rec.setBaseQualities(raw);
            }
        }
        return rec;
    }

    private SAMRecord buildSamRecord(WrittenGenomicRun run,
                                     SAMFileHeader header, int i) {
        OrderedTagRecord rec = new OrderedTagRecord(header);
        if (run.tags() != null) rec.setTagText(run.tags().get(i));

        String qnameRaw = run.readNames().get(i);
        rec.setReadName((qnameRaw == null || qnameRaw.isEmpty()) ? "*" : qnameRaw);

        rec.setFlags(run.flags()[i]);

        String rname = run.chromosomes().get(i);
        rec.setReferenceName((rname == null || rname.isEmpty()) ? "*" : rname);

        long pos = run.positions()[i];
        // htsjdk uses int alignment start; clamp out-of-range.
        rec.setAlignmentStart(pos > Integer.MAX_VALUE ? 0 : (int) pos);

        rec.setMappingQuality(run.mappingQualities()[i] & 0xFF);

        String cigar = run.cigars().get(i);
        rec.setCigarString((cigar == null || cigar.isEmpty()) ? "*" : cigar);

        // RNEXT collapse: only when mate matches and chromosome is not "*".
        String mateChromRaw = run.mateChromosomes().get(i);
        String mateChrom = (mateChromRaw == null || mateChromRaw.isEmpty())
            ? "*" : mateChromRaw;
        rec.setMateReferenceName(mateChrom);

        long matePos = run.matePositions()[i];
        rec.setMateAlignmentStart(matePos < 0 ? 0 : (int) matePos);

        rec.setInferredInsertSize(run.templateLengths()[i]);

        long offset = run.offsets()[i];
        int length = run.lengths()[i];
        if (length == 0) {
            rec.setReadBases(SAMRecord.NULL_SEQUENCE);
            rec.setBaseQualities(SAMRecord.NULL_QUALS);
        } else {
            int from = (int) offset;
            byte[] seqBytes = new byte[length];
            System.arraycopy(run.sequences(), from, seqBytes, 0, length);
            rec.setReadBases(seqBytes);

            // QUAL bytes are stored as ASCII Phred+33. htsjdk wants raw
            // Phred bytes; convert by subtracting 33. The all-0xFF
            // sentinel (M87 reader produces this for source QUAL '*'
            // with non-empty SEQ) maps back to NULL_QUALS so the round
            // trip canonicalises.
            byte[] qualBuf = run.qualities();
            boolean allFF = true;
            for (int k = from; k < from + length; k++) {
                if ((qualBuf[k] & 0xFF) != 0xFF) { allFF = false; break; }
            }
            if (allFF) {
                rec.setBaseQualities(SAMRecord.NULL_QUALS);
            } else {
                byte[] qualRaw = new byte[length];
                for (int k = 0; k < length; k++) {
                    int ascii = qualBuf[from + k] & 0xFF;
                    qualRaw[k] = (byte) (ascii - 33);
                }
                rec.setBaseQualities(qualRaw);
            }
        }

        return rec;
    }

    // ------------------------------------------------------------------
    // Writer factory (overridable by CramWriter)
    // ------------------------------------------------------------------

    /**
     * Build the SAMFileWriter for this output path + header. Default
     * implementation produces BAM; {@link CramWriter} overrides to
     * produce CRAM with a ReferenceSource.
     *
     * @param header   the header (sort order already set).
     * @param sort     whether output should be sorted on close.
     *                 (presorted=false + header SO=coordinate triggers
     *                 htsjdk's in-memory sort on close.)
     */
    protected SAMFileWriter makeWriter(SAMFileHeader header, boolean sort) {
        SAMFileWriterFactory factory = new SAMFileWriterFactory()
            .setCreateIndex(false)
            .setCreateMd5File(false);
        // presorted=true tells htsjdk "input is already sorted, don't
        // re-sort". When sort=false (output should be unsorted), we
        // also pass presorted=true so no sorting happens. When
        // sort=true, presorted=false so htsjdk sorts on close.
        boolean presorted = !sort;
        return factory.makeBAMWriter(header, presorted, path.toFile());
    }

    // ------------------------------------------------------------------
    // Backwards-compat: SAM text generation (test seam)
    // ------------------------------------------------------------------

    /**
     * Build the full SAM text (header + alignment lines) for
     * {@code run}, without involving htsjdk. Retained as a test seam:
     * the M88 cross-language conformance tests inspect the pre-write
     * SAM stream directly, and the legacy
     * {@code test_mate_collapse_to_equals} harness reads it through
     * this method.
     *
     * <p>Package-private visibility so {@link CramWriter} can reuse.</p>
     */
    String buildSamText(WrittenGenomicRun run,
                        List<ProvenanceRecord> provenance, boolean sort) {
        if (provenance == null) provenance = List.of();
        StringBuilder sb = new StringBuilder();
        appendHeader(sb, run, provenance, sort);
        appendAlignments(sb, run);
        return sb.toString();
    }

    private static void appendHeader(StringBuilder sb,
                                     WrittenGenomicRun run,
                                     List<ProvenanceRecord> provenance,
                                     boolean sort) {
        String so = sort ? "coordinate" : "unsorted";
        sb.append("@HD\tVN:1.6\tSO:").append(so).append('\n');

        Set<String> seen = new LinkedHashSet<>();
        for (String chrom : run.chromosomes()) {
            if (chrom == null || chrom.isEmpty() || "*".equals(chrom)) continue;
            if (!seen.add(chrom)) continue;
            sb.append("@SQ\tSN:").append(chrom)
              .append("\tLN:").append(DEFAULT_SQ_LENGTH).append('\n');
        }

        String sample = run.sampleName();
        String platform = run.platform();
        boolean hasSample = sample != null && !sample.isEmpty();
        boolean hasPlatform = platform != null && !platform.isEmpty();
        if (hasSample || hasPlatform) {
            sb.append("@RG\tID:rg1");
            if (hasSample) sb.append("\tSM:").append(sample);
            if (hasPlatform) sb.append("\tPL:").append(platform);
            sb.append('\n');
        }

        Set<String> usedIds = new HashSet<>();
        for (int idx = 0; idx < provenance.size(); idx++) {
            ProvenanceRecord prov = provenance.get(idx);
            String software = prov.software();
            String baseId = (software == null || software.isEmpty())
                ? ("pg" + idx) : software;
            String pgId = baseId;
            int n = 1;
            while (usedIds.contains(pgId)) {
                pgId = baseId + "." + n;
                n++;
            }
            usedIds.add(pgId);
            sb.append("@PG\tID:").append(pgId);
            sb.append("\tPN:").append(software == null ? "" : software);
            String cl = prov.parameters() != null
                ? prov.parameters().get("CL") : null;
            if (cl != null && !cl.isEmpty()) {
                sb.append("\tCL:").append(cl);
            }
            sb.append('\n');
        }
    }

    private static void appendAlignments(StringBuilder sb,
                                         WrittenGenomicRun run) {
        byte[] seqBuf = run.sequences();
        byte[] qualBuf = run.qualities();
        List<String> readNames = run.readNames();
        List<String> chromosomes = run.chromosomes();
        List<String> mateChromosomes = run.mateChromosomes();
        List<String> cigars = run.cigars();
        long[] positions = run.positions();
        long[] matePositions = run.matePositions();
        long[] offsets = run.offsets();
        int[] lengths = run.lengths();
        int[] flags = run.flags();
        int[] templateLengths = run.templateLengths();
        byte[] mappingQualities = run.mappingQualities();
        List<String> tags = run.tags();

        int n = readNames.size();
        for (int i = 0; i < n; i++) {
            String qnameRaw = readNames.get(i);
            String qname = (qnameRaw == null || qnameRaw.isEmpty()) ? "*" : qnameRaw;
            int flag = flags[i];
            String rnameRaw = chromosomes.get(i);
            String rname = (rnameRaw == null || rnameRaw.isEmpty()) ? "*" : rnameRaw;
            long pos = positions[i];
            int mapq = mappingQualities[i] & 0xFF;
            String cigarRaw = cigars.get(i);
            String cigar = (cigarRaw == null || cigarRaw.isEmpty()) ? "*" : cigarRaw;

            String mateChromRaw = mateChromosomes.get(i);
            String mateChrom = (mateChromRaw == null || mateChromRaw.isEmpty())
                ? "*" : mateChromRaw;
            String rnext;
            if (mateChrom.equals(rname) && !"*".equals(rname)) {
                rnext = "=";
            } else {
                rnext = mateChrom;
            }

            long matePos = matePositions[i];
            long pnext = matePos < 0 ? 0L : matePos;
            int tlen = templateLengths[i];
            long offset = offsets[i];
            int length = lengths[i];

            String seq;
            String qual;
            if (length == 0) {
                seq = "*";
                qual = "*";
            } else {
                int from = (int) offset;
                int to = from + length;
                seq = new String(seqBuf, from, length, StandardCharsets.US_ASCII);
                boolean allFF = true;
                for (int k = from; k < to; k++) {
                    if ((qualBuf[k] & 0xFF) != 0xFF) { allFF = false; break; }
                }
                if (allFF) {
                    qual = "*";
                } else {
                    qual = new String(qualBuf, from, length,
                        StandardCharsets.ISO_8859_1);
                }
            }

            sb.append(qname).append('\t')
              .append(flag).append('\t')
              .append(rname).append('\t')
              .append(pos).append('\t')
              .append(mapq).append('\t')
              .append(cigar).append('\t')
              .append(rnext).append('\t')
              .append(pnext).append('\t')
              .append(tlen).append('\t')
              .append(seq).append('\t')
              .append(qual);
            // M101: the optional fields follow column 11 verbatim.
            String tail = tags == null ? "" : tags.get(i);
            if (!tail.isEmpty()) sb.append('\t').append(tail);
            sb.append('\n');
        }
    }

    // ------------------------------------------------------------------
    // M101: SAM optional fields
    // ------------------------------------------------------------------

    /** A {@link SAMBinaryTagAndValue} that can be linked in a chosen
     *  order: htsjdk inserts attributes sorted by tag, which would
     *  reorder a read's optional fields on export. */
    static final class OrderedTag extends SAMBinaryTagAndValue {
        private static final long serialVersionUID = 1L;
        private final boolean unsignedArray;

        OrderedTag(String key, Object value, boolean unsignedArray) {
            super(SAMTag.makeBinaryTag(key), value);
            this.unsignedArray = unsignedArray;
        }

        @Override public boolean isUnsignedArray() { return unsignedArray; }

        void link(OrderedTag n) { this.next = n; }
    }

    /** A {@link SAMRecord} whose attributes keep the order of the
     *  stored tag text. */
    static final class OrderedTagRecord extends SAMRecord {
        private static final long serialVersionUID = 1L;

        OrderedTagRecord(SAMFileHeader header) { super(header); }

        /** Replace the attributes with the parsed tag text, in order. */
        void setTagText(String text) {
            List<OrderedTag> list = parseTagText(text);
            for (int k = 0; k + 1 < list.size(); k++) list.get(k).link(list.get(k + 1));
            setAttributes(list.isEmpty() ? null : list.get(0));
        }
    }

    /** Parse samtools-style tag text ({@code KEY:TYPE:VALUE}, tab-joined)
     *  into typed htsjdk values: {@code i} as Integer (Long past the
     *  int32 range; htsjdk picks the BAM width), {@code f} Float,
     *  {@code A} Character, {@code Z} String, {@code B} the matching
     *  primitive array (unsigned subtypes flagged). htsjdk has no
     *  {@code H} or {@code d} writer, so those are written as
     *  {@code Z}. */
    static List<OrderedTag> parseTagText(String text) {
        List<OrderedTag> out = new java.util.ArrayList<>();
        if (text == null || text.isEmpty()) return out;
        for (String f : text.split("\t", -1)) {
            if (f.length() < 5 || f.charAt(2) != ':' || f.charAt(4) != ':') {
                throw new IllegalArgumentException("malformed SAM optional field '" + f + "'");
            }
            String key = f.substring(0, 2);
            char type = f.charAt(3);
            String val = f.substring(5);
            Object v;
            boolean unsigned = false;
            switch (type) {
                case 'A' -> v = val.isEmpty() ? ' ' : val.charAt(0);
                case 'i' -> {
                    long l = Long.parseLong(val);
                    v = (l >= Integer.MIN_VALUE && l <= Integer.MAX_VALUE) ? (Object) (int) l : (Object) l;
                }
                case 'f' -> v = (float) Double.parseDouble(val);
                case 'B' -> {
                    String[] parts = val.split(",", -1);
                    char sub = parts[0].charAt(0);
                    int n = parts.length - 1;
                    unsigned = sub == 'C' || sub == 'S' || sub == 'I';
                    switch (sub) {
                        case 'c', 'C' -> {
                            byte[] a = new byte[n];
                            for (int k = 0; k < n; k++) a[k] = (byte) Long.parseLong(parts[k + 1]);
                            v = a;
                        }
                        case 's', 'S' -> {
                            short[] a = new short[n];
                            for (int k = 0; k < n; k++) a[k] = (short) Long.parseLong(parts[k + 1]);
                            v = a;
                        }
                        case 'i', 'I' -> {
                            int[] a = new int[n];
                            for (int k = 0; k < n; k++) a[k] = (int) Long.parseLong(parts[k + 1]);
                            v = a;
                        }
                        case 'f' -> {
                            float[] a = new float[n];
                            for (int k = 0; k < n; k++) a[k] = (float) Double.parseDouble(parts[k + 1]);
                            v = a;
                        }
                        default -> throw new IllegalArgumentException(
                            "SAM optional field '" + f + "': unknown B subtype '" + sub + "'");
                    }
                }
                default -> v = val;          // Z, and H / d written as Z
            }
            out.add(new OrderedTag(key, v, unsigned));
        }
        return out;
    }
}
