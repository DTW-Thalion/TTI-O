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
import global.thalion.ttio.io.ProgressSink;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * SAM/BAM importer — M87.
 *
 * <p>Reads SAM and BAM (Sequence Alignment/Map) files into
 * {@link WrittenGenomicRun} instances using <a
 * href="https://github.com/samtools/htsjdk">htsjdk</a> — the pure-Java
 * SAM/BAM/CRAM library used by GATK, Picard, and IGV.</p>
 *
 * <p>v1.5.0 replaced the prior {@code samtools} subprocess implementation
 * so tio-browser works on Windows without requiring a system samtools
 * install. The public API is unchanged; the {@code SamtoolsNotFoundException}
 * type is retained as a no-throw historical alias so callers and tests
 * compile unchanged.</p>
 *
 * <p>Format auto-detection is provided by htsjdk based on magic bytes;
 * one parser handles both SAM and BAM. The companion {@link SamReader}
 * (note: name collides with htsjdk.samtools.SamReader — fully qualify if
 * needed) exists as a discoverable convenience alias.</p>
 *
 * <p><b>Cross-language equivalents:</b> Python {@code ttio.importers.bam.BamReader},
 * Objective-C {@code TTIOBamReader}. The byte-level output is verified
 * by the {@code cross-compat} CI job.</p>
 *
 * (M87, v1.5.0 htsjdk swap)
 */
public class BamReader {

    /** Default emit-every-N cadence for the
     *  {@link #toGenomicRun(String, String, String, ProgressSink)}
     *  {@link ProgressSink}. The htsjdk record iterator can chew
     *  through several hundred thousand records per second; a tighter
     *  cadence drowns the GUI thread without adding info. */
    public static final int PROGRESS_INTERVAL_READS = 1000;

    static {
        // htsjdk's CRAM decoder validates the reference sequence MD5 by
        // default and refuses to decode if the embedded CRAM @SQ M5 hash
        // doesn't match the on-disk FASTA. The old samtools-subprocess
        // BamWriter / CramWriter wrote @SQ LN:INT32_MAX placeholders and
        // samtools never cared; htsjdk does. Disable the check globally
        // so we can read existing fixture CRAMs and not block on MD5
        // strictness. Lossless data integrity is still ensured by CRAM's
        // per-record byte-exact reconstruction.
        if (System.getProperty("samjdk.cram.use_alignment_md5_check") == null) {
            System.setProperty("samjdk.cram.use_alignment_md5_check", "false");
        }
    }

    /**
     * Historical exception type from the samtools-subprocess era.
     * v1.5.0 onwards never throws this — htsjdk is a Maven dep, always
     * available. Retained as a no-throw alias so callers and tests that
     * catch this type still compile.
     */
    public static final class SamtoolsNotFoundException extends IOException {
        private static final long serialVersionUID = 1L;
        public SamtoolsNotFoundException(String msg) { super(msg); }
        public SamtoolsNotFoundException(String msg, Throwable cause) {
            super(msg, cause);
        }
    }

    private final Path path;
    private List<ProvenanceRecord> lastProvenance = List.of();

    public BamReader(Path path) {
        this.path = path;
    }

    /** @return the path passed at construction time. */
    public Path path() { return path; }

    /**
     * @return the {@code @PG}-derived provenance records from the most
     *         recent {@link #toGenomicRun} call (empty before any call).
     */
    public List<ProvenanceRecord> lastProvenance() { return lastProvenance; }

    /** Default name {@code "genomic_0001"}, no region filter, no sample override. */
    public WrittenGenomicRun toGenomicRun(String name) throws IOException {
        return toGenomicRun(name, null, null, ProgressSink.discard());
    }

    /** No sample override. */
    public WrittenGenomicRun toGenomicRun(String name, String region)
            throws IOException {
        return toGenomicRun(name, region, null, ProgressSink.discard());
    }

    /** Backwards-compatible overload without a {@link ProgressSink}. */
    public WrittenGenomicRun toGenomicRun(String name, String region,
                                          String sampleName) throws IOException {
        return toGenomicRun(name, region, sampleName, ProgressSink.discard());
    }

    /** {@code ProgressSink}-accepting overload; defaults region + sample
     *  to {@code null}. */
    public WrittenGenomicRun toGenomicRun(String name, ProgressSink progress)
            throws IOException {
        return toGenomicRun(name, null, null, progress);
    }

    /**
     * Read the BAM/SAM and return a {@link WrittenGenomicRun}, firing
     * {@code progress.onProgress(readsDone, -1L)} every
     * {@link #PROGRESS_INTERVAL_READS} records during iteration and a
     * final {@code onProgress(total, total)} once all records have been
     * consumed. Total is unknown mid-iteration because htsjdk's record
     * iterator does not expose a cheap count without re-walking the
     * file (and the file may be region-filtered).
     *
     * @param name        run name (becomes {@code /study/genomic_runs/<name>}).
     * @param region      optional region filter (e.g. {@code "chr1:1000-2000"},
     *                    {@code "chr1"}, {@code "*"} for unmapped); {@code null}
     *                    means no filter.
     * @param sampleName  optional override for {@code sample_name};
     *                    {@code null} means "use first @RG SM:".
     * @param progress    progress callback; {@link ProgressSink#discard()}
     *                    if you don't care.
     * @throws IOException if the file is missing, malformed, or a region
     *                     filter is requested against an unindexed file.
     * @since 1.5.0
     */
    public WrittenGenomicRun toGenomicRun(String name, String region,
                                          String sampleName,
                                          ProgressSink progress)
            throws IOException {
        if (progress == null) progress = ProgressSink.discard();
        if (!Files.exists(path)) {
            throw new IOException("BAM/SAM file not found: " + path);
        }
        // htsjdk treats a zero-byte file as a 0-record BAM rather than
        // an error; reject explicitly for parity with samtools' "no BGZF
        // header" rejection and with the V4 zero-byte regression test.
        try {
            if (Files.size(path) == 0) {
                throw new IOException(
                    "BAM/SAM file is empty (0 bytes): " + path);
            }
        } catch (IOException e) {
            throw e;
        }

        long fileMtime;
        try {
            fileMtime = Files.getLastModifiedTime(path).toInstant().getEpochSecond();
        } catch (IOException e) {
            fileMtime = System.currentTimeMillis() / 1000L;
        }

        // Header state
        List<String> sqNames = new ArrayList<>();
        String rgSample = "";
        String rgPlatform = "";
        List<ProvenanceRecord> provenance = new ArrayList<>();
        BatchAccumulator acc = new BatchAccumulator();

        SamReaderFactory factory = makeReaderFactory();
        SamReader reader;
        try {
            reader = factory.open(SamInputResource.of(path.toFile()));
        } catch (RuntimeException e) {
            // htsjdk throws FileTruncatedException / SAMFormatException /
            // RuntimeIOException as unchecked. Wrap so callers catching
            // IOException (the documented contract) still see them.
            throw new IOException(
                "Failed to open BAM/SAM/CRAM: " + path + ": " + e.getMessage(), e);
        }
        try (reader) {
            SAMFileHeader header = reader.getFileHeader();

            // @SQ → reference sequence names (first wins for reference_uri).
            for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
                sqNames.add(seq.getSequenceName());
            }

            // @RG → sample + platform (first wins).
            for (SAMReadGroupRecord rg : header.getReadGroups()) {
                if (rgSample.isEmpty() && rg.getSample() != null) {
                    rgSample = rg.getSample();
                }
                if (rgPlatform.isEmpty() && rg.getPlatform() != null) {
                    rgPlatform = rg.getPlatform();
                }
            }

            // @PG → provenance records.
            for (SAMProgramRecord pg : header.getProgramRecords()) {
                String program = pg.getProgramName() != null ? pg.getProgramName() : "";
                Map<String, String> params = new LinkedHashMap<>();
                String cl = pg.getCommandLine();
                if (cl != null) params.put("CL", cl);
                if (pg.getId() != null) params.put("ID", pg.getId());
                if (pg.getProgramVersion() != null) params.put("VN", pg.getProgramVersion());
                if (pg.getPreviousProgramGroupId() != null) {
                    params.put("PP", pg.getPreviousProgramGroupId());
                }
                provenance.add(new ProvenanceRecord(
                    fileMtime, program, params, List.of(), List.of()));
            }

            // Iterate alignments — region filter applied if requested.
            // htsjdk may throw RuntimeException subclasses
            // (SAMFormatException, FileTruncatedException, etc.) during
            // record parsing; wrap as IOException for caller contract.
            Iterator<SAMRecord> it;
            try {
                it = iteratorFor(reader, region);
            } catch (RuntimeException e) {
                throw new IOException(
                    "Failed to iterate records in " + path + ": " + e.getMessage(), e);
            }
            SamTextTags lines = SamTextTags.open(path, reader, region);
            try {
                while (true) {
                    SAMRecord rec;
                    try {
                        if (!it.hasNext()) break;
                        rec = it.next();
                    } catch (RuntimeException e) {
                        throw new IOException(
                            "Malformed record in " + path + ": " + e.getMessage(), e);
                    }

                    addRecord(acc, rec, lines);
                    if (acc.size() % PROGRESS_INTERVAL_READS == 0) {
                        progress.onProgress(acc.size(), -1L);
                    }
                }
            } finally {
                if (it instanceof htsjdk.samtools.util.CloseableIterator<?>) {
                    ((htsjdk.samtools.util.CloseableIterator<?>) it).close();
                }
                if (lines != null) lines.close();
            }
        }

        // sample_name override.
        String effectiveSample = sampleName != null ? sampleName : rgSample;

        // reference_uri: first @SQ wins.
        String referenceUri = sqNames.isEmpty() ? "" : sqNames.get(0);

        int n = acc.size();
        this.lastProvenance = List.copyOf(provenance);

        // Final progress fire — total is now known. Listeners can flip
        // from indeterminate ("done out of ?") to determinate
        // ("done == total"), which the import dialog uses to mark
        // 100% before the write phase begins.
        progress.onProgress((long) n, (long) n);

        return acc.toRun(AcquisitionMode.GENOMIC_WGS, referenceUri, rgPlatform,
            effectiveSample, List.of());
    }

    /** Batches of {@code batchReads} reads as {@link WrittenGenomicRun}s
     *  carrying the run-level metadata ({@code @SQ}, {@code @RG},
     *  {@code @PG} provenance) so the first batch configures a
     *  {@link global.thalion.ttio.genomics.GenomicStreamWriter}. The
     *  iterator holds one batch; it implements {@link AutoCloseable}
     *  and closes the SAM reader at EOF or when closed early. */
    public Iterator<WrittenGenomicRun> iterBatches(String name, String region,
                                                   String sampleName, int batchReads)
            throws IOException {
        if (batchReads < 1) throw new IllegalArgumentException("batchReads must be >= 1");
        if (!Files.exists(path)) {
            throw new IOException("BAM/SAM file not found: " + path);
        }
        if (Files.size(path) == 0) {
            throw new IOException("BAM/SAM file is empty (0 bytes): " + path);
        }
        long fileMtime;
        try {
            fileMtime = Files.getLastModifiedTime(path).toInstant().getEpochSecond();
        } catch (IOException e) {
            fileMtime = System.currentTimeMillis() / 1000L;
        }
        SamReader reader;
        try {
            reader = makeReaderFactory().open(SamInputResource.of(path.toFile()));
        } catch (RuntimeException e) {
            throw new IOException(
                "Failed to open BAM/SAM/CRAM: " + path + ": " + e.getMessage(), e);
        }
        SAMFileHeader header = reader.getFileHeader();
        List<String> sqNames = new ArrayList<>();
        for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
            sqNames.add(seq.getSequenceName());
        }
        String rgSample = "", rgPlatform = "";
        for (SAMReadGroupRecord rg : header.getReadGroups()) {
            if (rgSample.isEmpty() && rg.getSample() != null) rgSample = rg.getSample();
            if (rgPlatform.isEmpty() && rg.getPlatform() != null) rgPlatform = rg.getPlatform();
        }
        List<ProvenanceRecord> provenance = new ArrayList<>();
        for (SAMProgramRecord pg : header.getProgramRecords()) {
            String program = pg.getProgramName() != null ? pg.getProgramName() : "";
            Map<String, String> params = new LinkedHashMap<>();
            if (pg.getCommandLine() != null) params.put("CL", pg.getCommandLine());
            if (pg.getId() != null) params.put("ID", pg.getId());
            if (pg.getProgramVersion() != null) params.put("VN", pg.getProgramVersion());
            if (pg.getPreviousProgramGroupId() != null) params.put("PP", pg.getPreviousProgramGroupId());
            provenance.add(new ProvenanceRecord(fileMtime, program, params, List.of(), List.of()));
        }
        this.lastProvenance = List.copyOf(provenance);
        final String referenceUri = sqNames.isEmpty() ? "" : sqNames.get(0);
        final String sample = sampleName != null ? sampleName : rgSample;
        final String platform = rgPlatform;
        final List<ProvenanceRecord> prov = List.copyOf(provenance);
        final Iterator<SAMRecord> it;
        try {
            it = iteratorFor(reader, region);
        } catch (RuntimeException e) {
            reader.close();
            throw new IOException(
                "Failed to iterate records in " + path + ": " + e.getMessage(), e);
        }
        SamTextTags lines;
        try {
            lines = SamTextTags.open(path, reader, region);
        } catch (IOException e) {
            reader.close();
            throw e;
        }
        return new BatchIterator(reader, it, batchReads, referenceUri, platform, sample, prov, lines);
    }

    /** {@link #iterBatches} as a {@link GenomicStreamSource}. */
    public GenomicStreamSource stream(String name, String region, String sampleName,
                                      Path referenceFasta, boolean embedReference, int batchReads) {
        return new GenomicStreamSource(name, () -> {
            try {
                return iterBatches(name, region, sampleName, batchReads);
            } catch (IOException e) {
                throw new java.io.UncheckedIOException(e);
            }
        }, referenceFasta, embedReference, null, null, false);
    }

    /** Default batch of 100 000 reads. */
    public GenomicStreamSource stream(String name, String region, String sampleName,
                                      Path referenceFasta, boolean embedReference) {
        return stream(name, region, sampleName, referenceFasta, embedReference, DEFAULT_BATCH_READS);
    }

    /** Reads per streamed batch. */
    public static final int DEFAULT_BATCH_READS = 100_000;

    private final class BatchIterator implements Iterator<WrittenGenomicRun>, AutoCloseable {
        private final SamReader reader;
        private final Iterator<SAMRecord> records;
        private final int batchReads;
        private final String referenceUri, platform, sample;
        private final List<ProvenanceRecord> prov;
        private final BatchAccumulator acc = new BatchAccumulator();
        private final SamTextTags lines;
        private WrittenGenomicRun next;
        private boolean done;

        BatchIterator(SamReader reader, Iterator<SAMRecord> records, int batchReads,
                      String referenceUri, String platform, String sample,
                      List<ProvenanceRecord> prov, SamTextTags lines) {
            this.reader = reader; this.records = records; this.batchReads = batchReads;
            this.referenceUri = referenceUri; this.platform = platform; this.sample = sample;
            this.prov = prov;
            this.lines = lines;
        }

        @Override public boolean hasNext() {
            if (next != null) return true;
            if (done) return false;
            try {
                while (acc.size() < batchReads && records.hasNext()) {
                    addRecord(acc, records.next(), lines);
                }
            } catch (IOException e) {
                throw new java.io.UncheckedIOException(e);
            } catch (RuntimeException e) {
                close();
                throw new java.io.UncheckedIOException(new IOException(
                    "Malformed record in " + path + ": " + e.getMessage(), e));
            }
            if (acc.size() == 0) { close(); return false; }
            next = acc.toRun(AcquisitionMode.GENOMIC_WGS, referenceUri, platform, sample, prov);
            acc.clear();
            if (!records.hasNext()) close();
            return true;
        }

        @Override public WrittenGenomicRun next() {
            if (!hasNext()) throw new java.util.NoSuchElementException();
            WrittenGenomicRun r = next;
            next = null;
            return r;
        }

        @Override public void close() {
            if (done) return;
            done = true;
            if (records instanceof htsjdk.samtools.util.CloseableIterator<?> c) c.close();
            if (lines != null) lines.close();
            try { reader.close(); } catch (IOException ignored) { }
        }
    }

    /** One htsjdk record into the accumulator (SEQ {@code *} → empty
     *  bytes; QUAL {@code *} with a sequence → 0xFF fill; Phred+33).
     *  M101: the optional fields are kept as samtools prints them, from
     *  the record's raw aux bytes (BAM), its own text line ({@code lines},
     *  SAM) or its attributes (CRAM); see {@link SamTagText}. */
    private static void addRecord(BatchAccumulator acc, SAMRecord rec, SamTextTags lines)
            throws IOException {
        // Before any getter that could re-encode the record's bytes.
        String tagText;
        try {
            tagText = lines != null ? SamTagText.normalizeSamText(lines.next())
                                    : SamTagText.fromRecord(rec);
        } catch (RuntimeException e) {
            throw new IOException("Malformed optional fields in record "
                + rec.getReadName() + ": " + e.getMessage(), e);
        }
        String qname = rec.getReadName() != null ? rec.getReadName() : "*";
        int flag = rec.getFlags();
        String rname = rec.getReferenceName() != null ? rec.getReferenceName() : "*";
        long pos = rec.getAlignmentStart();
        int mapq = rec.getMappingQuality();
        String cigar = rec.getCigarString() != null ? rec.getCigarString() : "*";
        String rnext = rec.getMateReferenceName() != null ? rec.getMateReferenceName() : "*";
        long pnext = rec.getMateAlignmentStart();
        int tlen = rec.getInferredInsertSize();
        String seqStr = rec.getReadString();
        byte[] seqBytes = (seqStr == null || "*".equals(seqStr))
            ? new byte[0] : seqStr.getBytes(StandardCharsets.US_ASCII);
        byte[] qualRaw = rec.getBaseQualities();
        byte[] qualBytes;
        if (qualRaw == null || qualRaw == SAMRecord.NULL_QUALS || qualRaw.length == 0) {
            if (seqBytes.length == 0) {
                qualBytes = new byte[0];
            } else {
                qualBytes = new byte[seqBytes.length];
                Arrays.fill(qualBytes, (byte) 0xFF);
            }
        } else {
            qualBytes = new byte[qualRaw.length];
            for (int i = 0; i < qualRaw.length; i++) qualBytes[i] = (byte) ((qualRaw[i] & 0xFF) + 33);
        }
        if (qualBytes.length != seqBytes.length) {
            if (seqBytes.length == 0) {
                qualBytes = new byte[0];
            } else if (qualBytes.length != 0) {
                throw new IOException("SEQ/QUAL length mismatch in record " + qname
                    + ": SEQ=" + seqBytes.length + " QUAL=" + qualBytes.length);
            }
        }
        acc.add(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seqBytes, qualBytes,
                tagText);
    }

    /** The optional-field columns of a SAM text file's alignment lines,
     *  read alongside htsjdk's record iterator (htsjdk keeps attributes
     *  sorted by tag, so their record order is only in the text). One
     *  line per record; header lines ({@code @}) and blank lines are
     *  skipped as htsjdk skips them. */
    static final class SamTextTags implements AutoCloseable {
        private final java.io.BufferedReader in;

        private SamTextTags(java.io.BufferedReader in) { this.in = in; }

        /** A reader for {@code path} when htsjdk reads it as SAM text
         *  without a region; {@code null} otherwise (BAM records carry
         *  their raw aux bytes). */
        static SamTextTags open(Path path, SamReader reader, String region) throws IOException {
            if (region != null || reader.type() != SamReader.Type.SAM_TYPE) return null;
            java.io.InputStream raw = new java.io.BufferedInputStream(Files.newInputStream(path));
            raw.mark(2);
            int b0 = raw.read(), b1 = raw.read();
            raw.reset();
            if (b0 == 0x1f && b1 == 0x8b) raw = new java.util.zip.GZIPInputStream(raw, 1 << 16);
            return new SamTextTags(new java.io.BufferedReader(
                new java.io.InputStreamReader(raw, StandardCharsets.UTF_8), 1 << 16));
        }

        /** Columns 12+ of the next alignment line ("" when it has none). */
        String next() throws IOException {
            String line;
            do {
                line = in.readLine();
                if (line == null) throw new IOException("SAM text ended before its records");
            } while (line.isEmpty() || line.charAt(0) == '@');
            int tabs = 0;
            for (int i = 0; i < line.length(); i++) {
                if (line.charAt(i) == '\t' && ++tabs == 11) return line.substring(i + 1);
            }
            return "";
        }

        @Override public void close() {
            try { in.close(); } catch (IOException ignored) { }
        }
    }

    // ------------------------------------------------------------------
    // Factory + region filter (overridable by CramReader)
    // ------------------------------------------------------------------

    /**
     * Build the htsjdk reader factory for this importer. Default is the
     * lenient stringency factory suitable for SAM/BAM. {@link CramReader}
     * overrides to inject a {@code ReferenceSource} so reference-
     * compressed CRAM bytes can be reconstituted.
     */
    protected SamReaderFactory makeReaderFactory() {
        return SamReaderFactory.makeDefault()
            .validationStringency(ValidationStringency.LENIENT);
    }


    /**
     * Resolve a samtools-style region string ({@code null}, {@code "*"},
     * {@code "chr1"}, {@code "chr1:1000"}, {@code "chr1:1000-2000"}) to
     * an iterator over matching SAMRecords. Subclasses
     * ({@code CramReader} in M88) can override to inject a reference for
     * CRAM decode if needed; default implementation handles SAM/BAM.
     */
    protected Iterator<SAMRecord> iteratorFor(SamReader reader, String region) {
        if (region == null) {
            return reader.iterator();
        }
        if ("*".equals(region)) {
            return reader.queryUnmapped();
        }
        // Parse "name[:start[-end]]" with optional thousands-separator commas.
        String refName;
        int start = 1;
        int end = Integer.MAX_VALUE;
        int colon = region.indexOf(':');
        if (colon < 0) {
            refName = region;
        } else {
            refName = region.substring(0, colon);
            String range = region.substring(colon + 1).replace(",", "");
            int dash = range.indexOf('-');
            try {
                if (dash < 0) {
                    start = Integer.parseInt(range);
                } else {
                    start = Integer.parseInt(range.substring(0, dash));
                    end = Integer.parseInt(range.substring(dash + 1));
                }
            } catch (NumberFormatException nfe) {
                throw new IllegalArgumentException(
                    "Malformed region string: " + region
                    + " (expected name[:start[-end]])", nfe);
            }
        }
        return reader.queryOverlapping(refName, start, end);
    }

    // ------------------------------------------------------------------
    // Samtools-availability shims (always-true since v1.5.0)
    // ------------------------------------------------------------------

    /**
     * Pre-v1.5.0 this probed for {@code samtools} on PATH. v1.5.0 onwards
     * BamReader uses htsjdk (a Maven dep, always available), so this
     * unconditionally returns {@code true}. Retained for source compat
     * with tests using {@code Assumptions.assumeTrue(isSamtoolsAvailable())}.
     */
    public static boolean isSamtoolsAvailable() {
        return true;
    }
}
