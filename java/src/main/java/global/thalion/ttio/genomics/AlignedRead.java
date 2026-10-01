/*
 * TTI-O Java Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
package global.thalion.ttio.genomics;

/**
 * One aligned sequencing read — the genomic analogue of
 * {@link global.thalion.ttio.MassSpectrum MassSpectrum}.
 *
 * <p>Immutable value record materialised by {@link GenomicRun#readAt(int)}
 * from the signal channels under
 * {@code /study/genomic_runs/<name>/signal_channels/}. No HDF5 I/O on
 * this class directly.</p>
 *
 * <p>Flag bits follow SAM convention:
 * 0x1 = paired, 0x4 = unmapped, 0x10 = reverse-strand,
 * 0x100 = secondary alignment, 0x800 = supplementary alignment.
 * Java's signed {@code byte} is sign-extended on read; callers should
 * mask with {@code & 0xFF} when interpreting raw base ASCII or quality
 * scores. The record's {@code qualities} field is {@code byte[]} for
 * compatibility with HDF5 UINT8 reads.</p>
 *
 * <p><b>API status:</b> Stable ().</p>
 *
 * <p><b>Cross-language equivalents:</b> Objective-C
 * {@code TTIOAlignedRead}, Python {@code ttio.aligned_read.AlignedRead}.</p>
 *
 * @param readName       Read identifier (e.g. {@code "read_000042"}).
 * @param chromosome     Reference sequence name (e.g. {@code "chr1"}).
 * @param position       0-based mapping position.
 * @param mappingQuality Phred-scaled mapping quality (0–93 typical, max 255).
 * @param cigar          CIGAR string (e.g. {@code "150M"}).
 * @param sequence       Base sequence in ACGTN ASCII.
 * @param qualities      Phred quality scores as raw bytes.
 * @param flags          SAM flags (UINT32 on disk; ints fit signed range).
 * @param mateChromosome {@code ""} if unpaired.
 * @param matePosition   {@code -1} if unpaired.
 * @param templateLength {@code 0} if unpaired.
 * @param tags           M101: the SAM optional fields (columns 12+),
 *                       tab-joined as samtools prints them; {@code ""}
 *                       when the read has none.
 */
public record AlignedRead(
    String readName,
    String chromosome,
    long   position,
    int    mappingQuality,
    String cigar,
    String sequence,
    byte[] qualities,
    int    flags,
    String mateChromosome,
    long   matePosition,
    int    templateLength,
    String tags
) {
    public AlignedRead {
        if (tags == null) tags = "";
    }

    /** Pre-M101 signature (11 components); no SAM tags. */
    public AlignedRead(String readName, String chromosome, long position,
                       int mappingQuality, String cigar, String sequence,
                       byte[] qualities, int flags, String mateChromosome,
                       long matePosition, int templateLength) {
        this(readName, chromosome, position, mappingQuality, cigar, sequence,
             qualities, flags, mateChromosome, matePosition, templateLength, "");
    }

    /** SAM flag 0x4: read failed to map. */
    public boolean isMapped()        { return (flags & 0x4) == 0; }
    /** SAM flag 0x1: read is part of a pair. */
    public boolean isPaired()        { return (flags & 0x1) != 0; }
    /** SAM flag 0x10: read is reverse-complemented relative to reference. */
    public boolean isReverse()       { return (flags & 0x10) != 0; }
    /** SAM flag 0x100: secondary alignment (multi-mapping). */
    public boolean isSecondary()     { return (flags & 0x100) != 0; }
    /** SAM flag 0x800: supplementary alignment (chimeric). */
    public boolean isSupplementary() { return (flags & 0x800) != 0; }
    /** Read length in bases (length of {@link #sequence}). */
    public int readLength()          { return sequence.length(); }
}
