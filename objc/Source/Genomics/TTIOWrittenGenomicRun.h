#ifndef TTIO_WRITTEN_GENOMIC_RUN_H
#define TTIO_WRITTEN_GENOMIC_RUN_H

#import <Foundation/Foundation.h>
#import "ValueClasses/TTIOEnums.h"

@class TTIOProvenanceRecord;

NS_ASSUME_NONNULL_BEGIN

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Genomics/TTIOWrittenGenomicRun.h</p>
 *
 * <p>Write-side container for a single genomic run, passed to
 * <code>+[TTIOSpectralDataset writeMinimalToPath:...genomicRuns:...]</code>.
 * Genomic analogue of <code>TTIOWrittenRun</code>. Pure data — no
 * methods beyond accessors and the designated initialisers.</p>
 *
 * <p>The class holds parallel C-array buffers (positions, mapping
 * qualities, flags, sequences, qualities, offsets, lengths, mate
 * positions, template lengths) plus per-read variable-length string
 * arrays (cigars, read names, mate chromosomes, chromosomes). The
 * writer materialises these into the on-disk channel layout
 * described in <code>docs/format-spec.md</code> §10.4-§10.10.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python:
 * <code>ttio.written_genomic_run.WrittenGenomicRun</code><br/>
 * Java:
 * <code>global.thalion.ttio.genomics.WrittenGenomicRun</code></p>
 */
@interface TTIOWrittenGenomicRun : NSObject

/** Acquisition mode. */
@property (readonly) TTIOAcquisitionMode acquisitionMode;

/** Reference genome URI (e.g. <code>@"GRCh38.p14"</code>). */
@property (readonly, copy) NSString *referenceUri;

/** Sequencing platform identifier. */
@property (readonly, copy) NSString *platform;

/** Sample identifier. */
@property (readonly, copy) NSString *sampleName;

/** int64 mapping positions. */
@property (readonly, copy) NSData *positionsData;

/** uint8 mapping qualities. */
@property (readonly, copy) NSData *mappingQualitiesData;

/** uint32 SAM flags. */
@property (readonly, copy) NSData *flagsData;

/** Concatenated uint8 read sequences (one ASCII byte per base). */
@property (readonly, copy) NSData *sequencesData;

/** Concatenated uint8 quality scores (Phred). */
@property (readonly, copy) NSData *qualitiesData;

/** uint64 per-read offsets into <code>sequencesData</code> /
 *  <code>qualitiesData</code>. */
@property (readonly, copy) NSData *offsetsData;

/** uint32 per-read base counts. */
@property (readonly, copy) NSData *lengthsData;

/** Per-read CIGAR strings. */
@property (readonly, copy) NSArray<NSString *> *cigars;

/** Per-read read names. */
@property (readonly, copy) NSArray<NSString *> *readNames;

/** Per-read mate chromosome names. */
@property (readonly, copy) NSArray<NSString *> *mateChromosomes;

/** int64 per-read mate positions; <code>-1</code> sentinel for
 *  unmapped mates. */
@property (readonly, copy) NSData *matePositionsData;

/** int32 per-read template lengths; <code>0</code> sentinel for
 *  unpaired reads. */
@property (readonly, copy) NSData *templateLengthsData;

/** Per-read chromosome names (for the genomic index). */
@property (readonly, copy) NSArray<NSString *> *chromosomes;

/** HDF5-filter compression codec applied to non-genomic-codec
 *  channels. Defaults to <code>TTIOCompressionZlib</code>. */
@property (readonly) TTIOCompression signalCompression;

/** Removes the V5 sequence-context strategies from the qualities
 *  auto-tune set (spec 2.4). Python:
 *  <code>opt_disable_qualities_v5</code>; Java:
 *  <code>optDisableQualitiesV5</code>. Defaults to NO. */
@property (nonatomic) BOOL optDisableQualitiesV5;

/** Role of the reads in an assembly context, persisted as the
 *  <code>@read_role</code> UTF-8 run attribute when non-nil.
 *  Recognised values: <code>hifi</code>, <code>ont_ul</code>,
 *  <code>hic_r1</code>, <code>hic_r2</code>,
 *  <code>parental_maternal</code>, <code>parental_paternal</code>,
 *  <code>illumina_wgs</code>; other strings are stored unchecked.
 *  Python: <code>read_role</code>; Java: <code>readRole</code>.
 *  Defaults to nil (attribute absent). */
@property (nonatomic, copy, nullable) NSString *readRole;

/** M101: SAM optional fields per read (columns 12+, tab-joined as
 *  samtools prints them; <code>@""</code> for a read without tags).
 *  nil, or an array whose entries are all empty, writes no
 *  <code>tags</code> channel. Stored with SAM_TAGS (codec id 18);
 *  MD:Z / NM:i are recomputed from the reference when the sequences
 *  channel is REF_DIFF_V2. Python: <code>tags</code>; Java:
 *  <code>tags</code>. Defaults to nil. */
@property (nonatomic, copy, nullable) NSArray<NSString *> *tags;

/** YES when <code>tags</code> holds at least one non-empty entry,
 *  i.e. the writer will emit a <code>tags</code> channel. */
- (BOOL)hasTags;

/** REF_DIFF_V2 slice byte budget: a slice closes before the read
 *  that would push it past this many bases (the 10,000-read cap
 *  still applies). 0 (the default) keeps the fixed-count rule.
 *  Writer policy only — the wire format and decoder are unchanged.
 *  Python: <code>ref_diff_slice_bytes</code>; Java:
 *  <code>refDiffSliceBytes</code>. */
@property (nonatomic) unsigned long long refDiffSliceBytes;

/**
 * Per-channel codec opt-in. Maps channel name
 * (<code>NSString *</code>) to a boxed
 * <code>TTIOCompression</code> value (<code>NSNumber *</code>).
 * Channels not listed use the
 * <code>signalCompression</code> path. Cross-language equivalent
 * of Python's
 * <code>WrittenGenomicRun.signal_codec_overrides</code>.
 */
@property (readonly, copy) NSDictionary<NSString *, NSNumber *> *signalCodecOverrides;

/** Per-run provenance records. Persisted under
 *  <code>&lt;run&gt;/provenance/steps</code> by the writer; read
 *  back via <code>-[TTIOGenomicRun provenanceChain]</code>.
 *  Defaults to an empty array. */
@property (nonatomic, copy) NSArray<TTIOProvenanceRecord *> *provenanceRecords;

/** When <code>YES</code> (default) and a context-aware codec is
 *  selected on the <code>sequences</code> channel, the writer
 *  embeds the chromosome sequences supplied in
 *  <code>referenceChromSeqs</code> at
 *  <code>/study/references/&lt;referenceUri&gt;/</code>. */
@property (nonatomic, assign) BOOL embedReference;

/** Map from chromosome name to uppercase ACGTN bytes. Required
 *  when <code>REF_DIFF_V2</code> is selected on
 *  <code>sequences</code> and <code>embedReference</code> is
 *  <code>YES</code>; otherwise the writer falls back silently to
 *  <code>BASE_PACK</code> on this channel. */
@property (nonatomic, copy, nullable) NSDictionary<NSString *, NSData *> *referenceChromSeqs;

/** External reference path stamped into file metadata for decoder
 *  fallback when the embedded reference is absent. The writer
 *  never reads this path; metadata only. */
@property (nonatomic, copy, nullable) NSString *externalReferencePath;

/** Number of reads in the run. */
@property (readonly) NSUInteger readCount;

/** Write this run in the pre-1.9 whole-channel layout instead of
 *  blocks_v1 (format-spec 10.12). Python:
 *  <code>opt_legacy_whole_channel</code>; Java:
 *  <code>optLegacyWholeChannel</code>. Defaults to NO. */
@property (nonatomic) BOOL optLegacyWholeChannel;

/** A copy of this run with a different per-channel codec map. Every
 *  other field, including the mutable options, is carried over. */
- (instancetype)copyWithSignalCodecOverrides:(NSDictionary<NSString *, NSNumber *> *)overrides;

/** A copy of this run with a different provenance chain. */
- (instancetype)copyWithProvenance:(NSArray<TTIOProvenanceRecord *> *)records;

/** A copy of this run with <code>optLegacyWholeChannel</code> set. */
- (instancetype)copyWithOptLegacyWholeChannel:(BOOL)legacy;

/// Phase 2c-T: verbatim v2 codec blobs for direct on-disk write,
/// bypassing the v2 codec encode step in the writer. Used by the
/// transport bulk-mode receiver. nil disables.
@property (nonatomic, strong, nullable) id /* TTIOBulkV2Blobs */ bulkV2Blobs;

/**
 * Convenience initialiser without per-channel codec overrides;
 * delegates to the designated initialiser with an empty overrides
 * dictionary.
 */
- (instancetype)initWithAcquisitionMode:(TTIOAcquisitionMode)mode
                           referenceUri:(NSString *)referenceUri
                               platform:(NSString *)platform
                             sampleName:(NSString *)sampleName
                              positions:(NSData *)positions
                       mappingQualities:(NSData *)mappingQualities
                                  flags:(NSData *)flags
                              sequences:(NSData *)sequences
                              qualities:(NSData *)qualities
                                offsets:(NSData *)offsets
                                lengths:(NSData *)lengths
                                 cigars:(NSArray<NSString *> *)cigars
                              readNames:(NSArray<NSString *> *)readNames
                        mateChromosomes:(NSArray<NSString *> *)mateChromosomes
                          matePositions:(NSData *)matePositions
                        templateLengths:(NSData *)templateLengths
                            chromosomes:(NSArray<NSString *> *)chromosomes
                      signalCompression:(TTIOCompression)signalCompression;

/**
 * Designated initialiser including the per-channel codec overrides
 * dictionary.
 *
 * Populates every channel of a writeable genomic run from parallel
 * arrays. `positions`, `mappingQualities`, `flags`, `offsets`,
 * `lengths`, `matePositions`, and `templateLengths` are typed
 * `NSData` buffers (int64, uint8, uint32, uint64, uint32, int64,
 * int32 respectively); `sequences` and `qualities` are concatenated
 * uint8 buffers sliced by `offsets` / `lengths`.
 *
 * @param mode                  Acquisition mode (`Genomic`, etc).
 * @param referenceUri          Reference URI string.
 * @param platform              Sequencing platform identifier.
 * @param sampleName            Sample name.
 * @param positions             Per-read int64 reference positions.
 * @param mappingQualities      Per-read uint8 MAPQ values.
 * @param flags                 Per-read uint32 SAM flags.
 * @param sequences             Concatenated ACGTN sequence bytes.
 * @param qualities             Concatenated Phred+33 quality bytes.
 * @param offsets               Per-read uint64 offsets into sequences/qualities.
 * @param lengths               Per-read uint32 read lengths.
 * @param cigars                Per-read CIGAR strings.
 * @param readNames             Per-read query names.
 * @param mateChromosomes       Per-read mate chromosome names.
 * @param matePositions         Per-read int64 mate positions.
 * @param templateLengths       Per-read int32 template lengths.
 * @param chromosomes           Per-read reference chromosome names.
 * @param signalCompression     Default compression for signal channels.
 * @param signalCodecOverrides  Per-channel-name codec id overrides
 *                              (empty dict = use defaults).
 * @return Initialised `TTIOWrittenGenomicRun`.
 */
- (instancetype)initWithAcquisitionMode:(TTIOAcquisitionMode)mode
                           referenceUri:(NSString *)referenceUri
                               platform:(NSString *)platform
                             sampleName:(NSString *)sampleName
                              positions:(NSData *)positions
                       mappingQualities:(NSData *)mappingQualities
                                  flags:(NSData *)flags
                              sequences:(NSData *)sequences
                              qualities:(NSData *)qualities
                                offsets:(NSData *)offsets
                                lengths:(NSData *)lengths
                                 cigars:(NSArray<NSString *> *)cigars
                              readNames:(NSArray<NSString *> *)readNames
                        mateChromosomes:(NSArray<NSString *> *)mateChromosomes
                          matePositions:(NSData *)matePositions
                        templateLengths:(NSData *)templateLengths
                            chromosomes:(NSArray<NSString *> *)chromosomes
                      signalCompression:(TTIOCompression)signalCompression
                   signalCodecOverrides:(NSDictionary<NSString *, NSNumber *> *)signalCodecOverrides;

@end

NS_ASSUME_NONNULL_END

#endif
