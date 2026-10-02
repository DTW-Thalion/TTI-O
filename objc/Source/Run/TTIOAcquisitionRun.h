#ifndef TTIO_ACQUISITION_RUN_H
#define TTIO_ACQUISITION_RUN_H

#import <Foundation/Foundation.h>
#import "Protocols/TTIOIndexable.h"
#import "Protocols/TTIOStreamable.h"
#import "Protocols/TTIOProvenanceable.h"
#import "Protocols/TTIOEncryptable.h"
#import "Protocols/TTIORun.h"
#import "ValueClasses/TTIOEnums.h"

#import "Providers/TTIOStorageProtocols.h"

@class TTIOSpectrum;
@class TTIOMassSpectrum;
@class TTIOChromatogram;
@class TTIOInstrumentConfig;
@class TTIOSpectrumIndex;
@class TTIOValueRange;
@class TTIOProvenanceRecord;
@class TTIOAccessPolicy;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Conforms To:</em> TTIOIndexable, TTIOStreamable,
 * TTIOProvenanceable, TTIOEncryptable, TTIORun</p>
 * <p><em>Declared In:</em> Run/TTIOAcquisitionRun.h</p>
 *
 * <p>An ordered run of spectra sharing an instrument configuration
 * and acquisition mode. The non-genomic counterpart of
 * <code>TTIOGenomicRun</code>; both conform to <code>TTIORun</code>
 * so cross-modality code can iterate uniformly.</p>
 *
 * <p>A run accepts any <code>TTIOSpectrum</code> subclass (mass
 * spectra, NMR spectra, Raman, IR, UV-Vis, ...) but every spectrum
 * within a single run must share a single subclass. Signal-channel
 * serialisation is name-driven, so an MS run writes
 * <code>mz_values</code> + <code>intensity_values</code> and an NMR
 * run writes <code>chemical_shift_values</code> +
 * <code>intensity_values</code>. The run group carries a
 * <code>spectrum_class</code> attribute identifying the subclass.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.acquisition_run.AcquisitionRun</code><br/>
 * Java: <code>global.thalion.ttio.AcquisitionRun</code></p>
 */
@interface TTIOAcquisitionRun : NSObject <TTIOIndexable,
                                          TTIOStreamable,
                                          TTIOProvenanceable,
                                          TTIOEncryptable,
                                          TTIORun>

/** Run identifier as stored in the .tio file (e.g.
 *  <code>@"run_0001"</code>). Defaults to the empty string for
 *  freshly constructed in-memory runs that have not yet been
 *  persisted. */
@property (readonly, copy) NSString *name;

/** Acquisition mode enum value identifying the protocol context. */
@property (readonly) TTIOAcquisitionMode acquisitionMode;

/** Instrument-configuration metadata. */
@property (readonly, strong) TTIOInstrumentConfig *instrumentConfig;

/** Per-spectrum offsets, lengths, and queryable scan metadata. */
@property (readonly, strong) TTIOSpectrumIndex *spectrumIndex;

/** Name of the dominant spectrum class for this run, e.g.
 *  <code>@"TTIOMassSpectrum"</code> or
 *  <code>@"TTIONMRSpectrum"</code>. */
@property (readonly, copy) NSString *spectrumClassName;

/** Omics modality this run carries. Storage attribute
 *  <code>@modality</code> (UTF-8). Defaults to
 *  <code>@"mass_spectrometry"</code>. */
@property (readonly, copy) NSString *modality;

/** Nucleus identifier for NMR runs (zero / <code>nil</code> for
 *  non-NMR). Propagated to every reconstructed
 *  <code>TTIONMRSpectrum</code>. */
@property (readonly, copy) NSString *nucleusType;

/** Spectrometer frequency in MHz for NMR runs. */
@property (readonly) double spectrometerFrequencyMHz;

/** Optional NMR solvent label (e.g. <code>@"CDCl3"</code>,
 *  <code>@"DMSO-d6"</code>). Empty string when not specified or when
 *  the run is not NMR. Stored as the <code>@solvent</code> string
 *  attribute on the run group.
 *
 *  <p>Cross-language equivalents: Java
 *  <code>AcquisitionRun.solvent()</code>, Python
 *  <code>AcquisitionRun.solvent</code>.</p>
 */
@property (readonly, copy) NSString *solvent;

/** Compression codec applied to signal-channel datasets when
 *  persisting this run. Defaults to <code>TTIOCompressionZlib</code>;
 *  writers may set <code>LZ4</code> or <code>NumpressDelta</code>
 *  explicitly before calling <code>-writeToGroup:</code>. */
@property (nonatomic) TTIOCompression signalCompression;

/** Phase 2 of the FLOAT_DELTA_ZSTD spec: MS float64 channels default
 *  to codec id 17 when <code>signalCompression</code> is left at
 *  <code>TTIOCompressionZlib</code>. Set to <code>YES</code> to keep
 *  the chunked-zlib layout — same opt-out pattern as
 *  <code>TTIOWrittenGenomicRun.optDisableInlineMateInfoV2</code>.
 *  Python: <code>WrittenRun.opt_disable_float_delta</code>; Java:
 *  <code>AcquisitionRun.setOptDisableFloatDelta</code>. */
@property (nonatomic) BOOL optDisableFloatDelta;

/** Chromatogram traces associated with this run (TIC / XIC / SRM).
 *  Empty by default. */
@property (readonly, copy) NSArray<TTIOChromatogram *> *chromatograms;

#pragma mark - In-memory construction

/**
 * Convenience initialiser without chromatograms.
 *
 * @param spectra Array of any single <code>TTIOSpectrum</code>
 *                subclass.
 * @param mode    Acquisition mode.
 * @param config  Instrument configuration.
 * @return An initialised run.
 */
- (instancetype)initWithSpectra:(NSArray *)spectra
                acquisitionMode:(TTIOAcquisitionMode)mode
               instrumentConfig:(TTIOInstrumentConfig *)config;

/**
 * Designated initialiser.
 *
 * @param spectra        Array of any single <code>TTIOSpectrum</code>
 *                       subclass.
 * @param chromatograms  Optional chromatograms; pass <code>nil</code>
 *                       or empty array for none.
 * @param mode           Acquisition mode.
 * @param config         Instrument configuration.
 * @return An initialised run.
 */
- (instancetype)initWithSpectra:(NSArray *)spectra
                  chromatograms:(NSArray<TTIOChromatogram *> *)chromatograms
                acquisitionMode:(TTIOAcquisitionMode)mode
               instrumentConfig:(TTIOInstrumentConfig *)config;

#pragma mark - Pixel coordinates (M102)

/**
 * Attaches imaging-grid positions to the run's spectrum index (each
 * <code>int32_t[spectrum count]</code>; three <code>nil</code>s drop
 * them). The index then writes the <code>pixel_x</code> /
 * <code>pixel_y</code> / <code>pixel_z</code> columns and the dataset
 * writer adds <code>opt_pixel_coordinates</code> (format-spec �4b).
 *
 * @return <code>NO</code> with <code>error</code> set when only some
 *         columns are given or a length does not match.
 */
- (BOOL)setPixelX:(NSData *)pixelX
           pixelY:(NSData *)pixelY
           pixelZ:(NSData *)pixelZ
            error:(NSError **)error;

#pragma mark - Storage round-trip

/**
 * Writes this run into a new sub-group named <code>name</code>
 * under <code>parent</code> via the
 * <code>TTIOStorageGroup</code> protocol.
 */
- (BOOL)writeToGroup:(id<TTIOStorageGroup>)parent
                name:(NSString *)name
               error:(NSError **)error;

/**
 * Reads a run from <code>parent/name</code>.
 */
+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)parent
                         name:(NSString *)name
                        error:(NSError **)error;

/**
 * Legacy alias for
 * <code>+readFromGroup:name:error:</code>; identical behaviour,
 * retained for source compatibility.
 */
+ (instancetype)readFromStorageGroup:(id)parent
                                name:(NSString *)name
                               error:(NSError **)error;

#pragma mark - Random access

/**
 * @param index Zero-based position; must satisfy
 *              <code>index &lt; count</code>.
 * @param error Out-parameter populated on failure.
 * @return The materialised spectrum (concrete subclass), or
 *         <code>nil</code> on failure.
 */
- (id)spectrumAtIndex:(NSUInteger)index error:(NSError **)error;

/**
 * Materialise all spectra in this run as an array. Convenience over
 * <code>-spectrumAtIndex:error:</code> + count for callers that want
 * to iterate or stream; rows that fail to materialise are skipped.
 *
 * <p>Cross-language equivalents: Java
 * <code>AcquisitionRun.spectra() : List&lt;Spectrum&gt;</code>, Python
 * <code>AcquisitionRun.spectra() -&gt; list[Spectrum]</code> (also
 * iterable via <code>__iter__</code>).</p>
 *
 * @return Immutable array of materialised spectra in index order.
 */
- (NSArray *)spectra;

/**
 * Element range <code>[start, start + count)</code> of a channel as
 * packed float64 bytes, read without materialising spectra: a decrypted
 * or already-cached column is sliced; a FLOAT_DELTA_ZSTD channel decodes
 * only the blocks the range covers (one block cached per channel); an
 * uncompressed channel is read as a hyperslab. Java
 * <code>AcquisitionRun.channelRange</code>.
 */
- (NSData *)channelRange:(NSString *)channelName
                   start:(NSUInteger)start
                   count:(NSUInteger)count
                   error:(NSError **)error;

/** As -channelRange:start:count:error: with FLOAT_DELTA_ZSTD blocks
 *  decoding on <code>threads</code> workers (0 = TTIO_THREADS; <= 1 is
 *  the serial read). Bytes are identical to the serial read's. */
- (NSData *)channelRange:(NSString *)channelName
                   start:(NSUInteger)start
                   count:(NSUInteger)count
                 threads:(NSUInteger)threads
                   error:(NSError **)error;

/**
 * Iterate the spectra in index order, reading the channels
 * <code>batch</code> spectra at a time through
 * <code>-channelRange:start:count:error:</code>, so a run is walked with
 * bounded memory. Set <code>*stop</code> to end early. Returns NO with
 * <code>error</code> when a read fails.
 */
- (BOOL)iterSpectraWithBatch:(NSUInteger)batch
                       error:(NSError **)error
                  usingBlock:(void (^)(id spectrum, NSUInteger index, BOOL *stop))block;

/** As -iterSpectraWithBatch:error:usingBlock: with each batch's channel
 *  reads decoding on <code>threads</code> workers. */
- (BOOL)iterSpectraWithBatch:(NSUInteger)batch
                     threads:(NSUInteger)threads
                       error:(NSError **)error
                  usingBlock:(void (^)(id spectrum, NSUInteger index, BOOL *stop))block;

/**
 * Visit the run's spectra one scheduling unit at a time.
 *
 * <p>The block runs on several threads at once and in no particular
 * order, so it must be safe to call that way. That relaxed ordering is
 * the whole difference from
 * <code>-iterSpectraWithBatch:threads:error:usingBlock:</code>, whose
 * ordering guarantee is what costs the parallelism.</p>
 *
 * <p>Spectrum <code>k</code> of <code>0 ..&lt; nSpectra</code> is
 * <code>[view spectrumAtIndex:viewStart + k]</code> and its run index
 * is <code>firstSpectrum + k</code>.</p>
 *
 * <p>Setting <code>*stop</code> stops further units being scheduled.
 * Units already in flight still run, so the block may be called after
 * it asks to stop.</p>
 *
 * <p>At one thread this is about 6 per cent slower than
 * <code>-iterSpectraWithBatch:threads:</code>, because a unit is
 * several times a batch and the channel codec's one-block cache serves
 * the smaller reads better. The two cross at two threads. Passing 0
 * resolves the count from <code>TTIO_THREADS</code>.</p>
 *
 * @param from    First spectrum index, inclusive.
 * @param to      Last spectrum index, exclusive.
 * @param threads 0 to resolve from <code>TTIO_THREADS</code>.
 * @param error   Populated on failure.
 * @param block   Visitor.
 * @return NO on a read failure, YES otherwise.
 */
- (BOOL)iterBlocksFrom:(NSUInteger)from
                    to:(NSUInteger)to
               threads:(NSUInteger)threads
                 error:(NSError **)error
            usingBlock:(void (^)(TTIOAcquisitionRun *view,
                                 NSUInteger viewStart,
                                 NSUInteger firstSpectrum,
                                 NSUInteger nSpectra,
                                 BOOL *stop))block;

/**
 * The unit's decoded float64 values for one channel, without building a
 * single spectrum.
 *
 * <p>Only a view handed to an
 * <code>-iterBlocksFrom:to:threads:error:usingBlock:</code> visitor has
 * columns; every other run returns <code>nil</code>.</p>
 *
 * <p>Constructing spectra is 73.7 per cent of the ordered reader's wall
 * time on an Orbitrap Exploris run, and the objects themselves are two
 * thirds of that, so a caller doing bulk arithmetic pays far more for
 * the objects than for the data. Reading the column directly measures
 * 820.7 MB/s against 215.5 through
 * <code>-iterSpectraWithBatch:</code>, both single-threaded, minima of
 * 5 interleaved passes.</p>
 *
 * <p>Spectrum <code>k</code> of the unit covers values
 * <code>[offsetAt:(firstSpectrum + k) - *outValueStart,
 * + lengthAt:(firstSpectrum + k))</code> of the returned buffer.</p>
 *
 * @param channel       Signal channel name.
 * @param outValueStart Run-global value index the buffer starts at;
 *                      may be NULL.
 * @return The unit's values, or nil when the receiver is not a view or
 *         has no such channel.
 */
- (nullable NSData *)unitColumnForChannel:(NSString *)channel
                            valueStart:(nullable unsigned long long *)outValueStart;

/** Write chromatograms under <code>runGroup/chromatograms/</code> in the
 *  layout <code>-writeToGroup:name:error:</code> uses. */
+ (BOOL)writeChromatograms:(NSArray<TTIOChromatogram *> *)chromatograms
                toRunGroup:(id<TTIOStorageGroup>)runGroup
                     error:(NSError **)error;

/** Write per-run provenance (compound <code>provenance/steps</code> on
 *  HDF5, always the <code>@provenance_json</code> mirror). */
+ (BOOL)writeProvenance:(NSArray<TTIOProvenanceRecord *> *)records
             toRunGroup:(id<TTIOStorageGroup>)runGroup
                  error:(NSError **)error;

/**
 * @param range Closed retention-time range in seconds.
 * @return Indices in ascending order whose retention time falls
 *         inside the range.
 */
- (NSArray<NSNumber *> *)indicesInRetentionTimeRange:(TTIOValueRange *)range;

#pragma mark - Persistence context

/**
 * Attaches file-path + run-name context after load so protocol
 * encryption methods can delegate to the in-place encryption
 * manager. Internal API.
 */
- (void)setPersistenceFilePath:(NSString *)path runName:(NSString *)runName;

/**
 * Releases all cached HDF5 handles (group + per-channel datasets).
 * After this call <code>-spectrumAtIndex:error:</code> fails; the
 * run keeps its index metadata so <code>count</code> / headers
 * remain queryable.
 */
- (void)releaseHDF5Handles;

#pragma mark - TTIOProvenanceable

/**
 * Append a processing-step record to the run's provenance chain.
 *
 * Concrete implementation of `-[TTIOProvenanceable
 * addProcessingStep:]` — see that protocol for the full contract.
 *
 * @param step  Record to append (input entities, activity, outputs).
 */
- (void)addProcessingStep:(TTIOProvenanceRecord *)step;

/**
 * Return the run's full provenance chain in insertion order.
 *
 * Concrete implementation of `-[TTIOProvenanceable provenanceChain]`.
 *
 * @return Records in insertion order; empty array if none recorded.
 */
- (NSArray<TTIOProvenanceRecord *> *)provenanceChain;

/**
 * Return every input-entity identifier referenced in the run's chain.
 *
 * Concrete implementation of `-[TTIOProvenanceable inputEntities]`.
 *
 * @return Input-entity identifiers in chain order; duplicates preserved.
 */
- (NSArray<NSString *> *)inputEntities;

/**
 * Return every output-entity identifier referenced in the run's chain.
 *
 * Concrete implementation of `-[TTIOProvenanceable outputEntities]`.
 *
 * @return Output-entity identifiers in chain order; duplicates preserved.
 */
- (NSArray<NSString *> *)outputEntities;

#pragma mark - TTIOEncryptable

/**
 * Encrypt the run in place at the requested granularity.
 *
 * Concrete implementation of `-[TTIOEncryptable
 * encryptWithKey:level:error:]`. Requires that
 * `-setPersistenceFilePath:runName:` has previously bound the run to
 * an on-disk container.
 *
 * @param key    32-byte AES-256-GCM key.
 * @param level  Granularity (`Channel`, `PerAU`, etc.) — see
 *               `TTIOEncryptionLevel`.
 * @param error  Out-error on missing persistence path, wrong key
 *               length, or HDF5 write failure.
 * @return YES on success, NO with `*error` set on failure.
 */
- (BOOL)encryptWithKey:(NSData *)key
                 level:(TTIOEncryptionLevel)level
                 error:(NSError **)error;

/**
 * Decrypt the run in place, restoring plaintext channels.
 *
 * Concrete implementation of `-[TTIOEncryptable
 * decryptWithKey:error:]`. Idempotent — returns YES with no error if
 * the run is already plaintext.
 *
 * @param key    32-byte AES-256-GCM key matching the one used to
 *               encrypt.
 * @param error  Out-error on GCM tag mismatch (wrong key), missing
 *               envelope, or HDF5 write failure.
 * @return YES on success, NO with `*error` set on failure.
 */
- (BOOL)decryptWithKey:(NSData *)key error:(NSError **)error;

/**
 * Return the access policy currently associated with the run.
 *
 * Concrete implementation of `-[TTIOEncryptable accessPolicy]`.
 *
 * @return The active `TTIOAccessPolicy`, or `nil` if none has been
 *         set.
 */
- (TTIOAccessPolicy *)accessPolicy;

/**
 * Replace the access policy associated with the run.
 *
 * Concrete implementation of `-[TTIOEncryptable setAccessPolicy:]`.
 *
 * @param policy  New access policy. Pass `nil` to clear an existing
 *                policy.
 */
- (void)setAccessPolicy:(TTIOAccessPolicy *)policy;

@end

#endif
