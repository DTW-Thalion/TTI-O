#ifndef TTIO_SPECTRUM_INDEX_H
#define TTIO_SPECTRUM_INDEX_H

#import <Foundation/Foundation.h>
#import "ValueClasses/TTIOEnums.h"

#import "Providers/TTIOStorageProtocols.h"

@class TTIOValueRange;
@class TTIOIsolationWindow;

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Run/TTIOSpectrumIndex.h</p>
 *
 * <p>Per-spectrum offsets, lengths, and queryable scan metadata for
 * one <code>TTIOAcquisitionRun</code>. Kept as parallel C arrays
 * inside <code>NSData</code> buffers for compact storage and fast
 * iteration. Persisted as parallel 1-D datasets under the run's
 * <code>spectrum_index/</code> sub-group.</p>
 *
 * <p>Range queries (RT, ms_level, polarity) operate on the
 * in-memory arrays and do not touch the signal channels — this is
 * the compressed-domain query property of the access-unit storage
 * model.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.acquisition_run.SpectrumIndex</code><br/>
 * Java: <code>global.thalion.ttio.SpectrumIndex</code></p>
 */
@interface TTIOSpectrumIndex : NSObject

/** Number of indexed spectra. */
@property (readonly) NSUInteger count;

/** @return Starting element index of spectrum <code>index</code>
 *          within <code>mz_values</code>. */
- (uint64_t)offsetAt:(NSUInteger)index;

/** @return Number of elements (peaks) in spectrum <code>index</code>. */
- (uint32_t)lengthAt:(NSUInteger)index;

/** @return Retention time in seconds. */
- (double)retentionTimeAt:(NSUInteger)index;

/** @return MS level (1, 2, 3, ...). */
- (uint8_t)msLevelAt:(NSUInteger)index;

/** @return Scan polarity. */
- (TTIOPolarity)polarityAt:(NSUInteger)index;

/** @return Precursor m/z; <code>0</code> for MS1. */
- (double)precursorMzAt:(NSUInteger)index;

/** @return Precursor charge state; <code>0</code> if unknown. */
- (uint8_t)precursorChargeAt:(NSUInteger)index;

/** @return Base-peak intensity. */
- (double)basePeakIntensityAt:(NSUInteger)index;

/**
 * @return Activation method for spectrum <code>index</code>;
 *         <code>TTIOActivationMethodNone</code> when the file was
 *         written without the
 *         <code>opt_ms2_activation_detail</code> feature flag.
 */
- (TTIOActivationMethod)activationMethodAt:(NSUInteger)index;

/**
 * @return Isolation window for spectrum <code>index</code>, or
 *         <code>nil</code> when the optional columns are absent or
 *         the stored target+offsets are all zero (MS1 sentinel).
 */
- (TTIOIsolationWindow *)isolationWindowAt:(NSUInteger)index;

/** <code>YES</code> when all four optional activation-detail
 *  parallel columns are present. */
@property (readonly) BOOL hasActivationDetail;

/**
 * @param range Closed retention-time range in seconds.
 * @return Indices whose retention time falls inside the range.
 */
- (NSIndexSet *)indicesInRetentionTimeRange:(TTIOValueRange *)range;

/**
 * @param msLevel MS level to filter on.
 * @return Indices whose MS level matches.
 */
- (NSIndexSet *)indicesForMsLevel:(uint8_t)msLevel;

#pragma mark - Construction

/**
 * Convenience initialiser without optional activation-detail
 * columns. Delegates to the designated initialiser with
 * <code>nil</code> trailing arguments.
 */
- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
            basePeakIntensities:(NSData *)basePeakIntensities;

/**
 * Designated initialiser. The four trailing optional columns
 * (<code>activationMethods</code>,
 * <code>isolationTargetMzs</code>,
 * <code>isolationLowerOffsets</code>,
 * <code>isolationUpperOffsets</code>) must all be <code>nil</code>
 * (legacy file) or all non-<code>nil</code> (activation-detail
 * columns present); mixed <code>nil</code> / non-<code>nil</code>
 * is rejected.
 */
- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
            basePeakIntensities:(NSData *)basePeakIntensities
              activationMethods:(NSData *)activationMethods
             isolationTargetMzs:(NSData *)isolationTargetMzs
          isolationLowerOffsets:(NSData *)isolationLowerOffsets
          isolationUpperOffsets:(NSData *)isolationUpperOffsets;

/**
 * Full initialiser including the optional <code>centroideds</code>
 * column. <code>centroideds</code> may be <code>nil</code> (legacy
 * file, no per-spectrum centroided info) or an <code>int32_t[count]</code>
 * with <code>0 = profile, 1 = centroided</code>.
 */
- (instancetype)initWithOffsets:(NSData *)offsets
                        lengths:(NSData *)lengths
                 retentionTimes:(NSData *)retentionTimes
                       msLevels:(NSData *)msLevels
                     polarities:(NSData *)polarities
                   precursorMzs:(NSData *)precursorMzs
               precursorCharges:(NSData *)precursorCharges
            basePeakIntensities:(NSData *)basePeakIntensities
              activationMethods:(NSData *)activationMethods
             isolationTargetMzs:(NSData *)isolationTargetMzs
          isolationLowerOffsets:(NSData *)isolationLowerOffsets
          isolationUpperOffsets:(NSData *)isolationUpperOffsets
                     centroideds:(NSData *)centroideds;

/**
 * @return Whether spectrum <code>i</code> is centroided (mzML
 * MS:1000127). Returns <code>NO</code> when the centroideds column is
 * absent — callers needing to distinguish "unknown" from "profile"
 * should also check <code>-hasCentroided</code>.
 */
- (BOOL)centroidedAt:(NSUInteger)index;

/** @return YES when the centroideds column is present. */
@property (readonly) BOOL hasCentroided;

#pragma mark - Pixel coordinates (M102)

/**
 * Imaging-grid x position of every spectrum as
 * <code>int32_t[count]</code>, or <code>nil</code> when the run
 * carries no pixel columns (format-spec §4b,
 * <code>opt_pixel_coordinates</code>). <code>pixelX</code>,
 * <code>pixelY</code> and <code>pixelZ</code> are present together
 * or not at all.
 */
@property (readonly, copy) NSData *pixelX;

/** Imaging-grid y positions, <code>int32_t[count]</code> or
 *  <code>nil</code>. */
@property (readonly, copy) NSData *pixelY;

/** Imaging-grid z positions, <code>int32_t[count]</code> or
 *  <code>nil</code>; <code>1</code> for a 2-D acquisition. */
@property (readonly, copy) NSData *pixelZ;

/** <code>YES</code> when the three pixel columns are present. */
@property (readonly) BOOL hasPixelCoordinates;

/** @return Pixel x of spectrum <code>index</code>; <code>0</code>
 *          when the pixel columns are absent. */
- (int32_t)pixelXAt:(NSUInteger)index;

/** @return Pixel y of spectrum <code>index</code>; <code>0</code>
 *          when the pixel columns are absent. */
- (int32_t)pixelYAt:(NSUInteger)index;

/** @return Pixel z of spectrum <code>index</code>; <code>0</code>
 *          when the pixel columns are absent. */
- (int32_t)pixelZAt:(NSUInteger)index;

/**
 * Returns a copy of this index carrying the given pixel columns
 * (each <code>int32_t[count]</code>) in place of any it had. Pass
 * three <code>nil</code>s to drop the columns.
 *
 * @return The new index, or <code>nil</code> with
 *         <code>error</code> set when only some of the columns are
 *         given or a column's length is not <code>count</code>.
 */
- (instancetype)indexWithPixelX:(NSData *)pixelX
                         pixelY:(NSData *)pixelY
                         pixelZ:(NSData *)pixelZ
                          error:(NSError **)error;

#pragma mark - Storage round-trip

/**
 * Writes the index under <code>parent/spectrum_index/</code>.
 *
 * @param parent Destination parent group.
 * @param error  Out-parameter populated on failure.
 * @return <code>YES</code> on success.
 */
- (BOOL)writeToGroup:(id<TTIOStorageGroup>)parent error:(NSError **)error;

/**
 * Reads the index from <code>parent/spectrum_index/</code>. A group
 * holding only some of <code>pixel_x</code> / <code>pixel_y</code> /
 * <code>pixel_z</code> is malformed and fails with
 * <code>TTIOErrorUnsupportedLayout</code>.
 *
 * @param parent Source parent group.
 * @param error  Out-parameter populated on failure.
 * @return The materialised index, or <code>nil</code> on failure.
 */
+ (instancetype)readFromGroup:(id<TTIOStorageGroup>)parent error:(NSError **)error;

/**
 * Legacy alias for <code>+readFromGroup:error:</code>; identical
 * behaviour, retained for source compatibility.
 */
+ (instancetype)readFromStorageGroup:(id)parent error:(NSError **)error;

@end

#endif
