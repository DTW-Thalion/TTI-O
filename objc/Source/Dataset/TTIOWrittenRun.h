/*
 * TTIOWrittenRun.h
 * TTI-O Objective-C Implementation
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#ifndef TTIO_WRITTEN_RUN_H
#define TTIO_WRITTEN_RUN_H

#import <Foundation/Foundation.h>
#import "ValueClasses/TTIOEnums.h"

@class TTIOProvenanceRecord;

NS_ASSUME_NONNULL_BEGIN

/**
 * <p><em>Inherits From:</em> NSObject</p>
 * <p><em>Declared In:</em> Dataset/TTIOWrittenRun.h</p>
 *
 * <p>Flat-buffer value object for
 * <code>+[TTIOSpectralDataset writeMinimalToPath:...]</code>.
 * Mirrors Python's
 * <code>ttio.spectral_dataset.WrittenRun</code> dataclass.</p>
 *
 * <p>The high-level <code>TTIOAcquisitionRun</code> /
 * <code>TTIOMassSpectrum</code> API requires one
 * <code>NSData</code> per spectrum per channel; at write time the
 * writer iterates every spectrum and memcpy-concatenates the
 * per-spectrum buffers into a flat channel. Callers with
 * already-flat arrays can skip both costs by building an
 * <code>TTIOWrittenRun</code> and calling
 * <code>+writeMinimalToPath:</code> on
 * <code>TTIOSpectralDataset</code>.</p>
 *
 * <p>All buffer fields must be little-endian and typed as
 * documented. No validation beyond length; the writer trusts the
 * caller.</p>
 *
 * <p><strong>API status:</strong> Stable.</p>
 *
 * <p><strong>Cross-language equivalents:</strong><br/>
 * Python: <code>ttio.spectral_dataset.WrittenRun</code><br/>
 * Java: <code>global.thalion.ttio.WrittenRun</code></p>
 */
@interface TTIOWrittenRun : NSObject

/** Spectrum class name (e.g. <code>@"TTIOMassSpectrum"</code> or
 *  <code>@"TTIONMRSpectrum"</code>). */
@property (nonatomic, copy, readonly) NSString *spectrumClassName;

/** Acquisition mode raw int value; see
 *  <code>TTIOAcquisitionMode</code>. */
@property (nonatomic, readonly) int64_t acquisitionMode;

/** Flat per-channel signal buffers, keyed by channel name. Each
 *  value is an <code>NSData</code> of float64 little-endian, already
 *  concatenated across spectra in index order. */
@property (nonatomic, copy, readonly) NSDictionary<NSString *, NSData *> *channelData;

/** int64 LE per-spectrum offsets into the channel buffers. */
@property (nonatomic, copy, readonly) NSData *offsets;

/** uint32 LE per-spectrum lengths. */
@property (nonatomic, copy, readonly) NSData *lengths;

/** float64 LE per-spectrum retention times in seconds. */
@property (nonatomic, copy, readonly) NSData *retentionTimes;

/** int32 LE per-spectrum MS levels. */
@property (nonatomic, copy, readonly) NSData *msLevels;

/** int32 LE per-spectrum polarities (cast from
 *  <code>TTIOPolarity</code>). */
@property (nonatomic, copy, readonly) NSData *polarities;

/** float64 LE per-spectrum precursor m/z values. */
@property (nonatomic, copy, readonly) NSData *precursorMzs;

/** int32 LE per-spectrum precursor charges. */
@property (nonatomic, copy, readonly) NSData *precursorCharges;

/** float64 LE per-spectrum base-peak intensities. */
@property (nonatomic, copy, readonly) NSData *basePeakIntensities;

/** Optional NMR run-level nucleus identifier; empty for MS runs. */
@property (nonatomic, copy) NSString *nucleusType;

/** Compression codec identifier for non-genomic-codec channels.
 *  One of <code>@"gzip"</code>, <code>@"none"</code>,
 *  <code>@"lz4"</code>, or <code>@"numpress_delta"</code>. Defaults
 *  to <code>@"gzip"</code>. */
@property (nonatomic, copy) NSString *signalCompression;

/** Phase 2 of the FLOAT_DELTA_ZSTD spec: MS float64 channels default
 *  to codec id 17 when <code>signalCompression</code> is left at
 *  <code>@"gzip"</code>. Set to <code>YES</code> to keep the
 *  chunked-zlib layout. Python:
 *  <code>WrittenRun.opt_disable_float_delta</code>; Java:
 *  <code>AcquisitionRun.setOptDisableFloatDelta</code>. */
@property (nonatomic) BOOL optDisableFloatDelta;

/** Per-run provenance records. Persisted under the run's group as
 *  <code>provenance/steps</code> (compound dataset) plus a legacy
 *  <code>@provenance_json</code> attribute mirror for backward
 *  compatibility. Defaults to an empty array. */
@property (nonatomic, copy) NSArray<TTIOProvenanceRecord *> *provenanceRecords;

/** M102 imaging-grid positions, each <code>int32_t[count]</code>,
 *  written as the <code>spectrum_index/pixel_x</code> /
 *  <code>pixel_y</code> / <code>pixel_z</code> columns (format-spec
 *  §4b). Set all three or none; default <code>nil</code> (no
 *  columns, no <code>opt_pixel_coordinates</code> flag). */
@property (nonatomic, copy, nullable) NSData *pixelX;
/** See <code>pixelX</code>. */
@property (nonatomic, copy, nullable) NSData *pixelY;
/** See <code>pixelX</code>. */
@property (nonatomic, copy, nullable) NSData *pixelZ;

/** <code>YES</code> when any of the pixel columns is set. */
@property (nonatomic, readonly) BOOL hasPixelCoordinates;

/**
 * Designated initialiser.
 *
 * @return An initialised written-run buffer.
 */
- (instancetype)initWithSpectrumClassName:(NSString *)spectrumClassName
                          acquisitionMode:(int64_t)acquisitionMode
                              channelData:(NSDictionary<NSString *, NSData *> *)channelData
                                  offsets:(NSData *)offsets
                                  lengths:(NSData *)lengths
                           retentionTimes:(NSData *)retentionTimes
                                 msLevels:(NSData *)msLevels
                               polarities:(NSData *)polarities
                             precursorMzs:(NSData *)precursorMzs
                         precursorCharges:(NSData *)precursorCharges
                      basePeakIntensities:(NSData *)basePeakIntensities;

@end

NS_ASSUME_NONNULL_END

#endif
