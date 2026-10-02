/*
 * TTIOWrittenRun.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOWrittenRun
 * Inherits From: NSObject
 * Declared In:   Dataset/TTIOWrittenRun.h
 *
 * Flat-buffer write-side value object for the
 * +[TTIOSpectralDataset writeMinimalToPath:...] fast paths.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOWrittenRun.h"

@implementation TTIOWrittenRun

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
                      basePeakIntensities:(NSData *)basePeakIntensities
{
    self = [super init];
    if (!self) return nil;
    _spectrumClassName = [spectrumClassName copy];
    _acquisitionMode = acquisitionMode;
    _channelData = [channelData copy];
    _offsets = [offsets copy];
    _lengths = [lengths copy];
    _retentionTimes = [retentionTimes copy];
    _msLevels = [msLevels copy];
    _polarities = [polarities copy];
    _precursorMzs = [precursorMzs copy];
    _precursorCharges = [precursorCharges copy];
    _basePeakIntensities = [basePeakIntensities copy];
    _nucleusType = @"";
    _signalCompression = @"gzip";
    _provenanceRecords = @[];
    return self;
}

- (BOOL)hasPixelCoordinates
{
    return _pixelX != nil || _pixelY != nil || _pixelZ != nil;
}

@end
