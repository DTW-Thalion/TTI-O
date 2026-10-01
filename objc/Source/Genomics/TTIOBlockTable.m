/*
 * TTI-O Objective-C Implementation
 * Copyright (c) 2026 The Thalion Initiative
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#import "Genomics/TTIOBlockTable.h"
#import "Genomics/TTIOGenomicBlocks.h"
#import "HDF5/TTIOHDF5Errors.h"

@implementation TTIOBlockTable {
    NSUInteger _n;
    BOOL _hasCodecs;
    unsigned long long *_readStart;
    NSUInteger *_nReads;
    unsigned long long *_baseStart;
    unsigned long long *_nBases;
    NSMutableDictionary<NSString *, NSMutableData *> *_off;
    NSMutableDictionary<NSString *, NSMutableData *> *_len;
    NSMutableDictionary<NSString *, NSMutableData *> *_codec;
    NSArray<NSString *> *_channels;
}

static BOOL ttioRowNum(NSDictionary *row, NSString *key, unsigned long long *out, NSError **error)
{
    id v = row[key];
    if (v == nil) {
        if (error) *error = TTIOMakeError(TTIOErrorDatasetRead, @"blocks/index: missing column %@", key);
        return NO;
    }
    *out = [v unsignedLongLongValue];
    return YES;
}

+ (instancetype)readFromRunGroup:(id<TTIOStorageGroup>)runGroup error:(NSError **)error
{
    id<TTIOStorageGroup> blocks = [runGroup openGroupNamed:@"blocks" error:error];
    if (!blocks) return nil;
    id<TTIOStorageDataset> ds = [blocks openDatasetNamed:@"index" error:error];
    if (!ds) return nil;
    NSArray *rows = [ds readRows:error];
    if (!rows) return nil;
    TTIOBlockTable *t = [[TTIOBlockTable alloc] init];
    if (![t _loadRows:rows error:error]) return nil;
    return t;
}

- (BOOL)_loadRows:(NSArray<NSDictionary *> *)rows error:(NSError **)error
{
    _n = rows.count;
    NSArray *channels = [TTIOGenomicBlocks blockChannels];
    NSArray *required = [TTIOGenomicBlocks requiredBlockChannels];
    _hasCodecs = _n > 0 && rows[0][[channels[0] stringByAppendingString:@"_codec"]] != nil;
    // An optional channel's columns are absent from runs that do not
    // carry it: every block then has an empty range (M101).
    NSMutableArray<NSString *> *present = [NSMutableArray arrayWithArray:required];
    for (NSString *ch in [TTIOGenomicBlocks optionalBlockChannels]) {
        if (_n > 0 && rows[0][[ch stringByAppendingString:@"_off"]] != nil) [present addObject:ch];
    }
    _channels = [present copy];
    _readStart = calloc(MAX(_n, (NSUInteger)1), sizeof(unsigned long long));
    _nReads = calloc(MAX(_n, (NSUInteger)1), sizeof(NSUInteger));
    _baseStart = calloc(MAX(_n, (NSUInteger)1), sizeof(unsigned long long));
    _nBases = calloc(MAX(_n, (NSUInteger)1), sizeof(unsigned long long));
    _off = [NSMutableDictionary dictionary];
    _len = [NSMutableDictionary dictionary];
    _codec = [NSMutableDictionary dictionary];
    for (NSString *ch in channels) {
        _off[ch] = [NSMutableData dataWithLength:_n * sizeof(unsigned long long)];
        _len[ch] = [NSMutableData dataWithLength:_n * sizeof(unsigned long long)];
        if (_hasCodecs) _codec[ch] = [NSMutableData dataWithLength:_n * sizeof(NSUInteger)];
    }
    for (NSUInteger i = 0; i < _n; i++) {
        NSDictionary *r = rows[i];
        unsigned long long v;
        if (!ttioRowNum(r, @"read_start", &v, error)) return NO; _readStart[i] = v;
        if (!ttioRowNum(r, @"n_reads", &v, error)) return NO; _nReads[i] = (NSUInteger)v;
        if (!ttioRowNum(r, @"base_start", &v, error)) return NO; _baseStart[i] = v;
        if (!ttioRowNum(r, @"n_bases", &v, error)) return NO; _nBases[i] = v;
        for (NSString *ch in channels) {
            if (![_channels containsObject:ch]) continue;
            if (!ttioRowNum(r, [ch stringByAppendingString:@"_off"], &v, error)) return NO;
            ((unsigned long long *)_off[ch].mutableBytes)[i] = v;
            if (!ttioRowNum(r, [ch stringByAppendingString:@"_len"], &v, error)) return NO;
            ((unsigned long long *)_len[ch].mutableBytes)[i] = v;
            if (_hasCodecs) {
                if (!ttioRowNum(r, [ch stringByAppendingString:@"_codec"], &v, error)) return NO;
                ((NSUInteger *)_codec[ch].mutableBytes)[i] = (NSUInteger)v;
            }
        }
    }
    return YES;
}

- (void)dealloc
{
    free(_readStart);
    free(_nReads);
    free(_baseStart);
    free(_nBases);
}

- (NSUInteger)count { return _n; }
- (BOOL)hasCodecs { return _hasCodecs; }
- (NSArray<NSString *> *)channels { return _channels ?: [TTIOGenomicBlocks requiredBlockChannels]; }
- (BOOL)hasTags { return [_channels containsObject:@"tags"]; }

- (unsigned long long)readCount
{
    return _n == 0 ? 0 : _readStart[_n - 1] + _nReads[_n - 1];
}

- (unsigned long long)readStartAt:(NSUInteger)b { return _readStart[b]; }
- (NSUInteger)nReadsAt:(NSUInteger)b { return _nReads[b]; }
- (unsigned long long)baseStartAt:(NSUInteger)b { return _baseStart[b]; }
- (unsigned long long)nBasesAt:(NSUInteger)b { return _nBases[b]; }

- (unsigned long long)offsetOf:(NSString *)channel at:(NSUInteger)b
{
    NSMutableData *d = _off[channel];
    return d ? ((const unsigned long long *)d.bytes)[b] : 0;
}

- (unsigned long long)lengthOf:(NSString *)channel at:(NSUInteger)b
{
    NSMutableData *d = _len[channel];
    return d ? ((const unsigned long long *)d.bytes)[b] : 0;
}

- (NSUInteger)codecOf:(NSString *)channel at:(NSUInteger)b
{
    NSMutableData *d = _codec[channel];
    return d ? ((const NSUInteger *)d.bytes)[b] : 0;
}

- (NSUInteger)blockForRead:(unsigned long long)i
{
    if (_n == 0) return NSNotFound;
    NSUInteger lo = 0, hi = _n;
    while (lo < hi) {
        NSUInteger mid = lo + (hi - lo) / 2;
        if (_readStart[mid] <= i) lo = mid + 1; else hi = mid;
    }
    if (lo == 0) return NSNotFound;
    NSUInteger b = lo - 1;
    if (i >= _readStart[b] + _nReads[b]) return NSNotFound;
    return b;
}

@end
