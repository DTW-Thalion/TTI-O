/*
 * TTIOAlignedRead.m
 * TTI-O Objective-C Implementation
 *
 * Class:         TTIOAlignedRead
 * Inherits From: NSObject
 * Conforms To:   NSCopying
 * Declared In:   Genomics/TTIOAlignedRead.h
 *
 * Immutable per-read aligned-sequencing value object. No storage
 * I/O on this class directly — TTIOGenomicRun materialises
 * instances from the signal channels.
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#import "TTIOAlignedRead.h"

@implementation TTIOAlignedRead

- (instancetype)initWithReadName:(NSString *)readName
                      chromosome:(NSString *)chromosome
                        position:(int64_t)position
                  mappingQuality:(uint8_t)mappingQuality
                           cigar:(NSString *)cigar
                        sequence:(NSString *)sequence
                       qualities:(NSData *)qualities
                           flags:(uint32_t)flags
                  mateChromosome:(NSString *)mateChromosome
                    matePosition:(int64_t)matePosition
                  templateLength:(int32_t)templateLength
{
    return [self initWithReadName:readName
                       chromosome:chromosome
                         position:position
                   mappingQuality:mappingQuality
                            cigar:cigar
                         sequence:sequence
                        qualities:qualities
                            flags:flags
                   mateChromosome:mateChromosome
                     matePosition:matePosition
                   templateLength:templateLength
                             tags:@""];
}

- (instancetype)initWithReadName:(NSString *)readName
                      chromosome:(NSString *)chromosome
                        position:(int64_t)position
                  mappingQuality:(uint8_t)mappingQuality
                           cigar:(NSString *)cigar
                        sequence:(NSString *)sequence
                       qualities:(NSData *)qualities
                           flags:(uint32_t)flags
                  mateChromosome:(NSString *)mateChromosome
                    matePosition:(int64_t)matePosition
                  templateLength:(int32_t)templateLength
                            tags:(NSString *)tags
{
    self = [super init];
    if (self) {
        _readName        = [readName copy];
        _chromosome      = [chromosome copy];
        _position        = position;
        _mappingQuality  = mappingQuality;
        _cigar           = [cigar copy];
        _sequence        = [sequence copy];
        _qualities       = [qualities copy];
        _flags           = flags;
        _mateChromosome  = [mateChromosome copy];
        _matePosition    = matePosition;
        _templateLength  = templateLength;
        _tags            = tags ? [tags copy] : @"";
    }
    return self;
}

- (id)copyWithZone:(NSZone *)zone
{
    return self;
}

- (BOOL)isEqual:(id)other
{
    if (self == other) return YES;
    if (![other isKindOfClass:[TTIOAlignedRead class]]) return NO;
    TTIOAlignedRead *o = (TTIOAlignedRead *)other;
    return [_readName isEqualToString:o.readName]
        && [_chromosome isEqualToString:o.chromosome]
        && _position == o.position
        && _mappingQuality == o.mappingQuality
        && [_cigar isEqualToString:o.cigar]
        && [_sequence isEqualToString:o.sequence]
        && [_qualities isEqualToData:o.qualities]
        && _flags == o.flags
        && [_mateChromosome isEqualToString:o.mateChromosome]
        && _matePosition == o.matePosition
        && _templateLength == o.templateLength
        && [_tags isEqualToString:o.tags];
}

- (NSUInteger)hash
{
    return _readName.hash ^ (NSUInteger)_position ^ (NSUInteger)_flags;
}

- (BOOL)isMapped         { return (_flags & 0x4) == 0; }
- (BOOL)isPaired         { return (_flags & 0x1) != 0; }
- (BOOL)isReverse        { return (_flags & 0x10) != 0; }
- (BOOL)isSecondary      { return (_flags & 0x100) != 0; }
- (BOOL)isSupplementary  { return (_flags & 0x800) != 0; }
- (NSUInteger)readLength { return _sequence.length; }

@end
