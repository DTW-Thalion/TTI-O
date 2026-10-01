#import <Foundation/Foundation.h>
#import "Codecs/TTIOReferenceResolver.h"
NS_ASSUME_NONNULL_BEGIN

/** Run-derived context for codecs. All fields optional/nullable; plain codecs ignore it. */
@interface TTIOCodecContext : NSObject
@property (nullable, strong) NSArray<NSNumber *> *readLengths;     // fqzcomp
@property (nullable, strong) NSArray<NSNumber *> *revcompFlags;    // fqzcomp
@property (nullable, strong) NSNumber *elementSize;               // delta encode
@property (nullable, strong) NSNumber *readCount;
@property (nullable, strong) NSData *positions;                  // int64-LE
@property (nullable, copy)   NSArray<NSString *> *(^cigarsProvider)(void);  // lazy thunk
@property (nullable, strong) NSNumber *totalBases;
@property (nullable, strong) NSArray<NSString *> *chromosomes;
@property (nullable, strong) NSData *ownChromIds;                // mate_info
@property (nullable, strong) NSData *ownPositions;               // mate_info
@property (nullable, strong) NSNumber *nRecords;
@property (nullable, strong) TTIOReferenceResolver *referenceResolver;
// encode-only (ref_diff):
@property (nullable, strong) NSData *offsets;
@property (nullable, strong) NSData *reference;
@property (nullable, strong) NSData *referenceMd5;
@property (nullable, strong) NSString *referenceUri;
@property (nullable, strong) NSNumber *readsPerSlice;
/** REF_DIFF_V2 byte budget: a slice closes before the read that
 *  would push it past this many bases (readsPerSlice still caps the
 *  read count). nil or 0 = the fixed reads-per-slice rule. */
@property (nullable, strong) NSNumber *sliceBytes;

// fqzcomp V5 (sequence context). Decode: lazy block returning the
// run's decoded sequences bytes, invoked only for version-5 streams.
// Encode: the flat base bytes, set by the writer when the run carries
// a base-parallel sequences channel and V5 is not opted out.
@property (nullable, copy)   NSData * _Nullable (^sequencesProvider)(void);
@property (nullable, strong) NSData *sequences;
/** fqzcomp encode strategy: nil/-1 auto, 0..4 V4 preset, 5/6 forced
 *  V5, TTIOM94ZHintV4Auto V4 with internal preset selection. */
@property (nullable, strong) NSNumber *qualStrategyHint;

// SAM_TAGS (M101): the reference bases per ownChromIds value, for
// MD/NM derivation (NSData, or NSNull where there is none). Encode:
// the array itself (nil/empty = no derivation). Decode: a lazy block,
// called once per decode. The other derivation inputs are
// sequences / sequencesProvider, offsets / readLengths, cigarsProvider,
// positions and ownChromIds.
@property (nullable, copy)   NSArray *tagReferences;
@property (nullable, copy)   NSArray * _Nullable (^tagReferencesProvider)(void);

+ (instancetype)emptyContext;
@end

NS_ASSUME_NONNULL_END
