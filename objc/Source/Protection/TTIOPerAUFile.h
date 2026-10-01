/*
 * TTIOPerAUFile.h
 * TTI-O Objective-C Implementation
 *
 * File-level per-AU encryption orchestrator. Reads plaintext
 * <channel>_values datasets, encrypts each spectrum independently
 * with TTIOPerAUEncryption, and rewrites the file's signal_channels
 * groups with the <channel>_segments compound layout from
 * docs/format-spec.md §9.1. Routes through the TTIOStorageProvider
 * / TTIOStorageGroup abstraction so any backend that supports
 * VL_BYTES compound fields (HDF5, Memory, ...) works.
 *
 * Cross-language equivalents:
 *   Python: ttio.encryption_per_au.encrypt_per_au / decrypt_per_au
 *   Java:   global.thalion.ttio.protection.PerAUFile
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * Copyright (c) 2026 The Thalion Initiative
 */
#ifndef TTIO_PER_AU_FILE_H
#define TTIO_PER_AU_FILE_H

#import <Foundation/Foundation.h>

NS_ASSUME_NONNULL_BEGIN

@interface TTIOPerAUFile : NSObject

/** Encrypt ``path`` in place with per-AU AES-256-GCM.
 *
 *  For each MS run, reads plaintext <channel>_values, slices by the
 *  run's spectrum_index/{offsets,lengths}, encrypts each row with a
 *  fresh IV + AAD = dataset_id || au_sequence || channel_name, and
 *  writes a <channel>_segments compound dataset. The plaintext
 *  <channel>_values dataset is deleted.
 *
 *  When ``encryptHeaders`` is YES, also encrypts the six semantic
 *  index arrays (retention_times, ms_levels, polarities,
 *  precursor_mzs, precursor_charges, base_peak_intensities) into
 *  spectrum_index/au_header_segments and deletes the plaintext
 *  arrays. Offsets + lengths remain plaintext (structural framing,
 *  not semantic PHI).
 *
 *  Sets the opt_per_au_encryption feature flag on the root group
 *  (and opt_encrypted_au_headers when header encryption is used).
 *
 *  ``providerName`` may be nil (scheme-based routing via
 *  TTIOProviderRegistry); passing "hdf5" forces the HDF5 provider.
 */
+ (BOOL)encryptFilePath:(NSString *)path
                     key:(NSData *)key
         encryptHeaders:(BOOL)encryptHeaders
            providerName:(nullable NSString *)providerName
                   error:(NSError * _Nullable *)error;

/** Read-only decrypt: return ``{run_name: {channel_name: NSData
 *  float64 LE}}`` for a per-AU-encrypted file. When
 *  opt_encrypted_au_headers is set, the run map also carries
 *  ``__au_headers__`` as an NSArray<TTIOAUHeaderPlaintext *>.
 *  A genomic run that carries SAM tags (M101) maps ``tags`` to an
 *  NSArray<NSString *>, one string per read.
 */
+ (nullable NSDictionary<NSString *, NSDictionary *> *)
    decryptFilePath:(NSString *)path
                key:(NSData *)key
       providerName:(nullable NSString *)providerName
              error:(NSError * _Nullable *)error;

/** Per-AU AES-256-GCM decrypt-in-place. The persist-to-disk counterpart to
 *  ``+encryptFilePath:`` -- mirrors the ``intensity_values_encrypted``
 *  legacy decrypt-in-place (``+[TTIOSpectralDataset decryptInPlaceAtPath:]``)
 *  but for the per-AU compound layout (``<channel>_segments``).
 *
 *  For each MS run with ``<channel>_segments`` under ``signal_channels``,
 *  decrypts each spectrum row with the per-AU GCM scheme (AAD =
 *  ``dataset_id || au_sequence || channel_name``), writes the concatenated
 *  plaintext back as ``<channel>_values`` (float64), deletes
 *  ``<channel>_segments``, and removes the ``<channel>_algorithm`` attribute.
 *
 *  Genomic runs are handled the same way: ``study/genomic_runs/<name>/
 *  signal_channels/<sequences|qualities>_segments`` is decrypted (uint8)
 *  and written back as the bare ``<sequences|qualities>`` dataset.
 *  ``dataset_id`` continues from where the MS loop left off, mirroring
 *  ``encryptFilePath:``'s AAD numbering exactly.
 *
 *  When the file carries ``opt_encrypted_au_headers``, the six plaintext
 *  index datasets (``retention_times``, ``ms_levels``, ``polarities``,
 *  ``precursor_mzs``, ``precursor_charges``, ``base_peak_intensities``)
 *  are restored from ``au_header_segments`` and the encrypted compound is
 *  removed.
 *
 *  Clears the ``opt_per_au_encryption`` / ``opt_encrypted_au_headers``
 *  feature flags and the root ``@encrypted`` attribute on completion, so
 *  subsequent readers see a fully unprotected file. Idempotent: returns
 *  YES with no error when the file is already plaintext.
 *
 *  Cross-language equivalents:
 *    Python: ``ttio.encryption_per_au.decrypt_per_au_in_place``
 *    Java:   ``PerAUFile.decryptFileInPlace``
 */
+ (BOOL)decryptFilePathInPlace:(NSString *)path
                            key:(NSData *)key
                   providerName:(nullable NSString *)providerName
                          error:(NSError * _Nullable *)error;

#pragma mark - Region-based per-AU encryption

/** Encrypt genomic signal channels with a per-chromosome key map.
 *
 *  Reads whose chromosome appears in ``keyMap`` are AES-256-GCM
 *  encrypted with the corresponding 32-byte key. Reads on chromosomes
 *  NOT in ``keyMap`` are stored as "clear segments" — the same
 *  ``<channel>_segments`` compound is reused, but with a length-0 IV
 *  / length-0 tag and the raw plaintext bytes stored in the
 *  ``ciphertext`` slot. The decoder branches on ``len(seg.iv)``, so
 *  uniform-IV files (every IV is exactly 12 bytes) decode unchanged
 *  under the same code path.
 *
 *  When ``keyMap`` contains the reserved key ``"_headers"``,
 *  the four genomic_index columns (chromosomes, positions,
 *  mapping_qualities, flags) are ALSO encrypted with that key —
 *  closing the gap where a reader without any signal-channel key
 *  could still see read locations + counts. ``offsets`` and
 *  ``lengths`` always stay plaintext (structural framing, not
 *  semantic PHI). The presence of the ``"_headers"`` entry is the
 *  opt-in signal for ``opt_encrypted_au_headers``.
 *
 *  MS runs are NOT touched — chromosome is a genomic concept. Use
 *  ``encryptFilePath:`` for MS encryption.
 *
 *  Sets the ``opt_per_au_encryption`` and
 *  ``opt_region_keyed_encryption`` feature flags on the root group
 *  (and ``opt_encrypted_au_headers`` when the ``"_headers"`` key is
 *  used). The per-channel ``<channel>_algorithm`` attribute is set
 *  to ``"aes-256-gcm-by-region"``.
 */
+ (BOOL)encryptFilePathByRegion:(NSString *)path
                          keyMap:(NSDictionary<NSString *, NSData *> *)keyMap
                    providerName:(nullable NSString *)providerName
                           error:(NSError * _Nullable *)error;

/** Decrypt a region-encrypted file using a per-chromosome key map.
 *
 *  Caller may supply only a subset of the keys used at encryption
 *  time. Clear segments (length-0 IV) decode without any key.
 *  Encrypted segments whose chromosome key isn't in ``keyMap``
 *  cause the call to fail (no key available).
 *
 *  Returns ``{run_name: {channel_name: NSData uint8}}`` like
 *  ``decryptFilePath:`` (genomic dtype is uint8, one ASCII byte per
 *  base). MS runs are not part of region-keyed encryption; this
 *  helper does not touch them.
 */
+ (nullable NSDictionary<NSString *, NSDictionary *> *)
    decryptFilePathByRegion:(NSString *)path
                      keyMap:(NSDictionary<NSString *, NSData *> *)keyMap
                providerName:(nullable NSString *)providerName
                       error:(NSError * _Nullable *)error;

@end

NS_ASSUME_NONNULL_END

#endif
