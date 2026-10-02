"""``FeatureFlags`` — the registry of format feature strings used by a file.

The actual on-disk read/write lives in :mod:`ttio._hdf5_io`; this module
provides a thin value object and the canonical list of reserved feature
names documented in ``docs/feature-flags.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

# v0.2 features (required unless prefixed with ``opt_``)
BASE_V1 = "base_v1"
COMPOUND_IDENTIFICATIONS = "compound_identifications"
COMPOUND_QUANTIFICATIONS = "compound_quantifications"
COMPOUND_PROVENANCE = "compound_provenance"
OPT_COMPOUND_HEADERS = "opt_compound_headers"
OPT_NATIVE_2D_NMR = "opt_native_2d_nmr"
OPT_NATIVE_MSIMAGE_CUBE = "opt_native_msimage_cube"
OPT_DATASET_ENCRYPTION = "opt_dataset_encryption"
OPT_DIGITAL_SIGNATURES = "opt_digital_signatures"

# Reserved for v0.3
OPT_CANONICAL_SIGNATURES = "opt_canonical_signatures"
COMPOUND_PER_RUN_PROVENANCE = "compound_per_run_provenance"

# v0.4 (M25): envelope encryption + rotatable KEK wrapping the DEK.
OPT_KEY_ROTATION = "opt_key_rotation"

# v0.4 (M28): file has been through the anonymization pipeline.
OPT_ANONYMIZED = "opt_anonymized"

# v0.8 (M49): file uses post-quantum crypto — ML-KEM-1024 key wrapping
# and/or ML-DSA-87 (v3:) signatures. Set whenever either primitive is
# used. Opt-flag (reader without PQC can still open the file and read
# unencrypted datasets; it just cannot verify v3: signatures or unwrap
# ML-KEM-wrapped DEKs).
OPT_PQC_PREVIEW = "opt_pqc_preview"

# per-AU encryption of signal channels — one AES-GCM op per
# spectrum per channel. Replaces the v0.x channel-granular
# encryption layout with ``<channel>_segments`` compound datasets.
# See ``docs/format-spec.md`` §9.1 and
# ``docs/transport-encryption-design.md``.
OPT_PER_AU_ENCRYPTION = "opt_per_au_encryption"

# additionally encrypts the AU semantic header fields
# (ms_level, polarity, RT, precursor m/z, etc.). Requires
# OPT_PER_AU_ENCRYPTION. Server-side filtering is disabled when
# this flag is set.
OPT_ENCRYPTED_AU_HEADERS = "opt_encrypted_au_headers"

# v0.11.1 (M73.1): native storage of a 2D-COS spectrum pair
# (synchronous + asynchronous correlation matrices) as two rank-2
# HDF5 datasets sharing the spectral-variable axis scale.
OPT_NATIVE_2D_COS = "opt_native_2d_cos"

# v0.12.0 (M74): spectrum_index carries four optional parallel columns
# recording MS/MS activation method and precursor isolation window
# (activation_methods, isolation_target_mzs, isolation_lower_offsets,
# isolation_upper_offsets). Files that set this flag are written with
# format version "1.3"; pre-M74 readers that ignore the flag still
# parse the base columns.
OPT_MS2_ACTIVATION_DETAIL = "opt_ms2_activation_detail"

# M102: spectrum_index carries the pixel_x / pixel_y / pixel_z int32
# columns of an imaging run (format-spec §4b). Added only when a run
# carries them.
OPT_PIXEL_COORDINATES = "opt_pixel_coordinates"

@dataclass(frozen=True, slots=True)
class FeatureFlags:
    """Immutable set of format feature strings with a version label."""

    version: str = "1.1"
    features: tuple[str, ...] = field(default_factory=tuple)

    @classmethod
    def from_iterable(cls, version: str, features: Iterable[str]) -> "FeatureFlags":
        """Build a :class:`FeatureFlags` from any iterable of feature strings.

        Parameters
        ----------
        version : str
            Format version label (e.g. ``"1.1"``, ``"1.5"``).
        features : Iterable[str]
            Feature flag names; consumed once and frozen into a tuple
            preserving iteration order.

        Returns
        -------
        FeatureFlags
            Immutable feature flag set with the given version label.
        """
        return cls(version=version, features=tuple(features))

    def has(self, name: str) -> bool:
        """Return whether a feature flag is set.

        Parameters
        ----------
        name : str
            Feature flag string to test (e.g.
            ``"opt_ms2_activation_detail"``).

        Returns
        -------
        bool
            True when ``name`` is present in :attr:`features`.
        """
        return name in self.features

    def required(self) -> tuple[str, ...]:
        """Features whose absence should cause a reader to refuse the file."""
        return tuple(f for f in self.features if not f.startswith("opt_"))
