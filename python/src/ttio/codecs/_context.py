"""Value objects for the codec registry (codec-registry refactor).

CodecContext carries everything any TTI-O codec might need; ChannelPayload
hides dataset-vs-group storage; DecodedChannel/EncodedChannel are closed
tagged unions so one Codec.decode/encode signature covers heterogeneous
codecs portably. No wire/format change — these only restructure dispatch.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from ..genomic.reference_resolver import ReferenceResolver
    from ..providers.base import StorageGroup


@dataclass(frozen=True)
class CodecContext:
    """Run-derived context for codecs. All fields optional; plain codecs ignore it."""
    read_lengths: "np.ndarray | None" = None        # fqzcomp; == index.lengths
    revcomp_flags: "np.ndarray | None" = None        # fqzcomp; (flags & 16) != 0
    # fqzcomp V5 (sequence context). Decode: lazy provider returning the
    # run's decoded sequences bytes, called only for version-5 streams.
    # Encode: the flat base bytes, populated by the writer when the run
    # carries a base-parallel sequences channel and V5 is not opted out.
    sequences_provider: "Callable[[], bytes] | None" = None
    sequences: "bytes | None" = None
    # fqzcomp encode strategy: -1 auto, 0..4 V4 preset, 5/6 forced V5,
    # 7 V4 with internal preset selection (fqzcomp_nx16_z.HINT_V4_AUTO).
    qual_strategy_hint: int = -1
    element_size: "int | None" = None                # delta_rans encode
    read_count: "int | None" = None                  # == index.count
    positions: "np.ndarray | None" = None            # ref_diff
    cigars_provider: "Callable[[], list[str]] | None" = None  # ref_diff (lazy)
    total_bases: "int | None" = None                 # ref_diff
    chromosomes: "list[str] | None" = None           # ref_diff
    own_chrom_ids: "np.ndarray | None" = None         # mate_info
    own_positions: "np.ndarray | None" = None         # mate_info
    n_records: "int | None" = None                    # mate_info
    reference_resolver: "ReferenceResolver | None" = None  # ref_diff (decode)
    # ref_diff ENCODE-only inputs (Task 5c). The decode path resolves
    # the reference lazily via reference_resolver + the blob header; the
    # encode path needs them up front because they are written *into*
    # the blob header. All optional — only the writer populates them.
    offsets: "np.ndarray | None" = None               # ref_diff encode; uint64, n_reads+1
    cigar_strings: "list[str] | None" = None          # ref_diff encode; per-read CIGAR
    reference: "bytes | None" = None                  # ref_diff encode; chrom sequence
    reference_md5: "bytes | None" = None              # ref_diff encode; 16 raw bytes
    reference_uri: "str | None" = None                # ref_diff encode
    reads_per_slice: "int | None" = None              # ref_diff encode; slice granularity
    # ref_diff encode; slice byte budget (M97): a slice closes before
    # the read that would push it past this many bases. None/0 = the
    # fixed reads_per_slice rule.
    slice_bytes: "int | None" = None
    # SAM_TAGS (M101): the reference bases per own_chrom_ids value, for
    # MD/NM derivation. Encode: the list itself (empty or None = no
    # derivation). Decode: a lazy provider, called once per decode.
    tag_references: "list[bytes | None] | None" = None
    tag_references_provider: "Callable[[], list[bytes | None]] | None" = None

    @staticmethod
    def empty() -> "CodecContext":
        return CodecContext()


@dataclass(frozen=True)
class ChannelPayload:
    """Encoded payload: either flat dataset bytes or a storage group (ref_diff)."""
    _bytes: "bytes | None" = None
    _group: "StorageGroup | None" = None

    @classmethod
    def of_bytes(cls, b: bytes) -> "ChannelPayload":
        return cls(_bytes=b)

    @classmethod
    def of_group(cls, g: "StorageGroup") -> "ChannelPayload":
        return cls(_group=g)

    def as_bytes(self) -> bytes:
        if self._bytes is None:
            raise TypeError("ChannelPayload holds a group, not bytes")
        return self._bytes

    def group(self) -> "StorageGroup":
        if self._group is None:
            raise TypeError("ChannelPayload holds bytes, not a group")
        return self._group


@dataclass(frozen=True)
class DecodedChannel:
    """Closed union of a decoded channel value: bytes | list[str] | mate-info dict."""
    _bytes: "bytes | None" = None
    _str_list: "list[str] | None" = None
    _mate_info: "dict[str, Any] | None" = None
    _kind: str = ""

    @classmethod
    def of_bytes(cls, b: bytes) -> "DecodedChannel":
        return cls(_bytes=b, _kind="bytes")

    @classmethod
    def of_str_list(cls, s: "list[str]") -> "DecodedChannel":
        return cls(_str_list=s, _kind="str_list")

    @classmethod
    def of_mate_info(cls, d: "dict[str, Any]") -> "DecodedChannel":
        return cls(_mate_info=d, _kind="mate_info")

    def as_bytes(self) -> bytes:
        if self._kind != "bytes":
            raise TypeError(f"DecodedChannel is {self._kind!r}, not bytes")
        return self._bytes  # type: ignore[return-value]

    def as_str_list(self) -> "list[str]":
        if self._kind != "str_list":
            raise TypeError(f"DecodedChannel is {self._kind!r}, not str_list")
        return self._str_list  # type: ignore[return-value]

    def as_mate_info(self) -> "dict[str, Any]":
        if self._kind != "mate_info":
            raise TypeError(f"DecodedChannel is {self._kind!r}, not mate_info")
        return self._mate_info  # type: ignore[return-value]


@dataclass(frozen=True)
class EncodedChannel:
    """Closed union of encode output: a dataset byte blob or a group layout (ref_diff)."""
    dataset_bytes: "bytes | None" = None
    group_children: "dict[str, bytes] | None" = None
    group_attrs: "dict[str, Any] | None" = None
    is_group: bool = False

    @classmethod
    def of_dataset(cls, b: bytes) -> "EncodedChannel":
        return cls(dataset_bytes=b, is_group=False)

    @classmethod
    def of_group(cls, children: "dict[str, bytes]", attrs: "dict[str, Any]") -> "EncodedChannel":
        return cls(group_children=children, group_attrs=attrs, is_group=True)
