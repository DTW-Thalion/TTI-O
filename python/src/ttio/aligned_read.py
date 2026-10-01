"""AlignedRead — one aligned sequencing read.

Genomic analogue of :class:`ttio.mass_spectrum.MassSpectrum`. Frozen
value object materialised by :meth:`ttio.genomic_run.GenomicRun.__getitem__`
from the signal channel arrays under
``/study/genomic_runs/<name>/signal_channels/``.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class AlignedRead:
    """One aligned sequencing read."""

    read_name: str
    chromosome: str
    position: int
    mapping_quality: int
    cigar: str
    sequence: str
    qualities: bytes
    flags: int
    mate_chromosome: str
    mate_position: int
    template_length: int
    # M101: the SAM optional fields (columns 12+), tab-joined as
    # samtools prints them; "" when the read has none.
    tags: str = ""

    @property
    def is_mapped(self) -> bool:
        """Whether the read is aligned to a reference position.

        Returns
        -------
        bool
            True when SAM flag 0x4 (``unmapped``) is cleared.
        """
        return not (self.flags & 0x4)

    @property
    def is_paired(self) -> bool:
        """Whether the read is part of a paired-end template.

        Returns
        -------
        bool
            True when SAM flag 0x1 (``paired``) is set.
        """
        return bool(self.flags & 0x1)

    @property
    def is_reverse(self) -> bool:
        """Whether the read aligns to the reverse strand.

        Returns
        -------
        bool
            True when SAM flag 0x10 (``reverse complement``) is set.
        """
        return bool(self.flags & 0x10)

    @property
    def is_secondary(self) -> bool:
        """Whether the alignment is a secondary alignment.

        Returns
        -------
        bool
            True when SAM flag 0x100 (``secondary alignment``) is set.
        """
        return bool(self.flags & 0x100)

    @property
    def is_supplementary(self) -> bool:
        """Whether the alignment is a supplementary alignment.

        Returns
        -------
        bool
            True when SAM flag 0x800 (``supplementary alignment``) is
            set; typically marks chimeric / split-read parts.
        """
        return bool(self.flags & 0x800)

    @property
    def read_length(self) -> int:
        """Length of the read's sequence in bases.

        Returns
        -------
        int
            ``len(self.sequence)``. Zero when the BAM record stored
            ``SEQ=*`` (sequence stripped).
        """
        return len(self.sequence)
