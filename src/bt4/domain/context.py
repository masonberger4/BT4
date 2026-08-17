"""The sequence around the CDS: 5'UTR, 3' flank, and vector backbone.

Until now BT4 optimized a coding sequence **in isolation**. Nothing in the engine
could see the 5'UTR that sits in front of the CDS in the real transcript, or the
vector backbone the CDS is cloned into. That made a whole class of real defects
invisible: a restriction site created at the insert/backbone junction, a uORF
opened by an ATG in the leader that reads in-frame into the CDS, an initiator
Kozak context that is *half* the user's to choose, and splice models fed a
nucleotide vacuum instead of real flanks.

:class:`ConstructContext` is the value object that carries that sequence. It is
deliberately **pure data**: it holds the flanks, says where the CDS starts, and
answers "what are the last ``k`` bases before the CDS?" -- it enforces nothing.
Constraints stay unchanged; the pipeline wraps them so their existing
``ok_suffix`` sees the junction (see :mod:`bt4.constraints.seeded`).

Three rules make it safe to hand to code that was written for a bare CDS:

* **Unknown bases end the flank.** Real constructs are often only partly known, and
  a placeholder ``N`` is not a base a constraint can reason about. Each flank is
  therefore truncated at the ``N`` nearest the CDS, so what reaches a constraint is
  always contiguous, real, and adjacent to the coding sequence. Nothing is silently
  invented and no ``N`` ever reaches a rule that would mis-handle it.
* **Coordinates are stated, never inferred.** :attr:`cds_offset` is the index of
  the CDS's first base within the assembled construct, so a rule that needs a
  reading frame (uORF pairing, internal-ATG position) can convert between construct
  and CDS coordinates instead of assuming the CDS starts at 0 -- the assumption that
  silently inverts every frame classification once a leader exists.
* **Repeats-by-construction are maskable.** AAV ITRs are a 145-nt palindrome and
  lentiviral LTRs are duplicated by design, so a whole-construct repeat audit over
  those two systems is pure noise unless the user can exclude them.
  :attr:`masked_spans` is that exclusion, in construct coordinates.

This module depends only on the standard library (``domain`` imports nothing).
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["CIRCULAR", "LINEAR", "ConstructContext"]

LINEAR = "linear"
"""A construct with two ends (a PCR product, an mRNA, a linearized vector)."""

CIRCULAR = "circular"
"""A closed construct (a plasmid): its 3' end is adjacent to its 5' end."""

_TOPOLOGIES = frozenset({LINEAR, CIRCULAR})

# Bases a flank may contain. ``N`` is accepted as an explicit "unknown" marker and
# truncates the flank (see the module docstring); every other IUPAC code is
# rejected, because a degenerate flank is not something a constraint can evaluate.
_FLANK_BASES = frozenset("ACGTN")


def _clean(sequence: str, *, keep: str) -> str:
    """Upper-case ``sequence``, validate its alphabet, and truncate it at ``N``.

    Args:
        sequence: Raw flank text (whitespace is ignored, case is not significant).
        keep: ``"suffix"`` to keep the part nearest a downstream CDS (an upstream
            flank), ``"prefix"`` to keep the part nearest an upstream CDS (a
            downstream flank).

    Returns:
        The contiguous, fully-known run of bases adjacent to the CDS.

    Raises:
        ValueError: If the flank contains a character outside ``{A,C,G,T,N}``.
    """
    seq = "".join(sequence.split()).upper()
    bad = sorted(set(seq) - _FLANK_BASES)
    if bad:
        raise ValueError(
            "construct context must be over {A,C,G,T} with N marking unknown bases; "
            f"got unsupported character(s): {bad}"
        )
    if "N" not in seq:
        return seq
    # Keep only the known run touching the CDS: everything past the nearest N is
    # unknown, and guessing across it would be exactly the kind of invented context
    # this object exists to avoid.
    return seq.rsplit("N", 1)[-1] if keep == "suffix" else seq.split("N", 1)[0]


@dataclass(frozen=True, slots=True)
class ConstructContext:
    """The known sequence flanking the CDS being designed.

    Attributes:
        upstream: Sequence immediately 5' of the CDS (5'UTR, and as much of the
            backbone before it as is known). Truncated at the ``N`` nearest the
            CDS, so only the contiguous known run adjacent to the start codon is
            kept.
        downstream: Sequence immediately 3' of the CDS (3'UTR / backbone),
            truncated at its nearest ``N`` the same way.
        topology: ``"linear"`` or ``"circular"``. Recorded so a whole-construct
            audit can treat a plasmid's ends as adjacent; it does not change how
            the CDS itself is optimized.
        masked_spans: Half-open ``(start, end)`` spans in **construct
            coordinates** to exclude from whole-construct repeat/motif reporting --
            for regions that are repeats by design (AAV ITRs, lentiviral LTRs).
        label: Optional free-text name for the construct, carried for display only.
    """

    upstream: str = ""
    downstream: str = ""
    topology: str = LINEAR
    masked_spans: tuple[tuple[int, int], ...] = ()
    label: str = ""
    _normalized: bool = field(default=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Normalize the flanks and validate topology and spans.

        Raises:
            ValueError: On an unsupported base, an unknown topology, or a
                malformed masked span.
        """
        if not self._normalized:
            object.__setattr__(self, "upstream", _clean(self.upstream, keep="suffix"))
            object.__setattr__(self, "downstream", _clean(self.downstream, keep="prefix"))
            object.__setattr__(self, "_normalized", True)
        if self.topology not in _TOPOLOGIES:
            raise ValueError(
                f"topology must be one of {sorted(_TOPOLOGIES)}, got {self.topology!r}"
            )
        for start, end in self.masked_spans:
            if start < 0 or end < start:
                raise ValueError(
                    f"masked span must be a half-open (start, end) with 0 <= start "
                    f"<= end, got ({start}, {end})"
                )

    @property
    def is_empty(self) -> bool:
        """True when no flanking sequence is known (equivalent to no context)."""
        return not self.upstream and not self.downstream

    @property
    def cds_offset(self) -> int:
        """Index of the CDS's first base within the assembled construct."""
        return len(self.upstream)

    def upstream_tail(self, k: int) -> str:
        """Return the last ``k`` known bases before the CDS (shorter if fewer exist).

        This is what a LOCAL constraint needs to judge the 5' junction: it is
        exactly the prefix the trellis would have had, had the CDS not started at
        index 0.

        Args:
            k: Number of trailing bases requested; ``<= 0`` yields ``""``.
        """
        if k <= 0:
            return ""
        return self.upstream[-k:]

    def downstream_head(self, k: int) -> str:
        """Return the first ``k`` known bases after the CDS (shorter if fewer exist)."""
        if k <= 0:
            return ""
        return self.downstream[:k]

    def assemble(self, cds: str) -> str:
        """Return the whole construct: ``upstream + cds + downstream``.

        Used by the whole-construct audit. For a circular construct the caller is
        responsible for handling wrap-around; this returns the linearized form with
        the CDS in the middle, so both junctions are interior to a linear scan.
        """
        return f"{self.upstream}{cds.upper()}{self.downstream}"

    def is_masked(self, start: int, end: int) -> bool:
        """True when ``[start, end)`` (construct coordinates) lies in a masked span."""
        return any(start >= lo and end <= hi for lo, hi in self.masked_spans)
