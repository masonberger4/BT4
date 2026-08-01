"""%MinMax codon-usage shaping — an honest split between a DP term and a metric.

%MinMax (Clarke & Clark 2008) is a sliding-window measure of how common or rare
the chosen synonymous codons are, relative to the max / min / average synonymous
usage of each residue. Over a window of ``W`` consecutive codons it compares the
actual codon frequencies (``Xactual``) to the sums of the most frequent
(``Xmax``), least frequent (``Xmin``), and average (``Xavg``) synonymous
frequencies, yielding a value in ``[-100, +100]``: positive when the window
leans on common codons, negative when it leans on rare ones.

That windowed statistic is genuinely **non-additive** per codon — a single
codon participates in up to ``W`` overlapping windows, and each window is
normalized by its own ``Xmax``/``Xmin`` — so folding it into the trellis DP as a
per-codon ``delta`` would repeat exactly the BT3 sin CLAUDE.md §5.4 / §10.3
forbids (a window statistic re-scored per codon, ``delta != score``). This
module therefore keeps two things strictly separate and honest about each:

* :func:`min_max_profile` is the **reporting metric**: the classic per-window
  %MinMax profile. It is a deterministic, read-only helper for display and
  auditing; it does **not** participate in the DP.
* :class:`MinMaxTerm` is the **DP-participating objective term**: a genuinely
  additive *codon-commonness prior*. Its per-codon contribution reads only the
  codon itself — ``f(codon) - f_avg(aa)``, how much more common than its
  synonymous average the codon is — so its ``delta`` sum equals its
  whole-sequence ``score`` exactly (invariant #4, "delta == score"). It is a
  cheap prior in the same spirit as %MinMax, not the windowed statistic itself.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from bt4.domain.genetic_code import CODON_TABLE, STOP, synonymous_codons
from bt4.domain.scope import Scope
from bt4.objectives.base import iter_codons

__all__ = ["MinMaxTerm", "min_max_profile"]

# Amino acids with a single codon: no coding choice, so no commonness signal.
_NON_DEGENERATE: frozenset[str] = frozenset({"M", "W"})

# Allowed orientations for the additive term.
_DIRECTIONS: frozenset[str] = frozenset({"max", "min"})


def _synonymous_averages(frequencies: Mapping[str, float]) -> dict[str, float]:
    """Return the mean synonymous frequency per amino acid.

    Codons are grouped by the amino acid (or stop ``'*'``) they encode and each
    group's frequencies are averaged. Only codons present in ``frequencies``
    contribute, so a caller may pass any subset that covers the codons it will
    score.

    Args:
        frequencies: Mapping ``codon -> frequency`` on any consistent scale.

    Returns:
        Mapping ``amino_acid -> mean synonymous frequency``.

    Raises:
        KeyError: If a key of ``frequencies`` is not a valid DNA codon.
    """
    groups: dict[str, list[float]] = {}
    for codon, value in frequencies.items():
        aa = CODON_TABLE[codon.upper()]
        groups.setdefault(aa, []).append(value)
    return {aa: sum(values) / len(values) for aa, values in groups.items()}


def _window_min_max(xactual: float, xavg: float, xmax: float, xmin: float) -> float:
    """Return the %MinMax value of one window from its four accumulated sums.

    Implements the Clarke & Clark piecewise definition, clamped so an all
    single-codon window (where ``Xmax == Xavg == Xmin``) contributes ``0.0``
    rather than dividing by zero.

    Args:
        xactual: Sum of the actual chosen codons' frequencies.
        xavg: Sum of per-residue average synonymous frequencies.
        xmax: Sum of per-residue maximum synonymous frequencies.
        xmin: Sum of per-residue minimum synonymous frequencies.

    Returns:
        The window's %MinMax value in ``[-100.0, 100.0]``.
    """
    if xactual > xavg:
        denominator = xmax - xavg
        if denominator <= 0.0:
            return 0.0
        return 100.0 * (xactual - xavg) / denominator
    denominator = xavg - xmin
    if denominator <= 0.0:
        return 0.0
    return -100.0 * (xavg - xactual) / denominator


def min_max_profile(
    dna: str,
    frequencies: Mapping[str, float],
    window: int = 18,
) -> list[float]:
    """Compute the classic sliding-window %MinMax profile (a reporting metric).

    For each window of ``window`` consecutive codons this sums, over the codons
    in the window, the maximum (``Xmax``), minimum (``Xmin``), average
    (``Xavg``), and actual (``Xactual``) synonymous frequencies, then returns the
    Clarke & Clark %MinMax value: ``100 * (Xactual - Xavg) / (Xmax - Xavg)`` when
    ``Xactual > Xavg`` (a positive "%Max"), else ``-100 * (Xavg - Xactual) /
    (Xavg - Xmin)`` (a negative "%Min"). Every value lies in ``[-100, 100]``:
    positive means common codons, negative means rare ones.

    This is a **reporting/profile** helper only — it is deterministic and does
    not participate in the optimizer. Per-residue max/min/average are taken over
    the synonymous codons *present in* ``frequencies``. A single-codon amino acid
    (Met, Trp) contributes its one frequency equally to all four sums, so it
    cancels out of both the numerator and denominator differences and thus has no
    effect on the %MinMax value (no special-casing needed); a window consisting
    entirely of such codons has a zero denominator and is reported as ``0.0``.

    ``frequencies`` may be on any consistent scale (per-thousand, fraction, …):
    because each window is normalized by its own ``Xmax``/``Xmin``, a global
    rescale leaves the profile unchanged. Normalizing so that each amino acid's
    synonymous frequencies sum to one (a fraction within each amino acid) is the
    cleanest input.

    Args:
        dna: Coding DNA whose length is a multiple of three.
        frequencies: Mapping ``codon -> frequency``; must contain every codon of
            ``dna``.
        window: Number of consecutive codons per window (classically 17-18).

    Returns:
        The per-window %MinMax values, of length ``max(0, n_codons - window +
        1)`` (empty when the sequence is shorter than one window).

    Raises:
        ValueError: If ``window`` is less than one, if ``len(dna)`` is not a
            multiple of three, or if a codon of ``dna`` is absent from
            ``frequencies``.
    """
    if window < 1:
        raise ValueError(f"window must be >= 1, got {window!r}")
    codons = [codon.upper() for _, codon in iter_codons(dna)]
    for codon in codons:
        if codon not in frequencies:
            raise ValueError(f"codon {codon!r} missing from frequencies mapping")
    n_codons = len(codons)
    if n_codons < window:
        return []

    # Precompute each codon's (actual, avg, max, min) synonymous frequencies once.
    per_codon: list[tuple[float, float, float, float]] = []
    for codon in codons:
        aa = CODON_TABLE[codon]
        present = [c for c in synonymous_codons(aa) if c in frequencies]
        freqs = [frequencies[c] for c in present]
        per_codon.append(
            (frequencies[codon], sum(freqs) / len(freqs), max(freqs), min(freqs))
        )

    profile: list[float] = []
    for start in range(n_codons - window + 1):
        xactual = xavg = xmax = xmin = 0.0
        for actual, avg, cmax, cmin in per_codon[start : start + window]:
            xactual += actual
            xavg += avg
            xmax += cmax
            xmin += cmin
        profile.append(_window_min_max(xactual, xavg, xmax, xmin))
    return profile


@dataclass(frozen=True, slots=True)
class MinMaxTerm:
    """Additive codon-commonness prior in the %MinMax spirit (larger is better).

    Each codon's contribution is ``f(codon) - f_avg(aa)`` — how much more common
    than its synonymous average the codon is — oriented by ``direction``:
    ``"max"`` rewards common codons (positive contribution), ``"min"`` rewards
    rare codons (the negation). Because :meth:`delta` reads only the codon (never
    the growing prefix or the position), the term is ``LOCAL`` with
    ``context_len == 0`` and its per-codon deltas sum exactly to :meth:`score`
    (invariant #4, "delta == score").

    This is deliberately **not** the windowed %MinMax statistic (that is
    non-additive; see :func:`min_max_profile` and the module docstring). It is a
    cheap additive prior over codon commonness that the DP can optimize honestly.

    Holds the frequency mapping directly (``codon -> f``) rather than a codon-
    usage table, so this pure objective term depends on nothing below ``domain``.
    Build it from a table via ``MinMaxTerm(table.frequency)``. Any consistent
    frequency scale works; normalizing to a fraction within each amino acid makes
    the per-codon contributions comparable across residues.

    Attributes:
        frequencies: Mapping ``codon -> frequency`` on any consistent scale.
        direction: ``"max"`` to reward common codons (sign ``+1``) or ``"min"``
            to reward rare codons (sign ``-1``).
    """

    frequencies: Mapping[str, float]
    direction: str = "max"
    _avg: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate ``direction`` and precompute per-residue synonymous averages.

        Raises:
            ValueError: If ``direction`` is not ``"max"`` or ``"min"``.
        """
        if self.direction not in _DIRECTIONS:
            raise ValueError(
                f"direction must be one of {sorted(_DIRECTIONS)}, got {self.direction!r}"
            )
        object.__setattr__(self, "_avg", _synonymous_averages(self.frequencies))

    @property
    def name(self) -> str:
        """Return the stable identifier, e.g. ``"minmax_max"``."""
        return f"minmax_{self.direction}"

    @property
    def sign(self) -> float:
        """Return ``+1.0`` for ``"max"`` (common) or ``-1.0`` for ``"min"`` (rare)."""
        return 1.0 if self.direction == "max" else -1.0

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``0`` - each codon is scored independently."""
        return 0

    def delta(self, prefix: str, codon: str, pos: int) -> float:
        """Return the codon-commonness contribution of ``codon``.

        The value is ``sign * (f(codon) - f_avg(aa(codon)))``: for ``"max"`` a
        codon more common than its synonymous average scores positive, for
        ``"min"`` the sign is flipped. Met, Trp, and stop codons carry no coding
        choice and score ``0.0``.

        Args:
            prefix: Unused (this term needs no context).
            codon: The 3-nt codon being placed.
            pos: Unused 0-based codon index.
        """
        up = codon.upper()
        aa = CODON_TABLE[up]
        if aa == STOP or aa in _NON_DEGENERATE:
            return 0.0
        return self.sign * (self.frequencies[up] - self._avg[aa])

    def score(self, dna: str) -> float:
        """Return the sum of :meth:`delta` over the codons of ``dna``.

        Equals the running sum of :meth:`delta` the DP accumulates (CLAUDE.md
        invariant #4), since the delta depends on nothing but the codon.

        Args:
            dna: The coding sequence to score.
        """
        return sum(self.delta("", codon, pos) for pos, codon in iter_codons(dna))
