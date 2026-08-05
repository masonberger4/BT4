"""Strong splice-consensus donor/acceptor motif constraint (CLAUDE.md §6, §9 Phase 3).

A designed coding sequence can accidentally spell out a **cryptic splice site** --
a donor (5' splice site) or acceptor (3' splice site) consensus that a spliceosome
might act on, excising part of the intended ORF. This is a genuine mRNA-design
hazard. :class:`SpliceSiteMotifConstraint` is a **LOCAL, structural heuristic**
that forbids the *strong consensus* motifs of the two splice sites, on the mRNA
**sense strand only** (splicing is strand-specific -- unlike a double-stranded
restriction site, there is no reverse-complement banning here).

**What this is, and pointedly what it is not.** This is a bounded-context motif
veto, *not* a splice-prediction model. It reduces only the most *obvious* risk --
a near-perfect textbook consensus. The real cryptic-splice adjudication is done
out-of-loop by the wrapped SpliceAI / Pangolin CNNs (CLAUDE.md §6; the
:mod:`bt4.biomodels.splice` backends), which read context far beyond any fixed
motif. So a sequence that passes this constraint is *not* certified splice-safe;
it is merely free of the strong consensus motifs, and its status is honestly a
heuristic prior, never a calibrated risk claim (§10.6). It is off by default.

**Why only the strong consensus (and never the bare dinucleotide).** The
invariant intronic donor ``GT`` and acceptor ``AG`` dinucleotides are ubiquitous
in any coding sequence; forbidding them outright would make almost every protein
infeasible and is expressly *not* what this does. Only the longer, specific
consensus can be banned (design-of-record governing rule 3). The bundled defaults
(Shapiro & Senapathy 1987; Zhang 1998):

* **Donor (5' splice site)** -- :data:`DEFAULT_DONOR_MOTIFS` = ``("GTRAGT",)``: the
  strong intronic donor core, intron positions +1..+6 of the ``MAG|GTRAGT``
  consensus -- the invariant ``GT`` plus the highly conserved ``RAGT``. The exonic
  ``MAG`` flank is deliberately *not* pinned: a strong donor does not require a
  fixed exon-side context, so pinning ``MAGGTRAGT`` (~1/65000) would *under*-forbid,
  missing genuine strong donors whose -3..-1 flank diverges. ``GTRAGT`` (~1/2048 in
  i.i.d. sequence) is the shortest pattern that is still unambiguously *strong donor
  consensus* rather than the bare ``GT``, and it leaves abundant synonymous escape
  codons.
* **Acceptor (3' splice site)** -- :data:`DEFAULT_ACCEPTOR_MOTIFS` =
  ``("YYYYYYNYAGG",)``: a 6-nt pyrimidine tract, the branch-proximal ``N``, then the
  conserved ``YAG`` ending the intron (``AG`` invariant at -2/-1, with ``Y`` at -3
  capturing both ``CAG`` and ``TAG`` acceptors) and the first exonic ``G`` (+1):
  ``YYYYYY N YAG|G``. The polypyrimidine run is what distinguishes a strong acceptor
  from the ubiquitous bare ``AG``. A 6-Y tract (~1/8192) is the defensible middle:
  longer (e.g. 10-Y) is so rare it essentially never fires, shorter over-fires on
  incidentally pyrimidine-rich coding windows.

Because both defaults are specific, this constraint fires only on strong
consensus; that conservatism is intentional and matches its role (obvious-risk
filter, with the CNN audit doing the real work). A user who wants a stricter or
custom filter may pass their own IUPAC ``donor_motifs`` / ``acceptor_motifs``.

**Honest scope limits (documented, not hidden).** (1) This is the human /
mammalian **major-spliceosome** (U2, ``GT-AG``) consensus only -- it does *not*
cover U12/minor-spliceosome (``AT-AC``) sites, non-canonical ``GC-AG`` donors, or
non-mammalian consensus, so it should be narrowed or disabled for other organisms
rather than assumed transferable. (2) A real 3' splice site is a *variable-length,
often-interrupted* polypyrimidine tract plus an upstream branch point plus the
terminal ``NYAG``; a fixed-window motif cannot capture that, so the acceptor arm
deliberately under-fires on long/interrupted tracts and reduces only *obvious*
acceptor consensus. (3) Collapsing the source position-weight matrices to a single
IUPAC string loses per-position information -- the degeneracy choices are heuristic
specificity/feasibility knobs, which is exactly why they are exposed as parameters.

**Feasibility.** As a hard veto in the exact DP, a pathological protein could in
principle be made infeasible if a forbidden core is pinned by synonymously
immovable codons; with these specific defaults that is rare (a designer almost
always has a synonymous codon avoiding a stray core). Like every other hard
constraint here it makes no ``relax()`` claim (that is not in the shipped
``Constraint`` protocol); an infeasible instance surfaces honestly as
``InfeasibleError``, and the user opts out by disabling this off-by-default
constraint or narrowing the patterns (the acceptor being the harder to satisfy).

**Why ``context_len`` suffices (invariant #3).** Exactly as
:class:`~bt4.constraints.rules.ForbiddenMotifConstraint`, an occurrence of an
``L``-length pattern can straddle the codon seam; catching one whose rightmost base
falls inside the incoming codon needs the ``L - 1`` preceding bases from the
prefix. With ``M = max`` pattern length, ``context_len == M - 1`` and the veto
inspects ``prefix[-(M-1):] + next_codon``, flagging only occurrences that overlap
``next_codon`` (a match fully inside the already-feasible prefix is not this
extension's fault). A sequence built respecting ``ok_suffix`` therefore carries
zero hard violations under ``validate`` -- property-tested.

This constraint is purely hard (``penalty`` is ``0.0``).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.constraints.iupac import find_iupac, is_iupac
from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = [
    "DEFAULT_ACCEPTOR_MOTIFS",
    "DEFAULT_DONOR_MOTIFS",
    "SpliceSiteMotifConstraint",
]

# Strong 5' splice-site (donor) consensus core, sense strand: intron +1..+6 of the
# MAG|GTRAGT consensus -- invariant GT plus the conserved RAGT (Shapiro & Senapathy
# 1987; Zhang 1998). The exon-side MAG is intentionally not pinned (it varies among
# strong donors). Never the bare intronic GT.
DEFAULT_DONOR_MOTIFS: tuple[str, ...] = ("GTRAGT",)

# Strong 3' splice-site (acceptor) consensus, sense strand: a 6-nt polypyrimidine
# tract + N + the conserved YAG ending the intron + the first exonic G (YAG | G).
# Y at -3 captures both CAG and TAG acceptors. Never the bare acceptor AG (the
# pyrimidine run is what makes it specific).
DEFAULT_ACCEPTOR_MOTIFS: tuple[str, ...] = ("YYYYYYNYAGG",)


@dataclass(frozen=True, slots=True)
class SpliceSiteMotifConstraint:
    """Forbid strong splice-consensus donor/acceptor motifs (sense strand, LOCAL).

    An honest structural heuristic: it bans the *strong* 5'ss / 3'ss consensus
    (IUPAC patterns), never the ubiquitous bare ``GT`` / ``AG``. It makes no
    calibrated splice-risk claim; the wrapped SpliceAI/Pangolin CNNs do the real
    audit out of loop (CLAUDE.md §6, §10.6).

    Attributes:
        donor_motifs: IUPAC patterns for the 5' splice-site (donor) consensus,
            sense strand. Defaults to :data:`DEFAULT_DONOR_MOTIFS`.
        acceptor_motifs: IUPAC patterns for the 3' splice-site (acceptor)
            consensus, sense strand. Defaults to :data:`DEFAULT_ACCEPTOR_MOTIFS`.
            Pass an empty tuple for either side to disable that arm; both empty
            makes the constraint an inert no-op.
    """

    donor_motifs: tuple[str, ...] = DEFAULT_DONOR_MOTIFS
    acceptor_motifs: tuple[str, ...] = DEFAULT_ACCEPTOR_MOTIFS
    name: str = field(default="splice_site", init=False)
    # (pattern, kind) pairs, deduplicated and sorted for deterministic auditing.
    _patterns: tuple[tuple[str, str], ...] = field(init=False, repr=False, compare=False)
    _maxlen: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate and normalize the IUPAC patterns.

        Raises:
            ValueError: If any motif is not a valid non-empty IUPAC string.
        """
        collected: dict[str, str] = {}
        for kind, motifs in (("donor", self.donor_motifs), ("acceptor", self.acceptor_motifs)):
            for motif in motifs:
                up = motif.upper()
                if not is_iupac(up):
                    raise ValueError(f"not a valid IUPAC {kind} motif: {motif!r}")
                # First kind to claim a pattern owns its label (donor before acceptor);
                # identical patterns are not double-counted.
                collected.setdefault(up, kind)
        patterns = tuple(sorted(collected.items()))
        maxlen = max((len(p) for p, _ in patterns), default=0)
        object.__setattr__(self, "_patterns", patterns)
        object.__setattr__(self, "_maxlen", maxlen)

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.LOCAL`."""
        return Scope.LOCAL

    def context_len(self) -> int:
        """Return ``max(0, maxlen - 1)`` -- enough for a motif to cross the seam."""
        return max(0, self._maxlen - 1)

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff a consensus motif occurrence overlaps ``next_codon``.

        The window is ``prefix[-(maxlen-1):] + next_codon``; an occurrence counts
        only when its exclusive end index lies beyond the start of ``next_codon``
        (so it actually involves the codon being appended). With no motifs the
        constraint is inert and this always returns ``True``.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        if not self._patterns:
            return True
        tail = self._maxlen - 1
        window = (prefix[-tail:] if tail > 0 else "") + next_codon
        seam = len(window) - len(next_codon)
        for pattern, _kind in self._patterns:
            for idx in find_iupac(window, pattern):
                if idx + len(pattern) > seam:
                    return False
        return True

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per consensus motif occurrence in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` for every occurrence of every
            donor/acceptor consensus motif, in a deterministic (pattern-sorted)
            order, each tagged with which splice site it matches.
        """
        seq = dna.upper()
        for pattern, kind in self._patterns:
            for idx in find_iupac(seq, pattern):
                yield Violation(
                    constraint="splice_site",
                    severity=Severity.HARD,
                    start=idx,
                    end=idx + len(pattern),
                    detail=f"strong {kind} splice-consensus motif {pattern} at {idx}",
                )
