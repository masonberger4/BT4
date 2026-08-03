"""uORF constraint: an out-of-frame internal ATG paired with a downstream stop.

BT3's ``almost-there`` branch shipped an ``out_of_frame_atg`` rule that flagged
an out-of-frame internal ATG which pairs with a downstream in-frame stop -- i.e.
a short **upstream/internal open reading frame (uORF)** embedded in the coding
sequence. uORFs are among the strongest *cis* dampeners of translation: a
ribosome that leaky-scans past the real start and initiates at the uORF's ATG is
diverted from the main ORF, lowering yield and producing a truncated peptide
(CLAUDE.md §6). Synonymous codon shuffling readily spawns such out-of-frame ATGs,
so a design tool should be able to suppress the ones it accidentally created.

This is the genuinely **non-local** half of internal-ATG handling, and it is kept
honestly separate from two neighbours:

* :class:`~bt4.constraints.kozak.InternalStartConstraint` (``avoid_internal_start``)
  is LOCAL: it bans an internal ATG in a *strong Kozak context* (purine at -3, G
  at +4), in any frame, within a bounded window. It says nothing about pairing.
* This constraint is ``Scope.GLOBAL``: a uORF's ATG and its closing in-frame stop
  can sit an arbitrary distance apart, so no bounded ``ok_suffix`` window sees
  both. It is therefore **never merged into the exact-DP trellis** (that would
  silently over-merge, CLAUDE.md §10.1); it is enforced in the refinement layer
  exactly like :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint`, and any
  residual uORFs are reported honestly.

What counts as a uORF here. For an ATG whose A-index ``a`` satisfies:

* ``a >= min_start`` (internal -- the real start codon at index 0 is skipped);
* ``a % 3 != 0`` (**out of the main reading frame** -- an in-frame internal ATG
  is just a Met codon forced by the protein and is not removable anyway);
* ``a < region_nt`` when a 5'-proximal scan window is set (uORFs matter most near
  the start, where leaky scanning happens -- and an unbounded ban would be
  infeasible, since out-of-frame ATG..stop pairs are common in any long CDS);

and for which the **first stop codon in the ATG's own frame** exists downstream,
the span ``ATG .. stop`` is a uORF. This is a purely structural, deterministic
rule: it makes **no** claim about how much expression the uORF costs (that would
require the validated expression model of Phase 4) -- it only suppresses a
known-repressive motif.

Honesty of ``ok_suffix`` <=> ``validate`` (invariant #3). Because a uORF is
closed by its *downstream* stop, appending a codon can only *complete* a uORF by
supplying that stop; the ATG is always already in the prefix. So ``ok_suffix`` is
defined directly in terms of :meth:`~UorfConstraint.validate`: it vetoes the
extension iff a uORF appears whose closing stop lands in the incoming codon. A
build respecting ``ok_suffix`` therefore never completes a uORF, so ``validate``
finds none -- the two agree by construction.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from dataclasses import dataclass, field

from bt4.domain.result import Severity, Violation
from bt4.domain.scope import Scope

__all__ = ["UorfConstraint"]

_STOPS = frozenset({"TAA", "TAG", "TGA"})


@dataclass(frozen=True, slots=True)
class UorfConstraint:
    """Ban an out-of-frame internal ATG that pairs with a downstream in-frame stop.

    A violation is any out-of-main-frame internal ATG (A-index ``a`` with
    ``a >= min_start``, ``a % 3 != 0``, and ``a < region_nt`` when a window is
    set) for which the first stop codon in the ATG's own frame exists downstream
    -- a bounded uORF ``ATG .. stop``. The ATG and stop can be far apart, so this
    is ``Scope.GLOBAL`` and is enforced in the refinement layer, never the trellis.

    Attributes:
        min_start: Only ATGs whose A-index is ``>= min_start`` are "internal". The
            default of ``3`` skips the first codon (the real start).
        region_nt: If not ``None``, only ATGs whose A-index is ``< region_nt`` are
            scanned -- the 5'-proximal window where leaky-scanning re-initiation
            matters (BT3 used 200; the default here is 100 = the first ~33 codons).
            ``None`` scans the whole sequence (rarely feasible for a long CDS).
    """

    min_start: int = 3
    region_nt: int | None = 100
    name: str = field(default="uorf", init=False)

    def __post_init__(self) -> None:
        """Validate the bounds.

        Raises:
            ValueError: If ``min_start`` is negative or ``region_nt`` is set and
                not positive.
        """
        if self.min_start < 0:
            raise ValueError(f"min_start must be >= 0, got {self.min_start!r}")
        if self.region_nt is not None and self.region_nt <= 0:
            raise ValueError(f"region_nt must be > 0 or None, got {self.region_nt!r}")

    def scope(self) -> Scope:
        """Return :attr:`~bt4.domain.scope.Scope.GLOBAL`.

        A uORF's ATG and its closing in-frame stop can lie any distance apart, so
        the veto depends on the whole sequence. The pipeline enforces it via
        :meth:`validate` in the refinement layer, never in the merged trellis
        (CLAUDE.md §10.1).
        """
        return Scope.GLOBAL

    def context_len(self) -> int:
        """Return ``sys.maxsize`` -- an unbounded sentinel.

        :meth:`ok_suffix` inspects the whole prefix (the paired ATG can be
        anywhere upstream), so no finite trailing window suffices; the sentinel
        states that honestly (mirroring
        :class:`~bt4.constraints.max_repeat.MaxRepeatConstraint`). It is safe: were
        this constraint ever handed to the trellis, a whole-sequence merge key
        would prevent *any* state merge (correct but exponential), never a silent
        over-merge -- which is why the pipeline routes it to refinement instead.
        """
        return sys.maxsize

    def _uorf_spans(self, seq: str) -> list[tuple[int, int]]:
        """Return the sorted ``(atg_start, uorf_end)`` spans of every uORF in ``seq``.

        For each qualifying out-of-frame internal ATG, finds the first stop codon
        in the ATG's own frame downstream; when one exists the pair is a uORF and
        its span ``(a, stop_index + 3)`` is emitted.

        Args:
            seq: The (upper-cased) sequence to scan.
        """
        n = len(seq)
        spans: list[tuple[int, int]] = []
        for a in range(self.min_start, n - 2):
            if self.region_nt is not None and a >= self.region_nt:
                break  # ATGs only open within the 5'-proximal window; a only grows.
            if a % 3 == 0:
                continue  # in the main frame -> a Met codon, not an out-of-frame uORF.
            if seq[a : a + 3] != "ATG":
                continue
            s = a + 3  # same frame as a, since (a + 3) % 3 == a % 3.
            while s + 3 <= n:
                if seq[s : s + 3] in _STOPS:
                    spans.append((a, s + 3))
                    break
                s += 3
        return spans

    def ok_suffix(self, prefix: str, next_codon: str) -> bool:
        """Return ``False`` iff appending ``next_codon`` completes a new uORF.

        Reads the whole prefix (this is a global constraint). A uORF is closed by
        its downstream stop, so a new one appears only when the incoming codon
        supplies that stop; such a uORF has ``uorf_end > len(prefix)``. ``prefix``
        is assumed already feasible.

        Args:
            prefix: The feasible DNA chosen so far.
            next_codon: The codon being appended.
        """
        full = (prefix + next_codon).upper()
        p = len(prefix)
        return not any(end > p for _a, end in self._uorf_spans(full))

    def penalty(self, prefix: str, next_codon: str) -> float:
        """Return ``0.0`` (this constraint is purely hard)."""
        return 0.0

    def validate(self, dna: str) -> Iterator[Violation]:
        """Yield one HARD violation per uORF in ``dna``.

        Args:
            dna: The whole coding sequence to audit.

        Yields:
            A :class:`~bt4.domain.result.Violation` spanning each ``ATG .. stop``
            uORF (out-of-frame internal ATG with its first in-frame stop), in
            ascending ATG order.
        """
        seq = dna.upper()
        for a, end in self._uorf_spans(seq):
            yield Violation(
                constraint="uorf",
                severity=Severity.HARD,
                start=a,
                end=end,
                detail=(
                    f"out-of-frame uORF: internal ATG at {a} (frame +{a % 3}) "
                    f"pairs with an in-frame stop ending at {end}"
                ),
            )
