"""Cross-backend splice agreement -- a first-class uncertainty signal.

CLAUDE.md section 6 makes backend *agreement* a design goal, not redundancy:
"Running *both* [SpliceAI and Pangolin] and reporting their agreement /
disagreement is a first-class uncertainty signal (section 8), not redundancy."
This module computes that signal for any set of
:class:`~bt4.biomodels.splice.base.SplicePredictor` backends -- the labeled PWM
baseline, Pangolin, and (later) SpliceAI -- over a panel of candidate sequences
scored against one reference.

The comparison is done at the **backend-agnostic** level -- pooled
Delta-splicing per candidate -- because backends differ in what per-position
scores mean (Pangolin reports one combined splice-site probability; the PWM
baseline reports separated donor / acceptor). Delta-splicing (larger is better)
and its ranking across candidates are directly comparable across backends, so
this module reports:

* each backend's Delta-splicing vector over the candidates;
* pairwise **Spearman rank correlation** of those vectors -- do the backends
  agree on *which* redesigns lower splice risk?
* the fraction of candidates on which all backends agree on the *sign* of
  Delta-splicing (did the redesign help or hurt?).

It **reports, it does not judge**: no backend is declared right, and a low
correlation is surfaced as uncertainty, not resolved. Pure standard library (no
numpy); it depends only on the :class:`SplicePredictor` protocol and
:mod:`bt4.domain` (transitively via the protocol).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt4.biomodels._stats import pearson, spearman
from bt4.biomodels.splice.base import DEFAULT_TOP_K, SplicePredictor, pooled_risk

__all__ = [
    "AgreementReport",
    "SiteCallAgreement",
    "agreement_from_deltas",
    "backend_agreement",
    "pearson",
    "site_call_agreement",
    "spearman",
]


@dataclass(frozen=True, slots=True)
class AgreementReport:
    """Cross-backend splice-risk agreement over a candidate panel.

    Attributes:
        backends: The backend names, in input order.
        delta_by_backend: For each backend, its Delta-splicing (larger is better)
            for each candidate vs the reference, aligned to the candidate order.
        rank_correlations: Pairwise Spearman rank correlation of the
            Delta-splicing vectors, keyed by ``(name_a, name_b)`` for each
            unordered backend pair. Empty when fewer than two candidates (rank
            correlation is undefined) or fewer than two backends.
        sign_agreement: Fraction of candidates on which *all* backends agree on
            the sign of Delta-splicing (treating a value within ``sign_epsilon``
            of zero as "no change"). ``1.0`` when there is only one backend.
        sign_epsilon: The magnitude below which a Delta-splicing is treated as
            zero for the sign-agreement count.
        n_candidates: Number of candidates scored.
    """

    backends: tuple[str, ...]
    delta_by_backend: dict[str, tuple[float, ...]]
    rank_correlations: dict[tuple[str, str], float]
    sign_agreement: float
    sign_epsilon: float
    n_candidates: int


def _sign(value: float, epsilon: float) -> int:
    """Return -1 / 0 / +1 for ``value``, with a dead-band of ``epsilon`` at zero."""
    if value > epsilon:
        return 1
    if value < -epsilon:
        return -1
    return 0


def agreement_from_deltas(
    delta_by_backend: dict[str, tuple[float, ...]],
    *,
    sign_epsilon: float = 1e-9,
) -> AgreementReport:
    """Build an :class:`AgreementReport` from **precomputed** Delta-splicing vectors.

    The rank/sign half of :func:`backend_agreement`, split out so a caller that has
    already scored every candidate (e.g.
    :func:`bt4.biomodels.splice.audit.audit_splice`) can reuse those Delta-splicing
    values instead of re-running the backends over the whole panel a second time
    (CLAUDE.md §7). Each vector must be aligned to the same candidate order and have
    the same length.

    Args:
        delta_by_backend: ``{backend name: Delta-splicing per candidate}`` (larger =
            better), in the desired backend order (dict insertion order is kept).
        sign_epsilon: Dead-band around zero for sign agreement.

    Returns:
        An :class:`AgreementReport` over the given deltas.
    """
    names = tuple(delta_by_backend)
    n_candidates = len(next(iter(delta_by_backend.values()))) if delta_by_backend else 0

    rank_correlations: dict[tuple[str, str], float] = {}
    if len(names) >= 2 and n_candidates >= 2:
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                rank_correlations[(names[i], names[j])] = spearman(
                    delta_by_backend[names[i]], delta_by_backend[names[j]]
                )

    if len(names) <= 1 or n_candidates == 0:
        sign_agreement = 1.0
    else:
        agree = sum(
            1
            for idx in range(n_candidates)
            if len({_sign(delta_by_backend[name][idx], sign_epsilon) for name in names}) == 1
        )
        sign_agreement = agree / n_candidates

    return AgreementReport(
        backends=names,
        delta_by_backend=dict(delta_by_backend),
        rank_correlations=rank_correlations,
        sign_agreement=sign_agreement,
        sign_epsilon=sign_epsilon,
        n_candidates=n_candidates,
    )


def backend_agreement(
    predictors: Sequence[SplicePredictor],
    candidates: Sequence[str],
    reference: str,
    *,
    top_k: int | None = None,
    sign_epsilon: float = 1e-9,
) -> AgreementReport:
    """Compare splice backends over ``candidates`` scored against ``reference``.

    For each backend and each candidate, computes Delta-splicing (larger is
    better) as the contract's ``pooled_risk(reference) - pooled_risk(candidate)``
    -- equal to :meth:`SplicePredictor.delta_splicing`, but scoring the reference
    once per backend rather than per candidate -- then reports pairwise Spearman
    rank agreement and sign agreement across backends.

    Args:
        predictors: The splice backends to compare. Names must be distinct.
        candidates: Candidate coding sequences (e.g. frontier members) to rank.
        reference: The reference sequence each candidate is compared against.
        top_k: Optional pooling depth override. When ``None`` each backend is
            pooled at its own depth (its ``top_k`` attribute, else
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`); when set, every
            backend is pooled at this ``top_k`` uniformly. Either way the reference
            is scored once per backend and reused across candidates.
        sign_epsilon: Dead-band around zero for sign agreement.

    Returns:
        An :class:`AgreementReport`.

    Raises:
        ValueError: If ``predictors`` or ``candidates`` is empty, or two
            predictors share a name.
    """
    if not predictors:
        raise ValueError("need at least one splice backend to compare")
    if not candidates:
        raise ValueError("need at least one candidate sequence to compare")
    names = tuple(p.name for p in predictors)
    if len(set(names)) != len(names):
        raise ValueError(f"backend names must be distinct, got {names}")

    delta_by_backend: dict[str, tuple[float, ...]] = {}
    for predictor in predictors:
        # Pool at `top_k` when overridden, else at the backend's own depth (its
        # `top_k` attribute, or DEFAULT_TOP_K if it exposes none). The contract
        # defines delta_splicing as pooled_risk(reference) - pooled_risk(designed),
        # so scoring the reference ONCE per backend and reusing it equals calling
        # delta_splicing per candidate -- while avoiding re-running the reference
        # through a heavy CNN C times (CLAUDE.md section 7, "everything incremental").
        depth = top_k if top_k is not None else getattr(predictor, "top_k", DEFAULT_TOP_K)
        ref_risk = pooled_risk(predictor.score_sequence(reference), depth)
        deltas = [
            ref_risk - pooled_risk(predictor.score_sequence(cand), depth)
            for cand in candidates
        ]
        delta_by_backend[predictor.name] = tuple(deltas)

    return agreement_from_deltas(delta_by_backend, sign_epsilon=sign_epsilon)


# --------------------------------------------------------------------------
# Site-panel agreement: do two backends call the SAME positions?


@dataclass(frozen=True, slots=True)
class SiteCallAgreement:
    """Whether two backends flag the same positions on a site-prediction panel.

    A different question from :class:`AgreementReport`, which ranks *candidate
    sequences* by Delta-splicing for the design flow. Here the panel is annotated
    genomic sequence and the question is positional: given the same window, do the
    two models point at the same bases?

    It is not answered by comparing their gate metrics. Two backends can each score
    skill 0.98 on the same panel while being confident about different positions --
    the metrics would look like agreement and the models would not agree.

    Attributes:
        backends: The two backend names, in input order.
        n_sites: Annotated sites across the panel -- also the per-window ``k``.
        both / only_first / only_second / neither: How many annotated sites each
            backend recovered in its own top-``k``. ``neither`` is a site both
            missed, which is a real and reportable outcome.
        jaccard: ``|A n B| / |A u B|`` over the *called* positions (each backend's
            top-``k`` per window), whether or not those calls are correct. This is
            the headline: it measures agreement about where the sites are, and is
            unaffected by prevalence because ``k`` is fixed by the annotation.
        spearman_on_called: Rank correlation of the two backends' scores, computed
            over the **union of their called positions** rather than every base.
            Over a whole panel the correlation is dominated by the ~99.96% of
            positions where both are ~0 and would read near 1.0 regardless of
            whether the models agree anywhere that matters.
        n_called_union: Size of that union, so the correlation has a denominator.
    """

    backends: tuple[str, str]
    n_sites: int
    both: int
    only_first: int
    only_second: int
    neither: int
    jaccard: float
    spearman_on_called: float
    n_called_union: int


def _top_k_indices(scores: Sequence[float], k: int) -> set[int]:
    """Return the indices of the ``k`` highest **positive** scores.

    A call requires a positive score. Taking the top ``k`` unconditionally pads the
    set with zero-scored positions whenever a backend has fewer than ``k`` peaks,
    and those padded entries are chosen by the index tie-break -- so two backends
    that found nothing in common still "agree" on positions 0, 1, 2 ... The error
    is in the flattering direction, which is the one to design against: a window
    where both models are silent would report partial agreement rather than none.

    Ties among genuinely positive scores are still broken by index, which keeps the
    result deterministic (invariant #7).
    """
    if k <= 0:
        return set()
    positive = [i for i, score in enumerate(scores) if score > 0.0]
    positive.sort(key=lambda i: (-scores[i], i))
    return set(positive[:k])


def site_call_agreement(
    first_tracks: Sequence[Sequence[float]],
    second_tracks: Sequence[Sequence[float]],
    site_indices: Sequence[Sequence[int]],
    names: tuple[str, str],
) -> SiteCallAgreement:
    """Compare where two backends call sites, window by window.

    Each backend's calls are its top-``k`` positions **within each window**, where
    ``k`` is that window's annotated-site count. Per-window rather than pooled: a
    long window would otherwise absorb every call and a short one none, which
    measures window length rather than model agreement.

    Args:
        first_tracks / second_tracks: One score track per window, already collapsed
            to a single per-position series (see ``pipeline.splice_gate._tracks``)
            and in the **same frame** as each other.
        site_indices: The annotated site positions per window, in the backends'
            frame -- i.e. already anchor-shifted. Both wrapped CNNs share the same
            anchors, so backend-vs-backend agreement needs no shift between them;
            the shift is only needed to compare either against the annotation.
        names: The two backend names.

    Returns:
        A :class:`SiteCallAgreement`.

    Raises:
        ValueError: On mismatched window counts, or two identical names (which
            would make the report unreadable rather than wrong).
    """
    if names[0] == names[1]:
        raise ValueError(f"backend names must differ, got {names[0]!r} twice")
    if not (len(first_tracks) == len(second_tracks) == len(site_indices)):
        raise ValueError(
            f"window counts differ: {len(first_tracks)} / {len(second_tracks)} tracks "
            f"for {len(site_indices)} annotated windows"
        )

    both = only_first = only_second = neither = 0
    intersection = union = 0
    paired_first: list[float] = []
    paired_second: list[float] = []

    for first, second, sites in zip(first_tracks, second_tracks, site_indices, strict=True):
        if len(first) != len(second):
            raise ValueError(
                f"track lengths differ within a window: {len(first)} vs {len(second)}"
            )
        k = len(set(sites))
        called_first = _top_k_indices(first, k)
        called_second = _top_k_indices(second, k)

        intersection += len(called_first & called_second)
        called_union = called_first | called_second
        union += len(called_union)
        for index in sorted(called_union):
            paired_first.append(first[index])
            paired_second.append(second[index])

        for site in set(sites):
            in_first = site in called_first
            in_second = site in called_second
            if in_first and in_second:
                both += 1
            elif in_first:
                only_first += 1
            elif in_second:
                only_second += 1
            else:
                neither += 1

    return SiteCallAgreement(
        backends=names,
        n_sites=both + only_first + only_second + neither,
        both=both,
        only_first=only_first,
        only_second=only_second,
        neither=neither,
        # An empty union means neither backend called anything, which is total
        # agreement about nothing -- 0.0 is the honest reading, not 1.0.
        jaccard=(intersection / union) if union else 0.0,
        # `spearman` raises below two points; a union that small carries no rank
        # information, so 0.0 is the honest report and matches `_stats`'s own
        # convention of returning 0.0 on degenerate input rather than raising.
        spearman_on_called=(
            spearman(paired_first, paired_second) if len(paired_first) >= 2 else 0.0
        ),
        n_called_union=union,
    )
