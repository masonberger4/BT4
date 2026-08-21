#!/usr/bin/env python3
"""Zero-data sanity checks on the wrapped RiboNN head, before any panel is bought.

RiboNN ships ``calibrated=False`` because it was trained and validated across
*different natural genes*, while BT4 asks a narrower question: given one protein and
one fixed UTR pair, which **synonymous** coding sequence is best? That experiment has
never been published (see ``docs/RESEARCH_codon_optimization_SOTA.md``), so promotion
needs a measured CDS-variant panel and a passing
:func:`bt4.biomodels.expression.verify_expression_gate`.

This script is what you run **first**, and it needs no measured data at all -- only a
RiboNN checkout and sequences already in this repository. Its job is to find out
whether a panel is even worth acquiring, and several of its outcomes are decisive on
their own:

``utr-control``
    Positive control. Score one CDS under two *different* UTR pairs; the scores must
    move. If they do not, the harness is broken and every "no effect" result below
    would be a wiring bug misread as biology. **Run this before believing anything.**

``cds-spread``
    The decisive check, and the confound check with it (one RiboNN invocation covers
    both). Holding the UTRs fixed, how much does the score move across *synonymous
    variants of the same protein*, compared with how much it moves *between*
    proteins? And is whatever movement exists just GC3 or CAI wearing a hat? A
    within-protein spread at the noise floor means RiboNN cannot rank BT4's
    candidates, no panel will change that, and the honest outcome is to record it.

``direction``
    Crude smoke test: over many proteins, does RiboNN prefer a max-CAI design to a
    min-CAI one more often than chance? Reported as an exact two-sided sign test. A
    coin flip is a strong warning; a clear preference proves nothing on its own,
    because optimized sequences also differ in GC.

``ladder``
    Dose-response along a real BT4 frontier for one protein. A jagged, noise-like
    response says the signal is not usable for ranking even if it is nonzero.

**Nothing here is a calibration gate and nothing here can promote a backend.** Every
number is in RiboNN's native CLR-residual TE units, carries ``calibrated=False``, and
is a diagnostic. A verdict of "responds to synonymous change" is a licence to go
looking for a panel -- not evidence that the response is *correct*.

Run it directly::

    # Real UTRs are required and are YOUR modelling choice -- they are recorded in
    # the report rather than defaulted, because a bundled UTR would be a hidden one.
    python scripts/ribonn_sensitivity.py --check utr-control \
        --utr5 utr5.fa --utr3 utr3.fa --utr5-alt alt5.fa --utr3-alt alt3.fa

    python scripts/ribonn_sensitivity.py --check cds-spread \
        --fasta scripts/data/ranaghan2021_tab4.fasta \
        --utr5 utr5.fa --utr3 utr3.fa --json > stage1_spread.json

On Windows use ``^`` for line continuations, and note that the adapter already
defaults to ``num_workers=0`` (required there) and ``batch_size=64``.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import statistics
from collections.abc import Mapping, Sequence
from typing import Any

from bt4 import api
from bt4._accel import gc_count
from bt4.biomodels._stats import spearman
from bt4.biomodels.codon.tables import CodonUsageTable, load_table
from bt4.biomodels.codon.tai import TaiTable, load_tai_table
from bt4.biomodels.expression import (
    BatchExpressionPredictor,
    ExpressionPredictor,
    resolve_backend,
)
from bt4.domain.genetic_code import STOP, translate
from bt4.io.fasta import read_fasta

__all__ = [
    "CHECKS",
    "build_report",
    "group_of",
    "main",
    "score_all",
    "sign_test_p_value",
]

CHECKS = ("utr-control", "cds-spread", "direction", "ladder")

# The in-tree CC BY 4.0 panel: 3 human proteins x (1 native + 3 algorithms x 10 runs).
# Sequence-only -- Ranaghan et al. measured expression for exactly one sequence, in
# E. coli, so this panel is usable for sensitivity work and NOT as a validation panel.
DEFAULT_PANEL = (
    pathlib.Path(__file__).resolve().parent / "data" / "ranaghan2021_tab4.fasta"
)

# Below this, two scores are the same number as far as float arithmetic goes; used
# only to say "did anything move at all", never as a biological threshold.
_NOISE_FLOOR = 1e-9


# --------------------------------------------------------------------------- #
# helpers


def group_of(header: str) -> str:
    """Return the grouping key for a FASTA header: the text before the first ``|``.

    The in-tree panel uses ``Protein|Source|runN|...``, so this yields the protein --
    which is the right grouping unit throughout: synonymous variants of one protein
    are a dependent cluster, and the whole point of ``cds-spread`` is to separate
    *within*-protein movement from *between*-protein movement. A header with no ``|``
    is its own group.
    """
    return header.split("|", 1)[0].strip() or header.strip()


def read_flank(value: str) -> str:
    """Return an ACGT flank from a literal sequence or a FASTA path.

    Mirrors the CLI's flank handling: a path is read and its records concatenated, a
    bare string is taken literally. Whitespace is stripped and the result upper-cased.

    Raises:
        ValueError: If the result is empty or contains a non-ACGT character. RiboNN
            refuses empty UTRs (its loader reads an all-empty column as NaN), and the
            UTRs carry most of its signal, so an empty flank is never acceptable here.
    """
    path = pathlib.Path(value)
    raw = "".join(seq for _h, seq in read_fasta(path)) if path.is_file() else value
    flank = "".join(raw.split()).upper()
    if not flank:
        raise ValueError("UTR context is empty; RiboNN requires non-empty utr5/utr3")
    bad = sorted({ch for ch in flank if ch not in "ACGT"})
    if bad:
        raise ValueError(f"UTR context has non-ACGT character(s): {bad}")
    return flank


def gc3_fraction(dna: str) -> float:
    """Return the GC fraction at codon third positions (0.0 for an empty sequence).

    Third-position GC is the compositional axis synonymous recoding moves most, so it
    is the confounder to check a "codon" signal against: a model whose within-protein
    response is explained by GC3 is a GC detector in this regime, whatever it is
    called.
    """
    thirds = dna.upper()[2::3]
    if not thirds:
        return 0.0
    return sum(1 for ch in thirds if ch in "GC") / len(thirds)


def sign_test_p_value(successes: int, trials: int) -> float:
    """Return the exact two-sided binomial p-value against ``p = 0.5``.

    Computed with :func:`math.comb` -- no scipy, matching BT4's dependency-free
    statistics (:mod:`bt4.biomodels._stats`). Two-sided by doubling the smaller tail
    and clamping at 1.0, the standard exact convention for a symmetric null.
    """
    if trials <= 0:
        return 1.0
    tail = min(successes, trials - successes)
    cumulative = math.fsum(math.comb(trials, k) for k in range(tail + 1))
    return min(1.0, 2.0 * cumulative / (2.0**trials))


def score_all(predictor: ExpressionPredictor, dnas: Sequence[str]) -> list[float]:
    """Score every sequence, using one batched invocation where the backend allows.

    RiboNN's cost is dominated by fixed per-invocation overhead, so a batched
    backend must be driven through :meth:`score_many`; falling back to per-sequence
    scoring would multiply the wall clock by ``len(dnas)``.
    """
    if not dnas:
        return []
    if isinstance(predictor, BatchExpressionPredictor):
        return [result.score for result in predictor.score_many(list(dnas))]
    return [predictor.score_sequence(dna).score for dna in dnas]


def _spread(values: Sequence[float]) -> dict[str, float]:
    """Return ``n`` / ``mean`` / ``sd`` / ``range`` for a series of scores."""
    n = len(values)
    return {
        "n": n,
        "mean": statistics.fmean(values) if n else 0.0,
        "sd": statistics.stdev(values) if n > 1 else 0.0,
        "range": (max(values) - min(values)) if n else 0.0,
    }


def _protein_of(dna: str) -> str:
    """Return the stop-free protein encoded by ``dna``."""
    return translate(dna).rstrip(STOP)


# --------------------------------------------------------------------------- #
# the checks


def check_utr_control(
    predictor: ExpressionPredictor,
    cds: str,
    utrs: tuple[str, str],
    alt_utrs: tuple[str, str],
) -> dict[str, Any]:
    """Positive control: the same CDS under two UTR pairs must score differently.

    This is the check that makes a *null* result elsewhere trustworthy. RiboNN's own
    attribution puts the densest per-nucleotide signal in the 5'UTR, so if swapping
    both UTRs leaves the score untouched, the sequences are not reaching the model and
    nothing else this script reports means anything.

    Returns:
        A report with both scores, their difference, and ``harness_ok`` -- ``True``
        only when the score actually moved.
    """
    import dataclasses

    # The check is defined by swapping the UTR context, so a backend that has no UTR
    # context cannot take it. Refuse with a clear message instead of letting
    # dataclasses.replace raise about an unexpected field.
    fields = (
        {f.name for f in dataclasses.fields(predictor)}
        if dataclasses.is_dataclass(predictor)
        else set()
    )
    if not {"utr5", "utr3"} <= fields:
        raise ValueError(
            f"the utr-control check needs a UTR-aware backend; {predictor.name} has no "
            "utr5/utr3 context to swap (the null placeholder cannot be controlled this "
            "way -- use --backend ribonn)"
        )
    alt = dataclasses.replace(predictor, utr5=alt_utrs[0], utr3=alt_utrs[1])  # type: ignore[type-var]
    (primary_score,) = score_all(predictor, [cds])
    (alt_score,) = score_all(alt, [cds])
    difference = abs(primary_score - alt_score)
    return {
        "check": "utr-control",
        "cds_length_nt": len(cds),
        "utr5_length_nt": len(utrs[0]),
        "utr3_length_nt": len(utrs[1]),
        "utr5_alt_length_nt": len(alt_utrs[0]),
        "utr3_alt_length_nt": len(alt_utrs[1]),
        "score_primary_utrs": primary_score,
        "score_alt_utrs": alt_score,
        "abs_difference": difference,
        "harness_ok": difference > _NOISE_FLOOR,
        "note": (
            "Positive control only. harness_ok=False means the UTR context is not "
            "reaching the model, so no other result in this script is interpretable. "
            "harness_ok=True proves plumbing, not skill."
        ),
    }


def check_cds_spread(
    predictor: ExpressionPredictor,
    records: Sequence[tuple[str, str]],
    table: CodonUsageTable,
    tai_table: TaiTable,
) -> dict[str, Any]:
    """The decisive check: within-protein vs between-protein movement, plus confounds.

    Scores every record in **one** invocation with the UTRs held fixed, then for each
    protein reports the spread of scores across its own synonymous variants and the
    rank correlation of those scores against CAI, tAI, GC, GC3 and length -- all
    recomputed here by BT4's own functions (``reported == computed``).

    The headline number is ``within_over_between``: the median within-protein standard
    deviation divided by the standard deviation of the per-protein means. Near zero
    means RiboNN is effectively blind to synonymous change under a fixed UTR and is
    reading gene identity instead -- which is exactly what it was trained to do, and
    exactly what BT4 cannot use.

    Returns:
        A report with per-group spreads, per-group confound correlations, and the
        pooled summary.
    """
    grouped: dict[str, list[tuple[str, str]]] = {}
    for header, dna in records:
        grouped.setdefault(group_of(header), []).append((header, dna))

    flat = [dna for members in grouped.values() for _header, dna in members]
    scores = score_all(predictor, flat)
    by_dna = dict(zip(flat, scores, strict=True))

    features = ("cai", "tai", "gc", "gc3", "length_nt")
    groups: list[dict[str, Any]] = []
    for name, members in sorted(grouped.items()):
        dnas = [dna for _header, dna in members]
        group_scores = [by_dna[dna] for dna in dnas]
        measured = {
            "cai": [table.cai(dna) for dna in dnas],
            "tai": [tai_table.tai(dna) for dna in dnas],
            "gc": [gc_count(dna) / len(dna) for dna in dnas],
            "gc3": [gc3_fraction(dna) for dna in dnas],
            "length_nt": [float(len(dna)) for dna in dnas],
        }
        groups.append(
            {
                "group": name,
                "protein_length_aa": len(_protein_of(dnas[0])),
                "score": _spread(group_scores),
                # Spearman is 0.0 by construction when a series has no variance
                # (bt4.biomodels._stats returns an honest 0.0 rather than NaN), which
                # is the right reading: no detectable monotonic relationship.
                "confound_spearman": {
                    feature: spearman(measured[feature], group_scores)
                    for feature in features
                },
            }
        )

    within_sds = [g["score"]["sd"] for g in groups]
    group_means = [g["score"]["mean"] for g in groups]
    between_sd = statistics.stdev(group_means) if len(group_means) > 1 else 0.0
    median_within = statistics.median(within_sds) if within_sds else 0.0
    ratio = (median_within / between_sd) if between_sd > _NOISE_FLOOR else None

    responds = median_within > _NOISE_FLOOR
    gc3_rhos = [abs(g["confound_spearman"]["gc3"]) for g in groups]
    return {
        "check": "cds-spread",
        "covers": ["Step 1.3 (synonymous sensitivity)", "Step 1.4 (GC confound)"],
        "n_records": len(flat),
        "n_groups": len(groups),
        "groups": groups,
        "median_within_group_sd": median_within,
        "between_group_sd": between_sd,
        "within_over_between": ratio,
        "responds_to_synonymous_change": responds,
        "median_abs_gc3_spearman": statistics.median(gc3_rhos) if gc3_rhos else 0.0,
        "note": (
            "responds_to_synonymous_change=False is a decisive result for the backend "
            "under test: it cannot rank synonymous variants under a fixed UTR, so no "
            "panel will promote it. (For the null placeholder that is expected and is "
            "the reference for what 'blind' looks like.) If True, read "
            "median_abs_gc3_spearman before celebrating -- a response that tracks GC3 "
            "is a GC detector, and BT4 already has GC as an objective term for free."
        ),
    }


def check_direction(
    predictor: ExpressionPredictor,
    proteins: Sequence[str],
    table: CodonUsageTable,
    organism: str,
    reference_set: str | None,
    seed: int,
) -> dict[str, Any]:
    """Sign test: does RiboNN prefer a max-CAI design to a min-CAI one?

    For each protein BT4 builds two variants -- ``cai_weight=+1`` and
    ``cai_weight=-1``, everything else off -- and both are scored in one invocation
    with the UTRs fixed. The result is an exact two-sided binomial test against
    chance.

    Deliberately crude. Codon-optimized and codon-deoptimized sequences differ in GC
    as well as in codon choice, so a clear preference here is consistent with RiboNN
    being a GC detector; only ``cds-spread``'s confound correlations can separate
    them. A coin flip, on the other hand, is a strong warning.

    Returns:
        A per-protein table plus ``successes`` / ``trials`` / ``p_value``.
    """
    pairs: list[tuple[str, str, str]] = []
    for protein in proteins:
        best = api.optimize(
            protein,
            api.OptimizeConfig(
                organism=organism,
                reference_set=reference_set,
                cai_weight=1.0,
                gc_weight=0.0,
                max_homopolymer=None,
                seed=seed,
            ),
        )
        worst = api.optimize(
            protein,
            api.OptimizeConfig(
                organism=organism,
                reference_set=reference_set,
                cai_weight=-1.0,
                gc_weight=0.0,
                max_homopolymer=None,
                seed=seed,
            ),
        )
        pairs.append((protein, best.dna, worst.dna))

    flat = [dna for _p, high, low in pairs for dna in (high, low)]
    scores = score_all(predictor, flat)

    rows: list[dict[str, Any]] = []
    successes = 0
    ties = 0
    for i, (protein, high, low) in enumerate(pairs):
        high_score, low_score = scores[2 * i], scores[2 * i + 1]
        # A tie is not a failure. The sign test is defined on non-tied pairs, so ties
        # are counted and dropped rather than silently scored against the optimized
        # design -- otherwise a backend that is simply blind (every score identical)
        # would report "0/N prefer the optimized design", which reads as a strong
        # preference for the deoptimized one.
        if abs(high_score - low_score) <= _NOISE_FLOOR:
            ties += 1
            preferred = None
        else:
            preferred = high_score > low_score
            successes += int(preferred)
        rows.append(
            {
                "protein_length_aa": len(protein),
                "cai_max_cai_design": table.cai(high),
                "cai_min_cai_design": table.cai(low),
                "gc3_max_cai_design": gc3_fraction(high),
                "gc3_min_cai_design": gc3_fraction(low),
                "score_max_cai": high_score,
                "score_min_cai": low_score,
                "prefers_optimized": preferred,
            }
        )
    trials = len(rows) - ties
    return {
        "check": "direction",
        "pairs": len(rows),
        "ties": ties,
        "trials": trials,
        "successes": successes,
        "p_value": sign_test_p_value(successes, trials),
        "proteins": rows,
        "note": (
            "A crude smoke test, not evidence of skill: max-CAI and min-CAI designs "
            "differ in GC as well as codon choice, so a significant result is equally "
            "consistent with a GC detector. Read it together with cds-spread's "
            "confound correlations. A coin flip is the informative failure, and "
            "ties=pairs means the backend is blind to the difference entirely (which "
            "is what the null placeholder does, by construction). Ties are excluded "
            "from trials, as the sign test requires."
        ),
    }


def check_ladder(
    predictor: ExpressionPredictor,
    protein: str,
    organism: str,
    reference_set: str | None,
    steps: int,
    seed: int,
) -> dict[str, Any]:
    """Dose-response along a real BT4 Pareto frontier for one protein.

    Every frontier point is an exact, proven-optimal solve at a different CAI/GC
    trade-off, so the set is a designed ladder rather than random noise. If RiboNN's
    score moves monotonically with CAI along it, the response is at least orderly; a
    jagged response says the signal is unusable for ranking even when it is nonzero.

    Returns:
        The per-rung table and the Spearman of score against CAI along the ladder.
    """
    frontier = api.frontier(
        protein,
        api.OptimizeConfig(
            organism=organism, reference_set=reference_set, gc_weight=1.0, seed=seed
        ),
        steps=steps,
    )
    dnas: list[str] = []
    cais: list[float] = []
    table = load_table(organism, reference_set=reference_set)
    for result in frontier.results:
        if result.dna not in dnas:
            dnas.append(result.dna)
            cais.append(table.cai(result.dna))
    scores = score_all(predictor, dnas)
    return {
        "check": "ladder",
        "protein_length_aa": len(protein),
        "n_rungs": len(dnas),
        "rungs": [
            {"cai": cai, "score": score}
            for cai, score in sorted(zip(cais, scores, strict=True))
        ],
        "spearman_score_vs_cai": spearman(cais, scores) if len(dnas) > 1 else 0.0,
        "score_spread": _spread(scores),
        "note": (
            "Orderliness, not correctness: a high |Spearman| here says RiboNN tracks "
            "CAI along a designed ladder, which is a coherence check. It is not "
            "evidence that RiboNN's ordering matches measured expression."
        ),
    }


# --------------------------------------------------------------------------- #
# driver


def build_report(
    check: str,
    predictor: ExpressionPredictor,
    *,
    records: Sequence[tuple[str, str]],
    utrs: tuple[str, str],
    alt_utrs: tuple[str, str] | None,
    organism: str,
    reference_set: str | None,
    steps: int,
    seed: int,
    max_proteins: int,
) -> dict[str, Any]:
    """Run one check and return a JSON-ready report (no printing, no file writes).

    Args:
        check: One of :data:`CHECKS`.
        predictor: The expression backend, already carrying the fixed UTR context.
        records: ``(header, dna)`` pairs from the panel FASTA.
        utrs: The ``(utr5, utr3)`` the predictor was built with, for the record.
        alt_utrs: The second UTR pair, required by ``utr-control`` only.
        organism / reference_set: Codon-table selection, threaded everywhere so a
            reported CAI always names the question it answers (CLAUDE.md §8).
        steps: Frontier steps for ``ladder``.
        seed: Solver seed (determinism, invariant #7).
        max_proteins: Cap on proteins used by ``direction``.

    Returns:
        The report dict, always carrying a ``backend`` block so an uncalibrated
        result can never be read as a calibrated one.

    Raises:
        ValueError: If ``check`` is unknown, or a check's inputs are missing.
    """
    table = load_table(organism, reference_set=reference_set)
    header: dict[str, Any] = {
        "backend": {
            "name": predictor.name,
            "calibrated": predictor.calibrated,
            "units": "RiboNN CLR-residual TE (never exponentiated)",
        },
        "organism": organism,
        "reference_set": reference_set or table.reference_set,
        "seed": seed,
        "utr5_sha_prefix": _short_hash(utrs[0]),
        "utr3_sha_prefix": _short_hash(utrs[1]),
        "honesty": (
            "Diagnostics only. calibrated is False, no result here can promote a "
            "backend, and every score is an uncalibrated model output in arbitrary "
            "CLR-residual units (CLAUDE.md §6/§10.6)."
        ),
    }

    if check == "utr-control":
        if alt_utrs is None:
            raise ValueError("utr-control needs --utr5-alt and --utr3-alt")
        if not records:
            raise ValueError("utr-control needs a CDS (--fasta or --cds)")
        body = check_utr_control(predictor, records[0][1], utrs, alt_utrs)
    elif check == "cds-spread":
        if not records:
            raise ValueError("cds-spread needs --fasta")
        tai_table = load_tai_table(organism)
        body = check_cds_spread(predictor, records, table, tai_table)
    elif check == "direction":
        proteins = _panel_proteins(records)[:max_proteins]
        if not proteins:
            raise ValueError("direction needs --fasta with at least one CDS")
        body = check_direction(predictor, proteins, table, organism, reference_set, seed)
    elif check == "ladder":
        proteins = _panel_proteins(records)
        if not proteins:
            raise ValueError("ladder needs --fasta with at least one CDS")
        body = check_ladder(predictor, proteins[0], organism, reference_set, steps, seed)
    else:
        raise ValueError(f"unknown check {check!r}; choose from {list(CHECKS)}")

    return {**header, **body}


def _short_hash(text: str) -> str:
    """Return a 12-hex-char content hash -- identifies a UTR without printing it."""
    import hashlib

    return hashlib.sha256(text.encode("ascii")).hexdigest()[:12]


def _panel_proteins(records: Sequence[tuple[str, str]]) -> list[str]:
    """Return one protein per group, translated from that group's first record.

    Grouping first means a 93-record panel of three proteins yields three proteins,
    not 93 near-duplicates -- the same dependent-cluster reasoning that makes the
    protein the grouping unit everywhere else in this work.
    """
    seen: dict[str, str] = {}
    for header, dna in records:
        seen.setdefault(group_of(header), dna)
    proteins: list[str] = []
    for dna in seen.values():
        protein = _protein_of(dna)
        if protein:
            proteins.append(protein)
    return proteins


def _render(report: Mapping[str, Any]) -> str:
    """Render a report as a compact text table (the ``--json`` alternative)."""
    backend = report["backend"]
    flag = "calibrated" if backend["calibrated"] else "UNCALIBRATED"
    lines = [
        f"check:    {report['check']}",
        f"backend:  {backend['name']}  [{flag}]",
        f"organism: {report['organism']} (reference set: {report['reference_set']})",
        "",
    ]
    check = report["check"]
    if check == "utr-control":
        lines += [
            f"score with primary UTRs : {report['score_primary_utrs']:+.6f}",
            f"score with alt UTRs     : {report['score_alt_utrs']:+.6f}",
            f"absolute difference     : {report['abs_difference']:.6g}",
            "",
            f"harness_ok: {report['harness_ok']}",
        ]
    elif check == "cds-spread":
        lines.append(
            f"{'group':<12}{'n':>4}{'score sd':>12}{'score range':>14}"
            f"{'rho(GC3)':>10}{'rho(CAI)':>10}"
        )
        for group in report["groups"]:
            score, rho = group["score"], group["confound_spearman"]
            lines.append(
                f"{group['group'][:11]:<12}{score['n']:>4}{score['sd']:>12.6f}"
                f"{score['range']:>14.6f}{rho['gc3']:>10.3f}{rho['cai']:>10.3f}"
            )
        ratio = report["within_over_between"]
        if ratio is not None:
            ratio_text = f"{ratio:.4g}"
        elif report["n_groups"] < 2:
            ratio_text = "n/a (only one group in the panel)"
        else:
            ratio_text = "n/a (between-group sd is zero: no per-protein differences)"
        lines += [
            "",
            f"median within-group sd : {report['median_within_group_sd']:.6g}",
            f"between-group sd       : {report['between_group_sd']:.6g}",
            f"within / between       : {ratio_text}",
            f"median |rho(GC3)|      : {report['median_abs_gc3_spearman']:.3f}",
            "",
            f"responds_to_synonymous_change: {report['responds_to_synonymous_change']}",
        ]
    elif check == "direction":
        lines += [
            f"prefers the max-CAI design in {report['successes']}/{report['trials']} "
            f"non-tied proteins ({report['ties']} tie(s) of {report['pairs']} pairs, "
            "excluded)",
            f"exact two-sided sign test p = {report['p_value']:.4g}",
        ]
    elif check == "ladder":
        lines.append(f"{'CAI':>10}{'score':>14}")
        for rung in report["rungs"]:
            lines.append(f"{rung['cai']:>10.4f}{rung['score']:>14.6f}")
        lines += [
            "",
            f"spearman(score, CAI) = {report['spearman_score_vs_cai']:.3f}",
        ]
    lines += ["", f"note: {report['note']}", f"honesty: {report['honesty']}"]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else "",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--check", required=True, choices=CHECKS)
    parser.add_argument(
        "--fasta",
        default=str(DEFAULT_PANEL),
        help="panel FASTA; group key is the header text before the first '|' "
        f"(default: {DEFAULT_PANEL.name}, the in-tree CC BY 4.0 sequence-only panel)",
    )
    parser.add_argument("--cds", help="a single literal CDS, instead of --fasta")
    parser.add_argument(
        "--utr5", required=True, help="5' UTR: a literal ACGT string or a FASTA path"
    )
    parser.add_argument("--utr3", required=True, help="3' UTR: literal or FASTA path")
    parser.add_argument("--utr5-alt", help="second 5' UTR (utr-control only)")
    parser.add_argument("--utr3-alt", help="second 3' UTR (utr-control only)")
    parser.add_argument("--backend", default="ribonn", help="expression backend name")
    parser.add_argument("--species", default="human", choices=("human", "mouse"))
    parser.add_argument(
        "--cell-type",
        action="append",
        default=None,
        dest="cell_types",
        help="restrict to a RiboNN cell type (repeatable); default averages all",
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="RiboNN dataloader workers; 0 is required on Windows/macOS (spawn)",
    )
    parser.add_argument("--organism", default="homo_sapiens")
    parser.add_argument("--reference-set", default=None)
    parser.add_argument("--steps", type=int, default=11, help="frontier steps (ladder)")
    parser.add_argument("--max-proteins", type=int, default=25, help="cap (direction)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--json", action="store_true", help="emit JSON, not a table")
    args = parser.parse_args(argv)

    try:
        utrs = (read_flank(args.utr5), read_flank(args.utr3))
        alt_utrs = (
            (read_flank(args.utr5_alt), read_flank(args.utr3_alt))
            if args.utr5_alt and args.utr3_alt
            else None
        )
    except ValueError as exc:
        parser.error(str(exc))

    records = (
        [("cli", "".join(args.cds.split()).upper())]
        if args.cds
        else read_fasta(args.fasta)
    )

    predictor = resolve_backend(
        args.backend,
        species=args.species,
        utr5=utrs[0],
        utr3=utrs[1],
        top_k=args.top_k,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        cell_types=tuple(args.cell_types or ()),
        # Never honour a standing $BT4_EXPRESSION_USE_ATTESTED: these are the checks
        # that decide whether a panel is worth acquiring, so a prior attestation must
        # not colour them -- the same rule the gate itself follows.
        use_attested=False,
    )

    try:
        report = build_report(
            args.check,
            predictor,
            records=records,
            utrs=utrs,
            alt_utrs=alt_utrs,
            organism=args.organism,
            reference_set=args.reference_set,
            steps=args.steps,
            seed=args.seed,
            max_proteins=args.max_proteins,
        )
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=__import__("sys").stderr)
        return 2

    print(json.dumps(report, indent=2, sort_keys=True) if args.json else _render(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
