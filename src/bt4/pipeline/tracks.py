"""Per-site risk tracks: sliding-window composition profiles of a sequence.

BT4 Studio's design (CLAUDE.md §6.6) calls for *per-site risk tracks* along the
delivered sequence so a user can see where composition risk concentrates rather
than only a whole-sequence number. This module assembles those tracks from the
honest, already-shipped **reporting** profiles -- it computes nothing new about
optimality and adds no objective; every value is a deterministic sliding-window
statistic recomputed from the DNA:

* ``gc_fraction`` -- GC fraction over a nucleotide window (occlusion / synthesis).
* ``cpg_density`` -- CpG density over a nucleotide window (innate-immune / stealth
  signal), from :func:`~bt4.objectives.dinuc_profile.dinucleotide_profile`.
* ``minmax`` -- the Clarke & Clark sliding-window %MinMax over a codon window
  (translational ramp / rare-codon clusters), from
  :func:`~bt4.objectives.minmax.min_max_profile`, using the organism's codon
  frequencies.
* ``splice_site`` -- per-nucleotide cryptic-splice-site score (the larger of the
  donor and acceptor score at each position), from whichever
  :class:`~bt4.biomodels.splice.SplicePredictor` backend is asked for. Unlike the
  three above this one comes from a *model*, so it carries that model's calibration
  status: with the default PWM baseline the numbers are **uncalibrated
  pseudo-scores in arbitrary units** and must be read as "where the consensus-like
  positions are", never as a splice probability (CLAUDE.md §10.6).

Each track carries its window size and unit so the caller can plot or audit it
without guessing. This is reporting, not optimization; nothing here feeds the
solver (see the profiles' own module docstrings and the %MinMax / CpG
additive-vs-reporting split).

There is deliberately **no folding track**. A sliding-window folding deltaG would
have to fold one window per position, and the dependency-free baseline is
``O(w^3)`` per window -- seconds to minutes for a real CDS -- so a track would
either be unusably slow or quietly computed at a coarser stride than its axis
implies. The 5' folding number the refinement pass actually optimizes is reported
directly in the result audit instead.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt4._accel import gc_count
from bt4.biomodels.codon.tables import load_table
from bt4.domain.context import ConstructContext
from bt4.domain.sequence import validate_dna
from bt4.objectives.dinuc_profile import dinucleotide_profile
from bt4.objectives.minmax import min_max_profile

__all__ = ["Track", "TracksResult", "run_tracks", "summarize"]


@dataclass(frozen=True, slots=True)
class Track:
    """One named sliding-window reporting profile over a sequence.

    Attributes:
        name: Stable track identifier (e.g. ``"gc_fraction"``).
        window: Window length; ``unit`` says whether it is nucleotides or codons.
        unit: The value's unit (e.g. ``"fraction"``, ``"cpg/slot"``, ``"%minmax"``).
        window_unit: ``"nt"`` or ``"codon"`` -- the stride/window basis.
        values: Per-window values; entry ``i`` describes the window starting at
            position ``i`` (in ``window_unit`` units). Empty when the sequence is
            shorter than one window.
    """

    name: str
    window: int
    unit: str
    window_unit: str
    values: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class TracksResult:
    """A bundle of per-site tracks recomputed from one coding sequence.

    Attributes:
        dna: The sequence the tracks were computed from.
        tracks: The named profiles.
        organism: The codon table's organism, or ``""`` when no table was read
            (the %MinMax track is omitted for a non-codon-aligned sequence).
        reference_set: Which of that organism's reference sets supplied the
            %MinMax frequencies, or ``""`` when no table was read. %MinMax is a
            *codon-commonness* profile, so the same sequence yields a different
            track under each reference set -- a result that did not carry the
            label would be a plot with no stated question, which is the failure
            this whole axis exists to prevent.
    """

    dna: str
    tracks: tuple[Track, ...]
    organism: str = ""
    reference_set: str = ""
    splice_model: str = ""
    splice_calibrated: bool = False

    def get(self, name: str) -> Track | None:
        """Return the track named ``name``, or ``None`` if absent."""
        return next((t for t in self.tracks if t.name == name), None)


def _gc_fraction_profile(dna: str, window: int) -> list[float]:
    """Per-window GC fraction over a ``window``-nt sliding window.

    Args:
        dna: An upper-cased ACGT sequence.
        window: Window length in nucleotides (``>= 1``).

    Returns:
        ``gc_count(dna[i : i + window]) / window`` for each window start ``i``;
        empty when ``len(dna) < window``.
    """
    n = len(dna)
    if window < 1 or n < window:
        return []
    return [gc_count(dna[i : i + window]) / window for i in range(n - window + 1)]


def run_tracks(
    dna: str,
    organism: str = "homo_sapiens",
    *,
    reference_set: str | None = None,
    nt_window: int = 50,
    codon_window: int = 18,
    splice: bool = False,
    splice_predictor: object | None = None,
    context: ConstructContext | None = None,
) -> TracksResult:
    """Compute the per-site composition tracks for a coding sequence.

    Args:
        dna: An ACGT coding sequence (case-insensitive). GC and CpG tracks accept
            any length; the %MinMax track needs a length that is a multiple of
            three (codon-aligned) and is omitted otherwise.
        organism: Codon-usage table key/alias whose frequencies define the
            %MinMax reference.
        reference_set: Which of that organism's reference sets to read the
            frequencies from (``None`` = its default). %MinMax is a *codon
            commonness* profile, so which genes the frequencies were counted
            over changes what the track means.
        nt_window: Sliding-window length (nucleotides) for the GC and CpG tracks.
        codon_window: Sliding-window length (codons) for the %MinMax track.
        splice: When ``True``, add the per-nucleotide ``splice_site`` track. It is
            opt-in because it runs a model rather than a composition statistic, and
            because the resulting numbers inherit that model's calibration status.
        splice_predictor: The :class:`~bt4.biomodels.splice.SplicePredictor` to
            score with; ``None`` uses :func:`bt4.biomodels.splice.default` (the
            uncalibrated PWM baseline). The result records which model ran and
            whether it is calibrated, so a plot can never imply more than the
            backend claims.

    Returns:
        A :class:`TracksResult` with the ``gc_fraction`` and ``cpg_density``
        tracks always, plus ``minmax`` when the sequence is codon-aligned and
        ``splice_site`` when ``splice`` is set.

    Raises:
        ValueError: On non-ACGT input, or a non-positive window.
    """
    d = validate_dna(dna)
    if nt_window < 1 or codon_window < 1:
        raise ValueError("window lengths must be >= 1")

    tracks: list[Track] = [
        Track(
            name="gc_fraction",
            window=nt_window,
            unit="fraction",
            window_unit="nt",
            values=tuple(_gc_fraction_profile(d, nt_window)),
        ),
        Track(
            name="cpg_density",
            window=nt_window,
            unit="cpg/slot",
            window_unit="nt",
            values=tuple(dinucleotide_profile(d, "CG", nt_window, density=True)),
        ),
    ]
    resolved_organism = ""
    resolved_reference_set = ""
    if len(d) % 3 == 0:
        table = load_table(organism, reference_set=reference_set)
        resolved_organism = table.organism
        resolved_reference_set = table.reference_set
        frequencies = table.frequency
        tracks.append(
            Track(
                name="minmax",
                window=codon_window,
                unit="%minmax",
                window_unit="codon",
                values=tuple(min_max_profile(d, frequencies, codon_window)),
            )
        )
    splice_model = ""
    splice_calibrated = False
    if splice:
        # Lazy: keeps `import bt4` light and the splice backends optional.
        from bt4.biomodels.splice import SplicePredictor, score_in_context
        from bt4.biomodels.splice import default as splice_default

        predictor = splice_predictor if splice_predictor is not None else splice_default()
        assert isinstance(predictor, SplicePredictor)
        # With known flanks the model reads real sequence around the CDS instead of
        # padding; the scores come back aligned to the CDS either way.
        scored = (
            score_in_context(predictor, d, context.upstream, context.downstream)
            if context is not None and not context.is_empty
            else predictor.score_sequence(d)
        )
        splice_model = scored.model_name
        splice_calibrated = scored.calibrated
        tracks.append(
            Track(
                name="splice_site",
                window=1,
                # The unit names the calibration status so a plot legend cannot
                # quietly upgrade a PWM pseudo-score into a probability.
                unit="p(site)" if scored.calibrated else "score (UNCALIBRATED)",
                window_unit="nt",
                values=tuple(
                    max(don, acc)
                    for don, acc in zip(scored.donor, scored.acceptor, strict=True)
                ),
            )
        )
    return TracksResult(
        dna=d,
        tracks=tuple(tracks),
        organism=resolved_organism,
        reference_set=resolved_reference_set,
        splice_model=splice_model,
        splice_calibrated=splice_calibrated,
    )


def summarize(tracks: Sequence[Track]) -> list[dict[str, object]]:
    """Return a compact per-track summary (count, min, max, mean) for display."""
    out: list[dict[str, object]] = []
    for track in tracks:
        values = track.values
        out.append(
            {
                "name": track.name,
                "window": track.window,
                "window_unit": track.window_unit,
                "unit": track.unit,
                "n": len(values),
                "min": min(values) if values else None,
                "max": max(values) if values else None,
                "mean": sum(values) / len(values) if values else None,
            }
        )
    return out
