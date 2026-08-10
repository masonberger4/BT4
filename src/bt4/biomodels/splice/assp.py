"""ASSP-backed splice cross-check -- an opt-in, out-of-loop **network** validator.

:class:`AsspSplicePredictor` wraps the online **ASSP** (Alternative Splice Site
Predictor; Wang & Marin, *Nucleic Acids Research* 2006,
doi:10.1093/nar/gkl556) web service behind the existing
:class:`~bt4.biomodels.splice.base.SplicePredictor` contract. ASSP classifies
putative *constitutive*, *alternative-isoform*, and *cryptic* donor (5') and
acceptor (3') splice sites of an input sequence, each with a strength ``score``
and a ``confidence`` in ``[0, 1]``.

BT3's original splice "model" **scraped this exact service in the optimizer's
inner loop as its only splice path** -- fragile, non-reproducible, network-bound,
with no offline model behind it (CLAUDE.md section 10.15, the tool's cautionary
origin story). BT4 keeps ASSP as a *kept, supported feature* but inverts every one
of those properties; the rules are made structural here:

* **Opt-in and out-of-the-inner-loop.** ASSP is never scored per optimizer move.
  It runs only as a *final audit / validation pass* on an already-delivered
  sequence (``bt4 validate --splice-backend assp``, ``bt4 optimize --check-splice
  assp``), gated behind an explicit flag *and* the ``bt4[assp]`` extra (httpx).
  The default, reproducible splice path stays the local PWM baseline (and, when
  installed, the wrapped SpliceAI / Pangolin CNNs) -- ASSP is *never* returned by
  :func:`bt4.biomodels.splice.default` or
  :func:`bt4.pipeline.splice_audit.available_splice_backends`.

* **Never silent, never blocking.** Every ASSP request is rate-limited and retried
  with exponential backoff, and responses are cached by sequence hash
  (:class:`CachingAsspTransport`) so reruns and overlapping designs are free. If
  the service is unreachable or returns an unparseable body the *raw* predictor
  raises an :class:`AsspError` (honest failure, exactly like the CNN adapters
  raise when their weights are missing) -- but the graceful-degradation wrapper
  :func:`bt4.pipeline.splice_crosscheck.run_splice_crosscheck` catches it and
  reports "unavailable" so an ASSP failure can **never** fail an optimization.

* **Labeled network-derived and non-reproducible.** :attr:`network_derived` is
  ``True`` and :attr:`calibrated` is ``False``. ASSP numbers are stamped
  network-derived and are **excluded from BT4's reproducible-from-manifest
  guarantee** (unlike the hash-pinned local models): the cross-check is reported
  as a separate, clearly-labeled section and is never folded into a
  :class:`~bt4.domain.result.Result` audit or provenance manifest.

**Honest scope of the wire format.** The live :class:`HttpAsspTransport` targets
the ASSP web form and the tabular site report ASSP documents, but its exact
request/response wire format is **unverified against the live service** (which was
unreachable during development). CI never makes a live call: it drives the adapter
entirely from committed **offline fixtures** (:class:`FixtureAsspTransport`,
selected via the ``BT4_ASSP_FIXTURE_DIR`` environment variable). Those fixtures
are *synthetic, illustrative ASSP-format reports* used to exercise the
parsing / caching / backoff / degradation logic deterministically -- **not** real
captured responses (capturing one needs live network access CI forbids). This
mirrors how the wrapped CNN backends ship with no bundled reference panel: the
promotion path is a maintainer confirming the live transport against the real
service, exactly as a maintainer records the CNNs' fidelity gate.

The network dependency (``httpx``) is imported **only inside methods**, never at
module load, so importing this module -- and ``import bt4`` -- stays lightweight
(CLAUDE.md section 3). This module depends only on :mod:`bt4.domain`, the standard
library, and -- lazily, inside methods -- the optional ``httpx`` package.
"""

from __future__ import annotations

import hashlib
import os
import re
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from bt4.biomodels.splice.base import (
    DEFAULT_TOP_K,
    SpliceResult,
    pooled_risk,
)
from bt4.domain.sequence import validate_dna

__all__ = [
    "ASSP_ENDPOINT",
    "FIXTURE_DIR_ENV_VAR",
    "AsspError",
    "AsspReportError",
    "AsspSite",
    "AsspSplicePredictor",
    "AsspTransport",
    "AsspUnavailableError",
    "CachingAsspTransport",
    "FixtureAsspTransport",
    "HttpAsspTransport",
    "cache_key",
    "default_assp_transport",
    "parse_assp_report",
]


ASSP_ENDPOINT: str = "http://wangcomputing.com/assp/"
"""Documented base URL of the ASSP web service (Wang & Marin 2006).

**Unverified wire detail.** The live :class:`HttpAsspTransport` POSTs the sequence
here, but the exact CGI path and form-field names are not confirmed against the
live service (which was unreachable during development). A maintainer bringing the
live path online should confirm the endpoint and payload against the real form;
CI never touches it (it uses committed offline fixtures).
"""

FIXTURE_DIR_ENV_VAR: str = "BT4_ASSP_FIXTURE_DIR"
"""Env var selecting an **offline** ASSP transport backed by stored fixtures.

When set, :func:`default_assp_transport` returns a :class:`FixtureAsspTransport`
reading committed fixture reports from that directory instead of a live
:class:`HttpAsspTransport` -- so CI (and any offline run) drives the adapter with
zero network calls. Each fixture file is named ``<cache_key(sequence)>.txt``.
"""


class AsspError(Exception):
    """Base class for every ASSP-adapter failure (all are caught as degradations)."""


class AsspUnavailableError(AsspError):
    """The ASSP service (or the ``httpx`` dependency) could not be reached / used.

    Raised by a transport when the network request cannot be completed (dependency
    missing, connection error, non-2xx status, retries exhausted, or a missing
    offline fixture). The cross-check wrapper catches it and degrades gracefully --
    an ASSP outage never fails an optimization (CLAUDE.md section 15/10.15).
    """


class AsspReportError(AsspError):
    """The ASSP response could not be parsed into splice sites.

    Raised when the report body has no recognizable site table, or a site row is
    malformed (unparseable score / confidence / site type, or a position outside
    the submitted sequence). Also caught by the cross-check wrapper, so a garbled
    response degrades gracefully rather than crashing a run.
    """


@dataclass(frozen=True, slots=True)
class AsspSite:
    """One splice site ASSP predicted for a sequence (immutable).

    Attributes:
        position: 0-based anchor index of the site in the submitted coding
            sequence (ASSP's 1-based position minus one).
        kind: ``"donor"`` (5' splice site) or ``"acceptor"`` (3' splice site).
        site_class: ASSP's classification -- ``"constitutive"``, ``"alternative"``
            (alternative-isoform), ``"cryptic"``, or ``"unknown"`` if the report did
            not label it.
        score: ASSP's splice-site strength score (unbounded; larger = stronger).
        confidence: ASSP's confidence in ``[0, 1]`` -- an **uncalibrated,
            network-derived** pseudo-probability (never a calibrated splice
            probability; see the module docstring).
    """

    position: int
    kind: str
    site_class: str
    score: float
    confidence: float


def cache_key(dna: str) -> str:
    """Return the content-addressed cache / fixture key for a sequence.

    The lowercase hex SHA-256 of the upper-cased, validated coding sequence. Used
    to name both cache entries (:class:`CachingAsspTransport`) and offline fixtures
    (:class:`FixtureAsspTransport`), so the same sequence always maps to the same
    key regardless of input casing.

    Args:
        dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

    Returns:
        The 64-character hex digest keying this sequence.

    Raises:
        ValueError: If ``dna`` is empty or contains non-ACGT characters.
    """
    return hashlib.sha256(validate_dna(dna).encode("ascii")).hexdigest()


# ---------------------------------------------------------------------------
# Report parsing (pure, standard library only -- fully testable offline).
# ---------------------------------------------------------------------------

# Header-column aliases (case-insensitive substring match) locating each field in
# ASSP's tabular site report. The parser is header-driven so it tolerates extra
# columns and column reordering.
_POSITION_ALIASES: tuple[str, ...] = ("position", "pos")
_SCORE_ALIASES: tuple[str, ...] = ("score",)
_CONFIDENCE_ALIASES: tuple[str, ...] = ("confidence", "conf")
_TYPE_ALIASES: tuple[str, ...] = ("site", "type", "prediction")


def _split_columns(line: str) -> list[str]:
    """Split a report line into columns on tabs, else on runs of 2+ spaces.

    ASSP's tabular report is tab-delimited; a space-padded fixed-width variant is
    also tolerated (runs of two or more spaces separate columns) so a single-spaced
    ``site type`` label is not split mid-token. On the tab path, **empty cells are
    preserved** (only stripped): the parser maps columns by header index, so
    dropping an empty field would shift every later column left and misalign the
    row against the header. The space-padded path cannot represent an empty middle
    column (a run of spaces is a single delimiter), so its stray leading/trailing
    empties are dropped.
    """
    if "\t" in line:
        return [cell.strip() for cell in line.split("\t")]
    return [cell.strip() for cell in re.split(r"\s{2,}", line.strip()) if cell.strip() != ""]


def _match_column(cells: Sequence[str], aliases: Sequence[str]) -> int | None:
    """Return the index of the first cell whose lowercase text contains an alias."""
    for idx, cell in enumerate(cells):
        low = cell.lower()
        if any(alias in low for alias in aliases):
            return idx
    return None


@dataclass(frozen=True, slots=True)
class _Header:
    """Resolved column indices for the site table (from the header row)."""

    position: int
    score: int
    confidence: int
    site_type: int


def _find_header(lines: Sequence[str]) -> tuple[int, _Header] | None:
    """Find the site-table header row and resolve its column indices.

    A header is the first line whose columns name (by alias) a position, a score,
    a confidence, and a site-type column. Returns ``(line_index, header)`` or
    ``None`` when no such row exists.
    """
    for line_index, line in enumerate(lines):
        cells = _split_columns(line)
        if len(cells) < 4:
            continue
        pos = _match_column(cells, _POSITION_ALIASES)
        score = _match_column(cells, _SCORE_ALIASES)
        conf = _match_column(cells, _CONFIDENCE_ALIASES)
        site = _match_column(cells, _TYPE_ALIASES)
        if None in (pos, score, conf, site):
            continue
        # mypy: the None check above proves these are ints.
        return line_index, _Header(position=pos, score=score, confidence=conf, site_type=site)  # type: ignore[arg-type]
    return None


def _classify(type_cell: str) -> tuple[str, str] | None:
    """Map a site-type cell to ``(kind, site_class)``, or ``None`` if unrecognized.

    ``kind`` is ``"donor"`` / ``"acceptor"`` (recognizing ASSP's ``5'`` / ``3'``
    aliases); ``site_class`` is ``constitutive`` / ``alternative`` / ``cryptic`` /
    ``unknown``.
    """
    low = type_cell.lower()
    if "donor" in low or "5'" in low or "5 '" in low or "5ss" in low:
        kind = "donor"
    elif "acceptor" in low or "3'" in low or "3 '" in low or "3ss" in low:
        kind = "acceptor"
    else:
        return None
    if "constitutive" in low:
        site_class = "constitutive"
    elif "alternative" in low or "isoform" in low:
        site_class = "alternative"
    elif "cryptic" in low:
        site_class = "cryptic"
    else:
        site_class = "unknown"
    return kind, site_class


def parse_assp_report(report: str, seq_len: int) -> tuple[AsspSite, ...]:
    """Parse an ASSP tabular site report into :class:`AsspSite` records.

    Header-driven and tolerant (see :func:`_split_columns`): the parser locates the
    site-table header (a row naming position / score / confidence / site-type
    columns), then reads each subsequent row whose position column is an integer.
    Rows whose position column is not an integer are treated as prose / footer and
    skipped; a row with an integer position but an unparseable score, confidence, or
    site type -- or a position outside ``[1, seq_len]`` -- is a genuine malformed
    row and raises.

    Args:
        report: The raw ASSP response body.
        seq_len: Length of the submitted sequence (for 1-based position bounds).

    Returns:
        The parsed sites, ordered by ``(position, kind)``. Empty when the report has
        a header but lists no sites (a legitimate "no predicted sites" result).

    Raises:
        AsspReportError: If no site-table header is found (an unrecognized /
            empty / error body), or a data row is malformed.
    """
    lines = report.splitlines()
    found = _find_header(lines)
    if found is None:
        raise AsspReportError(
            "ASSP response has no recognizable site table (expected a header row "
            "naming position / score / confidence / site-type columns)"
        )
    header_index, header = found
    sites: list[AsspSite] = []
    max_index = max(header.position, header.score, header.confidence, header.site_type)
    for line in lines[header_index + 1 :]:
        cells = _split_columns(line)
        if len(cells) <= max_index:
            continue  # prose / footer / short line -- not a full site row
        try:
            position_1based = int(cells[header.position])
        except ValueError:
            continue  # position column is not an integer -> not a site row
        if not 1 <= position_1based <= seq_len:
            raise AsspReportError(
                f"ASSP site position {position_1based} is outside the submitted "
                f"sequence [1, {seq_len}]"
            )
        try:
            score = float(cells[header.score])
            confidence = float(cells[header.confidence])
        except ValueError as exc:
            raise AsspReportError(f"ASSP site row has an unparseable number: {line!r}") from exc
        classified = _classify(cells[header.site_type])
        if classified is None:
            raise AsspReportError(f"ASSP site row has an unrecognized site type: {line!r}")
        kind, site_class = classified
        sites.append(
            AsspSite(
                position=position_1based - 1,
                kind=kind,
                site_class=site_class,
                score=score,
                confidence=confidence,
            )
        )
    sites.sort(key=lambda s: (s.position, s.kind))
    return tuple(sites)


# ---------------------------------------------------------------------------
# Transports: how a sequence turns into a raw ASSP report body.
# ---------------------------------------------------------------------------


@runtime_checkable
class AsspTransport(Protocol):
    """Maps a validated coding sequence to a raw ASSP report body.

    The seam between the ASSP adapter and the outside world, so the network can be
    swapped for offline fixtures in tests / CI. Implementations must raise
    :class:`AsspUnavailableError` when they cannot produce a body (network error,
    missing dependency, missing fixture) rather than returning junk.
    """

    def fetch(self, dna: str) -> str:
        """Return the raw ASSP report body for ``dna`` (or raise ``AsspError``)."""
        ...

    def available(self) -> bool:
        """Return whether this transport could run (cheap, never raises)."""
        ...


def _import_httpx() -> Any:
    """Import and return the ``httpx`` module (lazy, guarded).

    Raises:
        AsspUnavailableError: If ``httpx`` is not installed -- the ``bt4[assp]``
            extra provides it. Raised (not ``ModuleNotFoundError``) so the
            cross-check wrapper degrades gracefully.
    """
    try:
        import httpx  # type: ignore[import-not-found]

        return httpx
    except ImportError as exc:
        raise AsspUnavailableError(
            "httpx is not installed; install the 'bt4[assp]' extra to use the "
            "opt-in ASSP splice cross-check, or use an offline splice backend"
        ) from exc


def _throttle(
    state: dict[str, float],
    min_interval_s: float,
    monotonic: Callable[[], float],
    sleep: Callable[[float], None],
) -> None:
    """Sleep just long enough to keep requests at most one per ``min_interval_s``.

    Polite rate-limiting: records the last request time in ``state`` and sleeps the
    remaining interval before allowing the next. Pure of ``httpx``; the clock and
    sleep are injected so the throttle is unit-testable without real waiting.
    """
    if min_interval_s <= 0:
        return
    last = state.get("last")
    now = monotonic()
    if last is not None:
        wait = min_interval_s - (now - last)
        if wait > 0:
            sleep(wait)
    state["last"] = monotonic()


def _with_retries(
    attempt_fn: Callable[[], str],
    *,
    attempts: int,
    backoff_base_s: float,
    sleep: Callable[[float], None],
    retryable: tuple[type[BaseException], ...],
) -> str:
    """Call ``attempt_fn`` with exponential backoff on ``retryable`` exceptions.

    Retries up to ``attempts`` total tries, sleeping ``backoff_base_s * 2**i``
    between them; re-raises the last error as :class:`AsspUnavailableError` when all
    tries are exhausted. Pure of ``httpx`` (the retryable types and ``sleep`` are
    injected), so it is unit-testable without the network.

    Args:
        attempt_fn: The zero-arg request to retry; returns the report body.
        attempts: Maximum number of tries (``>= 1``).
        backoff_base_s: Base backoff; the delay before try ``i`` (0-based) is
            ``backoff_base_s * 2**i``.
        sleep: Injected sleep (``time.sleep`` in production).
        retryable: Exception types that trigger a retry (e.g. ``httpx.HTTPError``).

    Returns:
        The successful attempt's report body.

    Raises:
        ValueError: If ``attempts < 1``.
        AsspUnavailableError: If every try raised a ``retryable`` error.
    """
    if attempts < 1:
        raise ValueError(f"attempts must be >= 1, got {attempts}")
    last_exc: BaseException | None = None
    for i in range(attempts):
        try:
            return attempt_fn()
        except retryable as exc:
            last_exc = exc
            if i + 1 < attempts:
                sleep(backoff_base_s * (2**i))
    raise AsspUnavailableError(
        f"ASSP request failed after {attempts} attempt(s): {last_exc}"
    ) from last_exc


@dataclass(frozen=True, slots=True)
class HttpAsspTransport:
    """Live ASSP transport: POST the sequence to the web service via ``httpx``.

    Rate-limited (:func:`_throttle`) and retried with exponential backoff
    (:func:`_with_retries`); wraps every network failure in
    :class:`AsspUnavailableError`. See the module docstring for the honest
    "wire format unverified" caveat -- CI never exercises this path.

    Attributes:
        endpoint: The ASSP form URL (default :data:`ASSP_ENDPOINT`).
        timeout_s: Per-request timeout in seconds.
        max_attempts: Total tries before giving up (``>= 1``).
        backoff_base_s: Base for the exponential backoff between retries.
        min_interval_s: Minimum seconds between successive requests (polite
            rate-limiting).
    """

    endpoint: str = ASSP_ENDPOINT
    timeout_s: float = 30.0
    max_attempts: int = 4
    backoff_base_s: float = 2.0
    min_interval_s: float = 1.0
    # Mutable rate-limit clock state and injected clock/sleep -- excluded from
    # equality/repr so two identically-configured transports still compare equal.
    _state: dict[str, float] = field(default_factory=dict, compare=False, repr=False)
    _sleep: Callable[[float], None] = field(default=time.sleep, compare=False, repr=False)
    _monotonic: Callable[[], float] = field(default=time.monotonic, compare=False, repr=False)

    def available(self) -> bool:
        """Return whether ``httpx`` is importable (does not probe the network)."""
        try:
            _import_httpx()
        except AsspError:
            return False
        return True

    def _payload(self, dna: str) -> dict[str, str]:
        """Return the documented ASSP form payload for ``dna`` (unverified fields)."""
        return {"sequence": dna}

    def fetch(self, dna: str) -> str:
        """POST ``dna`` to ASSP and return the raw report body.

        Args:
            dna: A validated coding sequence over ``{A,C,G,T}``.

        Returns:
            The raw ASSP response text.

        Raises:
            AsspUnavailableError: If ``httpx`` is missing, or the request fails on
                every retry.
        """
        httpx = _import_httpx()
        _throttle(self._state, self.min_interval_s, self._monotonic, self._sleep)

        def attempt() -> str:
            response = httpx.post(self.endpoint, data=self._payload(dna), timeout=self.timeout_s)
            response.raise_for_status()
            text = response.text
            if not isinstance(text, str):  # pragma: no cover - defensive
                raise httpx.HTTPError("ASSP response body was not text")
            return text

        return _with_retries(
            attempt,
            attempts=self.max_attempts,
            backoff_base_s=self.backoff_base_s,
            sleep=self._sleep,
            retryable=(httpx.HTTPError,),
        )


@dataclass(frozen=True, slots=True)
class FixtureAsspTransport:
    """Offline ASSP transport: read stored report bodies from a fixtures directory.

    Reads ``<directory>/<cache_key(dna)>.txt`` -- so CI (and any offline run) drives
    the adapter with zero network calls. A missing fixture raises
    :class:`AsspUnavailableError`, so the graceful-degradation path is testable
    (point at an empty directory).

    Attributes:
        directory: The directory holding ``<cache_key>.txt`` fixture reports.
    """

    directory: str

    def available(self) -> bool:
        """Return whether the fixtures directory exists (never raises)."""
        return Path(self.directory).is_dir()

    def fetch(self, dna: str) -> str:
        """Return the stored report body for ``dna``.

        Raises:
            AsspUnavailableError: If no fixture file exists for ``dna``'s cache key.
        """
        path = Path(self.directory) / f"{cache_key(dna)}.txt"
        try:
            return path.read_text(encoding="utf-8")
        except OSError as exc:
            raise AsspUnavailableError(
                f"no offline ASSP fixture for this sequence at {path}"
            ) from exc


@dataclass(frozen=True, slots=True)
class CachingAsspTransport:
    """Wrap a transport with an in-memory cache keyed by sequence hash.

    "Responses cached by sequence hash" (CLAUDE.md section 6): a sequence is fetched
    from the wrapped transport at most once per process, so reruns and overlapping
    designs cost nothing. Only successful fetches are cached (a failure is not
    memoized, so a transient outage can recover on retry).

    Attributes:
        inner: The wrapped transport that actually produces report bodies.
    """

    inner: AsspTransport
    _cache: dict[str, str] = field(default_factory=dict, compare=False, repr=False)

    def available(self) -> bool:
        """Delegate to the wrapped transport (never raises)."""
        return self.inner.available()

    def fetch(self, dna: str) -> str:
        """Return the cached body for ``dna``, fetching (and caching) on a miss."""
        key = cache_key(dna)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        body = self.inner.fetch(dna)
        self._cache[key] = body
        return body


def default_assp_transport() -> AsspTransport:
    """Return the default ASSP transport (offline fixtures if configured, else live).

    If :data:`FIXTURE_DIR_ENV_VAR` is set, returns a cache-wrapped
    :class:`FixtureAsspTransport` over that directory -- the fully-offline path CI
    uses. Otherwise returns a cache-wrapped live :class:`HttpAsspTransport`. Either
    way caching by sequence hash is on.
    """
    fixture_dir = os.environ.get(FIXTURE_DIR_ENV_VAR)
    inner: AsspTransport = (
        FixtureAsspTransport(fixture_dir) if fixture_dir else HttpAsspTransport()
    )
    return CachingAsspTransport(inner)


# ---------------------------------------------------------------------------
# The predictor.
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AsspSplicePredictor:
    """Opt-in, out-of-loop ASSP splice cross-check behind the ``SplicePredictor`` contract.

    Turns ASSP's predicted donor / acceptor sites into per-position tracks: a site's
    :attr:`AsspSite.confidence` (in ``[0, 1]``) is placed at its position in the
    ``donor`` or ``acceptor`` track. See the module docstring for the opt-in /
    out-of-loop / network-derived / non-reproducible rules -- this backend is never
    the default and never runs in the optimizer loop.

    Attributes:
        transport: How a sequence becomes an ASSP report body. Defaults (via
            :func:`default_assp_transport`) to a cache-wrapped live HTTP transport,
            or a cache-wrapped offline fixture transport when
            :data:`FIXTURE_DIR_ENV_VAR` is set.
        top_k: Number of strongest sites summed by :meth:`delta_splicing`'s top-k /
            log-odds pooling. Defaults to
            :data:`~bt4.biomodels.splice.base.DEFAULT_TOP_K`.
    """

    transport: AsspTransport = field(default=None)  # type: ignore[assignment]
    top_k: int = DEFAULT_TOP_K

    def __post_init__(self) -> None:
        """Resolve the default transport and validate the pooling depth.

        Raises:
            ValueError: If ``top_k`` is not a positive integer.
        """
        if self.transport is None:
            object.__setattr__(self, "transport", default_assp_transport())
        if self.top_k <= 0:
            raise ValueError(f"top_k must be a positive integer, got {self.top_k}")

    @property
    def name(self) -> str:
        """Backend identifier."""
        return "assp"

    @property
    def calibrated(self) -> bool:
        """Always ``False``: ASSP is an external heuristic, uncalibrated in BT4's ledger."""
        return False

    @property
    def network_derived(self) -> bool:
        """Always ``True``: ASSP numbers are network-derived and excluded from the manifest.

        Not part of the :class:`~bt4.biomodels.splice.base.SplicePredictor` protocol
        (the local backends are all reproducible-from-manifest); the cross-check
        wrapper reads it to stamp ASSP results network-derived and keep them out of
        the reproducible-from-manifest guarantee (CLAUDE.md section 6/10.15).
        """
        return True

    def available(self) -> bool:
        """Return whether the transport could run (deps / fixtures present; never raises).

        A cheap check that does **not** probe the live service -- reachability is
        discovered at call time and degrades gracefully. For the default live
        transport this reports whether ``httpx`` is installed.
        """
        return self.transport.available()

    def sites(self, dna: str) -> tuple[AsspSite, ...]:
        """Return ASSP's predicted splice sites for ``dna``.

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            The parsed :class:`AsspSite` records, ordered by ``(position, kind)``.

        Raises:
            AsspUnavailableError: If the service / dependency is unavailable.
            AsspReportError: If the response cannot be parsed.
            ValueError: If ``dna`` is empty or contains non-ACGT characters.
        """
        seq = validate_dna(dna)
        report = self.transport.fetch(seq)
        return parse_assp_report(report, len(seq))

    def score_sequence(self, dna: str) -> SpliceResult:
        """Return per-position donor / acceptor site scores as a :class:`SpliceResult`.

        Each predicted site's confidence is placed at its position in the matching
        track (donor or acceptor); positions with no predicted site score ``0.0``.
        Scores are **uncalibrated, network-derived** pseudo-probabilities
        (:attr:`calibrated` is always ``False``).

        Args:
            dna: A coding sequence over ``{A,C,G,T}`` (case-insensitive).

        Returns:
            A :class:`SpliceResult` with ``model_name = "assp"`` and
            ``calibrated = False``.

        Raises:
            AsspUnavailableError / AsspReportError / ValueError: As on :meth:`sites`.
        """
        seq = validate_dna(dna)
        n = len(seq)
        donor = [0.0] * n
        acceptor = [0.0] * n
        for site in self.sites(seq):
            track = donor if site.kind == "donor" else acceptor
            # A position could carry both a donor and an acceptor call in different
            # rows; keep the stronger confidence per track (never sum -- that is not
            # a probability).
            if site.confidence > track[site.position]:
                track[site.position] = site.confidence
        return SpliceResult(
            donor=tuple(donor),
            acceptor=tuple(acceptor),
            model_name=self.name,
            calibrated=False,
        )

    def delta_splicing(self, designed_dna: str, reference_dna: str) -> float:
        """Return the negated added splice risk of ``designed`` vs ``reference``.

        See :meth:`bt4.biomodels.splice.base.SplicePredictor.delta_splicing` for the
        fixed *larger-is-better* orientation. Concretely returns
        ``pooled_risk(reference) - pooled_risk(designed)`` using top-k / log-odds
        pooling, so it is ``0.0`` for identical sequences, positive when the redesign
        lowers ASSP-predicted splice risk, and negative when it raises it.

        Raises:
            AsspUnavailableError / AsspReportError / ValueError: As on :meth:`sites`.
        """
        designed_risk = pooled_risk(self.score_sequence(designed_dna), self.top_k)
        reference_risk = pooled_risk(self.score_sequence(reference_dna), self.top_k)
        return reference_risk - designed_risk
