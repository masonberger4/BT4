"""Tests for the opt-in, out-of-loop ASSP splice cross-check backend.

The ASSP adapter talks to an online service; CI never makes a live call. These
tests drive it entirely from **synthetic offline fixtures** and injected
transports (mirroring how ``test_splice_pangolin.py`` never touches the real GPL
weights), covering:

* the contract surface and honesty flags (``calibrated is False``,
  ``network_derived is True``, ``name == "assp"``), and that ASSP is *never* the
  default and *never* auto-included in ``available_splice_backends``;
* the tolerant, header-driven report parser (tab- and space-delimited, donor /
  acceptor, constitutive / alternative / cryptic, empty-but-valid, and every
  malformed-report refusal);
* the transports -- fixture lookup by sequence hash, hash-keyed caching, missing
  fixture -> graceful ``AsspUnavailableError``, and the fixture-vs-live default
  selection;
* the pure retry/backoff and rate-limit helpers (no real waiting, no ``httpx``);
* end-to-end scoring against a committed fixture; and
* the guarantee that importing the module does not pull ``httpx``.
"""

from __future__ import annotations

import subprocess
import sys
import types
from pathlib import Path

import pytest

from bt4.biomodels.splice import assp as assp_mod
from bt4.biomodels.splice import default
from bt4.biomodels.splice.assp import (
    FIXTURE_DIR_ENV_VAR,
    AsspReportError,
    AsspSplicePredictor,
    AsspUnavailableError,
    CachingAsspTransport,
    FixtureAsspTransport,
    HttpAsspTransport,
    _throttle,
    _with_retries,
    cache_key,
    default_assp_transport,
    parse_assp_report,
)
from bt4.biomodels.splice.base import SplicePredictor, pooled_risk
from bt4.biomodels.splice.baseline import ConsensusPwmSplicePredictor
from bt4.pipeline.splice_audit import available_splice_backends

FIXTURES = Path(__file__).parent / "fixtures" / "assp"
SEQ_WITH_SITES = "ATGGCCGGCGATCGATCGATCGTAA"  # committed fixture: three sites
SEQ_NO_SITES = "ATGAAATTTGGGCCCTAA"  # committed fixture: header, no site rows


def _fixture_predictor() -> AsspSplicePredictor:
    """An ASSP predictor whose transport reads the committed offline fixtures."""
    return AsspSplicePredictor(transport=CachingAsspTransport(FixtureAsspTransport(str(FIXTURES))))


# --------------------------------------------------------------------------- #
# Contract surface & honesty
# --------------------------------------------------------------------------- #


def test_adapter_is_a_splice_predictor() -> None:
    model = _fixture_predictor()
    assert isinstance(model, SplicePredictor)
    assert model.name == "assp"
    assert isinstance(model.calibrated, bool)


def test_honesty_flags() -> None:
    model = _fixture_predictor()
    # ASSP is an external heuristic: never calibrated, always network-derived.
    assert model.calibrated is False
    assert model.network_derived is True


def test_assp_is_never_default_and_never_auto_included() -> None:
    # default() stays the labeled PWM baseline -- ASSP is opt-in, never the default.
    assert isinstance(default(), ConsensusPwmSplicePredictor)
    # available_splice_backends() never includes the network backend (it is opt-in,
    # requested explicitly by name, never auto-discovered).
    assert all(b.name != "assp" for b in available_splice_backends())


def test_top_k_validation() -> None:
    with pytest.raises(ValueError):
        AsspSplicePredictor(transport=FixtureAsspTransport(str(FIXTURES)), top_k=0)


# --------------------------------------------------------------------------- #
# Report parser
# --------------------------------------------------------------------------- #

_TAB_REPORT = (
    "ASSP report\n"
    "position\tsite type\tscore\tconfidence\n"
    "5\tconstitutive donor\t8.10\t0.95\n"
    "12\tcryptic acceptor\t6.20\t0.91\n"
)


def test_parse_tab_report() -> None:
    sites = parse_assp_report(_TAB_REPORT, seq_len=25)
    assert len(sites) == 2
    donor, acceptor = sites
    assert (donor.position, donor.kind, donor.site_class) == (4, "donor", "constitutive")
    assert donor.score == pytest.approx(8.10)
    assert donor.confidence == pytest.approx(0.95)
    assert (acceptor.position, acceptor.kind, acceptor.site_class) == (11, "acceptor", "cryptic")


def test_parse_space_padded_report() -> None:
    report = (
        "position    site type              score    confidence\n"
        "3           alternative isoform donor   5.5      0.90\n"
    )
    (site,) = parse_assp_report(report, seq_len=25)
    assert (site.position, site.kind, site.site_class) == (2, "donor", "alternative")


def test_parse_five_prime_three_prime_aliases() -> None:
    report = (
        "pos\ttype\tscore\tconf\n"
        "4\t5' splice site\t7.0\t0.8\n"
        "9\t3' splice site\t7.0\t0.8\n"
    )
    kinds = {s.kind for s in parse_assp_report(report, seq_len=25)}
    assert kinds == {"donor", "acceptor"}


def test_parse_header_but_no_sites_is_empty() -> None:
    report = "position\tsite type\tscore\tconfidence\n\n(no sites found)\n"
    assert parse_assp_report(report, seq_len=25) == ()


def test_parse_no_header_raises() -> None:
    with pytest.raises(AsspReportError):
        parse_assp_report("ASSP is temporarily unavailable, please try later.", seq_len=25)


def test_parse_malformed_number_raises() -> None:
    report = "position\tsite type\tscore\tconfidence\n5\tdonor\tNOTANUMBER\t0.9\n"
    with pytest.raises(AsspReportError):
        parse_assp_report(report, seq_len=25)


def test_parse_unrecognized_site_type_raises() -> None:
    report = "position\tsite type\tscore\tconfidence\n5\tenhancer\t8.0\t0.9\n"
    with pytest.raises(AsspReportError):
        parse_assp_report(report, seq_len=25)


def test_parse_position_out_of_range_raises() -> None:
    report = "position\tsite type\tscore\tconfidence\n999\tdonor\t8.0\t0.9\n"
    with pytest.raises(AsspReportError):
        parse_assp_report(report, seq_len=25)


def test_parse_tab_row_with_empty_interior_field_stays_aligned() -> None:
    # An empty INTERIOR tab field must not shift columns left: the parser maps by
    # header index, so dropping it would misalign score/confidence/site-type.
    report = (
        "position\tflag\tsite type\tscore\tconfidence\n"
        "5\t\tdonor\t8.0\t0.95\n"  # empty 'flag' column
    )
    (site,) = parse_assp_report(report, seq_len=25)
    assert (site.position, site.kind) == (4, "donor")
    assert site.score == pytest.approx(8.0)
    assert site.confidence == pytest.approx(0.95)


def test_parse_skips_prose_rows() -> None:
    # A footer line with the right column count but a non-integer position column
    # is prose, silently skipped -- not a malformed site row.
    report = (
        "position\tsite type\tscore\tconfidence\n"
        "5\tdonor\t8.0\t0.95\n"
        "Total\tsites: 1\t-\t-\n"
    )
    (site,) = parse_assp_report(report, seq_len=25)
    assert site.position == 4


# --------------------------------------------------------------------------- #
# Transports
# --------------------------------------------------------------------------- #


def test_fixture_transport_reads_by_hash() -> None:
    t = FixtureAsspTransport(str(FIXTURES))
    assert t.available() is True
    body = t.fetch(SEQ_WITH_SITES)
    assert "constitutive donor" in body


def test_fixture_transport_missing_degrades() -> None:
    t = FixtureAsspTransport(str(FIXTURES))
    # A sequence with no committed fixture -> AsspUnavailableError (graceful).
    with pytest.raises(AsspUnavailableError):
        t.fetch("ACGTACGTACGTACGTACGT")


def test_fixture_transport_available_false_for_missing_dir() -> None:
    assert FixtureAsspTransport("/no/such/dir").available() is False


class _CountingTransport:
    """A fake transport that records how many times ``fetch`` ran, per sequence."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch(self, dna: str) -> str:
        self.calls.append(dna)
        return "position\tsite type\tscore\tconfidence\n1\tdonor\t8.0\t0.9\n"

    def available(self) -> bool:
        return True


def test_caching_transport_caches_by_hash() -> None:
    inner = _CountingTransport()
    cache = CachingAsspTransport(inner)
    s1, s2 = "ATGAAATTT", "ATGCCCGGG"
    cache.fetch(s1)
    cache.fetch(s1)  # cache hit -- inner not called again
    cache.fetch(s1.lower())  # same sequence, different casing -> same key -> hit
    cache.fetch(s2)  # different sequence -> miss
    assert inner.calls == [s1, s2]


class _FlakyTransport:
    """A fake transport that fails its first fetch, then succeeds (per sequence)."""

    def __init__(self) -> None:
        self.calls = 0

    def fetch(self, dna: str) -> str:
        self.calls += 1
        if self.calls == 1:
            raise AsspUnavailableError("transient outage")
        return "position\tsite type\tscore\tconfidence\n1\tdonor\t8.0\t0.9\n"

    def available(self) -> bool:
        return True


def test_caching_does_not_memoize_failure() -> None:
    # A failed fetch must NOT be cached, or a transient outage would become a
    # permanent "unavailable" for that sequence (defeating retry recovery).
    inner = _FlakyTransport()
    cache = CachingAsspTransport(inner)
    with pytest.raises(AsspUnavailableError):
        cache.fetch("ATGAAATTT")  # call 1 raises, nothing cached
    body = cache.fetch("ATGAAATTT")  # call 2 re-invokes inner and succeeds
    assert "donor" in body
    cache.fetch("ATGAAATTT")  # call 3 is a cache hit -- inner not called again
    assert inner.calls == 2


def _fake_httpx(script: list[str]) -> tuple[object, dict[str, int]]:
    """A minimal stand-in for the ``httpx`` module for the live-transport tests.

    ``script`` is a list of outcomes consumed one per ``post`` call: the sentinel
    ``"RAISE"`` raises the fake ``HTTPError`` (a retryable transient), any other
    string is returned as the response body.
    """
    calls = {"n": 0}

    class HTTPError(Exception):
        pass

    def post(url: str, data: object = None, timeout: object = None) -> object:
        i = calls["n"]
        calls["n"] += 1
        outcome = script[min(i, len(script) - 1)]
        if outcome == "RAISE":
            raise HTTPError("transient")
        return types.SimpleNamespace(text=outcome, raise_for_status=lambda: None)

    fake = types.SimpleNamespace(post=post, HTTPError=HTTPError)
    return fake, calls


def test_http_transport_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    fake, calls = _fake_httpx(["position\tsite type\tscore\tconfidence\n1\tdonor\t8\t0.9\n"])
    monkeypatch.setattr(assp_mod, "_import_httpx", lambda: fake)
    t = HttpAsspTransport(min_interval_s=0.0, _sleep=lambda _s: None)
    assert t.available() is True
    body = t.fetch("ATGAAATTT")
    assert "donor" in body
    assert calls["n"] == 1


def test_http_transport_retries_then_succeeds(monkeypatch: pytest.MonkeyPatch) -> None:
    fake, calls = _fake_httpx(["RAISE", "RAISE", "recovered body"])
    monkeypatch.setattr(assp_mod, "_import_httpx", lambda: fake)
    slept: list[float] = []
    t = HttpAsspTransport(
        min_interval_s=0.0, max_attempts=4, backoff_base_s=2.0, _sleep=slept.append
    )
    body = t.fetch("ATGAAATTT")
    assert body == "recovered body"
    assert calls["n"] == 3  # two transient failures, third try succeeds
    assert slept == [2.0, 4.0]  # exponential backoff between the retries


def test_http_transport_exhausts_and_degrades(monkeypatch: pytest.MonkeyPatch) -> None:
    fake, _calls = _fake_httpx(["RAISE"])
    monkeypatch.setattr(assp_mod, "_import_httpx", lambda: fake)
    t = HttpAsspTransport(min_interval_s=0.0, max_attempts=2, _sleep=lambda _s: None)
    with pytest.raises(AsspUnavailableError):
        t.fetch("ATGAAATTT")


def test_default_transport_selects_fixture_when_env_set(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(FIXTURE_DIR_ENV_VAR, str(FIXTURES))
    t = default_assp_transport()
    assert isinstance(t, CachingAsspTransport)
    assert isinstance(t.inner, FixtureAsspTransport)


def test_default_transport_is_live_without_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(FIXTURE_DIR_ENV_VAR, raising=False)
    t = default_assp_transport()
    assert isinstance(t, CachingAsspTransport)
    assert isinstance(t.inner, HttpAsspTransport)


def test_http_transport_available_reflects_httpx() -> None:
    result = HttpAsspTransport().available()
    assert isinstance(result, bool)
    if "httpx" not in sys.modules:
        try:
            import httpx  # type: ignore[import-not-found]  # noqa: F401
        except ImportError:
            assert result is False


def test_http_transport_without_httpx_degrades() -> None:
    t = HttpAsspTransport()
    if t.available():
        pytest.skip("httpx is installed; the missing-dependency path is not exercised")
    with pytest.raises(AsspUnavailableError):
        t.fetch("ATGAAATTT")


# --------------------------------------------------------------------------- #
# Retry / backoff / rate-limit helpers (pure, no httpx, no real sleeping)
# --------------------------------------------------------------------------- #


class _Transient(Exception):
    pass


def test_with_retries_succeeds_first_try() -> None:
    sleeps: list[float] = []
    out = _with_retries(
        lambda: "ok", attempts=4, backoff_base_s=2.0, sleep=sleeps.append, retryable=(_Transient,)
    )
    assert out == "ok"
    assert sleeps == []  # no retries, no backoff


def test_with_retries_backs_off_then_succeeds() -> None:
    sleeps: list[float] = []
    state = {"n": 0}

    def attempt() -> str:
        state["n"] += 1
        if state["n"] < 3:
            raise _Transient("boom")
        return "ok"

    out = _with_retries(
        attempt, attempts=4, backoff_base_s=2.0, sleep=sleeps.append, retryable=(_Transient,)
    )
    assert out == "ok"
    assert sleeps == [2.0, 4.0]  # exponential: 2*2^0, 2*2^1


def test_with_retries_exhausts_and_raises() -> None:
    sleeps: list[float] = []

    def always_fail() -> str:
        raise _Transient("down")

    with pytest.raises(AsspUnavailableError):
        _with_retries(
            always_fail,
            attempts=3,
            backoff_base_s=1.0,
            sleep=sleeps.append,
            retryable=(_Transient,),
        )
    assert sleeps == [1.0, 2.0]  # slept before tries 2 and 3, not after the last


def test_with_retries_rejects_bad_attempts() -> None:
    with pytest.raises(ValueError):
        _with_retries(
            lambda: "x", attempts=0, backoff_base_s=1.0, sleep=lambda _s: None, retryable=()
        )


def test_throttle_first_call_no_wait_then_rate_limits() -> None:
    sleeps: list[float] = []
    clock = iter([100.0, 100.0, 100.3, 100.3])
    state: dict[str, float] = {}
    _throttle(state, 1.0, monotonic=lambda: next(clock), sleep=sleeps.append)
    assert sleeps == []  # first request: nothing to wait for
    _throttle(state, 1.0, monotonic=lambda: next(clock), sleep=sleeps.append)
    assert sleeps == [pytest.approx(0.7)]  # 0.3s elapsed of a 1.0s interval -> wait 0.7s


def test_throttle_disabled_when_interval_nonpositive() -> None:
    sleeps: list[float] = []
    _throttle({}, 0.0, monotonic=lambda: 0.0, sleep=sleeps.append)
    assert sleeps == []


# --------------------------------------------------------------------------- #
# End-to-end scoring against committed fixtures
# --------------------------------------------------------------------------- #


def test_cache_key_stable_and_case_insensitive() -> None:
    import hashlib

    key = cache_key(SEQ_WITH_SITES)
    assert key == hashlib.sha256(SEQ_WITH_SITES.encode()).hexdigest()
    assert cache_key(SEQ_WITH_SITES.lower()) == key


def test_score_sequence_places_confidence_at_sites() -> None:
    model = _fixture_predictor()
    result = model.score_sequence(SEQ_WITH_SITES)
    assert result.model_name == "assp"
    assert result.calibrated is False
    assert len(result.donor) == len(SEQ_WITH_SITES)
    # Fixture: constitutive donor @ pos 5 (0-based 4), cryptic acceptor @ pos 12
    # (0-based 11), alternative donor @ pos 19 (0-based 18).
    assert result.donor[4] == pytest.approx(0.95)
    assert result.donor[18] == pytest.approx(0.90)
    assert result.acceptor[11] == pytest.approx(0.91)
    # Everywhere else is zero.
    assert sum(1 for x in result.donor if x > 0) == 2
    assert sum(1 for x in result.acceptor if x > 0) == 1


def test_sites_returns_classified_records() -> None:
    sites = _fixture_predictor().sites(SEQ_WITH_SITES)
    assert [s.site_class for s in sites] == ["constitutive", "cryptic", "alternative"]


def test_score_sequence_empty_result_is_all_zero() -> None:
    result = _fixture_predictor().score_sequence(SEQ_NO_SITES)
    assert not any(result.donor)
    assert not any(result.acceptor)
    assert pooled_risk(result) == 0.0


def test_delta_splicing_is_zero_for_identical() -> None:
    assert _fixture_predictor().delta_splicing(SEQ_WITH_SITES, SEQ_WITH_SITES) == 0.0


def test_delta_splicing_orientation_is_larger_is_better() -> None:
    # The load-bearing orientation: pooled_risk(reference) - pooled_risk(designed).
    model = _fixture_predictor()
    # designed = NO_SITES (less risk) vs reference = WITH_SITES -> POSITIVE (better).
    assert model.delta_splicing(SEQ_NO_SITES, SEQ_WITH_SITES) > 0.0
    # designed = WITH_SITES (more risk) vs reference = NO_SITES -> NEGATIVE (worse).
    assert model.delta_splicing(SEQ_WITH_SITES, SEQ_NO_SITES) < 0.0


def test_raw_predictor_raises_when_unavailable() -> None:
    # The RAW predictor is honest: it RAISES when the service is unavailable
    # (graceful degradation is the cross-check wrapper's job, not the backend's).
    model = AsspSplicePredictor(transport=FixtureAsspTransport("/no/such/dir"))
    with pytest.raises(AsspUnavailableError):
        model.score_sequence(SEQ_WITH_SITES)


def test_importing_module_does_not_load_httpx() -> None:
    # Guard against a regression that adds a top-level httpx import, which would
    # make `import bt4` heavier (CLAUDE.md section 3).
    code = (
        "import bt4.biomodels.splice, bt4.biomodels.splice.assp, sys;"
        "bad=[m for m in ('httpx',) if m in sys.modules];"
        "print(bad); assert not bad, bad"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
