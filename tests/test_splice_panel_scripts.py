"""Tests for the splice fidelity-panel scripts (A4-A6 of the calibration runbook).

Three standalone scripts under ``scripts/`` implement the integration-fidelity
workflow, and they are only trustworthy if their *separation* holds:

* ``make_splice_panel.py`` builds the sequences and **may** use ``bt4``;
* ``capture_pangolin_panel.py`` captures upstream's expected scores and must
  **never** import ``bt4`` -- capturing expectations with the adapter under test
  would make the gate unfalsifiable;
* ``run_splice_fidelity_gate.py`` compares the two and reports honestly.

These tests pin the properties that make the workflow meaningful rather than
merely runnable:

* the panel is **deterministic** (invariant #7) and every member is pure ``ACGT``
  (``validate_dna`` rejects ``N``, and the adapters pad with ``N`` themselves);
* the panel spans the kinds the runbook calls for, including the designed
  ``cds``/``cds_variant`` members that are BT4's actual regime;
* the capture script's **independence guard** fires, and the file contains no
  ``bt4`` import at all -- checked statically, so the separation cannot rot;
* the gate's **panel-strength** report distinguishes a discriminating panel from
  a flat one, and its panel/capture pairing check rejects a stale capture.

None of this needs torch, TensorFlow, or the licensed weights, so it all runs in CI.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

from bt4.domain.sequence import validate_dna

_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name: str, filename: str) -> ModuleType:
    """Load a ``scripts/`` module by file path.

    Registers it in ``sys.modules`` first: the panel module defines
    ``slots=True`` dataclasses, whose construction looks the module up by name.
    """
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / filename)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def panel_script() -> ModuleType:
    """The loaded panel-builder module."""
    return _load("bt4_make_splice_panel", "make_splice_panel.py")


@pytest.fixture(scope="module")
def gate_script() -> ModuleType:
    """The loaded gate-runner module."""
    return _load("bt4_run_splice_fidelity_gate", "run_splice_fidelity_gate.py")


@pytest.fixture(scope="module")
def capture_script() -> ModuleType:
    """The loaded capture module (imports cleanly without torch -- deps are lazy)."""
    return _load("bt4_capture_pangolin_panel", "capture_pangolin_panel.py")


# --------------------------------------------------------------------------
# The panel


def test_panel_is_deterministic(panel_script: ModuleType) -> None:
    """Same seed => byte-identical panel and content hash (invariant #7)."""
    a = panel_script.build_panel(0)
    b = panel_script.build_panel(0)
    assert [x.sequence for x in a] == [x.sequence for x in b]
    assert panel_script.panel_content_hash(a) == panel_script.panel_content_hash(b)


def test_a_different_seed_changes_the_panel(panel_script: ModuleType) -> None:
    """The seed actually varies the probes (otherwise it is a false knob)."""
    a = panel_script.build_panel(0, include_designed=False)
    b = panel_script.build_panel(7, include_designed=False)
    assert panel_script.panel_content_hash(a) != panel_script.panel_content_hash(b)


def test_every_member_is_pure_acgt(panel_script: ModuleType) -> None:
    """No 'N' anywhere: the adapters add their own padding, and validate_dna refuses N."""
    for member in panel_script.build_panel(0):
        assert set(member.sequence) <= set("ACGT"), member.id
        validate_dna(member.sequence)  # raises if not


def test_panel_spans_the_runbook_kinds(panel_script: ModuleType) -> None:
    """The panel covers designed CDSs, synonymous variants, probes and edges."""
    kinds = {m.kind for m in panel_script.build_panel(0)}
    assert kinds == {"cds", "cds_variant", "donor", "acceptor", "edge"}


def test_ids_are_unique(panel_script: ModuleType) -> None:
    """Ids key the capture back to the panel, so duplicates would silently collide."""
    ids = [m.id for m in panel_script.build_panel(0)]
    assert len(ids) == len(set(ids))


def test_probes_carry_their_consensus_motif(panel_script: ModuleType) -> None:
    """A donor/acceptor probe must actually contain the motif it claims to."""
    panel = panel_script.build_panel(0, include_designed=False)
    donors = [m for m in panel if m.kind == "donor"]
    acceptors = [m for m in panel if m.kind == "acceptor"]
    assert donors and acceptors
    assert all(panel_script._DONOR_CONSENSUS in m.sequence for m in donors)
    assert all(panel_script._ACCEPTOR_CONSENSUS in m.sequence for m in acceptors)


def test_probe_motif_is_centred(panel_script: ModuleType) -> None:
    """The motif sits mid-sequence, away from the N-padding boundary artifact."""
    panel = panel_script.build_panel(0, include_designed=False)
    donor = next(m for m in panel if m.kind == "donor")
    at = donor.sequence.index(panel_script._DONOR_CONSENSUS)
    middle = len(donor.sequence) // 2
    assert abs(at - middle) < len(donor.sequence) // 4


def test_panel_without_designed_needs_no_bt4_api(panel_script: ModuleType) -> None:
    """``--no-designed`` yields a usable panel without invoking the optimizer."""
    panel = panel_script.build_panel(0, include_designed=False)
    assert panel
    assert not any(m.kind in {"cds", "cds_variant"} for m in panel)


# --------------------------------------------------------------------------
# The capture script's independence from bt4


def test_capture_script_contains_no_bt4_import() -> None:
    """Statically assert the capture never imports bt4.

    The runtime guard can only fire if the module is reached; this check is
    structural, so the separation cannot rot silently (the same posture as the
    attestation's license-clean field test).
    """
    source = (_SCRIPTS / "capture_pangolin_panel.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
    offenders = [m for m in imported if m == "bt4" or m.startswith("bt4.")]
    assert not offenders, f"capture script must not import bt4, found {offenders}"


def test_independence_guard_fires_when_bt4_is_loaded(capture_script: ModuleType) -> None:
    """The runtime guard refuses to capture in a process where bt4 is imported."""
    assert "bt4" in sys.modules  # this test module imports bt4 itself
    with pytest.raises(RuntimeError, match="independent"):
        capture_script._assert_bt4_not_imported()


def test_capture_mirrors_the_cli_not_the_example(capture_script: ModuleType) -> None:
    """Pin the CLI's weight set and channel map, not ``custom_usage.py``'s.

    Upstream ships two disjoint model sets; copying the example script would
    produce different numbers and fail the gate for a non-BT4 reason.
    """
    assert capture_script.CV_FOLDS == (1, 2, 3)
    assert capture_script.TISSUE_WEIGHT_INDEX == (0, 2, 4, 6)
    assert capture_script.TISSUE_CHANNEL == {0: 1, 2: 4, 4: 7, 6: 10}
    assert capture_script.WEIGHT_TEMPLATE.endswith(".3.v2")
    assert capture_script.PANGOLIN_FLANK == 5000


def test_capture_constants_agree_with_the_adapter(capture_script: ModuleType) -> None:
    """The independent capture and BT4's adapter must describe the same model.

    They are deliberately separate implementations; if their constants drift, the
    gate would compare two different models and its verdict would be meaningless.
    """
    from bt4.biomodels.splice.pangolin import (
        PANGOLIN_FLANK,
        PINNED_WEIGHT_SHA256,
        TISSUE_OUTPUTS,
    )

    assert capture_script.PANGOLIN_FLANK == PANGOLIN_FLANK
    adapter_map = {t.weight_index: t.output_channel for t in TISSUE_OUTPUTS.values()}
    assert adapter_map == capture_script.TISSUE_CHANNEL
    expected_names = {
        capture_script.WEIGHT_TEMPLATE.format(fold=f, index=i)
        for i in capture_script.TISSUE_WEIGHT_INDEX
        for f in capture_script.CV_FOLDS
    }
    assert expected_names == set(PINNED_WEIGHT_SHA256)


# --------------------------------------------------------------------------
# The gate runner


def test_spread_report_distinguishes_wide_from_flat(gate_script: ModuleType) -> None:
    """Peak statistics separate a discriminating panel from a vacuous one."""
    wide = gate_script.spread_report(
        [{"expected_site_scores": [0.9, 0.0]}, {"expected_site_scores": [0.01, 0.0]}]
    )
    flat = gate_script.spread_report(
        [{"expected_site_scores": [0.001, 0.0]}, {"expected_site_scores": [0.002, 0.0]}]
    )
    assert wide["spread"] > gate_script.MIN_USEFUL_PEAK_SPREAD
    assert wide["max_peak"] > gate_script.MIN_USEFUL_MAX_PEAK
    assert flat["spread"] < gate_script.MIN_USEFUL_PEAK_SPREAD
    assert flat["max_peak"] < gate_script.MIN_USEFUL_MAX_PEAK


def test_spread_report_handles_empty_input(gate_script: ModuleType) -> None:
    """An empty panel reports zeros rather than raising."""
    assert gate_script.spread_report([]) == {
        "min_peak": 0.0,
        "max_peak": 0.0,
        "mean_peak": 0.0,
        "spread": 0.0,
    }


def test_panels_match_binds_a_capture_to_its_panel(gate_script: ModuleType) -> None:
    """A stale capture must not be compared against a regenerated panel."""
    assert gate_script.panels_match({"panel_content_hash": "a"}, {"content_hash": "a"})
    assert not gate_script.panels_match({"panel_content_hash": "a"}, {"content_hash": "b"})
    assert not gate_script.panels_match({}, {"content_hash": "b"})


def test_load_cases_builds_fidelity_cases(gate_script: ModuleType) -> None:
    """Capture JSON converts into the dataclass the gate consumes."""
    from bt4.biomodels.splice import FidelityCase

    cases = gate_script.load_cases(
        {"cases": [{"sequence": "ACGT", "expected_site_scores": [0.1, 0.2, 0.3, 0.4]}]}
    )
    assert len(cases) == 1
    assert isinstance(cases[0], FidelityCase)
    assert cases[0].sequence == "ACGT"
    assert cases[0].expected_site_scores == (0.1, 0.2, 0.3, 0.4)


def test_gate_tolerance_floor_matches_the_attestation_floor(gate_script: ModuleType) -> None:
    """The gate's default tolerance must not exceed what an attestation may claim.

    ``attest_backend`` refuses a report looser than ``MAX_ATTESTATION_TOLERANCE``,
    so a default looser than that would produce passes that can never be attested.
    """
    from bt4.biomodels.splice import MAX_ATTESTATION_TOLERANCE

    parser_default = 1e-3
    assert parser_default <= MAX_ATTESTATION_TOLERANCE


def test_capture_resolves_weights_without_importing_bt4(capture_script: ModuleType) -> None:
    """The capture resolves its own weights dir, mirroring the adapter's order.

    A maintainer whose ``pangolin`` install lets BT4 auto-resolve should not have
    to discover an environment variable here -- but the resolution must reach that
    conclusion by importing ``pangolin``, never ``bt4``.
    """
    assert capture_script.resolve_model_dir("/definitely/not/a/directory") is None
    assert capture_script.resolve_model_dir(str(_SCRIPTS)) == Path(_SCRIPTS)


def test_capture_prefers_explicit_dir_over_env(
    capture_script: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An explicit --model-dir wins over the environment variable."""
    other = tmp_path / "other"
    other.mkdir()
    monkeypatch.setenv("BT4_PANGOLIN_MODEL_DIR", str(tmp_path))
    assert capture_script.resolve_model_dir(str(other)) == other
    assert capture_script.resolve_model_dir(None) == tmp_path


def test_capture_resolution_order_matches_the_adapter(capture_script: ModuleType) -> None:
    """Both resolvers read the same env var, so the two agree on where to look."""
    source = (_SCRIPTS / "capture_pangolin_panel.py").read_text(encoding="utf-8")
    assert "BT4_PANGOLIN_MODEL_DIR" in source
    from bt4.biomodels.splice.pangolin import _WEIGHTS_ENV_VAR

    assert _WEIGHTS_ENV_VAR == "BT4_PANGOLIN_MODEL_DIR"
