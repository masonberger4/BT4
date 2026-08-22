"""Shared pytest fixtures.

Currently one job: keep the process-wide splice-promotion switch from leaking
between tests.
"""

from __future__ import annotations

import os
from collections.abc import Iterator

import pytest

# Environment variables that a *code path under test* may set on ``os.environ``
# directly, rather than through ``monkeypatch``. Anything listed here is snapshotted
# before each test and restored after it.
_PROCESS_WIDE_SWITCHES = ("BT4_SPLICE_USE_ATTESTED",)


@pytest.fixture(autouse=True)
def _isolate_process_wide_switches() -> Iterator[None]:
    """Restore process-wide opt-in switches after every test.

    ``bt4.cli.__main__._enable_attested_splice`` implements ``--use-attested-splice``
    by assigning ``os.environ[BT4_SPLICE_USE_ATTESTED] = "1"``. That is correct for
    the CLI -- one switch has to govern every path (api, Studio, library) and a
    process-wide variable is how that is done -- but it means a test that drives
    ``main()`` with the flag mutates the environment of every test that runs after
    it, and ``monkeypatch.delenv`` does not undo an assignment it never made.

    Measured before this fixture existed: ``BT4_SPLICE_USE_ATTESTED`` escaped
    ``test_cli.py::test_use_attested_splice_flag_is_wired`` and stayed set for the
    remainder of the session, so every later test ran with attested promotion
    enabled. On CI that was invisible -- no CNN weights are installed, so the audit
    path's ``promote_if_attested`` is a no-op and the assertions still passed. On a
    machine holding the licensed weights it *would* promote Pangolin to
    ``calibrated=True`` across the suite -- the attestation ships and
    ``promote_if_attested`` honours the variable -- turning the ``all_calibrated is
    False`` assertions into environment-dependent tests. That consequence is inferred
    from the code path, not observed: it has not been run on a machine holding the
    weights. It is precisely the "promotion leaked into the default path" failure the
    opt-in design exists to prevent.

    Autouse and unconditional: the guarantee is worth more than the microseconds,
    and an opt-in isolation fixture only protects the tests that remember to ask.
    """
    saved = {name: os.environ.get(name) for name in _PROCESS_WIDE_SWITCHES}
    try:
        yield
    finally:
        for name, value in saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
