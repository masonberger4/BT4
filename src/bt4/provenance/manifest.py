"""Deterministic provenance — a result must be reproducible from its stamp alone.

BT3 hashed config *field names*, so a swapped custom codon table produced the
same ``config_hash`` — a false claim of identical provenance. BT4 hashes the
actual *contents* of every input that influences a result (config, codon/model
data, seed, code version), so any change that could change the output changes
the manifest.

Pure stdlib (``hashlib`` + ``json``). No timestamps by default — determinism is
the point.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

__all__ = [
    "Manifest",
    "build_manifest",
    "config_hash",
    "content_hash",
    "resolve_git_commit",
]


def content_hash(data: bytes | str) -> str:
    """Return the hex SHA-256 of ``data`` (UTF-8 encoded if a ``str``)."""
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


@lru_cache(maxsize=1)
def resolve_git_commit() -> str | None:
    """Return the current source commit SHA, or ``None`` when unavailable.

    Resolves ``git rev-parse HEAD`` in the directory containing this package.
    Returns ``None`` when git is absent or the source is not a checkout (e.g. an
    installed wheel), so provenance degrades gracefully rather than failing. The
    result is constant within a checkout, so it does not break determinism
    (invariant #7): two runs on the same commit stamp identically, while runs
    from different commits differ (invariant #9, "plus git SHA"). Cached so the
    subprocess runs at most once per process.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    sha = out.stdout.strip()
    return sha if out.returncode == 0 and sha else None


def _canonical_json(obj: object) -> str:
    """Serialize ``obj`` to canonical JSON (sorted keys, stable separators)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def config_hash(config: Mapping[str, object]) -> str:
    """Return a stable hex SHA-256 over the canonical JSON of ``config``.

    Two configs that serialize to the same canonical JSON hash identically;
    any difference in a value (not just a key name) changes the hash.
    """
    return content_hash(_canonical_json(config))


@dataclass(frozen=True, slots=True)
class Manifest:
    """A deterministic, content-addressed record of everything behind a result.

    Attributes:
        bt4_version: The single-sourced BT4 version string.
        config_hash: Hash over the run configuration's *values*.
        inputs: Map of logical input name -> content hash (codon tables, model
            weights, custom data files). Hashing contents, not names, is what
            makes a swapped table produce a different manifest.
        seed: The master RNG seed threaded through the run (``None`` if unseeded).
        git_commit: The source commit the run was produced from, when known.
        extra: Any additional deterministic key/value provenance.
    """

    bt4_version: str
    config_hash: str
    inputs: Mapping[str, str] = field(default_factory=dict)
    seed: int | None = None
    git_commit: str | None = None
    extra: Mapping[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-ready dict of this manifest."""
        return {
            "bt4_version": self.bt4_version,
            "config_hash": self.config_hash,
            "inputs": dict(sorted(self.inputs.items())),
            "seed": self.seed,
            "git_commit": self.git_commit,
            "extra": dict(sorted(self.extra.items())),
        }

    @property
    def stamp(self) -> str:
        """A single hex digest identifying this exact manifest.

        Reproducing a result means reproducing this stamp: same inputs, config,
        seed, and code version ⇒ same stamp.
        """
        return content_hash(_canonical_json(self.to_dict()))


def build_manifest(
    *,
    bt4_version: str,
    config: Mapping[str, object],
    inputs: Mapping[str, str] | None = None,
    seed: int | None = None,
    git_commit: str | None = None,
    extra: Mapping[str, str] | None = None,
) -> Manifest:
    """Assemble a :class:`Manifest`, hashing ``config`` into its ``config_hash``.

    Args:
        bt4_version: The single-sourced version string (``bt4.__version__``).
        config: The run configuration; hashed by value.
        inputs: Optional map of input name -> precomputed content hash.
        seed: Optional master RNG seed.
        git_commit: Optional source commit hash.
        extra: Optional additional deterministic provenance.

    Returns:
        A fully-populated, deterministic manifest.
    """
    return Manifest(
        bt4_version=bt4_version,
        config_hash=config_hash(config),
        inputs=dict(inputs or {}),
        seed=seed,
        git_commit=git_commit,
        extra=dict(extra or {}),
    )
