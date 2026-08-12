"""Tests for the REBASE-derived restriction-enzyme catalog.

The catalog used to be seventeen hand-typed pairs described only as
"textbook-correct". It is now derived from a version-pinned REBASE release by
``scripts/build_enzyme_catalog.py``, so these tests check the two things that
change makes checkable:

* **Provenance** -- the shipped bytes are content-hashed and carry a complete
  re-derivation trail (REBASE version, URL, source digest, selection rule), the
  same standard the recounted codon tables meet (CLAUDE.md §8).
* **No regression and no silent loss** -- every previously shipped enzyme must
  still resolve to the *same* site, and the Type IIS enzymes modern cloning
  depends on must be present (an earlier build silently dropped all of them by
  mishandling REBASE's two-strand listing).

Plus the governing principle (CLAUDE.md §1): a requested site is either enforced
or the run is honestly refused -- never silently ignored.
"""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from bt4 import api
from bt4.constraints.iupac import find_iupac, is_iupac, reverse_complement_iupac
from bt4.constraints.restriction import (
    ENZYMES,
    RestrictionSiteConstraint,
    available_enzymes,
    enzyme_provenance,
    resolve_enzyme,
)
from bt4.domain.sequence import validate_dna

# The exact catalog BT4 shipped before it was REBASE-derived. Every one of these
# must survive the switch with an identical site: this is what proves the new
# catalog agrees with the old hand-checked values rather than quietly redefining
# a site somebody's protocol depends on.
LEGACY_CATALOG: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG",
    "NdeI": "CATATG",
    "NcoI": "CCATGG",
    "EcoRV": "GATATC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "HinfI": "GANTC",
    "DraIII": "CACNNNGTG",
}

# Type IIS enzymes: the Golden Gate / MoClo workhorses. REBASE lists their
# asymmetric site once per strand, and a build that treats that as "two different
# sites" drops every one of them -- which is exactly what happened once.
TYPE_IIS: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
    "Esp3I": "CGTCTC",
    "AarI": "CACCTGC",
}


def test_catalog_is_substantial_and_sorted() -> None:
    names = available_enzymes()
    assert len(names) > 400, "the derived catalog should cover the commercial set"
    assert list(names) == sorted(names)
    assert len(set(names)) == len(names)


def test_every_site_is_valid_iupac() -> None:
    for name, site in ENZYMES.items():
        assert site, f"{name} has an empty site"
        assert is_iupac(site), f"{name} site {site!r} is not IUPAC"
        assert site == site.upper()
        assert 4 <= len(site) <= 12, f"{name} site {site!r} has implausible length"


def test_catalog_is_read_only() -> None:
    # Shipped data, not a mutable registry: a caller must not be able to redefine
    # a recognition site at runtime and silently invalidate the provenance stamp.
    with pytest.raises(TypeError):
        ENZYMES["EcoRI"] = "AAAAAA"  # type: ignore[index]


# --------------------------------------------------------------------------- #
# No regression, no silent loss.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("name", "site"), sorted(LEGACY_CATALOG.items()))
def test_legacy_enzymes_survive_unchanged(name: str, site: str) -> None:
    assert name in ENZYMES, f"{name} disappeared from the catalog"
    assert ENZYMES[name] == site, f"{name} site changed"


@pytest.mark.parametrize(("name", "site"), sorted(TYPE_IIS.items()))
def test_type_iis_enzymes_are_present(name: str, site: str) -> None:
    """REBASE's two-strand listing must not be mistaken for two sites."""
    assert name in ENZYMES, f"{name} missing: the asymmetric-site parse regressed"
    assert ENZYMES[name] == site


def test_isoschizomers_agree_on_their_shared_site() -> None:
    # Users name the enzyme they own; different names for one site must not
    # disagree about what that site is.
    assert ENZYMES["Acc65I"] == ENZYMES["KpnI"] == "GGTACC"
    assert ENZYMES["Esp3I"] == ENZYMES["BsmBI"] == "CGTCTC"


# --------------------------------------------------------------------------- #
# Name resolution.
# --------------------------------------------------------------------------- #


def test_resolve_is_case_insensitive() -> None:
    assert resolve_enzyme("EcoRI") == "GAATTC"
    assert resolve_enzyme("ecori") == "GAATTC"
    assert resolve_enzyme("  BsaI  ") == "GGTCTC"


def test_unknown_enzyme_suggests_instead_of_dumping_the_catalog() -> None:
    """A 500-name wall of text hides the answer; near misses give it."""
    with pytest.raises(ValueError) as excinfo:
        resolve_enzyme("EcoR1")
    message = str(excinfo.value)
    assert "EcoRI" in message, "the obvious correction should be offered"
    # The whole catalog must not be pasted into the error.
    assert message.count(",") < 10
    assert len(message) < 300


def test_unknown_enzyme_without_near_miss_still_explains() -> None:
    with pytest.raises(ValueError) as excinfo:
        resolve_enzyme("zzzzzzzz")
    assert "unknown enzyme" in str(excinfo.value)
    assert str(len(ENZYMES)) in str(excinfo.value)


# --------------------------------------------------------------------------- #
# Provenance.
# --------------------------------------------------------------------------- #


def test_provenance_hash_matches_the_shipped_catalog() -> None:
    raw = files("bt4.constraints.data").joinpath("rebase_enzymes.tsv").read_bytes()
    import hashlib

    assert enzyme_provenance()["sha256"] == hashlib.sha256(raw).hexdigest()


def test_provenance_carries_the_rederivation_trail() -> None:
    prov = enzyme_provenance()
    for key in (
        "source",
        "source_url",
        "source_sha256",
        "rebase_version",
        "build",
        "enzyme_count",
        "selection",
        "rebuild_command",
        "note",
        "license_note",
    ):
        assert key in prov, f"provenance is missing {key!r}"
    assert len(str(prov["source_sha256"])) == 64
    assert prov["enzyme_count"] == len(ENZYMES)
    # REBASE's terms are citation-gated, not a CC/public-domain grant, and the
    # sidecar must say so rather than implying a permissive licence.
    licence = str(prov["license_note"])
    assert "Roberts" in licence
    assert "not a" in licence.lower() or "NOT a" in licence
    # The honest scope limit: BT4 models the site, not the digest.
    assert "cut position" in str(prov["note"])


def test_selection_tally_is_monotonic() -> None:
    selection = enzyme_provenance()["selection"]
    assert isinstance(selection, dict)
    # Each filter can only narrow the previous stage.
    assert (
        selection["records_in_source"]
        >= selection["type_ii"]
        >= selection["commercially_available"]
        >= selection["with_usable_site"]
    )
    assert selection["shipped"] == len(ENZYMES)


def test_provenance_is_json_serializable() -> None:
    json.dumps(enzyme_provenance())


# --------------------------------------------------------------------------- #
# The governing principle: enforced, or honestly refused. Never ignored.
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["EcoRI", "BsaI", "SapI", "BslI", "MspJI"])
def test_requested_site_is_enforced_or_the_run_is_refused(name: str) -> None:
    """No third option: a delivered sequence never contains a banned site.

    ``BslI`` (``CCNNNNNNNGG``) and ``MspJI`` (``CNNR``) are real, highly
    degenerate REBASE sites. Banning a near-wildcard inside a coding sequence may
    be genuinely unsatisfiable -- in which case BT4 must say so, not hand back a
    sequence that still contains it (CLAUDE.md §1, §10).
    """
    site = resolve_enzyme(name)
    try:
        result = api.optimize("MAALKHETQW", api.OptimizeConfig(restriction_enzymes=(name,)))
    except api.InfeasibleError as exc:
        assert "restriction_site" in exc.constraints
        return
    assert list(find_iupac(result.dna, site)) == []
    # Sites are double-stranded, so the reverse complement is banned as well.
    assert list(find_iupac(result.dna, reverse_complement_iupac(site))) == []
    assert result.metrics.hard_violations == 0


def test_constraint_bans_both_strands_of_an_asymmetric_site() -> None:
    """An asymmetric Type IIS site must be banned in both orientations."""
    constraint = RestrictionSiteConstraint(enzymes=("BsaI",))
    forward = validate_dna("AAA" + "GGTCTC" + "AAA")
    reverse = validate_dna("AAA" + reverse_complement_iupac("GGTCTC") + "AAA")
    assert list(constraint.validate(forward)), "forward strand site not flagged"
    assert list(constraint.validate(reverse)), "reverse strand site not flagged"


def test_api_exposes_the_catalog_and_its_provenance() -> None:
    """Frontends reach the catalog through ``bt4.api`` alone (layering §3)."""
    assert api.available_enzymes() == available_enzymes()
    assert api.resolve_enzyme("bsai") == "GGTCTC"
    assert api.enzyme_provenance()["enzyme_count"] == len(ENZYMES)
