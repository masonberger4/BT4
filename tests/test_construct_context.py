"""Construct context: the sequence around the CDS (CLAUDE.md §1, §5 #3, #7, #9).

BT4 used to optimize a coding sequence in a vacuum, so a defect formed *between*
the user's 5'UTR/backbone and the first codon was structurally unreachable. These
tests pin the three things that make context safe to add:

* it changes the answer where it should (a junction-formed site is now avoided),
* it changes nothing where it should not (no context => byte-identical output),
* and it never leaks (the one outbound path can never see the backbone).
"""

from __future__ import annotations

import inspect

import pytest

from bt4 import api
from bt4.constraints.rules import HomopolymerConstraint
from bt4.constraints.seeded import SeededConstraint, seed_constraints
from bt4.constraints.uorf import UorfConstraint
from bt4.domain import ConstructContext
from bt4.domain.genetic_code import translate

# --------------------------------------------------------------------------- #
# The value object.
# --------------------------------------------------------------------------- #


def test_flanks_are_truncated_at_the_nearest_unknown_base() -> None:
    """Only the contiguous KNOWN run touching the CDS is kept; nothing is invented."""
    context = ConstructContext(upstream="AACCNNGGTT", downstream="AACCNNGGTT")
    assert context.upstream == "GGTT"  # everything past the last N is dropped
    assert context.downstream == "AACC"  # everything past the first N is dropped


def test_context_rejects_degenerate_bases() -> None:
    # N means "unknown" and truncates; any other IUPAC code is not something a
    # constraint can evaluate, so it is refused rather than silently mishandled.
    with pytest.raises(ValueError, match="A,C,G,T"):
        ConstructContext(upstream="ACGRT")


def test_cds_offset_and_assembly() -> None:
    context = ConstructContext(upstream="GCCACC", downstream="TGATAA")
    assert context.cds_offset == 6
    assert context.upstream_tail(4) == "CACC"
    assert context.upstream_tail(0) == ""
    assert context.assemble("ATGAAA") == "GCCACCATGAAATGATAA"
    assert ConstructContext().is_empty


def test_masked_spans_are_validated_and_queryable() -> None:
    context = ConstructContext(upstream="ACGT", masked_spans=((0, 4),))
    assert context.is_masked(1, 3)
    assert not context.is_masked(2, 9)
    with pytest.raises(ValueError, match="masked span"):
        ConstructContext(masked_spans=((5, 2),))


def test_topology_must_be_known() -> None:
    with pytest.raises(ValueError, match="topology"):
        ConstructContext(topology="knotted")


# --------------------------------------------------------------------------- #
# Seeding: the junction becomes visible without changing the Constraint protocol.
# --------------------------------------------------------------------------- #


def test_seeded_constraint_sees_the_junction() -> None:
    inner = HomopolymerConstraint(6)
    seeded = SeededConstraint(inner, "GCCACCAAAA")
    # A lysine AAA codon would complete a 7-A run across the junction. The bare
    # rule cannot see it (its prefix is empty); the seeded one can.
    assert inner.ok_suffix("", "AAA") is True
    assert seeded.ok_suffix("", "AAA") is False
    # Once the prefix is longer than the rule's context, the seed cannot matter.
    assert seeded.ok_suffix("CCCCCCCC", "AAA") is inner.ok_suffix("CCCCCCCC", "AAA")


def test_seeded_validate_reports_in_cds_coordinates() -> None:
    seeded = SeededConstraint(HomopolymerConstraint(4), "TTTTAAAA")
    violations = list(seeded.validate("AAAGGG"))
    assert violations, "a junction-spanning run must be reported"
    for v in violations:
        assert 0 <= v.start <= v.end <= len("AAAGGG")  # indexes the delivered CDS
    assert any("junction" in v.detail for v in violations)


def test_seeded_validate_drops_violations_wholly_upstream() -> None:
    # A run entirely inside the leader is not part of the designed sequence and no
    # codon choice can fix it, so it must not be charged to the design.
    seeded = SeededConstraint(HomopolymerConstraint(4), "AAAAAAAA")
    assert list(seeded.validate("GCGCGCGC")) == []


def test_seed_constraints_is_a_no_op_without_upstream() -> None:
    constraints = [HomopolymerConstraint(6)]
    assert seed_constraints(constraints, "") is constraints


# --------------------------------------------------------------------------- #
# End to end: it changes the answer, and only where it should.
# --------------------------------------------------------------------------- #


def test_no_context_is_byte_identical(  # invariant #7
) -> None:
    config = api.OptimizeConfig(max_homopolymer=5)
    with_none = api.OptimizeConfig(max_homopolymer=5, context=None)
    first = api.optimize("MAALKHETQWY", config)
    second = api.optimize("MAALKHETQWY", with_none)
    assert first.dna == second.dna
    assert first.audit["manifest"] == second.audit["manifest"]


def test_context_avoids_a_motif_formed_at_the_junction() -> None:
    """The headline capability: a defect that only exists in the assembled construct.

    Human Val prefers GTG. With a leader ending ``CCA``, choosing GTG completes a
    forbidden ``CCAGTG`` that lies half in the leader and half in the CDS -- so it
    is invisible to a CDS-only optimizer and unavoidable without context.
    """
    context = ConstructContext(upstream="GGCACCA")
    common = {
        "forbidden_motifs": ("CCAGTG",),
        "avoid_reverse_complement": False,
        "max_homopolymer": None,
    }
    blind = api.optimize("VK", api.OptimizeConfig(**common))  # type: ignore[arg-type]
    aware = api.optimize("VK", api.OptimizeConfig(**common, context=context))  # type: ignore[arg-type]

    assert "CCAGTG" in context.assemble(blind.dna), "premise: the blind run forms it"
    assert "CCAGTG" not in context.assemble(aware.dna)
    assert blind.dna != aware.dna
    # The design is still a valid back-translation of the protein.
    assert translate(aware.dna) == "VK*"


def test_context_changes_the_manifest() -> None:
    # Invariant #9: two runs that used different flanking sequence must not stamp
    # the same provenance, whichever provenance policy is in force.
    plain = api.optimize("MKV", api.OptimizeConfig(max_homopolymer=None))
    with_context = api.optimize(
        "MKV",
        api.OptimizeConfig(
            max_homopolymer=None, context=ConstructContext(upstream="GCCACCA")
        ),
    )
    assert plain.audit["manifest"] != with_context.audit["manifest"]


def test_context_provenance_policy_is_the_users_choice() -> None:
    """'omit' records the shape; 'hash' records identity. Both are honest, and differ."""
    context = ConstructContext(upstream="GCCACCAAAT")
    omitted = api.optimize(
        "MKV", api.OptimizeConfig(max_homopolymer=None, context=context)
    )
    hashed = api.optimize(
        "MKV",
        api.OptimizeConfig(
            max_homopolymer=None, context=context, context_provenance="hash"
        ),
    )
    assert omitted.audit["manifest"] != hashed.audit["manifest"]

    # Two different backbones of the SAME shape are distinguishable only under
    # "hash" -- which is exactly the trade-off the policy exists to expose.
    other = ConstructContext(upstream="TTTTTTTTTT")
    omitted_other = api.optimize(
        "MKV", api.OptimizeConfig(max_homopolymer=None, context=other)
    )
    hashed_other = api.optimize(
        "MKV",
        api.OptimizeConfig(
            max_homopolymer=None, context=other, context_provenance="hash"
        ),
    )
    assert omitted.audit["manifest"] == omitted_other.audit["manifest"]
    assert hashed.audit["manifest"] != hashed_other.audit["manifest"]


# --------------------------------------------------------------------------- #
# uORF pairing across the junction (the cheapest, best-evidenced construct win).
# --------------------------------------------------------------------------- #


def test_leader_atg_reading_into_the_cds_is_detected() -> None:
    leader, cds = "CATGCC", "ATGAAGTAA"
    constraint = UorfConstraint(region_nt=100, cds_offset=len(leader))
    spans = [(v.start, v.end) for v in constraint.validate(leader + cds)]
    assert spans == [(1, 10)], "the leader's out-of-frame ATG pairs with a stop in the CDS"
    # A CDS-only scan cannot see it at all -- that is the gap being closed.
    assert list(UorfConstraint(region_nt=100).validate(cds)) == []


def test_cds_offset_zero_preserves_the_old_behaviour() -> None:
    cds = "ATGAATGAAATAAGGGCCC"
    assert list(UorfConstraint(region_nt=100).validate(cds)) == list(
        UorfConstraint(region_nt=100, cds_offset=0).validate(cds)
    )


def test_uorf_frame_follows_the_cds_not_index_zero() -> None:
    """The frame is (a - cds_offset) % 3; a leader must not invert it.

    MUTATION THAT MUST FAIL THIS: use ``a % 3`` instead of ``(a - cds_offset) % 3``
    in ``_uorf_spans``.
    """
    # One sequence, two leader lengths. The ATG at index 3 is followed by TAA at
    # index 6, so whether it opens a uORF depends purely on the frame -- which
    # depends purely on where the CDS starts.
    sequence = "GCCATGTAAGGGCCC"
    # cds_offset 6: (3 - 6) % 3 == 0 -> the main frame, so this is an ordinary Met
    # codon, not a uORF opener.
    in_frame = list(UorfConstraint(region_nt=100, cds_offset=6).validate(sequence))
    assert [v for v in in_frame if v.start == 3] == []
    # cds_offset 7: (3 - 7) % 3 == 2 -> out of frame, and TAA closes it.
    out_of_frame = list(UorfConstraint(region_nt=100, cds_offset=7).validate(sequence))
    assert [(v.start, v.end) for v in out_of_frame if v.start == 3] == [(3, 9)]


def test_uorf_is_enforced_end_to_end_with_a_leader() -> None:
    context = ConstructContext(upstream="CATGCC")
    result = api.optimize(
        "MK", api.OptimizeConfig(avoid_uorf=True, max_homopolymer=None, context=context)
    )
    # Whatever it managed, it must report honestly and round-trip.
    assert result.audit["uorf_enforced"] in {"clean", "partial"}
    assert translate(result.dna) == "MK*"


# --------------------------------------------------------------------------- #
# The privacy rider: the one outbound path can never see the backbone.
# --------------------------------------------------------------------------- #


def test_the_network_crosscheck_cannot_be_handed_a_construct() -> None:
    """ASSP is the only outbound path; a backbone must never reach it.

    The guarantee is structural rather than a runtime filter: the cross-check
    takes a bare coding sequence and has no context parameter at all, so there is
    no argument through which flanking sequence could be transmitted.

    MUTATION THAT MUST FAIL THIS: add a ``context`` parameter to
    ``api.splice_crosscheck``.
    """
    signature = inspect.signature(api.splice_crosscheck)
    leaky = [p for p in signature.parameters if "context" in p or "backbone" in p]
    assert not leaky, f"the outbound cross-check must not accept context: {leaky}"


# --------------------------------------------------------------------------- #
# The whole-construct audit: "is the thing I am about to build clean?"
# --------------------------------------------------------------------------- #


def test_construct_audit_reports_enzyme_uniqueness() -> None:
    """The genuinely new capability: is this enzyme still a single cutter?

    A site can be absent from the insert and present in the backbone, so
    "unique in the plasmid" is a question only a whole-construct audit can answer.
    """
    context = ConstructContext(upstream="TTTGAATTCAAACCCGGGTTT", downstream="AAACCCGGG")
    config = api.OptimizeConfig(
        restriction_enzymes=("EcoRI",), max_homopolymer=None, context=context
    )
    result = api.optimize("MKVAA", config)
    audit = api.audit_construct(result.dna, config)

    (ecori,) = audit.enzymes
    assert ecori.enzyme == "EcoRI"
    assert ecori.site == "GAATTC"
    # The optimizer kept the insert free of the site, so the backbone's single
    # occurrence survives and the enzyme still cuts exactly once.
    assert ecori.in_cds == 0
    assert ecori.count == 1
    assert ecori.unique
    assert audit.cds_start == len(context.upstream)
    assert audit.cds_end == audit.cds_start + len(result.dna)


def test_construct_audit_attributes_a_junction_finding() -> None:
    context = ConstructContext(upstream="CCCCAAAA")
    config = api.OptimizeConfig(max_homopolymer=6, context=context)
    audit = api.audit_construct("AAAGGGCCC", config)
    assert audit.junction_violations == 1
    assert any("[junction]" in v.detail for v in audit.violations)
    assert not audit.is_clean


def test_construct_audit_masks_repeats_by_construction() -> None:
    """AAV ITRs / LVV LTRs are repeats by design; unmasked they drown the report."""
    itr = "GCGCGCTCGCTCGCTCGCTGGCTCGCTCGCTCGCTC"
    context = ConstructContext(upstream=itr, masked_spans=((0, len(itr)),))
    config = api.OptimizeConfig(
        max_repeat_length=6, max_homopolymer=None, context=context
    )
    audit = api.audit_construct("ATGAAGGTGTAA", config)
    assert audit.masked_violations > 0  # they were found...
    assert audit.violations == ()  # ...and excluded, not silently ignored
    assert audit.is_clean


def test_construct_audit_without_context_audits_the_bare_cds() -> None:
    config = api.OptimizeConfig(max_homopolymer=4)
    audit = api.audit_construct("AAAAAAGGG", config)
    assert audit.cds_start == 0
    assert audit.construct == "AAAAAAGGG"
    assert audit.junction_violations == 0
    assert any(v.constraint == "homopolymer" for v in audit.violations)


# --------------------------------------------------------------------------- #
# Junction folding: the region reported is the region optimized.
# --------------------------------------------------------------------------- #


def test_folding_region_is_shared_between_objective_and_audit() -> None:
    """One window function for both, so a reported/optimized mismatch cannot recur."""
    from bt4.biomodels.folding import default as folding_default
    from bt4.biomodels.folding import junction_window

    protein = "MAALKHETQWYCDEFGHIKLM"
    model = folding_default()

    plain = api.optimize(
        protein,
        api.OptimizeConfig(refine=True, refine_iterations=120, max_homopolymer=None),
    )
    assert plain.audit["folding_dg"] == model.five_prime_dg(
        junction_window("", plain.dna), None
    )
    assert plain.audit["folding_spans_junction"] is False

    context = ConstructContext(upstream="GGGCCCAAA" * 6)
    aware = api.optimize(
        protein,
        api.OptimizeConfig(
            refine=True, refine_iterations=120, max_homopolymer=None, context=context
        ),
    )
    assert aware.audit["folding_dg"] == model.five_prime_dg(
        junction_window(context.upstream, aware.dna), None
    )
    # With a leader the folded region really does span the junction, so it is
    # longer than the CDS-only window.
    assert aware.audit["folding_spans_junction"] is True
    assert aware.audit["folding_region_nt"] > plain.audit["folding_region_nt"]
