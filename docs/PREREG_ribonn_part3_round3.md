# Pre-registration — RiboNN Part 3, round 3

**Status: DRAFT until merged. Nothing has been run against this design.** Merging freezes
it. After that no threshold, protein, UTR, seed, variant recipe or analysis choice may
change, and a failure is recorded as a failure.

## 0. Where this sits

| round | outcome |
|---|---|
| 1 | **Stop.** Failed its own gate: median \|ρ(GC3)\| 0.7385 vs a pre-registered ≤ 0.70. |
| 2 | **Void.** Failed its sanity floor (`within_over_between` 0.180 < 0.20). See [`RESULT_ribonn_part3_round2.md`](RESULT_ribonn_part3_round2.md). |

Round 1's verdict stands and is reported alongside every later round. Round 2 answered
nothing — its variant set was too narrow to exercise the axis under test.

**Round 3 changes exactly two things**, both diagnosed from round 2's own data rather
than argued for:

1. **The variant recipe**, which caused the void.
2. **The adequacy check**, which caught the void but could not say *why* — because it was
   computed from the model's own output.

The gate itself is **unchanged**. It was validated against known-answer regimes
(`scripts/prereg_round2_selftest.py`) and nothing about round 2 impugns it. Reusing it
verbatim is deliberate: a gate that survives a void run unmodified is not a gate tuned to
an outcome.

## 1. What this round can and cannot establish

Unchanged from round 2, and restated because it is the easiest thing to lose:

**No measured expression data is used.** Every number is a RiboNN *prediction*. A moved
prediction is a statement about the model's output, **not** evidence that translation
efficiency changed. Whether RiboNN's ordering is *correct* is unanswerable here by
construction.

**"The score moved" is not a result** — RiboNN is deterministic, so a non-zero response is
guaranteed before the run starts. No endpoint below has that null.

The answerable question:

> For synonymous CDS variants of one protein under a fixed UTR — BT4's deployment regime —
> does RiboNN's score carry information that is **(a)** not available from the features BT4
> already optimizes for free, and **(b)** a stable property of the CDS rather than of the
> UTR it was scored in?

A pass is a **necessary condition** and a licence to look for a panel. It is not evidence
of skill, promotes nothing, and `calibrated` stays `False`.

## 2. Design

**Proteins — 16, freshly drawn.** Same rule as round 2 (MANE Select 1.5, 150–1200 aa,
round 1's genes excluded) with **`seed=20260822`**, via
`scripts/prereg_round3_draw_panel.py`. *Why redraw rather than reuse:* round 2's
per-protein gate values have been seen. The variants would differ, so reuse is defensible
— but a fresh draw removes the question for the cost of one script run.

> **Why a new draw script rather than round 2's.** Round 2's version validated the UTR-span
> constraint against the **full list of candidates it had collected** and then returned
> `contexts[:n]` — so truncation could discard the very context that satisfied the span,
> and the function reported success anyway. It did not bite in round 2 (that panel spans
> 2.300×, compliant), but it bit immediately here: the first round-3 draw returned 5′UTRs
> of **50 / 84 / 69 / 63 — a 1.68× span against a required ≥ 2×** — while exiting 0.
> `prereg_round3_draw_panel.py` validates *exactly what it returns*, and when the set is
> full but the span still fails it swaps a candidate into whichever **non-anchor** slot
> most improves the span, keeping it only if it does. Round 2's script is left **untouched
> at its recorded hash**: it is the code round 2 actually ran, and repairing it in place
> would destroy that run's verifiability to fix a run that is already void.

**Variants — 42 per protein, a temperature ladder.** Six rungs,
**T ∈ {0.4, 0.7, 1.0, 1.6, 2.6, 4.2}**, seven draws each, via `api.library` at the seed
recorded in the panel manifest.

The ladder is the round-2 fix, and it is chosen on measurement. Temperature does not
*widen* a sample — each rung spans only ~0.02–0.05 in GC — it **moves** it (ZNF286A
centres at GC 0.53 at T = 0.4 and 0.46 at T = 4.2). The union across rungs is what
reaches deployment-relevant amplitude:

| source | GC span | CAI span |
|---|---|---|
| round 1 (3 optimization algorithms — the panel that produced usable signal) | 0.062 – 0.145 | 0.091 – 0.170 |
| round 2 (single temperature 1.0) — **caused the void** | 0.020 – 0.055 | 0.045 – 0.072 |
| **round 3 ladder (measured on 5 proteins)** | **0.093 – 0.113** | **0.247 – 0.288** |

**No frontier stratum.** Round 2 pre-registered one (`api.candidates`, n = 24) as a
secondary. Measured, it yields **3–4 distinct sequences** per protein after
de-duplication — below the ≥ 12 the analysis requires. It is dropped rather than carried
forward as something that would fail silently.

**UTR contexts — 4**, drawn by the same rule as round 2 with the new seed, requiring the
five-prime UTR lengths to span ≥ 2×. Every variant is scored in all four.

**Cell type — `HEK293T` only** (round 2 §2 records why the all-cell-type read is not
free in BT4's adapter). **Model configuration — bound, not varied:** `species=human`,
`top_k=5`, `batch_size=64`, `num_workers=0`.

**Free-feature set** — GC, GC3, CAI, tAI, CpG density, CDS length.

Scale: 16 × 42 × 4 = **2,688 scorings** ≈ 2 h 50 m at the measured 3.8 s/seq.

## 3. Adequacy — checked on the inputs, not the model

This replaces round 2's `within_over_between` floor, and the reason is the sharpest
lesson of round 2:

> `within_over_between` is computed from **the model's own output**, so it cannot
> distinguish *"the panel is too narrow"* from *"the model is insensitive"*. Those imply
> opposite actions — void the run, or record the finding — and the statistic returns the
> same number for both.

Round 3 checks adequacy on the **sequences**, which is model-independent and cannot be
confounded with what is being measured. Per protein, over its 42 variants:

> **ADEQUATE iff** within-protein **GC span ≥ 0.060** *and* **CAI span ≥ 0.090**.

Both thresholds are the *minimum round 1 actually achieved* (GC 0.0622 on PDE3A, CAI
0.0905 on PDE3A) — i.e. anchored to the panel that produced interpretable signal, not
chosen for convenience. The ladder clears both by a wide margin on every protein measured
(≥ 0.093 / ≥ 0.247), so the floor is **known achievable before it is committed to** —
which is exactly the check round 2 skipped.

If **more than 2 of 16 proteins** fail adequacy, the run is **void**. If 1–2 fail, they
are dropped, the drop is reported, and the analysis proceeds on the rest.

**Harness positive control**, separately: scoring one fixed sequence under two different
UTR contexts must move the score by ≥ 0.01 in absolute value. This is the round-1 Step-10
check, and it is the *only* legitimate model-output sanity test here, because it has a
known expected direction independent of the hypothesis.

## 4. The gate — unchanged from round 2

Uncertainty is a **cluster bootstrap resampling whole proteins** (10,000 draws).

**Free-feature fit.** Within each protein × UTR context, regress RiboNN's score on the six
free features **rank-transformed plus their squares**. What remains is "the residual".
*Void condition:* if the across-protein median |ρ(residual, GC)| exceeds **0.10**, the fit
did not remove GC and the run is void.

**The gate.** Per protein, combine size and stability,

```
stable_non_free = residual_fraction x max(0, median cross-UTR residual rho)
```

and compare against the worst case over a family of **free-feature baseline scorers** run
through the identical pipeline: pure GC, pure CAI, and 24 seeded random linear blends of
all six features.

> **PASS iff** the cluster-bootstrap 95% CI for the across-protein median of
> `stable_non_free(RiboNN) − max over baselines` **excludes 0 from above.**

Both halves are required, and the known-answer regimes in
`scripts/prereg_round2_selftest.py` are why: stability alone passes a pure GC/CAI blend
(after the fit, such a model leaves deterministic dust, identical in every context, so it
correlates at exactly 1.0); size alone passes per-context noise; and size × stability
against *single-feature* baselines still passes a blend, because the fit is in rank space
and a monotone transform of a blend is not linear in the ranked features. The self-test
must pass before the real result is read.

**Stated limitation, unchanged:** this design cannot separate "information beyond the free
features" from "structure the fit failed to capture". The blend family bounds that
misspecification empirically; it does not eliminate it.

## 5. Detection is not sufficiency

At 16 × 42 a very small effect clears zero. The gate decides only **whether a non-free,
stable signal exists**. Whether it is worth a five-figure panel is an explicit **human
judgement on magnitudes, recorded on the record**, reported with CIs: cross-UTR residual
ρ, residual fraction, GC-stratified retention (0 = pure GC detector, 1 = GC-independent),
and the fit diagnostic.

The self-test gives the reference scale: a genuine stable CDS signal reads **+0.83**,
misspecification-only reads **−0.07**. A pass at +0.05 is real and almost certainly not
worth the money.

## 6. Decision rule

- **Void** (adequacy, fit contamination, or failed positive control) → neither pass nor
  fail; diagnose and re-pre-register.
- **Gate fails** → stop. Record it. RiboNN stays `calibrated=False` and `NEXT_SESSION.md`
  item 11 closes with reasons.
- **Gate passes** → record the magnitudes and make the spending judgement explicitly,
  carrying forward that no identified public panel reaches the size floor (PERSIST-seq
  ~77 usable rows / 4 groups; Mauger ~82 rows / 3 proteins). A pass licenses **building or
  hunting for** a panel; it does not make either rejected panel usable.

## 7. What is not allowed

- Changing any threshold, protein, UTR, seed, variant recipe or endpoint after any score
  is seen.
- Re-running with a different ladder, temperatures, contexts or `top_k` and reporting the
  better outcome. The design is the design, not a set of attempts.
- Re-analysing round 2's 2,560 recorded scores under this document's rules. They were
  produced under a voided design; reusing them here would be a search.
- Reporting the gate or the §5 magnitudes selectively.
- Treating a pass as evidence of predictive skill, or any score movement as evidence about
  translation efficiency.

## 8. Code and panel frozen with this document

**The analysis and its self-test are reused *verbatim*** from round 2 —
`scripts/prereg_round2_analyze.py` and `scripts/prereg_round2_selftest.py`, unmodified, at
the hashes recorded in `PREREG_ribonn_part3_round2.md` §7. That is the point rather than a
convenience: **the gate was not retuned after a void.** Anyone can confirm the deciding
code is byte-identical to what was frozen before round 2's data existed.

**Two genuinely new artifacts**, both because an earlier draft of this section claimed
reuse that turned out to be false:

- `scripts/prereg_round3_score.py` — round 2's scorer hardcodes `TEMPERATURE = 1.0` and
  cannot draw a ladder. The new one adds it and runs the §3 adequacy check and positive
  control **before** scoring, so an inadequate panel costs seconds rather than three hours.
- `scripts/prereg_round3_draw_panel.py` — round 2's draw script has the
  validate-then-truncate defect described in §2, which produced a non-compliant panel on
  its first round-3 invocation.

Both are committed in this PR before any scoring. Hashes are of **committed content**
(LF); verify with `git show HEAD:<path> | sha256sum` (round 2 §7 records why).

| file | sha256 |
|---|---|
| `scripts/prereg_round3_draw_panel.py` | `f4c76ec98b59ffc885113df7ce1ba9392f2adcfb5dd12796d365e3b906070d0f` |
| `scripts/prereg_round3_score.py` | `e6a6eb169f291473d58fcc03a7bf0258119d8ad7d63d15c9fdd50c2812fd9296` |
| `scripts/data/prereg_round3_panel.json` | `f0915a84679863429e5829921710c5078f35ede7125f7e9e31cfdfb4a8e55bb0` |

Reused unmodified, at the hashes recorded in `PREREG_ribonn_part3_round2.md` §7:
`scripts/prereg_round2_analyze.py`, `scripts/prereg_round2_selftest.py`.

The panel additionally carries its own `content_hash`
(`104a05c08e922769f0a5fdd04da1c05a06dbef27620dd479d5371fc4e3698386`) over the drawn
proteins, contexts and seeds, from the same pinned NCBI MANE 1.5
(`d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c`, 19,363 Select rows).

## 9. The design was verified achievable before it was frozen

Round 2 pre-registered a floor without checking that its own recipe could clear it, and
voided. Round 3 ran the §3 checks against **this exact panel** first
(`prereg_round3_score.py --check-only`, no scoring):

| check | threshold | observed | margin |
|---|---|---|---|
| within-protein GC span | ≥ 0.060 | 0.0917 – 0.1272 (16/16) | 1.5–2.1× |
| within-protein CAI span | ≥ 0.090 | 0.2464 – 0.3634 (16/16) | 2.7–4.0× |
| inadequate proteins | ≤ 2 | **0** | — |
| positive control \|Δ\| | ≥ 0.01 | 0.0738 | 7.4× |

**The panel, for the record** (16 proteins, 150–1012 aa): PROX1, EDEM2, RIBC1, TMEM79,
AXDND1, PTGFR, ACER1, STAM, MS4A7, PIAS4, GFUS, WBP11, TLCD3B, SLITRK5, CCL25, STAT1.
UTR contexts: HBB (5′ 50 / 3′ 134), ACTB (84 / 600), PLCD1 (103 / 277), ZNF644 (63 / 1494)
— 5′UTR span **2.060×**, satisfying §2's ≥ 2× requirement, checked against the returned set
rather than against a superset of it.

Passing these checks is **not** evidence about RiboNN. It only establishes that a null
result would mean something — that if the gate fails, it failed on a panel that genuinely
exercised the axis, rather than on one too narrow to have shown anything.
