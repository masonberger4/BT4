# Result — RiboNN Part 3, round 2: **VOID**

The pre-registered design is [`PREREG_ribonn_part3_round2.md`](PREREG_ribonn_part3_round2.md),
frozen at merge (#141) before any sequence was scored. This file records what happened.
Per §6 of that document, everything is reported whatever the verdict.

## Verdict

**The run is void — neither a pass nor a fail.** It failed the pre-registered sanity
floor:

```
within_over_between = 0.180  <  0.20        →  void
```

§3 binds that outcome in advance: failing the floor "voids the run … rather than making
it negative". No conclusion about RiboNN may be drawn from it.

**The gate came back strongly positive.** Recorded here because suppressing it would be
worse than stating it, and because this is precisely the configuration a pre-registration
exists to survive:

| statistic | median | 95% CI (cluster bootstrap over proteins) |
|---|---|---|
| **gate** (stable non-free vs. baseline floor) | **+0.4124** | [0.3292, 0.5534] |
| cross-UTR residual ρ | 0.9561 | [0.9393, 0.9672] |
| residual fraction | 0.6207 | [0.5249, 0.7842] |
| GC-stratified retention | 0.9594 | [0.9137, 0.9812] |
| ρ(residual, GC) — fit diagnostic | 0.0241 | [0.0157, 0.0274] |
| stable non-free, RiboNN | 0.5770 | — |
| stable non-free, baseline floor (26 baselines) | 0.1285 | — |

All 16 proteins were individually positive (+0.135 … +0.785). Absent the void this
passes. **It is not a result**, and neither is the `retention ≈ 0.96` that would appear
to cut against round 1's GC-detector reading. Both come from a void run in a regime
established below to carry 3.6× less within-protein signal than round 1.

## Why the floor failed — measured, not assumed

The first explanation offered was that the denominator grew, because 16 diverse MANE
proteins spread wider than round 1's three. **That was wrong, and testing it took
minutes.** The denominator shrank too; the numerator collapsed far harder:

| | round 1 (Ranaghan, 3 proteins) | round 2 (16 MANE proteins) | change |
|---|---|---|---|
| median within-protein score SD | 0.2524 | **0.0706** | −72% |
| between-protein score SD | 0.6213 | 0.3935 | −37% |
| ratio | 0.4062 | **0.1795** | fails |

The ratio is 0.1795–0.1919 across all four UTR contexts, so it is not one context's fluke.

**Root cause: the variant set was too narrow.** Round 2 drew 40 variants per protein from
`api.library` at **temperature 1.0** — samples from the natural codon distribution, which
cluster tightly. Round 1's 31 spellings per protein were outputs of *three different
codon-optimization algorithms*, which are extremal by construction. Measured
within-protein feature spans:

| source | GC span | CAI span |
|---|---|---|
| round 1 (3 optimization algorithms) | 0.062 – 0.145 | 0.091 – 0.170 |
| round 2 (`library`, T = 1.0) | **0.020 – 0.055** | 0.045 – 0.072 |

RiboNN's response was small because the design space it was shown was small. This is a
defect in the round-2 design, and an avoidable one: that document itself calls the
frontier "what BT4 actually delivers" and then makes temperature-1.0 sampling the primary
stratum anyway.

## What the floor got right, and what it got wrong

**Right:** it caught that the panel under-sampled the axis under test, which is exactly
what a sanity floor is for. Without it, round 2 would have reported a passing gate from a
regime that does not resemble deployment.

**Wrong, and this is the deeper defect:** `within_over_between` is computed **from the
model's own output**, so it *cannot distinguish* "the panel is too narrow" from "the model
is insensitive". Those are opposite conclusions — one voids the run, the other *is* the
finding — and the statistic returns the same number for both. Round 2 was lucky that the
input-side evidence (feature spans above) resolved the ambiguity after the fact.

Round 3 therefore replaces it with an adequacy check on the **inputs**: the within-protein
spread of GC and CAI, which is model-independent and cannot be confounded with the thing
being measured.

## Two further defects found while running it

- **The pre-registered frontier stratum could never have worked.** §2 specifies
  `api.candidates(protein, steps=11, n=24)` as a secondary stratum. Measured, it returns
  **3–4 distinct sequences** per protein after de-duplication — below the ≥12 the analysis
  requires per protein × context. It was bound as "never decisive", so nothing rests on
  it, but it would have failed silently had it been run.
- **The freeze hashes did not verify on Windows.** `core.autocrlf=true` with no
  `.gitattributes` meant checkout rewrote LF→CRLF, so `sha256sum` over the working tree
  hashed different bytes than the committed content — indistinguishable, to a reader, from
  a violated freeze. Fixed in #142; the committed blobs were always correct.

## Cost and provenance

2,560 scorings (16 proteins × 40 variants × 4 UTR contexts), 9,623 s ≈ 2 h 40 m at
3.49–3.88 s/seq. Scores: `round2_scores.jsonl`, sha256
`d973f403c6fc326e39354c83d23682e5e4c3d56628ba8c599807895552152226`, every row
`calibrated: false`. Panel `content_hash`
`b7da49ea1bd241994cd875b13101215259a7c81f055444461801e8f6916a5359`, drawn from NCBI MANE
1.5. The analysis ran from the frozen `scripts/prereg_round2_analyze.py`; its stale text
label (it prints the name of the *earlier, rejected* gate) is a display defect only — the
`--json` output above is the record, and the script was deliberately **not** edited
mid-run so the freeze stays checkable.

## Standing conclusion

**RiboNN remains `calibrated=False`.** Round 1's verdict — stop, on a GC3 confound of
0.7385 against a pre-registered 0.70 — is still the only completed verdict, and round 2
neither confirms nor overturns it. The open question moves to
[`PREREG_ribonn_part3_round3.md`](PREREG_ribonn_part3_round3.md).
