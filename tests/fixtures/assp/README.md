# Offline ASSP cross-check fixtures

These files are **synthetic, illustrative** ASSP-format report bodies used to drive
`bt4.biomodels.splice.assp.FixtureAsspTransport` (selected via the
`BT4_ASSP_FIXTURE_DIR` environment variable) so the ASSP splice cross-check is
exercised **entirely offline** in CI — no live network call is ever made (CLAUDE.md
§6, §10.15).

**They are NOT real captured ASSP responses.** The live ASSP web service
(http://wangcomputing.com/assp/) was unreachable during development, and capturing
a real response requires live network access CI forbids — exactly as no reference
panel ships for the wrapped SpliceAI / Pangolin CNNs. These fixtures reproduce the
*shape* of ASSP's tabular site report (position / site-type / score / confidence
columns, with constitutive / alternative-isoform / cryptic donor & acceptor
classifications) purely to test the adapter's parsing, caching, backoff, and
graceful-degradation logic deterministically. The site scores and confidences are
made-up, not measured, and are never presented as calibrated splice probabilities.

## File naming

Each file is named `<cache_key>.txt`, where `cache_key` is the lowercase hex
SHA-256 of the upper-cased, validated coding sequence
(`bt4.biomodels.splice.assp.cache_key`). Current fixtures:

| cache_key (file stem) | sequence | contents |
|---|---|---|
| `74c17699…e43b…1161c` | `ATGGCCGGCGATCGATCGATCGTAA` | three predicted sites (constitutive donor, cryptic acceptor, alternative donor) |
| `632320c4…4e43ab` | `ATGAAATTTGGGCCCTAA` | header present, **no** predicted sites (a legitimate empty result) |

To add a fixture, compute `cache_key(sequence)` and write the report to
`<key>.txt`.
