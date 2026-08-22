#!/usr/bin/env python
"""Draw the frozen protein / UTR panel for RiboNN Part 3 round 3.

This script exists so the panel is chosen by a **rule with a recorded seed**, never
by taste. Picking proteins after seeing scores -- or picking them because they look
representative -- is the failure mode the pre-registration
(`docs/PREREG_ribonn_part3_round3.md`) is written to prevent, so the draw is
committed *before* the first scoring run and never re-run to get a different set.

The rule, in full:

1. Download the **release-pinned** NCBI MANE summary (never a moving ``current``
   link, per CLAUDE.md 8) and record its own SHA-256 in the manifest.
2. Keep rows whose ``MANE_status`` is exactly ``MANE Select``, dropping the three
   genes already used in round 1 (BECN1, KRAS, PDE3A) so round 2 is independent.
3. Deterministically permute that list with ``seed``, then walk it **in order**,
   fetching each candidate's protein and UTRs from Ensembl, keeping the first
   ``--n-proteins`` that satisfy every length rule. Walking a seeded permutation
   in order is what keeps the draw reproducible: the accept/reject test never
   looks at anything downstream of the sequence itself.
4. Draw the UTR contexts the same way, with two fixed anchors (HBB, ACTB) plus
   further draws until the four 5'UTR lengths span at least ``--utr5-span``x. The
   span requirement is the point of the contexts -- four near-identical 5'UTRs
   would not test cross-context stability at all.

Length rules are RiboNN's own hard limits plus the pre-registered protein window:
its loader refuses a 5'UTR over 1381 nt or a CDS+3'UTR over 11937 nt, so a context
whose 3'UTR is long enough to push the *longest* drawn protein over that ceiling is
rejected here rather than silently dropping rows at scoring time.

The script **refuses rather than fabricates**: if a candidate cannot be fetched, or
its CDS is not found inside its cDNA, it is skipped and the skip is counted; if the
walk exhausts without filling the panel, it aborts instead of writing a short one.

Usage::

    python scripts/prereg_round3_draw_panel.py --out scripts/data/prereg_round3_panel.json
    python scripts/prereg_round3_draw_panel.py --verify scripts/data/prereg_round3_panel.json
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

MANE_RELEASE = "1.5"
MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
    f"release_{MANE_RELEASE}/MANE.GRCh38.v{MANE_RELEASE}.summary.txt.gz"
)
ENSEMBL = "https://rest.ensembl.org"

# Round 1's proteins, excluded so round 2 is an independent sample.
ROUND1_GENES = frozenset({"BECN1", "KRAS", "PDE3A"})

# Fixed UTR anchors: the pair already used in round 1, plus the alternate context
# from the Step 10 positive control. Keeping them makes round 2 comparable to
# round 1 rather than a wholly separate experiment.
UTR_ANCHORS = (("HBB", "ENST00000335295"), ("ACTB", "ENST00000646664"))

# RiboNN's own loader limits (bt4.biomodels.expression.ribonn).
RIBONN_MAX_UTR5 = 1381
RIBONN_MAX_CDS_PLUS_UTR3 = 11937

STOPS = ("TAA", "TGA", "TAG")


def _curl(url: str) -> bytes:
    """Fetch via the system ``curl``.

    Needed because this environment's Python SSL stack fails against
    ``ftp.ncbi.nlm.nih.gov`` with ``[ASN1: NOT_ENOUGH_DATA]`` while ``curl`` against
    the same URL returns 200. Rather than disable certificate verification -- which
    would silently accept any peer for a file we then hash and pin -- fall back to a
    client that completes the handshake properly. ``curl`` ships with Windows 10+,
    macOS and essentially every Linux.
    """
    out = subprocess.run(
        ["curl", "-sSL", "--fail", "--max-time", "300", url],
        capture_output=True,
        check=True,
    )
    return out.stdout


def _get(url: str, *, tries: int = 4) -> bytes:
    """Fetch ``url``, retrying politely. Raises on definitive failure."""
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "bt4-prereg/1"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                return bytes(fh.read())
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404):  # a definitive "no such thing"
                raise
            last = exc
        except Exception as exc:
            last = exc
            try:
                return _curl(url)
            except Exception as curl_exc:
                last = curl_exc
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {url}: {last}")


def _fasta_body(raw: bytes) -> str:
    return "".join(
        line.strip() for line in raw.decode().splitlines() if not line.startswith(">")
    ).upper()


def download_mane() -> tuple[list[dict[str, str]], str]:
    """Return (MANE Select rows, sha256 of the pinned summary file)."""
    blob = _get(MANE_URL)
    sha = hashlib.sha256(blob).hexdigest()
    text = gzip.decompress(blob).decode()
    lines = text.splitlines()
    header = lines[0].lstrip("#").split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        rec = dict(zip(header, line.split("\t"), strict=False))
        if rec.get("MANE_status") == "MANE Select":
            rows.append(rec)
    return rows, sha


def fetch_transcript(enst: str) -> dict[str, Any] | None:
    """Return protein / CDS / UTRs for one Ensembl transcript, or None if unusable."""
    base = enst.split(".")[0]
    try:
        fmt = "content-type=text/x-fasta"
        cdna = _fasta_body(_get(f"{ENSEMBL}/sequence/id/{base}?type=cdna;{fmt}"))
        cds = _fasta_body(_get(f"{ENSEMBL}/sequence/id/{base}?type=cds;{fmt}"))
    except Exception:
        return None
    if not cds or not cdna or set(cds) - set("ACGT") or set(cdna) - set("ACGT"):
        return None
    if len(cds) % 3 or cds[-3:] not in STOPS or not cds.startswith("ATG"):
        return None
    if cdna.count(cds) != 1:  # ambiguous or absent -> refuse, do not guess
        return None
    i = cdna.index(cds)
    utr5, utr3 = cdna[:i], cdna[i + len(cds) :]
    from bt4.domain.genetic_code import translate

    protein = translate(cds[:-3])
    if "*" in protein:  # an internal stop means the frame is not what we think
        return None
    return {
        "transcript": base,
        "protein": protein,
        "aa_len": len(protein),
        "cds_len": len(cds),
        "utr5": utr5,
        "utr3": utr3,
        "utr5_len": len(utr5),
        "utr3_len": len(utr3),
    }


def draw_proteins(
    rows: list[dict[str, str]], seed: int, n: int, min_aa: int, max_aa: int
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    eligible = [r for r in rows if r.get("symbol") and r["symbol"] not in ROUND1_GENES]
    order = list(range(len(eligible)))
    random.Random(seed).shuffle(order)
    picked: list[dict[str, Any]] = []
    seen_symbols: set[str] = set()
    skips = {"fetch_failed": 0, "length_window": 0, "duplicate_symbol": 0, "examined": 0}
    for idx in order:
        if len(picked) >= n:
            break
        rec = eligible[idx]
        symbol = rec["symbol"]
        if symbol in seen_symbols:
            skips["duplicate_symbol"] += 1
            continue
        skips["examined"] += 1
        info = fetch_transcript(rec["Ensembl_nuc"])
        if info is None:
            skips["fetch_failed"] += 1
            continue
        if not (min_aa <= info["aa_len"] <= max_aa):
            skips["length_window"] += 1
            continue
        seen_symbols.add(symbol)
        info["symbol"] = symbol
        info.pop("utr5", None)  # a protein's own UTRs are not used; contexts supply them
        info.pop("utr3", None)
        picked.append(info)
        time.sleep(0.12)  # be polite to Ensembl
    return picked, skips


def draw_utr_contexts(
    rows: list[dict[str, str]],
    seed: int,
    n: int,
    span: float,
    max_cds_len: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    contexts: list[dict[str, Any]] = []
    skips = {"fetch_failed": 0, "too_long_for_ribonn": 0, "empty_utr": 0, "examined": 0}
    for name, enst in UTR_ANCHORS:
        info = fetch_transcript(enst)
        if info is None:
            raise RuntimeError(f"fixed UTR anchor {name} ({enst}) could not be fetched")
        contexts.append({"name": name, **info})

    by_symbol = {r["symbol"]: r for r in rows if r.get("symbol")}
    order = sorted(by_symbol)
    random.Random(seed + 1).shuffle(order)
    n_anchor = len(UTR_ANCHORS)
    for symbol in order:
        if len(contexts) == n and _span_ok(contexts, span):
            break
        if symbol in ROUND1_GENES or any(c["name"] == symbol for c in contexts):
            continue
        skips["examined"] += 1
        info = fetch_transcript(by_symbol[symbol]["Ensembl_nuc"])
        if info is None:
            skips["fetch_failed"] += 1
            continue
        if info["utr5_len"] == 0 or info["utr3_len"] == 0:
            skips["empty_utr"] += 1
            continue
        if info["utr5_len"] > RIBONN_MAX_UTR5:
            skips["too_long_for_ribonn"] += 1
            continue
        if max_cds_len + info["utr3_len"] > RIBONN_MAX_CDS_PLUS_UTR3:
            skips["too_long_for_ribonn"] += 1
            continue
        candidate = {"name": symbol, **info}
        if len(contexts) < n:
            contexts.append(candidate)
        else:
            # Full but the span still fails. Swap the candidate into whichever
            # NON-ANCHOR slot most improves the span, and keep it only if it does.
            # (The anchors are fixed by design and are never displaced.)
            best_span, best_set = _span_of(contexts), None
            for i in range(n_anchor, len(contexts)):
                trial = [*contexts[:i], candidate, *contexts[i + 1 :]]
                if _span_of(trial) > best_span:
                    best_span, best_set = _span_of(trial), trial
            if best_set is not None:
                contexts = best_set
        time.sleep(0.12)
    # Validate exactly what is returned. The round-2 version validated the full
    # collected list and then returned `contexts[:n]`, so a truncation could drop the
    # very context that satisfied the span -- which is how the first round-3 draw
    # produced 5'UTRs of 50/84/69/63 (1.68x) while reporting success.
    if len(contexts) != n or not _span_ok(contexts, span):
        raise RuntimeError(
            f"could not assemble {n} UTR contexts spanning >= {span}x in 5'UTR length; "
            f"got {len(contexts)} ({[c['utr5_len'] for c in contexts]})"
        )
    return contexts, skips


def _span_of(contexts: list[dict[str, Any]]) -> float:
    """Ratio of longest to shortest 5'UTR across the set (0.0 if degenerate)."""
    lens = [c["utr5_len"] for c in contexts]
    if not lens or min(lens) <= 0:
        return 0.0
    return max(lens) / min(lens)


def _span_ok(contexts: list[dict[str, Any]], span: float) -> bool:
    return _span_of(contexts) >= span


def build(args: argparse.Namespace) -> dict[str, Any]:
    rows, mane_sha = download_mane()
    print(f"MANE Select rows: {len(rows)}  sha256={mane_sha}", file=sys.stderr)

    proteins, pskips = draw_proteins(rows, args.seed, args.n_proteins, args.min_aa, args.max_aa)
    if len(proteins) < args.n_proteins:
        raise RuntimeError(
            f"exhausted MANE without filling the panel: {len(proteins)}/{args.n_proteins}"
        )
    max_cds = max(p["cds_len"] for p in proteins)
    contexts, uskips = draw_utr_contexts(rows, args.seed, args.n_utr, args.utr5_span, max_cds)

    rng = random.Random(args.seed + 2)
    variant_seeds = {p["symbol"]: rng.randrange(2**31) for p in proteins}
    return {
        "schema": "bt4.prereg.round2.panel/1",
        "pre_registration": "docs/PREREG_ribonn_part3_round2.md",
        "seed": args.seed,
        "mane": {
            "release": MANE_RELEASE,
            "url": MANE_URL,
            "sha256": mane_sha,
            "n_mane_select": len(rows),
        },
        "filters": {
            "min_aa": args.min_aa,
            "max_aa": args.max_aa,
            "excluded_genes": sorted(ROUND1_GENES),
            "protein_draw_skips": pskips,
            "utr_draw_skips": uskips,
        },
        "proteins": proteins,
        "utr_contexts": contexts,
        "variant_seeds": variant_seeds,
    }


def content_hash(manifest: dict[str, Any]) -> str:
    """Stable hash of the drawn panel (not of the incidental skip tallies)."""
    core = {
        "seed": manifest["seed"],
        "mane_sha256": manifest["mane"]["sha256"],
        "proteins": [[p["symbol"], p["transcript"], p["protein"]] for p in manifest["proteins"]],
        "utr_contexts": [
            [c["name"], c["transcript"], c["utr5"], c["utr3"]] for c in manifest["utr_contexts"]
        ],
        "variant_seeds": manifest["variant_seeds"],
    }
    blob = json.dumps(core, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, help="write the manifest here")
    ap.add_argument("--verify", type=Path, help="re-draw and diff against this manifest")
    ap.add_argument("--seed", type=int, default=20260822)
    ap.add_argument("--n-proteins", type=int, default=16)
    ap.add_argument("--n-utr", type=int, default=4)
    ap.add_argument("--min-aa", type=int, default=150)
    ap.add_argument("--max-aa", type=int, default=1200)
    ap.add_argument("--utr5-span", type=float, default=2.0)
    args = ap.parse_args()
    if not args.out and not args.verify:
        ap.error("one of --out / --verify is required")

    manifest = build(args)
    manifest["content_hash"] = content_hash(manifest)

    if args.verify:
        old = json.loads(args.verify.read_text())
        same = old.get("content_hash") == manifest["content_hash"]
        print(f"committed: {old.get('content_hash')}")
        print(f"rebuilt  : {manifest['content_hash']}")
        print("MATCH" if same else "DIFFERENT -- the draw is not reproducible")
        return 0 if same else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", newline="\n")
    print(f"wrote {args.out}  content_hash={manifest['content_hash']}")
    print(f"proteins: {', '.join(p['symbol'] for p in manifest['proteins'])}")
    print(
        "UTR contexts: "
        + ", ".join(
            f"{c['name']}(5'={c['utr5_len']},3'={c['utr3_len']})"
            for c in manifest["utr_contexts"]
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
