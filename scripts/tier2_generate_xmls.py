#!/usr/bin/env python
"""
Generate Tier-2 timing benchmark XML variants for both example datasets, at
two grid baselines:

  * baseline=user1x  — user's published gridRateShifts and mascotshifts
                       (the practical drop-in target).
  * baseline=user5x  — both grids densified 5x. This shifts the converged
                       limit (the spline interpolation is closer to the cubic
                       spline) and lets us check whether the new path's
                       speed-up grows when we ask for finer accuracy.

For each baseline and dataset we emit:
  * one OLD config: the baseline itself, useMaxInterval=false.
  * eight NEW configs at maxInterval ∈ {1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4,
    1e-4, 1e-5}. The values bracket each XML's mascotshifts step (2e-4 for
    the 6_2 small XML, 5e-3 for SARS) so that "match user-grid logP" sits
    inside the swept range.

Naming: <dataset>__<baseline>__<mode>__<param>.xml

Tier2Runner.java parses the file name on double-underscore — we have THREE
labels here (baseline, mode, param) so the runner needs the matching update;
done in the same patch.

Usage:  conda run -n biopython_env python scripts/tier2_generate_xmls.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "sandbox" / "tier2"

DATASETS: list[tuple[str, Path]] = [
    ("6_2_small", REPO / "examples" / "6_2_simulation_datastreams.xml"),
    ("sarscov2_1000seq", REPO / "examples" /
     "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml"),
]

# (baseline_label, grid_factor)
BASELINES: list[tuple[str, int]] = [
    ("user1x", 1),
    ("user5x", 5),
]

# (mode, param, useMaxInterval, maxInterval)
# - "old" rows reproduce the baseline grid via the legacy path (no maxInterval).
# - "new" rows use the maxInterval branch with rates queried by time.
MAX_INTERVAL_VALUES = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 1e-5]


def _label(mi: float) -> str:
    """Stable filename-safe label for a maxInterval value, e.g. 5e-3 -> 'max5e-3'."""
    s = f"{mi:.0e}".replace("+0", "").replace("-0", "-").replace("+", "")
    # Drop leading zeros in mantissa: 1e-3 stays 1e-3, but format may give 1e-03.
    return "max" + s


VARIANTS: list[tuple[str, str, bool, float | None]] = (
    [("old", "grid", False, None)]
    + [("new", _label(mi)[3:] and _label(mi).removeprefix("max") and f"max{mi:g}", True, mi)
       for mi in MAX_INTERVAL_VALUES]
)
# The list-comprehension above is overly clever; rebuild it simply:
VARIANTS = [("old", "grid", False, None)] + [
    ("new", f"max{mi:g}", True, mi) for mi in MAX_INTERVAL_VALUES
]


def densify_one_grid(values: list[float], factor: int) -> list[float]:
    if factor <= 1:
        return values
    out: list[float] = []
    for i in range(len(values) - 1):
        a, b = values[i], values[i + 1]
        for k in range(factor):
            out.append(a + (b - a) * k / factor)
    out.append(values[-1])
    return out


# Match <gridRateShifts ...>body</gridRateShifts> (spline evaluation grid)
# AND <rateShifts id="mascotshifts" ...>body</rateShifts> (dynamics outer
# iteration grid for SARS-like setups). We do NOT touch SkygrowthRateShifts —
# those are spline knot times paired with a parameter whose dimension would
# stop matching. In the small XML the spline eval grid also serves as the
# dynamics outer grid (via <rateShifts idref=...>), so densifying it once is
# enough.
GRID_PATTERN = re.compile(
    r'(<gridRateShifts\b[^>]*>)([^<]*)(</gridRateShifts>)', re.DOTALL,
)
MASCOT_SHIFTS_PATTERN = re.compile(
    r'(<rateShifts\s+id="mascotshifts"[^>]*>)([^<]*)(</rateShifts>)', re.DOTALL,
)


def patch_grids(xml_text: str, factor: int) -> str:
    if factor <= 1:
        return xml_text

    def repl(m: re.Match[str]) -> str:
        open_tag, body, close_tag = m.group(1), m.group(2), m.group(3)
        values = [float(x) for x in body.split()]
        dense = densify_one_grid(values, factor)
        return open_tag + " ".join(f"{v:.6f}" for v in dense) + close_tag

    new_text, n_grid = GRID_PATTERN.subn(repl, xml_text)
    new_text, n_mascot = MASCOT_SHIFTS_PATTERN.subn(repl, new_text)
    if n_grid + n_mascot == 0:
        sys.exit("[tier2_generate_xmls] no gridRateShifts or mascotshifts elements found")
    return new_text


MASCOT_PATTERN = re.compile(
    r'(<distribution\s+id="Mascot\.t:SimDataset"\s+spec="mascotdatastreams\.distribution\.MascotLogPflag"[^>]*?)(>)',
)


def patch_mascot(xml_text: str, use_max_interval: bool, max_interval: float | None) -> str:
    m = MASCOT_PATTERN.search(xml_text)
    if not m:
        sys.exit("[tier2_generate_xmls] could not find MascotLogPflag <distribution> tag")
    head, close = m.group(1), m.group(2)
    extra = ""
    if use_max_interval:
        if max_interval is None:
            sys.exit("[tier2_generate_xmls] use_max_interval requires max_interval")
        extra = f' implementation="java" useMaxInterval="true" maxInterval="{max_interval}"'
    return xml_text[: m.start()] + head + extra + close + xml_text[m.end():]


COUPLED_RUN_PATTERN = re.compile(
    r'<run\s+id="mcmc"\s+spec="coupledMCMC\.CoupledMCMC"([^>]*?)>',
)


def patch_run_spec(xml_text: str) -> str:
    """Rewrite coupledMCMC.CoupledMCMC -> plain MCMC. The harness only walks
    the parsed BEASTObject graph — it never runs the chain — so this just
    avoids needing the CoupledMCMC package on classpath."""
    def repl(m: re.Match[str]) -> str:
        attrs = m.group(1)
        cl = re.search(r'chainLength="([^"]+)"', attrs)
        return f'<run id="mcmc" spec="MCMC" chainLength="{cl.group(1) if cl else "1"}">'

    return COUPLED_RUN_PATTERN.sub(repl, xml_text)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    n_written = 0
    for dataset, src in DATASETS:
        if not src.is_file():
            print(f"[tier2_generate_xmls] WARNING: missing {src} — skipping {dataset}")
            continue
        base_text = patch_run_spec(src.read_text())
        for baseline_label, factor in BASELINES:
            scaled_text = patch_grids(base_text, factor)
            for mode, param, use_mi, mi in VARIANTS:
                text = patch_mascot(scaled_text, use_mi, mi)
                out = OUT_DIR / f"{dataset}__{baseline_label}__{mode}__{param}.xml"
                out.write_text(text)
                print(f"  wrote {out.relative_to(REPO)}")
                n_written += 1
    print(f"[tier2_generate_xmls] {n_written} variants written under {OUT_DIR.relative_to(REPO)}")


if __name__ == "__main__":
    main()
