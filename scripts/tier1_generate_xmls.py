#!/usr/bin/env python
"""
Generate Tier-1 benchmark XML variants for the max-interval ODE-subdivision
correctness check.

Reads:   examples/6_2_simulation_datastreams.xml
Writes:  sandbox/tier1/<base>__<mode>__<param>.xml

Variants
--------
old    grid1x           original gridRateShifts, useMaxInterval=false (baseline)
old    grid5x           gridRateShifts densified 5x, useMaxInterval=false (fine reference)
new    max1e-2 ... 1e-5 useMaxInterval=true at four resolutions

The harness (Tier1Runner.java) parses the file name to label CSV rows; the
double-underscore separator is load-bearing.

Usage:  python scripts/tier1_generate_xmls.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC_XML = REPO / "examples" / "6_2_simulation_datastreams.xml"
OUT_DIR = REPO / "sandbox" / "tier1"
BASE = SRC_XML.stem  # "6_2_simulation_datastreams"

# Variants: (mode, param, gridDensify, useMaxInterval, maxInterval)
VARIANTS: list[tuple[str, str, int, bool, float | None]] = [
    ("old", "grid1x", 1, False, None),
    ("old", "grid5x", 5, False, None),
    ("new", "max1e-2", 1, True, 1e-2),
    ("new", "max1e-3", 1, True, 1e-3),
    ("new", "max1e-4", 1, True, 1e-4),
    ("new", "max1e-5", 1, True, 1e-5),
]


def densify_grid(values: list[float], factor: int) -> list[float]:
    """Insert (factor-1) equally-spaced points between each consecutive pair."""
    if factor <= 1:
        return values
    out: list[float] = []
    for i in range(len(values) - 1):
        a, b = values[i], values[i + 1]
        for k in range(factor):
            out.append(a + (b - a) * k / factor)
    out.append(values[-1])
    return out


def patch_grid(xml_text: str, factor: int) -> str:
    """Replace the body of <gridRateShifts id="SplineGridRateShifts" ...>."""
    if factor <= 1:
        return xml_text
    pattern = re.compile(
        r'(<gridRateShifts\s+id="SplineGridRateShifts"[^>]*>)([^<]*)(</gridRateShifts>)',
        flags=re.DOTALL,
    )
    m = pattern.search(xml_text)
    if not m:
        sys.exit("[tier1_generate_xmls] could not find gridRateShifts body to densify")
    open_tag, body, close_tag = m.group(1), m.group(2), m.group(3)
    values = [float(x) for x in body.split()]
    dense = densify_grid(values, factor)
    new_body = " ".join(f"{v:.6f}" for v in dense)
    return xml_text[: m.start()] + open_tag + new_body + close_tag + xml_text[m.end():]


def patch_mascot(xml_text: str, use_max_interval: bool, max_interval: float | None) -> str:
    """Add useMaxInterval / maxInterval attributes to the MascotLogPflag <distribution>."""
    pattern = re.compile(
        r'(<distribution\s+id="Mascot\.t:SimDataset"\s+spec="mascotdatastreams\.distribution\.MascotLogPflag"[^>]*?)(>)',
    )
    m = pattern.search(xml_text)
    if not m:
        sys.exit("[tier1_generate_xmls] could not find MascotLogPflag <distribution> tag")
    head, close = m.group(1), m.group(2)
    extra = ""
    if use_max_interval:
        if max_interval is None:
            sys.exit("[tier1_generate_xmls] use_max_interval=True requires max_interval")
        # implementation="java" is required by the validation guard in
        # MascotLogPflag.initAndValidate() when useMaxInterval=true.
        extra = f' implementation="java" useMaxInterval="true" maxInterval="{max_interval}"'
    return xml_text[: m.start()] + head + extra + close + xml_text[m.end():]


def main() -> None:
    if not SRC_XML.is_file():
        sys.exit(f"[tier1_generate_xmls] missing source XML: {SRC_XML}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base_text = SRC_XML.read_text()

    for mode, param, factor, use_mi, mi in VARIANTS:
        text = patch_grid(base_text, factor)
        text = patch_mascot(text, use_mi, mi)
        out = OUT_DIR / f"{BASE}__{mode}__{param}.xml"
        out.write_text(text)
        print(f"  wrote {out.relative_to(REPO)}")

    print(f"[tier1_generate_xmls] {len(VARIANTS)} variants written under {OUT_DIR.relative_to(REPO)}")


if __name__ == "__main__":
    main()
