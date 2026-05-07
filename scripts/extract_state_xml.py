#!/usr/bin/env python
"""
Build an XML that uses the LAST tree from a Tier-3 BEAST run as a fixed
TreeParser initialiser, so the IntervalDumper can be run on the
equilibrated-tree state without re-running MCMC.

Reads:
  - source XML        (existing tier3 input XML)
  - .trees file       (the BEAST `-SimDataset.trees` log; takes the last sample)
Writes:
  - modified XML in sandbox/tier3_state/

Replaces the <init spec="RandomTree" ...>...</init> block with
<init spec="TreeParser" ... newick="..."/> using the last tree from the run,
with numeric Newick labels expanded via the trees file's Translate block.

Note: parameter values are *not* updated. ODE step counts depend on the tree
and on grid/maxInterval timing, NOT on parameter values, so this is sufficient
for measuring n_doEuler at the equilibrated tree.

Usage:
    conda run -n biopython_env python scripts/extract_state_xml.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "sandbox" / "tier3_state"

CASES = [
    {
        "label": "old",
        "src_xml": REPO / "sandbox" / "tier3_old"
                  / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml",
        "trees": REPO / "sandbox" / "tier3_old"
                / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip-SimDataset.trees",
        "out_xml": OUT_DIR / "sars_state_old.xml",
    },
    {
        "label": "new",
        "src_xml": REPO / "sandbox" / "tier3_new"
                  / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip_interval.xml",
        "trees": REPO / "sandbox" / "tier3_new"
                / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip_interval-SimDataset.trees",
        "out_xml": OUT_DIR / "sars_state_new.xml",
    },
]


def parse_translate_and_last_tree(trees_path: Path) -> tuple[dict[int, str], int, str]:
    """Return (taxon_map, last_state_number, last_newick).

    taxon_map: int (numeric label in tree) -> str (taxon name)
    """
    text = trees_path.read_text()
    # The Translate block: `Translate\n\t1 name1,\n\t2 name2,\n\t...\n\t999 nameN\n;`
    m = re.search(r"Translate\s+(.*?)\s*;", text, re.DOTALL)
    if not m:
        sys.exit(f"[extract_state_xml] no Translate block in {trees_path}")
    taxon_map: dict[int, str] = {}
    for line in m.group(1).split(","):
        parts = line.strip().split(None, 1)
        if len(parts) != 2:
            continue
        idx, name = parts
        # Strip quotes if present.
        name = name.strip().strip("'").strip('"')
        taxon_map[int(idx)] = name

    # Last tree line: `tree STATE_<N> = (...);`
    last_state = -1
    last_newick = ""
    for m in re.finditer(r"tree\s+STATE_(\d+)\s*=\s*([^;]*);", text):
        s = int(m.group(1))
        if s > last_state:
            last_state = s
            last_newick = m.group(2).strip()
    if last_state < 0:
        sys.exit(f"[extract_state_xml] no STATE_* tree found in {trees_path}")
    return taxon_map, last_state, last_newick


def expand_numeric_labels(newick: str, taxon_map: dict[int, str]) -> str:
    """Replace numeric tip labels in Newick with the mapped taxon names."""
    # In Newick, a numeric label appears as `(` or `,` then a number then `:`.
    pat = re.compile(r"([(,])(\d+)(:|\)|,)")

    def repl(m: re.Match[str]) -> str:
        n = int(m.group(2))
        if n in taxon_map:
            return f"{m.group(1)}{taxon_map[n]}{m.group(3)}"
        return m.group(0)

    # Run twice — single pass misses overlapping matches at boundaries.
    out = pat.sub(repl, newick)
    out = pat.sub(repl, out)
    return out


# Multi-line match for the RandomTree block, including its nested
# <populationModel>...</populationModel>.
RANDOMTREE_PATTERN = re.compile(
    r'<init\s+id="RandomTree\.t:SimDataset"[^>]*>.*?</init>',
    re.DOTALL,
)


def replace_randomtree(xml_text: str, newick: str) -> str:
    """Replace the RandomTree <init> block with a TreeParser init using `newick`."""
    # XML attribute escaping: ampersands, quotes are the things to worry about
    # in Newick; standard newick has none of those, but we're cautious anyway.
    safe = (newick.replace("&", "&amp;").replace('"', "&quot;"))
    replacement = (
        f'<init id="NewickTree.t:SimDataset" '
        f'spec="beast.base.evolution.tree.TreeParser" '
        f'adjustTipHeights="false" '
        f'initial="@Tree.t:SimDataset" '
        f'taxa="@SimDataset" '
        f'IsLabelledNewick="true" '
        f'newick="{safe}"/>'
    )
    new_text, n = RANDOMTREE_PATTERN.subn(replacement, xml_text)
    if n == 0:
        sys.exit("[extract_state_xml] could not find RandomTree <init> block")
    return new_text


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        if not case["src_xml"].is_file():
            print(f"[skip {case['label']}] missing {case['src_xml']}")
            continue
        if not case["trees"].is_file():
            print(f"[skip {case['label']}] missing {case['trees']}")
            continue

        taxon_map, state_n, newick = parse_translate_and_last_tree(case["trees"])
        print(f"[{case['label']}] last tree = STATE_{state_n}, "
              f"{len(taxon_map)} taxa in translate block")
        expanded = expand_numeric_labels(newick, taxon_map)
        # Sanity: any remaining numeric leaf labels means the regex didn't catch them.
        leftover_numeric = re.findall(r"[(,](\d+)[:),]", expanded)
        if leftover_numeric:
            print(f"  WARNING: {len(leftover_numeric)} numeric labels still in newick "
                  f"after expansion (first few: {leftover_numeric[:5]})")

        xml_text = case["src_xml"].read_text()
        new_xml = replace_randomtree(xml_text, expanded)
        case["out_xml"].write_text(new_xml)
        print(f"  wrote {case['out_xml'].relative_to(REPO)}")


if __name__ == "__main__":
    main()
