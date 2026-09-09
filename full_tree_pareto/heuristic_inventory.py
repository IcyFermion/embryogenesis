"""Presentation scope and code-inspection provenance for full-tree heuristics.

This inventory is not a certification of historical caches. Only the two main
methods are regenerated with structural invariants by the publication pipeline.
"""

from pathlib import Path

import pandas as pd


INVENTORY = [
    {
        "heuristic": "Layerwise assignment", "presentation": "Main + supplement",
        "internal_positions": "Measured", "internal_proteins": "Measured",
        "roots": "Four measured roots fixed", "scored_edges": "1,000; validated",
        "constraints_and_status": "Exact Hungarian assignment within each contraction round; natural parent-slot multiplicities fixed. Not a global forest optimum.",
    },
    {
        "heuristic": "Degree-constrained spanning forest", "presentation": "Main + supplement",
        "internal_positions": "Measured", "internal_proteins": "Measured",
        "roots": "Four measured roots fixed", "scored_edges": "1,000; validated",
        "constraints_and_status": "Greedy constrained Kruskal plus terminal assignment; four components and at most two children per parent validated at every weight.",
    },
    {
        "heuristic": "Top-down rebuild", "presentation": "Supplement",
        "internal_positions": "Measured", "internal_proteins": "Measured",
        "roots": "Fixed in current code", "scored_edges": "1,000 intended; cache",
        "constraints_and_status": "Greedy rooted growth plus terminal assignment. Near-overlap with constrained forest; historical scalar-cost cache does not store reconstructed edges for revalidation.",
    },
    {
        "heuristic": "Bottom-up by layer", "presentation": "Supplement",
        "internal_positions": "Measured; identities reassigned", "internal_proteins": "Measured; identities reassigned",
        "roots": "Root slots fixed; identities may move", "scored_edges": "1,000 on fixed topology",
        "constraints_and_status": "Sequential assignment of internal identities on fixed topology; current pool includes root identities. Cached costs; broader identity freedom than the main methods.",
    },
    {
        "heuristic": "Paired bottom-up", "presentation": "Supplement",
        "internal_positions": "Measured", "internal_proteins": "Measured",
        "roots": "Biological root identities not enforced", "scored_edges": "1,000 intended; scope unverified",
        "constraints_and_status": "Pair children, then assign measured parents. Current code permits root identities in the parent pool; coverage and final components require a separate cache audit.",
    },
    {
        "heuristic": "Terminal-only rebuild", "presentation": "Held; table only",
        "internal_positions": "Constructed daughter midpoints", "internal_proteins": "Measured profiles chosen without replacement",
        "roots": "Measured root states not fixed", "scored_edges": "Cache unverified; diagnostic: 998",
        "constraints_and_status": "Final profile is removed before a break that skips scoring its pair. Valid deterministic pairings reproduce 499 scored merges (998 edges); cached optimized run not regenerated. Excluded from plots pending audit.",
    },
]


def write_inventory(heuristics: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    frame = pd.DataFrame(INVENTORY)
    if set(frame.heuristic) != set(heuristics.heuristic):
        raise AssertionError("Heuristic inventory does not cover every cached method")
    frame["cached_weight_samples"] = frame.heuristic.map(heuristics.groupby("heuristic").size())
    frame["cache_source"] = frame.heuristic.map(heuristics.groupby("heuristic").source_cache.first())
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "ce_full_tree_heuristic_inventory.csv", index=False)
    table = frame.drop(columns=["cached_weight_samples", "cache_source"]).rename(columns={
        "heuristic": "Method", "presentation": "Placement",
        "internal_positions": "Internal positions", "internal_proteins": "Internal protein states",
        "roots": "Root constraint", "scored_edges": "Scored edges",
        "constraints_and_status": "Algorithm and audit status",
    })
    columns = "".join(r">{\raggedright\arraybackslash}p{" + str(width) + "mm}"
                      for width in (30, 19, 25, 30, 30, 24, 65))
    table.to_latex(output_dir / "ce_full_tree_heuristic_inventory_rows.tex", index=False,
                   escape=True, column_format=columns)
    return frame
