"""Full C. elegans protein feature ranking from the linear baseline.

Extracted from ``elegans_protein_linear_baseline.ipynb`` (cells 2-7, 10, 14):
L2 multinomial logistic regression (C = 1, cell-type merge schema D) fit on
terminal cells of ``data/protein/aggregated_all/s3.csv``; global importance of
a protein is the maximum over classes of |coefficient|.  The notebook saved
only the top 20; this script regenerates the full ranking of all 210 proteins,
checks that its first 20 reproduce the saved ``top20_protein_names.csv`` in
order, and only then writes ``protein_feature_ranking.csv`` beside it.

Run from the repository root::

    python expression_embedding/elegans_protein_feature_ranking.py
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

BUNDLE_DIR = Path(__file__).resolve().parent
ROOT_DIR = BUNDLE_DIR.parent
DATA_DIR = ROOT_DIR / "data"
S3_PATH = DATA_DIR / "protein" / "aggregated_all" / "s3.csv"
CELL_LINEAGE_PATH = DATA_DIR / "cell_lineage.json"
CELL_TYPE_PATH = DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv"
RUN_DIR = BUNDLE_DIR / "results" / "elegans_protein_linear_baseline"
TOP20_PATH = RUN_DIR / "top20_protein_names.csv"
RANKING_PATH = RUN_DIR / "protein_feature_ranking.csv"

C_VALUE = 1.0

# Schema D (notebook cell 5).
MERGE_MAP_D = {
    "neuron": "neuron", "muscle": "muscle", "repro": "reproduction",
    "hypoderm": "epithelium", "epithelium": "epithelium",
    "sheath": "glial", "socket": "glial",
    "coelomocyte": "coelomocyte", "excretory": "excretory",
    "mesoderm": "mesoderm", "other": "other",
    "intestine": "alimentary", "valve": "alimentary",
    "marginal": "alimentary", "gland": "alimentary",
    "rectal": "alimentary",
    "tail": "programmed_death", "programmed_death": "programmed_death",
}
DEATH_LABEL = "programmed_death"


def map_names(did: str) -> str:
    if did == "P4a":
        return "Z3"
    if did == "P4p":
        return "Z2"
    if did == "P0a":
        return "AB"
    return did


def lineage_terminals() -> tuple[list[str], dict[str, list[str]]]:
    """Terminal nodes and, per internal node, its terminal descendants (cell 2)."""
    with CELL_LINEAGE_PATH.open("r", encoding="utf-8") as f:
        lineage = json.load(f)
    terminals: list[str] = []
    descendants: dict[str, list[str]] = defaultdict(list)

    def dfs(node, ancestors):
        name = map_names(node["did"])
        children = node.get("children", [])
        if not children:
            terminals.append(name)
            for ancestor in ancestors:
                descendants[ancestor].append(name)
        else:
            for child in children:
                dfs(child, ancestors + [name])

    dfs(lineage, [])
    return terminals, descendants


def terminal_cell_types(terminals: list[str]) -> dict[str, str]:
    """Schema-D cell type per terminal node (cells 3 and 5)."""
    key = pd.read_csv(CELL_TYPE_PATH)
    types: dict[str, str] = {}
    for node in terminals:
        rows = key[key["wormweb.lineage"] == node]
        if len(rows) == 0:
            continue
        labels = rows["wormweb.type"].dropna().unique()
        types[node] = labels[0] if len(labels) else DEATH_LABEL
    for node in terminals:
        types.setdefault(node, DEATH_LABEL)
    return {node: MERGE_MAP_D[t] for node, t in types.items()}


def class_order(cell_type_dict: dict[str, str]) -> list[str]:
    """Classes by descending frequency, programmed death last (cell 6)."""
    living = sorted(set(cell_type_dict.values()) - {DEATH_LABEL})
    counts = pd.Series(list(cell_type_dict.values())).value_counts()
    living = sorted(living, key=lambda c: counts.get(c, 0), reverse=True)
    return living + [DEATH_LABEL]


def fit_ranking() -> tuple[pd.DataFrame, dict]:
    terminals, descendants = lineage_terminals()
    cell_type_dict = terminal_cell_types(terminals)
    classes = class_order(cell_type_dict)
    class_index = {c: i for i, c in enumerate(classes)}

    # Cell 7: proteins x cells, drop all-zero cells, z-score each cell across
    # proteins, then standardise each protein across cells.
    protein_exp = pd.read_csv(S3_PATH, index_col=0).T.fillna(0)
    protein_exp = protein_exp[~(protein_exp == 0).all(axis=1)]
    protein_exp_zscore = protein_exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    X_full = StandardScaler().fit_transform(protein_exp_zscore.values.T.astype(np.float64))
    feature_names = np.array(protein_exp_zscore.index.tolist())
    cells = [map_names(n) for n in protein_exp_zscore.columns]
    is_terminal = np.array([len(descendants[c]) == 0 for c in cells])

    # Hard labels exist only for terminals; internal cells are not used.
    X_term = X_full[is_terminal]
    y_term = np.array([class_index[cell_type_dict[c]] for c, t in zip(cells, is_terminal) if t])

    # Cell 14.
    model = LogisticRegression(C=C_VALUE, max_iter=5000, tol=1e-4, random_state=42)
    model.fit(X_term, y_term)
    coef = model.coef_
    importance = np.abs(coef).max(axis=0)
    order = np.argsort(importance)[::-1]

    fitted_classes = [classes[i] for i in model.classes_]
    table = pd.DataFrame({
        "rank": np.arange(1, len(order) + 1),
        "protein": feature_names[order],
        "global_importance": importance[order],
        "top_class": [fitted_classes[j] for j in np.abs(coef[:, order]).argmax(axis=0)],
    })
    for i, name in enumerate(fitted_classes):
        table[f"{name}_coef"] = coef[i, order]

    info = {
        "n_cells": int(X_full.shape[0]),
        "n_proteins": int(X_full.shape[1]),
        "n_terminal": int(is_terminal.sum()),
        "classes": fitted_classes,
    }
    return table, info


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dry-run", action="store_true", help="check only; do not write")
    args = parser.parse_args(argv)

    table, info = fit_ranking()
    print(f"cells={info['n_cells']} proteins={info['n_proteins']} "
          f"terminal={info['n_terminal']} classes={len(info['classes'])}")

    saved = pd.read_csv(TOP20_PATH)["protein"].tolist()
    regenerated = table["protein"].head(len(saved)).tolist()
    if regenerated != saved:
        print("MISMATCH: regenerated top 20 differs from the saved list", file=sys.stderr)
        for r, (a, b) in enumerate(zip(saved, regenerated), start=1):
            flag = "" if a == b else "   <--"
            print(f"  {r:2d}. saved={a:22s} regenerated={b}{flag}", file=sys.stderr)
        return 1
    print(f"Top {len(saved)} reproduces {TOP20_PATH.name} exactly (same order).")

    imp = table["global_importance"].to_numpy()
    for k in (5, 10, 20, 30):
        print(f"  k={k:2d}: importance at rank k = {imp[k - 1]:.4f}, "
              f"rank k+1 = {imp[k]:.4f}, gap = {imp[k - 1] - imp[k]:.4f}")

    if args.dry_run:
        return 0
    table.to_csv(RANKING_PATH, index=False)
    print(f"Wrote {RANKING_PATH.relative_to(ROOT_DIR)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
