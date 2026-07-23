"""Data loader for terminal-only Pareto front analysis.

Loads C. elegans and C. briggsae tracking data, protein/RNA expression matrices,
lineage tree, and feature selection lists. All paths are relative to the repo root.
"""

import json
import numpy as np
import pandas as pd
import os as _os

# Repo root: two levels up from this file (terminal_pareto/data_loader.py)
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
def _p(relpath):
    """Resolve a data path relative to the repo root."""
    return _os.path.join(_REPO_ROOT, relpath)


# ── Constants ──
T_CE = 255          # C. elegans full cutoff (subtree & cell type analyses)
T_CE_EARLY = 230    # C. elegans early cutoff (cross-species parity with old briggsae)
T_CB = 155          # C. briggsae old cutoff (original Nadin-lab file)
T_CB_NEW = 143      # C. briggsae new cutoff (stage-matched to elegans T=255, ~293 terminals vs 299)

CE_TRACKS_PATH = _p("data/embryo1/tracks.txt")
CB_TRACKS_PATH = _p("data/c_briggsae/CD140715HLH1cbp1.csv")
CB_NEW_PATH = _p("data/c_briggsae/yiming/CD240731cbhis72p1.csv")

# ── Replicate file lists ──
CE_REPLICATES = [
    ('embryo1', _p("data/embryo1/tracks.txt"), 255),
    ('embryo2', _p("data/embryo2/tracks.txt"), 247),
    ('embryo3', _p("data/embryo3/tracks.txt"), 226),
]

CB_REPLICATES_OLD = [
    ('140715HLH1cbp1', _p("data/c_briggsae/CD140715HLH1cbp1.csv"), 155),
    ('140711HLH1cbp2', _p("data/c_briggsae/CD140711HLH1cbp2.csv"), 155),
    ('140711HLH1cbp3', _p("data/c_briggsae/CD140711HLH1cbp3.csv"), 150),
    ('140712HLH1cbp1', _p("data/c_briggsae/CD140712HLH1cbp1.csv"), 150),
    ('140712HLH1cbp2', _p("data/c_briggsae/CD140712HLH1cbp2.csv"), 150),
]

CB_REPLICATES_NEW = [
    ('240731cbhis72p1 (she-1)', _p("data/c_briggsae/yiming/CD240731cbhis72p1.csv"), 143),
    ('240731cbhis72p2 (she-1)', _p("data/c_briggsae/yiming/CD240731cbhis72p2.csv"), 146),
    ('240731cbhis72p3 (she-1)', _p("data/c_briggsae/yiming/CD240731cbhis72p3.csv"), 143),
    ('241202cbhis72p1 (she-1)', _p("data/c_briggsae/yiming/CD241202cbhis72p1.csv"), 138),
    ('241202cbhis72p2 (she-1)', _p("data/c_briggsae/yiming/CD241202cbhis72p2.csv"), 149),
    ('241202cbhis72p4 (she-1)', _p("data/c_briggsae/yiming/CD241202cbhis72p4.csv"), 158),
    ('210519ZZY0874p1 (AF16)',  _p("data/c_briggsae/yiming/CD210519ZZY0874p1.csv"), 148),
    ('210519ZZY0874p5 (AF16)',  _p("data/c_briggsae/yiming/CD210519ZZY0874p5.csv"), 156),
]


# ── Utility ──

def load_json(path):
    with open(path) as f:
        return json.load(f)


def map_names(did):
    if did == "P4a": return "Z3"
    elif did == "P4p": return "Z2"
    elif did == "P0a": return "AB"
    else: return did


# ── C. elegans tracking ──

def load_elegans_tracking(tcut, path=CE_TRACKS_PATH):
    """Load C. elegans 4D tracking, filter by time <= tcut.
    Positions: [z, y, x] * 0.1625 (µm conversion).
    Returns (xyz_map: dict[name -> np.array], valid_names: list[str])."""
    ts_ce = pd.read_csv(path, sep="\t")
    ts_ce = ts_ce[ts_ce["t"] <= tcut]
    xyz_map = {}
    valid = []
    for name in ts_ce["name"].unique():
        tp = ts_ce[ts_ce['name'] == name]["t"].values
        if len(tp) == 1 and tp[0] == tcut:
            continue  # cell only appears at cutoff — not tracked through development
        xyz_map[name] = ts_ce[ts_ce['name'] == name].values[-1][1:4] * 0.1625
        valid.append(name)
    return xyz_map, valid


# ── C. briggsae tracking ──

def load_briggsae_tracking(tcut, path=CB_TRACKS_PATH):
    """Load C. briggsae CSV tracking, filter by time <= tcut.
    Positions: [z, x, y] from columns 8:11.
    Returns (xyz_map: dict[name -> np.array], valid_names: list[str])."""
    ts_cb = pd.read_csv(path)
    ts_cb = ts_cb[ts_cb["time"] <= tcut]
    xyz_map = {}
    valid = []
    for name in ts_cb["cell"].unique():
        tp = ts_cb[ts_cb['cell'] == name]["time"].values
        if len(tp) == 1 and tp[0] == tcut:
            continue
        xyz_map[name] = ts_cb[ts_cb['cell'] == name].values[-1][8:11]
        valid.append(name)
    return xyz_map, valid


# ── Expression data ──

def load_protein_expression(path=None):
    if path is None: path = _p("data/protein/aggregated_all/s3_zscore.csv")
    return pd.read_csv(path, index_col=0).T


def load_ce_rna(path=None):
    if path is None: path = _p("data/c_briggsae/science.adu8249/c_elegans_tf.csv")
    return pd.read_csv(path, index_col=0).T


def load_cb_rna(path=None):
    if path is None: path = _p("data/c_briggsae/science.adu8249/c_briggsae_tf.csv")
    return pd.read_csv(path, index_col=0).T


def load_prot_sel(path=None):
    if path is None: path = _p("expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv")
    return pd.read_csv(path)['protein'].tolist()


def load_rna_sel(path=None):
    if path is None: path = _p("expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv")
    return pd.read_csv(path, names=['tf'], sep='\t')['tf'].tolist()


def load_apoptotic(path=None):
    if path is None: path = _p("data/apoptotic_cells.txt")
    return pd.read_csv(path, sep='\t', header=None)[0].tolist()


def load_cell_type_map(path=None):
    if path is None: path = _p("data/2023-06-29_entropy_cell_key_V2.csv")
    return pd.read_csv(path)


# ── Lineage tree traversal ──

def find_node(root, target):
    """DFS through lineage tree to find a node by its did (after map_names)."""
    t = map_names(target)

    def search(n):
        if map_names(n.get("did", "")) == t:
            return n
        for c in n.get("children", []):
            r = search(c)
            if r is not None:
                return r
        return None
    return search(root)


def collect_terminals(root, valid_names, subtree=None):
    """Collect terminal (leaf) cells and their parents from the lineage tree.
    Only includes cells where both child AND parent are in valid_names.
    Optionally restricts to a named subtree root.

    Returns (terminal_nodes: list[str], terminal_parents: list[str])."""
    start = find_node(root, subtree) if subtree else root
    tn, tp = [], []

    def dfs(node, parent):
        children = node.get("children", [])
        ln = map_names(node["did"])
        if len(children) == 0:
            pln = map_names(parent['did']) if parent else None
            if ln in valid_names and pln in valid_names:
                tn.append(ln)
                tp.append(pln)
        else:
            for c in children:
                dfs(c, node)

    dfs(start, None)
    return tn, tp


def collect_all_subtrees(root, valid_names, min_cells=10):
    """Enumerate all lineage subtrees with >= min_cells terminal cells.
    Returns list of (subtree_name, [(cell, parent), ...]) sorted by size descending."""
    terminal_descendants = {}

    def postorder_collect(node, parent=None):
        children = node.get("children", [])
        lookup_name = map_names(node["did"])
        my_terms = []
        if len(children) == 0:
            p_lookup_name = map_names(parent['did']) if parent else None
            if lookup_name in valid_names and p_lookup_name in valid_names:
                my_terms.append((lookup_name, p_lookup_name))
        else:
            for child in children:
                my_terms.extend(postorder_collect(child, node))
        terminal_descendants[lookup_name] = my_terms
        return my_terms

    postorder_collect(root)
    qualifying = [(n, t) for n, t in terminal_descendants.items() if len(t) >= min_cells]
    qualifying.sort(key=lambda x: -len(x[1]))
    return qualifying
