"""Data loading, preprocessing, cell type labels, sublineage split, batch sampling."""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler

# Import utilities from sibling modules
_here = Path(__file__).resolve().parents[0]
_bundle = _here.parent
if str(_bundle) not in sys.path:
    sys.path.insert(0, str(_bundle))

from timepoint_embedding import map_names, sublineage_split

# ---------------------------------------------------------------------------
# Schema D — cell type merging
# ---------------------------------------------------------------------------
SCHEMA_D = {
    "neuron": "neuron",
    "muscle": "muscle",
    "repro": "reproduction",
    "hypoderm": "epithelium",
    "epithelium": "epithelium",
    "sheath": "glial",
    "socket": "glial",
    "intestine": "alimentary",
    "valve": "alimentary",
    "marginal": "alimentary",
    "gland": "alimentary",
    "rectal": "alimentary",
    "coelomocyte": "coelomocyte",
    "excretory": "excretory",
    "mesoderm": "mesoderm",
    "other": "other",
    "tail": "programmed_death",
    "programmed_death": "programmed_death",
}


def _apply_schema_d(raw_type):
    return SCHEMA_D.get(raw_type, raw_type)


# ---------------------------------------------------------------------------
# Cell type label construction
# ---------------------------------------------------------------------------


def build_cell_type_labels(cell_names, lineage_path, cell_type_path=None):
    """Build soft labels for all cells following the notebook convention.

    Parameters
    ----------
    cell_names : list[str]
        All cell names in the dataset (both species combined).
    lineage_path : str
        Path to cell_lineage.json.
    cell_type_path : str or None
        Path to 2023-06-29_entropy_cell_key_V2.csv.  If None, types come
        exclusively from the lineage tree's ``data.type`` field.

    Returns
    -------
    y : np.ndarray (n_cells, n_classes) float32 — soft labels summing to 1
    hard_mask : np.ndarray (n_cells,) bool — True for effectively one-hot labels
    sample_weights : np.ndarray (n_cells,) float32 — entropy-based, mean=1
    class_names : list[str] — ordered class names (programmed_death last)
    """
    # 1. Load and flatten lineage tree
    with open(lineage_path) as f:
        root = json.load(f)

    terminal_nodes = []
    intermediate_nodes = []
    descendant_list_dict = defaultdict(list)

    def dfs(node, ancestors):
        children = node.get("children", [])
        lookup = map_names(node["did"])
        if not children:
            terminal_nodes.append(lookup)
            for anc in ancestors:
                descendant_list_dict[anc].append(lookup)
        else:
            intermediate_nodes.append(lookup)
            for child in children:
                dfs(child, ancestors + [lookup])

    dfs(root, [])

    # 2. Terminal type lookup — try lineage tree's own type annotations first
    # then fall back to the external CSV for terminal nodes not annotated in the tree.
    all_nodes_flat = {}

    def _flatten(node):
        did = node.get("did", "")
        if did:
            all_nodes_flat[map_names(did)] = node
        for child in node.get("children", []):
            _flatten(child)

    _flatten(root)

    # Build raw type dict from lineage tree's data.type
    terminal_type_raw = {}
    for tn in terminal_nodes:
        raw_t = (all_nodes_flat.get(tn, {}).get("data") or {}).get("type", "")
        terminal_type_raw[tn] = raw_t if raw_t else None

    # Fill gaps from external CSV if provided
    if cell_type_path is not None:
        cell_type_df = pd.read_csv(cell_type_path)
        for tn in terminal_nodes:
            if terminal_type_raw[tn] is not None:
                continue
            cur = cell_type_df[cell_type_df["wormweb.lineage"] == tn]
            cur_types = cur["wormweb.type"].dropna().unique()
            if len(cur_types) == 0:
                terminal_type_raw[tn] = "programmed_death"
            else:
                terminal_type_raw[tn] = cur_types[0]
    else:
        for tn in terminal_nodes:
            if terminal_type_raw[tn] is None:
                terminal_type_raw[tn] = "programmed_death"

    # 3. Apply Schema D merging
    terminal_type_merged = {tn: _apply_schema_d(t) for tn, t in terminal_type_raw.items()}

    # 4. Build class list (programmed_death always last)
    raw_class_set = set(terminal_type_merged.values())
    class_names = sorted(raw_class_set - {"programmed_death"})
    class_names.append("programmed_death")
    n_classes = len(class_names)
    class_to_idx = {c: i for i, c in enumerate(class_names)}

    # 5. One-hot for terminals
    cell_type_one_hot = {}
    for tn in terminal_nodes:
        vec = np.zeros(n_classes, dtype=np.float32)
        merged_t = terminal_type_merged[tn]
        if merged_t in class_to_idx:
            vec[class_to_idx[merged_t]] = 1.0
        cell_type_one_hot[tn] = vec

    # 6. Soft labels for intermediates
    for node in intermediate_nodes:
        desc = descendant_list_dict.get(node, [])
        if not desc:
            cell_type_one_hot[node] = np.zeros(n_classes, dtype=np.float32)
            continue
        summed = np.sum([cell_type_one_hot.get(d, np.zeros(n_classes, dtype=np.float32)) for d in desc], axis=0)
        total = summed.sum()
        cell_type_one_hot[node] = (summed / total).astype(np.float32) if total > 0 else summed

    # 7. Assemble labels for requested cell_names
    n_cells = len(cell_names)
    y = np.zeros((n_cells, n_classes), dtype=np.float32)
    for i, name in enumerate(cell_names):
        lbl = cell_type_one_hot.get(map_names(name))
        if lbl is None:
            y[i] = np.ones(n_classes, dtype=np.float32) / n_classes  # uniform fallback
        else:
            y[i] = lbl

    # 8. Hard mask (one-hot or near one-hot)
    hard_mask = y.max(axis=1) >= 0.999

    # 9. Entropy-based sample weights
    sample_weights = _entropy_based_weights(y, alpha=3.0)

    return y, hard_mask, sample_weights, class_names


def _entropy_based_weights(y_soft, alpha=3.0):
    y_clipped = np.clip(y_soft, 1e-10, 1.0)
    entropies = -np.sum(y_clipped * np.log(y_clipped), axis=1)
    max_ent = np.log(y_soft.shape[1])
    norm_ent = entropies / max_ent
    weights = np.exp(-alpha * norm_ent)
    weights = weights / weights.mean()
    return weights.astype(np.float32)


# ---------------------------------------------------------------------------
# Data loading & preprocessing
# ---------------------------------------------------------------------------


def load_cross_species_data(config):
    """Load and preprocess cross-species RNA data, build labels, split.

    Returns
    -------
    data : dict with keys:
        X_train, X_val, y_train, y_val, species_train, species_val,
        sample_weights_train, sample_weights_val, hard_mask_train, hard_mask_val,
        cell_names_train, cell_names_val, class_names, n_features, gene_names
    """
    rng = np.random.default_rng(config.seed)

    # ---- 1. Load RNA data ----
    ele_raw = pd.read_csv(config.ele_rna_path, index_col=0).T  # cells × genes
    bri_raw = pd.read_csv(config.bri_rna_path, index_col=0).T

    # ---- 2. Filter to shared TF orthologs ----
    shared_genes = sorted(set(ele_raw.columns) & set(bri_raw.columns))
    ele_df = ele_raw[shared_genes].copy()
    bri_df = bri_raw[shared_genes].copy()
    gene_names = np.array(shared_genes)
    n_features = len(gene_names)

    n_ele = len(ele_df)
    n_bri = len(bri_df)
    print(f"Elegans: {n_ele} cells × {n_features} TFs")
    print(f"Briggsae: {n_bri} cells × {n_features} TFs")

    # ---- 3. Preprocess ----
    ele_vals = ele_df.values.astype(np.float64)
    bri_vals = bri_df.values.astype(np.float64)

    if config.log_transform:
        ele_vals = np.log1p(ele_vals)
        bri_vals = np.log1p(bri_vals)

    # Per-cell L2 normalization
    for vals in (ele_vals, bri_vals):
        norms = np.linalg.norm(vals, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        vals /= norms

    # Per-species z-score
    scaler_ele = StandardScaler()
    scaler_bri = StandardScaler()
    ele_scaled = scaler_ele.fit_transform(ele_vals).astype(np.float32)
    bri_scaled = scaler_bri.fit_transform(bri_vals).astype(np.float32)

    # Verify per-species z-score
    assert np.allclose(ele_scaled.mean(axis=0), 0, atol=1e-5), "Elegans z-score failed"
    assert np.allclose(bri_scaled.mean(axis=0), 0, atol=1e-5), "Briggsae z-score failed"

    # ---- 4. Stack ----
    X = np.vstack([ele_scaled, bri_scaled])  # (2126, 237)
    species = np.array([0] * n_ele + [1] * n_bri, dtype=np.int64)
    cell_names_all = np.array(list(ele_df.index) + list(bri_df.index))

    # ---- 5. Build cell type labels ----
    y, hard_mask, sample_weights, class_names = build_cell_type_labels(
        cell_names_all.tolist(), config.lineage_path, config.cell_type_path
    )
    n_classes = len(class_names)
    print(f"Classes ({n_classes}): {class_names}")
    print(f"Hard (terminal) cells: {hard_mask.sum()} / {len(cell_names_all)}")

    # Exclude programmed_death terminals from training if configured
    death_idx = class_names.index("programmed_death") if "programmed_death" in class_names else -1
    train_mask = np.ones(len(cell_names_all), dtype=bool)
    if config.exclude_dead_cells and death_idx >= 0:
        death_terminal = hard_mask & (y.argmax(axis=1) == death_idx)
        train_mask &= ~death_terminal
        print(f"Excluded {death_terminal.sum()} programmed_death terminal cells from training")

    # ---- 6. Sublineage split (shared across species) ----
    # Cell names are identical in both species, so split once and apply to both.
    ele_idx = np.arange(n_ele)
    bri_idx = np.arange(n_ele, n_ele + n_bri)

    shared_meta = pd.DataFrame({"cell_name": ele_df.index.tolist()})
    train_loc, val_loc = sublineage_split(
        config.lineage_path, shared_meta, config.sublineage_depth,
        config.val_fraction, config.seed
    )

    ele_train_global = ele_idx[train_loc]
    ele_val_global = ele_idx[val_loc]
    bri_train_global = bri_idx[train_loc]
    bri_val_global = bri_idx[val_loc]

    train_idx = np.sort(np.concatenate([ele_train_global, bri_train_global]))
    val_idx = np.sort(np.concatenate([ele_val_global, bri_val_global]))

    # Sanity: no train/val overlap
    assert len(set(train_idx) & set(val_idx)) == 0, "Train/val overlap detected"

    print(f"Train: {len(train_idx)} ({np.sum(species[train_idx]==0)} ele + {np.sum(species[train_idx]==1)} bri)")
    print(f"Val:   {len(val_idx)} ({np.sum(species[val_idx]==0)} ele + {np.sum(species[val_idx]==1)} bri)")

    # ---- 7. Apply train_mask (exclude_dead_cells) only to training set ----
    train_idx = train_idx[train_mask[train_idx]]

    # ---- 8. Convert to tensors ----
    X_train_t = torch.FloatTensor(X[train_idx])
    X_val_t = torch.FloatTensor(X[val_idx])
    y_train_t = torch.FloatTensor(y[train_idx])
    y_val_t = torch.FloatTensor(y[val_idx])
    species_train = species[train_idx]
    species_val = species[val_idx]
    sw_train = torch.FloatTensor(sample_weights[train_idx])
    sw_val = torch.FloatTensor(sample_weights[val_idx])
    hm_train = hard_mask[train_idx]
    hm_val = hard_mask[val_idx]

    return {
        "X_train": X_train_t,
        "X_val": X_val_t,
        "y_train": y_train_t,
        "y_val": y_val_t,
        "species_train": species_train,
        "species_val": species_val,
        "sample_weights_train": sw_train,
        "sample_weights_val": sw_val,
        "hard_mask_train": hm_train,
        "hard_mask_val": hm_val,
        "train_idx": train_idx,
        "val_idx": val_idx,
        "cell_names_all": cell_names_all,
        "species_all": species,
        "class_names": class_names,
        "n_features": n_features,
        "gene_names": gene_names,
        "X_all": X,
        "y_all": y,
        "hard_mask_all": hard_mask,
        "sample_weights_all": sample_weights,
    }


# ---------------------------------------------------------------------------
# Mixed-species batch sampler
# ---------------------------------------------------------------------------


class MixedSpeciesBatchSampler:
    """Yields batches with both species in proportion to their dataset sizes.

    Each batch contains roughly the same species ratio as the full dataset.
    """

    def __init__(self, species, batch_size, shuffle=True, seed=42):
        self.species = np.asarray(species)
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.rng = np.random.default_rng(seed)

        self.ele_idx = np.where(self.species == 0)[0]
        self.bri_idx = np.where(self.species == 1)[0]
        self.n_ele = len(self.ele_idx)
        self.n_bri = len(self.bri_idx)

        self.ele_per_batch = max(1, int(batch_size * self.n_ele / len(species)))
        self.bri_per_batch = batch_size - self.ele_per_batch
        self.batches_per_epoch = max(self.n_ele // self.ele_per_batch,
                                     self.n_bri // self.bri_per_batch)

    def __iter__(self):
        ele_perm = self.rng.permutation(self.n_ele) if self.shuffle else np.arange(self.n_ele)
        bri_perm = self.rng.permutation(self.n_bri) if self.shuffle else np.arange(self.n_bri)

        ele_ptr = 0
        bri_ptr = 0

        while ele_ptr + self.ele_per_batch <= self.n_ele and bri_ptr + self.bri_per_batch <= self.n_bri:
            batch_ele = self.ele_idx[ele_perm[ele_ptr:ele_ptr + self.ele_per_batch]]
            batch_bri = self.bri_idx[bri_perm[bri_ptr:bri_ptr + self.bri_per_batch]]
            batch = np.concatenate([batch_ele, batch_bri])
            self.rng.shuffle(batch)
            yield batch
            ele_ptr += self.ele_per_batch
            bri_ptr += self.bri_per_batch

    def __len__(self):
        return min(self.n_ele // self.ele_per_batch, self.n_bri // self.bri_per_batch)
