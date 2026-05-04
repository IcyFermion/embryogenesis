"""
timepoint_embedding.py
======================
Per-timepoint expression embedding model for C. elegans protein atlas data.

Encoder: MLP that maps a single-timepoint expression vector to a fixed-dimensional
embedding. Training combines reconstruction, terminal-type classification with soft
labels, and a consecutive-timepoint smoothness loss. No lineage-relational supervision.

Train/validation split is by sub-lineage to avoid timepoint leakage.
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

mpl.rcParams["figure.dpi"] = 150

# ─────────────────────────────────────────────────────────────────────────────
# Path discovery (same pattern as other bundle scripts)
# ─────────────────────────────────────────────────────────────────────────────

BUNDLE_DIR = Path(__file__).resolve().parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = BUNDLE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class Config:
    # Paths
    h5ad_path: str = str(DATA_DIR / "protein" / "aggregated_scanpy.h5ad")
    lineage_path: str = str(DATA_DIR / "cell_lineage.json")
    cell_type_path: str = str(DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv")
    output_dir: str = str(RESULTS_DIR / "timepoint_embedding")

    # Feature set
    use_all_features: bool = False   # if True, use all 266 TFs; if False, use selected 25
    selected_tfs_path: str = str(RESULTS_DIR / "nn_selected_proteins_rev.tsv")

    # Model architecture (hidden_dims chosen automatically in main() based on feature count)
    n_features: int = 25  # set automatically from data
    hidden_dims: tuple = (64, 32)
    embed_dim: int = 32
    dropout: float = 0.1
    use_layer_norm: bool = True

    # Loss coefficients
    alpha: float = 1.0   # reconstruction weight
    beta: float = 1.0    # classification weight
    gamma: float = 0.1   # smoothness weight

    # Training
    batch_size: int = 256
    n_epochs: int = 200
    lr: float = 1e-3
    weight_decay: float = 1e-4
    patience: int = 25
    grad_clip: float = 1.0

    # Split
    sublineage_depth: int = 5
    val_fraction: float = 0.2

    # Smoothness pair sampling: how many (c,t)→(c,t+1) pairs to include per batch
    # as a fraction of batch_size (0 = none, 0.5 = half the batch)
    smoothness_pair_fraction: float = 0.5

    # Optional: per-sample label-confidence weight decay by distance from terminal time
    use_time_decay_weights: bool = False
    time_decay_tau: float = 50.0  # decay constant in timepoint units

    # Misc
    seed: int = 42
    device: str = "auto"


# ─────────────────────────────────────────────────────────────────────────────
# Name mapping (shared across the bundle)
# ─────────────────────────────────────────────────────────────────────────────

def map_names(did: str) -> str:
    if did == "P4a":
        return "Z3"
    elif did == "P4p":
        return "Z2"
    elif did == "P0a":
        return "AB"
    return did


# ─────────────────────────────────────────────────────────────────────────────
# Data loading & preprocessing
# ─────────────────────────────────────────────────────────────────────────────

def load_selected_tfs(path: str) -> list[str]:
    """Load the list of selected TF gene names (no strain suffix)."""
    tfs = pd.read_csv(path, sep="\t", header=None)[0].tolist()
    return [t.strip() for t in tfs if t.strip()]


def load_and_preprocess(config: Config):
    """Load h5ad, filter to selected TFs (or all), z-score normalize, build sample table.

    Returns
    -------
    X : ndarray (n_samples, n_features) float32
    sample_meta : DataFrame with columns [cell_name, time, orig_idx]
    feature_names : list[str]
    """
    import anndata as ad

    print(f"Loading {config.h5ad_path} ...")
    adata = ad.read_h5ad(config.h5ad_path)
    var_names = adata.var_names.tolist()

    if config.use_all_features:
        tf_names_found = var_names
        tf_indices = list(range(len(var_names)))
        print(f"  Using all {len(tf_names_found)} TFs")
    else:
        selected_tfs = load_selected_tfs(config.selected_tfs_path)
        tf_indices = []
        tf_names_found = []
        for tf in selected_tfs:
            if tf in var_names:
                tf_indices.append(var_names.index(tf))
                tf_names_found.append(tf)
            else:
                print(f"  Warning: {tf} not found in h5ad var names, skipping")
        print(f"  Using {len(tf_names_found)}/{len(selected_tfs)} selected TFs")

    # Extract expression matrix for selected TFs
    import scipy.sparse
    X_full = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.array(adata.X)
    X = X_full[:, tf_indices].astype(np.float64)

    # z-score normalize per TF across all samples
    feat_mean = X.mean(axis=0)
    feat_std = X.std(axis=0)
    feat_std[feat_std == 0] = 1.0
    X = (X - feat_mean) / feat_std
    X = np.clip(X, -5.0, 5.0).astype(np.float32)

    # Build sample metadata
    obs_df = adata.obs.copy()
    obs_df["cell_name"] = obs_df["Cell-name"].apply(map_names)
    obs_df["time"] = obs_df["Time"].values.astype(int)
    obs_df["orig_idx"] = np.arange(len(obs_df))

    # Sort by cell_name then time — critical for smoothness pair indexing
    obs_df = obs_df.sort_values(["cell_name", "time"]).reset_index(drop=True)
    X = X[obs_df["orig_idx"].values]

    # Recompute orig_idx after sorting
    obs_df["orig_idx"] = np.arange(len(obs_df))

    config.n_features = len(tf_names_found)
    config.hidden_dims = tuple(config.hidden_dims)

    print(f"  Samples: {X.shape[0]}, Features: {X.shape[1]}")
    print(f"  Unique cells: {obs_df['cell_name'].nunique()}")
    print(f"  Time range: [{obs_df['time'].min()}, {obs_df['time'].max()}]")

    return X, obs_df, tf_names_found, feat_mean, feat_std


# ─────────────────────────────────────────────────────────────────────────────
# Label construction
# ─────────────────────────────────────────────────────────────────────────────

def build_lineage_labels(lineage_path: str, cell_type_path: str):
    """Build hard (one-hot) and soft (distribution) labels for all lineage nodes.

    Reuses the logic from protein_feature_select.ipynb:
    - Terminal cells get the one-hot of their known cell type.
    - Intermediate cells get the average of their descendant terminals' one-hots.

    Returns
    -------
    cell_type_one_hot : dict str → ndarray (n_classes,) float32
    cell_types : list[str]  — ordered class names
    terminal_nodes : list[str]
    intermediate_nodes : list[str]
    descendant_list_dict : dict str → list[str]
    lineage_data : dict  — raw lineage tree
    """
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    terminal_nodes = []
    intermediate_nodes = []
    descendant_list_dict = defaultdict(list)

    def dfs(node, parent, ancestors=None):
        if ancestors is None:
            ancestors = []
        children = node.get("children", [])
        lookup_name = map_names(node["did"])
        if len(children) == 0:
            terminal_nodes.append(lookup_name)
            for ancestor in ancestors:
                descendant_list_dict[ancestor].append(lookup_name)
        else:
            intermediate_nodes.append(lookup_name)
            for child in children:
                dfs(child, node, ancestors + [lookup_name])

    dfs(lineage_data, None)

    # Load cell type assignments
    cell_type_df = pd.read_csv(cell_type_path)
    cell_type_dict = {}
    for node in terminal_nodes:
        cur = cell_type_df[cell_type_df["wormweb.lineage"] == node]
        cur_types = cur["wormweb.type"].dropna().unique()
        if len(cur_types) == 0:
            cell_type_dict[node] = "programmed_death"
        else:
            cell_type_dict[node] = cur_types[0]

    # Build ordered class list (programmed_death last)
    cell_types = sorted(set(cell_type_dict.values()))
    if "programmed_death" in cell_types:
        cell_types.remove("programmed_death")
    cell_types = sorted(cell_types, key=lambda x: sum(1 for v in cell_type_dict.values() if v == x), reverse=True)
    cell_types.append("programmed_death")

    cell_type_to_int = {ct: i for i, ct in enumerate(cell_types)}

    # One-hot for terminal nodes
    cell_type_one_hot = {}
    for node, ct in cell_type_dict.items():
        one_hot = np.zeros(len(cell_types), dtype=np.float32)
        one_hot[cell_type_to_int[ct]] = 1.0
        cell_type_one_hot[node] = one_hot

    # Soft labels for intermediate nodes: average of descendants' one-hots
    for node in intermediate_nodes:
        desc = descendant_list_dict.get(node, [])
        if len(desc) == 0:
            cell_type_one_hot[node] = np.zeros(len(cell_types), dtype=np.float32)
            continue
        summed = np.sum([cell_type_one_hot[d] for d in desc], axis=0)
        cell_type_one_hot[node] = (summed / summed.sum()).astype(np.float32)

    # Build name→type dict for diagnostic use (mapped names)
    name_to_type = {}
    for node in terminal_nodes:
        name_to_type[node] = cell_type_dict[node]
    for node in intermediate_nodes:
        # intermediate nodes: type is argmax of soft label
        soft = cell_type_one_hot[node]
        name_to_type[node] = cell_types[soft.argmax()]

    print(f"  Terminal nodes: {len(terminal_nodes)}, Intermediate: {len(intermediate_nodes)}")
    print(f"  Cell types ({len(cell_types)}): {cell_types}")

    return (cell_type_one_hot, cell_types, terminal_nodes,
            intermediate_nodes, descendant_list_dict, lineage_data, name_to_type)


def assign_labels_to_samples(sample_meta: pd.DataFrame,
                              cell_type_one_hot: dict,
                              n_classes: int,
                              config: Config):
    """Map per-sample labels from the per-cell label dictionary.

    Every timepoint of a cell inherits that cell's label.

    Returns
    -------
    y : ndarray (n_samples, n_classes) float32
    hard_mask : ndarray (n_samples,) bool  — True for terminal (hard-label) cells
    sample_weights : ndarray (n_samples,) float32
    """
    cell_names = sample_meta["cell_name"].values
    times = sample_meta["time"].values

    y = np.zeros((len(cell_names), n_classes), dtype=np.float32)
    hard_mask = np.zeros(len(cell_names), dtype=bool)
    n_descendants = np.zeros(len(cell_names), dtype=int)

    for i, cn in enumerate(cell_names):
        if cn in cell_type_one_hot:
            y[i] = cell_type_one_hot[cn]
            # Terminal cell = one-hot (exactly one class has value 1)
            is_hard = (y[i].max() >= 0.999)
            hard_mask[i] = is_hard
        else:
            y[i] = np.ones(n_classes, dtype=np.float32) / n_classes
            hard_mask[i] = False

    # Entropy-based sample weights (reuse pattern from protein_feature_select)
    y_clipped = np.clip(y, 1e-10, 1.0)
    entropies = -np.sum(y_clipped * np.log(y_clipped), axis=1)
    max_ent = np.log(n_classes)
    norm_ent = entropies / max_ent
    weights = np.exp(-3.0 * norm_ent)
    sample_weights = (weights / weights.mean()).astype(np.float32)

    # Optional: time-decay confidence weights
    if config.use_time_decay_weights:
        # For each cell, find its max time (terminal time). Distance from terminal
        # time reduces confidence.
        cell_max_t = sample_meta.groupby("cell_name")["time"].transform("max")
        dist_to_terminal = (cell_max_t - times).values.astype(np.float32)
        time_decay = np.exp(-dist_to_terminal / config.time_decay_tau)
        sample_weights = sample_weights * time_decay
        sample_weights = (sample_weights / sample_weights.mean()).astype(np.float32)

    return y, hard_mask, sample_weights


# ─────────────────────────────────────────────────────────────────────────────
# Sub-lineage split
# ─────────────────────────────────────────────────────────────────────────────

def sublineage_split(lineage_path: str, sample_meta: pd.DataFrame,
                     depth: int, val_fraction: float, seed: int):
    """Split samples by sub-lineage to avoid timepoint leakage.

    Why sub-lineage split: splitting by timepoint or by individual cell would
    leak information because different timepoints of the same cell (or closely
    related cells within the same sub-lineage) share developmental context.
    By holding out entire sub-lineages, we ensure validation cells share no
    recent common ancestor with training cells above the configured depth.
    """
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    # Collect all nodes at the given depth from root
    sublineage_roots = []

    def collect_at_depth(node, current_depth):
        lookup_name = map_names(node["did"])
        children = node.get("children", [])
        if current_depth == depth:
            sublineage_roots.append(lookup_name)
            return
        for child in children:
            collect_at_depth(child, current_depth + 1)

    collect_at_depth(lineage_data, 0)

    # For each sub-lineage root, collect all descendant cells
    descendant_cache = {}

    def collect_descendants(node, accum):
        name = map_names(node["did"])
        children = node.get("children", [])
        if len(children) == 0:
            accum.append(name)
        else:
            for child in children:
                collect_descendants(child, accum)

    def get_descendants(root_name, node):
        name = map_names(node["did"])
        if name == root_name:
            accum = []
            collect_descendants(node, accum)
            return accum
        for child in node.get("children", []):
            result = get_descendants(root_name, child)
            if result is not None:
                return result
        return None

    sublineage_cells = {}
    for root_name in sublineage_roots:
        cells = get_descendants(root_name, lineage_data)
        if cells:
            sublineage_cells[root_name] = set(cells)

    # Map each sample to its sub-lineage root
    sample_cells = sample_meta["cell_name"].values
    sample_to_sublineage = np.full(len(sample_cells), None, dtype=object)
    cell_to_sublineage = {}

    for sl_name, sl_cells in sublineage_cells.items():
        for c in sl_cells:
            cell_to_sublineage[c] = sl_name

    for i, c in enumerate(sample_cells):
        sample_to_sublineage[i] = cell_to_sublineage.get(c, None)

    # Unique sub-lineages with samples
    unique_sl = sorted(set(sl for sl in sample_to_sublineage if sl is not None))
    print(f"  Sub-lineages at depth {depth}: {len(unique_sl)} (with samples)")

    # Shuffle and split
    rng = np.random.RandomState(seed)
    rng.shuffle(unique_sl)
    n_val = max(1, int(len(unique_sl) * val_fraction))
    val_sl = set(unique_sl[:n_val])
    train_sl = set(unique_sl[n_val:])

    train_idx = np.where([sl in train_sl for sl in sample_to_sublineage])[0]
    val_idx = np.where([sl in val_sl for sl in sample_to_sublineage])[0]

    # Assign unassigned samples to train
    unassigned = np.where(np.array([sl is None for sl in sample_to_sublineage]))[0]
    if len(unassigned) > 0:
        print(f"  Note: {len(unassigned)} samples from {len(set(sample_cells[unassigned]))} "
              f"cells not assigned to any sub-lineage, adding to train")
        train_idx = np.concatenate([train_idx, unassigned])

    print(f"  Train samples: {len(train_idx)} ({len(set(sample_cells[train_idx]))} cells)")
    print(f"  Val samples:   {len(val_idx)} ({len(set(sample_cells[val_idx]))} cells)")

    return train_idx, val_idx


# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────

class TimepointEncoder(nn.Module):
    """MLP encoder: expression vector → embedding."""

    def __init__(self, n_features: int, hidden_dims: tuple, dropout: float = 0.1,
                 use_layer_norm: bool = True):
        super().__init__()
        layers = []
        in_dim = n_features
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            if use_layer_norm:
                layers.append(nn.LayerNorm(h))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
            in_dim = h
        self.encoder = nn.Sequential(*layers)
        self.embed_dim = in_dim

    def forward(self, x):
        return self.encoder(x)


class ReconstructionDecoder(nn.Module):
    """Mirror of the encoder: embedding → reconstructed input."""

    def __init__(self, n_features: int, hidden_dims: tuple, dropout: float = 0.1,
                 use_layer_norm: bool = True):
        super().__init__()
        layers = []
        in_dim = hidden_dims[-1] if hidden_dims else n_features
        for h in reversed(hidden_dims[:-1]):
            layers.append(nn.Linear(in_dim, h))
            if use_layer_norm:
                layers.append(nn.LayerNorm(h))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
            in_dim = h
        layers.append(nn.Linear(in_dim, n_features))
        self.decoder = nn.Sequential(*layers)

    def forward(self, z):
        return self.decoder(z)


class ClassifierHead(nn.Module):
    """Small classification head on the embedding."""

    def __init__(self, embed_dim: int, n_classes: int, hidden_dim: int = 32,
                 dropout: float = 0.1):
        super().__init__()
        self.head = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_classes),
        )

    def forward(self, z):
        return self.head(z)


class TimepointEmbeddingModel(nn.Module):
    """Full model: encoder + reconstruction head + classification head."""

    def __init__(self, n_features: int, n_classes: int, hidden_dims: tuple = (64, 32),
                 dropout: float = 0.1, use_layer_norm: bool = True):
        super().__init__()
        self.encoder = TimepointEncoder(n_features, hidden_dims, dropout, use_layer_norm)
        self.decoder = ReconstructionDecoder(n_features, hidden_dims, dropout, use_layer_norm)
        self.classifier = ClassifierHead(self.encoder.embed_dim, n_classes,
                                         hidden_dim=self.encoder.embed_dim, dropout=dropout)
        self.embed_dim = self.encoder.embed_dim

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        logits = self.classifier(z)
        return z, x_recon, logits


# ─────────────────────────────────────────────────────────────────────────────
# Dataset
# ─────────────────────────────────────────────────────────────────────────────

class TimepointDataset(Dataset):
    """Dataset of (cell, timepoint) samples."""

    def __init__(self, X: np.ndarray, y: np.ndarray, weights: np.ndarray,
                 next_t_idx: np.ndarray):
        self.X = torch.FloatTensor(X)
        self.y = torch.FloatTensor(y)
        self.weights = torch.FloatTensor(weights)
        self.next_t_idx = torch.LongTensor(next_t_idx)  # -1 if no next timepoint

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.weights[idx], idx, self.next_t_idx[idx]


# ─────────────────────────────────────────────────────────────────────────────
# Training
# ─────────────────────────────────────────────────────────────────────────────

def soft_cross_entropy(logits, soft_targets, weights=None):
    """Cross-entropy with soft targets and optional per-sample weights."""
    log_p = F.log_softmax(logits, dim=-1)
    per_sample = -(soft_targets * log_p).sum(dim=-1)
    if weights is not None:
        per_sample = per_sample * weights
    return per_sample.mean()


def _resolve_device(device_str: str) -> torch.device:
    if device_str == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        elif torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device_str)


def train(config: Config, X_train, y_train, w_train, next_t_train,
          X_val, y_val, w_val, next_t_val,
          hard_mask_train, hard_mask_val, n_classes: int):
    """Training loop with reconstruction + classification + smoothness losses."""

    device = _resolve_device(config.device)
    print(f"  Device: {device}")

    train_ds = TimepointDataset(X_train, y_train, w_train, next_t_train)
    val_ds = TimepointDataset(X_val, y_val, w_val, next_t_val)

    train_loader = DataLoader(train_ds, batch_size=config.batch_size, shuffle=True,
                              drop_last=False)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size * 2, shuffle=False,
                            drop_last=False)

    model = TimepointEmbeddingModel(
        n_features=config.n_features,
        n_classes=n_classes,
        hidden_dims=config.hidden_dims,
        dropout=config.dropout,
        use_layer_norm=config.use_layer_norm,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=config.lr,
                                  weight_decay=config.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,
                                                            T_max=config.n_epochs)

    history = {
        "train_total": [], "train_recon": [], "train_classify": [], "train_smooth": [],
        "train_acc": [],
        "val_total": [], "val_recon": [], "val_classify": [], "val_smooth": [],
        "val_acc": [],
    }

    best_val_loss = float("inf")
    best_state = None
    wait = 0

    # Pre-compute the full embedding lookup for smoothness pairs
    # We'll compute on-the-fly within batches for efficiency

    for epoch in range(config.n_epochs):
        # ── Train ──────────────────────────────────────────────────────────
        model.train()
        epoch_losses = {"total": 0.0, "recon": 0.0, "classify": 0.0, "smooth": 0.0}
        epoch_acc = 0.0
        n_batches = 0

        for xb, yb, wb, idx_b, next_b in train_loader:
            xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
            next_b = next_b.to(device)

            # Forward pass for all samples in batch
            z, x_recon, logits = model(xb)

            # Reconstruction loss
            loss_recon = F.mse_loss(x_recon, xb)

            # Classification loss
            loss_classify = soft_cross_entropy(logits, yb, wb)

            # Smoothness loss: for samples that have a t+1, compute embedding diff
            valid_mask = next_b >= 0
            loss_smooth = torch.tensor(0.0, device=device)
            if valid_mask.any():
                valid_idx_in_batch = torch.where(valid_mask)[0]
                next_global_idx = next_b[valid_idx_in_batch]

                # Get the t+1 samples
                x_next = torch.FloatTensor(X_train[next_global_idx.cpu().numpy()]).to(device)
                z_next = model.encoder(x_next)

                # ||z(c,t) - z(c,t+1)||^2
                z_curr = z[valid_idx_in_batch]
                loss_smooth = ((z_curr - z_next) ** 2).sum(dim=-1).mean()

            # Combined loss
            loss = (config.alpha * loss_recon +
                    config.beta * loss_classify +
                    config.gamma * loss_smooth)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
            optimizer.step()

            epoch_losses["total"] += loss.item()
            epoch_losses["recon"] += loss_recon.item()
            epoch_losses["classify"] += loss_classify.item()
            epoch_losses["smooth"] += loss_smooth.item()

            # Hard-label accuracy within batch
            b_hard = hard_mask_train[idx_b.cpu().numpy()]
            if b_hard.any():
                b_pred = logits.detach().cpu().argmax(dim=-1).numpy()
                b_true = yb.cpu().argmax(dim=-1).numpy()
                epoch_acc += (b_pred[b_hard] == b_true[b_hard]).mean()

            n_batches += 1

        scheduler.step()

        n_b = max(n_batches, 1)
        history["train_total"].append(epoch_losses["total"] / n_b)
        history["train_recon"].append(epoch_losses["recon"] / n_b)
        history["train_classify"].append(epoch_losses["classify"] / n_b)
        history["train_smooth"].append(epoch_losses["smooth"] / n_b)
        history["train_acc"].append(epoch_acc / n_b if n_batches > 0 else 0.0)

        # ── Validate ────────────────────────────────────────────────────────
        model.eval()
        val_losses = {"total": 0.0, "recon": 0.0, "classify": 0.0, "smooth": 0.0}
        val_acc = 0.0
        n_val_b = 0

        with torch.no_grad():
            for xb, yb, wb, idx_b, next_b in val_loader:
                xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
                next_b = next_b.to(device)

                z, x_recon, logits = model(xb)

                v_recon = F.mse_loss(x_recon, xb)
                v_classify = soft_cross_entropy(logits, yb, wb)

                v_smooth = torch.tensor(0.0, device=device)
                valid_mask = next_b >= 0
                if valid_mask.any():
                    valid_idx_in_batch = torch.where(valid_mask)[0]
                    next_global_idx = next_b[valid_idx_in_batch]
                    x_next = torch.FloatTensor(X_val[next_global_idx.cpu().numpy()]).to(device)
                    z_next = model.encoder(x_next)
                    z_curr = z[valid_idx_in_batch]
                    v_smooth = ((z_curr - z_next) ** 2).sum(dim=-1).mean()

                v_total = (config.alpha * v_recon +
                           config.beta * v_classify +
                           config.gamma * v_smooth)

                val_losses["total"] += v_total.item()
                val_losses["recon"] += v_recon.item()
                val_losses["classify"] += v_classify.item()
                val_losses["smooth"] += v_smooth.item()

                b_hard = hard_mask_val[idx_b.cpu().numpy()]
                if b_hard.any():
                    b_pred = logits.cpu().argmax(dim=-1).numpy()
                    b_true = yb.cpu().argmax(dim=-1).numpy()
                    val_acc += (b_pred[b_hard] == b_true[b_hard]).mean()

                n_val_b += 1

        n_vb = max(n_val_b, 1)
        history["val_total"].append(val_losses["total"] / n_vb)
        history["val_recon"].append(val_losses["recon"] / n_vb)
        history["val_classify"].append(val_losses["classify"] / n_vb)
        history["val_smooth"].append(val_losses["smooth"] / n_vb)
        history["val_acc"].append(val_acc / n_vb if n_val_b > 0 else 0.0)

        # Early stopping
        val_total = history["val_total"][-1]
        if val_total < best_val_loss:
            best_val_loss = val_total
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            wait = 0
        else:
            wait += 1

        if epoch % 20 == 0 or epoch == config.n_epochs - 1 or wait >= config.patience:
            print(f"  Epoch {epoch:3d} | "
                  f"train total={history['train_total'][-1]:.4f} "
                  f"acc={history['train_acc'][-1]:.3f} | "
                  f"val total={history['val_total'][-1]:.4f} "
                  f"acc={history['val_acc'][-1]:.3f}")

        if wait >= config.patience:
            print(f"  Early stopping at epoch {epoch}")
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model.cpu(), history


# ─────────────────────────────────────────────────────────────────────────────
# Embedding extraction
# ─────────────────────────────────────────────────────────────────────────────

@torch.no_grad()
def extract_all_embeddings(model: TimepointEmbeddingModel, X: np.ndarray,
                            device: str = "auto") -> np.ndarray:
    """Extract embedding vectors for all samples."""
    device = _resolve_device(device)
    model = model.to(device)
    model.eval()

    X_t = torch.FloatTensor(X).to(device)
    z = model.encoder(X_t)
    return z.cpu().numpy()


# ─────────────────────────────────────────────────────────────────────────────
# Sanity checks
# ─────────────────────────────────────────────────────────────────────────────

def run_sanity_checks(config, train_idx, val_idx, sample_meta, y, next_t_idx,
                       lineage_path, cell_types):
    """Run all required sanity checks before reporting results."""
    print("\n─── Sanity Checks ───")

    # 1. No lineage leakage: train and val cells are disjoint
    train_cells = set(sample_meta.iloc[train_idx]["cell_name"])
    val_cells = set(sample_meta.iloc[val_idx]["cell_name"])
    assert train_cells.isdisjoint(val_cells), \
        f"Train/val cell overlap: {len(train_cells & val_cells)} cells"
    print(f"  ✓ Train/val cell sets disjoint ({len(train_cells)} train, {len(val_cells)} val)")

    # 2. No timepoint leakage: all timepoints of val cells are in val
    for cell in val_cells:
        cell_samples = sample_meta[sample_meta["cell_name"] == cell].index.values
        assert set(cell_samples).issubset(set(val_idx)), \
            f"Cell {cell} has samples outside val set"
    print(f"  ✓ All timepoints of val cells are in val")

    # 3. Smoothness pairs are same-cell consecutive
    for i in range(len(next_t_idx)):
        if next_t_idx[i] >= 0:
            assert sample_meta.iloc[i]["cell_name"] == sample_meta.iloc[next_t_idx[i]]["cell_name"], \
                f"Smoothness pair crosses cell boundary at {i}"
            t_curr = sample_meta.iloc[i]["time"]
            t_next = sample_meta.iloc[next_t_idx[i]]["time"]
            assert t_next == t_curr + 1, \
                f"Smoothness pair not consecutive: t={t_curr} → t={t_next} at idx {i}"
    print(f"  ✓ All smoothness pairs are same-cell, consecutive timepoints "
          f"({(next_t_idx >= 0).sum()} pairs)")

    # 4. Soft labels sum to 1
    row_sums = y.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-5), \
        f"Soft labels don't sum to 1: range [{row_sums.min():.6f}, {row_sums.max():.6f}]"
    print(f"  ✓ Soft labels sum to 1 (range [{row_sums.min():.6f}, {row_sums.max():.6f}])")

    # 5. Hard-label samples exist
    n_hard = (y.max(axis=1) >= 0.999).sum()
    assert n_hard > 0, "No hard-label samples found"
    print(f"  ✓ Hard-label samples: {n_hard}/{len(y)}")

    print("─── All sanity checks passed ───\n")


# ─────────────────────────────────────────────────────────────────────────────
# Evaluation & diagnostics
# ─────────────────────────────────────────────────────────────────────────────

def plot_training_curves(history: dict, output_dir: str):
    """Save loss and accuracy curves."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    epochs = range(1, len(history["train_total"]) + 1)

    ax = axes[0]
    ax.plot(epochs, history["train_total"], label="Train", color="steelblue")
    ax.plot(epochs, history["val_total"], label="Val", color="tomato")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Total Loss")
    ax.set_title("Total Loss")
    ax.legend()

    ax = axes[1]
    ax.plot(epochs, history["train_recon"], label="Train Recon", color="steelblue",
            linestyle="-")
    ax.plot(epochs, history["val_recon"], label="Val Recon", color="tomato",
            linestyle="-")
    ax.plot(epochs, history["train_classify"], label="Train Classify",
            color="steelblue", linestyle="--")
    ax.plot(epochs, history["val_classify"], label="Val Classify",
            color="tomato", linestyle="--")
    ax.plot(epochs, history["train_smooth"], label="Train Smooth",
            color="steelblue", linestyle=":")
    ax.plot(epochs, history["val_smooth"], label="Val Smooth",
            color="tomato", linestyle=":")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss Component")
    ax.set_title("Loss Components")
    ax.legend(fontsize=7)
    ax.set_yscale("log")

    ax = axes[2]
    ax.plot(epochs, history["train_acc"], label="Train Hard Acc", color="steelblue")
    ax.plot(epochs, history["val_acc"], label="Val Hard Acc", color="tomato")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Accuracy")
    ax.set_title("Hard-Label Classification Accuracy")
    ax.legend()

    plt.tight_layout()
    path = os.path.join(output_dir, "training_curves.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved training curves → {path}")

    # Also save CSV
    csv_path = os.path.join(output_dir, "training_curves.csv")
    pd.DataFrame(history).to_csv(csv_path, index=False)
    print(f"  Saved training curves CSV → {csv_path}")


def compute_sibling_distances(embeddings: np.ndarray, sample_meta: pd.DataFrame,
                               lineage_path: str, name_to_type: dict,
                               output_dir: str):
    """Histograms of embedding distances: same-type vs different-type sibling pairs."""
    # Build parent→children map
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    children_dict = defaultdict(list)
    parent_dict = {}

    def build_family(node, parent_name=None):
        name = map_names(node["did"])
        if parent_name is not None:
            parent_dict[name] = parent_name
            children_dict[parent_name].append(name)
        for child in node.get("children", []):
            build_family(child, name)

    build_family(lineage_data)

    type_lookup = name_to_type

    # Get per-cell embedding (mean over timepoints)
    cell_emb = {}
    cell_names = sample_meta["cell_name"].values
    for i, cn in enumerate(cell_names):
        if cn not in cell_emb:
            cell_emb[cn] = []
        cell_emb[cn].append(embeddings[i])
    for cn in cell_emb:
        cell_emb[cn] = np.mean(cell_emb[cn], axis=0)

    # Find terminal sibling pairs
    terminal_set = set()
    for cn in cell_emb:
        ct = type_lookup.get(cn)
        if ct is not None and not pd.isna(ct):
            terminal_set.add(cn)

    same_pairs, diff_pairs = [], []
    for parent, children in children_dict.items():
        term_children = [c for c in children if c in terminal_set and c in cell_emb]
        for i in range(len(term_children)):
            for j in range(i + 1, len(term_children)):
                a, b = term_children[i], term_children[j]
                type_a = type_lookup.get(a)
                type_b = type_lookup.get(b)
                if type_a is None or type_b is None:
                    continue
                if type_a == type_b:
                    same_pairs.append((a, b))
                else:
                    diff_pairs.append((a, b))

    same_dists = np.array([np.linalg.norm(cell_emb[a] - cell_emb[b]) for a, b in same_pairs])
    diff_dists = np.array([np.linalg.norm(cell_emb[a] - cell_emb[b]) for a, b in diff_pairs])

    if len(same_pairs) == 0 or len(diff_pairs) == 0:
        print("  Not enough sibling pairs for distance analysis")
        return

    u_stat, p_mw = mannwhitneyu(same_dists, diff_dists, alternative="less")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    bins = np.linspace(0, max(same_dists.max(), diff_dists.max()) * 1.05, 40)
    axes[0].hist(same_dists, bins=bins, alpha=0.65, color="steelblue",
                 label=f"Same type (n={len(same_dists)})")
    axes[0].hist(diff_dists, bins=bins, alpha=0.65, color="tomato",
                 label=f"Diff type (n={len(diff_dists)})")
    axes[0].axvline(np.median(same_dists), color="steelblue", linestyle="--")
    axes[0].axvline(np.median(diff_dists), color="tomato", linestyle="--")
    axes[0].set_xlabel("L2 distance in embedding space")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Per-cell embedding: sibling distances")
    axes[0].legend()
    axes[0].text(0.97, 0.97, f"MW p={p_mw:.2e}", transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=9,
                 bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    violin_data = pd.DataFrame({
        "pair_type": ["Same type"] * len(same_dists) + ["Diff type"] * len(diff_dists),
        "distance": np.concatenate([same_dists, diff_dists]),
    })
    sns.violinplot(data=violin_data, x="pair_type", y="distance",
                   hue="pair_type", palette=["steelblue", "tomato"],
                   legend=False, ax=axes[1])
    axes[1].set_ylabel("L2 distance in embedding space")
    axes[1].set_title("Sibling embedding distances")

    plt.tight_layout()
    path = os.path.join(output_dir, "sibling_distances.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved sibling distance plot → {path}")


def compute_lineage_distance_scatter(embeddings: np.ndarray, sample_meta: pd.DataFrame,
                                      lineage_path: str, output_dir: str):
    """Scatter of embedding pairwise distance vs lineage-tree distance."""
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    # Build ancestry for each node
    node_ancestry = {}

    def dfs_ancestry(node, ancestors=None):
        if ancestors is None:
            ancestors = []
        name = map_names(node["did"])
        node_ancestry[name] = list(ancestors)
        for child in node.get("children", []):
            dfs_ancestry(child, ancestors + [name])

    dfs_ancestry(lineage_data)

    def lineage_distance(cell1, cell2):
        a1 = node_ancestry.get(cell1, [])
        a2 = node_ancestry.get(cell2, [])
        common = set(a1) & set(a2)
        if not common:
            return len(a1) + len(a2)
        lca = max(common, key=lambda x: a1.index(x))
        return (len(a1) - a1.index(lca) - 1) + (len(a2) - a2.index(lca) - 1)

    # Per-cell embeddings
    cell_emb = {}
    cell_names = sample_meta["cell_name"].values
    for i, cn in enumerate(cell_names):
        if cn not in cell_emb:
            cell_emb[cn] = []
        cell_emb[cn].append(embeddings[i])
    for cn in cell_emb:
        cell_emb[cn] = np.mean(cell_emb[cn], axis=0)

    cells = sorted(cell_emb.keys())
    n_cells = min(500, len(cells))  # subsample for pairwise computation
    rng = np.random.RandomState(42)
    sampled = rng.choice(cells, size=n_cells, replace=False)

    lineage_dists = []
    emb_dists = []
    for i in range(len(sampled)):
        for j in range(i + 1, len(sampled)):
            ld = lineage_distance(sampled[i], sampled[j])
            ed = np.linalg.norm(cell_emb[sampled[i]] - cell_emb[sampled[j]])
            lineage_dists.append(ld)
            emb_dists.append(ed)

    r_pearson, _ = pearsonr(lineage_dists, emb_dists)
    r_spearman, _ = spearmanr(lineage_dists, emb_dists)

    fig, ax = plt.subplots(figsize=(8, 6))
    lineage_emb_df = pd.DataFrame({
        "lineage_distance": lineage_dists,
        "embedding_distance": emb_dists,
    })
    sns.boxplot(x="lineage_distance", y="embedding_distance", data=lineage_emb_df,
                ax=ax, color="steelblue", fliersize=2)
    ax.set_xlabel("Lineage Tree Distance")
    ax.set_ylabel("Embedding L2 Distance")
    ax.set_title(f"Embedding distance vs lineage distance\n"
                 f"Pearson r={r_pearson:.3f}, Spearman r={r_spearman:.3f}")

    plt.tight_layout()
    path = os.path.join(output_dir, "lineage_distance_scatter.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved lineage distance scatter → {path}")


def plot_umap_embeddings(embeddings: np.ndarray, sample_meta: pd.DataFrame,
                          y_full: np.ndarray, hard_mask: np.ndarray,
                          cell_types: list, output_dir: str):
    """2D PCA of embeddings colored by terminal type and developmental time."""
    n_vis = min(10000, len(embeddings))
    rng = np.random.RandomState(42)
    vis_idx = rng.choice(len(embeddings), size=n_vis, replace=False)

    pca = PCA(n_components=2)
    emb_2d = pca.fit_transform(embeddings[vis_idx])

    times = sample_meta.iloc[vis_idx]["time"].values
    hard_vis = hard_mask[vis_idx]
    y_vis = y_full[vis_idx]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: colored by cell type (hard-labeled only)
    ax = axes[0]
    hard_idx = np.where(hard_vis)[0]
    if len(hard_idx) > 0:
        cmap = plt.get_cmap("tab20", len(cell_types))
        y_hard = y_vis[hard_idx].argmax(axis=1)
        for i, ct in enumerate(cell_types):
            mask = y_hard == i
            if mask.sum() > 0:
                ax.scatter(emb_2d[hard_idx][mask, 0], emb_2d[hard_idx][mask, 1],
                           c=[cmap(i)], s=8, alpha=0.6, label=ct)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title("PCA of embeddings – hard-labeled cells by type")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # Right: colored by developmental time
    ax = axes[1]
    sc = ax.scatter(emb_2d[:, 0], emb_2d[:, 1], c=times, cmap="viridis",
                    s=6, alpha=0.5)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title("PCA of embeddings – colored by developmental time")
    plt.colorbar(sc, ax=ax, label="Time")

    plt.tight_layout()
    path = os.path.join(output_dir, "embedding_pca.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved embedding PCA → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main(config: Config = None):
    if config is None:
        config = Config()

    # When using all features, use a wider architecture and separate output dir
    if config.use_all_features:
        config.hidden_dims = (128, 64, 32)
        config.output_dir = str(Path(config.output_dir).with_name("timepoint_embedding_all_features"))
    else:
        config.hidden_dims = (64, 32)

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Timepoint Embedding Experiment")
    print("=" * 60)
    feature_mode = "ALL 266 TFs" if config.use_all_features else "25 selected TFs"
    print(f"Feature mode: {feature_mode}")
    print(f"Output directory: {output_dir}")
    print(f"Config: alpha={config.alpha}, beta={config.beta}, gamma={config.gamma}")
    print(f"Hidden dims: {config.hidden_dims}, embed_dim: {config.embed_dim}")
    print(f"Sub-lineage depth: {config.sublineage_depth}, val frac: {config.val_fraction}")

    # ── 1. Load & preprocess data ──────────────────────────────────────────
    print("\n─── Data Loading ───")
    X, sample_meta, feature_names, feat_mean, feat_std = load_and_preprocess(config)

    # ── 2. Build labels ────────────────────────────────────────────────────
    print("\n─── Label Construction ───")
    (cell_type_one_hot, cell_types, terminal_nodes,
     intermediate_nodes, descendant_list_dict, lineage_data, name_to_type) = \
        build_lineage_labels(config.lineage_path, config.cell_type_path)

    n_classes = len(cell_types)
    y_full, hard_mask, sample_weights = assign_labels_to_samples(
        sample_meta, cell_type_one_hot, n_classes, config)

    # ── 3. Build smoothness pair index ─────────────────────────────────────
    # For each sample (c, t), find the index of (c, t+1) if it exists
    print("\n─── Building smoothness pairs ───")
    next_t_idx = np.full(len(sample_meta), -1, dtype=int)
    cell_to_indices = defaultdict(list)
    for i, (cn, t) in enumerate(zip(sample_meta["cell_name"], sample_meta["time"])):
        cell_to_indices[(cn, t)] = i

    for (cn, t), i in cell_to_indices.items():
        next_i = cell_to_indices.get((cn, t + 1), -1)
        if next_i >= 0:
            next_t_idx[i] = next_i

    n_pairs = (next_t_idx >= 0).sum()
    print(f"  Smoothness pairs: {n_pairs} / {len(sample_meta)} samples "
          f"({n_pairs / len(sample_meta) * 100:.1f}%)")

    # ── 4. Sub-lineage split ───────────────────────────────────────────────
    print("\n─── Train/Val Split ───")
    train_idx, val_idx = sublineage_split(
        config.lineage_path, sample_meta,
        config.sublineage_depth, config.val_fraction, config.seed)

    X_train, X_val = X[train_idx], X[val_idx]
    y_train, y_val = y_full[train_idx], y_full[val_idx]
    w_train, w_val = sample_weights[train_idx], sample_weights[val_idx]
    hard_train, hard_val = hard_mask[train_idx], hard_mask[val_idx]

    # Remap next_t_idx from global indices to local split indices
    def _remap_next_t(subset_idx, global_next_t):
        global_to_local = {g: l for l, g in enumerate(subset_idx)}
        local_next = np.full(len(subset_idx), -1, dtype=int)
        for local_i, global_i in enumerate(subset_idx):
            nxt = global_next_t[global_i]
            if nxt >= 0 and nxt in global_to_local:
                local_next[local_i] = global_to_local[nxt]
        return local_next

    next_t_train = _remap_next_t(train_idx, next_t_idx)
    next_t_val = _remap_next_t(val_idx, next_t_idx)

    # ── 5. Sanity checks ──────────────────────────────────────────────────
    run_sanity_checks(config, train_idx, val_idx, sample_meta, y_full,
                       next_t_idx, config.lineage_path, cell_types)

    # ── 6. Train ──────────────────────────────────────────────────────────
    print("\n─── Training ───")
    model, history = train(
        config, X_train, y_train, w_train, next_t_train,
        X_val, y_val, w_val, next_t_val,
        hard_train, hard_val, n_classes)

    # ── 7. Extract embeddings ──────────────────────────────────────────────
    print("\n─── Extracting Embeddings ───")
    embeddings = extract_all_embeddings(model, X, config.device)
    print(f"  Embeddings shape: {embeddings.shape}")

    # ── 8. Save outputs ────────────────────────────────────────────────────
    print("\n─── Saving Outputs ───")

    # 8a. Model checkpoint
    ckpt_path = output_dir / "model_checkpoint.pt"
    torch.save({
        "model_state_dict": model.state_dict(),
        "config": {k: v for k, v in config.__dict__.items()},
        "feature_names": feature_names,
        "cell_types": cell_types,
        "feat_mean": feat_mean,
        "feat_std": feat_std,
        "history": history,
    }, str(ckpt_path))
    print(f"  Saved checkpoint → {ckpt_path}")

    # 8b. Embeddings table: (cell_id, timepoint) → embedding vector
    emb_cols = [f"emb_{i}" for i in range(embeddings.shape[1])]
    emb_df = pd.DataFrame(embeddings, columns=emb_cols)
    emb_df.insert(0, "time", sample_meta["time"].values)
    emb_df.insert(0, "cell_name", sample_meta["cell_name"].values)
    emb_df.insert(0, "orig_idx", sample_meta["orig_idx"].values)
    emb_path = output_dir / "embeddings.csv"
    emb_df.to_csv(str(emb_path), index=False)
    print(f"  Saved embeddings → {emb_path} ({emb_df.shape})")

    # 8c. Per-cell embedding (mean over timepoints)
    cell_emb_means = emb_df.groupby("cell_name", observed=True)[emb_cols].mean()
    cell_emb_means_path = output_dir / "cell_embeddings_mean.csv"
    cell_emb_means.to_csv(str(cell_emb_means_path))
    print(f"  Saved per-cell mean embeddings → {cell_emb_means_path} ({cell_emb_means.shape})")

    # Also compute std
    cell_emb_stds = emb_df.groupby("cell_name", observed=True)[emb_cols].std()
    cell_emb_stds_path = output_dir / "cell_embeddings_std.csv"
    cell_emb_stds.to_csv(str(cell_emb_stds_path))
    print(f"  Saved per-cell std embeddings → {cell_emb_stds_path}")

    # 8d. Training curves
    plot_training_curves(history, str(output_dir))

    # ── 9. Diagnostic plots ────────────────────────────────────────────────
    print("\n─── Diagnostic Plots ───")

    # 9a. Sibling distances
    compute_sibling_distances(embeddings, sample_meta,
                               config.lineage_path, name_to_type,
                               str(output_dir))

    # 9b. Lineage distance scatter
    compute_lineage_distance_scatter(embeddings, sample_meta,
                                      config.lineage_path, str(output_dir))

    # 9c. PCA visualization
    plot_umap_embeddings(embeddings, sample_meta, y_full, hard_mask, cell_types,
                          str(output_dir))

    # ── 10. Final metrics ──────────────────────────────────────────────────
    print("\n─── Final Metrics ───")
    device = _resolve_device(config.device)
    model = model.to(device)
    model.eval()

    # Accuracy on hard-labeled samples
    hard_val_idx = val_idx[hard_mask[val_idx]]
    if len(hard_val_idx) > 0:
        X_hard_t = torch.FloatTensor(X[hard_val_idx]).to(device)
        y_hard_t = torch.FloatTensor(y_full[hard_val_idx]).to(device)
        with torch.no_grad():
            _, _, logits = model(X_hard_t)
            preds = logits.argmax(dim=-1).cpu().numpy()
            trues = y_hard_t.argmax(dim=-1).cpu().numpy()
            val_hard_acc = (preds == trues).mean()

        # Train hard accuracy
        hard_train_idx = train_idx[hard_mask[train_idx]]
        X_htrain_t = torch.FloatTensor(X[hard_train_idx]).to(device)
        y_htrain_t = torch.FloatTensor(y_full[hard_train_idx]).to(device)
        with torch.no_grad():
            _, _, logits = model(X_htrain_t)
            preds = logits.argmax(dim=-1).cpu().numpy()
            trues = y_htrain_t.argmax(dim=-1).cpu().numpy()
            train_hard_acc = (preds == trues).mean()

        print(f"  Train hard-label accuracy: {train_hard_acc:.4f}")
        print(f"  Val hard-label accuracy:   {val_hard_acc:.4f}")

    # Soft-label metrics on val
    X_val_t = torch.FloatTensor(X_val).to(device)
    y_val_t = torch.FloatTensor(y_val).to(device)
    with torch.no_grad():
        _, _, logits = model(X_val_t)
        probs = F.softmax(logits, dim=-1).cpu().numpy()
        expected_target_prob = (probs * y_val).sum(axis=-1).mean()
        soft_ce = -(y_val * np.log(np.clip(probs, 1e-10, 1.0))).sum(axis=-1).mean()
    print(f"  Val expected target prob:  {expected_target_prob:.4f}")
    print(f"  Val soft cross-entropy:    {soft_ce:.4f}")

    # Save metrics summary
    metrics = {
        "train_hard_accuracy": float(train_hard_acc),
        "val_hard_accuracy": float(val_hard_acc),
        "val_expected_target_prob": float(expected_target_prob),
        "val_soft_cross_entropy": float(soft_ce),
        "best_val_total_loss": float(min(history["val_total"])),
        "epochs_run": len(history["train_total"]),
        "n_train_samples": len(train_idx),
        "n_val_samples": len(val_idx),
        "n_train_cells": int(sample_meta.iloc[train_idx]["cell_name"].nunique()),
        "n_val_cells": int(sample_meta.iloc[val_idx]["cell_name"].nunique()),
        "n_features": config.n_features,
        "n_classes": n_classes,
    }
    pd.Series(metrics).to_csv(str(output_dir / "metrics.csv"), header=False)
    print(f"  Saved metrics → {output_dir / 'metrics.csv'}")

    print("\n" + "=" * 60)
    print("Done. All outputs in:", str(output_dir))
    print("=" * 60)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Timepoint embedding experiment")
    parser.add_argument("--all-features", action="store_true",
                        help="Use all 266 TFs instead of the selected 25")
    parser.add_argument("--epochs", type=int, default=None,
                        help="Override number of epochs")
    parser.add_argument("--output", type=str, default=None,
                        help="Override output directory")
    args = parser.parse_args()

    config = Config()
    if args.all_features:
        config.use_all_features = True
    if args.epochs is not None:
        config.n_epochs = args.epochs
    if args.output is not None:
        config.output_dir = args.output

    main(config)
