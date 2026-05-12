"""Three-way comparison on strict AB-branch split.

  1. Linear encoder baseline (single linear layer + L2-norm)
  2. MLP with biased near-pair sampling (original approach)
  3. MLP with cross-branch pair sampling

All evaluated on the same strict AB hold-out split.
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from expression_embedding.timepoint_embedding import map_names
from expression_embedding.rna_protein_align.config import Config
from expression_embedding.rna_protein_align.data_loader import (
    load_rna_protein_data,
    precompute_protein_distances,
    PairSampler,
)
from expression_embedding.rna_protein_align.trainer import (
    _resolve_device,
    train,
    make_fixed_val_pairs,
    pearson_on_pairs,
    compute_alignment_loss,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Strict-split utilities
# ═══════════════════════════════════════════════════════════════════════════════


def build_lineage_ancestry(lineage_path: str) -> dict:
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)
    node_ancestry = {}
    def dfs(node, ancestors=None):
        if ancestors is None:
            ancestors = []
        name = map_names(node["did"])
        node_ancestry[name] = list(ancestors)
        for child in node.get("children", []):
            dfs(child, ancestors + [name])
    dfs(lineage_data)
    return node_ancestry


def get_branch(cell_name: str) -> str:
    """Classify a cell into its major lineage branch."""
    if cell_name.startswith("AB"):
        return "AB"
    elif cell_name.startswith("MS"):
        return "MS"
    elif cell_name.startswith("E"):
        return "E"
    elif cell_name.startswith("C"):
        return "C"
    elif cell_name.startswith("D"):
        return "D"
    elif cell_name in ("Z2", "Z3"):
        return "Z"
    else:
        return "other_P1"


def strict_branch_split(all_cells, branch="AB"):
    """Hold out one branch as val, rest as train."""
    val_cells = [c for c in all_cells if get_branch(c) == branch]
    train_cells = [c for c in all_cells if get_branch(c) != branch]
    return train_cells, val_cells


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-branch pair sampler
# ═══════════════════════════════════════════════════════════════════════════════


class CrossBranchPairSampler:
    """Samples pairs where the two cells come from DIFFERENT major lineage branches.

    Optionally biases toward the nearest cross-branch pairs in protein space.
    """

    def __init__(
        self,
        cell_indices: np.ndarray,
        all_cells: list,
        full_prot_dist: np.ndarray,
        near_fraction: float = 0.5,
        seed: int = 42,
    ):
        self.cell_indices = cell_indices
        self.all_cells = all_cells
        self.full_prot_dist = full_prot_dist
        self.near_fraction = near_fraction
        self.rng = np.random.RandomState(seed)

        # Map local index → branch
        local_branches = [get_branch(all_cells[i]) for i in cell_indices]
        self.local_branches = np.array(local_branches)

        n = len(cell_indices)
        if n < 2:
            self.near_pairs = np.empty((0, 2), dtype=int)
            self.uniform_pairs = np.empty((0, 2), dtype=int)
            return

        # Build all cross-branch upper-triangular pairs
        triu_rows, triu_cols = np.triu_indices(n, k=1)
        cross_mask = self.local_branches[triu_rows] != self.local_branches[triu_cols]
        local_pairs = np.column_stack([triu_rows[cross_mask], triu_cols[cross_mask]])

        if len(local_pairs) == 0:
            # Fall back to all pairs if no cross-branch pairs exist
            local_pairs = np.column_stack([triu_rows, triu_cols])

        # Get protein distances for these pairs
        full_i = cell_indices[local_pairs[:, 0]]
        full_j = cell_indices[local_pairs[:, 1]]
        pair_dists = full_prot_dist[full_i, full_j]

        # Near pairs: bottom quartile among cross-branch pairs
        q1 = np.percentile(pair_dists, 25)
        near_mask = pair_dists <= q1
        self.near_pairs = local_pairs[near_mask]
        self.uniform_pairs = local_pairs

    def sample(self, n_pairs: int) -> tuple[np.ndarray, np.ndarray]:
        if len(self.uniform_pairs) == 0:
            return np.array([], dtype=int), np.array([], dtype=int)

        n_near = int(n_pairs * self.near_fraction)
        n_uniform = n_pairs - n_near

        if len(self.near_pairs) > 0:
            near_idx = self.rng.choice(len(self.near_pairs), size=n_near, replace=True)
            near_selected = self.near_pairs[near_idx]
        else:
            near_selected = np.empty((0, 2), dtype=int)

        uniform_idx = self.rng.choice(len(self.uniform_pairs), size=n_uniform, replace=True)
        uniform_selected = self.uniform_pairs[uniform_idx]

        if len(near_selected) > 0:
            selected = np.vstack([near_selected, uniform_selected])
        else:
            selected = uniform_selected

        self.rng.shuffle(selected)
        return selected[:, 0], selected[:, 1]


# ═══════════════════════════════════════════════════════════════════════════════
# Linear encoder
# ═══════════════════════════════════════════════════════════════════════════════


class LinearEncoder(nn.Module):
    """Single linear projection + L2 normalization. No hidden layers."""
    def __init__(self, n_features: int, embed_dim: int = 32):
        super().__init__()
        self.linear = nn.Linear(n_features, embed_dim)
        self.embed_dim = embed_dim

    def forward(self, x):
        return F.normalize(self.linear(x), p=2, dim=1)


class LinearDecoder(nn.Module):
    def __init__(self, n_features: int, embed_dim: int = 32):
        super().__init__()
        self.linear = nn.Linear(embed_dim, n_features)

    def forward(self, z):
        return self.linear(z)


class LinearModel(nn.Module):
    def __init__(self, n_features: int, embed_dim: int = 32):
        super().__init__()
        self.encoder = LinearEncoder(n_features, embed_dim)
        self.decoder = LinearDecoder(n_features, embed_dim)
        self.embed_dim = embed_dim

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return z, x_recon


# ═══════════════════════════════════════════════════════════════════════════════
# Training adaptation for custom pair samplers
# ═══════════════════════════════════════════════════════════════════════════════


def train_with_sampler(
    config, X_train, X_val, prot_train, prot_val,
    prot_dist_matrix, train_pair_sampler, val_pair_sampler,
    n_features, model,
):
    """Train with the given model and pair samplers. Returns model, history."""
    from expression_embedding.rna_protein_align.trainer import train_epoch, validate_epoch
    import copy
    from scipy.stats import pearsonr

    device = _resolve_device(config.device)
    model = model.to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.lr, weight_decay=config.weight_decay
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.n_epochs)

    X_train_t = torch.as_tensor(X_train, dtype=torch.float32, device=device)
    X_val_t = torch.as_tensor(X_val, dtype=torch.float32, device=device)

    fixed_a, fixed_b = make_fixed_val_pairs(len(X_val), config.val_n_pairs, seed=config.seed)

    history = {
        "train_total": [], "train_align": [], "train_recon": [],
        "val_total": [], "val_align": [], "val_recon": [], "val_pearson": [],
    }

    best_val_pearson = -float("inf")
    best_state = None
    wait = 0

    for epoch in range(config.n_epochs):
        train_losses = train_epoch(
            model, X_train_t, prot_dist_matrix,
            train_pair_sampler, config, optimizer, device,
        )
        scheduler.step()

        val_losses = validate_epoch(
            model, X_val_t, prot_dist_matrix,
            val_pair_sampler, config, device,
        )

        val_pearson, _, _ = pearson_on_pairs(
            model, X_val, prot_val, fixed_a, fixed_b, device,
        )

        history["train_total"].append(train_losses["total"])
        history["train_align"].append(train_losses["align"])
        history["train_recon"].append(train_losses["recon"])
        history["val_total"].append(val_losses["total"])
        history["val_align"].append(val_losses["align"])
        history["val_recon"].append(val_losses["recon"])
        history["val_pearson"].append(val_pearson)

        if val_pearson > best_val_pearson:
            best_val_pearson = val_pearson
            best_state = copy.deepcopy(model.state_dict())
            wait = 0
        else:
            wait += 1

        if epoch % 20 == 0 or epoch == config.n_epochs - 1 or wait >= config.patience:
            print(f"  Epoch {epoch:3d} | train total={train_losses['total']:.4f} "
                  f"pearson={val_pearson:.4f}")

        if wait >= config.patience:
            print(f"  Early stopping at epoch {epoch}")
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model.cpu(), history


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment runner
# ═══════════════════════════════════════════════════════════════════════════════


def run_experiment(name, config, X_train, X_val, prot_train, prot_val,
                   prot_dist, train_sampler, val_sampler, n_features, model):
    """Run one experiment, return (best_pearson, history)."""
    print(f"\n{'='*60}")
    print(f"Experiment: {name}")
    print(f"{'='*60}")
    print(f"  Train cells: {len(X_train)}, Val cells: {len(X_val)}")

    model_cpu, history = train_with_sampler(
        config, X_train, X_val, prot_train, prot_val,
        prot_dist, train_sampler, val_sampler,
        n_features, model,
    )

    best = max(history["val_pearson"])
    best_ep = history["val_pearson"].index(best)
    print(f"  Best val Pearson: {best:.4f} at epoch {best_ep}")
    return best, history


def main():
    import argparse
    from expression_embedding.rna_protein_align.model import RnaEncoderModel

    parser = argparse.ArgumentParser(description="Strict-split experiments")
    parser.add_argument("--branch", type=str, default="AB",
                        choices=["AB", "MS", "E", "C", "D"],
                        help="Branch to hold out as val")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory for results")
    args = parser.parse_args()

    # ── Config from the best sweep run ──────────────────────────────────
    config = Config()
    config.hidden_dims = (256, 128, 64)
    config.near_fraction = 0.3
    config.beta = 0.1
    config.dropout = 0.0
    config.lr = 0.002638
    config.n_pairs_per_epoch = 10000
    config.sublineage_depth = 5
    config.seed = 53002
    config.weight_decay = 5e-5
    config.n_epochs = 300
    config.patience = 50

    out_dir = Path(args.output) if args.output else \
              Path(config.output_dir).parent / f"strict_split_{args.branch}_exps"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output: {out_dir}")

    # ── Load full data (no sublineage split — we split by branch) ──────
    print("\n─── Loading data ───")
    data = load_rna_protein_data(
        rna_path=config.rna_path,
        protein_path=config.protein_emb_path,
        lineage_path=config.lineage_path,
        sublineage_depth=config.sublineage_depth,
        val_fraction=0.0,
        seed=config.seed,
        log_transform=config.log_transform,
    )

    all_cells = data["all_cells"]
    X_all = data["X_all"]
    prot_all = data["prot_all"]

    # Branch distribution
    branch_counts = defaultdict(int)
    for c in all_cells:
        branch_counts[get_branch(c)] += 1
    print(f"  Branch distribution: {dict(branch_counts)}")

    # Strict split
    train_cells, val_cells = strict_branch_split(all_cells, args.branch)
    cell_to_idx = {c: i for i, c in enumerate(all_cells)}
    train_idx = np.array([cell_to_idx[c] for c in train_cells], dtype=int)
    val_idx = np.array([cell_to_idx[c] for c in val_cells], dtype=int)

    X_train = X_all[train_idx]
    X_val = X_all[val_idx]
    prot_train = prot_all[train_idx]
    prot_val = prot_all[val_idx]

    print(f"  Strict '{args.branch}' split: {len(train_cells)} train, {len(val_cells)} val")

    prot_dist = precompute_protein_distances(prot_all)
    n_features = data["n_features"]

    # ── Experiment 1: Linear encoder ───────────────────────────────────
    linear_model = LinearModel(n_features, config.embed_dim)
    train_samp = PairSampler(train_idx, prot_dist, config.near_fraction, config.seed)
    val_samp = PairSampler(val_idx, prot_dist, 0.0, config.seed + 1)

    lin_best, lin_hist = run_experiment(
        "Linear encoder (strict split)",
        config, X_train, X_val, prot_train, prot_val,
        prot_dist, train_samp, val_samp, n_features, linear_model,
    )

    # ── Experiment 2: MLP with original biased sampling ─────────────────
    mlp_model = RnaEncoderModel(
        n_features, config.hidden_dims, config.embed_dim,
        config.dropout, config.use_layer_norm,
    )
    # Reuse same samplers (already created above)

    mlp_best, mlp_hist = run_experiment(
        "MLP + biased sampling (strict split)",
        config, X_train, X_val, prot_train, prot_val,
        prot_dist, train_samp, val_samp, n_features, mlp_model,
    )

    # ── Experiment 3: MLP with cross-branch sampling ────────────────────
    mlp_cb_model = RnaEncoderModel(
        n_features, config.hidden_dims, config.embed_dim,
        config.dropout, config.use_layer_norm,
    )
    cb_train_samp = CrossBranchPairSampler(
        train_idx, all_cells, prot_dist,
        near_fraction=config.near_fraction, seed=config.seed,
    )
    cb_val_samp = CrossBranchPairSampler(
        val_idx, all_cells, prot_dist,
        near_fraction=0.0, seed=config.seed + 1,
    )

    mlp_cb_best, mlp_cb_hist = run_experiment(
        "MLP + cross-branch sampling (strict split)",
        config, X_train, X_val, prot_train, prot_val,
        prot_dist, cb_train_samp, cb_val_samp, n_features, mlp_cb_model,
    )

    # ── Summary ────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"RESULTS SUMMARY (strict '{args.branch}' hold-out)")
    print(f"{'='*60}")
    print(f"  Linear encoder:              r = {lin_best:.4f}")
    print(f"  MLP + biased near-pair samp:  r = {mlp_best:.4f}")
    print(f"  MLP + cross-branch sampling:  r = {mlp_cb_best:.4f}")
    print(f"  Train: {len(train_cells)} cells, Val: {len(val_cells)} cells")
    print(f"{'='*60}")

    # ── Plot comparison ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(lin_hist["val_pearson"], color="gray", alpha=0.8,
            label=f"Linear encoder (best={lin_best:.3f})")
    ax.plot(mlp_hist["val_pearson"], color="tomato", alpha=0.8,
            label=f"MLP + biased sampling (best={mlp_best:.3f})")
    ax.plot(mlp_cb_hist["val_pearson"], color="steelblue", alpha=0.8,
            label=f"MLP + cross-branch sampling (best={mlp_cb_best:.3f})")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Val Pearson r")
    ax.set_title(f"Strict '{args.branch}' hold-out — three-way comparison")
    ax.legend()
    plt.tight_layout()
    comp_path = out_dir / "three_way_comparison.png"
    fig.savefig(str(comp_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {comp_path}")

    # Save histories
    import pandas as pd
    pd.DataFrame({
        "epoch": range(len(lin_hist["val_pearson"])),
        "linear": lin_hist["val_pearson"],
    }).to_csv(out_dir / "linear_history.csv", index=False)
    pd.DataFrame({
        "epoch": range(len(mlp_hist["val_pearson"])),
        "mlp_biased": mlp_hist["val_pearson"],
    }).to_csv(out_dir / "mlp_biased_history.csv", index=False)
    pd.DataFrame({
        "epoch": range(len(mlp_cb_hist["val_pearson"])),
        "mlp_cross_branch": mlp_cb_hist["val_pearson"],
    }).to_csv(out_dir / "mlp_cross_branch_history.csv", index=False)

    # Save summary
    with open(out_dir / "summary.txt", "w") as f:
        f.write(f"Strict '{args.branch}' hold-out\n")
        f.write(f"Train: {len(train_cells)}, Val: {len(val_cells)}\n")
        f.write(f"Linear encoder:              {lin_best:.4f}\n")
        f.write(f"MLP + biased near-pair:      {mlp_best:.4f}\n")
        f.write(f"MLP + cross-branch sampling: {mlp_cb_best:.4f}\n")

    print(f"\nAll results saved to {out_dir}")


if __name__ == "__main__":
    main()
