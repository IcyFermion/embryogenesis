"""Data loading, preprocessing, and pair sampling for mRNA-protein alignment."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure repo root is on sys.path for cross-bundle imports
BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from expression_embedding.timepoint_embedding import sublineage_split, map_names


def load_rna_protein_data(
    rna_path: str,
    protein_path: str,
    lineage_path: str,
    sublineage_depth: int = 5,
    val_fraction: float = 0.2,
    seed: int = 42,
    log_transform: bool = True,
) -> dict:
    """Load and align RNA-seq and protein embedding data.

    Args:
        log_transform: If True, apply log1p to mRNA values before L2-normalization.
            Recommended since mRNA count data is highly right-skewed.

    Returns a dict with keys:
        X_train, X_val: (n_cells, n_genes) float32, L2-normalized
        prot_train, prot_val: (n_cells, embed_dim) float32
        train_cells, val_cells: list[str]
        all_cells: list[str] — all overlapping cells
        X_all: (n_total, n_genes) float32
        prot_all: (n_total, embed_dim) float32
        gene_names: list[str]
        n_features: int
        train_indices, val_indices: indices into all_cells
    """
    # ── Load RNA-seq data ───────────────────────────────────────────────
    rna_df = pd.read_csv(rna_path, index_col=0)  # genes × cells
    rna_df = rna_df.T  # cells × genes
    rna_df.index.name = "cell_name"

    # ── Filter to shared TF orthologs (C. elegans ∩ C. briggsae) ────────
    # The shared set is frozen now; Stage 2 will reuse it for cross-species.
    shared_count = len(rna_df.columns)
    bri_path = str(Path(rna_path).parent / "c_briggsae_tf.csv")
    bri_genes = set(pd.read_csv(bri_path, index_col=0).index)
    ele_genes = set(rna_df.columns)
    shared_genes = sorted(ele_genes & bri_genes)
    excluded = sorted(ele_genes - bri_genes)
    rna_df = rna_df[shared_genes]
    if excluded:
        print(f"  Filtered to {len(shared_genes)}/{shared_count} shared orthologs "
              f"(excluded {len(excluded)} elegans-only: {', '.join(excluded[:6])}...)")

    rna_cells = set(rna_df.index)
    gene_names = list(rna_df.columns)
    n_genes = len(gene_names)
    print(f"  RNA data: {rna_df.shape[0]} cells × {n_genes} genes")

    # ── Load protein embeddings ─────────────────────────────────────────
    prot_df = pd.read_csv(protein_path, index_col=0)  # cells × embed_dim
    prot_cells = set(prot_df.index)
    embed_dim = prot_df.shape[1]
    print(f"  Protein embeddings: {prot_df.shape[0]} cells × {embed_dim}D")

    # ── Intersect cells ─────────────────────────────────────────────────
    shared_cells = sorted(rna_cells & prot_cells)
    if len(shared_cells) < 2:
        raise ValueError(
            f"Only {len(shared_cells)} cells in common between RNA and protein data"
        )
    print(f"  Shared cells: {len(shared_cells)}")

    rna_aligned = rna_df.loc[shared_cells].values.astype(np.float64)
    prot_aligned = prot_df.loc[shared_cells].values.astype(np.float32)

    # ── Check for NaN ───────────────────────────────────────────────────
    assert not np.any(np.isnan(rna_aligned)), "NaN in RNA data"
    assert not np.any(np.isnan(prot_aligned)), "NaN in protein embeddings"

    # ── Log-transform mRNA (highly skewed count data) ────────────────────
    if log_transform:
        rna_aligned = np.log1p(rna_aligned)
        print(f"  Applied log1p transform to mRNA data")

    # ── L2-normalize mRNA rows ──────────────────────────────────────────
    norms = np.linalg.norm(rna_aligned, axis=1, keepdims=True)
    zero_norm = norms[:, 0] < 1e-10
    if zero_norm.any():
        print(
            f"  Warning: {zero_norm.sum()} cells have zero expression; "
            f"using uniform vector"
        )
        rna_aligned[zero_norm] = 1.0 / np.sqrt(n_genes)
        norms[zero_norm] = 1.0
    X_all = (rna_aligned / norms).astype(np.float32)

    # ── Build sample_meta for sublineage split ──────────────────────────
    sample_meta = pd.DataFrame({"cell_name": shared_cells})

    # ── Sublineage split ────────────────────────────────────────────────
    print(f"\n  Sublineage split (depth={sublineage_depth}, val={val_fraction}):")
    train_idx, val_idx = sublineage_split(
        lineage_path, sample_meta, sublineage_depth, val_fraction, seed
    )

    train_cells = [shared_cells[i] for i in train_idx]
    val_cells = [shared_cells[i] for i in val_idx]
    print(f"  Train cells: {len(train_cells)}, Val cells: {len(val_cells)}")

    X_train = X_all[train_idx]
    X_val = X_all[val_idx]
    prot_train = prot_aligned[train_idx]
    prot_val = prot_aligned[val_idx]

    return {
        "X_train": X_train,
        "X_val": X_val,
        "prot_train": prot_train,
        "prot_val": prot_val,
        "train_cells": train_cells,
        "val_cells": val_cells,
        "all_cells": shared_cells,
        "X_all": X_all,
        "prot_all": prot_aligned,
        "gene_names": gene_names,
        "n_features": n_genes,
        "train_indices": train_idx,
        "val_indices": val_idx,
    }


def precompute_protein_distances(prot_all: np.ndarray) -> np.ndarray:
    """Compute full N×N cosine distance matrix for protein embeddings.

    L2-normalizes the protein embeddings first (they may not be unit-norm).
    Returns float32 array of shape (N, N) where d_ij = 1 - cos_sim(i, j).
    """
    norms = np.linalg.norm(prot_all, axis=1, keepdims=True)
    norms = np.clip(norms, 1e-10, None)
    prot_normed = prot_all / norms
    cos_sim = prot_normed @ prot_normed.T
    cos_dist = np.clip(1.0 - cos_sim, 0.0, 2.0)
    return cos_dist.astype(np.float32)


class PairSampler:
    """Efficient per-epoch pair sampler with bias toward near protein-space pairs.

    Pre-indexes all valid upper-triangular pairs within the split and sorts them
    by protein distance. Each call to sample() returns pair indices biased toward
    the bottom quartile of protein distances.
    """

    def __init__(
        self,
        cell_indices: np.ndarray,
        full_prot_dist: np.ndarray,
        near_fraction: float = 0.5,
        seed: int = 42,
    ):
        """
        Args:
            cell_indices: indices into the full protein distance matrix for
                          cells in this split (train or val).
            full_prot_dist: N_total×N_total precomputed cosine distance matrix.
            near_fraction: fraction of pairs to draw from bottom quartile.
            seed: random seed for reproducibility.
        """
        self.cell_indices = cell_indices
        self.full_prot_dist = full_prot_dist
        self.near_fraction = near_fraction
        self.rng = np.random.RandomState(seed)

        n = len(cell_indices)
        if n < 2:
            self.near_pairs = np.empty((0, 2), dtype=int)
            self.uniform_pairs = np.empty((0, 2), dtype=int)
            return

        # Build all upper-triangular pair indices (local)
        triu_rows, triu_cols = np.triu_indices(n, k=1)
        local_pairs = np.column_stack([triu_rows, triu_cols])  # (n_pairs, 2)

        # Map to full-matrix indices and get protein distances
        full_i = cell_indices[triu_rows]
        full_j = cell_indices[triu_cols]
        pair_dists = full_prot_dist[full_i, full_j]

        # Near pairs: bottom quartile
        q1 = np.percentile(pair_dists, 25)
        near_mask = pair_dists <= q1
        self.near_pairs = local_pairs[near_mask]  # (n_near, 2)
        self.uniform_pairs = local_pairs  # (n_total_pairs, 2)

    def sample(self, n_pairs: int) -> tuple[np.ndarray, np.ndarray]:
        """Sample n_pairs pairs.

        Returns (idx_a, idx_b) — each shape (n_pairs,), local indices into the
        split's cell array.
        """
        if len(self.uniform_pairs) == 0:
            return np.array([], dtype=int), np.array([], dtype=int)

        n_near = int(n_pairs * self.near_fraction)
        n_uniform = n_pairs - n_near

        near_idx = self.rng.choice(
            max(len(self.near_pairs), 1), size=n_near, replace=True
        )
        # Clamp to actual near pool size, then sample with replacement if needed
        if len(self.near_pairs) > 0:
            near_idx = self.rng.choice(
                len(self.near_pairs), size=n_near, replace=True
            )
            near_selected = self.near_pairs[near_idx]
        else:
            near_selected = np.empty((0, 2), dtype=int)

        uniform_idx = self.rng.choice(
            len(self.uniform_pairs), size=n_uniform, replace=True
        )
        uniform_selected = self.uniform_pairs[uniform_idx]

        if len(near_selected) > 0:
            selected = np.vstack([near_selected, uniform_selected])
        else:
            selected = uniform_selected

        self.rng.shuffle(selected)
        return selected[:, 0], selected[:, 1]
