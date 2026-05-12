"""Sanity checks, metrics, and diagnostic plots for mRNA-protein alignment."""

import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.decomposition import PCA

mpl.rcParams["figure.dpi"] = 150

# Import for name mapping reuse
import sys
BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
from expression_embedding.timepoint_embedding import map_names


# ═══════════════════════════════════════════════════════════════════════════════
# Sanity checks
# ═══════════════════════════════════════════════════════════════════════════════


def run_sanity_checks(
    X_train: np.ndarray,
    X_val: np.ndarray,
    prot_train: np.ndarray,
    prot_val: np.ndarray,
    train_cells: list,
    val_cells: list,
    model: torch.nn.Module,
    device: torch.device,
    config,
) -> dict:
    """Run all pre/post-training sanity checks. Returns baseline metrics dict."""
    print("\n" + "=" * 60)
    print("Sanity Checks")
    print("=" * 60)

    results = {}

    # 1. Train/val cell sets disjoint
    train_set = set(train_cells)
    val_set = set(val_cells)
    assert train_set.isdisjoint(val_set), (
        f"Train/val overlap: {len(train_set & val_set)} cells"
    )
    print(f"  [PASS] Train/val cells disjoint ({len(train_set)} train, {len(val_set)} val)")

    # 2. Encoder outputs L2-normalized
    model.eval()
    with torch.no_grad():
        X_t = torch.as_tensor(X_train[:min(64, len(X_train))], dtype=torch.float32, device=device)
        z = model.encoder(X_t)
        norms = z.norm(p=2, dim=1).cpu().numpy()
    assert np.allclose(norms, 1.0, atol=1e-4), (
        f"Encoder output norms not 1: [{norms.min():.6f}, {norms.max():.6f}]"
    )
    print(f"  [PASS] Encoder outputs L2-normalized (norms ~1.0, max dev={abs(norms - 1).max():.2e})")

    # 3. Linear probe baseline: raw mRNA → protein cosine distance R²
    print("  Computing linear probe baseline ...")
    r2_baseline = _linear_probe_baseline(X_train, prot_train)
    results["linear_probe_r2"] = r2_baseline
    print(f"  Linear probe R² = {r2_baseline:.4f} (expect ~0.33)")

    if r2_baseline < 0.2:
        print("  [WARN] Linear probe R² is low — check data normalization")
    elif r2_baseline > 0.5:
        print("  [WARN] Linear probe R² unusually high — possible data leakage?")

    print("=" * 60)
    return results


def _linear_probe_baseline(X: np.ndarray, prot: np.ndarray) -> float:
    """Compute R² between raw mRNA cosine distances and protein cosine distances.

    Uses a subset of pairs for efficiency (max ~50K pairs).
    """
    n = len(X)
    # Sample pairs
    max_pairs = 50000
    rng = np.random.RandomState(42)
    if n * (n - 1) // 2 > max_pairs:
        idx_a = rng.randint(0, n, size=max_pairs * 2)
        idx_b = rng.randint(0, n, size=max_pairs * 2)
        same = idx_a == idx_b
        idx_a = idx_a[~same][:max_pairs]
        idx_b = idx_b[~same][:max_pairs]
    else:
        triu_rows, triu_cols = np.triu_indices(n, k=1)
        idx_a = triu_rows
        idx_b = triu_cols

    # Raw mRNA cosine distance (X is already L2-normalized)
    cos_sim_rna = (X[idx_a] * X[idx_b]).sum(axis=1)
    cos_dist_rna = 1.0 - cos_sim_rna

    # Protein cosine distance (normalize first)
    prot_norms = np.linalg.norm(prot, axis=1, keepdims=True)
    prot_norms = np.clip(prot_norms, 1e-10, None)
    prot_n = prot / prot_norms
    cos_sim_prot = (prot_n[idx_a] * prot_n[idx_b]).sum(axis=1)
    cos_dist_prot = np.clip(1.0 - cos_sim_prot, 0.0, 2.0)

    r, _ = pearsonr(cos_dist_rna, cos_dist_prot)
    return r ** 2


# ═══════════════════════════════════════════════════════════════════════════════
# Diagnostic plots
# ═══════════════════════════════════════════════════════════════════════════════


def plot_training_curves(history: dict, output_dir: str):
    """Save loss and Pearson correlation curves."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    epochs = range(1, len(history["train_total"]) + 1)

    ax = axes[0]
    ax.plot(epochs, history["train_total"], label="Train total", color="steelblue")
    ax.plot(epochs, history["val_total"], label="Val total", color="tomato")
    ax.plot(
        epochs, history["train_align"], label="Train align",
        color="steelblue", linestyle="--", alpha=0.6,
    )
    ax.plot(
        epochs, history["val_align"], label="Val align",
        color="tomato", linestyle="--", alpha=0.6,
    )
    ax.plot(
        epochs, history["train_recon"], label="Train recon",
        color="steelblue", linestyle=":", alpha=0.6,
    )
    ax.plot(
        epochs, history["val_recon"], label="Val recon",
        color="tomato", linestyle=":", alpha=0.6,
    )
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title("Loss Components")
    ax.legend(fontsize=7)

    ax = axes[1]
    ax.plot(epochs, history["val_pearson"], color="darkgreen", marker="o", markersize=2)
    ax.axhline(y=0.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Pearson r")
    ax.set_title("Validation Alignment Correlation (headline metric)")

    plt.tight_layout()
    path = os.path.join(output_dir, "training_curves.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    csv_path = os.path.join(output_dir, "training_curves.csv")
    pd.DataFrame(history).to_csv(csv_path, index=False)
    print(f"  Saved training curves → {path}")


def plot_distance_scatter(
    mrna_dists: np.ndarray,
    protein_dists: np.ndarray,
    output_dir: str,
    title_suffix: str = "val",
):
    """Scatter of protein cosine distance vs mRNA-encoder cosine distance."""
    r, _ = pearsonr(mrna_dists, protein_dists)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        protein_dists, mrna_dists, s=3, alpha=0.3, color="steelblue", edgecolors="none"
    )
    lims = [0, max(protein_dists.max(), mrna_dists.max()) * 1.05]
    ax.plot(lims, lims, "k--", alpha=0.4, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Protein cosine distance")
    ax.set_ylabel("mRNA-encoder cosine distance")
    ax.set_title(f"Distance alignment ({title_suffix})\nPearson r={r:.4f}, R²={r**2:.4f}")
    ax.set_aspect("equal")

    plt.tight_layout()
    path = os.path.join(output_dir, f"distance_scatter_{title_suffix}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved distance scatter → {path}")


def plot_sibling_distances(
    embeddings: np.ndarray,
    cell_names: list,
    lineage_path: str,
    cell_type_path: str,
    output_dir: str,
):
    """Histograms of embedding distances: same-type vs different-type sibling pairs."""
    # Build cell type lookup
    cell_type_df = pd.read_csv(cell_type_path)
    name_to_type = {}
    for _, row in cell_type_df.iterrows():
        lineage = row["wormweb.lineage"]
        ct = row["wormweb.type"]
        if isinstance(ct, str):
            name_to_type[lineage] = ct

    # Build parent→children map from lineage tree
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    children_dict = defaultdict(list)

    def build_family(node, parent_name=None):
        name = map_names(node["did"])
        if parent_name is not None:
            children_dict[parent_name].append(name)
        for child in node.get("children", []):
            build_family(child, name)

    build_family(lineage_data)

    # Cell → embedding lookup
    cell_to_emb = dict(zip(cell_names, embeddings))

    # Find terminal sibling pairs
    terminal_set = set()
    for cn in cell_names:
        if cn in name_to_type:
            terminal_set.add(cn)

    same_pairs_dists = []
    diff_pairs_dists = []
    for parent, children in children_dict.items():
        term_children = [c for c in children if c in terminal_set and c in cell_to_emb]
        for i in range(len(term_children)):
            for j in range(i + 1, len(term_children)):
                a, b = term_children[i], term_children[j]
                type_a = name_to_type.get(a)
                type_b = name_to_type.get(b)
                if type_a is None or type_b is None:
                    continue
                dist = np.linalg.norm(cell_to_emb[a] - cell_to_emb[b])
                if type_a == type_b:
                    same_pairs_dists.append(dist)
                else:
                    diff_pairs_dists.append(dist)

    if len(same_pairs_dists) == 0 or len(diff_pairs_dists) == 0:
        print("  Not enough sibling pairs for distance analysis")
        return

    same_dists = np.array(same_pairs_dists)
    diff_dists = np.array(diff_pairs_dists)
    _, p_mw = mannwhitneyu(same_dists, diff_dists, alternative="less")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    bins = np.linspace(0, max(same_dists.max(), diff_dists.max()) * 1.05, 40)
    axes[0].hist(
        same_dists, bins=bins, alpha=0.65, color="steelblue",
        label=f"Same type (n={len(same_dists)})",
    )
    axes[0].hist(
        diff_dists, bins=bins, alpha=0.65, color="tomato",
        label=f"Diff type (n={len(diff_dists)})",
    )
    axes[0].axvline(np.median(same_dists), color="steelblue", linestyle="--")
    axes[0].axvline(np.median(diff_dists), color="tomato", linestyle="--")
    axes[0].set_xlabel("L2 distance in embedding space")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Per-cell mRNA embedding: sibling distances")
    axes[0].legend()
    axes[0].text(
        0.97, 0.97, f"MW p={p_mw:.2e}", transform=axes[0].transAxes,
        ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    violin_data = pd.DataFrame({
        "pair_type": ["Same type"] * len(same_dists) + ["Diff type"] * len(diff_dists),
        "distance": np.concatenate([same_dists, diff_dists]),
    })
    sns.violinplot(
        data=violin_data, x="pair_type", y="distance",
        hue="pair_type", palette=["steelblue", "tomato"],
        legend=False, ax=axes[1],
    )
    axes[1].set_ylabel("L2 distance in embedding space")
    axes[1].set_title("Sibling mRNA embedding distances")

    plt.tight_layout()
    path = os.path.join(output_dir, "sibling_distances.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved sibling distance plot → {path}")


def plot_lineage_distance_scatter(
    embeddings: np.ndarray,
    cell_names: list,
    lineage_path: str,
    output_dir: str,
):
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
        if not a1 or not a2:
            return -1
        common = set(a1) & set(a2)
        if not common:
            return len(a1) + len(a2)
        lca = max(common, key=lambda x: a1.index(x))
        return (len(a1) - a1.index(lca) - 1) + (len(a2) - a2.index(lca) - 1)

    cell_to_emb = dict(zip(cell_names, embeddings))
    cells_with_ancestry = [c for c in cell_names if c in node_ancestry]

    n_cells = min(500, len(cells_with_ancestry))
    rng = np.random.RandomState(42)
    sampled = rng.choice(cells_with_ancestry, size=n_cells, replace=False)

    lineage_dists = []
    emb_dists = []
    for i in range(len(sampled)):
        for j in range(i + 1, len(sampled)):
            ld = lineage_distance(sampled[i], sampled[j])
            if ld < 0:
                continue
            ed = np.linalg.norm(cell_to_emb[sampled[i]] - cell_to_emb[sampled[j]])
            lineage_dists.append(ld)
            emb_dists.append(ed)

    r_pearson, _ = pearsonr(lineage_dists, emb_dists)
    r_spearman, _ = spearmanr(lineage_dists, emb_dists)

    fig, ax = plt.subplots(figsize=(8, 6))
    df = pd.DataFrame({
        "lineage_distance": lineage_dists,
        "embedding_distance": emb_dists,
    })
    sns.boxplot(
        x="lineage_distance", y="embedding_distance", data=df,
        ax=ax, color="steelblue", fliersize=2,
    )
    ax.set_xlabel("Lineage Tree Distance")
    ax.set_ylabel("Embedding L2 Distance")
    ax.set_title(
        f"mRNA embedding distance vs lineage distance\n"
        f"Pearson r={r_pearson:.3f}, Spearman r={r_spearman:.3f}"
    )

    plt.tight_layout()
    path = os.path.join(output_dir, "lineage_distance_scatter.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved lineage distance scatter → {path}")


def plot_umap_embeddings(
    embeddings: np.ndarray,
    cell_names: list,
    lineage_path: str,
    cell_type_path: str,
    output_dir: str,
):
    """2D UMAP (or PCA fallback) of mRNA embeddings, colored by cell type and lineage depth."""
    # Cell type lookup
    cell_type_df = pd.read_csv(cell_type_path)
    name_to_type = {}
    for _, row in cell_type_df.iterrows():
        lineage = row["wormweb.lineage"]
        ct = row["wormweb.type"]
        if isinstance(ct, str):
            name_to_type[lineage] = ct

    # Lineage depth lookup
    with open(lineage_path, "r", encoding="utf-8") as f:
        lineage_data = json.load(f)

    node_depth = {}

    def dfs_depth(node, depth=0):
        name = map_names(node["did"])
        node_depth[name] = depth
        for child in node.get("children", []):
            dfs_depth(child, depth + 1)

    dfs_depth(lineage_data)

    # 2D reduction: UMAP if available, else PCA
    n_vis = min(2000, len(embeddings))
    rng = np.random.RandomState(42)
    vis_idx = rng.choice(len(embeddings), size=n_vis, replace=False)

    try:
        import umap
        reducer = umap.UMAP(n_components=2, random_state=42, verbose=False)
        emb_2d = reducer.fit_transform(embeddings[vis_idx])
        method = "UMAP"
    except ImportError:
        pca = PCA(n_components=2)
        emb_2d = pca.fit_transform(embeddings[vis_idx])
        method = "PCA"

    vis_cells = [cell_names[i] for i in vis_idx]
    cell_types = [name_to_type.get(c, "unknown") for c in vis_cells]
    depths = [node_depth.get(c, -1) for c in vis_cells]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: colored by cell type
    ax = axes[0]
    type_set = sorted(set(cell_types))
    cmap = plt.get_cmap("tab20", len(type_set))
    for i, ct in enumerate(type_set):
        mask = np.array([t == ct for t in cell_types])
        if mask.sum() > 0:
            ax.scatter(
                emb_2d[mask, 0], emb_2d[mask, 1],
                c=[cmap(i)], s=8, alpha=0.6, label=ct,
            )
    ax.set_xlabel(f"{method} 1")
    ax.set_ylabel(f"{method} 2")
    ax.set_title(f"{method} of mRNA embeddings — by cell type")
    ax.legend(fontsize=6, ncol=2, loc="upper right")

    # Right: colored by lineage depth (developmental time proxy)
    ax = axes[1]
    sc = ax.scatter(
        emb_2d[:, 0], emb_2d[:, 1], c=depths, cmap="viridis",
        s=6, alpha=0.5,
    )
    ax.set_xlabel(f"{method} 1")
    ax.set_ylabel(f"{method} 2")
    ax.set_title(f"{method} of mRNA embeddings — by lineage depth")
    plt.colorbar(sc, ax=ax, label="Lineage depth")

    plt.tight_layout()
    path = os.path.join(output_dir, "embedding_umap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved embedding {method} → {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Output saving
# ═══════════════════════════════════════════════════════════════════════════════


def save_outputs(
    model: torch.nn.Module,
    config,
    history: dict,
    embeddings: np.ndarray,
    all_cells: list,
    gene_names: list,
    baseline_metrics: dict,
    output_dir: str,
):
    """Save checkpoint, embeddings, and metrics."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Checkpoint
    ckpt_path = out / "model_checkpoint.pt"
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "config": {k: v for k, v in config.__dict__.items()},
            "gene_names": gene_names,
            "history": history,
        },
        str(ckpt_path),
    )
    print(f"  Saved checkpoint → {ckpt_path}")

    # Per-cell embeddings
    emb_cols = [f"emb_{i}" for i in range(embeddings.shape[1])]
    emb_df = pd.DataFrame(embeddings, index=all_cells, columns=emb_cols)
    emb_df.index.name = "cell_name"
    emb_path = out / "cell_embeddings_mrna.csv"
    emb_df.to_csv(str(emb_path))
    print(f"  Saved embeddings → {emb_path} ({emb_df.shape})")

    # Training curves
    plot_training_curves(history, str(out))

    # Metrics summary
    metrics = {
        "linear_probe_r2": baseline_metrics.get("linear_probe_r2", float("nan")),
        "best_val_pearson": float(max(history["val_pearson"])),
        "best_val_pearson_epoch": int(np.argmax(history["val_pearson"])),
        "final_val_pearson": float(history["val_pearson"][-1]),
        "best_val_total_loss": float(min(history["val_total"])),
        "epochs_run": len(history["train_total"]),
        "n_train_cells": baseline_metrics.get("n_train_cells", -1),
        "n_val_cells": baseline_metrics.get("n_val_cells", -1),
        "n_total_cells": baseline_metrics.get("n_total_cells", -1),
        "n_features": len(gene_names),
    }
    pd.Series(metrics).to_csv(str(out / "metrics.csv"), header=False)
    print(f"  Saved metrics → {out / 'metrics.csv'}")
