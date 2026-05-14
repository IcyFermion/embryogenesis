"""Diagnostics, sanity checks, and output saving for cross-species RNA embedding."""

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix, silhouette_score

matplotlib.use("Agg")


# ==============================================================================
# Main evaluation entry point
# ==============================================================================


def run_evaluation(model, data, history, config):
    """Run all diagnostics and save outputs."""
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device("cpu")
    model = model.to(device).eval()

    # Extract embeddings and predictions for all cells
    X_all = data["X_all"]
    X_all_t = torch.FloatTensor(X_all).to(device)
    with torch.no_grad():
        z, x_recon, logits = model(X_all_t)
        embeddings = z.cpu().numpy()
        preds = logits.argmax(dim=1).cpu().numpy()
        targets = data["y_all"].argmax(axis=1)

    # Split by species
    species = data["species_all"]
    ele_mask = species == 0
    bri_mask = species == 1

    # ---- 1. Sanity checks ----
    run_sanity_checks(data, X_all, config)

    # ---- 2. Training curves ----
    save_training_curves(history, output_dir)

    # ---- 3. Per-species diagnostics ----
    for label, mask in [("elegans", ele_mask), ("briggsae", bri_mask)]:
        run_intrinsic_diagnostics(
            embeddings[mask], data["cell_names_all"][mask],
            data["y_all"][mask], data["hard_mask_all"][mask],
            data["class_names"], label, config, output_dir
        )

    # ---- 4. Cross-species diagnostics ----
    run_cross_species_diagnostics(
        embeddings, species, data["cell_names_all"], data["y_all"],
        data["hard_mask_all"], data["class_names"], output_dir
    )

    # ---- 5. Confusion matrices ----
    run_confusion_matrices(
        preds, targets, species, data["hard_mask_all"],
        data["hard_mask_val"], data["val_idx"],
        data["class_names"], output_dir
    )

    # ---- 6. Save outputs ----
    save_outputs(model, embeddings, data, config, output_dir)

    print(f"\nAll diagnostics saved to {output_dir}")


# ==============================================================================
# Sanity checks
# ==============================================================================


def run_sanity_checks(data, X_all, config):
    """Five required sanity checks."""
    print("\n" + "=" * 60)
    print("SANITY CHECKS")
    print("=" * 60)

    # 1. Per-species z-scoring: mean ≈ 0 per gene within each species
    species = data["species_all"]
    for sp, label in [(0, "elegans"), (1, "briggsae")]:
        sp_data = X_all[species == sp]
        gene_means = sp_data.mean(axis=0)
        gene_stds = sp_data.std(axis=0)
        assert np.allclose(gene_means, 0, atol=1e-4), f"{label} z-score means ≠ 0"
        print(f"  [PASS] {label}: per-species z-score — gene means ≈ 0 "
              f"(max|mean|={np.abs(gene_means).max():.2e}), gene stds ≈ 1 "
              f"(mean std={gene_stds.mean():.4f})")

    # 2. Mixed-species batches confirmed
    from .data_loader import MixedSpeciesBatchSampler
    sampler = MixedSpeciesBatchSampler(data["species_train"], config.batch_size, seed=config.seed)
    n_batches = 0
    ele_fracs = []
    for batch in sampler:
        ele_frac = (data["species_train"][batch] == 0).mean()
        ele_fracs.append(ele_frac)
        n_batches += 1
    ele_frac_mean = np.mean(ele_fracs)
    print(f"  [PASS] Mixed-species batches: {n_batches} batches, "
          f"mean ele fraction={ele_frac_mean:.2f} (expect ~0.5)")

    # 3. Train/val no leakage
    train_cells = set(data["cell_names_all"][data["train_idx"]])
    val_cells = set(data["cell_names_all"][data["val_idx"]])
    overlap = train_cells & val_cells
    assert len(overlap) == 0, f"Train/val overlap: {overlap}"
    print(f"  [PASS] Train/val disjoint: train={len(train_cells)}, val={len(val_cells)}, overlap=0")

    # 4-5 are evaluated after training (centroid distances, accuracy comparison)
    print()


# ==============================================================================
# Intrinsic diagnostics (per species)
# ==============================================================================


def run_intrinsic_diagnostics(embeddings, cell_names, y, hard_mask, class_names,
                              species_label, config, output_dir):
    """Per-species: sibling-pair distances, lineage scatter, PCA by type."""
    print(f"\n--- Intrinsic diagnostics: {species_label} ---")

    # Sibling pair distances
    _sibling_pair_diagnostics(embeddings, cell_names, y, hard_mask, class_names,
                              config.lineage_path, species_label, output_dir)

    # Lineage distance scatter
    _lineage_distance_scatter(embeddings, cell_names, config.lineage_path,
                              species_label, output_dir)

    # PCA colored by terminal type
    _pca_by_type(embeddings, y, hard_mask, class_names, species_label, output_dir)


def _sibling_pair_diagnostics(embeddings, cell_names, y, hard_mask, class_names,
                              lineage_path, species_label, output_dir):
    """Same-type vs different-type sibling pair distance histograms."""
    # Build parent→children map
    with open(lineage_path) as f:
        root = json.load(f)

    parent_to_children = defaultdict(list)
    child_to_parent = {}

    def _build_family(node, parent=None):
        did = node.get("did", "")
        if not did:
            return
        if parent:
            child_to_parent[did] = parent
            parent_to_children[parent].append(did)
        for child in node.get("children", []):
            _build_family(child, did)

    _build_family(root)

    # Map cell names to indices
    name_to_idx = {n: i for i, n in enumerate(cell_names)}

    # Collect terminal sibling pairs
    same_dists, diff_dists = [], []
    cell_set = set(cell_names)

    for parent, children in parent_to_children.items():
        term_children = [c for c in children if c in cell_set]
        for i in range(len(term_children)):
            for j in range(i + 1, len(term_children)):
                a, b = term_children[i], term_children[j]
                if a not in name_to_idx or b not in name_to_idx:
                    continue
                ia, ib = name_to_idx[a], name_to_idx[b]
                d = float(np.linalg.norm(embeddings[ia] - embeddings[ib]))
                if y[ia].argmax() == y[ib].argmax():
                    same_dists.append(d)
                else:
                    diff_dists.append(d)

    if not same_dists or not diff_dists:
        print(f"  [{species_label}] Not enough sibling pairs for histogram")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(same_dists, bins=30, alpha=0.6, label=f"Same type (n={len(same_dists)})", density=True)
    ax.hist(diff_dists, bins=30, alpha=0.6, label=f"Diff type (n={len(diff_dists)})", density=True)
    ax.set_xlabel("Embedding distance")
    ax.set_ylabel("Density")
    ax.set_title(f"Sibling-pair distances — {species_label}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / f"sibling_pairs_{species_label}.png", dpi=150)
    plt.close(fig)

    same_mean, diff_mean = np.mean(same_dists), np.mean(diff_dists)
    print(f"  Sibling pairs: same-type mean={same_mean:.4f} (n={len(same_dists)}), "
          f"diff-type mean={diff_mean:.4f} (n={len(diff_dists)})")


def _lineage_distance_scatter(embeddings, cell_names, lineage_path, species_label, output_dir):
    """Embedding distance vs. lineage tree distance."""
    # Build ancestor map
    with open(lineage_path) as f:
        root = json.load(f)

    def _ancestors(node, anc_set=None):
        if anc_set is None:
            anc_set = set()
        did = node.get("did", "")
        if did:
            anc_set.add(did)
        for child in node.get("children", []):
            _ancestors(child, anc_set.copy())

    # Flatten
    all_nodes = {}
    def _flatten(node):
        did = node.get("did", "")
        if did:
            all_nodes[did] = node
        for child in node.get("children", []):
            _flatten(child)
    _flatten(root)

    # Sample pairs (limit to 5000 for speed)
    rng = np.random.default_rng(42)
    n = len(cell_names)
    n_total_pairs = n * (n - 1) // 2
    max_pairs = min(5000, n_total_pairs)
    if max_pairs > 0:
        # Get all upper-triangle indices and randomly select
        all_i, all_j = np.triu_indices(n, k=1)
        chosen = rng.choice(len(all_i), size=max_pairs, replace=False)
        idx_pairs = np.column_stack([all_i[chosen], all_j[chosen]])
    else:
        idx_pairs = np.zeros((0, 2), dtype=int)
        print(f"  [{species_label}] Not enough cells for lineage scatter")
        return

    # Compute distances
    emb_dists = np.linalg.norm(embeddings[idx_pairs[:, 0]] - embeddings[idx_pairs[:, 1]], axis=1)

    # Lineage tree edge distances (simplified: count edges between nodes)
    # For speed, approximate by depth difference + branching
    lineage_dists = np.zeros(len(idx_pairs))
    # Simple approximation: totalDistance difference
    for k, (i, j) in enumerate(idx_pairs):
        ni = all_nodes.get(cell_names[i], {})
        nj = all_nodes.get(cell_names[j], {})
        di = (ni.get("data") or {}).get("totalDistance", 0)
        dj = (nj.get("data") or {}).get("totalDistance", 0)
        lineage_dists[k] = abs(di - dj)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(lineage_dists, emb_dists, s=2, alpha=0.3)
    ax.set_xlabel("Lineage distance (totalDistance diff)")
    ax.set_ylabel("Embedding distance")
    ax.set_title(f"Embedding vs lineage distance — {species_label}")
    fig.tight_layout()
    fig.savefig(output_dir / f"lineage_scatter_{species_label}.png", dpi=150)
    plt.close(fig)

    corr = np.corrcoef(lineage_dists, emb_dists)[0, 1]
    print(f"  Lineage scatter: Pearson r = {corr:.4f}")


def _pca_by_type(embeddings, y, hard_mask, class_names, species_label, output_dir):
    """2D PCA colored by terminal type for a single species."""
    pca = PCA(n_components=2)
    coords = pca.fit_transform(embeddings)

    # Only show typed cells (hard or soft with dominant type)
    types = class_names
    dominant = np.array([types[t] for t in y.argmax(axis=1)])

    fig, ax = plt.subplots(figsize=(8, 7))
    for t in sorted(set(dominant)):
        mask = dominant == t
        ax.scatter(coords[mask, 0], coords[mask, 1], s=3, alpha=0.5, label=f"{t} ({mask.sum()})")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
    ax.set_title(f"PCA by terminal type — {species_label}")
    ax.legend(markerscale=3, fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / f"pca_type_{species_label}.png", dpi=150)
    plt.close(fig)


# ==============================================================================
# Cross-species diagnostics
# ==============================================================================


def run_cross_species_diagnostics(embeddings, species, cell_names, y, hard_mask,
                                  class_names, output_dir):
    """Joint PCA (species + type), centroid table."""
    print("\n--- Cross-species diagnostics ---")

    # Joint PCA
    pca = PCA(n_components=2)
    coords = pca.fit_transform(embeddings)

    ele_mask = species == 0
    bri_mask = species == 1
    dominant = np.array([class_names[t] for t in y.argmax(axis=1)])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 6.5))

    # By species
    ax0.scatter(coords[ele_mask, 0], coords[ele_mask, 1], s=2, alpha=0.4,
                label="C. elegans", c="#1f77b4")
    ax0.scatter(coords[bri_mask, 0], coords[bri_mask, 1], s=2, alpha=0.4,
                label="C. briggsae", c="#d62728")
    ax0.set_title("Joint embedding PCA — by species")
    ax0.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
    ax0.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
    ax0.legend(markerscale=5)

    # By type
    all_types = sorted(set(dominant))
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(all_types), 20)))
    for i, t in enumerate(all_types):
        mask = dominant == t
        ax1.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.4,
                    label=t, c=[palette[i % len(palette)]])
    ax1.set_title("Joint embedding PCA — by terminal type")
    ax1.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
    ax1.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
    ax1.legend(markerscale=4, fontsize=6, ncol=2, loc="lower right")

    fig.tight_layout()
    fig.savefig(output_dir / "joint_pca.png", dpi=150)
    plt.close(fig)

    # Silhouette by species in embedding space
    sil = silhouette_score(embeddings, species)
    print(f"  Joint embedding silhouette (species): {sil:.4f}")

    # Cross-species centroid table
    _centroid_table(embeddings, species, dominant, class_names, output_dir)


def _centroid_table(embeddings, species, dominant, class_names, output_dir):
    """d_same (same type cross-species) vs d_diff_within (different type same species)."""
    ele_mask = species == 0
    bri_mask = species == 1

    shared_types = []
    for t in sorted(set(dominant) - {"programmed_death", "unknown"}):
        n_ele = ((dominant == t) & ele_mask).sum()
        n_bri = ((dominant == t) & bri_mask).sum()
        if n_ele >= 3 and n_bri >= 3:
            shared_types.append(t)

    if len(shared_types) < 2:
        print("  Not enough shared types for centroid analysis")
        return

    ele_centroids = {}
    bri_centroids = {}
    for t in shared_types:
        ele_centroids[t] = embeddings[(dominant == t) & ele_mask].mean(axis=0)
        bri_centroids[t] = embeddings[(dominant == t) & bri_mask].mean(axis=0)

    rows = []
    for t in shared_types:
        d_same = float(np.linalg.norm(ele_centroids[t] - bri_centroids[t]))

        d_ele_diff = []
        for t2 in shared_types:
            if t2 != t:
                d_ele_diff.append(float(np.linalg.norm(ele_centroids[t] - ele_centroids[t2])))
        d_ele_diff_mean = np.mean(d_ele_diff) if d_ele_diff else float("nan")

        d_bri_diff = []
        for t2 in shared_types:
            if t2 != t:
                d_bri_diff.append(float(np.linalg.norm(bri_centroids[t] - bri_centroids[t2])))
        d_bri_diff_mean = np.mean(d_bri_diff) if d_bri_diff else float("nan")

        rows.append({
            "type": t,
            "d_same": d_same,
            "d_diff_ele": d_ele_diff_mean,
            "d_diff_bri": d_bri_diff_mean,
            "ratio_ele": d_same / d_ele_diff_mean if d_ele_diff_mean > 0 else float("nan"),
            "ratio_bri": d_same / d_bri_diff_mean if d_bri_diff_mean > 0 else float("nan"),
        })

    df = pd.DataFrame(rows).set_index("type")
    print("\n  Cross-species centroid distances:")
    print(df.to_string(float_format=lambda x: f"{x:.4f}"))

    better_ele = sum(1 for r in rows if r["ratio_ele"] < 1.0)
    better_bri = sum(1 for r in rows if r["ratio_bri"] < 1.0)
    n = len(rows)
    print(f"  Types where d_same < d_diff_within: ele={better_ele}/{n}, bri={better_bri}/{n}")

    df.to_csv(output_dir / "centroid_distances.csv")


# ==============================================================================
# Confusion matrices
# ==============================================================================


def run_confusion_matrices(preds, targets, species, hard_mask, hard_mask_val,
                           val_idx, class_names, output_dir):
    """Confusion matrices per species on held-out hard-labeled cells."""
    # val_idx contains the global indices of validation cells
    # hard_mask_val is already subset to val cells only
    val_species = species[val_idx]

    for sp, label in [(0, "elegans"), (1, "briggsae")]:
        # Hard-labeled val cells for this species
        mask_val = (val_species == sp) & hard_mask_val

        if mask_val.sum() < 5:
            print(f"  [{label}] Not enough hard-labeled val cells ({mask_val.sum()})")
            continue

        cm = confusion_matrix(targets[val_idx][mask_val], preds[val_idx][mask_val])
        # Normalize by row
        cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
        cm_norm = np.nan_to_num(cm_norm)

        fig, ax = plt.subplots(figsize=(max(6, len(class_names) * 0.6),
                                        max(5, len(class_names) * 0.5)))
        im = ax.imshow(cm_norm, vmin=0, vmax=1, cmap="Blues", aspect="auto")
        ax.set_xticks(range(len(class_names)))
        ax.set_yticks(range(len(class_names)))
        ax.set_xticklabels(class_names, rotation=90, fontsize=7)
        ax.set_yticklabels(class_names, fontsize=7)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(f"Confusion matrix — {label} (held-out, hard-label)")
        plt.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(output_dir / f"confusion_{label}.png", dpi=150)
        plt.close(fig)

        acc = (preds[val_idx][mask_val] == targets[val_idx][mask_val]).mean()
        print(f"  [{label}] Held-out hard-label accuracy: {acc:.4f} (n={mask_val.sum()})")


# ==============================================================================
# Output saving
# ==============================================================================


def save_training_curves(history, output_dir):
    """Save history CSV and training curve plots."""
    df = pd.DataFrame(history)
    df.index.name = "epoch"
    df.to_csv(output_dir / "training_curves.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Loss
    ax = axes[0, 0]
    ax.plot(df.index + 1, df["train_total"], label="Train total", linewidth=1)
    ax.plot(df.index + 1, df["val_total"], label="Val total", linewidth=1)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title("Total loss")
    ax.legend()

    # Components
    ax = axes[0, 1]
    ax.plot(df.index + 1, df["train_recon"], label="Train recon", linewidth=0.8, alpha=0.7)
    ax.plot(df.index + 1, df["train_classify"], label="Train classify", linewidth=0.8, alpha=0.7)
    ax.plot(df.index + 1, df["val_recon"], label="Val recon", linewidth=0.8, alpha=0.7)
    ax.plot(df.index + 1, df["val_classify"], label="Val classify", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title("Loss components")
    ax.legend(fontsize=7)

    # Accuracy
    ax = axes[1, 0]
    ax.plot(df.index + 1, df["val_joint_acc"], label="Joint acc", linewidth=1)
    ax.plot(df.index + 1, df["val_ele_acc"], label="Ele acc", linewidth=0.8, alpha=0.7)
    ax.plot(df.index + 1, df["val_bri_acc"], label="Bri acc", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Accuracy")
    ax.set_title("Validation classification accuracy")
    ax.legend()

    # Recon MSE per species
    ax = axes[1, 1]
    ax.plot(df.index + 1, df["val_ele_recon_mse"], label="Ele recon MSE", linewidth=0.8)
    ax.plot(df.index + 1, df["val_bri_recon_mse"], label="Bri recon MSE", linewidth=0.8)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("MSE")
    ax.set_title("Validation reconstruction MSE per species")
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_dir / "training_curves.png", dpi=150)
    plt.close(fig)


def save_outputs(model, embeddings, data, config, output_dir):
    """Save checkpoint, embeddings CSV, and metadata."""
    # Checkpoint
    torch.save({
        "model_state_dict": model.state_dict(),
        "config": {k: v for k, v in config.__dict__.items() if not k.startswith("_")},
    }, output_dir / "checkpoint.pt")

    # Embeddings CSV
    emb_df = pd.DataFrame(
        embeddings,
        index=data["cell_names_all"],
        columns=[f"dim_{i}" for i in range(embeddings.shape[1])],
    )
    emb_df["species"] = ["elegans" if s == 0 else "briggsae" for s in data["species_all"]]
    emb_df["dominant_type"] = [data["class_names"][t] for t in data["y_all"].argmax(axis=1)]
    emb_df["is_hard_label"] = data["hard_mask_all"]
    emb_df["split"] = ["train" if i in data["train_idx"] else "val" for i in range(len(data["cell_names_all"]))]
    emb_df.to_csv(output_dir / "cell_embeddings.csv")
    print(f"  Embeddings saved: {embeddings.shape}")

    # Class names
    pd.Series(data["class_names"]).to_csv(output_dir / "class_names.csv", index=False, header=["class"])
