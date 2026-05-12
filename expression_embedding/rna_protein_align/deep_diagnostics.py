"""Deep validation diagnostics for mRNA-to-protein alignment model.

Checks:
  1. Lineage-distance audit: distribution of lineage-tree distance from each
     val cell to its nearest train cell.
  2. Permutation test: shuffle mRNA→cell mapping for val cells, re-evaluate.
  3. Strict-branch split: hold out an entire major lineage branch, retrain.
  4. Pair-set check: verify no cell ID appears in both train and val pairs.
  5. Within-val vs cross-split evaluation: Pearson on (val,val) pairs vs
     (train,val) pairs.

Usage:
    python expression_embedding/rna_protein_align/deep_diagnostics.py \
        --output expression_embedding/results/rna_protein_align_final \
        [--skip-strict] [--strict-branch AB]
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

# Path setup
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
from expression_embedding.rna_protein_align.model import RnaEncoderModel
from expression_embedding.rna_protein_align.trainer import (
    _resolve_device,
    make_fixed_val_pairs,
    pearson_on_pairs,
    train,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Shared utilities
# ═══════════════════════════════════════════════════════════════════════════════


def load_model_and_data(output_dir: str):
    """Load the trained model checkpoint and data from an output directory."""
    out = Path(output_dir)
    ckpt = torch.load(str(out / "model_checkpoint.pt"), map_location="cpu",
                      weights_only=False)
    saved_cfg = ckpt["config"]

    # Reconstruct config
    config = Config()
    for k, v in saved_cfg.items():
        if hasattr(config, k):
            setattr(config, k, v)

    # Re-load data with the same seed
    data = load_rna_protein_data(
        rna_path=config.rna_path,
        protein_path=config.protein_emb_path,
        lineage_path=config.lineage_path,
        sublineage_depth=config.sublineage_depth,
        val_fraction=config.val_fraction,
        seed=config.seed,
        log_transform=config.log_transform,
    )

    # Rebuild model and load weights
    device = _resolve_device(config.device)
    model = RnaEncoderModel(
        n_features=data["n_features"],
        hidden_dims=config.hidden_dims,
        embed_dim=config.embed_dim,
        dropout=config.dropout,
        use_layer_norm=config.use_layer_norm,
    ).to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    return model, config, data, device, out


def build_lineage_ancestry(lineage_path: str) -> dict:
    """Return dict mapping cell_name → list of ancestor names (from root)."""
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


def lineage_distance(cell1: str, cell2: str, node_ancestry: dict) -> int:
    """Number of edges between two cells in the lineage tree."""
    a1 = node_ancestry.get(cell1)
    a2 = node_ancestry.get(cell2)
    if a1 is None or a2 is None:
        return -1
    common = set(a1) & set(a2)
    if not common:
        return len(a1) + len(a2)
    lca = max(common, key=lambda x: a1.index(x))  # deepest common ancestor
    return (len(a1) - a1.index(lca) - 1) + (len(a2) - a2.index(lca) - 1)


# ═══════════════════════════════════════════════════════════════════════════════
# Check 1: Lineage-distance audit
# ═══════════════════════════════════════════════════════════════════════════════


def check_lineage_audit(train_cells, val_cells, lineage_path, output_dir):
    """For each val cell, compute lineage distance to nearest train cell."""
    print("\n" + "=" * 60)
    print("Check 1: Lineage-distance audit")
    print("=" * 60)

    ancestry = build_lineage_ancestry(lineage_path)
    dists = []
    unmatched = 0
    for vc in val_cells:
        min_dist = float("inf")
        for tc in train_cells:
            d = lineage_distance(vc, tc, ancestry)
            if d >= 0 and d < min_dist:
                min_dist = d
        if min_dist == float("inf"):
            unmatched += 1
        else:
            dists.append(min_dist)

    dists = np.array(dists)
    print(f"  Val cells: {len(val_cells)}")
    print(f"  Unmatched (no common ancestor found): {unmatched}")
    print(f"  Nearest-train lineage distance:")
    print(f"    min={dists.min():.0f}  median={np.median(dists):.0f}  mean={dists.mean():.1f}  max={dists.max():.0f}")

    # Binned stats
    bins = [0, 2, 4, 6, 8, 10, 15, 20, 30, 50]
    print(f"  Distance distribution:")
    for i in range(len(bins) - 1):
        count = ((dists >= bins[i]) & (dists < bins[i + 1])).sum()
        pct = count / len(dists) * 100
        bar = "#" * int(pct / 2)
        print(f"    [{bins[i]:2d}, {bins[i+1]:2d}): {count:4d} ({pct:5.1f}%) {bar}")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(dists, bins=40, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(np.median(dists), color="tomato", linestyle="--",
               label=f"median = {np.median(dists):.0f}")
    ax.set_xlabel("Lineage-tree distance to nearest train cell")
    ax.set_ylabel("Count")
    ax.set_title("Check 1: Lineage-distance audit of val cells")
    ax.legend()
    plt.tight_layout()
    path = os.path.join(output_dir, "diagnostic_lineage_audit.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")

    # Verdict
    median = np.median(dists)
    if median < 4:
        print(f"  [WARN] Median distance {median:.0f} is small — val cells are close relatives of train")
    else:
        print(f"  [OK] Median distance {median:.0f} — val cells are reasonably separated")
    return dists


# ═══════════════════════════════════════════════════════════════════════════════
# Check 2: Permutation test
# ═══════════════════════════════════════════════════════════════════════════════


def check_permutation(model, X_val, prot_val, device, output_dir, n_perm=20):
    """Shuffle cell→mRNA mapping for val cells, re-evaluate Pearson."""
    print("\n" + "=" * 60)
    print("Check 2: Permutation test")
    print("=" * 60)

    # True correlation (from fixed pairs)
    n_pairs = 5000
    fixed_a, fixed_b = make_fixed_val_pairs(len(X_val), n_pairs, seed=42)
    true_r, _, _ = pearson_on_pairs(model, X_val, prot_val, fixed_a, fixed_b, device)
    print(f"  True Pearson (val-val pairs): {true_r:.4f}")

    # Permuted correlations
    perm_rs = []
    rng = np.random.RandomState(42)
    for i in range(n_perm):
        perm_idx = rng.permutation(len(X_val))
        X_perm = X_val[perm_idx]
        r, _, _ = pearson_on_pairs(model, X_perm, prot_val, fixed_a, fixed_b, device)
        perm_rs.append(r)

    perm_rs = np.array(perm_rs)
    print(f"  Permuted Pearson: mean={perm_rs.mean():.4f}  std={perm_rs.std():.4f}  "
          f"min={perm_rs.min():.4f}  max={perm_rs.max():.4f}")
    print(f"  Drop: {true_r - perm_rs.mean():.4f} (should be large — true >> permuted)")

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(perm_rs, bins=15, color="lightcoral", edgecolor="white", alpha=0.8,
            label=f"Permuted (n={n_perm})")
    ax.axvline(true_r, color="steelblue", linewidth=2.5,
               label=f"True (r={true_r:.4f})")
    ax.axvline(0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("Count")
    ax.set_title("Check 2: Permutation test — shuffled mRNA→cell mapping")
    ax.legend()
    plt.tight_layout()
    path = os.path.join(output_dir, "diagnostic_permutation.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")

    if perm_rs.mean() > 0.1:
        print(f"  [WARN] Permuted mean {perm_rs.mean():.4f} > 0.1 — possible leakage")
    else:
        print(f"  [OK] Permuted near zero — no cell-identity leakage")

    return true_r, perm_rs


# ═══════════════════════════════════════════════════════════════════════════════
# Check 3: Strict-branch split
# ═══════════════════════════════════════════════════════════════════════════════


def check_strict_branch(config, lineage_path, output_dir, branch="AB"):
    """Hold out an entire major lineage branch, retrain, report Pearson."""
    print("\n" + "=" * 60)
    print(f"Check 3: Strict-branch split (hold out '{branch}' descendants)")
    print("=" * 60)

    # Load the full data
    data = load_rna_protein_data(
        rna_path=config.rna_path,
        protein_path=config.protein_emb_path,
        lineage_path=config.lineage_path,
        sublineage_depth=config.sublineage_depth,
        val_fraction=0.0,  # we'll do our own split
        seed=config.seed,
        log_transform=config.log_transform,
    )

    # Identify branch descendants
    ancestry = build_lineage_ancestry(lineage_path)
    all_cells = data["all_cells"]

    if branch == "AB":
        # All cells under AB (names starting with "AB")
        val_cells = [c for c in all_cells if c.startswith("AB")]
    elif branch == "P1":
        # All cells under P1 (not under AB)
        # P1 descendants are everything except AB descendants
        # P1 → P2 → P3 → P4 → (MS, E, C, D, Z2, Z3)
        ab_cells = set(c for c in all_cells if c.startswith("AB"))
        val_cells = [c for c in all_cells if c not in ab_cells]
    else:
        raise ValueError(f"Unknown branch: {branch}. Use 'AB' or 'P1'.")

    val_set = set(val_cells)
    train_cells = [c for c in all_cells if c not in val_set]

    print(f"  Branch '{branch}': {len(val_cells)} val cells, {len(train_cells)} train cells")

    if len(val_cells) < 10:
        print("  [SKIP] Too few val cells")
        return None

    # Map cell names to indices in the all_cells array
    cell_to_idx = {c: i for i, c in enumerate(all_cells)}
    train_idx = np.array([cell_to_idx[c] for c in train_cells], dtype=int)
    val_idx = np.array([cell_to_idx[c] for c in val_cells], dtype=int)

    X_train = data["X_all"][train_idx]
    X_val = data["X_all"][val_idx]
    prot_train = data["prot_all"][train_idx]
    prot_val = data["prot_all"][val_idx]

    prot_dist = precompute_protein_distances(data["prot_all"])
    train_samp = PairSampler(train_idx, prot_dist, config.near_fraction, config.seed)
    val_samp = PairSampler(val_idx, prot_dist, 0.0, config.seed + 1)

    # Retrain with the same config
    print(f"  Retraining with strict split...")
    model, history = train(
        config, X_train, X_val, prot_train, prot_val,
        prot_dist, train_samp, val_samp, data["n_features"],
    )

    best = max(history["val_pearson"])
    best_ep = history["val_pearson"].index(best)
    print(f"  Strict-split best Pearson: {best:.4f} at epoch {best_ep}")

    # Save strict-split curves
    strict_dir = Path(output_dir) / f"strict_split_{branch}"
    strict_dir.mkdir(exist_ok=True)
    pd.DataFrame(history).to_csv(strict_dir / "training_curves.csv", index=False)

    # Save model
    torch.save(
        {"model_state_dict": model.state_dict(), "config": {k: v for k, v in config.__dict__.items()}},
        str(strict_dir / "model_checkpoint.pt"),
    )

    # Plot comparison
    # Load original history
    orig_hist = pd.read_csv(Path(output_dir) / "training_curves.csv")
    orig_best = orig_hist["val_pearson"].max()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(history["val_pearson"], color="tomato", alpha=0.8,
            label=f"Strict '{branch}' split (best={best:.3f})")
    ax.plot(orig_hist["val_pearson"], color="steelblue", alpha=0.8,
            label=f"Original sublineage split (best={orig_best:.3f})")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Val Pearson r")
    ax.set_title(f"Check 3: Strict-branch split comparison")
    ax.legend()
    plt.tight_layout()
    path = str(strict_dir / "split_comparison.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")

    if best < 0.3:
        print(f"  [WARN] Strict-split Pearson {best:.4f} < 0.3 — original split may have been easy")
    else:
        print(f"  [OK] Strict-split Pearson {best:.4f} — model generalizes across branches")

    return {"best_pearson": best, "best_epoch": best_ep, "branch": branch,
            "n_train": len(train_cells), "n_val": len(val_cells)}


# ═══════════════════════════════════════════════════════════════════════════════
# Check 4: Pair-set check
# ═══════════════════════════════════════════════════════════════════════════════


def check_pair_sets(config, data, prot_dist_matrix, output_dir):
    """Verify no cell ID appears in both train and val sampled pairs."""
    print("\n" + "=" * 60)
    print("Check 4: Pair-set leakage check")
    print("=" * 60)

    train_samp = PairSampler(data["train_indices"], prot_dist_matrix,
                              config.near_fraction, config.seed)
    val_samp = PairSampler(data["val_indices"], prot_dist_matrix,
                            0.0, config.seed + 1)

    all_cells = data["all_cells"]
    train_cells = data["train_cells"]
    val_cells = data["val_cells"]

    # Sample pairs
    tr_a, tr_b = train_samp.sample(100)
    va_a, va_b = val_samp.sample(100)

    # Map local indices to cell names
    train_pairs = [(train_cells[a], train_cells[b]) for a, b in zip(tr_a, tr_b)]
    val_pairs = [(val_cells[a], val_cells[b]) for a, b in zip(va_a, va_b)]

    train_cells_in_pairs = set()
    for a, b in train_pairs:
        train_cells_in_pairs.add(a)
        train_cells_in_pairs.add(b)

    val_cells_in_pairs = set()
    for a, b in val_pairs:
        val_cells_in_pairs.add(a)
        val_cells_in_pairs.add(b)

    overlap = train_cells_in_pairs & val_cells_in_pairs

    print(f"  Unique cells in 100 train pairs: {len(train_cells_in_pairs)}")
    print(f"  Unique cells in 100 val pairs:   {len(val_cells_in_pairs)}")
    print(f"  Overlap: {len(overlap)}")

    print(f"\n  First 10 train pairs:")
    for i, (a, b) in enumerate(train_pairs[:10]):
        print(f"    {a} ↔ {b}")

    print(f"\n  First 10 val pairs:")
    for i, (a, b) in enumerate(val_pairs[:10]):
        print(f"    {a} ↔ {b}")

    if len(overlap) > 0:
        print(f"\n  [FAIL] {len(overlap)} cells appear in both train and val pairs!")
        print(f"  Overlapping cells: {sorted(overlap)[:20]}")
    else:
        print(f"\n  [PASS] No cell appears in both train and val pairs")

    return len(overlap)


# ═══════════════════════════════════════════════════════════════════════════════
# Check 5: Within-val vs cross-split evaluation
# ═══════════════════════════════════════════════════════════════════════════════


def check_cross_split_eval(model, config, data, device, output_dir):
    """Compare Pearson on (val,val) pairs vs (train,val) pairs."""
    print("\n" + "=" * 60)
    print("Check 5: Within-val vs cross-split evaluation")
    print("=" * 60)

    X_train = data["X_train"]
    X_val = data["X_val"]
    prot_train = data["prot_train"]
    prot_val = data["prot_val"]
    train_cells = data["train_cells"]
    val_cells = data["val_cells"]

    n_pairs = 5000
    n_val = len(X_val)
    n_train = len(X_train)

    # Case (a): both cells in val — honest test
    va_a, va_b = make_fixed_val_pairs(n_val, n_pairs, seed=42)
    r_val_val, mrna_vv, prot_vv = pearson_on_pairs(
        model, X_val, prot_val, va_a, va_b, device,
    )

    # Case (b): one train, one val — partial leakage
    rng = np.random.RandomState(42)
    tv_a = rng.randint(0, n_train, size=n_pairs)
    tv_b = rng.randint(0, n_val, size=n_pairs)

    # Encode both sets
    model.eval()
    X_tr_t = torch.as_tensor(X_train, dtype=torch.float32, device=device)
    X_va_t = torch.as_tensor(X_val, dtype=torch.float32, device=device)
    z_train = model.encoder(X_tr_t)
    z_val = model.encoder(X_va_t)

    # mRNA cosine distances (train,val pairs)
    cos_sim_tv = (z_train[tv_a] * z_val[tv_b]).sum(dim=1)
    mrna_tv = (1.0 - cos_sim_tv).detach().cpu().numpy()

    # Protein cosine distances
    prot_a = prot_train[tv_a]
    prot_b = prot_val[tv_b]
    prot_a_n = prot_a / np.clip(np.linalg.norm(prot_a, axis=1, keepdims=True), 1e-10, None)
    prot_b_n = prot_b / np.clip(np.linalg.norm(prot_b, axis=1, keepdims=True), 1e-10, None)
    prot_tv = np.clip(1.0 - (prot_a_n * prot_b_n).sum(axis=1), 0.0, 2.0)

    from scipy.stats import pearsonr
    r_train_val, _ = pearsonr(mrna_tv, prot_tv)

    print(f"  (a) val-val pairs:   r = {r_val_val:.4f}  ← honest metric")
    print(f"  (b) train-val pairs:  r = {r_train_val:.4f}  ← partial leakage")
    print(f"  Difference:           Δ = {r_train_val - r_val_val:.4f}")

    # If cross-split correlation is much higher, the model may be exploiting
    # cell-lineage proximity rather than mRNA→protein mapping
    if r_train_val - r_val_val > 0.15:
        print(f"  [WARN] Cross-split r is substantially higher — "
              f"possible lineage proximity leakage")
    else:
        print(f"  [OK] val-val and train-val are comparable")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.scatter(prot_vv, mrna_vv, s=3, alpha=0.3, color="steelblue",
               edgecolors="none")
    ax.set_xlabel("Protein cosine distance")
    ax.set_ylabel("mRNA-encoder cosine distance")
    ax.set_title(f"val-val pairs (r={r_val_val:.4f})")
    ax.plot([0, 2], [0, 2], "k--", alpha=0.3, linewidth=0.8)

    ax = axes[1]
    ax.scatter(prot_tv, mrna_tv, s=3, alpha=0.3, color="tomato",
               edgecolors="none")
    ax.set_xlabel("Protein cosine distance")
    ax.set_ylabel("mRNA-encoder cosine distance")
    ax.set_title(f"train-val pairs (r={r_train_val:.4f})")
    ax.plot([0, 2], [0, 2], "k--", alpha=0.3, linewidth=0.8)

    plt.tight_layout()
    path = os.path.join(output_dir, "diagnostic_cross_split.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")

    return r_val_val, r_train_val


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Deep validation diagnostics for mRNA-protein alignment"
    )
    parser.add_argument("--output", type=str, required=True,
                        help="Path to the trained model output directory")
    parser.add_argument("--skip-strict", action="store_true",
                        help="Skip the strict-branch retraining (Check 3)")
    parser.add_argument("--strict-branch", type=str, default="AB",
                        choices=["AB", "P1"],
                        help="Which branch to hold out for Check 3")
    args = parser.parse_args()

    output_dir = args.output

    # ── Load model and data ───────────────────────────────────────────────
    model, config, data, device, out_dir = load_model_and_data(output_dir)
    print(f"Loaded model from {output_dir}")
    print(f"  Train cells: {len(data['train_cells'])}, "
          f"Val cells: {len(data['val_cells'])}")

    prot_dist = precompute_protein_distances(data["prot_all"])

    # ── Check 1: Lineage-distance audit ──────────────────────────────────
    check_lineage_audit(
        data["train_cells"], data["val_cells"],
        config.lineage_path, output_dir,
    )

    # ── Check 2: Permutation test ────────────────────────────────────────
    check_permutation(model, data["X_val"], data["prot_val"], device, output_dir)

    # ── Check 3: Strict-branch split (optional, slow) ────────────────────
    if not args.skip_strict:
        check_strict_branch(config, config.lineage_path, output_dir,
                            branch=args.strict_branch)
    else:
        print("\nCheck 3: [SKIPPED]")

    # ── Check 4: Pair-set check ──────────────────────────────────────────
    check_pair_sets(config, data, prot_dist, output_dir)

    # ── Check 5: Within-val vs cross-split evaluation ────────────────────
    check_cross_split_eval(model, config, data, device, output_dir)

    print("\n" + "=" * 60)
    print("All diagnostics complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
