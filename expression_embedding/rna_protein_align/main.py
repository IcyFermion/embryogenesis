"""Orchestration and CLI entry point for mRNA-to-protein embedding alignment."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

# Path setup (mirrors timepoint_embedding.py pattern)
BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from .config import Config
from .data_loader import (
    PairSampler,
    load_rna_protein_data,
    precompute_protein_distances,
)
from .evaluator import (
    plot_distance_scatter,
    plot_lineage_distance_scatter,
    plot_sibling_distances,
    plot_umap_embeddings,
    run_sanity_checks,
    save_outputs,
)
from .model import RnaEncoderModel
from .trainer import (
    _resolve_device,
    make_fixed_val_pairs,
    pearson_on_pairs,
    train,
)


def main(config: Config = None):
    if config is None:
        config = Config()

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("mRNA-to-Protein Embedding Alignment (Stage 1)")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Config: alpha={config.alpha}, beta={config.beta}")
    print(f"Hidden dims: {config.hidden_dims}, embed_dim: {config.embed_dim}")
    print(f"Sublineage depth: {config.sublineage_depth}, val frac: {config.val_fraction}")
    print(f"Pairs/epoch: {config.n_pairs_per_epoch}, near fraction: {config.near_fraction}")

    # ── 1. Load and preprocess data ──────────────────────────────────────
    print("\n─── Data Loading ───")
    data = load_rna_protein_data(
        rna_path=config.rna_path,
        protein_path=config.protein_emb_path,
        lineage_path=config.lineage_path,
        sublineage_depth=config.sublineage_depth,
        val_fraction=config.val_fraction,
        seed=config.seed,
        log_transform=config.log_transform,
    )

    X_train = data["X_train"]
    X_val = data["X_val"]
    prot_train = data["prot_train"]
    prot_val = data["prot_val"]
    train_cells = data["train_cells"]
    val_cells = data["val_cells"]
    all_cells = data["all_cells"]
    X_all = data["X_all"]
    prot_all = data["prot_all"]
    gene_names = data["gene_names"]
    n_features = data["n_features"]
    train_indices = data["train_indices"]
    val_indices = data["val_indices"]

    device = _resolve_device(config.device)
    print(f"  Device: {device}")

    # ── 2. Precompute protein cosine distance matrix ─────────────────────
    print("\n─── Precomputing protein distances ───")
    prot_dist_matrix = precompute_protein_distances(prot_all)
    print(f"  Full matrix: {prot_dist_matrix.shape}, "
          f"range [{prot_dist_matrix.min():.4f}, {prot_dist_matrix.max():.4f}]")

    # ── 3. Build pair samplers ───────────────────────────────────────────
    train_pair_sampler = PairSampler(
        cell_indices=train_indices,
        full_prot_dist=prot_dist_matrix,
        near_fraction=config.near_fraction,
        seed=config.seed,
    )
    val_pair_sampler = PairSampler(
        cell_indices=val_indices,
        full_prot_dist=prot_dist_matrix,
        near_fraction=0.0,  # uniform only for validation
        seed=config.seed + 1,
    )

    # ── 4. Instantiate model for sanity checks ───────────────────────────
    model = RnaEncoderModel(
        n_features=n_features,
        hidden_dims=config.hidden_dims,
        embed_dim=config.embed_dim,
        dropout=config.dropout,
        use_layer_norm=config.use_layer_norm,
    ).to(device)

    # ── 5. Pre-training sanity checks ────────────────────────────────────
    baseline_metrics = run_sanity_checks(
        X_train=X_train,
        X_val=X_val,
        prot_train=prot_train,
        prot_val=prot_val,
        train_cells=train_cells,
        val_cells=val_cells,
        model=model,
        device=device,
        config=config,
    )

    # ── 6. Train ─────────────────────────────────────────────────────────
    print("\n─── Training ───")
    model, history = train(
        config=config,
        X_train=X_train,
        X_val=X_val,
        prot_train=prot_train,
        prot_val=prot_val,
        prot_dist_matrix=prot_dist_matrix,
        train_pair_sampler=train_pair_sampler,
        val_pair_sampler=val_pair_sampler,
        n_features=n_features,
    )

    best_epoch = int(np.argmax(history["val_pearson"]))
    best_pearson = history["val_pearson"][best_epoch]
    print(f"\n  Best val Pearson: {best_pearson:.4f} at epoch {best_epoch}")

    # ── 7. Extract final embeddings for ALL overlapping cells ────────────
    print("\n─── Extracting embeddings ───")
    model = model.to(device)
    model.eval()
    with torch.no_grad():
        X_all_t = torch.as_tensor(X_all, dtype=torch.float32, device=device)
        z_all = model.encoder(X_all_t)
        all_embeddings = z_all.cpu().numpy()
    print(f"  Embeddings shape: {all_embeddings.shape}")

    # ── 8. Validate Pearson on final model (fresh held-out pairs) ───────
    fresh_a, fresh_b = make_fixed_val_pairs(
        len(X_val), config.val_n_pairs, seed=config.seed + 1000
    )
    val_pearson, mrna_dists, prot_dists = pearson_on_pairs(
        model, X_val, prot_val, fresh_a, fresh_b, device,
    )
    print(f"  Final val Pearson: {val_pearson:.4f}")

    # Pearson check: warn if out of expected range
    if val_pearson < 0.3:
        print("  [WARN] Val Pearson < 0.3 — investigate training or data")
    elif val_pearson > 0.6:
        print("  [WARN] Val Pearson > 0.6 — suspect data leakage")

    # ── 9. Save outputs ──────────────────────────────────────────────────
    print("\n─── Saving outputs ───")
    baseline_metrics["n_train_cells"] = len(train_cells)
    baseline_metrics["n_val_cells"] = len(val_cells)
    baseline_metrics["n_total_cells"] = len(all_cells)
    save_outputs(
        model=model,
        config=config,
        history=history,
        embeddings=all_embeddings,
        all_cells=all_cells,
        gene_names=gene_names,
        baseline_metrics=baseline_metrics,
        output_dir=str(output_dir),
    )

    # ── 10. Diagnostic plots ─────────────────────────────────────────────
    print("\n─── Diagnostic plots ───")

    # Distance scatter on held-out pairs
    plot_distance_scatter(mrna_dists, prot_dists, str(output_dir), title_suffix="val")

    # Sibling distances
    plot_sibling_distances(
        all_embeddings, all_cells,
        config.lineage_path, config.cell_type_path,
        str(output_dir),
    )

    # Lineage distance scatter
    plot_lineage_distance_scatter(
        all_embeddings, all_cells,
        config.lineage_path, str(output_dir),
    )

    # PCA/UMAP
    plot_umap_embeddings(
        all_embeddings, all_cells,
        config.lineage_path, config.cell_type_path,
        str(output_dir),
    )

    # ── 11. Final summary ────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Done!")
    print(f"  Linear probe baseline R²: {baseline_metrics['linear_probe_r2']:.4f}")
    print(f"  Best val Pearson:         {best_pearson:.4f} (epoch {best_epoch})")
    print(f"  Encoder improvement:       +{best_pearson - np.sqrt(baseline_metrics['linear_probe_r2']):.4f} r")
    print(f"  All outputs in: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="mRNA-to-protein embedding alignment (Stage 1)"
    )
    parser.add_argument("--epochs", type=int, default=None, help="Override n_epochs")
    parser.add_argument("--lr", type=float, default=None, help="Override learning rate")
    parser.add_argument(
        "--hidden-dims", type=str, default=None,
        help="Comma-separated hidden dims, e.g. '128,64'",
    )
    parser.add_argument(
        "--alpha", type=float, default=None, help="Alignment loss weight"
    )
    parser.add_argument(
        "--beta", type=float, default=None, help="Reconstruction loss weight"
    )
    parser.add_argument(
        "--near-fraction", type=float, default=None,
        help="Fraction of pairs from near protein-space quartile",
    )
    parser.add_argument(
        "--output", type=str, default=None, help="Override output directory"
    )
    parser.add_argument(
        "--depth", type=int, default=None, help="Sublineage depth for split"
    )
    parser.add_argument(
        "--no-log-transform", action="store_true",
        help="Disable log1p transform on mRNA data",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for data split and init"
    )
    parser.add_argument(
        "--dropout", type=float, default=None, help="Dropout rate"
    )
    parser.add_argument(
        "--weight-decay", type=float, default=None, help="Weight decay for AdamW"
    )
    parser.add_argument(
        "--patience", type=int, default=None, help="Early stopping patience"
    )
    args = parser.parse_args()

    config = Config()
    if args.epochs is not None:
        config.n_epochs = args.epochs
    if args.lr is not None:
        config.lr = args.lr
    if args.hidden_dims is not None:
        config.hidden_dims = tuple(int(x.strip()) for x in args.hidden_dims.split(","))
    if args.alpha is not None:
        config.alpha = args.alpha
    if args.beta is not None:
        config.beta = args.beta
    if args.near_fraction is not None:
        config.near_fraction = args.near_fraction
    if args.output is not None:
        config.output_dir = args.output
    if args.depth is not None:
        config.sublineage_depth = args.depth
    if args.no_log_transform:
        config.log_transform = False
    if args.seed is not None:
        config.seed = args.seed
    if args.dropout is not None:
        config.dropout = args.dropout
    if args.weight_decay is not None:
        config.weight_decay = args.weight_decay
    if args.patience is not None:
        config.patience = args.patience

    main(config)
