"""Random hyperparameter sweep for mRNA-to-protein embedding alignment.

Runs N randomly-sampled configurations sequentially, saving per-run metrics
to a summary CSV. Designed for ~7-8 hour overnight run.

Usage:
    python expression_embedding/rna_protein_align/sweep.py --n-runs 1500
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

# Path setup
BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from expression_embedding.rna_protein_align.config import Config
from expression_embedding.rna_protein_align.data_loader import (
    PairSampler,
    load_rna_protein_data,
    precompute_protein_distances,
)
from expression_embedding.rna_protein_align.model import RnaEncoderModel
from expression_embedding.rna_protein_align.trainer import (
    _resolve_device,
    make_fixed_val_pairs,
    pearson_on_pairs,
    train,
)


def sample_config(rng: np.random.RandomState) -> dict:
    """Randomly sample a configuration from the focused parameter space."""
    hidden_dims_opts = [
        (128, 64, 32),
        (256, 128, 64),
        (512, 256, 128),
        (256, 128),
        (512, 256),
    ]
    return {
        "hidden_dims": hidden_dims_opts[rng.randint(0, len(hidden_dims_opts))],
        "near_fraction": float(rng.choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])),
        "beta": float(rng.choice([0.02, 0.05, 0.1, 0.2, 0.5])),
        "sublineage_depth": int(rng.choice([5, 6, 7, 8])),
        "lr": float(10 ** rng.uniform(-3.3, -2.3)),  # log-uniform [5e-4, 5e-3]
        "dropout": float(rng.choice([0.0, 0.05, 0.1])),
        "n_pairs_per_epoch": int(rng.choice([10000, 15000, 20000])),
        "weight_decay": float(rng.choice([1e-5, 1e-4, 5e-4, 1e-3])),
        "seed": int(rng.randint(0, 100000)),
    }


def run_single(config_dict: dict, run_id: int, results_dir: Path) -> dict:
    """Run a single training with the given config overrides. Returns metrics dict."""
    config = Config()

    # Override from sweep config
    config.hidden_dims = config_dict["hidden_dims"]
    config.near_fraction = config_dict["near_fraction"]
    config.beta = config_dict["beta"]
    config.sublineage_depth = config_dict["sublineage_depth"]
    config.lr = config_dict["lr"]
    config.dropout = config_dict["dropout"]
    config.n_pairs_per_epoch = config_dict["n_pairs_per_epoch"]
    config.weight_decay = config_dict.get("weight_decay", 1e-4)
    config.seed = config_dict["seed"]
    config.val_n_pairs = 5000
    config.n_epochs = 300
    config.patience = 50

    run_dir = results_dir / f"run_{run_id:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    config.output_dir = str(run_dir)

    try:
        # Load data
        data = load_rna_protein_data(
            rna_path=config.rna_path,
            protein_path=config.protein_emb_path,
            lineage_path=config.lineage_path,
            sublineage_depth=config.sublineage_depth,
            val_fraction=config.val_fraction,
            seed=config.seed,
            log_transform=config.log_transform,
        )

        prot_dist = precompute_protein_distances(data["prot_all"])
        train_samp = PairSampler(
            data["train_indices"], prot_dist, config.near_fraction, config.seed,
        )
        val_samp = PairSampler(
            data["val_indices"], prot_dist, 0.0, config.seed + 1,
        )

        # Train
        model, history = train(
            config, data["X_train"], data["X_val"],
            data["prot_train"], data["prot_val"],
            prot_dist, train_samp, val_samp, data["n_features"],
        )

        best_pearson = float(max(history["val_pearson"]))
        best_epoch = int(np.argmax(history["val_pearson"]))
        final_pearson = float(history["val_pearson"][-1])
        epochs_run = len(history["train_total"])

        # Linear probe baseline on val data
        device = _resolve_device(config.device)
        fixed_a, fixed_b = make_fixed_val_pairs(len(data["X_val"]), 5000, seed=config.seed)
        baseline_pearson, _, _ = pearson_on_pairs(
            RnaEncoderModel(data["n_features"], config.hidden_dims, config.embed_dim,
                          config.dropout, config.use_layer_norm).to(device),
            data["X_val"], data["prot_val"], fixed_a, fixed_b, device,
        )

        # Save minimal outputs for this run
        run_metrics = {
            "run_id": run_id,
            "hidden_dims": str(config.hidden_dims),
            "near_fraction": config.near_fraction,
            "beta": config.beta,
            "sublineage_depth": config.sublineage_depth,
            "lr": config.lr,
            "dropout": config.dropout,
            "n_pairs_per_epoch": config.n_pairs_per_epoch,
            "seed": config.seed,
            "n_train_cells": len(data["train_cells"]),
            "n_val_cells": len(data["val_cells"]),
            "best_val_pearson": best_pearson,
            "best_epoch": best_epoch,
            "final_val_pearson": final_pearson,
            "epochs_run": epochs_run,
            "baseline_pearson": float(baseline_pearson),
            "best_val_total": float(min(history["val_total"])),
            "best_train_align": float(min(history["train_align"])),
            "best_val_align": float(min(history["val_align"])),
            "best_train_recon": float(min(history["train_recon"])),
            "best_val_recon": float(min(history["val_recon"])),
        }

        # Save per-run CSV
        pd.DataFrame(history).to_csv(run_dir / "training_curves.csv", index=False)

    except Exception as e:
        run_metrics = {
            "run_id": run_id,
            "error": str(e),
            **{k: config_dict.get(k, None) for k in [
                "hidden_dims", "near_fraction", "beta", "sublineage_depth",
                "lr", "dropout", "n_pairs_per_epoch", "seed",
            ]},
        }

    return run_metrics


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Hyperparameter sweep")
    parser.add_argument("--n-runs", type=int, default=1500,
                        help="Number of random configurations to try")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory for sweep results")
    args = parser.parse_args()

    default_out = BUNDLE_DIR / "results" / "rna_protein_align_sweep"
    results_dir = Path(args.output) if args.output else default_out
    results_dir.mkdir(parents=True, exist_ok=True)

    summary_path = results_dir / "summary.csv"

    # Determine start index from existing runs
    start_id = 0
    if summary_path.exists():
        existing = pd.read_csv(summary_path)
        if len(existing) > 0:
            start_id = int(existing["run_id"].max()) + 1
            print(f"Resuming from run {start_id} ({len(existing)} existing runs)")

    rng = np.random.RandomState(42 + start_id)

    n_runs = args.n_runs
    print(f"Starting sweep: {n_runs} runs from run_id={start_id}")
    print(f"Output: {summary_path}")
    print(f"Estimated time: ~{n_runs * 15 / 3600:.1f} hours")
    print("=" * 60)

    t_start = time.time()

    for i in range(n_runs):
        run_id = start_id + i
        t_run = time.time()

        config_dict = sample_config(rng)
        config_dict["seed"] = rng.randint(0, 100000)

        metrics = run_single(config_dict, run_id, results_dir)

        # Append to summary CSV
        write_header = not summary_path.exists() or (i == 0 and start_id == 0)
        with open(summary_path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=sorted(metrics.keys()))
            if write_header:
                writer.writeheader()
            writer.writerow(metrics)

        elapsed_run = time.time() - t_run
        elapsed_total = time.time() - t_start
        eta = (elapsed_total / (i + 1)) * (n_runs - i - 1)

        best = metrics.get("best_val_pearson", float("nan"))
        status = f"best_pearson={best:.4f}" if "error" not in metrics else f"ERROR: {metrics['error'][:60]}"

        print(f"  [{run_id:4d}] {status} | "
              f"{elapsed_run:.1f}s | "
              f"elapsed={elapsed_total/60:.0f}m | eta={eta/60:.0f}m")

    # Final summary
    elapsed_total = time.time() - t_start
    print("=" * 60)
    print(f"Sweep complete. {n_runs} runs in {elapsed_total/60:.1f} minutes.")

    if summary_path.exists():
        df = pd.read_csv(summary_path)
        if "best_val_pearson" not in df.columns:
            print("\nAll runs errored — no metrics collected.")
            return
        valid = df[df["best_val_pearson"].notna()].copy()
        if len(valid) > 0:
            valid = valid.sort_values("best_val_pearson", ascending=False)
            print(f"\nTop 10 configurations:")
            cols = ["run_id", "best_val_pearson", "best_epoch", "hidden_dims",
                    "near_fraction", "beta", "sublineage_depth", "lr", "dropout",
                    "n_pairs_per_epoch", "n_train_cells", "n_val_cells"]
            available = [c for c in cols if c in valid.columns]
            print(valid[available].head(10).to_string(index=False))

            top10_path = results_dir / "top10.csv"
            valid.head(10).to_csv(top10_path, index=False)
            print(f"\nTop 10 saved → {top10_path}")

    print(f"\nFull summary → {summary_path}")


if __name__ == "__main__":
    main()
