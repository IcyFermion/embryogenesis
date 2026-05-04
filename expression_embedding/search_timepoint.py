"""
Hyperparameter search for the timepoint embedding model (all 266 features).

Runs a focused set of configs, each with 80 epochs + early stopping.
Reports a ranked summary table and saves the best config's full outputs.
"""

import json
import os
import sys
from dataclasses import asdict, dataclass
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import torch

# Ensure we can import from the bundle
BUNDLE_DIR = Path(__file__).resolve().parent
if str(BUNDLE_DIR) not in sys.path:
    sys.path.insert(0, str(BUNDLE_DIR))

from timepoint_embedding import (  # noqa: E402
    Config,
    load_and_preprocess,
    build_lineage_labels,
    assign_labels_to_samples,
    sublineage_split,
    run_sanity_checks,
    train,
    extract_all_embeddings,
    _resolve_device,
)

RESULTS_DIR = BUNDLE_DIR / "results" / "timepoint_search"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Shared data loading (done once) ─────────────────────────────────────────

print("Loading data (shared across all configs)...")
base_config = Config(use_all_features=True, n_epochs=80, patience=20, seed=42)
X, sample_meta, feature_names, feat_mean, feat_std = load_and_preprocess(base_config)
(cell_type_one_hot, cell_types, terminal_nodes,
 intermediate_nodes, descendant_list_dict, lineage_data, name_to_type) = \
    build_lineage_labels(base_config.lineage_path, base_config.cell_type_path)
n_classes = len(cell_types)
y_full, hard_mask, sample_weights = assign_labels_to_samples(
    sample_meta, cell_type_one_hot, n_classes, base_config)

# Smoothness pairs
next_t_idx = np.full(len(sample_meta), -1, dtype=int)
cell_to_indices = {}
for i, (cn, t) in enumerate(zip(sample_meta["cell_name"], sample_meta["time"])):
    cell_to_indices[(cn, t)] = i
for (cn, t), i in cell_to_indices.items():
    next_i = cell_to_indices.get((cn, t + 1), -1)
    if next_i >= 0:
        next_t_idx[i] = next_i

# Train/val split (same for all configs)
train_idx, val_idx = sublineage_split(
    base_config.lineage_path, sample_meta,
    base_config.sublineage_depth, base_config.val_fraction, base_config.seed)

X_train, X_val = X[train_idx], X[val_idx]
y_train, y_val = y_full[train_idx], y_full[val_idx]
w_train, w_val = sample_weights[train_idx], sample_weights[val_idx]
hard_train, hard_val = hard_mask[train_idx], hard_mask[val_idx]

def _remap(subset_idx, global_next):
    g2l = {g: l for l, g in enumerate(subset_idx)}
    out = np.full(len(subset_idx), -1, dtype=int)
    for li, gi in enumerate(subset_idx):
        nxt = global_next[gi]
        if nxt >= 0 and nxt in g2l:
            out[li] = g2l[nxt]
    return out

next_train = _remap(train_idx, next_t_idx)
next_val = _remap(val_idx, next_t_idx)

print(f"Train: {len(train_idx)} samples, Val: {len(val_idx)} samples")
print(f"Features: {base_config.n_features}, Classes: {n_classes}")

# ── Search configs ──────────────────────────────────────────────────────────

@dataclass
class SearchEntry:
    name: str
    hidden_dims: tuple
    dropout: float
    alpha: float    # reconstruction
    beta: float     # classification
    gamma: float    # smoothness
    lr: float
    weight_decay: float
    batch_size: int = 256


SEARCH_CONFIGS = [
    # ── Baseline reference ──
    SearchEntry("baseline", (128, 64, 32), 0.1, 1.0, 1.0, 0.1, 1e-3, 1e-4),

    # ── Architecture sweep ──
    SearchEntry("deeper_4layer", (256, 128, 64, 32), 0.1, 1.0, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("wider_256", (256, 128, 64), 0.1, 1.0, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("narrow_2layer", (64, 32), 0.1, 1.0, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("bottleneck", (256, 128, 64, 64), 0.1, 1.0, 1.0, 0.1, 1e-3, 1e-4),

    # ── Regularization sweep ──
    SearchEntry("dropout_0.2", (128, 64, 32), 0.2, 1.0, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("dropout_0.3", (128, 64, 32), 0.3, 1.0, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("dropout_0.2_wd_1e-3", (128, 64, 32), 0.2, 1.0, 1.0, 0.1, 1e-3, 1e-3),

    # ── Loss coefficient sweep ──
    SearchEntry("gamma_0.5", (128, 64, 32), 0.1, 1.0, 1.0, 0.5, 1e-3, 1e-4),
    SearchEntry("gamma_1.0", (128, 64, 32), 0.1, 1.0, 1.0, 1.0, 1e-3, 1e-4),
    SearchEntry("alpha_0.3", (128, 64, 32), 0.1, 0.3, 1.0, 0.1, 1e-3, 1e-4),
    SearchEntry("beta_3.0", (128, 64, 32), 0.1, 1.0, 3.0, 0.1, 1e-3, 1e-4),
    SearchEntry("classify_focus", (128, 64, 32), 0.1, 0.1, 5.0, 0.3, 1e-3, 1e-4),

    # ── Learning rate sweep ──
    SearchEntry("lr_3e-4", (128, 64, 32), 0.1, 1.0, 1.0, 0.1, 3e-4, 1e-4),

    # ── Combined best guesses ──
    SearchEntry("combo_A", (256, 128, 64), 0.2, 1.0, 1.0, 0.5, 3e-4, 1e-3),
    SearchEntry("combo_B", (128, 64, 32), 0.2, 0.3, 2.0, 0.5, 1e-3, 1e-3),
    SearchEntry("combo_C", (256, 128, 64, 32), 0.2, 1.0, 1.0, 1.0, 3e-4, 1e-3),
]


# ── Run search ──────────────────────────────────────────────────────────────

results = []

for i, entry in enumerate(SEARCH_CONFIGS):
    print(f"\n{'=' * 50}")
    print(f"[{i+1}/{len(SEARCH_CONFIGS)}] {entry.name}")
    print(f"  hidden={entry.hidden_dims}, drop={entry.dropout}, "
          f"α={entry.alpha}, β={entry.beta}, γ={entry.gamma}, "
          f"lr={entry.lr}, wd={entry.weight_decay}")
    print(f"{'=' * 50}")

    cfg = Config(
        use_all_features=True,
        n_features=base_config.n_features,
        hidden_dims=entry.hidden_dims,
        dropout=entry.dropout,
        alpha=entry.alpha,
        beta=entry.beta,
        gamma=entry.gamma,
        lr=entry.lr,
        weight_decay=entry.weight_decay,
        batch_size=entry.batch_size,
        n_epochs=80,
        patience=20,
        seed=42,
    )

    try:
        model, history = train(
            cfg, X_train, y_train, w_train, next_train,
            X_val, y_val, w_val, next_val,
            hard_train, hard_val, n_classes,
        )

        best_epoch = np.argmin(history["val_total"])
        best_val_loss = history["val_total"][best_epoch]
        train_acc_best = history["train_acc"][best_epoch]
        val_acc_best = history["val_acc"][best_epoch]

        # Compute soft-label metrics on val set at best epoch
        device = _resolve_device(cfg.device)
        model = model.to(device)
        model.eval()
        X_val_t = torch.FloatTensor(X_val).to(device)  # noqa: F821
        y_val_t = torch.FloatTensor(y_val).to(device)  # noqa: F821
        with torch.no_grad():
            _, _, logits = model(X_val_t)
            probs = torch.nn.functional.softmax(logits, dim=-1).cpu().numpy()
            target_prob = (probs * y_val).sum(axis=-1).mean()
            soft_ce = -(y_val * np.log(np.clip(probs, 1e-10, 1.0))).sum(axis=-1).mean()
        model.cpu()

        result = {
            "name": entry.name,
            "hidden_dims": str(entry.hidden_dims),
            "dropout": entry.dropout,
            "alpha": entry.alpha,
            "beta": entry.beta,
            "gamma": entry.gamma,
            "lr": entry.lr,
            "weight_decay": entry.weight_decay,
            "epochs_run": len(history["train_total"]),
            "best_epoch": int(best_epoch),
            "train_total_loss": float(history["train_total"][best_epoch]),
            "val_total_loss": float(best_val_loss),
            "train_hard_acc": float(train_acc_best),
            "val_hard_acc": float(val_acc_best),
            "val_target_prob": float(target_prob),
            "val_soft_ce": float(soft_ce),
            "overfit_gap": float(train_acc_best - val_acc_best),
        }
        print(f"  → val_acc={val_acc_best:.4f}, train_acc={train_acc_best:.4f}, "
              f"gap={result['overfit_gap']:.4f}, epochs={result['epochs_run']}", flush=True)

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        result = {"name": entry.name, "error": str(e)}

    results.append(result)

# ── Rank and save ───────────────────────────────────────────────────────────

results_df = pd.DataFrame(results)
results_df = results_df.sort_values("val_hard_acc", ascending=False, na_position="last")

print(f"\n{'=' * 80}")
print("RANKED RESULTS (by val hard-label accuracy)")
print(f"{'=' * 80}")

display_cols = [
    "name", "val_hard_acc", "train_hard_acc", "overfit_gap",
    "val_target_prob", "val_soft_ce", "epochs_run",
    "hidden_dims", "dropout", "alpha", "beta", "gamma", "lr", "weight_decay",
]
available_cols = [c for c in display_cols if c in results_df.columns]
print(results_df[available_cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

csv_path = RESULTS_DIR / "search_results.csv"
results_df.to_csv(str(csv_path), index=False)
print(f"\nSaved → {csv_path}")

# ── Retrain best config with full epochs ────────────────────────────────────

best_name = results_df.iloc[0]["name"]
best_entry = next(e for e in SEARCH_CONFIGS if e.name == best_name)
print(f"\n{'=' * 60}")
print(f"Best config: {best_name}")
print(f"Retraining with 200 epochs + patience=30...")
print(f"{'=' * 60}")

best_cfg = Config(
    use_all_features=True,
    n_features=base_config.n_features,
    hidden_dims=best_entry.hidden_dims,
    dropout=best_entry.dropout,
    alpha=best_entry.alpha,
    beta=best_entry.beta,
    gamma=best_entry.gamma,
    lr=best_entry.lr,
    weight_decay=best_entry.weight_decay,
    batch_size=best_entry.batch_size,
    n_epochs=200,
    patience=30,
    seed=42,
    output_dir=str(RESULTS_DIR / f"best_{best_name}"),
)

model, history = train(
    best_cfg, X_train, y_train, w_train, next_train,
    X_val, y_val, w_val, next_val,
    hard_train, hard_val, n_classes,
)

best_epoch = np.argmin(history["val_total"])
print(f"  Best epoch: {best_epoch}, val_loss={history['val_total'][best_epoch]:.4f}")
print(f"  Train hard acc: {history['train_acc'][best_epoch]:.4f}")
print(f"  Val hard acc:   {history['val_acc'][best_epoch]:.4f}")

# Save best model
output_dir = Path(best_cfg.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
torch.save({
    "model_state_dict": model.state_dict(),
    "config": {k: v for k, v in best_cfg.__dict__.items()},
    "feature_names": feature_names,
    "cell_types": cell_types,
    "feat_mean": feat_mean,
    "feat_std": feat_std,
    "history": history,
}, str(output_dir / "model_checkpoint.pt"))

# Extract and save embeddings
embeddings = extract_all_embeddings(model, X, best_cfg.device)
emb_cols = [f"emb_{i}" for i in range(embeddings.shape[1])]
emb_df = pd.DataFrame(embeddings, columns=emb_cols)
emb_df.insert(0, "time", sample_meta["time"].values)
emb_df.insert(0, "cell_name", sample_meta["cell_name"].values)
emb_df.to_csv(str(output_dir / "embeddings.csv"), index=False)

cell_means = emb_df.groupby("cell_name", observed=True)[emb_cols].mean()
cell_means.to_csv(str(output_dir / "cell_embeddings_mean.csv"))

print(f"Saved best model + embeddings → {output_dir}")
print("\nDone.")
