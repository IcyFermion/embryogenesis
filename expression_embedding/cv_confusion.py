"""
cv_confusion.py
===============
Cross-validated confusion matrix for the best Phase 2 protein feature
selection configuration.  Reuses the exact data pipeline from
protein_feature_select.ipynb and collects hard-label predictions across all
CV folds to show where the model confuses cell types.

Usage:
    python3 expression_embedding/cv_confusion.py
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
import torch.nn.functional as F
from sklearn.metrics import confusion_matrix

mpl.rcParams["figure.dpi"] = 200

# ── Paths (same as protein_feature_select.ipynb) ────────────────────────────

BUNDLE_DIR = Path(__file__).resolve().parent
ROOT_DIR = BUNDLE_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
if str(BUNDLE_DIR) not in sys.path:
    sys.path.insert(0, str(BUNDLE_DIR))

RESULTS_DIR = BUNDLE_DIR / "results"
DATA_DIR = ROOT_DIR / "data"
PROTEIN_DATA_DIR = DATA_DIR / "protein" / "aggregated_all"
CELL_LINEAGE_PATH = DATA_DIR / "cell_lineage.json"
CELL_TYPE_PATH = DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv"
S3_PATH = PROTEIN_DATA_DIR / "s3.csv"

# ── Best config from Phase 2 CV (top row of phase2_cv_summary.csv) ──────────
# n_select=25, hidden_dims=(32,), dropout=0.1, dist_lambda=0.0
BEST_SELECTOR_CONFIG = {
    "l1_lambda": 0.002,
    "hidden_dims": (128, 64),
    "dropout": 0.2,
}
BEST_PHASE2_CONFIG = {
    "n_select": 25,
    "hidden_dims": (32,),
    "dropout": 0.1,
    "dist_lambda": 0.0,
}

N_FOLDS = 5
SELECTOR_EPOCHS = 300
FOCUSED_EPOCHS = 400
FOCUSED_BATCH_SIZE = 64
SEED = 42

# ── Device ──────────────────────────────────────────────────────────────────

device = (
    torch.device("mps")
    if torch.backends.mps.is_available()
    else torch.device("cuda")
    if torch.cuda.is_available()
    else torch.device("cpu")
)
print(f"Device: {device}")


# ── Reuse the name mapping ──────────────────────────────────────────────────

def map_names(did):
    if did == "P4a":    return "Z3"
    elif did == "P4p":  return "Z2"
    elif did == "P0a":  return "AB"
    return did


# ── 1. Load lineage & build labels (identical to protein_feature_select) ────

print("Building lineage labels ...")
with CELL_LINEAGE_PATH.open("r", encoding="utf-8") as f:
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
print(f"  {len(intermediate_nodes)} intermediate, {len(terminal_nodes)} terminal")

# Cell type mapping
cell_type_df = pd.read_csv(CELL_TYPE_PATH)
cell_type_dict = {}
for node in terminal_nodes:
    cur = cell_type_df[cell_type_df["wormweb.lineage"] == node]
    cur_types = cur["wormweb.type"].dropna().unique()
    cell_type_dict[node] = cur_types[0] if len(cur_types) > 0 else "programmed_death"

# Build ordered class list
cell_types = sorted(set(cell_type_dict.values()))
if "programmed_death" in cell_types:
    cell_types.remove("programmed_death")
cell_types = sorted(cell_types, key=lambda x: sum(1 for v in cell_type_dict.values() if v == x), reverse=True)
cell_types.append("programmed_death")
cell_type_to_int = {ct: i for i, ct in enumerate(cell_types)}

# One-hot for terminals, soft for intermediates
cell_type_one_hot = {}
for node, ct in cell_type_dict.items():
    oh = np.zeros(len(cell_types), dtype=np.float32)
    oh[cell_type_to_int[ct]] = 1.0
    cell_type_one_hot[node] = oh

for node in intermediate_nodes:
    desc = descendant_list_dict.get(node, [])
    summed = np.sum([cell_type_one_hot[d] for d in desc], axis=0) if desc else np.zeros(len(cell_types))
    cell_type_one_hot[node] = (summed / summed.sum()).astype(np.float32) if summed.sum() > 0 else np.ones(len(cell_types)) / len(cell_types)

# ── 2. Load expression data ─────────────────────────────────────────────────

print("Loading expression data ...")
protein_exp = pd.read_csv(S3_PATH, index_col=0).T.fillna(0)
# Drop all-zero genes
protein_exp = protein_exp[~(protein_exp == 0).all(axis=1)]
# z-score per gene
protein_exp_zscore = protein_exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

feature_names = np.array(protein_exp_zscore.index.tolist())
X = protein_exp_zscore.values.T.astype(np.float32)
y = np.array([cell_type_one_hot[map_names(node)] for node in protein_exp_zscore.columns], dtype=np.float32)

# Sample weights (entropy-based, same as notebook)
y_clipped = np.clip(y, 1e-10, 1.0)
entropies = -np.sum(y_clipped * np.log(y_clipped), axis=1)
norm_ent = entropies / np.log(len(cell_types))
sample_weights = np.exp(-3.0 * norm_ent)
sample_weights = (sample_weights / sample_weights.mean()).astype(np.float32)

terminal_mask = np.array([
    len(descendant_list_dict.get(map_names(node), [])) == 0
    for node in protein_exp_zscore.columns
], dtype=bool)

cell_names = np.array([map_names(n) for n in protein_exp_zscore.columns])

print(f"  X: {X.shape}, y: {y.shape}")
print(f"  Terminal (hard-label): {terminal_mask.sum()}, Intermediate (soft): {(~terminal_mask).sum()}")
print(f"  Classes: {len(cell_types)}")


# ── 3. Per-fold training & prediction ──────────────────────────────────────

print(f"\nRunning {N_FOLDS}-fold CV with best Phase 2 config ...")

# Import from the bundle copy (not the root duplicate)
import importlib
feature_selection = importlib.import_module("feature_selection")
train_one_pass = feature_selection.train_one_pass
train_focused = feature_selection.train_focused
_jaccard = feature_selection._jaccard

from sklearn.model_selection import StratifiedKFold

strat_labels = y.argmax(axis=1)
splitter = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
folds = list(splitter.split(X, strat_labels))

all_probs = np.zeros((len(X), len(cell_types)), dtype=np.float32)
fold_assignment = np.full(len(X), -1, dtype=int)
selected_features_per_fold = []
fold_accuracies = []

for fold_idx, (train_idx, val_idx) in enumerate(folds):
    print(f"\n  Fold {fold_idx + 1}/{N_FOLDS} "
          f"(train={len(train_idx)}, val={len(val_idx)}) ...")

    # Phase 1: train selector on train fold
    _, top_idx, gate_vals, _, _ = train_one_pass(
        X[train_idx], y[train_idx], sample_weights[train_idx],
        hidden_dims=BEST_SELECTOR_CONFIG["hidden_dims"],
        l1_lambda=BEST_SELECTOR_CONFIG["l1_lambda"],
        dropout=BEST_SELECTOR_CONFIG.get("dropout", 0.2),
        n_select=BEST_PHASE2_CONFIG["n_select"],
        n_epochs=SELECTOR_EPOCHS,
        seed=SEED + fold_idx,
        device=str(device),
    )

    selected = feature_names[top_idx].tolist()
    selected_features_per_fold.append(selected)
    X_tr_sel = X[train_idx][:, top_idx].astype(np.float32)
    X_val_sel = X[val_idx][:, top_idx].astype(np.float32)

    # Phase 2: train focused model on selected features
    model, _, losses = train_focused(
        X_tr_sel, y[train_idx], sample_weights[train_idx],
        hidden_dims=BEST_PHASE2_CONFIG["hidden_dims"],
        dropout=BEST_PHASE2_CONFIG.get("dropout", 0.1),
        n_epochs=FOCUSED_EPOCHS,
        batch_size=FOCUSED_BATCH_SIZE,
        seed=SEED + fold_idx,
        device=str(device),
        dist_lambda=BEST_PHASE2_CONFIG.get("dist_lambda", 0.0),
    )

    # Predict on val fold
    with torch.no_grad():
        logits_val, _ = model(torch.FloatTensor(X_val_sel).to(device))
    val_probs = F.softmax(logits_val, dim=-1).cpu().numpy()
    all_probs[val_idx] = val_probs
    fold_assignment[val_idx] = fold_idx

    # Hard-label accuracy on this fold
    val_hard = terminal_mask[val_idx]
    if val_hard.sum() > 0:
        fold_acc = (val_probs[val_hard].argmax(axis=1) == y[val_idx][val_hard].argmax(axis=1)).mean()
    else:
        fold_acc = float("nan")
    fold_accuracies.append(fold_acc)
    print(f"    Hard-label val acc: {fold_acc:.4f}  |  epochs: {len(losses)}")


# ── 4. Confusion matrix on hard-labeled (terminal) samples ──────────────────

print(f"\nOverall hard-label OOF accuracy: {np.nanmean(fold_accuracies):.4f} "
      f"± {np.nanstd(fold_accuracies):.4f}")

term_idx = np.where(terminal_mask)[0]
term_preds = all_probs[term_idx].argmax(axis=1)
term_trues = y[term_idx].argmax(axis=1)
term_names = cell_names[term_idx]

# Confusion matrix
present_classes = sorted(set(term_trues) | set(term_preds))
class_labels = [cell_types[i] for i in present_classes]
conf = confusion_matrix(term_trues, term_preds, labels=present_classes)
conf_norm = confusion_matrix(term_trues, term_preds, labels=present_classes, normalize="true")

# ── Per-class metrics ───────────────────────────────────────────────────────

per_class = defaultdict(list)
for p, t in zip(term_preds, term_trues):
    per_class[cell_types[t]].append(int(p == t))

print(f"\n{'Cell type':30s}  {'Acc':>6}  {'n':>5}")
print("-" * 45)
for ct, vals in sorted(per_class.items(), key=lambda x: -np.mean(x[1])):
    print(f"  {ct:28s}  {np.mean(vals):.4f}  {len(vals):5d}")

# Most common confusions
print(f"\nTop off-diagonal confusions:")
off_diag = []
for i, true_cls in enumerate(present_classes):
    for j, pred_cls in enumerate(present_classes):
        if i != j and conf[i, j] > 0:
            off_diag.append((cell_types[true_cls], cell_types[pred_cls], conf[i, j]))
off_diag.sort(key=lambda x: -x[2])
for true_cls, pred_cls, count in off_diag[:10]:
    print(f"  {true_cls:28s} → {pred_cls:28s}  ({count} samples)")


# ── 5. Feature stability across folds ───────────────────────────────────────

print(f"\nFeature stability across folds:")
for i in range(N_FOLDS):
    for j in range(i + 1, N_FOLDS):
        jac = _jaccard(selected_features_per_fold[i], selected_features_per_fold[j])
        print(f"  Fold {i+1} vs Fold {j+1}: Jaccard = {jac:.4f}")

# Features selected in ≥4 folds
feature_freq = defaultdict(int)
for sf in selected_features_per_fold:
    for f in set(sf):
        feature_freq[f] += 1
stable = sorted([(f, c) for f, c in feature_freq.items() if c >= 4], key=lambda x: -x[1])
print(f"\nFeatures selected in ≥4/5 folds ({len(stable)}):")
for f, c in stable:
    print(f"  {f:30s}  {c}/5")


# ── 6. Plots ────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# Absolute counts
sns.heatmap(conf, annot=True, fmt="d", cmap="Blues",
            xticklabels=class_labels, yticklabels=class_labels, ax=axes[0],
            linewidths=0.3, cbar_kws={"label": "Count"})
axes[0].set_title("Confusion matrix – terminal (hard-label) cells\n(5-fold CV, best Phase 2 config)")
axes[0].set_xlabel("Predicted")
axes[0].set_ylabel("True")
axes[0].tick_params(axis="x", rotation=45)
axes[0].tick_params(axis="y", rotation=0)

# Row-normalized (recall)
sns.heatmap(conf_norm, annot=True, fmt=".2f", cmap="YlOrRd", vmin=0, vmax=1,
            xticklabels=class_labels, yticklabels=class_labels, ax=axes[1],
            linewidths=0.3, cbar_kws={"label": "Recall"})
axes[1].set_title("Row-normalized confusion matrix (recall per true class)")
axes[1].set_xlabel("Predicted")
axes[1].set_ylabel("True")
axes[1].tick_params(axis="x", rotation=45)
axes[1].tick_params(axis="y", rotation=0)

plt.tight_layout()
out_path = RESULTS_DIR / "cv_confusion_matrix.png"
fig.savefig(str(out_path), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved confusion matrix → {out_path}")

# ── 7. Per-fold accuracy bar chart ──────────────────────────────────────────

fig, ax = plt.subplots(figsize=(7, 4))
bars = ax.bar(range(1, N_FOLDS + 1), fold_accuracies, color="steelblue", edgecolor="white")
ax.axhline(y=np.mean(fold_accuracies), color="tomato", linestyle="--",
           label=f"Mean = {np.mean(fold_accuracies):.3f}")
for bar, acc in zip(bars, fold_accuracies):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
            f"{acc:.3f}", ha="center", fontsize=9)
ax.set_xlabel("Fold")
ax.set_ylabel("Hard-label accuracy")
ax.set_title("Per-fold hard-label accuracy (terminal cells only)")
ax.set_ylim(0, 1)
ax.legend()
ax.set_xticks(range(1, N_FOLDS + 1))

plt.tight_layout()
out_path = RESULTS_DIR / "cv_fold_accuracy.png"
fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved fold accuracy → {out_path}")

# ── 8. Save per-sample predictions ──────────────────────────────────────────

pred_df = pd.DataFrame({
    "cell_name": cell_names,
    "is_terminal": terminal_mask,
    "true_class": [cell_types[y[i].argmax()] for i in range(len(y))],
    "pred_class": [cell_types[all_probs[i].argmax()] for i in range(len(all_probs))],
    "correct": [all_probs[i].argmax() == y[i].argmax() for i in range(len(all_probs))],
    "confidence": all_probs.max(axis=1),
    "fold": fold_assignment,
})
# Add per-class probabilities
for ci, ct in enumerate(cell_types):
    pred_df[f"prob_{ct}"] = all_probs[:, ci]

out_path = RESULTS_DIR / "cv_oof_predictions.csv"
pred_df.to_csv(str(out_path), index=False)
print(f"Saved OOF predictions → {out_path}")

print("\nDone.")
