#!/usr/bin/env python3
"""
phase2_alternative.py
─────────────────────────────────────────────────────────────────────────────
Fundamentally different approach to Phase 2 cell type prediction.

Key changes from the neural network approach:
1. Use ALL 210 features with built-in regularization (no Phase 1 selection)
2. Tree-based models (XGBoost, LightGBM, Random Forest) designed for small tabular data
3. Train on terminal (hard-label) samples only - soft labels add noise
4. Proper stratified CV with per-class analysis
5. Ensemble of diverse model families for robustness

Why this should work better:
- Tree models handle high-dim/low-sample naturally (built-in feature selection)
- No gradient-based optimization → no overfitting to batch noise
- Ensemble of different algorithms → better generalization
- Using all features lets the model discover interactions the gate missed
"""

import numpy as np
import pandas as pd
import json
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import (
    RandomForestClassifier,
    GradientBoostingClassifier,
    ExtraTreesClassifier,
    StackingClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, classification_report
import xgboost as xgb
import lightgbm as lgb
from tqdm import tqdm
import warnings

from bundle_paths import (
    CELL_LINEAGE_PATH,
    CELL_TYPE_CSV,
    PHASE2_ALTERNATIVE_RESULTS_CSV,
    S3_CSV,
    ensure_repo_on_path,
    ensure_results_dir,
)

ensure_repo_on_path()
ensure_results_dir()

warnings.filterwarnings('ignore')

print("="*70)
print("PHASE 2 ALTERNATIVE: Tree-Based & Ensemble Approach")
print("="*70)

# ─── Load Data ────────────────────────────────────────────────────────────────
print("\nLoading data...")

with CELL_LINEAGE_PATH.open('r', encoding='utf-8') as f:
    lineage_data = json.load(f)

def map_names(did):
    if   did == "P4a": return "Z3"
    elif did == "P4p": return "Z2"
    elif did == "P0a": return "AB"
    else: return did

terminal_nodes = []
intermediate_nodes = []
descendant_list_dict = defaultdict(list)

def dfs(node, parent, ancestors=[]):
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

cell_type_df = pd.read_csv(CELL_TYPE_CSV)
protein_exp = pd.read_csv(S3_CSV, index_col=0).T
protein_exp = protein_exp.fillna(0)
protein_exp = protein_exp[(protein_exp != 0).any(axis=1)]
protein_exp_zscore = protein_exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Build cell type mappings
cell_type_dict = {}
for node in terminal_nodes:
    cur_cell_type_df = cell_type_df[cell_type_df['wormweb.lineage'] == node]
    if len(cur_cell_type_df) == 0:
        continue
    cur_lineage_types = cur_cell_type_df["wormweb.type"]
    cur_lineage_types = cur_lineage_types[~cur_lineage_types.isna()].unique()
    if len(cur_lineage_types) == 0:
        cell_type_dict[node] = "programmed_death"
        continue
    cell_type_dict[node] = cur_lineage_types[0]

cell_types = list(set(cell_type_dict.values()))
cell_types = [ct for ct in cell_types if ct != "programmed_death"]
cell_types = sorted(cell_types, key=lambda x: sum([1 for v in cell_type_dict.values() if v == x]), reverse=True)
cell_types.append("programmed_death")
cell_type_to_int = {ct: i for i, ct in enumerate(cell_types)}

# Build one-hot labels
cell_type_one_hot = {}
for node, ct in cell_type_dict.items():
    one_hot = np.zeros(len(cell_types))
    one_hot[cell_type_to_int[ct]] = 1
    cell_type_one_hot[node] = one_hot

for node in intermediate_nodes:
    descendant_types = [cell_type_one_hot[desc] for desc in descendant_list_dict[node]]
    if len(descendant_types) == 0:
        cell_type_one_hot[node] = np.zeros(len(cell_types))
        continue
    summed_types = np.sum(descendant_types, axis=0)
    cell_type_one_hot[node] = summed_types / np.sum(summed_types)

# Full dataset
X_all = protein_exp_zscore.values.T
y_all = np.array([cell_type_one_hot[map_names(node)] for node in protein_exp_zscore.columns])
terminal_mask = np.array([
    len(descendant_list_dict[map_names(node)]) == 0
    for node in protein_exp_zscore.columns
])
feature_names = np.array(protein_exp_zscore.index.tolist())

# ─── KEY CHANGE: Use terminal (hard-label) samples only for training ──────────
X_term = X_all[terminal_mask]
y_term_onehot = y_all[terminal_mask]
y_term = y_term_onehot.argmax(axis=1)
cell_names_term = np.array([map_names(n) for n in protein_exp_zscore.columns])[terminal_mask]

print(f"\n✓ Data loaded")
print(f"  All samples:      {X_all.shape[0]} × {X_all.shape[1]} features")
print(f"  Terminal samples:  {X_term.shape[0]} (hard-label only)")
print(f"  Classes:           {len(cell_types)}")
print(f"  Features:          {X_all.shape[1]} (ALL proteins, no Phase 1 selection)")

# Class distribution
# Remap labels to be contiguous (some classes may not appear in terminal samples)
# Also drop very rare classes (< n_splits samples) since they can't be stratified
unique_classes_raw = np.unique(y_term)
class_counts = {c: np.sum(y_term == c) for c in unique_classes_raw}
min_samples_per_class = 5  # Need at least 5 for 5-fold CV
keep_classes = [c for c in unique_classes_raw if class_counts[c] >= min_samples_per_class]
keep_mask = np.isin(y_term, keep_classes)
X_term = X_term[keep_mask]
y_term = y_term[keep_mask]
cell_names_term = cell_names_term[keep_mask]

label_remap = {old: new for new, old in enumerate(sorted(np.unique(y_term)))}
y_term = np.array([label_remap[c] for c in y_term])
remap_cell_types = [cell_types[c] for c in sorted(label_remap.keys())]
n_actual_classes = len(np.unique(y_term))

dropped = len(unique_classes_raw) - len(keep_classes)
print(f"  Actual classes in terminal: {n_actual_classes} (dropped {dropped} rare classes with <{min_samples_per_class} samples)")
unique, counts = np.unique(y_term, return_counts=True)
print(f"\nClass distribution (terminal samples):")
for cls_idx, count in sorted(zip(unique, counts), key=lambda x: -x[1]):
    print(f"  {remap_cell_types[cls_idx]:30s}  n={count}")

# ─── Approach 1: Use ALL features (no Phase 1 selection) ─────────────────────
print("\n" + "="*70)
print("APPROACH 1: All 210 features with tree-based models")
print("="*70)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    'XGBoost': xgb.XGBClassifier(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.6,
        min_child_weight=3,
        reg_alpha=0.1,
        reg_lambda=1.0,
        random_state=42,
        tree_method='hist',
        device='cuda',
        verbosity=0,
        num_class=n_actual_classes,
        objective='multi:softprob',
    ),
    'LightGBM': lgb.LGBMClassifier(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.6,
        min_child_samples=5,
        reg_alpha=0.1,
        reg_lambda=1.0,
        random_state=42,
        verbose=-1,
        device='cpu',
    ),
    'RandomForest': RandomForestClassifier(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=3,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1,
    ),
    'ExtraTrees': ExtraTreesClassifier(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=3,
        max_features='sqrt',
        random_state=42,
        n_jobs=-1,
    ),
    'LogReg_L1': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegression(
            C=1.0, penalty='l1', solver='saga',
            max_iter=5000, random_state=42,
        )),
    ]),
    'LogReg_L2': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegression(
            C=1.0, penalty='l2', solver='lbfgs',
            max_iter=5000, random_state=42,
        )),
    ]),
}

results_all_features = {}

print(f"\nRunning 5-fold CV on {len(models)} models with ALL {X_term.shape[1]} features...\n")

for name, model in models.items():
    fold_accs = []
    fold_preds = np.zeros_like(y_term)

    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X_term, y_term)):
        try:
            model_clone = type(model)(**model.get_params()) if not isinstance(model, Pipeline) else Pipeline([
                (name, type(step).__call__(type(step), **step.get_params()) if hasattr(step, 'get_params') else step)
                for name, step in model.steps
            ])
        except Exception:
            from sklearn.base import clone
            model_clone = clone(model)

        X_tr, X_val = X_term[train_idx], X_term[val_idx]
        y_tr, y_val = y_term[train_idx], y_term[val_idx]

        model_clone.fit(X_tr, y_tr)
        preds = model_clone.predict(X_val)
        acc = accuracy_score(y_val, preds)
        fold_accs.append(acc)
        fold_preds[val_idx] = preds

    mean_acc = np.mean(fold_accs)
    std_acc = np.std(fold_accs)
    results_all_features[name] = {
        'mean_acc': mean_acc,
        'std_acc': std_acc,
        'fold_accs': fold_accs,
        'predictions': fold_preds,
    }
    print(f"  {name:20s}  {mean_acc:.4f} ± {std_acc:.4f}  folds: {[f'{a:.3f}' for a in fold_accs]}")

# ─── Approach 2: Use Phase-1-selected 25 features ────────────────────────────
print("\n" + "="*70)
print("APPROACH 2: Phase-1-selected 25 features with tree-based models")
print("="*70)

# Load the previously selected features
from feature_selection import _train_sparse_gate, _resolve_device
import torch
from sklearn.base import clone

device = _resolve_device('auto')
selector_model, _ = _train_sparse_gate(
    X_all, y_all,
    np.ones(len(X_all), dtype=np.float32),  # uniform weights for selection
    hidden_dims=(128, 64),
    l1_lambda=0.004,
    dropout=0.2,
    n_epochs=300,
    batch_size=128,
    seed=42,
    device=device,
)
gate_vals = selector_model.feature_importance()
top_25 = np.argsort(gate_vals)[::-1][:25]
X_term_25 = X_term[:, top_25]

print(f"Using {len(top_25)} selected features...")

results_25_features = {}

for name, model in models.items():
    fold_accs = []
    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X_term_25, y_term)):
        from sklearn.base import clone
        model_clone = clone(model)
        X_tr, X_val = X_term_25[train_idx], X_term_25[val_idx]
        y_tr, y_val = y_term[train_idx], y_term[val_idx]
        model_clone.fit(X_tr, y_tr)
        preds = model_clone.predict(X_val)
        fold_accs.append(accuracy_score(y_val, preds))

    mean_acc = np.mean(fold_accs)
    std_acc = np.std(fold_accs)
    results_25_features[name] = {'mean_acc': mean_acc, 'std_acc': std_acc, 'fold_accs': fold_accs}
    print(f"  {name:20s}  {mean_acc:.4f} ± {std_acc:.4f}  folds: {[f'{a:.3f}' for a in fold_accs]}")

# ─── Approach 3: Ensemble of best models ──────────────────────────────────────
print("\n" + "="*70)
print("APPROACH 3: Stacking Ensemble (all features)")
print("="*70)

# Sort by accuracy
sorted_models = sorted(results_all_features.items(), key=lambda x: -x[1]['mean_acc'])
top_3_names = [name for name, _ in sorted_models[:3]]
print(f"Top 3 base models: {top_3_names}")

# Simple voting ensemble
print("\nVoting ensemble (majority vote of top models)...")
top_preds = np.array([results_all_features[name]['predictions'] for name in top_3_names])
# Majority vote
from scipy.stats import mode
ensemble_preds = mode(top_preds, axis=0, keepdims=False).mode

ensemble_acc = accuracy_score(y_term, ensemble_preds)
print(f"  Voting ensemble accuracy: {ensemble_acc:.4f}")

# Stacking ensemble with CV
print("\nStacking ensemble (LogReg meta-learner)...")

base_estimators = []
for name in top_3_names:
    base_estimators.append((name, clone(models[name])))

stacking = StackingClassifier(
    estimators=base_estimators,
    final_estimator=LogisticRegression(
        C=1.0, max_iter=5000, random_state=42
    ),
    cv=5,
    n_jobs=-1,
)

stacking_fold_accs = []
for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X_term, y_term)):
    stacking_clone = clone(stacking)
    X_tr, X_val = X_term[train_idx], X_term[val_idx]
    y_tr, y_val = y_term[train_idx], y_term[val_idx]
    stacking_clone.fit(X_tr, y_tr)
    preds = stacking_clone.predict(X_val)
    stacking_fold_accs.append(accuracy_score(y_val, preds))

stacking_acc = np.mean(stacking_fold_accs)
stacking_std = np.std(stacking_fold_accs)
print(f"  Stacking ensemble: {stacking_acc:.4f} ± {stacking_std:.4f}  folds: {[f'{a:.3f}' for a in stacking_fold_accs]}")

# ─── Approach 4: XGBoost hyperparameter tuning ───────────────────────────────
print("\n" + "="*70)
print("APPROACH 4: XGBoost Hyperparameter Tuning (all features)")
print("="*70)

xgb_param_grid = [
    {'n_estimators': 300, 'max_depth': 4, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.5, 'min_child_weight': 1, 'reg_alpha': 0.0, 'reg_lambda': 1.0},
    {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_weight': 3, 'reg_alpha': 0.1, 'reg_lambda': 1.0},
    {'n_estimators': 800, 'max_depth': 8, 'learning_rate': 0.03, 'subsample': 0.7, 'colsample_bytree': 0.5, 'min_child_weight': 5, 'reg_alpha': 0.5, 'reg_lambda': 2.0},
    {'n_estimators': 500, 'max_depth': 4, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.7, 'min_child_weight': 1, 'reg_alpha': 0.0, 'reg_lambda': 0.5},
    {'n_estimators': 1000, 'max_depth': 6, 'learning_rate': 0.02, 'subsample': 0.8, 'colsample_bytree': 0.4, 'min_child_weight': 3, 'reg_alpha': 1.0, 'reg_lambda': 3.0},
    {'n_estimators': 500, 'max_depth': 3, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.8, 'min_child_weight': 1, 'reg_alpha': 0.0, 'reg_lambda': 1.0},
    {'n_estimators': 700, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.85, 'colsample_bytree': 0.6, 'min_child_weight': 2, 'reg_alpha': 0.3, 'reg_lambda': 1.5},
    {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_weight': 3, 'reg_alpha': 0.1, 'reg_lambda': 1.0, 'gamma': 0.1},
]

print(f"Testing {len(xgb_param_grid)} XGBoost configurations...")

xgb_results = []
for i, params in enumerate(xgb_param_grid):
    model = xgb.XGBClassifier(
        **params,
        random_state=42,
        tree_method='hist',
        device='cuda',
        verbosity=0,
        num_class=n_actual_classes,
        objective='multi:softprob',
    )

    fold_accs = []
    for train_idx, val_idx in skf.split(X_term, y_term):
        model_clone = clone(model)
        model_clone.fit(X_term[train_idx], y_term[train_idx])
        preds = model_clone.predict(X_term[val_idx])
        fold_accs.append(accuracy_score(y_term[val_idx], preds))

    mean_acc = np.mean(fold_accs)
    xgb_results.append({
        'config_idx': i,
        'mean_acc': mean_acc,
        'std_acc': np.std(fold_accs),
        'fold_accs': fold_accs,
        **params,
    })
    print(f"  Config {i}: {mean_acc:.4f} ± {np.std(fold_accs):.4f}")

xgb_results.sort(key=lambda x: -x['mean_acc'])
print(f"\nBest XGBoost config: {xgb_results[0]['mean_acc']:.4f}")

# ─── Approach 5: Train on ALL samples (terminal + intermediate with soft labels) ─
print("\n" + "="*70)
print("APPROACH 5: XGBoost with all samples (soft-label augmented)")
print("="*70)

# For tree models, use hard labels from argmax of soft distribution
y_all_hard = y_all.argmax(axis=1)
# Remap to contiguous labels matching terminal set
y_all_hard_remap = np.array([label_remap.get(c, -1) for c in y_all_hard])
# Only keep samples whose dominant class appears in terminal set
valid_mask = y_all_hard_remap >= 0

# Give higher sample weight to terminal (hard-label) nodes
from feature_selection import entropy_based_weights
sw = entropy_based_weights(y_all, alpha=3.0)

best_xgb_params = {k: v for k, v in xgb_results[0].items()
                   if k not in ('config_idx', 'mean_acc', 'std_acc', 'fold_accs')}

fold_accs_aug = []
X_valid = X_all[valid_mask]
y_valid = y_all_hard_remap[valid_mask]
sw_valid = sw[valid_mask]
term_in_valid = terminal_mask[valid_mask]

for train_idx, val_idx in skf.split(X_valid, y_valid):
    # Only evaluate on terminal samples in validation
    val_term_mask = term_in_valid[val_idx]
    if val_term_mask.sum() == 0:
        continue

    model = xgb.XGBClassifier(
        **best_xgb_params,
        random_state=42,
        tree_method='hist',
        device='cuda',
        verbosity=0,
    )
    model.fit(X_valid[train_idx], y_valid[train_idx], sample_weight=sw_valid[train_idx])
    preds = model.predict(X_valid[val_idx][val_term_mask])
    acc = accuracy_score(y_valid[val_idx][val_term_mask], preds)
    fold_accs_aug.append(acc)

aug_acc = np.mean(fold_accs_aug)
aug_std = np.std(fold_accs_aug)
print(f"  XGBoost (all samples, weighted): {aug_acc:.4f} ± {aug_std:.4f}")

# ─── Final Summary ────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

print(f"\n{'Method':50s}  {'Accuracy':>8}  {'±':>6}")
print("-"*70)

# Neural net baseline
print(f"{'NN baseline (Phase 2 CV, 25 features)':50s}  {'0.7574':>8}  {'0.0300':>6}")

# All features results
for name, res in sorted(results_all_features.items(), key=lambda x: -x[1]['mean_acc']):
    print(f"{f'{name} (all {X_term.shape[1]} features)':50s}  {res['mean_acc']:>8.4f}  {res['std_acc']:>6.4f}")

# 25 features results
for name, res in sorted(results_25_features.items(), key=lambda x: -x[1]['mean_acc']):
    print(f"{f'{name} (25 selected features)':50s}  {res['mean_acc']:>8.4f}  {res['std_acc']:>6.4f}")

# Ensemble
print(f"{'Voting ensemble (top 3, all features)':50s}  {ensemble_acc:>8.4f}  {'N/A':>6}")
print(f"{'Stacking ensemble (top 3, all features)':50s}  {stacking_acc:>8.4f}  {stacking_std:>6.4f}")

# XGBoost tuned
print(f"{'XGBoost tuned (all features)':50s}  {xgb_results[0]['mean_acc']:>8.4f}  {xgb_results[0]['std_acc']:>6.4f}")

# Augmented
print(f"{'XGBoost (all samples + soft-label weight)':50s}  {aug_acc:>8.4f}  {aug_std:>6.4f}")

# Find overall best
all_results = [
    ('NN baseline', 0.7574),
    ('Voting ensemble', ensemble_acc),
    ('Stacking ensemble', stacking_acc),
    ('XGBoost tuned', xgb_results[0]['mean_acc']),
    ('XGBoost augmented', aug_acc),
]
for name, res in results_all_features.items():
    all_results.append((f'{name} (all feat)', res['mean_acc']))
for name, res in results_25_features.items():
    all_results.append((f'{name} (25 feat)', res['mean_acc']))

all_results.sort(key=lambda x: -x[1])

print(f"\n{'='*70}")
print(f"BEST OVERALL: {all_results[0][0]} = {all_results[0][1]:.4f}")
improvement = all_results[0][1] - 0.7574
print(f"Improvement over NN baseline: +{improvement*100:.2f}%")
print(f"Target 80% reached: {all_results[0][1] >= 0.80}")
print(f"{'='*70}")

# Save results
results_df = pd.DataFrame([
    {'method': name, 'accuracy': acc}
    for name, acc in all_results
])
results_df.to_csv(PHASE2_ALTERNATIVE_RESULTS_CSV, index=False)
print(f"\nResults saved to: {PHASE2_ALTERNATIVE_RESULTS_CSV}")
