#!/usr/bin/env python3
"""
phase2_push_to_80.py - Aggressive push toward 80% accuracy
"""

import numpy as np
import pandas as pd
import json
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score
from sklearn.decomposition import PCA
from sklearn.base import clone
import xgboost as xgb
import lightgbm as lgb
import warnings

from bundle_paths import (
    CELL_LINEAGE_PATH,
    CELL_TYPE_CSV,
    PHASE2_PUSH_RESULTS_CSV,
    S3_CSV,
    ensure_results_dir,
)

ensure_results_dir()

warnings.filterwarnings('ignore')

print("="*70)
print("PUSHING TO 80%: Aggressive Optimization")
print("="*70)

# ─── Load Data (same as before) ──────────────────────────────────────────────
with CELL_LINEAGE_PATH.open('r', encoding='utf-8') as f:
    lineage_data = json.load(f)

def map_names(did):
    if   did == "P4a": return "Z3"
    elif did == "P4p": return "Z2"
    elif did == "P0a": return "AB"
    else: return did

terminal_nodes, intermediate_nodes = [], []
descendant_list_dict = defaultdict(list)

def dfs(node, parent, ancestors=[]):
    children = node.get("children", [])
    lookup_name = map_names(node["did"])
    if len(children) == 0:
        terminal_nodes.append(lookup_name)
        for a in ancestors:
            descendant_list_dict[a].append(lookup_name)
    else:
        intermediate_nodes.append(lookup_name)
        for c in children:
            dfs(c, node, ancestors + [lookup_name])

dfs(lineage_data, None)

cell_type_df = pd.read_csv(CELL_TYPE_CSV)
protein_exp = pd.read_csv(S3_CSV, index_col=0).T.fillna(0)
protein_exp = protein_exp[(protein_exp != 0).any(axis=1)]
protein_exp_zscore = protein_exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

cell_type_dict = {}
for node in terminal_nodes:
    df = cell_type_df[cell_type_df['wormweb.lineage'] == node]
    if len(df) == 0:
        continue
    types = df["wormweb.type"].dropna().unique()
    cell_type_dict[node] = types[0] if len(types) > 0 else "programmed_death"

cell_types = sorted(set(cell_type_dict.values()) - {"programmed_death"},
                    key=lambda x: -sum(1 for v in cell_type_dict.values() if v == x))
cell_types.append("programmed_death")
ct2i = {ct: i for i, ct in enumerate(cell_types)}

cell_type_one_hot = {}
for node, ct in cell_type_dict.items():
    oh = np.zeros(len(cell_types))
    oh[ct2i[ct]] = 1
    cell_type_one_hot[node] = oh

for node in intermediate_nodes:
    desc = [cell_type_one_hot[d] for d in descendant_list_dict[node]]
    if desc:
        s = np.sum(desc, axis=0)
        cell_type_one_hot[node] = s / s.sum()
    else:
        cell_type_one_hot[node] = np.zeros(len(cell_types))

X_all = protein_exp_zscore.values.T
y_all = np.array([cell_type_one_hot[map_names(n)] for n in protein_exp_zscore.columns])
terminal_mask = np.array([len(descendant_list_dict[map_names(n)]) == 0 for n in protein_exp_zscore.columns])

X_term = X_all[terminal_mask]
y_term_raw = y_all[terminal_mask].argmax(axis=1)

# Drop rare classes < 5 samples
counts = np.bincount(y_term_raw)
keep = np.where(counts >= 5)[0]
keep_mask = np.isin(y_term_raw, keep)
X_term = X_term[keep_mask]
y_term_raw = y_term_raw[keep_mask]
remap = {old: new for new, old in enumerate(sorted(np.unique(y_term_raw)))}
y_term = np.array([remap[c] for c in y_term_raw])
remap_types = [cell_types[c] for c in sorted(remap.keys())]
n_classes = len(np.unique(y_term))

print(f"\nData: {X_term.shape[0]} samples, {X_term.shape[1]} features, {n_classes} classes")
for i, ct in enumerate(remap_types):
    print(f"  {ct:25s}  n={np.sum(y_term == i)}")

# ─── Feature Engineering ──────────────────────────────────────────────────────
print("\n" + "="*70)
print("FEATURE ENGINEERING")
print("="*70)

# Strategy 1: PCA features as additional inputs
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_term)
pca = PCA(n_components=50)
X_pca = pca.fit_transform(X_scaled)
print(f"PCA: {X_pca.shape[1]} components explain {pca.explained_variance_ratio_.sum():.2%} variance")

# Strategy 2: Raw + PCA concatenation
X_augmented = np.hstack([X_term, X_pca])
print(f"Augmented features: {X_augmented.shape[1]} (raw {X_term.shape[1]} + PCA {X_pca.shape[1]})")

# Strategy 3: Top feature interactions
# Get top features from a quick random forest
rf_quick = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
rf_quick.fit(X_term, y_term)
top_feat_idx = np.argsort(rf_quick.feature_importances_)[::-1][:30]
print(f"Top 30 RF features identified")

# Add pairwise products of top features
X_interactions = []
for i in range(min(10, len(top_feat_idx))):
    for j in range(i+1, min(10, len(top_feat_idx))):
        X_interactions.append(X_term[:, top_feat_idx[i]] * X_term[:, top_feat_idx[j]])
X_interactions = np.array(X_interactions).T
X_full = np.hstack([X_term, X_pca, X_interactions])
print(f"Full feature set: {X_full.shape[1]} (raw + PCA + {X_interactions.shape[1]} interactions)")

# ─── Class weights ────────────────────────────────────────────────────────────
class_counts = np.bincount(y_term)
class_weight = {i: len(y_term) / (n_classes * count) for i, count in enumerate(class_counts) if count > 0}
print(f"\nClass weights: {[f'{w:.2f}' for w in class_weight.values()]}")

# ─── Massive XGBoost Grid Search ─────────────────────────────────────────────
print("\n" + "="*70)
print("MASSIVE XGBOOST GRID SEARCH")
print("="*70)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

feature_sets = {
    'raw_210': X_term,
    'pca_50': X_pca,
    'raw+pca': X_augmented,
    'raw+pca+interact': X_full,
}

xgb_configs = [
    # Shallow, heavily regularized
    {'n_estimators': 300, 'max_depth': 3, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.8, 'min_child_weight': 1, 'reg_alpha': 0, 'reg_lambda': 1},
    {'n_estimators': 500, 'max_depth': 3, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.7, 'min_child_weight': 1, 'reg_alpha': 0.1, 'reg_lambda': 2},
    # Medium depth
    {'n_estimators': 300, 'max_depth': 4, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.7, 'min_child_weight': 1, 'reg_alpha': 0, 'reg_lambda': 1},
    {'n_estimators': 500, 'max_depth': 4, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_weight': 2, 'reg_alpha': 0.1, 'reg_lambda': 1.5},
    {'n_estimators': 800, 'max_depth': 4, 'learning_rate': 0.03, 'subsample': 0.8, 'colsample_bytree': 0.5, 'min_child_weight': 3, 'reg_alpha': 0.3, 'reg_lambda': 2},
    # Standard depth
    {'n_estimators': 300, 'max_depth': 5, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.7, 'min_child_weight': 1, 'reg_alpha': 0, 'reg_lambda': 0.5},
    {'n_estimators': 500, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.85, 'colsample_bytree': 0.6, 'min_child_weight': 2, 'reg_alpha': 0.1, 'reg_lambda': 1},
    {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_weight': 3, 'reg_alpha': 0.1, 'reg_lambda': 1},
    # Deep
    {'n_estimators': 500, 'max_depth': 7, 'learning_rate': 0.05, 'subsample': 0.7, 'colsample_bytree': 0.5, 'min_child_weight': 5, 'reg_alpha': 0.5, 'reg_lambda': 3},
    {'n_estimators': 1000, 'max_depth': 6, 'learning_rate': 0.02, 'subsample': 0.7, 'colsample_bytree': 0.4, 'min_child_weight': 3, 'reg_alpha': 1, 'reg_lambda': 5},
    # Very regularized
    {'n_estimators': 500, 'max_depth': 4, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.5, 'min_child_weight': 5, 'reg_alpha': 1.0, 'reg_lambda': 5.0, 'gamma': 0.3},
    {'n_estimators': 300, 'max_depth': 3, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.9, 'min_child_weight': 1, 'reg_alpha': 0, 'reg_lambda': 0.5, 'gamma': 0},
]

all_results = []

for feat_name, X_feat in feature_sets.items():
    print(f"\nFeature set: {feat_name} ({X_feat.shape[1]} features)")

    for ci, params in enumerate(xgb_configs):
        # With and without class weights
        for use_weights in [False, True]:
            sw = None
            if use_weights:
                sw = np.array([class_weight[c] for c in y_term])

            fold_accs = []
            for train_idx, val_idx in skf.split(X_feat, y_term):
                model = xgb.XGBClassifier(
                    **params, random_state=42, tree_method='hist',
                    device='cuda', verbosity=0,
                    num_class=n_classes, objective='multi:softprob',
                )
                if sw is not None:
                    model.fit(X_feat[train_idx], y_term[train_idx], sample_weight=sw[train_idx])
                else:
                    model.fit(X_feat[train_idx], y_term[train_idx])
                preds = model.predict(X_feat[val_idx])
                fold_accs.append(accuracy_score(y_term[val_idx], preds))

            mean_acc = np.mean(fold_accs)
            all_results.append({
                'features': feat_name,
                'config_idx': ci,
                'weighted': use_weights,
                'mean_acc': mean_acc,
                'std_acc': np.std(fold_accs),
                'folds': fold_accs,
                **params,
            })

            if mean_acc >= 0.78:
                w_str = "+cw" if use_weights else ""
                print(f"  ★ Config {ci}{w_str}: {mean_acc:.4f} ± {np.std(fold_accs):.4f}")

# Sort and show top results
all_results.sort(key=lambda x: -x['mean_acc'])

print("\n" + "="*70)
print("TOP 15 CONFIGURATIONS")
print("="*70)

for i, r in enumerate(all_results[:15]):
    w_str = "+weighted" if r['weighted'] else ""
    print(f"  {i+1:2d}. {r['mean_acc']:.4f} ± {r['std_acc']:.4f}  "
          f"feat={r['features']:20s}  depth={r['max_depth']}  "
          f"n_est={r['n_estimators']}  lr={r['learning_rate']}  {w_str}")

# ─── LightGBM Grid Search ────────────────────────────────────────────────────
print("\n" + "="*70)
print("LIGHTGBM GRID SEARCH")
print("="*70)

lgb_configs = [
    {'n_estimators': 300, 'max_depth': 3, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.8, 'min_child_samples': 3, 'reg_alpha': 0, 'reg_lambda': 1, 'num_leaves': 15},
    {'n_estimators': 500, 'max_depth': 4, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_samples': 5, 'reg_alpha': 0.1, 'reg_lambda': 1, 'num_leaves': 31},
    {'n_estimators': 500, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.85, 'colsample_bytree': 0.6, 'min_child_samples': 3, 'reg_alpha': 0.1, 'reg_lambda': 1, 'num_leaves': 31},
    {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.6, 'min_child_samples': 5, 'reg_alpha': 0.1, 'reg_lambda': 1, 'num_leaves': 63},
    {'n_estimators': 800, 'max_depth': 4, 'learning_rate': 0.03, 'subsample': 0.8, 'colsample_bytree': 0.5, 'min_child_samples': 5, 'reg_alpha': 0.5, 'reg_lambda': 2, 'num_leaves': 15},
    {'n_estimators': 300, 'max_depth': -1, 'learning_rate': 0.1, 'subsample': 0.9, 'colsample_bytree': 0.8, 'min_child_samples': 10, 'reg_alpha': 0, 'reg_lambda': 1, 'num_leaves': 20},
]

lgb_results = []

for feat_name, X_feat in feature_sets.items():
    for ci, params in enumerate(lgb_configs):
        for use_weights in [False, True]:
            sw = np.array([class_weight[c] for c in y_term]) if use_weights else None

            fold_accs = []
            for train_idx, val_idx in skf.split(X_feat, y_term):
                model = lgb.LGBMClassifier(**params, random_state=42, verbose=-1, device='cpu')
                if sw is not None:
                    model.fit(X_feat[train_idx], y_term[train_idx], sample_weight=sw[train_idx])
                else:
                    model.fit(X_feat[train_idx], y_term[train_idx])
                preds = model.predict(X_feat[val_idx])
                fold_accs.append(accuracy_score(y_term[val_idx], preds))

            mean_acc = np.mean(fold_accs)
            lgb_results.append({
                'features': feat_name,
                'config_idx': ci,
                'weighted': use_weights,
                'mean_acc': mean_acc,
                'std_acc': np.std(fold_accs),
            })

            if mean_acc >= 0.78:
                w_str = "+cw" if use_weights else ""
                print(f"  ★ {feat_name} Config {ci}{w_str}: {mean_acc:.4f}")

lgb_results.sort(key=lambda x: -x['mean_acc'])
print(f"\nBest LightGBM: {lgb_results[0]['mean_acc']:.4f} ({lgb_results[0]['features']}, weighted={lgb_results[0]['weighted']})")

# ─── Multi-model Ensemble ────────────────────────────────────────────────────
print("\n" + "="*70)
print("MULTI-MODEL ENSEMBLE (Best configs)")
print("="*70)

# Get best feature set
best_feat = all_results[0]['features']
X_best = feature_sets[best_feat]
print(f"Using feature set: {best_feat}")

# Train top-5 XGBoost configs + top LightGBM + RF + ExtraTrees
ensemble_models = []

# Top 3 diverse XGBoost configs (different depths)
seen_depths = set()
for r in all_results:
    if r['features'] == best_feat and r['max_depth'] not in seen_depths:
        seen_depths.add(r['max_depth'])
        params = {k: v for k, v in r.items()
                  if k not in ('features', 'config_idx', 'weighted', 'mean_acc', 'std_acc', 'folds')}
        ensemble_models.append(('xgb_d' + str(r['max_depth']),
                                xgb.XGBClassifier(**params, random_state=42, tree_method='hist',
                                                  device='cuda', verbosity=0,
                                                  num_class=n_classes, objective='multi:softprob')))
        if len(seen_depths) >= 3:
            break

# Add LightGBM
ensemble_models.append(('lgb', lgb.LGBMClassifier(
    n_estimators=500, max_depth=5, learning_rate=0.05,
    subsample=0.85, colsample_bytree=0.6, min_child_samples=3,
    reg_alpha=0.1, reg_lambda=1, num_leaves=31,
    random_state=42, verbose=-1, device='cpu',
)))

# Add RF and ExtraTrees
ensemble_models.append(('rf', RandomForestClassifier(
    n_estimators=500, min_samples_leaf=2, max_features='sqrt', random_state=42, n_jobs=-1,
)))
ensemble_models.append(('et', ExtraTreesClassifier(
    n_estimators=500, min_samples_leaf=2, max_features='sqrt', random_state=42, n_jobs=-1,
)))

# Add L1 LogReg
ensemble_models.append(('lr', Pipeline([
    ('scaler', StandardScaler()),
    ('clf', LogisticRegression(C=1.0, penalty='l1', solver='saga', max_iter=5000, random_state=42)),
])))

print(f"Ensemble: {[name for name, _ in ensemble_models]}")

# Soft-voting ensemble via CV
from scipy.stats import mode

fold_accs_vote = []
fold_accs_avg = []

for train_idx, val_idx in skf.split(X_best, y_term):
    X_tr, X_val = X_best[train_idx], X_best[val_idx]
    y_tr, y_val = y_term[train_idx], y_term[val_idx]

    preds_list = []
    prob_list = []

    for name, model in ensemble_models:
        m = clone(model)
        m.fit(X_tr, y_tr)
        preds_list.append(m.predict(X_val))
        if hasattr(m, 'predict_proba'):
            prob_list.append(m.predict_proba(X_val))

    # Hard voting
    preds_stack = np.array(preds_list)
    vote_preds = mode(preds_stack, axis=0, keepdims=False).mode
    fold_accs_vote.append(accuracy_score(y_val, vote_preds))

    # Soft voting (probability averaging)
    if prob_list:
        avg_probs = np.mean(prob_list, axis=0)
        avg_preds = np.argmax(avg_probs, axis=1)
        fold_accs_avg.append(accuracy_score(y_val, avg_preds))

vote_acc = np.mean(fold_accs_vote)
avg_acc = np.mean(fold_accs_avg) if fold_accs_avg else 0

print(f"\nHard voting ensemble: {vote_acc:.4f} ± {np.std(fold_accs_vote):.4f}  folds: {[f'{a:.3f}' for a in fold_accs_vote]}")
print(f"Soft voting ensemble: {avg_acc:.4f} ± {np.std(fold_accs_avg):.4f}  folds: {[f'{a:.3f}' for a in fold_accs_avg]}")

# ─── Repeated CV for more stable estimates ────────────────────────────────────
print("\n" + "="*70)
print("REPEATED 5-FOLD CV (3 repeats) - BEST SINGLE MODEL")
print("="*70)

best_params = {k: v for k, v in all_results[0].items()
               if k not in ('features', 'config_idx', 'weighted', 'mean_acc', 'std_acc', 'folds')}
use_sw = all_results[0]['weighted']

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=42)
rep_accs = []

for train_idx, val_idx in rskf.split(X_best, y_term):
    model = xgb.XGBClassifier(**best_params, random_state=42, tree_method='hist',
                              device='cuda', verbosity=0,
                              num_class=n_classes, objective='multi:softprob')
    if use_sw:
        sw = np.array([class_weight[c] for c in y_term[train_idx]])
        model.fit(X_best[train_idx], y_term[train_idx], sample_weight=sw)
    else:
        model.fit(X_best[train_idx], y_term[train_idx])
    preds = model.predict(X_best[val_idx])
    rep_accs.append(accuracy_score(y_term[val_idx], preds))

print(f"Repeated CV: {np.mean(rep_accs):.4f} ± {np.std(rep_accs):.4f}")
print(f"Per-repeat means: {[f'{np.mean(rep_accs[i*5:(i+1)*5]):.4f}' for i in range(3)]}")

# ─── Final Summary ────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

best_single = all_results[0]['mean_acc']
best_ensemble = max(vote_acc, avg_acc)
best_overall = max(best_single, best_ensemble)

print(f"\nNN baseline (Phase 2 CV):     0.7574")
print(f"Best XGBoost single model:    {best_single:.4f}  (feat={all_results[0]['features']})")
print(f"Best LightGBM single model:   {lgb_results[0]['mean_acc']:.4f}")
print(f"Hard voting ensemble:         {vote_acc:.4f}")
print(f"Soft voting ensemble:         {avg_acc:.4f}")
print(f"Repeated CV (best model):     {np.mean(rep_accs):.4f} ± {np.std(rep_accs):.4f}")
print(f"\nBEST OVERALL:                 {best_overall:.4f}  (+{(best_overall-0.7574)*100:.2f}%)")
print(f"Target 80% reached:           {best_overall >= 0.80}")

# Save comprehensive results
results_df = pd.DataFrame(all_results[:20])
results_df.to_csv(PHASE2_PUSH_RESULTS_CSV, index=False)
print(f"\nSaved top 20 configs to: {PHASE2_PUSH_RESULTS_CSV}")
