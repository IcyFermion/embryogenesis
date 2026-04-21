#!/usr/bin/env python3
"""
phase2_full_optimization.py
Complete Phase 2 optimization pipeline with all stages.
Runs in dev conda environment.
"""

import os
import json
import sys
import numpy as np
import pandas as pd
import torch
import pickle
from collections import defaultdict
from tqdm import tqdm

from bundle_paths import (
    CELL_LINEAGE_PATH,
    CELL_TYPE_CSV,
    S3_CSV,
    STAGE1_RESULTS_CSV,
    STAGE1_RESULTS_PKL,
    STAGE2_RESULTS_CSV,
    STAGE2_RESULTS_PKL,
    ensure_repo_on_path,
    ensure_results_dir,
)

ensure_repo_on_path()
ensure_results_dir()

from feature_selection import (
    cross_validate_features,
    entropy_based_weights,
    eval_prediction_metrics,
    train_one_pass,
    train_focused,
)
from phase2_improved_cv import cross_validate_focused_improved

# Device setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Device: {device}\n')

# Load data
print("Loading data...")
import json
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

# Load protein expression
cell_type_df = pd.read_csv(CELL_TYPE_CSV)
protein_exp = pd.read_csv(S3_CSV, index_col=0).T
protein_exp = protein_exp.fillna(0)
protein_exp = protein_exp[(protein_exp != 0).any(axis=1)]
protein_exp_zscore = protein_exp.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

print(f'Protein expression shape: {protein_exp_zscore.shape}')
print(f'Terminal nodes: {len(terminal_nodes)}, Intermediate: {len(intermediate_nodes)}')

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
    cur_type = cur_lineage_types[0]
    cell_type_dict[node] = cur_type

cell_types = list(set(cell_type_dict.values()))
cell_types = [ct for ct in cell_types if ct != "programmed_death"]
cell_types = sorted(cell_types, key=lambda x: sum([1 for v in cell_type_dict.values() if v == x]), reverse=True)
cell_types.append("programmed_death")
cell_type_to_int = {ct: i for i, ct in enumerate(cell_types)}

# Build soft labels
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
    normalized_types = summed_types / np.sum(summed_types)
    cell_type_one_hot[node] = normalized_types

X = protein_exp_zscore.values.T
y = np.array([cell_type_one_hot[map_names(node)] for node in protein_exp_zscore.columns])
n_descendants = np.array([len(descendant_list_dict[map_names(node)]) for node in protein_exp_zscore.columns])

sample_weights_v2 = entropy_based_weights(y, alpha=3.0)
terminal_mask = np.array([len(descendant_list_dict[map_names(node)]) == 0 for node in protein_exp_zscore.columns])

print(f'Data shape: X={X.shape}, y={y.shape}')
print(f'Terminal (hard): {terminal_mask.sum()}, Intermediate (soft): {(~terminal_mask).sum()}')
print(f'Number of classes: {len(cell_types)}')

# Best Phase 1 config
best_selector_config = {
    'l1_lambda': 0.004,
    'hidden_dims': (128, 64),
    'dropout': 0.2,
}

print('\n' + '='*70)
print('STAGE 1: Architecture Comparison')
print('='*70)

STAGE1_PARAM_GRID = {
    'n_select': [25],
    'hidden_dims': [(128, 64)],
    'dropout': [0.2],
    'dist_lambda': [0.0],
    'model_type': ['focused', 'resnet', 'wide', 'attention'],
    'optimizer_type': ['adam'],
    'lr': [1e-3],
}

print(f'\nTesting {len(STAGE1_PARAM_GRID["model_type"])} architectures...')

try:
    stage1_results, stage1_best = cross_validate_focused_improved(
        X, y, sample_weights_v2, terminal_mask,
        selector_config=best_selector_config,
        param_grid=STAGE1_PARAM_GRID,
        n_splits=5,
        selector_epochs=250,
        focused_epochs=400,
        seed=42,
        device=str(device),
    )

    stage1_df = pd.DataFrame([
        {
            'model_type': r['config'].get('model_type', 'focused'),
            'mean_val_acc': r['summary']['mean_val_acc'],
            'mean_overfit_gap': r['summary']['mean_overfit_gap'],
            'soft_target_prob': r['summary']['mean_expected_target_probability'],
            'score': r['score'],
        }
        for r in stage1_results
    ])

    print('\nStage 1 Results:')
    print(stage1_df.to_string(index=False))
    print(f'\nBest architecture: {stage1_best["model_type"]}')
    print(f'  Validation accuracy: {stage1_results[0]["summary"]["mean_val_acc"]:.4f}')
    print(f'  Improvement from baseline (0.7574): +{(stage1_results[0]["summary"]["mean_val_acc"] - 0.7574)*100:.2f}%')

    # Save results
    with STAGE1_RESULTS_PKL.open('wb') as f:
        pickle.dump(stage1_results, f)
    stage1_df.to_csv(STAGE1_RESULTS_CSV, index=False)

except Exception as e:
    print(f'ERROR in Stage 1: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)

print('\n' + '='*70)
print('STAGE 2: Expanded Hyperparameter Grid Search')
print('='*70)

best_model_type = stage1_results[0]['config']['model_type']
print(f'Using best model type from Stage 1: {best_model_type}')

STAGE2_PARAM_GRID = {
    'n_select': [20, 25, 30],
    'hidden_dims': [(64, 64), (128, 64), (128, 128, 64), (256, 128, 64)],
    'dropout': [0.0, 0.1, 0.2],
    'dist_lambda': [0.0, 0.05, 0.1],
    'model_type': [best_model_type],
    'optimizer_type': ['adam', 'adamw'],
    'lr': [5e-4, 1e-3, 2e-3],
    'label_smoothing': [0.0, 0.1],
}

n_configs_stage2 = np.prod([len(v) for v in STAGE2_PARAM_GRID.values()])
print(f'\nGrid size: {n_configs_stage2} configurations')
print(f'Running 5-fold CV (this will take 30-60 minutes)...\n')

try:
    stage2_results, stage2_best = cross_validate_focused_improved(
        X, y, sample_weights_v2, terminal_mask,
        selector_config=best_selector_config,
        param_grid=STAGE2_PARAM_GRID,
        n_splits=5,
        selector_epochs=250,
        focused_epochs=400,
        seed=42,
        device=str(device),
    )

    print(f'\n=== Stage 2 Complete ===')
    best_stage2_acc = stage2_results[0]['summary']['mean_val_acc']
    print(f'Best accuracy: {best_stage2_acc:.4f}')
    print(f'Improvement from baseline (0.7574): +{(best_stage2_acc - 0.7574)*100:.2f}%')

    print(f'\nBest config:')
    for k, v in stage2_best.items():
        print(f'  {k}: {v}')

    # Show top 10
    stage2_df = pd.DataFrame([
        {
            'n_select': r['config']['n_select'],
            'hidden_dims': str(r['config']['hidden_dims']),
            'dropout': r['config']['dropout'],
            'lr': r['config']['lr'],
            'dist_lambda': r['config']['dist_lambda'],
            'label_smoothing': r['config'].get('label_smoothing', 0.0),
            'optimizer': r['config'].get('optimizer_type', 'adam'),
            'val_acc': r['summary']['mean_val_acc'],
            'overfit_gap': r['summary']['mean_overfit_gap'],
            'soft_target': r['summary']['mean_expected_target_probability'],
            'score': r['score'],
        }
        for r in stage2_results[:10]
    ])

    print('\nTop 10 Stage 2 Configurations:')
    print(stage2_df.to_string(index=False))

    # Save results
    with STAGE2_RESULTS_PKL.open('wb') as f:
        pickle.dump(stage2_results, f)
    stage2_df.to_csv(STAGE2_RESULTS_CSV, index=False)

except Exception as e:
    print(f'ERROR in Stage 2: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Summary
print('\n' + '='*70)
print('OPTIMIZATION SUMMARY')
print('='*70)
print(f'\nBaseline accuracy: 0.7574')
print(f'Stage 1 best: {stage1_results[0]["summary"]["mean_val_acc"]:.4f}')
print(f'Stage 2 best: {best_stage2_acc:.4f}')
print(f'\nTotal improvement: +{(best_stage2_acc - 0.7574)*100:.2f}%')
print(f'Target 80% reached: {best_stage2_acc >= 0.80}')

print('\n✓ Optimization complete! Results saved:')
print(f'  - {STAGE1_RESULTS_PKL}')
print(f'  - {STAGE1_RESULTS_CSV}')
print(f'  - {STAGE2_RESULTS_PKL}')
print(f'  - {STAGE2_RESULTS_CSV}')
