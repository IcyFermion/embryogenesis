"""
phase2_improved_cv.py
─────────────────────────────────────────────────────────────────────────────
Expanded Phase 2 hyperparameter search with new architectures and training strategies.

This module extends the standard cross_validate_focused() with:
- Multiple model architectures (focused, resnet, wide, attention)
- Multiple optimizers (adam, adamw, sgd)
- Learning rate tuning
- Early stopping
- Label smoothing
- Gradient clipping
"""

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import StratifiedKFold
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
from itertools import combinations

from bundle_paths import ensure_repo_on_path

ensure_repo_on_path()

from feature_selection import (
    _train_sparse_gate,
    build_focused_classifier,
    eval_prediction_metrics,
    soft_cross_entropy,
    _batch_dist_corr_loss,
    _count_parameters,
    _jaccard,
    _resolve_device,
    _build_param_grid,
)


def train_focused_improved(
    X_sel,
    y,
    sample_weights,
    hidden_dims=(64, 32),
    dropout=0.2,
    n_epochs=600,
    batch_size=64,
    seed=42,
    device='auto',
    dist_lambda=0.1,
    model_type='focused',
    optimizer_type='adam',
    lr=1e-3,
    label_smoothing=0.0,
    early_stopping_patience=None,
    gradient_clip=None,
):
    """
    Enhanced Phase 2 training with multiple architectures and optimizers.
    Wrapper around direct training for consistency with CV loops.
    """
    ensure_repo_on_path()
    from feature_selection import train_focused

    return train_focused(
        X_sel, y, sample_weights,
        hidden_dims=hidden_dims, dropout=dropout,
        n_epochs=n_epochs, batch_size=batch_size,
        seed=seed, device=device, dist_lambda=dist_lambda,
        model_type=model_type, optimizer_type=optimizer_type,
        lr=lr, label_smoothing=label_smoothing,
        early_stopping_patience=early_stopping_patience,
        gradient_clip=gradient_clip,
    )


def cross_validate_focused_improved(
    X,
    y,
    sample_weights,
    terminal_mask,
    selector_config,
    param_grid,
    n_splits=5,
    selector_epochs=250,
    focused_epochs=400,
    selector_batch_size=128,
    focused_batch_size=64,
    seed=42,
    device='auto',
    score_weights=None,
):
    """
    Expanded Phase 2 CV with support for multiple architectures and training strategies.

    Parameters
    ----------
    X              : ndarray (n_samples, n_features)
    y              : ndarray (n_samples, n_classes)  – soft labels
    sample_weights : ndarray (n_samples,)
    terminal_mask  : ndarray (n_samples,) bool – hard-label samples
    selector_config : dict – best Phase 1 configuration
    param_grid     : dict of lists
        May include:
        - n_select, hidden_dims, dropout, dist_lambda (original)
        - model_type: ['focused', 'resnet', 'wide', 'attention']
        - optimizer_type: ['adam', 'adamw', 'sgd']
        - lr: [1e-3, 2e-3, ...]
        - label_smoothing: [0.0, 0.1, ...]
        - early_stopping_patience: [None, 50, 100]
        - gradient_clip: [None, 1.0, 5.0]

    Returns
    -------
    results     : list[dict]  – sorted by composite score
    best_config : dict       – best configuration
    """
    device = _resolve_device(device)
    weights = {
        'val_acc': 0.55,
        'soft_target_probability': 0.25,
        'feature_stability': 0.10,
        'overfit_gap': 0.15,
        'param_count': 0.10,
    }
    if score_weights is not None:
        weights.update(score_weights)

    configs = _build_param_grid(param_grid)

    strat_labels = y.argmax(axis=1)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = list(skf.split(X, strat_labels))

    selector_hidden_dims = selector_config.get('hidden_dims', (128, 64))
    selector_l1_lambda = selector_config.get('l1_lambda', 0.004)
    selector_dropout = selector_config.get('dropout', 0.3)

    results = []

    for cfg_idx, cfg in enumerate(tqdm(configs, desc="Phase-2 CV configs")):
        # Extract config parameters
        n_select = cfg.get('n_select', 15)
        focused_hidden_dims = cfg.get('hidden_dims', (64, 32))
        focused_dropout = cfg.get('dropout', 0.2)
        dist_lambda = cfg.get('dist_lambda', 0.1)
        model_type = cfg.get('model_type', 'focused')
        optimizer_type = cfg.get('optimizer_type', 'adam')
        lr = cfg.get('lr', 1e-3)
        label_smoothing = cfg.get('label_smoothing', 0.0)
        early_stopping_patience = cfg.get('early_stopping_patience', None)
        gradient_clip = cfg.get('gradient_clip', None)

        fold_results = []
        for fold_idx, (train_idx, val_idx) in enumerate(folds):
            X_tr, X_val = X[train_idx], X[val_idx]
            y_tr, y_val = y[train_idx], y[val_idx]
            w_tr, w_val = sample_weights[train_idx], sample_weights[val_idx]
            term_tr = terminal_mask[train_idx]
            term_val = terminal_mask[val_idx]

            # Phase 1: select features
            selector_model, _ = _train_sparse_gate(
                X_tr,
                y_tr,
                w_tr,
                selector_hidden_dims,
                selector_l1_lambda,
                selector_dropout,
                selector_epochs,
                selector_batch_size,
                seed + fold_idx,
                device,
            )
            gate_vals = selector_model.feature_importance()
            top_k = np.argsort(gate_vals)[::-1][:n_select].tolist()

            X_tr_sel = X_tr[:, top_k].astype(np.float32)
            X_val_sel = X_val[:, top_k].astype(np.float32)

            # Phase 2: train focused model with new parameters
            model, _, _ = train_focused_improved(
                X_tr_sel,
                y_tr,
                w_tr,
                hidden_dims=focused_hidden_dims,
                dropout=focused_dropout,
                n_epochs=focused_epochs,
                batch_size=focused_batch_size,
                seed=seed + fold_idx,
                device=str(device),
                dist_lambda=dist_lambda,
                model_type=model_type,
                optimizer_type=optimizer_type,
                lr=lr,
                label_smoothing=label_smoothing,
                early_stopping_patience=early_stopping_patience,
                gradient_clip=gradient_clip,
            )

            X_tr_t = torch.FloatTensor(X_tr_sel)
            X_val_t = torch.FloatTensor(X_val_sel)
            y_tr_t = torch.FloatTensor(y_tr)
            y_val_t = torch.FloatTensor(y_val)

            train_metrics = eval_prediction_metrics(
                model,
                X_tr_t,
                y_tr_t,
                term_tr,
                device,
            )
            val_hard_metrics = eval_prediction_metrics(
                model,
                X_val_t,
                y_val_t,
                term_val,
                device,
            )
            val_soft_metrics = eval_prediction_metrics(
                model,
                X_val_t,
                y_val_t,
                None,
                device,
                sample_weights=w_val,
            )

            fold_results.append({
                'fold': fold_idx,
                'train_acc': train_metrics['argmax_accuracy'],
                'val_acc': val_hard_metrics['argmax_accuracy'],
                'val_expected_target_probability': val_soft_metrics['expected_target_probability'],
                'val_soft_cross_entropy': val_soft_metrics['soft_cross_entropy'],
                'top_k_indices': top_k,
                'n_parameters': _count_parameters(model),
            })

        # Compute stability metrics
        all_top_k = [fr['top_k_indices'] for fr in fold_results]
        pair_jaccards = [
            _jaccard(all_top_k[i], all_top_k[j])
            for i, j in combinations(range(n_splits), 2)
        ]

        n_features = X.shape[1]
        freq = np.zeros(n_features, dtype=int)
        for tk in all_top_k:
            freq[np.array(tk)] += 1
        consensus = np.where(freq >= int(np.ceil(n_splits * 0.8)))[0].tolist()

        val_accs = [fr['val_acc'] for fr in fold_results if not np.isnan(fr['val_acc'])]
        train_accs = [fr['train_acc'] for fr in fold_results]
        overfit_gaps = [tr - va for tr, va in zip(train_accs, val_accs)]
        expected_target_probs = [fr['val_expected_target_probability'] for fr in fold_results]
        val_soft_cross_entropies = [fr['val_soft_cross_entropy'] for fr in fold_results]
        n_parameters = fold_results[0]['n_parameters'] if fold_results else 0

        stability_summary = {
            'mean_jaccard': float(np.mean(pair_jaccards)) if pair_jaccards else float('nan'),
            'consensus_ratio': float(min(len(consensus), n_select) / n_select),
        }
        selection_stability = (
            0.7 * stability_summary['mean_jaccard']
            + 0.3 * stability_summary['consensus_ratio']
        )

        summary = {
            'mean_val_acc': float(np.mean(val_accs)) if val_accs else float('nan'),
            'std_val_acc': float(np.std(val_accs)) if val_accs else float('nan'),
            'mean_train_acc': float(np.mean(train_accs)) if train_accs else float('nan'),
            'mean_overfit_gap': float(np.mean(overfit_gaps)) if overfit_gaps else float('nan'),
            'mean_expected_target_probability': float(np.mean(expected_target_probs)),
            'mean_soft_cross_entropy': float(np.mean(val_soft_cross_entropies)),
            'mean_jaccard': stability_summary['mean_jaccard'],
            'consensus_ratio': stability_summary['consensus_ratio'],
            'selection_stability_score': selection_stability,
            'n_parameters': n_parameters,
        }

        results.append({
            'config': cfg,
            'fold_results': fold_results,
            'stability': {
                'pair_jaccards': pair_jaccards,
                'feature_frequency': freq,
                'consensus_features': consensus,
            },
            'summary': summary,
            'score': 0.0,
        })

    # Compute composite scores
    if results:
        param_counts = np.array([r['summary']['n_parameters'] for r in results], dtype=np.float64)
        if np.allclose(param_counts.max(), param_counts.min()):
            normalized_param_counts = np.zeros_like(param_counts)
        else:
            normalized_param_counts = (param_counts - param_counts.min()) / (param_counts.max() - param_counts.min())

        for result, param_penalty in zip(results, normalized_param_counts):
            summary = result['summary']
            mean_val_acc = summary['mean_val_acc']
            score = 0.0
            if not np.isnan(mean_val_acc):
                score = (
                    weights['val_acc'] * mean_val_acc
                    + weights['soft_target_probability'] * summary['mean_expected_target_probability']
                    + weights['feature_stability'] * summary['selection_stability_score']
                    - weights['overfit_gap'] * max(summary['mean_overfit_gap'], 0.0)
                    - weights['param_count'] * param_penalty
                )
            summary['normalized_param_count'] = float(param_penalty)
            result['score'] = score

    results.sort(key=lambda r: r['score'], reverse=True)
    best_config = results[0]['config']
    return results, best_config
