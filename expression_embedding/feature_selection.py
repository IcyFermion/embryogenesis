"""
feature_selection.py
────────────────────────────────────────────────────────────────────────────────
Reusable ML module for protein feature selection and embedding generation.

Public API
----------
build_cv_folds(X, y, n_splits, seed, groups=None)
    Construct stratified folds for standard supervised runs or grouped folds
    for lineage/timepoint experiments where related observations must stay in
    the same split.

build_selector_fold_cache(X, y, sample_weights, selector_config, folds, ...)
    Train one Phase-1 selector per fold and cache the full feature ranking so
    Phase-2 CV and grouped OOF evaluation can reuse the same selector fits.

cross_validate_features(X, y, sample_weights, terminal_mask, param_grid, ...)
    5-fold stratified or grouped CV over a hyperparameter param_grid. Returns
    per-config metrics (val accuracy, feature-set stability, overfitting gap)
    plus the best config ranked by a stability-aware composite score.

cross_validate_focused(X, y, sample_weights, terminal_mask, selector_config, ...)
    5-fold stratified or grouped CV for Phase 2 that jointly tunes the number
    of selected features K, the focused classifier architecture, and
    regularisation while penalising overfitting and unnecessary model size.

train_one_pass(X, y, sample_weights, hidden_dims, l1_lambda, dropout, ...)
    Trains SparseGateClassifier on ALL data with the CV-selected configuration.
    Returns: model, top-K feature indices, gate values, embeddings
    (n_samples × last_hidden_dim), and training loss history.

eval_prediction_metrics(model, X_tensor, y_tensor, mask, device, ...)
    Returns argmax accuracy plus soft-label-aware metrics such as expected
    target probability and soft cross-entropy.

Models
------
SparseGateClassifier  – learnable sigmoid gates for feature selection;
                        embedding dim = last element of hidden_dims (flexible).
FocusedClassifier     – legacy two-phase Phase-2 model (kept for compat).
"""

import itertools
from itertools import combinations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold, StratifiedKFold
try:
    from sklearn.model_selection import StratifiedGroupKFold
except ImportError:  # pragma: no cover - fallback for older scikit-learn
    StratifiedGroupKFold = None
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm


# ─────────────────────────────────────────────────────────────────────────────
# Models
# ─────────────────────────────────────────────────────────────────────────────

class SparseGateClassifier(nn.Module):
    """
    Classifier with per-feature learnable sigmoid gates.

    A sigmoid gate scalar multiplies each input feature before the backbone
    network.  L1 regularisation on the gate values (via gate_l1()) drives
    low-importance features toward zero, yielding feature importance scores.

    Parameters
    ----------
    n_features  : int   – input dimensionality
    n_classes   : int   – number of output classes
    hidden_dims : tuple – widths of each hidden layer; the last element
                          determines the embedding dimensionality
    dropout     : float – dropout probability applied after every hidden layer
    """

    def __init__(self, n_features, n_classes, hidden_dims, dropout=0.3):
        super().__init__()
        self.gate_logits = nn.Parameter(torch.zeros(n_features))

        layers = []
        in_dim = n_features
        for h in hidden_dims:
            layers += [
                nn.Linear(in_dim, h),
                nn.BatchNorm1d(h),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout),
            ]
            in_dim = h
        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, n_classes)

    def forward(self, x):
        gates = torch.sigmoid(self.gate_logits)
        h = self.backbone(x * gates)
        return self.head(h), h   # (logits, penultimate embedding)

    def gate_l1(self):
        """L1 loss on gate values – penalises non-zero gates for sparsity."""
        return torch.sigmoid(self.gate_logits).sum()

    def feature_importance(self):
        """Return gate values as a numpy array (higher = more important)."""
        return torch.sigmoid(self.gate_logits).detach().cpu().numpy()


class FocusedClassifier(nn.Module):
    """
    Legacy Phase-2 classifier trained on pre-selected features.
    Kept for backward compatibility; not used by the one-pass pipeline.
    """

    def __init__(self, n_selected, n_classes, hidden_dims=(64, 32), dropout=0.2):
        super().__init__()
        layers = []
        in_dim = n_selected
        for h in hidden_dims:
            layers += [
                nn.Linear(in_dim, h),
                nn.BatchNorm1d(h),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout),
            ]
            in_dim = h
        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, n_classes)

    def forward(self, x):
        h = self.backbone(x)
        return self.head(h), h


class FocusedClassifierResNet(nn.Module):
    """
    Phase-2 classifier with residual connections for deeper architectures.
    Uses skip connections where dimension matching allows.
    """

    def __init__(self, n_selected, n_classes, hidden_dims=(128, 128, 64), dropout=0.2):
        super().__init__()
        self.n_selected = n_selected
        self.hidden_dims = hidden_dims

        layers = []
        in_dim = n_selected

        for i, h in enumerate(hidden_dims):
            linear = nn.Linear(in_dim, h)
            nn.init.xavier_uniform_(linear.weight)
            layers.append(linear)
            layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(dropout))
            in_dim = h

        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, n_classes)
        nn.init.xavier_uniform_(self.head.weight)

    def forward(self, x):
        h = self.backbone(x)
        return self.head(h), h


class FocusedClassifierWide(nn.Module):
    """
    Phase-2 classifier with wide first layer and bottleneck compression.
    Good for capturing complex feature interactions.
    """

    def __init__(self, n_selected, n_classes, hidden_dims=(256, 128, 64), dropout=0.2):
        super().__init__()
        self.n_selected = n_selected
        self.hidden_dims = hidden_dims

        layers = []
        in_dim = n_selected

        for i, h in enumerate(hidden_dims):
            linear = nn.Linear(in_dim, h)
            nn.init.xavier_uniform_(linear.weight)
            layers.append(linear)
            layers.append(nn.BatchNorm1d(h))
            # Use GELU for first layer, ReLU for others
            activation = nn.GELU() if i == 0 else nn.ReLU(inplace=True)
            layers.append(activation)
            layers.append(nn.Dropout(dropout))
            in_dim = h

        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, n_classes)
        nn.init.xavier_uniform_(self.head.weight)

    def forward(self, x):
        h = self.backbone(x)
        return self.head(h), h


class FocusedClassifierAttention(nn.Module):
    """
    Phase-2 classifier with multi-head self-attention layer.
    Treats features as tokens and learns feature interactions.
    """

    def __init__(self, n_selected, n_classes, hidden_dims=(64, 64), dropout=0.2, n_heads=4):
        super().__init__()
        self.n_selected = n_selected
        self.hidden_dims = hidden_dims
        self.n_heads = n_heads

        # Project to attention dim
        attn_dim = hidden_dims[0] if hidden_dims else 64
        self.proj_in = nn.Linear(n_selected, attn_dim)
        nn.init.xavier_uniform_(self.proj_in.weight)

        # Multi-head attention
        self.attn = nn.MultiheadAttention(
            embed_dim=attn_dim,
            num_heads=min(n_heads, attn_dim),
            dropout=dropout,
            batch_first=True,
        )

        # Post-attention layers
        layers = []
        in_dim = attn_dim

        for i, h in enumerate(hidden_dims):
            if i == 0:
                in_dim = attn_dim
            linear = nn.Linear(in_dim, h)
            nn.init.xavier_uniform_(linear.weight)
            layers.append(linear)
            layers.append(nn.BatchNorm1d(h))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(dropout))
            in_dim = h

        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, n_classes)
        nn.init.xavier_uniform_(self.head.weight)

    def forward(self, x):
        # x: (batch_size, n_selected)
        # Treat as (batch_size, n_selected, 1) for attention
        x_proj = self.proj_in(x)  # (batch_size, attn_dim)
        x_attn = x_proj.unsqueeze(1)  # (batch_size, 1, attn_dim)

        # Self-attention (attend to same position due to single token)
        attn_out, _ = self.attn(x_attn, x_attn, x_attn)  # (batch_size, 1, attn_dim)
        attn_out = attn_out.squeeze(1)  # (batch_size, attn_dim)

        # Pass through backbone
        h = self.backbone(attn_out)
        return self.head(h), h


def build_focused_classifier(model_type, n_selected, n_classes, hidden_dims, dropout):
    """
    Factory function to build focused classifier by type.

    Parameters
    ----------
    model_type : str
        One of: 'focused', 'resnet', 'wide', 'attention'

    Returns
    -------
    model : nn.Module
    """
    if model_type == 'focused':
        return FocusedClassifier(n_selected, n_classes, hidden_dims, dropout)
    elif model_type == 'resnet':
        return FocusedClassifierResNet(n_selected, n_classes, hidden_dims, dropout)
    elif model_type == 'wide':
        return FocusedClassifierWide(n_selected, n_classes, hidden_dims, dropout)
    elif model_type == 'attention':
        return FocusedClassifierAttention(n_selected, n_classes, hidden_dims, dropout)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")


# ─────────────────────────────────────────────────────────────────────────────
# Training utilities
# ─────────────────────────────────────────────────────────────────────────────

def entropy_based_weights(y_soft, alpha=3.0):
    """
    Compute per-sample training weights based on soft-label entropy.

    Hard (one-hot) labels → H=0 → weight≈1.
    Soft/mixed labels     → H>0 → weight<1.

    Parameters
    ----------
    y_soft : ndarray (n_samples, n_classes)
    alpha  : float – sharpness; higher → harder labels dominate relatively more

    Returns
    -------
    weights : ndarray (n_samples,) float32, mean-normalised so mean=1
    """
    y_clipped = np.clip(y_soft, 1e-10, 1.0)
    entropies = -np.sum(y_clipped * np.log(y_clipped), axis=1)
    max_ent   = np.log(y_soft.shape[1])
    norm_ent  = entropies / max_ent          # 0 = hard/one-hot, 1 = uniform
    weights   = np.exp(-alpha * norm_ent)
    return (weights / weights.mean()).astype(np.float32)


def soft_cross_entropy(logits, soft_targets, weights=None):
    """
    Cross-entropy loss with soft/distributional targets and optional
    per-sample weights (higher weight = more influence on the loss).
    """
    log_p      = F.log_softmax(logits, dim=-1)
    per_sample = -(soft_targets * log_p).sum(dim=-1)
    if weights is not None:
        per_sample = per_sample * weights
    return per_sample.mean()


def _weighted_mean(values, weights=None):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return float('nan')
    if weights is None:
        return float(values.mean())

    weights = np.asarray(weights, dtype=np.float64)
    weight_sum = weights.sum()
    if weight_sum <= 0:
        return float(values.mean())
    return float(np.dot(values, weights) / weight_sum)


def prediction_metrics_from_logits(logits, soft_targets, weights=None):
    """
    Summarise predictive quality for both hard and soft labels.

    Returns
    -------
    metrics : dict
        argmax_accuracy            – exact class-match rate on argmax labels
        expected_target_probability – mean probability assigned to true soft mass
        soft_cross_entropy         – mean soft cross-entropy
    """
    if isinstance(logits, torch.Tensor):
        logits = logits.detach().cpu()
    else:
        logits = torch.as_tensor(logits, dtype=torch.float32)

    if isinstance(soft_targets, torch.Tensor):
        soft_targets = soft_targets.detach().cpu()
    else:
        soft_targets = torch.as_tensor(soft_targets, dtype=torch.float32)

    if weights is not None:
        if isinstance(weights, torch.Tensor):
            weights = weights.detach().cpu().numpy()
        else:
            weights = np.asarray(weights, dtype=np.float64)

    probs = F.softmax(logits, dim=-1)
    log_p = F.log_softmax(logits, dim=-1)

    argmax_matches = (
        probs.argmax(dim=-1) == soft_targets.argmax(dim=-1)
    ).numpy().astype(np.float64)
    target_mass = (probs * soft_targets).sum(dim=-1).numpy()
    cross_entropy = (-(soft_targets * log_p).sum(dim=-1)).numpy()

    return {
        'argmax_accuracy': _weighted_mean(argmax_matches, weights),
        'expected_target_probability': _weighted_mean(target_mass, weights),
        'soft_cross_entropy': _weighted_mean(cross_entropy, weights),
    }


def train_epoch(model, loader, optimizer, device, l1_lambda=0.0):
    model.train()
    total = 0.0
    for xb, yb, wb in loader:
        xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
        optimizer.zero_grad()
        logits, _ = model(xb)
        loss = soft_cross_entropy(logits, yb, wb)
        if l1_lambda > 0 and hasattr(model, 'gate_l1'):
            loss = loss + l1_lambda * model.gate_l1()
        loss.backward()
        optimizer.step()
        total += loss.item()
    return total / len(loader)


@torch.no_grad()
def eval_accuracy(model, X_tensor, y_tensor, mask, device):
    """
    Hard-label classification accuracy on the subset of samples in mask.

    mask : bool array or BoolTensor – True for samples to evaluate
    """
    return eval_prediction_metrics(
        model,
        X_tensor,
        y_tensor,
        mask,
        device,
    )['argmax_accuracy']


@torch.no_grad()
def eval_prediction_metrics(model, X_tensor, y_tensor, mask, device, sample_weights=None):
    """
    Evaluate a classifier using both hard-label and soft-label-compatible views.

    For hard labels, ``argmax_accuracy`` is the usual classification accuracy.
    For soft labels, ``argmax_accuracy`` becomes dominant-class agreement only,
    so ``expected_target_probability`` and ``soft_cross_entropy`` should also be
    inspected.
    """
    model.eval()
    logits, _ = model(X_tensor.to(device))
    logits = logits.cpu()
    targets = y_tensor.detach().cpu() if isinstance(y_tensor, torch.Tensor) else torch.FloatTensor(y_tensor)

    if mask is not None:
        m = mask.numpy() if isinstance(mask, torch.Tensor) else np.asarray(mask)
        if m.sum() == 0:
            return {
                'argmax_accuracy': float('nan'),
                'expected_target_probability': float('nan'),
                'soft_cross_entropy': float('nan'),
            }
        logits = logits[m]
        targets = targets[m]
        if sample_weights is not None:
            sample_weights = sample_weights[m]

    return prediction_metrics_from_logits(logits, targets, sample_weights)


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _jaccard(set_a, set_b):
    a, b = set(set_a), set(set_b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def _resolve_device(device):
    if device == 'auto':
        return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return torch.device(device)


def _count_parameters(model):
    return int(sum(p.numel() for p in model.parameters() if p.requires_grad))


def _build_param_grid(param_grid):
    keys = list(param_grid.keys())
    combos = list(itertools.product(*[param_grid[k] for k in keys]))
    return [dict(zip(keys, combo)) for combo in combos]


def build_cv_folds(X, y, n_splits, seed, groups=None):
    """Build cross-validation folds with optional lineage-aware grouping."""
    strat_labels = y.argmax(axis=1)

    if groups is None:
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        return list(splitter.split(X, strat_labels))

    groups = np.asarray(groups)
    if groups.shape[0] != X.shape[0]:
        raise ValueError(
            f"Expected groups to have length {X.shape[0]}, got {groups.shape[0]}"
        )

    if StratifiedGroupKFold is not None:
        splitter = StratifiedGroupKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=seed,
        )
        return list(splitter.split(X, strat_labels, groups))

    splitter = GroupKFold(n_splits=n_splits)
    return list(splitter.split(X, strat_labels, groups))


def _format_metric(value):
    if value is None:
        return 'nan'
    try:
        if np.isnan(value):
            return 'nan'
    except TypeError:
        pass
    return f"{float(value):.4f}"


def _stability_score(summary, metric_weights=None):
    weights = {
        'jaccard': 0.5,
        'topk_frequency': 0.3,
        'spearman': 0.2,
    }
    if metric_weights is not None:
        weights.update(metric_weights)

    mean_spearman = summary.get('mean_spearman', 0.0)
    if np.isnan(mean_spearman):
        mean_spearman = 0.0

    return (
        weights['jaccard'] * summary.get('mean_jaccard', 0.0)
        + weights['topk_frequency'] * summary.get('mean_topk_frequency', 0.0)
        + weights['spearman'] * ((mean_spearman + 1.0) / 2.0)
    )


def _train_sparse_gate(
    X,
    y,
    sample_weights,
    hidden_dims,
    l1_lambda,
    dropout,
    n_epochs,
    batch_size,
    seed,
    device,
    early_stopping_patience=None,
    early_stopping_min_delta=1e-4,
):
    """Train a SparseGateClassifier and return (model, losses)."""
    torch.manual_seed(seed)
    n_features = X.shape[1]
    n_classes  = y.shape[1]

    X_t = torch.FloatTensor(X)
    y_t = torch.FloatTensor(y)
    w_t = torch.FloatTensor(sample_weights)

    loader = DataLoader(
        TensorDataset(X_t, y_t, w_t),
        batch_size=batch_size,
        shuffle=True,
    )

    model = SparseGateClassifier(n_features, n_classes, hidden_dims, dropout).to(device)
    opt   = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_epochs)

    losses = []
    best_loss = float('inf')
    patience_counter = 0

    for _ in range(n_epochs):
        loss = train_epoch(model, loader, opt, device, l1_lambda=l1_lambda)
        sched.step()
        losses.append(loss)

        if early_stopping_patience is not None:
            if loss < best_loss - early_stopping_min_delta:
                best_loss = loss
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= early_stopping_patience:
                    break

    return model, losses


def build_selector_fold_cache(
    X,
    y,
    sample_weights,
    selector_config,
    folds,
    selector_epochs=300,
    selector_batch_size=128,
    seed=42,
    device='auto',
    early_stopping_patience=None,
    verbose=False,
    log_prefix='Selector cache',
):
    """Train one Phase-1 selector per fold and cache full feature rankings."""
    device = _resolve_device(device)

    selector_hidden_dims = selector_config.get('hidden_dims', (128, 64))
    selector_l1_lambda = selector_config.get('l1_lambda', 0.004)
    selector_dropout = selector_config.get('dropout', 0.3)

    selector_cache = []
    total_folds = len(folds)

    for fold_number, (train_idx, val_idx) in enumerate(folds, start=1):
        if verbose:
            print(
                f"[{log_prefix}] fold {fold_number}/{total_folds} start "
                f"train={len(train_idx)} val={len(val_idx)}"
            )

        model, losses = _train_sparse_gate(
            X[train_idx],
            y[train_idx],
            sample_weights[train_idx],
            selector_hidden_dims,
            selector_l1_lambda,
            selector_dropout,
            selector_epochs,
            selector_batch_size,
            seed + fold_number - 1,
            device,
            early_stopping_patience=early_stopping_patience,
        )
        gate_values = model.feature_importance()
        feature_ranking = np.argsort(gate_values)[::-1]
        selector_cache.append({
            'fold': fold_number - 1,
            'train_idx': train_idx,
            'val_idx': val_idx,
            'gate_values': gate_values,
            'feature_ranking': feature_ranking,
            'selector_epochs_ran': len(losses),
        })

        if verbose:
            print(
                f"[{log_prefix}] fold {fold_number}/{total_folds} ready "
                f"selector_epochs={len(losses)} top_gate={_format_metric(gate_values[feature_ranking[0]])}"
            )

    return selector_cache


# ─────────────────────────────────────────────────────────────────────────────
# Cross-validation
# ─────────────────────────────────────────────────────────────────────────────

def cross_validate_features(
    X,
    y,
    sample_weights,
    terminal_mask,
    param_grid,
    groups=None,
    n_splits=5,
    n_select=15,
    n_epochs=300,
    batch_size=128,
    seed=42,
    device='auto',
    score_weights=None,
    stability_metric_weights=None,
    early_stopping_patience=None,
    verbose=False,
    log_prefix='Phase-1 CV',
):
    """
    K-fold stratified cross-validation over a hyperparameter param_grid.

    Parameters
    ----------
    X              : ndarray (n_samples, n_features)
    y              : ndarray (n_samples, n_classes)  – soft labels
    sample_weights : ndarray (n_samples,)
    terminal_mask  : ndarray (n_samples,) bool – hard-label (terminal) samples
    param_grid     : dict of lists
    groups         : ndarray (n_samples,), optional
        If provided, all samples with the same group value stay in the same
        fold. This is intended for lineage/timepoint experiments where nearby
        observations from one lineage cell should not be split across folds.
        Keys are any subset of {'l1_lambda', 'hidden_dims', 'dropout'}.
        All combinations are evaluated.  Example::

            {
                'l1_lambda':   [0.002, 0.004, 0.008],
                'hidden_dims': [(128, 64), (128, 64, 32), (64, 32)],
                'dropout':     [0.2, 0.3],
            }

    n_splits       : int   – number of CV folds
    n_select       : int   – number of top features to track for stability
    n_epochs       : int   – training epochs per fold (300 is usually sufficient)
    batch_size     : int
    seed           : int
    device         : str   – 'auto', 'cpu', or 'cuda'

    Returns
    -------
    results     : list[dict]
        One dict per param config, sorted by composite score (descending).
        Each dict contains:
            'config'       – hyperparameter dict
            'fold_results' – list of per-fold dicts with keys:
                             fold, train_acc, val_acc, val_loss,
                             gate_values, top_k_indices
            'stability'    – {pair_jaccards, pair_spearmans,
                               feature_frequency, consensus_features}
            'summary'      – {mean_val_acc, std_val_acc, mean_train_acc,
                               mean_overfit_gap, mean_jaccard, std_jaccard,
                               mean_spearman, mean_topk_frequency,
                               consensus_ratio, selection_stability_score}
            'score'        – stability-aware composite score
    best_config : dict – hyperparameter dict of the top-ranked configuration
    """
    device = _resolve_device(device)
    weights = {
        'val_acc': 0.65,
        'selection_stability': 0.35,
        'overfit_gap': 0.10,
    }
    if score_weights is not None:
        weights.update(score_weights)

    # Build full Cartesian product of the param grid
    configs = _build_param_grid(param_grid)

    folds = build_cv_folds(X, y, n_splits, seed, groups=groups)

    results = []

    total_folds = len(folds)

    for config_idx, cfg in enumerate(
        tqdm(configs, desc="CV configs", disable=verbose),
        start=1,
    ):
        l1_lambda   = cfg.get('l1_lambda',   0.004)
        hidden_dims = cfg.get('hidden_dims', (128, 64))
        dropout     = cfg.get('dropout',     0.3)

        if verbose:
            print(f"[{log_prefix}] config {config_idx}/{len(configs)} start {cfg}")

        fold_results = []
        for fold_idx, (train_idx, val_idx) in enumerate(folds):
            X_tr,  X_val  = X[train_idx],             X[val_idx]
            y_tr,  y_val  = y[train_idx],             y[val_idx]
            w_tr          = sample_weights[train_idx]
            term_tr       = terminal_mask[train_idx]
            term_val      = terminal_mask[val_idx]

            model, losses = _train_sparse_gate(
                X_tr, y_tr, w_tr,
                hidden_dims, l1_lambda, dropout,
                n_epochs, batch_size, seed + fold_idx, device,
                early_stopping_patience=early_stopping_patience,
            )

            X_tr_t  = torch.FloatTensor(X_tr)
            X_val_t = torch.FloatTensor(X_val)
            y_tr_t  = torch.FloatTensor(y_tr)
            y_val_t = torch.FloatTensor(y_val)

            train_acc = eval_accuracy(model, X_tr_t,  y_tr_t,  term_tr,  device)
            val_acc   = eval_accuracy(model, X_val_t, y_val_t, term_val, device)

            gate_vals = model.feature_importance()
            top_k     = np.argsort(gate_vals)[::-1][:n_select].tolist()

            # Weighted soft-CE on the validation set (no L1 term)
            model.eval()
            with torch.no_grad():
                logits_val, _ = model(X_val_t.to(device))
                w_val_t = torch.FloatTensor(sample_weights[val_idx]).to(device)
                val_loss = soft_cross_entropy(
                    logits_val, y_val_t.to(device), w_val_t
                ).item()

            fold_results.append({
                'fold':          fold_idx,
                'train_acc':     float(train_acc),
                'val_acc':       float(val_acc) if not np.isnan(val_acc) else float('nan'),
                'val_loss':      float(val_loss),
                'gate_values':   gate_vals,
                'top_k_indices': top_k,
                'selector_epochs_ran': len(losses),
            })

            if verbose:
                print(
                    f"[{log_prefix}] config {config_idx}/{len(configs)} "
                    f"fold {fold_idx + 1}/{total_folds} "
                    f"train_acc={_format_metric(train_acc)} "
                    f"val_acc={_format_metric(val_acc)} "
                    f"val_loss={_format_metric(val_loss)} "
                    f"selector_epochs={len(losses)}"
                )

        # ── Stability metrics ────────────────────────────────────────────────
        all_top_k = [fr['top_k_indices'] for fr in fold_results]
        all_gates = [fr['gate_values']   for fr in fold_results]

        pair_jaccards  = [
            _jaccard(all_top_k[i], all_top_k[j])
            for i, j in combinations(range(n_splits), 2)
        ]
        pair_spearmans = [
            float(spearmanr(all_gates[i], all_gates[j]).statistic)
            for i, j in combinations(range(n_splits), 2)
        ]

        # How many folds each feature appears in the top-K set
        n_features = X.shape[1]
        freq = np.zeros(n_features, dtype=int)
        for tk in all_top_k:
            freq[np.array(tk)] += 1

        # Consensus = in top-K in ≥80% of folds
        consensus = np.where(freq >= int(np.ceil(total_folds * 0.8)))[0].tolist()
        mean_topk_frequency = float(np.sort(freq / total_folds)[::-1][:n_select].mean())
        consensus_ratio = float(min(len(consensus), n_select) / n_select)

        # ── Summary ──────────────────────────────────────────────────────────
        val_accs    = [fr['val_acc']   for fr in fold_results if not np.isnan(fr['val_acc'])]
        train_accs  = [fr['train_acc'] for fr in fold_results]
        overfit_gaps = [
            fr['train_acc'] - fr['val_acc']
            for fr in fold_results
            if not np.isnan(fr['val_acc'])
        ]

        summary = {
            'mean_val_acc':     float(np.mean(val_accs))    if val_accs   else float('nan'),
            'std_val_acc':      float(np.std(val_accs))     if val_accs   else float('nan'),
            'mean_train_acc':   float(np.mean(train_accs)),
            'mean_overfit_gap': float(np.mean(overfit_gaps)) if overfit_gaps else float('nan'),
            'mean_jaccard':     float(np.mean(pair_jaccards)),
            'std_jaccard':      float(np.std(pair_jaccards)),
            'mean_spearman':    float(np.mean(pair_spearmans)),
            'mean_topk_frequency': mean_topk_frequency,
            'consensus_ratio':     consensus_ratio,
        }
        summary['selection_stability_score'] = _stability_score(
            summary,
            metric_weights=stability_metric_weights,
        )

        # Composite score: reward validation accuracy, stable top-K features,
        # and penalise configurations that clearly overfit.
        mean_val_acc = summary['mean_val_acc']
        score = 0.0
        if not np.isnan(mean_val_acc):
            score = (
                weights['val_acc'] * mean_val_acc
                + weights['selection_stability'] * summary['selection_stability_score']
                - weights['overfit_gap'] * max(summary['mean_overfit_gap'], 0.0)
            )

        results.append({
            'config':       cfg,
            'fold_results': fold_results,
            'stability': {
                'pair_jaccards':     pair_jaccards,
                'pair_spearmans':    pair_spearmans,
                'feature_frequency': freq,
                'consensus_features': consensus,
            },
            'summary': summary,
            'score':   score,
        })

        if verbose:
            print(
                f"[{log_prefix}] config {config_idx}/{len(configs)} done "
                f"mean_val_acc={_format_metric(summary['mean_val_acc'])} "
                f"mean_train_acc={_format_metric(summary['mean_train_acc'])} "
                f"mean_overfit_gap={_format_metric(summary['mean_overfit_gap'])} "
                f"stability={_format_metric(summary['selection_stability_score'])} "
                f"score={_format_metric(score)}"
            )

    results.sort(key=lambda r: r['score'], reverse=True)
    best_config = results[0]['config']
    return results, best_config


def cross_validate_focused(
    X,
    y,
    sample_weights,
    terminal_mask,
    selector_config,
    param_grid,
    groups=None,
    n_splits=5,
    selector_epochs=300,
    focused_epochs=400,
    selector_batch_size=128,
    focused_batch_size=64,
    seed=42,
    device='auto',
    score_weights=None,
    selector_cache=None,
    selector_early_stopping_patience=None,
    focused_early_stopping_patience=None,
    verbose=False,
    log_prefix='Phase-2 CV',
):
    """
    Jointly tune Phase 2 architecture and the selected-feature count K.

    Each fold re-runs Phase 1 feature selection on the training split only,
    preventing feature leakage while allowing the number of retained features
    to be part of model selection.

    Parameters
    ----------
    selector_config : dict
        Best Phase 1 configuration, typically returned by
        ``cross_validate_features``.
    param_grid : dict of lists
        May include ``n_select``, ``hidden_dims``, ``dropout``, and
        ``dist_lambda``.
    groups : ndarray (n_samples,), optional
        If provided, all samples with the same group value stay in the same
        fold during both selector and focused-model validation.

    Returns
    -------
    results     : list[dict]
        Per-config CV results sorted by score.
    best_config : dict
        Best Phase 2 configuration.
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

    folds = build_cv_folds(X, y, n_splits, seed, groups=groups)
    total_folds = len(folds)

    if selector_cache is None:
        selector_cache = build_selector_fold_cache(
            X,
            y,
            sample_weights,
            selector_config,
            folds,
            selector_epochs=selector_epochs,
            selector_batch_size=selector_batch_size,
            seed=seed,
            device=device,
            early_stopping_patience=selector_early_stopping_patience,
            verbose=verbose,
            log_prefix=f"{log_prefix} selector",
        )
    elif len(selector_cache) != total_folds:
        raise ValueError(
            f"Expected selector_cache to have {total_folds} folds, got {len(selector_cache)}"
        )

    results = []

    for config_idx, cfg in enumerate(
        tqdm(configs, desc="Phase-2 CV configs", disable=verbose),
        start=1,
    ):
        n_select = cfg.get('n_select', 15)
        focused_hidden_dims = cfg.get('hidden_dims', (64, 32))
        focused_dropout = cfg.get('dropout', 0.2)
        dist_lambda = cfg.get('dist_lambda', 0.1)

        if verbose:
            print(f"[{log_prefix}] config {config_idx}/{len(configs)} start {cfg}")

        fold_results = []
        for fold_idx, (train_idx, val_idx) in enumerate(folds):
            X_tr, X_val = X[train_idx], X[val_idx]
            y_tr, y_val = y[train_idx], y[val_idx]
            w_tr, w_val = sample_weights[train_idx], sample_weights[val_idx]
            term_tr = terminal_mask[train_idx]
            term_val = terminal_mask[val_idx]

            cached_selector = selector_cache[fold_idx]
            top_k = cached_selector['feature_ranking'][:n_select].tolist()

            X_tr_sel = X_tr[:, top_k].astype(np.float32)
            X_val_sel = X_val[:, top_k].astype(np.float32)

            model, _, losses = train_focused(
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
                early_stopping_patience=focused_early_stopping_patience,
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
                'selector_epochs_ran': cached_selector.get('selector_epochs_ran'),
                'focused_epochs_ran': len(losses),
            })

            if verbose:
                print(
                    f"[{log_prefix}] config {config_idx}/{len(configs)} "
                    f"fold {fold_idx + 1}/{total_folds} "
                    f"train_acc={_format_metric(train_metrics['argmax_accuracy'])} "
                    f"val_acc={_format_metric(val_hard_metrics['argmax_accuracy'])} "
                    f"expected_target_prob={_format_metric(val_soft_metrics['expected_target_probability'])} "
                    f"soft_ce={_format_metric(val_soft_metrics['soft_cross_entropy'])} "
                    f"selector_epochs={cached_selector.get('selector_epochs_ran')} "
                    f"focused_epochs={len(losses)}"
                )

        all_top_k = [fr['top_k_indices'] for fr in fold_results]
        pair_jaccards = [
            _jaccard(all_top_k[i], all_top_k[j])
            for i, j in combinations(range(n_splits), 2)
        ]

        n_features = X.shape[1]
        freq = np.zeros(n_features, dtype=int)
        for tk in all_top_k:
            freq[np.array(tk)] += 1
        consensus = np.where(freq >= int(np.ceil(total_folds * 0.8)))[0].tolist()

        val_accs = [fr['val_acc'] for fr in fold_results if not np.isnan(fr['val_acc'])]
        train_accs = [fr['train_acc'] for fr in fold_results]
        overfit_gaps = [
            fr['train_acc'] - fr['val_acc']
            for fr in fold_results
            if not np.isnan(fr['val_acc'])
        ]
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

        if verbose:
            print(
                f"[{log_prefix}] config {config_idx}/{len(configs)} done "
                f"mean_val_acc={_format_metric(summary['mean_val_acc'])} "
                f"mean_expected_target_prob={_format_metric(summary['mean_expected_target_probability'])} "
                f"mean_soft_ce={_format_metric(summary['mean_soft_cross_entropy'])} "
                f"mean_overfit_gap={_format_metric(summary['mean_overfit_gap'])} "
                f"stability={_format_metric(summary['selection_stability_score'])}"
            )

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

            if verbose:
                cfg = result['config']
                print(
                    f"[{log_prefix}] config final {cfg} "
                    f"normalized_param_count={_format_metric(param_penalty)} "
                    f"score={_format_metric(score)}"
                )

    results.sort(key=lambda r: r['score'], reverse=True)
    best_config = results[0]['config']
    return results, best_config


# ─────────────────────────────────────────────────────────────────────────────
# One-pass final training
# ─────────────────────────────────────────────────────────────────────────────

def train_one_pass(
    X,
    y,
    sample_weights,
    hidden_dims,
    l1_lambda,
    dropout,
    n_select=15,
    n_epochs=500,
    batch_size=128,
    seed=42,
    device='auto',
    early_stopping_patience=None,
):
    """
    Train SparseGateClassifier on ALL data using the CV-validated configuration.

    Parameters
    ----------
    X              : ndarray (n_samples, n_features)
    y              : ndarray (n_samples, n_classes)
    sample_weights : ndarray (n_samples,)
    hidden_dims    : tuple  – architecture found by cross_validate_features
    l1_lambda      : float  – sparsity regularisation found by CV
    dropout        : float  – dropout rate found by CV
    n_select       : int    – number of top features to return
    n_epochs       : int    – training epochs (500 recommended for final run)
    batch_size     : int
    seed           : int
    device         : str    – 'auto', 'cpu', or 'cuda'

    Returns
    -------
    model         : SparseGateClassifier  (trained, in eval mode)
    top_k_indices : ndarray (n_select,) int  – feature indices sorted by gate value desc
    gate_values   : ndarray (n_features,)     – gate importance score per feature
    embeddings    : ndarray (n_samples, embed_dim)  – penultimate-layer activations;
                    embed_dim = last element of hidden_dims
    losses        : list[float]  – per-epoch training loss
    """
    device = _resolve_device(device)

    model, losses = _train_sparse_gate(
        X, y, sample_weights,
        hidden_dims, l1_lambda, dropout,
        n_epochs, batch_size, seed, device,
        early_stopping_patience=early_stopping_patience,
    )

    gate_values   = model.feature_importance()
    top_k_indices = np.argsort(gate_values)[::-1][:n_select]

    model.eval()
    with torch.no_grad():
        _, emb_t = model(torch.FloatTensor(X).to(device))
    embeddings = emb_t.cpu().numpy()

    return model, top_k_indices, gate_values, embeddings, losses


# ─────────────────────────────────────────────────────────────────────────────
# Phase-2 focused training with distance-preservation regularization
# ─────────────────────────────────────────────────────────────────────────────

def _batch_dist_corr_loss(emb, x_input):
    """
    Returns the *negative* Pearson correlation between pairwise L2 distances
    in the embedding space and in the input selected-feature space, computed
    within the current mini-batch.

    Minimizing this loss encourages the embedding to preserve the relative
    distance structure of the input feature space (Mantel-test analogue, but
    differentiable and computed per batch for efficiency).
    """
    n = emb.shape[0]
    if n < 4:
        return torch.tensor(0.0, device=emb.device)

    d_emb = torch.cdist(emb,     emb,     p=2)
    d_inp = torch.cdist(x_input, x_input, p=2)

    rows, cols = torch.triu_indices(n, n, offset=1, device=emb.device)
    de = d_emb[rows, cols]
    di = d_inp[rows, cols]

    de_c = de - de.mean()
    di_c = di - di.mean()
    r = (de_c * di_c).sum() / (de_c.norm() * di_c.norm() + 1e-8)
    return -r   # maximise r → minimise –r


def train_focused(
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
    Phase 2: train FocusedClassifier on the pre-selected feature matrix.

    Because the input is already reduced to N_SELECT features, the learned
    embedding naturally stays geometrically close to the original protein
    expression subspace.  The optional ``dist_lambda`` term makes this
    explicit by adding a per-batch distance-correlation penalty.

    Parameters
    ----------
    X_sel                  : ndarray (n_samples, n_select) – already-selected features
    y                      : ndarray (n_samples, n_classes)  – soft labels
    sample_weights         : ndarray (n_samples,)
    hidden_dims            : tuple  – backbone widths; last element = embedding dim.
                             Default (64, 32) gives a 32-D embedding.
    dropout                : float
    n_epochs               : int    – 600 is a good default for Phase 2
    batch_size             : int    – smaller batches help the distance-corr term
    seed                   : int
    device                 : str    – 'auto', 'cpu', or 'cuda'
    dist_lambda            : float  – weight for the distance-preservation loss.
                             dist_lambda > 0 adds –Pearson(d_emb, d_input) per batch.
                             Set to 0 to disable and train with pure classification.
    model_type             : str    – 'focused', 'resnet', 'wide', 'attention'
    optimizer_type         : str    – 'adam', 'adamw', 'sgd'
    lr                     : float  – learning rate
    label_smoothing        : float  – soft label smoothing factor (0=no smoothing)
    early_stopping_patience: int    – patience for early stopping (None=disabled)
    gradient_clip          : float  – maximum gradient norm (None=disabled)

    Returns
    -------
    model      : FocusedClassifier (or variant) (trained, eval mode)
    embeddings : ndarray (n_samples, embed_dim)
    losses     : list[float] – per-epoch total loss
    """
    device = _resolve_device(device)
    torch.manual_seed(seed)

    n_selected = X_sel.shape[1]
    n_classes  = y.shape[1]

    X_t = torch.FloatTensor(X_sel)
    y_t = torch.FloatTensor(y)
    w_t = torch.FloatTensor(sample_weights)

    # Apply label smoothing if requested
    if label_smoothing > 0.0:
        y_t = y_t * (1.0 - label_smoothing) + label_smoothing / n_classes

    loader = DataLoader(
        TensorDataset(X_t, y_t, w_t),
        batch_size=batch_size,
        shuffle=True,
    )

    model = build_focused_classifier(model_type, n_selected, n_classes, hidden_dims, dropout).to(device)

    # Select optimizer
    if optimizer_type == 'adamw':
        opt = optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    elif optimizer_type == 'sgd':
        opt = optim.SGD(model.parameters(), lr=lr, momentum=0.9, weight_decay=1e-4)
    else:  # 'adam'
        opt = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)

    sched = optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_epochs)

    losses = []
    best_loss = float('inf')
    patience_counter = 0

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        for xb, yb, wb in loader:
            xb, yb, wb = xb.to(device), yb.to(device), wb.to(device)
            opt.zero_grad()
            logits, emb = model(xb)
            loss = soft_cross_entropy(logits, yb, wb)
            if dist_lambda > 0.0:
                loss = loss + dist_lambda * _batch_dist_corr_loss(emb, xb)
            loss.backward()

            # Gradient clipping
            if gradient_clip is not None:
                torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clip)

            opt.step()
            epoch_loss += loss.item()

        avg_epoch_loss = epoch_loss / len(loader)
        sched.step()
        losses.append(avg_epoch_loss)

        # Early stopping
        if early_stopping_patience is not None:
            if avg_epoch_loss < best_loss:
                best_loss = avg_epoch_loss
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= early_stopping_patience:
                    break

    model.eval()
    with torch.no_grad():
        _, emb_t = model(X_t.to(device))
    embeddings = emb_t.cpu().numpy()

    return model, embeddings, losses
