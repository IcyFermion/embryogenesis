"""
feature_selection.py
────────────────────────────────────────────────────────────────────────────────
Reusable ML module for protein feature selection and embedding generation.

Public API
----------
cross_validate_features(X, y, sample_weights, terminal_mask, param_grid, ...)
    5-fold stratified CV over a hyperparameter param_grid.  Returns per-config
    metrics (val accuracy, feature stability/Jaccard, overfitting gap) plus the
    best config ranked by val_acc × jaccard composite score.

train_one_pass(X, y, sample_weights, hidden_dims, l1_lambda, dropout, ...)
    Trains SparseGateClassifier on ALL data with the CV-selected configuration.
    Returns: model, top-K feature indices, gate values, embeddings
    (n_samples × last_hidden_dim), and training loss history.

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
from sklearn.model_selection import StratifiedKFold
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
    model.eval()
    logits, _ = model(X_tensor.to(device))
    preds   = logits.argmax(dim=-1).cpu().numpy()
    targets = y_tensor.argmax(dim=-1).cpu().numpy()
    m = mask.numpy() if isinstance(mask, torch.Tensor) else mask
    if m.sum() == 0:
        return float('nan')
    return (preds[m] == targets[m]).mean()


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


def _train_sparse_gate(X, y, sample_weights, hidden_dims, l1_lambda, dropout,
                        n_epochs, batch_size, seed, device):
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
    for _ in range(n_epochs):
        loss = train_epoch(model, loader, opt, device, l1_lambda=l1_lambda)
        sched.step()
        losses.append(loss)

    return model, losses


# ─────────────────────────────────────────────────────────────────────────────
# Cross-validation
# ─────────────────────────────────────────────────────────────────────────────

def cross_validate_features(
    X,
    y,
    sample_weights,
    terminal_mask,
    param_grid,
    n_splits=5,
    n_select=15,
    n_epochs=300,
    batch_size=128,
    seed=42,
    device='auto',
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
                               mean_spearman}
            'score'        – composite score (val_acc × jaccard)
    best_config : dict – hyperparameter dict of the top-ranked configuration
    """
    device = _resolve_device(device)

    # Build full Cartesian product of the param grid
    keys   = list(param_grid.keys())
    combos = list(itertools.product(*[param_grid[k] for k in keys]))
    configs = [dict(zip(keys, combo)) for combo in combos]

    # Stratify on dominant cell type
    strat_labels = y.argmax(axis=1)
    skf   = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = list(skf.split(X, strat_labels))

    results = []

    for cfg in tqdm(configs, desc="CV configs"):
        l1_lambda   = cfg.get('l1_lambda',   0.004)
        hidden_dims = cfg.get('hidden_dims', (128, 64))
        dropout     = cfg.get('dropout',     0.3)

        fold_results = []
        for fold_idx, (train_idx, val_idx) in enumerate(folds):
            X_tr,  X_val  = X[train_idx],             X[val_idx]
            y_tr,  y_val  = y[train_idx],             y[val_idx]
            w_tr          = sample_weights[train_idx]
            term_tr       = terminal_mask[train_idx]
            term_val      = terminal_mask[val_idx]

            model, _ = _train_sparse_gate(
                X_tr, y_tr, w_tr,
                hidden_dims, l1_lambda, dropout,
                n_epochs, batch_size, seed + fold_idx, device,
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
            })

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
        consensus = np.where(freq >= int(np.ceil(n_splits * 0.8)))[0].tolist()

        # ── Summary ──────────────────────────────────────────────────────────
        val_accs    = [fr['val_acc']   for fr in fold_results if not np.isnan(fr['val_acc'])]
        train_accs  = [fr['train_acc'] for fr in fold_results]
        overfit_gaps = [tr - va for tr, va in zip(train_accs, val_accs)]

        summary = {
            'mean_val_acc':     float(np.mean(val_accs))    if val_accs   else float('nan'),
            'std_val_acc':      float(np.std(val_accs))     if val_accs   else float('nan'),
            'mean_train_acc':   float(np.mean(train_accs)),
            'mean_overfit_gap': float(np.mean(overfit_gaps)) if overfit_gaps else float('nan'),
            'mean_jaccard':     float(np.mean(pair_jaccards)),
            'std_jaccard':      float(np.std(pair_jaccards)),
            'mean_spearman':    float(np.mean(pair_spearmans)),
        }

        # Composite score: reward both high val accuracy and high feature stability
        score = summary['mean_val_acc'] * summary['mean_jaccard'] \
                if not np.isnan(summary['mean_val_acc']) else 0.0

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
):
    """
    Phase 2: train FocusedClassifier on the pre-selected feature matrix.

    Because the input is already reduced to N_SELECT features, the learned
    embedding naturally stays geometrically close to the original protein
    expression subspace.  The optional ``dist_lambda`` term makes this
    explicit by adding a per-batch distance-correlation penalty.

    Parameters
    ----------
    X_sel          : ndarray (n_samples, n_select) – already-selected features
    y              : ndarray (n_samples, n_classes)  – soft labels
    sample_weights : ndarray (n_samples,)
    hidden_dims    : tuple  – backbone widths; last element = embedding dim.
                     Default (64, 32) gives a 32-D embedding.
    dropout        : float
    n_epochs       : int    – 600 is a good default for Phase 2
    batch_size     : int    – smaller batches help the distance-corr term
    seed           : int
    device         : str    – 'auto', 'cpu', or 'cuda'
    dist_lambda    : float  – weight for the distance-preservation loss.
                     dist_lambda > 0 adds –Pearson(d_emb, d_input) per batch.
                     Set to 0 to disable and train with pure classification.

    Returns
    -------
    model      : FocusedClassifier (trained, eval mode)
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

    loader = DataLoader(
        TensorDataset(X_t, y_t, w_t),
        batch_size=batch_size,
        shuffle=True,
    )

    model = FocusedClassifier(n_selected, n_classes, hidden_dims, dropout).to(device)
    opt   = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_epochs)

    losses = []
    for _ in range(n_epochs):
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
            opt.step()
            epoch_loss += loss.item()
        sched.step()
        losses.append(epoch_loss / len(loader))

    model.eval()
    with torch.no_grad():
        _, emb_t = model(X_t.to(device))
    embeddings = emb_t.cpu().numpy()

    return model, embeddings, losses
