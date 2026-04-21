"""
autoregressive_embedding.py
────────────────────────────────────────────────────────────────────────────────
Unsupervised autoregressive GRU model for learning per-cell protein expression
embeddings from temporal sequences.

Public API
----------
load_temporal_sequences(h5ad_path, qc_threshold, min_seq_len)
    Load and preprocess the AnnData into per-cell temporal sequences.

train_autoregressive(sequences, ...)
    Train a 2-layer GRU to predict next-step expression.  Returns the trained
    model, training history, and copy-baseline MSE.

extract_embeddings(model, sequences, norm_params, method, embed_dim)
    Extract fixed-size per-cell embeddings from GRU hidden states via
    mean/last/weighted pooling + PCA.

Models
------
AutoregressiveGRU  – Linear projection → 2-layer GRU → linear prediction head.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence, pad_sequence
from torch.utils.data import DataLoader, Dataset
from sklearn.decomposition import PCA
from tqdm import tqdm
import anndata as ad


# ─────────────────────────────────────────────────────────────────────────────
# Data Loading
# ─────────────────────────────────────────────────────────────────────────────

def load_temporal_sequences(h5ad_path, qc_threshold=3.2, min_seq_len=3):
    """Load h5ad and build per-cell temporal sequences.

    Parameters
    ----------
    h5ad_path : str
        Path to the aggregated_scanpy.h5ad file.
    qc_threshold : float
        Values below this are treated as null signal and set to 0.
    min_seq_len : int
        Cells with fewer time points than this are dropped.

    Returns
    -------
    sequences : list of dict
        Each dict has keys 'cell_name' (str), 'times' (ndarray of int),
        'expression' (ndarray of shape (T, n_features), normalized).
    norm_params : dict
        'mean', 'std' arrays for inverse-transforming, plus 'qc_threshold'.
    feature_names : ndarray of str
        TF protein names (columns of the expression matrix).
    """
    adata = ad.read_h5ad(h5ad_path)
    import scipy.sparse
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else np.array(adata.X)
    X = X.astype(np.float64)

    # QC: values below threshold are null signal
    X[X < qc_threshold] = 0.0

    # Normalize: log1p then per-feature z-score, clip to [-5, 5]
    X = np.log1p(X)
    feat_mean = X.mean(axis=0)
    feat_std = X.std(axis=0)
    feat_std[feat_std == 0] = 1.0  # avoid division by zero
    X = (X - feat_mean) / feat_std
    X = np.clip(X, -5.0, 5.0)

    norm_params = {
        "mean": feat_mean,
        "std": feat_std,
        "qc_threshold": qc_threshold,
    }
    feature_names = np.array(adata.var.index)

    # Group by Cell-name, sort by Time
    cell_names_col = adata.obs["Cell-name"].values
    times_col = adata.obs["Time"].values

    unique_cells = np.unique(cell_names_col)
    sequences = []
    for cell in unique_cells:
        mask = cell_names_col == cell
        times = times_col[mask]
        expr = X[mask]
        order = np.argsort(times)
        times = times[order]
        expr = expr[order]
        if len(times) < min_seq_len:
            continue
        sequences.append({
            "cell_name": cell,
            "times": times,
            "expression": expr.astype(np.float32),
        })

    return sequences, norm_params, feature_names


# ─────────────────────────────────────────────────────────────────────────────
# Dataset
# ─────────────────────────────────────────────────────────────────────────────

class CellTemporalDataset(Dataset):
    """PyTorch dataset that returns (input_seq, target_seq, length) per cell."""

    def __init__(self, sequences):
        self.sequences = sequences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        expr = self.sequences[idx]["expression"]  # (T, n_feat)
        input_seq = torch.from_numpy(expr[:-1])   # (T-1, n_feat)
        target_seq = torch.from_numpy(expr[1:])    # (T-1, n_feat)
        length = input_seq.shape[0]
        return input_seq, target_seq, length


def collate_fn(batch):
    """Pad variable-length sequences and return (inputs, targets, lengths)."""
    inputs, targets, lengths = zip(*batch)
    lengths = torch.tensor(lengths, dtype=torch.long)
    # Sort by length descending for pack_padded_sequence
    sorted_idx = torch.argsort(lengths, descending=True)
    inputs = pad_sequence([inputs[i] for i in sorted_idx], batch_first=True)
    targets = pad_sequence([targets[i] for i in sorted_idx], batch_first=True)
    lengths = lengths[sorted_idx]
    return inputs, targets, lengths


# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────

class AutoregressiveGRU(nn.Module):
    """2-layer GRU with linear input projection and prediction head.

    Parameters
    ----------
    n_features : int
        Number of input features (TF proteins).
    d_model : int
        Hidden dimension for both the projection and GRU.
    n_layers : int
        Number of GRU layers.
    dropout : float
        Dropout between GRU layers.
    """

    def __init__(self, n_features, d_model=128, n_layers=2, dropout=0.1):
        super().__init__()
        self.input_proj = nn.Sequential(
            nn.Linear(n_features, d_model),
            nn.LayerNorm(d_model),
            nn.ReLU(),
        )
        self.gru = nn.GRU(
            input_size=d_model,
            hidden_size=d_model,
            num_layers=n_layers,
            dropout=dropout if n_layers > 1 else 0.0,
            batch_first=True,
        )
        self.pred_head = nn.Linear(d_model, n_features)

    def forward(self, x, lengths=None):
        """
        Parameters
        ----------
        x : Tensor (batch, seq_len, n_features)
        lengths : LongTensor (batch,), optional

        Returns
        -------
        pred : Tensor (batch, seq_len, n_features)
        hidden : Tensor (batch, seq_len, d_model)
        """
        h = self.input_proj(x)
        if lengths is not None:
            h = pack_padded_sequence(h, lengths.cpu(), batch_first=True, enforce_sorted=True)
        h, _ = self.gru(h)
        if lengths is not None:
            h, _ = pad_packed_sequence(h, batch_first=True)
        pred = self.pred_head(h)
        return pred, h


# ─────────────────────────────────────────────────────────────────────────────
# Training
# ─────────────────────────────────────────────────────────────────────────────

def _masked_mse_loss(pred, target, lengths):
    """Compute MSE only on non-padded positions."""
    batch, max_len, n_feat = pred.shape
    mask = torch.zeros(batch, max_len, device=pred.device, dtype=torch.bool)
    for i, l in enumerate(lengths):
        mask[i, :l] = True
    mask = mask.unsqueeze(-1).expand_as(pred)
    diff_sq = (pred - target) ** 2
    return (diff_sq * mask).sum() / mask.sum()


def _resolve_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def _compute_copy_baseline(sequences):
    """MSE of naive copy baseline: predict x_{t+1} = x_t."""
    total_sq_err = 0.0
    total_count = 0
    for seq in sequences:
        expr = seq["expression"]  # (T, F)
        diff = expr[1:] - expr[:-1]
        total_sq_err += (diff ** 2).sum()
        total_count += diff.size
    return total_sq_err / total_count


def train_autoregressive(
    sequences,
    n_features=None,
    d_model=128,
    n_layers=2,
    dropout=0.1,
    n_epochs=200,
    batch_size=32,
    lr=1e-3,
    weight_decay=1e-4,
    grad_clip=1.0,
    val_fraction=0.1,
    patience=20,
    seed=42,
    device=None,
):
    """Train the autoregressive GRU model.

    Parameters
    ----------
    sequences : list of dict
        Output of load_temporal_sequences.
    n_features : int, optional
        Inferred from data if None.
    d_model, n_layers, dropout : model hyperparameters
    n_epochs, batch_size, lr, weight_decay, grad_clip : training hyperparams
    val_fraction : float
        Fraction of cells held out for validation.
    patience : int
        Early stopping patience (epochs without val improvement).
    seed : int
        Random seed for reproducibility.
    device : torch.device, optional

    Returns
    -------
    model : AutoregressiveGRU (eval mode, on CPU)
    history : dict with 'train_loss', 'val_loss' lists
    copy_baseline : float
    """
    if device is None:
        device = _resolve_device()
    rng = np.random.RandomState(seed)

    if n_features is None:
        n_features = sequences[0]["expression"].shape[1]

    # Train/val split on cells
    n_val = max(1, int(len(sequences) * val_fraction))
    idx = rng.permutation(len(sequences))
    val_idx, train_idx = idx[:n_val], idx[n_val:]
    train_seqs = [sequences[i] for i in train_idx]
    val_seqs = [sequences[i] for i in val_idx]

    copy_baseline = _compute_copy_baseline(sequences)

    train_ds = CellTemporalDataset(train_seqs)
    val_ds = CellTemporalDataset(val_seqs)
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              collate_fn=collate_fn, drop_last=False)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            collate_fn=collate_fn, drop_last=False)

    model = AutoregressiveGRU(n_features, d_model, n_layers, dropout).to(device)
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    history = {"train_loss": [], "val_loss": []}
    best_val_loss = float("inf")
    best_state = None
    wait = 0

    for epoch in range(n_epochs):
        # --- Train ---
        model.train()
        total_loss, n_batches = 0.0, 0
        for inp, tgt, lengths in train_loader:
            inp, tgt = inp.to(device), tgt.to(device)
            optimizer.zero_grad()
            pred, _ = model(inp, lengths)
            loss = _masked_mse_loss(pred, tgt, lengths)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1
        scheduler.step()
        train_loss = total_loss / max(n_batches, 1)

        # --- Validate ---
        model.eval()
        total_val, n_val_b = 0.0, 0
        with torch.no_grad():
            for inp, tgt, lengths in val_loader:
                inp, tgt = inp.to(device), tgt.to(device)
                pred, _ = model(inp, lengths)
                total_val += _masked_mse_loss(pred, tgt, lengths).item()
                n_val_b += 1
        val_loss = total_val / max(n_val_b, 1)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)

        if epoch % 20 == 0 or epoch == n_epochs - 1:
            print(f"  Epoch {epoch:3d} | train {train_loss:.6f} | val {val_loss:.6f} | copy baseline {copy_baseline:.6f}")

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                print(f"  Early stopping at epoch {epoch} (patience={patience})")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    model.cpu()
    return model, history, copy_baseline


# ─────────────────────────────────────────────────────────────────────────────
# Embedding Extraction
# ─────────────────────────────────────────────────────────────────────────────

@torch.no_grad()
def extract_embeddings(model, sequences, method="mean", embed_dim=32, device=None):
    """Extract per-cell embeddings from GRU hidden states.

    Parameters
    ----------
    model : AutoregressiveGRU (eval mode)
    sequences : list of dict
    method : {'mean', 'last', 'weighted_mean'}
    embed_dim : int
        Final embedding dimensionality (via PCA).
    device : torch.device, optional

    Returns
    -------
    cell_names : list of str
    embeddings : ndarray (n_cells, embed_dim)
    pca : fitted PCA object
    """
    if device is None:
        device = _resolve_device()
    model = model.to(device)
    model.eval()

    all_hidden = []
    cell_names = []

    for seq_dict in sequences:
        expr = seq_dict["expression"]  # (T, F)
        x = torch.from_numpy(expr).unsqueeze(0).to(device)  # (1, T, F)
        lengths = torch.tensor([expr.shape[0]], dtype=torch.long)
        _, hidden_states = model(x, lengths)  # (1, T, d_model)
        h = hidden_states.squeeze(0).cpu().numpy()  # (T, d_model)

        if method == "mean":
            emb = h.mean(axis=0)
        elif method == "last":
            emb = h[-1]
        elif method == "weighted_mean":
            weights = np.linspace(0.5, 1.0, len(h))
            weights /= weights.sum()
            emb = (h * weights[:, None]).sum(axis=0)
        else:
            raise ValueError(f"Unknown method: {method}")

        all_hidden.append(emb)
        cell_names.append(seq_dict["cell_name"])

    embeddings_full = np.stack(all_hidden)  # (n_cells, d_model)
    pca = PCA(n_components=embed_dim)
    embeddings = pca.fit_transform(embeddings_full)

    model.cpu()
    return cell_names, embeddings, pca


# ─────────────────────────────────────────────────────────────────────────────
# Feature Importance
# ─────────────────────────────────────────────────────────────────────────────

def feature_importance(model):
    """Compute per-feature importance from the input projection weights.

    Returns
    -------
    importance : ndarray (n_features,)
        L2 norm of each input feature's column in the projection weight matrix.
    """
    weight = model.input_proj[0].weight.detach().cpu().numpy()  # (d_model, n_features)
    return np.linalg.norm(weight, axis=0)
