"""Training loop for mRNA-to-protein embedding alignment."""

import copy

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import pearsonr

from .config import Config
from .data_loader import PairSampler
from .model import RnaEncoderModel


def _resolve_device(device_str: str) -> torch.device:
    if device_str == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        elif torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device_str)


def compute_alignment_loss(
    embeddings: torch.Tensor,
    prot_dist_matrix: np.ndarray,
    pair_a: np.ndarray,
    pair_b: np.ndarray,
) -> torch.Tensor:
    """MSE between mRNA-encoder and protein cosine distances for sampled pairs.

    Protein distances come from numpy — no gradient flows into them.
    """
    idx_a = torch.as_tensor(pair_a, dtype=torch.long, device=embeddings.device)
    idx_b = torch.as_tensor(pair_b, dtype=torch.long, device=embeddings.device)

    cos_sim_mrna = (embeddings[idx_a] * embeddings[idx_b]).sum(dim=1)
    cos_dist_mrna = 1.0 - cos_sim_mrna

    cos_dist_protein = torch.as_tensor(
        prot_dist_matrix[pair_a, pair_b],
        dtype=torch.float32,
        device=embeddings.device,
    )

    return F.mse_loss(cos_dist_mrna, cos_dist_protein)


@torch.no_grad()
def pearson_on_pairs(
    model: RnaEncoderModel,
    X: np.ndarray,
    prot: np.ndarray,
    pair_a: np.ndarray,
    pair_b: np.ndarray,
    device: torch.device,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute Pearson r between mRNA-encoder and protein cosine distances.

    Uses pre-specified pair indices for consistent comparison across epochs.

    Returns (pearson_r, mrna_dists, prot_dists).
    """
    model.eval()
    X_t = torch.as_tensor(X, dtype=torch.float32, device=device)
    z = model.encoder(X_t)  # L2-normalized

    cos_sim = (z[pair_a] * z[pair_b]).sum(dim=1)
    mrna_dists = (1.0 - cos_sim).cpu().numpy()

    # Protein cosine distances
    prot_a = prot[pair_a]
    prot_b = prot[pair_b]
    prot_a_n = prot_a / np.clip(np.linalg.norm(prot_a, axis=1, keepdims=True), 1e-10, None)
    prot_b_n = prot_b / np.clip(np.linalg.norm(prot_b, axis=1, keepdims=True), 1e-10, None)
    prot_dists = np.clip(1.0 - (prot_a_n * prot_b_n).sum(axis=1), 0.0, 2.0)

    r, _ = pearsonr(mrna_dists, prot_dists)
    return r, mrna_dists, prot_dists


def make_fixed_val_pairs(n_val: int, n_pairs: int, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Generate fixed validation pair indices for consistent monitoring."""
    rng = np.random.RandomState(seed)
    idx_a = rng.randint(0, n_val, size=n_pairs * 2)
    idx_b = rng.randint(0, n_val, size=n_pairs * 2)
    same = idx_a == idx_b
    idx_a = idx_a[~same][:n_pairs]
    idx_b = idx_b[~same][:n_pairs]
    return idx_a, idx_b


def train_epoch(
    model: RnaEncoderModel,
    X_tensor: torch.Tensor,
    prot_dist_matrix: np.ndarray,
    pair_sampler: PairSampler,
    config: Config,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
) -> dict:
    model.train()

    z, x_recon = model(X_tensor)

    idx_a, idx_b = pair_sampler.sample(config.n_pairs_per_epoch)
    loss_align = compute_alignment_loss(z, prot_dist_matrix, idx_a, idx_b)
    loss_recon = F.mse_loss(x_recon, X_tensor)
    loss = config.alpha * loss_align + config.beta * loss_recon

    optimizer.zero_grad()
    loss.backward()
    nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
    optimizer.step()

    return {
        "total": loss.item(),
        "align": loss_align.item(),
        "recon": loss_recon.item(),
    }


@torch.no_grad()
def validate_epoch(
    model: RnaEncoderModel,
    X_val_tensor: torch.Tensor,
    prot_dist_matrix: np.ndarray,
    val_pair_sampler: PairSampler,
    config: Config,
    device: torch.device,
) -> dict:
    model.eval()

    z, x_recon = model(X_val_tensor)

    idx_a, idx_b = val_pair_sampler.sample(config.n_pairs_per_epoch)
    loss_align = compute_alignment_loss(z, prot_dist_matrix, idx_a, idx_b)
    loss_recon = F.mse_loss(x_recon, X_val_tensor)
    loss = config.alpha * loss_align + config.beta * loss_recon

    return {
        "total": loss.item(),
        "align": loss_align.item(),
        "recon": loss_recon.item(),
    }


def train(
    config: Config,
    X_train: np.ndarray,
    X_val: np.ndarray,
    prot_train: np.ndarray,
    prot_val: np.ndarray,
    prot_dist_matrix: np.ndarray,
    train_pair_sampler: PairSampler,
    val_pair_sampler: PairSampler,
    n_features: int,
) -> tuple[RnaEncoderModel, dict]:
    device = _resolve_device(config.device)
    print(f"  Device: {device}")

    model = RnaEncoderModel(
        n_features=n_features,
        hidden_dims=config.hidden_dims,
        embed_dim=config.embed_dim,
        dropout=config.dropout,
        use_layer_norm=config.use_layer_norm,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.lr, weight_decay=config.weight_decay
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.n_epochs
    )

    X_train_t = torch.as_tensor(X_train, dtype=torch.float32, device=device)
    X_val_t = torch.as_tensor(X_val, dtype=torch.float32, device=device)

    # Fixed validation pairs for consistent early stopping
    fixed_pairs_a, fixed_pairs_b = make_fixed_val_pairs(
        len(X_val), config.val_n_pairs, seed=config.seed
    )

    history = {
        "train_total": [],
        "train_align": [],
        "train_recon": [],
        "val_total": [],
        "val_align": [],
        "val_recon": [],
        "val_pearson": [],
    }

    best_val_pearson = -float("inf")
    best_state = None
    wait = 0

    for epoch in range(config.n_epochs):
        train_losses = train_epoch(
            model, X_train_t, prot_dist_matrix,
            train_pair_sampler, config, optimizer, device,
        )
        scheduler.step()

        val_losses = validate_epoch(
            model, X_val_t, prot_dist_matrix,
            val_pair_sampler, config, device,
        )

        val_pearson, _, _ = pearson_on_pairs(
            model, X_val, prot_val,
            fixed_pairs_a, fixed_pairs_b, device,
        )

        history["train_total"].append(train_losses["total"])
        history["train_align"].append(train_losses["align"])
        history["train_recon"].append(train_losses["recon"])
        history["val_total"].append(val_losses["total"])
        history["val_align"].append(val_losses["align"])
        history["val_recon"].append(val_losses["recon"])
        history["val_pearson"].append(val_pearson)

        if val_pearson > best_val_pearson:
            best_val_pearson = val_pearson
            best_state = copy.deepcopy(model.state_dict())
            wait = 0
        else:
            wait += 1

        if epoch % 20 == 0 or epoch == config.n_epochs - 1 or wait >= config.patience:
            print(
                f"  Epoch {epoch:3d} | "
                f"train total={train_losses['total']:.4f} "
                f"align={train_losses['align']:.4f} "
                f"recon={train_losses['recon']:.4f} | "
                f"val total={val_losses['total']:.4f} "
                f"pearson={val_pearson:.4f}"
            )

        if wait >= config.patience:
            print(f"  Early stopping at epoch {epoch}")
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model.cpu(), history
