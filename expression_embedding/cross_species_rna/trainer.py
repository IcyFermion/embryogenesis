"""Training loop with mixed-species batches and per-species metrics."""

import copy

import numpy as np
import torch
import torch.nn.functional as F
from torch import optim

from .data_loader import MixedSpeciesBatchSampler
from .model import CrossSpeciesModel


def soft_cross_entropy(logits, soft_targets, weights=None):
    """Cross-entropy accepting soft (probabilistic) targets.

    Parameters
    ----------
    logits : Tensor (N, C)
    soft_targets : Tensor (N, C) — each row sums to 1
    weights : Tensor (N,) or None — per-sample weight

    Returns
    -------
    scalar loss
    """
    log_p = F.log_softmax(logits, dim=-1)
    per_sample = -(soft_targets * log_p).sum(dim=-1)
    if weights is not None:
        per_sample = per_sample * weights
    return per_sample.mean()


def _compute_metrics(logits, y, hard_mask, x_recon, x, species):
    """Compute per-species classification accuracy, recon MSE, and losses."""
    device = logits.device
    n_total = len(y)
    ele_mask = (species == 0)
    bri_mask = (species == 1)

    results = {"n_total": n_total}

    for label, mask in [("ele", ele_mask), ("bri", bri_mask)]:
        n = mask.sum().item()
        if n == 0:
            for key in ["acc", "recon_mse", "n"]:
                results[f"{label}_{key}"] = 0.0
            continue

        # Accuracy on hard-labeled cells
        hm_sp = hard_mask & mask
        n_hard = hm_sp.sum().item()
        if n_hard > 0:
            preds = logits[hm_sp].argmax(dim=1)
            targets = y[hm_sp].argmax(dim=1)
            acc = (preds == targets).float().mean().item()
        else:
            acc = 0.0

        recon_mse = F.mse_loss(x_recon[mask], x[mask]).item()
        results[f"{label}_acc"] = acc
        results[f"{label}_recon_mse"] = recon_mse
        results[f"{label}_n"] = int(n)

    return results


def train_epoch(model, X, y, species, sample_weights, optimizer, scheduler, config):
    """One training epoch with mixed-species batches."""
    model.train()
    device = next(model.parameters()).device

    sampler = MixedSpeciesBatchSampler(species, config.batch_size, shuffle=True,
                                       seed=config.seed)
    total_loss = 0.0
    total_recon = 0.0
    total_class = 0.0
    n_batches = 0

    for batch_idx in sampler:
        X_batch = X[batch_idx].to(device)
        y_batch = y[batch_idx].to(device)
        sw_batch = sample_weights[batch_idx].to(device)

        optimizer.zero_grad()
        z, x_recon, logits = model(X_batch)

        loss_recon = F.mse_loss(x_recon, X_batch)
        loss_class = soft_cross_entropy(logits, y_batch, weights=sw_batch)
        loss = config.alpha * loss_recon + config.beta * loss_class

        loss.backward()
        if config.grad_clip > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()

        total_loss += loss.item()
        total_recon += loss_recon.item()
        total_class += loss_class.item()
        n_batches += 1

    if scheduler is not None:
        scheduler.step()

    return {
        "total": total_loss / max(n_batches, 1),
        "recon": total_recon / max(n_batches, 1),
        "classify": total_class / max(n_batches, 1),
    }


@torch.no_grad()
def validate_epoch(model, X, y, species, sample_weights, hard_mask, config):
    """Validation with full forward pass (no batching needed for ~400 cells)."""
    model.eval()
    device = next(model.parameters()).device

    X_t = X.to(device)
    y_t = y.to(device)
    hm_t = hard_mask.to(device) if isinstance(hard_mask, torch.Tensor) else torch.BoolTensor(hard_mask).to(device)
    sp_t = species.to(device) if isinstance(species, torch.Tensor) else torch.LongTensor(species).to(device)

    z, x_recon, logits = model(X_t)

    loss_recon = F.mse_loss(x_recon, X_t)
    loss_class = soft_cross_entropy(logits, y_t,
                                    weights=sample_weights.to(device))
    loss = config.alpha * loss_recon + config.beta * loss_class

    metrics = _compute_metrics(logits, y_t, hm_t, x_recon, X_t, sp_t)
    metrics["total"] = loss.item()
    metrics["recon"] = loss_recon.item()
    metrics["classify"] = loss_class.item()
    metrics["joint_acc"] = (metrics.get("ele_acc", 0) + metrics.get("bri_acc", 0)) / 2.0
    return metrics


def train(model, data, config):
    """Full training loop with early stopping by joint validation accuracy.

    Returns
    -------
    model : CrossSpeciesModel (best checkpoint, on CPU)
    history : dict[str, list[float]]
    """
    device = (
        torch.device("cuda" if torch.cuda.is_available() else
                     "mps" if torch.backends.mps.is_available() else "cpu")
        if config.device == "auto" else torch.device(config.device)
    )
    print(f"Device: {device}")
    model = model.to(device)

    optimizer = optim.AdamW(model.parameters(), lr=config.lr,
                            weight_decay=config.weight_decay)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.n_epochs)

    X_train = data["X_train"]
    y_train = data["y_train"]
    species_train = data["species_train"]
    sw_train = data["sample_weights_train"]

    X_val = data["X_val"]
    y_val = data["y_val"]
    species_val = data["species_val"]
    sw_val = data["sample_weights_val"]
    hm_val = data["hard_mask_val"]

    history = {
        "train_total": [], "train_recon": [], "train_classify": [],
        "val_total": [], "val_recon": [], "val_classify": [],
        "val_joint_acc": [],
        "val_ele_acc": [], "val_bri_acc": [],
        "val_ele_recon_mse": [], "val_bri_recon_mse": [],
    }

    best_joint_acc = -1.0
    best_state = None
    best_epoch = 0
    patience_counter = 0

    for epoch in range(1, config.n_epochs + 1):
        train_m = train_epoch(model, X_train, y_train, species_train, sw_train,
                              optimizer, scheduler, config)
        val_m = validate_epoch(model, X_val, y_val, species_val, sw_val, hm_val, config)

        # Record
        history["train_total"].append(train_m["total"])
        history["train_recon"].append(train_m["recon"])
        history["train_classify"].append(train_m["classify"])
        history["val_total"].append(val_m["total"])
        history["val_recon"].append(val_m["recon"])
        history["val_classify"].append(val_m["classify"])
        history["val_joint_acc"].append(val_m["joint_acc"])
        history["val_ele_acc"].append(val_m["ele_acc"])
        history["val_bri_acc"].append(val_m["bri_acc"])
        history["val_ele_recon_mse"].append(val_m["ele_recon_mse"])
        history["val_bri_recon_mse"].append(val_m["bri_recon_mse"])

        joint_acc = val_m["joint_acc"]

        if joint_acc > best_joint_acc:
            best_joint_acc = joint_acc
            best_state = copy.deepcopy(model.state_dict())
            best_epoch = epoch
            patience_counter = 0
        else:
            patience_counter += 1

        if epoch % 10 == 0 or epoch == 1 or patience_counter == 0:
            print(f"Epoch {epoch:3d} | "
                  f"train {train_m['total']:.4f} (R:{train_m['recon']:.4f} C:{train_m['classify']:.4f}) | "
                  f"val {val_m['total']:.4f} (R:{val_m['recon']:.4f} C:{val_m['classify']:.4f}) | "
                  f"joint_acc {joint_acc:.4f} (ele:{val_m['ele_acc']:.4f} bri:{val_m['bri_acc']:.4f})")

        if patience_counter >= config.patience:
            print(f"Early stopping at epoch {epoch} (best joint_acc={best_joint_acc:.4f} at epoch {best_epoch})")
            break

    model.load_state_dict(best_state)
    model = model.cpu()
    print(f"Best model: epoch {best_epoch}, joint_acc={best_joint_acc:.4f}")
    return model, history
