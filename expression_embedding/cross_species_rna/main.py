"""Cross-species RNA embedding pipeline — orchestration and CLI."""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch

# Ensure sibling imports work
_here = Path(__file__).resolve().parents[0]
_bundle = _here.parent
if str(_bundle) not in sys.path:
    sys.path.insert(0, str(_bundle))

from .config import Config
from .data_loader import load_cross_species_data
from .evaluator import run_evaluation
from .model import CrossSpeciesModel
from .trainer import train


def main():
    parser = argparse.ArgumentParser(description="Cross-species RNA embedding")
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--lr", type=float, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--hidden-dims", type=str, default=None,
                        help="Comma-separated, e.g. '128,64'")
    parser.add_argument("--embed-dim", type=int, default=None)
    parser.add_argument("--dropout", type=float, default=None)
    parser.add_argument("--classifier-hidden", type=int, default=None)
    parser.add_argument("--alpha", type=float, default=None)
    parser.add_argument("--beta", type=float, default=None)
    parser.add_argument("--weight-decay", type=float, default=None)
    parser.add_argument("--patience", type=int, default=None)
    parser.add_argument("--depth", type=int, default=None,
                        help="Sublineage depth for split")
    parser.add_argument("--val-fraction", type=float, default=None)
    parser.add_argument("--exclude-dead-cells", action="store_true", default=None)
    parser.add_argument("--no-log-transform", action="store_true", default=None)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--device", type=str, default=None)
    args = parser.parse_args()

    config = Config()

    # Apply CLI overrides
    overrides = {
        "n_epochs": args.epochs,
        "lr": args.lr,
        "batch_size": args.batch_size,
        "embed_dim": args.embed_dim,
        "dropout": args.dropout,
        "classifier_hidden": args.classifier_hidden,
        "alpha": args.alpha,
        "beta": args.beta,
        "weight_decay": args.weight_decay,
        "patience": args.patience,
        "sublineage_depth": args.depth,
        "val_fraction": args.val_fraction,
        "output_dir": args.output,
        "seed": args.seed,
        "device": args.device,
    }
    for key, val in overrides.items():
        if val is not None:
            setattr(config, key, val)

    if args.hidden_dims is not None:
        config.hidden_dims = tuple(int(x) for x in args.hidden_dims.split(","))
    if args.no_log_transform:
        config.log_transform = False
    if args.exclude_dead_cells:
        config.exclude_dead_cells = True

    # Set seeds
    random.seed(config.seed)
    np.random.seed(config.seed)
    torch.manual_seed(config.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(config.seed)

    print("=" * 60)
    print("Cross-Species RNA Embedding")
    print("=" * 60)
    for k, v in config.__dict__.items():
        if not k.startswith("_"):
            print(f"  {k}: {v}")
    print()

    # ---- Load data ----
    print("Loading data...")
    data = load_cross_species_data(config)

    # ---- Create model ----
    n_features = data["n_features"]
    n_classes = len(data["class_names"])
    model = CrossSpeciesModel(
        n_features=n_features,
        n_classes=n_classes,
        hidden_dims=config.hidden_dims,
        embed_dim=config.embed_dim,
        classifier_hidden=config.classifier_hidden,
        dropout=config.dropout,
    )
    n_params = sum(p.numel() for p in model.parameters())
    print(f"\nModel: {n_params:,} parameters")
    print(f"  hidden_dims={config.hidden_dims}, embed_dim={config.embed_dim}")
    print(f"  n_features={n_features}, n_classes={n_classes}")

    # ---- Train ----
    model, history = train(model, data, config)

    # ---- Evaluate ----
    run_evaluation(model, data, history, config)


if __name__ == "__main__":
    main()
