#!/usr/bin/env python3
"""
Quick Stage 1 test - verify architectures work
"""

import numpy as np
import pandas as pd
import torch

from bundle_paths import ensure_repo_on_path

ensure_repo_on_path()

from feature_selection import build_focused_classifier

print("Testing new architectures...")

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Device: {device}\n')

# Test parameters
n_selected = 25
n_classes = 18
n_samples = 32
hidden_dims_configs = [
    (32,),
    (128, 64),
    (128, 128, 64),
    (256, 128, 64),
]

model_types = ['focused', 'resnet', 'wide', 'attention']

print("Testing model architectures with dummy data:\n")

for model_type in model_types:
    for hidden_dims in hidden_dims_configs[:2]:  # Test first 2 configs per type
        try:
            # Create model
            model = build_focused_classifier(
                model_type, n_selected, n_classes, hidden_dims, dropout=0.2
            ).to(device)

            # Test forward pass
            X = torch.randn(n_samples, n_selected).to(device)
            logits, emb = model(X)

            n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

            status = "✓"
            print(f"{status} {model_type:12} | hidden={str(hidden_dims):20} | "
                  f"logits={tuple(logits.shape)} | emb_dim={emb.shape[-1]} | "
                  f"params={n_params:,}")
        except Exception as e:
            print(f"✗ {model_type:12} | hidden={str(hidden_dims):20} | Error: {str(e)[:50]}")

print("\n✓ Architecture test complete!")
