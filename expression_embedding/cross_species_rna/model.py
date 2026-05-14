"""Cross-species shared encoder, mirror decoder, and classifier head."""

import torch
import torch.nn as nn
import torch.nn.functional as F


class CrossSpeciesEncoder(nn.Module):
    """Shared MLP encoder with BatchNorm1d after input layer.

    Input:  (batch, n_features) — preprocessed mRNA vector
    Output: (batch, embed_dim)  — L2-normalized embedding
    """

    def __init__(self, n_features, hidden_dims=(128, 64), embed_dim=32, dropout=0.1):
        super().__init__()
        layers = []
        in_dim = n_features

        for i, hd in enumerate(hidden_dims):
            layers.append(nn.Linear(in_dim, hd))
            if i == 0:
                layers.append(nn.BatchNorm1d(hd))  # hedge against species shift
            else:
                layers.append(nn.BatchNorm1d(hd))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
            in_dim = hd

        layers.append(nn.Linear(in_dim, embed_dim))
        self.net = nn.Sequential(*layers)
        self.embed_dim = embed_dim

    def forward(self, x):
        h = self.net(x)
        return F.normalize(h, p=2, dim=1)


class MirrorDecoder(nn.Module):
    """Reconstructs input from embedding (reverse of encoder, no BatchNorm)."""

    def __init__(self, n_features, hidden_dims=(128, 64), embed_dim=32, dropout=0.1):
        super().__init__()
        reversed_dims = list(reversed(hidden_dims))
        layers = []
        in_dim = embed_dim
        for hd in reversed_dims:
            layers.append(nn.Linear(in_dim, hd))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
            in_dim = hd
        layers.append(nn.Linear(in_dim, n_features))
        self.net = nn.Sequential(*layers)

    def forward(self, z):
        return self.net(z)


class Classifier(nn.Module):
    """Small MLP classification head on embedding."""

    def __init__(self, embed_dim, n_classes, hidden_dim=64, dropout=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_classes),
        )

    def forward(self, z):
        return self.net(z)


class CrossSpeciesModel(nn.Module):
    """Shared encoder + decoder + classifier for cross-species RNA embedding."""

    def __init__(self, n_features, n_classes, hidden_dims=(128, 64), embed_dim=32,
                 classifier_hidden=64, dropout=0.1):
        super().__init__()
        self.encoder = CrossSpeciesEncoder(n_features, hidden_dims, embed_dim, dropout)
        self.decoder = MirrorDecoder(n_features, hidden_dims, embed_dim, dropout)
        self.classifier = Classifier(embed_dim, n_classes, classifier_hidden, dropout)

    def forward(self, x):
        z = self.encoder(x)
        x_recon = self.decoder(z)
        logits = self.classifier(z)
        return z, x_recon, logits
