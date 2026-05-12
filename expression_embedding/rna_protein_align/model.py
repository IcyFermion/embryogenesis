"""MLP encoder with mirror decoder for mRNA-to-protein alignment."""

import torch
import torch.nn as nn
import torch.nn.functional as F


class RnaEncoder(nn.Module):
    """MLP encoder: L2-normalized mRNA vector -> L2-normalized embedding."""

    def __init__(
        self,
        n_features: int,
        hidden_dims: tuple = (128, 64),
        embed_dim: int = 32,
        dropout: float = 0.1,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        layers = []
        in_dim = n_features
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            if use_layer_norm:
                layers.append(nn.LayerNorm(h))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
            in_dim = h
        layers.append(nn.Linear(in_dim, embed_dim))
        self.net = nn.Sequential(*layers)
        self.embed_dim = embed_dim

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.net(x)
        return F.normalize(h, p=2, dim=1)


class MirrorDecoder(nn.Module):
    """Mirror of encoder: embedding -> reconstructed mRNA vector."""

    def __init__(
        self,
        n_features: int,
        hidden_dims: tuple = (128, 64),
        embed_dim: int = 32,
        dropout: float = 0.1,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        layers = []
        rev_dims = (embed_dim,) + tuple(reversed(hidden_dims))
        for i in range(len(rev_dims) - 1):
            layers.append(nn.Linear(rev_dims[i], rev_dims[i + 1]))
            if use_layer_norm:
                layers.append(nn.LayerNorm(rev_dims[i + 1]))
            layers.append(nn.GELU())
            layers.append(nn.Dropout(dropout))
        layers.append(nn.Linear(rev_dims[-1], n_features))
        self.decoder = nn.Sequential(*layers)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)


class RnaEncoderModel(nn.Module):
    """Container: encoder + decoder."""

    def __init__(
        self,
        n_features: int,
        hidden_dims: tuple = (128, 64),
        embed_dim: int = 32,
        dropout: float = 0.1,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        self.encoder = RnaEncoder(
            n_features, hidden_dims, embed_dim, dropout, use_layer_norm
        )
        self.decoder = MirrorDecoder(
            n_features, hidden_dims, embed_dim, dropout, use_layer_norm
        )
        self.embed_dim = embed_dim

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        z = self.encoder(x)
        x_recon = self.decoder(z)
        return z, x_recon
