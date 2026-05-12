"""Config dataclass for mRNA-to-protein embedding alignment."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

BUNDLE_DIR = Path(__file__).resolve().parent.parent
ROOT_DIR = BUNDLE_DIR.parent
DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = BUNDLE_DIR / "results"


@dataclass
class Config:
    # Paths
    rna_path: str = str(
        DATA_DIR / "c_briggsae" / "science.adu8249" / "c_elegans_tf.csv"
    )
    protein_emb_path: str = str(
        RESULTS_DIR / "timepoint_embedding_all_features" / "cell_embeddings_mean.csv"
    )
    lineage_path: str = str(DATA_DIR / "cell_lineage.json")
    cell_type_path: str = str(DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv")
    output_dir: str = str(RESULTS_DIR / "rna_protein_align")

    # Model architecture
    hidden_dims: tuple = (128, 64)
    embed_dim: int = 32
    dropout: float = 0.1
    use_layer_norm: bool = True

    # Loss coefficients
    alpha: float = 1.0  # alignment weight
    beta: float = 0.1  # reconstruction weight

    # Training
    lr: float = 1e-3
    weight_decay: float = 1e-4
    n_epochs: int = 200
    patience: int = 25
    grad_clip: float = 1.0

    # Pair sampling
    n_pairs_per_epoch: int = 10000
    near_fraction: float = 0.5  # fraction from bottom quartile
    val_n_pairs: int = 5000  # fixed validation pairs per epoch

    # Split
    sublineage_depth: int = 5
    val_fraction: float = 0.2

    # Data preprocessing
    log_transform: bool = True  # log1p mRNA values before L2-normalization

    # Misc
    seed: int = 42
    device: str = "auto"
