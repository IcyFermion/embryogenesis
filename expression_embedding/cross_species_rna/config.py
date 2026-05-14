"""Configuration dataclass for cross-species RNA embedding training."""

from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
BUNDLE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BUNDLE_DIR / "results"


@dataclass
class Config:
    # ---- Paths ----
    ele_rna_path: str = str(DATA_DIR / "c_briggsae" / "science.adu8249" / "c_elegans_tf.csv")
    bri_rna_path: str = str(DATA_DIR / "c_briggsae" / "science.adu8249" / "c_briggsae_tf.csv")
    lineage_path: str = str(DATA_DIR / "cell_lineage.json")
    cell_type_path: str = str(DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv")
    output_dir: str = str(RESULTS_DIR / "cross_species_rna_embedding")

    # ---- Model architecture ----
    hidden_dims: tuple = (128, 64)
    embed_dim: int = 32
    dropout: float = 0.1
    classifier_hidden: int = 64

    # ---- Loss coefficients ----
    alpha: float = 1.0  # reconstruction weight
    beta: float = 1.0   # classification weight

    # ---- Training ----
    lr: float = 1e-3
    weight_decay: float = 1e-4
    n_epochs: int = 300
    batch_size: int = 128
    patience: int = 30
    grad_clip: float = 1.0

    # ---- Split ----
    sublineage_depth: int = 5
    val_fraction: float = 0.2

    # ---- Data preprocessing ----
    log_transform: bool = True
    exclude_dead_cells: bool = False  # include programmed_death cells in training

    # ---- Misc ----
    seed: int = 42
    device: str = "auto"  # "auto", "cuda", "mps", "cpu"

    # ---- Schema ----
    merge_schema: str = "D"  # cell type merge schema (only D implemented)
