### Task Description: Stage 1 — C. elegans mRNA Encoder with Protein-Embedding Alignment

### previous stage
We have done a dual head encoder based protein time point embedding with all c elegans protein atlas data, in ``expression_embedding/timepoint_embedding.py``. The next stage is too align rna-seq embeddings from two species to this protein-based embedding.

#### Objective
Train a per-cell mRNA encoder on C. elegans that produces embeddings whose pairwise cosine distances match those of the frozen per-timepoint protein embedding from the previous stage. No classification, no lineage-relational supervision.
#### Reference Artifacts

- Frozen protein-embedding model and its per-cell summarized embeddings from the previous stage. Use these as fixed alignment targets; do not update.
- C. elegans mRNA atlas (Large et al.), per-cell.
- TF ortholog list shared between C. elegans and C. briggsae. Freeze this feature set now; it will be reused in Stage 2.
- expression_comparison.ipynb for downstream sanity checks.

#### Data and Features

- c. elegans embedding results: ``expression_embedding/results/timepoint_embedding_all_features/cell_embeddings_mean.csv``
- c. elegans rna-seq data: ``data/c_briggsae/science.adu8249/c_elegans_tf.csv``
- c. briggsae rna-seq data: ``data/c_briggsae/science.adu8249/c_briggsae_tf.csv``, the TFs in briggsae are all present in elegans data, so we can just use TFs here as the shared RNA-seq feature across species.
- One sample = one C. elegans lineage cell with a per-cell mRNA vector restricted to shared TF orthologs.
- Pair samples: ~985 cells with both mRNA and a protein embedding.
- Apply L2 normalization to input mRNA vectors before the encoder.

#### Model

- Small MLP encoder, 2–3 hidden layers, embedding dim 32 (configurable).
- L2-normalize the encoder output so cosine distance is well-defined.
- Mirror decoder for reconstruction.
- All architecture hyperparameters configurable.

#### Losses
Total loss = α * L_align + β * L_recon, with defaults α=1.0, β=0.1.

- L_align: for sampled pairs (a, b), (d_cos_mRNA(a, b) - d_cos_protein(a, b))^2. Protein-side distances precomputed once and cached.
- L_recon: MSE between input mRNA vector and reconstruction.
- Both coefficients configurable for sweeps.

#### Pair Sampling

- Sample ~10K pairs per epoch (configurable).
- Bias toward near pairs in protein-embedding space: e.g., 50% of pairs sampled from the bottom quartile of protein distances, 50% uniform. Make the bias scheme configurable.
- Gradients must not flow into the protein embeddings.

#### Train/Validation Split

- Split by sub-lineage, same scheme as the previous stage. Reuse that utility.
- Held-out cells contribute neither to pair sampling nor to reconstruction.

#### Training Loop
Track per epoch:

- Total loss and each component (train + val).
- Pearson correlation between mRNA-encoder cosine distances and protein-embedding cosine distances on held-out pairs. This is the headline metric.
- Reconstruction MSE on held-out cells.

Save best model by validation alignment correlation, not total loss.
#### Outputs
In a configurable output directory under ``expression_embedding/results``:

- Trained encoder checkpoint.
Per-cell mRNA embeddings for all C. elegans cells (parquet or HDF5).
- Training curves (CSV + PNG): all loss components and the alignment correlation.
#### Diagnostic plots:

- Scatter of protein-distance vs. mRNA-encoder-distance on held-out pairs.
- Same-type vs. different-type sibling-pair distance histograms.
- Embedding-distance vs. lineage-tree-distance scatter.
- 2D UMAP of mRNA embeddings, colored by terminal type and developmental time.



#### Sanity Checks (Required)

- Linear-probe baseline replicated: confirm raw mRNA → protein-embedding R² ≈ 0.33 with cosine on the same data and split. Trained encoder should clearly exceed this.
- Encoder outputs are L2-normalized.
- Pair sampler produces no train/val leakage.
- Protein embeddings are frozen (gradient check).
- Validation alignment correlation lands in the 0.4–0.55 range. Below 0.3 → investigate. Above 0.6 → suspect leakage.

Out of Scope

- Classification heads.
- Any lineage-relational loss.
- C. briggsae data (Stage 2).
- Per-timepoint inputs (mRNA atlas is per-cell).
- Fine-tuning the protein embedding.

Code Organization
Self-contained script or notebook with config (YAML or dataclass) under ``expression_embedding/``. Separate modules for data loading, model, training, evaluation. Do not modify previous stage's code.