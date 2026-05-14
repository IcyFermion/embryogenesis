### Objective
Train a shared mRNA encoder on C. elegans and C. briggsae jointly, producing per-cell embeddings in a common latent space. No protein alignment, no lineage-relational supervision.
### Reference Artifacts

- C. elegans and C.briggsae mRNA atlas (Large et al.).
- scRNA-seq data with marker-gene-derived lineage labels. ``data/c_briggsae/science.adu8249/c_elegans_tf.csv`` and ``data/c_briggsae/science.adu8249/c_briggsae_tf.csv``
- Shared TF ortholog list, taken as the shared TFs from the scRNA-seq data.
- Terminal-type label set shared across species. ``data/2023-06-29_entropy_cell_key_V2.csv`
- smaller cell types should be merged together according to either schema A or D in ``expression_embedding/cell_type_merge.md``, default option should be schema D.
- look up ``expression_embedding/protein_feature_select.ipynb`` and ``expression_embedding/feature_selection.py`` as reference on how to do cell type encoding. Important things to follow are the hard one-hot encoding for terminal cells with types, soft labels and associated decaying sample weights for intermediate cells; option to include or exclude programmed death cells in training.

### Data and Features

- Sample = one lineage cell (per-cell, not per-timepoint) with a TF-ortholog mRNA vector.
- Input preprocessing: log-transform if not already applied, then per-species z-scoring (compute mean/std within each species independently, apply within-species).
- Species label retained as metadata for batching and evaluation, not used as model input.

### Model

- Shared MLP encoder, 2–3 hidden layers, embedding dim 32 (configurable).
- BatchNorm1d immediately after the input layer as a hedge against residual species shift.
- L2-normalize the encoder output.
- Mirror decoder for reconstruction.
- Classification head: linear or small MLP on the embedding, output dim = number of shared terminal types.
- No species-specific adapter. No alignment loss.

### Losses
Total = ``α * L_recon + β * L_classify``, defaults ``α=1.0, β=1.0``, both configurable.

- L_recon: MSE between input and reconstruction. Applied to both species.
- L_classify: cross-entropy against soft labels (terminal hard labels are one-hot soft labels; progenitors use the descendant-mixture soft labels constructed identically in both species). Applied to both species.
- No cross-species alignment loss. No lineage-relational loss.

### Batch Composition
Stratified mixed-species batches: each batch contains samples from both species in roughly proportional ratios to dataset sizes. This is required, not optional — species-pure batches would let BatchNorm and the encoder develop species-specific behavior.

### Train/Validation Split
Sub-lineage-aware, applied independently within each species. Reuse the branch-level hold-out utility from Stage 1. Track validation metrics per species.
### Training Loop
Track per epoch, split by species where applicable:

- Total loss and each component (train + val, per species).
- Classification accuracy on held-out hard-labeled cells (per species).
- Reconstruction MSE on held-out cells (per species).
- Cross-species same-type centroid distance (see diagnostics).

Save best model by joint validation classification accuracy averaged across species.

### Outputs
In configurable output directory (default to ``expression_embedding/results/cross_species_rna_embedding``):

- Trained encoder checkpoint.
- Per-cell embeddings for all C. elegans and C. briggsae cells, with species labels (parquet or HDF5).
- Training curves (CSV + PNG): all loss components and per-species metrics.
- Diagnostic plots and tables (see below).

### Diagnostics
Per species (intrinsic quality):

- Same-type vs. different-type sibling-pair distance histograms.
- Embedding distance vs. lineage-tree distance scatter.
- 2D PCA colored by terminal type and by developmental stage.

### Cross-species (the key checks):

- 2D PCA of joint embeddings, two-panel: colored by species, colored by terminal type. The type panel should show same-type cells from both species co-clustering.
- Per-type centroid table: for each terminal type, compute centroids in each species, report d_same (same type, cross-species) and compare against d_diff_within (different types, same species). Should preserve the 14/14 result from the input-level diagnostic.
- Confusion matrix of classifier on held-out cells, per species, against shared label set.

### Sanity Checks (Required)

1. Per-species z-scoring applied correctly (each species independently).
2. Mixed-species batches confirmed at the data-loader level.
3. Sub-lineage split has no train/val leakage in either species.
4. Cross-species same-type centroid distance is smaller than within-species different-type centroid distance for the majority of types (target: same as input-level result, 14/14 or close).
5. Validation classification accuracy on held-out C. elegans branches is comparable to (or better than) the per-timepoint protein model's accuracy. Far below → investigate. Far above → suspect C. briggsae labeling-circularity leakage into the joint training.

### Out of Scope

- Protein-embedding alignment loss.
- Species-specific input adapters (skipped per diagnostic).
- Lineage-relational losses.
- Per-timepoint inputs.
- Domain-adversarial training (not needed per diagnostic).

### Code Organization
Self-contained script or notebook with config. Separate modules for data loading, model, training, evaluation. Inherit utilities from ``expression_embedding`` as you see fit, don't change the existing experiments runtime code.