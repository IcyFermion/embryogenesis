# Feature Embedding Push

This directory contains the recent feature-embedding experiment bundle that was previously spread across the repository root.

## Scope

The bundle covers three connected experiment threads:

- Protein feature selection and Phase 2 classifier tuning
- Autoregressive temporal embedding generation
- Expression-space comparison notebooks

Shared source data still lives under `data/` at the repository root. Bundle-generated artifacts now live under `expression_embedding/results/`.

## Layout

- `protein_feature_select.ipynb`: feature selection, z-scored protein matrix export, Phase 2 CV summary, 32D embeddings, selected protein list
- `phase2_improved_cv.py`: expanded cross-validation utilities for the focused classifier family
- `phase2_full_optimization.py`: long-running two-stage architecture and hyperparameter search
- `phase2_alternative.py`: alternative tree/boosting baseline experiments
- `phase2_push_to_80.py`: aggressive search toward the original 80% target
- `phase2_hyperparameter_study.py`: notebook-template generator for structured hyperparameter studies
- `demo_improvements.py`: quick demo of the improved Phase 2 search space
- `test_architectures.py`: smoke test for architecture construction/import wiring
- `autoregressive_embedding.py`: reusable autoregressive embedding module
- `autoregressive_feature_embedding.ipynb`: temporal embedding training and export notebook
- `expression_comparison.ipynb`: downstream comparison notebook for bundle-generated embeddings and expression baselines
- `results/`: generated CSV, PKL, PNG, TSV, and checkpoint outputs for this bundle

## Inputs And Outputs

Read-only shared inputs:

- `data/cell_lineage.json`
- `data/2023-06-29_entropy_cell_key_V2.csv`
- `data/protein/aggregated_all/s3.csv`
- `data/protein/aggregated_all/lineage_binary_expression.csv`
- `data/protein/aggregated_all/s3_pca_10d.csv`
- `data/protein/aggregated_all/tissue_specific_tf.tsv`
- `data/protein/aggregated_scanpy.h5ad`
- `data/viscello/lin_sc_expr_190602.rds`
- `data/embryo1/tracks.txt`

Bundle-local generated outputs:

- `results/s3_zscore.csv`
- `results/phase2_cv_summary.csv`
- `results/embeddings_32d.csv`
- `results/nn_selected_proteins_rev.tsv`
- `results/ar_embeddings_32d.csv`
- `results/ar_feature_importance.tsv`
- `results/ar_model_checkpoint.pt`
- `results/stage1_results.csv`
- `results/stage1_results.pkl`
- `results/stage2_results.csv`
- `results/stage2_results.pkl`
- `results/phase2_alternative_results.csv`
- `results/phase2_push_results.csv`
- `results/phase2_stage2_results.png`

## Quick Start

From the repository root:

```bash
conda activate dev

python3 expression_embedding/cv_push/test_architectures.py
python3 expression_embedding/cv_push/demo_improvements.py
```

Longer optimization run:

```bash
python3 expression_embedding/cv_push/phase2_full_optimization.py
```

## Notebook Order

Recommended order if you want to regenerate the moved notebook artifacts:

1. `protein_feature_select.ipynb`
2. `autoregressive_feature_embedding.ipynb`
3. `expression_comparison.ipynb`

That order ensures the downstream comparison notebook can see both the protein-feature outputs and the autoregressive outputs in `results/`.

## Path Conventions

Python scripts inside this bundle use `bundle_paths.py` to do two things consistently:

- add the repository root to `sys.path` so shared modules like `feature_selection.py` still import cleanly
- keep generated outputs inside this bundle instead of writing back into the repository root

The moved notebooks now follow the same convention with notebook-local path setup cells.