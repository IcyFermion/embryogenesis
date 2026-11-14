# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment Setup

```bash
conda activate dev
```

All work runs in the `dev` conda environment (defined in `environment.yml`). Key dependencies: numpy, pandas, scipy, scikit-learn, torch, scanpy/anndata, networkx, matplotlib/seaborn, zss (tree edit distance).

## Running Scripts

```bash
# Smoke tests for expression_embedding bundle
python3 expression_embedding/cv_push/test_architectures.py
python3 expression_embedding/cv_push/demo_improvements.py

# Full optimization run
python3 expression_embedding/cv_push/phase2_full_optimization.py

# Pareto optimization on SLURM (pass num_threads as arg)
python3 pareto_slurm.py <num_threads>
```

Notebooks should be run in the order documented in `expression_embedding/README.md`.

## Architecture Overview

This is a computational biology research codebase for analyzing *C. elegans* embryogenesis. The project studies cell lineage trees and transcription factor (TF) protein expression patterns.

### Core Data Model

- **Cell lineage tree**: stored in `data/cell_lineage.json`. Nodes have `did` (lineage ID) and `children`. A handful of nodes use `name` instead of `did` — `utils.lineage_name_mapping()` handles these.
- **TF protein expression**: per-TF CSVs in `data/protein/` (named `<PROTEIN>_<STRAIN>.csv`). Aggregated forms: `data/protein/aggregated_all/s3.csv`, `s3_pca_10d.csv`, `lineage_binary_expression.csv`. AnnData form: `data/protein/aggregated_scanpy.h5ad`.
- **Cell tracking data**: `data/embryo1/tracks.txt` (and embryo2, embryo3). Tab-separated with a `t` (time) column used for temporal cutoffs.

### Key Python Modules

| File | Purpose |
|---|---|
| `utils.py` | `bidict` (bidirectional dict), `lineage_name_mapping()`, `load_json()` |
| `pareto_core.py` | `LineageTree` and `LineageOptimization` — the main tree data structure and optimization engine; uses multiprocessing via `process_map` |
| `lineage_optimization.py` | Duplicate/extracted version of `pareto_core.py` classes, created to fix multiprocessing pickling issues |
| `feature_selection.py` | `SparseGateClassifier` (sigmoid-gated feature selection + embedding), `FocusedClassifier`, CV utilities `cross_validate_features()` / `cross_validate_focused()` / `train_one_pass()` |
| `expression_embedding/autoregressive_embedding.py` | `AutoregressiveGRU` — 2-layer GRU for temporal protein expression; `load_temporal_sequences()`, `train_autoregressive()`, `extract_embeddings()` |
| `kdtree.py` | Vendored kd-tree implementation (from stefankoegl/kdtree) |
| `phylogeny.py` | `Node` / `PhyloTree` — phylogenetic tree dataclasses (separate from lineage tree) |
| `tree_visualizer.py` | Visualization utilities for lineage trees |
| `pareto_slurm.py` | Entry point for SLURM HPC runs of `LineageOptimization` |

### `LineageTree` Internals

Nodes are indexed by a compact integer `tree_id` (assigned in insertion order). Key mappings:
- `lineage_id_mapping[tree_id]` → lineage ID string
- `reverse_lineage_id_mapping[lineage_id]` → tree_id
- `children_list[tree_id]`, `parent_list[tree_id]` (root parent = -1)

### ML Pipeline (expression_embedding bundle)

Two parallel embedding strategies feed into comparison notebooks:

1. **Supervised (SparseGateClassifier)**: Phase 1 selects features via stratified CV over a hyperparameter grid; Phase 2 fine-tunes on the selected features. Outputs: selected protein list + 32D embeddings.
2. **Autoregressive (AutoregressiveGRU)**: Unsupervised; trains on temporal sequences from `aggregated_scanpy.h5ad`; extracts per-cell embeddings via mean/last/weighted pooling + PCA.

Recommended notebook execution order:
1. `expression_embedding/protein_feature_select.ipynb`
2. `expression_embedding/autoregressive_feature_embedding.ipynb`
3. `expression_embedding/timepoint_terminal_feature_select.ipynb`
4. `expression_embedding/expression_comparison.ipynb`

Results land in `expression_embedding/results/`. Scripts inside `expression_embedding/` use `bundle_paths.py` to add the repo root to `sys.path` so shared modules import cleanly.

### Pareto analyses and publication figures

- `terminal_pareto/` and `full_tree_pareto/` are numerical back ends (analyses, caches under `output/runs/*/*/analysis/`, validators). They do not draw figures.
- `publication/` is the front end: notation (`publication/notation.py` is the single source of displayed symbols), figures, captions and Table 1. Read `publication/README.md` before figure or notation work; the project overview is `README.md`.
- **Manuscript figures** are one free-standing set at `publication/output/manuscript/current`, every figure and Table 1 drawn by `publication/v2/` with `publication/palette.py`: `python -m publication.v2 --manuscript --output-dir publication/output/candidates/NAME`, then `python -m publication.v2.manuscript --v2-build ... --output-dir publication/output/manuscript/NAME`, then `--promote` that set (commit first); `--verify` checks the official set. Captions live in `publication/v2/captions.py`. The manuscript text is written in the authors' Overleaf document, outside this repository.
- Figures and analyses kept in reserve are recorded with their regeneration commands in `publication/BACKLOG.md`. The registry's figure families (`python -m publication build --family F --output-dir publication/output/candidates/NAME`) only regenerate backlog figures; releases into a production folder are retired (2026-11-13).
- The research tracks (A-E) are closed; their briefs (`plans/`) and the old handoff documents were retired on 2026-11-13 and remain in git (`git show 0b38600:PATH`). Extend the package READMEs rather than adding new handoff documents.

### Path Conventions

All scripts assume they are run from the **repository root**. Data paths are relative (`data/cell_lineage.json`, etc.). The `expression_embedding/` bundle writes outputs to `expression_embedding/results/` rather than back to root.
