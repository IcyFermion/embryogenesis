# Terminal-Only Pareto Front Analysis

Pareto-optimal assignment of terminal cell identities trading off **spatial proximity**
(L2 distance in µm) vs **expression similarity** (cosine distance on 20 selected proteins/TFs).
All costs are expressed as z-scores relative to a cousin-randomisation null — (0,0) = null
mean, 1 unit = 1σ.

## Structure

```
terminal_pareto/
├── README.md
├── __init__.py
├── data_loader.py       # All I/O: tracking, expression, lineage, feature lists
├── pareto_engine.py     # Pure computation: Pareto sweep, null models, perturbation
├── plot_style.py        # Shared publication palette, typography, and export helper
└── main.py              # Marimo notebook orchestrating the full analysis
```

## Quick start

```bash
conda activate dev
python terminal_pareto/main.py          # run headless (generates all figures/tables)
marimo edit terminal_pareto/main.py     # interactive notebook
```

## `data_loader.py`

Path resolution via `_p(relpath)` — all data paths are relative to the repo root
and resolved at import time from the file's location. Works regardless of cwd.

Key exports:
- `load_elegans_tracking(tcut, path)` / `load_briggsae_tracking(tcut, path)` — 4D embryo tracking
- `load_protein_expression()`, `load_ce_rna()`, `load_cb_rna()` — expression matrices
- `load_prot_sel()`, `load_rna_sel()` — selected feature lists
- `collect_terminals(root, valid_names, subtree)` — DFS lineage traversal
- Constants: `T_CE=255`, `T_CB_NEW=143`, replicate file lists

## `pareto_engine.py`

All functions are pure — input arrays/dicts, output arrays/dicts. No file I/O, no plotting.

| Category | Key functions |
|----------|--------------|
| Pareto | `pareto_sweep`, `lineage_edge_ratio` |
| Cost matrices | `build_cost_matrices` |
| Std normalisation | `compute_std_scaled_pareto`, `compute_cousin_random_stats`, `get_null_cloud`, `lineage_std_position` |
| Null models | `build_random_cost_matrices` (random-20-features), `run_z_noise_pareto` (Gaussian z-axis) |
| Perturbation | `edge_perturbation_choose2`, `edge_perturbation_choose3` |
| Lineage | `build_grandparent_map`, `build_cousin_groups`, `collect_all_subtrees` |
| Cell types | `type_shannon_entropy`, `decompose_by_cell_type` |

## `main.py` — Marimo notebook sections

| Section | Contents |
|---------|----------|
| **Data loading** | Tracking, expression, lineage, cost matrices, random stats |
| **1. Replicate consistency** | 3 C. elegans embryos + 13 C. briggsae imaging sessions |
| **2. Cross-species** | C. elegans protein (T=255) vs C. elegans RNA (T=255) vs C. briggsae RNA (T=143) |
| **2a. 2D XY-plane** | Same comparison with z-axis dropped (controls for anisotropic resolution) |
| **2b. Z-noise null** | XY-calibrated Gaussian noise replacing real z-coordinates |
| **3. Subtree analysis** | 5 major branches (Full, AB, ABa, ABp, P1) × 3 configs |
| **4. Edge perturbation** | C(n,2) and C(n,3) edge swaps testing local optimality |
| **5. Random-20-features null** | 200 random feature sets vs selected 20 (from `terminal_expression.ipynb`) |
| **6. Grand summary** | Combined results table |

## Key results (from `pareto_presentation.ipynb`)

| Config | N | Expr Opt ER | Max ER | Spatial Opt ER |
|--------|---|------------|--------|----------------|
| C. elegans Protein (T=255) | 299 | 0.204 | 0.819 | 0.599 |
| C. elegans RNA (T=255) | 231 | 0.074 | 0.805 | 0.684 |
| C. briggsae RNA (T=143) | 209 | 0.077 | 0.641 | 0.593 |

- **Z-axis**: Brigg. z-signal = +0.376 vs Eleg. = +0.50–0.55 (xy-calibrated noise null)
- **Perturbation**: Save Both% ≈ 0 for all configs — lineage is at a tight local Pareto optimum
- **Random features**: Selected 20 outperform random by ~2× (CE protein) to ~1.7× (CE/RNA)

## Dependencies

marimo ≥ 0.23, numpy, scipy, pandas, matplotlib, seaborn, joblib, tqdm.
All available in the `dev` conda environment (`environment.yml`).

## Figure output

Figures use a colour-vision-deficiency-safe palette, compact journal typography,
light grids, and editable TrueType text in vector files. The notebook writes
300-dpi PNG files to `terminal_pareto/output/`; `plot_style.save_figure` is
available when matching editable PDF output is needed for final assembly.

Cross-species figures remain directly in `output/`. Configuration-specific
figures are grouped under `output/ce_protein/`, `output/ce_rna/`, and
`output/cb_rna/`, including replicate, subtree, and lineage-proximity panels.
