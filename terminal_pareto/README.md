# Terminal-cell lineage Pareto analysis

This module tests whether terminal-cell parentage in nematode embryogenesis is
Pareto optimal with respect to two competing objectives:

- **travel distance**, derived from embryo cell tracking; and
- **cell-state distance**, derived from protein or RNA expression profiles.

Terminal-cell identities, final positions, and molecular states are held fixed.
Alternative lineages are reconstructed by assigning terminal cells to candidate
parents with minimum-cost bipartite matching across a sweep of objective weights.

The primary publication analysis uses *C. elegans* protein expression plus cell
tracking (`n = 299` terminal cells). Comparative analyses cover *C. elegans* RNA
and *C. briggsae* RNA, tracking replicates, subtrees, two-dimensional tracking,
z-noise controls, edge perturbations, and random-feature controls.

## Coordinate convention

Travel and cell-state costs are reported in standard-deviation units from the
first-cousin-shuffle null distribution. All displayed cost coordinates are then
translated so the natural lineage is at `(0, 0)`. This translation changes only
the presentation origin; it does not change assignments, the Pareto front, cost
rankings, or the relative-Pareto statistic.

## Final null models

The publication figures use progressively less constrained assignment shuffles:

1. first-cousin shuffle (shared grandparent);
2. second-cousin shuffle;
3. third-cousin shuffle; and
4. full-random assignment across candidate parents.

The superseded 2014 *C. briggsae* tracking series (`1407…`) is excluded from the
final replicate analysis.

## Repository structure

```text
terminal_pareto/
├── data_loader.py          # Tracking, expression, lineage, and feature I/O
├── pareto_engine.py        # Assignment, null-model, and perturbation methods
├── lineage_metrics.py      # Edge retention and lineage-tree distance metrics
├── plot_style.py           # Shared publication plotting style and export helpers
├── main.py                 # Marimo notebook for the complete analysis
├── publication_figures.py  # Final main and supporting publication panels
└── output/
    ├── ce_protein/         # Configuration-specific diagnostic figures
    ├── ce_rna/
    ├── cb_rna/
    └── publication/        # Final panels, TeX sources, and assembled PDFs
```

Generated output is ignored by Git. The three final LaTeX source files in
`output/publication/` are intentionally versioned because they contain the
publication captions and the TikZ null-model schematic.

## Running the analysis

From the repository root, with the `dev` environment active:

```bash
python terminal_pareto/main.py
```

This runs the complete Marimo analysis and refreshes the diagnostic figures in
`terminal_pareto/output/`. For an interactive session:

```bash
marimo edit terminal_pareto/main.py
```

The notebook is the analysis pipeline only; it no longer assembles publication
figures.

## Generating publication figures

Run the dedicated presentation pipeline:

```bash
python terminal_pareto/publication_figures.py
```

Supporting panel B is exported under both definitions by default:

- `terminal_support_B_er_td_changed_edges.pdf`; and
- `terminal_support_B_er_td_all_edges.pdf`.

The canonical `terminal_support_B_er_td.pdf` used by the supporting TeX figure
is the changed-edges version. To generate only one definition, run either:

```bash
python terminal_pareto/publication_figures.py --tree-distance-mode changed_edges
python terminal_pareto/publication_figures.py --tree-distance-mode all_edges
```

It writes these canonical panel assets:

- `terminal_main_pareto_front.pdf` — main Pareto-front figure;
- `terminal_support_B_er_td.pdf` — edge retention and tree distance;
- `terminal_support_C_structural_retention.pdf` — structural retention while moving toward
  the two single-objective optima.

The supporting schematic and figure wrappers are:

- `A_null_models_tikz.tex`;
- `ce_terminal_pareto_main.tex`; and
- `ce_terminal_pareto_supporting.tex`.

To assemble the final PDFs with Tectonic:

```bash
cd terminal_pareto/output/publication
tectonic ce_terminal_pareto_main.tex
tectonic ce_terminal_pareto_supporting.tex
```

## Principal metrics

- **Edge retention:** fraction of natural terminal parent–child edges preserved.
- **Mean lineage-tree distance:** mean tree separation between natural and
  reconstructed parent assignments.
- **TWER:** tree-distance-weighted edge retention.
- **Local perturbation tests:** pair and triple edge swaps that test whether both
  objectives can be improved near the natural lineage.

The final supporting figure emphasizes the structural price of pursuing the
last 5% of attainable distance reduction from the maximum-edge-retention
compromise toward either single-objective optimum.

## Dependencies

Python dependencies are provided by `environment.yml` and include Marimo,
NumPy, SciPy, pandas, Matplotlib, seaborn, joblib, and tqdm. Tectonic is used to
compile the standalone publication figures.
