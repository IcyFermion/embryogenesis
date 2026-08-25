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

## Publication color semantics

Recurring quantitative concepts use one shared visual vocabulary across the
manuscript figures:

- sequential blue: natural-edge retention;
- teal: travel direction or travel-minimizing assignment;
- vermillion: cell-state direction or cell-state-minimizing assignment;
- neutral gray: canonical front position `u` or tree geometry;
- purple: natural-lineage-to-front distance `d_LP`; and
- gold: the first-cousin null and its front distance `d_NP`.

Maximum-retention assignments use a distinctive marker whose fill remains on
the blue retention scale. Cell types and species retain local categorical
palettes with explicit legends.

## Repository structure

```text
terminal_pareto/
├── data_loader.py          # Tracking, expression, lineage, and feature I/O
├── pareto_engine.py        # Assignment, null-model, and perturbation methods
├── lineage_metrics.py      # Edge retention and lineage-tree distance metrics
├── plot_style.py           # Shared publication plotting style and export helpers
├── main.py                 # Marimo notebook for the complete analysis
├── fig2_fig3_ce_terminal_pareto.py       # Figures 2--3
├── fig4_ce_subtree_map.py                # Figure 4
├── fig5_table1_ce_canonical_metrics.py   # Figure 5 metrics + Table 1
├── fig5_figs1_ce_canonical_summary.py    # Figure 5 + Figure S1
├── fig6a_figs2_ce_cell_types.py          # Figure 6A + Figure S2
├── fig6bc_ce_within_type.py              # Figure 6B--C
├── FIGURES_4_6.md                        # Current science and workflow notes
└── output/
    ├── ce_protein/         # Configuration-specific diagnostic figures
    ├── ce_rna/
    ├── cb_rna/
    └── publication/        # Final panels, TeX sources, and assembled PDFs
```

Generated output is ignored by Git. The final LaTeX wrappers and TikZ
null-model schematic in `output/publication/` are intentionally versioned
because they contain publication captions and layout.

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
python terminal_pareto/fig2_fig3_ce_terminal_pareto.py
```

Supporting panel B uses the changed-edges tree-distance definition. The
all-edges alternative remains selectable but is written to the diagnostic
`output/ce_protein/` directory rather than the publication directory:

```bash
python terminal_pareto/fig2_fig3_ce_terminal_pareto.py --tree-distance-mode changed_edges
python terminal_pareto/fig2_fig3_ce_terminal_pareto.py --tree-distance-mode all_edges
```

It writes these canonical panel assets:

- `fig2_ce_terminal_pareto_front.pdf` — main Pareto-front panel;
- `fig3B_ce_edge_retention_tree_distance.pdf` — edge retention and tree distance;
- `fig3C_ce_structural_retention.pdf` — structural retention while moving toward
  the two single-objective optima.

The supporting schematic and figure wrappers are:

- `fig3A_ce_null_models.tex`;
- `fig2_ce_terminal_pareto_main.tex`; and
- `fig3_ce_terminal_pareto_supporting.tex`.

To assemble the final PDFs with Tectonic:

```bash
cd terminal_pareto/output/publication
tectonic fig2_ce_terminal_pareto_main.tex
tectonic fig3_ce_terminal_pareto_supporting.tex
```

The finalized C. elegans manuscript outputs are:

| Manuscript item | TeX/PDF stem | LaTeX label |
|---|---|---|
| Figure 2 | `fig2_ce_terminal_pareto_main` | `ce_terminal_pareto_main` |
| Figure 3 | `fig3_ce_terminal_pareto_supporting` | `ce_terminal_pareto_supporting` |
| Figure 4 | `fig4_ce_subtree_map` | `ce_subtree_map` |
| Figure 5 | `fig5_ce_canonical_summary` | `ce_canonical_summary` |
| Figure 6 | `fig6_ce_cell_types` | `ce_cell_types` |
| Figure S1 | `figS1_ce_canonical_summary_cousin_r` | `ce_canonical_summary_cousin_r` |
| Figure S2 | `figS2_ce_cell_type_cost_gain` | `ce_cell_type_cost_gain` |
| Table 1 | `table1_ce_subtree_statistics` | `ce_subtree_statistics` |

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
