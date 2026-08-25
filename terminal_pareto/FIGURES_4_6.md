# C. elegans Figures 4--6: analysis and production record

Updated 2026-08-25. This document replaces the former exploration and
handoff notes. It records the selected manuscript figures, the reproducible
scientific results behind them, the remaining caveats, and the commands needed
to regenerate the assets. Historical candidate layouts are intentionally not
documented here; Git history retains that development record.

## Status

Figures 4--6, Figures S1--S2, and Table 1 are selected, named consistently,
and compiling as one-page standalone PDFs. The finalized deliverables are:

| Item | Purpose | TeX/PDF stem | LaTeX label |
|---|---|---|---|
| Figure 4 | Investigated-subtree map | `fig4_ce_subtree_map` | `ce_subtree_map` |
| Figure 5 | Canonical subtree Pareto summary | `fig5_ce_canonical_summary` | `ce_canonical_summary` |
| Figure 6 | Cell-type retention and within-type optimization | `fig6_ce_cell_types` | `ce_cell_types` |
| Figure S1 | Cousin-null sensitivity summary | `figS1_ce_canonical_summary_cousin_r` | `ce_canonical_summary_cousin_r` |
| Figure S2 | Cell-type endpoint cost gains | `figS2_ce_cell_type_cost_gain` | `ce_cell_type_cost_gain` |
| Table 1 | Per-subtree statistics | `table1_ce_subtree_statistics` | `ce_subtree_statistics` |

All publication wrappers and their required component assets live in
`terminal_pareto/output/publication/`. Numerical caches live one level above,
in `terminal_pareto/output/`, so the publication directory contains only
manuscript deliverables.

## Analysis configuration

- Species and data: *C. elegans* embryo 1, three-dimensional tracking,
  `T <= 255`, top-20 protein features.
- Scope: terminal cells only.
- Pareto objectives: travel distance and cell-state distance.
- Main null: first-cousin assignment shuffles, seed 42.
- Display coordinates in the terminal pipeline: separate null-standardized
  axes translated so the natural lineage is at `(0, 0)`.
- Subtree inclusion threshold: at least 12 usable biological terminal cells,
  giving 43 investigated subtrees.
- Optimization sweep: 300 weights for publication output.

Null means, variances, and covariance are computed in closed form for every
subtree. Random draws are used only to display null clouds or estimate
optimality frequencies when exact assignment enumeration is infeasible. A
non-positive null standard deviation raises before division.

Cell types come from `data/2023-06-29_entropy_cell_key_V2.csv`, joined on
`wormweb.lineage` and consolidated with the repository's `MERGE_MAP`. The 32
terminal cells absent from that key exactly match `data/apoptotic_cells.txt`
and are typed as programmed death. The two separately annotated dying tail
hypodermal cells, `ABplppppppa` and `ABprppppppa`, remain distinct in the
analysis and enter the displayed `other (n<=4)` group.

## Figure 4: investigated-subtree map

Renderer: `fig4_ce_subtree_map.py`.

The figure places the 43 investigated subtrees in lineage order, with P0 at
the top. Each pie summarizes the fate composition of the exact usable terminal
set used in that subtree's Pareto analysis. Radius scales with `sqrt(n)`.
A solid circumference means the natural assignment occurs on the sampled
Pareto front; a dashed circumference means it does not. Edges connect each
displayed subtree to its nearest qualifying descendant; non-qualifying nodes
are omitted. Subtree names and cell counts use the shared sequential-blue
scale to encode maximum natural-edge retention. The final node layout has no
measured label collisions.

Selected component:
`fig4_ce_subtree_map_panel.{pdf,png}`.

## Figure 5: canonical subtree Pareto summary

Metric builder: `fig5_table1_ce_canonical_metrics.py`.

Renderer: `fig5_figs1_ce_canonical_summary.py`.

Each subtree has its own travel-minimizing assignment `T*`, cell-state-
minimizing assignment `S*`, sampled Pareto front, and natural assignment `L`.
The two objectives are normalized separately using only the endpoint spans:

\[
D_1(a)=\frac{T(a)-T(T^*)}{T(S^*)-T(T^*)},\qquad
D_2(a)=\frac{S(a)-S(S^*)}{S(T^*)-S(S^*)}.
\]

The null-independent natural-lineage distance is

\[
d_{LP}(L)=\min_{p\in\mathcal P_{\mathrm{sampled}}}
\left\|D(L)-D(p)\right\|_2.
\]

The minimum is evaluated only over actual sampled Pareto assignments, never
an interpolated segment. Coordinates are not clipped to the endpoint box;
the natural lineage is not guaranteed to lie inside it. The canonical
position `u_L` is normalized arc length along the front, from travel optimum
`u=0` to cell-state optimum `u=1`, evaluated at the same nearest assignment.

Figure 5 contains:

- a symbolic definition of `u`, `d_LP`, and `d_NP`;
- paired `d_LP` and `d_NP` values for P0, AB, ABa, ABp, and P1; and
- all 43 eligible subtrees in `(u, d_LP)` space, with the five major-subtree
  `d_NP` values overlaid as references.

Marker shape and direct labels identify subtrees. Neutral gray encodes `u`,
purple encodes `d_LP`, and gold encodes the first-cousin null mean and
`d_NP`. Maximum edge retention is deliberately omitted from Figure 5 so each
marker communicates canonical position and proximity only; it remains in
Figures 2, 3, 4, and 6 and in Table 1.

The older cousin-relative distance `r` remains in Table 1 and appears in
Figure S1 as a sensitivity analysis. Across the 43 subtrees, `d_LP` and `r`
have Spearman correlation 0.95, indicating that the broad proximity ranking
does not depend on the first-cousin normalization.

Table 1 also reports the endpoint-normalized distance from the first-cousin
null mean to the same sampled assignment used for `d_LP`:

\[
d_{NP}=\left\|D(N)-D(P^*)\right\|_2.
\]

This `d_NP` shares the Figure 5 geometry with `d_LP`. The supplementary `r`
retains its original first-cousin-null-standardized geometry, so it must not
be interpreted as `d_LP / d_NP`.

Selected components:

- `fig5A_ce_canonical_definition.{pdf,png}`;
- `fig5B_ce_canonical_major_subtrees.{pdf,png}`;
- `fig5C_ce_canonical_all_subtrees.{pdf,png}`;
- `figS1A_ce_cousin_r_definition.{pdf,png}`;
- `figS1B_ce_cousin_r_major_subtrees.{pdf,png}`;
- `figS1C_ce_cousin_r_all_subtrees.{pdf,png}`; and
- `table1_ce_subtree_statistics.tex`.

The manuscript wrappers assemble these component PDFs directly. Superseded
combined canvases and canonical-map prototypes are no longer generated.

## Figure 6: terminal cell types

Renderers:

- `fig6a_figs2_ce_cell_types.py` for Figure 6A and Figure S2;
- `fig6bc_ce_within_type.py` for Figure 6B--C.

Panel A shows natural-edge retention for each terminal fate across all 145
distinct assignments observed along the global Pareto sweep. Biological
parent identity defines retention, so swapping two assignment slots belonging
to the same parent retains both natural edges. Weighted cell-type retention is
asserted to recompose the global edge-retention curve exactly at every cached
assignment; the maximum-retention assignment has global retention 0.819398.
Irregular front positions are resampled to a uniform canonical-position grid
and narrowly averaged for display only. Exact assignments, keypoints, and
statistics remain unsmoothed.

Panels B--C repeat the Pareto optimization after forbidding assignments across
terminal cell-type groups. Six sufficiently large types receive individual
fronts in shared embryo-wide first-cousin-null standard-deviation units per
cell. This null supplies a common scale only; the Pareto assignments are
constrained by cell type, not cousin relationships. The aggregate
type-preserving front is then compared with the unrestricted front. At the
corresponding single-objective endpoints, type restriction retains 3.101
standard deviations more travel distance and 15.073 standard deviations more
cell-state distance. It realizes 16.8% of the unrestricted travel saving and
63.0% of the unrestricted cell-state saving.

Figure S2 preserves the complementary endpoint ledger: signed per-cell travel
and cell-state changes at the travel optimum, maximum-retention compromise,
and cell-state optimum. The two objective changes remain on separate axes and
are never combined linearly.

The shared semantic colors are retained here: teal marks travel optima,
vermillion marks cell-state optima, and blue marks edge retention or the
maximum-retention assignment. Cell-type curve colors are local categorical
encodings identified by their panel titles. The aggregate comparison uses
charcoal for unrestricted and orange for type-restricted assignments.

Selected components:

- `fig6A_ce_retention_heatmap.{pdf,png}`;
- `fig6B_ce_within_type_fronts.{pdf,png}`;
- `fig6C_ce_type_restricted_aggregate.{pdf,png}`; and
- `figS2_ce_cell_type_cost_gain_panel.{pdf,png}`.

The manuscript wrapper assembles panels B and C directly; the former combined
B--C canvas is no longer generated.

Analysis caches:

- `global_front_retention_by_cell_type.csv`;
- `global_front_by_cell_type.csv`;
- `within_type_fronts.csv`;
- `within_type_summary.csv`; and
- `type_preserving_aggregate_front.csv`.

## Reproducible biological findings

1. The natural terminal assignment is exactly sampled-Pareto-optimal in 16
   of 43 investigated subtrees, all with at most 35 terminal cells. Exact
   recovery becomes less frequent among larger subtrees, but the continuous
   size relationships are modest: size versus on-front status has Spearman
   rho `-0.35`, versus front proximity `+0.30`, versus cousin-relative
   distance `+0.18`, and versus maximum retention `-0.06`. This supports a
   frequency statement, not a smooth inefficiency-with-size gradient.
2. Exact recovery is not explained by small assignment spaces. Full-random
   assignments are on the front in 0/1000 draws for every subtree. The
   cousin-shuffle fraction is usually zero and reaches 33% only in the small,
   exactly enumerable MSpa assignment space.
3. Apparent relationships between attainable reductions and subtree depth or
   fate diversity disappear after per-cell normalization (`|rho| < 0.2`).
   Subtree size and median depth are strongly confounded (`rho = 0.90`).
4. True-sibling contrasts are at chance level: 8 of 14 informative pairs
   agree with the proposed direction, with three distance ties excluded.
5. Per-cell endpoint contributions differ by fate. Programmed-death cells
   have the largest cell-state improvement at the cell-state optimum
   (`-0.25 sigma/cell`); epithelium pays the largest travel price there
   (`+0.47 sigma/cell`); muscle pays the largest cell-state price at the
   travel optimum (`+0.22 sigma/cell`). Neurons dominate totals primarily
   because they are abundant.

## Interpretation and statistical caveats

- Nested subtrees are descriptive observations, not independent replicates.
  EMS/MS, ABplpp/ABplppp, and ABprpp/ABprppp can share identical cousin-group
  terminal sets and therefore identical null moments.
- Do not report regression p-values that treat the 43 subtree points as
  independent.
- `d_LP` is a geometric display statistic after endpoint normalization; it is
  not a biological weighted sum of travel and cell-state distances.
- Cousin-relative `r` is a two-dimensional Euclidean distance in separately
  null-standardized axes, not a conventional one-dimensional z-score.
- Equal weighting of standardized objectives must not be described as a
  biologically balanced cost because the objectives are not linearly
  comparable.
- Within-type cell-state gains are not uniformly negligible. Attributing them
  to expression noise would require an independent replicate or stability
  analysis.

## Regeneration workflow

Run from the repository root with the `dev` environment:

```bash
conda run -n dev python3 terminal_pareto/subtree_analysis.py --min-cells 12
conda run -n dev python3 terminal_pareto/fig4_ce_subtree_map.py --min-cells 12
conda run -n dev python3 terminal_pareto/fig5_table1_ce_canonical_metrics.py --min-cells 12
conda run -n dev python3 terminal_pareto/fig5_figs1_ce_canonical_summary.py
conda run -n dev python3 terminal_pareto/fig6a_figs2_ce_cell_types.py --iteration 300
conda run -n dev python3 terminal_pareto/fig6bc_ce_within_type.py --iteration 300
```

Compile the wrappers from the publication directory:

```bash
cd terminal_pareto/output/publication
conda run -n dev tectonic fig4_ce_subtree_map.tex
conda run -n dev tectonic fig5_ce_canonical_summary.tex
conda run -n dev tectonic fig6_ce_cell_types.tex
conda run -n dev tectonic figS1_ce_canonical_summary_cousin_r.tex
conda run -n dev tectonic figS2_ce_cell_type_cost_gain.tex
conda run -n dev tectonic table1_ce_subtree_statistics.tex
```

Run the canonical-metric validation checklist with:

```bash
conda run -n dev python3 terminal_pareto/fig5_table1_ce_canonical_metrics.py \
  --min-cells 12 --selftest
```

The self-test checks endpoint mapping, duplicate-solution invariance, sampled
front proximity, agreement with the subtree summary, on-front classification,
maximum-retention selection, threshold invariance, and output geometry.

## Supporting analysis code

- `subtree_analysis.py`: reusable per-subtree summary tables, exact null
  moments, and exact/Monte Carlo optimality frequencies with Wilson intervals.
- `subtree_explore.py`: global-front cell-type decomposition helpers used by
  Figure 6A and Figure S2 (despite the historical filename, it no longer
  contains exploratory plotting layouts).
- `ce_subtree_canonical_metrics.csv`: authoritative 43-row canonical metric
  cache used by Figure 5, Figure S1, and Table 1.
- `ce_subtree_canonical_curves.npz`: canonical curve cache used by validation.

## Output policy and remaining work

Only manuscript wrappers, their directly included panel assets, and the
authoritative numerical caches listed above belong in the active output loop.
Exploratory subtree scatters/fronts, alternative minimum-cell threshold CSVs,
combined pre-split canvases, and prototype canonical maps have been removed.
They remain recoverable from Git history or can be recomputed from the analysis
helpers if a new question requires them.

The figures are release candidates. Remaining work is editorial integration:
copy or include the finalized wrappers in the manuscript build, reconcile the
older planning outline in `task.md` with the accepted numbering, and ensure all
figure and table labels resolve in the complete paper. `paper/paperv2.tex` is a
reference manuscript and is not modified by this pipeline.
