# Track A: expression-feature calibration

**Question.** Is the Pareto geometry (front shape, natural-lineage position u,
distances d_LP and d_CP, maximum edge retention) a property of the chosen
markers or of the cost structure itself? Compare the selected top-k features
with random feature sets of the same size, as the tracking study compared
individual embryos with pooled travel.

## A0. Feature rankings (do first)

The selected sets come from L2 multinomial logistic regression (C = 1,
cell-type merge schema D) on terminal cells; global importance is
max over classes of |coefficient|.

- **RNA** (cross-species, 237 TFs shared by both species):
  `expression_embedding/results/cross_species_rna_linear/linear_coefficients.csv`
  holds coefficients for all 237. Ranking by max |coef| reproduces the saved
  top 20 (`rna_selected_features.tsv`) exactly (checked 2026-10-07). Top-5,
  10 and 30 can be read off directly.
- **Protein** (210 features, named protein_strain, e.g. `CND-1_SYS393`): only
  the top 20 were saved (`expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv`).
  The ranking was computed in `expression_embedding/elegans_protein_linear_baseline.ipynb`
  on `data/protein/aggregated_all/s3.csv` (not the z-scored table the Pareto
  analyses use). Extract the notebook's selection code into a script,
  regenerate the full ranking, confirm it reproduces the saved top 20, and save
  the full ranking beside it. If it does not reproduce, stop and report.

## A1. Feature-set survey

**Feature sets**, for k in {5, 10, 20, 30}:
- `top-k`: the first k of each ranking (k = 20 is the published set).
- `random-all-k`: k features drawn uniformly from the full pool.
- `random-hv-k`: k features drawn from highly variable features only.

**Highly variable pool.** Define it on non-z-scored values (every z-scored
feature has unit variance, so variance on `s3_zscore.csv` is meaningless):
protein on `s3.csv`, RNA on the stored TF values (check their scale first).
Use a dispersion or variance ranking over the cells in the analysed cohort and
choose the pool size after looking at the distribution; it must comfortably
exceed 30 (a starting proposal is the top half: about 105 proteins, 118 TFs).
Record how many top-k features fall in the pool.

**Draws.** Start with 200 random sets per (k, pool); fixed, recorded seeds.
RNA draws are one TF set used for both species (paired comparison).

**Cohorts and configurations** (terminal only at first):
1. Primary C. elegans protein cohort, 275 pooled edges (Figures 2-6).
2. Cross-species matched cohort, 187 edges: CE protein, CE RNA, CB RNA, each in
   3D and XY (Figures 8-9).

Travel matrices do not depend on features: load them from the existing caches
and recompute only the molecular matrix per feature set, with each pipeline's
current metric (cosine for terminal). Reuse
`terminal_pareto.cross_species_analysis.analyze` (matrix-generic: exact
first-cousin scaling, 301 weights, endpoint ties). Before the survey, confirm
that the published top-20 set reproduces the cached published metrics for both
cohorts; the primary cohort's metrics come from
`fig5_table1_ce_canonical_metrics.compute_front_landmarks` and the global
analysis, so check the two implementations agree.

**Record per feature set and configuration:** u, d_LP, d_CP (cache name
`d_NP`), maximum edge retention, retention at P*, whether the natural lineage
is exactly on the front, endpoint spans, natural-lineage costs in
first-cousin-null SD units, and the feature list itself.

**Implementation.** New add-on module (proposed
`terminal_pareto/feature_sensitivity.py`) writing a new run under
`terminal_pareto/output/runs/feature_sensitivity_v1/<run_id>/`, with a cache
identity and source hashes like the other add-ons, and tests. Do not edit
these hash-pinned modules: `cross_species_analysis.py`, `subtree_analysis.py`,
`front_coordinates.py`, `analysis_context.py`, `data_loader.py`,
`lineage_metrics.py`, `pareto_engine.py`, `global_analysis.py`,
`fig5_table1_ce_canonical_metrics.py`, `figS3_cross_geometry.py`, and in
`full_tree_pareto/`: `pooled_analysis.py`, `publication_analysis.py`,
`cousin_references.py`, `cross_species_analysis.py`. If a solver or sampler
is added, add it to the guard in `publication/tests/support.py`.

**Cost.** One sweep is 301 assignments on at most 275 x 275; expect seconds.
About 4 k x 3 set types x 200 draws x 7 configurations, parallelized: hours
locally.

**Sketches** (in `publication/sketches/`, branch `figure-ideas` after merge):
per metric and configuration, the random-all and random-hv distributions
against k with the top-k values marked; the author decides what to present.

## A2. Distance metric (after A1 settles)

The author intended cosine molecular distance everywhere. Terminal analyses use
cosine; the pooled full tree, the partial forest and the per-round references
use Euclidean. Compare both metrics on the terminal cohorts first (cheap),
then decide whether the full-tree analyses must be rerun with cosine. That is a
scientific change: new run identities for the pooled full tree and the partial
forest, and new Figures 7, 10, 11 and S4 (plus the round sketches).
