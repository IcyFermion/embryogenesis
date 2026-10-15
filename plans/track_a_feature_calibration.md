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
- **Done 2026-10-07.** `expression_embedding/elegans_protein_feature_ranking.py`
  reproduces the saved top 20 in order (importances match the notebook to 4
  decimals) and wrote `protein_feature_ranking.csv` (all 210, with per-class
  coefficients) beside it. Importance gaps at the k / k+1 boundary are small:
  protein 0.018, 0.012, 0.009, 0.012 and RNA 0.006, 0.004, 0.002, 0.0005 for
  k = 5, 10, 20, 30, so membership near each cutoff is fragile.

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

**Decisions while implementing (2026-10-07).**
- *C. elegans only for RNA* (author): C. briggsae tracking has larger z-axis
  uncertainty, so CB RNA is dropped. Without the CB match the RNA cohort is
  the **217** RNA-covered edges of the 275 (up from 187); no tracked RNA edge
  lies outside the 275, so 217 is the maximum. Survey configurations: primary
  protein (275, 3D) and CE protein / CE RNA on the 217 edges in 3D and XY.
  3D travel is the published pooled matrix restricted with fixed 275-edge
  denominators; XY is recomputed as in `cross_species_analysis` and equals the
  cached 187-edge XY matrix on those edges. The 187-edge caches serve only
  the top-20 reproduction check.
- *Highly variable pool*: variance of log1p values over cohort cells (raw
  variance tracks the mean, Spearman 0.96), top half: 105 proteins on `s3.csv`
  over the 275-edge cells, 118 TFs on CE RNA over the 217-edge cells. A quarter
  of the proteins are expressed in at most about 1% of cohort cells.
- *Zero vectors*: cosine is undefined for an all-zero cell; RNA draws that
  leave one (about 20% of random k = 5 sets) are redrawn and the number of
  rejections recorded. Protein z-scores are never zero.
- *Reproduction*: the published top-20 sets reproduce the primary P0 row and
  all four CE configurations of the 187-edge cache to 4e-15, and the
  cross-species `analyze` agrees with the Figure 5 landmark code on P0.
- *Worktrees*: the context cache key embeds the absolute tracking path, so
  `validate_existing_context_manifest` (and the publication terminal-primary
  tests) fail outside the main checkout. The survey checks matrix and source
  hashes instead. `data/embryo*` must be real copies, not links.

**Status 2026-10-07.** Survey run
`terminal_pareto/output/runs/feature_sensitivity_v1/survey_20261007/` (deleted 2026-10-14, superseded)
(analysis id `883aa81b…`): 3,208 feature sets, 8,020 sweeps, 4.4 min on 24
workers; validated (redraw, summary, 45 replays, 26 source hashes).
`analysis/summary.csv` places each top-k value in both random distributions.
11 random-all protein sweeps (k = 5, 10) are `degenerate_state`: mostly silent
panels give no first-cousin variance in cell-state cost. They are counted, not
summarized. Next: sketches, then A2.

**Sketches** (in `publication/sketches/`, branch `figure-ideas` after merge):
per metric and configuration, the random-all and random-hv distributions
against k with the top-k values marked; the author decides what to present.

**Bottom-k controls (2026-10-14).** Run `survey_20261014` (analysis id
`a93b23df…`) adds two fixed sets per k: the classifier's bottom k, and the
bottom k within the HV pool. Random draws are identical to `survey_20261007`,
which no longer validates against the current code and can be deleted. RNA
controls behave as negative controls; protein ones do not. The bottom-ranked
proteins are nearly silent, while the bottom HV proteins still carry cell-type
information. Max |coef| over 210 correlated proteins tracks single-protein
cell-type information weakly. CE RNA gives only 71 distinct profiles for the
217 cells, one per cell type, whereas protein gives one per cell.

## A1b. Choosing the cell-state representation (decision rule, 2026-10-14)

**Goal.** A simple, defensible choice of cell-state representation for the
Pareto analyses. **The front is the outcome, never a selection input.**
Choosing features by how close the natural lineage sits to the front would be
circular. The survey (A1) serves only as a robustness check.

**Roles of the modalities (author).**
- *Protein* is resolved per cell and is the C. elegans gold standard for all
  terminal, subtree and cell-type analyses.
- *RNA* is the only data shared by C. elegans and C. briggsae. Its role is
  limited to the basic terminal-only cross-species comparison (187-edge
  cohort, Figures 8-9). It is resolved by cell type, so it is not used for
  subtree or cell-type analyses. Its representation must be a single TF set
  used for both species.

**Candidates, per modality.** Classifier top-20 (published), the whole HV pool
(no classifier, no k), and the whole panel. All use cosine distance; A2 adds
Euclidean as a second axis.

**Primary criterion: fate discrimination.** Sibling separation, Cohen's d of
cosine distance between terminal sisters of different type and sisters of the
same type. It is reported for Schema D types and for fine types
(`wormweb.desc`), with bootstrap intervals over sister pairs. It is
independent of tracking and of the front.

**Validity check: lineage memory.** Among cells of the same fine type, the
Spearman correlation between expression distance and lineage path length. The
null distribution permutes profiles within each fine type. A positive excess
means the cost carries lineage history, which flatters the natural lineage by
construction. Such a representation is disqualified, or the effect is
reported.

**First results (2026-10-14, scratch, to be recorded as a run).**
- *Cohorts.* Protein on the primary 275 edges (71 fine types; 58/45
  same/different-type sisters, Schema D). RNA on the cross-species 187 edges
  (51 fine types; 46/22 sisters). The RNA HV pool is ranked over CE and CB
  jointly (118 TFs).
- *Fate discrimination* (Schema D / fine types; 95% bootstrap intervals are
  wide, about ±0.4):
  - Protein: top-20 0.73 / 0.70, HV pool 0.43 / 0.58, whole panel 0.16 / 0.35.
  - CE RNA: top-20 1.93 / 2.05, HV pool 1.19 / 2.06, whole panel 1.14 / 2.12.
  - The top-20 lead on Schema D is partly in-sample, since the classifier was
    trained on those labels. On fine types, protein top-20 and the HV pool
    overlap, and the three RNA candidates are equal.
- *Lineage memory* is significant for every candidate in both modalities and
  both species. The excess Spearman correlation over the within-type
  permutation null ranges from 0.19 to 0.35, with every p at the 1/1001
  minimum.
  - Protein: top-20 carries the least (0.195), the HV pool (0.345) and the
    whole panel (0.315) more.
  - It is a property of the data, not of one candidate, so it cannot simply
    disqualify a representation.
- *Fate-only control.* Each child is represented by its fine-type mean
  profile; parents keep their measured profiles. This removes within-type,
  lineage-correlated state while keeping fate differences.
  - Protein, natural-lineage state z (measured → fate-only): top-20
    −4.8 → −3.7, HV pool −7.4 → −4.1, whole panel −7.1 → −4.1.
  - Protein d_SD (measured → fate-only): 1.12 → 1.39, 0.23 → 1.04,
    0.66 → 1.07.
  - Under the fate-only cost the three protein candidates converge, and the
    lineage remains well below the cousin null.
  - RNA fronts barely change in either species, since RNA is already
    resolved by type.
  - The earlier "bottom HV proteins do well" result is lineage memory: broad
    panels carry more inherited signal.

**Recorded (2026-10-14).** `terminal_pareto/representation_criteria.py`, run
`representation_criteria_v1/criteria_metric_20261014` (analysis id `55be99f6…`; it supersedes `criteria_20261014`, now deleted),
validated by a full replay. It reproduces the scratch results above exactly;
only the bootstrap and null intervals moved with the seeds.
- *Decision:* keep the published top-20 for both modalities. Report the
  fate-only front as the control for lineage memory.
- *Methods and results:* `terminal_pareto/METHODS_cell_state_representation.md`,
  the basis for the methodology subsection.
- *Display:* d_SD needs a notation concept before any figure uses it.
- *Outside the pipelines:* a standalone review page lives in
  `artifacts/representation_audit/` (untracked).

**Then.** The choice is unchanged (the published top-20), so the published
fronts and figures stand. Still open: whether and how to show the fate-only
control and the candidate comparison (author), and A2 (distance metric).

## A2. Distance metric (after A1 settles)

**Terminal comparison (2026-10-14).** It is the metric axis of
`representation_criteria.py` (run `criteria_metric_20261014`); results are in
the methods doc. On z-scored protein, cosine and Euclidean fronts agree. On
stored RNA values, Euclidean nearly removes the lineage's cell-state advantage
(C. elegans 3D state z −0.25 against −3.14 for cosine). **Decision (author):**
cosine everywhere. The full tree and the partial forest are switched under new
run identities. A read-only migration map (2026-10-14):
- Switch points: `pooled_analysis.py:148` and
  `full_tree_pareto/cross_species_analysis.py:208`.
- The pooled Gaussian reference scores simulated states with Euclidean edge
  norms and needs a cosine scoring path.
- The partial-forest validator re-derives the state as Euclidean.
- Both modules are hash-pinned by every full-tree cache, so existing runs stop
  validating from this branch. Production keeps verifying, because it checks
  only file hashes.
- No full-tree cohort has a zero expression vector.
- The single-embryo-1 sensitivities (clock, leaf-conditioned, branch
  variability) are Euclidean by construction and outside the publication
  families.

**Status 2026-10-15.** Both full-tree analyses now use cosine:
- *Code:* `pooled_analysis.py` and `full_tree_pareto/cross_species_analysis.py`
  (versions 2), with cosine edge scoring for the Gaussian reference's protein
  block. `resume_paired.py` now creates `analysis/` in a fresh run.
- *Runs:* `pooled_full_tree_v1/terminal_clamped_cosine_20261014` (paired
  matching took 12 min on 24 workers) and
  `cross_species_layerwise_v1/terminal_anchored_cosine_20261014`, both
  replay-validated. Round references are rebuilt for both. The 50 full-tree
  tests pass.
- *Candidates:* `publication/output/candidates/full_tree_cosine_20261014`
  (Figures 7 and S4) and `.../full_tree_cross_species_cosine_20261014`
  (Figures 10 and 11).
- *Not released.* Production still holds the Euclidean figures from
  2026-10-04, so the Figure 7A pixel-parity test fails until a release (the
  author's call).
- *Downstream:* the round sketches and the `figures-v2` drafts in the main
  checkout read the round references, so they need rebuilding after merge.
- *Results:* in the methods doc (A2 section).

The author intended cosine molecular distance everywhere. Terminal analyses use
cosine; the pooled full tree, the partial forest and the per-round references
use Euclidean. Compare both metrics on the terminal cohorts first (cheap),
then decide whether the full-tree analyses must be rerun with cosine. That is a
scientific change: new run identities for the pooled full tree and the partial
forest, and new Figures 7, 10, 11 and S4 (plus the round sketches).
