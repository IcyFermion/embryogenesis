# Terminal-cell Pareto analysis

**Current release, 2026-09-23:** C. elegans protein analysis using pooled travel
from three embryos and endpoint-normalized Pareto displays. This replaces the
single-embryo publication set. [HANDOFF.md](HANDOFF.md) records current status,
validation, remaining manuscript work, and the short development history.
These are the module's two maintained Markdown documents. The resumed
cross-species terminal comparison is described below and uses its own run.

## Start here

- **Manuscript figures:** `output/publication/` (frozen published copy).
- **Published identity:** `output/publication/release_manifest.json` (profile,
  run, rendered-file hashes, and analysis-cache hashes).
- **Working run:** `output/runs/pooled_tracking_v1/migration_candidate_20260920/`.
  The historical run ID is retained to preserve cache identity; it is now the
  source of the release, not an unapproved scientific alternative.
- **Previous figures:** `output/legacy/embryo1/publication/`.
- **Older diagnostic plots:** `output/legacy/diagnostics/`.
- **Archived prototypes/pilots:** compact source and output archives described
  in the handoff. They are outside the active build.

Run commands from the repository root in the `dev` Conda environment.

```bash
# Verify the published bundle without rebuilding.
python terminal_pareto/publication_release.py --verify

# Rebuild layout from validated caches in the working run.
python terminal_pareto/publication_build.py --layout-only

# Validate, compile, then publish a revised layout; preserve the previous release.
python terminal_pareto/publication_build.py --layout-only --publish

# Full pooled pipeline (includes downstream analysis/rendering).
python terminal_pareto/publication_build.py

# Validate scientific results, source identity, and the legacy regression baseline.
python terminal_pareto/validate_pooled_migration.py
python terminal_pareto/test_analysis_context.py
python -m unittest terminal_pareto.test_publication_release
```

A build writes into its run directory. Publication is explicit via `--publish`
or `publication_release.py`; it validates first, checks compiled one-page
wrappers, stages a copy, preserves the previous publication, and writes the
release manifest. Editing working figures does not silently change the release.
Generated PDFs, arrays, and local output archives are Git-ignored and need the
project's external backup. Versioned wrapper templates and source archives do
not replace those data backups.

## Analysis conventions

Terminal identities, measured positions, molecular states, and parent-slot
multiplicities are fixed. Hungarian bipartite matching assigns children to
parent slots across 301 weights (300 intervals including both endpoints).
These are sampled optimal assignments, not enumeration of the complete discrete
Pareto set. Duplicate slots of one biological parent are interchangeable;
edge retention/agreement compares biological parent identity.

The primary cohort is the strict intersection of **275 natural terminal
parent–child edges** across embryo 1 (cutoff 255), embryo 2 (247), and embryo 3
(225). Both endpoints require tracking and the fixed protein representation.
Cutoffs are predefined stage matches, not equal elapsed developmental times.
At least 12 usable terminal cells gives **42 qualifying subtrees**.

Cell-state cost is cosine dissimilarity of the selected top-20 z-scored protein
features from `data/protein/aggregated_all/s3_zscore.csv`, with selection in
`expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv`.
Travel is Euclidean displacement between the retained parent/child positions,
not integrated trajectory length or measured metabolic energy.

For embryo r, let T_r(L) be the natural lineage's total travel on the global
matched cohort. Pooled edge cost is

\[
\bar d_{ij}=\frac{1}{3}\sum_{r=1}^{3}\frac{d^{(r)}_{ij}}{T_r(L)}.
\]

The same global denominators are used for every subtree and null assignment.
Pooled natural total is one. Coordinates are not averaged into a synthetic
embryo. Pooling remains additive over edges and optimizes the equal-weight
mean of the three observed normalized geometries.

Profiles remain explicit: `pooled_tracking_v1` (primary, 275 edges/42 subtrees),
`embryo1_matched` (coverage control, 275/42), and `embryo1_legacy` (299/43).
Low-level loaders retain historical defaults for reproducibility. The main
publication command defaults to the pooled profile. Historical no-profile
figure commands write to `output/legacy/rebuild/publication/`, not the current
release or the frozen old figures.

First-cousin assignment shuffling is the primary null, with seed 42 and 1,000
visualization draws. Figures 2–3 also use second-/third-cousin and full-random
shuffles. Pooled/matched subtree means, variances, and covariance use exact
first-cousin moments. Nonpositive null SDs and invalid endpoint spans are
handled explicitly. Optimizer scaling and display normalization are separate.

## Display coordinates and canonical metrics

Let A minimize travel and B minimize cell-state cost. For any assignment a,

\[
D_1(a)=\frac{T(a)-T(A)}{T(B)-T(A)},\qquad
D_2(a)=\frac{S(a)-S(B)}{S(A)-S(B)}.
\]

A maps to (0, 1), B to (1, 0). Natural lineage and null samples use the same
anchors and can fall outside the endpoint box; they are not clipped. Competing
assignment sets in one panel share that panel's reference endpoints. Separate
subtree/type endpoint scales support shape comparisons, not absolute cost
comparisons between groups. Percentage and null-SD views remain diagnostics;
Figure S1 and the Figure S2 ledger retain their specifically defined scales.

P* is the closest **attained sampled assignment** to natural lineage L in
endpoint coordinates. Canonical distance d_LP is its Euclidean distance to L;
d_NP measures from first-cousin-null mean N to that same P*. Position u is
normalized front arc length from A to P*. No interpolated lineage is used.
The maximum-edge-retention assignment highlighted in Figure 2 need not be P*.
The cousin-relative statistic r uses separately null-standardized geometry;
it is a sensitivity statistic, and d_NP multiplied by r need not equal d_LP.

## Published figures

All wrapper stems below are in `output/publication/` as PDF and TeX unless
specified otherwise. Panel PDF/PNG assets, S3 SVG and projection data, and
provenance accompany the wrappers.

| Item | Wrapper stem | Content |
|---|---|---|
| Figure 1 amendment | `fig1_ce_endpoint_normalization_amendment` | Endpoint/canonical schematic; co-author integration pending |
| Figure 2 | `fig2_ce_terminal_pareto_main` | Single pooled front, null clouds, natural lineage and maximum retention |
| Figure 3 | `fig3_ce_terminal_pareto_supporting` | Null schematic, changed-edge tree distance and structural trade-off |
| Figure 4 | `fig4_ce_subtree_map` | 42 subtrees, fate composition, on-front status and retention |
| Figure 5 | `fig5_ce_canonical_summary` | Major-subtree and all-subtree canonical summaries, A–B |
| Figure 6 | `fig6_ce_cell_types` | Fate retention, within-type fronts and restricted aggregate |
| Figure S1 | `figS1_ce_canonical_summary_cousin_r` | Cousin-relative sensitivity |
| Figure S2 | `figS2_ce_cell_type_cost_gain` | Separate per-cell cost changes at three selected assignments |
| Figure S3 | `figS3_ce_tracking_robustness` | Four travel geometries, full fronts with near-natural insets |
| Table 1 | `table1_ce_subtree_statistics.tex` | Canonical and null-relative metrics |

Figure 5 separates two questions. Panel A pairs natural-lineage and
first-cousin-null-mean distances to each major subtree's closest assignment
P*, using named rows and purple/gold points. Panel B shows all 42 subtrees
once in branch-grouped rows, with separate columns for canonical position
u (gray) and natural-lineage distance d_LP (purple). Column scales are shared
between groups; bold names link the five major subtrees across panels.
Branch names replace the former shape key, and every subtree is directly
labelled. The layout-only build regenerates these panels from the validated
canonical table without changing analysis caches or Figure S1.

S3 measures travel in embryos 1–3 (A–C) and the pooled geometry (D). Color
identifies the source of optimization, solid lines the native front, dashed
lines transferred assignments, and X natural lineage. Every panel evaluates
all four saved 301-solution sweeps unchanged; its inset uses the same
coordinates. The 4,816 projections are replay-validated. Overall front shapes
are similar; near natural lineage, pooled assignments form a compromise closer
to each native front than assignments transferred from another embryo.
The pooled optimum's advantage in its own objective alone is not independent
support for pooling.

## Cross-species terminal comparison (2026-09-30)

`cross_species_analysis.py` builds CE protein, CE RNA and CB AF16 RNA on the
same **187 terminal parent--child edges**, with identical biological-parent
capacities and endpoint displays. The six-panel figure follows Figure 2:
blue retention colors, natural-lineage crosses, outlined sampled-retention
maxima, three cousin-shuffle clouds, and full-random insets. Main axes and the
0--1 retention scale are shared across panels. Its top row is existing 3D
(A--C), its bottom row is 2D (XY; D--F), and columns are CE protein, CE RNA
and CB RNA. Natural lineage and maximum retention are marked separately,
without a connecting segment, as in Figure 2. Figure 9A--B compares front
shapes in paired 3D/XY panels with matching main axes, near-natural zoom
limits, species colors and line styles. Hollow diamonds mark the closest
attained point P*, and dashed segments connect natural lineage to P*, never
to maximum retention. Figure 9C follows the Figure 5B named-row style: three
configurations grouped under each geometry, with separate columns for
canonical position u, natural distance d_LP and analytic first-cousin-null-mean
distance d_NP. Both distances use the same P*. Each metric column has one
scale across all six rows. Each configuration retains its own endpoint anchors.

Working comparison:
`output/runs/cross_species_terminal_v1/pooled_comparison_20260930/`.
The author assigned the six-panel comparison to **Figure 8** and the paired
front overlays plus canonical-metric rows to **Figure 9**. Captioned standalone pages are in this run's
`publication/`; their sources and the September pilot checkpoint remain
independent of the existing 275-edge publication release.

- CE uses the published three-embryo pooled matrix restricted from 275 to
  187 edges, retaining the original 275-edge normalization totals. Protein
  and RNA use identical travel matrices. The CE RNA-covered three-embryo
  cohort has 217 edges before CB matching.
- CB pools the two AF16 embryos at cutoffs 148/156. Fixed denominators are
  their natural totals on the 196-edge shared tracking/RNA reference before
  CE matching. Pooling averages normalized candidate-edge distances.
- The frozen protein panel and shared RNA 20-TF panel use the existing
  representations and cosine distance. Annotation `Missing` flags remain
  descriptive; their measurement meaning is unresolved.
- Exact first-cousin moments are computed on each pooled matrix, preserving
  dependence induced by applying the same assignment across embryos. The
  optimizer uses those SDs and 301 weights. Endpoints use the original pilot's
  checked numerical secondary-objective tie handling. A nested 1,201-weight
  check is saved for all six pooled configuration/geometry combinations.
- XY recomputes pairwise distances and fixed reference totals on the same
  275/196 reference cohorts, then recomputes pooled matrices, nulls, optima,
  and endpoint coordinates. It is a separate projection sensitivity.

The working run contains `analysis/` for raw matrices, assignments, fronts,
metrics, exact null moments, null samples, coverage, normalization reference
edges, source hashes and transform metadata; `figures/` for PNG/PDF outputs,
caption and render hashes; and `validation/replay.json` for the replay report.

| Figure stem under `figures/` | Purpose |
|---|---|
| `terminal_cross_species_comparison` | Primary six-panel comparison: 3D above XY |
| `terminal_cross_species_overlay` | Figure 9: paired 3D/XY overlays (A--B) and canonical-metric rows (C) |
| `terminal_cross_species_tracking_sensitivity` | Individual-embryo optima evaluated under pooled costs and shared pooled anchors |

Individual row exports (`terminal_cross_species_existing_3d` and
`terminal_cross_species_xy`) and individual overlay exports
(`terminal_cross_species_overlay_existing_3d` and
`terminal_cross_species_overlay_xy`) remain available for reuse. The two
primary figures above supply Figures 8 and 9. The captioned build reuses the
same numerical checkpoint and changes no metrics.

| Figure | Wrapper stem under the run's `publication/` |
|---|---|
| Figure 8: six-panel 3D/XY comparison | `fig8_terminal_cross_species_comparison` |
| Figure 9: paired overlays and canonical metrics | `fig9_terminal_cross_species_overlays` |

Each wrapper has editable TeX and a compiled one-page PDF in the same style
as the earlier terminal publication figures. Panel PDF/PNG files accompany
them. `publication_manifest.json` records figure numbering, rendered-source
identity, artifact hashes and unchanged analysis hashes. Rebuilding preserves
the preceding assembled layout with verified hashes under `layout_history/`.
This does not promote or replace `output/publication/` or edit the manuscript.

From the repository root in `dev`:

```bash
# Generate a new numerical run (an existing completed run is never overwritten).
python terminal_pareto/cross_species_analysis.py --run-id YOUR_NEW_RUN_ID

# Rebuild presentation only from the current, identity-checked caches.
python terminal_pareto/cross_species_analysis.py --render-only

# Assemble and compile captioned Figures 8 and 9 from verified panels.
python terminal_pareto/cross_species_publication.py

# Replay every assignment and check source/cache hashes without rendering.
python terminal_pareto/cross_species_analysis.py --verify-only
python -m unittest terminal_pareto.test_cross_species_terminal terminal_pareto.test_cross_species_publication
```

The C. briggsae 3D caveat belongs in the figure captions, not plot titles:
traditional embryo-tracking techniques have limitations, particularly for
z-axis measurements. Keep this wording provisional pending experimental
collaborator input; do not frame it as an unverified axial scale. Stage
alignment and RNA measurement provenance remain unresolved. Endpoint displays
compare relative front shapes and proximity; raw costs and reference cohorts
remain explicit. Protein/RNA panels differ, and tracking embryos share
molecular data. Generated outputs are Git-ignored and require the same
external backup as other runs.

## Source map and reproducibility

| Files | Role |
|---|---|
| `publication_build.py`, `publication_release.py`, `publication_wrappers.py` | Build, validated promotion and wrapper templates |
| `analysis_context.py`, `global_analysis.py`, `front_coordinates.py` | Profiles, identity-checked caches, display transforms |
| `cross_species_analysis.py`, `fig_terminal_cross_species.py`, `test_cross_species_terminal.py` | Pooled molecular/species comparison, Figure 2 style rendering, scientific checks |
| `cross_species_publication.py`, `test_cross_species_publication.py` | Figures 8/9 caption templates, one-page compilation, numbering and archive checks |
| `fig1_endpoint_amendment.py`, `fig2_fig3_*.py`, `fig4_*.py`, `fig5_*.py`, `fig6*.py`, `figS3_cross_geometry.py` | Current figure renderers |
| `data_loader.py`, `pareto_engine.py`, `lineage_metrics.py`, `subtree_analysis.py`, `subtree_explore.py`, `plot_style.py` | Shared analysis/rendering support |
| `tracking_geometry_sensitivity.py`, `figS3_ce_tracking_robustness.py` | Saved tracking sweeps and historical tracking audit support |
| `main.py` | Historical exploratory CLI, not the primary publication coordinator |
| `assets/fig3A_ce_null_models.tex` | Versioned common null schematic; no dependency on retired figure outputs |
| `validate_pooled_migration.py`, `test_analysis_context.py`, `test_publication_release.py` | Scientific, cache, and promotion checks |
| `legacy/development_sources_20260923.tar.gz` and `.json` | Exact pre-consolidation sources/docs, hashes, original paths |

Required local inputs beyond source datasets are the working run's `analysis/`
cache, `output/tracking_geometry_sensitivity/` (S3 assignments/provenance), and
`output/tracking_metric_comparison/` (audited pooled matrices and matched-subtree
validation fixtures). `output/runs/baseline_legacy_20260920/` is the immutable
numerical regression baseline. Keep these when clearing generated output.
Cache reuse checks profile, cohort/order, context identity, sweep, schema and
hashes. To change a scientific configuration, use a new explicit run/profile;
do not rename a run directory or rewrite manifests to force reuse.
Loose root-level CSV/NPZ tables are historical no-profile analysis caches,
retained for those entrypoints; the primary pooled tables are in the run's
`analysis/` directory and are identified by the published release manifest.

## Interpretation

Distances are geometric and expression proxies, not direct energetic or
regulatory work. Three embryos mix biological and tracking variation; they do
not estimate population variability or isolate measurement error. Pooled travel
summarizes those three geometries. It does not establish denoising or
out-of-sample improvement: the earlier leave-one-geometry-out comparison
improved held-out travel in only 2 of 12 settings for each method.
Small cost margins can coexist with large changes in exact parent assignments.
Do not infer biological parent agreement from slot permutations or front shape.
Full-tree analysis has different costs and references and remains a separate
pipeline. Cross-species findings require acquisition/calibration and
expression-provenance review before scientific interpretation.
