# Terminal-cell Pareto analysis

*C. elegans* terminal-cell analysis using pooled travel from three embryos and
endpoint-normalized Pareto displays, plus the 187-edge cross-species
comparison. This package is numerical only: figures, captions and Table 1 are
drawn by [`publication/`](../publication/README.md), whose manuscript set
(`publication/output/manuscript/current`) holds the manuscript's figures.
Backup figures and parked analyses are in
[`publication/BACKLOG.md`](../publication/BACKLOG.md).

## Start here

- **Working runs (numerical):** `output/runs/pooled_tracking_v1/migration_candidate_20260920/analysis/`
  (primary 275-edge analysis; the S1 projection table is in
  `analysis/s3_cross_geometry/`) and
  `output/runs/cross_species_terminal_v1/pooled_comparison_20260930/analysis/`.
  Run IDs are historical and retained to preserve cache identity.
- **Add-on runs read by the manuscript set:**
  `representation_criteria_v1/criteria_metric_20261014` (Figure S2) and the
  Track D runs `subtree_symmetry_v1/compare_20261015`,
  `apoptosis_20261019`, `apoptosis_subtrees_20261021` and
  `apoptosis_single_20261022` (Figure 4 and its ABplp note; the adapter loads
  them as `TRACK_D_RUNS`). Never delete these.
- **Archives:** `legacy/development_sources_20260923.tar.gz` (+ `.json`) holds
  the code, TeX and docs from before the 2026-09 consolidation;
  `output/legacy/development_outputs_20260923.tar.gz` holds prototypes,
  pilots, consensus, tracking-sensitivity and minimax outputs, and
  `output/legacy/pre_promotion_run_20260923.tar.gz` the run before the
  cleanup. Read members with `tar -xOf ARCHIVE MEMBER`; resume pilots only in a
  separate checkout. Consolidation is not permission to discard scientific
  inputs.

Run commands from the repository root in the `dev` Conda environment.

```bash
# Numerical pipeline (global front, subtrees, canonical metrics, cell types,
# within-type fronts, S1 projections) into the run's analysis/ folder.
python terminal_pareto/analysis_pipeline.py

# Validate scientific results, source identity and the legacy regression baseline.
python terminal_pareto/validate_pooled_migration.py   # writes validation/migration_validation.{json,md}
python terminal_pareto/test_analysis_context.py       # about 87 s

# Figures: see publication/README.md (the manuscript set and the backlog families).
```

Generated outputs, caches and archives are Git-ignored and need the
project's external backup.

## Analysis conventions

Terminal identities, measured positions, molecular states, and parent-slot
multiplicities are fixed. Hungarian bipartite matching assigns children to
parent slots across 301 weights (300 intervals including both endpoints).
These are sampled optimal assignments, not enumeration of the complete discrete
Pareto set. Duplicate slots of one biological parent are interchangeable;
edge retention (lineage similarity in the figures) compares biological parent
identity.

The primary cohort is the strict intersection of **275 natural terminal
parent–child edges** across embryo 1 (cutoff 255), embryo 2 (247), and embryo 3
(225). Both endpoints require tracking and the fixed protein representation.
Cutoffs are predefined stage matches, not equal elapsed developmental times.
A terminal cell is a leaf of `data/cell_lineage.json` observed, with its
parent, in all three embryos at the cutoffs and with protein data (Methods
definition; 46 cells that end embryogenesis undivided but are internal in the
reference tree are therefore excluded, see "Cohort" below).

Cell-state cost is cosine dissimilarity of the selected top-20 z-scored protein
features from `data/protein/aggregated_all/s3_zscore.csv`, with selection in
`expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv`
(cosine everywhere since 2026-10-14; see Track A below).
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
pipeline defaults to the pooled profile. Historical no-profile figure commands
write to `output/legacy/rebuild/publication/`.

First-cousin assignment shuffling is the primary null, with seed 42 and 1,000
visualization draws. Figure 2 also shows second-/third-cousin and full-random
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
comparisons between groups. Endpoint normalization is for display; it inflates
d_LP and d_CP for near-point fronts (for example a *C. briggsae* d_LP of 0.76
in P2 comes from a raw travel gap under 2%).

P* is the closest **attained sampled assignment** to natural lineage L in
endpoint coordinates. Canonical distance d_LP is its Euclidean distance to L;
d_CP (cache field `d_NP`) measures from first-cousin-null mean C (`N`) to that
same P*. Position u is normalized front arc length from A to P*. No
interpolated lineage is used. The maximum-edge-retention assignment need not be
P*. The cousin-relative statistic r uses separately null-standardized
geometry; it is a sensitivity statistic (Table 1), and d_LP/d_CP need not equal
r (ABplp 0.20 against 0.24).

Readouts measured from the null rather than the front:

- d̃_LP (cache `d_SD`, notation `natural_front_distance_sd`): distance to the
  nearest front point in first-cousin-null SD units; that point can differ
  from P*.
- z_T (`natural_travel_z`) and z_S (`natural_state_z`): the natural lineage's
  travel and cell-state costs as z-scores against the first-cousin null.
- d_CP explodes as a subtree's front span collapses (ABplpa 64; leave-one-out
  range 50-358), one reason for the subtree floor below.

## Subtrees, apoptotic cells and the cohort (Track D)

**Subtree floor (author, 2026-11-02).** The back-end caches keep the 12-cell
floor (`PRIMARY_MIN_CELLS = 12`, 42 subtrees). The manuscript shows only
subtrees of at least 24 cells (`SUBTREE_FLOOR = 24` in
`publication/adapters/terminal.py`): P0, AB, ABa, ABp, P1, EMS, ABal, MS, ABpl,
ABpr, ABar, ABprp, ABplp, ABalp, ABara, MSa, P2 and MSp. A subtree is resolved
when its first-cousin null shuffles at least 12 cells and its front spans at
least 1 null SD on both axes (`TINY_SPAN_SD`); the floor is the smallest size
above which every subtree is resolved. The largest unresolved subtrees are
ABala (22 cells, null on 9) and ABprpa (21, state span 0.91 SD); no subtree
has 23 or 24 cells. The floor was not moved to drop MSp or ABplp. The build
checks it against `compare_20261015`'s `subtree_metrics.csv`
(`publication/tests/test_terminal_floor.py`).

Why: readouts measured from the front (u, d_LP, d̃_LP, d_CP) can be decided by
one cell at 12-36 edges, while z_T is stable and z_S is stable unless null
support is 11 cells or fewer. The twin correlation of u is 0.42 at 24 cells,
0.49 at 48 and 0.62 at 96, so the text compares no u values between subtrees.

**Claims.** Every subtree of at least 24 cells favours travel (u < 0.5),
including P2, which holds the C and D lineages: all 18 have u < 0.5 (highest
MSp 0.38; ABalp and P2 at 0); at the 12-cell floor 41 of 42 do. Six disjoint
subtrees (ABal, ABar, ABpl, ABpr, MS, P2) cover 273 of 275 cells. 11 of 18 are
off their front; ABplp is farthest (d_LP 0.061), then MSp (0.046) and ABpl
(0.027). d_LP < d_CP in all 18. C (14 cells, resolved, u 0, on its front) and
D (12 cells, null on 8, unresolved) need no exception.

**ABplp and one apoptotic cell.** ABplp's distance from its front comes from
the endpoint of one apoptotic cell, ABplpppapp: holding that endpoint at the
cell's birth position, or dropping the cell, puts ABplp on its front, while
the same change to its other deaths (ABplpappap, ABplppaaap) does not. Matched
living controls: 0 of 1,000 draws (p = 0.001, 119 distinct control sets;
matching is weak, controls are one division deeper and born 23-68 frames
later). Its endpoint lies farther than its mirror ABprpppapp's in every embryo
(8.9/5.7, 4.8/3.4 and 7.0/4.4 µm). Removing the living ABplpppaaap or
ABplppppapa also closes the gap. The paper must not call this an artifact.

**Apoptotic cells stay in the analysis** (author, 2026-10-21). 30 of the 275
cohort cells die. Dropping them takes d̃_LP from 1.11 to 0.57, against a
matched-control median of 0.59 (p 0.47); holding their endpoints at birth gives
0.49 (p 0.31). The departure from the front therefore sits in late-born, deep
cells in general, not in deaths. No death times exist; the `as_is` and
`birth_endpoint` variants bracket any death-time endpoint, and many endpoints
are corpse positions.

**Cohort (author, 2026-11-09).** 46 cells end embryogenesis undivided but are
internal in the reference tree (10 E, 11 each in ABar, ABpl and ABpr, 3 in
MS); `data/embryonic_terminals.txt` was never used. Adding them (321 edges)
gives d_LP 0.031 instead of 0.015 and the lineage stays far better than its
first-cousin null; the 275-edge cohort was kept. Numbers and run are in
`publication/BACKLOG.md`.

**Cell keys.** Our key calls ABplpapppap "RECT HYP" where WormAtlas says K′
(whether this changes the fate-only variant was not checked); ABprpappap gives
rise to RIS and DB4.

Left-right symmetry left the paper on 2026-11-02; its reviewer defense
(mirrored pairs, simulated twins, mirror geometry) is summarized in the
backlog.

## Figures in the manuscript set

| Figure | What this package supplies |
|---|---|
| 2 | Global front, nested cousin shuffles and full-random null of the primary run |
| 3 | Lineage similarity along the global sweep; rewiring at T, M and S (assignments replayed from the global analysis) |
| 4, Table 1 | Canonical metrics of the 18 subtrees of at least 24 cells, with the Track D runs |
| 5 | Cell-type lineage similarity, within-type fronts and the type-restricted aggregate (`fig6a_*`, `fig6bc_*` caches) |
| 7 | *C. briggsae* RNA, AF16 3D, from the cross-species run |
| S1 | The 4,816 S3 projections (four travel geometries) |
| S2 | Track A criteria run |

S1 measures travel in embryos 1–3 (A–C) and the pooled geometry (D). Color
identifies the source of optimization, solid lines the native front, dashed
lines transferred assignments, and X natural lineage. Every panel evaluates
all four saved 301-solution sweeps unchanged; its inset uses the same
coordinates. Overall front shapes are similar; near natural lineage, pooled
assignments form a compromise closer to each native front than assignments
transferred from another embryo. The pooled optimum's advantage in its own
objective alone is not independent support for pooling. Framing (author):
a shared global shape, local optima that depend on tracking and a pooled
compromise; keep the holdout result; claim neither isolated tracking error nor
a "correct" geometry. Fronts look alike while parent identities vary:
mixed-weight biological-parent agreement is about 0.76–0.79, and pure-travel
transfer penalties are much larger.

Development record (descriptive only; not the reason pooled travel was
adopted): the natural lineage lies on the sampled front in 16/43 subtrees
(legacy), 15/42 (matched) and 27/42 (pooled); global d_LP 0.0204, 0.0204 and
0.0150; maximum retention 0.8194, 0.8073 and 0.8327. At equal state cost,
pooled travel lowers mean travel by 2.984% but raises embryo 3 by 0.127%;
minimax lowers it by 2.334% with all three embryos lower (data only in the
development-outputs archive).

## Cross-species terminal comparison (2026-09-30)

`cross_species_analysis.py` builds CE protein, CE RNA and CB AF16 RNA on the
same **187 terminal parent--child edges**, with identical biological-parent
capacities and endpoint displays. Run:
`output/runs/cross_species_terminal_v1/pooled_comparison_20260930/`. The
manuscript set uses its *C. briggsae* RNA 3D front as Figure 7; the
comparisons with *C. elegans* (released Figures 8 and 9, family
`terminal-cross-species`) are in the backlog since 2026-11-10.

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
  optimizer uses those SDs and 301 weights. A nested 1,201-weight check is
  saved for all six configuration/geometry combinations: maximum retention is
  unchanged, distances change by < 5e-16 and u by < 1.2e-7. Tied optima are
  not enumerated.
- XY recomputes pairwise distances and fixed reference totals on the same
  275/196 reference cohorts, then recomputes pooled matrices, nulls, optima,
  and endpoint coordinates. It is a separate projection sensitivity.

The run's `analysis/` holds raw matrices, assignments
(`assignments_and_costs.npz`, including each single embryo's optima),
fronts, metrics, exact null moments, null samples, coverage, normalization
reference edges, source hashes and transform metadata; `validation/replay.json`
replays 13,828 assignments and checks 20 source hashes.

| d_LP, maximum retention | CE protein | CE RNA | CB RNA |
|---|---|---|---|
| 3D | 0.0099, 0.8717 | 0.0107, 0.9198 | 0.0311, 0.7380 |
| XY | 0.0196, 0.7594 | 0.0494, 0.6471 | 0.0411, 0.6524 |

XY reverses the order of CE RNA and CB, so the comparison supports no species
ranking and no modality claim.

`cross_species_partitions.py` (`python -m terminal_pareto.cross_species_partitions`)
writes `analysis/partitions/` of the same run for the backlog sketches: 25
subtrees (at least 12 matched edges) and 4 cell types (`MIN_TYPE` 15), with
retention by type along the cached sweeps and no new solves; P0 reproduces the
cached comparison exactly. Alimentary's within-type first-cousin null has
zero variance, so it uses full-random scaling and has no d_CP.

```bash
# Generate a new numerical run (an existing completed run is never overwritten).
python terminal_pareto/cross_species_analysis.py --run-id YOUR_NEW_RUN_ID

# Replay every assignment and check source/cache hashes.
python terminal_pareto/cross_species_analysis.py --verify-only
python -m unittest terminal_pareto.test_cross_species_terminal
```

The *C. briggsae* 3D caveat belongs in captions, not plot titles: traditional
embryo-tracking techniques have limitations, particularly for z-axis
measurements. Keep this wording provisional pending experimental collaborator
input; do not frame it as an unverified axial scale. Stage alignment and RNA
measurement provenance remain unresolved. Endpoint displays compare relative
front shapes and proximity; raw costs and reference cohorts remain explicit.
Protein/RNA panels differ, and tracking embryos share molecular data.

## Cell-state representation (Track A, 2026-10-14)

Two compute-only add-ons support the choice of expression representation
behind the cell-state cost. The methods, parameters and results are in
[`METHODS_cell_state_representation.md`](METHODS_cell_state_representation.md),
written as the basis of the publication's methodology subsection; Figure S2
shows the selection. The rule is that the front is the outcome, never a
selection input.

- `feature_sensitivity.py` (A1, run `output/runs/feature_sensitivity_v1/survey_20261014/`)
  places the classifier top-k, bottom-k and bottom-k-within-HV-pool sets
  among 200 + 200 random sets for k in {5, 10, 20, 30}. It covers protein on
  275 edges and CE protein / CE RNA on the 217 RNA-covered edges, in 3D and
  XY.
- `representation_criteria.py` (A1b/A2, run `output/runs/representation_criteria_v1/criteria_metric_20261014/`)
  scores top-20, HV-pool and whole-panel candidates without the front:
  - sibling separation (fate discrimination);
  - lineage memory (same-type similarity by lineage distance, against a
    within-type shuffle);
  - measured vs fate-only fronts (children as fine-type means).
  - every criterion and front for cosine and Euclidean (A2, paired seeds).

  Protein uses the 275 primary edges. CE and CB RNA use the 187
  cross-species edges.

Both runs start only after the published top-20 metrics are reproduced. The
criteria run also reproduces the notebook's sibling separation. The outcome
is to keep the published top-20 sets, and protein remains the *C. elegans* gold
standard. RNA, which is resolved by type, is limited to the terminal
cross-species comparison. The fate-only front is the control for lineage
memory. Cosine replaced Euclidean everywhere because with Euclidean a few
highly expressed TFs dominate raw RNA (CE z_S −0.25 against −3.14 with
cosine).

```bash
python terminal_pareto/feature_sensitivity.py --run-id NEW_ID        # survey; --verify-only to replay
python terminal_pareto/representation_criteria.py --run-id NEW_ID    # criteria; --check-only runs the gates
python -m unittest terminal_pareto.test_feature_sensitivity terminal_pareto.test_representation_criteria
```

In a git worktree the context cache key (an absolute tracking path) differs,
so both modules check matrix and source hashes instead of the context
manifest. `data/embryo*` must be real copies there, not links.

## Source map and reproducibility

Numerical code only; drawing, captions, notation (the cache field `d_NP` is
displayed as `d_CP`, null mean `C`) and Table 1 formatting are in
[`publication/`](../publication/README.md).

| Files | Role |
|---|---|
| `analysis_pipeline.py` | Compute-only orchestration of the pooled terminal caches |
| `analysis_context.py`, `global_analysis.py`, `front_coordinates.py` | Profiles, identity-checked caches, display transforms |
| `cross_species_analysis.py`, `test_cross_species_terminal.py` | Pooled molecular/species comparison and scientific checks (`fig_terminal_cross_species.py` only redirects `--render-only`) |
| `cross_species_partitions.py` | Cross-species comparison by subtree and cell type (backlog sketches) |
| `feature_sensitivity.py`, `representation_criteria.py` and their tests | Cell-state representation survey and criteria (compute-only; see `METHODS_cell_state_representation.md`) |
| `subtree_symmetry.py`, `subtree_symmetry_robustness.py` | Track D: mirrored pairs, subtree comparison and cell tables (runs `subtree_symmetry_v1/pairing_20261014`, `compare_20261015`, `cells_20261015`, `robustness_20261015`) |
| `apoptosis_sensitivity.py`, `apoptosis_subtree_control.py`, `apoptosis_single_deaths.py` | Track D4: apoptotic-cell variants and matched controls (runs `subtree_symmetry_v1/apoptosis_*`) |
| `cohort_sensitivity.py` | Track D7: the 321-edge cohort check (`cohort_sensitivity_v1/embryonic_terminals_20261108`) |
| `fig2_fig3_*.py`, `fig5_table1_*.py`, `fig6a_*.py`, `fig6bc_*.py`, `figS3_cross_geometry.py` | Analyses behind the figures (historical file names). `fig5_table1_*` and `figS3_cross_geometry.py` stay byte-identical because validators pin their source; they still write presentation side-files, which the pipeline discards |
| `data_loader.py`, `pareto_engine.py`, `lineage_metrics.py`, `subtree_analysis.py`, `subtree_explore.py` | Shared analysis support; `plot_style.py` re-exports `publication/style.py` for historical scripts |
| `tracking_geometry_sensitivity.py`, `figS3_ce_tracking_robustness.py` | Saved tracking sweeps and historical tracking audit support |
| `main.py` | Historical exploratory CLI, not the primary publication coordinator |
| `validate_pooled_migration.py`, `test_analysis_context.py` | Scientific and cache checks: 23 checks, including exact legacy reproduction and replay of all 4,816 S1 projections (max error 4.55e-13); figure checks are in `publication/tests` |

**Hash-pinned modules, never edit:** `analysis_context.py`, `data_loader.py`,
`lineage_metrics.py`, `pareto_engine.py`, `front_coordinates.py`,
`subtree_analysis.py` (it also pins the pooled full-tree caches),
`global_analysis.py`, `cross_species_analysis.py`,
`fig5_table1_ce_canonical_metrics.py`, `figS3_cross_geometry.py`,
`feature_sensitivity.py` and `representation_criteria.py`; each Track D add-on
pins its own run. Import them and write add-on modules with their own run
folder, source hashes, replay validation and tests (pattern:
`full_tree_pareto/cousin_references.py`). `solve_block_sweep` and
`exact_null_stats` are not pinned and are reusable. The cell-type code assumes
the 275-edge context.

Required local inputs beyond source datasets are the working run's `analysis/`
cache, `output/tracking_geometry_sensitivity/` (S1 assignments/provenance), and
`output/tracking_metric_comparison/` (audited pooled matrices and matched-subtree
validation fixtures). `output/runs/baseline_legacy_20260920/` is the immutable
numerical regression baseline. Keep these when clearing generated output.
Cache reuse checks profile, cohort/order, context identity, sweep, schema and
hashes. To change a scientific configuration, use a new explicit run/profile;
do not rename a run directory or rewrite manifests to force reuse.
Loose root-level CSV/NPZ tables are historical no-profile analysis caches,
retained for those entrypoints; the primary pooled tables are in the run's
`analysis/` directory.

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
