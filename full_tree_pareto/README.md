# Full-tree Pareto analysis

Updated 2026-09-16. This module extends the terminal-cell analysis to the
measured *C. elegans* lineage containing both internal and terminal cells. The
publication pipeline produces two candidate main figures:

| Item | Role | TeX/PDF stem | LaTeX label |
|---|---|---|---|
| Figure 7 | Bottom-up layerwise assignment and aggregate front | `fig7_ce_full_tree_layerwise` | `ce_full_tree_layerwise` |
| Figure 8 | Layerwise assignment versus constrained reconstruction | `fig8_ce_full_tree_collective` | `ce_full_tree_collective` |

The exploratory notebook remains useful development history, but
`publication_analysis.py` and `fig7_fig8_ce_full_tree_pareto.py` are the
authoritative analysis and presentation paths for these figures.

## Fixed publication configuration

- Organism and molecular data: *C. elegans*, top-20 z-scored protein features.
- Tracking: embryo 1 in three dimensions through `T <= 255`, scaled by 0.1625.
- Molecular distance: Euclidean distance across the 20 protein features.
- Spatial distance: Euclidean distance across final tracked `(x, y, z)`.
- Objective sweep: 300 intervals and 301 endpoint-inclusive travel weights.
- Early anchors: P0, AB, and P1 are required to construct the represented tree
  but have no measured coordinates or protein state and do not contribute
  scored edges.
- Evaluated forest: four earliest measured roots, 500 internal cells, 504
  terminal cells, and exactly 1,000 parent--child edges.

The represented tree contains 1,007 nodes including the three early anchors.
It is not the complete invariant lineage atlas: it is the complete internal and
terminal tree within the joint tracking/protein measurement space. Missing
intermediate atlas cells are skipped while their measured descendants remain
attached to the nearest represented ancestor.

Travel and cell-state objectives remain separate throughout. A weighted sum is
used only to select assignments at each sampled weight and must not be
interpreted as a biological efficiency score.

## Figure 7: bottom-up layerwise assignment

The layerwise analysis starts with every measured terminal cell and contracts
the represented tree toward the four roots. Because branches leave the
measurement space at different canonical depths, these are asynchronous
bottom-up contraction rounds, not uniform developmental layers.

| Round | Position | Reassigned edges |
|---:|---|---:|
| 1 | bottom | 504 |
| 2 | intermediate | 230 |
| 3 | intermediate | 126 |
| 4 | intermediate | 68 |
| 5 | intermediate | 38 |
| 6 | intermediate | 19 |
| 7 | intermediate | 10 |
| 8 | top | 5 |

The eight rounds partition all 1,000 evaluated edges exactly. Within a round,
the child identities are assigned by the Hungarian algorithm to the multiset
of their natural parent slots. The number of slots belonging to each parent is
fixed, but a child may move to any slot in the round.

Each round receives its own matched parent-slot null: 1,000 deterministic random
permutations of those same slots. Displayed travel and cell-state coordinates
are divided by their respective null standard deviations and translated so the
natural assignment is at `(0, 0)`. The aggregate front sums all eight optimized
rounds at the same weight. Its null sums independent per-round permutations.

Edge retention uses biological parent identity, not assignment-slot identity.
Two indistinguishable slots belonging to the same parent therefore retain both
natural edges when exchanged. This corrects the exploratory cache, which used
the stricter diagonal-slot comparison.

The natural assignment is sampled-dominated in rounds 1--5 and in the aggregate
assignment space, but not in rounds 6--8. The aggregate contains 202 unique
sampled solutions, and its maximum-retention solution preserves 49.9% of the
1,000 natural edges.

Accepted components:

- `fig7A_ce_full_tree_layerwise_rounds.{pdf,png}`; and
- `fig7B_ce_full_tree_layerwise_aggregate.{pdf,png}`.

## Figure 8: layerwise assignment versus constrained reconstruction

The presentation adopted on 2026-09-09 shows two complementary methods:

1. **Layerwise assignment:** exact Hungarian assignment within each contraction
   round, connecting terminal-only linear assignment and Figure 7. Natural
   parent-slot multiplicities are preserved.
2. **Degree-constrained spanning forest:** greedy reconstruction over measured
   internal states followed by terminal assignment; four measured roots,
   binary branching and 1,000 scored edges are validated at every weight.

Each curve contains its own non-dominated cost pairs from 301 sampled weights.
This compares relaxation of assignment constraints, not an exact global Pareto
front. Selection is based on distinct scientific roles rather than smoothness
or the number of sampled solutions. Top-down rebuilding nearly overlaps the
forest and is retained in the supplement as an algorithmic comparison.

Figure 8 retains first-cousin shuffle, the separate-clock Gaussian reference
and one broad null, random rebuild. All use ordinary continuous axes, without
an inset or axis breaks. Random rebuild preserves four roots and complements
the optimized forest; full-assignment shuffle is nearly redundant in cost
space (mean differences 0.39 travel and 1.12 protein units, versus cloud SDs
around 215 and 31). The distinct internal-layer shuffle and full-assignment
shuffle remain in the supplement. Coordinates use first-cousin-shuffle standard
deviations and place the natural lineage at zero. Plus signs use full
distribution means from the summary cache, not display-subset means.

### Supplementary comparison and six-method inventory

`figs_ce_full_tree_heuristics.tex` presents five cached strategy curves and
all five reference/null families on continuous axes. Solid curves are the main
methods; dashed curves add top-down, bottom-up-by-layer and paired bottom-up
comparisons. These additional caches have 1,001 sampled weights. Bottom-up by
layer and paired bottom-up have different root-identity freedoms and must not
be described as equally constrained, fully validated alternatives.

`table_ce_full_tree_heuristics.tex` lists all six methods, including the held
terminal-only method. `heuristic_inventory.py` generates the table rows and
CSV, recording measured versus constructed states, root constraints, edge
scope, algorithm and audit status. Only the two main methods are regenerated
with structural invariants by the publication pipeline; the table distinguishes
this from inspection of historical scalar-cost caches. Supplementary numbering
is provisional until manuscript integration.

### Held terminal-only comparison: scope and scoring issue

The earlier blanket description of all six curves as comparable 1,000-edge,
four-fixed-root reconstructions was too strong. Inspection on 2026-09-09 found
that `terminal_only_rebuild(..., use_kd_tree=True, exp_replacement=False)`:

- constructs internal spatial positions at daughter midpoints;
- chooses measured internal protein profiles without replacement, so
  "terminal-only" does not mean that it ignores internal protein observations;
- does not fix four measured root states; and
- removes the final available protein profile, then breaks before scoring its
  pair. Valid deterministic pairings yield 499 scored merges, or 998 edges.

This diagnostic did not regenerate the saved optimized curve.
`terminal_only_no_replacement.npz` is retained as historical data but withheld
from both main and supplementary plots pending a dedicated audit. The inventory
lists the method without presenting its cached costs as validated results.
Moving a method to the supplement does not repair a scope or scoring issue.
The historical `collective_nondominated` cache column remains for provenance;
figures and reports no longer use its six-method union to rank methods or
claim a common feasible set. The apparent terminal-only travel advantage
must not be interpreted as superiority under the main methods' constraints.

### Historical Brownian reference (superseded 2026-09-09)

The historical `phylo_bm.npz` cache is excluded from publication analysis for
two independent reasons. First, its cost scale and dimensionality identify it
as output from the older 3-spatial-plus-10-PC experiment; it cannot be combined
with the current 3-spatial-plus-20-z-scored-protein analysis.
Second, the exploratory implementation indexed measured arrays with tree IDs
instead of mapped biological lineage IDs. None of the 504 terminal tree IDs in
this configuration equals its biological data-row ID. The routine then sampled
every reconstructed node independently, destroying the parent--child
covariance that defines Brownian motion.

The initial publication replacement was an all-edge plug-in Brownian bootstrap. It
maps biological identities explicitly, fits one full 23-by-23 covariance to
the 1,000 measured increments after division by the square root of represented
branch time, fixes the four optimization-root states, and simulates every
scored descendant jointly down the same topology. Each draw is evaluated on
exactly the 1,000 edges used for the natural lineage. The pipeline generates
10,000 deterministic draws for numerical summaries and marks a deterministic
1,000-draw subset for display.

This cloud is intentionally described as a reference rather than a realistic
embryogenesis null. Zero drift, independent edge increments, variance linear
in branch time, one common covariance, and negligible observation error are
strong working assumptions. Its former role in Figure 8 was limited to probing where a
simple neutral diffusion process on the fixed lineage topology lies relative
to the natural lineage and the collective sampled Pareto envelope. Proximity
or separation must not be interpreted as evidence that embryogenesis follows
Brownian dynamics.

With seed 242, the 10,000-draw reference averages 4,831 travel units and 5,241
cell-state units, compared with 4,466 and 3,102 for the natural lineage. In the
Figure 8 coordinates these differences are approximately +9.42 and +385.27
first-cousin-null standard deviations. The extreme cell-state separation is
consistent with poor common-rate model adequacy rather than a precise measure
of biological inefficiency: the largest 5% of observed branch-standardized
cell-state increments contribute about 56% of the fitted squared rate, which a
homogeneous Gaussian generator redistributes across all branches. The
historical diagnostics cache records this concentration and related rate checks.
This audited sampler remains available through `parametric_brownian_reference()`
and `clock_sensitivity.py` for comparison, but is no longer a Figure 8 input.

### Audited spanning-forest replacement

The historical `mst_rebuild.npz` cache has been removed and is not an input to
the publication analysis. The saved cache did include both internal and
terminal costs, but its unconstrained MST is one unrooted component, does not
preserve the four biological roots, and permits internal degrees incompatible
with a binary lineage. The current legacy `mst_rebuild` implementation has a
separate early-return bug that reports only terminal-to-opening costs, so it
cannot reproducibly regenerate the old saved curve either.

The publication replacement preserves the intended spanning-tree idea while
enforcing the common optimization scope. At every weight it:

- holds four roots in separate components;
- uses greedy degree-constrained Kruskal selection over the 500 internal cells;
- prevents a component from acquiring two roots;
- limits every rooted internal cell to at most two children;
- exposes exactly 504 terminal openings and assigns every terminal once; and
- scores 496 internal plus 504 terminal edges, totaling 1,000.

Do not restore the legacy MST curve to Figure 8 without a different audited
replacement satisfying the same edge and topology constraints.

The replacement result cache is `output/ce_full_tree_spanning_forest.csv`.
It is intentionally not named `mst_rebuild.npz`: the publication pipeline
recomputes and validates the constrained forest before writing this CSV, which
keeps it distinct from the incompatible exploratory method and cache format.

Accepted component: `fig8_ce_full_tree_collective_panel.{pdf,png}`.

### Backlog: cell-state Brownian rate heterogeneity

The very large cell-state offset in the former Brownian reference motivated
the separate-clock replacement adopted below. Remaining protein heterogeneity
is an open model-adequacy question. A read-only sensitivity audit on
2026-09-02 found that covariance coupling is not its cause:

- the mean absolute travel--protein correlation in the fitted covariance is
  0.026, and the bootstrap correlation between total travel and cell-state
  costs is 0.024;
- setting the travel--protein covariance blocks to zero leaves both marginal
  cost distributions unchanged by construction;
- diagonalizing the 3-by-3 travel block changes the expected travel cost only
  from approximately 4,834 to 4,836; and
- treating all 20 proteins as independent increases the expected cell-state
  cost from approximately 5,242 to 5,273, rather than moving it toward the
  natural value of 3,102.

The stronger diagnosis is a mismatch between the homogeneous Brownian rate
fit and heterogeneous measured protein increments. The top 5% of
branch-standardized protein increments account for about 56% of the fitted
squared rate but only about 11% of the natural sum-of-norms cost. The weighted
mean standardized protein norm is 0.518 whereas its root-mean-square value is
0.896. A common Gaussian rate estimated from squared increments redistributes
the influence of rare large changes over all 1,000 branches, inflating the
sum of unsquared edge norms used in Figure 8.

Branch-duration scaling is also poorly supported for protein state. A log--log
regression of squared increment magnitude on represented branch duration gives
a slope of approximately 0.06 for protein state, compared with the Brownian
prediction of 1 and an observed travel slope of approximately 0.72. Adding a
single fitted drift vector does not resolve the discrepancy; its expected
cell-state cost remains approximately 5,262.

The initial audit proposed the following sensitivity models; clock mismatch
was subsequently prioritized and resolved by the separate-clock convention:

1. a block-diagonal travel/protein covariance for clearer interpretation,
   while recognizing that it cannot move the marginal null positions;
2. layer- or branch-class-specific rate matrices;
3. a heavy-tailed scale-mixture Brownian model that permits many quiet edges
   and a small number of large regulatory changes; and
4. a semiparametric residual bootstrap, stratified by branch duration or
   contraction round, to separate distributional shape from time-scaling
   failure.

Any further candidate must retain the four fixed roots, the measured topology, all
1,000 scored edges, and deterministic seeds. Compare its cost distribution and
edgewise residual diagnostics with the current reference before changing
Figure 8. Keep the adopted cloud labeled as a geometric
parametric reference and do not read its cell-state displacement as an
efficiency effect size.

### Official separate-clock Gaussian reference (adopted 2026-09-09)

**Manuscript decision reaffirmed 2026-09-16:** retain the root-fixed,
leaf-unconstrained separate-clock Gaussian. Both covariance blocks are fitted
to all 1,000 observed internal and terminal edges, then the fitted process
generates all descendants without requiring the observed leaf outcomes.
The branch-variability and leaf-conditioned alternatives remain reproducible
side analyses, not inputs to the main or supplementary manuscript figures.
Conditioning on leaves with an all-edge covariance estimate is mathematically
valid, but asks a different endpoint-constrained question; it is not rejected
as an implementation error. No model is selected merely for proximity to the
natural total, and none establishes biological efficiency.

Following the clock audit and author approval, Figure 8 now uses
`publication_analysis.separate_clock_reference()`, with comparisons evaluated
by `clock_sensitivity.py`. The decision is based on compatibility between the
measurement windows and the assumed clock, not on tuning the cloud toward
the observed costs. The former full-covariance shared-time Brownian model
remains an explicitly historical sensitivity comparison. The authoritative
methods specification is `methods/separate_clock_reference.tex`; the figure
caption records the model, cutoff rationale and remaining limitations.

The cutoff audit identified a clock/measurement compatibility concern beyond
generic rate heterogeneity. Of 504 terminal cells, 492 end at tracking cutoff
255. Of the 50 largest time-standardized protein increments, 49 end at that
cutoff, with median represented duration 5. The 118 edges with duration at
most 10 contribute 63.9% of the fitted protein squared rate. Protein inputs
are atlas aggregates, not endpoint measurements aligned to tracking embryo 1.
For example, `ABprpapppa -> ABprpapppap` is assigned tracking duration 5,
while its dominant reporter ZAG-1 observes the daughter at 67 distinct times
on a different timeline. A cutoff would not itself invalidate Brownian motion
if protein states were measured at the corresponding endpoints.

The model uses independent, zero-mean Gaussian edge increments with:

- spatial covariance `t_e * Sigma_xyz`, where
  `Sigma_xyz = mean(outer(delta_xyz, delta_xyz) / t_e)`;
- protein covariance `Sigma_protein`, where
  `Sigma_protein = mean(outer(delta_protein, delta_protein))`; and
- zero spatial/protein cross-covariance, retaining full covariance within
  each block. The displayed outer products refer to individual edge vectors.

Protein variance is per canonical transition, not per minute. All 1,000 edges
are validated to span one canonical generation. No drift, fitted time
exponent, extra rate class, or heavy-tail parameter is introduced. The shared
sampler fixes four observed roots, generates descendants along the same tree,
validates reconstructed parent--child increments, and scores the same edges.
The two blocks of the exported covariance have different clock units.

With 10,000 draws and seed 242:

| Objective | Observed | Original mean | Separate-clock mean | Separate-clock central 95% interval |
|---|---:|---:|---:|---:|
| Travel | 4,465.55 | 4,831.32 | 4,831.80 | 4,690.48--4,977.13 |
| Protein state | 3,102.41 | 5,240.78 | 3,603.77 | 3,554.64--3,653.72 |

In Figure 8 coordinates the adopted mean is approximately (+9.43, +90.33)
first-cousin-null standard deviations from the natural lineage; these are
plotting units, not standard deviations of the Gaussian reference itself.
The protein excess falls from 68.9% to 16.2%. The spatial marginal
distribution is unchanged by construction; finite simulation means differ.
A second 10,000-draw run with seed 243 gives separate-clock means 4,831.11
and 3,603.97, supporting numerical stability.

The revised model still fails useful adequacy checks. Observed protein edge
lengths have coefficient of variation 0.644, versus simulated mean 0.221.
The largest 5% of raw squared protein increments account for 30.6% in the
data versus 11.3% in the new model. Observed mean internal-edge length is
2.558 and terminal-edge length 3.638, whereas the new model predicts about
3.604 for either group. Thus removing protein time scaling exposes remaining
distributional and branch-class heterogeneity; it does not establish a
realistic developmental process or a biological efficiency effect.

The source audit records dominant-protein observation windows for the 12
largest original rate outliers. Three of the top 50 outlier edges touch at
least one missing selected-protein value in `s3.csv`, which was zero-filled
before z-scoring. These flags concern upstream provenance; the sensitivity
does not change preprocessing. Reporter `Time` and tracking `t` must not be
substituted for one another without alignment. Heavy-tailed or class-specific
extensions were initially deferred pending these diagnostics and held-out evaluation.
The one-parameter branch-variability sensitivity below now provides that first
evaluation; it does not replace the official reference.

### Branch-variability sensitivity (2026-09-10; not adopted)

`branch_variability_sensitivity.py` adds a mean-one lognormal protein variance
multiplier per edge, with one fitted variability parameter. Covariance remains
the training-edge second moment; the new parameter is fitted by conditional
20-dimensional increment likelihood, not by matching total cost. Spatial draws
and the 1,000-edge, four-root scope are unchanged.

The fitted parameter is 1.22525. With 10,000 draws, the protein mean falls from
3,604 to 2,988 versus 3,102 observed, but edge CV increases to 0.724 versus
0.644 observed and top-5% squared-change share to 36.3% versus 30.6%. Held-out
vector log density improves in all eight depth-three subtrees, while seven of
eight mixture total prediction intervals still miss the observation. A closer
pooled total therefore masks remaining lineage-specific discrepancies.

Figure 8 remains unchanged. See [BRANCH_VARIABILITY_SENSITIVITY.md](BRANCH_VARIABILITY_SENSITIVITY.md)
for model specification, results, source audit, held-out protocol, limitations
and reproduction commands. Outputs are isolated in `output/branch_variability/`.

### Leaf-conditioned Gaussian sensitivity (2026-09-10; not adopted)

`leaf_conditioned_reference.py` conditions the separate-clock Gaussian on all
504 measured leaf positions and protein profiles, in addition to the four
roots. The remaining 496 internal states are sampled jointly through the
weighted tree precision matrix. Covariance calibration is unchanged; this is
not a lognormal-mixture extension or a held-out parameter fit.

With 10,000 draws, expected travel falls to 4,065.82 versus 4,465.55 observed,
and protein cost to 2,935.67 versus 3,102.41 observed. A second seed reproduces
the result. Fixing endpoints reduces uncertainty and smooths internal states;
it does not necessarily increase scored travel. The model still fails edge-
heterogeneity and internal/terminal checks, and allows internal configurations
without developmental feasibility constraints. Neither these costs nor costs
at conditional mean states establish biological inefficiency or an exact optimum.

Figure 8 remains unchanged. See [LEAF_CONDITIONED_REFERENCE.md](LEAF_CONDITIONED_REFERENCE.md)
for the conditioning equations, results, limitations, checks and reproduction
commands. Outputs are isolated in `output/leaf_conditioned/`.

### Reproducing the original clock comparison

Run the standalone comparison without recomputing Pareto optimizations:

```bash
env OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache \
  XDG_CACHE_HOME=/tmp/embryogenesis_font_cache conda run -n dev \
  python full_tree_pareto/clock_sensitivity.py

# Optional numerical-only run or independently seeded comparison:
conda run -n dev python full_tree_pareto/clock_sensitivity.py \
  --analysis-only --seed 243 --output-dir full_tree_pareto/output/clock_sensitivity/seed243

conda run -n dev python -m unittest full_tree_pareto.test_clock_reference
```

Default outputs are in `output/clock_sensitivity/`:

- `draws.csv`, `summary.csv`: both models, seeds, totals, central intervals,
  and descriptive lower-tail counts (not calibrated p-values);
- `diagnostics.csv`: observed and simulated dispersion, tail concentration,
  and grouped mean edge lengths;
- `edge_quantiles.csv`: observed quantiles and within-replicate simulation
  quantiles with central 95% intervals;
- `covariances.csv`: fitted covariances with explicit clock labels;
- `edges.csv`, `observation_windows.csv`: identity, clock, increment, and
  source-measurement audit records; and
- `clock_comparison.{pdf,png}`: reference clouds, protein totals, edge
  distributions and overlapping branch-group comparisons.

The regression checks cover unchanged original simulation results,
determinism, scope/cost recomposition, protein-fit invariance to tracking-time
rescaling, and Gaussian second moments for both clocks. The full publication
pipeline also passes with the shared sampler. A publication integration check
verifies selection of the adopted model and rejection of historical generator
provenance. Both figure wrappers and the methods documents compile; Figure 8
is visually inspected as a one-page letter PDF with embedded fonts.

## Code and cache map

```text
full_tree_pareto/
├── README.md
├── publication_analysis.py              # Validated context, manifests, and caches
├── clock_sensitivity.py                 # Separate-clock reference and diagnostics
├── test_clock_reference.py              # Clock and sampler regression checks
├── methods/separate_clock_reference.tex # Authoritative Figure 8 model and rationale
├── fig7_fig8_ce_full_tree_pareto.py      # Publication renderer
├── heuristic_inventory.py               # Constraint/audit inventory and table generator
├── figs_ce_full_tree_heuristics.tex      # Supplementary comparison wrapper
├── table_ce_full_tree_heuristics.tex     # Six-method inventory wrapper
├── run_pipeline_multi.py                 # Parameterized exploratory runner
├── run_pipeline.py                       # Legacy single-configuration runner
├── full_tree_pareto.ipynb                # Exploratory notebook
└── output/
    ├── internal_opt/                     # Expensive notebook caches
    └── publication/                      # Panels, wrappers, and assembled PDFs
```

Publication-oriented numerical caches in `output/` are:

- `ce_full_tree_node_manifest.csv`;
- `ce_full_tree_layer_manifest.csv`;
- `ce_full_tree_layerwise_fronts.csv`;
- `ce_full_tree_layerwise_nulls.csv`;
- `ce_full_tree_spanning_forest.csv`;
- `ce_full_tree_collective_heuristics.csv`;
- `ce_full_tree_collective_nulls.csv`;
- `ce_full_tree_null_summary.csv`;
- `ce_full_tree_separate_clock_covariance.csv`; and
- `ce_full_tree_separate_clock_diagnostics.csv`.

The covariance cache labels each feature's clock, and the diagnostics cache
records the adopted generator, seed, both clocks, and block independence.
Old `ce_full_tree_brownian_*.csv` files, when present, are historical artifacts
and are not read or overwritten by the current publication pipeline.

The node manifest records measurement provenance, represented and canonical
parents/depths, terminal status, optimization roots, and contraction-round
membership. The layer manifest records parent-slot multiplicity, natural and
null costs, null scales, sampled-dominance status, unique-solution counts, and
maximum biological edge retention.

All generated output is ignored by Git. Publication TeX wrappers contain
captions and layout and must be explicitly force-added when versioned.

## Regeneration

From the repository root with the `dev` Conda environment:

```bash
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python full_tree_pareto/fig7_fig8_ce_full_tree_pareto.py
```

The full run reconstructs and validates the manifests, recomputes the
301-weight layerwise and spanning-forest analyses, recomputes the all-edge
fixed-topology separate-clock Gaussian reference, reuses the other audited notebook caches, and writes
both panel sets. To generate numerical caches only:

```bash
conda run -n dev python \
  full_tree_pareto/fig7_fig8_ce_full_tree_pareto.py --analysis-only
```

Compile the standalone figures:

```bash
cd full_tree_pareto/output/publication
conda run -n dev tectonic fig7_ce_full_tree_layerwise.tex
conda run -n dev tectonic fig8_ce_full_tree_collective.tex

# From the repository root, compile the authoritative methods specification:
# conda run -n dev tectonic full_tree_pareto/methods/separate_clock_reference.tex \
#   --outdir full_tree_pareto/output/publication
```

## Validation before handoff or commit

The renderer also writes the supplementary panel and inventory rows. From the
repository root, compile their wrappers with:

```bash
conda run -n dev tectonic full_tree_pareto/figs_ce_full_tree_heuristics.tex \
  --outdir full_tree_pareto/output/publication
conda run -n dev tectonic full_tree_pareto/table_ce_full_tree_heuristics.tex \
  --outdir full_tree_pareto/output/publication
```

Supplementary outputs are `output/publication/figs_ce_full_tree_heuristics_panel.{pdf,png}`,
`ce_full_tree_heuristic_inventory.csv`, and `ce_full_tree_heuristic_inventory_rows.tex`.
The wrapper sources live outside ignored output directories so captions remain
versionable. Compile both after generating the panel and table rows.

1. Compile `publication_analysis.py` and the renderer with `py_compile`.
2. Run the cache-only analysis and confirm every invariant passes.
3. Confirm that the eight contraction rounds contain exactly 1,000 unique
   evaluated edges and recompose both natural objective totals.
4. Confirm every spanning-forest solution has four components, at most two
   children per parent, and 1,000 scored edges.
5. Regenerate and visually inspect the main panels and supplementary comparison.
6. Compile figure and inventory wrappers; confirm one-page letter output
   (landscape for the inventory) with embedded fonts.
7. Run `git diff --check` and exclude unrelated or expensive generated caches
   from the staged set.

Shared visual semantics follow `terminal_pareto/plot_style.py`: sequential blue
encodes natural-edge retention, natural lineages use a dominant black cross,
nulls use muted local colors with explicit labels, axes use sentence case, and
PDF is canonical.
