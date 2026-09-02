# Full-tree Pareto analysis

Updated 2026-09-02. This module extends the terminal-cell analysis to the
measured *C. elegans* lineage containing both internal and terminal cells. The
publication pipeline produces two candidate main figures:

| Item | Role | TeX/PDF stem | LaTeX label |
|---|---|---|---|
| Figure 7 | Bottom-up layerwise assignment and aggregate front | `fig7_ce_full_tree_layerwise` | `ce_full_tree_layerwise` |
| Figure 8 | Collective heuristic front and notebook null models | `fig8_ce_full_tree_collective` | `ce_full_tree_collective` |

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

## Figure 8: collective heuristic envelope

Figure 8 compares six reconstruction strategies from the notebook's analysis
program:

1. bottom-up optimization by lineage layer;
2. independent layerwise parent assignment from Figure 7;
3. paired bottom-up rebuilding;
4. degree-constrained spanning-forest rebuilding;
5. top-down rebuilding; and
6. terminal-only rebuilding without reuse of measured internal protein states.

Each colored curve contains only the non-dominated cost pairs sampled for that
strategy. The strategies are plotted against one another rather than as a
single union, since the envelope is a sampled collection and not a proof of the
exact global Pareto front. The non-dominated solutions are dominated by the
degree-constrained spanning forest, top-down rebuilding, and terminal-only
rebuilding; the other three strategies remain scientifically informative
comparisons but rarely reach the shared envelope.

The degree-constrained spanning forest and the top-down rebuild are nearly
identical because both are greedy reconstructions over the same combined cost
matrix `alpha*xyz + (1-alpha)*exp` with the same structural constraints (four
fixed roots, at most two children per parent, and exactly 1,000 scored edges
with terminals assigned by linear-sum assignment). At matched weights the two
fronts differ by only a few travel and cell-state units out of thousands, so
they overlap almost exactly on the main axes.

All collective coordinates use the cached first-cousin terminal shuffle as the
common scale and are translated so the natural lineage is at `(0, 0)`. The five
notebook null models are retained:

- first-cousin terminal shuffle;
- internal-identity shuffle within represented lineage depths;
- full internal/terminal identity shuffle;
- random full-tree rebuild; and
- a fixed-topology parametric Brownian reference.

The first-cousin null and parametric Brownian reference remain on the main
axes. The broad random models use an inset, following the terminal-only
presentation convention for nulls outside the biologically relevant front
range. Endpoint symbols and repeated strategy-contribution dots are omitted:
the curve endpoints are already explicit, and the colored strategy curves
identify their own locations without a second point layer.

### Audited parametric Brownian reference

The historical `phylo_bm.npz` cache is excluded from publication analysis for
two independent reasons. First, its cost scale and dimensionality identify it
as output from the older 3-spatial-plus-10-PC experiment; it cannot be combined
with the current 3-spatial-plus-20-z-scored-protein analysis.
Second, the exploratory implementation indexed measured arrays with tree IDs
instead of mapped biological lineage IDs. None of the 504 terminal tree IDs in
this configuration equals its biological data-row ID. The routine then sampled
every reconstructed node independently, destroying the parent--child
covariance that defines Brownian motion.

The publication replacement is an all-edge plug-in parametric bootstrap. It
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
strong working assumptions. Its role in Figure 8 is limited to probing where a
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
diagnostics cache records this concentration and related rate checks.

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

The very large cell-state offset in the parametric Brownian reference remains
an open model-adequacy question and must be reviewed with the authors before
the reference is interpreted or revised. A read-only sensitivity audit on
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

If the authors decide that the reference should be revised, evaluate these as
explicit sensitivity models rather than silently replacing the current null:

1. a block-diagonal travel/protein covariance for clearer interpretation,
   while recognizing that it cannot move the marginal null positions;
2. layer- or branch-class-specific rate matrices;
3. a heavy-tailed scale-mixture Brownian model that permits many quiet edges
   and a small number of large regulatory changes; and
4. a semiparametric residual bootstrap, stratified by branch duration or
   contraction round, to separate distributional shape from time-scaling
   failure.

Any candidate must retain the four fixed roots, the measured topology, all
1,000 scored edges, and deterministic seeds. Compare its cost distribution and
edgewise residual diagnostics with the current reference before changing
Figure 8. Until that review, keep the current cloud labeled as a geometric
parametric reference and do not read its cell-state displacement as an
efficiency effect size.

## Code and cache map

```text
full_tree_pareto/
├── README.md
├── publication_analysis.py              # Validated context, manifests, and caches
├── fig7_fig8_ce_full_tree_pareto.py      # Publication renderer
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
- `ce_full_tree_brownian_covariance.csv`; and
- `ce_full_tree_brownian_diagnostics.csv`.

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
fixed-topology parametric Brownian reference, reuses the other audited notebook caches, and writes
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
```

## Validation before handoff or commit

1. Compile `publication_analysis.py` and the renderer with `py_compile`.
2. Run the cache-only analysis and confirm every invariant passes.
3. Confirm that the eight contraction rounds contain exactly 1,000 unique
   evaluated edges and recompose both natural objective totals.
4. Confirm every spanning-forest solution has four components, at most two
   children per parent, and 1,000 scored edges.
5. Regenerate and visually inspect all three component panels.
6. Compile both wrappers and confirm one-page letter output with embedded fonts.
7. Run `git diff --check` and exclude unrelated or expensive generated caches
   from the staged set.

Shared visual semantics follow `terminal_pareto/plot_style.py`: sequential blue
encodes natural-edge retention, natural lineages use a dominant black cross,
nulls use muted local colors with explicit labels, axes use sentence case, and
PDF is canonical.
