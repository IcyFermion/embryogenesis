# Publication-figure handoff

Updated 2026-09-23.

## Terminal figures: pooled release promoted

The author approved the pooled-travel, endpoint-normalized terminal figures.
The current release is `terminal_pareto/output/publication/`, including the
standalone Figure 1 amendment, single-panel Figure 2, and four-geometry S3
with insets. Previous single-embryo figures are in
`terminal_pareto/output/legacy/embryo1/publication/`.

Use [the terminal README](terminal_pareto/README.md) for analysis conventions,
figure paths, and build commands, and [the terminal handoff](terminal_pareto/HANDOFF.md)
for current validation, short development history, archived pilots/prototypes,
and remaining manuscript/co-author integration. These supersede the older
terminal sections of this repository-level document. The manuscript remains
unchanged by publication promotion. Full-tree work below remains separate.

## Full-tree Pareto analysis (candidate Figures 7--8)

The *C. elegans* protein configuration is now implemented for the represented
full tree: 500 internal cells, 504 terminal cells, four measured roots, and
1,000 evaluated edges. P0, AB, and P1 remain construction anchors with zeroed
unmeasured state and do not enter the scored forest. The implementation writes
cell- and layer-level manifests before rendering.

Figure 7 optimizes eight asynchronous bottom-up contraction rounds containing
504, 230, 126, 68, 38, 19, 10, and 5 edges. These rounds partition the 1,000
evaluated edges exactly. They must not be described as uniform developmental
depths because tracked branches end at different lineage depths. Each round
uses a matched parent-slot permutation null and biological-parent edge
retention. The eight fronts are then summed at common weights to form the
aggregate full-tree front.

Figure 8 now compares two methods: exact assignment within contraction rounds
and the greedy degree-constrained spanning forest. This emphasizes relaxation
of assignment constraints, with measured states, four roots and 1,000 edges
validated for both. Curves are sampled within-method fronts, not a global
optimum. First-cousin shuffle, the separate-clock Gaussian reference and one
broad random-rebuild null share continuous axes without inset or breaks.
Random rebuild preserves roots and nearly overlaps full-assignment shuffle;
the latter and the distinct internal-layer shuffle move to the supplement.

The supplementary figure (`full_tree_pareto/figs_ce_full_tree_heuristics.tex`)
contains five curves and all five null/reference families. The inventory table
(`full_tree_pareto/table_ce_full_tree_heuristics.tex`) lists all six methods,
including terminal-only, which is held out of all plots. Code inspection and
a controlled bookkeeping diagnostic found constructed internal spatial
midpoints, unfixed measured roots, and omission of the last pair (998 scored
edges) in its current generator. The saved optimized curve was not regenerated.
Earlier statements that all six methods shared validated scope were too strong;
moving a method to the supplement does not repair its scope. The table
distinguishes pipeline-validated main methods, code-inspected comparisons and
unverified historical caches. Six-method collective-nondominance flags remain
only as historical cache metadata.

The degree-constrained spanning forest and top-down rebuild trace nearly the
same front. Both are greedy reconstructions over the same combined cost matrix
`alpha*xyz + (1-alpha)*exp` with the same structural constraints, so they
differ by only a few cost units at matched weights and overlap on the main
axes.

The audit rejected and removed the notebook's historical MST cache. The saved
cache contains full internal-plus-terminal costs, but its unconstrained graph
forms one unrooted component, does not preserve the four biological roots, and
allows more than two children per parent. Separately, the current legacy
`mst_rebuild` implementation has an early-return bug and reports terminal costs
only, so it cannot regenerate that cache consistently. Figure 8 replaces it
with a four-root, degree-constrained spanning forest that validates 496 internal
edges, 504 terminal edges, and at most two children per parent at every sampled
weight. Its publication cache is `ce_full_tree_spanning_forest.csv`; do not
reintroduce `mst_rebuild.npz` into a publication comparison.

The audit also rejected `phylo_bm.npz`. It belongs to the older 10-PC/raw-cost
experiment and is on the scale visible in the earlier exploratory figure, not
the current top-20 z-scored protein analysis. The old implementation also used
tree IDs as measurement-row IDs—wrong for all 504 terminals here—and sampled
ancestral nodes independently, breaking Brownian parent--child covariance. The
initial replacement mapped biological identities explicitly and fitted a full
23-by-23 covariance to tracking-time-standardized increments. On 2026-09-09,
following the clock audit and author approval, the official Figure 8 reference
changed to separate clocks: spatial covariance per tracking-time unit and
protein covariance per canonical lineage transition, with full covariance
within blocks and zero cross-block covariance. All 1,000 edges span one
canonical generation, so this adds no fitted clock parameter. It
fixes the four optimization-root states and simulates every scored descendant
down the measured topology, producing 10,000 deterministic draws and a marked
1,000-draw display subset. This is a geometric reference for asking where
this simple Gaussian process on the fixed topology falls relative to the natural
lineage and sampled Pareto envelope. Its assumptions are not strongly grounded
in embryogenesis, so it must not be presented as a realistic biological null
or used for strong process-level inference.

The rationale is measurement-clock compatibility. Protein atlas aggregates
are not measured at the cutoff-truncated tracking endpoints; 49 of the top
50 former protein rate outliers end at the tracking cutoff, with median
duration 5. Dividing these differences by short tracking intervals inflates
the common protein rate. Adoption removes this unsupported scaling rather
than tuning the model to match the observed total.

The adopted 10,000-draw reference (seed 242) averages approximately 4,832
travel and 3,604 cell-state cost, versus 4,466 and 3,102 observed. The former
shared-time reference averaged 4,831 and 5,241; it is retained only as a
historical sensitivity comparison. The new model still misses observed
protein tail concentration and internal/terminal differences. Its displacement
must not be interpreted as a biological-efficiency effect size. Further
heterogeneity models require diagnostics and held-out evaluation. The decision,
diagnostics and regeneration commands are recorded
in `full_tree_pareto/README.md`, and the authoritative specification is
`full_tree_pareto/methods/separate_clock_reference.tex`. The Figure 8 caption
records both the cutoff rationale and remaining limitations. Publication
covariance and diagnostics now use `ce_full_tree_separate_clock_*.csv`.

Final reference review (2026-09-16): retain the **root-fixed, leaf-unconstrained
separate-clock Gaussian** in the manuscript. Its covariance is estimated from
all 1,000 internal and terminal edges, then used for forward simulation of all
descendants, without requiring the observed leaf outcomes. Two side analyses
are preserved in `full_tree_pareto/BRANCH_VARIABILITY_SENSITIVITY.md` and
`full_tree_pareto/LEAF_CONDITIONED_REFERENCE.md`, with their generators and tests.
The lognormal branch multiplier moves protein cost to about 2,988 but still
misses seven of eight held-out subtree total intervals. Jointly fixing root
and leaf states instead gives travel about 4,066 and protein about 2,936,
while permitting unconstrained internal configurations. Leaf conditioning is
mathematically valid with an all-edge covariance fit, but answers a different
endpoint-constrained question. Neither closer totals nor these lower costs
justify model adoption or biological-efficiency claims. Neither alternative
is an input to the main or supplementary manuscript figures; retain the
existing reference label, numerical cloud and plotting scales.

The next immediate step is author review of the candidate scientific story,
panel hierarchy, captions, and numbering. If the figures are accepted, update
the manuscript integration and treat the validated CSV caches and standalone
wrappers as the release path. Preserve the separation between full-tree and
terminal-only pipelines.

## Later work

Resolve the Figure 1 amendment with the co-author and revise the manuscript
for the pooled terminal release. Paused molecular/cross-species comparisons
remain archived exploratory work; consult the terminal handoff before resuming.
No terminal reorganization changes the full-tree analyses above.
