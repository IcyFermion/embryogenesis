# Publication-figure handoff

> **Layout change, 2026-10-04.** Figures, captions and releases moved to
> [`publication/`](publication/README.md); production is `publication/output/production/`.
> This package is numerical only. Figure material that this record places in
> run folders (`figures/`, `publication/`, `layout_history/`) or in the old
> `output/publication/` is archived in
> `publication/output/archive/migrated_20261004/`; figure history formerly in
> `output/legacy/` (release archives, embryo-1 figures, diagnostics) is in
> `publication/output/archive/legacy_20261004/` at the same relative paths. Commands named below that drew or released figures are retired;
> see the package README for current commands.

Updated 2026-10-04.

All figures, captions and Table 1 are built and released by the shared
presentation package; see [publication/README.md](publication/README.md).
Production is `publication/output/production/` (released 2026-10-04 with
`d_CP`/`C` notation and one universal style). Older paths below are historical.

## Terminal figures: pooled release promoted

The author approved the pooled-travel, endpoint-normalized terminal figures.
The current release is `terminal_pareto/output/publication/`, including the
standalone Figure 1 amendment, single-panel Figure 2, and four-geometry S3
with insets, now extended with cross-species Figures 8 and 9. Previous single-embryo figures are in
`terminal_pareto/output/legacy/embryo1/publication/`.

Use [the terminal README](terminal_pareto/README.md) for analysis conventions,
figure paths, and build commands, and [the terminal handoff](terminal_pareto/HANDOFF.md)
for current validation, short development history, archived pilots/prototypes,
and remaining manuscript/co-author integration. These supersede the older
terminal sections of this repository-level document. The manuscript remains
unchanged by publication promotion. Full-tree work below remains separate.

## Cross-species terminal Figures 8 and 9

On 2026-09-30, the author assigned Figure 8 to the pooled 187-edge six-panel
comparison (3D above XY; CE protein, CE RNA and CB RNA columns) and Figure 9 to
the separate paired 3D/XY front overlays, now supplemented by canonical-metric
panel C in the Figure 5B named-row style. It shows u, d_LP and d_NP for all
three configurations under each geometry; distances use the same closest
attained front point P*, not maximum retention. Figure 9A--B marks P* with
hollow diamonds and natural-to-P* dashed connectors. Figure 8 has no
natural-to-maximum-retention connectors. Both are assembled with captions
in the earlier terminal publication-wrapper style, compiled to one page and
visually checked. On 2026-10-01 their editable TeX, PDF pages and panel assets
were additively published in `terminal_pareto/output/publication/`:

- `fig8_terminal_cross_species_comparison.pdf`
- `fig9_terminal_cross_species_overlays.pdf`

The C. briggsae 3D tracking-method caveat is caption-only and provisional
pending experimental collaborator input. Numerical caches are unchanged;
figure numbers and source/artifact/analysis hashes are recorded in the build
manifest. The source assembly remains under
`terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930/publication/`.
The main release manifest now verifies 66 assets and separately records the
187-edge comparison cohort/provenance alongside the primary 275-edge analysis.
All 57 earlier production assets and nine comparison numerical files remain
unchanged. The preceding release is recoverable with verified hashes under
`terminal_pareto/output/legacy/releases/20261001T154243342749Z/publication/`.
Use `python terminal_pareto/publication_release.py --cross-species` to repeat
the explicit additive release and `--verify` to check production. Future
pooled-layout releases retain these comparison assets automatically; staged
failure rollback and retention are tested. See the terminal README/handoff
for reproduction and validation. The manuscript remains unchanged.
Older full-tree Figure 8 references below are historical;
the current full-tree specification includes Figure 7, its supplement and the
partial-forest cross-species Figures 10--11.

## Cross-species partial-forest Figures 10 and 11

A separate cross-species **partial-forest layerwise analysis** is now complete.
See [its handoff](full_tree_pareto/CROSS_SPECIES_HANDOFF.md). Starting from the
same 187 matched terminal edges as Figures 8--9, it ascends only through
available canonical ancestors, retaining 485 measured cells and 454 edges
with 31 fixed boundary roots. It has no imputation, skipped ancestors, other
reconstruction heuristics or Gaussian reference. This is not a complete tree.
New pooled-travel, Euclidean molecular-distance and endpoint-normalized caches
are isolated under
`full_tree_pareto/output/runs/cross_species_layerwise_v1/terminal_anchored_20260930/`.
The captioned six-panel comparison and overlays/canonical-metric companion
follow Figures 8/9's layouts and Figure 7B's reference inventory. The author
approved and released them as **Figures 10 and 11** on 2026-10-01, with
301-weight displays and separately cached
1,201-weight resolution checks; CE RNA XY has a modest closest-point
sensitivity recorded in the handoff. Captioned one-page PDFs, editable TeX
and panel PDF/PNG assets are in `full_tree_pareto/output/publication/`:

- `fig10_full_tree_cross_species_comparison.pdf`
- `fig11_full_tree_cross_species_overlays.pdf`

The additive release preserved all 15 unrelated production files, including
Figure 7 and its supplement, and all 29 partial-forest numerical artifacts.
Only the four obsolete `fig8_ce_full_tree_collective` PDF/TeX/panel PDF/panel PNG
files were removed from production; they remain recoverable in the hash-checked
pre-release archive
`full_tree_pareto/output/legacy/cross_species_releases/20261001T150611720032Z/publication/`.
Production assets, preserved files and archive hashes are recorded in
`cross_species_release_manifest.json`; verify with
`python -m full_tree_pareto.cross_species_publication --verify-production`.
Both PDFs were visually checked; all 54 full-tree and 13 terminal cross-species
tests passed. No pooled Figure 7 promotion or manuscript integration occurred.

## Full-tree Pareto analysis: pooled Figure 7 and supplement

**Released 2026-10-04:** pooled Figure 7 and provisional S4 are in
`full_tree_pareto/output/publication/` beside Figures 10/11, via the staged
mixed release in `publication/release.py`; the former directory is archived
under `full_tree_pareto/output/legacy/releases/20261004T024949483965Z/`. See
`full_tree_pareto/POOLED_HANDOFF.md`.

The author requested one main full-tree figure and one supplementary comparison,
replacing the two-main-figure arrangement. Current scientific specification and
build commands are in [the full-tree README](full_tree_pareto/README.md).
[The resumable handoff](full_tree_pareto/POOLED_HANDOFF.md) records live build
status, validated caches, tests and remaining checks; consult it before resuming.

- Figure 7A: eight bottom-up layerwise fronts. Figure 7B: aggregate layerwise
  front, first-/second-/third-cousin shuffles, with random rebuild in an inset
  matching the terminal layout. Gaussian was removed from the main panel at
  the author's request on 2026-09-29; its phylogenetic analogy is nonessential
  and its distinct internal-state inventory limits comparison to the fronts.
- Supplement: five regenerated heuristics, seven reference families (including
  Gaussian) and a
  constraint inventory. Terminal-only remains withheld.
- Travel: mean normalized pairwise distances from embryos 1--3 at cutoffs
  255/247/225, on 978 shared measured nodes and 974 full-tree edges. Each
  embryo uses its natural full-tree total as the global denominator.
- Display: terminal-style endpoint normalization, with each round's anchors
  in A and shared aggregate layerwise anchors for B and the supplement.
- Cousin shuffles: permute measured leaves (including cutoff/coverage leaves)
  within canonical ancestor groups two/three/four generations back; position
  and protein state move together. Internal states and four roots stay fixed;
  all 974 edges are scored. New 10,000-draw second-/third-cousin caches are
  independent add-ons; the expensive 301-weight reconstruction caches remain
  untouched. Existing 1,505 forests and all 20,000 new draws are replay-validated.
- Gaussian (supplement only): the author explicitly confirmed free forward internal simulation
  followed by observed terminal-state clamping, **not** joint conditioning.
  All measured leaves, including cutoff/coverage leaves, and four roots are
  fixed. Spatial processes are separate per embryo; protein uses a shared
  per-transition clock. Terminal reconnection costs are recomputed.
- Protein distance remains Euclidean top-20 z-scored features, unlike terminal
  cosine. This distinction and unequal heuristic/reference feasible sets are
  retained explicitly; no biological efficiency or calibrated-null claim.

The new build is isolated in
`full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_20260927/`.
Publication promotion preserves the former bundle under a hash-checked archive.
The old pipeline remains available and is not a source of pooled caches.
Manuscript integration is a separate step.

The author approved the revised Figure 7 and supplement and requested a
repository checkpoint on 2026-09-30. Source, caption templates, tests and
documentation are included; generated artifacts remain Git-ignored. This
approval does not itself replace `full_tree_pareto/output/publication/`.

## Historical full-tree decisions: pre-pooled Figures 7--8

**The remainder of this section records the former single-embryo release. It is
not the current figure specification.** In particular, the old decision to use
leaf-unconstrained simulation was superseded by the explicit terminal-clamping
request above. Old numerical totals refer to 1,000 edges and must not be mixed
with the pooled cohort. See also `full_tree_pareto/LEGACY_EMBRYO1.md`.

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

### Backlog: restore and validate 1,001-weight Pareto sweeps

Recorded with author approval on 2026-09-28; **deferred, not implemented**.
The original exploratory implementations use 1,000 intervals and 1,001
endpoint-inclusive weights. Current terminal publication configurations and
the new pooled full-tree build instead use 300 intervals and 301 weights.
Terminal publication code already specified 300 intervals in `ae557ab`;
full-tree publication code introduced that setting in `d493499`, and the
pooled terminal release `3b6c31f` retained it. The reason for the reduction is
not established; do not attribute it to runtime savings as a verified fact.

- Audit sweep settings across terminal global/subtree/cell-type/tracking
  analyses and all full-tree heuristics; restore the intended 1,001-weight
  convention consistently, including builders, cache identities and captions.
- Preserve the 301-weight caches and figures as a resolution baseline. Use
  new run identities for the denser sweeps, with per-weight checkpoints for
  expensive full-tree reconstruction. The 301- and 1,001-point grids are not
  nested, so do not assume all existing assignments can be reused.
- Compare attained fronts, retention and reported front-dependent metrics;
  verify endpoint anchors and raw-cost replay before adopting new figures.
  A denser sweep still does not enumerate the complete discrete Pareto set.
- Null distributions do not depend on sweep resolution. Reuse existing draws
  only after verifying unchanged cohort, objectives, reference model, seeds
  and source provenance; do not rewrite manifests to bypass cache checks.
- Update documentation and explicitly promote revised figures only after
  validation/review. This backlog entry changes no algorithms, caches,
  figures or publication assets.

### Other deferred work

Resolve the Figure 1 amendment with the co-author and revise the manuscript
for the pooled terminal release. The original molecular/cross-species pilot
remains archived. On 2026-09-30, the author requested a resumed terminal-only
CE protein / CE RNA / CB RNA comparison in Figure 2 style, subsequently
numbered Figures 8 and 9 and assembled with captions as recorded above.
Its pooled 187-edge numerical checkpoint and build are documented in
`terminal_pareto/README.md` and `HANDOFF.md`, under
`terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930/`.
No terminal reorganization changes the full-tree analyses above.
