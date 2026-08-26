# Publication-figure handoff

Updated 2026-08-26. This repository-level document hands completed publication
work from `terminal_pareto/` to the next analysis in `full_tree_pareto/`. The
older plan in `task.md` predates the accepted terminal-cell figure sequence;
use the deliverable names and numbering below unless the authors explicitly
revise them. The manuscript draft `paper/paperv2.tex` is reference context and
should not be edited as part of a figure-only task.

## Current status

The *C. elegans* terminal-cell result is represented by Figures 2--6, Figures
S1--S3, and Table 1. All selected figures have dedicated renderers, split panel
PDFs where appropriate, standalone LaTeX wrappers, and compiled one-page PDFs.
The current versions have passed editorial and visual review and should be
treated as release candidates rather than exploratory layouts.

| Item | Role | TeX/PDF stem | LaTeX label |
|---|---|---|---|
| Figure 2 | Terminal Pareto front and null models | `fig2_ce_terminal_pareto_main` | `ce_terminal_pareto_main` |
| Figure 3 | Null schematic and structural support | `fig3_ce_terminal_pareto_supporting` | `ce_terminal_pareto_supporting` |
| Figure 4 | Investigated-subtree lineage map | `fig4_ce_subtree_map` | `ce_subtree_map` |
| Figure 5 | Canonical subtree proximity summary | `fig5_ce_canonical_summary` | `ce_canonical_summary` |
| Figure 6 | Cell-type retention and restricted fronts | `fig6_ce_cell_types` | `ce_cell_types` |
| Figure S1 | Cousin-null-relative sensitivity | `figS1_ce_canonical_summary_cousin_r` | `ce_canonical_summary_cousin_r` |
| Figure S2 | Cell-type endpoint cost ledger | `figS2_ce_cell_type_cost_gain` | `ce_cell_type_cost_gain` |
| Figure S3 | Three-tracking-replicate robustness | `figS3_ce_tracking_robustness` | `ce_tracking_robustness` |
| Table 1 | Per-subtree statistics | `table1_ce_subtree_statistics` | `ce_subtree_statistics` |

The authoritative scientific and production record for the complete
terminal-cell analysis is `terminal_pareto/README.md`.

## Fixed Figures 2--6 analysis configuration

Figure S3 instead applies the matched three-replicate configuration documented
in the Figure S3 section of `terminal_pareto/README.md`.

- Organism/data: *C. elegans* embryo 1, 3D tracking through `T <= 255`, and
  the top-20 protein features.
- Scope: terminal cells only. Alternative assignments change direct terminal
  parentage while terminal identities, positions, and molecular states remain
  fixed.
- Objectives: travel distance and cell-state distance. Never add or describe
  them as a single biologically meaningful efficiency score.
- Main null: first-cousin assignment shuffling, seed 42. Figures 2--3 also use
  second-cousin, third-cousin, and full-random shuffles.
- Subtree set: at least 12 usable terminal cells, yielding 43 nested subtrees.
- Pareto sweep: 300 weights for publication output.
- Terminal display coordinates in Figures 2--3: each objective is scaled by
  its first-cousin-null standard deviation and translated so the natural
  lineage is at `(0, 0)`.

Subtree null means, variances, and covariance are computed in closed form.
Random draws are used only for displayed null clouds and for optimality
frequencies when exact enumeration is infeasible. Nested subtrees are not
independent biological replicates, so ordinary regression p-values across the
43 points are inappropriate.

## Canonical subtree metrics

For each subtree, travel and cell-state objectives are separately normalized
by their spans between the travel-minimizing assignment `T*` and the
cell-state-minimizing assignment `S*`. The primary proximity statistic is

`d_LP = min_{p in sampled Pareto front} ||D(L) - D(p)||_2`,

where `L` is the natural assignment. The closest point must be an actual
sampled assignment; an interpolated segment is not an attainable lineage.
`u_L` is normalized arc length from the travel optimum (`u=0`) to the
cell-state optimum (`u=1`) at that same closest assignment.

`d_NP` is the endpoint-normalized distance from the first-cousin null mean to
the same selected Pareto assignment. The older cousin-relative statistic `r`
uses separate null-standardized geometry. Therefore `r` is a sensitivity
measure, not a scale factor, and `d_NP * r` must not be asserted to equal
`d_LP`. Both metrics remain in Table 1; Figure 5 uses `d_LP`/`d_NP`, while
Figure S1 shows `r`.

## Figure-specific decisions

- Figure 2 colors front assignments by natural-edge retention and shows four
  increasingly permissive null models; the full-random cloud remains in an
  inset because it lies far from the biologically relevant range.
- Figure 3B uses mean lineage-tree distance among changed edges. The all-edge
  definition remains a selectable diagnostic mode, not a second publication
  panel. Figure 3C reports the structural price of the final 5% of attainable
  reduction toward each single-objective endpoint.
- Figure 4 maps all 43 qualifying subtrees. Pie area reports cell-type
  composition, radius scales with `sqrt(n)`, circumference style reports
  whether the natural assignment is on the sampled front, and label color
  reports maximum edge retention.
- Figure 5A defines `u`, `d_LP`, and `d_NP`; Figure 5B directly labels the five
  major subtrees; Figure 5C shows all subtrees with major-subtree `d_NP`
  references. Larger subtrees receive higher z-order so major lineages are not
  hidden by smaller overlapping markers.
- Figure 6A decomposes edge retention along one unrestricted global front.
  Figures 6B--C instead optimize within cell types by forbidding assignments
  across types. The first-cousin null supplies common scaling there; the
  assignments themselves are type-restricted, not cousin-restricted.

## Shared visual semantics

- sequential blue: natural-edge retention;
- teal: travel direction or travel-minimizing assignment;
- vermillion: cell-state direction or cell-state-minimizing assignment;
- charcoal/gray: canonical position `u` and tree geometry;
- purple: natural-lineage-to-front distance `d_LP`;
- gold: first-cousin null and `d_NP`;
- warm earth-tone sequence: ordered null-model family; and
- local categorical colors: cell types, always accompanied by a legend.

Use sentence case for axes, keep the natural-lineage cross visually dominant,
label every displayed null in the legend, and use PDF as the canonical output.

## Source and cache map

Production scripts:

- `terminal_pareto/fig2_fig3_ce_terminal_pareto.py` -- Figures 2--3;
- `terminal_pareto/subtree_analysis.py` -- reusable subtree summaries/nulls;
- `terminal_pareto/fig4_ce_subtree_map.py` -- Figure 4;
- `terminal_pareto/fig5_table1_ce_canonical_metrics.py` -- Figure 5 metrics,
  validation, and Table 1;
- `terminal_pareto/fig5_figs1_ce_canonical_summary.py` -- Figure 5 and S1
  component panels;
- `terminal_pareto/subtree_explore.py` -- cell-type decomposition helpers
  retained under its historical filename;
- `terminal_pareto/fig6a_figs2_ce_cell_types.py` -- Figure 6A and Figure S2;
- `terminal_pareto/fig6bc_ce_within_type.py` -- Figure 6B--C;
- `terminal_pareto/figS3_ce_tracking_robustness.py` -- matched three-replicate
  audit, analysis, and Figure S3; and
- `terminal_pareto/plot_style.py` -- shared visual semantics and export style.

Authoritative caches in `terminal_pareto/output/`:

- `subtree_summary_min12.csv`;
- `ce_subtree_canonical_metrics.csv`;
- `ce_subtree_canonical_curves.npz`;
- `global_front_retention_by_cell_type.csv`;
- `global_front_by_cell_type.csv`;
- `within_type_fronts.csv`;
- `within_type_summary.csv`;
- `type_preserving_aggregate_front.csv`;
- `ce_tracking_replicate_audit.csv`;
- `ce_tracking_replicate_cell_manifest.csv`;
- `ce_tracking_replicate_subtree_manifest.csv`;
- `ce_tracking_replicate_fronts.csv`;
- `ce_tracking_replicate_null_clouds.csv`; and
- `ce_tracking_replicate_major_metrics.csv`.

`terminal_pareto/output/publication/` should contain only accepted wrappers,
their directly included panel assets, and compiled manuscript figures. Do not
restore combined pre-split canvases, abandoned canonical-map prototypes,
exploratory subtree plots, or minimum-cell threshold test caches unless a new
analysis question explicitly needs them. Generated output is ignored by Git;
publication TeX files require `git add -f` when they are intentionally staged.

## Regeneration

From the repository root using the `dev` Conda environment:

```bash
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig2_fig3_ce_terminal_pareto.py
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/subtree_analysis.py --min-cells 12
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig4_ce_subtree_map.py --min-cells 12
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig5_table1_ce_canonical_metrics.py --min-cells 12
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig5_figs1_ce_canonical_summary.py
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig6a_figs2_ce_cell_types.py --iteration 300
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/fig6bc_ce_within_type.py --iteration 300
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/figS3_ce_tracking_robustness.py
```

Compile the wrappers from `terminal_pareto/output/publication/`:

```bash
for stem in \
  fig2_ce_terminal_pareto_main \
  fig3_ce_terminal_pareto_supporting \
  fig4_ce_subtree_map \
  fig5_ce_canonical_summary \
  fig6_ce_cell_types \
  figS1_ce_canonical_summary_cousin_r \
  figS2_ce_cell_type_cost_gain \
  figS3_ce_tracking_robustness \
  table1_ce_subtree_statistics
do
  conda run -n dev tectonic "${stem}.tex"
done
```

## Validation before handoff or commit

1. Run `python -m py_compile` on every changed Python source.
2. Run the canonical self-test:
   `conda run -n dev python terminal_pareto/fig5_table1_ce_canonical_metrics.py --min-cells 12 --selftest`.
3. Regenerate any renderer whose code changed and inspect its PNG/PDF output.
4. Compile every affected standalone wrapper with Tectonic.
5. Run `git diff --check` and inspect `git status --short`.
6. Confirm no unrelated worktree changes or expensive full-tree caches were
   staged.

## Tracking-replicate robustness (completed)

Figure S3 implements the planned three-replicate comparison in a dedicated
renderer without changing the accepted Figures 2--6 pipeline. The predefined
stage-matched cutoffs are embryo 1 `T=255`, embryo 2 `T=247`, and embryo 3
`T=226`; the numerical times differ because developmental rates differ. They
provide 299, 301, and 309 usable terminal edges, respectively. The strict
three-way intersection contains 279 edges, retaining 90.3--93.3% per embryo,
so all displayed comparisons use the same terminal identities and natural
parents. Forty-two subtrees remain eligible at `n >= 12`.

Panel A shows the matched P0 fronts and per-replicate first-cousin nulls in
coordinated null-SD axes. Panel B reports `u_L`, `d_LP`, and maximum edge
retention for P0, AB, ABa, ABp, and P1. The global fronts are closely aligned,
and P0 variation is modest (`u_L` range 0.019, `d_LP` range 0.012, maximum
retention range 0.057). ABa shows genuine local replicate dependence: embryo
3 places the natural assignment on the sampled front and reaches maximum
retention 1.0. This mixture of global robustness and local sensitivity is
better presented as a supplementary robustness result than as a new main-text
claim. See `terminal_pareto/README.md` for the complete matching decision,
analysis record, and regeneration workflow.

## Next work: full-tree Pareto analysis

The immediate next publication task moves to `full_tree_pareto/`. It will test
Pareto organization across the complete embryonic lineage, including internal
and terminal cells, rather than extending the terminal-only assignment model.
Begin with the *C. elegans* protein configuration and optimize the full tree
layer by layer. The resulting analysis is intended to support the next main
publication figure.

Start with a code-and-data audit before designing panels:

1. Inventory the existing full-tree runners, optimization variables, cached
   `embryo1_prot` results, and legacy assumptions in
   `full_tree_pareto/full_tree_pareto.ipynb`.
2. Define precisely what one optimization layer contains: the cells being
   reassigned, their admissible parents, the fixed state inherited from earlier
   layers, and how a layer contributes to whole-tree travel and cell-state
   costs.
3. Build a layer manifest containing cell identities, lineage depths,
   measurement availability, candidate-parent counts, and exclusions. Validate
   that each natural edge is admissible before running optimization.
4. Establish a reproducible *C. elegans* protein baseline, then generate and
   validate the layer-by-layer Pareto results before selecting the scientific
   story or figure layout.

Carry forward the norms that generalize from `terminal_pareto/`:

- keep travel and cell-state objectives separate and avoid describing a linear
  combination as a biological efficiency score;
- record the exact measurement space, optimization scope, null model, random
  seed, sweep resolution, and cache schema;
- distinguish attainable sampled assignments from interpolated display
  geometry;
- retain natural-lineage, objective-endpoint, and structural-retention
  landmarks where their full-tree definitions remain valid;
- use deterministic validation, manuscript color semantics, sentence-case
  labels, canonical PDF panels, and standalone LaTeX assembly; and
- place expensive generated results under `full_tree_pareto/output/` without
  coupling them to or modifying the accepted terminal-only pipeline.

Do not copy terminal-specific normalization, cousin-shuffle constraints, or
edge-retention definitions without checking that they remain meaningful for a
sequential full-tree reconstruction. Document any full-tree replacement next
to the implementation. The existing `full_tree_pareto/README.md` describes the
current runners and cache layout; update it as the audited design becomes the
new authoritative record for that module.

## Later work

After the full-tree *C. elegans* protein figure, later candidates include the
held Figure 1 overview, additional molecular configurations, and matched
*C. briggsae* comparisons. These analyses should reuse validated caches and
matched measurement spaces and should not be silently folded into the accepted
terminal-only pipeline.
