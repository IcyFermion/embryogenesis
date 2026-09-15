# Terminal-cell lineage Pareto analysis

Updated 2026-08-26. This directory contains the complete analysis and
publication pipeline for the terminal-cell portion of the project. This README
is the authoritative module record: it combines the analysis conventions,
figure decisions, reproducible findings, caveats, and production instructions
for Figures 2--6, Figures S1--S3, and Table 1.

The analysis tests whether terminal-cell parentage in nematode embryogenesis is
Pareto optimal with respect to two competing objectives:

- **travel distance**, derived from embryo cell tracking; and
- **cell-state distance**, derived from protein or RNA expression profiles.

Terminal-cell identities, final positions, and molecular states are held fixed.
Alternative lineages are reconstructed by assigning terminal cells to candidate
parents with minimum-cost bipartite matching across a sweep of objective
weights. The primary publication analysis uses *C. elegans* protein expression
and embryo-1 tracking. Comparative analyses cover RNA expression, tracking
replicates, subtrees, two-dimensional tracking, z-noise controls, edge
perturbations, and random-feature controls.

## Publication status

The selected manuscript outputs are release candidates. They have dedicated
renderers, split panel assets where appropriate, standalone LaTeX wrappers, and
compiled one-page PDFs.

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

The repository-level `PUBLICATION_FIGURES_HANDOFF.md` records the transition
from this completed module to the next publication analysis.

## Fixed analysis conventions

The accepted Figures 2--6 configuration is:

- *C. elegans* embryo 1, three-dimensional tracking through `T <= 255`;
- top-20 protein features with cosine cell-state distance;
- terminal cells only;
- first-cousin assignment shuffling with seed 42 as the main null;
- at least 12 usable biological terminal cells per subtree, giving 43 nested
  investigated subtrees; and
- a 300-weight Pareto sweep for the established publication pipeline.

Figures 2--3 additionally show second-cousin, third-cousin, and full-random
assignment shuffles. The superseded 2014 *C. briggsae* tracking series
(`1407...`) is excluded from final replicate comparisons.

Travel and cell-state costs are reported in standard-deviation units from the
first-cousin-shuffle null distribution where a common embryo-wide scale is
needed. In Figures 2--3, displayed coordinates are translated so the natural
lineage is at `(0, 0)`. Translation changes only the presentation origin; it
does not change assignments, fronts, rankings, or relative-Pareto statistics.

Subtree null means, variances, and covariance are computed in closed form.
Random draws are used only for displayed null clouds and for optimality
frequencies when exact enumeration is infeasible. A non-positive null standard
deviation raises before division.

## Canonical subtree metrics

For each subtree, let `T*` be its travel-minimizing assignment, `S*` its
cell-state-minimizing assignment, and `L` its natural assignment. The two
objectives are normalized separately using only the endpoint spans:

\[
D_1(a)=\frac{T(a)-T(T^*)}{T(S^*)-T(T^*)},\qquad
D_2(a)=\frac{S(a)-S(S^*)}{S(T^*)-S(S^*)}.
\]

The natural-lineage distance to the sampled front is

\[
d_{LP}(L)=\min_{p\in\mathcal P_{\mathrm{sampled}}}
\left\|D(L)-D(p)\right\|_2.
\]

The minimum is evaluated over actual sampled assignments, never an interpolated
segment. Coordinates are not clipped to the endpoint box. The canonical
position `u_L` is normalized arc length along the front, from the travel
optimum (`u=0`) to the cell-state optimum (`u=1`), evaluated at the same nearest
assignment.

The endpoint-normalized distance from the first-cousin null mean `N` to that
selected assignment `P*` is

\[
d_{NP}=\left\|D(N)-D(P^*)\right\|_2.
\]

The older cousin-relative statistic `r` uses separately null-standardized
geometry. It is a sensitivity measure, not a scale factor: `d_NP * r` must not
be asserted to equal `d_LP`. Table 1 retains both metrics; Figure 5 uses
`d_LP`/`d_NP`, while Figure S1 shows `r`.

## Figure records

### Figures 2--3: global terminal-cell analysis

Figure 2 shows the global Pareto front, natural assignment, and increasingly
permissive null models. The full-random cloud remains in an inset because it
lies far from the biologically relevant range.

Figure 3B uses mean lineage-tree distance among changed edges. The all-edge
definition remains a diagnostic mode, not a publication panel. Figure 3C
reports the structural price of pursuing the last 5% of attainable distance
reduction from the maximum-edge-retention compromise toward either
single-objective endpoint.

Principal structural metrics are:

- **edge retention:** fraction of natural terminal parent--child edges
  preserved;
- **mean lineage-tree distance:** mean tree separation between natural and
  reconstructed parent assignments;
- **TWER:** tree-distance-weighted edge retention; and
- **local perturbation tests:** pair and triple edge swaps that test whether
  both objectives can be improved near the natural lineage.

### Figure 4: investigated-subtree map

Figure 4 places all 43 qualifying subtrees in lineage order, with P0 at the
top. Each pie summarizes the fate composition of the exact usable terminal set
used in that subtree's Pareto analysis. Radius scales with `sqrt(n)`. A solid
circumference means the natural assignment occurs on the sampled front; a
dashed circumference means it does not. Edges connect each displayed subtree
to its nearest qualifying descendant, and non-qualifying nodes are omitted.
Subtree names and cell counts use sequential blue to encode maximum natural-edge
retention. The accepted layout has no measured label collisions.

### Figure 5, Figure S1, and Table 1: canonical subtree summary

Figure 5A defines `u`, `d_LP`, and `d_NP`. Figure 5B reports `d_LP` and `d_NP`
for P0, AB, ABa, ABp, and P1. Figure 5C places all 43 subtrees in `(u, d_LP)`
space and overlays the five major-subtree `d_NP` references. Larger subtrees
receive higher z-order so major lineages are not hidden by smaller markers.
Maximum edge retention is deliberately omitted so each marker communicates
canonical position and proximity only; retention remains in Figures 2--4 and
6 and in Table 1.

Figure S1 preserves cousin-relative `r` as a null-model sensitivity analysis.
Across the 43 subtrees, `d_LP` and `r` have Spearman correlation 0.95, so the
broad proximity ranking does not depend on the first-cousin normalization.

### Figure 6 and Figure S2: terminal cell types

Cell types are loaded from `data/2023-06-29_entropy_cell_key_V2.csv`, joined on
`wormweb.lineage`, and consolidated with `MERGE_MAP`. The 32 terminal cells
absent from that key exactly match `data/apoptotic_cells.txt` and are typed as
programmed death. The separately annotated dying tail hypodermal cells
`ABplppppppa` and `ABprppppppa` remain distinct and enter the displayed
`other (n<=4)` group.

Figure 6A shows natural-edge retention by terminal fate across all 145 distinct
assignments observed along the global Pareto sweep. Biological parent identity
defines retention, so swapping assignment slots belonging to the same parent
retains both natural edges. Weighted cell-type retention is asserted to
recompose the global curve exactly at every cached assignment. The accepted
maximum-retention assignment has global retention 0.819398. Irregular front
positions are resampled to a uniform canonical-position grid and narrowly
averaged for display only; exact assignments, keypoints, and statistics remain
unsmoothed.

Figures 6B--C repeat the optimization after forbidding assignments across
terminal cell-type groups. Six sufficiently large types receive individual
fronts in shared embryo-wide first-cousin-null standard-deviation units per
cell. This null supplies a common scale only; assignments are constrained by
cell type, not cousin relationships. The aggregate type-preserving front is
then compared with the unrestricted front. At the corresponding
single-objective endpoints, type restriction retains 3.101 standard deviations
more travel distance and 15.073 standard deviations more cell-state distance.
It realizes 16.8% of the unrestricted travel saving and 63.0% of the
unrestricted cell-state saving.

Figure S2 is the complementary endpoint ledger: signed per-cell travel and
cell-state changes at the travel optimum, maximum-retention compromise, and
cell-state optimum. The objectives remain on separate axes and are never
combined linearly.

### Figure S3: three-replicate tracking robustness

Figure S3 is generated independently of the accepted Figures 2--6 pipeline.
Numerical time points are not equated across embryos because developmental
timing differs. The analysis retains the predefined stage-matched cutoffs and
constructs a strict intersection of natural terminal parent--child edges. Both
the terminal and its natural parent must be tracked and have the fixed protein
representation.

| Replicate | Cutoff | Usable terminal edges | Strictly matched | Retained |
|---|---:|---:|---:|---:|
| embryo 1 | 255 | 299 | 275 | 92.0% |
| embryo 2 | 247 | 301 | 275 | 91.4% |
| embryo 3 | 225 | 298 | 275 | 92.3% |

The matched set leaves 42 qualifying subtrees at `n >= 12`. All replicates use
the same 275 edges, top-20 protein features, cosine cell-state distance,
first-cousin shuffle with seed 42, exact null moments, 1,000 displayed null
draws, and a 300-interval endpoint-inclusive sweep containing 301 weights.

Panel A overlays the matched P0 fronts in replicate-specific
first-cousin-null-SD coordinates, with each natural lineage translated to the
origin. Panel B reports `u_L`, `d_LP`, and maximum edge retention for P0, AB,
ABa, ABp, and P1. Equal-length horizontal marks are centered on each subtree's
categorical position, match the replicate colors in Panel A, and use uniform
transparency so coincident values visibly blend. Marker shape does not encode
replicate identity. `d_NP` remains cached but is omitted from the panel.

The fronts are closely aligned globally. For P0, replicate ranges are 0.009
for `u_L`, 0.016 for `d_LP`, and 0.044 for maximum edge retention. ABa shows
the largest local spread: 0.127 for `u_L`, 0.054 for `d_LP`, and 0.156 for
maximum edge retention. The largest `d_NP` spread is instead P1 at 0.042. The
result supports global tracking robustness while preserving genuine local
replicate dependence; it does not establish that every subtree is replicate
invariant.

The matching decision is recorded in:

- `output/ce_tracking_replicate_cell_manifest.csv`;
- `output/ce_tracking_replicate_audit.csv`; and
- `output/ce_tracking_replicate_subtree_manifest.csv`.

The remaining S3 caches are `ce_tracking_replicate_fronts.csv`,
`ce_tracking_replicate_null_clouds.csv`, and
`ce_tracking_replicate_major_metrics.csv`.

## Reproducible biological findings

1. The natural terminal assignment is exactly sampled-Pareto-optimal in 16 of
   43 investigated subtrees, all with at most 35 terminal cells. Exact recovery
   becomes less frequent among larger subtrees, but continuous size
   relationships are modest: size versus on-front status has Spearman rho
   `-0.35`, versus front proximity `+0.30`, versus cousin-relative distance
   `+0.18`, and versus maximum retention `-0.06`.
2. Exact recovery is not explained by small assignment spaces. Full-random
   assignments are on the front in 0/1000 draws for every subtree. The
   cousin-shuffle fraction is usually zero and reaches 33% only in the small,
   exactly enumerable MSpa assignment space.
3. Apparent relationships between attainable reductions and subtree depth or
   fate diversity disappear after per-cell normalization (`|rho| < 0.2`).
   Subtree size and median depth are strongly confounded (`rho = 0.90`).
4. True-sibling contrasts are at chance level: 8 of 14 informative pairs agree
   with the proposed direction, with three distance ties excluded.
5. Per-cell endpoint contributions differ by fate. Programmed-death cells have
   the largest cell-state improvement at the cell-state optimum
   (`-0.25 sigma/cell`); epithelium pays the largest travel price there
   (`+0.47 sigma/cell`); muscle pays the largest cell-state price at the travel
   optimum (`+0.22 sigma/cell`). Neurons dominate totals primarily because
   they are abundant.

## Interpretation and statistical caveats

- Nested subtrees are descriptive observations, not independent replicates.
  EMS/MS, ABplpp/ABplppp, and ABprpp/ABprppp can share identical cousin-group
  terminal sets and therefore identical null moments.
- Do not report regression p-values that treat the 43 subtree points as
  independent.
- `d_LP` is an endpoint-normalized geometric display statistic, not a
  biological weighted sum of travel and cell-state distances.
- Cousin-relative `r` is a two-dimensional Euclidean distance in separately
  null-standardized axes, not a conventional one-dimensional z-score.
- Equal weighting of standardized objectives must not be described as a
  biologically balanced cost; the objectives are not linearly comparable.
- Within-type cell-state gains are not uniformly negligible. Attributing them
  to expression noise requires an independent replicate or stability analysis.

## Shared visual semantics

- sequential blue: natural-edge retention;
- teal: travel direction or travel-minimizing assignment;
- vermillion: cell-state direction or cell-state-minimizing assignment;
- charcoal or neutral gray: canonical position `u` and tree geometry;
- purple: natural-lineage-to-front distance `d_LP`;
- gold: first-cousin null and `d_NP`;
- warm earth-tone sequence: ordered null-model family; and
- local categorical colors: cell types or replicates, always identified in the
  panel or legend.

Use sentence case for axes, keep the natural-lineage cross visually dominant,
label every displayed null, and use PDF as the canonical output.

## Code and output layout

```text
terminal_pareto/
├── README.md                              # This authoritative module record
├── data_loader.py                         # Tracking, expression, and lineage I/O
├── pareto_engine.py                       # Assignment and null-model methods
├── lineage_metrics.py                     # Retention and tree-distance metrics
├── plot_style.py                          # Shared publication style
├── main.py                                # Complete Marimo analysis notebook
├── subtree_analysis.py                    # Per-subtree summaries and exact nulls
├── subtree_explore.py                     # Cell-type decomposition helpers
├── fig2_fig3_ce_terminal_pareto.py        # Figures 2--3
├── fig4_ce_subtree_map.py                 # Figure 4
├── fig5_table1_ce_canonical_metrics.py    # Figure 5 metrics and Table 1
├── fig5_figs1_ce_canonical_summary.py     # Figure 5 and Figure S1
├── fig6a_figs2_ce_cell_types.py           # Figure 6A and Figure S2
├── fig6bc_ce_within_type.py               # Figure 6B--C
├── figS3_ce_tracking_robustness.py        # Figure S3
└── output/
    ├── ce_protein/                        # Configuration diagnostics
    ├── ce_rna/
    ├── cb_rna/
    └── publication/                       # Panels, wrappers, and final PDFs
```

Generated output is ignored by Git. Publication TeX wrappers are intentionally
versioned because they contain captions and layout, so a new ignored wrapper
requires `git add -f`. Do not restore abandoned layouts, combined pre-split
canvases, or threshold-test caches unless a new analysis question needs them.

Accepted component stems in `output/publication/` are:

- Figure 2: `fig2_ce_terminal_pareto_front`;
- Figure 3: `fig3A_ce_null_models`,
  `fig3B_ce_edge_retention_tree_distance`, and
  `fig3C_ce_structural_retention`;
- Figure 4: `fig4_ce_subtree_map_panel`;
- Figure 5: `fig5A_ce_canonical_definition`,
  `fig5B_ce_canonical_major_subtrees`, and
  `fig5C_ce_canonical_all_subtrees`;
- Figure 6: `fig6A_ce_retention_heatmap`,
  `fig6B_ce_within_type_fronts`, and
  `fig6C_ce_type_restricted_aggregate`;
- Figure S1: `figS1A_ce_cousin_r_definition`,
  `figS1B_ce_cousin_r_major_subtrees`, and
  `figS1C_ce_cousin_r_all_subtrees`;
- Figure S2: `figS2_ce_cell_type_cost_gain_panel`; and
- Figure S3: `figS3A_ce_tracking_replicate_fronts` and
  `figS3B_ce_tracking_replicate_metrics`.

Authoritative numerical caches in `output/` are:

- `subtree_summary_min12.csv`;
- `ce_subtree_canonical_metrics.csv` and
  `ce_subtree_canonical_curves.npz`;
- `global_front_retention_by_cell_type.csv` and
  `global_front_by_cell_type.csv`;
- `within_type_fronts.csv`, `within_type_summary.csv`, and
  `type_preserving_aggregate_front.csv`; and
- the six `ce_tracking_replicate_*.csv` files listed in the Figure S3 section.

## Regeneration

Run from the repository root using the `dev` Conda environment:

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

For the complete notebook analysis, run `python terminal_pareto/main.py`; use
`marimo edit terminal_pareto/main.py` for an interactive session.

Compile wrappers from `terminal_pareto/output/publication/`:

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
3. Regenerate each changed renderer and inspect its PNG/PDF output.
4. Compile each affected standalone wrapper with Tectonic.
5. Run `git diff --check` and inspect `git status --short`.
6. Confirm no unrelated files or expensive full-tree caches are staged.

Python dependencies are defined in `environment.yml`. Tectonic compiles the
standalone publication figures.

## Cross-dataset terminal pilot (2026-09-09)

`pilot_cross_dataset.py` implements an unnumbered comparison of C. elegans
protein, C. elegans RNA, and two AF16 C. briggsae RNA/tracking configurations
on 188 shared terminal edges. See `PILOT_CROSS_DATASET.md` for findings,
endpoint tie handling, regeneration, and calibration/provenance limitations.
The four-page review and caches are in `output/pilot_cross_dataset/`.
This pilot includes existing-3D and XY-only views plus coverage and cell-set
sensitivities; it is not an accepted publication figure.

Paused by author request on 2026-09-14. The pilot record includes the review
clarifications (coverage is not edge retention), the CE depth diagnostic,
validation status, and an ordered resume checklist. Preserve the current
ignored output directory before rerunning; generated plots/caches are local,
not included in the source commit. No publication promotion is authorized.
