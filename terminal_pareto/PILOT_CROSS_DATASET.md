# Terminal cross-dataset pilot

Implemented 2026-09-09. This is an unnumbered exploratory figure, not a release
candidate. The dedicated script is `pilot_cross_dataset.py`; outputs are in
`output/pilot_cross_dataset/`. Start with `pilot_review.pdf` (four pages) or
the two standalone `comparison_existing_3d.pdf` and `comparison_xy.pdf` files.
`pilot_cross_dataset.tex` assembles a separate captioned two-page comparison.

**Status: paused by author request on 2026-09-14.** Keep the current plots and
numerical caches as the exploratory checkpoint; no new sweep, geometry
correction, or publication promotion is part of this pause. Accepted terminal
Figures 2--6 and the separate full-tree pipeline are not changed by this pilot.
Resume with the checklist below, not by rerunning the analysis immediately.

## Configuration and matching

The primary set contains 188 natural terminal parent--child edges shared by
C. elegans embryo 1 protein, C. elegans embryo 1 RNA, and two AF16 C. briggsae
tracking embryos. Both child and parent require tracking and expression.
Terminal identities, natural parent identities, and the available parent-slot
multiplicities are identical. Candidate assignments are unrestricted within
this common terminal set; cousin relationships constrain the null only.

| Configuration | Tracking cutoff | Available edges | Shared edges |
|---|---:|---:|---:|
| C. elegans protein, embryo 1 | 255 | 299 | 188 |
| C. elegans RNA, embryo 1 | 255 | 231 | 188 |
| C. briggsae RNA, AF16 210519ZZY0874p1 | 148 | 210 | 188 |
| C. briggsae RNA, AF16 210519ZZY0874p5 | 156 | 213 | 188 |

The protein representation uses the established 20 selected z-scored protein
features. Both RNA matrices use the same frozen 20 selected TFs and cosine
distance on stored expression values, with no additional transform. Selected
vectors are checked for finite values and nonzero norms. Features and input
SHA-256 hashes are recorded in `provenance.json`.

The two RNA representations have the same row identities. Shared annotation
parents agree with the tree on all available edges checked. This checks data
joins, not independent cross-species biological homology or RNA measurement
provenance. The manifest retains the source annotation's `Missing` flag for
child and parent; the flag is not treated as proof that an expression vector
was measured or imputed. Predefined cutoffs are retained without claiming a
new independent validation of developmental alignment.

## Geometry and endpoints

The existing 3D convention uses the current tracking loaders. C. briggsae z
and xy units have not been established as commensurate: the metadata supplies
`0.44um/94slices` but no xy pixel size or definitive exported-coordinate
conversion. These panels explicitly say that the axial scale is unverified.
No anatomy-based rescaling is introduced. XY-only calculations drop the first
coordinate in both species and are displayed on a separate review page.

First-cousin null means, variances, and covariance are exact. There are 1,000
display draws per configuration, seed 42. Each weighted assignment problem
scales the two cost matrices by their respective null standard deviations.
The pilot uses 300 intervals (301 weights including endpoints), with a nested
1,200-interval check on all eight matched configuration/geometry combinations.

RNA endpoint assignments have ties: equally state-optimal assignments can
have different travel costs. The pilot resolves each endpoint by minimizing
the other objective, using diminishing perturbations and verifying that the
primary optimum is preserved to floating-point summation tolerance. It also
requires agreement of the secondary cost at two successive perturbations.
Primary gaps, secondary savings, perturbation size, and tolerance are cached.
This changes endpoint spans relative to arbitrary solver-selected endpoints;
it explains why these canonical distances differ from the earlier read-only
proposal. Interior ties are not exhaustively enumerated, and maximum retention
means the maximum observed among sampled solutions.

Each front is separately endpoint-normalized, travel optimum `(0,1)` to
state optimum `(1,0)`. `d_LP` selects an actual sampled assignment. `d_NP`
measures the distance from the exact cousin-null mean to that same assignment.
Canonical `u` is normalized arc length, with travel at 0 and state at 1.
Connecting lines are drawing aids, not extra feasible assignments. The zoom
uses its own axis limits; compare distances numerically, not by screen length.

## Pilot findings

| Matched configuration | d_LP, existing 3D | d_LP, XY | Max retention, existing 3D | Max retention, XY |
|---|---:|---:|---:|---:|
| C. elegans protein | 0.0083 | 0.0360 | 0.8989 | 0.6809 |
| C. elegans RNA | 0.0155 | 0.0848 | 0.8191 | 0.5691 |
| C. briggsae AF16 p1 RNA | 0.0463 | 0.0574 | 0.6596 | 0.6330 |
| C. briggsae AF16 p5 RNA | 0.0493 | 0.0631 | 0.6277 | 0.5638 |

All four natural lineages are closer to their selected front assignment than
their corresponding first-cousin null mean is. The apparent retention
separation in the existing 3D calculation is substantially reduced by XY
projection. The C. elegans RNA versus C. briggsae proximity ranking also changes.
Consequently these results do not yet establish a species ranking independent
of spatial geometry.

Matching is selective: only 1/32 programmed-death cells, 7/20 glial cells,
and 52/97 neurons in the original protein set remain, compared with 53/61
muscle cells. Protein maximum retention rises from 0.8194 on all 299 edges
to 0.8961 on the CE-only matched 231 and 0.8989 on the common 188. Much of
this change therefore occurs when restricting to RNA-covered CE cells.

Across the eight matched runs, the dense sweep leaves `d_LP`, `d_NP`, and
maximum retention unchanged to numerical precision; the largest absolute
change in `u_L` is below 2e-7. This is grid stability, not exhaustive recovery
of all nondominated discrete assignments or tied optima.

## Review clarifications at the pause

### Fourth review page, panel A: coverage, not optimized retention

The numbers in `coverage_and_set_sensitivity.pdf` panel A are **matched /
available terminal cells**, grouped by cell type. The numerator is membership
in the common 188-edge set across all four configurations; the denominator is
the corresponding type count in the original 299-edge CE protein set. Thus
neurons are 52/97, programmed-death cells 1/32, glia 7/20, and muscle 53/61.
These counts are determined before optimization. They are neither maxima of
one front nor maxima pooled across four fronts. Gray bars show availability;
blue bars show common-set membership. Panel B, separately, reports maximum
observed edge retention for each configuration and cohort.

On resumption, make the panel A wording explicitly say "Matched / available
terminal cells". This label improvement is deferred; the paused plots have
not been rerendered.

### Why including depth can increase CE edge retention

A read-only diagnostic on the cached matched CE travel matrices finds the
natural biological parent strictly closer than every other biological parent
for 127/188 children in existing 3D, versus 60/188 in XY. Seventy children
change from not strictly nearest in XY to strictly nearest in 3D; three change
in the reverse direction. Duplicate slots for the same biological parent are
excluded as alternatives, and the strict distance margin is greater than
1e-10 in the stored distance units. These counts were rechecked at the pause.

This local nearest-parent diagnostic is distinct from the globally
capacity-constrained assignment: the travel-optimal assignment retains
127/188 edges (67.6%) in 3D versus 64/188 (34.0%) in XY. The two CE modalities
share the same travel matrices. The agreement of the two 3D counts is not an
equivalence of algorithms. Their maximum retention over the sampled mixed
objectives is different again, as reported in the findings table above.

Depth disambiguation is a plausible geometric explanation: XY projection can
make unrelated parents appear close when they are separated in z. The
diagnostic supports that explanation in this CE geometry, but does not prove
that development optimizes travel or that adding a reliable axis must always
increase retention. It also does not validate CB coordinate calibration or
establish a species effect. A global distance rescaling cannot repair an
incorrect relative z-to-xy scale.

The nearest-parent counts are an exploratory diagnostic, not a separately
versioned pipeline output. To reproduce them without regenerating any cache,
run this Python snippet in the `dev` environment from the repository root:

```python
import numpy as np

with np.load("terminal_pareto/output/pilot_cross_dataset/assignments_and_costs.npz") as cache:
    nearest = {}
    for geometry in ("raw3d", "xy"):
        prefix = f"matched__{geometry}__ce_protein__"
        costs = cache[prefix + "xm"]  # parent slots by candidate children
        parents = cache[prefix + "parents"]
        other = parents[:, None] != parents[None, :]
        margin = np.where(other, costs, np.inf).min(axis=0) - np.diag(costs)
        nearest[geometry] = margin > 1e-10
        print(geometry, int(nearest[geometry].sum()), len(parents))
    print("XY to 3D", int((~nearest["xy"] & nearest["raw3d"]).sum()))
    print("3D to XY", int((~nearest["raw3d"] & nearest["xy"]).sum()))
```

## Files and reproducibility

- `cell_manifest.csv`: available edges, common membership, type, major region,
  and source annotation flags.
- `input_audit.csv`, `coverage_by_cell_type.csv`, `coverage_by_region.csv`:
  coverage, source cutoffs, and coordinate extents.
- `fronts.csv`, `metrics.csv`, `null_clouds.csv`: raw and normalized costs,
  exact null moments, landmarks, and display samples. Cohorts are `matched`,
  `matched_dense`, `available`, and `ce231`.
- `assignments_and_costs.npz`: ordered terminals/parents, cost matrices, and
  permutations. For a permutation `p`, parent slot `i` receives child `p[i]`;
  child `j` receives biological parent `parents[argsort(p)[j]]`. Dense runs
  reuse their matched base cost matrices.
- `convergence.csv`: dense-minus-base metric changes.
- `provenance.json`: settings, features, source hashes, and unresolved limits.
- `pilot_review.pdf`: endpoint comparison, XY comparison, null-SD companion,
  and coverage/set sensitivity. Individual PDF and PNG pages are also saved.

From the repository root:

```bash
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/pilot_cross_dataset.py
env MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache conda run -n dev \
  python terminal_pareto/pilot_cross_dataset.py --render-only
conda run -n dev tectonic terminal_pareto/pilot_cross_dataset.tex \
  --outdir terminal_pareto/output/pilot_cross_dataset
```

The renderer-only mode reads numerical caches; rerun the analysis after
changing data or numerical code. Generated artifacts remain ignored by Git.
The commit preserves the script, TeX wrapper, findings and handoff, not the
generated PDFs/PNGs, CSVs or NPZ. Those remain locally under the output path
above. A fresh checkout therefore needs the ignored source data and feature
selection files as well as regeneration; Git alone is not an artifact backup.

## Validation checkpoint

Implementation-stage validation passed: Python compilation, the canonical
self-test (17/17), toy checks of endpoint tie handling and biological-parent
retention, exact cousin moments against a small exhaustive enumeration, and
replay of all 15,628 cached permutations against both cost matrices and
biological-parent retention. Dense-sweep stability is summarized above.
The standalone plots and captioned wrapper were rendered, inspected, and
compiled successfully during implementation.

At the 2026-09-14 documentation pause, all 16 source hashes recorded in the
local `provenance.json` still matched, and the nearest-parent diagnostic was
rechecked. No new numerical sweep or figure rendering was needed for the
documentation-only pause. Prior validation is not a substitute for validation
after future input or numerical-code changes.

## How to resume

1. Read this record and `PUBLICATION_FIGURES_HANDOFF.md`, then inspect the
   four-page review and captioned comparison. Keep this pilot separate from
   accepted figures and from full-tree work; figure numbering is undecided.
2. Inspect the worktree and compare current inputs/code with the hashes in
   `provenance.json`. Preserve the current output directory before any rerun:
   the script defaults to the same path (use `--out` for a separate run).
   If caches are absent, restore
   the documented inputs before regenerating. Do not mix old caches with new
   numerical code. Retain all four configurations and the 188-edge baseline.
3. Resolve CB exported coordinate order/units and relative z-to-xy calibration
   from acquisition/export provenance, not by fitting anatomy or optimizing
   retention. Separately audit developmental alignment at the fixed cutoffs
   and measured/aggregated/imputed RNA-state provenance. Record unresolved
   items explicitly if they cannot be resolved.
4. Decide the scientific comparison with the authors. Preserve raw-3D and XY
   as distinct sensitivities; any calibrated geometry should be a new,
   separately labeled pilot result. Keep molecular panels and common parent
   capacities fixed when isolating geometry. Retain the available-set and
   CE-only 231-edge controls to distinguish coverage from geometry effects.
5. Make the deferred coverage-label clarification. If only labels change,
   use `--render-only` after preserving the checkpoint. For numerical changes,
   regenerate null moments, endpoints, sweeps and convergence checks together;
   replay assignments, inspect plots, and compile the wrapper. Promote the
   nearest-parent diagnostic to a tested output only if used in a figure or
   substantive quantitative claim.
6. Review the story, captions and limitations before proposing publication
   integration. Do not interpret tracking-replicate spread as RNA uncertainty,
   modality differences as a controlled modality effect, or sampled retention
   maxima as exhaustive maxima over all tied/nondominated assignments.

Before a publication release, resolve spatial calibration, stage alignment,
and measured/aggregated/imputed RNA-state provenance; decide which geometry
supports the intended claim. Protein/RNA feature panels differ, and tracking
replicates share molecular data, so neither a pure modality effect nor RNA
replicate uncertainty is estimated here.
