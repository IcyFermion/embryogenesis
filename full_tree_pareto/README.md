# Full-tree Pareto analysis

Updated 2026-09-29. The current rebuild produces **one main Figure 7** and
**one supplementary figure with a heuristic inventory**, using pooled travel
from all three *C. elegans* tracking embryos and endpoint-normalized displays.
See [POOLED_HANDOFF.md](POOLED_HANDOFF.md) for current execution/validation
status and resumption instructions. The former two-main-figure configuration
is preserved in [LEGACY_EMBRYO1.md](LEGACY_EMBRYO1.md).

## Build and outputs

Run from the repository root in the `dev` Conda environment:

```bash
# Build/resume scientific caches, render and compile both figures.
# On a fresh run, checkpoint the expensive paired sweep first:
python -m full_tree_pareto.resume_paired
python -m full_tree_pareto.publication_build

# Optional override for another machine (both commands default to 24 workers).
python -m full_tree_pareto.resume_paired --workers 8
python -m full_tree_pareto.publication_build --workers 8

# Replay-validate assignments and redraw without rerunning solvers.
python -m full_tree_pareto.publication_build --layout-only

# Explicit promotion; preserve the former publication in a hashed archive.
python -m full_tree_pareto.publication_build --layout-only --publish

python -m full_tree_pareto.publication_build --verify
python -m unittest discover -s full_tree_pareto -p 'test_*.py'
```

Working run: `output/runs/pooled_full_tree_v1/terminal_clamped_20260927/`.
`analysis/` holds source/cohort identity, matrices, node manifest, complete
parent-array sweeps, raw costs, null draws, Gaussian covariances and validation.
`analysis/cousin_shuffles/` adds independently hash-pinned second- and
third-cousin draws and leaf permutations, without changing the solver caches.
`publication/` holds panels, two compiled wrappers and display transforms.
Each redraw first saves a hash-checked preceding layout under `layout_history/`.
Wrapper stems are `fig7_ce_full_tree_layerwise` and
`figs_ce_full_tree_heuristics`; the latter is provisionally Figure S4, following
terminal supplements S1--S3.
Wrapper templates/captions are versioned in `publication_build.py`.

Promotion copies figures to `output/publication/` with a release manifest,
preserving the former bundle in `output/legacy/releases/`. A failed promotion
removes its staging folder and partial archive and restores the former bundle.
The manuscript is not automatically edited. Everything under `output/`,
including compiled wrappers and local archives, is Git-ignored and untracked
(the two historical wrappers were untracked on 2026-09-29), so it requires
external backup. `--run PATH` selects an isolated run;
`--workers N` controls parallel solver processes (default **24** in both
entrypoints; positive integers only). Use `--workers 1` for one solver process.
This is a process count, not CPU affinity: the OS schedules work across cores,
and 24 workers need not yield a fourfold speedup over six. Lower it on machines
with fewer cores or less memory. Keep `OPENBLAS_NUM_THREADS=1` to avoid nested
BLAS oversubscription. Completed methods resume
after identity/hash checks. Identity includes ordered cohort, source data,
algorithm code, sweep, draws and seed. Changed configurations require new
runs: never rewrite manifests to force reuse. Legacy scalar-cost caches are
not inputs. Historical renderer commands target `output/legacy/rebuild/`;
that folder does not exist until `fig7_fig8_ce_full_tree_pareto.py` is rerun,
so the historical `figs_`/`table_ce_full_tree_heuristics.tex` wrappers compile
only after such a rebuild.

Sweep weights, draws and seed are the `INTERVALS`, `DRAWS` and `SEED`
constants in `publication_build.py`; builds, `resume_paired`, captions and
display subsets all read them. `DEFAULT_WORKERS = 24` in the same module is
the shared runtime default for the build and checkpoint entrypoints; their
`--workers` options override it. Worker count is not part of scientific cache
identity. The frozen low-level `pooled_analysis.build` fallback remains four,
but both supported entrypoints pass the chosen count explicitly.
**Do not edit `pooled_analysis.py` casually:**
its source hash is part of every cache identity, so any change, even a
comment, invalidates all caches, including the expensive paired sweep.
`resume_paired` saves each finished weight separately under `checkpoints/`,
verifies those checkpoints on restart, and then assembles the standard paired
cache. Use the same `--run` path for it and the publication builder.

## Measurement and cohort

Tracking uses cutoffs **255, 247, 225**, matching the terminal pipeline.
Cells appearing only once exactly at cutoff are excluded. Positions use the
last retained row and scale 0.1625. The strict tracking/protein intersection
is an ancestor-complete binary forest: **978 measured cells, 487 internal,
491 leaves, four roots and 974 scored edges**. Roots are ABa, ABp, EMS and P2.
Unmeasured P0/AB/P1 anchors are not scored. Every edge spans one canonical
generation; there are no skipped ancestors or imputed positions.

The node manifest records canonical biological-terminal status, measured-tree
leaf status, last tracking times and cutoff flags in all embryos. A measured
leaf may result from biological termination, cutoff or cross-embryo coverage.
All are included in the terminal clamping below.

For embryo r, T_r(L) is its natural travel total on these **974 full-tree
edges**, not the terminal pipeline's 275-edge denominator:

\[
\bar d_{ij}=\frac{1}{3}\sum_{r=1}^{3}\frac{\lVert x_i^{(r)}-x_j^{(r)}\rVert}{T_r(L)}.
\]

The same denominators apply to every round, heuristic and reference. Pooled
natural travel is one. Pairwise distances are pooled; coordinates are never
averaged into a synthetic embryo. Travel is a retained-position displacement
proxy, not an integrated trajectory or metabolic cost.

Molecular distance remains **Euclidean distance across the selected 20
z-scored protein features**. Terminal figures use cosine dissimilarity;
harmonizing travel and displays does not silently change this distinction.

## Optimization and presentation

The eight contraction rounds contain **491, 223, 123, 66, 38, 18, 10, 5** edges.
They partition the forest and are asynchronous bottom-up rounds, not uniform
developmental stages. Each round assigns children to the multiset of natural
parent slots by exact Hungarian assignment. Retention compares biological
parents, so swaps of duplicate slots belonging to one parent remain retained.

Every method uses **301 endpoint-inclusive weights**. During optimization,
both objectives are divided by their **global exact first-cousin-null SDs**,
matching the terminal pipeline convention. These shared scales avoid arbitrary
unit imbalance after pooling and allow round costs at a common weight to sum
to an attained optimum in the product of round-wise assignment spaces. This
is not unrestricted global full-tree optimization. Optimizer scaling is
separate from endpoint display normalization.

- **Figure 7A:** eight round-specific fronts colored by parent retention,
  with natural lineage and maximum-retention markers.
- **Figure 7B:** aggregate layerwise front only, first-cousin shuffle,
  second- and third-cousin shuffles, and random rebuild in an inset; no
  spanning-forest front or Gaussian reference. The inset uses the same
  endpoint transform with separate axis limits, matching the terminal plot.
- **Supplement:** five regenerated heuristic curves and all seven reference
  families, plus a constraint/audit inventory.

Each curve shows its own sampled nondominated cost pairs, not a cross-method
global optimum or the complete discrete Pareto set. Bottom-up-by-layer and
paired bottom-up permit root identities to move; layerwise, top-down and
degree-constrained forest fix four observed roots. Every method saves full
parent arrays and validates cohort coverage, 974 edges, binary degree, four
components and absence of cycles. Different feasible sets still preclude
equally constrained rankings. The historical terminal-only algorithm remains
excluded due to unresolved constructed-state/scope/scoring issues.

For travel optimum A and state optimum B:

\[
D_1(a)=\frac{T(a)-T(A)}{T(B)-T(A)},\qquad
D_2(a)=\frac{S(a)-S(B)}{S(A)-S(B)}.
\]

A maps to (0,1), B to (1,0). Each round has its own anchors. Main collective
and supplementary comparisons **share aggregate layerwise anchors** for all
curves, observations and nulls. Values outside the endpoint box are not clipped.
The implementation directly reuses terminal `EndpointTransform`, records its
parameters in JSON, and rejects degenerate spans explicitly.

## Cousin-shuffle references

The main-panel reference set was revised at the author's request on
2026-09-29: remove Gaussian from Figure 7B and add second- and third-cousin
shuffles. The phylogenetic analogy is not central to the story, and the
Gaussian cloud is far from the observed lineage and reconstruction fronts;
its distinct state space also limits direct comparison. Retain it in the
supplement and scientific cache rather than discard the sensitivity record.

All three cousin models permute **measured leaf identities only** within
canonical ancestor groups two, three or four transitions back (first, second
and third cousins). This includes cutoff/coverage leaves, not just the 275
canonical biological terminals. Position across all embryos and protein state
move together. Internal identities, root identities, binary parent-slot counts
and internal edges stay fixed. Every draw scores all 974 edges; the unchanged
internal contribution is included, not omitted. Sibling exchanges may leave
the biological parent unchanged. These descriptive randomizations do not
represent all developmental constraints or provide calibrated significance.

The higher-cousin implementation reuses the terminal pipeline's canonical
ancestry/grouping functions, including unscored ancestors above the four
measured roots. Unshufflable singleton groups remain fixed. Each family has
10,000 draws (seeds 43 and 44); 1,000 are displayed and means use all draws.
Saved leaf permutations are checked for group membership, unique coverage and
score replay on every load. A separate source/cohort/seed/hash identity guards
the add-on cache, so `pooled_analysis.py` and existing reconstruction caches
remain untouched. Normal builds generate missing add-on caches; layout-only
builds require them already present. Promotion includes nested analysis hashes.

## Terminal-clamped Gaussian reference (supplement only)

Author-confirmed change, 2026-09-27: **freely simulate internal states from
observed roots, then clamp all measured leaves to observed values**. This is
not the earlier jointly leaf-conditioned Gaussian sensitivity.

For each embryo separately, fit zero-drift spatial covariance as the mean
outer product of observed displacements divided by that embryo's edge
duration. Fit one protein covariance as the mean outer product of observed
protein increments per canonical transition. All 974 edges enter these fits.
Keep full covariance within blocks, independent spatial/protein blocks and
independent spatial processes across embryos; no cross-embryo covariance
model is estimated from three replicates.

Internal descendants inherit simulated parent states plus Gaussian increments.
Roots stay observed. Every measured leaf's coordinates in each embryo and
its protein vector are copied exactly from observations. Terminal costs are
recomputed from the fixed leaf to its freely simulated parent. Thus terminal
increments are **not ordinary Brownian increments after clamping**, and
internal states do not respond to endpoint information as in Gaussian
conditioning. Omitting unused leaf innovations is equivalent to generating
them and then overwriting the leaves.

Each draw pools its three spatial totals using the observed denominators and
scores the shared protein realization once. Defaults: seed 242, 10,000 draws,
deterministic 1,000-draw display subset; plus markers use all draws. Tests
assert exact root/leaf values, free internal moments, determinism and pooled
score recomposition.

This holds observed outcomes fixed but **does not equalize feasible sets**:
Gaussian internal states are arbitrary, while heuristics use measured internal
states. Long terminal reconnections may dominate costs. Packing, mechanics,
intermediate developmental constraints and protein-rate heterogeneity remain
absent. Neither separation from natural lineage nor proximity to a front is
a calibrated significance test or evidence for selection on efficiency.
Independent pooling changes reference dispersion; it does not estimate
population variation from three embryos.

Historical root-fixed leaf-free, branch-variable and jointly leaf-conditioned
analyses remain in [LEGACY_EMBRYO1.md](LEGACY_EMBRYO1.md),
[BRANCH_VARIABILITY_SENSITIVITY.md](BRANCH_VARIABILITY_SENSITIVITY.md), and
[LEAF_CONDITIONED_REFERENCE.md](LEAF_CONDITIONED_REFERENCE.md). Their numerical
results apply to the old cohort and are not mixed with this build.

## Source map

| Source | Role |
|---|---|
| `pooled_analysis.py` | Cohort, matrices, five heuristics, references and replay-validated caches |
| `cousin_references.py` | Separate second-/third-cousin caches using terminal ancestry definitions and full-edge replay |
| `publication_build.py` | Current panels, captions, compilation and explicit archived promotion |
| `resume_paired.py` | Per-weight checkpointed paired sweep, assembled into the standard cache |
| `test_pooled_analysis.py` | Pooled scoring, solver invariants, endpoints, clamping and cache guards |
| `test_publication_build.py` | Endpoint display, release hashes, caption settings and promotion rollback |
| `test_cousin_references.py` | Canonical groups, terminal sampler agreement, forest invariants and cache/replay guards |
| `POOLED_HANDOFF.md` | Execution checkpoint and remaining work |
| `publication_analysis.py`, `fig7_fig8_ce_full_tree_pareto.py` | Historical single-embryo pipeline |
| `clock_sensitivity.py`, `branch_variability_sensitivity.py`, `leaf_conditioned_reference.py` | Historical reference sensitivities |
