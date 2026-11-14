# Full-tree Pareto analysis

Pooled travel from all three *C. elegans* tracking embryos over the measured
full tree, analysed layer by layer, plus the partial-forest cross-species
comparison. This package is numerical only: figures and captions are drawn by
[`publication/`](../publication/README.md). The manuscript set uses this
package for **Figure 6** (full tree by layer: 6A redraws every layer's front
from the pooled run's per-round references, 6B-C come from the Track E runs)
and for the backlog page 7C (*C. briggsae* partial forest by layer). The
released Figures 7, S4, S5, 10 and 11 are backlog figures
([`publication/BACKLOG.md`](../publication/BACKLOG.md)).

## Build and outputs

Run from the repository root in the `dev` Conda environment:

```bash
# Build/resume scientific caches. On a fresh run, checkpoint the expensive
# paired sweep first (both commands default to 24 workers):
python -m full_tree_pareto.resume_paired [--workers 8]
python -m full_tree_pareto.pooled_pipeline [--workers 8]

# Replay-validate existing caches without solving.
python -m full_tree_pareto.pooled_pipeline --layout-only
python -m unittest discover -s full_tree_pareto -p 'test_*.py'   # about 26 s

# Figures: see publication/README.md.
```

Working run: `output/runs/pooled_full_tree_v1/terminal_clamped_cosine_20261014/`
(cosine; the Euclidean `terminal_clamped_20260927/` is superseded).
`analysis/` holds source/cohort identity, matrices, node manifest, complete
parent-array sweeps, raw costs, null draws, Gaussian covariances and validation.
`analysis/cousin_shuffles/` adds independently hash-pinned second- and
third-cousin draws and leaf permutations, and `analysis/round_references/` the
per-round references (below), without changing the solver caches.
`checkpoints/` holds the per-weight paired sweep. Everything under `output/` is
Git-ignored and requires external backup. Caches are always generated fresh,
never rescaled from embryo-1 results (author decision).

`--run PATH` selects an isolated run; `--workers N` controls parallel solver
processes (default **24** in both entrypoints; positive integers only). Use
`--workers 1` for one solver process. This is a process count, not CPU
affinity: the OS schedules work across cores, and 24 workers need not yield a
fourfold speedup over six. Lower it on machines with fewer cores or less
memory. Keep `OPENBLAS_NUM_THREADS=1` to avoid nested BLAS oversubscription.
Run long sweeps detached (for example as a systemd user service such as
`full-tree-pareto-20260928.service`): interactive sessions died mid-sweep.
Completed methods resume after identity/hash checks. Identity includes ordered
cohort, source data, algorithm code, sweep, draws and seed. Changed
configurations require new runs: never rewrite manifests to force reuse.
Legacy scalar-cost caches are not inputs. Historical renderer commands target
`output/legacy/rebuild/`; that folder does not exist until
`fig7_fig8_ce_full_tree_pareto.py` is rerun.

Sweep weights, draws and seed are the `INTERVALS`, `DRAWS` and `SEED`
constants in `pooled_pipeline.py`; builds, `resume_paired` and the publication
captions all read them. `DEFAULT_WORKERS = 24` in the same module is
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
generation; there are no skipped ancestors or imputed positions. The 491
leaves include exactly the terminal pipeline's 275 cells; cutoff-flagged
leaves number 466, 464 and 466 in embryos 1-3.

The node manifest records canonical biological-terminal status, measured-tree
leaf status, last tracking times and cutoff flags in all embryos. A measured
leaf may result from biological termination, cutoff or cross-embryo coverage.
All are included in the terminal clamping below.

For embryo r, T_r(L) is its natural travel total on these **974 full-tree
edges** (4417.1253, 4407.4359 and 4403.2182), not the terminal pipeline's
275-edge denominator:

\[
\bar d_{ij}=\frac{1}{3}\sum_{r=1}^{3}\frac{\lVert x_i^{(r)}-x_j^{(r)}\rVert}{T_r(L)}.
\]

The same denominators apply to every round, heuristic and reference. Pooled
natural travel is one. Pairwise distances are pooled; coordinates are never
averaged into a synthetic embryo. Travel is a retained-position displacement
proxy, not an integrated trajectory or metabolic cost.

Molecular distance is **cosine distance across the selected 20 z-scored
protein features**, as in the terminal analyses (author decision, 2026-10-14,
Track A2; see `terminal_pareto/METHODS_cell_state_representation.md`). Runs
before that date used Euclidean distance and no longer validate against the
current code. The Gaussian reference scores spatial edges by Euclidean length
and protein edges by cosine distance. No full-tree cohort has a zero
expression vector.

Key numbers of the cosine run: natural cell-state total 389.0792 (first-cousin
state SD 1.39260; travel SD 0.00830310, unchanged from Euclidean); the
layerwise sweep has 259 distinct cost pairs and maximum retention 0.5298.
Reference means (travel/state): first cousin 1.10654/398.674, second cousin
1.26112/424.221, third cousin 1.42626/456.011, internal-layer 3.88086/737.281,
full assignment 3.99238/879.734, random rebuild 3.99332/879.880, clamped
Gaussian 2.47477/536.757.

## Optimization and presentation

The eight contraction rounds contain **491, 223, 123, 66, 38, 18, 10, 5** edges.
They partition the forest and are asynchronous bottom-up rounds, not uniform
developmental stages. Text and figures number them top-down as **layers**
(layer = 9 − round; layer 8 holds the leaves); code keeps rounds. Each round
assigns children to the multiset of natural parent slots by exact Hungarian
assignment. Retention compares biological parents, so swaps of duplicate
slots belonging to one parent remain retained.

Every method uses **301 endpoint-inclusive weights**. During optimization,
both objectives are divided by their **global exact first-cousin-null SDs**,
matching the terminal pipeline convention (kept for the per-round analysis,
author decision 2026-10-05). These shared scales avoid arbitrary unit
imbalance after pooling and allow round costs at a common weight to sum to an
attained optimum in the product of round-wise assignment spaces. This is not
unrestricted global full-tree optimization. Optimizer scaling is separate from
endpoint display normalization.

Released figures (backlog): Figure 7A, eight round-specific fronts colored by
parent retention (redrawn as manuscript Figure 6A); Figure 7B, the aggregate
layerwise front with first-, second- and third-cousin shuffles and a
random-rebuild inset; provisional S4, five regenerated heuristic curves and all
seven reference families, plus a constraint/audit inventory.

Each curve shows its own sampled nondominated cost pairs, not a cross-method
global optimum or the complete discrete Pareto set. Bottom-up-by-layer and
paired bottom-up permit root identities to move; layerwise, top-down and
degree-constrained forest fix four observed roots. Every method saves full
parent arrays and validates cohort coverage, 974 edges, binary degree, four
components and absence of cycles. Different feasible sets still preclude
equally constrained rankings. The historical terminal-only algorithm remains
excluded due to unresolved constructed-state/scope/scoring issues (998 edges,
roots not fixed).

For travel optimum A and state optimum B:

\[
D_1(a)=\frac{T(a)-T(A)}{T(B)-T(A)},\qquad
D_2(a)=\frac{S(a)-S(B)}{S(A)-S(B)}.
\]

A maps to (0,1), B to (1,0). Each round has its own anchors. Aggregate
comparisons **share aggregate layerwise anchors** for all curves, observations
and nulls. Values outside the endpoint box are not clipped. The implementation
directly reuses terminal `EndpointTransform`, records its parameters in JSON,
and rejects degenerate spans explicitly.

## Per-layer fronts (Track E; manuscript Figure 6)

`round_references.py` (`--target pooled|cross-species|all`) writes
`analysis/round_references/` in both cosine runs: per-round first-, second-
and third-cousin and full-random draws (2,000 each), exact first-cousin
moments and canonical metrics per round, on the cached sweeps (fronts are never
re-solved). A per-round cousin null is used for every round (author,
2026-10-05). Checks: round-1 moments equal the pooled leaf null; the
partial-forest round-1 null plus the fixed internal edges equals the cached
aggregate null; fronts equal `fronts.csv`.

Track E (closed 2026-11-05) asked first whether layer-specific features beat
the global top-20 (`layer_features.py`, `layer_ranking.py`, `layer_fronts.py`;
runs `layer_features_v1/e0_groupweight_20261015`, `e1_20261015`,
`e1_fronts_20261015`; `e0_20261015` is superseded): they do not, and the
global set stands. It then measured, layer by layer, how much the natural
lineage beats a cousin shuffle on travel and on cell state
(`layer_hourglass.py`, `layer_stages.py`; runs `hourglass_20261019` and
`stages_20261020`, which manuscript Figure 6B-C read). Methods, numbers and the
agreed claims are in [`LAYER_FEATURES.md`](LAYER_FEATURES.md) and
[`LAYERWISE_HOURGLASS.md`](LAYERWISE_HOURGLASS.md): the travel advantage is
smallest in the mid-embryo layers, in both species, while the cell-state
advantage stays flat; claims are limited to "consistent with" fate-based cell
sorting and "coincides with" the expression hourglass waist. Layers 1-2 (5 and
10 edges) are not drawn; claims start at layer 3.

Caveats: in pooled rounds 6-7 (18 and 10 edges) the natural lineage sits at
the travel optimum (u 0, d_LP 0, retention 1); partial-forest round 6 (2
edges) is degenerate for CE RNA and CB RNA in both geometries.

## Cousin-shuffle references

The reference set was revised at the author's request on 2026-09-29: remove
Gaussian from the main aggregate panel and add second- and third-cousin
shuffles. The phylogenetic analogy is not central to the story, and the
Gaussian cloud is far from the observed lineage and reconstruction fronts;
its distinct state space also limits direct comparison. It is retained in the
supplement and scientific cache rather than discarded.

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
measured roots. Unshufflable singleton groups remain fixed; second and third
cousins have 97 and 55 groups with more than one cell, covering 489 leaves.
Each family has 10,000 draws (seeds 43 and 44); 1,000 are displayed and means
use all draws. Saved leaf permutations are checked for group membership,
unique coverage and score replay on every load (`--layout-only` replays 1,505
forests and the 20,000 cousin draws). A separate source/cohort/seed/hash
identity guards the add-on cache, so `pooled_analysis.py` and existing
reconstruction caches remain untouched: editing `cousin_references.py`
invalidates only `analysis/cousin_shuffles/`. Normal builds generate missing
add-on caches; layout-only builds require them already present.

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
analyses are in [BRANCH_VARIABILITY_SENSITIVITY.md](BRANCH_VARIABILITY_SENSITIVITY.md)
and [LEAF_CONDITIONED_REFERENCE.md](LEAF_CONDITIONED_REFERENCE.md) (embryo 1,
see below). Their numerical results apply to the old cohort and are not mixed
with this build.

## Cross-species partial-forest comparison

Run: `output/runs/cross_species_layerwise_v1/terminal_anchored_cosine_20261014/`
(cosine, analysis id `d2c0850b…`; the Euclidean `terminal_anchored_20260930/`
and its id `4f5f7e9f…` are retired). It compares CE protein, CE RNA and CB AF16
RNA in 3D and XY, using only layerwise assignment.

Why a partial forest: the species share 803 nodes and 762 direct edges, but
with 41 disconnected roots, 7 upper islands and 64 unary parents, and both RNA
tables lack the four early roots and 29 internal nodes. The author therefore
approved a bottom-up partial-tree comparison rather than imputing states or
skipping ancestors. Starting from the same **187 matched biological terminal
edges** as the terminal comparison, each branch ascends until its first
unavailable canonical ancestor and never resumes above a gap (6, 17, 62, 86
and 16 terminals reach 2-6 transitions). This retains **485 measured cells,
454 scored edges, 31 boundary roots and 187 leaves**. Contraction rounds
contain **187, 119, 80, 48, 18, 2** edges. The 298 retained internal nodes
include 142 with one represented child and 156 with two; the observed
parent-slot capacities are preserved exactly. These are coverage boundaries,
not a claim of biologically unary divisions. This is not a complete or
near-complete embryonic tree.

Travel pools three CE embryos and two AF16 embryos (`210519ZZY0874p1` at cutoff
148, `p5` at 156; CB coordinates reordered from z,x,y to x,y,z; XY drops z)
using fixed natural totals on this 454-edge reference, separately for 3D and
XY. It does not use the terminal 275-edge denominators or reuse the full
tree's 974-edge caches. Molecular distance is **cosine**: frozen top-20
z-scored protein or shared 20 RNA TFs on stored values. The six aggregate
fronts use global exact first-cousin-null SD scaling, 301 weights, checked
endpoint tie handling and nested 1,201-weight checks (which change no maximum
retention; the largest shift is CE protein 3D, Δu −0.00145, Δd_CP −0.00043).
Each configuration has its own endpoint display spans. Cousin shuffles move
only the 187 terminals; random rebuild fixes the 31 roots and the 142/156
one-/two-child capacities and is descriptive, not uniform. Each reference has
10,000 draws (seeds 42-45), paired across the six cost spaces. Validation
replays 9,012 assignments and 240,000 reference realizations and checks 24
source hashes. Cite `analysis/metrics.csv` (for example 3D CE protein: u
0.292, d_LP 0.0828, d_CP 0.1287, retention 0.6476).

Limits: it answers a question about the lower lineage and is not a modality
test; replicates do not estimate molecular uncertainty; stage alignment and
RNA provenance are unresolved; there is no calibrated significance.

```bash
OPENBLAS_NUM_THREADS=1 python -m full_tree_pareto.cross_species_analysis
python -m full_tree_pareto.cross_species_analysis --verify-only
python -m unittest full_tree_pareto.test_cross_species
```

All commands run from the repository root in `dev`. They accept `--run-id`
for another isolated run. Completed case caches resume after identity/hash
checks. Changed scientific inputs or source require a new run; do not edit
manifests to force reuse.

## Legacy embryo-1 analysis

The first full-tree analysis used embryo 1 only, Euclidean distance and 1,000
edges: Figure 7 (maximum retention 49.9%) and Figure 8, which compared
layerwise assignment with a degree-constrained spanning forest against cousin,
separate-clock Gaussian and random-rebuild references. The pooled 974-edge
cohort replaced it on 2026-09-27; Figure 8 was retired on 2026-10-01, and the
single-embryo 7B, heuristic inventory and methods PDFs on 2026-10-04. Lessons
that still apply: retention counts biological parents, not slots; the
terminal-only rebuild stays withheld; `mst_rebuild.npz` and `phylo_bm.npz` are
invalid. Its caches stay loose in `output/` (`ce_full_tree_*.csv`,
`internal_opt/`, `clock_sensitivity/`, `branch_variability/`,
`leaf_conditioned/`); the code is `publication_analysis.py`,
`fig7_fig8_ce_full_tree_pareto.py` and the sensitivity modules below.

## Source map

Numerical code only; drawing and captions are in
[`publication/`](../publication/README.md).

| Source | Role |
|---|---|
| `pooled_analysis.py` | Cohort, matrices, five heuristics, references and replay-validated caches |
| `pooled_pipeline.py`, `test_pooled_pipeline.py` | Sweep settings and compute-only cache build/replay |
| `resume_paired.py` | Per-weight checkpointed paired sweep, assembled into the standard cache |
| `cousin_references.py`, `test_cousin_references.py` | Separate second-/third-cousin caches using terminal ancestry definitions and full-edge replay |
| `round_references.py`, `test_round_references.py` | Per-round cousin references and metrics for the pooled tree and the partial forest |
| `layer_features.py`, `layer_ranking.py`, `layer_fronts.py` and tests | Track E0-E1: layer-specific features (closed; the global top-20 stands) |
| `layer_hourglass.py`, `layer_stages.py` and tests | Track E3-E4: per-layer advantage and birth stages (manuscript Figure 6B-C) |
| `cross_species_analysis.py`, `test_cross_species.py` | Isolated terminal-anchored partial-forest cohort, six exact layerwise sweeps, references and replay guards (`fig_cross_species.py` only redirects `--render-only`) |
| `test_pooled_analysis.py` | Pooled scoring, solver invariants, endpoints, clamping and cache guards |
| `publication_analysis.py`, `fig7_fig8_ce_full_tree_pareto.py`, `heuristic_inventory.py` | Historical single-embryo pipeline |
| `clock_sensitivity.py`, `branch_variability_sensitivity.py`, `leaf_conditioned_reference.py` | Historical reference sensitivities |

**Hash-pinned modules, never edit:** `pooled_analysis.py`,
`publication_analysis.py`, `cousin_references.py`, `cross_species_analysis.py`
and `round_references.py`; `layer_features.py` and `layer_hourglass.py` pin
their runs. Import them and add new modules with their own caches, following
`cousin_references.py`.
