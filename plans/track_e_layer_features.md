# Track E: per-layer fronts in the full tree

**Revised 2026-10-19 (author).** Track E began as a question about
layer-specific features (E0-E2 below). E1 answered it: the global top-20
stands. Along the way the per-layer fronts showed a pattern worth finishing
instead: the natural lineage sits farthest from the front in the mid-embryo
layers, and the gap is spatial. The remaining work (E3-E5) turns that into a
small set of defensible claims and then stops.

**Question (revised).** Layer by layer, how much does the natural lineage beat
a cousin shuffle on travel and on cell state, and does that advantage change
with developmental stage?

**Layer convention.** Analysis code and run outputs keep contraction rounds
(`pooled_analysis.contraction_layers`), which count bottom-up: round 1 holds
the 491 measured leaves, round 8 holds 5 edges. Text and figures use layers,
counted top-down: **layer = 9 − round** (collaborator request, 2026-10-19;
matches the v2 Figure 6). Layer 1 holds the 5 cells born around the 6-cell
stage; layer 8 holds the measured leaves.

## Closed: E0-E2 (layer-specific features)

- **E0** (decided 2026-10-15): Schema D descendant-mixture labels and the
  equal-relatedness separation criterion. Run
  `layer_features_v1/e0_groupweight_20261015`.
- **E1**: per-round top-20 sets fail the held-out gate in strata 3 and 4
  (run `layer_features_v1/e1_20261015`). **The global top-20 stands.**
- **Front check** (run `layer_features_v1/e1_fronts_20261015`): per-round
  sets shift the lineage toward the cell-state end in rounds 3-5, which lineage
  memory alone would also produce. Descriptive only. The mid-layer hump in
  d_LP survives the per-round sets; this is E3's robustness check.
- **E2 is not run.** The scaling choice and the per-round fate-only control
  are dropped with it.

Details: `full_tree_pareto/LAYER_FEATURES.md`.

## Method (fixed; no new optimizer)

Every cell in the tree is measured, both position and protein. The bottom-up
tree therefore splits into independent per-layer assignment problems: a
layer's daughters are matched to its mothers' slots, with every other layer
fixed. Nothing is solved recursively.

- **Features and scaling** as in the pooled cosine run: global top-20,
  cosine cell-state distance, one global leaf first-cousin SD (Track A; E1).
- **Per layer:** the cached layerwise front and the per-round first-cousin
  shuffle (`round_references`, already run for the pooled tree and the
  cross-species partial forest).
- **Three numbers per layer:**
  - front distance d_LP (the existing canonical metric);
  - travel advantage = 1 − natural travel / first-cousin shuffle mean;
  - cell-state advantage, defined the same way on cell state.

Runs used:
- `pooled_full_tree_v1/terminal_clamped_cosine_20261014` (3 C. elegans embryos);
- `cross_species_layerwise_v1/terminal_anchored_cosine_20261014`, geometry
  raw3d: C. elegans protein and RNA, and C. briggsae (AF16, 2 embryos) RNA.

## What the existing metrics already show (2026-10-19)

Travel advantage over the first-cousin shuffle, layers 8 → 4:

| Series | L8 | L7 | L6 | L5 | L4 |
|---|---|---|---|---|---|
| Pooled C. elegans | 0.22 | 0.12 | 0.09 | 0.15 | 0.23 |
| Partial forest, C. elegans | 0.23 | 0.08 | 0.06 | 0.16 | 0.17 |
| Partial forest, C. briggsae | 0.23 | 0.08 | 0.07 | 0.18 | 0.14 |

- Against the second-cousin shuffle the pooled pattern is the same:
  0.41, 0.25, 0.22, 0.28, 0.43.
- Travel z against the first-cousin null is 3.6-12.8 in pooled layers 4-8, so
  lineage beats the shuffle in every layer.
- Cell-state advantage is small and flat (0.02-0.06).
- Layers 6-7 hold cells born at about 70-190 cells.
- **Layers 1-3** hold 5-18 edges and have degenerate front distances. They are
  plotted greyed out and carry no claims.
- **No uncertainty yet** on the advantages: this is E3.

## E3. Uncertainty for the per-layer advantages

New module `full_tree_pareto/layer_hourglass.py`; run
`layer_features_v1/hourglass_<date>`. It reuses the cached contexts and
`round_references` functions (`round_groups`, `round_blocks`) and solves
nothing.

- **Bootstrap unit: first-cousin groups within a layer.** For each group, the
  natural cost and the exact shuffle mean (each child's mean cost over its
  group's slots). Singleton groups count, with natural = shuffle mean.
  Resample groups with replacement, 2,000 draws, fixed seed; advantage =
  1 − Σ natural / Σ shuffle mean over the draw. Travel and cell state are
  resampled together.
- **Pre-specified contrast (fixed before results):**
  - mid = mean advantage of layers 6 and 7;
  - flanks = layer 8 on one side and layers 4-5 on the other, each the mean of
    its per-layer advantages;
  - dip = mid − flank, computed for each side separately.
- **Decision rule:** the dip is claimed only if the 95% interval of the dip
  excludes zero on both sides, in pooled C. elegans travel. Otherwise the
  hourglass framing is dropped and the claim shrinks to "lineage keeps cells
  close in space in every layer; its cell-state advantage is small."
- **Sensitivity**, reported but not gated:
  - the second-cousin shuffle;
  - the partial forest (C. elegans protein and RNA, C. briggsae RNA);
  - the cell-state dip, run through the same contrast.
- **Independence.** C. briggsae is an independent replicate. The C. elegans
  partial forest reuses the pooled embryos, so it is a different panel on the
  same embryos, not a replicate.
- **Reproduction gate:** with no resampling, the per-layer advantages equal
  the cached `round_references` metrics at 1e-12.

## E4. Write-up and figure

- **`full_tree_pareto/LAYERWISE_HOURGLASS.md`**, containing:
  - the method paragraph;
  - the E3 results;
  - the claims table below;
  - what is not claimed.

  `LAYER_FEATURES.md` gets a closing note pointing to it.
- **One candidate figure** in three panels:
  - front distance by layer;
  - travel and cell-state advantage by layer, both species, with intervals;
  - the birth-stage axis with the Schnabel sorting window and the Large et al.
    expression minimum marked.

  Built under `publication/output/candidates/`; production waits for the
  author.
- Update the collaborator review page to the trimmed claims.

### Claims (subject to E3)

| Claim | Wording |
|---|---|
| Lineage beats cousin shuffles on travel in layers 4-8 | Stated |
| It beats them by the least in layers 6-7 | Stated if E3 passes |
| Cell-state advantage is small and flat | Stated |
| The window matches fate-based cell sorting (Schnabel et al. 2006: 50 to ~400 cells) | "Consistent with" |
| It coincides with the expression hourglass waist (Large et al. 2025: 200 cells) | "Coincides with" |

**Not claimed:**
- the Spearman correlation with between-species JSD (7 points, both smooth in
  time);
- selection versus constraint, gene age (TAI), lethality and redundancy. These
  go in one Discussion paragraph citing the literature;
- anything about layers 1-3.

## E5. Close

Commit, open a PR and record decisions in `plans/README.md`. Track E ends here.

**Parked unless a reviewer asks:**
- fronts partitioned by division stage;
- a third, spatial-arrangement criterion;
- TAI per layer (Ma & Zheng 2023 phylostrata);
- reopening Track A with the shared criterion-ranked set;
- the mechanistic test from Schnabel et al.: each cell's distance to its fate
  region at the premorphogenetic stage.

## Implementation rules

- **Do not edit these hash-pinned modules:**
  - `pooled_analysis.py`, `publication_analysis.py`, `cousin_references.py`,
    `cross_species_analysis.py` and `round_references.py` in
    `full_tree_pareto/`;
  - every terminal module listed in Track D;
  - `representation_criteria.py`.

  Import them; do not change them.
- New modules write their own run folder, with source hashes, replay
  validation and tests.
- Add the bootstrap sampler to the guard in `publication/tests/support.py`.
- The output folders are links into the main checkout; never write one run
  folder from two sessions.
