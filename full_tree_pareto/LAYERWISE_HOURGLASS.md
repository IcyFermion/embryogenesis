# Per-layer lineage advantage and the mid-embryo travel dip (Track E, E3-E4)

Brief: `plans/track_e_layer_features.md`. Modules: `layer_hourglass.py` (E3) and
`layer_stages.py` (birth stages and the Large et al. JSD), with tests in
`test_layer_hourglass.py` and `test_layer_stages.py`. Figure candidates:
`publication/sketches/hourglass.py`.

Layers are numbered top-down in text and figures, **layer = 9 − round**. Code and
run outputs keep contraction rounds.

Notation follows `publication/notation.py`: d_LP is the endpoint-coordinate
distance from the lineage to P*. z_T is the natural travel cost minus the
first-cousin null mean, in null SDs (Track D), so negative values beat the
shuffle. The descriptive null-SD distance d̃_LP (persisted as `d_SD`) appears
only on the review page's background section and carries no claim.

## Method

Every cell of the full tree is measured, both its position and its protein
levels. The bottom-up tree therefore splits into independent per-layer
assignment problems: a layer's daughters are matched to its mothers' slots, and
every other layer stays fixed. Nothing is solved recursively, and E3 solves
nothing new.

- **Inputs.**
  - Pooled cosine run `pooled_full_tree_v1/terminal_clamped_cosine_20261014`:
    three C. elegans embryos, global top-20 proteins, 978 cells, 8 rounds.
  - Cross-species partial forest `cross_species_layerwise_v1/terminal_anchored_cosine_20261014`,
    raw3d geometry: C. elegans protein, C. elegans RNA and C. briggsae RNA
    (two AF16 embryos).
  - Both runs' per-round cousin references (`analysis/round_references`).
- **Advantage per layer**, computed separately for travel and for cell state:

  > advantage = 1 − natural cost / exact shuffle mean

  The shuffle permutes a layer's children within first-cousin groups (second
  cousins as a sensitivity check).
- **Intervals.** First-cousin groups are resampled with replacement within each
  layer, 2,000 draws. Children outside every group are singletons, whose
  natural cost equals their shuffle mean. Travel and cell state are resampled
  together. The partial forest's configurations share the resample indices.
- **Dip.** The mean advantage of layers 6-7 minus each flank. The flanks are
  layer 8, and the mean of layers 4 and 5.
- **Decision rule.** The dip is claimed only if both 95% intervals lie below
  zero, for pooled C. elegans travel against first cousins.
- **When it was fixed.** The contrast and the rule were fixed before any
  interval was computed, but **after** the point estimates had shown layers 6-7
  lowest in both species. A pass therefore shows that the dip is larger than
  resampling noise. It is not a confirmatory test, and only new embryos could
  confirm it.
- **Sensitivity**, reported but not gated:
  - the second-cousin shuffle;
  - the partial-forest configurations;
  - each tracked embryo's own travel matrix (stored in the cross-species run);
  - the cell-state dip.
- **Reproduction.** With no resampling, every per-round natural cost and
  first-cousin mean equals the cached round references at 1e-12.
  `--verify-only` replays every table and draw.

## Results (run `layer_features_v1/hourglass_20261019`, analysis id `a4ca3d85…`)

Travel advantage over the first-cousin shuffle, with 95% intervals:

| Layer (cells at birth, full tree) | Pooled C. elegans | Partial C. elegans | Partial C. briggsae |
|---|---|---|---|
| 8 (370) | 0.22 [0.19, 0.25] | 0.23 [0.18, 0.27] | 0.24 [0.18, 0.28] |
| 7 (192) | 0.12 [0.10, 0.15] | 0.09 [0.06, 0.12] | 0.08 [0.05, 0.12] |
| 6 (126) | 0.09 [0.06, 0.13] | 0.05 [0.03, 0.08] | 0.07 [0.04, 0.10] |
| 5 (73) | 0.15 [0.11, 0.20] | 0.16 [0.10, 0.23] | 0.19 [0.11, 0.26] |
| 4 (38) | 0.23 [0.16, 0.28] | 0.17 [0.07, 0.25] | 0.14 [0.04, 0.23] |

- **Gate passed.** The pooled C. elegans travel dip is −0.115 [−0.150, −0.077]
  against layer 8 and −0.082 [−0.126, −0.033] against layers 4-5.
- **Travel sensitivity.** The travel dip also excludes zero in every
  sensitivity series: the partial forest, C. briggsae, second cousins, and each
  of the five embryos on its own. The one exception is C. elegans embryo 1
  against layers 4-5, whose interval reaches zero.
- **Cell state.** The advantage is 0.02-0.06 in the pooled tree. No cell-state
  dip excludes zero in any series.
- **What changes in the gap** (figure candidates E and G; descriptive, no
  intervals). The advantage falls when the natural cost moves toward the
  shuffle cost: either cousins sit closer (the shuffle gets cheaper) or
  daughters end up farther from their own mother (the natural cost rises).
  - Travel to a random mother (full-random shuffle mean) is flat per edge,
    4.0-4.1 in layers 3-8 (C. briggsae 8.2-9.1), so per-edge costs compare
    across layers.
  - Against layer 8, the shuffle cost per edge is 67% and 38% higher at
    layers 6 and 7, so cousins are farther apart; the natural cost rises more
    (95% and 56%).
  - Against layers 4-5, layer 6's shuffle cost is unchanged (+1%) while the
    natural cost rises 13%. At layer 7 both fall, the shuffle more (−17%
    against −10%), so layer 7's gap to the early flank is partly cousins
    sitting closer.
  - C. briggsae repeats every contrast.

  So the dip is not explained by cousins sitting closer together, except in
  part for layer 7 against layers 4-5.
- **Front distance.** d_LP peaks in the same layers: 0.21 at layer 6 in the
  pooled tree. The hump survives the E1 per-round feature sets
  (`layer_features_v1/e1_fronts_20261015`).

Birth stages (run `layer_features_v1/stages_20261020`):
- Each layer's median number of cells alive at its cells' birth, from the
  Large et al. cell table: layers 8-4 are born at 370, 192, 126, 73 and 38
  cells.
- The same run recomputes the between-species JSD by stage. Its minimum is at
  the 200-cell stage, about 0.02-0.06 below the published values.

## Claims

| Claim | Wording |
|---|---|
| Lineage beats cousin shuffles on travel in layers 3-8 | Stated (z_T from −12.8 to −3.0, pooled; layer 3 added 2026-11-04) |
| It beats them by the least in layers 6-7 | Stated descriptively from the per-layer values and intervals (see "Presentation" below); the two contrasts appear once, in Methods, as a resampling check on layers chosen after looking |
| Cell-state advantage is small and flat | Stated (relative cost 0.94-0.98 in layers 4-8; layer 3's interval, 0.58-1.00, is wide) |
| The window matches fate-based cell sorting (Schnabel et al. 2006, 50 to ~400 cells) | "Consistent with" |
| It coincides with the expression hourglass waist (Large et al. 2025, 200 cells) | "Coincides with" |

Not claimed:
- the rank correlation with between-species JSD (seven points, both smooth in
  time);
- selection versus constraint, gene age (TAI), lethality and redundancy. These
  are background for one Discussion paragraph;
- anything about layers 1-2 (5 and 10 edges). Layer 3 (18 edges in 6
  first-cousin groups) is claimed from 2026-11-04; its front has only 5
  distinct points, so its d_LP of 0 is shown but not interpreted.

**Independence.** The C. briggsae embryos are independent, but their dip was
seen before the layers were fixed. The C. elegans partial forest reuses the
pooled embryos. So the single-embryo rows show consistency within the
available data, not replication.

## Figures (E4)

`python -m publication.sketches --output-dir publication/output/sketches/layer_hourglass_20261020 --only layer-hourglass`
builds five candidates in the production style:
- three main-text candidates:
  - A, four rows on a stage axis;
  - B, layer axis plus the dip;
  - C, a single headline panel;
- two supplement candidates:
  - D, every dip contrast;
  - E, the travel decomposition.

They are shown on the review page (https://claude.ai/artifact/U6vJK9ZwvpQGpkQJd6psHQ).

Compact set (2026-10-23, after discussion feedback and a collaborator review
the same day), from
`python -m publication.sketches --output-dir publication/output/sketches/layer_hourglass_20261023b --only layer-hourglass-compact`
(the version reviewed first is in `layer_hourglass_20261023`):
- F, main text: C. elegans full tree only, layers 3-8; d_LP and the travel
  and cell-state cost relative to the first-cousin shuffle (layer 3: 18
  cells in 6 cousin groups; its front has 5 distinct points and d_LP = 0);
- G, supplement: per-edge travel to the own mother, an aunt and a random
  mother, and the percent change in each dip contrast;
- H: approximate timeline against Large et al. Fig. 5E, with stages as
  birth-time windows and layers as birth distributions.

Compact review page: https://claude.ai/artifact/7hJxA3U7tP7TQJjaqiXBgj.

**Author decision (2026-11-04):** F goes to the main text and G to the
supplement; H is kept in reserve for a reviewer who asks about the timeline.
Candidates A-E are not used.

**Production (2026-11-04).** F and G were first built by the publication
family `full-tree-layers` as Figure 12 and supplement S5 (candidate
`track_e_layers_20261104`, not released).

**Author decision (2026-11-05): F becomes Figure 6B-C of the revised layout.**
The draft Figure 6 already showed d_LP per layer (old 6C), so a separate
Figure 12 would repeat it. Figure 6 keeps A (fronts of all eight layers;
layers 1-2 stay, to show why they are left out of the layer comparison) and
replaces the old B (lineage-similarity curves) and C (u and d_LP rows) with
two panels on the birth-size axis, layers 3-8:
- **6B, where the lineage sits:** d_LP, and the share of lineage edges
  rewired at the closest front point P* (1 − edge retention at P*; old 6B's
  quantity, one value per layer). Both peak at layer 6 (0.21; 0.82). The
  rewired share is descriptive, with no interval and no claim, and it is not
  independent of d_LP (same P*).
- **6C, lineage against the first-cousin shuffle:** the natural / shuffle
  cost for travel and cell state on one axis, with 95% bootstrap intervals.

The panels are drawn by `publication/figures/full_tree_layers.py`
(`compose_position`, `compose_cost`) from the adapter
`publication/adapters/full_tree.py:layers`. The adapter checks both runs'
file and source hashes and reproduces the per-layer costs against the
cached round references at 1e-12; it does not rerun the bootstrap. The
revised layout (`publication.v2`, branch figures-v2) places them. The family
`full-tree-layers` now holds only supplement S5 (G):
`python -m publication build --family full-tree-layers --output-dir publication/output/candidates/NAME`.

Caption draft for Figure 6B-C (to be merged into the Figure 6 caption):
"Layers 3-8, one point per layer at the median embryo size (cells alive)
when its cells are born (cell table of Large et al. 2025). Layers 1-2 (5 and
10 cells) are left out: their fronts have 2 and 6 distinct points (A).
(B) Distance d_LP from the lineage to the layer's front, and the share of
lineage edges rewired at the closest front point P*. Layer 3's front has
only 5 distinct points and the lineage lies on one of them, so both are 0.
(C) Natural travel and cell-state cost divided by the mean cost of the exact
first-cousin shuffle, which reassigns daughters among mothers that share a
grandmother; 1 means no better than the shuffle. Error bars are 95% bootstrap
percentile intervals from 2,000 resamples of first-cousin groups within each
layer (layer 3: 6 groups); B has no intervals."

### Presentation (collaborator review, 2026-10-23)

- **Orientation.** Panels b-c plot natural / first-cousin shuffle cost
  (1 − advantage), so they peak in layers 6-7 like d_LP; 1 means no better
  than the shuffle. 1/d_LP was ruled out because d_LP is 0 at layers 2-3. The
  JSD in H keeps its own direction; only timing is compared.
- **Error bars.** 95% bootstrap percentile intervals (2,000 resamples of
  first-cousin groups within a layer), not SD or SEM. d_LP has none. In H,
  horizontal bars and JSD boxes show spread (middle 50%), not uncertainty.
- **Layers 1-2** are not drawn.
- **Dip contrast.** Out of the figure and the main-text results. Without it,
  readers judge the early side from bars alone, and layer 5's interval
  overlaps layers 6-7 (layers 4 and 8 do not). So it stays as one Methods
  sentence, in relative-cost terms (mid minus flank): +0.115 [0.077, 0.150]
  against layer 8 and +0.082 [0.033, 0.126] against layers 4-5; the cell-state
  contrasts' intervals include zero. The main-text wording leans on what the
  bars show: layers 6-7 sit above layers 4 and 8, with layer 5 between.
