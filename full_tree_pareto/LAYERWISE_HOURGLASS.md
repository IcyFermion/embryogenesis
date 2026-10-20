# Per-layer lineage advantage and the mid-embryo travel dip (Track E, E3-E4)

Brief: `plans/track_e_layer_features.md`. Modules: `layer_hourglass.py` (E3) and
`layer_stages.py` (birth stages and the Large et al. JSD), with tests in
`test_layer_hourglass.py` and `test_layer_stages.py`. Figure candidates:
`publication/sketches/hourglass.py`.

Layers are numbered top-down in text and figures, **layer = 9 − round**. Code and
run outputs keep contraction rounds.

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
- **What changes in the gap.** Natural travel per edge against the shuffle mean
  (figure candidate E):
  - at layer 6 the natural lineage's daughters end up farther from their
    mother's last tracked position, while the shuffle cost barely changes;
  - at layer 7 both fall.

  So the dip is not only cousins sitting closer together.
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
| Lineage beats cousin shuffles on travel in layers 4-8 | Stated (travel z 3.6-12.8, pooled) |
| It beats them by the least in layers 6-7 | Stated as larger than resampling noise; layers chosen after looking |
| Cell-state advantage is small and flat | Stated |
| The window matches fate-based cell sorting (Schnabel et al. 2006, 50 to ~400 cells) | "Consistent with" |
| It coincides with the expression hourglass waist (Large et al. 2025, 200 cells) | "Coincides with" |

Not claimed:
- the rank correlation with between-species JSD (seven points, both smooth in
  time);
- selection versus constraint, gene age (TAI), lethality and redundancy. These
  are background for one Discussion paragraph;
- anything about layers 1-3 (5-18 edges, degenerate fronts).

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
Production waits for the author's choice.
