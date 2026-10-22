# Track E: layer-specific expression features in the full tree

**Question.** The pooled full-tree layerwise assignment uses one feature set,
the top-20 proteins chosen on terminal cells, in every contraction round.
Earlier cells may be told apart by different factors than terminal fates.
Is a different feature set per layer more reasonable, and does it change the
fronts?

**Layer convention.** Analysis code uses contraction rounds
(`pooled_analysis.contraction_layers`), which count bottom-up. Round 1 holds the
491 terminal edges and round 8 holds 5 edges. The v2 Figure 6 numbers layers
top-down. Keep rounds in all analysis output and convert only in figures.

## Principles (from Track A)

See `terminal_pareto/METHODS_cell_state_representation.md`.
- **The front never selects features.** Per-round features are chosen by a
  criterion that uses neither travel nor the front. Fronts are a robustness
  check only.
- **The single global top-20 is the baseline.** A per-round scheme must beat
  it on that criterion, and by more than random sets of the same size do. Keep
  the options few: one k for all rounds and one ranking method.
- **Cosine cell-state distance**, as decided in Track A2.

## Existing material

- **Pooled cosine run** `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_cosine_20261014/`:
  layerwise sweeps, plus per-round references (`analysis/round_references/`:
  fronts, metrics and draws per round).
- **Protein values for internal cells.** `data/protein/aggregated_all/s3_zscore.csv`
  covers internal cells (210 features); raw values are in `s3.csv`. Start with
  protein. RNA is resolved by cell type and covers only the 485-cell partial
  forest.
- **Earlier time-resolved feature work** in `expression_embedding/`:
  - `timepoint_terminal_feature_select.ipynb` and
    `stepped_avg_feature_select.ipynb`;
  - `timepoint_embedding.py` and `search_timepoint.py`.

  Check what they selected, and with which labels, before building anew.
- **Reusable machinery.** `terminal_pareto/representation_criteria.py` has
  sibling separation, lineage memory (within-type shuffle) and the fate-only
  cost.

## E0. Labels and criterion (do first)

- **Labels.** Internal cells have no fate. Proposed label: the composition of
  terminal cell types among a cell's descendants (fine or Schema D types).
- **Criterion per round.** A generalized sibling separation: do sisters whose
  descendant compositions differ have more distant expression than sisters
  whose compositions agree? Measure, for example, the rank correlation of
  expression distance with the Jaccard or Bray-Curtis distance between
  compositions. Agree the exact form with the author.
- **Expression coverage per round.** Report the fraction of cells with a
  variable reporter signal. Early rounds may carry nearly silent panels, which
  makes every criterion uninformative there.
- **Lineage memory per round.** Run the within-type-shuffle check per round,
  using descendant compositions. Daughters inherit their mother's protein, so
  memory is likely strongest in early rounds.

## E1. Per-round feature ranking

- For each round, rank features by the E0 criterion, or by a per-round
  classifier that predicts descendant composition, and take the top k.
- Compare each per-round set with the global top-20 and with random sets from
  that round's highly variable pool. Report overlap with the global set and
  stability under bootstrap.
- **Decision gate.** Proceed to E2 only if per-round sets beat the global set
  clearly on the criterion. Otherwise the global set stands (Occam).

## E2. Fronts (robustness)

- Rebuild the layerwise sweep with a per-round state matrix.
- **Scaling decision.** The pooled analysis scales the optimizer with one
  global first-cousin SD shared by all rounds. Per-round features change the
  state scale per round, so choose between that global scaling and per-round
  scaling; the per-round cousin nulls in `round_references` can supply it.
- Compare per-round metrics with the global-feature run, and include the
  fate-only control per round.
- **Cost.** Layerwise sweeps take seconds each; per-round references about a
  minute.

## Implementation

- Add a new module (proposed `full_tree_pareto/layer_features.py`) that
  imports the pooled functions. Write
  `full_tree_pareto/output/runs/layer_features_v1/<run_id>/`, with source
  hashes, replay validation and tests.
- **Do not edit these hash-pinned modules:**
  - `pooled_analysis.py`, `publication_analysis.py`, `cousin_references.py`,
    `cross_species_analysis.py` and `round_references.py` in `full_tree_pareto/`;
  - every terminal module listed in Track D;
  - `representation_criteria.py`, whose own run is pinned to it. Import it; do
    not change it.
- Add any new solver or sampler to the guard in `publication/tests/support.py`.

**Decisions for the author:** the descendant-composition label and the
criterion (E0); k, and whether per-round sets are allowed to differ in size;
and global against per-round scaling (E2).
