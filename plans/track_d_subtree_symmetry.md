# Track D: left-right subtree symmetry

**Question.** Bilaterally symmetric sublineages should behave alike in the
terminal Pareto analysis. Examples are ABplpa ↔ ABprpa and MSaa ↔ MSpa. Do
mirrored subtrees give similar front geometry and natural-lineage positions?
Where they differ, is the cause biology (known left-right asymmetric fates),
the data (coverage, tracking) or numerics (tiny endpoint spans)?

## Existing material

- **Per-subtree canonical metrics** for the 42 subtrees with at least 12 edges,
  in `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920/analysis/`:
  - `ce_subtree_canonical_metrics.csv`: u, d_LP, d_CP (`d_np`), maximum
    retention, endpoint spans, on-front flags;
  - `ce_subtree_canonical_curves.npz` and `subtree_manifest.csv`.
- **Candidate pairs already in that table:** ABpl/ABpr, ABpla/ABpra,
  ABplp/ABprp, ABplpa/ABprpa, ABplpp/ABprpp, ABplppp/ABprppp, MSa/MSp,
  MSaa/MSpa and MSap/MSpp. The ABal/ABar descendants need the pairing step
  below.
- **Overlap with Track C.** ABplpa (u = 1, d_CP 64) and ABprpa (d_CP 3.4) both
  have near-zero endpoint spans. Symmetry and outliers overlap, so share the
  Track C dossier tooling where it exists.
- **Cell annotation.** The cell key `data/2023-06-29_entropy_cell_key_V2.csv`
  has WormBase cell names with L/R suffixes (`wormweb.cell`) and fine types
  (`wormweb.desc`). `terminal_pareto/representation_criteria.py` loads both.
- **Robustness caches.** Per-embryo tracking sensitivity (`tracking_sensitivity/`,
  `s3_cross_geometry/`), and the 187-edge cross-species cohort (CE protein, CE
  RNA, CB RNA).
- **Conventions** from Track A (`terminal_pareto/METHODS_cell_state_representation.md`):
  - cosine cell-state distance and the published top-20 features;
  - d_SD (null-SD distance) as the sensitive readout;
  - the lineage-memory and fate-only machinery.

## D0. Pairing (do first)

Names are only a starting list. In the anterior AB lineage, left-right
counterparts come from branches that do not correspond by name.
1. Map each terminal cell to its bilateral counterpart from the WormBase cell
   names (…L/…R). Record unpaired cells, cells with left-right asymmetric fates
   (e.g. ASEL/ASER) and cell deaths.
2. Pair subtrees whose terminal descendants map onto each other. Report the
   matched fraction for every candidate pair, and keep pairs above a threshold
   agreed with the author.
3. Check cohort coverage. Both partners should have comparable numbers of
   tracked edges. Different tracking cutoffs or missing cells are a confound.

## D1. Compare paired subtrees

- **Subtree metrics.**
  - u, d_LP, d_CP, maximum retention and on-front status.
  - Natural-lineage travel and state z in each subtree's first-cousin null.
  - d_SD and the endpoint spans.
  - Raw costs.
  - Where a span is tiny, u and d_CP are unstable, so use the raw and null-SD
    metrics there.
- **Cell-level checks.**
  - Expression distance between counterpart cells, against same-type cells
    that are not counterparts.
  - Travel geometry: lengths of counterpart edges, and whether the embryo
    frame is mirrored.
- **What counts as similar.** The difference within a mirrored pair is judged
  against two references:
  - differences between random subtree pairs matched for edge count and depth;
  - differences between tracking embryos for the same subtree, which gives the
    measurement floor.

## D2. Explain the differences

Classify each asymmetric pair:
- **data:** coverage, cutoffs, z-axis uncertainty, one embryo only;
- **numerical:** endpoint spans;
- **biological:** known asymmetric fates or deaths.

Check robustness per embryo, in 3D against XY, with the fate-only cost, and
in the cross-species cohort where both partners survive the 12-edge
threshold. A literature pass on the biological candidates can be folded into
the Track C research session.

## Implementation

- Add a new module (proposed `terminal_pareto/subtree_symmetry.py`) writing
  `terminal_pareto/output/runs/subtree_symmetry_v1/<run_id>/`. Follow the other
  add-ons: write-once runs, source hashes, replay validation, tests.
- Read the canonical caches rather than recompute them. Any new front (e.g.
  for d_SD) goes through `cross_species_analysis.analyze`.
- **Do not edit these hash-pinned modules:**
  - `cross_species_analysis.py`, `subtree_analysis.py`, `front_coordinates.py`,
    `analysis_context.py`, `data_loader.py`, `lineage_metrics.py`,
    `pareto_engine.py`, `global_analysis.py`,
    `fig5_table1_ce_canonical_metrics.py` and `figS3_cross_geometry.py`;
  - `feature_sensitivity.py` and `representation_criteria.py`, which pin their
    own runs.
- Add any new solver or sampler to the guard in `publication/tests/support.py`.

**Decisions for the author:** the pairing threshold; whether to report d_SD
or endpoint metrics as primary for pairs with tiny spans; and whether the
results become a figure or a supplementary table.
