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

## Decisions (2026-10-15)

- **D0 result.** Run `subtree_symmetry_v1/pairing_20261014`
  (`terminal_pareto/subtree_symmetry.py`). Twelve candidate pairs: the nine
  brief pairs, plus ABa/ABp, ABalp/ABara and ABalpa/ABarap. The last two are
  the anterior counterparts; they pair by fate, not by lineage position.
- **Threshold (strict).** A pair is kept only if:
  - each subtree is the other's best fate match;
  - the left-right (sided) match is at least 0.85;
  - both partners have at least 12 cohort edges.
  This keeps ABplp/ABprp, ABplpa/ABprpa, ABplpp/ABprpp, ABplppp/ABprppp,
  MSa/MSp, MSaa/MSpa and MSap/MSpp.
- **Priority.** ABplpa/ABprpa and MSaa/MSpa come first. Experts in C. elegans
  embryo development suggested them.
- **Secondary.** The anterior pairs (ABalp/ABara, ABalpa/ABarap) are
  exploratory. Run them after the strict pairs, but they are not the main
  question.
- **Primary readout.** d_SD. Endpoint metrics (u, d_LP, d_CP) are reported
  alongside it.

## D1 results (2026-10-15)

Runs `subtree_symmetry_v1/compare_20261015` (fronts) and `cells_20261015`
(cell level), both replay-verified. Every pooled subtree front reproduces the
canonical table to 3e-12.

- **d_SD saturates.** The natural lineage lies on its own subtree front
  (d_SD = 0) in 26 of 42 subtrees, including all four priority subtrees. The
  continuous readout in the same null-SD frame is the lineage's position
  against first-cousin shuffles: travel z and state z.
- **Null support.** Many cohort cells have only their sister in the cohort, so
  their first-cousin group cannot change a cost. MSpa's null rests on one
  three-cell group (3 of 12 cells); MSaa's on 7 of 16.
- **ABplpa / ABprpa.** Similar.
  - Both have d_SD = 0. State z is −2.70 against −2.83, closer than 90% of
    size- and depth-matched random pairs. Travel z is −2.54 against −3.13,
    a typical difference.
  - u differs (1 against 0), but only because ABplpa's travel span is
    0.04 null SD. This is the expected numerical artifact.
  - Counterpart cells are the closest cell of their fine type in 87% of
    pairs.
- **MSaa / MSpa.** Different on the state axis: state z −1.70 against +0.92,
  at the 96th percentile of random pairs. The data cause this:
  - MSpa's four late-born pharyngeal cells are missing (born t 262–278);
  - MSpa's null rests on one group;
  - at the mirror position of MSaa's valve cells (VPI2V, VPI3V), MSpa makes
    two body-wall muscle cells, a real fate asymmetry. These are sister-only
    groups, so they do not move z.
  - Counterparts are no closer in expression than other pharyngeal muscle
    cells (closest of type in 11% of pairs).
- **Other kept pairs.** Within-pair state z differences sit at the 12th–47th
  percentile of random pairs. d_SD differences (ABplp 1.16, ABplpp 1.67) are
  within the between-embryo floor (1.13, 2.28).
- **Embryo frame.** In each embryo, one reflection plane maps the 55 sided
  cohort pairs onto each other. The residual is 2.0–2.5 µm, against
  11–13 µm left-right separation.
  - Counterpart edge lengths agree to within the between-embryo floor, except
    MSaa/MSpa (0.96 µm against 0.57).
  - Mirrored displacement angles are 20–54°, against about 76–94° for
    non-counterparts.
- **Expression symmetry by founder.** Posterior AB counterparts are nearly
  identical (closest of their fine type in 64–100% of pairs). MS counterparts
  are not (11–22%).

Caveat: state z has no tracking floor (one protein dataset), so its reference
is the random pairs only.

## D2 results (2026-10-15)

Run `subtree_symmetry_v1/robustness_20261015`
(`terminal_pareto/subtree_symmetry_robustness.py`), replay-verified.

**Variants.** Every subtree with ≥ 12 edges is recomputed:
- in XY;
- with the fate-only cost;
- in the 187-edge cohort: CE protein, CE RNA and CB RNA, each in 3D and XY;
- on a counterpart-matched cohort.

Each variant reproduces its cache at P0 to ≤ 1e-16. That covers the Track A
fate-only front, the cross-species metrics, and the XY matrix against the
187-edge cache.

**Classification** (rules in the module docstring):

| Pair | Tier | Class | Differs (3D) | Notes |
|---|---|---|---|---|
| ABplpa / ABprpa | priority | numerical | u | Travel span 0.04 null SD. The u gap shrinks in XY (0.61) and fate-only (0.38) and vanishes cross-species. State z matches in CE RNA (−0.70 / −0.73); CB RNA differs at an ordinary level (−0.62 / +0.07, 26th percentile). Travel z matches in C. briggsae (−2.17 / −2.23). |
| MSaa / MSpa | priority | data | state z | MSpa's null rests on one group (support 0.25), cohort ratio 0.75, late-born pharyngeal cells missing. Fate-only: MSpa's null is degenerate, and MSaa's state z goes from −1.70 to +1.00, so its measured advantage is lineage memory within pharyngeal muscle. Not in the 187 cohort. |
| ABplp / ABprp | strict | data | d_SD, d_LP | Gone in the matched cohort; late-born gap 9 cells. |
| ABplpp / ABprpp, ABplppp / ABprppp, MSa / MSp, MSap / MSpp | strict | symmetric | — | — |
| ABalp / ABara, ABalpa / ABarap | exploratory | symmetric | — | Under fate-only only: state z −2.04 / +1.58 and −1.69 / +1.38 (94th–97th percentile). A lead, not classified. |

**No pair classifies as biological.** The fate asymmetries are real but do not
move the front-level readouts:
- ABplpa: excretory cell, duct and gland, K, rectal cells, PVT.
- ABprpa: AVL, G1, RIH, RIR, RIS, AVG, VIRL/VIRR.
- MSaa/MSpa: VPI against body-wall muscle.
- Known asymmetric pairs: AWCL/AWCR, and ASEL in ABalp.

The per-pair lists in `classification.csv` are the input for the literature
pass in the Track C research session.
