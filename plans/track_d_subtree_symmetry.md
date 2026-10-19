# Track D: left-right subtree symmetry

**Question.** Bilaterally symmetric sublineages should behave alike in the
terminal Pareto analysis. Examples are ABplpa ↔ ABprpa and MSaa ↔ MSpa. Do
mirrored subtrees give similar front geometry and natural-lineage positions?
Where they differ, is the cause biology (known left-right asymmetric fates),
the data (coverage, tracking) or numerics (tiny endpoint spans)?

**Extension (2026-10-19).** The ABplpa/ABprpa follow-up (D3) traced the
travel span to an apoptotic cell. D4 asks how apoptotic cells should be
handled in the terminal analysis as a whole. It stays in this track at the
author's request, although it reaches beyond symmetry.

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

## D3. Trajectories of ABplpa and ABprpa (2026-10-19)

Requested after the collaborator review: is ABplpa's tiny travel span a
tracking defect? Scratch analysis on the stored runs; no run written yet.
Review page: https://claude.ai/artifact/UcSsS7P1DehMnxUUVMw1Cm (animation,
mirror-distance chart, swap diagram, track checks).

**Tracks are clean.**
- No missing or duplicated frames in either subtree, in any embryo.
- Steps above the 99.9th percentile are as common on the left as on the right.
- The only track ending before the cutoff is ABplpappap (apoptotic), in
  embryos 2 and 3.
- Reflected through a per-frame midline plane (fitted on the other ABpl/ABpr
  pairs), left cells fall within the reference residual. The exceptions repeat
  in all three embryos: the ABxpapp branch (different fates) and the early
  ABxpaa cells.
- Division timing, left minus right, at the last round: +4.5, −9.5 and −2
  frames in embryos 1–3. No side is consistently late. The embryo 1 lag pushes
  SIAD, SIBV, AMSO and CEMV past the cutoff.

**The travel span is one swap.**
- In ABplpa, natural = state optimum in every embryo.
- The travel optimum swaps the apoptotic ABplpappap with RMEV (EXC GL L in
  embryo 2), saving 1.83, 1.70 and 0.11 µm. That is the whole span.
- Without ABplpappap, the front is a single point (natural = both optima).
- In ABprpa, natural = travel optimum (pooled, embryo 2). Its state optimum is
  a near-tie in protein expression: a five-cell reshuffle changes summed
  cosine distance by −0.027 on the right and +0.49 on the left; SAAVR is the
  largest term (0.77 against SAAVL's 0.45).
- On the 15 name-matched cells, ABplpa is still a point and ABprpa keeps a
  0.18 SD state optimum.

**Conclusion.** No tracking defect. The u gap is numerical (confirms D2).
Report z for this pair.

**To do if used in the paper:** promote the analysis to a
`subtree_symmetry.py` stage with a write-once run (per-frame planes, mirror
distances, track checks, swap table), and decide whether the mirrored overlay
becomes a supplementary movie or a snapshot panel.

## D4. Apoptotic cells

**Question.** How should apoptotic cells enter the terminal analysis? Death
is an invariant, lineage-determined fate, so it is informative for cell state.
A dying cell has no destination, so its travel endpoint is poorly defined.

**Scratch evidence (2026-10-19, pooled canonical cohort).**
- 30 of the 275 cohort cells are apoptotic. 13 sit in first-cousin groups
  with another parent and so contribute to the null.
- Endpoints are inconsistent between embryos. 12 tracks end before the cutoff
  in some embryos but not others; 1 ends early in all three and 17 in none.
  Lifespan to the end of the track or the cutoff ranges from 1 to 92 frames.
- Natural travel is slightly longer than for living cells (median 2.6–3.4
  against 2.3 µm). Natural state cost is higher (cosine 0.53 against 0.36).
- Whole-cohort front, pooled:

  | Variant | d_SD | travel z | state z | retention at nearest |
  |---|---|---|---|---|
  | As is | 1.12 | −9.09 | −4.83 | 0.79 |
  | Drop deaths | 0.57 | −8.48 | −4.33 | 0.87 |
  | Drop 30 random living cells (median of 40) | 0.84 | −8.27 | −4.43 | 0.83 |
  | Deaths travel-neutral (constant column) | 1.89 | −8.86 | −4.83 | 0.75 |

  Dropping deaths brings natural closer to the front than 37 of 40 random
  drops. The travel-neutral variant is not neutral: it lets the optimum park
  costs on dead cells for free, so it is ruled out.

**Steps.**
1. **Characterize the 30 cells.**
   - Track length after the expected death: do step sizes change, which would
     mean corpse tracking?
   - How many frames each protein profile rests on.
   - Whether the same cells die in C. briggsae.
   - `data/apoptotic_cells.txt` lists names only (124). Ask the author whether
     death times exist; with them, the endpoint can be set at death.
2. **Variants**, each on the whole cohort, the 42 subtree fronts, the
   fate-only cost and the 187-edge cross-species cohort:
   - as is;
   - drop deaths;
   - endpoint at birth (deaths do not migrate; same duration in every embryo);
   - endpoint at death, if death times exist.
3. **Control.** About 1000 removals of living cells matched to the deaths on
   depth and birth time. The 40-draw result is a hint only.
4. **Decision rule (set before the run).**
   - If drop deaths beats the matched control on d_SD and retention
     (one-sided, 5%), report it as a finding: the natural lineage departs from
     the front mostly at cells that die.
   - Otherwise keep deaths in the canonical analysis and show the
     drop-deaths variant in the supplement.
   - In either case, check whether any D1/D2 classification changes.

**Implementation.** A new module `terminal_pareto/apoptosis_sensitivity.py`
writing `subtree_symmetry_v1/apoptosis_<date>/`. It reads the canonical
context, builds every front through `cross_species_analysis.analyze`, gates
the as-is variant against the canonical metrics, and adds its sampler to the
guard in `publication/tests/support.py`. The canonical runs stay unchanged
until the author decides.

## D4 results (2026-10-19)

Run `subtree_symmetry_v1/apoptosis_20261019`
(`terminal_pareto/apoptosis_sensitivity.py`), replay-verified.

**Gates.**
- The per-embryo travel matrices rebuilt from the tracks match the canonical
  components exactly.
- The as-is whole-cohort front reproduces the Track A criteria front
  (measured and fate-only, ≤ 5e-15).
- The as-is pair comparison reproduces D1.

**Death times.** None are in the data, and a literature search found no
downloadable per-cell table:
- Sulston et al. 1983 gives approximate times only in its lineage diagrams
  (most deaths 20–30 min after birth).
- The Richards et al. 2012 AceTree model file is offline (404).
- A 2026 Cell Death & Differentiation paper times corpse formation for 14
  deaths, mostly in a figure (ABalapapaa: 13 min after birth).

Nucleus radius is quantized and noisy, and the dying daughter is smaller than
its sister from birth, so it does not mark death either. The two variants
bracket any death-time endpoint: `birth_endpoint` (dies at once) and `as_is`
(end of the corpse track).

**The 30 deaths.**
- Tracks end before the cutoff for 4, 10 and 7 deaths in embryos 1–3.
- Median 28–35 frames from birth to the endpoint used. If a frame is about a
  minute, this is beyond the literature's 13–30 min, so many endpoints are
  corpse positions.
- The birth-to-endpoint shift is median 2.9 µm (0.2–9.9).
- Controls: every death's 8 nearest candidates have its depth; median
  birth-frame gap 1.75 (MSpaapp, born early, 22.5).

**Whole cohort, pooled, measured cost** (1000 matched draws):

| Variant | d_SD | Retention at nearest | Travel z | State z |
|---|---|---|---|---|
| As is | 1.11 | 0.79 | −9.09 | −4.83 |
| Drop deaths | 0.57 (control median 0.59, p = 0.47) | 0.87 (0.88, p = 0.67) | −8.48 | −4.33 |
| Endpoint at birth | 0.49 (0.57, p = 0.31) | 0.88 (0.88, p = 0.47) | −9.18 | −4.83 |

**Decision rule: no finding.** Deaths are not special. Removing 30 living
cells matched on depth and birth frame helps just as much, so the departure
of the natural lineage from the front sits in late-born, deep cells in
general. The unmatched 40-draw hint (37/40) was a birth-time confound.
Deaths stay in the canonical analysis; `drop_deaths` goes to the supplement.

**Subtrees** (pooled d_SD, as is → drop deaths / endpoint at birth):
- In the posterior ABpl subtrees, the distance to the front is entirely the
  deaths' endpoints. ABp 0.47, ABpl 0.54, ABplp 1.16, ABplpp 1.67 and ABplppp
  1.67 all go to 0 in both variants. (ABplppp drops below 12 edges without its
  death.)
- ABarp goes 0.78 → 0 when dropped but stays 0.73 with the birth endpoint.
- The whole-cohort rule above does not test subtrees; no matched control was
  run at subtree level.

**Pair comparison.**
- ABplp/ABprp: the D2 differences (d_SD, d_LP) vanish under both variants. Its
  "data" class now has a concrete cause, the apoptotic endpoints on the
  ABplp side.
- ABalpa/ABarap (exploratory): a d_NP difference appears. That is an endpoint
  metric, so it is not followed up.
- Nothing else changes status.

## Decisions (2026-10-19)

- Track D absorbs the trajectory follow-up (D3) and the apoptosis question
  (D4) rather than opening a new track.
- A symmetry part goes into the publication; its form (table, figure, movie)
  is still open.
