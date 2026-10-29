# Track D: left-right subtree symmetry

**Status (2026-10-27): reopened for D5** (readouts that one cell or a
collapsed span cannot decide; see D5 at the end). Branch `subtree-symmetry`,
recreated from `main` at `2542a03`.

D0–D4c were closed on 2026-10-22 and merged to `main` (PR #6). Each is a
replay-verified run under `subtree_symmetry_v1/`. Everything below up to D5
describes that state.

In the paper:
- Figure 5B numbers the mirrored pairs, and new panel 5C shows the pairs.
  The caption carries notes † ‡ §.
- Release candidate: `publication/output/candidates/fig5_symmetry_20261022`.
- Supplementary Note A is kept in reserve, outside the registry, at
  `publication/output/notes/noteA_20261022/`.

Handed off:
- **Figure 5 layout** (`figures-v2`):
  - The two-page layout lasts until panel A moves to Figure 1; then drop
    `PAGES` in `captions/terminal.py`.
  - v2's Figure 4B and 4C reuse `publication/figures/terminal_symmetry.py`.
- **z names:** with the author's colleagues; the change is one line in
  `notation.py`.
- **Literature pass:** fate-unmatched cells (`classification.csv`) go to
  Track C.
- **Not pursued:**
  - a D3 mirrored-trajectory movie (the analysis stays scratch; see the
    review page);
  - a later tracking cutoff for MSpa;
  - the anterior fate-only lead (ABalp/ABara).

**Question.** Bilaterally symmetric sublineages should behave alike in the
terminal Pareto analysis. Examples are ABplpa ↔ ABprpa and MSaa ↔ MSpa. Do
mirrored subtrees give similar front geometry and natural-lineage positions?
Where they differ, is the cause biology (known left-right asymmetric fates),
the data (coverage, tracking) or numerics (tiny endpoint spans)?

**Extension (2026-10-19).** The ABplpa/ABprpa follow-up (D3) traced the
travel span to an apoptotic cell. D4 asks how apoptotic cells should be
handled in the terminal analysis as a whole. It stays in this track at the
author's request, although it reaches beyond symmetry.

**Extension (2026-10-27).** D5 reopens the track. Figure 5C does not make
the case for symmetry at a glance, because the core readouts are fragile in
tightly bound subtrees. D5 looks for readouts that hold up there.

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
Deaths stay in the canonical analysis. The rule sent `drop_deaths` to the
supplement, but the author kept it to this brief (decision 2026-10-21).

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

## Decisions (2026-10-21)

- Apoptotic cells stay in the canonical analysis as before. The D4 variants
  are not added to the supplement. They are kept here so the manuscript can
  answer a reviewer's question about them.
- The subtree-level effect in posterior ABpl gets a matched control (D4b).
  After that comes the question of what to show in the paper.
- Figure sketch on the review page (https://claude.ai/artifact/UxCTYL2LAUf5cYCPLABDi1,
  "Figure sketch"). It has four panels: (a) pair coverage; (b) both fronts
  per pair in the null-SD frame; (c) left against right z; (d) reference
  percentiles. Main, Figure 4 row or supplement is still open.
- Table alternative, preferred by the author (same page, "Table sketch").
  Number the 18 paired subtrees in Figure 5B. Add a panel 5C with one row per
  pair: u, d_LP, d_SD, z_S, z_T, with left and right on the same row, cells
  at or above the 90th reference percentile shaded, and the D2 class. A side
  note covers the standouts:
  - † ABplpa/ABprpa u: one apoptotic-cell swap, travel span 0.04 null SD;
  - ‡ MSaa/MSpa u and z_S: missing late-born cells, one-group null;
  - § ABplp/ABplpp/ABplppp d_LP: the ABplpppapp corpse endpoint.
  Not yet built in `publication/`.
- Notation: the author wants d_SD shown as a decorated d_LP. Suggested
  `\tilde{d}_{LP}`, with "nearest front point" in the definition, because in
  null-SD units the nearest point can differ from P*. Not yet in
  `publication/notation.py`. The z names (sketch: z_T, z_S) go to the
  author's colleagues.

## D4b. Subtree-level control (set 2026-10-21, before the run)

**Question.** Inside the posterior ABpl subtrees, d_SD goes to 0 when the
deaths are dropped or held at birth. Is that because they are deaths, or
would removing any equally many late-born cells of the same subtree do the
same? Small subtrees have few cells near the front, so a large change from
removing a few cells is expected.

**Control.** As in D4, but within the subtree. Each draw picks one living cell
of that subtree per death, at random among the `POOL_SIZE` nearest free
candidates (same depth first, then birth frame), and applies the same
operation: drop, or endpoint at birth. 1000 draws per subtree and variant,
pooled travel, measured cost. Small subtrees have few distinct control sets;
their number is reported with each p-value.

**Subtrees tested.** Every canonical subtree with at least one death and a
pooled as-is d_SD above 1e-9, except P0 (done in D4). Ties in the p-values are
counted within the same 1e-9.
- **Primary:** ABplp, the left partner of the D2 pair ABplp/ABprp.
- **Supporting:** ABp, ABpl, ABplpp, ABplppp (nested with ABplp, so they
  share its deaths). ABplppp without its death drops below 12 edges, so only
  its birth-endpoint variant is tested.
- **Descriptive:** the rest, including ABarp.

**Decision rule.**
- **Finding:** in ABplp, drop deaths beats the matched control on d_SD
  (lower) and on retention at the nearest point (higher), each one-sided at
  5%. Then the ABplp/ABprp difference has a specific cause: the deaths'
  endpoints. Endpoint at birth is reported as the check.
- **Otherwise:** the departure is a small-subtree effect that any matched
  removal reproduces. ABplp/ABprp stays "data", with no cause specific to
  apoptosis.

**Implementation.** `terminal_pareto/apoptosis_subtree_control.py` reuses the
D4 context and sampler and writes `subtree_symmetry_v1/apoptosis_subtrees_<date>/`.
Its sampler joins the guard in `publication/tests/support.py`.

## D4b results (2026-10-21)

Run `subtree_symmetry_v1/apoptosis_subtrees_20261021`
(`terminal_pareto/apoptosis_subtree_control.py`), replay-verified. The as-is and variant
fronts reproduce D4's subtree table (d_SD ≤ 1e-16). 14 subtrees tested.

**Decision rule: finding.** In ABplp, dropping the three deaths takes d_SD
from 1.16 to 0 and retention from 0.91 to 1. None of 1000 matched draws (119
distinct sets) gets there: control d_SD median 1.22, 5th percentile 1.15;
p = 0.001 on both metrics. Endpoint at birth gives the same (p = 0.001).

| Subtree | Role | Drop deaths: d_SD p / control at 0 | Endpoint at birth: d_SD p / control at 0 |
|---|---|---|---|
| ABplp | primary | 0.001 / 0% | 0.001 / 0% |
| ABplpp | supporting | 0.24 / 24% | 0.001 / 0% |
| ABplppp | supporting | below 12 edges | 0.001 / 0% (8 distinct sets) |
| ABpl | supporting | 0.11 / 11% | 0.11 / 11% |
| ABp | supporting | 0.14 / 14% | 0.09 / 9% |
| ABarp | descriptive | 0.24 / 24% | 0.22 / 22% |

Retention p equals d_SD p in every row above.

**One death carries it** (scratch, same context). Dropping or holding
ABplpppapp at birth alone gives d_SD 0 in ABplp; dropping the other two deaths,
alone or together, leaves 1.16. Its mirror ABprpppapp also dies, but ABprp is
on its front. The left cell's endpoint is farther from its parent in every
embryo (8.9, 4.8 and 7.0 µm against 5.7, 3.4 and 4.4) and its state cost is
higher (cosine 1.21 against 0.97). Its tracks run 49–69 frames past birth,
past the literature's death times, so the endpoint is likely a corpse
position. Endpoint at birth removes it, so travel, not state, is the cause.

**Caveat: matching is weak inside ABplp.** It has no living depth-9 cohort
cells, so every candidate is one division deeper and born 23–68 frames later.
The control asks "any three other cells of the subtree", not "three cells
like the deaths".

**Reading.**
- The ABplp/ABprp difference (D2: d_SD, d_LP; class "data") has a specific
  cause: the corpse endpoint of ABplpppapp. The class stays "data".
- ABplpp and ABplppp, nested on the same cell, agree with endpoint at birth.
  Dropping 2 of 15 cells from ABplpp is too coarse a test (24% of controls
  also reach 0).
- In ABpl and ABp the effect is diluted: 9–14% of matched removals also bring
  the lineage onto the front.
- Descriptive: in MS, EMS, P1 and MSp the deaths bring the lineage closer to
  the front (dropping them raises d_SD, e.g. MS 0.20 → 0.61). In AB, holding
  deaths at birth beats the control (p = 0.016 / 0.014); this is not
  followed up.

## Decisions (2026-10-22): Figure 5 symmetry layer

- **Table route.** The four-panel figure sketch is dropped.
  - Figure 5B numbers the 18 paired subtrees, with a side dot after each name,
    and puts † ‡ § beside the readout each note explains.
  - New panel 5C has one row per pair, with columns u, d_LP, d̃_LP, z_T and
    z_S. A cell is shaded when the within-pair difference is at or above the
    90th reference percentile.
  - The notes go into the Figure 5 caption.
- **Notation** (in `publication/notation.py`):
  - `natural_front_distance_sd` is `\tilde{d}_{LP}`, defined to the nearest
    front point, which in this frame is not necessarily P*.
  - `natural_travel_position` is `z_T` and `natural_state_position` is `z_S`.
  - If the colleagues choose other names, the change is one line.
- **Code.**
  - `publication/figures/terminal_symmetry.py` holds the pair numbers, note
    marks, panel C and `check_notes`. `check_notes` fails the build if a
    caption claim stops holding in the stored runs.
  - `publication/adapters/terminal.py::_symmetry` reads the compare,
    robustness, D4 and D4b tables. It hash-checks each one against its
    provenance and records them as inputs of the terminal-primary family.
  - The helpers take plain axes, so v2 can reuse them when it adopts the
    layout (its Figure 4B rows, with panel C as 4C).
- **Every caption claim comes from a stored run.** The D3 detail (the RMEV
  swap) is not cited. The claim that ABplpppapp alone accounts for ABplp is
  replaced by run-backed statements: ABplppp has that cell as its only death,
  and holding birth endpoints puts all three subtrees on their fronts.
- **Page.** Panels A–C are taller than a page, so the caption is on a second
  page (`PAGES` in `captions/terminal.py`). This lasts until panel A moves to
  Figure 1 in the next layout iteration (author, 2026-10-22).
- **Builds** of terminal-primary must run from `/home/bingran/code/embryogenesis`.
  The analysis spec stores absolute tracking paths, so a worktree gets a
  different context cache key.

## D4c. One death at a time in the posterior ABpl family (set 2026-10-22, before the run)

**Question.** The Figure 5 caption says only that holding the deaths'
endpoints at birth puts ABplp, ABplpp and ABplppp on their fronts. The claim
that one cell, ABplpppapp, carries ABplp came from a scratch check. Is that
claim right, and does the same treatment ever move a right partner off its
front? D4c stores the answer in a run, so it can be cited if a reviewer asks.

**Subtrees.** Left: ABpl, ABplp, ABplpp, ABplppp. Right: ABpr, ABprp, ABprpp,
ABprppp. ABpl/ABpr is not a kept pair; it is included as context.

**Treatments**, for every death d of a subtree, under `drop_deaths` and
`birth_endpoint`:
- `only`: d alone;
- `all_but`: every other death of the subtree (skipped when d is the only one).

The subtree is also computed with no deaths treated and with all of them
treated. These rows must reproduce D4's subtree table to 1e-12.

**Endpoint table.** Every death in the family, paired with its mirror cell
(l ↔ r at the fourth letter). For each embryo: birth frame, frames from birth
to the endpoint, the parent-division-to-endpoint distance (the travel cost),
and the birth-to-endpoint distance.

**Readouts** (descriptive, no p-values). The tie tolerance is 1e-9, as in
D4b.
- **Single cell:** ABplpppapp `only` gives d_SD ≤ 1e-9 in ABplp under both
  variants, and `all_but` ABplpppapp leaves d_SD > 1e-9 under both.
- **Right partners:** every treatment leaves every right subtree at
  d_SD ≤ 1e-9.
- **Endpoint offset:** ABplpppapp's travel exceeds its mirror's in each of
  the three embryos.

**Implementation.** `terminal_pareto/apoptosis_single_deaths.py` reuses the
D4 context and the D4b subtree front, and writes
`subtree_symmetry_v1/apoptosis_single_<date>/`. Its computation joins the
guard in `publication/tests/support.py`.

## D4c results (2026-10-22)

Run `subtree_symmetry_v1/apoptosis_single_20261022`
(`terminal_pareto/apoptosis_single_deaths.py`), replay-verified. The
no-treatment and all-deaths rows reproduce D4's subtree table (d_SD ≤ 1e-12).
All three readouts hold.

- **Single cell.** Treating ABplpppapp alone puts every left subtree on its
  front, under both removal and birth endpoint: ABpl 0.54 → 0, ABplp
  1.16 → 0, ABplpp 1.67 → 0, ABplppp 1.67 → 0. ABplppp's removal is below 12
  edges, so it is tested with the birth endpoint only. Treating all the other
  deaths together leaves every left subtree off its front: ABplp 1.16 and
  1.16; ABpl 0.92 removed and 0.54 at birth.
- **Right partners.** None of the 44 right-subtree fronts moves off its front
  (max d_SD 1.6e-15).
- **Endpoint offset.** ABplpppapp's travel exceeds its mirror's in every
  embryo: 8.9 / 5.7, 4.8 / 3.4 and 7.0 / 4.4 µm. Both tracks run a similar
  time after birth (49–69 and 63–67 frames).
- **The other death.** ABplppaaap also ends farther than its mirror in 2 of 3
  embryos, but treating it moves no subtree onto its front.

**Reading.** The lineage is symmetric. The difference sits in where one dying
cell's track ends. Because the offset recurs in every embryo, it may be real
(corpse displacement or engulfment) rather than tracking noise; the data cannot
tell. The paper should not call it an artifact.

**Reserve note.** Supplementary Note A, a one-page supplementary-style page
for reviewer questions. Code: `publication/notes/apoptosis.py`; build with
`python -m publication.notes --output-dir publication/output/notes/<name>`.
Current build: `publication/output/notes/noteA_20261022/`. Every number comes
from the hash-checked D4, D4b, D4c and compare tables, and the build fails if
a sentence stops holding. It is outside the release registry, in keeping with
the decision to keep deaths out of the supplement; promoting it is a separate
decision.

## Decisions (2026-10-27)

- The author reopened Track D. The symmetry presentation (Figure 5C and its
  notes) does not convince: the narrative stays "mirrored subtrees behave
  alike", but the readouts must support it without footnotes.
- Tweaking the current metrics and adding new ones are both in scope.
- The released Figure 5B–C and Note A stay in production until a D5
  candidate replaces them.
- **Ground rule accepted** (author, same day): readouts are chosen on
  stability, continuity and resolution across all 42 subtrees before the
  pairs are read, with random pairs as the negative control.
- **Scope:** candidates are tested on the subtree analysis first. Whether a
  readout carries over to the rest of the paper (Figures 5 and 9) is decided
  after D5b, once it is known to work there.

## D5. Readouts that hold up in tight subtrees (set 2026-10-27)

**Why.** In the stored compare run, 5 of the 9 pairs in Figure 5C have a
shaded cell (within-pair difference ≥ 90th reference percentile). D2–D4c
traced every one of them to how the readout behaves, not to the lineages:
- **Collapsed span (u, d_LP, d_CP).** These are measured in the endpoint
  frame, which divides each axis by the span between the two optima. When a
  subtree's cells are tightly bound, the span collapses and one reassignment
  sets the value. ABplpa's u is 1 against ABprpa's 0 (97th percentile) on a
  travel span of 0.04 null SD, one swap of ABplpappap with RMEV (D3).
- **Zero inflation (d_LP, d̃_LP).** Both are 0 whenever the lineage is on its
  front (26 of 42 subtrees for d̃_LP), so one cell off the front becomes a
  whole-pair difference. ABplpppapp alone gives ABplp, ABplpp and ABplppp
  d̃_LP 1.16–1.67 against 0 (90th–95th percentile; D4c).
- **Thin null (z_T, z_S, and the d̃_LP scale).** The first-cousin null can
  rest on very few cells. MSpa's rests on one three-cell group, behind the
  MSaa/MSpa z_S difference (96th percentile).

In contrast, z_T is below the 72nd percentile in every pair, and z_S in every
pair except MSaa/MSpa.

**Goal.** A pair comparison whose readouts one cell, a collapsed span or a
thin null cannot decide. The claim tested stays the same: mirrored subtrees
behave alike.

**Ground rule (set before any pair is read).** Choose a readout only on
properties that do not depend on the pair outcome, scored on all 42 subtrees:
1. **Stability:** how far one cell's removal moves the value
   (leave-one-cell-out).
2. **Continuity:** where the span is wide and the null well supported, it
   agrees with the readout it replaces, so it measures the same thing as the
   rest of the paper.
3. **Resolution:** it still separates the natural lineage from the
   first-cousin null.

Only then read the pairs, against the D1 references: size- and depth-matched
random pairs, and the between-embryo floor. Random pairs are the negative
control. A readout that also makes random pairs look alike says nothing about
mirroring, and a reviewer will check that.

**Candidates** (to narrow with the author after D5a):
- **a. Stabilized frame:** u and d_LP in null-SD units, or with the endpoint
  span floored at one null SD (the `TINY_SPAN_SD` threshold), so that a
  collapsed front cannot amplify a difference.
- **b. Single-cell influence:** each readout with its leave-one-cell-out
  range, and a flag where one cell decides it. This turns the † and § notes
  into one stated rule.
- **c. Rank readouts:** the natural lineage's percentile among first-cousin
  null draws on each axis, and the share of null lineages that dominate it.
  They are bounded, scale-free and assume no normal null.
- **d. Shared ruler:** both partners on one scale, either the pair's pooled
  null SD or the counterpart-matched cells (as in D2's matched cohort), so
  that the difference has one unit.
- **e. Equivalence framing:** the left-right difference is no larger than the
  same subtree's difference between embryos (travel), with random pairs as
  the negative control. State has one protein dataset and no floor, so its
  reference stays the random pairs.

**Scope (decided 2026-10-27).** u and d_LP are the core readouts of Figures 5
and 9. D5 tests candidates on the subtree analysis only. Carrying a readout
over to the rest of the paper would re-release every terminal family; that
is decided after D5b.

**Steps.**
- **D5a. Diagnose (no decision).** For all 42 subtrees: leave-one-cell-out
  influence on u, d_LP, d̃_LP, d_CP, z_T and z_S, alongside span and null
  support. List which 5C differences survive the removal of any single cell.
  Scratch first; a write-once run if it feeds the paper.
- **D5b. Candidates.** Compute a–e on all 42 subtrees and score them on the
  ground rule. The author picks.
- **D5c. Pairs.** Compare the pairs with the chosen readouts against random
  pairs and the embryo floor. Build a new 5C candidate; v2's 4C follows
  through `terminal_symmetry`.

**Implementation.**
- New module (proposed `terminal_pareto/symmetry_readouts.py`) writing
  `subtree_symmetry_v1/readouts_<date>/`. Write-once, with source hashes,
  replay validation and tests.
- It reads the canonical context and the stored runs. Any new front goes
  through `cross_species_analysis.analyze`.
- The hash-pinned modules above stay untouched. D0–D4c runs are not
  rewritten.
- Builds and releases run from the main checkout.

## D5a results (2026-10-27, scratch)

For each of the 42 canonical subtrees, the pooled measured-cost front was
recomputed with each cohort cell left out once: 1599 fronts through
`apoptosis_sensitivity._front`. The unmodified fronts reproduce the compare
run (max deviation 4e-16). Three removals leave no front to measure:
- ABplpappap from ABplpa: degenerate span;
- ABalpaapppa from ABalpa: degenerate span;
- MSpapaapp from MSpa: no first-cousin null left.

Scripts and tables are in the session scratchpad (`d5a/loo.py`,
`summary.py`). No run has been written.

**Readouts measured from the front are extremes; those measured from the
null are averages.**
- **d_LP and d̃_LP are on/off.** "On the front" means the natural lineage
  equals one of the swept optimal assignments, which is a discrete event.
  - 11 of the 15 subtrees off their front reach it when a single cell is
    removed; 1–3 different cells do it in each.
  - None of the 27 on-front subtrees ever leaves it.
  - ABplp, ABplpp and ABplppp each reach their front by removing any one of
    three cells: the death ABplpppapp, or the living ABplpppaaap or
    ABplppppapa. D4c's statement still holds (treating ABplpppapp alone is
    enough), but the death is not the only cell that closes the gap.
- **d_CP explodes as the span collapses.** Single-cell influence follows the
  smaller span (Spearman −0.91). ABplpa is 64, with a leave-one-out range of
  50–358.
- **u jumps.** In 12 subtrees one cell moves u by at least the typical
  difference between unrelated subtrees (0.093), mostly in MS. This does
  not follow the span (Spearman 0.2): the nearest front point jumps between
  distant parts of the front. ABplpa's u stays 1 under every removal but
  one: removing ABplpappap collapses the front to a point, where u is
  undefined.
- **z_T is stable.** The largest single-cell shift has a median of 0.24–0.37
  (by null support), at most 0.52, about the typical between-subtree
  difference (0.51).
- **z_S is stable except where the null is thin.** Median shift 0.36–0.53;
  up to 1.65 where null support is at most 11 cells (ABpra, ABprpp/ABprppp,
  ABaraa, ABarp).

**Figure 5C, comparing like with like.** Each pair difference is recomputed
after the most favourable single removal in either partner. The reference
pairs get the same treatment.
- The d_LP and d̃_LP shading of ABplp, ABplpp and ABplppp drops to the
  37th–48th percentile.
- ABplpa/ABprpa u stays at the 97th percentile. The one removal that would
  change it makes u undefined, so the gap comes from the collapse itself, not
  from a removable cell.
- MSaa/MSpa z_S stays at the 98th percentile. This is not a single-cell
  effect but D2's data problem: missing late-born cells and a one-group
  null.

**Leads for D5b.**
- z_T can carry the comparison as it is.
- z_S needs a minimum null support, or an interval.
- A continuous replacement for d_LP and d̃_LP must not depend on the natural
  lineage hitting an optimum exactly. One option is the share of the
  null-to-front distance covered, in null-SD units.
- u needs a smoother position (soft nearest point) or should drop out of
  the pair comparison.
- Rank readouts (candidate c) do not fix a thin null.

## D5b. Candidate readouts (set 2026-10-27, before the run)

All candidates live in the null-SD frame: x = T/σ_T and y = S/σ_S, with σ
from the exact first-cousin null. They use the canonical sweep, 301 points
at α_k = k/300, where α is the weight on travel, so α = 0 is the state
optimum and α = 1 the travel optimum. With d_k the distance from the natural
lineage to sweep point k:
- **α̂, trade-off position (replaces u).**
  - α̂ = Σ w_k α_k / Σ w_k, with w_k = exp(−d_k² / 2τ²) and τ = 1 null SD.
    τ is the frame's unit, the same as `TINY_SPAN_SD`, and is fixed here.
  - τ = 0.5 and τ = 2 are sensitivity checks only.
  - Hard variant α̂₀: the mean α of the sweep points nearest the natural
    lineage. For a lineage on its front, this is the middle of its trade-off
    window, the α range where the natural lineage is optimal.
  - A front collapsed to a point gives α̂ = 0.5 instead of an arbitrary end.
- **ρ, proximity (replaces d_LP and d̃_LP).**
  - ρ = 1 − d̃_LP / R_N, where R_N is the distance from the null mean to the
    nearest sweep point.
  - ρ = 1 on the front and 0 as far from it as the null mean. It still has a
    ceiling at 1.
- **G, mean trade-off gap.**
  - g_k(P) = α_k (x_P − x_k) + (1 − α_k)(y_P − y_k), and
    G = mean_k g_k(natural) / mean_k g_k(null mean).
  - G = 0 when the lineage is optimal at every trade-off, and 1 when it is
    no better than the null mean on average.
  - It is continuous across the front, but it mixes position with distance.
- **z_T and z_S** are unchanged.
- **Jackknife SE** for every readout, over the single-cell removals that
  leave a front: SE = √((m−1)/m · Σ(θ_(i) − θ̄)²).

Candidate c (rank among null shuffles) is the shuffle share computed under
resolution below. Candidates d (a shared ruler) and e (an equivalence margin)
act on pairs, so they go to D5c.

**Scores**, over all 42 subtrees. The current readouts (u, d_LP, d̃_LP,
d_CP, z_T, z_S) are scored the same way.
1. **Stability.** For each subtree, the largest rank shift one removal
   causes: the perturbed value's midrank percentile among the other 41
   subtrees' values, minus the unperturbed one. Report the median, the 90th
   percentile and the count of subtrees shifted by 0.25 or more, plus the
   median jackknife SE over the across-subtree SD.
2. **Continuity.** Spearman correlation with the replaced readout, over the
   subtrees with both spans ≥ 1 null SD and null support ≥ 12 (α̂ against
   1 − u; ρ against −d̃_LP; G against d_LP).
3. **Resolution.**
   - Distance readouts (ρ, G, z): the natural lineage against 1000
     first-cousin shuffles (`pareto_engine.compute_group_shuffle_costs`,
     seed 42). Report the share of shuffles at least as good as the natural
     lineage, and the separation (natural − shuffle mean) / shuffle SD.
   - Position readouts (α̂, u): the across-subtree SD over the median
     jackknife SE.

**Not read in D5b:** pair differences of any candidate. They come in D5c,
after the author has picked the readouts.

**Implementation (scratch).** The script repeats `analyze`'s sweep up to the
costs, so that α̂, ρ and G exist where the endpoint frame does not. Gate:
the costs reproduce `analyze`'s curve on every non-degenerate front to
1e-12.

## D5b results (2026-10-27, scratch)

1641 fronts: 42 subtrees as they are, plus each single-cell removal.
- **Gates:** the sweep reproduces `analyze` exactly; d̃_LP, z_T and z_S
  reproduce to ≤ 3e-13; the unmodified fronts match the compare run
  (4e-16).
- **Fronts lost:** α̂, ρ and G are undefined only where the null vanishes
  (MSpa without MSpapaapp). The endpoint-frame readouts also lose the two
  degenerate spans.
- **Scripts:** `d5b/candidates.py`, `score.py` and `sigma_diag.py`, in the
  session scratchpad.
- **Pairs:** no pair difference was computed.
- **Review page:** https://claude.ai/artifact/3MnYp8RyUaD1tQt92j83mk. The
  18 mirrored-pair subtrees appear under random codes (M01–M18) until the
  readouts are picked.

| Readout | Rank shift, median / 90th pct | Subtrees shifted ≥ 0.25 | Jackknife SE / SD across subtrees | Continuity (Spearman) | Resolution |
|---|---|---|---|---|---|
| u | 0.20 / 0.46 | 16 | 0.36 | — | SD/SE 2.7 |
| d_LP | 0.13 / 0.58 | 17 | 0 (27 at 0) | — | — |
| d̃_LP | 0.12 / 0.56 | 17 | 0 (27 at 0) | — | — |
| d_CP | 0.18 / 0.41 | 13 | 0.03 | — | — |
| **z_T** | **0.07 / 0.23** | **4** | 0.30 | — | separation median 2.6, min 1.2 |
| **z_S** | **0.16 / 0.27** | **5** | 0.67 | — | separation median 1.6, min −0.9 |
| α̂ (τ = 1) | 0.24 / 0.41 | 19 | 0.79 | 0.31 against 1 − u | SD/SE 1.3 |
| α̂ (τ = 0.5 / 2 / hard) | 0.22–0.23 | 16–20 | 0.74–0.90 | 0.14–0.53 | SD/SE 1.1–1.4 |
| ρ | 0.01 / 0.58 | 15 | 0 (27 at 1) | 0.97 against −d̃_LP | separation median 2.6 |
| G | 0.18 / 0.39 | 14 | 0.86 | 0.17 against d_LP | separation median 2.4, min 0.5 |

Continuity is measured on the 23 subtrees with both spans ≥ 1 null SD and
null support ≥ 12.

**Reading.**
- **No new candidate passes the ground rule better than the readouts it would
  replace:**
  - ρ behaves like d̃_LP (it is a rescaling with the same ceiling).
  - G is continuous but no more stable than u, and it fails continuity: it
    mixes position with distance.
  - α̂ is less stable than u and resolves subtrees about half as well.
    Holding the null SDs at their unperturbed values improves it only
    partly (rank shift 0.24 → 0.18, 15 subtrees ≥ 0.25, SD/SE 1.6). The front
    itself reshapes when a cell is removed.
  - Candidate c (rank among shuffles) saturates. The natural lineage beats
    all 1000 shuffles on travel in 16 subtrees and on state in 6.
- **The readouts anchored on the null (z_T, and z_S except in thin nulls)
  are the stable ones.** Every readout anchored on the front, old or new, is
  fragile at these subtree sizes (12–30 edges for most pairs), because a
  front is an optimum, and an optimum can turn on one cell.
- **Caveat for D5c: noise.** z's single-cell noise is about the size of the
  typical between-subtree difference. The jackknife SE is median 0.43–0.52
  for z_T and 0.68–0.81 for z_S; the typical differences are 0.51 and 0.79.
  A per-pair percentile against random pairs therefore has little power, in
  both directions. A difference within noise is not evidence of symmetry on
  its own.

**Proposal for D5c (for the author).**
- Compare the pairs on z_T and z_S, each with its jackknife interval. Flag z_S
  where null support is ≤ 11.
- Test across the 9 pairs together, not pair by pair: are mirror-pair
  differences smaller than random-pair differences overall? The random
  pairs are the negative control.
- Show u and d_LP (and d̃_LP) as descriptive columns without shading, or
  drop them from 5C.
- MSaa/MSpa stays a stated coverage caveat (D2), whatever readout is used.

