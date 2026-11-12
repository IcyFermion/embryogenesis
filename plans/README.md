# Plan after the collaborator discussion (October 2026)

Feedback from the discussion on the figure-ideas review page, the response
agreed with the author on 2026-10-07, and one brief per track. Each track is
meant to start in its own session from its brief.

## Feedback (summary)

1. Terminal front: show how reconstructed lineages differ from the natural
   lineage at the two optima and at maximum edge retention; look for standout
   edge changes with biological meaning.
2. Find a more readable name for "maximum edge retention".
3. Figure 3 needs simplification and a stronger visual; 3B-C may be dropped,
   3A could move to Figure 1 or the supplement; show tree distance if 3B stays.
4. Figure 4 is visually busy for outsiders (no concrete idea yet).
5. Subtree outliers (ABplpa, ABarp, ABplppp): is there a biological story?
6. Calibrate the expression geometry: how do fronts behave with random
   protein/RNA feature sets, analogous to the earlier tracking-data study?

## Tracks

| Track | Items | Brief | Depends on | Branch |
|---|---|---|---|---|
| A. Feature calibration | 6, then the distance-metric check; since 2026-11-08 the feature figure (A3; Figure 2, Supplementary Figure S2 since 2026-11-10); since 2026-11-11 a backlog item, p-values for sibling separation (A4) | [track_a_feature_calibration.md](track_a_feature_calibration.md) | nothing | done, merged 2026-10-15 (PR #4); A3 on `figures-v2` |
| B. Figure redesign | 1, 2, 3, 4; since 2026-11-08 the third feedback round on the revised layout (B5: new numbering, C. briggsae 8B-C, Figure 5 order and r, backlog, palette); since 2026-11-10 the fourth (B6: feature figure to S2, main figures 1-7 again, d_CP beside d_LP in Figure 4, C. briggsae front alone with its comparisons in the backlog); since 2026-11-12 the layout is official (B7: the manuscript set at `publication/output/manuscript/current`, verified and promoted; text in the authors' Overleaf) | [track_b_figure_redesign.md](track_b_figure_redesign.md) | nothing; B3 waits for the name | `figure-ideas`; B5-B7 on `figures-v2` (PR to `main`) |
| C. Outlier subtrees | 5 | [track_c_outlier_subtrees.md](track_c_outlier_subtrees.md) | B1 edge-change tooling | `figure-ideas` (dossier), then a research session |
| D. Subtree symmetry | Left-right counterparts (ABplpa/ABprpa, MSaa/MSpa); since 2026-10-19 also apoptotic cells (D4); since 2026-10-27 readouts robust in tight subtrees (D5); since 2026-11-02 the 24-cell subtree floor, with symmetry dropped from the paper (D6); since 2026-11-08 the terminal cohort check (D7: cells that end embryogenesis undivided) | [track_d_subtree_symmetry.md](track_d_subtree_symmetry.md) | nothing; shares tooling with C | merged 2026-11-04 (PR #7: 24-cell floor, mirrored-pair layer retired; D0–D4c were PR #6); D7 decided 2026-11-09 (275-edge cohort kept; 321-edge check in the backlog) |
| E. Per-layer fronts (was: layer-specific features) | Per-layer travel and cell-state advantage in the full tree; the mid-embryo dip | [track_e_layer_features.md](track_e_layer_features.md) | Track A conventions | done 2026-11-05 (claims for layers 3-8; main-text panels become revised Figure 6B-C; S5 supplement; timeline H in reserve) |

Order: start A first, because it can change how much weight the geometric
metrics carry in every figure. Run B in parallel. Start C after B1.

Added 2026-10-15 (author): D and E. Start D first; it reuses the terminal
subtree caches and may feed C. E is larger. Its E1 decision gate decides
whether per-round features go any further.

## Decisions so far (2026-10-07)

- Overall plan approved; iterate track by track.
- Replacement name for "maximum edge retention": on hold while the author
  discusses candidates with colleagues ("edge agreement", "lineage agreement",
  "edge similarity").
- Track A: top-k selected features for k = 5, 10, 20, 30, and random draws of
  the same sizes from (i) the full panel and (ii) highly variable features
  only. Decide what to present after seeing the results.
- Distance metric (Track A2, decided 2026-10-14): cosine molecular distance
  everywhere. The full tree and the partial forest were rerun with cosine
  under new run identities; Figures 7, 10, 11 and S4 need re-release.
- Track E, E0 (decided 2026-10-15): the label is the Schema D mixture of a
  cell's terminal descendants in the complete lineage (programmed death kept
  as a class); rounds 5-8 are pooled; the criterion is equal-relatedness
  separation (Spearman of cosine against Bray-Curtis mixture distance within
  path-length groups, each path length weighted equally), with sisters-only
  beside it; lineage memory is measured on near-identical same-majority pairs
  (Bray-Curtis at most 0.2) against majority-class shuffles. Mixed generations
  within a round are a known caveat. Run
  `layer_features_v1/e0_groupweight_20261015`; write-up in
  `full_tree_pareto/LAYER_FEATURES.md`.
- Track E, E1 (settings decided 2026-10-15 before results): k = 20 everywhere,
  univariate ranking on the E0 criterion, held-out gate (beat the top-20 in
  at least 80% of halves and at least the 90th random percentile, in each of
  strata 2-4). Run `layer_features_v1/e1_20261015`: the gate fails (only
  strata 1-2 pass), so the global set stands and E2 is not run. A single
  criterion-ranked set shared by all rounds matches the per-round sets and
  beats the top-20 on held-out halves; whether that should reopen Track A is
  for the author.
- Track E, redirected (2026-10-19, author): per-layer features are closed (E2
  not run). The per-layer fronts show the natural lineage's travel advantage
  over cousin shuffles is smallest in layers 6-7 (layer = 9 − round; cells
  born at about 70-190 cells), in both species, while the cell-state advantage
  stays flat. Remaining work: a bootstrap over first-cousin groups with a
  dip contrast and decision rule fixed before the intervals (the layers were
  chosen after seeing the point estimates, so a pass shows only that the dip
  exceeds resampling noise) (E3), a write-up and one
  candidate figure with claims limited to "consistent with" fate-based cell
  sorting and "coincides with" the expression hourglass waist (E4), then
  close (E5). Layers are numbered top-down in text and figures; code keeps
  rounds.

## Running tracks in parallel

Use one git worktree per concurrent session; two sessions switching branches in
one checkout will collide. `data/`, `expression_embedding/results/` and every
`output/` folder are Git-ignored, so a fresh worktree has none of them. Link
them from the main checkout (new analyses write new run folders, so sharing is
safe; never let two sessions write the same run folder):

```bash
git worktree add ../embryogenesis-subtree-symmetry -b subtree-symmetry main
cd ../embryogenesis-subtree-symmetry
for p in expression_embedding/results terminal_pareto/output full_tree_pareto/output publication/output; do
  ln -s "$(realpath ../embryogenesis/$p)" "$p"
done
# data/ already exists (it has tracked files), so link its untracked entries one by one
for f in ../embryogenesis/data/*; do
  [ -e "data/$(basename "$f")" ] || ln -s "$(realpath "$f")" "data/$(basename "$f")"
done
# The tracking files must be real copies: the cache checks reject paths that resolve outside the worktree
for e in embryo1 embryo2 embryo3; do rm "data/$e" && cp -r "../embryogenesis/data/$e" "data/$e"; done
```

Because the output folders are links, every run a worktree writes lands in the
main checkout. Untracked work you want to keep (e.g. `artifacts/`) must be
moved out before `git worktree remove --force`. When the branch is merged,
remove the worktree from the main checkout, then `git worktree prune` and
`git branch -d BRANCH`. In a worktree the terminal context cache key differs
(absolute paths), so `publication.tests.test_terminal_primary` passes only in
the main checkout.

Never edit hash-pinned analysis modules (listed in each brief); add new modules
with their own caches, following `full_tree_pareto/cousin_references.py`.
