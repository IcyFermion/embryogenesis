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
| A. Feature calibration | 6, then the distance-metric check | [track_a_feature_calibration.md](track_a_feature_calibration.md) | nothing | done, merged 2026-10-15 (PR #4) |
| B. Figure redesign | 1, 2, 3, 4 | [track_b_figure_redesign.md](track_b_figure_redesign.md) | nothing; B3 waits for the name | `figure-ideas` |
| C. Outlier subtrees | 5 | [track_c_outlier_subtrees.md](track_c_outlier_subtrees.md) | B1 edge-change tooling | `figure-ideas` (dossier), then a research session |
| D. Subtree symmetry | Left-right counterparts (ABplpa/ABprpa, MSaa/MSpa) | [track_d_subtree_symmetry.md](track_d_subtree_symmetry.md) | nothing; shares tooling with C | `subtree-symmetry` |
| E. Layer-specific features | Per-round feature sets in the full-tree layerwise assignment | [track_e_layer_features.md](track_e_layer_features.md) | Track A conventions | `layer-features` |

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
