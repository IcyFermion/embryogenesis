# Terminal publication handoff

Updated 2026-09-23. The author approved switching to pooled travel and
endpoint-normalized fronts. The current publication is `output/publication/`;
its release manifest identifies the working run and records artifact/cache
hashes. [README.md](README.md) is the analysis and reproduction reference.

## Current state

- Primary: `pooled_tracking_v1`, 275 terminal edges, 42 eligible subtrees.
- Working run: `output/runs/pooled_tracking_v1/migration_candidate_20260920/`.
  The historical run ID is retained; there is no pending figure-promotion gate.
- Figure 1 amendment is standalone; Figure 2 is a single panel; Figure 3A is
  included; Figure 5 has two panels; S3 is one 2-by-2 page with insets.
- `figS3_cross_geometry.py` is the current S3 renderer. The prototype name and
  superseded plotting entrypoints have been retired.
- Standalone wrappers now use the intended figure numbers, including S1–S3;
  the independent schematic says “Figure 1 amendment”.
- Previous embryo-1 figures are frozen in `output/legacy/embryo1/publication/`.
  Historical no-profile commands write separately to `output/legacy/rebuild/`.
- The working-run analysis caches and saved assignment sets are unchanged by
  promotion. Rebuilds are staged; publication is explicit and preserves the
  previous release. Use `publication_release.py --verify` to check the snapshot.

## Validation and review

The prior 29-check migration suite passed, including exact legacy reproduction,
pooled reconstruction, profile/cache identity, endpoint inverse checks, wrapper
organization and replay of all 4,816 S3 projections (maximum error
\(4.55\times10^{-13}\)). The release promotion reruns this suite and requires
all nine compiled figure wrappers to be one page. The current report is in
`output/runs/pooled_tracking_v1/migration_candidate_20260920/validation/`.
After promotion, historical publication checks target the frozen legacy
bundle; they must not compare old hashes against the newly published figures.

Promotion completed with **29/29 migration checks, 13/13 profile/cache tests,
and 3/3 release-boundary tests passing**. The published manifest verifies
57 figure/data/provenance files and all recorded analysis files. The frozen
old publication matches its pre-organization hashes.

All figure wrappers were visually inspected during the organization review.
The publication pass also checks the new standalone numbering and S3 render.
The subsequent Figure 5 readability revision replaces shape-coded scatter
panels with a named five-subtree natural/null comparison and a grouped row
display of all 42 subtrees. Canonical position and natural distance have
separate columns with common scales; every subtree name is legible. The
compiled page was visually inspected. An artist-level audit confirms all ten
major-subtree distances and all 84 position/distance values match the cached
table exactly, with each of the 42 names present once in panel B. Figure S1
and the analysis caches are unchanged. The prior Figure 5 is preserved in
the working run's `validation/fig5_before_redesign_20260923/`, and the release
workflow archives the complete previous publication. Promotion tests cover
preserved legacy files, hash verification, and rollback on failure.

## Next manuscript iteration

1. Integrate the standalone amendment with the co-author's Figure 1; its final
   panel position is still the co-author's decision. Do not edit their artwork
   implicitly as part of rebuilding the terminal figures.
2. Update Methods and Results to the pooled definition and 275-edge/42-subtree
   scope. Use the published filenames and endpoint/canonical definitions in
   the README. Remove obsolete references to the schematic as Figure 2A/5A.
3. Frame S3 around shared global shape, tracking-dependent local optima, and
   the pooled compromise across A–C. Retain the limited holdout result in the
   sensitivity discussion. Avoid claims of isolated tracking error or a
   uniquely correct biological geometry.
4. Keep full-tree figures and the paused cross-species pilot separate. No new
   optimizations or pilot extensions are required for manuscript integration.

## Short development record

- Tracking sensitivity: similar front shapes coexist with variable parent
  identities and conserved neighborhoods. Mixed-weight biological-parent
  agreement was about 0.76–0.79; pure-travel transfer penalties were much larger
  than typical mixed-objective penalties. These are three-geometry findings.
- Pooled versus minimax: at no increase in state cost, the certified pooled
  budget solution reduced mean travel by 2.984%, but embryo 3 increased by
  0.127%. The minimax solution reduced mean travel by 2.334%, with reductions
  in all three embryos (2.409%, 2.316%, 2.276%). These budget solutions are
  separate from the 301-weight sweeps shown in S3.
- Pooled migration: explicit compatibility profiles separated the 299-edge
  legacy cohort from the matched 275-edge cohort. Endpoint normalization
  remained display-only. Natural lineage appeared on the sampled front in
  16/43 legacy, 15/42 matched embryo-1, and 27/42 pooled subtrees. This
  descriptive comparison was not the criterion for adopting pooled travel.
- Presentation: schematic moved from Figure 5 to Figure 2, then to a standalone
  Figure 1 amendment. Three-panel and two-panel S3 drafts were superseded by
  the four-geometry page with insets. September 23: author-approved promotion
  and consolidation into two maintained documents.
- Paused CE protein/RNA and CB AF16 pilot: 188 shared edges, raw-3D and XY
  comparisons, two AF16 embryos. Axial calibration and molecular provenance
  remain unresolved. Coverage was matched/available cells, not optimized
  retention. No species ranking is adopted.

## Legacy and intermediate material

No intermediate layout is a live publication alternative. The old loose
prototype directories/links, standalone amendment draft, prototype scripts,
exploratory consensus/metric-comparison scripts, and paused pilot scripts have
been cleared from the active directory after archive verification.

- `legacy/development_sources_20260923.tar.gz`: exact pre-consolidation Python,
  TeX, and seven Markdown records at their original repository-relative paths.
  Its adjacent JSON provides member hashes and the archive hash. This retains
  solver settings, pilot resume/calibration notes, detailed sensitivity results,
  and old reproduction entrypoints without maintaining parallel documents.
- `output/legacy/development_outputs_20260923.tar.gz`: verified output checkpoint
  of prototypes, prior housekeeping archive, paused pilots, consensus results,
  tracking sensitivity and pooled/minimax comparison. Adjacent JSON records
  hashes. Active S3/validation input directories remain available separately.
- `output/legacy/pre_promotion_run_20260923.tar.gz`: working run before the
  publication cleanup. Existing baseline archives remain the regression source.
- `output/legacy/embryo1/publication/`: old published figures, directly browsable.
- `output/legacy/diagnostics/`: 26 older diagnostic plots gathered from the
  top-level output and CE/CB plotting folders, with a relocation/hash manifest.
- Later promotions preserve replaced bundles under `output/legacy/releases/`.

To inspect an archived document or script, use `tar -xOf ARCHIVE MEMBER` with
its repository-relative member name. To resume a pilot, first extract the
source/output archives into a separate temporary checkout and read its original
record; do not unpack over the active module or silently rerun into the release.
Archives of generated output are local and Git-ignored, so they still require
external backup. Consolidation is not permission to discard scientific inputs.
