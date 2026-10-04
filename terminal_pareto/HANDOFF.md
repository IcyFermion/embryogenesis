# Terminal publication handoff

> **Layout change, 2026-10-04.** Figures, captions and releases moved to
> [`publication/`](../publication/README.md); production is `publication/output/production/`.
> This package is numerical only. Figure material that this record places in
> run folders (`figures/`, `publication/`, `layout_history/`) or in the old
> `output/publication/` is archived in
> `publication/output/archive/migrated_20261004/` and removed after the
> refactor PR. Commands named below that drew or released figures are retired;
> see the package README for current commands.

Updated 2026-10-01. The author approved switching to pooled travel and
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
- The additional cross-species figures are author-numbered Figures 8 and 9;
  their captioned PDFs, TeX and panel PDF/PNG assets are now in
  `output/publication/`. Their source builds remain in the separate comparison
  run described below, and a `cross_species` section of the release manifest
  records this distinct 187-edge analysis rather than relabeling the 275-edge run.
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
57 original figure/data/provenance files and all recorded analysis files. The
2026-10-01 Figure 8/9 addition extends verification to 66 assets. The frozen
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
4. Keep full-tree figures and the cross-species review figures separate. No new
   optimizations or pilot extensions are required for manuscript integration.

## Resumed cross-species comparison (2026-09-30)

The author requested a CE protein / CE RNA / CB RNA terminal-only figure in
the current Figure 2 style. The new implementation is
`cross_species_analysis.py`, with `fig_terminal_cross_species.py` for layout
and `test_cross_species_terminal.py` for scientific checks. The original
188-edge pilot remains in the archives below.

Current run: `output/runs/cross_species_terminal_v1/pooled_comparison_20260930/`.
The matched set is 187 edges; the original pilot's glial edge
`ABalppaapp -> ABalppaappp` is absent from the required CE tracking intersection.
CE protein and RNA share published pooled travel, restricted from the 275-edge
reference with its fixed denominators. CB pools the two AF16 embryos over a
196-edge shared tracking/RNA reference, using fixed natural-total denominators.
See the README for commands, generated-file stems and measurement conventions.

| Configuration | Existing-3D d_LP | Existing-3D max retention | XY d_LP | XY max retention |
|---|---:|---:|---:|---:|
| CE protein | 0.0099 | 0.8717 | 0.0196 | 0.7594 |
| CE RNA | 0.0107 | 0.9198 | 0.0494 | 0.6471 |
| CB AF16 RNA | 0.0311 | 0.7380 | 0.0411 | 0.6524 |

Six focused tests passed. All 13,828 saved base/dense assignments passed
permutation, raw travel/state cost, biological-parent retention, and endpoint
coordinate replay; 20 input/code hashes were checked. Base sweeps have 301
weights; the six pooled runs also have nested 1,201-weight checks. Maximum
retention is unchanged, distance changes are below 5e-16, and the largest
absolute change in canonical position is below 1.2e-7. This checks sampled
grid stability and does not enumerate all nondominated/tied optima.

The author's layout revision combines existing 3D above 2D (XY) in one
six-panel comparison (`terminal_cross_species_comparison`), with CE protein,
CE RNA and CB RNA in columns. A separate paired overlay figure
(`terminal_cross_species_overlay`) shows both geometries with identical main
axes and near-natural zoom limits. Species colors, line styles and landmarks
are consistent across overlays. Individual row and overlay exports, plus the
tracking transfer companion, remain available. All are standalone PNG/PDF
artifacts with captions and render provenance; the numerical checkpoint is
unchanged by this revision. At the author's request, the C. briggsae 3D caveat
is caption-only and describes limitations of traditional embryo-tracking
techniques, particularly for z-axis measurements, rather than an unverified
axial scale. Exact wording remains provisional pending experimental
collaborator input. Stage alignment and expression provenance still need
review; XY sensitivity changes the proximity ordering of CE RNA versus CB RNA.
No species ranking or controlled protein/RNA modality claim is adopted.
The initial comparison build included no promotion or 1,001-weight migration;
the subsequent approved Figure 8/9 release is recorded below.

The author subsequently assigned **Figure 8** to the six-panel comparison
and **Figure 9** to the paired overlays and requested captioned publication
pages. Figure 8 now omits natural-to-maximum-retention connectors. Figure 9
adds panel C below its 3D/XY overlays, following the Figure 5B named-row style.
The three configurations are grouped under each geometry, with separate
columns for u, d_LP and d_NP. Scales are common across the six rows within
each metric. Hollow diamonds and dashed natural-lineage connectors identify
the closest attained point P* in Figure 9A--B, not the maximum-retention point.
Both distances use this same P*, and u is normalized arc length from the
travel optimum to P*. All 18 displayed values are replay-checked against the
unchanged numerical cache, including analytic first-cousin-mean distances.
`cross_species_publication.py` assembles these from hash-checked panels
using the earlier terminal wrapper style. Both PDFs compile to one page,
carry the requested figure labels, and have been visually inspected with no
layout defects or TeX warnings. Their stems under the run's `publication/`
are `fig8_terminal_cross_species_comparison` and
`fig9_terminal_cross_species_overlays`; editable TeX and panel PDF/PNG assets
accompany them. The tracking caveat is caption-only, shared across captions,
and still explicitly provisional. Thirteen tests pass (six scientific and
seven publication/numbering/archive/artist checks), including connector
endpoints, the absence of Figure 8 connectors, and all 18 metric-panel points.
The build rechecked all 13,828
assignments and 20 source hashes; its manifest records eight assembled
artifacts and nine unchanged analysis files. Rebuilds archive the previous
assembled layout with hashes. Initial numbering did not promote the comparison
or revise historical full-tree Figure 8 assets. Subsequent production promotion
is explicit, and the manuscript remains untouched.

Presentation revision, 2026-10-01: Figure 8 now uses a shared heading style
with the partial-forest comparison. Each of its six panels displays the
species/molecular configuration, followed immediately by a bold `3D tracking`
or `2D (XY) tracking` label. Embryo counts no longer clutter titles; replicate
counts and pooling definitions remain in captions. The faint row-wide geometry
labels are removed and configuration names repeat in the XY row. The preceding
panel and captioned layouts are preserved under `layout_history/`; numerical
caches and all canonical metrics remain unchanged.

Production addition, 2026-10-01: the author requested Figures 8/9 alongside
the other terminal figures in `output/publication/`. Their eight assembled
assets were copied byte-for-byte without rebuilding, with a frozen assembly
manifest named `cross_species_publication_manifest.json`. The working run is
retained for provenance rather than moved away from its numerical caches.
`publication_release.py --cross-species` performs this additive, staged,
rollback-tested release. It extends the main release manifest with a distinct
comparison section, nine comparison analysis hashes and figure numbers.
All 57 preceding production assets and nine comparison numerical files match
the pre-release hashes. The entire former 58-file production directory is
recoverable in the verified archive
`output/legacy/releases/20261001T154243342749Z/publication/`.

Both pages were rendered and visually checked again. The release validates
all 13,828 assignments and 20 source hashes, and `publication_release.py
--verify` now checks 66 assets plus both runs' recorded analysis files. All
20 focused release/cross-species tests pass, including failed-promotion rollback,
source-run preservation and retention of Figures 8/9 through future pooled
layout promotion. The 13 profile/cache regression checks also pass. No
full-tree figures, numerical optimizations or manuscript files were changed.

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
