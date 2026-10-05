# Shared publication pipeline: implementation plan

Prepared 2026-10-02 against commit `6b464be`. Status: **plan only**.
No analysis, figure, cache, or production changes are authorized by this document
alone. The user requested a detailed plan for a subsequent implementation agent.

## Objective and boundary

The author's architectural requirement is **front-end/back-end separation**:

- **Back ends:** `terminal_pareto/` and `full_tree_pareto/` conduct experiments
  and expose validated scientific results.
- **Front end:** `publication/` consumes those results to produce figures,
  tables and TeX documents. It must not run experiments implicitly.
- **Interface:** explicit result contracts and provenance, with adapters for
  existing caches. A presentation change should not require editing both back
  ends; a scientific change produces a new validated result for the front end.

This is a separation of responsibilities, not merely a shared utility folder.
Active layouts and document production ultimately belong to the front end;
old back-end production entrypoints survive only as compatibility delegates.

Create a repository-level `publication/` package alongside `terminal_pareto/`
and `full_tree_pareto/`. It owns active publication presentation: notation,
visual conventions, reusable artists, figure layouts, table formatting,
caption templates, document assembly, rendering provenance, and release
coordination. The two scientific packages continue to own datasets, cohorts,
objectives, optimization, null/reference generation, scientific metrics,
numerical caches, and scientific validation.

The primary user story is: edit a canonical metric's displayed notation once,
then rebuild all affected figures, captions and tables from existing validated
caches. This must not rename persisted scientific fields or rerun optimization.
Changing the mathematical definition of a metric is a separate scientific change.

The first migration preserves current scientific results and appearance.
Do not combine it with a redesign, a 1,001-weight migration, new reference models,
manuscript edits, Figure 1 co-author integration, or pooled Figure 7 promotion.
Move production code, not numerical datasets or existing output directories.

## Authoritative scope and observed checkpoint

Read `PUBLICATION_FIGURES_HANDOFF.md`, both package READMEs,
`terminal_pareto/HANDOFF.md`, `full_tree_pareto/POOLED_HANDOFF.md`, and
`full_tree_pareto/CROSS_SPECIES_HANDOFF.md` in that order before implementation.
The repository-level handoff explicitly marks its older full-tree material as
historical; do not use those sections as the active specification.

| Family | Scientific input | Presentation/release status |
|---|---|---|
| Terminal primary | 275 edges, 42 subtrees, pooled travel, cosine molecular distance | Figures 1 amendment, 2–6, S1–S3 and Table 1 published |
| Terminal cross-species | 187 matched terminal edges; six configuration/geometry cases | Figures 8–9 published |
| Pooled full tree | 978 measured nodes, 974 edges; Euclidean molecular distance | Revised Figure 7 and provisional S4 approved in working run, **not promoted** |
| Partial-forest cross-species | 485 nodes, 454 edges, 31 boundary roots | Figures 10–11 published |
| Historical full tree | Former single-embryo Figure 7 and supplement | Still present in production beside Figures 10–11 |

Read-only inspection found a clean tracked checkout, terminal production with
67 files (66 assets recorded by its main manifest plus that manifest), full-tree
production with 25 files, and the pooled Figure 7 working publication with 15
files. These are inventory observations, not a fresh scientific validation.
The implementation must run baseline verifiers before changing code.

Preserve these run identities and paths:

- `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920/`
- `terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930/`
- `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_20260927/`
- `full_tree_pareto/output/runs/cross_species_layerwise_v1/terminal_anchored_20260930/`

## Findings that determine the design

1. Sharing already exists, but ownership is inverted. Full-tree renderers import
   `terminal_pareto.plot_style`, `EndpointTransform`, and artists from terminal
   figure scripts. `full_tree_pareto/fig_cross_species.py` imports
   `fig_terminal_cross_species` and Figure 2 helpers; full-tree cross-species
   assembly imports the private terminal wrapper helper `_write`.
2. Symbols are repeated across `fig1_endpoint_amendment.py`, Figure 5 rendering,
   the Table 1 writer, both cross-species captions, and cross-species metric
   artists. Palette sharing alone does not solve this.
3. `fig5_table1_ce_canonical_metrics.py` mixes scientific calculations, cache
   validation, table formatting and a CLI. Move its table formatting separately;
   do not relocate the whole file as if it were only presentation.
4. `full_tree_pareto/publication_build.py` mixes rendering, captions, compilation,
   analysis orchestration and promotion. `resume_paired` also imports its runtime
   and sweep defaults. Preserve these imports and values through compatibility
   exports; moving the renderer must not break computation entrypoints.
5. `front_coordinates.py` is scientific-cache hash-pinned by both cross-species
   analyses. `pooled_analysis.py` is pinned to expensive full-tree caches.
   Moving them or replacing them with shims changes their source hashes.
6. Compilation, layout archives, artifact hashing, staging and rollback are
   implemented repeatedly, with different manifests and release policies.
   Current full-tree Figure 7 promotion refuses to replace Figures 10–11.
   That guard must remain until a tested mixed-release replacement exists.
7. The legacy heuristic inventory describes historical scope. It cannot become
   the pooled supplement's inventory just because both are tables.

## Proposed package and dependency direction

Use a small package, not a general plotting framework. Create files only when
the corresponding responsibility is extracted:

```text
publication/
  __init__.py
  __main__.py              # cache-only presentation CLI
  README.md               # commands, ownership, notation-change recipe
  notation.py             # stable concept IDs, symbols, labels, descriptions
  style.py                # palette, typography, export conventions
  registry.py             # figure families, numbering, stems, dependencies
  contracts.py            # minimal validated presentation inputs
  artists.py              # fronts, markers, insets, comparison headings
  canonical_panels.py     # named metric rows and closest-point overlays
  tables.py               # Table 1 and inventory formatting
  assembly.py             # wrapper preamble, compilation, page/label checks
  provenance.py           # presentation dependency hashes and manifests
  release.py              # shared archive/stage/verify/rollback mechanics
  adapters/
    terminal.py           # existing terminal cache loaders and validators
    full_tree.py          # existing full-tree cache loaders and validators
  figures/
    terminal.py           # terminal-specific layouts (split if unwieldy)
    full_tree.py          # pooled Figure 7 and supplement
    cross_species.py      # common layout, explicit scope-specific inputs
  captions/
    terminal.py
    full_tree.py
    cross_species.py
  tests/
```

Dependency direction: scientific packages supply validated inputs to adapters;
adapters supply presentation records to layouts; layouts consume shared artists,
notation, and assembly. Shared artists must not import either scientific package
or import a figure-specific script. Adapters may import scientific code.
Old entrypoints delegate inward; avoid cycles through those compatibility shims.

An input contract needs only the fields actually used: analysis/run identity,
scope/cohort, ordered attained front points and assignment identifiers, anchors,
natural/null coordinates, canonical metrics, reference-family definitions,
validated table rows, and explicit caption metadata. Do not require every figure
to fit one universal dataframe. Scientific validators remain authoritative.

Keep terminal/full-tree scientific distinctions explicit in adapters and captions:
cosine versus Euclidean costs; 275/187/974/454-edge cohorts; per-round versus
aggregate anchors; analytic means versus sampled means; differing references;
partial forest versus complete represented forest. Never infer scope from a
filename, figure number, or number of dataframe rows.

## Single-source notation

Key concepts by semantic identity, not by current printed symbol. For example,
`natural_front_distance` maps to persisted `d_lp` or `d_LP` in the appropriate
adapter and currently displays the natural-to-front distance symbol. Similarly
map `canonical_position` to existing `u_L`/`u_lineage_lp` fields without renaming
them. Record definition, symbol body, readable label, and formatting precision
in one place; keep context-dependent scientific qualifications in captions.

Provide explicit Matplotlib mathtext, LaTeX and plain-text formatting functions.
Store one mathematical symbol body and add the appropriate delimiters per
backend; do not do broad string replacements on complete TeX documents.
Caption templates and table headings must request symbols from this registry.
Generated standalone TeX must remain self-contained (inline macro definitions
are acceptable); no dependency on the repository at compile time.

Cover endpoint axes, natural lineage, null mean, closest attained point,
canonical position, both canonical distances, retention and the distinct
cousin-relative statistic. Preserve existing definitions: both distances use
the same closest attained point, not an interpolated or maximum-retention point.
Narrative Markdown and the manuscript remain manually maintained; document that
boundary rather than claiming a notation change edits all prose in the repo.

## Source-to-destination migration map

| Current source | Action |
|---|---|
| `terminal_pareto/plot_style.py` | Extract implementation to `publication/style.py`; keep public compatibility exports after hash audit |
| `terminal_pareto/front_coordinates.py` | **Keep unchanged initially**; adapters reuse its transform; do not copy a second implementation |
| `terminal_pareto/fig2_fig3_ce_terminal_pareto.py` | Extract shared colormap, limits and front artists, then terminal layouts; retain CLI |
| Terminal `fig1_*`, `fig4_*`, `fig5_figs1_*`, `fig6*`, `figS3_cross_geometry.py` | Move presentation incrementally; keep analysis helpers/loaders in original scientific scope |
| `terminal_pareto/fig5_table1_ce_canonical_metrics.py` | Retain metric computation/cache API; extract `write_stats_table` formatting and delegate if hash-safe |
| `terminal_pareto/fig_terminal_cross_species.py` | Extract artists, headings, canonical rows and layout; retain compatible `render` entrypoint |
| `full_tree_pareto/fig_cross_species.py` | Consume common cross-species presentation without importing terminal figure scripts |
| `full_tree_pareto/publication_build.py` | Extract plots, wrapper templates, inventory formatting, compilation and release mechanics; preserve scientific build API/defaults |
| Both `cross_species_publication.py` files and terminal `publication_wrappers.py` | Move captions/assembly to shared package; retain command/API wrappers |
| Terminal `publication_build.py` / `publication_release.py` | Delegate presentation and release machinery; preserve explicit profile/run and legacy commands |
| Legacy `fig7_fig8_*`, `heuristic_inventory.py`, historical TeX wrappers | Keep historical behavior and output targets; do not silently relabel as current pooled publication |
| `assets/fig3A_ce_null_models.tex` | Audit notation; shared presentation asset may move with a stable compatibility path |

Before editing any listed file, inspect actual cached source manifests to learn
whether it is pinned. The table is a responsibility map, not permission to
invalidate caches. Retain a documented compatibility island where necessary.

## Cache and provenance policy

- Separate **scientific identity**, **presentation identity**, and **release
  identity**. Notation, colors, captions and layout belong only to presentation.
- Never edit old scientific manifests, disable source checks, or re-hash changed
  source into an existing scientific identity to make reuse pass.
- Keep currently pinned analysis files byte-identical during this extraction.
  If a desired extraction cannot meet that rule, leave that dependency in place
  and record a later versioned scientific migration; it is not a blocker to
  centralizing presentation.
- New presentation manifests include the full used shared dependency set,
  including notation, style, captions, templates, artists and adapter source.
  A stale shared dependency must invalidate affected rendered assets.
- Preserve historical manifests. Verification of a frozen release must remain
  possible after presentation source changes; distinguish artifact-integrity
  verification from readiness to assemble a new build with current source.
  Resolve that distinction explicitly where current verifiers check live source.
- Default new CLI builds are cache-only and require an explicit candidate output
  directory. Missing/stale numerical inputs fail with an actionable command;
  rendering must not fall through to solvers or null generation.
- Keep numerical run paths and scientific filenames unchanged. Permit a separate
  presentation destination even when the old CLI writes beside its input run.

## Phased implementation and checkpoints

### 0. Baseline and preservation

Record commit/status, source hashes, production inventories, scientific manifests,
and the four run inventories. Run current release verifiers. Save hash manifests
and recoverable copies of any ignored artifact directory before it is touched.
Inspect existing validation reports and reuse numerical caches. Missing cache
inputs should be reported; do not reconstruct them as part of cleanup.

### 1. Shared notation/style and one vertical slice

Create notation/style/assembly basics. Migrate Figures 9 and 11 together because
they already share overlay and canonical-row logic. Preserve values, colors,
titles, caveats, sizes, ordering, numbering, and output stems. Render to isolated
candidate directories. This proves terminal/full-tree independence and the
single-source notation path before moving the larger pipeline.

Checkpoint: both families build, numerical hashes are unchanged, current symbols
match the baseline, and a temporary notation override updates both families.

### 2. Complete active presentation extraction

Migrate Figures 8/10, then terminal Figures 1–6/S1–S3 and Table 1, then pooled
Figure 7/S4 and its inventory. Migrate small coherent groups and validate each
before continuing. Do not generate historical full-tree panels to stand in for
the approved pooled candidate. Keep compatibility CLIs/imports working.

Checkpoint: all active figures/tables obtain notation and shared conventions
from `publication/`; no full-tree presentation imports terminal figure scripts.
Document intentional legacy and hash-pinned exceptions.

### 3. Unified build registry and release coordination

Add stable family IDs, wrapper stems, numbering, input profiles, status, and
owned asset lists. Treat Figure 1 amendment and provisional S4 explicitly.
Expose proposed commands such as:

```text
python -m publication list
python -m publication build --family terminal-cross-species --output-dir PATH
python -m publication build --family full-tree-pooled --output-dir PATH
python -m publication verify --build PATH
```

Allow a build manifest to describe all families while retaining separate current
production roots. Unify low-level transactional mechanics, but preserve each
family's scientific validation and manifest semantics. Make release ownership
explicit: additions/replacements may alter only the selected family's assets.

For future Figure 7 promotion, stage a complete mixed full-tree release that
preserves Figures 10/11, updates all affected preserved-file records coherently,
verifies the staged result, archives the old release, and rolls back on failure.
Test this with temporary fixtures. Keep the current production guard until this
path passes. Implementing this capability does **not** promote the pooled build.

### 4. Documentation and implementation handoff

Update source maps, new/old command mappings, notation-change recipe, artifact
ownership and cache/release boundaries. Add a short link from the repository
handoff to the new package README. Preserve scientific conventions in the two
analysis READMEs. Record completed phases, candidate paths, checks, visual review
and any compatibility exceptions. Deliver a reviewable diff and candidates;
publication and manuscript integration remain separate actions.

## Validation and acceptance

Baseline commands from repository root, in `dev` (existing commands, not new APIs):

```bash
python terminal_pareto/publication_release.py --verify
python -m full_tree_pareto.cross_species_publication --verify-production
python -m unittest terminal_pareto.test_publication_release terminal_pareto.test_cross_species_publication
python -m unittest discover -s full_tree_pareto -p 'test_*.py'
```

Use `/home/bingran/miniconda3/envs/dev/bin/python` when Conda activation is
unavailable; set `OPENBLAS_NUM_THREADS=1` and a writable `MPLCONFIGDIR` for render
runs. Inspect each test's scope before running expensive checks. Reuse completed
scientific replay reports where appropriate, but do run affected replay/cache
validation and the final migration check (`validate_pooled_migration.py`).

Required evidence:

1. Pre/post numerical file hashes, identities, cohort order, assignments, anchors,
   metric values and null draws agree. Tests make solver/null-generation calls
   fail if reached from the new cache-only build path.
2. A temporary notation change propagates to every affected active artist,
   caption and table, including Figures 1/5/9/11 and Table 1. Stable CSV keys
   and numerical hashes remain unchanged. Restore the intended notation.
3. Existing artist-level checks still verify all canonical rows/values and
   connector endpoints. Preserve closest-attained versus maximum-retention
   distinctions and recorded front orientation/tie handling.
4. Compile wrappers; check expected pages, figure labels, missing assets and
   relevant TeX warnings. Compare rasterized pages/panels with baselines and
   inspect differences, not PDF byte equality (metadata can change). Inspect
   final-width labels, legends, insets and table fit visually.
5. Shared-dependency changes invalidate presentation manifests. Altered numerical
   caches still fail validation. Frozen production artifact verification works
   without silently blessing changed live source.
6. Test additive retention, family asset collisions, staged verification failure,
   archive integrity and rollback. Figure 7 replacement must retain Figures
   10/11 and both verification paths must agree on the resulting manifests.
7. Old documented commands/imports remain usable with the same defaults and
   output boundaries. Historical entrypoints do not write to current production.
8. Current production and frozen legacy archives remain byte-identical throughout
   the implementation. No manuscript or publication promotion is part of done.

No test suite or live scientific verifier was run while preparing this plan.
The planning pass read source, documentation and live manifest inventories only.

## Efficient execution instructions for the next agent

Work one phase at a time and keep a compact progress record in this document.
Start with the cited entrypoints and migration map instead of scanning the whole
repository. Do not read or regenerate large arrays unless a validation needs
them. Reuse existing tests and cached fixtures; prioritize the notation vertical
slice, provenance boundaries and release regressions over exhaustive rewrites.
Do not spawn additional agents unless the user requests them.

Suggested task prompt:

> Implement `PUBLICATION_PIPELINE_REFACTOR_PLAN.md`, beginning with baseline
> preservation and the Figures 9/11 vertical slice. Keep scientific source hashes
> and all numerical caches unchanged. Centralize active figure/table presentation
> in `publication/`, retain compatibility entrypoints, and validate each phase.
> Build isolated candidates; do not publish, promote pooled Figure 7, change
> scientific definitions, or edit the manuscript. Report the completed phases,
> tests, visual checks and any explicitly retained compatibility boundaries.

## Progress record

Branch `publication-refactor` (from `6b464be`). Author decisions on 2026-10-02:
midpoint review after phases 0–1; candidates under
`publication/output/candidates/` (ignored by the existing `output/` rule; ignore
policy to be revisited); visual equivalence at a glance suffices, not pixel
identity; long checks may run locally with timings recorded; display the
first-cousin-null-mean distance as `d_CP` instead of `d_NP`.

### Phase 0 — done (2026-10-02)

`publication/output/baseline_20261002/`: HEAD/status, SHA-256 of all 1,822
back-end output files and tracked back-end source, a 562 MB tar of both
`output/` trees, and timed logs from `run_baseline.sh` (about 3 min total):

| Check | Result | Time |
|---|---|---|
| `terminal_pareto/publication_release.py --verify` | 66 files verified | <1 s |
| `full_tree_pareto.cross_species_publication --verify-production` | pass | 3 s |
| terminal release/cross-species unittests | 20/20 | <1 s |
| `full_tree_pareto` unittest discover | 55/55 | 26 s |
| terminal / partial-forest `cross_species_analysis --verify-only` | 13,828 / 9,012 assignments replayed | <1 s / 3 s |
| `terminal_pareto/test_analysis_context.py` | 13/13 | 87 s |
| `terminal_pareto/validate_pooled_migration.py` | 29/29 | 49 s |
| `full_tree_pareto.publication_build --verify` | **pre-existing failure**: production has no `release_manifest.json` (pooled Figure 7 never promoted) | — |

No verifier changed any output file (hash diff empty).

### Phase 1 — done; midpoint review 2026-10-03

- `publication/` created (see its README): notation registry, style,
  contract, adapters, shared artists/canonical panels, Figure 9/11 layout and
  captions, assembly, provenance and a cache-only CLI.
- Parity: with the old `d_NP`/`N` notation, Figures 9 and 11 rebuild with
  byte-identical TeX and pixel-identical panel PNGs and compiled pages
  (`candidates/phase1_parity_20261002/`).
- Candidates with `d_CP`/`C`: `candidates/phase1_20261002/{terminal,full-tree}-cross-species/`;
  one page each, no TeX warnings, visually checked.
- 13 new tests: scope, replayed canonical values (all 18 per family), artist
  points and connector endpoints, notation override propagation to both
  families' artists and captions, numerical hashes unchanged, solver/null calls
  blocked, integrity vs readiness, destination guards, legacy parity.
- After the phase: tracked back-end source and all 1,822 output files
  unchanged; both production verifiers and the 20 + 55 legacy tests pass.
- Review decisions (2026-10-03): the null-mean point is labeled `C`; move
  the verifier change (integrity versus readiness) to the start of phase 2;
  mixed notation across not-yet-migrated figures is acceptable meanwhile.

### Phase 2 — done (2026-10-03)

- Verifier split first (`4cc8fdc`): production/promotion checks verify frozen
  artifacts, renders and numerical inputs strictly and report stale live
  presentation source; assembly keeps the strict current-source check.
- All active presentation moved into `publication/` (`78d34c2`, `1e80023`,
  `ffdd15f`): Figures 8-11; terminal 1 amendment-6, S1-S3, Table 1 and all
  wrapper captions; pooled Figure 7 and provisional S4. Legacy entrypoints
  delegate. Full-tree presentation imports no terminal figure script.
- Parity with the previous notation: Figures 8-11 and Figure 7/S4 identical
  (TeX bytes, panel and page pixels); terminal 48/50 identical, remaining two
  PNGs differ by antialiasing only (Figure 6C 81 px, S3 57 px; PDFs identical).
  Legacy `--render-only`, numbered assembly and `--layout-only` commands were
  run on scratch run copies and reproduced every historical file.
- Findings recorded for the author: published Figure 2 and S3 used Matplotlib
  defaults (reproduced via `style.matplotlib_defaults`); pre-existing TeX
  overflow warnings in Figure 7 (8.27 pt) and standalone Table 1 (3.75 pt)
  are explicitly tolerated, any other warning fails.
- Compatibility exceptions (validator-pinned): `figS3_cross_geometry.py`,
  `fig5_table1_ce_canonical_metrics.py`. Figure 3A asset moved to
  `publication/assets/`.

### Phase 3 — done, nothing promoted (2026-10-03, `b78646c`)

- Registry: families, production roots, printed numbers, owned assets
  (disjoint per root, tested); `build --all` with `build_set.json`; set
  verification.
- `release.py`: inventories, hash-checked archives, staged transaction with
  rollback; mixed full-tree release (pooled Figure 7/S4 in, five superseded
  historical files retired, methods PDFs and Figures 10/11 preserved, Figures
  10/11 preserved-file record rewritten with history). Fixture tests cover
  success with both real verifiers, foreign-asset refusal, staged and final
  failure rollback and tampered builds. Rehearsed on a scratch copy of real
  production with the real verifiers: 5 added, 8 changed, 5 removed;
  production unchanged. The CLI only rehearses; the legacy guard remains.
- Not unified yet: terminal and Figures 10/11 promoters keep their own tested
  transaction code.

### Phase 4 — done (2026-10-03)

- `publication/README.md` (families, commands, notation recipe, boundaries,
  open author items); source maps in both package READMEs; link from
  `PUBLICATION_FIGURES_HANDOFF.md`.
- Candidates: `publication/output/candidates/phase3_20261003/` (all families,
  `d_CP` notation, verified ready); parity builds `phase*_parity_*`.
- Final checks: 35 publication tests, 16 terminal and 53 full-tree tests,
  29/29 migration checks, 13/13 profile tests, both production verifiers,
  13,828 + 9,012 assignment replays; all 1,822 back-end output files and every
  validator-pinned source byte-identical to the phase 0 baseline.

### Post-review decisions and Figure 7 release (2026-10-03/04)

- Author decisions: one universal style for all figures (Figure 2 and S3
  restyled in builds); vertical overflow accepted (recorded, not fatal);
  phylogenetic methods PDFs dropped from full-tree production for now (kept in
  the archive); promote pooled Figure 7.
- Restyled terminal candidates: `publication/output/candidates/restyle_20261003/`
  (Figure 2 and S3 inspected; other assets unchanged). Terminal production not
  re-released.
- Figure 7/S4 release: candidate `candidates/release_fig7_20261003/` (identical
  to the approved build), rehearsed, then applied with `python -m publication
  release --family full-tree-pooled --apply` at `20261004T024949483965Z`.
  Archive: `full_tree_pareto/output/legacy/releases/20261004T024949483965Z/`.
  Post-release: `publication_build --verify` passes (it failed at baseline for
  lack of a release manifest), Figures 10/11 production and archive verify,
  terminal production verifies; only full-tree production and the new archive
  differ from the baseline hashes. Snapshot:
  `publication/output/baseline_20261002/output_sha256_after_fig7_release_20261004.txt`.

### Production and back-end separation (2026-10-04)

Author decisions: production moves to `publication/output/` (one flat folder);
back-end commands compute only; figure material in runs is archived and then
purged with the old production folders after the PR; open the PR.

- `publication/output/production/`: one release manifest with a section per
  family; releases refuse stale builds and foreign assets, archive to
  `publication/output/archive/releases/`, roll back on failure. First release
  of all four families (84 files, `d_CP`, unified style) at
  `20261004T180225456854Z`.
- Back ends compute only: `terminal_pareto/analysis_pipeline.py` (formerly
  `publication_build.py`) and `full_tree_pareto/pooled_pipeline.py`; figure
  scripts reduced to analyses or removed; render/assembly/release modules
  removed (`--render-only` hooks in the hash-pinned analysis CLIs now point to
  `publication`). Display helpers moved into `publication/`. The S3
  projection table moved to the pooled run's `analysis/s3_cross_geometry/`.
- Migration validator is numerical only (24 checks); its five figure checks
  are covered by `publication/tests` (compiled pages, inventories, wrapper
  organization).
- Figure material archived (30 paths, 418 files) in
  `publication/output/archive/migrated_20261004/`; purge after the PR with
  `python -m publication.migration --purge --confirm`.

### Post-merge cleanup (2026-10-04)

- PR #2 merged with author fixes (render hooks exit 0 after successful builds;
  tests read production). Batch `migrated_20261004` purged (30 paths), which
  also removed ten formerly tracked old wrapper `.tex` files from Git.
- Batch `legacy_20261004`: figure history in both `output/legacy/` folders
  (6 paths, 242 files) archived, verified and purged; the numerical tarballs
  stay in `terminal_pareto/output/legacy/`. The validator's archived-figure
  check moved to `publication/tests/test_migration.py` (validator: 23 checks).
- Intermediate candidates removed. Maintenance tools (`migration.py`,
  `parity.py`) no longer count as presentation source; production was rebuilt
  (figures identical) and re-released at `20261004T213512971027Z` from
  `candidates/production_20261004b/`, the only candidate kept.

### Per-family presentation tracking (2026-10-05)

Figure 6B-C z-order fix released as a single-family release (`63ee580`), which
exposed cross-family staleness noise. Builds now record per-family
dependencies; all families were rebuilt (figures identical) and re-released
at `20261005T035146633045Z`; `verify --production` reports nothing stale.
53 publication tests pass.
