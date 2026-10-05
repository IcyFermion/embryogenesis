# Publication front end

`publication/` turns validated back-end results into figures, captions, tables
and standalone TeX, and owns production releases. It never runs experiments:
`terminal_pareto/` and `full_tree_pareto/` own datasets, optimization,
null/reference generation, metrics, numerical caches and scientific
validation, and since 2026-10-04 they no longer draw figures at all.
Migration plan and progress:
[`PUBLICATION_PIPELINE_REFACTOR_PLAN.md`](../PUBLICATION_PIPELINE_REFACTOR_PLAN.md).

## Where things are

| Path | Contents |
|---|---|
| `publication/output/production/` | **Manuscript figures**: every family side by side, `release_manifest.json` (a section per family plus release history) and each family's `*.presentation_manifest.json` |
| `publication/output/candidates/<name>/` | Builds for review (`build --all` writes one folder per family plus `build_set.json`) |
| `publication/output/archive/releases/<stamp>/` | Previous production states, hash-checked, written by every release |
| `publication/output/archive/migrated_20261004/` | Figure material moved out of the back ends (old production folders, run `figures/`, `publication/`, `layout_history/`, render logs) |
| `publication/output/archive/legacy_20261004/` | Former `output/legacy/` figure history: old production release archives, the original embryo-1 figures, diagnostic plots |
| `publication/output/baseline_20261002/` | Pre-refactor backup tar and hash lists |
| `terminal_pareto/output/runs/*/*/analysis/`, `full_tree_pareto/output/runs/*/*/analysis/` | Numerical inputs (unchanged; read through adapters) |

Everything under `output/` is Git-ignored; keep an external backup.

## Families

| Family | Figures | Back-end input | Status |
|---|---|---|---|
| `terminal-primary` | 1 amendment, 2-6, S1-S3, Table 1 | `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920` | released |
| `terminal-cross-species` | 8, 9 | `.../cross_species_terminal_v1/pooled_comparison_20260930` | released |
| `full-tree-pooled` | 7, provisional S4 | `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_20260927` | released |
| `full-tree-cross-species` | 10, 11 | `.../cross_species_layerwise_v1/terminal_anchored_20260930` | released |

`python -m publication list` prints this from `registry.py`, which also holds
wrapper stems, printed numbers and each family's owned production assets
(disjoint, so one family's release cannot touch another's files).

## Commands

From the repository root in `dev` (or with the `dev` env's `bin/` on `PATH`
for `tectonic`), with `OPENBLAS_NUM_THREADS=1` and a writable `MPLCONFIGDIR`:

```bash
python -m publication list
python -m publication build --all --output-dir publication/output/candidates/NAME
python -m publication build --family terminal-primary --output-dir publication/output/candidates/NAME
python -m publication verify --build publication/output/candidates/NAME
python -m publication release --build publication/output/candidates/NAME --rehearse   # scratch copy
python -m publication release --build publication/output/candidates/NAME --apply      # production
python -m publication release --build SET --family terminal-primary --apply           # one family from a set
python -m publication verify --production
python -m publication.parity --candidate BUILD --reference publication/output/production
python -m unittest discover -s publication/tests -t .
```

Builds are cache-only: adapters call the back ends' validators/loaders
(replay, never solve; tests block solvers, null generation and RNG) and
write into a new or empty directory outside both back ends' `output/`. A
complete `--all` build takes about 30 s. Missing caches fail with the
back-end command that creates them. A release refuses stale builds (made
from older presentation source) and assets a family does not own, retires
assets a family dropped, archives the previous production and rolls back on
any failure.

Back-end numerical commands: `python terminal_pareto/analysis_pipeline.py`,
`python -m full_tree_pareto.pooled_pipeline`, both `cross_species_analysis`
modules and the validators (see the package READMEs).

## Layout

| Module | Responsibility |
|---|---|
| `notation.py` | Concept registry: symbol, plain form, label, definition, persisted field names |
| `style.py` | Palette, rcParams, retention colormap, export |
| `contracts.py` | Validated inputs per scope: `FrontComparisonInput`, `TerminalPrimaryInput`, `FullTreePooledInput` |
| `adapters/terminal.py`, `adapters/full_tree.py` | Only modules that import back ends; load caches after validation, prepare display coordinates |
| `artists.py`, `canonical_panels.py` | Shared front artists, P* replay checks, named-row metric panels |
| `figures/` | `cross_species.py` (8-11), `terminal*.py` (1-6, S1-S3), `full_tree.py` (7, S4) |
| `captions/` | `cross_species.py`, `terminal.py`, `full_tree.py` |
| `tables.py` | Table 1 |
| `assets/fig3A_ce_null_models.tex` | Figure 3A TikZ schematic |
| `assembly.py` | Wrapper preamble, numbering, compilation; page, label and TeX-warning checks |
| `provenance.py` | Presentation hashes and build manifests; integrity versus readiness |
| `release.py` | Production releases, rehearsal, verification, archives, transactions |
| `migration.py` | One-off archive and post-PR purge of figure material left in the back ends |
| `registry.py`, `__main__.py`, `parity.py` | Families and CLI; candidate/production comparison |

## Changing displayed notation

1. Edit the concept's `symbol`/`plain` (or `label`) in `notation.py`. Keys are
   semantic (`null_front_distance`), so cache fields such as `d_NP`/`d_np`
   are never renamed and numerical hashes do not change.
2. `build --all` into a new candidate folder and inspect.
3. `release --rehearse`, then `--apply`. Earlier builds still verify for
   integrity but are stale and cannot be released.

Symbols reach every figure label, caption, Table 1 and the Figure 1/5/S1
schematics. Narrative Markdown, handoffs and the manuscript are edited
manually. Changing a metric's *definition* is a scientific change made in the
back end under a new run identity.

**Current notation (author decision, 2026-10-02/03):** the distance from the
first-cousin null mean to P* is `d_CP` (formerly `d_NP`) and the null-mean
point is `C` (formerly `N`).

## Style and compilation policy (author decisions, 2026-10-03)

- One universal style (`style.configure`) for every figure.
- Vertical float overflow ("Float too large for page", e.g. Figure 7 by
  8.27 pt, standalone Table 1 by 3.75 pt) is recorded, not fatal; release can
  use two-page rendering. Every other TeX warning fails the build.

## Provenance boundaries

- *Scientific* identity: back-end caches and validators (unchanged).
- *Presentation* identity, per family: hashes of the `publication/` sources
  and assets that family's build depends on (derived from its registry entry
  points through `publication` imports; `Family.assets` for non-Python files)
  plus the notation snapshot (`presentation_manifest.json` per build). A
  change outside a family's dependencies does not make it stale.
- *Release* identity: `production/release_manifest.json` records, per family,
  the build, presentation identity, released files and the numerical input
  hashes; `verify --production` checks files and inputs strictly and reports,
  per family, presentation source that changed since release.

## Remaining back-end exceptions

- `terminal_pareto/figS3_cross_geometry.py` and
  `terminal_pareto/fig5_table1_ce_canonical_metrics.py` stay byte-identical
  because validators hash their source. They still write figure side-files
  when run; the numerical pipeline sends those to a discarded folder. S3 is
  drawn here from the projection table in the pooled run's
  `analysis/s3_cross_geometry/`.
- `terminal_pareto/plot_style.py` re-exports `style.py` for historical
  exploratory scripts. Historical single-embryo full-tree renderers
  (`fig7_fig8_*`, `heuristic_inventory.py`) are not part of any family; they
  only write to `output/legacy/rebuild/` if rerun. `terminal_pareto/output/legacy/`
  keeps the numerical pilot-output and pre-promotion run tarballs.
