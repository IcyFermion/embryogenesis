# Publication front end

`publication/` turns validated back-end results into figures, captions, tables
and standalone TeX. It never runs experiments: `terminal_pareto/` and
`full_tree_pareto/` own datasets, optimization, null/reference generation,
metrics, numerical caches and scientific validation. Migration plan and
progress: [`PUBLICATION_PIPELINE_REFACTOR_PLAN.md`](../PUBLICATION_PIPELINE_REFACTOR_PLAN.md).

## Families

| Family | Figures | Back-end input | Production root | Status |
|---|---|---|---|---|
| `terminal-primary` | 1 amendment, 2-6, S1-S3, Table 1 | `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920` | `terminal_pareto/output/publication` | published |
| `terminal-cross-species` | 8, 9 | `.../cross_species_terminal_v1/pooled_comparison_20260930` | `terminal_pareto/output/publication` | published |
| `full-tree-pooled` | 7, provisional S4 | `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_20260927` | `full_tree_pareto/output/publication` | approved, **not promoted** |
| `full-tree-cross-species` | 10, 11 | `.../cross_species_layerwise_v1/terminal_anchored_20260930` | `full_tree_pareto/output/publication` | published |

`python -m publication list` prints this from `registry.py`, which also holds
wrapper stems, printed numbers and each family's owned production assets.

## Commands

From the repository root in `dev` (or with the `dev` env's `bin/` on `PATH`
for `tectonic`), with `OPENBLAS_NUM_THREADS=1` and a writable `MPLCONFIGDIR`:

```bash
python -m publication list
python -m publication build --family terminal-primary --output-dir PATH   # one family
python -m publication build --all --output-dir PATH                       # all families + build_set.json
python -m publication verify --build PATH                                 # a family build or a build set
python -m publication release --family full-tree-pooled --build PATH --rehearse [--keep DIR]
python -m publication.parity --candidate BUILD --reference PRODUCTION     # TeX bytes, pixel diffs
python -m unittest discover -s publication/tests -t .
```

Builds are cache-only: adapters call the back ends' validators/loaders
(replay, never solve; tests block solvers, null generation and RNG) and
write into a new or empty directory outside both back ends' `output/`. A
complete `--all` build takes about 30 s. Missing caches fail with the command
that creates them. `release --rehearse` applies the mixed full-tree release
to a scratch copy of production with the real verifiers; actual promotion is
not exposed and remains a separate, approved step.

## Layout

| Module | Responsibility |
|---|---|
| `notation.py` | Concept registry: symbol, plain form, label, definition, persisted field names |
| `style.py` | Palette, rcParams, retention colormap, export, `matplotlib_defaults` |
| `contracts.py` | Validated inputs per scope: `FrontComparisonInput`, `TerminalPrimaryInput`, `FullTreePooledInput` |
| `adapters/terminal.py`, `adapters/full_tree.py` | Only modules that import back ends; load caches after validation |
| `artists.py`, `canonical_panels.py` | Shared front artists, P* replay checks, named-row metric panels |
| `figures/` | `cross_species.py` (8-11), `terminal*.py` (1-6, S1-S3), `full_tree.py` (7, S4) |
| `captions/` | `cross_species.py`, `terminal.py`, `full_tree.py` |
| `tables.py` | Table 1 |
| `assets/fig3A_ce_null_models.tex` | Figure 3A TikZ schematic (moved from `terminal_pareto/assets/`) |
| `assembly.py` | Wrapper preamble, numbering, compilation; page, label and TeX-warning checks |
| `provenance.py` | Presentation hashes and manifests; integrity versus readiness; `stale_files` |
| `release.py` | Inventories, hash-checked archives, staged transactions, mixed full-tree release |
| `registry.py`, `__main__.py`, `parity.py` | Families and CLI; candidate/production comparison |

## Changing displayed notation

1. Edit the concept's `symbol`/`plain` (or `label`) in `notation.py`. Keys are
   semantic (`null_front_distance`), so cache fields such as `d_NP`/`d_np`
   are never renamed and numerical hashes do not change.
2. `python -m publication build --all --output-dir NEW` and inspect.
3. Earlier builds still verify for integrity but report `ready: false`.

Symbols reach every figure label, caption, Table 1 and the Figure 1/5/S1
schematics. Narrative Markdown, handoffs and the manuscript are edited
manually. Changing a metric's *definition* is a scientific change made in the
back end under a new run identity.

**Current notation (author decision, 2026-10-02/03):** the distance from the
first-cousin null mean to P* is `d_CP` (formerly `d_NP`) and the null-mean
point is `C` (formerly `N`). Production assets still show `d_NP`/`N` until a
release; with the old symbols restored, builds reproduce production (below).

## Provenance and release boundaries

- *Scientific* identity: back-end caches and validators (unchanged).
- *Presentation* identity: hashes of every `publication/` source and asset
  plus the notation snapshot (`presentation_manifest.json` per build).
- *Release* identity: production manifests. Production verifiers now check
  frozen artifacts and numerical inputs strictly and **report** changed live
  presentation source rather than failing; assembly steps still require
  current source.
- Each family owns explicit production assets; a release refuses to change
  anything else. The mixed full-tree release retires the superseded historical
  Figure 7/supplement assets listed in `release.FULL_TREE_POOLED_RETIRES`,
  preserves the three methods PDFs and Figures 10/11, and rewrites the Figures
  10/11 preserved-file record (previous record kept in its history).

## Compatibility boundaries

- Legacy entrypoints keep their commands, defaults and output locations and
  import drawing/captions from here (`plot_style.py`, `publication_wrappers.py`,
  `fig*_*.py`, both `fig_*cross_species.py` and `cross_species_publication.py`,
  `full_tree_pareto/publication_build.py`). Under the old notation their
  outputs are byte/pixel-identical to the historical ones.
- Kept byte-identical because scientific validators hash their source:
  `terminal_pareto/figS3_cross_geometry.py` (S3 projections; its `render` is
  copied to `figures/terminal_s3.py`) and
  `terminal_pareto/fig5_table1_ce_canonical_metrics.py` (its Table 1 writer is
  superseded by `tables.py` but still used by the full legacy pipeline).
  Also untouched: `front_coordinates.py`, both `cross_species_analysis.py`,
  `pooled_analysis.py`, `cousin_references.py` and other pinned analysis code.
- Historical single-embryo full-tree renderers (`fig7_fig8_*`,
  `heuristic_inventory.py`) are unchanged and not part of any family.
- The terminal and Figures 10/11 promoters keep their own tested transaction
  code; only the new mixed full-tree release uses `release.py` so far.

## Known items for the author

- Published Figure 2 and S3 were drawn with Matplotlib defaults, not the house
  style (Figure 2 via the `--layout-only` path, S3 by its script). Builds
  reproduce this through `style.matplotlib_defaults`; restyling is a decision.
- Recorded, pre-existing TeX warnings (pages still one page): Figure 7 float
  8.27 pt too tall; Table 1 3.75 pt too tall when compiled standalone.
- Whether the mixed release should also retire the three methods PDFs.
