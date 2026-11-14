# Publication front end

`publication/` turns validated back-end results into the manuscript's figures,
captions and Table 1. It never runs experiments: `terminal_pareto/` and
`full_tree_pareto/` own datasets, optimization, null and reference generation,
metrics, numerical caches and scientific validation, and they draw no figures.

Interface rule (author, 2026-10-02): a presentation change never edits a back
end; a scientific change becomes a new validated run that the adapters read.

## The manuscript set

The manuscript's figures are one free-standing set (official since
2026-11-12; free-standing since 2026-11-13), at
`publication/output/manuscript/current`, a symlink to the promoted set;
`history.json` beside it lists promotions. One v2 build draws every figure and
Table 1, and nothing is taken from another figure set.

| Number | Figure | Caption |
|---|---|---|
| 1 | Schematics: the authors' files `paper/FIGS/Fig1_overview.pdf` (A-C) and `FigX_distances_schematic.pdf` (D), unchanged | draft |
| 2 | Terminal Pareto front with cousin clouds (A), nested null-model schematic (B) | draft |
| 3 | Lineage similarity along the front (A), rewiring at the travel optimum T, maximum similarity M and cell-state optimum S (B-D) | draft |
| 4 | 18 subtrees of at least 24 cells: fate pies, u, and d_LP with d_CP on one 0-1 axis; left-right counterparts adjacent | draft |
| 5 | Cell types: lineage similarity by fate class along the front (A), each class's own front (B), the cost of restricting assignments within cell types (C) | released |
| 6 | Full tree by layer: every layer's front (A), position (B) and cost (C) on the birth-size axis | agreed (2026-11-05) |
| 7 | *C. briggsae* terminal front | draft |
| S1 | Terminal assignments in four tracking geometries | released |
| S2 | Cell-state feature selection (Track A) | draft |
| Table 1 | Per-subtree statistics; r lives here | released, citing Figure 4 |

The set holds one standalone wrapper per figure (`figN_<stem>.tex` with caption
and label), its compiled page and the one graphic it includes
(`<stem>_panel.pdf`), `review.pdf` joining every page, `backlog/` with the v2
page 7B-C (*C. briggsae* against *C. elegans*), `README.md` (notes for the main
text, open items) and `manifest.json` (file hashes, the presentation sources
and external files the build recorded). Everything else kept as a backup is
recorded, with its regeneration command, in [`BACKLOG.md`](BACKLOG.md).

```bash
python -m publication.v2 --manuscript --output-dir publication/output/candidates/v2_release_YYYYMMDD
python -m publication.v2.manuscript --v2-build publication/output/candidates/v2_release_YYYYMMDD \
    --output-dir publication/output/manuscript/v2_release_YYYYMMDD
python -m publication.v2.manuscript --promote publication/output/manuscript/v2_release_YYYYMMDD
python -m publication.v2.manuscript --verify            # the official set (current)
```

- Commit first. The build records the hashes of its presentation sources
  (`v2.__main__.presentation_sources`: the v2 build, captions and manuscript
  modules and every `publication` module they import, plus `v2/assets/`) and
  of Figure 1's external files. Assembly refuses a build whose sources or
  external files changed since; `--promote` refuses a stale set, a set whose
  files changed and one assembled from uncommitted code.
- `--verify` fails on a missing, extra or modified file and reports staleness.
  Any change to a module in the source list makes the official set stale:
  rebuild, assemble and promote a new set.
- Each wrapper must compile to one page that prints its own number (read from
  the `.aux` label) with no TeX warning except a tolerated float overflow.
- Captions live in `v2/captions.py` (`AGREED`, `DRAFT`, `RELEASED` and
  `BACKLOG`); change them there and assemble a new set. A wrapper edited
  inside a set fails verification.
- `python -m publication.v2 --output-dir DIR` (without `--manuscript`) writes
  draft pages with titles, for review. The build also writes the Figure 3
  edge-change tables (`fig3_edge_changes_<point>_<weight>.csv`) and the
  palette sheet.
- The manuscript text is written in the authors' Overleaf document. How the
  set's figures and captions sync with it is open (author, 2026-11-12). The
  Git-ignored `paper/` holds a local `paperv2.tex` and Figure 1's source files
  in `paper/FIGS/`, which the build reads.

**Open items.** Captions of Figures 1-4, 7 and S2 are drafts in the v1 style
awaiting the author's review; Figures 5 and S1, redrawn in the v2 style on
2026-11-13, keep the released captions (Figure 5's says natural-edge retention).
Figure 1 is the authors' file and bypasses
`notation.py`: panels B-C read "expression" and 1D writes d_LP*, which need a
manual fix. The new name for maximum edge retention is pending (below). Five
fate classes share a hue with a measure (`palette.SHARED_HUES`); the author
decides at manuscript wrap-up.

## Where things are

| Path | Contents |
|---|---|
| `publication/output/manuscript/current` | The official manuscript set (above); `history.json` beside it |
| `publication/output/candidates/<name>/` | Builds: the v2 build behind a set, family builds that regenerate backlog figures |
| `publication/output/sketches/<name>/`, `publication/output/notes/<name>/` | Sketch and reserve-note builds, made on demand (`BACKLOG.md`) |
| `paper/` (Git-ignored) | Local `paperv2.tex`; `paper/FIGS/` with Figure 1's source files |
| `terminal_pareto/output/runs/*/*/analysis/`, `full_tree_pareto/output/runs/*/*/analysis/` | Numerical inputs, read through adapters |

Everything under `output/` is Git-ignored; keep an external backup. Older
figure output (the v1 production releases, release archives, the pre-refactor
backup and the 2026-10-04 migration archive) was deleted on 2026-11-13; the
ideas in it that still work with the current back end are in `BACKLOG.md`.

## Commands

From the repository root in `dev` (or with the `dev` env's `bin/` on `PATH`
for `tectonic`), with `OPENBLAS_NUM_THREADS=1` and a writable `MPLCONFIGDIR`:

```bash
python -m publication list                                  # figure families
python -m publication build --all --output-dir publication/output/candidates/NAME
python -m publication build --family terminal-primary --output-dir publication/output/candidates/NAME
python -m publication verify --build publication/output/candidates/NAME
python -m publication.sketches --output-dir publication/output/sketches/NAME [--only KEY ...]
python -m publication.notes --output-dir publication/output/notes/NAME
python -m unittest discover -s publication/tests -t .
```

Builds are cache-only: adapters call the back ends' validators and loaders
(replay, never solve; tests block solvers, null generation and RNG through
`tests/support.py::no_experiments`, where any new back-end solver or sampler
must be added) and write into a new or empty directory outside both back ends'
`output/`. A complete `--all` build takes about 30 s. Missing caches fail with
the back-end command that creates them. terminal-primary builds and
`publication.tests.test_terminal_primary` work only from the main checkout
(the terminal context cache key holds absolute paths).

Back-end numerical commands: `python terminal_pareto/analysis_pipeline.py`,
`python -m full_tree_pareto.pooled_pipeline`, both `cross_species_analysis`
modules and the validators (see the package READMEs).

## Layout

| Module | Responsibility |
|---|---|
| `v2/` | The manuscript set: `__main__.py` (build), `panels.py` (every figure's drawing), `data.py` (loaders on the adapters and add-on runs), `captions.py`, `manuscript.py` (assemble, verify, promote), `assets/null_models.tex` (Figure 2B). It imports nothing from `figures/` or `captions/` |
| `palette.py` | One colour per concept (measures, null shuffles, fate classes, feature-figure categories, neutrals) and the palette sheet; every v2 panel draws from it |
| `notation.py` | Concept registry: symbol, plain form, label, definition, persisted field names |
| `style.py` | rcParams (the set's base style too), retention colormap and export of the figure families |
| `contracts.py` | Validated inputs per scope: `FrontComparisonInput`, `TerminalPrimaryInput`, `FullTreePooledInput` |
| `adapters/terminal.py`, `adapters/full_tree.py` | Only modules that import back ends; load caches after validation, prepare display coordinates |
| `artists.py`, `canonical_panels.py` | Shared front artists, P* replay checks, named-row metric panels |
| `figures/`, `captions/` | Figure families' drawing and captions (backlog figures only) |
| `subtrees.py` | Subtree helpers shared by the adapter, the set and the families: lineage order, node layout, the ABplp caption note and its check |
| `tables.py` | Table 1 (the set's caption cites Figure 4; the family's keeps the released Figure 5) |
| `assets/fig3A_ce_null_models.tex` | Released Figure 3A TikZ schematic |
| `assembly.py` | Wrapper preamble, numbering, compilation; page, label and TeX-warning checks |
| `provenance.py` | Presentation hashes and build manifests; integrity versus readiness |
| `registry.py`, `__main__.py` | Figure families and the CLI (`list`, `build`, `verify --build`) |
| `sketches/` | Exploratory figures from the October review; outside the registry |
| `notes/` | Reserve notes for reviewer questions (Note A: the ABplp apoptotic-cell endpoint); statements checked against hash-checked runs |

## Figure families

The registry's families render the released figures (v1 numbering) from the
same caches as the set. Since 2026-11-13 they only regenerate backlog figures;
`python -m publication list` prints them.

| Family | Figures | Back-end input |
|---|---|---|
| `terminal-primary` | 1 amendment, 2-6, S1-S3, Table 1 | `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920` |
| `terminal-cross-species` | 8, 9 | `.../cross_species_terminal_v1/pooled_comparison_20260930` |
| `full-tree-pooled` | 7, provisional S4 | `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_cosine_20261014` |
| `full-tree-cross-species` | 10, 11 | `.../cross_species_layerwise_v1/terminal_anchored_cosine_20261014` |
| `full-tree-layers` | S5 | `full_tree_pareto/output/runs/layer_features_v1/` runs `hourglass_20261019` and `stages_20261020` |

The registry also holds wrapper stems, printed numbers and each family's
output files (disjoint across families) and non-Python assets.

## Changing displayed notation

1. Edit the concept's `symbol`/`plain` (or `label`) in `notation.py`. Keys are
   semantic (`null_front_distance`), so cache fields such as `d_NP`/`d_np` are
   never renamed and numerical hashes do not change.
2. Build, assemble and promote a new manuscript set and inspect it.

Symbols reach every figure label, caption and Table 1. Figure 1 (the authors'
files), narrative Markdown and the manuscript text are edited by hand.
Changing a metric's *definition* is a scientific change made in the back end
under a new run identity.

- **Notation (author, 2026-10-02/03):** the distance from the first-cousin
  null mean to P* is `d_CP` (formerly `d_NP`) and the null-mean point is `C`
  (formerly `N`).
- **Lineage similarity (open, Track B3):** the set relabels retention as
  "Lineage similarity" through `v2.__main__.NOTATION`, and v2 captions write
  it as literal text, while `notation.py`'s `retention` label and Table 1 still
  say edge retention. To adopt a name: change the label in `notation.py`, grep
  `publication/` for "edge retention", drop the override and assemble a new
  set.

## Conventions

- **Style (author, 2026-10-03):** one style for every figure. Vertical float
  overflow ("Float too large for page") is recorded, not fatal; every other TeX
  warning fails the build.
- **Comparisons (author, 2026-10-02):** visual equivalence at a glance, not
  pixel identity; compare rasterized pages, not PDF bytes (metadata changes).
- **Scope is stated, never inferred:** adapters and captions name the cohort
  (275, 187, 974 or 454 edges), per-round or aggregate anchors, analytic or
  sampled null means, the reference and partial or complete forests; never
  infer scope from a file name, figure number or row count.
- **Caption claims are checked at build time.** `terminal_canonical.NOTES`
  (ABplp's † in Figure 4) and `check_notes` fail the build when a claim stops
  holding; cite only statements backed by a stored run.
- **Figure 4 encoding (author, 2026-11-10):** d_LP (pink) and d_CP (gold, the
  first-cousin hue) share one axis, so r can be read from the two; d_LP/d_CP is
  close to r but not equal (r scales each cost by its null SD; ABplp 0.20 vs
  0.24). The build raises if the fate legend leaves the left half.

## Provenance boundaries

- *Scientific* identity: back-end caches and validators.
- *Presentation* identity of a family build: hashes of the `publication/`
  sources and assets that family depends on (derived from its registry entry
  points through `publication` imports; `Family.assets` for non-Python files)
  plus the notation snapshot (`presentation_manifest.json`). Reusing another
  family's drawing module makes it a dependency of both.
- *Manuscript set* identity: the set's `manifest.json` (above).

## Remaining back-end exceptions

- `terminal_pareto/figS3_cross_geometry.py` and
  `terminal_pareto/fig5_table1_ce_canonical_metrics.py` stay byte-identical
  because validators hash their source. They still write figure side-files
  when run; the numerical pipeline sends those to a discarded folder. S1 is
  drawn here from the projection table in the pooled run's
  `analysis/s3_cross_geometry/`.
- `terminal_pareto/plot_style.py` re-exports `style.py` for historical
  exploratory scripts. Historical single-embryo full-tree renderers
  (`fig7_fig8_*`, `heuristic_inventory.py`) are not part of any family; they
  only write to `output/legacy/rebuild/` if rerun. `terminal_pareto/output/legacy/`
  keeps the numerical pilot-output and pre-promotion run tarballs.
