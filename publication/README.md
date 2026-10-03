# Publication front end

`publication/` turns validated back-end results into figures, captions and
standalone TeX. It never runs experiments: `terminal_pareto/` and
`full_tree_pareto/` own datasets, optimization, nulls/references, metrics,
caches and scientific validation. See
[`PUBLICATION_PIPELINE_REFACTOR_PLAN.md`](../PUBLICATION_PIPELINE_REFACTOR_PLAN.md)
for the migration plan and progress.

**Status (phase 1):** Figures 9 and 11 (cross-species overlays and canonical
metrics) are built here. All other figures still use their back-end
entrypoints; production directories are unchanged and nothing is promoted.

## Commands

From the repository root in `dev` (or with the `dev` env's `bin/` on `PATH`,
for `tectonic`); set `OPENBLAS_NUM_THREADS=1` and a writable `MPLCONFIGDIR`.

```bash
python -m publication list
python -m publication build --family terminal-cross-species --output-dir PATH
python -m publication build --family full-tree-cross-species --output-dir PATH
python -m publication verify --build PATH
python -m unittest discover -s publication/tests -t .
```

`build` loads caches through the back end's validator (replay only; tests make
solver/null-generation calls fail), renders panels, writes and compiles
one-page captioned wrappers, and writes `presentation_manifest.json`. The
output directory must be new or empty and outside both back ends' `output/`.
Missing caches fail with the command that creates them.

## Layout

| Module | Responsibility |
|---|---|
| `notation.py` | Concept registry: symbol, plain form, label, definition, persisted field names |
| `style.py` | Palette, rcParams and export conventions |
| `contracts.py` | `FrontComparisonInput`: the validated fields layouts may use |
| `adapters/terminal.py`, `adapters/full_tree.py` | Read back-end caches after their validators pass |
| `artists.py`, `canonical_panels.py` | Front overlays, P* replay checks, named-row metric panels |
| `figures/cross_species.py`, `captions/cross_species.py` | Figure 9/11 layout and captions |
| `assembly.py` | Wrapper preamble, numbering, compilation, page/label/TeX-warning checks |
| `provenance.py` | Presentation hashes, build manifests, integrity vs readiness verification |
| `registry.py`, `__main__.py` | Families, figure numbers/stems, cache-only CLI |

Shared modules import neither back end; only `adapters/` do.

## Changing displayed notation

1. Edit the concept's `symbol`/`plain` (or `label`) in `notation.py`. Keys are
   semantic (`null_front_distance`), so cache fields such as `d_NP` are never
   renamed and numerical hashes do not change.
2. Rebuild affected families into a new directory and inspect them.
3. Changing a symbol changes the presentation identity: earlier builds still
   verify for integrity but report `ready: false`.

Changing a metric's *definition* is a scientific change made in the back end
with a new run identity. Narrative Markdown, handoffs and the manuscript are
edited manually; the registry does not rewrite prose.

Current notation change (author request, 2026-10-02): the first-cousin-null
mean to P* distance is displayed as `d_CP` (formerly `d_NP`), and the null mean
point as `C` (formerly `N`) so that the subscript reads "from C to P*". Until
phase 2, back-end-built figures (Figures 1, 5, S1, Table 1, and the published
9/11) still print `d_NP`/`N`.

## Compatibility boundaries (phase 1)

- Full-tree production verification (`cross_species_publication
  --verify-production`) re-hashes the *live* renderer and caption sources
  recorded in the Figures 10/11 manifests, including the verifier itself:
  `terminal_pareto/{plot_style,fig_terminal_cross_species,fig2_fig3_ce_terminal_pareto,cross_species_publication,publication_wrappers}.py`
  and `full_tree_pareto/{fig_cross_species,cross_species_publication}.py`.
  Turning these into delegates would fail that verification, so they are
  byte-identical copies for now. `tests/test_cross_species.py::LegacyParity`
  keeps the copied palette, rcParams, labels, caveats and preamble in step.
  Replacing them with shims waits for a release verifier that separates
  artifact integrity from live-source readiness (phase 3).
- Scientific modules pinned by numerical caches (for example
  `front_coordinates.py`, both `cross_species_analysis.py`, `pooled_analysis.py`)
  are not touched; adapters import them.
