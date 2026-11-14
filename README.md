# Embryogenesis: lineage Pareto analysis

Does the *C. elegans* embryonic lineage balance how far cells travel against how
much their state changes? This repository places the natural lineage on the
Pareto front of two costs, travel (pooled cell displacement from three tracked
embryos) and cell state (cosine distance of 20 transcription-factor reporter
proteins), against lineages rewired by cousin shuffles. It covers the terminal
divisions, their subtrees and cell types, the full measured tree layer by
layer, and *C. briggsae*.

- **Manuscript figures:** `publication/output/manuscript/current`, built by
  `publication/v2/` ([publication/README.md](publication/README.md)).
  Figures and analyses kept in reserve: [publication/BACKLOG.md](publication/BACKLOG.md).
- **Manuscript text:** the authors' Overleaf document (not in this repository).

## Packages

| Path | Contents |
|---|---|
| [`terminal_pareto/`](terminal_pareto/README.md) | Terminal-cell analysis: the 275-edge pooled cohort, subtrees, cell types, tracking robustness, the 187-edge cross-species comparison, the cell-state representation (Track A, [methods](terminal_pareto/METHODS_cell_state_representation.md)) and the subtree checks of Track D |
| [`full_tree_pareto/`](full_tree_pareto/README.md) | Full measured tree (974 edges) by contraction round and layer, cousin references, heuristics, the per-layer fronts of Track E ([features](full_tree_pareto/LAYER_FEATURES.md), [hourglass](full_tree_pareto/LAYERWISE_HOURGLASS.md)) and the partial-forest cross-species comparison |
| [`publication/`](publication/README.md) | Figures, captions, Table 1 and the manuscript set; reads the back ends' validated caches through adapters |
| [`expression_embedding/`](expression_embedding/README.md) | Reporter feature selection (the published top-20 classifier ranking) and expression embeddings |
| [`ap_axis_analysis/`](ap_axis_analysis/README.md) | Anterior-posterior axis of the tracking data |
| `data/` | Lineage tree, tracking embryos, protein and RNA tables (large inputs Git-ignored) |
| Root modules and notebooks | Earlier exploratory work (`pareto_core.py`, `lineage_optimization.py`, notebooks); see `CLAUDE.md` |

The back ends compute and validate numerical runs under
`*/output/runs/<family>/<run_id>/analysis/`; `publication/` only draws from
them. A changed configuration gets a new run id; never rewrite manifests to
force reuse.

## Runs behind the manuscript set

| Analysis | Run | Scope |
|---|---|---|
| Terminal primary | `terminal_pareto/output/runs/pooled_tracking_v1/migration_candidate_20260920` | 275 edges; 42 subtrees in the caches, 18 of at least 24 cells shown |
| Terminal cross-species | `terminal_pareto/output/runs/cross_species_terminal_v1/pooled_comparison_20260930` | 187 edges (Figure 7) |
| Cell-state representation | `terminal_pareto/output/runs/representation_criteria_v1/criteria_metric_20261014` | Figure S2 |
| Subtree checks (Track D) | `terminal_pareto/output/runs/subtree_symmetry_v1/` `compare_20261015`, `apoptosis_20261019`, `apoptosis_subtrees_20261021`, `apoptosis_single_20261022` | Figure 4 and its ABplp note |
| Full tree | `full_tree_pareto/output/runs/pooled_full_tree_v1/terminal_clamped_cosine_20261014` | 978 cells, 974 edges (Figure 6A) |
| Per-layer fronts (Track E) | `full_tree_pareto/output/runs/layer_features_v1/` `hourglass_20261019`, `stages_20261020` | Figure 6B-C |
| Partial forest | `full_tree_pareto/output/runs/cross_species_layerwise_v1/terminal_anchored_cosine_20261014` | 485 cells, 454 edges, 31 roots (backlog 7C) |

The Euclidean runs `terminal_clamped_20260927` and `terminal_anchored_20260930`
are superseded.

## Conventions

- **Cohort and tracking:** embryo cutoffs 255, 247 and 225 are stage matches.
  Travel is each embryo's displacement divided by its natural-lineage total on
  the cohort (275, 974 or 454 edges), then averaged; coordinates are never
  averaged. The terminal cohort is 275 edges; subtrees shown have at least 24
  cells.
- **Cell state:** cosine distance everywhere (author, 2026-10-14) on the
  published top-20 features of each modality; protein is the *C. elegans* gold
  standard and RNA serves only the cross-species comparisons. The front is the
  outcome, never a selection input.
- **Optimization and nulls:** 301 weights; exact first-cousin shuffle as the
  primary null; endpoint normalization is for display only.
- **Layers:** layer = 9 − contraction round, counted top-down in text and
  figures; code keeps rounds.
- **Discipline:** decision rules are fixed before results; runs are write-once
  with source hashes, replay validation and tests, and scratch work becomes a
  run if it feeds the paper. Never edit hash-pinned modules (listed in each
  package README); add a new module with its own cache, following
  `full_tree_pareto/cousin_references.py`. New back-end solvers or samplers go
  into `publication/tests/support.py::no_experiments`.

## Environment

Use the `dev` Conda environment (`environment.yml`; `tectonic` is inside it)
and run everything from the repository root:

```bash
conda activate dev
export OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache
python -m unittest discover -s publication/tests -t .
```

Some CPUs (for example under WSL2) also need `OPENBLAS_CORETYPE=Haswell`.
Everything under `output/` is Git-ignored and needs an external backup.

**Parallel sessions** need one git worktree each. Link the Git-ignored outputs
from the main checkout and copy the tracking files (the cache checks reject
paths that resolve outside the worktree):

```bash
git worktree add ../embryogenesis-NAME -b NAME main
cd ../embryogenesis-NAME
for p in expression_embedding/results terminal_pareto/output full_tree_pareto/output publication/output; do
  ln -s "$(realpath ../embryogenesis/$p)" "$p"
done
for f in ../embryogenesis/data/*; do
  [ -e "data/$(basename "$f")" ] || ln -s "$(realpath "$f")" "data/$(basename "$f")"
done
for e in embryo1 embryo2 embryo3; do rm "data/$e" && cp -r "../embryogenesis/data/$e" "data/$e"; done
```

Runs written in a worktree land in the main checkout, so never let two
sessions write one run folder. Move untracked work out before
`git worktree remove --force`, then `git worktree prune` and `git branch -d`.
terminal-primary figure builds and `publication.tests.test_terminal_primary`
work only in the main checkout.

## History

The research tracks (A feature calibration, B figure redesign, C outlier
subtrees, D subtree symmetry, E per-layer fronts) closed by 2026-11-13. Their
briefs (`plans/`) and the session handoffs (`PUBLICATION_FIGURES_HANDOFF.md`,
`PUBLICATION_PIPELINE_REFACTOR_PLAN.md`, `terminal_pareto/HANDOFF.md`,
`full_tree_pareto/POOLED_HANDOFF.md`, `CROSS_SPECIES_HANDOFF.md`,
`LEGACY_EMBRYO1.md`) were retired on that date; their durable content is in
the package READMEs and the backlog, and the originals remain in git
(`git show 0b38600:PATH`).
