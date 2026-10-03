# Pooled full-tree figure rebuild

Updated 2026-10-01. **Computation, rendering and visual review are complete.**
The author approved the revised figures and requested a repository checkpoint
on 2026-09-30. Publication promotion remains a separate, unperformed step.
[README.md](README.md) holds the scientific specification and build commands.

The subsequent cross-species layerwise comparison is a separate, completed
**terminal-anchored partial-forest analysis**, documented in
[CROSS_SPECIES_HANDOFF.md](CROSS_SPECIES_HANDOFF.md). It does not alter this
978-node Figure 7 build or its caches. RNA coverage limits the new matched
cohort to 485 cells/454 observed edges reached upward from 187 common terminal
edges, without imputation or skipped ancestors. The author approved these as
**Figures 10 and 11**, additively released to `output/publication/` on
2026-10-01. The original Figure 7/supplement remain unchanged. Four obsolete
full-tree Figure 8 files were retired and remain recoverable in the hashed
archive `output/legacy/cross_species_releases/20261001T150611720032Z/`.

## Current status

- Run: `output/runs/pooled_full_tree_v1/terminal_clamped_20260927/`.
- All five method caches, with 301 weights each, are saved. Layout-only builds
  replay-validate all **1,505 reconstructed forests**. The original five
  reference families and two new higher-cousin families each contain 10,000
  draws. All 20,000 new draws replay from saved leaf permutations.
- The final PDFs in that run's `publication/` were rendered and visually
  inspected on 2026-09-29 after the cousin-null revision. They are clean single-page layouts, and no TeX
  warnings were found:
  - `fig7_ce_full_tree_layerwise.pdf`: main Figure 7, with round-wise A and
    collective B.
  - `figs_ce_full_tree_heuristics.pdf`: provisional Figure S4, with five
    heuristics, all references and the constraint inventory.
- `output/publication/` retains the former single-embryo Figure 7/supplement,
  now alongside Figures 10/11 and their manifests; obsolete full-tree Figure 8
  assets are archived, not in production. The former bundle's two
  wrappers were untracked from Git on 2026-09-29, so promotion no longer
  produces tracked-file diffs.
- All 38 `full_tree_pareto` tests pass (19 pooled/build/cousin, 19 historical),
  including terminal-sampler agreement, forest/slot invariants, cached-draw
  replay, nested release hashes, archive preservation and worker options.
  Original scientific cache identities/hashes and all 1,505 parent arrays
  were rechecked; no reconstructions or original null families were rerun.
- The preceding Gaussian-main layout is preserved with hashes under
  `layout_history/20260930T033756463612Z/` (UTC stamp; local date September 29).
  A subsequent intermediate layout is also archived. Redraws now automatically
  archive the preceding working layout; Figure 7 production assets are untouched.

## Next action

Open review finding: `resume_paired` still needs to create `run/analysis/`
before assembling its final cache on a completely fresh run. The worker-default
change does not fix this separate issue; the completed existing run is unaffected.

When pooled Figure 7 publication promotion is explicitly requested, **first
adapt/coordinate promotion to preserve production Figures 10/11 and refresh
their preserved-file hashes**. The existing `publication_build --publish`
replaces the entire directory and now refuses a production folder containing
the Figures 10/11 release manifest. The historical command sequence below is
deferred until that change:

```bash
python -m full_tree_pareto.publication_build --layout-only --publish
python -m full_tree_pareto.publication_build --verify
```

Promotion archives the former bundle under `output/legacy/releases/<stamp>/`
with hashes. On failure it restores the former bundle and removes its staging
folder and partial archive. Then record the archive path here. The current
repository checkpoint includes source, caption templates, tests and handoffs;
generated figures, caches and local archives remain Git-ignored.

**Do not edit `pooled_analysis.py`** before promotion. Its source hash is part
of every cache identity, so any edit invalidates the expensive paired
sweep. Caption/layout changes belong in `publication_build.py`, which is not
hashed. Higher-cousin code lives in `cousin_references.py` with a separate
source-pinned cache under `analysis/cousin_shuffles/`. Changing it invalidates
only those add-on draws; do not rewrite its manifest to force stale reuse.
Normal builds generate missing add-on draws; layout-only requires existing
ones. Promotion verifies these draws and includes nested analysis hashes.

## Author-approved scope

- One main Figure 7. Panel A shows bottom-up contraction-round fronts. Panel B
  shows the aggregate layerwise front with first-/second-/third-cousin shuffles
  and random rebuild in an inset, following the terminal layout. It has no
  spanning-forest front or Gaussian. On September 29 the author requested
  Gaussian's removal from the main story: its phylogenetic analogy is
  nonessential, it lies far from the front/natural lineage, and it samples a
  different internal-state inventory. The supplement keeps five heuristics and
  all seven reference families, including Gaussian. Terminal-only remains excluded.
- Pool pairwise travel across embryo cutoffs 255, 247 and 225. Divide each
  embryo's distances by its natural full-tree total on the shared cohort, then
  average; never average coordinates. Molecular distance remains Euclidean
  top-20 z-scored protein distance, not terminal cosine.
- Endpoint displays reuse terminal `EndpointTransform`. Each round has its own
  anchors; the main and supplementary collective panels share aggregate
  layerwise anchors.
- Cousin shuffles use terminal-pipeline canonical ancestor groups two/three/four
  generations back, but include all 491 measured leaves. Position and protein
  state move together; internals and root identities remain fixed, with all
  974 edges scored. Second/third families have 97/55 nonsingleton groups,
  each covering 489 leaves; the other two remain fixed. Seeds are 43/44.
- The author explicitly confirmed **free forward internal simulation followed
  by terminal clamping**, NOT joint endpoint conditioning. All measured leaves
  and the four roots are fixed. This Gaussian specification is retained only
  in the supplementary comparison and record, not Figure 7B.
- Keep old artifacts and legacy code intact. New caches are generated, never
  rescaled from embryo-1 results.

## Key numbers

- Cohort: 978 measured nodes, with 487 internal, 491 leaves (275 canonical
  biological terminals), four roots and 974 edges. It is ancestor-complete and
  binary. Cutoff-flagged leaves number 466/464/466 by embryo.
- Rounds: 491, 223, 123, 66, 38, 18, 10 and 5 edges.
- Natural pooled travel is 1, and natural protein cost is 2991.6901.
  Denominators are about 4417.1253, 4407.4359 and 4403.2182. Exact
  first-cousin optimizer SDs are 0.00830310 and 5.55338249.
- Layerwise has 261 distinct sampled cost pairs and a maximum retention of
  0.5030800821.
- Reference means (travel, protein):
  - First cousin: 1.10653883, 3025.17559.
  - Second cousin: 1.26111536, 3118.05194.
  - Third cousin: 1.42626429, 3228.39034.
  - Internal-layer: 3.88086035, 4219.82283.
  - Full assignment: 3.99237622, 4688.17910.
  - Random rebuild: 3.99332308, 4688.74654.
  - Clamped Gaussian: 2.47477171, 6860.47526.
- Cross-check: the terminal pipeline's 275 cells exactly equal the
  biological-terminal subset of this tree.

## Operational notes

- Runtime default updated 2026-09-29: **24 worker processes**, shared through
  `publication_build.DEFAULT_WORKERS` by the build and checkpoint entrypoints.
  Both accept `--workers N` (positive integers). Completed computations used
  six workers; this update does not rerun them or change cache identities.
  The hash-pinned low-level build's fallback of four is unchanged, but the
  entrypoints always pass the configured count. No CPU affinity is imposed.
- The former 12-hour estimate may refer to an older 1,001-weight setup, as
  suggested by the author; it is not a verified benchmark for the current
  301-weight run. Do not use it to predict runtime on 24 workers.
- For fresh runs, use
  `python -m full_tree_pareto.resume_paired --run PATH`, which saves and
  replay-verifies each weight under `checkpoints/`. Earlier interactive
  sessions died mid-sweep; a detached user service
  (`full-tree-pareto-20260928.service`) completed it successfully.
- Environment: `OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache
  conda run --no-capture-output -n dev python ...` from the repository root.
