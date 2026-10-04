# Terminal-anchored partial-forest cross-species comparison

> **Layout change, 2026-10-04.** Figures, captions and releases moved to
> [`publication/`](../publication/README.md); production is `publication/output/production/`.
> This package is numerical only. Figure material that this record places in
> run folders (`figures/`, `publication/`, `layout_history/`) or in the old
> `output/publication/` is archived in
> `publication/output/archive/migrated_20261004/`; figure history formerly in
> `output/legacy/` (release archives, embryo-1 figures, diagnostics) is in
> `publication/output/archive/legacy_20261004/` at the same relative paths. Commands named below that drew or released figures are retired;
> see the package README for current commands.

Updated 2026-10-02. Numerical build, captioned assembly and visual review are
complete. The author approved production **Figures 10 and 11**, now released
in `output/publication/`. This analysis does not replace Figure 7 or its
supplement and does not edit the manuscript.

## Production checkpoint

- Figure 10: `output/publication/fig10_full_tree_cross_species_comparison.pdf`.
- Figure 11: `output/publication/fig11_full_tree_cross_species_overlays.pdf`.
- Both have editable TeX and panel PDF/PNG files, with the same caption-wrapper
  style as preceding publication figures. Figure 11 explicitly refers to
  Figure 10. Both are visually checked, single-page PDFs without TeX warnings.
- Promotion is additive and rollback-tested. It preserves all **15 unrelated
  production files** byte-for-byte, including Figure 7 and the supplement.
  Neither the pooled Figure 7 build nor terminal publication was promoted.
- Only these four obsolete full-tree Figure 8 assets were retired:
  `fig8_ce_full_tree_collective.pdf`, `fig8_ce_full_tree_collective.tex`,
  `fig8_ce_full_tree_collective_panel.pdf`,
  `fig8_ce_full_tree_collective_panel.png`.
- The entire preceding 19-file production folder, including those four files,
  is recoverable with verified hashes at
  `output/legacy/cross_species_releases/20261001T150611720032Z/publication/`.
  Production now has 25 files: 15 preserved, eight Figure 10/11 assets and two
  cross-species manifests. `cross_species_release_manifest.json` records the
  source run, new assets, preserved assets, retired assets and archive hashes.
- The prior unnumbered captioned assembly is preserved under the run's
  `layout_history/20261001T150500070345Z/`. All **29 numerical artifacts** remain
  byte-identical; numbering/promotion did not rerun assignments or null draws.
- **Pooled Figure 7 was promoted on 2026-10-04 without touching Figures 10/11**
  (see `POOLED_HANDOFF.md`); this release's preserved-file record was updated
  then, with the previous record kept in `preserved_files_history`. Historical note: The older
  `publication_build --publish` replaces the whole directory and now refuses
  a folder containing the Figures 10/11 release manifest. Preserve Figures 10/11 and refresh cross-species preserved-file
  hashes when replacing Figure 7 assets.

Presentation revision, 2026-10-01: every comparison panel now places a bold
`3D tracking` or `2D (XY) tracking` label directly beneath its species/molecular
configuration title. Configuration names repeat in the XY row; embryo counts
were removed from panel titles and remain in captions. The faint row-wide
geometry labels were removed. This same heading style is retrofitted to
terminal Figure 8. Prior layouts were archived with hashes; scientific caches,
endpoint anchors, assignments and metrics are unchanged.

## Author-approved scope and coverage decision

Brief pre-commit review, 2026-10-02: added and tested the Figure 7 promotion
guard above, preventing accidental replacement of the new production figures.
The scientific implementation and cached results required no changes. All 55
full-tree tests and 20 terminal release/cross-species tests pass; both production
bundles verify. The author requested a repository checkpoint of the source,
tests and handoffs. Generated figures, numerical caches and archives remain
local and Git-ignored.

Compare CE protein, CE RNA and CB AF16 RNA under 3D and XY tracking using
**only exact layerwise assignment**. The main comparison follows Figure 8's
two-row/three-column layout with Figure 7B's reference models. The companion
follows Figure 9: 3D/XY overlays in A--B and canonical metrics in C.

The first audit found 803 common measured nodes and 762 direct canonical
edges, but 41 disconnected roots, seven isolated upper islands and 64 unary
parents. Both RNA tables lack the four Figure 7 roots and 29 of its measured
internal nodes. The author explicitly allowed a bottom-up partial-tree
comparison, stopping when RNA coverage becomes unavailable, without
imputation or skipped ancestors. The implementation further anchors the
comparison to the **same 187 matched biological terminal edges as Figures
8--9**, rather than retaining arbitrary disconnected upper islands.

Ascend independently from each matched terminal through the shared measured
nodes, stopping at the **first missing canonical ancestor**. Do not resume
above a gap even when older ancestors are measured. Take the union of those
chains. Current cohort:

- 485 measured cells, 454 one-generation canonical edges, 31 boundary roots.
- 187 biological terminal leaves and 298 retained internal nodes.
- 142 observed one-child capacities and 156 two-child capacities. A missing
  sibling is not imputed or replaced by an invented slot.
- Six asynchronous contraction rounds: 187, 119, 80, 48, 18, 2 edges.
- Terminal chains have 2--6 scored ancestral transitions: 6 terminals reach
  two, 17 reach three, 62 reach four, 86 reach five and 16 reach six.

All configurations/geometries use this identical cohort, rounds and capacities.
Boundary-root identities are fixed; descendant/component membership may move
under permitted assignments. An assignment cannot cross contraction rounds.
Summing exact round optima at one common weight is an attained optimum in
their product space, **not unrestricted full-tree optimization**.

## Costs, pooling and references

- CE tracking embryos 1/2/3: cutoffs 255/247/225. CB AF16 p1/p5:
  `210519ZZY0874p1` at 148 and `210519ZZY0874p5` at 156.
- Last retained measured positions; cutoff-only single rows excluded by
  existing loaders. CE scale 0.1625; CB stored tracking values. Coordinate
  order is made explicitly x,y,z, converting CB loader z,x,y; XY drops z.
- Each embryo's pairwise distances are divided by its natural total on the
  frozen **454-edge partial forest**, then averaged equally within species.
  Natural pooled travel is one. Denominators remain fixed for every assignment
  and reference; 3D/XY have separate totals. No synthetic average coordinates.
- CE protein/RNA use identical travel matrices. Euclidean molecular distance
  uses frozen top-20 z-scored protein features or shared 20 RNA TFs on stored
  values, with no new feature selection, RNA rescaling or imputation. This
  differs from terminal cosine distance and retains the full-tree convention.
- Global analytic first-cousin-null SDs scale optimization in every round.
  Exact full-forest means include the fixed internal contribution; variances
  and pooled cross-embryo covariance come from shared terminal permutations.
- 301 endpoint-inclusive weights; nested 1,201-weight checks for all six
  cases. Primary endpoint optima use checked diminishing secondary-objective
  perturbations; no exhaustive tied-optimum enumeration claim.
- First/second/third cousins shuffle only the 187 terminal identities within
  canonical ancestor groups two/three/four transitions back. Position and
  molecular state move together; internal identities remain fixed; all 454
  edges are scored. Canonical ancestry may include unscored ancestors, as in
  Figure 7; this grouping does not create a scored skipped-ancestor edge.
- Random rebuild fixes 31 roots, the measured-state inventory and the
  observed one-/two-child capacities. Randomly attach nonroot internal nodes
  to available slots in a sampled order, then assign terminals to remaining
  slots. This is a descriptive sampler, not uniform enumeration of forests.
- 10,000 draws per reference, seeds 42/43/44/45; display deterministic 1,000
  subsets, mean markers use all draws except analytic first-cousin mean.
  Identical permutation/topology draws are paired across all six cost spaces.
- No spanning forest, top-down, bottom-up, paired reconstruction or Gaussian.

Endpoint normalization maps travel/state optima to (0,1)/(1,0) independently
per configuration. Observations and references share that configuration's
anchors; raw costs and assignments remain unchanged. Natural-lineage distance
and null-mean distance use the same closest **attained sampled assignment**,
not a line-segment projection or maximum retention. Canonical position is arc
length from the travel optimum divided by total sampled front arc length.
Retention compares biological parents over all 454 edges, not slot identity.

## Execution and cache protection

Run: `output/runs/cross_species_layerwise_v1/terminal_anchored_20260930/`.

```bash
env OPENBLAS_NUM_THREADS=1 conda run --no-capture-output -n dev python -m full_tree_pareto.cross_species_analysis
env OPENBLAS_NUM_THREADS=1 conda run --no-capture-output -n dev python -m full_tree_pareto.cross_species_analysis --verify-only
# Figures 10/11 (replaces --render-only, numbered assembly and --publish):
python -m publication build --family full-tree-cross-species --output-dir publication/output/candidates/NAME
python -m publication release --build publication/output/candidates/NAME --rehearse   # then --apply
```

Scientific caches pin cohort, data/source hashes, dependency versions, settings
and artifact hashes. Each completed case can resume; changed scientific code
requires a new run identity. Full parent arrays, leaf permutations and random
forests are saved and replay-validated. Plot/caption source identity is separate
from scientific identity. `fig_cross_species.py` reuses existing overlay and
canonical metric artists without modifying terminal presentation sources.

Output inventory:

- `analysis/inputs.json`, `nodes.csv`, `coverage.csv`, `matrices.npz`:
  frozen input identity, exclusion audit, states, geometry and pooled matrices.
- Six case NPZ/JSON caches and endpoint-tie audit JSON files.
- `references.npz`: paired cousin permutations and full random parent arrays.
- `fronts.csv`, `metrics.csv`, `null_clouds.csv`, `convergence.csv`,
  `provenance.json`: replayed results, canonical summaries and resolution check.
- `figures/partial_forest_cross_species_comparison.{png,pdf}` and
  `partial_forest_cross_species_overlay.{png,pdf}`, with a render manifest.
- `publication/`: numbered Figure 10/11 captioned one-page PDFs, editable TeX,
  copied panels and assembly manifest; the same eight assets are in production.
- `layout_history/`: hash-checked preceding working layouts, not a backup
  outside this checkout. All generated outputs/caches are Git-ignored.

## Completed validation checkpoint

- 16 focused tests pass: gap-stopping closure,
  live coverage/XY geometry, mixed-capacity rounds, brute-force product optimum,
  exact moments including fixed internal costs, cycles/roots/round guards,
  reproducible random rebuild, cache guards, six-panel and canonical artists,
  caption scope/numbered labels, hash-checked archives and live equivalence
  to the existing `pooled_analysis.layerwise` heuristic on the partial forest;
  additive release preservation, final-verification rollback and refusal of
  ambiguous Figure 8 retirement, pinned release conflicts or symlinks.
- A read-only endpoint probe confirmed positive travel and molecular endpoint
  spans for every configuration and geometry.
- All six 301- and 1,201-weight sweeps complete: **9,012 assignments** are
  replayed for raw objectives, biological-parent retention, round slots,
  roots/capacities/cycles, endpoint coordinates and canonical metrics.
  Analysis identity is
  `4f5f7e9fbafe6f3dcc5e94382f14b815affaba5af82eb1c6652d3725848ebcb9`;
  24 data/code source hashes are checked on every load.
- Four reference families each have 10,000 paired draws. All **240,000
  scored reference realizations** across six cost spaces replay from saved
  permutations or parent arrays. The underlying 10,000 random forests preserve
  the observed one-/two-child capacities and boundary-root identities.
- Both captioned PDFs compile to one page and were rendered and visually
  inspected. No TeX warnings, overfull/underfull boxes or clipped panels found.
- All **55 full-tree tests** pass (39 pooled/historical tests plus 16 new ones).
  All 13 terminal cross-species tests pass. At this checkpoint the accepted
  terminal release verified 57 published files plus its recorded analysis
  files; its subsequent Figure 8/9 addition verifies 66 assets, documented in
  the terminal README/handoff.
- Before production approval, a 154-file pre-build snapshot of existing
  full-tree/terminal analysis and publication artifacts was rehashed unchanged.
  The subsequent Figure 10/11 release writes only the approved cross-species
  production assets/manifests and retires the four named obsolete Figure 8
  assets. All 29 new-run numerical files and 15 unrelated production files
  remain unchanged against the separately captured pre-release hashes.

## Numerical summary and sweep-resolution sensitivity

The displayed 301-weight results are:

| Geometry | Configuration | u | d_LP | d_NP | Maximum parent retention |
|---|---|---:|---:|---:|---:|
| 3D | CE protein | 0.281 | 0.08479 | 0.12975 | 0.64537 |
| 3D | CE RNA | 0.134 | 0.11194 | 0.15838 | 0.59031 |
| 3D | CB RNA | 0.288 | 0.17928 | 0.22838 | 0.50881 |
| XY | CE protein | 0.313 | 0.11291 | 0.15229 | 0.56828 |
| XY | CE RNA | 0.153 | 0.16724 | 0.21112 | 0.47577 |
| XY | CB RNA | 0.323 | 0.19091 | 0.23915 | 0.46696 |

Consult `metrics.csv` for exact values. The dense sweeps contain all base-grid
weights and reproduce their raw costs. They do **not** make every front-dependent
summary identical:

- CE RNA XY: the denser closest attained assignment reduces d_LP by 0.002800
  and d_NP by 0.001861, with u increasing by 0.019486. This is the principal
  resolution sensitivity retained in the publication handoff.
- CE protein XY: d_LP decreases by 0.000016, d_NP by 0.000066 and u by 0.000326.
- CE RNA 3D: dense maximum retention increases by one edge (1/454 = 0.002203).
- Other distances and maximum retention are unchanged; their u changes are
  below 2e-7. The retained 301-weight display matches Figures 7--9's convention;
  the 1,201-weight results remain separate validation evidence, not a silently
  adopted resolution change. Neither sweep enumerates the discrete Pareto set.

The author approved the existing 301-weight presentation for production;
numbering and release are complete, and the source checkpoint was requested
on 2026-10-02. Manuscript integration remains a separate task. Generated artifacts/caches remain local and
Git-ignored; arrange external backup rather than relying solely on this checkout.

## Interpretation limits

The coverage-limited partial forest answers a lower-lineage assignment question,
not reconstruction of the full embryonic tree. Molecular representations and
data provenance differ, so this is not a controlled protein/RNA modality test.
Tracking replication does not estimate molecular uncertainty. Travel is a
retained-position displacement proxy, not integrated trajectory length.
Traditional embryo-tracking limitations, particularly CB z measurements, are
caption-only and provisionally worded pending experimental collaborators.
Developmental alignment and RNA measured/aggregated/imputed provenance remain
unresolved. No species efficiency ranking or calibrated significance claim.
