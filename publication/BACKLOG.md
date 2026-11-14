# Backlog: figures and analyses kept in reserve

Backups that are not in the manuscript set (`publication/output/manuscript/current`; see
[README.md](README.md)), ready if a reviewer or co-author asks. Figure files are not kept:
regenerate them on demand. Every command below was rerun against the current back end on
2026-11-13 (cosine cell-state distance, pooled tracking, current runs) and works. Released
figure numbers are the v1 numbering; the set's numbers are given where a figure moved.

## In the manuscript set's `backlog/`

- **Figure 7B-C**, *C. briggsae* counterparts of Figures 3A and 6B (the former 8B-C), with its
  caption. Moved out of Figure 7 on 2026-11-10: gaps in data resolution between the species make
  the comparisons hard to justify. Numbers: 7B lineage similarity at M 0.74, T 0.50, S 0.21; 7C
  d_LP 0.07, 0.04, 0.07, 0.09, 0.03 and lowest share of edges differing 22, 17, 52, 43 and 26%
  for layers 4-8. A *C. elegans* reference was considered and rejected (2026-11-09): on the same
  187 cells CE RNA gives d_LP 0.23 at layer 6 against 0.075, and the weaker *C. briggsae*
  tracking (z above all) makes the difference hard to claim.

## Released figures not in the set

`python -m publication build --family FAMILY --output-dir publication/output/candidates/NAME`
(cache-only, about 30 s for `--all`).

| Released | What it shows | Family | Why it is not in the set |
|---|---|---|---|
| 1 amendment | Endpoint-normalized coordinates | `terminal-primary` | Superseded by Figure 1D (the authors' file) |
| 2 | Terminal front, single panel | `terminal-primary` | Redrawn as Figure 2A |
| 3A | Nested null-model schematic | `terminal-primary` | Is Figure 2B |
| 3B | Tree distance of the edges that change along the front | `terminal-primary` | Not placed (once planned for the supplement) |
| 3C | Structural trade-off: retention along the front | `terminal-primary` | Not placed; Figure 3A shows lineage similarity |
| 4 | Subtree map, 42 subtrees at the 12-cell floor | `terminal-primary` | Superseded by Figure 4 (24-cell floor) |
| 5A-B | Canonical subtree summary | `terminal-primary` | Superseded by Figure 4 and Table 1 |
| 6, S3 | Cell types; four tracking geometries | `terminal-primary` | Redrawn in the v2 style as Figure 5 and S1 (same data and captions) |
| S1 | Cousin-relative distance r by subtree | `terminal-primary` | Retired: Figure 4 draws d_LP and d_CP, r stays in Table 1 |
| S2 | Per-cell travel and cell-state changes by cell type at T, M and S | `terminal-primary` | Not placed |
| 7A | Full tree, layerwise fronts by round | `full-tree-pooled` | Redrawn as Figure 6A |
| 7B | Aggregate layerwise front with cousin clouds and a random-rebuild inset | `full-tree-pooled` | Not placed |
| S4 | Provisional full-tree heuristics: five reconstructions, seven reference families (with the terminal-clamped Gaussian), inventory; terminal-only withheld | `full-tree-pooled` | Not placed |
| 8, 9 | Terminal cross-species comparison and overlays (187 edges) | `terminal-cross-species` | *C. briggsae* vs *C. elegans* (2026-11-10); the *C. briggsae* 3D caveat stays provisional |
| 10, 11 | Partial-forest cross-species comparison and overlays (454 edges, 31 roots; not a complete tree) | `full-tree-cross-species` | *C. briggsae* vs *C. elegans* (2026-11-10) |
| S5 | Which travel cost moves in the mid-embryo layers | `full-tree-layers` | Not placed (Figure 6B-C carry the claim) |

## Sketches

`python -m publication.sketches --output-dir publication/output/sketches/NAME --only KEY`. They
use round numbering (layer = 9 − round) and say "edge retention".

- `full-tree-rounds` (10 figures): the full tree as eight Figure-2-style round panels, a grid,
  and u, d_LP, d_CP and retention by round. Input: the pooled run's
  `analysis/round_references/`. Figure 6A covers the per-layer fronts.
- `cross-species-rounds` (6): the partial forest by round in the Figure 8/10 layout, with
  per-round metrics against the aggregate; the degenerate round 6 is left out.
- `cross-species-partitions` (8): the 187-edge terminal comparison by subtree (map and metrics,
  3D and XY) and by cell type (retention along the front, within-type fronts, metrics). Input:
  `terminal_pareto/cross_species_partitions.py`, written to
  `cross_species_terminal_v1/pooled_comparison_20260930/analysis/partitions/`. Author decisions
  (2026-10-05): *C. elegans* fates are assumed for matched *C. briggsae* cells (caption caveat);
  within-type fronts only for types of at least 15 cells. It uses the 12-edge subtree floor (25
  subtrees; the 24-cell floor would leave 12).
- `layer-hourglass` (5): Track E candidates A-E (stacked, by layer, headline, dip forest,
  decomposition). Figure 6B-C came out of this work.
- `layer-hourglass-compact` (3): F (compact), G (decomposition) and H, an approximate timeline
  against Large et al. 2025 Figure 5E, kept in reserve.

## Reserve note

- **Note A** (`python -m publication.notes --output-dir publication/output/notes/NAME`): ABplp's
  distance from its front comes from one apoptotic cell's track endpoint (ABplpppapp), with the
  mirror subtree and cells as controls. Its sentences are checked against the Track D runs at
  build time.

## Ideas with data but no drawing code

- **Tracking-transfer companion** for the terminal cross-species comparison: each single
  embryo's optima scored under pooled costs. The per-embryo permutations are in
  `pooled_comparison_20260930`'s `assignments_and_costs.npz`; the drawing code was retired on
  2026-10-04.
- **Track A feature-set survey** (top-k against random k-sets from the panel and the HV pool):
  `terminal_pareto/output/runs/feature_sensitivity_v1/survey_20261014/analysis/summary.csv` (`terminal_pareto/feature_sensitivity.py`).
- **Fate-only control front** (Track A1b; every child replaced by its fine type's mean profile):
  `criteria_metric_20261014`'s `fronts.csv`.
- **S2 panel E**, same-type sister pairs against the within-type permutation null (rho 0.36
  against 0.16; dropped on 2026-11-09, the measure stays in Methods): `memory_curves.csv` of the
  same run.
- **301 against 1,201 weights**, for reviewers: the `cohort=dense` rows of both cross-species
  `metrics.csv` files (1,201 weights change no maximum retention; the largest shift is CE
  protein 3D, Δu −0.00145).

## Not kept: they depend on retired inputs

- Embryo-1 full-tree Figures 7-8 (Euclidean, 1,000 edges; the spanning-forest comparison), the
  six-method heuristics table and the three phylogenetic-reference methods PDFs (sources in
  `full_tree_pareto/methods/*.tex`). `fig7_fig8_ce_full_tree_pareto.py` and
  `heuristic_inventory.py` still write to `output/legacy/rebuild/` if rerun; their numerical
  inputs stay in `full_tree_pareto/output/`.
- Euclidean builds (the `rounds_v1` sketches, the Euclidean releases of Figures 7, 10 and 11).
- Mirrored-pair panels (released 5C, symmetry sketches, pair tags): the mirrored-pair layer was
  retired on 2026-11-04 (Track D6); the code is at `4e3bec0^:publication/figures/terminal_symmetry.py`.
- Back-end figures and diagnostics from before the 2026-10-04 refactor. All of these files,
  with the v1 production releases, were deleted on 2026-11-13.

## Parked analyses

- **Significance of sibling separation (S2).** Report Cohen's d with a p-value? Exploratory
  numbers (a scratch script; redo in an add-on module): one-sided Mann-Whitney AUC 0.76 / 0.69
  / 0.60 for the top-20, HV pool and whole panel. Plain permutation p of d: 0.0009 / 0.004 /
  0.05. Permutation within founder lineages: 0.021 / 0.11 / 0.36. Founder confound: 57 of 63
  different-type pairs are in AB, while 24 of 35 same-type pairs are in MS, C or D. The top-20
  does not beat the HV pool on the same pairs (ΔAUC 0.07, 95% CI −0.04 to 0.18). 27 of 200
  random 20-sets reach its d. The classifier was trained on these cells' Schema D labels.
  Proposal: keep d in S2C and add the within-founder p ("the top-20 separates discordant from
  concordant sisters within founder lineages (p = 0.02), as well as the HV pool does"). Methods
  must state that `cohen_d` pools the two variances unweighted (the n−1-weighted form gives
  0.68 / 0.57 / 0.34).
- **Terminal cohort of 321 edges** (Track D7, kept at 275 on 2026-11-09). 46 cells end
  embryogenesis undivided but are internal in the reference lineage (10 E, 11 each in ABar,
  ABpl and ABpr, 3 in MS). With them: u 0.204 (0.243), d_LP 0.031 (0.015), d_CP 0.141 (0.136),
  maximum similarity 0.75 (0.83), z_T/z_S −9.9/−5.1 (−9.1/−4.8), d̃_LP 3.22 (1.11). Run
  `terminal_pareto/output/runs/cohort_sensitivity_v1/embryonic_terminals_20261108`. If raised:
  find which of the 46 cells move the lineage off the front, by founder group.
- **Outlier subtrees** (Track C; dossier and literature session never run). ABplpa (17 cells),
  ABarp (15) and ABplppp (12) are all below the 24-cell floor. ABplpa's near-zero travel span is
  one swap of the dying ABplpappap with RMEV, within the MS-induced left-right asymmetry (Hutter
  and Schnabel 1995; Moskowitz and Rothman 1996). ABplppp is off its front only through the
  ABplpppapp endpoint. ABarp is unexplained: dropping its deaths gives 0, but holding the
  endpoints at birth leaves 0.73.
- **Left-right symmetry, reviewer defense** (Track D0-D5; symmetry left the paper on
  2026-11-02). No mirrored pair classifies as biological. Simulated twins: at 12-34 edges the
  lineage sets only a third to a half of u. The coarse lean is not evidence (P 0.83). A
  front-free check of the travel lean, θ = atan2(−z_S, −z_T), is below 45° in 41 of 42
  subtrees (median 32°). Mirror geometry: one reflection plane per embryo maps the 55 sided
  pairs to within 2.0-2.5 µm, against 11-13 µm left-right separation. Runs
  `subtree_symmetry_v1/pairing_20261014`, `compare_20261015`, `cells_20261015`,
  `robustness_20261015` (`terminal_pareto/README.md`).
- **Per-layer fronts, parked unless a reviewer asks** (Track E): fronts by division stage; a
  spatial-arrangement criterion; TAI per layer (Ma and Zheng 2023 phylostrata); a shared
  criterion-ranked feature set (beats the top-20 on held-out halves but carries more lineage
  memory); the mechanistic test of Schnabel et al. (each cell's distance to its fate region at
  the premorphogenetic stage).
- **1,001-weight sweeps** (approved 2026-09-28, not implemented). The exploratory code used
  1,001 weights and publication uses 301 (terminal `ae557ab`, full tree `d493499`, kept in
  `3b6c31f`); the reason is not established, so do not attribute it to runtime. Audit the
  terminal global, subtree, cell-type and tracking analyses and every full-tree heuristic;
  restore 1,001 consistently in builders, cache identities and captions under new run ids, with
  per-weight checkpoints for the full tree; keep the 301-weight caches as the baseline (the
  grids are not nested). Compare fronts, retention and metrics, and check anchors and raw-cost
  replay; a denser sweep still does not enumerate the Pareto set. Null draws do not depend on
  the resolution and can be reused only after checking cohort, objectives, reference, seeds and
  provenance. Promote only after review.
- **Full-tree reconstructions:** an audit of the terminal-only reconstruction (withheld: 998
  edges, roots not fixed); protein-rate heterogeneity models (block-diagonal, class-specific
  rates, heavy-tailed, residual bootstrap), proposed and never run.
- **Provenance check:** missing values in `s3.csv` were zero-filled before z-scoring, touching 3
  of the 50 largest outlier edges; reporter `Time` is not tracking `t`. Check whether
  `s3_zscore.csv` inherits the zero-fill.
- **Cross-species, with the collaborators:** stage alignment, RNA provenance and the wording of
  the *C. briggsae* z-axis caveat.
- **Not pursued:** a later MSpa cutoff; the ABalp/ABara fate-only lead.
