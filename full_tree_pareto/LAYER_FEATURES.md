# Layer-specific expression features (Track E)

Brief: `plans/track_e_layer_features.md`. Module: `layer_features.py`, tests:
`test_layer_features.py`. Runs: `output/runs/layer_features_v1/<run_id>/`.
Principles follow `terminal_pareto/METHODS_cell_state_representation.md`:
features are chosen without travel or the front, the global top-20 is the
baseline, and distances are cosine.

## E0: labels, criterion, coverage and memory

Decided with the author on 2026-10-15.

**Cohort and strata.** The pooled cosine context (`pooled_analysis.build_context`,
978 cells, identical context key to `pooled_full_tree_v1/terminal_clamped_cosine_20261014`)
and its contraction rounds, counted bottom-up. The strata are rounds 1, 2, 3
and 4, plus rounds 5-8 pooled (71 cells; rounds 6-8 alone hold 18, 10 and 5).
Round 1 holds 491 measured leaves: 275 biological terminals and 216 cells cut
off by the tracking time limit.

**Label.** The Schema D mixture of a cell: the mean one-hot Schema D class over
its terminal descendants in the complete lineage (`data/cell_lineage.json`).
Leaves cut off by the time limit take the mixture of the cells they will
produce; biological terminals are one-hot and agree with the Track A labels.
`programmed_death` (Schema D's class for deaths and for terminals without a
WormBase type, 12% of terminals) is a class; `--exclude-death` renormalizes over
the ten living classes.

**Criterion: equal-relatedness separation.** For one stratum and one feature
set:
- Every pair of cells in the stratum gets the cosine distance between their
  z-scored profiles and the Bray-Curtis distance between their mixtures.
  Ancestor-descendant pairs are excluded; they occur only in the pooled
  stratum (143 pairs).
- Pairs are grouped by lineage path length (sisters 2, first cousins 4, ...).
  Within each group, the Spearman correlation of the two distances.
- The score is the mean over groups with equal weight per path length, so
  sisters count as much as the tens of thousands of distant pairs (author,
  2026-10-15; the pair-weighted mean is reported as
  `separation_pair_weighted`). A group enters only with at least 10 pairs and
  varying mixture distances, so every feature set is scored on the same pairs.
  Every group that enters holds at least 18 pairs. Constant expression
  distances score zero.

Comparing only equally related pairs keeps features that merely mark ancestry
(lineage memory) from scoring well. Reported beside it: the sisters-only
correlation (the Track A analogue) with a 2,000-resample bootstrap over sister
pairs, and the all-pairs correlation. The reference is 200 random 20-sets from
the Track A protein HV pool. Placement is the share of random sets below a
candidate, with ties counted half.

**Lineage memory.** Among non-ancestral pairs that share their majority Schema D
class and have near-identical mixtures (Bray-Curtis ≤ 0.2): the Spearman
correlation between expression distance and path length. The null shuffles
profiles 1,000 times among the stratum's cells of the same majority class.
Reported: ρ, null mean and 95% range, excess, z and a one-sided p. All
same-majority pairs serve as a sensitivity check.

**Coverage.** For each round and stratum, on raw intensities (`s3.csv`, NaN as
0): a protein is *expressed* if raw > 0 in at least 10% of the cells, and
*variable* if also in at most 90%. `protein_coverage.csv` holds every
protein's expressed fraction and log1p variance per round and stratum, for E1.

### Outputs (`analysis/`)

| File | Content |
|---|---|
| `labels.csv` | Cell, round, stratum, terminal descendants, majority class, mixture |
| `design.csv` | Pairs per stratum and path length; which enter the criterion and memory |
| `coverage.csv`, `protein_coverage.csv` | Reporter signal per feature set and per protein |
| `criteria.csv` | Separation, sisters-only (with interval), all-pairs; random placement |
| `path_criteria.csv` | Correlation per path length |
| `memory.csv` | Lineage memory, near-identical and same-majority pairs |
| `random_sets.csv`, `random_metrics.csv` | Random reference sets and their criteria |
| `reproduction.json`, `provenance.json` | Gates; settings, source and file hashes |

Gates: the pooled context key matches the cosine run, the Track A top-20
equals the pooled run's proteins and z-scored values exactly, and the Schema D
class of each of the 275 Track A terminals is reproduced. `--verify-only`
checks hashes, reruns the gates and recomputes every table at 1e-12.

### Results (run `e0_groupweight_20261015`)

Analysis id `d10b9a13…`, 10 s on 24 cores. It supersedes `e0_20261015`, which
weighted path lengths by pair count and no longer verifies against this code.
Memory and coverage are identical in both runs.

**Separation** (score; percentile among random 20-sets in brackets; then the
pair-weighted score):

| Stratum | top-20 | HV pool | whole panel | random 20, median [5-95%] |
|---|---|---|---|---|
| 1 | 0.205 (71); 0.204 | 0.250 (94); 0.249 | 0.248 (93); 0.243 | 0.183 [0.121, 0.255] |
| 2 | 0.089 (45); 0.045 | 0.168 (97); 0.180 | 0.173 (98); 0.178 | 0.096 [0.034, 0.163] |
| 3 | 0.067 (48); 0.021 | 0.123 (87); 0.125 | 0.123 (87); 0.152 | 0.068 [0.004, 0.153] |
| 4 | 0.078 (88); 0.154 | 0.018 (54); 0.118 | 0.176 (100); 0.200 | 0.011 [−0.077, 0.108] |
| 5-8 | −0.065 (53); −0.094 | −0.098 (34); −0.053 | −0.002 (84); 0.068 | −0.073 [−0.153, 0.050] |

- The global top-20 is at the random median in rounds 2, 3 and 5-8, and above
  it only in rounds 1 (71st) and 4 (88th). Pair weighting had put it below
  random in rounds 2-3 (7th and 16th), because its weakness lies in the
  numerous distant pairs.
- The whole panel is at the 84th-100th percentile in every stratum; the HV
  pool beats random only in rounds 1-3.
- Sisters only, the top-20 is strong in round 1 (0.41 [0.29, 0.52], 95th
  percentile) and round 2 (0.23 [0.01, 0.45], 89th). Above round 2 the sister
  intervals span zero for every set.
- **Mixed generations.** A round holds cells of different depths (round 1
  spans depths 5-10). Path lengths alternate: for the top-20 in round 1,
  even paths give 0.05-0.44 and odd paths, which always join cells of
  different generations, give 0.01-0.26. A cell compared with a cousin
  of another generation is compared against a mixture of a different
  resolution, so odd paths may understate separation. Keep in mind for E1
  and E2 (noted by the author, 2026-10-15).

**Lineage memory** (near-identical pairs; excess over the null, then p):

| Stratum | top-20 | HV pool | whole panel |
|---|---|---|---|
| 1 | 0.250, 0.001 | 0.294, 0.001 | 0.286, 0.001 |
| 2 | 0.176, 0.001 | 0.297, 0.001 | 0.287, 0.001 |
| 3 | 0.182, 0.001 | 0.376, 0.001 | 0.376, 0.001 |
| 4 | 0.100, 0.15 | 0.446, 0.001 | 0.590, 0.001 |
| 5-8 | 0.183, 0.05 | 0.207, 0.03 | 0.117, 0.13 |

- Memory is present in every stratum. The top-20 carries the least everywhere
  except 5-8, as at the terminal level (Track A). For the larger sets it peaks
  in rounds 3-4 rather than at the top. Same-majority pairs give the same
  picture.

**Coverage** (variable proteins, i.e. expressed in 10-90% of cells):

| Round | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| top-20 | 11 | 8 | 8 | 8 | 2 | 2 | 2 | 0 |
| HV pool | 69 | 50 | 39 | 32 | 17 | 14 | 9 | 4 |
| whole panel | 78 | 56 | 47 | 39 | 26 | 19 | 13 | 5 |

The top-20 is nearly silent from round 5 up; the whole panel keeps 22 variable
proteins in the pooled 5-8 stratum.

### Implications for E1

- There is room: larger global sets beat the top-20 on the criterion in
  rounds 2-3 and 5-8.
- The larger sets also carry much more memory. Each per-round set needs its
  memory reported beside its separation.
- Ranking features on a stratum's own pairs and scoring on the same pairs is
  in-sample. E1 needs held-out evaluation, for example ranking on half of the
  sister families or lineage clades and scoring on the rest.
