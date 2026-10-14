# Methods: choosing the cell-state representation

Source material for the publication's methodology subsection on the
cell-state cost of the terminal Pareto analysis. It covers two recorded
add-ons:

- `feature_sensitivity.py` (Track A1): feature-set survey, run
  `output/runs/feature_sensitivity_v1/survey_20261014/`.
- `representation_criteria.py` (Track A1b): front-independent criteria, lineage
  memory and the fate-only control, run
  `output/runs/representation_criteria_v1/criteria_20261014/`.

No figure or manuscript text exists yet. Displayed symbols must come from
`publication/notation.py`. `d_SD` (below) has no notation concept yet; add one
before any figure uses it. The cache field `d_NP` is displayed as d_CP.

## Principle

The terminal Pareto front compares the natural lineage with alternative
parent-child assignments under two costs: cell travel and cell-state
(expression) distance. How close the natural lineage lies to the front is the
result. The expression representation behind the cell-state cost is therefore
chosen only by criteria that use neither travel nor the front. Fronts computed
from different representations serve only as a robustness check.

## Data and cohorts

| Dataset | Values | Cohort | Role |
|---|---|---|---|
| C. elegans protein | z-scored reporter intensities (`s3_zscore.csv`, 210 protein-strain features) | 275 pooled terminal edges (Figures 2-6) | Gold standard for C. elegans; terminal, subtree and cell-type analyses |
| C. elegans RNA | TF expression (`c_elegans_tf.csv`), 237 TFs shared with C. briggsae | 187 cross-species edges (Figures 8-9) | Terminal cross-species comparison only |
| C. briggsae RNA | TF expression (`c_briggsae_tf.csv`), same 237 TFs | same 187 edges | Terminal cross-species comparison only |

The protein data resolve single cells (275 distinct profiles among the 275
cohort children). The RNA data resolve cell types: the 187 cross-species
children carry 63 distinct profiles in each species, and cells of one type
often share a profile across distant lineages. RNA is therefore not used for
subtree or cell-type analyses, and one TF set is used for both species.
C. briggsae RNA is excluded from the feature survey, because its tracking has
larger z-axis uncertainty.

Travel matrices are the published pooled matrices: three C. elegans embryos,
existing 3D and XY, with fixed 275-edge normalization. The C. briggsae
matrices pool two AF16 embryos. Cell-state cost is the cosine distance between
a parent slot's profile and a child's profile.

## Cell types and lineage distance

- **Schema D types** (11 classes) are the labels the expression classifiers
  were trained on (`expression_embedding/elegans_protein_feature_ranking.py`).
  Terminal cells without a WormBase type are `programmed_death`.
- **Fine types** are the WormBase cell descriptions (`wormweb.desc` in
  `data/2023-06-29_entropy_cell_key_V2.csv`). The cohorts contain 71
  (protein) and 51 (cross-species) fine types. They are independent of the
  classifier labels.
- **Lineage path length** between two cells is the number of divisions
  separating them through their last common ancestor: sisters 2, first
  cousins 4.

## Candidate representations

Each modality has three candidates:

1. **top-20**: the published set, the 20 features with the largest maximum
   |coefficient| of an L2 multinomial logistic regression (C = 1, Schema D)
   on terminal cells. The full rankings are in A0.
2. **HV pool**: the top half of the panel by variance of log1p values over the
   cohort cells. There is no classifier and no k. That is 105 proteins, ranked
   over the 275-edge cells on `s3.csv`. For RNA it is 118 TFs, ranked by the
   mean variance rank in C. elegans and C. briggsae over the 187-edge cells.
3. **whole panel**: all 210 proteins or all 237 TFs.

## Criterion: fate discrimination

Sibling separation is Cohen's d of the cosine distance between terminal
sisters of different type (set T) and sisters of the same type (set S):

    d = (mean_T - mean_S) / sqrt((var_S + var_T) / 2)

with sample variances. Only cohort sisters with both cells typed are used. It
is reported for Schema D and for fine types. The Schema D value is partly
in-sample for the top-20 sets, so the fine-type value is the primary
comparison.

The 95% intervals come from 2,000 bootstrap resamples, drawing S and T pairs
with replacement within each set. As a reference, 200 random 20-feature sets
are drawn from each modality's HV pool. RNA sets are shared by both species,
and sets leaving a zero vector are redrawn.

This criterion reproduces the original justification of the protein top-20
(`elegans_protein_linear_baseline.ipynb`, all terminal cells, notebook
representation): cosine d 0.684, Euclidean d 0.637.

## Check: lineage memory

Lineage memory is similarity between close relatives beyond what their
fates imply. If a daughter's expression partly reflects its mother's, for
example through a reporter protein that persists across division, the natural
lineage gains a cell-state advantage by construction.

- **Statistic.** For all pairs of cohort cells with the same fine type, ρ is
  the Spearman correlation between their cosine expression distance and their
  lineage path length.
- **Null.** Profiles are permuted among the cells of each fine type, 1,000
  times. Cells without a fine type stay fixed. This keeps the composition of
  every type, including how closely related its members are, and breaks only
  the alignment of profiles with the lineage.
- **Reported:** observed ρ; the null mean and 95% range; the excess
  (ρ − null mean); z; and a one-sided permutation p, (1 + #{null ≥ ρ}) /
  (1 + 1,000).
- **By lineage distance.** Mean distance of same-type pairs in path-length
  bins, against the same mean under the null (median and 5-95%). Consecutive
  path lengths are merged until a bin holds at least 30 pairs.
- **By cell type** (top-20 sets). Within each fine type with at least 6
  cohort cells, Spearman ρ between distance and path length, with a one-sided
  p from 500 permutations within the type.

Memory is not necessarily an artifact. Cells of one type that come from
different lineages can reach that fate through different regulatory states.
Because memory turned out to be present in every candidate, it is quantified
through the fate-only control rather than used to exclude a candidate.

## Control: fate-only cell-state cost

Each child's profile is replaced by the mean profile of its fine type over
the cohort children. A child without a fine type uses its Schema D class.
Parents keep their measured profiles. The cost then measures the distance from
a parent's state to the fates its children must reach. Within-type,
lineage-correlated variation is removed, and fate differences remain. Travel
is unchanged.

## Fronts and front metrics

Every front, measured or fate-only, uses the published sweep
(`cross_species_analysis.analyze` via `feature_sensitivity.evaluate`). That is
301 weights, each cost scaled by its exact first-cousin-null SD, Hungarian
assignment, and checked endpoint ties. Reported metrics:

- u, d_LP, d_CP (`d_NP`) and maximum edge retention, as in Figures 5 and 9.
- **d_SD**: the distance from L to the closest sampled assignment in null-SD
  coordinates centred on L, ((T − T_L)/σ_T, (S − S_L)/σ_S). Endpoint
  normalization absorbs much of a representation's effect on the front
  (shown with a noise ladder in A1), so d_SD is the sensitive readout.
- **Natural-lineage state z**: (S_L − mean_null S) / σ_S, the lineage's
  cell-state cost against first-cousin shuffles.

## Reproduction gates

A run starts only if all of the following hold:

1. The measured top-20 cell-state matrices and front metrics reproduce the
   published caches to 1e-9. That covers the primary P0 row (also checked
   against the Figure 5 landmark code) and all six 187-edge configurations
   (CE protein, CE RNA and CB RNA in 3D and XY).
2. The notebook sibling separation of the protein top-20 reproduces to 5e-4.

## Feature survey (A1)

- **Sets.** For k ∈ {5, 10, 20, 30}, 200 random sets each from the full panel
  and from the HV pool, plus three fixed classifier-ranked sets: top-k,
  bottom-k, and bottom-k within the HV pool.
- **Configurations.** Protein on 275 edges (3D), and CE protein and CE RNA on
  the 217 RNA-covered edges, in 3D and XY.
- **Placement.** Each fixed set's metrics are placed in the random
  distributions: the share of random sets below it, plus half the ties.
- **Degenerate fronts.** Mostly silent protein panels give no first-cousin
  variance in cell-state cost. Those sets are recorded as degenerate and
  counted, not summarized.

## Parameters

| Setting | Value |
|---|---|
| Seed | 20261014 (criteria); survey 20261007 |
| Bootstrap resamples | 2,000 |
| Within-type permutations | 1,000 (pooled); 500 (per type) |
| Random reference sets | 200 per modality, k = 20 |
| Path-length bins | ≥ 30 pairs per bin |
| Per-type minimum | 6 cohort cells |
| Sweep | 301 weights |
| Distance | cosine |

## Outputs (`analysis/` of a criteria run)

| File | Content |
|---|---|
| `candidates.csv` | Feature lists per dataset and candidate |
| `cell_types.csv` | Cohort cells with founder, Schema D class, fine type, fate-only group |
| `criteria.csv` | Sibling separation with intervals; lineage memory with null; profile counts |
| `memory_curves.csv` | Memory by path-length bin, observed and null |
| `memory_by_type.csv` | Memory within each fine type (top-20 sets) |
| `random_sets.csv`, `random_metrics.csv` | Random reference sets and their point criteria |
| `fronts.csv`, `front_curves.npz` | Measured and fate-only front metrics and sweeps |
| `reproduction.json`, `provenance.json` | Gate results; settings, source and file hashes |

`--verify-only` checks hashes and source identity, reruns the gates,
recomputes every table and compares at 1e-12. Random-set criteria are checked
on a seeded sample.

## Results (run criteria_20261014)

Analysis id `17f3a18e…`. Both gates passed: the published fronts reproduce to
≤ 6e-17, and the notebook values to 0.6842 and 0.6372.

**Sibling separation** (Cohen's d with 95% interval; fine types, then Schema D):

| Dataset | top-20 | HV pool | whole panel | random 20, median [5-95%] |
|---|---|---|---|---|
| CE protein | 0.70 [0.28, 1.17]; 0.73 | 0.58 [0.18, 1.10]; 0.43 | 0.35 [−0.04, 0.77]; 0.16 | 0.41 [0.01, 0.87]; 0.36 |
| CE RNA | 2.05 [1.50, 2.93]; 1.93 | 2.06 [1.33, 3.32]; 1.19 | 2.12 [1.40, 3.38]; 1.14 | 1.84 [0.95, 2.52]; 1.08 |
| CB RNA | 1.82 [1.24, 2.55]; 1.99 | 1.78 [1.17, 2.78]; 1.06 | 1.83 [1.18, 2.87]; 1.03 | 1.71 [0.90, 2.31]; 1.07 |

- Sister pairs with fine types (same / different): protein 35 / 63, RNA 32 / 36.
- The top-20 lead appears on Schema D labels only. On fine types, protein
  top-20 and the HV pool overlap, and the RNA candidates are equal.

**Lineage memory** (ρ, then null mean and 95% range, then excess):

| Dataset | top-20 | HV pool | whole panel |
|---|---|---|---|
| CE protein | 0.355; 0.160 [0.123, 0.203]; 0.195 | 0.396; 0.050 [0.007, 0.104]; 0.346 | 0.374; 0.059 [0.013, 0.117]; 0.315 |
| CE RNA | 0.306; 0.028 [−0.022, 0.098]; 0.279 | 0.097; −0.107 [−0.153, −0.044]; 0.204 | 0.113; −0.089 [−0.136, −0.022]; 0.202 |
| CB RNA | 0.288; 0.033 [−0.016, 0.106]; 0.255 | 0.093; −0.091 [−0.135, −0.023]; 0.184 | 0.139; −0.069 [−0.110, −0.010]; 0.209 |

- Every observed ρ exceeds all 1,000 permutations (p = 1/1001, z 6.4-14.2).
- **Protein memory is large and broad.** Same-type sisters are 0.30 closer in
  mean cosine distance than the null median (observed 0.14; null 5% bound 0.36),
  and first cousins 0.28 closer. The effect fades by path length 9.
- Protein memory is significant within pharyngeal muscle (ρ 0.48, p 0.002),
  body-wall muscle (0.15, p 0.008), head hypodermis (0.28, p 0.02) and the
  large hypodermal syncytium (0.14, p 0.04).
- **RNA memory is significant by rank but small in size** (sisters within
  0.02 of the null). It comes from the hypodermis, where the atlas profiles
  split by founder: large syncytium ρ 0.96 (CE) and 0.94 (CB), head hypodermis
  0.50 and 0.55. RNA shows none in muscle.

**Fate-only control** (measured → fate-only):

| Dataset, geometry | Candidate | d_SD | Lineage state z | u |
|---|---|---|---|---|
| CE protein, 3D | top-20 | 1.12 → 1.39 | −4.83 → −3.65 | 0.243 → 0.146 |
| | HV pool | 0.23 → 1.04 | −7.40 → −4.06 | 0.440 → 0.236 |
| | whole panel | 0.66 → 1.07 | −7.08 → −4.08 | 0.396 → 0.234 |
| CE RNA, 3D / XY | all three | 0.41-0.48 → 0.26-0.42 / 1.72-2.13 → 1.86-2.24 | −3.1 to −3.8 → −2.7 to −3.8 | small changes |
| CB RNA, 3D / XY | all three | 1.68-2.09 → 1.69-1.99 / 2.21-2.78 → 2.11-2.73 | −2.5 to −3.8 → −2.1 to −4.2 | small changes |

- **Protein.** Removing memory moves the HV pool and the whole panel far more
  than the top-20, and the three candidates converge. Under fate-only costs
  the lineage stays 3.7-4.1 null SDs better than its first cousins and
  1.0-1.4 null SDs from the front.
- **RNA** is resolved by type and barely changes in either species.

**Decision.** Keep the published top-20 for both modalities.
- It has the highest fate discrimination, clearly so only on Schema D labels.
- Among protein candidates it carries the least lineage memory, so its front
  changes least when memory is removed.
- Report the fate-only front as the control for memory.

## Interpretation and limits

- The criteria rank representations by fate information and inherited
  similarity. They do not validate the cost as a measure of regulatory work.
- Fine types are curated descriptions of uneven granularity. Some group
  several cell classes; programmed deaths have no description.
- Sister-pair counts are small (98 protein and 68 RNA pairs with fine types),
  so sibling-separation intervals are wide.
- RNA profiles are shared within types, so RNA distances have many exact ties.
  RNA memory is significant by rank but small in size.
