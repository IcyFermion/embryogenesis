# Separate-clock branch-variability sensitivity

Run and evaluated 2026-09-10. **Exploratory sensitivity only: Figure 8 and its
official separate-clock Gaussian reference are unchanged.**

Final manuscript review, 2026-09-16: retain the root-fixed separate-clock
Gaussian without the extra branch multiplier. Preserve this implementation,
tests and results as a research-record sensitivity; inclusion in a manuscript
supplement would require a separate decision. Better pooled agreement and
relative log density do not resolve the documented absolute-fit failures.

## Question and model

Does heterogeneous protein-change amplitude explain the remaining Gaussian
cost excess after removing inappropriate tracking-time scaling for protein?
The test keeps the same 20 z-scored proteins, four measured roots, measured
topology and 1,000 scored parent--child transitions. The spatial process is
unchanged. Protein increments are

\[
\Delta p_e=\sqrt{s_e}Lz_e,\quad z_e\sim N(0,I_{20}),\quad
\log s_e\sim N(-\tau^2/2,\tau^2),\quad LL^\top=\widehat\Sigma_p.
\]

The Gaussian vectors and multipliers are independent across edges and of each
other. Each edge receives one scalar shared across all proteins. Thus
`E[s] = 1`, preserving covariance, and `tau = 0` recovers the Gaussian. This is
an amplitude model, not a mechanism for gene regulation, a fitted clock, or a
lineage-specific rate model. Scalar mixing cannot change increment directions
or the fraction of an edge's squared change attributable to a given protein.

The covariance remains the uncentered empirical second moment of **training
increments**. With that covariance fixed, tau maximizes the marginal likelihood
of the full 20-dimensional training vectors, integrating out the lognormal
scale. This is a two-stage, conditional plug-in fit, not joint maximum
likelihood over covariance and tau. No total-distance or CV matching is used.

The search includes the Gaussian boundary and brackets the best point on a
0--3 grid. It rejects an upper-bound solution. Gaussian quadrature uses 512
nodes, checked against 1,024 nodes at every training and held-out fit. Singular
covariances are rejected without silently adding a ridge. The largest observed
full-data log-density quadrature discrepancy was below `5e-14`.

Independent edge increments on the validated forest uniquely specify all
descendant states when the four roots are fixed. The sensitivity directly
simulates their norms, an exact shortcut for these edge-additive costs; it does
not independently resample node states. The official sampler validates the
scope and supplies the unchanged travel draws. Gaussian and mixture protein
draws share Gaussian vectors, with a separate RNG stream for mixing scales.
These are distributionally equivalent to the publication Gaussian marginal,
not bitwise copies of its protein draws.

## Full-data results

10,000 draws, seed 242; fitted `tau = 1.22525`. The analytic mean-norm multiplier
is `exp(-tau^2/8) = 0.828901`. Independent quadrature gives Gaussian and mixture
expected protein totals of 3,603.921 and 2,987.293, respectively.

| Statistic | Observed | Gaussian simulation mean | Mixture simulation mean |
|---|---:|---:|---:|
| Protein total | 3,102.415 | 3,604.314 | 2,987.547 |
| Relative total excess | -- | +16.18% | -3.70% |
| Total squared protein change | 13,621.607 | 13,624.406 | 13,626.511 |
| Edge-length CV | 0.644 | 0.221 | 0.724 |
| Largest 5% share of squared change | 30.63% | 11.29% | 36.29% |
| Maximum edge length | 15.427 | 7.075 | 20.459 |
| Mean internal-edge length | 2.558 | 3.604 | 2.988 |
| Mean terminal-edge length | 3.638 | 3.604 | 2.987 |

The mixture's protein total central 95% simulation interval is
2,856.036--3,125.993, containing the observation. However, its CV interval is
0.663--0.805 and top-5% squared-share interval is 0.313--0.427, both slightly
above observation. Thus matching the total more closely does not establish
distributional adequacy. The mixture overcorrects dispersion and cannot
represent the different internal and terminal means with one common rate law.

The independent seed-243 run gives protein means 3,603.393 and 2,987.161;
mixture interval 2,856.766--3,125.065. The conclusion is numerically stable.
Travel is identical between the two models within each run by construction.

All intervals are conditional plug-in simulation intervals, not confidence
intervals incorporating fitting uncertainty and not calibrated significance
tests. The mean reduction is mathematical: preserving second moments while
increasing amplitude heterogeneity reduces expected unsquared distance.

## Held-out subtrees: better density, inadequate local totals

Eight subtrees rooted at canonical depth three are held out in turn: ABal,
ABar, ABpl, ABpr, E, MS, C and P3. The 992 edges internal to those subtrees are
tested exactly once. Their eight incoming stem edges are not tested. Each
fold also purges the held-out subtree's incoming edge from training, ensuring
that no training edge touches any held-out measured node. Both covariance and
tau are refitted using training edges only. Each fold has 2,000 simulated
held-out replicates, with its seed recorded in `fold_manifest.csv`.

Feature selection and atlas z-scoring remain fixed. This is conditional
within-atlas validation, not independent-embryo replication or validation of
the feature-selection pipeline. Biological and measurement dependence may
remain even after shared nodes are purged.

| Held-out subtree | Test edges | Mean log-density gain, nats/edge | Gaussian total error | Mixture total error |
|---|---:|---:|---:|---:|
| ABal | 192 | +4.89 | +7.1% | -12.3% |
| ABar | 164 | +20.09 | -3.1% | -20.5% |
| ABpl | 172 | +6.77 | -7.2% | -22.7% |
| ABpr | 192 | +5.79 | +8.2% | -10.4% |
| E | 30 | +8.35 | +133.8% | +95.2% |
| MS | 148 | +7.26 | +94.4% | +63.1% |
| C | 60 | +3.88 | +12.8% | -6.8% |
| P3 | 34 | +5.15 | +85.2% | +54.1% |

The mixture improves full-vector predictive log density in all eight folds,
averaging +8.31 nats per test edge. Fitted tau ranges from 1.183 to 1.262.
These gains are descriptive scores, not independent observations for a
significance test. Density improvement and total-distance accuracy are
different criteria: tail and near-origin density can improve even while the
mean norm is less accurate.

The mixture still misses seven of eight held-out total intervals (only C is
covered); the Gaussian also covers only one, ABar. It makes AB-subtree mean
totals less accurate while reducing, but not removing, overprediction in the
non-AB subtrees. Consequently the improved pooled total partly masks systematic
lineage-group errors. An exchangeable branch-amplitude distribution alone is
not sufficient for a generally adequate developmental reference.

## Source/outlier sensitivity

The audit ranks **raw protein norms**, not the historical time-standardized
outliers. The top 50 edges contribute 30.63% of squared change; their median
dominant-protein fraction is 74.33%, versus 50.04% across all edges. This flags
the importance of feature-specific changes, which a shared scalar does not
explicitly model; it is not by itself a test of directional Gaussianity.

Seventeen edges touch a selected-protein missing value in the upstream atlas;
three are among the top 50. Together the seventeen contribute only 2.58% of
raw squared change. Excluding them from **fitting**, while still predicting
the original 1,000-edge scope, gives tau 1.19467 and analytic mixture total
3,000.391 (Gaussian 3,586.393). The overall result is therefore not driven
solely by the known missing-source edges. This is a robustness check, not a
claim that missingness is random or that excluding these edges repairs the
atlas.

`raw_outlier_windows.csv` records parent and child reporter windows, raw
adjusted means and atlas values for the dominant protein on the top 20 edges.
All 40 parent/child entries have raw reporter rows. Raw adjusted means are
not identical to processed atlas values; no equality is assumed, and this
audit does not validate every aggregation/background-correction step. Reporter
times are not aligned to the spatial tracking timeline. No preprocessing was
changed and no raw outlier was declared either artifact or regulatory event.

## Decision and next step

Retain this as supplementary sensitivity, not an automatic replacement in
Figure 8. It supports branch-amplitude heterogeneity as an important source of
the Gaussian offset, but does not establish a biological mechanism or lineage
efficiency. Better relative predictive density is not adequate absolute fit.

Before adding more parameters, investigate the systematic AB versus non-AB
and internal versus terminal residuals and the feature-specific raw changes.
These diagnostics would justify or reject a small, biologically motivated
rate grouping more clearly than tuning tau to the pooled total. Do not add
post-hoc lineage rates solely to make every observed total fall inside a cloud.

## Reproduction and outputs

```bash
env OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache \
  conda run -n dev python -m full_tree_pareto.branch_variability_sensitivity

env OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache \
  conda run -n dev python -m full_tree_pareto.branch_variability_sensitivity \
  --seed 243 --analysis-only \
  --output-dir full_tree_pareto/output/branch_variability/seed243

conda run -n dev python -m unittest \
  full_tree_pareto.test_branch_variability full_tree_pareto.test_clock_reference
```

Outputs are isolated in `output/branch_variability/`: fits and covariance,
draws and summaries, edge-level diagnostics and quantile bands, fold manifest,
held-out log densities and predictive diagnostics, edge identities, source
audit and raw outlier windows, plus `branch_variability_comparison.png`.
Generated outputs follow the repository's existing ignored-output convention;
this report and the implementation/tests preserve the reproducible results.

The plot shows full-data protein totals, edge quantile bands, held-out total
ratios and held-out log-density gains. Bands and error bars are central 95%
conditional simulation intervals. The black line denotes the observed total
or edge distribution, and a held-out total ratio of one denotes exact agreement.
Positive log-density gains favor the mixture, not necessarily its predicted
total. Neither model is asserted to be a calibrated embryogenesis null.
