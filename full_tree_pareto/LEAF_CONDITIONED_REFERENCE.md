# Leaf-conditioned separate-clock Gaussian sensitivity

Evaluated 2026-09-10. **Side analysis only: Figure 8, its reference, and its
publication caches are unchanged.** This conditions the official separate-clock
Gaussian, not the historical shared-clock Brownian or lognormal-mixture model.

Final manuscript review, 2026-09-16: retain the root-fixed, leaf-unconstrained
Gaussian in the manuscript and preserve this alternative in the research
record. The all-edge fit characterizes increments over the entire measured
tree; root-only simulation retains the intended forward-generative comparison.
Leaf conditioning is mathematically valid with that same fitted covariance,
but adds endpoint information and changes the question. Its lower costs are
not grounds for adoption or evidence of biological inefficiency. This side
analysis is not automatically included in the manuscript supplement.

## Scope and question

The current root-fixed reference can end at arbitrary leaf states. This test
asks how the same Gaussian process behaves when it must reproduce every
observed measured leaf's position **and** 20-dimensional expression vector.

- Fix all 504 measured leaves and the four measured root states, exactly.
- Jointly sample the remaining 496 internal states.
- Preserve the same measured topology, identities and 1,000 scored edges.
- Keep the fitted spatial and protein covariance matrices unchanged.
- Use tracking duration for spatial increments and one canonical transition
  per protein increment; retain zero cross-block covariance and zero prior drift.

The leaf observations are the actual inputs to the publication analysis:
scaled final tracked coordinates and top-20 z-scored atlas protein values.
They are not a reconstructed final embryo. Some measured leaves end before
the tracking cutoff; atlas expression is not aligned to those tracking times.
The test fixes the same *observed leaf configuration*, not a synchronized final
developmental outcome.

The covariance estimates still use all observed increments, including internal
states. Thus this isolates endpoint conditioning at a fixed plug-in calibration;
it is not held-out prediction, a leaf-only parameter fit, or a significance test.

## Exact joint conditioning

For either feature block, let \(B\) be the signed edge--node incidence matrix,
with -1 at the parent and +1 at the child. Set
\(W=\operatorname{diag}(1/\ell_e)\), where \(\ell_e=t_e\) for position and
\(\ell_e=1\) for protein. Partition nodes into free internal nodes \(U\) and
fixed roots/leaves \(F\). Then

\[
\begin{aligned}
Q &= B_U^\top W B_U,\\
M_U &= -Q^{-1}B_U^\top W B_F X_F,\\
\operatorname{Cov}(X_i,X_j\mid X_F) &= (Q^{-1})_{ij}\Sigma,
\qquad i,j\in U
\end{aligned}
\]

Because every node in a feature block has the same feature covariance factor,
the mean is a clock-weighted harmonic interpolation of fixed states. Internal
nodes are correlated; no node is independently sampled from its marginal.

For efficient exact sampling, generate independent standard-normal entries in
an edge-by-feature matrix \(Z\), solve
\(QH=B_U^\top W^{1/2}Z\), and use \(X_U=M_U+HL^\top\), where
\(LL^\top=\Sigma\). The node covariance of \(H\) is exactly \(Q^{-1}\).
Sparse factorization is reused across draws. Fixed states are copied verbatim;
leaves are not snapped back after an unconstrained simulation.

Both blocks are sampled independently. Every draw is scored by recomputing
all 1,000 parent--child differences from the complete sampled state array.
Fixed values are checked exactly in every batch. Conditional edge variance
factors are also checked not to exceed their prior clock factors.

## Results

10,000 draws per model, seed 242:

| Objective | Observed | Root-fixed mean | Leaf-conditioned mean | Leaf-conditioned central 95% interval |
|---|---:|---:|---:|---:|
| Travel | 4,465.55 | 4,831.80 | 4,065.82 | 3,948.97--4,186.39 |
| Protein | 3,102.41 | 3,603.77 | 2,935.67 | 2,892.27--2,980.62 |

Conditioning reduces the mean, rather than increasing it: travel is 8.95%
below observation and protein 5.37% below. The observation lies above both
conditional central 95% intervals. These are simulation intervals at fixed
fitted covariance, not parameter-confidence intervals or calibrated tests.

An independent 10,000-draw seed-243 run gives leaf-conditioned travel 4,066.01
and protein 2,935.80, confirming the numerical result.

### Why the costs decrease

Conditioning removes much of the fluctuation possible under root-only sampling.
The mean ratio of conditional to prior marginal edge variance factors is
0.496 for each block. This follows from the trace of the projection associated
with the 496 free nodes over 1,000 edges; it is an average of edgewise ratios,
not a claim that every edge variance or total-cost variance is halved.

The conditional mean states smoothly interpolate between endpoints. Their
travel and protein costs are 2,321.39 and 1,314.31, respectively. Adding
conditional Gaussian fluctuations raises expected costs to about 4,066 and
2,936, but not above observation. Costs at mean states are not the same as
mean costs, because the norm is nonlinear.

The conditional mean minimizes the Gaussian quadratic energy, **not** the
sum of Euclidean edge lengths used in Figure 8. It is neither an exact Pareto
optimum nor a computed minimum-travel configuration. The plot marks it only
as a deterministic smoothing diagnostic.

### Remaining adequacy failures

| Edge statistic | Observed | Root-fixed Gaussian | Leaf-conditioned Gaussian |
|---|---:|---:|---:|
| Travel CV | 0.611 | 0.575 | 0.526 |
| Travel mean, internal edges | 5.862 | 5.309 | 4.657 |
| Travel mean, terminal edges | 3.091 | 4.362 | 3.484 |
| Protein CV | 0.644 | 0.221 | 0.282 |
| Protein largest-5% squared-change share | 30.63% | 11.29% | 14.88% |
| Protein maximum edge norm | 15.427 | 7.081 | 8.367 |
| Protein mean, internal edges | 2.558 | 3.604 | 3.064 |
| Protein mean, terminal edges | 3.638 | 3.604 | 2.810 |

The protein model still underrepresents heterogeneity and even reverses the
observed ordering of internal versus terminal mean edge lengths. Fixing all
leaves does not fix where along their ancestral paths expression changes occur.
The spatial model still overpredicts terminal-edge lengths while underpredicting
internal-edge lengths. Aggregate proximity should not hide these discrepancies.

## Interpretation and decision

These results do **not** establish that development is inefficient. The model
permits arbitrary internal coordinates and continuous expression profiles,
with no cell packing, excluded volume, mechanical interactions, allowable
division geometry, regulatory constraints or requirements to reproduce observed
intermediate stages. Its expression states need not be observed biological
profiles. Endpoint agreement alone does not make those histories feasible.

The travel score is a sum of parent--child final-coordinate displacements,
not integrated cell trajectories. The conditioned model is therefore not a
test of how much extra physical travel is required to reach a final embryo.

Retain this as a separate endpoint-conditioning sensitivity. The comparison
demonstrates that fixed endpoints can reduce uncertainty enough to lower both
costs. It does not supply an automatic replacement for Figure 8, or establish
an ordering between biological and conditioned-random costs in general.

## Verification and reproduction

Seven new tests cover exact chain-bridge means/covariances, branched conditioning
on both daughters, recovery of independent increments under root-only
conditioning, invariance to unused free-node observations at fixed covariance,
input checks, publication scope and exact boundaries in both blocks, covariance
reduction, and end-to-end plotting. Together with existing clock and branch-
variability tests, all 19 tests pass.

The analytic conditional expectation of summed squared edge norms is checked
against simulation. Seed 242 gives travel 21,115.24 versus analytic 21,106.32
(Monte Carlo SE 6.76), and protein 9,302.72 versus analytic 9,303.12 (SE 1.38).
The independent small-tree tests also verify off-diagonal joint covariances,
not just total-cost moments.

```bash
env OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache \
  conda run -n dev python -m full_tree_pareto.leaf_conditioned_reference

env OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/embryogenesis_mpl_cache \
  conda run -n dev python -m full_tree_pareto.leaf_conditioned_reference \
  --seed 243 --analysis-only --output-dir full_tree_pareto/output/leaf_conditioned/seed243

conda run -n dev python -m unittest full_tree_pareto.test_leaf_conditioned_reference \
  full_tree_pareto.test_clock_reference full_tree_pareto.test_branch_variability
```

Outputs are isolated in `output/leaf_conditioned/`: draws and summaries, edge
diagnostics and quantile bands, analytic moment checks, covariance, edge
identities/variance factors, conditional mean states with fixed-boundary flags,
and `leaf_conditioned_comparison.png`. Generated outputs follow the existing
ignored-output convention; this report and the source/tests preserve the
reproducible specification and principal results.

Plot caption: Root-only versus root-and-leaf-conditioned separate-clock
Gaussian reference, with observed values in black. The scatter shows 1,000
draws per model from the full 10,000-draw distributions; histograms use all
draws. The orange diamond is the cost of the conditional mean states, not a
minimum-cost solution. Protein edge-quantile bands show central 95% intervals
across simulated replicates, conditional on fitted covariance and endpoints.
