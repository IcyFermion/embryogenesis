# Tracking geometry and assignment stability

Investigation: 2026-09-14. Exploratory, global P0 only; no publication promotion.
Scope follows `README.md` and `PUBLICATION_FIGURES_HANDOFF.md`.

**Checkpoint: 2026-09-16, paused for manuscript writing at the author's request.**
The underlying CE tracking-variance investigation is sufficient for the current
drafting stage. Figure S3 presentation still needs work: the amendment is a
review artifact, not an accepted layout or settled manuscript narrative.
Keep the current numerical results and local artifacts as the checkpoint;
do not extend sweeps, resume AF16, or promote the candidate while merely
polishing prose. Revisit presentation and any necessary analysis once the
manuscript makes the scientific claim and required evidence clearer.

## Manuscript handoff

The evidence supports three distinct observations on the matched 275-edge set:

1. Global front shape and natural-lineage proximity are similar across the
   three observed tracking geometries, while exact parent assignments vary.
   Biological-parent agreement is higher than slot agreement; slots of the
   same parent must be treated as the same biological edge.
2. Candidate neighborhoods are substantially conserved. Nevertheless, the
   transferred pure-travel optima lose essentially all of their travel saving
   over nature. A small penalty relative to total cost need not be small
   relative to attainable improvement.
3. Joint optimization finds a more modest improvement shared across all three
   geometries. Requiring no increase in cell-state cost yields travel savings
   of 2.41%, 2.32%, and 2.28%; requiring at least 1% cell-state improvement
   still yields travel savings of 1.51%, 1.49%, and 1.47%.

Keep these statements conditional on this terminal assignment model, fixed
protein representation, matched set, and three measured geometries. This is
between-embryo sensitivity, not isolated tracking-error variance or a population
confidence interval. The common solution is optimal for the stated minimax
criterion within solver tolerances; it is not a unique biological optimum,
need not lie on every individual front, and has no held-out validation.
The variance analysis has not been extended across all 42 matched subtrees.

When returning to the figure:

1. Start from the drafted claim and decide which panels are needed. Rework S3
   hierarchy, density, labels, and caption; current rendering checks establish
   legibility, not author acceptance of the presentation.
2. Keep front stability, edge stability, and common improvement separate in
   the prose. Define the transfer direction and denominator wherever quoting
   a penalty; distinguish sampled fronts from complete discrete Pareto sets.
3. Preserve and verify the checkpoint before any rerun. Use the explicit
   common-solver commands below: only 0%, 0.5%, and 1% expression budgets were
   completed. The solver's default list also includes unfinished 2% and 5%
   budgets; its internal time limit is not a reliable wall-clock cutoff.
4. Request additional analysis only when needed to support the chosen claim;
   the remaining-scope list below is a backlog, not an active work plan.

### Checkpoint files and commit scope

The source commit contains the two analysis scripts, amendment renderer and
caption, this investigation record, and the README/publication handoffs.
Generated artifacts and input data remain Git-ignored and local; the commit
is not a backup of the following directories:

- `output/tracking_geometry_sensitivity/`: base and `refined/` sweeps, actual
  permutations/cost matrices, transfer/neighborhood/exclusion diagnostics,
  source hashes, and the original three-target review figure.
- `output/common_tracking_assignment/`: the three completed common solutions,
  their actual edges, solver bounds/statuses, and expression-tradeoff outputs.
- `output/tracking_s3_amendment/`: the two-page candidate, component panels,
  plotted coordinates, figure hashes, and isolated canonical-validation log.

Primary review artifact:
`output/tracking_s3_amendment/figS3_ce_tracking_amended_review.pdf`.
The accepted S3 assets under `output/publication/` are preserved. The stopped
AF16 attempt is outside this commit and has no reported optimization result.

## Main finding

Similar replicate-specific Pareto fronts coexist with appreciable differences
in biological parentage. Local candidate neighborhoods are more conserved than
the selected edges. However, transfer losses near the pure-travel endpoint are
large enough to erase the travel saving over the target natural lineage.
The evidence supports shared local geometry and a stable global trade-off
shape, but does not support a blanket claim that transferred solutions are
nearly optimal in a practically negligible sense.

`tracking_geometry_sensitivity.py` is an independent, deterministic diagnostic.
Its default results are in `output/tracking_geometry_sensitivity/`; the
fourfold denser run is in that directory's `refined/` subdirectory.
The generated PDF/PNG has three panels showing all source fronts evaluated in
each target geometry, followed by biological-parent agreement, transfer regret,
and candidate-parent neighborhood overlap. Connecting lines are visual guides;
they are not additional attainable assignments.

## Controlled comparison

- Same 275 terminal identities and natural edges as Figure S3, in reference
  lineage traversal order; 172 distinct biological parents with fixed slot
  multiplicities.
- Same predefined stage-matched cutoffs: 255, 247, 225.
- One shared top-20 protein cosine matrix. Only the spatial matrix changes.
- Parent-slot-to-child permutations preserve the original matching constraints.
  Comparisons invert the permutation and compare parent identity for each child.
- Exact first-cousin-null standard deviations, following S3.
- 301 endpoint-inclusive weights, repeated with 1,201 weights. The sweep is
  a sampling convention, not a biological distribution over trade-offs.

The loader takes each cell's last retained position. Thus the travel objective
is Euclidean displacement between these representative parent/child positions,
not integrated trajectory length, and positions are not all sampled at one time.
Between-embryo differences combine biological geometry, tracking error, and
residual developmental-stage differences. These data do not isolate measurement
noise or estimate population uncertainty from three embryos.

## Recomputed results

Numbers below use the 301-weight run unless stated otherwise. Pairwise ranges
are descriptive; six transfer directions are not six independent replicates.

| Quantity | Result |
|---|---|
| Spatial matrix Pearson correlation, distinct-parent rows | 0.9780–0.9791 |
| Median biological-parent agreement over equal weights | 75.6–78.9% |
| Biological-parent agreement at pure-travel optimum | 49.1–56.0% |
| Median agreement over 101 equally spaced canonical positions | 76.0–79.6% |
| Mean overlap of each child's 3 nearest distinct candidate parents | 67.3–70.7% |
| Mean overlap of 10 nearest distinct candidate parents | 78.7–79.5% |
| Mean overlap of 20 nearest distinct candidate parents | 83.8–84.8% |
| Median weighted transfer regret, as percent of target optimal total | 1.07–1.39% |
| Pure-travel transfer regret, as percent of target travel optimum | 16.06–23.45% |

For two independent uniformly chosen ten-parent sets from 172 parents, expected
overlap fraction is 10/172, or 5.8%. This is a combinatorial reference, not a
biologically realistic spatial null or a significance test.

The README's preliminary slot-agreement numbers should not be interpreted as
biological-edge agreement. Exchanging two slots of one parent changes the slot
permutation without changing any parent-child edges. This run gives median
slot agreement 49.5–52.4%, versus biological agreement 75.6–78.9%.
Slot tie choices can also depend on floating-point details.

At the expression endpoint, all three runs select the same biological assignment.
This follows from the shared expression matrix and solver tie handling; it is
not independent evidence of tracking robustness or uniqueness of that optimum.

### The denominator changes the interpretation of transfer loss

The small median percentages above divide by the complete, uncentered optimal
weighted cost. They do not mean only 1% of the improvement over nature was lost.

For the pure-travel endpoints:

| Source → target | Extra travel / target optimum | Extra travel / target natural-to-optimum saving |
|---|---:|---:|
| 1 → 2 | 23.45% | 122.63% |
| 1 → 3 | 19.44% | 163.66% |
| 2 → 1 | 20.01% | 124.19% |
| 2 → 3 | 16.06% | 135.16% |
| 3 → 1 | 19.95% | 123.80% |
| 3 → 2 | 19.14% | 100.05% |

All six transferred travel endpoints have greater travel cost than the target
natural lineage; 3 → 2 is only marginally worse (0.069 raw travel units).
This is a travel-axis statement, not a claim of two-objective dominance by
the natural lineage.

At matched expression budgets, bounds on median transfer travel loss also
show why the full-cost denominator can conceal substantial sensitivity:

| Source → target | Median travel loss / available natural travel saving, bounded |
|---|---:|
| 1 → 2 | 33.2–40.0% |
| 1 → 3 | 45.4–56.2% |
| 2 → 1 | 46.5–55.1% |
| 2 → 3 | 49.4–69.1% |
| 3 → 1 | 39.5–48.5% |
| 3 → 2 | 30.6–37.1% |

These are optimization bounds, not confidence intervals. The denominator is
the fixed target travel saving from its natural assignment to its unconstrained
travel optimum. It does not assert that nature satisfies each expression budget.
The 1,201-weight run tightens some lower bounds substantially (2 → 3 becomes
58.7–69.1%); exact constrained optima have not been solved. Thus no claim of
convergence of the complete front or all transfer bounds is made.

### How local are changed assignments?

For each changed child, measure separation between its two selected parents
in the target geometry, divided by the target parent's mean distance to its
five nearest other distinct parents. Pooling the six transfer directions, the
median is 1.04 local-spacing units at travel weights 0.9 and 1.0. At weight
0.5 it is 1.37, and at 0.1 it is 1.83. At the travel endpoint, the 90th
percentile is 1.63 local-spacing units.

Thus travel-endpoint replacements tend to be spatially local even though their
aggregate cost penalties are substantial. Expression-heavy solutions can exchange
more spatially distant parents. These are coordinate-neighborhood diagnostics;
they do not establish shared physical contacts or a developmental mechanism.

## Recommended rigorous measurements

### 1. Transfer at the same cell-state budget

For source assignment \(A\), let \(S(A)\) be its shared expression cost and
\(T_q(A)\) its travel cost in target embryo \(q\). Define

\[
g_q(s)=\min_{B:\,S(B)\leq s}T_q(B),\qquad
R_{r\to q}(A)=T_q(A)-g_q(S(A))
\]

This asks how much travel could be saved in the target while achieving at least
the source assignment's cell-state performance. It needs no alignment of weights
or replicate-specific arc lengths. Report raw travel, target null-SD units,
and fraction of the target's available natural travel saving; retain signs when
also reporting changes relative to nature.

The current script bounds this quantity without a mixed-integer solve. Every
exact weighted optimum supplies a lower bound on constrained minimum travel:

\[
g_q(s)\geq \sigma_{T,q}\max_{\alpha>0}
\frac{J_q^*(\alpha)-(1-\alpha)s/\sigma_S}{\alpha}
\]

Here \(J_q^*(\alpha)\) minimizes separately null-SD-scaled objectives.
Feasible assignments from the identical pooled candidate set (all three sweeps
plus nature), together with the transferred assignment itself, supply upper
bounds on \(g_q\). Subtract these bounds in reverse order to bound regret.
Numerical feasibility tolerance on summed expression is 1e-10.

A decisive next computation would solve the expression-constrained *integer*
assignment at selected shared budgets. Adding an expression-budget constraint
does not preserve the usual assignment relaxation's integrality. Weighted sums
can miss unsupported efficient solutions, so neither interpolation nor a denser
weight sweep certifies the complete discrete front. See
[Ozlen, Burton and MacRae (2014)](https://www.sciencedirect.com/science/article/pii/S0377221713006474).

### 2. Neighborhood stability and edge identifiability

Use biological-parent overlap at equal weight, equal canonical position, and
matched expression performance. Canonical positions must select actual sampled
assignments; record the achieved positions because sparse regions cannot be
matched exactly. Neighborhood overlap should count distinct parents and report
several neighborhood sizes, spatial ranks, and local replacement distances.

To test whether an exact edge is meaningful, compute its exclusion margin:

\[
m_{q,\alpha}(c)=
\min_{B:\,p_B(c)\ne p_{A^*}(c)}J_{q,\alpha}(B)
-J_{q,\alpha}(A^*)
\]

Forbid every slot of the selected biological parent for that child, then solve
the entire matching again. This accounts for compensating rearrangements and
capacity competition. A row's nearest-versus-second-nearest cost difference
does not account for these constraints.

The script computes 3,300 such solves: 275 children × three geometries × four
weights (0.1, 0.5, 0.9, 1.0). Dividing each margin by its geometry/weight's
natural-to-optimum weighted improvement gives pooled medians of 0.081%, 0.683%,
3.456%, and 1.013%, respectively. This is an algorithmic sensitivity scale,
not a biological combined efficiency score. Margins do not simply decrease
monotonically toward the travel endpoint.

At a declared additive tolerance \(\epsilon\), a selected edge is required
in every solution within \(\epsilon\) of the optimum exactly when
\(m(c)>\epsilon\), up to solver tolerance. Report stability curves across
tolerances instead of choosing one unexplained cutoff. Marginal alternatives
for different children need not coexist in a single feasible assignment.
General simultaneous cost-perturbation bounds are also available in
[Michael et al. (2020)](https://arxiv.org/abs/2005.11792).

### 3. A common solution with an explicit cost allowance

Independent majority votes for parentage can violate parent capacities.
Construct a feasible consensus instead: at each shared expression budget,
minimize worst-case scaled travel regret across geometries, or maximize retained
natural/shared biological edges subject to a specified regret allowance in each
geometry. Compare the structural gain as the allowance grows.

Evaluate mean-geometry or consensus methods leave-one-embryo-out: construct from
two geometries, evaluate on the third. This avoids presenting in-sample agreement
as transfer performance. Three folds remain descriptive, not population inference.

## Remaining scope for a full supplementary analysis

1. Separate embryo-1 cell-set effects (299 versus 275 edges) from geometry
   effects. Matching alters the feasible parent capacities as well as totals.
2. Extend the matched analysis to all 42 qualifying subtrees, with tree-distance,
   fate-specific retention, and canonical metrics; treat nested subtrees as
   dependent. This P0 diagnostic does not establish where local sensitivity
   concentrates, including the README's ABa observation.
3. Compare replicate-specific null-SD weighting with one fixed spatial scale
   across geometries. Equal-weight results currently include scale differences;
   the same-expression-budget definition avoids arbitrary weight matching.
4. Use expression-constrained integer solves to tighten transfer bounds at
   prespecified budgets, and implement tolerance-constrained feasible consensus.
5. If isolating tracking error is required, obtain localization/retracking
   information or explicitly model correlated coordinate perturbations after
   registration and stage checks. Do not independently jitter cost-matrix entries:
   moving one cell changes many related distances. With three embryos, a fitted
   unrestricted cell-by-cell covariance is unsupported.

## Regeneration and validation

Run from the repository root in the `dev` Conda environment:

```bash
python terminal_pareto/tracking_geometry_sensitivity.py --selftest
python terminal_pareto/tracking_geometry_sensitivity.py
python terminal_pareto/tracking_geometry_sensitivity.py --intervals 1200 \
  --out terminal_pareto/output/tracking_geometry_sensitivity/refined
```

Validation completed:

- Exhaustive enumeration of five small assignment problems verifies both
  constrained-regret bounds; parent-slot symmetry check passes.
- Vectorized spatial/expression construction checked against production builder;
  permutation validity, fixed expression costs, nonnegative optimality regrets,
  zero self-transfer regret, and ordered bounds asserted during both runs.
- Saved assignments replay both objectives; all 301 shared-weight assignments
  are identical in the 301- and 1,201-weight runs.
- Baseline native front costs reproduce existing S3 caches to 1e-10 null-SD units.
- The denser run discovers 180–187 distinct biological assignments versus
  136–147 initially. Maximum retention and nearest-front distance are unchanged;
  canonical natural position changes by less than 1e-7. This verifies those
  summaries, not exhaustive enumeration of assignments.
- The required canonical self-test passes 17/17 checks. That existing self-test
  rewrites its established canonical caches/table wrapper; no versioned
  publication wrapper diff resulted.
- New source compiles; `git diff --check` passes. Review PNG inspected for layout.

`assignments_and_costs.npz` stores actual permutations, matrices and cell identities;
CSV files contain transfer, neighborhood, locality, exclusion-margin, and matching
records. `provenance.json` records input/source hashes, features and numerical
library versions. Generated outputs are ignored by Git and must be preserved
separately if sharing this exploratory result.

## Follow-up: one common assignment can improve on nature in all three embryos

The additional joint optimization in `common_tracking_assignment.py` answers
the existence question affirmatively. It uses the same 275 matched children,
172 distinct candidate parents, and original per-parent capacities. A common
assignment means one biological parent for each child, used in all three
geometries; it is not restricted to edges selected on every individual sweep.

For expression budget \(b\), define the best common travel improvement by

\[
\max_A\min_{r\in\{1,2,3\}}
\frac{T_r(L)-T_r(A)}{T_r(L)}
\quad\text{subject to}\quad S(A)\leq b
\]

Equivalently, minimize an auxiliary variable \(z\) subject to
\(T_r(A)/T_r(L)\leq z\) for each embryo. Assignment variables are binary;
each child receives one parent and each parent retains its original number
of available child slots. Repeated parent rows are collapsed and checked for
exact equality before constructing the capacitated integer problem. This
optimizes the least-improved embryo, using percent of natural travel as the
explicit normalization. Other definitions of best could select another solution.

With expression cost constrained not to exceed nature, the optimum is
\(z=0.9772399253\), giving a guaranteed travel reduction of 2.2760%:

| Embryo | Natural travel | Common-assignment travel | Reduction |
|---|---:|---:|---:|
| 1 | 740.845103 | 723.001090 | 2.4086% |
| 2 | 780.637040 | 762.554955 | 2.3163% |
| 3 | 750.240285 | 733.164760 | 2.2760% |

The same assignment reduces cell-state cost from 128.978448 to 128.945038
(0.02590%) and retains 227/275 natural edges (82.55%; 48 changed parents).
Thus it strictly improves both objectives in each of the three geometries.
The expression improvement is small; this solution principally optimizes travel.
HiGHS returned optimal status and zero reported MIP gap; the objective lower
bound and verified feasible value differ by less than 1e-14. Optimality is
within the solver's numerical tolerances and does not imply uniqueness.

Requiring at least 1% expression reduction also yields a certified optimum:
actual expression reduction is 1.0254%, and travel reductions are 1.5118%,
1.4867%, and 1.4687%. This assignment retains 231/275 natural edges (84%).
Consequently, simultaneous advantage is not limited to the almost unchanged
expression cost of the first example.

Three completed budget solves illustrate the joint trade-off:

| Required expression reduction | Actual expression reduction | Least travel reduction across embryos |
|---|---:|---:|
| 0% | 0.0259% | 2.2760% |
| 0.5% | 0.5529% | 1.8968% |
| 1% | 1.0254% | 1.4687% |

All three report optimal status with zero MIP gap. A subsequent 2% budget solve
ran beyond its internal 45-second time limit and was manually interrupted;
the requested 5% budget was not attempted. Neither has a claimed result.

These common solutions need not lie exactly on any individual embryo's front.
The budget solves sample a joint minimax travel versus expression trade-off;
they do not enumerate the entire common Pareto set or apply a secondary
expression tie-break among equal minimax optima. Their advantage is evaluated
on the three geometries used for fitting, not on a held-out embryo. They show
that failure to transfer a single-embryo travel endpoint does not preclude
finding a more modest travel improvement shared across all three.

Results and actual child-to-parent assignments are stored in
`output/common_tracking_assignment/`. Additional expression budgets are in
its `expression_tradeoff/` subdirectory. The CSVs record solver bounds/status
so an incumbent obtained at a time limit is distinguishable from an optimum.

```bash
python terminal_pareto/common_tracking_assignment.py --selftest
python terminal_pareto/common_tracking_assignment.py --reductions 0 --seconds 45
python terminal_pareto/common_tracking_assignment.py --reductions 0.5 1 \
  --seconds 45 --out terminal_pareto/output/common_tracking_assignment/expression_tradeoff
```

Validation: exhaustive enumeration of three small capacity-constrained problems
matches the integer solver; all three saved common assignments replay all four
objective totals and parent capacities. Their combined verified summary is
`output/common_tracking_assignment/verified_common_tradeoff.csv`.
New source compilation and
`git diff --check` pass. No publication renderers or caches are regenerated
by this follow-up.

## CE-only Figure S3 amendment candidate (2026-09-15)

The author requested a presentable supplement and then chose to keep it
C. elegans only after confirming that the available AF16 set has two embryos.
The AF16 attempt was stopped; no AF16 result is included in this figure.

`output/tracking_s3_amendment/figS3_ce_tracking_amended_review.pdf` is the
two-page review artifact:

- Page 1 retains the global fronts and subtree comparisons from S3 panels
  A--B. The 2026-09-16 revision adds three zooms to candidate panel A, placing
  every common assignment in each embryo's own travel/cell-state space.
  Stars, squares, and triangles identify the common solutions requiring 0%,
  0.5%, and 1% cell-state improvement, respectively. All nine evaluated points
  have both costs below nature, as shown by the shaded lower-left region.
  The caption and standalone figure counter are updated in a generated review
  copy; accepted publication files remain intact.
- Page 2 adds C: biological-parent agreement; D: candidate-parent neighborhood
  overlap; E: directed transfer penalties relative to target weighted optima;
  F: pure-travel transfer penalties relative to natural-to-optimal travel
  savings; G: the three certified common assignments and retained natural edges.

This is an amendment candidate, not an adopted replacement. Retaining A--B
preserves the existing subtree sensitivity evidence. E and F deliberately use
different denominators, stated on the axes and in the caption. The complete
three-target projected-front view remains available in
`output/tracking_geometry_sensitivity/tracking_geometry_review.pdf`.

The renderer `figS3_ce_tracking_variance_amendment.py` uses the validated
301-weight caches and does not rerun optimization. It checks source hashes,
replays all nine source/target cost evaluations, and verifies all three common
assignments' costs, parent capacities, and solver statuses before plotting.
The page-2 caption source is `figS3_ce_tracking_variance_amendment.tex`.
The renderer creates the page-1 review wrapper from the accepted caption,
adding the common-assignment explanation. It writes the nine plotted
coordinates to `common_assignments_in_front_space.csv` and verifies that the
original native front costs agree with the shared optimization cache.

Regenerate from the repository root in the `dev` environment:

```bash
python terminal_pareto/figS3_ce_tracking_variance_amendment.py
cd terminal_pareto
tectonic figS3_ce_tracking_variance_amendment.tex --outdir output/tracking_s3_amendment
cd output/tracking_s3_amendment
tectonic figS3_existing_numbered.tex
pdfunite figS3_existing_numbered.pdf figS3_ce_tracking_variance_amendment.pdf \
  figS3_ce_tracking_amended_review.pdf
pdftoppm -scale-to 1700 -png figS3_ce_tracking_amended_review.pdf review
```

Both pages are rendered and visually checked. Python compilation and
`git diff --check` pass. The canonical 17-check self-test is run with its output
directories redirected to the amendment's `validation/` directory. An initial
attempt to redirect to `/tmp` hit the existing table writer's relative-path
assumption and rewrote the standard canonical NPZ through a captured default
argument; the final check explicitly redirects that writer too. Accepted
figure PDFs and panel assets are unchanged.
Figure inputs and their hashes are recorded in `figure_provenance.json` in
the amendment directory. Generated PDF/PNG/CSV artifacts remain ignored by Git.
The 2026-09-16 change is a rendering-only addition using the same verified
solutions; its cost replay, compilation, and visual checks were repeated.
