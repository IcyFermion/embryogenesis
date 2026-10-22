"""Captions for the terminal publication wrappers (Figures 1 amendment-6, S1-S3).

Moved from ``terminal_pareto/publication_wrappers.py``. Each entry returns the
complete ``figure`` environment; mathematical symbols come from the notation
registry (``$...$`` delimiters, as in the published wrappers).
"""

from pathlib import Path

from publication import notation as nt

# Versioned TikZ schematic for Figure 3A, \input by the Figure 3 wrapper.
NULL_SCHEMATIC = Path(__file__).resolve().parents[1] / "assets" / "fig3A_ce_null_models.tex"

SYMBOLS = ['canonical_position', 'cell_state_cost', 'closest_point', 'cousin_relative', 'natural_front_distance', 'natural_front_distance_sd', 'natural_lineage', 'natural_state_position', 'natural_travel_position', 'null_front_distance', 'null_mean', 'travel_cost']


def _m():
    return {key: nt.math(key) for key in SYMBOLS}


def fig1_ce_endpoint_normalization_amendment():
    M = _m()
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=0.82\textwidth]{{fig1_endpoint_normalization_amendment.pdf}}
\caption{{\textbf{{Endpoint coordinates and canonical metrics.}} The travel
optimum maps to $(0,1)$ and the cell-state optimum to $(1,0)$. Travel cost
{M['travel_cost']} and cell-state cost {M['cell_state_cost']} are normalized by the spans between these two
attained endpoint assignments. The natural lineage {M['natural_lineage']} and the first-cousin
null mean {M['null_mean']} use the same anchors. {M['closest_point']} is the closest attained sampled-front
assignment to {M['natural_lineage']} in these coordinates; {M['natural_front_distance']} and {M['null_front_distance']} are measured
from {M['natural_lineage']} and {M['null_mean']} to this same assignment. Canonical position {M['canonical_position']} is the
fraction of front arc length from the travel optimum to {M['closest_point']}. Each subtree
uses its own endpoints for shape summaries, while assignment comparisons
within a cost space share one reference. This display transformation leaves
the optimization and selected assignments unchanged. This independent
component awaits integration into the co-author's Figure~1.}}
\label{{fig:ce_endpoint_normalization_amendment}}
\end{{figure}}
"""


def fig2_ce_terminal_pareto_main():
    M = _m()
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=\textwidth]{{fig2_ce_terminal_pareto_front.pdf}}
\caption{{\textbf{{The pooled terminal lineage lies near the travel--cell-state
Pareto front.}} The highlighted circle marks the assignment with maximum
natural-edge retention; it need not be the closest assignment {M['closest_point']} defined
in the separate Figure~1 amendment.
The analysis uses the 275 natural terminal parent--child edges
shared by three stage-matched \emph{{C.\@ elegans}} tracking embryos. Travel is
the equal mean of the three edge-distance matrices after division by each
embryo's global natural travel total; cell-state distance is top-20 protein
cosine distance. Axes are fractions of the pooled front's endpoint cost spans,
with the travel optimum at $(0,1)$ and the cell-state optimum at $(1,0)$.
The natural lineage is shown at its transformed position rather than forced to
the origin. First-, second-, and third-cousin shuffles and the full-random
inset use the same pooled endpoint reference. Blue encodes natural-edge
retention.}}
\label{{fig:ce_terminal_pareto_main}}
\end{{figure}}
"""


def fig3_ce_terminal_pareto_supporting():
    return rf"""
\begin{{figure}}[p]
\centering
\resizebox{{\textwidth}}{{!}}{{\input{{fig3A_ce_null_models.tex}}}}\\[2mm]
\begin{{minipage}}{{0.49\textwidth}}\centering
\includegraphics[width=\linewidth]{{fig3B_ce_edge_retention_tree_distance.pdf}}
\end{{minipage}}\hfill
\begin{{minipage}}{{0.49\textwidth}}\centering
\includegraphics[width=\linewidth]{{fig3C_ce_structural_retention.pdf}}
\end{{minipage}}
\caption{{\textbf{{Null models and structural support along the pooled
terminal-cell Pareto front.}} \textbf{{(A)}} Nested cousin-shuffling null models
progressively relax the lineage neighborhood within which terminal parent
assignments are permuted, ending with the full-random assignment null.
Terminal identities, positions, and molecular states remain fixed.
\textbf{{(B)}} Natural-edge retention and mean lineage-tree distance
among changed edges are evaluated for the same pooled assignments as
Figure~2; the horizontal cost coordinate uses the pooled P0 endpoint
reference. \textbf{{(C)}} From the maximum-retention compromise, the final
portion of attainable improvement toward either single-objective endpoint
requires disproportionate loss of natural parentage. Retention and tree
distance keep their original definitions.}}
\label{{fig:ce_terminal_pareto_supporting}}
\end{{figure}}
"""


def fig4_ce_subtree_map():
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=\textwidth]{{fig4_ce_subtree_map_panel.pdf}}
\caption{{\textbf{{Pooled terminal Pareto behavior across the reference
lineage.}} The 42 displayed subtrees contain at least 12 of the 275 strictly
matched natural terminal edges. Pies report the fate composition of the exact
profile cohort; circle area scales with cell count. Solid outlines mark the 27
subtrees in which the natural assignment is an attained sampled-front point,
and dashed outlines mark the remaining subtrees. Label color reports maximum
natural-edge retention.}}
\label{{fig:ce_subtree_map}}
\end{{figure}}
"""


def fig5_ce_canonical_summary():
    M = _m()
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=\textwidth]{{fig5A_ce_canonical_major_subtrees.pdf}}\\[2mm]
\includegraphics[width=\textwidth]{{fig5B_ce_canonical_all_subtrees.pdf}}\\[2mm]
\includegraphics[width=\textwidth]{{fig5C_ce_mirrored_pairs.pdf}}
\end{{figure}}
\begin{{figure}}[p]
\caption{{\textbf{{Canonical position and endpoint-normalized proximity for the
pooled subtree fronts.}} Canonical metrics and the attained closest assignment
{M['closest_point']} are defined in the separate Figure~1 amendment. Each subtree uses its own
two single-objective endpoint spans. \textbf{{(A)}} For five major subtrees,
purple marks the natural-lineage distance {M['natural_front_distance']} and gold the distance
{M['null_front_distance']} from the first-cousin-null mean to the same {M['closest_point']}. Gray segments pair
the two distances for each named subtree; they are not uncertainty intervals.
\textbf{{(B)}} Every eligible subtree appears once, grouped into root/AB and
the ABa, ABp, and P1 branches. Gray dots show canonical position {M['canonical_position']}; purple
dots show natural-lineage distance {M['natural_front_distance']}. Column scales are shared across
groups. Bold names identify the five subtrees in A. A zero distance places
natural lineage on its sampled front. These shape-normalized coordinates
are not absolute performance comparisons between subtrees; nested subtrees
are not independent biological replicates. Numbers after a name mark the
mirrored pairs of C. \textbf{{(C)}} Left (blue) and right (orange) partners
of the seven left--right subtree pairs in which each subtree is the other's
best fate match (sided match $\geq 0.85$, at least 12 edges per side), and of
two exploratory anterior pairs (hollow). Besides {M['canonical_position']} and
{M['natural_front_distance']}, the columns give three readouts in units of the
first-cousin null SD: the distance {M['natural_front_distance_sd']} from {M['natural_lineage']} to
the nearest attained front point (in this frame not necessarily
{M['closest_point']}), and the travel and cell-state positions
{M['natural_travel_position']} and {M['natural_state_position']} of {M['natural_lineage']} relative
to the null mean. Shading marks a within-pair difference at or above the 90th
percentile of random subtree pairs matched for size and depth. MSa--MSp and
MSap--MSpp agree throughout; the other differences arise in the data or the
metric rather than in development. $\dagger$~ABplpa's front spans 0.04 null SD
of travel, all of it from its one apoptotic cell (without that cell the front
is a single point), so {M['canonical_position']} is not informative; both
partners lie on their fronts. $\ddagger$~MSpa has fewer cohort edges than MSaa
(12 against 16) because late-born cells divide after the tracking cutoff, and
its first-cousin null rests on a single cell group. $\S$~ABplp, ABplpp and
ABplppp lie off their fronts because of where the tracks of their apoptotic
cells end: with those endpoints held at the cells' birth positions all three
lie on their fronts. ABplppp contains one apoptotic cell, ABplpppapp; in
ABplp, the same change to three matched living cells never places the lineage
on its front (0 of 1000 draws).}}
\label{{fig:ce_canonical_summary}}
\end{{figure}}
"""


def fig6_ce_cell_types():
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=0.95\textwidth]{{fig6A_ce_retention_heatmap.pdf}}\\[1mm]
\includegraphics[width=0.95\textwidth]{{fig6B_ce_within_type_fronts.pdf}}\\[1mm]
\includegraphics[width=0.95\textwidth]{{fig6C_ce_type_restricted_aggregate.pdf}}
\caption{{\textbf{{Terminal fates differ in retained lineage structure and in
the cost of type restriction under pooled travel.}} \textbf{{(A)}} Fate-specific
natural-edge retention along the pooled unrestricted front, parameterized by
canonical arc length. \textbf{{(B)}} Separate within-type shape comparisons use
each type's own endpoint spans. All six displayed types have valid endpoint
spans; degenerate small-type cases remain explicitly flagged in the numerical
cache. \textbf{{(C)}} The type-restricted
and unrestricted aggregate fronts share the unrestricted pooled endpoint
reference, so their performance difference remains visible. Endpoint
penalties are also reported in the interpretable global null-SD units.}}
\label{{fig:ce_cell_types}}
\end{{figure}}
"""


def figS1_ce_canonical_summary_cousin_r():
    M = _m()
    return rf"""
\begin{{figure}}[p]
\centering
\begin{{minipage}}{{0.39\textwidth}}\centering
\includegraphics[width=\linewidth]{{figS1A_ce_cousin_r_definition.pdf}}
\end{{minipage}}\hfill
\begin{{minipage}}{{0.59\textwidth}}\centering
\includegraphics[width=\linewidth]{{figS1B_ce_cousin_r_major_subtrees.pdf}}
\end{{minipage}}\\[1mm]
\includegraphics[width=\textwidth]{{figS1C_ce_cousin_r_all_subtrees.pdf}}
\caption{{\textbf{{First-cousin-null-relative sensitivity for the pooled
subtree analysis.}} The statistic {M['cousin_relative']} retains its separately null-standardized
definition and is not used as an endpoint scale. The panels report the five
major subtrees and all 42 eligible subtrees.}}
\label{{fig:ce_canonical_summary_cousin_r}}
\end{{figure}}
"""


def figS2_ce_cell_type_cost_gain():
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=\textwidth]{{figS2_ce_cell_type_cost_gain_panel.pdf}}
\caption{{\textbf{{Fate-specific cost changes at three pooled global-front
assignments.}} Per-cell travel and cell-state changes are reported separately
in pooled embryo-wide first-cousin-null standard-deviation units relative to
the natural lineage. The travel optimum, maximum-retention compromise, and
cell-state optimum are exact assignments; the two objectives are not combined
into a biological cost.}}
\label{{fig:ce_cell_type_cost_gain}}
\end{{figure}}
"""


def figS3_ce_tracking_robustness():
    return rf"""
\begin{{figure}}[p]
\centering
\includegraphics[width=\textwidth]{{figS3_four_geometries_with_insets.pdf}}
\caption{{\textbf{{Geometry dependence of terminal assignments and the pooled
travel compromise.}} The same 275 terminal edges and top-20 protein cosine
costs are used throughout. Colors identify assignment sets optimized using
each of three embryos separately or pooled travel. Travel is measured in
embryo~1 (A), embryo~2 (B), embryo~3 (C), or the pooled geometry (D), with
assignments held fixed during transfer. Pooled travel is the equal mean of
embryo-specific travel costs divided by their respective natural-lineage
totals. Solid curves are optimized for the displayed geometry; dashed curves
show the other three sources. Each curve connects 301 weighted-sweep
solutions; connecting segments guide the eye and need not represent
attainable assignments. Each panel uses its native front's endpoint anchors
for all four curves and the natural lineage (X). Main axes show the full
sampled fronts; insets enlarge the marked near-natural regions in the same
coordinates. The pooled assignments provide a compromise between the
distinct native fronts of the three measured embryos. Variation in biology
and tracking is not separated here, and these embryos do not establish
performance in unmeasured embryos; the earlier leave-one-geometry-out check
improved held-out travel in only 2 of 12 settings.}}
\label{{fig:ce_tracking_robustness}}
\end{{figure}}
"""


# Panels A-C are taller than a page, so the caption follows on its own page
# until panel A moves to Figure 1 (author decision 2026-10-22).
PAGES = {"fig5_ce_canonical_summary": 2}

# Wrapper stem -> (printed figure number, caption body).
CAPTIONED_FIGURES = {
    "fig1_ce_endpoint_normalization_amendment": ("1 amendment", fig1_ce_endpoint_normalization_amendment),
    "fig2_ce_terminal_pareto_main": ("2", fig2_ce_terminal_pareto_main),
    "fig3_ce_terminal_pareto_supporting": ("3", fig3_ce_terminal_pareto_supporting),
    "fig4_ce_subtree_map": ("4", fig4_ce_subtree_map),
    "fig5_ce_canonical_summary": ("5", fig5_ce_canonical_summary),
    "fig6_ce_cell_types": ("6", fig6_ce_cell_types),
    "figS1_ce_canonical_summary_cousin_r": ("S1", figS1_ce_canonical_summary_cousin_r),
    "figS2_ce_cell_type_cost_gain": ("S2", figS2_ce_cell_type_cost_gain),
    "figS3_ce_tracking_robustness": ("S3", figS3_ce_tracking_robustness),
}
