"""Captions for the cross-species comparison figures (terminal 8/9, partial-forest 10/11).

Every mathematical symbol is requested from ``publication.notation``. Scope
statements (cohort, distance convention, analytic null mean) are explicit per
family because the scientific definitions differ.
"""

from __future__ import annotations

from publication import notation as nt

TRACKING_CAVEAT = (
    "All C. briggsae 3D panels are subject to limitations of traditional embryo-tracking techniques, "
    "particularly for z-axis measurements."
)
TRACKING_CAPTION_STATUS = (
    "This description is provisional and will be refined with experimental collaborators."
)


def _caveat():
    return TRACKING_CAVEAT.replace("C. briggsae", r"\emph{C.\@ briggsae}")


def _symbols():
    t = {key: nt.tex(key) for key in ("closest_point", "natural_front_distance", "null_front_distance",
                                      "null_mean", "natural_lineage", "canonical_position")}
    s = {key: nt.symbol(key) for key in ("closest_point", "natural_front_distance", "null_front_distance",
                                         "null_mean", "natural_lineage")}
    t["lp_definition"] = (rf"\({s['natural_front_distance']}=\lVert D({s['natural_lineage']})"
                          rf"-D({s['closest_point']})\rVert_2\)")
    t["cp_definition"] = (rf"\({s['null_front_distance']}=\lVert D({s['null_mean']})"
                          rf"-D({s['closest_point']})\rVert_2\)")
    return t


def terminal_overlay(meta: dict) -> str:
    """Terminal Figure 9."""
    n, S = meta["edges"], _symbols()
    return rf"""\textbf{{Endpoint-normalized terminal fronts and canonical metrics.}}
The matched {n}-edge cohort, pooled travel, molecular representations and
sampled assignment sweeps are the same as in Figure~8.
\textbf{{(A)}} Three-dimensional distances. \textbf{{(B)}} XY distances.
Blue solid lines show \emph{{C.\@ elegans}} protein, orange dashed lines
\emph{{C.\@ elegans}} RNA, and green dash-dotted lines \emph{{C.\@ briggsae}}
AF16 RNA. Colored crosses mark natural lineage; outlined circles mark
sampled maximum biological-parent retention. Hollow diamonds mark the closest
attained front assignment {S['closest_point']}; dashed connectors run from natural lineage
to {S['closest_point']}, not to maximum retention. Insets enlarge the region near
the natural lineages using the same coordinates as the main axes.
Both panels share main-axis and inset limits, colors and line styles.
Each dataset/geometry retains its own endpoint spans, with travel and
cell-state optima at \((0,1)\) and \((1,0)\), respectively.
\textbf{{(C)}} Named rows group the three configurations under 3D and XY,
following the separate-metric-column layout of Figure~5B. Gray shows canonical
position {S['canonical_position']}, purple natural-lineage distance {S['natural_front_distance']}, and gold distance
{S['null_front_distance']} from the analytic first-cousin null mean {S['null_mean']} to the same {S['closest_point']}.
Here {S['closest_point']} minimizes Euclidean distance to natural lineage {S['natural_lineage']} among the
attained sampled assignments in endpoint coordinates;
{S['lp_definition']} and
{S['cp_definition']}.
Position {S['canonical_position']} is front arc length from the travel optimum to {S['closest_point']}, divided
by total front arc length. Scales are shared across all six rows within each
metric; neither the closest point nor the distance is defined by maximum
retention. Overlays and metrics compare relative trade-off shapes and
natural-lineage proximity, not absolute travel or molecular costs.
Curves connect 301 attained weighted
solutions rather than enumerate the complete discrete Pareto set; connecting
segments need not be attainable assignments. Null clouds are omitted here
and shown in Figure~8. Protein/RNA differences are not a controlled modality
effect, and tracking replicates share molecular measurements.
{_caveat()} {TRACKING_CAPTION_STATUS}"""


def full_tree_overlay(meta: dict) -> str:
    """Partial-forest Figure 11."""
    n, nodes, S = meta["edges"], meta["nodes"], _symbols()
    return rf"""\textbf{{Partial-forest front overlays and canonical metrics.}}
The same {n}-edge, {nodes}-cell terminal-anchored cohort, pooled travel and
layerwise assignment spaces are used as in Figure~10.
\textbf{{(A)}} 3D tracking. \textbf{{(B)}} XY tracking. Blue solid lines show
\emph{{C.\@ elegans}} protein, orange dashed lines \emph{{C.\@ elegans}} RNA,
and green dash-dotted lines \emph{{C.\@ briggsae}} AF16 RNA. Crosses mark natural
lineage and outlined circles sampled maximum biological-parent retention.
Hollow diamonds mark the closest attained assignment {S['closest_point']}, and dashed
connectors join natural lineage to {S['closest_point']}, not maximum retention. Insets
enlarge the natural-lineage region without changing coordinates. Both panels
share main-axis and zoom limits; each configuration retains its own aggregate
layerwise endpoint spans. \textbf{{(C)}} Named rows group configurations under
3D and XY, using separate metric columns as in Figures~5B and 9C.
Gray shows position {S['canonical_position']}, purple distance {S['natural_front_distance']} from natural lineage
{S['natural_lineage']} to {S['closest_point']}, and gold distance {S['null_front_distance']} from the analytic partial-forest
first-cousin-null mean {S['null_mean']} to the same {S['closest_point']}.
{S['closest_point']} minimizes Euclidean distance to {S['natural_lineage']} among sampled attained
assignments in endpoint coordinates; {S['canonical_position']} is front arc length from the travel
optimum to {S['closest_point']}, divided by total arc length. Metric scales are shared
across all six rows within each column. Curves connect {meta['weights']} attained
weighted solutions, not the complete discrete Pareto set; connecting segments
need not be attainable. Separate {meta['dense_weights']:,}-weight
checks are retained in the numerical cache. These panels compare relative
trade-off shapes and lineage proximity, not absolute travel or molecular costs.
Coverage determines the partial-tree boundary; developmental alignment and
RNA measurement provenance remain unresolved. Protein/RNA differences are not
a controlled modality effect, and tracking replicates share molecular matrices.
{_caveat()} {TRACKING_CAPTION_STATUS}"""


def terminal_comparison(meta: dict) -> str:
    """Terminal Figure 8."""
    n, cb_reference = meta["edges"], meta["cb_reference_size"]
    return rf"""\textbf{{Terminal-only trade-offs across molecular representations and species.}}
The same {n} natural terminal parent--child edges and biological-parent
capacities are used throughout. Columns show \emph{{C.\@ elegans}} protein,
\emph{{C.\@ elegans}} RNA and \emph{{C.\@ briggsae}} AF16 RNA.
\textbf{{(A--C)}} Three-dimensional tracking distances.
\textbf{{(D--F)}} Two-dimensional XY distances, with z omitted and travel,
nulls, optima and endpoint coordinates recomputed.
Travel is the equal mean of per-embryo pairwise displacements divided by
fixed natural-lineage totals: three \emph{{C.\@ elegans}} embryos on the
275-edge reference and two AF16 embryos on their {cb_reference}-edge reference,
before restriction to the matched cohort. Displacement is not integrated
trajectory length. Cell-state cost is cosine dissimilarity across the frozen
20 z-scored protein reporters or the shared 20 RNA transcription factors
on stored values. Each sampled front connects 301 weighted assignment
solutions; connecting segments need not be attainable assignments.
Each panel uses its own endpoint cost spans, mapping the travel optimum to
\((0,1)\) and the cell-state optimum to \((1,0)\); natural lineage and nulls
share those anchors. All panels share axis limits and the retention scale.
Blue encodes biological-parent edge retention, black crosses mark natural
lineage, and outlined circles mark sampled maximum retention. Duplicate slots
of one parent do not count as different edges. First-, second- and third-cousin
shuffles are shown in the main axes; full-random shuffles use insets with the
same coordinate transform and separate limits. Plus signs mark the analytic
first-cousin mean or display-sample means for the other nulls.
These are descriptive comparisons, not calibrated significance tests or
absolute cost comparisons between datasets. {_caveat()} {TRACKING_CAPTION_STATUS}"""


def full_tree_comparison(meta: dict) -> str:
    """Partial-forest Figure 10."""
    n, nodes, terminals, roots = meta["edges"], meta["nodes"], meta["terminals"], meta["roots"]
    rounds = ", ".join(map(str, meta["round_edges"]))
    weights, draws = meta["weights"], meta["draws"]
    return rf"""\textbf{{Layerwise trade-offs on a matched terminal-anchored partial forest.}}
Columns show \emph{{C.\@ elegans}} protein, \emph{{C.\@ elegans}} RNA and
\emph{{C.\@ briggsae}} AF16 RNA; \textbf{{(A--C)}} use 3D tracking and
\textbf{{(D--F)}} use XY, recomputing distances, nulls and assignments after omitting z.
Starting at the same {terminals} matched biological terminal edges as Figures~8--9,
each branch ascends through observed canonical ancestors and stops at its first
measurement gap. The resulting {nodes} measured cells and {n} edges form a
partial forest with {roots} fixed boundary roots, not a complete embryonic tree.
No missing state is imputed and no ancestor is skipped. Six asynchronous
contraction rounds ({rounds} edges) preserve observed parent-slot multiplicities,
including one-child boundaries. Exact linear assignments at {weights} shared
weights use global analytic first-cousin-null SD scaling; aggregate costs sum
all rounds. This optimizes a product of round-wise assignment spaces, not
unrestricted tree reconstruction. Travel averages pairwise displacements from
three CE embryos (cutoffs 255/247/225) or two AF16 embryos (148/156), each divided
by its fixed natural total on these {n} edges. Coordinates are not averaged.
Cell-state distance is Euclidean across 20 z-scored protein reporters or the
shared 20 RNA TFs on stored values, unlike terminal cosine distance.
Each panel maps its own travel optimum to \((0,1)\) and state optimum to \((1,0)\);
lineage and nulls share those anchors, with common main-axis limits.
Blue indicates biological-parent retention; black crosses mark natural lineage
and outlined circles sampled maximum retention. Cousin shuffles permute only
the {terminals} terminal identities within canonical ancestor groups two, three
or four transitions back, keeping internal states fixed and scoring all {n}
edges. Random rebuild fixes roots and observed capacities and appears in insets
with the same transform but separate limits. Dots show 1,000 of {draws:,} draws;
plus signs show the analytic first-cousin mean or full-draw means for other
references. These descriptive comparisons are not calibrated significance
tests or absolute-cost species rankings. {_caveat()} {TRACKING_CAPTION_STATUS}"""
