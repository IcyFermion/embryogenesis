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
