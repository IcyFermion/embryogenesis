"""Captions of the revised layout's main figures, as TeX caption bodies (title sentence first).

Numbering of 2026-11-10: the cell-state feature figure is Supplementary Figure S2,
so the main figures are 1-7 again (as in draft 4). ``AGREED`` holds captions the
author has approved: Figure 6 (full tree by layer; agreed 2026-11-05, record:
``full_tree_pareto/LAYERWISE_HOURGLASS.md``). ``DRAFT`` holds captions written in
the style of the v1 baseline, for review:
- Figure 1 takes A-C from the manuscript draft and D from the released
  endpoint-coordinates amendment.
- Figures 2 and 4 adapt the released Figures 2-3A and 5 (4 draws d_LP and d_CP
  on one axis; r stays in Table 1).
- Figure 7 takes its facts from released Figure 8; it is the C. briggsae front
  alone (2026-11-10).
- Figure 3 is new; S2 takes its facts from Track A
  (``terminal_pareto/METHODS_cell_state_representation.md``).
``RELEASED`` holds the released captions of Figure 5 (released Figure 6, cell types)
and S1 (released S3, four tracking geometries), unchanged; they moved here on
2026-11-13, when the set stopped reusing released figures.
``BACKLOG`` holds the caption of the backlog page 7B-C (C. briggsae counterparts
of Figures 3A and 6B, the former 8B-C).
Counts are those of the released analyses; Figure 3's values were checked against the cached global sweep.
"""

from __future__ import annotations

from publication import notation as nt


def figure6() -> str:
    d = nt.tex("natural_front_distance")
    return rf"""\textbf{{Full-tree trade-offs by layer.}}
The tracked \emph{{C.\@ elegans}} lineage starts at the 4-cell stage (ABa, ABp, EMS and P2) and ends at 491
measured leaves (three embryos). Its 974 mother--daughter edges are split into eight layers, built bottom-up from
the leaves: the leaves form layer 8, and each layer above holds the mothers whose daughters have all been placed
below. A cell's layer is therefore set by its distance from the leaves, counted in generations along its longest
branch. Layer 1 holds five daughters of the 4-cell-stage cells (ABal, ABar, ABpl, ABpr and MS); E and P3 sit lower
because their branches are shorter. Within each layer, the daughters are reassigned to their mothers' slots with
every other layer fixed (Methods).
\textbf{{(A)}} Each layer's Pareto front (black) and the natural lineage (red cross), in that layer's endpoint
coordinates: travel optimum at \((0,1)\), cell-state optimum at \((1,0)\). The clouds permute daughters that share a
grandmother (first-cousin shuffle), a great-grandmother (second-cousin) or a great-great-grandmother (third-cousin)
among their mothers' slots. The insets permute all of the layer's daughters (full-random shuffle) on separately
scaled axes. Each cloud shows 600 of 2,000 draws. The fronts of layers 1 and 2 (5 and 10 cells) have only 2 and 6
distinct points, so B and C start at layer 3.
\textbf{{(B, C)}} One point per layer, at the median number of cells alive when the layer's cells are born (cell
table of Large et al.\ 2025); the top axis names the layer. \textbf{{(B)}} Top: distance {d} from the lineage to its
layer's front. Bottom: the lowest percentage of edges that differ from the lineage, across all assignments on the
front (100\% \(-\) maximum lineage similarity). Layer 3's front has only 5 distinct points and contains the lineage,
so both values are 0 there. \textbf{{(C)}} Natural travel and cell-state cost divided by the exact mean cost of the
first-cousin shuffle; 1 means no better than the shuffle. Error bars are 95\% bootstrap percentile intervals from
2,000 resamples of first-cousin groups within each layer (6 groups in layer 3). B has no intervals."""


def _symbols():
    keys = ("travel_optimum", "cell_state_optimum", "natural_lineage", "null_mean", "closest_point",
            "natural_front_distance", "null_front_distance", "canonical_position")
    return {key: nt.tex(key) for key in keys}


def figure1() -> str:
    m = _symbols()
    return rf"""\textbf{{Lineage structure creates a trade-off between cell positioning and fate specification.}}
\textbf{{(A)}} Schematic of a developmental lineage tree. A sequence of cell divisions generates terminal cells with
distinct positions and molecular identities that form different tissues and organs.
\textbf{{(B)}} Two alternative lineage organizations that produce the same terminal cell identities and final
positions. Left: daughters with similar molecular states are generated from the same parent, resulting in a small
cell-state distance but requiring some cells to travel farther to their final positions. Right: distinct cell types
are generated closer to their final positions, reducing travel distance but requiring larger parent-to-daughter
changes in molecular state.
\textbf{{(C)}} The space of all possible lineages, evaluated according to their total travel and cell-state distances.
The Pareto front contains the set of lineages for which neither objective can be improved without worsening the
other. Lineages away from the front are dominated by alternatives that achieve lower costs along both objectives and
would therefore be disfavored if these costs influence developmental performance or organismal fitness.
\textbf{{(D)}} Endpoint coordinates and canonical metrics. The travel optimum {m['travel_optimum']} maps to \((0,1)\) and
the cell-state optimum {m['cell_state_optimum']} to \((1,0)\); travel and cell-state costs are normalized by the spans
between these two attained endpoint assignments. The natural lineage {m['natural_lineage']} and the first-cousin null
mean {m['null_mean']} use the same anchors. {m['closest_point']} is the closest attained sampled-front assignment to
{m['natural_lineage']} in these coordinates; {m['natural_front_distance']} and {m['null_front_distance']} are measured
from {m['natural_lineage']} and {m['null_mean']} to this same assignment. Canonical position {m['canonical_position']}
is the fraction of front arc length from {m['travel_optimum']} to {m['closest_point']} (green arc; the orange arc runs
from {m['closest_point']} to {m['cell_state_optimum']}). Each subtree uses its own endpoints for shape summaries, while
assignment comparisons within a cost space share one reference. This display transformation leaves the optimization
and selected assignments unchanged."""


def figure_s2() -> str:
    return r"""\textbf{Choosing the reporter proteins of the cell-state cost.}
Cell-state distance is the cosine distance between protein-reporter profiles (z-scored intensities) of a parent slot
and a child, over 20 of the 210 reporters. The 20 are chosen by terminal fate alone, never by travel or by the
Pareto front. All panels use the 275-cell \emph{C.\@ elegans} terminal cohort.
\textbf{(A)} The 20 reporters with the largest coefficients in an L2-regularized multinomial logistic regression that
classifies terminal cells into fate classes; bar color gives the class of each reporter's largest coefficient.
Inset: importance of all 210 reporters by rank.
\textbf{(B)} Cosine distance between terminal sisters of the same fine cell type (35 pairs) and of different fine
types (63 pairs; WormBase cell descriptions), with the top-20 and with all 210 reporters. \(d\) is Cohen's \(d\) of
the difference, the sibling separation.
\textbf{(C)} Sibling separation with 95\% bootstrap intervals (2,000 resamples of sister pairs) for the top-20, the
more variable half of the panel (HV pool, 105 reporters) and all 210 reporters, against 200 random 20-reporter sets
drawn from the HV pool (violin; dot, median).
\textbf{(D)} Mean cosine distance between all pairs of cohort cells by lineage path length, the number of divisions
separating two cells through their last common ancestor (sisters 2, first cousins 4; lengths pooled in twos), for
the top-20, all 210 reporters and 200 random 20-reporter sets from the HV pool (band, 5--95\%)."""


def figure2() -> str:
    return r"""\textbf{The pooled terminal lineage lies near the travel--cell-state Pareto front.}
\textbf{(A)} The analysis uses the 275 natural terminal parent--child edges shared by three stage-matched
\emph{C.\@ elegans} tracking embryos. Travel is the equal mean of the three edge-distance matrices after division by
each embryo's global natural travel total; cell-state distance is top-20 protein cosine distance. The black curve
connects 301 weighted assignment solutions; connecting segments need not be attainable assignments. Axes are
fractions of the pooled front's endpoint cost spans, with the travel optimum at \((0,1)\) and the cell-state optimum
at \((1,0)\). The red cross marks the natural lineage, shown at its transformed position rather than forced to the
origin. First-, second- and third-cousin shuffles and the full-random inset use the same pooled endpoint reference;
the inset has its own axis limits.
\textbf{(B)} Nested cousin-shuffling null models progressively relax the lineage neighborhood within which terminal
parent assignments are permuted. Terminal cells (white circles) are permuted among the parent slots of cells that
share a grandparent (first-cousin shuffle), a great-grandparent (second-cousin) or a great-great-grandparent
(third-cousin); each shared ancestor takes the color of its shuffle, and brackets mark each shuffle's range. The
full-random shuffle permutes all terminal cells. Terminal identities, positions and molecular states remain fixed."""


def figure3() -> str:
    m = _symbols()
    return rf"""\textbf{{Lineage similarity and rewiring along the pooled terminal Pareto front.}}
Lineage similarity is the fraction of the 275 natural terminal parent--child edges that an assignment keeps.
\textbf{{(A)}} Lineage similarity along the pooled front of Figure~2A, against travel distance in the same endpoint
coordinates. The red cross marks the natural lineage (similarity 1). Squares mark the travel optimum
{m['travel_optimum']} and the cell-state optimum {m['cell_state_optimum']}; the circle marks the maximum-similarity
assignment M. \textbf{{(B--D)}} The natural lineage (top, founder lineages labelled) and one row each for
{m['travel_optimum']}, M and {m['cell_state_optimum']}. Each arc runs from a terminal cell that changes parent (black
dot, at its place in the natural lineage) to its assigned parent, drawn at the mean position of that parent's natural
daughters. Arc color gives the tree distance between the natural and the assigned parent, in divisions. Lineage
similarity is 0.71 at {m['travel_optimum']}, 0.83 at M and 0.19 at {m['cell_state_optimum']}, where 81, 46 and 223 of the
275 terminal cells change parent."""


def figure4() -> str:
    m = _symbols()
    return rf"""\textbf{{Canonical position and distances to the front for the pooled subtree fronts.}}
Each of the 18 subtrees with at least 24 of the 275 natural terminal edges appears once, the smallest size above
which every subtree's first-cousin null and front are resolved (Methods). Rows are grouped into the root and AB branch
and the ABa, ABp and P1 branches. Left-right counterparts sit in adjacent rows joined by a bracket: ABpl and ABpr,
ABplp and ABprp, MSa and MSp, and ABalp and ABara, which correspond by fate rather than by lineage position. Canonical
metrics and the attained closest assignment {m['closest_point']} are defined in Figure~1D; each subtree uses its own
two single-objective endpoint spans. Pies report the terminal fate composition of each subtree. Gray dots show
canonical position {m['canonical_position']}. On a shared axis, pink dots show the distance {m['natural_front_distance']}
from the natural lineage to {m['closest_point']} and gold dots the distance {m['null_front_distance']} from the
first-cousin null mean {m['null_mean']} to the same assignment; each dot sits on a stem from zero. In every subtree
the natural lineage lies closer to the front than its null mean. Column scales are shared across groups. A zero
distance places the natural lineage on its sampled front. These shape-normalized coordinates are not absolute
performance comparisons between subtrees; nested subtrees are not independent biological replicates.
\(\dagger\)~ABplp lies farthest from its front because of where the track of one apoptotic cell, ABplpppapp, ends: with
that endpoint held at the cell's birth position, ABplp lies on its front, whereas the same change to either of its
other two apoptotic cells, or to three living cells matched on depth and birth time (0 of 1000 draws), leaves it
off."""


def figure7() -> str:
    return r"""\textbf{Terminal trade-offs in \emph{C.\@ briggsae}.}
The natural \emph{C.\@ briggsae} AF16 lineage against its Pareto front, laid out as in Figure~2A. The analysis uses
the 187 natural terminal parent--child edges matched between \emph{C.\@ elegans} and \emph{C.\@ briggsae}. Travel is
the equal mean of per-embryo pairwise displacements in two AF16 embryos (3D tracking), each divided by its fixed
natural-lineage total on the 196-edge \emph{C.\@ briggsae} reference before restriction to the matched cohort.
Cell-state cost is cosine dissimilarity across the 20 shared RNA transcription factors. The black curve connects 301
weighted assignment solutions in the front's own endpoint coordinates, with the travel optimum at \((0,1)\) and the
cell-state optimum at \((1,0)\); the red cross marks the natural lineage, and the cousin shuffles and full-random
inset are defined as in Figure~2. Figure~2A uses all 275 \emph{C.\@ elegans} edges and protein cell state, so the two
fronts compare trade-off shapes, not absolute costs. All \emph{C.\@ briggsae} 3D tracking is subject to the
limitations of traditional embryo-tracking techniques, particularly for z-axis measurements."""


def backlog_figure7bc() -> str:
    """Backlog page (2026-11-10): the former Figure 8B-C, kept in reserve for reviewers."""
    return r"""\textbf{\emph{C.\@ briggsae} counterparts of Figures~3A and 6B (backlog).}
Kept in reserve: the gaps in data resolution between the two species (weaker \emph{C.\@ briggsae} tracking,
especially in z; RNA rather than protein cell state) make these comparisons hard to interpret.
\textbf{(B)} Lineage similarity along the \emph{C.\@ briggsae} front of Figure~7 against travel distance, as in
Figure~3A: the fraction of the 187 natural terminal edges that each assignment keeps. Squares mark the travel and
cell-state optima and the circle the maximum-similarity assignment M (similarity 0.74; 0.50 at the travel optimum and
0.21 at the cell-state optimum).
\textbf{(C)} The \emph{C.\@ briggsae} lineage by layer, as in Figure~6B, for the partial forest of the matched
cells: the 187 matched terminal cells and their ancestors, with layers numbered as in the \emph{C.\@ elegans} full
tree (layer 8 holds the leaves). Each layer's daughters are reassigned to their mothers' slots with the other layers
fixed, using AF16 3D travel and \emph{C.\@ briggsae} RNA cell state. Points sit at the median number of cells alive
when the matched \emph{C.\@ elegans} cells are born (Large et al.\ 2025). Layer 3 holds 2 cells and has no
resolved front, so C starts at layer 4. The distance axis is shared with Figure~6B."""


def figure5() -> str:
    """Released caption of the cell-type figure (released Figure 6), unchanged."""
    return r"""\textbf{Terminal fates differ in retained lineage structure and in
the cost of type restriction under pooled travel.} \textbf{(A)} Fate-specific
natural-edge retention along the pooled unrestricted front, parameterized by
canonical arc length. \textbf{(B)} Separate within-type shape comparisons use
each type's own endpoint spans. All six displayed types have valid endpoint
spans; degenerate small-type cases remain explicitly flagged in the numerical
cache. \textbf{(C)} The type-restricted
and unrestricted aggregate fronts share the unrestricted pooled endpoint
reference, so their performance difference remains visible. Endpoint
penalties are also reported in the interpretable global null-SD units."""


def figure_s1() -> str:
    """Released caption of the tracking-robustness figure (released S3), unchanged."""
    return r"""\textbf{Geometry dependence of terminal assignments and the pooled
travel compromise.} The same 275 terminal edges and top-20 protein cosine
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
improved held-out travel in only 2 of 12 settings."""


AGREED = {"6": figure6}
DRAFT = {"1": figure1, "2": figure2, "3": figure3, "4": figure4, "7": figure7, "S2": figure_s2}
RELEASED = {"5": figure5, "S1": figure_s1}
BACKLOG = {"7B-C": backlog_figure7bc}
