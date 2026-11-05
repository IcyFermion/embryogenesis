"""Captions agreed for the revised layout, as TeX caption bodies (title sentence first).

Only agreed captions are here; the other main figures are "to follow". Figure 6
was agreed with the author on 2026-11-05 (record and rationale:
``full_tree_pareto/LAYERWISE_HOURGLASS.md``). Its counts are those of the
released full-tree cohort and the Track E runs.
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


AGREED = {"6": figure6}
