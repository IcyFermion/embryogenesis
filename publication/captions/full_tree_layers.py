"""Caption for the per-layer supplement S5 (Track E).

Figure 6B-C of the revised layout is captioned with the rest of that layout.

Counts and ranges come from the validated input, never from the figure.
"""

from __future__ import annotations


def _range(values, digits=1):
    low, high = (f"{v:.{digits}f}" for v in values)
    return low if low == high else f"{low}--{high}"


def travel_decomposition(meta: dict) -> str:
    """Supplementary Figure S5."""
    return rf"""\textbf{{The mid-layer rise in relative travel cost comes mostly from the lineage's own travel.}}
\textbf{{(A)}} Travel per edge in the \emph{{C.\@ elegans}} full tree, from each daughter to its own mother
(natural lineage), to an aunt (mean of the first-cousin shuffle) and to a random mother in the same layer
(mean of the full shuffle). Travel to a random mother is flat across layers
({_range(meta['random'])}\(\times10^{{-3}}\) normalized units per edge), so per-edge costs compare across
layers. \textbf{{(B)}} Percent change in per-edge travel from each flank (layer 8, and the mean of layers
4 and 5) to layers 6 and 7, for the natural lineage and the first-cousin shuffle. Bars: \emph{{C.\@ elegans}} full tree. Triangles:
\emph{{C.\@ briggsae}} partial forest (two AF16 embryos, RNA cell state, raw 3D positions; layers
numbered the same way). Against layer 8, and at layer 6 against layers 4--5, the shuffle cost does not
fall and the natural cost rises more; only layer 7 against layers 4--5 is partly cousins sitting closer
together. Descriptive; no intervals."""
