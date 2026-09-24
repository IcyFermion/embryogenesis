"""Write pooled-profile LaTeX wrappers beside generated publication panels."""

from __future__ import annotations

from pathlib import Path
import shutil


PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[margin=0.55in]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{amsmath}
\usepackage{tikz}
\captionsetup{font=small,labelfont=bf}
\pagestyle{empty}
\begin{document}
"""


def _write(out_dir: Path, stem: str, body: str) -> Path:
    path = out_dir / f"{stem}.tex"
    # Standalone wrappers should carry the same labels as the manuscript.
    number = stem.split("_")[0].removeprefix("fig")
    if number == "1":
        numbering = r"\renewcommand{\thefigure}{1 amendment}" + "\n"
    elif number.startswith("S"):
        numbering = rf"\renewcommand{{\thefigure}}{{{number}}}" + "\n"
    else:
        numbering = rf"\setcounter{{figure}}{{{int(number)-1}}}" + "\n"
    path.write_text(PREAMBLE + numbering + body.strip() + "\n\\end{document}\n")
    return path


def write_pooled_wrappers(out_dir: Path | str) -> list[Path]:
    """Generate candidate wrappers without touching accepted legacy files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    null_schematic = (Path(__file__).resolve().parent / "assets"
                      / "fig3A_ce_null_models.tex")
    null_schematic_copy = out_dir / null_schematic.name
    if null_schematic.resolve() != null_schematic_copy.resolve():
        shutil.copy2(null_schematic, null_schematic_copy)
    wrappers = []
    wrappers.append(_write(out_dir, "fig1_ce_endpoint_normalization_amendment", r"""
\begin{figure}[p]
\centering
\includegraphics[width=0.82\textwidth]{fig1_endpoint_normalization_amendment.pdf}
\caption{\textbf{Endpoint coordinates and canonical metrics.} The travel
optimum maps to $(0,1)$ and the cell-state optimum to $(1,0)$. Travel cost
$T$ and cell-state cost $S$ are normalized by the spans between these two
attained endpoint assignments. The natural lineage $L$ and the first-cousin
null mean $N$ use the same anchors. $P^*$ is the closest attained sampled-front
assignment to $L$ in these coordinates; $d_{LP}$ and $d_{NP}$ are measured
from $L$ and $N$ to this same assignment. Canonical position $u$ is the
fraction of front arc length from the travel optimum to $P^*$. Each subtree
uses its own endpoints for shape summaries, while assignment comparisons
within a cost space share one reference. This display transformation leaves
the optimization and selected assignments unchanged. This independent
component awaits integration into the co-author's Figure~1.}
\label{fig:ce_endpoint_normalization_amendment}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "fig2_ce_terminal_pareto_main", r"""
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{fig2_ce_terminal_pareto_front.pdf}
\caption{\textbf{The pooled terminal lineage lies near the travel--cell-state
Pareto front.} The highlighted circle marks the assignment with maximum
natural-edge retention; it need not be the closest assignment $P^*$ defined
in the separate Figure~1 amendment.
The analysis uses the 275 natural terminal parent--child edges
shared by three stage-matched \emph{C.\@ elegans} tracking embryos. Travel is
the equal mean of the three edge-distance matrices after division by each
embryo's global natural travel total; cell-state distance is top-20 protein
cosine distance. Axes are fractions of the pooled front's endpoint cost spans,
with the travel optimum at $(0,1)$ and the cell-state optimum at $(1,0)$.
The natural lineage is shown at its transformed position rather than forced to
the origin. First-, second-, and third-cousin shuffles and the full-random
inset use the same pooled endpoint reference. Blue encodes natural-edge
retention.}
\label{fig:ce_terminal_pareto_main}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "fig3_ce_terminal_pareto_supporting", r"""
\begin{figure}[p]
\centering
\resizebox{\textwidth}{!}{\input{fig3A_ce_null_models.tex}}\\[2mm]
\begin{minipage}{0.49\textwidth}\centering
\includegraphics[width=\linewidth]{fig3B_ce_edge_retention_tree_distance.pdf}
\end{minipage}\hfill
\begin{minipage}{0.49\textwidth}\centering
\includegraphics[width=\linewidth]{fig3C_ce_structural_retention.pdf}
\end{minipage}
\caption{\textbf{Null models and structural support along the pooled
terminal-cell Pareto front.} \textbf{(A)} Nested cousin-shuffling null models
progressively relax the lineage neighborhood within which terminal parent
assignments are permuted, ending with the full-random assignment null.
Terminal identities, positions, and molecular states remain fixed.
\textbf{(B)} Natural-edge retention and mean lineage-tree distance
among changed edges are evaluated for the same pooled assignments as
Figure~2; the horizontal cost coordinate uses the pooled P0 endpoint
reference. \textbf{(C)} From the maximum-retention compromise, the final
portion of attainable improvement toward either single-objective endpoint
requires disproportionate loss of natural parentage. Retention and tree
distance keep their original definitions.}
\label{fig:ce_terminal_pareto_supporting}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "fig4_ce_subtree_map", r"""
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{fig4_ce_subtree_map_panel.pdf}
\caption{\textbf{Pooled terminal Pareto behavior across the reference
lineage.} The 42 displayed subtrees contain at least 12 of the 275 strictly
matched natural terminal edges. Pies report the fate composition of the exact
profile cohort; circle area scales with cell count. Solid outlines mark the 27
subtrees in which the natural assignment is an attained sampled-front point,
and dashed outlines mark the remaining subtrees. Label color reports maximum
natural-edge retention.}
\label{fig:ce_subtree_map}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "fig5_ce_canonical_summary", r"""
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{fig5A_ce_canonical_major_subtrees.pdf}\\[2mm]
\includegraphics[width=\textwidth]{fig5B_ce_canonical_all_subtrees.pdf}
\caption{\textbf{Canonical position and endpoint-normalized proximity for the
pooled subtree fronts.} Canonical metrics and the attained closest assignment
$P^*$ are defined in the separate Figure~1 amendment. Each subtree uses its own
two single-objective endpoint spans. \textbf{(A)} For five major subtrees,
purple marks the natural-lineage distance $d_{LP}$ and gold the distance
$d_{NP}$ from the first-cousin-null mean to the same $P^*$. Gray segments pair
the two distances for each named subtree; they are not uncertainty intervals.
\textbf{(B)} Every eligible subtree appears once, grouped into root/AB and
the ABa, ABp, and P1 branches. Gray dots show canonical position $u$; purple
dots show natural-lineage distance $d_{LP}$. Column scales are shared across
groups. Bold names identify the five subtrees in A. A zero distance places
natural lineage on its sampled front. These shape-normalized coordinates
are not absolute performance comparisons between subtrees; nested subtrees
are not independent biological replicates.}
\label{fig:ce_canonical_summary}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "fig6_ce_cell_types", r"""
\begin{figure}[p]
\centering
\includegraphics[width=0.95\textwidth]{fig6A_ce_retention_heatmap.pdf}\\[1mm]
\includegraphics[width=0.95\textwidth]{fig6B_ce_within_type_fronts.pdf}\\[1mm]
\includegraphics[width=0.95\textwidth]{fig6C_ce_type_restricted_aggregate.pdf}
\caption{\textbf{Terminal fates differ in retained lineage structure and in
the cost of type restriction under pooled travel.} \textbf{(A)} Fate-specific
natural-edge retention along the pooled unrestricted front, parameterized by
canonical arc length. \textbf{(B)} Separate within-type shape comparisons use
each type's own endpoint spans. All six displayed types have valid endpoint
spans; degenerate small-type cases remain explicitly flagged in the numerical
cache. \textbf{(C)} The type-restricted
and unrestricted aggregate fronts share the unrestricted pooled endpoint
reference, so their performance difference remains visible. Endpoint
penalties are also reported in the interpretable global null-SD units.}
\label{fig:ce_cell_types}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "figS1_ce_canonical_summary_cousin_r", r"""
\begin{figure}[p]
\centering
\begin{minipage}{0.39\textwidth}\centering
\includegraphics[width=\linewidth]{figS1A_ce_cousin_r_definition.pdf}
\end{minipage}\hfill
\begin{minipage}{0.59\textwidth}\centering
\includegraphics[width=\linewidth]{figS1B_ce_cousin_r_major_subtrees.pdf}
\end{minipage}\\[1mm]
\includegraphics[width=\textwidth]{figS1C_ce_cousin_r_all_subtrees.pdf}
\caption{\textbf{First-cousin-null-relative sensitivity for the pooled
subtree analysis.} The statistic $r$ retains its separately null-standardized
definition and is not used as an endpoint scale. The panels report the five
major subtrees and all 42 eligible subtrees.}
\label{fig:ce_canonical_summary_cousin_r}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "figS2_ce_cell_type_cost_gain", r"""
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figS2_ce_cell_type_cost_gain_panel.pdf}
\caption{\textbf{Fate-specific cost changes at three pooled global-front
assignments.} Per-cell travel and cell-state changes are reported separately
in pooled embryo-wide first-cousin-null standard-deviation units relative to
the natural lineage. The travel optimum, maximum-retention compromise, and
cell-state optimum are exact assignments; the two objectives are not combined
into a biological cost.}
\label{fig:ce_cell_type_cost_gain}
\end{figure}
"""))
    wrappers.append(_write(out_dir, "figS3_ce_tracking_robustness", r"""
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figS3_four_geometries_with_insets.pdf}
\caption{\textbf{Geometry dependence of terminal assignments and the pooled
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
improved held-out travel in only 2 of 12 settings.}
\label{fig:ce_tracking_robustness}
\end{figure}
"""))
    return wrappers
