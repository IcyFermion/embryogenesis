"""Captions and wrapper layout for pooled full-tree Figure 7 and supplement S4.

Moved from ``full_tree_pareto/publication_build.py``. Counts are explicit
caption metadata from the validated cohort; the supplement keeps its
constraint inventory table inside the figure environment.
"""

PREAMBLE = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=12mm]{geometry}
\usepackage{graphicx,caption,booktabs}
\pagestyle{empty}
\captionsetup{font=footnotesize,labelfont=bf}
\begin{document}
"""


FIGURE7_NUMBERING = "\\setcounter{figure}{6}\n"
SUPPLEMENT_NUMBERING = "\\renewcommand{\\thefigure}{S\\arabic{figure}}\n\\setcounter{figure}{3}\n"
NUMBER_WORDS = dict(enumerate(("Zero", "One", "Two", "Three", "Four", "Five", "Six",
                               "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve")))


def _fill(template, values):
    for key, value in values:
        template = template.replace(key, str(value))
    return template


def counts(meta):
    rounds = meta["round_edges"]
    return dict(NROUNDS=NUMBER_WORDS.get(len(rounds), str(len(rounds))), NWEIGHTS=f"{meta['weights']:,}",
                NSHOWN=f"{min(meta['display_draws'], meta['draws']):,}", NDRAWS=f"{meta['draws']:,}")


def figure7(meta):
    c = counts(meta)
    return _fill(r"""\begin{figure}[p]
\centering
\includegraphics[width=0.95\textwidth]{fig7A_ce_full_tree_layerwise_rounds.pdf}
\vspace{1mm}
\includegraphics[width=0.95\textwidth]{fig7B_ce_full_tree_collective.pdf}
\caption{\textbf{Full-tree trade-offs under layerwise assignment.}
The shared \emph{C. elegans} protein/tracking cohort contains INTERNAL internal
cells, LEAVES measured leaves and four roots, with EDGES scored edges.
Travel is the equal mean of pairwise displacements in three embryos (cutoffs
255, 247 and 225), each divided by its natural full-cohort travel total;
coordinates are not averaged. Cell-state distance is Euclidean distance
across the selected 20 z-scored protein features.
(A) NROUNDS asynchronous bottom-up contraction rounds (ROUNDS edges)
each reassign children to their natural parent slots using exact linear
assignment across NWEIGHTS shared weights after global first-cousin-null SD scaling.
Color indicates biological-parent
retention; outlined circles mark maximum-retention sampled assignments.
Rounds are not uniform developmental depths.
(B) Costs summed across rounds at each weight form the aggregate layerwise
front, compared with first-, second- and third-cousin shuffles. These permute
measured leaf identities within groups sharing an ancestor two, three or four
canonical generations back, respectively, including cutoff/coverage leaves;
all internal states, parent-slot multiplicities and four root identities stay
fixed. Position and protein state move together. All EDGES edges are scored,
including unchanged internal edges. Random rebuild appears in an inset, using
the same endpoint transform but separately scaled axes, as in the terminal plot.
Dots show NSHOWN of NDRAWS draws; plus signs show full-draw means. Black crosses
mark natural lineage. Coordinates use endpoint normalization: the travel
optimum is \((0,1)\) and the state optimum \((1,0)\), separately per round in A
and using aggregate layerwise endpoints in B. Observations are not clipped.
These are descriptive randomization references, not calibrated significance
tests; separation does not establish biological efficiency.}
\label{fig:ce_full_tree_layerwise}
\end{figure}
""", (("INTERNAL", meta["internal"]), ("LEAVES", meta["leaves"]), ("EDGES", meta["edges"]),
                    ("NROUNDS", c["NROUNDS"]), ("ROUNDS", ", ".join(map(str, meta["round_edges"]))),
                    ("NWEIGHTS", c["NWEIGHTS"]), ("NSHOWN", c["NSHOWN"]), ("NDRAWS", c["NDRAWS"])))


def supplement(meta):
    c = counts(meta)
    return _fill(r"""\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figs_ce_full_tree_heuristics_panel.pdf}
\caption{\textbf{Full-tree reconstruction strategies and reference models.}
Five strategies are regenerated on the same pooled-travel cohort and protein
representation as Figure 7, using NWEIGHTS endpoint-inclusive weights. Each curve
shows its own sampled non-dominated cost pairs, not a cross-method global
optimum. All curves, natural lineage and reference draws use the \emph{same}
aggregate layerwise endpoints from Figure 7B. These common display anchors do
not make the methods' feasible sets equivalent. All EDGES edges and all
measured nodes are scored once in each reconstruction; stored parent arrays
are replay-validated for binary degree, four components and absence of cycles.
The three cousin shuffles permute measured leaves within canonical ancestor
groups two, three or four generations back, with internal states and roots fixed;
internal-layer shuffle permutes internal identities within canonical depths
and holds leaves fixed; full assignment shuffle permutes internal and terminal
identities separately. The latter two permit root identities to move.
Random rebuild fixes the four roots and samples a binary topology over measured
states. The Gaussian reference, retained only here, fits separate-clock
covariance to all observed edges, freely simulates internal states from four
observed roots, then clamps every measured leaf to its observed position and
protein state. Independent spatial processes are fitted per embryo and their
costs pooled; protein variance is per canonical transition. Terminal edges
reconnect free internal states to observed leaves. This is not joint endpoint
conditioning or a calibrated significance test; Gaussian internal states are
not restricted to the measured-state inventory used by reconstruction methods.
Dots show deterministic NSHOWN-draw subsets; plus signs use all NDRAWS draws.
Axes are continuous and observations are not clipped.}
\label{fig:ce_full_tree_heuristics_supplement}
\vspace{4mm}
\footnotesize
\begin{tabular}{@{}p{0.30\textwidth}p{0.66\textwidth}@{}}
\toprule
Strategy & Constraints and algorithm \\
\midrule
Layerwise assignment & Exact assignment within fixed contraction rounds;
natural parent-slot multiplicities and four root identities preserved. \\
Degree-constrained forest & Greedy Kruskal internal forest with one fixed root
per component and binary capacity; terminal Hungarian assignment. \\
Top-down rebuild & Greedy attachment of internal states from four fixed roots;
terminal Hungarian assignment to remaining binary openings. \\
Bottom-up by layer & Fixed topology; internal identities reassigned by
depth from the bottom up, including movable root identities. Leaves fixed. \\
Paired bottom-up & Minimum-weight daughter pairing followed by parent
assignment; four root identities emerge, rather than being fixed. \\
Terminal-only (excluded) & Historical constructed-state algorithm has unresolved
scope/scoring problems; no pooled curve generated or displayed. \\
\bottomrule
\end{tabular}
\end{figure}
""", (("EDGES", meta["edges"]), ("NWEIGHTS", c["NWEIGHTS"]),
                    ("NSHOWN", c["NSHOWN"]), ("NDRAWS", c["NDRAWS"])))
