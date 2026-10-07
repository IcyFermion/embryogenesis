# Track B: figure redesign

Branch `figure-ideas`; sketches in `publication/sketches/`, shown on the
review page (see `PUBLICATION_FIGURES_HANDOFF.md` section 0). Production
figures change only through a normal family build and release.

## B1. Lineage rewiring view (feedback item 1)

Show, for the primary 275-edge terminal cohort, how the assignment differs from
the natural lineage at three points of the global front: the travel optimum T,
the cell-state optimum S, and the maximum-edge-retention assignment (also mark
P*). No new optimization: per-weight assignments are cached by the primary
global analysis (load through the publication terminal adapter) and, for the
cross-species cohort, in `assignments_and_costs.npz` ("parent slot i receives
child permutation[i]").

- Draw: terminals in lineage order; for each reassigned child, an arc from its
  natural parent to its assigned parent, colored by tree distance (path length
  between the two parents). This also shows what tree distance is, which
  covers Figure 3B's concept.
- Rank edge changes by tree distance and by cost gain; list the standouts.
- Hypothesis to test: swaps between left-right symmetric lineages producing
  equivalent fates (e.g. ABpl/ABpr descendants), and same-fate swaps. A
  homolog map needs a literature source, not name heuristics alone.
- Output an edge-change table per point (child, natural parent, assigned
  parent, tree distance, travel and state change, fates); Track C reuses it.

## B2. Figure 3 (item 3)

Proposal agreed in principle: B1 becomes the new Figure 3. 3A moves to the
schematic Figure 1 (co-author owned) or the supplement; 3B-C move to the
supplement. Decide after the B1 sketch.

## B3. Rename "maximum edge retention" (item 2)

On hold for the author's choice. When decided: change the `retention` concept
label in `publication/notation.py`, then the literal strings that bypass it
(grep for "edge retention" and "Maximum edge retention" in `publication/`),
rebuild all four families, rehearse and release.

## B4. Figure 4 (item 4)

Two or three alternatives to compare:
- icicle layout (lineage blocks, width proportional to terminal count, color =
  d_LP, on-front marked by a simple symbol);
- fate composition moved out of the tree into a thin strip;
- Figure 4 merged with Figure 5B into one annotated tree.
