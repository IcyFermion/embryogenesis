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

## B5. Third round of feedback: manuscript set of 2026-11-08 (draft 5)

Author and supervisor, on draft 4 (review page version 11). Branch
`figures-v2`; build `publication/output/candidates/v2_draft5`, manuscript set
`publication/output/manuscript/v2_20261108`.

**Numbering.** Figure 2 is new (cell-state features, Track A3); the former
Figures 2-7 are now 3-8: 3 terminal front, 4 lineage similarity, 5 subtrees,
6 cell types (released panels, placeholder), 7 full tree by layer, 8
*C. briggsae*. Captions in `publication/v2/captions.py` follow the new numbers.

Items and what was done:
1. **C. briggsae counterparts of 3A and 6B (now 4A and 7B).** New Figure 8
   row of two half-width panels. 8B: lineage similarity along the 187-edge
   *C. briggsae* front (M 0.74; T 0.50; S 0.21). 8C: the partial forest by
   layer (layers 4-8; layer 3 has 2 cells and no resolved front), *C. briggsae*
   RNA and AF16 3D travel, at the birth sizes of the matched *C. elegans*
   cells; the d_LP axis is shared with 7B. Loader `data.briggsae_layers`
   (validated round references and the hash-checked Track E stage table).
   d_LP 0.07, 0.04, 0.07, 0.09, 0.03 (L4-L8); lowest % differing 22, 17, 52,
   43, 26: the % peak sits at layer 6 as in *C. elegans*, the d_LP peak does not.
2. **Figure 5 (subtrees) order.** Left-right counterparts adjacent and
   bracketed: ABpl/ABpr, ABplp/ABprp, MSa/MSp, ABalp/ABara (by fate, Track D0).
   ABal/ABar are sisters, kept together above their pair.
   `panels.SUBTREE_ORDER`, `SUBTREE_PAIRS`.
3. **S1 metrics into Figure 5.** Third column, cousin-relative distance r
   (dashed reference at r = 1), from `canonical_display.relative_distance`.
   Two group columns now (Root & AB + ABa | ABp + P1); the fate legend fills
   the space under the left column. Released S1 is retired from the set.
4. **S2, S4, S5 to the backlog.** Copied to `backlog/` of the manuscript set
   with their captions at the end of `captions.tex`; not in `review.pdf`.
5. **S3 (now S1) note for the main text.** In the set's README: native-front
   d_LP 0.020, 0.029, 0.013 in embryos 1-3 (mean 0.021), pooled 0.015 (the
   global d_LP). Pooled is below two of three embryos and the mean, not the
   smallest. Values recorded by the build (`build.json`,
   `main_text_notes`).
6. **New Figure 2.** See Track A, A3.
7. **Common palette.** `publication/palette.py`: one colour per concept
   (MEASURE, NULL, FATE, SISTERS, FEATURE_SET, NEUTRAL) and `sheet()` (built
   as `palette.pdf`). All v2 panels draw from it; Figure 7B-C are now drawn in
   `publication.v2.panels` (`layer_position`, `layer_cost`) so they follow the
   palette. `style.py` is unchanged, because it is hashed into every family's
   build identity; the values are equal today. Five fate classes share a hue
   with a measure (`palette.SHARED_HUES`), shown with an asterisk on the sheet.
8. **E lineage in Figure 4's tree** (side question). Drawn correctly; the
   cohort rule is the cause. Track D, D7.

**Decisions (author, 2026-11-09).**
- Terminal cohort: keep the 275 edges; the 321-edge check is in the backlog
  (Track D, D7).
- Fate/measure shared hues: kept for now; the theme is decided when the
  whole manuscript wraps up.
- Figure 2E dropped: the memory measure is good (Methods) but the picture is
  not convincing. Figure 2 is A-D.
- Draft captions: good for now, revisited later.
- *C. elegans* reference in Figure 8: not added (it would have to go into 8A,
  8B and 8C alike). On the same 187 matched cells *C. elegans* RNA differs
  most at layer 6 (d_LP 0.23 against 0.075), and the weaker *C. briggsae*
  tracking, especially in z, makes that difference hard to claim. Figure 8
  stays *C. briggsae* alone.

**Release set (author request, 2026-11-09).** The manuscript set is now laid
out like `publication/output/production`: every figure and table has a
standalone captioned LaTeX wrapper (`figN_<stem>.tex`, with `\label`) next to
its compiled page and the graphics it includes. Main figures include the v2
manuscript-mode figure as `<stem>_panel.pdf`, at most 70% of the text height
so the caption fits. Figure 6 and S1 are the released wrappers, renumbered
(S1's files renamed `figS1_*`); the backlog keeps released wrappers unchanged.
`captions.tex` and `graphics/` are gone (each caption lives in its wrapper);
`review.tex` only joins the compiled pages. Set:
`publication/output/manuscript/v2_release_20261109`, from build
`candidates/v2_release_20261109`.
