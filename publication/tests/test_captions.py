"""Caption content and wrapper numbering (ported from the retired back-end presentation tests)."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from publication import assembly
from publication.captions import cross_species, full_tree, terminal

TERMINAL_META = dict(edges=187, cb_reference_size=196)
PARTIAL_META = dict(edges=454, nodes=485, terminals=187, roots=31, round_edges=[187, 119, 80, 48, 18, 2],
                    weights=301, dense_weights=1201, draws=10000)
POOLED_META = dict(edges=974, leaves=491, internal=487, round_edges=[491, 223, 123, 66, 38, 18, 10, 5],
                   weights=301, draws=10000, display_draws=1000)


class CrossSpeciesCaptions(unittest.TestCase):
    def test_terminal_figures_8_9(self):
        comparison, overlay = cross_species.terminal_comparison(TERMINAL_META), cross_species.terminal_overlay(TERMINAL_META)
        for text in (comparison, overlay):
            self.assertIn("187", text)
            self.assertIn("traditional embryo-tracking techniques", text)
            self.assertIn("provisional", text)
            self.assertNotIn("axial scale unverified", text)
        self.assertIn("196-edge reference", comparison)
        self.assertIn(r"\textbf{(C)}", overlay)
        self.assertIn(r"d_{LP}", overlay)
        self.assertIn(r"d_{CP}", overlay)

    def test_partial_forest_figures_10_11_state_their_scope(self):
        comparison, overlay = cross_species.full_tree_comparison(PARTIAL_META), cross_species.full_tree_overlay(PARTIAL_META)
        for text in (comparison, overlay):
            self.assertIn("454", text)
            self.assertIn("partial", text)
            self.assertNotIn("pilot", text)
            self.assertIn("traditional embryo-tracking", text)
        self.assertIn("No missing state is imputed", comparison)
        self.assertIn(r"d_{CP}", overlay)
        self.assertIn("Figure~10", overlay)


class FullTreeCaptions(unittest.TestCase):
    def test_figure7_and_s4_fill_every_count(self):
        figure7, supplement = full_tree.figure7(POOLED_META), full_tree.supplement(POOLED_META)
        for text in (figure7, supplement):
            for token in ("INTERNAL", "LEAVES", "EDGES", "ROUNDS", "NWEIGHTS", "NSHOWN", "NDRAWS"):
                self.assertNotIn(token, text)
            self.assertIn("301", text)
            self.assertIn("10,000", text)
        self.assertNotIn("Gaussian", figure7)
        self.assertIn("third-cousin", figure7)
        self.assertIn("Eight asynchronous", figure7)
        self.assertIn("retained only here", supplement)

    def test_reference_sets_keep_gaussian_in_the_supplement_only(self):
        from full_tree_pareto import pooled_analysis as pa
        from publication.adapters.full_tree import reference_sets
        main, supplement = reference_sets()
        self.assertNotIn(pa.REFERENCE, main)
        self.assertIn(pa.REFERENCE, supplement)


class TerminalWrappers(unittest.TestCase):
    def test_numbering_and_organization(self):
        with tempfile.TemporaryDirectory() as tmp:
            texts = {}
            for stem, (number, body) in terminal.CAPTIONED_FIGURES.items():
                texts[stem] = assembly.write_wrapper(Path(tmp), stem, number=number, body=body()).read_text()
                self.assertIn(r"\documentclass[10pt]{article}", texts[stem])
        self.assertIn(r"\renewcommand{\thefigure}{1 amendment}", texts["fig1_ce_endpoint_normalization_amendment"])
        self.assertIn(r"\renewcommand{\thefigure}{S3}", texts["figS3_ce_tracking_robustness"])
        self.assertIn(r"\setcounter{figure}{4}", texts["fig5_ce_canonical_summary"])
        self.assertIn(r"\input{fig3A_ce_null_models.tex}", texts["fig3_ce_terminal_pareto_supporting"])
        self.assertIn("Figure~1 amendment", texts["fig5_ce_canonical_summary"])
        self.assertIn("leave-one-geometry-out", texts["figS3_ce_tracking_robustness"])
        self.assertIn("fig5B_ce_canonical_all_subtrees.pdf", texts["fig5_ce_canonical_summary"])


if __name__ == "__main__":
    unittest.main()
