"""Partial-forest cohort, assignment, exact-null and presentation safeguards."""
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.spatial.distance import cdist

from full_tree_pareto import cross_species_analysis as ca
from full_tree_pareto import cousin_references as cr
from full_tree_pareto import fig_cross_species as figmod
from full_tree_pareto import cross_species_publication as pub


def fixture():
    xyz = np.array([[0., 0, 0], [1., 2, 0], [-2., 1, 1], [1., 4, 1],
                    [3., 2, 2], [-3., 3, 1], [-1., 5, 0]])
    state = np.array([[40., 20], [100., 0], [20., 40], [0., 0], [100., 0], [0., 0], [100., 0]])
    return ca.forest(list("abcdefg"), np.array([-1, 0, 0, 1, 1, 2, 5]), cdist(xyz, xyz), cdist(state, state))


class PartialForestTests(unittest.TestCase):
    def test_ancestor_closure_stops_at_first_gap(self):
        canonical = dict(a=None, b="a", c="b", d="c", e="d")
        kept, heights = ca.ancestor_closure(["e"], {"a", "c", "d", "e"}, canonical)
        self.assertEqual(kept, {"c", "d", "e"})
        self.assertEqual(heights, dict(e=2))
        self.assertNotIn("a", kept)  # Do not resume above missing b.
        with self.assertRaises(ValueError):
            ca.ancestor_closure(["z"], {"a"}, canonical)

    def test_contraction_rounds_partition_mixed_capacity_forest(self):
        ctx = fixture()
        np.testing.assert_array_equal(np.concatenate(ctx.layers), [3, 4, 6, 1, 5, 2])
        self.assertEqual(sum(map(len, ctx.layers)), len(ctx.edges))
        self.assertEqual(np.sum(ctx.counts == 1), 2)

    def test_full_forest_exact_moments_include_fixed_internal_cost(self):
        ctx = fixture()
        groups = [list(range(len(ctx.leaves)))]
        stats = ca.exact_null(ctx, groups)
        perms = np.array(list(itertools.permutations(range(len(ctx.leaves)))))
        scores = cr.score_permutations(ctx, perms)
        expected = [*scores.mean(axis=0), *scores.var(axis=0), np.cov(scores.T, ddof=0)[0, 1]]
        np.testing.assert_allclose(stats, expected, rtol=1e-11, atol=1e-11)
        self.assertGreater(stats[0], ctx.travel[ctx.leaves, ctx.parents[ctx.leaves]].sum())

    def test_layerwise_sweep_exact_product_optimum_and_cost_sum(self):
        ctx = fixture()
        stats = ca.exact_null(ctx, [list(range(len(ctx.leaves)))])
        parents, costs, ties = ca.solve_sweep(ctx, stats, 4)
        candidates = []
        for permutations in itertools.product(*(list(itertools.permutations(ctx.parents[ch])) for ch in ctx.layers)):
            parent = ctx.parents.copy()
            for children, slots in zip(ctx.layers, permutations):
                parent[children] = slots
            candidates.append(parent)
        all_costs = ctx.score(np.array(candidates))
        for i, alpha in enumerate(np.linspace(0, 1, 5)):
            weights = np.array([alpha, 1-alpha])/np.sqrt(stats[2:4])
            self.assertAlmostEqual(float(costs[i] @ weights), float((all_costs @ weights).min()), places=11)
            subtotal = np.sum([[ctx.travel[ch, parents[i, ch]].sum(), ctx.state[ch, parents[i, ch]].sum()]
                               for ch in ctx.layers], axis=0)
            np.testing.assert_allclose(subtotal, costs[i])
        self.assertEqual(len(ties), len(ctx.layers))
        frame, metric, transform = ca.summarize(ctx, parents, costs, stats, "toy")
        np.testing.assert_allclose(frame[["D1", "D2"]].iloc[[0, -1]], [[1, 0], [0, 1]])
        np.testing.assert_allclose(np.column_stack(transform.inverse(frame.D1, frame.D2)), costs)
        self.assertTrue(0 <= metric["u_L"] <= 1)

    def test_root_capacity_cycle_and_round_guards(self):
        ctx = fixture()
        for change in ({0: 1}, {6: 0}, {2: 5, 6: 0}):
            parent = ctx.parents.copy()
            for child, p in change.items():
                parent[child] = p
            with self.assertRaises(ValueError):
                ca.validate_forests(ctx, parent)
        other_round = ctx.parents.copy()
        other_round[[2, 3]] = other_round[[3, 2]]
        ca.validate_forests(ctx, other_round)
        with self.assertRaises(ValueError):
            ca.validate_forests(ctx, other_round, roundwise=True)

    def test_random_rebuild_preserves_unary_and_binary_capacity(self):
        ctx = fixture()
        parents = ca.random_rebuilds(ctx, 200, 45)
        ca.validate_forests(ctx, parents)
        np.testing.assert_array_equal(parents, ca.random_rebuilds(ctx, 200, 45))
        self.assertGreater(len(np.unique(parents, axis=0)), 1)

    def test_cache_identity_and_hash_are_checked(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "case.npz"
            ca.save_cache(path, dict(scope="toy"), parents=fixture().parents)
            np.testing.assert_array_equal(ca.load_cache(path, dict(scope="toy"))["parents"], fixture().parents)
            with self.assertRaises(ValueError):
                ca.load_cache(path, dict(scope="other"))
            with path.open("ab") as handle:
                handle.write(b"changed")
            with self.assertRaises(ValueError):
                ca.load_cache(path, dict(scope="toy"))

    def test_live_terminal_anchored_coverage_and_geometry(self):
        info, nodes, coverage, matrices, components, positions, expression = ca.load_inputs()
        ctx = ca.forest(info["names"], info["parents"], *matrices["raw3d", "ce_protein"])
        self.assertEqual((len(ctx.names), len(ctx.edges), len(ctx.leaves), len(ctx.roots)), (485, 454, 187, 31))
        self.assertEqual([len(ch) for ch in ctx.layers], [187, 119, 80, 48, 18, 2])
        self.assertEqual(set(np.array(ctx.names)[ctx.leaves]), set(info["terminal_seeds"]))
        self.assertTrue(nodes.loc[nodes.boundary_root, "unavailable_parent"].all())
        self.assertEqual(int(coverage.retained.sum()), 485)
        self.assertEqual((min(info["ancestor_steps"].values()), max(info["ancestor_steps"].values())), (2, 6))
        for geometry in ca.GEOMETRIES:
            np.testing.assert_array_equal(matrices[geometry, "ce_protein"][0], matrices[geometry, "ce_rna"][0])
            for config in ca.PRIMARY:
                forest = ca.forest(info["names"], info["parents"], *matrices[geometry, config])
                self.assertAlmostEqual(forest.natural[0], 1., places=12)
                species = "cb" if config == "cb_rna" else "ce"
                pooled = np.mean([components[f"{geometry}__{species}__{r}"] for r in positions[species]], axis=0)
                np.testing.assert_allclose(pooled, forest.travel)
                np.testing.assert_allclose(cdist(expression[config], expression[config]), forest.state)
        for species in positions:
            for replicate, xyz in positions[species].items():
                raw = cdist(xyz[:, :2], xyz[:, :2])
                total = raw[ctx.edges, ctx.parents[ctx.edges]].sum()
                np.testing.assert_allclose(raw/total, components[f"xy__{species}__{replicate}"])

    def test_live_solver_matches_existing_layerwise_heuristic_on_partial_forest(self):
        info, _, _, matrices, *_ = ca.load_inputs()
        from terminal_pareto.lineage_metrics import build_lineage_tree_index
        from terminal_pareto.pareto_engine import build_ancestor_groups
        index = build_lineage_tree_index(ca.dl.load_json(ca.ROOT / "data/cell_lineage.json"))
        ctx = ca.forest(info["names"], info["parents"], *matrices["raw3d", "ce_protein"])
        groups = build_ancestor_groups([ctx.names[i] for i in ctx.leaves], index, 2)
        stats = ca.exact_null(ctx, groups)
        parents, costs, _ = ca.solve_sweep(ctx, stats, 2)
        for row, alpha in enumerate((0., .5, 1.)):
            combined = alpha*ctx.travel/np.sqrt(stats[2]) + (1-alpha)*ctx.state/np.sqrt(stats[3])
            old_parent = ca.pa.layerwise(ctx, combined)
            ca.validate_forests(ctx, old_parent, roundwise=True)
            weights = np.array([alpha, 1-alpha])/np.sqrt(stats[2:4])
            np.testing.assert_allclose(costs[row] @ weights, ctx.score(old_parent) @ weights)
            if alpha == .5:
                np.testing.assert_array_equal(old_parent, parents[row])


class PresentationTests(unittest.TestCase):
    def test_numbered_figures_are_explicit_about_partial_scope(self):
        record = dict(edges=454, cohort_size=485, terminal_count=187, roots=31,
                      round_edges=[187, 119, 80, 48, 18, 2], settings=dict(intervals=300, dense_intervals=1200, draws=10000))
        with tempfile.TemporaryDirectory() as folder:
            wrappers = pub.write_wrappers(Path(folder), record)
            for wrapper in wrappers:
                content = wrapper.read_text()
                number = pub.FIGURES[wrapper.stem][0]
                self.assertIn(rf"\setcounter{{figure}}{{{number-1}}}", content)
                self.assertNotIn("pilot", content)
                self.assertIn("454", content)
                self.assertIn("partial", content)
                self.assertIn("traditional embryo-tracking", content)
                self.assertNotIn("axial scale unverified", content)
            self.assertIn("No missing state is imputed", wrappers[0].read_text())
            self.assertIn(r"d_{CP}", wrappers[1].read_text())
            self.assertIn("Figure~10", wrappers[1].read_text())

    def test_layout_archive_is_hash_checked_and_recoverable(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            (run / "figures").mkdir()
            old = run / "figures/old.pdf"
            old.write_bytes(b"previous layout")
            archive = figmod.archive_previous(run, "figures")
            self.assertEqual(ca.digest(old), ca.digest(archive / "figures/old.pdf"))
            self.assertEqual(json.loads((archive / "manifest.json").read_text())["files"]["old.pdf"], ca.digest(old))


class ProductionTests(unittest.TestCase):
    def prepare(self, folder):
        root = Path(folder)
        run, target = root / "run", root / "output/publication"
        (run / "publication").mkdir(parents=True)
        target.mkdir(parents=True)
        for name in pub.publication_files():
            (run / "publication" / name).write_bytes(f"numbered {name}".encode())
        (run / "publication/publication_manifest.json").write_text('{"analysis_id": "toy"}\n')
        for name in pub.LEGACY_FIG8:
            (target / name).write_bytes(f"legacy {name}".encode())
        (target / "fig7_ce_full_tree_layerwise.pdf").write_bytes(b"preserve Figure 7 exactly")
        (target / "table.pdf").write_bytes(b"preserve unrelated table")
        return run, target

    @patch.object(pub, "verify", return_value={"analysis_id": "toy"})
    def test_additive_release_preserves_others_and_archives_retired_figure8(self, _verify):
        with tempfile.TemporaryDirectory() as folder:
            run, target = self.prepare(folder)
            before = pub.inventory(target)
            pub.promote(run, production=target)
            release = pub.verify_release(target)
            self.assertEqual(release["figure_numbers"], {stem: number for stem, (number, _) in pub.FIGURES.items()})
            self.assertEqual(set(release["retired_figure8_files"]), set(pub.LEGACY_FIG8))
            self.assertEqual(set(pub.inventory(target)), pub.publication_files() | {
                pub.ASSEMBLY_COPY, pub.RELEASE_MANIFEST, "fig7_ce_full_tree_layerwise.pdf", "table.pdf"})
            self.assertEqual(pub.inventory(Path(release["previous_publication"]) / "publication"), before)
            for name, value in release["preserved_files"].items():
                self.assertEqual(ca.digest(target / name), before[name])
            self.assertFalse(any((target / name).exists() for name in pub.LEGACY_FIG8))

    @patch.object(pub, "verify", return_value={"analysis_id": "toy"})
    def test_failed_final_verification_restores_original_production(self, _verify):
        with tempfile.TemporaryDirectory() as folder:
            run, target = self.prepare(folder)
            before = pub.inventory(target)
            original_verify = pub.verify_release

            def fail_final(publication, **kwargs):
                if Path(publication) == target:
                    raise ValueError("injected final verification failure")
                return original_verify(publication, **kwargs)

            with patch.object(pub, "verify_release", side_effect=fail_final), self.assertRaisesRegex(ValueError, "injected"):
                pub.promote(run, production=target)
            self.assertEqual(pub.inventory(target), before)
            self.assertFalse(list(target.parent.glob(".cross-species-stage-*")))
            self.assertFalse(list((target.parent / "legacy/cross_species_releases").iterdir()))

    @patch.object(pub, "verify", return_value={"analysis_id": "toy"})
    def test_ambiguous_retirement_and_symlinks_are_refused(self, _verify):
        for kind in ("unknown_figure8", "pinned_release", "symlink", "directory_symlink"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                run, target = self.prepare(folder)
                if kind == "unknown_figure8":
                    (target / "fig8_unknown.pdf").write_bytes(b"not authorized")
                elif kind == "pinned_release":
                    (target / "release_manifest.json").write_text(json.dumps(dict(files={pub.LEGACY_FIG8[0]: "pinned"})))
                elif kind == "symlink":
                    (target / "alias.pdf").symlink_to(target / "table.pdf")
                else:
                    alias = target.parent / "alias"
                    alias.symlink_to(target, target_is_directory=True)
                before = pub.inventory(target)
                with self.assertRaises(ValueError):
                    pub.promote(run, production=alias if kind == "directory_symlink" else target)
                self.assertEqual(pub.inventory(target), before)


if __name__ == "__main__":
    unittest.main()
