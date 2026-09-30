import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from full_tree_pareto import cousin_references as cr
from full_tree_pareto import pooled_analysis as pa
from terminal_pareto.pareto_engine import compute_group_shuffle_costs


class CousinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ctx = pa.build_context()

    def test_first_cousin_group_membership_matches_existing_null(self):
        actual = {tuple(sorted(g)) for g in cr.groups_for_context(self.ctx, 2)}
        expected = {tuple(sorted(g)) for g in self.ctx.cousin_groups if len(g) > 1}
        self.assertEqual(actual, expected)

    def test_groups_and_permutations_preserve_forests_and_scores(self):
        ctx = self.ctx
        for steps in cr.ANCESTOR_STEPS:
            groups = cr.groups_for_context(ctx, steps)
            joined = np.concatenate(groups)
            self.assertEqual(len(joined), len(np.unique(joined)))
            perms = cr.sample_permutations(len(ctx.leaves), groups, 12, 42)
            cr.validate_permutations(perms, len(ctx.leaves), groups)
            scores = cr.score_permutations(ctx, perms)
            for perm, cost in zip(perms, scores):
                mapping = np.arange(len(ctx.names))
                mapping[ctx.leaves] = ctx.leaves[perm]
                parent = ctx.parents.copy()
                parent[mapping[ctx.edges]] = mapping[ctx.parents[ctx.edges]]
                np.testing.assert_array_equal(parent[ctx.internal], ctx.parents[ctx.internal])
                np.testing.assert_array_equal(np.bincount(parent[ctx.edges], minlength=len(ctx.names)),
                                              np.bincount(ctx.parents[ctx.edges], minlength=len(ctx.names)))
                np.testing.assert_allclose(pa.validate_assignment(ctx, parent), cost, atol=1e-12)
            # Exact shared RNG/group orientation agrees with terminal sampler.
            costs = compute_group_shuffle_costs(
                ctx.travel[np.ix_(ctx.parents[ctx.leaves], ctx.leaves)],
                ctx.state[np.ix_(ctx.parents[ctx.leaves], ctx.leaves)],
                groups, n_random=12, seed=42)
            constant = cr.score_permutations(ctx, np.arange(len(ctx.leaves))[None, :])[0]
            natural_leaves = np.array([ctx.travel[ctx.leaves, ctx.parents[ctx.leaves]].sum(),
                                       ctx.state[ctx.leaves, ctx.parents[ctx.leaves]].sum()])
            np.testing.assert_allclose(scores, np.column_stack(costs) + constant-natural_leaves, atol=1e-12)

    def test_invalid_permutations_rejected(self):
        groups = [np.array([0, 1]), np.array([2, 3])]
        for invalid in (np.array([[0, 0, 2, 3]]), np.array([[2, 1, 0, 3]]),
                        np.array([[0., 1., 2., 3.]])):
            with self.assertRaises(ValueError):
                cr.validate_permutations(invalid, 4, groups)
        valid = cr.sample_permutations(5, groups, 10, 4)
        cr.validate_permutations(valid, 5, groups)
        np.testing.assert_array_equal(valid[:, 4], 4)

    def test_cache_reuse_hash_identity_and_replay_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            with self.assertRaises(FileNotFoundError):
                cr.build(self.ctx, run, draws=12, seed=42, layout_only=True)
            first = cr.build(self.ctx, run, draws=12, seed=42)
            second = cr.build(self.ctx, run, draws=12, seed=42, layout_only=True)
            for name in cr.REFERENCES:
                np.testing.assert_array_equal(first[name], second[name])
            with self.assertRaises(ValueError):
                cr.build(self.ctx, run, draws=13, seed=42, layout_only=True)
            path = run / "analysis/cousin_shuffles/draws.npz"
            # Even a rehashed incorrect cost must fail independent replay.
            with np.load(path) as data:
                arrays = {key: data[key] for key in data.files}
            arrays["costs_0"][0, 0] += .1
            identity = cr.cache_identity(self.ctx, 12, 42)
            pa.save_cache(path, identity, **arrays)
            validation = path.parent / "validation.json"
            record = json.loads(validation.read_text())
            record["hashes"] = {p.name: pa.digest(p) for p in (path, path.with_suffix(".json"))}
            pa.write_json(validation, record)
            with self.assertRaises(AssertionError):
                cr.build(self.ctx, run, draws=12, seed=42, layout_only=True)


if __name__ == "__main__":
    unittest.main()
