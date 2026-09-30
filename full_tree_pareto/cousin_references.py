"""Auditable higher-cousin nulls without invalidating frozen solver caches."""
from pathlib import Path
import json

import numpy as np

from full_tree_pareto import pooled_analysis as pa
from terminal_pareto.lineage_metrics import build_lineage_tree_index
from terminal_pareto.pareto_engine import build_ancestor_groups
from utils import load_json

REFERENCES = ("Second-cousin shuffle", "Third-cousin shuffle")
ANCESTOR_STEPS = (3, 4)
VERSION = "measured-leaf-cousin-shuffles-1"


def groups_for_context(ctx, steps):
    """Use the terminal pipeline's canonical ancestry, including unscored anchors."""
    index = build_lineage_tree_index(load_json(pa.ROOT / "data/cell_lineage.json"))
    names = [ctx.names[i] for i in ctx.leaves]
    return [np.asarray(g, dtype=int) for g in build_ancestor_groups(names, index, steps)]


def sample_permutations(n, groups, draws, seed):
    rng = np.random.default_rng(seed)
    result = np.tile(np.arange(n, dtype=np.int32), (draws, 1))
    for perm in result:
        for group in groups:
            perm[group] = rng.permutation(group)
    return result


def validate_permutations(permutations, n, groups):
    if permutations.ndim != 2 or permutations.shape[1] != n:
        raise ValueError("Incomplete leaf permutation cache")
    if not np.issubdtype(permutations.dtype, np.integer):
        raise ValueError("Leaf permutations must be integers")
    if not np.array_equal(np.sort(permutations, axis=1),
                          np.broadcast_to(np.arange(n), permutations.shape)):
        raise ValueError("Leaf identities must occur exactly once per draw")
    # Singleton / ancestry-unavailable leaves each have a private fixed group.
    membership = np.arange(n)
    for group in groups:
        membership[group] = group[0]
    if not np.array_equal(membership[permutations],
                          np.broadcast_to(membership, permutations.shape)):
        raise ValueError("Permutation crosses canonical ancestor groups")


def score_permutations(ctx, permutations):
    """Reassign only measured leaves; score all edges, with internal costs fixed."""
    fixed = np.setdiff1d(ctx.edges, ctx.leaves)
    constant = np.array([ctx.travel[fixed, ctx.parents[fixed]].sum(),
                         ctx.state[fixed, ctx.parents[fixed]].sum()])
    parents = ctx.parents[ctx.leaves]
    costs = np.empty((len(permutations), 2))
    for start in range(0, len(permutations), 250):
        chunk = permutations[start:start+250]
        children = ctx.leaves[chunk]
        costs[start:start+len(chunk), 0] = constant[0] + ctx.travel[children, parents].sum(axis=1)
        costs[start:start+len(chunk), 1] = constant[1] + ctx.state[children, parents].sum(axis=1)
    return costs


def cache_identity(ctx, draws, seed):
    sources = (Path(__file__), pa.ROOT / "full_tree_pareto/pooled_analysis.py",
               pa.ROOT / "terminal_pareto/pareto_engine.py",
               pa.ROOT / "terminal_pareto/lineage_metrics.py",
               pa.ROOT / "terminal_pareto/data_loader.py", pa.ROOT / "utils.py")
    return dict(version=VERSION, context_key=ctx.identity["context_key"],
                draws=draws, seeds=[seed+1, seed+2], ancestor_steps=list(ANCESTOR_STEPS),
                numpy_version=np.__version__,
                source_hashes={str(p.relative_to(pa.ROOT)): pa.digest(p) for p in sources})


def build(ctx, run, *, draws, seed, layout_only=False):
    """Save draws and leaf permutations separately; replay every draw on reuse."""
    out = Path(run) / "analysis/cousin_shuffles"
    path = out / "draws.npz"
    identity = cache_identity(ctx, draws, seed)
    groups = [groups_for_context(ctx, steps) for steps in ANCESTOR_STEPS]
    if path.exists():
        pa.verify_cache(path, identity)
        record = json.loads((out / "validation.json").read_text())
        if record["identity"] != identity:
            raise ValueError("Cousin validation identity mismatch")
        for filename, sha in record["hashes"].items():
            if pa.digest(out / filename) != sha:
                raise ValueError(f"Changed cousin reference file: {filename}")
        with np.load(path) as data:
            arrays = {key: data[key] for key in data.files}
    else:
        if layout_only:
            raise FileNotFoundError("Build the higher-cousin caches before layout-only rendering")
        out.mkdir(parents=True, exist_ok=True)
        arrays = {}
        summaries = []
        for i, (name, steps, family_groups) in enumerate(zip(REFERENCES, ANCESTOR_STEPS, groups)):
            permutations = sample_permutations(len(ctx.leaves), family_groups, draws, seed+i+1)
            costs = score_permutations(ctx, permutations)
            arrays[f"permutations_{i}"] = permutations
            arrays[f"costs_{i}"] = costs
            summaries.append(dict(model=name, ancestor_steps=steps, seed=seed+i+1,
                                  groups=len(family_groups),
                                  group_sizes=[len(g) for g in family_groups],
                                  shuffled_leaves=sum(map(len, family_groups)),
                                  travel_mean=float(costs[:, 0].mean()), state_mean=float(costs[:, 1].mean()),
                                  travel_sd=float(costs[:, 0].std()), state_sd=float(costs[:, 1].std())))
        pa.save_cache(path, identity, **arrays)
        pa.write_json(out / "validation.json", dict(identity=identity, summaries=summaries,
                      nodes=len(ctx.names), edges=len(ctx.edges), leaves=len(ctx.leaves),
                      fixed_internal_states=True, fixed_root_identities=True,
                      hashes={p.name: pa.digest(p) for p in (path, path.with_suffix(".json"))}))
    result = {}
    for i, (name, family_groups) in enumerate(zip(REFERENCES, groups)):
        permutations, costs = arrays[f"permutations_{i}"], arrays[f"costs_{i}"]
        if permutations.shape != (draws, len(ctx.leaves)) or costs.shape != (draws, 2):
            raise ValueError("Incomplete higher-cousin draws")
        validate_permutations(permutations, len(ctx.leaves), family_groups)
        np.testing.assert_allclose(score_permutations(ctx, permutations), costs, rtol=1e-12, atol=1e-12)
        if not np.isfinite(costs).all():
            raise ValueError("Nonfinite higher-cousin scores")
        result[name] = costs
        print(f"Validated {name}: {draws} draws, mean {costs.mean(axis=0)}", flush=True)
    return result
