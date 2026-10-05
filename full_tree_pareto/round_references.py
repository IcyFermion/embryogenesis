"""Per-round cousin references and canonical metrics for layerwise fronts.

python -m full_tree_pareto.round_references [--target pooled|cross-species|all]

Each contraction round is its own assignment problem: the round's children
are matched to the round's parent slots. The per-round null shuffles only that
round's children over those slots, within groups sharing a canonical ancestor
two, three or four generations back (first-, second-, third-cousin), or over
all slots (full random). Internal states, roots and every other round stay
fixed, and only the round's edges are scored. In round 1 (the measured leaves)
this is the existing leaf-cousin shuffle restricted to leaf edges.

The null mean C is the exact first-cousin mean of the round. Canonical metrics
(u_L, d_LP, d_NP = displayed d_CP, retention at P*) follow
``full_tree_pareto.cross_species_analysis.summarize`` in each round's own
endpoint coordinates. Fronts are the cached layerwise sweeps (global leaf
first-cousin SD scaling, shared with the aggregate); nothing is re-solved.

An add-on: frozen caches and their hash-pinned modules are read, never
modified. Outputs live in ``<run>/analysis/round_references/``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from full_tree_pareto import cross_species_analysis as fcs
from full_tree_pareto import pooled_analysis as pa
from terminal_pareto.front_coordinates import DegenerateEndpointSpan, EndpointTransform
from terminal_pareto.lineage_metrics import build_lineage_tree_index
from terminal_pareto.pareto_engine import build_ancestor_groups
from terminal_pareto.subtree_analysis import exact_cousin_stats
from utils import load_json

VERSION = "per-round-cousin-references-1"
FAMILIES = ("first_cousin", "second_cousin", "third_cousin", "full_random")
ANCESTOR_STEPS = dict(first_cousin=2, second_cousin=3, third_cousin=4)
DRAWS = 2000
SEED = 4200
SUBDIR = "analysis/round_references"


def tree_index():
    return build_lineage_tree_index(load_json(pa.ROOT / "data/cell_lineage.json"))


def round_groups(names, children, index, family):
    """Groups as positions within ``children``; full random is one group of all."""
    if family == "full_random":
        return [np.arange(len(children))] if len(children) > 1 else []
    groups = build_ancestor_groups([names[i] for i in children], index, ANCESTOR_STEPS[family])
    return [np.asarray(g, dtype=int) for g in groups]


def sample_permutations(n, groups, draws, seed):
    """Within-group permutations of positions 0..n-1 (fixed outside groups)."""
    rng = np.random.default_rng(seed)
    result = np.tile(np.arange(n, dtype=np.int32), (draws, 1))
    for perm in result:
        for group in groups:
            perm[group] = rng.permutation(group)
    return result


def validate_permutations(permutations, n, groups):
    if permutations.ndim != 2 or permutations.shape[1] != n:
        raise ValueError("Incomplete round permutation cache")
    if not np.array_equal(np.sort(permutations, axis=1), np.broadcast_to(np.arange(n), permutations.shape)):
        raise ValueError("Round children must occur exactly once per draw")
    membership = np.arange(n)
    for group in groups:
        membership[group] = group[0]
    if not np.array_equal(membership[permutations], np.broadcast_to(membership, permutations.shape)):
        raise ValueError("Permutation crosses its round's ancestor groups")


def round_blocks(travel, state, parents, children):
    """Child-by-slot cost blocks; slot j is the natural parent of child j."""
    slots = parents[children]
    return travel[np.ix_(children, slots)], state[np.ix_(children, slots)]


def score_permutations(blocks, permutations):
    """Draw d sends child perm[d, j] to slot j; returns (draws, 2) round costs."""
    xm, em = blocks
    columns = np.arange(xm.shape[1])
    return np.column_stack((xm[permutations, columns].sum(axis=1), em[permutations, columns].sum(axis=1)))


def exact_moments(blocks, groups):
    """(mean travel, mean state, var travel, var state, cov) of the round's first-cousin null."""
    return np.array(exact_cousin_stats(*blocks, groups), dtype=float)


def landmarks(costs, natural, null_mean, retention, reference_id):
    """Canonical metrics of one sweep, as in ``fcs.summarize`` (travel optimum at the last index)."""
    try:
        transform = EndpointTransform.from_endpoints(
            reference_analysis_id=reference_id,
            travel_optimum_assignment_id=f"{reference_id}:{len(costs)-1}",
            state_optimum_assignment_id=f"{reference_id}:0",
            travel_optimum_costs=costs[-1], state_optimum_costs=costs[0])
    except DegenerateEndpointSpan:
        return None, dict(degenerate=True), None
    x, y = transform.transform(costs[:, 0], costs[:, 1])
    nx, ny = map(float, transform.transform(*natural))
    qx, qy = map(float, transform.transform(*null_mean))
    nearest = int(np.argmin(np.hypot(x-nx, y-ny)))
    arc = np.r_[0., np.cumsum(np.hypot(np.diff(x[::-1]), np.diff(y[::-1])))]
    u = (arc / arc[-1])[::-1]
    candidates = np.flatnonzero(retention == retention.max())
    maximum = int(candidates[np.argmin(np.abs(u[candidates]-u[nearest]))])
    frame = pd.DataFrame(dict(sweep_index=np.arange(len(costs)), alpha_travel=np.linspace(0, 1, len(costs)),
        travel=costs[:, 0], cell_state=costs[:, 1], D1=x, D2=y, u=u, edge_retention=retention,
        nearest_assignment=np.arange(len(costs)) == nearest, maximum_retention=np.arange(len(costs)) == maximum))
    metrics = dict(degenerate=False, natural_D1=nx, natural_D2=ny, null_D1=qx, null_D2=qy,
        nearest_index=nearest, maximum_index=maximum, u_L=float(u[nearest]),
        d_LP=float(np.hypot(x[nearest]-nx, y[nearest]-ny)), d_NP=float(np.hypot(x[nearest]-qx, y[nearest]-qy)),
        max_edge_retention=float(retention[maximum]), retention_at_nearest=float(retention[nearest]),
        retention_travel=float(retention[-1]), retention_state=float(retention[0]),
        distinct_front_points=int(len(np.unique(np.round(costs, 12), axis=0))))
    return frame, metrics, transform


def draw_references(names, parents, layers, index, draws, seed):
    """Permutations per round and family; shared by every configuration on the same forest."""
    arrays, groups = {}, {}
    for r, children in enumerate(layers, 1):
        for f, family in enumerate(FAMILIES):
            g = round_groups(names, children, index, family)
            groups[r, family] = g
            arrays[f"round_{r}__{family}"] = sample_permutations(len(children), g, draws, seed + 100*r + f)
    return arrays, groups


def analyze_forest(names, parents, layers, travel, state, sweep, permutations, groups, *, label, reference_id):
    """Per-round clouds, exact moments, fronts and metrics for one configuration."""
    costs_out, fronts, metrics, transforms = {}, [], [], {}
    for r, children in enumerate(layers, 1):
        blocks = round_blocks(travel, state, parents, children)
        natural = np.array([np.trace(blocks[0]), np.trace(blocks[1])])
        first = groups[r, "first_cousin"]
        stats = exact_moments(blocks, first) if first else np.r_[natural, 0., 0., 0.]
        row = dict(config=label, round=r, edges=len(children), natural_travel=natural[0], natural_state=natural[1],
                   null_mean_travel=stats[0], null_mean_state=stats[1], null_variance_travel=stats[2],
                   null_variance_state=stats[3], null_covariance=stats[4])
        for family in FAMILIES:
            perms = permutations[f"round_{r}__{family}"]
            validate_permutations(perms, len(children), groups[r, family])
            cloud = score_permutations(blocks, perms)
            costs_out[f"{label}__round_{r}__{family}"] = cloud
            row[f"{family}_groups"] = len(groups[r, family])
            row[f"{family}_shuffled"] = int(sum(map(len, groups[r, family])))
            row[f"{family}_travel_mean"], row[f"{family}_state_mean"] = cloud.mean(axis=0)
        round_costs = np.column_stack((travel[children[None, :], sweep[:, children]].sum(axis=1),
                                       state[children[None, :], sweep[:, children]].sum(axis=1)))
        retention = (sweep[:, children] == parents[children]).mean(axis=1)
        frame, metric, transform = landmarks(round_costs, natural, stats[:2], retention, f"{reference_id}:round_{r}")
        row.update(metric)
        metrics.append(row)
        if frame is not None:
            fronts.append(frame.assign(config=label, round=r))
            transforms[f"{label}__round_{r}"] = transform.metadata()
    return costs_out, fronts, metrics, transforms


def _identity(**fields):
    sources = [Path(__file__), pa.ROOT / "full_tree_pareto/pooled_analysis.py",
               pa.ROOT / "full_tree_pareto/cross_species_analysis.py",
               pa.ROOT / "terminal_pareto/front_coordinates.py", pa.ROOT / "terminal_pareto/pareto_engine.py",
               pa.ROOT / "terminal_pareto/subtree_analysis.py", pa.ROOT / "terminal_pareto/lineage_metrics.py",
               pa.ROOT / "utils.py"]
    return dict(version=VERSION, families=list(FAMILIES), ancestor_steps=ANCESTOR_STEPS,
                numpy_version=np.__version__,
                source_hashes={str(p.relative_to(pa.ROOT)): pa.digest(p) for p in sources}, **fields)


def _write_or_verify(out, identity, permutations, analyze):
    """Draw once; on reuse verify identity and file hashes, then replay every score."""
    path = out / "draws.npz"
    if path.exists():
        pa.verify_cache(path, identity)
        record = json.loads((out / "validation.json").read_text())
        if record["identity"] != identity:
            raise ValueError("Round-reference validation identity mismatch")
        for name, sha in record["hashes"].items():
            if pa.digest(out / name) != sha:
                raise ValueError(f"Changed round-reference file: {name}")
        with np.load(path) as data:
            saved = {k: data[k] for k in data.files}
        permutations = {k: v for k, v in saved.items() if "__costs" not in k}
        costs, fronts, metrics, transforms = analyze(permutations)
        for key, value in costs.items():
            np.testing.assert_allclose(value, saved[f"{key}__costs"], rtol=1e-12, atol=1e-12)
        for name, frame in (("fronts.csv", fronts), ("metrics.csv", metrics)):
            saved_frame = pd.read_csv(out / name)
            numeric = frame.select_dtypes(include="number").columns
            np.testing.assert_allclose(frame[numeric].to_numpy(float), saved_frame[numeric].to_numpy(float),
                                       rtol=1e-11, atol=1e-12, equal_nan=True)
        return costs, fronts, metrics, record["transforms"]
    out.mkdir(parents=True, exist_ok=True)
    costs, fronts, metrics, transforms = analyze(permutations)
    pa.save_cache(path, identity, **permutations, **{f"{k}__costs": v for k, v in costs.items()})
    fronts.to_csv(out / "fronts.csv", index=False)
    metrics.to_csv(out / "metrics.csv", index=False)
    pa.write_json(out / "validation.json", dict(identity=identity, transforms=transforms,
                  hashes={p.name: pa.digest(p) for p in (path, path.with_suffix(".json"),
                                                         out / "fronts.csv", out / "metrics.csv")}))
    return costs, fronts, metrics, transforms


def _concat(results):
    costs, fronts, metrics, transforms = {}, [], [], {}
    for c, f, m, t in results:
        costs.update(c)
        fronts += f
        metrics += m
        transforms.update(t)
    return costs, pd.concat(fronts, ignore_index=True), pd.DataFrame(metrics), transforms


def build_pooled(run=pa.DEFAULT_RUN, *, draws=DRAWS, seed=SEED):
    """CE protein pooled full tree (978 cells, 974 edges, eight rounds)."""
    from full_tree_pareto import pooled_pipeline as settings
    run = Path(run)
    ctx, fronts, _, _ = pa.build(run=run, layout_only=True, **settings.sweep_settings())
    with np.load(run / "analysis/layerwise_assignment.npz") as data:
        sweep = data["parents"]  # Hash-verified and replayed by pa.build above.
    index = tree_index()
    validation = json.loads((run / "analysis/validation.json").read_text())
    identity = _identity(target="pooled", draws=draws, seed=seed, context_key=ctx.identity["context_key"],
                         pooled_identity=validation["identity"], rounds=[len(x) for x in ctx.layers])
    permutations, groups = draw_references(ctx.names, ctx.parents, ctx.layers, index, draws, seed)
    if not np.array_equal(ctx.layers[0], ctx.leaves) or \
            {tuple(sorted(g)) for g in groups[1, "first_cousin"]} != {tuple(sorted(g)) for g in ctx.cousin_groups if len(g) > 1}:
        raise AssertionError("Round-1 first-cousin groups differ from the pooled leaf null")

    def analyze(perms):
        result = _concat([analyze_forest(ctx.names, ctx.parents, ctx.layers, ctx.travel, ctx.state, sweep, perms,
                                         groups, label="ce_protein", reference_id=f"{pa.PROFILE}:layerwise")])
        layerwise = fronts[fronts.method == pa.METHODS[0]]
        for r in range(1, len(ctx.layers) + 1):
            saved = layerwise[layerwise.scope == f"round_{r}"].sort_values("weight_index")
            mine = result[1][result[1]["round"] == r]
            np.testing.assert_allclose(mine[["travel", "cell_state", "edge_retention"]].to_numpy(),
                                       saved[["travel", "state", "retention"]].to_numpy(), rtol=1e-12, atol=1e-12)
        first = result[2].iloc[0]
        np.testing.assert_allclose([first.null_variance_travel, first.null_variance_state], ctx.scales**2, rtol=1e-10)
        return result
    return _write_or_verify(run / SUBDIR, identity, permutations, analyze)


def build_cross_species(run_id=fcs.DEFAULT_RUN, *, draws=DRAWS, seed=SEED):
    """Partial forest: six geometry x configuration sweeps over the same rounds; paired draws."""
    run = fcs.run_path(run_id)
    report = fcs.validate_run(run)
    analysis = run / "analysis"
    record = json.loads((analysis / "provenance.json").read_text())
    with np.load(analysis / "matrices.npz", allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    names, parents = record["names"], arrays["parents"]
    base = fcs.forest(names, parents, arrays["raw3d__ce_protein__travel"], arrays["raw3d__ce_protein__state"])
    index = tree_index()
    identity = _identity(target="cross_species", draws=draws, seed=seed, analysis_id=report["analysis_id"],
                         rounds=[len(x) for x in base.layers])
    permutations, groups = draw_references(names, base.parents, base.layers, index, draws, seed)

    def analyze(perms):
        results = []
        for geometry in fcs.GEOMETRIES:
            for config in fcs.PRIMARY:
                key = f"{geometry}__{config}"
                ctx = fcs.forest(names, parents, arrays[f"{key}__travel"], arrays[f"{key}__state"])
                cached = fcs.load_cache(analysis / f"{key}.npz", dict(analysis_id=report["analysis_id"], geometry=geometry,
                                        config=config, settings=record["settings"]))
                c, f, m, t = analyze_forest(names, ctx.parents, ctx.layers, ctx.travel, ctx.state, cached["parents"],
                                            perms, groups, label=key, reference_id=f"{report['analysis_id']}:{key}")
                # Round 1 + fixed internal edges must reproduce the cached full-forest first-cousin moments.
                stats = fcs.exact_null(ctx, build_ancestor_groups([names[i] for i in ctx.leaves], index, 2))
                fixed = stats[:2] - [m[0]["null_mean_travel"], m[0]["null_mean_state"]]
                internal = np.setdiff1d(ctx.edges, ctx.leaves)
                np.testing.assert_allclose(fixed, [ctx.travel[internal, ctx.parents[internal]].sum(),
                                                   ctx.state[internal, ctx.parents[internal]].sum()], rtol=1e-10)
                np.testing.assert_allclose(stats[2:], [m[0]["null_variance_travel"], m[0]["null_variance_state"],
                                                       m[0]["null_covariance"]], rtol=1e-10)
                for row in m:
                    row["geometry"], row["config"] = geometry, config
                results.append((c, [x.assign(geometry=geometry, config=config) for x in f], m, t))
        return _concat(results)
    return _write_or_verify(run / SUBDIR, identity, permutations, analyze)


def load(run, target):
    """Read-only: verify identities and hashes without solving or drawing."""
    out = Path(run) / SUBDIR
    record = json.loads((out / "validation.json").read_text())
    if record["identity"]["version"] != VERSION or record["identity"]["target"] != target:
        raise ValueError("Unsupported round-reference cache")
    for name, sha in record["hashes"].items():
        if pa.digest(out / name) != sha:
            raise ValueError(f"Changed round-reference file: {name}")
    for path, sha in record["identity"]["source_hashes"].items():
        if pa.digest(pa.ROOT / path) != sha:
            raise ValueError(f"Changed source since round references were drawn: {path}")
    with np.load(out / "draws.npz") as data:
        costs = {k[:-len("__costs")]: data[k] for k in data.files if k.endswith("__costs")}
    return dict(costs=costs, fronts=pd.read_csv(out / "fronts.csv"), metrics=pd.read_csv(out / "metrics.csv"),
                transforms=record["transforms"], identity=record["identity"],
                identity_sha=hashlib.sha256(json.dumps(record["identity"], sort_keys=True).encode()).hexdigest())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", choices=("pooled", "cross-species", "all"), default="all")
    args = parser.parse_args(argv)
    for target, build in (("pooled", build_pooled), ("cross-species", build_cross_species)):
        if args.target in (target, "all"):
            _, _, metrics, _ = build()
            columns = [c for c in ("geometry", "config", "round", "edges", "first_cousin_shuffled", "u_L", "d_LP",
                                   "d_NP", "max_edge_retention", "degenerate") if c in metrics]
            print(metrics[columns].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
