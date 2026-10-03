"""Terminal-anchored, observed-edge partial-forest species comparison.

No missing ancestor is bridged or imputed. Existing full-tree publication
caches and their hash-pinned implementation are never modified.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import cousin_references as cr
from terminal_pareto import data_loader as dl
from terminal_pareto.cross_species_analysis import endpoint_assignment, pool_matrices
from terminal_pareto.front_coordinates import EndpointTransform
from terminal_pareto.lineage_metrics import build_lineage_tree_index
from terminal_pareto.pareto_engine import build_ancestor_groups
from terminal_pareto.subtree_analysis import exact_cousin_stats

ROOT = pa.ROOT
VERSION = "terminal-anchored-partial-forest-cross-species-1"
RUN_ROOT = ROOT / "full_tree_pareto/output/runs/cross_species_layerwise_v1"
DEFAULT_RUN = "terminal_anchored_20260930"
PRIMARY = ("ce_protein", "ce_rna", "cb_rna")
GEOMETRIES = ("raw3d", "xy")
FAMILIES = ("first_cousin", "second_cousin", "third_cousin", "random_rebuild")
INTERVALS, DENSE_INTERVALS, DRAWS, SEED = 300, 1200, 10000, 42


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run_path(run_id):
    if not run_id or Path(run_id).name != run_id or run_id in (".", ".."):
        raise ValueError("Run identity must be one directory name")
    return RUN_ROOT / run_id


def ancestor_closure(seeds, available, canonical_parent):
    """Ascend independently from each seed and stop at the first unavailable node."""
    kept, heights = set(), {}
    for seed in seeds:
        if seed not in available:
            raise ValueError(f"Unavailable terminal seed: {seed}")
        current, steps, seen = seed, 0, set()
        while current in available:
            if current in seen:
                raise ValueError("Cycle in canonical ancestry")
            seen.add(current)
            kept.add(current)
            current = canonical_parent[current]
            steps += 1
        heights[seed] = steps - 1
    return kept, heights


@dataclass
class Forest:
    names: list[str]
    parents: np.ndarray
    travel: np.ndarray
    state: np.ndarray
    layers: list[np.ndarray]
    roots: np.ndarray
    leaves: np.ndarray
    internal: np.ndarray
    edges: np.ndarray
    counts: np.ndarray

    def score(self, parents):
        parents = np.asarray(parents)
        return np.stack((self.travel[self.edges, parents[..., self.edges]].sum(axis=-1),
                         self.state[self.edges, parents[..., self.edges]].sum(axis=-1)), axis=-1)

    @property
    def natural(self):
        return self.score(self.parents)


def forest(names, parents, travel, state):
    parents = np.asarray(parents, dtype=int)
    edges, roots = np.flatnonzero(parents >= 0), np.flatnonzero(parents < 0)
    if not len(edges) or np.any(parents[edges] >= edges):
        raise ValueError("Expected a nonempty canonical root-to-tip forest")
    counts = np.bincount(parents[edges], minlength=len(names))
    leaves, internal = np.flatnonzero(counts == 0), np.flatnonzero(counts > 0)
    if np.intersect1d(roots, leaves).size or np.any(counts > 2):
        raise ValueError("Isolated roots or noncanonical parent capacity")
    for matrix in (travel, state):
        if matrix.shape != (len(names), len(names)) or not np.isfinite(matrix).all():
            raise ValueError("Invalid candidate edge matrix")
        if np.any(matrix < -1e-12):
            raise ValueError("Negative edge distance")
    return Forest(list(names), parents, travel, state, pa.contraction_layers(parents),
                  roots, leaves, internal, edges, counts)


def validate_forests(ctx, parents, *, roundwise=False):
    """Replay full costs; check fixed roots, observed capacities and all cycles."""
    parents = np.asarray(parents)
    arrays = parents[None, :] if parents.ndim == 1 else parents
    n = len(ctx.names)
    if arrays.ndim != 2 or arrays.shape[1] != n or not np.issubdtype(arrays.dtype, np.integer):
        raise ValueError("Invalid parent array")
    if np.any(arrays < -1) or np.any(arrays >= n):
        raise ValueError("Parent outside cohort")
    if np.any(arrays[:, ctx.roots] != -1) or np.any(arrays[:, ctx.edges] < 0):
        raise ValueError("Changed boundary roots or scored edge scope")
    for start in range(0, len(arrays), 250):
        chunk = arrays[start:start+250]
        offsets = np.arange(len(chunk))[:, None] * n
        counts = np.bincount((chunk[:, ctx.edges] + offsets).ravel(),
                             minlength=len(chunk)*n).reshape(len(chunk), n)
        if not np.all(counts == ctx.counts):
            raise ValueError("Changed observed biological-parent capacities")
        ancestors = chunk.copy()
        for _ in range(n):
            if np.all(ancestors < 0):
                break
            ancestors = np.where(ancestors < 0, -1,
                                 np.take_along_axis(chunk, np.maximum(ancestors, 0), axis=1))
        else:
            raise ValueError("Cycle in reconstructed forest")
        if roundwise:
            for children in ctx.layers:
                if not np.all(np.sort(chunk[:, children], axis=1) == np.sort(ctx.parents[children])):
                    raise ValueError("Assignment crosses contraction-round parent slots")
    return ctx.score(parents)


def source_hashes():
    paths = ["data/cell_lineage.json", "data/protein/aggregated_all/s3_zscore.csv",
             "expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv",
             "data/c_briggsae/science.adu8249/c_elegans_tf.csv",
             "data/c_briggsae/science.adu8249/c_briggsae_tf.csv",
             "expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv",
             "data/c_briggsae/yiming/metadata.csv",
             "data/c_briggsae/science.adu8249/shared/CellTable_20241027.txt",
             "full_tree_pareto/cross_species_analysis.py", "full_tree_pareto/pooled_analysis.py",
             "full_tree_pareto/publication_analysis.py", "full_tree_pareto/cousin_references.py",
             "terminal_pareto/cross_species_analysis.py", "terminal_pareto/data_loader.py",
             "terminal_pareto/front_coordinates.py", "terminal_pareto/subtree_analysis.py",
             "terminal_pareto/pareto_engine.py", "terminal_pareto/lineage_metrics.py", "utils.py"]
    paths += [str(Path(path).relative_to(ROOT)) for _, path, _ in dl.CE_REPLICATES]
    paths += [str(Path(path).relative_to(ROOT)) for label, path, _ in dl.CB_REPLICATES_NEW if "AF16" in label]
    return {path: digest(ROOT / path) for path in paths}


def load_inputs():
    ce = pa.build_context()  # Read only; preserve its 978-node reference and caches.
    lineage = dl.load_json(ROOT / "data/cell_lineage.json")
    canonical_parent, depths = pa._canonical_maps(lineage)
    ce_rna, cb_rna, features = dl.load_ce_rna(), dl.load_cb_rna(), dl.load_rna_sel()
    metadata = pd.read_csv(ROOT / "data/c_briggsae/yiming/metadata.csv").set_index("ID")
    available = set(ce.names) & set(ce_rna.index) & set(cb_rna.index)
    cb_xyz, cb_times, tracking = {}, {}, {"ce": [], "cb": []}
    for label, path, cutoff in dl.CB_REPLICATES_NEW:
        if "AF16" not in label:
            continue
        replicate = label.split()[0]
        if metadata.loc[replicate, "Type"] != "AF16":
            raise ValueError("Expected AF16 tracking")
        xyz, valid = dl.load_briggsae_tracking(cutoff, path)
        available &= set(valid)
        # Standardize orders explicitly: the loader returns [z,x,y].
        cb_xyz[replicate] = {name: np.asarray(value, float)[[1, 2, 0]] for name, value in xyz.items()}
        table = pd.read_csv(path)
        cb_times[replicate] = table[table.time <= cutoff].groupby("cell", sort=False).time.last().to_dict()
        tracking["cb"].append(dict(label=replicate, path=str(Path(path).relative_to(ROOT)), cutoff=cutoff))
    for label, path, cutoff in dl.CE_REPLICATES:
        tracking["ce"].append(dict(label=label, path=str(Path(path).relative_to(ROOT)), cutoff=cutoff))
    biological_terminals = set(ce.nodes.loc[ce.nodes.biological_terminal, "cell"])
    seeds = [name for name in ce.names if name in biological_terminals & available
             and canonical_parent[name] in available]
    kept, heights = ancestor_closure(seeds, available, canonical_parent)
    names = [name for name in ce.names if name in kept]
    index = {name: i for i, name in enumerate(names)}
    parents = np.array([index.get(canonical_parent[name], -1) for name in names])
    edges = np.flatnonzero(parents >= 0)
    if (len(seeds), len(names), len(edges)) != (187, 485, 454):
        raise ValueError("Coverage changed: re-audit the terminal-anchored cohort")
    positions = dict(ce={label: ce.xyz[r, [ce.names.index(name) for name in names]]
                         for r, (label, _, _) in enumerate(dl.CE_REPLICATES)},
                     cb={label: np.array([xyz[name] for name in names], float) for label, xyz in cb_xyz.items()})
    expression = dict(ce_protein=ce.expression[[ce.names.index(name) for name in names]],
                      ce_rna=ce_rna.loc[names, features].to_numpy(float),
                      cb_rna=cb_rna.loc[names, features].to_numpy(float))
    if len(features) != 20 or len(set(features)) != 20:
        raise ValueError("Expected frozen shared 20 RNA features")
    states = {}
    for config, values in expression.items():
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite measured states: {config}")
        states[config] = cdist(values, values, metric="euclidean")
    matrices, totals, components = {}, {}, {}
    for geometry in GEOMETRIES:
        totals[geometry] = {}
        for species, replicates in positions.items():
            distances = {label: cdist(xyz if geometry == "raw3d" else xyz[:, :2],
                                      xyz if geometry == "raw3d" else xyz[:, :2])
                         for label, xyz in replicates.items()}
            denominators = {label: float(matrix[edges, parents[edges]].sum())
                            for label, matrix in distances.items()}
            totals[geometry][species] = denominators
            pooled = pool_matrices(distances, denominators)
            for label, matrix in distances.items():
                components[f"{geometry}__{species}__{label}"] = matrix / denominators[label]
            for config in ("cb_rna",) if species == "cb" else ("ce_protein", "ce_rna"):
                matrices[geometry, config] = (pooled, states[config])
    nodes = []
    for i, name in enumerate(names):
        row = dict(index=i, cell=name, canonical_parent=canonical_parent[name], parent_index=int(parents[i]),
                   boundary_root=parents[i] < 0, terminal_seed=name in seeds, canonical_depth=depths[name],
                   unavailable_parent=canonical_parent[name] not in available)
        ce_index = ce.names.index(name)
        for r, (label, _, cutoff) in enumerate(dl.CE_REPLICATES):
            row[f"{label}_last_time"] = int(ce.times[r, ce_index])
            row[f"{label}_at_cutoff"] = bool(ce.times[r, ce_index] == cutoff)
        for replicate, times in cb_times.items():
            row[f"{replicate}_last_time"] = int(times[name])
        nodes.append(row)
    coverage = []
    for name in ce.names:
        coverage.append(dict(cell=name, canonical_parent=canonical_parent[name],
            ce_rna_available=name in ce_rna.index, cb_rna_available=name in cb_rna.index,
            cb_tracking_available=all(name in xyz for xyz in cb_xyz.values()),
            available=name in available, retained=name in kept, terminal_seed=name in seeds,
            exclusion="" if name in kept else "measurement unavailable" if name not in available else
                      "not connected upward from a matched terminal without a gap"))
    info = dict(version=VERSION, names=names, parents=parents.tolist(), terminal_seeds=seeds,
        terminal_count=len(seeds), cohort_size=len(names), edges=len(edges), tracking=tracking,
        ancestor_steps=heights, normalization_totals=totals,
        normalization_reference="all 454 natural edges in this frozen terminal-anchored partial forest; fixed for every assignment and reference within each geometry",
        travel="equal mean of per-embryo pairwise distances divided by natural partial-forest totals; no averaged coordinates",
        coordinate_order="x,y,z; XY drops z (column 2) after explicit loader-order conversion",
        molecular="Euclidean distances, top20 z-scored protein or shared20 RNA TF stored values; no additional RNA rescaling",
        features=dict(ce_protein=ce.identity["proteins"], ce_rna=features, cb_rna=features),
        cohort_rule="start at matched biological terminal edges; ascend independently until first unavailable canonical ancestor; no skipped ancestors or imputation",
        sources=source_hashes(), numpy_version=np.__version__, scipy_version=scipy.__version__,
        limitations=["Partial measured forest, not the complete or near-complete embryonic tree",
                     "C. briggsae 3D subject to traditional embryo-tracking limitations, particularly z measurements; wording provisional",
                     "Developmental alignment and RNA measured/aggregated/imputed provenance remain unresolved",
                     "Protein/RNA differences are not a controlled modality effect",
                     "Tracking replicates share molecular matrices, not molecular replicate uncertainty",
                     "Sampled weighted optima do not enumerate the full discrete Pareto set"])
    return info, pd.DataFrame(nodes), pd.DataFrame(coverage), matrices, components, positions, expression


def exact_null(ctx, groups):
    stats = np.array(exact_cousin_stats(ctx.travel[np.ix_(ctx.leaves, ctx.parents[ctx.leaves])],
                                      ctx.state[np.ix_(ctx.leaves, ctx.parents[ctx.leaves])], groups))
    fixed = np.setdiff1d(ctx.edges, ctx.leaves)
    stats[:2] += [ctx.travel[fixed, ctx.parents[fixed]].sum(), ctx.state[fixed, ctx.parents[fixed]].sum()]
    if not np.isfinite(stats).all() or np.any(stats[2:4] <= 0):
        raise ValueError("Invalid exact full-forest first-cousin moments")
    return stats


def solve_sweep(ctx, stats, intervals):
    scales = np.sqrt(stats[2:4])
    result = np.broadcast_to(ctx.parents, (intervals+1, len(ctx.names))).copy()
    ties = []
    for round_index, children in enumerate(ctx.layers, 1):
        slots = ctx.parents[children]
        travel, state = (matrix[np.ix_(children, slots)] for matrix in (ctx.travel, ctx.state))
        state_perm, state_tie = endpoint_assignment(state, travel)
        travel_perm, travel_tie = endpoint_assignment(travel, state)
        rows = np.arange(len(children))
        for step in range(intervals+1):
            alpha = step / intervals
            chosen = (state_perm if step == 0 else travel_perm if step == intervals else
                      linear_sum_assignment(alpha*travel/scales[0] + (1-alpha)*state/scales[1])[1])
            result[step, children[rows]] = slots[chosen]
        ties.append(dict(round=round_index, edges=len(children), travel=travel_tie, state=state_tie))
    costs = validate_forests(ctx, result, roundwise=True)
    tolerance = 1e-10 * np.maximum(1, np.max(np.abs(costs), axis=0))
    if np.any(np.diff(costs[:, 0]) > tolerance[0]) or np.any(np.diff(costs[:, 1]) < -tolerance[1]):
        raise ValueError("Nonmonotone weighted objective sweep")
    return result, costs, ties


def summarize(ctx, parents, costs, stats, identity):
    transform = EndpointTransform.from_endpoints(reference_analysis_id=identity,
        travel_optimum_assignment_id=f"{identity}:{len(costs)-1}", state_optimum_assignment_id=f"{identity}:0",
        travel_optimum_costs=costs[-1], state_optimum_costs=costs[0])
    x, y = transform.transform(costs[:, 0], costs[:, 1])
    nx, ny = map(float, transform.transform(*ctx.natural))
    qx, qy = map(float, transform.transform(*stats[:2]))
    nearest = int(np.argmin(np.hypot(x-nx, y-ny)))
    arc = np.r_[0., np.cumsum(np.hypot(np.diff(x[::-1]), np.diff(y[::-1])))]
    u = (arc / arc[-1])[::-1]
    retention = (parents[:, ctx.edges] == ctx.parents[ctx.edges]).mean(axis=1)
    candidates = np.flatnonzero(retention == retention.max())
    maximum = int(candidates[np.argmin(np.abs(u[candidates]-u[nearest]))])
    frame = pd.DataFrame(dict(sweep_index=np.arange(len(costs)), alpha_travel=np.linspace(0, 1, len(costs)),
        travel=costs[:, 0], cell_state=costs[:, 1], D1=x, D2=y, u=u, edge_retention=retention,
        nearest_assignment=np.arange(len(costs)) == nearest, maximum_retention=np.arange(len(costs)) == maximum))
    metrics = dict(n=len(ctx.edges), nodes=len(ctx.names), roots=len(ctx.roots), leaves=len(ctx.leaves),
        internal=len(ctx.internal), unary=int(np.sum(ctx.counts == 1)), intervals=len(costs)-1,
        natural_travel=float(ctx.natural[0]), natural_cell_state=float(ctx.natural[1]), natural_D1=nx, natural_D2=ny,
        null_mean_travel=float(stats[0]), null_mean_state=float(stats[1]), null_D1=qx, null_D2=qy,
        null_variance_travel=float(stats[2]), null_variance_state=float(stats[3]), null_covariance=float(stats[4]),
        travel_sigma=float(np.sqrt(stats[2])), state_sigma=float(np.sqrt(stats[3])),
        nearest_index=nearest, maximum_index=maximum, u_L=float(u[nearest]),
        d_LP=float(np.hypot(x[nearest]-nx, y[nearest]-ny)), d_NP=float(np.hypot(x[nearest]-qx, y[nearest]-qy)),
        max_edge_retention=float(retention[maximum]), retention_at_nearest=float(retention[nearest]),
        retention_travel=float(retention[-1]), retention_state=float(retention[0]))
    return frame, metrics, transform


def random_rebuilds(ctx, draws, seed):
    """Random acyclic topology, fixed roots and observed 1/2-child capacities."""
    rng = np.random.default_rng(seed)
    result = np.full((draws, len(ctx.names)), -1, dtype=np.int32)
    nonroot_internal = np.setdiff1d(ctx.internal, ctx.roots)
    for parent in result:
        slots = list(np.repeat(ctx.roots, ctx.counts[ctx.roots]))
        for child in rng.permutation(nonroot_internal):
            parent[child] = slots.pop(int(rng.integers(len(slots))))
            slots.extend([int(child)] * int(ctx.counts[child]))
        if len(slots) != len(ctx.leaves):
            raise ValueError("Random rebuild changed measured capacity")
        parent[ctx.leaves] = rng.permutation(slots)
    validate_forests(ctx, result)
    return result


def save_cache(path, identity, **arrays):
    np.savez_compressed(path, **arrays)
    write_json(path.with_suffix(".json"), dict(identity=identity, sha256=digest(path)))


def load_cache(path, identity):
    record = json.loads(path.with_suffix(".json").read_text())
    if record["identity"] != identity or record["sha256"] != digest(path):
        raise ValueError(f"Changed cache identity/hash: {path}; choose a new run")
    with np.load(path, allow_pickle=False) as saved:
        return {key: saved[key] for key in saved.files}


def check_hashes(base, hashes):
    for name, expected in hashes.items():
        if digest(Path(base) / name) != expected:
            raise ValueError(f"Changed source/artifact: {name}")


def build(run, *, intervals=INTERVALS, dense_intervals=DENSE_INTERVALS, draws=DRAWS):
    if intervals <= 0 or dense_intervals <= intervals or dense_intervals % intervals or draws <= 0:
        raise ValueError("Use positive draws and a denser nested endpoint-inclusive grid")
    run, out = Path(run), Path(run) / "analysis"
    if (out / "provenance.json").exists():
        validate_run(run)
        saved_settings = json.loads((out / "provenance.json").read_text())["settings"]
        if (saved_settings["intervals"], saved_settings["dense_intervals"], saved_settings["draws"]) != (intervals, dense_intervals, draws):
            raise ValueError("Completed run has different settings; use a new run")
        return run
    info, nodes, coverage, matrices, components, positions, expression = load_inputs()
    settings = dict(intervals=intervals, dense_intervals=dense_intervals, draws=draws,
                    seeds=dict(zip(FAMILIES, (SEED, SEED+1, SEED+2, SEED+3))))
    identity = hashlib.sha256(json.dumps(dict(info=info, settings=settings), sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True, exist_ok=True)
    inputs_path = out / "inputs.json"
    if inputs_path.exists():
        previous = json.loads(inputs_path.read_text())
        if previous["analysis_id"] != identity:
            raise ValueError("Partial run has different inputs/settings; use a new run")
        check_hashes(out, previous["files"])
    else:
        nodes.to_csv(out / "nodes.csv", index=False)
        coverage.to_csv(out / "coverage.csv", index=False)
        arrays = dict(parents=np.array(info["parents"]), names=np.array(info["names"]))
        arrays.update({f"{g}__{c}__{axis}": m for (g, c), values in matrices.items()
                       for axis, m in zip(("travel", "state"), values)})
        arrays.update(components)
        arrays.update({f"coordinates__{s}__{r}": x for s, replicas in positions.items() for r, x in replicas.items()})
        arrays.update({f"expression__{c}": x for c, x in expression.items()})
        np.savez_compressed(out / "matrices.npz", **arrays)
        write_json(inputs_path, dict(info=info, settings=settings, analysis_id=identity,
            files={name: digest(out / name) for name in ("nodes.csv", "coverage.csv", "matrices.npz")}))
    ctx = forest(info["names"], info["parents"], *matrices["raw3d", "ce_protein"])
    tree_index = build_lineage_tree_index(dl.load_json(ROOT / "data/cell_lineage.json"))
    groups = {family: build_ancestor_groups([ctx.names[i] for i in ctx.leaves], tree_index, steps)
              for family, steps in zip(FAMILIES[:3], (2, 3, 4))}
    info.update(round_edges=[len(layer) for layer in ctx.layers], roots=len(ctx.roots), internal=len(ctx.internal),
                unary=int(np.sum(ctx.counts == 1)), group_sizes={f: [len(g) for g in gs] for f, gs in groups.items()})
    print(f"Partial forest: {len(ctx.names)} cells, {len(ctx.edges)} edges; rounds {info['round_edges']}", flush=True)
    shared_path = out / "references.npz"
    shared_identity = dict(analysis_id=identity, settings=settings)
    if shared_path.exists():
        shared = load_cache(shared_path, shared_identity)
    else:
        print(f"Building {draws} paired draws per reference", flush=True)
        shared = {f"{family}__permutations": cr.sample_permutations(len(ctx.leaves), groups[family], draws,
                    settings["seeds"][family]) for family in FAMILIES[:3]}
        shared["random_rebuild__parents"] = random_rebuilds(ctx, draws, settings["seeds"]["random_rebuild"])
        save_cache(shared_path, shared_identity, **shared)
    fronts, metrics, clouds, transforms, convergence, tie_records = [], [], [], {}, [], {}
    for geometry in GEOMETRIES:
        for config in PRIMARY:
            key = f"{geometry}__{config}"
            ctx = forest(info["names"], info["parents"], *matrices[geometry, config])
            stats = exact_null(ctx, groups["first_cousin"])
            path = out / f"{key}.npz"
            cache_id = dict(analysis_id=identity, geometry=geometry, config=config, settings=settings)
            if path.exists():
                cached = load_cache(path, cache_id)
                ties = json.loads(path.with_suffix(".ties.json").read_text())
            else:
                print(f"Solving {key}: {intervals+1} + {dense_intervals+1} weights", flush=True)
                parents, costs, ties = solve_sweep(ctx, stats, intervals)
                dense_parents, dense_costs, _ = solve_sweep(ctx, stats, dense_intervals)
                cached = dict(parents=parents, costs=costs, dense_parents=dense_parents, dense_costs=dense_costs)
                for family in FAMILIES[:3]:
                    cached[f"{family}__costs"] = cr.score_permutations(ctx, shared[f"{family}__permutations"])
                cached["random_rebuild__costs"] = ctx.score(shared["random_rebuild__parents"])
                save_cache(path, cache_id, **cached)
                write_json(path.with_suffix(".ties.json"), ties)
            tie_records[key] = ties
            for family in FAMILIES[:3]:
                cr.validate_permutations(shared[f"{family}__permutations"], len(ctx.leaves), groups[family])
            results = {}
            for cohort, prefix in (("base", ""), ("dense", "dense_")):
                parents, costs = cached[f"{prefix}parents"], cached[f"{prefix}costs"]
                np.testing.assert_allclose(validate_forests(ctx, parents, roundwise=True), costs, rtol=1e-12, atol=1e-9)
                frame, metric, transform = summarize(ctx, parents, costs, stats, f"{identity}:{key}")
                fronts.append(frame.assign(geometry=geometry, config=config, cohort=cohort))
                metrics.append(dict(geometry=geometry, config=config, cohort=cohort, **metric))
                results[cohort] = (frame, metric, transform)
            stride = dense_intervals // intervals
            np.testing.assert_allclose(cached["costs"], cached["dense_costs"][::stride], rtol=1e-12, atol=1e-9)
            transform = results["base"][2]
            transforms[key] = transform.metadata()
            convergence.append(dict(geometry=geometry, config=config, **{f"delta_{field}":
                results["dense"][1][field]-results["base"][1][field]
                for field in ("u_L", "d_LP", "d_NP", "max_edge_retention")}))
            for family in FAMILIES:
                costs = cached[f"{family}__costs"]
                x, y = transform.transform(costs[:, 0], costs[:, 1])
                clouds.append(pd.DataFrame(dict(geometry=geometry, config=config, family=family,
                    draw=np.arange(draws), travel=costs[:, 0], cell_state=costs[:, 1], D1=x, D2=y)))
            print(key, {f: round(results['base'][1][f], 5) for f in ("d_LP", "d_NP", "max_edge_retention")}, flush=True)
    pd.concat(fronts, ignore_index=True).to_csv(out / "fronts.csv", index=False)
    pd.DataFrame(metrics).to_csv(out / "metrics.csv", index=False)
    pd.concat(clouds, ignore_index=True).to_csv(out / "null_clouds.csv", index=False)
    pd.DataFrame(convergence).to_csv(out / "convergence.csv", index=False)
    write_json(out / "provenance.json", dict(**info, settings=settings, analysis_id=identity,
        endpoint_transforms=transforms, endpoint_ties=tie_records,
        files={p.name: digest(p) for p in out.iterdir() if p.is_file() and p.name != "provenance.json"}))
    validate_run(run)
    return run


def validate_run(run):
    out = Path(run) / "analysis"
    record = json.loads((out / "provenance.json").read_text())
    if record["version"] != VERSION:
        raise ValueError("Unsupported comparison cache")
    check_hashes(ROOT, record["sources"])
    check_hashes(out, record["files"])
    inputs = json.loads((out / "inputs.json").read_text())
    identity = hashlib.sha256(json.dumps(dict(info=inputs["info"], settings=inputs["settings"]), sort_keys=True).encode()).hexdigest()
    if identity != record["analysis_id"] or identity != inputs["analysis_id"]:
        raise ValueError("Analysis identity mismatch")
    with np.load(out / "matrices.npz", allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    np.testing.assert_array_equal(arrays["names"], record["names"])
    np.testing.assert_array_equal(arrays["parents"], record["parents"])
    canonical, depths = pa._canonical_maps(dl.load_json(ROOT / "data/cell_lineage.json"))
    for i, p in enumerate(arrays["parents"]):
        if p >= 0 and (canonical[record["names"][i]] != record["names"][p] or
                       depths[record["names"][i]]-depths[record["names"][p]] != 1):
            raise ValueError("Imputed or skipped canonical edge")
    settings = record["settings"]
    shared = load_cache(out / "references.npz", dict(analysis_id=identity, settings=settings))
    tree_index = build_lineage_tree_index(dl.load_json(ROOT / "data/cell_lineage.json"))
    fronts, metrics = pd.read_csv(out / "fronts.csv"), pd.read_csv(out / "metrics.csv")
    clouds = pd.read_csv(out / "null_clouds.csv")
    verified, null_draws = 0, 0
    for geometry in GEOMETRIES:
        np.testing.assert_array_equal(arrays[f"{geometry}__ce_protein__travel"], arrays[f"{geometry}__ce_rna__travel"])
        for config in PRIMARY:
            key = f"{geometry}__{config}"
            ctx = forest(record["names"], arrays["parents"], arrays[f"{key}__travel"], arrays[f"{key}__state"])
            species = "cb" if config == "cb_rna" else "ce"
            components = []
            for replicate in record["tracking"][species]:
                label = replicate["label"]
                xyz = arrays[f"coordinates__{species}__{label}"]
                coords = xyz if geometry == "raw3d" else xyz[:, :2]
                raw = cdist(coords, coords)
                denom = record["normalization_totals"][geometry][species][label]
                np.testing.assert_allclose(raw[ctx.edges, ctx.parents[ctx.edges]].sum(), denom, rtol=1e-12)
                np.testing.assert_allclose(raw/denom, arrays[f"{geometry}__{species}__{label}"], rtol=1e-12, atol=1e-14)
                components.append(raw/denom)
            np.testing.assert_allclose(np.mean(components, axis=0), ctx.travel, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(cdist(arrays[f"expression__{config}"], arrays[f"expression__{config}"]), ctx.state, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(ctx.natural[0], 1, rtol=1e-12)
            groups = {f: build_ancestor_groups([ctx.names[i] for i in ctx.leaves], tree_index, step)
                      for f, step in zip(FAMILIES[:3], (2, 3, 4))}
            stats = exact_null(ctx, groups["first_cousin"])
            cached = load_cache(out / f"{key}.npz", dict(analysis_id=identity, geometry=geometry, config=config, settings=settings))
            for cohort, prefix, intervals in (("base", "", settings["intervals"]), ("dense", "dense_", settings["dense_intervals"])):
                parents, costs = cached[f"{prefix}parents"], cached[f"{prefix}costs"]
                if parents.shape != (intervals+1, len(ctx.names)):
                    raise ValueError("Incomplete assignment sweep")
                np.testing.assert_allclose(validate_forests(ctx, parents, roundwise=True), costs, rtol=1e-12, atol=1e-9)
                frame, metric, transform = summarize(ctx, parents, costs, stats, f"{identity}:{key}")
                saved = fronts[(fronts.geometry == geometry) & (fronts.config == config) & (fronts.cohort == cohort)]
                np.testing.assert_allclose(frame.select_dtypes(include="number"), saved[frame.select_dtypes(include="number").columns], rtol=1e-11, atol=1e-9)
                m = metrics[(metrics.geometry == geometry) & (metrics.config == config) & (metrics.cohort == cohort)].iloc[0]
                np.testing.assert_allclose([metric[k] for k in metric], [m[k] for k in metric], rtol=1e-11, atol=1e-9)
                if cohort == "base":
                    if transform.metadata() != record["endpoint_transforms"][key]:
                        raise ValueError("Changed display transform")
                verified += len(parents)
            transform = EndpointTransform(**{k: v for k, v in record["endpoint_transforms"][key].items()
                                              if k not in ("display_mode", "clipping", "axis_limits")})
            for family in FAMILIES:
                if family == "random_rebuild":
                    if verified == settings["intervals"]+settings["dense_intervals"]+2:
                        validate_forests(ctx, shared["random_rebuild__parents"])
                    replay = ctx.score(shared["random_rebuild__parents"])
                else:
                    perms = shared[f"{family}__permutations"]
                    if perms.shape != (settings["draws"], len(ctx.leaves)):
                        raise ValueError("Incomplete cousin reference")
                    cr.validate_permutations(perms, len(ctx.leaves), groups[family])
                    replay = cr.score_permutations(ctx, perms)
                np.testing.assert_allclose(replay, cached[f"{family}__costs"], rtol=1e-12, atol=1e-9)
                saved = clouds[(clouds.geometry == geometry) & (clouds.config == config) & (clouds.family == family)]
                if len(saved) != settings["draws"]:
                    raise ValueError("Incomplete display/reference scores")
                np.testing.assert_allclose(replay, saved[["travel", "cell_state"]], rtol=1e-12, atol=1e-9)
                x, y = transform.transform(replay[:, 0], replay[:, 1])
                np.testing.assert_allclose(np.column_stack((x, y)), saved[["D1", "D2"]], rtol=1e-11, atol=1e-12)
                null_draws += len(replay)
    return dict(analysis_id=identity, nodes=record["cohort_size"], edges=record["edges"],
                terminal_seeds=record["terminal_count"], verified_assignments=verified,
                verified_reference_scores=null_draws, verified_source_hashes=len(record["sources"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=DEFAULT_RUN)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    run = run_path(args.run_id)
    if args.verify_only or args.render_only:
        print(json.dumps(validate_run(run), indent=2), flush=True)
    else:
        build(run)
    if not args.verify_only:
        from full_tree_pareto.fig_cross_species import render
        render(run)
        print(f"Partial-forest comparison ready: {run / 'figures'}", flush=True)


if __name__ == "__main__":
    main()
