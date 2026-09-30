"""Pooled full-tree analysis, isolated from the frozen embryo-1 pipeline.

All node indices are canonical breadth-first cohort indices. Caches store
complete parent arrays, not just scalar costs, for independent score replay.
"""
from __future__ import annotations

from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import hashlib
import heapq
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from full_tree_pareto.publication_analysis import _canonical_maps, _nondominated_mask
from terminal_pareto.subtree_analysis import exact_cousin_stats
from utils import load_json, lineage_name_mapping

ROOT = Path(__file__).resolve().parents[1]
PROFILE = "pooled_full_tree_v1"
DEFAULT_RUN = ROOT / "full_tree_pareto/output/runs" / PROFILE / "terminal_clamped_20260927"
CUTOFFS = (255, 247, 225)
REFERENCE = "Terminal-clamped Gaussian reference"
METHODS = ("Layerwise assignment", "Degree-constrained spanning forest",
           "Top-down rebuild", "Bottom-up by layer", "Paired bottom-up")
NULLS = ("First-cousin shuffle", "Internal-layer shuffle", "Full assignment shuffle",
         "Random rebuild", REFERENCE)
VERSION = "pooled-full-tree-1"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


@dataclass
class Context:
    names: list[str]
    parents: np.ndarray
    depths: np.ndarray
    xyz: np.ndarray
    times: np.ndarray
    expression: np.ndarray
    travel: np.ndarray
    state: np.ndarray
    denominators: np.ndarray
    roots: np.ndarray
    internal: np.ndarray
    leaves: np.ndarray
    edges: np.ndarray
    layers: list[np.ndarray]
    cousin_groups: list[np.ndarray]
    scales: np.ndarray
    identity: dict
    nodes: pd.DataFrame

    def score(self, parents):
        parents = np.asarray(parents)
        children = np.flatnonzero(parents >= 0)
        return np.array([self.travel[children, parents[children]].sum(),
                         self.state[children, parents[children]].sum()])

    @property
    def natural(self):
        return self.score(self.parents)

    def combined(self, alpha):
        return alpha * self.travel / self.scales[0] + (1-alpha) * self.state / self.scales[1]


def contraction_layers(parents):
    """Score each child once; promote a parent only after all daughters finish."""
    children = [np.flatnonzero(parents == i) for i in range(len(parents))]
    roots = set(np.flatnonzero(parents < 0))
    current = [i for i, ch in enumerate(children) if not len(ch) and i not in roots]
    done = set()
    layers = []
    while current:
        layers.append(np.array(sorted(current), dtype=int))
        done.update(current)
        current = sorted({int(parents[i]) for i in current} - roots - done)
        current = [i for i in current if set(children[i]).issubset(done)]
    flat = np.concatenate(layers)
    if set(flat) != set(np.flatnonzero(parents >= 0)) or len(flat) != len(set(flat)):
        raise AssertionError("Contraction rounds do not partition edges")
    return layers


def build_context():
    lineage_path = ROOT / "data/cell_lineage.json"
    expression_path = ROOT / "data/protein/aggregated_all/s3_zscore.csv"
    selected_path = ROOT / "expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv"
    lineage = load_json(lineage_path)
    canonical_parent, depth = _canonical_maps(lineage)
    expression = pd.read_csv(expression_path, index_col=0).T
    selected = pd.read_csv(selected_path)["protein"].tolist()
    expression = expression[selected].fillna(0)
    xyz_maps, time_maps, sources = [], [], [lineage_path, expression_path, selected_path]
    for embryo, cutoff in enumerate(CUTOFFS, 1):
        path = ROOT / f"data/embryo{embryo}/tracks.txt"
        sources.append(path)
        tracks = pd.read_csv(path, sep="\t")
        tracks = tracks[tracks.t <= cutoff]
        xyz, times = {}, {}
        for name, rows in tracks.groupby("name", sort=False):
            if name not in expression.index or (len(rows) == 1 and rows.t.iloc[0] == cutoff):
                continue
            xyz[name] = rows[["x", "y", "z"]].iloc[-1].to_numpy(float) * .1625
            times[name] = int(rows.t.iloc[-1])
        xyz_maps.append(xyz)
        time_maps.append(times)
    shared = set.intersection(*(set(x) for x in xyz_maps))
    names, queue = [], deque([lineage])
    while queue:
        node = queue.popleft()
        name = lineage_name_mapping(node["did"])
        if name in shared:
            names.append(name)
        queue.extend(node.get("children", []))
    index = {n: i for i, n in enumerate(names)}
    parents = np.array([index.get(canonical_parent[n], -1) for n in names])
    roots = np.flatnonzero(parents < 0)
    edges = np.flatnonzero(parents >= 0)
    counts = np.bincount(parents[edges], minlength=len(names))
    leaves, internal = np.flatnonzero(counts == 0), np.flatnonzero(counts > 0)
    if set(np.array(names)[roots]) != {"ABa", "ABp", "EMS", "P2"}:
        raise ValueError("Shared cohort is not an ancestor-complete four-root forest")
    if not np.all(counts[internal] == 2) or np.any(parents[edges] >= edges):
        raise ValueError("Expected a complete binary forest in root-to-tip order")
    xyz = np.array([[mapping[n] for n in names] for mapping in xyz_maps])
    times = np.array([[mapping[n] for n in names] for mapping in time_maps])
    if np.any(times[:, edges] - times[:, parents[edges]] <= 0):
        raise ValueError("Nonpositive represented spatial clock")
    components = np.array([cdist(x, x) for x in xyz])
    denominators = components[:, edges, parents[edges]].sum(axis=1)
    travel = (components / denominators[:, None, None]).mean(axis=0)
    protein = expression.loc[names].to_numpy(float)
    state = cdist(protein, protein)
    grouped = defaultdict(list)
    for j, leaf in enumerate(leaves):
        grouped[canonical_parent[canonical_parent[names[leaf]]]].append(j)
    groups = [np.array(g) for g in grouped.values()]
    stats = exact_cousin_stats(travel[np.ix_(leaves, parents[leaves])],
                               state[np.ix_(leaves, parents[leaves])], [g for g in groups if len(g) > 1])
    scales = np.sqrt(stats[2:4])
    if not np.isfinite(scales).all() or np.any(scales <= 0):
        raise ValueError("Nonpositive exact cousin variance")
    layers = contraction_layers(parents)
    round_map = {int(c): r+1 for r, layer in enumerate(layers) for c in layer}
    rows = []
    for i, n in enumerate(names):
        row = dict(index=i, cell=n, parent=canonical_parent[n], parent_index=int(parents[i]),
                   canonical_depth=depth[n], root=i in roots, terminal=i in leaves,
                   biological_terminal=not any(p == n for p in canonical_parent.values()),
                   contraction_round=round_map.get(i, 0))
        for r, cutoff in enumerate(CUTOFFS):
            row[f"embryo{r+1}_last_time"] = int(times[r, i])
            row[f"embryo{r+1}_at_cutoff"] = bool(times[r, i] == cutoff)
        rows.append(row)
    identity = dict(version=VERSION, profile=PROFILE, cutoffs=CUTOFFS,
                    names=names, parents=parents.tolist(), proteins=selected,
                    source_hashes={str(p.relative_to(ROOT)): digest(p) for p in sources},
                    travel="equal mean of per-embryo pairwise distances / full-cohort natural totals",
                    molecular="Euclidean top20 z-scored protein",
                    optimizer_scaling="global exact first-cousin null SDs, shared by every round and heuristic",
                    denominators=denominators.tolist(), scales=scales.tolist(),
                    reference="independent spatial processes per embryo; shared protein process; forward internal states then terminal clamp")
    identity["context_key"] = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    ctx = Context(names, parents, np.array([depth[n] for n in names]), xyz, times, protein,
                  travel, state, denominators, roots, internal, leaves, edges, layers,
                  groups, scales, identity, pd.DataFrame(rows))
    if not np.isclose(ctx.natural[0], 1):
        raise AssertionError("Pooled natural travel must equal one")
    return ctx


def validate_assignment(ctx, parent, *, fixed_roots=True):
    parent = np.asarray(parent)
    n = len(ctx.names)
    if parent.shape != (n,) or not np.issubdtype(parent.dtype, np.integer):
        raise ValueError("Invalid parent-array shape/type")
    if np.any(parent < -1) or np.any(parent >= n):
        raise ValueError("Parent outside cohort")
    edges, roots = np.flatnonzero(parent >= 0), np.flatnonzero(parent < 0)
    if len(edges) != len(ctx.edges) or len(roots) != len(ctx.roots):
        raise ValueError("Incorrect edge/component scope")
    if fixed_roots and set(roots) != set(ctx.roots):
        raise ValueError("Changed fixed roots")
    counts = np.bincount(parent[edges], minlength=n)
    if np.any(counts[ctx.leaves]) or not np.all(counts[ctx.internal] == 2):
        raise ValueError("Internal/terminal identities or binary capacity changed")
    for start in range(n):
        seen, node = set(), start
        while node >= 0:
            if node in seen:
                raise ValueError("Cycle in reconstructed forest")
            seen.add(node)
            node = int(parent[node])
    return ctx.score(parent)


def layerwise(ctx, cost):
    parent = ctx.parents.copy()
    for children in ctx.layers:
        slots = ctx.parents[children]
        r, c = linear_sum_assignment(cost[np.ix_(children, slots)])
        parent[children[r]] = slots[c]
    return parent


def finish_terminals(ctx, cost, parent):
    assigned = parent[parent >= 0]
    count = np.bincount(assigned, minlength=len(parent))
    slots = np.repeat(ctx.internal, 2-count[ctx.internal])
    if len(slots) != len(ctx.leaves):
        raise AssertionError("Wrong terminal capacity")
    r, c = linear_sum_assignment(cost[np.ix_(ctx.leaves, slots)])
    parent[ctx.leaves[r]] = slots[c]
    return parent


def spanning(ctx, cost):
    internal = ctx.internal
    roots = set(ctx.roots)
    u, v = np.triu_indices(len(internal), 1)
    order = np.argsort(cost[internal[u], internal[v]], kind="stable")
    uf = np.arange(len(internal))
    nr = np.array([int(i in roots) for i in internal])
    degrees = np.zeros(len(internal), dtype=int)
    graph = nx.Graph()
    graph.add_nodes_from(internal.tolist())
    def find(i):
        while uf[i] != i:
            uf[i] = uf[uf[i]]
            i = uf[i]
        return i
    for j in order:
        a, b = int(u[j]), int(v[j])
        if degrees[a] >= (2 if internal[a] in roots else 3) or degrees[b] >= (2 if internal[b] in roots else 3):
            continue
        ra, rb = find(a), find(b)
        if ra == rb or nr[ra]+nr[rb] > 1:
            continue
        uf[rb] = ra
        nr[ra] += nr[rb]
        degrees[a] += 1
        degrees[b] += 1
        graph.add_edge(int(internal[a]), int(internal[b]))
        if graph.number_of_edges() == len(internal)-len(roots):
            break
    if nx.number_connected_components(graph) != len(roots):
        raise AssertionError("Incomplete spanning forest")
    parent = np.full(len(ctx.names), -1, dtype=int)
    for root in sorted(roots):
        for a, b in nx.bfs_edges(graph, root):
            parent[b] = a
    return finish_terminals(ctx, cost, parent)


def top_down(ctx, cost):
    remaining = set(ctx.internal)-set(ctx.roots)
    queue = [(cost[p, c], int(p), int(c)) for p in ctx.roots for c in sorted(remaining)]
    heapq.heapify(queue)
    parent = np.full(len(ctx.names), -1, dtype=int)
    counts = np.zeros(len(ctx.names), dtype=int)
    while remaining:
        _, p, c = heapq.heappop(queue)
        if c not in remaining or counts[p] >= 2:
            continue
        parent[c] = p
        counts[p] += 1
        remaining.remove(c)
        for child in sorted(remaining):
            heapq.heappush(queue, (cost[c, child], c, child))
    return finish_terminals(ctx, cost, parent)


def bottom_up(ctx, cost):
    """Fixed topology, internal identities assigned bottom-up; roots may move."""
    mapping = np.arange(len(ctx.names))
    pool = list(ctx.internal)
    for depth in sorted(set(ctx.depths), reverse=True):
        slots = np.array([i for i in ctx.internal if ctx.depths[i] == depth], dtype=int)
        if not len(slots):
            continue
        columns = [cost[np.ix_(pool, mapping[np.flatnonzero(ctx.parents == p)])].sum(axis=1)
                   for p in slots]
        r, c = linear_sum_assignment(np.column_stack(columns))
        mapping[slots[c]] = np.array(pool)[r]
        used = set(np.array(pool)[r])
        pool = [i for i in pool if i not in used]
    parent = np.full(len(ctx.names), -1, dtype=int)
    parent[mapping[ctx.edges]] = mapping[ctx.parents[ctx.edges]]
    return parent


def paired(ctx, cost):
    """Minimum-weight blossom pairing then Hungarian parent assignment.

    Preserve every unmatched/unselected component at each round. Stop when all
    internal states are used; four roots emerge, with identities not fixed.
    """
    pool, bottom = list(ctx.internal), list(ctx.leaves)
    parent = np.full(len(ctx.names), -1, dtype=int)
    while pool:
        graph = nx.Graph()
        graph.add_nodes_from(bottom)
        graph.add_weighted_edges_from((a, b, float(cost[a, b]))
                                     for j, a in enumerate(bottom) for b in bottom[j+1:])
        pairs = sorted(tuple(sorted(p)) for p in nx.min_weight_matching(graph))
        if not pairs:
            raise AssertionError("Pairing stalled with unused internal states")
        pcost = cost[np.ix_(pool, [p[0] for p in pairs])] + cost[np.ix_(pool, [p[1] for p in pairs])]
        r, c = linear_sum_assignment(pcost)
        used_children, used_parents = set(), set()
        for i, j in zip(r, c):
            p = pool[i]
            a, b = pairs[j]
            parent[a] = parent[b] = p
            used_children.update((a, b))
            used_parents.add(p)
        bottom = sorted((set(bottom)-used_children) | used_parents)
        pool = [p for p in pool if p not in used_parents]
    return parent


SOLVERS = dict(zip(METHODS, (layerwise, spanning, top_down, bottom_up, paired)))


def _solve(payload):
    ctx, method, index, intervals = payload
    parent = SOLVERS[method](ctx, ctx.combined(index/intervals))
    validate_assignment(ctx, parent, fixed_roots=method in METHODS[:3])
    return parent


def random_nulls(ctx, draws=10000, seed=42):
    results = {}
    n = len(ctx.names)
    for model_index, model in enumerate(NULLS[:4]):
        rng = np.random.default_rng(seed+model_index)
        costs = np.empty((draws, 2))
        for draw in range(draws):
            mapping = np.arange(n)
            if model_index == 0:
                for group in ctx.cousin_groups:
                    ids = ctx.leaves[group]
                    mapping[ids] = rng.permutation(ids)
            elif model_index == 1:
                for depth in sorted(set(ctx.depths[ctx.internal])):
                    ids = ctx.internal[ctx.depths[ctx.internal] == depth]
                    mapping[ids] = rng.permutation(ids)
            elif model_index == 2:
                mapping[ctx.internal] = rng.permutation(ctx.internal)
                mapping[ctx.leaves] = rng.permutation(ctx.leaves)
            if model_index < 3:
                parent = np.full(n, -1, dtype=int)
                parent[mapping[ctx.edges]] = mapping[ctx.parents[ctx.edges]]
            else:
                parent = np.full(n, -1, dtype=int)
                slots = list(np.repeat(ctx.roots, 2))
                for child in rng.permutation(np.setdiff1d(ctx.internal, ctx.roots)):
                    parent[child] = slots.pop(int(rng.integers(len(slots))))
                    slots.extend((child, child))
                parent[ctx.leaves] = rng.permutation(slots)
            # Score the full mapped topology. Never call the legacy coordinate
            # scorer, which would silently revert to embryo 1.
            costs[draw] = ctx.score(parent)
            if draw < 3:
                validate_assignment(ctx, parent, fixed_roots=model_index in (0, 3))
        results[model] = costs
    return results


def forward_clamped(values, parents, clock, leaves, draws, rng, *, return_states=False):
    """Free Gaussian internal states, then observed leaves; NOT conditioning.

    Omitting unused leaf innovations is distributionally identical to forward
    simulating all nodes and overwriting the leaves afterwards.
    """
    edges = np.flatnonzero(parents >= 0)
    roots = np.flatnonzero(parents < 0)
    delta = values[edges]-values[parents[edges]]
    scaled = delta/np.sqrt(clock[edges, None])
    covariance = scaled.T @ scaled / len(edges)
    factor = np.linalg.cholesky(covariance)
    states = np.zeros((draws, len(parents), values.shape[1]))
    states[:, roots] = values[roots]
    terminal = set(leaves)
    for child in edges:
        if child in terminal:
            states[:, child] = values[child]
        else:
            states[:, child] = states[:, parents[child]] + np.sqrt(clock[child]) * (rng.normal(size=(draws, values.shape[1])) @ factor.T)
    if not np.array_equal(states[:, leaves], np.broadcast_to(values[leaves], states[:, leaves].shape)):
        raise AssertionError("Terminal clamping failed")
    if not np.array_equal(states[:, roots], np.broadcast_to(values[roots], states[:, roots].shape)):
        raise AssertionError("Root states changed")
    edge_norm = np.linalg.norm(states[:, edges]-states[:, parents[edges]], axis=2)
    return edge_norm.sum(axis=1), covariance, states if return_states else None


def gaussian(ctx, draws=10000, seed=242):
    rngs = [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(4)]
    clocks = []
    for times in ctx.times:
        clock = np.ones(len(ctx.names))
        clock[ctx.edges] = times[ctx.edges]-times[ctx.parents[ctx.edges]]
        clocks.append(clock)
    clocks.append(np.ones(len(ctx.names)))
    values = list(ctx.xyz) + [ctx.expression]
    components = np.empty((draws, 4))
    covariances = {}
    for start in range(0, draws, 100):
        stop = min(start+100, draws)
        for r in range(4):
            totals, cov, _ = forward_clamped(values[r], ctx.parents, clocks[r], ctx.leaves, stop-start, rngs[r])
            components[start:stop, r] = totals
            covariances[f"block_{r}"] = cov
    costs = np.column_stack(((components[:, :3]/ctx.denominators).mean(axis=1), components[:, 3]))
    return costs, components, covariances


def cache_identity(ctx, intervals, draws, seed):
    return dict(context_key=ctx.identity["context_key"], version=VERSION,
                intervals=intervals, draws=draws, seed=seed,
                analysis_source_sha256=digest(Path(__file__)),
                dependency_hashes={str(p.relative_to(ROOT)): digest(p) for p in (
                    ROOT / "terminal_pareto/subtree_analysis.py",
                    ROOT / "full_tree_pareto/publication_analysis.py", ROOT / "utils.py")},
                numpy_version=np.__version__, networkx_version=nx.__version__)


def verify_cache(path, identity):
    record = json.loads(path.with_suffix(".json").read_text())
    if record["identity"] != identity or record["sha256"] != digest(path):
        raise ValueError(f"Cache identity/hash mismatch: {path}; use a new run, not stale caches")


def save_cache(path, identity, **arrays):
    np.savez_compressed(path, **arrays)
    write_json(path.with_suffix(".json"), dict(identity=identity, sha256=digest(path)))


def build(run=DEFAULT_RUN, intervals=300, draws=10000, seed=42, workers=4, layout_only=False):
    ctx = build_context()
    run = Path(run)
    out = run / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    identity = cache_identity(ctx, intervals, draws, seed)
    context_path = out / "context.json"
    if context_path.exists():
        if json.loads(context_path.read_text()) != json.loads(json.dumps(ctx.identity)):
            raise ValueError("Context changed: select a new run directory")
    elif layout_only:
        raise FileNotFoundError("Layout-only requires a complete existing run")
    else:
        write_json(context_path, ctx.identity)
        ctx.nodes.to_csv(out / "nodes.csv", index=False)
        np.savez_compressed(out / "matrices.npz", travel=ctx.travel, state=ctx.state,
                            coordinates=ctx.xyz, expression=ctx.expression, parents=ctx.parents)
    print(f"Cohort: {len(ctx.names)} cells, {len(ctx.edges)} edges; rounds {[len(x) for x in ctx.layers]}", flush=True)
    assignments = {}
    for method in METHODS:
        path = out / (method.lower().replace(" ", "_")+".npz")
        if path.exists():
            verify_cache(path, identity)
            with np.load(path) as data:
                parent = data["parents"]
                costs = data["costs"]
        else:
            if layout_only:
                raise FileNotFoundError(path)
            print(f"Computing {method}: {intervals+1} weights", flush=True)
            payloads = ((ctx, method, i, intervals) for i in range(intervals+1))
            if workers > 1:
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    parent = np.array(list(pool.map(_solve, payloads, chunksize=1)))
            else:
                parent = np.array([_solve(p) for p in payloads])
            costs = np.array([ctx.score(p) for p in parent])
            save_cache(path, identity, parents=parent, costs=costs)
        if parent.shape != (intervals+1, len(ctx.names)):
            raise ValueError("Incomplete assignment sweep")
        replay = np.array([validate_assignment(ctx, p, fixed_roots=method in METHODS[:3]) for p in parent])
        np.testing.assert_allclose(replay, costs, atol=1e-12, rtol=1e-12)
        assignments[method] = parent
        print(f"Validated {method}", flush=True)
    path = out / "nulls.npz"
    if path.exists():
        verify_cache(path, identity)
        with np.load(path) as data:
            nulls = {name: data[f"null_{i}"] for i, name in enumerate(NULLS)}
    else:
        if layout_only:
            raise FileNotFoundError(path)
        print(f"Computing nulls and terminal-clamped Gaussian: {draws} draws each", flush=True)
        nulls = random_nulls(ctx, draws, seed)
        nulls[REFERENCE], components, covariances = gaussian(ctx, draws, seed+200)
        save_cache(path, identity, **{f"null_{i}": nulls[name] for i, name in enumerate(NULLS)},
                   gaussian_components=components, **covariances)
    for values in nulls.values():
        if values.shape != (draws, 2) or not np.isfinite(values).all():
            raise ValueError("Invalid null cache")
    fronts, layer_rows = [], []
    for method, parent_arrays in assignments.items():
        scopes = [("aggregate", ctx.edges)]
        if method == METHODS[0]:
            scopes += [(f"round_{i+1}", layer) for i, layer in enumerate(ctx.layers)]
        for scope, children in scopes:
            xy = np.array([[ctx.travel[children, p[children]].sum(), ctx.state[children, p[children]].sum()]
                           for p in parent_arrays]) if scope != "aggregate" else np.array([ctx.score(p) for p in parent_arrays])
            mask = _nondominated_mask(xy[:, 0], xy[:, 1])
            for i, (p, costs) in enumerate(zip(parent_arrays, xy)):
                fronts.append(dict(method=method, scope=scope, weight_index=i, alpha_travel=i/intervals,
                                   travel=costs[0], state=costs[1], nondominated=bool(mask[i]),
                                   retention=float(np.mean(p[ctx.edges] == ctx.parents[ctx.edges])) if scope == "aggregate"
                                   else float(np.mean(p[children] == ctx.parents[children]))))
            if method == METHODS[0]:
                natural = [ctx.travel[children, ctx.parents[children]].sum(), ctx.state[children, ctx.parents[children]].sum()]
                layer_rows.append(dict(scope=scope, edges=len(children), natural_travel=natural[0], natural_state=natural[1]))
    fronts = pd.DataFrame(fronts)
    layers = pd.DataFrame(layer_rows)
    summaries = pd.DataFrame([dict(model=name, travel_mean=x[:, 0].mean(), state_mean=x[:, 1].mean(),
                                  travel_sd=x[:, 0].std(), state_sd=x[:, 1].std(), draws=len(x)) for name, x in nulls.items()])
    if not layout_only:
        fronts.to_csv(out / "fronts.csv", index=False)
        layers.to_csv(out / "layers.csv", index=False)
        summaries.to_csv(out / "null_summary.csv", index=False)
        write_json(out / "validation.json", dict(identity=identity, assignments_replayed=(intervals+1)*len(METHODS),
                   nodes=len(ctx.names), edges=len(ctx.edges), leaves=len(ctx.leaves), roots=len(ctx.roots),
                   rounds=[len(x) for x in ctx.layers], natural=ctx.natural.tolist(),
                   hashes={p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "validation.json"}))
    else:
        record = json.loads((out / "validation.json").read_text())
        if record["identity"] != identity:
            raise ValueError("Validation identity mismatch")
        for name, sha in record["hashes"].items():
            if digest(out / name) != sha:
                raise ValueError(f"Analysis file changed: {name}")
    return ctx, fronts, layers, nulls
