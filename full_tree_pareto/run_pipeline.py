"""Standalone runner for the full-tree Pareto pipeline."""
import json, sys, numpy as np, pandas as pd, os
from pathlib import Path
from collections import defaultdict, deque

MODULE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = MODULE_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))
from utils import lineage_name_mapping as map_names, load_json
from pareto_core import LineageTree, LineageOptimization

OUT = str(MODULE_DIR / "output" / "internal_opt" / "celegans_prot" / "embryo1")
NULL_DIR = os.path.join(OUT, 'null')
HEUR_DIR = os.path.join(OUT, 'heuristics')
os.makedirs(NULL_DIR, exist_ok=True)
os.makedirs(HEUR_DIR, exist_ok=True)

RANDOM_ITER = 100_000
MAX_WORKERS = 24

# ═══ DATA LOADING (exact old pipeline from pareto_core.ipynb) ═══
print("=" * 60)
print("DATA LOADING")
print("=" * 60)

lineage_data = load_json(PROJECT_ROOT / 'data/cell_lineage.json')
lineage_exp_df = pd.read_csv(PROJECT_ROOT / "data/protein/aggregated_all/s3_pca_10d.csv", index_col=0)
lineage_exp_df = lineage_exp_df.fillna(0)
print(f"Expression: {lineage_exp_df.shape}")

tracking_df_all = pd.read_csv(PROJECT_ROOT / "data/embryo1/tracks.txt", sep="\t")
tracking_df = tracking_df_all[tracking_df_all["t"] <= 255]
S = 0.1625

# DFS
terminal_nodes = []; internal_nodes = []; parent_dict = {}
def dfs(node, parent, ancestors=[]):
    children = node.get("children", [])
    lookup_name = map_names(node["did"])
    if len(children) == 0:
        p_lookup_name = map_names(parent['did'])
        parent_dict[lookup_name] = p_lookup_name
        terminal_nodes.append(lookup_name)
    else:
        internal_nodes.append(lookup_name)
        if parent is not None:
            parent_dict[lookup_name] = map_names(parent['did'])
        for child in children:
            dfs(child, node, ancestors + [lookup_name])
dfs(lineage_data, None)
print(f"Terminals: {len(terminal_nodes)}, Internal: {len(internal_nodes)}")

# XYZ and tracking times
cell_names = tracking_df["name"].unique()
lineage_xyz_df = pd.DataFrame(columns=["x", "y", "z"])
lineage_name_to_last_tracking_time = {}
for name in cell_names:
    if name not in lineage_exp_df.index: continue
    lt = tracking_df[tracking_df['name'] == name]
    tp = lt["t"].values
    if len(tp) == 1 and tp[0] == 255: continue
    lineage_name_to_last_tracking_time[name] = int(tp[-1])
    lineage_xyz_df.loc[name] = lt[["x", "y", "z"]].values[-1] * S
lineage_name_to_last_tracking_time["AB"] = -3
lineage_name_to_last_tracking_time["P1"] = -2
lineage_name_to_last_tracking_time["P0"] = -20

# Untracked nodes + common_lineages (OLD PIPELINE: only 3 hardcoded)
untracked_nodes = ["P0", "AB", "P1"]
common_lineages = lineage_xyz_df.index
for name in untracked_nodes:
    lineage_xyz_df.loc[name] = np.repeat(np.nan, 3)
    lineage_exp_df.loc[name] = np.repeat(np.nan, lineage_exp_df.shape[1])
common_lineages = untracked_nodes + common_lineages.tolist()
print(f"Common lineages: {len(common_lineages)}")

lt_list = [-1] * len(common_lineages)
for idx, name in enumerate(common_lineages):
    lt_list[idx] = lineage_name_to_last_tracking_time.get(name, -1)

# Build tree
lineage_tree = LineageTree()
lineage_tree.add_node(0, -1); lineage_tree.root = 0
lineage_tree.add_node(1, 0); lineage_tree.add_node(2, 0)
bq = deque(); bq.append((lineage_data, -1))
while bq:
    node, pidx = bq.popleft()
    cn = map_names(node['did'])
    if cn in common_lineages:
        ci = common_lineages.index(cn)
        if ci >= 3:
            ptid = lineage_tree.reverse_lineage_id_mapping[pidx]
            lineage_tree.add_node(ci, ptid)
        for child in node.get('children', []):
            bq.append((child, ci))
    else:
        for child in node.get('children', []):
            bq.append((child, pidx))
lineage_tree.record_last_tracked_time(lt_list)
print(f"Tree size: {lineage_tree.size}")

lxm = np.nan_to_num(lineage_xyz_df.loc[common_lineages].values, nan=0.0)
lem = np.nan_to_num(lineage_exp_df.loc[common_lineages].values.astype(float), nan=0.0)

# Optimization
# Compute first_layer dynamically from tree topology
first_layer = []
for pid in [1, 2]:
    for cid in lineage_tree.children_list[pid]:
        first_layer.append((cid, 2))
print(f"first_internal_layer: {first_layer}")
tids = [lineage_tree.lineage_id_mapping[t] for t in range(lineage_tree.size) if len(lineage_tree.children_list[t]) == 0]
ccd = {common_lineages[lid]: 0 for lid in tids}

opt = LineageOptimization(lxm, lem, lineage_tree,
    first_internal_layer=first_layer, lineage_names=common_lineages,
    lineage_type_code_dict=ccd, exp_norm=2, max_workers=MAX_WORKERS)
print(f"Lineage: xyz={opt.lineage_xyz_cost:.1f}, exp={opt.lineage_exp_cost:.1f}")
print(f"Terms={len(opt.terminal_tree_ids)}, Internal={len(opt.internal_tree_ids)}")

# ═══ NULL MODELS ═══
print("\n" + "=" * 60)
print("NULL MODELS")
print("=" * 60)

null_data = {}

# Random full
path = os.path.join(NULL_DIR, 'random_full.npz')
if os.path.exists(path):
    print("Random full: loading cached")
    d = np.load(path); r_xyz, r_exp = d['xyz'], d['exp']
else:
    print(f"Random full: running ({RANDOM_ITER} iter)")
    res = opt.random_assignment_runner(iterations=RANDOM_ITER)
    r_xyz = np.array([c[0] for c in res]); r_exp = np.array([c[1] for c in res])
    np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
null_data['Random full'] = (r_xyz, r_exp)
print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

# Random by layer
path = os.path.join(NULL_DIR, 'random_by_layer.npz')
if os.path.exists(path):
    print("Random by layer: loading cached")
    d = np.load(path); r_xyz, r_exp = d['xyz'], d['exp']
else:
    print(f"Random by layer: running ({RANDOM_ITER} iter)")
    res = opt.random_assignment_by_layer_runner(iterations=RANDOM_ITER)
    r_xyz = np.array([c[0] for c in res]); r_exp = np.array([c[1] for c in res])
    np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
null_data['Random by layer'] = (r_xyz, r_exp)
print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

# Random cousin
path = os.path.join(NULL_DIR, 'random_cousin.npz')
if os.path.exists(path):
    print("Random cousin: loading cached")
    d = np.load(path); r_xyz, r_exp = d['xyz'], d['exp']
else:
    print(f"Random cousin: running ({RANDOM_ITER} iter)")
    res = opt.random_cousin_shuffle_runner(degree=1, iterations=RANDOM_ITER)
    r_xyz = np.array([c[0] for c in res]); r_exp = np.array([c[1] for c in res])
    np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
null_data['Random cousin'] = (r_xyz, r_exp)
print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

# Random rebuild
path = os.path.join(NULL_DIR, 'random_rebuild.npz')
if os.path.exists(path):
    print("Random rebuild: loading cached")
    d = np.load(path); r_xyz, r_exp = d['xyz'], d['exp']
else:
    print(f"Random rebuild: running ({RANDOM_ITER} iter)")
    res = opt.random_rebuild_runner(iterations=RANDOM_ITER)
    r_xyz = np.array([c[0] for c in res]); r_exp = np.array([c[1] for c in res])
    np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
null_data['Random rebuild'] = (r_xyz, r_exp)
print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

# ═══ HEURISTICS ═══
print("\n" + "=" * 60)
print("HEURISTICS")
print("=" * 60)

heur_data = {}

# Bottom-up by layer
path = os.path.join(HEUR_DIR, 'bottom_up_by_layer.npz')
if os.path.exists(path):
    print("Bottom-up by layer: loading cached")
    d = dict(np.load(path, allow_pickle=True)); h_xyz, h_exp = d['xyz'], d['exp']
else:
    print("Bottom-up by layer: running (1001 steps)")
    res = opt.bottom_up_by_layer_runner()
    h_xyz = np.array([c[0] for c in res]); h_exp = np.array([c[1] for c in res])
    np.savez(path, xyz=h_xyz, exp=h_exp)
heur_data['Bottom-up by layer'] = (h_xyz, h_exp)
print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

# Top-down rebuild
path = os.path.join(HEUR_DIR, 'top_down_rebuild.npz')
if os.path.exists(path):
    print("Top-down rebuild: loading cached")
    d = dict(np.load(path, allow_pickle=True)); h_xyz, h_exp = d['xyz'], d['exp']
else:
    print("Top-down rebuild: running (1001 steps)")
    res = opt.top_down_rebuild_runner()
    h_xyz = np.array([c[0] for c in res]); h_exp = np.array([c[1] for c in res])
    np.savez(path, xyz=h_xyz, exp=h_exp)
heur_data['Top-down rebuild'] = (h_xyz, h_exp)
print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

# Paired bottom-up
path = os.path.join(HEUR_DIR, 'paired_bottom_up.npz')
if os.path.exists(path):
    print("Paired bottom-up: loading cached")
    d = dict(np.load(path, allow_pickle=True)); h_xyz, h_exp = d['xyz'], d['exp']
else:
    print("Paired bottom-up: running (1001 steps, ~1-2 hr)")
    res = opt.paired_bottom_up_rebuild_runner()
    h_xyz = np.array([c[0] for c in res]); h_exp = np.array([c[1] for c in res])
    np.savez(path, xyz=h_xyz, exp=h_exp)
heur_data['Paired bottom-up'] = (h_xyz, h_exp)
print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

# Terminal-only (no kd)
path = os.path.join(HEUR_DIR, 'terminal_only_no_kd.npz')
if os.path.exists(path):
    print("Terminal-only (no kd): loading cached")
    d = dict(np.load(path, allow_pickle=True)); h_xyz, h_exp = d['xyz'], d['exp']
else:
    print("Terminal-only (no kd): running (1001 steps)")
    res = opt.terminal_only_rebuild_runner(use_kd_tree=False, exp_replacement=False)
    h_xyz = np.array([c[0] for c in res]); h_exp = np.array([c[1] for c in res])
    np.savez(path, xyz=h_xyz, exp=h_exp)
heur_data['Terminal-only (no kd)'] = (h_xyz, h_exp)
print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

# Terminal-only (with kd)
path = os.path.join(HEUR_DIR, 'terminal_only_kd.npz')
if os.path.exists(path):
    print("Terminal-only (with kd): loading cached")
    d = dict(np.load(path, allow_pickle=True)); h_xyz, h_exp = d['xyz'], d['exp']
else:
    print("Terminal-only (with kd): running (1001 steps)")
    res = opt.terminal_only_rebuild_runner(use_kd_tree=True, exp_replacement=True)
    h_xyz = np.array([c[0] for c in res]); h_exp = np.array([c[1] for c in res])
    np.savez(path, xyz=h_xyz, exp=h_exp)
heur_data['Terminal-only (with kd)'] = (h_xyz, h_exp)
print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

# ═══ COMPARE WITH OLD RESULTS ═══
print("\n" + "=" * 60)
print("COMPARISON WITH OLD RESULTS")
print("=" * 60)

old_csv = MODULE_DIR / "output/internal_opt/embryo1/10_pca_255_l2/pareto_comparison.csv"
if os.path.exists(old_csv):
    old_df = pd.read_csv(old_csv)
    for idx in [0, 250, 500, 750, 1000]:
        old_xyz = old_df['Bottom-up By Layer XYZ Cost'].iloc[idx]
        old_exp = old_df['Bottom-up By Layer Exp Cost'].iloc[idx]
        our_xyz = heur_data['Bottom-up by layer'][0][idx]
        our_exp = heur_data['Bottom-up by layer'][1][idx]
        xyz_ok = abs(our_xyz - old_xyz) / old_xyz < 0.05
        exp_ok = abs(our_exp - old_exp) / old_exp < 0.05
        print(f"alpha={idx/1000:.3f}: our=({our_xyz:.0f},{our_exp:.0f}) old=({old_xyz:.0f},{old_exp:.0f}) " + ("MATCH" if (xyz_ok and exp_ok) else "DIFF"))

# ═══ SUMMARY ═══
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"Data: 10D PCA, embryo1, T<=255")
print(f"Lineage: xyz={opt.lineage_xyz_cost:.1f}, exp={opt.lineage_exp_cost:.1f}")
print(f"Tree: {lineage_tree.size} nodes, {len(opt.terminal_tree_ids)} terminals")
print(f"\nNull models ({RANDOM_ITER} iterations each):")
for name, (xz, ex) in null_data.items():
    print(f"  {name:20s}: xyz={xz.mean():.0f}±{xz.std():.0f}, exp={ex.mean():.0f}±{ex.std():.0f}")
print(f"\nHeuristic minima:")
for name, (xz, ex) in heur_data.items():
    print(f"  {name:25s}: xyz={xz.min():.0f} ({(xz.min()-opt.lineage_xyz_cost)/opt.lineage_xyz_cost*100:+.1f}%), exp={ex.min():.0f} ({(ex.min()-opt.lineage_exp_cost)/opt.lineage_exp_cost*100:+.1f}%)")

print("\nDone!")
