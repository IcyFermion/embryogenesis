"""Parameterized full-tree Pareto pipeline runner. Accepts config name as argument."""
import json, sys, numpy as np, pandas as pd, os
from pathlib import Path
from collections import defaultdict, deque

MODULE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = MODULE_DIR.parent
CACHE_ROOT = MODULE_DIR / 'output' / 'internal_opt'
sys.path.insert(0, str(PROJECT_ROOT))
from utils import lineage_name_mapping as map_names, load_json
from pareto_core import LineageTree, LineageOptimization

# ═══ Configuration registry ═══
CONFIGS = {
    # ── Protein 10D PCA ──
    'embryo1_prot': {
        'tracking_file': './data/embryo1/tracks.txt',
        'tracking_sep': '\t',
        'tracking_time_col': 't',
        'tracking_cell_col': 'name',
        'tracking_xyz_cols': ['x', 'y', 'z'],
        'tracking_scale': 0.1625,
        'time_cutoff': 255,
        'exp_file': 'data/protein/aggregated_all/s3_zscore.csv',
        'exp_index_col': 0,
        'exp_transpose': True,  # proteins are rows, cells are columns
        'feature_selection_file': 'expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv',
        'feature_selection_col': 'protein',
        'manual_times': {'AB': -3, 'P1': -2, 'P0': -20},
        'output_dir': CACHE_ROOT / 'celegans_prot/embryo1',
    },
    'embryo2_prot': {
        'tracking_file': './data/embryo2/tracks.txt',
        'tracking_sep': '\t',
        'tracking_time_col': 't',
        'tracking_cell_col': 'name',
        'tracking_xyz_cols': ['x', 'y', 'z'],
        'tracking_scale': 0.1625,
        'time_cutoff': 247,
        'exp_file': 'data/protein/aggregated_all/s3_pca_10d.csv',
        'exp_index_col': 0,
        'exp_transpose': False,
        'manual_times': {'AB': -3, 'P1': -2, 'P0': -20},
        'output_dir': CACHE_ROOT / 'celegans_prot/embryo2',
    },
    'embryo3_prot': {
        'tracking_file': './data/embryo3/tracks.txt',
        'tracking_sep': '\t',
        'tracking_time_col': 't',
        'tracking_cell_col': 'name',
        'tracking_xyz_cols': ['x', 'y', 'z'],
        'tracking_scale': 0.1625,
        'time_cutoff': 226,
        'exp_file': 'data/protein/aggregated_all/s3_pca_10d.csv',
        'exp_index_col': 0,
        'exp_transpose': False,
        'manual_times': {'AB': -3, 'P1': -2, 'P0': -20},
        'output_dir': CACHE_ROOT / 'celegans_prot/embryo3',
    },
    # ── C. elegans RNA-seq (2D xy-only for cross-species comparison) ──
    'elegans_rna': {
        'tracking_file': './data/embryo1/tracks.txt',
        'tracking_sep': '\t',
        'tracking_time_col': 't',
        'tracking_cell_col': 'name',
        'tracking_xyz_cols': ['x', 'y', 'z'],
        'tracking_scale': 0.1625,
        'spatial_dims': ['x', 'y'],  # 2D xy-only for comparison with briggsae
        'time_cutoff': 230,  # early cutoff matching pareto_presentation.ipynb
        'exp_file': 'data/c_briggsae/science.adu8249/c_elegans_tf.csv',
        'exp_index_col': 0,
        'exp_transpose': True,  # TFs are rows, cells are columns → transpose
        'feature_selection_file': 'expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv',
        'feature_selection_col': 'tf',
        'manual_times': {'AB': -3, 'P1': -2, 'P0': -20},
        'output_dir': CACHE_ROOT / 'celegans_rna_2d/embryo1',
    },
    # ── C. briggsae RNA-seq (new tracking, 2D xy-only) ──
    'briggsae_rna': {
        'tracking_file': 'data/c_briggsae/yiming/CD210519ZZY0874p1.csv',
        'tracking_sep': ',',
        'tracking_time_col': 'time',
        'tracking_cell_col': 'cell',
        'tracking_xyz_cols': ['x', 'y', 'z'],
        'tracking_scale': 1.0,  # no scaling for C. briggsae
        'spatial_dims': ['x', 'y'],  # 2D xy-only: briggsae z-axis is less accurate
        'time_cutoff': 148,  # comparable developmental stage to elegans T=230
        'exp_file': 'data/c_briggsae/science.adu8249/c_briggsae_tf.csv',
        'exp_index_col': 0,
        'exp_transpose': True,
        'feature_selection_file': 'expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv',
        'feature_selection_col': 'tf',
        'manual_times': {'AB': -3, 'P1': -2, 'P0': -20},
        'output_dir': CACHE_ROOT / 'briggsae_rna_2d/embryo1',
    },
}

RANDOM_ITER = 100_000
MAX_WORKERS = 16


def run_pipeline(cfg_name):
    cfg = CONFIGS[cfg_name]
    OUT = str(cfg['output_dir'])
    NULL_DIR = os.path.join(OUT, 'null')
    HEUR_DIR = os.path.join(OUT, 'heuristics')
    os.makedirs(NULL_DIR, exist_ok=True)
    os.makedirs(HEUR_DIR, exist_ok=True)

    print("=" * 60)
    print(f"CONFIG: {cfg_name}")
    print(f"  Tracking: {cfg['tracking_file']}, T<={cfg['time_cutoff']}")
    print(f"  Expression: {cfg['exp_file']}")
    print(f"  Output: {OUT}")
    print("=" * 60)

    # ═══ DATA LOADING ═══
    lineage_data = load_json(PROJECT_ROOT / 'data/cell_lineage.json')

    # Expression
    lineage_exp_df = pd.read_csv(PROJECT_ROOT / cfg['exp_file'], index_col=cfg['exp_index_col'])
    if cfg['exp_transpose']:
        lineage_exp_df = lineage_exp_df.T
    lineage_exp_df = lineage_exp_df.fillna(0)
    # Feature selection (top20 proteins/TFs matching pareto_presentation.ipynb)
    if 'feature_selection_file' in cfg:
        sel_path = PROJECT_ROOT / cfg['feature_selection_file']
        sel_col = cfg['feature_selection_col']
        # Check if file has a header matching the expected column name
        with open(sel_path) as f:
            first = f.readline().strip()
        if first == sel_col:
            # Has header row
            sel_df = pd.read_csv(sel_path, sep='\t' if sel_path.suffix == '.tsv' else ',')
            sel_features = sel_df[sel_col].tolist()
        else:
            # No header — read as plain list of names
            sel_features = pd.read_csv(sel_path, header=None)[0].tolist()
        available = [f for f in sel_features if f in lineage_exp_df.columns]
        lineage_exp_df = lineage_exp_df[available]
        print(f"Expression: {lineage_exp_df.shape} (filtered to {len(available)}/{len(sel_features)} selected features)")
    else:
        print(f"Expression: {lineage_exp_df.shape}")

    # Tracking
    tracking_df_all = pd.read_csv(PROJECT_ROOT / cfg['tracking_file'], sep=cfg['tracking_sep'])
    time_col = cfg['tracking_time_col']
    cell_col = cfg['tracking_cell_col']
    xyz_cols = cfg['tracking_xyz_cols']
    spatial_dims = cfg.get('spatial_dims', ['x', 'y', 'z'])  # e.g. ['x','y'] for 2D RNA comparison
    tracking_df = tracking_df_all[tracking_df_all[time_col] <= cfg['time_cutoff']]
    S = cfg['tracking_scale']

    # DFS
    terminal_nodes = []
    internal_nodes = []
    parent_dict = {}

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
    cell_names = tracking_df[cell_col].unique()
    lineage_xyz_df = pd.DataFrame(columns=["x", "y", "z"])
    lineage_name_to_last_tracking_time = {}

    for name in cell_names:
        if name not in lineage_exp_df.index:
            continue
        lt = tracking_df[tracking_df[cell_col] == name]
        tp = lt[time_col].values
        if len(tp) == 1 and tp[0] == cfg['time_cutoff']:
            continue
        lineage_name_to_last_tracking_time[name] = int(tp[-1])
        lineage_xyz_df.loc[name] = lt[xyz_cols].values[-1] * S

    # Manual times
    for k, v in cfg['manual_times'].items():
        lineage_name_to_last_tracking_time[k] = v

    # Untracked nodes (hardcoded 3, matching old pipeline)
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
    lineage_tree.add_node(0, -1)
    lineage_tree.root = 0
    lineage_tree.add_node(1, 0)
    lineage_tree.add_node(2, 0)
    bq = deque()
    bq.append((lineage_data, -1))
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
            # Node has no expression/tracking data — skip it but continue
            # exploring descendants (using parent's index for correct linkage)
            for child in node.get('children', []):
                bq.append((child, pidx))
    lineage_tree.record_last_tracked_time(lt_list)
    print(f"Tree size: {lineage_tree.size}")

    # Apply spatial dimension filter (e.g. xy-only for cross-species RNA comparison)
    xyz_df_filtered = lineage_xyz_df[spatial_dims]
    lxm = np.nan_to_num(xyz_df_filtered.loc[common_lineages].values, nan=0.0)
    lem = np.nan_to_num(lineage_exp_df.loc[common_lineages].values.astype(float), nan=0.0)
    print(f'xyz_mat: {lxm.shape} (dims={spatial_dims}), exp_mat: {lem.shape}')

    # Optimization — compute first_layer dynamically from tree topology
    # AB (tree_id 1) and P1 (tree_id 2) are children of P0; their tracked
    # descendants form the first layer of internal nodes for traversal.
    first_layer = []
    for pid in [1, 2]:
        for cid in lineage_tree.children_list[pid]:
            first_layer.append((cid, 2))  # depth 2 from root
    print(f"first_internal_layer: {first_layer}")
    tids = [lineage_tree.lineage_id_mapping[t] for t in range(lineage_tree.size)
            if len(lineage_tree.children_list[t]) == 0]
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
        d = np.load(path)
        r_xyz, r_exp = d['xyz'], d['exp']
    else:
        print(f"Random full: running ({RANDOM_ITER} iter)")
        res = opt.random_assignment_runner(iterations=RANDOM_ITER)
        r_xyz = np.array([c[0] for c in res])
        r_exp = np.array([c[1] for c in res])
        np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
    null_data['Random full'] = (r_xyz, r_exp)
    print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

    # Random by layer
    path = os.path.join(NULL_DIR, 'random_by_layer.npz')
    if os.path.exists(path):
        print("Random by layer: loading cached")
        d = np.load(path)
        r_xyz, r_exp = d['xyz'], d['exp']
    else:
        print(f"Random by layer: running ({RANDOM_ITER} iter)")
        res = opt.random_assignment_by_layer_runner(iterations=RANDOM_ITER)
        r_xyz = np.array([c[0] for c in res])
        r_exp = np.array([c[1] for c in res])
        np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
    null_data['Random by layer'] = (r_xyz, r_exp)
    print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

    # Random cousin
    path = os.path.join(NULL_DIR, 'random_cousin.npz')
    if os.path.exists(path):
        print("Random cousin: loading cached")
        d = np.load(path)
        r_xyz, r_exp = d['xyz'], d['exp']
    else:
        print(f"Random cousin: running ({RANDOM_ITER} iter)")
        res = opt.random_cousin_shuffle_runner(degree=1, iterations=RANDOM_ITER)
        r_xyz = np.array([c[0] for c in res])
        r_exp = np.array([c[1] for c in res])
        np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
    null_data['Random cousin'] = (r_xyz, r_exp)
    print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

    # Random rebuild
    path = os.path.join(NULL_DIR, 'random_rebuild.npz')
    if os.path.exists(path):
        print("Random rebuild: loading cached")
        d = np.load(path)
        r_xyz, r_exp = d['xyz'], d['exp']
    else:
        print(f"Random rebuild: running ({RANDOM_ITER} iter)")
        res = opt.random_rebuild_runner(iterations=RANDOM_ITER)
        r_xyz = np.array([c[0] for c in res])
        r_exp = np.array([c[1] for c in res])
        np.savez_compressed(path, xyz=r_xyz, exp=r_exp)
    null_data['Random rebuild'] = (r_xyz, r_exp)
    print(f"  xyz={r_xyz.mean():.0f}±{r_xyz.std():.0f}, exp={r_exp.mean():.0f}±{r_exp.std():.0f}")

    # ═══ HEURISTICS ═══
    print("\n" + "=" * 60)
    print("HEURISTICS")
    print("=" * 60)

    skip = cfg.get('skip_heuristics', [])

    heur_data = {}

    # Bottom-up by layer
    path = os.path.join(HEUR_DIR, 'bottom_up_by_layer.npz')
    if os.path.exists(path):
        print("Bottom-up by layer: loading cached")
        d = dict(np.load(path, allow_pickle=True))
        h_xyz, h_exp = d['xyz'], d['exp']
    else:
        print("Bottom-up by layer: running (1001 steps)")
        res = opt.bottom_up_by_layer_runner()
        h_xyz = np.array([c[0] for c in res])
        h_exp = np.array([c[1] for c in res])
        np.savez(path, xyz=h_xyz, exp=h_exp)
    heur_data['Bottom-up by layer'] = (h_xyz, h_exp)
    print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

    # Layerwise Pareto assignment (per-layer Hungarian edge reassignment)
    path = os.path.join(HEUR_DIR, 'layerwise_pareto_assignment.npz')
    if os.path.exists(path):
        print("Layerwise Pareto assignment: loading cached")
        d = dict(np.load(path, allow_pickle=True))
        h_xyz, h_exp = d['xyz'], d['exp']
    else:
        print("Layerwise Pareto assignment: running (1001 steps)")
        res = opt.layerwise_pareto_assignment_runner()
        h_xyz = res['xyz']
        h_exp = res['exp']
        np.savez_compressed(path,
            xyz=h_xyz, exp=h_exp,
            per_layer_pareto_fronts=np.array(res['per_layer_pareto_fronts'], dtype=object),
            per_layer_edge_retention=np.array(res['per_layer_edge_retention'], dtype=object),
            per_layer_costs=np.array(res['per_layer_costs'], dtype=object),
            layer_sizes=np.array(res['layer_sizes']))
    heur_data['Layerwise Pareto assignment'] = (h_xyz, h_exp)
    print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

    # Top-down rebuild
    if 'top_down_rebuild' not in skip:
        path = os.path.join(HEUR_DIR, 'top_down_rebuild.npz')
        if os.path.exists(path):
            print("Top-down rebuild: loading cached")
            d = dict(np.load(path, allow_pickle=True))
            h_xyz, h_exp = d['xyz'], d['exp']
        else:
            print("Top-down rebuild: running (1001 steps)")
            res = opt.top_down_rebuild_runner()
            h_xyz = np.array([c[0] for c in res])
            h_exp = np.array([c[1] for c in res])
            np.savez(path, xyz=h_xyz, exp=h_exp)
        heur_data['Top-down rebuild'] = (h_xyz, h_exp)
        print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")
    else:
        print("Top-down rebuild: SKIPPED (known hang with skip-through tree)")

    # Paired bottom-up
    path = os.path.join(HEUR_DIR, 'paired_bottom_up.npz')
    if os.path.exists(path):
        print("Paired bottom-up: loading cached")
        d = dict(np.load(path, allow_pickle=True))
        h_xyz, h_exp = d['xyz'], d['exp']
    else:
        print("Paired bottom-up: running (1001 steps, ~1-2 hr)")
        res = opt.paired_bottom_up_rebuild_runner()
        h_xyz = np.array([c[0] for c in res])
        h_exp = np.array([c[1] for c in res])
        np.savez(path, xyz=h_xyz, exp=h_exp)
    heur_data['Paired bottom-up'] = (h_xyz, h_exp)
    print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

    # MST-based rebuild
    path = os.path.join(HEUR_DIR, 'mst_rebuild.npz')
    if os.path.exists(path):
        print("MST rebuild: loading cached")
        d = dict(np.load(path, allow_pickle=True))
        h_xyz, h_exp = d['xyz'], d['exp']
    else:
        print("MST rebuild: running (1001 steps)")
        res = opt.mst_rebuild_runner()
        h_xyz = np.array([c[0] for c in res])
        h_exp = np.array([c[1] for c in res])
        np.savez(path, xyz=h_xyz, exp=h_exp)
    heur_data['MST rebuild'] = (h_xyz, h_exp)
    print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

    # Terminal-only (no exp replacement)
    path = os.path.join(HEUR_DIR, 'terminal_only_no_replacement.npz')
    if os.path.exists(path):
        print("Terminal-only (no replacement): loading cached")
        d = dict(np.load(path, allow_pickle=True))
        h_xyz, h_exp = d['xyz'], d['exp']
    else:
        print("Terminal-only (no replacement): running (1001 steps)")
        res = opt.terminal_only_rebuild_runner(use_kd_tree=True, exp_replacement=False)
        h_xyz = np.array([c[0] for c in res])
        h_exp = np.array([c[1] for c in res])
        np.savez(path, xyz=h_xyz, exp=h_exp)
    heur_data['Terminal-only (no repl.)'] = (h_xyz, h_exp)
    print(f"  xyz=[{h_xyz.min():.0f}, {h_xyz.max():.0f}], exp=[{h_exp.min():.0f}, {h_exp.max():.0f}]")

    # ═══ SUMMARY ═══
    print("\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"Config: {cfg_name}")
    print(f"Tree: {lineage_tree.size} nodes, {len(opt.terminal_tree_ids)} terminals")
    print(f"Lineage: xyz={opt.lineage_xyz_cost:.1f}, exp={opt.lineage_exp_cost:.1f}")
    print(f"\nNull models ({RANDOM_ITER} iterations each):")
    for name, (xz, ex) in null_data.items():
        print(f"  {name:20s}: xyz={xz.mean():.0f}±{xz.std():.0f}, exp={ex.mean():.0f}±{ex.std():.0f}")
    print(f"\nHeuristic minima:")
    for name, (xz, ex) in heur_data.items():
        print(f"  {name:25s}: xyz={xz.min():.0f} ({(xz.min()-opt.lineage_xyz_cost)/opt.lineage_xyz_cost*100:+.1f}%), exp={ex.min():.0f} ({(ex.min()-opt.lineage_exp_cost)/opt.lineage_exp_cost*100:+.1f}%)")

    print("\nDone!")
    return opt, null_data, heur_data


if __name__ == '__main__':
    cfg = sys.argv[1] if len(sys.argv) > 1 else 'embryo1_prot'
    run_pipeline(cfg)
