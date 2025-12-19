import pandas as pd
import numpy as np
from collections import deque, defaultdict
from utils import lineage_name_mapping, load_json, bidict
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import seaborn as sns
import sys
from pareto_core import LineageOptimization, LineageTree

if __name__ == "__main__":
    num_threads = int(sys.argv[1])  # Read number of threads from SLURM export
    lineage_data = load_json('./data/cell_lineage.json')
    lineage_exp_df = pd.read_csv("data/protein/aggregated_all/s3_pca_10d.csv", index_col=0)
    lineage_exp_df = lineage_exp_df.fillna(0)
    tracking_df_all = pd.read_csv("./data/embryo3/tracks.txt", sep="\t")
    tracking_time_cutoff = 226
    tracking_scale = 0.1625
    tracking_df = tracking_df_all.loc[tracking_df_all["t"] <= tracking_time_cutoff]
    untracked_nodes = ["P0", "AB", "P1"]
    first_internal_layer = ["ABa", "ABp", "EMS", "P2"]

    def depth_cutoff_func(depth):
        if depth <= 9:
            return 1
        else:
            return float('inf')
        
    terminal_nodes = []
    internal_nodes = []
    parent_dict = {}
    terminal_parents_dict = defaultdict(list)
    terminal_ancestry = defaultdict(list)
    def dfs(node, parent, ancestors=[]):
        children = node.get("children", [])
        lookup_name = lineage_name_mapping(node["did"])
        if len(children) == 0:
            p_lookup_name = lineage_name_mapping(parent['did'])
            parent_dict[lookup_name] = p_lookup_name
            terminal_nodes.append(lookup_name)
            terminal_parents_dict[p_lookup_name].append(lookup_name)
            terminal_ancestry[lookup_name] = ancestors
        else:
            internal_nodes.append(lookup_name)
            if parent is not None:
                p_lookup_name = lineage_name_mapping(parent['did'])
                parent_dict[lookup_name] = p_lookup_name
            for child in children:
                dfs(child, node, ancestors + [lookup_name])

    dfs(lineage_data, None)

    cell_type_df = pd.read_csv("data/2023-06-29_entropy_cell_key_V2.csv")
    count = 0
    lineage_type_code_dict = defaultdict(str)
    type_code_dict = {}
    na_lineages = []
    for lineage in terminal_nodes:
        cur_cell_type_df = cell_type_df[cell_type_df['wormweb.lineage'] == lineage]
        if cur_cell_type_df.empty:
            # print(f"No cell types found for lineage: {lineage}")
            na_lineages.append(lineage)
            continue
        cur_lineage_types = cur_cell_type_df["wormweb.type"]
        # remove nan from cur_lineage_types
        cur_lineage_types = cur_lineage_types[~cur_lineage_types.isna()]
        cur_lineage_types = cur_lineage_types.unique()
        if len(cur_lineage_types) == 0:
            # print(f"No cell types found for lineage: {lineage}")
            na_lineages.append(lineage)
            continue
        cur_type = cur_lineage_types[0]
        if cur_type not in type_code_dict:
            # type_code_dict[cur_type] = chr(97 + len(type_code_dict))  # start from 'a'
            type_code_dict[cur_type] = len(type_code_dict)
        lineage_type_code_dict[lineage] = type_code_dict[cur_type]
    # dead/NA type
    na_type_code = len(type_code_dict)
    type_code_dict["dead/NA"] = na_type_code # unknown/dead type
    for lineage in na_lineages:
        lineage_type_code_dict[lineage] = na_type_code

    complexity_code_dict = {}
    cell_div_count = 0
    parent_code_mapping = {}
    for node in terminal_nodes:
        complexity_code_dict[node] = lineage_type_code_dict[node]
    bfs_queue = deque(terminal_nodes)
    while bfs_queue:
        cur_node = bfs_queue.popleft()
        cur_parent = parent_dict.get(cur_node, None)
        cur_node_code = complexity_code_dict.get(cur_node, None)
        if cur_node_code is None:
            print(f"Type code not found for node: {cur_node}")
            continue
        if cur_parent is None:
            print(f"Parent not found for node: {cur_node}")
            continue
        if cur_parent not in complexity_code_dict:
            complexity_code_dict[cur_parent] = cur_node_code
            bfs_queue.append(cur_parent)
        else:
            cur_parent_code = complexity_code_dict[cur_parent]
            if cur_parent_code > cur_node_code:
                parent_code_tuple = (cur_parent_code, cur_node_code)
            else:
                parent_code_tuple = (cur_node_code, cur_parent_code)
            if parent_code_tuple not in parent_code_mapping:
                parent_code_mapping[parent_code_tuple] = len(parent_code_mapping)+len(type_code_dict)
            complexity_code_dict[cur_parent] = parent_code_mapping[parent_code_tuple]
            cell_div_count += 1

    cell_names = tracking_df["name"].unique()
    lineage_xyz_df = pd.DataFrame(columns=["x", "y", "z"])
    lineage_name_to_last_tracking_time = {}
    for name in cell_names:
        time_points = tracking_df.loc[tracking_df['name'] == name]["t"].values
        if len(time_points) == 1 and time_points[0] == tracking_time_cutoff:
            continue
        lineage_tracks = tracking_df.loc[tracking_df['name'] == name]
        last_xyz_coordinates = lineage_tracks[["x", "y", "z"]].values[-1]
        last_tracking_time = lineage_tracks["t"].values[-1]
        lineage_name_to_last_tracking_time[name] = int(last_tracking_time)
        # appending last_xyz_coordinates to lineage_xyz_df as a new row with cell name as index
        lineage_xyz_df.loc[name] = last_xyz_coordinates * tracking_scale

    # manual estimates based on video and wormweb data for untracked nodes
    lineage_name_to_last_tracking_time["AB"] = -11
    lineage_name_to_last_tracking_time["P1"] = -10
    lineage_name_to_last_tracking_time["P0"] = -28

    common_lineages = lineage_xyz_df.index.intersection(lineage_exp_df.index)
    for name in untracked_nodes:
        lineage_xyz_df.loc[name] = np.repeat(np.nan, 3)
        lineage_exp_df.loc[name] = np.repeat(np.nan, lineage_exp_df.shape[1])
    common_lineages = untracked_nodes + common_lineages.tolist()
    lineage_id_to_last_tracking_time = [-1] * len(common_lineages)
    for idx, name in enumerate(common_lineages):
        lineage_id_to_last_tracking_time[idx] = lineage_name_to_last_tracking_time.get(name, -1)
        if lineage_name_to_last_tracking_time.get(name, -1) == -1:
            print(f"Warning: No last tracking time for lineage {name}")

    # pruning common lineages to include those only leads to terminal nodes
    pruned_common_lineages = set()
    bfs_queue = deque(terminal_nodes)
    while bfs_queue:
        cur_node = bfs_queue.popleft()
        if cur_node in common_lineages:
            pruned_common_lineages.add(cur_node)
            if cur_node in parent_dict:
                bfs_queue.append(parent_dict[cur_node])

    # using pruned common lineages, comment out to use all common lineages
    new_common_lineages = []
    for lineage in common_lineages:
        if lineage in pruned_common_lineages:
            new_common_lineages.append(lineage)
    old_commmon_lineages = common_lineages.copy()
    # common_lineages = new_common_lineages

    lineage_tree = LineageTree()
    # add first three nodes: ["AB", "P0", "P1"]
    missing_first_nodes = ["P0", "AB", "P1"]
    lineage_tree.add_node(0, -1)
    lineage_tree.root = 0
    lineage_tree.add_node(1, 0)
    lineage_tree.add_node(2, 0)
    bfs_queue = deque()
    bfs_queue.append((lineage_data, -1))
    while bfs_queue:
        node, parent_idx = bfs_queue.popleft()
        curr_name = lineage_name_mapping(node['did'])
        if curr_name not in common_lineages:
            continue
        curr_idx = common_lineages.index(curr_name)
        if curr_idx >= 3:
            parent_idx = lineage_tree.reverse_lineage_id_mapping[parent_idx]
            lineage_tree.add_node(curr_idx, parent_idx)
        children = node.get('children', [])
        for child_node in children:
            bfs_queue.append((child_node, curr_idx))
    lineage_tree.record_last_tracked_time(lineage_id_to_last_tracking_time)

    print("first layer of tracked nodes:", lineage_tree.children_list[1]+lineage_tree.children_list[2])

    lineage_xyz_mat = lineage_xyz_df.loc[common_lineages].values
    lineage_exp_mat = lineage_exp_df.loc[common_lineages].values

    first_layer = [(3,2), (4,2), (5,2), (6,2)]
    # first_layer = [(5,2)]

    opt = LineageOptimization(
        lineage_xyz_mat, 
        lineage_exp_mat,
        lineage_tree,
        first_internal_layer=first_layer,
        lineage_names=common_lineages, 
        lineage_type_code_dict=complexity_code_dict,
        exp_norm=2,
        max_workers=num_threads)
    
    # pareto_list_by_layer = opt.bottom_up_by_layer_runner()
    # pareto_list_top_down_rebuild = opt.top_down_rebuild_runner()
    # pareto_mst_list = opt.mst_rebuild_runner()
    # pareto_list_top_down_rebuild_depth_cutoff = opt.top_down_rebuild_runner(depth_weight_type=depth_cutoff_func)

    # pareto_xyz_cost_list_by_layer = [cost[0] for cost in pareto_list_by_layer]
    # pareto_exp_cost_list_by_layer = [cost[1] for cost in pareto_list_by_layer]
    # pareto_xyz_cost_list_top_down_rebuild = [cost[0] for cost in pareto_list_top_down_rebuild]
    # pareto_exp_cost_list_top_down_rebuild = [cost[1] for cost in pareto_list_top_down_rebuild]
    # pareto_xyz_cost_list_mst = [cost[0] for cost in pareto_mst_list]
    # pareto_exp_cost_list_mst = [cost[1] for cost in pareto_mst_list]
    # pareto_xyz_cost_list_top_down_rebuild_depth_cutoff = [cost[0] for cost in pareto_list_top_down_rebuild_depth_cutoff]
    # pareto_exp_cost_list_top_down_rebuild_depth_cutoff = [cost[1] for cost in pareto_list_top_down_rebuild_depth_cutoff]
    
    pareto_list_terminal_only_kd_internal = opt.terminal_only_rebuild_runner(use_kd_tree=True, exp_replacement=False)
    pareto_xyz_cost_list_terminal_only_kd_internal = [cost[0] for cost in pareto_list_terminal_only_kd_internal]
    pareto_exp_cost_list_terminal_only_kd_internal = [cost[1] for cost in pareto_list_terminal_only_kd_internal]

    # save pareto_list_terminal_only_kd_internal as a pandas dataframe
    pareto_df_terminal_only_kd_internal = pd.DataFrame(pareto_list_terminal_only_kd_internal, columns=['xyz_cost', 'exp_cost', 'internal_selection_count', 'complexity_score', 'min_consines', 'mean_consines'])
    pareto_df_terminal_only_kd_internal.to_csv('output/internal_opt/embryo3/10_pca_226_l2/pareto_terminal_only_kd_internal_no_exp_replacement.tsv', index=False, sep='\t')

    pareto_list_terminal_only_kd_internal = opt.terminal_only_rebuild_runner(use_kd_tree=True)
    pareto_xyz_cost_list_terminal_only_kd_internal = [cost[0] for cost in pareto_list_terminal_only_kd_internal]
    pareto_exp_cost_list_terminal_only_kd_internal = [cost[1] for cost in pareto_list_terminal_only_kd_internal]

    # save pareto_list_terminal_only_kd_internal as a pandas dataframe
    pareto_df_terminal_only_kd_internal = pd.DataFrame(pareto_list_terminal_only_kd_internal, columns=['xyz_cost', 'exp_cost', 'internal_selection_count', 'complexity_score', 'min_consines', 'mean_consines'])
    pareto_df_terminal_only_kd_internal.to_csv('output/internal_opt/embryo3/10_pca_226_l2/pareto_terminal_only_kd_internal.tsv', index=False, sep='\t')

    pareto_list_paired_bottom_up_rebuild = opt.paired_bottom_up_rebuild_runner()
    pareto_xyz_cost_list_paired_bottom_up_rebuild = [cost[0] for cost in pareto_list_paired_bottom_up_rebuild]
    pareto_exp_cost_list_paired_bottom_up_rebuild = [cost[1] for cost in pareto_list_paired_bottom_up_rebuild]

    # results_df = pd.DataFrame({
    #     'Bottom-up By Layer XYZ Cost': pareto_xyz_cost_list_by_layer,
    #     'Bottom-up By Layer Exp Cost': pareto_exp_cost_list_by_layer,
    #     'Top-Down Rebuild XYZ Cost': pareto_xyz_cost_list_top_down_rebuild,
    #     'Top-Down Rebuild Exp Cost': pareto_exp_cost_list_top_down_rebuild,
    #     'Paired Bottom-Up Rebuild XYZ Cost': pareto_xyz_cost_list_paired_bottom_up_rebuild,
    #     'Paired Bottom-Up Rebuild Exp Cost': pareto_exp_cost_list_paired_bottom_up_rebuild,
    #     'MST Based Rebuild XYZ Cost': pareto_xyz_cost_list_mst,
    #     'MST Based Rebuild Exp Cost': pareto_exp_cost_list_mst,
    #     'Top-Down Rebuild (Hard Depth Cutoff) XYZ Cost': pareto_xyz_cost_list_top_down_rebuild_depth_cutoff,
    #     'Top-Down Rebuild (Hard Depth Cutoff) Exp Cost': pareto_exp_cost_list_top_down_rebuild_depth_cutoff
    # })
    # results_df.to_csv('output/internal_opt/embryo3/10_pca_226_l2/pareto_comparison.csv', index=False)

    


        
        
    