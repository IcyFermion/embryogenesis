import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from collections import defaultdict, deque
from multiprocessing import Manager
from scipy.optimize import linear_sum_assignment, minimize
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm
import networkx as nx
import random
import heapq
import json
from copy import deepcopy
from utils import lineage_name_mapping, load_json, bidict
from zss import simple_distance, Node
from functools import partial
from itertools import combinations
from typing import Callable
from scipy.spatial import KDTree
import kdtree


class LineageTree:
    """"
    A binary tree class representing a lineage tree.
    Each node will have a unique internal id starting from 0
    When a node is created and added to the tree, it will be assigned the size of the tree as its id.
    Children information is stored in children_list, where children_list[i] is a list of the internal ids of the children of node i.
    Parent information is stored in parent_list, where parent_list[i] is the internal id of the parent of node i. The root node's parent is -1.
    The size of the tree is stored in size.
    Each node also carry an interger value representing the lineage identity.
    The mapping from internal id to interger value is stored in list lineage_id_list, where lineage_id_list[i] is the interger value of node i.
    """
    def __init__(self, size=0):
        self.size = size
        self.children_list = [[]] * size
        self.parent_list = [-1] * size
        self.lineage_id_mapping = []
        self.reverse_lineage_id_mapping = {}
        self.root = -1
    
    def add_node(self, lineage_id, parent_id=-1):
        current_id = self.size
        # lineage_id_mapping: tree internal id -> lineage id
        self.lineage_id_mapping.append(lineage_id)
        # reverse_lineage_id_mapping: lineage id -> tree internal id
        self.reverse_lineage_id_mapping[lineage_id] = current_id
        self.children_list.append([])
        self.parent_list.append(parent_id)
        if parent_id != -1:
            self.children_list[parent_id].append(current_id)
        self.size += 1
        return current_id
    
    def record_last_tracked_time(self, lineage_id_to_last_tracking_time: list[int]):
        # last_tracked_time_list: tree internal id -> last tracked time
        self.last_tracked_time_list = [-1] * self.size
        for tree_id in range(self.size):
            lineage_id = self.lineage_id_mapping[tree_id]
            if lineage_id != -1:
                self.last_tracked_time_list[tree_id] = lineage_id_to_last_tracking_time[lineage_id]
        # branch_time_length: tree internal id -> time length from parent to this node
        self.branch_time_length = [0] * self.size
        for tree_id in range(self.size):
            parent_id = self.parent_list[tree_id]
            if parent_id != -1:
                self.branch_time_length[tree_id] = self.last_tracked_time_list[tree_id] - self.last_tracked_time_list[parent_id]

    def path_from_root(self, tree_id: int):
        path = []
        current_id = tree_id
        while current_id != -1:
            path.append(current_id)
            current_id = self.parent_list[current_id]
        path.reverse()
        return path
    
    def shared_time_from_root(self, tree_id1: int, tree_id2: int):
        path1 = self.path_from_root(tree_id1)
        path2 = self.path_from_root(tree_id2)
        shared_time = 0
        for id1, id2 in zip(path1, path2):
            if id1 == id2:
                shared_time += self.branch_time_length[id1]
            else:
                break
        return shared_time
    
    def compute_vcv_matrix(self, terminal_tree_ids: list[int], sigma_sq=1.0):
        vcv_matrix = np.zeros((len(terminal_tree_ids), len(terminal_tree_ids)))
        for i in range(len(terminal_tree_ids)):
            for j in range(len(terminal_tree_ids)):
                tree_id_i = terminal_tree_ids[i]
                tree_id_j = terminal_tree_ids[j]
                vcv_matrix[i][j] = sigma_sq * self.shared_time_from_root(tree_id_i, tree_id_j)
                vcv_matrix[j][i] = vcv_matrix[i][j]
        return vcv_matrix
    

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        new_tree = LineageTree(self.size)
        memo[id(self)] = new_tree
        new_tree.children_list = [children.copy() for children in self.children_list]
        new_tree.parent_list = self.parent_list.copy()
        new_tree.lineage_id_mapping = self.lineage_id_mapping.copy()
        new_tree.reverse_lineage_id_mapping = self.reverse_lineage_id_mapping.copy()
        new_tree.root = self.root
        return new_tree


class LineageOptimization:
    def __init__(self, 
                 xyz_mat, 
                 exp_mat, 
                 lineage_tree: LineageTree, 
                 first_internal_layer: list[tuple[int, int]],
                 lineage_names: list[str],
                 lineage_type_code_dict: dict,
                 exp_norm = 2,
                 max_workers=10):
        self.xyz_mat = xyz_mat
        self.exp_mat = exp_mat
        self.lineage_tree = lineage_tree
        self.lineage_names = lineage_names
        self.first_internal_layer = first_internal_layer
        self.exp_norm = exp_norm
        self.max_workers = max_workers
        # gather untracked lineage ids
        untracked_tree_ids = set()
        for tree_id, _depth in first_internal_layer:
            cur_id = tree_id
            while cur_id != -1:
                parent_tree_id = self.lineage_tree.parent_list[cur_id]
                if parent_tree_id != -1:
                    untracked_tree_ids.add(parent_tree_id)
                cur_id = parent_tree_id
        self.untracked_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in untracked_tree_ids]

        self.tree_ids_by_depth, self.terminal_tree_ids, self.internal_tree_ids = self.lineage_traverse()
        self.lineage_xyz_cost, self.lineage_exp_cost = self.calc_lineage_cost(debugging=True)
        terminal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        # complexity_code_dict: terminal_lineage_id -> new type code starting from 0
        complexity_code_dict = {}
        new_type_converter = {}
        for lineage_id in terminal_lineage_ids:
            original_type_code = lineage_type_code_dict[self.lineage_names[lineage_id]]
            # some monkey fix to force the type of cutoff terminals to be the same
            original_type_code = min(original_type_code, 18)
            if original_type_code not in new_type_converter:
                new_type_converter[original_type_code] = len(new_type_converter)
            complexity_code_dict[lineage_id] = new_type_converter[original_type_code]
        self.terminal_type_code_dict = complexity_code_dict.copy()
        self.terminal_type_nums = len(new_type_converter)
        # compute complexity score for this lineage tree

        lineage_parent_list = [-1] * self.lineage_tree.size
        for tree_id in range(self.lineage_tree.size):
            parent_tree_id = self.lineage_tree.parent_list[tree_id]
            if parent_tree_id != -1:
                lineage_parent_list[self.lineage_tree.lineage_id_mapping[tree_id]] = self.lineage_tree.lineage_id_mapping[parent_tree_id]
        complexity_score = self.complexity_score(lineage_parent_list)
        print("Complexity Score of the original lineage tree:", complexity_score)


        self.xyz_cost_mat = np.zeros((self.lineage_tree.size, self.lineage_tree.size))
        self.exp_cost_mat = np.zeros((self.lineage_tree.size, self.lineage_tree.size))
        for i in range(self.lineage_tree.size):
            for j in range(self.lineage_tree.size):
                self.xyz_cost_mat[i, j] = np.linalg.norm(self.xyz_mat[i] - self.xyz_mat[j])
                self.exp_cost_mat[i, j] = np.linalg.norm(self.exp_mat[i] - self.exp_mat[j], ord=exp_norm)
        
        # internal_exp_map = self.exp_mat[[self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in self.internal_tree_ids]]
        # self.internal_exp_kd_tree = kdtree.create(internal_exp_map.tolist())

    def genetic_complexity_optimization(self, iterations=100000):
        optimization_cost_list = []
        lineage_ids = list(range(self.lineage_tree.size))
        lineage_parent_list = [-1] * self.lineage_tree.size
        for tree_id in range(self.lineage_tree.size):
            parent_tree_id = self.lineage_tree.parent_list[tree_id]
            cur_lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
            if parent_tree_id != -1:
                lineage_parent_list[cur_lineage_id] = self.lineage_tree.lineage_id_mapping[parent_tree_id]
        complexity_score = self.complexity_score(lineage_parent_list)
        print("Initial Complexity Score:", complexity_score)
        for it in tqdm(range(iterations)):
            # randomly swap two nodes
            lineage_id1, lineage_id2 = random.sample(lineage_ids, 2)
            if lineage_id1 in self.untracked_lineage_ids or lineage_id2 in self.untracked_lineage_ids:
                continue
            cur_lineage_parent_list = deepcopy(lineage_parent_list)
            parent_lineage_id1 = lineage_parent_list[lineage_id1]
            parent_lineage_id2 = lineage_parent_list[lineage_id2]
            cur_lineage_parent_list[lineage_id1] = parent_lineage_id2
            cur_lineage_parent_list[lineage_id2] = parent_lineage_id1
            new_complexity_score = self.complexity_score(cur_lineage_parent_list)
            if new_complexity_score < complexity_score:
                complexity_score = new_complexity_score
                lineage_parent_list = cur_lineage_parent_list
                # print(f"Iteration {it}: New best complexity score: {complexity_score}")
                optimization_cost_list.append(self.cost_with_lineage_parent_id_list(lineage_parent_list))

        return lineage_parent_list, optimization_cost_list, complexity_score

    def cost_with_lineage_parent_id_list(self, lineage_parent_list: list[int]):
        xyz_cost = 0
        exp_cost = 0
        for lineage_id, parent_lineage_id in enumerate(lineage_parent_list):
            if parent_lineage_id not in self.untracked_lineage_ids and parent_lineage_id != -1:
                xyz_cost += self.xyz_cost_mat[lineage_id][parent_lineage_id]
                exp_cost += self.exp_cost_mat[lineage_id][parent_lineage_id]
        return xyz_cost, exp_cost

    def complexity_score(self, lineage_parent_list: list[int]):
        terminal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        complexity_code_dict = self.terminal_type_code_dict.copy()
        cell_div_count = 0
        parent_code_mapping = {}
        bfs_queue = deque(terminal_lineage_ids)
        while bfs_queue:
            cur_node_lineage_id = bfs_queue.popleft()
            cur_parent_lineage_id = lineage_parent_list[cur_node_lineage_id] if cur_node_lineage_id != -1 else None
            cur_node_code = complexity_code_dict[cur_node_lineage_id]
            if cur_node_code is None:
                # print(f"Type code not found for node: {cur_node_lineage_id}")
                continue
            if cur_parent_lineage_id is None:
                # print(f"Parent not found for node: {cur_node_lineage_id}")
                continue
            if cur_parent_lineage_id not in complexity_code_dict:
                complexity_code_dict[cur_parent_lineage_id] = cur_node_code
                bfs_queue.append(cur_parent_lineage_id)
            else:
                cur_parent_code = complexity_code_dict[cur_parent_lineage_id]
                if cur_parent_code > cur_node_code:
                    parent_code_tuple = (cur_parent_code, cur_node_code)
                else:
                    parent_code_tuple = (cur_node_code, cur_parent_code)
                if parent_code_tuple not in parent_code_mapping:
                    parent_code_mapping[parent_code_tuple] = len(parent_code_mapping)+self.terminal_type_nums
                complexity_code_dict[cur_parent_lineage_id] = parent_code_mapping[parent_code_tuple]
                cell_div_count += 1

        return len(parent_code_mapping) / cell_div_count
    
    def random_complexity_score_runner(self, iterations=1000000):
        terminal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        random_complexity_scores = process_map(
            partial(self.random_complexity_score, terminal_lineage_ids), 
            range(iterations),
            max_workers=self.max_workers,
            chunksize=200,
            desc="Computing random complexity scores..."
        )
        return random_complexity_scores


    def random_complexity_score(self, terminal_lineage_ids: list[int], seed: int):
        random.seed(seed)
        type_code_list = [self.terminal_type_code_dict[lineage_id] for lineage_id in terminal_lineage_ids]
        random.shuffle(type_code_list)
        parent_code_mapping = {}
        bottom_layer = list(np.arange(len(type_code_list)))
        cell_div_count = 0
        while len(bottom_layer) > 1:
            next_bottom_layer = []
            if len(bottom_layer) % 2 != 0:
                next_bottom_layer.append(bottom_layer[-1])
            for i in range(0, len(bottom_layer), 2):
                if i + 1 >= len(bottom_layer):
                    break
                code1 = type_code_list[bottom_layer[i]]
                code2 = type_code_list[bottom_layer[i+1]]
                if code1 > code2:
                    parent_code_tuple = (code1, code2)
                else:
                    parent_code_tuple = (code2, code1)
                if parent_code_tuple not in parent_code_mapping:
                    parent_code_mapping[parent_code_tuple] = len(parent_code_mapping)+self.terminal_type_nums
                parent_code = parent_code_mapping[parent_code_tuple]
                next_bottom_layer.append(len(type_code_list))
                type_code_list.append(parent_code)
                cell_div_count += 1
            bottom_layer = next_bottom_layer
        return len(parent_code_mapping) / cell_div_count


    def lineage_traverse(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            first_internal_layer = self.first_internal_layer
        tree_ids_by_depth = defaultdict(list)
        terminal_tree_ids = []
        internal_tree_ids = []
        bfs_queue = deque(first_internal_layer)
            # tree_ids_by_depth[depth].append(tree_id)
        while bfs_queue:
            tree_id, depth = bfs_queue.popleft()
            tree_ids_by_depth[depth].append(tree_id)
            if len(self.lineage_tree.children_list[tree_id]) == 0:
                terminal_tree_ids.append(tree_id)
            else:
                internal_tree_ids.append(tree_id)
            for child in self.lineage_tree.children_list[tree_id]:
                bfs_queue.append((child, depth + 1))
        return tree_ids_by_depth, terminal_tree_ids, internal_tree_ids

    def calc_lineage_cost(self, first_internal_tree_ids: list[int]=None, lineage_id_mapping=None, debugging=False):
        if lineage_id_mapping is None:
            lineage_id_mapping = self.lineage_tree.lineage_id_mapping
        if first_internal_tree_ids is None:
            first_internal_tree_ids = [tree_id for tree_id, depth in self.first_internal_layer]
        xyz_cost = 0
        exp_cost = 0
        edge_count = 0
        bfs_queue = deque(first_internal_tree_ids)
        while bfs_queue:
            tree_id = bfs_queue.popleft()
            lineage_id = lineage_id_mapping[tree_id]
            for child in self.lineage_tree.children_list[tree_id]:
                child_lineage_id = lineage_id_mapping[child]
                xyz_cost += np.linalg.norm(self.xyz_mat[lineage_id] - self.xyz_mat[child_lineage_id])
                exp_cost += np.linalg.norm(self.exp_mat[lineage_id] - self.exp_mat[child_lineage_id], ord=self.exp_norm)
                bfs_queue.append(child)
                edge_count += 1
        if debugging:
            print(xyz_cost, exp_cost, edge_count)
        return xyz_cost, exp_cost
    
    def group_cousins_by_degree(self, degree=1):
        cousin_groups = defaultdict(list)
        for tree_id in self.terminal_tree_ids:
            path = self.lineage_tree.path_from_root(tree_id)
            if len(path) > degree + 1:
                cousin_ancestor_id = path[-(degree + 2)]
                cousin_groups[cousin_ancestor_id].append(tree_id)
        cousin_groups_list = list(cousin_groups.values())
        return cousin_groups_list
    
    def random_cousin_shuffle_runner(self, degree=1, iterations=1000000):
        cousin_groups = self.group_cousins_by_degree(degree)
        random_cousin_shuffle_costs = process_map(
            partial(self.random_cousin_shuffle, cousin_groups),
            range(iterations),
            max_workers=self.max_workers,
            chunksize=200,
            desc="Computing random cousin shuffle costs..."
        )
        return random_cousin_shuffle_costs

    def random_cousin_shuffle(self, cousin_groups, seed):
        random.seed(seed)
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        for group in cousin_groups:
            lineage_ids = [cur_lineage_id_mapping[tree_id] for tree_id in group]
            random.shuffle(lineage_ids)
            for tree_id, lineage_id in zip(group, lineage_ids):
                cur_lineage_id_mapping[tree_id] = lineage_id
        return self.calc_lineage_cost(lineage_id_mapping=cur_lineage_id_mapping)

    def random_assignment_runner(self, iterations=1000000):
        random_assignment_costs = process_map(
            self.random_assignment_cost,
            range(iterations),
            max_workers=self.max_workers,
            chunksize=200,
            desc="Computing random assignment costs..."
        )
        return random_assignment_costs

    def random_assignment_cost(self, seed):
        random.seed(seed)
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        internal_lineage_ids = [cur_lineage_id_mapping[tree_id] for tree_id in self.internal_tree_ids]
        terminal_lineage_ids = [cur_lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        random.shuffle(internal_lineage_ids)
        random.shuffle(terminal_lineage_ids)
        for tree_id, lineage_id in zip(self.internal_tree_ids, internal_lineage_ids):
            cur_lineage_id_mapping[tree_id] = lineage_id
        for tree_id, lineage_id in zip(self.terminal_tree_ids, terminal_lineage_ids):
            cur_lineage_id_mapping[tree_id] = lineage_id
        return self.calc_lineage_cost(lineage_id_mapping=cur_lineage_id_mapping)
    
    def random_assignment_by_layer_runner(self, iterations=1000000):
        internal_tree_ids_by_depth = defaultdict(list)
        for depth in sorted(self.tree_ids_by_depth.keys()):
            for tree_id in self.tree_ids_by_depth[depth]:
                if tree_id in self.internal_tree_ids:
                    internal_tree_ids_by_depth[depth].append(tree_id)
        random_assignment_costs = process_map(
            partial(self.random_assignment_by_layer, internal_tree_ids_by_depth),
            range(iterations),
            max_workers=self.max_workers,
            chunksize=200,
            desc="Computing random assignment by layer costs..."
        )
        return random_assignment_costs

    def random_assignment_by_layer(self,internal_tree_ids_by_depth, seed):
        random.seed(seed)
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        for depth, tree_ids in internal_tree_ids_by_depth.items():
            lineage_ids = [cur_lineage_id_mapping[tree_id] for tree_id in tree_ids]
            random.shuffle(lineage_ids)
            for tree_id, lineage_id in zip(tree_ids, lineage_ids):
                cur_lineage_id_mapping[tree_id] = lineage_id
        terminal_lineage_ids = [cur_lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        for tree_id, lineage_id in zip(self.terminal_tree_ids, terminal_lineage_ids):
            cur_lineage_id_mapping[tree_id] = lineage_id
        return self.calc_lineage_cost(lineage_id_mapping=cur_lineage_id_mapping)


    def random_rebuild_runner(self, iterations=1000000, depth_limit=float('inf')):
        first_layer_tree_ids = [tree_id for tree_id, depth in self.first_internal_layer]
        first_layer_lineage_id_with_depth = [(self.lineage_tree.lineage_id_mapping[tree_id], depth)  
                                             for tree_id, depth in self.first_internal_layer]
        # gather non first layer internal lineage ids
        internal_candidate_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] 
                                          for tree_id in self.internal_tree_ids if tree_id not in first_layer_tree_ids]
        terminal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in self.terminal_tree_ids]
        random_rebuild_costs = process_map(
            partial(self.random_rebuild_cost, depth_limit, first_layer_lineage_id_with_depth,
                    internal_candidate_lineage_ids, terminal_lineage_ids),
            range(iterations),
            max_workers=self.max_workers,
            chunksize=200,
            desc="Computing random rebuild costs..."
        )
        return random_rebuild_costs

    def random_rebuild_cost(self, depth_limit, first_layer_lineage_id_with_depth,
                             internal_candidate_lineage_ids, terminal_lineage_ids, seed):
        random.seed(seed)
        shuffled_internal_ids = internal_candidate_lineage_ids.copy()
        random.shuffle(shuffled_internal_ids)
        parent_list = []
        xyz_cost = 0
        exp_cost = 0
        for lineage_id, depth in first_layer_lineage_id_with_depth:
            parent_list.append((lineage_id, depth))
            parent_list.append((lineage_id, depth))
        for lineage_id in shuffled_internal_ids:
            # randomly pick a parent from parent_list
            parent_idx = random.randint(0, len(parent_list)-1)
            parent_lineage_id, parent_depth = parent_list[parent_idx]
            while parent_depth >= depth_limit:
                parent_idx = random.randint(0, len(parent_list)-1)
                parent_lineage_id, parent_depth = parent_list[parent_idx]
            xyz_cost += self.xyz_cost_mat[lineage_id][parent_lineage_id]
            exp_cost += self.exp_cost_mat[lineage_id][parent_lineage_id]
            parent_list.pop(parent_idx)
            parent_list.append((lineage_id, parent_depth + 1))
            parent_list.append((lineage_id, parent_depth + 1))
        random.shuffle(parent_list)
        for idx, lineage_id in enumerate(terminal_lineage_ids):
            parent_lineage_id, parent_depth = parent_list[idx]
            xyz_cost += self.xyz_cost_mat[lineage_id][parent_lineage_id]
            exp_cost += self.exp_cost_mat[lineage_id][parent_lineage_id]
        return xyz_cost, exp_cost
    
    def estimate_phylogenetic_evolution_rate(self):
        terminal_xyz_values = self.xyz_mat[self.terminal_tree_ids].copy()
        terminal_exp_values = self.exp_mat[self.terminal_tree_ids].copy()
        terminal_values = np.hstack((terminal_xyz_values, terminal_exp_values))
        # terminal_x_values = list(terminal_xyz_values[:, 2])
        n_terminals = len(self.terminal_tree_ids)
        vcv_matrix = self.lineage_tree.compute_vcv_matrix(self.terminal_tree_ids, sigma_sq=1.0)
        vcv_inv = np.linalg.inv(vcv_matrix)
        ones = np.ones(n_terminals)
        root_state = (ones @ vcv_inv @ terminal_values) / (ones @ vcv_inv @ ones)
        centered_values = terminal_values - root_state
        sigma_sq_mle = (centered_values.T @ vcv_inv @ centered_values) / n_terminals
        return sigma_sq_mle.diagonal()
    
    def phylogenetic_reconstruction(self, root_id=0):
        epsilon = 1e-8
        terminal_xyz_values = self.xyz_mat[self.terminal_tree_ids].copy()
        terminal_exp_values = self.exp_mat[self.terminal_tree_ids].copy()
        terminal_values = np.hstack((terminal_xyz_values, terminal_exp_values))
        n_terminals = len(self.terminal_tree_ids)
        xyz_sigma_sq = self.estimate_phylogenetic_evolution_rate()
        node_estimates = {}
        node_variances = {}
        for i, tree_id in enumerate(self.terminal_tree_ids):
            node_estimates[tree_id] = terminal_values[i]
            node_variances[tree_id] = np.zeros(3+self.exp_mat.shape[1])

        def post_order_traversal(node_id):
            if self.lineage_tree.children_list[node_id] == []:
                return node_estimates[node_id], node_variances[node_id]
            child_estimates = []
            for child in self.lineage_tree.children_list[node_id]:
                if child not in node_estimates:
                    child_val, child_var = post_order_traversal(child)
                    node_estimates[child] = child_val
                    node_variances[child] = child_var
                child_estimates.append({
                    'value': node_estimates[child],
                    'variance': node_variances[child],
                    'branch_length': self.lineage_tree.branch_time_length[child]
                })
            # Calculate weighted average for this node
            weights = [1/(c['branch_length'] * xyz_sigma_sq + c['variance']) for c in child_estimates]
            values = [c['value'] for c in child_estimates]
            weighted_avg = sum(w * v for w, v in zip(weights, values)) / sum(weights)
            variance = 1 / sum(weights)

            return weighted_avg, variance

        root_val, root_var = post_order_traversal(root_id)
        node_estimates[root_id] = root_val
        node_variances[root_id] = root_var

        def pre_order_traversal(node_id, parent_value=None):
            if node_id == root_id:
                return
            if len(self.lineage_tree.children_list[node_id]) == 0:
                return
            upward_est = node_estimates[node_id]
            upward_var = node_variances[node_id]
            # Prediction from parent
            parent_pred = parent_value  # Under BM, expected value = parent value
            parent_pred_var = self.lineage_tree.branch_time_length[node_id] * xyz_sigma_sq

            # Combine estimates (inverse-variance weighting)
            w_upward = 1 / (upward_var + epsilon)
            w_parent = 1 / (parent_pred_var + epsilon)

            if not np.isinf(w_upward).any():
                refined_value = (upward_est * w_upward + parent_pred * w_parent) / (w_upward + w_parent)
                refined_var = 1 / (w_upward + w_parent)
            else:
                refined_value = upward_est
                refined_var = upward_var
            node_estimates[node_id] = refined_value
            node_variances[node_id] = refined_var

            for child in self.lineage_tree.children_list[node_id]:
                pre_order_traversal(child, refined_value)
        
        for child in self.lineage_tree.children_list[root_id]:
            pre_order_traversal(child, root_val)

        return node_estimates, node_variances
    
    def phylogenetic_sampling(self, root_id=0, n_samples=10000):
        node_estimates, node_variances = self.phylogenetic_reconstruction(root_id=root_id)
        samples = {}
        for tree_id in tqdm(node_estimates.keys(), desc="Sampling from phylogenetic estimates..."):
            mean = node_estimates[tree_id]
            var = node_variances[tree_id]
            samples[tree_id] = np.random.multivariate_normal(mean, np.diag(var), size=n_samples)
        xyz_cost = np.zeros((n_samples))
        exp_cost = np.zeros((n_samples))
        mean_xyz_cost = 0
        mean_exp_cost = 0
        for tree_id in tqdm(node_estimates.keys(), desc="Calculating phylogenetic costs..."):
            cur_val = samples[tree_id]
            parent_id = self.lineage_tree.parent_list[tree_id]
            if tree_id == root_id:
                continue
            parent_val = samples[parent_id]
            diff = cur_val - parent_val
            xyz_cost += np.linalg.norm(diff[:, :3], axis=1)
            exp_cost += np.linalg.norm(diff[:, 3:], axis=1)
            mean_xyz_cost += np.linalg.norm(node_estimates[tree_id][:3] - node_estimates[parent_id][:3])
            mean_exp_cost += np.linalg.norm(node_estimates[tree_id][3:] - node_estimates[parent_id][3:])
        return xyz_cost, exp_cost, mean_xyz_cost, mean_exp_cost

    def layerwise_pareto_assignment(self):
        """Bottom-up Hungarian assignment at each layer independently.

        At each depth level, runs Hungarian to find the optimal child→parent
        edge assignment. The per-layer Pareto fronts are aggregated into
        1001 full-tree parent mappings (one per alpha), from which total
        xyz/exp costs are computed via calc_lineage_cost.

        Returns (list_of_pareto_fronts, list_of_layer_costs, lineage_id_to_parent_list,
                 list_of_layer_assignments).
          - list_of_pareto_fronts[layer][k] = (xyz_cost, exp_cost) for alpha=k/1000
          - list_of_layer_costs[layer] = (xyz_cost, exp_cost) of the original (diagonal) assignment
          - lineage_id_to_parent_list[k][child_lid] = assigned parent lineage_id at alpha=k
          - list_of_layer_assignments[layer][k] = col_indices array from Hungarian
        """
        list_of_pareto_fronts = []
        list_of_layer_costs = []
        list_of_layer_assignments = []  # NEW: col_indices per alpha per layer
        current_layer = self.terminal_tree_ids.copy()
        first_layer_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
        seen_tree_ids = set(current_layer)
        lineage_id_to_parent_list = []
        for i in range(1001):
            lineage_id_to_parent_list.append([-1] * self.lineage_tree.size)
        while True:
            # if any(tree_id in first_layer_tree_ids for tree_id in current_layer):
            #     break
            print(len(current_layer), "nodes in current layer")
            current_layer_cost = (0, 0)
            current_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in current_layer]
            current_parent_tree_ids = [self.lineage_tree.parent_list[tree_id] for tree_id in current_layer]
            current_parent_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in current_parent_tree_ids]
            xyz_cost_mat = np.zeros((len(current_lineage_ids), len(current_parent_lineage_ids)))
            exp_cost_mat = np.zeros((len(current_lineage_ids), len(current_parent_lineage_ids)))
            for i in tqdm(range(len(current_lineage_ids))):
                candidate_lineage_id = current_lineage_ids[i]
                for j in range(len(current_parent_lineage_ids)):
                    parent_lineage_id = current_parent_lineage_ids[j]
                    xyz_cost_mat[i][j] += np.linalg.norm(self.xyz_mat[candidate_lineage_id] - self.xyz_mat[parent_lineage_id], ord=2)
                    exp_cost_mat[i][j] += np.linalg.norm(self.exp_mat[candidate_lineage_id] - self.exp_mat[parent_lineage_id], ord=self.exp_norm)
                current_layer_cost = (current_layer_cost[0] + xyz_cost_mat[i][i], current_layer_cost[1] + exp_cost_mat[i][i])
            list_of_layer_costs.append(current_layer_cost)

            current_pareto_front = []
            current_layer_assignments = []  # NEW: col_indices at each alpha for this layer
            for k in tqdm(range(1001)):
                alpha = 0 + 0.001 * k  # weight for xyz cost
                cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
                row_indices, col_indices = linear_sum_assignment(cost_mat)
                current_layer_assignments.append(col_indices.copy())  # NEW: record assignment
                total_xyz_cost = 0
                total_exp_cost = 0
                for row, col in zip(row_indices, col_indices):
                    total_xyz_cost += xyz_cost_mat[row][col]
                    total_exp_cost += exp_cost_mat[row][col]
                    child_lineage_id = current_lineage_ids[row]
                    parent_lineage_id = current_parent_lineage_ids[col]
                    lineage_id_to_parent_list[k][child_lineage_id] = parent_lineage_id
                current_pareto_front.append((total_xyz_cost, total_exp_cost))
            list_of_pareto_fronts.append(current_pareto_front)
            list_of_layer_assignments.append(current_layer_assignments)  # NEW
            # move to next layer
            next_layer = []
            for tree_id in current_parent_tree_ids:
                if tree_id in first_layer_tree_ids:
                    continue
                if tree_id in seen_tree_ids:
                    next_layer.append(tree_id)
                else:
                    seen_tree_ids.add(tree_id)
            if len(next_layer) == 0:
                break
            current_layer = next_layer

        return list_of_pareto_fronts, list_of_layer_costs, lineage_id_to_parent_list, list_of_layer_assignments

    def layerwise_pareto_assignment_runner(self):
        """Pareto front from per-layer Hungarian edge assignment.

        Calls layerwise_pareto_assignment and aggregates the per-layer
        parent mappings into full-tree costs by summing edge costs directly.
        Also computes per-layer edge retention and records per-layer Pareto fronts.

        Returns a dict with keys:
          - 'xyz', 'exp': aggregated Pareto front arrays (length 1001)
          - 'per_layer_pareto_fronts': list of per-layer (xyz, exp) arrays
          - 'per_layer_edge_retention': list of per-layer ER arrays (length 1001)
          - 'per_layer_costs': list of per-layer diagonal (original) costs
          - 'layer_sizes': list of node counts per layer
        """
        list_of_pareto_fronts, list_of_layer_costs, lineage_id_to_parent_list, \
            list_of_layer_assignments = self.layerwise_pareto_assignment()

        # ── Aggregated Pareto front (sum edge costs across all layers) ──
        pareto_xyz = np.zeros(1001)
        pareto_exp = np.zeros(1001)
        for i in range(1001):
            xyz_cost, exp_cost = 0, 0
            cur_mapping = lineage_id_to_parent_list[i]
            for child_lid, parent_lid in enumerate(cur_mapping):
                if parent_lid == -1:
                    continue
                xyz_cost += self.xyz_cost_mat[child_lid][parent_lid]
                exp_cost += self.exp_cost_mat[child_lid][parent_lid]
            pareto_xyz[i] = xyz_cost
            pareto_exp[i] = exp_cost

        # ── Per-layer edge retention ──
        # Edge retention at layer L, alpha k = fraction of children assigned
        # to their original (biological) parent (col_indices[i] == i)
        per_layer_er = []
        layer_sizes = []
        for layer_idx, assignments in enumerate(list_of_layer_assignments):
            er_curve = np.zeros(1001)
            n = len(assignments[0])
            layer_sizes.append(n)
            for k, col_indices in enumerate(assignments):
                er_curve[k] = np.sum(col_indices == np.arange(n)) / n
            per_layer_er.append(er_curve)

        # ── Per-layer Pareto fronts as arrays ──
        per_layer_pf = []
        for layer_pf in list_of_pareto_fronts:
            pf_xyz = np.array([p[0] for p in layer_pf])
            pf_exp = np.array([p[1] for p in layer_pf])
            per_layer_pf.append((pf_xyz, pf_exp))

        return {
            'xyz': pareto_xyz,
            'exp': pareto_exp,
            'per_layer_pareto_fronts': per_layer_pf,
            'per_layer_edge_retention': per_layer_er,
            'per_layer_costs': list_of_layer_costs,
            'layer_sizes': layer_sizes,
        }

    # Backward-compatible alias
    bottom_up_pareto_at_each_layer = layerwise_pareto_assignment

    def free_xyz_distance(self):
        print(np.nan_to_num(self.xyz_mat).mean())

        def local_energy(point, neighbors):
            # energy function for internal relaxation
            return sum(np.linalg.norm(point - n) for n in neighbors)

        first_layer_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
        new_xyz_dict = {}
        neighbor_ids = defaultdict(list)
        internal_tree_ids = set()
        for tree_id in self.terminal_tree_ids:
            lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
            lineage_xyz = self.xyz_mat[lineage_id].copy()
            new_xyz_dict[tree_id] = lineage_xyz
        cur_tree_ids = self.terminal_tree_ids.copy()
        while cur_tree_ids:
            next_layer_tree_ids = []
            for tree_id in cur_tree_ids:
                if tree_id in first_layer_tree_ids: continue
                lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
                lineage_xyz = self.xyz_mat[lineage_id].copy()
                parent_id = self.lineage_tree.parent_list[tree_id]
                internal_tree_ids.add(parent_id)
                neighbor_ids[parent_id].append(tree_id)
                neighbor_ids[tree_id].append(parent_id)
                if parent_id not in new_xyz_dict:
                    new_xyz_dict[parent_id] = lineage_xyz
                else:
                    new_xyz_dict[parent_id] += lineage_xyz
                    new_xyz_dict[parent_id] /= 2
                    next_layer_tree_ids.append(parent_id)
            cur_tree_ids = next_layer_tree_ids
        internal_tree_ids = list(internal_tree_ids)
        bfs_queue = deque(first_layer_tree_ids)
        total_distance = 0
        while bfs_queue:
            tree_id = bfs_queue.popleft()
            for child in self.lineage_tree.children_list[tree_id]:
                total_distance += np.linalg.norm(new_xyz_dict[tree_id] - new_xyz_dict[child])
                bfs_queue.append(child)
        unrelaxed_total_distance = total_distance
        print("Total free xyz distance:", unrelaxed_total_distance)

        for i in range(10):
            max_shift = 0
            for node in internal_tree_ids:
                # Gather neighbor positions
                if len(neighbor_ids[node]) < 2:
                    continue
                neighbors = [new_xyz_dict[n] for n in neighbor_ids[node]]
                res = minimize(local_energy, new_xyz_dict[node], args=(neighbors,), method='L-BFGS-B')
                shift = np.linalg.norm(new_xyz_dict[node] - res.x)
                new_xyz_dict[node] = res.x
                if shift > max_shift:
                    max_shift = shift

            bfs_queue = deque(first_layer_tree_ids)
            total_distance = 0
            while bfs_queue:
                tree_id = bfs_queue.popleft()
                for child in self.lineage_tree.children_list[tree_id]:
                    total_distance += np.linalg.norm(new_xyz_dict[tree_id] - new_xyz_dict[child])
                    bfs_queue.append(child)

            print(f"Iteration {i}: max shift = {max_shift}, total distance = {total_distance}")
            
        return total_distance, unrelaxed_total_distance, new_xyz_dict

    def top_down_free_xyz_distance(self):
        def local_energy(point, neighbors):
            # energy function for internal relaxation
            return sum(np.linalg.norm(point - n) for n in neighbors)
        first_layer_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
        descendant_ids = defaultdict(list)
        internal_tree_ids = set()
        neighbor_ids = defaultdict(set)
        def gather_descendants(tree_id, terminal_id):
            parent_tree_id = self.lineage_tree.parent_list[tree_id]
            descendant_ids[tree_id].append(terminal_id)
            if tree_id not in first_layer_tree_ids:
                internal_tree_ids.add(parent_tree_id)
                neighbor_ids[parent_tree_id].add(tree_id)
                neighbor_ids[tree_id].add(parent_tree_id)
                gather_descendants(parent_tree_id, terminal_id)
        for tree_id in self.terminal_tree_ids:
            gather_descendants(tree_id, tree_id)

        new_xyz_dict = {}
        for tree_id, descendants in descendant_ids.items():
            avg_xyz = np.zeros(3)
            for desc_id in descendants:
                desc_lineage_id = self.lineage_tree.lineage_id_mapping[desc_id]
                desc_xyz = self.xyz_mat[desc_lineage_id].copy()
                avg_xyz += desc_xyz
            avg_xyz /= len(descendants)
            new_xyz_dict[tree_id] = avg_xyz

        internal_tree_ids = list(internal_tree_ids)
        
        bfs_queue = deque(first_layer_tree_ids)
        total_distance = 0
        while bfs_queue:
            tree_id = bfs_queue.popleft()
            for child in self.lineage_tree.children_list[tree_id]:
                total_distance += np.linalg.norm(new_xyz_dict[tree_id] - new_xyz_dict[child])
                bfs_queue.append(child)
        unrelaxed_total_distance = total_distance
        print("Total free xyz distance:", unrelaxed_total_distance)

        for i in range(10):
            max_shift = 0
            for node in internal_tree_ids:
                # Gather neighbor positions
                if len(neighbor_ids[node]) < 2:
                    continue
                neighbors = [new_xyz_dict[n] for n in neighbor_ids[node]]
                res = minimize(local_energy, new_xyz_dict[node], args=(neighbors,), method='L-BFGS-B')
                shift = np.linalg.norm(new_xyz_dict[node] - res.x)
                new_xyz_dict[node] = res.x
                if shift > max_shift:
                    max_shift = shift

            bfs_queue = deque(first_layer_tree_ids)
            total_distance = 0
            while bfs_queue:
                tree_id = bfs_queue.popleft()
                for child in self.lineage_tree.children_list[tree_id]:
                    total_distance += np.linalg.norm(new_xyz_dict[tree_id] - new_xyz_dict[child])
                    bfs_queue.append(child)

            print(f"Iteration {i}: max shift = {max_shift}, total distance = {total_distance}")


        return unrelaxed_total_distance, total_distance, descendant_ids

    def bottom_up_by_layer_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list_layerwise = process_map(
                partial(self.bottom_up_by_layer, top_internal_tree_ids, self.internal_tree_ids, self.tree_ids_by_depth),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list_layerwise = process_map(
                partial(self.bottom_up_by_layer, top_internal_tree_ids, internal_tree_ids, tree_ids_by_depth),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        return pareto_list_layerwise

    def bottom_up_by_layer(self, first_internal_tree_ids: list[int], internal_tree_ids: list[int], tree_ids_by_depth: dict[int, list[int]], idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        optimization_pool = deepcopy(internal_tree_ids)
        optimization_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in optimization_pool]
        for depth in sorted(tree_ids_by_depth.keys(), reverse=True):
            children_tree_ids = tree_ids_by_depth[depth]
            parents_tree_ids = set(self.lineage_tree.parent_list[tree_id] for tree_id in children_tree_ids)
            parents_tree_ids = list(parents_tree_ids)
            xyz_cost_mat = np.zeros((len(optimization_pool), len(parents_tree_ids)))
            exp_cost_mat = np.zeros((len(optimization_pool), len(parents_tree_ids)))
            for i in range(len(optimization_pool)):
                for j in range(len(parents_tree_ids)):
                    parent_tree_id = parents_tree_ids[j]
                    candidate_lineage_id = optimization_pool[i]
                    for child_tree_id in self.lineage_tree.children_list[parent_tree_id]:
                        child_lineage_id = cur_lineage_id_mapping[child_tree_id]
                        xyz_cost_mat[i][j] += np.linalg.norm(self.xyz_mat[candidate_lineage_id] - self.xyz_mat[child_lineage_id])
                        exp_cost_mat[i][j] += np.linalg.norm(self.exp_mat[candidate_lineage_id] - self.exp_mat[child_lineage_id], ord=self.exp_norm)
            cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
            row_indices, col_indices = linear_sum_assignment(cost_mat)
            for row, col in zip(row_indices, col_indices):
                cur_lineage_id_mapping[parents_tree_ids[col]] = optimization_pool[row]
            for index in reversed(row_indices):
                del optimization_pool[index]
        xyz_cost, exp_cost = self.calc_lineage_cost(first_internal_tree_ids=first_internal_tree_ids, lineage_id_mapping=cur_lineage_id_mapping)
        return xyz_cost, exp_cost
    
    def top_down_by_layer_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            internal_tree_ids = self.internal_tree_ids
            terminal_tree_ids = self.terminal_tree_ids
            tree_ids_by_depth = self.tree_ids_by_depth
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
        return self.top_down_by_layer(top_internal_tree_ids, internal_tree_ids, terminal_tree_ids, tree_ids_by_depth, 1001)
        # pareto_list_layerwise = process_map(
        #     partial(self.top_down_by_layer, top_internal_tree_ids, internal_tree_ids, terminal_tree_ids, tree_ids_by_depth),
        #     range(1001),
        #     max_workers=self.max_workers,
        #     chunksize=20,
        #     desc="Computing pareto costs"
        # )
        return pareto_list_layerwise

    def top_down_by_layer(self, first_internal_tree_ids: list[int], internal_tree_ids: list[int], 
                          terminal_tree_ids: list[int], tree_ids_by_depth: dict[int, list[int]], idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        optimization_pool = deepcopy(internal_tree_ids)
        optimization_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in optimization_pool if tree_id not in first_internal_tree_ids]
        cur_xyz_cost = 0
        cur_exp_cost = 0
        for depth in sorted(tree_ids_by_depth.keys()):
            parent_tree_ids = tree_ids_by_depth[depth]
            children_tree_ids = []
            for parent_tree_id in parent_tree_ids:
                children_tree_ids.extend(self.lineage_tree.children_list[parent_tree_id])
            children_tree_ids = [tree_id for tree_id in children_tree_ids if tree_id not in terminal_tree_ids]
            xyz_cost_mat = np.zeros((len(children_tree_ids), len(optimization_pool)))
            exp_cost_mat = np.zeros((len(children_tree_ids), len(optimization_pool)))
            for i in range(len(children_tree_ids)):
                for j in range(len(optimization_pool)):
                    parent_lineage_id = cur_lineage_id_mapping[self.lineage_tree.parent_list[children_tree_ids[i]]]
                    candidate_lineage_id = optimization_pool[j]
                    xyz_cost_mat[i][j] += np.linalg.norm(self.xyz_mat[candidate_lineage_id] - self.xyz_mat[parent_lineage_id])
                    exp_cost_mat[i][j] += np.linalg.norm(self.exp_mat[candidate_lineage_id] - self.exp_mat[parent_lineage_id], ord=self.exp_norm)
            cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
            row_indices, col_indices = linear_sum_assignment(cost_mat)
            for row, col in zip(row_indices, col_indices):
                cur_lineage_id_mapping[children_tree_ids[row]] = optimization_pool[col]
                cur_xyz_cost += xyz_cost_mat[row][col]
                cur_exp_cost += exp_cost_mat[row][col]
            for index in sorted(col_indices, reverse=True):
                del optimization_pool[index]
            if not optimization_pool:
                print("finished internal nodes optimization")
                break
        terminal_parent_lineage_ids = []
        for tree_id in terminal_tree_ids:
            parent_tree_id = self.lineage_tree.parent_list[tree_id]
            terminal_parent_lineage_ids.append(cur_lineage_id_mapping[parent_tree_id])
        terminal_xyz_cost_mat = np.zeros((len(terminal_tree_ids), len(terminal_parent_lineage_ids)))
        terminal_exp_cost_mat = np.zeros((len(terminal_tree_ids), len(terminal_parent_lineage_ids)))
        for i, tree_id in enumerate(terminal_tree_ids):
            lineage_id = cur_lineage_id_mapping[tree_id]
            for j, parent_lineage_id in enumerate(terminal_parent_lineage_ids):
                terminal_xyz_cost_mat[i][j] = np.linalg.norm(self.xyz_mat[lineage_id] - self.xyz_mat[parent_lineage_id])
                terminal_exp_cost_mat[i][j] = np.linalg.norm(self.exp_mat[lineage_id] - self.exp_mat[parent_lineage_id], ord=self.exp_norm)
        row_indices, col_indices = linear_sum_assignment(alpha * terminal_xyz_cost_mat + (1 - alpha) * terminal_exp_cost_mat)
        for row, col in zip(row_indices, col_indices):
            cur_xyz_cost += terminal_xyz_cost_mat[row][col]
            cur_exp_cost += terminal_exp_cost_mat[row][col]
        return cur_xyz_cost, cur_exp_cost

    def bottom_up_by_cell_runner(self, first_internal_layer: list[tuple[int, int]] = None):


        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list = process_map(
                partial(self.bottom_up_by_cell, top_internal_tree_ids, self.internal_tree_ids, self.terminal_tree_ids),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list = process_map(
                partial(self.bottom_up_by_cell, top_internal_tree_ids, internal_tree_ids, terminal_tree_ids),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        return pareto_list

    def bottom_up_by_cell(self, first_internal_tree_ids: list[int], internal_tree_ids: list[int], terminal_tree_ids: list[int], idx: int):
        alpha = 0 + 0.001 * idx
        cur_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        optimization_pool = deepcopy(internal_tree_ids)
        optimization_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in optimization_pool]
        is_optimized = [False] * self.lineage_tree.size
        for tree_id in terminal_tree_ids: is_optimized[tree_id] = True
        optimization_openings = set()
        for tree_id in terminal_tree_ids:
            parent_tree_id = self.lineage_tree.parent_list[tree_id]
            to_add = True
            for child in self.lineage_tree.children_list[parent_tree_id]:
                if not is_optimized[child]:
                    to_add = False
                    break
            if to_add: optimization_openings.add(parent_tree_id)
        optimization_openings = list(optimization_openings)
        
        while optimization_pool:
            cur_best = float('inf')
            cur_best_lineage_id = -1
            cur_best_lineage_id_idx = -1
            cur_best_tree_id = -1
            cur_best_tree_id_idx = -1
            for lineage_id_idx, lineage_id in enumerate(optimization_pool):
                for tree_id_idx, tree_id in enumerate(optimization_openings):
                    cur_cost = 0
                    for child_id in self.lineage_tree.children_list[tree_id]:
                        child_lineage_id = cur_lineage_id_mapping[child_id]
                        cur_cost += cur_cost_mat[lineage_id][child_lineage_id]
                    cur_cost = cur_cost / len(self.lineage_tree.children_list[tree_id])
                    if cur_cost < cur_best:
                        cur_best = cur_cost
                        cur_best_lineage_id = lineage_id
                        cur_best_tree_id = tree_id
                        cur_best_lineage_id_idx = lineage_id_idx
                        cur_best_tree_id_idx = tree_id_idx
            if cur_best_lineage_id != -1:
                cur_lineage_id_mapping[cur_best_tree_id] = cur_best_lineage_id
                is_optimized[cur_best_tree_id] = True
                del optimization_pool[cur_best_lineage_id_idx]
                del optimization_openings[cur_best_tree_id_idx]
                # no need to check new openings if we optimized one node at the top
                if cur_best_tree_id in first_internal_tree_ids: continue
                next_parent = self.lineage_tree.parent_list[cur_best_tree_id]
                to_add = True
                for child in self.lineage_tree.children_list[next_parent]:
                    if not is_optimized[child]:
                        to_add = False
                        break
                if to_add: optimization_openings.append(next_parent)
        return self.calc_lineage_cost(first_internal_tree_ids=first_internal_tree_ids, lineage_id_mapping=cur_lineage_id_mapping)

    def top_down_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None, depth_weight_type: str | Callable = "static"):
        if first_internal_layer is None:
            first_internal_layer = [(tree_id, depth) for tree_id, depth in self.first_internal_layer]
            internal_tree_ids = self.internal_tree_ids
            terminal_tree_ids = self.terminal_tree_ids
        else:
            _, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
        # Sequential loop instead of process_map — avoids OOM from pickling
        # the full LineageOptimization object in multiprocessing workers.
        pareto_list = []
        for idx in tqdm(range(1001), desc="Computing pareto costs"):
            pareto_list.append(
                self.top_down_rebuild(first_internal_layer, internal_tree_ids,
                                      terminal_tree_ids, depth_weight_type, idx))
        return pareto_list

    def top_down_rebuild(self, first_internal_layer: list[tuple[int, int]], internal_tree_ids: list[int], terminal_tree_ids: list[int], depth_weight_type: str | Callable, idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        cur_children_list = [[] for _ in range(self.lineage_tree.size)]
        cur_parent_list = [-1] * self.lineage_tree.size
        # use the same lineage id as tree id in the rebuilt tree
        cur_lineage_id_mapping = np.arange(self.lineage_tree.size).tolist()
        internal_pool = []
        internal_top_nodes = [(self.lineage_tree.lineage_id_mapping[tree_id], depth) for tree_id, depth in first_internal_layer]
        top_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
        for tree_id in internal_tree_ids:
            if tree_id not in top_tree_ids:
                internal_pool.append(self.lineage_tree.lineage_id_mapping[tree_id])
        terminal_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        cur_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        cur_xyz_cost = 0
        cur_exp_cost = 0
        
        def weight_factor_by_depth(depth):
            if depth_weight_type == "linear":
                return depth
            elif depth_weight_type == "quadratic":
                return depth * depth
            elif depth_weight_type == "cubic":
                return depth * depth * depth
            elif depth_weight_type == "sqrt":
                return np.sqrt(depth)
            elif depth_weight_type == "log":
                return np.log(depth + 1)
            elif callable(depth_weight_type):
                return depth_weight_type(depth)
            else:
                return 1

        cost_pq = []
        for i, depth in internal_top_nodes:
            depth_weight = weight_factor_by_depth(depth)
            for j in internal_pool:
                heapq.heappush(cost_pq, (cur_cost_mat[i][j] * depth_weight, i, j, depth))
        # Heap safety: max entries = |top|*|pool| + |pool|*(|pool|+1)/2
        # Use 5× this as a generous safety bound to catch pathological growth.
        n_pool = len(internal_pool)
        n_top = len(internal_top_nodes)
        max_heap = n_top * n_pool + n_pool * (n_pool + 1) // 2
        max_heap = max(50000, int(max_heap * 5))
        max_iter = max(200000, max_heap * 5)
        iter_count = 0
        while internal_pool:
            iter_count += 1
            if iter_count > max_iter:
                print(f"Warning: top_down_rebuild stuck at alpha={alpha:.3f}, "
                      f"{len(internal_pool)} nodes remain, {len(cost_pq)} in heap, breaking")
                break
            if len(cost_pq) > max_heap:
                print(f"Warning: heap overflow at alpha={alpha:.3f}, "
                      f"heap={len(cost_pq)} > max={max_heap}, breaking")
                break
            if len(cost_pq) == 0:
                print("Error: cost_pq empty with pool remaining")
                break
            cur_cost, parent_id, child_id, depth = heapq.heappop(cost_pq)
            if child_id not in internal_pool:
                continue
            if len(cur_children_list[parent_id]) >= 2:
                continue
            cur_children_list[parent_id].append(child_id)
            cur_parent_list[child_id] = parent_id
            internal_pool.remove(child_id)
            cur_xyz_cost += self.xyz_cost_mat[parent_id][child_id]
            cur_exp_cost += self.exp_cost_mat[parent_id][child_id]
            new_depth = depth + 1
            depth_weight = weight_factor_by_depth(new_depth)
            for next_child in internal_pool:
                heapq.heappush(cost_pq, (cur_cost_mat[child_id][next_child] * depth_weight,
                                         child_id, next_child, new_depth))
        terminal_parent_list = []
        for tree_id in internal_tree_ids:
            lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
            children_count = len(cur_children_list[lineage_id])
            while children_count < 2:
                terminal_parent_list.append(lineage_id)
                children_count += 1
        terminal_cost_mat = np.zeros((len(terminal_pool), len(terminal_parent_list)))
        terminal_xyz_cost_mat = np.zeros((len(terminal_pool), len(terminal_parent_list)))
        terminal_exp_cost_mat = np.zeros((len(terminal_pool), len(terminal_parent_list)))
        for i, lineage_id in enumerate(terminal_pool):
            for j, parent_lineage_id in enumerate(terminal_parent_list):
                terminal_cost_mat[i][j] = cur_cost_mat[parent_lineage_id][lineage_id]
                terminal_xyz_cost_mat[i][j] = self.xyz_cost_mat[parent_lineage_id][lineage_id]
                terminal_exp_cost_mat[i][j] = self.exp_cost_mat[parent_lineage_id][lineage_id]
        row_indices, col_indices = linear_sum_assignment(terminal_cost_mat)
        for row, col in zip(row_indices, col_indices):
            lineage_id = terminal_pool[row]
            parent_lineage_id = terminal_parent_list[col]
            cur_children_list[parent_lineage_id].append(lineage_id)
            cur_parent_list[lineage_id] = parent_lineage_id
            cur_xyz_cost += terminal_xyz_cost_mat[row][col]
            cur_exp_cost += terminal_exp_cost_mat[row][col]

        lineage_id_by_depth = defaultdict(list)
        bfs_queue = deque()
        visited = set()
        for tree_id, depth in self.first_internal_layer:
            lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
            bfs_queue.append((lineage_id, depth))
            lineage_id_by_depth[depth].append(lineage_id)
            visited.add(lineage_id)
        while bfs_queue:
            lineage_id, depth = bfs_queue.popleft()
            for child in cur_children_list[lineage_id]:
                if child in visited:
                    continue  # cycle detected in greedy assignment — skip
                visited.add(child)
                lineage_id_by_depth[depth + 1].append(child)
                bfs_queue.append((child, depth + 1))
        
        # return cur_xyz_cost, cur_exp_cost, cur_children_list, cur_parent_list, max(list(lineage_id_by_depth.keys()))
        complexity_score = self.complexity_score(cur_parent_list)
        return cur_xyz_cost, cur_exp_cost, complexity_score

    def top_down_balanced_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list = process_map(
                partial(self.top_down_balanced_rebuild, top_internal_tree_ids, self.internal_tree_ids, \
                         self.terminal_tree_ids, self.xyz_cost_mat, self.exp_cost_mat, self.lineage_tree.lineage_id_mapping, self.lineage_tree.size),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            _, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list = process_map(
                partial(self.top_down_balanced_rebuild, top_internal_tree_ids, internal_tree_ids, \
                         terminal_tree_ids, self.xyz_cost_mat, self.exp_cost_mat, self.lineage_tree.lineage_id_mapping, self.lineage_tree.size),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        return pareto_list

    def top_down_balanced_rebuild(self, top_tree_ids: list[int], internal_tree_ids: list[int], terminal_tree_ids: list[int], \
                                  xyz_cost_mat, exp_cost_mat, lineage_id_mapping, lineage_tree_size, idx: int):
        
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        all_cost_mat = xyz_cost_mat * alpha + exp_cost_mat * (1 - alpha)
        cur_xyz_cost = 0
        cur_exp_cost = 0
        cur_children_list = [[] for _ in range(lineage_tree_size)]
        cur_parent_list = [-1] * lineage_tree_size
        # use the same lineage id as tree id in the rebuilt tree
        cur_lineage_id_mapping = np.arange(lineage_tree_size).tolist()
        top_layer = []
        internal_pool = []
        for tree_id in internal_tree_ids:
            lineage_id = lineage_id_mapping[tree_id]
            if tree_id in top_tree_ids:
                top_layer.append(lineage_id)
                top_layer.append(lineage_id)
            else:
                internal_pool.append(lineage_id)
        while internal_pool:
            cur_cost_mat = np.zeros((len(top_layer), len(internal_pool)))
            for i, lineage_id_i in enumerate(top_layer):
                for j, lineage_id_j in enumerate(internal_pool):
                    cur_cost_mat[i][j] = all_cost_mat[lineage_id_i][lineage_id_j]
            row_indices, col_indices = linear_sum_assignment(cur_cost_mat)
            next_top_layer = []
            for row, col in zip(row_indices, col_indices):
                parent_lineage_id = top_layer[row]
                child_lineage_id = internal_pool[col]
                cur_children_list[parent_lineage_id].append(child_lineage_id)
                cur_parent_list[child_lineage_id] = parent_lineage_id
                cur_xyz_cost += xyz_cost_mat[parent_lineage_id][child_lineage_id]
                cur_exp_cost += exp_cost_mat[parent_lineage_id][child_lineage_id]
                next_top_layer.append(child_lineage_id)
                next_top_layer.append(child_lineage_id)
            if len(top_layer) > len(internal_pool):
                for idx in row_indices[::-1]:
                    del top_layer[idx]
                next_top_layer.extend(top_layer)
            for idx in sorted(col_indices, reverse=True):
                del internal_pool[idx]
            top_layer = next_top_layer
        terminal_nodes = [lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        terminal_cost_mat = np.zeros((len(top_layer), len(terminal_nodes)))
        for i, lineage_id_i in enumerate(top_layer):
            for j, lineage_id_j in enumerate(terminal_nodes):
                terminal_cost_mat[i][j] = all_cost_mat[lineage_id_i][lineage_id_j]
        row_indices, col_indices = linear_sum_assignment(terminal_cost_mat)
        for row, col in zip(row_indices, col_indices):
            parent_lineage_id = top_layer[row]
            child_lineage_id = terminal_nodes[col]
            cur_children_list[parent_lineage_id].append(child_lineage_id)
            cur_parent_list[child_lineage_id] = parent_lineage_id
            cur_xyz_cost += xyz_cost_mat[parent_lineage_id][child_lineage_id]
            cur_exp_cost += exp_cost_mat[parent_lineage_id][child_lineage_id]
        return cur_xyz_cost, cur_exp_cost, cur_children_list, cur_parent_list


    def paired_bottom_up_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        with Manager() as manager:
            shared_data = manager.dict()
            shared_data = {
                'xyz_cost_mat': self.xyz_cost_mat,
                'exp_cost_mat': self.exp_cost_mat,
                'lineage_id_mapping': self.lineage_tree.lineage_id_mapping,
                'lineage_tree_size': self.lineage_tree.size
            }
            if first_internal_layer is None:
                top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
                internal_tree_ids = self.internal_tree_ids
                terminal_tree_ids = self.terminal_tree_ids
                
            else:
                top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
                tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list_first_quarter = process_map(
                partial(self.paired_bottom_up_rebuild, internal_tree_ids, terminal_tree_ids, \
                            shared_data["xyz_cost_mat"], shared_data["exp_cost_mat"], shared_data["lineage_id_mapping"], shared_data["lineage_tree_size"]),
                range(250),
                max_workers=self.max_workers,
                chunksize=5,
                desc="Computing pareto costs 1/4"
            )
            pareto_list_second_quarter = process_map(
                partial(self.paired_bottom_up_rebuild, internal_tree_ids, terminal_tree_ids, \
                            shared_data["xyz_cost_mat"], shared_data["exp_cost_mat"], shared_data["lineage_id_mapping"], shared_data["lineage_tree_size"]),
                range(250, 500),
                max_workers=self.max_workers,
                chunksize=5,
                desc="Computing pareto costs 2/4"
            )
            pareto_list_third_quarter = process_map(
                partial(self.paired_bottom_up_rebuild, internal_tree_ids, terminal_tree_ids, \
                            shared_data["xyz_cost_mat"], shared_data["exp_cost_mat"], shared_data["lineage_id_mapping"], shared_data["lineage_tree_size"]),
                range(500, 750),
                max_workers=self.max_workers,
                chunksize=5,
                desc="Computing pareto costs 3/4"
            )
            pareto_list_fourth_quarter = process_map(
                partial(self.paired_bottom_up_rebuild, internal_tree_ids, terminal_tree_ids, \
                            shared_data["xyz_cost_mat"], shared_data["exp_cost_mat"], shared_data["lineage_id_mapping"], shared_data["lineage_tree_size"]),
                range(750, 1001),
                max_workers=self.max_workers,
                chunksize=5,
                desc="Computing pareto costs 4/4"
            )
        return pareto_list_first_quarter + pareto_list_second_quarter + pareto_list_third_quarter + pareto_list_fourth_quarter

    def paired_bottom_up_rebuild(self, internal_tree_ids: list[int], terminal_tree_ids: list[int], \
                                  xyz_cost_mat, exp_cost_mat, lineage_id_mapping, lineage_tree_size, idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        all_cost_mat = xyz_cost_mat * alpha + exp_cost_mat * (1 - alpha)
        cur_xyz_cost = 0
        cur_exp_cost = 0
        cur_children_list = [[] for _ in range(lineage_tree_size)]
        cur_parent_list = [-1] * lineage_tree_size
        # use the same lineage id as tree id in the rebuilt tree
        cur_lineage_id_mapping = np.arange(lineage_tree_size).tolist()
        internal_pool = [lineage_id_mapping[tree_id] for tree_id in internal_tree_ids]
        bottom_layer = [lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        # to compute complexity score
        complexity_code_dict = deepcopy(self.terminal_type_code_dict)
        cell_div_count = 0
        parent_code_mapping = {}
        num_terminal_types = len(set(complexity_code_dict.values()))

        while internal_pool and len(bottom_layer) > 1:
            # Pre-allocate for expected number of edges
            n_combinations = len(bottom_layer) * (len(bottom_layer) - 1) // 2
            bottom_edges = []
            
            for a, b in combinations(bottom_layer, 2):
                bottom_edges.append((a, b, all_cost_mat[a][b]))
            
            G = nx.Graph()
            G.add_weighted_edges_from(bottom_edges)
            bottom_pairs = list(nx.min_weight_matching(G))
            
            # Clean up graph and edges immediately
            G.clear()
            del G
            del bottom_edges
            
            next_bottom_layer = []
            if len(bottom_layer) % 2 == 1:
                paired_children = set()
                for pair in bottom_pairs:
                    paired_children.update(pair)
                unpaired_child = (set(bottom_layer) - paired_children).pop()
                next_bottom_layer.append(unpaired_child)
            
            # Use more memory-efficient matrix indexing
            cur_cost_mat = all_cost_mat[np.ix_(internal_pool, [p[0] for p in bottom_pairs])]
            cur_cost_mat += all_cost_mat[np.ix_(internal_pool, [p[1] for p in bottom_pairs])]
            
            row_indices, col_indices = linear_sum_assignment(cur_cost_mat)
            
            for row, col in zip(row_indices, col_indices):
                parent_lineage_id = internal_pool[row]
                child1, child2 = bottom_pairs[col]
                cur_xyz_cost += xyz_cost_mat[parent_lineage_id][child1] + xyz_cost_mat[parent_lineage_id][child2]
                cur_exp_cost += exp_cost_mat[parent_lineage_id][child1] + exp_cost_mat[parent_lineage_id][child2]
                cur_children_list[parent_lineage_id].append(child1)
                cur_children_list[parent_lineage_id].append(child2)
                cur_parent_list[child1] = parent_lineage_id
                cur_parent_list[child2] = parent_lineage_id
                child1_code = complexity_code_dict.get(child1, -1)
                child2_code = complexity_code_dict.get(child2, -1)
                if child1_code > child2_code:
                    parent_code_tuple = (child2_code, child1_code)
                else:
                    parent_code_tuple = (child1_code, child2_code)
                if parent_code_tuple not in parent_code_mapping:
                    parent_code_mapping[parent_code_tuple] = num_terminal_types + len(parent_code_mapping)
                complexity_code_dict[parent_lineage_id] = parent_code_mapping[parent_code_tuple]
                cell_div_count += 1
                
                next_bottom_layer.append(parent_lineage_id)
            
            # Clean up current cost matrix before next iteration
            del cur_cost_mat
            
            for i in row_indices[::-1]:
                del internal_pool[i]
            
            bottom_layer = next_bottom_layer
            del bottom_pairs  # Clean up pairs list
        
        complexity_score = len(parent_code_mapping)/cell_div_count

        # Final cleanup
        del all_cost_mat

        return cur_xyz_cost, cur_exp_cost, complexity_score, parent_code_mapping
        # return cur_xyz_cost, cur_exp_cost, cur_children_list, cur_parent_list, complexity_score
    
    def direct_bottom_up_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            pareto_list = process_map(
                partial(self.direct_bottom_up_rebuild, self.internal_tree_ids, self.terminal_tree_ids),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            _, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list = process_map(
                partial(self.direct_bottom_up_rebuild, internal_tree_ids, terminal_tree_ids),
                range(1001),
                max_workers=self.max_workers,
                chunksize=20,
                desc="Computing pareto costs"
            )
        return pareto_list

    def direct_bottom_up_rebuild(self, internal_tree_ids: list[int], terminal_tree_ids: list[int], idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        all_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        cur_children_list = [[] for _ in range(self.lineage_tree.size)]
        cur_parent_list = [-1] * self.lineage_tree.size
        # use the same lineage id as tree id in the rebuilt tree
        cur_lineage_id_mapping = np.arange(self.lineage_tree.size).tolist()
        internal_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in internal_tree_ids]
        internal_pool = internal_pool + internal_pool # binary tree, each internal node can have two children
        bottom_layer = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        cur_xyz_cost = 0
        cur_exp_cost = 0
        while bottom_layer and internal_pool:
            # if len(internal_pool) < len(bottom_layer):
            #     print("Error: internal pool is smaller than bottom layer")
            cur_cost_mat = all_cost_mat[np.ix_(internal_pool, bottom_layer)]
            row_indices, col_indices = linear_sum_assignment(cur_cost_mat)
            next_bottom_layer = set()
            for row, col in zip(row_indices, col_indices):
                parent_lineage_id = internal_pool[row]
                next_bottom_layer.add(parent_lineage_id)
                child_lineage_id = bottom_layer[col]
                cur_children_list[parent_lineage_id].append(child_lineage_id)
                cur_parent_list[child_lineage_id] = parent_lineage_id
                cur_xyz_cost += self.xyz_cost_mat[parent_lineage_id][child_lineage_id]
                cur_exp_cost += self.exp_cost_mat[parent_lineage_id][child_lineage_id]
            for i in row_indices[::-1]:
                del internal_pool[i]
            bottom_layer = list(next_bottom_layer)
        # if bottom_layer:
        #     print(len(bottom_layer), alpha)
        #     extra_xyz_cost, extra_exp_cost = self.terminal_only_rebuild(bottom_layer, idx)
        #     cur_xyz_cost += extra_xyz_cost
        #     cur_exp_cost += extra_exp_cost
        return cur_xyz_cost, cur_exp_cost, cur_children_list, cur_parent_list, len(bottom_layer)
    
    def terminal_only_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None, use_kd_tree: bool=True, exp_replacement: bool=True):
        # some of the rows in exp_mat can have nan values
        # exclude those rows from kd tree
        if first_internal_layer is not None:
            _, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
        else:
            terminal_tree_ids = self.terminal_tree_ids
            internal_tree_ids = self.internal_tree_ids
        terminal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        internal_lineage_ids = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in internal_tree_ids]
        pareto_list_1 = process_map(
            partial(self.terminal_only_rebuild, terminal_lineage_ids, internal_lineage_ids, use_kd_tree, exp_replacement),
            range(500),
            max_workers=self.max_workers,
            chunksize=20,
            desc="Computing pareto costs 1/2"
        )
        pareto_list_2 = process_map(
            partial(self.terminal_only_rebuild, terminal_lineage_ids, internal_lineage_ids, use_kd_tree, exp_replacement),
            range(500, 1001),
            max_workers=self.max_workers,
            chunksize=20,
            desc="Computing pareto costs 2/2"
        )
        return pareto_list_1 + pareto_list_2

    def terminal_only_rebuild(self, terminal_lineage_ids: list[int], internal_lineage_ids: list[int], use_kd_tree: bool = False, exp_replacement: bool = True, idx: int = 0):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        cur_xyz_cost = 0
        cur_exp_cost = 0
        # terminal only, so relabeling everything start from 0
        # and add new imaginary nodes to the end as we rebuild the tree
        cur_xyz_mat = self.xyz_mat[terminal_lineage_ids].copy()  # Explicit copy
        cur_exp_mat = self.exp_mat[terminal_lineage_ids].copy()  # Explicit copy
        bottom_layer = list(range(len(terminal_lineage_ids)))  # Use list() instead of tolist()
        internal_exp_mat = self.exp_mat[internal_lineage_ids].copy()
        internal_selection_count = np.zeros(len(internal_lineage_ids), dtype=int).tolist() if use_kd_tree else None
        parent_code_mapping = {}
        type_code_list = [self.terminal_type_code_dict[lineage_id] for lineage_id in terminal_lineage_ids]
        cell_div_count = 0
        # record ancestry path for each terminal node
        ancestry_paths = [[i] for i in range(len(terminal_lineage_ids))]
        offspring_list = [[i] for i in range(len(terminal_lineage_ids))]

        while len(bottom_layer) > 1 and len(internal_exp_mat):
            # Pre-allocate arrays for better memory management
            n_combinations = len(bottom_layer) * (len(bottom_layer) - 1) // 2
            bottom_edges = []
            # bottom_edges.reserve = n_combinations  # Hint for list allocation
            
            for a, b in combinations(bottom_layer, 2):
                pair_cost = alpha * np.linalg.norm(cur_xyz_mat[a] - cur_xyz_mat[b]) \
                      + np.linalg.norm(cur_exp_mat[a] - cur_exp_mat[b], ord=self.exp_norm) * (1 - alpha)
                bottom_edges.append((a, b, pair_cost))
            
            G = nx.Graph()
            G.add_weighted_edges_from(bottom_edges)
            bottom_pairs = list(nx.min_weight_matching(G))
            
            # Clear graph to free memory
            G.clear()
            del G
            del bottom_edges  # Explicitly delete large list
            
            next_bottom_layer = []
            next_ancestry_paths = []
            if len(bottom_layer) % 2 == 1:
                paired_children = set()
                for pair in bottom_pairs:
                    paired_children.update(pair)
                unpaired_child = (set(bottom_layer) - paired_children).pop()
                next_bottom_layer.append(unpaired_child)
            
            # Pre-calculate new matrix size and pre-allocate
            new_nodes_count = len(bottom_pairs)
            
            # Use more memory-efficient approach for matrix expansion
            if new_nodes_count > 0:
                new_xyz_rows = np.empty((new_nodes_count, cur_xyz_mat.shape[1]))
                new_exp_rows = np.empty((new_nodes_count, cur_exp_mat.shape[1]))
                
                for i, (child1, child2) in enumerate(bottom_pairs):
                    mid_xyz = (cur_xyz_mat[child1] + cur_xyz_mat[child2]) * 0.5  # More efficient than /2
                    mid_exp = (cur_exp_mat[child1] + cur_exp_mat[child2]) * 0.5
                    exp_1 = cur_exp_mat[child1]
                    exp_2 = cur_exp_mat[child2]
                    # Find closest internal expression profile
                    if use_kd_tree:
                        dist_1 = np.linalg.norm(internal_exp_mat - exp_1, ord=self.exp_norm, axis=1)
                        dist_2 = np.linalg.norm(internal_exp_mat - exp_2, ord=self.exp_norm, axis=1)
                        combined_dist = dist_1 + dist_2
                        min_idx = np.argmin(combined_dist)
                        mid_exp = internal_exp_mat[min_idx]
                        internal_selection_count[min_idx] += 1
                        if not exp_replacement:
                            internal_exp_mat = np.delete(internal_exp_mat, min_idx, axis=0)
                            if len(internal_exp_mat) <= 0:
                                break

                    # Calculate costs before storing
                    cur_xyz_cost += np.linalg.norm(cur_xyz_mat[child1] - mid_xyz) + np.linalg.norm(cur_xyz_mat[child2] - mid_xyz)
                    cur_exp_cost += np.linalg.norm(cur_exp_mat[child1] - mid_exp, ord=self.exp_norm) + np.linalg.norm(cur_exp_mat[child2] - mid_exp, ord=self.exp_norm)
                    
                    new_xyz_rows[i] = mid_xyz
                    new_exp_rows[i] = mid_exp
                    next_bottom_layer.append(len(cur_xyz_mat) + i)
                    type_1 = type_code_list[child1]
                    type_2 = type_code_list[child2]
                    combined_offspring = offspring_list[child1] + offspring_list[child2]
                    offspring_list.append(combined_offspring)
                    for child in combined_offspring:
                        ancestry_paths[child].append(len(cur_xyz_mat) + i)
                    if type_1 > type_2:
                        parent_code_tuple = (type_2, type_1)
                    else:
                        parent_code_tuple = (type_1, type_2)
                    if parent_code_tuple not in parent_code_mapping:
                        parent_code_mapping[parent_code_tuple] = len(parent_code_mapping) + self.terminal_type_nums
                    parent_type = parent_code_mapping[parent_code_tuple]
                    type_code_list.append(parent_type)
                    cell_div_count += 1
                
                # Concatenate once instead of multiple vstack calls
                cur_xyz_mat = np.concatenate([cur_xyz_mat, new_xyz_rows], axis=0)
                cur_exp_mat = np.concatenate([cur_exp_mat, new_exp_rows], axis=0)
                
                # Clean up temporary arrays
                del new_xyz_rows, new_exp_rows
            
            bottom_layer = next_bottom_layer
            del bottom_pairs  # Clean up
        min_consines = []
        mean_cosines = []
        for path in ancestry_paths:
            cos_angles = []
            for i in range(1, len(path) -1):
                vec1 = cur_exp_mat[path[i]] - cur_exp_mat[path[i -1]]
                vec2 = cur_exp_mat[path[i +1]] - cur_exp_mat[path[i]]
                cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2) + 1e-10)
                cos_angles.append(cos_angle)
            min_consines.append(np.min(cos_angles) if cos_angles else 1.0)
            mean_cosines.append(np.mean(cos_angles) if cos_angles else 1.0)
        del cur_xyz_mat, cur_exp_mat  # Final cleanup
        complexity_score = len(parent_code_mapping) / cell_div_count if cell_div_count > 0 else 0
        return cur_xyz_cost, cur_exp_cost, internal_selection_count, complexity_score, min_consines, mean_cosines
    
    def mst_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            internal_tree_ids = self.internal_tree_ids
            terminal_tree_ids = self.terminal_tree_ids
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
        # Sequential loop — avoids multiprocessing pickle/OOM issues
        mst_stat_list = []
        for idx in tqdm(range(1001), desc="Computing MST costs"):
            mst_stat_list.append(
                self.mst_rebuild(top_internal_tree_ids, internal_tree_ids,
                                 terminal_tree_ids, idx))
        return mst_stat_list

    def mst_rebuild(self, top_internal_tree_ids, internal_tree_ids, terminal_tree_ids, idx):
        edges = []
        internal_nodes = [self.lineage_tree.lineage_id_mapping[i] for i in internal_tree_ids if i not in top_internal_tree_ids]
        all_leaf_nodes = [self.lineage_tree.lineage_id_mapping[id] for id in terminal_tree_ids + top_internal_tree_ids]
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        cur_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        for i in range(len(internal_nodes)):
            for j in range(i + 1, len(internal_nodes)):
                edges.append((internal_nodes[i], internal_nodes[j], cur_cost_mat[internal_nodes[i]][internal_nodes[j]]))
        G = nx.Graph()
        G.add_weighted_edges_from(edges)
        mst = nx.minimum_spanning_tree(G)
        leaf_openings = []
        for node, degree in mst.degree():
            if degree == 1:
                leaf_openings.append(node)
                leaf_openings.append(node)
            if degree == 2:
                leaf_openings.append(node)
        if len(leaf_openings) < len(all_leaf_nodes):
            print("Error: leaf openings less than required")
        assignment_cost_mat = np.zeros((len(leaf_openings), len(all_leaf_nodes)))
        for i in range(len(leaf_openings)):
            for j in range(len(all_leaf_nodes)):
                assignment_cost_mat[i][j] = cur_cost_mat[leaf_openings[i]][all_leaf_nodes[j]]
        row_indices, col_indices = linear_sum_assignment(assignment_cost_mat)
        debug_xyz_cost = 0; debug_exp_cost = 0
        for row, col in zip(row_indices, col_indices):
            debug_xyz_cost += self.xyz_cost_mat[leaf_openings[row]][all_leaf_nodes[col]]
            debug_exp_cost += self.exp_cost_mat[leaf_openings[row]][all_leaf_nodes[col]]
            mst.add_edge(leaf_openings[row], all_leaf_nodes[col], weight=assignment_cost_mat[row][col])
        if alpha in [0, 0.5, 1]:
            print(f"Debug: MST rebuilt at alpha={alpha:.3f}, ")
        return debug_xyz_cost, debug_exp_cost

        lineage_id_by_depth = defaultdict(list)
        bfs_queue = deque()
        searched_dict = defaultdict(int)
        for tree_id, depth in self.first_internal_layer:
            lineage_id = self.lineage_tree.lineage_id_mapping[tree_id]
            bfs_queue.append((lineage_id, depth))
            lineage_id_by_depth[depth].append(lineage_id)
            searched_dict[lineage_id] = 1
        while bfs_queue:
            lineage_id, depth = bfs_queue.popleft()
            if lineage_id not in mst:
                continue  # top-level node not in MST (skip-through tree)
            for neighbor in mst.neighbors(lineage_id):
                if searched_dict[neighbor] == 0:
                    searched_dict[neighbor] = 1
                    bfs_queue.append((neighbor, depth + 1))
                    lineage_id_by_depth[depth + 1].append(neighbor)
        
        cur_xyz_cost = 0
        cur_exp_cost = 0
        for u, v in mst.edges():
            cur_xyz_cost += self.xyz_cost_mat[u][v]
            cur_exp_cost += self.exp_cost_mat[u][v]

        return cur_xyz_cost, cur_exp_cost, max(list(lineage_id_by_depth.keys()))


    def bottom_up_iterative_assignment(self, internal_tree_ids: list[int], terminal_tree_ids: list[int], idx: int):
        alpha = 0 + 0.001 * idx  # weight for xyz cost
        all_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        cur_xyz_cost = 0
        cur_exp_cost = 0
        cur_children_list = [[] for _ in range(self.lineage_tree.size)]
        cur_parent_list = [-1] * self.lineage_tree.size
        # use the same lineage id as tree id in the rebuilt tree
        cur_lineage_id_mapping = np.arange(self.lineage_tree.size).tolist()
        internal_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in internal_tree_ids]
        bottom_layer = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in terminal_tree_ids]
        while internal_pool and len(bottom_layer) > 1:
            # Pre-allocate for expected number of edges
            n_combinations = len(bottom_layer) * (len(bottom_layer) - 1) // 2
            bottom_edges = []
            
            for a, b in combinations(bottom_layer, 2):
                bottom_edges.append((a, b, all_cost_mat[a][b]))
            
            G = nx.Graph()
            G.add_weighted_edges_from(bottom_edges)
            bottom_pairs = list(nx.min_weight_matching(G))
            
            # Clean up graph and edges immediately
            G.clear()
            del G
            del bottom_edges
            
            next_bottom_layer = []
            if len(bottom_layer) % 2 == 1:
                paired_children = set()
                for pair in bottom_pairs:
                    paired_children.update(pair)
                unpaired_child = (set(bottom_layer) - paired_children).pop()
                next_bottom_layer.append(unpaired_child)
            
            # Use more memory-efficient matrix indexing