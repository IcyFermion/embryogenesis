"""
Lineage optimization classes for embryogenesis analysis.

This module contains the LineageTree and LineageOptimization classes
extracted from the Jupyter notebook to resolve pickling issues with multiprocessing.
"""

import pandas as pd
import numpy as np
from collections import defaultdict, deque
from multiprocessing import Manager
from scipy.optimize import linear_sum_assignment
from tqdm.contrib.concurrent import process_map
import networkx as nx
import random
import heapq
import json
from copy import deepcopy
from utils import lineage_name_mapping, load_json, bidict
from zss import simple_distance, Node
from matplotlib import pyplot as plt
from typing import List, Tuple, Dict, Any
from functools import partial
from itertools import combinations


class LineageTree:
    """
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
        self.lineage_id_mapping.append(lineage_id)
        self.reverse_lineage_id_mapping[lineage_id] = current_id
        self.children_list.append([])
        self.parent_list.append(parent_id)
        if parent_id != -1:
            self.children_list[parent_id].append(current_id)
        self.size += 1
        return current_id

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
                 lineage_names: list[str]):
        self.xyz_mat = xyz_mat
        self.exp_mat = exp_mat
        self.lineage_tree = lineage_tree
        self.lineage_names = lineage_names
        self.first_internal_layer = first_internal_layer
        self.tree_ids_by_depth, self.terminal_tree_ids, self.internal_tree_ids = self.lineage_traverse()
        self.lineage_xyz_cost, self.lineage_exp_cost = self.calc_lineage_cost(debugging=True)
        self.xyz_cost_mat = np.zeros((self.lineage_tree.size, self.lineage_tree.size))
        self.exp_cost_mat = np.zeros((self.lineage_tree.size, self.lineage_tree.size))
        for i in range(self.lineage_tree.size):
            for j in range(self.lineage_tree.size):
                self.xyz_cost_mat[i, j] = np.linalg.norm(self.xyz_mat[i] - self.xyz_mat[j])
                self.exp_cost_mat[i, j] = np.linalg.norm(self.exp_mat[i] - self.exp_mat[j])

    def lineage_traverse(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            first_internal_layer = self.first_internal_layer
        tree_ids_by_depth = defaultdict(list)
        terminal_tree_ids = []
        internal_tree_ids = []
        bfs_queue = deque()
        for tree_id, depth in first_internal_layer:
            bfs_queue.append((tree_id, depth))
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
                exp_cost += np.linalg.norm(self.exp_mat[lineage_id] - self.exp_mat[child_lineage_id])
                bfs_queue.append(child)
                edge_count += 1
        if debugging:
            print(xyz_cost, exp_cost, edge_count)
        return xyz_cost, exp_cost

    def bottom_up_by_layer_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list_layerwise = process_map(
                partial(self.bottom_up_by_layer, top_internal_tree_ids, self.internal_tree_ids, self.tree_ids_by_depth),
                range(1001),
                max_workers=10,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list_layerwise = process_map(
                partial(self.bottom_up_by_layer, top_internal_tree_ids, internal_tree_ids, tree_ids_by_depth),
                range(1001),
                max_workers=10,
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
                        exp_cost_mat[i][j] += np.linalg.norm(self.exp_mat[candidate_lineage_id] - self.exp_mat[child_lineage_id])
            cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
            row_indices, col_indices = linear_sum_assignment(cost_mat)
            for row, col in zip(row_indices, col_indices):
                cur_lineage_id_mapping[parents_tree_ids[col]] = optimization_pool[row]
            for index in reversed(row_indices):
                del optimization_pool[index]
        return self.calc_lineage_cost(first_internal_tree_ids=first_internal_tree_ids, lineage_id_mapping=cur_lineage_id_mapping)

    def bottom_up_by_cell_runner(self, first_internal_layer: list[tuple[int, int]] = None):

        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list = process_map(
                partial(self.bottom_up_by_cell, top_internal_tree_ids, self.internal_tree_ids, self.terminal_tree_ids),
                range(1001),
                max_workers=10,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            tree_ids_by_depth, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list = process_map(
                partial(self.bottom_up_by_cell, top_internal_tree_ids, internal_tree_ids, terminal_tree_ids),
                range(1001),
                max_workers=10,
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
    
    def top_down_rebuild_runner(self, first_internal_layer: list[tuple[int, int]] = None):
        # top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
        # return self.top_down_rebuild(top_internal_tree_ids, self.internal_tree_ids, self.terminal_tree_ids, 0)
        if first_internal_layer is None:
            top_internal_tree_ids = [tree_id for tree_id, _ in self.first_internal_layer]
            pareto_list = process_map(
                partial(self.top_down_rebuild, top_internal_tree_ids, self.internal_tree_ids, self.terminal_tree_ids),
                range(1001),
                max_workers=10,
                chunksize=20,
                desc="Computing pareto costs"
            )
        else:
            top_internal_tree_ids = [tree_id for tree_id, _ in first_internal_layer]
            _, terminal_tree_ids, internal_tree_ids = self.lineage_traverse(first_internal_layer)
            pareto_list = process_map(
                partial(self.top_down_rebuild, top_internal_tree_ids, internal_tree_ids, terminal_tree_ids),
                range(1001),
                max_workers=10,
                chunksize=20,
                desc="Computing pareto costs"
            )
        return pareto_list
    
    def top_down_rebuild(self, first_internal_tree_ids: list[int], internal_tree_ids: list[int], terminal_tree_ids: list[int], idx: int):
        alpha = 0 + 0.001 * idx
        cur_cost_mat = self.xyz_cost_mat * alpha + self.exp_cost_mat * (1 - alpha)
        cur_lineage_id_mapping = deepcopy(self.lineage_tree.lineage_id_mapping)
        optimization_pool = deepcopy(internal_tree_ids)
        optimization_pool = [self.lineage_tree.lineage_id_mapping[tree_id] for tree_id in optimization_pool]
        is_optimized = [False] * self.lineage_tree.size
        for tree_id in first_internal_tree_ids: is_optimized[tree_id] = True
        optimization_openings = set(first_internal_tree_ids)
        
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
                        if not is_optimized[child_id]:
                            child_lineage_id = cur_lineage_id_mapping[child_id]
                            cur_cost += cur_cost_mat[lineage_id][child_lineage_id]
                    child_count = 0
                    for child_id in self.lineage_tree.children_list[tree_id]:
                        if not is_optimized[child_id]:
                            child_count += 1
                    if child_count > 0:
                        cur_cost = cur_cost / child_count
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
                # append new optimizable children to the openings
                for child_id in self.lineage_tree.children_list[cur_best_tree_id]:
                    if child_id not in terminal_tree_ids:
                        optimization_openings.append(child_id)
        return self.calc_lineage_cost(first_internal_tree_ids=first_internal_tree_ids, lineage_id_mapping=cur_lineage_id_mapping)

    def mst_test_runner(self):
        mst_stat_list = process_map(
            self.mst_test,
            range(1001),
            max_workers=10,
            chunksize=20,
            desc="Computing mst test"
        )
        return mst_stat_list

    def mst_test(self, idx):
        np.random.seed(idx)
        # sample points uniformly in an unit cube
        n_points = len(self.terminal_tree_ids)
        random_points = np.random.uniform(0, 1, (n_points, 3))
        # compute the minimum spanning tree
        G = nx.Graph()
        for i in range(n_points):
            for j in range(i+1, n_points):
                weight = np.linalg.norm(random_points[i] - random_points[j])
                G.add_edge(i, j, weight=weight)
        mst = nx.minimum_spanning_tree(G)
        total_weight = sum([G[u][v]['weight'] for u, v in mst.edges()])
        return total_weight

    def pareto_approximation_runner(self, optimization_strategy, n_runs=1001):
        if optimization_strategy == 'layer':
            return self.bottom_up_by_layer_runner()
        elif optimization_strategy == 'cell':
            return self.bottom_up_by_cell_runner()
        elif optimization_strategy == 'topdown':
            return self.top_down_rebuild_runner()
        else:
            raise ValueError("Invalid optimization strategy")

    def generate_tree_edit_distance_data(self):
        # generate random trees and compute edit distance between them
        pass

    def tree_edit_distance_runner(self):
        pass

    def tree_edit_distance(self):
        pass