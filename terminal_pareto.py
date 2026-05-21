import numpy as np
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map  # or thread_map
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon, cosine, minkowski
import seaborn as sns
from functools import partial
import random
from collections import defaultdict
from heapq import heapify, heappop, heappush
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
default_colors = sns.color_palette(n_colors=10)


class ParetoFront:
    def __init__(self, terminal_nodes, terminal_parents, terminal_ancestry, xyz_coord_dict, lineage_exp, binary_lineage_exp=None, binary_lineage_end_exp=None, title="Pareto Front"):
        self.terminal_nodes = terminal_nodes
        self.terminal_parents = terminal_parents
        self.terminal_ancestry = terminal_ancestry
        self.title = title

        self.xyz_coord_dict = xyz_coord_dict
        self.exp_vector_dict = {}
        self.binary_exp_vector_dict  = {}
        self.binary_end_exp_vector_dict = {}
        self.lca_distance_mat = self.terminal_parent_pairwise_distance()
        for node, parent in zip(terminal_nodes, terminal_parents):
            self.exp_vector_dict[node] = lineage_exp.loc[:, node].values
            self.exp_vector_dict[parent] = lineage_exp.loc[:, parent].values
            if binary_lineage_exp is not None and binary_lineage_end_exp is not None:
                self.binary_exp_vector_dict[node] = lineage_exp.loc[:, node].values
                self.binary_exp_vector_dict[parent] = lineage_exp.loc[:, parent].values
                self.binary_end_exp_vector_dict[node] = lineage_exp.loc[:, node].values
                self.binary_end_exp_vector_dict[parent] = lineage_exp.loc[:, parent].values
            else:
                self.binary_exp_vector_dict[node] = binary_lineage_exp.loc[:, node].values
                self.binary_exp_vector_dict[parent] = binary_lineage_exp.loc[:, parent].values
                self.binary_end_exp_vector_dict[node] = binary_lineage_end_exp.loc[:, node].values
                self.binary_end_exp_vector_dict[parent] = binary_lineage_end_exp.loc[:, parent].values

    def compute_cost_matrices(self, norm="l2", normalize=False):
        """
        Compute cost matrices for xyz coordinates and expression vectors.
        
        :param norm: The norm to use for distance calculation. Options are "l2", "cosine", "l2+cosine".
        :return: Tuple of (xyz_cost_mat, exp_cost_mat).
        """
        xyz_cost_mat = np.zeros((len(self.terminal_parents), len(self.terminal_nodes)))
        exp_cost_mat = np.zeros((len(self.terminal_parents), len(self.terminal_nodes)))
        for i in tqdm(range(len(self.terminal_parents))):
            for j in range(len(self.terminal_nodes)):
                if norm == "l1":
                    xyz_cost_mat[i][j] = minkowski(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]], p=1)
                    exp_cost_mat[i][j] = minkowski(self.exp_vector_dict[self.terminal_parents[i]], self.exp_vector_dict[self.terminal_nodes[j]], p=1)
                elif norm == "l2":
                    xyz_cost_mat[i][j] = minkowski(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]], p=2)
                    exp_cost_mat[i][j] = minkowski(self.exp_vector_dict[self.terminal_parents[i]], self.exp_vector_dict[self.terminal_nodes[j]], p=2)
                elif norm == "cosine":
                    xyz_cost_mat[i][j] = cosine(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]])
                    exp_cost_mat[i][j] = cosine(self.exp_vector_dict[self.terminal_parents[i]], self.exp_vector_dict[self.terminal_nodes[j]])
                elif norm == "l2+cosine":
                    xyz_cost_mat[i][j] = minkowski(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]], p=2)
                    exp_cost_mat[i][j] = cosine(self.exp_vector_dict[self.terminal_parents[i]], self.exp_vector_dict[self.terminal_nodes[j]])
                elif norm == "l2+shannon":
                    xyz_cost_mat[i][j] = minkowski(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]], p=2)
                    exp_cost_mat[i][j] = jensenshannon(self.binary_exp_vector_dict[self.terminal_parents[i]], self.binary_exp_vector_dict[self.terminal_nodes[j]])
                elif norm == "l2+binary_state":
                    xyz_cost_mat[i][j] = minkowski(self.xyz_coord_dict[self.terminal_parents[i]], self.xyz_coord_dict[self.terminal_nodes[j]], p=2)
                    exp_cost_mat[i][j] = (~self.binary_end_exp_vector_dict[self.terminal_parents[i]] & self.binary_exp_vector_dict[self.terminal_nodes[j]]).sum()

        if normalize:
            # rowwise normalize each cost matrix to have mean 0 and std 1
            xyz_cost_mat = (xyz_cost_mat - xyz_cost_mat.mean(axis=1, keepdims=True)) / xyz_cost_mat.std(axis=1, keepdims=True)
            exp_cost_mat = (exp_cost_mat - exp_cost_mat.mean(axis=1, keepdims=True)) / exp_cost_mat.std(axis=1, keepdims=True)

        return xyz_cost_mat, exp_cost_mat
    
    def lca_distance(self, ancestry1: list, ancestry2: list):
        min_length = min(len(ancestry1), len(ancestry2))
        max_length = max(len(ancestry1), len(ancestry2))
        lca_index = 0
        for i in range(min_length):
            if ancestry1[i] == ancestry2[i]:
                lca_index = i
            else:
                break
        return max_length -  (lca_index + 1)

    def monte_carlo_simulation(self, xyz_cost_mat, exp_cost_mat, _):
        row_ind = range(len(xyz_cost_mat))
        random_permutation = np.random.permutation(len(row_ind))
        return xyz_cost_mat[row_ind, random_permutation].sum(), exp_cost_mat[row_ind, random_permutation].sum()

    def terminal_cousin_group(self, degree=2):
        lowest_common_ancestor_group = defaultdict(list)
        for terminal_name, ancestry in self.terminal_ancestry.items():
            if len(ancestry) < degree:
                continue
            lca = ancestry[-degree]
            lowest_common_ancestor_group[lca].append(terminal_name)
        self.cousin_group = lowest_common_ancestor_group
        return lowest_common_ancestor_group
    
    def terminal_parent_pairwise_distance(self):
        distance_mat = np.zeros((len(self.terminal_parents), len(self.terminal_parents)))
        for i in range(len(self.terminal_parents)):
            for j in range(len(self.terminal_parents)):
                # find the distance between terminal_i and terminal_j in the lineage tree using their ancestry
                ancestry_i = self.terminal_ancestry[i]
                ancestry_j = self.terminal_ancestry[j]
                distance_mat[i][j] = self.lca_distance(ancestry_i, ancestry_j)

        return distance_mat



    def compute_random_cost(self, xyz_cost_mat, exp_cost_mat, _):
        new_xyz_cost = 0
        new_exp_cost = 0
        for group in self.cousin_group.values():
            if len(group) < 2:
                new_xyz_cost += xyz_cost_mat[group, group]
                new_exp_cost += exp_cost_mat[group, group]
            else:
                shuffled_group = random.sample(group, len(group))
                new_xyz_cost += xyz_cost_mat[group, shuffled_group].sum()
                new_exp_cost += exp_cost_mat[group, shuffled_group].sum()
        return new_xyz_cost, new_exp_cost
    
    def calculate_neighbor_change_distances(self, neighbor_dict, terminal_parent_pairwise_distances, children_dict, parent_dict):
        total_distances = 0
        for terminal_cell in parent_dict.keys():
            terminal_neighbors_to_parent_distances = []
            parent_neighbors_to_terminal_distances = []
            parent_cell = parent_dict[terminal_cell]
            if parent_cell not in neighbor_dict["terminal_parents"]:
                continue
            parent_neighbors = neighbor_dict["terminal_parents"][parent_cell]
            terminal_neighbors = neighbor_dict["terminal"][terminal_cell]
            for neighbor in parent_neighbors:

                if neighbor not in self.terminal_nodes:
                    if neighbor not in children_dict:
                        continue
                    neighbor_list = children_dict[neighbor]
                else:
                    neighbor_list = [neighbor]
                for node in neighbor_list:
                    distance = np.linalg.norm(self.xyz_coord_dict[node]- self.xyz_coord_dict[terminal_cell])
                    parent_neighbors_to_terminal_distances.append(distance)
            for neighbor in terminal_neighbors:
                if neighbor not in terminal_parent_pairwise_distances[parent_cell]:
                    neighbor = parent_dict[neighbor]
                    if neighbor not in terminal_parent_pairwise_distances[parent_cell]:
                        continue
                distance = terminal_parent_pairwise_distances[parent_cell][neighbor]
                terminal_neighbors_to_parent_distances.append(distance)
            avg_distance = (np.mean(terminal_neighbors_to_parent_distances) if terminal_neighbors_to_parent_distances else 0
                            + np.mean(parent_neighbors_to_terminal_distances) if parent_neighbors_to_terminal_distances else 0) / 2
            total_distances += avg_distance
        return total_distances

    def pareto_front_with_alt_exp(self, alt_exp_vector_dict, norm="l2", step_size=0.001):
        xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm=norm)
        alt_exp_cost_mat = np.zeros((len(self.terminal_parents), len(self.terminal_nodes)))
        for i in tqdm(range(len(self.terminal_parents))):
            for j in range(len(self.terminal_nodes)):
                if norm == "l1":
                    alt_exp_cost_mat[i][j] = minkowski(alt_exp_vector_dict[self.terminal_parents[i]], alt_exp_vector_dict[self.terminal_nodes[j]], p=1)
                elif norm == "l2":
                    alt_exp_cost_mat[i][j] = minkowski(alt_exp_vector_dict[self.terminal_parents[i]], alt_exp_vector_dict[self.terminal_nodes[j]], p=2)
                elif norm == "cosine":
                    alt_exp_cost_mat[i][j] = cosine(alt_exp_vector_dict[self.terminal_parents[i]], alt_exp_vector_dict[self.terminal_nodes[j]])
        loop_len = int(1 / step_size) + 1
        pareto_xyz_cost_list = []
        pareto_exp_cost_list = []
        pareto_alt_exp_cost_list = []
        for i in range(loop_len):
            alpha = 0 + step_size * i
            cur_cost_mat = alpha * xyz_cost_mat + (1 - alpha) * (alt_exp_cost_mat)
            # # add small cost to neighboring edges to discourage their assignment
            # for edge in neighbor_edges:
            #     parent_ind = self.terminal_parents.index(edge[0])
            #     node_ind = self.terminal_nodes.index(edge[1])
            #     cur_cost_mat[parent_ind][node_ind] += 1e6  # large penalty
            cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
            cur_xyz_cost = xyz_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_exp_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_alt_exp_cost = alt_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            pareto_xyz_cost_list.append(cur_xyz_cost)
            pareto_exp_cost_list.append(cur_exp_cost)
            pareto_alt_exp_cost_list.append(cur_alt_exp_cost)
        return pareto_xyz_cost_list, pareto_exp_cost_list, pareto_alt_exp_cost_list

    def pareto_front_with_neighbors(self, neighbor_edges, neighbor_dict, terminal_parent_pairwise_distances, step_size=0.001):

        xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm="l2")
        loop_len = int(1 / step_size) + 1
        pareto_xyz_cost_list = []
        pareto_exp_cost_list = []
        lineage_exp_cost = exp_cost_mat.diagonal().sum()
        lineage_xyz_cost = xyz_cost_mat.diagonal().sum()
        lineage_neighbor_distance = 0
        pareto_neighbor_distance_list = []
        pareto_neighbor_change_distance_list = []
        for i, j in neighbor_edges:
            lineage_neighbor_distance += self.lca_distance(self.terminal_ancestry[i], self.terminal_ancestry[j])
        for i in range(loop_len):
            alpha = 0 + step_size * i
            cur_cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
            # # add small cost to neighboring edges to discourage their assignment
            # for edge in neighbor_edges:
            #     parent_ind = self.terminal_parents.index(edge[0])
            #     node_ind = self.terminal_nodes.index(edge[1])
            #     cur_cost_mat[parent_ind][node_ind] += 1e6  # large penalty
            cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
            cur_xyz_cost = xyz_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_exp_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_neighbor_distance = 0
            for i, j in neighbor_edges:
                new_i = cur_col_ind[i]
                new_j = cur_col_ind[j]
                cur_neighbor_distance += self.lca_distance(self.terminal_ancestry[new_i], self.terminal_ancestry[new_j])
            cur_children_dict = defaultdict(list)
            cur_parent_dict = {self.terminal_nodes[cur_col_ind[i]]: self.terminal_parents[i] for i in range(len(cur_col_ind))}
            for terminal, parent in cur_parent_dict.items():
                cur_children_dict[parent].append(terminal)
            cur_neighbor_change_distance = self.calculate_neighbor_change_distances(neighbor_dict, terminal_parent_pairwise_distances, cur_children_dict, cur_parent_dict)
            pareto_neighbor_change_distance_list.append(cur_neighbor_change_distance)
            pareto_neighbor_distance_list.append(cur_neighbor_distance)
            pareto_xyz_cost_list.append(cur_xyz_cost)
            pareto_exp_cost_list.append(cur_exp_cost)
        return pareto_xyz_cost_list, pareto_exp_cost_list, pareto_neighbor_distance_list, lineage_xyz_cost, lineage_exp_cost, lineage_neighbor_distance, pareto_neighbor_change_distance_list

    def compute_pareto_front(self, xyz_cost_mat, exp_cost_mat, monte_carlo=False, step_size=0.001, return_lca_distances=False):
        loop_len = int(1 / step_size) + 1
        if monte_carlo:
            repeat = 1000000
            # multiprocessing to speed up monte carlo simulation
            monte_carlo_list = process_map(partial(self.monte_carlo_simulation, xyz_cost_mat, exp_cost_mat), range(repeat), max_workers=20, chunksize=100, desc="Monte Carlo Simulation")
            xyz_monte_carlo_mean = np.mean([x[0] for x in monte_carlo_list])
            exp_monte_carlo_mean = np.mean([x[1] for x in monte_carlo_list])
            xyz_monte_carlo_std = np.std([x[0] for x in monte_carlo_list]) + 1e-6  # add small value to avoid division by zero
            exp_monte_carlo_std = np.std([x[1] for x in monte_carlo_list]) + 1e-6
        else:
            xyz_monte_carlo_mean = 0
            exp_monte_carlo_mean = 0
            xyz_monte_carlo_std = 1
            exp_monte_carlo_std = 1


        pareto_xyz_cost_list = []
        pareto_exp_cost_list = []
        avg_lca_distances = []
        unchanged_assignment_ratio_list = []
        lineage_exp_cost = exp_cost_mat.diagonal().sum()
        lineage_xyz_cost = xyz_cost_mat.diagonal().sum()
        for i in range(loop_len):
            alpha = 0 + step_size * i
            cur_cost_mat = alpha * xyz_cost_mat + (1 - alpha) * exp_cost_mat
            if monte_carlo:
                cur_cost_mat = alpha * (xyz_cost_mat/xyz_monte_carlo_std) + (1 - alpha) * (exp_cost_mat/exp_monte_carlo_std)
            cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
            cur_xyz_cost = xyz_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_exp_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            avg_lca_distances.append(self.lca_distance_mat[cur_row_ind, cur_col_ind].sum() / len(self.terminal_parents))
            cur_xyz_cost_dev = (cur_xyz_cost - xyz_monte_carlo_mean) / xyz_monte_carlo_std
            cur_exp_cost_dev = (cur_exp_cost - exp_monte_carlo_mean) / exp_monte_carlo_std
            pareto_xyz_cost_list.append(cur_xyz_cost_dev)
            pareto_exp_cost_list.append(cur_exp_cost_dev)
            unchanged_assignment_count = 0
            for i, j in zip(cur_row_ind, cur_col_ind):
                if self.terminal_parents[i] == self.terminal_parents[j]:
                    unchanged_assignment_count += 1
            unchanged_assignment_ratio_list.append(unchanged_assignment_count / len(cur_row_ind))
        # normalize pareto cost lists to have a range of 0 to 1
        normalized_pareto_xyz_cost_list = (np.array(pareto_xyz_cost_list) - np.min(pareto_xyz_cost_list)) / (np.max(pareto_xyz_cost_list) - np.min(pareto_xyz_cost_list))
        normalized_pareto_exp_cost_list = (np.array(pareto_exp_cost_list) - np.min(pareto_exp_cost_list)) / (np.max(pareto_exp_cost_list) - np.min(pareto_exp_cost_list))
        normalized_lineage_xyz_cost = (lineage_xyz_cost - np.min(pareto_xyz_cost_list)) / (np.max(pareto_xyz_cost_list) - np.min(pareto_xyz_cost_list))
        normalized_lineage_exp_cost = (lineage_exp_cost - np.min(pareto_exp_cost_list)) / (np.max(pareto_exp_cost_list) - np.min(pareto_exp_cost_list))
        # find best alpha fit for to minimize distance to lineage costs in normalized space
        best_fit = 0
        min_distance = float('inf')
        min_xyz_distance = float('inf')
        for i in range(loop_len):
            alpha = 0 + step_size * i
            cur_distance = np.sqrt((normalized_pareto_xyz_cost_list[i] - normalized_lineage_xyz_cost) ** 2 + (normalized_pareto_exp_cost_list[i] - normalized_lineage_exp_cost) ** 2)
            if cur_distance < min_distance:
                min_distance = cur_distance
                min_xyz_distance = pareto_xyz_cost_list[i]
                best_fit = alpha
        print(best_fit, min_xyz_distance)

        if return_lca_distances:
            return best_fit, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost, avg_lca_distances, unchanged_assignment_ratio_list

        return best_fit, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost 
        
    def pareto_rescale(self, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost):
        """
        Rescale the pareto cost lists and lineage costs to have a range of 0 to 1.
        """
        min_xyz = np.min(pareto_xyz_cost_list)
        max_xyz = np.max(pareto_xyz_cost_list)
        min_exp = np.min(pareto_exp_cost_list)
        max_exp = np.max(pareto_exp_cost_list)

        pareto_xyz_cost_list = (np.array(pareto_xyz_cost_list) - min_xyz) / (max_xyz - min_xyz)
        pareto_exp_cost_list = (np.array(pareto_exp_cost_list) - min_exp) / (max_exp - min_exp)
        
        lineage_xyz_cost = (lineage_xyz_cost - min_xyz) / (max_xyz - min_xyz)
        lineage_exp_cost = (lineage_exp_cost - min_exp) / (max_exp - min_exp)

        return pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost
    

    def plot_pareto_front(self, lineage_xyz_cost, lineage_exp_cost, 
                      pareto_xyz_cost_list, pareto_exp_cost_list, 
                      alpha, norm, rescale=False, color='blue'):
        plt.figure(figsize=(8, 6))
        if rescale:
            pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost = self.pareto_rescale(
                pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost
            )
        plt.plot(pareto_xyz_cost_list, pareto_exp_cost_list, marker='o', linestyle='-', color=color, markersize=3, label='Pareto Front')
        plt.xlabel('Motility Cost')
        plt.ylabel('Expression Cost')
        plt.scatter(lineage_xyz_cost, lineage_exp_cost, marker="D", color=color, label='Lineage Cost')
        plt.title(f'{self.title} (alpha={alpha}, norm={norm})')
        plt.grid(True)
        plt.legend()
        plt.show()
        return plt.figure

    def plot_pareto_front_series(self, pareto_front_series):
        """
        Plot multiple pareto fronts in the same figure.
        Each entry in pareto_front_series should be a tuple of the form:
        (lineage_xyz_cost, lineage_exp_cost, pareto_xyz_cost_list, pareto_exp_cost_list, alpha, norm)
        """
        plt.figure(figsize=(8, 6))
        line_colors = sns.color_palette(n_colors=len(pareto_front_series))
        lineage_colors = sns.color_palette(n_colors=len(pareto_front_series))
        for pareto_entry, line_color, lineage_color in zip(pareto_front_series, line_colors, lineage_colors):
            # Unpack the pareto entry
            lineage_xyz_cost, lineage_exp_cost, pareto_xyz_cost_list, pareto_exp_cost_list, alpha, norm = pareto_entry
            # Rescale the axis to have a range of 0 to 1
            pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost = self.pareto_rescale(
                pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost
            )
            plt.plot(pareto_xyz_cost_list, pareto_exp_cost_list, marker='o', linestyle='-', color=line_color, markersize=3, label=f'{norm} norm')
            plt.scatter(lineage_xyz_cost, lineage_exp_cost, marker="D" ,color=lineage_color, label=f'{norm} norm lineage')
        plt.xlabel('Motility Cost')
        plt.ylabel('Expression Cost')
        plt.title(f'{self.title}s for Different Norms')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_dual_axis_jointplot(self, x_data, y1_data, y2_data, x_label, y1_label, y2_label, title, x_name, y1_name, y2_name):
        """
        Create a joint plot with dual y-axes showing two comparisons against the same x variable.
        Uses seaborn styling for consistency with other plots.
        """
        # Convert to arrays for easier manipulation
        x_data = np.array(x_data)
        y1_data = np.array(y1_data)
        y2_data = np.array(y2_data)
        
        # Create figure with seaborn style
        sns.set_style("whitegrid")
        fig = plt.figure(figsize=(12, 8))
        
        # Create grid layout
        gs = fig.add_gridspec(3, 3, width_ratios=[1, 3, 1], height_ratios=[1, 3, 1], 
                             hspace=0.05, wspace=0.05)
        
        # Main plot with dual y-axes
        ax_main = fig.add_subplot(gs[1, 1], zorder=5)
        ax_main2 = ax_main.twinx()
        
        # Use seaborn color palette
        colors = sns.color_palette(n_colors=2)
        color1, color2 = colors[0], colors[1]
        
        # Scatter plots with seaborn styling
        ax_main.scatter(x_data, y1_data, alpha=0.5, color=color1, s=30, 
                       edgecolors='white', linewidth=0.5, label=y1_label)
        ax_main2.scatter(x_data, y2_data, alpha=0.5, color=color2, s=30, 
                        edgecolors='white', linewidth=0.5, label=y2_label)
        
        # Regression lines using seaborn's regplot approach
        sns.regplot(x=x_data, y=y1_data, ax=ax_main, scatter=False, 
                   color=color1, line_kws={'linewidth': 2, 'alpha': 0.8})
        sns.regplot(x=x_data, y=y2_data, ax=ax_main2, scatter=False, 
                   color=color2, line_kws={'linewidth': 2, 'alpha': 0.8})
        
        # Labels and styling
        ax_main.set_xlabel(x_label, fontsize=12)
        ax_main.set_ylabel(y1_label, color=color1, fontsize=12)
        ax_main2.set_ylabel(y2_label, color=color2, fontsize=12)
        ax_main2.yaxis.set_label_coords(0.95,0.5)
        ax_main.tick_params(labelleft=False, labelright=False)
        ax_main2.tick_params(labelleft=False, labelright=False)
        
        # Marginal distributions using seaborn
        ax_top = fig.add_subplot(gs[0, 1], sharex=ax_main)
        ax_right1 = fig.add_subplot(gs[1, 0], sharey=ax_main)
        ax_right2 = fig.add_subplot(gs[1, 2], sharey=ax_main2)
        
        # Top marginal (x distribution) - using seaborn
        sns.histplot(x=x_data, ax=ax_top, bins=40, color='gray', 
                    kde=True, alpha=0.6, stat='density')
        ax_top.set_ylabel('Density', fontsize=10)
        ax_top.set_xlabel('')
        ax_top.tick_params(labelbottom=False)
        
        # Right marginals (y distributions) - using seaborn
        sns.histplot(y=y1_data, ax=ax_right1, bins=40, color=color1, 
                    kde=True, alpha=0.6, stat='density')
        ax_right1.set_xlabel('Density', fontsize=10)
        ax_right1.set_ylabel('')
        ax_right1.tick_params(labelleft=False)
        
        sns.histplot(y=y2_data, ax=ax_right2, bins=40, color=color2, 
                    kde=True, alpha=0.6, stat='density')
        ax_right2.set_xlabel('Density', fontsize=10)
        ax_right2.set_ylabel('')
        ax_right2.tick_params(labelright=False)
        ax_right2.yaxis.set_label_position('left')
        ax_right2.yaxis.tick_left()
        
        # some stats:
        spearman_cur_comp1 = spearmanr(x_data, y1_data)
        spearman_cur_comp2 = spearmanr(x_data, y2_data)
        stat1 = (f" \n Rank cor between {x_name} and {y1_name}: {"{:.2}".format(spearman_cur_comp1.correlation)}, p-value: {"{:.2e}".format(spearman_cur_comp1.pvalue)}")
        stat2 = (f" \n Rank cor between {x_name} and {y2_name}: {"{:.2}".format(spearman_cur_comp2.correlation)}, p-value: {"{:.2e}".format(spearman_cur_comp2.pvalue)}")

        # Add title
        plt.suptitle(title+stat1+stat2, y=0.98, fontsize=12, fontweight='bold')
        
        # Add legends
        lines1, labels1 = ax_main.get_legend_handles_labels()
        lines2, labels2 = ax_main2.get_legend_handles_labels()
        ax_main.legend(lines1 + lines2, labels1 + labels2, loc='upper left', framealpha=0.9)
        
        plt.tight_layout()
        plt.show()


    def restricted_shuffle(self, norm="l2"):
        xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm=norm)
        self.terminal_cousin_group()
        # Use multiprocessing for speed up
        random_cost_list = process_map(
            partial(self.compute_random_cost, xyz_cost_mat, exp_cost_mat),
            range(100000), 
            max_workers=20, 
            chunksize=1000, 
            desc="Computing random costs"
        )

        # plotting the distribution of random costs in a 2d plot
        lineage_xyz_cost = xyz_cost_mat.diagonal().sum()
        lineage_exp_cost = exp_cost_mat.diagonal().sum()
        alpha, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost = self.compute_pareto_front(xyz_cost_mat, exp_cost_mat, monte_carlo=False)

        plt.figure(figsize=(8, 6))
        plt.plot(pareto_xyz_cost_list, pareto_exp_cost_list, marker='o', linestyle='-', color='blue', markersize=3, label='Pareto Front')
        plt.xlabel('Motility Cost')
        plt.ylabel('Expression Cost')
        plt.scatter(lineage_xyz_cost, lineage_exp_cost, marker='*', color='white', label='Real Lineage Cost')
        plt.scatter([x[0] for x in random_cost_list], [x[1] for x in random_cost_list], alpha=0.5, s=10, color=default_colors[8], label='Random Cousin Shuffles')
        plt.title(f'{self.title}')
        # plt.title(f'{self.title} (alpha={alpha}, norm={norm})')
        plt.grid(True)
        plt.legend()
        plt.show()

    def optimal_xyz_assignment(self, norm="l2"):
        xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm=norm)
        cur_cost_mat = xyz_cost_mat
        cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
        print(cur_row_ind, cur_col_ind)
        cur_xyz_cost = xyz_cost_mat[cur_row_ind, cur_col_ind].sum()
        cur_exp_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
        print(f"Optimal xyz assignment costs (norm={norm}): Motility Cost = {cur_xyz_cost}, Expression Cost = {cur_exp_cost}")
        return cur_row_ind, cur_col_ind, cur_xyz_cost, cur_exp_cost

    def pre_plot(self, norm="l2"):
        # plt.style.use('dark_background')
        plt.style.use('default')
        xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm=norm)
        
        # optimal xyz assignment
        cur_cost_mat = xyz_cost_mat
        cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
        cur_xyz_cost_list = xyz_cost_mat[cur_row_ind, cur_col_ind]
        cur_xyz_cost_sum = cur_xyz_cost_list.sum()
        lineage_xyz_cost_list = xyz_cost_mat[cur_row_ind, cur_row_ind]
        xyz_cost_diff_list = cur_xyz_cost_list - lineage_xyz_cost_list
        cur_exp_cost = exp_cost_mat[cur_row_ind, cur_col_ind]
        cur_exp_cost_sum = cur_exp_cost.sum()
        # plot the distribution of optimal xyz costs vs lineage xyz costs
        plt.figure(figsize=(8, 6))
        plt.scatter(lineage_xyz_cost_list, cur_xyz_cost_list, color='yellow', label='Optimal XYZ Assignment Costs', s=10)
        # draw y=x line for reference
        max_cost = max(max(lineage_xyz_cost_list), max(cur_xyz_cost_list))
        plt.plot([0, max_cost], [0, max_cost], color='red', linestyle='--', zorder=-1, label='y=x Reference Line')
        plt.xlabel('Lineage Displacement (um)')
        plt.ylabel('Optimally Assigned Displacement (um)')
        # plt.title(f'Optimal XYZ Assignment vs Lineage Costs (norm={norm})')
        plt.show()
        print(f"Optimal xyz assignment costs (norm={norm}): Motility Cost = {cur_xyz_cost_sum}, Expression Cost = {cur_exp_cost_sum}")
        print("number of optimally assigned displacements stayed the same:", sum(cur_xyz_cost_list == lineage_xyz_cost_list))
        top5_ind = np.argsort(cur_xyz_cost_list - lineage_xyz_cost_list)[:5]
        print("Top 5 Terminal Nodes with Largest Decrease in Displacement Cost:")
        for ind in top5_ind:
            terminal_node_ind = cur_col_ind[ind]
            print(f"Terminal Node: {self.terminal_nodes[terminal_node_ind]}, Lineage Parent: {self.terminal_parents[terminal_node_ind]}, Optimal Parent: {self.terminal_parents[ind]}, Displacement Difference: {xyz_cost_diff_list[ind]:.4f}")
    
        # optimal exp assignment
        cur_cost_mat = exp_cost_mat
        cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
        cur_exp_cost_list = exp_cost_mat[cur_row_ind, cur_col_ind]
        cur_exp_cost_sum = cur_exp_cost_list.sum()
        lineage_exp_cost_list = exp_cost_mat[cur_row_ind, cur_row_ind]
        exp_cost_diff_list = cur_exp_cost_list - lineage_exp_cost_list
        cur_xyz_cost = xyz_cost_mat[cur_row_ind, cur_col_ind]
        cur_xyz_cost_sum = cur_xyz_cost.sum()
        # plot the distribution of optimal xyz costs vs lineage xyz costs
        plt.figure(figsize=(8, 6))
        plt.scatter(lineage_exp_cost_list, cur_exp_cost_list, color='yellow', label='Optimal Expression Assignment Costs', s=10)
        # draw y=x line for reference
        max_cost = max(max(lineage_exp_cost_list), max(cur_exp_cost_list))
        plt.plot([0, max_cost], [0, max_cost], color='red', linestyle='--', zorder=-1, label='y=x Reference Line')
        plt.xlabel('Lineage Expression Cost')
        plt.ylabel('Optimally Assigned Expression Cost')
        # plt.title(f'Optimal Expression Assignment vs Lineage Costs (norm={norm})')
        plt.show()
        print(f"Optimal Expression assignment costs (norm={norm}): Expression Cost = {cur_exp_cost_sum}, Motility Cost = {cur_xyz_cost_sum}")
        print("number of optimally assigned displacements stayed the same:", sum(cur_exp_cost_list == lineage_exp_cost_list))
        top5_ind = np.argsort(exp_cost_diff_list)[:5]
        print("Top 5 Terminal Nodes with Largest Decrease in Expression Cost:")
        for ind in top5_ind:
            terminal_node_ind = cur_col_ind[ind]
            print(f"Terminal Node: {self.terminal_nodes[terminal_node_ind]}, Lineage Parent: {self.terminal_parents[terminal_node_ind]}, Optimal Parent: {self.terminal_parents[ind]}, Displacement Difference: {exp_cost_diff_list[ind]:.4f}")

        monte_carlo_list = process_map(partial(self.monte_carlo_simulation, xyz_cost_mat, exp_cost_mat), range(1000000), max_workers=20, chunksize=100, desc="Monte Carlo Simulation")
        self.terminal_cousin_group()
        # Use multiprocessing for speed up
        random_cost_list_1st_cousins = process_map(
            partial(self.compute_random_cost, xyz_cost_mat, exp_cost_mat),
            range(100000), 
            max_workers=20, 
            chunksize=1000, 
            desc="Computing random costs by 1st cousins shuffle"
        )

        self.terminal_cousin_group(degree=3)
        # Use multiprocessing for speed up
        random_cost_list_2nd_cousins = process_map(
            partial(self.compute_random_cost, xyz_cost_mat, exp_cost_mat),
            range(100000), 
            max_workers=20, 
            chunksize=1000, 
            desc="Computing random costs by 2nd cousins shuffle"
        )

        self.terminal_cousin_group(degree=4)
        # Use multiprocessing for speed up
        random_cost_list_3rd_cousins = process_map(
            partial(self.compute_random_cost, xyz_cost_mat, exp_cost_mat),
            range(100000), 
            max_workers=20, 
            chunksize=1000, 
            desc="Computing random costs by 3rd cousins shuffle"
        )

        # plotting the distribution of random costs in a 2d plot
        lineage_xyz_cost = xyz_cost_mat.diagonal().sum()
        lineage_exp_cost = exp_cost_mat.diagonal().sum()
        alpha, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, \
            lineage_exp_cost, avg_lca_distances, unchanged_assignment_ratio_list \
                = self.compute_pareto_front(xyz_cost_mat, exp_cost_mat, monte_carlo=False, return_lca_distances=True)

        print("lineage costs:", "xyz,", lineage_xyz_cost, "exp,", lineage_exp_cost)
        print("bottom right corner:", pareto_xyz_cost_list[0], pareto_exp_cost_list[0])
        print("top left corner:", pareto_xyz_cost_list[-1], pareto_exp_cost_list[-1])
        monte_carlo_xyz_list = [x[0] for x in monte_carlo_list]
        monte_carlo_exp_list = [x[1] for x in monte_carlo_list]
        print("monte carlo costs:", "xyz,", np.mean(monte_carlo_xyz_list), "exp,", np.mean(monte_carlo_exp_list))
        print("monte carlo costs std:", "xyz,", np.std(monte_carlo_xyz_list), "exp,", np.std(monte_carlo_exp_list))
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 10), sharex=True)
        ax1.plot(pareto_xyz_cost_list, pareto_exp_cost_list, marker='o', linestyle='-', color='black', markersize=3, label='Pareto Front')
        ax1.set_ylabel('Expression Cost')
        ax1.scatter(lineage_xyz_cost, lineage_exp_cost, marker='*', color='#FF0000', label='Real Lineage Cost')
        ax1.scatter([x[0] for x in random_cost_list_1st_cousins], [x[1] for x in random_cost_list_1st_cousins], alpha=0.5, s=10, color="#B59F00", label='Random 1st Cousin Shuffles')
        ax1.scatter([x[0] for x in random_cost_list_2nd_cousins], [x[1] for x in random_cost_list_2nd_cousins], alpha=0.5, s=10, color="#FF8C42", label='Random 2nd Cousin Shuffles')
        ax1.scatter([x[0] for x in random_cost_list_3rd_cousins], [x[1] for x in random_cost_list_3rd_cousins], alpha=0.5, s=10, color="#20B2AA", label='Random 3rd Cousin Shuffles')
        # ax1.scatter([x[0] for x in monte_carlo_list], [x[1] for x in monte_carlo_list], alpha=0.5, s=10, color="#00D060", label='Random Shuffles')
        # ax1.set_title(f'{self.title} (alpha={alpha}, norm={norm})')
        # ax1.grid(True)
        # ax1.legend()
        ax2.plot(pareto_xyz_cost_list, avg_lca_distances, marker='o', linestyle='-', color='black', markersize=3, label='Avg LCA Distance')
        # mark the xyz cost of the real lineage
        ax2.plot(lineage_xyz_cost, 0, color='#FF0000', marker='*', label='Real Lineage XYZ Cost')
        # mark the lowest average lca distance point on the plot
        min_lca_distance_index = np.argmin(avg_lca_distances)
        ax2.scatter(pareto_xyz_cost_list[min_lca_distance_index], avg_lca_distances[min_lca_distance_index], marker='d', color='#FFEA00', s=50, label='Lowest Avg LCA Distance Point', zorder=10)
        ax2.set_xlabel('Motility Cost')
        ax2.set_ylabel('Avg LCA Distance of Assigned Parents')
        ax2.set_ylim(0, 10)
        ax3.plot(pareto_xyz_cost_list, unchanged_assignment_ratio_list, marker='o', linestyle='-', color='black', markersize=3, label='Unchanged Assignment Ratio')
        # mark the xyz cost of the real lineage
        ax3.plot(lineage_xyz_cost, 1, color='#FF0000', marker='*', label='Real Lineage XYZ Cost')
        # mark the highest unchanged assignment ratio point on the plot
        max_unchanged_assignment_index = np.argmax(unchanged_assignment_ratio_list)
        ax3.scatter(pareto_xyz_cost_list[max_unchanged_assignment_index], unchanged_assignment_ratio_list[max_unchanged_assignment_index], marker='d', color='#FFEA00', s=50, label='Highest Unchanged Assignment Ratio Point', zorder=10)
        ax3.set_xlabel('Motility Cost')
        ax3.set_ylabel('Ratio of Unchanged Assignments')
        ax3.set_ylim(0, 1)
        plt.tight_layout()
        plt.show()

    def normal_run(self):
        pareto_plot_series = []
        colors = sns.color_palette(n_colors=4)
        for norm, color in zip(["l1", "l2", "cosine", "l2+cosine"], colors):
            print(f"Computing cost matrices with norm: {norm}")
            xyz_cost_mat, exp_cost_mat = self.compute_cost_matrices(norm=norm)
            alpha, pareto_xyz_cost_list, pareto_exp_cost_list, lineage_xyz_cost, lineage_exp_cost = self.compute_pareto_front(xyz_cost_mat, exp_cost_mat, monte_carlo=False)
            self.plot_pareto_front(lineage_xyz_cost, lineage_exp_cost, pareto_xyz_cost_list, pareto_exp_cost_list, alpha, norm=norm, rescale=False, color=color)
            pareto_plot_series.append((lineage_xyz_cost, lineage_exp_cost, pareto_xyz_cost_list, pareto_exp_cost_list, alpha, norm))
            # print(f"Alpha: {alpha}")
        self.plot_pareto_front_series(pareto_plot_series)

    def norm_comp(self):
        l1_xyz_cost_mat, l1_exp_cost_mat = self.compute_cost_matrices(norm="l1", normalize=False)
        l2_xyz_cost_mat, l2_exp_cost_mat = self.compute_cost_matrices(norm="l2", normalize=False)
        cosine_xyz_cost_mat, cosine_exp_cost_mat = self.compute_cost_matrices(norm="cosine", normalize=False)

        # fix xyz cost as l2 norm, compute pareto fronts for l1, l2 and cosine expression norms
        # then as the alpha varies, check how the expression cost changes between different norms
        # for example, for l2+l2 norms, see how the expression cost changes for l1 and cosine norms
        for exp_cost_mat, norm in zip([l1_exp_cost_mat, l2_exp_cost_mat, cosine_exp_cost_mat], ["l1", "l2", "cosine"]):
            print(f"Computing cost matrices with expression norm: {norm}")
            step_size = 0.001
            loop_len = int(1 / step_size) + 1
            pareto_xyz_cost_list = []
            pareto_l1_exp_cost_list = []
            pareto_l2_exp_cost_list = []
            pareto_cosine_exp_cost_list = []
            for i in range(loop_len):
                alpha = 0 + step_size * i
                alpha = 0.5 * (1 - np.cos(np.pi * alpha)) 
                # alpha = 3 * alpha ** 3 - 2 * alpha ** 3
                cur_cost_mat = alpha * l2_xyz_cost_mat + (1 - alpha) * exp_cost_mat
                cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
                cur_xyz_cost = l2_xyz_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_l1_exp_cost = l1_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_l2_exp_cost = l2_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_cosine_exp_cost = cosine_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                pareto_xyz_cost_list.append(cur_xyz_cost)
                pareto_l1_exp_cost_list.append(cur_l1_exp_cost)
                pareto_l2_exp_cost_list.append(cur_l2_exp_cost)
                pareto_cosine_exp_cost_list.append(cur_cosine_exp_cost)
            # plotting l1 expression cost vs l2 expression cost vs cosine expression cost
            # plot as 3 separate 2d plots
            if norm == "l1":
                cur_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_name = "l2"
                comp2_name = "cosine"
            elif norm == "l2":
                cur_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_name = "l1"
                comp2_name = "cosine"
            elif norm == "cosine":
                cur_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp1_name = "l1"
                comp2_name = "l2"

            # Create combined dual-axis joint plot
            self.plot_dual_axis_jointplot(
                x_data=cur_pareto_exp_cost_list,
                y1_data=comp1_pareto_exp_cost_list,
                y2_data=comp2_pareto_exp_cost_list,
                x_label=f"{norm} norm",
                y1_label=f"{comp1_name} norm",
                y2_label=f"{comp2_name} norm",
                x_name=norm,
                y1_name=comp1_name,
                y2_name=comp2_name,
                title=f"{norm} vs {comp1_name} and {comp2_name} Expression Costs (Linear Alpha Sampling)"
            )
            print(f"unique number of {norm} cost values: {len(np.unique(cur_pareto_exp_cost_list))}")
        return pareto_l2_exp_cost_list, pareto_cosine_exp_cost_list
    

    def norm_comp_with_adaptive_sampling(self):
        l1_xyz_cost_mat, l1_exp_cost_mat = self.compute_cost_matrices(norm="l1", normalize=False)
        l2_xyz_cost_mat, l2_exp_cost_mat = self.compute_cost_matrices(norm="l2", normalize=False)
        cosine_xyz_cost_mat, cosine_exp_cost_mat = self.compute_cost_matrices(norm="cosine", normalize=False)
        round_digits = 7
        # results dict use hash(round(alpha, round_digits)) as key
        # alpha ranges from 0 to 1 with an adaptive samping strategy
        # the main sampling target is target expression cost, and we want a granularity of smaller than 0.5% 
        # of the total target expression cost range
        for exp_cost_mat, norm in zip([l1_exp_cost_mat, l2_exp_cost_mat, cosine_exp_cost_mat], ["L1", "L2", "Cosine"]):
            print(f"Computing cost matrices with expression norm: {norm}")
            results_dict = {}
            pareto_exp_cost_list = []
            pareto_l1_exp_cost_list = []
            pareto_l2_exp_cost_list = []
            pareto_cosine_exp_cost_list = []
            cur_row_ind, cur_col_ind = linear_sum_assignment(exp_cost_mat)
            min_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_l1_exp_cost = l1_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_l2_exp_cost = l2_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_cosine_exp_cost = cosine_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            results_dict[hash(round(0.0, round_digits))]  = (0.0, min_cost, cur_l1_exp_cost, cur_l2_exp_cost, cur_cosine_exp_cost)
            cur_row_ind, cur_col_ind = linear_sum_assignment(l2_xyz_cost_mat)
            max_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_l1_exp_cost = l1_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_l2_exp_cost = l2_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            cur_cosine_exp_cost = cosine_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
            results_dict[hash(round(1.0, round_digits))]  = (1.0, max_cost, cur_l1_exp_cost, cur_l2_exp_cost, cur_cosine_exp_cost)
            target_granularity = 0.005 * (max_cost - min_cost)
            sampling_stack = [(0.0, 1.0)]
            while sampling_stack:
                start_alpha, end_alpha = sampling_stack.pop()
                mid_alpha = (start_alpha + end_alpha) / 2
                if hash(round(mid_alpha, round_digits)) in results_dict:
                    continue
                cur_cost_mat = mid_alpha * l2_xyz_cost_mat + (1 - mid_alpha) * exp_cost_mat
                cur_row_ind, cur_col_ind = linear_sum_assignment(cur_cost_mat)
                cur_xyz_cost = exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_l1_exp_cost = l1_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_l2_exp_cost = l2_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                cur_cosine_exp_cost = cosine_exp_cost_mat[cur_row_ind, cur_col_ind].sum()
                results_dict[hash(round(mid_alpha, round_digits))]  = (mid_alpha, cur_xyz_cost, cur_l1_exp_cost, cur_l2_exp_cost, cur_cosine_exp_cost)
                # check granularity
                start_exp_cost = results_dict[hash(round(start_alpha, round_digits))][1]
                end_exp_cost = results_dict[hash(round(end_alpha, round_digits))][1]
                if abs(end_exp_cost - start_exp_cost) > target_granularity:
                    sampling_stack.append((start_alpha, mid_alpha))
                    sampling_stack.append((mid_alpha, end_alpha))
            results_list = list(results_dict.values())
            heapify(results_list)
            recorded_values = set()
            while results_list:
                entry = heappop(results_list)
                if round(entry[1], round_digits) in recorded_values:
                    continue
                recorded_values.add(round(entry[1], round_digits))
                pareto_exp_cost_list.append(entry[1])
                pareto_l1_exp_cost_list.append(entry[2])
                pareto_l2_exp_cost_list.append(entry[3])
                pareto_cosine_exp_cost_list.append(entry[4])
            # plotting l1 expression cost vs l2 expression cost vs cosine expression cost
            # plot as 3 separate 2d plots
            if norm == "L1":
                cur_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_name = "L2"
                comp2_name = "Cosine"
            elif norm == "L2":
                cur_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_name = "L1"
                comp2_name = "Cosine"
            elif norm == "Cosine":
                cur_pareto_exp_cost_list = pareto_cosine_exp_cost_list
                comp1_pareto_exp_cost_list = pareto_l1_exp_cost_list
                comp2_pareto_exp_cost_list = pareto_l2_exp_cost_list
                comp1_name = "L1"
                comp2_name = "L2"

            # Create combined dual-axis joint plot
            self.plot_dual_axis_jointplot(
                x_data=cur_pareto_exp_cost_list,
                y1_data=comp1_pareto_exp_cost_list,
                y2_data=comp2_pareto_exp_cost_list,
                x_label=f"{norm} norm",
                y1_label=f"{comp1_name} norm",
                y2_label=f"{comp2_name} norm",
                x_name=norm,
                y1_name=comp1_name,
                y2_name=comp2_name,
                title=f"{norm} vs {comp1_name} and {comp2_name} Exp Costs (Adaptive Sampling)"
            )            
            print(f"unique number of {norm} cost values: {len(recorded_values)}")
        
        return results_dict

        