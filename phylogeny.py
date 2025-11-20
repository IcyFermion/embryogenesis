# %%
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Dict

@dataclass
class Node:
    """A node in the phylogenetic tree"""
    name: str
    branch_length: float = 0.0  # length from parent to this node
    parent: Optional['Node'] = None
    children: List['Node'] = field(default_factory=list)
    trait_value: Optional[float] = None  # for storing trait data
    time_from_root: float = 0.0  # cumulative time from root
    
    def add_child(self, child: 'Node'):
        """Add a child node"""
        child.parent = self
        self.children.append(child)
        child.time_from_root = self.time_from_root + child.branch_length
        
    def is_tip(self) -> bool:
        """Check if this is a tip (leaf) node"""
        return len(self.children) == 0
    
    def get_tips(self) -> List['Node']:
        """Recursively get all tips descended from this node"""
        if self.is_tip():
            return [self]
        tips = []
        for child in self.children:
            tips.extend(child.get_tips())
        return tips

class PhyloTree:
    """Simple phylogenetic tree implementation"""
    
    def __init__(self, root: Node):
        self.root = root
        self.tips = []
        self.internal_nodes = []
        self._catalog_nodes()
    
    def _catalog_nodes(self):
        """Catalog all tips and internal nodes"""
        self._recursive_catalog(self.root)
        
    def _recursive_catalog(self, node: Node):
        """Recursively catalog nodes"""
        if node.is_tip():
            self.tips.append(node)
        else:
            self.internal_nodes.append(node)
            for child in node.children:
                self._recursive_catalog(child)
    
    def compute_vcv_matrix(self, sigma_sq: float = 1.0) -> np.ndarray:
        """
        Compute variance-covariance matrix under Brownian motion
        
        This is the KEY method - let's understand what it does!
        """
        n_tips = len(self.tips)
        vcv = np.zeros((n_tips, n_tips))
        
        # For each pair of tips, find their shared evolutionary history
        for i, tip_i in enumerate(self.tips):
            for j, tip_j in enumerate(self.tips):
                if i <= j:  # matrix is symmetric
                    shared_time = self._shared_evolutionary_time(tip_i, tip_j)
                    vcv[i, j] = shared_time * sigma_sq
                    vcv[j, i] = shared_time * sigma_sq
        
        return vcv
    
    def _shared_evolutionary_time(self, node1: Node, node2: Node) -> float:
        """
        Calculate shared evolutionary time from root to MRCA of two nodes
        
        This is where the magic happens!
        """
        # Get paths from root to each node
        path1 = self._path_to_root(node1)
        path2 = self._path_to_root(node2)
        
        # Find MRCA (most recent common ancestor)
        mrca = None
        for n1 in path1:
            if n1 in path2:
                mrca = n1
                break
        
        # Shared time is time from root to MRCA
        return mrca.time_from_root if mrca else 0.0
    
    def _path_to_root(self, node: Node) -> List[Node]:
        """Get path from node to root"""
        path = []
        current = node
        while current is not None:
            path.append(current)
            current = current.parent
        return path
    
def create_example_tree():
    """Create a simple example tree with 30 tips"""
    # Create root
    root = Node("root")
    
    # First major split
    clade1 = Node("clade1", branch_length=5.0)
    clade2 = Node("clade2", branch_length=5.0)
    root.add_child(clade1)
    root.add_child(clade2)
    
    # Add more structure (simplified for clarity)
    # In practice, you'd have more complex branching
    # ... (we'd add the rest of the 30 tips here)
    
    return PhyloTree(root)


def create_balanced_cell_lineage_tree(n_generations=5):
    """
    Create a balanced binary tree representing cell divisions.
    n_generations=5 gives us 32 tips (2^5), close to our target of 30.
    
    For cell lineages, branch lengths might represent:
    - Time between divisions
    - Number of cell cycles
    - Developmental time
    """
    
    def build_subtree(parent_node, depth, max_depth, cell_id_counter):
        """Recursively build a balanced binary tree"""
        if depth == max_depth:
            return
        
        # Create two daughter cells
        for i in range(2):
            cell_id_counter[0] += 1
            
            if depth == max_depth - 1:
                # This will be a tip (terminal cell)
                child_name = f"cell_{cell_id_counter[0]}"
            else:
                # This is an intermediate division
                child_name = f"div_{cell_id_counter[0]}"
            
            # Branch length could represent time to next division
            # Let's make it slightly variable to be more realistic
            branch_length = 1.0 + np.random.normal(0, 0.1)  # ~1 time unit ± noise
            branch_length = max(0.1, branch_length)  # ensure positive
            
            child = Node(child_name, branch_length=branch_length)
            parent_node.add_child(child)
            
            # Recursively add children
            build_subtree(child, depth + 1, max_depth, cell_id_counter)
    
    # Create root (zygote or initial cell)
    root = Node("initial_cell", branch_length=0.0)
    cell_id_counter = [0]  # mutable counter
    
    # Build the tree
    build_subtree(root, 0, n_generations, cell_id_counter)
    
    # Update time_from_root for all nodes
    def update_times(node, time=0):
        node.time_from_root = time
        for child in node.children:
            update_times(child, time + child.branch_length)
    
    update_times(root)
    
    return PhyloTree(root)

# Create our cell lineage tree
tree = create_balanced_cell_lineage_tree(n_generations=5)
print(f"Created tree with {len(tree.tips)} cells")
print(f"Tree has {len(tree.internal_nodes)} internal nodes (cell divisions)")

# Let's simulate gene expression data for these cells
np.random.seed(42)  # for reproducibility

def simulate_gene_expression_on_tree(tree, ancestral_expression=10, sigma_sq=0.5):
    """
    Simulate gene expression evolution along the cell lineage tree.
    
    ancestral_expression: Expression level in the initial cell
    sigma_sq: Rate of expression change (variance per unit time)
    """
    
    def simulate_brownian_motion(node, parent_value):
        """Recursively simulate expression values"""
        if node.parent is None:
            # Root node
            node.trait_value = ancestral_expression
        else:
            # Expression change is normal with variance = branch_length * sigma_sq
            expression_change = np.random.normal(
                0, 
                np.sqrt(node.branch_length * sigma_sq)
            )
            node.trait_value = parent_value + expression_change
        
        # Simulate for children
        for child in node.children:
            simulate_brownian_motion(child, node.trait_value)
    
    # Start simulation from root
    simulate_brownian_motion(tree.root, None)
    
    # Extract tip values
    tip_values = np.array([tip.trait_value for tip in tree.tips])
    return tip_values

# Simulate expression data
expression_values = simulate_gene_expression_on_tree(tree)
print(f"\nExpression values range: {expression_values.min():.2f} to {expression_values.max():.2f}")
print(f"Mean expression: {expression_values.mean():.2f}")
# %%

vcv_matrix = tree.compute_vcv_matrix(sigma_sq=1.0)
print(f"\nVCV matrix shape: {vcv_matrix.shape}")


# %%
def reconstruct_all_ancestral_states(tree, tip_values, sigma_sq=0.5):
    """
    Reconstruct gene expression at ALL internal nodes in the tree.
    
    This uses the "local" or "three-point" algorithm which works 
    node by node from tips toward root, then back down.
    """
    
    # First, get the full VCV matrix for tips
    vcv = tree.compute_vcv_matrix(sigma_sq)
    vcv_inv = np.linalg.inv(vcv)
    
    # Store reconstructed values
    reconstructed_values = {}
    reconstruction_variance = {}
    
    # Step 1: Calculate contrasts and work up from tips to root
    # This is based on Felsenstein's independent contrasts algorithm
    
    def calculate_node_value(node, tip_values_dict):
        """
        Recursively calculate ancestral states from tips up to root.
        Uses weighted averaging based on branch lengths.
        """
        if node.is_tip():
            # For tips, we have observed values
            return node.trait_value, 0  # variance is 0 for observed values
        
        # For internal nodes, calculate weighted average of children
        child_values = []
        child_weights = []
        child_variances = []
        
        for child in node.children:
            if child.is_tip():
                val = tip_values_dict[child.name]
                var = 0
            else:
                val, var = calculate_node_value(child, tip_values_dict)
                reconstructed_values[child.name] = val
                reconstruction_variance[child.name] = var
            
            # Weight is inverse of branch length plus accumulated variance
            weight = 1.0 / (child.branch_length * sigma_sq + var)
            child_values.append(val)
            child_weights.append(weight)
            child_variances.append(var)
        
        # Weighted average
        total_weight = sum(child_weights)
        weighted_value = sum(v * w for v, w in zip(child_values, child_weights)) / total_weight
        
        # Variance of reconstruction
        node_variance = 1.0 / total_weight
        
        return weighted_value, node_variance
    
    # Create tip values dictionary
    tip_values_dict = {tip.name: val for tip, val in zip(tree.tips, tip_values)}
    
    # Reconstruct from tips to root
    root_value, root_var = calculate_node_value(tree.root, tip_values_dict)
    reconstructed_values[tree.root.name] = root_value
    reconstruction_variance[tree.root.name] = root_var
    
    return reconstructed_values, reconstruction_variance

# Now let's do a more sophisticated version using the full algorithm
def reconstruct_all_nodes_full_algorithm(tree: PhyloTree, tip_values, sigma_sq=0.5):
    """
    Full ancestral state reconstruction using the two-pass algorithm:
    1. Pass up: Calculate partial likelihoods from tips to root
    2. Pass down: Calculate final states from root to internal nodes
    
    This gives us both the ML estimates AND their uncertainty!
    """
    
    # Initialize storage for all nodes
    all_nodes = [tree.root] + tree.internal_nodes + tree.tips
    node_states = {}
    node_variances = {}
    
    # For tips, set observed values
    for tip, value in zip(tree.tips, tip_values):
        node_states[tip.name] = value
        node_variances[tip.name] = 0  # No uncertainty in observed values
    
    # Pass 1: Calculate preliminary estimates going up the tree
    def upward_pass(node):
        """Calculate preliminary estimates from children"""
        if node.is_tip():
            return node_states[node.name], 0
        
        child_estimates = []
        for child in node.children:
            if child.name not in node_states:
                child_val, child_var = upward_pass(child)
                node_states[child.name] = child_val
                node_variances[child.name] = child_var
            
            child_estimates.append({
                'value': node_states[child.name],
                'variance': node_variances[child.name],
                'branch_length': child.branch_length
            })
        
        # Calculate weighted average for this node
        weights = [1/(c['branch_length'] * sigma_sq + c['variance']) 
                  for c in child_estimates]
        values = [c['value'] for c in child_estimates]
        
        total_weight = sum(weights)
        weighted_mean = sum(v*w for v,w in zip(values, weights)) / total_weight
        variance = 1 / total_weight
        
        return weighted_mean, variance
    
    # Do upward pass
    root_val, root_var = upward_pass(tree.root)
    node_states[tree.root.name] = root_val
    node_variances[tree.root.name] = root_var
    
    # Pass 2: Refine estimates going down the tree (optional but more accurate)
    # This accounts for information from the whole tree, not just descendants
    
    def downward_pass(node, parent_value=None):
        """Refine estimates using parent information"""
        if node.parent is None:
            # Root node - already have best estimate
            return
        
        if node.is_tip():
            # Tips don't need refinement
            return
        
        # Use parent information to refine estimate
        # This is a weighted average between upward estimate and parent-based prediction
        upward_est = node_states[node.name]
        upward_var = node_variances[node.name]
        
        # Prediction from parent
        parent_pred = parent_value  # Under BM, expected value = parent value
        parent_pred_var = node.branch_length * sigma_sq
        
        # Combine estimates (inverse-variance weighting)
        w_upward = 1 / upward_var if upward_var > 0 else float('inf')
        w_parent = 1 / parent_pred_var
        
        if not np.isinf(w_upward):
            refined_value = (upward_est * w_upward + parent_pred * w_parent) / (w_upward + w_parent)
            refined_var = 1 / (w_upward + w_parent)
        else:
            refined_value = upward_est
            refined_var = upward_var
        
        node_states[node.name] = refined_value
        node_variances[node.name] = refined_var
        
        # Recursively refine children
        for child in node.children:
            downward_pass(child, refined_value)
    
    # Do downward pass
    for child in tree.root.children:
        downward_pass(child, root_val)
    
    return node_states, node_variances

# Let's test this on our tree!
reconstructed_states, variances = reconstruct_all_nodes_full_algorithm(
    tree, expression_values, sigma_sq=0.5
)

# Compare true vs reconstructed for internal nodes
print("\nComparing true vs reconstructed expression at internal nodes:")
print("-" * 60)
for node in tree.internal_nodes[:5]:  # Show first 5 for brevity
    true_val = node.trait_value
    recon_val = reconstructed_states[node.name]
    uncertainty = np.sqrt(variances[node.name])
    print(f"{node.name:10s}: True={true_val}, Reconstructed={recon_val} ± {uncertainty}")
# %%
def estimate_sigma_squared(tree: PhyloTree, tip_values):
    """
    Estimate sigma^2 (rate of evolution) using maximum likelihood.
    
    This finds the σ² that makes the observed tip values most likely
    given the tree structure and Brownian motion model.
    """
    
    # Method 1: Analytical MLE (most efficient)
    def analytical_mle():
        """
        Closed-form solution for σ² MLE
        """
        n = len(tip_values)
        
        # Build VCV matrix with σ² = 1 (we'll scale it)
        vcv_unit = tree.compute_vcv_matrix(sigma_sq=1.0)
        
        # Get inverse and determinant
        vcv_inv = np.linalg.inv(vcv_unit)
        
        # Estimate ancestral state at root first
        ones = np.ones(n)
        root_state = (ones @ vcv_inv @ tip_values) / (ones @ vcv_inv @ ones)
        
        # Center the data (subtract root state)
        centered_values = tip_values - root_state
        
        # MLE for σ²
        sigma_sq_mle = (centered_values @ vcv_inv @ centered_values) / n
        
        return sigma_sq_mle, root_state
    
    # Method 2: Using likelihood optimization (more flexible, can add constraints)
    def likelihood_optimization():
        """
        Find σ² by maximizing the likelihood function
        """
        from scipy.optimize import minimize_scalar
        
        def negative_log_likelihood(sigma_sq):
            """
            Calculate negative log-likelihood for given σ²
            """
            if sigma_sq <= 0:
                return np.inf
            
            vcv = tree.compute_vcv_matrix(sigma_sq)
            
            try:
                # Multivariate normal log-likelihood
                vcv_inv = np.linalg.inv(vcv)
                vcv_det = np.linalg.det(vcv)
                
                # Estimate root state for this σ²
                ones = np.ones(len(tip_values))
                root_state = (ones @ vcv_inv @ tip_values) / (ones @ vcv_inv @ ones)
                
                # Calculate likelihood
                diff = tip_values - root_state
                log_lik = -0.5 * (len(tip_values) * np.log(2 * np.pi) + 
                                  np.log(vcv_det) + 
                                  diff @ vcv_inv @ diff)
                
                return -log_lik  # Return negative for minimization
                
            except np.linalg.LinAlgError:
                return np.inf
        
        # Optimize
        result = minimize_scalar(negative_log_likelihood, 
                                bounds=(0.001, 100), 
                                method='bounded')
        
        sigma_sq_mle = result.x
        
        # Get root state at optimal σ²
        vcv = tree.compute_vcv_matrix(sigma_sq_mle)
        vcv_inv = np.linalg.inv(vcv)
        ones = np.ones(len(tip_values))
        root_state = (ones @ vcv_inv @ tip_values) / (ones @ vcv_inv @ ones)
        
        return sigma_sq_mle, root_state
    
    # Method 3: REML (Restricted Maximum Likelihood) - accounts for uncertainty in root
    def reml_estimate():
        """
        REML estimate - often preferred as it's less biased
        """
        n = len(tip_values)
        vcv_unit = tree.compute_vcv_matrix(sigma_sq=1.0)
        vcv_inv = np.linalg.inv(vcv_unit)
        
        ones = np.ones(n)
        root_state = (ones @ vcv_inv @ tip_values) / (ones @ vcv_inv @ ones)
        centered_values = tip_values - root_state
        
        # REML correction: divide by (n-1) instead of n
        sigma_sq_reml = (centered_values @ vcv_inv @ centered_values) / (n - 1)
        
        return sigma_sq_reml, root_state
    
    # Let's use the analytical MLE
    sigma_sq_est, root_est = analytical_mle()
    
    # Also calculate confidence interval using likelihood profile
    def get_confidence_interval(mle_value):
        """
        Calculate 95% CI using likelihood ratio test
        """
        from scipy.stats import chi2
        
        vcv_unit = tree.compute_vcv_matrix(sigma_sq=1.0)
        vcv_inv = np.linalg.inv(vcv_unit)
        n = len(tip_values)
        
        # Log-likelihood at MLE
        vcv_mle = vcv_unit * mle_value
        ones = np.ones(n)
        root_state = (ones @ vcv_inv @ tip_values) / (ones @ vcv_inv @ ones)
        diff = tip_values - root_state
        
        log_lik_mle = -0.5 * (n * np.log(2 * np.pi * mle_value) + 
                              np.log(np.linalg.det(vcv_unit)) + 
                              diff @ vcv_inv @ diff / mle_value)
        
        # Find values where log-likelihood drops by chi2(1, 0.95)/2
        threshold = log_lik_mle - chi2.ppf(0.95, df=1) / 2
        
        # Simple search for CI bounds (could be more sophisticated)
        lower = mle_value * 0.5
        upper = mle_value * 2.0
        
        return lower, upper
    
    return sigma_sq_est, root_est

# Test on our simulated data
estimated_sigma_sq, estimated_root = estimate_sigma_squared(tree, expression_values)

print(f"True σ² used in simulation: 0.5")
print(f"Estimated σ² from data: {estimated_sigma_sq:.3f}")
print(f"True root expression: 10.0")
print(f"Estimated root expression: {estimated_root:.2f}")
# %%
