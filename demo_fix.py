#!/usr/bin/env python3
"""
Demonstration script showing that the PicklingError has been resolved.

This script demonstrates that the LineageOptimization class can now be used
with multiprocessing without encountering the PicklingError that occurred
when the class was defined in a Jupyter notebook.
"""

import numpy as np
from lineage_optimization import LineageTree, LineageOptimization
from tqdm.contrib.concurrent import process_map

def main():
    print("🧬 Embryogenesis LineageOptimization PicklingError Fix Demo")
    print("=" * 60)
    
    print("\n1. Creating test data...")
    # Create realistic test data
    np.random.seed(42)  # For reproducible results
    n_cells = 15
    xyz_mat = np.random.rand(n_cells, 3) * 10  # 3D coordinates
    exp_mat = np.random.rand(n_cells, 8)       # Gene expression data
    
    print(f"   Created data for {n_cells} cells")
    print(f"   Spatial coordinates: {xyz_mat.shape}")
    print(f"   Expression data: {exp_mat.shape}")
    
    print("\n2. Building lineage tree...")
    # Create a lineage tree structure
    lineage_tree = LineageTree()
    
    # Add root nodes
    lineage_tree.add_node(0, -1)  # Root
    lineage_tree.root = 0
    lineage_tree.add_node(1, 0)   # First branch
    lineage_tree.add_node(2, 0)   # Second branch
    
    # Add more nodes to create a realistic tree
    for i in range(3, n_cells):
        parent = np.random.choice(range(max(1, i-3), i))
        lineage_tree.add_node(i, parent)
    
    first_internal_layer = [(1, 1), (2, 1)]
    lineage_names = [f"Cell_{i:02d}" for i in range(n_cells)]
    
    print(f"   Created tree with {lineage_tree.size} nodes")
    
    print("\n3. Creating LineageOptimization instance...")
    opt = LineageOptimization(
        xyz_mat=xyz_mat,
        exp_mat=exp_mat, 
        lineage_tree=lineage_tree,
        first_internal_layer=first_internal_layer,
        lineage_names=lineage_names
    )
    
    print(f"   Initial lineage costs:")
    print(f"   - Spatial cost: {opt.lineage_xyz_cost:.3f}")
    print(f"   - Expression cost: {opt.lineage_exp_cost:.3f}")
    
    print("\n4. Testing multiprocessing capabilities...")
    print("   Running mst_test with multiprocessing...")
    
    try:
        # This would previously fail with PicklingError
        mst_results = process_map(
            opt.mst_test,
            range(20),  # Small test set
            max_workers=4,
            chunksize=5,
            desc="   MST Tests"
        )
        
        print(f"   ✅ SUCCESS: Completed {len(mst_results)} MST tests")
        print(f"   Average MST cost: {np.mean(mst_results):.3f}")
        print(f"   Cost range: {np.min(mst_results):.3f} - {np.max(mst_results):.3f}")
        
    except Exception as e:
        print(f"   ❌ FAILED: {e}")
        return False
    
    print("\n5. Testing optimization methods...")
    try:
        # Test one of the optimization methods that uses multiprocessing
        print("   Running small-scale optimization test...")
        
        # Create a smaller test for faster execution
        small_results = []
        for i in range(3):
            result = opt.bottom_up_by_layer(
                first_internal_tree_ids=[1, 2],
                internal_tree_ids=opt.internal_tree_ids[:5],  # Limit for speed
                tree_ids_by_depth=opt.tree_ids_by_depth,
                idx=i
            )
            small_results.append(result)
        
        print(f"   ✅ Optimization methods work correctly")
        print(f"   Sample optimization results: {len(small_results)} computed")
        
    except Exception as e:
        print(f"   ❌ Optimization test failed: {e}")
        return False
    
    print("\n🎉 SUCCESS: All tests passed!")
    print("\nThe PicklingError has been successfully resolved by:")
    print("1. Extracting LineageTree and LineageOptimization classes to lineage_optimization.py")
    print("2. Updating the Jupyter notebook to import classes instead of defining them inline")
    print("3. Classes defined in proper Python modules can be pickled for multiprocessing")
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        print("\n❌ Demo failed - issue may not be fully resolved")
        exit(1)
    else:
        print("\n✨ Demo completed successfully!")