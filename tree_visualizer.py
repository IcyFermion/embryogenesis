class TreeVisualizer:
    def __init__(self, compact=True, max_depth=None, show_coordinates=True):
        """
        Initialize the tree visualizer.
        
        Args:
            compact: If True, uses more compact notation
            max_depth: Maximum depth to display (None for all levels)
            show_coordinates: Whether to show node coordinates
        """
        self.compact = compact
        self.max_depth = max_depth
        self.show_coordinates = show_coordinates
        self.stats = {
            'total_nodes': 0,
            'max_depth': 0,
            'nodes_per_level': {}
        }
    
    def visualize(self, root, collapse_after=None):
        """
        Generate a text visualization of the tree.
        
        Args:
            root: The root node of the tree
            collapse_after: Depth after which to show summary only
        """
        if not root:
            return "Empty tree"
        
        # First, gather statistics
        self._gather_stats(root)
        
        # Generate the visualization
        lines = []
        lines.append("=" * 60)
        lines.append(f"Tree Visualization (Total nodes: {self.stats['total_nodes']}, Depth: {self.stats['max_depth'] + 1})")
        lines.append("=" * 60)
        lines.append("")
        
        # Add legend
        lines.append("Legend: [O] = Optimized, [T] = Terminal, [*] = Regular node")
        lines.append("-" * 60)
        lines.append("")
        
        # Visualize the tree
        self._visualize_node(root, lines, 0, collapse_after)
        
        # Add summary statistics
        lines.append("")
        lines.append("-" * 60)
        lines.append("Level Statistics:")
        for level in sorted(self.stats['nodes_per_level'].keys()):
            count = self.stats['nodes_per_level'][level]
            lines.append(f"  Level {level:2d}: {count:3d} nodes")
        
        return "\n".join(lines)
    
    def _gather_stats(self, node, depth=0):
        """Gather statistics about the tree."""
        self.stats['total_nodes'] += 1
        self.stats['max_depth'] = max(self.stats['max_depth'], depth)
        
        if depth not in self.stats['nodes_per_level']:
            self.stats['nodes_per_level'][depth] = 0
        self.stats['nodes_per_level'][depth] += 1
        
        for child in node.children:
            self._gather_stats(child, depth + 1)
    
    def _visualize_node(self, node, lines, depth, collapse_after, is_last_child=True, prefix=""):
        """Recursively visualize a node and its children."""
        if self.max_depth is not None and depth > self.max_depth:
            return
        
        # Create the node representation
        if self.compact:
            indent = prefix
            if depth > 0:
                connector = "└── " if is_last_child else "├── "
                indent += connector
        else:
            indent = "  " * depth + "├─ "
            if depth > 0 and is_last_child:
                indent = "  " * depth + "└─ "
        
        # Node symbol
        if node.optimized and node.is_terminal:
            symbol = "[O,T]"
        elif node.optimized:
            symbol = "[O]"
        elif node.is_terminal:
            symbol = "[T]"
        else:
            symbol = "[*]"
        
        # Build node info
        node_info = f"{symbol} {node.key}"
        if self.show_coordinates and node.coordinate:
            node_info += f" @ {node.coordinate}"
        
        # Check if we should collapse
        if collapse_after is not None and depth >= collapse_after:
            child_count = self._count_descendants(node)
            if child_count > 0:
                node_info += f" ... ({child_count} descendants)"
            lines.append(indent + node_info)
            return
        
        lines.append(indent + node_info)
        
        # Process children
        child_count = len(node.children)
        for i, child in enumerate(node.children):
            is_last = (i == child_count - 1)
            
            if self.compact and depth > 0:
                # Update prefix for children
                if is_last_child:
                    child_prefix = prefix + "    "
                else:
                    child_prefix = prefix + "│   "
            else:
                child_prefix = prefix
            
            self._visualize_node(child, lines, depth + 1, collapse_after, is_last, child_prefix)
    
    def _count_descendants(self, node):
        """Count all descendants of a node."""
        count = len(node.children)
        for child in node.children:
            count += self._count_descendants(child)
        return count
    
    def visualize_breadth_first(self, root, show_level_separators=True):
        """
        Alternative visualization showing nodes level by level.
        Good for seeing patterns across levels.
        """
        if not root:
            return "Empty tree"
        
        from collections import deque
        
        lines = []
        lines.append("=" * 80)
        lines.append("Breadth-First Tree Visualization")
        lines.append("=" * 80)
        lines.append("")
        
        queue = deque([(root, 0)])
        current_level = -1
        level_nodes = []
        
        while queue:
            node, level = queue.popleft()
            
            if level != current_level:
                # Output previous level
                if current_level >= 0:
                    self._output_level(lines, current_level, level_nodes, show_level_separators)
                
                current_level = level
                level_nodes = []
                
                # Stop if we've reached max depth
                if self.max_depth is not None and level > self.max_depth:
                    lines.append(f"\n... (Truncated at depth {self.max_depth})")
                    break
            
            level_nodes.append(node)
            
            # Add children to queue
            for child in node.children:
                queue.append((child, level + 1))
        
        # Output last level
        if level_nodes and (self.max_depth is None or current_level <= self.max_depth):
            self._output_level(lines, current_level, level_nodes, show_level_separators)
        
        return "\n".join(lines)
    
    def _output_level(self, lines, level, nodes, show_separator):
        """Output all nodes at a given level."""
        if show_separator:
            lines.append(f"\n{'=' * 20} Level {level} ({len(nodes)} nodes) {'=' * 20}")
        else:
            lines.append(f"\nLevel {level} ({len(nodes)} nodes):")
        
        # Group nodes by parent for better readability
        node_strs = []
        for node in nodes:
            symbol = "[O]" if node.optimized else "[T]" if node.is_terminal else "[ ]"
            node_str = f"{symbol} {node.key}"
            if self.show_coordinates and node.coordinate:
                node_str += f"@{node.coordinate}"
            node_strs.append(node_str)
        
        # Output nodes in rows
        nodes_per_row = 10 if self.compact else 5
        for i in range(0, len(node_strs), nodes_per_row):
            row_nodes = node_strs[i:i + nodes_per_row]
            lines.append("  " + "  ".join(row_nodes))


# Example usage and testing
def create_test_tree(depth=5, children_per_node=3):
    """Create a test tree for demonstration."""
    class Node:
        def __init__(self, key, coordinate):
            self.key = key
            self.coordinate = coordinate
            self.children = []
            self.optimized = False
            self.parent = None
            self.is_terminal = False
    
    def build_tree(key_prefix, current_depth, max_depth):
        node = Node(key_prefix, (current_depth, len(key_prefix)))
        
        # Randomly set some properties
        import random
        node.optimized = random.random() > 0.7
        node.is_terminal = current_depth == max_depth or random.random() > 0.9
        
        if current_depth < max_depth and not node.is_terminal:
            num_children = random.randint(2, children_per_node)
            for i in range(num_children):
                child_key = f"{key_prefix}.{i}"
                child = build_tree(child_key, current_depth + 1, max_depth)
                child.parent = node
                node.children.append(child)
        
        return node
    
    return build_tree("root", 0, depth)


# Example of how to use the visualizer
if __name__ == "__main__":
    # Create a test tree
    root = create_test_tree(depth=6, children_per_node=4)
    
    # Create visualizer
    viz = TreeVisualizer(compact=True, show_coordinates=True)
    
    print("1. Standard Tree Visualization (Collapsed after level 3):")
    print(viz.visualize(root, collapse_after=3))
    
    print("\n\n2. Breadth-First Visualization (First 4 levels):")
    viz.max_depth = 3
    print(viz.visualize_breadth_first(root))
    
    # For very deep trees (like your 70-level tree), you might want:
    print("\n\n3. Ultra-Compact Mode for Deep Trees:")
    viz_deep = TreeVisualizer(compact=True, show_coordinates=False, max_depth=10)
    print(viz_deep.visualize(root, collapse_after=5))