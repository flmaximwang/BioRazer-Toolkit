from treelib import Tree
from pyrosetta.rosetta.core.kinematics import FoldTree


class FoldTreeVisualizer:

    def __init__(self, fold_tree: FoldTree):
        self.fold_tree = fold_tree
        self.visualizer = Tree()
        root = fold_tree.root()
        self.visualizer.create_node(f"{root}", root)
        edges_str = str(fold_tree).split("EDGE ")[1:]
        for edge_str in edges_str:
            edge_str_list = edge_str.split(" ")
            edge_parent = int(edge_str_list[0])
            edge_child = int(edge_str_list[1])
            edge_desc = " ".join(edge_str_list[2:])
            self.visualizer.create_node(
                f"{edge_child}: {edge_desc}", edge_child, parent=edge_parent
            )

    def show(self):
        self.visualizer.show()


def show_fold_tree(fold_tree: FoldTree):
    """
    Visualize the FoldTree structure using the Tree library.
    """
    visualizer = FoldTreeVisualizer(fold_tree)
    visualizer.show()
    return visualizer
