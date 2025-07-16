import argparse

from phylo_blocks import *

parser = argparse.ArgumentParser()
parser.add_argument('--tree', '-t', required=True, type=str)
parser.add_argument('--n-folds', '-n', required=True, type=int)
parser.add_argument('--output', '-o', required=True, type=str)

if __name__ == '__main__':
    args = parser.parse_args()

    tree = Tree(args.tree, format=1, quoted_node_names=True)

    # reformat the tree - unique-ify internal node names, calculate root dists
    # generating pre-processed root-node distance matrix
    # uniquely re-naming each node
    for i, node in enumerate(tree.traverse('levelorder')):
    
        if not node.is_leaf():
            node.name = f'c{i}'
            
    # compute distances to roots
    print('Computing root dists')
    root = tree.get_tree_root()
    d_root = {node.name: node.get_distance(root) for node in tree.traverse("preorder")}

    # pre-compute maps
    print('Pre-computing maps')
    desc_map = build_descendant_map(tree)
    lca_cache = compute_lca_cache_fast_v2(tree)

    print('get cophenetic matrix')
    matrix, leaves = fast_cophenetic_matrix(tree, d_root, lca_cache)

    print('Getting folds')
    fold_dict = cluster_cophenetic_matrix(matrix, leaves, args.n_folds)

    # Apply coloring
    color_branches_by_fold(tree, fold_dict)

    print(compute_min_interfold_distances(matrix, leaves, fold_dict))
    
    # Create and show tree style
    print('Drawing')
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
    ts = TreeStyle()
    ts.show_leaf_name = False
    tree.render(os.path.join(args.output, 'tree.png'), tree_style=ts, w=400)
    


    