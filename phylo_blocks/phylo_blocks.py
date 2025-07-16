import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from itertools import combinations
from ete3 import Tree, TreeStyle, NodeStyle
import os

from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from collections import defaultdict, Counter

# code for pre-computing lca cache
def build_descendant_map(tree):
    """Maps each node name to the set of its descendant leaves."""
    descendant_map = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            descendant_map[node.name] = {node.name}
        else:
            descendants = set()
            for child in node.children:
                descendants.update(descendant_map[child.name])
            descendant_map[node.name] = descendants
    return descendant_map

def compute_lca_cache_fast_v2(tree):
    '''
    code for pre-computing an lca cache for distance calculations
    '''
    
    lca_cache = {}
    seen = set()
    descendant_map = build_descendant_map(tree)
    leaves = tree.get_leaf_names()
    
    for leaf_num, leaf in enumerate(tree.iter_leaves()):
        current = leaf
        path = []
        while current:
            path.append(current)
            current = current.up

        for i, node in enumerate(path):
            new_descs = descendant_map[node.name] - descendant_map[path[i-1].name]

            for nd in new_descs:
                # add, sorting first
                key = (leaf.name, nd) if leaf.name < nd else (nd, leaf.name)
                if key not in lca_cache:
                    lca_cache[key] = node.name
            # seen.add(leaf.name)

        # Print progress only every 1% or every 10,000 steps
        if leaf_num % max(1, leaf_num // 1000) == 0:
            frac = leaf_num / len(leaves)
            bar = int(frac * 40)
            print(f'{frac:<20.6f} |' + '*' * bar + '-' * (40 - bar) + '|', end='\r')

    print()
    
    return lca_cache

def fast_cophenetic_matrix(tree: Tree, d_root : dict, lca_cache):
    leaves = tree.get_leaf_names()
    leaf_index = {name: i for i, name in enumerate(leaves)}
    n = len(leaves)

    matrix = np.zeros((n, n))

    for i, a in enumerate(leaves):
        for j in range(i + 1, n):
            b = leaves[j]
            lca = lca_cache[(a, b) if a < b else (b, a)]
            d = d_root[a] + d_root[b] - 2 * d_root[lca]
            matrix[i, j] = matrix[j, i] = d

        frac = (i * j) / (len(leaves) ** 2)
        count = frac
        frac = int(frac * 40)
        
        print( f'{count:<20}' + '|' + '*' * frac + '-' * (40-frac) + '|' , end='\r')

    return matrix, leaves

def cluster_cophenetic_matrix(D, labels, k=5):
    """
    D: patristic distance matrix (n x n)
    labels: list of leaf names
    k: number of folds
    method: 'mds'
    """
    print('performing embeddings')
    embed = MDS(n_components=5, dissimilarity='precomputed', random_state=12345, n_jobs=16)
    coords = embed.fit_transform(D)

    # KMeans on embedding
    print('performing kmeans')
    kmeans = KMeans(n_clusters=k, random_state=111, n_init='auto')
    cluster_labels = kmeans.fit_predict(coords)

    # Build fold_dict
    print('Performing dictionary')
    fold_dict = {name: int(fold) for name, fold in zip(labels, cluster_labels)}

    # Optional sanity check
    print("Cluster sizes:", Counter(fold_dict.values()))
    return fold_dict

# Function to color branches by fold
def color_branches_by_fold(tree, fold_dict):
    
    colors = ["red", "blue", "green", "orange", "purple"]

    for leaf in tree.iter_leaves():
        fold = fold_dict.get(leaf.name, None)
        if fold is None:
            continue  # skip leaves not in fold_dict

        color = colors[fold % len(colors)]

        # Color path from leaf to root
        node = leaf
        while node.up:
            ns = NodeStyle()
            ns["vt_line_color"] = color   # branch color
            ns["hz_line_color"] = color
            ns["size"] = 0          # hide node circles
            node.set_style(ns)
            node = node.up

def compute_min_interfold_distances(D, leaves, fold_dict):
    """
    D: n x n cophenetic distance matrix
    labels: list of leaf names corresponding to D rows/columns
    fold_dict: mapping from leaf name to fold ID
    Returns: 2D array [k x k] where entry (i, j) is min dist between fold i and j
    """
    n = len(leaves)
    folds = sorted(set(fold_dict.values()))
    k = len(folds)

    # Map fold ID to list of indices in D
    fold_indices = defaultdict(list)
    for i, leaf_name in enumerate(leaves):
        fold = fold_dict[leaf_name]
        fold_indices[fold].append(i)

    # Initialize result matrix
    min_dist = np.full((k, k), np.inf)

    for i in folds:
        for j in folds:
            if i == j:
                continue
            for a in fold_indices[i]:
                for b in fold_indices[j]:
                    d = D[a, b] if a < b else D[b, a]
                    if d < min_dist[i, j]:
                        min_dist[i, j] = d

    return min_dist


