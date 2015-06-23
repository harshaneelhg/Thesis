# This piece of code finds out minimum spanning tree for 
# Disease-Chemical-Gene association graph using 'Networkx'
# Python library.

import scipy.sparse
from networkx import minimum_spanning_tree
from networkx import from_scipy_sparse_matrix
from scipy.io import loadmat
import time

def find_min_spanning_tree(A):
	"""
		Input:
			A : Adjecency matrix in scipy.sparse format.
		Output:
			T : Minimum spanning tree.
			run_time : Total runtime to find minimum spanning tree 



	"""
	start = time.time()
	G = from_scipy_sparse_matrix(A)
	T = minimum_spanning_tree(G)
	run_time = time.time()-start
	return T, run_time

# Read Data from pre-processed data file
data = loadmat('../Data/preprocessed_dcg_sparse.mat')

# Separate Adjecency matrix
adj = data['adjecency_matrix']

# Get Minimum Spanning Tree
tree, run_time = find_min_spanning_tree(adj)

# *** USE T.edges TO ACCESS EDGES IN THE RETURNED MINIMUM SPANNING TREE.
# *** USE T.nodes TO ACCESS NODES IN THE RETURNED MINIMUM SPANNING TREE.