# This piece of code finds out minimum spanning tree for 
# Disease-Chemical-Gene association graph using 'Networkx'
# Python library.

import scipy.sparse
from networkx import minimum_spanning_tree
from networkx import from_scipy_sparse_matrix
from networkx.readwrite import json_graph
from scipy.io import loadmat
import time
import json
import os.path

def find_min_spanning_tree(A):
	"""
		Input:
			A : Adjecency matrix in scipy.sparse format.
		Output:
			T : Minimum spanning tree.
			run_time : Total runtime to find minimum spanning tree 

	"""
	# Record start time.
	start = time.time()

	# Check if graph is pre-processed, if yes then don't process it again.
	if os.path.exists('../Data/dcg_graph.json'):
		with open('../Data/dcg_graph.json') as data:
			d = json.load(data)
		G = json_graph.node_link_graph(d)

	# If graph is not preprocessed then convert it to a Graph and save it to a JSON file.
	else:
		G = from_scipy_sparse_matrix(A)
		data = json_graph.node_link_data(G)
		with open('../Data/dcg_graph.json', 'w') as outfile:
			json.dump(data, outfile)

	# Find MST.
	T = minimum_spanning_tree(G)

	#Record total Runtime
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

print 'Runtime: %s'%(run_time)