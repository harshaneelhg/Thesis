# This piece of code finds out shortest path between two nodes 
# of Disease-Chemical-Gene association graph using 'Networkx'
# Python library.

import scipy.sparse
from networkx import shortest_path
from networkx import from_scipy_sparse_matrix
from networkx.readwrite import json_graph
from scipy.io import loadmat
import time
import json
import os.path

def find_shortest_path(A, source, target):
	"""
		Input:
			A : Adjecency matrix in scipy.sparse format.
			source : Source node.
			target : Destination node.
		Output:
			P : Shortest path.
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

	# Find shortest path.
	P = shortest_path(G, source = source, target = target)

	#Record total Runtime
	run_time = time.time()-start
	return P, run_time

# Read Data from pre-processed data file
data = loadmat('../Data/preprocessed_dcg_sparse.mat')

# Separate Adjecency matrix
adj = data['adjecency_matrix']

source = 100
target = 5000

# Get Minimum Spanning Tree
path, run_time = find_shortest_path(adj, source, target)

print 'Path: ',
print path
print 'Runtime: %s'%(run_time)