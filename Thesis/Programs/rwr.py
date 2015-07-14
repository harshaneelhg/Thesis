#!/usr/bin/env python

import numpy as np
import scipy.sparse
from sklearn.preprocessing import normalize

#========================================== FUNCTION DEFINITIONS ==========================================
#==========================================================================================================

def rwr_algo(q, c, W):
	"""
		Input:
		q: Sparse query vector.
		c: Restart probabilities.
		W: Sparse adjecency matrix.

		Output:
		r: Sparse relevancy vector with every other node.

		This function implements basic version of RWR to output
		Relevancy of query node with other nodes in a graph represented
		by adjecency matrix W.
	"""

	#Column Normalize Adjecency Matrix.

	W = normalize(W, norm='l1', axis=0)

	#Basic Random Walk with Restarts Algorithm.

	r = q
	r1 = c*(r*W) + (1-c)*q
	i=0

	while (r1-r).dot((r1-r).T) > 1e-5 :
		r = r1
		r1 = c*(r*W) + (1-c)*q
		i= i+1

	#print 'Reached to the convergence in %s iterations.'%i
	return r1

#============================================== MAIN PROGRAM ==============================================
#==========================================================================================================

#========================================= Uncomment for test run==========================================
"""
W = [[0.000,1.000,1.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000],
	 [1.000,0.000,1.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000],
	 [1.000,1.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000],
	 [1.000,0.000,1.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000],
	 [0.000,0.000,0.000,1.000,0.000,1.000,1.000,1.000,0.000,0.000,0.000,0.000],
	 [0.000,0.000,0.000,0.000,1.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000],
	 [0.000,0.000,0.000,0.000,1.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000],
	 [0.000,1.000,0.000,0.000,1.000,0.000,0.000,0.000,1.000,0.000,1.000,0.000],
	 [0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,1.000,0.000,0.000],
	 [0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,1.000,1.000],
	 [0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,1.000,0.000,1.000],
	 [0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,1.000,0.000]]
W1 = scipy.sparse.coo_matrix(W)

q =  [0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000]
q1 = scipy.sparse.coo_matrix(q)

c = 0.9

scores = rwr_algo(q1,c,W1)
scores = np.array(scores.todense()).reshape(-1,).tolist()

# Adjust precison upto 3 decimal places.
scores = [float("%.3f"%scores[i]) for i in range(0,len(scores))]
print np.matrix(scores).T

"""