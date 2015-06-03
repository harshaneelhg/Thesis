#!/usr/bin/env python

import numpy as np
import scipy.io



#========================================== FUNCTION DEFINITIONS ==========================================
#==========================================================================================================

def cross_rank(e, c, a, A_norm, Y_norm):
	"""
		Input:
		e: Aggregated query vector for all the networks.
		c: Restart probability.
		a: Regularization parameter for cross-network consistancy
		A_norm: Aggregated digonal block matrix of symmetrically normalized matrices of every network.
		Y_norm: Aggregated digonal block matrix to track common nodes in multiple networks.

		Output:
		r: Relevancy with every other node.

		This function implements basic version of CrossRank algorithm.
	"""
	e = np.array(e)
	A_norm = np.matrix(A_norm)
	Y_norm = np.matrix(Y_norm)

	r = e
	#r1 = ((c/(1+2*a))*A_norm + (2*a/(1+2*a))*Y_norm).(dot(r)) + ((1-c)/(1+2*a))*e

	while np.linalg.norm(r1-r) > 1e-3 :
		r = np.array(r1)
		#r1 = ((c/(1+2*a))*A_norm + (2*a/(1+2*a))*Y_norm).(dot(r)) + ((1-c)/(1+2*a))*e

	return np.array(r1).reshape(-1,).tolist()


def cross_rank_preprocess_A(A):
	"""
		Input:
		A: Data of Networks of Networks in .mat format 
		  (Set of adjecency matrices of individual networks)

		Output:
		A_norm: Preprossed digonal matrix A_norm.
	"""
	networks=[]
	index = 0
	for a in A:
		d = sum(a.todense())
		D = np.array(d)**0.5
		D = np.diag(D[0].tolist())
		Ai = scipy.sparse.coo_matrix(D)*a*scipy.sparse.coo_matrix(D)
		networks.append(Ai)
		if index%20 == 0:
			print index+1
		index = index + 1
	#print networks
	scipy.io.savemat('networks.mat',{'networks':networks})
	A_norm = scipy.linalg.block_diag(*networks)
	return A_norm


def cross_rank_preprocess_Y(A_ID, G):
	"""
		Input:
		A_ID: The corresponding IDs of domain-specific networks in A.
		G: The adjacency matrix of the main network

		Output:
		Y_norm: Aggregated digonal block matrix to track common nodes in multiple networks.
	"""
	for a in A_ID:
		for b in A_ID:
			v = intersects(a,b)
			d = np.ones(len(v[1]))
			x = len(v[1])
			y = len(v[2])
			Oij = scipy.sparse.coo_matrix((d, (ia,ib)), shape=(x,y))
			


def intersects(A,B):
	"""
		Input:
		A and B: lists of numbers.

		Output:
		O = [C, ia, ib] where,
			C = Common data.
			ia = indices of A for common data.
			ib = indices of B for common data.
	"""
	C = set(A).intersection(B)
	C = list(C)

	ia = [A.index(k) for k in C]
	ib = [B.index(k) for k in C]

	return [C, ia, ib]


#============================================== MAIN PROGRAM ==============================================
#==========================================================================================================
pass
#data = scipy.io.loadmat('../Data/DBLP_NoN.mat')
#A = data['CoAuthorNets'][0]
#A_norm = cross_rank_preprocess_A(A)
#scipy.io.savemat('A_norm.mat',{'A_norm':A_norm})
