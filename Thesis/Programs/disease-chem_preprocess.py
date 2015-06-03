import pandas as pd
import numpy as np
import scipy.sparse
import scipy.io
import csv
from sklearn.preprocessing import normalize

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

	while (r1-r).dot((r1-r).T) > 1e-4 :
		r = r1
		r1 = c*(r*W) + (1-c)*q
		i= i+1

	print 'Reached to the convergence in %s iterations.'%i
	return r1

def preprocess_disease_chem(filename):
	"""
		Input:
			filename: Name of the CSV file containing disease chemical association.
		Output:
			Sparse matrix of disease_chem interaction.

		Example:
			Input-
			ID	Name Chem
			1	abc1 chm_1
			2	abc2 chm_2
			3	abc3 chm_3
			:	:	 :	
			n 	abcn chm_n

			Output-
			Sparse matrix:
			(0,0) chm_1
			(0,1) chm_2
			:	  :
	"""
	# Read Dataset
	data = pd.read_csv('../Data/disease_chemical_associations.csv')
	# Separate required columns
	doid= data['DOID']
	diseases = data['Disease_Name']
	chid = data['Chemical_ID']
	chemicals = data['Chemical_Name']

	# Find Unique Diseases and Chemicals

	diseases_unique = np.unique(doid).tolist()
	chemicals_unique = np.unique(chid).tolist()

	n_diseases = len(diseases_unique)
	n_chemicals = len(chemicals_unique)
	n = n_diseases+n_chemicals

	dc = scipy.sparse.lil_matrix((n,n),dtype=float)

	for i in range(0, len(doid)):
		id1 = diseases_unique.index(doid[i])
		id2 = chemicals_unique.index(chid[i])
		dc[id1, n_diseases+id2] = 1.0
		dc[n_diseases+id2, id1] = 1.0

	dc = dc.tocoo()
	return n, dc


n, A = preprocess_disease_chem('temp')
q = scipy.sparse.lil_matrix((1,n),dtype=float)
q[0,0]=1.0

r1 = rwr_algo(q, 0.9, A)
print r1
scipy.io.savemat('ranking.mat',{'ranking_vector':r1})