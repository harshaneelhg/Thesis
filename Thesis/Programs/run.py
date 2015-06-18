# This program is has code to run the algorithms after preprocessing.
# This progrm takes input for the query node, creates query vector and
# calls RWR algorithm to find Disease-Gene_chemical associations.

import numpy as np 
import scipy.io
import scipy.sparse
import json
from sklearn.preprocessing import normalize
import time

###################################### Function Definitions ########################################


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

############################################ Main Program ##########################################

# Record start time.
start = time.time()

# Read preprocessed data.
data = scipy.io.loadmat('../Data/preprocessed_dcg_sparse.mat')

A = data['adjecency_matrix']
n = int(data['number_of_nodes'])
n_diseases = int(data['number_of_diseases'])
n_chemicals = int(data['number_of_chemicals'])
n_genes = int(data['number_of_genes'])
with open('../Data/diseases_dict.json') as f:
    diseases_dict = dict(json.load(f))
with open('../Data/chem_dict.json') as f:
    chem_dict = dict(json.load(f))
with open('../Data/disease_id.json') as f:
    disease_id = dict(json.load(f))
with open('../Data/chem_id.json') as f:
    chem_id = dict(json.load(f))
with open('../Data/gene_id.json') as f:
    gene_id = dict(json.load(f))

# Initialize query vector.
q = scipy.sparse.lil_matrix((1,n),dtype=float)

# Input query node 'k'.
k = 1000

# Construct sparse query vector.
q[0,k]=1.0
print 'Executing query...'
r1 = rwr_algo(q, 0.9, A)

# Save Ranking results
scipy.io.savemat('../Data/ranking.mat',{'ranking_vector':r1})

# Display query information.
if k < n_diseases:
	print '\n***You are querying for a disease : %s, %s\n' % (disease_id[str(k)], diseases_dict[disease_id[str(k)]])
elif k >= n_diseases and k < n_diseases+n_chemicals:
	print '\n***You are qurying for a chemical : %s, %s\n' % (chem_id[str(k)], chem_dict[chem_id[str(k)]])
else:
	print '\n***You are qurying for a gene : %s\n' % gene_id[k]

# Separate ranking vector into diseases, chemicals and genes.
r1 = r1.todense()
r1 =r1.tolist()[0]
d = r1[0:n_diseases]
c = r1[n_diseases:n_diseases+n_chemicals]
g = r1[n_diseases+n_chemicals:n]

# Sort Ranking results.
sorted_d = [[i[0],i[1]] for i in sorted(enumerate(d), key=lambda x:x[1])]
sorted_c = [[i[0],i[1]] for i in sorted(enumerate(c), key=lambda x:x[1])]
sorted_g = [[i[0],i[1]] for i in sorted(enumerate(g), key=lambda x:x[1])]

# Print most relevant results for Diseases, Genes and Chemicals.
print 'Five most closely associated Diseases are:'
last = len(sorted_d)-1
for i in range(0,5):
	print ('    '+str(i+1)+': %s, %s\n       Score: %s') % (disease_id[str(sorted_d[last-i][0])], diseases_dict[disease_id[str(sorted_d[last-i][0])]], sorted_d[last-i][1])

print 'Five most closely associated Chemicals are:'
last = len(sorted_c)-1
for i in range(0,5):
	print ('    '+str(i+1)+': %s, %s\n       Score: %s') % (chem_id[str(n_diseases+sorted_c[last-i][0])], chem_dict[chem_id[str(n_diseases+sorted_c[last-i][0])]], sorted_c[last-i][1])

print 'Five most closely associated Genes are:'
last = len(sorted_g)-1
for i in range(0,5):
	print ('    '+str(i+1)+': %s\n       Score: %s') % (gene_id[str(n_chemicals+n_diseases+sorted_g[last-i][0])], sorted_g[last-i][1])

# Print elasped time.
print '\nElapsed time: %s seconds' % (time.time()-start) 