#!/usr/bin/env python

# This piece of code implements pre-processing steps to generate
# sparse matrix from disease-chemical and disease-gene asscociation
# data storeed in two separate csv files.

import pandas as pd
import numpy as np
import scipy.sparse
import scipy.io
import csv
import json
from sklearn.preprocessing import normalize
import time


###################################### Function Definitions ########################################

def preprocess_disease_chem():
	"""
		Input:
			filename: Name of the CSV file containing disease-chemical association
			and disease-gene association.
		Output:
			n : Size of adjecency matrix.
			n_diseases : Number of unique diseases.
			n_chemicals : Number of unique chemicals.
			n_genes : Number of unique genes.
			A : Adjecency matrix of disease-chemical-gene association.
			diseases_dict : Dictionary of disease names and their DOIDs.
			chem_dict : Dictionary of chemical names and their CHIDs.
			disease_id : Dictionary of column/row numbers of diseases and their DOID.
			chem_id : Dictionary of column/row numbers of chemicals and their CHID.
			gene_id : Dictionary of column/row numbers of genes and gene names.

		Block matrix structure for Disease-Chemical-Genes association:
		(0 indicates block of zeros in block matrix)

				D 		C 		G 
			-------------------------
		    |	 	|		|		|
		D	|	0	|		|		|
			|		|		|		|
			-------------------------
		 	|		|		|		|
		C	|		|	0	|	0	|
			|		|		|		|
			-------------------------
		 	|		|		|		|
		G	|		|	0	|	0	|
			|		|		|		|
			-------------------------
	"""
	
	# Read Dataset.
	print 'reading dataset...'
	data_dc = pd.read_csv('../Data/disease_chemical_associations.csv')
	data_dg = pd.read_csv('../Data/disease_gene_associations.csv')
	
	# Separate required columns.

	doid_dc= data_dc['DOID'].tolist()
	diseases_dc = data_dc['Disease_Name'].tolist()
	chid = data_dc['Chemical_ID'].tolist()
	chemicals = data_dc['Chemical_Name'].tolist()

	doid_dg= data_dg['DOID'].tolist()
	diseases_dg = data_dg['Disease_Name'].tolist()
	genes = data_dg['Gene_Symbol'].tolist()
	
	# Find Unique Diseases,
	# Create disease-name dictionary,
	# Create disease-id dictionary.

	diseases_unique = np.unique(doid_dc+doid_dg).tolist()
	n_diseases = len(diseases_unique)
	temp_diseases = diseases_dc+diseases_dg
	temp_doid = doid_dc+doid_dg
	diseases_dict = {}
	disease_id = {}
	for d in diseases_unique:
		i = temp_doid.index(d)
		diseases_dict[temp_doid[i]] = temp_diseases[i] 
		disease_id[diseases_unique.index(d)] = d

	# Find Unique Chemicals,
	# Create chemical-name dictionary,
	# Create chemical-id dictionary.

	chemicals_unique = np.unique(chid).tolist()
	n_chemicals = len(chemicals_unique)
	chem_dict= {}
	chem_id = {}
	for c in chemicals_unique:
		i = chid.index(c)
		chem_dict[chid[i]]=chemicals[i]
		chem_id[chemicals_unique.index(c)+n_diseases] = c

	# Find Unique Genes,
	# Create gene-id dictionary.

	genes_unique = np.unique(genes).tolist()
	n_genes = len(genes_unique)
	gene_id = {}
	for g in genes_unique:
		gene_id[genes_unique.index(g)+n_diseases+n_chemicals] = g

	n = n_diseases+n_chemicals+n_genes

	# Create sparse matrix of disease-chemical-gene association.
	print 'Constructing sparse matrix...'
	dcg = scipy.sparse.lil_matrix((n,n),dtype=float)

	for i in range(0, len(doid_dc)):
		id1 = diseases_unique.index(doid_dc[i])
		id2 =  n_diseases+chemicals_unique.index(chid[i])
		dcg[id1, id2] = 1.0
		dcg[id2, id1] = 1.0

	for i in range(0, len(doid_dg)):
		id1 = diseases_unique.index(doid_dg[i])
		id2 = n_diseases+n_chemicals+genes_unique.index(genes[i])
		dcg[id1, id2] = 1.0
		dcg[id2, id1] = 1.0

	dcg = dcg.tocoo()
	return n, n_diseases, n_chemicals, n_genes, dcg, diseases_dict, chem_dict, disease_id, chem_id, gene_id

############################################ Main Program ##########################################

# Record start time.
start = time.time()

# Pre-process data.
print 'Pre-processing data...'
n, n_diseases, n_chemicals, n_genes, A, diseases_dict, chem_dict, disease_id, chem_id, gene_id = preprocess_disease_chem()

save = {
		'adjecency_matrix': A,
		'number_of_nodes': n,
		'number_of_diseases': n_diseases,
		'number_of_chemicals': n_chemicals,
		'number_of_genes': n_genes,
		}
scipy.io.savemat('../Data/preprocessed_dcg_sparse.mat', save)

# Save dictionaries to JSON files.
with open('../Data/diseases_dict.json', 'w') as outfile:
     json.dump(diseases_dict, outfile)
with open('../Data/chem_dict.json', 'w') as outfile:
     json.dump(chem_dict, outfile)
with open('../Data/disease_id.json', 'w') as outfile:
     json.dump(disease_id, outfile)
with open('../Data/chem_id.json', 'w') as outfile:
     json.dump(chem_id, outfile)
with open('../Data/gene_id.json', 'w') as outfile:
     json.dump(gene_id, outfile)

# Print elasped time.
print '\nElapsed time: %s seconds' % (time.time()-start)