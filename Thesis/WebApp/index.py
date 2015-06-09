import flask, flask.views
import json
import scipy.io
import scipy.sparse
import imp
import time

app = flask.Flask(__name__)
app.secret_key= 'RandomWalk'
session = []
data = scipy.io.loadmat('../Data/preprocessed_dcg_sparse.mat')

r = imp.load_source('rwr_algo', '../Programs/rwr.py')

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

diseases = []
chemicals = []
genes = []

for d in disease_id:
	diseases.append((d,disease_id[d]+' | '+diseases_dict[disease_id[d]]))
for c in chem_id:
	chemicals.append((c,chem_id[c]+' | '+chem_dict[chem_id[c]]))
for g in gene_id:
	genes.append((g,gene_id[g]))

class Home(flask.views.MethodView):
	
	def get(self):
		if 'time' in session:
			session.pop()
		return flask.render_template('index.html', genes = genes, diseases = diseases, chemicals = chemicals, session= session)

	def post(self):
		start = time.time()
		index = None
		if 'disease_query' in flask.request.form :
			index = flask.request.form['disease']
		elif 'chemical_query' in flask.request.form :
			index = flask.request.form['chemical']
		else:
			index = flask.request.form['gene']
		q = scipy.sparse.lil_matrix((1,n),dtype=float)
		q[0,int(index)]=1.0
		r1 = r.rwr_algo(q, 0.9, A)
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

		disease_display=[]
		last = len(sorted_d)-1
		for i in range(0,5):
			disease_display.append([disease_id[str(sorted_d[last-i][0])], diseases_dict[disease_id[str(sorted_d[last-i][0])]], sorted_d[last-i][1]])

		chem_display=[]
		last = len(sorted_c)-1
		for i in range(0,5):
			chem_display.append([chem_id[str(n_diseases+sorted_c[last-i][0])], chem_dict[chem_id[str(n_diseases+sorted_c[last-i][0])]], sorted_c[last-i][1]])

		gene_display=[]
		last = len(sorted_g)-1
		for i in range(0,5):
			gene_display.append([gene_id[str(n_chemicals+n_diseases+sorted_g[last-i][0])], sorted_g[last-i][1]])

		session.append('time')

		return flask.render_template('index.html', 
									  genes = genes, 
									  diseases = diseases, 
									  chemicals = chemicals,
									  disease_display = disease_display,
									  chem_display = chem_display,
									  gene_display = gene_display,
									  time = time.time()-start,
									  session = session
									  )

app.add_url_rule(
				 '/',
				 view_func=Home.as_view('index'),
				 methods=['GET', 'POST']
				)
app.debug = True
app.run()