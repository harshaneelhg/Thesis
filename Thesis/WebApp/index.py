import flask, flask.views
import json
import scipy.io

app = flask.Flask(__name__)

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
		return flask.render_template('index.html', genes = genes, diseases = diseases, chemicals = chemicals)

	def post(self):
		return flask.render_template('index.html', genes = genes, diseases = diseases, chemicals = chemicals)
		if 'disease-query' in flask.request.form :
			index = flask.request.form['disease']

app.add_url_rule(
				 '/',
				 view_func=Home.as_view('index'),
				 methods=['GET', 'POST']
				)
app.debug = True
app.run()