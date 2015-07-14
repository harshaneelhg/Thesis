import flask, flask.views
import json
import scipy.io
import scipy.sparse
import imp
import time
import networkx as nx
from networkx.readwrite import json_graph
import pdb
import random

app = flask.Flask(__name__)

# Create secret key for POST responses. This key must be random and must be secret.
app.secret_key= 'RandomWalk'

# List of session variables.
session = []

# Load processed data. Done once when application starts.
data = scipy.io.loadmat('../Data/preprocessed_dcg_sparse_new.mat')
r = imp.load_source('rwr_algo', '../Programs/rwr.py')

# Separate members from .mat file.
A = data['adjecency_matrix']
n = int(data['number_of_nodes'])
n_diseases = int(data['number_of_diseases'])
n_chemicals = int(data['number_of_chemicals'])
n_genes = int(data['number_of_genes'])

# Load Graph
print 'Loading graph...'
G = ""
G=nx.read_gpickle("../Data/dcg_graph.gpickle")
print 'Graph loaded...'

# Read dictionaries which keep track of IDs and names.
with open('../Data/diseases_dict_new.json') as f:
    diseases_dict = dict(json.load(f))
with open('../Data/chem_dict_new.json') as f:
    chem_dict = dict(json.load(f))
with open('../Data/gene_dict_new.json') as f:
    gene_dict = dict(json.load(f))
with open('../Data/disease_id_new.json') as f:
    disease_id = dict(json.load(f))
with open('../Data/chem_id_new.json') as f:
    chem_id = dict(json.load(f))
with open('../Data/gene_id_new.json') as f:
    gene_id = dict(json.load(f))
with open('../Data/disease_rev_lookup.json') as f:
    disease_rev_lookup = dict(json.load(f))
with open('../Data/chem_rev_lookup.json') as f:
    chem_rev_lookup = dict(json.load(f))
with open('../Data/gene_rev_lookup.json') as f:
    gene_rev_lookup = dict(json.load(f))

# Prepare lists for dropdown menus for users to select.
diseases = []
chemicals = []
genes = []

for d in disease_id:
	diseases.append((d,disease_id[d]+' | '+diseases_dict[disease_id[d]]))
for c in chem_id:
	chemicals.append((c,chem_id[c]+' | '+chem_dict[chem_id[c]]))
for g in gene_id:
	genes.append((g,str(gene_id[g])+' | '+gene_dict[str(gene_id[g])]))

def get_name(node_id):
	if str(node_id) in disease_id:
		return 'Disease: '+str(disease_id[str(node_id)])
	elif str(node_id) in chem_id:
		return 'Chemical: '+str(chem_id[str(node_id)])
	elif str(node_id) in gene_id:
		return 'Gene: '+str(gene_id[str(node_id)])
	else:
		return 'error'

def get_json(tree, node, visited):
	if tree.degree(node) == 1:
		return {"name" : get_name(node)}
	else:
		sub_tree = {"name" : get_name(node), "children":[]}
		neighbors = tree.neighbors(node)
		for neighbor in neighbors:
			if neighbor not in visited:
				visited.append(neighbor)
				child_tree=get_json(tree, neighbor, visited)
				sub_tree["children"].append(child_tree)
		return sub_tree

def get_shortest_path(source, target):
	P = nx.shortest_path(G, source = source, target = target) 
	return P

def get_spanning_tree(graph):
	T = nx.minimum_spanning_tree(G)
	return T

def get_connections(graph):
	return get_spanning_tree(graph)

def path_to_graph(path):
	G = nx.Graph()
	G.add_path(path)
	return G

def get_links(st_links):
	edges = nx.edges(st_links)
	l = []
	for edge in edges:
		l.append({'source': str(get_name(edge[0])), 'target':str(get_name(edge[1]))})
	return l

# Create a view for homepage.
class Home(flask.views.MethodView):
	
	# Define GET method.
	def get(self):
		if 'time' in session:
			session.pop()
		return flask.render_template('index.html', genes = genes, diseases = diseases, chemicals = chemicals, session= session)

	# Define POST method.
	def post(self):

		# Record start time.
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

		# Find and store most relevant results for Diseases, Genes and Chemicals.
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
			gene_display.append([gene_id[str(n_chemicals+n_diseases+sorted_g[last-i][0])], gene_dict[str(gene_id[str(n_chemicals+n_diseases+sorted_g[last-i][0])])], sorted_g[last-i][1]])

		end = time.time()

		# Add Query-Time in session.
		session.append('time')

		# Return calculated values to web-page for rendering.
		return flask.render_template('index.html', 
									  genes = genes, 
									  diseases = diseases, 
									  chemicals = chemicals,
									  disease_display = disease_display,
									  chem_display = chem_display,
									  gene_display = gene_display,
									  time = end-start,
									  session = session
									  )

# Create a view for connection analysis page.
class Connections(flask.views.MethodView):
	
	# Define GET method.
	def get(self):
		if 'path' in session:
			session.remove('path')
		return flask.render_template('connections.html', session= session)

	def post(self):
		if 'find_path' in flask.request.form:
			session = []
			session.append("path")

			source_node = flask.request.form['source_node']
			target_node = flask.request.form['target_node']
			if str(source_node) in disease_rev_lookup:
				source_id = int(disease_rev_lookup[source_node])
			elif str(source_node) in chem_rev_lookup:
				source_id = int(chem_rev_lookup[source_node])
			elif str(source_node) in gene_rev_lookup:
				source_id = int(gene_rev_lookup[source_node])
			else:
				session.remove("path")
				session.append("error")
				return flask.render_template('connections.html', session= session)
			if str(target_node) in disease_rev_lookup:
				target_id = int(disease_rev_lookup[target_node])
			elif str(target_node) in chem_rev_lookup:
				target_id = int(chem_rev_lookup[target_node])
			elif str(target_node) in gene_rev_lookup:
				target_id = int(gene_rev_lookup[target_node])
			else:
				session.remove("path")
				session.append("error")
				return flask.render_template('connections.html', session= session)
			path = get_shortest_path(source_id, target_id)
			tree = nx.Graph()
			tree.add_path(path)
			#nodes = tree.nodes()
			#pdb.set_trace()
			#root = nodes[len(path)/2]
			#while tree.degree(root) == 1:
			#	x = random.randint(0,len(path)-1)
			#	root = nodes[x]
			#H = get_json(tree, root, [root])
			#json_tree = json.dumps(H).encode('utf8')
			return flask.render_template('connections.html', session= session, tree_links=get_links(tree))
		#pdb.set_trace()
		if 'find_connections' in flask.request.form:
			session = []
			session.append('connection')

			ids = str(flask.request.form['input_ids'])

			ids = ids.split('\n')
			#pdb.set_trace()
			list_ids = []
			for i in ids:
				if str(i.strip()) in disease_rev_lookup:
					list_ids.append(int(disease_rev_lookup[str(i.strip())]))
				elif str(i.strip()) in chem_rev_lookup:
					list_ids.append(int(chem_rev_lookup[str(i.strip())]))
				elif str(i.strip()) in gene_rev_lookup:
					list_ids.append(int(gene_rev_lookup[str(i.strip())]))
				else:
					pdb.set_trace()
					session.remove('connection')
					session.append('error')
					return flask.render_template('connections.html', session= session)
			#pdb.set_trace()
			st_graph = nx.Graph()
			#p1= [1000,1001,1002,1003,1004]
			#p2 = [2000,1001,2002,2004,1003]
			#st_graph.add_path(p1)
			#st_graph.add_path(p2)
			for a in list_ids:
				for b in list_ids:
					P = get_shortest_path(a,b)
					st_graph.add_path(P)
			nodes = st_graph.nodes()
			#pdb.set_trace()
			root = nodes[len(nodes)/2]
			while st_graph.degree(root) == 1:
				x = random.randint(0,len(nodes)-1)
				root = nodes[x]
			steiner_tree = nx.minimum_spanning_tree(st_graph)
			H = get_json(st_graph, root, [root])
			#json_tree = json.dumps(H).encode('utf8')
			st_links = get_links(steiner_tree)
			return flask.render_template('connections.html', session= session, tree_links=st_links) 

# Add links in the context of web application.
app.add_url_rule(
				 '/random_walk_simulator',
				 view_func=Home.as_view('index'),
				 methods=['GET', 'POST']
				)
app.add_url_rule(
				 '/connection_analysis',
				 view_func=Connections.as_view('connections'),
				 methods=['GET', 'POST']
				)
app.debug = True
app.run()