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

# Load Graph
print 'Loading graph...'
G = ""
G=nx.read_gpickle("../Data/dcg_graph.gpickle")
print 'Graph loaded...'

# Get name of the disease, chemical or gene for given index.
def get_name(node_id):
	if str(node_id) in disease_id:
		return 'Disease: '+str(disease_id[str(node_id)])
	elif str(node_id) in chem_id:
		return 'Chemical: '+str(chem_id[str(node_id)])
	elif str(node_id) in gene_id:
		return 'Gene: '+str(gene_id[str(node_id)])
	else:
		return 'error'

# BFS implementation to get json version of graph.
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

# Returns list of links in a graph.
def get_links(st_links, link_type, n_path):
	if link_type == 'path':
		edges = nx.edges(st_links)
		l = []
		for edge in edges:
			l.append({'source': str(get_name(edge[0]))+'('+str(n_path)+')', 'target':str(get_name(edge[1]))+'('+str(n_path)+')'})
		return l
	else:
		edges = nx.edges(st_links)
		l = []
		for edge in edges:
			l.append({'source': str(get_name(edge[0])), 'target':str(get_name(edge[1]))})
		return l

# Returns IDs for Diseases, Chemicals or Genes.
def get_key(query,query_type):
	if query_type == "disease":
		keys = diseases_dict.keys()
		for k in keys:
			if diseases_dict[k] == query:
				return k
	elif query_type == "chem":
		keys = chem_dict.keys()
		for k in keys:
			if chem_dict[k] == query:
				return k
	else:
		keys = gene_dict.keys()
		for k in keys:
			if gene_dict[k] == query:
				return k

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

	# Define post method.
	def post(self):
		# Handle request for shortest path.
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
			if nx.has_path(G, source_id, target_id):	
				path = nx.all_shortest_paths(G, source_id, target_id)
			else:
				path = [[source_id,source_id],[target_id,target_id]]
			selected = list()
			all_links = list()
			n_paths = 0
			for p in path:
				tree = nx.Graph()
				tree.add_path(p)
				l = get_links(tree,'path',n_paths+1)
				all_links = all_links + l
				n_paths = n_paths + 1
				selected.append(get_name(source_id)+'('+str(n_paths)+')')
				selected.append(get_name(target_id)+'('+str(n_paths)+')')
				if n_paths >= 3:
					break
			return flask.render_template('connections.html', session= session, tree_links=all_links, selected=selected)
		# Handle request for Steiner tree visualization.
		if 'find_connections' in flask.request.form:
			session = []
			session.append('connection')
			# Separate IDs from input textarea.
			ids = str(flask.request.form['input_ids'])
			ids = ids.strip()
			ids = ids.split('\n')
			selected = []
			list_ids = []
			for i in ids:
				# Reverse lookup dictionary to find indices of IDs in the graph.
				if str(i.strip()) in disease_rev_lookup:
					list_ids.append(int(disease_rev_lookup[str(i.strip())]))
				elif str(i.strip()) in chem_rev_lookup:
					list_ids.append(int(chem_rev_lookup[str(i.strip())]))
				elif str(i.strip()) in gene_rev_lookup:
					list_ids.append(int(gene_rev_lookup[str(i.strip())]))
				# Throw error if ID/s are inappropriate. 
				else:
					session.remove('connection')
					session.append('error')
					return flask.render_template('connections.html', session= session)
			# Initialize an empty graph.
			st_graph = nx.Graph()
			# Add shortest path to graph between every pair of node.
			isolated_paths = []
			for a in list_ids:
				for b in list_ids:
					if nx.has_path(G,a,b):
						P = nx.shortest_path(G,a,b)
						st_graph.add_path(P)
					else:
						isolated_paths.append([a,a])
						isolated_paths.append([b,b])
				selected.append(get_name(a))
			nodes = st_graph.nodes()
			# Find spanning tree of subgraph created by adding shortest paths every pair of nodes.
			steiner_tree = nx.minimum_spanning_tree(st_graph)
			for ip in isolated_paths:
				steiner_tree.add_path(ip)
			st_links = get_links(steiner_tree,'conn',0)
			return flask.render_template('connections.html', session= session, tree_links=st_links, selected=selected) 

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

# Add a route for AJAX call for search query.
@app.route('/get_info', methods=['GET'])
def get_info():
	query = flask.request.args['search']
	if query in diseases_dict:
		message = "Type : Disease\n"+"Name : "+str(diseases_dict[str(query)])+"\nID : "+str(query)
	elif query in diseases_dict.values():
		message = "Type : Disease\n"+"Name : "+str(query)+"\nID : "+ get_key(query,"disease")
	elif query in chem_dict:
		message = "Type : Chemical\n"+"Name : "+str(chem_dict[str(query)])+"\nID : "+str(query)
	elif query in chem_dict.values():
		message = "Type : Chemical\n"+"Name : "+str(query)+"\nID : "+ get_key(query,"chem")
	elif query in gene_dict:
		message = "Type : Gene\n"+"Name : "+str(gene_dict[str(query)])+"\nID : "+str(query)
	elif query in gene_dict.values():
		message = "Type : Gene\n"+"Name : "+str(query)+"\nID : "+ get_key(query,"gene")
	else:
		message = "Incorrect query!"
	return json.dumps({'message': message})

# Add route to populate autocomplete options for AJAX call.
@app.route('/get_available_options', methods=['GET'])
def get_available_options():
	option = flask.request.args['type']
	if option == "all":
		return json.dumps({'names' : diseases_dict.values()+chem_dict.values()+gene_dict.values(), 'ids' : diseases_dict.keys()+chem_dict.keys()+gene_dict.keys()})
	else:
		return json.dumps({'error' : "Invalid query option..."})


app.debug = True
app.run(port=6063)