import flask
import networkx as nx
import json
from networkx.readwrite import json_graph

G = nx.barbell_graph(6,3)
# this d3 example uses the name attribute for the mouse-hover value,
# so add a name to each node
for n in G:
    G.node[n]['name'] = n
# write json formatted data
d = json_graph.node_link_data(G) # node-link format to serialize
# write json
json.dump(d, open('force/force.json','w'))
print('Wrote node-link JSON data to force/force.json')


app = flask.Flask(__name__, static_folder = 'force')

@app.route('/<path:path>')
def static_proxy(path):
    return app.send_static_file(path)
print('Go to http://localhost:8000/force.html')
app.run(port=8000)
