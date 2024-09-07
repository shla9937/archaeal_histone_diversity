from netgraph import ArcDiagram
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap

# Create a graph
G = nx.Graph()

# Add nodes and edges
for edge, pairs in pairs_dict.items():
    if -1 in edge: 
        continue
    elif len(edge) == 1:  # Self-loop case
        if pairs:  # Check if there are data for the self-loop
            G.add_edge(edge[0], edge[0], weight=len(pairs)*.1)
    else:
        G.add_edge(edge[0], edge[1], weight=len(pairs)*.1)

# Create node labels
node_labels = {node: f"Clust {node}\n{len(proteins[proteins['Cluster'] == str(node)])}" for node in G.nodes()}
#node_sizes = {node: len(proteins[proteins['Cluster'] == str(node)])*0.0015 for node in G.nodes()}

edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]

# Add labels for edges
edge_labels = {}
for edge in G.edges():
    if edge[0] == edge[1]:
        if (edge[0],) in pairs_dict:
            edge_labels[edge] = len(pairs_dict[(edge[0],)])  # Handle self-loops with data
    else:
        try:
            edge_labels[edge] = len(pairs_dict[edge])
        except: 
            None
# Define the colors
light_gray = mcolors.to_rgba('#eeeeee')
black = mcolors.to_rgba('#000000')
min_weight = min(weights)
max_weight = max(weights)
scaled_weights = [(w - min_weight) / (max_weight - min_weight) for w in weights]
colors = [mcolors.to_hex(tuple((1 - gray) * np.array(light_gray) + gray * np.array(black))) for gray in scaled_weights]
edge_colors = {k: v for k, v in zip(edges, colors)}

# Convert rgb to hex
def parse_rgb(rgb_string):
    match = re.match(r'rgb\((\d+), (\d+), (\d+)\)', rgb_string)
    if match:
        try:
            # Convert each RGB value to float and divide by 255
            r, g, b = map(lambda x: int(x) / 255, match.groups())
            return r, g, b
        except:
            None
    else:
        raise ValueError(f"Invalid RGB string: {rgb_string}")

# Convert color strings to RGB tuples and then to hexadecimal format
hex_colors = [mcolors.to_hex(parse_rgb(rgb)) for rgb in tab_20_pastel]
node_colors = {k: hex_colors[k+1] for k in G.nodes()}

# Create an arc diagram with node and edge labels
arc_diagram = ArcDiagram(G, node_labels=node_labels, edge_labels=edge_labels, node_color=node_colors,
                         node_order=set(node_labels), edge_layout_kwargs=dict(selfloop_radius=-.05),
                         edge_color=edge_colors, edge_width=1, node_size=5.5,
                         edge_label_fontdict={'size':8, 'bbox': {'boxstyle': 'round', 'visible': False},
                                              'horizontalalignment': 'center', 
                                              'verticalalignment': 'bottom'},
                         node_label_fontdict={'size':8})

# Show the plot
plt.title('Co-occurrence of archaeal histones')
plt.axis('off')
#plt.savefig('figures/co-occur_2024_archaeal_histones.png', dpi=300)
plt.show()