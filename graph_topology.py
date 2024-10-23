#written by Eli Niktab

import igraph as ig
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np

def parse_gfa_to_igraph(gfa_file):
    """
    Parses a GFA file and returns an igraph Graph object.
    """
    # Initialize a graph
    G = ig.Graph(directed=False)  # Change to directed=True if your GFA graph is directed
    nodes_set = set()  # To store unique nodes
    edges = []  # To store edges before adding them to the graph

    with open(gfa_file, 'r') as file:
        for line in file:
            if line.startswith('S'):  # Segment line (node)
                parts = line.strip().split('\t')
                node_id = parts[1]  # Extract node ID
                nodes_set.add(node_id)
            elif line.startswith('L'):  # Link line (edge)
                parts = line.strip().split('\t')
                source = parts[1]  # Source node
                target = parts[3]  # Target node
                nodes_set.add(source)
                nodes_set.add(target)
                edges.append((source, target))

    # Add all unique nodes to the graph
    G.add_vertices(list(nodes_set))

    # Map node IDs to indices
    node_id_map = {node_id: idx for idx, node_id in enumerate(G.vs['name'])}

    # Add all the edges using node indices
    edge_list = [(node_id_map[source], node_id_map[target]) for source, target in edges]
    G.add_edges(edge_list)

    return G

def calculate_centrality_measures(graph):
    """
    Calculates centrality measures for a graph and returns them as a dictionary.
    Handles NaN values by replacing them with zeros.
    """
    measures = {}
    # Degree centrality
    measures['degree'] = graph.degree()

    # Betweenness centrality
    betweenness = graph.betweenness()
    # Replace NaNs with zeros
    betweenness = [0 if np.isnan(b) else b for b in betweenness]
    measures['betweenness'] = betweenness

    # Closeness centrality
    closeness = graph.closeness()
    # Replace NaNs with zeros
    closeness = [0 if np.isnan(c) else c for c in closeness]
    measures['closeness'] = closeness

    # Eigenvector centrality
    eigenvector = graph.eigenvector_centrality()
    # Replace NaNs with zeros
    eigenvector = [0 if np.isnan(e) else e for e in eigenvector]
    measures['eigenvector'] = eigenvector

    return measures

def compare_centrality(measure1, measure2, measure_name):
    """
    Compares two centrality measures using the Mann-Whitney U test and prints the results.
    Handles NaN values by filtering them out.
    """
    # Convert measures to numpy arrays for easier handling
    measure1 = np.array(measure1)
    measure2 = np.array(measure2)

    # Filter out NaN values
    measure1 = measure1[~np.isnan(measure1)]
    measure2 = measure2[~np.isnan(measure2)]

    # Check if we have enough data to perform the test
    if len(measure1) == 0 or len(measure2) == 0:
        print(f"\nCannot perform statistical test for {measure_name} centrality due to insufficient data.")
        return

    # Perform Mann-Whitney U test
    stat, p_value = mannwhitneyu(measure1, measure2, alternative='two-sided')

    print(f"\nComparing {measure_name.capitalize()} Centrality:")
    print(f"Statistic: {stat:.4f}")
    print(f"P-value: {p_value:.4f}")

    if p_value < 0.05:
        print(f"The difference in {measure_name} centrality distributions is statistically significant.")
    else:
        print(f"No significant difference found in {measure_name} centrality distributions.")


def main():
    # Paths to your GFA files
    gfa_file1 = '/home/niktabel/workspace/media/gadi_g_te53/sj2852/blgrp/analyses/pangenome/analysis-v4/communities/pggb-corrected/A4GALT.pggb.out/A4GALT.fasta.c0e6cc8.11fba48.f0903f4.smooth.final.gfa'
    gfa_file2 = '/home/niktabel/workspace/media/gadi_g_te53/sj2852/blgrp/analyses/pangenome/analysis-v4/communities/pggb-corrected/ABCB6.pggb.out/ABCB6.fasta.c0e6cc8.11fba48.8ac1809.smooth.final.gfa'

    # Parse both GFA files into graphs
    graph1 = parse_gfa_to_igraph(gfa_file1)
    graph2 = parse_gfa_to_igraph(gfa_file2)

    # Calculate centrality measures for both graphs
    centrality_measures1 = calculate_centrality_measures(graph1)
    centrality_measures2 = calculate_centrality_measures(graph2)

    # Compare centrality measures
    #for measure_name in ['degree', 'betweenness', 'closeness', 'eigenvector']:
        #measure1 = centrality_measures1[measure_name]
        #measure2 = centrality_measures2[measure_name]
        #compare_centrality(measure1, measure2, measure_name)


if __name__ == '__main__':
    main()


