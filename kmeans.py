from matplotlib import colors, cm
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict

from sklearn.cluster import KMeans, SpectralClustering

def write_network():
    proteins = {}
    data = {}
    i = 0
    with open("4932.protein.links.v11.5.txt") as f:
        f.readline() # First empty line
        while True:
            line = f.readline().strip()
            if line == "":
                break
            protein_1 = line.split()[0][5:]
            protein_2 = line.split()[1][5:]
            if (protein_1 == "YJR104C" and protein_1 not in proteins) or (protein_2 == "YJR104C" and protein_2 not in proteins):
                sod1_index = i
                print("SOD1's ID: " + str(i))
            strength = line.split()[2]
            if protein_1 not in proteins:
                proteins[protein_1] = i
                i += 1
            if protein_2 not in proteins:
                proteins[protein_2] = i
                i += 1
            data[(proteins[protein_1], proteins[protein_2])] = strength
    f.close()

    with open("network_cleaned.txt", "w") as f:
        for key in data:
            f.write(str(key[0]) + " " + str(key[1]) + " " + data[key] + "\n")
    f.close()

    return sod1_index, proteins

def clean_network(filename, cleaned, node_type):
    # Read in network
    if node_type == int:
        G = nx.read_weighted_edgelist(filename, comments="#", nodetype=int)
    elif node_type == str:
        G = nx.read_weighted_edgelist(filename, comments="#", nodetype=str)

    if not cleaned:
        return G

    # Delete low-confidence links
    threshold_score = 990
    for edge in G.edges:
        weight = list(G.get_edge_data(edge[0], edge[1]).values())
        if(weight[0] <= threshold_score):
            G.remove_edge(edge[0], edge[1])
        else:
            G.remove_edge(edge[0], edge[1])
            G.add_edge(edge[0], edge[1]) # To set value to 1 instead of confidence link

    # Read in essential proteins
    essential_proteins = []
    with open("essential_nodes.txt") as f:
        for line in f:
            essential_proteins.append(line.strip())
    f.close()
    print("No. of essential proteins removed: " + str(len(essential_proteins)))

    # Get network w/o essential proteins
    G_prime = nx.Graph.copy(G)
    for p in essential_proteins:
        p = "4932." + p
        if p in G_prime.nodes:
            G_prime.remove_node(p)
    G = G_prime

    # Take only largest connected component
    largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc)

    print("No. of nodes: " + str(G.number_of_nodes()))
    print("No. of edges: " + str(G.number_of_edges()))

    return G

def cluster_alg(G, alg_type, num_clusters, cleaned, print_results):
    adj_matrix = nx.adjacency_matrix(G)
    if not cleaned:
        adj_matrix /= 1000 # Normalize for faster processing
    # if print_results:
    #     print("Adjacency matrix:\n" + str(adj_matrix))

    k_clusters = num_clusters
    results = []
    algorithm = {}

    if alg_type == "kmeans":
        algorithm['kmeans'] = KMeans(n_clusters=k_clusters)
    elif alg_type == "spectral":
        algorithm['spectral'] = SpectralClustering(n_clusters=k_clusters, affinity="precomputed")  # no init
    for model in algorithm.values():
        model.fit(adj_matrix)
        results.append(list(model.labels_))
    if print_results:
        if alg_type == "kmeans":
            print("K-means clustering algorithm (which community each protein is in):\n" + str(results[0]))
        elif alg_type == "spectral":
            print("Spectral clustering algorithm (which community each protein is in):\n" + str(results[0]))

    return results[0], algorithm

def draw_communities(G, membership, pos, num_clusters):
    fig, ax = plt.subplots(figsize=(16, 9))

    # Convert membership list to a dict where key=club, value=list of students in club
    club_dict = defaultdict(list)
    for student, club in enumerate(membership):
        club_dict[club].append(student)

    # Normalize number of clubs for choosing a color
    norm = colors.Normalize(vmin=0, vmax=len(club_dict.keys()))

    for club, members in club_dict.items():
        nx.draw_networkx_nodes(G, pos,
                               nodelist=members,
                               node_color=cm.jet(norm(club)),
                               node_size=500,
                               alpha=0.8,
                               ax=ax)

    # Draw edges (social connections) and show final plot
    plt.title("Graph with " + str(num_clusters) + " clusters")
    nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax)
    plt.show()

def run_experiment(G, sod1_index, proteins, node_frequencies, alg_type, num_clusters, cleaned, print_results):
    # Use k-means or spectral clustering algorithm
    results, algorithms = cluster_alg(G, alg_type, num_clusters, cleaned, print_results)

    # Make results able to look up specific proteins
    indices = [i for i in range(len(results))]
    results = list(zip(results, indices))

    # Sort nodes by cluster
    clusters = {}
    for tup in results:
        cluster = tup[0]
        index = tup[1]
        if cluster not in clusters:
            clusters[cluster] = [index]
        else:
            clusters[cluster].append(index)
    if print_results:
        print("Clusters (which proteins are in each community):\n" + str(clusters))
    assert len(clusters) == num_clusters

    # Find the cluster which SOD1 is in
    for c in clusters:
        if clusters[c].count(sod1_index) > 0:
            sod1_community = c
            sod1_cluster = clusters[c]
        # if print_results:
        #     print("Community " + str(c) + " has " + str(len(clusters[c])) + " proteins")
    if print_results:
        print("SOD1 is in community no. " + str(sod1_community))

    # Get nodes in same community as SOD1
    same_community = []
    for node in sod1_cluster:
        protein = list(proteins.keys())[node]
        same_community.append("4932." + protein)
    if print_results:
        print("Proteins in community no. " + str(sod1_community) + ": " + str(same_community))
        print("Size of community no. " + str(sod1_community) + ": " + str(len(sod1_cluster)) + "\n")

    for node in same_community:
        if node not in node_frequencies:
            node_frequencies[node] = 1
        else:
            node_frequencies[node] += 1

    return node_frequencies

if __name__ == "__main__":
    # Get user input
    while True:
        user = input("Do you want to clean the network? (Y/N) ")
        if user == "Y":
            cleaned = True
            break
        elif user == "N":
            cleaned = False
            break
    while True:
        try:
            num_clusters = int(input("How many clusters do you want? (k>=2) "))
            break
        except ValueError:
            pass
    while True:
        alg_type = input("Which algorithm would you like to use? (kmeans/spectral) ")
        if alg_type == "kmeans" or alg_type == "spectral":
            break
    while True:
        try:
            num_experiments = int(input("How many experiments do you want to perform? "))
            break
        except ValueError:
            pass
    while True:
        user = input("Do you want to keep track of the clustering? (Y/N) ")
        if user == "Y":
            print_results = True
            break
        elif user == "N":
            print_results = False
            break

    # Rewrite network and possibly clean it
    print("##### Preprocessing network... #####")
    sod1_index, proteins = write_network()
    print("##### Cleaning network... #####")
    G = clean_network("network_cleaned.txt", cleaned, int)
    # G = nx.karate_club_graph()

    # Keep count of no. of frequencies each node in SOD1's community appears after all experiments
    node_frequencies = {}
    for i in range(num_experiments):
        if i == 0 or (i + 1) % 10 == 0 or i == num_experiments - 1:
            print("##### Running experiment " + str(i + 1) + "/" + str(num_experiments) + "... #####")
        node_frequencies = run_experiment(G, sod1_index, proteins, node_frequencies, alg_type, num_clusters, cleaned, print_results)

    # Find threshold to keep track of most frequent nodes in SOD1's community
    print("##### Results #####")
    print("Frequencies of nodes in SOD1's community after " + str(num_experiments) + " experiments: " + str(node_frequencies))
    print("Community of proteins: " + str(list(node_frequencies.keys())))
    print("Size of community: " + str(len(list(node_frequencies.keys()))))

    lengths = []
    for i in range(num_experiments):
        freq_dict_filtered = dict((k, v) for k, v in node_frequencies.items() if v > i)
        lengths.append(len(freq_dict_filtered))
    plt.plot(lengths)
    plt.title("Lengths:")
    plt.show()

    freq_filtered = {}
    for f in node_frequencies:
        if node_frequencies[f] >= num_experiments * (3 / 4):
            freq_filtered[f] = node_frequencies[f]
    print("Frequencies of proteins in SOD1's community after " + str(num_experiments) + " experiments (filtered): " + str(freq_filtered))
    print("Proteins in SOD1's community (filtered): " + str(list(freq_filtered.keys())))
    print("Size of community (filtered): " + str(len(list(freq_filtered.keys()))) + "\n")

    # Draw clusters (takes a while, and you CANNOT clean the network beforehand)
    # pos = nx.spring_layout(G)
    # draw_communities(G, algorithms['kmeans'].labels_, pos, num_clusters)
    # draw_communities(G, algorithms[''].labels_, pos, num_clusters)