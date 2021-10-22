import numpy as np

from kmeans import write_network, clean_network
import matplotlib.pyplot as plt
import networkx as nx

while True:
    user = input("Do you want to clean the network? (Y/N) ")
    if user == "Y":
        cleaned = True
        break
    elif user == "N":
        cleaned = False
        break

# Preprocess and retrieve network
G = clean_network("4932.protein.links.v11.5.txt", cleaned, str)

# Get most important nodes according to Page Rank
measure = nx.pagerank(G)
top10 = [[measure[k],k] for k in measure.keys()] # Choose the 10 largest values
top10.sort(reverse=True)
print("\n Centrality Measure: Page Rank")
for idx,pair in enumerate(top10):
    if pair[1] == "4932.YJR104C":
        print(str(idx+1),": \t is node ",pair[1]," with value:\t",pair[0])
    # else:
    #     print("Not SOD1", end=", ")
for idx,pair in enumerate(top10[:10]):
    print(str(idx+1),": \t is node ",pair[1]," with value:\t",pair[0])

# generate array of the centrality measure
G_plot = G
PageRankS = list(nx.pagerank(G_plot).values())
plt.title("Centrality measure for all nodes")
plt.plot(PageRankS)
plt.show()

# compute degree sequence
degS=[G.degree()[node] for node in list(G.nodes())]
degS.sort()
degS=np.array(degS)
mean_degree = np.average(degS)
print('\nMean degree:',mean_degree)

# subset of high degree nodes
highDegreePageRank = []
# threshold = mean_degree
threshold = 150
for i in range(len(degS)):
    if degS[i] >= threshold:
        highDegreePageRank.append(PageRankS[i])

print("high degrees: " + str(highDegreePageRank))
