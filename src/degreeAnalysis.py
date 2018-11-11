import snap
import numpy as np
import matplotlib.pyplot as plt 


def genErdosRenyi(G):
    """
    G: Graph to generate Erdos Renyi graph based on its properties

    return type: snap.PUNGraph
    return: Erdos-Renyi graph with N nodes and E edges
    """
    N = G.GetNodes()
    E = G.GetEdges()

    Graph = snap.TUNGraph.New()
    for i in range(N):
        Graph.AddNode(i)
        
    i = 0
    while i < E:
        edge = (np.random.randint(0, E), np.random.randint(0, E))
        if Graph.IsNode(edge[0]) and Graph.IsNode(edge[1]):
            if not Graph.IsEdge(edge[0], edge[1]) and edge[0] != edge[1]:
                Graph.AddEdge(edge[0], edge[1])
                i += 1

    # print "Erdos edges: {}".format(Graph.GetEdges())
    # print "Erdos nodes: {}".format(Graph.GetNodes())
    return Graph

def getDataPointsToPlot(Graph):
    """
    :param - Graph: snap.PUNGraph object representing an undirected graph

    return values:
    X: list of degrees
    Y: list of frequencies: Y[i] = fraction of nodes with degree X[i]
    """
    X, Y = [], []
    deg = [node.GetDeg() for node in Graph.Nodes()]
    X = [i for i in range(max(deg))]
    Y = [float(deg.count(i))/Graph.GetNodes() for i in X]
    X,Y = zip(*[(x,y) for x,y in zip(X,Y) if y > 0])
    print X, Y
    return X, Y

G_LS174t = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/LS174T_Edgelist_spatialGraph_RIN.txt", 0, 1)
G_SW1222 = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/SW122_Edgelist_spatialGraph_RIN.txt", 0, 1)
X_LS174T, Y_LS174T = getDataPointsToPlot(G_LS174t)
X_SW1222, Y_SW1222 = getDataPointsToPlot(G_SW1222)

#Erdos graph for LS174t
G_ER_LS = genErdosRenyi(G_LS174t)
X_ER_LS, Y_ER_LS = getDataPointsToPlot(G_ER_LS)

#Erdos graph for SW1222
G_ER_LS = genErdosRenyi(G_SW1222)
X_ER_SW, Y_ER_SW = getDataPointsToPlot(G_ER_LS)

# plt.scatter(X_LS174T, Y_LS174T, color = 'y', label = 'LS174T')
# plt.scatter(X_ER_LS, Y_ER_LS, color = 'b', label = 'LS174T - Erdos')

plt.scatter(X_SW1222, Y_SW1222, color = 'y', label = 'SW1222')
plt.scatter(X_ER_SW, Y_ER_SW, color = 'b', label = 'SW1222 - Erdos')


# plt.scatter(X_LS174T, [np.log(x) for x in Y_LS174T], color = 'y', label = 'LS174T')
# plt.scatter(X_ER_LS, [np.log(x) for x in Y_ER_LS], color = 'b', label = 'LS174T - Erdos')
# plt.scatter(X_SW1222, [np.log(x) for x in Y_SW1222], color = 'r', label = 'SW1222')

#plt.yscale('log')

plt.xlabel('Node Degree (log)')
plt.ylabel('Count of Nodes with a Given Degree')
plt.title('Degree Distribution Different Tumours')
plt.legend()
plt.show()







# print zip(X_LS174T, Y_LS174T)
# G_LS174t = snap.LoadEdgeList(snap.PNEANET, "../data/Edgelist/LS174T_Edgelist_spatialGraph_RIN.txt", 0, 1)
# X_LS174T, Y_LS174T = getDataPointsToPlot(G_LS174t)
# print zip(X_LS174T, Y_LS174T)
