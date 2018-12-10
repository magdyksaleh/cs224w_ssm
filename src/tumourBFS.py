import snap
import collections
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

#Running BFS on the graph


def getEfficiency(G):
    NIdToDistH = snap.TIntH()
    shortestPaths = 0
    for i in tqdm(range(G.GetNodes())):
        snap.GetShortPath(G, i, NIdToDistH)
        shortestPaths += sum([NIdToDistH[item] for item in NIdToDistH])
        # if i%100==0: print (i, 1.0/shortestPaths)
    return 1.0/shortestPaths


def genErdosRenyi(N, E):
    """
    :param - N: number of nodes
    :param - E: number of edges

    return type: snap.PUNGraph
    return: Erdos-Renyi graph with N nodes and E edges
    """
    ############################################################################
    Graph = snap.TNGraph.New()
    for i in range(N):
        Graph.AddNode(i)
        
    i = 0
    while i < E:
        edge = (np.random.randint(0, E), np.random.randint(0, E))
        if Graph.IsNode(edge[0]) and Graph.IsNode(edge[1]):
            if not Graph.IsEdge(edge[0], edge[1]) and edge[0] != edge[1]:
                Graph.AddEdge(edge[0], edge[1])
                i += 1

    ############################################################################
    return Graph


def genCircle(N=5242):
    """
    :param - N: number of nodes

    return type: snap.PUNGraph
    return: Circle graph with N nodes and N edges. Imagine the nodes form a
        circle and each node is connected to its two direct neighbors.
    """
    ############################################################################
    # TODO: Your code here!
    Graph = snap.TNGraph.New()
    for i in range(N):
        Graph.AddNode(i)
    
    for i in range(N):
        Graph.AddEdge(i, (i+1)%N)

    ############################################################################
    return Graph


def connectNbrOfNbr(Graph, N=5242):
    """
    :param - Graph: snap.PUNGraph object representing a circle graph on N nodes
    :param - N: number of nodes

    return type: snap.PUNGraph
    return: Graph object with additional N edges added by connecting each node
        to the neighbors of its neighbors
    """
    ############################################################################
    for i in range(N):
        Graph.AddEdge(i, (i+2)%N)
    ############################################################################
    return Graph


def connectRandomNodes(Graph, M=4000):
    """
    :param - Graph: snap.PUNGraph object representing an undirected graph
    :param - M: number of edges to be added

    return type: snap.PUNGraph
    return: Graph object with additional M edges added by connecting M randomly
        selected pairs of nodes not already connected.
    """
    ############################################################################
    # TODO: Your code here!
    N = Graph.GetNodes()
    for i in range(M):
        edge = (np.random.randint(0, N), np.random.randint(0, N))
        if not Graph.IsEdge(edge[0], edge[1]) and edge[0] != edge[1]:
            Graph.AddEdge(edge[0], edge[1])
    ############################################################################
    return Graph


def genSmallWorld(N=5242, E=14484):
    """
    :param - N: number of nodes
    :param - E: number of edges

    return type: snap.PUNGraph
    return: Small-World graph with N nodes and E edges
    """
    Graph = genCircle(N)
    Graph = connectNbrOfNbr(Graph, N)
    Graph = connectRandomNodes(Graph, 4000)
    return Graph

def genNullModelBfs(h):
    totalNodes = 3*(2**h)-2
    x_cords, y_cords = [0], [0]
    for i in range(h+1): #negative loop
        frac = 2.0**(i)/totalNodes
        reach = i+1
        x_cords.append(x_cords[-1])
        y_cords.append(reach)
        x_cords.append(frac+x_cords[-1])
        y_cords.append(reach)
    
    for i in range(1, h+1): #pos loop
        frac = 2.0**(h-i)/totalNodes
        reach = 3*(2**i)+h-i-2
        x_cords.append(x_cords[-1])
        y_cords.append(reach)
        x_cords.append(frac+x_cords[-1])
        y_cords.append(reach)
        
    return (x_cords, y_cords)

def generatePlots(G, name):
    inList, inList_GNP, outList, outList_GNP = [], [], [], []
    GNP = genErdosRenyi(G.GetNodes(), G.GetEdges())
    numSamples = 5000
    for i in tqdm(range(numSamples)):
        nodeId = G.GetRndNId()
        nodeGNPId = GNP.GetRndNId()
        outList.append(snap.GetBfsTree(G, nodeId, True, False).GetNodes())
        outList_GNP.append(snap.GetBfsTree(GNP, nodeGNPId, True, False).GetNodes())
        inList.append(snap.GetBfsTree(G, nodeId, False, True).GetNodes())
        inList_GNP.append(snap.GetBfsTree(GNP, nodeGNPId, False, True).GetNodes())
    inList.sort()
    outList.sort()
    inList_GNP.sort()
    outList_GNP.sort()

    x = [float(i)/numSamples for i in range(numSamples)]
    plt.subplot(2, 1, 1)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.plot(x, inList, label=name)
    plt.plot(x, inList_GNP, label='Erdos Reini')
    x_cords13, y_cords13 = genNullModelBfs(13)
    plt.plot(x_cords13, y_cords13, label='Null Model')
    plt.title('Reachability using Inlinks - {}'.format(name))
    plt.xlabel('Frac. of starting node (%)')
    plt.ylabel('number of nodes reached')
    plt.legend()

    plt.subplot(2, 1, 2)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.plot(x, outList, label=name)
    plt.plot(x, outList_GNP, label='Erdos Reini')
    x_cords12, y_cords12 = genNullModelBfs(12)
    plt.plot(x_cords12, y_cords12, label='Null Model')
    plt.title('Reachability using outlinks - {}'.format(name))
    plt.xlabel('Frac. of starting node (%)')
    plt.ylabel('number of nodes reached')
    plt.legend()
    plt.tight_layout()
    plt.savefig(name+'_BFS.pdf')
    plt.clf()



G_LS174t = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist_v2/LS174T_clean_EdgesList.txt", 0, 1)
generatePlots(G_LS174t, 'LS174T')
# print "Efficiency of LS174T: {}".format( getEfficiency(G_LS174t))

G_SW122 = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist_v2/SW1222_clean_EdgesList.txt", 0, 1)
generatePlots(G_SW122, 'SW1222')
# print "Efficiency of SW1222: {}".format( getEfficiency(G_SW122))