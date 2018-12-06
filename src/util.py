import matplotlib.pyplot as plt
import snap
import copy
import util_data
import collections
from scipy.stats import pearsonr

#tree+inverted tree nullmodel
def genTreeInvTreeNullModelPlot(h):
    """
     Generate a tree-inverted tree null model with height of tree/inverted tree is h

     /  
    /\ 
    \/ 
     \ 
    """
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

def genTreeInvTreeNullModel(h):
    G = snap.TNGraph.New()
    [G.AddNode(i) for i in range(1, 2**(h))]
    [G.AddEdge(i, 2*i) for i in range(1, 2**(h-1))]
    [G.AddEdge(i, 2*i+1) for i in range(1, 2**(h-1))]

    lastNode =  2**(h)-1

    [G.AddNode(lastNode+i) for i in range(1, 2**(h+1))]
    [G.AddEdge(lastNode+2*i, lastNode+i) for i in range(1, 2**(h))]
    [G.AddEdge(lastNode+2*i+1, lastNode+i) for i in range(1, 2**(h))]
    
    for i in range(2**(h-1), 2**h):
        G.AddEdge(i, 2*i+lastNode)
        G.AddEdge(i, 2*i+1+ lastNode)
    return G


def getNodeSplit(G, nodes=[]):
    """
    Finds all the nodes where arteries split
    return: (divergent, convergent) - number of nodes in each

              out
              ^
              |
              |
    in --------> ------> out
    """

    divNodes = 0
    convNodes = 0

    if len(nodes) == 0:
        for node in G.Nodes():
            if node.GetInDeg() > node.GetOutDeg():
                convNodes += 1
            if node.GetInDeg() < node.GetOutDeg():
                divNodes += 1 
    else:
        for node in nodes:
            node = G.GetNI(node)
            if node.GetInDeg() > node.GetOutDeg():
                convNodes += 1
            if node.GetInDeg() < node.GetOutDeg():
                divNodes += 1 
    return (divNodes, convNodes)




def getInletIds(G):
    return [node.GetId() for node in G.Nodes() if node.GetInDeg() == 0]

def getOutletIds(G):
    return [node.GetId() for node in G.Nodes() if node.GetOutDeg() == 0]


def BFS_mod(G, nodeIds):
    """
    Outlink following BFS with printing every time the queue is empty
    G - tumour graph
    nodeId - node to start from
    return (outBFS, propFrontDivLen, propFrontConvLen, propFrontLens)
    """
    # node = G.GetNI(nodeId)
    print "Inlets: {}".format(nodeIds)
    queue = set(nodeIds)
    propFront = set()
    propFrontDivLen, propFrontConvLen, propFrontLens = [], [], []
    outBFS = copy.deepcopy(queue)
    seen = set()
    while len(queue) != 0 or len(propFront) != 0:
        if len(queue) == 0:
            queue = queue.union(propFront)
            outBFS = outBFS.union(propFront)
            divNodes, convNodes = getNodeSplit(G, propFront)
            propFrontDivLen.append(divNodes)              
            propFrontConvLen.append(convNodes)               
            propFrontLens.append(len(propFront))
            print "\ntotal nodes: {}".format(len(propFront))
            print "divergent nodes: {}".format(divNodes)
            print "convergent nodes: {}".format(convNodes)
            propFront.clear()
        newNode = queue.pop()
        seen.add(newNode)
        propFront = propFront.union(set([G.GetNI(newNode).GetOutNId(i) for i in range(G.GetNI(newNode).GetOutDeg()) if G.GetNI(newNode).GetOutNId(i) not in seen]))
    return (outBFS, propFrontDivLen, propFrontConvLen, propFrontLens)
        

def plotBFSProp(Graph, name):
    inlets = getInletIds(Graph)
    outBFS, propFrontDivLen, propFrontConvLen, propFrontLens = BFS_mod(Graph, inlets)
    plt.plot(propFrontConvLen, label='Convergent Nodes')
    plt.plot(propFrontDivLen, label='Diveregent Nodes')
    # plt.plot(propFrontLens, label='Total Nodes')
    plt.title(name)
    plt.xlabel('BFS iteration')
    plt.ylabel('Number of Nodes')
    plt.legend()
    plt.savefig('Propagation_'+name+'.pdf')
    plt.clf()
    # plt.show()

def getInletHist(G, name):
    """
    Run an outward BFS from each inlet and plot the size of the outset
    """

    inlets = getInletIds(G)
    BFS_out = []  
    for inlet in inlets:
        BFS_out.append(snap.GetBfsTree(G, inlet, True, False).GetNodes())
    # print len(BFS_out)
    # print sorted(BFS_out, key=lambda args: args[1])
    plt.hist(BFS_out)
    plt.xlabel('Size of Outset')
    plt.ylabel('Number of Nodes')
    plt.savefig("Inlet_BFS_"+name+"_Hist.pdf")
    plt.clf()

def getInletScatter(G, name):
    """
    Run an outward BFS from each inlet and plot the size of the outset
    """

    inlets = getInletIds(G)
    BFS_out = []  
    for inlet in inlets:
        BFS_out.append(snap.GetBfsTree(G, inlet, True, False).GetNodes())
    # print len(BFS_out)
    # print sorted(BFS_out, key=lambda args: args[1])
    plt.scatter(*zip(*collections.Counter(BFS_out).items()))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Size of Outset')
    plt.ylabel('Number of Nodes')
    plt.savefig("Inlet_BFS_"+name+"_Scatter.pdf")
    plt.clf()

def getInletRadiusPlot(G, name, filename):
    """
    Run an outward BFS from each inlet and plot the size of the outset vs node "Radius"
    """

    inlets = getInletIds(G)
    BFS_out, x_cords = [], []
    BFS_mean_rad = {}
    if name == 'Healthy':
        radii = util_data.loadNodeAttr(filename, small=True)
    else:
        radii = util_data.loadNodeAttr(filename)

    for inlet in inlets:
        BFS_out.append((snap.GetBfsTree(G, inlet, True, False).GetNodes(), radii[inlet]))
    
    for elem in BFS_out:
        if elem[0] in BFS_mean_rad:
            BFS_mean_rad[elem[0]].append(elem[1])
        else: 
            BFS_mean_rad[elem[0]] = [elem[1]]

    for elem in BFS_mean_rad:
        BFS_mean_rad[elem] = sum(BFS_mean_rad[elem])/len(BFS_mean_rad[elem])
    
    # y_cords, x_cords = zip(*BFS_out)
    print "Pearson Coeff for {}: {}".format(name, pearsonr(BFS_mean_rad.values(), BFS_mean_rad.keys()))
    plt.plot(BFS_mean_rad.values(), BFS_mean_rad.keys(), 'o')
    plt.xlabel('Mean Radius $(\mu m)$', usetex=True)
    plt.ylabel('Size of Outset')
    plt.savefig("Inlet_BFS_Radius_"+name+".pdf")
    plt.clf()
