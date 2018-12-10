import matplotlib.pyplot as plt
import snap
import copy
import util_data
import collections
import bisect
import numpy as np
from collections import OrderedDict
from scipy.stats import pearsonr
from tqdm import tqdm
from data_cleaning.snair.clean_am import *

#tree+inverted tree nullmodel
def genTreeInvTreeNullModelPlot(h):
    """
     Generate a tree-inverted tree null model with height of tree/inverted tree is h

     | 
    /\ 
    \/ 
     | 
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
    #print "Inlets: {}".format(nodeIds)
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
            #print "\ntotal nodes: {}".format(len(propFront))
            #print "divergent nodes: {}".format(divNodes)
            #print "convergent nodes: {}".format(convNodes)
            propFront.clear()
        newNode = queue.pop()
        seen.add(newNode)
        propFront = propFront.union(set([i for i in G.GetNI(newNode).GetOutEdges() if i not in seen]))
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

def influence_maximisation(G, possible_nodes, k=5):
    # Take in graph G, nodes from which to select top k, and returns optimal set of k
    # nodes to maximise influence as well as the set of influenced nodes
    possible_nodes = list(possible_nodes)
    influence_dict = {node:get_influence_set(G, node) for node in possible_nodes}
    influence_dict = OrderedDict(sorted(influence_dict.items(), key=lambda t: len(t[1]), reverse = True))
    optimal_set = set()
    current_influence = set()
    list_of_lengths = []
    for i in tqdm(range(k)):
        # print "Beginning selection %d"%i
        top_node, influence_set = influence_dict.popitem(False)
        optimal_set.add(top_node)
        current_influence = current_influence | influence_set
        if i==k-1:
            break
        influence_list = influence_dict.items()
        calculated_set = set()
        # pbar = tqdm(total=G.GetNodes())
        while True:
            if influence_list[0][0] in calculated_set:
                break
            checked_node, checked_influence = influence_list[0]
            checked_influence = get_influence_set(G, checked_node, current_influence)
            calculated_set.add(checked_node)
            _, influences = zip(*influence_list)
            influences = [len(inf_set) for inf_set in influences]
            insert_idx = bisect.bisect_left(influences[::-1][:-1], len(checked_influence))
            if insert_idx == len(influences):
                break
            else:
                del influence_list[0]
                influence_list.insert(len(influence_list)-insert_idx, (checked_node, checked_influence))
        #     pbar.update(1)
        # pbar.close()
        influence_dict = OrderedDict(influence_list)
        list_of_lengths.append(len(current_influence))
    return optimal_set, current_influence, list_of_lengths


def get_influence_set(G, x, S = set([])):
    #S is set of already influenced nodes
    x_influence = set([node.GetId() for node in snap.GetBfsTree(G, x, True, False).Nodes()])
    return x_influence - S


def test_influence_maximisation():
    G = snap.TNGraph.New()
    for i in range(7):
        G.AddNode(i)
    G.AddEdge(0,1)
    G.AddEdge(1,2)
    G.AddEdge(2,3)
    G.AddEdge(4,5)
    G.AddEdge(5,6)
    possible_nodes = [node.GetId() for node in G.Nodes()]
    assert influence_maximisation(G, possible_nodes, k=2) == set([0,4])


def get_basic_feature(graph, node_id):
    NI = graph.GetNI(node_id)
    v_1 = NI.GetDeg()
    nbrs = NI.GetOutEdges()
    nbr_vec = snap.TIntV()
    nbrs = [nbr_vec.Add(nbr) for nbr in nbrs]
    nbr_vec.Add(node_id)
    subgraph = snap.GetSubGraph(graph, nbr_vec)
    v_2 = subgraph.GetEdges()
    total_edges = 0
    for node in subgraph.Nodes():
        orig_NI = graph.GetNI(node.GetId())
        total_edges += orig_NI.GetDeg()
    v_3 = total_edges - 2*v_2
    feature = snap.TFltV()
    feature.Add(v_1)
    feature.Add(v_2)
    feature.Add(v_3)

    return feature

def calc_similarity(feature_1, feature_2):
    norm_1 = np.linalg.norm(np.asarray(feature_1))
    norm_2 = np.linalg.norm(np.asarray(feature_2))
    if norm_1 == 0 or norm_2 == 0:
        return 0
    return np.sum(np.asarray(feature_1) * np.asarray(feature_2))/(norm_1 * norm_2)



def run_rolx_iteration(graph, features):
    new_features = snap.TIntFltVH()
    iterator =  features.BegI()
    while not iterator.IsEnd():
        node_id = iterator.GetKey()
        nodeI = graph.GetNI(node_id)
        node_nbrs = [nbr for nbr in nodeI.GetOutEdges()]
        sum_aggregate = np.zeros(features[0].Len())
        num_nbrs = len(node_nbrs)
        for node_nbr in node_nbrs:
            sum_aggregate += np.asarray(features[node_nbr], dtype=np.float32)
        new_feature = snap.TFltV(features[node_id])
        if num_nbrs == 0:
            mean_aggregate = np.zeros(features[0].Len())
            sum_aggregate = np.zeros(features[0].Len())
        else:
            mean_aggregate = sum_aggregate/len(node_nbrs)
        for i in range(mean_aggregate.size):
            new_feature.Add(mean_aggregate[i])
        for i in range(sum_aggregate.size):
            new_feature.Add(sum_aggregate[i])
        new_features[node_id] = new_feature
        iterator.Next()
    return new_features

def get_pos(txt_file, small = False):
    if small:
        coordinates = np.loadtxt(txt_file)[:,:2]
        return coordinates
    LS = read_am(txt_file)
    LS_edge_prop = edgeProps(LS, small = small)
#     assert(len(LS_edge_prop)==len(LS[2]))
    LS_clean = cleanEdgelist(LS, LS_edge_prop)
    coordinates = {}
    for edge, attr in LS_clean.items():
        coordinates[edge[0]] = np.asarray(attr['srcCords'])[:2]
        coordinates[edge[1]] = np.asarray(attr['dstCords'])[:2]
#     return dict(zip(range(coordinates.shape[0]), coordinates))
    return coordinates
