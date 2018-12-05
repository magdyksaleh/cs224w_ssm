import matplotlib.pyplot as plt
import snap

#tree+inverted tree nullmodel
def genTreeInvTreeNullModel(h):
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

def getNodeSplit(G):
    """
    Finds all the nodes where arteries split
    return: (divergent, convergent)

              out
              ^
              |
              |
    in --------> ------> out
    """



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
    queue = [
        G.GetNI(nodeId).GetOutNId(i) for nodeId in nodeIds for i in range(G.GetNI(nodeId).GetOutDeg())
    ]
    propFront = []
    propFrontDivLen = []
    propFrontConvLen = []
    propFrontLens = []
    outBFS = queue[:]
    while len(queue) != 0 or len(propFront) != 0:
        if len(queue) == 0:
            queue += propFront
            outBFS += propFront
            divNodes, convNodes = getNodeSplit(G, propFront)
            propFrontDivLen.append(divNodes)              
            propFrontConvLen.append(convNodes)               
            propFrontLens.append(len(propFront))
            print "total nodes: {}".format(propFront[-1])
            print "divergent nodes: {}".format(propFrontDivLen[-1])
            print "convergent nodes: {}".format(propFrontConvLen[-1])
            propFront = []
        newNode = queue.pop()
        propFront += [G.GetNI(newNode).GetOutNId(i) for i in range(G.GetNI(newNode).GetOutDeg())]
    return (outBFS, propFrontDivLen, propFrontConvLen, propFrontLens)
        


