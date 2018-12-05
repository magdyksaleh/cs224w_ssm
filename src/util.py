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

              out
              ^
              |
              |
    in --------> ------> out
    """
    divNodes = 0
    convNodes = 0
    for node in G.Nodes():
        if node.GetInDeg() > node.GetOutDeg():
            convNodes += 1
        if node.GetInDeg() < node.GetOutDeg():
            divNodes += 1 
    return (divNodes, convNodes)


