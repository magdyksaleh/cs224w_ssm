import snap
import numpy as np


def loadNodeAttr(filename, small=False):
    """ 
    Loads the node attributes from the amira files
    """
    #get number of edges points
    f = open(filename, 'r')
    edgeFlag, edgePointFlag, radiusFlag = False, False, False
    edges, radii = [], []
    edgePointDict, nodeDict = {}, {}
    edgePointCntr = 0

    for line in f:
        #grab edges
        if len(line.strip()) == 0: continue
        if edgeFlag and ("@3" in line): 
            edgeFlag = False
            edgePointFlag = True
            continue
        
        if edgePointFlag and ("@4" in line): 
            edgePointFlag = False
            continue
        
        if ("@5" in line) and not ("{" in line): 
            radiusFlag = True
            continue
        
        if ("@6" in line) and not ("{" in line): 
            radiusFlag = False
            break

        if radiusFlag:
            radii.append(float(line.strip()))
        
        if edgeFlag:
            [srd, dst] = line.split()
            edges.append((int(srd), int(dst)))
        if ("@2" in line) and not ("{" in line): 
            edgeFlag = True

        
        #grab edgepoints
        if edgePointFlag and not small:
            edgePointDict[edges[edgePointCntr]] = int(line.strip())
            edgePointCntr += 1
        

    cntr = 0
    for edge in edges:
        numPoints = edgePointDict[edge] if not small else 2
        src, dst = edge

        if src in nodeDict:
            nodeDict[src].append(radii[cntr])
        else: nodeDict[src] = [radii[cntr]]
        
        cntr += numPoints - 1
        if dst in nodeDict:
            nodeDict[dst].append(radii[cntr])
        else: nodeDict[dst] = [radii[cntr]]
        cntr += 1
    
    for node in nodeDict:
        nodeDict[node] = sum(nodeDict[node])/len(nodeDict[node])

    return nodeDict