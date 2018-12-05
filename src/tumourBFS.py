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
    inList, outList = [], []
    numSamples = 5000
    for i in tqdm(range(numSamples)):
        nodeId = G.GetRndNId()
        outList.append(snap.GetBfsTree(G, nodeId, True, False).GetNodes())
        inList.append(snap.GetBfsTree(G, nodeId, False, True).GetNodes())
    inList.sort()
    outList.sort()

    x = [float(i)/numSamples for i in range(numSamples)]
    plt.subplot(2, 1, 1)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.plot(x, inList, label='LS174T')
    x_cords13, y_cords13 = genNullModelBfs(13)
    plt.plot(x_cords13, y_cords13, label='Null Model')
    plt.title('Reachability using Inlinks - {}'.format(name))
    plt.xlabel('Frac. of starting node (%)')
    plt.ylabel('number of nodes reached')
    plt.legend()

    plt.subplot(2, 1, 2)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.plot(x, outList, label='LS174T')
    x_cords12, y_cords12 = genNullModelBfs(12)
    plt.plot(x_cords12, y_cords12, label='Null Model')
    plt.title('Reachability using outlinks - {}'.format(name))
    plt.xlabel('Frac. of starting node (%)')
    plt.ylabel('number of nodes reached')
    plt.legend()
    plt.tight_layout()
    plt.savefig(name+'_BFS.pdf')
    plt.show()
    plt.clf()



G_LS174t = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/LS174T_clean_EdgesList.txt", 0, 1)
generatePlots(G_LS174t, 'LS174T')
# print "Efficiency of LS174T: {}".format( getEfficiency(G_LS174t))

G_SW122 = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/SW1222_clean_EdgesList.txt", 0, 1)
# generatePlots(G_LS174t, 'SW1222')
print "Efficiency of SW1222: {}".format( getEfficiency(G_SW122))