import snap
from util import *
import matplotlib.pyplot as plt
import numpy as np


G_LS174t = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/LS174T_clean_EdgesList.txt", 0, 1)
G_SW122 = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/SW1222_clean_EdgesList.txt", 0, 1)

G_Mesentery = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/Mesentery_clean_EdgeList.txt", 0, 1)


div, conv = getNodeSplit(G_LS174t)

print "converegent frac: {}".format(float(conv)/(div+conv))
print "diveregent frac: {}".format(float(div)/(div+conv))

print "Total Nodes: {}".format(G_LS174t.GetNodes())
print "Nodes w/ equal in and out: {}".format(G_LS174t.GetNodes()-div-conv)

inlets = getInletIds(G_Mesentery)
print inlets[1]

outBFS, propFrontDivLen, propFrontConvLen, propFrontLens = BFS_mod(G_Mesentery, inlets)

print len(outBFS)

plt.plot(propFrontConvLen, label='Convergent Nodes')
plt.plot(propFrontDivLen, label='Diveregent Nodes')
# plt.plot(propFrontLens, label='Total Nodes')
plt.xlabel('BFS iteration')
plt.ylabel('Number of Nodes')
plt.legend()
plt.show()
