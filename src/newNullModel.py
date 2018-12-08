import snap
from util import *
import matplotlib.pyplot as plt
import numpy as np
import util_data




G_LS174t = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/LS174T_clean_EdgesList.txt", 0, 1)
G_SW1222 = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/SW1222_clean_EdgesList.txt", 0, 1)
G_Mesentery = snap.LoadEdgeList(snap.PNEANet, "../data/Edgelist/Mesentery_clean_EdgeList.txt", 0, 1)
div, conv = getNodeSplit(G_LS174t)


G_NullTit = genTreeInvTreeNullModel(13)
# for edge in G_NullTit.Edges():
#     print edge.GetSrcNId(), edge.GetDstNId()
# plotBFSProp(G_NullTit, 'Null_h_3') 

print(snap.GetMxWcc(G_Mesentery).GetNodes())
getInletHist(G_Mesentery, 'Healty')
getInletScatter(G_Mesentery, 'Healthy')
getInletRadiusPlot(G_Mesentery, 'Healthy',
    '/Users/magdy/Desktop/Stanford/Fall18/224w/project/data/og_files/Flow2Amira(Small).txt')



print(snap.GetMxWcc(G_LS174t).GetNodes())
getInletHist(G_LS174t, 'LS174T')
getInletScatter(G_LS174t, 'LS174T')
getInletRadiusPlot(G_LS174t, 'LS174T',
    '/Users/magdy/Desktop/Stanford/Fall18/224w/project/data/og_files/spatialGraph_RIN.am')


print(snap.GetMxWcc(G_SW1222).GetNodes())
getInletHist(G_SW1222, 'SW1222')
getInletScatter(G_SW1222, 'SW1222') 

getInletRadiusPlot(G_SW1222, 'SW1222',
    '/Users/magdy/Desktop/Stanford/Fall18/224w/project/data/og_files/SW122_spatialGraph_RIN.txt')