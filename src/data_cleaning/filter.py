import numpy as np
import scipy.io as sio
import sys

LS174T_data = sio.loadmat('SW1222_data.mat')
edges = []
for elem in LS174T_data['seg_filtered'][0]:
    edgeRef = elem[0][0]
    [edges.append(tuple(edgeRef[:-1])) for i in range(edgeRef[2]+1)]

outFile = open("../../data/SW1222__clean_EdgesList.txt", "w")
for edge in edges:
    outFile.write(str(edge[0]) + "\t" + str(edge[1])+ "\n")
outFile.close()

    
