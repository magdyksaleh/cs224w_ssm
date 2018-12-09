from collections import defaultdict
import copy

def read_am(fname):
    f = open(fname)
    d = [x.strip() for x in f]
    f.close()
    dat = {}
    d = d[d.index('@1'):]
    cur = -1
    for x in d:
        if x:
            if x[0]=='@':
                cur = int(x[1:])
                dat[cur] = []
            else:
                dat[cur].append([float(y) for y in x.split(' ')])
    return dat

def edgeProps(dat, cord_idx = 1, edge_idx = 2, numedgepts_idx = 3, pres_idx = 6, flow_idx = 7, radii_idx = 5):
    edgeProps            = defaultdict(lambda: defaultdict(list))  #edge indexed dict
    edgePointCntr      = 0 

    for i, x in enumerate(dat[numedgepts_idx]):
        curEdge = tuple(dat[edge_idx][i])
        x = x[0]

        srcCords, dstCords  = dat[cord_idx][int(curEdge[0])], dat[cord_idx][int(curEdge[1])]
        srcPreassure, dstPreassure  = dat[pres_idx][edgePointCntr ][0], dat[pres_idx][edgePointCntr+int(x)-1 ][0]
        srcFlow, dstFlow            = dat[flow_idx][edgePointCntr ][0], dat[flow_idx][edgePointCntr+int(x)-1 ][0]
        srcRadius, dstRadius        = dat[radii_idx][edgePointCntr][0], dat[radii_idx][edgePointCntr+int(x)-1][0]        
        
        edgeProps[curEdge]['srcCords'     ].append(srcCords       )
        edgeProps[curEdge]['dstCords'     ].append(dstCords       )
        edgeProps[curEdge]['srcPreassure' ].append(srcPreassure   )
        edgeProps[curEdge]['dstPreassure' ].append(dstPreassure   )
        edgeProps[curEdge]['srcFlow'      ].append(srcFlow        )
        edgeProps[curEdge]['dstFlow'      ].append(dstFlow        )
        edgeProps[curEdge]['srcRadius'    ].append(srcRadius      )
        edgeProps[curEdge]['dstRadius'    ].append(dstRadius      )
        


        
        if srcPreassure > dstPreassure:
            edgeProps[curEdge]['pIdx'] = 1
        elif  srcPreassure < dstPreassure:
            edgeProps[curEdge]['pIdx'] = -1
        else:
            edgeProps[curEdge]['pIdx'] = 0
        edgePointCntr += int(x)
    return edgeProps

def cleanEdgelist(dat, edgeProps, edgelist_idx=2):
    for x in copy.deepcopy(edgeProps):
        curEdge = tuple(x)
        if x[0] == x[1]: #delete self edges
            del edgeProps[curEdge]
            continue
        
        #average values
        edgeProps[curEdge]['srcCords'     ] = edgeProps[curEdge]['srcCords'     ][0]
        edgeProps[curEdge]['dstCords'     ] = edgeProps[curEdge]['dstCords'     ][0]
        edgeProps[curEdge]['srcPreassure' ] = edgeProps[curEdge]['srcPreassure' ][0]
        edgeProps[curEdge]['dstPreassure' ] = edgeProps[curEdge]['dstPreassure' ][0]
        edgeProps[curEdge]['srcFlow'      ] = sum(edgeProps[curEdge]['srcFlow'      ])
        edgeProps[curEdge]['dstFlow'      ] = sum(edgeProps[curEdge]['dstFlow'      ])
        edgeProps[curEdge]['srcRadius'    ] = float(sum(edgeProps[curEdge]['srcRadius'    ]))/len(edgeProps[curEdge]['srcRadius'    ])
        edgeProps[curEdge]['dstRadius'    ] = float(sum(edgeProps[curEdge]['dstRadius'    ]))/len(edgeProps[curEdge]['dstRadius'    ])

        if edgeProps[curEdge]['pIdx'] == -1: #swap values
            print "old: ", edgeProps[curEdge]
            curEdgeProps = edgeProps[curEdge]
            srcCords, dstCords          = curEdgeProps['srcCords'],     curEdgeProps['dstCords']
            srcPreassure, dstPreassure  = curEdgeProps['srcPreassure'], curEdgeProps['dstPreassure']
            srcFlow, dstFlow            = curEdgeProps['srcFlow'],      curEdgeProps['dstFlow']
            srcRadius, dstRadius        = curEdgeProps['srcRadius'],    curEdgeProps['dstRadius']
            
            newEdge = (x[1], x[0])

            edgeProps[newEdge]['srcCords'     ] = (dstCords       )
            edgeProps[newEdge]['dstCords'     ] = (srcCords       )
            edgeProps[newEdge]['srcPreassure' ] = (dstPreassure   )
            edgeProps[newEdge]['dstPreassure' ] = (srcPreassure   )
            edgeProps[newEdge]['srcFlow'      ] = (dstFlow        )
            edgeProps[newEdge]['dstFlow'      ] = (srcFlow        )
            edgeProps[newEdge]['srcRadius'    ] = (dstRadius      )
            edgeProps[newEdge]['dstRadius'    ] = (srcRadius      )
            print "new: ", edgeProps[newEdge]
            del edgeProps[curEdge]

    return edgeProps
    
LS = read_am('../../../data/og_files/LS174T_spatialGraph_RIN.txt')
LS_edge_prop = edgeProps(LS)
# assert(len(LS_edge_prop)==len(LS[2]))
LS_clean = cleanEdgelist(LS, LS_edge_prop)

# write edgelist to file
# f = open('data/Edgelist_v2/LS174T_clean_EdgesList.txt', 'w')
# for x in LS_clean:
#     f.write('%d\t%d\n'%(x[0],x[1]))
# f.close()