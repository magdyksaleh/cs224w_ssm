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

def edge_dir(dat, numedgepts_idx=3, pres_idx=6):
    dirs = []
    i=0
    for x in dat[numedgepts_idx]:
        if dat[pres_idx][i]>dat[pres_idx][i+int(x)-1]:
            dirs.append(1)
        elif dat[pres_idx][i]<dat[pres_idx][i+int(x)-1]:
            dirs.append(-1)
        else:
            dirs.append(0)
        i+=int(x[0])  
    return dirs

def clean_edgelist(dat, edge_dirs, edgelist_idx=2):
    clean = []
    seen = set()
    for i,x in enumerate(dat[edgelist_idx]):
        if edge_dirs[i]!=0 and x[0]!=x[1]:
            if edge_dirs[i]==1:
                e = (int(x[0]),int(x[1]))
            else:
                e = (int(x[1]), int(x[0]))
            if e not in seen:
                clean.append(e)
                seen.add(e)
    return clean
    
LS = read_am('data/og_files/LS174T_spatialGraph_RIN.txt')
LS_edge_dir = edge_dir(LS)
assert(len(LS_edge_dir)==len(LS[2]))
LS_clean = clean_edgelist(LS, LS_edge_dir)

# write edgelist to file
# f = open('data/Edgelist_v2/LS174T_clean_EdgesList.txt', 'w')
# for x in LS_clean:
#     f.write('%d\t%d\n'%(x[0],x[1]))
# f.close()