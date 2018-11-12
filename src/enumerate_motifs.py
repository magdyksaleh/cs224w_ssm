import snap
import sys, os, glob
import itertools
import networkx as nx
from matplotlib import pyplot as plt
import tqdm
sys.path.append('../')

class MotifCounter(object):
    def __init__(self, num_nodes, subgraph_path = "../data/subgraphs"):
        self.node_motifs = {2:2, 3:13, 4:199}
        self.subgraph_path = os.path.join(subgraph_path, str(num_nodes))
        if not os.path.exists(self.subgraph_path):
            os.makedirs(self.subgraph_path)
        self.num_nodes = num_nodes
        if len(glob.glob(os.path.join(self.subgraph_path, '*.txt')))==0:
            self.motifs = self.create_motifs()
            [snap.SaveEdgeList(graph, os.path.join(self.subgraph_path,"{}.txt".format(i))) for i,graph in enumerate(self.motifs)]
        else:
            self.motifs = [snap.LoadEdgeList(snap.PNGraph, os.path.join(self.subgraph_path,"{}.txt".format(i)), 0, 1) for i in range(self.node_motifs[num_nodes])]
    def match(self, G1, G2):
        if G1.GetEdges() > G2.GetEdges():
            return False
        else:
            G = G2
            H = G1
            
        for p in itertools.permutations(range(self.num_nodes)):
            edge = G.BegEI()
            matches = True
            while edge < G.EndEI():
                if not H.IsEdge(p[edge.GetSrcNId()], p[edge.GetDstNId()]):
                    matches = False
                    break
                edge.Next()
            if matches:
                break
        return matches
    
    def create_motifs(self, draw=False):
        #print("hello")
        graph_list = []
        num_found = 0
        for graph in tqdm.tqdm(self.enumerate_graphs(self.num_nodes), total=2**(self.num_nodes*(self.num_nodes-1))):            
            if num_found==0:
                graph_list.append(graph)
                if draw:
                    self.draw_graph(graph, os.path.join(self.subgraph_path, str(num_found)+'.png') )
                num_found += 1
            else:
                for generated_subgraph in graph_list:
                    isomorphic = self.match(graph, generated_subgraph)
                    if isomorphic:
                        break
                if not isomorphic:
                    if draw:
                        self.draw_graph(graph, os.path.join(self.subgraph_path, str(num_found)+'.png') )
                    graph_list.append(graph)
                    num_found+=1
        return graph_list
    
    def enumerate_graphs(self, k):
        for seq in itertools.product("01", repeat=k*(k-1)):
            g = snap.TNGraph.New()
            for i in range(k): g.AddNode(i)
            for i,e in enumerate(seq):
                if e=='1':
                    start_node = i/(k-1)
                    end_node = i % (k-1)
                    if end_node >= start_node:
                        end_node += 1
                    g.AddEdge(start_node, end_node)
            if snap.GetMxWcc(g).GetNodes()==k:
                yield g 
            
    def draw_graph(self, g, outname):
        plt.clf()
        G = nx.DiGraph()
        for n in g.Nodes():
            G.add_node(n.GetId())
        for e in g.Edges():
            G.add_edge(e.GetSrcNId(), e.GetDstNId())
        nx.draw(G)
        plt.savefig(outname)

    def count_motifs(self, G):
        
        self.counts = [0]*self.node_motifs[num_nodes]
        for node in tqdm(G.Nodes(), total = G.GetNodes()):
            v = node.GetId()
            v_extension = set([nbr for nbr in node.GetOutEdges() if nbr > v])
            v_extension.update([nbr for nbr in node.GetInEdges() if nbr > v])
            self.extend_subgraph(G, self.num_nodes, [v], v_extension, v, verbose)
        return self.counts

    def extend_subgraph(self, G, k, sg, v_ext, node_id, verbose=False):

        if len(sg) is k:
            count_iso(G, sg, verbose)
            return
        sg_nbrs = copy.deepcopy(v_ext)

        while len(v_ext) != 0:
            w = v_ext.pop()
            w_nodeI = G.GetNI(w)
            v_new_ext = copy.deepcopy(v_ext)
            sg.append(w)
            v_new_ext.update([nbr for nbr in w_nodeI.GetOutEdges() if (nbr > node_id and nbr not in sg and nbr not in sg_nbrs)])
            v_new_ext.update([nbr for nbr in w_nodeI.GetInEdges() if (nbr > node_id and nbr not in sg and nbr not in sg_nbrs)])
            extend_subgraph(G,k, sg, v_new_ext, node_id, verbose)
            sg.remove(w)
    
    def count_iso(self, G, sg, verbose=False):
    if verbose:
        print(sg)
    nodes = snap.TIntV()
    for NId in sg:
        nodes.Add(NId)
    # This call requires latest version of snap (4.1.0)
    SG = snap.GetSubGraphRenumber(G, nodes)
    for i in range(len(self.motifs)):
        if self.match(self.motifs[i], SG):
            self.counts[i] += 1