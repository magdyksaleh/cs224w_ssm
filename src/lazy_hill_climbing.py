import snap
import numpy as np
import bisect
from collections import OrderedDict

count = 0
def influence_maximisation(G, possible_nodes, k=5):
	possible_nodes = list(possible_nodes)
	influence_dict = {node:get_influence_set(G, node) for node in possible_nodes}
	influence_dict = OrderedDict(sorted(influence_dict.items(), key=lambda t: len(t[1]), reverse = True))
	optimal_set = set()
	current_influence = set()
	for i in range(k):
		top_node, influence_set = influence_dict.popitem(False)
		optimal_set.add(top_node)
		current_influence = current_influence | influence_set
		if i==k-1:
			break
		influence_list = influence_dict.items()
		calculated_set = set()
		while True:
			if influence_list[0][0] in calculated_set:
				break
			checked_node, checked_influence = influence_list[0]
			checked_influence = get_influence_set(G, checked_node, current_influence)
			calculated_set.add(checked_node)
			_, influences = zip(*influence_list)
			influences = [len(inf_set) for inf_set in influences]
			insert_idx = bisect.bisect_left(influences[::-1][:-1], len(checked_influence))
			if insert_idx == len(influences):
				break
			else:
				del influence_list[0]
				influence_list.insert(len(influence_list)-insert_idx, (checked_node, checked_influence))
		influence_dict = OrderedDict(influence_list)
	return optimal_set


def get_influence_set(G, x, S = set([])):
	#S is set of already influenced nodes
	x_influence = set([node.GetId() for node in snap.GetBfsTree(G, x, True, False).Nodes()])
	return x_influence - S

if __name__ =='__main__':
	G = snap.TNGraph.New()
	for i in range(7):
		G.AddNode(i)
	G.AddEdge(0,1)
	G.AddEdge(1,2)
	G.AddEdge(2,3)
	G.AddEdge(4,5)
	G.AddEdge(5,6)
	possible_nodes = [node.GetId() for node in G.Nodes()]
	print influence_maximisation(G, possible_nodes, k=2)
	print count