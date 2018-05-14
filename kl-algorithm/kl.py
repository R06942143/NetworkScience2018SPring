from __future__ import division
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import islice
from operator import itemgetter
import random
import numpy as np
import networkx as nx
from networkx.utils import not_implemented_for
from networkx.algorithms.community.community_utils import is_partition
import time
__all__ = ['kernighan_lin_bisection']


def _compute_delta(G, A, B, weight):
    # helper to compute initial swap deltas for a pass
    delta = defaultdict(float)
    for u, v, d in G.edges(data=True):
        w = d.get(weight, 1)
        if u in A:
            if v in A:
                delta[u] -= w
                delta[v] -= w
            elif v in B:
                delta[u] += w
                delta[v] += w
        elif u in B:
            if v in A:
                delta[u] += w
                delta[v] += w
            elif v in B:
                delta[u] -= w
                delta[v] -= w
    return delta


def _update_delta(delta, G, A, B, u, v, weight):
    # helper to update swap deltas during single pass
    for _, nbr, d in G.edges(u, data=True):
        w = d.get(weight, 1)
        if nbr in A:
            delta[nbr] += 2 * w
        if nbr in B:
            delta[nbr] -= 2 * w
    for _, nbr, d in G.edges(v, data=True):
        w = d.get(weight, 1)
        if nbr in A:
            delta[nbr] -= 2 * w
        if nbr in B:
            delta[nbr] += 2 * w
    return delta


def _kernighan_lin_pass(G, A, B, weight):
    # do a single iteration of Kernighan–Lin algorithm
    # returns list of  (g_i,u_i,v_i) for i node pairs u_i,v_i
    multigraph = G.is_multigraph()
    delta = _compute_delta(G, A, B, weight)
    swapped = set()
    gains = []
    while len(swapped) < len(G):
        gain = []
        for u in A - swapped:
            for v in B - swapped:
                try:
                    if multigraph:
                        w = sum(d.get(weight, 1) for d in G[u][v].values())
                    else:
                        w = G[u][v].get(weight, 1)
                except KeyError:
                    w = 0
                gain.append((delta[u] + delta[v] - 2 * w, u, v))
        if len(gain) == 0:
            break
        maxg, u, v = max(gain, key=itemgetter(0))
        swapped |= {u, v}
        gains.append((maxg, u, v))
        delta = _update_delta(delta, G, A - swapped, B - swapped, u, v, weight)
    return gains


# [docs]@not_implemented_for('directed')
def kernighan_lin_bisection(G, partition=None, max_iter=100, weight='weight'):
    """Partition a graph into two blocks using the Kernighan–Lin
    algorithm.

    This algorithm paritions a network into two sets by iteratively
    swapping pairs of nodes to reduce the edge cut between the two sets.

    Parameters
    ----------
    G : graph

    partition : tuple
        Pair of iterables containing an initial partition. If not
        specified, a random balanced partition is used.

    max_iter : int
        Maximum number of times to attempt swaps to find an
        improvemement before giving up.

    weight : key
        Edge data key to use as weight. If None, the weights are all
        set to one.

    Returns
    -------
    partition : tuple
        A pair of sets of nodes representing the bipartition.

    Raises
    -------
    NetworkXError
        If partition is not a valid partition of the nodes of the graph.

    References
    ----------
    .. [1] Kernighan, B. W.; Lin, Shen (1970).
       "An efficient heuristic procedure for partitioning graphs."
       *Bell Systems Technical Journal* 49: 291--307.
       Oxford University Press 2011.

    """
    # If no partition is provided, split the nodes randomly into a
    # balanced partition.
    for div in range(2,8):
        if partition is None:
            nodes = list(G)
            random.shuffle(nodes)
            h = len(nodes) // div
            partition = (nodes[:h], nodes[h:])
        # Make a copy of the partition as a pair of sets.
        try:
            A, B = set(partition[0]), set(partition[1])
        except:
            raise ValueError('partition must be two sets')
        if not is_partition(G, (A, B)):
            raise nx.NetworkXError('partition invalid')
        for i in range(max_iter):
            # `gains` is a list of triples of the form (g, u, v) for each
            # node pair (u, v), where `g` is the gain of that node pair.
            gains = _kernighan_lin_pass(G, A, B, weight)
            csum = list(nx.utils.accumulate(g for g, u, v in gains))
            max_cgain = max(csum)
            if max_cgain <= 0:
                break
            # Get the node pairs up to the index of the maximum cumulative
            # gain, and collect each `u` into `anodes` and each `v` into
            # `bnodes`, for each pair `(u, v)`.
            index = csum.index(max_cgain)
            nodesets = islice(zip(*gains[:index + 1]), 1, 3)
            anodes, bnodes = (set(s) for s in nodesets)
            A |= bnodes
            A -= anodes
            B |= anodes
            B -= bnodes
            print(str(i)+'/'+str((max_iter)))
            color = np.zeros(len(input_nodes))
            for q in range(len(np.array(list(A)))):
                color[np.where(input_nodes == np.array(list(A))[q])] = aa
            nx.draw_networkx(g,with_labels = True,node_color = color,
                                        pos =p,node_size = 100,font_size = 3,font_color = 'w')
            tEnd = time.time()
            plt.title('ratio:'+str(div)+'    epoch:'+str(i)+'     time:'+str(int(tEnd -tStart)))
            plt.savefig('H:/master/code/python/networkScience/week10/pic_kl/{:03d}{}.png'.format(i,aa),format = 'png')
            plt.clf()
            # plt.show()
        return A, B


a = np.genfromtxt(r'H:\master\code\python\networkScience\week10\facebook\0.edges')
g = nx.Graph()
input_nodes = []
for i in a[:,0]:
    if(not (i in input_nodes)):
        input_nodes = np.append(input_nodes,i)
for i in a[:,1]:
    if(not (i in input_nodes)):
        input_nodes = np.append(input_nodes,i)
g.add_nodes_from(input_nodes)
g.add_edges_from(a)
p = nx.kamada_kawai_layout(g)

tStart = time.time()
aa = 1
A,B = kernighan_lin_bisection(g)


gg = nx.Graph()
gg.add_nodes_from(np.array(list(A)))
gg.add_edges_from(g.edges(np.array(list(A))))
aa = 2
A,B = kernighan_lin_bisection(gg)

