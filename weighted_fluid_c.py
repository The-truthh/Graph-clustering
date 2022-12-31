from collections import Counter
import random

import networkx as nx

__all__ = ['asyn_fluidc']


# Optional to fix the random seed
# random.seed(123)

def asyn_fluidc_new(G, alpha, num_neigbor_depth=1, weight_update=1, max_iter=100):
    if nx.is_connected(G) == True:
        result = asyn_fluidc(G, alpha, num_neigbor_depth, weight_update, max_iter)
        return result
    else:
        result = []
        for comp in nx.connected_components(G):
            sub_G = G.subgraph(comp)
            if len(sub_G) <= 10:
                result.append(list(sub_G))
            else:
                result.extend(asyn_fluidc(sub_G, alpha, num_neigbor_depth, weight_update, max_iter))
        return result


def asyn_fluidc(G, alpha, num_neigbor_depth=1, weight_update=1, max_iter=100):
    """
    Fluid Communities: A Competitive and Highly Scalable Community Detection Algorithm.
    Args:
        - G: Graph to run the algorithm into.
            + type: networkx.Graph
        - k: Number of communities to search.
            + type: int
        - max_iter: Number of maximum iterations allowed.
            + type: int
    Return:
        - List of communities, where each community is a list of vertex ID.
          Each vertex ID can be either an int or str.
            + type: list(list(int or str))
    """
    # Initialization
    max_density = 1.0
    vertices = list(G)
    random.shuffle(vertices)

    degrees = [val for (node, val) in G.degree()]
    degree_thres = max(degrees)  * alpha
    c_neighbors = set()
    core = []
    for vertice in vertices:
        if vertice in c_neighbors:
            continue

        if G.degree[vertice] > degree_thres:
            core.append(vertice)
            for i in range(1, num_neigbor_depth + 1):
                c_neighbors.update(get_neigbors(G, vertice, num_neigbor_depth)[i])
    communities = {n: i for i, n in enumerate(core)}


    density = {}
    com_to_numvertices = {}
    for vertex in communities.keys():
        com_to_numvertices[communities[vertex]] = 1
        density[communities[vertex]] = max_density
    # Set up control variables and start iterating
    iter_count = 0
    cont = True
    while cont:
        cont = False
        iter_count += 1
        # Loop over all vertices in graph in a random order
        vertices = list(G)
        random.shuffle(vertices)
        for vertex in vertices:
            # Updating rule
            com_counter = Counter()
            weight_counter = Counter()
            # Take into account self vertex community
            try:
                com_counter.update({communities[vertex]: density[communities[vertex]]})
            except KeyError:
                pass

            # Gather neighbour vertex communities
            for v in G[vertex]:
                try:
                    weight_temp = G[vertex][v]['weight']
                    p = round(weight_temp, 2) * 100
                    temp = random.randint(1, 100)
                    if temp < p:
                        indicator = 1
                    else:
                        indicator = 0
                    if indicator == 1:
                        com_counter.update({communities[v]: density[communities[v]]})
                        weight_counter.update({communities[v]: weight_temp * weight_update})
                except KeyError:
                    continue

            # Check which is the community with highest density
            new_com = -1
            if len(com_counter.keys()) > 0:
                com_counter.update(weight_counter)
                max_freq = max(com_counter.values())
                best_communities = [com for com, freq in com_counter.items()
                                    if (max_freq - freq) < 0.0001]
                # If actual vertex com in best communities, it is preserved
                try:
                    if communities[vertex] in best_communities:
                        new_com = communities[vertex]
                except KeyError:
                    pass
                # If vertex community changes...
                if new_com == -1:
                    # Set flag of non-convergence
                    cont = True
                    # Randomly chose a new community from candidates
                    new_com = random.choice(best_communities)
                    # Update previous community status
                    try:
                        density[communities[vertex]] = max_density / \
                                                       com_to_numvertices[communities[vertex]]
                    except KeyError:
                        pass
                    # Update new community status
                    communities[vertex] = new_com
                    com_to_numvertices[communities[vertex]] += 1
                    density[communities[vertex]] = max_density / \
                                                       com_to_numvertices[communities[vertex]]
        # If maximum iterations reached --> output actual results
        if iter_count > max_iter:
            print ("Exiting by max iterations!")
            break
    # Return results by grouping communities as list of vertices
    return list(_invert_dict(communities).values())


def _invert_dict(orig_dict):
    """
    Inverting Python dictionary keys and values: Many to one --> One to many
    Args:
        - orig_dict: Dictionary desired to invert.
            + type: dict
    Return:
        - Inverted dictionary
            + type: dict
    """
    return_dict = {}
    for v, k in orig_dict.items():
        try:
            return_dict[k].append(v)
        except KeyError:
            return_dict[k] = [v]
    return return_dict


def get_neigbors(g, node, depth=1):
    output = {}
    layers = dict(nx.bfs_successors(g, source=node, depth_limit=depth))
    nodes = [node]
    for i in range(1,depth+1):
        output[i] = []
        for x in nodes:
            output[i].extend(layers.get(x,[]))
        nodes = output[i]
    return output
