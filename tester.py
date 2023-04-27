import networkx as nx
import matplotlib.pyplot as plt
from random import randint, uniform, sample
from struct import pack, unpack
import subprocess
from time import time

G = nx.Graph()
graph_dict = dict()
# community_sizes = [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]
community_sizes = [30, 30, 30, 30, 30, 20, 30, 30, 50, 50, 30, 30, 30, 30, 50, 50, 30, 30, 30, 30, 50, 50, 30, 30, 30]
# community_sizes = [10,10,10,10,10,10]
# community_sizes = [20,20,20,20,20,20,20,20,20]


# set edges between nodes in same community
# G - nx graph
# community_nodes - list of nodes in community
# probability - float from (0,1) interval
def commune(G, community_nodes, probability):
    for u in range(len(community_nodes)):
        for v in range(u+1, len(community_nodes)):
            if uniform(0, 1) < probability:
                try:
                    graph_dict[community_nodes[u]].add(community_nodes[v])
                except:
                    graph_dict[community_nodes[u]] = {community_nodes[v]}

                try:
                    graph_dict[community_nodes[v]].add(community_nodes[u])
                except:
                    graph_dict[community_nodes[v]] = {community_nodes[u]}

                G.add_edge(community_nodes[u],community_nodes[v])


# set edges between nodes in different communities
# G - nx graph
# comm - list of lists of community nodes
# probability - number between 0 and 1
def connect_communities(G, comm, probability):
    for c1 in range(len(comm)):
        for c2 in range(c1, len(comm)):
            s = sample(comm[c1], int(len(comm[c1]) * probability))
            d = sample(comm[c2], int(len(comm[c2]) * probability))
            for w in s:
                for q in d:
                    if uniform(0,1) < probability:
                        # print("(" + str(w) + "," + str(q) + ")")
                        try:
                            graph_dict[w].add(q)
                        except:
                            graph_dict[w] = {q}

                        try:
                            graph_dict[q].add(w)
                        except:
                            graph_dict[q] = {w}

                        G.add_edge(w, q)


# set the participating nodes in each community
# community_sizes - list of integers of number of nodes in each community
# Returns: list of lists of community nodes
def set_communities(community_sizes):
    num_of_nodes = sum(community_sizes)
    avail_nodes = {x for x in range(num_of_nodes)}
    communities = []

    for size in community_sizes:
        current = set(sample(avail_nodes, size))
        avail_nodes = avail_nodes.difference(current)
        current = list(current)

        with open('test.in.txt', 'a') as f:
            i = 0
            for item in current:
                f.write("%s" % item)
                if i<len(current)-1:
                    f.write(", ")
                i += 1
            f.write("\n")

        communities.append(current)

    return communities


def read_text(file_name):
    with open(file_name, 'r') as f:
        lines = [line.rstrip().split(", ") for line in f]
        print(sorted(lines))
    return lines


def run_c_program():
    subprocess.call("final.exe test.in out")


def graph():
    G = nx.Graph()

    number_of_nodes = 0

    is_first = True
    count = 0
    curr_node = -1

    with open("test.in", "rb") as f:
        for byte in iter(lambda: f.read(4), ''):
            if byte == b'':
                break

            val = unpack('i', byte)[0]

            if is_first:  # If first byte then store total number of nodes
                number_of_nodes = val
                is_first = False
                continue

            if count == 0:  # If moved to next node
                count = val
                curr_node += 1
            else:  # If in current node's edges
                G.add_edge(curr_node, val)
                count -= 1

    is_first = True
    count = 0
    curr_group = 0

    color_map = [None for x in range(number_of_nodes)]
    color_dict = {1: 'lightcoral',
                  2: 'lightgreen',
                  3: 'skyblue',
                  4: 'wheat',
                  5: 'orchid',
                  6: 'plum',
                  7: 'red',
                  8: 'green',
                  9: 'blue',
                  10: 'lime',
                  11: 'aqua',
                  12: 'gold',
                  13: 'darkred',
                  14: 'darkseagreen',  #
                  15: 'olive',
                  16: 'lightsalmon',
                  17: 'deepskyblue',
                  18: 'teal',
                  19: 'khaki',
                  20: 'magenta',
                  21: 'springgreen',
                  22: 'yellow',
                  23: 'dodgerblue',
                  24: 'silver',
                  25: 'crimson',
                  26: 'powderblue',
                  27: 'navy',
                  28: 'dimgrey',
                  29: 'black',
                  30: 'tan',
                  31: 'black'}

    with open("out", "rb") as f:
        for byte in iter(lambda: f.read(4), ''):
            if byte == b'':
                break

            val = unpack('i', byte)[0]

            if is_first:  # If first byte then store total number of groups
                number_of_groups = val
                is_first = False
                continue

            if count == 0:  # If moved to next group
                count = val
                curr_group += 1
            else:  # If in current group
                color_map[val] = curr_group
                count -= 1

    color_map = [color_dict[color_map[x]] for x in range(len(color_map))]
    nx.draw_networkx(G, nodelist=[x for x in range(number_of_nodes)], node_color=color_map, with_labels=True)
    plt.show()


def check_partitioning():
    output_groups = []
    is_first = True
    count = 0
    curr_group = []

    # get partitions from algorithms output file
    with open("out", 'rb') as f:
        for byte in iter(lambda: f.read(4), ''):
            if byte == b'':
                if len(curr_group) > 0:
                    output_groups.append(curr_group)

                break

            val = unpack('i', byte)[0]

            if is_first:
                is_first = False
                continue

            if count == 0:
                count = val

                if len(curr_group) > 0:
                    output_groups.append(curr_group)
                    curr_group = []

            else:
                curr_group.append(val)
                count -= 1

    diff_nodes = set()

    # compare with original generator
    for group in output_groups:
        min_diff = None

        for community in communities:
            curr_diff = set(group).symmetric_difference(set(community))

            if min_diff is None:
                min_diff = curr_group
            if len(min_diff) > len(curr_diff):
                min_diff = curr_diff

        diff_nodes.update(min_diff)

    print(str(len(diff_nodes)) + " diffs: ", end='')
    print(sorted(diff_nodes))


start = time()

with open('test.in.txt', 'w') as f:
    pass

communities = set_communities(community_sizes)

for community in communities:
    commune(G, community, 0.8)

connect_communities(G, communities, 0.2)
lines = read_text("test.in.txt")

with open("test.in", "wb") as f:
    s = pack('i', len(graph_dict))
    f.write(s)
    for key in sorted(graph_dict.keys()):
        s = pack('i', len(graph_dict[key]))
        f.write(s)
        for neighbor in sorted(graph_dict[key], key=int):
            s = pack('i', int(neighbor))
            f.write(s)

run_c_program()
end = time()

check_partitioning()

print("elapsed time: " + str(end-start))

# graph()

# nx.draw(G, with_labels=True)
# plt.show()
