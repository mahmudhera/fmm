import time

def read_unitigs(unitigs_file):
    unitigs = set()
    with open(unitigs_file) as f:
        for line in f:
            if line[0] == '>':
                continue
            else:
                unitigs.add(line.strip())
    return unitigs

def build_de_bruijn_graph(unitigs, k=21):
    graph = {}
    children = {}
    parent = {}

    unitigs_list = list(unitigs)

    for i in range(len(unitigs)):
        unitig1 = unitigs_list[i]
        for j in range(i+1, len(unitigs)):
            unitig2 = unitigs_list[j]
            if unitig1[-k+1:] == unitig2[:k-1]:
                if unitig1 not in graph:
                    graph[unitig1] = []
                graph[unitig1].append(unitig2)
                if unitig1 not in children:
                    children[unitig1] = []
                children[unitig1].append(unitig2)
                if unitig2 not in parent:
                    parent[unitig2] = []
                parent[unitig2].append(unitig1)
            if unitig2[-k+1:] == unitig1[:k-1]:
                if unitig2 not in graph:
                    graph[unitig2] = []
                graph[unitig2].append(unitig1)
                if unitig2 not in children:
                    children[unitig2] = []
                children[unitig2].append(unitig1)
                if unitig1 not in parent:
                    parent[unitig1] = []
                parent[unitig1].append(unitig2)
    return graph, children, parent


def main():
    unitigs_file = 'cdbg.fa'
    
    tick = time.time()
    unitigs = read_unitigs(unitigs_file)
    print('Time to read unitigs:', time.time() - tick)

    tick = time.time()
    graph, children, parent = build_de_bruijn_graph(unitigs)
    print('Time to build graph:', time.time() - tick)

    print('Number of unitigs:', len(unitigs))
    print('Number of edges:', sum([len(graph[u]) for u in graph]))
    print('Number of children:', sum([len(children[u]) for u in children]))
    print('Number of parents:', sum([len(parent[u]) for u in parent]))


if __name__ == '__main__':
    main()