import time
from Bio import SeqIO

"""
This script takes the following arguments:
    - genome_file: the path to the genome file
"""

def read_genome(genome_file):
    """
    Reads a genome file and returns the genome as a string
    """
    genome = ""
    for record in SeqIO.parse(genome_file, "fasta"):
        genome += str(record.seq)
    return genome

def clean_genome_string(genome_string):
    """
    Removes all non-alphabet characters from a genome string
    """
    alphabet = set('ACGT')
    return ''.join(filter(alphabet.__contains__, genome_string))
                   
def get_kmers(genome_string, k):
    """
    Returns a list of all k-mers in a genome string
    """
    kmers = []
    for i in range(len(genome_string)-k+1):
        kmers.append(genome_string[i:i+k])
    return kmers

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
    orig_sequence_filename = 'ndl_orig.fasta'
    mutated_sequence_filename = 'ndl_mutated.fasta'
    k = 21

    tick = time.time()

    orig_genome = clean_genome_string(read_genome(orig_sequence_filename))
    mutated_genome = clean_genome_string(read_genome(mutated_sequence_filename))

    kmers_orig = get_kmers(orig_genome, k)
    kmers_mutated = get_kmers(mutated_genome, k)

    kmers_shared_set = set(kmers_orig).intersection(set(kmers_mutated))
    kmers_only_in_orig = set(kmers_orig).difference(set(kmers_mutated))
    kmers_only_in_mutated = set(kmers_mutated).difference(set(kmers_orig))

    tock = time.time()
    print('Time to process genome files:', tock - tick)
    
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