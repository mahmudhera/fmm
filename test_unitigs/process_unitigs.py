import time
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import multiprocessing

alphabet = set('ACGT')

def reverse_complement(kmer):
    """
    Returns the reverse complement of a k-mer
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(kmer))


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


def assign_colors_to_unitigs(unitigs, kmers_shared_set, kmers_only_in_orig, kmers_only_in_mutated, k):
    colors = {}
    unitigs_in_both = []
    unitigs_in_orig = []
    unitigs_in_mutated = []
    for unitig in unitigs:
        unitig_one_kmet = unitig[:k]
        revcomp_kmer = reverse_complement(unitig_one_kmet)
        if unitig_one_kmet in kmers_shared_set or revcomp_kmer in kmers_shared_set:
            colors[unitig] = 'both'
            unitigs_in_both.append(unitig)
        elif unitig_one_kmet in kmers_only_in_orig or revcomp_kmer in kmers_only_in_orig:
            colors[unitig] = 'orig'
            unitigs_in_orig.append(unitig)
        elif unitig_one_kmet in kmers_only_in_mutated or revcomp_kmer in kmers_only_in_mutated:
            colors[unitig] = 'mutated'
            unitigs_in_mutated.append(unitig)
        else:
            colors[unitig] = 'none'
    return colors, unitigs_in_both, unitigs_in_orig, unitigs_in_mutated


def find_bridging_unitigs(unitigs_in_orig, unitigs_in_mutated, graph, children, parent):
    bridging_unitig_pairs = []
    for unitig_orig in unitigs_in_orig:
        for unitig_mutated in unitigs_in_mutated:
            try:
                unitig_orig_parents = parent[unitig_orig]
                unitig_mutated_parents = parent[unitig_mutated]
                unitig_orig_children = children[unitig_orig]
                unitig_mutated_children = children[unitig_mutated]

                if len(set(unitig_orig_parents).intersection(set(unitig_mutated_parents))) > 0 and len(set(unitig_mutated_children).intersection(set(unitig_orig_children))) > 0:
                    bridging_unitig_pairs.append((unitig_orig, unitig_mutated))
            except KeyError:
                continue
    return bridging_unitig_pairs


def count_kmers_single_subst_delt_insert(bridging_unitig_pairs, k, kmers_orig_count):
    num_single_subst = 0
    num_single_delt = 0
    num_single_insert = 0

    alphabel = set('ACGT')

    for unitig_orig, unitig_mutated in bridging_unitig_pairs:
        alignment1 = pairwise2.align.globalms(unitig_orig, unitig_mutated, 3, -1, -1, -1)[0]
        alignment2 = pairwise2.align.globalms(unitig_orig, reverse_complement(unitig_mutated), 3, -1, -1, -1)[0]

        if alignment1.score > alignment2.score:
            alignment = alignment1
        else:
            alignment = alignment2

        seqA = alignment.seqA
        seqB = alignment.seqB
        
        num_chars = len(seqA)
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqA[i] in alphabet and seqB[i] in alphabet:
                    in_numbers[i] = 1

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                kmer_of_interest = seqA[i:].replace('-', '')[:k]
                if kmer_of_interest in kmers_orig_count:
                    num_single_subst += kmers_orig_count[kmer_of_interest]
                else:
                    try:
                        kmer_of_interest = reverse_complement(kmer_of_interest)
                        num_single_subst += kmers_orig_count[kmer_of_interest]
                    except:
                        print('Error in subst.:', kmer_of_interest)
                        num_single_subst += 1

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqA[i] == '-':
                    in_numbers[i] = 1

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                kmer_of_interest = seqA[i:].replace('-', '')[:k]
                if kmer_of_interest in kmers_orig_count:
                    num_single_insert += kmers_orig_count[kmer_of_interest]
                else:
                    try:
                        kmer_of_interest = reverse_complement(kmer_of_interest)
                        num_single_insert += kmers_orig_count[kmer_of_interest]
                    except:
                        print('Error in insert:', kmer_of_interest)
                        num_single_insert += 1

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqB[i] == '-':
                    in_numbers[i] = 1

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                kmer_of_interest = seqA[i:].replace('-', '')[:k]
                if kmer_of_interest in kmers_orig_count:
                    num_single_delt += kmers_orig_count[kmer_of_interest]
                else:
                    try:
                        kmer_of_interest = reverse_complement(kmer_of_interest)
                        num_single_delt += kmers_orig_count[kmer_of_interest]
                    except:
                        print('Error in delt:', kmer_of_interest)
                        num_single_delt += 1

    return num_single_subst, num_single_delt, num_single_insert

        
        
            
            
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
    colors, unitigs_in_both, unitigs_in_orig, unitigs_in_mutated = assign_colors_to_unitigs(unitigs, kmers_shared_set, kmers_only_in_orig, kmers_only_in_mutated, k)
    print('Time to assign colors:', time.time() - tick)
    # show five random unitigs
    print('Some sample colors:')
    for i in range(5):
        print(list(unitigs)[i], colors[list(unitigs)[i]])

    print('Number of unitigs in both:', len(unitigs_in_both))
    print('Number of unitigs in orig only:', len(unitigs_in_orig))
    print('Number of unitigs in mutated only:', len(unitigs_in_mutated))

    tick = time.time()
    graph, children, parent = build_de_bruijn_graph(unitigs)
    print('Time to build graph:', time.time() - tick)

    print('Number of unitigs:', len(unitigs))
    print('Number of edges:', sum([len(graph[u]) for u in graph]))
    print('Number of children:', sum([len(children[u]) for u in children]))
    print('Number of parents:', sum([len(parent[u]) for u in parent]))

    # show five random graph entries
    for i in range(5):
        print(list(graph.keys())[i], graph[list(graph.keys())[i]])

    tick = time.time()
    bridging_unitig_pairs = find_bridging_unitigs(unitigs_in_orig, unitigs_in_mutated, graph, children, parent)
    print('Time to find bridging unitigs:', time.time() - tick)
    print('Number of bridging unitig pairs:', len(bridging_unitig_pairs))
    # show first five bridging unitig pairs
    print('Some sample bridging unitig pairs:')
    for i in range(5):
        print(bridging_unitig_pairs[i])

    # get count of kmers in kmers_orig
    kmers_orig_count = {}
    for kmer in kmers_orig:
        if kmer in kmers_orig_count:
            kmers_orig_count[kmer] += 1
        else:
            kmers_orig_count[kmer] = 1

    tick = time.time()
    num_single_subst, num_single_delt, num_single_insert = count_kmers_single_subst_delt_insert(bridging_unitig_pairs, k, kmers_orig_count)
    print('Time to count single subst/delt/insert:', time.time() - tick)
    print('Number of single substitutions:', num_single_subst)
    print('Number of single deletions:', num_single_delt)
    print('Number of single insertions:', num_single_insert)


# returns index such that sorted_list[index] <= query < sorted_list[index+1]
def find_index(sorted_list, query):
    low, high = 0, len(sorted_list) - 1
    
    # If the query is out of bounds, handle edge cases
    if query < sorted_list[0]:
        return -1
    if query > sorted_list[-1]:
        return len(sorted_list)
    
    while low <= high:
        mid = (low + high) // 2
        
        if sorted_list[mid] == query:
            return mid
        elif sorted_list[mid] < query:
            low = mid + 1
        else:
            high = mid - 1
    
    return high



def main2():
    unitigs_orig_filename = 'cdbg_orig.fa'
    unitigs_mutated_filename = 'cdbg_mutated.fa'
    k = 21
    multiplier = 3.0

    unitigs_orig = read_unitigs(unitigs_orig_filename)
    unitigs_mutated = read_unitigs(unitigs_mutated_filename)

    # sort the unitigs by length, small to large
    unitigs_orig = sorted(list(unitigs_orig), key=lambda x: len(x))
    unitigs_mutated = sorted(list(unitigs_mutated), key=lambda x: len(x))

    list_of_unitig_lengths_orig = [len(unitig) for unitig in unitigs_orig]
    list_of_unitig_lengths_mutated = [len(unitig) for unitig in unitigs_mutated]

    num_kmers_single_subst = 0
    num_kmers_single_delt = 0

    f = open('log_alignments', 'w')
    f2 = open('log_lengths', 'w')

    num_completed = 0
    print('Number of unitigs:', len(unitigs_orig))
    for unitig1 in unitigs_orig:
        best_match_score = -999999999
        best_match_alignment = None
        best_match_length = -1

        unitig_2_length_low = len(unitig1) / multiplier
        unitig_2_length_high = len(unitig1) * multiplier

        low_index = find_index(list_of_unitig_lengths_mutated, unitig_2_length_low)
        high_index = find_index(list_of_unitig_lengths_mutated, unitig_2_length_high)

        low_index = max(0, low_index)
        high_index = min(len(unitigs_mutated)-1, high_index)

        for unitig2 in unitigs_mutated[low_index:high_index+1]:
            alignment = pairwise2.align.globalms(unitig1, unitig2, 3, -1, -1, -1)[0]
            if alignment.score > best_match_score:
                best_match_score = alignment.score
                best_match_alignment = alignment
                best_match_length = len(unitig2)
            unitig2 = reverse_complement(unitig2)
            alignment = pairwise2.align.globalms(unitig1, unitig2, 3, -1, -1, -1)[0]
            if alignment.score > best_match_score:
                best_match_score = alignment.score
                best_match_alignment = alignment
                best_match_length = len(unitig2)
        
        alignment = best_match_alignment
        seqA = alignment.seqA
        seqB = alignment.seqB
        
        num_chars = len(seqA)
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqA[i] in alphabet and seqB[i] in alphabet:
                    in_numbers[i] = 1
                else:
                    in_numbers[i] = 2
        
        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_subst += 1

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] == '-' and seqA[i] in alphabet:
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_delt += 1

        f.write(seqA + '\n')
        f.write(seqB + '\n')
        f.write('--\n')

        f2.write(str(len(unitig1)) + ',' + str(best_match_length) + '\n')


        # show progress
        num_completed += 1
        print('Completed {num_completed} unitigs'.format(num_completed=num_completed), end='\r')

    f.close()
    f2.close()

    print()
    print(num_kmers_single_subst)
    print(num_kmers_single_delt)



def process_unitigs(unitigs_orig_subset, unitigs_mutated, list_of_unitig_lengths_mutated, k, multiplier):
    num_kmers_single_subst = 0
    num_kmers_single_delt = 0
    results = []

    for unitig1 in unitigs_orig_subset:
        best_match_score = -999999999
        best_match_alignment = None

        unitig_2_length_low = len(unitig1) / multiplier
        unitig_2_length_high = len(unitig1) * multiplier

        low_index = find_index(list_of_unitig_lengths_mutated, unitig_2_length_low)
        high_index = find_index(list_of_unitig_lengths_mutated, unitig_2_length_high)

        low_index = max(0, low_index-1)
        high_index = min(len(unitigs_mutated)-1, high_index+1)

        for unitig2 in unitigs_mutated[low_index:high_index+1]:
            alignment = pairwise2.align.globalms(unitig1, unitig2, 3, -1, -1, -1)[0]
            if alignment.score > best_match_score:
                best_match_score = alignment.score
                best_match_alignment = alignment

            unitig2 = reverse_complement(unitig2)
            alignment = pairwise2.align.globalms(unitig1, unitig2, 3, -1, -1, -1)[0]
            if alignment.score > best_match_score:
                best_match_score = alignment.score
                best_match_alignment = alignment
        
        alignment = best_match_alignment
        seqA = alignment.seqA
        seqB = alignment.seqB
        
        num_chars = len(seqA)
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqA[i] in alphabet and seqB[i] in alphabet:
                    in_numbers[i] = 1
                else:
                    in_numbers[i] = 2
        
        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_subst += 1

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] == '-' and seqA[i] in alphabet:
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_delt += 1
        
        results.append((num_kmers_single_subst, num_kmers_single_delt))

    return results



def main3():
    unitigs_orig_filename = 'cdbg_orig.fa'
    unitigs_mutated_filename = 'cdbg_mutated.fa'
    k = 21
    multiplier = 3.0

    unitigs_orig = read_unitigs(unitigs_orig_filename)
    unitigs_mutated = read_unitigs(unitigs_mutated_filename)

    # Sort the unitigs by length, small to large
    unitigs_mutated = sorted(list(unitigs_mutated), key=lambda x: len(x))
    list_of_unitig_lengths_mutated = [len(unitig) for unitig in unitigs_mutated]

    # Divide the work among multiple processes
    num_cores = multiprocessing.cpu_count()
    chunk_size = len(unitigs_orig) // num_cores
    unitigs_orig_list = list(unitigs_orig)
    chunks = [unitigs_orig_list[i:i + chunk_size] for i in range(0, len(unitigs_orig), chunk_size)]

    pool = multiprocessing.Pool(num_cores)
    results = pool.starmap(process_unitigs, [(chunk, unitigs_mutated, list_of_unitig_lengths_mutated, k, multiplier) for chunk in chunks])

    # Aggregate results
    num_kmers_single_subst = 0
    num_kmers_single_delt = 0
    for result in results:
        for subst, delt in result:
            num_kmers_single_subst += subst
            num_kmers_single_delt += delt

    print(num_kmers_single_subst)
    print(num_kmers_single_delt)




if __name__ == '__main__':
    main3()