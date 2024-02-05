from collections import Counter
from scipy.sparse import csr_matrix 
from scipy.sparse.csgraph import maximum_bipartite_matching


def get_num_kmers_single_subst_delt_insert_shared(kmers_orig, kmers_mutated):
    """
    Returns the number of k-mers that are single substitutions, single insertion, single deletion
    """

    # print the number of kmers in the original string and the mutated string
    # print(len(kmers_orig), len(kmers_mutated))

    num_kmers_single_substitution = 0
    num_kmers_single_deletion = 0
    num_kmers_single_insertion = 0
    num_kmers_shared = 0

    # get num of shared kmers
    mutated_kmers_set = set(kmers_mutated)
    orig_kmers_set = set(kmers_orig)
    shared_kmers_set = orig_kmers_set.intersection(mutated_kmers_set)
    num_kmers_shared = len(shared_kmers_set)
    
    # create a hash set of all k-1 -mers in the original string
    all_k_minus_1_mers_orig = set()
    for kmer in kmers_orig:
        all_k_minus_1_mers_orig.add( kmer[1:] )
        all_k_minus_1_mers_orig.add( kmer[:-1] )

    # create a hash set of all k-1 -mers in the mutated string
    all_k_minus_1_mers_mutated = set()
    for kmer in kmers_mutated:
        all_k_minus_1_mers_mutated.add( kmer[1:] )
        all_k_minus_1_mers_mutated.add( kmer[:-1] )

    # create a hash set of all k-1 -mers in the mutated string that is 1 deletion away from a k-mer in the mutated string
    all_k_minus_1_mers_mutated_1_deletion = set()
    for kmer in kmers_mutated:
        for i in range(len(kmer)):
            all_k_minus_1_mers_mutated_1_deletion.add( kmer[:i] + kmer[i+1:] )

    # iterate through all k-mers in the original string, check if it has a 1 delt k-1 mer in all_k_minus_1_mers_mutated_1_deletion
    # if it does, then it is single substitution
    # sets needed: all_k_minus_1_mers_mutated_1_deletion
    # for each k-mer in the original string, check if it has a 1 delt k-1 mer in all_k_minus_1_mers_mutated
    # if it does, then it is single deletion
    # sets needed: all_k_minus_1_mers_mutated
    for kmer in orig_kmers_set:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_mutated_1_deletion and kmer not in shared_kmers_set:
                num_kmers_single_substitution += 1
                break

    print('Here..')
    num_kmers_single_substitution_dict = Counter()
    new_kmers_single_substitution_dict = Counter()
    for kmer in orig_kmers_set:
        #if kmer in shared_kmers_set:
        #    continue
        # generate all kmers that are 1 substitution away from kmer
        for i in range(len(kmer)):
            for base in ['A', 'C', 'G', 'T']:
                if base == kmer[i]:
                    continue
                new_kmer = kmer[:i] + base + kmer[i+1:]
                if new_kmer in mutated_kmers_set:
                    #if new_kmer in shared_kmers_set:
                    #    continue
                    if kmer in num_kmers_single_substitution_dict:
                        num_kmers_single_substitution_dict[kmer] += 1
                        new_kmers_single_substitution_dict[new_kmer] += 1
                    else:
                        num_kmers_single_substitution_dict[kmer] = 1
                        new_kmers_single_substitution_dict[new_kmer] = 1

    # print the top 10 kmers with the most number of single substitutions
    print ('Mostly occuring 10 kmers marked for subst. (from the set of orig kmers:)')
    for kmer, num in num_kmers_single_substitution_dict.most_common(10):
        print(kmer, num)

    print ('Mostly occuring 10 kmers marked for subst. (from the set of new kmers:)')
    for kmer, num in new_kmers_single_substitution_dict.most_common(10):
        print(kmer, num)

    # print the number of kmers in num_kmers_single_substitution_dict
    print( 'number of kmers marked for single char subst. in the orig and the new kmers, respectively: ' )
    print(len(num_kmers_single_substitution_dict), len(new_kmers_single_substitution_dict))


    ####################
    # maximal bipartite matching start
    ####################

    # create a graph with the orig kmers as the left nodes and the new kmers as the right nodes
    # the edges are between the orig kmer and the new kmer if the new kmer is 1 substitution away from the orig kmer
    # the weight of the edge is the number of times the new kmer is 1 substitution away from the orig kmer
    # the goal is to find the maximum matching of the graph
    # the maximum matching will give us the number of single substitutions
    # the maximum matching is the number of edges in the graph that are not overlapping

    print('Constructing the graph..')

    left_kmers_to_indices = {}
    right_kmers_to_indices = {}
    double_indices_to_edges = {}

    left_kmers = list(num_kmers_single_substitution_dict.keys())
    for i, kmer in enumerate(left_kmers):
        left_kmers_to_indices[kmer] = i

    right_kmers = list(new_kmers_single_substitution_dict.keys())
    for i, kmer in enumerate(right_kmers):
        right_kmers_to_indices[kmer] = i

    left_kmers_set = set(left_kmers)
    right_kmers_set = set(right_kmers)

    for kmer in left_kmers_set:
        # generate all kmers that are 1 substitution away from kmer
        for i in range(len(kmer)):
            for base in ['A', 'C', 'G', 'T']:
                if base == kmer[i]:
                    continue
                new_kmer = kmer[:i] + base + kmer[i+1:]
                if new_kmer in right_kmers_set:
                    left_kmer_index = left_kmers_to_indices[kmer]
                    right_kmer_index = right_kmers_to_indices[new_kmer]
                    if (left_kmer_index, right_kmer_index) in double_indices_to_edges:
                        double_indices_to_edges[(left_kmer_index, right_kmer_index)] += 1
                    else:
                        double_indices_to_edges[(left_kmer_index, right_kmer_index)] = 1

    # get the left indices, right indices, and the weights of the edges
    left_indices = []
    right_indices = []
    weights = []

    for (left_index, right_index), weight in double_indices_to_edges.items():
        left_indices.append(left_index)
        right_indices.append(right_index)
        weights.append(weight)

    # create a csr matrix
    num_left_nodes = len(left_kmers)
    num_right_nodes = len(right_kmers)
    graph_matrix = csr_matrix((weights, (left_indices, right_indices)), shape=(num_left_nodes, num_right_nodes))

    print ('Graph created., now finding the maximum matching..')

    # get the maximum matching
    max_matching = maximum_bipartite_matching(graph_matrix, perm_type='column')

    # get the number of values that are not -1
    num_edges_max_matching = 0
    for value in max_matching:
        if value != -1:
            num_edges_max_matching += 1

    print('Number of edges in the maximum matching: ', num_edges_max_matching)

    
    ####################
    # maximal bipartite matching end
    ####################

    for kmer in orig_kmers_set:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_mutated and kmer not in shared_kmers_set:
                num_kmers_single_deletion += 1
                break
    
    # for each k-mer in the mutated string, check if it has a 1 delt k-1 mer in all_k_minus_1_mers_orig
    # if it does, then it is single insertion
    # sets needed: all_k_minus_1_mers_orig
    for kmer in orig_kmers_set:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_orig and kmer not in shared_kmers_set:
                num_kmers_single_insertion += 1
                break

    return num_kmers_single_substitution, num_kmers_single_deletion, num_kmers_single_insertion, num_kmers_shared
