from collections import Counter

def get_num_kmers_single_subst_delt_insert_shared(kmers_orig, kmers_mutated):
    """
    Returns the number of k-mers that are single substitutions, single insertion, single deletion
    """

    # print the number of kmers in the original string and the mutated string
    print(len(kmers_orig), len(kmers_mutated))

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
    for kmer in kmers_orig:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_mutated_1_deletion and kmer not in shared_kmers_set:
                num_kmers_single_substitution += 1
                break

    print('Here..')
    num_kmers_single_substitution_dict = Counter()
    for kmer in kmers_orig:
        # generate all kmers that are 1 substitution away from kmer
        for i in range(len(kmer)):
            for base in ['A', 'C', 'G', 'T']:
                if base == kmer[i]:
                    continue
                if kmer[:i] + base + kmer[i+1:] in mutated_kmers_set:
                    if kmer in num_kmers_single_substitution_dict:
                        num_kmers_single_substitution_dict[kmer] += 1
                    else:
                        num_kmers_single_substitution_dict[kmer] = 1

    # print the top 10 kmers with the most number of single substitutions
    for kmer, num in num_kmers_single_substitution_dict.most_common(10):
        print(kmer, num)

    for kmer in kmers_orig:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_mutated and kmer not in shared_kmers_set:
                num_kmers_single_deletion += 1
                break
    
    # for each k-mer in the mutated string, check if it has a 1 delt k-1 mer in all_k_minus_1_mers_orig
    # if it does, then it is single insertion
    # sets needed: all_k_minus_1_mers_orig
    for kmer in kmers_mutated:
        for i in range(len(kmer)):
            if kmer[:i] + kmer[i+1:] in all_k_minus_1_mers_orig and kmer not in shared_kmers_set:
                num_kmers_single_insertion += 1
                break

    return num_kmers_single_substitution, num_kmers_single_deletion, num_kmers_single_insertion, num_kmers_shared
