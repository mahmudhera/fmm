import sys

from read_genome import read_genome
from read_genome import clean_genome_string
from read_genome import get_kmers


def estimate_S_D(genome_string, mutated_genome_string, k):
    # get k-mers in the original string
    long_kmers_orig = get_kmers(genome_string, k+2)
    kmers_orig_set = set(get_kmers(genome_string, k))

    orig_kmer_to_long_kmer_where_kmer_prefix = {}
    orig_kmer_to_long_kmer_where_kmer_suffix = {}

    for long_kmer in long_kmers_orig:
        orig_kmer_to_long_kmer_where_kmer_prefix[long_kmer[:k]] = long_kmer
        orig_kmer_to_long_kmer_where_kmer_suffix[long_kmer[1:]] = long_kmer

    # get k-mers in the mutated string
    kmers_mutated = get_kmers(mutated_genome_string, k)
    kmers_mutated_set = set(kmers_mutated)

    mut_kmer_to_long_kmer_where_kmer_prefix = {}
    mut_kmer_to_long_kmer_where_kmer_suffix = {}

    for long_kmer in get_kmers(mutated_genome_string, k+1):
        mut_kmer_to_long_kmer_where_kmer_prefix[long_kmer[:k]] = long_kmer
        mut_kmer_to_long_kmer_where_kmer_suffix[long_kmer[2:]] = long_kmer

    num_kmers_single_subst = 0
    for kmer_orig in kmers_orig_set:
        # generate new kmer by substituting one character in kmer_orig
        for i in range(k):
            for base in ['A', 'C', 'G', 'T']:
                if base == kmer_orig[i]:
                    continue
                new_kmer = kmer_orig[:i] + base + kmer_orig[i+1:]
                if new_kmer in kmers_orig_set:
                    continue
                if new_kmer in kmers_mutated_set:
                    if i != 0 and i != k-1:
                        num_kmers_single_subst += 1
                    elif i == 0:
                        if new_kmer not in mut_kmer_to_long_kmer_where_kmer_suffix or new_kmer not in orig_kmer_to_long_kmer_where_kmer_suffix:
                            continue
                        long_kmer_orig = orig_kmer_to_long_kmer_where_kmer_suffix[kmer_orig]
                        long_kmer_mut = mut_kmer_to_long_kmer_where_kmer_suffix[new_kmer]
                        if long_kmer_orig[0:1] == long_kmer_mut[0:1]:
                            num_kmers_single_subst += 1
                    else:
                        if new_kmer not in mut_kmer_to_long_kmer_where_kmer_prefix or new_kmer not in orig_kmer_to_long_kmer_where_kmer_prefix:
                            continue
                        long_kmer_orig = orig_kmer_to_long_kmer_where_kmer_prefix[kmer_orig]
                        long_kmer_mut = mut_kmer_to_long_kmer_where_kmer_prefix[new_kmer]
                        if long_kmer_orig[-2:-1] == long_kmer_mut[-2:-1]:
                            num_kmers_single_subst += 1

    num_kmers_single_delt = 0
    k_plus_one_mers_orig = get_kmers(genome_string, k+1)
    for k_plus_one_mer in k_plus_one_mers_orig:
        for i in range(k+1):
            new_kmer = k_plus_one_mer[:i] + k_plus_one_mer[i+1:]
            if new_kmer in kmers_orig_set:
                continue
            if new_kmer in kmers_mutated_set:
                num_kmers_single_delt += 1

    print(num_kmers_single_subst)
    print(num_kmers_single_delt)
            

    

if __name__ == '__main__':
    genome_name = 'ndl.fasta'
    mutated_genome_name = sys.argv[1]
    k = 21

    # read the genome and mutated genome
    genome_string = read_genome(genome_name)
    genome_string = clean_genome_string(genome_string)
    mutated_genome_string = read_genome(mutated_genome_name)
    mutated_genome_string = clean_genome_string(mutated_genome_string)

    estimate_S_D(genome_string, mutated_genome_string, k)