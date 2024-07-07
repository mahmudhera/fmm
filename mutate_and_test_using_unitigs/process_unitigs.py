"""
This script takes two fasta files as inputs
Builds the bdg graph for these two inputs
Goes over all unitigs
Compares unitig vs unitig
And writes num of kmers with single subst and delt
in an output file
"""


import time
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import multiprocessing
import random
import argparse
import os

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




def process_unitigs(unitigs_orig_subset, unitigs_mutated, list_of_unitig_lengths_mutated, k, multiplier):
    num_kmers_single_subst = 0
    num_kmers_single_delt = 0
    results = []
    all_alignments = []
    kmers_marked_for_subst = []

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
        #for unitig2 in unitigs_mutated:
            alignment = pairwise2.align.localms(unitig1, unitig2, 3, -1, -1, -1)[0]
            if alignment.score > best_match_score:
                best_match_score = alignment.score
                best_match_alignment = alignment

            unitig2 = reverse_complement(unitig2)
            alignment = pairwise2.align.localms(unitig1, unitig2, 3, -1, -1, -1)[0]
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
                kmers_marked_for_subst.append(seqA[i:i+k])

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] == '-' and seqA[i] in alphabet:
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_delt += 1

        all_alignments.append((unitig1, format_alignment(*alignment)))

    results.append((num_kmers_single_subst, num_kmers_single_delt, all_alignments, kmers_marked_for_subst))

    return results



def main3():

    # parse command line arguments for filenames and k
    parser = argparse.ArgumentParser()
    parser.add_argument('orig_filename', type=str, help='Filename of original string')
    parser.add_argument('mut_filename', type=str, help='Filename of mutated string')
    parser.add_argument('k', type=int, help='Length of k-mers')
    parser.add_argument('o', type=str, help='Output filename')

    args = parser.parse_args()

    orig_filename = args.orig_filename
    mut_filename = args.mut_filename
    k = args.k
    output_filename = args.o
    multiplier = 3.0

    # invoke cuttlefish2 with the files and generate the unitigs
    # command: cuttlefish build -s <filename> -k <k> -t <thread_count> -o <out_filename> -w . --ref
    unitigs_orig_filename = orig_filename + '_unitigs'
    unitigs_mutated_filename = mut_filename + '_unitigs'

    # remove files if they exist
    if os.path.exists(unitigs_orig_filename+'.fa'):
        os.system('rm ' + unitigs_orig_filename + '*')
    if os.path.exists(unitigs_mutated_filename+'.fa'):
        os.system('rm ' + unitigs_mutated_filename + '*')

    cmd1 = 'cuttlefish build -s {orig_filename} -k {k} -t 128 -o {unitigs_orig_filename} -w . --ref'.format(orig_filename=orig_filename, k=k, unitigs_orig_filename=unitigs_orig_filename)
    cmd2 = 'cuttlefish build -s {mut_filename} -k {k} -t 128 -o {unitigs_mutated_filename} -w . --ref'.format(mut_filename=mut_filename, k=k, unitigs_mutated_filename=unitigs_mutated_filename)
    os.system(cmd1)
    os.system(cmd2)

    unitigs_orig = read_unitigs(unitigs_orig_filename+'.fa')
    unitigs_mutated = read_unitigs(unitigs_mutated_filename+'.fa')

    # Sort the unitigs by length, small to large
    unitigs_mutated = sorted(list(unitigs_mutated), key=lambda x: len(x))
    list_of_unitig_lengths_mutated = [len(unitig) for unitig in unitigs_mutated]

    # Divide the work among multiple processes
    num_cores = multiprocessing.cpu_count()
    num_cores = min(num_cores, len(unitigs_orig))

    chunk_size = len(unitigs_orig) // num_cores
    unitigs_orig_list = list(unitigs_orig)
    random.shuffle(unitigs_orig_list)
    chunks = [unitigs_orig_list[i:i + chunk_size] for i in range(0, len(unitigs_orig), chunk_size)]

    pool = multiprocessing.Pool(num_cores)
    args = []
    for i, chunk in enumerate(chunks):
        args.append((chunk, unitigs_mutated, list_of_unitig_lengths_mutated, k, multiplier))
    results = pool.starmap(process_unitigs, args)

    print('Done, now aggregating results'
          ' and writing to output file')

    # Aggregate results
    f = open('alignments', 'w')
    f2 = open('substitution_kmers', 'w')
    num_kmers_single_subst = 0
    num_kmers_single_delt = 0
    for result in results:
        for subst, delt, alignment, kmers_subst in result:
            num_kmers_single_subst += subst
            num_kmers_single_delt += delt
            for u, a in alignment:
                f.write(u)
                f.write('\n')
                f.write(a)
                f.write('\n\n')
            for kmer in kmers_subst:
                f2.write(kmer)
                f2.write('\n')
    f.close()
    f2.close()

    # Write results to output file
    with open(output_filename, 'w') as f:
        f.write(str(num_kmers_single_subst) + '\n')
        f.write(str(num_kmers_single_delt) + '\n')



if __name__ == '__main__':
    main3()