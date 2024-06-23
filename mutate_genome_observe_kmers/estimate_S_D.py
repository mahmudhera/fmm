import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import random
import os
import subprocess
import time
import itertools

from read_genome import read_genome
from read_genome import clean_genome_string
from read_genome import get_kmers


def estimate_S_D(genome_string, mutated_genome_string, k):
    # get k-mers in the original string
    kmers_orig = get_kmers(genome_string, k)
    long_kmers_orig = get_kmers(genome_string, k+1)

    kmers_orig_set = set(kmers_orig)
    long_kmers_orig_set = set(long_kmers_orig)

    # get k-mers in the mutated string
    kmers_mutated = get_kmers(mutated_genome_string, k)

    kmers_mutated_set = set(kmers_mutated)

    num_kmers_single_subst = 0
    for kmer_orig in kmers_orig:
        for kmer_mutated in kmers_mutated:
            # generate new kmer by substituting one character in kmer_orig
            for i in range(k):
                for base in ['A', 'C', 'G', 'T']:
                    if base == kmer_orig[i]:
                        continue
                    new_kmer = kmer_orig[:i] + base + kmer_orig[i+1:]
                    if new_kmer in kmers_orig_set:
                        continue
                    if new_kmer == kmer_mutated:
                        num_kmers_single_subst += 1
                        break

    print(num_kmers_single_subst)
            

    

if __name__ == '__main__':
    genome_name = 'ndl.fasta'
    mutated_genome_name = 
    k = 21

    # read the genome and mutated genome
    genome_string = read_genome(genome_name)
    genome_string = clean_genome_string(genome_string)
    mutated_genome_string = read_genome(mutated_genome_name)
    mutated_genome_string = clean_genome_string(mutated_genome_string)

    estimate_S_D(genome_string, mutated_genome_string, k)