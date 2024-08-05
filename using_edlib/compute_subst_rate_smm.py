"""
In this script: we will compute the subst rate only, assuming the
simple mutation model
"""

import argparse

from itertools import product

from genome_readers import read_genome
from genome_readers import get_kmers

def parse_args():
    # arguments: input_file (fasta), output_file (txt)
    parser = argparse.ArgumentParser(description="Compute the substitution rate")
    parser.add_argument("input_file", type=str, help="Input file (fasta)")
    parser.add_argument("output_file", type=str, help="Output file (txt)")
    parser.add_argument("--k", type=int, default=21, help="k-mer size")
    parser.add_argument("--num_sim", type=int, default=5, help="Number of simulations")
    return parser.parse_args()

def compute_subst_rate(kmers_orig, kmers_mut, L1, L2):
    #L = (L1+L2)/2
    L = L1
    num_intersection = len(set(kmers_orig).intersection(set(kmers_mut)))
    #num_intersection = (L - k + 1) * (1 - ps)^k -- for the simple mutation model
    one_minus_ps = (num_intersection / (L - k + 1))**(1.0/k)
    return 1 - one_minus_ps

if __name__ == '__main__':
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file
    k = args.k
    num_sim = args.num_sim

    genome_string = read_genome(args.input_file)
    L1 = len(genome_string)
    kmers_orig = get_kmers(genome_string, k)
    genome_file_prefix = input_file.split(".")[0]

    mut_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    for ps, pd, d in product(mut_rates, repeat=3):
        for i in range(num_sim):
            mutated_filename = genome_file_prefix + "_mutated_" + str(ps) + "_" + str(pd) + "_" + str(d) + "_" + str(i) + ".fasta"
            mutated_genome_string = read_genome(mutated_filename)
            L2 = len(mutated_genome_string)
            kmers_mut = get_kmers(mutated_genome_string, k)

            # compute the substitution rate
            subst_rate = compute_subst_rate(kmers_orig, kmers_mut, L1, L2)
            with open(output_file, "a") as f:
                f.write(str(ps) + "\t" + str(pd) + "\t" + str(d) + "\t" + str(i) + "\t" + str(subst_rate) + "\n")
            print(ps, pd, d, i, subst_rate)
