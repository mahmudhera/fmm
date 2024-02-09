from read_genome import read_genome
import argparse
from read_genome import clean_genome_string
from mutation_model_simulator import mutation_model
from read_genome import get_kmers
from compare_sets_of_kmers import get_num_kmers_single_subst_delt_insert_shared
import os
import subprocess
import itertools
import time

# in this script, we will read a genome, vary our mutation rates, and then mutated genome using the varying mutation rates.

if __name__ == '__main__':
    # parse arguments

    parser = argparse.ArgumentParser(description='Take a genpme name, read the contents, then mutate the genome.')

    # make genome a required argument
    parser.add_argument('genome', type=str, help='the genome to mutate')

    # add optional arguments: number of simulations to run, kmer size, flag indicating only simulation
    parser.add_argument('--num_sim', type=int, default=100, help='number of simulations to run')
    parser.add_argument('--k', type=int, default=21, help='kmer size')
    parser.add_argument('--only_sim', action='store_true', help='only simulate, do not compare kmers', default=False)
    parser.add_argument('--run_once', action='store_true', help='only run the first iteration', default=False)

    # parse the arguments
    args = parser.parse_args()
    genome_name = args.genome
    num_simulations = args.num_sim
    k = args.k
    only_sim = args.only_sim
    run_once = args.run_once

    # read the genome
    genome_string = read_genome(genome_name)
    genome_string = clean_genome_string(genome_string)

    # get the length of the genome
    L1 = len(genome_string)

    # vary the mutation rates
    mutation_rates = [0.01, 0.05, 0.1]

    # vary p_s, p_d, d using the mutation rates
    num_completed = 0
    for seed, p_d, d, p_s in list( itertools.product(range(num_simulations), mutation_rates, mutation_rates, mutation_rates) ):

        # measure time taken
        start_time = time.time()
        
        # mutated filename format: genome_name_mutated_p_s_p_d_d_seed.fasta
        # if mutated file already exists, then skip this simulation
        mutated_filename = genome_name + '_mutated_' + str(p_s) + '_' + str(p_d) + '_' + str(d) + '_' + str(seed) + '.fasta'
        if not os.path.exists(mutated_filename):
            # mutate the genome file using command: ./mutate_genome seed, genome_filename, p_s, p_d, d, output_filename, kmer_size
            command = './mutate_genome ' + str(seed) + ' ' + genome_name + ' ' + str(p_s) + ' ' + str(p_d) + ' ' + str(d) + ' ' + mutated_filename + ' ' + str(k)
            print (command)
            subprocess.call(command, shell=True)

        # read the file for S, D, I, N_sh
        with open(mutated_filename, 'r') as f:
            first_line = f.readline()
            S, D, I, N_sh = first_line[1:].strip().split('_')[1:]
            S, D, I, N_sh = int(S), int(D), int(I), int(N_sh)
        
        # do not process any further if only_sim is True
        if only_sim:
            continue

        # get the mutated string
        mutated_string = read_genome(mutated_filename)
        mutated_string = clean_genome_string(mutated_string)

        # L2 = length of the mutated string
        L2 = len(mutated_string)

        # get k-mers in the original string
        kmers_orig = get_kmers(genome_string, k)
        long_kmers_orig = get_kmers(genome_string, k+1)

        kmers_orig_to_count = {}
        for kmer in kmers_orig:
            if kmer in kmers_orig_to_count:
                kmers_orig_to_count[kmer] += 1
            else:
                kmers_orig_to_count[kmer] = 1

        # get k-mers in the mutated string
        kmers_mutated = get_kmers(mutated_string, k)
        long_kmers_mutated = get_kmers(mutated_string, k+1)

        # test num of shared kmers by two ways
        shared_kmers_set = set(kmers_orig).intersection(set(kmers_mutated))
        num_shared_kmers = sum([ kmers_orig_to_count[kmer] for kmer in shared_kmers_set ])
        print(len(num_shared_kmers), N_sh)
        continue

        # get the observations using the two sets of k-mers: S_calc, D_calc, I_calc, N_sh_calc
        S_calc, D_calc, I_calc, N_sh_calc = get_num_kmers_single_subst_delt_insert_shared(kmers_orig, kmers_mutated, long_kmers_orig, long_kmers_mutated)

        # print everything: p_s, p_d, d, L1, L2, S, D, I, N_sh, S_calc, D_calc, I_calc, N_sh_calc
        print(p_s, p_d, d, S, S_calc, D, D_calc, I, I_calc, N_sh, N_sh_calc)

        # increment num_completed
        num_completed += 1
        print('Completed: ' + str(num_completed) + ' out of ' + str(len(mutation_rates)**3 * num_simulations))

        # print time taken
        print('Time taken: ' + str(time.time() - start_time))

        # break if run_once is True
        if run_once:
            break