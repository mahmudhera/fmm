from read_genome import read_genome
import argparse
from read_genome import clean_genome_string
from mutation_model_simulator import mutation_model
from read_genome import get_kmers
from compare_sets_of_kmers import get_num_shared_kmers
from compare_sets_of_kmers import get_num_kmers_single_subst_delt_insert
import os

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

    # parse the arguments
    args = parser.parse_args()
    genome_name = args.genome
    num_simulations = args.num_sim
    k = args.k
    only_sim = args.only_sim

    # read the genome
    genome_string = read_genome(genome_name)
    genome_string = clean_genome_string(genome_string)

    # get the length of the genome
    L1 = len(genome_string)

    # vary the mutation rates
    mutation_rates = [0.01, 0.05, 0.1]

    # vary p_s, p_d, d using the mutation rates
    num_completed = 0
    for p_s in mutation_rates:
        for p_d in mutation_rates:
            for d in mutation_rates:
                for seed in range(num_simulations):
                    num_completed += 1
                    print('Completed: ' + str(num_completed) + ' out of ' + str(len(mutation_rates)**3 * num_simulations))
                    
                    # mutated filename format: genome_name_mutated_p_s_p_d_d_seed.fasta
                    # if mutated file already exists, then skip this simulation
                    mutated_filename = genome_name + '_mutated_' + str(p_s) + '_' + str(p_d) + '_' + str(d) + '_' + str(seed) + '.fasta'
                    if not os.path.exists(mutated_filename):
                        # create a mutation model using the original string, the seed, and the mutation rates
                        mm = mutation_model(seed, genome_string, p_s, p_d, d)

                        # mutate the string, get the mutated string and observations: S, D, I, N_sh
                        mutated_string, S, I, D, N_sh = mm.mutate_string()

                        # write the mutated string in a file, filename format: genome_name_mutated_p_s_p_d_d_seed.fasta
                        # in the first line of the fasta, write S, D, I, N_sh. Format: "> mutated_string_S_D_I_N_sh"
                        with open(mutated_filename, 'w') as f:
                            f.write('> ' + mutated_string + '_' + str(S) + '_' + str(D) + '_' + str(I) + '_' + str(N_sh) + '\n')
                            f.write(mutated_string)

                    # read the file for S, D, I, N_sh
                    with open(mutated_filename, 'r') as f:
                        mutated_string = f.readline().strip().split('_')[0][1:]
                        S = int(f.readline().strip().split('_')[1])
                        D = int(f.readline().strip().split('_')[2])
                        I = int(f.readline().strip().split('_')[3])
                        N_sh = int(f.readline().strip().split('_')[4])
                    
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

                    # get k-mers in the mutated string
                    kmers_mutated = get_kmers(mutated_string, k)

                    # get the observations using the two sets of k-mers: S_calc, D_calc, I_calc, N_sh_calc
                    S_calc, D_calc, I_calc = get_num_kmers_single_subst_delt_insert(kmers_orig, kmers_mutated)
                    N_sh_calc = get_num_shared_kmers(kmers_orig, kmers_mutated)

                    # print everything: p_s, p_d, d, L1, L2, S, D, I, N_sh, S_calc, D_calc, I_calc, N_sh_calc
                    print(p_s, p_d, d, L1, L2, S, D, I, N_sh, S_calc, D_calc, I_calc, N_sh_calc)