import argparse
from create_random import create_random_genome
from itertools import product
from tqdm import tqdm
from mutate_genome import mutate as mutate_genome
from genome_readers import read_true_SDIN_values, reverse_complement, read_unitigs, read_genome
import numpy as np
import os
from run_cuttlefish import run_cuttlefish
from multiprocessing import Pool

import edlib



def compute_S_D_I_N(u1, unitig_set_mutd, k):
    
    num_kmers_single_subst, num_kmers_single_delt, num_kmers_no_mutation = 0, 0, 0
    num_kmers_single_insertion = 0
    k = 21

    for u2 in unitig_set_mutd:
        alignment, distance, st1, st2 = None, 9999999999, None, None
        
        r1 = edlib.align(u1, u2, mode = "HW", task = "path")
        r2 = edlib.align(u2, u1, mode = "HW", task = "path")
        
        u3 = reverse_complement(u1)
        r3 = edlib.align(u3, u2, mode = "HW", task = "path")
        r4 = edlib.align(u2, u3, mode = "HW", task = "path")
        
        for i, r in enumerate([r1, r2, r3, r4]):
            if r['editDistance'] < distance:
                alignment, distance = r, r['editDistance']
                if i == 0:
                    st1, st2 = u1, u2
                    flip = False
                elif i == 1:
                    st1, st2 = u2, u1
                    flip = True
                elif i == 2:
                    st1, st2 = u3, u2
                    flip = False
                else:
                    st1, st2 = u2, u3
                    flip = True
        
        nice = edlib.getNiceAlignment(alignment, st1, st2)
        seqA, seqB = nice['query_aligned'], nice['target_aligned']
        assert len(seqA) == len(seqB)
        
        if flip:
            seqB, seqA = seqA, seqB
            
        alphabet = set('ACGT')
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
            if sum(in_numbers[i:i+k]) == 0:
                num_kmers_no_mutation += 1
        
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] in alphabet and seqA[i] == '-':
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_insertion += 1
                    
                    
    return num_kmers_single_subst, num_kmers_single_delt, num_kmers_single_insertion, num_kmers_no_mutation



def compute_S_D_I_N_all(unitig_set_orig, unitig_set_mutd, k, num_threads = 64):
    # call compute_S_D_I_N using a multiprocessing pool
    # return the sum of all the values
    pool = Pool(num_threads)
    arg_list = [(u1, unitig_set_mutd, k) for u1 in unitig_set_orig]
    results = pool.starmap(compute_S_D_I_N, arg_list)
    pool.close()

    S, D, I, N = 0, 0, 0, 0
    for S_, D_, I_, N_ in results:
        S += S_
        D += D_
        I += I_
        N += N_

    return S, D, I, N



def estimate_rates(L, L2, S, D, fA, fA_mut):
    val1 = 3.0 * (fA_mut - 1.0*L2/4.0) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  )
    val2 = 3.0 * (4.0 * S * fA_mut / (4.0*S + 3.0*D) - L2 * S / (4.0*S + 3.0*D) ) / (L - 4.0 * fA)
    assert abs(val1 - val2) < 1e-6

    val1 = 3.0 * (- fA + 1.0*L/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  )
    val2 = 3.0 * S / (4.0*S + 3.0*D)
    assert abs(val1 - val2) < 1e-6

    val1 = 3.0 * (fA_mut - fA + 1.0*L/4 - 1.0*L2/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  )
    val2 = 3.0 * (4.0 * S * fA_mut / (4.0*S + 3.0*D) - L2 * S / (4.0*S + 3.0*D) ) / (L - 4.0 * fA)  + 3.0 * S / (4.0*S + 3.0*D)
    assert abs(val1 - val2) < 1e-6

    subst_rate = 3.0 * (fA_mut - fA + 1.0*L/4 - 1.0*L2/4) / (  (L-4.0*fA) * (1 + 3.0*D/(4.0*S))  )
    del_rate = 1.0 * D * subst_rate / S
    ins_rate = 1.0 * L2 / L - 1.0 + del_rate

    return subst_rate, del_rate, ins_rate



def perform_one_iteration(genome_file_prefix, ps, pd, d, i, args, unitigs_file_orig, L, fA):
    # create the mutated genome file using these mutation rates. args = genome_filename, ps, pd, d, seed, output_filename, k
    mutated_filename = genome_file_prefix + "_mutated_" + str(ps) + "_" + str(pd) + "_" + str(d) + "_" + str(i) + ".fasta"
    
    # if this file already exists, skip this simulation
    if not os.path.exists(mutated_filename):
        mutate_genome(args.genome_filename, ps, pd, d, i, mutated_filename, args.k)

    # read the mutated genome file
    mutated_string = read_genome(mutated_filename)
    L2 = len(mutated_string)
    fA_mut = mutated_string.count('A')

    # read true S D I N values from the mutated genome file
    S, D, I, N = read_true_SDIN_values(mutated_filename)

    mutated_cuttlefish_prefix = mutated_filename + '_unitigs'
    mutated_unitigs_file = mutated_cuttlefish_prefix + '.fa'
    
    # run cuttlefish on the mutated genome file to generate the unitigs file, store the name of the unitigs file
    if not os.path.exists(mutated_unitigs_file):
        run_cuttlefish(mutated_filename, args.k, 64, mutated_cuttlefish_prefix)

    assert os.path.exists(mutated_unitigs_file), f"Mutated unitigs file {mutated_unitigs_file} not found"

    # read two sets of unitigs
    unitigs_orig, unitigs_mut = read_unitigs(unitigs_file_orig), read_unitigs(mutated_unitigs_file)

    # run the alignment based approach to get an estimate of S D I N
    S_est, D_est, I_est, N_est = compute_S_D_I_N_all(unitigs_orig, unitigs_mut, args.k)

    # estimate the mutation rates
    subst_rate, del_rate, ins_rate = estimate_rates(L, L2, S, D, fA, fA_mut)

    # estimate the mutation rates using estimated S and D
    subst_rate_est, del_rate_est, ins_rate_est = estimate_rates(L, L2, S_est, D_est, fA, fA_mut)

    return ps, pd, d, i, S, D, I, N, S_est, D_est, I_est, N_est, subst_rate, del_rate, ins_rate, subst_rate_est, del_rate_est, ins_rate_est



def parse_args():
    parser = argparse.ArgumentParser(description="Check if by using alignments, we can get a good estimate of S D I N")
    parser.add_argument("genome_filename", help="The genome file name")
    parser.add_argument("k", type=int, help="The k-mer size")
    
    # argument indicating genome file exists (if not, we have to create it randomly)
    parser.add_argument("--genome_exists", action="store_true", help="Indicates that the genome file already exists")

    # arguments for creating a random genome
    parser.add_argument("--length", type=int, help="The length of the genome")
    parser.add_argument("--chunk_size", type=int, help="The size of each chunk")
    parser.add_argument("--num_chunks", type=int, help="The number of chunks")
    parser.add_argument("--seed", type=int, help="The seed for the random number generator")

    # arguments for mutation rate start, end, and step
    parser.add_argument("--ps_start", type=float, help="The start value of ps", default=0.01)
    parser.add_argument("--ps_end", type=float, help="The end value of ps", default=0.1)
    parser.add_argument("--ps_step", type=float, help="The step value of ps", default=0.01)

    # argument for the number of simulations
    parser.add_argument("--num_simulations", type=int, help="The number of simulations", default=10)

    # argument for output filename
    parser.add_argument("--output_observations", help="The output filename containing S D I N observations")
    parser.add_argument("--output_rates", help="The output filename containing estimated rates")

    return parser.parse_args()



def main():
    # parse the arguments
    args = parse_args()

    # if genome file does not exist, create it
    if not args.genome_exists:
        create_random_genome(args.genome_filename, args.length, args.chunk_size, args.num_chunks, args.seed)

    cuttlefish_prefix = args.genome_filename + '_unitigs'
    # run cuttlefish on the genome file to generate the unitigs file, args: genome_filename, k, num_threads, outoput_prefix
    run_cuttlefish(args.genome_filename, args.k, 64, cuttlefish_prefix)
    unitigs_file = cuttlefish_prefix + '.fa'

    assert os.path.exists(unitigs_file), f"Unitigs file {unitigs_file} not found"

    mutation_rates = [round(x, 2) for x in list( np.arange(args.ps_start, args.ps_end, args.ps_step) )]
    num_simulations = args.num_simulations
    genome_file_prefix = args.genome_filename.split('.')[0]

    genome_string = read_genome(args.genome_filename)
    L = len(genome_string)
    fA = genome_string.count('A') 

    # open the output files
    f = open(args.output_observations, "w")
    f.write("ps pd d i S D I N S_est D_est I_est N_est\n")
    f2 = open(args.output_rates, "w")
    f2.write("ps pd d i subst_rate del_rate ins_rate subst_rate_est del_rate_est ins_rate_est\n")  

    num_threads = 128
    pool_main = Pool(num_threads)  

    arg_list = []
    for ps, pd, d in product(mutation_rates, repeat=3):
        for i in range(num_simulations):
            arg_list.append((genome_file_prefix, ps, pd, d, i, args, unitigs_file, L, fA))

    results = pool_main.starmap(perform_one_iteration, tqdm(arg_list))
    pool_main.close()

    for res in results:
        ps, pd, d, i, S, D, I, N, S_est, D_est, I_est, N_est, subst_rate, del_rate, ins_rate, subst_rate_est, del_rate_est, ins_rate_est = res
        # write these values to the output file
        f.write(f"{ps} {pd} {d} {i} {S} {D} {I} {N} {S_est} {D_est} {I_est} {N_est}\n")
        print("True S D I N values:", ps, pd, d, i, S, D, I, N, S_est, D_est, I_est, N_est)

        # write these values to the output file
        f2.write(f"{ps} {pd} {d} {i} {subst_rate} {del_rate} {ins_rate} {subst_rate_est} {del_rate_est} {ins_rate_est}\n")
        print("Estimated rates:", ps, pd, d, i, subst_rate, del_rate, ins_rate, subst_rate_est, del_rate_est, ins_rate_est)

        f.flush()
        f2.flush()


if __name__ == "__main__":
    main()