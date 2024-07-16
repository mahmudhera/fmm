import argparse
from create_random import create_random_genome
from itertools import product
from tqdm import tqdm
from mutate_genome import mutate as mutate_genome
from genome_readers import read_true_SDIN_values
import numpy as np

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

    # argument indicating unitigs file exists (if not, we have to create it using cuttlefish)
    parser.add_argument("--unitigs_exists", action="store_true", help="Indicates that the unitigs file already exists")

    # arguments for mutation rate start, end, and step
    parser.add_argument("--ps_start", type=float, help="The start value of ps", default=0.01)
    parser.add_argument("--ps_end", type=float, help="The end value of ps", default=0.1)
    parser.add_argument("--ps_step", type=float, help="The step value of ps", default=0.01)

    # argument for the number of simulations
    parser.add_argument("--num_simulations", type=int, help="The number of simulations", default=10)

    return parser.parse_args()

def main():
    # parse the arguments
    args = parse_args()

    # if genome file does not exist, create it
    if not args.genome_exists:
        create_random_genome(args.genome_filename, args.length, args.chunk_size, args.num_chunks, args.seed)

    genome_file_prefix = args.genome_filename.split('.')[0]

    mutation_rates = [round(x, 2) for x in list( np.arange(args.ps_start, args.ps_end, args.ps_step) )]
    num_simulations = args.num_simulations

    for ps, pd, d in tqdm(list(product(mutation_rates, repeat=3))):
        for i in range(num_simulations):
            # create the mutated genome file using these mutation rates. args = genome_filename, ps, pd, d, seed, output_filename, k
            mutated_filename = genome_file_prefix + "_mutated_" + str(ps) + "_" + str(pd) + "_" + str(d) + "_" + str(i) + ".fasta"
            mutate_genome(args.genome_filename, ps, pd, d, i, mutated_filename, args.k)

            # read true S D I N values from the mutated genome file
            S, D, I, N = read_true_SDIN_values(mutated_filename)

            # show these values
            print("True S D I N values:", ps, pd, d, i, S, D, I, N)


if __name__ == "__main__":
    main()