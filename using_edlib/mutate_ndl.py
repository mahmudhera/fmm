from itertools import product
from mutate_genome import mutate
from tqdm import tqdm
from run_cuttlefish import run_cuttlefish

mutation_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
k = 21
for ps, pd, d, seed in tqdm(list(product(mutation_rates, mutation_rates, mutation_rates, range(5)))):
    mutated_filename = "data/ndl_mutated_" + str(ps) + "_" + str(pd) + "_" + str(d) + "_" + str(seed) + ".fasta"
    cuttlefish_prefix = mutated_filename + '_unitigs'
    run_cuttlefish(mutated_filename, k, 64, cuttlefish_prefix)