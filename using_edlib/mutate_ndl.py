from itertools import product
from mutate_genome import mutate
from tqdm import tqdm

mutation_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
k = 21
for ps, pd, d, seed in tqdm(list(product(mutation_rates, mutation_rates, mutation_rates, range(10)))):
    mutated_filename = "data/ndl_mutated_" + str(ps) + "_" + str(pd) + "_" + str(d) + "_" + str(seed) + ".fasta"
    mutate("data/ndl.fasta", ps, pd, d, seed, mutated_filename, k)