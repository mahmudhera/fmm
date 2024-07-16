import os

def mutate(genome_filename, ps, pd, d, seed, output_filename, k):
    """
    Mutate the genome by changing each base with the probability of mutation_rate
    """
    # command = ./a.out 0 random.fasta 0.00 0.05 0.07 random_mutated.fasta 21
    command = "./a.out " + str(seed) + " " + genome_filename + " " + str(ps) + " " + str(pd) + " " + str(d) + " " + output_filename + " " + str(k)
    os.system(command)