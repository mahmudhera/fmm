import random

def create_random_genome(genome_filename, length, chunk_size, num_chunks, seed):
    """
    Create a random genome of length length
    """
    random.seed(seed)
    alphabet = "ACGT"
    
    rand_str_list = []
    for i in range(num_chunks):
        rand_str = ''.join( random.choice(alphabet) for j in range(chunk_size) )
        rand_str_list.append(rand_str)
        
    num_samples = length // chunk_size

    rand_str = ""
    for i in range(num_samples):
        index = random.choice(list(range(num_chunks)))
        rand_str += rand_str_list[index]
    
    with open(genome_filename, "w") as f:
        f.write('> random\n')
        f.write(rand_str)