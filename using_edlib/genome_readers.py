from Bio import SeqIO

def read_genome(genome_file):
    """
    Reads a genome file and returns the genome as a string
    """
    genome = ""
    for record in SeqIO.parse(genome_file, "fasta"):
        genome += str(record.seq)
    return genome

def clean_genome_string(genome_string):
    """
    Removes all non-alphabet characters from a genome string
    """
    alphabet = set('ACGT')
    return ''.join(filter(alphabet.__contains__, genome_string))
                   
def get_kmers(genome_string, k):
    """
    Returns a list of all k-mers in a genome string
    """
    kmers = []
    for i in range(len(genome_string)-k+1):
        kmers.append(genome_string[i:i+k])
    return kmers

def read_unitigs(unitigs_file):
    unitigs = set()
    with open(unitigs_file) as f:
        for line in f:
            if line[0] == '>':
                continue
            else:
                unitigs.add(line.strip())
    return list(unitigs)


def read_true_SDIN_values(mutated_genome_file):
    with open(mutated_genome_file) as f:
        # read the first line, which looks like <someString>_<S>_<D>_<I>_<N>
        line = f.readline().strip()
        parts = line.split('_')
        S = int(parts[1])
        D = int(parts[2])
        I = int(parts[3])
        N = int(parts[4])
    return S, D, I, N