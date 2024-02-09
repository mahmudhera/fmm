"""
In this script, we will generate a random string of length L over the alphabet {A, C, G, T}.
After that, we will write the string to a fasta file.
"""

import argparse
import random
import string
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_string(L):
    """
    Generates a random string of length L over the alphabet {A, C, G, T}
    """
    return ''.join(random.choice('ACGT') for _ in range(L))

def write_genome(genome_string, genome_name):
    """
    Writes a genome string to a fasta file
    """
    record = SeqRecord(Seq(genome_string), id=genome_name, description="")
    SeqIO.write(record, genome_name, "fasta")

if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser(description='Generate a random genome string and write it to a fasta file')
    parser.add_argument('genome_name', type=str, help='the name of the genome file')
    parser.add_argument('L', type=int, help='the length of the genome string')
    args = parser.parse_args()
    genome_name = args.genome_name
    L = args.L

    # generate the genome string
    genome_string = generate_random_string(L)

    # write the genome string to a fasta file
    write_genome(genome_string, genome_name)
