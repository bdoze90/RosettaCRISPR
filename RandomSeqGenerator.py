"""This file generates a set random DNA sequence for use in modeling"""

import random

no_of_sequences = 10000
sequence_size = 20

dna_nts = ['A', 'T', 'C', 'G']


def create_sequences():
    for s in range(no_of_sequences):
        seq = str()
        for i in range(sequence_size):
            nt = dna_nts[random.randint(0, 3)]
            seq += nt
        print(seq)

# create a function to generate off targets for a set of sequences


create_sequences()
