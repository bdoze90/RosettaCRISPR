"""This file contains the SeqMutator class.  It allows one to input a sequence with ."""

import os, subprocess
import shutil


class SeqMutator:

    def __init__(self, inpt, oupt):
        self.input_directory = inpt
        self.output_directory = oupt
        self.base_pdb = str()  # base_pdb should be set to the pdb that is missing the respective DNA or RNA
        self.new_seqs = list()
        os.chdir(self.input_directory)

    def stitching(self,new_nucleic_acid):
        mut_pdb = "mut_file.pdb"
        shutil.copy2(self.base_pdb, mut_pdb)
        f = open(mut_pdb, 'a')
        f.write(new_nucleic_acid)


    def read_in_seqs(self,seq_file):
        f = open(self.input_directory + seq_file)
        for line in f:
            line[:-1]
    def write_modified_pdbs(self):

