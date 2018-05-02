"""This file contains the SeqMutator class.  It allows one to input a sequence with ."""

from tempfile import mkstemp
from shutil import move
import os
from RosettaSub import RosettaSingleProcess


class SeqMutator:

    def __init__(self, base_pdb, inpt, base):
        self.seq_list = inpt
        self.base_dir = base  # base_pdb should be set to the pdb that is missing the respective DNA or RNA
        self.new_seqs = list()
        self.nuc_type = base_pdb[0:3]
        self.base_pdb = base_pdb
        self.output_directory = self.nuc_type + "_MUT/"

        self.run_pool = list()

        self.read_in_seqs()  # STEP 1: get the RNA/DNA mutation sequences and put them in the SOLO folder
        self.change_chain_name()

    def read_in_seqs(self):
        os.chdir(self.base_dir)
        f = open(self.seq_list)
        for line in f:
            # use the line concatenator to set you sequence to the subseq you need for your base PDB
            s = line[:-1].split("\t")  # index 0: sequence_id; index 1: sequence for rosetta
            if self.nuc_type == "DNA":
                seq = 'CA' + self.revcom(s[1][4:])
            else:
                seq = s[1]
            print(seq)
            rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")
            rr.set_inputs(["-s", self.base_pdb, "-seq", seq.lower(), "-o", self.output_directory + s[0] + ".pdb"])
            rr.run_process()
        f.close()

    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self):
        new_chain_char = self.base_pdb[-5]
        os.chdir(self.base_dir + self.output_directory)
        for p_file in os.listdir(os.curdir):
            # Create temp file
            if p_file[:2] == "r_" or p_file[:2] == "d_":
                print("changing file: " + p_file)
                fh, abs_path = mkstemp()
                with os.fdopen(fh, 'w') as new_file:
                    with open(p_file) as old_file:
                        for line in old_file:
                            if line.find("TER") == -1:
                                new_line = line[:21] + new_chain_char + line[22:]  # chain char at position 21
                                new_file.write(new_line)
                            else:
                                new_file.write(line)
                # Remove original file
                os.remove(p_file)
                # Move new file
                move(abs_path, p_file)

    def revcom(self, sequence, complement=False):
        retseq = ""
        change = {'A': 'T',
                  'T': 'A',
                  'G': 'C',
                  'C': 'G'}
        for nt in sequence:
            rnt = change[nt]
            if complement:
                retseq += rnt
            else:
                retseq = rnt + retseq
        return retseq


SeqMutator(base_pdb="DNA_chainC.pdb", inpt="/home/trinhlab/Documents/RosettaCRISPR/4UN3/dna_seqs.txt",
           base="/home/trinhlab/Documents/RosettaCRISPR/4UN3/")

