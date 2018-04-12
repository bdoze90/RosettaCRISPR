"""This file contains the SeqMutator class.  It allows one to input a sequence with ."""

from tempfile import mkstemp
from shutil import move
import os
from RosettaSub import RosettaSubprocess


class SeqMutator:

    def __init__(self, nuc_type, inpt, oupt, base):
        self.seq_list = inpt
        self.output_directory = oupt
        self.base_dir = base  # base_pdb should be set to the pdb that is missing the respective DNA or RNA
        self.new_seqs = list()
        self.nuc_type = nuc_type

        #self.read_in_seqs()  # STEP 1: get the RNA/DNA mutation sequences and put them in the SOLO folder
        self.change_chain_name("D") # STEP 2: change the chain name on the files so that they don't conflict when added to the full PDB

    def read_in_seqs(self):
        os.chdir(self.base_dir)
        base_pdb = "compDNAonly.pdb"
        f = open(self.seq_list)
        for line in f:
            # use the line concatenator to set you sequence to the subseq you need for your base PDB
            s = line[:-1].split("\t")  # index 0: sequence_id; index 1: sequence for rosetta
            if self.nuc_type == "DNA":
                seq = self.revcom(s[1][-9:-6]) + 'TGGTATTG'
            else:
                seq = s[1]
            print(seq)
            rr = RosettaSubprocess("rna_thread.default.macosclangrelease")
            rr.set_inputs(["-s", base_pdb, "-seq", seq.lower(), "-o", self.output_directory + s[0] + ".pdb"])
            rr.run_process()
        f.close()

    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, new_chain_char):
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


SeqMutator(nuc_type="DNA", inpt="/Users/brianmendoza/Desktop/RosettaCRISPR/dna_seqs.txt",
           oupt="COMP_DNA_MUT_SOLO/",
           base="/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/")

