"""This file contains the SeqMutator class.  It allows one to input a sequence with ."""

from tempfile import mkstemp
from shutil import move
import os
from RosettaSub import RosettaSingleProcess


class SeqMutator:

    def __init__(self, base):
        self.seq_lists = {'r': "rna_seqs.txt", 'd': "dna_seqs.txt"}
        self.base_dir = base  # base_pdb should be set to the pdb that is missing the respective DNA or RNA
        self.new_seqs = list()

        self.cs_dict = {"4UN3/ChainA": ('r', 0, 81, '', ''),
                        "4UN3/ChainC": ('d', 0, 'n', 'rc','TGGTATTG'),
                        "4UN3/ChainD": ('d', 17, 'n', '','TGGTATTG'),
                        "4UN4/ChainA": ('r', 0, 81, '',''),
                        "4UN4/ChainC": ('d', 17, 'n', 'rc','TGGTATTG'),
                        "4UN4/ChainD": ('d', 18, 'n', '','TGGTATTG'),
                        "4UN4/ChainE": ('d', 0, 17, 'rc',''),
                        "4UN5/ChainA": ('r', 0, 81, '',''),
                        "4UN5/ChainC": ('d', 20, 'n', 'rc','TGGTATTG'),
                        "4UN5/ChainD": ('d', 18, 'n', '','TGGTATTG'),
                        "4UN5/ChainE": ('d', 0, 17, 'rc',''),
                        "5FQ5/ChainA": ('r', 0, 81, '',''),
                        "5FQ5/ChainC": ('d', 19, 'n', 'rc','TGGTATTG'),
                        "5FQ5/ChainD": ('d', 18, 'n', '','TGGTATTG'),
                        "5FQ5/ChainE": ('d', 0, 17, 'rc',''),
                        "4OO8ABC/ChainB": ('r', 0, 'n', '',''),
                        "4OO8ABC/ChainC": ('d', 0, 'n', 'rc','')
                        }

        for chain in self.cs_dict:
            self.read_in_seqs(chain)  # STEP 1: get the RNA/DNA mutation sequences and put them in the SOLO folder

    def read_in_seqs(self, chain_name):
        base_pdb = self.base_dir + chain_name + ".pdb"
        output_directory = self.base_dir + chain_name + "_MUT" + "/"
        seq_file = self.seq_lists[self.cs_dict[chain_name][0]]
        f = open(self.base_dir + seq_file)
        for line in f:
            s = line[:-1].split("\t")  # index 0: sequence_id; index 1: sequence for rosetta
            ix = self.cs_dict[chain_name]
            # check for addendum:
            qs = s[1] + ix[4]
            # check if the 'n' exists then just use the first index:
            if ix[2] == 'n':
                seq = qs[ix[1]:]
                # check for revcomness:
                if ix[3] == 'rc':
                    seq = self.revcom(seq)
            else:
                seq = qs[ix[1]:ix[2]]  # grabs second index and the length of the sequence for the structure
                if ix[3] == 'rc':
                    seq = self.revcom(seq)
            print(seq)
            rr = RosettaSingleProcess("rna_thread.default.macosclangrelease")
            rr.set_inputs(
                ["-s", base_pdb, "-seq", seq.lower(), "-o", output_directory + s[0] + ".pdb"])
            rr.run_process()
        f.close()
        self.change_chain_name(output_directory)


    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, out_dir):
        new_chain_char = out_dir[-2]
        os.chdir(out_dir)
        for p_file in os.listdir(os.curdir):
            # Create temp file
            if p_file[:2] == "r_" or p_file[:2] == "d_":
                print("changing file: " + p_file)
                fh, abs_path = mkstemp()
                with os.fdopen(fh, 'w') as new_file:
                    with open(p_file) as old_file:
                        for line in old_file:
                            if line.find("TER") != -1 or line.find("HET") != -1:
                                new_file.write(line)
                            else:
                                new_line = line[:21] + new_chain_char + line[22:]  # chain char at position 21
                                new_file.write(new_line)
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
            rnt = change[nt.upper()]
            if complement:
                retseq += rnt
            else:
                retseq = rnt + retseq
        return retseq


SeqMutator(base="/Users/brianmendoza/Desktop/RosettaCRISPR_Relaxed2/")

