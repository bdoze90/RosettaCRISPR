"""This file contains the SeqMutator class.  It creates the base seq mutations that create the FULL_MUT_PDB files from
which on-target data can be pulled and off-target sequences can be relaxed and then generated."""

from tempfile import mkstemp
from shutil import move
import os
from RosettaSub import RosettaSingleProcess
from PDBparse import PDB


class SeqMutator:

    def __init__(self, base, structure,Step1=True,Step2=True):
        self.rSequences = ["EMPTY"]
        self.dSequences = ["EMPTY"]
        self.grab_seqs(base)

        self.structureID = structure
        self.base_dir = base + "/" + structure  # base_pdb should be set to the pdb that is missing the respective DNA or RNA

        # List of all the appropriate indexes for the mutated sequences and crystal structures
        self.cs_dict = {"4UN3": {"ChainA": ('r', 0, 81, '', ''),
                                 "ChainB": ('protein', '', '', '', ''),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGTATTG'),
                                 "ChainD": ('d', 17, 'n', '', 'TGGTATTG')},

                        "4UN4": {"ChainA": ('r', 0, 81, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 17, 'n', 'rc', 'TGGTATTG'),
                                 "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                 "ChainE": ('d', 0, 17, 'rc', '')},

                        "4UN5": {"ChainA": ('r', 0, 81, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 20, 'n', 'rc', 'TGGTATTG'),
                                 "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                 "ChainE": ('d', 0, 17, 'rc', '')},

                        "5FQ5": {"ChainA": ('r', 0, 81, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 19, 'n', 'rc', 'TGGTATTG'),
                                 "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                 "ChainE": ('d', 0, 17, 'rc', '')},

                        "4OO8ABC": {"ChainB": ('r', 0, 'n', '', ''),
                                    "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                    "ChainC": ('d', 0, 'n', 'rc', '')}
                        }

        for i in range(1,6):
            ensemble_dir = self.base_dir + "Ensemble_" + str(i) + "/OFF_TARGET/"
            # First thing: Pull apart all the Off-target relaxed files:
            if Step1:
                self.pull_apart(ensemble_dir, i)
            # Second thing: Make the full mutant crystal structures that either need to be scored or need to be minimized
            if Step2:
                self.full_mut_pdbs(ensemble_dir)

                # Gets all the sequences from the RNA and DNA text files and puts them into the folders

    def grab_seqs(self, base_directory):
        R = open(base_directory + "/rna_seqs.txt")
        for line in R:
            self.rSequences.append(line[:-1].split("\t")[1].upper())
        R.close()
        D = open(base_directory + "/dna_seqs.txt")
        for line in D:
            self.dSequences.append(line[:-1].split("\t")[1].upper())
        D.close()

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

    def full_mut_pdbs(self, edir):
        # Iterate over every "on-target" directory in the OFF_TARGET folder:
        os.chdir(edir)
        # set the chain names for the dna and rna and protein:
        rchain = str()  # there will only be one chain
        dchain = list()  # list of chains
        protein_string = str()
        for chain in self.cs_dict[self.structureID]:
            if self.cs_dict[self.structureID][chain][0] == 'r':
                rchain = chain
            if self.cs_dict[self.structureID][chain][0] == 'd':
                dchain.append(chain)
            # if you find the protein, then create the protein string
            if self.cs_dict[self.structureID][chain][0] == 'protein':
                f = open(os.getcwd() + "/" + chain + ".pdb")
                for line in f:
                    protein_string += line
                f.close()

        rna_pdb_string = str()
        # Load the rna_pdb_string:
        f = open(os.getcwd() + "/" + rchain + ".pdb")
        for line in f:
            rna_pdb_string += line
        f.close()
        dna_pdb_string = ""
        # Load the dna_pdb_string:
        for dna in dchain:
            f = open(os.getcwd() + "/" + dna + ".pdb")
            for line in f:
                rna_pdb_string += line
            f.close()
        full_out_pdb = open(
            os.getcwd() + "/" + "FULL_MUT_PDBs/" + "ON_" + str(i) + ".pdb", 'w')
        full_out_pdb.write(rna_pdb_string)
        full_out_pdb.write(protein_string)
        full_out_pdb.write(dna_pdb_string)
        full_out_pdb.close()
        print("All rna mutants for " + edir + " created.")


    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, out_dir):
        new_chain_char = out_dir[-6]
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


SeqMutator("/Users/brianmendoza/Desktop/RosettaCRISPR/","4UN3",Step1=False)


