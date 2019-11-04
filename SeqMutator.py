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
        print(len(self.rSequences))
        print(len(self.dSequences))

        self.structureID = structure
        self.base_dir = base + structure  # base_pdb should be set to the pdb that is missing the respective DNA or RNA

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

                        # Crystal structure has lost a portion of the sequence and therefore only good for short guides
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

                        # Needs rna_seqs2.txt b/c of different scaffold RNA
                        "5F9R": {"ChainA": ('r', 0, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},


                        # Beginning of the saCas9 structures
                        "5CZZ": {"ChainA": ('r', 0, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5AXW": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        # Beginning of the Cas12 structures
                        "5XUU": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5XUS": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5XUT": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5XH6": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5XH7": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')}
                        }

        for i in range(1,2):
            ensemble_dir = self.base_dir + "/Ensemble_" + str(i) + "/"
            # First thing: Pull apart all the Off-target relaxed files:
            if Step1:
                self.pull_apart(ensemble_dir, i)
            # Second thing: Make the full mutant crystal structures that either need to be scored or need to be minimized
            if Step2:
                self.full_mut_pdbs(ensemble_dir)

                # Gets all the sequences from the RNA and DNA text files and puts them into the folders

    # This function pulls apart the relaxed pdb and puts its chains in the ensemble directory.
    def pull_apart(self, edir, ens):
        mypdb = PDB
        os.chdir(edir)
        for item in os.listdir(os.curdir):
            if item.endswith("_rel.pdb"):
                mypdb = PDB(os.getcwd()+ "/" + item)
        for chain in self.cs_dict[self.structureID]:
            f = open(edir + "/" + chain + ".pdb", 'w')
            f.write(mypdb.return_chain(chain[5]))
            f.close()
            # Now mutate the chain to all of its necessary off-target iterations:
            # make sure you aren't trying to parse the protein:
            if self.cs_dict[self.structureID][chain][0] != "protein":
                self.read_in_seqs(edir, chain)

    def grab_seqs(self, base_directory):
        R = open(base_directory + "rna_seqs_2.txt")
        for line in R:
            self.rSequences.append(line[:-1].split("\t")[1].upper())
        R.close()
        D = open(base_directory + "dna_seqs.txt")
        for line in D:
            self.dSequences.append(line[:-1].split("\t")[1].upper())
        D.close()

    def read_in_seqs(self, edir, chain_name):
        base_pdb = edir + chain_name + ".pdb"
        specs = self.cs_dict[self.structureID][chain_name]
        # check to see whether the chain is DNA or RNA, then set filename and gather appropriate sequence:
        if specs[0] == 'd':
            for i in range(1,len(self.dSequences)):
                sequence = self.dSequences[i]
                num = str(i)
                print(num,sequence)
                while len(num) < 6:
                    num = "0" + num
                outfilename = 'd_' + num + ".pdb"
                if specs[2] == 'n':
                    mysequence = sequence[specs[1]:] + specs[4]
                else:
                    mysequence = sequence[specs[1]:specs[2]] + specs[4]
                if specs[3] == 'rc':
                    mysequence = self.revcom(mysequence)
                self.rna_mutation_run(base_pdb, edir, chain_name, mysequence, outfilename)

        elif specs[0] == 'r':
            for i in range(1,len(self.rSequences)):
                num = str(i)
                sequence = self.rSequences[i]
                print(num,sequence)
                while len(num) < 6:
                    num = "0" + num
                outfilename = 'r_' + num + ".pdb"
                if specs[2] == 'n':
                    mysequence = sequence[specs[1]:] + specs[4]
                else:
                    mysequence = sequence[specs[1]:specs[2]] + specs[4]
                if specs[3] == 'rc':
                    mysequence = self.revcom(mysequence)
                self.rna_mutation_run(base_pdb, edir, chain_name, mysequence, outfilename)


    def rna_mutation_run(self,base, edir, chain_name, myseq, outfile):
        # Run the Rosetta Subprocess for each sequence:
        rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")
        rr.set_inputs(
            ["-s", base, "-seq", myseq.lower(), "-o",
             edir + "/" + chain_name + "_MUT/" + outfile])
        rr.run_process()
        print(outfile)


    def full_mut_pdbs(self, edir):
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

        # Get all the mutant dnas:
        for i in range(1, len(self.rSequences)):
            num = str(i)
            while len(num) < 6:
                num = "0" + num
            dna_pdb_string = ""
            for mydna in dchain:
                f = open(os.getcwd() + "/" + mydna + "_MUT/" + "d_" + num + ".pdb")
                for line in f:
                    dna_pdb_string += line
                f.close()
            rna_pdb_string = str()
            # Load the rna_pdb_string:
            f = open(os.getcwd() + "/" + rchain + "_MUT/" + "r_" + num + ".pdb")
            for line in f:
                rna_pdb_string += line
            f.close()
            full_out_pdb = open(
                os.getcwd() + "/" + "FULL_MUT_PDBs/" + "ON_" + num + ".pdb", 'w')
            full_out_pdb.write(rna_pdb_string)
            full_out_pdb.write(protein_string)
            full_out_pdb.write(dna_pdb_string)
            full_out_pdb.close()
            print("All mutants for " + edir + " created.")


    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, out_dir):
        new_chain_char = out_dir[-5]
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


SeqMutator("/home/trinhlab/Desktop/RosettaCRISPR/","5F9R")


