"""This file contains the SeqMutator class.  It creates the base seq mutations that create the FULL_MUT_PDB files from
which on-target data can be pulled and off-target sequences can be relaxed and then generated."""

from tempfile import mkstemp
from shutil import move
import os
from RosettaSub import RosettaSingleProcess
from PDBparse import PDB


class SeqMutator:

    def __init__(self, base, structure,Step1=True,Step2=True,Step3=True):
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
                        "5F9R": {"ChainA": ('r', 1, 116, '', ''),
                                 "ChainB": ('protein', 'n','n','n','n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},


                        # Beginning of the saCas9 structures
                        "5CZZ": {"ChainA": ('r', 0, 'n', '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TGGCGATTAG'),
                                 "ChainD": ('d', 11, 'n', '', 'TGGCGATTAG')},

                        "5AXW": {"ChainA": ('r', 0, 'n', '', ''),
                                 "ChainB": ('protein', 'n', 'n', 'n', 'n'),
                                 "ChainC": ('d', 0, 'n', 'rc', 'TTGGGTAG'),
                                 "ChainD": ('d', 'n', 'n', '', 'TTGGGTAG')},


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


        for i in range(1,6):
            ensemble_dir = self.base_dir + "/Ensemble_" + str(i)
            # First thing: Pull apart all the Off-target relaxed files:
            if Step1:
                self.read_in_seqs(ensemble_dir)
            # Second thing: Make the full mutant crystal structures that either need to be scored or need to be minimized
            if Step2:
                self.full_mut_pdbs(ensemble_dir, Step3)
            # Third thing: Assemble the NA only crystal structures into the NA_MUT_PDBs folder

                # Gets all the sequences from the RNA and DNA text files and puts them into the folders

    def grab_seqs(self, base_directory):
        R = open(base_directory + "/rna_seqs2.txt")
        for line in R:
            self.rSequences.append(line[:-1].split("\t")[1].upper())
        R.close()
        D = open(base_directory + "/dna_seqs.txt")
        for line in D:
            self.dSequences.append(line[:-1].split("\t")[1].upper())
        D.close()

    # Function: Generates the new mutated nucleic acid pdb and puts the file into the ChainX_MUT directory
    def read_in_seqs(self, edir):
        os.chdir(edir)

        # Generate the subChainPDBs:
        base_pdb = edir + "/" + self.structureID + "_rel.pdb"
        P = PDB(base_pdb)
        for chain in self.cs_dict[self.structureID]:
            f = open(os.getcwd() + "/" + chain + ".pdb",'w')
            f.write(P.return_chain(chain[5]))
            f.close()

        # Generate the library of sequences of subfiles (MUT directories)
        for chain in self.cs_dict[self.structureID]:
            if not self.cs_dict[self.structureID][chain][0] == "protein":
                self.create_new_rna_dnas(chain, edir)

    # Called from above!
    # Function is called from pull_apart and creates the new sub-pdb chain files for each of the desired mutations
    def create_new_rna_dnas(self, chain, working_dir):
        # Iterate across all the sequence ids for on-targets:
        for i in range(1,len(self.rSequences)):
            # Process the extra zeros for the filename:
            outnumber = str(i)
            for x in range(6-len(str(i))):
                outnumber = '0' + outnumber

            mysequence = str()
            outfilename = str()
            # get the tuple responsible for all chain identifiers:
            specs = self.cs_dict[self.structureID][chain]
            # check to see whether the chain is DNA or RNA, then set filename and gather appropriate sequence:
            if specs[0] == 'd':
                outfilename = 'd_' + outnumber + ".pdb"
                if specs[2] == 'n':
                    mysequence = self.dSequences[i][specs[1]:] + specs[4]
                else:
                    mysequence = self.dSequences[i][specs[1]:specs[2]] + specs[4]
            elif specs[0] == 'r':
                outfilename = 'r_' + outnumber + ".pdb"
                if specs[2] == 'n':
                    mysequence = self.rSequences[i][specs[1]:] + specs[4]
                else:
                    mysequence = self.rSequences[i][specs[1]:specs[2]] + specs[4]
            # check to see if you need to run the sequence through revcom algorithm:
            if specs[3] == 'rc':
                mysequence = self.revcom(mysequence)

            # Run the Rosetta Subprocess for each sequence:
            rr = RosettaSingleProcess("rna_thread.default.macosclangrelease")
            rr.set_inputs(
                ["-s", working_dir + "/" + chain + ".pdb", "-seq", mysequence.lower(), "-o",
                 working_dir + "/" + chain + "_MUT/" + outfilename])
            rr.run_process()
            print(outfilename, i)
        self.change_chain_name(working_dir + "/" + chain + "_MUT")

    # Function: Compiles the Chains into a single pdb that can be passed into Rosetta
    def full_mut_pdbs(self, edir, Step3):
        # Iterate over every sequence in the rna/dna folder:
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

        for i in range(1,len(self.rSequences)):
            # Process the extra zeros for the filename:
            outnumber = str(i)
            for x in range(6 - len(str(i))):
                outnumber = '0' + outnumber

            rna_pdb_string = str()
            # Get to the rna directory:
            os.chdir(edir + "/" + rchain + "_MUT/")
            # Load the rna_pdb_string:
            f = open(os.getcwd() + "/" + "r_" + outnumber + ".pdb")
            for line in f:
                rna_pdb_string += line
            f.close()
            dna_pdb_string = ""
            # Load the dna_pdb_string:
            for dna in dchain:
                os.chdir(edir + "/" + dna + "_MUT/")
                f = open(os.getcwd() + "/" + "d_" + outnumber + ".pdb")
                for line in f:
                    dna_pdb_string += line
                f.close()
            full_out_pdb = open(
                edir + "/" + "FULL_MUT_PDBs/" + "ON_" + outnumber + ".pdb", 'w')
            full_out_pdb.write(rna_pdb_string)
            full_out_pdb.write(protein_string)
            full_out_pdb.write(dna_pdb_string)
            full_out_pdb.close()
            if Step3:
                na_out_pdb = open(
                    edir + "/" + "NA_MUT_PDBs/" + "ON_" + outnumber + ".pdb", 'w')
                na_out_pdb.write(rna_pdb_string)
                na_out_pdb.write(dna_pdb_string)
        print("All rna mutants for " + edir + " created.")


    # -------------------------ACCESSORY FUNCTIONS-------------------------------- #

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


SeqMutator("/Volumes/Seagate_Drive/RosettaCRISPR","5F9R")


