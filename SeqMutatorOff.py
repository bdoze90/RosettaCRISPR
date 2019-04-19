"""This file creates the mutant pdbs (including stitching) for off-targets and includes them in the same directory as
the on-target bases (struct_output).  This allows multiple structure outputs to be run and scored."""

from tempfile import mkstemp
import shutil
import os
from RosettaSub import RosettaSingleProcess
from PDBparse import PDB


class OffMutator:

    def __init__(self, base, structureID, NA_stat=False, Step1=True, Step2=True, Step3=False):

        self.base_dir = base + structureID + "/"  # This should be the structure directory of the pdb file
        self.structureID = structureID
        self.NA_only = NA_stat

        # containers for the mutated sequences.  The first position is the "EMPTY" string to get indexing right.
        self.rSequences = ["EMPTY"]
        self.dSequences = ["EMPTY"]

        # List of the off-target sequences needed for the mutation process
        self.off_combos = {1957: (1900, 1956, 0, 0),
                           2015: (1958, 2014, 0, 0),
                           1899: (1842, 1898, 0, 0),
                           2073: (2016, 2072, 0, 0),
                           2074: (2075, 2084, 0, 0),
                           2090: (2091, 2102, 0, 0),
                           2103: (2104, 2136, 0, 0),
                           2137: (2138, 2143, 2252, 2257),
                           2144: (2145, 2162, 0, 0),
                           2163: (2164, 2202, 0, 0),
                           2203: (2204, 2209, 2258, 2263),
                           2210: (2211, 2222, 0, 0),
                           2223: (2224, 2241, 2264, 2272),
                           2242: (2243, 2251, 2273, 2281),
                           2282: (2283, 2297, 0, 0)
                           }

        # List of all the appropriate indexes for the mutated sequences and crystal structures
        self.cs_dict = {"4UN3": {"ChainA": ('r', 0, 81, '', ''),
                                 "ChainB": ('protein','','','',''),
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

                        "4OO8ABC":{"ChainB": ('r', 0, 'n', '', ''),
                                   "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                   "ChainC": ('d', 0, 'n', 'rc', '')}
                        }

        # Before doing anything else: pull in the rna and dna sequences into the RNA and DNA containers (lists):
        self.grab_seqs(base)

        # For each of the ensembles, the following functions need to be performed:
        # self.pull_apart()
        # self.full_mut_pdbs()
        # self.run_Rosetta()
        for i in range(1,3):
            ensemble_dir = self.base_dir + "Ensemble_" + str(i) + "/OFF_TARGET/"
            # First thing: Pull apart all the Off-target relaxed files:
            if Step1:
                self.pull_apart(ensemble_dir,i)
            # Second thing: Make the full mutant crystal structures that either need to be scored or need to be minimized
            if Step2:
                self.full_mut_pdbs(ensemble_dir)
            if Step3:
                self.truncation(ensemble_dir)


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

# --------------------FUNCTIONS CALLED IN THE ENSEMBLE LOOP---------------------------- #

    # FUNCTION # 1
    # This function pulls apart the individual chains of the perfectly matched pdb and sorts the chains into
    # the correct folder.  It will do this for every ensemble and every ON_ pdb in the OFF_TARGET folder
    def pull_apart(self, edir, ens):
        os.chdir(edir)
        # Iterate over all the ON_00xxxx directories:
        for directory in os.listdir(os.curdir):
            # Check to make sure it is not a .DS_Store file
            if directory.startswith("ON_"):
                os.chdir(edir + "/" + directory)
                # get the number of the on-target:
                onid = int(directory[3:])
                # dissect the file:
                myfile =  os.getcwd() + "/" + directory + "_relaxed_000" + str(ens) + ".pdb"
                P = PDB(myfile)
                for chain in self.cs_dict[self.structureID]:
                    wdir = edir + directory
                    f = open(wdir + "/" + chain + ".pdb",'w')
                    f.write(P.return_chain(chain[5]))
                    f.close()
                    # Now mutate the chain to all of its necessary off-target iterations:
                    # make sure you aren't trying to parse the protein:
                    if self.cs_dict[self.structureID][chain][0] != "protein":
                        self.create_new_rna_dnas(onid, chain, wdir)

    # FUNCTION # 1b
    # Function is called from pull_apart and creates the new sub-pdb chain files for each of the desired mutations
    def create_new_rna_dnas(self, ontarget_id, chain, working_dir):
        # Iterate across all the sequence ids for the off-target mutations:
        for i in range(self.off_combos[ontarget_id][0],self.off_combos[ontarget_id][1]+1):  # need the plus 1 to make it inclusive
            mysequence = str()
            outfilename = str()
            # get the tuple responsible for all chain identifiers:
            specs = self.cs_dict[self.structureID][chain]
            # check to see whether the chain is DNA or RNA, then set filename and gather appropriate sequence:
            if specs[0] == 'd':
                outfilename = 'd_00' + str(i) + ".pdb"
                if specs[2] == 'n':
                    mysequence = self.dSequences[i][specs[1]:] + specs[4]
                else:
                    mysequence = self.dSequences[i][specs[1]:specs[2]] + specs[4]
            elif specs[0] == 'r':
                outfilename = 'r_00' + str(i) + ".pdb"
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
                ["-s", working_dir + "/" + chain + ".pdb", "-seq", mysequence.lower(), "-o",  working_dir + "/" + chain + "_MUT/" + outfilename])
            rr.run_process()
            print(outfilename,ontarget_id,i)
        self.change_chain_name(working_dir + "/" + chain + "_MUT")


    # FUNCTION #2
    # This function is the second call in the ensemble loop and creates the mutant structures and the appropriate
    # directory for the structures to be scored in.
    def full_mut_pdbs(self, edir):
        # Iterate over every "on-target" directory in the OFF_TARGET folder:
        os.chdir(edir)
        for directory in os.listdir(os.curdir):
            # Check to make sure it is not a .DS_Store file
            if directory.startswith("ON_"):
                os.chdir(edir + "/" + directory)
            else:
                continue  # skip all directories that aren't on target directories
            # make the full_mut_directory (remove if old one is there):
            if os.path.isdir("full_mut_pdbs"):
                shutil.rmtree("full_mut_pdbs")
            os.mkdir("full_mut_pdbs")

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

            # Create the mutant pdbs by iterating over the off_combos container of the target that matches the directory:
            target = int(directory[3:])
            # Create the rna-on dna-off files:
            rna_pdb_string = str()
            # Load the rna_pdb_string:
            f = open(os.getcwd() + "/" + rchain + ".pdb")
            for line in f:
                rna_pdb_string += line
            f.close()
            # Get all the mutant dnas:
            for i in range(self.off_combos[target][0],self.off_combos[target][1]+1):
                dna_pdb_string = ""
                for mydna in dchain:
                    f = open(os.getcwd() + "/" + mydna + "_MUT/" + "d_00" + str(i) + ".pdb")
                    for line in f:
                        dna_pdb_string += line
                    f.close()
                # Consolidate all the chains into one file:
                full_out_pdb = open(os.getcwd() + "/" + "full_mut_pdbs/" + "r_00" + str(target) + "_d_00" + str(i) + ".pdb",'w')
                full_out_pdb.write(rna_pdb_string)
                full_out_pdb.write(protein_string)
                full_out_pdb.write(dna_pdb_string)
                full_out_pdb.close()
            print("All dna mutants for " + str(target) + " rna created.")

            # Create the rna-off dna-on files:
            dna_pdb_string = ""
            # Load the dna_pdb_string:
            for dna in dchain:
                f = open(os.getcwd() + "/" + dna + ".pdb")
                for line in f:
                    dna_pdb_string += line
                f.close()
            # Get all the mutant rnas:
            for i in range(self.off_combos[target][0], self.off_combos[target][1] + 1):
                rna_pdb_string = ""
                f = open(os.getcwd() + "/" + rchain + "_MUT/" + "r_00" + str(i) + ".pdb")
                for line in f:
                    rna_pdb_string += line
                f.close()
                # Consolidate all the chains into one file:
                full_out_pdb = open(
                    os.getcwd() + "/" + "full_mut_pdbs/" + "d_00" + str(target) + "_r_00" + str(i) + ".pdb", 'w')
                full_out_pdb.write(rna_pdb_string)
                if not self.NA_only:
                    full_out_pdb.write(protein_string)
                full_out_pdb.write(dna_pdb_string)
                full_out_pdb.close()
            print("All rna mutants for " + str(target) + " dna created.")

    # FUNCTION #3
    # This function creates the library of truncated RNAs for the investigation of such scores.
    # Truncated DNAs are more complicated to create and if the truncation trick works, we can look into this




    # -------------------------ACCESSORY FUNCTIONS-------------------------------- #


    # Changes the chain name to be consistent with the new PDB file:
    def change_chain_name(self, out_dir):
        new_chain_char = out_dir[-5]
        os.chdir(out_dir)
        for p_file in os.listdir(os.curdir):
            print(p_file)
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
                shutil.move(abs_path, p_file)

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


O = OffMutator("/Users/brianmendoza/Desktop/RosettaCRISPR/","4UN3",NA_stat=True, Step1=False)