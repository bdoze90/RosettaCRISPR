"""This file creates the mutant pdbs (including stitching) for off-targets and includes them in the same directory as
the on-target bases (struct_output).  This allows multiple structure outputs to be run and scored."""

from tempfile import mkstemp
from shutil import move
import os
from CasRosetta import RosettaBatch
from RosettaSub import RosettaSingleProcess
from PDBparse import PDB


class OffMutator:

    def __init__(self, base, structureID):
        self.base_dir = base + structureID + "/"  # This should be the structure directory of the pdb file
        self.structureID = structureID

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
        self.cs_dict = {"4UN3":{"ChainA": ('r', 0, 81, '', ''),
                                "ChainB": ('protein','n','n','n','n')
                                "ChainC": ('d', 0, 'n', 'rc', 'TGGTATTG'),
                                "ChainD": ('d', 17, 'n', '', 'TGGTATTG')},

                        "4UN4":{"ChainA": ('r', 0, 81, '', ''),
                                "ChainC": ('d', 17, 'n', 'rc', 'TGGTATTG'),
                                "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                "ChainE": ('d', 0, 17, 'rc', '')},

                        "4UN5":{"ChainA": ('r', 0, 81, '', ''),
                                "ChainC": ('d', 20, 'n', 'rc', 'TGGTATTG'),
                                "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                "ChainE": ('d', 0, 17, 'rc', '')},

                        "5FQ5":{"ChainA": ('r', 0, 81, '', ''),
                                "ChainC": ('d', 19, 'n', 'rc', 'TGGTATTG'),
                                "ChainD": ('d', 18, 'n', '', 'TGGTATTG'),
                                "ChainE": ('d', 0, 17, 'rc', '')},

                        "4OO8ABC":{"ChainB": ('r', 0, 'n', '', ''),
                                   "ChainC": ('d', 0, 'n', 'rc', '')}
                        }

        # First thing: Pull apart all the Off-target relaxed files:
        self.pull_apart()

        # Second thing: Make the mutant chains for all the required off-target mutations
        self.chain_mutations()

        self.RC = RosettaBatch()
        #RC.stitch_OFF(struct=b_pdb)


    # This function pulls apart the individual chains of the perfectly matched pdb and sorts the chains into
    # the correct folder.  It will do this for every ensemble and every ON_ pdb in the OFF_TARGET folder
    def pull_apart(self):
        for i in range(5):
            ensemble_dir = self.base_dir + "Ensemble_" + str(i) + "/OFF_TARGET"
            os.chdir(ensemble_dir)
            # Iterate over all the ON_00xxxx directories:
            for directory in os.listdir(os.curdir):
                # Check to make sure it is not a .DS_Store file
                if directory.startswith("ON_"):
                    os.chdir(ensemble_dir + "/" + directory)
                    # dissect the file:
                    myfile = directory + ".pdb"
                    P = PDB(myfile)
                    for chain in self.cs_dict[self.structureID]:
                        f = open(self.base_dir + "Ensemble_" + str(i) + "/OFF_TARGET/" + directory + "/" + chain + ".pdb")
                        f.write(P.return_chain(chain[5]))
                        f.close()




    # This function takes every file in the structure output and pulls apart its protein lines and gets the new chain
    def dissect_base_files(self, nucleotides):
        # Go through all of the base on-target files for the off-target investigation
        os.chdir(self.input_dir)
        for on_base in os.listdir(self.input_dir):
            base_file = PDB(on_base)

            if nucleotides == 'new':
                # get every mutant in the range listed for the off_combos with the key of the on_base:
                for i in range(self.off_combos[int(on_base[5:9])][0],self.off_combos[int(on_base[5:9])][1]):
                    # identify the sequence needed:
                    'd_00' + str(i)
                    'r_00' + str(i)



                self.create_new_rna_dnas(base_file.return_chain("A"), sequences)
            else:
                self.use_old_nucleotides()

    def create_new_rna_dnas(self, chain, mutations):
        # Create a temporary file for the chain to be turned into a pdb:
        tf = open("temp_pdb_file.pdb",'w')
        tf.writelines(chain)
        tf.close()
        # Run the Rosetta Subprocess for each sequence:
        for seq in mutations:
            rr = RosettaSingleProcess("rna_thread.default.macosclangrelease")
            rr.set_inputs(
                ["-s", "temp_pdb_file.pdb", "-seq", seq.lower(), "-o", output_directory + s[0] + ".pdb"])
            rr.run_process()
        os.remove("temp_pdb_file.pdb")


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
            rr = RosettaSingleProcess("rna_thread.default.linuxgccrelease")
            rr.set_inputs(
                ["-s", base_pdb, "-seq", seq.lower(), "-o", output_directory + s[0] + ".pdb"])
            rr.run_process()
        f.close()
        self.change_chain_name(output_directory)

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


