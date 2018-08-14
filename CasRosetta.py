"""This file runs the designed Rosetta protocol for a group of PDB files, each with different mutations."""

import os
from RosettaSub import RosettaSubprocess


class RosettaBatch:

    def __init__(self):
        # "/home/trinhlab/Documents/RosettaCRISPR/"
        self.base_dir = "/Users/brianmendoza/Desktop/RosettaCRISPR/"
        self.output_dir = "FULL_MUT_PDBs/"

        # KEY: structure  VALUE: tuple, first is the protein file, second is the RNA directory, third is the list of
        # DNA directories
        self.pdbs_and_dirs = {"4UN3/": ("ChainB.pdb", "ChainA_MUT/", ["ChainC_MUT/", "ChainD_MUT"]),
                              "4UN4/": ("ChainB.pdb", "ChainA_MUT/", ["ChainC_MUT/", "ChainD_MUT", "ChainE_MUT"]),
                              "4UN5/": ("ChainB.pdb", "ChainA_MUT/", ["ChainC_MUT/", "ChainD_MUT", "ChainE_MUT"]),
                              "5FQ5/": ("ChainB.pdb", "ChainA_MUT/", ["ChainC_MUT/", "ChainD_MUT", "ChainE_MUT"]),
                              "4OO8ABC/": ("ChainA.pdb", "ChainB_MUT/", ["ChainC_MUT/"]),
                              }

        self.off_combos = {1957: (1900, 1956, 0, 0),
                           2015: (1958, 2014, 0, 0),
                           1899: (1842, 1898, 0, 0),
                           2073: (2016, 2072, 0, 0),
                           2074: (2075, 2084, 0, 0),
                           2090: (2091, 2102, 0, 0),
                           2103: (2104, 2136, 0, 0),
                           2137: (2138, 2143, 2252, 2257),
                           2144: (2145, 2162, 0, 0),
                           2163: (2164, 2202, 0,0),
                           2203: (2204, 2209, 2258, 2263),
                           2210: (2211, 2222, 0, 0),
                           2223: (2224, 2241, 2264, 2272),
                           2242: (2243, 2251, 2273, 2281),
                           2282: (2283, 2297, 0, 0)
                           }

    def stitch_ON(self, all=True, struct=None):
        for crystal in self.pdbs_and_dirs:
            print("Working on on targets for crystal structure: " + crystal)  # Output tracking
            # Iterate across every rna file in the crystal directory
            for rna_file in os.listdir(self.base_dir + crystal + self.pdbs_and_dirs[crystal][1]):
                if rna_file.startswith(".DS_Store"):
                    continue
                dna_file_id = "d_" + rna_file[2:]
                self.single_stitch(crystal, rna_file, dna_file_id, off_flag=False)

    def stitch_OFF(self, all=True, struct=None):
        # Loop to go across each crystal
        for crystal in self.pdbs_and_dirs:
            print("Working on off targets for crystal structure: " + crystal)  # Output tracking
            # Loop through all the on-target sequences to generate the rna_files
            for target_sequence in self.off_combos:
                rna_file = "r_00" + str(target_sequence)
                # Loop through the range of id's to get the dna file ids for off targets
                for i in range(self.off_combos[1][0],self.off_combos[1][1]):
                    dna_file_id = "d_00" + str(target_sequence)
                    self.single_stitch(crystal, rna_file, dna_file_id, off_flag=True)
                # If there is a second range of off-target dna id's then go through another loop of this range
                if self.off_combos[1][2] != 0:
                    for i in range(self.off_combos[1][2], self.off_combos[1][3]):
                        dna_file_id = "d_00" + str(target_sequence)
                        self.single_stitch(crystal, rna_file, dna_file_id, off_flag=True)

    # Send in a single structure, r_file, and d_file identification to stitch together
    def single_stitch(self, structure, r_file, d_file, off_flag):
        # initiate the file for the final structure:
        if off_flag:
            out_file = open(self.base_dir + structure + self.output_dir +
                            "OFF_" + r_file[2:] + "_" + d_file[2:] + ".pdb", "w")
        else:
            out_file = open(self.base_dir + structure + self.output_dir + 'ON_' + r_file[2:], "w")
        # Change directory to get to the location of the protein chain file
        os.chdir(self.base_dir + structure)
        # Open and copy the lines from the protein to the output file
        pf = open(self.pdbs_and_dirs[structure][0])
        for line in pf:
            out_file.write(line)
        pf.close()
        # Change directory to get to the location of the rna files
        os.chdir(self.base_dir + structure + self.pdbs_and_dirs[structure][1])
        # Open and copy the lines form the rna file passed through to the output file
        rf = open(r_file)
        for line in rf:
            out_file.write(line)
        rf.close()
        # Iterate through the dna directories:
        for dna_directory in self.pdbs_and_dirs[structure][2]:
            os.chdir(self.base_dir + structure + dna_directory)
            df = open(d_file)
            for line in df:
                out_file.write(line)
            df.close()
        # print("Created file with the RNA " + r_file + " and the DNA " + d_file)
        out_file.close()

    def run_pdbs_in_directory(self, file_flag):
        os.chdir(self.base_dir + self.output_dir)
        for file in os.listdir(os.curdir):
            if file.startswith(file_flag):
                rr = RosettaSubprocess("score_jd2.default.macosclangrelease")
                rr.set_inputs(["-in:file:s", file, "-out:file:scorefile", "raw_score_output/" + file + ".sc"])
                rr.run_process()


# Code Execution
rc = RosettaBatch()
rc.stitch_ON()
#rc.run_pdbs_in_directory("ON_")
