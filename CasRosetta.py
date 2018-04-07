"""This file runs the designed Rosetta protocol for a group of PDB files, each with different mutations."""

import os
from RosettaSub import RosettaSubprocess


class RosettaBatch:

    def __init__(self):
        self.base_dir = "/Users/brianmendoza/Desktop/RosettaCRISPR/4OO8ABC/"
        self.output_dir = "FULL_MUT_PDBs/"

        self.base_file_string = str()

        self.rna_file_string = str()
        self.dna_file_string = str()

    def stitching(self, off=False, dna_id="none", rna_id="none"):
        os.chdir(self.base_dir)
        base_pdb = "PROTonly.pdb"
        rna_dir = "RNA_MUT_SOLO/"
        dna_dir = "DNA_MUT_SOLO/"

        # Create the base_file_string:
        bf = open(base_pdb)
        end_line = ""
        for line in bf:
            if line.find("END") != -1:
                end_line = line
            else:
                self.base_file_string += line

        if off:
            dna_file = "d_" + dna_id + ".pdb"
            rna_file = "r_" + rna_id + ".pdb"

        # Put the files into collective strings for input
        for rna_file in os.listdir(self.base_dir + rna_dir):
            os.chdir(self.base_dir + rna_dir)
            if rna_file.startswith('r_'):
                # get sequence ID for comparison lookup
                seq_id = rna_file[2:-4]
                print("working on sequence " + seq_id)
                cf = open(rna_file)
                for line in cf:
                    self.rna_file_string += line
                cf.close()
                # grabs the dna pdb from the same seq_id
                dna_file = "d_" + seq_id + ".pdb"
                ocf = open(self.base_dir + dna_dir + dna_file)
                for line in ocf:
                    self.dna_file_string += line
                ocf.close()
                # Creates the new PDB in the FULL_MUT_PDBs directory
                nf = open(self.base_dir + self.output_dir + "ON_" + seq_id + ".pdb", "w")
                nf.write(self.base_file_string + self.rna_file_string + self.dna_file_string + end_line)
                self.dna_file_string = ""
                self.rna_file_string = ""

    def run_pdbs_in_directory(self, file_flag):
        os.chdir(self.base_dir + self.output_dir)
        for file in os.listdir(os.curdir):
            if file.startswith(file_flag):
                rr = RosettaSubprocess("score_jd2.default.macosclangrelease")
                rr.set_inputs(["-in:file:s", file, "-out:file:scorefile", "raw_score_output/" + file + ".sc"])
                rr.run_process()


# Code Execution
rc = RosettaBatch()
rc.run_pdbs_in_directory("ON_")
