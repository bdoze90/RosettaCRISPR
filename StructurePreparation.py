"""This file contains the StructurePrep class which takes in a base pdb structure and either minimizes or relaxes it.
The output directories are named as such.  Please make sure that your flags file is named appropriately and the
appropriate directories have been created before running this file."""

import os
from RosettaSub import RosettaSubprocess


class StructurePrep:

    def __init__(self, base_pdb,runtype,base_directory):
        self.base_pdb = base_pdb
        self.base_dir = base_directory

        if runtype == 'relax':
            self.run_relaxation("hello")
        else:
            self.run_minimization()

    def run_relaxation(self, file):
        os.chdir(self.base_dir)
        rr = RosettaSubprocess("relax.default.linuxgccrelease")
        rr.set_inputs(["-in:file:s", file, "-out:file:scorefile", "raw_score_output/" + file + ".sc"])
        rr.run_process()

    def run_minimization(self):
        os.chdir(self.base_dir)
        cycles = 10  # change this if you want to but this is the standard based on the protocol we developed
        for i in range(cycles):
            rr = RosettaSubprocess("minimize.default.macosclangrelease")
            rr.set_inputs(["-s", self.base_pdb, "-run:min_tolerance", "0.00001", "-out:suffix", str(i)])
            rr.run_process()


S = StructurePrep("4un3.pdb", "minimize", "/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/")

