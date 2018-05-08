"""This protocol combines other elements of the library to relax a set of structures for on-targeting research"""

import os
from RosettaSub import RosettaSubprocess

class StructurePrep:

    def __init__(self, base_pdb,runtype,base_directory):
        self.base_pdb = base_pdb
        self.base_dir = base_directory
        self.runtype = runtype

        self.batch = list()

        # Setting the files to be processed:
        for file in os.listdir(self.base_dir+"FULL_MUT/"):
            if file.startswith("ON_"):
                self.batch.append(self.base_dir+"FULL_MUT/"+file)

        if self.runtype == 'relax':
            self.run_relaxation("hello")
        else:
            self.run_minimization("hello")

    def run_relaxation(self, filler):
        os.chdir(self.base_dir+"on_target_relaxed/")
        rr = RosettaSubprocess("relax.default.linuxgccrelease", max_p=16, process_list=self.batch)
        rr.set_inputs(["-s", filler, "-out:suffix _relaxed", "@general_relax_flags"])
        rr.run_batch()

    def run_minimization(self, filler):
        os.chdir(self.base_dir+"on_target_relaxed/")
        cycles = 1  # change this if you want to but this is the standard based on the protocol we developed
        for i in range(cycles):
            rr = RosettaSubprocess("minimize.default.linuxgccrelease", max_p=16, process_list=self.batch)
            rr.set_inputs(["-s", filler, "@general_minimize_flags"])
            rr.run_batch()


S = StructurePrep("4un3_min_relaxed_0001.pdb", "minimize", "/home/trinhlab/Documents/RosettaCRISPR/4UN4/")


