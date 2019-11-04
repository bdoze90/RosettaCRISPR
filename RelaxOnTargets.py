"""This protocol relaxes the on-target sequences that need to be used for the off-target assay."""

"""Setting to run minimization can also be used to relax all on-targets if necessary."""

import os
from RosettaSub import RosettaSubprocess

class StructurePrep:

    def __init__(self, runtype, base_directory):
        self.base_dir = base_directory
        self.runtype = runtype

        self.batch = list()

        self.off_base_pdbs = [2015, 1899, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]

        if self.runtype == 'relax':
            for offbaseid in self.off_base_pdbs:
                self.batch.append(self.base_dir + "FULL_MUT_PDBs/ON_00" + str(offbaseid) + ".pdb")
            self.run_relaxation("hello")
        elif self.runtype == "minimize_full_muts":
            # Setting the files to be processed for all relaxed on-targets:
            for offbaseid in self.off_base_pdbs:
                for file in os.listdir(self.base_dir+ "OFF_TARGET/ON_00" + str(offbaseid) + "/full_mut_pdbs"):
                    if file.endswith(".pdb"):
                        self.batch.append(self.base_dir + "OFF_TARGET/ON_00" + str(offbaseid) + "/full_mut_pdbs/" + file)
                self.run_minimization(self.base_dir + "OFF_TARGET/ON_00" + str(offbaseid))
        else:
            # Setting the files to be processed for all relaxed on-targets:
            for file in os.listdir(self.base_dir + "on_target_relaxed/struct_output"):
                if file.endswith(".pdb"):
                    self.batch.append(self.base_dir + "on_target_relaxed/struct_output/" + file)
            self.run_minimization("hello")

    def run_relaxation(self, filler):
        os.chdir(self.base_dir+"on_target_relaxed/")
        rr = RosettaSubprocess("relax.default.linuxgccrelease", max_p=10, process_list=self.batch)
        rr.set_inputs(["-s", filler, "-out:suffix _relaxed", "@general_relax_flags"])
        rr.run_batch()

    def run_minimization(self, subdirectory):
        os.chdir(subdirectory)
        cycles = 1  # change this if you want to but this is the standard based on the protocol we developed
        for i in range(cycles):
            rr = RosettaSubprocess("minimize.default.linuxgccrelease", max_p=16, process_list=self.batch)
            rr.set_inputs(["-s", "hello", "@general_minimize_flags"])
            rr.run_batch()
            self.batch = []

"""for i in range(3,6):
    # Runs all the minimizations for all the off-target groups in the given ensemble
    base_directory = "/home/trinhlab/Desktop/RosettaCRISPR/4UN4/Ensemble_" + str(i) + "/OFF_TARGET/"
    for directory in os.listdir(base_directory):
        if directory.startswith("ON_"):  # to avoid .DS_Store hidden files
            my_rosetta = StructurePrep("minimize", base_directory+directory +"/")"""


#Sm = StructurePrep("minimize_full_muts", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_4/")
Smaq = StructurePrep("minimize_full_muts", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_5/")
#Smaqo = StructurePrep("minimize_full_muts", "/home/trinhlab/Desktop/RosettaCRISPR/4UN3/Ensemble_3/")