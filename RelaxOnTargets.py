"""This protocol relaxes the on-target sequences that need to be used for the off-target assay."""

"""Setting to run minimization can also be used to relax all on-targets if necessary."""

import os
from RosettaSub import RosettaSubprocess

class StructurePrep:

    def __init__(self, runtype, base_directory):
        self.base_dir = base_directory
        self.runtype = runtype

        self.batch = list()

        self.off_base_pdbs = [1957, 2015, 1899, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]

        if self.runtype == 'relax':
            for offbaseid in self.off_base_pdbs:
                self.batch.append(self.base_dir + "FULL_MUT_PDBs/ON_00" + str(offbaseid) + ".pdb")
            self.run_relaxation("hello")
        else:
            # Setting the files to be processed for all relaxed on-targets:
            for file in os.listdir(self.base_dir + "on_target_relaxed/struct_output/"):
                if file.startswith("ON_"):
                    self.batch.append(self.base_dir + "on_target_relaxed/struct_output/" + file)
            self.run_minimization("hello")

    def run_relaxation(self, filler):
        os.chdir(self.base_dir+"on_target_relaxed/")
        rr = RosettaSubprocess("relax.default.linuxgccrelease", max_p=16, process_list=self.batch)
        rr.set_inputs(["-s", filler, "-out:suffix _relaxed", "@general_relax_flags"])
        rr.run_batch()

    def run_minimization(self, filler):
        os.chdir(self.base_dir+"on_target_relaxed/struct_output/")
        cycles = 1  # change this if you want to but this is the standard based on the protocol we developed
        for i in range(cycles):
            rr = RosettaSubprocess("minimize.default.linuxgccrelease", max_p=16, process_list=self.batch)
            rr.set_inputs(["-s", filler, "@general_minimize_flags"])
            rr.run_batch()


S = StructurePrep("relax", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_1/")
Xp = StructurePrep("relax", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_2/")
Xw = StructurePrep("relax", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_3/")
Soda = StructurePrep("relax", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_4/")
Pop = StructurePrep("relax", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_5/")

Sa = StructurePrep("minimize", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_1/")
Xpa = StructurePrep("minimize", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_2/")
Xwa = StructurePrep("minimize", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_3/")
Sodaa = StructurePrep("minimize", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_4/")
Popam = StructurePrep("minimize", "/home/trinhlab/Desktop/RosettaCRISPR/5F9R/Ensemble_5/")