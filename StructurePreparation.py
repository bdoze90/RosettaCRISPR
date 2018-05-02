"""This file contains the StructurePrep class which takes in a base pdb structure and either minimizes or relaxes it.
The output directories are named as such.  Please make sure that your flags file is named appropriately and the
appropriate directories have been created before running this file."""

import os
from RosettaSub import RosettaSubprocess
from multiprocessing.pool import ThreadPool


class StructurePrep:

    def __init__(self,runtype, base_directory):
        self.base_dir = base_directory

        if runtype == 'relax':
            self.run_relaxation()
        else:
            self.run_minimization("base_pdb")

    def run_relaxation(self):
        run_pool = list()
        tp = ThreadPool(3)
        os.chdir(self.base_dir + "on_target_relaxed/")
        mut_dir = self.base_dir + "FULL_MUT/"
        structures = os.listdir(mut_dir)
        fullstruct = list()
        for structure in structures:
            fullstruct.append(mut_dir + structure)
        structures.clear()
        for i in range(15):
            rr = RosettaSubprocess("relax.default.linuxgccrelease")
            rr.set_inputs(["-s", fullstruct[i], "-out:suffix _relaxed", "@general_relax_flags"])
            run_pool.append(rr)
        for item in run_pool:
            tp.apply_async(item.run_process)
        tp.close()
        tp.join()

    def run_minimization(self, base_pdb):
        os.chdir(self.base_dir)
        cycles = 10  # change this if you want to but this is the standard based on the protocol we developed
        for i in range(cycles):
            rr = RosettaSubprocess("minimize.default.linuxgccrelease")
            rr.set_inputs(["-s", base_pdb, "-run:min_tolerance", "0.00001", "-out:suffix", str(i)])
            rr.run_process()


S = StructurePrep("relax", "/home/trinhlab/Documents/RosettaCRISPR/4UN3/")

