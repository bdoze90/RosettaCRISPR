"""This file is used for the ACF computer to run a set of relaxations or minimizations for the on/off targets"""

import RosettaSub

# This time only do 5F9R and the Ensemble 1 and set it to that (Do ensembles 3, 4, and 5 later)
base_dir = "/lustre/haven/proj/UTK0073/RosettaCRISPR/5F9R/Ensemble_1/"
base_list = [1899, 1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]
# For every full-mut-pdb in the list of on targets, generate 10 structures on it with 5 iterations.
my_proc_list = list()
for item in base_list:
    full_on_pdb = base_dir + "FULL_MUT_PDBs/ON_00" + str(item) + ".pdb"
    my_proc_list.append(full_on_pdb)

rr = RosettaSub.RosettaSubprocess("relax.default.linuxgccrelease",16,my_proc_list)
rr.set_inputs(["-s", 'filler', "nstruct", "10", "relax:default_repeats", "5", "-out:path:pdb", base_dir + "on_target_relaxed","-out:suffix", "_rel"])
rr.run_batch()



# Put the results of the relaxations in their off-target folders with the suffix rel and the number.