"""This file runs the designated Rosetta protocol for a group of PDB files."""

import os
from RosettaSub import RosettaSubprocess, RosettaSingleProcess


class RosettaBatch:

    def __init__(self,base,structure,algorithm,batch=True):
        # "/home/trinhlab/Documents/RosettaCRISPR/"
        self.base_dir = base + "/" + structure

        self.batchmode = batch

        # Algorithm function and inputs dictionary:
        self.algorithms = {"Scoring": ["score_jd2.default.macosclangrelease",["-in:file:s",'filler', "-out:pdb"]],
                           "Minimization": ["minimize.default.macosclangrelease",
        ["-s", 'filler', "-out:suffix", "_min", "-run:min_tolerance", "0.001"]],
                           "Relaxation": ["relax.default.linuxgccrelease",
        ["-s", 'filler', "nstruct", "1", "relax:default_repeats", "5", "-out:suffix", "_rel"]]}

        self.selected_algorithm = self.algorithms[algorithm]

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
                                 "ChainB": ('protein', '', '', '', ''),
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

                        "4OO8ABC": {"ChainB": ('r', 0, 'n', '', ''),
                                    "ChainA": ('protein', 'n', 'n', 'n', 'n'),
                                    "ChainC": ('d', 0, 'n', 'rc', '')}
                        }

    def run_ensemble(self,id):
        os.chdir(self.base_dir + "/Ensemble_" + str(id) + "/OFF_TARGET")
        for ontarget in os.listdir(os.curdir):
            if ontarget.startswith("ON_"):  # check to make sure it is not .DS_Store
                os.chdir(os.getcwd()+ "/" + ontarget + "/full_mut_pdbs")
                self.run_pdbs_in_directory(os.getcwd())


    def run_structure(self):
        for ensemble in os.listdir(self.base_dir):
            if ensemble.startswith("Ensemble"):
                self.run_ensemble(int(ensemble[:-1]))


    def run_onTarget(self,targetID,ens_num):
        os.chdir(self.base_dir + "/" + "Ensemble_" + str(ens_num) + "/OFF_TARGET/" + targetID + "/full_mut_pdbs")
        self.run_pdbs_in_directory(os.getcwd())

    def run_truncation(self,targetID,ens_num):
        os.chdir(self.base_dir + "/" + "Ensemble_" + str(ens_num) + "/OFF_TARGET/" + targetID + "/full_mut_pdbs/truncs")
        self.run_pdbs_in_directory(os.getcwd())


    def run_individual(self, pdb_dir_path, pdb_local_name):
        self.run_pdbs_in_directory(pdb_dir_path,singlePDB=pdb_local_name)



    # Main function that calls either single or batch mode in RosettaSub to process the pdb files
    def run_pdbs_in_directory(self, directory, singlePDB=""):
        if singlePDB:
            rr = RosettaSingleProcess(self.selected_algorithm[0])
            self.selected_algorithm[1][1] = singlePDB
            rr.set_inputs(self.selected_algorithm[1])
            rr.run_process()
            return

        pdb_list = list()
        for file in os.listdir(directory):
            if file.endswith(".pdb"):
                pdb_list.append(file)
        if self.batchmode:
            rr = RosettaSubprocess(self.selected_algorithm[0],4,pdb_list)
            # fix this so that the batch is correct
            rr.set_inputs(self.selected_algorithm[1])
            rr.run_batch()
        else:
            for pdb in pdb_list:
                rr = RosettaSingleProcess(self.selected_algorithm[0])
                self.selected_algorithm[1][1] = pdb
                rr.set_inputs(self.selected_algorithm[1])
                rr.run_process()



# Code Execution
rc = RosettaBatch("/Users/brianmendoza/Desktop/RosettaCRISPR","4UN3","Scoring", batch=False)
rc.run_truncation("ON_001899",1)