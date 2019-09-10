"""This file runs the process for the Truncation exercise, which evaluates the PDB as an unstable transition-state
process, whereby each nucleotide is added into the structure and scored."""

import RosettaSub
import os

class Score_Truncate:

    def __init__(self, structure, ensemble):
        self.base_dir = "/Volumes/Seagate_Drive/RosettaCRISPR/" + structure + "/Ensemble_" + str(ensemble) + "/OFF_TARGET/"

        self.offbases = [1899, 1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]
        self.off_index_Hsu1 = [(1842, 1899), (1900, 1957), (1958, 2015), (2016, 2073), (2075, 2085), (2091, 2103),
                          (2104, 2137), (2138, 2144),
                          (2145, 2163), (2164, 2203), (2204, 2210), (2211, 2223), (2224, 2242), (2243, 2252),
                          (2283, 2297)]



    def run_truncation_scoring(self):
        R = RosettaSub.RosettaSingleProcess("score_jd2.default.macosclangrelease")
        for item in self.offbases:
            file_dir = self.base_dir + "ON_00" + str(item) + "/full_mut_pdbs/truncs_from_min/"
            for file in os.listdir(file_dir):
                if file.endswith(".pdb"):
                    R.set_inputs(["-in:file:s", file_dir+file, "-out:pdb"])
                    R.run_process()


S = Score_Truncate("4UN4",1)
S.run_truncation_scoring()