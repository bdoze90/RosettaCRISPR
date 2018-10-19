"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import sklearn, scipy, numpy
import os
from Analysis.DataImport import PoseData


def get_dataset(directory,suffix):
    outfile = "/Users/brianmendoza/Desktop/Rosetta_analysisTruncRNA.txt"
    f = open(outfile,'w')
    os.chdir(directory)
    mylabels = True
    for pdb in os.listdir(os.curdir):
        if pdb.endswith(suffix):
            print("Processing ",pdb)
            P = PoseData(os.path.abspath(pdb))
            if mylabels:
                for item in P.labels:
                    f.write(item + " ")
                    mylabels = False
                f.write("\n")
            f.write(pdb + " ")
            for item in P.return_cum_pose_values("RNA"):
                f.write(str(round(item,4)) + " ")
            f.write("\n")
    f.close()



class Analysis:

    def __init__(self, base_dir):
        self.basedirectory = base_dir

    def get_pearson_correlation(self, subset, subsetid, scores):
        array1 = list()
        array2 = list()


    def define_subset(self,subset_type,subset_id):
        if subset_type == 'all':
            return self.basedirectory
        elif subset_type == 'ensemble':
            return self.basedirectory + "/Ensemble_" + str(subset_id)
        elif subset_type == 'onid':
            return self.basedirectory  # need to iterate over ensembles as well


get_dataset("/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/Ensemble_1/OFF_TARGET/ON_001899/full_mut_pdbs/truncs/","scr.pdb")




