"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import sklearn, scipy, numpy


class Analysis:

    def __init__(self, base_dir):
        self.basedirectory = base_dir


    def comparePose(self):

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







