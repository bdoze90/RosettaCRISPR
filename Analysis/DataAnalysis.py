"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import sklearn, scipy, numpy
import os
from Analysis.DataImport import PoseData


class Analysis:

    def __init__(self, base_dir,structure):
        self.basedirectory = base_dir
        self.struct = structure

        self.experimental_dict = dict()

        # Get the score correlations

    def get_scores(self):
        f = open(self.basedirectory + "HsuScores.txt")
        for line in f:
            expcor = line[:-1].split("\t")
            self.experimental_dict[expcor[0]] = float(expcor[1])
            

    def get_pearson_correlation(self, subset, subsetid, scores):
        array1 = list()
        array2 = list()




A = Analysis("/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/","4UN4")



