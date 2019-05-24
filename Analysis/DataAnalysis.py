"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import sklearn, scipy, numpy
import os
from Analysis.DataImport import PoseData


class Analysis:

    def __init__(self, base_dir,structure):
        self.basedirectory = base_dir
        self.struct = structure

    def get_pearson_correlation(self, subset, subsetid, scores):
        array1 = list()
        array2 = list()


    

A = Analysis("/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/","4UN4")



