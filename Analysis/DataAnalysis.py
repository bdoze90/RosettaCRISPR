"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import sklearn, scipy, numpy


class Analysis:

    def __init__(self, base_dir):
        self.basedirectory = base_dir


    def get_pearson_correlation(self, compareID):



