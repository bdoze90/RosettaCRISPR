"""This file lays out the classes for storing the data in appropriate objects for querying in the DataImport
and DataAnalysis files."""

class StructureScores:

    def __init__(self):
        self.structure_name = str()
        self.Ensembles = list()


class Ensemble:

    def __init__(self):
        self.id = int()
        self.Targets = dict()  # dictionary of Target class objects with their id number as their key

    def return_score_matrix(self, onid):
        self.Targets[onid].returnallscores()


class Target:

    def __init__(self):
        self.sequence = str()
        self.on_score_DNA = float()
        self.on_score_RNA = float()
        self.on_score_min_DNA = float()
        self.on_score_min_RNA = float()

    def returnallscores(self):
        for item in 


class Score:

    def __init__(self):
        self.type = str() #DNA or RNA
        self.min = bool()
        self.total = float()

        self.subscores = list()  # see Pose for the position of each subscore

    def ret_subscore(self, index):
        return self.subscores[index]