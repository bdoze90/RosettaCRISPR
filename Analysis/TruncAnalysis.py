"""Python file for analyzing the truncation graphs and processing them into objects and vectors for machine learning.
NOTE: This file is currently structured for import from excel data charts and only for a single on-target structure.
This will eventually change once the actual data that is appropriate to import is decided and the truncation computation
is coming in in a rapid manner."""

import statistics, numpy

# Data storage structures where each item is a TruncStruct:
Structures = dict()


def import_data(importdir, struct, ensemble):
    f = open(importdir +struct+ ensemble + "TrimmedBaseTruncationScores.txt")
    for line in f:
        if line.startswith("EXP"):
            T = TruncStruct()
            T.experimental_value = float(line.split()[2])
            T.ID = line.split()[1]
            Structures[T.ID] = T
        else:
            datapoint = line.split()
            sequenceID = datapoint[0]
            score = datapoint[10:14]
            # min ddG, min dG, no min ddG, no min dG
            Structures[sequenceID].abs_dG = datapoint[10]
            Structures[sequenceID].min_ddG = datapoint[11]

    f.close()
    y = open(importdir+struct+ensemble+"TrimmedTruncationScores.txt")
    for line in f:
        myline = line.split("\t")
        myline
        Structures[int(myline[0])].basestruct[]

class TruncStruct:

    def __init__(self):
        self.ID = str()

        self.abs_dG = list()
        self.min_ddG = list()
        self.min_dG = list()
        self.nomin_dG = list()
        self.nomin_ddg = list()

        self.basestruct = list()

        # Attach the experimental value so that it can be tracked and added when doing analysis
        self.experimental_value = float()

    def import_from_trimmer(self,offid):


    # returns the slope of the
    def get_extreme(self,min_or_max, score_type):
        myscorevec = list()
        if score_type == "ABS":
            myscorevec = self.abs_dG[1:-1]
        elif score_type == "Min_ddG":
            myscorevec = self.min_ddG[1:-1]
        elif score_type == "NoMin_ddG":
            myscorevec = self.nomin_ddg[1:-1]
        elif score_type == "Min_dG":
            myscorevec = self.min_dG[1:-1]
        elif score_type == "NoMin_dG":
            myscorevec = self.nomin_dG[1:-1]

        # Now take the chosen vector and return the selected extrema
        if min_or_max == "MAX":
            return max(myscorevec)
        else:
            return min(myscorevec)

    def get_slope(self,score_type,slope_type):
        return "poop"


class TruncProcessor:

    def __init__(self,Structures,Expdata):
        self.ExpData = Expdata

        self.ProcessedStructs = dict()

        self.process_structures(Structures)



    def zscores(self,vector,threshold):
        return_list = list()
        for i in range(len(vector)):
            zs = (vector[i] - statistics.mean(vector))/statistics.stdev(vector)
            if zs > threshold:
                return_list.append((vector[i],i))
        return return_list

    def all_data(self):
        for item in self.ProcessedStructs:
            print(item)
            print(self.ProcessedStructs[item].return_array())

import_data("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmedScores/","5F9R", "Ensemble_1")


C = TruncProcessor(ComputationScores,ExperimentalData)
C.all_data()








