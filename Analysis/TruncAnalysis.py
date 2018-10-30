"""Python file for analyzing the truncation graphs and processing them into objects and vectors for machine learning.
NOTE: This file is currently structured for import from excel data charts and only for a single on-target structure.
This will eventually change once the actual data that is appropriate to import is decided and the truncation computation
is coming in in a rapid manner."""

# Data storage structures:
ExperimentalData = dict()
ComputationScores = dict()


def import_data(importfile):
    exp_data = False
    f = open(importfile)
    for line in f:
        if line.startswith("Exp"):
            exp_data = True
        elif exp_data:
            seqID = line.split()[0]
            expval = line.split()[1]
            ExperimentalData[seqID] = expval
        else:
            datapoint = line.split()
            sequenceID = datapoint[0]
            score = datapoint[10:14]
            if sequenceID in ComputationScores:
                for i in range(4):
                    ComputationScores[sequenceID][i] = score[i]
            else:
                ComputationScores[sequenceID] = ([score[0]], [score[1]], [score[2]], [score[3]])
    f.close()


class TruncProcessor:

    def __init__(self,Structures,Expdata):
        self.process_structures(Structures)
        self.ExpData = Expdata

        self.ProcessedStructs = dict()

        # need to create a struct out of this
        self.Parameters = {"StartValue" : 0,
                           "EndValue" : 0,
                           "MaxValue" : 0,
                           "MinValue" : 0,
                           "SigDeltas" : list()}


    # Each column in the structures needs to be processed to find the desired features
    def process_structures(self,Structures):
        # Get each structure and send to the processed structures dictionary object
        for item in Structures:
            # Initialize the processed structures dictionary:
            self.ProcessedStructs[item] = self.Parameters
            # Process each score in the preprocessed dictionary
            for score_vec in Structures[item]:
                self.ProcessedStructs[item]["StartValue"] = score_vec[1]  # 1 because the very first is 0
                self.ProcessedStructs[item]["EndValue"] = score_vec[-1]
                # Max and min need to be processed such that they are not the start or end (hence the shortened vectors)
                self.ProcessedStructs[item]["MaxValue"] = score_vec[2:-1]

















    "/Users/brianmendoza/Dropbox/Rosetta/importfile.txt"