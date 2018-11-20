"""Python file for analyzing the truncation graphs and processing them into objects and vectors for machine learning.
NOTE: This file is currently structured for import from excel data charts and only for a single on-target structure.
This will eventually change once the actual data that is appropriate to import is decided and the truncation computation
is coming in in a rapid manner."""

import statistics, numpy

# Data storage structures:
ExperimentalData = dict()
ComputationScores = dict()


def import_data(importfile):
    exp_data = False
    f = open(importfile)
    for line in f:
        if line.startswith("EXP"):
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
                    ComputationScores[sequenceID][i].append(float(score[i]))
            else:
                # min ddG, min dG, no min ddG, no min dG
                ComputationScores[sequenceID] = ([float(score[0])], [float(score[1])], [float(score[2])], [float(score[3])])
    f.close()

class Parameters:

    def __init__(self):
        self.StartValues = float()
        self.EndValues = float()
        self.MaxValues = tuple()
        self.MinValues = tuple()
        self.ddGSigDeltas = list()

    def return_array(self):
        ret_arr = list()
        ret_arr.append(self.StartValue)
        ret_arr.append(self.EndValue)
        ret_arr.append(self.MaxValue[0])
        ret_arr.append(self.MinValue[0])
        ret_arr.append(self.MaxValue[0]*self.MaxValue[1])
        ret_arr.append(self.MaxValue[0]/self.MaxValue[1])
        ret_arr.append(self.MinValue[0] * self.MinValue[1])
        ret_arr.append(self.MinValue[0] / self.MinValue[1])
        for item in self.SigDeltas:
            ret_arr.append(item[0])
            ret_arr.append(item[0]*item[1])
        return ret_arr


class TruncProcessor:

    def __init__(self,Structures,Expdata):
        self.ExpData = Expdata

        self.ProcessedStructs = dict()

        self.process_structures(Structures)


    # Each column in the structures needs to be processed to find the desired features
    def process_structures(self,Structures):
        # Get each structure and send to the processed structures dictionary object
        for item in Structures:
            # Initialize the processed structures dictionary:
            self.ProcessedStructs[item] = Parameters()
            # Process each score in the preprocessed dictionary
            for score_vec in Structures[item]:
                self.ProcessedStructs[item].StartValue = score_vec[1]  # 1 because the very first is 0
                self.ProcessedStructs[item].EndValue = score_vec[-1]
                # Max and min need to be processed such that they are not the start or end (hence the shortened vectors)
                self.ProcessedStructs[item].MaxValue = (max(score_vec[2:-1]), numpy.argmax(score_vec[2:-1]))
                self.ProcessedStructs[item].MinValue = (min(score_vec[2:-1]), numpy.argmin(score_vec[2:-1]))
                self.ProcessedStructs[item].SigDeltas = self.zscores(score_vec,2)



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

import_data("/Users/brianmendoza/Desktop/trunc_rawdata.txt")

C = TruncProcessor(ComputationScores,ExperimentalData)
C.all_data()








