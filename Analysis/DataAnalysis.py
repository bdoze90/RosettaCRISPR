"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import scipy.stats, numpy
import os
from Analysis.DataImport import PoseData
import MasterEnvVariables


class Ensemble:

    def __init__(self,number):
        self.numID = number
        self.dna_tot_struct_scores = dict()  # Dictionary with the score as value and the off-target as the key
        self.dna_trunc_struct_scores = dict()  # Dictionary with the off-target as key and a list of the truncation scores as value
        self.rna_tot_struct_scores = dict()
        self.rna_trunc_struct_scores = dict()
        self.revdna_tot_struct_scores = dict()
        self.revrna_tot_struct_scores = dict()

    def add_score(self,totalscoreline):
        putat_score = totalscoreline.split(", ")
        if putat_score[2] == "DNA":
            if putat_score[1].find("trunc") != -1:
                # Check to see if there is already a list at the position:
                index = int(putat_score[1][0:4])
                if index not in self.dna_trunc_struct_scores:
                    self.dna_trunc_struct_scores[index] = [None]*20
                # need to make sure to append the appropriate value to the list in trunc scores
                else:
                    self.dna_trunc_struct_scores[index][int(putat_score[1][putat_score[1].find("_")+1:])-1] = putat_score[3:]
            elif int(putat_score[1]) in MasterEnvVariables.OFFBASES_HSU_SPCAS9:
                self.revdna_tot_struct_scores[int(putat_score[0])] = putat_score[3:]
            else:
                self.dna_tot_struct_scores[int(putat_score[1])] = putat_score[3:]

        else:
            if putat_score[1].find("trunc") != -1:
                # Check to see if there is already a list at the position:
                index = int(putat_score[1][0:4])
                if index not in self.rna_trunc_struct_scores:
                    self.rna_trunc_struct_scores[index] = [None] * 20
                # need to make sure to append the appropriate value to the list in trunc scores
                else:
                    self.rna_trunc_struct_scores[index][
                       int(putat_score[1][putat_score[1].find("_") + 1:])-1] = putat_score[3:]
            elif int(putat_score[1]) in MasterEnvVariables.OFFBASES_HSU_SPCAS9:
                self.revrna_tot_struct_scores[int(putat_score[0])] = putat_score[3:]
            else:
                self.rna_tot_struct_scores[int(putat_score[1])] = putat_score[3:]

    # Get values for all the scores and return a list of integers or floats for the onid in question
    def processed_total_scores(self, moltype):
        retdict = dict()
        if moltype == "RNA":
            mydict = self.rna_tot_struct_scores
        else:
            mydict = self.dna_tot_struct_scores
        for scoreset in mydict:
            retscores = list()
            for subscr in mydict[scoreset]:
                retscores.append(float(subscr))
            retdict[scoreset] = retscores
        return retdict

    def processed_truncation_scores(self, moltype):




    def return_scores(self, moltype, scoretype):
        if moltype == 'RNA':
            if scoretype == "truncation":
                return self.rna_trunc_struct_scores
            else:
                return self.rna_tot_struct_scores
        else:
            if scoretype == 'truncation':
                return self.dna_trunc_struct_scores
            else:
                return self.dna_tot_struct_scores




class Analysis:

    def __init__(self, base_dir,structure):
        # Information on location of data and structure
        self.basedirectory = base_dir
        self.struct = structure

        # Containers for experimental data and the Rosetta generated data
        self.experimental_dict = dict()
        self.Ensemble_data = list()

        # Import the scores and their experimental values
        self.get_exp()
        self.get_rosetta_scores()

    # FUNCTIONS FOR IMPORTING DATA #
    def get_exp(self):
        f = open(self.basedirectory + "HsuScores.txt")
        for line in f:
            expcor = line[:-1].split("\t")
            self.experimental_dict[int(expcor[0])] = float(expcor[1])
        f.close()

    def get_rosetta_scores(self):
        # iterate through all of the ensembles
        for i in range(1,6):
            # Create an ensemble object to get the scores
            E = Ensemble(i)
            f = open(self.basedirectory + self.struct + "Ensemble_" + str(i) + "TrimmedScores.txt")
            for line in f:
                E.add_score(line[:-1])
            f.close()
            self.Ensemble_data.append(E)
    # END IMPORT DATA FUNCTIONS #

    # Get basic statistics for the total scores
    def align_scores_to_exp(self,scr_index,startrange=1800,endrange=3000):
        arrayexp = list()
        arrayscr = list()
        scores_set = self.Ensemble_data[0].processed_total_scores("RNA")
        for offtarget in scores_set:
            # check to make sure the offtarget is in the range specified:
            if startrange < offtarget < endrange:
                # check to make sure that both dictionary accessions do not return a key error
                if offtarget in self.experimental_dict:
                    arrayexp.append(self.experimental_dict[offtarget])
                    arrayscr.append(scores_set[offtarget][scr_index])
        print(scipy.stats.pearsonr(arrayexp,arrayscr))


    # Get the correlations for the truncations
    def correlation_of_truncation(self, scr_index, start=0, end=3000):
        arrayexp = list()
        arraytruncs = list()





A = Analysis("/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/","4UN4")
A.align_scores_to_exp(0,startrange=1899,endrange=1958)

