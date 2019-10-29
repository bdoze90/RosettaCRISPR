"""File takes in PoseData objects from DataImport and analyzes the data for corrleation."""

import scipy.stats, numpy
import os
from Analysis.DataImport import PoseData
import MasterEnvVariables

# ENSEMBLE CLASS:
# Stores the information for a single ensemble.  The number of the ensemble is stored in the numID variable
# Stores information regarding the pose structure scores in dictionaries keyed of the on-target key

class Ensemble:

    def __init__(self,number, struct):
        self.numID = number
        self.structure = struct
        # Dictionary of dictionaries
        # key1: ON ID
        # key2: OFF ID
        # Value: List of lists of truncation scores
        self.dna_off_scores = dict()
        self.rna_off_scores = dict()
        # Baseline Dictionary
        # key1: ON ID
        # Value: List of lists of truncation scores
        self.dna_on_scores = dict()
        self.rna_on_scores = dict()


    # Get all the scores from the four files of the Ensemble
    def gather_scores(self):
        base_dir_file_pre = "/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/" + self.structure + "Ensemble_" + str(self.numID)
        f = open(base_dir_file_pre + "TrimmedTruncationScores.txt")
        for line in f:
            putat_score = line.split(", ")
            if putat_score[2] == "DNA":
                self.add_off_score(self.dna_off_scores, putat_score)
            else:
                self.add_off_score(self.rna_off_scores, putat_score)
        f.close()
        # Gather the full structure score and add it to the end of the off scores dictionary
        y = open(base_dir_file_pre + "TrimmedTotalScores.txt")
        for line in y:
            putat_score = line.split(", ")
            if putat_score[2] == "DNA":
                self.add_off_score(self.dna_off_scores, putat_score, full=True)
            else:
                self.add_off_score(self.rna_off_scores, putat_score, full=True)
        y.close()
        # Gather the base structure truncations and add them to the on scores dictionary
        z = open(base_dir_file_pre + "TrimmedBaseTruncationScores.txt")
        for line in z:
            putat_score = line.split(", ")
            if putat_score[2] == "DNA":
                self.add_on_score(self.dna_on_scores, putat_score)
            else:
                self.add_on_score(self.rna_on_scores, putat_score)
        # Gather the full structure score for the base structure and add them to the on scores dictionary
        t = open(base_dir_file_pre + "TrimmedBaseTotalScores.txt")
        for line in t:
            putat_score = line.split(", ")
            if putat_score[2] == "DNA":
                self.add_on_score(self.dna_on_scores, putat_score, full=True)
            else:
                self.add_on_score(self.rna_on_scores, putat_score, full=True)



    # Sub-function for above to process the decision tree of where to add the score
    def add_off_score(self, score_dict, putat_score, full=False):
        offindex = putat_score[1][0:4]  # this removes the 00 from the front
        if full:
            offindex = str(int(putat_score[1]))
        onindex = str(int(putat_score[0]))
        # Check to see if the ON ID has been initialized:
        if onindex in score_dict:
            # Check to see if the OFF ID has been initialized:
            if offindex not in score_dict[onindex]:
                score_dict[onindex][offindex] = [None] * 20
                score_dict[onindex][offindex][int(putat_score[1][putat_score[1].find("_") + 1:]) - 1] = putat_score[3:]
            # place the list of scores in the correct index
            else:
                if full:
                    score_dict[onindex][offindex].append(putat_score[3:])
                else:
                    score_dict[onindex][offindex][int(putat_score[1][putat_score[1].find("_") + 1:]) - 1] = putat_score[3:]
        elif not full:
            score_dict[onindex] = dict()
            score_dict[onindex][offindex] = [None] * 20
            score_dict[onindex][offindex][int(putat_score[1][putat_score[1].find("_") + 1:]) - 1] = putat_score[3:]

    # Subfunction for above gather_scores function
    def add_on_score(self, score_dict, putat_score, full=False):
        onindex = str(int(putat_score[0]))
        # Check to see if the ON ID has been initialized:
        if onindex in score_dict:
            # Check to see if the OFF ID has been initialized:
            if onindex not in score_dict:
                score_dict[onindex] = [None] * 20  # Not 20 cuz the last comes from the tot scores file
                score_dict[onindex][int(putat_score[1][putat_score[1].find("trunc_") + 6:]) - 1] = putat_score[3:]
            # place the list of scores in the correct index
            else:
                if full:
                    score_dict[onindex].append(putat_score[3:])
                else:
                    score_dict[onindex][int(putat_score[1][putat_score[1].find("trunc_") + 6:]) - 1] = putat_score[3:]
        elif not full:
            score_dict[onindex] = dict()
            score_dict[onindex] = [None] * 20
            score_dict[onindex][int(putat_score[1][putat_score[1].find("trunc_") + 6:]) - 1] = putat_score[3:]



    # Get values for all the scores and return a list of integers or floats for the onid in question
    """def processed_total_scores(self, moltype):
        retdict = dict()
        if moltype == "RNA":
            mydict = self.rna_off_scores
        else:
            mydict = self.dna_off_scores
        for scoreset in mydict:
            retscores = list()
            for subscr in mydict[scoreset]:
                retscores.append(float(subscr))
            retdict[scoreset] = retscores
        return retdict"""

    # This function takes the truncation scores in the selected moltype and returns the absolute change, dG, ddG
    def processed_truncation_scores(self, moltype, onindex, score_index):
        normalized = dict()  # A dictionary with key being the off-target id and value is a list of length 20 with all the dG values
        if moltype == "DNA":
            base_score_matrix = self.dna_on_scores[onindex]
            for off_matrix in self.dna_off_scores[onindex]:
                normalized[off_matrix] = list()  # initializes the list at the site
                for i in range(len(self.dna_off_scores[onindex][off_matrix])):  # this should give you 20 for spCas9
                    normalized[off_matrix].append(float(self.dna_off_scores[onindex][off_matrix][i][score_index]) - float(base_score_matrix[i][score_index]))
        ddG_normalized = dict()
        for item in normalized:
            ddG_normalized[item] = [0]
            for i in range(1,len(normalized[item])):
                ddG_normalized[item].append(normalized[item][i] - normalized[item][i-1])
        """f = open("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/" + onindex + ".txt", 'w')
        for offindex in normalized:
            f.write(offindex)
            for item in normalized[offindex]:
                f.write("\t" + str(item))
            f.write("\n" + offindex)
            for item in ddG_normalized[offindex]:
                f.write("\t" + str(item))
            f.write("\n")
        f.close()"""

        return normalized, ddG_normalized




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
            self.experimental_dict[expcor[0]] = float(expcor[1])
        f.close()

    def get_rosetta_scores(self):
        # iterate through all of the ensembles
        for i in range(4,5):  # need to fix the other files so that the ensemble comes in correctly
            # Create an ensemble object to get the scores
            E = Ensemble(i,"5F9R")
            E.gather_scores()
            self.Ensemble_data.append(E)
    ##### END IMPORT DATA FUNCTIONS #######

    def trunc_information(self,num_ensemble, infotype1, moltype, onindex, infotype2="ddG"):
        unaveraged_data = dict()  # Dictionary containing the return data (ddG, position) keyed by offindex
        score_set = self.Ensemble_data[num_ensemble -1].processed_truncation_scores(moltype, onindex, -1)  # -1 hard-coded for now to get the total
        for offindex in score_set[1]:
            if infotype1 == "max":
                ddGmax = max(score_set[1][offindex])
                dGmax = max(score_set[0][offindex])
                ddGmin = min(score_set[1][offindex])
                unaveraged_data[offindex] = (ddGmax, ddGmin, score_set[1][offindex].index(ddGmax) , score_set[1][offindex].index(ddGmin),
                                             dGmax, score_set[0][offindex][0], score_set[0][offindex].index(dGmax))
        return unaveraged_data





    # Get basic statistics for the total scores
    def align_scores_to_exp(self,infotype, moltype):
        for ontargetgroup in MasterEnvVariables.OFF_COMBOS_HSU_SPCAS9:
            arrayexp = []
            arrayscr = []
            rosetta_matrix = self.trunc_information(infotype, moltype, str(ontargetgroup) )  # This is a dictionary with averaged tuples
            for offtarget in rosetta_matrix:
                if offtarget and str(ontargetgroup) in self.experimental_dict:  #not all are currently in Hsu scores so check to make sure
                    arrayexp.append(self.experimental_dict[offtarget]/self.experimental_dict[str(ontargetgroup)])
                    arrayscr.append(rosetta_matrix[offtarget][0]*rosetta_matrix[offtarget][1])
            print(ontargetgroup, scipy.stats.pearsonr(arrayexp,arrayscr))


    # Get the correlations for the truncations
    def correlation_of_truncation(self, scr_index, start=0, end=3000):
        arrayexp = list()
        arraytruncs = list()





A = Analysis("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/","5F9R")
for item in MasterEnvVariables.OFF_COMBOS_HSU_SPCAS9:
    myprintdict = A.trunc_information(1,"max","DNA",str(item))
    for index in myprintdict:
        print(index, myprintdict[index][0], myprintdict[index][1],
              myprintdict[index][2], myprintdict[index][3], myprintdict[index][4], myprintdict[index][5], myprintdict[index][6])
#A.align_scores_to_exp("max","DNA")

