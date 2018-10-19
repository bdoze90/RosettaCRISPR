"""This file runs the process for the Truncation exercise, which evaluates the PDB as an unstable transition-state
process, whereby each nucleotide is added into the structure and scored."""

from PDBparse import PDB

class Truncate:

    def __init__(self, structure, ensemble, targetIDs):
        self.base_dir = "/Users/brianmendoza/Desktop/RosettaCRISPR" + structure + "/Ensemble_" + str(ensemble) + "/OFF_TARGET"
        self.targetIDs = targetIDs



    def min_truncation(self,pdb):
