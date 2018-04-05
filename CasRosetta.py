"""This file runs the designed Rosetta protocol for a group of PDB files, each with different mutations."""

import os

class MutPDB:

    def __init__(self):
        self.CasType = str()


def run_pdbs_in_directory(directory):
    os.chdir(directory)
    for file in os.listdir(directory):
        # run subprocesses that generate scores for all the files
