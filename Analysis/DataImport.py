"""File for getting information from the Pose into Python data structures."""


import os
import pandas, quandl
import sklearn


class PoseData:

    def __init__(self, pdb_name):
        self.mydir = "/Users/brianmendoza/Desktop/RosettaCRISPRMAC/" + pdb_name + "/"
        self.PoseTables = dict()

        self.fill_pose_tables()

    def fill_pose_tables(self):
        os.chdir(self.mydir)


