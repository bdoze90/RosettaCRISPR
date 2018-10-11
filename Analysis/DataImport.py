"""File for getting information from the Pose into Python data structures."""


import os

class PoseData:

    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        # All the data on what sequence the PDB has
        self.get_ID_from_filename()
        self.rnaID = int()
        self.dnaID = int()

        # The pose table in somewhat query-able form
        self.PoseTable = dict()
        # Storage for the labels, weights, and full_pose stats
        self.labels = list()
        self.weights = list()
        self.pose = list()

        self.fill_pose_tables()

    # Sets the ID # of the DNA and RNA sequences in the structure
    def get_ID_from_filename(self):
        # Takes care of any base files that find their way into the program
        if self.pdb_path.find("base"):
            self.rnaID = 0
            self.dnaID = 0
        else:
            nameArray = self.pdb_path.split("_")
            if nameArray.startswith('OFF'):
                self.rnaID = int(nameArray[nameArray.index("r") + 1])
                self.dnaID = int(nameArray[nameArray.index("d") + 1])
            else:
                self.rnaID = int(nameArray[1])
                self.dnaID = int(nameArray[1])

    # Function for parsing the Pose table inside the file
    def fill_pose_tables(self):
        in_pose = False
        f = open(self.pdb_path)
        for line in f:
            if in_pose:
                myarray = line.split()
                if myarray[0] == 'label':
                    self.labels = myarray[1:]
                elif myarray[0] == 'weights':
                    for item in myarray[1:-1]:
                        self.weights.append(float(item))
                elif myarray[0] == 'pose':
                    for item in myarray[1:]:
                        self.pose.append(float(item))
                elif line.startswith('#END'):
                    break
                else:
                    self.PoseTable[myarray[0]] = list()
                    for item in myarray[1:]:
                        self.PoseTable[myarray[0]].append(float(item))
            elif line.startswith("#BEGIN"):  # Start of Pose file
                in_pose = True
            else:
                continue
        f.close()

