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
        self.ProteinPoseTable = dict()
        self.RNAPoseTable = dict()
        self.DNAPoseTable = dict()
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
        f = open(self.pdb_path)
        while True:
            line = f.readline()
            if line.startswith("#BEGIN"):
                # This chunk gets the first three main lines
                self.labels = f.readline().split()[1:]
                self.weights = f.readline().split()[1:-1]
                self.pose = f.readline().split()[1:]
                for i in range(len(self.weights)):
                    self.weights[i] = float(self.weights[i])
                for i in range(len(self.pose)):
                    self.pose[i] = float(self.pose[i])

                # This chunk gets the lines for the actual pose
                while True:
                    line = f.readline()
                    if line.find("LowerRNA") != -1:
                        self.add_line(line,"RNA")
                        while line.find("UpperRNA") == -1:
                            line = f.readline()
                            self.add_line(line,"RNA")
                    if line.find("LowerDNA") != -1:
                        self.add_line(line,"DNA")
                        while line.find("UpperDNA") == -1:
                            line = f.readline()
                            self.add_line(line,"DNA")
                    if line.find("NtermProtein") != -1:
                        self.add_line(line,"Protein")
                        while line.find("CtermProtein") == -1:
                            line = f.readline()
                            self.add_line(line,"Protein")
                    if line.startswith("#END_POSE"):
                        break
                break
            else:
                continue
        f.close()

    #  Moves to the function the ability to take a pure string line and add as floats to the dictionary
    def add_line(self, line, dictionary):
        id = line.split()[0]
        numarray = list()
        for item in line.split()[1:]:
            numarray.append(float(item))
        if dictionary == "RNA":
            self.RNAPoseTable[id] = numarray
        elif dictionary == "DNA":
            self.DNAPoseTable[id] = numarray
        elif dictionary == "Protein":
            self.ProteinPoseTable[id] = numarray

    def return_cum_pose_values(self, pose):
        output = list()
        select_pose = dict()
        if pose == "RNA":
            select_pose = self.RNAPoseTable
        elif pose == "DNA":
            select_pose = self.DNAPoseTable
        elif pose == "Protein":
            select_pose = self.ProteinPoseTable
        for i in range(len(self.labels)):
            tot_sum_score = 0
            for residue in select_pose:
                tot_sum_score += select_pose[residue][i]
            output.append(tot_sum_score)
        return output

