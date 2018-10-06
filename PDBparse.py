"""Has all the necessary functions for parsing PDB files.  When you load a PDB file it will automatically break down
the chains and the name.  Note: It ignores the Pose output from Rosetta, that is in another class/file."""


class PDB:

    def __init__(self, filename):
        # Container for the chains and their strings, none have any more than E chains
        self.Chain = {'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        # Container for the heteroatom chains of individual chains
        self.hetero_hold = {' ': "", 'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        # opens the file and calls the line parser
        self.load_file(filename)

        # DO NOT DELETE ANY SPACES, THESE ARE IMPORTANT FOR PDB FILE STRUCTURE!
        self.ter = "TER                                                                             \n"

    def load_file(self, file):
        f = open(file)
        for line in f:
            # If you get to the pose, end the process
            if line.startswith("#"):
                break
            # Parse the line by adding it to the containers
            self.parse_line(line)
        f.close()

    def parse_line(self, l):
        chain = l[21]
        if l.startswith("ATOM"):
            self.Chain[chain] += l
        elif l.startswith("HETNAM"):
            self.hetero_hold[chain] += l

    def return_chain(self,chain_id):
        return self.Chain[chain_id] + self.ter

"""def ATOM_parse(l):
    name = l[:6]
    serial = int(l[6:11])
    a_name = l[12:16]
    altloc = l[16]
    resName = l[17:20]
    chainID = l[21]
    resSeq = l[22:26]
    iCode = l[26]
    x = l[30:38]
    y = l[38:46]
    z = l[46:54]
    occupancy = l[54:60]"""

#p = PDB("/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/Ensemble_1/4un3_min_relaxed_0001.pdb")
#print(p.return_chain("A"))