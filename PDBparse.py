"""Has all the necessary functions for parsing PDB files.  When you load a PDB file it will automatically break down
the chains and the name.  Note: It ignores the Pose output from Rosetta, that is in another class/file."""


class PDB:

    def __init__(self, filename):
        # Container for the chains and their strings, none have any more than E chains
        self.Chain = {'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        # Container for the heteroatom chains of individual chains
        self.hetero_hold = {' ': "", 'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        self.Residues = list()

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

    def residue_breakdown(self,chain_ID):
        atoms = self.Chain[chain_ID].split("\n")[:-1]
        first_atom_position = int(atoms[0][22:26])
        cur_position = first_atom_position
        curresstring = str()
        for atom in atoms:
            if cur_position == int(atom[22:26]):
                curresstring += atom + "\n"
            else:
                self.Residues.append(curresstring)
                curresstring = atom + "\n"
                cur_position = int(atom[22:26])
        self.Residues.append(curresstring)


    # Returns a list of all the truncated chains (20bp only) in string format
    def auto_truncate(self,chain_id):
        trunc_list = list()
        for i in range(1,21):
            self.residue_breakdown(chain_id)
            ret_string = ""
            for item in self.Residues[:-i]:
                ret_string += item
            trunc_list.append(ret_string)
        return trunc_list


#p = PDB("/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/Ensemble_1/4un3_min_relaxed_0001.pdb")
#p.auto_truncate("A")