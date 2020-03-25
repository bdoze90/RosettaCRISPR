"""Has all the necessary functions for parsing PDB files.  When you load a PDB file it will automatically break down
the chains and the name.  Note: It ignores the Pose output from Rosetta, that is in another class/file.
NOTE: This file has the necessary information to make truncations."""

import os

class PDB:

    def __init__(self, directory, filename):
        # Container for the chains and their strings, none have any more than E chains
        self.Chain = {'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        # Container for the heteroatom chains of individual chains
        self.hetero_hold = {' ': "", 'A': "", 'B': "", 'C': "", 'D': "", 'E': ""}

        self.Residues = list()
        self.directory = directory
        self.filename = filename

        # opens the file and calls the line parser
        self.load_file(directory + filename)

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

    # Helper function for the auto_truncate function
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
    def auto_truncate(self,chain_id, length, direction):
        trunc_list = list()
        for i in range(1,length+1):
            self.residue_breakdown(chain_id)
            ret_string = ""
            if direction:
                for item in self.Residues[i:]:
                    ret_string += item
            else:
                for item in self.Residues[:-i]:
                    ret_string += item
            trunc_list.append(ret_string)
            self.Residues = []
        return trunc_list


    # This function reassembles the imported PDB either with a missing chain or with a truncation reassembly
    def reassemble(self, outid, seqlen, trunkchain, trunc_dir=True, removechain="E", truncation=True):
        if truncation:
            trunkdir = self.directory + "truncs_from_min/"
            if not os.path.isdir(trunkdir):
                os.mkdir(trunkdir)
            truncs = self.auto_truncate(trunkchain, seqlen, trunc_dir)
            for i in range(len(truncs)):
                file_string = ""
                for chain in self.Chain:
                    if chain == trunkchain:
                        file_string += truncs[i]
                    elif chain == removechain:
                        poo = 1
                    else:
                        file_string += self.Chain[chain]
                f = open(trunkdir + str(outid) + "trunc_" + str(seqlen-i) + ".pdb","w")
                f.write(file_string + "\n")
                f.close()

    def create_sub_file(self, outdirectory, chains):
        for chain in chains:
            mydir = outdirectory + "Chain" + chain + ".pdb"
            f = open(mydir, 'w')
            f.write(self.hetero_hold[chain])
            f.write(self.Chain[chain])
            f.close()

    def create_na_file(self, outdirectory, chains):
        f = open(outdirectory + self.filename[:-4] + "_NA.pdb", 'w')
        for chain in chains:
            f.write(self.Chain[chain])
            f.write(self.ter)
        f.close()

"""# List of tuples for off-target sequence IDs in the Hsu1 dataset:
offbases = [1899, 1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]
off_index_Hsu1 = [(1842,1899),(1900,1957),(1958,2015),(2016,2073),(2075,2085),(2091,2103),(2104,2137),(2138,2144),
                  (2145,2163),(2164,2203),(2204,2210),(2211,2223),(2224,2242),(2243,2252),(2283,2297)]
# use this for loop for the truncation of the off-target structures
for i in range(0,15):
    ob = str(offbases[i])
    for sid in range(off_index_Hsu1[i][0],off_index_Hsu1[i][1]):
        p = PDB("/Volumes/Seagate_Drive/RosettaCRISPR/4UN4/Ensemble_5/OFF_TARGET/ON_00" + str(ob) + "/full_mut_pdbs/",
                "d_00" + str(ob) + "_r_00" + str(sid) + "_min.pdb")
        p.reassemble(sid)
        print("completed " + ob + str(sid))"""

#for file in os.listdir("/home/trinhlab/Desktop/RosettaCRISPR/5XUS/Ensemble_1/FULL_MUT_PDBs"):
#p = PDB("/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/Ensemble_1/OFF_TARGET/ON_001957/ON_001957_relaxed_0001.pdb")
#p.reassemble("/Users/brianmendoza/Desktop/RosettaCRISPR/4UN3/Ensemble_1/OFF_TARGET/ON_001957/Truncs/",1957)
"""directory = "/Volumes/Seagate_Drive/RosettaCRISPR/4UN3/Ensemble_2/OFF_TARGET/ON_001899/full_mut_pdbs/"
for pdb_file in os.listdir(directory):
    if pdb_file.endswith(".pdb"):
        print(pdb_file)
        p = PDB(directory, pdb_file)
        p.create_na_file("/Volumes/Seagate_Drive/RosettaCRISPR/4UN3/Ensemble_2/OFF_TARGET/ON_001899/full_mut_pdbs/",["A","C","D"])"""