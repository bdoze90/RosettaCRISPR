"""Makes directories according to the structure assigned"""

import os

offbases = [1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]

# Set the directory to RosettaCRISPR
base_dir = "/Users/brianmendoza/Dropbox/RosettaCRISPR/"


def make_master_dirs():
    os.chdir(base_dir)
    # Make the structure folders:
    structures = ["4UN3/","4UN4/","4UN5/","5FQ5/","4OO8ABC/"]
    for structure in structures:
        os.mkdir(structure)

    # Make ensemble folders:
    for structure in structures:
        os.chdir(base_dir + structure)
        i = 0
        for i in range(5):
            name = "Ensemble_" + str(i + 1)
            os.mkdir(name)

    # Make off_target folders:
    for structure in structures:
        i = 0
        for i in range(5):
            name = "Ensemble_" + str(i + 1) + "/"
            os.chdir(base_dir+structure+name)
            os.mkdir("OFF_TARGET")

    # Make all the off-target base directories:
    for structure in structures:
        i = 0
        for i in range(5):
            name = "Ensemble_" + str(i + 1)
            os.chdir(base_dir+structure+name+"/OFF_TARGET")
            for offbase in offbases:
                os.mkdir("ON_00" + str(offbase))


def make_struct_specific_dirs(struct):
    os.chdir(base_dir + struct)
