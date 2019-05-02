"""Makes directories according to the structure assigned"""

import os

offbases = [1899, 1957, 2015, 2073, 2074, 2090, 2103, 2137, 2144, 2163, 2203, 2210, 2223, 2242, 2282]

# Set the directory to RosettaCRISPR
base_dir = "/Volumes/Seagate_Drive/RosettaCRISPR/"
structure = "4UN4"
ensemble = "Ensemble_5"


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

def rename_files_after_processing():
    for off_base in offbases:
        os.chdir(base_dir + structure + "/" + ensemble + "/OFF_TARGET/" + "ON_00" + str(off_base) + "/full_mut_pdbs/truncs_from_min/")
        for pdb in os.listdir(os.curdir):
            if pdb.endswith("0001.pdb") or pdb.endswith("_scr.pdb"):
                os.rename(pdb,pdb[:-9]+"_scr.pdb")
            else:
                os.remove(pdb)
        print("Completed " + str(off_base))

def delete_non_min_files():
    for off_base in offbases:
        os.chdir(base_dir + structure + "/" + ensemble + "/OFF_TARGET/" + "ON_00" + str(off_base) + "/full_mut_pdbs/")
        for pdb in os.listdir(os.curdir):
            if pdb.endswith("0001.pdb"):
                os.rename(pdb,pdb[:-9] + ".pdb")
            elif pdb.endswith(".pdb") and not pdb.endswith("min.pdb"):
                os.remove(pdb)

def make_off_target_dirs(pdbn, off_list):
    for i in range(1,6):
        for off_target in off_list:
            os.mkdir("/Volumes/Seagate_Drive/RosettaCRISPR/" + pdbn + "/Ensemble_" + str(i) + "/OFF_TARGET/" + "ON_00" + str(off_target))

#make_off_target_dirs("5CZZ",offbases)
rename_files_after_processing()
#delete_non_min_files()