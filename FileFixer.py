"""Scripts for fixing files that don't conform to the file structure put out by Rosetta."""

import os, random

def fix_protein_signal_pose(directory,protsignal):
    for file in os.listdir(directory):
        if file.endswith("min.pdb"):
            f = open(directory + "/" + file)
            y = open("tempfile.pdb",'w')
            for line in f:
                if line.startswith(protsignal):
                    y.write(line[:3] + "CtermProtein" + line[3:])
                else:
                    y.write(line)
            f.close()
            y.close()
            os.rename("tempfile.pdb", directory + "/" + file)

def library_corrector():
    f = open("/Users/brianmendoza/Desktop/mylibfile.txt")
    for line in f:
        myint = int(line[:-1])*100
        intrange = random.randint(0,300)
        myint = random.randint(-intrange,intrange) + myint
        print(abs(myint))

def file_weights_fixer(directory):
    for file in os.listdir(directory):
        if file.endswith("relmin_0001.pdb"):
            f = open(directory + "/" + file)
            y = open("tempfile.pdb", 'w')
            while True:
                line = f.readline()
                if line.startswith("#BEGIN"):
                    topline = line[:-1]
                    bottomline = f.readline()
                    y.write(topline+bottomline)
                else:
                    y.write(line)
                if line.startswith("#END_POSE"):
                    break
            f.close()
            y.close()
            os.rename("tempfile.pdb", directory + "/" + file)

def file_weights_fixer2(directory):
    for file in os.listdir(directory):
        if file.endswith("relmin_0001.pdb"):
            f = open(directory + "/" + file)
            y = open("tempfile.pdb", 'w')
            while True:
                line = f.readline()
                if line.startswith("# All scores"):
                    topline = f.readline()
                    sl = topline.find("pdblabel") + 3
                    y.write("#BEGIN_POSE_ENERGIES_TABLE " + topline[:sl] + "\n")
                    y.write(topline[sl:])
                else:
                    y.write(line)
                if line.startswith("#END_POSE"):
                    break
            f.close()
            y.close()
            os.rename("tempfile.pdb", directory + "/" + file)



#fix_protein_signal_pose("/Users/brianmendoza/Desktop/RosettaCRISPR/5XUS/Ensemble_1/OFF_TARGET/r_100018/","VAL_1206")
#library_corrector()
file_weights_fixer2("/Users/brianmendoza/Desktop/RosettaCRISPR/5XUS/Ensemble_1/FULL_MUT_PDBs")