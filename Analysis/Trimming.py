"""Creates a single file out of the truncated pdbs, deleting all information but the relevant RNA/DNA scores."""

import os
import MasterEnvVariables
from Analysis.DataImport import PoseData


def write_data(myfile, target, mutdir):
    print(myfile)
    retstr = str()
    outid = str(target) + ", " + str(myfile)[:-8]
    P = PoseData(mutdir + myfile)
    outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
    retstr += outstr
    outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
    retstr += outstr
    return retstr

# Iterate over all the off-target directories in the ensemble off target directory:
def trim_total_scores():
    Ensemble = 1
    Structure = "4UN4"
    for i in range(0,5):
        directory = "/Volumes/Seagate_Drive/RosettaCRISPR/" + Structure + "/Ensemble_" + str(Ensemble) + "/OFF_TARGET/"
        outputfile = "/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/" + Structure + "Ensemble_" + str(Ensemble) + "TrimmedTotalScores.txt"
        f = open(outputfile, "w")
        for target in MasterEnvVariables.OFFBASES_HSU_SPCAS9:
            fullmutdir = directory + "ON_00" + str(target) + "/full_mut_pdbs/"
            truncmutdir = fullmutdir + "truncs_from_min/"
            for file in os.listdir(fullmutdir):
                if file.endswith(".pdb"):
                    print(file)
                    did = str(file).find("d")
                    rid = str(file).find("r")
                    dnaid = str(file)[did+2:did+8]
                    rnaid = str(file)[rid+2:rid+8]
                    outid = dnaid + ", " + rnaid
                    print(outid)
                    P = PoseData(fullmutdir + file)
                    outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
                    f.write(outstr)
                    outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
                    f.write(outstr)
            """for file in os.listdir(truncmutdir):
                if file.endswith(".pdb"):
                    f.write(write_data(file, target, truncmutdir))"""
        f.close()
        Ensemble += 1

def fix_4un4_trim():
    mydir = "/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/"
    for file in os.listdir(mydir):
        if file.endswith("TrimmedScores.txt"):
            f = open(mydir + "Ensemble_" + file[13] + "TrimmedTotalScores.txt", 'w')
            y = open(mydir + file)
            for line in y:
                if line.find("trunc") == -1:
                    f.write(line)
            y.close()
            f.close()

fix_4un4_trim()