"""Creates a single file out of the truncated pdbs, deleting all information but the relevant RNA/DNA scores."""

import os
import MasterEnvVariables
from Analysis.DataImport import PoseData

Structure = "4UN4"
Ensemble = 1


# Iterate over all the off-target directories in the ensemble off target directory:
for i in range(0,5):
    directory = "/Volumes/Seagate_Drive/RosettaCRISPR/" + Structure + "/Ensemble_" + str(Ensemble) + "/OFF_TARGET/"
    outputfile = "/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/" + Structure + "Ensemble_" + str(Ensemble) + "TrimmedScores.txt"
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
        for file in os.listdir(truncmutdir):
            if file.endswith(".pdb"):
                print(file)
                outid = str(target) + ", " + str(file)[:-8]
                P = PoseData(truncmutdir + file)
                outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
                f.write(outstr)
                outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
                f.write(outstr)
    f.close()
    Ensemble += 1
