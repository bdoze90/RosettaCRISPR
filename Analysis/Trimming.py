"""Creates a single file out of the truncated pdbs, deleting all information but the relevant RNA/DNA scores."""

import os
import MasterEnvVariables
from Analysis.DataImport import PoseData


def write_data(myfile, target, mutdir):
    print(myfile)
    retstr = str()
    outid = str(target) + ", " + str(target)
    P = PoseData(mutdir + myfile)
    outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
    retstr += outstr
    outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
    retstr += outstr
    return retstr

# Iterate over all the off-target directories in the ensemble off target directory:
def trim_total_scores(Structure, offbasegroup):
    Ensemble = 1
    for i in range(0,5):
        directory = "/Users/brianmendoza/Desktop/RosettaCRISPR/" + Structure + "/Ensemble_" + str(Ensemble) + "/OFF_TARGET/"
        outputfile = "/Users/brianmendoza/Dropbox/Rosetta/TrimmedScores/" + Structure + "Ensemble_" + str(Ensemble) + "TrimmedTotalScores.txt"
        f = open(outputfile, "w")
        for target in offbasegroup:
            fullmutdir = directory + str(target) + "/full_mut_pdbs/"
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
            for myfile in os.listdir(truncmutdir):
                if myfile.endswith(".pdb"):
                    f.write(write_data(myfile, target, truncmutdir))
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


structure = "5F9R"
offbaselist = MasterEnvVariables.OFFBASES["lbCas12"]
#trim_total_scores(structure, offbaselist)
mylist = [1,3,4,5]
for i in mylist:
    f = open("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/" + structure + "Ensemble_" + str(i) + "TrimmedBaseTotalScores.txt", 'w')
    for file in os.listdir("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/" + structure + "_on_bases/Ensemble_" + str(i)):
        #print(file)
        if file.endswith("min.pdb"):
            f.write(write_data(file, file[5:9], "/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/" + structure + "_on_bases/Ensemble_" + str(i) + "/"))
    f.close()
