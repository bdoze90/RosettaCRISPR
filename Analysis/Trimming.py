"""Creates a single file out of the truncated pdbs, deleting all information but the relevant RNA/DNA scores.
Files generated and their descriptions:
TrimmedTotalScores - Contains the full structure pose information (non-truncated) for all the off-target combinations.
TrimmedBaseTotalScores - Contains the base (on-target) pose data that is non-truncated.
TrimmedBaseTruncationScores - Contains the base (on-target) pose data for all the truncated iterations of the base structure.
TrimmedTruncationScores - Contains the pose information for each truncated structure for each off-target combination."""

import os
import MasterEnvVariables
from Analysis.DataImport import PoseData


def write_data(myfile, mutdir, off=False):
    retstr = str()
    rid = str(myfile).find("r")
    did = str(myfile).find("d")
    if off:  # getting rid of on target contaminants in off-target file
        if rid == did:
            return ""
    outid = myfile[rid+2:rid+8] + ", " + myfile[did+2:did+8] + ", " + myfile[myfile.find("trunc_")+6:myfile.find("_scr")]
    P = PoseData(mutdir + myfile)
    outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
    retstr += outstr
    outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
    retstr += outstr
    return retstr

# Iterate over all the off-target directories in the ensemble off target directory:
def trim_total_scores(Structure, ensemble_num):
    directory = "/home/trinhlab/Desktop/RosettaCRISPR/" + Structure + "/Ensemble_" + str(ensemble_num)
    output_base = "/home/trinhlab/Desktop/RosettaCRISPRTrimmed/" + Structure + "Ensemble_" + str(ensemble_num)
    outputfile1 = output_base + "TrimmedTotalScores.txt"
    outputfile2 = output_base + "TrimmedTruncationScores.txt"
    outputfile3 = output_base + "TrimmedBaseTotalScores.txt"
    outputfile4 = output_base + "TrimmedBaseTruncationScores.txt"

    # This gets the base scores:
    f = open(outputfile3, 'w')
    fullmutdir = directory + "/FULL_MUT_PDBs/"
    for file in os.listdir(fullmutdir):
        if file.endswith("rel_min.pdb"):
            print(file)
            did = str(file).find("d")
            rid = str(file).find("r")
            dnaid = str(file)[did+2:did+8]
            rnaid = str(file)[rid+2:rid+8]
            outid = rnaid + ", " + dnaid
            print(outid)
            P = PoseData(fullmutdir + file)
            outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
            f.write(outstr)
            outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
            f.write(outstr)
    f.close()
    f = open(outputfile4, 'w')
    truncmutdir = fullmutdir + "truncs_from_min/"
    for myfile in os.listdir(truncmutdir):
        if myfile.endswith("scr.pdb"):
            print(myfile)
            f.write(write_data(myfile, truncmutdir))
    f.close()

    # This gets the off-target (non-base) scores
    f = open(outputfile1, 'w')
    fullmutdir = directory + "/OFF_TARGET/"
    for directory in os.listdir(fullmutdir):
        if directory.startswith("r_"):
            for file in os.listdir(fullmutdir + directory):
                if file.endswith("min.pdb"):
                    print(file)
                    did = str(file).find("d")
                    rid = str(file).find("r")
                    dnaid = str(file)[did + 2:did + 8]
                    rnaid = str(file)[rid + 2:rid + 8]
                    if dnaid != rnaid:  # ignores some input contamination that recreates the on target
                        outid = rnaid + ", " + dnaid
                        print(outid)
                        P = PoseData(fullmutdir + directory + "/" + file)
                        outstr = outid + ", DNA, " + str(P.return_cum_pose_values("DNA"))[1:-1] + "\n"
                        f.write(outstr)
                        outstr = outid + ", RNA, " + str(P.return_cum_pose_values("RNA"))[1:-1] + "\n"
                        f.write(outstr)
    f.close()
    f = open(outputfile2, 'w')
    for directory in os.listdir(fullmutdir):
        if directory.startswith("r_"):
            truncmutdir = fullmutdir + directory + "/truncs_from_min/"
            # check to make sure there even exists any off targets:
            if os.path.isdir(truncmutdir):
                print(directory)
                for myfile in os.listdir(truncmutdir):
                    if myfile.endswith("scr.pdb"):
                        f.write(write_data(myfile, truncmutdir,off=True))
    f.close()

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


structure = "5XUS"
mylist = [1,2,3,4,5]  # The ensembles that you want to investigate
for i in mylist:
    trim_total_scores(structure, i)
