"""Creates a single file out of the truncated pdbs, deleting all information but the relevant RNA/DNA scores."""

Structure = "4UN4"
Ensemble = 1

# Iterate over all the off-target directories in the ensemble off target directory:
directory = "/Volumes/Seagate_Drive/RosettaCRISPR/" + Structure + "/Ensemble_" + str(i)
fullmutdir = directory
for file in fullmutdir:
    if file.endswith(".pdb"):
        f = open(file)

        f.close()
