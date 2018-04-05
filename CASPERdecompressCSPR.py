"""File for decompressing .cspr files"""

from SeqTranslate import SeqTranslate

file = "/Users/brianmendoza/Desktop/pdespCas9.cspr"
f = open(file)
y = open("uncompressedpdespCas9.csv","w")

S = SeqTranslate()
# there needs to be a feature here for adding the PAMs by looking up information


in_repeats = False
for line in f:
    if line.find("REPEATS") != -1:
        break
    if line.find("CHROMOSOME") == -1:
        y.write(S.decompress_csf_tuple(line) + "\n")
    else:
        y.write(line)


