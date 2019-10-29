"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess
import numpy
import math
import re
import random

"""replace_nt = {"A":"G","G":"A","C":"T","T":"C"}

f = open("/Users/brianmendoza/Desktop/multiseqsinfo.txt")
for line in f:
    linetable = line[:-1].split("\t")
    basesequence = linetable[0]
    offlocations = linetable[3]
    actualoffseq = ""
    for i in range(len(offlocations)):
        if offlocations[i] == 'x':
            actualoffseq += replace_nt[basesequence[i]]
        else:
            actualoffseq += basesequence[i]
    print(actualoffseq)


f.close()"""
for i in range(1000):
    print(random.randrange(0,20000))
