"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess
import numpy
import math
import re

f = open("/Volumes/Seagate_Drive/RosettaCRISPR/rna_seqs2.txt")
for line in f:
    print(line)

f.close()
