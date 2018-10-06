"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess


os.chdir("/Users/brianmendoza/Desktop/RosettaCRISPR/")
print(os.listdir(os.curdir))
os.chdir(os.curdir + "4UN3")