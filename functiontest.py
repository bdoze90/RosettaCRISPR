"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess


f = open("/Users/brianmendoza/Desktop/ecospCas9.cspr")
while True:
    line = f.readline()
    if line.startswith("REPEAT"):
        break
    sequence = line.split(",")[1]
    if len(sequence) != 10:
        print("broke " + sequence)
f.close()