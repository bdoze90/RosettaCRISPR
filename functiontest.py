"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess


class TeClass:

    def __init__(self):
        self.blank_item = str()
        self.exec_dir = "/Users/brianmendoza/Desktop/CASPERnative"

    def run_subprocess(self):
        subprocess.run(self.exec_dir)



def testing_threadpool():
    tp = ThreadPool(3)






replace("/Users/brianmendoza/Desktop/testfile.txt")