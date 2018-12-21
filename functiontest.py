"""File to test simple functions"""


import os
from multiprocessing.pool import ThreadPool
import subprocess
import numpy
import math
import re

mystr = "asloiwe4.12asvei3.02awlef8a"
print(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",mystr))
