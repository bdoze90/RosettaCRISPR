"""File for analyzing the ddG data for correlation"""

import numpy as np
import scipy.stats
import statsmodels.api as sm

y = list()
x = [list(),list()]
f = open("/Users/brianmendoza/Dropbox/RosettaCRISPRTrimmed/MachLearn.txt")
for line in f:
    ml = line[:-1].split("\t")
    print(ml)
    y.append(float(ml[0]))
    x[0].append(float(ml[1]))
    x[1].append(float(ml[2]))

def reg_m(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

print(reg_m(y, x).summary())
print(scipy.stats.spearmanr(y,x[1]))