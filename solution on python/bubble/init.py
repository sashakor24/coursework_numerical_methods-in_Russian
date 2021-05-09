import bubble
import numpy as np
R = np.array([0.01, 0.1, 0.2, 0.4])
N = np.array([1000, 10000, 100000, 600000])
eu = np.array([False, True, True, True]) # True/False показывает, нужно ли считать при данном N[i]
rk = np.array([True,True,True, True])
lib = np.array([False,True,True,True])

for R0 in R:
    for i in range(len(N)):
        bubble.do(N[i],R0,eu[i],rk[i],lib[i])
