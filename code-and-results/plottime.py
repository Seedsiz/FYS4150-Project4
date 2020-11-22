### Programme which makes a plot of runtime
# with and without paralellization

import numpy as np
import matplotlib.pyplot as plt
import math

# input files
infile1 = open("./Results/cycles/timeO2.txt")
infile2 = open("./Results/cycles/timeO3.txt")
infile3 = open("./Results/cycles/timeOfast.txt")
infile4 = open("./Results/cycles/timeunp.txt")

# loadtext to get matrices
input1 = np.loadtxt(infile1.readlines())
input2 = np.loadtxt(infile2.readlines())
input3 = np.loadtxt(infile3.readlines())
input4 = np.loadtxt(infile4.readlines())

# transform to log scale with base 10
a1 = np.log10(input1[:,1])
a2 = np.log10(input1[:,0])

b1 = np.log10(input2[:,1])
b2 = np.log10(input2[:,0])

c1 = np.log10(input3[:,1])
c2 = np.log10(input3[:,0])

d1 = np.log10(input4[:,1])
d2 = np.log10(input4[:,0])

# plot
plt.plot(a1, a2, 'rx', label = "-O2")
plt.plot(b1, b2, 'bx', label = "-O3")
plt.plot(c1, c2, 'cx', label = "-Ofast")
plt.plot(d1, d2, 'gx', label = "Unparallellized w/out flags")

plt.xlabel('MC cycles $\log_{10}$(N)')
plt.ylabel('Time $\log_{10}$(t)')
plt.legend()
plt.show()
