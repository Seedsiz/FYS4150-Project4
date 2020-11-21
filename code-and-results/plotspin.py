### making spinplots ###
"""
Programme to read off txt files and
plot last spin after all cycles for
one temperature
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

Ncycles = int(input("Input number of Monte Carlo cycles:"))
L = int(input("Input number of spins for a given axis:"))
rank = int(input("Input rank number:"))


infile = open("./Results/spinmatrices/spinmatrix" + str(Ncycles) + \
                "-" + str(L) + "by" + str(L) + "rank" + str(rank) +".txt", 'r')

line = infile.readline()
T = line.split()[1]

spinvec = [];

for line in infile:
    numbers = line.split()
    spinvec.append(float(numbers[0]))

infile.close()

spinvec = np.array(spinvec)

spinmatrix = np.zeros((L,L))

for i in range(L): # rows
    for j in range(L): # colums
        spinmatrix[i,j] = spinvec[i*L + j];

colvec = np.linspace(0,L-1,L)
colvec = ["%d" % x for x in colvec]

spindf = pd.DataFrame(spinmatrix, columns = [colvec])

# Input the vector and make a heatmap
plt.title("T:{:.2f}, MC: {:d}".format(float(T),Ncycles), fontsize = 14)
sns.heatmap(spindf,vmin = -1.0, vmax = 1.0, cmap = "coolwarm")
plt.tight_layout()
plt.savefig("./Results/Figures/Spins/spinmatrix-{:d}by{:d}spinsystem-T-{:.2f}.png".format(L,L,float(T)))
plt.show()
