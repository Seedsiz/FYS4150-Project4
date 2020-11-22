# Programme that make use of the expectation values
# to plot energy, magnetic moment, heat capacity etc
# as a function of temperature

import numpy as np
import matplotlib.pyplot as plt


Ncycles = 100000
L = np.array((60, 80, 100))

rank0 = np.zeros(30)
rank1 = np.zeros(30)
rank2 = np.zeros(30)
rank3 = np.zeros(30)

T = np.zeros(4*30*len(L))
E = np.zeros(4*30*len(L))
Cv = np.zeros(4*30*len(L))
S = np.zeros(4*30*len(L))
Amm = np.zeros(4*30*len(L))

infile = open("./Results/exp_values/expvalsfinal.txt")
input = np.loadtxt(infile.readlines())
"""
        T = thread0[:,0] , thread1[:,0] , thread2[:,0] , thread3[:,0]
        E = thread0[:,3] , thread1[:,3] , thread2[:,3] , thread3[:,3]
        Amm = thread0[:,5] , thread1[:,5] , thread2[:,5] , thread3[:,5]
        Cv = thread0[:,6] , thread1[:,6] , thread2[:,6] , thread3[:,6]
        S = thread0[:,7] , thread1[:,7] , thread2[:,7] , thread3[:,7]
"""

def lattice(L):
    T = input[:,0]
    Cv = input[:,6]
    if(L == 60):
        T = T[0:30*4]
        Cv = Cv[0:30*4]

    if(L == 80):
        T = T[30*4:(30*4)*2]
        Cv = Cv[30*4:(30*4)*2]
    if(L == 100):
        T = T[(30*4)*2:(30*4)*3]
        Cv = Cv[(30*4)*2:(30*4)*3]
    return T, Cv

for i in L:
    results = lattice(i)
    plt.plot(results[0], results[1], label = "L = " + str(i))

plt.legend()
plt.show()
