"""
Code to read text files and make histogram
with distribution of energies for a steady state
for one particular temperature
PS: need to make textfile for wanted temperature before running
as of now, only handles one temperature
- If you want to plot energies or exp values, change arrays in isingmodel
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

Ncycles = int(input("Enter number of Monte Carlo cycles:"))
L = int(input("Enter number of spins for a given axis:"))
rank = int(input("Enter rank number (assumes only one rank used (1-4)):"))
print("Press 1 to plot energy histogram")
print("Press 2 to plot energy expectation value histogram")
do = int(input("Enter number:"))

# read the energies over cycles (after calibration)
infile = open("./Results/cycles/EMcycles" + str(Ncycles) + \
                "-" + str(L) + "by" + str(L) + "rank" + rank + ".txt", 'r')
infile.readline()

expE_cycles = []   # energies,  could add histogram for M if wanted

for line in infile:
    numbers = line.split()
    expE_cycles.append(float(numbers[4]))
    cycles = int(numbers[1])
infile.close()

expE_cycles = np.array(expE_cycles)

# read off the expectation value and variance
infile = open("./Results/exp_values/expvaluescycles" + str(Ncycles) + \
              "-" + str(L) + "by" + str(L) + ".txt");
infile.readline()

expE = [] # expecation values for different temperatures
stdE = [] # standard deviation calculated with Bessel correction(N-1)

for line in infile:
    numbers = line.split()
    T = float(numbers[0])
    cycles = int(numbers[1])        # cycles after calibration
    expE.append(float(numbers[3]))
    stdE_scalar = np.sqrt(1.0/(float(cycles)-1.0)*float(numbers[8]))  # std with Bessel correction
    stdE.append(stdE_scalar)
infile.close()

T = np.array(T)
expE = np.array(expE)
stdE = np.array(stdE_scalar)
expE_cycles = expE_cycles[Ncycles-cycles:Ncycles] # take out the expectation values that are not included (0 since not filled in)

if do == 1:
    sns.set_style('darkgrid')
    sns.displot(expE_cycles, height=3, aspect=1)

if do == 2:
    # decide upon number of bins
    width_interval = 1e-4 # after 1 cycle, the end point energies can have in each
    range_ =  np.max(expE_cycles) - np.min(expE_cycles) # range of lowest energies and min energies
    N_bins = int(range_/width_interval) + 1 # get number of bins

    sns.set_style('darkgrid')
    sns.displot(expE_cycles,bins = N_bins, height=3, aspect=1) # look into this one: A lot easier to use and super nice


plt.title("{:d}x{:d}, MC cycles:{:.1e}, Calibration Cycles: {:.1e}".format(L,L,Ncycles,Ncycles-cycles))
labels, counts = np.unique(expE_cycles, return_counts=True)
plt.grid(axis='y', alpha=0.75)
if do == 1:
    plt.xlabel('$Energy/L^{2}$')
elif do == 2:
    plt.xlabel('$\mathbb{E}[E]/L^{2}$')

plt.ylabel('Frequency')
plt.yscale('log')
plt.legend(["$\mu=%.3e, \sigma=%.1e$"%(expE, stdE)])
plt.savefig("./Results/Figures/Histograms/hist-{:d}by{:d}spinsystem-T-{:f}-option-{:d}.png".format(L,L,T,do))
plt.show()












#
