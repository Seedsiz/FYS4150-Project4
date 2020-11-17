# Programme to read off txt files and plot E,M and accepted flips as
# functions of Monte Carlo cycles
# PS: Make sure you ran Monte Carlo with the temperature(s) you want
import numpy as np
import matplotlib.pyplot as plt

Ncycles = int(input("Input number of Monte Carlo cycles:"))
n_temps = int(input("1 or 2 temperatures to be evaluated?"))
if n_temps == 1:
    TT = float(input("Set temperature:"))

print("Press 1 to plot E and M over cycles")
print("Press 2 to plot number of accepted flips over cycles")
do = int(input("Enter number:"))

infile = open("./Results/exp_values/EMcycles" + str(Ncycles) + ".txt", 'r')
infile.readline()

expE_cycles = []   # energies
expM_cycles= []    # magnetic moment
accepted = []      # number of accepted flips
temperatures = []


for line in infile:
    numbers = line.split()
    temperatures.append(float(numbers[0]))
    Lspins = int(numbers[2])
    accepted.append(float(numbers[3]))
    expE_cycles.append(float(numbers[4]))
    expM_cycles.append(float(numbers[5]))

expE_cycles = np.array(expE_cycles)
expM_cycles = np.array(expM_cycles)
accepted = np.array(accepted)

cycles = np.linspace(0,Ncycles-1,Ncycles)+1;

# to handle two temperatures in same figure (case with 20x20 lattice)
Ts = 1.0; Ts2 = 2.4;
if n_temps == 1:
    Ts = TT
    T1 = np.zeros(len(cycles)) + Ts
    T2 = np.zeros(len(cycles)) + Ts2

elif n_temps == 2:
    T1 = np.zeros(2*len(cycles)) + Ts
    T2 = np.zeros(2*len(cycles)) + Ts2


# get temperature indices
temp_points = np.equal(temperatures,T1)
temp_indices = np.where(temp_points == True)
temp_points2p4 = np.equal(temperatures,T2)
temp_indices2p4 = np.where(temp_points2p4 == True)
t1i = temp_indices[0]; t2i = temp_indices2p4[0]

# splitting up arrays
# for temperature 1
expE_cycles1 = np.array(expE_cycles)[t1i]
expM_cycles1 = np.array(expM_cycles)[t1i]
accepted1 = np.array(accepted)[t1i]

# for temperature 2.4
expE_cycles2 = np.array(expE_cycles)[t2i]
expM_cycles2 = np.array(expM_cycles)[t2i]
accepted2 = np.array(accepted)[t2i]

def get_ana(T,ncycles):
    ## analytical values for a 2x2 matrix
    L = 2;
    beta = float(1/T);
    z = 12. + 4*np.cosh(8*beta);
    exp_E = -(32./z)*np.sinh(8*beta)/L**2;
    mean_abs_M = (8./z)*(np.exp(8*beta) + 2)/L**2;
    ## make arrays
    ana_E = np.zeros(ncycles) + exp_E
    ana_M = np.zeros(ncycles) + mean_abs_M
    return ana_E, ana_M

L = int(np.sqrt(Lspins))

if do == 1:
    print("Press 1 to make plot for 2x2")
    print("Press 2 to make for other spin system")
    choice = int(input("Enter number:"))
    plt.figure()
    if choice == 1:
        ana_E1, ana_M1 = get_ana(Ts,len(cycles))
        plt.plot(cycles,ana_E1,"-", color = "red", label = "ana, T:"+str(Ts)+"", alpha = 0.7)
        if n_temps == 2:
            ana_E2, ana_M2 = get_ana(Ts2,len(cycles))
            plt.plot(cycles,ana_E2,"-", color = "purple", label = "ana, T:"+str(Ts2)+"", alpha = 0.7)

    ### plot energies ###
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,expE_cycles1,'--o',label = "num, T:"+str(Ts)+"",markersize = 1.5)
    if n_temps == 2:
        plt.plot(cycles,expE_cycles2,'--o',label = "num, T:"+str(Ts2)+"",markersize = 1.5)
    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("$\mathbb{E}[E]/L^{2}$",fontsize = 13)
    plt.legend(loc = "best",fontsize = 14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/Ecycles-{:d}by{:d}spinsystem-{:d}temps-T1-{:f}.png".format(L,L,n_temps,Ts))
    plt.show()

    ### plot absolute magnetic moment expectation ###
    plt.figure()
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,expM_cycles1,'--o',label = "num, T:"+str(Ts)+"",markersize = 1.5)
    if n_temps == 2:
        plt.plot(cycles,expM_cycles2,'--o',label = "num, T:"+str(Ts2)+"",markersize = 1.5)

    if choice == 1:
        plt.plot(cycles,ana_M1,"-", color = "red", label = "ana, T:"+str(Ts)+"", alpha = 0.7)
        if n_temps == 2:
            plt.plot(cycles,ana_M2,"-", color = "purple", label = "ana, T:"+str(Ts2)+"", alpha = 0.7)

    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("$\mathbb{E}[|M|]/L^{2}$",fontsize = 13)
    plt.legend(loc = "best",fontsize = 14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/Mcycles-{:d}by{:d}spinsystem-{:d}temps-T1-{:f}.png".format(L,L,n_temps,Ts))
    plt.show()

elif do == 2:
    ### plot accepted draws. PS the point 0,0 is not included ###
    plt.figure()
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,accepted1,'--o',label = "T:"+str(Ts)+"",markersize = 1.5)
    if n_temps == 2:
        plt.plot(cycles,accepted2,'--o',label = "T:"+str(Ts2)+"",markersize = 1.5)
    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("accepted draws",fontsize = 13)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc = "best",fontsize = 14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/accepted_draws-{:d}by{:d}spinsystem-{:d}temps-T1-{:f}.png".format(L,L,n_temps,Ts))
    plt.show()

#print out which cycle which cycle that gets close to convergence
if choice == 1:# for 2 by 2  for ts (first temperature = 1 if otherwise stated)
    tol = 1e-3
    epsE = expE_cycles1 - ana_E1
    Eindex = np.where(np.abs(epsE) < tol)[0][0] # get first elem
    epsM = expM_cycles1 - ana_M1;
    Mindex = np.where(np.abs(epsM) < tol)[0][0]
    print("MC cycle close to analytical - Eindex: {:d}. Mindex: {:d}".format(Eindex,Mindex))

#print out which cycle which cycle that gets close to convergence
if choice == 2:
    # eg. 20 by 20 grid
    tol = 1e-5 # need smaller tolerance because comparing between neighbouring points
    epsE = np.zeros(len(cycles)-1)
    epsM = np.zeros(len(cycles)-1)
    for i in range(len(cycles)-1):
        epsE[i] = expE_cycles1[i+1] - expE_cycles1[i]
        epsM[i] = expM_cycles1[i+1] - expM_cycles1[i]

    # find indexes of convergence
    Eindex = np.where(np.abs(epsE) < tol)[0][0]
    Mindex = np.where(np.abs(epsM) < tol)[0][0]

    print("MC cycle close to analytical - Eindex: {:d}. Mindex: {:d}".format(Eindex+1,Mindex+1))
