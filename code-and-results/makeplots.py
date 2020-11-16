# Programme to read off txt files and plot E,M and accepted flips as
# functions of Monte Carlo cycles
import numpy as np
import matplotlib.pyplot as plt

Ncycles = int(input("Input number of Monte Carlo cycles:"))
print("Press 1 to plot E and M over cycles")
print("Press 2 to plot number of accepted flips over cycles")
do = int(input("Enter number:"))

infile = open("./Results/exp_values/EMcycles" + str(Ncycles) + ".txt", 'r')
infile.readline()

expE_cycles = []   # energies
expM_cycles= []    # magnetic moment
accepted = []      # number of accepted flips

for line in infile:
    numbers = line.split()
    accepted.append(float(numbers[3]))
    expE_cycles.append(float(numbers[4]))
    expM_cycles.append(float(numbers[5]))
    Lspins = int(numbers[2])


expE_cycles = np.array(expE_cycles)
expM_cycles = np.array(expM_cycles)
accepted = np.array(accepted)
cycles = np.linspace(0,Ncycles-1,Ncycles)+1;

L = int(np.sqrt(Lspins))

if do == 1:

    print("Press 1 to make plot for 2x2")
    print("Press 2 to make for other spin system")
    choice = int(input("Enter number:"))

    plt.figure()
    if choice == 1:
        ## analytical values
        T = 1.0;
        L = 2;

        beta = float(1/T);
        z = 12. + 4*np.cosh(8*T);
        exp_E = -(32./z)*np.sinh(8*beta)/L**2;
        mean_abs_M = (8./z)*(np.exp(8*beta) + 2)/L**2;


        ## make arrays
        ana_E = np.zeros(Ncycles) + exp_E
        ana_M = np.zeros(Ncycles) + mean_abs_M

        plt.plot(cycles,ana_E,"-", color = "orange", label = "analytical", alpha = 0.7)

    ### plot energies ###
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,expE_cycles,'--o',label = "numerical <E>",markersize = 1.5)

    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("$\mathbb{E}[E]/L^{2}$",fontsize = 13)
    plt.legend(loc = "upper right",fontsize = 15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/Ecycles-{:d}by{:d}spinsystem.png".format(L,L))
    plt.show()

    ### plot absolute magnetic moment expectation ###
    plt.figure()
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,expM_cycles,'--o',label = "numerical <|M|>",markersize = 1.5)

    if choice == 1:
        plt.plot(cycles,ana_M,"-", label = "analytical", alpha = 0.7)

    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("$\mathbb{E}[|M|]/L^{2}$",fontsize = 13)
    plt.legend(loc = "upper right",fontsize = 15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/Mcycles-{:d}by{:d}spinsystem.png".format(L,L))
    plt.show()

elif do == 2:
    ### plot accepted draws. PS the point 0,0 is not included ###
    plt.figure()
    plt.title("$LxL$:{:d}x{:d}".format(L,L),fontsize = 14)
    plt.plot(cycles,accepted,'--o',markersize = 1.5)

    plt.xlabel("MC cycles",fontsize = 13)
    plt.ylabel("accepted draws",fontsize = 13)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("./Results/Figures/accepted_draws-{:d}by{:d}spinsystem.png".format(L,L))
    plt.show()
