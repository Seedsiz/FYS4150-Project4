#solving task f) and g) in the prosject description


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


Ncycles = 1000000      #number of monte carlo cycles
Tnum = 20              #number of temperature points for each core
Tpoints = 4*Tnum       #number of temperature points for the region T = [2.2,2.4]

L = np.array((40,60,80,100))    #different values of lattice size L

#storing temperatures, energy, mean magnetization, heat capacity and Susceptibility
#stored for each value of L
T = np.zeros((4,int(Tpoints)))
E = np.zeros((4,int(Tpoints)))
Mabs = np.zeros((4,int(Tpoints)))
Cv = np.zeros((4,int(Tpoints)))
xi = np.zeros((4,int(Tpoints)))


#Reading the files produced by the 4 different cores, rank0,rank1,rank2,rank3
for i in range(len(L)):
    infile_rank0 = open("./Results/exp_values/expvaluescycles" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(0) + ".txt", "r")
    infile_rank0.readline()
    rank0 = np.loadtxt(infile_rank0)
    infile_rank1 = open("./Results/exp_values/expvaluescycles" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(1) + ".txt", "r")
    infile_rank1.readline()
    rank1 = np.loadtxt(infile_rank1)
    infile_rank2 = open("./Results/exp_values/expvaluescycles" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(2) + ".txt", "r")
    infile_rank2.readline()
    rank2 = np.loadtxt(infile_rank2)
    infile_rank3 = open("./Results/exp_values/expvaluescycles" + str(Ncycles) + \
    "-" + str(L[i]) + "by" + str(L[i]) + "rank" + str(3) + ".txt", "r")
    infile_rank3.readline()
    rank3 = np.loadtxt(infile_rank3)
    #rank = np.array((rank0[:,0],rank1[:,0],rank2[:,0],rank3[:,0]))
    for a in range(Tnum):
        T[i,a] = rank0[a,0]
        E[i,a] = rank0[a,3]         #Filling values for exp values from the first core
        Mabs[i,a]= rank0[a,5]
        Cv[i,a] = rank0[a,6]
        xi[i,a] = rank0[a,7]

    a += 1
    for b in range(Tnum):
        T[i,a+b] = rank1[b,0]
        E[i,a+b] = rank1[b,3]        #Filling values for exp values from the second core
        Mabs[i,a+b] = rank1[b,5]
        Cv[i,a+b] = rank1[b,6]
        xi[i,a+b] = rank1[b,7]
    b += 1
    for c in range(Tnum):
        T[i,a+b+c] = rank2[c,0]
        E[i,a+b+c] = rank2[c,3]
        Mabs[i,a+b+c] = rank2[c,3]      #Filling values for exp values from the third core
        Cv[i,a+b+c] = rank2[c,6]
        xi[i,a+b+c] = rank2[c,7]
    c += 1
    for d in range(Tnum):
        T[i,a+b+c+d] = rank3[d,0]
        E[i,a+b+c+d] = rank3[d,3]      #Filling values for exp values from the fourth core
        Mabs[i,a+b+c+d] = rank3[d,5]
        Cv[i,a+b+c+d] = rank3[d,6]
        xi[i,a+b+c+d] = rank3[d,7]



task = str(input("Solve task f or g? (f/g): "))   #solve either task f) or g) ?

if task == "f": #solve task f)
    Lstr = ["L = 40", "L = 60", "L = 80", "L = 100"]

    for i in range(len(L)):    #plotting the various expectation values for different values of L
        plt.plot(T[i,:],xi[i,:], label = Lstr[i])
        plt.title('Phase transitions for suscptibility $\chi$ for L = 40, 60, 80 and 100')
        plt.xlabel('Temperature T')
        plt.ylabel('Susceptibility $\chi$')
        plt.legend(loc = "best",fontsize = 13)
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        plt.tight_layout()
        plt.savefig("./Results/Figures/Calibration-Plots/phasetransitions-xi-L=40-60-80-100.png")

    plt.show()

if task == "g": #solve task g)
    max_x = np.zeros(4)                   #storing the 4 T_C values that corresponds
        for i in range(len(L)):           #to the 4 peaks in the plot for Cv
        Cv_hat = savgol_filter(Cv[i,:], 51, 2)  #making a curvefit to the Cv values
        max_y = max(Cv_hat)                    #maximum Cv value for a given L
        max_x[i] = T[i,:][Cv_hat.argmax()]     #corresponding temp Tc to max_y

    def least_squares(x,y):  #calculating the uncertainty of estimated Tc(infinity)
        n = len(x)
        D = np.sum(x**2) - (1/n)*(np.sum(x)**2)
        E = np.sum(x*y)-(1/n)*np.sum(x)*np.sum(y)
        F = np.sum(y**2)-(1/n)*(np.sum(y)**2)
        delta_alfa = np.sqrt((1/(n-2))*(D*F - E**2)/D**2)
        return delta_alfa

    error = least_squares(1/L,max_x)   #the uncertainty of estimated Tc(infinity)
    print(error)

    m, b = np.polyfit(1/L, max_x, 1)  #making a curvefit to the 4 Tc values for L = 40,60,80,100

    print(b)

    #plotting the curve fitting the 4 Tc values for L = 40,60,80,100

    plt.plot(1/40, max_x[0], 'o', label = '$T_C(40) \\approx 2.2868$')
    plt.plot(1/60, max_x[1], 'o', label = '$T_C(60) \\approx 2.2815$')
    plt.plot(1/80, max_x[2], 'o', label = '$T_C(80) \\approx 2.2789$')
    plt.plot(1/100, max_x[3], 'o', label = '$T_C(100) \\approx 2.2736$')
    plt.plot(1/L, m*(1/L) + b, label = 'Critical $T_{C}(\infty) \\approx 2.267 \pm 0.1565$')
    plt.title('Curvefit for $T_{C}(L) = \\frac{1}{L} + T_{C}(\infty)$ between 4 points')
    plt.xlabel('1/L')
    plt.ylabel('$T_{C}(L)$')
    plt.legend(loc = "best",fontsize = 13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.tight_layout()
    plt.savefig("./Results/Figures/Calibration-Plots/TCcurve-L=40-60-80-100.png")
    plt.show()
