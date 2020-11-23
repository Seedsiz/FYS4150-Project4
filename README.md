# FYS4150-Project4
Git for Project 4 in Computational Physics (FYS4150).

### Main overview
* The programs in this repository aim at finding expectation values for a 2D Ising Model and studying phase transition in a magnetic system, as described in the : [Project 4 - Studies of phase transitions in magnetic systems ](https://github.com/Seedsiz/FYS4150-Project4/blob/main/report/Project%20description.pdf). The final report can be found at: [title](link).

* The main challenge was to apply Monte Carlo simulations and the metropolis sampling algorithm to obtain expectation values and study phase transition for a 2D system. To do so a random number generator was used (see link below). In addition the algorithm was parallelized (over temperatures) using OpenMP. Each thread was regarded as a separate experiment, thus each thread has its own unique seed. The spin matrix is initialized once for each thread. For temperatures within that thread, the first Monte Carlo (MC) cycle of the next temperature uses the last spin configuration of the previous temperature.

* Another central task was to verify the algorithm for a 2 by 2 matrix (numerical vs. analytical results), and check method performance with/without parallelization.

* Textfiles and figures can be found in the folder Results.

* The calculations are performed with the energy coefficent J, and boltzmann factor, k, set to 1. With this scaling the critical temperature is approximately 2.269 in the termodynamical limit L -> infinity (infinite number of spins in 2D system).

### Code: Link and description of programmes
- [main.cpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/main.cpp) : Runs the other programmes and provide user options through terminal.

 - [makefile](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/makefile) : Compiles and executes cpp files, and provides plot options of histograms, spin states for variables, expectation values and physical attributes (heath capacity, magnetic susceptibility etc)  

-  [montecarlo.hpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/montecarlo.hpp) : Headerfile for the superclass MonteCarlo, with subclasses IsingModel2D.

- [montecarlo.cpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/montecarlo.cpp) : Provides the superclass method draw_acceptance, which draws a random number between [0,1). Is used to initialize the spin matrix S.
- [isingmodel.cpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/isingmodel.cpp) : The subclass which attains the expectation values for a 2D spin matrix at equilibrium for a given temperature T. Subclass methods provided are given in following order
  1. init: Sets up member parameters and vectors, and initializes the flattened spin matrix.
  2. setup_boltzmann_ratio: Sets up the boltzmann ratio which is used as acceptance criteria for a specific temperature
  3. magnetic_moment: Get magnetic moment of the initial spin state.
  4. energy: Get energy of the initial spin state.
  5. find_deltaE: Find the energy difference between suggested state and present spin configuration. Also gives the boltzmann ratio used in montecarlo sampling.
  6. metropolis: Accepts all flips with dE < 0. For other dE compare boltzmann ratio with random number [0,1). If boltzmann ratio larger, accept.
  7. solve: Gives the energies and magnetic moment, and calculates the expectation values of the 2D Isingmodel for a given number of calibration cycles and Monte Carlo simulations. Writes to file.
  8. The other methods provided are write to file methods.

- [test.cpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/test.cpp): Unit tests for a 2 by 2 system with temperature 1. Has headerfile [catch.hpp](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/catch.hpp)

- [plotcycles.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plotcycles.py) plots the expectation values or energy and magnetic momentum values as a function of MC cycles.
- [plothist.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plothist.py) visualizes the energy distribution or energy expectation distribution for a given temperature as histograms.
- [plotspin.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plotspin.py) makes a heatmap of the final spin configuration for a given temperature.
- [plotexpvalues.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plotexpvalues.py) plots the heat capacity, magnetic susceptibility and other parameters of the system.
- [plottime.py](https://github.com/Seedsiz/FYS4150-Project4/blob/main/code-and-results/plottime.py) makes a plot for different runtime with and without flags/parallelization.

The files can be compiled with "make all". And plots can be generated with "make argument."

How to run the programmes to reproduce the results discussed in the article: The menu provides user input on the following variables
  1. Integer number of spin particles for an axis, eg. 2,20,40,60,80,100
  2. Start point temperature, eg. 1
  3. Endpoint temperature, eg. 2.4
  4. Integer number of temperature points to be evaluated within a thread, eg. 10
  5. integer number of MC cycles, eg. 100 thousands or 1-10 million to get a good result
  6. Enter integer number of calibration cycles, 20 000 to be on the safe side.
  7. Enter integer number of threads (1 if not parallelization wanted), a value between 1-4.

After experimentation with the MC cycles, the authors decided to use 20 000 calibration cycles, to get more accurate expectation values for less number of cycles. This ensures also good results for larger spin systems and higher temperatures, which needs more calibration cycles compared to lower temperatures and smaller spin system. If one want to plot more than just expectation values, one must go in the isingmodel.cpp to decomment variables m_accepted, E_cycles and M_cycles. Then you must set the bool objects to true.

### Links and packages
- The Mersenne Twister (pseudo)random number generator was used in generating uniform distribution to draw indices and acceptance criteria. Documentation on the class mt19937_64 can be found [here.](https://www.cplusplus.com/reference/random/mt19937_64/)

- Documentation for Matplotlib from python from [here](https://matplotlib.org/)

- Documentation on seaborn displot [here](https://seaborn.pydata.org/generated/seaborn.displot.html#seaborn.displot) and heatmaps [here](https://seaborn.pydata.org/generated/seaborn.heatmap.html)

- Documentation on pandas dataframes can be found [here](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)

- Documentation on parallelization with OpenMP can be found [here](https://www.openmp.org/wp-content/uploads/OpenMP-4.5-1115-CPP-web.pdf) or for more versions [here](https://www.openmp.org/resources/refguides/)

- Documentation on armadillo can be found [here](http://arma.sourceforge.net/docs.html)
