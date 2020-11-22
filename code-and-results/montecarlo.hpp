#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <armadillo>
#include <iostream>
#include <chrono>
#include <random> // To get access to the mersenne twister random generator

using namespace arma;
using namespace std;

class MonteCarlo{

protected:
  int m_calibration; // number of calibration cycles
  int m_nT; // number of temperatures to loop over
  int m_L; // number of grid points (spin particles) // move these two maybe
  int m_L2; // The number of spins in total
  vec m_accepted; // cumulative number of accepted spins as function of cycles
  int m_cumulative_accept;  //total cumulative number of accepted spins over many cycles
  vec m_T; // vector with temperatures to loop over
  int m_MC; //number of monte carlo cycles
  int m_rand_i; // the random index to draw
  int m_rand_j; // the random index to draw
  double m_check; // the random r in [0,1] acceptance criteria
  vec m_map; // a mapping instead of ghost cells
  vec getBoltzmann; // A vector containing boltzmann (deltaEs)
  double m_deltaE; // actual change in energy
  double m_w; // The boltzmann ratio gotten from getBoltzmann
  double m_beta; // m_beta = 1/m_T
  bool m_cont; // continue factor in metropolis algo
  double m_Energy; // total energy of the system
  double m_MagneticMoment; // magnetic moment M
  vec S;     //lattice converted to a flat array
  double exp_val_E, exp_val_E2;  //expectation values for energy and energy squared
  double exp_val_M, exp_val_M2;  //expectation values for magnetization and magnetization squared
  double exp_val_Mabs;   //expectation value for mean absolute value of magnetization
  double m_Cv, m_xi;    //specific heat and susceptibility
  /*random number generator, seeded for each rank once */
  mt19937_64 m_gen;     // seeded with sd
  uniform_real_distribution<double> m_distribution;  // creates [0,1)
  ofstream m_file_emcyc; // cycles to file,  to get access
  ofstream m_file_expv; // expectation values file, to close
  int m_rank;     // The given rank

public:
  void initialize();      // empty initializer
  void draw_acceptance(); // rnd get r [0,1] acceptance criteria,
                          //used in initating the spins system
};

class IsingModel2D: public MonteCarlo{

protected:
  vec S; // A vector containing all spins; must be initalized in a random state

public:
  void init(int L, double T_start, double T_end, int n_T, int MC, int rank); // initiates the spin system
  void setup_boltzmann_ratio(int tempi); // sets up the boltzmann ratio used in metropolis algo
  int magnetic_moment(); // Calculates magnetic moment after initialization
  void metropolis(double w); // accept or not accept drawn spin
  void energy(); // Calculates energy of the system after initialization
  void find_deltaE(int tempi, int flip_i, int flip_j); // Get deltaE and boltzmann ratio of suggested flip
  vec solve(bool save_cycles, int calibration); // Solves the system for MC cycles
  void open_exp_vals_to_file(ofstream&file); // Open exp_vals
  void write_exp_vals_to_file(vec expval,ofstream &file, int temp, double varE, double varM); // write expectation vaules to file
  void open_EM_cycles_to_file(ofstream&file); // Open cycles to file
  void write_EM_cycles_to_file(ofstream&file, vec E, vec M, int temp); // write E and M from each cycle
  void write_spin_to_file(bool check); // write spin to file if check is true
};
#endif
