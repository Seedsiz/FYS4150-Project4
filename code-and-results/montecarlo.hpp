#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

class MonteCarlo{

protected:
  int m_L; // number of grid points (spin particles) // move these two maybe
  vec m_accepted; // cumulative number of accepted spins
  double m_T; //temperature
  int m_MC; //number of monte carlo cycles
  int m_rand_i; // the random index to draw
  int m_rand_j; // the random index to draw
  double m_check; // the random r in [0,1] acceptance criteria
  vec m_map; // a mapping instead of ghost cells
  vec getBoltzmann; // A vector containing boltzmann (deltaEs)
  vec getdeltaE;    // Used to add delta E (get m_deltaE) to Energy in monte carlo cycle
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

public:
  void initialize(int L, double T);
  void draw_index(); // random number generator: get flip index
  void draw_acceptance(); // rnd get r [0,1] acceptance criteria
  void metropolis(double w); // sampling rule;
  void monte_carlo(vec S); // calculates one MC cycle only
  void expectation_values();             // Get expectation_values
};

class IsingModel2D: public MonteCarlo{

protected:
  vec S; // A vector containing all spins; must be initalized in a random state

public:
  void init(int L, double T, int MC);
  int magnetic_moment();
  void expectation_values();
  void energy();
  void find_deltaE(int flip_i, int flip_j);
  void specHeat();
  vec solve();
};
#endif
