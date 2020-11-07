#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

class MonteCarlo{

protected:
  int m_L; // number of grid points (spin particles) // move these two maybe
  int m_T; //temperature
  int m_rand; // the random index to draw
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

public:
  void initialize(int L, double T);
  void draw_index(); // random number generator: get flip index
  void draw_acceptance(); // rnd get r [0,1] acceptance criteria
  void metropolis(double w); // sampling rule;
  void monte_carlo(vec S, int flip_i, int flip_j); // calculates one MC cycle only
  void expectation_values();             // Get expectation_values
};

class IsingModel2D: public MonteCarlo{

protected:
  vec S; // A vector containing all spins; must be initalized in a random state

public:
  void init(int L, double temp);
  void magnetization();
  void energy();
  void find_deltaE(int flip_i, int flip_j);
  void specHeat();
  void solve();
};
#endif
