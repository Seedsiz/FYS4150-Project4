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
  vec S;     //lattice converted to a flat array
  vec getBoltzmann;
  double m_w;

public:
  void initialize(int L, double T);
  void draw_index(); // random number generator: get flip index
  void draw_acceptance(); // rnd get r [0,1] acceptance criteria
  void metropolis(); // sampling rule;
};

class IsingModel2D: public MonteCarlo{

protected:
  double E; // total energy of the system
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
