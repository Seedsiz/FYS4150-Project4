#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

class MonteCarlo {

public:
  void initialize();
  void RNG(); // random number generator;
};

class IsingModel : public MonteCarlo{

protected:
  double energy;

public:
  void init();
  void magnetization();
  void energy();
  void specHeat();
  vec spins;
};
#endif
