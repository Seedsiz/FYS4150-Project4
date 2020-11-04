#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

class MonteCarlo {
protected:
  double energy;

public:
  void initialize();
  void magnetization();
  void energy();
  void specHeat();
  vec spins;
};


#endif
