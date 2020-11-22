#include "montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

void MonteCarlo::initialize(){
}

void MonteCarlo::draw_acceptance(){
  /* Mersenne twister random generator suggest
  random number between [0,1);
  */
  int sd = chrono::high_resolution_clock::now().time_since_epoch().count(); // Used to obtain seed
  mt19937_64 gen(sd);                                                       // seeded with sd
  uniform_real_distribution<double> distribution(0.0,1.0);                  // creates [0,1)
  m_check =  distribution(gen);                                             // draw acceptance criteria
};
