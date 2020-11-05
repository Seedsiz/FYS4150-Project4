#include "montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

void MonteCarlo::initialize(int L, double T){
  m_L = L;
  m_T = T;
  }

void MonteCarlo::draw_index(){ // random number generator;
    /* Mersenne twister random generator suggest
    flipping of spin with random index;
    */
  int rd = chrono::high_resolution_clock::now().time_since_epoch().count(); //+ rank <--  for parallellization;
  mt19937_64 gen(rd);      // seeded with rd
  uniform_int_distribution<> distribution(0, (m_L-1)*(m_L-1)); // Choose uniform distr. with range m_L^2 (unsigned integer)
  m_rand =  distribution(gen); // Draw index, flip this
  cout <<  m_rand;
};

void MonteCarlo::draw_acceptance(){
  /* Mersenne twister random generator suggest
  random number between [0,1);
  */
  int sd = chrono::high_resolution_clock::now().time_since_epoch().count(); // Used to obtain seed
  mt19937_64 gen(sd);                                                       // seeded with sd
  uniform_real_distribution<double> distribution(0.0,1.0);                  // creates [0,1)
  m_check =  distribution(gen);                                             // draw acceptance criteria
  cout <<  m_check << endl;
};

void MonteCarlo::metropolis(){
// sampling rule for montecarlo method. Choose if suggested flip should be accepted
draw_acceptance();
double beta = 1/m_T;
double deltaE = 1;
double w = exp(-beta*deltaE);

if (m_check <= w){
// write something here to accept
// Calculate magnetization and energy
  }
}
