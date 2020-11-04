#include "montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
using namespace std;

void MonteCarlo::initialize(int L, double T){
  m_L = L;
  m_T = T;
  }

void MonteCarlo::draw_index(){ // random number generator;
    /* Mersenne twister random generator suggest
    flipping of spin with random index;
    */
  /* Where to place initializer to only seed once? use typedefs?*/
  // Initialize random number generator with a seed (only initialize one time)
  random_device rd;          // Used to obtain seed
  mt19937_64 gen(rd());      // seeded with rd
  uniform_int_distribution<> distribution(0, (m_L-1)*(m_L-1)); // Choose uniform distr. with range m_L^2 (unsigned integer)
  //  Could use something like
  // std::uniform_int_distribution<unsigned long long> dis (above)
  m_rand =  distribution(gen); // Draw index, flip this
  cout <<  m_rand << endl;
};

void MonteCarlo::draw_acceptance(){
  random_device sd;          // Used to obtain seed
  mt19937_64 gen(sd());      // seeded with rd
  uniform_real_distribution<double> distribution(0.0,1.0); // creates [0,1)
  m_check =  distribution(gen); // draw acceptance criteria
  cout <<  m_check << endl;
}

void MonteCarlo::metropolis(){
// sampling rule for montecarlo method. Choose if suggested flip should be accepted
// double w = exp(-beta*deltaE
//double Aji = min(1,exp(-beta*deltaE))
// if (m_rand =< Aji){
// Accept move move = True;
//}

}
