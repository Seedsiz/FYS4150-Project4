#include "montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

void MonteCarlo::initialize(int L, double T){

}

void MonteCarlo::draw_index(){ // random number generator;
  /* Mersenne twister random generator suggest
    flipping of spin with random index PS: indices needs to be mapped;
    */
  /*int rd = chrono::high_resolution_clock::now().time_since_epoch().count(); //+ rank <--  for parallellization;
  mt19937_64 gen_i(rd);      // seeded with rd
  uniform_int_distribution<> distribution_i(1, (m_L)); // Choose uniform distr. with range 1,(m_L-1) (unsigned integer)
  m_rand_i =  distribution_i(gen_i); // Draw index, flip this  (PS, needs to call for both i and j values)
  int sd = chrono::high_resolution_clock::now().time_since_epoch().count(); //+ rank <--  for parallellization;
  mt19937_64 gen_j(sd);     // seeded with sd
  uniform_int_distribution<> distribution_j(1, (m_L)); // Choose uniform distr. with range 1,(m_L-1) (
  m_rand_j =  distribution_j(gen_j);
  //cout <<  m_rand;
  */
};

void MonteCarlo::draw_acceptance(){
  /* Mersenne twister random generator suggest
  random number between [0,1);
  */
  int sd = chrono::high_resolution_clock::now().time_since_epoch().count(); // Used to obtain seed
  mt19937_64 gen(sd);                                                       // seeded with sd
  uniform_real_distribution<double> distribution(0.0,1.0);                  // creates [0,1)
  m_check =  distribution(gen);                                             // draw acceptance criteria
  //cout << m_check << endl;
};


void MonteCarlo::monte_carlo(vec S){

};


void MonteCarlo::expectation_values(){

}
