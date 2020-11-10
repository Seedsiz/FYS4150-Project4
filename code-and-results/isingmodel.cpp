#include "montecarlo.hpp"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;


void IsingModel2D::init(int L, double temp){
  // Create mapping vector so that physical mesh points are
  // connected to ghost cells (check that these are int!!!)
  m_map = vec(m_L+2);
  m_map(0) = m_L-1;
  m_map(m_L+1) = 0;
  m_T = temp;
  for (int i = 0; i <= m_L; i++){
    m_map(i+1) = i;
  }
  // remember to cout to see if this is correct
  m_beta = 1/m_T;
  getBoltzmann = vec(9); // 5 different energy states: Scale J to 1;
  getdeltaE = vec(9);
  getBoltzmann(0) = exp(-8.0);
  getBoltzmann(2) = exp(-4.0);
  getBoltzmann(4) = 1.0;       //exp(0)
  getBoltzmann(8) = exp(8.0);
  getdeltaE(0) = -8.0;
  getdeltaE(2) = -4.0;
  getdeltaE(4) = 0.0;
  getdeltaE(8) = 8.0;

  S = vec(L*L);   //Setting up lattice of L*L elements
  draw_acceptance();    //Getting random number
  if(m_T >= 1) {        //Temperature check
    for(int i = 0; i < L*L; i++) {    //If the temperature is greater than 1,
      if(m_check < 0.5) {             //the lattice is filled with random spins.
        S(i) = -1;
      }else {
        S(i) = 1;
      }
      draw_acceptance();
    }
  }else {
    for(int i = 0; i < L*L; i++) {    //If the temperature is smaller than 1,
      if(m_check < 0.5) {             //the mattice is filled with either only
        S(i) = -1;                    //positive spins, or only negative.
      }else {
        S(i) = 1;
      }
    }
  }
}

void IsingModel2D::magnetization(){
  /* Code for magnetization for one specific
  state with periodic boundary conditions (2D

  Calculating total magnetization, by summing over all spins
  for one specific state */
  m_magnetization = 0;
  for (int i = 0; i < m_L*m_L; i++){
    m_magnetization += S(i);
  }
}

void IsingModel2D::energy(){
/* Code for energy for one specific
state with periodic boundary conditions (2D) */
int i_p; int j_p;
// Calculating total energy by multiplying below and to the right
for (int i = 1; i <= m_L+1; i++){
  for (int j = 1; j <=m_L+1; j++){
    i_p = m_map(i); j_p = m_map(j); // mapping to physical mesh points
    m_Energy += S(i_p*m_L+j_p)*S(i_p*m_L+j_p+1) + \
        S(i_p*m_L + j_p)*S((i_p+1)*m_L + j_p);
    }
  }
}

void IsingModel2D::find_deltaE(int i, int j){
  // take in suggested random indices
  // use the sum of spins, then map to a deltaE;
  int S1 =  S(m_map(i-1)*m_L + m_map(j));
  int S2 =  S(m_map(i+1)*m_L + m_map(j));
  int S3 =  S(m_map(i)*m_L + m_map(j-1));
  int S4 =  S(m_map(i)*m_L + m_map(j+1));
  int spin_sum = S1 + S2 + S3 + S4;
  int mapping = spin_sum + 4;
  m_w = getBoltzmann(mapping); //
  m_deltaE = getdeltaE(mapping);
}


void MonteCarlo::expectation_values(){

  m_MagneticMoment = 0;
  for (int i = 0; i < L*L; i++){
    magnetization();



}

void IsingModel2D::specHeat() {
}


void IsingModel2D::solve(){
}
