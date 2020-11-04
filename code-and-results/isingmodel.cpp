#include "montecarlo.hpp"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

void init(){

}

void IsingModel2D::magnetization(){

}

void IsingModel2D::energy(){
/* Code for energy for one specific
state with periodic boundary conditions (2D) */

// Calculating total energy by multiplying below and to the right
int k = m_L + 1; // So that index is L + 1;
for (int i = 0; i < m_L; i++){
  for (int j = 0; j < m_L; j++){
    E += S(i*k+j)*S(i*k+j+1) + S(i*k + j)*S((i+1)*k + j);
    }
  }
}

void IsingModel2D::specHeat() {
}

void IsingModel2D::solve(){

}
