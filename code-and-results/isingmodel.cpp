#include "montecarlo.hpp"

#include <armadillo>
#include <iostream>
#include <iomanip>

using namespace arma;
using namespace std;


void IsingModel2D::init(int L, double T, int MC){
  // Create mapping vector so that physical mesh points are
  // connected to ghost cells (check that these are int!!!)
  //initialize(L,temp);
  m_T = T; // temperature
  m_beta = 1/m_T;
  m_L = L; // number of spins along a given axis (square 2D system)
  m_MC = MC; // number of Monte carlo cycles
  m_map = vec(m_L+2);
  m_map(0) = m_L-1; // first index in map --> last index in physical mesh
  m_map(m_L+1) = 0; // last index in map --> first index in physical mesh


  for (int i = 0; i < m_L; i++){
    m_map(i+1) = i; // mapping is correct
  }

  m_beta = 1/m_T;
  getBoltzmann = vec(9); // 5 different energy states: Scale J to 1;
  getdeltaE = vec(9);
  getBoltzmann(0) = exp(m_beta*8.0);
  getBoltzmann(2) = exp(m_beta*4.0);
  getBoltzmann(4) = 1.0;       //exp(0)
  getBoltzmann(6) = exp(-m_beta*4.0);
  getBoltzmann(8) = exp(-m_beta*8.0);
  getdeltaE(0) = -8.0;
  getdeltaE(2) = -4.0;
  getdeltaE(4) = 0.0;
  getdeltaE(6) = 4.0;
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
      if(m_check < 0.5) {             //the lattice is filled with either only
        S(i) = -1;                    //positive spins, or only negative.
      }else {
        S(i) = 1;
      }
    }
  }
  cout << S << "\n";
}


int IsingModel2D::magnetic_moment(){
  /* Code for magnetization for one specific
  state with periodic boundary conditions (2D

  Calculating total magnetization, by summing over all spins
  for one specific state */
  m_MagneticMoment = 0;
  for (int i = 0; i < m_L*m_L; i++){
    m_MagneticMoment += S(i);
  }
  return m_MagneticMoment;
}

void IsingModel2D::energy(){
/* Code for energy for one specific
state with periodic boundary conditions (2D) */
int i_p; int j_p;
int i_pp; int j_pp;   //Same as i_p and j_p only i+1 and j+1
// Calculating total energy by multiplying below and to the right
m_Energy = 0.0;
for (int i = 1; i < m_L+1; i++){
  for (int j = 1; j < m_L+1; j++){
    i_p = m_map(i); j_p = m_map(j); // mapping to physical mesh points
    i_pp = m_map(i+1); j_pp = m_map(j+1); // mapping to physical mesh points
    m_Energy += S(i_p*m_L+j_p)*S(i_p*m_L+j_pp) + \
        S(i_p*m_L + j_p)*S((i_pp)*m_L + j_p);
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
  //cout <<"map:"<< mapping << "\n";
  //cout << "dE: "<< m_deltaE << "\n";
}


void IsingModel2D::expectation_values(){

}

void IsingModel2D::specHeat() {
}


vec IsingModel2D::solve(){
  // calculates one cycle only
  // sends in the indices suggested if metropolis gives true
  // update expectation values and flip
  int m_L2 = m_L*m_L;
  vec exp_E_cycles = vec(m_MC); //Store mean energy for each cycle
  vec exp_M_cycles = vec(m_MC); // Store mean magnetic moment each cycle
  vec exp_values = vec(5); // The final expectation values
  m_accepted = vec(m_MC);

  // Initialize to zero
  m_accepted.fill(0);
  exp_val_E = 0; exp_val_E2 = 0;
  exp_val_M = 0; exp_val_M2 = 9;
  exp_val_Mabs = 0;
  int cumulative_accept = 0;

  energy(); // calculate initial energy
  magnetic_moment(); // calculate initial magnetic moment
  for (int c = 0; c < m_MC; c++){
    for (int i = 0; i < m_L*m_L; i++){
      draw_index();                     //drawing a random index i and j from the lattice S
      find_deltaE(m_rand_i, m_rand_j);  //calculating deltaE and m_w for flip of the random indices
      metropolis(m_w);                 //returns true or false, given deltaE
      if (m_cont == true){  // !!! Uncertain if this if statement needs to be here (could have only one)
    	  m_Energy += S(m_map(m_rand_i)*m_L + m_map(m_rand_j))*m_deltaE; // Calculating mean expectation value of cycle, get right sign
        cumulative_accept += 1;
        S(m_map(m_rand_i)*m_L + m_map(m_rand_j)) *= -1.0;    // if true, flip one spin and accept new spin config
    	  m_MagneticMoment += 2*S(m_map(m_rand_i)*m_L + m_map(m_rand_j)); // check why this is like this
        //cout << m_deltaE << "\n";
        //cout << m_Energy << "\n"; <--- Something is wrong here: only increases negatively

        }
      }
    // first store expectation value for this cycle
    exp_E_cycles(c) = m_Energy;
    exp_M_cycles(c) = m_MagneticMoment;

    //adding expectation values from each cycle
    exp_val_E += m_Energy;
    exp_val_E2 += m_Energy*m_Energy;
    exp_val_M += m_MagneticMoment;
    exp_val_M2 += m_MagneticMoment*m_MagneticMoment;
    exp_val_Mabs += fabs(m_MagneticMoment);

    // store accepted flips:
    m_accepted(c) = cumulative_accept;
  }
  //Get final expectation value over all cycles: Dividing the sum with number
  //of MC cycles m_MC to get expectation values.
  exp_values(0) = exp_val_E/((double) m_MC);
  exp_values(1) = exp_val_M/((double) m_MC);
  exp_values(2) = exp_val_Mabs/((double) m_MC);

  exp_val_E2 = exp_val_E2/((double) m_MC);
  exp_val_M2 = exp_val_M2/((double) m_MC);

  //Calculating variance for energy and magnetization
  //Finding specific heat m_Cv and suceptibility m_xi
  double varianceE = exp_val_E2 - exp_values(0)*exp_values(0);
  double varianceM = exp_val_M2 - exp_values(1)*exp_values(1);
  m_Cv = varianceE/(m_T*m_T);
  m_xi = varianceM/m_T;
  exp_values(3) = m_Cv;
  exp_values(4) = m_xi;

  //cout << setprecision(15) << ((double) 1/m_L2)*exp_values << "\n";
  //cout << m_accepted;
  cout << exp_E_cycles; //
  //cout << exp_M_cycles;

  // Scaling by L^2 spins because m_L is an intrinsic parameter
  exp_E_cycles = ((double) 1/m_L2)*exp_E_cycles;
  exp_M_cycles = ((double) 1/m_L2)*exp_M_cycles;
  return ((double) 1/m_L2)*exp_values;
}
