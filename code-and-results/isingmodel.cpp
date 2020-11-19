#include "montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace arma;
using namespace std;


void IsingModel2D::init(int L, double T_start, double T_end, int n_T, int MC, int rank){
  // Create mapping vector so that physical mesh points are
  // connected to ghost cells (check that these are int!!!)
  //initialize(L,temp);
  m_rank = rank;
  m_nT = n_T; // number of temperatures to loop over
  m_T = linspace<vec>(T_start, T_end, n_T);  // temperature vec to loop over
  m_L = L; // number of spins along a given axis (square 2D system)
  m_MC = MC; // number of Monte carlo cycles
  m_map = vec(m_L+2);
  m_map(0) = m_L-1; // first index in map --> last index in physical mesh
  m_map(m_L+1) = 0; // last index in map --> first index in physical mesh


  for (int i = 0; i < m_L; i++){
    m_map(i+1) = i; // mapping is correct
  }

  // Setting up boltzmann ratio vector
  getBoltzmann = zeros<vec>(17);

  S = vec(L*L);   //Setting up lattice of L*L elements
  draw_acceptance();    //Getting random number
  if(m_T(0) >= 3) {        //Temperature check
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
  //cout << S << "\n";
}

void IsingModel2D::setup_boltzmann_ratio(int temp_i){ // function to call for each temp
  // 5 different energy states: Scale J to 1;
    m_beta = 1/((double) m_T(temp_i));
    getBoltzmann(0) = exp(m_beta*8.0); // boltzmann for deltaE = -8
    getBoltzmann(4) = exp(m_beta*4.0); // boltzmann for deltaE = -4
    getBoltzmann(8) = 1.0; //deltaE = 0 exp^0 = 1 // boltzmann for deltaE = 0
    getBoltzmann(12) = exp(-m_beta*4.0); // boltzmann for deltaE = 4
    getBoltzmann(16) = exp(-m_beta*8.0); // boltzmann for deltaE = 8
}

int IsingModel2D::magnetic_moment(){
  /* Code for magnetization for one specific
  state with periodic boundary conditions (2D)
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
for (int i = 1; i <= m_L; i++){
  for (int j = 1; j <= m_L; j++){
    i_p = m_map(i); j_p = m_map(j); // mapping to physical mesh points
    i_pp = m_map(i+1); j_pp = m_map(j+1); // mapping to physical mesh points
    m_Energy -= S(i_p*m_L+j_p)*S(i_p*m_L+j_pp) + \
        S(i_p*m_L + j_p)*S((i_pp)*m_L + j_p);
    }
  }
}

void IsingModel2D::find_deltaE(int temp_i,int i, int j){
  // take in suggested random indices
  // use the sum of sorrounding spins
  int S_candid =  S(m_map(i)*m_L + m_map(j)); // spin suggested flipped
  int S1 =  S(m_map(i-1)*m_L + m_map(j));
  int S2 =  S(m_map(i+1)*m_L + m_map(j));
  int S3 =  S(m_map(i)*m_L + m_map(j-1));
  int S4 =  S(m_map(i)*m_L + m_map(j+1));
  m_deltaE = 2*S_candid*(S1 + S2 + S3 + S4);
  int mapping = m_deltaE + 8;
  m_w = getBoltzmann(mapping); //
  //cout << "Sc   " << S_candid << "\n";
  //cout << "sum S  " << (S1 + S2 + S3 + S4) << "\n";
  //cout << m_deltaE;
  //cout <<"map:"<< mapping << "\n";
}


void IsingModel2D::metropolis(double w){
// sampling rule for montecarlo method. Choose if suggested flip should be accepted
  if (m_deltaE < 0.0){
    m_Energy += m_deltaE; // Calculating value of cycle
    S(m_map(m_rand_i)*m_L + m_map(m_rand_j)) *= -1.0;    // if true, flip one spin and accept new spin config
    m_MagneticMoment += 2*S(m_map(m_rand_i)*m_L + m_map(m_rand_j)); // check why this is like this
    m_cumulative_accept += 1;
  }

  else if (m_distribution(m_gen) <= w){
    m_Energy += m_deltaE; // Calculating value of cycle
    S(m_map(m_rand_i)*m_L + m_map(m_rand_j)) *= -1.0;    // if true, flip one spin and accept new spin config
    m_MagneticMoment += 2*S(m_map(m_rand_i)*m_L + m_map(m_rand_j)); // check why this is like this
    m_cumulative_accept += 1;
    //m_check =  distribution(gen);
    //cout << "dE" << m_deltaE << "\n";
  }
}


vec IsingModel2D::solve(bool save_cycles, int calibration){ // calibration: number of calibration cycles
  // calculates all cycles
  // sends in the indices suggested if metropolis gives true
  // update expectation values and flip
  m_calibration = calibration;
  m_L2 = m_L*m_L;
  vec E_cycles = vec(m_MC); //Store last energy for each cycle
  vec M_cycles = vec(m_MC); // Store last magnetic moment each cycle
  vec exp_values = vec(5); // The final expectation values
  m_accepted = vec(m_MC);

  /* Mersenne twister random generator suggest
  flipping of spin with random index. PS: indices are thereafter mapped;
  */ // seed once in the beginning
  int rd = chrono::high_resolution_clock::now().time_since_epoch().count()+ m_rank; // <--  for parallellization;
  mt19937_64 gen_i(rd);      // seeded with rd
  uniform_int_distribution<> distribution_i(1, (m_L)); // Choose uniform distr. with range 1,(m_L) (unsigned integer)

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //  <--  for parallellization;
  mt19937_64 gen_j(sd);     // seeded with sd
  uniform_int_distribution<> distribution_j(1, (m_L)); // Choose uniform distr. with range 1,(m_L)

  int td = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //<--  for parallellization;
  m_gen.seed(td);

  open_exp_vals_to_file(m_file_expv); // opens file to be written

  if (save_cycles == true){
    open_EM_cycles_to_file(m_file_emcyc);
  }

  for (int temp = 0; temp < m_nT; temp++){
    // Initialize to zero for each temperature
    // state of S at in last MC cycle for previous temp
    m_accepted.fill(0.0);
    exp_val_E = 0.0; exp_val_E2 = 0.0;
    exp_val_M = 0.0; exp_val_M2 = 0.0;
    exp_val_Mabs = 0.0;
    m_cumulative_accept = 0.0;

    energy(); // calculate initial energy
    magnetic_moment(); // calculate initial magnetic moment

    setup_boltzmann_ratio(temp); // get right beta = 1/T

    // Calibration cycles: run through some percentage of samples before
    // adding energies and magnetic moment to expectation values and variances
    for (int c = 0; c < m_calibration; c++){
      for (int i = 0; i < m_L*m_L; i++){
        m_rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        m_rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        find_deltaE(temp, m_rand_i, m_rand_j);  //calculating deltaE and m_w for flip of the random indices
        //m_check =  distribution(gen);
        metropolis(m_w);       // draw acceptance criteria, flip or not
      }
      // store accepted flips:
      //m_accepted(c) = m_cumulative_accept;
    }

    // cycles contributing to mean and variance
    for (int c = m_calibration; c < m_MC; c++){
      for (int i = 0; i < m_L*m_L; i++){
        m_rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        m_rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        find_deltaE(temp, m_rand_i, m_rand_j);  //calculating deltaE and m_w for flip of the random indices
        //m_check =  distribution(gen);
        metropolis(m_w);       // draw acceptance criteria, flip or not
      }
      //adding expectation values from each cycle
      exp_val_E += m_Energy;
      exp_val_E2 += m_Energy*m_Energy;
      exp_val_M += m_MagneticMoment;
      exp_val_M2 += m_MagneticMoment*m_MagneticMoment;
      exp_val_Mabs += fabs(m_MagneticMoment);
      // store accepted flips:
      //m_accepted(c) = m_cumulative_accept;

      // first store endpoint energy value for this cycle
      // to get histograms (could also do this for plotcycles)
      //E_cycles(c) = m_Energy;
      //M_cycles(c) = m_MagneticMoment;

      // or one could store expectation values (for plotcycles)
      //E_cycles(c) = exp_val_E/((double) c - m_calibration+1);
      //M_cycles(c) = exp_val_Mabs/((double) c - m_calibration+1);
    }
    //Get final expectation value over all cycles for this temperature: Dividing the sum with number
    //of MC cycles m_MC to get expectation values.

    int c_contrib = m_MC - m_calibration; // cycles contributing to expectation values

    exp_values(0) = exp_val_E/((double) c_contrib);
    exp_values(1) = exp_val_M/((double) c_contrib);
    exp_values(2) = exp_val_Mabs/((double) c_contrib);

    exp_val_E2 = exp_val_E2/((double) c_contrib);
    exp_val_M2 = exp_val_M2/((double) c_contrib);

    //Calculating variance for energy and magnetization
    //Finding specific heat m_Cv and suceptibility m_xi
    double varianceE = exp_val_E2 - exp_values(0)*exp_values(0);
    double varianceM = exp_val_M2 - exp_values(1)*exp_values(1);
    m_Cv = varianceE/(m_T(temp)*m_T(temp));
    m_xi = varianceM/m_T(temp);
    exp_values(3) = m_Cv;
    exp_values(4) = m_xi;

    //cout << setprecision(15) << ((double) 1/m_L2)*exp_values << "\n";
    //cout << E_cycles; //
    //cout << M_cycles;

    // Scaling by L^2 spins because m_L is an intrinsic parameter
    E_cycles = ((double) 1/m_L2)*E_cycles;
    M_cycles = ((double) 1/m_L2)*M_cycles;
    exp_values = ((double) 1/m_L2)*exp_values;
    varianceE = ((double) 1/m_L2)*varianceE;
    varianceM = ((double) 1/m_L2)*varianceM;

    write_exp_vals_to_file(exp_values,m_file_expv,temp, varianceE,varianceM);
    if (save_cycles == true){
      write_EM_cycles_to_file(m_file_emcyc, E_cycles, M_cycles, temp);
    }
  }
  m_file_expv.close();
  if (save_cycles == true){
    m_file_emcyc.close();
  }
  return exp_values; // return something/or just write to file above?
}

void IsingModel2D::open_EM_cycles_to_file(ofstream&file){
  string filename("./Results/cycles/EMcycles" + to_string(m_MC) + \
                "-" + to_string(m_L) + "by" + to_string(m_L) + ".txt");
  file.open(filename);
  file  << "T" << setw(20) << "MC_cycles-ac" \
        << setw(20) << "N_spins" << setw(20) << "accepted flips" << setw(20) \
        << "<E>/N" << setw(20) << "<|M|>/N";
  file << "\n";
}

void IsingModel2D::write_EM_cycles_to_file(ofstream&file, vec E, vec M, int temp){ // write E,M for each cycle
  for (int cycle = 0; cycle < m_MC; cycle++){
      file    << setprecision(15) <<  m_T(temp) << setw(20) << m_MC-m_calibration << setw(20) \
              << m_L2 << setw(20) << m_accepted(cycle) << setw(20) << E(cycle) \
              << setw(20) << M(cycle);
      file << "\n";
      }
}


void IsingModel2D::open_exp_vals_to_file(ofstream&file){ // write expectation values after all cycles
  string filename("./Results/exp_values/expvaluescycles" + to_string(m_MC) + \
                  "-" + to_string(m_L) + "by" + to_string(m_L) + "rank" + to_string(m_rank) + ".txt");
  file.open(filename);
  file    << "T" << setw(25) << "MC_cycles-ac" << setw(25) << "N_spins" << setw(25)\
          << "<E>/N" << setw(25) << "<M>/N" << setw(25) <<  "<|M|>/N" << setw(25) \
          << "Cv" << setw(25) << "Xi" << setw(25) << "varE" << setw(25) << "varM";
  file << "\n";
}

void IsingModel2D::write_exp_vals_to_file(vec expval,ofstream&file, int temp, double varE, double varM){
  // write energies, magnetization, number of MC cycles,  to file
  // write in the end expectation values to file
  // post-process these in python.
  file    << setprecision(15) << m_T(temp) << setw(25) << m_MC-m_calibration << setw(25) << m_L2 << setw(25)\
          << expval(0) << setw(25) << expval(1) << setw(25) <<  expval(2) << setw(25) \
          << expval(3) << setw(25) << expval(4) << setw(25) << varE << setw(25) << varM;
  file << "\n";
}
