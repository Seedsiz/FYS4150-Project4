#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <armadillo>
#include <omp.h>
#include <stdio.h>

using namespace std;
using namespace arma;


void menu();

int main(int argc, char const *argv[]){
  menu();
  return 0;
}

void menu(){
  int L; int MC;
  double T_start, T_end;
  int n_T;
  int numthreads;
  bool save_over_cycles = false;
  int calib;

  cout << "Enter integer number of spin particles for each axis:" << " ";
  cin >> L;
  cout << "Enter start point temperature:"  << " ";
  cin >> T_start;
  cout << "Enter an endpoint temperature:"  << " ";
  cin >> T_end;
  cout << "Enter integer number of temperature points to be evaluated:"  << " ";
  cin >> n_T;
  cout << "Enter integer number of MC cycles:"  << " ";
  cin >> MC;
  cout << "Enter integer number of calibration cycles:"  << " ";
  cin >> calib; // eg. 20 000
  cout << "Enter integer number of threads:"  << " ";
  cin >> numthreads;


  //Tryout random generator
  //MonteCarlo mysolver;
  //mysolver.initialize(L,T);
  //mysolver.draw_index();
  //mysolver.draw_acceptance();

  // Start parallelization
  int temps_i;
  vec T_startvec;
  vec T_endvec;
  vec sol;

  // set up T_vec for start and end for each node
  // temperature vector with T_start, T_end for nodes
  double dist = (T_end-T_start)/numthreads;
  T_startvec = linspace<vec>(T_start, T_end-dist, numthreads);
  T_endvec =  linspace<vec>(T_start+dist, T_end, numthreads);
  //cout << T_startvec <<"\n";
  //cout << T_endvec << "\n";

  int rank;
  double start;
  double end;
  start = omp_get_wtime();
  omp_set_num_threads(numthreads);
  #pragma omp parallel for default(shared) num_threads(numthreads) private(temps_i)
  for (temps_i = 0; temps_i < numthreads; temps_i++){
    //T_start = T_startvec(temps_i);
    //cout << T_start << "\n";
    //T_end = T_endvec(temps_i);
    IsingModel2D model; // initate class object;
    model.init(L, T_startvec(temps_i),T_endvec(temps_i), n_T, MC,omp_get_thread_num());
    model.solve(save_over_cycles,calib);
    printf("Thread rank: %d\n", omp_get_thread_num());
    model.close_exp_vals_to_file();
  }
  end = omp_get_wtime();
  printf("Work took %f seconds\n", end - start);
  //IsingModel2D model;
  //model.init(L, T_start,T_end, n_T, MC);
  //model.solve();

  //Catch::Session().run();
}
