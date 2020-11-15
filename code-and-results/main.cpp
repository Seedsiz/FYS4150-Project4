#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <armadillo>
#include <omp.h>

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
  int num_threads;
  bool save_over_cycles = false;

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
  cout << "Enter integer number of threads:"  << " ";
  cin >> num_threads;

  //Tryout random generator
  //MonteCarlo mysolver;
  //mysolver.initialize(L,T);
  //mysolver.draw_index();
  //mysolver.draw_acceptance();

  // Start parallelization
  int temps_i;
  vec T_vec;
  vec sol;

  // set up T_vec for start and end for each node
  // temperature vector with T_start, T_end for nodes
  T_vec = linspace<vec>(T_start, T_end, num_threads+1);
  IsingModel2D model; // initate class object;

  omp_set_num_threads(num_threads);
  # pragma omp parallel for default(shared) private (temps_i);

  for (temps_i = 0; temps_i < num_threads;temps_i++){
    T_start = T_vec(temps_i);
    T_end = T_vec(temps_i+1);
    model.init(L, T_start,T_end, n_T, MC);
    sol = model.solve(save_over_cycles);
    cout << sol;
  }

  //IsingModel2D model;
  //model.init(L, T_start,T_end, n_T, MC);
  //model.solve();


  //Catch::Session().run();

}
