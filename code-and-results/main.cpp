#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "montecarlo.hpp"
#include <iostream>
#include <armadillo>

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
  /*
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

  //Tryout random generator
  //MonteCarlo mysolver;
  //mysolver.initialize(L,T);
  //mysolver.draw_index();
  //mysolver.draw_acceptance();
  IsingModel2D model;
  model.init(L, T_start,T_end, n_T, MC);
  model.solve();
  */

  Catch::Session().run();

}
