#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
#include "catch.hpp"
#include "montecarlo.hpp"
#include <iostream>

using namespace std;


void menu();

int main(int argc, char const *argv[]){
  menu();
  return 0;
}

void menu(){
  int L;
  int MC;
  double T;
  cout << "Enter integer number of spin particles for each axis:" << " ";
  cin >> L;
  cout << "Enter number of MC cycles:"  << " ";
  cin >> MC;
  cout << "Choose a temperature T:"  << " ";
  cin >> T;

  //Tryout random generator
  //MonteCarlo mysolver;
  //mysolver.initialize(L,T);
  //mysolver.draw_index();
  //mysolver.draw_acceptance();
  IsingModel2D model;
  model.init(L, T, MC);
  //model.solve();
}
