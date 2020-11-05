#include "montecarlo.hpp"
#include <iostream>

using namespace std;


void menu();

int main(int argc, char const *argv[]){
  menu();
  return 0;
}

void menu(){
  int L; double T;
  cout << "Enter integer number of spin particles for each axis:" << " ";
  cin >> L;
  cout << "Choose a temperature T:"  << " ";
  cin >> T;

  //Tryout random generator
  MonteCarlo mysolver;
  mysolver.initialize(L,T);
  mysolver.draw_index();
}
