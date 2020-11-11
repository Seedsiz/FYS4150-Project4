#include "catch.hpp"
#include "montecarlo.hpp"
#include<iostream>

TEST_CASE("Testing expectation values") {
  IsingModel2D mysolver;
  int L = 2;
  int MC = 100;
  double T = 1.0;
  mysolver.init(L, T, MC);
  vec exp_val;
  exp_val = mysolver.solve();

  double exp_E_num = exp_val(0);
  double exp_M_num = exp_val(1);
  double mean_abs_M_num = exp_val(2);
  double Cv_num = exp_val(0);
  double xi_num = exp_val(0);

  vec num_val = vec(4);
  num_val(0) = exp_E_num;
  num_val(1) = mean_abs_M_num;
  num_val(2) = Cv_num;
  num_val(3) = xi_num;

  double z = 12 + 4*cosh((double) (8/T));
  double exp_E = -((double) (32/z))*sinh((double) 8/T);
  double mean_abs_M = (double) (3/2);
  double Cv = ((double) 256/(T*T*z))*(cosh((double) 8/T) - (4/z)*sinh((double) 8/T)*sinh((double) 8/T));
  double xi = ((double) 32/(T*z))*(exp((double) 8/T) + 1);

  vec exact = vec(4);
  exact(0) = exp_E;
  exact(1) = mean_abs_M;
  exact(2) = Cv;
  exact(3) = xi;

  double tol = 1E-04;

  for(int i = 0; i < 4; i++) {
    cout << num_val[i] << " " << exact[i] << "\n";
    //REQUIRE((num_val[i] - exact[i]) < tol);
  }

}
