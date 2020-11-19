#include "catch.hpp"
#include "montecarlo.hpp"
#include<iostream>

TEST_CASE("Testing expectation values") {
  IsingModel2D mysolver;
  int calibration = 20000;
  int L = 2;
  int MC = 1e7;
  double T = 1.0;
  double T2 = 1.0;
  int n = 1;
  double B = 1./T;
  int rank = 0;
  mysolver.init(L, T,T2, n, MC, rank);
  bool save_over_cycles = false;
  vec exp_val;
  exp_val = mysolver.solve(save_over_cycles, calibration);

  double exp_E_num = exp_val(0);
  double exp_M_num = exp_val(1);
  double mean_abs_M_num = exp_val(2);
  double Cv_num = exp_val(3);
  double xi_num = exp_val(4);

  vec num_val = vec(4);
  num_val(0) = exp_E_num;
  num_val(1) = mean_abs_M_num;
  num_val(2) = Cv_num;
  num_val(3) = xi_num;

  double z = 12. + 4*cosh(8*B);
  double exp_E = -(32./z)*sinh(8*B);
  double mean_abs_M = (8./z)*(exp(8*B) + 2);
  double Cv = ((double) 256/(T*T*z))*(cosh(8*B) - (4/z)*sinh(8*B)*sinh(8*B));
  double xi = ((double) 32/(T*z))*(exp(8*B) + 1);

  vec exact = vec(4);
  exact(0) = exp_E/(L*L);
  exact(1) = mean_abs_M/(L*L);
  exact(2) = Cv/(L*L);
  exact(3) = xi/(L*L);

  double tol = 1E-02;
  for(int i = 0; i < 4; i++) {
    cout << num_val[i] << " " << exact[i] << "\n";
    REQUIRE((num_val[i] - exact[i]) < tol);
  }

}
