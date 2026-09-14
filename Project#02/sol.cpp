#include <bits/stdc++.h>
#include "molecule.h"

int main()
{
  Molecule M("./input/h2o_geom.txt");
  M.load_Hessian("./input/h2o_hessian.txt");
  cout << M.geom << endl;
  cout << "original hessian:" << endl << M.hessian << endl;
  M.weight_Hessian();
  cout << "weighted hessian" << endl << M.hessian << endl;
  auto solver = M.solve_Hessian();
  cout << solver.eigenvalues() << endl;

  for (auto a: solver.eigenvalues())
  {
    if (a<1e-8) a=0;
    a = std::sqrt(a) * 5140.5;
    cout << a << std::endl;
  }
};