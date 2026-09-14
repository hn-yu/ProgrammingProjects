#include <bits/stdc++.h>
#include "molecule.h"

int main()
{
  Molecule M("./input/h2o_geom.txt");
  M.load_Hessian("./input/h2o_hessian.txt");
  cout << M.geom << endl;
  cout << M.hessian << endl;
};