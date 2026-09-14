#include <eigen3/Eigen/Dense>
#include <bits/stdc++.h>
using namespace std;

class Molecule
{
public:
  int n_atoms;
  Eigen::VectorXd z;
  Eigen::MatrixXd geom;

  Molecule(string filename)
  {
    ifstream file(filename);
    file >> n_atoms;
    geom.resize(n_atoms, 3);
    z.resize(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
    {
      file >> z[i] >> geom(i,0) >> geom(i,1) >> geom(i,2);
    }
  }

  void load_Hessian(string filename)
  {
    ifstream file(filename);
    file >> n_atoms;
    Eigen::MatrixXd hessian;

    for (int i = 0; i < n_atoms; ++i)
    {
      file >> hessian(i,0) >> hessian(i,1) >> hessian(i,2);
    }
  }
};

int main()
{
  Molecule M("./input/h2o_geom.txt");
  cout << M.geom << endl;
  cout << M.geom(0,2) << endl;
  cout << M.z << endl;
};