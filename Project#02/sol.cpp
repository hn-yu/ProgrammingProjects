#include <eigen3/Eigen/Dense>
#include <bits/stdc++.h>
using namespace std;

class Molecule
{
public:
  int n_atoms;
  Eigen::VectorXd z;
  Eigen::MatrixXd geom;
  Eigen::MatrixXd hessian;

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
    hessian.resize(3*n_atoms, 3*n_atoms);

    for (int i = 0; i < 3*n_atoms; ++i)
    {
      for (int j = 0; j < 3*n_atoms; ++j)
      {
        file >> hessian(i,j);
      }
    }
  }
};

int main()
{
  Molecule M("./input/h2o_geom.txt");
  cout << M.geom << endl;
  cout << M.geom(0,2) << endl;
  cout << M.z << endl;
  cout << M.hessian << endl; 
};