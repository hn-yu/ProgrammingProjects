#include "molecule.h"
namespace fs = std::filesystem;
typedef Eigen::MatrixXd mat;
using namespace std;

mat read_one_elec_integrals(const fs::path &p, const int n_orbitals)
{
  ifstream file(p);
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n_orbitals, n_orbitals);
  int i, j;
  double val;
  while (file >> i >> j >> val)
  {
    m(i - 1, j - 1) = val;
    m(j - 1, i - 1) = val;
  }
  return m;
}

int pair_index(int i, int j)
{
  if (j > i)
    swap(i, j);
  return (i * (i + 1)) / 2 + j;
}

int eri_index(int i, int j, int k, int l)
{
  return pair_index(pair_index(i, j), pair_index(k, l));
}

vector<double> read_two_elec_integrals(const fs::path &p, const int n_orbitals)
{
  ifstream file(p);
  const int n_pairs = n_orbitals * (n_orbitals + 1) / 2;
  const int n = n_pairs * (n_pairs + 1) / 2;
  vector<double> eri(n, 0.0);

  int i, j, k, l;
  double val;
  while (file >> i >> j >> k >> l >> val)
  {
    i--;
    j--;
    k--;
    l--;
    eri[eri_index(i, j, k, l)] = val;
  }
  return eri;
}

// make S^{1/2} matrix
mat inv_s_1_2(mat S)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
  Eigen::VectorXd val = solver.eigenvalues();
  for (int i = 0; i < val.size(); i++)
  {
    val[i] = 1 / sqrt(val[i]);
  }
  mat S_inv_1_2 = solver.eigenvectors() * val.asDiagonal() * solver.eigenvectors().transpose();
  return S_inv_1_2;
}

int main()
{
  Molecule M("./input/h2o/STO-3G/geom.dat");
  ifstream nuc("./input/h2o/STO-3G/s.dat");
  int n_basis = 7;
  double E_nuc;
  nuc >> E_nuc;
  mat S = read_one_elec_integrals("./input/h2o/STO-3G/s.dat", n_basis);
  mat T = read_one_elec_integrals("./input/h2o/STO-3G/t.dat", n_basis);
  mat V = read_one_elec_integrals("./input/h2o/STO-3G/v.dat", n_basis);
  mat H = T + V;
  cout << H << endl;
  auto ERI = read_two_elec_integrals("./input/h2o/STO-3G/eri.dat", 7);
  cout << ERI[0] << endl;
  auto X = inv_s_1_2(S);
  cout << X.transpose() * S * X << endl;

  // Step 5：Build the Initial Guess Density
  mat F = H;
  mat F_prime = X.transpose() * H * X;
  Eigen::SelfAdjointEigenSolver<mat> solver(F_prime);
  auto eps = solver.eigenvalues();
  auto C_prime = solver.eigenvectors();
  auto C = X * C_prime;
  cout << eps << endl;
  int n_electrons = M.z.sum();
  mat C_occ = C.leftCols(n_electrons / 2);
  mat P = C_occ * C_occ.transpose();

  // Step 6 : initial SCF energy
  // E_total = E_elec + E_nuc
  double E_elec = 0;
  for (int mu = 0; mu < n_basis; mu++)
  {
    for (int nu = 0; nu < n_basis; nu++)
    {
      E_elec += P(mu, nu) * H(mu, nu);
    }
  }
  E_elec *= 2;
  // =−125.8420774
  cout << E_elec << endl;

  // step 7: build Fock Matrix

  for (int mu = 0; mu < n_basis; ++mu)
  {
    for (int nu = 0; nu < n_basis; ++nu)
    {
      for (int lambda = 0; lambda < n_basis; ++lambda)
      {
        for (int sigma = 0; sigma < n_basis; ++sigma)
        {
          F(mu, nu) += P(lambda, sigma) * (2 * ERI[eri_index(mu, nu, lambda, sigma)] - ERI[eri_index(mu, lambda, nu, sigma)]) ;
        }
      }
    }
  }
  
  cout << F(0, 0) << endl;
}