#include "molecule.h"
#include <eigen3/Eigen/Core>
namespace fs = std::filesystem;
typedef Eigen::MatrixXd mat;
using namespace std;

mat read_one_elec_integrals(const fs::path &p, const int n_orbitals) {
  ifstream file(p);
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n_orbitals, n_orbitals);
  int i, j;
  double val;
  while (file >> i >> j >> val) {
    m(i - 1, j - 1) = val;
    m(j - 1, i - 1) = val;
  }
  return m;
}

struct OrbitalSolution {
  Eigen::VectorXd eps;
  mat C;
};

int pair_index(int i, int j) {
  if (j > i)
    swap(i, j);
  return (i * (i + 1)) / 2 + j;
}

int eri_index(int i, int j, int k, int l) {
  return pair_index(pair_index(i, j), pair_index(k, l));
}

vector<double> read_two_elec_integrals(const fs::path &p,
                                       const int n_orbitals) {
  ifstream file(p);
  const int n_pairs = n_orbitals * (n_orbitals + 1) / 2;
  const int n = n_pairs * (n_pairs + 1) / 2;
  vector<double> eri(n, 0.0);

  int i, j, k, l;
  double val;
  while (file >> i >> j >> k >> l >> val) {
    i--;
    j--;
    k--;
    l--;
    eri[eri_index(i, j, k, l)] = val;
  }
  return eri;
}

// make S^{1/2} matrix
mat inv_s_1_2(mat S) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
  Eigen::VectorXd val = solver.eigenvalues();
  for (int i = 0; i < val.size(); i++) {
    val[i] = 1 / sqrt(val[i]);
  }
  mat S_inv_1_2 = solver.eigenvectors() * val.asDiagonal() *
                  solver.eigenvectors().transpose();
  return S_inv_1_2;
}

mat build_density(const mat &C, int n_occ) {
  mat C_occ = C.leftCols(n_occ);
  return C_occ * C_occ.transpose();
}

OrbitalSolution solve_fock(mat F, mat X) {
  mat F_prime = X.transpose() * F * X;
  Eigen::SelfAdjointEigenSolver<mat> solver(F_prime);
  mat C = X * solver.eigenvectors();
  return {solver.eigenvalues(), C};
}

mat build_fock(const mat &H, const mat &P, const vector<double> &ERI) {
  mat F = H;
  int n_basis = H.cols();
  for (int mu = 0; mu < n_basis; ++mu) {
    for (int nu = 0; nu < n_basis; ++nu) {
      for (int lambda = 0; lambda < n_basis; ++lambda) {
        for (int sigma = 0; sigma < n_basis; ++sigma) {
          F(mu, nu) +=
              P(lambda, sigma) * (2 * ERI[eri_index(mu, nu, lambda, sigma)] -
                                  ERI[eri_index(mu, lambda, nu, sigma)]);
        }
      }
    }
  }
  return F;
};

int main() {
  Molecule M("./input/h2o/STO-3G/geom.dat");
  ifstream nuc("./input/h2o/STO-3G/enuc.dat");
  int n_basis = 7;
  double E_nuc;
  nuc >> E_nuc;
  mat S = read_one_elec_integrals("./input/h2o/STO-3G/s.dat", n_basis);
  mat T = read_one_elec_integrals("./input/h2o/STO-3G/t.dat", n_basis);
  mat V = read_one_elec_integrals("./input/h2o/STO-3G/v.dat", n_basis);
  mat H = T + V;
  auto ERI = read_two_elec_integrals("./input/h2o/STO-3G/eri.dat", 7);
  auto X = inv_s_1_2(S);
  // cout << X.transpose() * S * X << endl;

  // Step 5：Build the Initial Guess Density

  auto solution = solve_fock(H, X);
  auto eps = solution.eps;
  auto C = solution.C;
  int n_electrons = M.z.sum();


  // we need a loop to iter to convergence

  cout << "iter    E_elec      E_total      ||dP||" << endl;
  mat P = build_density(C, n_electrons / 2);

  for (int i=0; i<100; i++) {
    mat F = build_fock(H, P, ERI);
    double E_elec = 0;
    for (int mu = 0; mu < n_basis; mu++) {
      for (int nu = 0; nu < n_basis; nu++) {
        E_elec += P(mu, nu) * (H(mu, nu) + F(mu, nu));
      }
    }

    auto solution = solve_fock(F, X);
    mat C_new = solution.C;
    mat P_new = build_density(C_new, n_electrons/2);
    double delta = (P-P_new).norm();
    if (delta < 1e-8) {
    break;
    }
    P = P_new;
    cout << i << "    "
     << E_elec << "    "
     << E_elec + E_nuc << "    "
     << delta << endl;
  }
}
