// molecule.h
#include <eigen3/Eigen/Dense>
namespace fs = std::filesystem;
using namespace std;

constexpr std::array<double, 51> masses = {
    0.000000000, // 0  dummy

    1.007825032,  // 1  H-1
    4.002603254,  // 2  He-4
    7.016003437,  // 3  Li-7
    9.012183065,  // 4  Be-9
    11.00930536,  // 5  B-11
    12.000000000, // 6  C-12
    14.003074004, // 7  N-14
    15.994914620, // 8  O-16
    18.998403163, // 9  F-19
    19.992440176, // 10 Ne-20

    22.989769282, // 11 Na-23
    23.985041697, // 12 Mg-24
    26.981538530, // 13 Al-27
    27.976926535, // 14 Si-28
    30.973761998, // 15 P-31
    31.972071174, // 16 S-32
    34.968852683, // 17 Cl-35
    39.962383123, // 18 Ar-40
    38.963706486, // 19 K-39
    39.962590864, // 20 Ca-40

    44.955908280, // 21 Sc-45
    47.947941980, // 22 Ti-48
    50.943957040, // 23 V-51
    51.940506230, // 24 Cr-52
    54.938043910, // 25 Mn-55
    55.934936330, // 26 Fe-56
    58.933194290, // 27 Co-59
    57.935342410, // 28 Ni-58
    62.929597720, // 29 Cu-63
    63.929142010, // 30 Zn-64

    68.925573500, // 31 Ga-69
    73.921177760, // 32 Ge-74
    74.921594570, // 33 As-75
    79.916521800, // 34 Se-80
    78.918337600, // 35 Br-79
    83.911497728, // 36 Kr-84
    84.911789738, // 37 Rb-85
    87.905612500, // 38 Sr-88
    88.905840300, // 39 Y-89
    89.904697700, // 40 Zr-90

    92.906373000,  // 41 Nb-93
    97.905404820,  // 42 Mo-98
    96.906366700,  // 43 Tc-97
    101.904344100, // 44 Ru-102
    102.905498000, // 45 Rh-103
    105.903480400, // 46 Pd-106
    106.905091600, // 47 Ag-107
    113.903365100, // 48 Cd-114
    114.903878800, // 49 In-115
    119.902201600, // 50 Sn-120

};

class Molecule
{
public:
  int n_atoms;
  Eigen::VectorXi z;
  Eigen::MatrixXd geom;
  Eigen::MatrixXd hessian;

  Molecule(fs::path filename)
  {
    ifstream file(filename);
    file >> n_atoms;
    geom.resize(n_atoms, 3);
    z.resize(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
    {
      double d;
      file >> d >> geom(i, 0) >> geom(i, 1) >> geom(i, 2);
      z[i] = (int) d;
    }
  }

  void load_Hessian(fs::path filename)
  {
    ifstream file(filename);
    int n;
    file >> n;
    assert(n == n_atoms);
    hessian.resize(3 * n_atoms, 3 * n_atoms);

    for (int i = 0; i < 3 * n_atoms; ++i)
    {
      for (int j = 0; j < 3 * n_atoms; ++j)
      {
        file >> hessian(i, j);
      }
    }
  }

  void weight_Hessian()
  {
    for (int i = 0; i < 3 * n_atoms; ++i)
    {
      for (int j = 0; j < 3 * n_atoms; ++j)
      {
        hessian(i, j) /= sqrt(masses[z[i/3]] * masses[z[j/3]]);
      }
    }
  }
};