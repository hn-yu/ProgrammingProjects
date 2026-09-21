#include <bits/stdc++.h>
using namespace std;

struct Atom
{
  int Z;
  double x, y, z;
};

double masses[] = {
    0.0000000, // 0: dummy
    1.007825,  // H-1
    4.002603,  // He-4
    7.016004,  // Li-7
    9.012183,  // Be-9
    11.009305, // B-11
    12.000000, // C-12
    14.003074, // N-14
    15.994915, // O-16
    18.998403, // F-19
    19.992440, // Ne-20
    22.989770, // Na-23
    23.985042, // Mg-24
    26.981538, // Al-27
    27.976927, // Si-28
    30.973762, // P-31
    31.972071, // S-32
    34.968853, // Cl-35
    39.962383, // Ar-40
    38.963707, // K-39
    39.962591  // Ca-40
};

class Molecule
{
public:
  vector<Atom> va;
  int atom_number = 0;
  vector<double> pair_dists;
  vector<double> bond_angles;

  Molecule(const filesystem::path &p)
  {

    ifstream file(p);
    file >> atom_number;
    for (int i = 0; i < atom_number; i++)
    {
      Atom a;
      file >> a.Z >> a.x >> a.y >> a.z;
      va.push_back(a);
    };
    pair_dists.resize(atom_number * atom_number);
    bond_angles.resize(atom_number * atom_number * atom_number);
  };

  void calc_bond_length()
  {
    for (int i = 0; i < atom_number; i++)
    {
      for (int j = 0; j < i; j++)
      {
        float dx = va[i].x - va[j].x;
        float dy = va[i].y - va[j].y;
        float dz = va[i].z - va[j].z;
        float atomic_length = sqrt(dx * dx + dy * dy + dz * dz);
        cout << i << " " << j << " " << atomic_length << '\n';
        pair_dists[i * atom_number + j] = atomic_length;
        pair_dists[j * atom_number + i] = atomic_length;
      }
    }
  }

  void calc_bond_angle()
  {
    for (int j = 0; j < atom_number; j++)
    {
      for (int k = 0; k < atom_number; k++)
      {
        if (pair_dists[j * atom_number + k] > 4.0 || j == k)
          continue;
        for (int i = 0; i < k; i++)
        {
          if (pair_dists[i * atom_number + j] > 4.0 || i == j)
            continue;
          double r_ji = pair_dists[i * atom_number + j];
          double r_jk = pair_dists[j * atom_number + k];
          double cos_angle_1 = ((va[j].x - va[i].x) / r_ji) * ((va[j].x - va[k].x) / r_jk);
          double cos_angle_2 = ((va[j].y - va[i].y) / r_ji) * ((va[j].y - va[k].y) / r_jk);
          double cos_angle_3 = ((va[j].z - va[i].z) / r_ji) * ((va[j].z - va[k].z) / r_jk);
          double angle = acos(cos_angle_1 + cos_angle_2 + cos_angle_3) * 180.0 / numbers::pi;
          cout << i << " " << j << " " << k << " " << angle << "\n";
          bond_angles[i * atom_number * atom_number + j * atom_number + k] = angle;
        }
      }
    }
  }

  double get_bond_angle(int i, int j, int k)
  {
    if (i > k)
      swap(i, k);
    return bond_angles[i * atom_number * atom_number + j * atom_number + k];
  };

  void calc_oop_angle()
  {
    // treat k as central atom
    //   i
    //   |
    //   k
    //  / \
    // j   l
    for (int k = 0; k < atom_number; k++)
    {
      for (int j = 0; j < atom_number; j++)
      {
        if (j == k || pair_dists[j * atom_number + k] > 4.0)
          continue;
        for (int l = 0; l < j; l++)
        {
          if (l == k || pair_dists[l * atom_number + k] > 4.0)
            continue;
          for (int i = 0; i < atom_number; i++)
          {
            if (i == j || i == k || i == l || pair_dists[i * atom_number + k] > 4.0)
              continue;
            double kj[3] = {
                (va[k].x - va[j].x) / pair_dists[k * atom_number + j],
                (va[k].y - va[j].y) / pair_dists[k * atom_number + j],
                (va[k].z - va[j].z) / pair_dists[k * atom_number + j]};
            double kl[3] = {
                (va[k].x - va[l].x) / pair_dists[k * atom_number + l],
                (va[k].y - va[l].y) / pair_dists[k * atom_number + l],
                (va[k].z - va[l].z) / pair_dists[k * atom_number + l]};
            double ki[3] = {
                (va[k].x - va[i].x) / pair_dists[k * atom_number + i],
                (va[k].y - va[i].y) / pair_dists[k * atom_number + i],
                (va[k].z - va[i].z) / pair_dists[k * atom_number + i]};
            double sin_jkl = sin(get_bond_angle(j, k, l) / 180.0 * numbers::pi);
            if (abs(sin_jkl) < 1e-12)
              continue;
            double s1 = kj[1] * kl[2] - kj[2] * kl[1];
            double s2 = kj[2] * kl[0] - kj[0] * kl[2];
            double s3 = kj[0] * kl[1] - kj[1] * kl[0];
            double sin_ijkl = (s1 * ki[0] + s2 * ki[1] + s3 * ki[2]) / sin_jkl;
            sin_ijkl = clamp(sin_ijkl, -1.0, 1.0);
            double angle = asin(sin_ijkl) * 180.0 / numbers::pi;
            cout << i << " " << j << " " << k << " " << l << " " << angle << "\n";
          }
        }
      }
    }
  };

  void calc_dihedral_angles()
  // i
  //  \
  //   j ===== k
  //            \
  //             l
  {
    for (int j = 0; j < atom_number; j++)
    {
      for (int k = j + 1; k < atom_number; k++)
      {
        if (pair_dists[j * atom_number + k] > 4.0)
          continue;
        for (int i = 0; i < atom_number; i++)
        {
          if (i == j || i == k || pair_dists[i * atom_number + j] > 4.0)
            continue;
          for (int l = 0; l < atom_number; l++)
          {
            if (l == i || l == j || l == k || pair_dists[l * atom_number + k] > 4.0)
            {
              continue;
            }
            // for all valid 4-tuple ijkl
            double ij[3] = {
                (va[i].x - va[j].x) / pair_dists[i * atom_number + j],
                (va[i].y - va[j].y) / pair_dists[i * atom_number + j],
                (va[i].z - va[j].z) / pair_dists[i * atom_number + j]};
            double jk[3] = {
                (va[j].x - va[k].x) / pair_dists[k * atom_number + j],
                (va[j].y - va[k].y) / pair_dists[k * atom_number + j],
                (va[j].z - va[k].z) / pair_dists[k * atom_number + j]};
            double kl[3] = {
                (va[k].x - va[l].x) / pair_dists[k * atom_number + l],
                (va[k].y - va[l].y) / pair_dists[k * atom_number + l],
                (va[k].z - va[l].z) / pair_dists[k * atom_number + l]};
            double cross_eij_ejk[3] =
                {
                    ij[1] * jk[2] - ij[2] * jk[1],
                    ij[2] * jk[0] - ij[0] * jk[2],
                    ij[0] * jk[1] - ij[1] * jk[0]};
            double cross_ejk_ekl[3] =
                {
                    jk[1] * kl[2] - jk[2] * kl[1],
                    jk[2] * kl[0] - jk[0] * kl[2],
                    jk[0] * kl[1] - jk[1] * kl[0]};
            double top = cross_eij_ejk[0] * cross_ejk_ekl[0] + cross_eij_ejk[1] * cross_ejk_ekl[1] + cross_eij_ejk[2] * cross_ejk_ekl[2];
            double sin_ijk = sin(get_bond_angle(i, j, k) / 180.0 * numbers::pi);
            double sin_jkl = sin(get_bond_angle(j, k, l) / 180.0 * numbers::pi);
            if (sin_ijk < 1e-12 || sin_jkl < 1e-12)
              continue;
            double cos_angle = top / (sin_ijk * sin_jkl);
            cos_angle = clamp(cos_angle, -1.0, 1.0);
            double angle = acos(cos_angle) * 180.0 / numbers::pi;
            cout << i << " " << j << " " << k << " " << l << " " << angle << endl;
          }
        }
      }
    }
  };

  void com_translation()
  {
    double mtotal = 0.0, mx = 0.0, my = 0.0, mz = 0.0;
    for (Atom &atom : va)
    {
      mx += masses[atom.Z] * atom.x;
      my += masses[atom.Z] * atom.y;
      mz += masses[atom.Z] * atom.z;
      mtotal += masses[atom.Z];
    }
    double x_cm = mx / mtotal, y_cm = my / mtotal, z_cm = mz / mtotal;
    cout << "Moment of inertia tensor: " << x_cm << " " << y_cm << " " << z_cm << endl;
    for (Atom &atom : va)
    {
      atom.x -= x_cm;
      atom.y -= y_cm;
      atom.z -= z_cm;
    }
  };
  
  void calc_moments_of_inertia()
  {
    ;
  }
  
};

int main()
{
  Molecule M("input/acetaldehyde.dat");
  M.calc_bond_length();
  M.calc_bond_angle();
  cout << "get the 0,6,5 bond angle from 5,6,0:\n";
  cout << M.get_bond_angle(5, 6, 0) << "\n";
  cout << "Out-of-plane angles: \n";
  M.calc_oop_angle();
  cout << "Torsional angles:" << endl;
  M.calc_dihedral_angles();
  for (Atom &atom : M.va)
  {
    cout << atom.x << " " << atom.y << " " << atom.z << " " << endl;
  }
  M.com_translation();
  for (Atom &atom : M.va)
  {
    cout << atom.x << " " << atom.y << " " << atom.z << " " << endl;
  }
}