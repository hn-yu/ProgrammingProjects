Now we've added a new member function to the Molecule class:
```c++
// Computes the angle between planes a-b-c and b-c-d
double Molecule::torsion(int a, int b, int c, int d)
{
  double eabc_x = (unit(1,b,a)*unit(2,b,c) - unit(2,b,a)*unit(1,b,c));
  double eabc_y = (unit(2,b,a)*unit(0,b,c) - unit(0,b,a)*unit(2,b,c));
  double eabc_z = (unit(0,b,a)*unit(1,b,c) - unit(1,b,a)*unit(0,b,c));

  double ebcd_x = (unit(1,c,b)*unit(2,c,d) - unit(2,c,b)*unit(1,c,d));
  double ebcd_y = (unit(2,c,b)*unit(0,c,d) - unit(0,c,b)*unit(2,c,d));
  double ebcd_z = (unit(0,c,b)*unit(1,c,d) - unit(1,c,b)*unit(0,c,d));

  double exx = eabc_x * ebcd_x;
  double eyy = eabc_y * ebcd_y;
  double ezz = eabc_z * ebcd_z;

  double tau = (exx + eyy + ezz)/(sin(angle(a,b,c)) * sin(angle(b,c,d)));

  if(tau < -1.0) tau = acos(-1.0);
  else if(tau > 1.0) tau = acos(1.0);
  else tau = acos(tau);

  // Compute the sign of the torsion 
  double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
  double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
  double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
  double norm = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
  double dot = 0.0;
  if(norm > 1e-14) {
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;
    dot = cross_x*unit(0,b,c)+cross_y*unit(1,b,c)+cross_z*unit(2,b,c);
  }
  double sign = 1.0;
  if(dot < 0.0) sign = -1.0;

  return tau*sign;
}
```

And we use the new function in the code as follows:

```c++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include "molecule.h"

using namespace std;

int main()
{
  Molecule mol("geom.dat", 0);

  cout << "Number of atoms: " << mol.natom << endl;
  cout << "Input Cartesian coordinates:\n";
  mol.print_geom();

  cout << "Interatomic distances (bohr):\n";
  for(int i=0; i < mol.natom; i++)
    for(int j=0; j < i; j++)
      printf("%d %d %8.5f\n", i, j, mol.bond(i,j));

  constexpr double distance_cutoff = 4.0;
  cout << "\nBond angles (both adjacent distances < 4.0 bohr):\n";
  // j is the central atom.  i < k removes only the i-j-k / k-j-i duplicate;
  // it does not constrain which atom may be the center.
  for(int j=0; j < mol.natom; j++) {
    for(int i=0; i < mol.natom; i++) {
      if(i == j) continue;
      for(int k=i+1; k < mol.natom; k++) {
        if(k == j) continue;
        if(mol.bond(i,j) < distance_cutoff && mol.bond(j,k) < distance_cutoff)
          printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i,j,k)*(180.0/acos(-1.0)));
      }
    }
  }

  cout << "\nOut-of-Plane angles:\n";
  for(int i=0; i < mol.natom; i++) {
    for(int k=0; k < mol.natom; k++) {
      for(int j=0; j < mol.natom; j++) {
        for(int l=0; l < j; l++) {
          if(i!=j && i!=k && i!=l && j!=k && k!=l && mol.bond(i,k) < 4.0 && mol.bond(k,j) < 4.0 && mol.bond(k,l) < 4.0)
              printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.oop(i,j,k,l)*(180.0/acos(-1.0)));
        }
      }
    }
  }

  cout << "\nTorsional angles:\n\n";
  // Canonicalize only the central bond (j < k).  A global index ordering such
  // as i > j > k > l would make the result depend on arbitrary atom numbering.
  for(int j=0; j < mol.natom; j++) {
    for(int k=j+1; k < mol.natom; k++) {
      if(mol.bond(j,k) >= distance_cutoff) continue;
      for(int i=0; i < mol.natom; i++) {
        if(i == j || i == k || mol.bond(i,j) >= distance_cutoff) continue;
        for(int l=0; l < mol.natom; l++) {
          if(l == i || l == j || l == k || mol.bond(k,l) >= distance_cutoff) continue;
          // A torsion is undefined when either plane normal vanishes.
          if(fabs(sin(mol.angle(i,j,k))) < 1e-12 ||
             fabs(sin(mol.angle(j,k,l))) < 1e-12) continue;
          printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
        }
      }
    }
  }

  return 0;
}
```

The above code produces the following output when applied to the acetaldehyde test case:

```
Number of atoms: 7
Input Cartesian coordinates:
6       0.000000000000       0.000000000000       0.000000000000
6       0.000000000000       0.000000000000       2.845112131228
8       1.899115961744       0.000000000000       4.139062527233
1      -1.894048308506       0.000000000000       3.747688672216
1       1.942500819960       0.000000000000      -0.701145981971
1      -1.007295466862      -1.669971842687      -0.705916966833
1      -1.007295466862       1.669971842687      -0.705916966833
Interatomic distances (bohr):
1 0  2.84511
2 0  4.55395
2 1  2.29803
3 0  4.19912
3 1  2.09811
3 2  3.81330
4 0  2.06517
4 1  4.04342
4 2  4.84040
4 3  5.87463
5 0  2.07407
5 1  4.05133
5 2  5.89151
5 3  4.83836
5 4  3.38971
6 0  2.07407
6 1  4.05133
6 2  5.89151
6 3  4.83836
6 4  3.38971
6 5  3.33994

Bond angles (both adjacent distances < 4.0 bohr):
 1- 0- 4 109.847056
 1- 0- 5 109.898406
 1- 0- 6 109.898406
 4- 0- 5 109.953682
 4- 0- 6 109.953682
 5- 0- 6 107.252646
 0- 1- 2 124.268308
 0- 1- 3 115.479341
 2- 1- 3 120.252351
 1- 2- 3  28.377448
 1- 3- 2  31.370201
 0- 4- 5  35.109529
 0- 4- 6  35.109529
 5- 4- 6  59.031048
 0- 5- 4  34.936789
 0- 5- 6  36.373677
 4- 5- 6  60.484476
 0- 6- 4  34.936789
 0- 6- 5  36.373677
 4- 6- 5  60.484476

Out-of-plane angles:
 0- 3- 1- 2  -0.000000
 0- 6- 4- 5  19.939726
 0- 6- 5- 4 -19.850523
 0- 5- 6- 4  19.850523
 1- 5- 0- 4  53.678778
 1- 6- 0- 4 -53.678778
 1- 6- 0- 5  54.977064
 2- 3- 1- 0   0.000000
 3- 2- 1- 0  -0.000000
 4- 5- 0- 1 -53.651534
 4- 6- 0- 1  53.651534
 4- 6- 0- 5 -54.869992
 4- 6- 5- 0  29.885677
 4- 5- 6- 0 -29.885677
 5- 4- 0- 1  53.626323
 5- 6- 0- 1 -56.277112
 5- 6- 0- 4  56.194621
 5- 6- 4- 0 -30.558964
 5- 4- 6- 0  31.064344
 6- 4- 0- 1 -53.626323
 6- 5- 0- 1  56.277112
 6- 5- 0- 4 -56.194621
 6- 5- 4- 0  30.558964
 6- 4- 5- 0 -31.064344

Torsional angles:

 4- 0- 1- 2   0.000000
 4- 0- 1- 3 180.000000
 5- 0- 1- 2 121.097586
 5- 0- 1- 3 -58.902414
 6- 0- 1- 2 -121.097586
 6- 0- 1- 3  58.902414
 1- 0- 4- 5 121.064344
 1- 0- 4- 6 -121.064344
 5- 0- 4- 6 117.871313
 6- 0- 4- 5 -117.871313
 1- 0- 5- 4 -121.033513
 1- 0- 5- 6 119.434473
 4- 0- 5- 6 -119.532014
 6- 0- 5- 4 119.532014
 1- 0- 6- 4 121.033513
 1- 0- 6- 5 -119.434473
 4- 0- 6- 5 119.532014
 5- 0- 6- 4 -119.532014
 0- 1- 2- 3 180.000000
 0- 1- 3- 2 180.000000
 0- 4- 5- 6  36.366799
 6- 4- 5- 0 -36.366799
 0- 4- 6- 5 -36.366799
 5- 4- 6- 0  36.366799
 0- 5- 6- 4  34.930266
 4- 5- 6- 0 -34.930266

```
