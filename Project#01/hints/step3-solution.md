Here again we've extended the Molecule class to compute the bond angles without storing them.  This also required the addition of a new function to compute the Cartesian unit vectors on the fly, whose declaration we add to the member functions of molecule.h:
```c++
#include <string>

using namespace std;

class Molecule
{
  public:
    int natom;
    int charge;
    int *zvals;
    double **geom;
    string point_group;

    void print_geom();
    void rotate(double phi);
    void translate(double x, double y, double z);
    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double torsion(int atom1, int atom2, int atom3, int atom4);
    double unit(int cart, int atom1, int atom2);

    Molecule(const char *filename, int q);
    ~Molecule();
};
```

Now the two new functions added to molecule.cc:
```c++
// Returns the value of the unit vector between atoms a and b
// in the cart direction (cart=0=x, cart=1=y, cart=2=z)
double Molecule::unit(int cart, int a, int b)
{
  return -(geom[a][cart]-geom[b][cart])/bond(a,b);
}

// Returns the angle between atoms a, b, and c in radians
double Molecule::angle(int a, int b, int c)
{
  double cosine = unit(0,b,a) * unit(0,b,c)
                + unit(1,b,a) * unit(1,b,c)
                + unit(2,b,a) * unit(2,b,c);
  if(cosine < -1.0) cosine = -1.0;
  if(cosine >  1.0) cosine =  1.0;
  return acos(cosine);
}
```

And finally, the code that makes use of the above: 

```c++
#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

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

  return 0;
}
```
Note that the value of &pi; is obtained with machine precision using `acos(-1.0)`.

The above code produces the following output for the acetaldehyde test case:

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
```
