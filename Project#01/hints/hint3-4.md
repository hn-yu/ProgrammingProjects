Printing every central-angle combination can produce many geometrically uninteresting angles.  The examples in this project therefore use a **4.0 bohr proximity cutoff** for the two distances adjacent to the central atom.  This cutoff is only a simple printing heuristic; it is not a general chemical-bond criterion.

The important indexing rule is that the central atom must be allowed to have **any** index.  Do not enforce `i < j < k` (or `i > j > k`) while also treating `j` as the center, because that makes the result depend on atom numbering and omits valid central atoms.

A clean loop structure is:

```c++
constexpr double distance_cutoff = 4.0;

for(int j=0; j < mol.natom; j++) {          // central atom
  for(int i=0; i < mol.natom; i++) {
    if(i == j) continue;
    for(int k=i+1; k < mol.natom; k++) {    // i < k avoids i-j-k / k-j-i duplicates
      if(k == j) continue;
      if(R[i][j] < distance_cutoff && R[j][k] < distance_cutoff) {
        ...
      }
    }
  }
}
```

If you want **chemical** bond angles rather than cutoff-defined angle candidates, first construct a molecular connectivity graph using an appropriate bonding criterion, then choose pairs of neighbors around each central atom `j`.
