To print each torsional angle only once, canonicalize the **central bond** rather than imposing a global ordering on all four atom indices. A rule such as `j < i`, `k < j`, `l < k` depends on the arbitrary atom numbering and can omit valid i-j-k-l chains.

One robust pattern is to require `j < k` for the central bond, then choose any distinct atom `i` close to `j` and any distinct atom `l` close to `k`:

```c++
for(int j=0; j < mol.natom; j++) {
  for(int k=j+1; k < mol.natom; k++) {
    if(mol.bond(j,k) >= distance_cutoff) continue;
    for(int i=0; i < mol.natom; i++) {
      if(i == j || i == k || mol.bond(i,j) >= distance_cutoff) continue;
      for(int l=0; l < mol.natom; l++) {
        if(l == i || l == j || l == k || mol.bond(k,l) >= distance_cutoff) continue;
        // i-j-k-l is now a unique candidate for this central bond.
      }
    }
  }
}
```

As in the bond-angle exercise, the 4.0 bohr cutoff is only a simple proximity filter, not a general chemical bond definition. Also skip cases where either adjacent bond angle is 0 or 180 degrees, because the corresponding plane normal and torsion are undefined.
