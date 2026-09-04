Previously the [Hessian](../../Project%2302/hints/hint1.md) ([Project #2](../../Project%2302/)) and the one-electron quantities so far in [Project 3](../../Project%2303) have been provided in an order that makes reading them in using a loop structure trivial. 

```c++
  ...
  
  FILE *input;
  int a, b;
  double **S;
  
  ...
  
  input = fopen("s.dat", "r");
  for(int i=0; i < nao; i++) {
    for(int j=0; j <= i; j++) {
      fscanf(input, "%d %d %lf", &a, &b, &S[i][j]);
      S[j][i] = S[i][j];
    }
  }
  fclose(input);
```

However the two-electron integrals are not provided in this convenient ordering, and the file omits integrals whose value is zero. Initialize your packed TEI storage to zero, then use the explicit indices on each line to place every non-zero value. This avoids hard-wired loops that assume a specific ordering or line count.

```c++
  ...
  
  FILE *input;
  int i, j, k, l;
  double val;
  double *TEI;
  
  ...
  
  input = fopen("eri.dat", "r");
  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) == 5) {
  
  ...
  
    TEI[ijkl] = val;
  }
  fclose(input);
```
`fscanf` returns the number of successfully converted fields, so `== 5` accepts only complete ERI records and stops cleanly at end-of-file or malformed input.
