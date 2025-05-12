/**
# Global polynomial reconstruction

We find the polynomial that describes the variation of the data long
the height ($z$).

$$f(z) = \sum_{n = 0}^{nl - 1} a_n z^n$$

For `nl` layers with depth $h_n$ ranging from $z_{n}$ to $z_{n + 1}$
and averaged scalar value $F_n$,  we have `nl` equations for the `nl`
unknows in the $\mathbf{a}_i$ array of $a_n$ values.

$$F_{\mathrm{old}} = A \mathbf{a_i},$$

$$\left(  \begin{array}{c} 
F_0 h_0 & \\ 
F_1 h_1&\\ 
\vdots &\\ 
 F_{nl-1} h_{nl - 1} & 
\end{array} \right)
 = 
\left( \begin{array}{cccc} 
z_1 - z_0 & \frac{z_1^2 - z_0^2}{2} & ... & \frac{z_1^{nl - 1} - z_0^{nl -1}}{nl - 1} \\ 
z_2 - z_1 & \frac{z_2^2 - z_1^2}{2} & ... & \frac{z_2^{nl - 1} - z_1^{nl -1}}{nl - 1} \\  
\vdots & & \ddots &\\
z_{nl} - z_{nl - 1} & \frac{z_{nl}^2 - z_{nl - 1}^2}{2} & ... & \frac{z_{nl}^{nl - 1} - z_{nl - 1}^{nl - 1}}{nl - 1}  
\end{array}\right) 
\left( \begin{array}{c}
a_0 & \\
a_1 &\\
\vdots & \\
a_{n-1} &
\end{array}\right)
$$
We can reintegrate these to any new mapping using these coefficients.

$$F_{\mathrm{new}} = B \mathbf{a}_i = B A^{-1} F_{\mathrm{old}}.$$

Where $B$ is the matrix based on the desired heights after mapping ($zn$),

$$ 
B = \left( \begin{array}{cccc} 
zn_1 - zn_0 & \frac{zn_1^2 - zn_0^2}{2} & ... & \frac{zn_1^{nl - 1} - zn_0^{nl -1}}{nl - 1} \\ 
zn_2 - zn_1 & \frac{zn_2^2 - zn_1^2}{2} & ... & \frac{zn_2^{nl - 1} - zn_1^{nl -1}}{nl - 1} \\  
\vdots & & \ddots &\\
zn_{nl} - zn_{nl - 1} & \frac{zn_{nl}^2 - zn_{nl - 1}^2}{2} & ... & \frac{zn_{nl}^{nl - 1} - zn_{nl - 1}^{nl - 1}}{nl - 1}  
\end{array}\right) 
$$

As such, we need a matrix invertor. Using Gauss--Jordan elimination.
 */
void inverse (int n, double A[n][n], double Ai[n][n]) {
  // Start with identity
  for (int i = 0; i < n ; i++)
    for (int j = 0; j < n; j++) 
      Ai[i][j] = i == j ? 1 : 0;
  // Repeated elimination of the Collumns
  for(int i = 0 ; i < n; i++) {
    assert(A[i][i] != 0);
    for(int j = 0; j < n; j++) {
      if(i != j) {
	double ratio = A[j][i]/A[i][i];
	for (int k = 0; k < n; k++) {
	  A[j][k]  -= ratio*A[i][k];
	  Ai[j][k] -= ratio*Ai[i][k];
	}
      }
    }
  }
  //scale rows
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      Ai[i][j] = Ai[i][j]/A[i][i];
  }
}
/**
A square matrix multiplication function for $AB = C$.
 */
void matmul (int n, double A[n][n], double B[n][n], double C[n][n]) {
  for (int i = 0; i < n; i++) 
    for (int j = 0; j < n; j++) {
      C[i][j] = 0.;
      for (int k = 0; k < n; k++) 
	C[i][j] += A[k][i]*B[j][k];
    }
}
/**
And a integer power-raising function for positive `i`*/
double raiser (double a, int i) {
  return i <= 0 ? 1 : a*raiser (a, i - 1);
}
/**
## Implementation

The desired heigh fractions are stored in `beta`.
*/
double * beta = NULL;

event defaults (i = 0) {
  beta = malloc (nl*sizeof(double));
  for (int l = 0; l < nl; l++)
    beta[l] = 1./nl;
}
/**
Because remapping can be expensive ($\mathcal{O}(nl^3)$), we only do
the work if the relative change of a cell's height is above a certain
fraction (`rem_frac`). The function returns the number of remapped
cells.
 */
double rem_frac = 0.1;

int remapper (scalar h, scalar * tracers) {
  int nvar = list_len(tracers);
  int remaps = 0;
  foreach(reduction (+:remaps)) {
    double H = 0.;
    foreach_layer()
      H += h[];
    // Do we bother to continue? 
    bool co = false;
    foreach_layer()
      if (fabs(beta[point.l]*H - h[])/h[] > rem_frac)
	co = true;
    if (co) {
      remaps++;
      int npos = nl + 1;
      double A[nl][nl], Ainv[nl][nl], B[nl][nl];
      double zpos[npos], znew[npos];
      double F[nl][nvar];
      zpos[0] = znew[0] = 0.;
      // Compute new heights, find matrices A, B and array F
      foreach_layer() {
	zpos[point.l + 1] = zpos[point.l] + h[];
	int it = 0;
	for (scalar s in tracers)
	  F[point.l][it++]     = s[]*h[];
	h[] = H*beta[point.l];
	znew[point.l + 1] = znew[point.l] + h[];
	for (int i = 0; i < nl; i++) {
	  A[i][point.l] = (raiser(zpos[point.l + 1], i + 1.) -
			   raiser(zpos[point.l], i+ 1.))/(i + 1.);
	  B[i][point.l] = (raiser(znew[point.l + 1], i + 1.) -
			   raiser(znew[point.l], i + 1.))/(i + 1.);
	}
      }
      // Invert and multiply with $B$
      inverse (nl, A, Ainv);
      matmul (nl, B, Ainv, A);
      // Recompte data.
      foreach_layer() {
	int it = 0;
	for (scalar s in tracers) {
	  double fa = 0;
	  for (int i = 0; i < nl; i++) 
	    fa += A[point.l][i]*F[i][it];
	  s[] = fa/h[];
	  it++;
	}
      }
    }
  }
  boundary ({h});
  boundary (tracers);
  return remaps;
}

/**
The remapping maybe applied every timestep
 */
int nr_remaps; 
event remap (i++) {
  if (nl > 1)
    nr_remaps = remapper (h, tracers);
}

event cleanup (t = end) {
  free (beta); beta = NULL;
}

/**
### Notes: 

* The remapping is conservative.
* Fitting a polynomial is likely to introduce osscilations at the edged
 of the domain.
* Smooth variation of the layers' depths is not enforced 

 */
