/**
# A solver for Hessenberg systems

An [Hessenberg
matrix](https://en.wikipedia.org/wiki/Hessenberg_matrix) is an "almost
triangular" matrix i.e. the sum of a triangular matrix and a
tridiagonal matrix.

The function below solves $Hx=b$ where $H$ is an upper Hessenberg
matrix of rank $n$. The right-hand side $b$ is given as vector $x$ and
is replaced by the solution. $H$ is given as a one dimensional array
where each matrix element is indexed as $H_{ij} = H[in+j]$.

## References

~~~bib
@book{henry1994,
  title={The shifted Hessenberg system solve computation},
  author={Henry, Greg},
  year={1994},
  publisher={Cornell Theory Center, Cornell University},
  url={https://pdfs.semanticscholar.org/df75/8d16317f246ac4049a1569b6f56510a4add7.pdf}
}
~~~
*/

static inline void givens (double x, double y, double * c, double * s)
{
#if 0
  #define sign2(a,b) ((b) > 0. ? ((a) > 0. ? (a) : -(a)) : ((a) > 0. ? -(a) : (a)))

  if (x == 0. && y == 0.)
    *c = 1., *s = 0.;
  else if (fabs(y) > fabs(x)) {
    double t = x/y;
    x = sqrt(1. + t*t);
    *s = - sign2(1./x, y);
    *c = t*(*s);
  }
  else {
    double t = y/x;
    y = sqrt2(1. + t*t);
    *c = sign2(1./y, x);
    *s = - t*(*c);
  }
#else
  double t = sqrt (sq(x) + sq(y));
  *c = x/t, *s = -y/t;
#endif
}

void solve_hessenberg (double * H, double * x, int n)
{
  double v[n], c[n], s[n];
  for (int i = 0; i < n; i++)
    v[i] = H[n*(i + 1) - 1];
  for (int k = n - 1; k >= 1; k--) {
    double a = H[k*n + k - 1];
    givens (v[k], a, &c[k], &s[k]);
    x[k] /= c[k]*v[k] - s[k]*a;
    double ykck = x[k]*c[k], yksk = x[k]*s[k];
    for (int l = 0; l <= k - 2; l++) {
      a = H[l*n + k - 1];
      x[l] -= ykck*v[l] - yksk*a;
      v[l] = c[k]*a + s[k]*v[l];
    }
    a = H[(k - 1)*n + k - 1];
    x[k-1] -= ykck*v[k-1] - yksk*a;
    v[k-1] = c[k]*a + s[k]*v[k-1];
  }
  double tau1 = x[0]/v[0];
  for (int k = 1; k < n; k++) {
    double tau2 = x[k];
    x[k-1] = c[k]*tau1 - s[k]*tau2;
    tau1 = c[k]*tau2 + s[k]*tau1;
  }
  x[n-1] = tau1;
}
