/**
# Harmonic decomposition 

The `harmonic_decomposition()` function performs a continuous,
least-square fit of the coefficients $Z$ and $(a_i,b_i)$ of the harmonic
fonction
$$
Z + \sum_i a_i \cos(\omega_i t) + b_i \sin(\omega_i t)
$$
to the scalar field `s`, with the $\omega_i$ given as input. 

If the optional argument `e` is given, it is used to store the
residual (i.e. the standard deviation) of the fit. */

attribute {
  struct {
    double * omega;
    scalar * A, * B, Z, E;
    double ** M, ** Mn;
    bool invertible; 
  } harmonic;
}

static double de (int n, double * ha, double * hx, double ** M)
{
  double xm = ha[2*n];
  double e = xm*(M[2*n][2*n]*xm - 2.*hx[2*n]);

  for (int i = 0; i < n; i++) {
    e += 2.*(ha[i]*(xm*M[i][2*n] - hx[i]) +
	     ha[n + i]*(xm*M[n + i][2*n] - hx[n + i]));
    for (int j = 0; j < n; j++)
      e += (ha[i]*ha[j]*M[j][i] + 
	    ha[n + i]*ha[n + j]*M[n + j][n + i] +
	    2.*ha[i]*ha[n + j]*M[n + j][i]);
  }
  return e;
}

void harmonic_decomposition (scalar s, double t, double * omega, scalar e = {-1})
{
  if (!s.harmonic.omega) {
    int n = 0;
    for (double * o = omega; *o > 0.; o++, n++);
    s.harmonic.omega = malloc(n*sizeof(double));
    memcpy (s.harmonic.omega, omega, n*sizeof(double));
    s.harmonic.M = (double **) matrix_new (2*n + 1, 2*n + 1, sizeof (double));
    s.harmonic.Mn = (double **) matrix_new (2*n + 1, 2*n + 1, sizeof (double));
    for (int i = 0; i < 2*n + 1; i++)  
      for (int j = 0; j < 2*n + 1; j++) {
	s.harmonic.M[i][j] = 0.;
	s.harmonic.Mn[i][j] = (i == j);
      }
    s.harmonic.E = e;
    scalar Z = new scalar;
    s.harmonic.Z = Z;
    s.harmonic.A = s.harmonic.B = NULL;
    for (int i = 0; i < n; i++) {
      scalar a = new scalar, b = new scalar;
      char name[80];
      snprintf (name, 79, "%s%da", s.name, i);
      free (a.name); a.name = strdup (name);
      snprintf (name, 79, "%s%db", s.name, i);
      free (b.name); b.name = strdup (name);      
      s.harmonic.A = list_append (s.harmonic.A, a);
      s.harmonic.B = list_append (s.harmonic.B, b);
    }
    reset (s.harmonic.A, 0.);
    reset (s.harmonic.B, 0.);
    s.harmonic.invertible = false;
  }
  
  scalar * A = s.harmonic.A, * B = s.harmonic.B;
  scalar Z = s.harmonic.Z, E = s.harmonic.E;
  double ** Mn = s.harmonic.Mn;
  int n = 0;
  for (scalar a in A)
    n++;
  
  double vsin[n], vcos[n];
  for (int i = 0; i < n; i++) {
    vsin[i] = sin (s.harmonic.omega[i]*t);
    vcos[i] = cos (s.harmonic.omega[i]*t);
  }

  double ** M = s.harmonic.M;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[i][j]         += vcos[j]*vcos[i];
      M[i][n + j]     += vsin[j]*vcos[i];
      M[n + i][j]     += vcos[j]*vsin[i];
      M[n + i][n + j] += vsin[j]*vsin[i];
    }
    M[i][2*n]     += vcos[i];
    M[n + i][2*n] += vsin[i];
  }
  for (int j = 0; j < n; j++) {
    M[2*n][j]     += vcos[j];
    M[2*n][n + j] += vsin[j];
  }
  M[2*n][2*n] += 1.;

  double ** iM = (double **) matrix_new (2*n + 1, 2*n + 1, sizeof (double));
  for (int i = 0; i < 2*n + 1; i++)
    for (int j = 0; j < 2*n + 1; j++)
      iM[i][j] = M[i][j];
  if (!matrix_inverse (iM, 2*n + 1, 1e-6)) {
    assert (!s.harmonic.invertible);
    foreach() {
      double x = s[];
      scalar a, b;
      int i = 0;
      for (a,b in A,B) {
	a[] += x*vcos[i];
	b[] += x*vsin[i];
	i++;
      }
      Z[] += x;
      if (E.i >= 0)
	E[] += x*x;
    }
  }
  else {
    foreach() {
      double x = s[], sx2 = 0.;

      /* A^n */
      double ha[2*n + 1];
      scalar a, b;
      int i = 0;
      for (a,b in A,B) {
	ha[i] = a[];
	ha[i + n] = b[];
	i++;
      }
      ha[2*n] = Z[];

      /* X^n = M^n.A^n */
      double hx[2*n + 1];
      for (int i = 0; i < 2*n + 1; i++) {
	hx[i] = 0.;
	for (int j = 0; j < 2*n + 1; j++)
	  hx[i] += Mn[i][j]*ha[j];
      }

      if (E.i >= 0) {
	if (s.harmonic.invertible)
	  sx2 = x*x + Mn[2*n][2*n]*E[] - de (n, ha, hx, Mn);
	else
	  sx2 = x*x + E[];
      }
  
      /* X^n+1 = X^n + Delta^n */
      for (int i = 0; i < n; i++) {
	hx[i]     += x*vcos[i];
	hx[i + n] += x*vsin[i];
      }
      hx[2*n] += x;

      /* A^n+1 = (M^n+1)^-1.X^n+1 */
      for (int i = 0; i < 2*n + 1; i++) {
	ha[i] = 0.;
	for (int j = 0; j < 2*n + 1; j++)
	  ha[i] += iM[i][j]*hx[j];
      }

      i = 0;
      for (a,b in A,B) {
	a[] = ha[i];
	b[] = ha[i + n];
	i++;
      }
      Z[] = ha[2*n];

      if (E.i >= 0)
	E[] = (sx2 + de (n, ha, hx, M))/M[2*n][2*n];      
    }
    s.harmonic.invertible = true;
    for (int i = 0; i < 2*n + 1; i++)
      for (int j = 0; j < 2*n + 1; j++)
	Mn[i][j] = M[i][j];
  }
  matrix_free (iM);
}

event cleanup (t = end)
{
  for (scalar s in all)
    if (s.harmonic.omega) {
      free (s.harmonic.omega); s.harmonic.omega = NULL;
      matrix_free (s.harmonic.M);
      matrix_free (s.harmonic.Mn);
      scalar Z = s.harmonic.Z;
      delete ({Z});
      delete (s.harmonic.A);
      free (s.harmonic.A);
      delete (s.harmonic.B);
      free (s.harmonic.B);
    }
}
