static void eigsrt (double d[dimension],
		    double v[dimension][dimension])
{
  int k, j, i;
  double p;

  for (i = 0; i < dimension - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < dimension; j++)
      if (d[j] >= p) 
	p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < dimension; j++) {
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }
  }
}

#define ROTATE(a,i,j,k,l) {\
    g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

/**
 * eigenvalues:
 * @a: a symmetric matrix.
 * @d: a vector.
 * @v: another matrix.
 *
 * Fills @d (resp. @v) with the eigenvalues (resp. eigenvectors) of
 * matrix @a.
 */
void eigenvalues (double a[dimension][dimension],
		  double d[dimension],
		  double v[dimension][dimension])
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, b[dimension], z[dimension];

  for (ip = 0; ip < dimension; ip++) {
    for (iq = 0; iq < dimension; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < dimension; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < dimension - 1; ip++) {
      for (iq = ip + 1; iq < dimension; iq++)
	sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(dimension*dimension);
    else
      tresh = 0.0;
    for (ip = 0; ip < dimension - 1; ip++) {
      for (iq = ip + 1; iq < dimension; iq++) {
	g = 100.0*fabs (a[ip][iq]);
	if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
	    fabs(d[iq]) + g == fabs(d[iq]))
	  a[ip][iq] = 0.0;
	else if (fabs (a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if (fabs(h) + g == fabs(h))
	    t = a[ip][iq]/h;
	  else {
	    theta = 0.5*h/a[ip][iq];
	    t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt (1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for (j = 0; j <= ip - 1; j++)
	    ROTATE (a, j, ip, j, iq);
	  for (j = ip + 1; j <= iq - 1; j++)
	    ROTATE (a, ip, j, j, iq);
	  for (j = iq + 1; j < dimension; j++)
	    ROTATE(a, ip, j, iq, j);
	  for (j = 0; j < dimension; j++)
	    ROTATE(v, j, ip, j, iq);
	}
      }
    }
    for (ip = 0; ip < dimension; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  /* Too many iterations */
  for (i = 0; i < dimension; i++) {
    for (j = 0; j < dimension; j++)
      fprintf (stderr, "%10.3g ", a[i][j]);
    fprintf (stderr, "\n");
  }
  assert (false);
}

void lambda2 (const vector u, scalar l2)
{
  foreach() {
    double JJ[dimension][dimension];
    scalar s = u.x;
    int i = 0;
    foreach_dimension()
      JJ[0][i++] = center_gradient (s);
    s = u.y; i = 0;
    foreach_dimension()
      JJ[1][i++] = center_gradient (s);
    s = u.z; i = 0;
    foreach_dimension()
      JJ[2][i++] = center_gradient (s);
    double S2O2[dimension][dimension];
    for (int i = 0; i < dimension; i++)
      for (int j = 0; j < dimension; j++) {
	S2O2[i][j] = 0.;
	for (int k = 0; k < dimension; k++)
	  S2O2[i][j] += JJ[i][k]*JJ[k][j] + JJ[k][i]*JJ[j][k];
      }
    double lambda[dimension], ev[dimension][dimension];
    eigenvalues (S2O2, lambda, ev);
    l2[] = lambda[1]/2.;
  }
}
