/**
#Splitted Riemann solver functions
*/

void flux (const double * state, double * flux, double * eigenvalue);

void riemann (const double * right, const double * left, double * f, int len)
{
  double fr[len], fl[len], er[2], el[2];
  flux (right, fr, er);
  flux (left,  fl, el);
  double ap = max(er[1], el[1]); ap = max(ap, 0.);
  double am = min(er[0], el[0]); am = min(am, 0.);
  double a = max(ap, -am); 

  if (a > 0.) {
    for (int i = 0; i < len; i++)
      f[i] = (ap*fl[i] - am*fr[i] + ap*am*(right[i] - left[i]))/(ap - am);
  }
  else
    for (int i = 0; i < len; i++)
      f[i] = 0.;
}

foreach_dimension ()
void get_fluxes_x (scalar * scalars, scalar * slopes, scalar * lflux,  Point point)
{

  int len = list_len (scalars);

  double r[len], l[len]; // right/left Riemann states
  double f[len];         // fluxes for each conserved quantity
  double dx = Delta/2.;
  int i = 0;
  scalar s,g;
  for (s,g in scalars,slopes) {
    r[i] = s[] - dx*g[];
    l[i++] = s[-1] + dx*g[-1];
  }

  riemann (r, l, f, len);
  i = 0;
  for (scalar fs in lflux)
    fs[] = fm.x[]*f[i++];

}