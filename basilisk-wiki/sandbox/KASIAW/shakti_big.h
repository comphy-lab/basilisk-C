/**
Create fields of gap height (f), current time (ct), and the pointers to them
*/

scalar f[];
scalar ct[];
scalar * tracers = {f,ct};

/**Impose atmospheric pressure at the outflow*/

p[left] = dirichlet(0);


/**Run on a 1km square domain, with initial conditions imposed on maximum resolution grid*/
int main()
{
  size(1.e3);
  init_grid (1 << MAXLEVEL);
  run();
}

/**Randomly seeded initial perturbation around 1mm gap height */
event init (i = 0) {
  foreach() {
    f[] = 0.9e-3*(1+(1e-1* noise()));
    }
}

/**Output current time in days, and statistics */
event logfile (i++)
{
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g\n", 
	   i, t/(24.*3600.), mgp.i, s.sum, s.min, s.max);
}

/**Make movies of gap height and pressure head evolution*/
event movies (t += 2.592e6)
{
  output_ppm (f, min = 0.0007, max = 0.001, file = "f.mp4");
  
  output_ppm (p, min = 0, max = 120, file = "p.mp4");

  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, min = 0, max = MAXLEVEL, linear = false, file = "level.mp4");
}

event output (t += 2.592e6; t <= 3.1536e7)
{
  static int nf = 0;
  printf ("file: field-%d\n", nf++);
  output_field ({p,f}, stdout, N);
}

/**Create coefficients for non-linear Poisson equation, as described in equation (11) of [https://gmd.copernicus.org/articles/11/2955/2018/#section2]()*/
event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    double Re = om*norm(u)/nu;
    beta.x[] = clamp(-ff*ff*ff/(1.+Re)*g/12./nu,-1000,-1.e-10);
  }
  foreach(){
    double uxx = (u.x[]+u.x[1,0])/2;
    double uyy = (u.y[]+u.y[0,1])/2;
    double q = sqrt(sq(uxx)+sq(uyy));
    double N = rhoi*g*hi + rhow*g*(dzb*x - p[]);
    double col = A*N*N*N*f[];
    double mdiff = ((m[1,0]*f[1,0]+m[]*f[])*(f[1,0]-f[])-((m[-1,0]*f[-1,0]+m[]*f[])*(f[]-f[-1,0])))/sq(Delta)/2 + ((m[0,1]*f[0,1]+m[]*f[])*(f[0,1]-f[])-((m[0,-1]*f[0,-1]+m[]*f[])*(f[]-f[0,-1])))/sq(Delta)/2;
    m[] = (G + 12.*nu*rhow*sq(q)/(f[]*f[]*f[])*(1.+om*q/nu))/L+mdiff;
    zeta[] = (1./rhow-1./rhoi)*m[] + col + ieb;
    opening[] = m[]/rhoi - col;
    ones[] = 1.;
    }
}

#if TREE
event adapt (i++) {
  double tolerance = 0.000005;
  adapt_wavelet ({f}, &tolerance, MAXLEVEL);
  event ("coefficients");
}
#endif
