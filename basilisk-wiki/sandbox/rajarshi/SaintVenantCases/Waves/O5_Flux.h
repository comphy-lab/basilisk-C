/**
# O5 FLUX COMPUTATION HEADER FILE
*/

#define BGHOSTS 2

foreach_dimension() {

  static double weno5_left_x (Point point, scalar s, double gradL) {
    return (1./3.*(s[-1] - 2.*Delta*gradL) - 13./6.*s[-2] +
	    47./6.*s[-1] + 9./2.*s[0] - 1./2.*s[1])/10.;
  }

  static double weno5_right_x (Point point, scalar s, double gradR) {
    return (-1./2.*s[-2] + 27./6.*s[-1] + 47./6.*s[0] -
	    13./6.*s[1] + 1./3.*(s[0] + 2.*Delta*gradR))/10.;
  }

  #define SQRT15 3.87298334620742

#if dimension > 1  
  static void quadrature3p_2D_y (Point point, face vector avg, double * q) {
    q[0] = ((- 9. - 22.*SQRT15)*avg.y[2] + (116. + 164.*SQRT15)*avg.y[1] +
	    2186.*avg.y[] +
	    (- 9. + 22.*SQRT15)*avg.y[-2] + (116. - 164.*SQRT15)*avg.y[-1])
      /2400.;

    q[1] = (9.*avg.y[2] - 116.*avg.y[1] +
	    2134.*avg.y[] +
	    9.*avg.y[-2] - 116.*avg.y[-1])/1920.;

    q[2] = ((- 9. + 22.*SQRT15)*avg.y[2] + (116. - 164.*SQRT15)*avg.y[1] +
	    2186.*avg.y[] +
	    (- 9. - 22.*SQRT15)*avg.y[-2] + (116. + 164.*SQRT15)*avg.y[-1])
      /2400.;
  }
#endif

}

static void fluxes (scalar f, face vector u, face vector flux)
{
  vector grad[];
  foreach()
    foreach_dimension()
    grad.x[] = (f[1] - f[-1])/(2.*Delta);
  boundary ((scalar *){grad});

  foreach_face() {
    if (u.x[] >= 0)
      flux.x[] = weno5_left_x (point, f, grad.x[-2]);
    else
      flux.x[] = weno5_right_x (point, f, grad.x[1]);
  }
}


void Tracer_fluxes (scalar f, face vector u, face vector flux) {

#if dimension == 1 
    fluxes (f, u, flux);
    foreach_face()
      flux.x[] *= u.x[];

#elif dimension == 2
    face vector line_avg[];
    fluxes (f, u, line_avg);
    boundary ((scalar *) {line_avg});
    
    foreach_face() {
      double fq[3];
      quadrature3p_2D_x (point, line_avg, fq);
      double uq[3];
      quadrature3p_2D_x (point, u, uq);
      flux.x[] = (5.*fq[0]*uq[0] + 8.*fq[1]*uq[1] + 5.*fq[2]*uq[2])/18.;
    }
#endif

    boundary_flux ({flux});

}

struct Advection {
  scalar * tracers;
  face vector u;
  double dt;
};

void advection (struct Advection p){
 
  for (scalar f in p.tracers) {
      face vector flux[];
      Tracer_fluxes(f,p.u,flux);
      foreach()
         foreach_dimension()
            f[] += p.dt*(flux.x[]-flux.x[1])/Delta;
   }
   boundary(p.tracers);
}
