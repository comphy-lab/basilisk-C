#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta)))	\
		      + ((neumann (0))*					\
			 ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))

double  heaviside ( double x, double minx1, double maxx1, double minx2, double maxx2, double minx3, double maxx3 ) {
  double a = 0;
  if ((x > minx1 && x < maxx1) || (x > minx2 && x < maxx2) || (x > minx3 && x < maxx3))
    a = 1;
  return a;
}

double T = 12*24;
double Pi1 = -6, Pi2 = 500, Pi4 = 2500, Re = 500;
double LAMBDA,LAMBDA_s, LAMBDA_w, sqN, NU, scaling_e, b0, L, Nbv, w, d;
double x1, x2, x3, x4, L0, QN, QN1, QN2, b_s, b_w;

int maxlevel = 9, minlevel = 7;

scalar b[], * tracers = {b};
b[bottom] = (robin(LAMBDA_s, NU, QN) * heaviside(x,X0,x1,x2,x3,x4,L0)) + (robin(LAMBDA_w, NU, QN) * heaviside(x,x1,x2,x3,x4,0,0));
b[top] = dirichlet (sqN*y);

FILE * gp;
face vector av[];

int main() {
  Nbv     = 1;
  sqN     = sq(Nbv);
  b0      = 1;
  L0      = 10;
  L       = b0/sqN;
  QN      = 1.2*10^-2; 
  LAMBDA  = 6*10^-3; 
  LAMBDA_s = LAMBDA;
  LAMBDA_w = LAMBDA/10;
  w       = 0.5;
  d       = 2;
  NU      = b0/(sq(Nbv)*Nbv)/Re;
  X0      = -L0/2;
  b_s     = QN/LAMBDA_s;
  b_w     = QN/LAMBDA_w;
  
  x1      = -(d/2) - (w/2);
  x2      = x1 + w;
  x3      = x1 + d;
  x4      = x3 + w;;

  scaling_e = sqN/sq(b_w)/w/L; 

  
  printf("sqN = %f\n",sqN);
  printf("QN = %f\n", QN);
  printf("Lambda = %f\n", LAMBDA);
  printf("b0 = %f\n", b0);
  printf("L0 = %f\n", L0);
  printf("NU = %f\n", NU);
  printf("b_s = %f\n", b_s);
  printf("b_w = %f\n", b_w);
  periodic (left);
  u.t[bottom] = dirichlet (0.);
#if (dimension == 3)
  u.r[bottom] = dirichlet (0.);
  periodic (back);
#endif
  const face vector muc[] = {NU, NU, NU};
  mu = muc;
  a = av;

  p.prolongation = p.refine = refine_linear; //3rd-order interpolation (vertical)
  foreach_dimension()
    u.x.refine = refine_linear;
  N = 1 << minlevel; // = bitwise 2^minlevel
  
  run();
}

event init (t = 0) {
  DT = 1./(sqrt(sqN)*10.);
  foreach() {
    b[] = sqN*y;
    foreach_dimension()
      u.x[] = noise()/1000.;
  }
  // Set initial buoyancy to surface:
  // b[bottom] = dirichlet (b0*heaviside(x,x1,x2,x3,x4,0,0));
  boundary ({b, u});
  gp = popen ("gnuplot", "w");
  fprintf(gp,
	  "set term pngcairo size 500, 500\n"
	  "set xr [%g: %g]\n"
	  "set yr [0 : %g]\n"
	  "set grid\n"
	  "set size square"
	  "set key topleft"
	  "set xlabel 'b'\n"
	  "set ylabel 'y/L_0'\n",
	  b0/LAMBDA, sqN*L0, L0);
}

event acceleration (i++) {
  coord grav = {0, 1, 0};
  boundary({b});
  foreach_face()
    av.x[] = grav.x*face_value(b, 0);
}

event tracer_diffusion(i++) 
  diffusion (b, dt, mu);

event adapt (i++) {
  double be = b0/150., ue = 0.01; 
  adapt_wavelet ({b, u}, (double[]){be, ue, ue, ue}, maxlevel, minlevel);
}

event damp (i++) {
  foreach() {
    if (y > 2.*L0/3.) {
      foreach_dimension()
	u.x[] *= exp (-dt*(y - 2*L0/3));
      b[] -= (b[] - sqN*y)*(1. - exp (-dt*(y - 2*L0/3)));
    }
  }
  boundary ({b, u});
}

event time_series (i += 25) {
  double e = 0, diss = 0, sbflx = 0, e_dim = 0;
  foreach(reduction(+:diss) reduction(+:e) reduction(+:sbflx)) {
    foreach_dimension() {
      e += sq(u.x[])*dv();
      diss += dv()*(sq(u.x[1] - u.x[-1]) +
		    sq(u.x[0,1] - u.x[0,-1]) +
		    sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
    if (Delta > y)
      sbflx += (b[0,-1] - b[])*dv()/sq(Delta);
  }
  diss *= -NU;
  sbflx *= NU;
  e /= 2.;
  e_dim += e*scaling_e;
  
  static FILE * fp = fopen ("timeseries", "w");
  if (i == 0)
    fprintf (fp, "t\ti\tn\twct\tspeed\te\te_dim\tdiss\tsbflx\n");
  fprintf (fp, "%g\t%d\t%ld\t%g\t%g\t%g\t%g\t%g\t%g\n",
	   t, i, grid->tn, perf.t, perf.speed, e, e_dim, diss, sbflx);
  fflush (fp);
}

#include "profile5c.h"
int frame = 0;
event movies (t += T/1000.) {
  scalar db[];
  foreach() {
    db[] = 0;
    foreach_dimension() 
      db[] += sq((b[1] - b[-1])/(2*Delta));
    db[] = db[] > 0 ? log(sqrt(db[]) + 1.) : 0;
  }
  boundary ({db});
  output_ppm (b, file = "b.mp4", linear = true,
	      n = 500, box = {{X0,0},{(X0+L0),L0}}, min = 0, max = (sqN*L0/4));
  output_ppm (db, file = "db.mp4", linear = true,
	      n = 500, box = {{X0,0},{(X0+L0),L0}}, min = 0, max = 3*log(sqN + 1));
  double by, yp = 1e-6;
  fprintf (gp,
	   "set output 'plot%d.png'\n"
	   "set title 'Time since sunrise:  %02d:%02d (HH:MM)'\n"	
	   "plot '-' w l lw 3 t 'buoyancy'\n",
	   frame, (int)(t/60.), (int)floor(fmod(t, 60.)));
  while (yp < L0) {
    by = 0;
    double dy = average_over_yp ({b}, &by, yp);
    fprintf (gp,"%g %g\n", by, yp);
    yp += dy;
  }
  fprintf (gp,"e\n");
  frame++;
}

event stop (t = T) {
  pclose (gp);
  system ("rm mov.mp4");
  system ("ffmpeg -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  sleep (10);
  system ("ffmpeg -y -i b.mp4 -i mov.mp4 -filter_complex hstack output.mp4");
  return 1;
}
