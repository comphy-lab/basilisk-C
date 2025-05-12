/**
#Couette flow of a Bingham fluid between rotating cylinders

We solve the couette flow of a Bingham material between two rotating cylinders. The inner cylinder is rotating, while the outer cylinder is still.

##Code
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

int main()
{
  origin (-0.5, -0.5);
  
//  stokes = true;
  TOLERANCE = 1e-4;
  init_grid (1 << 6);
    run();
}


double tEnd = 50.01;

#define WIDTH 0.5

face vector alphav[], muv[];
scalar un[];
scalar rhov[];
scalar eta[], shear[], yielded[], shear_norm[];
tensor shear1[];


event init (t = 0) {
    if (!restore (file = "restart"))
  {
  const face vector g[] = {1.,1.};
  a = g;
  mu = muv;

  vertex scalar phi[];
  foreach_vertex()
      phi[] = difference (sq(0.5) - sq(x) - sq(y), sq(0.25) - sq(x) - sq(y));
  boundary ({phi});
  fractions (phi, cs, fs);

  u.n[embed] = dirichlet (x*x + y*y > 0.14 ? 0. : - y);
  u.t[embed] = dirichlet (x*x + y*y > 0.14 ? 0. :   x);

  for (scalar s in {u})
    s.third = true;

    foreach()
    un[] = u.y[];
  }
}


event properties(i++){
  double epsilon = 1e-5;

  foreach()
  {
    shear1.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear1.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear1.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear1.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

    shear_norm[] = sqrt(0.5)*sqrt(sq(shear1.x.x[])+sq(shear1.x.y[])+sq(shear1.x.y[])+sq(shear1.y.y[]));
    eta[] = (0.1/(shear_norm[] + epsilon) + 0.01);   // tau_y = 0.1 and mu_p = 0.01
  }
  boundary ({eta, shear, shear_norm});


  foreach_face(){
    muv.x[] = fm.x[]*((eta[] + eta[-1,0])/2.);   // Non-Newtonian
//    muv.x[] = fm.x[];   // Newtonian
  }
  boundary((scalar  *){muv});

  foreach()
  {
      if (shear_norm[]*eta[] > 0.1)
      {
        yielded[] = 1.;
      }
      else
      {
        yielded[] = 0.; 
      }
  }
  boundary((scalar  *){yielded});
}


event snapshot (t = 0; t <= tEnd; t += 0.05)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}


event logfile(i+=1)
{
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}


event stop (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1;
}


event fields (t = 0.; t <= tEnd; t += 0.05)
//event fields (t = end)
 {
  scalar utheta[];
  foreach() 
  {
    double theta = atan2(y, x);
    if (cs[] > 0.)
    {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
    }
    else
      p[] = utheta[] = 0.;
  }

  char ux[80];
  sprintf (ux, "u.x-%.2f.png", t);
  draw_vof ("cs", "fs");
  squares ("utheta", linear = true, spread = -1);
  save (ux);

  char yi[80];
  sprintf (yi, "yield-%.2f.png", t);
  draw_vof ("cs", "fs");
  squares ("yielded", linear = true, spread = -1);
  save (yi);

  char rad[50];
  sprintf (rad, "rad-%.2f.txt", t);
  FILE * radd = fopen (rad, "w");
  for (double x = -0.5 + 1./pow(2,8); x <= -0.25 - 1./pow(2,8); x += 1./pow(2,8))
  {
    double theta = atan2(0, x), r = sqrt(x*x + 0*0);
    fprintf (radd, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", x, r, theta, interpolate(u.x, x, 0), interpolate(u.y, x, 0), interpolate(utheta, x, 0), interpolate(p, x, 0), interpolate(shear_norm, x, 0), interpolate(yielded, x, 0));
  }
  fclose (radd);

/*  char psy[50];
  sprintf (psy, "yprofile-%.2f.txt", t);
  FILE * fpy = fopen (psy, "w");
  for (double x = -0.5 ; x <= 0.; x += 1./pow(2,9))
  {
    double y = -x;
    fprintf(fpy, "%.5f %.5f %.5f %.5f %.5f %.5f\n", x, y, interpolate(u.x, x, y), interpolate(u.y, x, y), interpolate(shear_norm, x, y), interpolate(yielded, x, y));
  }
  fclose (fpy);


  char psx[50];
  sprintf (psx, "xprofile-%.2f.txt", t);
  FILE * fpx = fopen (psx, "w");
  for (double x = -0.5 + 1./pow(2,8); x < 0.5; x += 1./pow(2,8))
  {
    double y = x + 0.25;
    fprintf(fpx, "%.5f %.5f %.5f\n", x, y, interpolate(p, x, y));
  }
  fclose (fpx);

  char namei1[80];
  sprintf (namei1, "interface-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (cs, fp1);
  fclose (fp1);

  char field2[80];
  sprintf (field2, "fields-output-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  output_field ((scalar *){u, utheta, p, cs, yielded, shear_norm}, fld2, n = 500, linear = true,  box = {{-0.5,-0.5},{0.5,0.5}});
  fclose (fld2);*/
}


event end (t = tEnd) {}


/** 
## Results
Velocity Field
![Velocity field](rotating_couette/u.x-4.00.png)
Yielded Region
![Yielded(red)/Unyielded(blue) regions](rotating_couette/yield-4.00.png)

~~~gnuplot $\tau_y/\mu_p = 10$
set xlabel "Radius"
set ylabel "Velocity"
set xtics 0,0.05,1
set ytics 0,0.05,1
set format x "%.2f"
set format y "%.2f"
set size square

bingham(r,Rl)=(r > Rl ? 0. : r*sqrt(2.)*10./4.*((Rl/r)**2-2.*log(Rl/r)-1.))

plot [0.25:0.5][0:0.25] bingham(x,0.35) t "Theoretical", 'rad-4.00.txt' u 2:6 t "Numerical"
~~~

*/