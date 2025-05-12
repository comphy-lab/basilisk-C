/**
# Poiseuille flow of a yield-stress fluid using mask.

The channel width is 0.50 + ESP. ESP takes the values of 1e-2, 5e-3 and 2.5e-3.
*/

#include "embed.h"
//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"


double EPS;
double tau_y = 0.2;  // yield stress
double mu_p = 1.0;   // plastic viscosity
double epsilon;  // regularization parameters
scalar shear_norm[], eta[];
scalar un[];
scalar e[];
face vector muv[];
tensor shear[];

int main() 
{
  periodic (top);
  size (1.);
  
  stokes = true;
  TOLERANCE = 1e-5;

  u.n[right] = dirichlet(0);
  u.t[right] = dirichlet(0);
  
  epsilon = 1e-4;
  for (EPS = 1e-2; EPS >= 1e-3; EPS *= 5e-1)
  {
    for (N = 16; N <= 128; N *= 2)
    {
      run();
    }
  }
}



event init (t = 0) {
  const face vector g[] = {0.,1.};
  a = g;
  mu = muv;
 
  mask (x > 0.5 + EPS ? right : none);
}


event properties(i++)
{
  foreach()
  {
    shear.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

    shear_norm[] = sqrt(0.5)*sqrt(sq(shear.x.x[])+sq(shear.x.y[])+sq(shear.x.y[])+sq(shear.y.y[]));
    eta[] = (tau_y/(shear_norm[] + epsilon) + mu_p);
  }
  boundary ({eta, shear, shear_norm});

  foreach_face(){
    muv.x[] = fm.x[]*((eta[] + eta[-1,0])/2.);  // viscosity
  }
  boundary((scalar  *){muv});
}


double analytical (double x)
{ 
  double h = 0.5 + EPS;
  double dpdy = 1.;
  double X = tau_y/dpdy;
  double U = (1/mu_p)*((-dpdy/2)*(sq(X) - sq(h)) + tau_y*(X - h));
  return (x > X ? (1/mu_p)*((-dpdy/2)*(sq(x) - sq(h)) + tau_y*(x - h)) : U);
}


event logfile (t += 0.0001; i <= 10000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1;
  printf ("i = %d t = %.3g N = %d EPS = %.3e ep = %.2e du = %.2e\n", i,t,N, EPS, epsilon, du);
  fflush(stdout);
}


event profile (t = end) {
  char uprof[50];
  sprintf (uprof, "uprofile-%d-%.3e.txt", N, EPS);
  FILE * up = fopen (uprof, "w");
  foreach()
    fprintf (up, "%g %g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[], x - Delta/2., (u.y[] - u.y[-1])/Delta);
  printf ("\n");
  fclose (up);


  foreach(){
    e[] = u.y[] - analytical(x);
  }
  norm n = normf (e);
  fprintf (ferr, "%d %g %g %g %g %g\n", N, epsilon, EPS, n.avg, n.rms, n.max);
}





/**

## Results

The channel width is 0.5 + EPS.

~~~gnuplot Velocity profile: comparison between numerical and analytical solutions for EPS = 1.25e-3
reset
set xlabel 'x'
set ylabel 'u.y'
h = 0.50125
dpdy = 1.
tau_y = 0.2
mu = 1.0
X = tau_y/dpdy
U = (1/mu)*((-dpdy/2)*(X**2 - h**2) + tau_y*(X - h))
u(x) = x > X ? (1/mu)*((-dpdy/2)*(x**2 - h**2) + tau_y*(x - h)) : U
set xrange [0:0.51]
set yrange [0:0.06]
plot 'uprofile-16-1.250e-03.txt' u 1:4 w p t 'N = 16', 'uprofile-32-1.250e-03.txt' u 1:4 w p t 'N = 32', 'uprofile-64-1.250e-03.txt' u 1:4 w p t 'N = 64', 'uprofile-128-1.250e-03.txt' u 1:4 w p t 'N = 128', u(x) w l t 'analytical' lw 2
~~~

For EPS = 1e-2, the simulation did not converge for N = 64. Thus, we do not show the error for EPS = 1e-2 here

~~~gnuplot Convergence of the error on u.y with EPS = 5e-3
reset
set xlabel 'Resolution'
set ylabel 'Error norms'
set key bottom left
set logscale
set cbrange [1:2]
set xtics 8,2,128
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit f2(x) 'log' every ::8::11 u (log($1)):(log($5)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::8::11 u (log($1)):(log($6)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::8::11 u 1:6 t '|h|_{max}' ps 1.5, 'log' every ::8::11 u 1:5 t '|h|_2' ps 1.5
~~~


~~~gnuplot Convergence of the error on u.y with EPS = 2.5e-3
reset
set xlabel 'Resolution'
set ylabel 'Error norms'
set key bottom left
set logscale
set cbrange [1:2]
set xtics 8,2,128
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit f2(x) 'log' every ::12::15 u (log($1)):(log($5)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::12::15 u (log($1)):(log($6)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::12::15 u 1:6 t '|h|_{max}' ps 1.5, 'log' every ::12::15 u 1:5 t '|h|_2' ps 1.5
~~~



~~~gnuplot Convergence of the error on u.y with EPS = 1.25e-3
reset
set xlabel 'Resolution'
set ylabel 'Error norms'
set key bottom left
set logscale
set cbrange [1:2]
set xtics 8,2,128
set grid
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit f2(x) 'log' every ::16::19 u (log($1)):(log($5)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::16::19 u (log($1)):(log($6)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::16::19 u 1:6 t '|h|_{max}' ps 1.5, 'log' every ::16::19 u 1:5 t '|h|_2' ps 1.5
~~~

*/
