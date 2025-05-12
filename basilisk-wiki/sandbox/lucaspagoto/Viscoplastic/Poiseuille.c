/**
# Poiseuille flow of a yield-stress fluid
*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

double tau_y = 0.2;  // yield stress
double mu_p = 1.0;   // plastic viscosity
double epsilon;  // regularization parameters
scalar shear_norm[], eta[];
face vector muv[];
tensor shear[];

int main() 
{
  periodic (top);
  size (0.5);
  
  stokes = true;
  TOLERANCE = 1e-5;
  
  u.t[right] = dirichlet(0);

  const face vector g[] = {0.,1.};
  a = g;
  mu = muv;

  for (epsilon = 1e-2; epsilon >= 1e-4; epsilon *= 1e-1)
  {
    for (N = 8; N <= 64; N *= 2)
    {
      run();
    }
  }
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
  double h = 0.5;
  double dpdy = 1.;
  double X = tau_y/dpdy;
  double U = (1/mu_p)*((-dpdy/2)*(sq(X) - sq(h)) + tau_y*(X - h));
  return (x > X ? (1/mu_p)*((-dpdy/2)*(sq(x) - sq(h)) + tau_y*(x - h)) : U);
}


scalar un[];

event logfile (t += 0.001; i <= 10000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
  //  fprintf (stderr, "%g %d %g %d %d\n", t, i, du, mgp.i, mgu.i);
  printf ("i = %d t = %.3g N = %d ep = %.e du = %.e\n", i,t,N, epsilon, du);
  fflush(stdout);
}


event profile (t = end) {
  char uprof[50];
  sprintf (uprof, "uprofile-%d-%.0e.txt", N, epsilon);
  FILE * up = fopen (uprof, "w");
  foreach()
    fprintf (up, "%g %g %g %g %g %g %g\n", x, y, u.x[], u.y[], p[], x - Delta/2., (u.y[] - u.y[-1])/Delta);
  printf ("\n");
  fclose (up);

  scalar e[];
  foreach(){
    e[] = u.y[] - analytical(x);
  }
  norm n = normf (e);
  fprintf (ferr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}





/**
~~~gnuplot Velocity profile: comparison between numerical and analytical solutions
reset
set title "N = 64"
set xlabel 'x'
set ylabel 'u.y'
h = 0.5
dpdy = 1.
tau_y = 0.2
mu = 1.0
X = tau_y/dpdy
U = (1/mu)*((-dpdy/2)*(X**2 - h**2) + tau_y*(X - h))
u(x) = x > X ? (1/mu)*((-dpdy/2)*(x**2 - h**2) + tau_y*(x - h)) : U
set xrange [0:0.5]
set yrange [0:0.06]
plot 'uprofile-64-1e-02.txt' u 1:4 w p t '1e-2', 'uprofile-64-1e-03.txt' u 1:4 w p t '1e-3', 'uprofile-64-1e-04.txt' u 1:4 w p lc 7 t '1e-4', u(x) w l lw 2 lc -1 t 'analytical'
~~~




~~~gnuplot Convergence of the error on u.y with $\epsilon = 1e-2$
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
fit f2(x) 'log' every ::0::3 u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::0::3 u (log($1)):(log($4)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::0::3 u 1:4 t '|h|_{max}' ps 1.5, 'log' every ::0::3 u 1:3 t '|h|_2' ps 1.5
~~~


~~~gnuplot Convergence of the error on u.y with $\epsilon = 1e-3$
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
fit f2(x) 'log' every ::4::7 u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::4::7 u (log($1)):(log($4)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::4::7 u 1:4 t '|h|_{max}' ps 1.5, 'log' every ::4::7u 1:3 t '|h|_2' ps 1.5
~~~


~~~gnuplot Convergence of the error on u.y with $\epsilon = 1e-4$
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
fit f2(x) 'log' every ::8::11 u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' every ::8::11 u (log($1)):(log($4)) via am,bm
plot exp (fm(log(x))) t ftitle(am,bm), exp (f2(log(x))) t ftitle(a2,b2), 'log' every ::8::11 u 1:4 t '|h|_{max}' ps 1.5, 'log' every ::8::11 u 1:3 t '|h|_2' ps 1.5
~~~


*/