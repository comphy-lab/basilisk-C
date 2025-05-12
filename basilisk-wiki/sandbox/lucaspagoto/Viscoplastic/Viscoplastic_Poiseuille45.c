/**
# Poiseuille flow in a periodic channel inclined at 45 degrees

We solve the flow of a Bingham material between parallel plates dirve by gravity using embedded boundaries and periodic boundary conditions.

## Code
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

int main()
{
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  
  stokes = true;
  TOLERANCE = 1e-4;
  init_grid (1 << 7);
  
  run();
  run();
}

double tEnd = 20.01;
int cc = 1;
#define WIDTH 0.5
#define EPS 1e-14

face vector alphav[], muv[];
vector un[];
scalar rhov[];
scalar eta[], shear[], yielded[], shear_norm[];
tensor shear1[];

event init (t = 0) {
  const face vector g[] = {1.,1.};
  a = g;
  mu = muv;

  vertex scalar phi[];
  foreach_vertex()
    phi[] = difference (union (y - x - EPS, x - y - 0.3 + EPS),
			y - x - 0.7 + EPS);
  boundary ({phi});
  fractions (phi, cs, fs);

  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);

  for (scalar s in {u})
    s.third = true;
}



event properties(i++){
  double epsilon = 1e-4;
  double tau_y = 0.05;
  double mu_p = 0.02;

  foreach()
  {
    shear1.x.x[] = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear1.x.y[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear1.y.x[] = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear1.y.y[] = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);

    shear_norm[] = sqrt(0.5)*sqrt(sq(shear1.x.x[])+sq(shear1.x.y[])+sq(shear1.x.y[])+sq(shear1.y.y[]));
    eta[] = (tau_y/(shear_norm[] + epsilon) + mu_p);
  }
  boundary ({eta, shear, shear_norm});


  foreach_face(){
    if (cc == 1)
      muv.x[] = mu_p*fm.x[];  // Newtonian
    else
      muv.x[] = fm.x[]*((eta[] + eta[-1,0])/2.);  //Non-Newtonian
  }

  boundary((scalar  *){muv});

  foreach()
  {
      if (shear_norm[]*eta[] > tau_y)
      {
        yielded[] = 1;
      }
      else
      {
        yielded[] = 0; 
      }
  }
  boundary((scalar  *){yielded});
}




event logfile(i+=1)
{
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
}



event fields (t = 0.; t <= tEnd; t += 0.5)
 {
  char ux[80];
  sprintf (ux, "u.x-%d.png", cc);
  draw_vof ("cs", "fs");
  squares ("u.x", linear = true, spread = -1);
  save (ux);

  char psy[50];
  sprintf (psy, "yprofile-%.2f-%d.txt", t, cc);
  FILE * fpy = fopen (psy, "w");
  for (double x = -0.5 + 1./pow(2,9); x < 0.; x += 1./pow(2,9))
  {
    double y = -x;
    double uu = sqrt(sq(interpolate(u.x, x, y))+sq(interpolate(u.y, x, y)));
    fprintf(fpy, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n", x, y, interpolate(u.x, x, y), interpolate(u.y, x, y), uu, interpolate(shear_norm, x, y), interpolate(yielded, x, y));
  }
  fclose (fpy);

/*
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
  output_field ((scalar *){u, p, cs, yielded}, fld2, n = 500, linear = true,  box = {{-0.5,-0.5},{0.5,0.5}});
  fclose (fld2);*/
}


event end (t = tEnd) {cc = 2;}

/**
![Velocity field Newtonian](Viscoplastic_Poiseuille45/u.x-1.png)
![Velocity field Non-Newtonian](Viscoplastic_Poiseuille45/u.x-2.png)

~~~gnuplot Velocity profile Newtonian
set xlabel 'Width'
set ylabel 'Velocity'
set size square
h = 0.247487
mu_p = 0.02
G = 2**0.5
f(x) = -(1/mu_p)*(G/2)*(x**2 - h**2)

plot [0:0.25][0:2.5] 'yprofile-20.00-1.txt' u ((-0.247487+(($1)**2+($2)**2)**0.5)):(((($3)**2+($4)**2)**0.5)) w l lw 3 lc 7 t "Numerical", f(x) w l lw 4 dt 3 lc -1 t "Theoretical"
~~~


~~~gnuplot Velocity profile Non-Newtonian
set xlabel 'Width'
set ylabel 'Velocity'
set size square
h = 0.247487
mu_p = 0.02
G = 2**0.5
tau_y = 0.05
ap = tau_y/G
U = ((1/mu_p)*(-(G/2)*(ap**2-h**2)+tau_y*(ap-h)))
g(x) = x > ap ? ((1/mu_p)*(-(G/2)*(x**2-h**2)+tau_y*(x-h))) : U

plot [0.:0.25][0:2] 'yprofile-20.00-2.txt' u ((-0.247487+(($1)**2+($2)**2)**0.5)):(((($3)**2+($4)**2)**0.5)) w l lw 3 lc 7 t "Numerical", g(x) w l lw 4 dt 3 lc -1 t "Theoretical"
~~~

*/