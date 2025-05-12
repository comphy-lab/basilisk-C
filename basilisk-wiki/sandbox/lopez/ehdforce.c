/**
# Embed EHD stress test
*/
#define W 5.0
#define A 0.25

#define NX 1.
#define NY 1.

#include "embed.h"
#include "src/embed_utils.h"
#include "run.h"

scalar phi[], rho = unity;
face vector a[];
(const) face vector alpha = unityf, epsilon = unityf;

#include "src/stress.h"

static double exact (double x, double y)
{ 
  return A*sin(2.*NX*pi/W*x)*sin(2.*NY*pi/W*y);
}

#define FX(x,y) -4.*sq(A)*NX*(sq(NX)+sq(NY))*sq(sin(2.*NY*pi/W*y))*	\
			 sin(4.*NX*pi/W*x)*pow(pi/W,3.)
#define FY(x,y) -4.*sq(A)*NY*(sq(NX)+sq(NY))*sq(sin(2.*NX*pi/W*x))*	\
			 sin(4.*NY*pi/W*y)*pow(pi/W,3.)


phi[embed] = dirichlet(exact(x,y));
phi[left] = exact(x-Delta/2.,y);
phi[right] = exact(x+Delta/2.,y);
phi[top] = exact(x,y + Delta/2.);
phi[bottom] = exact(x,y - Delta/2.);


int main()
{
  size(W);
  origin(-W/2., -W/2.);
  for (int lvl = 5; lvl <= 8; lvl++) {
    N = 1 << lvl;
    run();
  }
}  

event defaults (i = 0) {
  vertex scalar phiv[];
  foreach_vertex()
    phiv[] = sq(x) + sq(y) - sq(1.1);
  
  fractions (phiv, cs, fs);

  epsilon = fs;
  alpha = fs;
  rho = cs;
  foreach() {
    double xc = x, yc = y;
#if 0    
      if (cs[] > 0. && cs[] < 1.) {
	coord n = facet_normal (point, cs, fs), p;
	double alpha = plane_alpha (cs[], n);
	line_center (n, alpha, cs[], &p);
	xc += p.x*Delta, yc += p.y*Delta;
      }
#endif      
    phi[] = cs[] > 0. ? exact(xc, yc) : nodata;
  }

  //  output_cells(stdout);
  //  output_facets(cs,stdout,fs);
}

event dibujar (i = 0) {
  face vector err[];

  foreach_face (x) {
    err.x[] = a.x[] - FX(x,y);
    printf ("x: %g %g %g %g %g  \n", x, y, a.x[], FX(x,y), err.x[]);
  }
     
  foreach_face (y) {
    err.y[] = a.y[] - FY(x,y);
    printf ("y: %g %g %g %g %g \n", x, y, a.y[], FY(x,y),err.y[]);
  }
  
  norm nx = normf (err.x), ny = normf (err.y);
     
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N,
	   nx.avg, nx.rms, nx.max,
	   ny.avg, ny.rms, ny.max);

}

/**
## Results

#### Errors

~~~gnuplot Average error convergence
reset
set terminal svg font ",16"
set key bottom left spacing 1.1
set xtics 8,2,512
set grid ytics
set ytics format "%.0e" 1.e-12,100,1.e2
set xlabel 'N'
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
fx(x)=a+b*x
fit fx(x) 'log' u (log($1)):(log($4)) via a,b
fx2(x)=a2+b2*x
fit fx2(x) 'log' u (log($1)):(log($2)) via a2,b2
fy(x)=c+d*x
fit fy(x) 'log' u (log($1)):(log($7)) via c,d
fy2(x)=c2+d2*x
fit fy2(x) 'log' u (log($1)):(log($5)) via c2,d2
fxO(x)=f+g*x
set xlabel 'N'
set ylabel 'Error'
set logscale
set yrange [:]
plot 'log' u 1:4 t 'Fx:max', exp(fx(log(x))) lw 2 t ftitle(a,b), \
     'log' u 1:2 t 'Fx:norm1', exp(fx2(log(x))) lw 2 t ftitle(a2,b2)
~~~

*/
