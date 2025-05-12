/**
# Advection of a scalar field - Uniform grid simulation

We use either the second-order Bell--Colella--Glaz scheme or the
5th-order WENO scheme

The Tracer function is chosen as the following volume averaged compact function

$$ Tracer(x,y) = \frac{1}{\Delta^2} \int_{y_c - \frac{\Delta}{2}}^
                 { y_c + \frac{\Delta}{2} } \int_{x_c - \frac{\Delta}{2}}
                 ^{ x_c + \frac{\Delta}{2} } \left(1 - \left( \frac{x - x_0}{a}
                 \right)^2 \right)^7  \left(1 - \left( \frac{y - y_0}{a} 
                 \right)^2 \right)^7  dx dy $$

*/


static double PointValues (double x, double y, double offsetx, double offsety) {

  if(sq(x-offsetx)<=0.0225 && sq(y-offsety)<=0.0225) 
      return ( pow(1-sq((x-offsetx)/0.15), 7) * pow(1-sq((y-offsety)/0.15), 7) );

  else
      return 0;

}

static double VolumeAvg (double x, double y, double Delta, double offsetx,
		         double offsety){

  double sum=0.;
  double q = (Delta/2.)*sqrt(3./5.);
  
  sum += (5./18.)*( (5./18.)*PointValues(x-q,y-q,offsetx,offsety) +
		    (8./18.)*PointValues(x,y-q,offsetx,offsety) +
		    (5./18.)*PointValues(x+q,y-q,offsetx,offsety));
  
  sum += (8./18.)*( (5./18.)*PointValues(x-q,y,offsetx,offsety) +
		    (8./18.)*PointValues(x,y,offsetx,offsety) +
		    (5./18.)*PointValues(x+q,y,offsetx,offsety));

  sum += (5./18.)*( (5./18.)*PointValues(x-q,y+q,offsetx,offsety) +
		    (8./18.)*PointValues(x,y+q,offsetx,offsety) +
		    (5./18.)*PointValues(x+q,y+q,offsetx,offsety));

  return (sum);

}

#include "grid/multigrid.h"
#define dimension 2
#include "../Header_Files/advection_r.h"
#include "view.h"

void face_velocity (face vector u, double t)
{
  foreach_face(x)
    u.x[] =  ( 1.5/(pi*Delta)) * sin(2.*pi*t/5.) * cos(pi*x)
             *( cos(pi*(y+Delta/2.)) - cos(pi*(y-Delta/2.)) );

  foreach_face(y)
    u.y[] = -( 1.5/(pi*Delta)) * sin(2.*pi*t/5.) * cos(pi*y)
             *( cos(pi*(x+Delta/2.)) - cos(pi*(x-Delta/2.)) );
  
  boundary((scalar *){u});

}

#if !WENO
event velocity (i++) {
  face_velocity (u, t);
}
#endif

scalar f[];
scalar * tracers = {f};

double cmax;
clock_t start,end;

int main()
{
  system("rm -f ErrorvsGrid.dat");
  origin (-0.5, -0.5);
  DT = 0.1;
  CFL = 0.8;

  FILE * fp2 = fopen ("TimeComputing.dat","w");
  for (N=128;N<=512;N*=2){
      start=clock();
      run();
      end=clock();
      fprintf(fp2,"%g \n",(end-start)/((double) CLOCKS_PER_SEC));
    }
  fclose(fp2);
  system("paste ErrorvsGrid.dat TimeComputing.dat
          > ErrorvsTime.dat");
}

scalar Error[];

event init (i = 0)
{ 
   foreach()
      f[] = VolumeAvg (x, y, Delta, -0.2, -0.236338);
   boundary({f});
}


event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %ld %.12f %g %g\n", t,
	   grid->tn, s.sum, s.min, s.max);
}


event deformed (t = 2.5) {
  view (fov = 19.35, width = 400, height = 400);
  clear();
  squares ("f", linear = true);
  save ("f.png");
}

event field (t = 5) {

  scalar e[];
  foreach()
    e[] = f[] - Neumann(x,y,Delta,-0.2,-0.236338);
  norm n = normf (e);

  FILE * fp1 = fopen("ErrorvsGrid.dat","a");
  fprintf (fp1, "%g \t %g \t %g \t %g \t %g\n",cmax,
	   sqrt(grid->tn), n.avg, n.rms, n.max);
  fclose(fp1);

  view (fov = 19.35, width = 400, height = 400);
  clear();
  squares ("f", linear = true);
  save ("e.png");

}

/**

# BCG SCHEME RESULTS

![Tracer field at $t=2.5$](advection/f.png)

![Error field at $t=5$ (BCG)](advection/e.png)

~~~gnuplot Convergence with spatial resolution (BCG)
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'ErrorvsGrid.dat' u (log($2)):(log($5)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'ErrorvsGrid.dat' u (log($2)):(log($3)) via a2,b2
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [64:1024]
set xtics 64,2,1024
set grid ytics
set cbrange [1:2]
plot 'ErrorvsGrid.dat' u 2:5 t 'max', 'ErrorvsGrid.dat' u 2:3 t 'norm1', \
      exp(f(log(x))) t ftitle(a,b), exp(f2(log(x))) t ftitle(a2,b2)
~~~ 

# WENO SCHEME RESULTS

![Error field at $t=5$ for $N=128$ (WENO)](advection-weno/e.png)

![Convergence with spatial resolution (WENO)](advection-weno/_plot0.png)

# Computation time analysis 

~~~gnuplot Error vs Time 
set xlabel 'Computing Time'
set ylabel 'Maximum-Error'
set logscale
set xrange [2:2048]
set xtics 2,2,2048
set grid ytics
plot 'ErrorvsTime.dat' u 6:5 w lp t 'Time-BCG' ps 2, \
     '../advection-weno/ErrorvsTime.dat' u 6:5 w lp t 'Time-WENO' ps 2 
~~~ 

*/
