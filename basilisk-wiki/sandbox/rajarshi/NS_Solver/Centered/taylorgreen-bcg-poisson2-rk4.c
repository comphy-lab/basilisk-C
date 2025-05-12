/**
#Taylor green vortex - bcg advection + Order 2 viscosity and projection + RK4 time advection
*/
#include "grid/multigrid.h"
#include "centered-bcg-poisson2-rk4.h"

scalar ke[];

int main() {

  origin(-0.5,-0.5);
  foreach_dimension()
     periodic(right);
  for (N=32;N<=256;N*=2)
     run();

}

event init (i=0) {

   foreach() {
     u.x[] = - cos(2.*pi*x)*sin(2.*pi*y); 
     u.y[] =   sin(2.*pi*x)*cos(2.*pi*y);
     p[]   = -(cos(4.*pi*x) + cos(4.*pi*y))/4.;
   }
   boundary({p});
   foreach()
     foreach_dimension()
        g.x[] = - (p[1]-p[-1])/(2.*Delta);
   boundary((scalar *){g});

}

event logfile (i++) {

   scalar div[];
   foreach() {
     div[] = (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
     ke[] = (sq(u.x[])+sq(u.y[]))/2.;
   }
   printf("%d %d %g %g %g \n",N,i,t,normf(div).max,statsf(ke).sum);

}

event error (t=2) {

   scalar e[];
   foreach(){
     double u0 = - cos(2.*pi*x)*sin(2.*pi*y);
     double v0 =   sin(2.*pi*x)*cos(2.*pi*y);
     e[] = norm(u) - sqrt(sq(u0) + sq(v0));
   }
   norm n = normf(e);
   fprintf(stderr,"%d %g %g %g\n",N,n.avg,n.rms,n.max);

}

event animation (i++,last) {
    output_ppm(p,min=-0.5,max=0.5,file="Pressure.mp4");
    output_ppm(ke,min=0,max=0.5,file="Kinetic-Energy.mp4");
}

/**
#We produce animations of the vorticity, pressure and kinetic-energy fields 

![Animation of the pressure field.](taylorgreen-bcg-poisson2-rk4/Pressure.mp4)

![Animation of the kinetic energy field.](taylorgreen-bcg-poisson2-rk4/Kinetic-Energy.mp4) 


#For this particular case, the BCG scheme converges at order 3

~~~gnuplot Accuracy of the solution as a function of the level of refinement
set xlabel 'Spatial Resolution'
set ylabel 'Error norm'
set cbrange [1:1]
set logscale
set xtics 16,2,256
ftitle(a,b) = sprintf("order %4.2f",-b)
f2(x) = a2+b2*x
fit [4:] f2(x) 'log' u (log($1)):(log($3)) via a2,b2
set xrange [16:512]
plot exp (f2(log(x))) t ftitle(a2,b2), \
     'log' u 1:3 t '|e|_2' ps 1.5, \
     'log' u 1:4 t '|e|_{max}' ps 1.5 lc 0
~~~ 

#The divergence of the centered velocity field is well behaved.

~~~gnuplot Evolution of the maximum divergence of the centered velocity field
reset
set xlabel 'Time'
set ylabel 'Maximum divergence'
set cbrange [1:1]
set logscale y
set xrange [0:2]
set yrange [1e-4:]
plot '< grep "^32 " out' u 3:4 w l t '32^2', \
     '< grep "^64 " out' u 3:4 w l t '64^2', \
     '< grep "^128 " out' u 3:4 w l t '128^2', \
     '< grep "^256 " out' u 3:4 w l t '256^2'
~~~

#By fitting the decrease of kinetic energy we get an estimate of the numerical viscosity

~~~gnuplot Equivalent Reynolds number as a function of resolution
reset
set xlabel 'Resolution'
set ylabel 'Equivalent Reynolds number'
set cbrange[1:1]
set logscale
set xtics 16,2,256
set xrange [16:512]
f(x) = a*exp(-b*x)
Re(b) = 1./(b/(4.*(2.*pi)**2))

set print "Re"
fit [1:] f(x) '< grep "^32 " out' u 3:5 via a,b
print 32,Re(b)
fit [1:] f(x) '< grep "^64 " out' u 3:5 via a,b
print 64,Re(b)
fit [1:] f(x) '< grep "^128 " out' u 3:5 via a,b
print 128,Re(b)
fit [1:] f(x) '< grep "^256 " out' u 3:5 via a,b
print 256,Re(b)

plot 'Re' w lp t 'centered'
~~~
*/
