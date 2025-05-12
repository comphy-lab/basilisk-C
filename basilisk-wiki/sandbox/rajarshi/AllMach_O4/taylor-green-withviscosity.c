/**
#Taylor Green Vortex - WENO5 + RK4 + O4 Projection + O4 Viscosity.
*/

#include "grid/multigrid.h"
#include "allmach_weno_poisson4_rk4.h"
#define u q

scalar ke[];

int main() {

  origin(-0.5,-0.5);
  foreach_dimension()
     periodic(right);
  const face vector muc[] = {0.0001,0.0001};
  mu = muc;
  for (N=32;N<=256;N*=2)
     run();

}

event init (i = 0) {
 
  foreach() {
    u.x[] =   (sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(cos(2.*pi*(y+Delta/2.))-cos(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));
    u.y[] = - (cos(2.*pi*(x+Delta/2.))-cos(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));
    p[]   = - (sin(4.*pi*(x+Delta/2.))-sin(4.*pi*(x-Delta/2.)) + sin(4.*pi*(y+Delta/2.))-sin(4.*pi*(y-Delta/2.)))/(16.*pi*Delta);
  }

  boundary ({p});
  boundary ((scalar *){u});
  foreach()
    foreach_dimension()
      g.x[] = - (p[-2] -8.*p[-1] + 8.*p[1] - p[2])/(48.*Delta);
  boundary ((scalar *){g});
}

event logfile (i++) {

   scalar div[];
   scalar divface[];
   foreach() {
     div[] = (u.x[-2,0] -8.*u.x[-1,0] + 8.*u.x[1,0] - u.x[2,0])/(48.*Delta) + (u.y[0,-2] -8.*u.y[0,-1] + 8.*u.y[0,1] - u.y[0,2])/(48.*Delta);
     ke[] = sq(u.x[])+sq(u.y[]);
     divface[] = 0.;
     foreach_dimension()
        divface[] = (uf.x[1] - uf.x[])/Delta;
   }
   printf("%d %d %g %g %.12g %g \n",N,i,t,normf(div).max,statsf(ke).sum,statsf(divface).sum);

}

event error (t=2) {

   scalar e[];
   foreach(){
     double u0 =    exp(-8.*pi*pi*mu.x[]*t)*(sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(cos(2.*pi*(y+Delta/2.))-cos(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));
     double v0 = -  exp(-8.*pi*pi*mu.x[]*t)*(cos(2.*pi*(x+Delta/2.))-cos(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/(sq(2.*pi*Delta));  
     e[] = norm(u) - sqrt(sq(u0) + sq(v0));
   }
   norm n = normf(e);
   fprintf(stderr,"%d %g %g %g\n",N,n.avg,n.rms,n.max);
   

}

event animation (i++,last) {
 if(N==256){
    output_ppm(p,min=-0.5,max=0.5,file="Pressure.mp4");
    output_ppm(ke,min=0,max=1,file="Kinetic-Energy.mp4");
 }
}

/**

#We produce animations of the vorticity, pressure and kinetic-energy fields 

![Animation of the pressure field.](taylorgreen_weno_poisson4_rk4/Pressure.mp4)

![Animation of the kinetic energy field.](taylorgreen_weno_poisson4_rk4/Kinetic-Energy.mp4) 


For this particular case, the Weno advection scheme
converges at fifth-order.

~~~gnuplot Accuracy of the solution as a function of the level of refinement
set xlabel 'Spatial resolution'
set ylabel 'Error norms'
set cbrange [1:1]
set logscale
set xtics 16,2,256
ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit [4:] f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit [4:] fm(x) 'log' u (log($1)):(log($4)) via am,bm
set xrange [16:512]
plot exp (f2(log(x))) t ftitle(a2,b2), \
     exp (fm(log(x))) t ftitle(am,bm),  \
     'log' u 1:3 t '|e|_2' ps 1.5, \
     'log' u 1:4 t '|e|_{max}' ps 1.5 lc 0 
~~~

The divergence of the centered velocity field is well-behaved.

~~~gnuplot Evolution of the maximum divergence of the centered velocity field
reset
set xlabel 'Time'
set ylabel 'Maximum divergence'
set cbrange [1:1]
set logscale y
set xrange [0:2]
set yrange [1e-9:0.1]
plot '< grep "^32 " out' u 3:4 w l t '32^2', \
     '< grep "^64 " out' u 3:4 w l t '64^2', \
     '< grep "^128 " out' u 3:4 w l t '128^2', \
     '< grep "^256 " out' u 3:4 w l t '256^2'   
~~~

By fitting the decrease of the kinetic energy, we get an estimate of
the numerical viscosity.

~~~gnuplot Equivalent Reynolds number as a function of resolution
reset
set xlabel 'Resolution'
set ylabel 'Equivalent Reynolds number'
set cbrange [1:1]
set logscale
set xtics 16,2,256
set xrange [16:512]
f(x)=a*exp(-b*x)
Re(b)=1./(b/(4.*(2.*pi)**2))

set print "Re"
fit [1:] f(x) '< grep "^32 " out' u 3:5 via a,b
print 32,Re(b)
fit [1:] f(x) '< grep "^64 " out' u 3:5 via a,b
print 64,Re(b)
fit [1:] f(x) '< grep "^128 " out' u 3:5 via a,b
print 128,Re(b)
fit [1:] f(x) '< grep "^256 " out' u 3:5 via a,b
print 256,Re(b)

plot 'Re' w lp t 'O2'

~~~
*/