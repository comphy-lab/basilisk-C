/**
# Taylor Green vortices solver !
*/

#include "grid/multigrid.h"

#if ALL_MACH
#include "all-mach.h"
# include "bcg.h"
# define u q
event tracer_advection (i++)
  advection ((scalar *){q}, uf, dt, (scalar *){g});

#else
#include "navier-stokes/centered.h"
#endif


int main(){

   origin (-0.5,-0.5);
   foreach_dimension()
      periodic(right);
   for(N=32;N<=64;N*=2)
      run();
}

event init(i=0) {

   foreach(){
      u.x[] = -cos(2.*pi*x)*sin(2.*pi*y);
      u.y[] =  sin(2.*pi*x)*cos(2.*pi*y);
      p[]   = -(cos(4.*pi*x) + cos(4.*pi*y))/4.;
      }
   boundary({p});
   foreach()
      foreach_dimension()
          g.x[] = -(p[-2] - 15.*p[-1] + 15.*p[] - p[1])/(12.*Delta);
   boundary((scalar *){g});

}

event logfile (i++){
  
  scalar div[],ke[];
  foreach(){
    div[] = (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
    ke[]  = sq(u.x[]) + sq(u.y[]);
  }
  printf("%d %d %g %g %g \n",N,i,t,normf(div).max,statsf(ke).sum);
}

event error (t=2){

   scalar e[];
   foreach(){
      double u0 = -cos(2.*pi*x)*sin(2.*pi*y);
      double v0 =  sin(2.*pi*x)*cos(2.*pi*y);
      e[] = norm(u) - sqrt(sq(u0)+sq(v0));
   }
   norm n = normf(e);
   fprintf(stderr, "%d %g %g %g \n",N,n.avg,n.rms,n.max);
}

/**
For this particular case, the Bell--Collela--Glaz advection scheme
converges at third-order.

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
     'log' u 1:4 t '|e|_{max}' ps 1.5 lc 0, \
     '../taylor-green-all-mach/log' u 1:3 t '|e|_2 (all Mach)' ps 1.5, \
     '' u 1:4 t '|e|_{max} (all Mach)' ps 1.5
~~~

The divergence of the centered velocity field is well-behaved.

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
     '< grep "^256 " out' u 3:4 w l t '256^2', \
     '< grep "^32 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '32^2 (all Mach)',		  \
     '< grep "^64 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '64^2 (all Mach)',		   \
     '< grep "^128 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '128^2 (all Mach)',		   \
     '< grep "^256 " ../taylor-green-all-mach/out' \
     u 3:4 w l t '256^2 (all Mach)'
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

set print "Re-all-mach"
fit [1:] f(x) '< grep "^32 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 32,Re(b)
fit [1:] f(x) '< grep "^64 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 64,Re(b)
fit [1:] f(x) '< grep "^128 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 128,Re(b)
fit [1:] f(x) '< grep "^256 " ../taylor-green-all-mach/out' u 3:5 via a,b
print 256,Re(b)

plot 'Re' w lp t 'centered', 'Re-all-mach' w lp t 'all Mach'
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/reynolds.html)
*/