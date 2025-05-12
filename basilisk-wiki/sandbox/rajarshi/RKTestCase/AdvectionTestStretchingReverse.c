/**
# Advection of a scalar field

Using either the second-order Bell--Colella--Glaz scheme or the
5th-order WENO scheme. */

static double FuncValues (double x, double y, double offsetx, double offsety) {
  if(sq(x-offsetx)<=0.01 && sq(y-offsety)<=0.01) 
    return (pow(1-sq((x-offsetx)/0.1),7)*pow(1-sq((y-offsety)/0.1),7));
  else
    return 0;
}

static double Neumann (double x, double y, double Delta, double offsetx, double offsety) {
  double sum=0.;
  double q = (Delta/2.)*sqrt(3./5.);
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y-q,offsetx,offsety) + (8./18.)*FuncValues(x,y-q,offsetx,offsety) + (5./18.)*FuncValues(x+q,y-q,offsetx,offsety) );
  sum += (8./18.)*( (5./18.)*FuncValues(x-q,y,offsetx,offsety) + (8./18.)*FuncValues(x,y,offsetx,offsety) + (5./18.)*FuncValues(x+q,y,offsetx,offsety) );
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y+q,offsetx,offsety) + (8./18.)*FuncValues(x,y+q,offsetx,offsety) + (5./18.)*FuncValues(x+q,y+q,offsetx,offsety) );
  return (sum);
}

#include "grid/multigrid.h"
#define dimension 2
#define BGHOSTS 2

#include "../Header_File/O5_Flux.h"
#include "../Header_File/runge-kutta-new.h"

#include "utils.h"

face vector u[];
void face_velocity (face vector u, double t)
{
  vertex scalar psi[];
  foreach_vertex()
    psi[] = - 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;
  trash ({u});
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[])/Delta;
  boundary ((scalar *){u}); 
}


void convergence(int depth){

  L0=1.;
  origin(-0.5,-0.5);
  init_grid(1<<depth);
  foreach_dimension()
      periodic(left);

  FILE *fp;
  scalar s[],error[];

  foreach()
    s[] = Neumann(x,y,Delta,-0.2,-0.236338);
  boundary({s});

  face_velocity(u,t); 

  double t,dt,tfinal;
  t=0.;
  tfinal = 5;
  dt = 5e-4;

  int i = 0;
  while(t<tfinal){
    runge_kutta ({s},t,dt,tracer_fluxes,4);
    t+=dt;
    i++;
    fprintf(stderr,"%g \t\t %g \t\t %g \n",pow(2,depth),t,statsf(s).max);
    if(depth==8)
       if(i%10==0)
          output_ppm(s,min=0,max=1,file="Tracer.mp4");
  }
 
  
  foreach()
     error[] = s[] - Neumann(x,y,Delta,-0.2,-0.236338);

  fp = fopen("ErrorvsGrid.dat","a");
  fprintf(fp,"%g %g \n",pow(2,depth),normf(error).max);
  fclose(fp);
  
}

int main(){
  system ("rm -f ErrorvsGrid.dat");
  for(int depth = 6; depth <= 8; depth++)
     convergence(depth);
}


/**

We produce animations of the Tracer field

![Animation of the Tracer
 field.](AdvectionTestStretchingReverse/Tracer.mp4)

~~~gnuplot Convergence with spatial resolution weno
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'ErrorvsGrid.dat' u (log($1)):(log($2)) via a,b
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [32:512]
set xtics 32,2,512
set format y "10^{%L}"
set cbrange [1:2]
plot 'ErrorvsGrid.dat' u 1:2 t 'Error' ps 1.5 lc 0, exp(f(log(x))) t ftitle(a,b) lc 1
~~~
*/
