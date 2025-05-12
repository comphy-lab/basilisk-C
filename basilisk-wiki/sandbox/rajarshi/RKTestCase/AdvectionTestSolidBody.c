/**
#Rigid Body advection test case
*/

#include "grid/multigrid.h"
#define dimension 2
#define BGHOSTS 2

#include "../Header_File/O5_Flux.h"
#include "../Header_File/runge-kutta-new.h"
#include "utils.h"

#define omega 10.*pi

face vector u[];
void face_velocity(face vector u, double time){
  foreach_face(x)
     u.x[] = -omega*y;
  foreach_face(y)
     u.y[] = omega*x;     
  boundary((scalar *){u});
} 

static double FuncValues (double x, double y, double offset) {
  if(sq(x-offset)<=0.01 && sq(y-offset)<=0.01) 
    return (pow(1-sq((x-offset)/0.1),7)*pow(1-sq((y-offset)/0.1),7));
  else
    return 0;
}

static double Neumann (double x, double y, double Delta, double offset) {
  double sum=0.;
  double q = (Delta/2.)*sqrt(3./5.);
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y-q,offset) + (8./18.)*FuncValues(x,y-q,offset) + (5./18.)*FuncValues(x+q,y-q,offset) );
  sum += (8./18.)*( (5./18.)*FuncValues(x-q,y,offset) + (8./18.)*FuncValues(x,y,offset) + (5./18.)*FuncValues(x+q,y,offset) );
  sum += (5./18.)*( (5./18.)*FuncValues(x-q,y+q,offset) + (8./18.)*FuncValues(x,y+q,offset) + (5./18.)*FuncValues(x+q,y+q,offset) );
  return (sum);
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
    s[] = Neumann(x,y,Delta,0.25);
  

  boundary({s});

  face_velocity(u,t); 

  double t,dt,tfinal;
  t=0.;
  tfinal = 0.2;
  dt = 1e-4;

  int i = 0;
  while(t<tfinal-dt){
    runge_kutta ({s},t,dt,tracer_fluxes,4);
    t+=dt;
    i++;
    fprintf(stderr,"%g \t\t %g \t\t %g \n",pow(2,depth),t,statsf(s).max);
    if(depth==9)
       if(i%5==0)
          output_ppm(s,min=0,max=1,file="Tracer.mp4");
  }
 
  
  foreach()
     error[] = s[] - Neumann(x,y,Delta,0.25);

  fp = fopen("ErrorvsGrid.dat","a");
  fprintf(fp,"%g %g \n",pow(2,depth),normf(error).max);
  fclose(fp);
  
}

int main(){
  system ("rm -f ErrorvsGrid.dat");
  for(int depth = 7; depth <= 9; depth++)
     convergence(depth);
}


/**

We produce animations of the pressure and kinetic-energy fields 

![Animation of the Tracer
 field.](AdvectionTestSolidBody/Tracer.mp4)

~~~gnuplot Convergence with spatial resolution weno - periodic
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'ErrorvsGrid.dat' u (log($1)):(log($2)) via a,b
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [64:1024]
set xtics 64,2,1024
set format y "10^{%L}"
set cbrange [1:2]
plot 'ErrorvsGrid.dat' u 1:2 t 'Error' ps 1.5 lc 0, exp(f(log(x))) t ftitle(a,b) lc 1
~~~
*/
