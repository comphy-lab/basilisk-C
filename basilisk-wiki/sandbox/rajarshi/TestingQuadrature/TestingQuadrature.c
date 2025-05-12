/**
#Testing Weno schemes - with and without quadrature formulation.
*/

#include "grid/multigrid.h"
#define dimension 2
#define BGHOSTS 2

#include "utils.h"

foreach_dimension(){

static double weno5_left_x (Point point, scalar X, double gradL){
  return ( (1./30.)*(X[-1,0]-2.*Delta*gradL) - (13./60.)*X[-2,0] + (47./60.)*X[-1,0] + (9./20.)*X[0,0] - (1./20.)*X[1,0] );
}

#if Quad
static double * quadrature_2D_x (Point point, face vector WENO_Avg ){

  static double QuadValues[3];

  QuadValues[0] = ( (-9.-22.*sqrt(15.))*WENO_Avg.x[0,2] + (116.+164.*sqrt(15.))*WENO_Avg.x[0,1] + 2186.*WENO_Avg.x[] + (116.-164.*sqrt(15.))*WENO_Avg.x[0,-1] + (-9.+22.*sqrt(15.))*WENO_Avg.x[0,-2] )/2400.;

  QuadValues[1] = (9*WENO_Avg.x[0,2] - 116.*WENO_Avg.x[0,1] + 2134.*WENO_Avg.x[] - 116.*WENO_Avg.x[0,-1] + 9.*WENO_Avg.x[0,-2])/1920.;

  QuadValues[2] = ( (-9.+22.*sqrt(15.))*WENO_Avg.x[0,2] + (116.-164.*sqrt(15.))*WENO_Avg.x[0,1] + 2186.*WENO_Avg.x[] + (116.+164.*sqrt(15.))*WENO_Avg.x[0,-1] + (-9.-22.*sqrt(15.))*WENO_Avg.x[0,-2] )/2400.;

  return(QuadValues);

}
#endif

}


static void flux_analytical (face vector flux_an){
  foreach_face(x)
     flux_an.x[] = sin(2.*pi*x)*cos(2.*pi*x)*(cos(4.*pi*(y-Delta/2.)) - cos(4.*pi*(y+Delta/2.)))/(8.*pi*Delta);
  foreach_face(y)
     flux_an.y[] = sin(2.*pi*y)*cos(2.*pi*y)*(cos(4.*pi*(x-Delta/2.)) - cos(4.*pi*(x+Delta/2.)))/(8.*pi*Delta); 
  boundary((scalar *){flux_an});  
}


void convergence(int depth){

  L0=1;
  origin(-0.5,-0.5);
  init_grid(1<<depth);
  foreach_dimension()
     periodic(left);
  

  scalar f[];
  vector u[];
  
  foreach(){
    f[] = (cos(2.*pi*(x+Delta/2.))-cos(2.*pi*(x-Delta/2.)))*(cos(2.*pi*(y+Delta/2.))-cos(2.*pi*(y-Delta/2.)))/sq(2.*pi*Delta);
    u.x[] = (sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/sq(2.*pi*Delta);
    u.y[] = (sin(2.*pi*(x+Delta/2.))-sin(2.*pi*(x-Delta/2.)))*(sin(2.*pi*(y+Delta/2.))-sin(2.*pi*(y-Delta/2.)))/sq(2.*pi*Delta);
  }
  boundary({f});
  boundary((scalar *){u});
  
  /* Interpolation of the velocity field to compute line average at cell boundaries */

  tensor gradu[];
  foreach()
    foreach_dimension()
       gradu.x.x[] = (u.x[1]-u.x[-1])/(2.*Delta); 
  boundary ((scalar *) {gradu});
  
  face vector u_avg[];
  foreach_face()
       u_avg.x[] = weno5_left_x  (point,u.x,gradu.x.x[-2]);
  boundary ((scalar *){u_avg});
    

  /* Interpolation of the tracer fields to compute line averages at cell boundaries */
  
  vector gradf[];
  foreach()
    foreach_dimension()
       gradf.x[] = (f[1]-f[-1])/(2.*Delta); 
  boundary ((scalar *) {gradf});
  
  face vector f_avg[];
  foreach_face()
       f_avg.x[] = weno5_left_x  (point,f,gradf.x[-2]);
  boundary ((scalar *) {f_avg});


  /** Computing the flux */

  face vector flux[];

#if Quad

  double *temp;
  double fq[3],uq[3];

  foreach_face(){
    temp = quadrature_2D_x (point,f_avg);
    fq[0] = *(temp);
    fq[1] = *(temp+1);
    fq[2] = *(temp+2);
    temp = quadrature_2D_x (point,u_avg);
    uq[0] = *(temp);
    uq[1] = *(temp+1);
    uq[2] = *(temp+2); 
    flux.x[] = (5.*fq[0]*uq[0] + 8.*fq[1]*uq[1] + 5.*fq[2]*uq[2])/18.;
  }
  
#else
    
  foreach_face()
    flux.x[] = u_avg.x[]*f_avg.x[];

#endif

  boundary ((scalar *){flux});

  
  face vector flux_an[];
  flux_analytical(flux_an);

  face vector Error[];
  double max = 0;
  foreach_face (){
    Error.x[] = flux.x[] - flux_an.x[];
    if(max<=fabs(Error.x[]))
      max = fabs(Error.x[]);
  }
  boundary((scalar *){Error});

  FILE * fp = fopen("Error.dat","w");
  foreach_face()
     fprintf(fp,"%g %g %g \n",x,y,Error.x[]);
  fclose(fp);

  fp = fopen("ErrorvsGrid.dat","a");
  fprintf(fp,"%g %g \n",pow(2,depth),max);
  fclose (fp);

}

int main(){
  system ("rm -f ErrorvsGrid.dat");
  for(int depth=6; depth <=9; depth++)
      convergence(depth);
}


/**

~~~gnuplot Error
splot 'Error.dat' u 1:2:3 w p t 'Error'
~~~ 

~~~gnuplot Convergence with spatial resolution weno - periodic
ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'ErrorvsGrid.dat' u (log($1)):(log($2)) via a,b
f2(x)=a2+b2*x
fit f2(x) '../TestingQuadratureYes/ErrorvsGrid.dat' u (log($1)):(log($2)) via a2,b2
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [32:1024]
set xtics 32,2,1024
set format y "10^{%L}"
set cbrange [1:2]
plot 'ErrorvsGrid.dat' u 1:2 t 'No-Quadrature' ps 1.5 lc 0, exp(f(log(x))) t ftitle(a,b), \
     '../TestingQuadratureYes/ErrorvsGrid.dat' u 1:2 t 'Quadrature' ps 1.5 lc 2, exp(f2(log(x))) t ftitle(a2,b2)      
~~~
*/
