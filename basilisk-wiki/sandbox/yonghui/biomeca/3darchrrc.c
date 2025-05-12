/**
#real 3d arch with KV model
This is a rough attempt to merge all the solvers together, plz understande the code by reading those simple exemples.
[RRC](biomeca/2dworrc.c) \&
[Womersley](biomeca/2dwoflow.c) \&
[ReadSTL](biomeca/readstl.c) \&
[MassFlow](biomeca/qtest.c) \&

*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"

#define MIN_LEVEL 5
#define LEVEL 7
#define MAX_LEVEL 10
#define eq_r 1.
#define ADAPT 0 //1 = use


#define tmax   1*2.*M_PI
#define alpha_w  10.

#define  ORR 0.1
#define  ORC 0.1
#define  OCC 10

FILE * fp1;
FILE * fp2;
FILE * fp3;
double	Pout;
double	Pold;
double	DEBIT;
double	Qold;
double	PPPP;




/**
##Importing the geometry */
void fraction_from_stl (scalar f, FILE * fp)
{
  /**
We read the STL file and compute the bounding box of the model. */
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
    maxl = max.x - min.x;
    scalar d[];
  distance (d, p);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f);
}

/**
## volume fraction field declaire */
scalar f0[];
int main()
{
  init_grid (8);
  size (8.);
  origin (-1*L0/4 ,-2*L0/5, 0.);
  fp1 = fopen ("testpoint", "w");
  fp2 = fopen ("pressure.dat", "w");
  fp3 = fopen ("uuu", "w"); 
  run();
}

/**
##Boundary Conditions*/
u.n[back] = ( y <= 0)? dirichlet(10.*f0[]*cos(t)):neumann(0);
u.t[back] = ( y <= 0)? dirichlet(0):neumann(0);
u.r[back] = ( y <= 0)? dirichlet(0):neumann(0);

p[front] = neumann(0) ;
pf[front] = neumann(0) ;

/**
## Initial Event*/
event init (t = 0) {
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("aortic_arch.stl", "r");
    char name2[10];
    for (int ii=MIN_LEVEL; ii<=LEVEL; ii++){	
      fraction_from_stl (f0, fp);
      adapt_wavelet ({f0}, (double[]){0.00001}, ii);	
      //dump file
      sprintf (name2, "dump_%d",ii);  	
      dump(file = name2);
    }
    fraction_from_stl (f0, fp);
    view (fov = 19.9677, quat = {0.614306,0.390179,0.155159,0.66807}, 
          tx = -0.120804, ty = -0.340923, width = 1920, height = 1080);
    draw_vof ("f0");
    squares ("f0", min = -10, max =10);
    save ("init.png");
    double viscosity =  sq(eq_r) / sq(alpha_w);
    const face vector muc[] = {viscosity, viscosity, viscosity};
    mu = muc;
  }
  DEBIT = 0.;
  Pout =0.;
  Pold =0.;
  Qold =0.;
}


event bc (t <= tmax; i += 1) {
  //BC : lateral velocity = 0
 foreach()
    foreach_dimension()
    u.x[] = f0[]*u.x[];
    boundary ((scalar *){u , p});
 
  /**
  inlet surface position [-1:3],[-3:0.5] in xy plan, cut z=0.4 seems good, you can also choose a plan by yourself (Recommande try it first in bview3D). */
  DEBIT = 0.;
  double zz = 0.4;
  for (double xx = -1.; xx < 3.; xx += 0.05 ){
    for (double yy = -3.; yy < 0.5; yy += 0.05 ){
      DEBIT += interpolate(u.z, xx, yy, zz)*sq(0.05);
     }	
  }


  /**
  RRC model , in this case, it seems still has some problem*/
  double Reee = 0.5;
  double alpha = (Reee + ORR) / ORC + 1.;
  double beta =  1. - Reee/(Reee + ORR +ORC);
  double QQQ = beta * (1.-exp( -alpha*t/((Reee+ORR)*OCC)))*DEBIT;
  double term1 = ORC * OCC * Pold / dt;
  double term2 = (ORR + ORC) * DEBIT ;
  double term3 = ORR * ORC * OCC * (DEBIT - Qold) / dt;
  double deno = 1. + ORC * OCC / dt ;
  //out puressure calculate  //RRC model
 
  Pout =  (term1 + term2 + term3) / deno;
  //Pout = ORR * DEBIT;  
  //Pout=0.	;
  fprintf(fp1,"%g %g %g %g %g\n" , t, DEBIT, Qold, Pout, Pold);

  //position du test point
  double px = 3.5;
  double py = 0.8;
  double pz = 2.;
  fprintf(stderr, "%d %g %g %g \n", i, t, dt,interpolate(u.x , px, py, pz));

  p[back] = ( y >= 0)? dirichlet(Pout):neumann(0) ;
  pf[back] = ( y >= 0)? dirichlet(Pout):neumann(0) ;

  //we give the P Q to the Pold Qold
  Pold = Pout;
  Qold = DEBIT;
}
/**
## Output 
you can delete some commande, it woundn't change the simulation*/
event snapshot (i += 1000)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  scalar l2[];
  lambda2 (u,l2);
  dump (file = name);
}
 
event video (t += tmax/300.){
  scalar l2[];
  lambda2 (u,l2);
  /**
input and output view*/
  view (fov = 20, quat = {0,0,0,1}, tx = -0.19684, ty = 0.00981146,  bg = {0.3,0.4,0.6}, width = 640, height = 640);
  // draw_vof("f0"); 	
  box();
  squares ("u.z", alpha = 0.4);
  cells( alpha = 0.4);
  save ("input_uz.mp4");
  /**
main chanel view*/
  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("u.z",n = {1,-1,0},alpha =2.4 );
  squares("u.z",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_uz.mp4");

  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("l2",n = {1,-1,0},alpha =2.4 );
  squares("l2",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_l2.mp4");

  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("u.x",n = {1,-1,0},alpha =2.4 );
  squares("u.x",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_ux.mp4");

  /**
cross section view*/
  view (fov = 11.1239, quat = {0.658889,-0.274942,-0.26525,0.648006}, tx = -0.207846, ty = -0.306235, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //  draw_vof("f0");
  squares("u.x",n = {-1,-1,0},alpha =-2.4);
  cells(n = {-1,-1,0},alpha =-2.4);
  save ("section_ux.mp4");

  squares("l2",n = {-1,-1,0},alpha =-2.4);
  cells(n = {-1,-1,0},alpha =-2.4);
  save ("section_l2.mp4");

}

/**
Adaption */

#if ADAPT
event adapt (i++) {
  double uemax = 0.01;
  adapt_wavelet ({u}, (double[]){uemax,uemax,uemax}, MAX_LEVEL, 4);
}
#endif




/**
# Results
~~~gnuplot testpoint
plot 'testpoint' u 1:2 w l t'debit' ,\
'testpoint' u 1:3 w l t'Q' ,\
'testpoint' u 1:4 w l t'Pressure' ,\
'testpoint' u 1:5 w l t'Pold'
~~~
 
~~~gnuplot pressure at y = L0/2
plot 'pressure.dat' u 3:4 w l t'pressure' ,\
'pressure.dat' u 3:5 w l t'pressure out'
~~~
*/ 
