//#include "grid/octree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "my_contact.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "utils.h"
#include "navier-stokes/perfs.h"
//#include "display.h"
//#include "view.h"
//#include "reduced.h"

//Dimensionless quantities

#define Oh 2.185e-03  // MU1/sqrt(Rho1*SIGMA*Radius)  // Ohnesorge number for water/air
#define Rratio 1000.0  // Rho1/Rho2
#define MUratio 100.0  //MU1/MU2

//System properties (SI)
#define Rho1 998.0  // liquid
#define MU1 0.00102

#define Rho2 Rho1/Rratio  // gas
#define MU2 MU1/MUratio

#define Radius 3.e-03 // Radius of the hemispherical droplet (3mm) V=56.5µL // Req = 0.79*Radius

#define SIGMA (sq(MU1)/(Rho1*Radius*sq(Oh))) // 72.8 mN/m for water/air
#define A (bo*SIGMA/(Rho1*0.79*Radius*0.79*Radius)) // acceleration of the plate in the frame of reference
#define Tc sqrt(Rho1*0.5*Radius*Radius*Radius/SIGMA) // capillary time 13.6ms


//System characteristics (dimensionless)
//We rescale everything according to Radius=1, Rho2_nd=Rho2 and Tc constant

#define Radius_nd (Radius/Radius)
#define Rho2_nd Rho2
#define Rho1_nd (Rratio*Rho2_nd)

#define SIGMA_nd (Rho1_nd*0.5*pow(Radius_nd,3)/pow(Tc,2))

#define A_nd (bo*SIGMA_nd/(Rho1_nd*0.79*Radius_nd*0.79*Radius_nd))

#define MU1_nd (Oh*sqrt(Rho1_nd*SIGMA_nd*Radius_nd))
#define MU2_nd (MU1_nd/MUratio)

#define theta_d  60.*(pi/180.)
#define tend 5*Tc //time upto which simulation runs or movie is generated

#define radius_underground Radius_nd/pow(  (1 - cos(theta_d) -0.5*sin(theta_d)*sin(theta_d)*cos(theta_d) )    ,(1./3.))

vector h[];
h.t[left] = contact_angle (theta_d); // impose a constant contact angle
u.t[left] = dirichlet(0); // no slip (comment these lines for free slip BC)
u.n[left] = dirichlet(0);
u.n[right] = u.n[] > 0 ? neumann(0) : 0;

int LEVEL = 14;

double bo;

int main (int argc, char * argv[])
{
    if (argc > 1)
       bo = atof(argv[1]);
     else {
       fprintf(stderr, "Error: bo value not specified.\n");
       return 1;  // Exit the program with an error code
     }

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */
  size (80*Radius_nd);
  f.height = h;
  rho1 = Rho1_nd;
  mu1 = MU1_nd;
  rho2 = Rho2_nd;
  mu2 = MU2_nd;
  f.sigma = SIGMA_nd;
  N = 1 << 8;
  init_grid(N);
  run();
}


event init (t = 0)
{
  if (!restore (file = "restart")){
  refine (sq(x) + sq(y) < 1.5*radius_underground && level < (LEVEL+1));      
  fraction (f, - (sq(y) + sq(x+radius_underground*cos(theta_d)) - sq(radius_underground)));


  
      fprintf(ferr,"Bo \t%e\t%e\t%e\n",bo,Rho1*A*Radius*Radius/SIGMA,Rho1_nd*A_nd*Radius_nd*Radius_nd/SIGMA_nd);
  fprintf(ferr,"Oh \t%e\t%e\t%e\n",Oh,MU1/sqrt(Rho1*SIGMA*Radius),MU1_nd/sqrt(Rho1_nd*Radius_nd*SIGMA_nd));
  fprintf(ferr,"Tc \t%e\t%e\t%e\n",Tc, sqrt(Rho1*pow(Radius,3)/SIGMA),sqrt(Rho1_nd*pow(Radius_nd,3)/SIGMA_nd));
  fprintf(ferr,"Val \t%e\t%e\t%e\n",A_nd,SIGMA_nd,MU1_nd);
  fprintf(ferr,"Volume in microliters \t%e\n",(1e9)*(4/3)*(pi)*pow(Radius*0.79 , 3));
  fprintf(ferr,"Actual Volume in microliters \t %e",(1e9)*(4/3)*(pi)*pow(Radius , 3)*0.5);  
  fprintf(ferr,"\n---------------------------------------\n---------------------------------------LEVEL = %d \n Number of grd pt per rad = %g \n",LEVEL , 10/(pow(2,LEVEL)));
}
  fprintf(ferr,"Bo \t%e\t%e\t%e\n",bo,Rho1*A*Radius*Radius/SIGMA,Rho1_nd*A_nd*Radius_nd*Radius_nd/SIGMA_nd);
  fprintf(ferr,"Oh \t%e\t%e\t%e\n",Oh,MU1/sqrt(Rho1*SIGMA*Radius),MU1_nd/sqrt(Rho1_nd*Radius_nd*SIGMA_nd));
  fprintf(ferr,"Tc \t%e\t%e\t%e\n",Tc, sqrt(Rho1*pow(Radius,3)/SIGMA),sqrt(Rho1_nd*pow(Radius_nd,3)/SIGMA_nd));
  fprintf(ferr,"Val \t%e\t%e\t%e\n",A_nd,SIGMA_nd,MU1_nd);
  fprintf(ferr,"Volume in microliters \t%e\n",(1e9)*(4/3)*(pi)*pow(Radius*0.79 , 3));
  fprintf(ferr,"Actual Volume in microliters \t %e",(1e9)*(4/3)*(pi)*pow(Radius , 3)*0.5);

  fprintf(ferr,"\n---------------------------------------\n---------------------------------------LEVEL = %d \n Number of grd pt per rad = %g \n",LEVEL , (pow(2,LEVEL))/10);

    boundary({f});
	dump("dump-initial");	
		
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] += A_nd;      
}


event droplets (i++){
heights(f,h);
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);
  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
 char name1[80];
 sprintf (name1, "tag_vol-%g", bo);
  static FILE * fp = fopen (name1,"w");
  for (int j = 0; j < n; j++)
    fprintf (fp, "%d %g %g %d %g %g %g\n", i, t/Tc , bo,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);

}

event end (t = tend){
  printf ("i= %d t = %g\n", i, t);
}

event adapt (i++){
  adapt_wavelet ((scalar*){f,u}, (double[]){0., 0.001 ,0.001}, LEVEL);
}

event snapshots (t += 0.1*Tc ; t <=tend)
{
  char name[80];
  sprintf (name, "s-%g-%d", t/Tc,LEVEL);
  scalar omega[];
  vorticity (u, omega);
  dump (name);

}

