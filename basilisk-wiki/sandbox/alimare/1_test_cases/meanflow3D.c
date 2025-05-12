/**
# Evolution by mean curvature flow

Solve
$$
\mathbf{n}\cdot\partial_\tau\mathbf{X} = -\kappa
$$
See e.g [Deckelnick et al,
2005](https://homepages.warwick.ac.uk/~masgat/PAPERS/DecDziEll05_ActaNumerica.pdf).

![Initial Interface is a deformed sphere](meanflow3D/cs.png)

![Evolution of 0-level-set](meanflow3D/cs.mp4)
*/
#define QUADRATIC 1
#define BGHOSTS 2
#include "profiling.h"
#include "grid/octree.h"
#include "embed.h"
#include "advection.h"
#include "curvature.h"

scalar * tracers = NULL, * interfaces = NULL;
#include "../simple_discretization.h"
#include "../LS_reinit.h"
#include "../LS_recons.h"
#include "../LS_curvature.h"
#include "view.h"


double circle (double x, double y, double z, double R0, double A, double m)
{ 
  double eps = 1.e-15;
  double r= sqrt(sq(x) + sq(y)), R = sqrt(sq(x) + sq(y) + sq(z) );
  double theta = atan2(y,x);


  double psi = acos(z/(R+eps));
  double modifiedR = R*(1.+ A*R0*(cos(m*psi) +  //x
                                  cos(m*psi) +  //y
                                  sin(m*psi)) );                     //z
  return R0 - modifiedR;
}

int main()
{
  origin(-0.5,-0.5,-0.5);
  int MAXLEVEL = 6;
  // this is the diffusive timestep limit (with diffusion coeff unity)
  double DT2 = 0.5*sq(L0/(1<<MAXLEVEL));
  fprintf(stderr,"DT2 %g\n",DT2);
  init_grid(1<< (MAXLEVEL-2));

  vertex scalar distn[];

  double R0  = 0.25;
  double A = 0.2;
  double m = 8;
  double NB_width = 0.08;

  foreach_vertex() {
      distn[] = clamp(circle(x, y, z, R0, A, m), -NB_width, NB_width);
  }
  boundary({distn});
  fractions(distn,cs,fs);
  
  scalar curve[];
  curvature(cs,curve);

// we refine the mesh
  int count = 4, n;
  for (int k = 0; k < count; k++){
    adapt_wavelet({cs},
    (double[]){1.e-2},MAXLEVEL, 4);

    foreach_vertex()
      distn[] = clamp(circle(x, y, z, R0, A, m), -NB_width, NB_width);
    boundary({distn});
    fractions(distn,cs);
    n= 0;
    foreach(reduction(+:n)){
      n++;
    }
    fprintf(stderr, "prerefine nbcells %d\n", n);

  }
  view (fov = 15.4043, quat = {-0.174473,-0.480791,-0.0947139,0.854064});
  draw_vof("cs");
  save("cs.png");


  scalar dist[];
  foreach(){
    dist[] = clamp(circle(x, y, z, R0, A, m), -NB_width, NB_width);
  }
  boundary({dist});


/**
Advection loop
*/
  count = 100;
  
  for (int k = 0; k < count; k++){
    // in interfacial cellsspeed = -(curve - curve0) where curve0 is equilibrium
   // value
    n= 0;
    foreach(reduction(+:n)){
      n++;
    }
    fprintf(stderr, "nbcells%d\n", n);
    fprintf(stderr, "iteration %d\n", k);
    cell2node(dist,distn);
    fractions(distn,cs,fs);
    curvature_LS(dist,curve);

    stats s =statsf(curve);
    fprintf(stderr, "moyenne curve %g %g\n", s.min, s.max);

    vector vpc[];
    foreach(){
      if(interfacial(point,cs)){
        coord n       = facet_normal( point, cs ,fs);
        normalize(&n);
  // 
        foreach_dimension(){
            vpc.x[] = -n.x*(-15+curve[]);
        }
      }
    }
    boundary((scalar *){vpc});
    restriction((scalar *){vpc});
    stats s2 = statsf(vpc.z);
    fprintf(stderr, "stats vpcz %g %g \n", s2.min, s2.max);

    // recons_speed
    scalar * speed_recons  = {vpc.x,vpc.y,vpc.z};
    double err = 0.;
    recons_speed(dist, DT2, speed_recons,
		 1.e-6, &err, 
		 20, 
		 cs, fs);

    // advect LS
    RK3(dist,vpc,DT2, 0.2);  
    draw_vof("cs");
    cells();
    cells(n = {0,1,0});
    save("cs.mp4");
    LS_reinit(dist);
    foreach(){
      dist[] = clamp(dist[], -NB_width, NB_width);
    }
    boundary({dist});
    restriction({dist});
    adapt_wavelet({cs},(double[]){1.e-2},
      maxlevel = MAXLEVEL, minlevel =4,{dist});
  }
}
