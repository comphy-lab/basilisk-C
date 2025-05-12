/**
# Imposing $speed = n$ on a sphere

We simply check that the initial sphere stays a sphere as it grows

![Evolution of the 0-level-set](sphere_imposed/cs.mp4)

~~~gnuplot Curvature
set yrange[-60:0]
f(x) = -1/(0.1*x)
set xlabel 'time'
set ylabel 'Curvature'
plot 'log' u 1:2 w l t 'min', '' u 1:3 w l t 'max', f(x) t '1/t'
~~~


*/
#define QUADRATIC 1
#define BGHOSTS 2
#define DEBUG 1
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
#include "../basic_geom.h"

int main()
{
  origin(-0.5,-0.5,-0.5);
  int MAXLEVEL = 7;
  // this is the diffusive timestep limit (with diffusion coeff unity)
  init_grid(1<< 6);

  vertex scalar distn[];
  coord center = {0.,0.,0.};

  double NB_width = 4*L0/(1<<MAXLEVEL);
  double size = L0/60.;
  foreach_vertex() {
      distn[] = clamp(sphere(x,y,z,center, size), -NB_width, NB_width);
  }
  boundary({distn});
  fractions(distn,cs,fs);
  
  scalar curve[];
  curvature(cs,curve);

// we refine the mesh
  int count = 4;
  for (int k = 0; k < count; k++){
    adapt_wavelet({cs},
    (double[]){1.e-2},MAXLEVEL, 4);
    foreach_vertex()
      distn[] = clamp(sphere(x,y,z,center, size), -NB_width, NB_width);
    boundary({distn});
    fractions(distn,cs);
  }
  


  scalar dist[];
  foreach(){
    dist[] = clamp(sphere(x,y,z,center, size), -NB_width, NB_width);
  }
  boundary({dist});


/**
Advection loop
*/
  count = 200;
  double myt = 0., t_fin = 1.;
  double mydt = 0.45*L0/(1 << MAXLEVEL);
  int k=0;   
  char name[80];
  view (fov = 15.4043, quat = {-0.174473,-0.480791,-0.0947139,0.854064});
  draw_vof("cs");
  cells();
  cells(n = {0,1,0});
  sprintf (name, "cs%g.png", myt);
  save(name);

  while (myt<t_fin){
    k++;
    myt+= mydt;
    // in interfacial cellsspeed = -(curve - curve0) where curve0 is equilibrium
   // value
    cell2node(dist,distn);
    fractions(distn,cs,fs);
    curvature_LS(dist,curve);

    stats s =statsf(curve);
    fprintf(stderr, "%g %g %g\n", myt, s.min, s.max);

    vector vpc[];
    foreach(){
      if(interfacial(point,cs)){
        coord n       = facet_normal( point, cs ,fs);
        normalize(&n);
  // 
        foreach_dimension(){
            vpc.x[] = -0.3*n.x;
        }
      }
      else
        foreach_dimension()
          vpc.x[] = 0.;
    }
    boundary((scalar *){vpc});
    restriction((scalar *){vpc});

    // recons_speed
    vector vpcr[];
    foreach(){
      foreach_dimension(){
        if(interfacial(point,cs))vpcr.x[] = vpc.x[];
        else vpcr.x[] = 0.;
      }
    }
    boundary((scalar * ){vpcr});
    restriction((scalar * ){vpcr});
    scalar * speed_recons  = {vpcr.x,vpcr.y,vpcr.z};
    double err = 0.;
    recons_speed(dist, deltat = 0.45*L0/(1<<MAXLEVEL), speed_recons,
     tolerance = 1.e-6, &err, 
     nb_iter = 30, 
     cs, fs,NB_width);

    // advect LS
    RK3(dist,vpcr,mydt, NB_width);

    if(k%130 ==0){    
      draw_vof("cs");
      cells();
      cells(n = {0,1,0});
      sprintf (name, "cs%g.png", myt);
      save(name);
    }
    LS_reinit(dist);
    foreach(){
      dist[] = clamp(dist[], -NB_width, NB_width);
    }
    boundary({dist});
    restriction({dist});
    scalar normvpc[];
    foreach(){
    if(fabs(dist[]< NB_width)){
      coord temp;
      foreach_dimension()
        temp.x = vpc.x[];

      normvpc[] = norm2(temp);
    }
    else{
      normvpc[] = 0.;
    }
  }

  boundary({normvpc});
  restriction({normvpc});

    adapt_wavelet({cs,normvpc},(double[]){1.e-2,1.e-2},
      maxlevel = MAXLEVEL, minlevel =4,{dist});
  }
}
