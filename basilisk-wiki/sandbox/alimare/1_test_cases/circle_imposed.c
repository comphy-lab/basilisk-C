/**
# Impose $speed = n$ on a circle

We simply check that the initial circle stays a circle as it grows by a ratio
250.

~~~gnuplot Curvature
plot 'log' u 1:2 w l t 'min', '' u 1:3 w l t 'max'
~~~

~~~gnuplot Interface
set size ratio -1
set object 1 circle front at 0.,0. size 0.45 fillcolor rgb "black" lw 1 dt 2
plot 'out' w l t ''
~~~



*/
#define QUADRATIC 1
#define BGHOSTS 2
#include "profiling.h"
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

FILE * fp1;
char filename [100];

int main()
{
  origin(-0.5,-0.5);
  int MAXLEVEL = 8;
  // this is the diffusive timestep limit (with diffusion coeff unity)
  init_grid(1<< (MAXLEVEL-2));

  vertex scalar distn[];
  coord center = {0.,0.};

  double NB_width = (1<<2) * L0/(1<< MAXLEVEL);
  double size = L0/119;

  foreach_vertex() {
    distn[] = clamp(circle(x,y,center, size), -NB_width, NB_width);
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
    distn[] = clamp(circle(x,y,center, size), -NB_width, NB_width);
    boundary({distn});
    fractions(distn,cs);
  }

  scalar dist[];
  foreach(){
    dist[] = clamp(circle(x,y,center, size), -NB_width, NB_width);
  }
  boundary({dist});
  curvature(cs,curve); 
/**
Advection loop
*/
  double myt = 0., t_fin = 4.5;
  double mydt = 0.5*L0/(1 << MAXLEVEL);
  int k=0;

  // fopen(fp1, "init_interface","w");
  // output_facets(cs,fp1);
  // fclose(fp1)

  while (myt<t_fin)
  { // mimics a time loop

    myt+= mydt;

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
          vpc.x[] = -0.1*n.x;
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
    scalar * speed_recons  = {vpcr.x,vpcr.y};
    double err = 0.;

    recons_speed(dist, deltat = 0.45*L0/(1<<MAXLEVEL), speed_recons,
     tolerance = 1.e-6, &err, 
     nb_iter = 30, 
     cs, fs,NB_width);

  // advect LS
    RK3(dist,vpcr,mydt, NB_width);
    if(k%75 ==0){    
      output_facets(cs,stdout);
    }
    k++;
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
