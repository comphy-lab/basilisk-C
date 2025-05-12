/**
#Test case for the LS_recons() + LS_reinit() + curvature_LS()

Given a value of a field on interfacial centroids we want to obtain a
cell-centered field such that the biquadratic interpolation of the cell-centered
field on the interfacial centroids re-gives the initial value.

Here we check the backbone of the hybrid level-set embedded method, with the
redistancing, the reconstruction off the interface and we do a comparison of the
curvature calculated using the level-set function and the VOF function.
This test case is similar to the test cases of Popinet 
[meanflow.c](/sandbox/popinet/meanflow.c) and Magdelaine 
[curvature_equilibrium.c](/sandbox/qmagdelaine/equilibrium_shapes/curvature_equilibrium.c)
we solve
$$
\mathbf{n}\cdot\partial_\tau\mathbf{X} = -\kappa
$$
taken from [Deckelnick et al,
2005](https://homepages.warwick.ac.uk/~masgat/PAPERS/DecDziEll05_ActaNumerica.pdf).

~~~gnuplot Evolution of the interface (zoom)
set term svg size 800,800
set output 'interface.svg'
set size ratio -1
set key top right
unset xlabel
unset xtics
unset ytics
unset border
unset margin
unset ylabel
plot 'out0' w l lw 2 t 'curve VOF',\
  'out1' w l lw 2 t 'curve LS'
~~~
*/

#define BIQUADRATIC 1
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846
#define LS_face_value(a,i)      ((a[i] + a[i-1])/2.)

#include "embed.h"
#include "advection.h"

#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../alex_functions.h"
#include "../LS_recons.h"
#include "../LS_curvature.h"
#include "../LS_advection.h"

scalar * tracers = NULL;
scalar dist[];
scalar * level_set = {dist};

/**
The field to be reconstructed is `reconsfield`.
*/
vector reconsfield[];
scalar * advectionfield   = {reconsfield.x, reconsfield.y};

int     nb_cell_NB;
double  NB_width ;    // length of the NB

char filename [100];
FILE * fp1;

double geometry(double x, double y) {
  double theta = atan2(y,x);
  double r = sqrt(sq(x) + sq(y));
  return r - 0.4*(1. + 0.2*cos(8.*theta));

}

int main() {
  periodic(right);

  L0 = 1.;
  CFL = 0.5;
  origin (-0.5, -0.5);
  
  int MAXLEVEL = 6 ;
  N = 1 << MAXLEVEL;
  init_grid (N);
  nb_cell_NB = 1<<3;
  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  boundary ({dist});
  restriction({dist});
  for (int j = 0; j<2; j++){
    foreach() {
    dist[] = -geometry(x,y);
    }
    boundary({dist});
    snprintf(filename, 100,  "out%d", j);
    fp1 = fopen (filename,"w");
    for (int i = 0; i<400; i++){
      vertex scalar dist_n[];
      cell2node(dist,dist_n);
  
      fractions (dist_n, cs, fs);
      fractions_cleanup(cs,fs);
      boundary({cs,fs});
      restriction({cs,fs});
  
      foreach_face(){
        reconsfield.x[] = 0.;
      }
      boundary((scalar *){reconsfield});
  
      scalar curve[];
      if(j==0)
        curvature(cs,curve);
      else
        curvature_LS(dist,curve);
  
      boundary({curve});
      foreach(){
        if(interfacial(point, cs)){
          coord n       = facet_normal( point, cs ,fs);
          normalize(&n);
  // 
          foreach_dimension(){
            reconsfield.x[] = -n.x*(/*4.+*/curve[]);
          }
        }
      }
      boundary((scalar *){reconsfield});
      restriction((scalar *){reconsfield});
  
      int err = 0;
    
      recons_speed(dist, dt= 0.5*L0 / (1 << MAXLEVEL), 
        advectionfield, 1.e-5, &err, 60,  cs ,fs);
  
      face vector vpcf[];
      foreach_face()
        vpcf.x[] = LS_face_value(reconsfield.x,0);
      boundary((scalar *){vpcf});
      restriction((scalar *){vpcf});
  /**
  Hardcoded timestep.
  */
      double dt = 1.8e-5;
      advection_no_metric (level_set, vpcf, dt);    
      LS_reinit(dist);
  
      if(i%30 == 0)output_facets (cs, fp1);
    }
    fclose(fp1);
  }
}

