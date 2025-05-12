/**
#Test case for the LS_recons() 


Given a value of a field on interfacial centroids we want to obtain a
cell-centered field such that the bilinear interpolation of the cell-centered
field on the interfacial centroids re-gives the initial value.



*/

#define BICUBIC 1
#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846

#include "embed.h"
#include "fractions.h"
#include "curvature.h"

#include "view.h"
#include "../alex_functions.h"
#include "../LS_recons.h"

#define T_eq          0.

#define DEBUG_LS 1


scalar dist[];
scalar * level_set = {dist};

/**
The field to be reconstructed is `reconsfield`.
*/
scalar reconsfield[];
scalar * LS_speed   = {reconsfield};

int 		nb_cell_NB;
double  NB_width ;    // length of the NB


/**
Initial geometry definition. Here the interface equation is :

$$
r\left(1+ 0.2 *cos(6\theta) \right) - R
$$
where $r = \sqrt{x^2 + y^2}$, $\theta = \arctan(x,y)$ and $R = \frac{L_0}{5}$

*/
double geometry(double x, double y, double Radius) {

  coord center;
  center.x = 0.5;
  center.y = 0.5;

  double theta = atan2 (y-center.y, x-center.x);
  double R2  =  sq(x - center.x) + sq (y - center.y) ;
  double s = ( Radius - sqrt(R2)*(1.+0.25*cos(12*theta)) );

  return s;
}

int main() {
  periodic(right);

  L0 = 1.;
  CFL = 0.5;
  origin (0., 0.);
  
  int MAXLEVEL = 7 ;
  N = 1 << MAXLEVEL;
  init_grid (N);
  nb_cell_NB = 1<<3;
  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  // vertex scalar dist_n[];
  // foreach_vertex() {
  //   dist_n[] = clamp(-geometry(x,y,L0/4.),-1.5*NB_width,1.5*NB_width);
  // }

  // boundary ({dist_n});
  // restriction({dist_n});

  // fractions (dist_n, cs, fs);
  // boundary({cs,fs});
  // restriction({cs,fs});

  foreach() {
    dist[] = clamp(-geometry(x,y,L0/4.),-1.5*NB_width,1.5*NB_width);
  }
  boundary ({dist});
  restriction({dist});
  vertex scalar dist_n[];
  cell2node(dist,dist_n);

  fractions (dist_n, cs, fs);
  fractions_cleanup(cs,fs);
  boundary({cs,fs});
  restriction({cs,fs});


  foreach_face(){
    reconsfield[] = 0.;
  }
  boundary({reconsfield});

/**
The field to be reconstructed is the normal of the interface projected in the
x axis $\mathbf{n_x}$.
*/

	foreach(){
		if(interfacial(point, cs)){
			coord n       = facet_normal( point, cs ,fs);
      normalize(&n);
      foreach_dimension(){
	      reconsfield[] = n.x;
      }
		}
	}
  boundary({reconsfield});
  int err = 0;

/**
To obtain a convergence to $eps = 1\times 10^{-12}$ we remove the limitation
activated by default by setting $k_{limit}$ = 1.
*/
  int k_limit = 1; 
  double deltat  = 0.5*L0 / (1 << MAXLEVEL);  // Delta

  recons_speed(dist, deltat, LS_speed, 
  	k_limit, 6.e-5, &err, 1000, 0.3, cs ,fs);
  view (fov = 24, tx = -0.487267, ty = -0.444062, 
        width = 800, height = 800);
  dump();
  draw_vof("cs");
  squares("reconsfield");
  save ("reconsfield.png");
  if(err==1){
  	fprintf(stderr, "NOT CONVERGED\n" );
  }
  exit(1);

  

}

/**
![The obtained reconstructed field](test_recons/reconsfield.png)


*/