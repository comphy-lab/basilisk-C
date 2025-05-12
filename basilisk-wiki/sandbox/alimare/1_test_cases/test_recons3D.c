/**
#Test case for the LS_recons() 


Given a value of a field on interfacial centroids we want to obtain a
cell-centered field such that the bilinear interpolation of the cell-centered
field on the interfacial centroids re-gives the initial value.



*/

#define DEBUG 1
#define QUADRATIC 1


#define DOUBLE_EMBED 1
#define LevelSet     1
#define Gibbs_Thomson 1
#define Pi 3.14159265358979323846

#include "grid/octree.h"
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

int         nb_cell_NB;
double  NB_width ;    // length of the NB


/**
Initial geometry definition. Here the interface equation is :

$$
r\left(1+ 0.2 *cos(6\theta) \right) - R
$$
where $r = \sqrt{x^2 + y^2}$, $\theta = \arctan(x,y)$ and $R = \frac{L_0}{5}$

*/
double geometry(Point point, double Radius) {

// spherical coordinates r theta phi
  coord center = {0.5,0.5,0.5};
  double theta = atan2 (sqrt(sq(x-center.x)+sq(y-center.y)),z-center.z);
  double psi   = atan2 (y-center.y, x-center.x);
  double R2 = sqrt(sq (x - center.x) +  sq (y - center.y) + sq (z - center.z));
  // double psi   = acos (z-center.z/R2+eps);
  if(x-center.x == 0)exit(1);
  // return  Radius - sq(x*cos(2*theta)*cos(2*psi))/3. + sq(y)/2. + sq(z)/1.;
  return  Radius - R2*(1.-0.3*cos(4*theta));
}

void draw_isolines(scalar s, double smin, double smax, int niso, 
  int w){
  scalar vdist[];
  cell2node(s,vdist);

  boundary ({vdist});
  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso){
    isoline ("vdist", sval, lw = w);
  }
}


int main() {
  periodic(right);

  L0 = 1.;
  CFL = 0.5;
  origin (0., 0. , 0.);
  for(int kk=0;kk<=0;kk++){ 
  int MAXLEVEL = 6+kk ;
  N = 1 << MAXLEVEL;
  init_grid (N);
  nb_cell_NB = 1<<4;
  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  foreach() {
    dist[] = -geometry(point,L0/4.);
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
      reconsfield[] = n.x;
        }
    }
  boundary({reconsfield});
  double tolerance = 1.e-6, mymax = 0.;
  stats s = statsf(reconsfield);
  mymax = max(mymax,max(fabs(s.min),s.max)) ;
  double err = 0.;
  tolerance = tolerance*mymax;

/**
To obtain a convergence to $eps = 1\times 10^{-12}$ we remove the limitation
activated by default by setting $k_{limit}$ = 1.
*/
  int k_limit = 0; 
  double deltat  = 0.45*L0 / (1 << MAXLEVEL);  // Delta

  recons_speed(dist, deltat, LS_speed, 
    tolerance, &err, 60, cs ,fs);
  fprintf(stderr, "%d %g\n",1<<MAXLEVEL, err);
  view (fov = 30, tx = -0.487267, ty = -0.444062, 
        width = 1400, height = 1400);
  dump();
  // draw_isolines(reconsfield, -0.999, 0.999, 30, 1);
  char filename[100];
  snprintf(filename, 100,  "reconsfield%d.png", kk);
  squares("reconsfield", min = -1, max = 1);
  save (filename);
  if(err==1){
    fprintf(stderr, "NOT CONVERGED\n" );
  }
  dump();
}
  exit(1);
}

/**
![The obtained reconstructed field](test_recons/reconsfield.png)


*/