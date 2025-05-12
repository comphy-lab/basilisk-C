/**
Impact imposing two different contact angles. One while spreading and another one while receding
*/

#define theta_adv (2.0944)
#define theta_rec (pi/4)
#define theta_0 (pi/2)
#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof_angle.h"
#include "tension_angle.h"
#include "tag.h"
#include "droplet_stat.h"

# define WE 32 //based in the droplet radius
# define RHO_r 0.001 
# define RE 1673 
# define MU_r 0.018

int LEVEL = 8;

scalar f[];
scalar * interfaces = {f};
face vector alphav[];
scalar rhov[];
face vector muv[];

u.n[right]  = neumann(0);
uf.n[right] = neumann(0);
p[right] = neumann(0);

u.t[left]   = dirichlet(0);
u.n[left]   = dirichlet(0);
uf.t[left]  = dirichlet(0);

u.n[top]  = neumann(0);
uf.n[top] = neumann(0);
p[top] = dirichlet(0);

uf.n[bottom] = 0;

#define rho(f) (clamp((f), 0, 1)*(1-RHO_r) + RHO_r)
#define mu(f)  (((f) + (1. - (f))*MU_r)/RE)
//#define mu(f) (1./((f) + (1. - (f))/MU_r)/RE)



int main() {
  size(4);
  init_grid (1 << 4);
  alpha = alphav;
  rho = rhov;
  mu = muv;
  f.sigma = 1.0/WE;
  run();
}

#if 1
#define ff f
#else
scalar ff[];
#endif

event init (t = 0) {

  refine (sq(x-1. - 0.2) + sq(y) - sq(1 + 0.02) < 0 && 
	 (sq(x-1. - 0.2) + sq(y) - sq(1 - 0.02) > 0 || x < 0.1) &&
	  level < LEVEL);
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = -sq(x-1 - 0.2) - sq(y) + 1.;
  fractions (phi, f);
  foreach() {
    u.x[] = -f[];
#ifndef ff
    ff[] = f[];
#endif
  }
}

event properties (i++) {

#if TREE  
  f.prolongation = refine_bilinear;
  boundary ({f});
#endif

#ifndef ff
  foreach()
    ff[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
  boundary ({ff});
#endif
  foreach_face() {
    double T = (ff[] + ff[-1])/2.;
    alphav.x[] = fm.x[]/rho(T);
    muv.x[] = fm.x[]*mu(T);
  }
  foreach()
    rhov[] = cm[]*rho(ff[]);

#if TREE  
  f.prolongation = fraction_refine;
  boundary ({f});
#endif
  boundary((scalar *){alphav, muv, rhov});
  
}

event adapt (i += 1) {
  adapt_wavelet ({f, p, u.x, u.y}, (double[]){1e-3, 1e-3, 1e-3, 1e-3},
		 maxlevel = LEVEL);
}

event drop_remove (i += 5) {
  remove_droplets (f, 1, 0);
}

event angle(i++){
  foreach_boundary(left)
	{
    	if (f[] < 1. && f[] > 0.)
      		{
      		if (u.y[] > 0.01){
      		        theta = theta_adv;
        
                        if (theta > 2.35619)
                                theta = 2.35619;
                        else if (theta < pi/4.)
                                theta = pi/4;
                                
                                }
                else if (u.y[] < 0.01){
                      		        
                       theta = theta_rec;
        
                        if (theta > 2.35619)
                                theta = 2.35619;
                        else if (theta < pi/4.)
                                theta = pi/4;
                                }
                else{
                        theta = theta_0;
                    }         
                }
          }
}

event profile(t = 5e-2; t < 20; t+= 5e-2) {
  char name[80];
  sprintf(name,"L8profile_t%f.dat", t);
  FILE * fp1 = fopen(name,"w");
  output_facets (f, fp1);
}

#if 0
event gfsview (i += 4) {
  static FILE * fp = popen("gfsview2D", "w");
  output_gfs (fp, t = t);
}
#endif

#if 0
event movies (i += 2) {
  char name[80];
  sprintf(name,"gfsview-batch2D detail.gfv | ppm2mpeg > impact.mpg");
  static FILE * fp = popen (name, "w");
  output_gfs (fp, t = t);
  fprintf (fp, "Save stdout { format = PPM width = 950 height = 410}\n");
}
#endif

