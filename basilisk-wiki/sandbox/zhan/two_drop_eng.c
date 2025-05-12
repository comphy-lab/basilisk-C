/**
# Two drops engulfing critical test case

The test case is taken from [Yu et al, 2019](#yu_versatile_2019)

![Time evolution of the interface with $\sigma_{12}$ = $\sigma_{31}$ = 0.526](two_drop_eng/movie.mp4)
*/

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

#define R2 			0.5		// Radius of the bubbles
#define R3 			0.666667
#define L 			2.5 	//size of the square domain
#define OH 			0.1 	// Ohnesorge number	
#define MAXLEVEL 	6
#define MINLEVEL 	3

vertex scalar phi[]; 				//distance function used to define the fractions


int main(){
  size (L);
  periodic(left);
  periodic(top);
#if AXI
#else
  origin (-L/2.0, -L/2.0);
#endif
  N = 1 << MAXLEVEL;
  init_grid (N);
  double sigma12 = 1.0;
  double sigma31 = 0.526;
  double sigma23 = sigma31;
  f1.sigma 	= 0.5 * (sigma12 + sigma31 - sigma23);
  f2.sigma 	= 0.5 * (sigma12 + sigma23 - sigma31);
  f3.sigma 	= 0.5 * (sigma23 + sigma31 - sigma12);
  rho1 		= 1.0;
  rho2 		= 1.0;
  rho3 		= 1.0;
  mu1 		= OH*1.0;
  mu2			= OH*1.0;
  mu3 		= OH*1.0;
  run();
  // movie generation is broken by the return below. This is a bug of Basilisk.
  //  return 0; 
}

event init (t = 0){
  if(!restore("dumpfile")){
    foreach_vertex()
      phi[] = ((sq(x) + sq(y) + sq(z))) - sq(R3) ;
    boundary({phi});
    restriction({phi});
    fractions(phi, f1);

    foreach_vertex()
      phi[] = sq(R2) - ((sq(x) + sq(y-R3+R2) + sq(z)));
    boundary({phi});
    restriction({phi});
    fractions(phi, f2);
    foreach()
      f3[] = 1. - f1[] - f2[];
    boundary({f1, f2, f3});
    restriction({f1, f2, f3});
  }
  else{}
}

event adapt(i++){
  adapt_wavelet({f1, f2, f3, u.x, u.y}, (double[]) {1.0e-9, 1.0e-9, 1.0e-9, 1.0e-3, 1.0e-3}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event movie (t += 0.1; t <= 30.)
{
  clear();
  box();
  cells();
  draw_vof("f1");
  draw_vof("f2", filled = 1, fc = {0.9, 0.4, 0.2});  
  draw_vof("f2");
  draw_vof("f3", filled = 1, fc = {0.2, 0.4, 0.9});
  draw_vof("f3");
  save("movie.mp4");
}

/**
## References

~~~bib
@article{yu_versatile_2019,
   title = {A versatile lattice {Boltzmann} model for immiscible ternary fluid flows},
   volume = {31},
   issn = {1070-6631, 1089-7666},
   url = {http://aip.scitation.org/doi/10.1063/1.5056765},
   doi = {10.1063/1.5056765},
   abstract = {We propose a lattice Boltzmann color-gradient model for immiscible ternary ﬂuid ﬂows, which is applicable to the ﬂuids with a full range of interfacial tensions, especially in near-critical and critical states. An interfacial force for N-phase systems is derived and then introduced into the model using a body force scheme, which helps reduce spurious velocities. A generalized recoloring algorithm is applied to produce phase segregation and ensure immiscibility of three different ﬂuids, where an enhanced form of segregation parameters is derived by considering the existence of Neumann’s triangle and the effect of the equilibrium contact angle in a three-phase junction. The proposed model is ﬁrst validated by two typical examples, namely, the Young-Laplace test for a compound droplet and the spreading of a droplet between two stratiﬁed ﬂuids. It is then used to study the structure and stability of double droplets in a static matrix. Consistent with the theoretical stability diagram, seven possible equilibrium morphologies are successfully reproduced by adjusting the interfacial tension ratio. By simulating near-critical and critical states of double droplets where the outcomes are very sensitive to the model accuracy, we show that the present model is advantageous to three-phase ﬂow simulations and allows for accurate simulation of near-critical and critical states. Finally, we investigate the inﬂuence of interfacial tension ratio on the behavior of a compound droplet in a three-dimensional shear ﬂow, and four different deformation and breakup modes are observed.},
   language = {en},
   number = {1},
   urldate = {2021-07-10},
   journal = {Physics of Fluids},
   author = {Yu, Yuan and Liu, Haihu and Liang, Dong and Zhang, Yonghao},
   month = jan,
   year = {2019},
   pages = {012108},
   file = {2019 - A versatile lattice Boltzmann model for immiscible.pdf:/home/xue/Zotero/storage/YCLPHDJW/2019 - A versatile lattice Boltzmann model for immiscible.pdf:application/pdf},
   }
~~~
*/
