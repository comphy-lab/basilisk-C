/**
 * We simulate the impact of a water drop on thin oil film. This work is in collaboration with Alidad Amirfazli from Surface Engineering and Instrumentation Lab (SEiL).
*/


#include "navier-stokes/centered.h"
#include "axi.h"
#include "three-phase.h"
#include "tension.h"
#include "tracer.h"
#include "conserving3.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "tag.h"
  
#define FILTERED
// default : box size define equals to 1
#define length 1.
#define radius 0.25 // numerical raidus equals to 25% of the box size
// Define all physical parameters

#define Rayon_physique 1.5e-3
#define hauteur_physique 2.e-3
#define hauteur_eau 0. // don't put the drop outside of the box !
#define hauteur_huile 3.e-4
// We choose a dimensionless speed equals to 1
#define V 0.5

#define Mu_w 0.002 // Pa.s
#define Mu_g 1.6e-4 // Pa.s
#define Mu_o 0.02
#define rho_w 1000.
#define rho_g 2.
#define rho_o 900.

#define Re (rho_w*V*Rayon_physique*2./Mu_w)

#define vrho_w 1.
#define vrho_g (rho_g/rho_w)
#define vrho_o (rho_o/rho_w)
#define vmu_w (2.*radius/Re)
#define vmu_g (2.*radius/Re*(Mu_g/Mu_w))
#define vmu_o (2.*radius/Re*(Mu_o/Mu_w))


#define We (rho_w*sq(V)*Rayon_physique*2./sigma12)

#define sigma12 7.e-2
#define sigma13 3.2e-2
#define sigma23 2.e-2


#define SIGMA12 2.*radius/We
#define SIGMA13 2.*radius/We*(sigma13/sigma12)
#define SIGMA23 2.*radius/We*(sigma23/sigma12)




// define convert ratio between physical and numerical values

#define ratio radius/Rayon_physique // convert all geometrical lengths in numerical value
// the phisical domain of calcul is define by length/ratio



// Convert all real values in numerical values

#define hnum hauteur_physique*ratio
#define hwater hauteur_eau*ratio
#define hoil hauteur_huile*ratio

#define water(x,y) (-sq(x-hnum)-sq(y)+sq(radius))
#define pool(x,y) hwater+hoil-x
#define domain(x,y) max(water(x,y),pool(x,y))
#define epsilon 0.95e-4*ratio
#define initwater(x,y) (-sq(x-hnum)-sq(y)+sq(radius+epsilon))
#define rfwater(x,y) (sq(radius-epsilon)<(sq(x-hnum)+sq(y)) && (sq(x-hnum)+sq(y)<sq(radius+epsilon)))

#define gad 10/(ratio*sq(V))


u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);


int maxlevel = 9;
int minlevel = 6;

double dtuser = 0.001;
double tmax = 0.25; //0.7 (simulation is too long);
double uem = 0.01;
double fem = 0.001;
double volref[3];
face vector av[];

scalar b1[],b2[],b3[],f[];
scalar * tracers = {f}; 

char comm1[80],comm2[80],comm3[80],tempo[160];
int h, min, day, mois, an,min,h;
time_t now;

int main(){
	origin(0.,0.);
	init_grid(256);
	rho1=vrho_g;
	rho2=vrho_w;
	rho3=vrho_o;
	mu1=vmu_g;
	mu2=vmu_w;
	mu3=vmu_o;
	f2.sigma = (SIGMA12 + SIGMA23 - SIGMA13)/2.; //eau
        f3.sigma = (SIGMA23-SIGMA12+SIGMA13)/2.; // oil
        f1.sigma = (SIGMA13 + SIGMA12 - SIGMA23)/2.; //gaz
	a=av;
	size(length);
	run();		

}

event init(i=0){
	refine((rfwater(x,y)) && (level<maxlevel));
	scalar f4[],f6[],f5[];

  fraction(f4,water(x,y));
  fraction(f5,initwater(x,y));
  fraction(f6,pool(x,y));
  foreach(){
    f2[]=f4[];
    f3[]=f6[];
    f1[]=clamp(1.-f2[]-f3[],0.,1.);
    rhov[]=rho(f2[],f3[]);
    f[]=f4[]+5*f6[];
  }
  boundary ((scalar *){f1,f2,f3,rhov});

	foreach(){
		u.x[]=-1*f5[];
  }
	
}


event acceleration(i++){
  foreach_face(x){
    av.x[]=-gad;
  }
  boundary((scalar*) {av});
}

event remove_drop(i+=5){
  remove_droplets(f1,4);
  remove_droplets(f2,4);
  remove_droplets(f3,4);
}

event adapt (i++) {
	
        adapt_wavelet ({u,f1,f2,f3}, (double[]){uem,uem,fem,fem,fem}, maxlevel,minlevel);
}



event movie (t =0.; t +=dtuser; t<= tmax) {
  clear();
  view (fov=20, quat = {0,0,-0.707,0.707}, tx =0.0, ty =-0.5, width = 2000, height = 2000);
  draw_vof ("f1",filled=1, fc = {0.,0.8,0.8}, lw = 2);
  draw_vof ("f2",filled=1, fc = {1.,0.,0.}, lw = 2);
  draw_vof ("f3",filled=1, fc = {0.7,0.7,0.}, lw = 2);
  mirror({0,-1}) {
  draw_vof ("f1", lw = 2);
  draw_vof ("f2", lw = 2);
  draw_vof ("f3", lw = 2);
  cells();
  }
  save("drop_impact_level9.mp4");

}
/**
 * The solution presented below is not well resolve, I recommand a maxlevel=11 to kill all the non physical effects we can observe here.
 * If the code is running in order to do some experiment a maxlevel of 13 is recommanded.
*/
/**
![Animation of the splash -level 9](drop_on_oil/drop_impact_level9.mp4)
*/