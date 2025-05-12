/**
# Engulfing droplet : Instable triple point

![Basic three-phase VOF](instable/movie_basic.mp4)




# Engulfing droplet : Instable triple point

![Corrected advection scheme](instable/movie_corrected.mp4)


 */



//#define Corrected 1

#include "axi.h"
#include "navier-stokes/centered.h"
#if Corrected
#include "headers/three-phase_corrected.h"
#else
#include "headers/three-phase.h"
#endif
#include "tension.h"
#include "headers/reduced3.h"
#include "view.h"
#include "tag.h"
#include "headers/norm-3p.h"
//#include "navier-stokes/perfs.h"

//#define FILTERED

double sigma12 = 0.075;  //water/gas
double sigma13 = 0.035;  //oil/gas
double sigma23 = 0.025;  //water/oil 

#define beta acos((sq(sigma13)-sq(sigma12)-sq(sigma23))/(2*sigma12*sigma23))/2

#define SIGMA12 1.
#define SIGMA13 (sigma13/sigma12)
#define SIGMA23 (sigma23/sigma12)

#define vol (4./3)*pi*pow(radius,3)
#define den ((4./3)*pi*pow(radius,3))/((2*pi*sq(1.-cos(beta))*(2.+cos(beta)))/(24.*pow(sin(beta),3)))
#define D pow(den,(1./3))

// default : box size define equals to 1
#define length 1.
#define radius 0.15
#define radius_exact D/(2*sin(beta)) // numerical raidus equals to 25% of the box size
// Define all physical parameters

#define Rayon_physique 1.5e-3
#define hauteur_physique Rayon_physique/radius*length/2
#define hauteur_eau 0. // don't put the drop outside of the box !
#define hauteur_huile 5.e-3
// We choose a dimensionless speed equals to 1


#define Mu_w 1./60 // Pa.s
#define Mu_g 1./60 // Pa.s
#define Mu_o 1./60
//double Mu_o=1/60;
#define rho_w 1.
#define rho_g 1.
#define rho_o 1.

#define Bond (rho_w*9.81*sq(2.*Rayon_physique)/sigma12)
#define Oh (Mu_w/sqrt(sigma12*rho_w*2.*Rayon_physique))

#define vrho_w 1.
#define vrho_g (rho_g/rho_w)
#define vrho_o (rho_o/rho_w)
#define vmu_w (sqrt(2*radius)*Oh)
#define vmu_g (sqrt(2*radius)*Oh*(Mu_g/Mu_w))
#define vmu_o (sqrt(2*radius)*Oh*(Mu_o/Mu_w))

// define convert ratio between physical and numerical values

#define ratio radius/Rayon_physique // convert all geometrical lengths in numerical value
// the phisical domain of calcul is define by length/ratio



// Convert all real values in numerical values

#define hnum hauteur_physique*ratio
#define hwater hauteur_eau*ratio
#define hoil hauteur_huile*ratio

#if AXI
#define water1_exact(x,y) (-sq(x-0.5+radius_exact*(-cos(beta)))-sq(y)+sq(radius_exact))
#define water2_exact(x,y) (-sq(x-0.5+radius_exact*(+cos(beta)))-sq(y)+sq(radius_exact))
#define water(x,y) (-sq(x-0.5)-sq(y)+sq(radius))
#define water_exact(x,y) min(water1_exact(x,y),water2_exact(x,y))
#else
#define water(x,y) (-sq(x-0.5)-sq(y-0.5)+sq(radius))
#endif
#define pool(x,y) hwater+hoil-x
#define domain(x,y) min(water(x,y),pool(x,y))
#define epsilon 0.95e-4*ratio
#define rfwater(x,y) (sq(radius-epsilon)<(sq(x-hnum)+sq(y)) && (sq(x-hnum)+sq(y)<sq(radius+epsilon)))


#if AXI
u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);
#endif


int maxlevel = 7;
int minlevel = 5;

double dtuser = 0.001;
double tmax = .7;
double uem = 0.001;
double fem = 0.01;

bool reboot=false;
scalar f1e[], f2e[], f3e[],f6[];
double erreur=0.;

stats inits1;
stats inits2;
stats inits3;
stats inits4;

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
  G.x=-Bond/sq(2*radius);
  size(length);
  run();     
}

event init(i=0){
  if (!restore (file = "restart")){
    reboot=true;
    refine((rfwater(x,y)) && (level<maxlevel));
    fraction(f2,water(x,y));
    #if AXI
    foreach(){
      if(x <= length/2.){
        f3[] = 1. - f2[];
        f1[] = 0.;
      }
      else{
        f1[] = 1. - f2[];
        f3[] = 0.;
      }
    }
    #else
    foreach(){
      if(y <= length/2.){
        f3[] = 1. - f2[];
        f1[] = 0.;
      }
      else{
        f1[] = 1. - f2[];
        f3[] = 0.;
      }
    }
    #endif
    //boundary ((scalar *){f1,f2,f3,rhov});
  }
  inits1=statsf(f1);
  inits2=statsf(f2);
  inits3=statsf(f3);
}


event remove_drop(i+=5){
  remove_droplets(f1,4);
  remove_droplets(f2,4);
  remove_droplets(f3,4);
}

event adapt (i++) {
  adapt_wavelet ({f1,f2,f3}, (double[]){fem,fem,fem}, maxlevel,minlevel);
}




event movie (t =0.; t +=dtuser; t<= tmax) {
  clear();
  box();
  view (fov=10, quat = {0,0,-0.707,0.707}, tx =.27, ty =-0.4);
  draw_vof ("f1", filled=1, fc = {0.,0.8,0.8}, lw = 2);
  draw_vof ("f2",filled=1, fc = {1.,0.,0.}, lw = 2);
  draw_vof ("f3",filled=1, fc = {0.7,0.7,0.}, lw = 2);
  mirror({0,-1}){
    draw_vof ("f1", filled=1, fc = {0.,0.8,0.8}, lw = 2);
    draw_vof ("f2",filled=1, fc = {1.,0.,0.}, lw = 2);
    draw_vof ("f3",filled=1, fc = {1.,1.,0.}, lw = 2);
    #if Corrected
    save ("movie_corrected.mp4");
    #else
    save("movie_basic.mp4");
    #endif
  }
}

event logfile(t=0.; t +=dtuser; t<= tmax){
  double ke1=0.,ke2=0.,ke3=0.;
  stats s1=statsf(f1);
  stats s2=statsf(f2);
  stats s3=statsf(f3);


  foreach(reduction(+:ke1) reduction(+:ke2) reduction(+:ke3)){
    foreach_dimension(){
      ke1+=dv()*f1[]*rho_w*sq(u.x[]);
      ke2+=dv()*f2[]*rho_g*sq(u.x[]);
      ke3+=dv()*f3[]*rho_o*sq(u.x[]);
    }
  }
  fraction(f2e,water_exact(x,y));
  foreach(){
    if(x <= length/2.){
      f3e[] = 1. - f2e[];
      f1e[] = 0.;
    }
    else{
      f1e[] = 1. - f2e[];
      f3e[] = 0.;
    }
    erreur+=dv()*(fabs(f1[]-f1e[])+fabs(f2[]-f2e[])+fabs(f3[]-f3e[]));
  }

  if (t==0){
    FILE * fp500 = fopen("facets_exact","w"); 
    output_facets(f1e, fp500);
    output_facets(f2e, fp500);
    output_facets(f3e, fp500);
    fclose(fp500);
  }
  FILE * fp200 = fopen("facets","w"); 
  FILE * fp300 = fopen("cells","w"); 
  FILE * fp100 = fopen("data","a");
  coord p1;
  struct jonction j1 = locate_triple_point(&p1, f1, f2, f3, 0., length, 0., length, 1.1);
  output_cells(fp300);
  output_facets(f1, fp200);
  output_facets(f2, fp200);
  output_facets(f3, fp200);

  fprintf(fp100, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", t,j1.pt12x,j1.pt12y,j1.pt13x,j1.pt13y,
  j1.pt23x,j1.pt23y, p1.x,p1.y, fractions_error(f1, f2, f3, sq(L0)), shape_error(f1, f2, f3, f1e, f2e, f3e, sq(L0)),s1.sum/inits1.sum,s2.sum/inits2.sum,s3.sum/inits3.sum,(s1.sum+s2.sum+s3.sum)/(inits1.sum+inits2.sum+inits3.sum));
  fclose(fp100);
  fclose(fp200);
  fclose(fp300);
}
