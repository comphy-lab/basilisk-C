/**
Impact of a drop of water on an oil film above a deep water pool
*/

/**
This simulation aim to observe if a jet is fomed after the impact and compare \n
it's velocity with the one observe during the experiment. All the experimental work has been done by Guy-Jean Michon at Institut d'Alembert. 
*/


#include "navier-stokes/centered.h"
#include "axi.h"
#include "three-phase.h"
#include "tension.h"
#include "conserving3.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "tag.h"


#define FILTERED

#define radius 0.2
#define Dre 2.8e-3
#define Vre 1.935
#define Mu_re2 2.0e-3
#define Mu_re3 5.e-2
#define Rho2_re 1000.
#define Rho3_re 900.
#define Rho1_re 1.21
#define Re (Rho2_re*Vre*Dre/Mu_re2)
#define Sig_re2 72e-3
#define Sig_re3 22e-3
#define Sig_re32 50e-3
#define We2 (Rho2_re*Vre*Vre*Dre/Sig_re2)
#define We3 (Rho2_re*Vre*Vre*Dre/Sig_re3)
#define We32 (Rho2_re*Vre*Vre*Dre/Sig_re32)
#define SIGMA2 (2.*radius/We2)/2.
#define SIGMA3 (2.*radius/We3)/2.
#define SIGMA32 (2.*radius/We32)/2.
#define Lref (Dre/(2.*radius))   
#define gad (Lref*9.81/(Vre*Vre))

#define vrho2 (1.)
#define vrho3 (Rho3_re/Rho2_re)
#define vrho1 (Rho1_re/Rho2_re)
#define vmu2 (2.*radius/Re)
#define vmu3 (2.*radius/Re*(Mu_re3/Mu_re2))
#define vmu1 (2.*radius/Re*(1.81e-4/Mu_re2))

#define ycc 0.
#define hpool (4.*radius*2.)
#define hlame (0.707*radius*2.*1.)
#define xcc0 0.
#define xcc2 (xcc0 + hpool)
#define xcc3 (xcc0 + hpool + hlame)
#define hgout (0.5*radius*2.)//(0.25*radius*2.)
#define xcc (xcc3 + hgout + radius)

#define goutte(x,y) (sq(x - xcc) + sq(y - ycc))
#define brd(x,y) (sq(radius) - goutte(x,y))
#define init_goutte_speed(x,y) (sq(radius+(2*radius*0.5)/2)-goutte(x,y));
#define water(x) max(xcc2-x,0)
#define oil(x) max(xcc3-x,0)
#define pool(x) intersection(oil(x),water(x))
#define rfpool(x,y) (sq(radius-0.005)<goutte(x,y) && sq(radius+0.005)>goutte(x,y))

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);


int maxlevel = 8; // 12.;
int minlevel = 6;
double uem = 0.01;
double fem = 0.01;
double dtuser = 0.03;
double tmax = 4.; //8.; Use tmax=8 in order to see the complete evolution of the impact
face vector av[];


int main(int argc, char * argv[]){
  origin(0,0);
  size(4);
  init_grid(256);
  if (argc > 1)
    maxlevel = atoi (argv[1]);
    if (argc > 2)
    uem = atof (argv[2]);
    rho1 = vrho1; //gaz
    rho2 = vrho2; // eau 
    rho3 = vrho3; // oil
    mu1 = vmu1+(4*vmu1*exp(-t/(dtuser)));
    mu2 = vmu2;
    mu3 =vmu3;
    a=av;

  f2.sigma = (SIGMA2 + SIGMA32 - SIGMA3)/2.;
  f3.sigma = (SIGMA3 + SIGMA32 - SIGMA2)/2.;
  f1.sigma = (SIGMA3 + SIGMA2 - SIGMA32)/2.;
  run();
}

event init(t=0) {

    refine((rfpool(x,y)) && (level<maxlevel));


  scalar f4[],f5[],f6[],f7[],f8[];

  fraction(f4,brd(x,y));
  fraction(f5,pool(x));
  fraction(f6,oil(x));
  fraction(f7,water(x));
  fraction(f8,init_goutte_speed(x,y)); //init 


  foreach(){
    f2[]=f5[]+f4[];
    f3[]=f6[]-f7[];
    f1[]=clamp(1.-f2[]-f3[],0.,1.);
    
      #ifdef FILTERED
        rhov[]=rho(f2[],f3[]);   
      #else
        rhov[]=rho(f2[],f3[]);
      #endif
    }
    
    
    boundary ((scalar *){f2,f3,rhov}); 

    foreach() {
      if (f8[]>0){
          u.x[]=-1*f8[];

      }
  }
    boundary({u.x});


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
        adapt_wavelet ({u,f2,f3}, (double[]){uem,uem,fem,fem}, maxlevel,minlevel);
}




/*event movie (t =0;  t+=dtuser; t<= tmax) {
  clear();
  view (fov=15, quat = {0,0,-0.707,0.707}, tx =0.0, ty =-0.5, width = 2000, height = 2000);
  draw_vof ("f1",lw = 2);
  draw_vof ("f3",filled=1, fc = {0.7,0.7,0.} ,lw = 2);
  draw_vof ("f2",filled=1, fc = {0.,0.7,0.7},lw = 2);
  mirror({0,-1}) {
  draw_vof ("f1",lw = 2);
  draw_vof ("f3",lw = 2);
  draw_vof ("f2",lw = 2);  
  cells();
  }
  save("splash3p_9.mp4");

}*/


/**
The following movie has been done in level 13, I recommand to use at least a level 11 to be accurate
*/
/**
<center>
<video width="1024" height="1024" controls>
<source src="splash3p.mp4" type="video/mp4">
</video>
*/
