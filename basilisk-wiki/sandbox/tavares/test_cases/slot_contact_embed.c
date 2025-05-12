/**
# Interface through an embed converging-diverging 2D solid slot
This test case has been initially run by Lopes Herrera to validate the coupling between vof and embed. I've added the contact angle computation between the vof interface and the embed solid. This a graphical test case, there is no quantitave calculation for validations.
 */

//#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact_embed.h"

/**
Geometrical parameters of the slot. */

#define A 0.4
#define B 0.1

u.n[left] = dirichlet(0);
u.t[left] = dirichlet(0);
u.n[right] = dirichlet(0);
u.t[right] = dirichlet(0);


f[bottom] = 1.;
u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0.);

// u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);

int LEVEL = 6;
int angle = 90;

scalar tag[], kappa[], triple_cell[], f1[], Alpha, app_angle[];
vector ns[], nf[], nf1[]; 
double volume_vof_init;

int main()
{
  origin(-0.5,0.);

  init_grid (1 << 6);
  
  rho1 = 1.;
  mu1 = 0.1;


  /**
  The gas phase is less denser and viscous */
  
  rho2 = 1.;
  mu2 = mu1;
  
  f.sigma=0.5;

  display_control (angle, 0, 180, "Angle");

  run();
}

#define T 15.


event init (i = 0)
{

  vertex scalar phi[];
  foreach_vertex() {
    double cone1 = -(x -(y-A));
    double cone2 = (x +(y-A));
    double slot1 = -(x-B);
    double slot2 = (x+B);


    double temp1 = intersection(slot1, slot2);
    double temp2 = intersection(cone1, cone2);

    phi[] = max(temp1,temp2);
  }

  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, 0.68 - y);
}


#if 1
event contact(i++; t<=T){
  /**Identifier les cellules coupées et leurs voisins dans le solide*/
  triple_cell_dectection(angle, f, cs, ns, nf1, nf);
  /**Reconstruire la fraction volumique dans les ghost cell*/
  fraction_reconstruction(f, nf);
}
#endif


#if 0
event apparent_angle(i++;t<=T){
  foreach(){
    app_angle[] = 0.;
    if (triple_cell[]){
      //coord ntemp1 = interface_normal (point, f);
      foreach_dimension() app_angle[] += nf1.x[]*ns.x[] ;
      app_angle[]= (acos(app_angle[])*180)/pi;
    } 
  }
}
#endif