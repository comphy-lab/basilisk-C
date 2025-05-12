/* This header file can be used to run in Large Eddy Simulation (LES) mode. 
By default, the eddy viscosity is based on the Vreman Eddy viscosity model. 
In the .c-script, this file should be included after the Navier-Stokes/centered.h 
module is included but either before all other events or after the grid 
adaptation event.
 */

#include "tracer.h"
#include "diffusion.h"
//#include "Vreman.h"
#include "Dynamic_Smagorinsky.h"


/* Some global variables are initialized, Kh and Km are the diffusivity 
for heat and momentum, respectively, Pr is the turbulent Prandlt number, 
defined as Pr=Kh−νKm−ν, with ν being the molecular viscosity (molvis), 
Csmag is the classical Smagorinsky constant and tracers is a list of 
diffusive flow tracers, e.g. the Buoyancy field and total water vapor 
field. */


face vector mu_t[]; // turbulent viscosity
scalar Evis[]; // Cell Centered diffusivity
double Csmag;
scalar * tracers;

/* Since Km∝Kh∝Δ2, proper care is required for evaluating corresponding 
   diffusivities at the different levels of resolution boundaries. */

static inline void Evisprol(Point point,scalar s){
  foreach_child()
    Evis[]=bilinear(point,Evis)/4.; 
}

static inline void Evisres(Point point,scalar s){
  double sum = 0.;
  foreach_child()
    sum += s[];
  s[] = sum/2.;
}

/* We set some default values for the parameters that may be overwritten 
by the users’ init() event. */

event defaults(i=0){
  if (dimension!=3) //Allow to run, but give a warning
    fprintf(stdout,"Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n",dimension);
  Csmag=0.12;
  Evis.prolongation=Evisprol;
  Evis.restriction=Evisres;

/* On tree grids we do not directly care about the diffusivities on 
refined and/or coarsend cells and faces. These should be properly 
reevaluated before they appear in any computation */

#if TREE  
  Evis.refine=no_restriction;
  Evis.coarsen=no_restriction;
  foreach_dimension(){
    mu_t.x.coarsen=no_restriction;
  }
#endif
}

/* The centered eddyviscosity is evaluated and translated into the 
mandadory face-field diffusivity for the usage in other parts of the 
solver. */


//event Eddyvis(i++){
//  eddyviscosity(Csmag,u,Evis); 
//  boundary({Evis});
//  foreach_face(){
//    Km.x[]=(Evis[]+Evis[-1])/2;  // Face center
//  }
  

  /* In 3D, there are 4 finer faces per coarser face. So consistency 
     with Km,Kh∝Δ2 is conviniently achieved by applying the Boundary_
     flux() function for these face fields */

//  boundary_flux({Km}); 
//}


// overload viscosity event
event viscous_term(i++){

  correction (dt); // better approximation of solenoidal velocity field
  // compute eddy viscosity
  //eddyviscosity(Csmag,u,Evis);
  eddyviscosity(u,Evis); // Dynamic Smagorinsky
  boundary({Evis});
  foreach_face() {
    mu_t.x[] = 0.5*(rho[]*Evis[]+rho[-1]*Evis[-1]);  // Face center.  Multiplied by rho to convert to dynamic viscosity
    }
  correction (-dt);
  boundary_flux({mu_t}); 

  // add to dynamic viscosity
  face vector muv = mu;

 
  foreach_face()
    muv.x[] = mu.x[] + mu_t.x[];
  }


