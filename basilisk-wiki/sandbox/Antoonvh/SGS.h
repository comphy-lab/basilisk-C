/** 
# Apply the Vreman Sub-Grid Scale Model
This header file can be used to run with a sub-grid scale model. By default, the eddy viscosity is based on the [Vreman Eddy viscosity model](vreman.h). This file should included after the Navier-Stokes/centered.h module is included but before grid adaptation is carried out. This module also applies diffusion to the list of tracers.  
 */
#include "tracer.h"
#include "diffusion.h"
#include "vreman.h"
/**
We initialize global variables: $K_h$ and $K_m$ are the diffusivity for heat and momentum, respectively. *Pr* is the turbulent Prandlt number, defined as $Pr=\frac{K_h-\nu}{K_m-\nu}$, with $\nu$ being the molecular viscosity (*molvis*). *Csmag* is the classical Smagorinsky constant and *tracers* is a list of diffusive flow tracers, e.g. the Buoyancy field.     
*/
face vector Km[],Kh[];
(const) face vector Pr;
scalar Evis[]; //Cell Centered diffusivity
double molvis;  
double Csmag;
scalar * tracers;
/**
Since $K_m\propto K_h \propto \Delta^2$, proper care is required for evaluating corresponding diffusivities at the different levels of resolution boundaries. 
*/

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
/**
We set some default values and for the parameters that may be overwritten by the users' *init* event. 
*/

event defaults(i=0){
  if (dimension!=3) //Allow to run, but give a warning
    fprintf(stdout,"Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n",dimension);
  mu=Kh;
  Pr=unityf;
  molvis=0.;
  Csmag=0.12;
  Evis.prolongation=Evisprol;
  Evis.restriction=Evisres;
/** On tree grids we do not care about the diffusivities on refined and/or coarsend cells and faces. These should be properly reevaluated before they apear in any computation*/
#if TREE  
  Evis.refine=no_restriction;
  Evis.coarsen=no_restriction;
  foreach_dimension(){
    Kh.x.refine=no_restriction;
    Km.x.coarsen=no_restriction;
  }
#endif
}
/**
We calculate the centered eddyviscosity and translate it into the required face-field diffusivity. 
*/
event Eddyvis(i++){
  eddyviscosity(Csmag,u,molvis,Evis); 
  boundary({Evis});
  foreach_face(){
    Km.x[]=(Evis[]+Evis[-1])/2;
    Kh.x[]=(Pr.x[]*(Km.x[]-molvis))+molvis;
  }
  /** In 3D, there are 4 finer faces per coarser face. So consistency with $K_m,K_h\propto \Delta^2$ is conviniently achieved by applying the *Boundary_flux* function for these face fields */ 
  boundary_flux({Km,Kh}); 
}
/**
   For some user-friendliness, the SGS-mixing of the *tracers* is hard-coded into this header file. 
*/

event tracer_diffusion(i++){
  for (scalar s in tracers)
    diffusion(s,dt,Kh);
}
