/**
## Periodic simulations of buoyant emulsions
*/

#define dimension 2
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "no-coalescence.h" 
#include "tag.h"

/** My own header files <br/>
algebra.h --> include operator for coord and tens struct <br/>
parameters.h --> Define all the parameters <br/>
MISC.h --> definition of miscellaneous functions <br/>
CA.h --> continuous average functions <br/>
Drop.h --> compute the properties of the droplets <br/>
PA.h --> compute the particle average closure terms <br/>
CW.h --> compute the average layer by layer and the nearest averaged data <br/>
*/ 

#include "parameters.h"
#include "algebra.h"
#include "print_header.h"
#include "MISC.h"

#define compute_output 1
#if compute_output
#define threshold_volume pow(Dx*3,dimension)

/** declaration of the main instance containing the statistics  and drops 
  properties */

#include "Drop.h"
#include "CA.h"
#include "PA.h"
#include "CW.h"
Drop Bub[Nb+Nbsup]; // Creation of a global Drop instance (defined in Drop.h)
Drop BubTmp[Nb+Nbsup]; 
// event track_bub(t=0; t <= TMAX; t += dtprint){
event track_bub(i=0; i += dtprint; t <= TMAX){
  int tagi = 0;
  int nbubble;
  int shift;
  for(scalar s in interfaces) {
    scalar tag_level[];
    foreach() tag_level[] = s[]> EPS;
    nbubble = tag(tag_level);
     // compute each drops properties 
    Compute_drops(s,tag_level,nbubble,tagi,Bub);     
    #if COAL
      // mark each neigboring scalar 
      neighboring_scalar(s,tag_level,nbubble,tagi,Bub); 
    #endif
    shift = 0;
    for(int j=0;j<nbubble;j++) { 
      if(Bub[tagi+j].vol < threshold_volume) {     //Fragment suppression
        foreach()                              
          if(tag_level[]==j+1+shift)
            s[]=0;
        fprintf(stdout,"remove frag tag %d vol %g\n",Bub[tagi+j].tag,Bub[tagi+j].vol);
        // remove the fragment from the list by replacing the others
        for(int k = j; k < nbubble-1; k++)
          Bub[tagi+k] = Bub[tagi+k + 1];

        nbubble -= 1;
        j       -= 1;
        shift   += 1;
      }
    }
    for(int k=0;k<nbubble;k++) {
      Bub[tagi].tag = tagi;// should be realtag 
      Bub[tagi].si = s.i;
      tagi++;
    }
  }

  for (int j = 0; j < tagi; j++) {
    Bub[j].tagmax = tagi;
    Bub[j].time = t;
  }
  
  // the first step has no sense has BubTmp isn't def  
  assign_tags_from_last_step(BubTmp,Bub);

  nearest_dist(BubTmp,Bub);


  #if COAL
  if(t>1 && !shift) coalescence(Bub,tagi); 
  #endif

  /** Now that we have all the information on the drops we can compute 
   * closure terms using continuous and particle averages */

  compute_mv(Bub,BubTmp);// compute drop change of momentum (force) 
  CA ca = calcul_CA();
  PA pa = calcul_PA(Bub);
  compute_alpha_layer();
  
  #if BI_DISPERSE
  calcul_PAs(Bub,i);
  #endif

  // print all the data in files 
  if(main_process()){
    print_Drop(Bub,i); 
    print_CA(ca,i);
    print_PA(pa,"fPA.csv",i);
  }
  
  memcpy(BubTmp, Bub, sizeof(Drop)*(Nb+Nbsup));// copy the arry

  // Compute eulerian nearest statistics 
  Eulerian_nearest_stat(Bub,i);
}

#endif // end compute outpute

/**Initializes the domain with zero velocity */
event init (t=0){
  if(restore(file = "../dump/dump-last")){
    restore_old_csv();
  }else{
    coord Positions[Nb];
    double Radii[Nb]; 
    for (int j = 0; j < Nb; j++){
      if(j % 2 == 0) Radii[j] = R1; 
      else Radii[j] = R2; 
    }

    generate_drops_pos(Positions,Radii);

    foreach(){
      // changes the max level value to obtain a better initialization of the vof fields
      f[] = refine_frac((coord) {x,y,z},Positions,Radii,Delta,1,1);  
      u.x[]=u.y[]=u.z[]=0.;
    }
    init_and_print_header();
  }
}

/** Overloads the acceleration event to apply the effect of gravity. 
 * Due to the periodic boundary conditions the acceleration needs to be reduced by
 *  $\frac{\rho_{av}}{\rho}g$ rho is not available directly as a face centered vector
 *  so the stagger of f[] is adjusted to calculate it */

event acceleration (i++) {
  double rhoav=avg_rho();
  face vector av = a;
  foreach_face(y)
    av.y[] -= (1-rhoav/((f[]+f[0,-1])/2*(rho1-rho2)+rho2))*G; 
  #if CORREC
  coord Uavg = avg_U();
  foreach_face()
    av.x[] -= Uavg.x/dt;
  #endif
}

/** Outputs videos of the velocity, acceleration, Pressure, and tag fields

![Vertical velocity](RS/uy.mp4)

![Vertical Acceleration](RS/P.mp4)

![Pressure](RS/a.mp4)

![VoF tracers](RS/inter.mp4)
*/

event movies(t=0; t<=MOVIES*1200; t+=MOVIES){
// event movies(t=0; t<=MOVIES*200; i+=MOVIES){
  #if dimension == 3
  view (fov=30,
       camera="front",
       width = 800, 
       height = 800);
  

  clear();
  for (scalar s in interfaces){
    draw_vof(s.name, lw = 2.);
  }
  squares(color = "u.y", n = {0,0,1}, alpha = -Ls/2);
  save("P.mp4");
  
  clear();
  view (fov=30,
       camera="iso",
       width = 800, 
       height = 800);
    for (scalar s in interfaces){
    draw_vof(s.name, lw = 2.);
  }
  squares(color = "u.y", n = {0,0,1} ,alpha = -Ls/2.);
  squares(color = "u.y", n = {0,1,0} ,alpha = -Ls/2.);
  squares(color = "u.y", n = {1,0,0} ,alpha = -Ls/2.);
  save("U.mp4");

  #elif dimension == 2
  box();
    for (scalar s in interfaces){
    draw_vof(s.name, fc = {0.13,0.47,0.77}, min = 0, max = 3, lw = 2.);
  }
  squares("u.y");
  save("uy.mp4");
  clear();

  box();
  for (scalar s in interfaces){
    draw_vof(s.name, fc = {0.13,0.47,0.77}, min = 0, max = 3, lw = 2.);
  }
  squares("p");
  save("P.mp4");
  clear();

  box();
  squares("a.y");
  save("a.mp4");
  clear();

  box();
  for (scalar s in interfaces){
    draw_vof(s.name, fc = {0.13,0.47,0.77}, min = 0, max = 3, lw = 2.);
  }
  squares("g.y");
  save("g.mp4");
  clear();

  box();
  int index = 0;
  float * colors[] = {(float[]){1,0,0},(float[]){0,0,1},(float[]){0,1,1},(float[]){0,1,1},(float[]){1,1,1}};
  for (scalar s in interfaces)
    draw_vof (s.name, fc = colors[index++], filled = 1, min = 0, max = 3, lw = 2.);
  save ("inter.mp4");
  clear();

  #endif
}

// LOI is the initial number of intefaces
// It is used to recover a simulation
#ifdef LOI
event defaults (i=0)
{
  length_of_interfaces = LOI;
}
#endif

/** The flow parameters are calculated based on the previously defined non-dimensional parameters. 
 * The tolerance is reduced to 1e-5 and the boundaries are set to be periodic*/
int main(int argc, char** argv){  
  // make the problem dimensionless 
  L0 = Ls [0];
  DT = HUGE [0];
  size(Ls);
  init_grid(1<<(LEVEL));
  origin(-Ls/2.,-Ls/2.,-Ls/2.);
  rho1      =rho_d;
  mu1       =mu_d;
  rho2      =rho_f;
  mu2       =mu_f;
  f.sigma   =sig;
  TOLERANCE =1e-5 [*];
  foreach_dimension()
    periodic(right);
  run();
}

event stop(t=TMAX){
  dump (file = "dump-last");
  int nvof = list_len(interfaces);
  FILE * para = fopen("number_of_vof.txt","w");
  fprintf(para,"%d",nvof);
  fclose(para);
  return 1;
}
