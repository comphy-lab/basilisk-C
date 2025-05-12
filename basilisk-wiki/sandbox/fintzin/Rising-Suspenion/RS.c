/**
## Tri/Bi-periodic simulations of rising droplets aiming to close the Navier-Stokes averaged equations.

Warning : Before running this file apply this 
[patch](http://basilisk.fr/sandbox/fintzin/Patches/einstein_sum2.patch)
with darcs apply file.patch on your basilisk installation. 
*/

#define dimension 3
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "output.h"
#include "no-coalescence.h" 
#include "tag.h"



/** My own header files <br/>
algebra.h --> include operator for coord and tens struct <br/>
parameters.h --> Define all the parameters <br/>
MISC.h --> definition of miscellaneous functions <br/>
CA.h --> continuous average functions <br/>
Drop.h --> compute the properties of tehs droplets <br/>
PA.h --> compute the particular average closure terms <br/>
CW.h --> compute the continuity waves data <br/>
*/ 

#include "parameters.h"
#include "algebra.h"
#include "MISC.h"

#define compute_output 1
#if compute_output
/** declaration of the main instance containing the statistics  and drops properties */
#include "Drop.h"
#include "CA.h"
#include "PA.h"
#include "CW.h"
Drop Bub[Nb+Nbsup]; // Creation of a global Drop instance (defined in Drop.h)
Drop BubTmp[Nb+Nbsup]; 
event track_bub(t=0; t <= TMAX; t += dtprint){
  int tagi = 0;
  int nbubble;
  int shift;
  for(scalar s in interfaces) {
    scalar tag_level[];
    foreach() tag_level[] = s[]> EPS;
    nbubble = tag(tag_level);
    Compute_drops(s,tag_level,nbubble,tagi,Bub);      // compute each drops properties 
    #if COAL
      neighboring_scalar(s,tag_level,nbubble,tagi,Bub); // mark each neigboring scalar 
    #endif
    shift = 0;
    for(int j=0;j<nbubble;j++) { 
      if(Bub[tagi+j].vol < pow(Dx*3,dimension)) {     //Fragment suppression
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

  Bub[0].tagmax = tagi;
  if(t==FIRSTSTEP) memcpy(BubTmp, Bub, sizeof(Drop)*(Nb+Nbsup));// copy the arry
  
  assign_tags_from_last_step(BubTmp,Bub,tagi);

  if(main_process())
    print_pair_dist(BubTmp,Bub,tagi,i,t);


  #if COAL
  if(t>1 && !shift) coalescence(Bub,tagi); 
  #endif

  /** Now that we have all the information on the drops we can compute 
   * closure terms using continuous and particular averages */
  CA ca;   // Creation of a global CA instance (defined in CA.h)
  PA pa;   // Creation of a global PA instance (defined in PA.h)

  calcul_CA(&ca);
  calcul_PA(&pa,Bub,BubTmp);
  compute_alpha_layer();

  // save the datas
  if(main_process()){
    print_Drop(BubTmp,Bub,tagi,i,t); 
    print_CA(ca,i);
    print_PA(pa,i);
  }

  dump (file = "dump-last");
  int nvof = list_len(interfaces);
  FILE * para = fopen("number_of_vof.txt","w");
  fprintf(para,"%d",nvof);
  fclose(para);

  memcpy(BubTmp, Bub, sizeof(Drop)*(Nb+Nbsup));// copy the arry
}

#endif // end compute outpute

/**Initializes the domain with zero velocity and outputs the flow parameters*/
event init (t=0){
  if(restore(file = "../dump/dump-last")){
    restore_old_csv();
  }else{
    foreach(){
      f[] = geometry(x,y,z);  
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
    av.y[] -= (1-rhoav/((f[]+f[0,-1])/2*(rho1-rho2)+rho2))*g; 
  #if CORREC
  coord Uavg = avg_U();
  foreach_face()
    av.x[] -= Uavg.x/dt;
  #endif
}

/** Outputs videos of the velocity, Pressure, vorticity, and tag fields

![Pressure face view](RS/P.mp4)

![Velocity iso view ](RS/U.mp4)
*/

event movies(t=0; t<=TMAX; t+=MOVIES){
  #if dimension == 3
  view (fov=30,
       camera="front",
       width = 800, 
       height = 800);
  clear();

  box();
  draw_vof("f", lw = 2.);
  squares(color = "u.y", n = {0,0,1}, alpha = -Ls/2);
  save("P.mp4");
  
  clear();
  view (fov=30,
       camera="iso",
       width = 800, 
       height = 800);
  box();
  draw_vof("f", lw = 2.);
  squares(color = "u.y", n = {0,0,1} ,alpha = -Ls/2.);
  squares(color = "u.y", n = {0,1,0} ,alpha = -Ls/2.);
  squares(color = "u.y", n = {1,0,0} ,alpha = -Ls/2.);
  save("U.mp4");

  #elif dimension == 2
  box();
  draw_vof("f", lw = 2.);
  squares("u.y");
  save("uy.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("p");
  save("P.mp4");
  clear();

  box();
  squares("a.y");
  save("a.mp4");
  clear();

  scalar inter[];
  for(scalar s in interfaces)
    foreach()
      inter[] += (s.i)* (s[]>EPS);

  box();
  squares("inter",map=randomap);
  save("inter.mp4");
  clear();
  #endif
}


#ifdef LOI
event defaults (i=0)
{
  length_of_interfaces = LOI;
}
#endif

/** The flow parameters are calculated based on the previously defined non-dimensional parameters. 
 * The tolerance is reduced to 1e-4 and the boundaries are set to be periodic*/
int main(int argc, char** argv){  
  size(Ls);
  init_grid(1<<(LEVEL));
  origin(-Ls/2.,-Ls/2.,-Ls/2.);
  rho1      =r1;
  mu1       =mu_d;
  rho2      =rho_f;
  mu2       =mu_f;
  f.sigma   =sig;
  TOLERANCE =1e-5;
  foreach_dimension()
    periodic(right);
  run();
}

event stop(t=TMAX){
  return 1;
}

/**
## References

~~~bib
@article{batchelor1970stress,
  title={The stress system in a suspension of force-free particles},
  author={Batchelor, GK},
  journal={Journal of fluid mechanics},
  volume={41},
  number={3},
  pages={545--570},
  year={1970},
  publisher={Cambridge University Press}
}
~~~
*/