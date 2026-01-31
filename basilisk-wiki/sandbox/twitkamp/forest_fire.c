/**
# Forest fire model 

The forest fire model is a simple cellular automaton that mimics the spreading of fire in a (dense) forest. Four simple rules are executed simultaneously. 


 1. Burning cell turns empty cell 
 
 2. Tree burns if any neighbour is burning
 
 3. Tree ignites with prob $f$ even if no neighbor is burning
 
 4. empty space becomes tree with prob $p$ 

The ratio of $f$ and $p$ dictates the dynamics of the system. 
*/

#include "run.h"
#include "utils.h"
#include "output.h"
#include "view.h"

double Niter, f, p;
int main(){ 
  L0 = 1.;
  N = 128;
  periodic(top);
  periodic(right);

  Niter = 1e3;
  
  run();
}

scalar tree_field[], tree_field_step[];
event init(i = 0){
  foreach(){
    tree_field[] = 0.;
  }
  boundary({tree_field});
}

event movie(i = 0; i += 1., i <= Niter){
  char mov[80];
  sprintf(mov, "tree_%g_%g.mp4", f, p);
  //output_ppm(tree_field, file = mov, min = -1, max = 1, n = 2*N); 
  
  view(fov=20, tx = 0, ty = -0.5);
  //view(fov=20, quat = {0,0,-0.707,0.707}, tx = 0, ty = -0.5);
  draw_vof("tree_field");
  save("tree.mp4");
}


//three states 
// Tree: 1
// Empty: 0      
// Burning: -1
event tree_flip(i++){
  f = 0.001;
  p = 0.01;

  foreach(){
    double prob_p = (noise() + 1.)/2.;
    double prob_f = (noise() + 1.)/2.;
    
    /** # Rule 1 */
    if (tree_field[] == -1.){
      tree_field_step[] = 0.;
    }
    
    /** # Rule 2 and 3 */
    else if (tree_field[] == 1. && 
             (tree_field[0, 1] == -1. || tree_field[0, -1] == -1. || tree_field[-1, 0] == -1. || tree_field[1, 0] == -1. || 
             prob_f < f)){
      tree_field_step[] = -1.;
    }
    
        
    /** # Rule 4 */
    else if (prob_p < p && tree_field[] == 0.){
      tree_field_step[] = 1.;
    }
  }
  foreach(){
    tree_field[] = tree_field_step[];
  }
}

/**
![Forest fire spreading for $f$ = 0.001 and $p$ = 0.01.](forest_fire/tree.mp4)(width="60%")
*/