/**
# Recursive initialization of volume fraction

This script allows to recursively refine around the interface to
increase the accuracy of the initialization providing a list of 
center and radius. This could be extended to other topology as 
done in the [built-in initialization algorithm](/src/fractions.h).
*/

/**
Chose an level of refinement. It means that it will refine 
at a level $Le = Le_{max} + Le_{init}$.
*/
int init_level = 7;

/** Recursive function which returns the volume fraction of cell after subdivising it into 8 finer subcells. */
double refine_frac(coord center, coord * pos, double * radius, int nb, double step, int level){
  // First check if the center of the cell lies in between a layer of Delta around the droplet radius
  int Rpd=0, Rmd=0;
  for (int j = 0; j < nb; j++){
      coord r;
      foreach_dimension() r.x = center.x-pos[j].x;
      #if dimension == 2
      Rpd +=  (sq(r.x) + sq(r.y) < sq(radius[j] + 0.5*sqrt(2)*step));
      Rmd +=  (sq(r.x) + sq(r.y) < sq(radius[j] - 0.5*sqrt(2)*step));
      #else
      Rpd +=  (sq(r.x) + sq(r.y) + sq(r.z) < sq(radius[j] + 0.5*sqrt(3)*step));
      Rmd +=  (sq(r.x) + sq(r.y) + sq(r.z) < sq(radius[j] - 0.5*sqrt(3)*step));
      #endif
  }
  // If the cell is in the inner radius, then the cell is full
  if (Rmd >= 1) return 1;
  // If the cell is not in the outer radius, then the cell is empty
  if (Rpd == 0) return 0;

  // Otherwise compute the fraction in the cell
  double floc = 0;
  // If the deepest level is reached, just take the value at the center of the cell
  if (level == init_level){ 
    for (int j = 0; j < nb; j++){
      coord r;
      foreach_dimension() r.x = center.x-pos[j].x;
      #if dimension == 2
      floc += (- sq(r.x) - sq(r.y) + sq(radius[j]) > 0);
      #else
      floc += (- sq(r.x) - sq(r.y) - sq(r.z) + sq(radius[j]) > 0);
      #endif
    }
    return floc;
  }
  // Othewise, subdivide into 8 subcells and recall this function
  else{ 
    coord subcenter;
    for (int icx = 0; icx < 2; icx++){
      subcenter.x = center.x + 0.5 * (icx - 0.5)*step;
      for (int icy = 0; icy < 2; icy++){
        subcenter.y = center.y + 0.5 * (icy - 0.5)*step;
        #if dimension == 2
        floc += refine_frac(subcenter, pos, radius, nb, 0.5*step, level+1);
        #else
        for (int icz = 0; icz < 2; icz++){
          subcenter.z = center.z + 0.5 * (icz - 0.5)*step;
          // Recursive call until level max or empty/full is reached
          floc += refine_frac(subcenter, pos, radius, nb, 0.5*step, level+1);
        }
        #endif
      }
    }
    return floc/pow(2,dimension);
  }
}