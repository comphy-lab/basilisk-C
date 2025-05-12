/**
# Equidistant Grid Diagnosis 
*/
int diaglevel = 6;
/**
Here is a function that should interpolate the equidistant-grid-diagosis points that lay within a coarse cell. For now it only returns the number of required interpolations. 
*/
int upsample_points(Point point){
  return  1 << ((diaglevel - point.level) * (dimension));
}

scalar s[];
int main(){
  X0 = Y0 = -L0/2.;
  /**
     There are cells at level 5, 6 and 7 and we aim to diagnose the
     solution at a resolution corresponding to level 6
  */
  init_grid(1 << 5);
  refine (x + y < -0.5 && level < 6);
  refine (sq(x - 0.25) + sq(y - 0.25) < sq(0.25) && level < 7);
  
  boundary({s});
  int cells = 0;
  for (int l = 1; l <= diaglevel; l++){ // This is Vincent's idea
    foreach_level(l){
      /**
	 The tree cells may be coarser or equal in size with respect to the
	 resolution of diagnosis. We do not iterate over the finer ones.
      */
      if (l < diaglevel){
	if (is_leaf(cell)){ // A Coarse leaf
	  cells += upsample_points(point);
	  printf("%d points should be interpolated\n", upsample_points(point));
      	}
      }else{               // l == diaglevel, we can directly access s[]. 
	  printf("%g %g %g\n", x, y, s[]);
	  cells++;
      }
    }
  }
  printf("We should have diagnosed %d cells out off %d that were requested.\n", cells,
	 (1 << (diaglevel * (dimension)))); 
}