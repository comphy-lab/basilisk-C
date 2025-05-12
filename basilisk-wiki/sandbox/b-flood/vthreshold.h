/**
# Velocity threshold for saint-venant
*/


double threshold = 16; // max value of u
scalar thresh[];

event threshold_vel( i++ ){
  foreach(){	  
    if( norm(u) > threshold ){
    	thresh[] = 1;
    	foreach_dimension()
    		u.x[] /= norm(u)/threshold;		
     }
  }
}

event init( i = 0 )
  fprintf(stderr,"# velocity threshold activated , threshold = %lf\n Check that speed normalization remains exceptional thanks to raster thresh[]",threshold);

/**
## Link to the homepage
* [Homepage](Readme)
*/