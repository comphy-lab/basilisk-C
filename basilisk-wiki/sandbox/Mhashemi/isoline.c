#include "grid/octree.h" 
#include "utils.h"
#include "fractions.h"
#include "view.h"
#include "iso3D.h"


#define d0 8e-4

double init=0.0049;
double end=0.0059;

char name[80];
char name2[80];
char name3[80];
char name4[80];

char image[80];
char image2[80];

int main(int argc, char * argv[]){

if (argc > 1)
    init = atof (argv[1]);
if (argc > 2)
    end = atof (argv[2]);

 
for (double t=init; t <= end; t+=2e-5){

    sprintf (name, "snapshot-%g", t);
    sprintf (name2, "reduced_size_dump-%g", t);
    sprintf (name3, "jet_column-new%g", t);
    sprintf (name4, "gfsfield-%g", t);


   init_grid(2);
   if (restore (file = name2)) {
   
      vector u[];
      scalar f[];
      f.prolongation = fraction_refine; //for volume fraction field

      scalar * interesting_scalars =  {f} ; 

      restore (name2, interesting_scalars);

  
      boundary(interesting_scalars);
      fprintf(fout, "restore is done for t-%g\n", t);

      double  z;
      double  z_d0;

      for (int i=0; i <=11; i+=1){
        
	z=d0*i*0.625;
	z_d0=i*0.625;
	
	if (i==11) {
            
	   z=d0*7;
           z_d0=7;

	} 
        //sprintf (image, "cross-section-%g-z/d0-%g.ppm", t, z_d0);

        sprintf (image, "cross-section-%g-%g.ppm", t, z_d0);
	sprintf (image2, "contour-%g-%g.ppm", t, z_d0);


        view (fov =  2.5, quat = {0,0.707107,0,0.707107}, tx = 0.0,
	ty = 0.0, bg = {1.0, 1.0, 1.0},   width = 800, height = 716, 
	samples = 1);
    
   
       /** isoline2 ("f", val = 0.0 ,
	      np = {1, 0, 0}, alpha = z, lc= {0.001,0,0}, lw = 2.5 ); */
        cross_section("f", np = {1, 0, 0}, alpha = z, lc= {0.001,0,0}, lw = 2.5 );
        save(image);
	clear();

	squares("f", n = {1, 0, 0}, alpha = z, linear= true);
	save(image2);
	clear();
     }


     delete (interesting_scalars);
     fprintf(fout, "******postprocess is done for t=%g\n", t);
    }
  }
  
}












