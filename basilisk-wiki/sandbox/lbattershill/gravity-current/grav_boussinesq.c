#include "grid/quadtree.h"
#include "embed.h"
#include "convection_boussinesq_grav.h"
#include "output_vtu_foreach.h"
#include "utils.h"

#define MINLEVEL 4
#define maxlevel 13
#define gr 6e6
#define sc 1

double ue = 1e-3;


int main() {
	  L0 = 32.;
	  origin (0, 0);
	  Gr = gr; Sc = sc; 
	  DT = 0.001;
	  TOLERANCE = 1e-3;
	  run(); 
}


u.n[top] = neumann (0);
u.t[top] = neumann (0);
uf.n[top] = neumann(0);
u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
uf.n[bottom] = dirichlet(0);
u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
uf.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);
u.t[right] = dirichlet(0);
uf.n[right] = dirichlet(0);
u.n[embed] =  dirichlet(0.);
u.t[embed] =  dirichlet(0.);
uf.n[embed] = dirichlet(0);

event init (t = 0)
{
	//refine(y > 0.8 && y < 1.2 && x > 9 && x < 11 && level < maxlevel);
	  vertex scalar phi[];
	    foreach_vertex() {
		        phi[] = y < 1;
			      }
	  boundary ({phi});
	  fractions (phi, cs, fs);

	  foreach() {
		  if (x <= 6.667 && y >= 0.)
			  f[] = 1;
		  else
			  f[] = 0;
		  foreach_dimension()
			  u.x[] = 0;
		}
	  boundary ({f,u});

}

#include "view.h"
event logfile (i++)
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);


//Output .txt files for every timestep (x, y, f[])

int time_step = 0.;
double shift_control = 1.;
event output_files  (t += 0.03; t <35.){ 
		
		        char one[80], two[80], three[80], four[80];;
			 
			//First two outputs show JPG snapshots for FOV 
			sprintf (one, "13_PNG_H_150_linear_false_fov_Gr_6e6_Sc_1-%d.png",time_step );
			//sprintf(two, "TXT_linear_false_fov_Gr_9.9e7_Sc_700-%d.txt", time_step);
			output_ppm(f, file = one, min = 0, max = 1, linear = true, n = 1280, box = {{11.67,0},{15,1}});
			

			//The next two outputs are for a shiftting box
  
			/*
			if (time_step > 3500) {
				sprintf(three, "PNG_linear_false_full_Gr_1e8_Sc_10-%d.png", time_step);
				double x_shift = 8.75- 0.00242*shift_control*0.45;
				double y_shift = 15.75- 0.00242*shift_control*0.45;
				
				output_ppm(f, file = three, min = 0, max = 1, linear = false, n = 1792, box = {{x_shift,0},{y_shift,1}});
				sprintf(four, "PNG_linear_true_full_Gr_1e8_Sc_10-%d.png", time_step);
				output_ppm(f, file = four, min = 0., max = 1., linear = true, n = 1792, box = {{x_shift,0},{y_shift,1}});

				
				shift_control += 1.;

			}	
		
			  FILE * fp = fopen (s3, "w");
			  foreach()
				if (x > 8.75 && x < 11.25 && y < 1.)    
			    		fprintf (fp,"%g %g %g\n",x,y,f[]);
			    fclose (fp);
 			*/
		 	   //output_ppm(f, file = s2, box = {{8.75,0},{11.25,1}});
			 time_step += 1.;
}






event adapt (i++)
{
	scalar omega[];
	vorticity(u,omega);
	adapt_wavelet((scalar *){cs,omega,f},(double[]){ue,ue,ue,0.001},maxlevel);
				        
}

/** 2D time averaged & compared against experiments: */
/**
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/16/Time_avg.svg/1280px-Time_avg.svg.png" width=750px\>
*/