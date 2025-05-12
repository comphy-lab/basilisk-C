/**
Simple small test case to checking vtp/htg/vti export.

*/
/** To test HTG Export */
//~ #define HTG 1
//~ #include "grid/octree.h"
//~ #include "grid/quadtree.h"
//~ #include "../../output_htg.h"


/** To test VTI Export */
#define VTI 1
#include "grid/multigrid3D.h"
//~ #include "grid/multigrid.h"
#include "output_vti_pv590_mpi.h"

/** To test VTP Export */
#define VTP 1
// needed, because there is no flow
#include "navier-stokes/centered.h"
// use "new_tracer_particles (0);" in init()

#include "tracer-particles.h"
// use "new_inertial_particles (0);" in init()    
#include "stokes-particles.h" 
#include "output_vtp_pv590.h"

Particles flow, heavy;
int np = 100;
double ti = 0.;

/** 
# Test
*/	


int main() {

	#if HTG
	system("mkdir -p htg");
	#endif
	#if VTI 
	system("mkdir -p vti");
	#endif
	#if VTP
	system("mkdir -p vtp");
	#endif
  
	L0 = 2;
	#if dimension == 3
	X0 = Y0 = Z0 = -L0/4.;
	#elif dimension == 2
	Z0 = Y0 = -L0/4.;
	#endif
	
	
	
	#if VTI || VTP
		init_grid (16); 
	
	#elif HTG
		init_grid (1<<1); // Initialize a 2 x 2 grid
			#if dimension == 2
				refine   ((x >= 12) && (y >= 12) && (level <  3)); // Refine to top right corner  
				unrefine ((x <  8 ) && (y <  8 ) && (level >= 1)); // Coarsen the bottom left corner
			#elif dimension == 3
				refine   ((x >= 12) && (y >= 12) && (z >= 12) && (level <  3)); // Refine to top right corner  
				unrefine ((x <   8) && (y <   8) && (z <   8) && (level >= 1)); // Coarsen the bottom left corner
			#endif	
	#endif

	
	#if VTP
	if (pid()==0){
		int j=0;			  

		flow = new_tracer_particles (np);
		heavy= new_inertial_particles (np);
		 while (j < np){
			  double xp=noise()*L0/2.;
			  double yp=noise()*L0/2.;			
			  #if dimension > 2
			  double zp=noise()*L0/2.;			
			  #endif
			  
				pl[flow][j].x  = xp;
				pl[flow][j].y  = yp;
				#if dimension > 2
				pl[flow][j].z  = zp;
				#endif
				xp=noise()*L0/2.;
				yp=noise()*L0/2.;	
				pl[heavy][j].x = xp;
				pl[heavy][j].y = yp;
				
				#if dimension > 2
				zp=noise()*L0/2.;		
				pl[heavy][j].z = zp;
				#endif
				pl[heavy][j].u2.x = 2000; // Density
				pl[heavy][j].u2.y = 0.05e-3 + 0.01e-3*noise(); // Radius
				//~ pl[heavy][j].u2.z = 1; // Timescale taup
				//~ printf("Partikel platziert:%i\n", j);
				j++;
			
		}

	  } else {
			flow =  new_tracer_particles (0);
			heavy = new_inertial_particles (0);
	  }
	  particle_boundary (flow);  
	  particle_boundary (heavy);  
	  set_particle_attributes (flow);   	  
	  set_particle_attributes (heavy);   
	#endif
	#if HTG || VTP || VTI
	
	scalar d[],mpi_color[];
	foreach(){
		d[]=x*x + y*y + z*z;
		mpi_color[]=pid();
	}
	
	boundary({d,mpi_color});
	scalar l_ref[];
	vector vectt[];
	foreach(){
		l_ref[] = level;
		
		//~ vectt.x[] = x;
		//~ vectt.y[] = y*y;
		//~ vectt.z[] = z*z*z;
		vectt.x[] = 1;
		vectt.y[] = 2;
		vectt.z[] = 3;
	}
	boundary(all);
	#endif


	int i = 0;
	double t = 0.;
	#if HTG
	{	  
		char path[]="htg/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_htg((scalar *) {l_ref,d,mpi_color},(vector *){vectt}, path, prefix, i, t);
	}
	#endif
	#if VTI
	{
		char path[]="vti/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_vti((scalar *) {l_ref,d,mpi_color},(vector *){vectt}, path, prefix, i, t);
	}
	
	#endif
	#if VTP
	{
		char path[]="vtp/";
		char prefix[80];
		sprintf(prefix, "flow_%06d", i); // underscore is necessary! 
		output_vtp(flow, false, path, prefix, i, t, ti);
		sprintf(prefix, "heavy_%06d", i);  
		output_vtp(heavy, true, path, prefix, i, t, ti);
	}
	#endif
}
