/**
Huge testcase for checking HTG/VTP/VTI-Export.
This test is based on "sphere.c"
*/
/** Test *vti? */
//~ #include "grid/multigrid.h"
//~ #include "grid/multigrid3D.h"
//~ #include "output_vti_pv590_mpi.h"

/** Test *.htg? */
//~ #include "grid/quadtree.h"	
#include "grid/octree.h"
#include "output_htg.h"

#include "embed.h"
#include "navier-stokes/centered.h"

/** Test *.vtp? */
//~ // use "new_tracer_particles (0);" in init()
#include "tracer-particles.h"
//~ // use "new_inertial_particles (0);" in init()    
#include "stokes-particles.h" 
#include "output_vtp_pv590.h"

#include "navier-stokes/perfs.h"

#if dimension==3
#include "lambda2.h"
#endif

int maxlevel = 6;
face vector muv[];

Particles flow, heavy;
int np = 100; // number particle
double ti = 0.5; // time inject

/**
# main() */
int main()
{
	#if TREE
		system("mkdir -p htg");

	#else
		system("mkdir -p vti");
	#endif
	
		system("mkdir -p vtp");
	
	
	periodic(bottom);
	#if dimension==3
		periodic(front);
	#endif


	#if MULTIGRID
	init_grid (1<<maxlevel);  
	#else
	init_grid (1<<5);
	#endif


	size (4.);
	origin (-L0/2, -L0/2., -L0/2.);

	mu = muv;
	run();
}

/**
The viscosity is just $1/Re$, because we chose a sphere of diameter
unity and an unit inflow velocity. */

event properties (i++)
{
	foreach_face()
		muv.x[] = fm.x[]/300.;
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
#if dimension==3
u.r[embed] = dirichlet(0.);
#endif
event init (t = 0) {

	/**
	We initially refine only in a sphere, slightly larger than the solid
	sphere. */
  
	#if TREE
		#if dimension == 2
			refine ( sq(x) + sq(y)  < sq(0.6) && level < maxlevel);
		#elif dimension == 3
			refine (sq(x-1.) + sq(y-1.) + sq(z-1.) < sq(0.6) && level < maxlevel);
		#endif
  #endif
  /**
  We define the unit sphere. */

	vertex scalar phi[];
	foreach_vertex(){
		#if dimension == 2
			phi[] = sq(x) + sq(y) - sq(0.5);  
		#elif dimension == 3
			phi[] = x*x + y*y + z*z - sq(0.5);
		#endif
	}
	boundary ({phi});
	fractions (phi, cs, fs);
  

	/**
	We set the initially horizontal velocity to unity everywhere
	(outside the sphere). */

	foreach(){
		u.x[] = cs[];
	}
	
	new_tracer_particles (0);
	new_inertial_particles (0);
	
}

event add_particle(t=ti){
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
}


event report(i+=10)
{
	scalar mpi_color[];
	foreach(){
		mpi_color[] = pid();
	}
	boundary(all);
	restriction (all);
	
	
  	#if TREE
	{	  
		char path[]="htg/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_htg((scalar *) {p,cs,mpi_color},(vector *){u}, path, prefix, i, t);
	}
	
	#else
	{
		char path[]="vti/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_vti((scalar *) {p,cs,mpi_color},(vector *){u}, path, prefix, i, t);
	}
	
	#endif

	if(t>=ti){
		char path[]="vtp/";
		char prefix[80];
		sprintf(prefix, "flow_%06d", i); // underscore is necessary! 
		output_vtp(flow, false, path, prefix, i, t, ti);
		sprintf(prefix, "heavy_%06d", i);  
		output_vtp(heavy, true, path, prefix, i, t, ti);
	}

}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */
event logfile (i++)
	fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);


#if TREE
event adapt (i++) {
	astats s = adapt_wavelet ({cs,u}, (double[]){1e-5,2e-3,2e-3,2e-3},
				maxlevel, 4);
	fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif

event ending(t=60){
}

