/**
  Simulation von versetzten, quadratischen Zylindern
  
  interne Variable
  mu = dyn. Viskosität
  rho = Dichte
  a = Beschleunigungvektor
  
  cm: skalare Größe, die die Fläche[2D]/ Volumen[3D] einer Zelle angibt (Anteil)
  fm: face vector, welcher die Länge[2D], Fläche[3D] einer Seite einer Zelle angibt
  Delta: Größe der Stencil-Zelle
  
  dv() = (cube(Delta)*cm[])
  
  Bei Embed gilt: cm[] = cs[] -> Innerhalb der Filterstruktur ist cm[]=0
 */


//~ #include "grid/octree.h"
#include "grid/quadtree.h"
#include "output_htg_pv590_mpi.h"

#include "embed.h"
#include "navier-stokes/centered.h"

#include "navier-stokes/perfs.h"
#include "maxruntime.h"

#define PHI 0.65


/** default level */
int LEVEL;

/** Distance Field */
scalar d[];
face vector muc[];

/** Main function */
int main(int argc, char * argv[])
{  
	maxruntime (&argc, argv);
	if (argc == 2) {
		LEVEL = atoi(argv[1]);
	} else {
		if(pid()==0){		
			printf("Specify Level\n");
			printf("Usage: mpirun -np X ./turbulent_flow... LEVEL\n");
		}
		return 1;
	}
			
	system("mkdir -p htg");
	
	init_grid (1<<(LEVEL));
	size (1e-2);
	origin(-L0/2.,-L0/2.,-L0/2.);
	
	/**
	* Periodische Randbedingung in allen Raumrichtungen
	*/

	foreach_dimension()
		periodic(right);
	
	const face vector g[] = {1.,0.,0.};
	a = g;	
	
	mu = muc;
	const scalar rhov[] = 1.;
	rho = rhov;
	
	//~ NITERMAX=250;
	TOLERANCE=1e-4; // Default Tolerance: 1e-3
	/** Limit Timstep */
	DT=1e-1; 
	run();
}

/**
 * Set BC on the filter surface */
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
#if dimension==3
u.r[embed] = dirichlet(0.);
#endif

scalar un[];

/** Initialiesierungszeitschritt */ 
event init (t = 0) {	
	/** Calculate Distancefield on vertices (average)*/
	vertex scalar phi[];
	
	double H = L0/2.;
	double D = sqrt(1-PHI)*H;
			
	//~ double eps = D/10.;
	//~ refine(	(fabs(x)<(D/2.+eps)) && (fabs(x)>(D/2.-eps)) &&  \
			//~ (fabs(y)<(H/2.-D/2.+eps)) && (fabs(y)>(H/2.-D/2.-eps)) && \
			//~ (fabs(y)<(H/2.+D/2.+eps)) && (fabs(y)>(H/2.+D/2.-eps)) && \
			//~ level < LEVEL); 
	
	foreach_vertex(){
		// 2 in the middle			
		phi[] = 					-2*( (fabs(x)<(D/2.)) && (fabs(y)>(H/2.-D/2.)) && (fabs(y)<(H/2.+D/2.)) ) + 1;
		
		// 4 corners
		phi[] = intersection(phi[],	-2*( (fabs(x)>(H-D/2.)) && (fabs(y)>(H-D/2.)) ) + 1);
				
		// left and right
		phi[] = intersection(phi[],	-2*( (fabs(x)>(H-D/2.)) && (fabs(y)<(D/2.)) ) + 1);
		    
	}
	boundary ({phi});


	/** Calucalte Fractions of the embed Geometry) */
	fractions (phi,cs,fs);	 		
	fractions_cleanup(cs,fs);
	boundary({cs,fs});
	restriction({cs,fs});
	
	foreach_face()
		muc.x[]=fm.x[]*1e-6;
	boundary ((scalar *){muc});
	
	//~ foreach()
		//~ u.x[]= cs[] >= 1.?0.01:0.;
	boundary((scalar*){u});
	foreach()
		un[] = u.x[];

}

/**
 * Laufzeitstatistiken
 */
event logfile (i++){
	//~ double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
	double du = 0.;
	fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
	fprintf (stderr, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
			LEVEL, i,
			mgp.i, mgp.nrelax, mgp.minlevel,
			mgu.i, mgu.nrelax, mgu.minlevel,
			du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
	fprintf(stderr, "dt: %g\n", dt);
	/** Stop the simulation */
	//~ if (i > 100 && (avg < 1e-9 || du < 0.0001)) {
		//~ if (pid()==0) printf("Simulation ist Stationär!");
		//~ return 1;
	//~ }
}

event disturb(t=0.1){
	foreach(){
		u.x[]+=cs[]*un[]*(0.1*noise());
		u.y[]+=cs[]*un[]*(0.1*noise());
	}
	boundary((scalar*){u});
}

event snapshot (t+=0.005)
//~ event snapshot (i+=1)
{	
	boundary(all);
	restriction(all);
	
	scalar mpi_color[], l_ref[];
	foreach_cell(){
		l_ref[] = level;
		mpi_color[] = pid();
	  }
	 	  
	char path[]="htg/";
	char prefix[80];
	sprintf(prefix, "data_%06d", i);  
	output_htg((scalar *) {cs, d,p, mpi_color},(vector *){u}, path, prefix, i, t);
	
}


event ending(t=60)
	printf("Simulation beendet");

/**
 * AMR-Refinement
 * 
 * s.nc: Anzahl der vergröberten Zellen
 * s.nf: Anzahl der verfeinerten Zellen
 * 
 * Codebaustein aus http://basilisk.fr/sandbox/Antoonvh/bootstrapbubble.c
 *   scalar ff[];
 *   foreach()
 *   	ff[] = filtergeom[];
 *   boundary({ff});
 *   adapt_wavelet ({filtergeom,ff,u},...);
 * forciert gleiches RF entlang der Filtergeometrie
 */

event adapt (i++) {
  scalar ff[];
  foreach()
    ff[] = cs[];
  boundary({ff});
  astats s = adapt_wavelet ({cs,ff,u}, (double[]){1e-2,1e-1,1e-3,1e-3,1e-3},LEVEL, 4);
  
  //~ astats s = adapt_wavelet ({cs,u}, (double[]){1e-6,1e-3,1e-3},LEVEL, 5);
  //~ astats s = adapt_wavelet ({cs}, (double[]){1e-4},LEVEL, 5);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

