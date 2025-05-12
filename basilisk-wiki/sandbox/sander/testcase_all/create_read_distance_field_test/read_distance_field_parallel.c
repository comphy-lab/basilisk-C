/**
This program reads the binary distance field and outputs a paraview-
compatible file.
*/


//~ #define AMR
//~ #include "grid/octree.h"
//~ #include "output_htg.h"

#define MG
#include "grid/multigrid3D.h"
#include "output_vti_pv590_mpi.h"

#include "embed.h"
//~ #include "navier-stokes/centered.h"

#include "readxyz.h"

int LEVEL = 7;


int main(int argc, char * argv[])
{  
	if (argc > 1)
		LEVEL = atoi(argv[1]);

	#ifdef AMR
		system("mkdir -p htg");
		printf("Prozess:%04i: Erzeuge FractionField aus DistanceField für AMR auf Level %i\n", pid(),LEVEL);
	#endif
	#ifdef MG
		system("mkdir -p vti");
		printf("Prozess:%04i: Erzeuge FractionField aus DistanceField für MG auf Level %i\n",pid(), LEVEL);
	#endif
	init_grid (1<<(LEVEL));
	size (1.);

	foreach_dimension()
		periodic(right);

	FILE * fp ; 
	scalar d[];
	#ifdef MG
		fp = fopen("kelvin08-distance.raw.mg", "rb");
	#endif
	#ifdef AMR
		fp = fopen("kelvin08-distance.raw.amr", "rb");
	#endif
	if (fp == NULL){
		printf("DistanceField not found\n");
		exit(1);
	}
	#ifdef MG
	read_xyz_float_v2(fp, d, LEVEL);
	#endif	
	#ifdef AMR
	read_foreach_float_v2(fp, d, LEVEL);
	#endif

	fclose(fp);
	  
	vertex scalar phi[];
	foreach_vertex(){
		phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
		d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
	}
	  
	boundary ({phi});		
	fractions (phi, cs);
	boundary ({cs});
	foreach()
		cs[]=1.-cs[];	
	boundary ({cs});
	
	int i = 0;
	double t = 0.;
	
	scalar mpi_color[];
	foreach(){
		mpi_color[]=pid();		
	}
	boundary(all);
	#ifdef AMR
	{	  
		char path[]="htg/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_htg((scalar *) {cs, mpi_color},(vector *){NULL}, path, prefix, i, t);
	}
	#endif
	#ifdef MG
	{
		char path[]="vti/";
		char prefix[80];
		sprintf(prefix, "data_%06d", i);  
		output_vti((scalar *) {cs, mpi_color},(vector *){NULL}, path, prefix, i, t);
	}	
	#endif

}
