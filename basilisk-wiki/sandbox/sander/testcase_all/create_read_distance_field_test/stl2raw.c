/**
 This Programm reads a STL-File and outputs a distance field 
 in N-order (foreach() operator)
 
 
 */
//================================
#define MG
#ifdef MG
#include "grid/multigrid3D.h"
#endif

//~ #define AMR
//~ #ifdef AMR
//~ #include "grid/octree.h"
//~ #endif
//================================

#include "distance.h"
int LEVEL = 7;


int main(int argc, char * argv[])
{ 
	if(argc > 1){
		LEVEL = atoi(argv[1]);
	}   
        #ifdef MG
        printf("Prozess:%04i: Erzeuge Distanzfeld aus STL für Multigrid auf Level %i\n", pid(), LEVEL);
	#endif
        #ifdef AMR
	printf("Prozess:%04i: Erzeuge Distanzfeld aus STL für AMR auf Level %i\n", pid(), LEVEL);
	#endif
	unsigned long long nvox = (1 << LEVEL);	
	init_grid (nvox);
	size (1.);

	foreach_dimension()
		periodic(right);

	FILE * fp ; 

	fp = fopen("../../../../../kelvinstl/kelvin08.stl", "r");
	if (fp == NULL){
		printf("STL not found\n");
		exit(1);
	}
	coord * p = input_stl (fp);
	fclose(fp);

	coord min, max;
	bounding_box (p, &min, &max);
	double maxl = -HUGE;
	foreach_dimension()
		if (max.x - min.x > maxl)
			maxl = max.x - min.x;

	scalar d[];
	distance (d, p);

	int i=0;
	float* data = (float*) malloc(sizeof(float)*nvox*nvox*nvox);
	foreach(){
		float value[1] = {d[]};
		data[i]=value[0];
		i++;
	}

	#ifdef MG
	fp = fopen("kelvin08-distance.raw.mg", "wb");
	#endif
	#ifdef AMR
	fp = fopen("kelvin08-distance.raw.amr", "wb");
	#endif
	fwrite(data, sizeof(float), nvox*nvox*nvox, fp);	
	fclose(fp);  

	free(data);
}
