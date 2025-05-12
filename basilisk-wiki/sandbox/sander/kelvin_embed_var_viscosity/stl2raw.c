/**
 * This Programm reads a STL-File and outputs a distance field 
 * x,y,z Order

 */

#include "grid/octree.h"
#include "distance.h"

int LEVEL = 7;

int main(int argc, char * argv[])
{ 
	if(argc > 1){
		LEVEL = atoi(argv[1]);
	}
  assert(npe() == 1);
   
	init_grid (1 << LEVEL);
	size (1.);
	/**
	* Periodische Randbedingung in allen Raumrichtungen
	*/
	foreach_dimension()
		periodic(right);

	FILE * fp ; 

	fp = fopen("kelvin08.stl", "r");
	if (!fp){ printf("STL not found\n"); exit(1); }
	coord * p = input_stl (fp);
	fclose(fp);

	coord min, max;
	bounding_box (p, &min, &max);

	scalar d[];
	distance (d, p);
  
  //~ restriction({d});
  unsigned long long cell_count = 0, offset=0;
  for (int lvl=0;lvl<=LEVEL;++lvl)
    cell_count += cube(1<<lvl);
    
  float* data = qmalloc(cell_count,float);
  
  for (int lvl=0;lvl<=LEVEL;++lvl) {
    fprintf(stderr,"Write Data on Level %i\n", lvl);
    offset += cube(1<<(lvl-1));
    fprintf(stderr,"Offset %llu\n", offset);
    foreach_level(lvl,serial) {
      int cells_on_axis = 1 << lvl;            
      int i = (x - X0)*cells_on_axis/L0,\
          j = (y - Y0)*cells_on_axis/L0,\
          k = (z - Z0)*cells_on_axis/L0;
      unsigned long long index = offset \
                              + (cells_on_axis*k + j)*cells_on_axis + i;
      #if DEBUG
      fprintf(stderr,"stl2raw.c, lvl:%i i:%i j:%i k:%i index:%llu cell_count:%llu\n",\
            lvl,i,j,k,index,cell_count);
      fflush(stderr);
      #endif
      float value[1] = {d[]};
      data[index]=value[0];            
    }  
	}
  
	fp = fopen("kelvin08-distance.raw.tree_structured", "wb");
	if (!fp) { fprintf(stderr,"File could not be written to!\n"); exit(1);}
	fwrite(&data[0], sizeof(float), cell_count, fp);	
	fclose(fp);  

  fp = fopen("kelvin08-distance.raw.last_level", "wb");
  if (!fp) { fprintf(stderr,"File could not be written to!\n"); exit(1);}   
  /** Set array pointer to only write the largest level*/
  fwrite(&data[offset], sizeof(float), cell_count - offset, fp);	
  fclose(fp);  
	free(data);
}
