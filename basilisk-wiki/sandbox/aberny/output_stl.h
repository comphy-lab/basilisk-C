/**
# STL output

This function output a scalar field into an stl file by using the facets
algoritm. The function write the data into a binary file so that the resuting
STL file can be read by any stl-reader (such as meshlab for example).

The input parameters  are a volume fraction field c and an optional file name.
By default the name will be "out.stl"

This output works in serial, with openmp and with MPI
*/
#if dimension >=3
#include "fractions.h"

trace
void stl_output_binary (scalar c, const char * name = "out.stl")
{

  /**
  For the MPI output, we open as many file as there is process.
  We will use a cat command at the end and then delete those external file*/

#if _MPI
  char outMPI[80];
  sprintf(outMPI, "_outSTL-%d.stl", pid() );
  FILE * fp = fopen(outMPI, "wb");
#else
  FILE * fp = fopen(name, "wb");
#endif

  /** We write the data into a binary file. The header contains a string of 80
chars that does not start with "solid". We put an empty one.*/

  char title_string[80];
  face vector s;
  s.x.i = -1;
  int nfacets = 0;

  if (pid()==0){ //Write in the 0 file
    fwrite( title_string, 1, sizeof(title_string), fp); //Title
    fwrite( &nfacets, sizeof(int), 1, fp); //Number of facets - we will update this later
  }

  foreach(reduction(+:nfacets)){
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (c[], n);
      coord v[12];
      int m = facets (n, alpha, v, 1.); //m is the number of coordinates
      //Assume all values in v belong to the same plane.
      for (int j = 0; j < m - 2; j++){
        /** 
        In an STL file, there is first the normal of a surface, then the
        cordinate of the triangle and then a short int equal to 0*/
        float word[12] = { n.x, n.y, n.z,
           x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta, 
           x + v[j+1].x*Delta, y + v[j+1].y*Delta, z + v[j+1].z*Delta, 
           x + v[j+2].x*Delta, y + v[j+2].y*Delta, z + v[j+2].z*Delta};
        unsigned short *trmpt = (unsigned short*)(&(word[12]) + 1);
        *trmpt = 0;
        fwrite(word, 12*sizeof(float) + sizeof(unsigned short), 1, fp);      
        nfacets += 1;
      } 
    }
  }
#if _MPI
  if (pid() == 0) {
#endif
    // fprintf(stderr, "nfacets=%d\n", nfacets);
    // fprintf(stderr, "Rewriting headers...\n");

    /**
      Update the number of facetes at the top of the file*/
    fseek(fp, 80, SEEK_SET);
    fwrite( &nfacets, sizeof(int), 1, fp);
#if _MPI
  }
  fclose(fp);
  /**
    Concatenate all the file into a single file*/
  char copyCat[99];
  sprintf(copyCat, "cat _outSTL-*.stl > %s", p.name);
  if (pid() == 0){
    system(copyCat);
    system("rm -f _outSTL-*.stl");
  }
#else
  /**
  Don't forget to close the file*/
  fclose(fp);
#endif
}
#endif

/**
An example using this output function can be found here: [csgBool.c](csgBool.c)

The stl file can be open with meshlab, for example.
*/