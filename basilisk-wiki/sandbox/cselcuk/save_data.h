/** 
# Wrapper for output-functions with Paraview
*/

# include "output_vtu_foreach.h"

void save_data(scalar * list, vector * vlist, int const iter, double const time,  char * name){
 
  FILE * fpvtk;
  char filename[80];
  char suffix[80];

  sprintf (filename,"%s", name);
  
#if  _MPI
  sprintf (suffix, "-%03d_n%3.3d.vtu", iter, pid());
  strcat(filename, suffix);
#else
  sprintf (suffix, "-%03d.vtu", iter);
  strcat(filename, suffix);
#endif
 
  fpvtk =  fopen (filename, "w");
  /* N = 1 << maxlevel; */ 
  output_vtu_bin_foreach (list, vlist, fpvtk, false);
  fclose(fpvtk);
   
#if _MPI
  if (pid() ==0 ) {
    sprintf (filename, "%s",name);
    sprintf(suffix, "-%03d.pvtu", iter);
    strcat(filename, suffix);
    fpvtk = fopen(filename, "w");
    
    sprintf (filename, "%s",name);
    sprintf(suffix, "-%03d", iter);
    strcat(filename, suffix);
    output_pvtu_bin (list, vlist, fpvtk, filename);
    fclose (fpvtk);
  }
#endif
  
  
}
