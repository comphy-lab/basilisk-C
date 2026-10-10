/** This is almost the same file as in cselcuk sandbox [save_data.h](https://basilisk.fr/sandbox/cselcuk/save_data.h). \
There is just the possibility to store the data in a subdirectory.
*/
                                                                           
#include "output_vtu_foreach.h"

void save_data(scalar * list, vector * vlist, int const iter, double const time, char * prefix, char * directory){

  FILE * fpvtk;
  char filename[80];
  char filepref[80];
  char suffix[80];

  sprintf (filepref,"%s", prefix);
  sprintf (filename,"%s", directory);
  strcat(filename, filepref);

#if  _MPI
  sprintf (suffix, "-%03d_n%3.3d.vtu", iter, pid());
  strcat(filename, suffix);
#else
  sprintf (suffix, "-%03d.vtu", iter);
  strcat(filename, suffix);
#endif

  fpvtk =  fopen (filename, "w");
  output_vtu_bin_foreach (list, vlist, fpvtk, false);
  fclose(fpvtk);

#if _MPI
  if (pid() ==0 ) {
    sprintf (filepref,"%s", prefix);
    sprintf (filename,"%s", directory);
    strcat(filename, filepref);
    sprintf(suffix, "-%03d.pvtu", iter);
    strcat(filename, suffix);
    fpvtk = fopen(filename, "w");

    sprintf (filename,"%s", prefix);
    sprintf(suffix, "-%03d", iter);
    strcat(filename, suffix);
    output_pvtu_bin (list, vlist, fpvtk, filename);
    fclose (fpvtk);
  }
#endif
}