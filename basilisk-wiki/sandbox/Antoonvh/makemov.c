/**
# Generate a movie from turbulence data
   
3D flow data is fetched from the [Johns Hopkins turbulence
database](http://turbulence.pha.jhu.edu/). The u,v and w data are
assigned to the r,g and b channels of a `.mp4` movie.

<div class="figure">
<video controls loop preload="metadata" width="600px">
<source src="https://surfdrive.surf.nl/files/index.php/s/cR4D3acmqxV0Q4N/download" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Music composed by Eric Satie, performed by Reinbert de Leeuw.
</p>
</div>
*/
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "turblib.h"

#define N 1024

int main() {
  char * authtoken = "---"; //Your token goes here
  char * dataset = "isotropic1024coarse";
  
  soapinit();
  turblibSetExitOnError(1);
  float * datu = (float*)malloc (3*N*N*sizeof(float));
  /**
     Convert data to .mp4 asap with avconv (or ffmpeg) using Sephane
     Popinet's Bash script.
     
     see www.basilisk.fr/src/ppm2mp4
  */
  FILE * fp = popen ("ppm2mp4 turbulence.mp4", "w");
  for (int z = 1; z <= N; z++) {
    getCutout (authtoken, dataset, "u", 1000, 1, 1, z,
	       N, N, z, 1, 1, 1, 1, &datu[0]);
    fprintf (fp, "P6\n%d %d\n%d\n", N, N, 255); // .ppm header 
    for (int j = 0; j < N*N*3 ; j += 3) {
      unsigned char u[3]; // 3 byte pixel
      for (int i = 0; i < 3; i++)
	u[i] =  (unsigned char)fmin(fabs(datu[j + i])*255, 255); // irriversible
      //u[i] =  (unsigned char)fmax(fmin(datu[j + i]*128 + 128, 255), 0);   
      fwrite (&u[0], 1, 3, fp);
    }
    printf ("frame %d\n", z);
  }
  /**
     clean up the mess
  */
  fclose (fp);
  soapdestroy();
  free (datu);
  return 0;
}