#include <stdio.h>


// faire un coup de for s in {a, b, c} comme dans Basilisk_C
void output_gnuplot_heatmap (scalar field, int N,
			     char filename[],
			     double min, double max, double box[2][2],
			     int masked, scalar mask,
			     char xlabel[], char ylabel[], char cblabel[], char title[])
{
  if ((min == 0.) && (max == 0.)) { // If min and max are 0, we compute the min and max values of the field to define the colormap
    if (masked) { // If a mask is used, we don't use the masked values to calibrate the colormap
      min = 1e10;
      max = -1e10;
      foreach (reduction(min:min) reduction(max:max)) {
	if (mask[] >= 0.) {
	  min = min (min, field[]);
	  max = max (max, field[]);
	}
      }
    } else {
      min = 1e10;
      max = -1e10;
      foreach (reduction(min:min) reduction(max:max)) {
	min = min (min, field[]);
	max = max (max, field[]);
      }
    }
  }

  // Write data
  FILE * fp = fopen ("plot.data", "w");
  FILE * fpp = fopen ("mask", "w");
  int ii = 0;
  foreach (serial) {
    double color = max (min, min (field[], max));
    if ((masked) && (mask[] < 0.))
      fprintf (fpp, "%g %g 0x00FFFFFF\n", x, y);
    else
      fprintf (fpp, "%g %g 0xFF000000\n", x, y);
    fprintf (fp, "%g %g %g %g\n", x, y, field[], color);
    ii += 1;
    if (ii % N == 0) {
      fprintf (fp, "\n");
    }
  }
  fprintf (fp, "\n");
  fclose (fp);
  fclose (fpp);

  // Write the gnuplot script
  fp = fopen ("plot.plot", "w");
  fprintf (fp,
	   "reset session\n"
	   "min = %g\n"
	   "max = %g\n"
	   "set term pngcairo size 1024,1024\n"
	   "set output '%s.png'\n"
	   "set view map scale 1\n"
	   "set size square\n"
	   "set pm3d implicit\n"
	   //"set palette defined (0 '#000022', 1 '#32369c', 2 '#00b9a1', 3 '#3dd872', 4 '#b9f18b', 5 '#fefe98', 6 '#c0ae77', 7 '#825f56', 8 '#926f66')\n"
	   "set palette defined (0 '#000088', 1./8 '#0000ff', 3./8 '#00ffff', 5./8 '#ffff00', 7./8 '#ff0000', 1 '#880000')\n"
	   "unset key\n"
	   "unset title\n"
	   "set xlabel '%s'\n"
	   "set ylabel '%s'\n"
	   "set cblabel '%s'\n"
	   "set title '%s'\n"
	   "set autoscale noextend\n"
	   "set cbrange [min:max]\n"
	   "file = 'plot.data'\n"
	   "splot file index (0) using 1:2:3:4 with pm3d palette",
	   min, max, filename, xlabel, ylabel, cblabel, title);
  if (masked)
    fprintf (fp, ", 'mask' using 1:2:($0):3 with points pointtype 5 pointsize 1. lc rgb var\n");
  fclose (fp);
  
  system ("gnuplot plot.plot");
}
