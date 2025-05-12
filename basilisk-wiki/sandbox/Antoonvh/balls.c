/**
# Citrus-slicing dream

I dreamt about slicing citrus fruit.

![It looked exactly like this](balls/ballz.mp4)(loop)
 */
#include "grid/multigrid3D.h"
#include "bwatch.h"

#define LEVEL 8

vector v[], v2[];
scalar s[];
unsigned char mov [1 << LEVEL][1 << LEVEL][1 << LEVEL][3];
int main() {
  N = 1 << LEVEL;
  X0 = Y0 = Z0 = -L0/2;
  init_grid (N);
  /**
## Load another movie and image data
 */
  system ("wget https://antoonvanhooft.nl/media/subsiding.mp4");
  system ("mv subsiding.mp4 m.mp4");
  char cmd[999];
  sprintf (cmd, "ffmpeg -i m.mp4 -vf scale=%d:%d -y r.mp4", N, N);
  system (cmd);
  system ("ffmpeg -r 1 -i r.mp4 -r 1 frame%04d.ppm");
  
  char line_one[2];
  int height, width, max_color;
  char fn[99];
  for (int ti = 1; ti <= N; ti++) {
    sprintf (fn, "frame%04d.ppm", ti);
    FILE * fp = fopen (fn, "rb");
    fscanf(fp, "%s %d %d %d\n", line_one, &height, &width, &max_color);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	fread(&mov[i][j][ti - 1][0], 1, 3, fp);
      }
    }
    fclose (fp);
  }
  system ("rm -f r.mp4 m.mp4 frame*");
  foreach() {
    s[] = x + y + z;
    v.x[] = mov[point.i - GHOSTS][point.j-GHOSTS][point.k-GHOSTS][0];
    v.y[] = mov[point.i - GHOSTS][point.j-GHOSTS][point.k-GHOSTS][1];
    v.z[] = mov[point.i - GHOSTS][point.j-GHOSTS][point.k-GHOSTS][2];
  }
  system ("wget -O img.jpg \"https://www.seriouseats.com/thmb/dIqfr4qwvFYvTZEVB8meWTL64g0=/610x0/filters:no_upscale():max_bytes(150000):strip_icc():format(webp)/__opt__aboutcom__coeus__resources__content_migration__serious_eats__seriouseats.com__images__2014__04__20140421-knife-skills-citrus-10-c7bac5a6bd53499fa64e3b65f21a052c.jpg\"");
  sprintf (cmd, "convert img.jpg -resize %dx%d^ \
                 -gravity center -extent %dx%d img.ppm", N, N, N, N);
  system (cmd);
  FILE * fpi = fopen ("img.ppm", "rb");
  fscanf(fpi, "%s %d %d %d\n", line_one, &height, &width, &max_color);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      fread(&mov[j][255 - i][0][0], 1, 3, fpi);
    }
  }
  fclose (fpi);
  foreach() {
    v2.x[] = mov[point.i - GHOSTS][point.j-GHOSTS][0][0];
    v2.y[] = mov[point.i - GHOSTS][point.j-GHOSTS][0][1];
    v2.z[] = mov[point.i - GHOSTS][point.j-GHOSTS][0][2];
  }
  boundary ({s, v});
  watch (fov = 1.2, nx = 1024, ny = 768);
  FILE * fp = popen ("ppm2mp4 ballz.mp4", "w");
  for (double j = 0; j < 2*pi; j += pi/100. + 1e-6) {
    sphere (C = {0.05 + 0.25 *sin(j), 0.  + 0.05*cos(3*j), 0.1 - 0.3 *cos(j)}, R = 0.1,
	    mat = {.v = v, .min = 0, .max = 255, .linear = true});
    sphere (C = {0.05 + 0.3 *sin(j - 1.5), 0.  + 0.05*cos(2*(j - 1)), 0.1 - 0.25 *cos(j - 1.5)}, R = 0.1,
	    mat = {.T = 0.7});
    sphere (C = {0.02 + 0.25*sin(j - 4), 0. + 0.1*cos(4*(j - 2)), 0.1 - 0.2*cos(j - 4)}, R = 0.1,
	    mat = {.col = {100, 100, 100}, .R = 0.5});
    sphere (C = {0.25*cos(j)           ,  0.3, -0.24*sin(j)},     R = 0.1,
	    mat = {.s = s, .R = 0.4, .min = 0, .max = 1.4});
    sphere (C = {-0.25*cos(j)           ,  -0.3, -0.24*sin(j)},     R = 0.10,
	    mat = {.v = v, .R = 0.2, .min = 0, .max = 255, .linear = true});
    sphere (C = { -0.24*sin(j - 3)                  ,  -0.25*cos(j - 3), 0.4},     R = 0.10,
	    mat = {.ind = 1.3});
    sphere (C = { -0.24*sin(-j + 2)                  ,  0.25*cos(j - 2), -0.4},     R = 0.10,
	    mat = {.R = 1});
    sphere (C = {0, -0.3*cos(j)}, R = 0.1,
	    mat = {.ind = 1.3});
    quadriangles (alpha = -0.49,
		  mat = {.v = v2, .min = 0, .max = 255});
    sphere (6, mat = {.dull = true});
    store (fp);
    plain();
    printf ("frame %g\n", j);
  }
}
