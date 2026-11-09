/**
# GPX tracks to heat map

Making a heatmap from GPX data tracks

![Strava in the wiki](https://www.antoonvanhooft.nl/media/heatmap_all.png)


*/
#include "distance.h" // for input_xy()
#include "utils.h" // output_ppm

// Color bar
void blue_body (double cmap[NCMAP][3])
{
  int ref_points = 4;
  static double basemap[4][3] = {
    {0.5, 0.5, 1},
    {0.1, 0, 1},
    {0.1, 0.0, 0.5},
    {0.1, 0.0, 0}
  };
  
  for (int i = 0; i < NCMAP; i++) {
    double x = i*((ref_points-1) - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

/**
## Segmet cell intersections
   
The  foreach_segment() macro is too [buggy](https://basilisk.fr/sandbox/bugs/segment_finder.c), but I can't see what is wrong. So we modify some bwatch methods, where finite rays are used to represent segments. Furthermore, the multigrid acceleration reduces computation time significantly. 
 */

typedef struct ray {
  coord O;
  coord dir;
  double length;
} ray;

static inline double ray_cell_intersect (ray r, Point point, coord a[2]) {
  coord cc = {x,y,z};
  double in = -HUGE, out = HUGE;
  foreach_dimension() {
    in  = max(in,  r.dir.x != 0 ? ((cc.x - (sign(r.dir.x)*Delta/2.)) - r.O.x)/r.dir.x : -HUGE);
    out = min(out, r.dir.x != 0 ? ((cc.x + (sign(r.dir.x)*Delta/2.)) - r.O.x)/r.dir.x :  HUGE);
  }
  if (in >= out || out < 0 || in > r.length) 
    return HUGE;
  in = in < 0 ? 0 : in; //The origin is in the cell.
  foreach_dimension() {
    a[0].x = r.O.x +  in*r.dir.x;
    a[1].x = r.O.x + out*r.dir.x;
  }
  return in;
}

macro foreach_ray_cell_intersection (ray r) {
  foreach_cell() { //MG acceleration
    coord _a[2] = {0};
    double _dist = ray_cell_intersect (r, point, _a);
    if (_dist < HUGE) {
      if (is_leaf(cell))
	{...}
    } else // Ray missed, skip children
      continue;
  }
}

scalar s[];


int main(int argc, char ** argv) {
  int counter = 1;
  init_grid (1024);
  
  // We find the bounds of our canvas
  X0 = HUGE;
  Y0 = HUGE;
  double X1 = -HUGE;
  double Y1 = -HUGE;
  
  while (counter++ < argc) {
    char srt[999];
    sprintf (srt, "./gpx2gnu.sh %s", argv[counter-1]);
    system (srt);
    coord * p = input_xy (fopen ("track3.gnu", "r"));
    coord minp, maxp;
    
    bounding_box (p, &minp, &maxp);
    X0 = min(X0, minp.x);
    Y0 = min(Y0, minp.y);
    X1 = max(X1, maxp.x);
    Y1 = max(Y1, maxp.y);
  }
  double YM = (Y0 + Y1)/2; // central lattitude
  double fac = cos(YM*pi/180.); // mapping factor
  double margin = 1.1;
  X1 *= fac;
  X0 *= fac;
  L0 = margin* max(X1 - X0, Y1 - Y0);
  X0 = (X1 + X0)/2 - (1/2.*L0);
  Y0 = YM - (1/2.*L0);
  
  printf ("%g %g %g\n", L0, X0, Y0);
  counter = 1;
  while (counter++ < argc) {
    char srt[999];
    sprintf (srt, "./gpx2gnu.sh %s", argv[counter-1]);
    printf ("%s \n \n ", srt);

    system (srt);
    coord * p = input_xy (fopen ("track3.gnu", "r"));
    // Read the polyline data
    while (p->x != nodata) {
      coord S[2], P[2];
      S[0] = *p;
      S[0].x *= fac;
      p++;
      if (p->x != nodata) {
	S[1] = *p;
	S[1].x *= fac;
	ray r;
	r.O = S[0];
	coord ta = {0, 0, 0};
	double len = 0;
	foreach_dimension() {
	  ta.x = S[1].x - S[0].x;
	  len += sq(ta.x);
	}
	r.length = sqrt(len);
	// We only want to plot relevant segments
	if (len > 0 && len < 1) { // corrupt data
	  normalize(&ta);
	  r.dir =  ta;
	  double cells = 0;
	  // Look for errouneous ray data
	  foreach_ray_cell_intersection(r) {
	    if (_dist >= 0) 
	      cells++;
	  }
	  if (cells < N/2.) { // exclude errouneous rays
	    foreach_ray_cell_intersection(r) {
	      // To prevent double counting we do not increase the starting cell
	      if (_dist > 0) { 
		// double len = 0;
		// foreach_dimension()
		//  len += sq(_a[0].x - _a[1].x);
		s[]++; //= sqrt(len)/Delta;
	      }
	    }
	  }
	}
      }
    }
  }
  // The often single-cell tracks need to me more wide
  scalar h[];
  foreach() {
    h[] = 0;
    if (s[] < 0.1) {
      double sum = 0;
      int n = 0;
      foreach_neighbor(1)
	if (s[] >= 0.1) {
	  sum += s[];
	  n++;
	}
      if (n)
	h[] = sum/n;
    }
  }
  // Widening itself and masking.
  foreach() {
    if (h[])
      s[] = h[];
    if (s[] < 0.1)
      s[] = nodata;
  }
  // The 1-2-3-2-1 triangle-filter kernel
  foreach() {
    int n = 0;
    double sum = 0;
    if (s[] != nodata) {
      sum = s[];
      n = 1;
      foreach_neighbor(1) {
	if (s[] != nodata) {
	  sum += s[];
	  n++;
	}
      }
      foreach_neighbor(2) {
	if (s[]  != nodata) {
	  sum += s[];
	  n++;
	}
      }
    }
    if (n)
      h[] = sum/n;
    else
      h[] = nodata;
  }
  // apply transformation for color scale.
  double maxv = 2;
  foreach() {
    if (h[] != nodata)
      h[] = min(log(h[]), maxv);
  }
  output_ppm (h, file = "track.png" , min = -0.1, max = maxv + 0.1, map = blue_body);
  // add transparancy
  system ("convert track.png -transparent black heatmap.png");
  char cmd[999];
  // We need to un-map the x coordinate for the osm python script (see below)
  sprintf (cmd, "./venv/bin/python osm.py %g %g %g %g\n", X0/fac, Y0, (X0 + L0)/fac, Y0 + L0);
  system(cmd);
}
/**
## Supporting scripts

### gpx to poly line

~~~bash
   #! /usr/bin/bash
   awk 'match($0, /lat="([^"]+)".*lon="([^"]+)"/, a) {
                      lat = a[1]
      		     lon = a[2]
		     print lon, lat }' $1 > track3.gnu
~~~

### OSM interface
   
This was vibe coded, well done chatGPT! 
   
~~~python
   #!/usr/bin/env python3
import math
from io import BytesIO
from PIL import Image
import mercantile
import requests
import sys


# Define your bounding box: [min_lon, min_lat, max_lon, max_lat]
bbox = (float(sys.argv[1]), float(sys.argv[2]) , float(sys.argv[3]), float(sys.argv[4]))
zoom = 12  # Choose your zoom level

# Get all tile coordinates intersecting the bounding box at this zoom level
tiles = list(mercantile.tiles(bbox[0], bbox[1], bbox[2], bbox[3], zoom))

min_x = min(t.x for t in tiles)
max_x = max(t.x for t in tiles)
min_y = min(t.y for t in tiles)
max_y = max(t.y for t in tiles)

# Calculate dimensions of the final stitched image
width = (max_x - min_x + 1) * 256
height = (max_y - min_y + 1) * 256

canvas = Image.new("RGB", (width, height))

# Headers required by OSM tile usage policy (must include a valid User-Agent)
headers = {"User-Agent": "heatmap script"
           "contact: your@mail.com}

for t in tiles:
  url = f"https://tile.openstreetmap.org/{zoom}/{t.x}/{t.y}.png"
  response = requests.get(url, headers=headers)
  if response.status_code == 200:
    tile_img = Image.open(BytesIO(response.content))
    # Paste tile into the correct relative position on the canvas
    x_offset = (t.x - min_x) * 256
    y_offset = (t.y - min_y) * 256
    canvas.paste(tile_img, (x_offset, y_offset))
  else:
    print(f"Failed to fetch tile {zoom}/{t.x}/{t.y}")



# 1. Get the pixel coordinates of your exact bounding box relative to the top-left tile (min_x, min_y)
ul_web = mercantile.xy(bbox[0], bbox[3])  # Top-left (min_lon, max_lat)
lr_web = mercantile.xy(bbox[2], bbox[1])  # Bottom-right (max_lon, min_lat)

# Convert web mercator coordinates to pixel coordinates at the given zoom level
# (tile size is 256px)
world_size = 256 * (2**zoom)


def xy_to_pixel(x, y):
  px = (x + 20037508.342789244) * world_size / 40075016.685578488
  py = (20037508.342789244 - y) * world_size / 40075016.685578488
  return px, py


ul_px, ul_py = xy_to_pixel(ul_web[0], ul_web[1])
lr_px, lr_py = xy_to_pixel(lr_web[0], lr_web[1])

# Calculate offset relative to the minimum tile's top-left corner
min_tile_origin_x = min_x * 256
min_tile_origin_y = min_y * 256

crop_left = int(ul_px - min_tile_origin_x)
crop_top = int(ul_py - min_tile_origin_y)
crop_right = int(lr_px - min_tile_origin_x)
crop_bottom = int(lr_py - min_tile_origin_y)

# 2. Crop the canvas to the exact bounds
exact_canvas = canvas.crop((crop_left, crop_top, crop_right, crop_bottom))
exact_canvas.save("osm_exact_output.png")
print("Saved cropped map as osm_exact_output.png!")
~~~

We can now resize and overlay the heatmap with the streetmap, using

   ~~~bash
w=$(identify -ping -format "%w" osm_exact_output.png)
h=$(identify -ping -format "%h" osm_exact_output.png)
convert osm_exact_output.png heatmap.png -resize ${w}x${h}\! -composite output.png
~~~
*/							  

							  

